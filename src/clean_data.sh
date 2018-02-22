
#clean_data.sh
#Ajay Vibhute,
#Mayuri Shinde
#19 Sept 2017
#Script takes level2 file, mkf file, bunch file, saaconfig, noise config and generates 
#level2 cleaned file
#Steps :
#1.cztdataselection- for removing saa
#2.cztnoisypixclean- for removing noisy pixels
#3.cztsuperbunchclean - for removing super bunches 
#4.cztheavybunchclean - for removing heavy bunches 
#5.cztflickpixclean - for removing flickering pixels
#6.cztdatasel
#7.cztevtclean
#8.cztbindata


#datadir="/home/cztipoc/Mayuri/cztipoc/noise_reduction/data_for_testing/grbtest/orbit_04229"
#mkffile=$datadir/AS1G05_240T01_9000000536_04229czt_level2.mkf
#saaconfig=$datadir/saaThreshold
#eventfile=$datadir/AS1G05_240T01_9000000536_04229cztM0_level2.evt
#caldbbadpix=$datadir/AS1cztbadpix20160809v01.fits
#bunchfile=$datadir/AS1G05_240T01_9000000536_04229cztM0_level2_bunch.fits
#noiseconfigfile=$datadir/noiseReductionThreshold
#caldblld=$datadir/AS1cztlld20160517v01.fits

mkffile=$1
saaconfig=$2
eventfile=$3
caldbbadpix=$4
bunchfile=$5
noiseconfigfile=$6
caldblld=$7
bunchcleanfile=$8

echo "Executing data selection, will create multiple files";
echo "cztdataselection $mkffile $saaconfig $eventfile"


cztbunchclean par_infile=$eventfile par_bunchfile=$bunchfile par_outfile=$bunchcleanfile par_livetimefile=livetime_bc.fits par_skipT1=0.0 par_skipT2=0.0 par_skipT3=0.0 par_bunchdeftime=20 par_bunch_length_thresh=20 par_livetime_binsize=1.0 clobber=yes history=yes

cztdataselection $mkffile $saaconfig $bunchcleanfile

temp=`echo $bunchcleanfile| sed 's/.evt//'`
lsstr="$temp""_""*[0-9].evt"

echo "Executing cztnoisypix clean"
for f in `ls $lsstr`;
do
	echo "Running cztnoisypix clean for $f"
	echo "cztnoisypixclean $f $caldbbadpix $bunchfile $noiseconfigfile"

	cztnoisypixclean $f $caldbbadpix $noiseconfigfile
	noisecleanedfile=`echo $f|sed 's/.evt/_nc.evt/'`
	noisebadpixfile=`echo $f|sed 's/.evt/_badpix_nc.fits/'`
	noiselivetime=`echo $f|sed 's/.evt/_nc_livetime.fits/'`
	
	orbitno=`echo $noisecleanedfile |cut -d'_' -f 4|sed 's/cztM0//'`
	echo "Running cztsuperbunchclean for $noisecleanedfile"
	cztsuperbunchclean $noisecleanedfile $bunchfile livetime_bc.fits $orbitno 

	superbunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc.evt/'`
	superbunchlivetime=`echo $f|sed 's/.evt/_nc_sbc_livetime.fits/'`
	echo "Running cztheavybunchclean for $superbunchcleanevt"
	cztheavybunchclean $superbunchcleanevt $bunchfile $caldblld


	heavybunchcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc.evt/'`
	echo "Running cztflickpixclean for $heavybunchcleanevt"
	cztflickpixclean $heavybunchcleanevt $noisebadpixfile	 

	flickcleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc.evt/'`
	flickbadpixfile=`echo $f|sed 's/.evt/_badpix_nc_fc.fits/'`
	echo "Running czteventseperation for $flickcleanevt"
	czteventsep $flickcleanevt

	eventsepevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single.evt/'`
	dscleanevt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single_ds.evt/'`
	
	gtievt=`echo $f|sed 's/.evt/_nc_sbc_hbc_fc_single.gti/'`
	cztgtigen $eventsepevt mkffile=$mkffile outfile=$gtievt mkfThresholds.txt usergtifile=-
	
	echo "Running cztdatasel"
	cztdatasel infile=$eventsepevt gtifile=$gtievt gtitype="quad" outfile=$dscleanevt clobber="y" history="y"


	cleanevt=`echo $f|sed 's/.evt/_quad.evt/'`
	echo "Running cztevt clean"
	
	cztevtclean infile=$dscleanevt outfile=$cleanevt alphaval="0" vetorange="0-0" clobber="y" history="y" 
		
	echo "Running cztbindata"
	j=0
	for i in 0-1000
	do
		bindataout=`echo $f|sed 's/.evt/_/'`
		bindataout=$bindataout"E"$i"_quad_clean"
		wtevt=`echo $f|sed 's/.evt/.wevt/'`
		echo $bindataout $wtevt 

		cztbindata inevtfile=$cleanevt mkffile=$mkffile  badpixfile=$flickbadpixfile livetimefile=$superbunchlivetime outfile=$bindataout outevtfile=$wtevt maskWeight="no" rasrc="1" decsrc="1" badpixThreshold=0 outputtype="lc" timebinsize="1" energyrange=$i clobber="y"

	done
done	




