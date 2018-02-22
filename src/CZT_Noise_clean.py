import os
import numpy as np
import numpy.ma as ma
from astropy.io import fits,ascii



#creating light curves
def light_curve(time_UT,energy, bin_size,energy_1,energy_2):
	bin_data = np.arange(UT,UT_last,bin_size)
	energy_index =[(energy<energy_1)|(energy>energy_2)][0]
	time_UT = ma.masked_array(time_UT, mask = energy_index)
	light_curve, bin_e =  np.array(np.histogram(time_UT.compressed(), bins = bin_data))
	light_curve = light_curve/bin_size
	print len(time_UT.compressed())
	return light_curve, bin_data
	

outfile_1= open("CZT_raw_lc.txt","w")
outfile_2= open("CZT_cleaned_lc.txt","w")


evt_file = fits.open("AS1G05_234T02_9000000474_03632cztM0_level2.evt")
bunch_file = fits.open("AS1G05_234T02_9000000474_03632cztM0_level2_bunch.fits")
Q_n = 4

#############################Reading Event file
time_evt = evt_file[Q_n].data['TIME']
energy_evt = evt_file[Q_n].data['ENERGY']
det_id = evt_file[Q_n].data['DetID']
pix_id = evt_file[Q_n].data['pixID']
pixel_id = det_id*256+pix_id

UT = int(time_evt[0])
UT_last = int(time_evt[len(time_evt)-1])
#print str(UT)+'\t'+str(UT_last)+'\n'
#print len(time_evt)

#time_evt,energy_evt,det_id,pix_id,pixel_id = zip(*sorted(zip(time_evt,energy_evt,det_id,pix_id,pixel_id)))
time_evt = np.array(time_evt)
energy_evt = np.array(energy_evt)
det_id = np.array(det_id)
pix_id = np.array(pix_id)
pixel_id = np.array(pixel_id)

################################Reading Bunch_file
time_bun = bunch_file[Q_n].data['TIME']
no_of_bun_events = bunch_file[Q_n].data['NumEvent']
bun_det_id_1 = bunch_file[Q_n].data['DetId1']
bun_det_id_2 = bunch_file[Q_n].data['DetId2']
bun_det_id_3 = bunch_file[Q_n].data['DetId3']
bun_det_id_4 = bunch_file[Q_n].data['DetId4']


#time_bun,no_of_bun_events,bun_det_id_1,bun_det_id_2,bun_det_id_3,bun_det_id_4 = zip(*sorted(zip(time_bun,no_of_bun_events,bun_det_id_1,bun_det_id_2,bun_det_id_3,bun_det_id_4)))
no_of_bun_events = np.array(no_of_bun_events)
time_bun = np.array(time_bun)
bun_det_id_1 = np.array(bun_det_id_1)
bun_det_id_2 = np.array(bun_det_id_2)
bun_det_id_3 = np.array(bun_det_id_3)
bun_det_id_4  = np.array(bun_det_id_4 )



print("Total events = "+str(len(time_evt)))

###############################Noisy pixel flagging
total_pix_counts = [None]*4096
total_pix_counts = np.array(total_pix_counts)
for i in range(0,4096):
	total_pix_counts[i] = 0

for i in range(0,len(time_evt)):
	#outfile.write(str(energy_evt[i])+"\n")
	total_pix_counts[pixel_id[i]] = total_pix_counts[pixel_id[i]]+1



max_bun_det = [None]*len(time_bun)
for i in range(0,len(time_bun)):
	det = []
	det.append(bun_det_id_1[i])
	det.append(bun_det_id_2[i])
	det.append(bun_det_id_3[i])
	det.append(bun_det_id_4[i])
	det = np.array(det)
	g	= np.bincount(det)
	max_bun_det[i]	=	np.argmax(g)
	
		
print "Finished Reading"	

flagged_pixels = np.array([total_pix_counts==0])
#print flagged_pixels
#print sum(flagged_pixels[0])
total_pix_counts = ma.masked_array(total_pix_counts, mask=flagged_pixels)
total_mean = np.mean(total_pix_counts)
total_sigma = np.std(total_pix_counts)
total_mean_temp = total_mean
total_sigma_temp = total_sigma


while True:
	temp = np.array([total_pix_counts > (total_mean+5*total_sigma)])
	#print temp
	#print sum(temp[0])
	flagged_pixels = np.add(flagged_pixels,temp)
	#print sum(flagged_pixels[0])
	total_pix_counts = ma.masked_array(total_pix_counts, mask=flagged_pixels)
	total_mean = np.mean(total_pix_counts)
	total_sigma = np.std(total_pix_counts)
	#print total_mean
	if (total_mean_temp-total_mean) < 0.1 and (total_sigma_temp-total_sigma)<0.1:
		break
	total_mean_temp = total_mean
	total_sigma_temp = total_sigma	

#print total_mean
#print sum(flagged_pixels[0])
#print len(flagged_pixels[0])

#print [pixel_id == 10000]
flagged_pixels_index_noisy = np.array([pixel_id == 10000])
#print flagged_pixels_index[0]
#print len(flagged_pixels_index)
#print sum(flagged_pixels_index[0])
for i in range(0,len(flagged_pixels[0])):
	if flagged_pixels[0][i]==True:
		temp = np.array([pixel_id == i])
		flagged_pixels_index_noisy = np.add(flagged_pixels_index_noisy,temp)
#print flagged_pixels_index[0]
print ("Noise flagged events = "+str(sum(flagged_pixels_index_noisy[0])))
#print len(flagged_pixels_index[0])

print "Noise flagging finished"





############################Bunch_Noise_flagging#####################################
flagged_bun_evt =np.array([pixel_id == 10000])
temp_time = int(time_evt[0])
time_bun_curr = np.array(time_bun[np.nonzero((time_bun>int(time_evt[0]))&(time_bun<int(time_evt[0])+2))[0]])
#print time_bun_curr
#time_evt = ma.masked_array(time_evt, mask=flagged_pixels_index_noisy[0]) 
for i in range(0,len(time_evt)):
	if flagged_pixels_index_noisy[0][i] ==True:
		continue
	
	#print str(i)+"\t"+str(len(time_evt))
	#print str(i)
	
	if energy_evt[i]<30:

		if temp_time != int(time_evt[i]):
			time_bun_curr = np.array(time_bun[np.nonzero((time_bun>int(time_evt[i]))&(time_bun<int(time_evt[i])+2))[0]])
			temp_time = int((time_evt[i]))
		
		Bunch_time_sub = np.array(np.add(time_evt[i],-time_bun_curr))
		Bunch_time_sub =  Bunch_time_sub[Bunch_time_sub>0]
		Bunch_near_event = Bunch_time_sub[Bunch_time_sub<100000E-6]
		
		if len(Bunch_near_event)>0:
			#min_value = min(Bunch_near_event)
			#bunch_index = np.where(time_bun==(time_evt[i]-min_value))[0][0]
			#print bunch_index
			#if len(Bunch_near_event)>0 or max_bun_det[bunch_index]==detid[i]:
			flagged_bun_evt[0][i] = True


			
print "Finished Bunch cleaning"			
print ("bunch related events = "+str(sum(flagged_bun_evt[0])))








'''
##############Bunch related Noise cleaning1 ###################
flagged_evt = []
for i in range(0,len(time_evt)):
	if flagged_pixels_index[0][i] == True:
		continue
	l = len(time_bun)
	#print l
	y=[]
	for j in range(0,l):
		if (time_evt[i]-time_bun[j]) > 100E-6:
			#print "true"
			y.append(j)
			continue
			
		elif 0 <= (time_evt[i]-time_bun[j]) < 100E-6:
			# and energy_evt[i] < 30:
			flagged_evt.append(i)
			break
		elif (time_evt[i]-time_bun[j])<0:
			break
	#print len(y)	
	time_bun = np.delete(time_bun,y)			


print len(time_evt)
print len(flagged_pixels_index[0])
print len(flagged_evt)

##############################################################













##########################Bunch related events flagging2###############
event_30keV_index = [energy_evt < 30]
event_above_30keV_index = [energy_evt > 30]

time_evt_30 = ma.masked_array(time_evt, mask = event_above_30keV_index)
time_evt_30_above = ma.masked_array(time_evt, mask = event_30keV_index)
energy_evt_30 = ma.masked_array(energy_evt, mask = event_above_30keV_index)
energy_evt_30_above = ma.masked_array(energy_evt, mask = event_30keV_index)
det_id_30 = ma.masked_array(det_id, mask = event_above_30keV_index)
det_id_30_above = ma.masked_array(det_id, mask = event_30keV_index)
pixel_id_30 = ma.masked_array(pixel_id, mask = event_above_30keV_index)
pixel_id_30_above =ma.masked_array(pixel_id, mask = event_30keV_index)

time_evt_30 = ma.masked_array(time_evt, mask = flagged_pixels_index)
time_evt_30_above = ma.masked_array(time_evt, mask = flagged_pixels_index)
energy_evt_30 = ma.masked_array(energy_evt, mask = flagged_pixels_index)
energy_evt_30_above = ma.masked_array(energy_evt, mask = flagged_pixels_index)
det_id_30 = ma.masked_array(det_id, mask = flagged_pixels_index)
det_id_30_above = ma.masked_array(det_id, mask = flagged_pixels_index)
pixel_id_30 = ma.masked_array(pixel_id, mask = flagged_pixels_index)
pixel_id_30_above =ma.masked_array(pixel_id, mask = flagged_pixels_index)
print "y"

events_near_bunch_30keV_index = np.array([abs(time_evt_30-time_bun[0])<75E-3])
events_near_bunch_above_30keV_index = np.array([abs(time_evt_30_above-time_bun[0])<1E-3])
print "y"
#print sum(det_id_30)

for i in range(1,len(time_bun)):
	events_near_bunch_30keV_index = np.add(events_near_bunch_30keV_index, np.array([abs(time_evt_30-time_bun[i])<75E-3]))
	events_near_bunch_above_30keV_index = np.add(events_near_bunch_30keV_index, np.array([abs(time_evt_30_above-time_bun[i])<1E-3]))

time_evt_30 = ma.masked_array(time_evt, mask = events_near_bunch_30keV_index)
time_evt_30_above = ma.masked_array(time_evt, mask = events_near_bunch_above_30keV_index)
energy_evt_30 = ma.masked_array(energy_evt, mask = events_near_bunch_30keV_index)
energy_evt_30_above = ma.masked_array(energy_evt, mask = events_near_bunch_above_30keV_index)
det_id_30 = ma.masked_array(det_id, mask = events_near_bunch_30keV_index)
det_id_30_above = ma.masked_array(det_id, mask = events_near_bunch_above_30keV_index)
pixel_id_30 = ma.masked_array(pixel_id, mask = events_near_bunch_30keV_index)
pixel_id_30_above = ma.masked_array(pixel_id, mask = events_near_bunch_above_30keV_index)

print "y"
#print sum(det_id_30)
#######################################################################################
'''


##########################################Flagging flickering pixels####################
flagged_pixels_flickering = np.array([False]*4096)
bin_bound = np.arange(UT,UT_last,100)
event_50keV_index = [energy_evt < 50][0]
event_above_50keV_index = [energy_evt > 50][0]
for i in range(0,4096):
	pixel_index = [pixel_id != i]
	
	pixel_time_50 = ma.masked_array(time_evt, mask = np.add( np.add(flagged_bun_evt[0],pixel_index),np.add(event_above_50keV_index,flagged_pixels_index_noisy[0])))
	pixel_time_above_50 = ma.masked_array(time_evt, mask = np.add( np.add(flagged_bun_evt[0],pixel_index),np.add(event_50keV_index,flagged_pixels_index_noisy[0])))
	
	#pixel_time_30 = ma.masked_array(time_evt, mask = event_above_30keV_index)
	#pixel_time_above_30 = ma.masked_array(time_evt, mask =  event_30keV_index)
	
	pixel_lc_50 , bin_e = np.array( np.histogram(pixel_time_50.compressed(), bins = bin_bound))
	pixel_lc_above_50, bin_e = np.array(np.histogram(pixel_time_above_50.compressed(), bins = bin_bound))
	
	#pixel_lc_above_30 = np.array(pixel_lc_above_30)
	pixel_lc_50 = pixel_lc_50[pixel_lc_50!=0]
	pixel_lc_above_50 = pixel_lc_above_50[pixel_lc_above_50!=0]
	#print pixel_lc_30
	#outfile.write(str(pixel_lc_30)+"\n")
	#print len(pixel_lc_30)
	

		
	
	pixel_mean_50 = np.mean(pixel_lc_50)
	pixel_mean_above_50 = np.mean(pixel_lc_above_50)
	pixel_std_50 = np.std(pixel_lc_50)
	pixel_std_above_50 = np.std(pixel_lc_above_50)
	
	outliers_50 = [pixel_lc_50 > (pixel_mean_50+5*pixel_std_50)][0]
	outliers_above_50 = [pixel_lc_50 > (pixel_mean_above_50+5*pixel_std_above_50)][0]
	#print outliers_30
	no_of_outliers_50 = sum(outliers_50)
	no_of_outliers_above_50 = sum(outliers_above_50)
	
	
	if no_of_outliers_50 >1:
		# or no_of_outliers_above_50 >2:
		#print pixel_lc_50
		#print pixel_lc_above_50
		flagged_pixels_flickering[i] = True
		
		
flagged_pixels_index_flickering = np.array([pixel_id == 10000])
for i in range(0,len(flagged_pixels_flickering)):
	if flagged_pixels_flickering[i]==True:
		temp = np.array([pixel_id == i])
		flagged_pixels_index_flickering  = np.add(flagged_pixels_index_flickering,temp)
		
print sum(flagged_pixels_index_flickering[0])

		
print "No of pixels flagged as flickering = "+"\t"+str(sum(flagged_pixels_flickering))+"\n"




lc_1 , UT_lc = light_curve(time_evt, energy_evt,1,0,50)
lc_2 , UT_lc = light_curve(time_evt, energy_evt,1,50,100)
lc_3 , UT_lc = light_curve(time_evt, energy_evt,1,100,150)
lc_4 , UT_lc = light_curve(time_evt, energy_evt,1,150,10000)

q1=0
q2=0
q3=0
q4=0

for i in range(0,len(lc_1)):
	q1=q1+lc_1[i]
	q2=q2+lc_2[i]
	q3=q3+lc_3[i]
	q4=q4+lc_4[i]
	
	outfile_1.write(str(UT_lc[i]-UT)+"\t"+str(lc_1[i])+"\t"+str(lc_2[i])+"\t"+str(lc_3[i])+"\t"+str(lc_4[i])+"\n")


	

time_evt = ma.masked_array(time_evt, mask = np.add(flagged_pixels_index_noisy[0],np.add(flagged_bun_evt[0],flagged_pixels_index_flickering[0])))
energy_evt = ma.masked_array(energy_evt, mask = np.add(flagged_pixels_index_noisy[0],np.add(flagged_bun_evt[0],flagged_pixels_index_flickering[0])))
det_id = ma.masked_array(det_id, mask = np.add(flagged_pixels_index_noisy[0],np.add(flagged_bun_evt[0],flagged_pixels_index_flickering[0])))
pixel_id = ma.masked_array(pixel_id, mask = np.add(flagged_pixels_index_noisy[0],np.add(flagged_bun_evt[0],flagged_pixels_index_flickering[0])))




lc_1 , UT_lc = light_curve(time_evt, energy_evt,1,0,50)
lc_2 , UT_lc = light_curve(time_evt, energy_evt,1,50,100)
lc_3 , UT_lc = light_curve(time_evt, energy_evt,1,100,150)
lc_4 , UT_lc = light_curve(time_evt, energy_evt,1,150,10000)

q1=0
q2=0
q3=0
q4=0

for i in range(0,len(lc_1)):
	q1=q1+lc_1[i]
	q2=q2+lc_2[i]
	q3=q3+lc_3[i]
	q4=q4+lc_4[i]
	
	outfile_2.write(str(UT_lc[i]-UT)+"\t"+str(lc_1[i])+"\t"+str(lc_2[i])+"\t"+str(lc_3[i])+"\t"+str(lc_4[i])+"\n")
	
#print len(time_evt)	
print "Total genuine events"+str(q1+q2+q3+q4)+"\t"+str(q1)+"\t"+str(q2)+"\t"+str(q3)+"\t"+str(q4)+"\n"
