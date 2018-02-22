from __future__ import division
import os
import numpy as np
import numpy.ma as ma
from astropy.io import fits,ascii
import math
import sys
from matplotlib import pyplot as plt
import scipy.stats.distributions as ps

path_to_peak_file = "<path>/<parameter_name>/<parameter_value>/results_peaks.txt"
path_to_dip_file = "<path>/<parameter_name>/<parameter_value>/results_dips.txt"
path_to_orbit_report_file = "<path>/<parameter_name>/<parameter_value>/orbit_report.txt"
path_to_plots = "<path>/<parameter_name>/<parameter_value>/plots/" #make the plots directory manually
binsize = 1##########change binsize here################



def poissonian_limit(average,rms,sample_size):
	j = average + 3*rms 
	i = average - 3*rms
	while True:
	
		c = ps.poisson.pmf(int(average),int(j))*sample_size
		#print "High "+str(c)+"\t"+str(average)+"\t"+str(j)
		if c < 1.0:
			break 
		j = j+0.1*rms
		
	while True:
		c = ps.poisson.pmf(int(average),int(i))*sample_size
		#print "Low "+str(c)+"\t"+str(average)+"\t"+str(i)
		if c < 1.0:
			break 
		i = i-0.1*rms
		
	return j,i


peak_file = open(path_to_peak_file,"a+")
dip_file = open(path_to_dip_file,"a+")
orbit_report_file = open(path_to_orbit_report_file, "a+")



for j in range(1,5):
	file_lc_Q_E = fits.open(sys.argv[j])
	
	file_name = sys.argv[j]
	quad_no =j # int((j-1)/3)
	energy_range = sys.argv[5]#j%3

	###################################Reading light curves##################################
	#Q_E_time, Q_E_rate, Q_E_error, Q_E_fracexp = file_lc_Q_E[i].data["TIME"], file_lc_Q_E[i].data["RATE"], file_lc_Q_E[i].data["ERROR"], file_lc_Q_E[i].data["FRACEXP"]
	Q_E_time, Q_E_rate, Q_E_fracexp = file_lc_Q_E[1].data["TIME"], file_lc_Q_E[1].data["RATE"],  file_lc_Q_E[1].data["FRACEXP"]
	index_0 = [Q_E_rate!=0][0]
	Q_E_time = Q_E_time[index_0]
	Q_E_rate = Q_E_rate[index_0]
	#Q0_E1_error = Q0_E1_error[index_0]
	Q_E_fracexp = Q_E_fracexp[index_0]
	Q_E_rate = Q_E_rate/Q_E_fracexp #comment this line if live time correction is already done
	
	





	##################################Dividing into chunks####################################
	lc_len_Q_E = len(Q_E_rate)

	chunks = int((Q_E_time[len(Q_E_time)-1]-Q_E_time[0])/(binsize*1000))
	chunk_residual = int(Q_E_time[len(Q_E_time)-1]-Q_E_time[0])%(binsize*1000)
	b1 = np.empty(shape=[chunks+1])
	b2 = np.empty(shape=[chunks+1])
	sample_size = np.empty(shape=[chunks+1])
	if chunks == 0:
		b1[0] = Q_E_time[0]-binsize/2.0
		b2[0] = Q_E_time[len(Q_E_time)-1]+binsize/2.0
		
	for i in range(0,int(chunks)+1):
		if i == chunks:
			if chunk_residual != 0:
				b2[i] = Q_E_time[len(Q_E_time)-1]+binsize/2.0
				b1[i] = b2[i] - (binsize*1000) - binsize
				break
			else:
				b1 = np.delete(b1,len(b1)-1)
				b2 = np.delete(b1,len(b2)-1)
				
		b1[i] = Q_E_time[0]-binsize/2.0+i*(binsize*1000)
		b2[i] = b1[i]+ (binsize*1000) + binsize

	sample_size = (b2-b1-binsize)/binsize

	#current_index = [(Q_E_time>b1[1]) & (Q_E_time<b2[1])][0]
	###########################################################################################
	

	for i in range(0,len(b1)):
	
		current_index = [(Q_E_time>b1[i]) & (Q_E_time<b2[i])][0]
		sample_size = (b2[i]-b1[i]-binsize)/binsize
		current_rate = Q_E_rate[current_index]
		current_time = Q_E_time[current_index]
		#current_error = Q_E_error[current_index]
		
		if i==0:
			current_time = np.delete(current_time,0)
			current_rate = np.delete(current_rate,0)
		if i==len(b1)-1:
			current_time = np.delete(current_time,len(b1)-1)
			current_rate = np.delete(current_rate,len(b1)-1)
			
		tstart = current_time[0]
		tstop = current_time[len(current_time)-1]
		
		a1, c1 = np.polyfit(current_time,current_rate,1)
		rate_linear = a1*current_time+c1
		
		a2,c2,d2 = np.polyfit(current_time,current_rate,2)
		rate_quad = a2*current_time**2+c2*current_time+d2

		mean = np.mean(current_rate)
		rms = np.std(current_rate)
		rms_linear  = np.std(current_rate-rate_linear)
		rms_quad = np.std(current_rate-rate_quad)
		

		#plt.plot(current_time,current_rate-rate_quad)
		#plt.plot(current_time,rate_linear)
		#plt.plot(current_time,rate_quad)
		#plt.legend()
		#plt.show()
		
		
		#print rms
		#print rms_linear
		#print rms_quad
		
		cq_mean = np.sum((current_rate-mean)**2/mean)/(len(current_time)-1)
		cq_linear = np.sum((current_rate-rate_linear)**2/rate_linear)/(len(rate_linear)-1)
		cq_quad = np.sum((current_rate-rate_quad)**2/rate_quad)/(len(rate_quad)-1)
		
		#print cq_mean
		#print cq_linear 
		#print cq_quad
		#print str(rms)+"\t"+str(mean)+"\t"+str(sample_size[i])
		thresh_high, thresh_low = poissonian_limit(mean,rms,sample_size) 
		
		#print str(thresh_high)+"\t"+str(thresh_low)+"\t"+str((thresh_high-mean)/rms)+"\t"+str((mean-thresh_low)/rms)
		#print np.sum([(current_rate-rate_linear) < (thresh_low-mean)][0])
		#print np.sum([(current_rate-rate_linear) > (thresh_high-mean)][0])
		#print np.sum([(current_rate-rate_quad) < (thresh_low-mean)][0])
		#print np.sum([(current_rate-rate_quad) > (thresh_high-mean)][0])
		
		
		if cq_mean<cq_linear and cq_mean<cq_quad:
			high_index =[current_rate > thresh_high][0]
			low_index =[current_rate < thresh_low][0]
			time_high = current_time[high_index]
			time_low = current_time[low_index]
			counts_high = current_rate[high_index]
			counts_low = current_rate[low_index]
			
			np1 =  np.sum(high_index)
			np2 =  np.sum(low_index)
			if np1>0:
				for k in range(0,np1):
					peak_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(0)+"\t"+str(i)+"\t"+str(time_high[k])+"\t"+str(counts_high[k])+"\t"+str(float((counts_high[k]-mean)/rms))+"\n")
			if np2>0:
				for k in range(0,np2):
					dip_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(0)+"\t"+str(i)+"\t"+str(time_low[k])+"\t"+str(counts_low[k])+"\t"+str(float((mean-counts_low[k])/rms))+"\n")
			if np1>0 or np2>0:
				plt.plot(current_time,current_rate)
				plt.xlabel("Time")
				plt.ylabel("Counts/sec")
				#plt.label("Count rate = "+str(mean/binsize))
				plt.title(str(file_name)+"_"+str(i)+"   "+str(mean/binsize)+"  meanfit")
				plt.savefig(path_to_plots+str(file_name)+"_"+str(i)+"_plot.png")
				#plt.close(fig)
				plt.clf()
				
				corrected_current_time = current_time[(~np.add(high_index,low_index))]
				corrected_count_rate = current_rate[(~np.add(high_index,low_index))]
				corrected_mean = np.mean(corrected_count_rate)
				corrected_rms = np.std(corrected_count_rate)
				corrected_cq_mean = np.sum((corrected_count_rate-corrected_mean)**2/corrected_mean)/(len(corrected_current_time)-1)
			else:
				corrected_cq_mean = cq_mean
				
			orbit_report_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(i)+"\t"+str(len(b1))+"\t"+str(tstart)+"\t"+str(tstop)+"\t"+str(sample_size)+"\t"+"\t"+str(mean)+"\t"+str(rms)+"\t"+str(cq_mean)+"\t"+str(cq_linear)+"\t"+str(cq_quad)+"\t"+str(0)+"\t"+str(np1)+"\t"+str(np2)+"\t"+str(corrected_cq_mean)+"\n")
			
		elif cq_linear<cq_mean and cq_linear<cq_quad:
			high_index =[(current_rate-rate_linear) > (thresh_high-mean)][0]
			low_index =[(current_rate-rate_linear) < (thresh_low-mean)][0]
			time_high = current_time[high_index]
			time_low = current_time[low_index]
			counts_high = current_rate[high_index]
			counts_low = current_rate[low_index]
			fit_high = rate_linear[high_index]
			fit_low = rate_linear[low_index]
			
			
			
			
			np1 =  np.sum(high_index)
			np2 =  np.sum(low_index)
			if np1>0:
				for k in range(0,np1):
					peak_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(1)+"\t"+str(i)+"\t"+str(time_high[k])+"\t"+str(counts_high[k])+"\t"+str(float((counts_high[k]-fit_high[k])/rms_linear))+"\n")
		
			if np2>0:
				for k in range(0,np2):
					dip_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(1)+"\t"+str(i)+"\t"+str(time_low[k])+"\t"+str(counts_low[k])+"\t"+str(float((fit_low[k]-counts_high[k])/rms_linear))+"\n")
				
			if np1>0 or np2>0:
				plt.plot(current_time,current_rate-rate_linear)
				plt.xlabel("Time")
				plt.ylabel("Counts/sec")
			#	plt.label("Count rate = "+str(mean/binsize))
				plt.title(str(file_name)+"_"+str(i)+"   "+str(mean/binsize)+"  linearfit")
				plt.savefig(path_to_plots+str(file_name)+"_"+str(i)+"_plot.png")
				#plt.close(fig)
				plt.clf()
				corrected_current_time = current_time[(~np.add(high_index,low_index))]
				corrected_count_rate = current_rate[(~np.add(high_index,low_index))]
				a1, c1 = np.polyfit(corrected_current_time,corrected_count_rate,1)
				corrected_rate_linear = a1*corrected_current_time+c1
				corrected_rms_linear  = np.std(corrected_count_rate-corrected_rate_linear)
				corrected_cq_linear  = np.sum((corrected_count_rate-corrected_rate_linear)**2/corrected_rate_linear)/(len(corrected_rate_linear)-1)
			else:
				corrected_cq_linear = cq_linear
			orbit_report_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(i)+"\t"+str(len(b1))+"\t"+str(tstart)+"\t"+str(tstop)+"\t"+str(sample_size)+"\t"+"\t"+str(mean)+"\t"+str(rms)+"\t"+str(cq_mean)+"\t"+str(cq_linear)+"\t"+str(cq_quad)+"\t"+str(1)+"\t"+str(np1)+"\t"+str(np2)+"\t"+str(corrected_cq_linear)+"\n")
		
			
		elif cq_quad<cq_mean and cq_quad<cq_linear:
			high_index = [(current_rate-rate_quad) > (thresh_high-mean)][0]
			low_index = [(current_rate-rate_quad) < (thresh_low-mean)][0]
			time_high = current_time[high_index]
			time_low = current_time[low_index]
			np1 =  np.sum(high_index)
			np2 =  np.sum(low_index)
			counts_high = current_rate[high_index]
			counts_low = current_rate[low_index]
			fit_high = rate_quad[high_index]
			fit_low = rate_quad[low_index]
			
			
			
			if np1>0:
				for k in range(0,np1):
					peak_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(2)+"\t"+str(i)+"\t"+str(time_high[k])+"\t"+str(counts_high[k])+"\t"+str(float((counts_high[k]-fit_high[k])/rms_quad))+"\n")
		
			if np2>0:
				for k in range(0,np2):
					dip_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(2)+"\t"+str(i)+"\t"+str(time_low[k])+"\t"+str(counts_low[k])+"\t"+str(float((fit_low[k]-counts_low[k])/rms_quad))+"\n")
			if np1>0 or np2>0:
				plt.plot(current_time,current_rate-rate_quad)
				#plt.plot(current_time,current_rate)
				plt.xlabel("Time")
				plt.ylabel("Counts/sec")
				#plt.label("Count rate = "+str(mean/binsize))
				plt.title(str(file_name)+"_"+str(i)+"   "+str(mean/binsize)+"  quadfit")
				plt.savefig(path_to_plots+str(file_name)+"_"+str(i)+"_plot.png")
				#plt.close(fig)
				plt.clf()
				corrected_current_time = current_time[(~np.add(high_index,low_index))]
				corrected_count_rate = current_rate[(~np.add(high_index,low_index))]
				a2,c2,d2 = np.polyfit(corrected_current_time,corrected_count_rate,2)
				corrected_rate_quad = a2*corrected_current_time**2+c2*corrected_current_time+d2
				corrected_rms_quad  = np.std(corrected_count_rate-corrected_rate_quad)
				corrected_cq_quad  = np.sum((corrected_count_rate-corrected_rate_quad)**2/corrected_rate_quad)/(len(corrected_rate_quad)-1)
			else:
				corrected_cq_quad = cq_quad
				
				
	
			orbit_report_file.write(str(file_name)+"\t"+str(quad_no)+"\t"+str(energy_range)+"\t"+str(i)+"\t"+str(len(b1))+"\t"+str(tstart)+"\t"+str(tstop)+"\t"+str(sample_size)+"\t"+"\t"+str(mean)+"\t"+str(rms)+"\t"+str(cq_mean)+"\t"+str(cq_linear)+"\t"+str(cq_quad)+"\t"+str(2)+"\t"+str(np1)+"\t"+str(np2)+"\t"+str(corrected_cq_quad)+"\n")
		
