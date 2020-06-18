import math
import numpy as np
import sys
import os
import re
import time
from random import randint

#dataset_name = sys.argv[1]
#f_in = open(dataset_name, 'r')

#datasets = ['mushroom', 'connect','retail', 'bms-pos', 'accidents', 'bms-web1', 'bms-web2','chess','pumb-star','T10I4D100K','T40I10D100K'];
#datasets = [ 'retail', 'T10I4D100K','bms-pos'];
datasets = ['kosarac','chess','mushroom','connect', 'retail', 'bms-pos', 'bms-web1', 'bms-web2','T10I4D100K','T40I10D100K'];

results_time = open('results_time.csv', 'a')
results_space = open('results_space.csv', 'a')

jp = 10000;
alpha = 0.05;
running_datasets = datasets
flag=""

if len(sys.argv) > 1:
	jp = sys.argv[1]
if len(sys.argv) > 2:
	flag = sys.argv[2]

while (len(running_datasets) > 0):

	# check if experiment concluded and collect results
	removed_something = 0
	for dat_name in running_datasets[:]:

		out_file_name = 'experimentnew_'+str(dat_name)+'_'+str(jp)+flag;
		outname = "datasets/" + dat_name + "/"+out_file_name+"_timing.txt"
		cmd = "tail -250 "+outname+" > tail_script.txt"
		os.system(cmd)
		#print cmd

		f_in = open('tail_script.txt', 'r')

		ok = 1
		found = 0
		while (ok==1) :
			line=f_in.readline()
			if len(line) > 0:

				searchObj = re.match( r'Total execution time: .*', line, re.M|re.I)
				if searchObj:
					time_txt = searchObj.group()
					time_txt = time_txt[22:]
					time_txt = time_txt[:(len(time_txt)-5)]
					time_value = float(time_txt)
					print dat_name+" - Running time: " + str(time_value)
					running_datasets.remove(dat_name)
					removed_something = 1
					found = 1
					results_time.write(dat_name+";"+str(jp)+";"+flag+";"+str(time_value)+";\n")

				#searchObj = re.match( r'Number of significant itemsets: .*', line, re.M|re.I)
				#if searchObj:
				#	res_txt = searchObj.group()
				#	res_txt = res_txt[32:]
				#	res_value = int(res_txt)
				#	print dat_name + " finished"
				#	print "Results: " + str(res_value)

				searchObj = re.match( r'Peak memory consumption: .*', line, re.M|re.I)
				if searchObj:
					memory_txt = searchObj.group()
					print memory_txt
					memory_txt = memory_txt[25:]
					#print memory_txt
					memory_txt = memory_txt[:len(memory_txt)-29]
					#print memory_txt
					memory_value = float(memory_txt)
					print dat_name+" - Memory usage: " + str(memory_value)
					results_space.write(dat_name+";"+str(jp)+";"+flag+";"+str(memory_value)+";\n")

			else:
				#if found == 0:
				#	print dat_name + " still running"
				ok = 0


	if removed_something == 1:
		removed_something = 0
		print running_datasets
	time.sleep( 10 )


f_in.close()
results_space.close()
results_time.close()
