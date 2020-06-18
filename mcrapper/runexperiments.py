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
datasets = ['kosarac','pumb-star','chess','mushroom','connect', 'retail', 'bms-pos', 'bms-web1', 'bms-web2','T10I4D100K','T40I10D100K'];
datasets = ['susy']
datasets = ['mushroom','retail','bms-web1', 'bms-web2','T10I4D100K','pumb-star','chess','mushroom','connect', 'retail', 'bms-pos', 'bms-web1', 'bms-web2','T10I4D100K','T40I10D100K','susy']

datasets = ['cod-rna'];

jp = 10000;
alpha = 0.01;
flag = ""

if len(sys.argv) > 1:
	flag = sys.argv[1]


for dat_name in datasets:
	out_file_name = 'experimentnew_'+str(dat_name)+'_'+str(jp)+flag;
	# launch the experiment
	seed = randint(0, 999)
	cmd = "cd datasets/" + dat_name + "/ && ./fim_closed " + out_file_name + " "+str(jp)+" "+str(alpha)+" " + dat_name + ".labels " + dat_name + ".dat " +str(seed)+" 2>&1 &"
	os.system(cmd)
	print cmd
	time.sleep( 60 )

running_datasets = datasets


while (len(running_datasets) > 0):

	time.sleep( 100 )
	print len(running_datasets)
	# check if experiment concluded and collect results
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
					time_value = float(time_txt)
					print "Running time: " + str(time_value)
					running_datasets.remove(dat_name)
					found = 1

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
					memory_txt = memory_txt[25:]
					memory_txt = memory_txt[:6]
					memory_value = float(memory_txt)
					print "Memory usage: " + str(memory_value)

			else:
				#if found == 0:
				#	print dat_name + " still running"
				ok = 0

print "done!"
f_in.close()
time.sleep(180)
