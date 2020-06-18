#!/usr/bin/python
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
#datasets = ['mushroom','bms-web2','bms-web1', 'retail', 'T10I4D100K'];
datasets = ['susy','kosarac','accidents','pumb-star','chess','mushroom','connect', 'retail', 'bms-pos', 'bms-web1', 'bms-web2','T10I4D100K','T40I10D100K'];
#datasets =['connect','chess']#['bms-pos','bms-web2','accidents'];
#datasets = ['chess','bms-pos','bms-web1']
#datasets = ['accidents','pumb-star','connect','bms-pos','susy']
datasets = ['kosarac','susy']
#datasets = ['mushroom','bms-pos','bms-web2']
#datasets = ['mushroom']
datasets = ['susy','accidents','pumb-star','chess','connect', 'bms-pos', 'T40I10D100K','mushroom','retail', 'bms-web1', 'bms-web2','T10I4D100K','svmguide1','cod-rna','covtype','breast-cancer','ijcnn1','phishing','real-sim','rcv1','gisette','a1a','svmguide3','a9a'];
datasets = ['retail','mushroom','bms-web1','bms-web2','cod-rna','covtype','breast-cancer','ijcnn1','phishing','accidents','a9a','svmguide3','pumb-star','connect','bms-pos','T10I4D100K','T40I10D100K','susy','chess'];


jp = 10000;
alpha = 0.05;
maxram = 50000;
flag = "final"
runs = 10


datasets = ['retail','mushroom','bms-web1','bms-web2','cod-rna','covtype','breast-cancer','ijcnn1','a9a','svmguide3','T10I4D100K'];


times = dict()
memory = dict()
tested = dict()
runs_count = dict()

for i in datasets:
	times[i] = np.zeros(runs)
	memory[i] = np.zeros(runs)
	tested[i] = np.zeros(runs)
	runs_count[i] = 0

print times



def launch_experiment(dat_name , jp , alpha , flag , run_id):
	print "Run new experiment for " + dat_name
	print "jp = "+str(jp)+" alpha = "+str(alpha)+" flag = "+str(flag)
	out_file_name = 'experimentnew_'+str(dat_name)+'_'+str(jp)+"_run"+str(run_id)+"_"+flag;
	# launch the experiment
	seed = randint(0, 999)
	cmd = "cd datasets/" + dat_name + "/ && ./fim_closed " + out_file_name + " "+str(jp)+" "+str(alpha)+" " + dat_name + "_new.labels " + dat_name + "_new.dat " +str(seed)+" 2>&1 &"
	os.system(cmd)
	print cmd


# launch the first experiments
for dat_name in datasets:
	launch_experiment(dat_name , jp , alpha , flag , 0)
	time.sleep(randint(10, 20))

#running_datasets = datasets
running_datasets = list()
for i in datasets:
	running_datasets.append(i)

running = 1
while (running == 1):


	time.sleep( 10 )

	# check if experiment concluded and collect results
	for dat_name in running_datasets[:]:

		out_file_name = 'experimentnew_'+str(dat_name)+'_'+str(jp)+"_run"+str(runs_count[dat_name])+"_"+flag+"_timing.txt";
		outname = "datasets/" + dat_name + "/"+out_file_name
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
					time_txt = time_txt[22:(len(time_txt)-4)]
					time_value = -1.0
					try:
						time_value = float(time_txt)
					except ValueError:
						time_value = -1.0

					print "Running time: " + str(time_value)

					results = times[dat_name]
					results[runs_count[dat_name]] = time_value
					times[dat_name] = results

					running_datasets.remove(dat_name)
					found = 1

				searchObj = re.match( r'Peak memory consumption: .*', line, re.M|re.I)
				if searchObj:
					memory_txt = searchObj.group()
					memory_txt = memory_txt[25:]
					memory_txt = memory_txt[:6]

					memory_value = -1.0
					try:
						memory_value = float(memory_txt)
					except ValueError:
						memory_value = -1.0

					results = memory[dat_name]
					results[runs_count[dat_name]] = memory_value
					memory[dat_name] = results

					print "Memory usage: " + str(memory_value)

			else:
				#if found == 0:
				#	print dat_name + " still running"
				ok = 0


	# check if all experiments concluded
	# this happens when all datasets performed all the runs
	if len(datasets) == 0:
		running = 0

	# experiments which are not running and not have done all the runs, start a new run
	for i in datasets:
		if i not in running_datasets:
			runs_count[i] = runs_count[i] + 1
			if runs_count[i] < runs:
				time.sleep(randint(10, 20))
				launch_experiment(i , jp , alpha , flag , runs_count[i])
				running_datasets.append(i)

				# print results of this run on file
				file_out_log = open(flag+"_"+str(jp)+".txt",'w')
				file_out_log.write(str(jp)+";"+str(alpha)+";"+i+";"+str(times[i][runs_count[i]-1])+";"+str(memory[i][runs_count[i]-1])+";"+str(runs_count[i]-1)+";\n")
				file_out_log.close()

			else:
				print " all runs done "
				# all runs concluded for this dataset
				# compute average of metrics and variance and print it to file
				b = times[i]
				times[i] = b[b>0]
				mean_time = np.mean(times[i])
				variance_time = np.var(times[i])

				b = memory[i]
				memory[i] = b[b>0]
				mean_memory = np.mean(memory[i])
				variance_memory = np.var(memory[i])

				all_results_file = open("all_results.csv",'a')
				print i+" done with mean "+str(mean_time)+" and variance "+str(variance_time)
				print i+" done with mean "+str(mean_memory)+" and variance "+str(variance_memory)
				all_results_file.write(str(jp)+";"+str(alpha)+";"+i+";"+str(mean_time)+";"+str(variance_time)+";"+str(mean_memory)+";"+str(variance_memory)+";\n")
				all_results_file.close()
				datasets.remove(i)

f_in.close()
