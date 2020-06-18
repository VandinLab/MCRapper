import math
import numpy as np
import sys
import os
import re
import time
import os.path


datasets = ['covtype','retail','mushroom','bms-web1','bms-web2','cod-rna',
        'breast-cancer','ijcnn1','a9a','svmguide3','T10I4D100K','chess','bms-pos',
        'T40I10D100K','phishing','connect','accidents','pumb-star','susy'];
run_light = 0
ks = [1000000 , 10, 100000 , 100, 10000 , 1000]
jp = 10000
runs = 10
alpha = 0.05
#small_wait = 10
#medium_wait = 60
#long_wait = 600
small_wait = 10
medium_wait = 60
long_wait = 600

def get_file_suffix(light , j , dataset , k , run_id):

    suffix_ = ""
    if light == 1:
        suffix_ = suffix_ + "light_"

    suffix_ = suffix_ + str(dataset)+"_j10e"+str(int(math.log10(j)))

    if k > 0:
        suffix_ = suffix_+"_k10e"+str(int(math.log10(k)))

    suffix_ = suffix_+"_r"+str(run_id)+".txt"

    return suffix_

def get_cmd(light , dataset , j , k, run_id):

    if light == 1:
        cmd_ = "./fim_closed_light"
    else:
        cmd_ = "./fim_closed"

    suff = get_file_suffix(light , j , dataset , k , run_id)
    suff_ = "results/" + str(suff)

    cmd_ = cmd_+ " "+str(suff_)

    cmd_ = cmd_+ " "+str(j)

    cmd_ = cmd_+ " "+str(alpha)

    cmd_ = cmd_+ " ../TopKWY/datasets/"+str(dataset)+"/"+str(dataset)+"_new.labels"
    cmd_ = cmd_+ " ../TopKWY/datasets/"+str(dataset)+"/"+str(dataset)+"_new.dat"

    cmd_ = cmd_+ " "+str(np.random.randint(1 , 999))

    if k > 0:
        cmd_ = cmd_+ " "+str(k)

    cmd_ = cmd_+ " > results/stats_"+str(suff)+" 2>&1 &"

    return cmd_


def check_results(light , dataset , j , k, run_id):

    if light == 1 and run_light == 0:
        return 1

    suff = get_file_suffix(light , j , dataset , k , run_id)
    stats_file_path = "results/stats_"+str(suff)
    if not os.path.exists(stats_file_path):
        return 0
    cmd = "tail -50 "+stats_file_path+" > temp_stats.txt"
    os.system(cmd)

    #store results
    results = -1
    running_time = 0.0
    wy_time = 0.0
    peak_memory = 0

    f_in = open("temp_stats.txt" , 'r')
    for line in f_in:

        #print line.replace('\n','')

        if light == 0:
            searchObj = re.match( r'Number of significant patterns: .*', line, re.M|re.I)
            if searchObj:
                temp_ = searchObj.group()
                temp_ = temp_[len("Number of significant patterns: "):]
                results = int(float(temp_))
                print str(dataset)+" finished run " + str(run_id) + " "+str(stats_file_path)
                print "  Results : " + str(results)

        searchObj = re.match( r'Runtime for correction .*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Runtime for correction (s): "):]
            running_time = float(temp_)
            if light == 1:
                print str(dataset)+" finished run " + str(run_id) + " "+str(stats_file_path)
                results = 0
            print "  Running time : " + str(running_time)

        searchObj = re.match( r'Runtime for WY computations .*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Runtime for WY computations (s): "):]
            wy_time = float(temp_)
            print "  Running time (WY) : " + str(wy_time)

        searchObj = re.match( r'Peak Memory.*', line, re.M|re.I)
        if searchObj:
            temp_ = searchObj.group()
            temp_ = temp_[len("Peak Memory (MB): "):]
            peak_memory = int(float(temp_))
            print "  Memory : " + str(peak_memory)

    if results >= 0:
        #print results for this run to file
        f_out = open("all_runs.csv" , 'a')
        f_out.write(str(light)+";"+str(k)+";"+str(j)+";"+str(dataset)+";"+str(run_id)+";"+str(results)+";"+str(running_time)+";"+str(wy_time)+";"+str(peak_memory)+";\n")
        return 1
    return 0


# main script

flag = 1
if len(sys.argv) > 1:
    flag = int(sys.argv[1])


dataset_completed_runs_light = dict()
for dataset in datasets:
    dataset_completed_runs_light[dataset] = 0

dataset_completed_runs_topk_inf = dict()
for dataset in datasets:
    dataset_completed_runs_topk_inf[dataset] = 0

dataset_completed_runs_topk = dict()
for dataset in datasets:
    temp = dict()
    for k in ks:
        temp[k] = 0
    dataset_completed_runs_topk[dataset] = temp

running = 1
running_datasets_light = set()
running_datasets_topk_inf = set()
running_datasets_topk = set()
running_datasets_topk_value = dict()

while running == 1:

    #print "new experiments "

    for dataset in datasets:
        if dataset not in running_datasets_light:
            running_datasets_light.add(dataset)
            # experiment for wy light
            cmd = get_cmd(1 , dataset , jp , 0 , dataset_completed_runs_light[dataset])
            if run_light == 1:
                os.system(cmd)
                print cmd
                time.sleep(np.random.randint(small_wait,2*small_wait))

        if dataset not in running_datasets_topk_inf:
            running_datasets_topk_inf.add(dataset)
            # experiment for topkwy with k=inf
            cmd = get_cmd(0 , dataset , jp , 0 , dataset_completed_runs_topk_inf[dataset])
            os.system(cmd)
            print cmd
            time.sleep(np.random.randint(small_wait,2*small_wait))

        min_run = runs
        min_run_k = -1
        if dataset not in running_datasets_topk:
            for k in ks:
                if dataset_completed_runs_topk[dataset][k] < min_run:
                    min_run_k = k
                    min_run = dataset_completed_runs_topk[dataset][k]

            if min_run_k != -1:
                running_datasets_topk.add(dataset)
                time.sleep(np.random.randint(small_wait,2*small_wait))
                # experiment for topkwy and k < inf
                k = min_run_k
                running_datasets_topk_value[dataset] = k
                cmd = get_cmd(0 , dataset , jp , k , dataset_completed_runs_topk[dataset][k])
                os.system(cmd)
                print cmd

    #print "checking "

    # check all datasets
    for dataset in datasets:

        # check result for wy light
        if dataset_completed_runs_light[dataset] < runs and check_results(1 , dataset , jp , 0 , dataset_completed_runs_light[dataset]) == 1:
            time.sleep(np.random.randint(small_wait,2*small_wait))
            running_datasets_light.discard(dataset)
            dataset_completed_runs_light[dataset] = dataset_completed_runs_light[dataset] + 1
        # check result for topkwy with k=inf
        if dataset_completed_runs_topk_inf[dataset] < runs and check_results(0 , dataset , jp , 0 , dataset_completed_runs_topk_inf[dataset]) == 1:
            time.sleep(np.random.randint(small_wait,2*small_wait))
            running_datasets_topk_inf.discard(dataset)
            dataset_completed_runs_topk_inf[dataset] = dataset_completed_runs_topk_inf[dataset] + 1

        # check result for topkwy with k < inf
        k = running_datasets_topk_value[dataset]
        if dataset in running_datasets_topk and check_results(0 , dataset , jp , k , dataset_completed_runs_topk[dataset][k]) == 1:
            time.sleep(np.random.randint(small_wait,2*small_wait))
            running_datasets_topk.discard(dataset)
            dataset_completed_runs_topk[dataset][k] = dataset_completed_runs_topk[dataset][k] + 1


    running = 0
    for dataset in datasets:
        if dataset_completed_runs_light[dataset] < runs or dataset_completed_runs_topk_inf[dataset] < runs:
            running = 1
        for k in ks:
            if dataset_completed_runs_topk[dataset][k] < runs:
                running = 1

    if running == 1:

        time.sleep(np.random.randint(medium_wait,2*medium_wait))
