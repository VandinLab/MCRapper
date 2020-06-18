import math
import numpy as np
import sys
import os
import re
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-dbs","--datasets",type=int, default=2, help="0 = only light , 1 = only heavy , 2 = both")
args = parser.parse_args()

db_sizes = {"svmguide3" : 1243,
"chess" : 3196,
"mushroom" : 8124,
"phishing" : 11055,
"breast-cancer" : 12773,
"a9a" : 32561,
"pumb-star" : 49046,
"bms-web1" : 58136,
"connect" : 67557,
"bms-web2" : 77158,
"retail" : 88162,
"ijcnn1" : 91701,
"T10I4D100K" : 100000,
"T40I10D100K" : 100000,
"cod-rna" : 271617,
"accidents" : 340183,
"bms-pos" : 515597,
"covtype" : 581012,
"susy" : 5000000}

heavy_datasets = ["chess","connect","pumb-star","accidents","svmguide3","phishing"]
sup_heavy_datasets = ["connect","accidents","phishing","pumb-star","susy"]
sup_heavy_datasets2 = ["chess"]
tfp_sup_heavy = 0.05
tfp_sup_heavy2 = 0.05
tfp_f_sup_heavy = 0.0
tfp_t_sup_heavy = 0.05
tfp_thr = 0.05
light_datasets = ["T10I4D100K","a9a","bms-pos","bms-web2","mushroom","susy","T40I10D100K","bms-web1","breast-cancer","covtype","ijcnn1","retail"]
if args.datasets == 2:
    datasets = set(light_datasets) | set(heavy_datasets)
if args.datasets == 0:
    datasets = set(light_datasets) - set(heavy_datasets)
if args.datasets == 1:
    datasets = set(heavy_datasets) - set(light_datasets)

#datasets = ["a9a","mushroom","breast-cancer","retail","ijcnn1","svmguide3","cod-rna","covtype"]

#sample_sizes = [1000, 10000, 100000, 1000000]
sample_sizes = [10000.]#np.logspace(3. , 6. , 6)
small_wait = 1
medium_wait = 2
long_wait = 3
delta = 0.1

run_fisher = 1
parallel = 1

runs = 10

#ks = [1 , 10, 100]
ks = [10]

print(datasets)
time.sleep(medium_wait)

def run_experiment(db , sz , k , run_id , g_id):
    sz = int(sz)
    time.sleep(np.random.randint(small_wait))
    uid = str(args.datasets)+"_"+str(g_id)
    cmd = "python run_radest.py -db "+str(db)+" -sz "+str(sz)+" -jp "+str(k)+" -d "+str(delta)+" -run "+str(run_id)+" -uid "+str(uid)
    min_t = 0.0
    tfp_thr_ = tfp_thr
    if db in sup_heavy_datasets:
        tfp_thr_ = tfp_sup_heavy
        min_t = tfp_t_sup_heavy
    if db in sup_heavy_datasets2:
        tfp_thr_ = tfp_sup_heavy2
        min_t = tfp_t_sup_heavy
    cmd = cmd + " -tfp "+str(tfp_thr_)
    cmd = cmd + " -t "+str(min_t)
    print(cmd)
    os.system(cmd)
    time.sleep(np.random.randint(small_wait))

if parallel == 1:
    ids = list()
    from multiprocessing import Pool
    pool = Pool()
    g_id = 0
    for run_id in range(runs):
        for sz in sample_sizes:
            for db in datasets:
                for k in ks:
                    #sz = min(db_sizes[db] , sz)
                    if db_sizes[db] >= sz or 1==1:
                        ids.append( pool.apply_async(run_experiment , [db , sz , k , run_id , g_id]))
                        g_id = g_id + 1
    print("Launched "+str(g_id)+" experiments")
    todo_total = len(ids)
    done_ = 0
    for id in ids:
        id.get()
        done_ = done_ + 1
        print("DONE "+str(done_)+"/"+str(todo_total))
else:
    for sz in sample_sizes:
        for db in datasets:
            sz = db_sizes[db]
            for run_id in range(runs):
                run_experiment(db , sz , run_id)
