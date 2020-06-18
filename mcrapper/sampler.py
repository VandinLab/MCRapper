# create random sample from dataset and launch rademacher estimator to compute the maximum epsilon
import sys
import os
import time

if len(sys.argv) < 7:
    print "error with input parameters; dataset, labels, word, sample size, jp, delta."
    exit()

dataset = sys.argv[1]
labels = sys.argv[2]
word = sys.argv[3]
sample_size = int(sys.argv[4])
jp = int(sys.argv[5])
delta = float(sys.argv[6])

fin = open(dataset,'r')

if sample_size < 1:
    print "error with k"
    exit()
if delta < 0 or delta > 1:
    print "error with delta"
    exit()


import math
import numpy as np

# count number of transactions
trans_num = 0
for line in fin:
    trans_num = trans_num + 1

random_sample_indexes = np.random.randint(0 , high = trans_num , size = sample_size)

# sort indexes
random_sample_indexes = sorted(random_sample_indexes)
#print random_sample_indexes

# reopen dataset and output file
fin.close()
fin = open(dataset,'r')
sample_path_labels = dataset+"_samplelabels.dat"
sample_path = dataset+"_sample.dat"
fout = open(sample_path,'w')
fout_labels = open(sample_path_labels,'w')

trans_id = -1
j = 0
written = 0
for line in fin:
    trans_id = trans_id + 1
    while j < len(random_sample_indexes) and random_sample_indexes[j] < trans_id:
        j = j + 1
    while j < len(random_sample_indexes) and random_sample_indexes[j] == trans_id:
        j = j + 1
        fout.write(line)
        fout_labels.write(str(np.random.randint(0,2))+"\n")
        written = written + 1

print "printed "+str(written)+" transactions"

fout_labels.close()
fout.close()

rand_seed = np.random.randint(0 , 1000)
cmd = "./fim_closed test_sampler_"+word+" "+str(jp)+" 0.05 "+sample_path_labels+" "+sample_path+" "+str(rand_seed)
print cmd
time.sleep(1)
os.system(cmd)
