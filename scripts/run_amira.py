import math
import os
import numpy as np
import sys
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db","--dataset", help="dataset name")
parser.add_argument("-sz","--samplesize", type=float ,help="sample size (>0).",default=10000.)
parser.add_argument("-g","--gamma", type=float ,help="confidence for Amira (in (0,alpha))",default=0.1)
parser.add_argument("-v","--verbose", help="increase output verbosity (def. false)")
parser.add_argument("-wd", help="path to the work directory ",default="work_dir/")
parser.add_argument("-run","--run",type=int ,help="run id",default=0)
parser.add_argument("-f",type=float ,help="assumed maximum frequency of the items",default=1.0)
parser.add_argument("-uid", help="unique id for experiments",default="0")
parser.add_argument("-keep", help="0= remove files for experiments",default="1")
args = parser.parse_args()

wait_time = 1
path_amira_folder = "../amira/"
datasets_path = "../datasets/"

def get_result(pattern , path ,  verbose=1):
    fin = open(path,'r')
    for line in fin:
        if pattern in line:
            line = line.replace('\n','')
            if verbose == 1:
                print(line)
            return line[len(pattern):]
    fin.close()

if not args.samplesize:
    print("Argument samplesize is needed")
    parser.print_help(sys.stderr)
    exit()
if not args.dataset:
    print("dataset name is needed!")
    parser.print_help(sys.stderr)
    exit()

print("dataset = "+str(args.dataset))
print("sample size = "+str(args.samplesize))
print("gamma = "+str(args.gamma))



def run_amira():
    global sample_size
    amira_path = path_amira_folder+"amira"
    work_dir_path = args.wd
    if not os.path.exists(work_dir_path):
        os.system("mkdir "+str(work_dir_path))
    temp_file_path = work_dir_path+"out_amira_"+str(args.uid)+".txt"

    trans_path = datasets_path+args.dataset+"/"+args.dataset+".dat"

    if not os.path.isfile(trans_path):
        print("path to dataset not correct! "+trans_path)
        exit()

    time.sleep(wait_time)

    # compute epsilon
    min_freq = 0.9
    sample_path = work_dir_path+"sample_amira_"+str(args.uid)+".dat"
    cmd = amira_path+" -f -n -p -s "+sample_path
    if args.f < 1.0:
        cmd = cmd+" -r "+str(args.f)
    cmd = cmd+" "+str(args.gamma)+" "+str(min_freq)+" "+str(args.samplesize)+" "+str(trans_path)+" > "+temp_file_path
    print("creating sample and computing epsilon for "+str(args.dataset)+"...")
    if args.verbose:
        print(cmd)
        time.sleep(wait_time)
    os.system(cmd)
    time.sleep(wait_time)

    #remove db with labels
    if args.keep == 0:
        os.system("rm "+sample_path)

    results_patterns = ("omega1: " , "rho1: ", "total: ","create_sample: ")
    results = list()
    for pattern in results_patterns:
        results.append(float(get_result(pattern , temp_file_path)))

    epsilon = results[1]
    rade = results[0]
    running_time = results[2]
    print(epsilon)
    fout = open(temp_file_path , 'w')
    i = 0
    for pattern in results_patterns:
        fout.write(pattern+str(results[i])+"\n")
        i = i + 1


run_amira()
