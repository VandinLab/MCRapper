import math
import os
import numpy as np
import sys
import time
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-db","--dataset", help="dataset name")
parser.add_argument("-sz","--samplesize", type=float ,help="sample size (>0).",default=10000.)
parser.add_argument("-d","--delta", type=float ,help="confidence delta",default=0.1)
parser.add_argument("-j", type=int ,help="number of draws of Rademacher r.v.s",default=1)
parser.add_argument("-e","--eta", type=float ,help="Hybrid algorithm: confidence eta",default=0.01)
parser.add_argument("-t", type=float ,help="Hybrid algorithm: minimum frequency of patterns to enumerate.",default=0.0)
parser.add_argument("-ub", type=float ,help="Hybrid algorithm: union bound term (don't needed as it is automatically computed)",default=0)
parser.add_argument("-tfp", type=float ,help="theta parameter for TFP algorithm.",default=1.0)
parser.add_argument("-r","--results",help="path where to write results",default="results_radest.csv")
parser.add_argument("-run","--run",type=int ,help="run id",default=0)
parser.add_argument("-k", type=int ,help="enumerate only the top-k frequent patterns. (currently disabled)",default=-1)
parser.add_argument("-f", type=float ,help="maximum frequency of patterns to process.",default=1.0)
parser.add_argument("-fa", type=float ,help="maximum assumed frequency of items for amira.",default=1.0)
parser.add_argument("-var", type=float ,help="maximum assumed variance of patterns.",default=0.25)
parser.add_argument("-v","--verbose", help="increase output verbosity (1 to enable, def. 0)")
parser.add_argument("-uid", help="unique id for experiments",default="0")
parser.add_argument("-keep", help="0= remove files for experiments",default="1")
parser.add_argument("-bfs", type=int ,help="=1 use bfs search instead of dfs",default=0)
args = parser.parse_args()

wait_time = 1
run_only_amira = 0

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
print("delta = "+str(args.delta))
print("k = "+str(args.k))

def checkinput_transactions(path):
    fin = open(path , 'r')
    t_id = -1
    for line in fin:
        line = line.replace('\n','')
        t_id = t_id + 1
        items = line.split(' ')
        trans_items = set()
        for item in items:
            if len(item) > 0:
                item = int(item)
                if item not in trans_items:
                    trans_items.add(item)
                else:
                    print("Error! item "+str(item)+" is at least doubled in transaction "+str(t_id))
                    print("Fix it with fix_db.py")
                    exit()
    print("checked "+path)

seed_ = np.random.randint(99999)

def compute_correction_factor_hybrid(trans_path):
    debug_corr_factor = 0
    corr_factor = 0.
    items_freq = dict()
    fin = open(trans_path , 'r')
    t_id = -1
    inv_sz = 1./float(args.samplesize)

    max_length = 0
    num_transactions_per_length = np.zeros(101)
    # first pass on sample to compute the frequencies of the items, and to count
    # the lengths of the transactions
    for line in fin:
        line = line.replace('\n','')
        t_id = t_id + 1
        items = line.split(' ')
        trans_items = set()
        for item in items:
            if len(item) > 0:
                item = int(item)
                trans_items.add(item)
                if item not in items_freq:
                    items_freq[item] = inv_sz
                else:
                    items_freq[item] = items_freq[item]+inv_sz
        max_length = max(max_length , len(trans_items))
        if num_transactions_per_length.shape[0] < len(trans_items) + 1:
            temp_array = np.zeros((len(trans_items) + 1)*2)
            for i in range(num_transactions_per_length.shape[0]):
                temp_array[i] = num_transactions_per_length[i]
            num_transactions_per_length = temp_array.copy()
        num_transactions_per_length[len(trans_items)] += 1
    # compute first naive correction factor
    corr_factor += max_length * math.log(2.)
    sum_temp = 0.
    for i in range(1,max_length+1):
        sum_temp += num_transactions_per_length[i]*2**(i-max_length)
        if debug_corr_factor == 1:
            if num_transactions_per_length[i] > 0:
                print(" i ",i," num_transactions_per_length[i] ",num_transactions_per_length[i])
    corr_factor += math.log(sum_temp)
    corr_factor += math.log(1./args.eta)
    fin.close()
    if debug_corr_factor == 1:
        #print(num_transactions_per_length)
        print(items_freq)
        print("checked "+trans_path)

    print("corr_factor ",corr_factor)
    print("max_length ",max_length)
    return corr_factor

def run_radest(run_amira):
    global sample_size
    correct_path = "../mcrapper/fim_closed"
    if args.bfs == 1:
        correct_path = "../TopKWY/src/topkwy"
    work_dir_path = "work_dir/"
    if run_only_amira == 1:
        work_dir_path = "work_dir_amira_only/"
    if args.tfp < 1.0:
        work_dir_path = "work_dir_tfp/"
    if not os.path.exists(work_dir_path):
        os.system("mkdir "+str(work_dir_path))
    temp_file_path = work_dir_path+"out_"+str(args.uid)

    # run amira to create the sample
    cmd = "python run_amira.py -db "+str(args.dataset)+" -sz "+str(args.samplesize)+" -g "+str(args.delta)+" -wd "+str(work_dir_path)+" -uid "+str(args.uid)
    if args.fa < 1.0:
        cmd = cmd+" -f "+str(args.fa)
    if args.verbose == 1:
        cmd = cmd+" -v "+str(args.verbose)
    if run_amira == 1:
        print(cmd)
        os.system(cmd)
        time.sleep(wait_time)

    results_patterns = ("omega1: " , "rho1: ", "total: ","create_sample: ")
    amira_out_path = work_dir_path+"out_amira_"+str(args.uid)+".txt"
    results = list()
    for pattern in results_patterns:
        results.append(float(get_result(pattern , amira_out_path)))

    epsilon_amira = results[1]
    running_time_amira = (results[2]-results[3]) / 1000.
    #rade_amira = results[3]
    rade_amira = results[0]

    trans_path = work_dir_path+"sample_amira_"+str(args.uid)+".dat"

    checkinput_transactions(trans_path)
    time.sleep(wait_time)

    if not os.path.isfile(trans_path):
        print("path to dataset not correct! "+trans_path)
        exit()

    # correct
    out_path = temp_file_path+"_output.txt"
    if args.bfs == 0:
        if args.t > 0.:
            # reduce theta since we are using the hybrid algorighm
            delta_ = args.delta - args.eta
        else:
            delta_ = args.delta
        cmd = correct_path+" -b "+temp_file_path+" -j "+str(args.j)+" -d "+str(delta_)+" -i "+str(trans_path)
        cmd = cmd+" -s "+str(seed_)
        if args.k > 0:
            cmd = cmd+" -k "+str(args.k)
        if args.t > 0.:
            cmd = cmd+" -t "+str(args.t)
            args.delta = args.delta - args.eta
            if args.ub == 0.:
                # we have to compute the correction factor
                args.ub = compute_correction_factor_hybrid(trans_path)
            # otherwise, use the given one
            if args.ub > 0.:
                cmd = cmd+" -u "+str(args.ub)
        if args.var < 0.25:
            cmd = cmd+" -v "+str(args.var)
        if args.f < 1.0:
            cmd = cmd+" -f "+str(args.f)
        cmd = cmd+" > "+out_path+" 2>&1"

    if args.bfs == 1:
        trans_path_nodat = trans_path.replace(".dat","")
        trans_path_input = trans_path_nodat+"_new.spec"
        cmd = "python ../TopKWY/scripts/compute_spec_file.py "+trans_path_nodat
        os.system(cmd)
        cmd = correct_path+" -s "+str(trans_path_input)
        if args.t > 0.:
            cmd = cmd+" -t "+str(args.t)
        #cmd = cmd+" -d "+str(args.delta)
        cmd = cmd+" > "+out_path+" 2>&1"

    print("correcting "+str(args.dataset)+"...")
    if args.verbose:
        print(cmd)
        time.sleep(wait_time)
    if run_only_amira == 0:
        os.system(cmd)
        time.sleep(wait_time)

    results_patterns = ["Epsilon (n=1) " , "Epsilon (SBF bound) ", "Estimated Rademacher: " , "Runtime for correction (s): " , "Number of explored itemsets: ", "Epsilon (Var bound) ", "Epsilon (Hyb bound) ", "Epsilon (Hyb bound n=1) "]

    res = list()
    if run_only_amira == 0:
        for pattern in results_patterns:
            res.append(float(get_result(pattern , out_path)))
    else:
        for pattern in results_patterns:
            res.append(0.0)

    epsilon = res[0]
    epsilon_sbf = res[1]
    epsilon_kera = res[2]
    rade_est = res[3]
    time_to_correct = res[4]
    explored_patterns = res[5]
    epsilon_var = res[6]

    res = list()
    res.append(rade_est)
    res.append(epsilon)
    res.append(epsilon_sbf)
    res.append(epsilon_kera)
    res.append(rade_amira)
    res.append(epsilon_amira)
    res.append(time_to_correct)
    res.append(running_time_amira)
    res.append(explored_patterns)
    res.append(epsilon_var)
    res.append(args.tfp)
    res.append(args.t)

    if run_only_amira == 1:
        args.results = "results_radest_onlyamira.csv"
    if args.tfp < 1.0:
        args.results = "results_radest_tfp.csv"
    if args.bfs == 1:
        args.results = "results_radest_bfs.csv"

    print(args.results)
    fout = open(args.results , 'a')
    out = args.dataset+";"+str(args.samplesize)+";"+str(args.delta)+";"+str(args.j)+";"+str(seed_)+";"
    for res_ in res:
        out = out+str(res_)+";"
    out = out+"\n"
    fout.write(out)

    #remove amira sample
    if args.keep == 0:
        os.system("rm "+trans_path)
        os.system("rm "+labels_path)
        os.system("rm "+out_path)
        os.system("rm "+amira_out_path)

    if run_only_amira == 1:
        return epsilon_amira

    return min(epsilon_sbf , epsilon_var)

if args.tfp == 1.0:
    run_radest(1)
else:
    i = 0
    # don't remove intermediate files
    args.keep = 1
    max_iterations = 10
    eps = list()
    eps.append(1.0)
    epsilon_ = 1.0
    deltas = list()
    delta = args.delta
    delta_i = args.delta#/2.
    deltas.append(delta_i)
    # create sample and do first iteration of algorithm
    args.delta = delta_i
    args.var = min(0.25 , args.tfp*(1.-args.tfp))
    #frequent_only = 0
    #if args.t > 0.0:
    #    args.fa = args.t
    #    frequent_only = 1
    # compute epsilon with amira, that is later on used by radest to prune patterns
    run_only_amira = 1
    args.delta = delta_i
    eps_i = run_radest(1)
    args.f = 0. # only count number of true frequent patterns according to amira
    temp_t = args.t
    args.t = args.tfp + eps_i
    run_only_amira = 0
    run_radest(0)
    args.delta = delta
    args.t = temp_t
    #run_only_amira = 0
    # iterations of the algorithm
    min_eps = eps_i
    while(eps_i <= eps[len(eps)-1] and i <= max_iterations):
        i = i+1
        args.f = args.tfp + eps_i
        #args.t = max(args.t , args.tfp - eps_i)
        delta_i = delta_i#/2.
        eps.append(eps_i)
        args.delta = delta_i
        eps_i = run_radest(0)
        min_eps = min(min_eps , eps_i)

    args.f = 0. # only count number of true frequent patterns according to radest
    args.t = args.tfp + eps_i
    run_radest(0)
