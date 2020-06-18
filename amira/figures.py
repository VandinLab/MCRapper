#! /usr/bin/env python3
#
# Copyright 2014,2018,2019 Matteo Riondato <riondato@acm.org>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import json
import math
import sys

import matplotlib.pyplot as plt
import numpy as np

def appendToDict(d):
    d['min'].append(1e26) # very big number
    d['max'].append(0.0)
    d['avg'].append(0.0)


def initializeDict(d):
    d['min'] = []
    d['max'] = []
    d['avg'] = []


def updateDict(d, v):
    if d['min'][-1] > v:
        d['min'][-1] = v
    if d['max'][-1] < v:
        d['max'][-1] = v
    d['avg'][-1] += v


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: {} varfile\n".format(sys.argv[0]))
        sys.exit(1)
    # "exec" the vars file to have all the variables in memory
    with open(sys.argv[1], 'rt') as varF:
        for line in varF:
            exec(line, globals())
    try:
        reps = int(REPETITIONS)
    except ValueError as e:
            sys.stderr.write(
                    "Error parsing the number of repetitions: {}\n".format(e))
            sys.exit(1)
    try:
        dbound = int(DBOUND)
    except ValueError as e:
            sys.stderr.write("Error parsing the dbound: {}\n".format(e))
            sys.exit(1)
    try:
        delta = float("0." + DELTAS)
    except ValueError as e:
            sys.stderr.write("Error parsing the delta: {}\n".format(e))
            sys.exit(1)
    logdelta = math.log(delta)
    try:
        theta = float("0." + THETAS)
    except ValueError as e:
            sys.stderr.write("Error parsing the theta: {}\n".format(e))
            sys.exit(1)
    try:
        exacttime = int(EXACT_TIME)
    except ValueError as e:
            sys.stderr.write("Error parsing the exact time: {}\n".format(e))
            sys.exit(1)
    sizesStr = SAMPLE_SIZES.split()
    sizes = []

    vceps = []
    eps = dict()
    initializeDict(eps)
    # For the errs, we compute the max of the max, the min of the min, the min
    # of the 1st quartile, the max of the 3rd quartile, the median of the
    # medians. Compute these quantities requires a special treatment.
    errs = []
    fps = dict()
    fps['fp'] = []
    fps['possfp'] = []
    fps['tp'] = []
    # The times dict contains only averages
    times = dict()
    times['sam'] = []
    times['rho1'] = []
    times['mine'] = []
    times['rho2'] = []
    times['prune'] = []
    tottimes = dict()
    initializeDict(tottimes)

    resFN_dir="../sampleres/"
    resFN_base = "_".join((DATASET, DELTAS, THETAS))
    resFN_suffix = "-res.json"

    for sizeStr in sizesStr:
        try:
            size = int(sizeStr)
        except ValueError as e:
            sys.stderr.write("Error parsing the size: {}\n".format(e))
            sys.exit(1)
        vceps.append(math.sqrt((dbound - logdelta) / size) * 2)
        sizes.append(size / 1000000)
        appendToDict(eps)
        errs.append([1.0, 1.0, 0.0, 0.0, 0.0])
        meds = []
        fps['fp'].append(0)
        fps['possfp'].append(0)
        fps['tp'].append(0)
        times['sam'].append(0.0)
        times['rho1'].append(0.0)
        times['mine'].append(0.0)
        times['rho2'].append(0.0)
        times['prune'].append(0.0)
        appendToDict(tottimes)
        tps = 0
        for rep in range(1, reps + 1):
            resFN = "_".join((resFN_base, sizeStr, str(rep))) + resFN_suffix
            with open(resFN_dir + resFN, 'tr') as resF:
                try:
                    res = json.load(resF)
                except json.JSONDecodeError as e:
                    sys.stderr.write("Error parsing the JSON file '{}': "
                            "{}\n".format(resFN ,e))
                    sys.exit(1)
            # XXX: should we also get rho1 and rho2?
            updateDict(eps, res['mine']['run']['eps'])
            meds.append(res['comp']['err_med'])
            if errs[-1][0] > res['comp']['err_min']:
                errs[-1][0] = res['comp']['err_min']
            if errs[-1][1] > res['comp']['err_1q']:
                errs[-1][1] = res['comp']['err_1q']
            if errs[-1][3] < res['comp']['err_3q']:
                errs[-1][3] = res['comp']['err_3q']
            if errs[-1][4] > res['comp']['err_max']:
                errs[-1][4] = res['comp']['err_max']
            fps['fp'][-1] += res['comp']['fp']
            fps['possfp'][-1] += res['comp']['possfp']
            fps['tp'][-1] += res['comp']['tp']
            if res['comp']['tp'] != res['comp']['orig']:
                sys.stderr.write(
                        "Error: missed some FIs in {}\n".format(resFN))
                sys.exit(1)
            if res['comp']['nafp'] != 0:
                sys.stderr.write(
                        "Error: got not-acceptable FPs in {}\n".format(resFN))
                sys.exit(1)
            if tps == 0:
                tps = res['comp']['tp']
            if tps != res['comp']['tp']:
                sys.stderr.write(
                        "Error: Unexpected number of TPs in {}\n".format(resFN))
                sys.exit(1)
            times['sam'][-1] += res['mine']['run']['runtimes']['create_sample']
            times['rho1'][-1] += res['mine']['run']['runtimes']['get_rho1']
            times['mine'][-1] += res['mine']['run']['runtimes']['mine']
            times['rho2'][-1] += res['mine']['run']['runtimes']['get_rho2']
            times['prune'][-1] += res['mine']['run']['runtimes']['prune']
            updateDict(tottimes, res['mine']['run']['runtimes']['total'])

        # Compute the averages from the sums
        eps['avg'][-1] /= reps
        meds.sort()
        errs[-1][2]=meds[reps // 2] # XXX: Assuming reps is odd.
        fps['fp'][-1] /= reps
        fps['possfp'][-1] /= reps
        fps['tp'][-1] /= reps
        times['sam'][-1] /= reps
        times['rho1'][-1] /= reps
        times['mine'][-1] /= reps
        times['rho2'][-1] /= reps
        times['prune'][-1] /= reps
        tottimes['avg'][-1] /= reps
    #We now have the statistics for all sample sizes
    sizesStr = [str(x) for x in sizes]
    pos = np.arange(1, len(sizes) + 1)

    #plt.rcParams['font.size'] = 12
    plt.rcParams['font.weight'] = 'bold'

    plt.figure(1)
    plt.plot(sizes, vceps, label='vc', ls='-', marker='s', lw=2.5)
    plt.plot(sizes, eps['max'], label='max', ls='-.', marker='v', lw=2.5)
    plt.plot(sizes, eps['avg'], label='avg', ls='-', marker='o', lw=2.5)
    plt.plot(sizes, eps['min'], label='min', ls=':', marker='^', lw=2.5)
    plt.ylabel(r'accuracy $\varepsilon$', fontweight='bold')
    plt.xlabel(r'sample size ($\times 10^6$)', fontweight='bold')
    plt.title(DATASET + r'--- $\theta={}$, $\delta={}$'.format(theta, delta),
            fontweight='bold')
    plt.gca().yaxis.grid(which="major", color='k', ls=':', lw=0.7)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.tight_layout()

    plt.figure(2)
    plt.plot(sizes, fps['fp'], label='FPs', ls='-.', marker='v', lw=2.5)
    plt.plot(sizes, fps['possfp'], label='Poss. FPs', ls='-', marker='o',
            lw=2.5)
    plt.plot(sizes, fps['tp'], label='TPs', ls=':', marker='^', lw=2.5)
    plt.ylabel('output itemsets', fontweight='bold')
    plt.xlabel(r'sample size ($\times 10^6$)', fontweight='bold')
    plt.title(DATASET + r'--- $\theta={}$, $\delta={}$'.format(theta, delta),
            fontweight='bold')
    plt.gca().yaxis.grid(which="major", color='k', ls=':', lw=0.7)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.tight_layout()

    plt.figure(3)
    boxLabels = list(x[0] + "\n(" + str(x[1]) + ")"
            for x in zip(sizesStr, eps['min']))
    plt.boxplot(errs, labels=boxLabels, whis='range')
    plt.ylabel('absolute error', fontweight='bold')
    plt.xlabel(r'sample size ($\times 10^6$)', fontweight='bold')
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.title(DATASET + r'--- $\theta={}$, $\delta={}$'.format(theta, delta),
            fontweight='bold')
    plt.tight_layout()

    # TODO: get runtimes for total and for VC.
    plt.figure(4)
    plt.plot(sizes, [exacttime] * len(sizes), label='exact', ls='-', marker='',
            lw=2.5, color='k')
    plt.plot(sizes, tottimes['max'], label='max', ls='-.', marker='v', lw=2.5)
    plt.plot(sizes, tottimes['avg'], label='avg', ls='-', marker='o', lw=2.5)
    plt.plot(sizes, tottimes['min'], label='min', ls=':', marker='^', lw=2.5)
    plt.ylabel('total runtime (ms)', fontweight='bold')
    plt.xlabel(r'sample size ($\times 10^6$)', fontweight='bold')
    plt.title(DATASET + r'--- $\theta={}$, $\delta={}$'.format(theta, delta),
            fontweight='bold')
    plt.gca().yaxis.grid(which="major", color='k', ls=':', lw=0.7)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend()
    plt.tight_layout()

    # Stacked barplot. See
    # https://python-graph-gallery.com/12-stacked-barplot-with-matplotlib/
    plt.figure(5)
    bottomMine = np.add(times['sam'], times['rho1']).tolist()
    bottomRho2 = np.add(bottomMine, times['mine']).tolist()
    bottomPrune = np.add(bottomRho2, times['rho2']).tolist()
    samBar = plt.bar(pos, times['sam'])
    rho1Bar = plt.bar(pos, times['rho1'], bottom=times['sam'], edgecolor='white')
    mineBar = plt.bar(pos, times['mine'], bottom=bottomMine, edgecolor='white')
    bars = [samBar, rho1Bar, mineBar]
    barsNames = ['sample', r'$\rho_1$', 'mine']
    if times['rho2'][0] != 0.0:
        rho2Bar = plt.bar(pos, times['rho2'], bottom=bottomRho2, edgecolor='white')
        pruneBar = plt.bar(pos, times['prune'], bottom=bottomPrune, edgecolor='white')
        bars += [rho2Bar, pruneBar]
        barNames += [r'$\rho_2', 'prune']
    plt.ylabel('total runtime (ms)', fontweight='bold')
    plt.xlabel(r'sample size ($\times 10^6$)', fontweight='bold')
    plt.xticks(pos, sizesStr)
    plt.title(DATASET + r'--- $\theta={}$, $\delta={}$'.format(theta, delta),
            fontweight='bold')
    plt.gca().yaxis.grid(which="major", color='k', ls=':', lw=0.7)
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.legend(bars, barsNames)
    plt.tight_layout()

    plt.show()

if __name__ == "__main__":
    main()
