# Collect results of an experiment by computing various statistics for the
# desired keywords in the full log
#
# Copyright 2014,2019 Matteo Riondato <matteo@cs.brown.edu>
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

import os.path
import statistics
import sys


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(
            os.path.basename(sys.argv[0]) +
            ": collect results of an experiment by computing various" +
            " statistics for the desired keywords in the logs\n")
        sys.stderr.write(
            "USAGE: {} key1[,key2[,...]] log1 [log2 [...]]\n".format(
                os.path.basename(sys.argv[0])))
        return 1

    nonstat_keys = [
        "dataset", "epsilon", "delta", "gamma", "topk",
        "stop_cond_standard_sat", "stop_cond_union_sat",
        "stop_cond_shatter_sat", "initial_sample_size",
        "geometric_schedule", "linear_schedule",
        "geometric_schedule_multiplier", "linear_schedule_addend",
        "use_lengths_CI_method", "use_items_CI_method",
        "use_size_CI_method", "use_eVC_CI_method",
        "use_refined_CI_method", "use_standard_stop_cond",
        "use_new_stop_cond", "use_refined_stop_cond",
        "use_optimize_stop_cond", "theta", "min_freq", "large_res",
        "sample_res", "orig_FIs", "d_index", "opt_param"]

    stat_keys = [
        "sample_size_avg", "sample_size_max", "sample_size_min",
        "sample_size_stdev", "iteration_avg", "iteration_max", "iteration_min",
        "iteration_stdev", "sampling_time_avg", "sampling_time_max",
        "sampling_time_min", "sampling_time_stdev", "stopcond_time_avg",
        "stopcond_time_max", "stopcond_time_min", "stopcond_time_stdev",
        "runtime_avg", "runtime_max", "runtime_min", "runtime_stdev",
        "mining_time_avg", "mining_time_max", "mining_time_min",
        "mining_time_stdev", "d_index_time_avg", "d_index_time_max",
        "d_index_time_min", "d_index_time_stdev", "total_time_avg",
        "total_time_max", "total_time_min", "total_time_stdev",
        "intersection_avg", "intersection_max", "intersection_min",
        "intersection_stdev", "false_negs_avg", "false_negs_max",
        "false_negs_min", "false_negs_stdev", "false_pos_avg", "false_pos_max",
        "false_pos_min", "false_pos_stdev", "non_acceptable_false_pos_avg",
        "non_acceptable_false_pos_max", "non_acceptable_false_pos_min",
        "non_acceptable_false_pos_stdev", "jaccard_avg", "jaccard_max",
        "jaccard_min", "jaccard_stdev", "wrong_eps_avg", "wrong_eps_max",
        "wrong_eps_min", "wrong_eps_stdev", "max_abs_err_avg",
        "max_abs_err_max", "max_abs_err_min", "max_abs_err_stdev",
        "avg_abs_err_avg", "avg_abs_err_max", "avg_abs_err_min",
        "avg_abs_err_stdev", "avg_rel_err_avg", "avg_rel_err_max",
        "avg_rel_err_min", "avg_rel_err_stdev"]

    stat_keys_prefixes = [
        "sample_size", "iteration", "sampling_time", "stopcond_time",
        "runtime", "mining_time", "d_index_time", "total_time", "intersection",
        "false_negs", "false_pos", "non_acceptable_false_pos", "jaccard",
        "wrong_eps", "max_abs_err", "avg_abs_err", "avg_rel_err",
        "opt_minimizer"]

    # Validate input keys
    input_keys = sys.argv[1].split(",")
    for key in input_keys:
        if key not in nonstat_keys and key not in stat_keys:
            sys.stderr.write("key '{}' not recognized.\n".format(key))
            return 1

    values = dict()
    for key in nonstat_keys + stat_keys_prefixes:
        values[key] = []

    for log_file in sys.argv[2:]:
        with open(log_file, "rt") as log_FILE:
            keys = log_FILE.readline().strip().split(",")
            log_values = log_FILE.readline().strip().split(",")
            file_dict = dict(zip(keys, log_values))
            for key in file_dict:
                if key in stat_keys_prefixes:
                    values[key].append(float(file_dict[key]))
                else:
                    values[key].append(file_dict[key])

    values_to_print = []
    for key in input_keys:
        if key in nonstat_keys:
            if len(values[key]) > 0:
                values_to_print.append(values[key][0])
            else:
                sys.stderr.write("Key {} not found in results, using 'not found' as value.\n".format(key))
                values_to_print.append("not found")
        else:
            last_underscore_index = key.rindex("_")
            prefix = key[:last_underscore_index]
            suffix = key[last_underscore_index + 1:]
            if suffix == "avg":
                values_to_print.append(str(statistics.mean(values[prefix])))
            elif suffix == "max":
                values_to_print.append(str(max(values[prefix])))
            elif suffix == "min":
                values_to_print.append(str(min(values[prefix])))
            elif suffix == "stdev":
                values_to_print.append(str(statistics.stdev(values[prefix])))

    print(",".join(input_keys))
    print(",".join(values_to_print))

    return 0

if __name__ == "__main__":
    main()
