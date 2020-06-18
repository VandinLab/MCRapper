#! /usr/bin/env bash
#
# Run an experiment or collect the results
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

#TODO 20190206: Update code to collect the results.

DATA_DIR=`readlink -f ../datasets`
LOG_DIR=`readlink -f ../logs`
SAMPLE_RES_DIR=`readlink -f ../sampleres`

if [ $# -lt 2 ]; then
    echo "$0: run an experiment or collect the results" >&2
    echo "USAGE: $0 {run | collect {key1[,key2[,...]]}} variables_file" >&2
    exit 1
fi

COMMAND="$1"

if [ "${COMMAND}" = "run" ]; then
    if [ $# -ne 2 ]; then
        echo "Error: wrong number of arguments" >&2
        echo "USAGE: $0 {run | collect {key1[,key2[,...]]}} variables_file" >&2
        exit 1
    fi
    VARIABLES_FILE=`readlink -f "$2"`
    if [ ! -r "${VARIABLES_FILE}" ]; then
        echo "Error: Variables file '$2' not readable" >&2
        exit 1
    fi
elif [ "${COMMAND}" = "collect" ]; then
    if [ $# -ne 3 ]; then
        echo "USAGE: $0 {run | collect {key1[,key2[,...]]}} variables_file" >&2
        exit 1
    fi

    COLLECT_KEYS="$2"
    VARIABLES_FILE=`readlink -f $3`
    if [ ! -r "${VARIABLES_FILE}" ]; then
        echo "Error: Variables file '$3' not readable" >&2
        exit 1

    fi
else
    echo "Error: invalid or no command specified" >&2
    echo "USAGE: $0 {run | collect {key1[,key2[,...]]}} variables_file" >&2
    exit 1
fi

. "${VARIABLES_FILE}"

DATASET_FILE="${DATA_DIR}/${DATASET}.dat"

if [ ${COMMAND} = "run" -a ! -r "${DATASET_FILE}" ]; then
    echo "Error: Dataset '${DATASET_FILE}' not readable" >&2
    exit 1
fi

if [ ${COMMAND} = "run" -a ! -r "${EXACT}" ]; then
    echo "Error: Exact results file '${EXACT}' not readable" >&2
    exit 1
fi

for DELTA in ${DELTAS}; do
    for THETA in ${THETAS}; do
        SAMPLE_BASE="${DATASET}_${DELTA}_${THETA}"
        if [ ${COMMAND} = "run" ]; then
            for SIZE in ${SAMPLE_SIZES}; do
                for REPETITION in $(seq 1 $REPETITIONS); do
                    SAMPLE="${SAMPLE_BASE}_${SIZE}_${REPETITION}"
                    echo ${SAMPLE} && \
                    ./amira -j -f -v -d ${DATASET_SIZE} ${ADDIT_FLAGS} 0.${DELTA} 0.${THETA} ${SIZE} "${DATASET_FILE}" > "${SAMPLE_RES_DIR}/${SAMPLE}-mine.json" 2> "${LOG_DIR}/${SAMPLE}.log" && \
                    printf "{\n\"mine\":" > "${SAMPLE_RES_DIR}/${SAMPLE}-res.json"
                    # 32 and 3 are a magic numbers depending on the
                    # output of amira.
                    head -32 "${SAMPLE_RES_DIR}/${SAMPLE}-mine.json" >> "${SAMPLE_RES_DIR}/${SAMPLE}-res.json" && \
                    tail -3 "${SAMPLE_RES_DIR}/${SAMPLE}-mine.json" >> "${SAMPLE_RES_DIR}/${SAMPLE}-res.json" && \
                    # We do really need to sort them.
                    ./sort_fis -j -v "${SAMPLE_RES_DIR}/${SAMPLE}-mine.json" > "${SAMPLE_RES_DIR}/${SAMPLE}-sort.json" 2>> "${LOG_DIR}/${SAMPLE}.log" && \
                    rm "${SAMPLE_RES_DIR}/${SAMPLE}-mine.json" && \
                    ./compare_fis -jJo -v "${EXACT}" "${SAMPLE_RES_DIR}/${SAMPLE}-sort.json" > "${SAMPLE_RES_DIR}/${SAMPLE}-comp.json" 2>> "${LOG_DIR}/${SAMPLE}.log" && \
                    printf ",\n\"comp\": " >> "${SAMPLE_RES_DIR}/${SAMPLE}-res.json" && \
                    cat "${SAMPLE_RES_DIR}/${SAMPLE}-comp.json" >> "${SAMPLE_RES_DIR}/${SAMPLE}-res.json"
                    printf "}\n" >> "${SAMPLE_RES_DIR}/${SAMPLE}-res.json"
                    sudo umount / 2> /dev/null
                done
            done
        elif [ ${COMMAND} = "collect" ]; then
            python3 collect_results.py ${COLLECT_KEYS} "${SAMPLE_RES_DIR}"/${SAMPLE_BASE}_*-res.json > "${SAMPLE_RES_DIR}"/${SAMPLE_BASE}_${THETA}-collect.csv
        else # unreached
            echo "ERROR: you reached code supposed to be unreachable. How?" >&2
            exit 1
        fi
    done # for THETA
done # for DELTA
