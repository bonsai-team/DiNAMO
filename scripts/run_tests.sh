#!/bin/bash
set -ev

RED='\033[1;32m'
GREEN='\033[1;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

CROSS_MARK="\xe2\x9c\x97\x0a"
CHECK_MARK="\xe2\x9c\x93\x0a "

TEST_DIR='test'

DATASETS=''
DATASETS+='GATA1 '
DATASETS+='KLF1 '
DATASETS+='OCT4 '
DATASETS+='SOX2 '
DATASETS+='STAT3 '

CTRL_FILE='ctrl.fa'
SIG_FILE='peaks.fa'

if ! [ -f bin/dinamo ]
then
    echo -e "${RED}Error : executable could not be found${NC}"
    exit 1
fi

for DATASET in ${DATASETS}
do
    echo -e "Now testing with ${DATASET} dataset :"
    # ChIP-Seq
    for (( l=2; l<=5; l++ ))
    do
        for (( d=0; d<=$l; d++ ))
        do
            echo -en "\t Testing ChIP-Seq mode with parameters -l ${l}, -d ${d}..."
            bin/dinamo -pf ${TEST_DIR}/${DATASET}/${SIG_FILE} -nf ${TEST_DIR}/${DATASET}/${CTRL_FILE} -l ${l} -d ${d} -o  ${TEST_DIR}/${DATASET}/chipseq_l${l}_d${d}_test.meme 2> /dev/null 1> /dev/null
            diff ${TEST_DIR}/${DATASET}/chipseq_l${l}_d${d}_test.meme  ${TEST_DIR}/${DATASET}/chipseq_l${l}_d${d}.meme
            if [ $? -ne 0 ]; then
                echo -en "${RED}${CROSS_MARK}${NC}\t Test failed; the program did not output the expected file." >&2
                exit 1
            fi
            echo -en "\t${GREEN}${CHECK_MARK}${NC}"
            rm -f ${TEST_DIR}/${DATASET}/chipseq_l${l}_d${d}_test.meme
        done
    done
    # Fixed pos
    for (( l=2; l<=5; l++ ))
    do
        for (( d=0; d<=$l; d++ ))
        do
            for p in 1 2 3 4
            do
                echo -en "\t Testing fixed position mode with parameters -l ${l}, -d ${d}, -p ${p}..."
                bin/dinamo -pf ${TEST_DIR}/${DATASET}/${SIG_FILE}  -nf ${TEST_DIR}/${DATASET}/${SIG_FILE}  -l ${l} -d ${d} -p ${p} --no-log >  ${TEST_DIR}/${DATASET}/position_l${l}_d${d}_p${p}_test.res
                diff ${TEST_DIR}/${DATASET}/position_l${l}_d${d}_p${p}_test.res ${TEST_DIR}/${DATASET}/position_l${l}_d${d}_p${p}.res
                if [ $? -ne 0 ]; then
                    echo -en "${RED}${CROSS_MARK}${NC}\t Test failed; the program did not output the expected file." >&2
                    exit 1
                fi
                echo -en "\t${GREEN}${CHECK_MARK}${NC}"
                rm -f ${TEST_DIR}/${DATASET}/position_l${l}_d${d}_p${p}_test.res
            done
        done
    done
done
