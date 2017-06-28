#!/bin/bash
set -ev
if [ "$TRAVIS_BRANCH" = "unstable" ]; then

	RED='\033[1;32m'
	GREEN='\033[1;32m'
	YELLOW='\033[1;33m'
	NC='\033[0m' # No Color

	CHECK_MARK='\u2714'
	CROSS='\u2716'

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

	for DATASET in $DATASETS
	do
		echo -e "Now testing with $DATASET dataset :"
		echo -e "\t Testing ChIP-Seq mode with parameters -l 6, -d 6..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 6 -d 6 -o temp.meme --no-log > /dev/null
		diff temp.meme $TEST_DIR/$DATASET/chipseq_6_6.meme
		if [ $? -ne 0 ]; then
        	echo -e "${RED}${CROSS}${NC}\t Test failed; the program did not output the expected file." >&2
			exit 1
    	fi

		echo -e "\t Testing ChIP-Seq mode with parameters -l 7, -d 3..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 7 -d 3 -o temp.meme --no-log > /dev/null
		diff temp.meme $TEST_DIR/$DATASET/chipseq_7_3.meme
		if [ $? -ne 0 ]; then
        	echo -e "${RED}${CROSS}${NC}\t Test failed; the program did not output the expected file." >&2
			exit 1
    	fi

		echo -e "\t Testing fixed position mode with parameters -l 5, -d 5, -p 1..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 5 -d 5 -p 1 --no-log > temp.res
		diff temp.res $TEST_DIR/$DATASET/position_5_5_1.results
		if [ $? -ne 0 ]; then
        	echo -e "${RED}${CROSS}${NC}\t Test failed; the program did not output the expected file." >&2
			exit 1
    	fi

		echo -e "\t Testing fixed position mode with parameters -l 7, -d 7, -p 4..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 7 -d 7 -p 4 --no-log > temp.res
		diff temp.res $TEST_DIR/$DATASET/position_7_7_4.results
		if [ $? -ne 0 ]; then
        	echo -e "${RED}${CROSS}${NC}\t Test failed; the program did not output the expected file." >&2
			exit 1
    	fi
		echo -e "${GREEN}$DATASET tests were successful !${NC}\n"
	done
	rm temp.meme
fi
