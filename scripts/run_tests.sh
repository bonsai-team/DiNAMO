#!/bin/bash
set -ev
if [ "$TRAVIS_BRANCH" = "unstable" ]; then
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
		echo "Error : executable could not be found"
	fi

	for DATASET in $DATASETS
	do
		echo "Testing ChIP-Seq mode on $DATASET with -l 6 -d 6..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 6 -d 6 -o temp.meme --no-log > /dev/null
		diff temp.meme $TEST_DIR/$DATASET/dinamo_6_6.meme
		if [ $? -ne 0 ]; then
        	echo "Test failed; the program did not output the expected file." >&2
			exit 1;
		else
			echo "Test passed, no difference found."
    	fi
			echo "Testing ChIP-Seq mode on $DATASET with -l 7 -d 3..."
		bin/dinamo -pf $TEST_DIR/$DATASET/$SIG_FILE -nf $TEST_DIR/$DATASET/$CTRL_FILE -l 7 -d 3 -o temp.meme --no-log > /dev/null
		diff temp.meme $TEST_DIR/$DATASET/dinamo_7_3.meme
		if [ $? -ne 0 ]; then
        	echo "Test failed; the program did not output the expected file." >&2
			exit 1;
		else
			echo "Test passed, no difference found."
    	fi
	done
	rm temp.meme
fi
