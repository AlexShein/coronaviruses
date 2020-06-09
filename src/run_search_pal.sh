#!/bin/bash

# First arg: target path

COUNTER=0
LEN=$(ls $1/*fasta | wc -l)
for filename in $(ls $1 | head -n$LEN); do
    ./src/dnapunctuation $1/$filename 10 50 3 15 5
    COUNTER=$[COUNTER + 1]
    if (($(($COUNTER%5)) == 0))
    then
        echo -ne "\r#Processing $COUNTER out of $LEN"
    fi
done;
echo -ne "\r#Done, processed $COUNTER files\n"