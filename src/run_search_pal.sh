#!/bin/bash

# First arg: target path
# Example: src/run_search_pal.sh processed_fasta

COUNTER=0
LEN=$(ls $1/*fasta | wc -l)
for filename in $(ls $1 | head -n$LEN); do
    # ./src/dnapunctuation $1/$filename 15 30 3 10 3
    ./src/dnapunctuation $1/$filename 10 20 0 8 2
    COUNTER=$[COUNTER + 1]
    if (($(($COUNTER%5)) == 0))
    then
        echo -ne "\r#Processing $COUNTER out of $LEN"
    fi
done;
echo -ne "\r#Done, processed $COUNTER files\n"