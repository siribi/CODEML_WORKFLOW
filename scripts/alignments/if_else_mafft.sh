#!/bin/sh
FALIST=`ls -1 *.fasta`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.fasta}
        echo $FILEBASE
        if mafft --thread -1 --auto $FILEBASE.fasta > $FILEBASE.aligned.fa | grep -m 1 "Illegal"; then
        mv $FILEBASE.aligned.fa $FILEBASE.ERROR.fa        
	continue
        fi
        done
