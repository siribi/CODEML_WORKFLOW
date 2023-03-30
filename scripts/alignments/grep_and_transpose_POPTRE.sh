#!/bin/bash
# Job name: grep_commands
FALIST=`ls -1 *.dist.sorted`
i=1
for FAFILE in $FALIST
        do
        FILEBASE=${FAFILE%.dist.sorted}
        echo $FILEBASE
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Poptre") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Popides") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Salpur") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Salsuc") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Avecar") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt      
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cansat") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Ficmic") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Distri") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Erijap") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Malbac") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Prunavi") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Prunper") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Pyrbet") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cucsat") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cucarg") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Momcha") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Ampedge") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Cercan") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Chafas") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Dalodo") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Faialb") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Glymax") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Phavul") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Vigung") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Alnglu") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Betpen") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Carfan") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Casgla") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Casmol") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Fagsyl") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Querob") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Carill") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Jugreg") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Jatcur") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Merann") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Linusi") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Kanobo") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Aratha") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
        awk '{for(i = 1; i <= NF; i++) {if($i ~ "Vitvin") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk '{for(i = 1; i <= NF; i++) {if($i ~ "Macint") {print $i; exit}}}' $FILEBASE.dist.sorted >> $FILEBASE.oneofeach.txt
	awk  '{ printf( "%s ", $1 ); } END { printf( "\n" ); }' $FILEBASE.oneofeach.txt > $FILEBASE.oneofeach.transposed.txt
	done
