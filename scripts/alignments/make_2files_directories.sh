#!/bin/sh
i=0; 
for f in *; 
do 
    d=dir.$(printf %1d $((i/2+1))); 
    mkdir -p $d; 
    mv "$f" $d; 
    let i++; 
done
