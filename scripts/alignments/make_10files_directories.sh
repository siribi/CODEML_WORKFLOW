#!/bin/sh
i=0; 
for f in *; 
do 
    d=dir.$(printf %1d $((i/10+1))); 
    mkdir -p $d; 
    mv "$f" $d; 
    let i++; 
done