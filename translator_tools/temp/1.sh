#!/bin/bash

path=.
files=$(ls $path)
for filename in $files
do
        # echo $path/$filename 
        hfile=$(ls $path/$filename)
        doxfile=$(ls $path/$filename/*.md)
        #echo $doxfile
        for doxfilename in $doxfile
        do
                echo $doxfilename  >> merge.txt
                cat $doxfilename  >>  merge.txt
                echo -e "\n"  >> merge.txt
        done
done




# 合并
#cd examples_collect
#cat * > examples_collect_all.h