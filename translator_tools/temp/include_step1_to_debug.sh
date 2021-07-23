#!/bin/bash
# coding=UTF-8


cd ../translator_factory
mkdir include_collect
translator_file=include_collect
cp -rf ../Origin_file/deal.II include_collect
path=include_collect/deal.II

files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        # step2-main code
        hfile=$(ls $path/$filename/*.h)
        for hfilename in $hfile
        do
                ../translator_tools/include_to.py $hfilename # debug
        done

done

