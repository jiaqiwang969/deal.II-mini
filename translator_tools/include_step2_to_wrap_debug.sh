#!/bin/bash
# coding=UTF-8



cd ../translator_factory
rm -rf include_collect/deal.II
mkdir include_collect
translator_file=include_collect
cp -rf ../Origin_file/deal.II include_collect
path=include_collect/deal.II

files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        hfile=$(ls $path/$filename/*.h)
        for hfilename in $hfile
        do
                ../translator_tools/include_to.py $hfilename 
                txtfile_0=`echo $hfilename | sed 's/.h$/\_0.txt/'`
                ../translator_tools/include_from.py $txtfile_0  
                hfile_0=`echo $hfilename | sed 's/.h$/\_0.h/'`
                hfile_0_0=`echo $hfilename | sed 's/.h$/\_0_0.h/'`
                ../translator_tools/include_wrap_debug.py  $hfile_0 > $hfile_0_0
        done

done

