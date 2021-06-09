#!/bin/bash

rm ../translator_factory/examples_collect
cd ../translator_factory
mkdir examples_collect
translator_file=examples_collect
path=../Origin_file/examples
files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        hfile=$(ls $path/$filename/doc)
        doxfile=$(ls $path/$filename/doc/*.dox)
        echo $doxfile
        for doxfilename in $doxfile
        do
                echo $doxfilename  >> $translator_file/$filename
                cat $doxfilename  >>  $translator_file/$filename
                echo -e "\n"  >> $translator_file/$filename
        done
done




# 合并
cd examples_collect
cat * > examples_collect_all.h

# to
../../translator_tools/examples_to.py examples_collect_all.h

# 等分 （翻译大小限制）
split -b 800k examples_collect_all_0.txt




