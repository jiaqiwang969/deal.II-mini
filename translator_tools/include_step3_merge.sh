#!/bin/bash



cd ../translator_factory
translator_file=include_collect/collection
mkdir $translator_file
path=../translator_factory/include_collect/deal.II
files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        txtfile_0=$(ls $path/$filename/*_0.txt)
        for txtfilename_0 in $txtfile_0
        do
                echo $txtfilename_0  >> $translator_file/$filename
                cat $txtfilename_0  >>  $translator_file/$filename
                echo -e "\n"  >> $translator_file/$filename
        done
done


pwd

# # # 合并
cd include_collect/collection
cat * > include_collect_all.txt


# # # 等分 （翻译大小限制）
split -b 800k include_collect_all.txt








