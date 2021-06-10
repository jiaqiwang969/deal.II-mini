#!/bin/bash
#rm ../Origin_file/examples
#cp -rf ../dealii/examples ../Origin_file/
cd ../translator_factory
rm examples_collect
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
cd ../../translator_tools
ls
pwd
./examples_to.py ../translator_factory/examples_collect/examples_collect_all.h

# # 等分 （翻译大小限制）
# #split -b 800k examples_collect_all_0.txt


# # from
./examples_from.py ../translator_factory/examples_collect/examples_collect_all_0.txt

## delete 
for filename in $files
do
        echo $path/$filename 
        hfile=$(ls $path/$filename/doc)
        rm -rf $path/$filename/doc/*.dox
done

# ## 
cd ../translator_factory
awk '/..\/Origin_file\//{close(out); out = $1 } {print  > out}' examples_collect/examples_collect_all_0.h

