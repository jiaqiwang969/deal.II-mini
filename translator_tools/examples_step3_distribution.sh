#!/bin/bash

# # ## step-1
# # # 删除deal.II-translator下所有文件，但保留文件夹和comm
# cd include
# rm -rf deal.II-translator
# cp -rf ../../include/deal.II-translator .

# path=deal.II-translator
# files=$(ls $path)
# for filename in $files
# do
#         echo $path/$filename 
#         rm -rf $path/$filename/*.h # 删除h文件，保留comment文件
#         rm -rf $path/$filename/*.txt # 删除txt文件，保留comment文件
# done


path=examples
files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        hfile=$(ls $path/$filename/doc)
        doxfile=$(ls $path/$filename/doc/*.dox)
        echo $doxfile
        for doxfilename in $doxfile
        do
                rm -rf $doxfilename
        done
done





# # ## step-2 分离翻译文件
# cd ..
# awk '/include\//{close(out); out = $1 } {print  > out}' post/A_T.txt

# ## step-3 转化
# path=include/deal.II-translator
# files=$(ls $path)
# for filename in $files
# do
#         #echo $path/$filename 
#         # step2-main code
#         txtfile=$(ls $path/$filename/*0.txt)
#         for txtfilename in $txtfile
#         do
#                 # 跳过翻译，合并文档
#                 echo $txtfilename
#                 hfile=`echo $txtfilename | sed 's/.txt$/.h/'`
#                 hfile1=`echo $txtfilename | sed 's/\_0.txt$/.h/'`
#                 ../contrib/translator/from_T.py $txtfilename  # 输出hfile0 (未修复版)
#                 #echo '//' > $hfile1
#                 ## 开启debug模式
#                 #../contrib/translator/wrap.py  $hfile >> $hfile1
#                 ../contrib/translator/wrap.py  $hfile >> $hfile1
#                 hfile2=`echo $hfile | sed 's/\_0.h$/.tempT/'`
#                 mv $hfile $hfile2 # bak the file
#         done
# done