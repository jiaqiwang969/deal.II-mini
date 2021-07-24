#!/bin/bash
# coding=UTF-8
cd ../

path=Translator_file/deal.II/design_pattern
rm  -rf temp/design_pattern
mkdir temp/design_pattern
fixpath=temp/design_pattern
files=$(ls $path)
for filename in $files
do
        #echo $path/$filename 
        #hfile= `echo $filename  | sed 's/.h$/\_0.txt/'`
        ./translator_tools/wrapcomments.py $path/$filename > $fixpath/$filename 

        # # 生成 "_0.txt" 分解文档
        # ./contrib/translator/to.py $path/$filename/$hfilename
        # # 跳过翻译，合并文档
        # txtfile_0=`echo $hfilename | sed 's/.h$/\_0.txt/'`
        # hfile0=`echo $hfilename | sed 's/.h$/\_0.h/'`
        # hfile1=`echo $hfilename | sed 's/.h$/.origin/'`
        # hfile2=`echo $hfilename | sed 's/.h$/.temp/'`
        # ./contrib/translator/from.py $path/$filename/$txtfile_0  # 输出hfile0 (未修复版)
        # mv $path/$filename/$hfilename $path/$filename/$hfile1 # bak the origin file
        # ./contrib/translator/wrap.py  $path/$filename/$hfile0 >> $path/$filename/$hfilename
        # mv $path/$filename/$hfile0 $path/$filename/$hfile2 # bak the file

done

cp -rf temp/design_pattern/* Translator_file/deal.II/design_pattern

