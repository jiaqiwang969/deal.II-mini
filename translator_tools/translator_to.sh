#!/bin/bash
# coding=UTF-8
mkdir translator_file
rm -rf translator_file/*
path=include/deal.II-translator
rm -rf include/deal.II-translator
cp -rf include/deal.II-origin include/deal.II-translator
files=$(ls $path)
for filename in $files
do
        echo $path/$filename 
        # step2-main code
        hfile=$(ls $path/$filename)
        for hfilename in $hfile
        do
                # 生成 "_0.txt" 分解文档
                ./contrib/translator/to.py $path/$filename/$hfilename
                # 跳过翻译，合并文档
                txtfile_0=`echo $hfilename | sed 's/.h$/\_0.txt/'`
                hfile0=`echo $hfilename | sed 's/.h$/\_0.h/'`
                hfile1=`echo $hfilename | sed 's/.h$/.origin/'`
                hfile2=`echo $hfilename | sed 's/.h$/.temp/'`
                ./contrib/translator/from.py $path/$filename/$txtfile_0  # 输出hfile0 (未修复版)
                mv $path/$filename/$hfilename $path/$filename/$hfile1 # bak the origin file
                ./contrib/translator/wrap.py  $path/$filename/$hfile0 >> $path/$filename/$hfilename
                mv $path/$filename/$hfile0 $path/$filename/$hfile2 # bak the file


        done
        # setp3-合并h文件，方便一次性翻译
        txtfile=$(ls $path/$filename/*.txt)
        for txtfilename in $txtfile
        do
                echo $txtfilename  >> translator_file/$filename
                cat $txtfilename  >>  translator_file/$filename
                echo -e "\n"  >> translator_file/$filename

        done


done

# for filename in $files
# do

#         # # # step1-初始化clear
#         # rm -rf $path/$filename/*_postcommands
#         # rm -rf $path/$filename/*_commands
#         # rm -rf $path/$filename/*_latex
#         # rm -rf $path/$filename/*_comment
#         # rm -rf $path/$filename/*.txt
#         # rm -rf $path/$filename/*.temp
#         # rm -rf $path/$filename/*.origin

# done



# # step3-将翻译合并成一个文件
# result_path=translator_file
# translate_files=$(ls $result_path)
# for translate_filename in $translate_files
# do
#         cat $result_path/$translate_filename >> $result_path/all.txt
#         echo $result_path/$translate_filename
# done


