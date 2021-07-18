#!/bin/bash
# # 删除所有comm
#path=deal.II-translator
files=$(ls $pwd)
for filename in $files
do
        echo $filename 
        rm -rf $filename/*_commands 
        rm -rf $filename/*_comment
        rm -rf $filename/*_latex
        rm -rf $filename/*_postcommands
        rm -rf $filename/*.txt
        rm -rf $filename/*.origin
        rm -rf $filename/*.temp
        rm -rf $filename/*.tempT	

done

