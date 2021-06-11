#!/bin/bash

cd ..
path=Translator_file/examples
files=$(ls $path)
for filename in $files
do
     #echo $filename.cc
     echo $path/$filename/step-*.cc
     rm -rf $path/$filename/step-*.cc
     cp -rf CCTest/CCTest_file/$filename.cc $path/$filename

     
#      rm -rf $path$/filename/step-*.cc
done


