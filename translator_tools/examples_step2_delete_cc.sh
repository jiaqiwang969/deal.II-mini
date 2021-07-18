#!/bin/bash

cd ../CCTest
rm -rf CCTest_file/*
awk '/CCTest_file\//{close(out); out = $1 } {print  > out}' merge_0_T.cc
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


