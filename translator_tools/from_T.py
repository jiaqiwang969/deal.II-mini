#!/usr/bin/python3
#-*- coding: UTF-8 -*- 
# -----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
# -----------------------------------------
# 目的： 为了解决header文库翻译内容为注释内容，但非注释内容不翻译的问题
# 多加了一个comment数据库，储存comment


import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.txt$', args.filename) == None):
    sys.exit('The input should be .txt file. Exit.')

print('Input file:', args.filename)


filebase = re.sub('_0.txt$', '', args.filename)
doxygen_comment = filebase+'_doxygen_comment'
doxygen_latex = filebase+'_doxygen_latex'
doxygen_commands = filebase+'_doxygen_commands'
# pre change to post, 最先应用，最后解析，嵌套
doxygen_postcommands = filebase+'_doxygen_postcommands'


# Load LaTeX data from binary files
with open(args.filename, 'r') as fin:
    source = fin.read()
with open(doxygen_comment, 'rb') as fp:
    comment = pickle.load(fp)
with open(doxygen_commands, 'rb') as fp:
    commands = pickle.load(fp)
with open(doxygen_latex, 'rb') as fp:
    latex = pickle.load(fp)
with open(doxygen_postcommands, 'rb') as fp:
    postcommands = pickle.load(fp)

# Replace weird characters introduced by translation
trtext = re.sub('\u200B', ' ', source)

# Fix spacing
trtext = re.sub(r'\\ ', r'\\', trtext)
trtext = re.sub(' ~ ', '~', trtext)
trtext = re.sub(' {', '{', trtext)

# Restore LaTeX and formulas
here = 0
newtext = ''
corrupted = []
for m in re.finditer('\[{0,1}[0124][\.\, ][xX][\.\,][0-9]+\]{0,1}', trtext):
    print(m)
    t = int(
        re.search('(?<=\[{0})[0124](?=[\.\,][xX][\.\,])', m.group()).group())
    n = int(
        re.search('(?<=[\.\,][xX][\.\,])[0-9]+(?=\]{0})', m.group()).group())

    if(t == 1):
        newtext += trtext[here:m.start()] + latex[n]
    elif (t == 0):
        newtext += trtext[here:m.start()] + comment[n]
    elif (t == 4):
        newtext += trtext[here:m.start()] + postcommands[n]
    elif(t == 2):
        newtext += trtext[here:m.start()] + " " + commands[n] + " "
    here = m.end()
newtext += trtext[here:]
trtext = newtext


# 再解析一遍，防止嵌套关系，有遗漏
here = 0
newtext = ''
corrupted = []
for m in re.finditer('\[{0,1}[0124][\.\, ][xX][\.\,][0-9]+\]{0,1}', trtext):
    print(m)
    t = int(
        re.search('(?<=\[{0})[0124](?=[\.\,][xX][\.\,])', m.group()).group())
    n = int(
        re.search('(?<=[\.\,][xX][\.\,])[0-9]+(?=\]{0})', m.group()).group())

    if(t == 1):
        newtext += trtext[here:m.start()] + latex[n]
    elif (t == 0):
        newtext += trtext[here:m.start()] + comment[n]
    elif (t == 4):
        newtext += trtext[here:m.start()] + postcommands[n]
    elif(t == 2):
        newtext += trtext[here:m.start()] + " " + commands[n] + " "
    here = m.end()
newtext += trtext[here:]
trtext = newtext

# # 最后再进行一些修剪的工作
regex = r"\/\*\。$|\/\*$|^\*\/"
subst = ""
trtext = re.sub(regex, subst, trtext, 0, re.MULTILINE)
trtext = '//'+trtext


# Save the processed output to .tex file
output_filename = re.sub('.txt$', '.h', args.filename)
with open(output_filename, 'w') as translation_file:
    translation_file.write(trtext)
print('Output file:', output_filename)

# Report the corrupted tokens
if(corrupted == []):
    print('No corrupted tokens. The translation is ready.')
else:
    print('Corrupted tokens detected:', end=' ')
    for c in corrupted:
        print(c, end=' ')
    print()
    print('To improve the output manually change the corrupted tokens in file',
          args.filename, 'and run from.py again.')
