#!/usr/bin/python
# coding=utf8
# -----------------------------------------
# (?<=0.x.)([0-9]+)(?=\] \n )|
# (?<=\] \n)( )(?=\n\[)


import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

# only work for merge_0_AUTODEBUG.txt
if(re.search('.txt$', args.filename) == None):
    sys.exit('The input should be .h file. Exit.')

filebase = re.sub('.txt$', '', args.filename)
with open(args.filename, 'r') as source_file:
    text = source_file.read()



# -------- 操作块 ---------
# 一个空行填充
regex = r"(?<=0.x.)([0-999999]+\] \n )(?=)"  # 插入地方的关键字
matches = re.finditer(regex, text, re.MULTILINE)

for matchNum, match in enumerate(matches, start=1):

    # 替换的内容，得先进行数字提取和转化
    regex0 = '] \n '
    subst = ""
    mat=match.group()
    num = int(re.sub(regex0, subst, mat, 0 ))
    if (num ==28289):
        a=1
    num2=str(num+1)
    subst1= mat[:-1]+"[0.x." + num2 +"] "
    text = re.sub(mat, subst1, text, 0 )

# # 两个空行填充
# regex = r"(?<=0.x.)([0-9]+\] \n )(?=)"  # 插入地方的关键字
# matches = re.finditer(regex, text, re.MULTILINE)

# for matchNum, match in enumerate(matches, start=1):

#     # 替换的内容，得先进行数字提取和转化
#     regex0 = '] \n '
#     subst = ""
#     mat=match.group()
#     num = int(re.sub(regex0, subst, mat, 0 ))
#     subst1= mat[:-1]+"[0.x." + str(num+1) +"] "
#     text = re.sub(mat, subst1, text, 0 )



# Save the processed output to .txt file
limit = 30000000  # Estimated Google Translate character limit
start = 0
npart = 0
for m in re.finditer(r'\.\n', text):
    if(m.end()-start < limit):
        end = m.end()
    else:
        output_filename = filebase+'_%d.cc' % npart
        npart += 1
        with open(output_filename, 'w') as txt_file:
            txt_file.write(text[start:end])
        print('Output file:', output_filename)
        start = end
        end = m.end()
output_filename = filebase+'_%d.cc' % npart
with open(output_filename, 'w') as txt_file:
    txt_file.write(text[start:])
print('Output file:', output_filename)
print('Supply the output file(s) to Google Translate')
