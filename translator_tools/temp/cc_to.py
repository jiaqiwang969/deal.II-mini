#!/usr/bin/python
#-*- coding: UTF-8 -*- 
# -----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
# -----------------------------------------
# 目的： 为了解决header文库翻译内容为注释内容，但非注释内容不翻译的问题
# 思路： 先对全部内容进行注释，//，然后，去除"// *"" 的注释内容, 对不翻译内容，打包处理【】
# 翻译完成以后（仅仅翻译注释内容），
# 再进行注释操作(前面对目的，是为了分离注释和非注释内容)

# 另外一条思路：逆向思维 （可行）
# */ comment /* 标记识别comment，然后将这部分提取为注释内容【】
# comment /* */ comment /*  */  comment
# 提取在上下添加辅助 得到： */ comment /* */ comment /*  */  comment /*
import re
import sys
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

if(re.search('.cc$', args.filename) == None):
    sys.exit('The input should be .h file. Exit.')

print('LaTeX file:', args.filename)
filebase = re.sub('.cc$', '', args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()

text = source



# 代码块掩盖
regex = r"^ {0,80}[^\/\/\n ].*?$|^\/\*.*"
commands0 = []
for m in re.finditer(regex, text, re.MULTILINE):
    commands0.append(m.group())
nc = 0

def repl_f(obj):
    global nc
    nc += 1
    return ' [0.x.%d] ' % (nc-1)

text = re.sub(regex, repl_f, text, 0, re.MULTILINE) #取代 text = recommand.sub(repl_f, text)



# 去掉//号
regex = r" {0,10}\/\/"
subst = ""
text = re.sub(regex, subst, text, 0, re.MULTILINE)


# 首先将\r\r的扩充几行，以免被”吃掉“
regex0 = r"\n{2}|\n\s\n|[\n ](?=\[0.x.)"
regex1 = r"(?<=\n\n)|(?<=\] )[\n]"
subst0 = "\n\n\n"  # "* " 不会被wrapcomments给剔除，但“*”会
text = re.sub(regex0, subst0, text, 0)
text = re.sub(regex1, subst0, text, 0)





# 关键字转化-01
commands1 = []
start_values = []
end_values = []
for m in re.finditer(r'@dot|@htmlonly|[@\\]f{align\*}|[@\\]f{align}|[@\\]f{eqnarray}|[@\\]f{eqnarray\*}|@f{equation\*}|@f{equation}|@f{multline\*}|@f{gather\*}|@code|@verbatim|[@\\]f\[|\<a[ \n]|\<h[1-9]\>|\<[bi]\>', text):
    start_values.append(m.start())
for m in re.finditer(r'@enddot|@endhtmlonly|@f}|@endcode|@endverbatim|\\f}|[@\\]f\]|\<\/a\>|\<\/h[1-9]\>|<\/[bi]\>', text):
    end_values.append(m.end())
nitems = len(start_values)
assert len(end_values) == nitems, filebase
if(nitems > 0):
    newtext = text[:start_values[0]]
    for neq in range(nitems-1):
        commands1.append(text[start_values[neq]:end_values[neq]])
        newtext += '[1.x.%d]' % (len(commands1)-1) + \
            text[end_values[neq]:start_values[neq+1]]
    commands1.append(text[start_values[nitems-1]:end_values[nitems-1]])
    newtext += '[1.x.%d]' % (len(commands1)-1) + text[end_values[nitems-1]:]
    text = newtext


# 关键字转化-02
recommand = re.compile(
    r'\@sect.*|\@relatesalso.*|\@name.*|\@addtogroup.*|\@defgroup.*|\@ingroup.*|\@brief.*|\/\*\*|\*\/|@ref\s\w+|@cite\s\w+|@note|@dealiiTutorialDOI.*|@dealiiVideoLecture.*|@image.*|@include.*|[sS]tep-\d{0,2}|@page.*|<div\s\w{5}=\"[^>]*?>[\s\S]*?<\/div>[\s\S]<\/div>|<div\s\w{5}=[^>]*?>[\s\S]*?<\/div>|(\$+)(?:(?!\1)[\s\S])*\1|(\<code\>)[\s\S].*(\<\/code\>+)|\<img.*[\s\S]\>|\<em\>|\<\/em\>|\S+\:\:\S+|@p \S+|<table\s\w{5}=\"[^>]*?>[\s\S]*?<\/table>|<table>[\s\S]*?<\/table>|<p\s\w{5}=\"[^>]*?>[\s\S]*?<\/p>|\<ol\>|\<li\>|\<\/li\>|\<\/ol\>|\<ul\>|\<\/ul\>|\<br\>|\<p\>|\<\/p\>|\@\S+|\<dl\>|\<\/dl\>|\<dd\>|\<\/dd\>')
commands2 = []
for m in recommand.finditer(text):
    commands2.append(m.group())
nc = 0

def repl_f(obj):
    global nc
    nc += 1
    return ' [2.x.%d] ' % (nc-1)

text = recommand.sub(repl_f, text)






# save doxygen commands file
doxygen_commands0 = filebase+'_doxygen_commands0'
with open(doxygen_commands0, 'wb') as fp:
    pickle.dump(commands0, fp)
doxygen_commands1 = filebase+'_doxygen_commands1'
with open(doxygen_commands1, 'wb') as fp:
    pickle.dump(commands1, fp)
doxygen_commands2 = filebase+'_doxygen_commands2'
with open(doxygen_commands2, 'wb') as fp:
    pickle.dump(commands2, fp)





# # # 最后将换行转化为空格
regex = r"(?<=[^\n\*])(\n)(?=[^-\*])"
subst = ""
text = re.sub(regex, subst, text, 0, re.MULTILINE)

# 删去过多的空格
regex = r"\] \n{0,20} {0,10}\["
subst = "] \\n["
text = re.sub(regex, subst, text, 0, re.MULTILINE)

regex = r"\n{4,120}"
subst = "\\n\\n"
text = re.sub(regex, subst, text, 0, re.MULTILINE)

# list fix
regex = r" - "
subst = "\\n\\n - "
text = re.sub(regex, subst, text, 0, re.MULTILINE)



# 重新注释
regex = r"^"
subst = "//"
text = re.sub(regex, subst, text, 0, re.MULTILINE)

regex = r"^\/\/\[0.x."
subst = "[0.x."
text = re.sub(regex, subst, text, 0, re.MULTILINE)



# # # 最后对list修复，可被warpcomments.py识别
# regex = r"^ {0,5}\-#{0,1}[^\>]"
# subst = "*  - "
# text = re.sub(regex, subst, text, 0, re.MULTILINE)


# Save the processed output to .txt file
limit = 300000000  # Estimated Google Translate character limit
start = 0
npart = 0
for m in re.finditer(r'\.\n', text):
    if(m.end()-start < limit):
        end = m.end()
    else:
        output_filename = filebase+'_%d.txt' % npart
        npart += 1
        with open(output_filename, 'w') as txt_file:
            txt_file.write(text[start:end])
        print('Output file:', output_filename)
        start = end
        end = m.end()
output_filename = filebase+'_%d.txt' % npart
with open(output_filename, 'w') as txt_file:
    txt_file.write(text[start:])
print('Output file:', output_filename)
print('Supply the output file(s) to Google Translate')
