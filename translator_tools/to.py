#!/usr/bin/python3
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

if(re.search('.h$', args.filename) == None):
    sys.exit('The input should be .h file. Exit.')

print('LaTeX file:', args.filename)
filebase = re.sub('.h$', '', args.filename)
with open(args.filename, 'r') as source_file:
    source = source_file.read()

# Search for possible token conflicts
conflicts = re.findall('\[ *[012][\.][0-9]+\]', source)
if(conflicts != []):
    print('Token conflicts detected: ', conflicts)
    sys.exit(
        'Tokens may overlap with the content. Change tokens or remove the source of conflict.')
else:
    print('No token conflicts detected. Proceeding.')


text = '*/'+source+'/*'
postamble = []

# # Step0 全局注释
# regex = r"^"
# subst = "##"
# # You can manually specify the number of replacements by changing the 4th argument
# text = re.sub(regex, subst, text, 0, re.MULTILINE)


# Step0
postrecommand = re.compile(r'\/\*.*\n {0,3}\*\/|\/\*.*\*\/')  # fix bug with "/* \n */" type
postcommands = []
for m in postrecommand.finditer(text):
    postcommands.append(m.group())
nc = 0


def prerepl_f(obj):
    global nc
    nc += 1
    return ' [4.x.%d] ' % (nc-1)


text = postrecommand.sub(prerepl_f, text)


# Step1 掩盖特殊字符区块，第一类型
# Hide LaTeX constructs \begin{...} ... \end{...}
# Step4 掩盖特殊字符区块，第一类型
# Hide LaTeX constructs \begin{...} ... \end{...}
comment = []
start_values = []
end_values = []
for m in re.finditer(r'\*\/', text):
    start_values.append(m.start())
for m in re.finditer(r'\/\*', text):
    end_values.append(m.end())
nitems = len(start_values)
assert len(end_values) == nitems,filebase
if(nitems > 0):
    newtext = text[:start_values[0]]
    for neq in range(nitems-1):
        comment.append(text[start_values[neq]:end_values[neq]])
        newtext += '[0.x.%d]' % (len(comment)-1) + \
            text[end_values[neq]:start_values[neq+1]]
    comment.append(text[start_values[nitems-1]:end_values[nitems-1]])
    newtext += '[0.x.%d]' % (len(comment)-1) + text[end_values[nitems-1]:]
    text = newtext


# Step3 取消注释部分对注释
# 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
regex = r" \*"
subst = ""
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)


# Step2 整理行间距
# 首先将\r\r的扩充几行，以免被”吃掉“
regex0 = r"\n{2}|\n\s\n|[\n ](?=\[0.x.|This tutorial depends on|@dot|@note|@htmlonly|@verbatim|@code|@f{align.*}|\s+-|-|\*\/|@defgroup|@ingroup|@brief|@section|@subsection)"
regex1 = r"(?<=\n\n)|(?<=@endhtmlonly)[\n]|(?<=@endverbatim)[\n]|(?<=@endcode)[\n]|(?<=@f})[\n]"
subst0 = "\n* \n"  # "* " 不会被wrapcomments给剔除，但“*”会
text = re.sub(regex0, subst0, text, 0)
text = re.sub(regex1, subst0, text, 0)


# Step4 掩盖特殊字符区块，第一类型
# Hide LaTeX constructs \begin{...} ... \end{...}
latex = []
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
        latex.append(text[start_values[neq]:end_values[neq]])
        newtext += '[1.x.%d]' % (len(latex)-1) + \
            text[end_values[neq]:start_values[neq+1]]
    latex.append(text[start_values[nitems-1]:end_values[nitems-1]])
    newtext += '[1.x.%d]' % (len(latex)-1) + text[end_values[nitems-1]:]
    text = newtext


# Step5 掩盖特殊字符区块，第2类型
# 将“!!.*” 非翻译的部分掩盖，翻译完之后，再替换成\n
# Replace LaTeX commands, formulas and comments by tokens
recommand = re.compile(
    r'\@relatesalso.*|\@name.*|\@addtogroup.*|\@defgroup.*|\@ingroup.*|\@brief.*|\/\*\*|\*\/|@ref\s\w+|@cite\s\w+|@note|@dealiiTutorialDOI.*|@dealiiVideoLecture.*|@image.*|@include.*|[sS]tep-\d{0,2}|@page.*|<div\s\w{5}=\"[^>]*?>[\s\S]*?<\/div>[\s\S]<\/div>|<div\s\w{5}=[^>]*?>[\s\S]*?<\/div>|(\$+)(?:(?!\1)[\s\S])*\1|(\<code\>)[\s\S].*(\<\/code\>+)|\<img.*[\s\S]\>|\<em\>|\<\/em\>|\S+\:\:\S+|@p \S+|<table\s\w{5}=\"[^>]*?>[\s\S]*?<\/table>|<table>[\s\S]*?<\/table>|<p\s\w{5}=\"[^>]*?>[\s\S]*?<\/p>|\<ol\>|\<li\>|\<\/li\>|\<\/ol\>|\<ul\>|\<\/ul\>|\<br\>|\<p\>|\<\/p\>|\@\S+|\<dl\>|\<\/dl\>|\<dd\>|\<\/dd\>')
commands = []
for m in recommand.finditer(text):
    commands.append(m.group())
nc = 0


def repl_f(obj):
    global nc
    nc += 1
    return ' [2.x.%d] ' % (nc-1)


text = recommand.sub(repl_f, text)


# save doxygen commands file
doxygen_postcommands = filebase+'_doxygen_postcommands'
doxygen_comment = filebase+'_doxygen_comment'
doxygen_latex = filebase+'_doxygen_latex'
doxygen_commands = filebase+'_doxygen_commands'
with open(doxygen_postcommands, 'wb') as fp:
    pickle.dump(postcommands, fp)
with open(doxygen_comment, 'wb') as fp:
    pickle.dump(comment, fp)
with open(doxygen_latex, 'wb') as fp:
    pickle.dump(latex, fp)
with open(doxygen_commands, 'wb') as fp:
    pickle.dump(commands, fp)


# 最后将换行户转化为空格
regex = r"(?<=[^\n\*])(\n)(?=[^-\*])"
subst = " "
text = re.sub(regex, subst, text, 0, re.MULTILINE)


# # 最后对list修复，可被warpcomments.py识别
regex = r"^ {0,5}\-#{0,1}[^\>]"
subst = "*  - "
text = re.sub(regex, subst, text, 0, re.MULTILINE)


# Save the processed output to .txt file
limit = 300000  # Estimated Google Translate character limit
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
