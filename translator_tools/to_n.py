#!/usr/bin/python3
#-*- coding: UTF-8 -*- 

# -----------------------------------------
# deal.ii translator for doxygen documents
# Jiaqi-Knight 2021
# based on Dmitry R. Gulevich 2020 getxfix
# -----------------------------------------
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


text = source
postamble = []

# # 删除 * 注释，该作用是对标号降级，稍微会有一些影响，但无妨
# regex = r"^\s\*|^\/\/.*"
# subst = ""
# # You can manually specify the number of replacements by changing the 4th argument
# text = re.sub(regex, subst, text, 0, re.MULTILINE)

# 首先将\r\r的扩充几行，以免被”吃掉“
regex0 = r"\n{2}|\n\s\n|[\n](?=[@\\]f{align\*}|[@\\]f{align}|[@\\]f{eqnarray}|[@\\]f{eqnarray\*}|@f{equation\*}|@f{multline\*}|@f{gather\*}|@code|@verbatim|[@\\]f\[|\@f|\<h[0-9]\>|@dot|@note|@htmlonly|@verbatim|@code|@f{align.*}|\s+-|-|\*\/)"
regex1 = r"(?<=.dox)[\n]|(?<=@endhtmlonly)[\n]|(?<=@endverbatim)[\n]|(?<=@endcode)[\n]|(?<=@f[}\]])[\n]"
subst0 = "\n\n\n"
text = re.sub(regex0, subst0, text, 0)
text = re.sub(regex1, subst0, text, 0)

# Hide LaTeX constructs \begin{...} ... \end{...}
latex = []
start_values = []
end_values = []
for m in re.finditer(r'@dot|@htmlonly|[@\\]f{align\*}|[@\\]f{align}|[@\\]f{eqnarray}|[@\\]f{eqnarray\*}|@f{equation\*}|@f{equation}|@f{multline\*}|@f{gather\*}|@code|@verbatim|[@\\]f\[|\<a[ \n]|\<h[1-9]\>|\<[bi]\>', text):
    start_values.append(m.start())
for m in re.finditer(r'@enddot|@endhtmlonly|@f}|@endcode|@endverbatim|\\f}|[@\\]f\]|\<\/a\>|\<\/h[1-9]\>|<\/[bi]\>', text):
    end_values.append(m.end())
nitems = len(start_values)
assert(len(end_values) == nitems)
if(nitems > 0):
    newtext = text[:start_values[0]]
    for neq in range(nitems-1):
        latex.append(text[start_values[neq]:end_values[neq]])
        newtext += '[1.x.%d]' % (len(latex)-1) + \
            text[end_values[neq]:start_values[neq+1]]
    latex.append(text[start_values[nitems-1]:end_values[nitems-1]])
    newtext += '[1.x.%d]' % (len(latex)-1) + text[end_values[nitems-1]:]
    text = newtext

# Replace LaTeX commands, formulas and comments by tokens
recommand = re.compile(
    r'\/\*\*|\*\/|@ref\s\w+|@cite\s\w+|@note|@dealiiTutorialDOI.*|@dealiiVideoLecture.*|@image.*|@include.*|@page.*|(\$+)(?:(?!\1)[\s\S])*\1|(\<code\>)[\s\S].*(\<\/code\>+)|\<img.*[\s\S]\>|\<em\>|\<\/em\>|\S+\:\:\S+|@p \S+|<table\s\w{5}=\"[^>]*?>[\s\S]*?<\/table>|<table>[\s\S]*?<\/table>|\<li\>|\<\/li\>|\<\/ol\>|\<ul\>|\<\/ul\>|\<br\>|\<p\>|\<\/p\>|\@\S+')
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

doxygen_latex = filebase+'_doxygen_latex'
doxygen_commands = filebase+'_doxygen_commands'
with open(doxygen_latex, 'wb') as fp:
    pickle.dump(latex, fp)
with open(doxygen_commands, 'wb') as fp:
    pickle.dump(commands, fp)


# 最后将换行户转化为空格，为了更好的翻译目的
# 接着，开始替换
regex = r"(?<=[^\n])(\n)(?=[^-])"
#regex = r"ssssssss"
subst = " "
# You can manually specify the number of replacements by changing the 4th argument
text = re.sub(regex, subst, text, 0, re.MULTILINE)


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
