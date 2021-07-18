# deal.II-translator
Using Doxygen format to translate packages into different languages-https://www.dealii.org/

- 首先，我们利用"./to.py filename" 生成剔除特殊字符的中转文件；
- 然后，利用deepl或者google翻译；
- 然后，利用"./from.py filename_0_T.txt";
- 然后，得到的替换文件回代到source中，再一次进行编译。注意编译需要开启mathjax，否则公式可能乱码"cmake -DDEAL_II_COMPONENT_DOCUMENTATION=ON -DDEAL_II_DOXYGEN_USE_MATHJAX=ON",然后再输入"make documentation"，若有改动，只需输入后者，无需重头开始编译。

- 可能用到的批量操作代码：
  - 拆分一次性翻译的文档：
awk '/\[2.x.0\]/{n++;w=1} n&&w{print >"step-"n"_0_T.txt"}' T.txt
  - 重新批量命名
for file in `ls *.h`;do mv $file `echo $file|sed 's/_0_T//g'`;done;
  - 批量允许python代码
for filename in headers/*.h; do
        ./to.py "$filename" 
done
for filename in examples/*_T.txt; do
        ./from.py ${filename} 
done
-在当前目录下，将所有aaaModule都替换为bbbName
grep -rl 'aaaModule' ./  | xargs sed -i "" "s/aaaModule/bbbName/g"

-r 表示搜索子目录
-l 表示输出匹配的文件名


# 将h文件全部翻译，并进行CI-测试代码 
使用说明：
- step1: terminal 运行 `translator_debug.sh > log1.txt` 用来测试，如果有错误，需要进行bug调试
- step2: terminal 运行 `translator_to.sh > log2.txt` 进行元素提取，整理打包文件在transltor_file 找到。 运行 `cat * > merge.txt` 合并为一个。
- step3: 在Testlab文件夹中，替换已经翻译好的merge.txt， 然后运行  `bash_T.sh` .
- step4: 最后进行整理，替换。 其中，tutorial有一些文件在修复wrapcomments过程中，出错为空，暂且替换为未修复版本。

# bug备忘录
- bug-01: 翻译误差
  - 案例：base/polunomials_p.h,fe/fe_simple_p.h 文件的代码块曾出现重复现象，，暂未解决，通过报错找到，并人工修复。 具体位置见commit。2021-06-07 
  - 原因: [,x,]在翻译的时候，deepl在短句上面，会重复翻译一次，导致的！
  - 解决办法: 调试会出错，根据出错信息定位，手动修复。
  


  - bug02: 输出为空。
    - 原因01：举例  `/* Equation 9: Complementary slackness */`  不要换行，否则无法识别。在过程中，会报错，进行修复。
    - 原因02：举例  `@code/* Equation 9: Complementary slackness */@endcode`  不要换行，否则无法识别。在过程中，会报错，进行修复。目前策略，修改为`@code // Equation 9: Complementary slackness  @endcode`。 提示：vscode文字颜色变化会提示错误

  - bug03: deepl 错误翻译
    - 案例1: one (which we do by using the second template argument that isusually defaulted to equal the first; here, we will create objects [2.x.100] dimensional cells in a  [2.x.101]  dimensional space).An element  [2.x.102]  of  [2.x.103]  is uniquelyidentified by the vector  [2.x.104]  of its coefficients [2.x.105] , that is:[1.x.71]where summation  is implied over repeated indexes. Note that we could usediscontinuous elements here &mdash; in fact, there is no real reason to usecontinuous ones since the integral formulation does notimply any derivatives on our trial functions so continuity is unnecessary,and often in the literature only piecewise constant elements are used.
    - 翻译：我们将有限维空间[2.x.98]定义为[1.x.70]的基函数[2.x.99]，我们将使用通常的FE_Qfinite元素，但这次它被定义在一个一维的流形上（我们通过使用第二个模板参数，通常默认为等于第一个；这里，我们将创建对象[2.x.100 2.x.103]的一个元素[2.x.102]被其系数[2.x.105]的向量[2.x.104]所唯一识别，即：[1.x.71]，其中的求和隐含在重复索引中。请注意，我们可以在这里使用非连续元素；事实上，没有真正的理由使用非连续元素，因为积分公式并不意味着我们的试探函数的任何导数，所以连续性是不必要的，而且在文献中通常只使用片断常数元素。

    - 操作： 通过regx `\[{0,1}[0124][\.\, ][xX][\.\,][0-9]+[. ]x[. ][0-9]+\]{0,1} `将其在翻译文件中找到，并人工将其分开。[2.x.100]。 [2.x.103]






  # 目前已完成的问题
  - CI测试通过，说明，include代码没有问题，翻译完毕。
  
  
  # 仍在测试
  - doc header 头文件 和 tutorial 头文件 。
