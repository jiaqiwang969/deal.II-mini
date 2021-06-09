//include/deal.II-translator/A-headers/automatic_and_symbolic_differentiation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

/**
 *    @defgroup auto_symb_diff Automatic and symbolic differentiation
 *
 * @brief A module dedicated to the implementation of functions and classes
 * that relate 对自动和符号微分。
 * 下面我们将非常简要地介绍什么是自动微分和符号微分，这些计算/数字方案有哪些变化，以及它们是如何被整合到deal.II的框架中的。所有这些方案的目的是在不想手工计算的情况下自动计算函数的导数或其近似值。常见的例子是，在有限元背景下，人们想解决一个非线性问题，要求一些残差
 * $F(u,\nabla u)=0$ ，其中 $F$
 * 是一个复杂的函数，需要进行微分以应用牛顿方法；以及人们得到一个与参数有关的问题
 * ${\cal A}(q,u,\nabla u) = f$  并希望对参数 $q$
 * 形成导数的情况，例如，为了优化 $q$
 * 的输出函数，或为了进行 $q$ 的敏感性分析。我们应该把
 * $q$
 * 看作是设计参数：例如，机翼的宽度或形状，选择用来建造物体的材料的硬度系数，发送到设备的功率，发送到燃烧器的气体的化学成分。在所有这些情况下，人们应该把
 * $F$ 和 $\cal A$ 看作是<i>complicated</i>，并且繁琐地加以区分
 *
 * - 至少在手工操作时是这样。 step-15 中显示了一个相对简单的非线性问题的案例，它已经突出了手工计算导数的繁琐。然而，在现实中，人们可能会想到诸如化学反应流等问题，其中流体方程的系数，如密度和粘度，强烈地、非线性地依赖于流体在每一点的化学成分、温度和压力；化学物种根据反应系数相互反应，而反应系数也非线性地、复杂地依赖于化学成分、温度和压力。在许多情况下，所有这些系数的精确公式可能要花几行字才能写出来，可能包括几个非线性项的指数和（谐波或几何）平均值，和/或可能包含数据点之间的表格查询和内插。仅仅把这些项写对就已经很困难了；计算这些项的导数在大多数应用中是不切实际的，在现实中是不可能写对的。如果没有计算机的帮助，更高一级的导数就更不可能做到。自动或符号微分是解决这个问题的一个方法。我们只需要实现一次计算这些系数的输入的函数，就可以得到（正确的！）导数，而不需要进一步的编码工作（尽管在运行时、编译时或两者都有不可忽略的计算成本）。
 *
 *
 * @section  auto_diff_1 自动微分法 <a
 * href="https://en.wikipedia.org/wiki/Automatic_differentiation">Automatic
 * differentiation </a>
 * （通常也被称为算法微分），是一种数字方法，可以用来
 * "自动
 * "计算函数对一个或多个输入变量的一阶，或许还有高阶导数。虽然这需要一定的计算成本，但使用这种工具的好处可能是巨大的。如果使用得当，通常复杂的函数的导数可以被计算得非常精确。尽管这些框架所能达到的确切精度在很大程度上取决于它们的基础数学公式，但一些实现方式的计算精度与机器精度相当。请注意，这与经典的数值微分（例如，通过在不同的点上对一个函数进行评估，使用有限差分近似）不同，后者的精度取决于扰动的大小以及所选择的有限差分方案；这些方法的误差要明显大于表述良好的自动微分方法。
 * 那么，在有限元背景下使用自动微分的三个实际例子将是
 *
 *
 * - 快速建立一个新的非线性公式的原型，而不需要手工计算线性化本身。
 *
 *
 * - 在复杂的多物理学框架内自动线性化有限元残差加法形成，以及
 *
 * - 验证用户对基于单元的计算（如残差）和基于连续点的计算（如非线性构成法的切线）的线性化实施。
 * 有相当多的自动微分数字的实现方法。它们主要分为两大类，即
 * <em>  源代码转换  </em>  和  <em>  操作符重载  </em>
 * 。第一种方法是根据一些输入函数生成新的、可编译的代码，当执行时，返回输入函数的导数。第二种方法是利用<tt>C++</tt>操作符定义的能力，为自定义类的类型进行重载。因此，一个代表这样一个可自动微分的数字的类，在对它进行每一次数学运算后，原则上可以评估和跟踪它的值以及它的方向性导数。由于专门实现
 * <em> 源代码转换 </em>
 * 方法的库共同描述了高度专业化的工具，作为函数预处理器使用，它们在deal.II本身没有直接支持。然而，后者代表了专门的数字类型，可以通过在适当的环境中使用模板元编程来支持。鉴于上面的例子，这意味着FEValues类（和朋友），以及Tensor和SymmetricTensor类应该支持用这些专门的数字进行计算。理论上，整个程序都可以被做成可微调的。例如，这在解决方案对输入参数的敏感性分析中可能是有用的。然而，到目前为止，这还没有被测试过）。)
 * 基于 <em> 运算符重载 </em>
 * 的专门框架的实现通常分为三类之一。在每一类中，一些定制的数据类代表了被评估函数的浮点值和它的导数，由
 *
 *
 * - 利用 <em> 对偶 </em> / <em> 复杂步骤 </em> / <em> 超对偶 </em> 公式（有时也称为 <em> 无磁带 </em> 方法）。
 *
 *
 * - 那些利用 <em> 拍打 </em> 的策略，以及
 *
 *
 * - 那些通过 <em> 表达式模板 </em> 使用的编译时优化。
 * 为了初步了解这些不同的实现方式在实践中可能是什么样的，我们对这些方法做了如下的一般性总结。
 *
 *
 * - 上面列出的前两种 <em> 无龙头 </em> 方法（对偶数和复数步法）使用截断泰勒级数的某种变化，以及对扰动参数定义的特殊选择，用基于有限差分的方法来计算函数导数。双重 "数构成了在函数值被评估时同时计算的累积方向导数；在复数步法中，虚值有效地起到了这个作用。扰动参数的选择决定了方案的数值质量，例如泰勒方案截断的影响；对偶数在其一阶导数中不包含任何高阶项，而对于复数阶方法，这些存在的高阶项被忽略了。可以证明，这两种方法都不受减法取消误差的影响，而且在其有限差分方案中，它们在数值上对为数值扰动选择的内部 step- 大小不敏感。因此，对偶数方法产生精确的一阶导数，而复数步骤近似法则不能。然而，对偶数的标准实现不能产生精确的二阶导数值。超对偶数对这一想法有不同的看法，数字以类似于四元数的形式表示（即带有额外的非实数成分），导数由泰勒级数所有四个成分的高阶截断来计算。其结果是，通过适当的实现，第一和第二导数都可以准确计算出来。
 *
 * - 使用 <em> 磁带 </em> 方法，选择一个指定的代码子区域，对该区域所有用活动（标记）输入变量执行的操作进行跟踪并记录在一个被称为磁带的数据结构中。在磁带区域结束时，可以通过用不同的输入变量集 "重放 "磁带来重新评估所记录的函数，而不是直接重新计算该函数。假设录音区域代表一个平滑的函数，那么该函数的任意高阶导数可以通过参考跟踪和存储在磁带上的代码路径来计算。    例如，这可以通过对感兴趣点周围的函数进行评估来实现）。有一些策略可以处理这样的情况：在被评估的点上，所录取的函数是不平滑的，或者它不是分析的。此外，我们可能需要考虑分支函数的情况，即磁带不再是连续的，而是在不同的评估路径上分叉，而这是由于最初记录的输入。
 *
 *
 * - 基于<a href="https://en.wikipedia.org/wiki/Expression_templates">expression templates</a>的方法利用了由抽象语法树（AST）构建的计算图（在这种情况下是<a href="https://en.wikipedia.org/wiki/Directed_acyclic_graph">directed acyclic graph (DAG)</a>），该图从其输入值中解析出函数输出。    树上最外层的叶子代表独立变量或常数，并由单数运算符转换，由二进制运算符连接（在最简单的情况下）。因此，在编译时，对函数输入进行的操作是已知的，与此相关的导数操作也可以同时使用众所周知的计算运算导数的规则（如加减法下导数的关联性、乘积规则和连锁规则）来定义。这个运算器返回的编译输出类型不需要是通用的，而是可以根据DAG顶点上给该特定运算器的特定输入（可能携带微分历史）进行专业化。通过这种方式，可以为用于评估依赖函数的每个中间结果的非常专门的个别操作生成一套编译时优化的指令。
 * 当然，这些方法中的每一种都有其优点和缺点，对于要解决的特定问题，其中一种可能比另一种更合适。由于上述实施细节（以及其他未讨论的细节）可能对用户是隐藏的，因此了解使用这些
 * "黑盒子
 * "自动差分数字中的任何一个的影响、运行时间成本和潜在的限制可能仍然很重要。
 * 除了所提供的链接文章外，用于提供这里所提供的细节的资源包括。
 *
 * @code{.bib}
 * @InProceedings{Fike2011a,
 * author    = {Fike, Jeffrey A and Alonso, Juan J},
 * title     = {The Development of Hyper-Dual Numbers for Exact Second-Derivative Calculations},
 * booktitle = {49th {AIAA} Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition},
 * year      = {2011},
 * volume    = {886},
 * pages     = {124},
 * month     = {jan},
 * publisher = {American Institute of Aeronautics and Astronautics},
 * doi       = {10.2514/6.2011-886},
 * }
 * @endcode
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 * title     = {Getting Started with ADOL-C},
 * author    = {Walther, Andrea and Griewank, Andreas},
 * year      = {2009},
 * booktitle = {Combinatorial scientific computing},
 * doi       = {10.1.1.210.4834},
 * pages     = {181--202}
 * }
 * @endcode
 *
 * # ###利用链式定律
 * 在最实际的意义上，上述任何一类都是利用链式规则来计算复合函数的总导数。为了执行这个动作，它们通常使用两种机制中的一种来计算导数，具体是
 *
 *
 * -  <em> 正向模式 </em> （或 <em> 正向积累 </em> ）自动分化，或
 *
 *
 * -  <em> 反向模式 </em> （或 <em> 反向积累 </em> ）自动差分。
 * 值得注意的是， <em> 最优雅各布积算 </em>
 * ，执行一组最小的计算，位于这两种极限情况之间。它对一般复合函数的计算仍然是图论中的一个公开问题。
 * 借助于下面的图表（它和一些列出的细节由这个<a
 * href="https://en.wikipedia.org/wiki/Automatic_differentiation">Wikipedia
 * article</a>提供），让我们思考一下函数 $f (\mathbf{x}) = \sin
 * (x_{1}) + x_{1} x_{2}$ 及其导数的计算的代表性。 <div
 * class="twocolumn" style="width: 80%"> <div class="parent"> <div class="img"
 * align="center"> <img
 * src="https://upload.wikimedia.org/wikipedia/commons/a/a4/ForwardAccumulationAutomaticDifferentiation.png"
 * alt="Forward mode automatic differentiation" width="400"> </div>   <div
 * class="text" align="center"> Forward mode automatic differentiation </div>
 * </div>  <div class="parent"> <div class="img" align="center"> <img
 * src="https://upload.wikimedia.org/wikipedia/commons/a/a0/ReverseaccumulationAD.png"
 * alt="Reverse mode automatic differentiation" width="400"> </div>   <div
 * class="text" align="center"> Reverse mode automatic differentiation </div>
 * </div> </div>
 * 具体来说，我们将简要介绍什么是正向和反向的自动差分。请注意，在图中，沿着文字图形的边缘是函数
 * $w$  相对于  $i$  -个变量的定向导数，用符号  $\dot{w} =
 * \dfrac{d w}{d x_{i}}$
 * 表示。在这个例子中，用于呈现函数值及其方向导数的具体计算方法在<a
 * href="https://en.wikipedia.org/wiki/Automatic_differentiation">source
 * article</a>中列出。对于第二个说明性例子，我们请感兴趣的读者参考<a
 * href="http://www.columbia.edu/~ahd2125/post/2015/12/5/">this article</a>。
 * 首先考虑任何复合函数 $f(x)$ ，这里表示为有两个独立变量，可以剖析为其基本函数@f[
 * f (\mathbf{x})
 * = f_{0} \circ f_{1} \circ f_{2} \circ \ldots \circ f_{n} (\mathbf{x})
 * \quad .
 * @f]的组合。 如前所述，如果每个原始运算 $f_{n}$ 都是平滑的和可微的，那么可以普遍采用链式规则来计算 $f$ 的总导数，即 $\dfrac{d f(x)}{d \mathbf{x}}$ 。区分 "正向 "和 "反向 "模式的是链规的评估方式，但最终都是计算总导数@f[
 * \dfrac{d f (\mathbf{x})}{d \mathbf{x}} = \dfrac{d f_{0}}{d f_{1}} \dfrac{d
 * f_{1}}{d f_{2}} \dfrac{d f_{2}}{d f_{3}} \ldots \dfrac{d f_{n}
 * (\mathbf{x})}{d \mathbf{x}} \quad . @f] 。
 * 在正向模式下，链式规则是自然地从 "内向外 "计算的。因此自变量是固定的，每个子函数 $f'_{i} \vert_{f'_{i+1}}$ 被递归计算，其结果作为输入返回给父函数。使用圆括号封装和固定操作顺序，这意味着我们计算@f[
 * \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 * = \dfrac{d f_{0}}{d f_{1}} \left( \dfrac{d f_{1}}{d f_{2}} \left(\dfrac{d f_{2}}{d f_{3}} \left(\ldots \left( \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)\right)\right)\right)
 * \quad .
 * @f] 前向扫频的计算复杂性与输入函数的计算复杂性成正比。然而，对于每一个需要计算的方向性导数，都需要对计算图进行一次扫频。
 * 在反向模式下，链式规则的计算有点不自然地从 "外向内 "进行。因变量的值首先被计算和固定，然后前面的微分运算被评估并与之前的结果从左到右连续相乘。同样，如果我们用括号封装和固定运算顺序，这意味着反向计算是由@f[
 * \dfrac{d f (\mathbf{x})}{d \mathbf{x}}
 * = \left( \left( \left( \left( \left( \dfrac{d f_{0}}{d f_{1}} \right) \dfrac{d f_{1}}{d f_{2}} \right) \dfrac{d f_{2}}{d f_{3}} \right) \ldots \right) \dfrac{d f_{n} (\mathbf{x})}{d \mathbf{x}} \right)
 * \quad .
 * @f]进行的。中间值 $\dfrac{d f_{i-1}}{d f_{i}}$ 被称为 <em> 邻接 </em> ，必须在计算图被遍历时计算和储存。然而，对于每个从属标量函数来说，计算图的一次扫描就能一次性呈现所有的方向性导数。
 * 总的来说，每种模式的效率是由独立（输入）变量和从属（输出）变量的数量决定的。如果输出的数量大大超过输入，那么可以证明正向模式比反向模式的效率更高。反之，当输入变量的数量大大超过输出变量的数量时，也是如此。这一点可以用来帮助告知哪种数字类型最适合于使用自动微分进行的哪组操作。例如，在许多需要计算二阶导数的应用中，结合反向和正向模式是合适的。前者通常用于计算第一导数，后者用于计算第二导数。
 *
 *
 * @subsection  auto_diff_1_1 支持的自动区分库
 * 我们目前有以下数字类型和组合的验证实现。
 *
 *
 *
 * - 录音的ADOL-C（理论上是可微分的，但将实现多达二阶导数的内部驱动）。
 *
 *
 *
 * - 无带子的ADOL-C（一旦可区分）。
 *
 *
 *
 * - 使用表达式模板进行动态内存分配的前向模式Sacado（一旦可区分）。
 *
 *
 *
 * - 使用表达式模板的嵌套正向模式Sacado（可两次微分）。
 *
 *
 *
 * - 反向模式的萨卡多（一旦可微调）。
 *
 *
 *
 * - 嵌套的反向和动态分配的正向模式Sacado（两次可微分，但结果是内存泄漏，在 Differentiation::AD::NumberTypes) 中描述的
 * 注意，在上面，"动态内存分配
 * "是指在编译时不需要指定独立变量的数量。 <a
 * href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">ADOL-C
 * user manual</a>中的
 *
 * @code{.bib}
 * @Manual{Walther2009a,
 * title     = {Getting Started with ADOL-C},
 * author    = {Walther, Andrea and Griewank, Andreas},
 * year      = {2009},
 * booktitle = {Combinatorial scientific computing},
 * doi       = {10.1.1.210.4834},
 * pages     = {181--202},
 * url       = {https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf}
 * }
 * @endcode
 *
 * 提供了对其有带和无带实现的原则性见解，以及如何将ADOL-C纳入用户代码中。为了解ADOL-C的实现，以及如何在数字代码中使用它的可能性，一些进一步的有用资源包括。
 *
 * @code{.bib}
 * @Article{Griewank1996a,
 * author    = {Griewank, Andreas and Juedes, David and Utke, Jean},
 * title     = {Algorithm 755: {ADOL-C}: a package for the automatic differentiation of algorithms written in {C/C++}},
 * journal   = {ACM Transactions on Mathematical Software (TOMS)},
 * year      = {1996},
 * volume    = {22},
 * number    = {2},
 * pages     = {131--167},
 * doi       = {10.1145/229473.229474},
 * publisher = {ACM}
 * }
 * @endcode
 * @code{.bib}
 * @InCollection{Bischof2008a,
 * author =    {Bischof, Christian and Guertler, Niels and Kowarz, Andreas and Walther, Andrea},
 * title =     {Parallel reverse mode automatic differentiation for OpenMP programs with ADOL-C},
 * booktitle = {Advances in Automatic Differentiation},
 * publisher = {Springer},
 * year =      {2008},
 * pages =     {163--173}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Kulshreshtha2012a,
 * chapter   = {Computing Derivatives in a Meshless Simulation Using Permutations in {ADOL}-C},
 * pages     = {321--331},
 * title     = {Recent Advances in Algorithmic Differentiation},
 * publisher = {Springer Berlin Heidelberg},
 * year      = {2012},
 * author    = {Kshitij Kulshreshtha and Jan Marburger},
 * editor    = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 * series    = {Lecture Notes in Computational Science and Engineering},
 * doi       = {10.1007/978-3-642-30023-3_29},
 * }
 * @endcode
 * @code{.bib}
 * @InProceedings{Kulshreshtha2013a,
 * author    = {Kulshreshtha, Kshitij and Koniaeva, Alina},
 * title     = {Vectorizing the forward mode of ADOL-C on a GPU using CUDA},
 * booktitle = {13th European AD Workshop},
 * year      = {2013},
 * month     = jun
 * }
 * @endcode
 *
 * 同样，为了解Sacado数字类型的实现（特别是表达式模板的使用和利用），选择了一些有用的资源，包括：*。
 *
 * @code{.bib}
 * @InCollection{Bartlett2006a,
 * author        = {Bartlett, R. A. and Gay, D. M. and Phipps, E. T.},
 * title         = {Automatic Differentiation of C++ Codes for Large-Scale Scientific Computing},
 * booktitle     = {International Conference on Computational Science {\textendash} {ICCS} 2006},
 * publisher     = {Springer Berlin Heidelberg},
 * year          = {2006},
 * editor        = {Alexandrov, V.N. and van Albada, G.D. and Sloot, P.M.A. amd Dongarra, J.},
 * pages         = {525--532},
 * doi           = {10.1007/11758549_73},
 * organization  = {Springer}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Gay2012a,
 * chapter   = {Using expression graphs in optimization algorithms},
 * pages     = {247--262},
 * title     = {Mixed Integer Nonlinear Programming},
 * publisher = {Springer New York},
 * year      = {2012},
 * author    = {Gay, D. M.},
 * editor    = {Lee, J. and Leyffer, S.},
 * isbn      = {978-1-4614-1927-3},
 * doi       = {10.1007/978-1-4614-1927-3_8}
 * }
 * @endcode
 * @code{.bib}
 * @InBook{Phipps2012a,
 * chapter     = {Efficient Expression Templates for Operator Overloading-based Automatic Differentiation},
 * pages       = {309--319},
 * title       = {Recent Advances in Algorithmic Differentiation},
 * publisher   = {Springer},
 * year        = {2012},
 * author      = {Eric Phipps and Roger Pawlowski},
 * editor      = {Forth S. and Hovland P. and Phipps E. and Utke J. and Walther A.},
 * series      = {Lecture Notes in Computational Science and Engineering},
 * volume      = {73},
 * date        = {2012-05-15},
 * doi         = {10.1007/978-3-642-30023-3_28},
 * eprint      = {1205.3506v1},
 * eprintclass = {cs.MS},
 * eprinttype  = {arXiv}
 * }
 * @endcode
 *
 * 正向和反向模式的Sacado数字的实现是相当复杂的。从Trilinos
 * 12.12开始，数学运算的实现涉及大量的预处理指令和宏编程。因此，代码可能难以理解，而且这些类也不存在有意义的配套文档。因此，理解这些数字的原理实现的有用资源可以在<a
 * href="https://trilinos.org/docs/dev/packages/sacado/doc/html/classSacado_1_1Fad_1_1SimpleFad.html">this
 * link for the Sacado::Fad::SimpleFad
 * class</a>中找到，它概述了一个不使用表达式模板的正向模式自动差分数字的参考（尽管据说效率很低）实现。虽然没有明确说明，但
 * Sacado::Fad::SimpleFad 类似乎是按照双数的精神实现的）。
 *
 *
 * @subsection  auto_diff_1_2 自动微分是如何集成到deal.II的？
 * 由于每个自动分化库的接口都有很大的不同，所以在不久的将来会建立一个统一的内部接口来连接每个数字。目标将是让一些驱动类（提供核心功能，以后将在下一节介绍）有一个一致的机制与不同的自动区分库互动。具体来说，它们需要能够正确地初始化和最终确定将被解释为公式的因变量和自变量的数据。
 * 实现支持的自动差分数字的接口的文件摘要如下。
 *
 *
 * - ad_drivers.h:提供作为内部支持的自动区分库接口的驱动的类。这些在内部被用作我们提供的帮助类的中介。
 *
 *
 * - ad_helpers.h:提供了一系列的类来帮助在一些不同的情况下进行自动区分。这些详见  @ref auto_diff_1_3  。
 *
 *
 * - ad_number_types.h:引入了一个枚举（称为类型代码），用于驱动类将支持的可自动区分的数字组合。   下面将讨论使用这种有点限制性的机制的理由。
 *
 *
 * - ad_number_traits.h:声明一些为每个自动区分库和/或数字类型专用的内部类。这些随后被用来通过NumberTraits和ADNumberTraits类提供一个统一的接口，这些类在整个驱动中被广泛使用。我们还提供了一些机制来轻松查询这些数字的选择属性，即一些类型特征。
 *
 *
 * - adolc_math.h:对ADOL-C数学运算的扩展，使这些数字在整个库中得到一致使用。
 *
 *
 * - adolc_number_types.h:实现内部类，定义我们如何使用ADOL-C数字。
 *
 *
 * - adolc_product_types.h:定义了一些乘积和标量类型，允许与Tensor和SymmetricTensor类一起使用ADOL-C数字。
 *
 *
 * - sacado_math.h:对Sacado数学运算的扩展，使这些数字在整个库中被一致使用。
 *
 *
 * - sacado_number_types.h:实现内部类，定义我们如何使用支持的Sacado数字。
 *
 *
 * - sacado_product_types.h:定义了一些乘积和标量类型，允许与Tensor和SymmetricTensor类一起使用支持的Sacado数。
 * 通过对每个支持的数字类型使用类型代码，我们人为地限制了可以在库中使用的可自动区分的数字类型。这种设计选择是由于确保每个数字类型被正确初始化，以及所有嵌套（模板化）类型的组合对库中执行的所有操作保持有效并不是一件小事。此外，库中还有一些冗长的函数，这些函数是为支持的数字类型而实例化的，并且有内部检查，只有在使用库所知道的可自动区分的数字时才会满足。这再次确保了所有计算的完整性得到维护。最后，使用一个简单的枚举作为类模板参数，最终使得在生产代码中切换使用的类型变得非常容易，几乎不需要对用户代码进行进一步的修改。
 * @subsubsection  auto_diff_1_3 自动分化库的用户接口
 * deal.II库为我们支持的自动区分库提供了一个统一的接口。到目前为止，已经为以下情况开发了帮助类。
 *
 *
 * - 设计用于在正交点水平（或任何一般连续点）操作的类。
 *
 *
 *
 * -  Differentiation::AD::ScalarFunction:  %一个标量值函数的微分。       一个典型的用途是直接从应变能量函数中开发构成法。在 step-71 中给出了这种确切使用情况的一个例子。
 *
 *
 *
 *
 * -  Differentiation::AD::VectorFunction:  %对矢量值函数进行微分。       这可用于构成法的运动变量的线性化，或协助解决局部内部变量的进化方程。
 *
 *
 * - 设计用于在细胞水平上操作的类。
 *
 *
 *
 *
 * -  Differentiation::AD::EnergyFunctional:  %标量值能量函数的微分，如可能产生于变分公式。这个类别的一个例子是在  step-72  中使用的。
 *
 *
 *
 *
 * -  Differentiation::AD::ResidualLinearization:  矢量值的有限元残差的差异化，导致其一致的线性化。   step-72 也提供了一个关于如何使用这个类的演示。
 * 当然，用户也可以自己管理初始化和导数计算。
 * 关于如何使用ADOL-C的最新例子可以在以下文献中找到
 *
 *
 * - 他们的<a href="https://projects.coin-or.org/ADOL-C/browser/trunk/ADOL-C/doc/adolc-manual.pdf?format=raw">user manual</a>。
 *
 *
 * - 他们的<a href="https://gitlab.com/adol-c/adol-c/tree/master/ADOL-C/examples">development repository</a>，和
 *
 *
 * - 我们的<a href="https://github.com/dealii/dealii/tree/master/tests/adolc">test-suite</a>。
 * 而对于萨卡多来说，说明性的例子可以在下面找到
 *
 *
 * - 他们的<a href="https://github.com/trilinos/Trilinos/tree/master/packages/sacado/example">development repository</a>。
 *
 *
 * - 一个<a href="https://github.com/dealii/code-gallery/tree/master/Quasi_static_Finite_strain_Compressible_Elasticity">code-gallery example</a>，和
 *
 *
 * - 我们的<a href="https://github.com/dealii/dealii/tree/master/tests/sacado">test-suite</a>。
 *
 *
 * @section  symb_diff_1 符号表达式和微分法 <a
 * href="https://en.wikipedia.org/wiki/Symbolic_differentiation">Symbolic
 * differentiation</a>就其设计和使用而言，与自动微分完全不同。任何符号库的基础都是一个计算机代数系统（CAS），它实现了一种语言和一系列的算法来操作符号（或
 * "类字符串"）表达。从哲学的角度来看，这与手工进行代数运算的方式最为相似。
 * 为了帮助更好地区分符号微分和自动微分等数值方法，让我们考虑一个非常简单的例子。假设函数
 * $f(x,y) = [2x+1]^{y}$  ，其中  $x$  和  $y$
 * 是相互独立的变量。通过应用连锁法则，这个函数的导数只是
 * $\dfrac{d f(x,y)}{d x} = 2y[2x+1]^{y-1}$  和  $\dfrac{d f(x,y)}{d y} =
 * [2x+1]^{y} \ln(2x+1)$
 * 。这些正是你在定义符号变量`x`和`y`，定义符号表达式`f
 * = pow(2x+1, y)`并计算导数`diff(f, x)`和`diff(f,
 * y)`后从CAS得到的结果。在这一点上，没有假设`x'和`y'代表什么；它们以后可能被解释为普通（标量）数、复数或其他一些幂函数和自然对数函数被定义好的东西。很明显，这意味着也没有关于哪一点来评估表达式或其导数的假设。我们可以很容易地将
 * $\dfrac{d f(x, y)}{d x}$ 的表达式在 $x=1, y=2.5$
 * 处评估，然后再在不重新计算导数表达式本身的情况下，在
 * $x=3.25, y=-6$
 * 处评估。事实上，任何符号变量或表达式的解释，以及变量之间的相互依赖关系，都可以在它们的操作过程中的任何时候被定义或重新定义；这导致了计算的灵活性，这是自动差分所不能比拟的。例如，我们可以执行永久替换
 * $g(x) = \dfrac{d f(x, y)}{d x} \vert_{y=1}$ ，然后针对 $x$
 * 的几个不同值重新计算 $g(x)$
 * 。我们也可以事后表达`x`和`y`之间的相互依赖关系，如 $y
 * \rightarrow y(x) := 2x$
 * 。对于这种情况，这意味着最初计算的导数 $\dfrac{d f(x,
 * y)}{d x} \rightarrow \dfrac{\partial f(x, y(x))}{\partial x} = 2y(x)
 * [2x+1]^{y(x)-1} = 4x[2x+1]^{2x-1}$ 和 $\dfrac{d f(x, y)}{d y} \rightarrow
 * \dfrac{\partial f(x, y(x))}{\partial y} = [2x+1]^{y(x)} \ln(2x+1) =
 * [2x+1]^{2x} \ln(2x+1)$
 * 真正代表了部分导数而不是总导数。当然，如果在计算导数
 * $\dfrac{d f(x, y(x))}{d x}$ 和 $\dfrac{d f(x, y(x))}{d y}$
 * 之前明确定义了这样的相互依赖关系，那么这可能对应于总导数（这也是自动差分在这个例子中唯一能够实现的结果）。
 * 由于复杂的CAS构成了符号操作的基础，操作类型不一定只限于微分，而是可能跨越与离散微分计算、纯数学主题等相关的操作范围。<a
 * href="https://www.sympy.org/en/index.html">SymPy</a>库的文档给出了大量的例子，突出了一个成熟的CAS的能力。通过
 * Differentiation::SD::Expression 类和 Differentiation::SD
 * 命名空间中的相关函数，我们为高性能的<a
 * href="https://github.com/symengine/symengine">SymEngine</a>符号操作库提供了一个封装，它具有丰富的运算符重载和一致的界面，使其易于使用且
 * "自然"。事实上，在许多情况下，这个类可以作为算术类型的
 * "滴入
 * "替代品，将操作从数字性质转变为符号性质；当类在底层数字类型上被模板化时，这就变得特别容易。由于专注于PDEs的数值模拟，在deal.II中暴露的CAS功能主要集中在符号表达式的创建、操作和微分上。
 * 对SymEngine功能的方便包装主要集中在仅涉及基于字典（即让人联想到
 * "基于字符串
 * "的东西）的操作。尽管SymEngine以一种有效的方式执行这些操作，但它们仍然是已知的计算昂贵的，特别是当这些操作是在大的表达式上执行时。因此，应该预计到在生产代码中使用时，执行微分、符号替换等
 * @b
 * 的代码部分的性能可能是一个限制因素。因此，deal.II提供了一个接口，通过
 * @p BatchOptimizer
 * 类（本身经常利用SymEngine提供的功能）加速对冗长的符号表达的评估。特别是，
 * @p BatchOptimizer
 * 同时使用常见的子表达式消除（CSE）等方法来优化符号表达式的集合，以及通过使用自定义生成的
 * `std::function` 或使用LLVM
 * JIT编译器来生成高性能代码路径来评估这些表达式。
 * Differentiation::SD::BatchOptimizer 类的用法在 step-71
 * 中有所例证。
 * 作为最后的说明，必须认识到deal.II目前实现的支持符号库的接口还有很大的缺陷。目前实现的功能水平有效地限制了符号代数在传统使用情况下的使用（即标量和张量代数，因为可能对定义构成关系或复杂函数作为边界条件或源项的应用有用）。事实上，
 * step-71
 * 展示了如何使用它来实现具有挑战性的构成模型。在未来，我们还将本着与
 * Differentiation::AD
 * 命名空间相同的精神，实现类来协助执行装配操作。
 * 实现支持的符号可微调数的接口的文件摘要如下。
 *
 *
 * - symengine_math.h:数学运算的实现，使实现符号表达式的类在整个库和用户代码中的使用一致。   它为标准命名空间中的许多数学函数提供对应的定义。
 *
 *
 * - symengine_number_traits.h:提供一些机制来轻松查询符号数的选择属性，即一些类型特征。
 *
 *
 * - symengine_number_types.h: Differentiation::SD::Expression 类的实现，可用于表示标量符号变量、标量符号表达式等。   这个表达式类被赋予了一整套重载的运算符，用于SymEngine库所支持的所有数学和逻辑运算，并且被认为在数值建模的背景下是有用的。
 *
 *
 * - symengine_optimizer.h: Differentiation::SD::BatchOptimizer 类的实现，可用于使用各种技术加速（在某些情况下，显著）符号表达式的评估。
 *
 *
 * - symengine_product_types.h:定义了一些积和标量类型，允许与Tensor和SymmetricTensor类一起使用符号表达式。
 *
 *
 * - symengine_scalar_operations.h:定义了许多可以对标量符号表达式或变量进行的操作。   这包括（但不限于）创建标量符号，对标量进行微分，以及在标量表达式中进行符号替换。
 *
 *
 * - symengine_tensor_operations.h:定义了可以对符号表达式或变量的张量进行的众多操作。   这包括（但不限于）创建符号的张量，对符号的张量进行微分，对符号的张量进行微分，以及在张量表达式中进行符号替换。
 *
 *
 * - symengine_types.h。为一些在符号计算范围内常用的类型提供别名。
 *
 *
 * - symengine_utilities.h:提供一些在符号计算范围内有用的实用函数。
 *
 */


