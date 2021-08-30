/**
@page step_1 The step-1 tutorial program
@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Aboutthetutorial"> About the tutorial </a>
        <li><a href="#Videolecturesontutorialprograms"> Video lectures on tutorial programs </a>
        <li><a href="#Whatthisprogramdoes"> What this program does </a>
        <li><a href="#Aboutscientificcomputingingeneral"> About scientific computing in general </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Creatingthefirstmesh">Creating the first mesh</a>
        <li><a href="#Creatingthesecondmesh">Creating the second mesh</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Differentadaptiverefinementstrategies"> Different adaptive refinement strategies </a>
        <li><a href="#Differentgeometries"> Different geometries </a>
        <li><a href="#Commentsaboutprogramminganddebugging"> Comments about programming and debugging </a>
        <li><a href="#Moreaboutgraphicaloutput"> More about graphical output </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-1/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Aboutthetutorial"></a><h3> About the tutorial </h3>


由于这是第一个教程程序，让我们首先评论一下这个教程和deal.II的其他文档应该如何工作。deal.II的文档基本上有三个不同的层次。

- 该教程。这是一个程序集，展示了deal.II在实践中的应用。它通常不在单个参数的层面上讨论单个函数，而是希望给出事物如何共同工作的大画面。换句话说，它讨论的是 "概念"：什么是deal.II的构件，它们如何在有限元程序中一起使用。

- 该手册。这是deal.II中每一个类和每一个（成员）函数的文档。例如，如果你点击本页顶部的 "主页 "或 "类 "选项卡，你就可以看到。在这里你可以查到 Triangulation::create_triangulation_compatibility 的第二个参数是什么意思，这只是一个略显晦涩的例子。当你知道你要做什么，但忘记了函数到底是怎么命名的，它的参数是什么，或者它的返回值是什么时，你就需要这种级别的文档。请注意，当你读完教程并点击任何一个类或函数名称时，你也会进入手册，也就是说，当你需要对某个函数或类进行更详细的描述时，教程包含了大量进入手册的链接。另一方面，手册并不是学习deal.II的好地方，因为它只给你一个微观的观点，而没有告诉你一个函数是如何融入大局的。

- 模块。这些是一起工作或具有相关功能的类和函数组。如果你点击本页顶部的 "模块 "标签，你就会进入一个列出许多此类组的页面。每个模块都讨论了这些类的基本原理；例如， @ref Sparsity 模块讨论了与存储矩阵的稀疏模式有关的各种不同问题。这就是中级水平的文档：它们给你一个特定领域的概述。例如，当你想知道存在哪些有限元类时，你会看一下 @ref fe 模块。当然，这些模块也与手册（有时也与教程）有交叉链接；如果你点击一个类的名字，比如说三角法，如果你想了解更多关于这个类的背景，也会在类名的右上方得到一个指向这个类所属模块的链接。

让我们回到教程中来，因为你正在看的是它的第一个程序（或 "步骤"）。每个教程的程序都被细分为以下几个部分。<ol>  <li>  <b>Introduction:</b> 这是讨论程序的作用，包括数学模型，以及与以前的教程程序相比有哪些新的编程技术。     <li>  <b>The commented program:</b> 广泛地记录了源代码的清单。在这里，我们经常记录个别行或代码块，并讨论它们做什么，如何做，以及为什么。评论中经常提到介绍，也就是说，你必须先了解<i>what</i>程序想要达到的目标（介绍中讨论的目标），然后才能了解<i>how</i>它想要达到的目标。     <li>  <b>Results:</b> 程序的输出，包括注释和解释。这一部分也经常有一个小节，给出如何在不同方向上扩展程序的建议；在早期的程序中，这是为了给你提供小实验的方向，旨在使你熟悉deal.II，而在后来的程序中，更多的是关于如何使用更高级的数值技术。     <li>  <b>The plain program:</b> 剥去所有注释的源代码。如果你想看到代码的 "全貌"，这很有用，因为程序的注释版本中间有很多文字，往往很难在屏幕上一次看到单个函数的全部代码。   </ol> 

教程不仅意味着是静态文档，而且你应该玩玩它们。为此，进入 <code>example/step-1</code> 目录（或任何你感兴趣的教程的编号），然后输入

@code
  cmake .
  make
  make run
@endcode

第一条命令设置了描述本教程程序所依赖的包含文件、如何编译以及如何运行的文件。这条命令应该也能找到已安装的deal.II库，这些库是在你按照<a href="../../readme.html" target="body">README</a>文件中描述的方式编译和安装一切时产生的。如果这个命令不能找到deal.II库，那么你需要用命令提供安装的路径

@code
  cmake -DDEAL_II_DIR=/path/to/installed/deal.II .
@endcode

而不是。

上述命令中的第二条将源代码编译成可执行文件，而最后一条则是执行它（严格来说，如果可执行文件还不存在， <code>make run</code> 也会编译代码，所以如果你想的话，你可以跳过第二条命令）。这就是运行代码和产生输出所需的全部内容，在教程程序的 "结果 "部分讨论。这个顺序需要在你想玩的所有教程目录中重复。

当学习这个库时，你需要玩玩它，看看会发生什么。为此，用你喜欢的编辑器打开 <code>example/step-1/step-1.cc</code> 的源文件，并以某种方式进行修改，保存它并按上述方式运行它。在这个程序的结果部分的末尾给出了一些可能的修改建议，在那里我们还提供了一些其他有用信息的链接。




<a name="Videolecturesontutorialprograms"></a><h3> Video lectures on tutorial programs </h3>


在关于deal.II和计算科学的<a
href="http://www.math.colostate.edu/~bangerth/videos.html">Wolfgang
Bangerth video lectures</a>中也讨论和演示了这个和其他几个教程程序。特别是，你可以看到他为运行这个和其他程序所执行的步骤，你会对可以用来处理deal.II的工具有一个更好的了解。特别是，第2和第4讲概述了deal.II和任何有限元代码的构建模块。(  @dealiiVideoLectureSeeAlso{2,4}) 

如果你还不熟悉使用Linux和在命令行上运行东西，你可能会有兴趣观看讲座2.9和2.91。(  @dealiiVideoLectureSeeAlso{2.9,2.91}) These give overviews over the command 行和关于编译程序时发生的事情，分别。

请注意，deal.II正在积极开发，在开发过程中，我们偶尔会对这些视频讲座中仍然引用的函数或类进行重新命名或废弃。  例如，视频讲座5中的步骤1代码使用了一个HyperShellBoundary类，后来被SphericalManifold类取代。此外，从deal.II 9.0版本开始， GridGenerator::hyper_shell() 现在自动将SphericalManifold附加到Triangulation上。否则，讲座材料的其余部分都是相关的。

<a name="Whatthisprogramdoes"></a><h3> What this program does </h3>


让我们回到步骤1，即当前的程序。在这第一个例子中，我们实际上并没有做很多事情，而是展示了两种技术：生成三角形对象的语法是什么，以及所有单元格上简单循环的一些元素。我们创建了两个网格，一个是有规律地细化的正方形（不是很刺激，但对于一些问题来说是常见的起始网格），还有一个是更多的几何尝试：一个环形域，向内边缘细化。通过这些，你将了解到每一个有限元程序都必须有的三样东西。一个用于网格的Triangulation类型的对象；对GridGenerator函数的调用以生成网格；以及涉及迭代器的所有单元的循环（迭代器是指针的泛化，在C++标准库中经常使用；在deal.II的背景下， @ref Iterators 模块谈到了它们）。

该程序在其他方面足够小，不需要大量的介绍。

 @dealiiVideoLecture{5,6} 




<a name="Aboutscientificcomputingingeneral"></a><h3> About scientific computing in general </h3>


如果你正在阅读这个教程程序，很可能你有兴趣继续使用deal.II来完成你自己的项目。因此，你即将开始一个使用大规模科学计算库的编程练习。除非你已经是大规模编程方法的资深用户，否则这对你来说可能是一个新的领域；伴随着所有的新规则，比如你将不得不处理别人写的代码，你可能不得不考虑记录自己的代码，因为你可能在一年后不记得它到底在做什么（或者因为别人也会使用它），或者想出一些方法来测试你的程序是否在做正确的事情。这些都不是我们通常训练数学家、工程师或科学家的东西，但当你开始编写超过几百行的软件时，这些就很重要了。请记住。制作软件并不等同于仅仅编写代码。

为了使你在这一旅程中生活得更轻松，让我们指出一些资源，这些资源在你开始任何大规模的编程之前是值得浏览的。

- 在<a href="https://github.com/dealii/dealii/wiki/Frequently-Asked-Questions">
  deal.II FAQ</a>中，有大量关于deal.II特定方面问题的答案，但也有一些更普遍的问题，如 "我如何调试科学计算代码？"或 "我能否训练自己写出错误更少的代码？"。

- 你将从成为一个更好的程序员中受益。为此，一个很好的资源是Steve McConnell的[Code Complete](https://en.wikipedia.org/wiki/Code_Complete)  @cite CodeComplete  。这本书已经有几年的历史了，最后一版是在2004年出版的，但它作为良好的编程实践指南的吸引力丝毫不减，一些主要的开发者把它作为他们研究小组每一代成员的集体阅读项目。

- <a href="http://software-carpentry.org/">Software Carpentry project</a>，提供了处理软件的许多重要主题的介绍，如版本控制、make文件、测试等。它是专门为科学家和工程师编写的，而不是为计算机科学家编写的，并以简短的、实用的课程为重点。

- <a href="https://bssw.io/">Better Scientific Software
  project</a>有大量的资源（和有趣的博文），涵盖了编写科学软件的许多方面。

- <a href="https://ideas-productivity.org/">IDEAS
  project</a>也有关于软件开发的资源，特别是用于并行计算。在该网站的 "活动 "部分有录制的教程和网络研讨会，涉及许多有趣的主题。

- 一篇关于<a href="http://arxiv.org/abs/1210.0530">Best
  Practices for Scientific Computing</a>的文章，介绍了许多方法，通过这些方法，你可以确保你是一个高效的程序员，写出的程序可以正常工作。

作为一个一般性建议。如果你期望在未来花几天时间来编写软件，请帮你自己的忙，学习能够使你的生活更有效率的工具，特别是调试器和集成开发环境。(  @dealiiVideoLectureSeeAlso{7,8,8.01,25})  你会发现，通过提高工作效率，你很快就会把学习这些工具的时间拿回来几倍!上面提到的几个视频讲座展示了如何使用集成开发环境或调试器等工具。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * 库中最基本的类是Triangulation类，它在这里声明。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/tria.h> 
 * 
 * @endcode
 * 
 * 这里有一些生成标准网格的函数。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * @endcode
 * 
 * 输出各种图形格式的网格。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_out.h> 
 * 
 * @endcode
 * 
 * 这对于C++输出来说是需要的。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * 
 * @endcode
 * 
 * 这是对 `std::sqrt` 和 `std::fabs` 函数声明的说明。
 * 

 * 
 * 
 * @code
 * #include <cmath> 
 * 
 * @endcode
 * 
 * 导入deal.II的最后一步是这样的。所有deal.II的函数和类都在一个命名空间 <code>dealii</code> 中，以确保它们不会与你可能想和deal.II一起使用的其他库的符号发生冲突。我们可以在使用这些函数和类时，在每个名字前加上 <code>dealii::</code> 的前缀，但这很快就会变得繁琐和令人厌烦。相反，我们只是简单地导入整个deal.II的名字空间，以供一般使用。
 * 

 * 
 * 
 * @code
 * using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Creatingthefirstmesh"></a> 
 * <h3>Creating the first mesh</h3>
 * 

 * 
 * 在下面的第一个函数中，我们简单地使用单位方格作为域，并从中产生一个全局细化网格。
 * 

 * 
 * 
 * @code
 * void first_grid() 
 * { 
 * 
 * @endcode
 * 
 * 首先要做的是为二维域的三角化定义一个对象。
 * 

 * 
 * 
 * @code
 *   Triangulation<2> triangulation; 
 * 
 * @endcode
 * 
 * 在这里和下面的许多情况下，类名后面的字符串"<2>"表示这是一个在两个空间维度上工作的对象。同样，也有一些三角形类的版本是在一个（"<1>"）和三个（"<3>"）空间维度上工作的。这种工作方式是通过一些模板魔法实现的，我们将在后面的示例程序中详细研究；在那里，我们也将看到如何以一种基本独立于维度的方式编写程序。
 * 

 * 
 * 接下来，我们要用一个正方形领域的单个单元来填充三角结构。三角形被细化了四次，总共得到 $4^4=256$ 个单元。
 * 

 * 
 * 
 * @code
 *   GridGenerator::hyper_cube(triangulation); 
 *   triangulation.refine_global(4); 
 * 
 * @endcode
 * 
 * 现在我们要将网格的图形表示写到输出文件中。deal.II的GridOut类可以用多种不同的输出格式来实现；在这里，我们选择可扩展矢量图（SVG）格式，你可以用你选择的网络浏览器来进行可视化。
 * 

 * 
 * 
 * @code
 *   std::ofstream out("grid-1.svg"); 
 *   GridOut       grid_out; 
 *   grid_out.write_svg(triangulation, out); 
 *   std::cout << "Grid written to grid-1.svg" << std::endl; 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Creatingthesecondmesh"></a> 
 * <h3>Creating the second mesh</h3>
 * 

 * 
 * 下面第二个函数中的网格略微复杂一些，因为我们使用了一个环形域，并对结果进行了一次全局细化。
 * 

 * 
 * 
 * @code
 * void second_grid() 
 * { 
 * 
 * @endcode
 * 
 * 我们再次开始定义一个二维域的三角化对象。
 * 

 * 
 * 
 * @code
 *   Triangulation<2> triangulation; 
 * 
 * @endcode
 * 
 * 然后我们用一个环形域来填充它。环的中心应是点(1,0)，内半径和外半径应是0.5和1。圆周单元的数量可以由这个函数自动调整，但我们选择在最后一个参数中明确设置为10。
 * 

 * 
 * 
 * @code
 *   const Point<2> center(1, 0); 
 *   const double   inner_radius = 0.5, outer_radius = 1.0; 
 *   GridGenerator::hyper_shell( 
 *     triangulation, center, inner_radius, outer_radius, 10); 
 * 
 * @endcode
 * 
 * 默认情况下，三角测量假定所有边界都是直线，所有单元都是双线性四边形或三线性六边形，并且它们是由粗略网格（我们刚刚创建的）的单元定义的。除非我们做一些特别的事情，否则当需要引入新的点时，域被假定为由粗网格的直线划定，而新的点将简单地位于周围的中间。然而，在这里，我们知道领域是弯曲的，我们想让三角法根据底层的几何形状来放置新的点。幸运的是，一些优秀的灵魂实现了一个描述球状域的对象，而环是球状域的一个部分；它只需要环的中心，并自动计算出如何指示三角计算在哪里放置新的点。这在deal.II中的工作方式是，你用一个通常被称为 "流形指标 "的数字来标记你想要弯曲的三角形部分，然后告诉三角形在所有有这个流形指标的地方使用一个特定的 "流形对象"。具体如何操作在此并不重要（你可以在 step-53 和 @ref manifold 中阅读）。GridGenerator中的函数在大多数情况下为我们处理这个问题：它们将正确的流形附加到一个域上，这样当三角形被细化时，新的单元就会被放置在正确的位置上。在目前的情况下， GridGenerator::hyper_shell 为所有的单元格附加了一个球形流形：这将导致单元格在球面坐标的计算下被细化（因此新的单元格的边缘要么是径向的，要么是位于原点周围的同心圆）。
 * 

 * 
 * 默认情况下（即对于手工创建的三角图或未调用GridGenerator函数（如 GridGenerator::hyper_shell 或 GridGenerator::hyper_ball), ），三角图的所有单元格和面都将其manifold_id设置为 numbers::flat_manifold_id, ，如果您想要一个产生直线边缘的流形，这是默认的，但您可以为个别单元格和面改变这个数字。在这种情况下，因此与数字0相关的曲面流形将不适用于那些流形指标为非零的部分，但其他流形描述对象可以与这些非零指标相关联。如果没有流形描述与特定的流形指标相关联，则暗示产生直角边缘的流形。(流形指标是一个略微复杂的话题；如果你对这里到底发生了什么感到困惑，你可能想看看 @ref GlossManifoldIndicator "关于这个话题的词汇表条目")。既然 GridGenerator::hyper_shell 选择的默认值是合理的，我们就不去管它。
 * 

 * 
 * 为了演示如何在所有单元格上写一个循环，我们将分五个步骤向域的内圈细化网格。
 * 

 * 
 * 
 * @code
 *   for (unsigned int step = 0; step < 5; ++step) 
 *     { 
 * 
 * @endcode
 * 
 * 接下来，我们需要对三角形的活动单元进行循环。你可以把三角形看作一个单元格的集合。如果它是一个数组，你只需要得到一个指针，用操作符`++`从一个元素递增到下一个元素。三角形的单元不是作为一个简单的数组来存储的，但是<i>iterator</i>的概念将指针的工作方式概括为任意的对象集合（更多信息见<a href= "http:en.wikipedia.org/wiki/Iterator#C.2B.2B">wikipedia</a>）。通常情况下，C++中的任何容器类型都会返回一个迭代器，指向集合的开始，方法称为`begin'，而迭代器则指向集合结束后的1，方法称为`end'。我们可以用操作符`++it`来增加一个迭代器`it`，用`*it`来解除引用以获得底层数据，并通过比较`it != collection.end()`来检查我们是否完成。
 * 

 * 
 * 第二个重要的部分是我们只需要活动单元。活动单元是那些没有被进一步细化的单元，也是唯一可以被标记为进一步细化的单元。deal.II提供了迭代器类别，允许我们在<i>all</i>单元（包括活动单元的父单元）或只在活动单元上迭代。因为我们要的是后者，所以我们需要调用方法 Triangulation::active_cell_iterators().  。
 * 

 * 
 * 把所有这些放在一起，我们可以用
 * <div class=CodeFragmentInTutorialComment>
 * @code{.cpp}
 *      for (auto it = triangulation.active_cell_iterators().begin();
 *           it != triangulation.active_cell_iterators().end();
 *           ++it)
 *        {
 *          auto cell = *it;
 *  //Then a miracle occurs...
 *        }
 *  @endcode
 * </div>
 * 在一个三角形的所有活动单元上循环。 在这个循环的初始化器中，我们使用了`auto`关键字作为迭代器`it`的类型。`auto`关键字意味着被声明的对象的类型将从上下文中推断出来。当实际的类型名称很长，甚至可能是多余的时候，这个关键字很有用。如果你不确定类型是什么，想查一下结果支持什么操作，你可以去看方法的文档  Triangulation::active_cell_iterators().  在这个例子中，`it`的类型是  `Triangulation::active_cell_iterator`.  
 * 

 * 
 * 虽然`auto`关键字可以让我们不用输入长长的数据类型名称，但我们仍然要输入大量冗余的关于开始和结束迭代器以及如何递增的声明。与其这样，我们不如使用<a href="http:en.cppreference.com/w/cpp/language/range-for">range-based for loops</a>，它将上面显示的所有语法包成一个更短的形式。
 * 

 * 
 * 
 * @code
 *       for (auto &cell : triangulation.active_cell_iterators()) 
 *         { 
 * @endcode
 * 
 * @note  关于deal.II中使用的迭代器类的更多信息，见 @ref Iterators ，关于基于范围的for循环和`auto`关键字的更多信息，见 @ref CPP11 。
 * 

 * 
 * 接下来，我们在单元格的所有顶点上循环。为此，我们查询一个顶点索引的迭代器（在2D中，这是一个包含元素`{0,1,2,3}`的数组，但是由于`cell->vertex_indices()`知道单元格所处的维度，因此返回的数组在所有维度上都是正确的，这使得无论我们在2D还是3D中运行这段代码都是正确的，也就是说，它实现了 "维度无关的编程" - 我们将在  step-4  中讨论一个重要部分）。
 * 

 * 
 * 
 * @code
 *           for (const auto v : cell->vertex_indices()) 
 *             { 
 * 
 * @endcode
 * 
 * 如果这个单元格位于内边界，那么它至少有一个顶点必须位于内环上，因此与中心的径向距离正好是0.5，达到浮点精度。所以我们计算这个距离，如果我们发现一个顶点具有这个属性，我们就标记这个单元，以便以后进行细化。然后我们也可以打破所有顶点的循环，转到下一个单元。        因为离中心的距离是以浮点数计算的，所以我们必须期望我们所计算的东西只能精确到[round-off](https:en.wikipedia.org/wiki/Round-off_error)以内。因此，我们永远不能指望通过平等的方式来比较距离和内半径。诸如 "if (distance_from_center == inner_radius) "这样的语句将会失败，除非我们运气特别好。相反，我们需要以一定的容忍度进行比较，通常的方法是写成`if  (std::abs(distance_from_center  。
 * 

 * 
 * - inner_radius) <= tolerance)`，其中`tolerance'是比四舍五入大的某个小数字。问题是如何选择它。我们可以直接选择，比如说，`1e-10'，但这只适合于我们比较的对象是大小为1的情况。如果我们创建了一个单元大小为`1e+10'的网格，那么`1e-10'将远远低于四舍五入，就像以前一样，只有在我们特别幸运的情况下，比较才会成功。相反，使公差*相对于被比较对象的典型 "比例 "几乎总是有用的。在这里，"尺度 "是指内半径，或者是细胞的直径。我们选择前者，并将公差设置为 $10^{-6}$ 倍环形物的内半径。
 * 

 * 
 * 
 * @code
 *               const double distance_from_center = 
 *                 center.distance(cell->vertex(v)); 
 * 
 *               if (std::fabs(distance_from_center - inner_radius) <= 
 *                   1e-6 * inner_radius) 
 *                 { 
 *                   cell->set_refine_flag(); 
 *                   break; 
 *                 } 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 现在我们已经标记了所有我们想要细化的单元格，我们让三角化实际做这个细化。这样做的函数的名字很长，因为我们也可以标记单元格进行粗化，该函数一次完成粗化和细化。
 * 

 * 
 * 
 * @code
 *       triangulation.execute_coarsening_and_refinement(); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，在这五次细化迭代之后，我们要再次将得到的网格写入文件，同样是SVG格式。这和上面的工作一样。
 * 

 * 
 * 
 * @code
 *   std::ofstream out("grid-2.svg"); 
 *   GridOut       grid_out; 
 *   grid_out.write_svg(triangulation, out); 
 * 
 *   std::cout << "Grid written to grid-2.svg" << std::endl; 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * 最后是主函数。这里没有什么可做的，只是调用两个子函数，产生两个网格。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   first_grid(); 
 *   second_grid(); 
 * } 
 * 
 * 
 * 
 * @endcode
examples/step-1/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序会产生两个网格的图形（grid-1.svg和grid-2.svg）。你可以用大多数网络浏览器打开它们--在最简单的情况下，只要在文件系统资源管理器中打开当前目录，然后点击文件。如果你喜欢在命令行上工作，你可以用该文件调用你的网络浏览器。`firefox grid-1.svg`，`google-chrome grid-1.svg`，或者任何你的浏览器的名字。如果你这样做，这两个网格应该是这样的。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-1-r9.2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-1.grid-2-r9.2.png" alt="">
    </td>
  </tr>
</table> 

左边那个，嗯，不是很刺激。右边的是&mdash；至少是&mdash；非传统的。这些图片对每个单元的 "细化水平 "进行了颜色编码。一个粗略的网格单元要被细分多少次才能得到给定的单元。在左图中，这是无聊的，因为网格被全局细化了若干次，也就是说，<i>every</i>单元被细化的次数相同。

(虽然第二个网状结构完全是人为捏造的，当然在应用中也不太实用，但令大家惊讶的是，它已经进入了文献：见  @cite Mu05  。显然，它至少对某些事情是有好处的）。)




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="Differentadaptiverefinementstrategies"></a><h4> Different adaptive refinement strategies </h4>


这个程序显然没有太多的功能，但特别是 <code>second_grid</code> 函数有一堆你可以玩弄它的地方。例如，你可以修改我们决定细化哪些单元格的标准。一个例子是把条件改成这样。

@code
      for (auto &cell: triangulation.active_cell_iterators())
        if (cell->center()[1] > 0)
          cell->set_refine_flag ();
@endcode

这将细化所有单元中心的 $y$ 坐标大于零的单元（我们通过解除引用 <code>cell</code> 迭代器调用的 <code>TriaAccessor::center</code> 函数返回一个Point<2>对象；下标 <code>[0]</code> 将得到 $x$ 坐标，下标 <code>[1]</code> 得到 $y$  坐标）。通过查看TriaAccessor提供的函数，你也可以使用更复杂的标准进行细化。

一般来说，你能用`cell->something()`形式的操作做什么，在文档中有点困难，因为`cell`不是一个指针，而是一个迭代器。你可以在单元格上调用的函数可以在`TriaAccessor'类的文档中找到（它的函数也可以在单元格的面或更普遍的、出现在三角形中的各种几何对象上调用），以及`CellAccessor'（它增加了一些专门针对*单元格的函数）。

对整个迭代器概念的更彻底的描述可以在 @ref Iterators 文档模块中找到。




<a name="Differentgeometries"></a><h4> Different geometries </h4>


另一种可能性是生成完全不同几何形状的网格。虽然对于复杂的几何体来说，使用网格生成器获得的网格是没有办法的，但是有大量的几何体，deal.II可以使用GridGenerator命名空间的函数来创建网格。许多这样的几何体（如本例程序中使用的几何体）包含有弯曲面的单元：换句话说，我们希望放置在边界上的新顶点位于一个圆上。deal.II通过Manifold类（以及从它继承的类）处理复杂的几何体；尤其是GridGenerator中对应于非笛卡尔网格的函数（如 GridGenerator::hyper_shell 或 GridGenerator::truncated_cone) 将一个Manifold对象附加到三角网格中应该是曲线的部分（分别为SphericalManifold和CylindricalManifold），并在应该是平面的部分使用另一个Manifold（FlatManifold）。关于这些类的设计理念和接口的描述，请参见Manifold的文档或 @ref manifold "manifold模块"。看看它们提供了什么，看看如何在这样的程序中使用它们。

我们还在第49步中讨论了其他各种创建和操作网格的方法（并描述了附加Manifolds的过程）。




<a name="Commentsaboutprogramminganddebugging"></a><h4> Comments about programming and debugging </h4>


最后，我们对用deal.II修改或编写程序做一个总体的评论。当你开始使用教程程序或你自己的应用程序时，你会发现错误会发生：你的程序会包含一些代码，这些代码要么是立即中止程序，要么是一些简单地导致错误结果的bug。无论哪种情况，你都会发现知道如何使用调试器是非常有帮助的：你可能会通过把调试输出放到你的程序中，编译它，然后运行它来应付一段时间，但最终用调试器寻找错误会更快，更方便，更可靠，因为你不必总是重新编译程序，而且你可以检查变量的值和它们的变化。

与其推迟学习如何使用调试器，直到你真的看不到任何其他方法来发现一个错误，这里是我们将在这个项目中提供的一个建议：尽快学习如何使用调试器。这将是很好的时间投资。( 从顶层的<a
href="http://www.dealii.org/">deal.II webpage</a>链接到的 @dealiiVideoLectureSeeAlso{25}) The deal.II Frequently Asked 问题（FAQ）页面也提供了大量关于调试deal.II程序的提示。




<a name="Moreaboutgraphicaloutput"></a><h4> More about graphical output </h4>


在你的论文或出版物中包含网格往往是有用的。为此，按细化级别对单元格进行颜色编码，并在每个单元格上打印单元格号，可能不是很有用。但这并不意味着一定要这样做 -- GridOut类允许为每种可能的输出格式设置标志（参见GridOutFlags命名空间中的类），以控制网格的具体绘制方式。当然，你也可以选择其他的输出文件格式，如VTK或VTU；这对三维网格特别有用，因为二维格式如SVG并不特别有用，因为它固定了三维物体的特定视角。因此，你可能想探索GridOut类中的其他选项。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-1.cc"
*/
