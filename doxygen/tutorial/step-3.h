/**
@page step_3 The step-3 tutorial program
This tutorial depends on step-2.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thebasicsetupoffiniteelementmethods">The basic set up of finite element methods</a>
        <li><a href="#Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"> Should we multiply by a test function from the left or from the right? </a>
        <li><a href="#Computingthematrixandrighthandsidevector"> Computing the matrix and right hand side vector </a>
        <li><a href="#Abouttheimplementation">About the implementation</a>
        <li><a href="#Anoteontypes"> A note on types </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Manynewincludefiles">Many new include files</a>
        <li><a href="#ThecodeStep3codeclass">The <code>Step3</code> class</a>
      <ul>
        <li><a href="#Step3Step3">Step3::Step3</a>
        <li><a href="#Step3make_grid">Step3::make_grid</a>
        <li><a href="#Step3setup_system">Step3::setup_system</a>
        <li><a href="#Step3assemble_system">Step3::assemble_system</a>
        <li><a href="#Step3solve">Step3::solve</a>
        <li><a href="#Step3output_results">Step3::output_results</a>
        <li><a href="#Step3run">Step3::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
        <li><a href="#UsingHDF5tooutputthesolutionandadditionaldata">Using HDF5 to output the solution and additional data</a>
      <ul>
        <li><a href="#Changingtheoutputtoh5"> Changing the output to .h5</a>
        <li><a href="#Addingthepointvalueandthemeanseeextensionaboveintotheh5file"> Adding the point value and the mean (see extension above) into the .h5 file</a>
      </ul>
        <li><a href="#UsingRandggplot2togenerateplots"> Using R and ggplot2 to generate plots</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-3/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{10} 

<a name="Thebasicsetupoffiniteelementmethods"></a><h3>The basic set up of finite element methods</h3>


这是第一个我们实际使用有限元来计算的例子。我们将解决一个简单的泊松方程，其边界值为零，但右手边非零。

@f{align*}


  -\Delta u &= f \qquad\qquad & \text{in}\ \Omega,
  \\
  u &= 0 \qquad\qquad & \text{on}\ \partial\Omega.


@f}

我们将在正方形 $\Omega=[-1,1]^2$ 上求解这个方程，你已经在步骤1和步骤2中学习了如何生成网格。在这个程序中，我们也将只考虑 $f(\mathbf x)=1$ 这个特殊情况，并在下一个教程程序中再来讨论如何实现更一般的情况，即步骤4。

如果你学过有限元方法的基本知识，你会记得我们需要采取的步骤，用有限维度的近似方法来近似解 $u$ 。具体来说，我们首先需要推导出上述方程的弱形式，通过将方程乘以测试函数 $\varphi$ <i>from the left</i>（我们将在下面回到从左而非从右相乘的原因）并在域 $\Omega$ 上积分得到。

@f{align*}


  -\int_\Omega \varphi \Delta u = \int_\Omega \varphi f.


@f}

这可以通过部件进行整合。

@f{align*}
  \int_\Omega \nabla\varphi \cdot \nabla u


  -
  \int_{\partial\Omega} \varphi \mathbf{n}\cdot \nabla u
   = \int_\Omega \varphi f.


@f}

测试函数 $\varphi$ 必须满足同样的边界条件（用数学术语来说：它需要来自我们寻求解决方案的集合的切线空间），因此在边界上 $\varphi=0$ ，因此我们正在寻找的弱形式为

@f{align*}
  (\nabla\varphi, \nabla u)
   = (\varphi, f),


@f}

其中我们使用了常用的符号  $(a,b)=\int_\Omega a\; b$  。然后，问题要求从适当的空间（这里是空间 $H^1$ ）中找出一个函数 $u$ ，对于该函数，这一声明对于所有测试函数 $\varphi$ 都是真的。

当然，在一般情况下，我们无法在计算机上找到这样的函数，而是寻求一个近似值 $u_h(\mathbf x)=\sum_j U_j \varphi_j(\mathbf
x)$  ，其中 $U_j$ 是我们需要确定的未知膨胀系数（这个问题的 "自由度"）， $\varphi_i(\mathbf x)$ 是我们将使用的有限元形状函数。为了定义这些形状函数，我们需要以下内容。

- 一个用来定义形状函数的网格。你已经看到如何在步骤1和步骤2中生成和操作描述网格的对象。

- 一个描述我们想在参考单元上使用的形状函数的有限元（在deal.II中总是单位间隔 $[0,1]$ 、单位正方形 $[0,1]^2$ 或单位立方体 $[0,1]^3$ ，取决于你在哪个空间维度工作）。在步骤2中，我们已经使用了FE_Q<2>类型的对象，它表示通常的拉格朗日元素，通过对支持点的插值来定义形状函数。最简单的是FE_Q<2>(1)，它使用1度的多项式。在2D中，这些通常被称为<i>bilinear</i>，因为它们在参考单元的两个坐标中都是线性的。(在1d中，它们是<i>linear</i>，在3d中是<i>tri-linear</i>；然而，在deal.II文档中，我们经常不做这种区分，而总是简单地称这些函数为 "线性"。)

- 一个DoFHandler对象，以有限元对象提供的参考单元描述为基础，枚举网格上的所有自由度。你也已经在步骤2中看到了如何做到这一点。

- 一个映射，告诉你如何从参考单元上的有限元类定义的形状函数中获得实数单元上的形状函数。默认情况下，除非你明确说明，否则deal.II将使用（双，三）线性映射，所以在大多数情况下，你不必担心这个步骤。

通过这些步骤，我们现在有一组函数 $\varphi_i$ ，我们可以定义离散问题的弱形式：找到一个函数 $u_h$ ，即找到上面提到的扩展系数 $U_j$ ，以便

@f{align*}
  (\nabla\varphi_i, \nabla u_h)
   = (\varphi_i, f),
   \qquad\qquad
   i=0\ldots N-1.


@f}

请注意，我们在此遵循惯例，即一切从零开始计算，这在C和C++中很常见。如果你插入表示法 $u_h(\mathbf x)=\sum_j U_j
\varphi_j(\mathbf x)$ ，这个方程可以重写为一个线性系统，然后观察到

@f{align*}{
  (\nabla\varphi_i, \nabla u_h)
  &= \left(\nabla\varphi_i, \nabla \Bigl[\sum_j U_j \varphi_j\Bigr]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \left[U_j \varphi_j\right]\right)
\\
  &= \sum_j \left(\nabla\varphi_i, \nabla \varphi_j \right) U_j.


@f}

有了这个，问题就成了。找到一个向量 $U$ ，以便

@f{align*}{
  A U = F,


@f}

其中矩阵 $A$ 和右手边 $F$ 定义为

@f{align*}
  A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j),
  \\
  F_i &= (\varphi_i, f).


@f}






<a name="Shouldwemultiplybyatestfunctionfromtheleftorfromtheright"></a><h3> Should we multiply by a test function from the left or from the right? </h3>


在我们继续描述如何计算这些数量之前，请注意，如果我们从<i>right</i>乘以一个测试函数而不是从左边乘以原方程，那么我们将得到一个形式为的线性系统

@f{align*}
  U^T A = F^T


@f}

有一个行向量  $F^T$  。通过转置这个系统，这当然等同于解决了

@f{align*}
  A^T U = F


@f}

这里与上面的 $A=A^T$ 相同。但一般来说不是，为了避免任何形式的混淆，经验表明，只要养成从左边而不是从右边乘方程的习惯（正如数学文献中经常做的那样），就可以避免一类常见的错误，因为在比较理论和实现时，矩阵会自动正确，不需要转置。本教程的第一个例子见第9步，我们有一个非对称的双线性方程，对于这个方程，我们从右面还是从左面相乘是有区别的。




<a name="Computingthematrixandrighthandsidevector"></a><h3> Computing the matrix and right hand side vector </h3>


现在我们知道我们需要什么（即：持有矩阵和向量的对象，以及计算 $A_{ij},F_i$ 的方法），我们可以看看需要什么来实现这一点。

-  $A$ 的对象是SparseMatrix类型，而 $U$ 和 $F$ 的对象则是Vector类型。我们将在下面的程序中看到哪些类是用来解决线性系统的。

- 我们需要一种方法来形成积分。在有限元方法中，最常见的是使用正交法，也就是说，积分被每个单元上的一组点的加权和所取代。也就是说，我们首先将 $\Omega$ 的积分分成所有单元的积分，@f{align*}
    A_{ij} &= (\nabla\varphi_i, \nabla \varphi_j)
    = \sum_{K \in {\mathbb T}} \int_K \nabla\varphi_i \cdot \nabla \varphi_j,
    \\
    F_i &= (\varphi_i, f)
    = \sum_{K \in {\mathbb T}} \int_K \varphi_i f,
  @f}

  然后用正交法对每个单元的贡献进行近似。   @f{align*}
    A^K_{ij} &=
    \int_K \nabla\varphi_i \cdot \nabla \varphi_j
    \approx
    \sum_q \nabla\varphi_i(\mathbf x^K_q) \cdot \nabla
    \varphi_j(\mathbf x^K_q) w_q^K,
    \\
    F^K_i &=
    \int_K \varphi_i f
    \approx
    \sum_q \varphi_i(\mathbf x^K_q) f(\mathbf x^K_q) w^K_q,
  @f}

  其中 $\mathbf x^K_q$ 是 $q$ 单元上的第三个正交点 $K$ ， $w^K_q$ 是 $q$ 的正交权。这样做需要有不同的部分，接下来我们将依次讨论它们。

- 首先，我们需要一种方法来描述正交点的位置  $\mathbf x_q^K$  和它们的权重  $w^K_q$  。它们通常以与形状函数相同的方式从参考单元映射出来，即隐含地使用MappingQ1类，或者，如果你明确地说，通过从Mapping派生的其他类之一。参考单元上的位置和权重由派生自正交基类的对象来描述。通常，人们选择一个正交公式（即一组点和权重），使正交正好等于矩阵中的积分；这可以实现，因为积分中的所有因子都是多项式，由高斯正交公式完成，在QGauss类中实现。

- 然后我们需要一些东西来帮助我们在 $K$ 单元上评估 $\varphi_i(\mathbf x^K_q)$ 。这就是FEValues类的作用：它需要一个有限元对象来描述参考单元上的 $\varphi$ ，一个正交对象来描述正交点和权重，以及一个映射对象（或隐含地采用MappingQ1类），并在位于 $K$ 的正交点上提供形状函数的值和导数，以及积分所需的各种其他信息。

FEValues确实是装配过程中的核心类。你可以这样看待它。FiniteElement和派生类描述了形状<i>functions</i>，即无限维度的对象：函数在每一点都有值。由于理论上的原因，我们需要这样做，因为我们想用函数的积分来进行分析。然而，对于计算机来说，这是一个非常困难的概念，因为它们一般只能处理有限的信息量，所以我们用正交点上的和来代替积分，我们通过使用定义在参考单元（正交对象）上的点映射（映射对象）到真实单元上的点来获得。实质上，我们将问题简化为我们只需要有限的信息，即形状函数值和导数、正交权重、法向量等，只需要在有限的点集合上。FEValues类就是将这三个部分结合在一起，并在一个特定的单元上提供这个有限的信息集  $K$  。当我们组装下面的线性系统时，你会看到它的作用。

值得注意的是，如果你只是在应用程序中自己创建这三个对象，并自己处理这些信息，那么所有这些也都可以实现。然而，这样做既不简单（FEValues类提供的正是你实际需要的信息），也不快：FEValues类经过高度优化，只在每个单元中计算你需要的特定信息；如果有任何东西可以从上一个单元中重复使用，那么它就会这样做，而且该类中有很多代码可以确保在有利的地方进行缓存。

这个介绍的最后一块是要提到，在得到一个线性系统后，要用迭代求解器进行求解，然后进行后处理：我们用DataOut类创建一个输出文件，然后可以用一个常见的可视化程序进行可视化。

 @note  前面对任何有限元实现的所有重要步骤的概述，在deal.II中也有对应的内容：该库可以自然地归纳为若干 "模块"，涵盖刚才概述的基本概念。你可以通过本页面顶部的标签访问这些模块。在<a href="index.html">front page of the deal.II manual</a>上也有对最基本的概念组的概述。




<a name="Abouttheimplementation"></a><h3>About the implementation</h3>


虽然这是你能用有限元方法解决的最简单的方程，但这个程序显示了大多数有限元程序的基本结构，也是几乎所有下面的程序基本上都会遵循的模板。具体来说，这个程序的主类看起来像这样。

@code
class Step3
{
  public:
    Step3 ();
    void run ();


  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;


    Triangulation<2>     triangulation;
    FE_Q<2>              fe;
    DoFHandler<2>        dof_handler;


    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    Vector<double>       solution;
    Vector<double>       system_rhs;
};
@endcode



这遵循了<a
href="http://en.wikipedia.org/wiki/Encapsulation_(object-oriented_programming)">data
encapsulation</a>的面向对象编程口号，也就是说，我们尽力将这个类的几乎所有内部细节隐藏在外部无法访问的私有成员中。

让我们从成员变量开始。这些遵循我们在上面的要点中所概述的构建模块，即我们需要一个三角形和一个DoFHandler对象，以及一个描述我们想要使用的各种形状函数的有限元对象。第二组对象与线性代数有关：系统矩阵和右手边以及解向量，还有一个描述矩阵稀疏模式的对象。这就是这个类所需要的全部内容（也是任何静止PDE的求解器所需要的基本内容），并且需要在整个程序中存活。与此相反，我们在装配时需要的FEValues对象只在整个装配过程中需要，因此我们在进行装配的函数中把它作为一个局部对象来创建，并在结束时再次销毁它。

其次，让我们来看看成员函数。这些，也已经构成了几乎所有下面的教程程序都会使用的共同结构。   <ul>   <li>   <code>make_grid()</code>  : 这就是人们所说的<i>preprocessing function</i>。顾名思义，它设置了存储三角图的对象。在以后的例子中，它还可以处理边界条件、几何形状等。     <li>   <code>setup_system()</code>  : 这是一个函数，其中设置了解决问题所需的所有其他数据结构。特别是，它将初始化DoFHandler对象并正确确定与线性代数有关的各种对象的大小。这个函数通常与上面的预处理函数分开，因为在一个与时间相关的程序中，每当网格被自适应细化时（我们将在步骤6中看到如何做），它可能至少每隔几个时间步就会被调用。另一方面，在上面的预处理函数中，设置网格本身只在程序开始时进行一次，因此，它被分离成自己的函数。     <li>   <code>assemble_system()</code>  : 这就是计算矩阵和右手边的内容的地方，在上面的介绍中已经详细讨论过。由于对这个线性系统进行处理在概念上与计算其条目有很大不同，我们将其与以下函数分开。     <li>   <code>solve()</code>  : 这就是我们计算线性系统 $U$ 的解的函数。在当前的程序中，这是一个简单的任务，因为矩阵是如此简单，但只要问题不再那么微不足道，它就会成为程序规模的重要部分（例如，一旦你对库有了更多的了解，请参阅步骤20，步骤22，或步骤31）。     <li>   <code>output_results()</code>  : 最后，当你计算出一个解决方案后，你可能想用它做一些事情。例如，你可能想以可视化的格式输出它，或者你可能想计算你感兴趣的量：例如，热交换器中的热通量、机翼的空气摩擦系数、最大桥梁载荷，或者仅仅是某一点上的数值解的值。因此，这个函数是对你的解进行后处理的地方。   </ul> 所有这些都是由单一的公共函数（除构造函数外），即 <code>run()</code> 函数来支撑的。它是在创建这种类型的对象的地方被调用的，它是按正确顺序调用所有其他函数的函数。把这个操作封装到 <code>run()</code> 函数中，而不是从 <code>main()</code> 中调用所有其他函数，确保你可以改变这个类中的关注点分离的实现方式。例如，如果其中一个函数变得太大了，你可以把它分成两个，而你唯一需要关注的地方就是这个类中的变化，而不是其他地方。

如上所述，你会看到这种一般的结构&mdash；有时在函数名称的拼写上会有一些变化，但基本上是按照这种功能分离的顺序&mdash；在下面的许多教程程序中也是如此。




<a name="Anoteontypes"></a><h3> A note on types </h3>


deal.II通过命名空间 dealii::types. 中的别名定义了一些积分%类型（在前一句中，"积分 "一词被用作与名词 "整数 "相对应的<i>adjective</i>。它不应该与表示曲线或曲面下的面积或体积的<i>noun</i>"积分 "混淆起来。形容词 "积分 "在C++世界中被广泛使用，如 "积分类型"、"积分常数 "等。）特别是，在这个程序中，你会在几个地方看到 types::global_dof_index ：一个整数类型，用来表示自由度的<i>global</i>索引，即在定义在三角形之上的DoFHandler对象中特定自由度的索引（而不是特定单元中的特定自由度的索引）。对于当前的程序（以及几乎所有的教程程序），你将有几千个到几百万个全局未知数（而且，对于 $Q_1$ 元素，你将有4个<i>locally on each cell</i>的2D和8个3D）。因此，允许为全局DoF指数存储足够大的数字的数据类型是 <code>unsigned int</code> ，因为它允许存储0到略高于40亿的数字（在大多数系统中，整数是32位的）。事实上，这就是 types::global_dof_index 的作用。

那么，为什么不马上使用 <code>unsigned int</code> 呢？deal.II在7.3版本之前一直是这样做的。然而，deal.II支持非常大的计算（通过步骤40中讨论的框架），当分布在几千个处理器上时，可能有超过40亿个未知数。因此，有些情况下 <code>unsigned int</code> 不够大，我们需要一个64位的无符号积分类型。为了实现这一点，我们引入了 types::global_dof_index ，它默认被定义为<code>unsigned int</code>，而如果有必要，可以通过在配置过程中传递一个特定的标志，将其定义为<code>unsigned long long int</code>（见ReadMe文件）。

这涵盖了技术方面。但是还有一个文档的目的：在图书馆和建立在它之上的代码中，如果你看到一个地方使用数据类型 types::global_dof_index, ，你就会立即知道被引用的数量实际上是一个全局dof指数。如果我们只是使用 <code>unsigned int</code> （它也可能是一个局部索引，一个边界指示器，一个材料ID，等等），就不会有这样的意义了。立即知道一个变量指的是什么也有助于避免错误：如果你看到一个 types::global_dof_index 类型的对象被分配给 types::subdomain_id, 类型的变量，这很明显一定有一个错误，尽管它们都是用无符号整数表示，因此编译器不会抱怨。

在更实际的情况下，这种类型的存在意味着在装配过程中，我们创建一个 $4\times 4$ 矩阵（在2d中，使用 $Q_1$ 元素）来表示我们当前所在单元的贡献，然后我们需要将这个矩阵的元素添加到全局（系统）矩阵的相应元素中。为此，我们需要获得当前单元的局部自由度的全局指数，为此我们将始终使用下面这段代码。

@code
  cell->get_dof_indices (local_dof_indices);
@endcode

其中 <code>local_dof_indices</code> 被声明为

@code
  std::vector<types::global_dof_index> local_dof_indices (fe.n_dofs_per_cell());
@endcode

这个变量的名字可能有点名不副实--它代表 "在当前单元上局部定义的那些自由度的全局指数"--但持有这种信息的变量在整个库中普遍是这样命名的。

 @note   types::global_dof_index  并不是这个命名空间中定义的唯一类型。相反，有一整个系列，包括 types::subdomain_id,  types::boundary_id, 和 types::material_id. 所有这些都是整数数据类型的别名，但正如上面所解释的，它们被用于整个库，以便（i）变量的意图变得更容易辨别，以及（ii）如果有必要，可以将实际类型改为一个更大的类型，而不必翻阅整个库，找出 <code>unsigned int</code> 的特定使用是否对应于，例如，一个材料指标。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Manynewincludefiles"></a> 
 * <h3>Many new include files</h3>
 * 

 * 
 * 这些包含文件已经为你所知。它们声明了处理三角形和自由度枚举的类。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * 
 * @endcode
 * 
 * 在这个文件中声明了创建网格的函数。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * @endcode
 * 
 * 这个文件包含了对拉格朗日插值有限元的描述。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_q.h> 
 * 
 * @endcode
 * 
 * 而这个文件是创建稀疏矩阵的稀疏模式所需要的，如前面的例子中所示。
 * 

 * 
 * 
 * @code
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * @endcode
 * 
 * 接下来的两个文件是在每个单元上使用正交法组装矩阵所需要的。下面将对其中声明的类进行解释。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * 
 * @endcode
 * 
 * 以下是我们在处理边界值时需要的三个包含文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * 
 * @endcode
 * 
 * 我们现在几乎到了终点。第二组到最后一组include文件是用于线性代数的，我们用它来解决拉普拉斯方程的有限元离散化所产生的方程组。我们将使用向量和全矩阵在每个单元中组装方程组，并将结果转移到稀疏矩阵中。然后我们将使用共轭梯度求解器来解决这个问题，为此我们需要一个预处理程序（在这个程序中，我们使用身份预处理程序，它没有任何作用，但我们还是需要包括这个文件）。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * 
 * @endcode
 * 
 * 最后，这是为了输出到文件和控制台。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/data_out.h> 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * ...这是为了将deal.II命名空间导入到全局范围。
 * 

 * 
 * 
 * @code
 * using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep3codeclass"></a> 
 * <h3>The <code>Step3</code> class</h3>
 * 

 * 
 * 在这个程序中，我们没有采用以前例子中的程序化编程，而是将所有东西都封装到一个类中。这个类由一些函数组成，这些函数分别执行有限元程序的某些方面，一个`main`函数控制先做什么和后做什么，还有一个成员变量列表。
 * 

 * 
 * 该类的公共部分相当简短：它有一个构造函数和一个从外部调用的函数`run`，其作用类似于`main`函数：它协调该类的哪些操作应以何种顺序运行。该类中的其他东西，即所有真正做事情的函数，都在该类的私有部分。
 * 

 * 
 * 
 * @code
 * class Step3 
 * { 
 * public: 
 *   Step3(); 
 * 
 *   void run(); 
 * 
 * @endcode
 * 
 * 然后，还有一些成员函数，它们主要是做它们名字所暗示的事情，在介绍中已经讨论过了。由于它们不需要从外部调用，所以它们是本类的私有函数。
 * 

 * 
 * 
 * @code
 * private: 
 *   void make_grid(); 
 *   void setup_system(); 
 *   void assemble_system(); 
 *   void solve(); 
 *   void output_results() const; 
 * 
 * @endcode
 * 
 * 最后我们还有一些成员变量。有一些变量描述了三角形和自由度的全局编号（我们将在这个类的构造函数中指定有限元的确切多项式程度）...
 * 

 * 
 * 
 * @code
 *   Triangulation<2> triangulation; 
 *   FE_Q<2>          fe; 
 *   DoFHandler<2>    dof_handler; 
 * 
 * @endcode
 * 
 * ...拉普拉斯方程离散化产生的系统矩阵的稀疏模式和数值的变量...
 * 

 * 
 * 
 * @code
 *   SparsityPattern      sparsity_pattern; 
 *   SparseMatrix<double> system_matrix; 
 * 
 * @endcode
 * 
 * .......以及用于保存右手边和解决方案向量的变量。
 * 

 * 
 * 
 * @code
 *   Vector<double> solution; 
 *   Vector<double> system_rhs; 
 * }; 
 * @endcode
 * 
 * 
 * <a name="Step3Step3"></a> 
 * <h4>Step3::Step3</h4>
 * 

 * 
 * 这里是构造函数。它除了首先指定我们需要双线性元素（由有限元对象的参数表示，它表示多项式的程度），并将dof_handler变量与我们使用的三角形相关联之外，没有做更多的工作。(注意，目前三角结构并没有设置网格，但是DoFHandler并不关心：它只想知道它将与哪个三角结构相关联，只有当你使用distribution_dofs()函数试图在网格上分布自由度时，它才开始关心实际的网格。) Step3类的所有其他成员变量都有一个默认的构造函数，它可以完成我们想要的一切。
 * 

 * 
 * 
 * @code
 * Step3::Step3() 
 *   : fe(1) 
 *   , dof_handler(triangulation) 
 * {} 
 * @endcode
 * 
 * 
 * <a name="Step3make_grid"></a> 
 * <h4>Step3::make_grid</h4>
 * 

 * 
 * 现在，我们要做的第一件事是生成我们想在其上进行计算的三角形，并对每个顶点进行自由度编号。我们之前在 step-1 和 step-2 中分别看到过这两个步骤。
 * 

 * 
 * 这个函数做的是第一部分，创建网格。 我们创建网格并对所有单元格进行五次细化。由于初始网格（也就是正方形 $[-1,1] \times [-1,1]$ ）只由一个单元组成，所以最终的网格有32乘以32个单元，总共是1024个。
 * 

 * 
 * 不确定1024是否是正确的数字？我们可以通过使用三角形上的 <code>n_active_cells()</code> 函数输出单元格的数量来检查。
 * 

 * 
 * 
 * @code
 * void Step3::make_grid() 
 * { 
 *   GridGenerator::hyper_cube(triangulation, -1, 1); 
 *   triangulation.refine_global(5); 
 * 
 *   std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *             << std::endl; 
 * } 
 * @endcode
 * 
 * @note  我们调用 Triangulation::n_active_cells() 函数，而不是 Triangulation::n_cells(). 这里，<i>active</i>指的是没有进一步提炼的单元。我们强调 "活跃 "这个形容词，因为还有更多的单元，即最细的单元的父单元，它们的父单元等等，直到构成初始网格的一个单元为止。当然，在下一个更粗的层次上，单元格的数量是最细层次上的单元格的四分之一，即256，然后是64、16、4和1。如果你在上面的代码中调用 <code>triangulation.n_cells()</code> ，你会因此得到一个1365的值。另一方面，单元格的数量（相对于活动单元格的数量）通常没有什么意义，所以没有很好的理由去打印它。
 * 

 * 
 * 
 * <a name="Step3setup_system"></a> 
 * <h4>Step3::setup_system</h4>
 * 

 * 
 * 接下来我们列举所有的自由度，并建立矩阵和向量对象来保存系统数据。枚举是通过使用 DoFHandler::distribute_dofs(), 来完成的，我们在 step-2 的例子中已经看到了。由于我们使用了FE_Q类，并且在构造函数中设置了多项式的度数为1，即双线性元素，这就将一个自由度与每个顶点联系起来。当我们在生成输出时，让我们也看看有多少自由度被生成。
 * 

 * 
 * 
 * @code
 * void Step3::setup_system() 
 * { 
 *   dof_handler.distribute_dofs(fe); 
 *   std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *             << std::endl; 
 * 
 * @endcode
 * 
 * 每个顶点应该有一个DoF。因为我们有一个32乘以32的网格，所以DoFs的数量应该是33乘以33，即1089。
 * 

 * 
 * 正如我们在前面的例子中所看到的，我们通过首先创建一个临时结构，标记那些可能为非零的条目，然后将数据复制到SparsityPattern对象中，然后可以被系统矩阵使用，来设置一个稀疏模式。
 * 

 * 
 * 
 * @code
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *   sparsity_pattern.copy_from(dsp); 
 * 
 * @endcode
 * 
 * 注意，SparsityPattern对象并不保存矩阵的值，它只保存条目所在的位置。条目本身存储在SparseMatrix类型的对象中，我们的变量system_matrix就是其中之一。
 * 

 * 
 * 稀疏模式和矩阵之间的区别是为了让几个矩阵使用相同的稀疏模式。这在这里似乎并不重要，但是当你考虑到矩阵的大小，以及建立稀疏模式可能需要一些时间时，如果你必须在程序中存储几个矩阵，这在大规模问题中就变得很重要了。
 * 

 * 
 * 
 * @code
 *   system_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 在这个函数中要做的最后一件事是将右侧向量和解向量的大小设置为正确的值。
 * 

 * 
 * 
 * @code
 *   solution.reinit(dof_handler.n_dofs()); 
 *   system_rhs.reinit(dof_handler.n_dofs()); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step3assemble_system"></a> 
 * <h4>Step3::assemble_system</h4>
 * 

 * 
 * 下一步是计算形成线性系统的矩阵和右手边的条目，我们从中计算出解决方案。这是每一个有限元程序的核心功能，我们在介绍中已经讨论了主要步骤。
 * 

 * 
 * 组装矩阵和向量的一般方法是在所有单元上循环，并在每个单元上通过正交计算该单元对全局矩阵和右侧的贡献。现在要认识到的一点是，我们需要实心单元上正交点位置的形状函数值。然而，有限元形状函数和正交点都只定义在参考单元上。因此，它们对我们帮助不大，事实上，我们几乎不会直接从这些对象中查询有关有限元形状函数或正交点的信息。
 * 

 * 
 * 相反，我们需要的是一种将这些数据从参考单元映射到实际单元的方法。能够做到这一点的类都是由Mapping类派生出来的，尽管人们常常不必直接与它们打交道：库中的许多函数都可以将映射对象作为参数，但当它被省略时，它们只是简单地诉诸于标准的双线性Q1映射。我们将走这条路，暂时不打扰它（我们将在 step-10 、 step-11 和 step-12 中再讨论这个问题）。
 * 

 * 
 * 所以我们现在有三个类的集合来处理：有限元、正交、和映射对象。这就太多了，所以有一种类型的类可以协调这三者之间的信息交流：FEValues类。如果给这三个对象各一个实例（或两个，以及一个隐式线性映射），它就能为你提供实心单元上正交点的形状函数值和梯度的信息。
 * 

 * 
 * 利用所有这些，我们将把这个问题的线性系统组装在以下函数中。
 * 

 * 
 * 
 * @code
 * void Step3::assemble_system() 
 * { 
 * 
 * @endcode
 * 
 * 好的，我们开始吧：我们需要一个正交公式来计算每个单元格的积分。让我们采用一个高斯公式，每个方向有两个正交点，即总共有四个点，因为我们是在二维。这个正交公式可以准确地积分三度以下的多项式（在一维）。很容易检查出，这对目前的问题来说是足够的。
 * 

 * 
 * 
 * @code
 *   QGauss<2> quadrature_formula(fe.degree + 1); 
 * 
 * @endcode
 * 
 * 然后我们初始化我们在上面简单谈及的对象。它需要被告知我们要使用哪个有限元，以及正交点和它们的权重（由一个正交对象共同描述）。如前所述，我们使用隐含的Q1映射，而不是自己明确指定一个。最后，我们必须告诉它我们希望它在每个单元上计算什么：我们需要正交点的形状函数值（对于右手 $(\varphi_i,f)$ ），它们的梯度（对于矩阵条目 $(\nabla \varphi_i, \nabla \varphi_j)$ ），以及正交点的权重和从参考单元到实际单元的雅各布变换的行列式。
 * 

 * 
 * 我们实际需要的信息列表是作为FEValues构造函数的第三个参数的标志集合给出的。由于这些值必须重新计算，或者说更新，每次我们进入一个新的单元时，所有这些标志都以前缀 <code>update_</code> 开始，然后指出我们想要更新的实际内容。如果我们想要计算形状函数的值，那么给出的标志是#update_values；对于梯度，它是#update_gradients。雅各布的行列式和正交权重总是一起使用的，所以只计算乘积（雅各布乘以权重，或者简称 <code>JxW</code> ）；由于我们需要它们，我们必须同时列出#update_JxW_values。
 * 

 * 
 * 
 * @code
 *   FEValues<2> fe_values(fe, 
 *                         quadrature_formula, 
 *                         update_values | update_gradients | update_JxW_values); 
 * 
 * @endcode
 * 
 * 这种方法的优点是，我们可以指定每个单元上究竟需要什么样的信息。很容易理解的是，这种方法可以大大加快有限元计算的速度，相比之下，所有的东西，包括二阶导数、单元的法向量等都在每个单元上计算，不管是否需要它们。
 * 

 * 
 * @note  <code>update_values | update_gradients | update_JxW_values</code>的语法对于那些不习惯用C语言编程多年的位操作的人来说不是很明显。首先， <code>operator|</code> 是<i>bitwise or operator</i>，也就是说，它接受两个整数参数，这些参数被解释为比特模式，并返回一个整数，其中每个比特都被设置，因为在两个参数中至少有一个的对应比特被设置。例如，考虑操作 <code>9|10</code>. In binary, <code>9=0b1001</code> （其中前缀 <code>0b</code> 表示该数字将被解释为二进制数字）和 <code>10=0b1010</code>  。通过每个比特，看它是否在其中一个参数中被设置，我们得出 <code>0b1001|0b1010=0b1011</code> ，或者用十进制符号表示， <code>9|10=11</code>  。你需要知道的第二个信息是，各种 <code>update_*</code> 标志都是有<i>exactly one bit set</i>的整数。例如，假设  <code>update_values=0b00001=1</code>  ,  <code>update_gradients=0b00010=2</code>  ,  <code>update_JxW_values=0b10000=16</code>  。那么<code>update_values | update_gradients | update_JxW_values = 0b10011 = 19</code>。换句话说，我们得到一个数字，即<i>encodes a binary mask representing all of the operations you want to happen</i>，其中每个操作正好对应于整数中的一个位，如果等于1，意味着每个单元格上应该更新一个特定的片断，如果是0，意味着我们不需要计算它。换句话说，即使 <code>operator|</code> 是<i>bitwise OR operation</i>，它真正代表的是<i>I want this AND that AND the other</i>。这样的二进制掩码在C语言编程中很常见，但在C++这样的高级语言中也许不是这样，但对当前的目的有很好的作用。
 * 

 * 
 * 为了在下文中进一步使用，我们为一个将被频繁使用的值定义了一个快捷方式。也就是每个单元的自由度数的缩写（因为我们是在二维，自由度只与顶点相关，所以这个数字是4，但是我们更希望在写这个变量的定义时，不妨碍我们以后选择不同的有限元，每个单元有不同的自由度数，或者在不同的空间维度工作）。
 * 

 * 
 * 一般来说，使用符号名称而不是硬编码这些数字是个好主意，即使你知道它们，因为例如，你可能想在某个时候改变有限元。改变元素就必须在不同的函数中进行，而且很容易忘记在程序的另一部分做相应的改变。最好不要依赖自己的计算，而是向正确的对象索取信息。在这里，我们要求有限元告诉我们每个单元的自由度数，无论我们在程序中的其他地方选择什么样的空间尺寸或多项式程度，我们都会得到正确的数字。
 * 

 * 
 * 这里定义的快捷方式主要是为了讨论基本概念，而不是因为它节省了大量的输入，然后会使下面的循环更容易阅读。在大型程序中，你会在很多地方看到这样的快捷方式，`dofs_per_cell`就是一个或多或少是这类对象的传统名称。
 * 

 * 
 * 
 * @code
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 * 
 * @endcode
 * 
 * 现在，我们说我们想逐个单元地组装全局矩阵和向量。我们可以将结果直接写入全局矩阵，但是这样做的效率并不高，因为对稀疏矩阵元素的访问是很慢的。相反，我们首先在一个小矩阵中计算每个单元的贡献，并在这个单元的计算结束后将其转移到全局矩阵中。我们对右手边的向量也是这样做的。所以我们首先分配这些对象（这些是局部对象，所有的自由度都与所有其他的自由度耦合，我们应该使用一个完整的矩阵对象，而不是一个用于局部操作的稀疏矩阵；以后所有的东西都将转移到全局的稀疏矩阵中）。
 * 

 * 
 * 
 * @code
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *   Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 在集合每个单元的贡献时，我们用自由度的局部编号（即从零到dofs_per_cell-1的编号）来做。然而，当我们将结果转移到全局矩阵时，我们必须知道自由度的全局编号。当我们查询它们时，我们需要为这些数字建立一个从头开始的（临时）数组（关于这里使用的类型， types::global_dof_index, ，见介绍末尾的讨论）。
 * 

 * 
 * 
 * @code
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 现在是所有单元格的循环。我们之前已经看到这对一个三角形是如何工作的。DoFHandler的单元格迭代器与Triangulation的迭代器完全类似，但有关于你所使用的有限元的自由度的额外信息。在自由度处理程序的活动单元上进行循环操作的方法与三角法相同。
 * 

 * 
 * 注意，这次我们将单元的类型声明为`const auto &`，而不是`auto`。在第1步中，我们通过用细化指标标记来修改三角形的单元。在这里，我们只检查单元格而不修改它们，所以把`cell`声明为`const`是很好的做法，以便执行这个不变性。
 * 

 * 
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators()) 
 *     { 
 * 
 * @endcode
 * 
 * 我们现在坐在一个单元上，我们希望计算形状函数的值和梯度，以及参考单元和真实单元之间映射的雅各布矩阵的行列式，在正交点上。由于所有这些值都取决于单元格的几何形状，我们必须让FEValues对象在每个单元格上重新计算它们。
 * 

 * 
 * 
 * @code
 *       fe_values.reinit(cell); 
 * 
 * @endcode
 * 
 * 接下来，在我们填充之前，将本地单元对全局矩阵和全局右手边的贡献重置为零。
 * 

 * 
 * 
 * @code
 *       cell_matrix = 0; 
 *       cell_rhs    = 0; 
 * 
 * @endcode
 * 
 * 现在是时候开始对单元进行积分了，我们通过对所有的正交点进行循环来完成，我们将用q_index来编号。
 * 

 * 
 * 
 * @code
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
 *         { 
 * 
 * @endcode
 * 
 * 首先组装矩阵。对于拉普拉斯问题，每个单元格上的矩阵是形状函数i和j的梯度的积分。由于我们不进行积分，而是使用正交，所以这是在所有正交点的积分之和乘以正交点的雅各布矩阵的行列式乘以这个正交点的权重。你可以通过使用 <code>fe_values.shape_grad(i,q_index)</code> 得到形状函数 $i$ 在数字q_index的正交点上的梯度；这个梯度是一个二维向量（事实上它是张量 @<1,dim@>, 类型，这里dim=2），两个这样的向量的乘积是标量乘积，即两个shape_grad函数调用的积是点乘。这又要乘以雅各布行列式和正交点权重（通过调用 FEValues::JxW() 得到）。最后，对所有形状函数 $i$ 和 $j$ 重复上述操作。
 * 

 * 
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices()) 
 *             for (const unsigned int j : fe_values.dof_indices()) 
 *               cell_matrix(i, j) += 
 *                 (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
 *                  fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
 *                  fe_values.JxW(q_index));           // dx 
 * 
 * @endcode
 * 
 * 然后我们对右手边做同样的事情。在这里，积分是对形状函数i乘以右手边的函数，我们选择的是常值为1的函数（更有趣的例子将在下面的程序中考虑）。
 * 

 * 
 * 
 * @code
 *           for (const unsigned int i : fe_values.dof_indices()) 
 *             cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
 *                             1. *                                // f(x_q) 
 *                             fe_values.JxW(q_index));            // dx 
 *         } 
 * 
 * @endcode
 * 
 * 现在我们有了这个单元的贡献，我们必须把它转移到全局矩阵和右手边。为此，我们首先要找出这个单元上的自由度有哪些全局数字。让我们简单地询问该单元的信息。
 * 

 * 
 * 
 * @code
 *       cell->get_dof_indices(local_dof_indices); 
 * 
 * @endcode
 * 
 * 然后再次循环所有形状函数i和j，并将局部元素转移到全局矩阵中。全局数字可以用local_dof_indices[i]获得。
 * 

 * 
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices()) 
 *         for (const unsigned int j : fe_values.dof_indices()) 
 *           system_matrix.add(local_dof_indices[i], 
 *                             local_dof_indices[j], 
 *                             cell_matrix(i, j)); 
 * 
 * @endcode
 * 
 * 再来，我们对右边的向量做同样的事情。
 * 

 * 
 * 
 * @code
 *       for (const unsigned int i : fe_values.dof_indices()) 
 *         system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *     } 
 * 
 * @endcode
 * 
 * 现在，几乎所有的东西都为离散系统的求解做好了准备。然而，我们还没有照顾到边界值（事实上，没有迪里切特边界值的拉普拉斯方程甚至不是唯一可解的，因为你可以在离散解中加入一个任意的常数）。因此，我们必须对这种情况做一些处理。
 * 

 * 
 * 为此，我们首先获得边界上的自由度列表以及形状函数在那里的值。为了简单起见，我们只对边界值函数进行插值，而不是将其投影到边界上。库中有一个函数正是这样做的。  VectorTools::interpolate_boundary_values(). 它的参数是（省略存在默认值而我们不关心的参数）：DoFHandler对象，用于获取边界上自由度的全局数字；边界上边界值应被内插的部分；边界值函数本身；以及输出对象。
 * 

 * 
 * 边界分量的含义如下：在很多情况下，你可能只想在边界的一部分施加某些边界值。例如，在流体力学中，你可能有流入和流出的边界，或者在身体变形计算中，身体的夹紧和自由部分。那么你就想用指标来表示边界的这些不同部分，并告诉interpolate_boundary_values函数只计算边界的某一部分（例如夹住的部分，或流入的边界）的边界值。默认情况下，所有的边界都有一个0的边界指标，除非另有规定。如果边界的部分有不同的边界条件，你必须用不同的边界指示器为这些部分编号。然后，下面的函数调用将只确定那些边界指标实际上是作为第二个参数指定的0的边界部分的边界值。
 * 

 * 
 * 描述边界值的函数是一个Function类型的对象或一个派生类的对象。其中一个派生类是 Functions::ZeroFunction, ，它描述了一个到处都是零的函数（并不意外）。我们就地创建这样一个对象，并将其传递给 VectorTools::interpolate_boundary_values() 函数。
 * 

 * 
 * 最后，输出对象是一对全局自由度数（即边界上的自由度数）和它们的边界值（这里所有条目都是零）的列表。这种自由度数到边界值的映射是由 <code>std::map</code> 类完成的。
 * 

 * 
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values; 
 *   VectorTools::interpolate_boundary_values(dof_handler, 
 *                                            0, 
 *                                            Functions::ZeroFunction<2>(), 
 *                                            boundary_values); 
 * 
 * @endcode
 * 
 * 现在我们得到了边界DoF的列表和它们各自的边界值，让我们用它们来相应地修改方程组。这可以通过以下函数调用来实现。
 * 

 * 
 * 
 * @code
 *   MatrixTools::apply_boundary_values(boundary_values, 
 *                                      system_matrix, 
 *                                      solution, 
 *                                      system_rhs); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step3solve"></a> 
 * <h4>Step3::solve</h4>
 * 

 * 
 * 下面的函数简单地求解了离散化的方程。由于该系统对于高斯消除或LU分解等直接求解器来说是一个相当大的系统，我们使用共轭梯度算法。你应该记住，这里的变量数量（只有1089个）对于有限元计算来说是一个非常小的数字，而100.000是一个比较常见的数字。 对于这个数量的变量，直接方法已经不能使用了，你不得不使用CG这样的方法。
 * 

 * 
 * 
 * @code
 * void Step3::solve() 
 * { 
 * 
 * @endcode
 * 
 * 首先，我们需要有一个对象，知道如何告诉CG算法何时停止。这是通过使用SolverControl对象来实现的，作为停止标准，我们说：在最多1000次迭代后停止（这远远超过了1089个变量的需要；见结果部分以了解真正使用了多少次），如果残差的规范低于 $10^{-12}$ 就停止。在实践中，后一个标准将是停止迭代的一个标准。
 * 

 * 
 * 
 * @code
 *   SolverControl solver_control(1000, 1e-12); 
 * 
 * @endcode
 * 
 * 然后，我们需要解算器本身。SolverCG类的模板参数是向量的类型，留下空的角括号将表明我们采取的是默认参数（即 <code>Vector@<double@></code>  ）。然而，我们明确地提到了模板参数。
 * 

 * 
 * 
 * @code
 *   SolverCG<Vector<double>> solver(solver_control); 
 * 
 * @endcode
 * 
 * 现在求解方程组。CG求解器的第四个参数是一个预处理程序。我们觉得还没有准备好深入研究这个问题，所以我们告诉它使用身份运算作为预处理。
 * 

 * 
 * 
 * @code
 *   solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
 * 
 * @endcode
 * 
 * 现在求解器已经完成了它的工作，求解变量包含了求解函数的结点值。
 * 

 * 
 * 
 * @code
 * } 
 * @endcode
 * 
 * 
 * <a name="Step3output_results"></a> 
 * <h4>Step3::output_results</h4>
 * 

 * 
 * 典型的有限元程序的最后一部分是输出结果，也许会做一些后处理（例如计算边界处的最大应力值，或者计算整个流出物的平均通量，等等）。我们这里没有这样的后处理，但是我们想把解决方案写到一个文件里。
 * 

 * 
 * 
 * @code
 * void Step3::output_results() const 
 * { 
 * 
 * @endcode
 * 
 * 为了将输出写入文件，我们需要一个知道输出格式等的对象。这就是DataOut类，我们需要一个该类型的对象。
 * 

 * 
 * 
 * @code
 *   DataOut<2> data_out; 
 * 
 * @endcode
 * 
 * 现在我们必须告诉它从哪里获取它要写的值。我们告诉它使用哪个DoFHandler对象，以及求解向量（以及求解变量在输出文件中的名称）。如果我们有不止一个我们想在输出中查看的向量（例如右手边，每个单元格的错误，等等），我们也要把它们加进去。
 * 

 * 
 * 
 * @code
 *   data_out.attach_dof_handler(dof_handler); 
 *   data_out.add_data_vector(solution, "solution"); 
 * 
 * @endcode
 * 
 * 在DataOut对象知道它要处理哪些数据后，我们必须告诉它把它们处理成后端可以处理的数据。原因是我们将前端（知道如何处理DoFHandler对象和数据向量）与后端（知道许多不同的输出格式）分开，使用一种中间数据格式将数据从前端传输到后端。数据通过以下函数转换为这种中间格式。
 * 

 * 
 * 
 * @code
 *   data_out.build_patches(); 
 * 
 * @endcode
 * 
 * 现在我们已经为实际输出做好了一切准备。只要打开一个文件，用VTK格式把数据写进去就可以了（在我们这里使用的DataOut类中还有很多其他函数，可以把数据写成postscript、AVS、GMV、Gnuplot或其他一些文件格式）。
 * 

 * 
 * 
 * @code
 *   std::ofstream output("solution.vtk"); 
 *   data_out.write_vtk(output); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step3run"></a> 
 * <h4>Step3::run</h4>
 * 

 * 
 * 最后，这个类的最后一个函数是主函数，调用 <code>Step3</code> 类的所有其他函数。这样做的顺序类似于大多数有限元程序的工作顺序。由于这些名字大多是不言自明的，所以没有什么可评论的。
 * 

 * 
 * 
 * @code
 * void Step3::run() 
 * { 
 *   make_grid(); 
 *   setup_system(); 
 *   assemble_system(); 
 *   solve(); 
 *   output_results(); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 这是程序的主函数。由于主函数的概念大多是C++编程之前的面向对象时代的遗留物，所以它通常不做更多的事情，只是创建一个顶层类的对象并调用其原理函数。
 * 

 * 
 * 最后，函数的第一行是用来启用deal.II可以生成的一些诊断程序的输出。  @p deallog 变量（代表deal-log，而不是de-allog）代表一个流，库的某些部分将输出写入其中。例如，迭代求解器将产生诊断程序（起始残差、求解器步骤数、最终残差），在运行这个教程程序时可以看到。
 * 

 * 
 * @p deallog 的输出可以写到控制台，也可以写到文件，或者两者都写。两者在默认情况下都是禁用的，因为多年来我们已经知道，一个程序只应该在用户明确要求的时候才产生输出。但这是可以改变的，为了解释如何做到这一点，我们需要解释 @p deallog 是如何工作的。当库的个别部分想要记录输出时，它们会打开一个 "上下文 "或 "部分"，这个输出将被放入其中。在想要写输出的部分结束时，人们再次退出这个部分。由于一个函数可以在这个输出部分打开的范围内调用另一个函数，所以输出实际上可以分层嵌套到这些部分。LogStream类（ @p deallog 是一个变量）将这些部分中的每一个称为 "前缀"，因为所有的输出都以这个前缀打印在行的左端，前缀由冒号分隔。总是有一个默认的前缀叫做 "DEAL"（暗示了deal.II的历史，它是以前一个叫做 "DEAL "的库的继承者，LogStream类是被带入deal.II的少数代码之一）。
 * 

 * 
 * 默认情况下， @p logstream 只输出前缀为零的行--也就是说，所有的输出都是禁用的，因为默认的 "DEAL "前缀总是存在的。但人们可以为应该输出的行设置不同的最大前缀数，以达到更大的效果，事实上在这里我们通过调用 LogStream::depth_console(). 将其设置为两个。这意味着对于所有的屏幕输出，在默认的 "DEAL "之外再推一个前缀的上下文被允许将其输出打印到屏幕上（"控制台"），而所有进一步嵌套的部分将有三个或更多的前缀被激活，会写到 @p deallog, ，但 @p deallog 并不转发这个输出到屏幕。因此，运行这个例子（或者看 "结果 "部分），你会看到解算器的统计数据前缀为 "DEAL:CG"，这是两个前缀。这对于当前程序的上下文来说已经足够了，但是你将在以后看到一些例子（例如，在 step-22 中），其中求解器嵌套得更深，你可能通过设置更高的深度来获得有用的信息。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   deallog.depth_console(2); 
 * 
 *   Step3 laplace_problem; 
 *   laplace_problem.run(); 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-3/doc/results.dox



<a name="Results"></a><h1>Results</h1>


程序的输出看起来如下。

@code
Number of active cells: 1024
Number of degrees of freedom: 1089
DEAL:cg::Starting value 0.121094
DEAL:cg::Convergence step 48 value 5.33692e-13
@endcode



前两行是我们写给  <code>cout</code>  的内容。最后两行是CG求解器在没有我们的干预下生成的。前两行说明了迭代开始时的残差，而最后一行告诉我们求解器需要47次迭代才能使残差的规范值达到5.3e-13，即低于我们在 "solve "函数中设置的阈值1e-12。我们将在下一个程序中展示如何抑制这种输出，这种输出有时对调试很有用，但往往会使屏幕显示变得混乱。

除了上面显示的输出，该程序还生成了文件 <code>solution.vtk</code> ，该文件为VTK格式，被当今许多可视化程序广泛使用--包括两个重量级的<a href="https://www.llnl.gov/visit">VisIt</a>和<a href="https://www.paraview.org">Paraview</a>，是当今最常使用的程序。

使用VisIt，生成一张像这样的解决方案的图片并不是很困难。   <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-3.solution-3.png" alt="Visualization of the solution of step-3">
    </td>
  </tr>
</table>  它同时显示了解和网格，根据每一点的解的值提升到 $x$  -  $y$ 平面之上。当然，这里的解并不特别令人兴奋，但这是拉普拉斯方程所代表的内容和我们为这个程序选择的右手边 $f(\mathbf x)=1$ 的结果。拉普拉斯方程描述了（在许多其他用途中）受外部（也是垂直）力作用的膜的垂直变形。在目前的例子中，膜的边界被夹在一个没有垂直变化的方形框架上；因此，一个恒定的力密度将直观地导致膜简单地向上隆起--就像上图所示。

VisIt和Paraview都允许玩各种可视化的解决方案。几个视频讲座展示了如何使用这些程序。   @dealiiVideoLectureSeeAlso{11,32} 




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


如果你想用这个程序玩一玩，这里有几个建议。   </p> 

 <ul>   <li>  改变几何图形和网格。在程序中，我们通过使用 <code>GridGenerator::hyper_cube</code> 函数生成了一个方形域和网格。然而， <code>GridGenerator</code> 也有大量的其他函数。试试L形域，环形域，或其他你在那里找到的域。     </li> 

    <li>  改变边界条件。代码使用 Functions::ZeroFunction 函数来生成零边界条件。然而，你可能想用 <code>ConstantFunction&lt;2&gt;(1)</code> 而不是 <code>ZeroFunction&lt;2&gt;()</code> 尝试非零常数边界值，以获得单位Dirichlet边界值。在函数命名空间的文档中描述了更多的奇异函数，你可以选择一个来描述你的特定边界值。     </li> 

    <li>  修改边界条件的类型。目前，发生的情况是，我们在周围使用迪里希特边界值，因为默认情况是所有边界部分的边界指标为零，然后我们告诉 VectorTools::interpolate_boundary_values() 函数，在所有指标为零的边界部分上将边界值插值为零。    <p>  如果我们给边界的部分分配不同的指标，我们可以改变这种行为。例如，在调用 GridGenerator::hyper_cube(): @code
  triangulation.begin_active()->face(0)->set_boundary_id(1);
  @endcode后立即尝试这样做。



  这样做的目的是，首先要求三角剖分返回一个迭代器，指向第一个活动单元。当然，由于这是一个正方形的三角测量的粗略网格，此刻三角测量只有一个单元，而且它是活动的。接下来，我们要求单元格返回它的第一个面的迭代器，然后我们要求面将该面的边界指标重置为1。接下来的事情就是这样。当网格被细化时，子单元的面会继承其父母的边界指示器，也就是说，即使在最细的网格上，广场一侧的面的边界指示器为1。稍后，当我们要插值边界条件时， VectorTools::interpolate_boundary_values() 调用将只为那些边界指标为零的面产生边界值，而对那些具有不同边界指标的面则不予理会。这样做的目的是对前者施加Dirichlet边界条件，而对后者施加同质的Neumann条件（即解的法向导数为零，除非在变分等式的右侧添加额外的条款来处理潜在的非零Neumann条件）。如果你运行该程序，你会看到这一点。

  另一种改变边界指标的方法是根据面中心的笛卡尔坐标来标注边界。   例如，我们可以通过检查单元格中心的y坐标是否在-1和1的公差（这里是1e-12）范围内，将沿上下边界的所有单元格标记为边界指示器1。在调用 GridGenerator::hyper_cube(), 后，像以前一样立即尝试这样做。   @code
  for (auto &face : triangulation.active_face_iterators())
    if (std::fabs(face->center()(1) - (-1.0)) < 1e-12 ||
        std::fabs(face->center()(1) - (1.0)) < 1e-12)
      face->set_boundary_id(1);
  @endcode

  虽然这段代码比以前长了一些，但它对复杂的几何形状很有用，因为它不需要脸部标签的知识。

    <li> 最后一点的一个小变化是像上面那样设置不同的边界值，但随后为边界指标一使用不同的边界值函数。在实践中，你要做的是为边界指标一增加对 <code>interpolate_boundary_values</code> 的第二次调用。   @code
  VectorTools::interpolate_boundary_values(dof_handler,
					   1,
					   ConstantFunction<2>(1.),
					   boundary_values);
  @endcode

  如果你在这个函数的第一个调用之后立即进行这个调用，那么它将把边界指标为1的面的边界值内插到单位值，并将这些内插值与之前计算的边界指标为0的值合并。

    <li>  观察收敛情况。我们将只讨论第7步中规范的计算误差，但很容易检查计算在这里已经收敛了。例如，我们可以在一个点上评估解的值，并比较不同%的全局细化的值（全局细化的步骤数在上面的 <code>LaplaceProblem::make_grid</code> 中设定）。为了评估某个点的解决方案，例如在 $(\frac 13, \frac 13)$ ，我们可以在 <code>LaplaceProblem::output_results</code> 函数中加入以下代码。   @code
    std::cout << "Solution at (1/3,1/3): "
              << VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3))
              << std::endl;
  @endcode

  对于1到9个全局细化步骤，我们就会得到以下的点值序列。     <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$u_h(\frac 13,\frac13)$</th> </tr>
    <tr> <td>1</td> <td>0.166667</td> </tr>
    <tr> <td>2</td> <td>0.227381</td> </tr>
    <tr> <td>3</td> <td>0.237375</td> </tr>
    <tr> <td>4</td> <td>0.240435</td> </tr>
    <tr> <td>5</td> <td>0.241140</td> </tr>
    <tr> <td>6</td> <td>0.241324</td> </tr>
    <tr> <td>7</td> <td>0.241369</td> </tr>
    <tr> <td>8</td> <td>0.241380</td> </tr>
    <tr> <td>9</td> <td>0.241383</td> </tr>
  </table>  通过注意到每两个连续值之间的差异减少了大约4倍，我们可以猜测 "正确 "的值可能是 $u(\frac 13, \frac 13)\approx 0.241384$  。事实上，如果我们假设这是正确的值，我们可以证明上面的序列确实显示了 ${\cal
  O}(h^2)$ 的收敛&mdash；理论上，收敛顺序应该是 ${\cal O}(h^2 |\log h|)$ ，但是领域和网格的对称性可能导致了观察到的更好的收敛顺序。

  这方面的一个小变种是用二次元重复测试。你需要做的就是在构造函数中把有限元的多项式程度设置为2  <code>LaplaceProblem::LaplaceProblem</code>  。

    <li>  平均值的收敛。一个不同的方法是计算解的平均数，以了解解是否真的收敛了（收敛到什么程度&mdash；我们无法判断它是否真的是正确的值！）。为此，在 <code>LaplaceProblem::output_results</code> 中添加以下代码：@code
    std::cout << "Mean value: "
              << VectorTools::compute_mean_value (dof_handler,
						  QGauss<2>(fe.degree + 1),
						  solution,
						  0)
              << std::endl;
  @endcode

  该函数的文档解释了第二和第四个参数的含义，而第一和第三个参数应该是很明显的。再次做同样的研究，我们改变了全局细化步骤的数量，我们得到以下结果。     <table align="center" class="doxtable">
    <tr> <th># of refinements</th> <th>$\int_\Omega u_h(x)\; dx$</th> </tr>
    <tr> <td>0</td> <td>0.09375000</td> </tr>
    <tr> <td>1</td> <td>0.12790179</td> </tr>
    <tr> <td>2</td> <td>0.13733440</td> </tr>
    <tr> <td>3</td> <td>0.13976069</td> </tr>
    <tr> <td>4</td> <td>0.14037251</td> </tr>
    <tr> <td>5</td> <td>0.14052586</td> </tr>
    <tr> <td>6</td> <td>0.14056422</td> </tr>
    <tr> <td>7</td> <td>0.14057382</td> </tr>
    <tr> <td>8</td> <td>0.14057622</td> </tr>
  </table>  同样，两个相邻值之间的差异下降了约四倍，表明收敛为  ${\cal O}(h^2)$  。   </ul> 




<a name="UsingHDF5tooutputthesolutionandadditionaldata"></a><h3>Using %HDF5 to output the solution and additional data</h3>


%HDF5是一种常用的格式，可以被许多脚本语言（如R或Python）读取。让deal.II产生一些%HDF5文件并不困难，然后可以在外部脚本中使用，对该程序产生的一些数据进行后处理。这里有一些关于可能的想法。




<a name="Changingtheoutputtoh5"></a><h4> Changing the output to .h5</h4>


为了充分利用自动化，我们首先需要为全局细化步骤的数量引入一个私有变量 <code>unsigned int n_refinement_steps </code> ，它将被用于输出文件名。在 <code>make_grid()</code> we then replace <code>triangulation.refine_global(5);</code> 中用

@code
n_refinement_steps = 5;
triangulation.refine_global(n_refinement_steps);
@endcode

deal.II库有两个不同的%HDF5绑定，一个在HDF5命名空间（用于对接通用数据文件），另一个在DataOut（专门用于为解决方案的可视化写文件）。尽管HDF5 deal.II绑定支持串行和MPI，但%HDF5 DataOut绑定只支持并行输出。由于这个原因，我们需要初始化一个只有一个处理器的MPI通信器。这可以通过添加以下代码来实现。

@code
int main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  ...
}
@endcode

接下来我们改变 `Step3::output_results()` 的输出例程，如DataOutBase命名空间文档中所述。

@code
const std::string filename_h5 = "solution_" + std::to_string(n_refinement_steps) + ".h5";
DataOutBase::DataOutFilterFlags flags(true, true);
DataOutBase::DataOutFilter data_filter(flags);
data_out.write_filtered_data(data_filter);
data_out.write_hdf5_parallel(data_filter, filename_h5, MPI_COMM_WORLD);
@endcode

然后，产生的文件可以被可视化，就像教程的原始版本产生的VTK文件一样；但是，由于%HDF5是一种更通用的文件格式，它也可以很容易地用脚本语言处理，用于其他目的。




<a name="Addingthepointvalueandthemeanseeextensionaboveintotheh5file"></a><h4> Adding the point value and the mean (see extension above) into the .h5 file</h4>


在输出解决方案后，可以再次打开该文件以包括更多的数据集。  这使得我们可以将实验的所有必要信息保存在一个结果文件中，然后可以由一些后处理脚本来读取和处理。关于可能的输出选项，请看 HDF5::Group::write_dataset() 的进一步信息）。

为了实现这一点，我们首先将必要的头文件纳入我们的文件。

@code
#include <deal.II/base/hdf5.h>
@endcode

在我们的输出例程的末尾添加以下几行，将关于某一点的解的值，以及解的平均值的信息添加到我们的%HDF5文件中。

@code
HDF5::File data_file(filename_h5, HDF5::File::FileAccessMode::open, MPI_COMM_WORLD);
Vector<double> point_value(1);
point_value[0] = VectorTools::point_value(dof_handler, solution,
                                          Point<2>(1./3, 1./3));
data_file.write_dataset("point_value", point_value);
Vector<double> mean_value(1);
mean_value[0] = VectorTools::compute_mean_value(dof_handler,
                                                QGauss<2>(fe.degree + 1),
                                                solution, 0);
data_file.write_dataset("mean_value",mean_value);
@endcode






<a name="UsingRandggplot2togenerateplots"></a><h3> Using R and ggplot2 to generate plots</h3>


上述放入%HDF5文件的数据，然后可以从脚本语言中使用，进行进一步的后处理。在下文中，让我们展示一下，特别是如何用<a href="https://en.wikipedia.org/wiki/R_(programming_language)">R
programming language</a>这个在统计数据分析中广泛使用的语言来完成。(例如，类似的事情也可以在Python中完成。)如果你不熟悉R和ggplot2，你可以看看R的数据木工课程<a href="https://datacarpentry.org/R-ecology-lesson/index.html">here</a>。此外，由于大多数搜索引擎对 "R+主题 "这种形式的搜索很吃力，我们建议使用专门的服务<a
href="http://rseek.org">RSeek </a>来代替。

R和其他语言最突出的区别是，赋值运算符（`a = 5`）通常被写成`a <- 5`。由于后者被认为是标准的，我们将在我们的例子中也使用它。要在R语言中打开`.h5`文件，你必须安装<a href="https://bioconductor.org/packages/release/bioc/html/rhdf5.html">rhdf5</a>包，它是Bioconductor软件包的一部分。

首先，我们将包括所有必要的包，并看看我们文件中的数据是如何结构化的。

@code{.r}
library(rhdf5)     # library for handling HDF5 files
library(ggplot2)   # main plotting library
library(grDevices) # needed for output to PDF
library(viridis)   # contains good colormaps for sequential data


refinement <- 5
h5f <- H5Fopen(paste("solution_",refinement,".h5",sep=""))
print(h5f)
@endcode

这给出了以下输出

@code{.unparsed}
HDF5 FILE
   name /
filename


    name       otype  dclass     dim
0 cells       H5I_DATASET INTEGER  x 1024
1 mean_value  H5I_DATASET FLOAT   1
2 nodes       H5I_DATASET FLOAT    x 1089
3 point_value H5I_DATASET FLOAT   1
4 solution    H5I_DATASET FLOAT    x 1089
@endcode

数据集可以通过  <code>h5f\$name</code>  访问。函数  <code>dim(h5f\$cells)</code>  给我们提供了用于存储我们单元格的矩阵的尺寸。我们可以看到以下三个矩阵，以及我们添加的两个额外数据点。   <ul>   <li>   <code>cells</code>  ：一个4x1024的矩阵，存储每个单元的（C++）顶点指数  <li>   <code>nodes</code>  ：一个2x1089的矩阵，存储我们单元顶点的位置值（x，y）  <li>   <code>solution</code>  : 一个1x1089的矩阵，存储我们的解决方案在每个顶点的值  </ul>  现在我们可以使用这些数据来生成各种图表。用ggplot2作图通常分为两步。首先，数据需要被处理并添加到一个  <code>data.frame</code>  。之后，构建一个 <code>ggplot</code> 对象，并通过向其添加绘图元素来进行操作。

 <code>nodes</code> and <code>cells</code> 包含我们绘制网格所需的所有信息。下面的代码将所有的数据打包成一个数据框架，用于绘制我们的网格。

@code{.r}
# Counting in R starts at 1 instead of 0, so we need to increment all
# vertex indices by one:
cell_ids <- h5f$cells+1


# Store the x and y positions of each vertex in one big vector in a
# cell by cell fashion (every 4 entries belong to one cell):
cells_x <- h5f$nodes[1,][cell_ids]
cells_y <- h5f$nodes[2,][cell_ids]


# Construct a vector that stores the matching cell by cell grouping
# (1,1,1,1,2,2,2,2,...):
groups <- rep(1:ncol(cell_ids),each=4)


# Finally put everything into one dataframe:
meshdata <- data.frame(x = cells_x, y = cells_y, id = groups)
@endcode



有了完成的数据框架，我们就有了绘制网格所需的一切。

@code{.r}
pdf (paste("grid_",refinement,".pdf",sep=""),width = 5,height = 5) # Open new PDF file
plt <- ggplot(meshdata,aes(x=x,y=y,group=id))                      # Construction of our plot
                                                                   # object, at first only data


plt <- plt + geom_polygon(fill="white",colour="black")             # Actual plotting of the grid as polygons
plt <- plt + ggtitle(paste("grid at refinement level #",refinement))


print(plt)                                                         # Show the current state of the plot/add it to the pdf
dev.off()                                                          # Close PDF file
@endcode



这个文件的内容看起来如下（不是很令人兴奋，但你会明白的）。   <table width="60%" align="center">
  <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.grid_5.png" alt="Grid after 5 refinement steps of step-3">
   </td>
  </tr>
</table> 

我们还可以将解决方案本身可视化，这看起来会更有趣。为了给我们的解决方案做一个二维伪色图，我们将使用  <code>geom_raster</code>  。这个函数需要一个结构化的网格，即在x和y方向上是均匀的。幸运的是，我们在这一点上的数据是以正确的方式结构化的。下面的代码将我们的曲面的伪彩色表示法绘制成一个新的PDF。

@code{.r}
pdf (paste("pseudocolor_",refinement,".pdf",sep=""),width = 5,height = 4.2) # Open new PDF file
colordata <- data.frame(x = h5f$nodes[1,],y = h5f$nodes[2,] , solution = h5f$solution[1,])
plt <- ggplot(colordata,aes(x=x,y=y,fill=solution))
plt <- plt + geom_raster(interpolate=TRUE)
plt <- plt + scale_fill_viridis()
plt <- plt + ggtitle(paste("solution at refinement level #",refinement))


print(plt)
dev.off()
H5Fclose(h5f) # Close the HDF5 file
@endcode

现在的情况是这样的。   <table width="60%" align="center">
 <tr>
   <td align="center">
     <img src="https://www.dealii.org/images/steps/developer/step-3.extensions.pseudocolor_5.png" alt="Solution after 5 refinement steps of step-3">
   </td>
 </tr>
</table> 

为了绘制收敛曲线，我们需要从1开始用不同的 <code>n_refinement_steps</code> 值多次重新运行C++代码。由于每个文件只包含一个数据点，我们需要对它们进行循环，并将结果串联成一个矢量。

@code{.r}
n_ref <- 8   # Maximum refinement level for which results are existing


# First we initiate all vectors with the results of the first level
h5f   <- H5Fopen("solution_1.h5")
dofs  <- dim(h5f$solution)[2]
mean  <- h5f$mean_value
point <- h5f$point_value
H5Fclose(h5f)


for (reflevel in 2:n_ref)
{
   h5f   <- H5Fopen(paste("solution_",reflevel,".h5",sep=""))
   dofs  <- c(dofs,dim(h5f\$solution)[2])
   mean  <- c(mean,h5f\$mean_value)
   point <- c(point,h5f\$point_value)
   H5Fclose(h5f)
}
@endcode

由于我们对数值本身不感兴趣，而是对与 "精确 "解决方案相比的误差感兴趣，我们将假设我们的最高细化水平是该解决方案，并从数据中省略它。

@code{.r}
# Calculate the error w.r.t. our maximum refinement step
mean_error  <- abs(mean[1:n_ref-1]-mean[n_ref])
point_error <- abs(point[1:n_ref-1]-point[n_ref])


# Remove the highest value from our DoF data
dofs     <- dofs[1:n_ref-1]
convdata <- data.frame(dofs = dofs, mean_value= mean_error, point_value = point_error)
@endcode

现在我们有所有的数据可以用来生成我们的图。在对数尺度上绘制误差往往是有用的，这在下面的代码中可以实现。

@code
pdf (paste("convergence.pdf",sep=""),width = 5,height = 4.2)
plt <- ggplot(convdata,mapping=aes(x = dofs, y = mean_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "mean value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)


plt <- ggplot(convdata,mapping=aes(x = dofs, y = point_value))
plt <- plt+geom_line()
plt <- plt+labs(x="#DoFs",y = "point value error")
plt <- plt+scale_x_log10()+scale_y_log10()
print(plt)


dev.off()
@endcode

这就产生了下面的图，显示了均值和所选点的解值的误差如何很好地收敛到零。   <table style="width:50%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_mean.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-3.extensions.convergence_point.png" alt=""></td>
  </tr>
</table> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-3.cc"
*/
