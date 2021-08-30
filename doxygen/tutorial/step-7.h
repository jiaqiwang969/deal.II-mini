/**
@page step_7 The step-7 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Verificationofcorrectness">Verification of correctness</a>
        <li><a href="#NonhomogeneousNeumannboundaryconditions">Non-homogeneous Neumann boundary conditions</a>
        <li><a href="#Anoteongoodprogrammingpractice">A note on good programming practice</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#TheHelmholtzsolverclass">The Helmholtz solver class</a>
        <li><a href="#TheHelmholtzProblemclassimplementation">The HelmholtzProblem class implementation</a>
      <ul>
        <li><a href="#HelmholtzProblemHelmholtzProblemconstructor">HelmholtzProblem::HelmholtzProblem constructor</a>
        <li><a href="#HelmholtzProblemsetup_system">HelmholtzProblem::setup_system</a>
        <li><a href="#HelmholtzProblemassemble_system">HelmholtzProblem::assemble_system</a>
        <li><a href="#HelmholtzProblemsolve">HelmholtzProblem::solve</a>
        <li><a href="#HelmholtzProblemrefine_grid">HelmholtzProblem::refine_grid</a>
        <li><a href="#HelmholtzProblemprocess_solution">HelmholtzProblem::process_solution</a>
        <li><a href="#HelmholtzProblemrun">HelmholtzProblem::run</a>
      <ul>
        <li><a href="#Outputofgraphicaldata">Output of graphical data</a>
        <li><a href="#Outputofconvergencetables">Output of convergence tables</a>
        <li><a href="#Furthertablemanipulations">Further table manipulations</a>
      </ul>
      </ul>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
      <ul>
        <li><a href="#Whenistheerrorsmall"> When is the error "small"? </a>
      </ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#HigherOrderElements"> Higher Order Elements </a>
        <li><a href="#ConvergenceComparison"> Convergence Comparison </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-7/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


在这个程序中，我们将主要考虑两个方面。<ol>  <li>  验证程序的正确性和生成收敛表；  <li>  亥姆霍兹方程的非均质诺伊曼边界条件。   </ol>  除了这些主题，还将再次展示各种改进和技巧。




<a name="Verificationofcorrectness"></a><h3>Verification of correctness</h3>


可能从来没有一个非琐碎的有限元程序从一开始就能正常工作。因此有必要找到方法来验证计算出的解是否正确。通常，这是通过选择模拟的设置来完成的，我们知道准确的连续解，并评估连续解和计算的离散解之间的差异。如果这个差值以正确的收敛顺序收敛为零，这已经是正确性的一个很好的指示，尽管可能还有其他的误差源持续存在，而这些误差对总误差的贡献很小，或者是更高阶的。在有限元模拟的背景下，这种通过选择适当的右手边和边界条件来选择解决方案的技术通常被称为<i>Method of Manufactured Solution</i>。

在这个例子中，我们不会去研究系统性软件验证的理论，这是一个非常复杂的问题。相反，我们将展示deal.II在这方面所能提供的工具。这基本上是围绕着一个单一的函数的功能， VectorTools::integrate_difference(). 这个函数计算一个给定的连续函数和每个单元上不同规范的有限元场之间的差值。当然，像其他积分一样，我们只能用正交公式来评估这些规范；因此，选择正确的正交公式对准确评估误差至关重要。这一点对于 $L_\infty$ 规范尤其成立，我们只在正交点评估数值和精确解的最大偏差；那么就不应该尝试使用正交规则，其评估只发生在[超级收敛](https://en.wikipedia.org/wiki/Superconvergence)可能发生的点上，例如最低阶高斯正交公式的高斯点，对于该矩阵的装配中的积分是正确的（例如，对于线性元素，不要使用QGauss（2）正交公式）。事实上，这通常也是对其他规范的好建议：如果你的正交点偶然选择在由于超收敛而恰好误差特别小的位置，那么计算出来的误差看起来会比实际情况小得多，甚至可能表明收敛阶数更高。因此，我们将为这些误差准则的积分选择不同的正交公式，而不是线性系统的装配。

函数 VectorTools::integrate_difference() 对三角形的每个单元 $K$ 进行所需规范的评估，并返回一个持有每个单元的这些值的向量。从局部值中，我们可以得到全局误差。例如，如果所有单元 $K$ 的元素为 $e_K$ 的向量 $\mathbf e$ 包含本地 $L_2$ 规范 $\|u-u_h\|_K$ ，那么

@f[
  E = \| {\mathbf e} \| = \left( \sum_K e_K^2 \right)^{1/2}


@f]

是全球  $L_2$  误差  $E=\|u-u_h\|_\Omega$  。

在程序中，我们将展示如何评估和使用这些量，并且我们将监测它们在网格细化下的值。当然，我们必须选择手头的问题，使我们能够明确地说明解决方案及其导数，但由于我们要评估程序的正确性，这才是合理的。如果我们知道程序对一个（或者，如果想真正确定：许多）特别选择的右手边产生了正确的解决方案，我们可以相当有信心，它也将计算出我们不知道确切数值的问题的正确解决方案。

除了简单地计算这些数量之外，我们将展示如何从这个程序产生的数据中生成格式良好的表格，自动计算收敛率等。此外，我们还将比较不同的网格细化策略。




<a name="NonhomogeneousNeumannboundaryconditions"></a><h3>Non-homogeneous Neumann boundary conditions</h3>


本例程序的第二个完全不相关的主题是使用非均质边界条件。这些条件包括在使用边界积分的变异形式中，我们必须在组装右手边的矢量时进行数值评估。

在我们进入编程之前，让我们简单看看数学公式。我们在这里要解决的方程是亥姆霍兹方程，"带有漂亮的符号"。

@f[


  -\Delta u + \alpha u = f,


@f]

在 $[-1,1]^2$ 与 $\alpha=1$ 的正方形上，用迪里希特边界条件来增强。

@f[
  u = g_1


@f]

在边界 $\Gamma_1$ 的某些部分 $\Gamma$ ，以及诺伊曼条件

@f[
  {\mathbf n}\cdot \nabla u = g_2


@f]

对其余的 $\Gamma_2 = \Gamma \backslash \Gamma_1$ 。在我们特定的测试案例中，我们将使用 $\Gamma_1=\Gamma \cap\{\{x=1\}
\cup \{y=1\}\}$  。 我们说这个方程有 "漂亮的符号"，因为算子 $-\Delta + \alpha I$ 与身份 $I$ 和 $\alpha>0$ 是一个正定算子；<a
href="https://en.wikipedia.org/wiki/Helmholtz_equation">equation with
the "bad sign"</a>是 $-\Delta u - \alpha u$ ，由时间谐波过程建模产生。如果 $\alpha>0$ 很大的话，该算子就不是正定的，这就导致了我们在这里不需要讨论的各种问题。如果 $\alpha$ 恰好是 $-\Delta$ 的特征值之一，那么算子也可能不是可逆的--即方程没有唯一的解--)

因为我们想验证我们的数字解  $u_h$  的收敛性，我们想要一个设置，以便我们知道精确的解  $u$  。这就是 "人造解法 "的作用。为此，让我们选择一个函数

@f[
  \bar u(x) = \sum_{i=1}^3 \exp\left(-\frac{|x-x_i|^2}{\sigma^2}\right)


@f]

其中指数的中心 $x_i$ 为 $x_1=(-\frac 12,\frac 12)$ ， $x_2=(-\frac 12,-\frac 12)$ ，和 $x_3=(\frac 12,-\frac 12)$ ，半宽设为 $\sigma=\frac {1}{8}$ 。然后，制造解的方法说：选择

@f{align*}
  f &= -\Delta \bar u + \bar u, \\
  g_1 &= \bar u|_{\Gamma_1}, \\
  g_2 &= {\mathbf n}\cdot \nabla\bar u|_{\Gamma_2}.


@f}

有了这个特殊的选择，我们就可以推断出原始问题的解恰好是 $u=\bar u$  。换句话说，通过以特定的方式选择方程的右边和边界条件，我们已经为自己制造了一个我们知道其解决方案的问题。这使我们能够计算出我们的数字解决方案的误差。在下面的代码中，我们用 <code>Solution</code> 类来表示 $\bar u$ ，其他类将用于表示 $\bar u|_{\Gamma_1}=g_1$ 和 ${\mathbf n}\cdot \nabla\bar u|_{\Gamma_2}=g_2$  。

利用上述定义，我们可以陈述方程的弱表述，即：找到 $u\in H^1_g=\{v\in H^1: v|_{\Gamma_1}=g_1\}$ ，以便

@f[
  {(\nabla v, \nabla u)}_\Omega + {(v,u)}_\Omega
  =
  {(v,f)}_\Omega + {(v,g_2)}_{\Gamma_2}


@f]

为所有测试函数  $v\in H^1_0=\{v\in H^1: v|_{\Gamma_1}=0\}$  。边界项 ${(v,g_2)}_{\Gamma_2}$ 是通过部分积分和使用 $\partial_n u=g_2$ 对 $\Gamma_2$ 和 $v=0$ 的方式出现的。因此，在离散公式中，我们用来建立全局矩阵和右侧向量的单元矩阵和向量看起来是这样的。

@f{eqnarray*}
  A_{ij}^K &=& \left(\nabla \varphi_i, \nabla \varphi_j\right)_K
              +\left(\varphi_i, \varphi_j\right)_K,
  \\
  F_i^K &=& \left(\varphi_i, f\right)_K
           +\left(\varphi_i, g_2\right)_{\partial K\cap \Gamma_2}.


@f}

由于域积分的生成已经在前面的例子中多次展示过，这里只对轮廓积分的生成感兴趣。它的工作原理如下：对于域积分，我们有一个 <code>FEValues</code> 类，它提供了形状值和梯度，以及雅各布行列式和其他信息，并指定了单元格中的正交点；同样，还有一个 <code>FEFaceValues</code> 类，为单元格面上的积分执行这些任务。我们向它提供一个流形的正交公式，该流形的维数比域的维数少一，以及我们想在其上进行积分的单元格和其面的数量。然后，该类将计算该面的正交点的值、梯度、法向量、权重等，然后我们可以以与域积分相同的方式使用这些值。下面的程序显示了这样做的细节。




<a name="Anoteongoodprogrammingpractice"></a><h3>A note on good programming practice</h3>


除了上面概述的数学主题外，我们还想用这个程序来说明良好编程实践的一个方面，即命名空间的使用。在为deal.II库编程时，我们非常注意不要为类和全局函数使用过于通用的名称，例如 <code>f(), sz(), rhs()</code> 等。此外，我们把所有的东西都放到了名字空间  <code>dealii</code>  中。但当一个人写的应用程序不是为了给别人使用时，就不一定会注意这么多。如果你遵循步骤1到步骤6的编程风格，这些函数就会最终进入全局命名空间，不幸的是，很多其他的东西也在那里（基本上是C语言提供的所有东西，以及你通过头文件从操作系统得到的所有东西）。更糟糕的是，C语言的设计者在避免使用通用名称方面也不总是很小心；例如，符号<code>j1, jn</code>在C头文件中被定义（它们表示贝塞尔函数）。

为了避免不同函数或变量的名称相撞时产生的问题（经常出现混乱的错误信息），把你所做的一切放到<a
href="http://en.wikipedia.org/wiki/Namespace_(computer_science)">namespace</a>中是一个好的做法。按照这种风格，我们将在程序的顶部打开一个名字空间 <code>Step7</code> ，将deal.II名字空间导入其中，将本程序特有的一切（除了 <code>main()</code> ，它必须在全局名字空间中）放入其中，并且只在文件的底部关闭它。换句话说，这个程序的结构是这样的

@code
  #includes ...


  namespace Step7
  {
    using namespace dealii;


    ...everything to do with the program...
  }


  int main ()
  {
    ...do whatever main() does...
  }
@endcode

我们将在整个交易.II教程的剩余部分遵循这一方案。


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
 * 这些第一个包含文件在前面的例子中都已经处理过了，所以我们不再解释其中的内容。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * @endcode
 * 
 * 在这个例子中，我们将不使用DoFHandler类默认使用的编号方案，而是使用Cuthill-McKee算法对其进行重新编号。正如在 step-2 中已经解释过的，必要的函数被声明在以下文件中。
 * 

 * 
 * 
 * @code
 * #include <deal.II/dofs/dof_renumbering.h> 
 * 
 * @endcode
 * 
 * 然后我们将展示一个小技巧，如何确保对象在仍在使用时不被删除。为此，deal.II有一个SmartPointer辅助类，它被声明在这个文件中。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/smartpointer.h> 
 * 
 * @endcode
 * 
 * 接下来，我们要使用介绍中提到的函数 VectorTools::integrate_difference() ，我们要使用一个ConvergenceTable，在运行过程中收集所有重要的数据，并在最后以表格形式打印出来。这些来自于以下两个文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/base/convergence_table.h> 
 * 
 * @endcode
 * 
 * 最后，我们需要使用FEFaceValues类，它与FEValues类在同一个文件中声明。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <array> 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 在我们继续实际执行之前的最后一步是打开一个命名空间 <code>Step7</code> ，我们将把所有的东西放进去，正如在介绍的最后所讨论的，并把命名空间 <code>dealii</code> 的成员导入其中。
 * 

 * 
 * 
 * @code
 * namespace Step7 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 在实现实际求解的类之前，我们首先声明和定义一些代表右手边和求解类的函数类。由于我们要将数值得到的解与精确的连续解进行比较，我们需要一个代表连续解的函数对象。另一方面，我们需要右手边的函数，而这个函数当然与解共享一些特征。为了减少如果我们必须同时改变两个类中的某些东西而产生的依赖性，我们将两个函数的共同特征移到一个基类中。
 * 

 * 
 * 解（正如介绍中所解释的，我们选择三个指数之和）和右手边的共同特征是：指数的数量，它们的中心，以及它们的半宽。我们在以下类别中声明它们。由于指数的数量是一个编译时的常数，我们使用一个固定长度的 <code>std::array</code> 来存储中心点。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class SolutionBase 
 *   { 
 *   protected: 
 *     static const std::array<Point<dim>, 3> source_centers; 
 *     static const double                    width; 
 *   }; 
 * 
 * @endcode
 * 
 * 表示指数中心和宽度的变量刚刚被声明，现在我们还需要给它们赋值。在这里，我们可以展示另一个小小的模板魔法，即我们如何根据维度给这些变量分配不同的值。我们将在程序中只使用2维的情况，但我们展示1维的情况是为了说明一个有用的技术。
 * 

 * 
 * 首先我们为1d情况下的中心赋值，我们将中心等距离地放在-1/3、0和1/3处。这个定义的<code>template &lt;&gt;</code>头显示了一个明确的专业化。这意味着，这个变量属于一个模板，但是我们并没有向编译器提供一个模板，让它通过用一些具体的值来替代 <code>dim</code> 来专门化一个具体的变量，而是自己提供一个专门化，在这个例子中是 <code>dim=1</code>  。如果编译器在模板参数等于1的地方看到了对这个变量的引用，它就知道它不需要通过替换 <code>dim</code> 从模板中生成这个变量，而是可以立即使用下面的定义。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   const std::array<Point<1>, 3> SolutionBase<1>::source_centers = { 
 *     {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}}; 
 * 
 * @endcode
 * 
 * 同样地，我们可以为 <code>dim=2</code> 提供一个明确的特殊化。我们将2d情况下的中心放置如下。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   const std::array<Point<2>, 3> SolutionBase<2>::source_centers = { 
 *     {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}}; 
 * 
 * @endcode
 * 
 * 还需要给指数的半宽指定一个值。我们希望对所有维度使用相同的数值。在这种情况下，我们只需向编译器提供一个模板，它可以通过用一个具体的值替换 <code>dim</code> 来生成一个具体的实例。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   const double SolutionBase<dim>::width = 1. / 8.; 
 * 
 * @endcode
 * 
 * 在声明和定义了解和右手的特征后，我们可以声明代表这两者的类。它们都代表连续函数，所以它们都派生于Function&lt;dim&gt;基类，它们也继承了SolutionBase类中定义的特征。
 * 

 * 
 * 实际的类是在下面声明的。请注意，为了计算数值解与连续解在L2和H1（半）准则下的误差，我们必须提供精确解的值和梯度。这比我们在以前的例子中所做的要多，在以前的例子中，我们所提供的只是一个或一列点的值。幸运的是，Function类也有用于梯度的虚拟函数，所以我们可以简单地重载Function基类中各自的虚拟成员函数。请注意，一个函数在 <code>dim</code> 空间维度上的梯度是一个大小为 <code>dim</code> 的向量，即一个等级为1、维度为 <code>dim</code> 的张量。就像其他很多东西一样，该库提供了一个合适的类。这个类的一个新特点是，它明确地使用了张量对象，之前在  step-3  和  step-4  中作为中间词出现。张量是标量（等级为零的张量）、向量（等级为一的张量）和矩阵（等级为二的张量）以及高维对象的概括。张量类需要两个模板参数：张量等级和张量维度。例如，在这里我们使用等级为一的张量（向量），维度为 <code>dim</code> (so they have <code>dim</code> 项）。虽然这比使用Vector的灵活性要差一些，但当编译时知道向量的长度时，编译器可以生成更快的代码。此外，指定一个秩为1、维数为 <code>dim</code> 的张量，可以保证张量具有正确的形状（因为它是内置于对象本身的类型中的），所以编译器可以为我们抓住大多数与尺寸有关的错误。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solution : public Function<dim>, protected SolutionBase<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *     virtual Tensor<1, dim> 
 *     gradient(const Point<dim> & p, 
 *              const unsigned int component = 0) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 精确解类的值和梯度的实际定义是根据其数学定义，不需要过多解释。
 * 

 * 
 * 唯一值得一提的是，如果我们访问一个依赖模板的基类的元素（在本例中是SolutionBase&lt;dim&gt;的元素），那么C++语言会强迫我们写  <code>this-&gt;source_centers</code>  ，对于基类的其他成员也是如此。如果基类不依赖模板，C++就不需要 <code>this-&gt;</code> 的限定。这一点的原因很复杂，C++书籍会在<i>two-stage (name) lookup</i>这句话下进行解释，在deal.II FAQs中也有很长的描述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double Solution<dim>::value(const Point<dim> &p, const unsigned int) const 
 *   { 
 *     double return_value = 0; 
 *     for (const auto &center : this->source_centers) 
 *       { 
 *         const Tensor<1, dim> x_minus_xi = p - center; 
 *         return_value += 
 *           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)); 
 *       } 
 * 
 *     return return_value; 
 *   } 
 * 
 * @endcode
 * 
 * 同样，这也是对解的梯度的计算。 为了从指数的贡献中积累梯度，我们分配了一个对象  <code>return_value</code>  ，它表示秩  <code>1</code>  和维  <code>dim</code>  的张量的数学量。它的默认构造函数将其设置为只包含零的向量，所以我们不需要明确关心它的初始化。
 * 

 * 
 * 注意，我们也可以把对象的类型定为Point&lt;dim&gt;，而不是Tensor&lt;1,dim&gt;。等级1的张量和点几乎是可以交换的，而且只有非常细微的数学含义不同。事实上，Point&lt;dim&gt;类是由Tensor&lt;1,dim&gt;类派生出来的，这就弥补了它们的相互交换能力。它们的主要区别在于它们在逻辑上的含义：点是空间中的点，比如我们要评估一个函数的位置（例如，见这个函数的第一个参数的类型）。另一方面，秩1的张量具有相同的变换属性，例如，当我们改变坐标系时，它们需要以某种方式旋转；然而，它们不具有点所具有的相同内涵，只是比坐标方向所跨越的空间更抽象的对象。事实上，梯度生活在 "对等 "的空间中，因为它们的分量的维度不是长度，而是长度上的一个）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   Tensor<1, dim> Solution<dim>::gradient(const Point<dim> &p, 
 *                                          const unsigned int) const 
 *   { 
 *     Tensor<1, dim> return_value; 
 * 
 *     for (const auto &center : this->source_centers) 
 *       { 
 *         const Tensor<1, dim> x_minus_xi = p - center; 
 * 
 * @endcode
 * 
 * 对于梯度，注意它的方向是沿着（x-x_i），所以我们把这个距离向量的倍数加起来，其中的因子是由指数给出。
 * 

 * 
 * 
 * @code
 *         return_value += 
 *           (-2. / (this->width * this->width) * 
 *            std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) * 
 *            x_minus_xi); 
 *       } 
 * 
 *     return return_value; 
 *   } 
 * 
 * @endcode
 * 
 * 除了代表精确解的函数外，我们还需要一个函数，在组装离散方程的线性系统时，我们可以将其作为右手。这可以通过下面的类和其函数的定义来实现。请注意，这里我们只需要函数的值，而不是它的梯度或高阶导数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim>, protected SolutionBase<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 右手边的值是由解的负拉普拉斯加上解本身给出的，因为我们要解决亥姆霍兹方程的问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> &p, 
 *                                    const unsigned int) const 
 *   { 
 *     double return_value = 0; 
 *     for (const auto &center : this->source_centers) 
 *       { 
 *         const Tensor<1, dim> x_minus_xi = p - center; 
 * 
 * @endcode
 * 
 * 第一个贡献是拉普拉斯的。
 * 

 * 
 * 
 * @code
 *         return_value += 
 *           ((2. * dim - 
 *             4. * x_minus_xi.norm_square() / (this->width * this->width)) / 
 *            (this->width * this->width) * 
 *            std::exp(-x_minus_xi.norm_square() / (this->width * this->width))); 
 * 
 * @endcode
 * 
 * 而第二个是解决方案本身。
 * 

 * 
 * 
 * @code
 *         return_value += 
 *           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)); 
 *       } 
 * 
 *     return return_value; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TheHelmholtzsolverclass"></a> 
 * <h3>The Helmholtz solver class</h3>
 * 

 * 
 * 然后我们需要做所有工作的类。除了它的名字，它的接口与前面的例子基本相同。
 * 

 * 
 * 其中一个不同点是，我们将在几种模式下使用这个类：用于不同的有限元，以及用于自适应细化和全局细化。全局细化还是自适应细化的决定是通过在类的顶部声明的枚举类型传达给该类的构造函数的。构造函数接收一个有限元对象和细化模式作为参数。
 * 

 * 
 * 除了 <code>process_solution</code> 函数外，其余的成员函数与之前一样。在解被计算出来后，我们对它进行一些分析，比如计算各种规范的误差。为了实现一些输出，它需要细化周期的编号，因此得到它作为一个参数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class HelmholtzProblem 
 *   { 
 *   public: 
 *     enum RefinementMode 
 *     { 
 *       global_refinement, 
 *       adaptive_refinement 
 *     }; 
 * 
 *     HelmholtzProblem(const FiniteElement<dim> &fe, 
 *                      const RefinementMode      refinement_mode); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void process_solution(const unsigned int cycle); 
 * 
 * @endcode
 * 
 * 现在是这个类的数据元素。在我们以前的例子中已经使用过的变量中，只有有限元对象不同。这个类的对象所操作的有限元被传递给这个类的构造函数。它必须存储一个指向有限元的指针，供成员函数使用。现在，对于本类来说，这没有什么大不了的，但由于我们想在这些程序中展示技术而不是解决方案，我们将在这里指出一个经常出现的问题--当然也包括正确的解决方案。
 * 

 * 
 * 考虑以下在所有示例程序中出现的情况：我们有一个三角形对象，我们有一个有限元对象，我们还有一个DoFHandler类型的对象，它同时使用前两个对象。这三个对象的寿命与其他大多数对象相比都相当长：它们基本上是在程序开始时或外循环时设置的，并在最后被销毁。问题是：我们能否保证DoFHandler使用的两个对象的寿命至少与它们被使用的时间相同？这意味着DoFHandler必须对其他对象的销毁情况有一定的了解。
 * 

 * 
 * 我们将在这里展示库如何设法找出对一个对象仍有活动的引用，并且从使用对象的角度来看，该对象仍然活着。基本上，该方法是沿着以下思路进行的：所有受到这种潜在危险的指针的对象都来自一个叫做Subscriptor的类。例如，Triangulation、DoFHandler和FiniteElement类的一个基类都派生于Subscriptor。后面这个类并没有提供太多的功能，但是它有一个内置的计数器，我们可以订阅这个计数器，因此这个类的名字就叫 "订阅器"。每当我们初始化一个指向该对象的指针时，我们可以增加它的使用计数器，而当我们移开指针或不再需要它时，我们再减少计数器。这样，我们就可以随时检查有多少个对象还在使用该对象。此外，该类需要知道一个指针，它可以用来告诉订阅对象它的无效性。
 * 

 * 
 * 如果一个从Subscriptor类派生出来的对象被销毁，它也必须调用Subscriptor类的析构函数。在这个析构器中，我们使用存储的指针告诉所有订阅的对象该对象的无效性。当对象出现在移动表达式的右侧时，也会发生同样的情况，也就是说，在操作后它将不再包含有效的内容。在试图访问被订阅的对象之前，订阅类应该检查存储在其相应指针中的值。
 * 

 * 
 * 这正是SmartPointer类正在做的事情。它基本上就像一个指针一样，也就是说，它可以被取消引用，可以被分配给其他指针，等等。除此之外，当我们试图解除引用这个类所代表的指针时，它使用上面描述的机制来找出这个指针是否是悬空的。在这种情况下，会抛出一个异常。
 * 

 * 
 * 在本例程序中，我们希望保护有限元对象，避免因某种原因导致所指向的有限元在使用中被破坏。因此，我们使用了一个指向有限元对象的SmartPointer；由于有限元对象在我们的计算中实际上从未改变，我们传递了一个const FiniteElement&lt;dim&gt;作为SmartPointer类的模板参数。请注意，这样声明的指针是在构造求解对象时被分配的，并在销毁时被销毁，所以对有限元对象销毁的锁定贯穿了这个HelmholtzProblem对象的生命周期。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     SmartPointer<const FiniteElement<dim>> fe; 
 * 
 *     AffineConstraints<double> hanging_node_constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 * 
 * @endcode
 * 
 * 倒数第二个变量存储了传递给构造函数的细化模式。由于它只在构造函数中设置，我们可以声明这个变量为常数，以避免有人不由自主地设置它（例如在一个 "if "语句中，==偶然被写成=）。
 * 

 * 
 * 
 * @code
 *     const RefinementMode refinement_mode; 
 * 
 * @endcode
 * 
 * 对于每个细化级别，一些数据（比如单元格的数量，或者数值解的L2误差）将被生成，并在之后打印出来。TableHandler可以用来收集所有这些数据，并在运行结束后以简单文本或LaTeX格式的表格输出。这里我们不仅使用TableHandler，还使用了派生类ConvergenceTable，它还可以评估收敛率。
 * 

 * 
 * 
 * @code
 *     ConvergenceTable convergence_table; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="TheHelmholtzProblemclassimplementation"></a> 
 * <h3>The HelmholtzProblem class implementation</h3>
 * 
 * <a name="HelmholtzProblemHelmholtzProblemconstructor"></a> 
 * <h4>HelmholtzProblem::HelmholtzProblem constructor</h4>
 * 

 * 
 * 在这个类的构造函数中，我们只设置作为参数传递的变量，并将DoF处理程序对象与三角形（不过目前是空的）相关联。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   HelmholtzProblem<dim>::HelmholtzProblem(const FiniteElement<dim> &fe, 
 *                                           const RefinementMode refinement_mode) 
 *     : dof_handler(triangulation) 
 *     , fe(&fe) 
 *     , refinement_mode(refinement_mode) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemsetup_system"></a> 
 * <h4>HelmholtzProblem::setup_system</h4>
 * 

 * 
 * 下面的函数设置了自由度、矩阵和向量的大小等。它的大部分功能在前面的例子中已经展示过了，唯一不同的是在第一次分配自由度后立即进行重新编号的步骤。
 * 

 * 
 * 重编自由度并不难，只要你使用库中的一种算法。它只需要一行代码。这方面的更多信息可以在  step-2  中找到。
 * 

 * 
 * 但是请注意，当你对自由度进行重新编号时，你必须在分配自由度后立即进行，因为诸如悬空节点、稀疏模式等都取决于重新编号后的绝对数。
 * 

 * 
 * 我们在这里介绍重新编号的原因是，这是一个相对便宜的操作，但往往有一个有利的效果。虽然CG迭代本身与自由度的实际排序无关，但我们将使用SSOR作为一个预处理程序。SSOR会经过所有的自由度，并做一些取决于之前发生的操作；因此，SSOR操作并不独立于自由度的编号，而且众所周知，它的性能会通过使用重新编号技术得到改善。一个小实验表明，确实如此，例如，用这里使用的Q1程序进行自适应细化的第五个细化周期的CG迭代次数，在没有重编号的情况下为40次，而在重编号的情况下为36次。对于这个程序中的所有计算，一般都可以观察到类似的节省。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(*fe); 
 *     DoFRenumbering::Cuthill_McKee(dof_handler); 
 * 
 *     hanging_node_constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     hanging_node_constraints.condense(dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemassemble_system"></a> 
 * <h4>HelmholtzProblem::assemble_system</h4>
 * 

 * 
 * 为手头的问题组装方程组，主要是像之前的例子程序一样。然而，无论如何，有些东西已经改变了，所以我们对这个函数进行了相当广泛的评论。
 * 

 * 
 * 在该函数的顶部，你会发现通常的各种变量声明。与以前的程序相比，重要的是我们希望解决的问题也是双二次元的，因此必须使用足够精确的正交公式。此外，我们需要计算面的积分，即 <code>dim-1</code> 维的对象。那么，面的正交公式的声明就很直接了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::assemble_system() 
 *   { 
 *     QGauss<dim>     quadrature_formula(fe->degree + 1); 
 *     QGauss<dim - 1> face_quadrature_formula(fe->degree + 1); 
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     const unsigned int dofs_per_cell = fe->n_dofs_per_cell(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 然后我们需要一些对象来评估正交点上的形状函数的值、梯度等。虽然看起来用一个对象来做域积分和面积分应该是可行的，但是有一个微妙的区别，因为域积分的权重包括域中单元的度量，而面积分的正交需要低维流形中面的度量。在内部，这两个类都根植于一个共同的基类，它完成了大部分工作，并为域积分和面积分提供了相同的接口。
 * 

 * 
 * 对于亥姆霍兹方程的双线性形式的域积分，我们需要计算值和梯度，以及正交点的权重。此外，我们需要实细胞上的正交点（而不是单位细胞上的正交点）来评估右手边的函数。我们用来获取这些信息的对象是之前讨论过的FEValues类。
 * 

 * 
 * 对于面积分，我们只需要形状函数的值以及权重。我们还需要实心单元上的法向量和正交点，因为我们要从精确解对象中确定Neumann值（见下文）。给我们提供这些信息的类被称为FEFaceValues。
 * 

 * 
 * 
 * @code
 *     FEValues<dim> fe_values(*fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_face_values(*fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_quadrature_points | 
 *                                        update_normal_vectors | 
 *                                        update_JxW_values); 
 * 
 * @endcode
 * 
 * 然后我们需要一些从以前的例子中已经知道的对象。一个表示右侧函数的对象，它在单元格上正交点的值，单元格矩阵和右侧，以及单元格上自由度的指数。
 * 

 * 
 * 请注意，我们对右手边对象的操作只是查询数据，绝不会改变该对象。因此我们可以声明它  <code>const</code>  。
 * 

 * 
 * 
 * @code
 *     const RightHandSide<dim> right_hand_side; 
 *     std::vector<double>      rhs_values(n_q_points); 
 * 
 * @endcode
 * 
 * 最后我们定义一个表示精确解函数的对象。我们将用它来计算边界上的诺伊曼值。通常情况下，我们当然会使用一个单独的对象来计算，特别是由于精确解通常是未知的，而诺伊曼值是规定的。然而，我们将有点偷懒，使用我们已经有的信息。当然，现实生活中的程序会在这里采取其他方式。
 * 

 * 
 * 
 * @code
 *     Solution<dim> exact_solution; 
 * 
 * @endcode
 * 
 * 现在是所有单元格的主循环。这与之前的例子基本没有变化，所以我们只对有变化的地方进行评论。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix = 0.; 
 *         cell_rhs    = 0.; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(), 
 *                                    rhs_values); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 * 
 * @endcode
 * 
 * 第一件改变的事情是双线性形式。它现在包含了亥姆霍兹方程的附加项。
 * 

 * 
 * 
 * @code
 *                 cell_matrix(i, j) += 
 *                   ((fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q) 
 *                       fe_values.shape_grad(j, q_point)     // grad phi_j(x_q) 
 *                     +                                      // 
 *                     fe_values.shape_value(i, q_point) *    // phi_i(x_q) 
 *                       fe_values.shape_value(j, q_point)) * // phi_j(x_q) 
 *                    fe_values.JxW(q_point));                // dx 
 * 
 *               cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q) 
 *                               rhs_values[q_point] *               // f(x_q) 
 *                               fe_values.JxW(q_point));            // dx 
 *             } 
 * 
 * @endcode
 * 
 * 然后是右手边的第二项，即等高线积分。首先我们要找出这个单元格的面与边界部分Gamma2的交点是否为非零。为此，我们对所有面进行循环，检查其边界指示器是否等于 <code>1</code> ，这是我们在下面的 <code>run()</code> 函数中为组成Gamma2的边界部分指定的值。(边界指示器的默认值是 <code>0</code> ，所以只有在我们明确设置的情况下，面的指示器才能等于 <code>1</code> 。)
 * 

 * 
 * 
 * @code
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary() && (face->boundary_id() == 1)) 
 *             { 
 * 
 * @endcode
 * 
 * 如果我们来到这里，那么我们已经找到了一个属于Gamma2的外部面。接下来，我们必须计算形状函数的值和其他数量，这些都是我们在计算轮廓积分时需要的。这是用 <code>reinit</code> 函数完成的，我们已经从FEValue类中知道了。
 * 

 * 
 * 
 * @code
 *               fe_face_values.reinit(cell, face); 
 * 
 * @endcode
 * 
 * 然后，我们可以通过在所有的正交点上进行循环来进行积分。        在每个正交点上，我们首先计算法线导数的值。我们使用精确解的梯度和从 <code>fe_face_values</code> 对象中获得的当前正交点处的面的法向量来进行计算。然后用它来计算这个面对右手边的额外贡献。
 * 

 * 
 * 
 * @code
 *               for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                    ++q_point) 
 *                 { 
 *                   const double neumann_value = 
 *                     (exact_solution.gradient( 
 *                        fe_face_values.quadrature_point(q_point)) * 
 *                      fe_face_values.normal_vector(q_point)); 
 * 
 *                   for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                     cell_rhs(i) += 
 *                       (fe_face_values.shape_value(i, q_point) * // phi_i(x_q) 
 *                        neumann_value *                          // g(x_q) 
 *                        fe_face_values.JxW(q_point));            // dx 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 现在我们有了本单元的贡献，我们可以把它转移到全局矩阵和右手边的向量，就像之前的例子一样。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           { 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               system_matrix.add(local_dof_indices[i], 
 *                                 local_dof_indices[j], 
 *                                 cell_matrix(i, j)); 
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 同样，对边界值的消除和处理也在前面显示过。
 * 

 * 
 * 然而，我们注意到，现在我们插值边界值的边界指标（由 <code>interpolate_boundary_values</code> 的第二个参数表示）不再代表整个边界了。相反，它是我们没有指定其他指标的那部分边界（见下文）。因此，边界上不属于Gamma1的自由度被排除在边界值的插值之外，就像我们希望的那样。
 * 

 * 
 * 
 * @code
 *     hanging_node_constraints.condense(system_matrix); 
 *     hanging_node_constraints.condense(system_rhs); 
 * 
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              Solution<dim>(), 
 *                                              boundary_values); 
 *     MatrixTools::apply_boundary_values(boundary_values, 
 *                                        system_matrix, 
 *                                        solution, 
 *                                        system_rhs); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemsolve"></a> 
 * <h4>HelmholtzProblem::solve</h4>
 * 

 * 
 * 解方程组的方法与之前一样。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::solve() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-12); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *     hanging_node_constraints.distribute(solution); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemrefine_grid"></a> 
 * <h4>HelmholtzProblem::refine_grid</h4>
 * 

 * 
 * 现在是做网格细化的函数。根据传递给构造函数的细化模式，我们进行全局或适应性细化。
 * 

 * 
 * 全局细化很简单，所以没有什么可评论的。 在适应性细化的情况下，我们使用的函数和类与前面的例子程序相同。请注意，我们可以将诺伊曼边界与迪里切特边界区别对待，事实上在这里也应该这样做，因为我们在部分边界上有诺伊曼边界条件，但是由于我们在这里没有描述诺伊曼值的函数（我们只是在组装矩阵时从精确解中构造这些值），我们省略了这个细节，尽管以严格正确的方式做这些并不难添加。
 * 

 * 
 * 在开关的最后，我们有一个看起来稍微有点奇怪的默认情况：一个 <code>Assert</code> statement with a <code>false</code> 条件。由于 <code>Assert</code> 宏在条件为假的时候会引发一个错误，这意味着只要我们碰到这个语句，程序就会被中止。这是故意的。现在我们只实现了两种细化策略（全局性和适应性），但有人可能想增加第三种策略（例如，具有不同细化标准的适应性），并在决定细化模式的枚举中增加第三个成员。如果不是switch语句的默认情况，这个函数会简单地运行到结束而不做任何事情。这很可能不是原意。因此，在deal.II库中，你会发现一个防御性的编程技术，那就是总是有默认的中止案例，以确保在switch语句中列出案例时没有考虑的值最终被抓住，并迫使程序员添加代码来处理它们。我们还将在下面的其他地方使用同样的技术。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::refine_grid() 
 *   { 
 *     switch (refinement_mode) 
 *       { 
 *         case global_refinement: 
 *           { 
 *             triangulation.refine_global(1); 
 *             break; 
 *           } 
 * 
 *         case adaptive_refinement: 
 *           { 
 *             Vector<float> estimated_error_per_cell( 
 *               triangulation.n_active_cells()); 
 * 
 *             KellyErrorEstimator<dim>::estimate( 
 *               dof_handler, 
 *               QGauss<dim - 1>(fe->degree + 1), 
 *               std::map<types::boundary_id, const Function<dim> *>(), 
 *               solution, 
 *               estimated_error_per_cell); 
 * 
 *             GridRefinement::refine_and_coarsen_fixed_number( 
 *               triangulation, estimated_error_per_cell, 0.3, 0.03); 
 * 
 *             triangulation.execute_coarsening_and_refinement(); 
 * 
 *             break; 
 *           } 
 * 
 *         default: 
 *           { 
 *             Assert(false, ExcNotImplemented()); 
 *           } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemprocess_solution"></a> 
 * <h4>HelmholtzProblem::process_solution</h4>
 * 

 * 
 * 最后，我们想在计算出解决方案后对其进行处理。为此，我们用各种（半）准则对误差进行积分，并生成表格，这些表格以后将被用来以漂亮的格式显示对连续解的收敛情况。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::process_solution(const unsigned int cycle) 
 *   { 
 * 
 * @endcode
 * 
 * 我们的第一个任务是计算误差准则。为了整合计算出的数值解和连续解之间的差异（由本文件顶部定义的Solution类描述），我们首先需要一个向量来保存每个单元的误差准则。由于16位数的精度对这些数量来说并不那么重要，我们通过使用 <code>float</code> 而不是 <code>double</code> 值来节省一些内存。
 * 

 * 
 * 下一步是使用库中的一个函数来计算每个单元的L2准则的误差。 我们必须将DoF处理程序对象、保存数值解的节点值的向量、作为函数对象的连续解、它应将每个单元上的误差规范放入的向量、计算该规范的正交规则，以及要使用的规范类型传递给它。这里，我们使用高斯公式，在每个空间方向上有三个点，并计算L2规范。
 * 

 * 
 * 最后，我们想得到全局L2准则。这当然可以通过对每个单元格上的规范的平方求和，然后取该值的平方根来得到。这相当于取每个单元格上的规范向量的l2（小写 <code>l</code>  ）规范。
 * 

 * 
 * 
 * @code
 *     Vector<float> difference_per_cell(triangulation.n_active_cells()); 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe->degree + 1), 
 *                                       VectorTools::L2_norm); 
 *     const double L2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 * @endcode
 * 
 * 通过同样的程序，我们可以得到H1半正态。我们重新使用 <code>difference_per_cell</code> 向量，因为在计算了上面的 <code>L2_error</code> 变量后，它不再被使用。全局 $H^1$ 半正态误差的计算方法是：取每个单元格上的误差的平方和，然后取其平方根--这个操作由 VectorTools::compute_global_error. 方便地执行。
 * 
 * @code
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe->degree + 1), 
 *                                       VectorTools::H1_seminorm); 
 *     const double H1_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::H1_seminorm); 
 * 
 * @endcode
 * 
 * 最后，我们计算出最大法线。当然，我们实际上不能计算域中*所有*点上的真正的最大误差，而只能计算有限的评估点上的最大误差，为了方便起见，我们仍然称之为 "正交点"，并用一个正交类型的对象来表示，尽管我们实际上没有进行任何积分。
 * 

 * 
 * 然后是我们想在哪些点上精确地进行评估的问题。事实证明，我们得到的结果相当敏感地取决于所使用的 "正交 "点。还有一个超融合的问题。在某些网格上，对于多项式程度 $k\ge 2$ ，有限元解决方案在节点点以及Gauss-Lobatto点上特别精确，比随机选择的点要精确得多。(参见 @cite Li2019 和第1.2节的讨论和参考文献，以了解更多这方面的信息)。换句话说，如果我们有兴趣找到最大的差值 $u(\mathbf x)-u_h(\mathbf x)$ ，那么我们应该看一下 $\mathbf x$ ，这些点特别不属于这种 "特殊 "的点，而且我们特别不应该用`QGauss(fe->degree+1)`来定义我们评估的地方。相反，我们使用一个特殊的正交规则，该规则是通过梯形规则迭代有限元的度数乘以2再加上每个空间方向的1而得到的。请注意，QIterated类的构造函数需要一个一维正交规则和一个数字，这个数字告诉它在每个空间方向重复这个规则的频率。
 * 

 * 
 * 使用这个特殊的正交规则，我们就可以尝试找到每个单元的最大误差。最后，我们通过调用 VectorTools::compute_global_error. 来计算每个单元上的L无穷大误差的全局L无穷大误差。
 * 
 * @code
 *     const QTrapezoid<1>  q_trapez; 
 *     const QIterated<dim> q_iterated(q_trapez, fe->degree * 2 + 1); 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       q_iterated, 
 *                                       VectorTools::Linfty_norm); 
 *     const double Linfty_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::Linfty_norm); 
 * 
 * @endcode
 * 
 * 在所有这些错误被计算出来之后，我们最终写出一些输出。此外，我们通过指定列的键和值将重要的数据添加到TableHandler中。 注意，没有必要事先定义列的键 -- 只需添加值即可，列将按照第一次添加值的顺序被引入到表中。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_active_cells = triangulation.n_active_cells(); 
 *     const unsigned int n_dofs         = dof_handler.n_dofs(); 
 * 
 *     std::cout << "Cycle " << cycle << ':' << std::endl 
 *               << "   Number of active cells:       " << n_active_cells 
 *               << std::endl 
 *               << "   Number of degrees of freedom: " << n_dofs << std::endl; 
 * 
 *     convergence_table.add_value("cycle", cycle); 
 *     convergence_table.add_value("cells", n_active_cells); 
 *     convergence_table.add_value("dofs", n_dofs); 
 *     convergence_table.add_value("L2", L2_error); 
 *     convergence_table.add_value("H1", H1_error); 
 *     convergence_table.add_value("Linfty", Linfty_error); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HelmholtzProblemrun"></a> 
 * <h4>HelmholtzProblem::run</h4>
 * 

 * 
 * 和前面的例子程序一样， <code>run</code> 函数控制执行的流程。基本布局与前面的例子一样：在连续细化的网格上有一个外循环，在这个循环中首先是问题的设置，组装线性系统，求解，和后处理。
 * 

 * 
 * 主循环的第一个任务是创建和细化网格。这和前面的例子一样，唯一的区别是我们想把边界的一部分标记为诺伊曼型，而不是迪里希型。
 * 

 * 
 * 为此，我们将使用以下惯例。属于Gamma1的面将有边界指示器 <code>0</code> （这是默认的，所以我们不需要明确设置），属于Gamma2的面将使用 <code>1</code> 作为边界指示器。 为了设置这些值，我们在所有单元格上循环，然后在给定单元格的所有面上循环，检查它是否是我们想用Gamma2表示的边界的一部分，如果是，则将其边界指示器设置为 <code>1</code>  。在本程序中，我们认为左边和底部的边界是Gamma2。我们通过询问一个面的中点的x或y坐标（即向量分量0和1）是否等于-1来确定一个面是否是该边界的一部分，但我们必须给出一些小的回旋余地，因为比较在中间计算中会有四舍五入的浮点数是不稳定的。
 * 

 * 
 * 值得注意的是，我们必须在这里对所有的单元格进行循环，而不仅仅是活动单元格。原因是在细化时，新创建的面会继承其父面的边界指标。如果我们现在只设置活动面的边界指示器，粗化一些单元并在以后细化它们，它们将再次拥有我们没有修改的父单元的边界指示器，而不是我们想要的那个。因此，我们必须改变Gamma2上所有单元的面的边界指标，无论它们是否处于活动状态。另外，我们当然也可以在最粗的网格上完成这项工作（即在第一个细化步骤之前），之后才细化网格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HelmholtzProblem<dim>::run() 
 *   { 
 *     const unsigned int n_cycles = 
 *       (refinement_mode == global_refinement) ? 5 : 9; 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
 *       { 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation, -1., 1.); 
 *             triangulation.refine_global(3); 
 * 
 *             for (const auto &cell : triangulation.cell_iterators()) 
 *               for (const auto &face : cell->face_iterators()) 
 *                 { 
 *                   const auto center = face->center(); 
 *                   if ((std::fabs(center(0) - (-1.0)) < 1e-12) || 
 *                       (std::fabs(center(1) - (-1.0)) < 1e-12)) 
 *                     face->set_boundary_id(1); 
 *                 } 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 * @endcode
 * 
 * 接下来的步骤在前面的例子中已经知道了。这主要是每个有限元程序的基本设置。
 * 

 * 
 * 
 * @code
 *         setup_system(); 
 * 
 *         assemble_system(); 
 *         solve(); 
 * 
 * @endcode
 * 
 * 在这一连串的函数调用中，最后一步通常是对自己感兴趣的数量的计算解进行评估。这在下面的函数中完成。由于该函数产生的输出显示了当前细化步骤的编号，我们将这个编号作为一个参数传递。
 * 

 * 
 * 
 * @code
 *         process_solution(cycle); 
 *       } 
 * @endcode
 * 
 * 
 * <a name="Outputofgraphicaldata"></a> 
 * <h5>Output of graphical data</h5>
 * 

 * 
 * 在最后一次迭代后，我们在最细的网格上输出解决方案。这是用下面的语句序列完成的，我们在以前的例子中已经讨论过了。第一步是生成一个合适的文件名（这里称为 <code>vtk_filename</code> ，因为我们想以VTK格式输出数据；我们添加前缀以区分该文件名与下面其他输出文件的文件名）。在这里，我们通过网格细化算法来增加名称，和上面一样，我们要确保在增加了另一种细化方法而没有通过下面的switch语句来处理的情况下，中止程序。
 * 

 * 
 * 
 * @code
 *     std::string vtk_filename; 
 *     switch (refinement_mode) 
 *       { 
 *         case global_refinement: 
 *           vtk_filename = "solution-global"; 
 *           break; 
 *         case adaptive_refinement: 
 *           vtk_filename = "solution-adaptive"; 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 * @endcode
 * 
 * 我们用一个后缀来增加文件名，表示我们在计算中使用的有限元。为此，有限元基类将每个坐标变量中形状函数的最大多项式程度存储为一个变量 <code>degree</code> ，我们在切换语句中使用（注意，双线性形状函数的多项式程度实际上是2，因为它们包含术语 <code>x*y</code> ；但是，每个坐标变量的多项式程度仍然只有1）。我们再次使用同样的防御性编程技术来防止多项式阶数具有意外值的情况，在switch语句的默认分支中使用 <code>Assert (false, ExcNotImplemented())</code> 这个成语。
 * 

 * 
 * 
 * @code
 *     switch (fe->degree) 
 *       { 
 *         case 1: 
 *           vtk_filename += "-q1"; 
 *           break; 
 *         case 2: 
 *           vtk_filename += "-q2"; 
 *           break; 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 * @endcode
 * 
 * 一旦我们有了输出文件的基本名称，我们就为VTK输出添加一个合适的扩展名，打开一个文件，并将解决方案的向量添加到将进行实际输出的对象中。
 * 

 * 
 * 
 * @code
 *     vtk_filename += ".vtk"; 
 *     std::ofstream output(vtk_filename); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 * 
 * @endcode
 * 
 * 现在像以前一样建立中间格式是下一步。我们在这里再介绍一下deal.II的一个特点。其背景如下：在这个函数的一些运行中，我们使用了双二次元的有限元。然而，由于几乎所有的输出格式都只支持双线性数据，所以数据只写成了双线性，信息因此而丢失。 当然，我们不能改变图形程序接受其输入的格式，但我们可以用不同的方式来写数据，这样我们就能更接近于四次方近似中的信息。例如，我们可以把每个单元写成四个子单元，每个子单元都有双线数据，这样我们在三角图中的每个单元都有九个数据点。当然，图形程序显示的这些数据仍然只是双线性的，但至少我们又给出了一些我们拥有的信息。
 * 

 * 
 * 为了允许在每个实际单元中写入多个子单元， <code>build_patches</code> 函数接受一个参数（默认为 <code>1</code>  ，这就是为什么你在之前的例子中没有看到这个参数）。这个参数表示每个空间方向上的每个单元应被细分为多少个子单元来输出。例如，如果你给出  <code>2</code>  ，这将导致二维的4个单元和三维的8个单元。对于二次元元素，每个空间方向的两个子单元显然是正确的选择，所以这就是我们所选择的。一般来说，对于多项式阶的元素 <code>q</code>, we use <code>q</code> 细分，元素的顺序也是按照上述方式确定的。
 * 

 * 
 * 有了这样生成的中间格式，我们就可以实际写入图形输出了。
 * 

 * 
 * 
 * @code
 *     data_out.build_patches(fe->degree); 
 *     data_out.write_vtk(output); 
 * @endcode
 * 
 * 
 * <a name="Outputofconvergencetables"></a> 
 * <h5>Output of convergence tables</h5>
 * 

 * 
 * 在图形输出之后，我们还想从我们在  <code>process_solution</code>  中进行的误差计算中生成表格。在那里，我们用每个细化步骤的单元格数量以及不同规范的误差来填充一个表格对象。
 * 

 * 
 * 为了使这些数据有更好的文本输出，我们可能想设置输出时写入数值的精度。我们使用3位数，这对误差规范来说通常是足够的。默认情况下，数据是以定点符号写入的。然而，对于人们想看到的科学符号的列，另一个函数调用设置了 <code>scientific_flag</code> to <code>true</code>  ，导致数字的浮点表示。
 * 

 * 
 * 
 * @code
 *     convergence_table.set_precision("L2", 3); 
 *     convergence_table.set_precision("H1", 3); 
 *     convergence_table.set_precision("Linfty", 3); 
 * 
 *     convergence_table.set_scientific("L2", true); 
 *     convergence_table.set_scientific("H1", true); 
 *     convergence_table.set_scientific("Linfty", true); 
 * 
 * @endcode
 * 
 * 对于输出到LaTeX文件的表格，默认的列的标题是作为参数给 <code>add_value</code> 函数的键。要想拥有不同于默认的TeX标题，你可以通过以下函数调用来指定它们。注意，`\\'被编译器简化为`\'，这样，真正的TeX标题就是，例如，` $L^\infty$  -error'。
 * 

 * 
 * 
 * @code
 *     convergence_table.set_tex_caption("cells", "\\# cells"); 
 *     convergence_table.set_tex_caption("dofs", "\\# dofs"); 
 *     convergence_table.set_tex_caption("L2", "L^2-error"); 
 *     convergence_table.set_tex_caption("H1", "H^1-error"); 
 *     convergence_table.set_tex_caption("Linfty", "L^\\infty-error"); 
 * 
 * @endcode
 * 
 * 最后，表格中每一列的默认LaTeX格式是`c'（居中）。要指定一个不同的（如`右'），可以使用以下函数。
 * 

 * 
 * 
 * @code
 *     convergence_table.set_tex_format("cells", "r"); 
 *     convergence_table.set_tex_format("dofs", "r"); 
 * 
 * @endcode
 * 
 * 在这之后，我们终于可以把表写到标准输出流 <code>std::cout</code> （在多写一行空行之后，使事情看起来更漂亮）。请注意，文本格式的输出是非常简单的，标题可能不会直接打印在特定的列上面。
 * 

 * 
 * 
 * @code
 *     std::cout << std::endl; 
 *     convergence_table.write_text(std::cout); 
 * 
 * @endcode
 * 
 * 该表也可以写成LaTeX文件。 在调用 "latex filename "和例如 "xdvi filename "后，可以查看（很好的）格式化的表格，其中filename是我们现在要写入输出的文件名。我们构建文件名的方法和以前一样，但有一个不同的前缀 "error"。
 * 

 * 
 * 
 * @code
 *     std::string error_filename = "error"; 
 *     switch (refinement_mode) 
 *       { 
 *         case global_refinement: 
 *           error_filename += "-global"; 
 *           break; 
 *         case adaptive_refinement: 
 *           error_filename += "-adaptive"; 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     switch (fe->degree) 
 *       { 
 *         case 1: 
 *           error_filename += "-q1"; 
 *           break; 
 *         case 2: 
 *           error_filename += "-q2"; 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     error_filename += ".tex"; 
 *     std::ofstream error_table_file(error_filename); 
 * 
 *     convergence_table.write_tex(error_table_file); 
 * @endcode
 * 
 * 
 * <a name="Furthertablemanipulations"></a> 
 * <h5>Further table manipulations</h5>
 * 

 * 
 * 在全局细化的情况下，输出收敛率也可能是有意义的。这可以通过ConvergenceTable提供的比常规TableHandler的功能来实现。然而，我们只为全局细化做这件事，因为对于自适应细化来说，确定像收敛顺序这样的事情是比较麻烦的。在此，我们还展示了一些可以用表来做的其他事情。
 * 

 * 
 * 
 * @code
 *     if (refinement_mode == global_refinement) 
 *       { 
 * 
 * @endcode
 * 
 * 第一件事是，人们可以将单个列组合在一起，形成所谓的超级列。从本质上讲，这些列保持不变，但被分组的那些列将得到一个贯穿一组中所有列的标题。例如，让我们把 "周期 "和 "单元格 "两列合并成一个名为 "n单元格 "的超级列。
 * 

 * 
 * 
 * @code
 *         convergence_table.add_column_to_supercolumn("cycle", "n cells"); 
 *         convergence_table.add_column_to_supercolumn("cells", "n cells"); 
 * 
 * @endcode
 * 
 * 接下来，没有必要总是输出所有的列，或者按照它们在运行过程中最初添加的顺序。选择和重新排列列的工作方式如下（注意，这包括超级列）。
 * 

 * 
 * 
 * @code
 *         std::vector<std::string> new_order; 
 *         new_order.emplace_back("n cells"); 
 *         new_order.emplace_back("H1"); 
 *         new_order.emplace_back("L2"); 
 *         convergence_table.set_column_order(new_order); 
 * 
 * @endcode
 * 
 * 对于在这之前发生在ConvergenceTable上的一切，使用一个简单的TableHandler就足够了。事实上，ConvergenceTable是由TableHandler派生出来的，但它提供了自动评估收敛率的额外功能。例如，下面是我们如何让表计算减少率和收敛率（收敛率是减少率的二进制对数）。
 * 

 * 
 * 
 * @code
 *         convergence_table.evaluate_convergence_rates( 
 *           "L2", ConvergenceTable::reduction_rate); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "L2", ConvergenceTable::reduction_rate_log2); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "H1", ConvergenceTable::reduction_rate); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "H1", ConvergenceTable::reduction_rate_log2); 
 * 
 * @endcode
 * 
 * 这些函数的每一次调用都会产生一个额外的列，与原来的列（在我们的例子中是 "L2 "和 "H1 "列）合并成一个超级列。
 * 

 * 
 * 最后，我们想再次写下这个收敛图，首先写到屏幕上，然后以LaTeX格式写到磁盘上。文件名还是按照上面的方法构建。
 * 

 * 
 * 
 * @code
 *         std::cout << std::endl; 
 *         convergence_table.write_text(std::cout); 
 * 
 *         std::string conv_filename = "convergence"; 
 *         switch (refinement_mode) 
 *           { 
 *             case global_refinement: 
 *               conv_filename += "-global"; 
 *               break; 
 *             case adaptive_refinement: 
 *               conv_filename += "-adaptive"; 
 *               break; 
 *             default: 
 *               Assert(false, ExcNotImplemented()); 
 *           } 
 *         switch (fe->degree) 
 *           { 
 *             case 1: 
 *               conv_filename += "-q1"; 
 *               break; 
 *             case 2: 
 *               conv_filename += "-q2"; 
 *               break; 
 *             default: 
 *               Assert(false, ExcNotImplemented()); 
 *           } 
 *         conv_filename += ".tex"; 
 * 
 *         std::ofstream table_file(conv_filename); 
 *         convergence_table.write_tex(table_file); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 在进入 <code>main()</code> 之前的最后一步是关闭命名空间 <code>Step7</code> ，我们已经把这个程序所需要的一切都放在这个命名空间里。
 * 

 * 
 * 
 * @code
 * } // namespace Step7 
 * @endcode
 * 
 * 
 * <a name="Mainfunction"></a> 
 * <h3>Main function</h3>
 * 

 * 
 * 主函数主要和以前一样。唯一不同的是，我们解了三次，一次是Q1和适应性细化，一次是Q1元素和全局细化，一次是Q2元素和全局细化。
 * 

 * 
 * 由于我们在下面为两个空间维度实例化了几个模板类，我们通过在函数的开头声明一个常数来表示空间维度的数量，使之更加通用。如果你想在1d或2d中运行程序，那么你只需要改变这个实例，而不是下面的所有用法。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   const unsigned int dim = 2; 
 * 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step7; 
 * 
 * @endcode
 * 
 * 现在是对主类的三次调用。每个调用都被封锁在大括号中，以便在区块结束时和我们进入下一个运行之前销毁各自的对象（即有限元和HelmholtzProblem对象）。这就避免了变量名称的冲突，也确保了在三次运行中的一次运行结束后立即释放内存，而不是只在 <code>try</code> 块的末尾释放。
 * 

 * 
 * 
 * @code
 *       { 
 *         std::cout << "Solving with Q1 elements, adaptive refinement" 
 *                   << std::endl 
 *                   << "=============================================" 
 *                   << std::endl 
 *                   << std::endl; 
 * 
 *         FE_Q<dim>             fe(1); 
 *         HelmholtzProblem<dim> helmholtz_problem_2d( 
 *           fe, HelmholtzProblem<dim>::adaptive_refinement); 
 * 
 *         helmholtz_problem_2d.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 * 
 *       { 
 *         std::cout << "Solving with Q1 elements, global refinement" << std::endl 
 *                   << "===========================================" << std::endl 
 *                   << std::endl; 
 * 
 *         FE_Q<dim>             fe(1); 
 *         HelmholtzProblem<dim> helmholtz_problem_2d( 
 *           fe, HelmholtzProblem<dim>::global_refinement); 
 * 
 *         helmholtz_problem_2d.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 * 
 *       { 
 *         std::cout << "Solving with Q2 elements, global refinement" << std::endl 
 *                   << "===========================================" << std::endl 
 *                   << std::endl; 
 * 
 *         FE_Q<dim>             fe(2); 
 *         HelmholtzProblem<dim> helmholtz_problem_2d( 
 *           fe, HelmholtzProblem<dim>::global_refinement); 
 * 
 *         helmholtz_problem_2d.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *       { 
 *         std::cout << "Solving with Q2 elements, adaptive refinement" 
 *                   << std::endl 
 *                   << "===========================================" << std::endl 
 *                   << std::endl; 
 * 
 *         FE_Q<dim>             fe(2); 
 *         HelmholtzProblem<dim> helmholtz_problem_2d( 
 *           fe, HelmholtzProblem<dim>::adaptive_refinement); 
 * 
 *         helmholtz_problem_2d.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *     } 
 *   catch (std::exception &exc) 
 *     { 
 *       std::cerr << std::endl 
 *                 << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       std::cerr << "Exception on processing: " << std::endl 
 *                 << exc.what() << std::endl 
 *                 << "Aborting!" << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       return 1; 
 *     } 
 *   catch (...) 
 *     { 
 *       std::cerr << std::endl 
 *                 << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       std::cerr << "Unknown exception!" << std::endl 
 *                 << "Aborting!" << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       return 1; 
 *     } 
 * 
 *   return 0; 
 * } 
 * 
 * @endcode
examples/step-7/doc/results.dox



<a name="Results"></a><h1>Results</h1>



该程序产生了两种输出。第一种是输出文件  <code>solution-adaptive-q1.vtk</code>  ,  <code>solution-global-q1.vtk</code>  , 和  <code>solution-global-q2.vtk</code>  。我们在此展示后者的三维视图。


 <img src="https://www.dealii.org/images/steps/developer/step-7.solution.png" alt=""> 





其次，该程序不仅将表格写入磁盘，而且在运行时也写入屏幕。输出结果如下（记得标为" <code>H1</code> "的列实际上显示的是错误的 $H^1$ <i>semi-</i>规范，而不是完整的 $H^1$ 规范）。




@code
examples/\step-7> make run
Solving with Q1 elements, adaptive refinement
=============================================


Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 81
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 157
Cycle 2:
   Number of active cells:       280
   Number of degrees of freedom: 341
Cycle 3:
   Number of active cells:       577
   Number of degrees of freedom: 690
Cycle 4:
   Number of active cells:       1099
   Number of degrees of freedom: 1264
Cycle 5:
   Number of active cells:       2191
   Number of degrees of freedom: 2452
Cycle 6:
   Number of active cells:       4165
   Number of degrees of freedom: 4510
Cycle 7:
   Number of active cells:       7915
   Number of degrees of freedom: 8440
Cycle 8:
   Number of active cells:       15196
   Number of degrees of freedom: 15912


cycle cells dofs     L2        H1      Linfty
    0    64    81 1.840e+00 2.858e+00 1.835e+00
    1   124   157 5.190e-02 1.200e+00 1.344e-01
    2   280   341 1.439e-02 7.892e-01 7.554e-02
    3   577   690 8.627e-03 5.061e-01 2.805e-02
    4  1099  1264 3.217e-03 3.030e-01 1.073e-02
    5  2191  2452 1.445e-03 2.097e-01 5.073e-03
    6  4165  4510 8.387e-04 1.460e-01 2.013e-03
    7  7915  8440 7.051e-04 1.053e-01 1.804e-03
    8 15196 15912 2.774e-04 7.463e-02 6.911e-04


Solving with Q1 elements, global refinement
===========================================


Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 81
Cycle 1:
   Number of active cells:       256
   Number of degrees of freedom: 289
Cycle 2:
   Number of active cells:       1024
   Number of degrees of freedom: 1089
Cycle 3:
   Number of active cells:       4096
   Number of degrees of freedom: 4225
Cycle 4:
   Number of active cells:       16384
   Number of degrees of freedom: 16641


cycle cells dofs     L2        H1      Linfty
    0    64    81 1.840e+00 2.858e+00 1.835e+00
    1   256   289 3.570e-02 1.199e+00 1.307e-01
    2  1024  1089 1.192e-02 7.565e-01 7.168e-02
    3  4096  4225 3.047e-03 3.823e-01 2.128e-02
    4 16384 16641 7.660e-04 1.917e-01 5.554e-03


n cells         H1                   L2
0    64 2.858e+00    -    - 1.840e+00     -    -
1   256 1.199e+00 2.38 1.25 3.570e-02 51.54 5.69
2  1024 7.565e-01 1.58 0.66 1.192e-02  2.99 1.58
3  4096 3.823e-01 1.98 0.98 3.047e-03  3.91 1.97
4 16384 1.917e-01 1.99 1.00 7.660e-04  3.98 1.99


Solving with Q2 elements, global refinement
===========================================


Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 289
Cycle 1:
   Number of active cells:       256
   Number of degrees of freedom: 1089
Cycle 2:
   Number of active cells:       1024
   Number of degrees of freedom: 4225
Cycle 3:
   Number of active cells:       4096
   Number of degrees of freedom: 16641
Cycle 4:
   Number of active cells:       16384
   Number of degrees of freedom: 66049


cycle cells dofs     L2        H1      Linfty
    0    64   289 1.606e-01 1.278e+00 3.029e-01
    1   256  1089 7.638e-03 5.248e-01 4.816e-02
    2  1024  4225 8.601e-04 1.086e-01 4.827e-03
    3  4096 16641 1.107e-04 2.756e-02 7.802e-04
    4 16384 66049 1.393e-05 6.915e-03 9.971e-05


n cells         H1                   L2
0    64 1.278e+00    -    - 1.606e-01     -    -
1   256 5.248e-01 2.43 1.28 7.638e-03 21.03 4.39
2  1024 1.086e-01 4.83 2.27 8.601e-04  8.88 3.15
3  4096 2.756e-02 3.94 1.98 1.107e-04  7.77 2.96
4 16384 6.915e-03 3.99 1.99 1.393e-05  7.94 2.99


Solving with Q2 elements, adaptive refinement
===========================================


Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 289
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 577
Cycle 2:
   Number of active cells:       289
   Number of degrees of freedom: 1353
Cycle 3:
   Number of active cells:       547
   Number of degrees of freedom: 2531
Cycle 4:
   Number of active cells:       1057
   Number of degrees of freedom: 4919
Cycle 5:
   Number of active cells:       2059
   Number of degrees of freedom: 9223
Cycle 6:
   Number of active cells:       3913
   Number of degrees of freedom: 17887
Cycle 7:
   Number of active cells:       7441
   Number of degrees of freedom: 33807
Cycle 8:
   Number of active cells:       14212
   Number of degrees of freedom: 64731


cycle cells dofs     L2        H1      Linfty
    0    64   289 1.606e-01 1.278e+00 3.029e-01
    1   124   577 7.891e-03 5.256e-01 4.852e-02
    2   289  1353 1.070e-03 1.155e-01 4.868e-03
    3   547  2531 5.962e-04 5.101e-02 1.876e-03
    4  1057  4919 1.977e-04 3.094e-02 7.923e-04
    5  2059  9223 7.738e-05 1.974e-02 7.270e-04
    6  3913 17887 2.925e-05 8.772e-03 1.463e-04
    7  7441 33807 1.024e-05 4.121e-03 8.567e-05
    8 14212 64731 3.761e-06 2.108e-03 2.167e-05
@endcode




我们可以看到网格细化后的误差减少，对于进行全局细化的情况，也可以看到收敛率。可以清楚地看到Q1和Q2元素在 $H^1$ 半规范下的线性和二次收敛率，以及 $L_2$ 规范下的二次和三次收敛率。





最后，该程序还生成了LaTeX版本的表格（此处未显示），这些表格被写入一个文件中，以便可以复制粘贴到LaTeX文件中。




<a name="Whenistheerrorsmall"></a><h4> When is the error "small"? </h4>


我们上面展示的是如何在一些不同的规范中确定误差 $\|u-u_h\|$ 的大小。我们这样做主要是因为我们对测试我们的解决方案是否*融合感兴趣。但是从工程的角度来看，这个问题往往更实际：我的网格要做得多细才能使误差 "足够小"？换句话说，如果在上表中， $H^1$ 的半规范已经减少到`4.121e-03`，这是否足以让我在蓝图上签字，宣布我们的数值模拟显示桥梁足够坚固？

在实践中，我们很少遇到这种情况，因为我通常不能在重要的情况下将数值解 $u_h$ 与精确解 $u$ 进行比较--如果我知道 $u$ ，我就不必计算 $u_h$  。但即使我可以，一般情况下要问的问题是。`4.121e-03`是什么？解决方案将有物理单位，例如公斤-米-平方，我在域上积分一个单位为上述平方的函数，然后取其平方根。因此，如果域是二维的， $\|u-u_h\|_{L_2}$ 的单位是公斤-米-立方。那么问题来了。 $4.121\times 10^{-3}$ 的单位是kg-times-mubed小吗？这取决于你要模拟的是什么。如果你是一个天文学家，习惯于以太阳质量为单位的质量和以光年为单位的距离，那么是的，这是一个小得惊人的数字。但是如果你是做原子物理的，那就不是：这并不小，而且你的误差肯定也不够小；你需要一个更细的网格。

换句话说，当我们看这些数字的时候，我们通常需要与一个 "尺度 "进行比较。做到这一点的一个方法是，不看任何规范中的*绝对*误差 $\|u-u_h\|$ ，而是看*相对*误差 $\|u-u_h\|/\|u\|$  。如果这个比率是 $10^{-5}$ ，那么你就知道，*平均而言， $u$ 和 $u_h$ 之间的差异是0.001%--对于工程而言可能足够小。

我们如何计算 $\|u\|$ ？我们只需要在所有的单元格上做一个积分循环，在这些单元格上做正交点，然后把东西加起来，最后取平方根。但有一个更简单的方法经常使用。你可以调用

@code
    Vector<double> zero_vector (dof_handler.n_dofs());
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      zero_vector,
                                      Solution<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
@endcode

来计算  $\|u-0\|_{L_2}$  。另外，如果你特别懒，不喜欢创建`零_向量'，你可以使用，如果网格不是太粗，那么 $\|u\| \approx \|u_h\|$  ，我们可以通过调用 $\|u\| \approx \|u_h\|=\|0-u_h\|$ 来计算。

@code
    Vector<float> norm_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      ZeroFunction<dim>(),
                                      norm_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
@endcode

在这两种情况下，我们只需要像在程序中一样，通过调用以下命令将单元格准则的向量合并为一个全局准则

@code
    const double L2_norm =
      VectorTools::compute_global_error(triangulation,
                                        norm_per_cell,
                                        VectorTools::L2_norm);
@endcode






<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="HigherOrderElements"></a><h4> Higher Order Elements </h4>


继续运行具有高阶元素的程序（ $Q_3$  ,  $Q_4$  , ...）。你会注意到，代码的几个部分的断言将被触发（例如，在为数据输出生成文件名时）。你可能必须解决这些问题，但要使程序正常工作应该不是很难

<a name="ConvergenceComparison"></a><h4> Convergence Comparison </h4>


Q1或Q2更好吗？自适应细化与全局细化又如何呢？比较它们的一个（有些不公平，但很典型）指标是看误差与未知数的关系。

要看到这一点，以对数风格创建一个图，未知数在 $x$ 轴上， $L_2$ 误差在 $y$ 轴上。你可以为 $h^2=N^{-1}$ 和 $h^3=N^{-3/2}$ 添加参考线，并检查全局和适应性细化是否遵循这些。如果我们做一个（并非完全不合理的）假设，即对于一个好的线性求解器，计算工作量与未知数 $N$ 成正比，那么很明显， ${\cal O}(N^{-3/2})$ 的误差减少比 ${\cal O}(N^{-1})$ 形式的减少要好得多：也就是说，与使用全局细化相比，自适应细化以更少的计算工作量给我们带来了期望的误差等级。这并不是一个特别令人惊讶的结论，但是值得在实践中检查这类假设。

当然，更公平的比较是在 $x$ 轴上绘制运行时间（先切换到发布模式！）而不是未知数的数量。如果你通过对每个细化步骤进行计时（例如，使用Timer类）来绘制运行时间与未知数数量的关系，你会发现线性求解器并不完美--其运行时间的增长速度超过了与线性系统大小成正比的速度--挑选一个更好的线性求解器可能适合这种比较。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-7.cc"
*/
