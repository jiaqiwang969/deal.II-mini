/**
@page step_29 The step-29 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemsetting">Problem setting</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeDirichletBoundaryValuescodeclass">The <code>DirichletBoundaryValues</code> class</a>
        <li><a href="#ThecodeParameterReadercodeclass">The <code>ParameterReader</code> class</a>
      <ul>
        <li><a href="#codeParameterReaderdeclare_parameterscode"><code>ParameterReader::declare_parameters</code></a>
        <li><a href="#codeParameterReaderread_parameterscode"><code>ParameterReader::read_parameters</code></a>
      </ul>
        <li><a href="#ThecodeComputeIntensitycodeclass">The <code>ComputeIntensity</code> class</a>
        <li><a href="#ThecodeUltrasoundProblemcodeclass">The <code>UltrasoundProblem</code> class</a>
      <ul>
        <li><a href="#codeUltrasoundProblemmake_gridcode"><code>UltrasoundProblem::make_grid</code></a>
        <li><a href="#codeUltrasoundProblemsetup_systemcode"><code>UltrasoundProblem::setup_system</code></a>
        <li><a href="#codeUltrasoundProblemassemble_systemcode"><code>UltrasoundProblem::assemble_system</code></a>
        <li><a href="#codeUltrasoundProblemsolvecode"><code>UltrasoundProblem::solve</code></a>
        <li><a href="#codeUltrasoundProblemoutput_resultscode"><code>UltrasoundProblem::output_results</code></a>
        <li><a href="#codeUltrasoundProblemruncode"><code>UltrasoundProblem::run</code></a>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-29/doc/intro.dox

 <br> 

<i>
This program was contributed by Moritz Allmaras at Texas A&amp;M
University. Some of the work on this tutorial program has been funded
by NSF under grant DMS-0604778.
</i>

<b>Note:</b> 为了运行这个程序，deal.II必须被配置为使用UMFPACK稀疏直接求解器。请参考<a
href="../../readme.html#umfpack">ReadMe</a>中的说明如何做到这一点。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



一个经常出现的问题是如何用deal.II解决涉及复值函数的问题。对于许多问题，与其直接使用复值有限元，不如将复值函数分成实部和虚部，并使用单独的标量有限元场来离散它们中的每一个，这样做往往更方便。基本上，这相当于把一个复值方程看作是两个实值方程的系统。这个简短的例子演示了如何在deal.II中通过使用 <code>FE_system</code> 对象来堆叠两个代表实部和虚部的有限元场来实现。(相反的方法，保持所有的复值，在另一个教程程序中演示：见步骤58)。当分成实部和虚部时，这里涉及的方程属于矢量值问题的范畴。这个主题的顶层概述可以在  @ref vector_valued  模块中找到。

除了这些讨论，我们还讨论了ParameterHandler类，它提供了一种方便的方法，在运行时从配置文件中读取参数，而不需要重新编译程序代码。




<a name="Problemsetting"></a><h3>Problem setting</h3>


这个程序的最初目的是模拟由一个几何形状可变的换能器透镜产生的超声波的聚焦特性。最近在医学成像方面的应用，不仅将超声波用于成像，而且还能激发材料中的某些局部效应，如光学特性的变化，然后可以用其他成像技术来测量。这些方法的一个重要因素是能够将超声波的强度集中在材料的一个特定部分，最好是一个点，以便能够检查该特定位置的材料特性。

为了推导出这个问题的模型，我们把超声波看成是由波浪方程支配的压力波。

@f[
	\frac{\partial^2 U}{\partial t^2}	-	c^2 \Delta U = 0


@f]

其中 $c$ 是波速（为简单起见，我们假设为常数）， $U
= U(x,t),\;x \in \Omega,\;t\in\mathrm{R}$  。边界 $\Gamma=\partial\Omega$ 分为两部分 $\Gamma_1$ 和 $\Gamma_2=\Gamma\setminus\Gamma_1$ ，其中 $\Gamma_1$ 代表换能器透镜， $\Gamma_2$ 是吸收边界（也就是说，我们希望在 $\Gamma_2$ 上选择边界条件，使其模仿一个更大的域）。在 $\Gamma_1$ 上，换能器产生一个恒定频率 ${\omega}>0$ 和恒定振幅（我们在这里选择为1）的波。

@f[
U(x,t) = \cos{\omega t}, \qquad x\in \Gamma_1


@f]



如果没有其他（内部或边界）源，并且由于唯一的源具有频率 $\omega$ ，那么解决方案可以接受形式为 $U(x,t) = \textrm{Re}\left(u(x)\,e^{i\omega
t})\right)$ 的变量分离。复值函数 $u(x)$ 描述了频率为 ${\omega}$ 的波的振幅和相位（相对于源）的空间依赖性，其中振幅是我们感兴趣的量。通过将这种形式的解决方案插入波浪方程，我们看到，对于 $u$ ，我们有

@f{eqnarray*}


-\omega^2 u(x) - c^2\Delta u(x) &=& 0, \qquad x\in\Omega,\\
u(x) &=& 1,  \qquad x\in\Gamma_1.


@f}



为了在 $\Gamma_2$ 上找到模拟吸收边界的合适条件，考虑一个频率为 ${\omega}$ 的波在 $k\in {\mathrm{R}^2}$ 方向上行驶。为了使 $V$ 能够解决波浪方程， $|k|={\frac{\omega}{c}}$ 必须成立。假设这个波以直角击中 $x_0\in\Gamma_2$ 的边界，即 $n=\frac{k}{|k|}$ ， $n$ 表示 $\Omega$ 在 $x_0$ 的外单位法线。然后在 $x_0$ ，这个波满足方程式

@f[
c (n\cdot\nabla V) + \frac{\partial V}{\partial t} = (i\, c\, |k| - i\, \omega) V = 0.


@f]

因此，通过强制执行边界条件

@f[
c (n\cdot\nabla U) + \frac{\partial U}{\partial t} = 0, \qquad x\in\Gamma_2,


@f]

以直角撞击边界 $\Gamma_2$ 的波将被完全吸收。另一方面，那些没有以直角撞击边界的波场部分不满足这个条件，将其作为边界条件强制执行会产生部分反射，即只有部分波会通过边界，就像它不在这里一样，而剩下的部分波会被反射回域中。

如果我们愿意接受这一点作为吸收边界的充分近似，我们最终得出以下问题，对于 $u$  。

@f{eqnarray*}


-\omega^2 u - c^2\Delta u &=& 0, \qquad x\in\Omega,\\
c (n\cdot\nabla u) + i\,\omega\,u &=&0, \qquad x\in\Gamma_2,\\
u &=& 1,  \qquad x\in\Gamma_1.


@f}

这是一个亥姆霍兹方程（类似于步骤7中的方程，但这次有''坏符号''），在 $\Gamma_1$ 上有Dirichlet数据，在 $\Gamma_2$ 上有混合边界条件。由于 $\Gamma_2$ 上的条件，我们不能只是分别处理 $u$ 的实部和虚部方程。然而，我们可以把 $u$ 的PDE看作是 $u$ 的实部和虚部的两个PDE系统， $\Gamma_2$ 的边界条件代表系统中两个部分之间的耦合项。这是按以下思路进行的。让 $v=\textrm{Re}\;u,\; w=\textrm{Im}\;u$ ，然后在 $v$ 和 $w$ 方面，我们有以下系统。

@f{eqnarray*}
  \left.\begin{array}{ccc}


    -\omega^2 v - c^2\Delta v &=& 0 \quad\\


    -\omega^2 w - c^2\Delta w &=& 0 \quad
  \end{array}\right\} &\;& x\in\Omega,
	\\
  \left.\begin{array}{ccc}
    c (n\cdot\nabla v) - \omega\,w &=& 0 \quad\\
    c (n\cdot\nabla w) + \omega\,v &=& 0 \quad
  \end{array}\right\} &\;& x\in\Gamma_2,
	\\
	\left.\begin{array}{ccc}
    v &=& 1 \quad\\
    w &=& 0 \quad
  \end{array}\right\} &\;& x\in\Gamma_1.


@f}



对于 $\phi,\psi$ 与 $\phi|_{\Gamma_1}=\psi|_{\Gamma_1}=0$ 的测试函数，经过通常的乘法，在 $\Omega$ 上的积分和应用部分积分，我们得到弱的表述

@f{eqnarray*}


-\omega^2 \langle \phi, v \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi, \nabla v \rangle_{\mathrm{L}^2(\Omega)}


- c \omega \langle \phi, w \rangle_{\mathrm{L}^2(\Gamma_2)} &=& 0, \\


-\omega^2 \langle \psi, w \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi, \nabla w \rangle_{\mathrm{L}^2(\Omega)}
+ c \omega \langle \psi, v \rangle_{\mathrm{L}^2(\Gamma_2)} &=& 0.


@f}



我们选择有限元空间 $V_h$ 和 $W_h$ ，基数为 $\{\phi_j\}_{j=1}^n,
\{\psi_j\}_{j=1}^n$ ，寻找近似解

@f[
v_h = \sum_{j=1}^n \alpha_j \phi_j, \;\; w_h = \sum_{j=1}^n \beta_j \psi_j.


@f]

将其插入变异形式中，可以得到方程组

@f[
\renewcommand{\arraystretch}{2.0}
\left.\begin{array}{ccc}
\sum_{j=1}^n
\left(


-\omega^2 \langle \phi_i, \phi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi_i, \nabla \phi_j \rangle_{\mathrm{L}^2(\Omega)}
\right)
\alpha_j


- \left(
c\omega \langle \phi_i,\psi_j\rangle_{\mathrm{L}^2(\Gamma_2)}\right)\beta_j
&=& 0 \\
\sum_{j=1}^n
\left(


-\omega^2 \langle \psi_i, \psi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi_i, \nabla \psi_j \rangle_{\mathrm{L}^2(\Omega)}
\right)\beta_j
+ \left(
c\omega \langle
\psi_i,\phi_j\rangle_{\mathrm{L}^2(\Gamma_2)}
\right)\alpha_j
&=& 0
\end{array}\right\}\;\;\forall\; i =1,\ldots,n.


@f]

用矩阵符号表示。

@f[
\renewcommand{\arraystretch}{2.0}
\left(
\begin{array}{cc}


-\omega^2 \langle \phi_i, \phi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \phi_i, \nabla \phi_j \rangle_{\mathrm{L}^2(\Omega)}
& -c\omega \langle \phi_i,\psi_j\rangle_{\mathrm{L}^2(\Gamma_2)} \\
c\omega \langle \psi_i,\phi_j\rangle_{\mathrm{L}^2(\Gamma_2)}
& -\omega^2 \langle \psi_{i}, \psi_j \rangle_{\mathrm{L}^2(\Omega)}
+ c^2 \langle \nabla \psi_{i}, \nabla \psi_j  \rangle_{\mathrm{L}^2(\Omega)}
\end{array}
\right)
\left(
\begin{array}{c}
\alpha \\ \beta
\end{array}
\right)
=
\left(
\begin{array}{c}
0 \\ 0
\end{array}
\right)


@f]

(不要被这里的右手边为零所迷惑，那是因为我们还没有包括Dirichlet边界数据)。由于非对角线区块的交替符号，我们已经可以看到这个系统是非对称的，事实上它甚至是不确定的。当然，没有必要选择空间 $V_h$ 和 $W_h$ 是相同的。然而，我们期望解的实部和虚部具有类似的性质，因此在实现中确实会采取 $V_h=W_h$ ，并且也会对两个空间使用相同的基函数 $\phi_i = \psi_i$ 。使用不同符号的原因只是让我们能够区分 $v$ 和 $w$ 的形状函数，因为这种区分在实施中起着重要作用。




<a name="Thetestcase"></a><h3>The test case</h3>


在计算中，我们将考虑波在单位方阵中的传播，超声由换能器透镜产生，透镜的形状是圆的一段，中心在 $(0.5, d)$ ，半径略大于 $d$ ；这种形状应该导致声波在圆中心的聚焦。改变 $d$ 会改变透镜的 "焦点"，并影响 $u$ 强度的空间分布，我们主要关注的是 $|u|=\sqrt{v^2+w^2}$ 的聚焦效果如何。

在下面的程序中，我们将使用实部和虚部分裂的公式来实现复值亥姆霍兹方程。我们还将讨论如何生成一个看起来像正方形并带有轻微隆起的模拟换能器的域（在 <code>UltrasoundProblem<dim>::make_grid()</code> 函数中），以及如何生成不仅包含解分量 $v$ 和 $w$ ，而且直接在输出文件中包含幅值 $\sqrt{v^2+w^2}$ 的图形输出（在 <code>UltrasoundProblem<dim>::output_results()</code> ）。最后，我们使用ParameterHandler类来轻松读取参数，如焦距 $d$ 、波速 $c$ 、频率 $\omega$ ，以及在运行时从输入文件中读取其他一些参数，而不是在源代码中固定这些参数，因为每次我们想改变参数时，都必须重新编译。


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
 * 下面的头文件以前都讨论过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * 
 * @endcode
 * 
 * 这个头文件包含了ParameterHandler类的必要声明，我们将用它来从配置文件中读取我们的参数。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/parameter_handler.h> 
 * 
 * @endcode
 * 
 * 为了解决线性系统，我们将使用UMFPACK提供的稀疏LU分解（见SparseDirectUMFPACK类），为此需要以下头文件。 请注意，为了编译这个教程程序，deal.II-library需要在UMFPACK支持下构建，默认情况下是启用的。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * @endcode
 * 
 * FESystem类允许我们将多个FE对象堆叠成一个复合的、矢量值的有限元场。该类的必要声明在该头文件中提供。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_system.h> 
 * 
 * @endcode
 * 
 * 最后，包括声明定时器类的头文件，我们将用它来确定我们程序的每个操作需要多少时间。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/timer.h> 
 * 
 * @endcode
 * 
 * 作为本程序开始时的最后一步，我们将本程序中的所有内容放入其命名空间，并在其中使deal.II命名空间中的所有内容全局可用，不需要在所有内容前加上 <code>dealii</code><code>::</code>  。
 * 

 * 
 * 
 * @code
 * namespace Step29 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeDirichletBoundaryValuescodeclass"></a> 
 * <h3>The <code>DirichletBoundaryValues</code> class</h3>
 * 

 * 
 * 首先我们为代表Dirichlet边界值的函数定义一个类。这在以前已经做过很多次了，因此不需要过多解释。
 * 

 * 
 * 由于有两个值 $v$ 和 $w$ 需要在边界处规定，我们必须告诉基类这是一个有两个分量的向量值函数， <code>vector_value</code> 函数和它的表亲 <code>vector_value_list</code> 必须返回有两个条目的向量。在我们的例子中，这个函数非常简单，它只是对实部 $v$ 返回1，对虚部 $w$ 返回0，而不管它在哪个点被评估。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class DirichletBoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     DirichletBoundaryValues() 
 *       : Function<dim>(2) 
 *     {} 
 * 
 *     virtual void vector_value(const Point<dim> & /*p*/, 
 *                               Vector<double> &values) const override 
 *     { 
 *       Assert(values.size() == 2, ExcDimensionMismatch(values.size(), 2)); 
 * 
 *       values(0) = 1; 
 *       values(1) = 0; 
 *     } 
 * 
 *     virtual void 
 *     vector_value_list(const std::vector<Point<dim>> &points, 
 *                       std::vector<Vector<double>> &  value_list) const override 
 *     { 
 *       Assert(value_list.size() == points.size(), 
 *              ExcDimensionMismatch(value_list.size(), points.size())); 
 * 
 *       for (unsigned int p = 0; p < points.size(); ++p) 
 *         DirichletBoundaryValues<dim>::vector_value(points[p], value_list[p]); 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ThecodeParameterReadercodeclass"></a> 
 * <h3>The <code>ParameterReader</code> class</h3>
 * 

 * 
 * 下一个类负责准备ParameterHandler对象并从输入文件中读取参数。 它包括一个声明所有必要参数的函数 <code>declare_parameters</code> 和一个从外部调用的 <code>read_parameters</code> 函数，以启动参数读取过程。
 * 

 * 
 * 
 * @code
 *   class ParameterReader : public Subscriptor 
 *   { 
 *   public: 
 *     ParameterReader(ParameterHandler &); 
 *     void read_parameters(const std::string &); 
 * 
 *   private: 
 *     void              declare_parameters(); 
 *     ParameterHandler &prm; 
 *   }; 
 * 
 * @endcode
 * 
 * 构造函数存储了一个传递给它的ParameterHandler对象的引用。
 * 

 * 
 * 
 * @code
 *   ParameterReader::ParameterReader(ParameterHandler &paramhandler) 
 *     : prm(paramhandler) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="codeParameterReaderdeclare_parameterscode"></a> 
 * <h4><code>ParameterReader::declare_parameters</code></h4>
 * 

 * 
 * <code>declare_parameters</code> 函数声明了我们的ParameterHandler对象能够从输入文件中读取的所有参数，以及它们的类型、范围条件和它们出现在哪个分段。我们将用一对大括号包住所有进入一个部分的条目，以迫使编辑器将它们缩进一级，使其更容易阅读哪些条目共同构成一个部分。
 * 

 * 
 * 
 * @code
 *   void ParameterReader::declare_parameters() 
 *   { 
 * 
 * @endcode
 * 
 * 网格和几何参数包括应用于初始粗略网格的全局细化步数和换能器镜头的焦距 $d$ 。对于细化步数，我们允许在 $[0,\infty)$ 范围内的整数，其中 Patterns::Integer 对象的第二个参数被省略，表示半开区间。 对于焦距，任何大于零的数字都可以接受。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Mesh & geometry parameters"); 
 *     { 
 *       prm.declare_entry("Number of refinements", 
 *                         "6", 
 *                         Patterns::Integer(0), 
 *                         "Number of global mesh refinement steps " 
 *                         "applied to initial coarse grid"); 
 * 
 *       prm.declare_entry("Focal distance", 
 *                         "0.3", 
 *                         Patterns::Double(0), 
 *                         "Distance of the focal point of the lens " 
 *                         "to the x-axis"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 下一小节专门讨论方程中出现的物理参数，它们是频率  $\omega$  和波速  $c$  。同样，两者都需要位于半开区间 $[0,\infty)$ 内，通过调用 Patterns::Double 类，仅以左端点为参数来表示。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Physical constants"); 
 *     { 
 *       prm.declare_entry("c", "1.5e5", Patterns::Double(0), "Wave speed"); 
 * 
 *       prm.declare_entry("omega", "5.0e7", Patterns::Double(0), "Frequency"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 最后但并非最不重要的是，我们希望能够通过配置文件中的条目来改变输出的一些属性，如文件名和格式，这就是最后一小节的目的。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Output parameters"); 
 *     { 
 *       prm.declare_entry("Output filename", 
 *                         "solution", 
 *                         Patterns::Anything(), 
 *                         "Name of the output file (without extension)"); 
 * 
 * @endcode
 * 
 * 由于不同的输出格式在生成输出时可能需要不同的参数（例如，postscript输出需要视角角度、线宽、颜色等），如果我们必须为库中支持的每一种可能的输出格式手工声明所有这些参数，那就太麻烦了。相反，每种输出格式都有一个 <code>FormatFlags::declare_parameters</code> 函数，它在自己的小节中声明了该格式的所有特定参数。下面调用 DataOutInterface<1>::declare_parameters 为所有可用的输出格式执行 <code>declare_parameters</code> ，这样就为每一种格式创建了一个自己的小节，为该特定的输出格式声明参数。(上面 <code>@<1@></code> 的调用中，模板参数的实际值在这里并不重要：该函数做了同样的工作，与维度无关，但恰好是在一个依赖模板参数的类中。)  要想知道哪种输出格式有哪些参数，你可以查阅DataOutBase类的文档，或者干脆在没有参数文件的情况下运行这个程序。然后它将创建一个文件，其中所有声明的参数都设置为默认值，这可以方便地作为一个起点，将参数设置为你想要的值。
 * 

 * 
 * 
 * @code
 *       DataOutInterface<1>::declare_parameters(prm); 
 *     } 
 *     prm.leave_subsection(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeParameterReaderread_parameterscode"></a> 
 * <h4><code>ParameterReader::read_parameters</code></h4>
 * 

 * 
 * 这是ParameterReader类中的主函数。 它从外部被调用，首先声明所有的参数，然后从输入文件中读取参数，文件名由调用者提供。对这个函数的调用完成后，可以用 <code>prm</code> 对象来检索从文件中读入的参数值。
 * 

 * 
 * 
 * @code
 *   void ParameterReader::read_parameters(const std::string &parameter_file) 
 *   { 
 *     declare_parameters(); 
 * 
 *     prm.parse_input(parameter_file); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeComputeIntensitycodeclass"></a> 
 * <h3>The <code>ComputeIntensity</code> class</h3>
 * 

 * 
 * 正如介绍中所提到的，我们真正追求的量是超声波强度的空间分布，它对应于  $|u|=\sqrt{v^2+w^2}$  。现在我们可以只满足于在输出中拥有 $v$ 和 $w$ ，并使用合适的可视化或后处理工具从我们计算的解决方案中得出 $|u|$ 。然而，也有一种方法可以输出从deal.II中的解决方案中得出的数据，我们在这里要利用这个机制。
 * 

 * 
 * 到目前为止，我们一直使用 DataOut::add_data_vector 函数将包含输出数据的向量添加到一个DataOut对象中。 这个函数有一个特殊的版本，除了数据向量之外，还有一个额外的参数类型为DataPostprocessor。当这个函数用于输出时，在每个要生成输出数据的点上，指定的DataPostprocessor对象的 DataPostprocessor::evaluate_scalar_field() 或 DataPostprocessor::evaluate_vector_field() 函数被调用，从数据向量代表的有限元函数的值、梯度和二阶导数计算输出量（在面相关数据的情况下，法向量也是可用的）。因此，这使我们可以输出任何可以从解的值及其导数中局部导出的数量。 当然，超声强度 $|u|$ 就是这样一个量，它的计算甚至不涉及任何 $v$ 或 $w$ 的导数。
 * 

 * 
 * 在实践中，DataPostprocessor类只提供了这个功能的接口，我们需要从它派生出我们自己的类，以实现接口所指定的功能。在最一般的情况下，我们必须实现几个成员函数，但是如果输出量是一个单一的标量，那么其中的一些模板代码可以由一个更专业的类DataPostprocessorScalar来处理，我们可以从这个类派生。这就是 <code>ComputeIntensity</code> 类的作用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ComputeIntensity : public DataPostprocessorScalar<dim> 
 *   { 
 *   public: 
 *     ComputeIntensity(); 
 * 
 *     virtual void evaluate_vector_field( 
 *       const DataPostprocessorInputs::Vector<dim> &inputs, 
 *       std::vector<Vector<double>> &computed_quantities) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 在构造函数中，我们需要用两个参数调用基类的构造函数。第一个参数表示由该类计算的单一标量在输出文件中应表示的名称。在我们的例子中，后处理程序有 $|u|$ 作为输出，所以我们使用 "Intensity"。
 * 

 * 
 * 第二个参数是一组标志，表示后处理程序需要哪些数据来计算输出量。 这可以是update_values、update_gradients和update_hessians（如果是脸部数据，也可以是update_normal_vector）的任何一个子集，这些都在UpdateFlags中记录。 当然，导数的计算需要额外的资源，所以这里只应该给出真正需要的数据的标志，就像我们使用FEValues对象时一样。在我们的例子中，只有  $v$  和  $w$  的函数值需要用来计算  $|u|$  ，所以我们用 update_values 标志就可以了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ComputeIntensity<dim>::ComputeIntensity() 
 *     : DataPostprocessorScalar<dim>("Intensity", update_values) 
 *   {} 
 * 
 * @endcode
 * 
 * 实际的后处理发生在下面这个函数中。它的输入是一个存储函数值的对象（这里是向量值），代表给 DataOut::add_data_vector, 的数据向量在我们产生输出的所有评估点的评估值，以及一些代表导数的张量对象（我们在这里没有使用，因为 $|u|$ 只是从 $v$ 和 $w$ 计算出来）。派生量在 <code>computed_quantities</code> 向量中返回。请记住，这个函数只能使用由  <code>get_needed_update_flags</code>  指定的各自更新标志的数据。例如，我们可能不会在这里使用导数，因为我们对  <code>get_needed_update_flags</code>  的实现要求只提供函数值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ComputeIntensity<dim>::evaluate_vector_field( 
 *     const DataPostprocessorInputs::Vector<dim> &inputs, 
 *     std::vector<Vector<double>> &               computed_quantities) const 
 *   { 
 *     Assert(computed_quantities.size() == inputs.solution_values.size(), 
 *            ExcDimensionMismatch(computed_quantities.size(), 
 *                                 inputs.solution_values.size())); 
 * 
 * @endcode
 * 
 * 计算本身是很简单的。我们遍历输出向量中的每个条目，并从 $v$ 和 $w$ 的相应值中计算出 $|u|$ 。我们通过创建一个复数 $u$ ，然后对结果调用 `std::abs()` 来实现。(我们可能想调用 `std::norm()`, ，但是在一个历史的怪圈中，C++委员会决定 `std::norm()` 应该返回绝对值的<i>square</i>--从而不满足数学家对所谓 "规范 "的属性要求。)
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < computed_quantities.size(); i++) 
 *       { 
 *         Assert(computed_quantities[i].size() == 1, 
 *                ExcDimensionMismatch(computed_quantities[i].size(), 1)); 
 *         Assert(inputs.solution_values[i].size() == 2, 
 *                ExcDimensionMismatch(inputs.solution_values[i].size(), 2)); 
 * 
 *         const std::complex<double> u(inputs.solution_values[i](0), 
 *                                      inputs.solution_values[i](1)); 
 * 
 *         computed_quantities[i](0) = std::abs(u); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ThecodeUltrasoundProblemcodeclass"></a> 
 * <h3>The <code>UltrasoundProblem</code> class</h3>
 * 

 * 
 * 最后这里是这个程序的主类。 它的成员函数与前面的例子非常相似，特别是 step-4 ，成员变量的列表也没有什么大的惊喜。传递给构造函数的ParameterHandler对象被存储为一个引用，以便于从该类的所有函数中访问参数。 由于我们正在使用矢量值的有限元，我们使用的FE对象是FESystem类型的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class UltrasoundProblem 
 *   { 
 *   public: 
 *     UltrasoundProblem(ParameterHandler &); 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void output_results() const; 
 * 
 *     ParameterHandler &prm; 
 * 
 *     Triangulation<dim> triangulation; 
 *     DoFHandler<dim>    dof_handler; 
 *     FESystem<dim>      fe; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 *     Vector<double>       solution, system_rhs; 
 *   }; 
 * 
 * @endcode
 * 
 * 构造函数接收ParameterHandler对象并将其存储在一个引用中。它还初始化了DoF-Handler和有限元系统，该系统由标量Q1场的两个副本组成，一个用于 $v$ ，一个用于 $w$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   UltrasoundProblem<dim>::UltrasoundProblem(ParameterHandler &param) 
 *     : prm(param) 
 *     , dof_handler(triangulation) 
 *     , fe(FE_Q<dim>(1), 2) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemmake_gridcode"></a> 
 * <h4><code>UltrasoundProblem::make_grid</code></h4>
 * 

 * 
 * 这里我们为我们的领域设置网格。 正如论述中所提到的，这个几何体只是一个单位正方形（2d），其边界部分代表换能器透镜，由一个圆的扇形代替。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::make_grid() 
 *   { 
 * 
 * @endcode
 * 
 * 首先我们生成一些日志输出，并启动一个定时器，这样我们就可以在这个函数完成后计算出执行时间。
 * 

 * 
 * 
 * @code
 *     std::cout << "Generating grid... "; 
 *     Timer timer; 
 * 
 * @endcode
 * 
 * 然后，我们从ParameterHandler对象中查询换能器镜头的焦距和网格细化步数的值。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Mesh & geometry parameters"); 
 * 
 *     const double       focal_distance = prm.get_double("Focal distance"); 
 *     const unsigned int n_refinements = prm.get_integer("Number of refinements"); 
 * 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 接下来，为换能器镜头的位置和焦点定义了两个点，也就是圆的中心，其线段将形成边界的换能器部分。注意，这是程序中唯一一个在二维和三维中略有不同的地方。尽管本教程只涉及二维情况，但要使这个程序在三维中发挥作用，必要的补充是非常少的，所以我们选择包括它们。
 * 

 * 
 * 
 * @code
 *     const Point<dim> transducer = 
 *       (dim == 2) ? Point<dim>(0.5, 0.0) : Point<dim>(0.5, 0.5, 0.0); 
 *     const Point<dim> focal_point = (dim == 2) ? 
 *                                      Point<dim>(0.5, focal_distance) : 
 *                                      Point<dim>(0.5, 0.5, focal_distance); 
 * 
 * @endcode
 * 
 * 作为初始粗网格，我们采用一个简单的单位正方形，每个方向上有5个细分。分区的数量是这样选择的：我们想指定为传感器边界的线段 $[0.4,0.6]$ 是由一个面来跨越的。然后，我们通过所有的单元格，找到换能器所在的面，事实上，这只是X轴上从0.4到0.6的一条边。这是我们希望根据圆环形边界进行细化的地方，所以我们用不同的流形指标来标记这个边缘。由于我们要在换能器上设置迪里希特边界条件，所以我们也要改变其边界指标。
 * 

 * 
 * 
 * @code
 *     GridGenerator::subdivided_hyper_cube(triangulation, 5, 0, 1); 
 * 
 *     for (auto &cell : triangulation.cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->at_boundary() && 
 *             ((face->center() - transducer).norm_square() < 0.01)) 
 *           { 
 *             face->set_boundary_id(1); 
 *             face->set_manifold_id(1); 
 *           } 
 * 
 * @endcode
 * 
 * 对于换能器镜头的圆形部分，使用了一个SphericalManifold对象（当然，在2D中只是代表一个圆），中心的计算方法如上。
 * 

 * 
 * 
 * @code
 *     triangulation.set_manifold(1, SphericalManifold<dim>(focal_point)); 
 * 
 * @endcode
 * 
 * 现在，全局细化被执行。靠近换能器位置的单元格将根据换能器透镜的圆形边界被自动细化。
 * 

 * 
 * 
 * @code
 *     triangulation.refine_global(n_refinements); 
 * 
 * @endcode
 * 
 * 最后，我们再生成一些日志输出。我们停止定时器，并查询自函数开始以来所经过的CPU秒数。
 * 

 * 
 * 
 * @code
 *     timer.stop(); 
 *     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
 * 
 *     std::cout << "  Number of active cells:  " << triangulation.n_active_cells() 
 *               << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemsetup_systemcode"></a> 
 * <h4><code>UltrasoundProblem::setup_system</code></h4>
 * 

 * 
 * 系统矩阵的初始化、稀疏模式和向量与前面的例子相同，因此不需要进一步评论。和前面的函数一样，我们也输出我们在这里所做的运行时间。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::setup_system() 
 *   { 
 *     std::cout << "Setting up system... "; 
 *     Timer timer; 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *     solution.reinit(dof_handler.n_dofs()); 
 * 
 *     timer.stop(); 
 *     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
 * 
 *     std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemassemble_systemcode"></a> 
 * <h4><code>UltrasoundProblem::assemble_system</code></h4>
 * 

 * 
 * 和以前一样，这个函数负责组装系统矩阵和右手边的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::assemble_system() 
 *   { 
 *     std::cout << "Assembling system matrix... "; 
 *     Timer timer; 
 * 
 * @endcode
 * 
 * 首先我们从ParameterHandler对象中查询波速和频率，并将其存储在本地变量中，因为它们将在本函数中频繁使用。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Physical constants"); 
 * 
 *     const double omega = prm.get_double("omega"), c = prm.get_double("c"); 
 * 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 像往常一样，计算积分时使用普通的高斯正交规则。由于我们的双线性形式涉及到 $\Gamma_2$ 上的边界积分，所以我们还需要一个正交法则来计算面的积分，这些面是 $dim-1$ 维的。
 * 

 * 
 * 
 * @code
 *     QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *     QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size(), 
 *                        n_face_q_points = face_quadrature_formula.size(), 
 *                        dofs_per_cell   = fe.n_dofs_per_cell(); 
 * 
 * @endcode
 * 
 * FEValues对象将为我们评估形状函数。 对于涉及到 $\Omega$ 上的积分的双线性形式的部分，我们需要形状函数的值和梯度，当然还有正交权重。 对于涉及边界积分的条款，只需要形状函数值和正交权重。
 * 

 * 
 * 
 * @code
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_face_values(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_JxW_values); 
 * 
 * @endcode
 * 
 * 像往常一样，系统矩阵是逐个单元组装的，我们需要一个矩阵来存储本地单元的贡献，以及一个索引向量来将单元的贡献转移到全局系统矩阵中的适当位置，然后。
 * 

 * 
 * 
 * @code
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 * 
 * @endcode
 * 
 * 在每个单元，我们首先需要重置本地贡献矩阵，并请求FEValues对象计算当前单元的形状函数。
 * 

 * 
 * 
 * @code
 *         cell_matrix = 0; 
 *         fe_values.reinit(cell); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           { 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               { 
 * 
 * @endcode
 * 
 * 在这一点上，重要的是要记住，我们所处理的是一个有两个组成部分的有限元系统。由于我们构造这个FESystem的方式，即作为两个标量有限元场的笛卡尔乘积，每个形状函数只有一个非零分量（用deal.II的行话说，它们是 @ref  GlossPrimitive "原始"）。 因此，每个形状函数可以被看作是引言中的 $\phi$ 's或 $\psi$ 's之一，同样，相应的自由度也可以归属于 $\alpha$ 或 $\beta$  。      然而，当我们遍历当前单元上的所有自由度时，它们并不以任何特定的顺序出现，因此我们无法立即决定索引为 $i$ 和 $j$ 的自由度是属于我们解决方案的实部还是虚部。 另一方面，如果你看一下介绍中的系统矩阵的形式，这个区别是至关重要的，因为它将决定当前一对DoF的贡献将归入系统矩阵的哪个块，因此我们需要从给定的两个形状函数中计算哪个数量。 幸运的是，FESystem对象可以为我们提供这些信息，即它有一个函数 FESystem::system_to_component_index, ，为每个局部的DoF索引返回一对整数，其中第一个表示该DoF属于系统的哪个组成部分。这对整数中的第二个整数表示该DoF在标量基有限元场中的索引，但这一信息在这里并不相关。如果你想知道更多关于这个函数和原始向量值元素背后的基本方案，可以看看 step-8 或 @ref vector_valued 模块，那里对这些主题有深入的解释。
 * 

 * 
 * 
 * @code
 *                 if (fe.system_to_component_index(i).first == 
 *                     fe.system_to_component_index(j).first) 
 *                   { 
 * 
 * @endcode
 * 
 * 如果DoF $i$ 和 $j$ 都属于同一个分量，即它们的形状函数都是 $\phi$ 的，或者都是 $\psi$ 的，贡献将最终出现在我们系统矩阵的一个对角块中，由于相应的条目是由同一个公式计算的，我们不必理会它们实际上是 $\phi$ 还是 $\psi$ 形状函数。我们可以简单地通过遍历所有正交点并将其贡献相加来计算条目，其中形状函数的值和梯度由我们的FEValues对象提供。
 * 

 * 
 * 
 * @code
 *                     for (unsigned int q_point = 0; q_point < n_q_points; 
 *                          ++q_point) 
 *                       cell_matrix(i, j) += 
 *                         (((fe_values.shape_value(i, q_point) * 
 *                            fe_values.shape_value(j, q_point)) * 
 *                             (-omega * omega) + 
 *                           (fe_values.shape_grad(i, q_point) * 
 *                            fe_values.shape_grad(j, q_point)) * 
 *                             c * c) * 
 *                          fe_values.JxW(q_point)); 
 * 
 * @endcode
 * 
 * 你可能认为我们在向FEValues对象请求形状函数值或梯度时，必须指定我们想评估的形状函数的哪个分量。然而，由于形状函数是原始的，它们只有一个非零分量，而且FEValues类足够聪明，它知道我们肯定对这一个非零分量感兴趣。
 * 

 * 
 * 
 * @code
 *                   } 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 我们还必须增加由于边界项的贡献。为此，我们对当前单元格的所有面进行循环，首先看它是否在边界上，其次看它是否有与 $\Gamma_2$ 相关的正确的边界指标，即我们有吸收边界条件的那部分边界。
 * 

 * 
 * 
 * @code
 *         for (const auto face_no : cell->face_indices()) 
 *           if (cell->face(face_no)->at_boundary() && 
 *               (cell->face(face_no)->boundary_id() == 0)) 
 *             { 
 * 
 * @endcode
 * 
 * 这些面肯定会对系统矩阵的非对角线块作出贡献，所以我们要求FEFaceValues对象为我们提供这个面上的形状函数值。
 * 

 * 
 * 
 * @code
 *               fe_face_values.reinit(cell, face_no); 
 * 
 * @endcode
 * 
 * 接下来，我们循环浏览当前单元格的所有DoF，找到属于不同组件且都支持当前面号的一对。
 * 

 * 
 * 
 * @code
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                   if ((fe.system_to_component_index(i).first != 
 *                        fe.system_to_component_index(j).first) && 
 *                       fe.has_support_on_face(i, face_no) && 
 *                       fe.has_support_on_face(j, face_no)) 
 * 
 * @endcode
 * 
 * 检查形状函数在一个面上是否有支持并不是严格必要的：如果我们不检查它，我们会简单地把本地单元矩阵的项加起来，这些项碰巧是零，因为至少有一个形状函数碰巧是零。然而，我们可以通过添加上面的检查来节省这项工作。
 * 

 * 
 * 在任何一种情况下，这些DoFs都会对系统矩阵的对角线区块的边界积分作出贡献。为了计算积分，我们在面的所有正交点上进行循环，然后用面的正交规则提供的正交权重对贡献进行加权求和。 与对角线块上的条目不同，这里的形状函数哪一个是 $\psi$ ，哪一个是 $\phi$ ，确实很重要，因为这将决定该条目的符号。 我们通过一个简单的条件语句来说明这一点，它决定了正确的符号。由于我们已经检查了DoF  $i$  和  $j$  属于不同的组件，这里只需测试其中一个组件属于哪个组件即可。
 * 

 * 
 * 
 * @code
 *                     for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                          ++q_point) 
 *                       cell_matrix(i, j) += 
 *                         ((fe.system_to_component_index(i).first == 0) ? -1 : 
 *                                                                         1) * 
 *                         fe_face_values.shape_value(i, q_point) * 
 *                         fe_face_values.shape_value(j, q_point) * c * omega * 
 *                         fe_face_values.JxW(q_point); 
 *             } 
 * 
 * @endcode
 * 
 * 现在，我们已经完成了这个单元，必须将其贡献从本地系统矩阵转移到全局系统矩阵。为此，我们首先得到这个单元的全局指数列表......
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 * @endcode
 * 
 * ...然后将这些条目逐一添加到系统矩阵中。
 * 

 * 
 * 
 * @code
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             system_matrix.add(local_dof_indices[i], 
 *                               local_dof_indices[j], 
 *                               cell_matrix(i, j)); 
 *       } 
 * 
 * @endcode
 * 
 * 唯一剩下的是 $\Gamma_1$ 上的迪里希特边界值，其特征是边界指标1。Dirichlet值由我们上面定义的 <code>DirichletBoundaryValues</code> 类提供。
 * 

 * 
 * 
 * @code
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              1, 
 *                                              DirichletBoundaryValues<dim>(), 
 *                                              boundary_values); 
 * 
 *     MatrixTools::apply_boundary_values(boundary_values, 
 *                                        system_matrix, 
 *                                        solution, 
 *                                        system_rhs); 
 * 
 *     timer.stop(); 
 *     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemsolvecode"></a> 
 * <h4><code>UltrasoundProblem::solve</code></h4>
 * 

 * 
 * 正如介绍中已经提到的，系统矩阵既不是对称的，也不是确定的，因此，如何提出一个迭代求解器和预处理器来很好地处理这个矩阵并不是很明显。 我们选择了另一种方式，用UMFPACK提供的稀疏LU分解来解决线性系统。这通常是二维问题的一个很好的首选，即使对于大量的DoF也能很好地工作。 SparseDirectUMFPACK类提供了UMFPACK的deal.II接口，它非常容易使用，使我们仅用3行代码就能解决我们的线性系统。
 * 

 * 
 * 再次注意，为了编译这个例子程序，你需要有支持UMFPACK的deal.II库。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::solve() 
 *   { 
 *     std::cout << "Solving linear system... "; 
 *     Timer timer; 
 * 
 * @endcode
 * 
 * 解决线性系统的代码很短：首先，我们分配一个正确类型的对象。下面的 <code>initialize</code> 调用提供了我们想要反转的矩阵给SparseDirectUMFPACK对象，并同时启动了LU分解。因此，这也是这个程序中大部分计算工作发生的地方。
 * 

 * 
 * 
 * @code
 *     SparseDirectUMFPACK A_direct; 
 *     A_direct.initialize(system_matrix); 
 * 
 * @endcode
 * 
 * 分解之后，我们可以把 <code>A_direct</code> 当作代表我们系统矩阵的逆矩阵来使用，所以要计算出解决方案，我们只需要与右边的向量相乘。
 * 

 * 
 * 
 * @code
 *     A_direct.vmult(solution, system_rhs); 
 * 
 *     timer.stop(); 
 *     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemoutput_resultscode"></a> 
 * <h4><code>UltrasoundProblem::output_results</code></h4>
 * 

 * 
 * 这里我们按照参数文件中指定的格式输出我们的解 $v$ 和 $w$ 以及导出的量 $|u|$ 。从 $v$ 和 $w$ 导出 $|u|$ 的大部分工作已经在 <code>ComputeIntensity</code> 类的实现中完成，因此输出程序相当简单，与前面教程中的内容非常相似。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::output_results() const 
 *   { 
 *     std::cout << "Generating output... "; 
 *     Timer timer; 
 * 
 * @endcode
 * 
 * 定义我们的  <code>ComputeIntensity</code>  类的对象和一个 DataOut 对象。
 * 

 * 
 * 
 * @code
 *     ComputeIntensity<dim> intensities; 
 *     DataOut<dim>          data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 * @endcode
 * 
 * 接下来我们从ParameterHandler查询输出相关的参数。 DataOut::parse_parameters  调用作为  <code>ParameterReader::declare_parameters</code>  中  DataOutInterface<1>::declare_parameters  调用的对应方。它从ParameterHandler收集所有与输出格式相关的参数，并相应地设置DataOut对象的相应属性。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Output parameters"); 
 * 
 *     const std::string output_filename = prm.get("Output filename"); 
 *     data_out.parse_parameters(prm); 
 * 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 现在，我们从ParameterHandler提供的基本名称和DataOut类提供的后缀（默认后缀被设置为正确的类型，与.prm文件中通过parse_parameters()设置的后缀相匹配）来组合文件名。
 * 

 * 
 * 
 * @code
 *     const std::string filename = output_filename + data_out.default_suffix(); 
 * 
 *     std::ofstream output(filename); 
 * 
 * @endcode
 * 
 * 解向量 $v$ 和 $w$ 以常规方式添加到DataOut对象中。
 * 

 * 
 * 
 * @code
 *     std::vector<std::string> solution_names; 
 *     solution_names.emplace_back("Re_u"); 
 *     solution_names.emplace_back("Im_u"); 
 * 
 *     data_out.add_data_vector(solution, solution_names); 
 * 
 * @endcode
 * 
 * 对于强度，我们只是再次调用 <code>add_data_vector</code> ，但这次用我们的 <code>ComputeIntensity</code> 对象作为第二个参数，这实际上是将 $|u|$ 加入到输出数据中。
 * 

 * 
 * 
 * @code
 *     data_out.add_data_vector(solution, intensities); 
 * 
 * @endcode
 * 
 * 最后的步骤和以前一样。注意，现在实际的输出格式是由输入文件中的内容决定的，也就是说，人们可以改变输出格式而不必重新编译这个程序。
 * 

 * 
 * 
 * @code
 *     data_out.build_patches(); 
 *     data_out.write(output); 
 * 
 *     timer.stop(); 
 *     std::cout << "done (" << timer.cpu_time() << "s)" << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeUltrasoundProblemruncode"></a> 
 * <h4><code>UltrasoundProblem::run</code></h4>
 * 

 * 
 * 这里我们简单地一个接一个地执行我们的函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void UltrasoundProblem<dim>::run() 
 *   { 
 *     make_grid(); 
 *     setup_system(); 
 *     assemble_system(); 
 *     solve(); 
 *     output_results(); 
 *   } 
 * } // namespace Step29 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h4>The <code>main</code> function</h4>
 * 

 * 
 * 最后是该程序的 <code>main</code> 功能。它的结构与其他几乎所有的教程程序相同。唯一的例外是，我们定义了ParameterHandler和 <code>ParameterReader</code> 对象，并让后者从一个名为 <code>step-29.prm</code> 的文本文件中读入参数值。这样读到的值会被交给UltrasoundProblem类的一个实例。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step29; 
 * 
 *       ParameterHandler prm; 
 *       ParameterReader  param(prm); 
 *       param.read_parameters("step-29.prm"); 
 * 
 *       UltrasoundProblem<2> ultrasound_problem(prm); 
 *       ultrasound_problem.run(); 
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
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-29/doc/results.dox

<a name="Results"></a>

<a name="Results"></a><h1>Results</h1>


当前程序从一个名为 <code>\step-29.prm</code> 的输入文件中读取其运行时参数，该文件看起来像这样。

@code
subsection Mesh & geometry parameters
  # Distance of the focal point of the lens to the x-axis
  set Focal distance        = 0.3


  # Number of global mesh refinement steps applied to initial coarse grid
  set Number of refinements = 5
end



subsection Physical constants
  # Wave speed
  set c     = 1.5e5


  # Frequency
  set omega = 3.0e7
end



subsection Output parameters
  # Name of the output file (without extension)
  set Output file   = solution


  # A name for the output format to be used
  set Output format = vtu
end
@endcode



可以看出，我们设置了 $d=0.3$  ，相当于换能器镜头的焦点在 $x=0.5$  ， $y=0.3$  。粗略的网格被细化了5次，结果是160x160个单元，输出结果以vtu格式写入。参数读取器可以理解更多的参数，特别是与输出的生成有关的参数，但是我们在这里不需要这些参数，因此坚持使用其默认值。

这是调试模式下程序的控制台输出。

@code
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Debug configuration
Generating grid... done (0.820449s)
  Number of active cells:  25600
Setting up system... done (1.18392s)
  Number of degrees of freedom: 51842
Assembling system matrix... done (2.33291s)
Solving linear system... done (1.34837s)
Generating output... done (2.05782s)
[100%] Built target run
@endcode



(当然，如果你在本地运行该程序，执行时间会有所不同。)事实上，大部分时间花在组装系统矩阵和生成输出上是由于在调试模式下需要检查许多断言。在发布模式下，程序的这些部分运行得更快，而求解线性系统的速度几乎没有加快。

@code
> make run
[ 66%] Built target step-29
Scanning dependencies of target run
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0144960s)
DEAL::  Number of active cells:  25600
DEAL::Setting up system... done (0.0356880s)
DEAL::  Number of degrees of freedom: 51842
DEAL::Assembling system matrix... done (0.0436570s)
DEAL::Solving linear system... done (1.54733s)
DEAL::Generating output... done (0.720528s)
[100%] Built target run
@endcode



程序的图形输出看起来如下。


 <table align="center" class="doxtable">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.v.png" alt="v = Re(u)">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.w.png" alt="w = Im(u)">
    </td>
  </tr>
  <tr>
    <td colspan="2">
      <img src="https://www.dealii.org/images/steps/developer/step-29.intensity.png" alt="|u|">
    </td>
  </tr>
</table> 

前两张图片显示了 $u$ 的实部和虚部，而最后一张显示了强度 $|u|$ 。我们可以清楚地看到，强度集中在镜头的焦点周围（0.5，0.3），焦点在 $x$ -方向上相当尖锐，但在 $y$ -方向上更加模糊，这是聚焦镜头的几何形状、其有限孔径和问题的波性的结果。

因为五颜六色的图形总是很有趣，而且为了进一步强调聚焦效果，这里还有一组图片，强调了强度在 $x$ -方向上的实际聚焦效果。

 <table align="center" class="doxtable">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.surface.png" alt="|u|">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-29.contours.png" alt="|u|">
    </td>
  </tr>
</table> 


最后，程序的结构使我们很容易确定程序的哪些部分可以随着网格的细化而很好地扩展，哪些部分不可以。下面是5、6、7次全局细化的运行时间。

@code
> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0135260s)
DEAL::  Number of active cells:  25600
DEAL::Setting up system... done (0.0213910s)
DEAL::  Number of degrees of freedom: 51842
DEAL::Assembling system matrix... done (0.0414300s)
DEAL::Solving linear system... done (1.56621s)
DEAL::Generating output... done (0.729605s)
[100%] Built target run


> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.0668490s)
DEAL::  Number of active cells:  102400
DEAL::Setting up system... done (0.109694s)
DEAL::  Number of degrees of freedom: 206082
DEAL::Assembling system matrix... done (0.160784s)
DEAL::Solving linear system... done (7.86577s)
DEAL::Generating output... done (2.89320s)
[100%] Built target run


> make run
[ 66%] Built target step-29
[100%] Run step-29 with Release configuration
DEAL::Generating grid... done (0.293154s)
DEAL::  Number of active cells:  409600
DEAL::Setting up system... done (0.491301s)
DEAL::  Number of degrees of freedom: 821762
DEAL::Assembling system matrix... done (0.605386s)
DEAL::Solving linear system... done (45.1989s)
DEAL::Generating output... done (11.2292s)
[100%] Built target run
@endcode



每次我们细化一次网格，所以每一步的单元和自由度的数量大约是四倍。可以看出，生成网格、设置自由度、组装线性系统和生成输出的规模相当接近于线性，而求解线性系统的操作，自由度的数量每增加4倍，就需要8倍的时间，也就是说，它是 ${\cal O}(N^{3/2})$  。这可以解释为（使用最优排序）有限元矩阵的带宽是  $B={\cal O}(N^{(dim-1)/dim})$  ，而使用LU分解解决带状线性系统的努力是  ${\cal O}(BN)$  。这也解释了为什么该程序也能在3D中运行（在改变了 <code>UltrasoundProblem</code> 对象的维度后），但其扩展性很差，需要非常耐心才能完成对具有明显分辨率的网格上的线性系统的求解，尽管该程序的其他部分的扩展性非常好。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序的一个明显的可能的扩展是在3D中运行它&mdash；毕竟，我们周围的世界是三维的，而超声束在三维介质中传播。你可以通过简单地改变 <code>main()</code> 中主类的模板参数并运行它来尝试。但这不会让你走得很远：当然，如果你按照参数文件中的设置做5个全局细化步骤，就更不会了。你的内存会耗尽，因为网格（含 $(2^5)^3 \cdot 5^3=2^{15}\cdot 125 \approx 4\cdot 10^6$ 单元），特别是稀疏直接求解器会占用太多的内存。然而，如果你有时间的话，你可以用3个全局细化步骤来求解：在2011年初，直接求解大约需要半个小时。然而，你会注意到，这个解是完全错误的：网格大小根本不够小，不能准确地解决解的波浪，你可以在解的图中看到这一点。因此，在这种情况下，如果你不想在这个问题上扔一个更大的（估计是%并行的）机器，那么自适应性是必不可少的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-29.cc"
*/
