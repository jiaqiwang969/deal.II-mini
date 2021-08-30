/**
@page step_33 The step-33 tutorial program
This tutorial depends on step-12, step-71.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Eulerflow">Euler flow</a>
        <li><a href="#Discretization">Discretization</a>
        <li><a href="#Automaticdifferentiation"> Automatic differentiation </a>
        <li><a href="#Trilinossolvers"> Trilinos solvers </a>
        <li><a href="#Adaptivity"> Adaptivity </a>
        <li><a href="#Inputdeckinitialandboundaryconditions">Input deck, initial and boundary conditions</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Eulerequationspecifics">Euler equation specifics</a>
      <ul>
        <li><a href="#Componentdescription">Component description</a>
        <li><a href="#Transformationsbetweenvariables">Transformations between variables</a>
        <li><a href="#EulerEquationscompute_flux_matrix">EulerEquations::compute_flux_matrix</a>
        <li><a href="#EulerEquationscompute_normal_flux">EulerEquations::compute_normal_flux</a>
        <li><a href="#EulerEquationscompute_forcing_vector">EulerEquations::compute_forcing_vector</a>
        <li><a href="#Dealingwithboundaryconditions">Dealing with boundary conditions</a>
        <li><a href="#EulerEquationscompute_refinement_indicators">EulerEquations::compute_refinement_indicators</a>
        <li><a href="#EulerEquationsPostprocessor">EulerEquations::Postprocessor</a>
      </ul>
        <li><a href="#Runtimeparameterhandling">Run time parameter handling</a>
      <ul>
        <li><a href="#ParametersSolver">Parameters::Solver</a>
        <li><a href="#ParametersRefinement">Parameters::Refinement</a>
        <li><a href="#ParametersFlux">Parameters::Flux</a>
        <li><a href="#ParametersOutput">Parameters::Output</a>
        <li><a href="#ParametersAllParameters">Parameters::AllParameters</a>
      </ul>
        <li><a href="#Conservationlawclass">Conservation law class</a>
      <ul>
        <li><a href="#ConservationLawConservationLaw">ConservationLaw::ConservationLaw</a>
        <li><a href="#ConservationLawsetup_system">ConservationLaw::setup_system</a>
        <li><a href="#ConservationLawassemble_system">ConservationLaw::assemble_system</a>
        <li><a href="#ConservationLawassemble_cell_term">ConservationLaw::assemble_cell_term</a>
        <li><a href="#ConservationLawassemble_face_term">ConservationLaw::assemble_face_term</a>
        <li><a href="#ConservationLawsolve">ConservationLaw::solve</a>
        <li><a href="#ConservationLawcompute_refinement_indicators">ConservationLaw::compute_refinement_indicators</a>
        <li><a href="#ConservationLawrefine_grid">ConservationLaw::refine_grid</a>
        <li><a href="#ConservationLawoutput_results">ConservationLaw::output_results</a>
        <li><a href="#ConservationLawrun">ConservationLaw::run</a>
      </ul>
        <li><a href="#main">main()</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Stabilization">Stabilization</a>
        <li><a href="#Betterlinearsolvers">Better linear solvers</a>
        <li><a href="#Cachetheexplicitpartofresidual">Cache the explicit part of residual</a>
        <li><a href="#Otherconservationlaws">Other conservation laws</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-33/doc/intro.dox

 <br> 

<i>
This program was written for fun by David Neckels (NCAR) while working
at Sandia (on the Wyoming Express bus to and from Corrales each day).
The main purpose was to better understand Euler flow.
The code solves the basic Euler equations of gas dynamics, by using a
fully implicit Newton iteration (inspired by Sandia's Aria code).  The
code may be configured by an input file to run different simulations
on different meshes, with differing boundary conditions.
<br>
The original code and documentation was later slightly modified by Wolfgang
Bangerth to make it more modular and allow replacing the parts that are
specific to the Euler equations by other hyperbolic conservation laws without
too much trouble.
</i>

 @note  程序使用<a
href="http://trilinos.org">Trilinos</a>线性求解器（这些可以在Trilinos的Aztec/Amesos包中找到）和一个自动微分包，Sacado，也是Trilinos的一部分。deal.II必须被配置为使用Trilinos。请参考<a
href="../../readme.html#trilinos">ReadMe</a>文件以了解如何做到这一点。

 @note  虽然这个程序很好地展示了自动微分的使用，但它并没有表达欧拉方程求解器的技术水平。对于这个方程有更快、更准确的方法，你应该看看步骤67和步骤69，看看这个方程如何更有效地得到解决。




<a name="Introduction"></a><a name="Intro"></a><h1>Introduction</h1> 。


<a name="Eulerflow"></a><h3>Euler flow</h3>


描述可压缩、无粘性气体运动的方程（所谓的气体动力学欧拉方程）是一个基本的守恒定律系统。在空间维度 $d$ 中，其内容为

@f[
\partial_t \mathbf{w} + \nabla \cdot \mathbf{F}(\mathbf{w}) =
\mathbf{G}(\mathbf w),


@f]

解 $\mathbf{w}=(\rho v_1,\ldots,\rho v_d,\rho,
E)^{\top}$ 包括 $\rho$ 流体密度， ${\mathbf v}=(v_1,\ldots v_d)^T$ 流速（因此 $\rho\mathbf v$ 是线性动量密度），和 $E$ 气体的能量密度。我们将上述方程式解释为  $\partial_t \mathbf{w}_i + \nabla \cdot \mathbf{F}_i(\mathbf{w}) = \mathbf
G_i(\mathbf w)$  ,  $i=1,\ldots,dim+2$  。

对于欧拉方程，通量矩阵 $\mathbf F$ （或通量函数系统）被定义为（这里显示的是情况 $d=3$ ）。

@f{eqnarray*}
  \mathbf F(\mathbf w)
  =
  \left(
  \begin{array}{ccc}
    \rho v_1^2+p & \rho v_2v_1  & \rho v_3v_1 \\
    \rho v_1v_2  & \rho v_2^2+p & \rho v_3v_2 \\
    \rho v_1v_3  & \rho v_2v_3  & \rho v_3^2+p \\
    \rho v_1 & \rho v_2 & \rho v_3 \\
    (E+p) v_1 & (E+p) v_2 & (E+p) v_3
  \end{array}
  \right),


@f}

我们将只选择重力的影响作为特定的右手边强制力，用以下方式描述

@f{eqnarray*}
  \mathbf G(\mathbf w)
  =
  \left(
  \begin{array}{c}
    g_1\rho \\
    g_2\rho \\
    g_3\rho \\
    0 \\
    \rho \mathbf g \cdot \mathbf v
  \end{array}
  \right),


@f}

其中 $\mathbf g=(g_1,g_2,g_3)^T$ 表示重力矢量。有了这个，整个方程组就变成了：

@f{eqnarray*}
  \partial_t (\rho v_i) + \sum_{s=1}^d \frac{\partial(\rho v_i v_s +
  \delta_{is} p)}{\partial x_s} &=& g_i \rho, \qquad i=1,\dots,d, \\
  \partial_t \rho + \sum_{s=1}^d \frac{\partial(\rho v_s)}{\partial x_s} &=& 0,  \\
  \partial_t E + \sum_{s=1}^d \frac{\partial((E+p)v_s)}{\partial x_s} &=&
  \rho \mathbf g \cdot \mathbf v.


@f}

这些方程分别描述了动量、质量和能量的守恒。该系统被一个定义压力的关系所封闭。   $p =
(\gamma -1)(E-\frac{1}{2} \rho |\mathbf v|^2)$  .对于空气（主要是氮气和氧气）和其他双原子气体的成分，其比热比为  $\gamma=1.4$  。

这个问题显然属于矢量值问题的范畴。关于如何在deal.II中处理这些问题的一般概述可以在 @ref vector_valued 模块中找到。

<a name="Discretization"></a><h3>Discretization</h3>


考虑到这是一个双曲问题，与步骤12中讨论的简单问题的风格相同，以通常的方式进行微调：我们选择一个有限元空间 $V_h$ ，并针对我们的（矢量值）测试函数 $\mathbf{z} \in V_h$ 积分我们的守恒法。  然后我们通过部分积分，用<i> numerical </i>通量 $\mathbf{H}$ 来近似边界通量。

@f{eqnarray*}
&&\int_{\Omega} (\partial_t \mathbf{w}, \mathbf{z}) + (\nabla \cdot \mathbf{F}(\mathbf{w}), \mathbf{z}) \\
&\approx &\int_{\Omega} (\partial_t \mathbf{w}, \mathbf{z}) - (\mathbf{F}(\mathbf{w}), \nabla \mathbf{z}) + h^{\eta}(\nabla \mathbf{w} , \nabla \mathbf{z}) + \int_{\partial \Omega} (\mathbf{H}(\mathbf{w}^+, \mathbf{w}^-, \mathbf{n}), \mathbf{z}^+),


@f}

其中上标 $+$ 表示一个函数的内部轨迹， $-$ 表示外部轨迹。扩散项 $h^{\eta}(\nabla \mathbf{w} , \nabla \mathbf{z})$ 是严格为了稳定而引入的，其中 $h$ 是网格大小， $\eta$ 是一个参数，规定要增加多少扩散。

在边界上，我们必须说清楚外痕 $\mathbf{w}^-$ 是什么。根据边界条件，我们规定以下两种情况。   <ul>   <li>  流入边界： $\mathbf{w}^-$ 被规定为理想值。   <li>  超音速流出边界： $\mathbf{w}^- = \mathbf{w}^+$   <li>  亚音速流出边界： $\mathbf{w}^- = \mathbf{w}^+$  除了能量变量被修改为支持规定的压力 $p_o$  ，即 $\mathbf{w}^- =(\rho^+, \rho v_1^+, \dots, \rho v_d^+, p_o/(\gamma -1) + 0.5 \rho |\mathbf{v}^+|^2)$   <li>  反射边界：我们设定 $\mathbf{w}^-$  ，使 $(\mathbf{v}^+ + \mathbf{v}^-) \cdot \mathbf{n} = 0$  和 $\rho^- = \rho^+,E^-=E^+$  。   </ul> 

关于这些问题的更多信息可以在Ralf Hartmann的博士论文中找到（"Adaptive Finite Element Methods for the Compressible Euler Equations"，博士论文，海德堡大学，2002）。

我们使用时间步长方案来替代上述方程中的时间导数。为了简单起见，我们将 $ \mathbf{B}({\mathbf{w}_{n}})(\mathbf z) $ 定义为时间步长 $n$ 的空间残差。

@f{eqnarray*}
 \mathbf{B}(\mathbf{w}_{n})(\mathbf z)  &=&


- \int_{\Omega} \left(\mathbf{F}(\mathbf{w}_n),
\nabla\mathbf{z}\right) +  h^{\eta}(\nabla \mathbf{w}_n , \nabla \mathbf{z}) \\
&& +
\int_{\partial \Omega} \left(\mathbf{H}(\mathbf{w}_n^+,
\mathbf{w}^-(\mathbf{w}_n^+), \mathbf{n}), \mathbf{z}\right)


-
\int_{\Omega} \left(\mathbf{G}(\mathbf{w}_n),
\mathbf{z}\right) .


@f}



因此，在每个时间步骤，我们的完全离散化是，应用于任何测试函数 $\mathbf z$ 的残差等于零。

@f{eqnarray*}
R(\mathbf{W}_{n+1})(\mathbf z) &=&
\int_{\Omega} \left(\frac{{\mathbf w}_{n+1} - \mathbf{w}_n}{\delta t},
\mathbf{z}\right)+
\theta \mathbf{B}({\mathbf{w}}_{n+1}) +  (1-\theta) \mathbf{B}({\mathbf w}_{n}) \\
&=& 0


@f}

其中 $ \theta \in [0,1] $ 和 $\mathbf{w}_i = \sum_k \mathbf{W}_i^k \mathbf{\phi}_k$  。选择 $\theta=0$ 的结果是显式（正向）欧拉方案， $\theta=1$ 是稳定的隐式（反向）欧拉方案，而 $\theta=\frac 12$ 是克拉克-尼克尔森方案。

在下面的实现中，我们选择Lax-Friedrichs通量的函数 $\mathbf H$ ，即 $\mathbf{H}(\mathbf{a},\mathbf{b},\mathbf{n}) =
\frac{1}{2}(\mathbf{F}(\mathbf{a})\cdot \mathbf{n} +
\mathbf{F}(\mathbf{b})\cdot \mathbf{n} + \alpha (\mathbf{a} - \mathbf{b}))$ ，其中 $\alpha$ 是输入文件中指定的一个固定数字，或者 $\alpha$ 是一个与网格有关的值。在后一种情况下，它被选为 $\frac{h}{2\delta T}$ ， $h$ 是施加磁通量的面的直径，而 $\delta T$ 是当前的时间步长。

有了这些选择，将残差等同于零就会产生一个非线性方程组  $R(\mathbf{W}_{n+1})=0$  。我们通过牛顿迭代来解决这个非线性系统（与步骤15中解释的方法相同），即通过迭代

@f{eqnarray*}
R'(\mathbf{W}^k_{n+1},\delta \mathbf{W}_{n+1}^k)(\mathbf z) & = & -
R(\mathbf{W}^{k}_{n+1})(\mathbf z) \qquad \qquad \forall \mathbf z\in V_h \\
\mathbf{W}^{k+1}_{n+1} &=& \mathbf{W}^k_{n+1} + \delta \mathbf{W}^k_{n+1},


@f}

直到 $|R(\mathbf{W}^k_{n+1})|$ （残差）足够小。通过用有限元空间的节点基础而不是所有的 $\mathbf z$ 进行测试，我们得出了一个 $\delta \mathbf W$ 的线性系统。

@f{eqnarray*}
\mathbf R'(\mathbf{W}^k_{n+1})\delta \mathbf{W}^k_{n+1} & = & -
\mathbf R(\mathbf{W}^{k}_{n+1}).


@f}

一般来说，这个线性系统既不是对称的，也没有任何特定的确定性属性。我们将使用直接求解器或Trilinos的GMRES实现来解决它。从<a href="#Results">results shown below</a>中可以看出，这种全隐式迭代收敛速度非常快（通常为3步），并具有牛顿方法所期望的二次收敛顺序。




<a name="Automaticdifferentiation"></a><h3> Automatic differentiation </h3>


由于计算雅各布矩阵 $\mathbf R'(\mathbf W^k)$ 是一个可怕的野兽，我们使用一个自动微分包，Sacado，来做这个。  Sacado是<a
href="http://trilinos.org" target="_top">Trilinos</a>框架内的一个包，提供了一个C++模板类 <code>Sacado::Fad::DFad</code> （ <code>Fad</code> 代表 "前向自动微分"），支持基本算术运算符和函数，如 <code> sqrt, sin, cos, pow, </code> 等。为了使用这个功能，人们声明一个这种类型的变量集合，然后将这个集合中的一些变量表示为自由度，其余的变量是独立变量的函数。  这些变量在算法中被使用，随着变量的使用，它们对自由度的敏感度被持续更新。

可以想象，对于整个雅各布矩阵来说，这可能是非常昂贵的：自变量的数量是 $\mathbf W^k$ ，因变量是向量 $\mathbf
R(\mathbf W^k)$ 的元素。这两个向量很容易有几万个元素或更多。  然而，需要注意的是，并非 $\mathbf R$ 的所有元素都依赖于 $\mathbf W^k$ 的所有元素：事实上， $\mathbf R$ 中的一个条目只依赖于 $\mathbf W^k$ 的一个元素，如果两个相应的形状函数重叠并以弱形式耦合。

具体来说，定义当前单元上的残差可能依赖的最小独立AD变量集是明智的：在每个元素上，我们定义那些对应于定义在这个单元上的自由度的独立变量（或者，如果我们必须计算单元之间的跳转项，则对应于定义在两个相邻单元上的自由度），而因变量是本地残差向量的元素。如果不这样做，即把<i>all</i>和 $\mathbf W^k$ 的元素定义为独立的，将导致大量零的计算非常昂贵：局部残差向量的元素几乎独立于解向量的所有元素，因此它们的导数为零；然而，试图计算这些零可以轻易地占用整个程序90%甚至更多的计算时间，正如这个程序首次编写几年后，一个学生无意中做的实验所示。


回到自动计算雅各布系数的问题上。作者将这种方法与手工编码的雅各布式并列使用，用于不可压缩的Navier-Stokes问题，发现Sacado方法与使用手工编码的雅各布式一样快，但无限简单，而且不容易出错。由于使用自动差分只需要编码残差 $R(\mathbf{W})$ ，确保代码的正确性和维护代码变得非常简单--雅各布矩阵 $\mathbf R'$ 基本上是由计算残差 $\mathbf
R$ 的同一代码计算的。

说了这么多，这里有一个非常简单的例子，显示Sacado如何被使用。

@code
#include <Sacado.hpp>
#include <iostream>


using fad_double = Sacado::Fad::DFad<double>;


main() {


  fad_double a,b,c;


  a = 1; b = 2;


  a.diff(0,2);  // Set a to be dof 0, in a 2-dof system.


  b.diff(1,2);  // Set b to be dof 1, in a 2-dof system.


  c = 2*a+cos(a*b);


  double *derivs = &c.fastAccessDx(0); // Access derivatives


  std::cout << "dc/da = " << derivs[0] << ", dc/db=" << derivs[1] << std::endl;


}
@endcode



输出的是 $c(a,b)=2a+\cos(ab)$ 在 $a=1,b=2$ 的导数 $\frac{\partial c(a,b)}{\partial a},
\frac{\partial c(a,b)}{\partial b}$  。

应该注意的是，Sacado提供了更多的自动差分功能，而不是本程序中使用的小子集。  然而，理解上面的例子就足以理解Sacado在这个欧拉流程序中的使用。

<a name="Trilinossolvers"></a><h3> Trilinos solvers </h3> 该程序使用Aztec迭代求解器或Amesos稀疏直接求解器，两者均由Trilinos包提供。  这个软件包本身就是为了用于并行程序而设计的，然而，它也可以像这里一样，轻松地用于串行程序。  Epetra软件包是基本的矢量/矩阵库，解算器是在此基础上建立的。  这个非常强大的包可以用来描述矢量的平行分布，并定义对这些矢量进行操作的稀疏矩阵。  请查看注释代码，了解更多关于这些求解器在例子中的使用细节。


<a name="Adaptivity"></a><h3> Adaptivity </h3> 这个例子使用了一个特别的细化指标，该指标在冲击类问题中显示出一定的作用，在包括下坡流的例子中也是如此。  我们根据密度的平方梯度进行细化。悬空节点的处理是通过计算不同细化水平的单元的数值通量来实现的，而不是像迄今为止的所有其他教程程序那样使用AffineConstraints类。  通过这种方式，这个例子结合了连续和DG的方法论。它还简化了Jacobian的生成，因为我们不必通过用于计算自由度的自动微分来跟踪受限自由度。


 @note  而这个程序是在2008年写的，我们不知道有什么出版物会真正使用这种方法。然而，A. Dedner、R. Kl&ouml;fkorn和M. Kr&auml;nkel最近的一篇论文（"Continuous Finite-Elements on Non-Conforming Grids Using Discontinuous Galerkin Stabilization", Proceedings of Finite Volumes for Complex Applications VII - Methods and Theoretical Aspects, Springer, 2014）接近。

此外，我们强制规定了细化水平的最大数量，以控制细化的程度。  根据作者的经验，对于与时间有关的问题的适应性，如果不注意的话，细化很容易导致仿真戛然而止，因为时间步长的限制，如果网格在领域的任何部分变得太细的话。  在这个例子中，细化的数量被限制，让用户指定在网格的任何地方出现的最大细化水平。  这样一来，细化就不会使模拟速度减慢到停滞不前。  当然，这纯粹是一种启发式的策略，如果作者的顾问听说了，作者可能会被永远放逐出有限元误差估计界。

<a name="Inputdeckinitialandboundaryconditions"></a><h3>Input deck, initial and boundary conditions</h3>


我们使用一个输入文件平台来驱动仿真。  通过这种方式，我们可以改变边界条件和其他重要的模拟属性，而不必重新编译。  关于格式的更多信息，请看<a href="#Results">results section</a>，在那里我们更详细地描述了一个输入文件的例子。

在以前的例子程序中，我们通常对初始和边界条件进行硬编码。在这个程序中，我们改用表达式解析器类FunctionParser，这样我们就可以在输入文件中指定一个通用表达式，并在运行时对其进行解析&mdash；这样，我们就可以改变初始条件而不需要重新编译程序。因此，在下面的程序中不会声明名为InitialConditions或BoundaryConditions的类。




<a name="Implementation"></a><h3>Implementation</h3>


这个程序的实现被分成三个基本部分。   <ul>   <li>   <code>EulerEquations</code> 类，封装了完全描述欧拉方程具体内容的一切。这包括通量矩阵 $\mathbf F(\mathbf W)$ 、数值通量 $\mathbf F(\mathbf
  W^+,\mathbf W^-,\mathbf n)$ 、右手边 $\mathbf G(\mathbf W)$ 、边界条件、细化指标、输出的后处理，以及需要了解解向量和方程的各个组成部分的含义的类似事情。

    <li>  一个命名空间，处理与运行时参数有关的一切。

    <li>   <code>ConservationLaw</code>  处理时间步进、外部非线性和内部线性求解、组装线性系统以及驱动所有这些的顶层逻辑的类。   </ul> 

这种方法的原因是它将程序中的各种问题分开： <code>ConservationLaw</code> 是以这样一种方式编写的，即相对简单地将其适用于不同的方程组。人们只需为其他双曲方程重新实现 <code>EulerEquations</code> 类的成员，或者用额外的方程来增加现有的方程（例如通过添加额外的变量，或者通过添加化学成分等）。然而，这种修改不会影响到时间步进，或者非线性求解器，如果正确的话，因此 <code>ConservationLaw</code> 中的任何内容都不必修改。

同样，如果我们想改进线性或非线性求解器，或者改进时间步进方案（正如在<a
href="#Results">results section</a>的末尾所暗示的），那么这根本不需要对 <code>EulerEquations</code> 进行修改。


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
 * 首先是一套标准的deal.II包括。这里没有什么特别需要评论的。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/base/function_parser.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/conditional_ostream.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/grid_in.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 *  
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/mapping_q1.h> 
 * #include <deal.II/fe/fe_q.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * @endcode
 * 
 * 然后，正如介绍中提到的，我们使用各种Trilinos软件包作为线性求解器以及自动微分。这些都在以下的包含文件中。
 * 

 * 
 * 由于deal.II提供了基本的Trilinos矩阵、预处理程序和求解器的接口，我们把它们作为deal.II线性代数结构类似地包括在内。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * #include <deal.II/lac/trilinos_solver.h> 
 * 
 * @endcode
 * 
 * Sacado是Trilinos中的自动微分包，用于寻找全隐式牛顿迭代的雅各布系数。
 * 

 * 
 * 
 * @code
 * #include <Sacado.hpp> 
 * 
 * @endcode
 * 
 * 这又是C++语言。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * #include <vector> 
 * #include <memory> 
 * #include <array> 
 * 
 * @endcode
 * 
 * 在本节结束时，将dealii库中的所有内容引入本程序内容将进入的命名空间。
 * 

 * 
 * 
 * @code
 * namespace Step33 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Eulerequationspecifics"></a> 
 * <h3>Euler equation specifics</h3>
 * 

 * 
 * 这里我们定义了这个特定的守恒定律系统的通量函数，以及几乎所有其他的气体动力学欧拉方程所特有的东西，原因在介绍中讨论过。我们将所有这些归入一个结构，该结构定义了所有与通量有关的东西。这个结构的所有成员都是静态的，也就是说，这个结构没有由实例成员变量指定的实际状态。更好的方法是使用命名空间，而不是一个拥有所有静态成员的结构--但是命名空间不能被模板化，而且我们希望结构中的一些成员变量取决于空间维度，我们以通常的方式用模板参数来引入。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct EulerEquations 
 *   { 
 * @endcode
 * 
 * 
 * <a name="Componentdescription"></a> 
 * <h4>Component description</h4>
 * 

 * 
 * 首先是几个变量，它们以一种通用的方式描述了我们的解向量的各个组成部分。这包括系统中分量的数量（欧拉方程中每个空间方向的动量都有一个条目，加上能量和密度分量，总共有 <code>dim+2</code> 个分量），以及描述第一个动量分量、密度分量和能量密度分量在解向量中的索引的函数。请注意，所有这些%数都取决于空间维度；以通用的方式定义它们（而不是以隐含的惯例）使我们的代码更加灵活，并使以后的扩展更加容易，例如，在方程中加入更多的分量。
 * 

 * 
 * 
 * @code
 *     static const unsigned int n_components             = dim + 2; 
 *     static const unsigned int first_momentum_component = 0; 
 *     static const unsigned int density_component        = dim; 
 *     static const unsigned int energy_component         = dim + 1; 
 * 
 * @endcode
 * 
 * 在这个程序中一路生成图形输出时，我们需要指定解变量的名称，以及各种成分如何分组为矢量和标量场。我们可以在这里进行描述，但是为了使与欧拉方程有关的事情在这里得到解决，并使程序的其他部分尽可能地通用，我们在以下两个函数中提供了这类信息。
 * 

 * 
 * 
 * @code
 *     static std::vector<std::string> component_names() 
 *     { 
 *       std::vector<std::string> names(dim, "momentum"); 
 *       names.emplace_back("density"); 
 *       names.emplace_back("energy_density"); 
 * 
 *       return names; 
 *     } 
 * 
 *     static std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *     component_interpretation() 
 *     { 
 *       std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *         data_component_interpretation( 
 *           dim, DataComponentInterpretation::component_is_part_of_vector); 
 *       data_component_interpretation.push_back( 
 *         DataComponentInterpretation::component_is_scalar); 
 *       data_component_interpretation.push_back( 
 *         DataComponentInterpretation::component_is_scalar); 
 * 
 *       return data_component_interpretation; 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Transformationsbetweenvariables"></a> 
 * <h4>Transformations between variables</h4>
 * 

 * 
 * 接下来，我们定义气体常数。我们将在紧接着这个类的声明之后的定义中把它设置为1.4（与整数变量不同，比如上面的变量，静态常量浮点成员变量不能在C++的类声明中被初始化）。这个1.4的值代表了由两个原子组成的分子的气体，比如空气，它几乎完全由 $N_2$ 和 $O_2$ 组成，痕迹很小。
 * 

 * 
 * 
 * @code
 *     static const double gas_gamma; 
 * 
 * @endcode
 * 
 * 在下文中，我们将需要从保守变量的矢量中计算动能和压力。我们可以根据能量密度和动能 $\frac 12 \rho |\mathbf v|^2= \frac{|\rho \mathbf v|^2}{2\rho}$ 来做这件事（注意，独立变量包含动量分量 $\rho v_i$ ，而不是速度 $v_i$ ）。
 * 

 * 
 * 
 * @code
 *     template <typename InputVector> 
 *     static typename InputVector::value_type 
 *     compute_kinetic_energy(const InputVector &W) 
 *     { 
 *       typename InputVector::value_type kinetic_energy = 0; 
 *       for (unsigned int d = 0; d < dim; ++d) 
 *         kinetic_energy += 
 *           W[first_momentum_component + d] * W[first_momentum_component + d]; 
 *       kinetic_energy *= 1. / (2 * W[density_component]); 
 * 
 *       return kinetic_energy; 
 *     } 
 * 
 *     template <typename InputVector> 
 *     static typename InputVector::value_type 
 *     compute_pressure(const InputVector &W) 
 *     { 
 *       return ((gas_gamma - 1.0) * 
 *               (W[energy_component] - compute_kinetic_energy(W))); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="EulerEquationscompute_flux_matrix"></a> 
 * <h4>EulerEquations::compute_flux_matrix</h4>
 * 

 * 
 * 我们将通量函数 $F(W)$ 定义为一个大矩阵。 这个矩阵的每一行都代表了该行成分的标量守恒定律。 这个矩阵的确切形式在介绍中给出。请注意，我们知道这个矩阵的大小：它的行数与系统的分量一样多， <code>dim</code> 列数一样多；我们没有为这样的矩阵使用FullMatrix对象（它的行数和列数是可变的，因此每次创建这样的矩阵时必须在堆上分配内存），而是马上使用一个矩形的数字阵列。
 * 

 * 
 * 我们将通量函数的数值类型模板化，这样我们就可以在这里使用自动微分类型。 同样地，我们将用不同的输入矢量数据类型来调用该函数，所以我们也对其进行模板化。
 * 

 * 
 * 
 * @code
 *     template <typename InputVector> 
 *     static void compute_flux_matrix(const InputVector &W, 
 *                                     ndarray<typename InputVector::value_type, 
 *                                             EulerEquations<dim>::n_components, 
 *                                             dim> &     flux) 
 *     { 
 * 
 * @endcode
 * 
 * 首先计算出现在通量矩阵中的压力，然后计算矩阵中对应于动量项的前 <code>dim</code> 列。
 * 

 * 
 * 
 * @code
 *       const typename InputVector::value_type pressure = compute_pressure(W); 
 * 
 *       for (unsigned int d = 0; d < dim; ++d) 
 *         { 
 *           for (unsigned int e = 0; e < dim; ++e) 
 *             flux[first_momentum_component + d][e] = 
 *               W[first_momentum_component + d] * 
 *               W[first_momentum_component + e] / W[density_component]; 
 * 
 *           flux[first_momentum_component + d][d] += pressure; 
 *         } 
 * 
 * @endcode
 * 
 * 然后是密度（即质量守恒）的条款，最后是能量守恒。
 * 

 * 
 * 
 * @code
 *       for (unsigned int d = 0; d < dim; ++d) 
 *         flux[density_component][d] = W[first_momentum_component + d]; 
 * 
 *       for (unsigned int d = 0; d < dim; ++d) 
 *         flux[energy_component][d] = W[first_momentum_component + d] / 
 *                                     W[density_component] * 
 *                                     (W[energy_component] + pressure); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="EulerEquationscompute_normal_flux"></a> 
 * <h4>EulerEquations::compute_normal_flux</h4>
 * 

 * 
 * 在域的边界和跨挂节点上，我们使用一个数值通量函数来强制执行边界条件。 这个程序是基本的Lax-Friedrich的通量，有一个稳定的参数  $\alpha$  。它的形式也已经在介绍中给出。
 * 

 * 
 * 
 * @code
 *     template <typename InputVector> 
 *     static void numerical_normal_flux( 
 *       const Tensor<1, dim> &                                      normal, 
 *       const InputVector &                                         Wplus, 
 *       const InputVector &                                         Wminus, 
 *       const double                                                alpha, 
 *       std::array<typename InputVector::value_type, n_components> &normal_flux) 
 *     { 
 *       ndarray<typename InputVector::value_type, 
 *               EulerEquations<dim>::n_components, 
 *               dim> 
 *         iflux, oflux; 
 * 
 *       compute_flux_matrix(Wplus, iflux); 
 *       compute_flux_matrix(Wminus, oflux); 
 * 
 *       for (unsigned int di = 0; di < n_components; ++di) 
 *         { 
 *           normal_flux[di] = 0; 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d]; 
 * 
 *           normal_flux[di] += 0.5 * alpha * (Wplus[di] - Wminus[di]); 
 *         } 
 *     } 
 * @endcode
 * 
 * 
 * <a name="EulerEquationscompute_forcing_vector"></a> 
 * <h4>EulerEquations::compute_forcing_vector</h4>
 * 

 * 
 * 与描述通量函数 $\mathbf F(\mathbf w)$ 的方式相同，我们也需要有一种方法来描述右侧的强迫项。正如介绍中提到的，我们在这里只考虑重力，这导致了具体的形式 $\mathbf G(\mathbf w) = \left( g_1\rho, g_2\rho, g_3\rho, 0, \rho \mathbf g \cdot \mathbf v \right)^T$ ，这里显示的是三维情况。更具体地说，我们将只考虑三维的 $\mathbf g=(0,0,-1)^T$ ，或二维的 $\mathbf g=(0,-1)^T$ 。这自然导致了以下函数。
 * 

 * 
 * 
 * @code
 *     template <typename InputVector> 
 *     static void compute_forcing_vector( 
 *       const InputVector &                                         W, 
 *       std::array<typename InputVector::value_type, n_components> &forcing) 
 *     { 
 *       const double gravity = -1.0; 
 * 
 *       for (unsigned int c = 0; c < n_components; ++c) 
 *         switch (c) 
 *           { 
 *             case first_momentum_component + dim - 1: 
 *               forcing[c] = gravity * W[density_component]; 
 *               break; 
 *             case energy_component: 
 *               forcing[c] = gravity * W[first_momentum_component + dim - 1]; 
 *               break; 
 *             default: 
 *               forcing[c] = 0; 
 *           } 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Dealingwithboundaryconditions"></a> 
 * <h4>Dealing with boundary conditions</h4>
 * 

 * 
 * 我们必须处理的另一件事是边界条件。为此，让我们首先定义一下我们目前知道如何处理的各种边界条件。
 * 

 * 
 * 
 * @code
 *     enum BoundaryKind 
 *     { 
 *       inflow_boundary, 
 *       outflow_boundary, 
 *       no_penetration_boundary, 
 *       pressure_boundary 
 *     }; 
 * 
 * @endcode
 * 
 * 接下来的部分是实际决定在每一种边界上做什么。为此，请记住，从介绍中可以看出，边界条件是通过在给定的不均匀性 $\mathbf j$ 的边界外侧选择一个值 $\mathbf w^-$ ，以及可能在内部选择解的值 $\mathbf w^+$ 来指定的。然后，两者都被传递给数值通量 $\mathbf H(\mathbf{w}^+, \mathbf{w}^-, \mathbf{n})$ ，以定义边界对双线性形式的贡献。
 * 

 * 
 * 边界条件在某些情况下可以为解矢量的每个分量独立指定。例如，如果分量 $c$ 被标记为流入，那么 $w^-_c = j_c$  。如果是流出，那么 $w^-_c = w^+_c$  。这两种简单的情况在下面的函数中首先得到处理。
 * 

 * 
 * 有一个小插曲，从C++语言的角度来看，这个函数是不愉快的。输出向量  <code>Wminus</code>  当然会被修改，所以它不应该是  <code>const</code>  的参数。然而，在下面的实现中，它却成为了参数，而且为了使代码能够编译，它必须成为参数。原因是我们在 <code>Wminus</code> 类型为 <code>Table@<2,Sacado::Fad::DFad@<double@> @></code> 的地方调用这个函数，这是一个2d表，其指数分别代表正交点和向量分量。我们用 <code>Wminus[q]</code> 作为最后一个参数来调用这个函数；对2d表进行下标会产生一个代表1d向量的临时访问器对象，这正是我们在这里想要的。问题是，根据C++ 1998和2003标准，临时访问器对象不能被绑定到一个函数的非静态引用参数上，就像我们在这里希望的那样（这个问题将在下一个标准中以rvalue引用的形式得到解决）。 我们在这里把输出参数变成常量，是因为<i>accessor</i>对象是常量，而不是它所指向的表：那个表仍然可以被写到。然而，这个黑客是不愉快的，因为它限制了可以作为这个函数的模板参数的数据类型：一个普通的向量是不行的，因为当标记为  <code>const</code>  时，不能被写入。由于目前没有好的解决方案，我们将采用这里显示的务实的，甚至是不漂亮的解决方案。
 * 

 * 
 * 
 * @code
 *     template <typename DataVector> 
 *     static void 
 *     compute_Wminus(const std::array<BoundaryKind, n_components> &boundary_kind, 
 *                    const Tensor<1, dim> &                        normal_vector, 
 *                    const DataVector &                            Wplus, 
 *                    const Vector<double> &boundary_values, 
 *                    const DataVector &    Wminus) 
 *     { 
 *       for (unsigned int c = 0; c < n_components; c++) 
 *         switch (boundary_kind[c]) 
 *           { 
 *             case inflow_boundary: 
 *               { 
 *                 Wminus[c] = boundary_values(c); 
 *                 break; 
 *               } 
 * 
 *             case outflow_boundary: 
 *               { 
 *                 Wminus[c] = Wplus[c]; 
 *                 break; 
 *               } 
 * 
 * @endcode
 * 
 * 规定的压力边界条件有点复杂，因为即使压力是规定的，我们在这里真正设定的是能量分量，它将取决于速度和压力。因此，尽管这似乎是一个Dirichlet类型的边界条件，但我们得到了能量对速度和密度的敏感性（除非这些也被规定了）。
 * 

 * 
 * 
 * @code
 *             case pressure_boundary: 
 *               { 
 *                 const typename DataVector::value_type density = 
 *                   (boundary_kind[density_component] == inflow_boundary ? 
 *                      boundary_values(density_component) : 
 *                      Wplus[density_component]); 
 * 
 *                 typename DataVector::value_type kinetic_energy = 0; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   if (boundary_kind[d] == inflow_boundary) 
 *                     kinetic_energy += boundary_values(d) * boundary_values(d); 
 *                   else 
 *                     kinetic_energy += Wplus[d] * Wplus[d]; 
 *                 kinetic_energy *= 1. / 2. / density; 
 * 
 *                 Wminus[c] = 
 *                   boundary_values(c) / (gas_gamma - 1.0) + kinetic_energy; 
 * 
 *                 break; 
 *               } 
 * 
 *             case no_penetration_boundary: 
 *               { 
 * 
 * @endcode
 * 
 * 我们规定了速度（我们在这里处理的是一个特定的分量，所以速度的平均值是与表面法线正交的。 这就形成了整个速度分量的敏感度。
 * 

 * 
 * 
 * @code
 *                 typename DataVector::value_type vdotn = 0; 
 *                 for (unsigned int d = 0; d < dim; d++) 
 *                   { 
 *                     vdotn += Wplus[d] * normal_vector[d]; 
 *                   } 
 * 
 *                 Wminus[c] = Wplus[c] - 2.0 * vdotn * normal_vector[c]; 
 *                 break; 
 *               } 
 * 
 *             default: 
 *               Assert(false, ExcNotImplemented()); 
 *           } 
 *     } 
 * @endcode
 * 
 * 
 * <a name="EulerEquationscompute_refinement_indicators"></a> 
 * <h4>EulerEquations::compute_refinement_indicators</h4>
 * 

 * 
 * 在这个类中，我们也要指定如何细化网格。这个类 <code>ConservationLaw</code> 将使用我们在 <code>EulerEquation</code> 类中提供的所有信息，对于它所求解的特定守恒定律是不可知的：它甚至不关心一个求解向量有多少个分量。因此，它不可能知道合理的细化指标是什么。另一方面，在这里我们知道，或者至少我们可以想出一个合理的选择：我们简单地看一下密度的梯度，然后计算  $\eta_K=\log\left(1+|\nabla\rho(x_K)|\right)$  ，其中  $x_K$  是单元格  $K$  的中心。
 * 

 * 
 * 当然也有很多同样合理的细化指标，但这个指标确实如此，而且很容易计算。
 * 

 * 
 * 
 * @code
 *     static void 
 *     compute_refinement_indicators(const DoFHandler<dim> &dof_handler, 
 *                                   const Mapping<dim> &   mapping, 
 *                                   const Vector<double> & solution, 
 *                                   Vector<double> &       refinement_indicators) 
 *     { 
 *       const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
 *       std::vector<unsigned int> dofs(dofs_per_cell); 
 * 
 *       const QMidpoint<dim> quadrature_formula; 
 *       const UpdateFlags    update_flags = update_gradients; 
 *       FEValues<dim>        fe_v(mapping, 
 *                          dof_handler.get_fe(), 
 *                          quadrature_formula, 
 *                          update_flags); 
 * 
 *       std::vector<std::vector<Tensor<1, dim>>> dU( 
 *         1, std::vector<Tensor<1, dim>>(n_components)); 
 * 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         { 
 *           const unsigned int cell_no = cell->active_cell_index(); 
 *           fe_v.reinit(cell); 
 *           fe_v.get_function_gradients(solution, dU); 
 * 
 *           refinement_indicators(cell_no) = std::log( 
 *             1 + std::sqrt(dU[0][density_component] * dU[0][density_component])); 
 *         } 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="EulerEquationsPostprocessor"></a> 
 * <h4>EulerEquations::Postprocessor</h4>
 * 

 * 
 * 最后，我们声明一个实现数据组件后处理的类。这个类解决的问题是，我们使用的欧拉方程的表述中的变量是保守的而不是物理形式的：它们是动量密度  $\mathbf m=\rho\mathbf v$  ，密度  $\rho$  ，和能量密度  $E$  。我们还想把速度  $\mathbf v=\frac{\mathbf m}{\rho}$  和压力  $p=(\gamma-1)(E-\frac{1}{2} \rho |\mathbf v|^2)$  放入我们的输出文件中。
 * 

 * 
 * 此外，我们还想增加生成Schlieren图的可能性。Schlieren图是一种将冲击和其他尖锐界面可视化的方法。Schlieren "这个词是一个德语单词，可以翻译成 "条纹"--不过，用一个例子来解释可能更简单：比如说，当你把高浓度的酒精或透明的盐水倒入水中时，你会看到schlieren；这两种物质的颜色相同，但它们的折射率不同，所以在它们完全混合之前，光线会沿着弯曲的光线穿过混合物，如果你看它，会导致亮度变化。这就是 "分光"。类似的效果发生在可压缩流中，因为折射率取决于气体的压力（以及因此的密度）。
 * 

 * 
 * 这个词的起源是指三维体积的二维投影（我们看到的是三维流体的二维图片）。在计算流体力学中，我们可以通过考虑其原因来了解这种效应：密度变化。因此，Schlieren图是通过绘制 $s=|\nabla \rho|^2$ 产生的；显然， $s$ 在冲击和其他高度动态的地方很大。如果用户需要（通过在输入文件中指定），我们希望除了上面列出的其他派生量之外，还能生成这些希里伦图。
 * 

 * 
 * 从解决我们问题的数量中计算出派生数量，并将其输出到数据文件中的算法的实现依赖于DataPostprocessor类。它有大量的文档，该类的其他用途也可以在  step-29  中找到。因此，我们避免了大量的评论。
 * 

 * 
 * 
 * @code
 *     class Postprocessor : public DataPostprocessor<dim> 
 *     { 
 *     public: 
 *       Postprocessor(const bool do_schlieren_plot); 
 * 
 *       virtual void evaluate_vector_field( 
 *         const DataPostprocessorInputs::Vector<dim> &inputs, 
 *         std::vector<Vector<double>> &computed_quantities) const override; 
 * 
 *       virtual std::vector<std::string> get_names() const override; 
 * 
 *       virtual std::vector< 
 *         DataComponentInterpretation::DataComponentInterpretation> 
 *       get_data_component_interpretation() const override; 
 * 
 *       virtual UpdateFlags get_needed_update_flags() const override; 
 * 
 *     private: 
 *       const bool do_schlieren_plot; 
 *     }; 
 *   }; 
 * 
 *   template <int dim> 
 *   const double EulerEquations<dim>::gas_gamma = 1.4; 
 * 
 *   template <int dim> 
 *   EulerEquations<dim>::Postprocessor::Postprocessor( 
 *     const bool do_schlieren_plot) 
 *     : do_schlieren_plot(do_schlieren_plot) 
 *   {} 
 * 
 * @endcode
 * 
 * 这是唯一值得评论的函数。在生成图形输出时，DataOut和相关的类将在每个单元格上调用这个函数，以获取每个正交点的值、梯度、Hessians和法向量（如果我们在处理面）。请注意，每个正交点的数据本身就是矢量值，即保守变量。我们在这里要做的是计算每个正交点上我们感兴趣的量。注意，为此我们可以忽略Hessians（"inputs.solution_hessians"）和法向量（"inputs.normals"）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EulerEquations<dim>::Postprocessor::evaluate_vector_field( 
 *     const DataPostprocessorInputs::Vector<dim> &inputs, 
 *     std::vector<Vector<double>> &               computed_quantities) const 
 *   { 
 * 
 * @endcode
 * 
 * 在函数的开始，让我们确保所有的变量都有正确的大小，这样我们就可以访问各个向量元素，而不必怀疑我们是否可能读或写无效的元素；我们还检查 <code>solution_gradients</code> 向量只包含我们真正需要的数据（系统知道这个，因为我们在下面的 <code>get_needed_update_flags()</code> 函数中这样说）。对于内向量，我们检查至少外向量的第一个元素有正确的内部大小。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_quadrature_points = inputs.solution_values.size(); 
 * 
 *     if (do_schlieren_plot == true) 
 *       Assert(inputs.solution_gradients.size() == n_quadrature_points, 
 *              ExcInternalError()); 
 * 
 *     Assert(computed_quantities.size() == n_quadrature_points, 
 *            ExcInternalError()); 
 * 
 *     Assert(inputs.solution_values[0].size() == n_components, 
 *            ExcInternalError()); 
 * 
 *     if (do_schlieren_plot == true) 
 *       { 
 *         Assert(computed_quantities[0].size() == dim + 2, ExcInternalError()); 
 *       } 
 *     else 
 *       { 
 *         Assert(computed_quantities[0].size() == dim + 1, ExcInternalError()); 
 *       } 
 * 
 * @endcode
 * 
 * 然后在所有的正交点上循环，在那里做我们的工作。这段代码应该是不言自明的。输出变量的顺序首先是 <code>dim</code> 速度，然后是压力，如果需要的话，还可以是SCHLIEREN图。请注意，我们尝试使用 <code>first_momentum_component</code> 和 <code>density_component</code> 的信息，对输入向量中的变量顺序进行通用处理。
 * 

 * 
 * 
 * @code
 *     for (unsigned int q = 0; q < n_quadrature_points; ++q) 
 *       { 
 *         const double density = inputs.solution_values[q](density_component); 
 * 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           computed_quantities[q](d) = 
 *             inputs.solution_values[q](first_momentum_component + d) / density; 
 * 
 *         computed_quantities[q](dim) = 
 *           compute_pressure(inputs.solution_values[q]); 
 * 
 *         if (do_schlieren_plot == true) 
 *           computed_quantities[q](dim + 1) = 
 *             inputs.solution_gradients[q][density_component] * 
 *             inputs.solution_gradients[q][density_component]; 
 *       } 
 *   } 
 * 
 *   template <int dim> 
 *   std::vector<std::string> EulerEquations<dim>::Postprocessor::get_names() const 
 *   { 
 *     std::vector<std::string> names; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       names.emplace_back("velocity"); 
 *     names.emplace_back("pressure"); 
 * 
 *     if (do_schlieren_plot == true) 
 *       names.emplace_back("schlieren_plot"); 
 * 
 *     return names; 
 *   } 
 * 
 *   template <int dim> 
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *   EulerEquations<dim>::Postprocessor::get_data_component_interpretation() const 
 *   { 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       interpretation(dim, 
 *                      DataComponentInterpretation::component_is_part_of_vector); 
 * 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 * 
 *     if (do_schlieren_plot == true) 
 *       interpretation.push_back( 
 *         DataComponentInterpretation::component_is_scalar); 
 * 
 *     return interpretation; 
 *   } 
 * 
 *   template <int dim> 
 *   UpdateFlags 
 *   EulerEquations<dim>::Postprocessor::get_needed_update_flags() const 
 *   { 
 *     if (do_schlieren_plot == true) 
 *       return update_values | update_gradients; 
 *     else 
 *       return update_values; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameterhandling"></a> 
 * <h3>Run time parameter handling</h3>
 * 

 * 
 * 我们接下来的工作是定义一些包含运行时参数的类（例如求解器的公差、迭代次数、稳定参数等等）。我们可以在主类中完成这项工作，但我们将其与主类分开，以使程序更加模块化和易于阅读。所有与运行时参数有关的东西都在以下命名空间中，而程序逻辑则在主类中。
 * 

 * 
 * 我们将把运行时参数分成几个独立的结构，我们将把这些结构全部放在一个命名空间  <code>Parameters</code>  中。在这些类中，有几个类将参数分组，用于单独的组，比如用于求解器、网格细化或输出。这些类中的每一个都有函数  <code>declare_parameters()</code>  和  <code>parse_parameters()</code>  ，分别在ParameterHandler对象中声明参数子段和条目，并从这样的对象中检索实际参数值。这些类在ParameterHandler的子段中声明它们的所有参数。
 * 

 * 
 * 以下命名空间的最后一个类结合了前面所有的类，从它们派生出来，并负责处理输入文件顶层的一些条目，以及其他一些奇特的条目，这些条目在子段中太短了，不值得本身有一个结构。
 * 

 * 
 * 这里值得指出的是一件事。下面这些类中没有一个构造函数可以初始化各种成员变量。不过这不是问题，因为我们将从输入文件中读取这些类中声明的所有变量（或者间接地：一个ParameterHandler对象将从那里读取，我们将从这个对象中获取数值），它们将以这种方式被初始化。如果输入文件中根本没有指定某个变量，这也不是问题。在这种情况下，ParameterHandler类将简单地采取默认值，这个默认值是在声明下面这些类的 <code>declare_parameters()</code> 函数中的一个条目时指定的。
 * 

 * 
 * 
 * @code
 *   namespace Parameters 
 *   { 
 * @endcode
 * 
 * 
 * <a name="ParametersSolver"></a> 
 * <h4>Parameters::Solver</h4>
 * 

 * 
 * 这些类中的第一个是关于线性内部求解器的参数。它提供的参数表明使用哪种求解器（GMRES作为一般非对称不定式系统的求解器，或稀疏直接求解器），要产生的输出量，以及各种调整阈值不完全LU分解（ILUT）的参数，我们使用它作为GMRES的预处理器。
 * 

 * 
 * 特别是，ILUT需要以下参数。
 * 

 * 
 * - ilut_fill：形成ILU分解时要增加的额外条目数
 * 

 * 
 * - ilut_atol, ilut_rtol: 在形成预处理程序时，对于某些问题，不好的条件（或者只是运气不好）会导致预处理程序的条件很差。 因此，将对角线扰动添加到原始矩阵中，并为这个稍好的矩阵形成预处理程序会有帮助。 ATOL是一个绝对扰动，在形成预处理之前加到对角线上，RTOL是一个比例因子  $rtol \geq 1$  。
 * 

 * 
 * - ilut_drop: ILUT将放弃任何幅度小于此值的数值。 这是一种管理该预处理程序所使用的内存量的方法。
 * 

 * 
 * 每个参数的含义在  ParameterHandler::declare_entry  调用的第三个参数中也有简要说明  <code>declare_parameters()</code>  。
 * 

 * 
 * 
 * @code
 *     struct Solver 
 *     { 
 *       enum SolverType 
 *       { 
 *         gmres, 
 *         direct 
 *       }; 
 *       SolverType solver; 
 * 
 *       enum OutputType 
 *       { 
 *         quiet, 
 *         verbose 
 *       }; 
 *       OutputType output; 
 * 
 *       double linear_residual; 
 *       int    max_iterations; 
 * 
 *       double ilut_fill; 
 *       double ilut_atol; 
 *       double ilut_rtol; 
 *       double ilut_drop; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Solver::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("linear solver"); 
 *       { 
 *         prm.declare_entry( 
 *           "output", 
 *           "quiet", 
 *           Patterns::Selection("quiet|verbose"), 
 *           "State whether output from solver runs should be printed. " 
 *           "Choices are <quiet|verbose>."); 
 *         prm.declare_entry("method", 
 *                           "gmres", 
 *                           Patterns::Selection("gmres|direct"), 
 *                           "The kind of solver for the linear system. " 
 *                           "Choices are <gmres|direct>."); 
 *         prm.declare_entry("residual", 
 *                           "1e-10", 
 *                           Patterns::Double(), 
 *                           "Linear solver residual"); 
 *         prm.declare_entry("max iters", 
 *                           "300", 
 *                           Patterns::Integer(), 
 *                           "Maximum solver iterations"); 
 *         prm.declare_entry("ilut fill", 
 *                           "2", 
 *                           Patterns::Double(), 
 *                           "Ilut preconditioner fill"); 
 *         prm.declare_entry("ilut absolute tolerance", 
 *                           "1e-9", 
 *                           Patterns::Double(), 
 *                           "Ilut preconditioner tolerance"); 
 *         prm.declare_entry("ilut relative tolerance", 
 *                           "1.1", 
 *                           Patterns::Double(), 
 *                           "Ilut relative tolerance"); 
 *         prm.declare_entry("ilut drop tolerance", 
 *                           "1e-10", 
 *                           Patterns::Double(), 
 *                           "Ilut drop tolerance"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Solver::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("linear solver"); 
 *       { 
 *         const std::string op = prm.get("output"); 
 *         if (op == "verbose") 
 *           output = verbose; 
 *         if (op == "quiet") 
 *           output = quiet; 
 * 
 *         const std::string sv = prm.get("method"); 
 *         if (sv == "direct") 
 *           solver = direct; 
 *         else if (sv == "gmres") 
 *           solver = gmres; 
 * 
 *         linear_residual = prm.get_double("residual"); 
 *         max_iterations  = prm.get_integer("max iters"); 
 *         ilut_fill       = prm.get_double("ilut fill"); 
 *         ilut_atol       = prm.get_double("ilut absolute tolerance"); 
 *         ilut_rtol       = prm.get_double("ilut relative tolerance"); 
 *         ilut_drop       = prm.get_double("ilut drop tolerance"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="ParametersRefinement"></a> 
 * <h4>Parameters::Refinement</h4>
 * 

 * 
 * 同样的，这里有几个参数决定了网格如何被细化（以及是否要被细化）。关于冲击参数的具体作用，请看下面的网格细化函数。
 * 

 * 
 * 
 * @code
 *     struct Refinement 
 *     { 
 *       bool   do_refine; 
 *       double shock_val; 
 *       double shock_levels; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Refinement::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("refinement"); 
 *       { 
 *         prm.declare_entry("refinement", 
 *                           "true", 
 *                           Patterns::Bool(), 
 *                           "Whether to perform mesh refinement or not"); 
 *         prm.declare_entry("refinement fraction", 
 *                           "0.1", 
 *                           Patterns::Double(), 
 *                           "Fraction of high refinement"); 
 *         prm.declare_entry("unrefinement fraction", 
 *                           "0.1", 
 *                           Patterns::Double(), 
 *                           "Fraction of low unrefinement"); 
 *         prm.declare_entry("max elements", 
 *                           "1000000", 
 *                           Patterns::Double(), 
 *                           "maximum number of elements"); 
 *         prm.declare_entry("shock value", 
 *                           "4.0", 
 *                           Patterns::Double(), 
 *                           "value for shock indicator"); 
 *         prm.declare_entry("shock levels", 
 *                           "3.0", 
 *                           Patterns::Double(), 
 *                           "number of shock refinement levels"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Refinement::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("refinement"); 
 *       { 
 *         do_refine    = prm.get_bool("refinement"); 
 *         shock_val    = prm.get_double("shock value"); 
 *         shock_levels = prm.get_double("shock levels"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="ParametersFlux"></a> 
 * <h4>Parameters::Flux</h4>
 * 

 * 
 * 接下来是关于通量修改的部分，使其更加稳定。特别是提供了两个选项来稳定Lax-Friedrichs通量：要么选择 $\mathbf{H}(\mathbf{a},\mathbf{b},\mathbf{n}) = \frac{1}{2}(\mathbf{F}(\mathbf{a})\cdot \mathbf{n} + \mathbf{F}(\mathbf{b})\cdot \mathbf{n} + \alpha (\mathbf{a} - \mathbf{b}))$ ，其中 $\alpha$ 是在输入文件中指定的一个固定数字，要么 $\alpha$ 是一个与网格有关的值。在后一种情况下，它被选择为 $\frac{h}{2\delta T}$ ，其中 $h$ 是施加流量的面的直径， $\delta T$ 是当前的时间步长。
 * 

 * 
 * 
 * @code
 *     struct Flux 
 *     { 
 *       enum StabilizationKind 
 *       { 
 *         constant, 
 *         mesh_dependent 
 *       }; 
 *       StabilizationKind stabilization_kind; 
 * 
 *       double stabilization_value; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Flux::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("flux"); 
 *       { 
 *         prm.declare_entry( 
 *           "stab", 
 *           "mesh", 
 *           Patterns::Selection("constant|mesh"), 
 *           "Whether to use a constant stabilization parameter or " 
 *           "a mesh-dependent one"); 
 *         prm.declare_entry("stab value", 
 *                           "1", 
 *                           Patterns::Double(), 
 *                           "alpha stabilization"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Flux::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("flux"); 
 *       { 
 *         const std::string stab = prm.get("stab"); 
 *         if (stab == "constant") 
 *           stabilization_kind = constant; 
 *         else if (stab == "mesh") 
 *           stabilization_kind = mesh_dependent; 
 *         else 
 *           AssertThrow(false, ExcNotImplemented()); 
 * 
 *         stabilization_value = prm.get_double("stab value"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="ParametersOutput"></a> 
 * <h4>Parameters::Output</h4>
 * 

 * 
 * 然后是关于输出参数的部分。我们提供产生Schlieren图（密度的平方梯度，一种可视化冲击前沿的工具），以及图形输出的时间间隔，以防我们不希望每个时间步骤都有输出文件。
 * 

 * 
 * 
 * @code
 *     struct Output 
 *     { 
 *       bool   schlieren_plot; 
 *       double output_step; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Output::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("output"); 
 *       { 
 *         prm.declare_entry("schlieren plot", 
 *                           "true", 
 *                           Patterns::Bool(), 
 *                           "Whether or not to produce schlieren plots"); 
 *         prm.declare_entry("step", 
 *                           "-1", 
 *                           Patterns::Double(), 
 *                           "Output once per this period"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Output::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("output"); 
 *       { 
 *         schlieren_plot = prm.get_bool("schlieren plot"); 
 *         output_step    = prm.get_double("step"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="ParametersAllParameters"></a> 
 * <h4>Parameters::AllParameters</h4>
 * 

 * 
 * 最后是将这一切结合起来的类。它自己声明了一些参数，主要是参数文件顶层的参数，以及一些太小的部分，以至于不值得有自己的类。它还包含了所有实际上与空间维度有关的东西，比如初始或边界条件。
 * 

 * 
 * 因为这个类是由上面所有的类派生出来的，所以 <code>declare_parameters()</code> and <code>parse_parameters()</code> 函数也会调用基类的相应函数。
 * 

 * 
 * 注意这个类也处理输入文件中指定的初始和边界条件的声明。为此，在这两种情况下，都有像 "w_0值 "这样的条目，它代表了 $x,y,z$ 方面的表达式，将初始或边界条件描述为一个公式，随后将由FunctionParser类来解析。类似的表达方式还有 "w_1"、"w_2 "等，表示欧拉系统的 <code>dim+2</code> 守恒变量。同样，我们允许在输入文件中最多使用 <code>max_n_boundaries</code> 个边界指标，这些边界指标中的每一个都可以与流入、流出或压力边界条件相关联，同质的边界条件要分别为每个组件和每个边界指标指定。
 * 

 * 
 * 用来存储边界指标的数据结构有点复杂。它是一个 <code>max_n_boundaries</code> 元素的数组，表示将被接受的边界指标的范围。对于这个数组中的每个条目，我们在 <code>BoundaryCondition</code> 结构中存储一对数据：首先是一个大小为 <code>n_components</code> 的数组，对于解向量的每个分量，它表明它是流入、流出还是其他类型的边界，其次是一个FunctionParser对象，它一次描述了这个边界ID的解向量的所有分量。
 * 

 * 
 * <code>BoundaryCondition</code> 结构需要一个构造器，因为我们需要在构造时告诉函数解析器对象它要描述多少个向量分量。因此，这个初始化不能等到我们在后面的 <code>AllParameters::parse_parameters()</code> 中实际设置FunctionParser对象所代表的公式。
 * 

 * 
 * 由于必须在构造时告诉Function对象其向量大小的同样原因，我们必须有一个 <code>AllParameters</code> 类的构造函数，至少要初始化另一个FunctionParser对象，即描述初始条件的对象。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     struct AllParameters : public Solver, 
 *                            public Refinement, 
 *                            public Flux, 
 *                            public Output 
 *     { 
 *       static const unsigned int max_n_boundaries = 10; 
 * 
 *       struct BoundaryConditions 
 *       { 
 *         std::array<typename EulerEquations<dim>::BoundaryKind, 
 *                    EulerEquations<dim>::n_components> 
 *           kind; 
 * 
 *         FunctionParser<dim> values; 
 * 
 *         BoundaryConditions(); 
 *       }; 
 * 
 *       AllParameters(); 
 * 
 *       double diffusion_power; 
 * 
 *       double time_step, final_time; 
 *       double theta; 
 *       bool   is_stationary; 
 * 
 *       std::string mesh_filename; 
 * 
 *  
 *       BoundaryConditions  boundary_conditions[max_n_boundaries]; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     template <int dim> 
 *     AllParameters<dim>::BoundaryConditions::BoundaryConditions() 
 *       : values(EulerEquations<dim>::n_components) 
 *     { 
 *       std::fill(kind.begin(), 
 *                 kind.end(), 
 *                 EulerEquations<dim>::no_penetration_boundary); 
 *     } 
 * 
 *     template <int dim> 
 *     AllParameters<dim>::AllParameters() 
 *       : diffusion_power(0.) 
 *       , time_step(1.) 
 *       , final_time(1.) 
 *       , theta(.5) 
 *       , is_stationary(true) 
 *       , initial_conditions(EulerEquations<dim>::n_components) 
 *     {} 
 * 
 *     template <int dim> 
 *     void AllParameters<dim>::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.declare_entry("mesh", 
 *                         "grid.inp", 
 *                         Patterns::Anything(), 
 *                         "input file name"); 
 * 
 *       prm.declare_entry("diffusion power", 
 *                         "2.0", 
 *                         Patterns::Double(), 
 *                         "power of mesh size for diffusion"); 
 * 
 *       prm.enter_subsection("time stepping"); 
 *       { 
 *         prm.declare_entry("time step", 
 *                           "0.1", 
 *                           Patterns::Double(0), 
 *                           "simulation time step"); 
 *         prm.declare_entry("final time", 
 *                           "10.0", 
 *                           Patterns::Double(0), 
 *                           "simulation end time"); 
 *         prm.declare_entry("theta scheme value", 
 *                           "0.5", 
 *                           Patterns::Double(0, 1), 
 *                           "value for theta that interpolated between explicit " 
 *                           "Euler (theta=0), Crank-Nicolson (theta=0.5), and " 
 *                           "implicit Euler (theta=1)."); 
 *       } 
 *       prm.leave_subsection(); 
 * 
 *       for (unsigned int b = 0; b < max_n_boundaries; ++b) 
 *         { 
 *           prm.enter_subsection("boundary_" + Utilities::int_to_string(b)); 
 *           { 
 *             prm.declare_entry("no penetration", 
 *                               "false", 
 *                               Patterns::Bool(), 
 *                               "whether the named boundary allows gas to " 
 *                               "penetrate or is a rigid wall"); 
 * 
 *             for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
 *                  ++di) 
 *               { 
 *                 prm.declare_entry("w_" + Utilities::int_to_string(di), 
 *                                   "outflow", 
 *                                   Patterns::Selection( 
 *                                     "inflow|outflow|pressure"), 
 *                                   "<inflow|outflow|pressure>"); 
 * 
 *                 prm.declare_entry("w_" + Utilities::int_to_string(di) + 
 *                                     " value", 
 *                                   "0.0", 
 *                                   Patterns::Anything(), 
 *                                   "expression in x,y,z"); 
 *               } 
 *           } 
 *           prm.leave_subsection(); 
 *         } 
 * 
 *       prm.enter_subsection("initial condition"); 
 *       { 
 *         for (unsigned int di = 0; di < EulerEquations<dim>::n_components; ++di) 
 *           prm.declare_entry("w_" + Utilities::int_to_string(di) + " value", 
 *                             "0.0", 
 *                             Patterns::Anything(), 
 *                             "expression in x,y,z"); 
 *       } 
 *       prm.leave_subsection(); 
 * 
 *       Parameters::Solver::declare_parameters(prm); 
 *       Parameters::Refinement::declare_parameters(prm); 
 *       Parameters::Flux::declare_parameters(prm); 
 *       Parameters::Output::declare_parameters(prm); 
 *     } 
 * 
 *     template <int dim> 
 *     void AllParameters<dim>::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       mesh_filename   = prm.get("mesh"); 
 *       diffusion_power = prm.get_double("diffusion power"); 
 * 
 *       prm.enter_subsection("time stepping"); 
 *       { 
 *         time_step = prm.get_double("time step"); 
 *         if (time_step == 0) 
 *           { 
 *             is_stationary = true; 
 *             time_step     = 1.0; 
 *             final_time    = 1.0; 
 *           } 
 *         else 
 *           is_stationary = false; 
 * 
 *         final_time = prm.get_double("final time"); 
 *         theta      = prm.get_double("theta scheme value"); 
 *       } 
 *       prm.leave_subsection(); 
 * 
 *       for (unsigned int boundary_id = 0; boundary_id < max_n_boundaries; 
 *            ++boundary_id) 
 *         { 
 *           prm.enter_subsection("boundary_" + 
 *                                Utilities::int_to_string(boundary_id)); 
 *           { 
 *             std::vector<std::string> expressions( 
 *               EulerEquations<dim>::n_components, "0.0"); 
 * 
 *             const bool no_penetration = prm.get_bool("no penetration"); 
 * 
 *             for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
 *                  ++di) 
 *               { 
 *                 const std::string boundary_type = 
 *                   prm.get("w_" + Utilities::int_to_string(di)); 
 * 
 *                 if ((di < dim) && (no_penetration == true)) 
 *                   boundary_conditions[boundary_id].kind[di] = 
 *                     EulerEquations<dim>::no_penetration_boundary; 
 *                 else if (boundary_type == "inflow") 
 *                   boundary_conditions[boundary_id].kind[di] = 
 *                     EulerEquations<dim>::inflow_boundary; 
 *                 else if (boundary_type == "pressure") 
 *                   boundary_conditions[boundary_id].kind[di] = 
 *                     EulerEquations<dim>::pressure_boundary; 
 *                 else if (boundary_type == "outflow") 
 *                   boundary_conditions[boundary_id].kind[di] = 
 *                     EulerEquations<dim>::outflow_boundary; 
 *                 else 
 *                   AssertThrow(false, ExcNotImplemented()); 
 * 
 *                 expressions[di] = 
 *                   prm.get("w_" + Utilities::int_to_string(di) + " value"); 
 *               } 
 * 
 *             boundary_conditions[boundary_id].values.initialize( 
 *               FunctionParser<dim>::default_variable_names(), 
 *               expressions, 
 *               std::map<std::string, double>()); 
 *           } 
 *           prm.leave_subsection(); 
 *         } 
 * 
 *       prm.enter_subsection("initial condition"); 
 *       { 
 *         std::vector<std::string> expressions(EulerEquations<dim>::n_components, 
 *                                              "0.0"); 
 *         for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) 
 *           expressions[di] = 
 *             prm.get("w_" + Utilities::int_to_string(di) + " value"); 
 *         initial_conditions.initialize( 
 *           FunctionParser<dim>::default_variable_names(), 
 *           expressions, 
 *           std::map<std::string, double>()); 
 *       } 
 *       prm.leave_subsection(); 
 * 
 *       Parameters::Solver::parse_parameters(prm); 
 *       Parameters::Refinement::parse_parameters(prm); 
 *       Parameters::Flux::parse_parameters(prm); 
 *       Parameters::Output::parse_parameters(prm); 
 *     } 
 *   } // namespace Parameters 
 * 
 * @endcode
 * 
 * 
 * <a name="Conservationlawclass"></a> 
 * <h3>Conservation law class</h3>
 * 

 * 
 * 这里终于出现了一个类，它实际上是对我们上面定义的所有欧拉方程和参数的具体内容做了一些事情。公共接口与以往基本相同（构造函数现在需要一个文件名来读取参数，这个文件名在命令行中传递）。私有函数接口也与通常的安排非常相似， <code>assemble_system</code> 函数被分成三个部分：一个包含所有单元的主循环，然后分别调用另外两个单元和面的积分。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ConservationLaw 
 *   { 
 *   public: 
 *     ConservationLaw(const char *input_filename); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 * 
 *     void assemble_system(); 
 *     void assemble_cell_term(const FEValues<dim> &                       fe_v, 
 *                             const std::vector<types::global_dof_index> &dofs); 
 *     void assemble_face_term( 
 *       const unsigned int                          face_no, 
 *       const FEFaceValuesBase<dim> &               fe_v, 
 *       const FEFaceValuesBase<dim> &               fe_v_neighbor, 
 *       const std::vector<types::global_dof_index> &dofs, 
 *       const std::vector<types::global_dof_index> &dofs_neighbor, 
 *       const bool                                  external_face, 
 *       const unsigned int                          boundary_id, 
 *       const double                                face_diameter); 
 * 
 *     std::pair<unsigned int, double> solve(Vector<double> &solution); 
 * 
 *     void compute_refinement_indicators(Vector<double> &indicator) const; 
 *     void refine_grid(const Vector<double> &indicator); 
 * 
 *     void output_results() const; 
 * 
 * @endcode
 * 
 * 前面的几个成员变量也相当标准。请注意，我们定义了一个映射对象，在整个程序中组装术语时使用（我们将把它交给每个FEValues和FEFaceValues对象）；我们使用的映射只是标准的 $Q_1$ 映射--换句话说，没有什么花哨的东西--但是在这里声明一个映射并在整个程序中使用它将使以后在有必要时改变它更加简单。事实上，这一点相当重要：众所周知，对于欧拉方程的跨音速模拟，如果边界近似没有足够高的阶数，计算就不会收敛，即使像 $h\rightarrow 0$ 那样。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim>   triangulation; 
 *     const MappingQ1<dim> mapping; 
 * 
 *     const FESystem<dim> fe; 
 *     DoFHandler<dim>     dof_handler; 
 * 
 *     const QGauss<dim>     quadrature; 
 *     const QGauss<dim - 1> face_quadrature; 
 * 
 * @endcode
 * 
 * 接下来是一些数据向量，对应于前一个时间步骤的解决方案（ <code>old_solution</code> ），当前解决方案的最佳猜测（ <code>current_solution</code> ；我们说<i>guess</i>是因为计算它的牛顿迭代可能还没有收敛，而 <code>old_solution</code> 是指前一个时间步骤的完全收敛的最终结果），以及下一个时间步骤的解决方案的预测器，通过将当前和之前的解决方案推算到未来一个时间步骤计算。
 * 

 * 
 * 
 * @code
 *     Vector<double> old_solution; 
 *     Vector<double> current_solution; 
 *     Vector<double> predictor; 
 * 
 *     Vector<double> right_hand_side; 
 * 
 * @endcode
 * 
 * 这一组最后的成员变量（除了最下面的保存所有运行时参数的对象和一个屏幕输出流，它只在要求verbose输出的情况下打印一些东西）处理我们在这个程序中与Trilinos库的接口，该库为我们提供了线性求解器。与在 step-17 和 step-18 中包括PETSc矩阵类似，我们需要做的是创建一个Trilinos稀疏矩阵而不是标准的deal.II类。该系统矩阵在每个牛顿步骤中被用于雅各布系数。由于我们不打算并行运行这个程序（不过用Trilinos数据结构也不难），所以我们不必考虑其他的事情，比如分配自由度。
 * 

 * 
 * 
 * @code
 *     TrilinosWrappers::SparseMatrix system_matrix; 
 * 
 *     Parameters::AllParameters<dim> parameters; 
 *     ConditionalOStream             verbose_cout; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ConservationLawConservationLaw"></a> 
 * <h4>ConservationLaw::ConservationLaw</h4>
 * 

 * 
 * 关于构造函数没有什么可说的。基本上，它读取输入文件并将解析后的值填充到参数对象中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ConservationLaw<dim>::ConservationLaw(const char *input_filename) 
 *     : mapping() 
 *     , fe(FE_Q<dim>(1), EulerEquations<dim>::n_components) 
 *     , dof_handler(triangulation) 
 *     , quadrature(fe.degree + 1) 
 *     , face_quadrature(fe.degree + 1) 
 *     , verbose_cout(std::cout, false) 
 *   { 
 *     ParameterHandler prm; 
 *     Parameters::AllParameters<dim>::declare_parameters(prm); 
 * 
 *     prm.parse_input(input_filename); 
 *     parameters.parse_parameters(prm); 
 * 
 *     verbose_cout.set_condition(parameters.output == 
 *                                Parameters::Solver::verbose); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ConservationLawsetup_system"></a> 
 * <h4>ConservationLaw::setup_system</h4>
 * 

 * 
 * 每次改变网格时都会调用下面这个（简单的）函数。它所做的就是根据我们在之前所有的教程程序中生成的稀疏模式来调整特里诺斯矩阵的大小。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::setup_system() 
 *   { 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 * 
 *     system_matrix.reinit(dsp); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawassemble_system"></a> 
 * <h4>ConservationLaw::assemble_system</h4>
 * 

 * 
 * 这个和下面两个函数是这个程序的核心。它们将牛顿方法应用于非线性守恒方程组所产生的线性系统组合起来。
 * 

 * 
 * 第一个函数将所有的装配部件放在一个例行程序中，为每个单元格/面分配正确的部件。 对这些对象的装配的实际实现是在以下函数中完成的。
 * 

 * 
 * 在函数的顶部，我们做了常规的内务处理：分配FEValues、FEFaceValues和FESubfaceValues对象，这些对象对单元、面和子面（在不同细化级别的相邻单元的情况下）进行积分。请注意，我们并不需要所有这些对象的所有信息（如值、梯度或正交点的实际位置），所以我们只让FEValues类通过指定最小的UpdateFlags集来获得实际需要的信息。例如，当使用邻接单元的FEFaceValues对象时，我们只需要形状值。给定一个特定的面，正交点和 <code>JxW</code> 值与当前单元格相同，法向量已知为当前单元格的法向量的负值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::assemble_system() 
 *   { 
 *     const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
 * 
 *     std::vector<types::global_dof_index> dof_indices(dofs_per_cell); 
 *     std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell); 
 * 
 *     const UpdateFlags update_flags = update_values | update_gradients | 
 *                                      update_quadrature_points | 
 *                                      update_JxW_values, 
 *                       face_update_flags = 
 *                         update_values | update_quadrature_points | 
 *                         update_JxW_values | update_normal_vectors, 
 *                       neighbor_face_update_flags = update_values; 
 * 
 *     FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags); 
 *     FEFaceValues<dim>    fe_v_face(mapping, 
 *                                 fe, 
 *                                 face_quadrature, 
 *                                 face_update_flags); 
 *     FESubfaceValues<dim> fe_v_subface(mapping, 
 *                                       fe, 
 *                                       face_quadrature, 
 *                                       face_update_flags); 
 *     FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
 *                                          fe, 
 *                                          face_quadrature, 
 *                                          neighbor_face_update_flags); 
 *     FESubfaceValues<dim> fe_v_subface_neighbor(mapping, 
 *                                                fe, 
 *                                                face_quadrature, 
 *                                                neighbor_face_update_flags); 
 * 
 * @endcode
 * 
 * 然后循环所有单元，初始化当前单元的FEValues对象，并调用在此单元上组装问题的函数。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_v.reinit(cell); 
 *         cell->get_dof_indices(dof_indices); 
 * 
 *         assemble_cell_term(fe_v, dof_indices); 
 * 
 * @endcode
 * 
 * 然后在这个单元的所有面上循环。 如果一个面是外部边界的一部分，那么就在那里集合边界条件（ <code>assemble_face_terms</code> 的第五个参数表示我们是在外部面还是内部面工作；如果是外部面，表示邻居自由度指数的第四个参数被忽略，所以我们传递一个空向量）。
 * 

 * 
 * 
 * @code
 *         for (const auto face_no : cell->face_indices()) 
 *           if (cell->at_boundary(face_no)) 
 *             { 
 *               fe_v_face.reinit(cell, face_no); 
 *               assemble_face_term(face_no, 
 *                                  fe_v_face, 
 *                                  fe_v_face, 
 *                                  dof_indices, 
 *                                  std::vector<types::global_dof_index>(), 
 *                                  true, 
 *                                  cell->face(face_no)->boundary_id(), 
 *                                  cell->face(face_no)->diameter()); 
 *             } 
 * 
 * @endcode
 * 
 * 另一种情况是，我们正在处理一个内部面。我们需要区分两种情况：这是在同一细化水平的两个单元之间的正常面，和在不同细化水平的两个单元之间的面。
 * 

 * 
 * 在第一种情况下，我们不需要做什么：我们使用的是连续有限元，在这种情况下，面条款不会出现在双线性表格中。第二种情况通常也不会导致面条款，如果我们强烈地执行悬挂节点约束的话（就像到目前为止，只要我们使用连续有限元的所有教程程序一样--这种执行是由AffineConstraints类和 DoFTools::make_hanging_node_constraints). 一起完成的）。 然而，在当前程序中，我们选择在不同细化水平的单元之间的面弱地执行连续性，原因有二。(i)因为我们可以，更重要的是(ii)因为我们必须通过AffineConstraints类的操作，将我们用来计算牛顿矩阵元素的自动微分穿起来。这是有可能的，但不是微不足道的，所以我们选择了这种替代方法。
 * 

 * 
 * 需要决定的是我们坐在两个不同细化水平的单元之间的接口的哪一边。
 * 

 * 
 * 让我们先来看看邻居更精细的情况。然后，我们必须在当前单元格的面的子代上循环，并在每个子代上进行整合。我们在代码中加入了几个断言，以确保我们试图找出邻居的哪个子面与当前单元格的某个子面重合的推理是正确的--有点防御性的编程永远不会有坏处。
 * 

 * 
 * 然后我们调用对面进行整合的函数；由于这是一个内部面，第五个参数是假的，第六个参数被忽略了，所以我们再次传递一个无效的值。
 * 

 * 
 * 
 * @code
 *           else 
 *             { 
 *               if (cell->neighbor(face_no)->has_children()) 
 *                 { 
 *                   const unsigned int neighbor2 = 
 *                     cell->neighbor_of_neighbor(face_no); 
 * 
 *                   for (unsigned int subface_no = 0; 
 *                        subface_no < cell->face(face_no)->n_children(); 
 *                        ++subface_no) 
 *                     { 
 *                       const typename DoFHandler<dim>::active_cell_iterator 
 *                         neighbor_child = 
 *                           cell->neighbor_child_on_subface(face_no, subface_no); 
 * 
 *                       Assert(neighbor_child->face(neighbor2) == 
 *                                cell->face(face_no)->child(subface_no), 
 *                              ExcInternalError()); 
 *                       Assert(neighbor_child->is_active(), ExcInternalError()); 
 * 
 *                       fe_v_subface.reinit(cell, face_no, subface_no); 
 *                       fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 
 * 
 *                       neighbor_child->get_dof_indices(dof_indices_neighbor); 
 * 
 *                       assemble_face_term( 
 *                         face_no, 
 *                         fe_v_subface, 
 *                         fe_v_face_neighbor, 
 *                         dof_indices, 
 *                         dof_indices_neighbor, 
 *                         false, 
 *                         numbers::invalid_unsigned_int, 
 *                         neighbor_child->face(neighbor2)->diameter()); 
 *                     } 
 *                 } 
 * 
 * @endcode
 * 
 * 我们必须关注的另一种可能性是邻居是否比当前单元更粗（特别是，由于每个面只有一个悬挂节点的通常限制，邻居必须正好比当前单元更粗一级，这是我们用断言检查的）。同样，我们在这个接口上进行整合。
 * 

 * 
 * 
 * @code
 *               else if (cell->neighbor(face_no)->level() != cell->level()) 
 *                 { 
 *                   const typename DoFHandler<dim>::cell_iterator neighbor = 
 *                     cell->neighbor(face_no); 
 *                   Assert(neighbor->level() == cell->level() - 1, 
 *                          ExcInternalError()); 
 * 
 *                   neighbor->get_dof_indices(dof_indices_neighbor); 
 * 
 *                   const std::pair<unsigned int, unsigned int> faceno_subfaceno = 
 *                     cell->neighbor_of_coarser_neighbor(face_no); 
 *                   const unsigned int neighbor_face_no = faceno_subfaceno.first, 
 *                                      neighbor_subface_no = 
 *                                        faceno_subfaceno.second; 
 * 
 *                   Assert(neighbor->neighbor_child_on_subface( 
 *                            neighbor_face_no, neighbor_subface_no) == cell, 
 *                          ExcInternalError()); 
 * 
 *                   fe_v_face.reinit(cell, face_no); 
 *                   fe_v_subface_neighbor.reinit(neighbor, 
 *                                                neighbor_face_no, 
 *                                                neighbor_subface_no); 
 * 
 *                   assemble_face_term(face_no, 
 *                                      fe_v_face, 
 *                                      fe_v_subface_neighbor, 
 *                                      dof_indices, 
 *                                      dof_indices_neighbor, 
 *                                      false, 
 *                                      numbers::invalid_unsigned_int, 
 *                                      cell->face(face_no)->diameter()); 
 *                 } 
 *             } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawassemble_cell_term"></a> 
 * <h4>ConservationLaw::assemble_cell_term</h4>
 * 

 * 
 * 这个函数通过计算残差的单元部分来组装单元项，将其负数加到右手边的向量上，并将其相对于局部变量的导数加到雅各布系数（即牛顿矩阵）上。回顾一下，单元格对残差的贡献为 $R_i = \left(\frac{\mathbf{w}^{k}_{n+1} - \mathbf{w}_n}{\delta t} , \mathbf{z}_i \right)_K $ 。
 * $ + \theta \mathbf{B}(\mathbf{w}^{k}_{n+1})(\mathbf{z}_i)_K $  
 * $ + (1-\theta) \mathbf{B}(\mathbf{w}_{n}) (\mathbf{z}_i)_K $ ，其中 $\mathbf{B}(\mathbf{w})(\mathbf{z}_i)_K = - \left(\mathbf{F}(\mathbf{w}),\nabla\mathbf{z}_i\right)_K $  。
 * $ + h^{\eta}(\nabla \mathbf{w} , \nabla \mathbf{z}_i)_K $  
 * $ - (\mathbf{G}(\mathbf {w}), \mathbf{z}_i)_K $ 为 $\mathbf{w} = \mathbf{w}^k_{n+1}$ 和 $\mathbf{w} = \mathbf{w}_{n}$  ， $\mathbf{z}_i$ 为 $i$ 的第1个向量值测试函数。  此外，标量积 $\left(\mathbf{F}(\mathbf{w}), \nabla\mathbf{z}_i\right)_K$ 可以理解为 $\int_K \sum_{c=1}^{\text{n\_components}}  \sum_{d=1}^{\text{dim}} \mathbf{F}(\mathbf{w})_{cd} \frac{\partial z^c_i}{x_d}$ ，其中 $z^c_i$ 是 $i$ 第1个测试函数的 $c$ 分量。
 * 

 * 
 * 在这个函数的顶部，我们做了一些常规的内务工作，即分配一些我们以后需要的局部变量。特别是，我们将分配一些变量来保存 $k$ 次牛顿迭代后的当前解 $W_{n+1}^k$ （变量 <code>W</code> ）和上一时间步长的解 $W_{n}$ （变量 <code>W_old</code> ）的值。
 * 

 * 
 * 除此以外，我们还需要当前变量的梯度。 我们必须计算这些是有点遗憾的，我们几乎不需要。 一个简单的守恒定律的好处是，通量一般不涉及任何梯度。 然而，我们确实需要这些梯度，用于扩散稳定化。
 * 

 * 
 * 我们存储这些变量的实际格式需要一些解释。首先，我们需要解向量的 <code>EulerEquations::n_components</code> 分量在每个正交点的数值。这就构成了一个二维表，我们使用deal.II的表类（这比 <code>std::vector@<std::vector@<T@> @></code> 更有效，因为它只需要分配一次内存，而不是为外向量的每个元素分配一次）。同样地，梯度是一个三维表，Table类也支持。
 * 

 * 
 * 其次，我们想使用自动微分。为此，我们使用 Sacado::Fad::DFad 模板来计算所有我们想计算导数的变量。这包括当前解和正交点的梯度（是自由度的线性组合），以及由它们计算出来的所有东西，如残差，但不包括前一个时间步长的解。这些变量都可以在函数的第一部分找到，同时还有一个变量，我们将用它来存储残差的一个分量的导数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::assemble_cell_term( 
 *     const FEValues<dim> &                       fe_v, 
 *     const std::vector<types::global_dof_index> &dof_indices) 
 *   { 
 *     const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 
 *     const unsigned int n_q_points    = fe_v.n_quadrature_points; 
 * 
 *     Table<2, Sacado::Fad::DFad<double>> W(n_q_points, 
 *                                           EulerEquations<dim>::n_components); 
 * 
 *     Table<2, double> W_old(n_q_points, EulerEquations<dim>::n_components); 
 * 
 *     Table<3, Sacado::Fad::DFad<double>> grad_W( 
 *       n_q_points, EulerEquations<dim>::n_components, dim); 
 * 
 *     Table<3, double> grad_W_old(n_q_points, 
 *                                 EulerEquations<dim>::n_components, 
 *                                 dim); 
 * 
 *     std::vector<double> residual_derivatives(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来，我们必须定义自变量，我们将尝试通过解决一个牛顿步骤来确定自变量。这些自变量是局部自由度的值，我们在这里提取。
 * 

 * 
 * 
 * @code
 *     std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
 *       dofs_per_cell); 
 *     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *       independent_local_dof_values[i] = current_solution(dof_indices[i]); 
 * 
 * @endcode
 * 
 * 下一步包含了所有的魔法：我们宣布自分变量的一个子集为独立自由度，而所有其他的变量仍然是依存函数。这些正是刚刚提取的局部自由度。所有引用它们的计算（无论是直接还是间接）都将积累与这些变量有关的敏感度。
 * 

 * 
 * 为了将这些变量标记为独立变量，下面的方法可以起到作用，将 <code>independent_local_dof_values[i]</code> 标记为总共 <code>dofs_per_cell</code> 中的 $i$ 个独立变量。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *       independent_local_dof_values[i].diff(i, dofs_per_cell); 
 * 
 * @endcode
 * 
 * 在所有这些声明之后，让我们实际计算一些东西。首先， <code>W</code>, <code>W_old</code>, <code>grad_W</code> 和 <code>grad_W_old</code> 的值，我们可以通过使用公式 $W(x_q)=\sum_i \mathbf W_i \Phi_i(x_q)$ 从局部DoF值计算出来，其中 $\mathbf W_i$ 是解向量（局部部分）的第 $i$ 项，而 $\Phi_i(x_q)$ 是在正交点 $x_q$ 评估的第 $i$ 个矢量值的形状函数的值。梯度可以用类似的方法来计算。
 * 

 * 
 * 理想情况下，我们可以通过调用类似 FEValues::get_function_values 和 FEValues::get_function_gradients, 的东西来计算这些信息，但是由于（i）我们必须为此扩展FEValues类，以及（ii）我们不想让整个 <code>old_solution</code> 矢量fad类型，只有局部单元变量，我们明确编码上面的循环。在这之前，我们增加一个循环，将所有的fad变量初始化为零。
 * 

 * 
 * 
 * @code
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       for (unsigned int c = 0; c < EulerEquations<dim>::n_components; ++c) 
 *         { 
 *           W[q][c]     = 0; 
 *           W_old[q][c] = 0; 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             { 
 *               grad_W[q][c][d]     = 0; 
 *               grad_W_old[q][c][d] = 0; 
 *             } 
 *         } 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *         { 
 *           const unsigned int c = 
 *             fe_v.get_fe().system_to_component_index(i).first; 
 * 
 *           W[q][c] += independent_local_dof_values[i] * 
 *                      fe_v.shape_value_component(i, q, c); 
 *           W_old[q][c] += 
 *             old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, c); 
 * 
 *           for (unsigned int d = 0; d < dim; d++) 
 *             { 
 *               grad_W[q][c][d] += independent_local_dof_values[i] * 
 *                                  fe_v.shape_grad_component(i, q, c)[d]; 
 *               grad_W_old[q][c][d] += old_solution(dof_indices[i]) * 
 *                                      fe_v.shape_grad_component(i, q, c)[d]; 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 接下来，为了计算单元贡献，我们需要在所有正交点评估 $\mathbf{F}({\mathbf w}^k_{n+1})$  ,  $\mathbf{G}({\mathbf w}^k_{n+1})$  和  $\mathbf{F}({\mathbf w}_n)$  ,  $\mathbf{G}({\mathbf w}_n)$  。为了存储这些，我们还需要分配一点内存。请注意，我们以自分变量的方式计算通量矩阵和右手边，这样以后就可以很容易地从中计算出雅各布贡献。
 * 

 * 
 * 
 * @code
 *     std::vector<ndarray<Sacado::Fad::DFad<double>, 
 *                         EulerEquations<dim>::n_components, 
 *                         dim>> 
 *       flux(n_q_points); 
 * 
 *     std::vector<ndarray<double, EulerEquations<dim>::n_components, dim>> 
 *       flux_old(n_q_points); 
 * 
 *  
 *       std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
 *       forcing(n_q_points); 
 * 
 *     std::vector<std::array<double, EulerEquations<dim>::n_components>> 
 *       forcing_old(n_q_points); 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         EulerEquations<dim>::compute_flux_matrix(W_old[q], flux_old[q]); 
 *         EulerEquations<dim>::compute_forcing_vector(W_old[q], forcing_old[q]); 
 *         EulerEquations<dim>::compute_flux_matrix(W[q], flux[q]); 
 *         EulerEquations<dim>::compute_forcing_vector(W[q], forcing[q]); 
 *       } 
 * 
 * @endcode
 * 
 * 我们现在已经有了所有的部件，所以进行组装。 我们有一个通过系统组件的外循环，和一个通过正交点的内循环，在那里我们积累了对 $i$ 的残差 $R_i$ 的贡献。这个残差的一般公式在引言和本函数的顶部给出。然而，考虑到  $i$  第三个（矢量值）测试函数  $\mathbf{z}_i$  实际上只有一个非零分量（关于这个主题的更多信息可以在  @ref  矢量值模块中找到），我们可以把它简化一下。它将由下面的变量 <code>component_i</code> 表示。有了这个，残差项可以重新写成
 * @f{eqnarray*}
 * R_i &=&
 * \left(\frac{(\mathbf{w}_{n+1} -
 * \mathbf{w}_n)_{\text{component\_i}}}{\delta
 * t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K
 * \\ &-& \sum_{d=1}^{\text{dim}} \left(  \theta \mathbf{F}
 * ({\mathbf{w}^k_{n+1}})_{\text{component\_i},d} + (1-\theta)
 * \mathbf{F} ({\mathbf{w}_{n}})_{\text{component\_i},d}  ,
 * \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}} {\partial
 * x_d}\right)_K
 * \\ &+& \sum_{d=1}^{\text{dim}} h^{\eta} \left( \theta \frac{\partial
 * (\mathbf{w}^k_{n+1})_{\text{component\_i}}}{\partial x_d} + (1-\theta)
 * \frac{\partial (\mathbf{w}_n)_{\text{component\_i}}}{\partial x_d} ,
 * \frac{\partial (\mathbf{z}_i)_{\text{component\_i}}}{\partial x_d}
 * \right)_K
 * \\ &-& \left( \theta\mathbf{G}({\mathbf{w}^k_n+1} )_{\text{component\_i}}
 * + (1-\theta)\mathbf{G}({\mathbf{w}_n})_{\text{component\_i}} ,
 * (\mathbf{z}_i)_{\text{component\_i}} \right)_K ,
 * @f}
 * ，其中积分可以理解为通过对正交点求和来评估。
 * 

 * 
 * 我们最初对残差的所有贡献进行正向求和，这样我们就不需要对雅各布项进行负数。 然后，当我们对 <code>right_hand_side</code> 矢量进行求和时，我们就否定了这个残差。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
 *       { 
 *         Sacado::Fad::DFad<double> R_i = 0; 
 * 
 *         const unsigned int component_i = 
 *           fe_v.get_fe().system_to_component_index(i).first; 
 * 
 * @endcode
 * 
 * 每一行（i）的残差将被累积到这个fad变量中。 在这一行的装配结束时，我们将查询这个变量的敏感度，并将其加入到雅各布系数中。
 * 

 * 
 * 
 * @code
 *         for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
 *           { 
 *             if (parameters.is_stationary == false) 
 *               R_i += 1.0 / parameters.time_step * 
 *                      (W[point][component_i] - W_old[point][component_i]) * 
 *                      fe_v.shape_value_component(i, point, component_i) * 
 *                      fe_v.JxW(point); 
 * 
 *             for (unsigned int d = 0; d < dim; d++) 
 *               R_i -= 
 *                 (parameters.theta * flux[point][component_i][d] + 
 *                  (1.0 - parameters.theta) * flux_old[point][component_i][d]) * 
 *                 fe_v.shape_grad_component(i, point, component_i)[d] * 
 *                 fe_v.JxW(point); 
 * 
 *             for (unsigned int d = 0; d < dim; d++) 
 *               R_i += 
 *                 1.0 * 
 *                 std::pow(fe_v.get_cell()->diameter(), 
 *                          parameters.diffusion_power) * 
 *                 (parameters.theta * grad_W[point][component_i][d] + 
 *                  (1.0 - parameters.theta) * grad_W_old[point][component_i][d]) * 
 *                 fe_v.shape_grad_component(i, point, component_i)[d] * 
 *                 fe_v.JxW(point); 
 * 
 *             R_i -= 
 *               (parameters.theta * forcing[point][component_i] + 
 *                (1.0 - parameters.theta) * forcing_old[point][component_i]) * 
 *               fe_v.shape_value_component(i, point, component_i) * 
 *               fe_v.JxW(point); 
 *           } 
 * 
 * @endcode
 * 
 * 在循环结束时，我们必须将敏感度加到矩阵上，并从右手边减去残差。Trilinos FAD数据类型让我们可以使用  <code>R_i.fastAccessDx(k)</code>  访问导数，所以我们将数据存储在一个临时数组中。然后，这些关于整行本地道夫的信息被一次性添加到特里诺斯矩阵中（支持我们选择的数据类型）。
 * 

 * 
 * 
 * @code
 *         for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *           residual_derivatives[k] = R_i.fastAccessDx(k); 
 *         system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 
 *         right_hand_side(dof_indices[i]) -= R_i.val(); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawassemble_face_term"></a> 
 * <h4>ConservationLaw::assemble_face_term</h4>
 * 

 * 
 * 在这里，我们做的事情与前面的函数基本相同。在顶部，我们引入自变量。因为如果我们在两个单元格之间的内部面上工作，也会使用当前的函数，所以自变量不仅是当前单元格上的自由度，而且在内部面上的情况下，也是邻近单元格上的自由度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::assemble_face_term( 
 *     const unsigned int                          face_no, 
 *     const FEFaceValuesBase<dim> &               fe_v, 
 *     const FEFaceValuesBase<dim> &               fe_v_neighbor, 
 *     const std::vector<types::global_dof_index> &dof_indices, 
 *     const std::vector<types::global_dof_index> &dof_indices_neighbor, 
 *     const bool                                  external_face, 
 *     const unsigned int                          boundary_id, 
 *     const double                                face_diameter) 
 *   { 
 *     const unsigned int n_q_points    = fe_v.n_quadrature_points; 
 *     const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 
 * 
 *     std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
 *       dofs_per_cell), 
 *       independent_neighbor_dof_values(external_face == false ? dofs_per_cell : 
 *                                                                0); 
 * 
 *     const unsigned int n_independent_variables = 
 *       (external_face == false ? 2 * dofs_per_cell : dofs_per_cell); 
 * 
 *     for (unsigned int i = 0; i < dofs_per_cell; i++) 
 *       { 
 *         independent_local_dof_values[i] = current_solution(dof_indices[i]); 
 *         independent_local_dof_values[i].diff(i, n_independent_variables); 
 *       } 
 * 
 *     if (external_face == false) 
 *       for (unsigned int i = 0; i < dofs_per_cell; i++) 
 *         { 
 *           independent_neighbor_dof_values[i] = 
 *             current_solution(dof_indices_neighbor[i]); 
 *           independent_neighbor_dof_values[i].diff(i + dofs_per_cell, 
 *                                                   n_independent_variables); 
 *         } 
 * 
 * @endcode
 * 
 * 接下来，我们需要定义保守变量  ${\mathbf W}$  在面的这一侧（  $ {\mathbf W}^+$  ）和另一侧（  ${\mathbf W}^-$  ）的值，对于  ${\mathbf W} = {\mathbf W}^k_{n+1}$  和  ${\mathbf W} = {\mathbf W}_n$  。"这一边 "的值可以用与前一个函数完全相同的方式计算，但注意 <code>fe_v</code> 变量现在是FEFaceValues或FESubfaceValues的类型。
 * 

 * 
 * 
 * @code
 *     Table<2, Sacado::Fad::DFad<double>> Wplus( 
 *       n_q_points, EulerEquations<dim>::n_components), 
 *       Wminus(n_q_points, EulerEquations<dim>::n_components); 
 *     Table<2, double> Wplus_old(n_q_points, EulerEquations<dim>::n_components), 
 *       Wminus_old(n_q_points, EulerEquations<dim>::n_components); 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *         { 
 *           const unsigned int component_i = 
 *             fe_v.get_fe().system_to_component_index(i).first; 
 *           Wplus[q][component_i] += 
 *             independent_local_dof_values[i] * 
 *             fe_v.shape_value_component(i, q, component_i); 
 *           Wplus_old[q][component_i] += 
 *             old_solution(dof_indices[i]) * 
 *             fe_v.shape_value_component(i, q, component_i); 
 *         } 
 * 
 * @endcode
 * 
 * 计算 "对立面 "就比较复杂了。如果这是一个内部面，我们可以像上面那样，简单地使用邻居的独立变量来计算它。
 * 

 * 
 * 
 * @code
 *     if (external_face == false) 
 *       { 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const unsigned int component_i = 
 *                 fe_v_neighbor.get_fe().system_to_component_index(i).first; 
 *               Wminus[q][component_i] += 
 *                 independent_neighbor_dof_values[i] * 
 *                 fe_v_neighbor.shape_value_component(i, q, component_i); 
 *               Wminus_old[q][component_i] += 
 *                 old_solution(dof_indices_neighbor[i]) * 
 *                 fe_v_neighbor.shape_value_component(i, q, component_i); 
 *             } 
 *       } 
 * 
 * @endcode
 * 
 * 另一方面，如果这是一个外部边界面，那么 $\mathbf{W}^-$ 的值将是 $\mathbf{W}^+$ 的函数，或者它们将是规定的，这取决于这里施加的边界条件的种类。
 * 

 * 
 * 为了开始评估，让我们确保为这个边界指定的边界ID是我们在参数对象中实际有数据的一个。接下来，我们对不均匀性的函数对象进行评估。 这有点棘手：一个给定的边界可能同时有规定的和隐含的值。 如果一个特定的成分没有被规定，那么这些值就会被评估为零，并在下面被忽略。
 * 

 * 
 * 剩下的部分由一个实际了解欧拉方程边界条件具体内容的函数完成。请注意，由于我们在这里使用的是fad变量，敏感度将被适当地更新，否则这个过程将是非常复杂的。
 * 

 * 
 * 
 * @code
 *     else 
 *       { 
 *         Assert(boundary_id < Parameters::AllParameters<dim>::max_n_boundaries, 
 *                ExcIndexRange(boundary_id, 
 *                              0, 
 *                              Parameters::AllParameters<dim>::max_n_boundaries)); 
 * 
 *         std::vector<Vector<double>> boundary_values( 
 *           n_q_points, Vector<double>(EulerEquations<dim>::n_components)); 
 *         parameters.boundary_conditions[boundary_id].values.vector_value_list( 
 *           fe_v.get_quadrature_points(), boundary_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; q++) 
 *           { 
 *             EulerEquations<dim>::compute_Wminus( 
 *               parameters.boundary_conditions[boundary_id].kind, 
 *               fe_v.normal_vector(q), 
 *               Wplus[q], 
 *               boundary_values[q], 
 *               Wminus[q]); 
 * 
 * @endcode
 * 
 * 这里我们假设边界类型、边界法向量和边界数据值在时间推进中保持不变。
 * 

 * 
 * 
 * @code
 *             EulerEquations<dim>::compute_Wminus( 
 *               parameters.boundary_conditions[boundary_id].kind, 
 *               fe_v.normal_vector(q), 
 *               Wplus_old[q], 
 *               boundary_values[q], 
 *               Wminus_old[q]); 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 现在我们有了 $\mathbf w^+$ 和 $\mathbf w^-$ ，我们可以去计算每个正交点的数值通量函数 $\mathbf H(\mathbf w^+,\mathbf w^-, \mathbf n)$ 。在调用这个函数之前，我们还需要确定Lax-Friedrich的稳定性参数。
 * 

 * 
 * 
 * @code
 *     std::vector< 
 *       std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
 *       normal_fluxes(n_q_points); 
 *     std::vector<std::array<double, EulerEquations<dim>::n_components>> 
 *       normal_fluxes_old(n_q_points); 
 * 
 *     double alpha; 
 * 
 *  
 *       { 
 *         case Parameters::Flux::constant: 
 *           alpha = parameters.stabilization_value; 
 *           break; 
 *         case Parameters::Flux::mesh_dependent: 
 *           alpha = face_diameter / (2.0 * parameters.time_step); 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *           alpha = 1; 
 *       } 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         EulerEquations<dim>::numerical_normal_flux( 
 *           fe_v.normal_vector(q), Wplus[q], Wminus[q], alpha, normal_fluxes[q]); 
 *         EulerEquations<dim>::numerical_normal_flux(fe_v.normal_vector(q), 
 *                                                    Wplus_old[q], 
 *                                                    Wminus_old[q], 
 *                                                    alpha, 
 *                                                    normal_fluxes_old[q]); 
 *       } 
 * 
 * @endcode
 * 
 * 现在以与前面函数中的单元格贡献完全相同的方式组装面项。唯一不同的是，如果这是一个内部面，我们还必须考虑到剩余贡献对相邻单元自由度的敏感性。
 * 

 * 
 * 
 * @code
 *     std::vector<double> residual_derivatives(dofs_per_cell); 
 *     for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
 *       if (fe_v.get_fe().has_support_on_face(i, face_no) == true) 
 *         { 
 *           Sacado::Fad::DFad<double> R_i = 0; 
 * 
 *           for (unsigned int point = 0; point < n_q_points; ++point) 
 *             { 
 *               const unsigned int component_i = 
 *                 fe_v.get_fe().system_to_component_index(i).first; 
 * 
 *               R_i += (parameters.theta * normal_fluxes[point][component_i] + 
 *                       (1.0 - parameters.theta) * 
 *                         normal_fluxes_old[point][component_i]) * 
 *                      fe_v.shape_value_component(i, point, component_i) * 
 *                      fe_v.JxW(point); 
 *             } 
 * 
 *           for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *             residual_derivatives[k] = R_i.fastAccessDx(k); 
 *           system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 
 * 
 *           if (external_face == false) 
 *             { 
 *               for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                 residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell + k); 
 *               system_matrix.add(dof_indices[i], 
 *                                 dof_indices_neighbor, 
 *                                 residual_derivatives); 
 *             } 
 * 
 *           right_hand_side(dof_indices[i]) -= R_i.val(); 
 *         } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawsolve"></a> 
 * <h4>ConservationLaw::solve</h4>
 * 

 * 
 * 在这里，我们实际解决线性系统，使用Trilinos的Aztec或Amesos线性求解器。计算的结果将被写入传递给这个函数的参数向量中。其结果是一对迭代次数和最终的线性残差。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<unsigned int, double> 
 *   ConservationLaw<dim>::solve(Vector<double> &newton_update) 
 *   { 
 *     switch (parameters.solver) 
 *       { 
 * 
 * @endcode
 * 
 * 如果参数文件指定要使用直接求解器，那么我们就到这里。这个过程很简单，因为deal.II在Trilinos中为Amesos直接求解器提供了一个封装类。我们所要做的就是创建一个求解器控制对象（这里只是一个虚拟对象，因为我们不会进行任何迭代），然后创建直接求解器对象。在实际进行求解时，注意我们没有传递一个预处理程序。无论如何，这对直接求解器来说没有什么意义。 最后我们返回求解器的控制统计信息&mdash;它将告诉我们没有进行任何迭代，并且最终的线性残差为零，这里没有任何可能提供的更好的信息。
 * 

 * 
 * 
 * @code
 *         case Parameters::Solver::direct: 
 *           { 
 *             SolverControl                                  solver_control(1, 0); 
 *             TrilinosWrappers::SolverDirect::AdditionalData data( 
 *               parameters.output == Parameters::Solver::verbose); 
 *             TrilinosWrappers::SolverDirect direct(solver_control, data); 
 * 
 *             direct.solve(system_matrix, newton_update, right_hand_side); 
 * 
 *             return {solver_control.last_step(), solver_control.last_value()}; 
 *           } 
 * 
 * @endcode
 * 
 * 同样地，如果我们要使用一个迭代求解器，我们使用Aztec的GMRES求解器。我们也可以在这里使用Trilinos的迭代求解器和预处理类，但是我们选择直接使用Aztec的求解器。对于给定的问题，Aztec的内部预处理实现优于deal.II的包装类，所以我们在AztecOO求解器中使用ILU-T预处理，并设置了一堆可以从参数文件中修改的选项。
 * 

 * 
 * 还有两个实际问题。由于我们将右手边和求解向量建立为deal.II向量对象（而不是矩阵，它是一个Trilinos对象），我们必须将Trilinos Epetra向量交给求解器。 幸运的是，他们支持 "视图 "的概念，所以我们只需发送一个指向deal.II向量的指针。我们必须为设置平行分布的向量提供一个Epetra_Map，这只是一个串行的假对象。最简单的方法是要求矩阵提供它的地图，我们要用它为矩阵-向量乘积做好准备。
 * 

 * 
 * 其次，Aztec求解器希望我们传入一个Trilinos Epetra_CrsMatrix，而不是 deal.II包装类本身。所以我们通过trilinos_matrix()命令来访问Trilinos包装类中的实际Trilinos矩阵。Trilinos希望矩阵是非常量的，所以我们必须使用const_cast手动删除常量。
 * 

 * 
 * 
 * @code
 *         case Parameters::Solver::gmres: 
 *           { 
 *             Epetra_Vector x(View, 
 *                             system_matrix.trilinos_matrix().DomainMap(), 
 *                             newton_update.begin()); 
 *             Epetra_Vector b(View, 
 *                             system_matrix.trilinos_matrix().RangeMap(), 
 *                             right_hand_side.begin()); 
 * 
 *             AztecOO solver; 
 *             solver.SetAztecOption( 
 *               AZ_output, 
 *               (parameters.output == Parameters::Solver::quiet ? AZ_none : 
 *                                                                 AZ_all)); 
 *             solver.SetAztecOption(AZ_solver, AZ_gmres); 
 *             solver.SetRHS(&b); 
 *             solver.SetLHS(&x); 
 * 
 *             solver.SetAztecOption(AZ_precond, AZ_dom_decomp); 
 *             solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
 *             solver.SetAztecOption(AZ_overlap, 0); 
 *             solver.SetAztecOption(AZ_reorder, 0); 
 * 
 *             solver.SetAztecParam(AZ_drop, parameters.ilut_drop); 
 *             solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill); 
 *             solver.SetAztecParam(AZ_athresh, parameters.ilut_atol); 
 *             solver.SetAztecParam(AZ_rthresh, parameters.ilut_rtol); 
 * 
 *             solver.SetUserMatrix( 
 *               const_cast<Epetra_CrsMatrix *>(&system_matrix.trilinos_matrix())); 
 * 
 *             solver.Iterate(parameters.max_iterations, 
 *                            parameters.linear_residual); 
 * 
 *             return {solver.NumIters(), solver.TrueResidual()}; 
 *           } 
 *       } 
 * 
 *     Assert(false, ExcNotImplemented()); 
 *     return {0, 0}; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawcompute_refinement_indicators"></a> 
 * <h4>ConservationLaw::compute_refinement_indicators</h4>
 * 

 * 
 * 这个函数是真正的简单。我们在这里并不假装知道一个好的细化指标会是什么。相反，我们认为 <code>EulerEquation</code> 类会知道这个问题，所以我们只是简单地服从于我们在那里实现的相应函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::compute_refinement_indicators( 
 *     Vector<double> &refinement_indicators) const 
 *   { 
 *     EulerEquations<dim>::compute_refinement_indicators(dof_handler, 
 *                                                        mapping, 
 *                                                        predictor, 
 *                                                        refinement_indicators); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ConservationLawrefine_grid"></a> 
 * <h4>ConservationLaw::refine_grid</h4>
 * 

 * 
 * 在这里，我们使用之前计算的细化指标来细化网格。在开始的时候，我们在所有的单元格上循环，并标记那些我们认为应该被细化的单元格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   ConservationLaw<dim>::refine_grid(const Vector<double> &refinement_indicators) 
 *   { 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         const unsigned int cell_no = cell->active_cell_index(); 
 *         cell->clear_coarsen_flag(); 
 *         cell->clear_refine_flag(); 
 * 
 *         if ((cell->level() < parameters.shock_levels) && 
 *             (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val)) 
 *           cell->set_refine_flag(); 
 *         else if ((cell->level() > 0) && 
 *                  (std::fabs(refinement_indicators(cell_no)) < 
 *                   0.75 * parameters.shock_val)) 
 *           cell->set_coarsen_flag(); 
 *       } 
 * 
 * @endcode
 * 
 * 然后，我们需要在进行细化的同时，将各种解决方案向量从旧网格转移到新网格。SolutionTransfer类是我们的朋友；它有相当丰富的文档，包括例子，所以我们不会对下面的代码做太多评论。最后三行只是把其他一些向量的大小重新设置为现在的正确大小。
 * 

 * 
 * 
 * @code
 *     std::vector<Vector<double>> transfer_in; 
 *     std::vector<Vector<double>> transfer_out; 
 * 
 *     transfer_in.push_back(old_solution); 
 *     transfer_in.push_back(predictor); 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 * 
 *     SolutionTransfer<dim> soltrans(dof_handler); 
 *     soltrans.prepare_for_coarsening_and_refinement(transfer_in); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 * 
 *     dof_handler.clear(); 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     { 
 *       Vector<double> new_old_solution(1); 
 *       Vector<double> new_predictor(1); 
 * 
 *       transfer_out.push_back(new_old_solution); 
 *       transfer_out.push_back(new_predictor); 
 *       transfer_out[0].reinit(dof_handler.n_dofs()); 
 *       transfer_out[1].reinit(dof_handler.n_dofs()); 
 *     } 
 * 
 *     soltrans.interpolate(transfer_in, transfer_out); 
 * 
 *  
 *     old_solution = transfer_out[0]; 
 * 
 *     predictor.reinit(transfer_out[1].size()); 
 *     predictor = transfer_out[1]; 
 * 
 *     current_solution.reinit(dof_handler.n_dofs()); 
 *     current_solution = old_solution; 
 *     right_hand_side.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConservationLawoutput_results"></a> 
 * <h4>ConservationLaw::output_results</h4>
 * 

 * 
 * 现在的这个函数是相当直接的。所有的魔法，包括将数据从保守变量转化为物理变量，都已经被抽象化，并被移到EulerEquations类中，这样在我们想要解决其他双曲守恒定律时就可以被替换。
 * 

 * 
 * 请注意，输出文件的数量是通过保持一个静态变量形式的计数器来确定的，这个计数器在我们第一次来到这个函数时被设置为零，并在每次调用结束时被增加一。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::output_results() const 
 *   { 
 *     typename EulerEquations<dim>::Postprocessor postprocessor( 
 *       parameters.schlieren_plot); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 *     data_out.add_data_vector(current_solution, 
 *                              EulerEquations<dim>::component_names(), 
 *                              DataOut<dim>::type_dof_data, 
 *                              EulerEquations<dim>::component_interpretation()); 
 * 
 *     data_out.add_data_vector(current_solution, postprocessor); 
 * 
 *     data_out.build_patches(); 
 * 
 *     static unsigned int output_file_number = 0; 
 *     std::string         filename = 
 *       "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtk(output); 
 * 
 *     ++output_file_number; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ConservationLawrun"></a> 
 * <h4>ConservationLaw::run</h4>
 * 

 * 
 * 这个函数包含了这个程序的顶层逻辑：初始化，时间循环，以及牛顿内部迭代。
 * 

 * 
 * 在开始时，我们读取参数文件指定的网格文件，设置DoFHandler和各种向量，然后在这个网格上插值给定的初始条件。然后我们在初始条件的基础上进行一系列的网格细化，以获得一个已经很适应起始解的网格。在这个过程结束时，我们输出初始解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConservationLaw<dim>::run() 
 *   { 
 *     { 
 *       GridIn<dim> grid_in; 
 *       grid_in.attach_triangulation(triangulation); 
 * 
 *       std::ifstream input_file(parameters.mesh_filename); 
 *       Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str())); 
 * 
 *       grid_in.read_ucd(input_file); 
 *     } 
 * 
 *     dof_handler.clear(); 
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 所有字段的大小。
 * 

 * 
 * 
 * @code
 *     old_solution.reinit(dof_handler.n_dofs()); 
 *     current_solution.reinit(dof_handler.n_dofs()); 
 *     predictor.reinit(dof_handler.n_dofs()); 
 *     right_hand_side.reinit(dof_handler.n_dofs()); 
 * 
 *     setup_system(); 
 * 
 *     VectorTools::interpolate(dof_handler, 
 *                              parameters.initial_conditions, 
 *                              old_solution); 
 *     current_solution = old_solution; 
 *     predictor        = old_solution; 
 * 
 *     if (parameters.do_refine == true) 
 *       for (unsigned int i = 0; i < parameters.shock_levels; ++i) 
 *         { 
 *           Vector<double> refinement_indicators(triangulation.n_active_cells()); 
 * 
 *           compute_refinement_indicators(refinement_indicators); 
 *           refine_grid(refinement_indicators); 
 * 
 *           setup_system(); 
 * 
 *           VectorTools::interpolate(dof_handler, 
 *                                    parameters.initial_conditions, 
 *                                    old_solution); 
 *           current_solution = old_solution; 
 *           predictor        = old_solution; 
 *         } 
 * 
 *     output_results(); 
 * 
 * @endcode
 * 
 * 然后我们进入主时间步进循环。在顶部，我们简单地输出一些状态信息，这样就可以跟踪计算的位置，以及显示非线性内部迭代进展的表格的标题。
 * 

 * 
 * 
 * @code
 *     Vector<double> newton_update(dof_handler.n_dofs()); 
 * 
 *     double time        = 0; 
 *     double next_output = time + parameters.output_step; 
 * 
 *     predictor = old_solution; 
 *     while (time < parameters.final_time) 
 *       { 
 *         std::cout << "T=" << time << std::endl 
 *                   << "   Number of active cells:       " 
 *                   << triangulation.n_active_cells() << std::endl 
 *                   << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *                   << std::endl 
 *                   << std::endl; 
 * 
 *         std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl 
 *                   << "   _____________________________________" << std::endl; 
 * 
 * @endcode
 * 
 * 然后是内牛顿迭代，解决每个时间步长的非线性问题。它的工作方式是将矩阵和右手边重置为零，然后组装线性系统。如果右手边的规范足够小，那么我们就宣布牛顿迭代已经收敛了。否则，我们求解线性系统，用牛顿增量更新当前解，并输出收敛信息。最后，我们检查牛顿迭代的次数是否超过了10次的限制--如果超过了，就说明迭代有可能出现了发散，继续迭代也没有什么好处。如果发生这种情况，我们就抛出一个异常，这个异常将在 <code>main()</code> 中被捕获，并在程序终止前显示状态信息。
 * 

 * 
 * 注意，我们写AssertThrow宏的方式基本上等同于写<code>if (!(nonlin_iter  @<=  10)) throw ExcMessage ("No convergence in nonlinear solver");</code>这样的话。唯一显著的区别是，AssertThrow还确保被抛出的异常带有它产生的位置（文件名和行号）的信息。这在这里不是太关键，因为只有一个地方可能发生这种异常；然而，当人们想找出错误发生的地方时，它通常是一个非常有用的工具。
 * 

 * 
 * 
 * @code
 *         unsigned int nonlin_iter = 0; 
 *         current_solution         = predictor; 
 *         while (true) 
 *           { 
 *             system_matrix = 0; 
 * 
 *             right_hand_side = 0; 
 *             assemble_system(); 
 * 
 *             const double res_norm = right_hand_side.l2_norm(); 
 *             if (std::fabs(res_norm) < 1e-10) 
 *               { 
 *                 std::printf("   %-16.3e (converged)\n\n", res_norm); 
 *                 break; 
 *               } 
 *             else 
 *               { 
 *                 newton_update = 0; 
 * 
 *                 std::pair<unsigned int, double> convergence = 
 *                   solve(newton_update); 
 * 
 *                 current_solution += newton_update; 
 * 
 *                 std::printf("   %-16.3e %04d        %-5.2e\n", 
 *                             res_norm, 
 *                             convergence.first, 
 *                             convergence.second); 
 *               } 
 * 
 *             ++nonlin_iter; 
 *             AssertThrow(nonlin_iter <= 10, 
 *                         ExcMessage("No convergence in nonlinear solver")); 
 *           } 
 * 
 * @endcode
 * 
 * 只有在牛顿迭代已经收敛的情况下，我们才会到达这一点，所以在这里做各种收敛后的任务。
 * 

 * 
 * 首先，我们更新时间，如果需要的话，产生图形输出。然后，我们通过近似 $\mathbf w^{n+1}\approx \mathbf w^n + \delta t \frac{\partial \mathbf w}{\partial t} \approx \mathbf w^n + \delta t \; \frac{\mathbf w^n-\mathbf w^{n-1}}{\delta t} = 2 \mathbf w^n - \mathbf w^{n-1}$ 来更新下一个时间步长的解决方案的预测器，以尝试使适应性更好地工作。 我们的想法是尝试在前面进行细化，而不是步入一个粗略的元素集并抹去旧的解决方案。 这个简单的时间推断器可以完成这个工作。有了这个，如果用户需要的话，我们就可以对网格进行细化，最后继续进行下一个时间步骤。
 * 

 * 
 * 
 * @code
 *         time += parameters.time_step; 
 * 
 *         if (parameters.output_step < 0) 
 *           output_results(); 
 *         else if (time >= next_output) 
 *           { 
 *             output_results(); 
 *             next_output += parameters.output_step; 
 *           } 
 * 
 *         predictor = current_solution; 
 *         predictor.sadd(2.0, -1.0, old_solution); 
 * 
 *         old_solution = current_solution; 
 * 
 *         if (parameters.do_refine == true) 
 *           { 
 *             Vector<double> refinement_indicators( 
 *               triangulation.n_active_cells()); 
 *             compute_refinement_indicators(refinement_indicators); 
 * 
 *             refine_grid(refinement_indicators); 
 *             setup_system(); 
 * 
 *             newton_update.reinit(dof_handler.n_dofs()); 
 *           } 
 *       } 
 *   } 
 * } // namespace Step33 
 * @endcode
 * 
 * 
 * <a name="main"></a> 
 * <h3>main()</h3>
 * 

 * 
 * 下面的``main''函数与前面的例子类似，不需要进行注释。请注意，如果在命令行上没有给出输入文件名，程序就会中止。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step33; 
 * 
 *       if (argc != 2) 
 *         { 
 *           std::cout << "Usage:" << argv[0] << " input_file" << std::endl; 
 *           std::exit(1); 
 *         } 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization( 
 *         argc, argv, dealii::numbers::invalid_unsigned_int); 
 * 
 *       ConservationLaw<2> cons(argv[1]); 
 *       cons.run(); 
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
 *     }; 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-33/doc/results.dox

<a name="Results"></a>

<a name="Results"></a><h1>Results</h1>


我们用网格 <code>slide.inp</code> （该文件与本程序的源代码在同一目录下）和以下的输入牌（在同一目录下有 <code>input.prm</code> ）运行该问题。

@verbatim
# Listing of Parameters
# ---------------------


# The input grid
set mesh = slide.inp


# Stabilization parameter
set diffusion power = 2.0


# --------------------------------------------------
# Boundary conditions
# We may specify boundary conditions for up to MAX_BD boundaries.
# Your .inp file should have these boundaries designated.
subsection boundary_1
  set no penetration = true # reflective boundary condition
end


subsection boundary_2
  # outflow boundary
  # set w_2 = pressure
  # set w_2 value = 1.5 - y
end


subsection boundary_3
  set no penetration = true # reflective
  # set w_3 = pressure
  # set w_3 value = 1.0
end


subsection boundary_4
  set no penetration = true #reflective
end


# --------------------------------------------------
# Initial Conditions
# We set the initial conditions of the conservative variables.  These lines
# are passed to the expression parsing function.  You should use x,y,z for
# the coordinate variables.


subsection initial condition
  set w_0 value = 0
  set w_1 value = 0
  set w_2 value = 10*(x<-0.7)*(y> 0.3)*(y< 0.45) + (1-(x<-0.7)*(y> 0.3)*(y< 0.45))*1.0
  set w_3 value = (1.5-(1.0*1.0*y))/0.4
end


# --------------------------------------------------
# Time stepping control
subsection time stepping
  set final time = 10.0 # simulation end time
  set time step  = 0.02 # simulation time step
  set theta scheme value = 0.5
end


subsection linear solver
  set output         = quiet
  set method         = gmres
  set ilut fill      = 1.5
  set ilut drop tolerance = 1e-6
  set ilut absolute tolerance = 1e-6
  set ilut relative tolerance = 1.0
end


# --------------------------------------------------
# Output frequency and kind
subsection output
  set step           = 0.01
  set schlieren plot = true
end


# --------------------------------------------------
# Refinement control
subsection refinement
  set refinement = true # none only other option
  set shock value = 1.5
  set shock levels = 1 # how many levels of refinement to allow
end


# --------------------------------------------------
# Flux parameters
subsection flux
 set stab = constant
 #set stab value = 1.0
end
@endverbatim



当我们运行该程序时，我们会得到以下的输出。

@verbatim
...
T=0.14
   Number of active cells:       1807
   Number of degrees of freedom: 7696


   NonLin Res     Lin Iter       Lin Res
   _____________________________________
   7.015e-03        0008        3.39e-13
   2.150e-05        0008        1.56e-15
   2.628e-09        0008        5.09e-20
   5.243e-16        (converged)


T=0.16
   Number of active cells:       1807
   Number of degrees of freedom: 7696


   NonLin Res     Lin Iter       Lin Res
   _____________________________________
   7.145e-03        0008        3.80e-13
   2.548e-05        0008        7.20e-16
   4.063e-09        0008        2.49e-19
   5.970e-16        (converged)


T=0.18
   Number of active cells:       1807
   Number of degrees of freedom: 7696


   NonLin Res     Lin Iter       Lin Res
   _____________________________________
   7.395e-03        0008        6.69e-13
   2.867e-05        0008        1.33e-15
   4.091e-09        0008        3.35e-19
   5.617e-16        (converged)
...
@endverbatim



这个输出报告了牛顿迭代的进度和时间步长。请注意，我们对牛顿迭代的实现确实显示了预期的二次收敛顺序：每一步的非线性残差的规范大致是上一步的规范的平方。这导致了我们在这里可以看到的非常快速的收敛。这种情况一直保持到 $t=1.9$ 时，这时非线性迭代报告缺乏收敛。

@verbatim
...


T=1.88
   Number of active cells:       2119
   Number of degrees of freedom: 9096


   NonLin Res     Lin Iter       Lin Res
   _____________________________________
   2.251e-01        0012        9.78e-12
   5.698e-03        0012        2.04e-13
   3.896e-05        0012        1.48e-15
   3.915e-09        0012        1.94e-19
   8.800e-16        (converged)


T=1.9
   Number of active cells:       2140
   Number of degrees of freedom: 9184


   NonLin Res     Lin Iter       Lin Res
   _____________________________________
   2.320e-01        0013        3.94e-12
   1.235e-01        0016        6.62e-12
   8.494e-02        0016        6.05e-12
   1.199e+01        0026        5.72e-10
   1.198e+03        0002        1.20e+03
   7.030e+03        0001        nan
   7.030e+03        0001        nan
   7.030e+03        0001        nan
   7.030e+03        0001        nan
   7.030e+03        0001        nan
   7.030e+03        0001        nan





----------------------------------------------------
Exception on processing:


--------------------------------------------------------
An error occurred in line <2476> of file <\step-33.cc> in function
    void Step33::ConservationLaw<dim>::run() [with int dim = 2]
The violated condition was:
    nonlin_iter <= 10
The name and call sequence of the exception was:
    ExcMessage ("No convergence in nonlinear solver")
Additional Information:
No convergence in nonlinear solver


--------------------------------------------------------


Aborting!


----------------------------------------------------
@endverbatim



我们可以通过查看解决方案的动画来找出原因和可能的补救措施。

运行这些计算的结果是一堆输出文件，我们可以将其传递给我们选择的可视化程序。当我们把它们整理成一个电影时，过去几个时间步骤的结果看起来是这样的。

 <img src="https://www.dealii.org/images/steps/developer/step-33.oscillation.gif " alt="" height="300"> 

正如我们所看到的，当大质量的流体碰到左下角时，会发生一些振荡，导致迭代的发散。解决这个问题的一个懒办法是添加更多的粘性。如果我们将扩散功率设置为 $\eta = 1.5$ ，而不是 $2.0$ ，模拟就能度过这一危机。那么，结果就会是这样。


 <img src="https://www.dealii.org/images/steps/developer/step-33.slide.ed2.gif " alt="" height="300"> 

沉重的流体在重力作用下被拉下斜坡，在那里与滑雪屋相撞，并被抛向空中！希望每个人都能逃脱。  希望每个人都能逃过一劫!还有，我们可以看到由于人为的粘性，重质和轻质之间的界限很快就模糊了。

我们还可以直观地看到自适应细化网格的演变。

 <img src="https://www.dealii.org/images/steps/developer/step-33.slide.adapt.ed2.gif " alt="" height="300"> 

根据上面讨论的启发式细化方案，自适应性跟随并先于流动模式。





<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Stabilization"></a><h4>Stabilization</h4>


我们选择的数值方案在人工粘度小的时候不是特别稳定，而在人工粘度大的时候则过于扩散。此外，众所周知，还有一些更先进的技术来稳定解决方案，例如流线扩散、最小二乘法稳定条款、熵粘性。




<a name="Betterlinearsolvers"></a><h4>Better linear solvers</h4>


虽然作为非线性求解器的牛顿方法在时间步长足够小的情况下似乎效果很好，但线性求解器是可以改进的。例如，在目前的方案中，只要我们使用迭代求解器，每个牛顿步骤都要重新计算ILU；同样，对于直接求解器，每个步骤都要计算牛顿矩阵的LU分解。这显然是一种浪费：从一个牛顿步骤到另一个牛顿步骤，可能还有不同的时间步骤，牛顿矩阵不会发生根本性的变化：一个牛顿步骤的ILU或稀疏LU分解可能仍然是下一个牛顿或时间步骤的非常好的预处理。因此，避免重新计算将是减少计算时间的一个好办法。

我们可以更进一步：由于接近收敛时，牛顿矩阵只发生一点变化，我们也许可以定义一个准牛顿方案，即在每次牛顿迭代中我们只重新计算残差（即右手边的向量），并重新使用牛顿矩阵。由此产生的方案很可能不是二次收敛的，我们必须期望多做几次非线性迭代；然而，鉴于我们不必每次都花时间建立牛顿矩阵，由此产生的方案很可能更快。




<a name="Cachetheexplicitpartofresidual"></a><h4>Cache the explicit part of residual</h4>


在 ConservationLaw::assemble_cell_term 函数中计算的残差为 $R_i = \left(\frac{\mathbf{w}^{k}_{n+1} - \mathbf{w}_n}{\delta t}
    , \mathbf{z}_i \right)_K  +
      \theta \mathbf{B}({\mathbf{w}^{k}_{n+1}})(\mathbf{z}_i)_K +
      (1-\theta) \mathbf{B}({\mathbf{w}_{n}}) (\mathbf{z}_i)_K $ 这意味着我们在一个牛顿迭代步骤中计算了两次空间残差：一次是关于当前解 $\mathbf{w}^{k}_{n+1}$ ，另一次是关于最后一个时间步长的解 $\mathbf{w}_{n}$ ，在一个时间步长的所有牛顿迭代中保持相同。在牛顿迭代过程中缓存残差 $ \mathbf{B}({\mathbf{w}_{n}}) (\mathbf{z}_i)_K$ 的显式部分将节省大量的人力。




<a name="Otherconservationlaws"></a><h4>Other conservation laws</h4>


最后，作为超越欧拉方程直接求解的一个方向，本程序非常努力地将所有专门针对欧拉方程的实现分离到一个类中（ <code>EulerEquation</code> 类），而将所有专门用于组装矩阵和向量、非线性和线性求解器以及一般顶层逻辑的实现分离到另一个类中（ <code>ConservationLaw</code> 类）。

通过替换该类中通量矩阵和数值通量的定义，以及那里定义的其他各种部分，应该也可以将 <code>ConservationLaw</code> 类应用于其他双曲守恒定律。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-33.cc"
*/
