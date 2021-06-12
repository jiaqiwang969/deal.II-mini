/**
@page step_31 The step-31 tutorial program
This tutorial depends on step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#TheBoussinesqequations">The Boussinesq equations</a>
        <li><a href="#Boundaryandinitialconditions">Boundary and initial conditions</a>
        <li><a href="#Solutionapproach">Solution approach</a>
      <ul>
        <li><a href="#Timestepping">Time stepping</a>
        <li><a href="#WeakformandspacediscretizationfortheStokespart">Weak form and space discretization for the Stokes part</a>
        <li><a href="#Stabilizationweakformandspacediscretizationforthetemperatureequation">Stabilization, weak form and space discretization for the temperature equation</a>
        <li><a href="#Linearsolvers">Linear solvers</a>
      <ul>
        <li><a href="#LinearsolversfortheStokesproblem">Linear solvers for the Stokes problem</a>
        <li><a href="#Linearsolversforthetemperatureequation">Linear solvers for the temperature equation</a>
      </ul>
      </ul>
        <li><a href="#Implementationdetails">Implementation details</a>
      <ul>
        <li><a href="#UsingdifferentDoFHandlerobjects">Using different DoFHandler objects</a>
        <li><a href="#UsingTrilinos">Using Trilinos</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
      <ul>
        <li><a href="#ThecodeInverseMatrixcodeclasstemplate">The <code>InverseMatrix</code> class template</a>
        <li><a href="#Schurcomplementpreconditioner">Schur complement preconditioner</a>
      </ul>
        <li><a href="#ThecodeBoussinesqFlowProblemcodeclasstemplate">The <code>BoussinesqFlowProblem</code> class template</a>
        <li><a href="#BoussinesqFlowProblemclassimplementation">BoussinesqFlowProblem class implementation</a>
      <ul>
        <li><a href="#BoussinesqFlowProblemBoussinesqFlowProblem">BoussinesqFlowProblem::BoussinesqFlowProblem</a>
        <li><a href="#BoussinesqFlowProblemget_maximal_velocity">BoussinesqFlowProblem::get_maximal_velocity</a>
        <li><a href="#BoussinesqFlowProblemget_extrapolated_temperature_range">BoussinesqFlowProblem::get_extrapolated_temperature_range</a>
        <li><a href="#BoussinesqFlowProblemcompute_viscosity">BoussinesqFlowProblem::compute_viscosity</a>
        <li><a href="#BoussinesqFlowProblemsetup_dofs">BoussinesqFlowProblem::setup_dofs</a>
        <li><a href="#BoussinesqFlowProblemassemble_stokes_preconditioner">BoussinesqFlowProblem::assemble_stokes_preconditioner</a>
        <li><a href="#BoussinesqFlowProblembuild_stokes_preconditioner">BoussinesqFlowProblem::build_stokes_preconditioner</a>
        <li><a href="#BoussinesqFlowProblemassemble_stokes_system">BoussinesqFlowProblem::assemble_stokes_system</a>
        <li><a href="#BoussinesqFlowProblemassemble_temperature_matrix">BoussinesqFlowProblem::assemble_temperature_matrix</a>
        <li><a href="#BoussinesqFlowProblemassemble_temperature_system">BoussinesqFlowProblem::assemble_temperature_system</a>
        <li><a href="#BoussinesqFlowProblemsolve">BoussinesqFlowProblem::solve</a>
        <li><a href="#BoussinesqFlowProblemoutput_results">BoussinesqFlowProblem::output_results</a>
        <li><a href="#BoussinesqFlowProblemrefine_mesh">BoussinesqFlowProblem::refine_mesh</a>
        <li><a href="#BoussinesqFlowProblemrun">BoussinesqFlowProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Resultsin2d"> Results in 2d </a>
        <li><a href="#Resultsin3d"> Results in 3d </a>
        <li><a href="#Numericalexperimentstodetermineoptimalparameters"> Numerical experiments to determine optimal parameters </a>
      <ul>
        <li><a href="#Choosingicsubksubiicsubksubiandbeta"> Choosing <i>c<sub>k</sub></i><i>c<sub>k</sub></i> and beta </a>
      <ul>
        <li><a href="#ResultsforQsub1subelements">Results for Q<sub>1</sub> elements</a>
        <li><a href="#ResultsforQsub2subelements">Results for Q<sub>2</sub> elements</a>
        <li><a href="#Resultsfor3d">Results for 3d</a>
        <li><a href="#Conclusions">Conclusions</a>
      </ul>
      </ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-31/doc/intro.dox

 <br> 

<i>This program was contributed by Martin Kronbichler and Wolfgang
Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation or of The
California Institute of Technology.
</i>


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="TheBoussinesqequations"></a><h3>The Boussinesq equations</h3>


这个程序涉及一个有趣的物理问题：如果流体（即液体或气体）遇到由温度差异引起的浮力差异，它是如何表现的？很明显，流体中温度较高（因此较轻）的部分会上升，温度较低（密度较大）的部分会在重力作用下下沉。

在流体运动速度足够慢，以至于惯性效应可以被忽略的情况下，描述这种行为的方程是布西尼斯克方程，其内容如下。

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=&


  -\rho\; \beta \; T\; \mathbf{g},
  \\
  \nabla \cdot {\mathbf u} &=& 0,
  \\
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma.


@f}

这些方程属于矢量值问题的范畴（这个主题的顶层概述可以在 @ref vector_valued 模块中找到）。这里， $\mathbf u$ 是速度场， $p$ 是压力， $T$ 是流体的温度。   $\varepsilon ({\mathbf u}) = \frac 12
[(\nabla{\mathbf u}) + (\nabla {\mathbf u})^T]$ 是速度的对称梯度。可以看出，速度和压力解决了描述不可压缩流体运动的斯托克斯方程，这个方程我们以前在步骤22中考虑过；我们将广泛借鉴在该程序中获得的经验，特别是关于高效线性斯托克斯求解器的经验。

流体运动的强制项是流体的浮力，表示为密度 $\rho$ 、热膨胀系数 $\beta$ 、温度 $T$ 和指向下方的重力矢量 $\mathbf{g}$ 的积。在第32步的介绍中给出了为什么右手边看起来像它的推导）。前两个方程描述了流体如何通过移动对温差做出反应，第三个方程说明了流体运动如何影响温度场：它是一个平流扩散方程，即温度附着在流体颗粒上，并在流场中平流，还有一个额外的扩散（热传导）项。在许多应用中，扩散系数相当小，温度方程实际上是传输的，而不是扩散主导的，因此其特征是双曲而不是椭圆；我们在开发一个稳定的离散化时必须考虑到这一点。

在上述方程中，右侧的 $\gamma$ 项表示热源，可能是一个空间和时间上的变化函数。   $\eta$ 和 $\kappa$ 表示粘度和扩散系数，在本教程程序中我们假定这两个系数为常数。当 $\eta$ 取决于温度时，更普遍的情况是物理应用中的一个重要因素。大多数材料随着温度的升高而变得更加流动（即 $\eta$ 随着 $T$ 的降低而降低）；有时，如在温度接近熔点的岩石矿物的情况下， $\eta$ 可能在典型的温度范围内发生数量级的变化。

我们注意到，上述斯托克斯方程可以通过引入<a target="_top"
href="http://en.wikipedia.org/wiki/Rayleigh_number">Rayleigh
number</a>  $\mathrm{Ra}=\frac{\|\mathbf{g}\| \beta \rho}{\eta \kappa} \delta T L^3$ 来实现非维度化，使用的是典型长度尺度 $L$ 、典型温差 $\delta T$ 、密度 $\rho$ 、热扩散率 $\eta$ 和热导率 $\kappa$  。   $\mathrm{Ra}$ 是一个无尺寸的数字，它描述了由温差引起的浮力变化导致的热传输和热扩散导致的热传输的比率。一个小的瑞利数意味着浮力相对于粘度来说并不强，流体运动 $\mathbf{u}$ 足够慢，因此热扩散 $\kappa\nabla T$ 是主要的热传输项。另一方面，高瑞利数的流体将显示出主导热传导的强烈对流。

对于我们感兴趣的计算热对流的大多数流体，瑞利数是非常大的，通常是 $10^6$ 或更大。从方程的结构中，我们看到这将导致大的压力差和大的速度。因此， $T$ 的对流-扩散方程中的对流项也将非常大，这个方程的精确解将要求我们选择小的时间步长。因此，具有大雷利数的问题很难用数值来解决，其原因与<a
href="http://en.wikipedia.org/wiki/Navier-stokes_equations">Navier-Stokes
equations</a>大时难以解决<a
href="http://en.wikipedia.org/wiki/Reynolds_number">Reynolds number
$\mathrm{Re}$</a>的问题相似。

请注意，大的瑞利数不一定涉及大的绝对速度。例如，地幔中的瑞利数大于 $10^6$  。然而，速度却很小：该材料实际上是固体岩石，但它是如此之热，而且处于压力之下，它可以非常缓慢地流动，每年最多只有几厘米的速度。然而，这可以导致在数百万年的时间尺度上的混合，这个时间尺度比相同数量的热量通过热传导分布要短得多，而且这个时间尺度与影响地球内部和表面结构的演变有关。

 @note 如果你对使用该程序作为你自己实验的基础感兴趣，你也会想看看它在step-32中的延续。此外，step-32后来被发展成更大的开放源代码ASPECT（见https://aspect.geodynamics.org/），它可以解决现实的问题，在试图将step-31变形为可以解决任何你想解决的问题之前，你可能想研究一下它。




<a name="Boundaryandinitialconditions"></a><h3>Boundary and initial conditions</h3>


由于Boussinesq方程是在流体运动的惯性不起作用的假设下推导出来的，所以流场在每个时间段完全由该时间段的浮力差决定，而不是由以前的流场决定。这反映在上面的前两个方程是不包含时间导数的稳态斯托克斯方程的事实。因此，我们不需要速度或压力的初始条件。另一方面，温度场确实满足一个有时间导数的方程，所以我们需要初始条件 $T$  。

至于边界条件：如果 $\kappa>0$ ，那么温度满足一个二阶微分方程，需要边界周围所有时间的边界数据。这些数据可以是规定的边界温度 $T|_{\partial\Omega}=T_b$ （Dirichlet边界条件），也可以是规定的热通量 $\mathbf{n}\cdot\kappa\nabla
T|_{\partial\Omega}=\phi$ ；在这个程序中，我们将使用一个绝缘的边界条件，即规定没有热通量。   $\phi=0$  .

同样地，速度场要求我们提出边界条件。这些条件可以是 $\mathbf{u}=0$ 上的无滑移无通量条件 $\partial\Omega$ ，如果流体粘在边界上，或者无正常通量条件 $\mathbf n \cdot \mathbf
u = 0$ ，如果流体可以沿边界流动但不能穿过边界，或者任何数量的其他物理上合理的条件。在这个程序中，我们将使用无正常通量条件。




<a name="Solutionapproach"></a><h3>Solution approach</h3>


与步骤21中解决的方程一样，我们这里有一个微分代数方程（DAE）系统：就时间变量而言，只有温度方程是微分方程，而 $\mathbf{u}$ 和 $p$ 的斯托克斯系统没有时间导数，因此属于必须在每个时间瞬间保持的那种代数约束。与第21步的主要区别是，那里的代数约束是一个混合拉普拉斯系统，其形式为

@f{eqnarray*}
  \mathbf u + {\mathbf K}\lambda \nabla p &=& 0, \\
  \nabla\cdot \mathbf u &=& f,


@f}

现在我们有一个斯托克斯系统

@f{eqnarray*}


  -\nabla \cdot (2 \eta \varepsilon ({\mathbf u})) + \nabla p &=& f, \\
  \nabla\cdot \mathbf u &=& 0,


@f}

其中 $\nabla \cdot \eta \varepsilon (\cdot)$ 是一个类似于拉普拉斯 $\Delta$ 的算子，适用于一个矢量场。

鉴于与我们在步骤21中所做的相似，我们选择类似的方法可能并不令人惊讶，尽管我们将不得不对微分算子左上角的算子变化进行调整。




<a name="Timestepping"></a><h4>Time stepping</h4>


作为DAE的问题结构允许我们使用与我们在步骤21中已经使用的相同的策略，即我们使用一个时间滞后方案：我们首先解决温度方程（使用外推的速度场），然后将新的温度解插入速度方程的右侧。不过，我们在代码中实现这一方案的方式是从一个稍微不同的角度来看问题。我们首先使用前一个时间步长的温度场来求解速度和压力的斯托克斯方程，这意味着我们得到前一个时间步长的速度。换句话说，我们首先求解时间步长 $n - 1$ 的斯托克斯系统，即

@f{eqnarray*}


  -\nabla \cdot (2\eta \varepsilon ({\mathbf u}^{n-1})) + \nabla p^{n-1} &=&


  -\rho\; \beta \; T^{n-1} \mathbf{g},
  \\
  \nabla \cdot {\mathbf u}^{n-1} &=& 0,


@f}

然后用外推速度场的温度方程到时间  $n$  。

与第21步相比，我们在这里将使用一个高阶时间步进方案，即用（单边）差分商 $\frac{\frac 32 T^{n}-2T^{n-1}+\frac 12 T^{n-2}}{k}$ 取代时间导数 $\frac{\partial T}{\partial t}$ ， $k$ 为时间步长。这就得到了离散化的时间温度方程

@f{eqnarray*}
  \frac 32 T^n


  -
  k\nabla \cdot \kappa \nabla T^n
  &=&
  2 T^{n-1}


  -
  \frac 12 T^{n-2}


  -
  k(2{\mathbf u}^{n-1} - {\mathbf u}^{n-2} ) \cdot \nabla (2T^{n-1}-T^{n-2})
  +
  k\gamma.


@f}

请注意温度方程是如何被半显式解决的：扩散被隐式处理，而平流被显式处理，使用温度和速度的外推法（或前推法），包括刚刚计算的速度  ${\mathbf u}^{n-1}$  。对当前时间水平的正向投影  $n$  是由泰勒扩展得出的，  $T^n
\approx T^{n-1} + k_n \frac{\partial T}{\partial t} \approx T^{n-1} + k_n
\frac{T^{n-1}-T^{n-2}}{k_n} = 2T^{n-1}-T^{n-2}$  。我们需要这个投影来保持BDF-2方案的精度。换句话说，我们在显式右手边使用的温度场是当前温度场的二阶近似值&mdash；不完全是显式时间步进方案，但从特征上看也不会太远。

温度外推的引入将时间步长限制在<a href="http://en.wikipedia.org/wiki/Courant–Friedrichs–Lewy_condition">
Courant-Friedrichs-Lewy (CFL) condition</a>，就像在 @ref step_21 "步骤-21 "中一样。(如果我们隐含地处理平流项，我们就不会有这个稳定条件，因为BDF-2方案是A级稳定的，代价是我们需要在每个时间步长建立一个新的温度矩阵。)我们将在<a href="#Results">results
section</a>中讨论时间步长的确切选择，但目前重要的是，这个CFL条件意味着时间步长 $k$ 可能在不同的时间步长中发生变化，我们必须稍微修改上述公式。如果 $k_n,k_{n-1}$ 是当前和前一个时间步长的时间步长，那么我们使用近似值

@f{align*}{
\frac{\partial T}{\partial t} \approx
 \frac 1{k_n}
 \left(
       \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^{n}


       -
       \frac{k_n+k_{n-1}}{k_{n-1}}T^{n-1}
       +
       \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}
 \right)
 @f}

和

@f{align*}{
T^n \approx
   T^{n-1} + k_n \frac{\partial T}{\partial t}
   \approx
   T^{n-1} + k_n
   \frac{T^{n-1}-T^{n-2}}{k_{n-1}}
   =
   \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2},


@f}

并将上述方程概括如下。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^n


  -
  k_n\nabla \cdot \kappa \nabla T^n
  &=&
  \frac{k_n+k_{n-1}}{k_{n-1}} T^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}


  -
  k_n{\mathbf u}^{*,n} \cdot \nabla T^{*,n}
  +
  k_n\gamma,


@f}



其中 ${(\cdot)}^{*,n} = \left(1+\frac{k_n}{k_{n-1}}\right)(\cdot)^{n-1} -
\frac{k_n}{k_{n-1}}(\cdot)^{n-2}$ 表示速度 $\mathbf u$ 和温度 $T$ 外推到时间级别 $n$ ，使用前两个时间步骤的数值。这不是一个容易读懂的方程，但会为我们提供所需的高阶精度。作为一致性检查，很容易验证，如果 $k_n=k_{n-1}$  ，它可以还原成与上面相同的方程。

最后我们注意到，选择高阶时间步进方案当然会迫使我们在内存中保留更多的时间步进；特别是，我们在这里需要保留 $T^{n-2}$ ，这是一个我们以前可以抛弃的向量。这似乎是一个麻烦，我们以前可以通过使用一阶时间步进方案来避免，但是正如我们在下面讨论稳定化问题时看到的那样，我们无论如何都需要这个向量，因此在时间离散化中保留它基本上是免费的，并给我们提供了使用高阶方案的机会。




<a name="WeakformandspacediscretizationfortheStokespart"></a><h4>Weak form and space discretization for the Stokes part</h4>


像解决混合拉普拉斯方程一样，解决斯托克斯方程需要我们为速度和压力变量选择特定的有限元对。因为这在步骤22中已经讨论过了，所以我们只简单介绍一下这个话题。这里，我们使用稳定对 $Q_{p+1}^d \times Q_p, p\ge 1$  。这些都是连续元素，所以我们可以通过部分积分和用离散函数替代连续函数来形成斯托克斯方程的弱形式，没有问题。

@f{eqnarray*}
  (\nabla {\mathbf v}_h, 2\eta \varepsilon ({\mathbf u}^{n-1}_h))


  -
  (\nabla \cdot {\mathbf v}_h, p^{n-1}_h)
  &=&


  -({\mathbf v}_h, \rho\; \beta \; T^{n-1}_h \mathbf{g}),
  \\
  (q_h, \nabla \cdot {\mathbf u}^{n-1}_h) &=& 0,


@f}

为所有测试函数  $\mathbf v_h, q_h$  。第一个方程的第一项被认为是张量之间的内积，即 $(\nabla {\mathbf v}_h, \eta \varepsilon ({\mathbf u}^{n-1}_h))_\Omega
 = \int_\Omega \sum_{i,j=1}^d [\nabla {\mathbf v}_h]_{ij}
           \eta [\varepsilon ({\mathbf u}^{n-1}_h)]_{ij}\, dx$  。因为这个乘积中的第二个张量是对称的，所以 $\nabla {\mathbf v}_h$ 的反对称分量不起作用，如果我们用 $\mathbf v_h$ 的对称梯度代替，会导致完全一样的形式。因此，我们考虑并实施的表述是

@f{eqnarray*}
  (\varepsilon({\mathbf v}_h), 2\eta \varepsilon ({\mathbf u}^{n-1}_h))


  -
  (\nabla \cdot {\mathbf v}_h, p^{n-1}_h)
  &=&


  -({\mathbf v}_h, \rho\; \beta \; T^{n-1}_h \mathbf{g}),
  \\
  (q_h, \nabla \cdot {\mathbf u}^{n-1}_h) &=& 0.


@f}



这与我们在第22步中已经讨论过的完全一样，这里就不多说了。




<a name="Stabilizationweakformandspacediscretizationforthetemperatureequation"></a><h4>Stabilization, weak form and space discretization for the temperature equation</h4>


更有趣的问题是如何处理温度平流-扩散方程。默认情况下，并不是所有这个方程的离散化都是同样稳定的，除非我们要么做一些像上卷、稳定化，或者所有这些的事情。实现这一点的方法之一是使用不连续元素（即我们在步骤12中离散传输方程或在步骤20和步骤21中离散压力时使用的FE_DGQ类），并在单元间的界面上定义一个考虑到上卷的流量。如果我们有一个纯粹的平流问题，这可能是最简单的方法。然而，这里我们也有一些扩散，用不连续元素对拉普拉斯算子进行离散化是很麻烦的，因为有大量的附加项需要在单元间的每个面上进行积分。不连续元素还有一个缺点，即使用数值通量会带来额外的数值扩散，这种扩散无处不在，而我们真的希望将数值扩散的影响降到最低，只在需要稳定方案的地方应用它。

因此，一个更好的选择是在模型中加入一些非线性粘度。从本质上讲，这样做的目的是将温度方程的形式从

@f{eqnarray*}
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T &=& \gamma


@f}

到类似于

@f{eqnarray*}
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot (\kappa+\nu(T)) \nabla T &=& \gamma,


@f}

其中 $\nu(T)$ 是一个额外的粘度（扩散）项，只在冲击和其他不连续点附近发挥作用。   $\nu(T)$ 的选择方式是，如果 $T$ 满足原始方程，则额外的粘性为零。

为了实现这一点，文献中包含了许多方法。我们在这里将遵循Guermond和Popov开发的一种方法，它建立在一个适当定义的残差和一个额外粘度的极限程序之上。为此，让我们定义一个残差 $R_\alpha(T)$ 如下。

@f{eqnarray*}
  R_\alpha(T)
  =
  \left(
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T - \gamma
  \right)
  T^{\alpha-1}


@f}

其中，我们以后将从 $[1,2]$ 范围内选择稳定指数 $\alpha$ 。请注意，如果 $T$ 满足温度方程， $R_\alpha(T)$ 将为零，因为此时括号内的项将为零。将条款相乘，我们得到以下完全等同的形式。

@f{eqnarray*}
  R_\alpha(T)
  =
  \frac 1\alpha
  \frac{\partial (T^\alpha)}{\partial t}
  +
  \frac 1\alpha
  {\mathbf u} \cdot \nabla (T^\alpha)


  -
  \frac 1\alpha
  \nabla \cdot \kappa \nabla (T^\alpha)
  +
  \kappa(\alpha-1)
  T^{\alpha-2} |\nabla T|^2


  -
  \gamma
  T^{\alpha-1}


@f}



有了这个残差，我们现在可以把人工黏度定义为一个片状常数函数，在直径为 $K$ 的每个单元上分别定义如下。

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
  \min\left\{
    h_K,
    h_K^\alpha
    \frac{\|R_\alpha(T)\|_{L^\infty(K)}}{c(\mathbf{u},T)}
  \right\}


@f}



这里， $\beta$ 是一个稳定常数（通过维度分析发现它是无单位的，因此与比例无关；我们将在<a href="#Results">results section</a>中讨论其选择）， $c(\mathbf{u},T)$ 是一个归一化常数，其单位必须是 $\frac{m^{\alpha-1}K^\alpha}{s}$  。我们将选择它作为 $c(\mathbf{u},T) =
 c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
 \ |\mathrm{diam}(\Omega)|^{\alpha-2}$  ，其中 $\mathrm{var}(T)=\max_\Omega T - \min_\Omega T$ 是目前温度值的范围（记住，浮力是由温度变化驱动的，而不是绝对温度）， $c_R$ 是一个无尺寸常数。为了理解这个方法为什么有效，请考虑这个问题。如果在一个特定的单元 $K$ 上，温度场是平滑的，那么我们希望那里的残差很小（事实上是在 ${\cal O}(h_K)$ 的数量级上），注入人工扩散的稳定项在那里的大小将是 $h_K^{\alpha+1}$ &mdash；也就是说，相当小，就像我们希望它在没有必要进行额外扩散时那样。另一方面，如果我们处于或接近温度场的不连续性，那么残差将很大； $\nu_\alpha(T)$ 定义中的最小操作将确保稳定项的大小为 $h_K$ &mdash；这是确保方案稳定的最佳人工粘性量。

这种方案是否真的有效是个好问题。Guermond和Popov的计算表明，这种形式的稳定方案实际上比其他大多数稳定方案（例如流线扩散，仅举最简单的一种）表现得更好。此外，对于 $\alpha\in
[1,2)$ ，他们甚至可以证明，对于线性传输方程，它比流线扩散产生更好的收敛阶数。对于 $\alpha=2$ ，目前还没有理论结果，但数值测试表明，其结果比 $\alpha=1$ 好得多。

一个更实际的问题是如何将这种人工扩散引入我们想要解决的方程。请注意，数值粘度 $\nu(T)$ 是随温度变化的，所以我们要解决的方程在 $T$ 中是非线性的&mdash；这不是人们对稳定方程的简单方法的期望，如果我们意识到 $\nu(T)$ 在 $T$ 中是不可分的，那就更不可能了。然而，我们没有理由绝望：我们仍然要在时间上进行离散，我们可以明确地处理这个术语。

在稳定参数的定义中，我们用  $\frac{\partial T}{\partial t} \approx
\frac{T^{n-1}-T^{n-2}}{k^{n-1}}$  对时间导数进行近似。这种近似只利用了可用的时间数据，这就是我们需要存储前两个时间步骤的数据的原因（这使我们能够使用BDF-2方案而不需要额外的存储成本）。我们现在可以简单地在 $t_{n-1}$ 处评估其余的项，但这样一来，离散残差无非是一个向后的欧拉近似，它只有一阶精度。因此，在平滑解的情况下，尽管外部BDF-2方案和空间FE离散化的时间精度为二阶，但残差仍为 $h$ 阶。这当然不是我们想要的（事实上，我们希望在解决方案表现良好的区域有较小的残差），所以需要更谨慎一些。这个问题的关键是观察我们构造的第一导数实际上是以 $t_{n-\frac{3}{2}}$ 为中心的。如果我们通过使用近似值 $\frac 12 T^{n-1}+\frac 12 T^{n-2}$ 来评估 $t_{n-\frac{3}{2}}$ 处的所有空间项，我们就可以得到所需的二阶精确残差计算，这意味着我们将非线性粘度计算为这个中间温度的函数， $\nu_\alpha =
\nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)$  。请注意，这种对残差的评估无非是一个Crank-Nicholson方案，所以我们可以肯定，现在一切正常了。人们可能会想，现在的数值粘度没有在时间 $n$ 进行评估（相对于方程的其余部分），这是否是一个问题。然而，这种偏移是不严谨的。对于平滑解， $\nu_\alpha$ 将连续变化，所以时间偏移的误差比非线性粘度本身要小 $k$ 倍，也就是说，它是被遗漏的一个小的高阶贡献。这很好，因为该项本身已经达到了光滑区域的离散化误差水平。

使用上面介绍的BDF-2方案，这就得到了更简单的大小为 $k$ 的均匀时间步长的情况。

@f{eqnarray*}
  \frac 32 T^n


  -
  k\nabla \cdot \kappa \nabla T^n
  &=&
  2 T^{n-1}


  -
  \frac 12 T^{n-2}
  \\
  &&
  +
  k\nabla \cdot
  \left[
    \nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)
    \ \nabla (2T^{n-1}-T^{n-2})
  \right]
  \\
  &&


  -
  k(2{\mathbf u}^{n-1}-{\mathbf u}^{n-2}) \cdot \nabla (2T^{n-1}-T^{n-2})
  \\
  &&
  +
  k\gamma.


@f}

在这个方程的左侧仍然是来自时间导数的项和我们隐含处理的原始（物理）扩散（这实际上是一个很好的项：从左侧产生的矩阵是质量矩阵和拉普拉斯矩阵的倍数&mdash；两者都是正定的，如果时间步长 $k$ 很小，和很容易反转）。在右侧，第一行的条款是时间导数的结果；第二行是时间 $t_{n-\frac
32}$ 的人工扩散；第三行包含平流条款，第四行是来源。请注意，人工扩散对当前时间的外推温度的作用，与我们在时间步进一节中讨论的平流作用相同。

我们在现实中必须使用的非均匀时间步长的形式要复杂一些（这就是为什么我们先展示了上面的简单形式），其内容为：。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} T^n


  -
  k_n\nabla \cdot \kappa \nabla T^n
  &=&
  \frac{k_n+k_{n-1}}{k_{n-1}} T^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T^{n-2}
  \\
  &&
  +
  k_n\nabla \cdot
  \left[
    \nu_\alpha\left(\frac 12 T^{n-1}+\frac 12 T^{n-2}\right)
    \ \nabla  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \right]
  \\
  &&


  -
  k_n
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right){\mathbf u}^{n-1} -
    \frac{k_n}{k_{n-1}}{\mathbf u}^{n-2}
  \right]
  \cdot \nabla
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \\
  &&
  +
  k_n\gamma.


@f}



在解决了所有这些问题之后，弱形式自然而然地从最后一个方程中显示的强形式中产生，我们立即得出了离散化方程的弱形式。

@f{eqnarray*}
  \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} (\tau_h,T_h^n)
  +
  k_n (\nabla \tau_h, \kappa \nabla T_h^n)
  &=&
  \biggl(\tau_h,
  \frac{k_n+k_{n-1}}{k_{n-1}} T_h^{n-1}


  -
  \frac{k_n^2}{k_{n-1}(k_n+k_{n-1})} T_h^{n-2}
  \\
  &&\qquad


  -
  k_n
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right){\mathbf u}^{n-1} -
    \frac{k_n}{k_{n-1}}{\mathbf u}^{n-2}
  \right]
  \cdot \nabla
  \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  +
  k_n\gamma \biggr)
  \\
  &&


  -
  k_n \left(\nabla \tau_h,
    \nu_\alpha\left(\frac 12 T_h^{n-1}+\frac 12 T_h^{n-2}\right)
    \ \nabla \left[
    \left(1+\frac{k_n}{k_{n-1}}\right)T^{n-1}-\frac{k_n}{k_{n-1}}T^{n-2}
  \right]
  \right)


@f}

为所有离散测试函数  $\tau_h$  。在这里，扩散项已经被部分整合，我们已经使用，我们将施加没有热通量，  $\mathbf{n}\cdot\kappa\nabla T|_{\partial\Omega}=0$  。

这就产生了一个矩阵方程，其形式为

@f{eqnarray*}
  \left( \frac{2k_n+k_{n-1}}{k_n+k_{n-1}} M+k_n A_T\right) T_h^n
  = F(U_h^{n-1}, U_h^{n-2},T_h^{n-1},T_h^{n-2}),


@f}

考虑到左边的矩阵结构（两个正定矩阵之和），使用共轭梯度法很容易解决这个问题。




<a name="Linearsolvers"></a><h4>Linear solvers</h4>


如上所述，我们解决速度/压力和温度的联合系统的方法是使用算子分割，我们首先用旧的温度场解决速度和压力的斯托克斯系统，然后用刚刚计算的速度场解决新的温度场。关于算子分割方法的更广泛的讨论可以在步骤58中找到）。




<a name="LinearsolversfortheStokesproblem"></a><h5>Linear solvers for the Stokes problem</h5>


解决来自斯托克斯系统的线性方程已经在步骤22中进行了详细的讨论。特别是在该程序的结果部分，我们讨论了一些替代的线性求解器策略，结果发现这些策略比原来的方法更有效。在那里确定的最佳替代方案是使用一个由涉及舒尔补码的块状矩阵预处理的GMRES求解器。具体来说，斯托克斯算子导致了一个块状结构的矩阵

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)


@f}

正如那里所讨论的，一个好的预处理程序是

@f{eqnarray*}
  P
  =
  \left(\begin{array}{cc}
    A & 0 \\ B & -S
  \end{array}\right),
  \qquad
  \text{or equivalently}
  \qquad
  P^{-1}
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
  \end{array}\right)


@f}

其中 $S$ 是斯托克斯算子的舒尔补 $S=B^TA^{-1}B$  。当然，这个预处理程序是没有用的，因为我们不能形成矩阵的各种倒数，但我们可以用下面的方法作为预处理程序。

@f{eqnarray*}
  \tilde P^{-1}
  =
  \left(\begin{array}{cc}
    \tilde A^{-1} & 0 \\ \tilde S^{-1} B \tilde A^{-1} & -\tilde S^{-1}
  \end{array}\right)


@f}

其中 $\tilde A^{-1},\tilde S^{-1}$ 是反矩阵的近似值。特别是，事实证明 $S$ 在光谱上等同于质量矩阵，因此，用适用于压力空间上的质量矩阵的CG求解器取代 $\tilde
S^{-1}$ 是一个不错的选择。与步骤22稍有不同的是，我们在这里的动量方程中有一个系数 $\eta$ ，通过与那里相同的推导，我们应该得出结论，我们应该使用的是具有条目 $\tilde S_{ij}=(\eta^{-1}\varphi_i,\varphi_j)$ 的加权质量矩阵。

想出一个好的替代方案 $\tilde
A^{-1}$ 更为复杂，它对应于矢量值速度场的离散化对称拉普拉斯，即 $A_{ij} = (\varepsilon {\mathbf v}_i, 2\eta \varepsilon ({\mathbf
v}_j))$  。在步骤22中，我们用 $A$ 的稀疏LU分解（使用SparseDirectUMFPACK类）来代替 $\tilde A^{-1}$ &mdash; 完美的前置条件&mdash; 在2D中，但对于3D来说，内存和计算时间通常不足以实际计算这个分解；因此，我们在3D中只使用不完全LU分解（ILU，使用稀疏ILU类）。

对于这个项目，我们想走得更远一点。为此，请注意，矢量场上的对称化双线性形式 $(\varepsilon {\mathbf v}_i, 2 \eta \varepsilon ({\mathbf v}_j))$ 与非对称化版本 $(\nabla {\mathbf v}_i, \eta \nabla {\mathbf v}_j)
= \sum_{k,l=1}^d
  (\partial_k ({\mathbf v}_i)_l, \eta \partial_k ({\mathbf v}_j)_l)
$ 相差不大（请注意，在这个形式中因子2已经消失了）。然而，后者的优点是测试函数的 <code>dim</code> 矢量分量不是耦合的（好吧，几乎是，见下文），也就是说，得到的矩阵是块对角线的：每个矢量分量有一个块，这些块中的每个都等于这个矢量分量的拉普拉斯矩阵。因此，假设我们以这样的方式排列自由度，即首先对速度的所有 $x$ 分量进行编号，然后是 $y$ 分量，然后是 $z$ 分量，那么与这种稍有不同的双线性形式相关的矩阵 $\hat A$ 具有如下形式

@f{eqnarray*}
  \hat A =
  \left(\begin{array}{ccc}
    A_s & 0 & 0 \\ 0 & A_s & 0 \\ 0 & 0 & A_s
  \end{array}\right)


@f}

其中 $A_s$ 是一个拉普拉斯矩阵，其大小等于与矢量值速度的每个分量相关的形状函数数量。有了这个矩阵，我们就可以对速度矩阵 $A$ 的预处理进行如下定义。

@f{eqnarray*}
  \tilde A^{-1} =
  \left(\begin{array}{ccc}
    \tilde A_s^{-1} & 0 & 0 \\
    0 & \tilde A_s^{-1} & 0 \\
    0 & 0 & \tilde A_s^{-1}
  \end{array}\right),


@f}

其中 $\tilde A_s^{-1}$ 是拉普拉斯矩阵的预处理程序&mdash;我们非常清楚如何建立良好的预处理程序!

在现实中，故事并不那么简单。为了使矩阵 $\tilde A$ 确定，我们需要通过应用边界条件使各个块 $\tilde
A_s$ 确定。我们可以尝试通过在边界周围应用狄氏边界条件来做到这一点，然后，如果后者的矩阵是由斯托克斯问题产生的，我们在领域周围的速度分量上也有狄氏边界条件，即如果我们执行 $\mathbf{u} =
0$ ，那么如此定义的前置条件 $\tilde A^{-1}$ 就变成了 $A$ 的良好前置条件。

不幸的是，这个 "如果 "是 "如果且仅是如果"：在下面的程序中，我们将希望使用 $\mathbf u
\cdot \mathbf n = 0$ 形式的无流量边界条件（即允许与边界平行的流量%，但没有通过边界的流量）。在这种情况下，事实证明，上面定义的块状对角线矩阵不是一个好的预处理程序，因为它忽略了边界上的成分耦合。因此，更好的方法是如果我们将矩阵 $\hat A$ 建立为矢量拉普拉斯矩阵 $\hat A_{ij} = (\nabla {\mathbf v}_i,
\eta \nabla {\mathbf v}_j)$ ，然后应用与我们应用于 $A$ 相同的边界条件。如果这是一个围绕域的迪里希特边界条件， $\hat A$ 将像上面那样解耦为三个对角线块，如果边界条件是 $\mathbf u
\cdot \mathbf n = 0$ 的形式，那么这将在边界引入自由度的耦合，但只在那里。事实上，这被证明是一个比上面介绍的更好的预处理程序，而且几乎具有我们希望得到的所有好处。


总结这整个故事，我们可以看到。   <ul>   <li>  与我们在步骤22中从对称梯度产生的原始矩阵 $A$ 建立一个预处理程序相比，我们不得不期待基于拉普拉斯双线性形式的预处理程序表现得更差，因为它没有考虑到向量分量之间的耦合。

    <li> 另一方面，拉普拉斯矩阵的预处理程序通常比矢量问题的预处理程序更成熟，性能更好。例如，在写这篇文章的时候，代数%多重网格（AMG）算法对于标量问题已经非常成熟，但对于矢量问题却不是如此。

    <li> 在建立这个预处理程序时，我们将不得不建立矩阵 $\hat A$ 及其预处理程序。虽然这意味着我们必须存储一个之前不需要的额外矩阵，但与存储耦合矩阵 $A$ 的预处理程序相比，预处理程序 $\tilde A_s^{-1}$ 可能需要的内存要少得多。这是因为矩阵 $A_s$ 每行只有三分之一的条目对应于内部自由度，并且只在边界条件引入耦合的部分包含向量分量之间的耦合。因此，存储该矩阵是比较便宜的，我们可以预期，计算和存储预处理程序 $\tilde A_s$ 也将比为完全耦合的矩阵做这些事情便宜得多。   </ul> 




<a name="Linearsolversforthetemperatureequation"></a><h5>Linear solvers for the temperature equation</h5>


这是最容易的部分。温度方程的矩阵具有 $\alpha M + \beta A$ 的形式，其中 $M,A$ 是温度空间上的质量和刚度矩阵， $\alpha,\beta$ 是与时间步进方案以及当前和前一个时间步进有关的常数。这是一个对称正定和一个对称正半定矩阵之和，其结果也是对称正定的。此外， $\frac\beta\alpha$ 是一个与时间步长成正比的数字，因此只要网格很细就会变小，从而阻尼当时条件不好的刚度矩阵的影响。

因此，用共轭梯度算法反转这个矩阵，使用一个简单的预处理程序，与反转斯托克斯矩阵相比是微不足道和非常便宜的。




<a name="Implementationdetails"></a><h3>Implementation details</h3>


<a name="UsingdifferentDoFHandlerobjects"></a><h4>Using different DoFHandler objects</h4>


关于下面的程序，值得事先解释的一件事是使用了两个不同的DoFHandler对象。如果看一下上述方程的结构和它们的求解方案，就会发现几乎没有什么共同点能使斯托克斯部分和温度部分保持一致。在我们以前讨论 @ref
vector_valued "矢量值问题 "的所有教程程序中，我们总是只使用一个具有几个矢量分量的单一有限元，以及一个DoFHandler对象。有时，我们将得到的矩阵分解成若干块，以方便特定的求解器方案；例如，在目前程序所依据的斯托克斯方程的第22步程序中就是如此。

当然，我们在这里也可以这样做。我们将得到的线性系统看起来像这样。

@f{eqnarray*}
  \left(\begin{array}{ccc}
    A & B^T & 0 \\ B & 0 &0 \\ C & 0 & K
  \end{array}\right)
  \left(\begin{array}{ccc}
    U^{n-1} \\ P^{n-1} \\ T^n
  \end{array}\right)
  =
  \left(\begin{array}{ccc}
    F_U(T^{n-1}) \\ 0 \\ F_T(U^{n-1},U^{n-2},T^{n-1},T^{n-2})
  \end{array}\right).


@f}

这方面的问题是。我们从未同时使用整个矩阵。事实上，它从未真正同时存在。如上所述， $K$ 和 $F_T$ 依赖于已经计算出的解 $U^n$ ，在第一种情况下，通过时间步长（这依赖于 $U^n$ ，因为它必须满足CFL条件）。所以我们只有在已经解决了左上角 $2\times 2$ 块斯托克斯系统后才能组装它，一旦我们转向温度方程，我们就不再需要斯托克斯部分了；我们为一个在任何时候都不会以整体存在于内存中的矩阵建立一个对象，这导致我们在步骤21中跳了一些圈套，所以我们不要重复这类错误。此外，我们实际上并没有建立矩阵 $C$ ：因为当我们进入温度方程时，我们已经知道了 $U^n$ ，而且因为我们必须在这个时候组装右手边的 $F_T$ ，我们只是将项 $CU^n$ 移到右手边，并将其与所有其他项组装在一起。这意味着矩阵中不存在温度变量和斯托克斯变量耦合的部分，因此所有自由度的全局列举不再重要：如果我们有所有斯托克斯自由度的列举，以及所有温度自由度的独立列举就足够了。

从本质上讲，将<i>everything</i>放入一个块状矩阵中并没有什么用处（当然，对于 $2\times 2$ 斯托克斯部分，也有同样好的理由这样做），或者，就这一点而言，将所有东西放入同一个DoFHandler对象。

但这样做是否有<i>downsides</i>的好处？这些问题是存在的，尽管它们一开始可能并不明显。主要问题是，如果我们需要创建一个包含速度、压力和温度形状函数的全局有限元，并使用它来初始化DoFHandler。但是我们也用这个有限元对象来初始化我们使用的所有FEValues或FEFaceValues对象。这可能看起来不是什么大问题，但是想象一下，例如，当我们评估我们需要计算人工粘度 $
  R_\alpha(T)
  =
  \left(
  \frac{\partial T}{\partial t}
  +
  {\mathbf u} \cdot \nabla T


  -
  \nabla \cdot \kappa \nabla T - \gamma
  \right)
  T^{\alpha-1}
$ 的残差 $\nu_\alpha(T)|_K$ 时会发生什么。  为此，我们需要温度的拉普拉斯，我们使用形状函数的二阶导数（Hessians）张量来计算（为此我们必须给FEValues对象加上 <code>update_hessians</code> 标志）。现在，如果我们有一个包含速度、压力和温度的形状函数的有限性，这意味着我们必须计算<i>all</i>形状函数的Hessians，包括速度的许多高阶形状函数。这是很多我们不需要的计算，事实上，如果一个人要这样做（就像我们在程序的早期版本中那样），组装右手边需要大约四分之一的整体计算时间。

所以我们要做的是使用两个不同的有限元对象，一个用于斯托克斯成分，一个用于温度。这样就有两个不同的DoFHandlers，两个稀疏模式和两个用于斯托克斯和温度部分的矩阵，等等。每当我们要组装包含温度和斯托克斯形状函数的东西时（特别是斯托克斯和温度方程的右侧），我们就使用两个FEValues对象，用两个单元格迭代器进行初始化，通过与同一三角化对象相关的两个DoFHandler对象进行平行行走。对于这两个FEValues对象，我们当然使用相同的正交对象，这样我们就可以在同一组正交点上进行迭代，但是每个FEValues对象将只根据它实际需要计算的内容来获得更新标志。特别是，当我们像上面那样计算残差时，我们只要求得到斯托克斯形状函数的值，但也要求得到温度形状函数的Hessians &mdash；确实便宜得多，而且事实证明：组装温度方程的右手边现在是程序中几乎无法测量的一个组成部分。

有了这些变化，对程序进行计时，可以得出只有以下操作与整个运行时间有关。   <ul>   <li>  解决斯托克斯系统：72%的运行时间。     <li>  组装斯托克斯预处理程序，并使用Trilinos ML包计算代数多网格层次结构：占运行时间的11%。     <li>  函数  <code>BoussinesqFlowProblem::setup_dofs</code>  : 占整体运行时间的7%。     <li>  组装斯托克斯和温度右侧向量以及组装矩阵。7%.   </ul>  实质上这意味着除了代数多重网格之外，所有的瓶颈都已经被移除。




<a name="UsingTrilinos"></a><h4>Using Trilinos</h4>


与我们在第17步和第18步中使用PETSc来支持我们的线性代数需求一样，我们在这个程序中使用了<a
href="http://trilinos.org">Trilinos</a>库的接口（安装说明见deal.II README文件）。Trilinos是一个非常大的集合，包括与线性和非线性代数有关的所有东西，以及围绕这些东西的各种工具（看起来它在未来也会向许多其他方向发展）。

使用Trilinos的主要原因，类似于我们探索的PETSc，是它是一个非常强大的库，比deal.II自己的线性代数库提供了很多工具。这尤其包括在集群上以%parallel方式工作的能力，使用MPI，以及更多种类的前置条件器。在后一类中，最有趣的能力之一是Trilinos ML包的存在，它实现了代数多栅（AMG）方法。我们将使用这个预处理程序对动量方程的二阶算子部分进行预处理。在步骤32中，我们将使用与这里讨论的相同的问题，探索以%并行方式解决问题的能力。

我们在第17步和第18步中使用的PETSc无疑是一个强大的库，它提供了大量处理矩阵、向量、迭代求解器和预处理器的函数，还有很多其他的东西，其中大部分在%parallel中运行得相当好。然而，它比Trilinos早了几年，是用C语言编写的，而且一般来说不像其他一些库那样容易使用。因此，deal.II也获得了与Trilinos的接口，Trilinos与PETSc有很多相同的功能。然而，它是一个年轻了好几年的项目，是用C++编写的，其作者一般都非常重视软件设计。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们在这里要解决的情况如下：我们用 $\kappa=10^{-6}, \eta=1, \rho=1, \beta=10$ 来解决上述的Boussinesq方程，即一个相对缓慢运动的流体，它几乎没有热扩散传导性，主要通过对流来传输热量。在边界上，我们将要求速度（ $\mathrm{n}\cdot\mathrm{u}=0$ ）和温度（ $\mathrm{n}\cdot\nabla T=0$ ）没有正态流量。这是在步骤22的介绍中讨论的情况之一，它固定了速度的一个分量，同时允许流动与边界%平行。还有 <code>dim-1</code> 分量需要固定，即法向应力的切向分量；对于这些分量，我们选择同质条件，这意味着我们不需要任何特殊条件。初始条件只对温度场是必要的，我们选择它为恒定的零。

然后，问题的演变完全由温度方程的右手边 $\gamma(\mathrm{x},t)$ 驱动，即由热源和汇驱动。在这里，我们选择了一个在圣诞讲座前发明的设置：美国的教室里当然禁止使用真实的蜡烛，但允许使用虚拟的蜡烛。因此，我们选择了三个球形的热源，不等距地靠近领域的底部，模仿三个蜡烛的样子。位于这些热源处的流体，最初处于静止状态，然后被加热，随着温度的升高，获得浮力，上升；更多的流体被拖上来，穿过热源，导致三个热羽上升，直到它们被外面下沉的流体循环所捕获，取代了因加热而上升的空气。


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
 * 像往常一样，第一步是包括这些著名的deal.II库文件和一些C++头文件的功能。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/block_sparsity_pattern.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * @endcode
 * 
 * 然后我们需要包括一些头文件，这些文件提供了矢量、矩阵和预处理类，这些类实现了各自Trilinos类的接口。特别是，我们将需要基于Trilinos的矩阵和向量类以及Trilinos预处理程序的接口。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_vector.h> 
 * #include <deal.II/lac/trilinos_parallel_block_vector.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * 
 * @endcode
 * 
 * 最后，这里有几个C++头文件还没有被上述头文件中的某个文件所包含。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * #include <memory> 
 * #include <limits> 
 * 
 * @endcode
 * 
 * 在这个顶层事项的最后，我们将所有deal.II的名字导入到全局命名空间。
 * 

 * 
 * 
 * @code
 * namespace Step31 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 同样，程序的下一阶段是定义方程数据，即各种边界条件、右手边和初始条件（记住，我们要解决的是一个时间依赖型系统）。这个定义的基本策略与  step-22  中的相同。不过关于细节，还是有一些区别。
 * 

 * 
 * 首先，我们没有在速度上设置任何不均匀的边界条件，因为正如介绍中所解释的，我们将使用无流条件  $\mathbf{n}\cdot\mathbf{u}=0$  。所以剩下的是应力张量法线分量的切向部分的条件 <code>dim-1</code> ， $\textbf{n} \cdot [p \textbf{1} - \eta\varepsilon(\textbf{u})]$ ；我们假定这些分量的值是同质的，也就是说，一个自然的边界条件，不需要具体的动作（它作为零项出现在弱形式的右边）。
 * 

 * 
 * 对于温度  $T$  ，我们假设没有热能通量，即  $\mathbf{n} \cdot \kappa \nabla T=0$  。这也是一个边界条件，不需要我们做任何特别的事情。
 * 

 * 
 * 第二，我们必须设定温度的初始条件（速度和压力不需要初始条件，因为我们在这里考虑的准稳态情况下的斯托克斯方程没有速度或压力的时间导数）。在这里，我们选择一个非常简单的测试案例，即初始温度为零，所有的动力学都由温度的右手边驱动。
 * 

 * 
 * 第三，我们需要定义温度方程的右边。我们选择它在域的底部某处的三个圆（或三维球）内为常数，如介绍中所解释的那样，而在域外为零。
 * 

 * 
 * 最后，或者说首先，在这个命名空间的顶部，我们定义我们需要的各种材料常数（ $\eta,\kappa$ ，密度 $\rho$ 和热膨胀系数 $\beta$ ）。
 * 

 * 
 * 
 * @code
 *   namespace EquationData 
 *   { 
 *     constexpr double eta     = 1; 
 *     constexpr double kappa   = 1e-6; 
 *     constexpr double beta    = 10; 
 *     constexpr double density = 1; 
 * 
 *     template <int dim> 
 *     class TemperatureInitialValues : public Function<dim> 
 *     { 
 *     public: 
 *       TemperatureInitialValues() 
 *         : Function<dim>(1) 
 *       {} 
 * 
 *       virtual double value(const Point<dim> & /*p*/, 
 *                            const unsigned int /*component*/ = 0) const override 
 *       { 
 *         return 0; 
 *       } 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  value) const override 
 *       { 
 *         for (unsigned int c = 0; c < this->n_components; ++c) 
 *           value(c) = TemperatureInitialValues<dim>::value(p, c); 
 *       } 
 *     }; 
 * 
 *     template <int dim> 
 *     class TemperatureRightHandSide : public Function<dim> 
 *     { 
 *     public: 
 *       TemperatureRightHandSide() 
 *         : Function<dim>(1) 
 *       {} 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override 
 *       { 
 *         (void)component; 
 *         Assert(component == 0, 
 *                ExcMessage("Invalid operation for a scalar function.")); 
 * 
 *         Assert((dim == 2) || (dim == 3), ExcNotImplemented()); 
 * 
 *         static const Point<dim> source_centers[3] = { 
 *           (dim == 2 ? Point<dim>(.3, .1) : Point<dim>(.3, .5, .1)), 
 *           (dim == 2 ? Point<dim>(.45, .1) : Point<dim>(.45, .5, .1)), 
 *           (dim == 2 ? Point<dim>(.75, .1) : Point<dim>(.75, .5, .1))}; 
 *         static const double source_radius = (dim == 2 ? 1. / 32 : 1. / 8); 
 * 
 *         return ((source_centers[0].distance(p) < source_radius) || 
 *                     (source_centers[1].distance(p) < source_radius) || 
 *                     (source_centers[2].distance(p) < source_radius) ? 
 *                   1 : 
 *                   0); 
 *       } 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  value) const override 
 *       { 
 *         for (unsigned int c = 0; c < this->n_components; ++c) 
 *           value(c) = TemperatureRightHandSide<dim>::value(p, c); 
 *       } 
 *     }; 
 *   } // namespace EquationData 
 * 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * 本节介绍了一些用于求解斯托克斯系统线性方程的对象，我们需要在每个时间步长中求解。这里使用的许多想法与 step-20 相同，其中介绍了基于Schur补的预处理程序和求解器，实际接口来自 step-22 （特别是 step-22 中 "结果 "部分的讨论，其中我们介绍了直接Schur补方法的替代品）。但是请注意，在这里我们不使用Schur补数来解决Stokes方程，尽管预处理程序中出现了一个近似的Schur补数（压力空间的质量矩阵）。
 * 

 * 
 * 
 * @code
 *   namespace LinearSolvers 
 *   { 
 * @endcode
 * 
 * 
 * <a name="ThecodeInverseMatrixcodeclasstemplate"></a> 
 * <h4>The <code>InverseMatrix</code> class template</h4>
 * 

 * 
 * 这个类是一个接口，用于计算 "倒置 "矩阵对向量的作用（使用 <code>vmult</code> 操作），其方式与 step-22 中的相应类相同：当请求这个类的对象的乘积时，我们使用CG方法解决与该矩阵有关的线性方程组，通过（模板化） <code>PreconditionerType</code> 类的预处理器加速。
 * 

 * 
 * 与 step-22 中同一类别的实现略有不同，我们让 <code>vmult</code> 函数接受任何类型的向量类型（但是，如果矩阵不允许与这种向量进行矩阵-向量乘积，它将产生编译器错误）。
 * 

 * 
 * 第二，我们捕捉解算器可能抛出的任何异常。原因如下。在调试这样的程序时，偶尔会犯一个错误，即把一个不确定或不对称的矩阵或预处理程序传递给当前的类。在这种情况下，求解器将不能收敛并抛出一个运行时异常。如果在这里没有被捕捉到，它就会在调用堆栈中传播，最后可能会在 <code>main()</code> 中出现，在那里我们会输出一个错误信息，说CG求解器失败。那么问题来了。哪个CG求解器？倒置质量矩阵的那个？用拉普拉斯算子反转左上角块的那个？还是在当前代码中我们使用线性求解器的其他几个嵌套位置中的一个CG求解器？在运行时异常中没有这方面的指示，因为它没有存储我们到达产生异常的地方的调用栈。
 * 所以
 * 与其让异常自由传播到 <code>main()</code> ，不如意识到如果内部求解器失败，外部函数能做的很少，不如将运行时异常转化为一个断言，该断言失败后会触发对 <code>abort()</code> 的调用，允许我们在调试器中追溯我们如何到达当前位置。
 * 

 * 
 * 
 * @code
 *     template <class MatrixType, class PreconditionerType> 
 *     class InverseMatrix : public Subscriptor 
 *     { 
 *     public: 
 *       InverseMatrix(const MatrixType &        m, 
 *                     const PreconditionerType &preconditioner); 
 * 
 *       template <typename VectorType> 
 *       void vmult(VectorType &dst, const VectorType &src) const; 
 * 
 *     private: 
 *       const SmartPointer<const MatrixType> matrix; 
 *       const PreconditionerType &           preconditioner; 
 *     }; 
 * 
 *     template <class MatrixType, class PreconditionerType> 
 *     InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
 *       const MatrixType &        m, 
 *       const PreconditionerType &preconditioner) 
 *       : matrix(&m) 
 *       , preconditioner(preconditioner) 
 *     {} 
 * 
 *     template <class MatrixType, class PreconditionerType> 
 *     template <typename VectorType> 
 *     void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
 *       VectorType &      dst, 
 *       const VectorType &src) const 
 *     { 
 *       SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm()); 
 *       SolverCG<VectorType> cg(solver_control); 
 * 
 *       dst = 0; 
 * 
 *       try 
 *         { 
 *           cg.solve(*matrix, dst, src, preconditioner); 
 *         } 
 *       catch (std::exception &e) 
 *         { 
 *           Assert(false, ExcMessage(e.what())); 
 *         } 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Schurcomplementpreconditioner"></a> 
 * <h4>Schur complement preconditioner</h4>
 * 

 * 
 * 这是在介绍中详细描述的舒尔补码预处理程序的实现。与 step-20 和 step-22 相反，我们使用GMRES一次性解决块系统，并使用块结构矩阵的Schur补码来建立一个良好的预处理程序。
 * 

 * 
 * 让我们看看介绍中描述的理想预处理矩阵  $P=\left(\begin{array}{cc} A & 0 \\ B & -S \end{array}\right)$  。如果我们在线性系统的求解中应用这个矩阵，迭代式GMRES求解器的收敛性将受矩阵
 * @f{eqnarray*} P^{-1}\left(\begin{array}{cc} A &
 * B^T \\ B & 0 \end{array}\right) = \left(\begin{array}{cc} I & A^{-1}
 * B^T \\ 0 & I \end{array}\right), 
 * @f}
 * 的制约，这确实非常简单。基于精确矩阵的GMRES求解器将在一次迭代中收敛，因为所有的特征值都是相等的（任何Krylov方法最多需要多少次迭代就有多少个不同的特征值）。Silvester和Wathen提出了这样一个用于受阻斯托克斯系统的预处理程序（"稳定的斯托克斯系统的快速迭代解第二部分。 Using general block preconditioners", SIAM J. Numer. Anal., 31 (1994), pp.1352-1367）。)
 * 

 * 
 * 用 $\tilde{P}$ 代替 $P$ 可以保持这种精神：乘积 $P^{-1} A$ 仍将接近于特征值为1的矩阵，其分布不取决于问题大小。这让我们希望能够得到一个与问题规模无关的GMRES迭代次数。
 * 

 * 
 * 已经通过 step-20 和 step-22 教程的deal.II用户当然可以想象我们将如何实现这一点。 我们用一些由InverseMatrix类构建的近似逆矩阵取代 $P^{-1}$ 中的精确逆矩阵，逆舒尔补码将由压力质量矩阵 $M_p$ 近似（如介绍中提到的由 $\eta^{-1}$ 加权）。正如在 step-22 的结果部分所指出的，我们可以通过应用一个预处理程序来取代 $A$ 的精确逆，在这种情况下，如介绍中所解释的那样，在一个矢量拉普拉斯矩阵上。这确实增加了（外部）GMRES的迭代次数，但仍然比精确的逆运算便宜得多，因为 <em> 的每个 </em> 外部求解器步骤（使用AMG预处理程序）需要20到35次CG迭代。
 * 

 * 
 * 考虑到上述解释，我们定义了一个具有 <code>vmult</code> 功能的预处理类，这就是我们在程序代码中进一步与通常的求解器函数交互所需要的。
 * 

 * 
 * 首先是声明。这与 step-20 中Schur补码的定义相似，不同的是我们在构造函数中需要更多的预处理程序，而且我们在这里使用的矩阵是建立在Trilinos之上的。
 * 

 * 
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     class BlockSchurPreconditioner : public Subscriptor 
 *     { 
 *     public: 
 *       BlockSchurPreconditioner( 
 *         const TrilinosWrappers::BlockSparseMatrix &S, 
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                             PreconditionerTypeMp> &Mpinv, 
 *         const PreconditionerTypeA &                Apreconditioner); 
 * 
 *       void vmult(TrilinosWrappers::MPI::BlockVector &      dst, 
 *                  const TrilinosWrappers::MPI::BlockVector &src) const; 
 * 
 *     private: 
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> 
 *         stokes_matrix; 
 *       const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                                              PreconditionerTypeMp>> 
 *                                  m_inverse; 
 *       const PreconditionerTypeA &a_preconditioner; 
 * 
 *       mutable TrilinosWrappers::MPI::Vector tmp; 
 *     }; 
 * 
 * @endcode
 * 
 * 当使用 TrilinosWrappers::MPI::Vector 或 TrilinosWrappers::MPI::BlockVector, 时，Vector被使用IndexSet初始化。IndexSet不仅用于调整 TrilinosWrappers::MPI::Vector 的大小，而且还将 TrilinosWrappers::MPI::Vector 中的一个索引与一个自由度联系起来（更详细的解释见 step-40 ）。函数complete_index_set()创建了一个IndexSet，每个有效的索引都是这个集合的一部分。请注意，这个程序只能按顺序运行，如果并行使用，将抛出一个异常。
 * 

 * 
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>:: 
 *       BlockSchurPreconditioner( 
 *         const TrilinosWrappers::BlockSparseMatrix &S, 
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                             PreconditionerTypeMp> &Mpinv, 
 *         const PreconditionerTypeA &                Apreconditioner) 
 *       : stokes_matrix(&S) 
 *       , m_inverse(&Mpinv) 
 *       , a_preconditioner(Apreconditioner) 
 *       , tmp(complete_index_set(stokes_matrix->block(1, 1).m())) 
 *     {} 
 * 
 * @endcode
 * 
 * 接下来是 <code>vmult</code> 函数。我们以三个连续的步骤实现上述 $P^{-1}$ 的动作。 在公式中，我们要计算 $Y=P^{-1}X$ ，其中 $X,Y$ 都是有两个块成分的向量。
 * 

 * 
 * 第一步用矩阵 $A$ 的预处理乘以矢量的速度部分，即计算 $Y_0={\tilde A}^{-1}X_0$  。 然后将得到的速度矢量乘以 $B$ 并减去压力，即我们要计算 $X_1-BY_0$  。这第二步只作用于压力向量，由我们矩阵类的残差函数完成，只是符号不对。因此，我们改变临时压力向量中的符号，最后乘以反压力质量矩阵，得到最终的压力向量，完成我们对斯托克斯预处理的工作。
 * 

 * 
 * 
 * @code
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     void 
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult( 
 *       TrilinosWrappers::MPI::BlockVector &      dst, 
 *       const TrilinosWrappers::MPI::BlockVector &src) const 
 *     { 
 *       a_preconditioner.vmult(dst.block(0), src.block(0)); 
 *       stokes_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1)); 
 *       tmp *= -1; 
 *       m_inverse->vmult(dst.block(1), tmp); 
 *     } 
 *   } // namespace LinearSolvers 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeBoussinesqFlowProblemcodeclasstemplate"></a> 
 * <h3>The <code>BoussinesqFlowProblem</code> class template</h3>
 * 

 * 
 * 定义了解决随时间变化的Boussinesq问题的顶层逻辑的类的定义主要是基于 step-22 的教程程序。主要的区别在于，现在我们还必须求解温度方程，这迫使我们为温度变量准备第二个DoFHandler对象，以及当前和之前的时间步骤的矩阵、右手边和求解向量。正如介绍中提到的，所有的线性代数对象都将使用相应的Trilinos功能的包装器。
 * 

 * 
 * 这个类的成员函数让人想起 step-21 ，在那里我们也使用了一个交错的方案，首先解决流动方程（这里是斯托克斯方程， step-21 是达西流），然后更新平流量（这里是温度，那里是饱和度）。新的函数主要涉及到确定时间步长，以及人工粘性稳定的适当大小。
 * 

 * 
 * 最后三个变量表示在下次调用相应的建立函数时，是否需要重建各种矩阵或预处理程序。这使得我们可以将相应的 <code>if</code> 移到相应的函数中，从而使我们的主 <code>run()</code> 函数保持干净，易于阅读。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BoussinesqFlowProblem 
 *   { 
 *   public: 
 *     BoussinesqFlowProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void   setup_dofs(); 
 *     void   assemble_stokes_preconditioner(); 
 *     void   build_stokes_preconditioner(); 
 *     void   assemble_stokes_system(); 
 *     void   assemble_temperature_system(const double maximal_velocity); 
 *     void   assemble_temperature_matrix(); 
 *     double get_maximal_velocity() const; 
 *     std::pair<double, double> get_extrapolated_temperature_range() const; 
 *     void                      solve(); 
 *     void                      output_results() const; 
 *     void                      refine_mesh(const unsigned int max_grid_level); 
 * 
 *     double compute_viscosity( 
 *       const std::vector<double> &        old_temperature, 
 *       const std::vector<double> &        old_old_temperature, 
 *       const std::vector<Tensor<1, dim>> &old_temperature_grads, 
 *       const std::vector<Tensor<1, dim>> &old_old_temperature_grads, 
 *       const std::vector<double> &        old_temperature_laplacians, 
 *       const std::vector<double> &        old_old_temperature_laplacians, 
 *       const std::vector<Tensor<1, dim>> &old_velocity_values, 
 *       const std::vector<Tensor<1, dim>> &old_old_velocity_values, 
 *       const std::vector<double> &        gamma_values, 
 *       const double                       global_u_infty, 
 *       const double                       global_T_variation, 
 *       const double                       cell_diameter) const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     double             global_Omega_diameter; 
 * 
 *     const unsigned int        stokes_degree; 
 *     FESystem<dim>             stokes_fe; 
 *     DoFHandler<dim>           stokes_dof_handler; 
 *     AffineConstraints<double> stokes_constraints; 
 * 
 *     std::vector<IndexSet>               stokes_partitioning; 
 *     TrilinosWrappers::BlockSparseMatrix stokes_matrix; 
 *     TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix; 
 * 
 *     TrilinosWrappers::MPI::BlockVector stokes_solution; 
 *     TrilinosWrappers::MPI::BlockVector old_stokes_solution; 
 *     TrilinosWrappers::MPI::BlockVector stokes_rhs; 
 * 
 *     const unsigned int        temperature_degree; 
 *     FE_Q<dim>                 temperature_fe; 
 *     DoFHandler<dim>           temperature_dof_handler; 
 *     AffineConstraints<double> temperature_constraints; 
 * 
 *     TrilinosWrappers::SparseMatrix temperature_mass_matrix; 
 *     TrilinosWrappers::SparseMatrix temperature_stiffness_matrix; 
 *     TrilinosWrappers::SparseMatrix temperature_matrix; 
 * 
 *     TrilinosWrappers::MPI::Vector temperature_solution; 
 *     TrilinosWrappers::MPI::Vector old_temperature_solution; 
 *     TrilinosWrappers::MPI::Vector old_old_temperature_solution; 
 *     TrilinosWrappers::MPI::Vector temperature_rhs; 
 * 
 *     double       time_step; 
 *     double       old_time_step; 
 *     unsigned int timestep_number; 
 * 
 *     std::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner; 
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC>  Mp_preconditioner; 
 * 
 *     bool rebuild_stokes_matrix; 
 *     bool rebuild_temperature_matrices; 
 *     bool rebuild_stokes_preconditioner; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemclassimplementation"></a> 
 * <h3>BoussinesqFlowProblem class implementation</h3>
 * 
 * <a name="BoussinesqFlowProblemBoussinesqFlowProblem"></a> 
 * <h4>BoussinesqFlowProblem::BoussinesqFlowProblem</h4>
 * 

 * 
 * 这个类的构造函数是对  step-22  中的构造函数的扩展。我们需要添加涉及温度的各种变量。正如介绍中所讨论的，我们将再次使用 $Q_2\times Q_1$ （Taylor-Hood）元素来表示斯托克斯部分，并使用 $Q_2$ 元素表示温度。然而，通过使用存储斯托克斯和温度有限元的多项式程度的变量，可以很容易地持续修改这些元素的程度以及下游使用的所有正交公式。此外，我们还初始化了时间步长以及矩阵组合和预处理的选项。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BoussinesqFlowProblem<dim>::BoussinesqFlowProblem() 
 *     : triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()) 
 *     , stokes_degree(1) 
 *     , stokes_fe(FE_Q<dim>(stokes_degree + 1), dim, FE_Q<dim>(stokes_degree), 1) 
 *     , stokes_dof_handler(triangulation) 
 *     , 
 * 
 *     temperature_degree(2) 
 *     , temperature_fe(temperature_degree) 
 *     , temperature_dof_handler(triangulation) 
 *     , 
 * 
 *     time_step(0) 
 *     , old_time_step(0) 
 *     , timestep_number(0) 
 *     , rebuild_stokes_matrix(true) 
 *     , rebuild_temperature_matrices(true) 
 *     , rebuild_stokes_preconditioner(true) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_maximal_velocity"></a> 
 * <h4>BoussinesqFlowProblem::get_maximal_velocity</h4>
 * 

 * 
 * 开始这个类的真正功能是一个辅助函数，确定域内（事实上是正交点）的最大（ $L_\infty$  ）速度。它是如何工作的，对所有已经达到本教程这一点的人来说应该是比较明显的。请注意，由于我们只对速度感兴趣，我们不使用 <code>stokes_fe_values.get_function_values</code> 来获取整个斯托克斯解的值（速度和压力），而是使用 <code>stokes_fe_values[velocities].get_function_values</code> 来提取速度部分。这样做的额外好处是，我们得到的是张量<1,dim>，而不是向量<double>中的一些分量，这样我们就可以马上用 <code>norm()</code> 函数来处理它，得到速度的大小。
 * 

 * 
 * 唯一值得思考的一点是如何选择我们在这里使用的正交点。由于这个函数的目标是通过查看每个单元格上的正交点来寻找域内的最大速度。所以我们应该问，我们应该如何最好地选择每个单元上的这些正交点。为此，回顾一下，如果我们有一个单一的 $Q_1$ 场（而不是高阶的矢量值场），那么最大值将在网格的一个顶点达到。换句话说，我们应该使用QTrapezoid类，它的正交点只在单元的顶点。
 * 

 * 
 * 对于高阶形状函数，情况更为复杂：最大值和最小值可能在形状函数的支持点之间达到（对于通常的 $Q_p$ 元素，支持点是等距的Lagrange插值点）；此外，由于我们正在寻找一个矢量值的最大幅值，我们更不能肯定地说潜在的最大点集合在哪里。然而，从直觉上讲，即使不能证明，拉格朗日插值点似乎也是比高斯点更好的选择。
 * 

 * 
 * 现在有不同的方法来产生一个正交公式，其正交点等于有限元的插值点。一种选择是使用 FiniteElement::get_unit_support_points() 函数，将输出减少到一组唯一的点以避免重复的函数评估，并使用这些点创建一个正交对象。另一个选择，这里选择的是使用QTrapezoid类，并将其与QIterated类相结合，该类在每个坐标方向的若干子单元上重复QTrapezoid公式。为了覆盖所有的支持点，我们需要对其进行 <code>stokes_degree+1</code> 次迭代，因为这是使用中的斯托克斯元素的多项式程度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double BoussinesqFlowProblem<dim>::get_maximal_velocity() const 
 *   { 
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(), stokes_degree + 1); 
 *     const unsigned int   n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(stokes_fe, quadrature_formula, update_values); 
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points); 
 *     double                      max_velocity = 0; 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         fe_values[velocities].get_function_values(stokes_solution, 
 *                                                   velocity_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           max_velocity = std::max(max_velocity, velocity_values[q].norm()); 
 *       } 
 * 
 *     return max_velocity; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemget_extrapolated_temperature_range"></a> 
 * <h4>BoussinesqFlowProblem::get_extrapolated_temperature_range</h4>
 * 

 * 
 * 接下来是一个函数，确定从前两个时间步长推算到当前步长时， $\Omega$ 内正交点的最低和最高温度。我们在计算人工粘性参数 $\nu$ 时需要这个信息，正如在介绍中所讨论的那样。
 * 

 * 
 * 外推温度的公式是  $\left(1+\frac{k_n}{k_{n-1}} \right)T^{n-1} + \frac{k_n}{k_{n-1}} T^{n-2}$  。计算的方法是在所有正交点上循环，如果当前值比前一个值大/小，则更新最大和最小值。在对所有正交点进行循环之前，我们将存储最大和最小值的变量初始化为可表示为双数的最小和最大的数字。这样我们就知道它比最小/最大值大/小，并且所有正交点的循环最终会用正确的值来更新初始值。
 * 

 * 
 * 这里唯一值得一提的复杂情况是，在第一个时间步骤中， $T^{k-2}$ 当然还不能使用。在这种情况下，我们只能使用 $T^{k-1}$ ，这是我们从初始温度得到的。作为正交点，我们使用与前一个函数相同的选择，但不同的是，现在重复的数量由温度场的多项式程度决定。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<double, double> 
 *   BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range() const 
 *   { 
 *     const QIterated<dim> quadrature_formula(QTrapezoid<1>(), 
 *                                             temperature_degree); 
 *     const unsigned int   n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(temperature_fe, quadrature_formula, update_values); 
 *     std::vector<double> old_temperature_values(n_q_points); 
 *     std::vector<double> old_old_temperature_values(n_q_points); 
 * 
 *     if (timestep_number != 0) 
 *       { 
 *         double min_temperature = std::numeric_limits<double>::max(), 
 *                max_temperature = -std::numeric_limits<double>::max(); 
 * 
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
 *           { 
 *             fe_values.reinit(cell); 
 *             fe_values.get_function_values(old_temperature_solution, 
 *                                           old_temperature_values); 
 *             fe_values.get_function_values(old_old_temperature_solution, 
 *                                           old_old_temperature_values); 
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q) 
 *               { 
 *                 const double temperature = 
 *                   (1. + time_step / old_time_step) * old_temperature_values[q] - 
 *                   time_step / old_time_step * old_old_temperature_values[q]; 
 * 
 *                 min_temperature = std::min(min_temperature, temperature); 
 *                 max_temperature = std::max(max_temperature, temperature); 
 *               } 
 *           } 
 * 
 *         return std::make_pair(min_temperature, max_temperature); 
 *       } 
 *     else 
 *       { 
 *         double min_temperature = std::numeric_limits<double>::max(), 
 *                max_temperature = -std::numeric_limits<double>::max(); 
 * 
 *         for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
 *           { 
 *             fe_values.reinit(cell); 
 *             fe_values.get_function_values(old_temperature_solution, 
 *                                           old_temperature_values); 
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q) 
 *               { 
 *                 const double temperature = old_temperature_values[q]; 
 * 
 *                 min_temperature = std::min(min_temperature, temperature); 
 *                 max_temperature = std::max(max_temperature, temperature); 
 *               } 
 *           } 
 * 
 *         return std::make_pair(min_temperature, max_temperature); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemcompute_viscosity"></a> 
 * <h4>BoussinesqFlowProblem::compute_viscosity</h4>
 * 

 * 
 * 最后一个工具函数计算单元 $\nu|_K$ 上的人工粘度参数 $K$ ，作为外推温度、其梯度和Hessian（二阶导数）、速度、当前单元正交点上的所有右手 $\gamma$ 和其他各种参数的函数，在介绍中已详细说明。
 * 

 * 
 * 这里有一些值得一提的通用常数。首先，我们需要固定 $\beta$ ；我们选择 $\beta=0.017\cdot dim$ ，这个选择在本教程程序的结果部分有详细讨论。其次是指数 $\alpha$ ； $\alpha=1$ 对于目前的程序似乎很好用，尽管选择 $\alpha = 2$ 可能会有一些额外的好处。最后，有一件事需要特别说明。在第一个时间步骤中，速度等于零， $\nu|_K$ 的公式没有定义。在这种情况下，我们返回 $\nu|_K=5\cdot 10^3 \cdot h_K$ ，这个选择无疑更多的是出于启发式的考虑（不过，它与第二个时间步骤中大多数单元的返回值处于同一数量级）。
 * 

 * 
 * 根据介绍中讨论的材料，该函数的其余部分应该是显而易见的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double BoussinesqFlowProblem<dim>::compute_viscosity( 
 *     const std::vector<double> &        old_temperature, 
 *     const std::vector<double> &        old_old_temperature, 
 *     const std::vector<Tensor<1, dim>> &old_temperature_grads, 
 *     const std::vector<Tensor<1, dim>> &old_old_temperature_grads, 
 *     const std::vector<double> &        old_temperature_laplacians, 
 *     const std::vector<double> &        old_old_temperature_laplacians, 
 *     const std::vector<Tensor<1, dim>> &old_velocity_values, 
 *     const std::vector<Tensor<1, dim>> &old_old_velocity_values, 
 *     const std::vector<double> &        gamma_values, 
 *     const double                       global_u_infty, 
 *     const double                       global_T_variation, 
 *     const double                       cell_diameter) const 
 *   { 
 *     constexpr double beta  = 0.017 * dim; 
 *     constexpr double alpha = 1.0; 
 * 
 *     if (global_u_infty == 0) 
 *       return 5e-3 * cell_diameter; 
 * 
 *     const unsigned int n_q_points = old_temperature.size(); 
 * 
 *     double max_residual = 0; 
 *     double max_velocity = 0; 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         const Tensor<1, dim> u = 
 *           (old_velocity_values[q] + old_old_velocity_values[q]) / 2; 
 * 
 *         const double dT_dt = 
 *           (old_temperature[q] - old_old_temperature[q]) / old_time_step; 
 *         const double u_grad_T = 
 *           u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2; 
 * 
 *         const double kappa_Delta_T = 
 *           EquationData::kappa * 
 *           (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) / 
 *           2; 
 * 
 *         const double residual = 
 *           std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) * 
 *                    std::pow((old_temperature[q] + old_old_temperature[q]) / 2, 
 *                             alpha - 1.)); 
 * 
 *         max_residual = std::max(residual, max_residual); 
 *         max_velocity = std::max(std::sqrt(u * u), max_velocity); 
 *       } 
 * 
 *     const double c_R            = std::pow(2., (4. - 2 * alpha) / dim); 
 *     const double global_scaling = c_R * global_u_infty * global_T_variation * 
 *                                   std::pow(global_Omega_diameter, alpha - 2.); 
 * 
 *     return ( 
 *       beta * max_velocity * 
 *       std::min(cell_diameter, 
 *                std::pow(cell_diameter, alpha) * max_residual / global_scaling)); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemsetup_dofs"></a> 
 * <h4>BoussinesqFlowProblem::setup_dofs</h4>
 * 

 * 
 * 这是一个函数，用于设置我们这里的DoFHandler对象（一个用于斯托克斯部分，一个用于温度部分），以及将本程序中线性代数所需的各种对象设置为合适的尺寸。它的基本操作与我们在  step-22  中的操作类似。
 * 

 * 
 * 该函数的主体首先列举了斯托克斯和温度系统的所有自由度。对于斯托克斯部分，自由度被排序，以确保速度优先于压力自由度，这样我们就可以将斯托克斯矩阵划分为一个 $2\times 2$ 矩阵。作为与 step-22 的区别，我们不进行任何额外的DoF重新编号。在那个程序中，它得到了回报，因为我们的求解器严重依赖ILU，而我们在这里使用AMG，它对DoF编号不敏感。用于压力质量矩阵反演的IC预处理程序当然会利用类似Cuthill-McKee的重新编号，但是与速度部分相比，其成本很低，所以额外的工作并没有得到回报。
 * 

 * 
 * 然后，我们继续生成悬挂的节点约束，这些约束来自两个DoFHandler对象的自适应网格细化。对于速度，我们通过向已经存储了悬挂节点约束矩阵的对象添加约束来施加无流边界条件 $\mathbf{u}\cdot \mathbf{n}=0$ 。函数中的第二个参数描述了总dof向量中的第一个速度分量，这里是零。变量 <code>no_normal_flux_boundaries</code> 表示要设置无通量边界条件的边界指标；这里是边界指标0。
 * 

 * 
 * 做完这些后，我们计算各块中的自由度数量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::setup_dofs() 
 *   { 
 *     std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0); 
 *     stokes_sub_blocks[dim] = 1; 
 * 
 *     { 
 *       stokes_dof_handler.distribute_dofs(stokes_fe); 
 *       DoFRenumbering::component_wise(stokes_dof_handler, stokes_sub_blocks); 
 * 
 *       stokes_constraints.clear(); 
 *       DoFTools::make_hanging_node_constraints(stokes_dof_handler, 
 *                                               stokes_constraints); 
 *       std::set<types::boundary_id> no_normal_flux_boundaries; 
 *       no_normal_flux_boundaries.insert(0); 
 *       VectorTools::compute_no_normal_flux_constraints(stokes_dof_handler, 
 *                                                       0, 
 *                                                       no_normal_flux_boundaries, 
 *                                                       stokes_constraints); 
 *       stokes_constraints.close(); 
 *     } 
 *     { 
 *       temperature_dof_handler.distribute_dofs(temperature_fe); 
 * 
 *       temperature_constraints.clear(); 
 *       DoFTools::make_hanging_node_constraints(temperature_dof_handler, 
 *                                               temperature_constraints); 
 *       temperature_constraints.close(); 
 *     } 
 * 
 *     const std::vector<types::global_dof_index> stokes_dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(stokes_dof_handler, stokes_sub_blocks); 
 * 
 *     const unsigned int n_u = stokes_dofs_per_block[0], 
 *                        n_p = stokes_dofs_per_block[1], 
 *                        n_T = temperature_dof_handler.n_dofs(); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << " (on " << triangulation.n_levels() << " levels)" << std::endl 
 *               << "Number of degrees of freedom: " << n_u + n_p + n_T << " (" 
 *               << n_u << '+' << n_p << '+' << n_T << ')' << std::endl 
 *               << std::endl; 
 * 
 * @endcode
 * 
 * 下一步是创建斯托克斯和温度系统矩阵的稀疏模式，以及建立斯托克斯预处理矩阵的预处理。如同在 step-22 中一样，我们选择使用DynamicSparsityPattern的封锁版本来创建模式。
 * 

 * 
 * 因此，我们首先释放存储在矩阵中的内存，然后建立一个BlockDynamicSparsityPattern类型的对象，该对象由 $2\times 2$ 块（用于斯托克斯系统矩阵和预处理器）或DynamicSparsityPattern（用于温度部分）组成。然后我们用非零模式填充这些对象，考虑到对于斯托克斯系统矩阵，在压力-压力块中没有条目（但所有速度矢量分量相互耦合并与压力耦合）。同样，在斯托克斯预处理矩阵中，只有对角线块是非零的，因为我们使用了介绍中讨论的矢量拉普拉斯。这个算子只把拉普拉斯的每个矢量分量与它自己联系起来，而不是与其他矢量分量联系起来。然而，应用无流量边界条件产生的约束条件将在边界处再次耦合向量分量）。
 * 

 * 
 * 在生成稀疏模式时，我们直接应用悬挂节点和无流边界条件的约束。这种方法在 step-27 中已经使用过了，但与早期教程中的方法不同，在早期教程中我们先建立原始的稀疏模式，然后才加入约束条件产生的条目。这样做的原因是，在以后的装配过程中，我们要在将本地道夫转移到全局道夫时立即分配约束。因此，在受限自由度的位置不会有数据写入，所以我们可以通过将最后一个布尔标志设置为 <code>false</code> ，让 DoFTools::make_sparsity_pattern 函数省略这些条目。一旦稀疏性模式准备好了，我们就可以用它来初始化特里诺斯矩阵。由于Trilinos矩阵在内部存储了稀疏模式，所以在初始化矩阵之后，没有必要再保留稀疏模式。
 * 

 * 
 * 
 * @code
 *     stokes_partitioning.resize(2); 
 *     stokes_partitioning[0] = complete_index_set(n_u); 
 *     stokes_partitioning[1] = complete_index_set(n_p); 
 *     { 
 *       stokes_matrix.clear(); 
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2); 
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u); 
 *       dsp.block(0, 1).reinit(n_u, n_p); 
 *       dsp.block(1, 0).reinit(n_p, n_u); 
 *       dsp.block(1, 1).reinit(n_p, n_p); 
 * 
 *       dsp.collect_sizes(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 * 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (!((c == dim) && (d == dim))) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         stokes_dof_handler, coupling, dsp, stokes_constraints, false); 
 * 
 *       stokes_matrix.reinit(dsp); 
 *     } 
 * 
 *     { 
 *       Amg_preconditioner.reset(); 
 *       Mp_preconditioner.reset(); 
 *       stokes_preconditioner_matrix.clear(); 
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2); 
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u); 
 *       dsp.block(0, 1).reinit(n_u, n_p); 
 *       dsp.block(1, 0).reinit(n_p, n_u); 
 *       dsp.block(1, 1).reinit(n_p, n_p); 
 * 
 *       dsp.collect_sizes(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (c == d) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         stokes_dof_handler, coupling, dsp, stokes_constraints, false); 
 * 
 *       stokes_preconditioner_matrix.reinit(dsp); 
 *     } 
 * 
 * @endcode
 * 
 * 温度矩阵（或者说是矩阵，因为我们提供了一个温度质量矩阵和一个温度刚度矩阵，它们将在时间离散化中被加在一起）的创建与斯托克斯矩阵的生成相同；只是在这里要简单得多，因为我们不需要照顾任何块或组件之间的耦合。注意我们是如何初始化三个温度矩阵的。我们只使用稀疏模式对第一个矩阵进行再初始化，而对其余两个再初始化则使用先前生成的矩阵。这样做的原因是，从一个已经生成的矩阵进行重新初始化，可以让Trilinos重新使用稀疏模式，而不是为每个副本生成一个新的模式。这样可以节省一些时间和内存。
 * 

 * 
 * 
 * @code
 *     { 
 *       temperature_mass_matrix.clear(); 
 *       temperature_stiffness_matrix.clear(); 
 *       temperature_matrix.clear(); 
 * 
 *       DynamicSparsityPattern dsp(n_T, n_T); 
 *       DoFTools::make_sparsity_pattern(temperature_dof_handler, 
 *                                       dsp, 
 *                                       temperature_constraints, 
 *                                       false); 
 * 
 *       temperature_matrix.reinit(dsp); 
 *       temperature_mass_matrix.reinit(temperature_matrix); 
 *       temperature_stiffness_matrix.reinit(temperature_matrix); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，我们将斯托克斯解的向量 $\mathbf u^{n-1}$ 和 $\mathbf u^{n-2}$ ，以及温度 $T^{n}$ 、 $T^{n-1}$ 和 $T^{n-2}$ （时间步进所需）和所有系统的右手边设置为正确的大小和块结构。
 * 

 * 
 * 
 * @code
 *     IndexSet temperature_partitioning = complete_index_set(n_T); 
 *     stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD); 
 *     old_stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD); 
 *     stokes_rhs.reinit(stokes_partitioning, MPI_COMM_WORLD); 
 * 
 *     temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD); 
 *     old_temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD); 
 *     old_old_temperature_solution.reinit(temperature_partitioning, 
 *                                         MPI_COMM_WORLD); 
 * 
 *     temperature_rhs.reinit(temperature_partitioning, MPI_COMM_WORLD); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_stokes_preconditioner"></a> 
 * <h4>BoussinesqFlowProblem::assemble_stokes_preconditioner</h4>
 * 

 * 
 * 这个函数组装了我们用于预处理斯托克斯系统的矩阵。我们需要的是速度分量上的矢量拉普拉斯矩阵和压力分量上的质量矩阵，并以 $\eta^{-1}$ 加权。我们首先生成一个适当阶数的正交对象，即FEValues对象，它可以给出正交点的值和梯度（连同正交权重）。接下来我们为单元格矩阵和局部与全局DoF之间的关系创建数据结构。向量 <code>grad_phi_u</code> and <code>phi_p</code> 将保存基函数的值，以便更快地建立局部矩阵，正如在 step-22 中已经完成的那样。在我们开始对所有活动单元进行循环之前，我们必须指定哪些成分是压力，哪些是速度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner() 
 *   { 
 *     stokes_preconditioner_matrix = 0; 
 * 
 *     const QGauss<dim> quadrature_formula(stokes_degree + 2); 
 *     FEValues<dim>     stokes_fe_values(stokes_fe, 
 *                                    quadrature_formula, 
 *                                    update_JxW_values | update_values | 
 *                                      update_gradients); 
 * 
 *     const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
 *     std::vector<double>         phi_p(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 *     for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
 *       { 
 *         stokes_fe_values.reinit(cell); 
 *         local_matrix = 0; 
 * 
 * @endcode
 * 
 * 本地矩阵的创建相当简单。只有一个拉普拉斯项（关于速度）和一个由 $\eta^{-1}$ 加权的质量矩阵需要生成，所以本地矩阵的创建在两行中完成。一旦本地矩阵准备好了（在每个正交点上循环查看本地矩阵的行和列），我们就可以得到本地的DoF指数，并将本地信息写入全局矩阵中。我们像在 step-27 中那样做，也就是说，我们直接应用本地悬挂节点的约束。这样做，我们就不必事后再做，而且我们也不会在消除约束时将矩阵的条目写成实际上将再次设置为零。
 * 

 * 
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 grad_phi_u[k] = stokes_fe_values[velocities].gradient(k, q); 
 *                 phi_p[k]      = stokes_fe_values[pressure].value(k, q); 
 *               } 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 local_matrix(i, j) += 
 *                   (EquationData::eta * 
 *                      scalar_product(grad_phi_u[i], grad_phi_u[j]) + 
 *                    (1. / EquationData::eta) * phi_p[i] * phi_p[j]) * 
 *                   stokes_fe_values.JxW(q); 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         stokes_constraints.distribute_local_to_global( 
 *           local_matrix, local_dof_indices, stokes_preconditioner_matrix); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblembuild_stokes_preconditioner"></a> 
 * <h4>BoussinesqFlowProblem::build_stokes_preconditioner</h4>
 * 

 * 
 * 这个函数生成将用于Schur互补块预处理的内部预处理。由于只有当矩阵发生变化时才需要重新生成预处理程序，因此在矩阵没有变化的情况下，该函数不需要做任何事情（即标志 <code>rebuild_stokes_preconditioner</code> 的值为 <code>false</code> ）。否则，它的第一个任务是调用 <code>assemble_stokes_preconditioner</code> 来生成预处理矩阵。
 * 

 * 
 * 接下来，我们为速度-速度矩阵  $A$  设置预处理程序。正如介绍中所解释的，我们将使用基于矢量拉普拉斯矩阵 $\hat{A}$ 的AMG预处理器（它在频谱上与斯托克斯矩阵 $A$ 接近）。通常， TrilinosWrappers::PreconditionAMG 类可以被看作是一个好的黑箱预处理程序，不需要任何特殊的知识。然而，在这种情况下，我们必须小心：因为我们为一个矢量问题建立了一个AMG，我们必须告诉预处理程序设置哪个道夫属于哪个矢量成分。我们使用 DoFTools::extract_constant_modes, 函数来做这件事，该函数生成一组 <code>dim</code> 向量，其中每个向量在向量问题的相应分量中为1，在其他地方为0。因此，这些是每个分量上的常数模式，这解释了变量的名称。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::build_stokes_preconditioner() 
 *   { 
 *     if (rebuild_stokes_preconditioner == false) 
 *       return; 
 * 
 *     std::cout << "   Rebuilding Stokes preconditioner..." << std::flush; 
 * 
 *     assemble_stokes_preconditioner(); 
 * 
 *     Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>(); 
 * 
 *     std::vector<std::vector<bool>> constant_modes; 
 *     FEValuesExtractors::Vector     velocity_components(0); 
 *     DoFTools::extract_constant_modes(stokes_dof_handler, 
 *                                      stokes_fe.component_mask( 
 *                                        velocity_components), 
 *                                      constant_modes); 
 *     TrilinosWrappers::PreconditionAMG::AdditionalData amg_data; 
 *     amg_data.constant_modes = constant_modes; 
 * 
 * @endcode
 * 
 * 接下来，我们再设置一些AMG预处理程序的选项。特别是，我们需要告诉AMG设置，我们对速度矩阵使用二次基函数（这意味着矩阵中有更多的非零元素，因此需要在内部选择一种更稳健的算法）。此外，我们希望能够控制粗化结构的建立方式。Trilinos平滑聚合AMG的方法是寻找哪些矩阵条目与对角线条目大小相似，以便代数式地建立一个粗网格结构。通过将参数 <code>aggregation_threshold</code> 设置为0.02，我们指定所有尺寸超过该行中一些对角线枢轴的百分之二的条目应该形成一个粗网格点。这个参数是比较特别的，对它进行一些微调会影响预处理程序的性能。根据经验，较大的 <code>aggregation_threshold</code> 值会减少迭代次数，但增加每次迭代的成本。看一下Trilinos的文档会提供更多关于这些参数的信息。有了这个数据集，我们就用我们想要的矩阵来初始化预处理程序。
 * 

 * 
 * 最后，我们也初始化预处理程序以反转压力质量矩阵。这个矩阵是对称的，表现良好，所以我们可以选择一个简单的预处理程序。我们坚持使用不完全Cholesky（IC）因子化预处理器，它是为对称矩阵设计的。我们也可以选择SSOR预处理器，其松弛系数约为1.2，但IC对我们的例子来说更便宜。我们把预处理程序包成一个 <code>std::shared_ptr</code> 指针，这使得下次重新创建预处理程序更加容易，因为我们不必关心破坏以前使用的对象。
 * 

 * 
 * 
 * @code
 *     amg_data.elliptic              = true; 
 *     amg_data.higher_order_elements = true; 
 *     amg_data.smoother_sweeps       = 2; 
 *     amg_data.aggregation_threshold = 0.02; 
 *     Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0, 0), 
 *                                    amg_data); 
 * 
 *     Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
 *     Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1, 1)); 
 * 
 *     std::cout << std::endl; 
 * 
 *     rebuild_stokes_preconditioner = false; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_stokes_system"></a> 
 * <h4>BoussinesqFlowProblem::assemble_stokes_system</h4>
 * 

 * 
 * 我们用于推进耦合的斯托克斯-温度系统的时滞方案迫使我们将装配（以及线性系统的解）分成两步。第一步是创建斯托克斯系统的矩阵和右手边，第二步是创建温度道夫的矩阵和右手边，这取决于速度的线性系统的结果。
 * 

 * 
 * 该函数在每个时间步长的开始时被调用。在第一个时间步骤中，或者如果网格已经改变，由 <code>rebuild_stokes_matrix</code> 表示，我们需要组装斯托克斯矩阵；另一方面，如果网格没有改变，矩阵已经有了，这就没有必要了，我们需要做的就是组装右手边的向量，它在每个时间步骤中都会改变。
 * 

 * 
 * 关于实现的技术细节，与  step-22  相比没有太大变化。我们重置矩阵和向量，在单元格上创建正交公式，然后创建相应的FEValues对象。对于更新标志，我们只在完全装配的情况下需要基函数导数，因为右手边不需要它们；像往常一样，根据当前需要选择最小的标志集，使程序中进一步调用  FEValues::reinit  的效率更高。
 * 

 * 
 * 有一件事需要评论&ndash；因为我们有一个单独的有限元和DoFHandler来处理温度问题，所以我们需要生成第二个FEValues对象来正确评估温度解决方案。要实现这一点并不复杂：只需使用温度结构，并为我们需要用于评估温度解决方案的基函数值设置一个更新标志。这里需要记住的唯一重要部分是，两个FEValues对象使用相同的正交公式，以确保我们在循环计算两个对象的正交点时得到匹配的信息。
 * 

 * 
 * 声明的过程中，有一些关于数组大小的快捷方式，本地矩阵和右手的创建，以及与全局系统相比，本地道夫的索引的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::assemble_stokes_system() 
 *   { 
 *     std::cout << "   Assembling..." << std::flush; 
 * 
 *     if (rebuild_stokes_matrix == true) 
 *       stokes_matrix = 0; 
 * 
 *     stokes_rhs = 0; 
 * 
 *     const QGauss<dim> quadrature_formula(stokes_degree + 2); 
 *     FEValues<dim>     stokes_fe_values( 
 *       stokes_fe, 
 *       quadrature_formula, 
 *       update_values | update_quadrature_points | update_JxW_values | 
 *         (rebuild_stokes_matrix == true ? update_gradients : UpdateFlags(0))); 
 * 
 *     FEValues<dim> temperature_fe_values(temperature_fe, 
 *                                         quadrature_formula, 
 *                                         update_values); 
 * 
 *     const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来我们需要一个向量，它将包含前一个时间层的温度解在正交点的值，以组装动量方程右侧的源项。让我们把这个向量称为  <code>old_solution_values</code>  。
 * 

 * 
 * 我们接下来创建的向量集包含了基函数的评估以及它们的梯度和对称梯度，将用于创建矩阵。将这些放到自己的数组中，而不是每次都向FEValues对象索取这些信息，是为了加速装配过程的优化，详情请参见 step-22 。
 * 

 * 
 * 最后两个声明是用来从整个FE系统中提取各个块（速度、压力、温度）的。
 * 

 * 
 * 
 * @code
 *     std::vector<double> old_temperature_values(n_q_points); 
 * 
 *     std::vector<Tensor<1, dim>>          phi_u(dofs_per_cell); 
 *     std::vector<SymmetricTensor<2, dim>> grads_phi_u(dofs_per_cell); 
 *     std::vector<double>                  div_phi_u(dofs_per_cell); 
 *     std::vector<double>                  phi_p(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 * @endcode
 * 
 * 现在开始对问题中的所有单元格进行循环。我们正在为这个装配例程处理两个不同的DoFHandlers，所以我们必须为使用中的两个对象设置两个不同的单元格迭代器。这可能看起来有点奇怪，因为斯托克斯系统和温度系统都使用相同的网格，但这是保持自由度同步的唯一方法。循环中的第一条语句也是非常熟悉的，按照更新标志的规定对有限元数据进行更新，将局部数组清零，并在正交点处获得旧解的值。然后我们准备在单元格上的正交点上循环。
 * 

 * 
 * 
 * @code
 *     auto       cell             = stokes_dof_handler.begin_active(); 
 *     const auto endc             = stokes_dof_handler.end(); 
 *     auto       temperature_cell = temperature_dof_handler.begin_active(); 
 * 
 *     for (; cell != endc; ++cell, ++temperature_cell) 
 *       { 
 *         stokes_fe_values.reinit(cell); 
 *         temperature_fe_values.reinit(temperature_cell); 
 * 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 *         temperature_fe_values.get_function_values(old_temperature_solution, 
 *                                                   old_temperature_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             const double old_temperature = old_temperature_values[q]; 
 * 
 * @endcode
 * 
 * 接下来我们提取与内积中的条款相关的基础函数的值和梯度。如 step-22 所示，这有助于加速装配。    一旦完成，我们开始在本地矩阵的行和列上进行循环，并将相关的乘积送入矩阵。右手边是由温度驱动的重力方向（在我们的例子中是垂直方向）的强迫项。 请注意，右手边的项总是生成的，而矩阵的贡献只有在 <code>rebuild_matrices</code> 标志要求时才会更新。
 * 

 * 
 * 
 * @code
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 phi_u[k] = stokes_fe_values[velocities].value(k, q); 
 *                 if (rebuild_stokes_matrix) 
 *                   { 
 *                     grads_phi_u[k] = 
 *                       stokes_fe_values[velocities].symmetric_gradient(k, q); 
 *                     div_phi_u[k] = 
 *                       stokes_fe_values[velocities].divergence(k, q); 
 *                     phi_p[k] = stokes_fe_values[pressure].value(k, q); 
 *                   } 
 *               } 
 * 
 *             if (rebuild_stokes_matrix) 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                   local_matrix(i, j) += 
 *                     (EquationData::eta * 2 * (grads_phi_u[i] * grads_phi_u[j]) - 
 *                      div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
 *                     stokes_fe_values.JxW(q); 
 * 
 *             const Point<dim> gravity = 
 *               -((dim == 2) ? (Point<dim>(0, 1)) : (Point<dim>(0, 0, 1))); 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               local_rhs(i) += (-EquationData::density * EquationData::beta * 
 *                                gravity * phi_u[i] * old_temperature) * 
 *                               stokes_fe_values.JxW(q); 
 *           } 
 * 
 * @endcode
 * 
 * 循环所有单元的最后一步是将局部贡献输入到全局矩阵和向量结构中，并将其输入到  <code>local_dof_indices</code>  指定的位置。 同样，我们让AffineConstraints类来完成将单元格矩阵元素插入全局矩阵的工作，这已经浓缩了悬挂的节点约束。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         if (rebuild_stokes_matrix == true) 
 *           stokes_constraints.distribute_local_to_global(local_matrix, 
 *                                                         local_rhs, 
 *                                                         local_dof_indices, 
 *                                                         stokes_matrix, 
 *                                                         stokes_rhs); 
 *         else 
 *           stokes_constraints.distribute_local_to_global(local_rhs, 
 *                                                         local_dof_indices, 
 *                                                         stokes_rhs); 
 *       } 
 * 
 *     rebuild_stokes_matrix = false; 
 * 
 *     std::cout << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_temperature_matrix"></a> 
 * <h4>BoussinesqFlowProblem::assemble_temperature_matrix</h4>
 * 

 * 
 * 这个函数组装温度方程中的矩阵。温度矩阵由两部分组成，质量矩阵和时间步长乘以刚度矩阵，由拉普拉斯项乘以扩散量给出。由于该矩阵取决于时间步长（从一个步长到另一个步长），温度矩阵需要在每个时间步长进行更新。我们可以简单地在每个时间步长中重新生成矩阵，但这并不真正有效，因为质量和拉普拉斯矩阵只有在我们改变网格时才会改变。因此，我们通过在这个函数中生成两个单独的矩阵，一个是质量矩阵，一个是刚度（扩散）矩阵，这样做更有效率。一旦我们知道了实际的时间步长，我们将把这个矩阵加上刚度矩阵乘以时间步长的总和。
 * 

 * 
 * 所以这第一步的细节非常简单。为了防止我们需要重建矩阵（即网格发生了变化），我们将数据结构归零，得到一个正交公式和一个FEValues对象，并为基函数创建局部矩阵、局部dof指数和评估结构。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::assemble_temperature_matrix() 
 *   { 
 *     if (rebuild_temperature_matrices == false) 
 *       return; 
 * 
 *     temperature_mass_matrix      = 0; 
 *     temperature_stiffness_matrix = 0; 
 * 
 *     QGauss<dim>   quadrature_formula(temperature_degree + 2); 
 *     FEValues<dim> temperature_fe_values(temperature_fe, 
 *                                         quadrature_formula, 
 *                                         update_values | update_gradients | 
 *                                           update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> local_stiffness_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     std::vector<double>         phi_T(dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 现在，让我们开始在三角结构中的所有单元上进行循环。我们需要将局部矩阵清零，更新有限元评估，然后在每个正交点上循环矩阵的行和列，然后我们创建质量矩阵和刚度矩阵（拉普拉斯项乘以扩散  <code>EquationData::kappa</code>  。最后，我们让约束对象将这些值插入全局矩阵中，并直接将约束条件浓缩到矩阵中。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
 *       { 
 *         local_mass_matrix      = 0; 
 *         local_stiffness_matrix = 0; 
 * 
 *         temperature_fe_values.reinit(cell); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 grad_phi_T[k] = temperature_fe_values.shape_grad(k, q); 
 *                 phi_T[k]      = temperature_fe_values.shape_value(k, q); 
 *               } 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   local_mass_matrix(i, j) += 
 *                     (phi_T[i] * phi_T[j] * temperature_fe_values.JxW(q)); 
 *                   local_stiffness_matrix(i, j) += 
 *                     (EquationData::kappa * grad_phi_T[i] * grad_phi_T[j] * 
 *                      temperature_fe_values.JxW(q)); 
 *                 } 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         temperature_constraints.distribute_local_to_global( 
 *           local_mass_matrix, local_dof_indices, temperature_mass_matrix); 
 *         temperature_constraints.distribute_local_to_global( 
 *           local_stiffness_matrix, 
 *           local_dof_indices, 
 *           temperature_stiffness_matrix); 
 *       } 
 * 
 *     rebuild_temperature_matrices = false; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemassemble_temperature_system"></a> 
 * <h4>BoussinesqFlowProblem::assemble_temperature_system</h4>
 * 

 * 
 * 这个函数对温度矩阵进行第二部分的装配工作，实际添加压力质量和刚度矩阵（时间步长在这里起作用），以及创建依赖于速度的右手边。这个函数中的右侧装配的声明与其他装配例程中使用的声明基本相同，只是这次我们把自己限制在矢量上。我们将计算温度系统的残差，这意味着我们必须评估二阶导数，由更新标志 <code>update_hessians</code> 指定。
 * 

 * 
 * 温度方程通过流体速度与斯托克斯系统相耦合。解决方案的这两部分与不同的DoFHandlers相关联，因此我们需要再次创建第二个FEValues对象来评估正交点的速度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::assemble_temperature_system( 
 *     const double maximal_velocity) 
 *   { 
 *     const bool use_bdf2_scheme = (timestep_number != 0); 
 * 
 *     if (use_bdf2_scheme == true) 
 *       { 
 *         temperature_matrix.copy_from(temperature_mass_matrix); 
 *         temperature_matrix *= 
 *           (2 * time_step + old_time_step) / (time_step + old_time_step); 
 *         temperature_matrix.add(time_step, temperature_stiffness_matrix); 
 *       } 
 *     else 
 *       { 
 *         temperature_matrix.copy_from(temperature_mass_matrix); 
 *         temperature_matrix.add(time_step, temperature_stiffness_matrix); 
 *       } 
 * 
 *     temperature_rhs = 0; 
 * 
 *     const QGauss<dim> quadrature_formula(temperature_degree + 2); 
 *     FEValues<dim>     temperature_fe_values(temperature_fe, 
 *                                         quadrature_formula, 
 *                                         update_values | update_gradients | 
 *                                           update_hessians | 
 *                                           update_quadrature_points | 
 *                                           update_JxW_values); 
 *     FEValues<dim>     stokes_fe_values(stokes_fe, 
 *                                    quadrature_formula, 
 *                                    update_values); 
 * 
 *     const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     Vector<double> local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来是向量的声明，用来保存旧的和更早的解决方案的值（分别作为时间级别 $n-1$ 和 $n-2$ 的符号）和当前单元的正交点的梯度。我们还声明了一个对象来保存温度的右侧值（ <code>gamma_values</code> ），并且我们再次使用温度基函数的快捷方式。最终，我们需要找到温度极值和计算域的直径，这将用于稳定参数的定义（我们得到了最大速度作为这个函数的输入）。
 * 

 * 
 * 
 * @code
 *     std::vector<Tensor<1, dim>> old_velocity_values(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_old_velocity_values(n_q_points); 
 *     std::vector<double>         old_temperature_values(n_q_points); 
 *     std::vector<double>         old_old_temperature_values(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_temperature_grads(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_old_temperature_grads(n_q_points); 
 *     std::vector<double>         old_temperature_laplacians(n_q_points); 
 *     std::vector<double>         old_old_temperature_laplacians(n_q_points); 
 * 
 *     EquationData::TemperatureRightHandSide<dim> temperature_right_hand_side; 
 *     std::vector<double>                         gamma_values(n_q_points); 
 * 
 *     std::vector<double>         phi_T(dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell); 
 * 
 *     const std::pair<double, double> global_T_range = 
 *       get_extrapolated_temperature_range(); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 * 
 * @endcode
 * 
 * 现在，让我们开始在三角结构中的所有单元格上进行循环。同样，我们需要两个单元格迭代器，平行走过两个参与的DoFHandler对象的单元格，用于斯托克斯和温度部分。在这个循环中，我们首先将局部rhs设置为零，然后在正交点上获得旧的解函数的值和导数，因为它们将被用于稳定参数的定义和作为方程中的系数，分别需要。请注意，由于温度有自己的DoFHandler和FEValues对象，我们在正交点得到整个解（反正只有标量温度场），而对于斯托克斯部分，我们仅限于通过使用 <code>stokes_fe_values[velocities].get_function_values</code> 提取速度部分（而忽略压力部分）。
 * 

 * 
 * 
 * @code
 *     auto       cell        = temperature_dof_handler.begin_active(); 
 *     const auto endc        = temperature_dof_handler.end(); 
 *     auto       stokes_cell = stokes_dof_handler.begin_active(); 
 * 
 *     for (; cell != endc; ++cell, ++stokes_cell) 
 *       { 
 *         local_rhs = 0; 
 * 
 *         temperature_fe_values.reinit(cell); 
 *         stokes_fe_values.reinit(stokes_cell); 
 * 
 *         temperature_fe_values.get_function_values(old_temperature_solution, 
 *                                                   old_temperature_values); 
 *         temperature_fe_values.get_function_values(old_old_temperature_solution, 
 *                                                   old_old_temperature_values); 
 * 
 *         temperature_fe_values.get_function_gradients(old_temperature_solution, 
 *                                                      old_temperature_grads); 
 *         temperature_fe_values.get_function_gradients( 
 *           old_old_temperature_solution, old_old_temperature_grads); 
 * 
 *         temperature_fe_values.get_function_laplacians( 
 *           old_temperature_solution, old_temperature_laplacians); 
 *         temperature_fe_values.get_function_laplacians( 
 *           old_old_temperature_solution, old_old_temperature_laplacians); 
 * 
 *         temperature_right_hand_side.value_list( 
 *           temperature_fe_values.get_quadrature_points(), gamma_values); 
 * 
 *         stokes_fe_values[velocities].get_function_values(stokes_solution, 
 *                                                          old_velocity_values); 
 *         stokes_fe_values[velocities].get_function_values( 
 *           old_stokes_solution, old_old_velocity_values); 
 * 
 * @endcode
 * 
 * 接下来，我们根据介绍中的讨论，使用专用函数计算用于稳定的人工粘性。有了这个，我们就可以进入正交点和局部rhs矢量分量的循环了。这里的术语相当冗长，但其定义遵循本方案介绍中开发的时间-离散系统。BDF-2方案比用于第一时间步的后向欧拉方案多需要一个旧时间步的术语（并且涉及更复杂的因素）。当所有这些都完成后，我们将局部向量分配到全局向量中（包括悬挂节点约束）。
 * 

 * 
 * 
 * @code
 *         const double nu = 
 *           compute_viscosity(old_temperature_values, 
 *                             old_old_temperature_values, 
 *                             old_temperature_grads, 
 *                             old_old_temperature_grads, 
 *                             old_temperature_laplacians, 
 *                             old_old_temperature_laplacians, 
 *                             old_velocity_values, 
 *                             old_old_velocity_values, 
 *                             gamma_values, 
 *                             maximal_velocity, 
 *                             global_T_range.second - global_T_range.first, 
 *                             cell->diameter()); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 grad_phi_T[k] = temperature_fe_values.shape_grad(k, q); 
 *                 phi_T[k]      = temperature_fe_values.shape_value(k, q); 
 *               } 
 * 
 *             const double T_term_for_rhs = 
 *               (use_bdf2_scheme ? 
 *                  (old_temperature_values[q] * (1 + time_step / old_time_step) - 
 *                   old_old_temperature_values[q] * (time_step * time_step) / 
 *                     (old_time_step * (time_step + old_time_step))) : 
 *                  old_temperature_values[q]); 
 * 
 *             const Tensor<1, dim> ext_grad_T = 
 *               (use_bdf2_scheme ? 
 *                  (old_temperature_grads[q] * (1 + time_step / old_time_step) - 
 *                   old_old_temperature_grads[q] * time_step / old_time_step) : 
 *                  old_temperature_grads[q]); 
 * 
 *             const Tensor<1, dim> extrapolated_u = 
 *               (use_bdf2_scheme ? 
 *                  (old_velocity_values[q] * (1 + time_step / old_time_step) - 
 *                   old_old_velocity_values[q] * time_step / old_time_step) : 
 *                  old_velocity_values[q]); 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               local_rhs(i) += 
 *                 (T_term_for_rhs * phi_T[i] - 
 *                  time_step * extrapolated_u * ext_grad_T * phi_T[i] - 
 *                  time_step * nu * ext_grad_T * grad_phi_T[i] + 
 *                  time_step * gamma_values[q] * phi_T[i]) * 
 *                 temperature_fe_values.JxW(q); 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         temperature_constraints.distribute_local_to_global(local_rhs, 
 *                                                            local_dof_indices, 
 *                                                            temperature_rhs); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemsolve"></a> 
 * <h4>BoussinesqFlowProblem::solve</h4>
 * 

 * 
 * 这个函数可以解决线性方程组的问题。在介绍之后，我们从斯托克斯系统开始，在这里我们需要生成我们的块状舒尔预处理器。由于所有相关的动作都在类 <code>BlockSchurPreconditioner</code> 中实现，我们所要做的就是适当地初始化这个类。我们需要传递的是一个用于压力质量矩阵的 <code>InverseMatrix</code> 对象，我们使用相应的类和我们已经生成的IC预处理器以及用于速度-速度矩阵的AMG预处理器一起设置。注意， <code>Mp_preconditioner</code> 和 <code>Amg_preconditioner</code> 都只是指针，所以我们用 <code>*</code> 来传递实际的预处理对象。
 * 

 * 
 * 一旦预处理程序准备好了，我们就为该块系统创建一个GMRES求解器。由于我们使用的是Trilinos数据结构，我们必须在求解器中设置相应的模板参数。GMRES需要在内部存储每次迭代的临时向量（见 step-22 的结果部分的讨论）&ndash；它可以使用的向量越多，一般来说性能越好。为了控制内存需求，我们将向量的数量设置为100。这意味着在求解器的100次迭代中，每个临时向量都可以被存储。如果求解器需要更频繁地迭代以获得指定的容忍度，它将通过每100次迭代重新开始，在一个减少的向量集上工作。
 * 

 * 
 * 有了这些设置，我们求解系统并在斯托克斯系统中分配约束条件，即悬挂节点和无流体边界条件，以便即使在受约束的道夫下也有适当的解值。最后，我们把迭代次数写到屏幕上。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::solve() 
 *   { 
 *     std::cout << "   Solving..." << std::endl; 
 * 
 *     { 
 *       const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                                          TrilinosWrappers::PreconditionIC> 
 *         mp_inverse(stokes_preconditioner_matrix.block(1, 1), 
 *                    *Mp_preconditioner); 
 * 
 *       const LinearSolvers::BlockSchurPreconditioner< 
 *         TrilinosWrappers::PreconditionAMG, 
 *         TrilinosWrappers::PreconditionIC> 
 *         preconditioner(stokes_matrix, mp_inverse, *Amg_preconditioner); 
 * 
 *       SolverControl solver_control(stokes_matrix.m(), 
 *                                    1e-6 * stokes_rhs.l2_norm()); 
 * 
 *       SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres( 
 *         solver_control, 
 *         SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(100)); 
 * 
 *       for (unsigned int i = 0; i < stokes_solution.size(); ++i) 
 *         if (stokes_constraints.is_constrained(i)) 
 *           stokes_solution(i) = 0; 
 * 
 *       gmres.solve(stokes_matrix, stokes_solution, stokes_rhs, preconditioner); 
 * 
 *       stokes_constraints.distribute(stokes_solution); 
 * 
 *       std::cout << "   " << solver_control.last_step() 
 *                 << " GMRES iterations for Stokes subsystem." << std::endl; 
 *     } 
 * 
 * @endcode
 * 
 * 一旦我们知道了斯托克斯解，我们就可以根据最大速度确定新的时间步长。我们必须这样做以满足CFL条件，因为对流项在温度方程中得到了明确的处理，正如在介绍中所讨论的那样。这里使用的时间步长公式的确切形式将在本程序的结果部分讨论。
 * 

 * 
 * 这里有一个插曲。该公式包含了对速度最大值的除法。然而，在计算开始时，我们有一个恒定的温度场（我们以恒定的温度开始，只有在源作用的第一个时间步长后，它才会变成非恒定的）。恒定温度意味着没有浮力作用，所以速度为零。除以它不可能得到什么好结果。
 * 

 * 
 * 为了避免产生无限的时间步长，我们问最大速度是否非常小（特别是小于我们在接下来的任何时间步长中遇到的值），如果是，我们就不除以零，而是除以一个小值，从而产生一个大的但有限的时间步长。
 * 

 * 
 * 
 * @code
 *     old_time_step                 = time_step; 
 *     const double maximal_velocity = get_maximal_velocity(); 
 * 
 *     if (maximal_velocity >= 0.01) 
 *       time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree * 
 *                   GridTools::minimal_cell_diameter(triangulation) / 
 *                   maximal_velocity; 
 *     else 
 *       time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree * 
 *                   GridTools::minimal_cell_diameter(triangulation) / .01; 
 * 
 *     std::cout << "   " 
 *               << "Time step: " << time_step << std::endl; 
 * 
 *     temperature_solution = old_temperature_solution; 
 * 
 * @endcode
 * 
 * 接下来我们用函数  <code>assemble_temperature_system()</code>  设置温度系统和右手边。 知道了温度方程的矩阵和右手边，我们设置了一个预处理程序和一个求解器。温度矩阵是一个质量矩阵（特征值在1左右）加上一个拉普拉斯矩阵（特征值在0和 $ch^{-2}$ 之间）乘以一个与时间步长成正比的小数字  $k_n$  。因此，产生的对称和正定矩阵的特征值在 $[1,1+k_nh^{-2}]$ 范围内（至于常数）。这个矩阵即使对于小的网格尺寸也只是适度的条件不良，我们通过简单的方法得到一个相当好的预处理，例如用一个不完全的Cholesky分解预处理（IC），我们也用它来预处理压力质量矩阵求解器。作为一个求解器，我们选择共轭梯度法CG。和以前一样，我们通过模板参数 <code>TrilinosWrappers::MPI::Vector</code> 告诉求解器使用Trilinos向量。最后，我们求解，分配悬挂节点约束，并写出迭代次数。
 * 

 * 
 * 
 * @code
 *     assemble_temperature_system(maximal_velocity); 
 *     { 
 *       SolverControl solver_control(temperature_matrix.m(), 
 *                                    1e-8 * temperature_rhs.l2_norm()); 
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 
 * 
 *       TrilinosWrappers::PreconditionIC preconditioner; 
 *       preconditioner.initialize(temperature_matrix); 
 * 
 *       cg.solve(temperature_matrix, 
 *                temperature_solution, 
 *                temperature_rhs, 
 *                preconditioner); 
 * 
 *       temperature_constraints.distribute(temperature_solution); 
 * 
 *       std::cout << "   " << solver_control.last_step() 
 *                 << " CG iterations for temperature." << std::endl; 
 * 
 * @endcode
 * 
 * 在这个函数的结尾，我们在向量中步进并读出最大和最小的温度值，我们也想输出这些值。在本程序的结果部分讨论的确定时间步长的正确常数时，这将非常有用。
 * 

 * 
 * 
 * @code
 *       double min_temperature = temperature_solution(0), 
 *              max_temperature = temperature_solution(0); 
 *       for (unsigned int i = 0; i < temperature_solution.size(); ++i) 
 *         { 
 *           min_temperature = 
 *             std::min<double>(min_temperature, temperature_solution(i)); 
 *           max_temperature = 
 *             std::max<double>(max_temperature, temperature_solution(i)); 
 *         } 
 * 
 *       std::cout << "   Temperature range: " << min_temperature << ' ' 
 *                 << max_temperature << std::endl; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemoutput_results"></a> 
 * <h4>BoussinesqFlowProblem::output_results</h4>
 * 

 * 
 * 该函数将解决方案写入VTK输出文件，用于可视化，每隔10个时间步长就会完成。这通常是一个相当简单的任务，因为deal.II库提供的函数几乎为我们完成了所有的工作。与以前的例子相比，有一个新的函数。我们想把斯托克斯解和温度都看作一个数据集，但是我们已经根据两个不同的DoFHandler对象完成了所有的计算。幸运的是，DataOut类已经准备好处理这个问题。我们所要做的就是不要在一开始就附加一个单一的DoFHandler，然后将其用于所有添加的向量，而是为每个向量分别指定DoFHandler。剩下的就像  step-22  中所做的那样。我们创建解决方案的名称（这些名称将出现在各个组件的可视化程序中）。第一个 <code>dim</code> 分量是矢量速度，然后我们有斯托克斯部分的压力，而温度是标量。这些信息是用DataComponentInterpretation辅助类读出来的。接下来，我们将数据向量与它们的DoFHandler对象连接起来，根据自由度建立补丁，这些补丁是描述可视化程序数据的（子）元素。最后，我们打开一个文件（包括时间步数）并将vtk数据写入其中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::output_results() const 
 *   { 
 *     if (timestep_number % 10 != 0) 
 *       return; 
 * 
 *     std::vector<std::string> stokes_names(dim, "velocity"); 
 *     stokes_names.emplace_back("p"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       stokes_component_interpretation( 
 *         dim + 1, DataComponentInterpretation::component_is_scalar); 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       stokes_component_interpretation[i] = 
 *         DataComponentInterpretation::component_is_part_of_vector; 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.add_data_vector(stokes_dof_handler, 
 *                              stokes_solution, 
 *                              stokes_names, 
 *                              stokes_component_interpretation); 
 *     data_out.add_data_vector(temperature_dof_handler, 
 *                              temperature_solution, 
 *                              "T"); 
 *     data_out.build_patches(std::min(stokes_degree, temperature_degree)); 
 * 
 *     std::ofstream output("solution-" + 
 *                          Utilities::int_to_string(timestep_number, 4) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemrefine_mesh"></a> 
 * <h4>BoussinesqFlowProblem::refine_mesh</h4>
 * 

 * 
 * 这个函数负责处理自适应网格细化。这个函数执行的三个任务是：首先找出需要细化/粗化的单元，然后实际进行细化，并最终在两个不同的网格之间传输解向量。第一个任务是通过对温度使用成熟的凯利误差估计器来实现的（对于这个程序，我们主要关注的是温度，我们需要在高温度梯度的区域保持精确，同时也要避免有太多的数值扩散）。第二项任务是实际进行再塑形。这也只涉及到基本函数，例如 <code>refine_and_coarsen_fixed_fraction</code> ，它可以细化那些具有最大估计误差的单元，这些误差合计占80%，并粗化那些具有最小误差的单元，这些误差合计占10%。
 * 

 * 
 * 如果像这样实施，我们会得到一个不会有太大进展的程序。请记住，我们期望的温度场几乎是不连续的（扩散率 $\kappa$ 毕竟非常小），因此我们可以预期，一个自由适应的网格会越来越细化到大梯度的区域。网格大小的减少将伴随着时间步长的减少，需要大量的时间步长来解决给定的最终时间。这也会导致在几个网格细化周期后，网格的不连续性解决得比开始时好得多。
 * 

 * 
 * 特别是为了防止时间步长的减少和相应的大量时间步长，我们限制了网格的最大细化深度。为此，在细化指标应用于单元格后，我们简单地在最细层的所有单元格上循环，如果它们会导致网格层次过高，则取消对它们的细化选择。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   BoussinesqFlowProblem<dim>::refine_mesh(const unsigned int max_grid_level) 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     KellyErrorEstimator<dim>::estimate(temperature_dof_handler, 
 *                                        QGauss<dim - 1>(temperature_degree + 1), 
 *                                        {}, 
 *                                        temperature_solution, 
 *                                        estimated_error_per_cell); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, 
 *                                                       estimated_error_per_cell, 
 *                                                       0.8, 
 *                                                       0.1); 
 *     if (triangulation.n_levels() > max_grid_level) 
 *       for (auto &cell : 
 *            triangulation.active_cell_iterators_on_level(max_grid_level)) 
 *         cell->clear_refine_flag(); 
 * 
 * @endcode
 * 
 * 作为网格细化的一部分，我们需要将旧的网格中的解决方案向量转移到新的网格中。为此，我们使用SolutionTransfer类，我们必须准备好需要转移到新网格的解向量（一旦完成细化，我们将失去旧的网格，所以转移必须与细化同时发生）。我们肯定需要的是当前温度和旧温度（BDF-2时间步长需要两个旧的解决方案）。由于SolutionTransfer对象只支持在每个dof处理程序中传输一个对象，我们需要在一个数据结构中收集两个温度解决方案。此外，我们也选择转移斯托克斯解，因为我们需要前两个时间步长的速度，其中只有一个是在飞行中计算的。
 * 

 * 
 * 因此，我们为斯托克斯和温度的DoFHandler对象初始化了两个SolutionTransfer对象，将它们附加到旧的dof处理程序中。有了这个，我们就可以准备三角测量和数据向量的细化了（按这个顺序）。
 * 

 * 
 * 
 * @code
 *     std::vector<TrilinosWrappers::MPI::Vector> x_temperature(2); 
 *     x_temperature[0]                            = temperature_solution; 
 *     x_temperature[1]                            = old_temperature_solution; 
 *     TrilinosWrappers::MPI::BlockVector x_stokes = stokes_solution; 
 * 
 *     SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> temperature_trans( 
 *       temperature_dof_handler); 
 *     SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> stokes_trans( 
 *       stokes_dof_handler); 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 *     temperature_trans.prepare_for_coarsening_and_refinement(x_temperature); 
 *     stokes_trans.prepare_for_coarsening_and_refinement(x_stokes); 
 * 
 * @endcode
 * 
 * 现在一切都准备好了，所以进行细化，在新的网格上重新创建dof结构，并初始化矩阵结构和 <code>setup_dofs</code> 函数中的新向量。接下来，我们实际执行网格之间的插值解。我们为温度创建另一份临时向量（现在与新网格相对应），并让插值函数完成这项工作。然后，产生的向量数组被写入各自的向量成员变量中。
 * 

 * 
 * 记住，约束集将在setup_dofs()调用中为新的三角结构进行更新。
 * 

 * 
 * 
 * @code
 *     triangulation.execute_coarsening_and_refinement(); 
 *     setup_dofs(); 
 * 
 *     std::vector<TrilinosWrappers::MPI::Vector> tmp(2); 
 *     tmp[0].reinit(temperature_solution); 
 *     tmp[1].reinit(temperature_solution); 
 *     temperature_trans.interpolate(x_temperature, tmp); 
 * 
 *     temperature_solution     = tmp[0]; 
 *     old_temperature_solution = tmp[1]; 
 * 
 * @endcode
 * 
 * 在解决方案被转移后，我们再对被转移的解决方案实施约束。
 * 

 * 
 * 
 * @code
 *     temperature_constraints.distribute(temperature_solution); 
 *     temperature_constraints.distribute(old_temperature_solution); 
 * 
 * @endcode
 * 
 * 对于斯托克斯矢量，一切都一样&ndash;除了我们不需要另一个临时矢量，因为我们只是插值了一个矢量。最后，我们必须告诉程序，矩阵和预处理程序需要重新生成，因为网格已经改变。
 * 

 * 
 * 
 * @code
 *     stokes_trans.interpolate(x_stokes, stokes_solution); 
 * 
 *     stokes_constraints.distribute(stokes_solution); 
 * 
 *     rebuild_stokes_matrix         = true; 
 *     rebuild_temperature_matrices  = true; 
 *     rebuild_stokes_preconditioner = true; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BoussinesqFlowProblemrun"></a> 
 * <h4>BoussinesqFlowProblem::run</h4>
 * 

 * 
 * 这个函数执行Boussinesq程序中的所有基本步骤。它首先设置一个网格（根据空间维度，我们选择一些不同级别的初始细化和额外的自适应细化步骤，然后在 <code>dim</code> 维度上创建一个立方体，并首次设置了道夫。由于我们想用一个自适应细化的网格开始时间步进，我们执行一些预细化步骤，包括所有的装配、求解和细化，但实际上没有在时间上推进。相反，我们使用被人诟病的 <code>goto</code> 语句，在网格细化后立即跳出时间循环，从 <code>start_time_iteration</code> 标签开始的新网格上重新开始。( <code>goto</code> 的使用将在 step-26 中讨论) 。
 * 

 * 
 * 在我们开始之前，我们将初始值投影到网格上，并获得 <code>old_temperature_solution</code> 矢量的第一个数据。然后，我们初始化时间步数和时间步长，开始时间循环。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoussinesqFlowProblem<dim>::run() 
 *   { 
 *     const unsigned int initial_refinement     = (dim == 2 ? 4 : 2); 
 *     const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3); 
 * 
 *     GridGenerator::hyper_cube(triangulation); 
 *     global_Omega_diameter = GridTools::diameter(triangulation); 
 * 
 *     triangulation.refine_global(initial_refinement); 
 * 
 *     setup_dofs(); 
 * 
 *     unsigned int pre_refinement_step = 0; 
 * 
 *   start_time_iteration: 
 * 
 *     VectorTools::project(temperature_dof_handler, 
 *                          temperature_constraints, 
 *                          QGauss<dim>(temperature_degree + 2), 
 *                          EquationData::TemperatureInitialValues<dim>(), 
 *                          old_temperature_solution); 
 * 
 *     timestep_number = 0; 
 *     time_step = old_time_step = 0; 
 * 
 *     double time = 0; 
 * 
 *     do 
 *       { 
 *         std::cout << "Timestep " << timestep_number << ":  t=" << time 
 *                   << std::endl; 
 * 
 * @endcode
 * 
 * 时间循环的第一步都是显而易见的；我们组装斯托克斯系统、预处理程序、温度矩阵（矩阵和预处理程序实际上只在我们之前重新处理的情况下发生变化），然后进行求解。在继续下一个时间步骤之前，我们必须检查我们是否应该首先完成预精炼步骤，或者是否应该重新啮合（每五个时间步骤），精炼到一个与初始精炼和预精炼步骤一致的水平。循环的最后一个步骤是推进解，即把解复制到下一个 "较早 "的时间层。
 * 

 * 
 * 
 * @code
 *         assemble_stokes_system(); 
 *         build_stokes_preconditioner(); 
 *         assemble_temperature_matrix(); 
 * 
 *         solve(); 
 * 
 *         output_results(); 
 * 
 *         std::cout << std::endl; 
 * 
 *         if ((timestep_number == 0) && 
 *             (pre_refinement_step < n_pre_refinement_steps)) 
 *           { 
 *             refine_mesh(initial_refinement + n_pre_refinement_steps); 
 *             ++pre_refinement_step; 
 *             goto start_time_iteration; 
 *           } 
 *         else if ((timestep_number > 0) && (timestep_number % 5 == 0)) 
 *           refine_mesh(initial_refinement + n_pre_refinement_steps); 
 * 
 *         time += time_step; 
 *         ++timestep_number; 
 * 
 *         old_stokes_solution          = stokes_solution; 
 *         old_old_temperature_solution = old_temperature_solution; 
 *         old_temperature_solution     = temperature_solution; 
 *       } 
 * 
 * @endcode
 * 
 * 做以上所有的工作，直到我们到达时间100。
 * 

 * 
 * 
 * @code
 *     while (time <= 100); 
 *   } 
 * } // namespace Step31 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 主函数看起来与所有其他程序几乎一样。
 * 

 * 
 * 有一个区别是我们必须要注意的。这个程序使用了Trilinos，而通常情况下，Trilinos被配置为可以使用MPI在%parallel中运行。这并不意味着它<i>has</i>可以在%parallel中运行，事实上这个程序（不像 step-32 ）根本没有尝试使用MPI在%parallel中做任何事情。然而，Trilinos希望MPI系统被初始化。我们通过创建一个类型为 Utilities::MPI::MPI_InitFinalize 的对象来做到这一点，该对象使用给main()的参数（即 <code>argc</code> 和 <code>argv</code> ）初始化MPI（如果可用的话），并在对象超出范围时再次去初始化它。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step31; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization( 
 *         argc, argv, numbers::invalid_unsigned_int); 
 * 
 * @endcode
 * 
 * 这个程序只能在串行中运行。否则，将抛出一个异常。
 * 

 * 
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
 *                   ExcMessage( 
 *                     "This program can only be run in serial, use ./step-31")); 
 * 
 *       BoussinesqFlowProblem<2> flow_problem; 
 *       flow_problem.run(); 
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
 * 
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
 * 
 * @endcode
examples/step-31/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Resultsin2d"></a><h3> Results in 2d </h3>


当你在2D中运行该程序时，输出将看起来像这样。<code> <pre> 活动单元的数量：256（在5层） 自由度的数量：3556（2178+289+1089)

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.919118温度的9次CG迭代。    温度范围：-0.16687 1.30011

活动单元的数量：280（在6层） 自由度的数量：4062（2490+327+1245）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.459559温度的9次CG迭代。    温度范围：-0.0982971 0.598503

活动单元的数量：520（在7个层面上） 自由度的数量：7432（4562+589+2281）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.229779 温度的9次CG迭代。    温度范围：-0.0551098 0.294493

活动单元的数量：1072（在8层） 自由度的数量：15294（9398+1197+4699）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.11489 温度的9次CG迭代。    温度范围：-0.0273524 0.156861

活动单元的数量：2116（在9层） 自由度的数量：30114（18518+2337+9259）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.0574449温度的9次CG迭代。    温度范围：-0.014993 0.0738328

时间步骤1：t=0.0574449 装配...    解决...    斯托克斯子系统的56次GMRES迭代。    时间步长：0.0574449 温度的9次CG迭代。    温度范围：-0.0273934 0.14488

...</pre> </code>

在开始的时候，我们自适应地细化了几次网格，并总是返回到时间步长为零的新细化的网格上重新开始。只有这样，我们才开始实际的时间迭代。

程序运行了一段时间。时间步数为0、500、1000、1500、2000、3000、4000和5000的温度字段看起来是这样的（注意温度使用的色标并不总是相同）。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.solution.07.png" alt="">
    </td>
  </tr>
</table> 

这里显示的视觉效果是使用实例的一个版本生成的，该版本在传输网格后没有强制执行约束。

可以看出，我们有三个加热流体的热源，因此产生了一个浮力效应，让流体的热袋上升并旋转起来。通过烟囱效应，这三股气流被来自外部并想加入上升气流的流体压在一起。请注意，由于流体最初处于静止状态，那些最初在源头上的流体部分比后来被充分发展的流场拖到源头上的流体获得更长的加热时间。因此，它更热，这一事实可以从三个羽流的红色尖端看出。还要注意流场的相对精细的特征，这是我们选择的温度方程的复杂传输稳定的结果。

除了上面的图片外，下面的图片显示了自适应网格和同一时间步长的流场。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.2d.grid.07.png" alt="">
    </td>
  </tr>
</table> 




<a name="Resultsin3d"></a><h3> Results in 3d </h3>


当然，同样的事情也可以在3D中完成，将 <code>main()</code> 中的BoussinesqFlowProblem对象的模板参数从2改为3，这样，现在的输出看起来如下。

<code> <pre> 活动单元的数量：64（在3层） 自由度的数量：3041（2187+125+729）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：2.45098 温度的9次CG迭代。    温度范围：-0.675683 4.94725

活动单元的数量：288（在4层） 自由度的数量：12379（8943+455+2981）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：1.22549 温度的9次CG迭代。    温度范围：-0.527701 2.25764

活动单元的数量：1296（在5层） 自由度的数量：51497（37305+1757+12435）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.612745温度的10次CG迭代。    温度范围：-0.496942 0.847395

活动单元的数量：5048（在6层） 自由度的数量：192425（139569+6333+46523）。

时间步数0: t=0 正在组装...    重建斯托克斯预处理程序...    解算...    0次GMRES迭代，用于斯托克斯子系统。    时间步长：0.306373 温度的10次CG迭代。    温度范围：-0.267683 0.497739

时间步数1：t=0.306373 正在组装...    解决...    斯托克斯子系统的27次GMRES迭代。    时间步长：0.306373 温度的10次CG迭代。    温度范围：-0.461787 0.958679

...</pre> </code>

在时间步数为0、50、100、150、200、300、400、500、600、700和800的情况下，将温度等值线可视化，得到以下图示。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.00.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.02.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.04.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.05.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.06.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.07.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.08.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.09.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.3d.solution.10.png" alt="">
    </td>
    <td>
    </td>
  </tr>
</table> 

第一张图片看起来像三只刺猬，这是因为我们的方案基本上是将源乘以第一时间步长投射到网格上，以获得第一时间步的温度场。由于源函数是不连续的，我们需要期待这个项目的过冲和欠冲。这就是事实上发生的情况（在2d中更容易检查），并导致等值面的皱缩外观。  这里显示的视觉效果是使用例子的一个版本生成的，该版本在传输网格后没有强制执行约束。




<a name="Numericalexperimentstodetermineoptimalparameters"></a><h3> Numerical experiments to determine optimal parameters </h3>


现在的程序有三个参数，我们在理论上并没有掌握如何以最佳方式进行选择。这三个参数是。   <ul>   <li>  时间步骤必须满足CFL条件  $k\le \min_K \frac{c_kh_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。这里， $c_k$ 是无量纲的，但什么是正确的值？     <li>  在计算人工黏度时。

@f{eqnarray*}
  \nu_\alpha(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)}
  \min\left\{
    h_K,
    h_K^\alpha
    \frac{\|R_\alpha(T)\|_{L^\infty(K)}}{c(\mathbf{u},T)}
  \right\},


@f}

      与 $c(\mathbf{u},T) =
      c_R\ \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
      \ |\mathrm{diam}(\Omega)|^{\alpha-2}$  。       这里，无量纲%数 $\beta,c_R$ 的选择是有意义的。   </ul>  在所有这些情况下，我们将不得不期望每个值的正确选择取决于其他值的正确选择，而且很可能也取决于用于温度的有限元的空间尺寸和多项式程度。下面我们将讨论一些数值实验来选择常数  $c_k$  和  $\beta$  。

下面，我们将不讨论 $c_R$ 的选择问题。在程序中，我们将其设定为 $c_R=2^{\frac{4-2\alpha}{d}}$  。这个值的原因有点复杂，与程序的历史而不是推理有关：虽然全局缩放参数 $c(\mathbf{u},T)$ 的正确公式如上所示，但程序（包括与deal.II 6.2一起出厂的版本）最初有一个错误，即我们计算的是 $c(\mathbf{u},T) =
      \|\mathbf{u}\|_{L^\infty(\Omega)} \ \mathrm{var}(T)
      \ \frac{1}{|\mathrm{diam}(\Omega)|^{\alpha-2}}$ ，而在这里我们将缩放参数设置为1。由于我们只在 $\mathrm{diam}(\Omega)=2^{1/d}$ 的单位平方/立方体上进行计算，这完全等同于使用 $c_R=\left(2^{1/d}\right)^{4-2\alpha}=2^{\frac{4-2\alpha}{d}}$ 的正确公式。由于 $c_R$ 的这个值对于当前的程序来说似乎很好用，我们在程序中修正了公式，并将 $c_R$ 设置为一个值，正好再现了我们之前的结果。不过，我们将在第32步中再次审视这个问题。

然而，现在回到讨论 $c_k$ 和 $\beta$ 的什么值来选择。




<a name="Choosingicsubksubiicsubksubiandbeta"></a><h4> Choosing <i>c<sub>k</sub></i><i>c<sub>k</sub></i> and beta </h4> 。


这两个常数肯定在某种程度上有联系。原因很容易看出来。在纯平流问题的情况下， $\frac{\partial T}{\partial t} + \mathbf{u}\cdot\nabla T = \gamma$ ，任何显式方案都必须满足形式为 $k\le \min_K \frac{c_k^a h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 的CFL条件。另一方面，对于纯扩散问题， $\frac{\partial T}{\partial t} + \nu \Delta T = \gamma$ ，显式方案需要满足一个条件 $k\le \min_K \frac{c_k^d h_K^2}{\nu}$ 。因此，鉴于上述 $\nu$ 的形式，像我们这里要解决的平流扩散问题将导致一个 $
k\le \min_K \min \left\{
  \frac{c_k^a h_K}{\|\mathbf{u}\|_{L^\infty(K)}},
  \frac{c_k^d h_K^2}{\beta \|\mathbf{u}\|_{L^\infty(K)} h_K}\right\}
  =
  \min_K \left( \min \left\{
  c_k^a,
  \frac{c_k^d}{\beta}\right\}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}} \right)
$ 的条件。因此，我们必须面对这样一个事实：我们可能想选择 $\beta$ 大一些，以提高数值方案的稳定性（通过增加人工扩散量），但我们必须以更小的、因而更多的时间步骤为代价。因此，在实践中，人们希望尽可能地选择 $\beta$ ，以保持传输问题的充分稳定，同时尽量选择大的时间步长，以减少总体工作量。

要找到正确的平衡，唯一的办法是做一些计算实验。下面是我们的做法。我们稍微修改了程序，允许更少的网格细化（所以我们不一定要等那么久），并选择 $
  \nu(T)|_K
  =
  \beta
  \|\mathbf{u}\|_{L^\infty(K)} h_K
$ 来消除常数 $c_R$ 的影响（我们知道通过使用这个版本的 $\nu(T)$ 作为人工粘度，解决方案是稳定的，但我们可以通过使用这个人工粘度的更复杂的公式来改善情况--即使解决方案更清晰）。然后我们对不同的值 $c_k,\beta$ 运行程序，观察域中的最大和最小温度。我们期望看到的情况是这样的。如果我们选择的时间步长过大（即选择一个比理论上允许的大的 $c_k$ ），那么我们将得到温度的指数式增长。如果我们选择 $\beta$ 太小，那么传输稳定变得不充分，解决方案将显示出明显的振荡，但不是指数级增长。




<a name="ResultsforQsub1subelements"></a><h5>Results for Q<sub>1</sub> elements</h5>


下面是我们对 $\beta=0.01, \beta=0.1$ ，和 $\beta=0.5$ ， $c_k$ 的不同选择，以及2d的双线性元素（ <code>temperature_degree=1</code> ）得到的结果。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.1.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q1.beta=0.5.png" alt="">
    </td>
  </tr>
</table> 

解释这些图表的方法是这样的：对于 $\beta=0.01$ 和 $c_k=\frac 12,\frac 14$ ，我们看到指数增长或至少是大的变化，但如果我们选择 $k=\frac 18\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小，那么这个方案虽然有点摇摆不定，但还是稳定的。对于更多的人工扩散，我们可以选择 $k=\frac 14\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小的 $\beta=0.03$ ， $k=\frac 13\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 或更小的 $\beta=0.1$ ，并再次需要 $k=\frac 1{15}\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 的 $\beta=0.5$ （这次是因为许多扩散需要一个小的时间步长）。

那么该如何选择呢？如果我们只是对大时间步长感兴趣，那么我们会选择 $\beta=0.1$ 和 $k=\frac 13\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。另一方面，我们也对准确性感兴趣，在这里，实际调查这些曲线所显示的内容可能会有兴趣。为此，请注意，我们从零温度开始，我们的来源是正的&mdash；所以我们会直观地期望温度永远不会降到零以下。但它确实如此，这是使用连续元素来近似不连续的解决方案时，吉布现象的结果。因此，我们可以看到，选择 $\beta$ 太小是不好的：太少的人工扩散会导致没有扩散掉的过冲和欠冲。另一方面，对于大的 $\beta$ ，最低温度在开始时下降到零以下，但随后迅速扩散回零。

另一方面，我们也来看看最高温度。观察溶液的电影，我们看到最初流体处于静止状态。源头不断加热相同体积的流体，其温度在开始时呈线性增长，直到其浮力能够使其向上移动。因此，流体中最热的部分被带离了溶液，取而代之的流体只被加热了很短的时间就被移出了源区，因此仍然比初始气泡要冷。如果 $\kappa=0$ （在程序中是非零的，但非常小），那么流体中最热的部分应该随着流动而平移，其温度不变。这就是我们在最小的 $\beta$ 图中可以看到的：一旦达到最高温度，它就几乎不再变化。另一方面，人工扩散越大，热点的扩散就越多。请注意，对于这个标准，时间步长的大小并不发挥重要作用。

因此，总结起来，可能最好的选择似乎是 $\beta=0.03$ 和 $k=\frac 14\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$  。曲线有点摇摆不定，但总的来说，图片看起来相当合理，除了由于吉布现象而在接近开始时间时出现一些过冲和欠冲的情况。




<a name="ResultsforQsub2subelements"></a><h5>Results for Q<sub>2</sub> elements</h5>


我们也可以对高阶元素重复同样的实验序列。这里是温度的双二次方形状函数（ <code>temperature_degree=2</code> ）的图形，同时我们保留了斯托克斯系统的 $Q_2/Q_1$ 稳定泰勒-胡德元素。

 <table align="center" class="doxtable">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.01.png" alt="">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-31.timestep.q2.beta=0.1.png" alt="">
    </td>
  </tr>
</table> 

同样， $\beta$ 的小值会导致较少的扩散，但我们必须选择非常小的时间步长来保持事情的控制。太大的 $\beta$ 值会导致更多的扩散，但同样需要小的时间步骤。最佳值似乎是 $\beta=0.03$ ，和 $Q_1$ 元素一样，然后我们必须选择 $k=\frac 18\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ &mdash；正好是 $Q_1$ 元素的一半大小。]元素，如果我们把CFL条件说成是要求时间步长足够小，以便运输在每个时间步长中的移动距离不超过一个<i>grid point</i>距离（对于 $Q_1$ 元素是 $h_K$ ，但对于 $Q_2$ 元素是 $h_K/2$ ），这个事实可能并不令人惊讶。事实证明， $\beta$ 需要稍微大一点，以便在模拟后期获得稳定的结果，时间大于60，所以我们实际上在代码中选择它作为 $\beta = 0.034$ 。




<a name="Resultsfor3d"></a><h5>Results for 3d</h5>


我们可以在3D中重复这些实验，找到每个 $\beta$ 值的最佳时间步骤，并找到 $\beta$ 的最佳值。人们发现，对于2d中已经使用的相同的 $\beta$ ，时间步长需要小一点，大约是1.2倍左右。这很容易解释：时间步长的限制是 $k=\min_K \frac{ch_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ ，其中 $h_K$ 是单元的<i>diameter</i>。然而，真正需要的是网格点之间的距离，它是 $\frac{h_K}{\sqrt{d}}$  。所以更合适的形式是  $k=\min_K \frac{ch_K}{\|\mathbf{u}\|_{L^\infty(K)}\sqrt{d}}$  。

第二个发现是，需要把 $\beta$ 选得稍微大一点（大约 $\beta=0.05$ 左右）。这就再次减少了我们可以采取的时间步骤。







<a name="Conclusions"></a><h5>Conclusions</h5>


总之，从上面的简单计算来看， $\beta=0.034$ 似乎是2D中稳定参数的一个好选择，而 $\beta=0.05$ 则是3D中的稳定参数。以独立于维度的方式，我们可以将其建模为 $\beta=0.017d$  。如果在更细的网格上做更长时间的计算（几千个时间步长），就会意识到时间步长还不够小，为了稳定，就必须把上述数值再降低一些（大约是 $\frac 78$ 的一个系数）。

因此，调和2D、3D和可变多项式程度并考虑到所有因素的公式如下。

@f{eqnarray*}
  k =
  \frac 1{2 \cdot 1.7} \frac 1{\sqrt{d}}
  \frac 2d
  \frac 1{q_T}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}
  =
  \frac 1{1.7 d\sqrt{d}}
  \frac 1{q_T}
  \frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}.


@f}

在第一种形式中（方程中心）， $\frac
1{2 \cdot 1.7}$ 是一个通用常数， $\frac 1{\sqrt{d}}$ 是说明单元直径和网格点间距的因素， $\frac 2d$ 说明 $\beta$ 随着空间尺寸的增加而增加， $\frac 1{q_T}$ 说明高阶元素的网格点之间的距离， $\frac{h_K}{\|\mathbf{u}\|_{L^\infty(K)}}$ 说明相对于单元尺寸的局部传输速度。这就是我们在程序中使用的公式。

至于对温度使用 $Q_1$ 或 $Q_2$ 元素的问题，以下考虑可能是有用的。首先，解决温度方程在整个方案中几乎不是一个因素，因为几乎所有的计算时间都用于解决每个时间步骤中的斯托克斯系统。因此，温度方程的高阶元素并不是一个重要的缺点。另一方面，如果比较一下由于不连续的源描述而产生的过冲和欠冲的大小，我们会注意到，对于上述 $\beta$ 和 $k$ 的选择， $Q_1$ 的解决方案下降到 $-0.47$ 左右，而 $Q_2$ 的解决方案只到 $-0.13$ （记住，精确解决方案根本不应该变成负数。这意味着 $Q_2$ 解明显更准确；因此程序使用这些高阶元素，尽管我们在较小的时间步长方面付出了代价。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


有各种方法来扩展当前的程序。当然，特别感兴趣的是使其更快和/或提高程序的分辨率，特别是在3D方面。这就是step-32教程程序的主题，它将实现在集群上以%并行方式解决这个问题的策略。它也是更大的开放源代码ASPECT（见https://aspect.geodynamics.org/）的基础，它可以解决现实问题，并构成step-32的进一步发展。

另一个方向是使流体流动更加真实。这个程序最初是为了模拟各种情况，模拟地幔中的物质对流，即外地核和固体地壳之间的区域：在那里，物质从下面被加热，从上面被冷却，导致热对流。然而，这种流体的物理学要比这个程序中显示的复杂得多。地幔材料的粘度与温度有很大的关系，即 $\eta=\eta(T)$ ，这种关系经常被模拟为粘度随温度升高而呈指数下降。其次，地幔的大部分动态是由化学反应决定的，主要是构成地幔的各种晶体的相变；然后，斯托克斯方程右边的浮力项不仅取决于温度，而且还取决于某个特定位置的化学成分，这些化学成分被流场平流，但也作为压力和温度的函数而变化。我们将在以后的教程程序中也研究其中的一些影响。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-31.cc"
*/
