/**
@page step_67 The step-67 tutorial program
This tutorial depends on step-33, step-48, step-59.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#TheEulerequations">The Euler equations</a>
        <li><a href="#HighorderdiscontinuousGalerkindiscretization">High-order discontinuous Galerkin discretization</a>
        <li><a href="#Explicittimeintegration">Explicit time integration</a>
        <li><a href="#Fastevaluationofintegralsbymatrixfreetechniques">Fast evaluation of integrals by matrix-free techniques</a>
        <li><a href="#Evaluationoftheinversemassmatrixwithmatrixfreetechniques">Evaluation of the inverse mass matrix with matrix-free techniques</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#LowstorageexplicitRungeKuttatimeintegrators">Low-storage explicit Runge&mdash;Kutta time integrators</a>
        <li><a href="#ImplementationofpointwiseoperationsoftheEulerequations">Implementation of point-wise operations of the Euler equations</a>
        <li><a href="#TheEulerOperationclass">The EulerOperation class</a>
      <ul>
        <li><a href="#Localevaluators">Local evaluators</a>
        <li><a href="#Theapplyandrelatedfunctions">The apply() and related functions</a>
      </ul>
        <li><a href="#TheEulerProblemclass">The EulerProblem class</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
        <li><a href="#Convergenceratesfortheanalyticaltestcase">Convergence rates for the analytical test case</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin2D">Results for flow in channel around cylinder in 2D</a>
        <li><a href="#Resultsforflowinchannelaroundcylinderin3D">Results for flow in channel around cylinder in 3D</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Moreadvancednumericalfluxfunctionsandskewsymmetricformulations">More advanced numerical flux functions and skew-symmetric formulations</a>
        <li><a href="#Equippingthecodeforsupersoniccalculations">Equipping the code for supersonic calculations</a>
        <li><a href="#ExtensiontothelinearizedEulerequations">Extension to the linearized Euler equations</a>
        <li><a href="#ExtensiontothecompressibleNavierStokesequations">Extension to the compressible Navier-Stokes equations</a>
        <li><a href="#Usingcellcentricloopsandsharedmemory">Using cell-centric loops and shared memory</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-67/doc/intro.dox



 <br> 

<i>
This program was contributed by Martin Kronbichler. Many ideas presented here
are the result of common code development with Niklas Fehn, Katharina Kormann,
Peter Munch, and Svenja Schoeder.


This work was partly supported by the German Research Foundation (DFG) through
the project "High-order discontinuous Galerkin for the exa-scale" (ExaDG)
within the priority program "Software for Exascale Computing" (SPPEXA).
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序使用显式时间积分器求解流体力学的欧拉方程，其无矩阵框架应用于空间的高阶非连续Galerkin离散化。关于欧拉系统的细节和另一种隐式方法，我们也参考了第33步教程程序。你可能还想看看第69步，看看解决这些方程的另一种方法。




<a name="TheEulerequations"></a><h3>The Euler equations</h3>


欧拉方程是一个守恒定律，描述了一个可压缩的无粘性气体的运动。

@f[
\frac{\partial \mathbf{w}}{\partial t} + \nabla \cdot \mathbf{F}(\mathbf{w}) =
\mathbf{G}(\mathbf w),


@f]

其中解向量的 $d+2$ 分量为 $\mathbf{w}=(\rho, \rho
u_1,\ldots,\rho u_d,E)^{\mathrm T}$  。这里， $\rho$  表示流体密度， ${\mathbf u}=(u_1,\ldots, u_d)^\mathrm T$  表示流体速度， $E$  表示气体的能量密度。速度不直接求解，而是用变量 $\rho \mathbf{u}$ ，即线性动量（因为这是一个守恒量）。

欧拉通量函数是一个 $(d+2)\times d$ 矩阵，定义为

@f[
  \mathbf F(\mathbf w)
  =
  \begin{pmatrix}
  \rho \mathbf{u}\\
  \rho \mathbf{u} \otimes \mathbf{u} + \mathbb{I}p\\
  (E+p)\mathbf{u}
  \end{pmatrix}


@f]

其中 $\mathbb{I}$ 为 $d\times d$ 身份矩阵， $\otimes$ 为外积；其组成部分分别表示质量、动量和能量通量。右手边的强制力由以下公式给出

@f[
  \mathbf G(\mathbf w)
  =
  \begin{pmatrix}
  0\\
  \rho\mathbf{g}\\
  \rho \mathbf{u} \cdot \mathbf{g}
  \end{pmatrix},


@f]

其中矢量 $\mathbf g$ 表示重力的方向和大小。然而，它也可以表示作用于流体的任何其他单位质量的外力。例如，想想外部电场对带电粒子所施加的静电力）。

这三块方程，第二块涉及 $d$ 成分，描述了质量、动量和能量的守恒。压力不是一个解决方案的变量，但需要通过其他变量的 "闭合关系 "来表达；我们在此选择适合由两个原子组成的分子的气体的关系，在中等温度下，由 $p=(\gamma - 1) \left(E-\frac 12 \rho
\mathbf{u}\cdot \mathbf{u}\right)$ 和常数 $\gamma = 1.4$ 给出。




<a name="HighorderdiscontinuousGalerkindiscretization"></a><h3>High-order discontinuous Galerkin discretization</h3>


对于空间离散化，我们使用高阶非连续加勒金（DG）离散化，使用的解扩展形式为

@f[
\mathbf{w}_h(\mathbf{x}, t) =
\sum_{j=1}^{n_\mathbf{dofs}} \boldsymbol{\varphi}_j(\mathbf{x}) {w}_j(t).


@f]

这里， $\boldsymbol{\varphi}_j$ 表示第 $j$ 个基函数，以矢量形式写出不同成分的独立形状函数，让 $w_j(t)$ 分别通过密度、动量和能量变量。在这种形式下，空间依赖性包含在形状函数中，时间依赖性包含在未知系数中  $w_j$  。与连续有限元方法中一些形状函数跨越元素边界不同，在DG方法中，形状函数是单个元素的局部，从一个元素到下一个元素是不连续的。从一个单元到其相邻单元的解的连接是由下面规定的数值通量来实现的。这允许一些额外的灵活性，例如，在数值方法中引入方向性，例如，上卷。

DG方法是解决传输特性问题的流行方法，因为它们结合了低分散误差和勉强解决的尺度上的可控耗散。这使得它们在流体动力学领域的模拟中特别有吸引力，因为在这个领域中，需要代表广泛的活动尺度，不充分解决的特征很容易干扰重要的良好解决的特征。此外，高阶DG方法非常适用于现代硬件的正确实施。同时，DG方法也不是万能的。特别是当解出现不连续（冲击）时，就像欧拉方程在某些流态下的典型情况一样，高阶DG方法容易出现振荡解，就像所有不使用通量或坡度限制器的高阶方法一样。这是<a
href="https://en.wikipedia.org/wiki/Godunov%27s_theorem">Godunov's theorem</a>的结果，即任何线性的总变差（TVD）方案（如基本的DG离散化）最多只能达到一阶精度。换句话说，由于DG方法的目标是高阶精度，因此它们不可能对出现冲击的解进行TVD。尽管有些人声称DG方法中的数值通量可以控制耗散，但除非问题中的<b>all</b>冲击与单元边界对齐，否则这一点的价值有限。任何穿过单元内部的冲击都会因为高阶多项式而再次产生振荡分量。在有限元和DG界，存在许多不同的方法来处理冲击，例如在有问题的单元上引入人工扩散（使用基于解的模态分解等的有问题单元指标），在子网格上转换为耗散性低阶有限体积方法，或者增加一些限制程序。考虑到这种情况下的大量可能性，再加上相当大的实施努力，我们在这里不考虑带有明显冲击的欧拉方程系统，而是集中在带有波浪状现象的亚音速流动系统。对于一个能很好地处理冲击的方法（但每个未知数的成本较高），我们可以参考step-69教程程序。

对于DG公式的推导，我们将欧拉方程与测试函数 $\mathbf{v}$ 相乘，并对单个单元进行积分 $K$ ，从而得到

@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}
+ \left(\mathbf{v}, \nabla \cdot \mathbf{F}(\mathbf{w})\right)_{K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.


@f]



然后我们对第二项进行分项积分，将分歧从解槽移到测试函数槽，并产生一个元素边界上的积分。

@f[
\left(\mathbf{v}, \frac{\partial \mathbf{w}}{\partial t}\right)_{K}


- \left(\nabla \mathbf{v}, \mathbf{F}(\mathbf{w})\right)_{K}
+ \left<\mathbf{v}, \mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})
\right>_{\partial K} =
\left(\mathbf{v},\mathbf{G}(\mathbf w)\right)_{K}.


@f]

在表面积分中，我们用术语 $\widehat{\mathbf{F}}(\mathbf w)$ 代替了术语 $\mathbf{F}(\mathbf w)$ ，即数值通量。数字通量的作用是连接相邻元素上的解，并弱化解的连续性。这保证了PDE的全局耦合反映在离散化中，尽管单元上有独立的基函数。通过将数值通量定义为来自内部面两侧的解的函数 $\widehat{\mathbf{F}}(\mathbf w^-,
\mathbf w^+)$ 和 $\mathbf w^+$ ，包括与邻居的连接。我们要求的一个基本属性是，数值通量需要是<b>conservative</b>。也就是说，我们希望所有的信息（即质量、动量和能量）在一个面上离开一个单元时，都能完整地进入邻近的单元，反之亦然。这可以表示为 $\widehat{\mathbf{F}}(\mathbf w^-, \mathbf w^+) =
\widehat{\mathbf{F}}(\mathbf w^+, \mathbf w^-)$ ，也就是说，数值通量从任何一边都评估为相同的结果。结合数值通量与所考虑的面的单位外法向量相乘的事实，即从两边指向相反的方向，我们看到守恒被满足了。数值通量的另一个观点是作为一个单值的中间状态，从两边微弱地连接解决方案。

有大量的数值通量函数可用，也称为黎曼解算器。对于欧拉方程，存在所谓的精确黎曼求解器--意味着来自双方的状态以一种与欧拉方程沿线不连续的方式结合起来--以及近似黎曼求解器，它违反了一些物理特性，并依靠其他机制来使方案总体上准确。近似黎曼求解器的优点是计算起来比较便宜。大多数通量函数都起源于有限体积界，它们类似于单元（称为体积）内的多项式0度的DG方法。由于欧拉算子 $\mathbf{F}$ 的体积积分对于恒定解和检验函数会消失，所以数值通量必须完全代表物理算子，这也解释了为什么该界有大量的研究。对于DG方法，一致性是由单元内的高阶多项式保证的，这使得数值通量不再是一个问题，通常只影响收敛率，例如，对于度数为 $\mathcal O(h^p)$ 的多项式，解是否收敛为 $\mathcal O(h^{p+1/2})$ 或 $\mathcal
O(h^{p+1})$ 的准则。因此，数值通量可以被看作是一种机制，用于选择更有利的耗散/分散特性或关于离散化和线性化算子的极值特征，这影响到显式时间积分器中最大的可接受的时间步长。

在这个教程程序中，我们实现了两种通量的变体，可以通过程序中的开关来控制（当然，要使它们成为通过输入文件控制的运行时参数也很容易）。第一个通量是本地的Lax--Friedrichs通量

@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
   \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
   \mathbf{n^-}.


@f]



在Lax--Friedrichs通量的原始定义中，使用了一个系数 $\lambda =
\max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ （对应于信息在界面两边移动的最大速度），说明两个状态之间的差异， $[\![\mathbf{w}]\!]$ 被欧拉通量中的最大特征值惩罚，即 $\|\mathbf{u}\|+c$  ，其中 $c=\sqrt{\gamma p / \rho}$  是音速。在下面的实现中，我们对惩罚项进行了一些修改，因为无论如何惩罚都是近似的。我们使用

@f{align*}{
\lambda
&=
\frac{1}{2}\max\left(\sqrt{\|\mathbf{u^-}\|^2+(c^-)^2},
                     \sqrt{\|\mathbf{u}^+\|^2+(c^+)^2}\right)
\\
&=
\frac{1}{2}\sqrt{\max\left(\|\mathbf{u^-}\|^2+(c^-)^2,
                           \|\mathbf{u}^+\|^2+(c^+)^2\right)}.


@f}

额外的因子 $\frac 12$ 降低了惩罚强度（这导致特征值的负实部减少，从而增加了可接受的时间步长）。使用和内的平方允许我们减少昂贵的平方根操作的数量，对于原始的Lax--Friedrichs定义是4个，现在只需要一个。这种简化导致参数 $\lambda$ 的减少最多为2倍，因为 $\|\mathbf{u}\|^2+c^2 \leq
\|\mathbf{u}\|^2+2 c |\mathbf{u}\| + c^2 = \left(\|\mathbf{u}\|+c\right)^2
\leq 2 \left(\|\mathbf{u}\|^2+c^2\right)$ ，最后一个不等式来自杨氏不等式。

第二个数值通量是由Harten、Lax和van Leer提出的，称为HLL通量。它考虑到欧拉方程的不同传播方向，取决于声速。它利用一些中间状态  $\bar{\mathbf{u}}$  和  $\bar{c}$  来定义两个分支  $s^\mathrm{p} = \max\left(0, \bar{\mathbf{u}}\cdot \mathbf{n} +
\bar{c}\right)$  和  $s^\mathrm{n} = \min\left(0, \bar{\mathbf{u}}\cdot
\mathbf{n} - \bar{c}\right)$  。从这些分支中，人们再定义出通量

@f[
\hat{\mathbf{F}}(\mathbf{w}^-,\mathbf{w}^+) =
\frac{s^\mathrm{p} \mathbf{F}(\mathbf{w}^-)-s^\mathrm{n} \mathbf{F}(\mathbf{w}^+)}
                   {s^\mathrm p - s^\mathrm{n} } +
\frac{s^\mathrm{p} s^\mathrm{n}}{s^\mathrm{p}-s^\mathrm{n}}
\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes \mathbf{n^-}.


@f]

关于中间状态的定义  $\bar{\mathbf{u}}$  和  $\bar{c}$  ，已经提出了几个变种。最初提出的变体使用密度平均的速度定义，  $\bar{\mathbf{u}}
= \frac{\sqrt{\rho^-} \mathbf{u}^- + \sqrt{\rho^+}\mathbf{u}^+}{\sqrt{\rho^-}
+ \sqrt{\rho^+}}$  。由于我们考虑的是没有冲击的欧拉方程，因此在本教程程序中，我们简单地使用算术平均值， $\bar{\mathbf{u}} = \frac{\mathbf{u}^- +
\mathbf{u}^+}{2}$ 和 $\bar{c} = \frac{c^- + c^+}{2}$ ，与 $c^{\pm} =
\sqrt{\gamma p^{\pm} / \rho^{\pm}}$ ，而把其他变体留给可能的扩展。我们还注意到，HLL通量在文献中被扩展为所谓的HLLC通量，其中C代表表示接触不连续的能力。

在没有邻接状态 $\mathbf{w}^+$ 的边界上，通常的做法是从边界条件中推导出合适的外部值（详见关于DG方法的一般文献）。在这个教程程序中，我们考虑三种类型的边界条件，即<b>inflow boundary conditions</b>，其中所有分量都是规定的。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho_\mathrm{D}(t)\\
(\rho \mathbf u)_{\mathrm D}(t) \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(Dirichlet)},


@f]

<b>subsonic outflow boundaries</b>，在这里我们不规定外部解，因为流场要离开域，而使用内部值；我们仍然需要规定能量，因为欧拉通量中还有一个传入特性。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E_\mathrm{D}(t)\end{pmatrix} \quad
 \text{(mixed Neumann/Dirichlet)},


@f]

和<b>wall boundary condition</b>，它们描述了无渗透配置。

@f[
\mathbf{w}^+ = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- - 2 [(\rho \mathbf u)^-\cdot \mathbf n] \mathbf{n}
 \\ E^-\end{pmatrix}.


@f]



解的多项式展开最后被插入到弱形式，测试函数被基函数取代。这就得到了一个空间上离散、时间上连续的非线性系统，其未知系数的数量有限  $w_j$  ,  $j=1,\ldots,n_\text{dofs}$  。关于DG方法中多项式度数的选择，截至2019年，文献中并没有关于什么多项式度数最有效的共识，决定取决于问题。高阶多项式可以确保更好的收敛率，因此对于中等到高精确度要求的<b>smooth</b>解来说，高阶多项式更有优势。同时，自由度所在的体积与表面的比率，随着高阶的增加而增加，这使得数值通量的影响变弱，通常会减少耗散。然而，在大多数情况下，解决方案是不平滑的，至少与可以承受的分辨率相比是不平滑的。例如，在不可压缩流体力学、可压缩流体力学以及与之相关的波浪传播课题中都是如此。在这个前渐进制度中，误差大约与数值分辨率成正比，而其他因素，如分散误差或耗散行为变得更加重要。非常高阶的方法往往被排除在外，因为它们带有根据未知数衡量的更多限制性的CFL条件，而且当涉及到表示复杂几何形状时，它们也不那么灵活。因此，2到6的多项式度数在实践中是最受欢迎的，例如见 @cite FehnWallKronbichler2019 中的效率评估和其中引用的参考文献。

<a name="Explicittimeintegration"></a><h3>Explicit time integration</h3>


为了进行时间离散化，我们稍微重新排列了弱的形式，并在所有单元上求和。

@f[
\sum_{K \in \mathcal T_h} \left(\boldsymbol{\varphi}_i,
\frac{\partial \mathbf{w}}{\partial t}\right)_{K}
=
\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w})\right)_{K}


-\left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w})\right>_{\partial K} +
\left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w)\right)_{K}
\right],


@f]

其中 $\boldsymbol{\varphi}_i$ 贯穿了从1到 $n_\text{dofs}$ 的所有基函数。

我们现在用 $\mathcal M$ 表示质量矩阵，其条目为 $\mathcal M_{ij} =
\sum_{K} \left(\boldsymbol{\varphi}_i,
\boldsymbol{\varphi}_j\right)_K$ ，并用

@f[
\mathcal L_h(t,\mathbf{w}_h) = \left[\sum_{K\in \mathcal T_h}
\left[
\left(\nabla \boldsymbol{\varphi}_i, \mathbf{F}(\mathbf{w}_h)\right)_{K}


- \left<\boldsymbol{\varphi}_i,
\mathbf{n} \cdot \widehat{\mathbf{F}}(\mathbf{w}_h)\right>_{\partial K}
+ \left(\boldsymbol{\varphi}_i,\mathbf{G}(\mathbf w_h)\right)_{K}
\right]\right]_{i=1,\ldots,n_\text{dofs}}.


@f]

给定一个与全局未知数矢量和使用中的有限元相关的函数 $\mathbf{w}_h$ ，对欧拉算子的右手边进行评估的算子。这个函数 $\mathcal L_h$ 是明确随时间变化的，因为在边界上评估的数值通量将涉及边界某些部分的随时间变化的数据 $\rho_\mathrm{D}$ 、 $(\rho \mathbf{u})_\mathrm{D}$ 和 $E_\mathbf{D}$ ，取决于边界条件的分配。有了这个符号，我们可以把空间上的离散、时间上的连续系统紧凑地写为

@f[
\mathcal M \frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal L_h(t, \mathbf{w}_h),


@f]

其中我们冒昧地用 $\mathbf{w}_h$ 表示全局解矢量（除了相应的有限元函数外）。等价地，上述系统的形式为

@f[
\frac{\partial \mathbf{w}_h}{\partial t} =
\mathcal M^{-1} \mathcal L_h(t, \mathbf{w}_h).


@f]



对于用高阶非连续Galerkin方法离散的双曲系统，该系统的显式时间积分非常流行。这是由于质量矩阵 $\mathcal M$ 是块对角线的（每个块只对应于定义在同一单元上的同类变量），因此很容易倒置。在每个时间步长--或Runge-Kutta方案的阶段--我们只需要用给定的数据评估一次微分算子，然后应用质量矩阵的逆。另一方面，对于隐式时间步进，人们首先必须将方程线性化，然后迭代解决线性系统，这涉及到几个残差评估和至少十几个线性化算子的应用，正如在步骤33教程程序中所展示的那样。

当然，显式时间步长的简单性是有代价的，即由于所谓的Courant-Friedrichs-Lewy（CFL）条件而产生的条件稳定性。它指出，时间步长不能大于离散微分算子的最快信息传播速度。用更现代的术语来说，传播速度对应于离散算子的最大特征值，反过来又取决于网格大小、多项式程度 $p$ 和欧拉算子的物理学，即 $\mathbf F(\mathbf w)$ 相对于 $\mathbf{w}$ 的线性化的特征值。在这个程序中，我们设定的时间步长如下。

@f[
\Delta t = \frac{\mathrm{Cr}}{p^{1.5}}\left(\frac{1}
           {\max\left[\frac{\|\mathbf{u}\|}{h_u} + \frac{c}{h_c}\right]}\right),


@f]



在所有正交点和所有单元中取最大值。无量纲数 $\mathrm{Cr}$ 表示库朗数，可以选择最大稳定数 $\mathrm{Cr}_\text{max}$ ，其值取决于所选择的时间步进方法及其稳定性。用于多项式缩放的幂 $p^{1.5}$ 是启发式的，代表1到8之间的多项式度数最接近，例如，见 @cite SchoederKormann2018 。在更高的度数限制下， $p>10$ ， $p^2$ 的比例更准确，与通常用于内部惩罚方法的逆向估计有关。关于公式中使用的<i>effective</i>网格尺寸 $h_u$ 和 $h_c$ ，我们注意到对流传输是定向的。因此，一个合适的比例是使用速度方向的元素长度  $\mathbf u$  。下面的代码从参考单元到实际单元的雅各布系数的倒数得出这个比例，也就是说，我们近似于  $\frac{\|\mathbf{u}\|}{h_u} \approx \|J^{-1} \mathbf
u\|_{\infty}$  。相反，声波具有各向同性的特点，这就是为什么我们使用最小的特征尺寸，由 $J$ 的最小奇异值代表，用于声学缩放  $h_c$  。最后，我们需要增加对流和声学限制，因为欧拉方程可以以速度传输信息  $\|\mathbf{u}\|+c$  。

在这个教程程序中，我们使用<a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods">explicit
Runge--Kutta methods</a>的一个特定变体，一般来说，它使用以下更新程序，从时间 $t^n$ 的状态 $\mathbf{w}_h^{n}$ 到新时间 $t^{n+1}$ 的 $\Delta t = t^{n+1}-t^n$  。

@f[
\begin{aligned}
\mathbf{k}_1 &= \mathcal M^{-1} \mathcal L_h\left(t^n, \mathbf{w}_h^n\right),
\\
\mathbf{k}_2 &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_2\Delta t,
                       \mathbf{w}_h^n + a_{21} \Delta t \mathbf{k}_1\right),
\\
&\vdots \\
\mathbf{k}_s &= \mathcal M^{-1} \mathcal L_h\left(t^n+c_s\Delta t,
  \mathbf{w}_h^n + \sum_{j=1}^{s-1} a_{sj} \Delta t \mathbf{k}_j\right),
\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^n + \Delta t\left(b_1 \mathbf{k}_1 +
b_2 \mathbf{k}_2 + \ldots + b_s \mathbf{k}_s\right).
\end{aligned}


@f]

在 $\mathbf{k}_i$ 、 $i=1,\ldots,s$ 的阶段性方案中，向量 $s$ 是算子在某个中间状态下的评价，并通过某种线性组合用于定义阶段性结束值 $\mathbf{w}_h^{n+1}$ 。该方案中的标量系数 $c_i$ 、 $a_{ij}$ 和 $b_j$ 的定义，使得高阶方案满足某些条件，最基本的是 $c_i = \sum_{j=1}^{i-1}a_{ij}$  。参数通常以所谓的<a
href="https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Explicit_Runge%E2%80%93Kutta_methods">Butcher
tableau</a>的形式收集，它收集了定义该方案的所有系数。对于一个五级方案，它看起来是这样的。

@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_{21} \\
c_3 & a_{31} & a_{32} \\
c_4 & a_{41} & a_{42} & a_{43} \\
c_5 & a_{51} & a_{52} & a_{53} & a_{54} \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}


@f]



在这个教程程序中，我们使用显式Runge--Kutta方法的一个子集，即所谓的低存储Runge--Kutta方法（LSRK），它假定了系数的额外结构。在参考文献 @cite KennedyCarpenterLewis2000 所使用的变体中，假设使用的是Butcher tableaus的形式

@f[
\begin{array}{c|ccccc}
0 \\
c_2 & a_1 \\
c_3 & b_1 & a_2 \\
c_4 & b_1 & b_2 & a_3 \\
c_5 & b_1 & b_2 & b_3 & a_4 \\
\hline
& b_1 & b_2 & b_3 & b_4 & b_5
\end{array}


@f]

有了这样的定义，对  $\mathbf{w}_h^n$  的更新与中间值  $\mathbf{k}_i$  的信息共享存储。从 $\mathbf{w}^{n+1}=\mathbf{w}^n$ 和 $\mathbf{r}_1 = \mathbf{w}^n$ 开始，每个 $s$ 阶段的更新都简化为

@f[
\begin{aligned}
\mathbf{k}_i &=
\mathcal M^{-1} \mathcal L_h\left(t^n+c_i\Delta t, \mathbf{r}_{i} \right),\\
\mathbf{r}_{i+1} &= \mathbf{w}_h^{n+1} + \Delta t \, a_i \mathbf{k}_i,\\
\mathbf{w}_h^{n+1} &= \mathbf{w}_h^{n+1} + \Delta t \, b_i \mathbf{k}_i.
\end{aligned}


@f]

除了连续更新的向量 $\mathbf w_h^{n+1}$ ，这个方案只需要两个辅助向量，即保存微分算子的评估的向量 $\mathbf{k}_i$ ，以及保存微分算子应用的右手边的向量 $\mathbf{r}_i$ 。在后续阶段  $i$  ，值  $\mathbf{k}_i$  和  $\mathbf{r}_i$  可以使用相同的存储。

低存储量变体的主要优点是一方面减少了内存消耗（如果必须在内存中装入非常多的未知数，持有所有的 $\mathbf{k}_i$ 来计算随后的更新，对于 $s$ 来说已经是一个极限，在5到8之间--记得我们使用的是显式方案，所以我们不需要存储任何通常比几个向量大很多的矩阵），另一方面是减少内存访问。在这个程序中，我们对后一个方面特别感兴趣。由于运算符评估的成本只是简单地从内存中流转输入和输出向量的一小部分，我们必须考虑向量更新的成本，而低存储的变体可以提供传统显式Runge--Kutta方法两倍的吞吐量，原因就在于此，例如，见 @cite SchoederKormann2018 中的分析。

除了参考文献 @cite KennedyCarpenterLewis2000 中的三阶、四阶和五阶精度的三个变体外，我们还使用了一个四阶精度的七级变体，该变体是为声学设置而优化的 @cite TseliosSimos2007  。声学问题是欧拉方程的亚音速制度的有趣方面之一，其中可压缩性导致了声波的传播；通常，人们使用围绕背景状态的线性化欧拉方程的进一步简化，或围绕固定框架的声波方程。




<a name="Fastevaluationofintegralsbymatrixfreetechniques"></a><h3>Fast evaluation of integrals by matrix-free techniques</h3>


这个程序中使用的主要成分是我们用来评估算子  $\mathcal L_h$  和反质量矩阵  $\mathcal M$  的快速无矩阵技术。实际上，<i>matrix-free</i>这个术语有点名不副实，因为我们是在处理一个非线性算子，并没有将反过来可以用矩阵表示的算子线性化。然而，作为稀疏矩阵-向量乘积的替代品，积分的快速评估已经变得很流行，如步骤-37和步骤-59所示，为此我们在交易二中创造了这个基础设施<i>matrix-free functionality</i>。此外，反质量矩阵确实是以无矩阵的方式应用的，详见下文。

无矩阵基础设施使我们能够快速评估弱形式的积分。其成分是将解系数快速插值为正交点的值和导数，在正交点进行逐点运算（在这里我们实现了上述的微分算子），以及与所有测试函数相乘和对正交点求和。第一和第三部分利用了和因子化，并在步骤37的单元积分教程和步骤59的面积分教程中进行了广泛的讨论。唯一的区别是，我们现在处理的是一个 $d+2$ 分量的系统，而不是以前教程程序中的标量系统。在代码中，所有的变化是FEEvaluation和FEFaceEvaluation类的一个模板参数，即设置分量的数量。对向量的访问和以前一样，都由评价器透明地处理。我们还注意到，下面的代码中选择的带有单一评价器的变体并不是唯一的选择--我们也可以为单独的组件 $\rho$ 、 $\rho \mathbf u$ 和 $E$ 使用单独的评价器；鉴于我们对所有组件的处理是类似的（也反映在我们把方程作为一个矢量系统的方式），这里会更复杂。和以前一样，FEEvaluation类通过结合对几个单元（和面）的操作来提供显式的矢量化，涉及的数据类型称为VectorizedArray。由于这种类型的算术运算都是重载的，所以我们不必为它费心，除了通过函数接口对函数进行评估，我们需要同时为几个正交点的位置提供特殊的<i>vectorized</i>评估。

这个程序中更大的变化是在正交点的操作。在这里，多分量评估器为我们提供了之前没有讨论过的返回类型。 FEEvaluation::get_value() 将为第37步的拉普拉斯返回一个标量（更准确地说，由于跨单元的矢量化，是一个VectorizedArray类型），现在它返回的类型是`Tensor<1,dim+2,VectorizedArray<Number>'。同样，梯度类型现在是`张量<1,dim+2,张量<1,dim,矢量化数组<Number>>`，其中外部张量收集了欧拉系统的`dim+2'分量，内部张量是各个方向的偏导数。例如，欧拉系统的通量 $\mathbf{F}(\mathbf{w})$ 就属于这种类型。为了减少我们为拼出这些类型而写的代码量，我们尽可能使用C++的`自动'关键字。

从实施的角度来看，非线性并不是一个很大的困难。它是在我们表达欧拉弱形式的条款时自然引入的，例如以动量条款的形式  $\rho \mathbf{u}
\otimes \mathbf{u}$  。为了得到这个表达式，我们首先从动量变量  $\rho \mathbf{u}$  推导出速度  $\mathbf{u}$  。鉴于 $\rho
\mathbf{u}$ 和 $\rho$ 一样被表示为 $p$ 度的多项式，速度 $\mathbf{u}$ 是参考坐标 $\hat{\mathbf{x}}$ 的一个有理表达。当我们进行乘法 $(\rho
\mathbf{u})\otimes \mathbf{u}$ 时，我们得到一个表达式，它是两个多项式的比值，分子中的多项式程度 $2p$ 和分母中的多项式程度 $p$ 。结合测试函数的梯度，分子中的积分度为 $3p$ ，分母中的积分度为 $p$ ，对于仿生单元，即平行四边形/平行四边形，已经有了积分。对于弧形单元，当积分乘以映射的雅各布系数时，会出现额外的多项式和有理表达式。在这一点上，人们通常需要放弃坚持精确的积分，而采取高斯（更确切地说，高斯-勒格伦德）正交提供的任何精度。这时的情况与拉普拉斯方程的情况类似，积分项包含非affince单元上的有理表达式，也只能进行近似积分。由于这些公式只对多项式进行精确积分，我们不得不以积分错误的形式忍受<a
href="https://mathoverflow.net/questions/26018/what-are-variational-crimes-and-who-coined-the-term">variational
crime</a>的影响。

虽然对于椭圆问题来说，不精确的积分通常是可以容忍的，但对于双曲问题来说，不精确的积分会引起一些令人头痛的效应，这种效应称为<b>aliasing</b>。这个术语来自于信号处理，表达了不适当的、过于粗糙的采样情况。就正交而言，不适当的采样意味着我们使用的正交点与准确采样变系数积分所需的点相比太少。在DG文献中已经表明，别离误差会在<i>barely</i>解析模拟的数值解中引入非物理性的振荡。别名主要影响到粗略的分辨率--而采用相同方案的更细的网格则工作良好--这一事实并不令人惊讶，因为分辨率高的模拟往往在一个单元的长度尺度上是平滑的（即，它们在较高的多项式程度上有小的系数，由于正交点太少而被遗漏，而在较低的多项式程度上的主要解贡献仍然被很好地捕获--这只是泰勒定理的一个结果）。为了解决这个问题，DG文献中提出了各种方法。一种技术是过滤，它可以抑制与高次多项式度数有关的解成分。由于所选择的节点基不是分层的，这就意味着要从节点基转化为分层基（例如，基于Legendre多项式的模态基），其中单元内的贡献是按多项式程度划分的。在这个基础上，我们可以将与高度数相关的求解系数乘以一个小数，保持低度数不变（以避免破坏一致性），然后再转换回节点基础。然而，过滤器会降低该方法的准确性。另一个在某种意义上更简单的策略是使用更多的正交点来更准确地捕捉非线性项。每个坐标方向使用超过 $p+1$ 个正交点有时被称为过度积分或一致积分。后者在不可压缩的Navier-Stokes方程中最为常见，其中 $\mathbf{u}\otimes \mathbf{u}$ 非线性导致 $3p$ 度的多项式积分（当同时考虑测试函数时），只要元素的几何形状是仿生的，每个方向的 $\textrm{floor}\left(\frac{3p}{2}\right)+1$ 正交点就可以精确积分。在非多项式积分的欧拉方程的背景下，选择就不那么明确了。根据各种变量的变化， $\textrm{floor}\left(\frac{3p}{2}\right)+1$ 或 $2p+1$ 点（分别精确积分度为 $3p$ 或 $4p$ 的多项式）都很常见。

为了反映程序中正交选择的这种可变性，我们把正交点的数量作为一个变量来指定，就像多项式的度数一样，并注意到人们会根据流量配置做出不同的选择。默认选择是 $p+2$ 点--比最小可能的 $p+1$ 点多一点。FEEvaluation和FEFaceEvaluation类允许通过模板参数无缝地改变点的数量，这样程序就不会因此而变得更复杂。




<a name="Evaluationoftheinversemassmatrixwithmatrixfreetechniques"></a><h3>Evaluation of the inverse mass matrix with matrix-free techniques</h3>


最后一个要素是反质量矩阵的评估  $\mathcal
M^{-1}$  。在具有显式时间积分的DG方法中，质量矩阵是块状对角线，因此很容易反转--人们只需要反转对角线块。然而，考虑到无矩阵的积分评估在成本上更接近于只访问向量，即使应用块对角矩阵（例如通过LU因子数组）也会比评估 $\mathcal L_h$ 贵几倍，仅仅是因为对于高阶有限元来说，仅仅存储和加载大小为`dofs_per_cell`x`dofs_per_cell`的矩阵是昂贵的。由于这显然是不可取的，部分社区已经转移到质量矩阵是对角线的基础，例如<i>L<sub>2</sub></i>正交Legendre基础，使用分层多项式或高斯四分法点上的拉格朗日多项式（这只是利用Legendre信息的另一种方式）。虽然对角线属性对于变形元素来说是失效的，但通过采取对角线质量矩阵而忽略其余部分（质量包络的变种，尽管不是步骤-48中利用的具有额外积分误差的变种）所产生的误差已被证明不会改变离散化精度。高斯正交点中的拉格朗日基础有时也被称为同位设置，因为多项式的结点与正交点重合（="同位"），避免了一些内插操作。鉴于我们想在 $\mathcal L_h$ 中对非线性项使用更多的正交点，然而，拼合属性就失去了。(更确切地说，在改变基础后，它仍然用于FEEvaluation和FEFaceEvaluation，见无矩阵论文  @cite KronbichlerKormann2019  。)

在这个教程程序中，我们使用拼合思想来应用反质量矩阵，但有一个小的转折。与其在高斯四分法的点上通过拉格朗日多项式使用配位，我们更倾向于在高斯-洛巴托点上使用传统的拉格朗日基础，因为那些使面积分的评估变得便宜。这是因为对于高斯-洛巴托点来说，一些节点点位于单元格的面上，而且不难证明，在任何给定的面上，唯一具有非零值的形状函数正是其节点点实际上位于该面上的那些。当然，我们也可以像步骤48那样使用高斯-洛巴托正交（有一些额外的积分误差），但我们不想牺牲精度，因为这些正交公式通常比一般的高斯正交公式的阶数低。相反，我们使用参考文献 @cite KronbichlerSchoeder2016 中描述的一个想法，其中提出为了应用反质量矩阵而改变基础。让我们用 $S$ 表示在正交点评价的形状函数矩阵，形状函数在矩阵的行中，正交点在列中。那么，单元格 $K$ 上的质量矩阵由以下公式给出

@f[
\mathcal M^K = S J^K S^\mathrm T.


@f]

这里， $J^K$ 是以雅各布系数乘以正交权重（JxW）的行列式作为条目的对角矩阵。矩阵 $S$ 被构造为一维矩阵的克朗克积（张量积），例如，在三维中为

@f[
S = S_{\text{1D}}\otimes S_{\text{1D}}\otimes S_{\text{1D}},


@f]

这是基函数是一维形状函数的张量积，正交公式是一维正交公式的张量积的结果。对于多项式的数量等于正交点的数量的情况， $S J^K S^\mathrm T$ 中的所有矩阵都是方形的，同样，克朗克积中的 $S$ 的成分也是方形的。因此，人们可以对每个矩阵进行反转，形成整体的逆。

@f[
\left(\mathcal M^K\right)^{-1} = S_{\text{1D}}^{-\mathrm T}\otimes
S_{\text{1D}}^{-\mathrm T}\otimes S_{\text{1D}}^{-\mathrm T}
\left(J^K\right)^{-1}
S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}\otimes S_{\text{1D}}^{-1}.


@f]

这个公式的结构与用和因子化技术对积分进行正向评价的步骤完全相同（即交易.II的FEEvaluation和MatrixFree框架）。因此，我们可以利用相同的代码路径，采用不同的插值矩阵， $S_{\mathrm{1D}}^{-\mathrm{T}}$ 而不是 $S_{\mathrm{1D}}$  。

类 MatrixFreeOperators::CellwiseInverseMassMatrix 实现了这个操作。它从有限元中包含的基（在这里是FE_DGQ）改变为高斯正交点中的拉格朗日基。在这里，可以评估对角线质量矩阵的逆值，这只是`JxW`因子的逆值（即正交权重乘以从参考坐标到实坐标的雅各布系数）。一旦这样做了，我们就可以变回标准的节点高斯-洛巴托基础。

这种应用反质量矩阵的特殊方式的优点是成本类似于质量矩阵的正向应用，这比用超积分和面积分评估空间算子 $\mathcal L_h$ 更便宜。(我们将在<a href="#Results">results section</a>中用详细的时间信息证明这一点)。事实上，它是如此便宜，以至于在大多数现代架构上，它被读取源向量、读取对角线和写入目的向量的带宽所限制。用于结果部分的硬件可以使计算的速度至少比从内存流向量的速度快一倍。




<a name="Thetestcase"></a><h3>The test case</h3>


在这个教程程序中，我们实现了两个测试案例。第一个案例是限于两个空间维度的收敛性测试。它运行一个所谓的等熵涡旋，它通过一个背景流场进行传输。第二个案例使用了一个更令人兴奋的设置。我们从一个浸在通道中的圆柱体开始，使用 GridGenerator::channel_with_cylinder() 函数。在这里，我们强加一个马赫数为 $\mathrm{Ma}=0.307$ 的亚音速初始场，在 $x$ 方向上速度不变。在顶壁和底壁以及圆柱体上，我们施加了一个无穿透（即切向流动）的条件。与初始条件相比，这种设置迫使气流重新定向，从而导致大的声波从圆柱体上传播出去。在上游方向，波的传播速度较慢（因为它必须逆着迎面而来的气体移动），包括密度和压力的不连续。在下游方向，由于声音的传播和流体的流动方向相同，传输速度较快，这在一定程度上抹去了不连续性。一旦声波碰到上下壁，声音就会被反射回来，形成一些漂亮的形状，如下图<a href="#Results">results section</a>所示。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 包含文件与之前的无矩阵教程程序 step-37 、 step-48 和 step-59 相似。
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/time_stepping.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/vectorization.h> 
 * 
 * #include <deal.II/distributed/tria.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/fe/fe_system.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/tria.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/la_parallel_vector.h> 
 * 
 * #include <deal.II/matrix_free/fe_evaluation.h> 
 * #include <deal.II/matrix_free/matrix_free.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iomanip> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 下面的文件包括CellwiseInverseMassMatrix数据结构，我们将在质量矩阵反演中使用它，这是本教程程序中唯一的新包含文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/matrix_free/operators.h> 
 * 
 * namespace Euler_DG 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 与其他无矩阵教程程序类似，我们在文件的顶部收集所有控制程序执行的参数。除了我们想要运行的维度和多项式程度，我们还指定了我们想要用于欧拉方程中非线性项的高斯正交公式的点数。此外，我们指定了随时间变化的问题的时间间隔，并实现了两个不同的测试案例。第一个是二维的分析解，而第二个是介绍中描述的围绕圆柱体的通道流。根据测试案例，我们还改变了运行模拟的最终时间，以及一个变量`output_tick`，它指定了我们要在哪个时间间隔内写入输出（假设tick大于时间步长）。
 * 

 * 
 * 
 * @code
 *   constexpr unsigned int testcase             = 0; 
 *   constexpr unsigned int dimension            = 2; 
 *   constexpr unsigned int n_global_refinements = 3; 
 *   constexpr unsigned int fe_degree            = 5; 
 *   constexpr unsigned int n_q_points_1d        = fe_degree + 2; 
 * 
 *   using Number = double; 
 * 
 *   constexpr double gamma       = 1.4; 
 *   constexpr double final_time  = testcase == 0 ? 10 : 2.0; 
 *   constexpr double output_tick = testcase == 0 ? 1 : 0.05; 
 * 
 * @endcode
 * 
 * 接下来是时间积分器的一些细节，即用公式 $\Delta t =
 * \text{Cr} n_\text{stages} \frac{h}{(p+1)^{1.5} (\|\mathbf{u} +
 * c)_\text{max}}$ 来衡量时间步长的库朗数，以及选择一些低存储量的Runge--Kutta方法。我们指定Runge--Kutta方案每级的Courant数，因为这对不同级数的方案给出了一个更实际的数值成本表达。
 * 

 * 
 * 
 * @code
 *   const double courant_number = 0.15 / std::pow(fe_degree, 1.5); 
 *   enum LowStorageRungeKuttaScheme 
 *   { 
 *     stage_3_order_3, /* Kennedy, Carpenter, Lewis, 2000 */ 
 * 
 * 
 *     stage_5_order_4, /* Kennedy, Carpenter, Lewis, 2000 */ 
 * 
 * 
 *     stage_7_order_4, /* Tselios, Simos, 2007 */ 
 * 
 * 
 *     stage_9_order_5, /* Kennedy, Carpenter, Lewis, 2000 */ 
 * 
 * 
 *   }; 
 *   constexpr LowStorageRungeKuttaScheme lsrk_scheme = stage_5_order_4; 
 * 
 * @endcode
 * 
 * 最终，我们选择了空间离散化的一个细节，即单元间面的数值通量（黎曼求解器）。在这个程序中，我们实现了Lax--Friedrichs通量和Harten--Lax--van Leer(HLL)通量的一个改进版本。
 * 

 * 
 * 
 * @code
 *   enum EulerNumericalFlux 
 *   { 
 *     lax_friedrichs_modified, 
 *     harten_lax_vanleer, 
 *   }; 
 *   constexpr EulerNumericalFlux numerical_flux_type = lax_friedrichs_modified; 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 我们现在定义了一个带有测试情况0的精确解的类和一个带有测试情况1的通道背景流场的类。鉴于欧拉方程是一个在 $d$ 维度上有 $d+2$ 个方程的问题，我们需要告诉函数基类正确的分量数量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ExactSolution : public Function<dim> 
 *   { 
 *   public: 
 *     ExactSolution(const double time) 
 *       : Function<dim>(dim + 2, time) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 就实际实现的函数而言，分析性测试案例是一个等熵涡旋案例（例如参见Hesthaven和Warburton的书，第209页第6.6节中的例6.1），它满足欧拉方程，右侧的力项为零。考虑到这个定义，我们返回密度、动量或能量，这取决于所要求的成分。请注意，密度的原始定义涉及一些表达式的 $\frac{1}{\gamma -1}$ -次方。由于 `std::pow()` 在某些系统上的实现相当慢，我们用对数和指数（以2为底）来代替它，这在数学上是等价的，但通常优化得更好。与 `std::pow()`, 相比，对于非常小的数字，这个公式可能会在最后一位数字上失去准确性，但我们还是很高兴，因为小数字映射为接近1的数据。
 * 

 * 
 * 对于通道测试案例，我们简单地选择密度为1， $x$ 方向的速度为0.4，其他方向的速度为0，以及对应于背景速度场测量的1.3声速的能量，根据关系 $E = \frac{c^2}{\gamma (\gamma -1)} + \frac 12 \rho \|u\|^2$ 计算得出。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double ExactSolution<dim>::value(const Point<dim> & x, 
 *                                    const unsigned int component) const 
 *   { 
 *     const double t = this->get_time(); 
 * 
 *     switch (testcase) 
 *       { 
 *         case 0: 
 *           { 
 *             Assert(dim == 2, ExcNotImplemented()); 
 *             const double beta = 5; 
 * 
 *             Point<dim> x0; 
 *             x0[0] = 5.; 
 *             const double radius_sqr = 
 *               (x - x0).norm_square() - 2. * (x[0] - x0[0]) * t + t * t; 
 *             const double factor = 
 *               beta / (numbers::PI * 2) * std::exp(1. - radius_sqr); 
 *             const double density_log = std::log2( 
 *               std::abs(1. - (gamma - 1.) / gamma * 0.25 * factor * factor)); 
 *             const double density = std::exp2(density_log * (1. / (gamma - 1.))); 
 *             const double u       = 1. - factor * (x[1] - x0[1]); 
 *             const double v       = factor * (x[0] - t - x0[0]); 
 * 
 *             if (component == 0) 
 *               return density; 
 *             else if (component == 1) 
 *               return density * u; 
 *             else if (component == 2) 
 *               return density * v; 
 *             else 
 *               { 
 *                 const double pressure = 
 *                   std::exp2(density_log * (gamma / (gamma - 1.))); 
 *                 return pressure / (gamma - 1.) + 
 *                        0.5 * (density * u * u + density * v * v); 
 *               } 
 *           } 
 * 
 *         case 1: 
 *           { 
 *             if (component == 0) 
 *               return 1.; 
 *             else if (component == 1) 
 *               return 0.4; 
 *             else if (component == dim + 1) 
 *               return 3.097857142857143; 
 *             else 
 *               return 0.; 
 *           } 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *           return 0.; 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LowstorageexplicitRungeKuttatimeintegrators"></a> 
 * <h3>Low-storage explicit Runge--Kutta time integrators</h3>
 * 

 * 
 * 接下来的几行实现了一些低存储量的Runge--Kutta方法的变体。这些方法有特定的布彻表，系数为 $b_i$ 和 $a_i$ ，如介绍中所示。如同Runge--Kutta方法的惯例，我们可以从这些系数中推导出时间步骤 $c_i = \sum_{j=1}^{i-2} b_i + a_{i-1}$ 。这种方案的主要优点是每个阶段只需要两个向量，即解的累积部分 $\mathbf{w}$ （在最后一个阶段后的新时间 $t^{n+1}$ 保持解 $\mathbf{w}^{n+1}$ ），在各阶段被评估的更新向量 $\mathbf{r}_i$ ，加上一个向量 $\mathbf{k}_i$ 来保持算子评估。这样的Runge--Kutta设置减少了内存存储和内存访问。由于内存带宽通常是现代硬件上的性能限制因素，当微分算子的评估得到很好的优化时，性能可以比标准的时间积分器得到改善。考虑到传统的Runge--Kutta方案可能允许稍大的时间步长，因为更多的自由参数可以获得更好的稳定性，这一点也是真实的。
 * 

 * 
 * 在本教程中，我们集中讨论Kennedy, Carpenter和Lewis(2000)文章中定义的低存储方案的几个变体，以及Tselios和Simos(2007)描述的一个变体。还有一大系列的其他方案，可以通过额外的系数集或稍微不同的更新公式来解决。
 * 

 * 
 * 我们为这四种积分器定义了一个单一的类，用上述的枚举来区分。对每个方案，我们再将 $b_i$ 和 $a_i$ 的向量填充到类中的给定变量。
 * 

 * 
 * 
 * @code
 *   class LowStorageRungeKuttaIntegrator 
 *   { 
 *   public: 
 *     LowStorageRungeKuttaIntegrator(const LowStorageRungeKuttaScheme scheme) 
 *     { 
 *       TimeStepping::runge_kutta_method lsrk; 
 * 
 * @endcode
 * 
 * 首先是Kennedy等人（2000）提出的三阶方案。虽然它的稳定区域比其他方案小得多，但它只涉及三个阶段，所以在每个阶段的工作方面很有竞争力。
 * 

 * 
 * 
 * @code
 *       switch (scheme) 
 *         { 
 *           case stage_3_order_3: 
 *             { 
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE3_ORDER3; 
 *               break; 
 *             } 
 * 
 * @endcode
 * 
 * 下一个方案是四阶的五级方案，同样在Kennedy等人（2000）的论文中定义。
 * 

 * 
 * 
 * @code
 *           case stage_5_order_4: 
 *             { 
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE5_ORDER4; 
 *               break; 
 *             } 
 * 
 * @endcode
 * 
 * 下面这个七级和四阶的方案已经明确地推导出用于声学问题。它在四阶方案中兼顾了虚特征值的精度，并结合了一个大的稳定区域。由于DG方案在最高频率之间是耗散的，这不一定转化为每级可能的最高时间步长。在本教程方案的背景下，数值通量在耗散中起着至关重要的作用，因此也是最大的稳定时间步长。对于修改后的Lax--Friedrichs通量，如果只考虑稳定性，该方案在每级步长方面与`stage_5_order_4`方案相似，但对于HLL通量来说，效率稍低。
 * 

 * 
 * 
 * @code
 *           case stage_7_order_4: 
 *             { 
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE7_ORDER4; 
 *               break; 
 *             } 
 * 
 * @endcode
 * 
 * 这里包括的最后一个方案是Kennedy等人（2000）的五阶九级方案。它是这里使用的方案中最精确的，但是较高的精度牺牲了一些稳定性，所以每级的归一化步长比四阶方案要小。
 * 

 * 
 * 
 * @code
 *           case stage_9_order_5: 
 *             { 
 *               lsrk = TimeStepping::LOW_STORAGE_RK_STAGE9_ORDER5; 
 *               break; 
 *             } 
 * 
 *           default: 
 *             AssertThrow(false, ExcNotImplemented()); 
 *         } 
 *       TimeStepping::LowStorageRungeKutta< 
 *         LinearAlgebra::distributed::Vector<Number>> 
 *         rk_integrator(lsrk); 
 *       rk_integrator.get_coefficients(ai, bi, ci); 
 *     } 
 * 
 *     unsigned int n_stages() const 
 *     { 
 *       return bi.size(); 
 *     } 
 * 
 * @endcode
 * 
 * 时间积分器的主要功能是通过阶段，评估算子，为下一次评估准备  $\mathbf{r}_i$  矢量，并更新解决方案矢量  $\mathbf{w}$  。我们把工作交给所涉及的`pde_operator`，以便能够把Runge--Kutta设置的矢量操作与微分算子的评估合并起来，以获得更好的性能，所以我们在这里所做的就是委托矢量和系数。
 * 

 * 
 * 我们单独调用第一阶段的算子，因为我们需要稍微修改一下那里的参数。我们从旧的解决方案 $\mathbf{w}^n$ 而不是 $\mathbf r_i$ 向量中评估解决方案，所以第一个参数是`solution`。我们在这里让阶段向量 $\mathbf{r}_i$ 也持有评估的临时结果，因为它在其他情况下不会被使用。对于所有后续阶段，我们使用向量`vec_ki`作为第二个向量参数来存储运算符的求值结果。最后，当我们到了最后一个阶段，我们必须跳过对向量 $\mathbf{r}_{s+1}$ 的计算，因为没有系数 $a_s$ 可用（也不会用到）。
 * 

 * 
 * 
 * @code
 *     template <typename VectorType, typename Operator> 
 *     void perform_time_step(const Operator &pde_operator, 
 *                            const double    current_time, 
 *                            const double    time_step, 
 *                            VectorType &    solution, 
 *                            VectorType &    vec_ri, 
 *                            VectorType &    vec_ki) const 
 *     { 
 *       AssertDimension(ai.size() + 1, bi.size()); 
 * 
 *       pde_operator.perform_stage(current_time, 
 *                                  bi[0] * time_step, 
 *                                  ai[0] * time_step, 
 *                                  solution, 
 *                                  vec_ri, 
 *                                  solution, 
 *                                  vec_ri); 
 * 
 *       for (unsigned int stage = 1; stage < bi.size(); ++stage) 
 *         { 
 *           const double c_i = ci[stage]; 
 *           pde_operator.perform_stage(current_time + c_i * time_step, 
 *                                      bi[stage] * time_step, 
 *                                      (stage == bi.size() - 1 ? 
 *                                         0 : 
 *                                         ai[stage] * time_step), 
 *                                      vec_ri, 
 *                                      vec_ki, 
 *                                      solution, 
 *                                      vec_ri); 
 *         } 
 *     } 
 * 
 *   private: 
 *     std::vector<double> bi; 
 *     std::vector<double> ai; 
 *     std::vector<double> ci; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofpointwiseoperationsoftheEulerequations"></a> 
 * <h3>Implementation of point-wise operations of the Euler equations</h3>
 * 

 * 
 * 在下面的函数中，我们实现了与欧拉方程有关的各种特定问题的运算。每个函数都作用于我们在解向量中持有的守恒变量向量 $[\rho, \rho\mathbf{u}, E]$ ，并计算各种派生量。
 * 

 * 
 * 首先是速度的计算，我们从动量变量 $\rho \mathbf{u}$ 除以 $\rho$ 得出。这里需要注意的是，我们用关键字`DEAL_II_ALWAYS_INLINE`来装饰所有这些函数。这是一个特殊的宏，映射到一个编译器专用的关键字，告诉编译器永远不要为这些函数创建一个函数调用，而是将实现<a href="https:en.wikipedia.org/wiki/Inline_function">inline</a>移到它们被调用的地方。这对性能至关重要，因为我们对其中一些函数的调用达到了几百万甚至几十亿次。例如，我们既使用速度来计算通量，也使用速度来计算压力，而这两个地方都要在每个单元的每个正交点进行评估。确保这些函数是内联的，不仅可以确保处理器不必执行跳转指令进入函数（以及相应的返回跳转），而且编译器可以在调用函数的地方之后的代码中重新使用一个函数的上下文的中间信息。(我们注意到，编译器通常很善于自己找出哪些函数要内联。这里有一个地方，编译器可能是自己想出来的，也可能不是，但我们可以肯定的是，内联是一种胜利。)
 * 

 * 
 * 我们应用的另一个技巧是为反密度设置一个单独的变量  $\frac{1}{\rho}$  。这使得编译器只对通量进行一次除法，尽管除法在多个地方使用。由于除法的费用大约是乘法或加法的10到20倍，避免多余的除法对性能至关重要。我们注意到，由于四舍五入的影响，在浮点运算中，先取反数，后与之相乘并不等同于除法，所以编译器不允许用标准的优化标志来交换一种方式。然而，以正确的方式编写代码也不是特别困难。
 * 

 * 
 * 总而言之，所选择的总是内联和仔细定义昂贵的算术运算的策略使我们能够写出紧凑的代码，而不需要将所有的中间结果传递出去，尽管要确保代码映射到优秀的机器码。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim, Number> 
 *     euler_velocity(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Number inverse_density = Number(1.) / conserved_variables[0]; 
 * 
 *     Tensor<1, dim, Number> velocity; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       velocity[d] = conserved_variables[1 + d] * inverse_density; 
 * 
 *     return velocity; 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数从保守变量的矢量中计算压力，使用公式  $p = (\gamma - 1) \left(E - \frac 12 \rho \mathbf{u}\cdot \mathbf{u}\right)$  。如上所述，我们使用来自`euler_velocity()`函数的速度。注意，我们需要在这里指定第一个模板参数`dim`，因为编译器无法从张量的参数中推导出它，而第二个参数（数字类型）可以自动推导出来。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Number 
 *     euler_pressure(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Tensor<1, dim, Number> velocity = 
 *       euler_velocity<dim>(conserved_variables); 
 * 
 *     Number rho_u_dot_u = conserved_variables[1] * velocity[0]; 
 *     for (unsigned int d = 1; d < dim; ++d) 
 *       rho_u_dot_u += conserved_variables[1 + d] * velocity[d]; 
 * 
 *     return (gamma - 1.) * (conserved_variables[dim + 1] - 0.5 * rho_u_dot_u); 
 *   } 
 * 
 * @endcode
 * 
 * 这里是欧拉通量函数的定义，也就是实际方程的定义。考虑到速度和压力（编译器的优化将确保只做一次），考虑到介绍中所说的方程，这是直截了当的。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> 
 *     euler_flux(const Tensor<1, dim + 2, Number> &conserved_variables) 
 *   { 
 *     const Tensor<1, dim, Number> velocity = 
 *       euler_velocity<dim>(conserved_variables); 
 *     const Number pressure = euler_pressure<dim>(conserved_variables); 
 * 
 *     Tensor<1, dim + 2, Tensor<1, dim, Number>> flux; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       { 
 *         flux[0][d] = conserved_variables[1 + d]; 
 *         for (unsigned int e = 0; e < dim; ++e) 
 *           flux[e + 1][d] = conserved_variables[e + 1] * velocity[d]; 
 *         flux[d + 1][d] += pressure; 
 *         flux[dim + 1][d] = 
 *           velocity[d] * (conserved_variables[dim + 1] + pressure); 
 *       } 
 * 
 *     return flux; 
 *   } 
 * 
 * @endcode
 * 
 * 接下来的这个函数是一个简化数值通量实现的助手，它实现了一个张量的张量（具有大小为`dim + 2`的非标准外维，所以deal.II的张量类提供的标准重载在此不适用）与另一个相同内维的张量的作用，即一个矩阵-向量积。
 * 

 * 
 * 
 * @code
 *   template <int n_components, int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, n_components, Number> 
 *     operator*(const Tensor<1, n_components, Tensor<1, dim, Number>> &matrix, 
 *               const Tensor<1, dim, Number> &                         vector) 
 *   { 
 *     Tensor<1, n_components, Number> result; 
 *     for (unsigned int d = 0; d < n_components; ++d) 
 *       result[d] = matrix[d] * vector; 
 *     return result; 
 *   } 
 * 
 * @endcode
 * 
 * 这个函数实现了数值通量（黎曼求解器）。它从一个界面的两边获得状态，并获得法向量，从解的一边  $\mathbf{w}^-$  向解  $\mathbf{w}^+$  的方向。在依赖片断恒定数据的有限体积方法中，数值通量是核心成分，因为它是唯一输入物理信息的地方。在DG方法中，由于元素内部的多项式和那里使用的物理通量，数值通量就不那么核心了。由于在连续解的极限中，两边的数值一致的高阶插值，数值通量可以被看作是对两边解的跳跃的控制，以弱化连续性。必须认识到，在存在冲击的情况下，仅靠数值通量是无法稳定高阶DG方法的，因此任何DG方法都必须与进一步的冲击捕捉技术相结合，以处理这些情况。在本教程中，我们将重点讨论欧拉方程在没有强不连续的亚声速体系中的波状解，我们的基本方案已经足够了。
 * 

 * 
 * 尽管如此，数值通量对整个方案的数值耗散起着决定性作用，并影响到显式Runge-Kutta方法的可接受的时间步长。我们考虑两种选择，一种是改良的Lax-Friedrichs方案，另一种是广泛使用的Harten-Lax-van Leer（HLL）通量。对于这两种方案，我们首先需要得到界面两边的速度和压力，并评估物理欧拉通量。
 * 

 * 
 * 对于局部Lax--Friedrichs通量，其定义是 $\hat{\mathbf{F}}
 * =\frac{\mathbf{F}(\mathbf{w}^-)+\mathbf{F}(\mathbf{w}^+)}{2} +
 * \frac{\lambda}{2}\left[\mathbf{w}^--\mathbf{w}^+\right]\otimes
 * \mathbf{n^-}$  ，其中因子 $\lambda =
 * \max\left(\|\mathbf{u}^-\|+c^-, \|\mathbf{u}^+\|+c^+\right)$ 给出了最大波速， $c = \sqrt{\gamma p / \rho}$ 是音速。在这里，考虑到通量对解的影响很小，为了计算效率的原因，我们选择了该表达式的两个修改。对于上述因子 $\lambda$ 的定义，我们需要取四个平方根，两个用于两个速度规范，两个用于两侧的声速。因此，第一个修改是宁可使用 $\sqrt{\|\mathbf{u}\|^2+c^2}$ 作为最大速度的估计（如介绍中所示，它与实际最大速度最多相差2倍）。这使我们能够从最大速度中提取平方根，并且只需进行一次平方根计算就可以了。第二个修改是进一步放宽参数 $\lambda$ --它越小，耗散系数就越小（与 $\mathbf{w}$ 的跳跃相乘，最终可能导致耗散变小或变大）。这使得我们可以用更大的时间步长将频谱纳入显式Runge--Kutta积分器的稳定区域。然而，我们不能使耗散太小，因为否则假想的特征值会越来越大。最后，目前的保守公式在 $\lambda\to 0$ 的极限中不是能量稳定的，因为它不是偏斜对称的，在这种情况下需要额外的措施，如分裂形式的DG方案。
 * 

 * 
 * 对于HLL通量，我们遵循文献中的公式，通过一个参数 $s$ 引入Lax--Friedrichs的两个状态的额外加权。它是由欧拉方程的物理传输方向得出的，以当前的速度方向和声速为准。对于速度，我们在此选择一个简单的算术平均数，这对危险情况和材料参数的适度跳跃是足够的。
 * 

 * 
 * 由于数值通量在弱形式下是与法向量相乘的，因此我们对方程中的所有项都用法向量来乘以结果。在这些乘法中，上面定义的 "操作符*"可以实现类似于数学定义的紧凑符号。
 * 

 * 
 * 在这个函数和下面的函数中，我们使用变量后缀`_m`和`_p`来表示从 $\mathbf{w}^-$ 和 $\mathbf{w}^+$ 得出的量，即在观察相邻单元时相对于当前单元的 "这里 "和 "那里 "的数值。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename Number> 
 *   inline DEAL_II_ALWAYS_INLINE // 
 *     Tensor<1, dim + 2, Number> 
 *     euler_numerical_flux(const Tensor<1, dim + 2, Number> &u_m, 
 *                          const Tensor<1, dim + 2, Number> &u_p, 
 *                          const Tensor<1, dim, Number> &    normal) 
 *   { 
 *     const auto velocity_m = euler_velocity<dim>(u_m); 
 *     const auto velocity_p = euler_velocity<dim>(u_p); 
 * 
 *     const auto pressure_m = euler_pressure<dim>(u_m); 
 *     const auto pressure_p = euler_pressure<dim>(u_p); 
 * 
 *     const auto flux_m = euler_flux<dim>(u_m); 
 *     const auto flux_p = euler_flux<dim>(u_p); 
 * 
 *     switch (numerical_flux_type) 
 *       { 
 *         case lax_friedrichs_modified: 
 *           { 
 *             const auto lambda = 
 *               0.5 * std::sqrt(std::max(velocity_p.norm_square() + 
 *                                          gamma * pressure_p * (1. / u_p[0]), 
 *                                        velocity_m.norm_square() + 
 *                                          gamma * pressure_m * (1. / u_m[0]))); 
 * 
 *             return 0.5 * (flux_m * normal + flux_p * normal) + 
 *                    0.5 * lambda * (u_m - u_p); 
 *           } 
 * 
 *         case harten_lax_vanleer: 
 *           { 
 *             const auto avg_velocity_normal = 
 *               0.5 * ((velocity_m + velocity_p) * normal); 
 *             const auto   avg_c = std::sqrt(std::abs( 
 *               0.5 * gamma * 
 *               (pressure_p * (1. / u_p[0]) + pressure_m * (1. / u_m[0])))); 
 *             const Number s_pos = 
 *               std::max(Number(), avg_velocity_normal + avg_c); 
 *             const Number s_neg = 
 *               std::min(Number(), avg_velocity_normal - avg_c); 
 *             const Number inverse_s = Number(1.) / (s_pos - s_neg); 
 * 
 *             return inverse_s * 
 *                    ((s_pos * (flux_m * normal) - s_neg * (flux_p * normal)) - 
 *                     s_pos * s_neg * (u_m - u_p)); 
 *           } 
 * 
 *         default: 
 *           { 
 *             Assert(false, ExcNotImplemented()); 
 *             return {}; 
 *           } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 这个函数和下一个函数是辅助函数，提供紧凑的评估调用，因为多个点通过VectorizedArray参数被分批放在一起（详见 step-37 教程）。这个函数用于亚音速外流边界条件，我们需要将能量分量设置为一个规定值。下一个函数请求所有分量上的解，用于流入边界，其中解的所有分量都被设置。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename Number> 
 *   VectorizedArray<Number> 
 *   evaluate_function(const Function<dim> &                      function, 
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized, 
 *                     const unsigned int                         component) 
 *   { 
 *     VectorizedArray<Number> result; 
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
 *       { 
 *         Point<dim> p; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           p[d] = p_vectorized[d][v]; 
 *         result[v] = function.value(p, component); 
 *       } 
 *     return result; 
 *   } 
 * 
 *   template <int dim, typename Number, int n_components = dim + 2> 
 *   Tensor<1, n_components, VectorizedArray<Number>> 
 *   evaluate_function(const Function<dim> &                      function, 
 *                     const Point<dim, VectorizedArray<Number>> &p_vectorized) 
 *   { 
 *     AssertDimension(function.n_components, n_components); 
 *     Tensor<1, n_components, VectorizedArray<Number>> result; 
 *     for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
 *       { 
 *         Point<dim> p; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           p[d] = p_vectorized[d][v]; 
 *         for (unsigned int d = 0; d < n_components; ++d) 
 *           result[d][v] = function.value(p, d); 
 *       } 
 *     return result; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerOperationclass"></a> 
 * <h3>The EulerOperation class</h3>
 * 

 * 
 * 这个类实现了欧拉问题的评估器，类似于  step-37  或  step-59  的 `LaplaceOperator` 类。由于本算子是非线性的，不需要矩阵接口（交给预处理程序），我们跳过了无矩阵算子中的各种`vmult`函数，只实现了`apply`函数以及`apply`与上述低存储Runge-Kutta时间积分器所需的矢量更新的组合（称为`perform_stage`）。此外，我们还增加了三个涉及无矩阵例程的额外函数，即一个是根据元素中的速度和声速计算时间步长的估计值（与实际时间步长的Courant数相结合），一个是解的投影（专门针对DG情况的 VectorTools::project() ），还有一个是计算与可能的分析解或与某些背景状态的规范的误差。
 * 

 * 
 * 该课的其余部分与其他无矩阵教程相似。正如介绍中所讨论的，我们提供了几个函数，允许用户在由 types::boundary_id 变量标记的领域边界的不同部分传递各种形式的边界条件，以及可能的体力。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   class EulerOperator 
 *   { 
 *   public: 
 *     static constexpr unsigned int n_quadrature_points_1d = n_points_1d; 
 * 
 *     EulerOperator(TimerOutput &timer_output); 
 * 
 *     void reinit(const Mapping<dim> &   mapping, 
 *                 const DoFHandler<dim> &dof_handler); 
 * 
 *     void set_inflow_boundary(const types::boundary_id       boundary_id, 
 *                              std::unique_ptr<Function<dim>> inflow_function); 
 * 
 *     void set_subsonic_outflow_boundary( 
 *       const types::boundary_id       boundary_id, 
 *       std::unique_ptr<Function<dim>> outflow_energy); 
 * 
 *     void set_wall_boundary(const types::boundary_id boundary_id); 
 * 
 *     void set_body_force(std::unique_ptr<Function<dim>> body_force); 
 * 
 *     void apply(const double                                      current_time, 
 *                const LinearAlgebra::distributed::Vector<Number> &src, 
 *                LinearAlgebra::distributed::Vector<Number> &      dst) const; 
 * 
 *     void 
 *     perform_stage(const Number cur_time, 
 *                   const Number factor_solution, 
 *                   const Number factor_ai, 
 *                   const LinearAlgebra::distributed::Vector<Number> &current_ri, 
 *                   LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
 *                   LinearAlgebra::distributed::Vector<Number> &      solution, 
 *                   LinearAlgebra::distributed::Vector<Number> &next_ri) const; 
 * 
 *     void project(const Function<dim> &                       function, 
 *                  LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     std::array<double, 3> compute_errors( 
 *       const Function<dim> &                             function, 
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     double compute_cell_transport_speed( 
 *       const LinearAlgebra::distributed::Vector<Number> &solution) const; 
 * 
 *     void 
 *     initialize_vector(LinearAlgebra::distributed::Vector<Number> &vector) const; 
 * 
 *   private: 
 *     MatrixFree<dim, Number> data; 
 * 
 *     TimerOutput &timer; 
 * 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
 *       inflow_boundaries; 
 *     std::map<types::boundary_id, std::unique_ptr<Function<dim>>> 
 *                                    subsonic_outflow_boundaries; 
 *     std::set<types::boundary_id>   wall_boundaries; 
 *     std::unique_ptr<Function<dim>> body_force; 
 * 
 *     void local_apply_inverse_mass_matrix( 
 *       const MatrixFree<dim, Number> &                   data, 
 *       LinearAlgebra::distributed::Vector<Number> &      dst, 
 *       const LinearAlgebra::distributed::Vector<Number> &src, 
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const; 
 * 
 *     void local_apply_cell( 
 *       const MatrixFree<dim, Number> &                   data, 
 *       LinearAlgebra::distributed::Vector<Number> &      dst, 
 *       const LinearAlgebra::distributed::Vector<Number> &src, 
 *       const std::pair<unsigned int, unsigned int> &     cell_range) const; 
 * 
 *     void local_apply_face( 
 *       const MatrixFree<dim, Number> &                   data, 
 *       LinearAlgebra::distributed::Vector<Number> &      dst, 
 *       const LinearAlgebra::distributed::Vector<Number> &src, 
 *       const std::pair<unsigned int, unsigned int> &     face_range) const; 
 * 
 *     void local_apply_boundary_face( 
 *       const MatrixFree<dim, Number> &                   data, 
 *       LinearAlgebra::distributed::Vector<Number> &      dst, 
 *       const LinearAlgebra::distributed::Vector<Number> &src, 
 *       const std::pair<unsigned int, unsigned int> &     face_range) const; 
 *   }; 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   EulerOperator<dim, degree, n_points_1d>::EulerOperator(TimerOutput &timer) 
 *     : timer(timer) 
 *   {} 
 * 
 * @endcode
 * 
 * 对于欧拉算子的初始化，我们设置了类中包含的MatrixFree变量。这可以通过给定一个描述可能的弯曲边界的映射以及一个描述自由度的DoFHandler对象来完成。由于我们在这个教程程序中使用的是不连续的Galerkin离散化，没有对解场施加强烈的约束，所以我们不需要传入AffineConstraints对象，而是使用一个假的来构造。关于正交，我们要选择两种不同的方式来计算基础积分。第一种是灵活的，基于模板参数`n_points_1d`（将被分配到本文件顶部指定的`n_q_points_1d`值）。更精确的积分是必要的，以避免由于欧拉算子中的可变系数而产生的混叠问题。第二个不太精确的正交公式是一个基于`fe_degree+1`的严密公式，需要用于反质量矩阵。虽然该公式只在仿生元素形状上提供了精确的反，而在变形元素上则没有，但它可以通过张量积技术快速反转质量矩阵，这对于确保整体的最佳计算效率是必要的。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::reinit( 
 *     const Mapping<dim> &   mapping, 
 *     const DoFHandler<dim> &dof_handler) 
 *   { 
 *     const std::vector<const DoFHandler<dim> *> dof_handlers = {&dof_handler}; 
 *     const AffineConstraints<double>            dummy; 
 *     const std::vector<const AffineConstraints<double> *> constraints = {&dummy}; 
 *     const std::vector<Quadrature<1>> quadratures = {QGauss<1>(n_q_points_1d), 
 *                                                     QGauss<1>(fe_degree + 1)}; 
 * 
 *     typename MatrixFree<dim, Number>::AdditionalData additional_data; 
 *     additional_data.mapping_update_flags = 
 *       (update_gradients | update_JxW_values | update_quadrature_points | 
 *        update_values); 
 *     additional_data.mapping_update_flags_inner_faces = 
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors | 
 *        update_values); 
 *     additional_data.mapping_update_flags_boundary_faces = 
 *       (update_JxW_values | update_quadrature_points | update_normal_vectors | 
 *        update_values); 
 *     additional_data.tasks_parallel_scheme = 
 *       MatrixFree<dim, Number>::AdditionalData::none; 
 * 
 *     data.reinit( 
 *       mapping, dof_handlers, constraints, quadratures, additional_data); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::initialize_vector( 
 *     LinearAlgebra::distributed::Vector<Number> &vector) const 
 *   { 
 *     data.initialize_dof_vector(vector); 
 *   } 
 * 
 * @endcode
 * 
 * 随后的四个成员函数是必须从外部调用的，以指定各种类型的边界。对于一个流入的边界，我们必须以密度  $\rho$  、动量  $\rho \mathbf{u}$  和能量  $E$  来指定所有成分。考虑到这些信息，我们将函数与各自的边界ID一起存储在这个类的地图成员变量中。同样，我们对亚音速外流边界（我们也要求一个函数，用来检索能量）和壁面（无穿透）边界进行处理，在壁面上我们施加零法线速度（不需要函数，所以我们只要求边界ID）。对于目前的DG代码来说，边界条件只作为弱形式的一部分被应用（在时间积分期间），设置边界条件的调用可以出现在对这个类的`reinit()`调用之前或之后。这与连续有限元代码不同，在连续有限元代码中，边界条件决定了被送入MatrixFree初始化的AffineConstraints对象的内容，因此需要在无矩阵数据结构的初始化之前设置。
 * 

 * 
 * 在四个函数中的每一个中添加的检查是用来确保边界条件在边界的各个部分是相互排斥的，也就是说，用户不会意外地将一个边界既指定为流入边界，又指定为亚声速流出边界。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_inflow_boundary( 
 *     const types::boundary_id       boundary_id, 
 *     std::unique_ptr<Function<dim>> inflow_function) 
 *   { 
 *     AssertThrow(subsonic_outflow_boundaries.find(boundary_id) == 
 *                     subsonic_outflow_boundaries.end() && 
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
 *  
 *  
 *  
 *  
 *  
 *                 ExcMessage("Expected function with dim+2 components")); 
 * 
 *     inflow_boundaries[boundary_id] = std::move(inflow_function); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_subsonic_outflow_boundary( 
 *     const types::boundary_id       boundary_id, 
 *     std::unique_ptr<Function<dim>> outflow_function) 
 *   { 
 *     AssertThrow(inflow_boundaries.find(boundary_id) == 
 *                     inflow_boundaries.end() && 
 *                   wall_boundaries.find(boundary_id) == wall_boundaries.end(), 
 *                 ExcMessage("You already set the boundary with id " + 
 *                            std::to_string(static_cast<int>(boundary_id)) + 
 *                            " to another type of boundary before now setting " + 
 *                            "it as subsonic outflow")); 
 *     AssertThrow(outflow_function->n_components == dim + 2, 
 *                 ExcMessage("Expected function with dim+2 components")); 
 * 
 *     subsonic_outflow_boundaries[boundary_id] = std::move(outflow_function); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_wall_boundary( 
 *     const types::boundary_id boundary_id) 
 *   { 
 *     AssertThrow(inflow_boundaries.find(boundary_id) == 
 *                     inflow_boundaries.end() && 
 *                   subsonic_outflow_boundaries.find(boundary_id) == 
 *                     subsonic_outflow_boundaries.end(), 
 *                 ExcMessage("You already set the boundary with id " + 
 *                            std::to_string(static_cast<int>(boundary_id)) + 
 *                            " to another type of boundary before now setting " + 
 *                            "it as wall boundary")); 
 * 
 *     wall_boundaries.insert(boundary_id); 
 *   } 
 * 
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::set_body_force( 
 *     std::unique_ptr<Function<dim>> body_force) 
 *   { 
 *     AssertDimension(body_force->n_components, dim); 
 * 
 *     this->body_force = std::move(body_force); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Localevaluators"></a> 
 * <h4>Local evaluators</h4>
 * 

 * 
 * 现在我们开始研究欧拉问题的局部评估器。评估器相对简单，遵循  step-37  、  step-48  或  step-59  中提出的内容。第一个显著的区别是，我们使用的是具有非标准正交点数量的FEE评估。以前我们总是将正交点的数量设置为等于多项式度数加1（确保在仿生元素形状上的精确积分），现在我们将正交点的数量设置为一个单独的变量（例如多项式度数加多项式度数的二分之一或三分之一），以更准确地处理非线性项。由于评估器通过模板参数输入了适当的循环长度，并在变量 FEEvaluation::n_q_points, 中保留了整个单元格的正交点数量，所以我们现在自动操作更精确的公式，而无需进一步修改。
 * 

 * 
 * 第二个区别是由于我们现在评估的是一个多分量系统，而不是之前考虑的标量系统。无矩阵框架提供了几种方法来处理多成分的情况。这里显示的变体是利用一个嵌入了多个分量的FEEvaluation对象，由第四个模板参数`dim + 2`指定欧拉系统中的分量。因此， FEEvaluation::get_value() 的返回类型不再是一个标量（这将返回一个VectorizedArray类型，收集几个元素的数据），而是一个`dim+2`组件的张量。该功能与标量的情况类似；它由一个基类的模板专业化处理，称为FEEvaluationAccess。另一个变体是使用几个FEEvaluation对象，一个标量对象用于密度，一个带`dim`分量的矢量值对象用于动量，另一个标量评价器用于能量。为了确保这些分量指向解决方案的正确部分，FEEvaluation的构造函数在所需的MatrixFree字段之后需要三个可选的整数参数，即多DoFHandler系统的DoFHandler编号（默认取第一个），如果有多个Quadrature对象，则取正交点的编号（见下文），以及作为第三个参数的矢量系统中的分量。由于我们有一个单一的矢量来表示所有的分量，我们将使用第三个参数，并将其设置为`0`表示密度，`1`表示矢量值的动量，`dim+1`表示能量槽。然后FEEvaluation在 FEEvaluationBase::read_dof_values() 和 FEEvaluation::distributed_local_to_global() 或更紧凑的 FEEvaluation::gather_evaluate() 和 FEEvaluation::integrate_scatter() 调用中挑选适当的解矢量子范围。
 * 

 * 
 * 当涉及到身体力向量的评估时，为了效率，我们区分了两种情况。如果我们有一个常数函数（源自 Functions::ConstantFunction), ），我们可以在正交点的循环外预先计算出数值，并简单地在所有地方使用该数值。对于一个更通用的函数，我们反而需要调用我们上面提供的`evaluate_function()`方法；这个路径更昂贵，因为我们需要访问与正交点数据有关的内存。
 * 

 * 
 * 其余部分沿用其他教程的程序。由于我们已经在单独的`euler_flux()`函数中实现了欧拉方程的所有物理学，我们在这里所要做的就是给定在正交点评估的当前解，由`phi.get_value(q)`返回，并告诉FEEvaluation对象，通过形状函数的梯度（这是一个外部`dim+2`分量的张量，每个张量持有一个`dim`分量的 $x,y,z$  ] 欧拉通量的分量）。) 最后值得一提的是，在我们得到一个外部函数的情况下，我们通过测试函数`phi.submit_value()`的值来排队测试数据的顺序。我们必须在调用`phi.get_value(q)'之后进行，因为`get_value()'（读取解决方案）和`submit_value()'（排队等待测试函数的乘法和正交点的求和）访问同一个底层数据域。这里很容易实现没有临时变量`w_q`，因为值和梯度之间没有混合。对于更复杂的设置，必须首先复制出例如正交点的值和梯度，然后通过 FEEvaluationBase::submit_value() 和 FEEvaluationBase::submit_gradient(). 再次排列结果。
 * 

 * 
 * 作为最后的说明，我们提到我们没有使用这个函数的第一个MatrixFree参数，这是一个来自 MatrixFree::loop(). 的回调，接口规定了现在的参数列表，但是由于我们在一个成员函数中，MatrixFree对象已经可以作为`data`变量，我们坚持使用，以避免混淆。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_cell( 
 *     const MatrixFree<dim, Number> &, 
 *     LinearAlgebra::distributed::Vector<Number> &      dst, 
 *     const LinearAlgebra::distributed::Vector<Number> &src, 
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const 
 *   { 
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data); 
 * 
 *     Tensor<1, dim, VectorizedArray<Number>> constant_body_force; 
 *     const Functions::ConstantFunction<dim> *constant_function = 
 *       dynamic_cast<Functions::ConstantFunction<dim> *>(body_force.get()); 
 * 
 *     if (constant_function) 
 *       constant_body_force = evaluate_function<dim, Number, dim>( 
 *         *constant_function, Point<dim, VectorizedArray<Number>>()); 
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.gather_evaluate(src, EvaluationFlags::values); 
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto w_q = phi.get_value(q); 
 *             phi.submit_gradient(euler_flux<dim>(w_q), q); 
 *             if (body_force.get() != nullptr) 
 *               { 
 *                 const Tensor<1, dim, VectorizedArray<Number>> force = 
 *                   constant_function ? constant_body_force : 
 *                                       evaluate_function<dim, Number, dim>( 
 *                                         *body_force, phi.quadrature_point(q)); 
 * 
 *                 Tensor<1, dim + 2, VectorizedArray<Number>> forcing; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   forcing[d + 1] = w_q[0] * force[d]; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   forcing[dim + 1] += force[d] * w_q[d + 1]; 
 * 
 *                 phi.submit_value(forcing, q); 
 *               } 
 *           } 
 * 
 *         phi.integrate_scatter(((body_force.get() != nullptr) ? 
 *                                  EvaluationFlags::values : 
 *                                  EvaluationFlags::nothing) | 
 *                                 EvaluationFlags::gradients, 
 *                               dst); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数涉及到内部面的积分计算，在这里我们需要与面相邻的两个单元的评估器。我们将变量`phi_m`与解分量 $\mathbf{w}^-$ 相关联，将变量`phi_p`与解分量 $\mathbf{w}^+$ 相关联。我们在FEFaceEvaluation的构造函数中通过第二个参数来区分两边，`true`表示内侧，`false`表示外侧，内侧和外侧表示相对于法向量的方向。
 * 

 * 
 * 注意调用 FEFaceEvaluation::gather_evaluate() 和 FEFaceEvaluation::integrate_scatter() 结合了对向量的访问和因式分解部分。这种合并操作不仅节省了一行代码，而且还包含了一个重要的优化。鉴于我们在Gauss-Lobatto正交公式的点上使用拉格朗日多项式的节点基础，在每个面上只有 $(p+1)^{d-1}$ 的基础函数评估为非零。因此，评估器只访问了向量中的必要数据，而跳过了乘以零的部分。如果我们首先读取向量，我们就需要从向量中加载所有的数据，因为孤立的调用不知道后续操作中需要哪些数据。如果随后的 FEFaceEvaluation::evaluate() 调用要求数值和导数，确实需要每个分量的所有 $(p+1)^d$ 向量条目，因为所有基函数的法向导数都是非零的。
 * 

 * 
 * 评价器的参数以及程序与单元评价相似。由于非线性项的存在，我们再次使用更精确的（过度）积分方案，指定为列表中第三个模板参数。在正交点上，我们再去找我们的自由函数来计算数值通量。它从两边（即 $\mathbf{w}^-$ 和 $\mathbf{w}^+$ ）接收在正交点评估的解决方案，以及到减去一边的法向量。正如上面所解释的，数值通量已经乘以来自减法侧的法向量了。我们需要转换符号，因为在引言中得出的弱形式中，边界项带有一个减号。然后，通量被排队在减号和加号上进行测试，由于加号上的法向量与减号上的法向量正好相反，所以要调换符号。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_face( 
 *     const MatrixFree<dim, Number> &, 
 *     LinearAlgebra::distributed::Vector<Number> &      dst, 
 *     const LinearAlgebra::distributed::Vector<Number> &src, 
 *     const std::pair<unsigned int, unsigned int> &     face_range) const 
 *   { 
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_m(data, 
 *                                                                       true); 
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi_p(data, 
 *                                                                       false); 
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face) 
 *       { 
 *         phi_p.reinit(face); 
 *         phi_p.gather_evaluate(src, EvaluationFlags::values); 
 * 
 *         phi_m.reinit(face); 
 *         phi_m.gather_evaluate(src, EvaluationFlags::values); 
 * 
 *         for (unsigned int q = 0; q < phi_m.n_q_points; ++q) 
 *           { 
 *             const auto numerical_flux = 
 *               euler_numerical_flux<dim>(phi_m.get_value(q), 
 *                                         phi_p.get_value(q), 
 *                                         phi_m.get_normal_vector(q)); 
 *             phi_m.submit_value(-numerical_flux, q); 
 *             phi_p.submit_value(numerical_flux, q); 
 *           } 
 * 
 *         phi_p.integrate_scatter(EvaluationFlags::values, dst); 
 *         phi_m.integrate_scatter(EvaluationFlags::values, dst); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 对于位于边界的面，我们需要施加适当的边界条件。在这个教程程序中，我们实现了上述的四种情况。第五种情况，即超音速流出条件，将在下面的 "结果 "部分讨论）。不连续的Galerkin方法对边界条件的施加不是作为约束条件，而只是弱化。因此，各种条件是通过找到一个适当的<i>exterior</i>量 $\mathbf{w}^+$ 来施加的，然后将其交给也用于内部面的数值通量函数。实质上，我们在域外 "假装 "一个状态，如果那是现实，PDE的解将满足我们想要的边界条件。
 * 

 * 
 * 对于墙的边界，我们需要对动量变量施加一个无正态通量的条件，而对于密度和能量，我们使用的是诺伊曼条件  $\rho^+ = \rho^-$  和  $E^+ = E^-$  。为了实现无正态通量条件，我们将外部数值设定为内部数值，并减去墙面法线方向，即法线矢量方向上的速度的2倍。
 * 

 * 
 * 对于流入边界，我们简单地将给定的Dirichlet数据 $\mathbf{w}_\mathrm{D}$ 作为边界值。另一种方法是使用 $\mathbf{w}^+ = -\mathbf{w}^- + 2 \mathbf{w}_\mathrm{D}$  ，即所谓的镜像原理。
 * 

 * 
 * 强加外流本质上是一个诺伊曼条件，即设定  $\mathbf{w}^+ = \mathbf{w}^-$  。对于亚声速流出的情况，我们仍然需要强加一个能量值，我们从各自的函数中得出这个值。对于<i>backflow</i>的情况，即在Neumann部分有动量通入域的情况，需要一个特殊的步骤。根据文献（这一事实可以通过适当的能量论证得出），我们必须切换到流入部分的通量的另一个变体，见Gravemeier, Comerford, Yoshihara, Ismail, Wall, "A novel formulation for Neumann inflow conditions in biomechanics", Int. J. Numer. Meth. 生物医学。Eng., vol. 28 (2012). 这里，动量项需要再次添加，这相当于去除动量变量上的通量贡献。我们在后处理步骤中这样做，而且只适用于我们都处于外流边界且法向量与动量（或等同于速度）之间的点积为负的情况。由于我们在SIMD矢量化中一次处理多个正交点的数据，这里需要明确地在SIMD数组的条目上循环。
 * 

 * 
 * 在下面的实现中，我们在正交点的层面上检查各种类型的边界。当然，我们也可以将决定权移出正交点循环，将整个面孔视为同类，这就避免了在正交点的内循环中进行一些地图/集合的查找。然而，效率的损失并不明显，所以我们在这里选择了更简单的代码。还要注意的是，最后的 "else "子句会捕捉到这样的情况，即边界的某些部分没有通过 `EulerOperator::set_..._boundary(...)`. 分配任何边界条件。
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_boundary_face( 
 *     const MatrixFree<dim, Number> &, 
 *     LinearAlgebra::distributed::Vector<Number> &      dst, 
 *     const LinearAlgebra::distributed::Vector<Number> &src, 
 *     const std::pair<unsigned int, unsigned int> &     face_range) const 
 *   { 
 *     FEFaceEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, true); 
 * 
 *     for (unsigned int face = face_range.first; face < face_range.second; ++face) 
 *       { 
 *         phi.reinit(face); 
 *         phi.gather_evaluate(src, EvaluationFlags::values); 
 * 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto w_m    = phi.get_value(q); 
 *             const auto normal = phi.get_normal_vector(q); 
 * 
 *             auto rho_u_dot_n = w_m[1] * normal[0]; 
 *             for (unsigned int d = 1; d < dim; ++d) 
 *               rho_u_dot_n += w_m[1 + d] * normal[d]; 
 * 
 *             bool at_outflow = false; 
 * 
 *             Tensor<1, dim + 2, VectorizedArray<Number>> w_p; 
 *             const auto boundary_id = data.get_boundary_id(face); 
 *             if (wall_boundaries.find(boundary_id) != wall_boundaries.end()) 
 *               { 
 *                 w_p[0] = w_m[0]; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   w_p[d + 1] = w_m[d + 1] - 2. * rho_u_dot_n * normal[d]; 
 *                 w_p[dim + 1] = w_m[dim + 1]; 
 *               } 
 *             else if (inflow_boundaries.find(boundary_id) != 
 *                      inflow_boundaries.end()) 
 *               w_p = 
 *                 evaluate_function(*inflow_boundaries.find(boundary_id)->second, 
 *                                   phi.quadrature_point(q)); 
 *             else if (subsonic_outflow_boundaries.find(boundary_id) != 
 *                      subsonic_outflow_boundaries.end()) 
 *               { 
 *                 w_p          = w_m; 
 *                 w_p[dim + 1] = evaluate_function( 
 *                   *subsonic_outflow_boundaries.find(boundary_id)->second, 
 *                   phi.quadrature_point(q), 
 *                   dim + 1); 
 *                 at_outflow = true; 
 *               } 
 *             else 
 *               AssertThrow(false, 
 *                           ExcMessage("Unknown boundary id, did " 
 *                                      "you set a boundary condition for " 
 *                                      "this part of the domain boundary?")); 
 * 
 *             auto flux = euler_numerical_flux<dim>(w_m, w_p, normal); 
 * 
 *             if (at_outflow) 
 *               for (unsigned int v = 0; v < VectorizedArray<Number>::size(); ++v) 
 *                 { 
 *                   if (rho_u_dot_n[v] < -1e-12) 
 *                     for (unsigned int d = 0; d < dim; ++d) 
 *                       flux[d + 1][v] = 0.; 
 *                 } 
 * 
 *             phi.submit_value(-flux, q); 
 *           } 
 * 
 *         phi.integrate_scatter(EvaluationFlags::values, dst); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数实现了质量矩阵的逆运算。在介绍中已经广泛讨论了算法和原理，所以我们在这里只讨论 MatrixFreeOperators::CellwiseInverseMassMatrix 类的技术问题。它所做的操作与质量矩阵的正向评估类似，只是使用了不同的插值矩阵，代表逆 $S^{-1}$ 因子。这些代表了从指定的基础（在这种情况下，高斯--洛巴托正交公式点中的拉格朗日基础）到高斯正交公式点中的拉格朗日基础的改变。在后者的基础上，我们可以应用点的逆向`JxW`因子，即正交权重乘以从参考坐标到实坐标的映射的雅各布系数。一旦完成了这一操作，基数将再次变回节点高斯-洛巴托基数。所有这些操作都由下面的 "apply() "函数完成。我们需要提供的是要操作的局部场（我们通过一个FEEvaluation对象从全局向量中提取），并将结果写回质量矩阵操作的目标向量。
 * 

 * 
 * 需要注意的一点是，我们在FEEvaluation的构造函数中添加了两个整数参数（可选），第一个是0（在多DoFHandler系统中选择DoFHandler；在这里，我们只有一个），第二个是1，用于进行正交公式选择。由于我们将正交公式0用于非线性项的过度积分，我们使用公式1与默认的 $p+1$ （或变量名称中的`fe_degree+1`）点用于质量矩阵。这导致了对质量矩阵的平方贡献，并确保了精确的积分，正如介绍中所解释的。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::local_apply_inverse_mass_matrix( 
 *     const MatrixFree<dim, Number> &, 
 *     LinearAlgebra::distributed::Vector<Number> &      dst, 
 *     const LinearAlgebra::distributed::Vector<Number> &src, 
 *     const std::pair<unsigned int, unsigned int> &     cell_range) const 
 *   { 
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number> 
 *       inverse(phi); 
 * 
 *     for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.read_dof_values(src); 
 * 
 *         inverse.apply(phi.begin_dof_values(), phi.begin_dof_values()); 
 * 
 *         phi.set_dof_values(dst); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Theapplyandrelatedfunctions"></a> 
 * <h4>The apply() and related functions</h4>
 * 

 * 
 * 我们现在来到实现欧拉算子整体评估的函数，即 $\mathcal M^{-1} \mathcal L(t, \mathbf{w})$  ，调用上面介绍的局部评估器。这些步骤在前面的代码中应该是清楚的。需要注意的一点是，我们需要调整与边界各部分相关的函数中的时间，以便在边界数据与时间相关的情况下与方程一致。然后，我们调用 MatrixFree::loop() 来执行单元和面的积分，包括在`src`向量中进行必要的ghost数据交换。该函数的第七个参数，"true"，指定我们要在开始向其累积积分之前，将 "dst "向量作为循环的一部分归零。这个变体比在循环之前明确调用`dst = 0.;`要好，因为归零操作是在矢量的子范围内完成的，其部分是由附近的积分写入的。这加强了数据的定位，并允许缓存，节省了向量数据到主内存的一次往返，提高了性能。循环的最后两个参数决定了哪些数据被交换：由于我们只访问一个面的形状函数的值，这是典型的一阶双曲问题，并且由于我们有一个节点基础，节点位于参考元素表面，我们只需要交换这些部分。这又节省了宝贵的内存带宽。
 * 

 * 
 * 一旦应用了空间算子 $\mathcal L$ ，我们需要进行第二轮操作，应用反质量矩阵。这里，我们调用 MatrixFree::cell_loop() ，因为只有单元格积分出现。单元循环比全循环更便宜，因为只访问与本地拥有的单元相关的自由度，这只是DG离散化的本地拥有的自由度。因此，这里不需要鬼魂交换。
 * 

 * 
 * 在所有这些函数的周围，我们设置了定时器范围来记录计算时间，以统计各部分的贡献。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::apply( 
 *     const double                                      current_time, 
 *     const LinearAlgebra::distributed::Vector<Number> &src, 
 *     LinearAlgebra::distributed::Vector<Number> &      dst) const 
 *   { 
 *     { 
 *       TimerOutput::Scope t(timer, "apply - integrals"); 
 * 
 *       for (auto &i : inflow_boundaries) 
 *         i.second->set_time(current_time); 
 *       for (auto &i : subsonic_outflow_boundaries) 
 *         i.second->set_time(current_time); 
 * 
 *       data.loop(&EulerOperator::local_apply_cell, 
 *                 &EulerOperator::local_apply_face, 
 *                 &EulerOperator::local_apply_boundary_face, 
 *                 this, 
 *                 dst, 
 *                 src, 
 *                 true, 
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values, 
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values); 
 *     } 
 * 
 *     { 
 *       TimerOutput::Scope t(timer, "apply - inverse mass"); 
 * 
 *       data.cell_loop(&EulerOperator::local_apply_inverse_mass_matrix, 
 *                      this, 
 *                      dst, 
 *                      dst); 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 让我们转到做Runge--Kutta更新的整个阶段的函数。它调用 EulerOperator::apply() ，然后对向量进行一些更新，即`next_ri = solution + factor_ai * k_i`和`solution += factor_solution * k_i`。与其通过向量接口执行这些步骤，我们在这里提出了一个替代策略，在基于缓存的架构上速度更快。由于向量所消耗的内存往往比缓存所能容纳的要大得多，因此数据必须有效地来自缓慢的RAM内存。这种情况可以通过循环融合来改善，即在一次扫描中对`next_ki`和`solution`进行更新。在这种情况下，我们将读取两个向量`rhs`和`solution`并写入`next_ki`和`solution`，而在基线情况下，至少有4次读取和两次写入。在这里，我们更进一步，当质量矩阵反转在向量的某一部分完成后，立即执行循环。  MatrixFree::cell_loop() 提供了一种机制，在单元格的循环第一次接触到一个向量条目之前，附加一个 `std::function` （我们在这里没有使用，但用于例如向量的归零），以及在循环最后接触到一个条目之后，调用第二个 `std::function` 。回调的形式是给定向量上的一个范围（就MPI宇宙中的本地索引编号而言），可以由`local_element()`函数来处理。
 * 

 * 
 * 对于这个第二个回调，我们创建一个lambda，在一个范围内工作，并在这个范围内写入相应的更新。理想情况下，我们会在本地循环之前添加`DEAL_II_OPENMP_SIMD_PRAGMA`，以建议编译器对这个循环进行SIMD并行化（这意味着在实践中我们要确保在循环内部使用的指针的索引范围之间没有重叠，也称为别名）。事实证明，在写这篇文章的时候，GCC 7.2无法编译lambda函数中的OpenMP pragma，所以我们在下面注释了这个pragma。如果你的编译器比较新，你应该可以再次取消注释这些行。
 * 

 * 
 * 注意，当我们不需要更新`next_ri`向量时，我们为最后的Runge--Kutta阶段选择不同的代码路径。这个策略带来了相当大的速度提升。在40核机器上，默认矢量更新时，逆质量矩阵和矢量更新需要60%以上的计算时间，而在更优化的变体中，这一比例约为35%。换句话说，这是一个大约三分之一的速度提升。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::perform_stage( 
 *     const Number                                      current_time, 
 *     const Number                                      factor_solution, 
 *     const Number                                      factor_ai, 
 *     const LinearAlgebra::distributed::Vector<Number> &current_ri, 
 *     LinearAlgebra::distributed::Vector<Number> &      vec_ki, 
 *     LinearAlgebra::distributed::Vector<Number> &      solution, 
 *     LinearAlgebra::distributed::Vector<Number> &      next_ri) const 
 *   { 
 *     { 
 *       TimerOutput::Scope t(timer, "rk_stage - integrals L_h"); 
 * 
 *       for (auto &i : inflow_boundaries) 
 *         i.second->set_time(current_time); 
 *       for (auto &i : subsonic_outflow_boundaries) 
 *         i.second->set_time(current_time); 
 * 
 *       data.loop(&EulerOperator::local_apply_cell, 
 *                 &EulerOperator::local_apply_face, 
 *                 &EulerOperator::local_apply_boundary_face, 
 *                 this, 
 *                 vec_ki, 
 *                 current_ri, 
 *                 true, 
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values, 
 *                 MatrixFree<dim, Number>::DataAccessOnFaces::values); 
 *     } 
 * 
 *     { 
 *       TimerOutput::Scope t(timer, "rk_stage - inv mass + vec upd"); 
 *       data.cell_loop( 
 *         &EulerOperator::local_apply_inverse_mass_matrix, 
 *         this, 
 *         next_ri, 
 *         vec_ki, 
 *         std::function<void(const unsigned int, const unsigned int)>(), 
 *         [&](const unsigned int start_range, const unsigned int end_range) { 
 *           const Number ai = factor_ai; 
 *           const Number bi = factor_solution; 
 *           if (ai == Number()) 
 *             { 
 * 
 *           /* DEAL_II_OPENMP_SIMD_PRAGMA  */ 
 *               for (unsigned int i = start_range; i < end_range; ++i) 
 *                 { 
 *                   const Number k_i          = next_ri.local_element(i); 
 *                   const Number sol_i        = solution.local_element(i); 
 *                   solution.local_element(i) = sol_i + bi * k_i; 
 *                 } 
 *             } 
 *           else 
 *             { 
 * 
 *               /* DEAL_II_OPENMP_SIMD_PRAGMA  */ 
 *               for (unsigned int i = start_range; i < end_range; ++i) 
 *                 { 
 *                   const Number k_i          = next_ri.local_element(i); 
 *                   const Number sol_i        = solution.local_element(i); 
 *                   solution.local_element(i) = sol_i + bi * k_i; 
 *                   next_ri.local_element(i)  = sol_i + ai * k_i; 
 *                 } 
 *             } 
 *         }); 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 在讨论了将解提前一个时间步长的函数的实现后，现在让我们来看看实现其他辅助性操作的函数。具体来说，这些是计算投影、评估误差和计算单元上信息传输速度的函数。
 * 

 * 
 * 这些函数中的第一个基本上等同于 VectorTools::project(), ，只是速度快得多，因为它是专门针对DG元素的，不需要设置和解决线性系统，因为每个元素都有独立的基函数。我们在这里展示代码的原因，除了这个非关键操作的小幅提速之外，还因为它显示了 MatrixFreeOperators::CellwiseInverseMassMatrix. 提供的额外功能。
 * 

 * 
 * 投影操作的工作原理如下。如果我们用 $S$ 表示在正交点评估的形状函数矩阵，那么在单元格 $K$ 上的投影是一个形式为 $\underbrace{S J^K S^\mathrm T}_{\mathcal M^K} \mathbf{w}^K = S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ 的操作，其中 $J^K$ 是包含雅各布系数乘以正交权重（JxW）的对角矩阵， $\mathcal M^K$ 是单元格的质量矩阵， $\tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$ 是要投影到正交点的领域评估。实际上，矩阵 $S$ 通过张量积有额外的结构，如介绍中所解释的）。这个系统现在可以等效地写成 $\mathbf{w}^K = \left(S J^K S^\mathrm T\right)^{-1} S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q} = S^{-\mathrm T} \left(J^K\right)^{-1} S^{-1} S J^K \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$  。现在，项 $S^{-1} S$ 和 $\left(J^K\right)^{-1} J^K$ 相抵消，导致最后的表达式 $\mathbf{w}^K = S^{-\mathrm T} \tilde{\mathbf{w}}(\mathbf{x}_q)_{q=1:n_q}$  。这个操作由 MatrixFreeOperators::CellwiseInverseMassMatrix::transform_from_q_points_to_basis(). 实现。这个名字来自于这个投影只是乘以 $S^{-\mathrm T}$ ，一个从高斯正交点的节点基到给定的有限元基的基数变化。请注意，我们调用 FEEvaluation::set_dof_values() 将结果写入矢量，覆盖之前的内容，而不是像典型的积分任务那样累积结果--我们可以这样做，因为对于不连续的Galerkin离散，每个矢量条目都只有一个单元的贡献。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   void EulerOperator<dim, degree, n_points_1d>::project( 
 *     const Function<dim> &                       function, 
 *     LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 
 *     MatrixFreeOperators::CellwiseInverseMassMatrix<dim, degree, dim + 2, Number> 
 *       inverse(phi); 
 *     solution.zero_out_ghost_values(); 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           phi.submit_dof_value(evaluate_function(function, 
 *                                                  phi.quadrature_point(q)), 
 *                                q); 
 *         inverse.transform_from_q_points_to_basis(dim + 2, 
 *                                                  phi.begin_dof_values(), 
 *                                                  phi.begin_dof_values()); 
 *         phi.set_dof_values(solution); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数再次重复了同样由deal.II库提供的功能，即 VectorTools::integrate_difference(). 我们在这里展示了明确的代码，以强调跨几个单元的矢量化是如何工作的，以及如何通过该接口累积结果。回顾一下，每个<i>lane</i>的矢量化数组持有来自不同单元的数据。通过对当前MPI进程所拥有的所有单元批的循环，我们就可以填充一个结果的VectorizedArray；为了得到一个全局的总和，我们需要进一步去对SIMD阵列中的条目进行求和。然而，这样的程序并不稳定，因为SIMD数组事实上可能并不持有其所有通道的有效数据。当本地拥有的单元的数量不是SIMD宽度的倍数时，就会发生这种情况。为了避免无效数据，我们必须在访问数据时明确地跳过那些无效的通道。虽然人们可以想象，我们可以通过简单地将空车道设置为零（从而不对总和做出贡献）来使其工作，但情况比这更复杂。如果我们要从动量中计算出一个速度呢？那么，我们就需要除以密度，而密度是零--结果就会是NaN，并污染结果。当我们在单元格批次中循环时，使用函数 MatrixFree::n_active_entries_per_cell_batch() 给我们提供有效数据的通道数，累积有效SIMD范围内的结果，就可以避免这种陷阱。它在大多数单元上等于 VectorizedArray::size() ，但如果单元数与SIMD宽度相比有余数，则在最后一个单元批上可能会更少。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   std::array<double, 3> EulerOperator<dim, degree, n_points_1d>::compute_errors( 
 *     const Function<dim> &                             function, 
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     TimerOutput::Scope t(timer, "compute errors"); 
 *     double             errors_squared[3] = {}; 
 *     FEEvaluation<dim, degree, n_points_1d, dim + 2, Number> phi(data, 0, 0); 
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.gather_evaluate(solution, EvaluationFlags::values); 
 *         VectorizedArray<Number> local_errors_squared[3] = {}; 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto error = 
 *               evaluate_function(function, phi.quadrature_point(q)) - 
 *               phi.get_value(q); 
 *             const auto JxW = phi.JxW(q); 
 * 
 *             local_errors_squared[0] += error[0] * error[0] * JxW; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               local_errors_squared[1] += (error[d + 1] * error[d + 1]) * JxW; 
 *             local_errors_squared[2] += (error[dim + 1] * error[dim + 1]) * JxW; 
 *           } 
 *         for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell); 
 *              ++v) 
 *           for (unsigned int d = 0; d < 3; ++d) 
 *             errors_squared[d] += local_errors_squared[d][v]; 
 *       } 
 * 
 *     Utilities::MPI::sum(errors_squared, MPI_COMM_WORLD, errors_squared); 
 * 
 *     std::array<double, 3> errors; 
 *     for (unsigned int d = 0; d < 3; ++d) 
 *       errors[d] = std::sqrt(errors_squared[d]); 
 * 
 *     return errors; 
 *   } 
 * 
 * @endcode
 * 
 * EulerOperator类的最后一个函数是用来估计传输速度的，由网格大小缩放，这与设置显式时间积分器的时间步长有关。在欧拉方程中，有两种传输速度，即对流速度 $\mathbf{u}$ 和相对于以速度 $\mathbf u$ 运动的介质而言，声波的传播速度 $c = \sqrt{\gamma p/\rho}$  。
 * 

 * 
 * 在时间步长的公式中，我们感兴趣的不是这些绝对速度，而是信息穿过一个单元所需的时间量。对于与介质一起传输的信息， $\mathbf u$ 是由网格大小缩放的，所以最大速度的估计可以通过计算 $\|J^{-\mathrm T} \mathbf{u}\|_\infty$  得到，其中 $J$ 是实域到参考域的转换的雅各布。请注意， FEEvaluationBase::inverse_jacobian() 返回的是反转和转置的雅各布，代表从实数到参考坐标的度量项，所以我们不需要再次转置。我们在下面的代码中把这个极限存储在变量`convective_limit`中。
 * 

 * 
 * 声音的传播是各向同性的，所以我们需要考虑到任何方向的网格尺寸。然后，适当的网格大小比例由 $J$ 的最小奇异值给出，或者，等同于 $J^{-1}$ 的最大奇异值。请注意，当忽略弯曲的单元时，可以用单元顶点之间的最小距离来近似这个量。为了得到Jacobian的最大奇异值，一般的策略是使用一些LAPACK函数。由于我们在这里需要的只是一个估计值，所以我们可以避免将一个向量数组的张量分解成几个矩阵的麻烦，并在没有向量的情况下进入一个（昂贵的）特征值函数，而是使用应用于 $J^{-1}J^{-\mathrm T}$ 的幂方法进行几次迭代（在下面的代码中为五次）。这种方法的收敛速度取决于最大特征值与次大特征值的比率以及初始猜测，即所有1的矢量。这可能表明，我们在接近立方体形状的单元上得到缓慢的收敛，在这种情况下，所有的长度几乎都是一样的。然而，这种缓慢的收敛意味着结果将位于两个最大的奇异值之间，而这两个奇异值无论如何都是接近最大值的。在所有其他情况下，收敛将是快速的。因此，我们可以只在这里硬编码5次迭代，并确信结果是好的。
 * 

 * 
 * 
 * @code
 *   template <int dim, int degree, int n_points_1d> 
 *   double EulerOperator<dim, degree, n_points_1d>::compute_cell_transport_speed( 
 *     const LinearAlgebra::distributed::Vector<Number> &solution) const 
 *   { 
 *     TimerOutput::Scope t(timer, "compute transport speed"); 
 *     Number             max_transport = 0; 
 *     FEEvaluation<dim, degree, degree + 1, dim + 2, Number> phi(data, 0, 1); 
 * 
 *     for (unsigned int cell = 0; cell < data.n_cell_batches(); ++cell) 
 *       { 
 *         phi.reinit(cell); 
 *         phi.gather_evaluate(solution, EvaluationFlags::values); 
 *         VectorizedArray<Number> local_max = 0.; 
 *         for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *           { 
 *             const auto solution = phi.get_value(q); 
 *             const auto velocity = euler_velocity<dim>(solution); 
 *             const auto pressure = euler_pressure<dim>(solution); 
 * 
 *             const auto inverse_jacobian = phi.inverse_jacobian(q); 
 *             const auto convective_speed = inverse_jacobian * velocity; 
 *             VectorizedArray<Number> convective_limit = 0.; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               convective_limit = 
 *                 std::max(convective_limit, std::abs(convective_speed[d])); 
 * 
 *             const auto speed_of_sound = 
 *               std::sqrt(gamma * pressure * (1. / solution[0])); 
 * 
 *             Tensor<1, dim, VectorizedArray<Number>> eigenvector; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               eigenvector[d] = 1.; 
 *             for (unsigned int i = 0; i < 5; ++i) 
 *               { 
 *                 eigenvector = transpose(inverse_jacobian) * 
 *                               (inverse_jacobian * eigenvector); 
 *                 VectorizedArray<Number> eigenvector_norm = 0.; 
 *                 for (unsigned int d = 0; d < dim; ++d) 
 *                   eigenvector_norm = 
 *                     std::max(eigenvector_norm, std::abs(eigenvector[d])); 
 *                 eigenvector /= eigenvector_norm; 
 *               } 
 *             const auto jac_times_ev   = inverse_jacobian * eigenvector; 
 *             const auto max_eigenvalue = std::sqrt( 
 *               (jac_times_ev * jac_times_ev) / (eigenvector * eigenvector)); 
 *             local_max = 
 *               std::max(local_max, 
 *                        max_eigenvalue * speed_of_sound + convective_limit); 
 *           } 
 * 
 * @endcode
 * 
 * 与前面的函数类似，我们必须确保只在一个单元格批次的有效单元格上积累速度。
 * 

 * 
 * 
 * @code
 *      for (unsigned int v = 0; v < data.n_active_entries_per_cell_batch(cell);
 *              ++v) 
 *           for (unsigned int d = 0; d < 3; ++d) 
 *             max_transport = std::max(max_transport, local_max[v]); 
 *       } 
 * 
 *     max_transport = Utilities::MPI::max(max_transport, MPI_COMM_WORLD); 
 * 
 *     return max_transport; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TheEulerProblemclass"></a> 
 * <h3>The EulerProblem class</h3>
 * 

 * 
 * 该类将EulerOperator类与时间积分器和通常的全局数据结构（如FiniteElement和DoFHandler）相结合，以实际运行Euler问题的模拟。
 * 

 * 
 * 成员变量是一个三角形、一个有限元、一个映射（用于创建高阶曲面，见 step-10 ），以及一个描述自由度的DoFHandler。此外，我们还保留了上面描述的EulerOperator的实例，它将完成所有积分方面的繁重工作，以及一些时间积分的参数，如当前时间或时间步长。
 * 

 * 
 * 此外，我们使用一个PostProcessor实例来向输出文件写入一些额外的信息，这与  step-33  中的做法类似。DataPostprocessor类的接口很直观，要求我们提供关于需要评估的信息（通常只有解决方案的值，除了Schlieren图，我们只在二维中启用它是有意义的），以及被评估的东西的名称。请注意，也可以通过可视化程序（如ParaView）中的计算器工具来提取大部分信息，但在写输出时就已经做了，这要方便得多。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class EulerProblem 
 *   { 
 *   public: 
 *     EulerProblem(); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid_and_dofs(); 
 * 
 *     void output_results(const unsigned int result_number); 
 * 
 *     LinearAlgebra::distributed::Vector<Number> solution; 
 * 
 *     ConditionalOStream pcout; 
 * 
 * #ifdef DEAL_II_WITH_P4EST 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 * #else 
 *     Triangulation<dim> triangulation; 
 * #endif 
 * 
 *     FESystem<dim>        fe; 
 *     MappingQGeneric<dim> mapping; 
 *     DoFHandler<dim>      dof_handler; 
 * 
 *     TimerOutput timer; 
 * 
 *     EulerOperator<dim, fe_degree, n_q_points_1d> euler_operator; 
 * 
 *     double time, time_step; 
 * 
 *     class Postprocessor : public DataPostprocessor<dim> 
 *     { 
 *     public: 
 *       Postprocessor(); 
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
 *   EulerProblem<dim>::Postprocessor::Postprocessor() 
 *     : do_schlieren_plot(dim == 2) 
 *   {} 
 * 
 * @endcode
 * 
 * 对于字段变量的主要评估，我们首先检查数组的长度是否等于预期值（长度`2*dim+4`或`2*dim+5`来自我们在下面get_names()函数中指定的名字的大小）。然后我们在所有的评估点上循环，填充相应的信息。首先，我们填写密度 $\rho$ 、动量 $\rho \mathbf{u}$ 和能量 $E$ 的原始解变量，然后我们计算得出速度 $\mathbf u$ 、压力 $p$ 、声速 $c=\sqrt{\gamma p / \rho}$ ，以及显示 $s = |\nabla \rho|^2$ 的Schlieren图，如果它被启用。参见 step-69 中另一个创建Schlieren图的例子）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EulerProblem<dim>::Postprocessor::evaluate_vector_field( 
 *     const DataPostprocessorInputs::Vector<dim> &inputs, 
 *     std::vector<Vector<double>> &               computed_quantities) const 
 *   { 
 *     const unsigned int n_evaluation_points = inputs.solution_values.size(); 
 * 
 *     if (do_schlieren_plot == true) 
 *       Assert(inputs.solution_gradients.size() == n_evaluation_points, 
 *              ExcInternalError()); 
 * 
 *     Assert(computed_quantities.size() == n_evaluation_points, 
 *            ExcInternalError()); 
 *     Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError()); 
 *     Assert(computed_quantities[0].size() == 
 *              dim + 2 + (do_schlieren_plot == true ? 1 : 0), 
 *            ExcInternalError()); 
 * 
 *     for (unsigned int q = 0; q < n_evaluation_points; ++q) 
 *       { 
 *         Tensor<1, dim + 2> solution; 
 *         for (unsigned int d = 0; d < dim + 2; ++d) 
 *           solution[d] = inputs.solution_values[q](d); 
 * 
 *         const double         density  = solution[0]; 
 *         const Tensor<1, dim> velocity = euler_velocity<dim>(solution); 
 *         const double         pressure = euler_pressure<dim>(solution); 
 * 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           computed_quantities[q](d) = velocity[d]; 
 *         computed_quantities[q](dim)     = pressure; 
 *         computed_quantities[q](dim + 1) = std::sqrt(gamma * pressure / density); 
 * 
 *         if (do_schlieren_plot == true) 
 *           computed_quantities[q](dim + 2) = 
 *             inputs.solution_gradients[q][0] * inputs.solution_gradients[q][0]; 
 *       } 
 *   } 
 * 
 *   template <int dim> 
 *   std::vector<std::string> EulerProblem<dim>::Postprocessor::get_names() const 
 *   { 
 *     std::vector<std::string> names; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       names.emplace_back("velocity"); 
 *     names.emplace_back("pressure"); 
 *     names.emplace_back("speed_of_sound"); 
 * 
 *     if (do_schlieren_plot == true) 
 *       names.emplace_back("schlieren_plot"); 
 * 
 *     return names; 
 *   } 
 * 
 * @endcode
 * 
 * 对于量的解释，我们有标量密度、能量、压力、声速和Schlieren图，以及动量和速度的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *   EulerProblem<dim>::Postprocessor::get_data_component_interpretation() const 
 *   { 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       interpretation; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       interpretation.push_back( 
 *         DataComponentInterpretation::component_is_part_of_vector); 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 * 
 *     if (do_schlieren_plot == true) 
 *       interpretation.push_back( 
 *         DataComponentInterpretation::component_is_scalar); 
 * 
 *     return interpretation; 
 *   } 
 * 
 * @endcode
 * 
 * 关于必要的更新标志，我们只需要所有数量的值，但Schlieren图除外，它是基于密度梯度的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   UpdateFlags EulerProblem<dim>::Postprocessor::get_needed_update_flags() const 
 *   { 
 *     if (do_schlieren_plot == true) 
 *       return update_values | update_gradients; 
 *     else 
 *       return update_values; 
 *   } 
 * 
 * @endcode
 * 
 * 这个类的构造函数并不令人惊讶。我们设置了一个基于 "MPI_COMM_WORLD "通信器的平行三角形，一个具有 "dim+2 "分量的密度、动量和能量的矢量有限元，一个与底层有限元相同程度的高阶映射，并将时间和时间步长初始化为零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   EulerProblem<dim>::EulerProblem() 
 *     : pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
 * #ifdef DEAL_II_WITH_P4EST 
 *     , triangulation(MPI_COMM_WORLD) 
 * #endif 
 *     , fe(FE_DGQ<dim>(fe_degree), dim + 2) 
 *     , mapping(fe_degree) 
 *     , dof_handler(triangulation) 
 *     , timer(pcout, TimerOutput::never, TimerOutput::wall_times) 
 *     , euler_operator(timer) 
 *     , time(0) 
 *     , time_step(0) 
 *   {} 
 * 
 * @endcode
 * 
 * 作为一个网格，本教程程序实现了两种选择，取决于全局变量`testcase`。对于分析型变量（`testcase==0`），域是 $(0, 10) \times (-5, 5)$ ，域的四周都有迪里希特边界条件（流入）。对于 "testcase==1"，我们将域设置为矩形箱中的圆柱体，源自Sch&auml;fer和Turek（1996）对不可压缩的粘性流动的圆柱体的流动测试案例。在这里，我们有更多种类的边界。通道左侧的流入部分是给定的流入类型，为此我们选择了一个恒定的流入轮廓，而我们在右侧设置了一个亚声速的流出。对于圆柱体周围的边界（边界id等于2）以及通道壁（边界id等于3），我们使用壁的边界类型，即无正态流。此外，对于三维圆柱体，我们还在垂直方向上增加了一个重力。有了基础网格（包括由 GridGenerator::channel_with_cylinder()), 设置的流形），我们就可以执行指定数量的全局细化，从DoFHandler创建未知的编号，并将DoFHandler和Mapping对象交给EulerOperator的初始化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EulerProblem<dim>::make_grid_and_dofs() 
 *   { 
 *     switch (testcase) 
 *       { 
 *         case 0: 
 *           { 
 *             Point<dim> lower_left; 
 *             for (unsigned int d = 1; d < dim; ++d) 
 *               lower_left[d] = -5; 
 * 
 *             Point<dim> upper_right; 
 *             upper_right[0] = 10; 
 *             for (unsigned int d = 1; d < dim; ++d) 
 *               upper_right[d] = 5; 
 * 
 *             GridGenerator::hyper_rectangle(triangulation, 
 *                                            lower_left, 
 *                                            upper_right); 
 *             triangulation.refine_global(2); 
 * 
 *             euler_operator.set_inflow_boundary( 
 *               0, std::make_unique<ExactSolution<dim>>(0)); 
 * 
 *             break; 
 *           } 
 * 
 *         case 1: 
 *           { 
 *             GridGenerator::channel_with_cylinder( 
 *               triangulation, 0.03, 1, 0, true); 
 * 
 *             euler_operator.set_inflow_boundary( 
 *               0, std::make_unique<ExactSolution<dim>>(0)); 
 *             euler_operator.set_subsonic_outflow_boundary( 
 *               1, std::make_unique<ExactSolution<dim>>(0)); 
 * 
 *             euler_operator.set_wall_boundary(2); 
 *             euler_operator.set_wall_boundary(3); 
 * 
 *             if (dim == 3) 
 *               euler_operator.set_body_force( 
 *                 std::make_unique<Functions::ConstantFunction<dim>>( 
 *                   std::vector<double>({0., 0., -0.2}))); 
 * 
 *             break; 
 *           } 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     triangulation.refine_global(n_global_refinements); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     euler_operator.reinit(mapping, dof_handler); 
 *     euler_operator.initialize_vector(solution); 
 * 
 * @endcode
 * 
 * 在下文中，我们输出一些关于问题的统计数据。因为我们经常会出现相当多的单元格或自由度，所以我们希望用逗号来分隔每一组的三位数来打印它们。这可以通过 "locales "来实现，尽管这种工作方式不是特别直观。  step-32 对此有稍微详细的解释。
 * 

 * 
 * 
 * @code
 *     std::locale s = pcout.get_stream().getloc(); 
 *     pcout.get_stream().imbue(std::locale("")); 
 *     pcout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *           << " ( = " << (dim + 2) << " [vars] x " 
 *           << triangulation.n_global_active_cells() << " [cells] x " 
 *           << Utilities::pow(fe_degree + 1, dim) << " [dofs/cell/var] )" 
 *           << std::endl; 
 *     pcout.get_stream().imbue(s); 
 *   } 
 * 
 * @endcode
 * 
 * 对于输出，我们首先让欧拉算子计算出数值结果的误差。更确切地说，对于分析解的情况，我们计算与分析结果的误差，而对于第二个测试情况，我们计算与密度和能量恒定的背景场以及 $x$ 方向的恒定速度的偏差。
 * 

 * 
 * 下一步是创建输出。这与 step-33 中的做法类似：我们让上面定义的后处理器控制大部分的输出，除了我们直接写的原始场。对于分析解的测试案例，我们还对分析解进行了另一次投影，并打印出该场和数值解之间的差异。一旦我们定义了所有要写的量，我们就建立输出的补丁。与 step-65 类似，我们通过设置适当的标志来创建一个高阶VTK输出，这使我们能够可视化高多项式度的场。最后，我们调用 `DataOutInterface::write_vtu_in_parallel()` 函数，将结果写入给定的文件名。这个函数使用了特殊的MPI并行写设施，与其他大多数教程程序中使用的标准库的 `std::ofstream` 变体相比，它通常对并行文件系统更加优化。`write_vtu_in_parallel()`函数的一个特别好的特点是，它可以将所有MPI行列的输出合并到一个文件中，使得没有必要有一个所有此类文件的中央记录（即 "pvtu "文件）。
 * 

 * 
 * 对于并行程序来说，看一下单元在处理器之间的划分往往是有启发的。为此，我们可以向 DataOut::add_data_vector() 传递一个数字向量，其中包含与当前处理器拥有的活动单元一样多的条目；然后这些数字应该是拥有这些单元的处理器的等级。例如，这样一个向量可以从 GridTools::get_subdomain_association(). 中获得。另一方面，在每个MPI进程中，DataOut将只读取那些对应于本地拥有的单元的条目，这些条目当然都有相同的值：即当前进程的等级。矢量的其余条目中的内容实际上并不重要，因此我们可以用一个廉价的技巧逃脱。我们只是把我们给 DataOut::add_data_vector() 的向量的所有*值都填上当前MPI进程的等级。关键是在每个进程中，只有对应于本地拥有的单元格的条目会被读取，而忽略其他条目中的（错误）值。事实上，每个进程提交的向量中的条目子集是正确的，这就足够了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EulerProblem<dim>::output_results(const unsigned int result_number) 
 *   { 
 *     const std::array<double, 3> errors = 
 *       euler_operator.compute_errors(ExactSolution<dim>(time), solution); 
 *     const std::string quantity_name = testcase == 0 ? "error" : "norm"; 
 * 
 *     pcout << "Time:" << std::setw(8) << std::setprecision(3) << time 
 *           << ", dt: " << std::setw(8) << std::setprecision(2) << time_step 
 *           << ", " << quantity_name << " rho: " << std::setprecision(4) 
 *           << std::setw(10) << errors[0] << ", rho * u: " << std::setprecision(4) 
 *           << std::setw(10) << errors[1] << ", energy:" << std::setprecision(4) 
 *           << std::setw(10) << errors[2] << std::endl; 
 * 
 *     { 
 *       TimerOutput::Scope t(timer, "output"); 
 * 
 *       Postprocessor postprocessor; 
 *       DataOut<dim>  data_out; 
 * 
 *       DataOutBase::VtkFlags flags; 
 *       flags.write_higher_order_cells = true; 
 *       data_out.set_flags(flags); 
 * 
 *       data_out.attach_dof_handler(dof_handler); 
 *       { 
 *         std::vector<std::string> names; 
 *         names.emplace_back("density"); 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           names.emplace_back("momentum"); 
 *         names.emplace_back("energy"); 
 * 
 *         std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *           interpretation; 
 *         interpretation.push_back( 
 *           DataComponentInterpretation::component_is_scalar); 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_part_of_vector); 
 *         interpretation.push_back( 
 *           DataComponentInterpretation::component_is_scalar); 
 * 
 *         data_out.add_data_vector(dof_handler, solution, names, interpretation); 
 *       } 
 *       data_out.add_data_vector(solution, postprocessor); 
 * 
 *       LinearAlgebra::distributed::Vector<Number> reference; 
 *       if (testcase == 0 && dim == 2) 
 *         { 
 *           reference.reinit(solution); 
 *           euler_operator.project(ExactSolution<dim>(time), reference); 
 *           reference.sadd(-1., 1, solution); 
 *           std::vector<std::string> names; 
 *           names.emplace_back("error_density"); 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             names.emplace_back("error_momentum"); 
 *           names.emplace_back("error_energy"); 
 * 
 *           std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *             interpretation; 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_scalar); 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             interpretation.push_back( 
 *               DataComponentInterpretation::component_is_part_of_vector); 
 *           interpretation.push_back( 
 *             DataComponentInterpretation::component_is_scalar); 
 * 
 *           data_out.add_data_vector(dof_handler, 
 *                                    reference, 
 *                                    names, 
 *                                    interpretation); 
 *         } 
 * 
 *       Vector<double> mpi_owner(triangulation.n_active_cells()); 
 *       mpi_owner = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD); 
 *       data_out.add_data_vector(mpi_owner, "owner"); 
 * 
 *       data_out.build_patches(mapping, 
 *                              fe.degree, 
 *                              DataOut<dim>::curved_inner_cells); 
 * 
 *       const std::string filename = 
 *         "solution_" + Utilities::int_to_string(result_number, 3) + ".vtu"; 
 *       data_out.write_vtu_in_parallel(filename, MPI_COMM_WORLD); 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * EulerProblem::run() 函数将所有的部分组合起来。它首先调用创建网格和设置数据结构的函数，然后初始化时间积分器和低存储积分器的两个临时向量。我们称这些向量为`rk_register_1`和`rk_register_2`，并使用第一个向量表示 $\mathbf{r}_i$ ，第二个向量表示 $\mathbf{k}_i$ ，在介绍中概述的Runge--Kutta方案的公式。在我们开始时间循环之前，我们通过 `EulerOperator::compute_cell_transport_speed()` 函数计算时间步长。为了便于比较，我们将那里得到的结果与最小网格尺寸进行比较，并将它们打印到屏幕上。对于像本教程程序中接近于统一的声速和速度，预测的有效网格尺寸将是接近的，但如果缩放比例不同，它们可能会有变化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EulerProblem<dim>::run() 
 *   { 
 *     { 
 *       const unsigned int n_vect_number = VectorizedArray<Number>::size(); 
 *       const unsigned int n_vect_bits   = 8 * sizeof(Number) * n_vect_number; 
 * 
 *       pcout << "Running with " 
 *             << Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) 
 *             << " MPI processes" << std::endl; 
 *       pcout << "Vectorization over " << n_vect_number << " " 
 *             << (std::is_same<Number, double>::value ? "doubles" : "floats") 
 *             << " = " << n_vect_bits << " bits (" 
 *             << Utilities::System::get_current_vectorization_level() << ")" 
 *             << std::endl; 
 *     } 
 * 
 *     make_grid_and_dofs(); 
 * 
 *     const LowStorageRungeKuttaIntegrator integrator(lsrk_scheme); 
 * 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_1; 
 *     LinearAlgebra::distributed::Vector<Number> rk_register_2; 
 *     rk_register_1.reinit(solution); 
 *     rk_register_2.reinit(solution); 
 * 
 *     euler_operator.project(ExactSolution<dim>(time), solution); 
 * 
 *     double min_vertex_distance = std::numeric_limits<double>::max(); 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         min_vertex_distance = 
 *           std::min(min_vertex_distance, cell->minimum_vertex_distance()); 
 *     min_vertex_distance = 
 *       Utilities::MPI::min(min_vertex_distance, MPI_COMM_WORLD); 
 * 
 *     time_step = courant_number * integrator.n_stages() / 
 *                 euler_operator.compute_cell_transport_speed(solution); 
 *     pcout << "Time step size: " << time_step 
 *           << ", minimal h: " << min_vertex_distance 
 *           << ", initial transport scaling: " 
 *           << 1. / euler_operator.compute_cell_transport_speed(solution) 
 *           << std::endl 
 *           << std::endl; 
 * 
 *     output_results(0); 
 * 
 * @endcode
 * 
 * 现在我们准备开始时间循环，我们一直运行到时间达到预期的结束时间。每隔5个时间步长，我们就计算一个新的时间步长估计值--由于解决方案是非线性的，在模拟过程中调整这个值是最有效的。如果Courant数选择得过于激进，模拟通常会在时间步数为NaN时爆炸，所以在这里很容易发现。有一点需要注意的是，由于不同的时间步长选择的相互作用，四舍五入的误差可能会传播到前几位数，从而导致略有不同的解决方案。为了降低这种敏感性，通常的做法是将时间步长四舍五入或截断到几位数，例如在这种情况下是3。如果当前时间接近规定的输出 "刻度 "值（如0.02），我们也会写出输出。在时间循环结束后，我们通过打印一些统计数据来总结计算，这主要由 TimerOutput::print_wall_time_statistics() 函数完成。
 * 

 * 
 * 
 * @code
 *     unsigned int timestep_number = 0; 
 * 
 *     while (time < final_time - 1e-12) 
 *       { 
 *         ++timestep_number; 
 *         if (timestep_number % 5 == 0) 
 *           time_step = 
 *             courant_number * integrator.n_stages() / 
 *             Utilities::truncate_to_n_digits( 
 *               euler_operator.compute_cell_transport_speed(solution), 3); 
 * 
 *         { 
 *           TimerOutput::Scope t(timer, "rk time stepping total"); 
 *           integrator.perform_time_step(euler_operator, 
 *                                        time, 
 *                                        time_step, 
 *                                        solution, 
 *                                        rk_register_1, 
 *                                        rk_register_2); 
 *         } 
 * 
 *         time += time_step; 
 * 
 *         if (static_cast<int>(time / output_tick) != 
 *               static_cast<int>((time - time_step) / output_tick) || 
 *             time >= final_time - 1e-12) 
 *           output_results( 
 *             static_cast<unsigned int>(std::round(time / output_tick))); 
 *       } 
 * 
 *     timer.print_wall_time_statistics(MPI_COMM_WORLD); 
 *     pcout << std::endl; 
 *   } 
 * 
 * } // namespace Euler_DG 
 * 
 * @endcode
 * 
 * main()函数并不令人惊讶，它遵循了以前所有MPI程序中的做法。当我们运行一个MPI程序时，我们需要调用`MPI_Init()`和`MPI_Finalize()`，我们通过 Utilities::MPI::MPI_InitFinalize 数据结构来完成。请注意，我们只用MPI来运行程序，并将线程数设置为1。
 * 

 * 
 * 
 * @code
 * int main(int argc, char **argv) 
 * { 
 *   using namespace Euler_DG; 
 *   using namespace dealii; 
 * 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *   try 
 *     { 
 *       deallog.depth_console(0); 
 * 
 *       EulerProblem<dimension> euler_problem; 
 *       euler_problem.run(); 
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
examples/step-67/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


在一台有40个进程的机器上以默认设置运行该程序，会产生以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 147,456 ( = 4 [vars] x 1,024 [cells] x 36 [dofs/cell/var] )
Time step size: 0.00689325, minimal h: 0.3125, initial transport scaling: 0.102759


Time:       0, dt:   0.0069, error rho:   2.76e-07, rho * u:  1.259e-06, energy: 2.987e-06
Time:    1.01, dt:   0.0069, error rho:   1.37e-06, rho * u:  2.252e-06, energy: 4.153e-06
Time:    2.01, dt:   0.0069, error rho:  1.561e-06, rho * u:   2.43e-06, energy: 4.493e-06
Time:    3.01, dt:   0.0069, error rho:  1.714e-06, rho * u:  2.591e-06, energy: 4.762e-06
Time:    4.01, dt:   0.0069, error rho:  1.843e-06, rho * u:  2.625e-06, energy: 4.985e-06
Time:    5.01, dt:   0.0069, error rho:  1.496e-06, rho * u:  1.961e-06, energy: 4.142e-06
Time:       6, dt:   0.0083, error rho:  1.007e-06, rho * u:  7.119e-07, energy: 2.972e-06
Time:       7, dt:   0.0095, error rho:  9.096e-07, rho * u:  3.786e-07, energy: 2.626e-06
Time:       8, dt:   0.0096, error rho:  8.439e-07, rho * u:  3.338e-07, energy:  2.43e-06
Time:       9, dt:   0.0096, error rho:  7.822e-07, rho * u:  2.984e-07, energy: 2.248e-06
Time:      10, dt:   0.0096, error rho:  7.231e-07, rho * u:  2.666e-07, energy: 2.074e-06


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     2.249s    30 |     2.249s |     2.249s     8 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |  0.008066s    13 |   0.00952s |   0.01041s    20 |
| compute transport speed       |       258 |   0.01012s    13 |   0.05392s |   0.08574s    25 |
| output                        |        11 |    0.9597s    13 |    0.9613s |    0.9623s     6 |
| rk time stepping total        |      1283 |    0.9827s    25 |     1.015s |      1.06s    13 |
| rk_stage - integrals L_h      |      6415 |    0.8803s    26 |    0.9198s |    0.9619s    14 |
| rk_stage - inv mass + vec upd |      6415 |   0.05677s    15 |   0.06487s |   0.07597s    13 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



程序输出显示，所有的误差都很小。这是由于我们使用了一个相对较细的 $32^2$ 单元的网格，用5度的多项式来求得一个平滑的解决方案。一个有趣的模式显示在时间步长上：虽然在时间5之前是0.0069，但在后来的时间里增加到0.0096。在时间5和6.5之间，一旦在声速之上有一些运动的旋涡（因此传播速度更快）离开计算域，步长就会增加。在这之后，气流只是在同一方向上是均匀的，与之前均匀速度被漩涡覆盖的状态相比，气体的最大速度有所下降。我们的时间步长公式认识到了这种影响。

最后一块输出显示了关于程序各个部分时间的详细信息；它通过显示最快和最慢的处理器所花费的时间以及平均时间将其分解开来--这在非常大的计算中通常很有用，可以发现是否有处理器持续过热（并因此节制其时钟速度）或因其他原因持续过慢。总结显示，在1.02秒内完成了1283个时间步骤（看所有MPI进程的平均时间），而11个文件的输出又花了0.96秒。将每个时间步数和五个Runge--Kutta阶段分解开来，每次评估的计算时间为0.16毫秒。这种高性能是无矩阵评估器的典型表现，也是显式时间积分对隐式求解器非常有竞争力的原因，特别是对于大规模模拟。程序运行结束时的计算时间细分显示， $\mathcal L_h$ 中的积分评估贡献了大约0.92秒，反质量矩阵的应用贡献了0.06秒。此外，对时间步长计算的运输速度的估计又贡献了0.05秒的计算时间。

如果我们再使用三个级别的全局细化和总共940万个DoF，最终的统计数据如下（对于修改后的Lax--Friedrichs通量， $p=5$  ，和同一系统的40个核心的双插槽Intel Xeon Gold 6230）。

@code
+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     244.9s    12 |     244.9s |     244.9s    34 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        11 |    0.4239s    12 |    0.4318s |    0.4408s     9 |
| compute transport speed       |      2053 |     3.962s    12 |     6.727s |     10.12s     7 |
| output                        |        11 |     30.35s    12 |     30.36s |     30.37s     9 |
| rk time stepping total        |     10258 |     201.7s     7 |     205.1s |     207.8s    12 |
| rk_stage - integrals L_h      |     51290 |     121.3s     6 |     126.6s |     136.3s    16 |
| rk_stage - inv mass + vec upd |     51290 |     66.19s    16 |     77.52s |     81.84s    10 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



每个时间步长，求解器现在需要0.02秒，大约是147k未知数的小问题的25倍。鉴于该问题涉及64倍的未知数，计算时间的增加并不令人惊讶。由于我们也做了8倍的时间步数，计算时间在理论上应该增加512倍。实际增加的时间是205秒/1.02秒=202。这是因为由于通信开销的原因，小问题的规模不能充分利用40个核心。如果我们研究一下每个时间步长所做操作的细节，这一点就很清楚了。带有近邻通信的微分算子 $\mathcal L_h$ 的评估时间从0.92秒到127秒，也就是说，它增加了138倍。另一方面，应用反质量矩阵和向量更新的成本，完全不需要在MPI进程之间通信，增加了1195倍。这一增长超过了理论上的512倍，因为对于较大的尺寸，操作受限于RAM内存的带宽，而对于较小的尺寸，所有的矢量都适合于CPU的缓存。数字显示，尽管使用了低存储量的Runge-Kutta积分器和合并矢量操作，但质量矩阵评估和矢量更新部分几乎消耗了Runge-Kutta阶段所花费的40%的时间。而且尽管对 $\mathcal L_h$ 算子使用了过度积分。对于更简单的微分算子和更昂贵的时间积分器，花费在质量矩阵和矢量更新部分的比例也可以达到70%。如果我们以每秒处理的DoFs和Runge--Kutta阶段计算一个吞吐量数字，我们得到@f[ \text{throughput} =
\frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} = \frac{10258 \cdot 5 \cdot
9.4\,\text{MDoFs}}{205s} = 2360\, \text{MDoFs/s} @f]这个吞吐量数字非常高，因为简单地将一个向量复制到另一个向量的运行速度只有大约10,000 MDoFs/s。

如果我们进入下一个更大的规模，有3770万个DoF，总的模拟时间是2196秒，其中1978秒用于时间步进。L_h算子的运行时间增加了9.3倍（1179秒对127秒），反质量矩阵和向量更新增加了10.3倍（797秒对77.5秒）。运行时间非最佳增长的原因可以追溯到给定硬件上的缓存效应（有40MB的二级缓存和55MB的三级缓存）。虽然不是所有的相关数据都适合940万DoF的缓存（一个向量需要75MB，我们有三个向量加上MatrixFree中的一些额外数据），但还是有能力满足一个半向量的需求。考虑到现代的缓存比天真的最近使用的策略更复杂（在这种情况下，我们几乎没有重复使用，因为数据是以类似流的方式使用的），我们可以假设，在940万DoFs的情况下，确实有相当一部分数据可以从缓存中交付。在更大的情况下，即使有最佳的缓存，也只有不到10%的数据可以放入缓存中，而且会有相关的性能损失。




<a name="Convergenceratesfortheanalyticaltestcase"></a><h3>Convergence rates for the analytical test case</h3>


对于修改后的Lax--Friedrichs通量和测量动量变量的误差，我们得到以下收敛表（密度和能量变量的速率非常相似）。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.373e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.130e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.899e-03</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">5.577e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">7.381e-03</td>
    <td>3.64</td>
    <td align="right">36,864</td>
    <td align="center">2.082e-04</td>
    <td>5.42</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.724e-03</td>
    <td>3.56</td>
    <td align="right">65,536</td>
    <td align="center">3.072e-04</td>
    <td>4.59</td>
    <td align="right">147,456</td>
    <td align="center">2.625e-06</td>
    <td>6.31</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">6.205e-04</td>
    <td>2.92</td>
    <td align="right">262,144</td>
    <td align="center">1.880e-05</td>
    <td>4.03</td>
    <td align="right">589,824</td>
    <td align="center">3.268e-08</td>
    <td>6.33</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">8.279e-05</td>
    <td>2.91</td>
    <td align="right">1,048,576</td>
    <td align="center">1.224e-06</td>
    <td>3.94</td>
    <td align="right">2,359,296</td>
    <td align="center">9.252e-10</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.105e-05</td>
    <td>2.91</td>
    <td align="right">4,194,304</td>
    <td align="center">7.871e-08</td>
    <td>3.96</td>
    <td align="right">9,437,184</td>
    <td align="center">1.369e-10</td>
    <td>2.77</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">1.615e-06</td>
    <td>2.77</td>
    <td align="right">16,777,216</td>
    <td align="center">4.961e-09</td>
    <td>3.99</td>
    <td align="right">37,748,736</td>
    <td align="center">7.091e-11</td>
    <td>0.95</td>
  </tr>
</table> 

如果我们改用Harten-Lax-van Leer通量，结果如下。   <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th colspan="3"><i>p</i>=2</th>
    <th colspan="3"><i>p</i>=3</th>
    <th colspan="3"><i>p</i>=5</th>
  </tr>
  <tr>
    <th>n_cells</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
    <th>n_dofs</th>
    <th>error mom</th>
    <th>rate</th>
  </tr>
  <tr>
    <td align="right">16</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">2,304</td>
    <td align="center">1.339e-01</td>
    <td>&nbsp;</td>
  </tr>
  <tr>
    <td align="right">64</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td>&nbsp;</td>
    <td align="right">4,096</td>
    <td align="center">9.037e-02</td>
    <td>&nbsp;</td>
    <td align="right">9,216</td>
    <td align="center">8.849e-03</td>
    <td>3.92</td>
  </tr>
  <tr>
    <td align="right">256</td>
    <td align="right">9,216</td>
    <td align="center">4.204e-02</td>
    <td>&nbsp;</td>
    <td align="right">16,384</td>
    <td align="center">9.143e-03</td>
    <td>3.31</td>
    <td align="right">36,864</td>
    <td align="center">2.501e-04</td>
    <td>5.14</td>
  </tr>
  <tr>
    <td align="right">1024</td>
    <td align="right">36,864</td>
    <td align="center">4.913e-03</td>
    <td>3.09</td>
    <td align="right">65,536</td>
    <td align="center">3.257e-04</td>
    <td>4.81</td>
    <td align="right">147,456</td>
    <td align="center">3.260e-06</td>
    <td>6.26</td>
  </tr>
  <tr>
    <td align="right">4096</td>
    <td align="right">147,456</td>
    <td align="center">7.862e-04</td>
    <td>2.64</td>
    <td align="right">262,144</td>
    <td align="center">1.588e-05</td>
    <td>4.36</td>
    <td align="right">589,824</td>
    <td align="center">2.953e-08</td>
    <td>6.79</td>
  </tr>
  <tr>
    <td align="right">16,384</td>
    <td align="right">589,824</td>
    <td align="center">1.137e-04</td>
    <td>2.79</td>
    <td align="right">1,048,576</td>
    <td align="center">9.400e-07</td>
    <td>4.08</td>
    <td align="right">2,359,296</td>
    <td align="center">4.286e-10</td>
    <td>6.11</td>
  </tr>
  <tr>
    <td align="right">65,536</td>
    <td align="right">2,359,296</td>
    <td align="center">1.476e-05</td>
    <td>2.95</td>
    <td align="right">4,194,304</td>
    <td align="center">5.799e-08</td>
    <td>4.02</td>
    <td align="right">9,437,184</td>
    <td align="center">2.789e-11</td>
    <td>3.94</td>
  </tr>
  <tr>
    <td align="right">262,144</td>
    <td align="right">9,437,184</td>
    <td align="center">2.038e-06</td>
    <td>2.86</td>
    <td align="right">16,777,216</td>
    <td align="center">3.609e-09</td>
    <td>4.01</td>
    <td align="right">37,748,736</td>
    <td align="center">5.730e-11</td>
    <td>-1.04</td>
  </tr>
</table> 

表中显示，我们对两种数值通量都得到了最佳的 $\mathcal O\left(h^{p+1}\right)$ 收敛率。对于 $p=2$ 的Lax--Friedrichs通量，误差略小，但对于 $p=3$ 的情况则相反；在任何情况下，这个测试案例的差异都相对较小。

对于 $p=5$ ，我们在最细的网格上用两种通量达到了 $10^{-11}$ 的舍入精度。还要注意的是，误差是绝对的，域长为 $10^2$ ，所以相对误差低于 $10^{-12}$ 。HLL通量对于最高度数来说要好一些，这是由于Lax--Friedrichs通量的轻微不准确造成的。Lax--Friedrichs通量对离开域的解设置了一个Dirichlet条件，这导致了一个小的人工反射，这在Lax--Friedrichs通量中被凸显出来。除此之外，我们看到数值通量的影响很小，因为元素内部的多项式部分是引起反射的主要动力。当试图用高阶DG设置来接近更具挑战性的设置时，通量的有限影响也会产生影响。以第33步的参数和网格为例，一旦高质部分接近边界，我们就会在两种通量下得到振荡（这反过来会使密度为负值，并使解决方案爆炸），这与低阶有限体积情况不同（ $p=0$ ）。因此，任何导致溶液中出现冲击的情况都需要某种形式的限制性或人工耗散。对于另一种选择，请参见step-69教程程序。




<a name="Resultsforflowinchannelaroundcylinderin2D"></a><h3>Results for flow in channel around cylinder in 2D</h3>


对于渠道中圆柱体周围的流动测试案例，我们需要将第一行代码改为

@code
  constexpr unsigned int testcase = 1;
@endcode

这个测试案例从一个马赫数为0.31的恒定速度和恒定的初始密度的背景场开始；气流必须绕过一个圆柱体形式的障碍物。由于我们对圆柱体壁施加了一个无穿透的条件，最初迎面撞上圆柱体的气流必须重新排列，这就产生了一个大的声波。下面的图片显示了二维情况下5级全局细化时0.1、0.25、0.5和1.0（左上至右下）的压力，使用了102,400个单元，多项式程度为5，所有4个求解变量的自由度为1470万。我们清楚地看到，在时间0.1的第一个快照中，不连续现象在上游方向传播缓慢，在下游方向传播较快。在时间0.25，声波已经到达顶部和底部的墙壁并反射到内部。从下壁和上壁反射波的不同距离，我们可以看到以 GridGenerator::channel_with_cylinder() 为代表的Sch&auml;fer-Turek试验案例的轻微不对称性，圆柱体上方的空间与下方相比要多一些。在后来的时间里，画面更加混乱，到处都是许多声波。

 <table align="center" class="doxtable" style="width:85%">
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_010.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_025.png" alt="" width="100%">
    </td>
  </tr>
  <tr>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_050.png" alt="" width="100%">
    </td>
    <td>
        <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_100.png" alt="" width="100%">
    </td>
  </tr>
</table> 

下一张图片显示了在相同分辨率下，从通道入口向出口看，时间为1.0时的压力仰角图--在这里，我们可以看到大量的反射。在该图中，可以看到两种类型的波。较大振幅的波对应于初始不连续物撞击墙壁时发生的各种反射，而与元素大小相似的小振幅波则对应于数值伪影。它们起源于方案的有限分辨率，并在不连续面通过高阶多项式的元素时出现。这种效应可以通过提高分辨率来治愈。除了这种效应之外，丰富的波浪结构是高阶DG方法的传输精度的结果。

 <img src="https://www.dealii.org/images/steps/developer/step-67.pressure_elevated.jpg" alt="" width="40%"> 

通过2级全局细化，1,600个单元，网格及其在40个MPI进程上的划分情况如下。

 <img src="https://www.dealii.org/images/steps/developer/step-67.grid-owner.png" alt="" width="70%"> 

当我们在40个核心上运行具有4级全局细化的代码时，我们得到以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 3,686,400 ( = 4 [vars] x 25,600 [cells] x 36 [dofs/cell/var] )
Time step size: 7.39876e-05, minimal h: 0.001875, initial transport scaling: 0.00110294


Time:       0, dt:  7.4e-05, norm rho:   4.17e-16, rho * u:  1.629e-16, energy: 1.381e-15
Time:    0.05, dt:  6.3e-05, norm rho:    0.02075, rho * u:    0.03801, energy:   0.08772
Time:     0.1, dt:  5.9e-05, norm rho:    0.02211, rho * u:    0.04515, energy:   0.08953
Time:    0.15, dt:  5.7e-05, norm rho:    0.02261, rho * u:    0.04592, energy:   0.08967
Time:     0.2, dt:  5.8e-05, norm rho:    0.02058, rho * u:    0.04361, energy:   0.08222
Time:    0.25, dt:  5.9e-05, norm rho:    0.01695, rho * u:    0.04203, energy:   0.06873
Time:     0.3, dt:  5.9e-05, norm rho:    0.01653, rho * u:     0.0401, energy:   0.06604
Time:    0.35, dt:  5.7e-05, norm rho:    0.01774, rho * u:    0.04264, energy:    0.0706


...


Time:    1.95, dt:  5.8e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05185
Time:       2, dt:  5.7e-05, norm rho:    0.01432, rho * u:    0.03969, energy:   0.04889


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |     273.6s    13 |     273.6s |     273.6s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.01112s    35 |    0.0672s |    0.1337s     0 |
| compute transport speed       |      6914 |     5.422s    35 |     15.96s |     29.99s     1 |
| output                        |        41 |     37.24s    35 |      37.3s |     37.37s     0 |
| rk time stepping total        |     34564 |     205.4s     1 |     219.5s |     230.1s    35 |
| rk_stage - integrals L_h      |    172820 |     153.6s     1 |     164.9s |     175.6s    27 |
| rk_stage - inv mass + vec upd |    172820 |     47.13s    13 |     53.09s |     64.05s    33 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



这里显示的各种数量的规范是对背景场（即初始条件）的偏差 $\rho'$ 、 $(\rho u)'$ 和 $E'$ 。运行时间的分布总体上与之前的测试案例相似。唯一略有不同的是，与反质量矩阵和矢量更新相比，在 $\mathcal L_h$ 中花费的时间比例较大。这是因为几何体是变形的，无矩阵框架需要从内存中加载额外的几何体阵列，这些阵列在仿生网格的情况下是被压缩的。

将全局细化的数量增加到5，输出就变成了。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 14,745,600 ( = 4 [vars] x 102,400 [cells] x 36 [dofs/cell/var] )


...


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              |      2693s    32 |      2693s |      2693s    23 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |   0.04537s    32 |     0.173s |    0.3489s     0 |
| compute transport speed       |     13858 |     40.75s    32 |     85.99s |     149.8s     0 |
| output                        |        41 |     153.8s    32 |     153.9s |     154.1s     0 |
| rk time stepping total        |     69284 |      2386s     0 |      2450s |      2496s    32 |
| rk_stage - integrals L_h      |    346420 |      1365s    32 |      1574s |      1718s    19 |
| rk_stage - inv mass + vec upd |    346420 |     722.5s    10 |     870.7s |      1125s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



对性能的影响与分析性测试案例相似--理论上，计算时间应该增加8倍，但我们实际上看到时间步骤增加了11倍（219.5秒对2450秒）。这可以追溯到缓存，小的情况下大多适合缓存。一个有趣的效果，是典型的本地通信（积分 $\mathcal L_h$ ）和全局通信（计算运输速度）混合的程序，有一些负载不平衡，可以通过查看分别遇到不同阶段的最小和最大时间的MPI等级来观察。级别0报告了 "rk时间步进总数 "部分的最快吞吐量。同时，对于 "计算传输速度 "部分，它似乎是最慢的，几乎比平均水平慢了2倍，与较快的等级相比几乎是4倍。由于后者涉及到全局通信，我们可以将这部分的缓慢归因于本地Runge--Kutta阶段在这个等级上推进得更快，需要等到其他处理器跟上。在这一点上，人们可以怀疑这种不平衡的原因。在所有的MPI进程中，单元格的数量几乎是相同的。然而，无矩阵框架在位于通道出口处的仿生和笛卡尔单元上速度更快，较低的MPI等级被分配到这些单元。另一方面，报告Runga--Kutta阶段最高运行时间的等级32拥有靠近圆柱体的弯曲单元，对于这些单元不可能有数据压缩。为了提高吞吐量，我们可以在划分 parallel::distributed::Triangulation 对象时给不同的单元类型分配不同的权重，甚至可以测量几个时间步骤的运行时间，然后尝试重新平衡。

对于1470万DoFs的测试案例，在346000个Runge--Kutta阶段中，每个Runge--Kutta阶段的吞吐量可以计算到2085 MDoFs/s，比上面报告的2360 MDoFs/s的笛卡尔网格吞吐量略慢。

最后，如果我们增加一个额外的细化，我们会记录以下输出。

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 58,982,400 ( = 4 [vars] x 409,600 [cells] x 36 [dofs/cell/var] )


...


Time:    1.95, dt:  1.4e-05, norm rho:    0.01488, rho * u:    0.03923, energy:   0.05183
Time:       2, dt:  1.4e-05, norm rho:    0.01431, rho * u:    0.03969, energy:   0.04887


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 2.166e+04s    26 | 2.166e+04s | 2.166e+04s    24 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.1758s    30 |     0.672s |     1.376s     1 |
| compute transport speed       |     27748 |     321.3s    34 |     678.8s |      1202s     1 |
| output                        |        41 |     616.3s    32 |     616.4s |     616.4s    34 |
| rk time stepping total        |    138733 | 1.983e+04s     1 | 2.036e+04s | 2.072e+04s    34 |
| rk_stage - integrals L_h      |    693665 | 1.052e+04s    32 | 1.248e+04s | 1.387e+04s    19 |
| rk_stage - inv mass + vec upd |    693665 |      6404s    10 |      7868s | 1.018e+04s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



rk时间步数总数 "部分对应的是2010 MDoFs/s的吞吐量。执行139k时间步长的总体运行时间是20k秒（5.7小时）或每秒7个时间步长--对于有近6000万个未知数来说还不错。通过在计算中添加更多的内核，可以实现更多的吞吐量。




<a name="Resultsforflowinchannelaroundcylinderin3D"></a><h3>Results for flow in channel around cylinder in 3D</h3>


将通道测试案例切换到3D，并进行3次全局细化，输出结果是

@code
Running with 40 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 221,184,000 ( = 5 [vars] x 204,800 [cells] x 216 [dofs/cell/var] )


...


Time:    1.95, dt:  0.00011, norm rho:    0.01131, rho * u:    0.03056, energy:   0.04091
Time:       2, dt:  0.00011, norm rho:     0.0119, rho * u:    0.03142, energy:   0.04425


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 1.734e+04s     4 | 1.734e+04s | 1.734e+04s    38 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |    0.6551s    34 |     3.216s |     7.281s     0 |
| compute transport speed       |      3546 |       160s    34 |     393.2s |     776.9s     0 |
| output                        |        41 |      1350s    34 |      1353s |      1357s     0 |
| rk time stepping total        |     17723 | 1.519e+04s     0 | 1.558e+04s | 1.582e+04s    34 |
| rk_stage - integrals L_h      |     88615 | 1.005e+04s    32 | 1.126e+04s |  1.23e+04s    11 |
| rk_stage - inv mass + vec upd |     88615 |      3056s    11 |      4322s |      5759s    32 |
+-------------------------------------------+------------------+------------+------------------+
@endcode



物理原理与二维情况类似，由于引力的作用，在Z方向有轻微的运动。在这种情况下，每个Runge-Kutta阶段的吞吐量为

@f[
\text{throughput} = \frac{n_\mathrm{time steps} n_\mathrm{stages}
n_\mathrm{dofs}}{t_\mathrm{compute}} =
\frac{17723 \cdot 5 \cdot 221.2\,\text{M}}{15580s} = 1258\, \text{MDoFs/s}.


@f]



吞吐量低于二维，因为 $\mathcal L_h$ 项的计算更加昂贵。这是由于 "度+2 "点的过度积分和较大比例的面积分（更差的体积-表面比率），以及更昂贵的通量计算。如果我们只考虑反质量矩阵和矢量更新部分，我们记录到等熵涡旋的二维案例的吞吐量为4857 MDoFs/s，有3770万个未知数，而三维案例的运行速度为4535 MDoFs/s。性能是相似的，因为这两种情况实际上都受到内存带宽的限制。

如果我们进行四级全局细化，我们需要增加进程的数量以在内存中容纳所有的东西--在这种情况下，计算需要大约350GB的RAM内存。另外，通过增加额外的资源，完成35k个时间步骤所需的时间也变得更容易忍受。因此，我们使用了6个节点，每个节点有40个核心，从而使计算有240个MPI进程。

@code
Running with 240 MPI processes
Vectorization over 8 doubles = 512 bits (AVX512)
Number of degrees of freedom: 1,769,472,000 ( = 5 [vars] x 1,638,400 [cells] x 216 [dofs/cell/var] )


...


Time:    1.95, dt:  5.6e-05, norm rho:    0.01129, rho * u:     0.0306, energy:   0.04086
Time:       2, dt:  5.6e-05, norm rho:    0.01189, rho * u:    0.03145, energy:   0.04417


+-------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed              | 5.396e+04s   151 | 5.396e+04s | 5.396e+04s     0 |
|                                           |                  |                               |
| Section                       | no. calls |   min time  rank |   avg time |   max time  rank |
+-------------------------------------------+------------------+------------+------------------+
| compute errors                |        41 |     2.632s   178 |     7.221s |     16.56s     0 |
| compute transport speed       |      7072 |       714s   193 |      1553s |      3351s     0 |
| output                        |        41 |      8065s   176 |      8070s |      8079s     0 |
| rk time stepping total        |     35350 |  4.25e+04s     0 |  4.43e+04s | 4.515e+04s   193 |
| rk_stage - integrals L_h      |    176750 | 2.936e+04s   134 | 3.222e+04s |  3.67e+04s    99 |
| rk_stage - inv mass + vec upd |    176750 |      7004s    99 | 1.207e+04s |  1.55e+04s   132 |
+-------------------------------------------+------------------+------------+------------------+
@endcode

这个模拟有近20亿个未知数--确实是一个相当大的计算量，而每个时间步长仍然只需要大约1.5秒。




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这里介绍的代码可以直接扩展到自适应网格，给定适当的指标来设置细化标志。在声波方程的背景下，类似求解器的大规模适应性已经由<a href="https://github.com/kronbichler/exwave">exwave
project</a>实现。然而，在目前的情况下，自适应性的好处往往只限于靠近声波起源的早期时间和效果，因为波最终会反射和衍射。这就导致了到处都是陡峭的梯度，类似于湍流，以及或多或少的全局细化网格。

我们在结果部分没有讨论的另一个话题是不同时间积分方案的比较。该程序提供了四种低存储量的Runga--Kutta积分器的变体，每一种都有轻微不同的精度和稳定性行为。在这里实现的方案中，高阶方案提供了额外的精度，但在违反CFL条件之前，每级步长的效率略低。一个有趣的扩展是将这里提出的低存储变体与标准的Runge--Kutta积分器进行比较，或者使用与质量矩阵运算分开运行的矢量运算，并比较性能。




<a name="Moreadvancednumericalfluxfunctionsandskewsymmetricformulations"></a><h4>More advanced numerical flux functions and skew-symmetric formulations</h4>


正如介绍中提到的，本程序中采用的修改的Lax--Friedrichs通量和HLL通量只是文献中关于欧拉方程的大量数值通量中的两个变种。一个例子是HLLC通量（Harten-Lax-van Leer-Contact）通量，它增加了HLL通量或Roe通量中缺少的稀疏波效应。正如介绍中提到的，数值通量对高阶DG方案的影响是有争议的（与低阶离散的情况不同）。

为了提高求解器的稳定性，一个相关的改进是也要考虑空间积分项。上面使用的相当幼稚的实现方式的一个缺点是，原始欧拉方程的能量守恒（在没有冲击的情况下）只适用于离散化误差。如果解决方案的分辨率不足，离散化误差会引起数值能量的增加，并最终导致离散化的不稳定。这是因为欧拉方程中的项的不精确数值积分，其中包含有理非线性和弯曲单元的高阶内容。摆脱这种困境的方法是所谓的倾斜对称公式，见 @cite Gassner2013 的一个简单变体。倾斜对称意味着在弱式中切换解 $\mathbf{w}$ 和检验函数 $\mathbf{v}$ 的作用，除了一些边界项外，产生原始量的精确负值。在离散设置中，挑战在于当积分只被近似计算时也要保持这种倾斜对称性（在连续情况下，倾斜对称性是部分积分的结果）。偏斜对称的数值方案平衡了保守形式的空间导数  $(\nabla \mathbf v, \mathbf{F}(\mathbf w))_{K}$  和对流形式的贡献  $(\mathbf v, \tilde{\mathbf{F}}(\mathbf w)\nabla
\mathbf{w})_{K}$  ，对于某些  $\tilde{\mathbf{F}}$  。准确的条款取决于方程和积分公式，在某些情况下可以通过特殊的倾斜对称有限差分方案来理解。

要想开始，有兴趣的读者可以看看https://github.com/kronbichler/advection_miniapp，其中用deal.II对一个简单的平流方程实现了倾斜对称的DG公式。

<a name="Equippingthecodeforsupersoniccalculations"></a><h4>Equipping the code for supersonic calculations</h4>


正如介绍中提到的，欧拉方程的解随着马赫数的增加而产生冲击，这需要额外的机制来稳定方案，例如限制器的形式。除了实际实施限制器或人工粘性方法外，主要的挑战是如何平衡计算，因为在有问题的单元中限制震荡所涉及的额外计算会使它们比没有限制的普通DG单元更昂贵。此外，更好地应对不连续情况的额外数值通量也是一种选择。

对于超音速流动来说，有一个因素也是必要的，那就是适当的边界条件。与介绍中讨论的并在程序中实现的亚音速流出边界相反，所有的特性都是超音速流出边界的外在表现，所以我们不想规定任何外部数据。

@f[
\mathbf{w}^+ = \mathbf{w}^- = \begin{pmatrix} \rho^-\\
(\rho \mathbf u)^- \\ E^-\end{pmatrix} \quad
 \text{(Neumann)}.


@f]



在代码中，我们将简单地添加额外的语句

@code
            else if (supersonic_outflow_boundaries.find(boundary_id) !=
                     supersonic_outflow_boundaries.end())
              {
                w_p        = w_m;
                at_outflow = true;
              }
@endcode

在 "local_apply_boundary_face() "函数中。

<a name="ExtensiontothelinearizedEulerequations"></a><h4>Extension to the linearized Euler equations</h4>


当对欧拉解的兴趣主要在于声波的传播时，围绕一个背景状态，即一个给定的密度、速度和能量（或压力）场，将欧拉方程线性化，只计算针对这些场的变化，往往是合理的。这就是航空声学的广泛领域的设置。即使有时分辨率要求大大降低，但由于线性化引起了额外的条款，实施起来就变得有些复杂了。从代码的角度来看，在算子评估中，我们还需要为代码配备要线性化的状态。这一信息可以由分析函数（根据正交点的位置进行评估）或由类似于解决方案的矢量提供。基于该矢量，我们将创建一个额外的FEEvaluation对象，从中读取并提供正交点的场值。如果背景速度为零，密度为常数，线性化的欧拉方程进一步简化，可以等效地写成声波方程的形式。

在声音传播的背景下，一个挑战往往是边界条件的定义，因为计算域需要是有限的，而实际模拟往往跨越无限的（或至少大得多）物理域。传统的Dirichlet或Neumann边界条件会引起声波的反射，最终传播到感兴趣的区域，破坏了解决方案。因此，各种非反射边界条件或海绵层的变体，通常以<a
href="https://en.wikipedia.org/wiki/Perfectly_matched_layer">perfectly
matched layers</a>的形式出现--其中解决方案被阻尼，没有反射

--是很常见的。




<a name="ExtensiontothecompressibleNavierStokesequations"></a><h4>Extension to the compressible Navier-Stokes equations</h4>


如 @cite FehnWallKronbichler2019 所述，本教程程序中的求解器也可以通过添加粘性项扩展到可压缩的Navier-Stokes方程。为了尽量保持这里获得的性能，尽管有额外的椭圆项的成本，例如通过内部惩罚方法，我们可以像步骤59的教程程序一样，将基础从FE_DGQ切换到FE_DGQHermite。




<a name="Usingcellcentricloopsandsharedmemory"></a><h4>Using cell-centric loops and shared memory</h4>


在本教程中，我们使用了以面为中心的循环。在这里，单元和面的积分在不同的循环中处理，导致对结果向量的多次写入访问，这在现代硬件上是比较昂贵的，因为写入操作通常也会导致隐含的读操作。另一方面，以元素为中心的循环是在处理一个单元的同时直接处理其所有的2d面。虽然这种循环意味着通量必须计算两次（对于一个内部面的每一面），但结果向量只需访问一次的事实--以及由此产生的算法没有竞赛条件，因此完全适合共享内存的事实--已经带来了性能的提升。如果你对这些高级主题感兴趣，你可以看一下步骤76，在那里我们对本教程进行了修改，以便我们能够使用这些功能。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-67.cc"
*/
