/**
@page step_70 The step-70 tutorial program
This tutorial depends on step-19, step-32, step-60.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Massivelyparallelnonmatchinggridsimulationsoffluidstructureinteractionproblems">Massively parallel non-matching grid simulations of fluid structure interaction problems</a>
      <ul>
        <li><a href="#Codimensiononecase">Co-dimension one case</a>
        <li><a href="#Codimensionzerocase">Co-dimension zero case</a>
        <li><a href="#Representationofand">Representation of Ω and Γ</a>
        <li><a href="#Usingparticlestotrack">Using particles to track Γ</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
        <li><a href="#Morereferences"> More references</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Runtimeparameterhandling">Run-time parameter handling</a>
        <li><a href="#TheStokesImmersedProblemclassdeclaration">The StokesImmersedProblem class declaration</a>
        <li><a href="#TheStokesImmersedProblemclassimplementation">The StokesImmersedProblem class implementation</a>
      <ul>
        <li><a href="#Objectconstructionandmeshinitializationfunctions">Object construction and mesh initialization functions</a>
        <li><a href="#Particleinitializationfunctions">Particle initialization functions</a>
        <li><a href="#DoFinitializationfunctions">DoF initialization functions</a>
        <li><a href="#Assemblyfunctions">Assembly functions</a>
        <li><a href="#Solvingthelinearsystem">Solving the linear system</a>
        <li><a href="#Meshrefinement">Mesh refinement</a>
        <li><a href="#Creatingoutputforvisualization">Creating output for visualization</a>
        <li><a href="#Therunfunction">The "run" function</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Twodimensionaltestcase"> Two dimensional test case </a>
        <li><a href="#Threedimensionaltestcase"> Three dimensional test case </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-70/doc/intro.dox

 <br> 

<i>This program was contributed by Luca Heltai (International School for
Advanced Studies, Trieste), Bruno Blais (Polytechnique Montréal),
and Rene Gassmöller (University of California Davis)
</i>

 @dealiiTutorialDOI{10.5281/zenodo.3829064,https://zenodo.org/badge/DOI/10.5281/zenodo.3829064.svg} 




<a name="Introduction"></a><h1>Introduction</h1>


<a name="Massivelyparallelnonmatchinggridsimulationsoffluidstructureinteractionproblems"></a><h3>Massively parallel non-matching grid simulations of fluid structure interaction problems</h3>


在本教程中，我们考虑了层流体系中的混合问题。这类问题出现在从化学工程到发电（如涡轮机械）等广泛的应用中。混合问题特别难以用数值来解决，因为它们通常涉及一个容器（有固定的边界，可能还有复杂的几何形状，如挡板），由域 $\Omega$ 表示，和一个（或多个）浸入和旋转的叶轮（由域 $\Omega^{\text{imp}}$ 表示）。我们希望解决流动方程的域是两个域之间的（与时间有关的）差值，即。   $\Omega\setminus\Omega^{\text{imp}}$  .

对于旋转叶轮，使用任意拉格朗日欧拉公式（其中流体域--连同网格！）是不可能的，除非只考虑小时间（即小的流体域变形）。-- 是不可能的，除非只考虑小时间（即小的流域变形）。如果想跟踪叶轮多次旋转时的流动演变，所产生的变形网格就会过于扭曲而无用。

在这种情况下，一个可行的替代策略是使用非匹配方法（类似于我们在step-60中所做的），其中一个背景固定网格（可能在时间上进行局部细化以更好地捕捉实体运动）与一个旋转的、独立的网格相耦合。

为了保持步骤60中使用的相同符号，我们使用 $\Omega$ 来表示 ${\mathbb R}^{\text{spacedim}}$ 中的域，代表流体和叶轮的容器，我们使用 $\Gamma$ 在 ${\mathbb R}^{\text{dim}}$ 来表示整个叶轮（当它的`spacedim`度量非负值时，也就是说，当我们可以把它表示为维数`dim`等于`spacedim`的网格时），薄叶轮的同维度表示，或者只是整个叶轮的边界。

域 $\Gamma$ 被嵌入到 $\Omega$ （ $\Gamma \subseteq \Omega$ ）中，它是不匹配的：一般来说，它不与任何体积网格的特征对齐。我们在 $\Omega$ 上求解一个偏微分方程，通过一些惩罚技术在嵌入域 $\Gamma$ 上强制执行一些问题的解决条件。在当前情况下，条件是流体在 $\Gamma$ 上各点的速度等于固体叶轮在该点的速度。

我们在此描述的技术在文献中使用了许多名称之一：<b>immersed finite element method</b>和<b>fictitious boundary
method</b>等。  其主要原理是两个网格的离散化保持完全独立。在本教程中，这种方法被用来求解由斯托克斯方程描述的粘性流体的运动，该流体被一个刚性的非变形叶轮搅动。

因此， $\Omega$ 中求解的方程是蠕动流的斯托克斯方程（即 $\text{Re}\rightarrow 0$ ），并且在与叶轮相关的移动*嵌入域* $\Gamma$ 上应用无滑动边界条件。然而，这个教程可以很容易地扩展到其他方程（例如，纳维-斯托克斯方程、线性弹性方程等）。它可以被看作是Step-60的一个自然扩展，它可以通过MPI使用分布式并行计算架构解决大型问题。

然而，与第60步相反， $\Gamma$ 上的迪里希特边界条件是弱加的，而不是通过使用拉格朗日乘法器，而且我们集中处理两个完全分布的三角形的耦合（这种组合在第60步的实施中是不可能的）。

当人们想在嵌入域上执行条件时，有两种有趣的情况发生  $\Gamma$  。

- 嵌入域 $\Gamma$ 的几何维度`dim`与域 $\Omega$ 相同（`spacedim`），也就是说， $\Gamma$ 的spacedim-维度不为零。在这种情况下，对 $\Gamma$ 施加Dirichlet边界的边界条件是通过体积惩罚完成的。如果施加的惩罚只取决于速度，这通常被称为 $\mathcal{L}^2$ 惩罚，而如果惩罚同时取决于速度及其梯度，则是 $\mathcal{H}^1$ 惩罚。 $\mathcal{L}^2$  惩罚的情况与Darcy型方法非常相似。对 $\mathcal{L}^2$ 和 $\mathcal{H}^1$ 两种惩罚方法都进行了广泛的分析（例如，见 @cite Angot1999 ）。

- 嵌入域 $\Gamma$ 的内在维度`dim`小于 $\Omega$ 的维度（`spacedim`），因此其spacedim维度为零；例如，它是一条嵌入二维域的曲线，或一个嵌入三维域的表面。当然，这在物理上是不可能的，但是如果金属片的厚度可以忽略不计的话，我们可以把在流体中运动的非常薄的金属片视为本质上的低维。在这种情况下，通过应用<a href="https://en.wikipedia.org/wiki/Joachim_Nitsche">Nitsche</a>方法（见 @cite Freund1995 ）对 $\Gamma$ 施加弱边界条件。

这两种方法都有非常相似的要求，并导致高度相似的公式。因此，我们几乎以同样的方式对待它们。

在本教程中，我们对 $\Gamma$ 的进一步细节不感兴趣：我们假设嵌入域的尺寸（`dim`）总是比嵌入域的尺寸 $\Omega$ （`spacedim`）小一或相等。

我们要解决以下微分问题：给定 $g$ 上的一个足够规则的函数 $\Gamma$ ，找到 $(\textbf{u},p)$ 的解。

@f{eqnarray*}


  -\Delta \mathbf{u} + \nabla p &=& 0,\\


  -\nabla \cdot \textbf{u} &=& 0,\\
  \textbf{u} &=& \textbf{g}  \text{ in } \Gamma,\\
  \textbf{u} &=& 0 \text{ on } \partial\Omega.


@f}



这个方程，我们通过缩放时间单位的方式将其规范化，使粘度的数值为1，描述了缓慢的粘性流动，如蜂蜜或岩浆。本教程的主要目的是展示如何用惩罚方法，以弱的方式将速度场条件 $\mathbf{u} = \mathbf{g}$ 强加于非匹配的 $\Gamma$ 。关于斯托克斯问题的更广泛的讨论，包括体力、不同的边界条件和解决策略，可以在步骤22中找到。

让我们开始单独考虑整个域 $\Omega$ 中的斯托克斯问题。我们寻找一个速度场 $\mathbf{u}$ 和一个压力场 $p$ ，满足斯托克斯方程和 $\partial\Omega$ 上的同质边界条件。

斯托克斯方程的微弱形式首先通过将其写成矢量形式而得到

@f{eqnarray*}
  \begin{pmatrix}
    {-\Delta \textbf{u} + \nabla p}
    \\
    {-\textrm{div}\;\textbf{u}}
  \end{pmatrix}
  =
  \begin{pmatrix}
  0
  \\
  0
  \end{pmatrix},


@f}

从左边开始与一个矢量值测试函数 $\phi = \begin{pmatrix}\textbf{v} \\ q\end{pmatrix}$ 形成点积，并在域 $\Omega$ 上进行积分，得到以下一组方程。

@f{eqnarray*}
  (\mathrm v,


   -\Delta \textbf{u} + \nabla p)_{\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  0


@f}

这对所有的测试函数都必须成立  $\phi = \begin{pmatrix}\textbf{v}
\\ q\end{pmatrix}$  。


通过部分积分并利用 $\partial\Omega$ 的边界条件，我们得到以下变分问题。

@f{eqnarray*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega} - (\textrm{div}\; \textbf{v}, p)_{\Omega}


 - (q, \textrm{div}\; \textbf{u})_{\Omega}&=& 0


@f}



其中 $(\cdot, \cdot)_{\Omega}$ 代表 $L^2$ 标量积。这也是步骤22中使用的变异形式。

这个变分公式没有考虑到嵌入域。与step-60相反，我们并不强行执行 $\textbf{u}$ 对 $\Gamma$ 的约束，而是通过惩罚项弱行执行这些约束。

对这种弱强加边界条件的分析取决于 $\Gamma$ 的spacedim-dimensional度量是正的（如果`dim`等于`spacedim`）或零（如果`dim`小于`spacedim`）。我们讨论这两种情况。




<a name="Codimensiononecase"></a><h4>Co-dimension one case</h4>


在这种情况下，我们假设 $\Gamma$ 是实际叶轮的边界，即嵌入二维域的封闭曲线或三维域的封闭表面。这种方法的思路首先是考虑在 $\Gamma$ 上弱加迪里切特边界条件，遵循尼采方法。这是通过在流体域上使用以下修改后的公式来实现的，其中没有对 $\Gamma$ 上的测试函数施加强条件。

@f{multline*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega\setminus\Omega^{\text{imp}}} - (\textrm{div}\;  \textbf{v}, p)_{\Omega\setminus\Omega^{\text{imp}}}


  - (q, \textrm{div}\; \textbf{u})_{\Omega\setminus\Omega^{\text{imp}}} \\


  - (\textbf{v},\nabla \textbf{u} \cdot \textbf{n})_{\Gamma}
  + (\textbf{v}\cdot \textbf{n},p)_{\Gamma} \\


 -  (\nabla\textbf{v}\cdot \textbf{n},\textbf{u})_{\Gamma}
 + (q, \textbf{u} \cdot \textbf{n})_{\Gamma}
 + \beta (\textbf{v},\textbf{u})_{\Gamma} \\
=  - (\nabla\textbf{v}\cdot \textbf{n},\textbf{g})_{\Gamma} + (q, \textbf{g} \cdot \textbf{n})_{\Gamma}
 + \beta (\textbf{v},\textbf{g})_{\Gamma}.


@f}



过 $\Gamma$ 的积分是低维积分。可以证明（见 @cite Freund1995 ），存在一个正的常数 $C_1$ ，所以如果 $\beta > C_1$ ，边界的弱强加将是一致和稳定的。在 $\Gamma$ 上的前两个附加积分（上式中的第二行）在通过部分积分后自然出现，此时我们不假设 $\mathbf{v}$ 在 $\Gamma$ 上是零。

上述方程中的第三行包含两个项，是为了确保弱形式的一致性而添加的，还有一个稳定项，是为了强制执行边界条件，其误差与近似误差一致。一致性项和稳定项是用实际的边界数据添加到右手边的  $\mathbf{g}$  。

当 $\mathbf{u}$ 满足 $\Gamma$ 上的条件 $\mathbf{u}=\mathbf{g}$ 时， $\Gamma$ 上的所有一致性和稳定性积分都被抵消，就剩下斯托克斯流的通常弱形式，也就是说，上述表述是一致的。

我们注意到，可以使用另一种（非对称的）表述方式。

@f{multline*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega\setminus\Omega^{\text{imp}}} -  (\textrm{div}\;  \textbf{v}, p)_{\Omega\setminus\Omega^{\text{imp}}}


  - (q, \textrm{div}\; \textbf{u})_{\Omega\setminus\Omega^{\text{imp}}} \\


  -(\textbf{v},\nabla \textbf{u} \cdot \textbf{n})_{\Gamma}
  + (\textbf{v}\cdot \textbf{n},p)_{\Gamma} \\
   +(\nabla\textbf{v}\cdot \textbf{n},\textbf{u})_{\Gamma}


 - (q, \textbf{u} \cdot \textbf{n})_{\Gamma}
 + \beta (\textbf{v},\textbf{u})_{\Gamma} \\
=   (\nabla\textbf{v}\cdot \textbf{n},\textbf{g})_{\Gamma} - (q, \textbf{g} \cdot \textbf{n})_{\Gamma}
 + \beta (\textbf{v},\textbf{g})_{\Gamma}.


@f}

注意第三行和第四行的第一项的不同符号。在这种情况下，稳定性和一致性条件成为 $\beta > 0$  。在对称情况下， $\beta$ 的值取决于 $h$ ，一般来说，它被选择为 $\beta = C h^{-1} $ ， $h$ 是衡量被整合面的大小， $C$ 是一个常数，以便 $1 \leq C \leq 10$  。这就像人们通常使用Nitsche惩罚方法来执行Dirichlet边界条件一样。

另一方面，非对称方法与非连续Galerkin方法的非对称内部惩罚方法（"NIPG "方法 @cite Riviere1999 ）的连续性的执行方式有关。即使非对称情况在稳定参数的可能选择方面似乎更有优势，我们还是选择了对称离散化，因为在这种情况下，可以证明对偶问题也是一致的，导致解决方案不仅能量准则以正确的顺序收敛，而且其 $L^2$ 准则也是如此。此外，得到的矩阵仍然是对称的。

上述表述是在假设领域被精确离散的情况下进行的。然而，如果叶轮的变形是一个刚体运动，就有可能人为地将斯托克斯问题的解扩展到螺旋桨本身，因为刚体运动也是斯托克斯问题的解。我们的想法是在 $\Omega^{\text{imp}}$ 内解决同样的问题，在 $\Gamma$ 上施加同样的边界条件，使用同样的惩罚技术，并用在 $\Omega$ 上全局连续的测试函数 $\mathbf{v}$ 来测试。

这导致了以下（中间）配方。

@f{multline*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega} - (\textrm{div}\;  \textbf{v}, p)_{\Omega}


  - (q, \textrm{div}\; \textbf{u})_{\Omega} \\


  - (\textbf{v},  \lbrack \nabla \textbf{u} \rbrack \cdot \textbf{n})_{\Gamma}
  + (\textbf{v}\cdot \textbf{n},\lbrack p \rbrack )_{\Gamma} \\


 -  (\lbrack \nabla\textbf{v} \rbrack \cdot \textbf{n},\textbf{u})_{\Gamma}
 + (\lbrack q \rbrack, \textbf{u} \cdot n)_{\Gamma}
 + 2\beta (\textbf{v},\textbf{u})_{\Gamma} \\
=  - (\lbrack \nabla\textbf{v}\rbrack\cdot \textbf{n},\textbf{g})_{\Gamma} + (\lbrack q\rbrack, \textbf{g} \cdot n)_{\Gamma}
 + 2\beta (\textbf{v},\textbf{g})_{\Gamma},


@f}

其中跳跃项，用 $\lbrack \cdot \rbrack$ 表示，是相对于法向量 $\textbf{n}$ 的一个固定方向计算的。2的因子出现在 $\beta$ 前面，因为我们看到 $\Gamma$ 的每一部分两次，一次来自流体内部，一次来自在其中移动的障碍物。对于 $\Gamma$ 上的所有其他积分，我们对 $\Gamma$ 的每一部分都访问了两次，但符号相反，因此得到的是跳跃项）。

这里我们注意到，与不连续的Galerkin方法不同，测试和试验函数在 $\Gamma$ 中是连续的。此外，如果 $\Gamma$ 不与单元边界对齐，所有的跳跃项也是零，因为一般来说，有限元函数空间在每个单元内都是平滑的，如果 $\Gamma$ 只在有限的几个点上切过一个单元与它的边界相交，除了稳定化的贡献外， $\Gamma$ 上的所有贡献都可以从公式中忽略掉，导致以下变量公式的最终形式。

@f{multline*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega} - (\textrm{div}\;  \textbf{v}, p)_{\Omega}


  - (q, \textrm{div}\; \textbf{u})_{\Omega}  + 2\beta (\textbf{v},\textbf{u})_{\Gamma} \\
=  2\beta (\textbf{v},\textbf{g})_{\Gamma}.


@f}



在step-60中，约束条件的施加需要以拉格朗日乘数的形式增加新的变量。本教程程序不存在这种情况。使用Nitsche方法施加边界条件只修改了系统矩阵和右手边，没有增加额外的未知数。然而，嵌入域上的速度矢量 $\textbf{u}$ 不会与规定的速度 $\textbf{g}$ 完全匹配，而只是达到一个数值误差，这个误差与有限元方法的插值误差相同。此外，与第60步一样，我们仍然需要在不匹配的嵌入网格上进行积分，以构建对 $\Gamma$ 施加边界条件的必要边界项。




<a name="Codimensionzerocase"></a><h4>Co-dimension zero case</h4>


在这种情况下， $\Gamma$ 具有相同的尺寸，但被嵌入到 $\Omega$ 中。我们可以把它看作是一个在流体中移动的厚物体。在 $\mathcal{L}^2$ 惩罚的情况下，额外的惩罚项可以被解释为 $\Gamma$ 内的达西项，结果是。

@f{eqnarray*}
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega} - & (\textrm{div}\;  \textbf{v}, p)_{\Omega}


  - (q, \textrm{div}\; \textbf{u})_{\Omega}  + \beta (\textbf{v},\textbf{u})_{\Gamma}
=  \beta (\textbf{v},\textbf{g})_{\Gamma}.


@f}



这里，对 $\Gamma$ 的积分只是对部分体积的积分。因此， $\mathcal{L}^2$ 的惩罚包括增加一个体积项，约束流体的速度与 $\Gamma$ 内刚体的速度保持一致。在这种情况下， $\beta$ 必须被选择得足够大，以确保 $\Gamma$ 中的迪里希特边界条件得到充分尊重，但也不能太高，以保持系统矩阵的适当调节。

一个 $\mathcal{H}^1$ 的惩罚可以用类似的方式构建，在惩罚中加入一个粘性成分，以阻尼 $\Gamma$ 内的速度梯度。

@f{eqnarray*}{
(\nabla \textbf{v}, \nabla \textbf{u})_{\Omega} - & (\textrm{div}\;  \textbf{v}, p)_{\Omega}


  - (q, \textrm{div}\; \textbf{u})_{\Omega}
  + \beta_1 (\textbf{v},\textbf{u})_{\Gamma}
  + \beta_2 (\nabla \textbf{v}, \nabla \textbf{u})_{\Gamma}
=  \beta_1 (\textbf{v},\textbf{g})_{\Gamma}
+ \beta_2 (\nabla \textbf{v}, \nabla \textbf{g})_{\Gamma}.


@f}



请注意， $L^2$ 的惩罚（`dim`等于`spacedim`）和Nitsche的惩罚（`dim`等于`spacedim-1`）导致了完全相同的数值实现，这要感谢deal.II的独立维度能力。




<a name="Representationofand"></a><h4>Representation of Ω and Γ</h4>


在本教程中，嵌入网格 $\Gamma$ 和嵌入网格都是用 parallel::distributed::Triangulation. 来描述的。这两个三角形可以通过GridGenerator命名空间中的函数来建立，或者通过读取其他应用程序（例如GMSH，见步骤-49的讨论）产生的网格文件来建立。这比之前在第60步中的做法略微通用了一些。

无论是在 "dim=spacedim "还是 "dim<spacedim "的情况下，增加沉没边界法，只是在系统矩阵和系统的右手边引入了额外的项，这些项是在 $\Gamma$ 上积分的结果。这并没有改变必须解决的问题的变量数量。因此，挑战与必须进行的积分有关  $\Gamma$  。

在有限元中，我们将这个积分分成来自用于离散化 $\Gamma$ 的所有单元的贡献，我们将 $K$ 上的积分转换为参考元素 $\hat K$ 上的积分，其中 $F_{K}$ 是 $\hat K$ 到 $K$ 的映射，并使用正交公式计算 $\hat K$ 上的积分。比如说。

\f[
\beta (\textbf{v},\textbf{u})_{\Gamma} =  \sum_{K\in \Gamma} \int_{\hat K}
\hat{\textbf{u}}(\hat x) (\textbf{v} \circ F_{K}) (\hat x) J_K (\hat x) \mathrm{d} \hat x =
\sum_{K\in \Gamma} \sum_{i=1}^{n_q}  \big(\hat{\textbf{u}}(\hat x_i)  (\textbf{v} \circ F_{K}) (\hat x_i) J_K (\hat x_i) w_i \big)
\f]

计算这个和是不容易的，因为我们必须评估 $(v_j \circ F_{K})
(\hat x_i)$  。一般来说，如果 $\Gamma$ 和 $\Omega$ 没有对齐，那么 $y_i = F_{K}(\hat x_i)$ 这个点相对于 $\Omega$ 来说是完全任意的，除非我们想出一个办法，将 $V_h(\Omega)$ 的所有基函数插在 $\Omega$ 上的一个任意点上，否则我们无法计算出需要的积分。


要评估 $(v_j \circ F_{K}) (\hat x_i)$ ，需要采取以下步骤（如下图所示）。

- 对于 $\Gamma$ 中的一个给定单元 $K$ ，计算实点 $y_i \dealcoloneq F_{K} (\hat
x_i)$ ，其中 $x_i$ 是用于 $K
\subseteq \Gamma$ 上的积分的正交点之一。这是最容易的部分。   FEValues::quadrature_point() 给了我们所有正交点的实空间位置。

- 找到 $\Omega$ 中 $y_i$ 所在的单元。我们将称这个元素为 $T$  。

- 找到 $T$ 内 $y_i$ 的参考坐标。为此，我们需要将参考元素 $\hat T$ 转换为元素 $T$ ： $\hat y_i = G^{-1}_{T} (y_i)$ 的映射 $G_T$ 的逆映射。

- 评估  $v_j$  网格在此点  $\hat y_i$  的基函数  $\Omega$  。这也是比较简单的，使用FEValues。


<p align="center"> <img src="https://www.dealii.org/images/steps/developer/step-60.C_interpolation.png" alt="">  </p> 

在步骤60中，上述第二至第四步是通过依次调用来计算的。

-  GridTools::find_active_cell_around_point(),  后面是

-  Mapping::transform_real_to_unit_cell().  然后我们

- 构建一个自定义的正交公式，包含参考单元格中的点，然后

- 构建一个FEValues对象，具有给定的正交公式，并以第一步中获得的单元格为初始化。

虽然这种方法对目前的情况是可行的，但它并不适合于使用分布式三角形的平行模拟。事实上，由于嵌入域 $\Gamma$ 单元上的正交点的位置与嵌入三角形的位置不一致，而且 $\Gamma$ 是不断移动的，这就要求代表 $\Gamma$ 的三角形被完整地存储在所有处理器中。随着处理器的数量和 $\Gamma$ 中单元格数量的增加，这将导致内存方面的严重瓶颈。因此，在这一步骤中寻求一种替代策略。




<a name="Usingparticlestotrack"></a><h4>Using particles to track Γ</h4>


请记住，对于惩罚法（ $\mathcal{L}^2$ 或 $\mathcal{H}^1$ ）和尼采法，我们要计算的是由正交近似的积分。也就是说，我们需要计算

\f[
\beta (\textbf{v},\textbf{u})_{\Gamma} =
\sum_{K\in \Gamma} \sum_{i=1}^{n_q}  \big(\hat{\textbf{u}}(\hat x_i)  (\textbf{v} \circ F_{K}) (\hat x_i) J_K (\hat x_i) w_i \big)
\f] 如果你跟随上面的讨论，那么你会记得  $\textbf{u}$  和  $\textbf{v}$  是定义在流体网格上的形状函数。唯一定义在实体网格上的东西是。   $F_K(\hat x_i)$  ，是实体单元上正交点的位置，是 $\Gamma$ 的一部分， $J_K$ 是其雅各布系数的行列式， $w_i$ 是相应的正交权值。

现在要认识到的重要部分是这样的。   $w_i$ 是正交公式的一个属性，不随时间变化。此外， $F_K$ 的雅各布矩阵本身随着固体障碍物在流体中的移动而变化，但由于固体被认为是非变形的（它只是平移和旋转，但不扩张），雅各布矩阵的行列式保持不变。因此，乘积 $J_K(\hat x_i) w_i$ （我们通常用`JxW`表示）在每个正交点上都保持不变。因此，我们唯一需要跟踪的是位置 $x_i=F_K(\hat x_i)$ --但这些位置随着实体域的速度移动。

换句话说，我们实际上根本不需要保留实体网格。我们所需要的只是位置 $x_i(t)$ 和相应的`JxW`值。由于这两个属性都是附着在实体材料上的点属性（或点向量），它们可以被理想化为一组不相连的无限小的 "粒子"，它们随着实体的运动携带所需的`JxW`信息。deal.II有能力以ParticleHandler类的形式在大规模并行计算中分配和存储这样一组粒子（关于实现的细节见 @cite GLHPW2018  ），我们将在本教程中使用这一功能。

因此，本步骤采取的方法如下。

- 为域名  $\Gamma$  创建一个  parallel::distributed::Triangulation  。

- 在 Particles::Particle 上的正交点位置创建 $\Gamma$  。

- 调用 Particles::ParticleHandler::insert_global_particles() 函数，将粒子分配到各个处理器上，*遵循实体三角形*的做法。

- 将 "JxW "值作为一个 "属性 "附加到每个 Particles::Particle 对象。

这种结构的生成相对来说比较昂贵，但是每次模拟必须只生成一次。一旦 Particles::ParticleHandler 被生成，并且所需的信息被附加到粒子上，就可以利用粒子在ParticleHandler内按单元分组的事实，对 $\Gamma$ 进行积分，使我们能够。

- 在 $\Omega$ 中至少包含一个粒子的所有单元格上循环操作

- 循环处理给定单元中的所有粒子

- 计算积分并填充全局矩阵。

由于 Particles::ParticleHandler 可以管理粒子从一个处理器到另一个处理器的交换，嵌入的三角形可以通过位移粒子而被移动或变形。与这种位移相关的唯一约束是，颗粒的位移距离不应大于一个单元的大小。这是因为这是 Particles::ParticleHandler 能够追踪离开当前单元的粒子现在所处的单元的极限。

一旦整个问题（斯托克斯问题和沉没边界施加）被集合起来，最后的鞍点问题由迭代求解器解决，应用于舒尔补数 $S$ （其构造例如在步骤22中描述），我们使用LinearOperator类构造 $S$ 。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们在这里解决的问题是对斯托克斯流的时间可逆性的证明。这在科学教育实验中经常用泰勒-库伊特流和染料液滴来说明，在流体以周期性的方式位移后，染料液滴又恢复到原来的形状。

@htmlonly


<iframe width="560" height="315" src="https://www.youtube.com/embed/p08_KlTKP50" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>


@endhtmlonly



在目前的问题中，一个非常粘稠的流体被一个叶轮的旋转所搅动，在二维中，叶轮被一个矩形网格所模拟。叶轮旋转了一定的圈数，之后流动被逆转，从而在相反的方向上进行相同圈数的旋转。我们回顾一下，由于斯托克斯方程是自交的，蠕动流是可逆的。因此，如果叶轮运动在相反的方向上被逆转，流体应该回到其原来的位置。在本例中，我们通过插入一圈被动示踪剂颗粒来说明这一点，这些颗粒被流体平移并返回到原来的位置，从而证明了流动的时间可逆性。




<a name="Morereferences"></a><h3> More references</h3>


本教程程序使用了一些关于对流体内部的非匹配界面施加速度条件的技术。要了解更多的背景材料，你可能要查阅以下参考资料。   @cite Freund1995  ,  @cite Angot1999  ,  @cite Glowinski1999  ,  @cite Boffi2008  ,  @cite Heltai2012  。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name=""></a> 
 * @sect3{Include files}  其中大部分已经在其他地方介绍过了，我们只对新的部分进行评论。靠近顶部的开关允许在 PETSc 和 Trilinos 线性代数功能之间进行选择，这与  step-40  和  step-50  中的开关类似。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/timer.h> 
 * 
 * #include <deal.II/lac/block_linear_operator.h> 
 * #include <deal.II/lac/generic_linear_algebra.h> 
 * #include <deal.II/lac/linear_operator.h> 
 * #include <deal.II/lac/linear_operator_tools.h> 
 * 
 * #define FORCE_USE_OF_TRILINOS 
 * 
 * namespace LA 
 * { 
 * #if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \ 
 *   !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS)) 
 *   using namespace dealii::LinearAlgebraPETSc; 
 * #  define USE_PETSC_LA 
 * #elif defined(DEAL_II_WITH_TRILINOS) 
 *   using namespace dealii::LinearAlgebraTrilinos; 
 * #else 
 * #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required 
 * #endif 
 * } // namespace LA 
 * 
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/base/parameter_acceptor.h> 
 * #include <deal.II/base/parsed_function.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/distributed/grid_refinement.h> 
 * #include <deal.II/distributed/solution_transfer.h> 
 * #include <deal.II/distributed/tria.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_nothing.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_fe_field.h> 
 * #include <deal.II/fe/mapping_q.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_in.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/petsc_precondition.h> 
 * #include <deal.II/lac/petsc_solver.h> 
 * #include <deal.II/lac/petsc_sparse_matrix.h> 
 * #include <deal.II/lac/petsc_vector.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/solver_minres.h> 
 * #include <deal.II/lac/sparsity_tools.h> 
 * #include <deal.II/lac/vector.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * @endcode
 * 
 * 这些是关于  step-60  的唯一新的包含文件。在本教程中，实体和流体之间的非匹配耦合是通过一个中间数据结构来计算的，该结构记录了实体的正交点在流体网格中的位置如何演变。这个数据结构需要跟踪描述实体域的每个单元上的正交点的位置、正交权重，如果实体域是同维度的，还需要跟踪每个点的法向量。
 * 

 * 
 * Deal.II通过ParticleHandler类在Particles命名空间中提供这些设施。ParticleHandler是一个允许你管理粒子集合的类（类型为 Particles::Particle), 的对象，代表具有一些附加属性（如id）的点的集合，漂浮在一个 parallel::distributed::Triangulation. 命名空间中的方法和类允许人们轻松实现Particle-In-Cell方法和在分布式三角形上的粒子追踪。
 * 

 * 
 * 我们 "滥用 "这个数据结构来存储嵌入周围流体网格中的实体正交点的位置信息，包括积分权重，以及可能的表面法线。我们之所以使用这个额外的数据结构，是因为实体网格和流体网格可能是不重叠的，如果我们使用两个独立的三角计算对象，那么它们将独立地分布在并行进程中。
 * 

 * 
 * 为了耦合这两个问题，我们依靠ParticleHandler类，在每个粒子中存储一个实体正交点的位置（一般来说，它不与任何流体正交点对齐），它的权重，以及耦合这两个问题可能需要的任何其他信息。这些位置然后与固体叶轮的（规定）速度一起传播。
 * 

 * 
 * 固体正交点的所有权最初是从固体网格本身的MPI分区中继承的。这样产生的粒子后来通过ParticleHandler类的方法分配到流体网格中。这允许MPI进程之间透明地交换关于流体单元和实体正交点之间的重叠模式的信息。
 * 

 * 
 * 
 * @code
 * #include <deal.II/particles/data_out.h> 
 * #include <deal.II/particles/generators.h> 
 * #include <deal.II/particles/particle_handler.h> 
 * #include <deal.II/particles/utilities.h> 
 * 
 * @endcode
 * 
 * 在生成网格时，我们允许从文件中读取它，如果deal.II已经建立了OpenCASCADE支持，我们也允许读取CAD文件，并将它们作为网格的流形描述符（参见 step-54 对OpenCASCADE命名空间中的各种流形描述符的详细描述）。
 * 

 * 
 * 
 * @code
 * #include <deal.II/opencascade/manifold_lib.h> 
 * #include <deal.II/opencascade/utilities.h> 
 * #ifdef DEAL_II_WITH_OPENCASCADE 
 * #  include <TopoDS.hxx> 
 * #endif 
 * 
 * #include <cmath> 
 * #include <fstream> 
 * #include <iostream> 
 * #include <memory> 
 * 
 * namespace Step70 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameterhandling"></a> 
 * <h3>Run-time parameter handling</h3>
 * 

 * 
 * 与我们在 step-60 中所做的类似，我们建立了一个持有我们问题的所有参数的类，并从ParameterAcceptor类中派生出来以简化参数文件的管理和创建。
 * 

 * 
 * ParameterAcceptor范式要求所有参数都可以被ParameterAcceptor方法写入。为了避免出现很难追踪的错误（比如写成`time = 0`而不是`time == 0`），我们在一个外部类中声明所有的参数，该类在实际的`StokesImmersedProblem`类之前被初始化，并将其作为`const`引用传递给主类。
 * 

 * 
 * 该类的构造函数负责该类的成员与ParameterHandler中的相应条目之间的连接。由于使用了 ParameterHandler::add_parameter() 方法，这种连接是微不足道的，但要求这个类的所有成员都是可写的。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim = dim> 
 *   class StokesImmersedProblemParameters : public ParameterAcceptor 
 *   { 
 *   public: 
 *     StokesImmersedProblemParameters(); 
 * 
 * @endcode
 * 
 * 然而，由于这个类将作为一个`const`引用传递给StokesImmersedProblem类，我们必须确保我们仍然可以在这里定义的Function类派生的对象中正确设置时间。为了做到这一点，我们声明 `StokesImmersedProblemParameters::rhs` 和 `StokesImmersedProblemParameters::angular_velocity` 成员都是 "可变 "的，并定义以下的小辅助方法，将它们的时间设置为正确的值。
 * 

 * 
 * 
 * @code
 *     void set_time(const double &time) const 
 *     { 
 *       rhs.set_time(time); 
 *       angular_velocity.set_time(time); 
 *     } 
 * 
 * @endcode
 * 
 * 该类的其余部分主要由描述模拟及其离散化细节的成员变量组成。下面的参数是关于输出的位置、空间和时间离散化（默认是 $Q_2\times Q_1$ Taylor-Hood离散化，它使用2度的多项式来计算速度），以及在我们再次生成图形输出之前应该经过多少时间步长。
 * 

 * 
 * 
 * @code
 *     std::string output_directory = "."; 
 * 
 *     unsigned int velocity_degree = 2; 
 * 
 *     unsigned int number_of_time_steps = 501; 
 *     double       final_time           = 1.0; 
 * 
 *     unsigned int output_frequency = 1; 
 * 
 * @endcode
 * 
 * 我们允许每个网格独立地被细化。在本教程中，固体网格上没有解决物理问题，其速度被作为基准点给出。然而，在本教程中加入一些弹性模型，并将其转化为一个完全成熟的FSI求解器是相对简单的。
 * 

 * 
 * 
 * @code
 *     unsigned int initial_fluid_refinement      = 5; 
 *     unsigned int initial_solid_refinement      = 5; 
 *     unsigned int particle_insertion_refinement = 3; 
 * 
 * @endcode
 * 
 * 为了提供对流体领域的粗略描述，我们使用extract_rtree_level()方法，该方法适用于流体三角结构中每个局部拥有的单元的边界盒树。树的级别越高，提取的边界盒数量就越多，对流体领域的描述也就越准确。然而，大量的边界盒也意味着巨大的通信成本，因为边界盒的收集是由所有进程收集的。
 * 

 * 
 * 
 * @code
 *     unsigned int fluid_rtree_extraction_level = 1; 
 * 
 * @endcode
 * 
 * 方程中使用的唯一两个数值参数是流体的粘度，以及Nitsche公式中使用的惩罚项 $\beta$ 。
 * 

 * 
 * 
 * @code
 *     double viscosity    = 1.0; 
 *     double penalty_term = 100; 
 * 
 * @endcode
 * 
 * 默认情况下，我们创建一个没有着色的hyper_cube，并且我们使用同质的Dirichlet边界条件。在这个集合中，我们存储了设置边界条件时要使用的边界ID。
 * 

 * 
 * 
 * @code
 *     std::list<types::boundary_id> homogeneous_dirichlet_ids{0}; 
 * 
 * @endcode
 * 
 * 我们在此说明另一种从参数文件中创建三角形的方法，使用 GridGenerator::generate_from_name_and_arguments(), ，该方法接收GridGenerator命名空间中的函数名称，其参数为一个字符串，代表参数的元组。
 * 

 * 
 * 在 Patterns::Tools::Convert 类中详细解释了将参数从字符串解析成字符串的机制，该类用于将字符串翻译成大多数基本STL类型（向量、映射、图元）和基本deal.II类型（点、张量、BoundingBox等）。
 * 

 * 
 * 一般来说，可以用等级1的统一元素表示的对象（即 std::vector<double>,  Point<dim>,  std::set<int>,  等）是用逗号分开的。额外的等级采取分号，允许你将字符串解析为 `std::vector<std::vector<double>>`, 或例如 `std::vector<Point<dim>>`, 类型的对象，如`0.0, 0.1; 0.1, 0.2`。这个字符串可以被解释为两个Point对象的向量，或者一个双数向量的向量。
 * 

 * 
 * 当条目不统一时，比如在元组的情况下，我们用冒号来分隔各个条目。例如，像`5: 0.1, 0.2`这样的字符串可以用来解析一个类型为 `std::pair<int,  Point<2>>的对象或者一个 `std::tuple<int,  的对象。
 * std::vector<double>>`.  
 * 

 * 
 * 在我们的例子中，大多数参数是点对象（代表中心、角、细分元素等）、整数值（细分数量）、双倍值（半径、长度等）或布尔选项（如许多GridGenerator函数采取的`colorize`选项）。
 * 

 * 
 * 在下面的例子中，我们设置了合理的默认值，但这些值可以在运行时通过选择GridGenerator命名空间的任何其他支持的函数来改变。如果GridGenerator函数失败，本程序将把网格的名称解释为vtk网格文件名，把参数解释为从manifold_id到描述域的几何形状的CAD文件的映射。每个CAD文件都将被分析，并根据CAD文件本身的内容生成OpenCASCADE命名空间的Manifold。
 * 

 * 
 * 为了尽可能的通用，我们对每个生成的网格都这样做：流体网格、固体网格，但也包括使用三角法生成的示踪粒子。
 * 

 * 
 * 
 * @code
 *     std::string name_of_fluid_grid       = "hyper_cube"; 
 *     std::string arguments_for_fluid_grid = "-1: 1: false"; 
 *     std::string name_of_solid_grid       = "hyper_rectangle"; 
 *     std::string arguments_for_solid_grid = spacedim == 2 ? 
 *                                              "-.5, -.1: .5, .1: false" : 
 *                                              "-.5, -.1, -.1: .5, .1, .1: false"; 
 *     std::string name_of_particle_grid = "hyper_ball"; 
 *     std::string arguments_for_particle_grid = 
 *       spacedim == 2 ? "0.3, 0.3: 0.1: false" : "0.3, 0.3, 0.3 : 0.1: false"; 
 * 
 * @endcode
 * 
 * 同样地，我们允许不同的局部细化策略。特别是，我们限制了细化水平的最大数量，以控制流体网格的最小尺寸，并保证它与实体网格兼容。细化级数的最小值也得到了控制，以确保在流动的大部分地区有足够的精度。此外，我们根据流体速度场的标准误差估计器进行局部细化。
 * 

 * 
 * 我们允许用户选择两种最常见的细化策略，即 "fixed_number "或 "fixed_fraction"，这两种策略参考了 GridRefinement::refine_and_coarsen_fixed_fraction() 和 GridRefinement::refine_and_coarsen_fixed_number(). 方法。
 * 

 * 
 * 细化可以每隔几个时间步骤进行一次，而不是连续进行，我们通过`细化_频率`参数来控制这个值。
 * 

 * 
 * 
 * @code
 *     int          max_level_refinement = 8; 
 *     int          min_level_refinement = 5; 
 *     std::string  refinement_strategy  = "fixed_fraction"; 
 *     double       coarsening_fraction  = 0.3; 
 *     double       refinement_fraction  = 0.3; 
 *     unsigned int max_cells            = 20000; 
 *     int          refinement_frequency = 5; 
 * 
 * @endcode
 * 
 * 最后，以下两个函数对象被用来控制斯托克斯流的源项和我们移动固体体的角速度。在一个更现实的模拟中，实体速度或其变形将来自于实体域上的辅助问题的解决。在这个例子中，我们把这部分放在一边，只是在浸没的固体上沿Z轴施加一个固定的旋转速度场，由一个可以在参数文件中指定的函数来控制。
 * 

 * 
 * 
 * @code
 *     mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> rhs; 
 *     mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> 
 *       angular_velocity; 
 *   }; 
 * 
 * @endcode
 * 
 * 还有一个任务就是声明我们在输入文件中可以接受哪些运行时参数。我们将这些参数分成不同的类别，把它们放在ParameterHandler类的不同部分。我们首先在全局范围内声明StokesImmersedProblem使用的所有全局参数。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   StokesImmersedProblemParameters<dim, 
 *                                   spacedim>::StokesImmersedProblemParameters() 
 *     : ParameterAcceptor("Stokes Immersed Problem/") 
 *     , rhs("Right hand side", spacedim + 1) 
 *     , angular_velocity("Angular velocity") 
 *   { 
 *     add_parameter( 
 *       "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1)); 
 * 
 *     add_parameter("Number of time steps", number_of_time_steps); 
 *     add_parameter("Output frequency", output_frequency); 
 * 
 *     add_parameter("Output directory", output_directory); 
 * 
 *     add_parameter("Final time", final_time); 
 * 
 *     add_parameter("Viscosity", viscosity); 
 * 
 *     add_parameter("Nitsche penalty term", penalty_term); 
 * 
 *     add_parameter("Initial fluid refinement", 
 *                   initial_fluid_refinement, 
 *                   "Initial mesh refinement used for the fluid domain Omega"); 
 * 
 *     add_parameter("Initial solid refinement", 
 *                   initial_solid_refinement, 
 *                   "Initial mesh refinement used for the solid domain Gamma"); 
 * 
 *     add_parameter("Fluid bounding boxes extraction level", 
 *                   fluid_rtree_extraction_level, 
 *                   "Extraction level of the rtree used to construct global " 
 *                   "bounding boxes"); 
 * 
 *     add_parameter( 
 *       "Particle insertion refinement", 
 *       particle_insertion_refinement, 
 *       "Refinement of the volumetric mesh used to insert the particles"); 
 * 
 *     add_parameter( 
 *       "Homogeneous Dirichlet boundary ids", 
 *       homogeneous_dirichlet_ids, 
 *       "Boundary Ids over which homogeneous Dirichlet boundary conditions are applied"); 
 * 
 * @endcode
 * 
 * 下一节专门介绍用于创建各种网格的参数。我们将需要三种不同的三角形。流体网格 "用于定义流体领域，"固体网格 "用于定义固体领域，"粒子网格 "用于分布一些示踪粒子，这些粒子随速度漂移，只作为被动示踪物使用。
 * 

 * 
 * 
 * @code
 *     enter_subsection("Grid generation"); 
 *     { 
 *       add_parameter("Fluid grid generator", name_of_fluid_grid); 
 *       add_parameter("Fluid grid generator arguments", arguments_for_fluid_grid); 
 * 
 *       add_parameter("Solid grid generator", name_of_solid_grid); 
 *       add_parameter("Solid grid generator arguments", arguments_for_solid_grid); 
 * 
 *  
 *       add_parameter("Particle grid generator arguments", 
 *                     arguments_for_particle_grid); 
 *     } 
 *     leave_subsection(); 
 * 
 *     enter_subsection("Refinement and remeshing"); 
 *     { 
 *       add_parameter("Refinement step frequency", refinement_frequency); 
 *       add_parameter("Refinement maximal level", max_level_refinement); 
 *       add_parameter("Refinement minimal level", min_level_refinement); 
 *       add_parameter("Refinement strategy", 
 *                     refinement_strategy, 
 *                     "", 
 *                     this->prm, 
 *                     Patterns::Selection("fixed_fraction|fixed_number")); 
 *       add_parameter("Refinement coarsening fraction", coarsening_fraction); 
 *       add_parameter("Refinement fraction", refinement_fraction); 
 *       add_parameter("Maximum number of cells", max_cells); 
 *     } 
 *     leave_subsection(); 
 * 
 * @endcode
 * 
 * 最后的任务是修正右侧函数的默认尺寸，并定义一个有意义的默认角速度而不是零。
 * 

 * 
 * 
 * @code
 *     rhs.declare_parameters_call_back.connect([&]() { 
 *       Functions::ParsedFunction<spacedim>::declare_parameters(this->prm, 
 *                                                               spacedim + 1); 
 *     }); 
 *     angular_velocity.declare_parameters_call_back.connect([&]() { 
 *       this->prm.set("Function expression", 
 *                     "t < .500001 ? 6.283185 : -6.283185"); 
 *     }); 
 *   } 
 * 
 * @endcode
 * 
 * 一旦角速度被提供为一个函数对象，我们就通过下面这个派生自函数类的类来重建点状实体速度。它通过假设实体以给定的角速度绕原点（或3D中的 $z$ 轴）旋转，提供实体在给定位置的速度值。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   class SolidVelocity : public Function<spacedim> 
 *   { 
 *   public: 
 *     static_assert(spacedim > 1, 
 *                   "Cannot instantiate SolidVelocity for spacedim == 1"); 
 * 
 *     SolidVelocity(const Functions::ParsedFunction<spacedim> &angular_velocity) 
 *       : angular_velocity(angular_velocity) 
 *     {} 
 * 
 *     virtual double value(const Point<spacedim> &p, 
 *                          unsigned int           component = 0) const override 
 *     { 
 *       Tensor<1, spacedim> velocity; 
 * 
 * @endcode
 * 
 * 我们假设角速度是沿Z轴方向的，也就是说，我们把实际的角速度模拟成二维旋转，而不考虑`spacedim`的实际值。
 * 

 * 
 * 
 * @code
 *       const double omega = angular_velocity.value(p); 
 *       velocity[0]        = -omega * p[1]; 
 *       velocity[1]        = omega * p[0]; 
 * 
 *       return velocity[component]; 
 *     } 
 * 
 *   private: 
 *     const Functions::ParsedFunction<spacedim> &angular_velocity; 
 *   }; 
 * 
 * @endcode
 * 
 * 同样地，我们假设固体的位置可以在每个时间步长明确地计算出来，利用角速度的知识。我们计算固体粒子的确切位置，假定固体的旋转量等于时间步长乘以在`p`点计算的角速度。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   class SolidPosition : public Function<spacedim> 
 *   { 
 *   public: 
 *     static_assert(spacedim > 1, 
 *                   "Cannot instantiate SolidPosition for spacedim == 1"); 
 * 
 *     SolidPosition(const Functions::ParsedFunction<spacedim> &angular_velocity, 
 *                   const double                               time_step) 
 *       : Function<spacedim>(spacedim) 
 *       , angular_velocity(angular_velocity) 
 *       , time_step(time_step) 
 *     {} 
 * 
 *     virtual double value(const Point<spacedim> &p, 
 *                          unsigned int           component = 0) const override 
 *     { 
 *       Point<spacedim> new_position = p; 
 * 
 *       double dtheta = angular_velocity.value(p) * time_step; 
 * 
 *       new_position[0] = std::cos(dtheta) * p[0] - std::sin(dtheta) * p[1]; 
 *       new_position[1] = std::sin(dtheta) * p[0] + std::cos(dtheta) * p[1]; 
 * 
 *       return new_position[component]; 
 *     } 
 * 
 *     void set_time_step(const double new_time_step) 
 *     { 
 *       time_step = new_time_step; 
 *     } 
 * 
 *   private: 
 *     const Functions::ParsedFunction<spacedim> &angular_velocity; 
 *     double                                     time_step; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="TheStokesImmersedProblemclassdeclaration"></a> 
 * <h3>The StokesImmersedProblem class declaration</h3>
 * 

 * 
 * 我们现在准备介绍我们的教程程序的主类。像往常一样，除了构造函数外，我们只留下一个公共入口：`run()`方法。其他的都是 "私有 "的，并通过run方法本身进行访问。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim = dim> 
 *   class StokesImmersedProblem 
 *   { 
 *   public: 
 *     StokesImmersedProblem( 
 *       const StokesImmersedProblemParameters<dim, spacedim> &par); 
 * 
 *     void run(); 
 * 
 * @endcode
 * 
 * 接下来的部分包含了该类的`private`成员。第一个方法类似于前一个例子中的方法。然而，它不仅负责生成流体的网格，而且还负责生成固体的网格。第二个方法是计算最大的时间步长，保证每个粒子最多移动一个单元。这对于确保 Particles::ParticleHandler 能够找到粒子最终所在的单元是非常重要的，因为它只能从一个单元看向它的近邻（因为在并行设置中，每个MPI进程只知道它拥有的单元以及它们的近邻）。
 * 

 * 
 * 
 * @code
 *   private: 
 *     void make_grid(); 
 * 
 *     double compute_time_step() const; 
 * 
 * @endcode
 * 
 * 接下来的两个函数将初始化这个类中使用的 Particles::ParticleHandler 对象。我们有两个这样的对象。一个代表被动追踪器，用于绘制流体粒子的轨迹，而另一个代表固体的材料粒子，它们被放置在固体网格的正交点上。
 * 

 * 
 * 
 * @code
 *     void setup_tracer_particles(); 
 *     void setup_solid_particles(); 
 * 
 * @endcode
 * 
 * 剩下的设置分为两部分。以下两个函数中的第一个创建了每次模拟需要的所有对象，而另一个则设置了所有需要在每个细化步骤中重新初始化的对象。
 * 

 * 
 * 
 * @code
 *     void initial_setup(); 
 *     void setup_dofs(); 
 * 
 * @endcode
 * 
 * 装配例程与其他斯托克斯装配例程非常相似，但Nitsche限制部分除外，它利用其中一个粒子处理程序在流体域的非匹配部分进行积分，对应于固体的位置。我们将这两部分分成两个独立的函数。
 * 

 * 
 * 
 * @code
 *     void assemble_stokes_system(); 
 *     void assemble_nitsche_restriction(); 
 * 
 * @endcode
 * 
 * 其余的函数求解线性系统（看起来与 step-60 中的线性系统几乎相同），然后对解进行后处理。refine_and_transfer()方法仅在每一个`refinement_frequency`步骤中被调用，以适应网格，并确保所有在细化前的时间步骤中计算的场都正确地转移到新的网格中。这包括矢量场，以及粒子信息。同样地，我们每隔`output_frequency`步就会调用两个输出方法。
 * 

 * 
 * 
 * @code
 *     void solve(); 
 * 
 *     void refine_and_transfer(); 
 * 
 *     void output_results(const unsigned int cycle, const double time) const; 
 *     void output_particles(const Particles::ParticleHandler<spacedim> &particles, 
 *                           std::string                                 fprefix, 
 *                           const unsigned int                          iter, 
 *                           const double time) const; 
 * 
 * @endcode
 * 
 * 接下来让我们来看看这个类的成员函数。第一个是处理从参数文件中读取的运行时参数。如前所述，我们通过使其成为一个`const`引用来确保我们不能从这个类中修改这个对象。
 * 

 * 
 * 
 * @code
 *     const StokesImmersedProblemParameters<dim, spacedim> &par; 
 * 
 * @endcode
 * 
 * 然后还有MPI通信器对象，如果程序是并行运行的，我们将用它来让进程在网络上发送信息，还有`pcout`对象和定时器信息，也被 step-40 采用，例如。
 * 

 * 
 * 
 * @code
 *     MPI_Comm mpi_communicator; 
 * 
 *     ConditionalOStream pcout; 
 * 
 *     mutable TimerOutput computing_timer; 
 * 
 * @endcode
 * 
 * 接下来是关于  step-60  的主要创新点之一。这里我们假设固体和流体都是完全分布的三角形。这使得问题可以扩展到非常大的自由度，代价是要沟通所有非匹配三角形之间的重叠区域。这一点特别棘手，因为我们没有对两个三角形的各个子域的相对位置或分布做出假设。特别是，我们假设每个进程只拥有 "solid_tria "的一部分，以及 "fluid_tria "的一部分，不一定在同一个物理区域，也不一定重叠。
 * 

 * 
 * 我们原则上可以尝试创建初始分区，使每个过程的子域在固体和流体区域之间重叠。然而，这种重叠在模拟过程中会被破坏，我们将不得不一次又一次地重新分配DoF。我们在本教程中采用的方法更加灵活，而且成本也不高。我们在模拟开始时进行两次全对全的通信，以交换每个处理器的几何占用信息（近似）（通过包围盒的集合完成）。
 * 

 * 
 * 这个信息被 Particles::ParticleHandler 类用来交换（使用某对某的通信模式）所有的粒子，因此每个进程都知道生活在它所拥有的流体子域所占区域的粒子。
 * 

 * 
 * 为了把重叠的区域连接起来，我们利用了ParticleHandler类中实现的设施。
 * 

 * 
 * 
 * @code
 *     parallel::distributed::Triangulation<spacedim>      fluid_tria; 
 *     parallel::distributed::Triangulation<dim, spacedim> solid_tria; 
 * 
 * @endcode
 * 
 * 接下来是对所使用的有限元的描述，以及适当的正交公式和相应的DoFHandler对象。在目前的实现中，只有`fluid_fe`是真正必要的。为了完整起见，并便于扩展，我们还保留了`solid_fe`，但它被初始化为一个FE_Nothing有限元空间，即没有自由度的空间。
 * 

 * 
 * 我们将这两个有限元空间声明为 `std::unique_ptr` 对象，而不是普通的成员变量，以便在`StokesImmersedProblemParameters'被初始化后生成它们。特别是，它们将在`initial_setup()`方法中被初始化。
 * 

 * 
 * 
 * @code
 *     std::unique_ptr<FiniteElement<spacedim>>      fluid_fe; 
 *     std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe; 
 * 
 *     std::unique_ptr<Quadrature<spacedim>> fluid_quadrature_formula; 
 *     std::unique_ptr<Quadrature<dim>>      solid_quadrature_formula; 
 * 
 *     DoFHandler<spacedim>      fluid_dh; 
 *     DoFHandler<dim, spacedim> solid_dh; 
 * 
 *     std::unique_ptr<MappingFEField<dim, spacedim>> solid_mapping; 
 * 
 * @endcode
 * 
 * 与 step-22 中的做法类似，我们使用一个块状系统来处理问题的斯托克斯部分，并非常密切地遵循那里的做法。
 * 

 * 
 * 
 * @code
 *     std::vector<IndexSet> fluid_owned_dofs; 
 *     std::vector<IndexSet> solid_owned_dofs; 
 * 
 *     std::vector<IndexSet> fluid_relevant_dofs; 
 *     std::vector<IndexSet> solid_relevant_dofs; 
 * 
 * @endcode
 * 
 * 利用这种自由度的划分，我们就可以定义所有必要的对象来描述有关的线性系统。
 * 

 * 
 * 
 * @code
 *     AffineConstraints<double> constraints; 
 * 
 *     LA::MPI::BlockSparseMatrix system_matrix; 
 *     LA::MPI::BlockSparseMatrix preconditioner_matrix; 
 * 
 *     LA::MPI::BlockVector solution; 
 *     LA::MPI::BlockVector locally_relevant_solution; 
 *     LA::MPI::BlockVector system_rhs; 
 * 
 * @endcode
 * 
 * 让我们转到这个程序的粒子方面。有两个 Particles::ParticleHandler 对象用于耦合固体和流体，以及描述被动追踪器。在许多方面，这些对象的作用类似于离散化中使用的DoFHandler类，也就是说，它们提供了粒子的枚举，并允许查询每个粒子的信息。
 * 

 * 
 * 
 * @code
 *     Particles::ParticleHandler<spacedim> tracer_particle_handler; 
 *     Particles::ParticleHandler<spacedim> solid_particle_handler; 
 * 
 * @endcode
 * 
 * 对于每个追踪器粒子，我们需要计算其当前位置的速度场，并使用离散时间步进方案更新其位置。我们使用分布式线性代数对象来做这件事，这些对象存储了每个粒子的位置或速度的坐标。也就是说，这些向量有`tracer_particle_handler.n_global_particles() * spacedim`项，我们将以一种方式来存储这些向量的一部分，以便在所有进程中进行划分。(隐含地，我们在此假设每个粒子的`spacedim'坐标被存储在向量的连续条目中)。因此，我们需要确定每个向量条目的所有者是谁。我们将这个所有者设定为等于在时间 $t=0$ 产生该粒子的进程。这个信息对每一个进程都存储在`locally_owned_tracer_particle_coordinates`索引集里。
 * 

 * 
 * 一旦粒子被分配到与拥有粒子所在区域的进程相匹配，我们将需要从该进程读取相应的速度场。我们通过填充一个只读的速度矢量场来实现这一目标，该矢量场包含了Ghost条目中的相关信息。这是通过`locally_relevant_tracer_particle_coordinates`索引集实现的，该索引集记录了模拟过程中的变化情况，也就是说，它记录了当前进程拥有的粒子最终出现在哪里，以及谁拥有最终出现在我的子域的粒子。
 * 

 * 
 * 虽然这不是最有效的策略，但我们保持这种方式是为了说明事情在真实的流固耦合（FSI）问题中是如何运作的。如果一个粒子与一个特定的固体自由度相联系，我们就不能自由选择谁拥有它，我们必须把这个信息传达给周围的人。我们在这里说明了这一点，并表明通信模式是点对点的，就算法的总成本而言可以忽略不计。
 * 

 * 
 * 然后，基于这些细分定义的向量被用来存储粒子的速度（只读，有幽灵条目）和它们的位移（读/写，没有幽灵条目）。
 * 

 * 
 * 
 * @code
 *     IndexSet locally_owned_tracer_particle_coordinates; 
 *     IndexSet locally_relevant_tracer_particle_coordinates; 
 * 
 *     LA::MPI::Vector tracer_particle_velocities; 
 *     LA::MPI::Vector relevant_tracer_particle_displacements; 
 * 
 * @endcode
 * 
 * 本教程程序的关键点之一是两个独立的 parallel::distributed::Triangulation 对象之间的耦合，其中一个对象可能相对于另一个对象移动和变形（可能有较大的变形）。当流体和实体的三角形都是 parallel::distributed::Triangulation, 类型时，每个进程只能访问这两个三角形中每个单元的局部拥有的部分。如上所述，一般情况下，本地拥有的域是不重叠的。
 * 

 * 
 * 为了允许在不重叠的 parallel::distributed::Triangulation 对象之间有效地交换信息，该库的一些算法要求用户提供三角形的本地拥有部分所占区域的粗略描述，其形式是每个进程的轴对齐的边界盒集合，这些边界盒提供了域的本地拥有部分的完整覆盖。这种信息就可以用于这样的情况：人们需要向已知位置周围的单元格的所有者发送信息，而不知道这个所有者实际上是谁。但是，如果我们知道每个进程拥有的几何区域或体积的边界盒集合，那么我们就可以确定可能拥有该位置所在单元的所有进程的一个子集：即其边界盒包含该点的所有进程。与其向所有进程发送与该位置相关的信息，不如只向具有点对点通信基元的一小部分进程发送信息。你会注意到，这也允许典型的时间与内存的权衡：我们愿意存储的关于每个进程拥有的区域的数据越多--以更精细的边界框信息的形式--我们必须执行的通信就越少）。
 * 

 * 
 * 我们通过收集一个向量（长度为 Utilities::MPI::n_mpi_processes()) 的BoundingBox对象的向量）来构建这些信息。我们用extract_rtree_level()函数填充这个向量，并允许用户选择要提取的树的哪一级。这个 "级别 "对应的是与边界框重叠的区域应该有多粗/多细。
 * 

 * 
 * 作为一个例子，这是由extract_rtree_level()函数对一个分布在三个过程中的二维超球所提取的结果。每张图片中，绿色显示的是与每个进程上的三角形的本地所有单元相关的边界框，紫色显示的是从rtree中提取的边界框。
 * 

 * 
 * @image html rtree-process-0.png  
 * @image html rtree-process-1.png  
 * @image html rtree-process-2.png  
 * 

 * 
 * 我们将这些盒子存储在一个全局成员变量中，在每个细化步骤中都会更新。
 * 

 * 
 * 
 * @code
 *     std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="TheStokesImmersedProblemclassimplementation"></a> 
 * <h3>The StokesImmersedProblem class implementation</h3>
 * 
 * <a name="Objectconstructionandmeshinitializationfunctions"></a> 
 * <h4>Object construction and mesh initialization functions</h4>
 * 

 * 
 * 在构造函数中，我们创建了mpi_communicator，以及流体和实体的三角计算和dof_handler。通过使用mpi_communicator，我们构建了ConditionalOStream和TimerOutput对象。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   StokesImmersedProblem<dim, spacedim>::StokesImmersedProblem( 
 *     const StokesImmersedProblemParameters<dim, spacedim> &par) 
 *     : par(par) 
 *     , mpi_communicator(MPI_COMM_WORLD) 
 *     , pcout(std::cout, 
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
 *     , computing_timer(mpi_communicator, 
 *                       pcout, 
 *                       TimerOutput::summary, 
 *                       TimerOutput::wall_times) 
 *     , fluid_tria(mpi_communicator, 
 *                  typename Triangulation<spacedim>::MeshSmoothing( 
 *                    Triangulation<spacedim>::smoothing_on_refinement | 
 *                    Triangulation<spacedim>::smoothing_on_coarsening)) 
 *     , solid_tria(mpi_communicator, 
 *                  typename Triangulation<dim, spacedim>::MeshSmoothing( 
 *                    Triangulation<dim, spacedim>::smoothing_on_refinement | 
 *                    Triangulation<dim, spacedim>::smoothing_on_coarsening)) 
 *     , fluid_dh(fluid_tria) 
 *     , solid_dh(solid_tria) 
 *   {} 
 * 
 * @endcode
 * 
 * 为了生成网格，我们首先尝试使用deal.II GridGenerator命名空间中的函数，通过利用 GridGenerator::generate_from_name_and_argument().  如果这个函数失败，那么我们使用以下方法，名称被解释为文件名，参数被解释为从流形ID到CAD文件的映射，并使用OpenCASCADE命名空间设施转换为流形描述符。在顶部，我们把文件读成一个三角图。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void read_grid_and_cad_files(const std::string &grid_file_name, 
 *                                const std::string &ids_and_cad_file_names, 
 *                                Triangulation<dim, spacedim> &tria) 
 *   { 
 *     GridIn<dim, spacedim> grid_in; 
 *     grid_in.attach_triangulation(tria); 
 *     grid_in.read(grid_file_name); 
 * 
 * @endcode
 * 
 * 如果我们走到这一步，那么三角图已经被读取，我们已经准备好将正确的流形描述附加到它上面。只有在deal.II支持OpenCASCADE的情况下，我们才会执行接下来的几行代码。对于地图中的每个条目，我们尝试打开相应的CAD文件，分析它，并根据其内容，选择一个 OpenCASCADE::ArcLengthProjectionLineManifold （如果CAD文件包含一个`TopoDS_Edge'或一个`TopoDS_Wire'）或一个 OpenCASCADE::NURBSPatchManifold, ，如果文件包含一个面。请注意，如果CAD文件不包含单一的线、边或面，在生成Manifold时将会抛出一个断言。
 * 

 * 
 * 我们使用 Patterns::Tools::Convert 类来完成从字符串到歧管ID和文件名之间的映射的转换。
 * 

 * 
 * 
 * @code
 * #ifdef DEAL_II_WITH_OPENCASCADE 
 *     using map_type  = std::map<types::manifold_id, std::string>; 
 *     using Converter = Patterns::Tools::Convert<map_type>; 
 * 
 *     for (const auto &pair : Converter::to_value(ids_and_cad_file_names)) 
 *       { 
 *         const auto &manifold_id   = pair.first; 
 *         const auto &cad_file_name = pair.second; 
 * 
 *         const auto extension = boost::algorithm::to_lower_copy( 
 *           cad_file_name.substr(cad_file_name.find_last_of('.') + 1)); 
 * 
 *         TopoDS_Shape shape; 
 *         if (extension == "iges" || extension == "igs") 
 *           shape = OpenCASCADE::read_IGES(cad_file_name); 
 *         else if (extension == "step" || extension == "stp") 
 *           shape = OpenCASCADE::read_STEP(cad_file_name); 
 *         else 
 *           AssertThrow(false, 
 *                       ExcNotImplemented("We found an extension that we " 
 *                                         "do not recognize as a CAD file " 
 *                                         "extension. Bailing out.")); 
 * 
 * @endcode
 * 
 * 现在我们检查一下这个 "形状 "中包含了多少个面。OpenCASCADE本质上是三维的，所以如果这个数字是零，我们就把它解释为线状流形，否则就解释为`spacedim`=3中的 OpenCASCADE::NormalToMeshProjectionManifold ，或者`spacedim`=2中的 OpenCASCADE::NURBSPatchManifold 。
 * 

 * 
 * 
 * @code
 *         const auto n_elements = OpenCASCADE::count_elements(shape); 
 *         if ((std::get<0>(n_elements) == 0)) 
 *           tria.set_manifold( 
 *             manifold_id, 
 *             OpenCASCADE::ArclengthProjectionLineManifold<dim, spacedim>(shape)); 
 *         else if (spacedim == 3) 
 *           { 
 * 
 * @endcode
 * 
 * 我们使用这个技巧，因为 OpenCASCADE::NormalToMeshProjectionManifold 只在spacedim = 3的情况下实现。上面的检查保证了事情的实际运作是正确的。
 * 

 * 
 * 
 * @code
 *             const auto t = reinterpret_cast<Triangulation<dim, 3> *>(&tria); 
 *             t->set_manifold(manifold_id, 
 *                             OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>( 
 *                               shape)); 
 *           } 
 *         else 
 * 
 * @endcode
 * 
 * 我们也允许基于单个NURBS补丁的二维空间的曲面描述。要做到这一点，CAD文件必须包含一个单一的`TopoDS_Face`。
 * 

 * 
 * 
 * @code
 *           tria.set_manifold(manifold_id, 
 *                             OpenCASCADE::NURBSPatchManifold<dim, spacedim>( 
 *                               TopoDS::Face(shape))); 
 *       } 
 * #else 
 *     (void)ids_and_cad_file_names; 
 *     AssertThrow(false, ExcNotImplemented("Generation of the grid failed.")); 
 * #endif 
 *   } 
 * 
 * @endcode
 * 
 * 现在让我们把东西放在一起，并制作所有必要的网格。如上所述，我们首先尝试在内部生成网格，如果我们失败了（即如果我们最终进入了`catch'子句），那么我们就继续执行上述函数。
 * 

 * 
 * 我们对流体和固体网格都重复这个模式。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::make_grid() 
 *   { 
 *     try 
 *       { 
 *         GridGenerator::generate_from_name_and_arguments( 
 *           fluid_tria, par.name_of_fluid_grid, par.arguments_for_fluid_grid); 
 *       } 
 *     catch (...) 
 *       { 
 *         pcout << "Generating from name and argument failed." << std::endl 
 *               << "Trying to read from file name." << std::endl; 
 *         read_grid_and_cad_files(par.name_of_fluid_grid, 
 *                                 par.arguments_for_fluid_grid, 
 *                                 fluid_tria); 
 *       } 
 *     fluid_tria.refine_global(par.initial_fluid_refinement); 
 * 
 *     try 
 *       { 
 *         GridGenerator::generate_from_name_and_arguments( 
 *           solid_tria, par.name_of_solid_grid, par.arguments_for_solid_grid); 
 *       } 
 *     catch (...) 
 *       { 
 *         read_grid_and_cad_files(par.name_of_solid_grid, 
 *                                 par.arguments_for_solid_grid, 
 *                                 solid_tria); 
 *       } 
 * 
 *     solid_tria.refine_global(par.initial_solid_refinement); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Particleinitializationfunctions"></a> 
 * <h4>Particle initialization functions</h4>
 * 

 * 
 * 一旦固体和流体网格被创建，我们就开始填充 Particles::ParticleHandler 对象。我们要处理的第一个对象是用来跟踪流体中的被动追踪器的对象。这些东西只是沿途传送，从某种意义上说，它们的位置并不重要：我们只是想用它们来观察流体被传送的位置。我们可以使用任何我们选择的方式来确定它们的初始位置。一个方便的方法是将初始位置创建为我们所选择的形状的网格顶点，这个选择由参数文件中的一个运行时参数决定。
 * 

 * 
 * 在这个实现中，我们使用FE_Q有限元空间的支持点来创建追踪器，这些支持点定义在一个临时网格上，然后被丢弃。在这个网格中，我们只保留与支撑点相关的 Particles::Particle 对象（存储在 Particles::ParticleHandler 类中）。
 * 

 * 
 * Particles::ParticleHandler 类提供了插入一组粒子的可能性，这些粒子实际生活在活动过程所拥有的域的一部分。然而，在这种情况下，这个功能是不够的。作为任意网格（与流体网格不匹配）上的FE_Q对象的本地拥有的支持点所产生的粒子没有理由位于流体网格的本地拥有的子域的同一物理区域内。事实上，这种情况几乎不会发生，尤其是我们想要跟踪粒子本身发生了什么。
 * 

 * 
 * 在粒子入室方法（PIC）中，人们通常习惯于将粒子的所有权分配给粒子所在的过程。在本教程中，我们说明了一种不同的方法，如果想跟踪与粒子有关的信息，这种方法是很有用的（例如，如果一个粒子与一个特定的自由度有关，而这个自由度是由一个特定的过程所拥有的，不一定是在任何特定时间拥有该粒子所在的流体单元的同一个过程）。在这里使用的方法中，粒子的所有权在开始时被分配一次，每当原始所有者需要从拥有粒子所在单元的进程中获得信息时，就会发生一对一的通信。我们确保使用初始粒子分布来设置粒子的所有权，并在程序的整个执行过程中保持相同的所有权。
 * 

 * 
 * 有了这个概述，让我们看看这个函数做什么。在顶部，我们创建了一个临时的三角形和DoFHandler对象，我们将从中获取初始粒子位置的节点位置。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::setup_tracer_particles() 
 *   { 
 *     parallel::distributed::Triangulation<spacedim> particle_insert_tria( 
 *       mpi_communicator); 
 *     GridGenerator::generate_from_name_and_arguments( 
 *       particle_insert_tria, 
 *       par.name_of_particle_grid, 
 *       par.arguments_for_particle_grid); 
 *     particle_insert_tria.refine_global(par.particle_insertion_refinement); 
 * 
 *     FE_Q<spacedim>       particles_fe(1); 
 *     DoFHandler<spacedim> particles_dof_handler(particle_insert_tria); 
 *     particles_dof_handler.distribute_dofs(particles_fe); 
 * 
 * @endcode
 * 
 * 这就是事情开始变得复杂的地方。由于我们可能会在并行环境中运行这个程序，每个并行进程现在都会创建这些临时三角形和DoFHandlers。但是，在完全分布式三角形中，活动进程只知道本地拥有的单元，而不知道其他进程是如何分布自己的单元的。这对于上面创建的临时三角形以及我们想嵌入粒子的流体三角形都是如此。另一方面，一般来说，这两个三角形的局部已知部分不会重合。也就是说，我们将从临时网格的节点位置创建的粒子的位置是任意的，并且可能落在当前进程无法访问的流体三角结构的区域内（即流体领域中细胞是人工的区域）。为了了解将这些粒子发送给谁，我们需要对流体网格在处理器中的分布有一个（粗略的）概念。
 * 

 * 
 * 我们通过以下方式来构建这一信息：首先建立一个以本地拥有的单元为边界的盒子的索引树，然后提取该树的第一层中的一个。
 * 

 * 
 * 
 * @code
 *     std::vector<BoundingBox<spacedim>> all_boxes; 
 *     all_boxes.reserve(fluid_tria.n_locally_owned_active_cells()); 
 *     for (const auto &cell : fluid_tria.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         all_boxes.emplace_back(cell->bounding_box()); 
 * 
 *     const auto tree = pack_rtree(all_boxes); 
 *     const auto local_boxes = 
 *       extract_rtree_level(tree, par.fluid_rtree_extraction_level); 
 * 
 * @endcode
 * 
 * 每个进程现在都有一个完全包围所有本地拥有的进程的边界盒集合（但可能与其他进程的边界盒相重叠）。然后我们在所有参与的进程之间交换这些信息，这样每个进程都知道所有其他进程的边界盒。
 * 

 * 
 * 有了这些信息，我们就可以将`tracer_particle_handler`初始化到流体网格，并从（临时）tracer particles triangulation的支持点生成粒子。这个函数调用使用了我们刚刚构建的`global_bounding_boxes`对象，以确定将位置来自`particles_dof_handler`的本地拥有部分的粒子发送到何处。在这个调用结束时，每个粒子将被分配到正确的进程（即拥有粒子所在的流体单元的进程）。在这一点上，我们也将他们的编号输出到屏幕上。
 * 

 * 
 * 
 * @code
 *     global_fluid_bounding_boxes = 
 *       Utilities::MPI::all_gather(mpi_communicator, local_boxes); 
 * 
 *     tracer_particle_handler.initialize(fluid_tria, 
 *                                        StaticMappingQ1<spacedim>::mapping); 
 * 
 *     Particles::Generators::dof_support_points(particles_dof_handler, 
 *                                               global_fluid_bounding_boxes, 
 *                                               tracer_particle_handler); 
 * 
 *     pcout << "Tracer particles: " 
 *           << tracer_particle_handler.n_global_particles() << std::endl; 
 * 
 * @endcode
 * 
 * 这样创建的每个粒子都有一个唯一的ID。在下面的算法中的某个时刻，我们将需要包含每个粒子的位置和速度信息的向量。这个向量的大小为`n_particles * // spacedim`，我们需要为每个粒子提供位置和速度信息。
 * spacedim`，我们将不得不以一种方式来存储这个向量的元素，以便每个并行进程 "拥有 "与它拥有的粒子的坐标相对应的那些元素。换句话说，我们必须在所有进程中划分0和`n_particles * spacedim`之间的索引空间。我们可以通过查询`tracer_particle_handler`的本地相关粒子的ID来做到这一点，并构建需要的索引，将所有粒子的位置和速度存储在一个（平行分布的）矢量中，其中我们隐含地假设我们将每个位置或速度的坐标存储在`spacedim`连续的矢量元素中（这就是 IndexSet::tensor_priduct() 函数的作用）。
 * 

 * 
 * 
 * @code
 *     locally_owned_tracer_particle_coordinates = 
 *       tracer_particle_handler.locally_owned_particle_ids().tensor_product( 
 *         complete_index_set(spacedim)); 
 * 
 * @endcode
 * 
 * 在模拟开始时，所有粒子都在它们的原始位置。当粒子移动时，它们可能会穿越到另一个进程所拥有的领域的某个部分。如果发生这种情况，当前进程会正式保持对粒子的 "所有权"，但可能需要从粒子落地的进程中读取访问。我们将这一信息保存在另一个索引集中，该索引集存储了当前进程的子域中的所有粒子的索引，不管它们是否一直在这里。
 * 

 * 
 * 保留这个索引集使我们能够利用线性代数类来进行有关粒子位置和速度的所有通信。这模拟了在固体域中解决另一个问题的情况下会发生的情况（如在流体-结构相互作用中。在后一种情况下，实体域上的额外DOFs将被耦合到流体域中发生的情况。
 * 

 * 
 * 
 * @code
 *     locally_relevant_tracer_particle_coordinates = 
 *       locally_owned_tracer_particle_coordinates; 
 * 
 * @endcode
 * 
 * 最后，我们要确保在细化时，粒子被正确转移。在进行局部细化或粗化时，粒子会落在另一个单元中。原则上，我们可以在细化后重新分配所有的粒子，但是这将是非常昂贵的。
 * 

 * 
 * Particles::ParticleHandler 类有一种方法可以在细化时将信息从一个单元转移到它的子单元或它的父单元，而不需要重构整个数据结构。这是通过向三角结构注册两个回调函数来实现的。这些函数将在细化即将发生和刚刚发生时收到一个信号，并将以最小的计算成本将所有信息转移到新的细化网格中。
 * 

 * 
 * 
 * @code
 *     fluid_tria.signals.pre_distributed_refinement.connect( 
 *       [&]() { tracer_particle_handler.register_store_callback_function(); }); 
 * 
 *     fluid_tria.signals.post_distributed_refinement.connect([&]() { 
 *       tracer_particle_handler.register_load_callback_function(false); 
 *     }); 
 *   } 
 * 
 * @endcode
 * 
 * 与我们对被动追踪器所做的类似，我们接下来设置追踪实体网格的正交点的粒子。这里的主要区别是，我们还想给每个粒子附加一个权重值（正交点的 "JxW "值），这样我们就可以在不直接访问原始实体网格的情况下计算积分。
 * 

 * 
 * 这是通过利用 Particles::Particle 类的 "属性 "概念实现的。它可以（以一种有效的内存方式）在一个 Particles::ParticleHandler 对象内为每个 Particles::Particle 对象存储任意数量的`双`数字。我们利用这种可能性来存储实体网格的正交点的JxW值。
 * 

 * 
 * 在我们的例子中，我们只需要为每个粒子存储一个属性：实体网格上的积分的JxW值。这将在构造时作为最后一个参数传递给solid_particle_handler对象。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::setup_solid_particles() 
 *   { 
 *     QGauss<dim> quadrature(fluid_fe->degree + 1); 
 * 
 *     const unsigned int n_properties = 1; 
 *     solid_particle_handler.initialize(fluid_tria, 
 *                                       StaticMappingQ1<spacedim>::mapping, 
 *                                       n_properties); 
 * 
 * @endcode
 * 
 * 我们在本地生成的粒子数等于本地拥有的单元总数乘以每个单元中使用的正交点的数量。我们将所有这些点存储在一个向量中，并将其相应的属性存储在一个向量的向量中。
 * 

 * 
 * 
 * @code
 *     std::vector<Point<spacedim>> quadrature_points_vec; 
 *     quadrature_points_vec.reserve(quadrature.size() * 
 *                                   solid_tria.n_locally_owned_active_cells()); 
 * 
 *     std::vector<std::vector<double>> properties; 
 *     properties.reserve(quadrature.size() * 
 *                        solid_tria.n_locally_owned_active_cells()); 
 * 
 *     FEValues<dim, spacedim> fe_v(*solid_fe, 
 *                                  quadrature, 
 *                                  update_JxW_values | update_quadrature_points); 
 *     for (const auto &cell : solid_dh.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           fe_v.reinit(cell); 
 *           const auto &points = fe_v.get_quadrature_points(); 
 *           const auto &JxW    = fe_v.get_JxW_values(); 
 * 
 *           for (unsigned int q = 0; q < points.size(); ++q) 
 *             { 
 *               quadrature_points_vec.emplace_back(points[q]); 
 *               properties.emplace_back( 
 *                 std::vector<double>(n_properties, JxW[q])); 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 我们以处理示踪粒子的同样方式进行，重新使用计算出的边界盒。然而，我们首先检查`global_fluid_bounding_boxes`对象是否已经被填充。这里当然应该是这样的，因为这个方法是在初始化示踪粒子的方法之后调用的。然而，我们要确保，如果将来有人决定（无论出于什么原因）先初始化固体粒子处理程序，或者只复制教程的这一部分，当事情没有按照预期进行时，会抛出一个有意义的异常。
 * 

 * 
 * 由于我们已经存储了正交点的位置，我们可以使用这些位置来直接使用`solid_particle_handler`插入粒子，而不必通过 Particles::Generators 函数。
 * 

 * 
 * 
 * @code
 *     Assert(!global_fluid_bounding_boxes.empty(), 
 *            ExcInternalError( 
 *              "I was expecting the " 
 *              "global_fluid_bounding_boxes to be filled at this stage. " 
 *              "Make sure you fill this vector before trying to use it " 
 *              "here. Bailing out.")); 
 * 
 *     solid_particle_handler.insert_global_particles(quadrature_points_vec, 
 *                                                    global_fluid_bounding_boxes, 
 *                                                    properties); 
 * 
 * @endcode
 * 
 * 和前面的函数一样，我们最后要确保在细化时，粒子被正确转移。
 * 

 * 
 * 
 * @code
 *     fluid_tria.signals.pre_distributed_refinement.connect( 
 *       [&]() { solid_particle_handler.register_store_callback_function(); }); 
 * 
 *     fluid_tria.signals.post_distributed_refinement.connect( 
 *       [&]() { solid_particle_handler.register_load_callback_function(false); }); 
 * 
 *     pcout << "Solid particles: " << solid_particle_handler.n_global_particles() 
 *           << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="DoFinitializationfunctions"></a> 
 * <h4>DoF initialization functions</h4>
 * 

 * 
 * 我们设置了有限元空间和整个步骤中使用的正交公式。对于流体，我们使用Taylor-Hood元素（例如 $Q_k \times Q_{k-1}$  ）。由于我们没有解决固体领域的任何方程，所以产生了一个空的有限元空间。这个程序的一个自然扩展是解决流体结构的相互作用问题，这就要求`solid_fe`使用更有用的FiniteElement类。
 * 

 * 
 * 和其他许多函数一样，我们在这里存储了进行操作所需的时间。当前的函数把它的时间信息放到一个标签为 "初始设置 "的部分。在不同的函数中对这个定时器进行了许多其他的调用。它们允许监测每个单独函数的绝对和相对成本，以确定瓶颈。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::initial_setup() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Initial setup"); 
 * 
 *     fluid_fe = 
 *       std::make_unique<FESystem<spacedim>>(FE_Q<spacedim>(par.velocity_degree), 
 *                                            spacedim, 
 *                                            FE_Q<spacedim>(par.velocity_degree - 
 *                                                           1), 
 *                                            1); 
 * 
 *     solid_fe = std::make_unique<FE_Nothing<dim, spacedim>>(); 
 *     solid_dh.distribute_dofs(*solid_fe); 
 * 
 *     fluid_quadrature_formula = 
 *       std::make_unique<QGauss<spacedim>>(par.velocity_degree + 1); 
 *     solid_quadrature_formula = 
 *       std::make_unique<QGauss<dim>>(par.velocity_degree + 1); 
 *   } 
 * 
 * @endcode
 * 
 * 接下来我们构建分布式块状矩阵和向量，用于解决问题中出现的线性方程。这个函数改编自 step-55 ，我们参考这个步骤进行全面解释。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::setup_dofs() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Setup dofs"); 
 * 
 *     fluid_dh.distribute_dofs(*fluid_fe); 
 * 
 *     std::vector<unsigned int> stokes_sub_blocks(spacedim + 1, 0); 
 *     stokes_sub_blocks[spacedim] = 1; 
 *     DoFRenumbering::component_wise(fluid_dh, stokes_sub_blocks); 
 * 
 *     auto dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(fluid_dh, stokes_sub_blocks); 
 * 
 *     const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1]; 
 * 
 *     pcout << "   Number of degrees of freedom: " << fluid_dh.n_dofs() << " (" 
 *           << n_u << '+' << n_p << " -- " 
 *           << solid_particle_handler.n_global_particles() << '+' 
 *           << tracer_particle_handler.n_global_particles() << ')' << std::endl; 
 * 
 *     fluid_owned_dofs.resize(2); 
 *     fluid_owned_dofs[0] = fluid_dh.locally_owned_dofs().get_view(0, n_u); 
 *     fluid_owned_dofs[1] = 
 *       fluid_dh.locally_owned_dofs().get_view(n_u, n_u + n_p); 
 * 
 *     IndexSet locally_relevant_dofs; 
 *     DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs); 
 *     fluid_relevant_dofs.resize(2); 
 *     fluid_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u); 
 *     fluid_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p); 
 * 
 *     { 
 *       constraints.reinit(locally_relevant_dofs); 
 * 
 *       FEValuesExtractors::Vector velocities(0); 
 *       DoFTools::make_hanging_node_constraints(fluid_dh, constraints); 
 *       VectorTools::interpolate_boundary_values( 
 *         fluid_dh, 
 *         0, 
 *         Functions::ZeroFunction<spacedim>(spacedim + 1), 
 *         constraints, 
 *         fluid_fe->component_mask(velocities)); 
 *       constraints.close(); 
 *     } 
 * 
 *     auto locally_owned_dofs_per_processor = 
 *       Utilities::MPI::all_gather(mpi_communicator, 
 *                                  fluid_dh.locally_owned_dofs()); 
 *     { 
 *       system_matrix.clear(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1); 
 *       for (unsigned int c = 0; c < spacedim + 1; ++c) 
 *         for (unsigned int d = 0; d < spacedim + 1; ++d) 
 *           if (c == spacedim && d == spacedim) 
 *             coupling[c][d] = DoFTools::none; 
 *           else if (c == spacedim || d == spacedim || c == d) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         fluid_dh, coupling, dsp, constraints, false); 
 * 
 *       SparsityTools::distribute_sparsity_pattern( 
 *         dsp, 
 *         locally_owned_dofs_per_processor, 
 *         mpi_communicator, 
 *         locally_relevant_dofs); 
 * 
 *       system_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator); 
 *     } 
 * 
 *     { 
 *       preconditioner_matrix.clear(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1); 
 *       for (unsigned int c = 0; c < spacedim + 1; ++c) 
 *         for (unsigned int d = 0; d < spacedim + 1; ++d) 
 *           if (c == spacedim && d == spacedim) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         fluid_dh, coupling, dsp, constraints, false); 
 *       SparsityTools::distribute_sparsity_pattern( 
 *         dsp, 
 *         locally_owned_dofs_per_processor, 
 *         mpi_communicator, 
 *         locally_relevant_dofs); 
 *       preconditioner_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator); 
 *     } 
 * 
 *     locally_relevant_solution.reinit(fluid_owned_dofs, 
 *                                      fluid_relevant_dofs, 
 *                                      mpi_communicator); 
 *     system_rhs.reinit(fluid_owned_dofs, mpi_communicator); 
 *     solution.reinit(fluid_owned_dofs, mpi_communicator); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Assemblyfunctions"></a> 
 * <h4>Assembly functions</h4>
 * 

 * 
 * 我们将系统矩阵、预处理矩阵和右手边组合起来。这段代码改编自 step-55 ，基本上是 step-27 的内容，如果你知道斯托克斯方程是什么样子的，就会觉得很标准。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::assemble_stokes_system() 
 *   { 
 *     system_matrix         = 0; 
 *     preconditioner_matrix = 0; 
 *     system_rhs            = 0; 
 * 
 *     TimerOutput::Scope t(computing_timer, "Assemble Stokes terms"); 
 * 
 *     FEValues<spacedim> fe_values(*fluid_fe, 
 *                                  *fluid_quadrature_formula, 
 *                                  update_values | update_gradients | 
 *                                    update_quadrature_points | 
 *                                    update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = fluid_quadrature_formula->size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<Vector<double>> rhs_values(n_q_points, 
 *                                            Vector<double>(spacedim + 1)); 
 * 
 *     std::vector<Tensor<2, spacedim>> grad_phi_u(dofs_per_cell); 
 *     std::vector<double>              div_phi_u(dofs_per_cell); 
 *     std::vector<double>              phi_p(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 *     const FEValuesExtractors::Vector     velocities(0); 
 *     const FEValuesExtractors::Scalar     pressure(spacedim); 
 * 
 *     for (const auto &cell : fluid_dh.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           cell_matrix  = 0; 
 *           cell_matrix2 = 0; 
 *           cell_rhs     = 0; 
 * 
 *           fe_values.reinit(cell); 
 *           par.rhs.vector_value_list(fe_values.get_quadrature_points(), 
 *                                     rhs_values); 
 *           for (unsigned int q = 0; q < n_q_points; ++q) 
 *             { 
 *               for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                 { 
 *                   grad_phi_u[k] = fe_values[velocities].gradient(k, q); 
 *                   div_phi_u[k]  = fe_values[velocities].divergence(k, q); 
 *                   phi_p[k]      = fe_values[pressure].value(k, q); 
 *                 } 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     { 
 *                       cell_matrix(i, j) += 
 *                         (par.viscosity * 
 *                            scalar_product(grad_phi_u[i], grad_phi_u[j]) - 
 *                          div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
 *                         fe_values.JxW(q); 
 * 
 *                       cell_matrix2(i, j) += 1.0 / par.viscosity * phi_p[i] * 
 *                                             phi_p[j] * fe_values.JxW(q); 
 *                     } 
 * 
 *                   const unsigned int component_i = 
 *                     fluid_fe->system_to_component_index(i).first; 
 *                   cell_rhs(i) += fe_values.shape_value(i, q) * 
 *                                  rhs_values[q](component_i) * fe_values.JxW(q); 
 *                 } 
 *             } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           constraints.distribute_local_to_global(cell_matrix, 
 *                                                  cell_rhs, 
 *                                                  local_dof_indices, 
 *                                                  system_matrix, 
 *  
 * 
 *           constraints.distribute_local_to_global(cell_matrix2, 
 *                                                  local_dof_indices, 
 *                                                  preconditioner_matrix); 
 *         } 
 * 
 *     system_matrix.compress(VectorOperation::add); 
 *     preconditioner_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 *   } 
 * 
 * @endcode
 * 
 * 下面的方法是处理因对叶轮施加速度而产生的惩罚项。从某种意义上说，它是本教程的核心，但它相对简单。这里我们利用`solid_particle_handler`来计算Nitsche限制或嵌入域中的惩罚。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::assemble_nitsche_restriction() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Assemble Nitsche terms"); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(spacedim); 
 * 
 *     SolidVelocity<spacedim> solid_velocity(par.angular_velocity); 
 * 
 *     std::vector<types::global_dof_index> fluid_dof_indices( 
 *       fluid_fe->n_dofs_per_cell()); 
 * 
 *     FullMatrix<double>     local_matrix(fluid_fe->n_dofs_per_cell(), 
 *                                     fluid_fe->n_dofs_per_cell()); 
 *     dealii::Vector<double> local_rhs(fluid_fe->n_dofs_per_cell()); 
 * 
 *     const auto penalty_parameter = 
 *       1.0 / GridTools::minimal_cell_diameter(fluid_tria); 
 * 
 * @endcode
 * 
 * 我们在所有的本地粒子上循环。虽然这可以直接通过循环所有的单元格来实现，但这将迫使我们循环许多不包含粒子的单元格。因此，我们在所有的粒子上循环，但是，我们得到粒子所在的单元格的参考，然后在该单元格中循环所有的粒子。这使得我们能够跳过不包含粒子的单元格，但又能集合每个单元格的局部矩阵和rhs来应用Nitsche的限制。一旦我们完成了一个单元格上的所有粒子，我们就将`粒子`迭代器推进到当前单元格上的粒子的末端（这是`while`循环体的最后一行）。
 * 

 * 
 * 
 * @code
 *     auto particle = solid_particle_handler.begin(); 
 *     while (particle != solid_particle_handler.end()) 
 *       { 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 * @endcode
 * 
 * 我们从粒子本身得到一个通往粒子所在单元的迭代器。然后，我们就可以像通常那样在系统矩阵和右手边组装附加项了。
 * 

 * 
 * 
 * @code
 *         const auto &cell = particle->get_surrounding_cell(fluid_tria); 
 *         const auto &dh_cell = 
 *           typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh); 
 *         dh_cell->get_dof_indices(fluid_dof_indices); 
 * @endcode
 * 
 * 所以
 * 然后让我们得到位于这个单元格上的单元格集合，并对它们进行迭代。从每个粒子中，我们收集该粒子的位置和参考位置，以及附加在该粒子上的额外信息。在本例中，这些信息是用于生成粒子的正交点的 "JxW"。
 * 

 * 
 * 利用这些信息，我们可以将正交点的贡献加入到local_matrix和local_rhs中。我们可以利用每个粒子的参考位置，轻松地评估其位置上的形状函数值。
 * 

 * 
 * 
 * @code
 *         const auto pic = solid_particle_handler.particles_in_cell(cell); 
 *         Assert(pic.begin() == particle, ExcInternalError()); 
 *         for (const auto &p : pic) 
 *           { 
 *             const auto &ref_q  = p.get_reference_location(); 
 *             const auto &real_q = p.get_location(); 
 *             const auto &JxW    = p.get_properties()[0]; 
 * 
 *             for (unsigned int i = 0; i < fluid_fe->n_dofs_per_cell(); ++i) 
 *               { 
 *                 const auto comp_i = 
 *                   fluid_fe->system_to_component_index(i).first; 
 *                 if (comp_i < spacedim) 
 *                   { 
 *                     for (unsigned int j = 0; j < fluid_fe->n_dofs_per_cell(); 
 *                          ++j) 
 *                       { 
 *                         const auto comp_j = 
 *                           fluid_fe->system_to_component_index(j).first; 
 *                         if (comp_i == comp_j) 
 *                           local_matrix(i, j) += 
 *                             penalty_parameter * par.penalty_term * 
 *                             fluid_fe->shape_value(i, ref_q) * 
 *                             fluid_fe->shape_value(j, ref_q) * JxW; 
 *                       } 
 *                     local_rhs(i) += penalty_parameter * par.penalty_term * 
 *                                     solid_velocity.value(real_q, comp_i) * 
 *                                     fluid_fe->shape_value(i, ref_q) * JxW; 
 *                   } 
 *               } 
 *           } 
 * 
 *         constraints.distribute_local_to_global(local_matrix, 
 *                                                local_rhs, 
 *                                                fluid_dof_indices, 
 *                                                system_matrix, 
 *                                                system_rhs); 
 *         particle = pic.end(); 
 *       } 
 * 
 *     system_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solvingthelinearsystem"></a> 
 * <h4>Solving the linear system</h4>
 * 

 * 
 * 这个函数用FGMRES求解线性系统，有一个对角线块的预处理和一个对角线块的代数多重网格（AMG）方法。该预处理程序对 $(0,0)$ （即速度-速度）块应用V循环，对 $(1,1)$ 块应用质量矩阵的CG（这是我们对舒尔补码的近似值：上面组装的压力质量矩阵）。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::solve() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Solve"); 
 * 
 *     LA::MPI::PreconditionAMG prec_A; 
 *     { 
 *       LA::MPI::PreconditionAMG::AdditionalData data; 
 * 
 * #ifdef USE_PETSC_LA 
 *       data.symmetric_operator = true; 
 * #endif 
 *       prec_A.initialize(system_matrix.block(0, 0), data); 
 *     } 
 * 
 *     LA::MPI::PreconditionAMG prec_S; 
 *     { 
 *       LA::MPI::PreconditionAMG::AdditionalData data; 
 * 
 * #ifdef USE_PETSC_LA 
 *       data.symmetric_operator = true; 
 * #endif 
 *       prec_S.initialize(preconditioner_matrix.block(1, 1), data); 
 *     } 
 * 
 *     const auto A = linear_operator<LA::MPI::Vector>(system_matrix.block(0, 0)); 
 *     const auto amgA = linear_operator(A, prec_A); 
 * 
 *     const auto S = 
 *       linear_operator<LA::MPI::Vector>(preconditioner_matrix.block(1, 1)); 
 *     const auto amgS = linear_operator(S, prec_S); 
 * 
 *     ReductionControl          inner_solver_control(100, 
 *                                           1e-8 * system_rhs.l2_norm(), 
 *                                           1.e-2); 
 *     SolverCG<LA::MPI::Vector> cg(inner_solver_control); 
 * 
 *     const auto invS = inverse_operator(S, cg, amgS); 
 * 
 *     const auto P = block_diagonal_operator<2, LA::MPI::BlockVector>( 
 *       std::array< 
 *         dealii::LinearOperator<typename LA::MPI::BlockVector::BlockType>, 
 *         2>{{amgA, amgS}}); 
 * 
 *     SolverControl solver_control(system_matrix.m(), 
 *                                  1e-10 * system_rhs.l2_norm()); 
 * 
 *     SolverFGMRES<LA::MPI::BlockVector> solver(solver_control); 
 * 
 *     constraints.set_zero(solution); 
 * 
 *  
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations." 
 *           << std::endl; 
 * 
 *     constraints.distribute(solution); 
 * 
 *     locally_relevant_solution = solution; 
 *     const double mean_pressure = 
 *       VectorTools::compute_mean_value(fluid_dh, 
 *                                       QGauss<spacedim>(par.velocity_degree + 2), 
 *                                       locally_relevant_solution, 
 *                                       spacedim); 
 *     solution.block(1).add(-mean_pressure); 
 *     locally_relevant_solution.block(1) = solution.block(1); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Meshrefinement"></a> 
 * <h4>Mesh refinement</h4>
 * 

 * 
 * 我们以一种完全标准的方式处理网格细化问题。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::refine_and_transfer() 
 *   { 
 *     TimerOutput::Scope               t(computing_timer, "Refine"); 
 *     const FEValuesExtractors::Vector velocity(0); 
 * 
 *     Vector<float> error_per_cell(fluid_tria.n_active_cells()); 
 *     KellyErrorEstimator<spacedim>::estimate(fluid_dh, 
 *                                             QGauss<spacedim - 1>( 
 *                                               par.velocity_degree + 1), 
 *                                             {}, 
 *                                             locally_relevant_solution, 
 *                                             error_per_cell, 
 *                                             fluid_fe->component_mask(velocity)); 
 * 
 *     if (par.refinement_strategy == "fixed_fraction") 
 *       { 
 *         parallel::distributed::GridRefinement:: 
 *           refine_and_coarsen_fixed_fraction(fluid_tria, 
 *                                             error_per_cell, 
 *                                             par.refinement_fraction, 
 *                                             par.coarsening_fraction); 
 *       } 
 *     else if (par.refinement_strategy == "fixed_number") 
 *       { 
 *         parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
 *           fluid_tria, 
 *           error_per_cell, 
 *           par.refinement_fraction, 
 *           par.coarsening_fraction, 
 *           par.max_cells); 
 *       } 
 * 
 *     for (const auto &cell : fluid_tria.active_cell_iterators()) 
 *       { 
 *         if (cell->refine_flag_set() && 
 *             cell->level() == par.max_level_refinement) 
 *           cell->clear_refine_flag(); 
 *         if (cell->coarsen_flag_set() && 
 *             cell->level() == par.min_level_refinement) 
 *           cell->clear_coarsen_flag(); 
 *       } 
 * 
 *     parallel::distributed::SolutionTransfer<spacedim, LA::MPI::BlockVector> 
 *       transfer(fluid_dh); 
 *     fluid_tria.prepare_coarsening_and_refinement(); 
 *     transfer.prepare_for_coarsening_and_refinement(locally_relevant_solution); 
 *     fluid_tria.execute_coarsening_and_refinement(); 
 * 
 *     setup_dofs(); 
 * 
 *     transfer.interpolate(solution); 
 *     constraints.distribute(solution); 
 *     locally_relevant_solution = solution; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Creatingoutputforvisualization"></a> 
 * <h4>Creating output for visualization</h4>
 * 

 * 
 * 我们使用deal.II的标准并行功能在流体域上输出结果（速度和压力）。编写一个压缩的vtu文件，将所有处理器的信息聚集在一起。另外写一个`.pvd`记录，将物理时间与vtu文件联系起来。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void 
 *   StokesImmersedProblem<dim, spacedim>::output_results(const unsigned int cycle, 
 *                                                        double time) const 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Output fluid"); 
 * 
 *     std::vector<std::string> solution_names(spacedim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         spacedim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<spacedim> data_out; 
 *     data_out.attach_dof_handler(fluid_dh); 
 *     data_out.add_data_vector(locally_relevant_solution, 
 *                              solution_names, 
 *                              DataOut<spacedim>::type_dof_data, 
 *                              data_component_interpretation); 
 * 
 *     Vector<float> subdomain(fluid_tria.n_active_cells()); 
 *     for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *       subdomain(i) = fluid_tria.locally_owned_subdomain(); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(cycle) + ".vtu"; 
 *     data_out.write_vtu_in_parallel(par.output_directory + "/" + filename, 
 *                                    mpi_communicator); 
 * 
 *     static std::vector<std::pair<double, std::string>> times_and_names; 
 *     times_and_names.push_back(std::make_pair(time, filename)); 
 *     std::ofstream ofile(par.output_directory + "/" + "solution.pvd"); 
 *     DataOutBase::write_pvd_record(ofile, times_and_names); 
 *   } 
 * 
 * @endcode
 * 
 * 同样地，我们通过 Particles::DataOut 对象将粒子（无论是来自实体还是追踪器）写成一个单一的压缩vtu文件。这个简单的对象并不写作为 "属性 "附加到粒子上的额外信息，而只写它们的id--但是，无论如何，我们并不关心这些粒子位置的 "JxW "值，所以我们可能想要可视化的信息并没有丢失。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::output_particles( 
 *     const Particles::ParticleHandler<spacedim> &particles, 
 *     std::string                                 fprefix, 
 *     const unsigned int                          iter, 
 *     const double                                time) const 
 *   { 
 *     Particles::DataOut<spacedim> particles_out; 
 *     particles_out.build_patches(particles); 
 *     const std::string filename = 
 *       (fprefix + "-" + Utilities::int_to_string(iter) + ".vtu"); 
 *     particles_out.write_vtu_in_parallel(par.output_directory + "/" + filename, 
 *                                         mpi_communicator); 
 * 
 *     static std::map<std::string, std::vector<std::pair<double, std::string>>> 
 *       times_and_names; 
 *     if (times_and_names.find(fprefix) != times_and_names.end()) 
 *       times_and_names[fprefix].push_back(std::make_pair(time, filename)); 
 *     else 
 *       times_and_names[fprefix] = {std::make_pair(time, filename)}; 
 *     std::ofstream ofile(par.output_directory + "/" + fprefix + ".pvd"); 
 *     DataOutBase::write_pvd_record(ofile, times_and_names[fprefix]); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Therunfunction"></a> 
 * <h4>The "run" function</h4>
 * 

 * 
 * 这个函数现在负责协调整个模拟过程。它与其他时间相关的教程程序非常相似--以 step-21 或 step-26 为例。在开始的时候，我们会输出一些状态信息，同时将所有的当前参数保存到输出目录下的文件中，以利于重现。
 * 

 * 
 * 
 * @code
 *   template <int dim, int spacedim> 
 *   void StokesImmersedProblem<dim, spacedim>::run() 
 *   { 
 * #ifdef USE_PETSC_LA 
 *     pcout << "Running StokesImmersedProblem<" 
 *           << Utilities::dim_string(dim, spacedim) << "> using PETSc." 
 *           << std::endl; 
 * #else 
 *     pcout << "Running StokesImmersedProblem<" 
 *           << Utilities::dim_string(dim, spacedim) << "> using Trilinos." 
 *           << std::endl; 
 * #endif 
 *     par.prm.print_parameters(par.output_directory + "/" + "used_parameters_" + 
 *                                std::to_string(dim) + std::to_string(spacedim) + 
 *                                ".prm", 
 *                              ParameterHandler::Short); 
 * 
 * @endcode
 * 
 * 然后我们开始时间循环。我们在第一个循环中初始化模拟的所有元素
 * 

 * 
 * 
 * @code
 *     const double time_step    = par.final_time / (par.number_of_time_steps - 1); 
 *     double       time         = 0; 
 *     unsigned int output_cycle = 0; 
 * 
 *     for (unsigned int cycle = 0; cycle < par.number_of_time_steps; 
 *          ++cycle, time += time_step) 
 *       { 
 *         par.set_time(time); 
 *         pcout << "Cycle " << cycle << ':' << std::endl 
 *               << "Time : " << time << ", time step: " << time_step << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             make_grid(); 
 *             initial_setup(); 
 *             setup_dofs(); 
 *             setup_tracer_particles(); 
 *             setup_solid_particles(); 
 *             tracer_particle_velocities.reinit( 
 *               locally_owned_tracer_particle_coordinates, mpi_communicator); 
 *             output_results(output_cycle, time); 
 *             { 
 *               TimerOutput::Scope t(computing_timer, "Output tracer particles"); 
 *               output_particles(tracer_particle_handler, 
 *                                "tracer", 
 *                                output_cycle, 
 *                                time); 
 *             } 
 *             { 
 *               TimerOutput::Scope t(computing_timer, "Output solid particles"); 
 *               output_particles(solid_particle_handler, 
 *                                "solid", 
 *                                output_cycle, 
 *                                time); 
 *             } 
 *           } 
 * 
 * @endcode
 * 
 * 在第一个时间步长之后，我们在每个时间步长的开始时对实体进行位移，以考虑到它已经移动的事实。
 * 

 * 
 * 
 * @code
 *         else 
 *           { 
 *             TimerOutput::Scope t(computing_timer, 
 *                                  "Set solid particle position"); 
 * 
 *             SolidPosition<spacedim> solid_position(par.angular_velocity, 
 *                                                    time_step); 
 *             solid_particle_handler.set_particle_positions(solid_position, 
 *                                                           false); 
 *           } 
 * 
 * @endcode
 * 
 * 为了更新系统的状态，我们首先对示踪粒子位置的流体速度进行插值，并采用天真的显式欧拉方案对无质量示踪粒子进行漂移。
 * 

 * 
 * 
 * @code
 *         { 
 *           TimerOutput::Scope t(computing_timer, "Set tracer particle motion"); 
 *           Particles::Utilities::interpolate_field_on_particles( 
 *             fluid_dh, 
 *             tracer_particle_handler, 
 *             locally_relevant_solution, 
 *             tracer_particle_velocities, 
 *             fluid_fe->component_mask(FEValuesExtractors::Vector(0))); 
 * 
 *           tracer_particle_velocities *= time_step; 
 * 
 *           locally_relevant_tracer_particle_coordinates = 
 *             tracer_particle_handler.locally_owned_particle_ids().tensor_product( 
 *               complete_index_set(spacedim)); 
 * 
 *           relevant_tracer_particle_displacements.reinit( 
 *             locally_owned_tracer_particle_coordinates, 
 *             locally_relevant_tracer_particle_coordinates, 
 *             mpi_communicator); 
 * 
 *           relevant_tracer_particle_displacements = tracer_particle_velocities; 
 * 
 *           tracer_particle_handler.set_particle_positions( 
 *             relevant_tracer_particle_displacements); 
 *         } 
 * 
 * @endcode
 * 
 * 利用这些新的位置，我们就可以组装斯托克斯系统，并解决它。
 * 

 * 
 * 
 * @code
 *         assemble_stokes_system(); 
 *         assemble_nitsche_restriction(); 
 *         solve(); 
 * 
 * @endcode
 * 
 * 在适当的频率下，我们再将固体粒子、示踪粒子和流体领域的信息写入文件，以便进行可视化，并通过适应网格来结束时间步骤。
 * 

 * 
 * 
 * @code
 *         if (cycle % par.output_frequency == 0) 
 *           { 
 *             output_results(output_cycle, time); 
 *             { 
 *               TimerOutput::Scope t(computing_timer, "Output tracer particles"); 
 *               output_particles(tracer_particle_handler, 
 *                                "tracer", 
 *                                output_cycle, 
 *                                time); 
 *             } 
 *             { 
 *               TimerOutput::Scope t(computing_timer, "Output solid particles"); 
 *               output_particles(solid_particle_handler, 
 *                                "solid", 
 *                                output_cycle, 
 *                                time); 
 *             } 
 *             ++output_cycle; 
 *           } 
 *         if (cycle % par.refinement_frequency == 0 && 
 *             cycle != par.number_of_time_steps - 1) 
 *           refine_and_transfer(); 
 *       } 
 *   } 
 * 
 * } // namespace Step70 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 代码的其余部分，即`main()`函数，是标准的，除了对输入参数文件的处理。我们允许用户指定一个可选的参数文件作为程序的参数。如果没有指定，我们就使用默认文件 "parameters.prm"，如果不存在，我们就创建这个文件。文件名首先被扫描为字符串 "23"，然后是 "3"。如果文件名包含字符串 "23"，问题类将分别以模板参数2和3进行实例化。如果只找到 "3 "这个字符串，那么两个模板参数都被设置为3，否则都被设置为2。
 * 

 * 
 * 如果程序被调用时没有任何命令行参数（即`argc==1`），那么我们就默认使用 "参数.prm"。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   using namespace Step70; 
 *   using namespace dealii; 
 *   deallog.depth_console(1); 
 *   try 
 *     { 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       std::string prm_file; 
 *       if (argc > 1) 
 *         prm_file = argv[1]; 
 *       else 
 *         prm_file = "parameters.prm"; 
 * 
 *       if (prm_file.find("23") != std::string::npos) 
 *         { 
 *           StokesImmersedProblemParameters<2, 3> par; 
 *           ParameterAcceptor::initialize(prm_file); 
 * 
 *           StokesImmersedProblem<2, 3> problem(par); 
 *           problem.run(); 
 *         } 
 *       else if (prm_file.find("3") != std::string::npos) 
 *         { 
 *           StokesImmersedProblemParameters<3> par; 
 *           ParameterAcceptor::initialize(prm_file); 
 * 
 *           StokesImmersedProblem<3> problem(par); 
 *           problem.run(); 
 *         } 
 *       else 
 *         { 
 *           StokesImmersedProblemParameters<2> par; 
 *           ParameterAcceptor::initialize(prm_file); 
 * 
 *           StokesImmersedProblem<2> problem(par); 
 *           problem.run(); 
 *         } 
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
examples/step-70/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序的目录中包含一些样本参数文件，你可以用它们来重现本节介绍的结果。如果你没有在命令行中指定参数文件作为参数，程序将默认尝试读取文件"`parameters.prm`"，并执行二维版本的代码。正如在源代码的讨论中所解释的那样，如果你的文件名包含字符串 "23"，那么程序将运行一个三维问题，即共维度为1的沉入式实体。如果文件名包含字符串 "3"，它将运行一个三维问题，同维度的沉浸实体为零，否则它将运行一个二维问题，同维度的沉浸实体为零。

无论具体的参数文件名是什么，如果指定的文件不存在，当你执行程序时，你会得到一个异常，即找不到这样的文件。

@code


----------------------------------------------------
Exception on processing:


--------------------------------------------------------
An error occurred in line <74> of file <../source/base/parameter_acceptor.cc> in function
    static void dealii::ParameterAcceptor::initialize(const std::string &, const std::string &, const ParameterHandler::OutputStyle, dealii::ParameterHandler &)
The violated condition was:
    false
Additional information:
    You specified <parameters.prm> as input parameter file, but it does not exist. We created it for you.


--------------------------------------------------------


Aborting!


----------------------------------------------------
@endcode



然而，正如错误信息已经指出的，触发异常的代码也将生成指定的文件（"`parameters.prm`"在这种情况下），该文件仅仅包含这个程序关心的所有参数的默认值（对于正确的尺寸和辅助尺寸，根据文件名中是否包含字符串 "23 "或 "3"）。通过检查默认参数文件，我们看到以下内容。

@code
# Listing of Parameters
# ---------------------
subsection Stokes Immersed Problem
  set Final time                            = 1
  # Extraction level of the rtree used to construct global bounding boxes
  set Fluid bounding boxes extraction level = 1


  # Boundary Ids over which homogeneous Dirichlet boundary conditions are
  # applied
  set Homogeneous Dirichlet boundary ids    = 0


  # Initial mesh refinement used for the fluid domain Omega
  set Initial fluid refinement              = 5


  # Initial mesh refinement used for the solid domain Gamma
  set Initial solid refinement              = 5
  set Nitsche penalty term                  = 100
  set Number of time steps                  = 501
  set Output directory                      = .
  set Output frequency                      = 1


  # Refinement of the volumetric mesh used to insert the particles
  set Particle insertion refinement         = 3
  set Velocity degree                       = 2
  set Viscosity                             = 1



  subsection Angular velocity
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  =


    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = t < .500001 ? 6.283185 : -6.283185 # default: 0


    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end


  subsection Grid generation
    set Fluid grid generator              = hyper_cube
    set Fluid grid generator arguments    = -1: 1: false
    set Particle grid generator           = hyper_ball
    set Particle grid generator arguments = 0.3, 0.3: 0.1: false
    set Solid grid generator              = hyper_rectangle
    set Solid grid generator arguments    = -.5, -.1: .5, .1: false
  end


  subsection Refinement and remeshing
    set Maximum number of cells        = 20000
    set Refinement coarsening fraction = 0.3
    set Refinement fraction            = 0.3
    set Refinement maximal level       = 8
    set Refinement minimal level       = 5
    set Refinement step frequency      = 5
    set Refinement strategy            = fixed_fraction
  end


  subsection Right hand side
    # Sometimes it is convenient to use symbolic constants in the expression
    # that describes the function, rather than having to use its numeric value
    # everywhere the constant appears. These values can be defined using this
    # parameter, in the form `var1=value1, var2=value2, ...'.
    #
    # A typical example would be to set this runtime parameter to
    # `pi=3.1415926536' and then use `pi' in the expression of the actual
    # formula. (That said, for convenience this class actually defines both
    # `pi' and `Pi' by default, but you get the idea.)
    set Function constants  =


    # The formula that denotes the function you want to evaluate for
    # particular values of the independent variables. This expression may
    # contain any of the usual operations such as addition or multiplication,
    # as well as all of the common functions such as `sin' or `cos'. In
    # addition, it may contain expressions like `if(x>0, 1, -1)' where the
    # expression evaluates to the second argument if the first argument is
    # true, and to the third argument otherwise. For a full overview of
    # possible expressions accepted see the documentation of the muparser
    # library at http://muparser.beltoforion.de/.
    #
    # If the function you are describing represents a vector-valued function
    # with multiple components, then separate the expressions for individual
    # components by a semicolon.
    set Function expression = 0; 0; 0


    # The names of the variables as they will be used in the function,
    # separated by commas. By default, the names of variables at which the
    # function will be evaluated are `x' (in 1d), `x,y' (in 2d) or `x,y,z' (in
    # 3d) for spatial coordinates and `t' for time. You can then use these
    # variable names in your function expression and they will be replaced by
    # the values of these variables at which the function is currently
    # evaluated. However, you can also choose a different set of names for the
    # independent variables at which to evaluate your function expression. For
    # example, if you work in spherical coordinates, you may wish to set this
    # input parameter to `r,phi,theta,t' and then use these variable names in
    # your function expression.
    set Variable names      = x,y,t
  end


end
@endcode



如果你现在运行该程序，你会在参数`Output directory`（默认为当前目录）指定的目录下得到一个名为`parameters_22.prm`的文件，其中包含上述参数的简短版本（没有注释和文档），记录了所有用于运行程序的参数。

@code
subsection Stokes Immersed Problem
  set Final time                            = 1
  set Fluid bounding boxes extraction level = 1
  set Homogeneous Dirichlet boundary ids    = 0
  set Initial fluid refinement              = 5
  set Initial solid refinement              = 5
  set Nitsche penalty term                  = 100
  set Number of time steps                  = 501
  set Output directory                      = .
  set Output frequency                      = 1
  set Particle insertion refinement         = 3
  set Velocity degree                       = 2
  set Viscosity                             = 1
  subsection Angular velocity
    set Function constants  =
    set Function expression = t < .500001 ? 6.283185 : -6.283185
    set Variable names      = x,y,t
  end
  subsection Grid generation
    set Fluid grid generator              = hyper_cube
    set Fluid grid generator arguments    = -1: 1: false
    set Particle grid generator           = hyper_ball
    set Particle grid generator arguments = 0.3, 0.3: 0.1: false
    set Solid grid generator              = hyper_rectangle
    set Solid grid generator arguments    = -.5, -.1: .5, .1: false
  end
  subsection Refinement and remeshing
    set Maximum number of cells        = 20000
    set Refinement coarsening fraction = 0.3
    set Refinement fraction            = 0.3
    set Refinement maximal level       = 8
    set Refinement minimal level       = 5
    set Refinement step frequency      = 5
    set Refinement strategy            = fixed_fraction
  end
  subsection Right hand side
    set Function constants  =
    set Function expression = 0; 0; 0
    set Variable names      = x,y,t
  end
end
@endcode



首先创建 "parameters.prm "文件（程序第一次运行时），然后创建 "output/parameters_22.prm "文件（每次使用现有的输入文件运行程序时），这是因为你可能想把大多数参数保留为默认值，只修改其中的一小部分，同时仍然能够重现结果，检查特定模拟使用了哪些参数。一般来说，将用于模拟的参数文件与模拟输出一起保存起来是很好的科学做法，这样你就可以在以后的时间里重复相同的运行。

另一个原因是输入文件可能只包含那些与默认值不同的参数。例如，你可以在本教程程序中使用以下（完全有效的）参数文件。

@code
subsection Stokes Immersed Problem
  set Final time                         = 1
  set Nitsche penalty term               = 10
  set Number of time steps               = 101
  set Velocity degree                    = 3
end
@endcode

你将使用Q3/Q2 Taylor-Hood有限元运行程序，进行101步，使用Nitsche惩罚为`10`，并将所有其他参数保持为默认值。输出目录不仅包含了这些参数的记录，而且包含了仿真中使用的所有参数。你可以在生成的文件`parameters_22.prm`中查看所有其他参数。




<a name="Twodimensionaltestcase"></a><h3> Two dimensional test case </h3>


默认问题产生了一个同维度的零叶轮，由一个旋转的矩形网格组成，在一个方向上旋转半个时间单位，在相反方向上旋转半个时间单位，恒定的角速度等于 $\approx 2\pi \frac{\text{rad}}{\text{time unit}}$  。因此，叶轮做了半个旋转，并返回到原来的位置。下面的动画显示了速度的大小，固体叶轮和示踪粒子的运动。


<p align="center"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-70.2d_tracing.gif" alt="" width="500"> </div>  </p> 

在一个核心上，程序的输出将看起来像下面这样。

@code
bash$ mpirun -np 1 ./step-70 test.prm
Running StokesImmersedProblem<2> using Trilinos.
Cycle 0:
Time : 0, time step: 0.002
   Number of degrees of freedom: 9539 (8450+1089 -- 0+0)
Tracer particles: 337
Solid particles: 9216
   Solved in 158 iterations.
   Number of degrees of freedom: 9845 (8722+1123 -- 9216+337)
Cycle 1:
Time : 0.002, time step: 0.002
   Solved in 142 iterations.
Cycle 2:
Time : 0.004, time step: 0.002
   Solved in 121 iterations.
Cycle 3:
Time : 0.006, time step: 0.002
   Solved in 121 iterations.


...


Cycle 499:
Time : 0.998, time step: 0.002
   Solved in 199 iterations.
Cycle 500:
Time : 1, time step: 0.002
   Solved in 196 iterations.


+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |       302s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble Nitsche terms          |       501 |      43.3s |        14% |
| Assemble Stokes terms           |       501 |      21.5s |       7.1% |
| Initial setup                   |         1 |  0.000792s |         0% |
| Output fluid                    |       502 |      31.8s |        11% |
| Output solid particles          |       502 |      32.2s |        11% |
| Output tracer particles         |       502 |      0.61s |       0.2% |
| Refine                          |       100 |      4.68s |       1.5% |
| Set solid particle position     |       500 |      3.34s |       1.1% |
| Set tracer particle motion      |       501 |     0.729s |      0.24% |
| Setup dofs                      |       101 |       2.2s |      0.73% |
| Solve                           |       501 |       164s |        54% |
+---------------------------------+-----------+------------+------------+
@endcode



你可能会注意到，组装耦合系统比组装斯托克斯部分更昂贵。这在很大程度上取决于用于应用Nitsche限制的高斯点（固体粒子）的数量。在目前的情况下，所使用的示踪粒子的数量相对较少。因此，跟踪它们的运动是相对便宜的。

下面的影片显示了解决方案随时间的演变。

@htmlonly
<p align="center">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/y4Gypj2jpXw"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly



这部电影显示了灰色的旋转障碍物（实际上是用足够大的点绘制的固体粒子的叠加，使它们重叠），浅色的<a
href="https://en.wikipedia.org/wiki/Streamlines,_streaklines,_and_pathlines">streamlines
of the fluid flow</a>（包括在模拟过程中特定时间形成的角顶点），以及蓝色色调的示踪粒子。

模拟结果显示，在结束的时候，示踪剂颗粒已经在一定程度上回到了原来的位置，尽管它们已经被流场扭曲了。下面的图片比较了粒子在一个时间单位的流动后的初始和最终位置。

<p align="center"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-70.tracer_comparison.png" alt="" width="500"> </div>  </p> 

在这种情况下，我们看到在叶轮扫过的体积之外的示踪剂颗粒已经非常接近它们的初始位置，而在扫过的体积内的示踪剂颗粒的变形略大。这种变形是非物理性的。它是由用于平移粒子的显式欧拉方案引起的数值误差、由虚构领域引起的精度损失以及最后由斯托克斯方程的离散化误差引起的。前两个错误是造成这种变形的主要原因，它们可以通过使用更细的网格和更小的时间步长来缓解。




<a name="Threedimensionaltestcase"></a><h3> Three dimensional test case </h3>


为了玩一玩，我们将虚构的领域复杂化（取自https://grabcad.com/library/lungstors-blower-1），并在三个空间维度上运行共维一模拟，使用以下"`参数_23.prm`"文件。

@code
subsection Stokes Immersed Problem
  set Final time                            = 1
  set Homogeneous Dirichlet boundary ids    = 0
  set Fluid bounding boxes extraction level = 1
  set Initial fluid refinement              = 3
  set Initial solid refinement              = 0
  set Nitsche penalty term                  = 10
  set Number of time steps                  = 101
  set Output frequency                      = 1
  set Particle insertion refinement         = 3
  set Velocity degree                       = 2
  set Viscosity                             = 1
  subsection Angular velocity
    set Function constants  =
    set Function expression = t < .500001 ? 5 : -5
    set Variable names      = x,y,z,t
  end
  subsection Grid generation
    set Fluid grid generator              = hyper_rectangle
    set Fluid grid generator arguments    = -50,-50, -10: 50, 50, 40: false
    set Solid grid generator              = impeller.vtk
    set Solid grid generator arguments    = 1:impeller.step
    set Particle grid generator           = hyper_ball
    set Particle grid generator arguments = 30, 30, 20: 10: false
  end
  subsection Refinement and remeshing
    set Maximum number of cells        = 100000
    set Refinement coarsening fraction = 0.3
    set Refinement fraction            = 0.3
    set Refinement maximal level       = 6
    set Refinement step frequency      = 5
    set Refinement strategy            = fixed_fraction
  end
  subsection Right hand side
    set Function constants  =
    set Function expression = 0; 0; 0; 0
    set Variable names      = x,y,z,t
  end
end
@endcode



在这种情况下，定时输出有点不同。

@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |  5.54e+03s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble Nitsche terms          |       101 |       111s |         2% |
| Assemble Stokes terms           |       101 |       208s |       3.8% |
| Initial setup                   |         1 |   0.00187s |         0% |
| Output fluid                    |       102 |      15.5s |      0.28% |
| Output solid particles          |       102 |      2.63s |         0% |
| Output tracer particles         |       102 |      2.49s |         0% |
| Refine                          |        20 |      18.4s |      0.33% |
| Set solid particle position     |       100 |       6.1s |      0.11% |
| Set tracer particle motion      |       101 |      10.8s |       0.2% |
| Setup dofs                      |        21 |      13.9s |      0.25% |
| Solve                           |       101 |  5.16e+03s |        93% |
+---------------------------------+-----------+------------+------------+
@endcode



现在，求解器在三维空间中占用了大部分的求解时间，就运行时间而言，粒子运动和Nitsche装配仍然相对不重要。




@htmlonly
<p align="center">
  <iframe width="560" height="315" src="https://www.youtube.com/embed/Srwq7zyR9mg"
   frameborder="0"
   allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture"
   allowfullscreen></iframe>
 </p>
@endhtmlonly




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


目前的教程程序显示了流体和固体之间的单向耦合，其中固体运动是强加的（而不是求解的），并通过利用固体正交点的位置和权重在固体域中读取。

代码的结构已经允许人们通过利用读取实体网格正交点上流体速度值的可能性来实现双向耦合。为了提高MPI通信模式的效率，我们应该将正交点的所有权保持在实体处理器上，该处理器拥有创建这些正交点的单元。在目前的代码中，通过使用实体分区而不是初始流体分区来定义用于交换正交点信息的向量索引集就足够了。

这使得本教程程序中使用的技术与教程步骤-60中提出的技术相结合，以解决带有分布式拉格朗日乘数的流体结构交互问题，在 parallel::distributed::Triangulation 对象上。

上面的时间显示，目前的预处理策略对Nitsche惩罚的效果并不好，如果我们想瞄准更大的问题，我们应该想出一个更好的预处理方法。此外，应该实施检查点重启策略，以允许较长的模拟被中断和恢复，例如在step-69教程中就是这样做的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-70.cc"
*/
