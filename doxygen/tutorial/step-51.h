/**
@page step_51 The step-51 tutorial program
This tutorial depends on step-7, step-9, step-61.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#HybridizablediscontinuousGalerkinmethods"> Hybridizable discontinuous Galerkin methods </a>
      <ul>
        <li><a href="#Reducingthesizeofthelinearsystem"> Reducing the size of the linear system </a>
        <li><a href="#RelationwithStaticCondensation"> Relation with Static Condensation </a>
        <li><a href="#Solutionqualityandratesofconvergence"> Solution quality and rates of convergence</a>
        <li><a href="#Alternativeapproaches"> Alternative approaches </a>
      </ul>
        <li><a href="#HDGappliedtotheconvectiondiffusionproblem"> HDG applied to the convection-diffusion problem </a>
      <ul>
        <li><a href="#Postprocessingandsuperconvergence"> Post-processing and super-convergence </a>
      </ul>
        <li><a href="#Problemspecificdata"> Problem specific data </a>
        <li><a href="#Implementation"> Implementation </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#TheHDGsolverclass">The HDG solver class</a>
        <li><a href="#TheHDGclassimplementation">The HDG class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#HDGsetup_system">HDG::setup_system</a>
        <li><a href="#HDGPerTaskDatapPerTaskDataScratchDataptrace_reconstructHDGDCA1">HDG::PerTaskData}  接下来是定义并行装配的本地数据结构。第一个结构 @p PerTaskData 包含了被写入全局矩阵的本地向量和矩阵，而ScratchData包含了我们在本地装配中需要的所有数据。这里有一个变量值得注意，即布尔变量 @p  trace_reconstruct。正如介绍中提到的，我们分两步解决HDG系统。首先，我们为骨架系统创建一个线性系统，通过舒尔补码 $D-CA^{-1</a>
        <li><a href="#HDGScratchData">HDG::ScratchData</a>
        <li><a href="#HDGPostProcessScratchData">HDG::PostProcessScratchData</a>
        <li><a href="#HDGassemble_system">HDG::assemble_system</a>
        <li><a href="#HDGassemble_system_one_cell">HDG::assemble_system_one_cell</a>
        <li><a href="#HDGcopy_local_to_global">HDG::copy_local_to_global</a>
        <li><a href="#HDGsolve">HDG::solve</a>
        <li><a href="#HDGpostprocess">HDG::postprocess</a>
        <li><a href="#HDGpostprocess_one_cell">HDG::postprocess_one_cell</a>
        <li><a href="#HDGoutput_results">HDG::output_results</a>
        <li><a href="#HDGrefine_grid">HDG::refine_grid</a>
        <li><a href="#HDGrun">HDG::run</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Programoutput">Program output</a>
      <ul>
        <li><a href="#Convergencetables">Convergence tables</a>
      </ul>
        <li><a href="#Comparisonwithcontinuousfiniteelements">Comparison with continuous finite elements</a>
      <ul>
        <li><a href="#Resultsfor2D">Results for 2D</a>
        <li><a href="#Resultsfor3D">Results for 3D</a>
      </ul>
        <li><a href="#Possibilitiesforimprovements">Possibilities for improvements</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-51/doc/intro.dox

 <br> 

<i>
This program was contributed by Martin Kronbichler and Scott Miller.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序介绍了对流-扩散方程的可混合的不连续Galkerin方法的实现。

<a name="HybridizablediscontinuousGalerkinmethods"></a><h3> Hybridizable discontinuous Galerkin methods </h3>


反对使用非连续Galerkin元素的一个常见论点是，在隐式系统中必须解决大量的全局耦合自由度。  这是因为，与连续有限元不同，在典型的非连续元中，每个顶点有一个自由度<i>for each of the adjacent elements</i>，而不是只有一个，对边和面也是如此。  作为未知数增长速度的例子，请考虑FE_DGPMonomial基础：每个标量解分量都由度数为 $p$ 的多项式表示，每个元素有 $(1/\text{dim}!) \prod_{i=1}^{\text{dim}}(p+i)$ 个自由度。通常，一个元素的所有自由度都与相邻元素的所有自由度相耦合。  由此产生的离散方程会很快产生非常大的线性系统，特别是对于2或3维的方程系统。

<a name="Reducingthesizeofthelinearsystem"></a><h4> Reducing the size of the linear system </h4>为了减轻解决这种大型线性系统的计算成本，Cockburn和同事们引入了可混合的非连续Galerkin（HDG）方法（见Nguyen和Peraire最近发表的HDG概述文章中的参考资料 @cite Ngu2012 ）。


HDG方法通过使用Dirichlet-to-Neumann映射对数学问题进行表述来实现这一目标。  偏微分方程首先被写成一阶系统，然后每个场通过DG方法进行离散。  在这一点上，网格骨架上的单值 "跟踪 "值，即元素面，被视为独立的未知量。这就产生了离散公式中的未知数，这些未知数分为两类。

- 面的未知数，只与面的两边的单元格未知数耦合。

- 单元未知数只与同一单元内定义的单元和面未知数相耦合。最重要的是，一个单元格的内部自由度不会与另一个单元格的任何内部自由度相耦合。

然后，Dirichlet-to-Neumann地图的概念允许以下解决程序。<ol>  <li>  使用局部元素内部数据来强制执行三角形骨架上的Neumann条件。  然后，全局问题是求解轨迹值，这是唯一全局耦合的未知数。     <li>  使用已知的骨架值作为Dirichlet数据来求解局部元素级的解决方案。  这被称为 "局部求解器"，是一个<i>embarrassingly parallel</i>逐个元素的求解过程。   </ol> 

<a name="RelationwithStaticCondensation"></a><h4> Relation with Static Condensation </h4>上述程序也有线性代数的解释--被称为<i>static condensation</i>--被Guyan在连续有限元的背景下 @cite G65 ，以及被Fraeijs de Veubeke用于混合方法 @cite F65 ，用来减少全局线性系统的大小。在后一种情况下（混合公式），系统的减少是通过使用不连续的通量，结合引入一个额外的辅助变量<i>hybrid</i>来实现的，该变量在每个元素的边界上近似于未知数的轨迹。这个过程被称为混合化，并且通过类比，这也是为什么Cockburn、Gopalakrishnan和Lazarov在2009年引入的局部非连续Galerkin方法 @cite CGL2009 ，以及随后由他们的合作者开发，最终被称为<i>hybridizable discontinuous Galerkin</i>（HDG）方法的原因。


让我们把与HDG问题相关的完整线性系统写成一个块状系统，离散DG（单元内部）变量 $U$ 为第一块，骨架（面）变量 $\Lambda$ 为第二块。

@f{eqnarray*}
\begin{pmatrix} A & B \\ C & D \end{pmatrix}
\begin{pmatrix} U \\ \Lambda \end{pmatrix}
=
\begin{pmatrix} F \\ G \end{pmatrix}.


@f}

我们现在的目的是用类似于步骤20的舒尔补码方法来消除 $U$ 块，这导致了以下两个步骤。

@f{eqnarray*}
(D - C A^{-1} B) \Lambda &=& G - C A^{-1} F, \\
A U &=& F - B \Lambda.


@f}

关键是 $A^{-1}$ 的存在不是问题，因为 $A$ 是一个块对角线矩阵，每个块对应一个单元，因此足够容易反转。与其他单元的耦合是由骨架变量上的矩阵 $B$ 和 $C$ 引入的。 $A$ 的块对角性以及 $B$ 和 $C$ 的结构使我们能够逐元反转矩阵 $A$ （迪里希特问题的局部解），并从 $D$ 中减去 $CA^{-1}B$ 。因此，迪里切特到诺曼映射概念的步骤对应于<ol>  <li>  构建舒尔补码矩阵 $D-C A^{-1} B$ 和右手边 $G - C A^{-1} F$  <i>locally on each cell</i>并以通常的方式将贡献插入全局跟踪矩阵， <li> 求解舒尔互补系统 $\Lambda$ ， <li> 使用第二个方程求解 $U$ ，给出 $\Lambda$  。   </ol> 




<a name="Solutionqualityandratesofconvergence"></a><h4> Solution quality and rates of convergence</h4> 对传统DG方法的另一个批评是，近似通量的收敛是次优的。  局部HDG解可以被证明是收敛的 $\mathcal{O}(h^{p+1})$  ，即以最优顺序收敛。  此外，还可以利用超级收敛特性对新的近似解进行后处理，使其以  $\mathcal{O}(h^{p+2})$  的速率收敛。





<a name="Alternativeapproaches"></a><h4> Alternative approaches </h4>


可混合的非连续Galerkin方法只是解决非连续Galerkin方法的问题的一种方法。另一个想法是所谓的 "弱Galerkin "方法。它在步骤61中进行了探讨。




<a name="HDGappliedtotheconvectiondiffusionproblem"></a><h3> HDG applied to the convection-diffusion problem </h3>


本例中使用的HDG公式取自 <br> <b>
  N.C. Nguyen, J. Peraire, B. Cockburn:
  <i>An implicit high-order hybridizable discontinuous Galerkin method
  for linear convection–diffusion equations</i><i>An implicit high-order hybridizable discontinuous Galerkin method
  for linear convection–diffusion equations</i>,
  Journal of Computational Physics, 2009, 228:9, 3232-3254.
  <a href="http://dx.doi.org/10.1016/j.jcp.2009.01.030">[DOI]</a><a href="http://dx.doi.org/10.1016/j.jcp.2009.01.030">[DOI]</a>
</b>。

我们考虑域 $\Omega$ 上的对流-扩散方程，该方程具有迪里切特边界 $\partial \Omega_D$ 和诺伊曼边界 $\partial \Omega_N$  。

@f{eqnarray*}
	\nabla \cdot (\mathbf{c} u) - \nabla \cdot (\kappa \nabla u) &=& f,
	\quad \text{ in } \Omega, \\
	u &=& g_D, \quad \text{ on } \partial \Omega_D, \\
	(\mathbf{c} u - \kappa \nabla u)\cdot \mathbf{n} &=& g_N,
	\quad \text{ on }  \partial \Omega_N.


@f}



引入辅助变量 $\mathbf{q}=-\kappa \nabla u$ ，将上述方程改写为一阶系统。

@f{eqnarray*}
  \mathbf{q} + \kappa \nabla u &=& 0, \quad \text{ in } \Omega, \\
  \nabla \cdot (\mathbf{c} u + \mathbf{q}) &=& f, \quad \text{ in } \Omega, \\
  u &=& g_D, \quad \text{ on } \partial \Omega_D, \\
  (\mathbf{q} + \mathbf{c}u)\cdot\mathbf{n}  &=& g_N,
	\quad \text{ on }  \partial \Omega_N.


@f}



我们将这些方程乘以权重函数 $\mathbf{v}, w$ ，并对每个元素 $K$ 进行分项积分，得到。

@f{eqnarray*}
  (\mathbf{v}, \kappa^{-1} \mathbf{q})_K - (\nabla\cdot\mathbf{v}, u)_K
    + \left<\mathbf{v}\cdot\mathbf{n}, {\hat{u}}\right>_{\partial K} &=& 0, \\


  - (\nabla w, \mathbf{c} u + \mathbf{q})_K
    + \left<w, (\widehat{\mathbf{c} u}+{\hat{\mathbf{q}}})\cdot\mathbf{n}\right>_{\partial K}
    &=& (w,f)_K.


@f}



带帽子的术语表示数值轨迹（通常也被称为数值通量）。  它们是对元素边界上的内部值的近似。  为了确保守恒，这些项在任何给定的元素边上都必须是单值的 $\partial K$ ，尽管对于不连续的形状函数，当然可能有多个值来自界面附近的单元。我们通过使用以下形式的跟踪来消除数字跟踪 $\hat{\mathbf{q}}$ 。

@f{eqnarray*}
  \widehat{\mathbf{c} u}+\hat{\mathbf{q}} = \mathbf{c}\hat{u} + \mathbf{q}
  + \tau(u - \hat{u})\mathbf{n} \quad \text{ on } \partial K.


@f}



变量 $\hat {u}$ 作为一个额外的自变量被引入，是我们最终建立一个全局耦合线性系统的变量。如上所述，它被定义在元素面上，并且在面与面的交汇处（2D中的顶点，3D中的边缘和顶点）不连续。数值跟踪函数中出现的 $u$ 和 $\mathbf{q}$ 的值被认为是限制在边界 $\partial K$ 的单元内部解。

局部稳定参数 $\tau$ 对HDG解决方案的稳定性和准确性有影响；进一步的讨论见文献。据报道，稳定参数为1是给出最佳结果的选择。趋向于无穷大的稳定参数 $\tau$ 禁止解在元素边界上的跳跃，使HDG解接近连续有限元素的近似值。在下面的程序中，我们选择稳定参数为

@f{eqnarray*}
  \tau = \frac{\kappa}{\ell} + |\mathbf{c} \cdot \mathbf{n}|


@f}

其中我们设定扩散 $\kappa=1$ 和扩散长度尺度为 $\ell = \frac{1}{5}$  。

HDG方法中的轨迹/骨架变量在元素面上是单值的。  因此，它们必须强烈地代表 $\partial\Omega_D$ 上的迪里希特数据。  这意味着

@f{equation*}
  \hat{u}|_{\partial \Omega_D} = g_D,


@f}

其中等号实际上是指边界函数 $L_2$ 对脸部变量空间的 $g$ 投射（例如脸部的线性函数）。然后，这个约束被应用于骨架变量 $\hat{u}$ ，使用非均质约束的方法 VectorTools::project_boundary_values. 。

将三角结构中所有元素的贡献相加，强制执行数值通量的法线分量，并对由 $w$ 加权的方程进行分项积分，我们得出问题的最终形式：找到 $(\mathbf{q}_h, u_h, \hat{u}_h) \in
\mathcal{V}_h^p \times \mathcal{W}_h^p \times \mathcal{M}_h^p$ ，以便

@f{align*}
  (\mathbf{v}, \kappa^{-1} \mathbf{q}_h)_{\mathcal{T}}


    - ( \nabla\cdot\mathbf{v}, u_h)_{\mathcal{T}}
    + \left<\mathbf{v}\cdot\mathbf{n}, \hat{u}_h\right>_{\partial\mathcal{T}}
    &= 0,
    \quad &&\forall \mathbf{v} \in \mathcal{V}_h^p,
\\


   - (\nabla w, \mathbf{c} u_h)_{\mathcal{T}}
   + (w, \nabla \cdot \mathbf{q}_h)_{\mathcal{T}}
   + (w, (\mathbf{c}\cdot\mathbf{n}) \hat{u}_h)_{\partial \mathcal{T}}
    + \left<w, \tau (u_h - \hat{u}_h)\right>_{\partial \mathcal{T}}
    &=
    (w, f)_{\mathcal{T}},
    \quad &&\forall w \in \mathcal{W}_h^p,
\\
  \left< \mu, \hat{u}_h\mathbf{c} \cdot \mathbf{n}
  		+ \mathbf{q}_h\cdot \mathbf{n}
  	    + \tau (u_h - \hat{u}_h)\right>_{\partial \mathcal{T}}
    &=
    \left<\mu, g_N\right>_{\partial\Omega_N},
    \quad &&\forall \mu \in \mathcal{M}_h^p.


@f}



未知数 $(\mathbf{q}_h, u_h)$ 被称为局部变量；它们被表示为标准的DG变量。  未知数 $\hat{u}_h$ 是骨架变量，在网格的一维表面（面）上有支持。

我们用符号 $(\cdot, \cdot)_{\mathcal{T}} = \sum_K (\cdot, \cdot)_K$ 表示所有单元的积分之和， $\left<\cdot,
\cdot\right>_{\partial \mathcal{T}} = \sum_K \left<\cdot,
\cdot\right>_{\partial K}$ 表示所有单元的所有面的积分，也就是说，内部面被访问两次，一次来自每侧，并有相应的法向量。当结合共享一个面的两个元素的贡献时，上述方程产生了DG方法中熟悉的条款，解在单元格边界上有跳跃性。

在上述方程中，标量变量  $\mathcal {W}_h^{p}$  的空间  $u_h$  被定义为在每个单元上为张量积多项式  $p$  且在元素边界上不连续的函数空间  $\mathcal Q_{-p}$  ，即由  <code>FE_DGQ<dim>(p)</code>  描述的空间。梯度或通量变量的空间  $\mathbf{q}_i$  是一个矢量元素空间，其中每个分量是局部多项式且不连续  $\mathcal Q_{-p}$  。在下面的代码中，我们将这两个局部部分收集在一个FESystem中，其中第一个 @p dim组件表示梯度部分，最后一个标量组件对应标量变量。对于骨架部分 $\hat{u}_h$ ，我们定义了一个由住在元素面上的不连续张量乘积多项式组成的空间，在deal.II中由FE_FaceQ类实现。这个空间在其他方面与FE_DGQ相似，即解函数在两个相邻的面之间不连续，也可参见下面的结果部分进行说明。

在上面给出的弱形式中，我们可以注意到以下的耦合模式。<ol>  <li>  矩阵 $A$  由局部-局部耦合项组成。  当局部加权函数 $(\mathbf{v}, w)$ 与局部求解项 $(\mathbf{q}_h, u_h)$ 相乘时就会产生这些耦合项。因为这些元素是不连续的， $A$ 是块对角线。     <li>  矩阵 $B$ 代表局部面的耦合。  这些是具有加权函数 $(\mathbf{v}, w)$ 的条款，乘以骨架变量 $\hat{u}_h$  。     <li>  矩阵 $C$ 代表面-本地耦合，它涉及加权函数 $\mu$ 乘以本地解 $(\mathbf{q}_h, u_h)$  。     <li>  矩阵  $D$  是面-面耦合；条款涉及  $\mu$  和  $\hat{u}_h$  。   </ol> 

<a name="Postprocessingandsuperconvergence"></a><h4> Post-processing and super-convergence </h4>


HDG方法的一个特点是，它们通常允许构建一个丰富的解决方案，以提高精度。这种后处理方法以逐个元素的方式获取HDG解决方案，并将其结合起来，从而在使用度数为 $p$ 的多项式时可以获得 $\mathcal O(h^{p+2})$ 的精度。要做到这一点，有两个必要的成分。<ol>  <li>  计算的解梯度  $\mathbf{q}_h$  以最佳速度收敛，即  $\mathcal{O}(h^{p+1})$  。     <li>  解的标量部分的单元平均数， $\frac{(1,u_h)_K}{\text{vol}(K)}$  ，以 $\mathcal{O}(h^{p+2})$  的速度超级收敛。   </ol> 

我们现在引入一个新的变量 $u_h^* \in \mathcal{V}_h^{p+1}$ ，我们通过在约束条件 $\left(1, u_h^*\right)_K = \left(1,
u_h\right)_K$ 下对单元格 $K$ 的表达式进行最小化来找到它。这个约束是必要的，因为最小化函数并不能确定 $u_h^*$ 的常数部分。这就转化为以下方程组。

@f{eqnarray*}
\left(1, u_h^*\right)_K &=& \left(1, u_h\right)_K\\
\left(\nabla w_h^*, \kappa \nabla u_h^*\right)_K &=&


-\left(\nabla w_h^*, \mathbf{q}_h\right)_K
\quad \text{for all } w_h^* \in \mathcal Q^{p+1}.


@f}



由于我们在第二组方程中用度数为 $p+1$ 的张量积多项式空间中的整组基函数进行测试，这是一个过度确定的系统，方程比未知数多一个。我们在下面的代码中通过省略其中一个方程来解决这个问题（因为拉普拉斯的行在代表一个常数函数时是线性依赖的）。正如我们将在下面看到的，这种形式的后处理给出了所需的超级收敛结果，速率为 $\mathcal {O}(h^{p+2})$  。  应该指出的是，在构建 $u_h^*$ 时有一定的自由度，这种从梯度中提取信息的最小化方法不是唯一的方法。特别是，这里定义的后处理方案在任何意义上都不满足对流-扩散方程。作为替代方案，上面引用的Nguyen、Peraire和Cockburn的论文提出了另一个有点复杂的对流-扩散公式，该公式也可以将通量变量后处理为 $H(\Omega,\mathrm{div})$ -符合的变体，并且在扩散较小时更好地表示局部对流-扩散算子。我们把更复杂的后处理的实现作为一个可能的扩展留给感兴趣的读者。

请注意，对于矢量值的问题，后处理的工作原理是类似的。我们只需为每个向量分量的平均值分别设置约束，并将梯度作为主要信息来源。

<a name="Problemspecificdata"></a><h3> Problem specific data </h3>


在这个教程程序中，我们考虑的测试案例与步骤7中的几乎相同。计算域是 $\Omega \dealcoloneq [-1,1]^d$ ，精确的解决方案与步骤7中的解决方案相对应，除了一个缩放比例。我们使用以下源中心 $x_i$ 作为指数 <ul>   <li>  1D： $\{x_i\}^1 = \{ -\frac{1}{3}, 0, \frac{1}{3} \}$  ,  <li>  2D： $\{\mathbf{x}_i\}^2 = \{ (-\frac{1}{2},\frac{1}{2}),
                        		 (-\frac{1}{2},-\frac{1}{2}),
  					 (\frac{1}{2},-\frac{1}{2})
  				   \}$  ,  <li>  3D： $\{\mathbf{x}_i\}^3 = \{ (-\frac{1}{2},\frac{1}{2}, \frac{1}{4}),
  				      (-\frac{3}{5},-\frac{1}{2}, -\frac{1}{8}),
  				      (\frac{1}{2},-\frac{1}{2}, \frac{1}{2})
  				   \}$  。   </ul> 

有了精确的解决方案，我们就可以选择右手边的强制力和诺伊曼边界条件，从而得到这个解决方案（制造的解决方案技术）。在这个例子中，我们选择扩散等于1，对流为

\f[
\mathbf{c} = \begin{cases}
1, & \textrm{dim}=1 \\
(y, -x), & \textrm{dim}=2 \\
(y, -x, 1), & \textrm{dim}=3
\end{cases}
\f] 注意，对流是无发散的，  $\nabla \cdot c = 0$  。

<a name="Implementation"></a><h3> Implementation </h3>


除了实现上述方程，下面的实现还提供了以下功能。   <ul>   <li>  WorkStream来并行化本地求解器。在步骤9中已经详细介绍了Workstream。     <li>  从跟踪中重构本地DG解。     <li>  对解进行后处理以实现超融合。     <li>  用于直接输出全局骨架解的DataOutFaces。   </ul> 


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
 * 大多数deal.II的include文件已经在前面的例子中涉及到了，没有注释。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/tensor_function.h> 
 * #include <deal.II/base/exceptions.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/work_stream.h> 
 * #include <deal.II/base/convergence_table.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_bicgstab.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * @endcode
 * 
 * 然而，我们确实有一些新的包括在这个例子中。第一个定义了三角形面的有限元空间，我们把它称为 "骨架"。这些有限元在元素内部没有任何支持，它们代表的是在每个模数一的表面上有一个单一的值的多项式，但在模数二的表面上允许有不连续。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_face.h> 
 * 
 * @endcode
 * 
 * 我们包含的第二个新文件定义了一种新的稀疏矩阵类型。 常规的 <code>SparseMatrix</code> 类型存储了所有非零条目的索引。  <code>ChunkSparseMatrix</code> 则是利用了DG解的耦合性。 它存储了一个指定大小的矩阵子块的索引。 在HDG背景下，这个子块大小实际上是由骨架解场定义的每个面的自由度数量。这使得矩阵的内存消耗减少了三分之一，并且在求解器中使用矩阵时也会有类似的速度提升。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/chunk_sparse_matrix.h> 
 * 
 * @endcode
 * 
 * 这个例子的最后一个新的包括涉及到数据输出。 由于我们在网格的骨架上定义了一个有限元场，我们希望能够直观地看到这个解决方案的实际情况。DataOutFaces正是这样做的；它的接口与我们熟悉的DataOut几乎一样，但输出的数据只有模拟的二维1数据。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/data_out_faces.h> 
 * 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 我们首先将所有的类放入自己的命名空间。
 * 

 * 
 * 
 * @code
 * namespace Step51 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 分析解的结构与 step-7 中相同。有两个例外情况。首先，我们也为3D情况创建了一个解决方案，其次，我们对解决方案进行了缩放，使其在解决方案的所有宽度值上的规范是统一的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class SolutionBase 
 *   { 
 *   protected: 
 *     static const unsigned int n_source_centers = 3; 
 *     static const Point<dim>   source_centers[n_source_centers]; 
 *     static const double       width; 
 *   }; 
 * 
 *   template <> 
 *   const Point<1> 
 *     SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers] = 
 *       {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}; 
 * 
 *   template <> 
 *   const Point<2> 
 *     SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers] = 
 *       {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}; 
 * 
 *   template <> 
 *   const Point<3> 
 *     SolutionBase<3>::source_centers[SolutionBase<3>::n_source_centers] = { 
 *       Point<3>(-0.5, +0.5, 0.25), 
 *       Point<3>(-0.6, -0.5, -0.125), 
 *       Point<3>(+0.5, -0.5, 0.5)}; 
 * 
 *   template <int dim> 
 *   const double SolutionBase<dim>::width = 1. / 5.; 
 * 
 *   template <int dim> 
 *   class Solution : public Function<dim>, protected SolutionBase<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/ = 0) const override 
 *     { 
 *       double sum = 0; 
 *       for (unsigned int i = 0; i < this->n_source_centers; ++i) 
 *         { 
 *           const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 
 *           sum += 
 *             std::exp(-x_minus_xi.norm_square() / (this->width * this->width)); 
 *         } 
 * 
 *       return sum / 
 *              std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
 *     } 
 * 
 *     virtual Tensor<1, dim> 
 *     gradient(const Point<dim> &p, 
 *              const unsigned int /*component*/ = 0) const override 
 *     { 
 *       Tensor<1, dim> sum; 
 *       for (unsigned int i = 0; i < this->n_source_centers; ++i) 
 *         { 
 *           const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 
 * 
 *           sum += 
 *             (-2 / (this->width * this->width) * 
 *              std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) * 
 *              x_minus_xi); 
 *         } 
 * 
 *       return sum / 
 *              std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 这个类实现了一个函数，标量解和它的负梯度被收集在一起。这个函数在计算HDG近似的误差时使用，它的实现是简单地调用Solution类的值和梯度函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class SolutionAndGradient : public Function<dim>, protected SolutionBase<dim> 
 *   { 
 *   public: 
 *     SolutionAndGradient() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  v) const override 
 *     { 
 *       AssertDimension(v.size(), dim + 1); 
 *       Solution<dim>  solution; 
 *       Tensor<1, dim> grad = solution.gradient(p); 
 *       for (unsigned int d = 0; d < dim; ++d) 
 *         v[d] = -grad[d]; 
 *       v[dim] = solution.value(p); 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 接下来是对流速度的实现。如介绍中所述，我们选择的速度场在二维是 $(y, -x)$ ，在三维是 $(y, -x, 1)$ 。这就得到了一个无发散的速度场。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ConvectionVelocity : public TensorFunction<1, dim> 
 *   { 
 *   public: 
 *     ConvectionVelocity() 
 *       : TensorFunction<1, dim>() 
 *     {} 
 * 
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override 
 *     { 
 *       Tensor<1, dim> convection; 
 *       switch (dim) 
 *         { 
 *           case 1: 
 *             convection[0] = 1; 
 *             break; 
 *           case 2: 
 *             convection[0] = p[1]; 
 *             convection[1] = -p[0]; 
 *             break; 
 *           case 3: 
 *             convection[0] = p[1]; 
 *             convection[1] = -p[0]; 
 *             convection[2] = 1; 
 *             break; 
 *           default: 
 *             Assert(false, ExcNotImplemented()); 
 *         } 
 *       return convection; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 我们实现的最后一个函数是用于制造解决方案的右手边。它与 step-7 非常相似，不同的是我们现在有一个对流项而不是反应项。由于速度场是不可压缩的，即 $\nabla \cdot \mathbf{c} =0$ ，对流项简单读作  $\mathbf{c} \nabla u$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim>, protected SolutionBase<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/ = 0) const override 
 *     { 
 *       ConvectionVelocity<dim> convection_velocity; 
 *       Tensor<1, dim>          convection = convection_velocity.value(p); 
 *       double                  sum        = 0; 
 *       for (unsigned int i = 0; i < this->n_source_centers; ++i) 
 *         { 
 *           const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 
 * 
 *           sum += 
 *             ((2 * dim - 2 * convection * x_minus_xi - 
 *               4 * x_minus_xi.norm_square() / (this->width * this->width)) / 
 *              (this->width * this->width) * 
 *              std::exp(-x_minus_xi.norm_square() / (this->width * this->width))); 
 *         } 
 * 
 *       return sum / 
 *              std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="TheHDGsolverclass"></a> 
 * <h3>The HDG solver class</h3>
 * 

 * 
 * HDG的求解过程与  step-7  的求解过程非常相似。主要区别在于使用了三套不同的DoFHandler和FE对象，以及ChunkSparseMatrix和相应的解决方案向量。我们还使用WorkStream来实现多线程的本地求解过程，该过程利用了本地求解器的尴尬的并行性质。对于WorkStream，我们定义了对单元格的本地操作和复制到全局矩阵和向量的函数。我们这样做既是为了装配（装配要运行两次，一次是在我们生成系统矩阵时，另一次是在我们从骨架值计算元素内部解时），也是为了后处理，在后处理中我们提取一个在高阶收敛的解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class HDG 
 *   { 
 *   public: 
 *     enum RefinementMode 
 *     { 
 *       global_refinement, 
 *       adaptive_refinement 
 *     }; 
 * 
 *     HDG(const unsigned int degree, const RefinementMode refinement_mode); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_system(const bool reconstruct_trace = false); 
 *     void solve(); 
 *     void postprocess(); 
 *     void refine_grid(const unsigned int cycle); 
 *     void output_results(const unsigned int cycle); 
 * 
 * @endcode
 * 
 * 用于组装和解决原始变量的数据。
 * 

 * 
 * 
 * @code
 *     struct PerTaskData; 
 *     struct ScratchData; 
 * 
 * @endcode
 * 
 * 对解决方案进行后处理以获得  $u^*$  是一个逐个元素的过程；因此，我们不需要组装任何全局数据，也不需要声明任何 "任务数据 "供WorkStream使用。
 * 

 * 
 * 
 * @code
 *     struct PostProcessScratchData; 
 * 
 * @endcode
 * 
 * 以下三个函数被 WorkStream 用来完成程序的实际工作。
 * 

 * 
 * 
 * @code
 *     void assemble_system_one_cell( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       ScratchData &                                         scratch, 
 *       PerTaskData &                                         task_data); 
 * 
 *     void copy_local_to_global(const PerTaskData &data); 
 * 
 *     void postprocess_one_cell( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       PostProcessScratchData &                              scratch, 
 *       unsigned int &                                        empty_data); 
 * 
 *     Triangulation<dim> triangulation; 
 * 
 * @endcode
 * 
 * "局部 "解是每个元素的内部。 这些代表了原始解场  $u$  以及辅助场  $\mathbf{q}$  。
 * 

 * 
 * 
 * @code
 *     FESystem<dim>   fe_local; 
 *     DoFHandler<dim> dof_handler_local; 
 *     Vector<double>  solution_local; 
 * 
 * @endcode
 * 
 * 新的有限元类型和相应的 <code>DoFHandler</code> 被用于耦合元素级局部解的全局骨架解。
 * 

 * 
 * 
 * @code
 *     FE_FaceQ<dim>   fe; 
 *     DoFHandler<dim> dof_handler; 
 *     Vector<double>  solution; 
 *     Vector<double>  system_rhs; 
 * 
 * @endcode
 * 
 * 如介绍中所述，HDG解可以通过后处理达到  $\mathcal{O}(h^{p+2})$  的超收敛率。 后处理的解是一个不连续的有限元解，代表每个单元内部的原始变量。 我们定义了一个程度为 $p+1$ 的FE类型来表示这个后处理的解，我们只在构造后用于输出。
 * 

 * 
 * 
 * @code
 *     FE_DGQ<dim>     fe_u_post; 
 *     DoFHandler<dim> dof_handler_u_post; 
 *     Vector<double>  solution_u_post; 
 * 
 * @endcode
 * 
 * 与骨架相对应的自由度强烈地执行Dirichlet边界条件，就像在连续Galerkin有限元方法中一样。我们可以通过AffineConstraints对象以类似的方式强制执行边界条件。此外，悬挂节点的处理方式与连续有限元的处理方式相同。对于只在面定义自由度的面元素，这个过程将精炼面的解设置为与粗略面的表示相吻合。
 * 

 * 
 * 请注意，对于HDG来说，消除悬空节点并不是唯一的可能性，就HDG理论而言，我们也可以使用精炼侧的未知数，通过精炼侧的跟踪值来表达粗略侧的局部解。然而，这样的设置在deal.II循环方面并不容易实现，因此没有进一步分析。
 * 

 * 
 * 
 * @code
 *     AffineConstraints<double> constraints; 
 * 
 * @endcode
 * 
 * ChunkSparseMatrix类的用法与通常的稀疏矩阵类似。你需要一个ChunkSparsityPattern类型的稀疏模式和实际的矩阵对象。在创建稀疏模式时，我们只需要额外传递局部块的大小。
 * 

 * 
 * 
 * @code
 *     ChunkSparsityPattern      sparsity_pattern; 
 *     ChunkSparseMatrix<double> system_matrix; 
 * 
 * @endcode
 * 
 * 与  step-7  相同。
 * 

 * 
 * 
 * @code
 *     const RefinementMode refinement_mode; 
 *     ConvergenceTable     convergence_table; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="TheHDGclassimplementation"></a> 
 * <h3>The HDG class implementation</h3>
 * 
 * <a name=""></a> 
 * @sect4{Constructor}  该构造函数与其他例子中的构造函数类似，除了处理多个DoFHandler和FiniteElement对象。请注意，我们为局部DG部分创建了一个有限元系统，包括梯度/通量部分和标量部分。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   HDG<dim>::HDG(const unsigned int degree, const RefinementMode refinement_mode) 
 *     : fe_local(FE_DGQ<dim>(degree), dim, FE_DGQ<dim>(degree), 1) 
 *     , dof_handler_local(triangulation) 
 *     , fe(degree) 
 *     , dof_handler(triangulation) 
 *     , fe_u_post(degree + 1) 
 *     , dof_handler_u_post(triangulation) 
 *     , refinement_mode(refinement_mode) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::setup_system}  HDG解决方案的系统是以类似于其他大多数教程程序的方式设置的。 我们小心翼翼地用我们所有的DoFHandler对象来分配道夫。  @p solution 和 @p system_matrix 对象与全局骨架解决方案一起。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::setup_system() 
 *   { 
 *     dof_handler_local.distribute_dofs(fe_local); 
 *     dof_handler.distribute_dofs(fe); 
 *     dof_handler_u_post.distribute_dofs(fe_u_post); 
 * 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     solution_local.reinit(dof_handler_local.n_dofs()); 
 *     solution_u_post.reinit(dof_handler_u_post.n_dofs()); 
 * 
 *     constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *     std::map<types::boundary_id, const Function<dim> *> boundary_functions; 
 *     Solution<dim>                                       solution_function; 
 *     boundary_functions[0] = &solution_function; 
 *     VectorTools::project_boundary_values(dof_handler, 
 *                                          boundary_functions, 
 *                                          QGauss<dim - 1>(fe.degree + 1), 
 *                                          constraints); 
 *     constraints.close(); 
 * 
 * @endcode
 * 
 * 在创建块状稀疏模式时，我们首先创建通常的动态稀疏模式，然后设置块状大小，该大小等于一个面的道夫数，当把它复制到最终的稀疏模式时。
 * 

 * 
 * 
 * @code
 *     { 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
 *       sparsity_pattern.copy_from(dsp, fe.n_dofs_per_face()); 
 *     } 
 *     system_matrix.reinit(sparsity_pattern); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::PerTaskData}  接下来是定义并行装配的本地数据结构。第一个结构 @p PerTaskData 包含了被写入全局矩阵的本地向量和矩阵，而ScratchData包含了我们在本地装配中需要的所有数据。这里有一个变量值得注意，即布尔变量 @p  trace_reconstruct。正如介绍中提到的，我们分两步解决HDG系统。首先，我们为骨架系统创建一个线性系统，通过舒尔补码 $D-CA^{-1}B$  将局部部分浓缩到其中。然后，我们用骨架的解来解决局部部分。对于这两个步骤，我们需要两次元素上的相同矩阵，我们希望通过两个装配步骤来计算。由于大部分的代码是相似的，我们用相同的函数来做这件事，但只是根据我们在开始装配时设置的一个标志在两者之间切换。因为我们需要把这个信息传递给本地的工作程序，所以我们把它存储在任务数据中一次。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct HDG<dim>::PerTaskData 
 *   { 
 *     FullMatrix<double>                   cell_matrix; 
 *     Vector<double>                       cell_vector; 
 *     std::vector<types::global_dof_index> dof_indices; 
 * 
 *     bool trace_reconstruct; 
 * 
 *     PerTaskData(const unsigned int n_dofs, const bool trace_reconstruct) 
 *       : cell_matrix(n_dofs, n_dofs) 
 *       , cell_vector(n_dofs) 
 *       , dof_indices(n_dofs) 
 *       , trace_reconstruct(trace_reconstruct) 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="HDGScratchData"></a> 
 * <h4>HDG::ScratchData</h4>
 * @p ScratchData  包含WorkStream中每个线程的持久化数据。 FEValues、矩阵和矢量对象现在应该很熟悉了。 有两个对象需要讨论。  `std::vector<std::vector<unsigned  int> > fe_local_support_on_face` 和  `std::vector<std::vector<unsigned  int> > fe_support_on_face`。 这些用于指示所选择的有限元是否在与 @p fe_local 相关的局部部分和骨架部分 @p fe. 的参考单元的特定面上有支持（非零值）。 我们在构造函数中提取这一信息，并为我们工作的所有单元存储一次。 如果我们不存储这一信息，我们将被迫在每个单元上装配大量的零项，这将大大降低程序的速度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct HDG<dim>::ScratchData 
 *   { 
 *     FEValues<dim>     fe_values_local; 
 *     FEFaceValues<dim> fe_face_values_local; 
 *     FEFaceValues<dim> fe_face_values; 
 * 
 *     FullMatrix<double> ll_matrix; 
 *     FullMatrix<double> lf_matrix; 
 *     FullMatrix<double> fl_matrix; 
 *     FullMatrix<double> tmp_matrix; 
 *     Vector<double>     l_rhs; 
 *     Vector<double>     tmp_rhs; 
 * 
 *     std::vector<Tensor<1, dim>> q_phi; 
 *     std::vector<double>         q_phi_div; 
 *     std::vector<double>         u_phi; 
 *     std::vector<Tensor<1, dim>> u_phi_grad; 
 *     std::vector<double>         tr_phi; 
 *     std::vector<double>         trace_values; 
 * 
 *     std::vector<std::vector<unsigned int>> fe_local_support_on_face; 
 *     std::vector<std::vector<unsigned int>> fe_support_on_face; 
 * 
 *     ConvectionVelocity<dim> convection_velocity; 
 *     RightHandSide<dim>      right_hand_side; 
 *     const Solution<dim>     exact_solution; 
 * 
 *     ScratchData(const FiniteElement<dim> &fe, 
 *                 const FiniteElement<dim> &fe_local, 
 *                 const QGauss<dim> &       quadrature_formula, 
 *                 const QGauss<dim - 1> &   face_quadrature_formula, 
 *                 const UpdateFlags         local_flags, 
 *                 const UpdateFlags         local_face_flags, 
 *                 const UpdateFlags         flags) 
 *       : fe_values_local(fe_local, quadrature_formula, local_flags) 
 *       , fe_face_values_local(fe_local, 
 *                              face_quadrature_formula, 
 *                              local_face_flags) 
 *       , fe_face_values(fe, face_quadrature_formula, flags) 
 *       , ll_matrix(fe_local.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
 *       , lf_matrix(fe_local.n_dofs_per_cell(), fe.n_dofs_per_cell()) 
 *       , fl_matrix(fe.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
 *       , tmp_matrix(fe.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
 *       , l_rhs(fe_local.n_dofs_per_cell()) 
 *       , tmp_rhs(fe_local.n_dofs_per_cell()) 
 *       , q_phi(fe_local.n_dofs_per_cell()) 
 *       , q_phi_div(fe_local.n_dofs_per_cell()) 
 *       , u_phi(fe_local.n_dofs_per_cell()) 
 *       , u_phi_grad(fe_local.n_dofs_per_cell()) 
 *       , tr_phi(fe.n_dofs_per_cell()) 
 *       , trace_values(face_quadrature_formula.size()) 
 *       , fe_local_support_on_face(GeometryInfo<dim>::faces_per_cell) 
 *       , fe_support_on_face(GeometryInfo<dim>::faces_per_cell) 
 *       , exact_solution() 
 *     { 
 *       for (unsigned int face_no : GeometryInfo<dim>::face_indices()) 
 *         for (unsigned int i = 0; i < fe_local.n_dofs_per_cell(); ++i) 
 *           { 
 *             if (fe_local.has_support_on_face(i, face_no)) 
 *               fe_local_support_on_face[face_no].push_back(i); 
 *           } 
 * 
 *       for (unsigned int face_no : GeometryInfo<dim>::face_indices()) 
 *         for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
 *           { 
 *             if (fe.has_support_on_face(i, face_no)) 
 *               fe_support_on_face[face_no].push_back(i); 
 *           } 
 *     } 
 * 
 *     ScratchData(const ScratchData &sd) 
 *       : fe_values_local(sd.fe_values_local.get_fe(), 
 *                         sd.fe_values_local.get_quadrature(), 
 *                         sd.fe_values_local.get_update_flags()) 
 *       , fe_face_values_local(sd.fe_face_values_local.get_fe(), 
 *                              sd.fe_face_values_local.get_quadrature(), 
 *                              sd.fe_face_values_local.get_update_flags()) 
 *       , fe_face_values(sd.fe_face_values.get_fe(), 
 *                        sd.fe_face_values.get_quadrature(), 
 *                        sd.fe_face_values.get_update_flags()) 
 *       , ll_matrix(sd.ll_matrix) 
 *       , lf_matrix(sd.lf_matrix) 
 *       , fl_matrix(sd.fl_matrix) 
 *       , tmp_matrix(sd.tmp_matrix) 
 *       , l_rhs(sd.l_rhs) 
 *       , tmp_rhs(sd.tmp_rhs) 
 *       , q_phi(sd.q_phi) 
 *       , q_phi_div(sd.q_phi_div) 
 *       , u_phi(sd.u_phi) 
 *       , u_phi_grad(sd.u_phi_grad) 
 *       , tr_phi(sd.tr_phi) 
 *       , trace_values(sd.trace_values) 
 *       , fe_local_support_on_face(sd.fe_local_support_on_face) 
 *       , fe_support_on_face(sd.fe_support_on_face) 
 *       , exact_solution() 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="HDGPostProcessScratchData"></a> 
 * <h4>HDG::PostProcessScratchData</h4>
 * @p PostProcessScratchData  包含WorkStream在对本地解决方案进行后处理时使用的数据  $u^*$  。 它与  @p ScratchData.  类似，但要简单得多。
 * 
 * @code
 *   template <int dim> 
 *   struct HDG<dim>::PostProcessScratchData 
 *   { 
 *     FEValues<dim> fe_values_local; 
 *     FEValues<dim> fe_values; 
 * 
 *     std::vector<double>         u_values; 
 *     std::vector<Tensor<1, dim>> u_gradients; 
 *     FullMatrix<double>          cell_matrix; 
 * 
 *     Vector<double> cell_rhs; 
 *     Vector<double> cell_sol; 
 * 
 *     PostProcessScratchData(const FiniteElement<dim> &fe, 
 *                            const FiniteElement<dim> &fe_local, 
 *                            const QGauss<dim> &       quadrature_formula, 
 *                            const UpdateFlags         local_flags, 
 *                            const UpdateFlags         flags) 
 *       : fe_values_local(fe_local, quadrature_formula, local_flags) 
 *       , fe_values(fe, quadrature_formula, flags) 
 *       , u_values(quadrature_formula.size()) 
 *       , u_gradients(quadrature_formula.size()) 
 *       , cell_matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell()) 
 *       , cell_rhs(fe.n_dofs_per_cell()) 
 *       , cell_sol(fe.n_dofs_per_cell()) 
 *     {} 
 * 
 *     PostProcessScratchData(const PostProcessScratchData &sd) 
 *       : fe_values_local(sd.fe_values_local.get_fe(), 
 *                         sd.fe_values_local.get_quadrature(), 
 *                         sd.fe_values_local.get_update_flags()) 
 *       , fe_values(sd.fe_values.get_fe(), 
 *                   sd.fe_values.get_quadrature(), 
 *                   sd.fe_values.get_update_flags()) 
 *       , u_values(sd.u_values) 
 *       , u_gradients(sd.u_gradients) 
 *       , cell_matrix(sd.cell_matrix) 
 *       , cell_rhs(sd.cell_rhs) 
 *       , cell_sol(sd.cell_sol) 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::assemble_system}   @p assemble_system 函数与 Step-32 上的函数类似，其中正交公式和更新标志被设置，然后 <code>WorkStream</code> 被用来以多线程的方式进行工作。  @p trace_reconstruct  输入参数用于决定我们是求全局骨架解（false）还是局部解（true）。
 * 

 * 
 * 对于汇编的多线程执行，有一点值得注意的是，`assemble_system_one_cell()`中的局部计算会调用BLAS和LAPACK函数，如果这些函数在deal.II中可用。因此，底层的BLAS/LAPACK库必须支持同时来自多个线程的调用。大多数实现都支持这一点，但有些库需要以特定方式构建以避免问题。例如，在BLAS/LAPACK调用内部没有多线程的情况下编译的OpenBLAS需要在构建时将一个名为`USE_LOCKING'的标志设置为true。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::assemble_system(const bool trace_reconstruct) 
 *   { 
 *     const QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *     const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 * 
 *     const UpdateFlags local_flags(update_values | update_gradients | 
 *                                   update_JxW_values | update_quadrature_points); 
 * 
 *     const UpdateFlags local_face_flags(update_values); 
 * 
 *     const UpdateFlags flags(update_values | update_normal_vectors | 
 *                             update_quadrature_points | update_JxW_values); 
 * 
 *     PerTaskData task_data(fe.n_dofs_per_cell(), trace_reconstruct); 
 *     ScratchData scratch(fe, 
 *                         fe_local, 
 *                         quadrature_formula, 
 *                         face_quadrature_formula, 
 *                         local_flags, 
 *                         local_face_flags, 
 *                         flags); 
 * 
 *     WorkStream::run(dof_handler.begin_active(), 
 *                     dof_handler.end(), 
 *                     *this, 
 *                     &HDG<dim>::assemble_system_one_cell, 
 *                     &HDG<dim>::copy_local_to_global, 
 *                     scratch, 
 *                     task_data); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::assemble_system_one_cell}  HDG程序的实际工作由  @p assemble_system_one_cell.  组装局部矩阵  $A, B, C$  在这里完成，同时还有全局矩阵的局部贡献  $D$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::assemble_system_one_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     ScratchData &                                         scratch, 
 *     PerTaskData &                                         task_data) 
 *   { 
 * 
 * @endcode
 * 
 * 为Dof_handler_local构建迭代器，用于FEValues的reinit函数。
 * 

 * 
 * 
 * @code
 *     typename DoFHandler<dim>::active_cell_iterator loc_cell(&triangulation, 
 *                                                             cell->level(), 
 *                                                             cell->index(), 
 *                                                             &dof_handler_local); 
 * 
 *     const unsigned int n_q_points = 
 *       scratch.fe_values_local.get_quadrature().size(); 
 *     const unsigned int n_face_q_points = 
 *       scratch.fe_face_values_local.get_quadrature().size(); 
 * 
 *     const unsigned int loc_dofs_per_cell = 
 *       scratch.fe_values_local.get_fe().n_dofs_per_cell(); 
 * 
 *     const FEValuesExtractors::Vector fluxes(0); 
 *     const FEValuesExtractors::Scalar scalar(dim); 
 * 
 *     scratch.ll_matrix = 0; 
 *     scratch.l_rhs     = 0; 
 *     if (!task_data.trace_reconstruct) 
 *       { 
 *         scratch.lf_matrix     = 0; 
 *         scratch.fl_matrix     = 0; 
 *         task_data.cell_matrix = 0; 
 *         task_data.cell_vector = 0; 
 *       } 
 *     scratch.fe_values_local.reinit(loc_cell); 
 * 
 * @endcode
 * 
 * 我们首先计算对应于局部-局部耦合的 @p ll_matrix 矩阵（在介绍中称为矩阵 $A$ ）的单元内部贡献，以及局部右手向量。 我们在每个正交点存储基函数、右手边值和对流速度的值，以便快速访问这些场。
 * 

 * 
 * 
 * @code
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         const double rhs_value = scratch.right_hand_side.value( 
 *           scratch.fe_values_local.quadrature_point(q)); 
 *         const Tensor<1, dim> convection = scratch.convection_velocity.value( 
 *           scratch.fe_values_local.quadrature_point(q)); 
 *         const double JxW = scratch.fe_values_local.JxW(q); 
 *         for (unsigned int k = 0; k < loc_dofs_per_cell; ++k) 
 *           { 
 *             scratch.q_phi[k] = scratch.fe_values_local[fluxes].value(k, q); 
 *             scratch.q_phi_div[k] = 
 *               scratch.fe_values_local[fluxes].divergence(k, q); 
 *             scratch.u_phi[k] = scratch.fe_values_local[scalar].value(k, q); 
 *             scratch.u_phi_grad[k] = 
 *               scratch.fe_values_local[scalar].gradient(k, q); 
 *           } 
 *         for (unsigned int i = 0; i < loc_dofs_per_cell; ++i) 
 *           { 
 *             for (unsigned int j = 0; j < loc_dofs_per_cell; ++j) 
 *               scratch.ll_matrix(i, j) += 
 *                 (scratch.q_phi[i] * scratch.q_phi[j] - 
 *                  scratch.q_phi_div[i] * scratch.u_phi[j] + 
 *                  scratch.u_phi[i] * scratch.q_phi_div[j] - 
 *                  (scratch.u_phi_grad[i] * convection) * scratch.u_phi[j]) * 
 *                 JxW; 
 *             scratch.l_rhs(i) += scratch.u_phi[i] * rhs_value * JxW; 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 脸部条款是在所有元素的所有面上集合起来的。这与更传统的DG方法相反，在组装过程中，每个面只被访问一次。
 * 

 * 
 * 
 * @code
 *     for (const auto face_no : cell->face_indices()) 
 *       { 
 *         scratch.fe_face_values_local.reinit(loc_cell, face_no); 
 *         scratch.fe_face_values.reinit(cell, face_no); 
 * 
 * @endcode
 * 
 * 在求解局部变量时需要已经得到的  $\hat{u}$  值。
 * 

 * 
 * 
 * @code
 *         if (task_data.trace_reconstruct) 
 *           scratch.fe_face_values.get_function_values(solution, 
 *                                                      scratch.trace_values); 
 * 
 *         for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *           { 
 *             const double     JxW = scratch.fe_face_values.JxW(q); 
 *             const Point<dim> quadrature_point = 
 *               scratch.fe_face_values.quadrature_point(q); 
 *             const Tensor<1, dim> normal = 
 *               scratch.fe_face_values.normal_vector(q); 
 *             const Tensor<1, dim> convection = 
 *               scratch.convection_velocity.value(quadrature_point); 
 * 
 * @endcode
 * 
 * 这里我们计算介绍中讨论的稳定参数：由于扩散是1，并且扩散长度尺度被设定为1/5，它只是导致扩散部分的贡献为5，而对流部分的贡献是通过元素边界的居中方案中的对流大小。
 * 

 * 
 * 
 * @code
 *             const double tau_stab = (5. + std::abs(convection * normal)); 
 * 
 * @endcode
 * 
 * 我们存储非零通量和标量值，利用我们在 @p ScratchData. 中创建的 support_on_face 信息。
 * 
 * @code
 *             for (unsigned int k = 0; 
 *                  k < scratch.fe_local_support_on_face[face_no].size(); 
 *                  ++k) 
 *               { 
 *                 const unsigned int kk = 
 *                   scratch.fe_local_support_on_face[face_no][k]; 
 *                 scratch.q_phi[k] = 
 *                   scratch.fe_face_values_local[fluxes].value(kk, q); 
 *                 scratch.u_phi[k] = 
 *                   scratch.fe_face_values_local[scalar].value(kk, q); 
 *               } 
 * 
 * @endcode
 * 
 * 当  @p trace_reconstruct=false,  我们准备为骨架变量  $\hat{u}$  组装系统。如果是这种情况，我们必须组装所有与问题相关的局部矩阵：局部-局部、局部-面部、面部-局部和面部-面部。 面-面矩阵被存储为 @p TaskData::cell_matrix, ，这样就可以通过 @p copy_local_to_global将其组装到全局系统中。
 * 

 * 
 * 
 * @code
 *             if (!task_data.trace_reconstruct) 
 *               { 
 *                 for (unsigned int k = 0; 
 *                      k < scratch.fe_support_on_face[face_no].size(); 
 *                      ++k) 
 *                   scratch.tr_phi[k] = scratch.fe_face_values.shape_value( 
 *                     scratch.fe_support_on_face[face_no][k], q); 
 *                 for (unsigned int i = 0; 
 *                      i < scratch.fe_local_support_on_face[face_no].size(); 
 *                      ++i) 
 *                   for (unsigned int j = 0; 
 *                        j < scratch.fe_support_on_face[face_no].size(); 
 *                        ++j) 
 *                     { 
 *                       const unsigned int ii = 
 *                         scratch.fe_local_support_on_face[face_no][i]; 
 *                       const unsigned int jj = 
 *                         scratch.fe_support_on_face[face_no][j]; 
 *                       scratch.lf_matrix(ii, jj) += 
 *                         ((scratch.q_phi[i] * normal + 
 *                           (convection * normal - tau_stab) * scratch.u_phi[i]) * 
 *                          scratch.tr_phi[j]) * 
 *                         JxW; 
 * 
 * @endcode
 * 
 * 注意face_no-local矩阵的符号。 我们在组装时否定了这个符号，这样我们就可以在计算舒尔补时使用 FullMatrix::mmult 的加法。
 * 

 * 
 * 
 * @code
 *                       scratch.fl_matrix(jj, ii) -= 
 *                         ((scratch.q_phi[i] * normal + 
 *                           tau_stab * scratch.u_phi[i]) * 
 *                          scratch.tr_phi[j]) * 
 *                         JxW; 
 *                     } 
 * 
 *                 for (unsigned int i = 0; 
 *                      i < scratch.fe_support_on_face[face_no].size(); 
 *                      ++i) 
 *                   for (unsigned int j = 0; 
 *                        j < scratch.fe_support_on_face[face_no].size(); 
 *                        ++j) 
 *                     { 
 *                       const unsigned int ii = 
 *                         scratch.fe_support_on_face[face_no][i]; 
 *                       const unsigned int jj = 
 *                         scratch.fe_support_on_face[face_no][j]; 
 *                       task_data.cell_matrix(ii, jj) += 
 *                         ((convection * normal - tau_stab) * scratch.tr_phi[i] * 
 *                          scratch.tr_phi[j]) * 
 *                         JxW; 
 *                     } 
 * 
 *                 if (cell->face(face_no)->at_boundary() && 
 *                     (cell->face(face_no)->boundary_id() == 1)) 
 *                   { 
 *                     const double neumann_value = 
 *                       -scratch.exact_solution.gradient(quadrature_point) * 
 *                         normal + 
 *                       convection * normal * 
 *                         scratch.exact_solution.value(quadrature_point); 
 *                     for (unsigned int i = 0; 
 *                          i < scratch.fe_support_on_face[face_no].size(); 
 *                          ++i) 
 *                       { 
 *                         const unsigned int ii = 
 *                           scratch.fe_support_on_face[face_no][i]; 
 *                         task_data.cell_vector(ii) += 
 *                           scratch.tr_phi[i] * neumann_value * JxW; 
 *                       } 
 *                   } 
 *               } 
 * 
 * @endcode
 * 
 * 这最后一个项将 $\left<w,\tau u_h\right>_{\partial \mathcal T}$ 项的贡献加入到本地矩阵中。相对于上面的脸部矩阵，我们在两个装配阶段都需要它。
 * 

 * 
 * 
 * @code
 *             for (unsigned int i = 0; 
 *                  i < scratch.fe_local_support_on_face[face_no].size(); 
 *                  ++i) 
 *               for (unsigned int j = 0; 
 *                    j < scratch.fe_local_support_on_face[face_no].size(); 
 *                    ++j) 
 *                 { 
 *                   const unsigned int ii = 
 *                     scratch.fe_local_support_on_face[face_no][i]; 
 *                   const unsigned int jj = 
 *                     scratch.fe_local_support_on_face[face_no][j]; 
 *                   scratch.ll_matrix(ii, jj) += 
 *                     tau_stab * scratch.u_phi[i] * scratch.u_phi[j] * JxW; 
 *                 } 
 * 
 * @endcode
 * 
 * 当 @p trace_reconstruct=true, 时，我们在逐个元素的基础上求解局部解。 局部右手边的计算是通过用计算值 @p trace_values替换 @p 计算中的基函数 @p tr_phi。 当然，现在矩阵的符号是减号，因为我们已经把所有的东西移到了方程的另一边。
 * 

 * 
 * 
 * @code
 *             if (task_data.trace_reconstruct) 
 *               for (unsigned int i = 0; 
 *                    i < scratch.fe_local_support_on_face[face_no].size(); 
 *                    ++i) 
 *                 { 
 *                   const unsigned int ii = 
 *                     scratch.fe_local_support_on_face[face_no][i]; 
 *                   scratch.l_rhs(ii) -= 
 *                     (scratch.q_phi[i] * normal + 
 *                      scratch.u_phi[i] * (convection * normal - tau_stab)) * 
 *                     scratch.trace_values[q] * JxW; 
 *                 } 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 一旦完成所有局部贡献的组装，我们必须：（1）组装全局系统；（2）计算局部贡献。(1)组装全局系统，或者(2)计算局部解值并保存。无论哪种情况，第一步都是对局部-局部矩阵进行反转。
 * 

 * 
 * 
 * @code
 *     scratch.ll_matrix.gauss_jordan(); 
 * 
 * @endcode
 * 
 * 对于(1)，我们计算舒尔补码，并将其添加到 @p  cell_matrix，介绍中的矩阵 $D$ 。
 * 

 * 
 * 
 * @code
 *     if (task_data.trace_reconstruct == false) 
 *       { 
 *         scratch.fl_matrix.mmult(scratch.tmp_matrix, scratch.ll_matrix); 
 *         scratch.tmp_matrix.vmult_add(task_data.cell_vector, scratch.l_rhs); 
 *         scratch.tmp_matrix.mmult(task_data.cell_matrix, 
 *                                  scratch.lf_matrix, 
 *                                  true); 
 *         cell->get_dof_indices(task_data.dof_indices); 
 *       } 
 * 
 * @endcode
 * 
 * 对于(2)，我们只是求解(ll_matrix). (solution_local) = (l_rhs)。因此，我们用 @p l_rhs 乘以我们已经倒置的局部-局部矩阵，并用 <code>set_dof_values</code> 函数来存储结果。
 * 

 * 
 * 
 * @code
 *     else 
 *       { 
 *         scratch.ll_matrix.vmult(scratch.tmp_rhs, scratch.l_rhs); 
 *         loc_cell->set_dof_values(scratch.tmp_rhs, solution_local); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * 如果我们处于解题的第一步，即 @sect4{HDG::copy_local_to_global} ，那么我们就把局部矩阵组装到全局系统中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::copy_local_to_global(const PerTaskData &data) 
 *   { 
 *     if (data.trace_reconstruct == false) 
 *       constraints.distribute_local_to_global(data.cell_matrix, 
 *                                              data.cell_vector, 
 *                                              data.dof_indices, 
 *                                              system_matrix, 
 *                                              system_rhs); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::solve}  骨架解是通过使用带有身份预处理程序的BiCGStab求解器来解决的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::solve() 
 *   { 
 *     SolverControl                  solver_control(system_matrix.m() * 10, 
 *                                  1e-11 * system_rhs.l2_norm()); 
 *     SolverBicgstab<Vector<double>> solver(solver_control); 
 *     solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
 * 
 *     std::cout << "   Number of BiCGStab iterations: " 
 *               << solver_control.last_step() << std::endl; 
 * 
 *     system_matrix.clear(); 
 *     sparsity_pattern.reinit(0, 0, 0, 1); 
 * 
 *     constraints.distribute(solution); 
 * 
 * @endcode
 * 
 * 一旦我们求出了骨架解，我们就可以以逐个元素的方式求出局部解。 我们通过重新使用相同的 @p assemble_system 函数来做到这一点，但将 @p trace_reconstruct 切换为真。
 * 

 * 
 * 
 * @code
 *     assemble_system(true); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="HDGpostprocess"></a> 
 * <h4>HDG::postprocess</h4>
 * 

 * 
 * 后处理方法有两个目的。首先，我们要在度数为 $p+1$ 的元素空间中构造一个后处理的标量变量，我们希望它能在阶 $p+2$ 上收敛。这也是一个逐个元素的过程，只涉及标量解以及局部单元上的梯度。为了做到这一点，我们引入了已经定义好的从头开始的数据以及一些更新标志，并运行工作流来并行地完成这一工作。
 * 

 * 
 * 第二，我们要计算离散化误差，就像我们在  step-7  中做的那样。整个过程与调用 VectorTools::integrate_difference. 相似，区别在于我们如何计算标量变量和梯度变量的误差。在 step-7 中，我们通过计算 @p L2_norm 或 @p H1_seminorm 的贡献来做到这一点。在这里，我们有一个DoFHandler，计算了这两个贡献，并按其矢量分量排序， <code>[0, dim)</code> 为梯度， @p dim 为标量。为了计算它们的值，我们用一个ComponentSelectFunction来计算它们中的任何一个，再加上上面介绍的 @p SolutionAndGradient类，它包含了它们中任何一个的分析部分。最终，我们还计算了后处理的解决方案的L2-误差，并将结果添加到收敛表中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::postprocess() 
 *   { 
 *     { 
 *       const QGauss<dim> quadrature_formula(fe_u_post.degree + 1); 
 *       const UpdateFlags local_flags(update_values); 
 *       const UpdateFlags flags(update_values | update_gradients | 
 *                               update_JxW_values); 
 * 
 *       PostProcessScratchData scratch( 
 *         fe_u_post, fe_local, quadrature_formula, local_flags, flags); 
 * 
 *       WorkStream::run( 
 *         dof_handler_u_post.begin_active(), 
 *         dof_handler_u_post.end(), 
 *         [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *                PostProcessScratchData &                              scratch, 
 *                unsigned int &                                        data) { 
 *           this->postprocess_one_cell(cell, scratch, data); 
 *         }, 
 *         std::function<void(const unsigned int &)>(), 
 *         scratch, 
 *         0U); 
 *     } 
 * 
 *     Vector<float> difference_per_cell(triangulation.n_active_cells()); 
 * 
 *     ComponentSelectFunction<dim> value_select(dim, dim + 1); 
 *     VectorTools::integrate_difference(dof_handler_local, 
 *                                       solution_local, 
 *                                       SolutionAndGradient<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe.degree + 2), 
 *                                       VectorTools::L2_norm, 
 *                                       &value_select); 
 *     const double L2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 *     ComponentSelectFunction<dim> gradient_select( 
 *       std::pair<unsigned int, unsigned int>(0, dim), dim + 1); 
 *     VectorTools::integrate_difference(dof_handler_local, 
 *                                       solution_local, 
 *                                       SolutionAndGradient<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe.degree + 2), 
 *                                       VectorTools::L2_norm, 
 *                                       &gradient_select); 
 *     const double grad_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 *     VectorTools::integrate_difference(dof_handler_u_post, 
 *                                       solution_u_post, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe.degree + 3), 
 *                                       VectorTools::L2_norm); 
 *     const double post_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 *     convergence_table.add_value("cells", triangulation.n_active_cells()); 
 *     convergence_table.add_value("dofs", dof_handler.n_dofs()); 
 * 
 *     convergence_table.add_value("val L2", L2_error); 
 *     convergence_table.set_scientific("val L2", true); 
 *     convergence_table.set_precision("val L2", 3); 
 * 
 *     convergence_table.add_value("grad L2", grad_error); 
 *     convergence_table.set_scientific("grad L2", true); 
 *     convergence_table.set_precision("grad L2", 3); 
 * 
 *     convergence_table.add_value("val L2-post", post_error); 
 *     convergence_table.set_scientific("val L2-post", true); 
 *     convergence_table.set_precision("val L2-post", 3); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="HDGpostprocess_one_cell"></a> 
 * <h4>HDG::postprocess_one_cell</h4>
 * 

 * 
 * 这是为后处理所做的实际工作。根据介绍中的讨论，我们需要建立一个系统，将DG解的梯度部分投影到后处理变量的梯度上。此外，我们还需要将新的后处理变量的平均值设置为等于标量DG解在单元上的平均值。
 * 

 * 
 * 从技术上讲，梯度的投影是一个有可能填满我们的 @p dofs_per_cell 乘以 @p dofs_per_cell 矩阵的系统，但它是单数（所有行的总和为零，因为常数函数的梯度为零）。因此，我们拿掉一行，用它来强加标量值的平均值。我们为标量部分挑选第一行，尽管我们可以为 $\mathcal Q_{-p}$ 元素挑选任何一行。然而，如果我们使用FE_DGP元素，第一行将对应常数部分，删除例如最后一行将得到一个奇异系统。这样一来，我们的程序也可以用于这些元素。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::postprocess_one_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     PostProcessScratchData &                              scratch, 
 *     unsigned int &) 
 *   { 
 *     typename DoFHandler<dim>::active_cell_iterator loc_cell(&triangulation, 
 *                                                             cell->level(), 
 *                                                             cell->index(), 
 *                                                             &dof_handler_local); 
 * 
 *     scratch.fe_values_local.reinit(loc_cell); 
 *     scratch.fe_values.reinit(cell); 
 * 
 *     FEValuesExtractors::Vector fluxes(0); 
 *     FEValuesExtractors::Scalar scalar(dim); 
 * 
 *     const unsigned int n_q_points = scratch.fe_values.get_quadrature().size(); 
 *     const unsigned int dofs_per_cell = scratch.fe_values.dofs_per_cell; 
 * 
 *     scratch.fe_values_local[scalar].get_function_values(solution_local, 
 *                                                         scratch.u_values); 
 *     scratch.fe_values_local[fluxes].get_function_values(solution_local, 
 *                                                         scratch.u_gradients); 
 * 
 *     double sum = 0; 
 *     for (unsigned int i = 1; i < dofs_per_cell; ++i) 
 *       { 
 *         for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *           { 
 *             sum = 0; 
 *             for (unsigned int q = 0; q < n_q_points; ++q) 
 *               sum += (scratch.fe_values.shape_grad(i, q) * 
 *                       scratch.fe_values.shape_grad(j, q)) * 
 *                      scratch.fe_values.JxW(q); 
 *             scratch.cell_matrix(i, j) = sum; 
 *           } 
 * 
 *         sum = 0; 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           sum -= (scratch.fe_values.shape_grad(i, q) * scratch.u_gradients[q]) * 
 *                  scratch.fe_values.JxW(q); 
 *         scratch.cell_rhs(i) = sum; 
 *       } 
 *     for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *       { 
 *         sum = 0; 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           sum += scratch.fe_values.shape_value(j, q) * scratch.fe_values.JxW(q); 
 *         scratch.cell_matrix(0, j) = sum; 
 *       } 
 *     { 
 *       sum = 0; 
 *       for (unsigned int q = 0; q < n_q_points; ++q) 
 *         sum += scratch.u_values[q] * scratch.fe_values.JxW(q); 
 *       scratch.cell_rhs(0) = sum; 
 *     } 
 * 
 * @endcode
 * 
 * 集合了所有条款后，我们又可以继续解决这个线性系统。我们对矩阵进行反转，然后将反转结果乘以右手边。另一种方法（数字上更稳定）是只对矩阵进行因式分解，然后应用因式分解。
 * 

 * 
 * 
 * @code
 *     scratch.cell_matrix.gauss_jordan(); 
 *     scratch.cell_matrix.vmult(scratch.cell_sol, scratch.cell_rhs); 
 *     cell->distribute_local_to_global(scratch.cell_sol, solution_u_post); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::output_results}  我们有三组我们想输出的结果：局部解决方案，后处理的局部解决方案，以及骨架解决方案。前两个结果都 "活 "在元素体积上，而后者则活在三角形的一维表面上。 我们的 @p output_results 函数将所有的局部解决方案写入同一个vtk文件，尽管它们对应于不同的DoFHandler对象。 骨架变量的图形输出是通过使用DataOutFaces类完成的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::output_results(const unsigned int cycle) 
 *   { 
 *     std::string filename; 
 *     switch (refinement_mode) 
 *       { 
 *         case global_refinement: 
 *           filename = "solution-global"; 
 *           break; 
 *         case adaptive_refinement: 
 *           filename = "solution-adaptive"; 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     std::string face_out(filename); 
 *     face_out += "-face"; 
 * 
 *     filename += "-q" + Utilities::int_to_string(fe.degree, 1); 
 *     filename += "-" + Utilities::int_to_string(cycle, 2); 
 *     filename += ".vtk"; 
 *     std::ofstream output(filename); 
 * 
 *     DataOut<dim> data_out; 
 * 
 * @endcode
 * 
 * 我们首先定义本地解决方案的名称和类型，并将数据添加到  @p data_out.  中。
 * 
 * @code
 *     std::vector<std::string> names(dim, "gradient"); 
 *     names.emplace_back("solution"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       component_interpretation( 
 *         dim + 1, DataComponentInterpretation::component_is_part_of_vector); 
 *     component_interpretation[dim] = 
 *       DataComponentInterpretation::component_is_scalar; 
 *     data_out.add_data_vector(dof_handler_local, 
 *                              solution_local, 
 *                              names, 
 *                              component_interpretation); 
 * 
 * @endcode
 * 
 * 我们添加的第二个数据项是后处理的解决方案。在这种情况下，它是一个属于不同DoFHandler的单一标量变量。
 * 

 * 
 * 
 * @code
 *     std::vector<std::string> post_name(1, "u_post"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       post_comp_type(1, DataComponentInterpretation::component_is_scalar); 
 *     data_out.add_data_vector(dof_handler_u_post, 
 *                              solution_u_post, 
 *                              post_name, 
 *                              post_comp_type); 
 * 
 *     data_out.build_patches(fe.degree); 
 *     data_out.write_vtk(output); 
 * 
 *     face_out += "-q" + Utilities::int_to_string(fe.degree, 1); 
 *     face_out += "-" + Utilities::int_to_string(cycle, 2); 
 *     face_out += ".vtk"; 
 *     std::ofstream face_output(face_out); 
 * 
 * @endcode
 * 
 * <code>DataOutFaces</code> 类的工作原理与 <code>DataOut</code> class when we have a <code>DoFHandler</code> 类似，后者定义了三角形骨架上的解决方案。 我们在此将其视为如此，代码与上面类似。
 * 

 * 
 * 
 * @code
 *     DataOutFaces<dim>        data_out_face(false); 
 *     std::vector<std::string> face_name(1, "u_hat"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       face_component_type(1, DataComponentInterpretation::component_is_scalar); 
 * 
 *     data_out_face.add_data_vector(dof_handler, 
 *                                   solution, 
 *                                   face_name, 
 *                                   face_component_type); 
 * 
 *     data_out_face.build_patches(fe.degree); 
 *     data_out_face.write_vtk(face_output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="HDGrefine_grid"></a> 
 * <h4>HDG::refine_grid</h4>
 * 

 * 
 * 我们为HDG实现了两种不同的细化情况，就像在 <code>Step-7</code> 中一样：adaptive_refinement和global_refinement。 global_refinement选项每次都会重新创建整个三角形。这是因为我们想使用比一个细化步骤更细的网格序列，即每个方向2、3、4、6、8、12、16...个元素。
 * 

 * 
 * adaptive_refinement模式使用 <code>KellyErrorEstimator</code> 对标量局部解中的非规则区域给出一个体面的指示。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::refine_grid(const unsigned int cycle) 
 *   { 
 *     if (cycle == 0) 
 *       { 
 *         GridGenerator::subdivided_hyper_cube(triangulation, 2, -1, 1); 
 *         triangulation.refine_global(3 - dim); 
 *       } 
 *     else 
 *       switch (refinement_mode) 
 *         { 
 *           case global_refinement: 
 *             { 
 *               triangulation.clear(); 
 *               GridGenerator::subdivided_hyper_cube(triangulation, 
 *                                                    2 + (cycle % 2), 
 *                                                    -1, 
 *                                                    1); 
 *               triangulation.refine_global(3 - dim + cycle / 2); 
 *               break; 
 *             } 
 * 
 *           case adaptive_refinement: 
 *             { 
 *               Vector<float> estimated_error_per_cell( 
 *                 triangulation.n_active_cells()); 
 * 
 *               FEValuesExtractors::Scalar scalar(dim); 
 *               std::map<types::boundary_id, const Function<dim> *> 
 *                 neumann_boundary; 
 *               KellyErrorEstimator<dim>::estimate(dof_handler_local, 
 *                                                  QGauss<dim - 1>(fe.degree + 1), 
 *                                                  neumann_boundary, 
 *                                                  solution_local, 
 *                                                  estimated_error_per_cell, 
 *                                                  fe_local.component_mask( 
 *                                                    scalar)); 
 * 
 *               GridRefinement::refine_and_coarsen_fixed_number( 
 *                 triangulation, estimated_error_per_cell, 0.3, 0.); 
 * 
 *               triangulation.execute_coarsening_and_refinement(); 
 * 
 *               break; 
 *             } 
 * 
 *           default: 
 *             { 
 *               Assert(false, ExcNotImplemented()); 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 就像在 step-7 中一样，我们将其中两个面的边界指标设置为1，在这里我们要指定诺伊曼边界条件而不是迪里希特条件。由于我们每次都会为全局细化重新创建三角形，所以在每个细化步骤中都会设置标志，而不仅仅是在开始时。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : triangulation.cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->at_boundary()) 
 *           if ((std::fabs(face->center()(0) - (-1)) < 1e-12) || 
 *               (std::fabs(face->center()(1) - (-1)) < 1e-12)) 
 *             face->set_boundary_id(1); 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{HDG::run}  这里的功能与 <code>Step-7</code>  基本相同。我们在10个周期中循环，在每个周期中细化网格。 在最后，收敛表被创建。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HDG<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 10; ++cycle) 
 *       { 
 *         std::cout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         refine_grid(cycle); 
 *         setup_system(); 
 *         assemble_system(false); 
 *         solve(); 
 *         postprocess(); 
 *         output_results(cycle); 
 *       } 
 * 
 * @endcode
 * 
 * 与 step-7 相比，收敛表有一个微小的变化：由于我们没有在每个周期内以2的系数细化我们的网格（而是使用2，3，4，6，8，12，...的序列），我们需要告诉收敛率评估这一点。我们通过设置单元格数量作为参考列，并额外指定问题的维度来实现这一目的，这为单元格数量和网格大小之间的关系提供了必要的信息。
 * 

 * 
 * 
 * @code
 *     if (refinement_mode == global_refinement) 
 *       { 
 *         convergence_table.evaluate_convergence_rates( 
 *           "val L2", "cells", ConvergenceTable::reduction_rate_log2, dim); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "grad L2", "cells", ConvergenceTable::reduction_rate_log2, dim); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "val L2-post", "cells", ConvergenceTable::reduction_rate_log2, dim); 
 *       } 
 *     convergence_table.write_text(std::cout); 
 *   } 
 * 
 * } // end of namespace Step51 
 * 
 * int main() 
 * { 
 *   const unsigned int dim = 2; 
 * 
 *   try 
 *     { 
 * 
 * @endcode
 * 
 * 现在是对主类的三次调用，完全类似于  step-7  。
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
 *         Step51::HDG<dim> hdg_problem(1, Step51::HDG<dim>::adaptive_refinement); 
 *         hdg_problem.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 * 
 *       { 
 *         std::cout << "Solving with Q1 elements, global refinement" << std::endl 
 *                   << "===========================================" << std::endl 
 *                   << std::endl; 
 * 
 *         Step51::HDG<dim> hdg_problem(1, Step51::HDG<dim>::global_refinement); 
 *         hdg_problem.run(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 * 
 *       { 
 *         std::cout << "Solving with Q3 elements, global refinement" << std::endl 
 *                   << "===========================================" << std::endl 
 *                   << std::endl; 
 * 
 *         Step51::HDG<dim> hdg_problem(3, Step51::HDG<dim>::global_refinement); 
 *         hdg_problem.run(); 
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
examples/step-51/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Programoutput"></a><h3>Program output</h3>


我们首先看一下程序在二维运行时产生的输出。在下面的四张图片中，我们展示了多项式度数 $p=1$ 和程序的周期2、3、4和8的解决方案。在图中，我们将从内部数据（DG部分）产生的数据与骨架部分（ $\hat{u}$ ）叠加到同一图中。我们不得不生成两个不同的数据集，因为单元格和面孔代表不同的几何实体，它们的组合（在同一个文件中）在VTK输出的deal.II中不被支持。

这些图像显示了HDG的明显特征：细胞的解决方案（彩色的表面）在细胞之间是不连续的。骨架变量上的解位于面的位置，并将局部部分联系起来。骨架解决方案在面与面之间的顶点上是不连续的，尽管它的值沿着同一坐标方向的线相当接近。骨架解可以被解释为两边之间的橡胶弹簧，它可以平衡解的跳跃（或者说，通量 $\kappa \nabla u
+ \mathbf{c} u$ ）。从左上角的图片可以看出，大体解经常出现过冲和欠冲，而骨架变量确实是对精确解更好的近似；这解释了为什么我们可以通过后处理步骤得到更好的解。

随着网格的细化，单元之间的跳跃变得很小（我们代表一个平滑的解决方案），骨架解决方案接近内部部分。对于第8周期，两个变量没有明显的区别。我们还看到边界条件是如何弱化实施的，内部变量并不完全满足边界条件。在下部和左侧边界，我们设置了诺伊曼边界条件，而在右侧和顶部边界，我们设置了迪里希特条件。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.sol_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.sol_3.png" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.sol_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.sol_8.png" alt=""></td>
  </tr>
</table> 

接下来，我们看一下后处理的解决方案，还是在周期2、3、4和8。这是一个不连续的解决方案，局部由二阶多项式描述。虽然在第二周期的网格上，解决方案看起来不是很好，但在第三和第四周期，它看起来好得多。正如下面的收敛表所示，我们发现它也更快地收敛到了分析解。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.post_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.post_3.png" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.post_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.post_8.png" alt=""></td>
  </tr>
</table> 

最后，我们看一下 $p=3$ 在第二周期的解。尽管网格较粗，只有64个单元，但经过后处理的解在质量上与第8周期4,096个单元的线性解（未经过后处理）相似。这清楚地表明了高阶方法对于平滑解的优越性。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.sol_q3_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.post_q3_2.png" alt=""></td>
  </tr>
</table> 

<a name="Convergencetables"></a><h4>Convergence tables</h4>


当程序运行时，它还会输出各自的步骤和收敛表的信息，并在最后输出各部分的误差。在2D中，收敛表看起来如下。

@code
Q1 elements, adaptive refinement:
cells dofs   val L2    grad L2  val L2-post
   16    80 1.804e+01 2.207e+01   1.798e+01
   31   170 9.874e+00 1.322e+01   9.798e+00
   61   314 7.452e-01 3.793e+00   4.891e-01
  121   634 3.240e-01 1.511e+00   2.616e-01
  238  1198 8.585e-02 8.212e-01   1.808e-02
  454  2290 4.802e-02 5.178e-01   2.195e-02
  898  4378 2.561e-02 2.947e-01   4.318e-03
 1720  7864 1.306e-02 1.664e-01   2.978e-03
 3271 14638 7.025e-03 9.815e-02   1.075e-03
 6217 27214 4.119e-03 6.407e-02   9.975e-04


Q1 elements, global refinement:
cells dofs      val L2        grad L2      val L2-post
   16    80 1.804e+01    - 2.207e+01    - 1.798e+01    -
   36   168 6.125e+00 2.66 9.472e+00 2.09 6.084e+00 2.67
   64   288 9.785e-01 6.38 4.260e+00 2.78 7.102e-01 7.47
  144   624 2.730e-01 3.15 1.866e+00 2.04 6.115e-02 6.05
  256  1088 1.493e-01 2.10 1.046e+00 2.01 2.880e-02 2.62
  576  2400 6.965e-02 1.88 4.846e-01 1.90 9.204e-03 2.81
 1024  4224 4.018e-02 1.91 2.784e-01 1.93 4.027e-03 2.87
 2304  9408 1.831e-02 1.94 1.264e-01 1.95 1.236e-03 2.91
 4096 16640 1.043e-02 1.96 7.185e-02 1.96 5.306e-04 2.94
 9216 37248 4.690e-03 1.97 3.228e-02 1.97 1.599e-04 2.96


Q3 elements, global refinement:
cells dofs      val L2        grad L2      val L2-post
   16   160 3.613e-01    - 1.891e+00    - 3.020e-01    -
   36   336 6.411e-02 4.26 5.081e-01 3.24 3.238e-02 5.51
   64   576 3.480e-02 2.12 2.533e-01 2.42 5.277e-03 6.31
  144  1248 8.297e-03 3.54 5.924e-02 3.58 6.330e-04 5.23
  256  2176 2.254e-03 4.53 1.636e-02 4.47 1.403e-04 5.24
  576  4800 4.558e-04 3.94 3.277e-03 3.96 1.844e-05 5.01
 1024  8448 1.471e-04 3.93 1.052e-03 3.95 4.378e-06 5.00
 2304 18816 2.956e-05 3.96 2.104e-04 3.97 5.750e-07 5.01
 4096 33280 9.428e-06 3.97 6.697e-05 3.98 1.362e-07 5.01
 9216 74496 1.876e-06 3.98 1.330e-05 3.99 1.788e-08 5.01
@endcode




我们可以看到网格细化后的误差减少，对于进行全局细化的情况，也可以看到收敛率。在 $L_2$ 准则下，标量变量和梯度变量的Q1元素的二次收敛率很明显，在 $L_2$ 准则下，后处理的标量变量的三次收敛率也是如此。注意HDG解决方案的这一明显特征。在典型的连续有限元中，阶 $p$ 的解的梯度收敛率只有 $p$ ，与实际解的 $p+1$ 相反。即使有限元的超收敛结果也是可用的（例如Zienkiewicz和Zhu首次提出的超收敛补丁恢复），但这些通常只限于结构化网格和其他特殊情况。对于Q3 HDG变量，标量变量和梯度在四阶收敛，后处理的标量变量在五阶收敛。

在3D中观察到相同的收敛率。

@code
Q1 elements, adaptive refinement:
cells   dofs    val L2    grad L2  val L2-post
     8     144 7.122e+00 1.941e+01   6.102e+00
    29     500 3.309e+00 1.023e+01   2.145e+00
   113    1792 2.204e+00 1.023e+01   1.912e+00
   379    5732 6.085e-01 5.008e+00   2.233e-01
  1317   19412 1.543e-01 1.464e+00   4.196e-02
  4579   64768 5.058e-02 5.611e-01   9.521e-03
 14596  199552 2.129e-02 3.122e-01   4.569e-03
 46180  611400 1.033e-02 1.622e-01   1.684e-03
144859 1864212 5.007e-03 8.371e-02   7.364e-04
451060 5684508 2.518e-03 4.562e-02   3.070e-04


Q1 elements, global refinement:
cells   dofs       val L2          grad L2       val L2-post
     8     144 7.122e+00    - 1.941e+01     - 6.102e+00    -
    27     432 5.491e+00 0.64 2.184e+01 -0.29 4.448e+00 0.78
    64     960 3.646e+00 1.42 1.299e+01  1.81 3.306e+00 1.03
   216    3024 1.595e+00 2.04 8.550e+00  1.03 1.441e+00 2.05
   512    6912 6.922e-01 2.90 5.306e+00  1.66 2.511e-01 6.07
  1728   22464 2.915e-01 2.13 2.490e+00  1.87 8.588e-02 2.65
  4096   52224 1.684e-01 1.91 1.453e+00  1.87 4.055e-02 2.61
 13824  172800 7.972e-02 1.84 6.861e-01  1.85 1.335e-02 2.74
 32768  405504 4.637e-02 1.88 3.984e-01  1.89 5.932e-03 2.82
110592 1354752 2.133e-02 1.92 1.830e-01  1.92 1.851e-03 2.87


Q3 elements, global refinement:
cells   dofs       val L2        grad L2      val L2-post
     8     576 5.670e+00    - 1.868e+01    - 5.462e+00    -
    27    1728 1.048e+00 4.16 6.988e+00 2.42 8.011e-01 4.73
    64    3840 2.831e-01 4.55 2.710e+00 3.29 1.363e-01 6.16
   216   12096 7.883e-02 3.15 7.721e-01 3.10 2.158e-02 4.55
   512   27648 3.642e-02 2.68 3.305e-01 2.95 5.231e-03 4.93
  1728   89856 8.546e-03 3.58 7.581e-02 3.63 7.640e-04 4.74
  4096  208896 2.598e-03 4.14 2.313e-02 4.13 1.783e-04 5.06
 13824  691200 5.314e-04 3.91 4.697e-03 3.93 2.355e-05 4.99
 32768 1622016 1.723e-04 3.91 1.517e-03 3.93 5.602e-06 4.99
110592 5419008 3.482e-05 3.94 3.055e-04 3.95 7.374e-07 5.00
@endcode



<a name="Comparisonwithcontinuousfiniteelements"></a><h3>Comparison with continuous finite elements</h3>


<a name="Resultsfor2D"></a><h4>Results for 2D</h4>


收敛表验证了介绍中所述的预期收敛率。现在，我们想在本教程的问题上展示一下HDG方法与普通有限元（连续Galkerin）方法相比的计算效率的快速比较。当然，与连续有限元相比，HDG方法对于传输为主的问题的稳定性方面在实践中也很重要，这是在平滑分析解的问题上看不到的方面。在下面的图片中，我们比较了 $L_2$ 误差作为自由度数的函数（左）和线性求解器中花费的计算时间（右），连续有限元（CG）和本教程中介绍的混合非连续Galerkin方法的两个空间维度。相对于教程中我们只使用无条件的BiCGStab，下面的数字中显示的时间使用了 TrilinosWrappers::PreconditionAMG. 中的Trilinos代数多网格预处理器 对于HDG部分，为了利用最细级别的矩阵中的块结构，我们使用了ChunkSparseMatrix周围的跟踪变量的包装器。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2d_plain.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2dt_plain.png" width="400" alt=""></td>
  </tr>
</table> 

图中的结果显示，HDG方法在 $p=1$ 时比连续有限元慢，对立方体元素的速度差不多，对六阶元素的速度更快。然而，我们在上面已经看到，HDG方法实际上产生的解比原始变量所表示的更准确。因此，在下面两幅图中，我们转而显示HDG的后处理解的误差（例如用 $p=1^*$ 表示）。现在我们看到，对于 $p=3$ 和 $p=6$ ，在相同的工作量下，HDG有明显的优势，而对于 $p=1$ ，质量也差不多。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2d_post.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2dt_post.png" width="400" alt=""></td>
  </tr>
</table> 

由于HDG方法实际产生的结果收敛为 $h^{p+2}$ ，我们应该将其与具有相同渐进收敛行为的连续Galerkin解决方案进行比较，即程度为 $p+1$ 的FE_Q。如果我们这样做，我们会得到下面的收敛曲线。我们看到，用二阶多项式的CG又明显优于用线型的HDG。然而，HDG对高阶的优势依然存在。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2d_postb.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.2dt_postb.png" width="400" alt=""></td>
  </tr>
</table> 

这些结果与一般的DG方法的特性是一致的。最佳性能通常不是在线性元素上实现的，而是在更高的阶数上，通常在 $p=3$ 左右。这是因为不连续解的体积-表面效应，有太多的解存在于表面，因此当元素为线性时，会出现重复工作。换句话说，尽管DG方法专注于不连续（因此看起来精度不高）的解的表示，但在相对高的阶数上使用时往往是最有效的。

<a name="Resultsfor3D"></a><h4>Results for 3D</h4>


我们现在展示了同样的三维数字：第一行显示了自由度数和计算时间与标量变量 $L_2$ 中的 $u$ 误差的关系，在 $p$ 阶的CG和HDG，第二行显示了后处理的HDG方案，而不是原始方案，第三行比较了后处理的HDG方案与阶 $p+1$ 的CG。在三维中，体积-表面效应使得HDG的成本更高，对于线型来说，CG的解决方案显然比HDG更好。对于立方体，HDG和CG的质量相似，而HDG对于六阶多项式又更有效率。我们也可以使用FE_DGP和FE_FaceP的组合来代替（FE_DGQ, FE_FaceQ），它们不使用 $p$ 度的张量积多项式，而是<i>complete</i>度的Legendre多项式。在给定的网格尺寸下，FE_FaceP的骨架变量的自由度较少，但求解质量（误差与自由度数量的关系）与FE_FaceQ的结果非常相似。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3d_plain.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3dt_plain.png" width="400" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3d_post.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3dt_post.png" width="400" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3d_postb.png" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-51.3dt_postb.png" width="400" alt=""></td>
  </tr>
</table> 

关于效率比较的最后一点说明。我们试图使用通用的稀疏矩阵结构和类似的求解器（两者的最佳AMG预处理器，没有对任何一个AMG参数进行特别的调整），在一个玩具的例子上，对两种方法的成本与精度进行公平的描述。然而，应该注意的是，连续有限元的几何多网格（GMG）对于 $p=3$ 和 $p=6$ 来说，大约快了四到五个系数。截至2019年，HDG的最优复杂度迭代求解器仍在研究界开发中。另外，还有其他方面的CG的实现，如步骤37所示的快速无矩阵方法，使高阶连续元素更具竞争力。同样，本教程的作者也不清楚是否可以对HDG做出类似的改进。我们参考<a href="https://dx.doi.org/10.1137/16M110455X">Kronbichler
and Wall (2018)</a>，了解最近的效率评估。




<a name="Possibilitiesforimprovements"></a><h3>Possibilities for improvements</h3>


正如在介绍中已经提到的，一种可能性是实施文献中讨论的另一种后处理技术。

第二项没有做得很好的是与这个程序的性能有关，这当然是实际应用中的一个问题（也要权衡(H)DG方法对传输为主的问题有更好的解决质量）。让我们来看看这个教程程序的计算时间和各个部分的份额。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th>&nbsp;</th>
    <th>Setup</th>
    <th>Assemble</th>
    <th>Solve</th>
    <th>Trace reconstruct</th>
    <th>Post-processing</th>
    <th>Output</th>
  </tr>
  <tr>
    <th>&nbsp;</th>
    <th>Total time</th>
    <th colspan="6">Relative share</th>
  </tr>
  <tr>
    <td align="left">2D, Q1, cycle 9, 37,248 dofs</td>
    <td align="center">5.34s</td>
    <td align="center">0.7%</td>
    <td align="center">1.2%</td>
    <td align="center">89.5%</td>
    <td align="center">0.9%</td>
    <td align="center">2.3%</td>
    <td align="center">5.4%</td>
  </tr>
  <tr>
    <td align="left">2D, Q3, cycle 9, 74,496 dofs</td>
    <td align="center">22.2s</td>
    <td align="center">0.4%</td>
    <td align="center">4.3%</td>
    <td align="center">84.1%</td>
    <td align="center">4.1%</td>
    <td align="center">3.5%</td>
    <td align="center">3.6%</td>
  </tr>
  <tr>
    <td align="left">3D, Q1, cycle 7, 172,800 dofs</td>
    <td align="center">9.06s</td>
    <td align="center">3.1%</td>
    <td align="center">8.9%</td>
    <td align="center">42.7%</td>
    <td align="center">7.0%</td>
    <td align="center">20.6%</td>
    <td align="center">17.7%</td>
  </tr>
  <tr>
    <td align="left">3D, Q3, cycle 7, 691,200 dofs</td>
    <td align="center">516s</td>
    <td align="center">0.6%</td>
    <td align="center">34.5%</td>
    <td align="center">13.4%</td>
    <td align="center">32.8%</td>
    <td align="center">17.1%</td>
    <td align="center">1.5%</td>
  </tr>
</table> 

从表中可以看出，解算器和汇编调用在程序的运行时间中占主导地位。这也清楚地表明，在哪些方面的改进是最有意义的。

<ol>  <li>  更好的线性求解器。我们使用的是BiCGStab迭代求解器，没有预处理程序，迭代次数随着问题大小的增加而增加（Q1元素和全局细化的迭代次数在小尺寸时从35次开始，但在最大尺寸时增加到701次）。为了做得更好，例如可以使用Trilinos的代数多网格预处理程序，或者像<a
  href="https://dx.doi.org/10.1137/16M110455X">Kronbichler and Wall
  (2018)</a>中讨论的一些更高级的变体。对于以扩散为主的问题，比如目前的问题，只要我们不与MPI并行工作，就可以设计这样一个求解器，在最细的层次上使用更高效的ChunkSparseMatrix的矩阵-向量积。对于MPI并行化的计算，可以使用一个标准 TrilinosWrappers::SparseMatrix 。

    <li>  通过预先组装那些不从一个单元改变到另一个单元的部件（那些既不包含可变系数也不包含依赖映射的项）来加快组装速度。   </ol> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-51.cc"
*/
