/**
@page step_61 The step-61 tutorial program
This tutorial depends on step-51.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#WeakGalerkinfiniteelementmethods"> Weak Galerkin finite element methods </a>
        <li><a href="#Theequationtosolve"> The equation to solve </a>
        <li><a href="#WeakGalerkinscheme"> Weak Galerkin scheme </a>
        <li><a href="#Representingtheweakgradient"> Representing the weak gradient </a>
        <li><a href="#Assemblingthelinearsystem"> Assembling the linear system </a>
        <li><a href="#PostprocessingandiLsub2subiiLsub2subierrors"> Post-processing and <i>L<sub>2</sub></i><i>L<sub>2</sub></i>-errors </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#TheWGDarcyEquationclasstemplate">The WGDarcyEquation class template</a>
        <li><a href="#Righthandsideboundaryvaluesandexactsolution">Right hand side, boundary values, and exact solution</a>
        <li><a href="#WGDarcyEquationclassimplementation">WGDarcyEquation class implementation</a>
      <ul>
        <li><a href="#WGDarcyEquationWGDarcyEquation">WGDarcyEquation::WGDarcyEquation</a>
        <li><a href="#WGDarcyEquationmake_grid">WGDarcyEquation::make_grid</a>
        <li><a href="#WGDarcyEquationsetup_system">WGDarcyEquation::setup_system</a>
        <li><a href="#WGDarcyEquationassemble_system">WGDarcyEquation::assemble_system</a>
        <li><a href="#WGDarcyEquationdimsolve">WGDarcyEquation<dim>::solve</a>
        <li><a href="#WGDarcyEquationdimcompute_postprocessed_velocity">WGDarcyEquation<dim>::compute_postprocessed_velocity</a>
        <li><a href="#WGDarcyEquationdimcompute_pressure_error">WGDarcyEquation<dim>::compute_pressure_error</a>
        <li><a href="#WGDarcyEquationdimcompute_velocity_error">WGDarcyEquation<dim>::compute_velocity_error</a>
        <li><a href="#WGDarcyEquationoutput_results">WGDarcyEquation::output_results</a>
        <li><a href="#WGDarcyEquationrun">WGDarcyEquation::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestresultsoniWGQsub0subQsub0subRTsub0subiiWGQsub0subQsub0subRTsub0subi">Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i><i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik0iik0i">Convergence table for <i>k=0</i><i>k=0</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub1subQsub1subRTsub1subiiWGQsub1subQsub1subRTsub1subi">Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i><i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik1iik1i">Convergence table for <i>k=1</i><i>k=1</i></a>
      </ul>
        <li><a href="#TestresultsoniWGQsub2subQsub2subRTsub2subiiWGQsub2subQsub2subRTsub2subi">Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i><i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></a>
      <ul>
        <li><a href="#Convergencetableforik2iik2i">Convergence table for <i>k=2</i><i>k=2</i></a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-61/doc/intro.dox

 <br> 

<i>
This program was contributed by Zhuoran Wang.
Some more information about this program, as well as more numerical
results, are presented in @cite Wang2019 .
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程程序介绍了泊松方程的 "弱加勒金 "有限元方法的实现。从某种意义上说，考虑这种方法的动机与步骤51中的动机相同：我们想考虑不连续的形状函数，但又需要解决这样一个事实：与通常的连续Galerkin方法相比，所产生的问题有更多的自由度（因为，例如，每个顶点携带的自由度与相邻单元一样多）。我们还必须解决这样一个事实：与连续Galerkin方法不同，<i>every</i>一个单元上的自由度与它的每个面邻单元上的所有自由度相耦合。因此，从 "传统的 "非连续Galerkin方法得到的矩阵既大又相对密集。

step-51中的混合非连续Galerkin方法（HDG）和本教程中的弱Galerkin（WG）方法都是通过引入额外的自由度来解决耦合问题的，这些自由度的形状函数只存在于单元间的一个面上（即网格的 "骨架 "上），因此它们将相邻单元上的自由度相互 "隔离"：单元自由度只与同一单元上的其他单元自由度以及面自由度耦合，而与相邻单元上的单元自由度不耦合。因此，这些细胞自由度的形状函数的耦合确实正好耦合在一个细胞和定义在其面上的自由度上。

对于一个给定的方程，例如二阶泊松方程，HDG和WG方法的区别在于如何精确地制定连接所有这些不同形状函数的问题。事实上，对于某些WG和HDG的表述，有可能表明它们是等价的）。HDG的做法是用一阶方程系统重新表述二阶问题，然后在概念上把面的自由度看作是这个一阶系统的 "通量"。相比之下，WG方法保持二阶形式，并将面的自由度视为与主解变量相同的类型，只是限制在低维的面。为了方程的目的，在定义对其应用微分算子的含义时，人们需要以某种方式将这些形状函数 "扩展 "到单元的内部。与HDG相比，该方法的优势在于它不会因为将方程重写为一阶系统而导致未知数的增加，但它也不太容易实现。然而，正如我们在下文中所看到的，这种额外的努力并不可怕。




<a name="WeakGalerkinfiniteelementmethods"></a><h3> Weak Galerkin finite element methods </h3>


弱加勒金有限元方法（WGFEMs）使用离散的弱函数来近似标量未知数，使用离散的弱梯度来近似经典梯度。该方法最初是由王俊平和叶秀在论文<a href="https://doi.org/10.1016/j.cam.2012.10.003">
<i>A weak Galerkin finite element method for second order elliptic problems</i><i>A weak Galerkin finite element method for second order elliptic problems</i>,
J. Comput. Appl. Math., 103-115, 2013</a>中提出。与连续Galerkin方法相比，弱Galerkin方法满足重要的物理特性，即局部质量守恒和体法通量连续。它的结果是一个SPD线性系统，并且通过网格细化可以获得最佳收敛率。




<a name="Theequationtosolve"></a><h3> The equation to solve </h3> 该程序使用弱加尔金有限元法求解泊松方程。


@f{align*}{
  \nabla \cdot \left( -\mathbf{K} \nabla p \right)
    &= f,
    \qquad \mathbf{x} \in \Omega, \\
  p &=  p_D,\qquad \mathbf{x} \in \Gamma^D, \\
  \mathbf{u} \cdot \mathbf{n} &= u_N,
  \qquad \mathbf{x} \in \Gamma^N,


@f}

其中 $\Omega \subset \mathbb{R}^n (n=2,3)$ 是一个有界域。在流体流经多孔介质的背景下， $p$ 是压力， $\mathbf{K}$ 是渗透性张量， $f$ 是源项， $p_D, u_N$ 代表Dirichlet和Neumann边界条件。我们可以引入一个通量， $\mathbf{u} = -\mathbf{K} \nabla p$ ，对应于达西速度（以我们在步骤20中的方式），这个变量在下面的考虑中很重要。

在这个程序中，我们将考虑一个测试案例，即在单位平方域上的确切压力为 $p = \sin \left( \pi x\right)\sin\left(\pi y \right)$ ，具有同质Dirichelet边界条件和 $\mathbf{K}$ 身份矩阵。然后我们将计算压力、速度和通量的 $L_2$ 误差。




<a name="WeakGalerkinscheme"></a><h3> Weak Galerkin scheme </h3>


上面的泊松方程有一个解 $p$ ，需要满足问题的弱表述。

@f{equation*}
\mathcal{A}\left(p,q \right) = \mathcal{F} \left(q \right),


@f}

为所有测试函数  $q$  ，其中

@f{equation*}
\mathcal{A}\left(p,q\right)
  \dealcoloneq \int_\Omega \left(\mathbf{K} \nabla p\right) \cdot \nabla q \;\mathrm{d}x,


@f}

和

@f{equation*}
\mathcal{F}\left(q\right)
  \dealcoloneq \int_\Omega f \, q \;\mathrm{d}x


  - \int_{\Gamma^N} u_N q \; \mathrm{d}x.


@f}

在这里，我们以双线性形式进行了部分积分，我们在内部评估 $p,p$ 的梯度，在域的边界评估 $q$ 的值。所有这些都是很好的定义，因为我们假设解是在 $H^1$ 中，对它来说，取梯度和评估边界值是有效的操作。

弱Galerkin方法的想法是用一个<i>discontinuous function</i> $p_h$ 来近似精确的 $p$ 解。这个函数可能只在单元格之间的界面上不连续，由于我们也想沿着界面评估这个函数，我们不仅要规定它在单元格内部应该有什么值，还要规定它在界面上的值。我们通过说 $p_h$ 实际上是一个元组， $p_h=(p^\circ,p^\partial)$ ，尽管它实际上只是一个单一的函数，它要么等于 $p^\circ(x)$ ，要么等于 $p^\partial(x)$ ，这取决于它是在位于细胞内部还是在细胞界面的某一点 $x$ 上被评估。

然后我们想把这个近似值简单地贴到上面的双线性表格中。这适用于我们必须在边界上评估测试函数 $q_h$ 的情况（我们只需取其界面部分 $q_h^\partial$ ），但我们必须小心处理梯度，因为它只在单元格内部定义。因此，泊松方程的弱Galerkin方案被定义为

@f{equation*}
\mathcal{A}_h\left(p_h,q \right) = \mathcal{F} \left(q_h \right),


@f}

对于所有离散测试函数  $q_h$  ，其中

@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,


@f}

和

@f{equation*}
\mathcal{F}\left(q_h\right)
  \dealcoloneq \sum_{K \in \mathbb{T}} \int_K f \, q_h^\circ \;\mathrm{d}x


  - \sum_{\gamma \in \Gamma_h^N} \int_\gamma u_N q_h^\partial \;\mathrm{d}x,


@f}

关键的一点是，在这里，我们用<i>discrete weak gradient</i>算子 $\nabla_{w,d} p_h$ 代替了梯度 $\nabla p_h$ ，这对于我们特殊定义的近似 $p_h$ 是有意义的。

那么问题是该算子如何工作。为此，让我们首先说说我们是如何看待压力的离散近似值 $p_h$ 的。如上所述，"函数" $p_h$ 实际上由两部分组成：单元内部的值 $p_h^\circ$ 和界面上的 $p_h^\partial$ 。我们必须为这两部分定义离散的（有限维）函数空间；在这个程序中，我们将用FE_DGQ来表示 $p_h^\circ$ 作为细胞内部的空间（在每个细胞上定义，但一般沿界面是不连续的），用FE_FaceQ表示 $p_h^\partial$ 作为界面上的空间。

那么让我们只考虑一个单元（因为上面的积分都是逐个单元定义的，而且弱离散梯度是逐个单元定义的）。 $p_h$ 对 $K$ , $p_h|_K$ 的限制由一对 $(p_h^\circ|_K,p_h^\partial|_{\partial K})$ 组成。从本质上讲，我们可以认为 $\nabla_{w,d} p_h$ 是定义在 $K$ 上的某个函数，它近似于梯度；特别是，如果 $p_h|_K$ 是一个可微函数的限制（对 $K$ 的内部和边界--这将使它在内部和边界之间连续），那么 $\nabla_{w,d} p_h$  将只是精确梯度 $\nabla p_h$  。但是，由于 $p_h|_K$ 在 $K$ 的内部和边界之间不连续，我们需要一个更一般的定义；此外，我们不能处理任意函数，因此要求 $\nabla_{w,d} p_h$ 也在一个有限元空间中（由于梯度是一个矢量，必须是矢量值，而且由于弱梯度是在每个单元上单独定义的，因此在单元之间也将是不连续的）。

这样做的方法是以下列方式定义这个弱梯度算子 $\nabla_{w,d}|_K :
DGQ_k(K) \times DGQ_r(\partial K) \rightarrow RT_s(K)$ （其中 $RT_s(K)$ 是单元格 $K$ 上阶为 $s$ 的矢量值Raviart-Thomas空间）。

@f{equation*}{
  \int_K \mathbf v_h \cdot (\nabla_{w,d} p_h)
  =


  -\int_K (\nabla \cdot \mathbf v_h) p_h^\circ
  +\int_{\partial K} (\mathbf v_h \cdot \mathbf n) p_h^\partial,


@f}

为所有测试函数  $\mathbf v_h \in RT_s(K)$  。从本质上讲，这只是一个逐部积分公式的应用。换句话说，对于一个给定的 $p_h=(p^\circ_h,p^\partial_h)$ ，我们需要把 $\nabla_{w,d} p_h|_K$ 看作是度数为 $s$ 的Raviart-Thomas函数，对于这个函数，左手边和右手边在所有测试函数中是相等的。

那么，需要说明的一个关键点是以下几点。通常的梯度 $\nabla$ 是一个*本地*算子，它仅仅根据一个函数在某一点及其（无限小）邻域的值来计算导数，而弱离散梯度 $\nabla_{w,d}$ 却没有这个特性。它取决于它所应用的函数在整个单元上的值，包括单元的边界。然而，两者都是线性算子，从上面 $\nabla_{w,d}$ 的定义可以看出，这将允许我们在下面的讨论中通过矩阵来表示 $\nabla_{w,d}$ 。

 @note  值得指出的是，虽然弱的离散梯度是Raviart-Thomas空间 $RT_s(K)$ 在每个单元 $K$ 的一个元素，但它在单元之间是不连续的。另一方面，定义在整个网格上并由FE_RaviartThomas类实现的Raviart-Thomas空间 $RT_s=RT_s({\mathbb T})$ 代表在单元间界面上具有连续法线分量的函数。这意味着<i>globally</i>,  $\nabla_{w,d} p_h$ 不在 $RT_s$ 中，尽管它在 $K$ 中的每个单元上。   相反，它是在一个 "破碎的 "拉维-托马斯空间中，下面我们将用符号 $DGRT_s$ 来表示。 这里的术语 "破碎 "指的是 "把东西打碎 "的过程，而不是表达 "没有功能 "的同义词。因此，人们可能会（理所当然地）争辩说，在弱加尔金文献中使用的符号有点误导，但这往往取决于使用某种符号的背景--在目前的背景下，对Raviart-Thomas空间或元素的提及总是被理解为对 "破碎 "空间的提及。

 @note  deal.II恰好有一个实现了这个破碎的Raviart-Thomas空间。FE_DGRT类。因此，在本教程中，我们将简单地一直使用FE_DGRT类，尽管在所有那些我们必须计算单元格本地矩阵和向量的地方，它没有任何区别。




<a name="Representingtheweakgradient"></a><h3> Representing the weak gradient </h3>


由于 $p_h$ 是有限元空间的一个元素，我们可以像往常一样在一个基础上展开它，也就是说，我们可以写出

@f{equation*}{
  p_h(\mathbf x) = \sum_j P_j \varphi_j(\mathbf x).


@f}

这里，由于 $p_h$ 有两个分量（内部分量和界面分量），对于基函数 $\varphi_j(\mathbf x)$ 也必须如此，我们可以写成 $\varphi_j = (\varphi_j^\circ,\varphi_j^\partial)$  。如果你按照步骤8、步骤20和 @ref vector_valued "向量值问题文件模块 "中的描述，就不会感到奇怪，对于 $j$ 的某些值， $\varphi_j^\circ$ 将为零，而对于 $j$ 的其他值， $\varphi_j^\partial$ 将为零--也就是说，形状函数将是一种或另一种类型。然而，这在这里并不重要。重要的是，我们需要思考如何表示 $\nabla_{w,d} \varphi_j$ ，因为当我们想实现双线性形式时，这显然是问题中会出现的东西

@f{equation*}
\mathcal{A}_h\left(p_h,q_h\right)
  = \sum_{K \in \mathbb{T}}
    \int_K \mathbf{K} \nabla_{w,d} p_h \cdot \nabla_{w,d} q_h \;\mathrm{d}x,


@f}



关键的一点是，已知 $\nabla_{w,d} \varphi_j$ 是 "破碎的 "Raviart-Thomas空间 $DGRT_s$ 的一个成员。这意味着我们可以（在每个单元 $K$ 上分别表示

@f{equation*}
\nabla_{w,d} \varphi_j|_K
  = \sum_k C_{jk}^K \mathbf v_k|_K


@f}

其中，函数 $\mathbf v_k \in DGRT_s$ ，以及 $C^K$ 是一个维数的矩阵

@f{align*}{
 \text{dim}\left(DGQ_k(K) \times DGQ_r(K)\right) &\times \text{dim}\left(RT_s(K)\right)
  \\
 &=
 \left(\text{dim}(DGQ_k(K)) + \text{dim}(DGQ_r(K))\right) \times \text{dim}\left(RT_s(K)\right).


@f}

弱离散梯度可以被表示为一个矩阵，这不应该是一个惊喜：它是一个从一个有限维空间到另一个有限维空间的线性算子。如果为这两个空间都选择基数，那么<i>every linear operator</i>当然可以写成一个矩阵，将与算子的域空间的基数有关的扩展系数向量映射到与图像空间的基数有关的扩展系数向量）。)

利用这个扩展，我们可以很容易地使用上面的弱离散梯度的定义来定义矩阵要做什么。

@f{equation*}{
  \int_K \mathbf v_i \cdot \left(\sum_k C_{jk}^K \mathbf v_k\right)
  =


  -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
  +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,


@f}

对于所有的测试功能  $\mathbf v_i \in DGRT_s$  。

这显然导致了一个线性系统，其形式为

@f{equation*}{
  \sum_k M_{ik}^K C_{jk}^K
  =
  G_{ij}^K


@f}

与

@f{equation*}{
  M_{ik}^K = \int_K \mathbf v_i \cdot \mathbf v_k,
  \qquad\qquad
  G_{ij}^K = -\int_K (\nabla \cdot \mathbf v_i) \varphi_j^\circ
             +\int_{\partial K} (\mathbf v_i \cdot \mathbf n) \varphi_j^\partial,


@f}

因此

@f{equation*}{
  \left(C^K\right)^T = \left(M^K\right)^{-1} G^K.


@f}

(在这最后一步中，我们假设指数 $i,j,k$ 只涉及在单元 $K$ 上活动的自由度，从而确保空间 $RT_s(K)$ 上的质量矩阵是可逆的。)等价地，利用矩阵 $M$ 的对称性，我们可以看到

@f{equation*}{
  C^K = \left(G^K\right)^{T} \left(M^K\right)^{-1}.


@f}

另外值得指出的是，矩阵 $C^K$ 和 $G^K$ 当然不是正方形而是长方形。




<a name="Assemblingthelinearsystem"></a><h3> Assembling the linear system </h3>


在解释了弱离散梯度是如何定义的之后，我们现在可以回到有关方程的线性系统应该如何组装的问题上。具体来说，利用上面显示的双线性形式 ${\cal A}_h$ 的定义，我们就需要计算局部对全局矩阵的贡献元素。

@f{equation*}{
  A^K_{ij} = \int_K \left({\mathbf K} \nabla_{w,d} \varphi_i\right) \cdot \nabla_{w,d} \varphi_j.


@f}

如上所述，我们可以用Raviart-Thomas基础在每个单元格上展开 $\nabla_{w,d} \varphi_i$ ，同样，对于 $\nabla_{w,d} \varphi_j$ 也是如此。

@f{equation*}{
  A^K_{ij} = \int_K
    \left(
      {\mathbf K}
      \sum_k C_{ik}^K \mathbf v_k|_K
    \right)
    \cdot
    \sum_l C_{jl}^K \mathbf v_l|_K.


@f}

通过重新排列和，可以得到以下表达式。

@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
     \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K.


@f}

因此，如果我们有每个单元格 $K$ 的矩阵 $C^K$ ，那么我们可以很容易地计算出单元格 $K$ 对矩阵 $A$ 的贡献 $A^K$ ，如下所示。

@f{equation*}{
  A^K_{ij} =
    \sum_k \sum_l C_{ik}^K C_{jl}^K
    H^K_{kl}
    =
    \sum_k \sum_l C_{ik}^K H^K_{kl} C_{jl}^K
    =
    \left(C^K H^K (C^K)^T \right)_{ij}.


@f}

在这里。

@f{equation*}{
  H^K_{kl} =
  \int_K
    \left(
      {\mathbf K}
      \mathbf v_k|_K
    \right)
    \cdot
    \mathbf v_l|_K,


@f}

这实际上只是单元 $K$ 上的质量矩阵，使用Raviart-Thomas基础并通过渗透性张量 $\mathbf K$ 加权。这里的推导表明，弱加尔金法实际上只需要我们计算每个单元 $C^K$ 和 $H^K$ 的矩阵，然后再计算 $A^K = C^K H^K (C^K)^T$ ，这很容易计算出来。下面要显示的代码正是这样做的。

在计算出单元格 $A^K$ 对全局矩阵的贡献后，我们要做的就是将这些局部贡献 "分配 "到全局矩阵中。如何做到这一点，首先显示在步骤3和步骤4中。在目前的程序中，这将通过调用 AffineConstraints::distribute_local_to_global(). 来促进。

一个线性系统当然也需要一个右手边。除了我们只需要对每个形状函数 $\varphi_i^\circ$ 使用单元格内部部分外，这里没有与计算右手边有关的困难。




<a name="PostprocessingandiLsub2subiiLsub2subierrors"></a><h3> Post-processing and <i>L<sub>2</sub></i><i>L<sub>2</sub></i>-errors </h3> 。


前面几节的讨论已经给了我们一个线性系统，我们可以求解数值压力 $p_h$  。我们可以用它来计算变量 $\mathbf u = -{\mathbf K}\nabla p$ 的近似值，如果这是我们要解决的模型，它对应于介质在多孔介质中的流动速度。这种步骤--从离散问题的解中计算一个派生量--通常被称为 "后处理"。

这里，我们不使用 $p_h$ 的精确梯度，而是使用 $p_h$ 的离散弱梯度来计算每个元素上的速度。如上所述，在每个元素上，数值压力 $\nabla p$ 的梯度可以用离散弱梯度 $ \nabla_{w,d}\phi_i$ 来近似。

@f{equation*}
\nabla_{w,d} p_h
= \nabla_{w,d} \left(\sum_{i} P_i \phi_i\right)
= \sum_{i} P_i \nabla_{w,d}\phi_i.


@f}



在单元格 $K$ 上，数值速度 $ \mathbf{u}_h = -\mathbf{K} \nabla_{w,d}p_h$ 可写为

@f{align*}{
  \mathbf{u}_h
  &= -\mathbf{K} \nabla_{w,d} p_h
   = -\mathbf{K}\sum_{i} \sum_{j} P_i C^K_{ij}\mathbf{v}_j,


@f}

其中 $C^K$ 是上面的扩展矩阵， $\mathbf{v}_j$ 是 $RT$ 空间在一个单元上的基函数。

不幸的是， $\mathbf{K} \mathbf{v}_j$ 可能不在 $RT$ 空间中（当然，除非如果 $\mathbf K$ 是常数乘以身份矩阵）。因此，为了在有限元程序中表示它，我们需要把它投射回我们可以处理的有限维空间。在这里，我们将使用 $L_2$ 投影法将其投影回（破碎的） $RT$ 空间。

我们将每个单元格  $K$  上的投影定义为  $ \mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right) =
\sum_{k} d_{jk}\mathbf{v}_k$  。对于任何  $j$  ,  $\left( \mathbf{Q}_h \left( \mathbf{Kv}_j \right),\mathbf{v}_k \right)_K =
\left( \mathbf{Kv}_j,\mathbf{v}_k \right)_K.$  所以，与其说是上面的公式，不如说是  $K$  单元上的数字速度变成了

@f{equation*}
\mathbf{u}_h = \mathbf{Q}_h \left( -\mathbf{K}\nabla_{w,d}p_h \right) =


-\sum_i \sum_j P_i B^K_{ij}\mathbf{Q}_h \left( \mathbf{K}\mathbf{v}_j \right),


@f}

我们有以下系统来解决系数问题  $d_{jk}$  。

@f{equation*}
 \sum_j
  \left(\mathbf{v}_i,\mathbf{v}_j\right)
   d_{jk}
   =
    \left( \mathbf{Kv}_j,\mathbf{v}_k \right).


@f}

在下面的实现中，元素为 $
   d_{jk}
$ 的矩阵被称为 <code>cell_matrix_D</code>  ，而元素为 $
      \left( \mathbf{Kv}_j,\mathbf{v}_k \right)
$ 的矩阵被称为 <code>cell_matrix_E</code> 。

那么元素速度为

@f{equation*}
\mathbf{u}_h = -\sum_{i} \sum_{j}P_ic_{ij}\sum_{k}d_{jk}\mathbf{v}_k =
\sum_{k}- \left(\sum_{j} \sum_{i} P_ic_{ij}d_{jk} \right)\mathbf{v}_k,


@f}

其中 $-\sum_{j} \sum_{i} P_ic_{ij}d_{jk}$ 在代码中被称为`细胞速度'。

利用这个通过 "后处理 "得到的速度，我们可以通过以下公式定义压力、速度和通量的 $L_2$ 误差。

@f{align*}{
\|p-p_h^\circ\|^2
  &= \sum_{K \in \mathbb{T}} \|p-p_h^\circ\|_{L_2(K)}^2, \\
 \|\mathbf{u}-\mathbf{u}_h\|^2
  &= \sum_{K \in \mathbb{T}} \|\mathbf{u}-\mathbf{u}_h\|_{L_2(K)^2}^d,\\
\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|^2
  &= \sum_{K \in \mathbb{T}} \sum_{\gamma \subset \partial K}
    \frac{|K|}{|\gamma|} \|\mathbf{u} \cdot \mathbf{n} - \mathbf{u}_h \cdot \mathbf{n}\|_{L_2(\gamma)}^2,


@f}

其中 $| K |$ 为元素的面积， $\gamma$ 为元素的面， $\mathbf{n}$ 为每个面的单位法向量。这些规范中的最后一条衡量了网格单元之间界面上速度向量的法向分量的精度。缩放因子 $|K|/|\gamma|$ 的选择是为了随着网格大小的变化，缩放出界面集合的长度（或面积）的差异。

上面的第一个错误很容易用 VectorTools::integrate_difference. 计算出来，其他的需要多做一些工作，在下面的代码中实现。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name=""></a> 
 * @sect3{Include files}  这个程序是基于 step-7  、 step-20  和  step-51  ，所以下面的头文件大部分是熟悉的。我们需要以下文件，其中只有导入FE_DGRaviartThomas类的文件（即`deal.II/fe/fe_dg_vector.h`）是真正的新文件；FE_DGRaviartThomas实现了介绍中讨论的 "破碎 "Raviart-Thomas空间。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/tensor_function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/point.h> 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/fe/fe_raviart_thomas.h> 
 * #include <deal.II/fe/fe_dg_vector.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_face.h> 
 * #include <deal.II/fe/component_mask.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/data_out_faces.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 我们的第一步，像往常一样，是把所有与本教程程序有关的东西放到自己的命名空间中。
 * 

 * 
 * 
 * @code
 * namespace Step61 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="TheWGDarcyEquationclasstemplate"></a> 
 * <h3>The WGDarcyEquation class template</h3>
 * 

 * 
 * 这是本程序的主类。我们将使用弱加勒金（WG）方法求解内部和面上的数值压力，并计算出压力的 $L_2$ 误差。在后处理步骤中，我们还将计算速度和通量的 $L_2$  误差。
 * 

 * 
 * 该类的结构与以前的教程程序没有根本的不同，所以除了一个例外，没有必要对细节进行评论。该类有一个成员变量`fe_dgrt`，对应于介绍中提到的 "破碎 "的Raviart-Thomas空间。还有一个与之匹配的`dof_handler_dgrt`，表示从这个元素创建的有限元场的全局枚举，还有一个向量`darcy_velocity`，用于保持这个场的节点值。在求解压力后，我们将使用这三个变量来计算一个后处理的速度场，然后我们可以对其进行误差评估，并将其输出用于可视化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class WGDarcyEquation 
 *   { 
 *   public: 
 *     WGDarcyEquation(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void compute_postprocessed_velocity(); 
 *     void compute_velocity_errors(); 
 *     void compute_pressure_error(); 
 *     void output_results() const; 
 * 
 *     Triangulation<dim> triangulation; 
 * 
 *     FESystem<dim>   fe; 
 *     DoFHandler<dim> dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 * 
 *     FE_DGRaviartThomas<dim> fe_dgrt; 
 *     DoFHandler<dim>         dof_handler_dgrt; 
 *     Vector<double>          darcy_velocity; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandexactsolution"></a> 
 * <h3>Right hand side, boundary values, and exact solution</h3>
 * 

 * 
 * 接下来，我们定义系数矩阵 $\mathbf{K}$ （这里是身份矩阵），迪里希特边界条件，右手边 $f = 2\pi^2 \sin(\pi x) \sin(\pi y)$  ，以及与这些选择相对应的 $K$ 和 $f$ 的精确解，即 $p = \sin(\pi x) \sin(\pi y)$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Coefficient : public TensorFunction<2, dim> 
 *   { 
 *   public: 
 *     Coefficient() 
 *       : TensorFunction<2, dim>() 
 *     {} 
 * 
 *     virtual void value_list(const std::vector<Point<dim>> &points, 
 *                             std::vector<Tensor<2, dim>> &values) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   void Coefficient<dim>::value_list(const std::vector<Point<dim>> &points, 
 *                                     std::vector<Tensor<2, dim>> &  values) const 
 *   { 
 *     Assert(points.size() == values.size(), 
 *            ExcDimensionMismatch(points.size(), values.size())); 
 *     for (unsigned int p = 0; p < points.size(); ++p) 
 *       values[p] = unit_symmetric_tensor<dim>(); 
 *   } 
 * 
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     BoundaryValues() 
 *       : Function<dim>(2) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
 *                                     const unsigned int /*component*/) const 
 *   { 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   };
 * 
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> &p, 
 *                                    const unsigned int /*component*/) const 
 *   { 
 *     return (2 * numbers::PI * numbers::PI * std::sin(numbers::PI * p[0]) * 
 *             std::sin(numbers::PI * p[1])); 
 *   } 
 * 
 * @endcode
 * 
 * 实现精确压力解决方案的类有一个奇怪的地方，我们把它作为一个有两个分量的向量值来实现。(我们在构造函数中说它有两个分量，在这里我们调用基函数类的构造函数)。在`value()`函数中，我们不测试`component`参数的值，这意味着我们为向量值函数的两个分量返回相同的值。我们这样做是因为我们将本程序中使用的有限元描述为一个包含内部和界面压力的矢量值系统，当我们计算误差时，我们希望使用相同的压力解来测试这两个分量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ExactPressure : public Function<dim> 
 *   { 
 *   public: 
 *     ExactPressure() 
 *       : Function<dim>(2) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double ExactPressure<dim>::value(const Point<dim> &p, 
 *                                    const unsigned int /*component*/) const 
 *   { 
 *     return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
 *   } 
 * 
 *   template <int dim> 
 *   class ExactVelocity : public TensorFunction<1, dim> 
 *   { 
 *   public: 
 *     ExactVelocity() 
 *       : TensorFunction<1, dim>() 
 *     {} 
 * 
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   Tensor<1, dim> ExactVelocity<dim>::value(const Point<dim> &p) const 
 *   { 
 *     Tensor<1, dim> return_value; 
 *     return_value[0] = -numbers::PI * std::cos(numbers::PI * p[0]) * 
 *                       std::sin(numbers::PI * p[1]); 
 *     return_value[1] = -numbers::PI * std::sin(numbers::PI * p[0]) * 
 *                       std::cos(numbers::PI * p[1]); 
 *     return return_value; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationclassimplementation"></a> 
 * <h3>WGDarcyEquation class implementation</h3>
 * 
 * <a name="WGDarcyEquationWGDarcyEquation"></a> 
 * <h4>WGDarcyEquation::WGDarcyEquation</h4>
 * 

 * 
 * 在这个构造函数中，我们创建了一个矢量值函数的有限元空间，这里将包括用于内部和界面压力的函数， $p^\circ$  和  $p^\partial$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   WGDarcyEquation<dim>::WGDarcyEquation(const unsigned int degree) 
 *     : fe(FE_DGQ<dim>(degree), 1, FE_FaceQ<dim>(degree), 1) 
 *     , dof_handler(triangulation) 
 *     , fe_dgrt(degree) 
 *     , dof_handler_dgrt(triangulation) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationmake_grid"></a> 
 * <h4>WGDarcyEquation::make_grid</h4>
 * 

 * 
 * 我们在单位平方域上生成一个网格并对其进行细化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::make_grid() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, 0, 1); 
 *     triangulation.refine_global(5); 
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "   Total number of cells: " << triangulation.n_cells() 
 *               << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationsetup_system"></a> 
 * <h4>WGDarcyEquation::setup_system</h4>
 * 

 * 
 * 在我们创建了上面的网格后，我们分配自由度并调整矩阵和向量的大小。这个函数中唯一值得关注的部分是我们如何插值压力的边界值。由于压力由内部和界面分量组成，我们需要确保我们只插值到矢量值解空间中与界面压力相对应的分量上（因为这些分量是唯一定义在域的边界上的）。我们通过一个只针对界面压力的分量屏蔽对象来做到这一点。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     dof_handler_dgrt.distribute_dofs(fe_dgrt); 
 * 
 *     std::cout << "   Number of pressure degrees of freedom: " 
 *               << dof_handler.n_dofs() << std::endl; 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     { 
 *       constraints.clear(); 
 *       const FEValuesExtractors::Scalar interface_pressure(1); 
 *       const ComponentMask              interface_pressure_mask = 
 *         fe.component_mask(interface_pressure); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                interface_pressure_mask); 
 *       constraints.close(); 
 *     } 
 * 
 * @endcode
 * 
 * 在双线性形式中，在两个相邻单元之间的面上没有积分项，所以我们可以直接使用 <code>DoFTools::make_sparsity_pattern</code> 来计算稀疏矩阵。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationassemble_system"></a> 
 * <h4>WGDarcyEquation::assemble_system</h4>
 * 

 * 
 * 这个函数比较有趣。正如介绍中所详述的，线性系统的装配要求我们评估形状函数的弱梯度，这是Raviart-Thomas空间的一个元素。因此，我们需要定义一个Raviart-Thomas有限元对象，并有FEValues对象在正交点评估它。然后我们需要计算每个单元 $K$ 上的矩阵 $C^K$ ，为此我们需要介绍中提到的矩阵 $M^K$ 和 $G^K$ 。
 * 

 * 
 * 有一点可能不是很明显，在之前所有的教程程序中，我们总是用DoFHandler的单元格迭代器来调用 FEValues::reinit() 。这样就可以调用诸如 FEValuesBase::get_function_values() 这样的函数，在单元格的正交点上提取有限元函数的值（用DoF值的矢量表示）。为了使这种操作发挥作用，人们需要知道哪些向量元素对应于给定单元上的自由度--也就是说，正是DoFHandler类所提供的那种信息和操作。
 * 

 * 
 * 我们可以为 "破碎的 "Raviart-Thomas空间创建一个DoFHandler对象（使用FE_DGRT类），但是我们在这里真的不想这样做。至少在当前函数中，我们不需要任何与这个破碎空间相关的全局定义的自由度，而只需要引用当前单元上的这种空间的形状函数。因此，我们利用这样一个事实，即人们也可以用单元格迭代器来调用 FEValues::reinit() 的Triangulation对象（而不是DoFHandler对象）。在这种情况下，FEValues当然只能为我们提供只引用单元格的信息，而不是这些单元格上列举的自由度。所以我们不能使用 FEValuesBase::get_function_values(), ，但我们可以使用 FEValues::shape_value() 来获取当前单元上正交点的形状函数值。下面我们要利用的就是这种功能。下面给我们提供Raviart-Thomas函数信息的变量是`fe_values_rt`（和相应的`fe_face_values_rt`）对象。
 * 

 * 
 * 鉴于上述介绍，下面的声明应该是非常明显的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::assemble_system() 
 *   { 
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 
 * 
 *     FEValues<dim>     fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values); 
 *     FEFaceValues<dim> fe_face_values(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_normal_vectors | 
 *                                        update_quadrature_points | 
 *                                        update_JxW_values); 
 * 
 *     FEValues<dim>     fe_values_dgrt(fe_dgrt, 
 *                                  quadrature_formula, 
 *                                  update_values | update_gradients | 
 *                                    update_quadrature_points | 
 *                                    update_JxW_values); 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
 *                                           face_quadrature_formula, 
 *                                           update_values | 
 *                                             update_normal_vectors | 
 *                                             update_quadrature_points | 
 *                                             update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 
 * 
 *     RightHandSide<dim>  right_hand_side; 
 *     std::vector<double> right_hand_side_values(n_q_points); 
 * 
 *     const Coefficient<dim>      coefficient; 
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来，让我们声明介绍中讨论的各种单元格矩阵。
 * 

 * 
 * 
 * @code
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 *     Vector<double>     cell_solution(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 我们需要  <code>FEValuesExtractors</code>  来访问形状函数的  @p interior  和  @p face  部分。
 * 

 * 
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure_interior(0); 
 *     const FEValuesExtractors::Scalar pressure_face(1); 
 * 
 * @endcode
 * 
 * 这最终让我们在所有单元格上进行循环。在每个单元中，我们将首先计算用于构建局部矩阵的各种单元矩阵--因为它们取决于相关的单元，所以它们需要在每个单元中重新计算。我们还需要Raviart-Thomas空间的形状函数，为此我们需要首先创建一个通往三角化单元的迭代器，我们可以通过从指向DoFHandler的单元中的赋值来获得。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 * 
 *         const typename Triangulation<dim>::active_cell_iterator cell_dgrt = 
 *           cell; 
 *         fe_values_dgrt.reinit(cell_dgrt); 
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(), 
 *                                    right_hand_side_values); 
 *         coefficient.value_list(fe_values.get_quadrature_points(), 
 *                                coefficient_values); 
 * 
 * @endcode
 * 
 * 我们要计算的第一个单元矩阵是拉维-托马斯空间的质量矩阵。 因此，我们需要循环计算速度FEValues对象的所有正交点。
 * 

 * 
 * 
 * @code
 *         cell_matrix_M = 0; 
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *             { 
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
 *                 { 
 *                   const Tensor<1, dim> v_k = 
 *                     fe_values_dgrt[velocities].value(k, q); 
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 接下来我们通过使用 FullMatrix::gauss_jordan(). 对这个矩阵进行求逆 它将被用来计算后面的系数矩阵 $C^K$ 。值得一提的是，后面的 "cell_matrix_M "实际上包含了*的逆*。
 * 在这个调用之后的 $M^K$ 的*逆。
 * 

 * 
 * 
 * @code
 *         cell_matrix_M.gauss_jordan(); 
 * 
 * @endcode
 * 
 * 从介绍中，我们知道定义 $C^K$ 的方程的右边 $G^K$ 是面积分和单元积分的区别。在这里，我们对内部的贡献的负值进行了近似。这个矩阵的每个分量都是多项式空间的一个基函数与拉维-托马斯空间的一个基函数的发散之间的乘积的积分。这些基函数是在内部定义的。
 * 

 * 
 * 
 * @code
 *         cell_matrix_G = 0; 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *             { 
 *               const double div_v_i = 
 *                 fe_values_dgrt[velocities].divergence(i, q); 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const double phi_j_interior = 
 *                     fe_values[pressure_interior].value(j, q); 
 * 
 *                   cell_matrix_G(i, j) -= 
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q)); 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 接下来，我们用正交法对面的积分进行近似。每个分量都是多项式空间的基函数与Raviart-Thomas空间的基函数与法向量的点积的积分。所以我们在元素的所有面上循环，得到法向量。
 * 

 * 
 * 
 * @code
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             fe_face_values.reinit(cell, face); 
 *             fe_face_values_dgrt.reinit(cell_dgrt, face); 
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *               { 
 *                 const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *                   { 
 *                     const Tensor<1, dim> v_i = 
 *                       fe_face_values_dgrt[velocities].value(i, q); 
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                       { 
 *                         const double phi_j_face = 
 *                           fe_face_values[pressure_face].value(j, q); 
 * 
 *                         cell_matrix_G(i, j) += 
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
 *                       } 
 *                   } 
 *               } 
 *           } 
 * @endcode
 * 
 * @p cell_matrix_C 是 $G^K$ 的转置与质量矩阵的逆之间的矩阵乘积（该逆存储在 @p cell_matrix_M): 中）。
 * 
 * @code
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 
 * 
 * @endcode
 * 
 * 最后我们可以计算出本地矩阵  $A^K$  。 元素  $A^K_{ij}$  由  $\int_{E} \sum_{k,l} C_{ik} C_{jl} (\mathbf{K} \mathbf{v}_k) \cdot \mathbf{v}_l \mathrm{d}x$  得到。我们在上一步已经计算了系数 $C$ ，因此在适当地重新排列循环后得到以下结果。
 * 

 * 
 * 
 * @code
 *         local_matrix = 0; 
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
 *               { 
 *                 const Tensor<1, dim> v_k = 
 *                   fe_values_dgrt[velocities].value(k, q); 
 *                 for (unsigned int l = 0; l < dofs_per_cell_dgrt; ++l) 
 *                   { 
 *                     const Tensor<1, dim> v_l = 
 *                       fe_values_dgrt[velocities].value(l, q); 
 * 
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                       for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                         local_matrix(i, j) += 
 *                           (coefficient_values[q] * cell_matrix_C[i][k] * v_k) * 
 *                           cell_matrix_C[j][l] * v_l * fe_values_dgrt.JxW(q); 
 *                   } 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 接下来，我们计算右手边， $\int_{K} f q \mathrm{d}x$  。
 * 

 * 
 * 
 * @code
 *         cell_rhs = 0; 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               cell_rhs(i) += (fe_values[pressure_interior].value(i, q) * 
 *                               right_hand_side_values[q] * fe_values.JxW(q)); 
 *             } 
 * 
 * @endcode
 * 
 * 最后一步是将本地矩阵的组件分配到系统矩阵中，并将单元格右侧的组件转移到系统右侧。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global( 
 *           local_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimsolve"></a> 
 * <h4>WGDarcyEquation<dim>::solve</h4>
 * 

 * 
 * 这一步相当琐碎，与之前的许多教程程序相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::solve() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> solver(solver_control); 
 *     solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 
 *     constraints.distribute(solution); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_postprocessed_velocity"></a> 
 * <h4>WGDarcyEquation<dim>::compute_postprocessed_velocity</h4>
 * 

 * 
 * 在这个函数中，根据之前计算的压力解计算出速度场。速度被定义为 $\mathbf{u}_h = \mathbf{Q}_h \left(-\mathbf{K}\nabla_{w,d}p_h \right)$ ，这需要我们计算许多与系统矩阵组装相同的项。还有一些矩阵 $E^K,D^K$ 我们也需要组装（见介绍），但它们实际上只是遵循相同的模式。
 * 

 * 
 * 在这里计算与我们在`assemble_system()`函数中已经完成的相同的矩阵，当然是浪费CPU时间的。同样地，我们把那里的一些代码复制到这个函数中，这通常也是一个糟糕的主意。一个更好的实现可能会提供一个函数来封装这些重复的代码。我们也可以考虑使用计算效率和内存效率之间的经典权衡，在装配过程中每个单元只计算一次 $C^K$ 矩阵，把它们存储在边上的某个地方，然后在这里重新使用它们。例如， step-51 就是这样做的，`assemble_system()`函数需要一个参数来决定是否重新计算本地矩阵，类似的方法--也许是将本地矩阵存储在其他地方--可以适用于当前的程序）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::compute_postprocessed_velocity() 
 *   { 
 *     darcy_velocity.reinit(dof_handler_dgrt.n_dofs()); 
 * 
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_face_values(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_normal_vectors | 
 *                                        update_quadrature_points | 
 *                                        update_JxW_values); 
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt, 
 *                                  quadrature_formula, 
 *                                  update_values | update_gradients | 
 *                                    update_quadrature_points | 
 *                                    update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
 *                                           face_quadrature_formula, 
 *                                           update_values | 
 *                                             update_normal_vectors | 
 *                                             update_quadrature_points | 
 *                                             update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell      = fe.n_dofs_per_cell(); 
 *     const unsigned int dofs_per_cell_dgrt = fe_dgrt.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points      = fe_values.get_quadrature().size(); 
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 
 * 
 *     const unsigned int n_face_q_points = fe_face_values.get_quadrature().size(); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices_dgrt( 
 *       dofs_per_cell_dgrt); 
 * 
 *     FullMatrix<double> cell_matrix_M(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
 *     FullMatrix<double> cell_matrix_G(dofs_per_cell_dgrt, dofs_per_cell); 
 *     FullMatrix<double> cell_matrix_C(dofs_per_cell, dofs_per_cell_dgrt); 
 *     FullMatrix<double> cell_matrix_D(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
 *     FullMatrix<double> cell_matrix_E(dofs_per_cell_dgrt, dofs_per_cell_dgrt); 
 * 
 *     Vector<double> cell_solution(dofs_per_cell); 
 *     Vector<double> cell_velocity(dofs_per_cell_dgrt); 
 * 
 *     const Coefficient<dim>      coefficient; 
 *     std::vector<Tensor<2, dim>> coefficient_values(n_q_points_dgrt); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure_interior(0); 
 *     const FEValuesExtractors::Scalar pressure_face(1); 
 * 
 * @endcode
 * 
 * 在介绍中，我们解释了如何计算单元上的数值速度。我们需要每个单元上的压力解值、格拉姆矩阵的系数和 $L_2$ 投影的系数。我们已经计算了全局解，所以我们将从全局解中提取单元解。格拉姆矩阵的系数在我们计算压力的系统矩阵时已经计算过了。我们在这里也要这样做。对于投影的系数，我们做矩阵乘法，即用格拉姆矩阵的倒数乘以 $(\mathbf{K} \mathbf{w}, \mathbf{w})$ 的矩阵作为组成部分。然后，我们将所有这些系数相乘，称之为β。数值速度是贝塔和拉维尔特-托马斯空间的基础函数的乘积。
 * 

 * 
 * 
 * @code
 *     typename DoFHandler<dim>::active_cell_iterator 
 *       cell = dof_handler.begin_active(), 
 *       endc = dof_handler.end(), cell_dgrt = dof_handler_dgrt.begin_active(); 
 *     for (; cell != endc; ++cell, ++cell_dgrt) 
 *       { 
 *         fe_values.reinit(cell); 
 *         fe_values_dgrt.reinit(cell_dgrt); 
 * 
 *         coefficient.value_list(fe_values_dgrt.get_quadrature_points(), 
 *                                coefficient_values); 
 * 
 * @endcode
 * 
 * 这个 <code>cell_matrix_E</code> 的分量是 $(\mathbf{K} \mathbf{w}, \mathbf{w})$ 的积分。  <code>cell_matrix_M</code> 是格拉姆矩阵。
 * 

 * 
 * 
 * @code
 *         cell_matrix_M = 0; 
 *         cell_matrix_E = 0; 
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *             { 
 *               const Tensor<1, dim> v_i = fe_values_dgrt[velocities].value(i, q); 
 *               for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
 *                 { 
 *                   const Tensor<1, dim> v_k = 
 *                     fe_values_dgrt[velocities].value(k, q); 
 * 
 *                   cell_matrix_E(i, k) += 
 *                     (coefficient_values[q] * v_i * v_k * fe_values_dgrt.JxW(q)); 
 * 
 *                   cell_matrix_M(i, k) += (v_i * v_k * fe_values_dgrt.JxW(q)); 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 为了计算介绍中提到的矩阵 $D$ ，我们就需要按照介绍中的解释来评估 $D=M^{-1}E$ 。
 * 

 * 
 * 
 * @code
 *         cell_matrix_M.gauss_jordan(); 
 *         cell_matrix_M.mmult(cell_matrix_D, cell_matrix_E); 
 * 
 * @endcode
 * 
 * 然后，我们还需要再次计算矩阵 $C$ ，用于评估弱离散梯度。这与组装系统矩阵时使用的代码完全相同，所以我们只需从那里复制它。
 * 

 * 
 * 
 * @code
 *         cell_matrix_G = 0; 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *             { 
 *               const double div_v_i = 
 *                 fe_values_dgrt[velocities].divergence(i, q); 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const double phi_j_interior = 
 *                     fe_values[pressure_interior].value(j, q); 
 * 
 *                   cell_matrix_G(i, j) -= 
 *                     (div_v_i * phi_j_interior * fe_values.JxW(q)); 
 *                 } 
 *             } 
 * 
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             fe_face_values.reinit(cell, face); 
 *             fe_face_values_dgrt.reinit(cell_dgrt, face); 
 * 
 *             for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *               { 
 *                 const Tensor<1, dim> &normal = fe_face_values.normal_vector(q); 
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell_dgrt; ++i) 
 *                   { 
 *                     const Tensor<1, dim> v_i = 
 *                       fe_face_values_dgrt[velocities].value(i, q); 
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                       { 
 *                         const double phi_j_face = 
 *                           fe_face_values[pressure_face].value(j, q); 
 * 
 *                         cell_matrix_G(i, j) += 
 *                           ((v_i * normal) * phi_j_face * fe_face_values.JxW(q)); 
 *                       } 
 *                   } 
 *               } 
 *           } 
 *         cell_matrix_G.Tmmult(cell_matrix_C, cell_matrix_M); 
 * 
 * @endcode
 * 
 * 最后，我们需要提取对应于当前单元的压力未知数。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_values(solution, cell_solution); 
 * 
 * @endcode
 * 
 * 我们现在可以计算当地的速度未知数（相对于我们将 $-\mathbf K \nabla_{w,d} p_h$ 项投影到的Raviart-Thomas空间而言）。
 * 

 * 
 * 
 * @code
 *         cell_velocity = 0; 
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               cell_velocity(k) += 
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 
 * 
 * @endcode
 * 
 * 我们计算达西速度。这与cell_velocity相同，但用于绘制Darcy速度图。
 * 

 * 
 * 
 * @code
 *         cell_dgrt->get_dof_indices(local_dof_indices_dgrt); 
 *         for (unsigned int k = 0; k < dofs_per_cell_dgrt; ++k) 
 *           for (unsigned int j = 0; j < dofs_per_cell_dgrt; ++j) 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               darcy_velocity(local_dof_indices_dgrt[k]) += 
 *                 -(cell_solution(i) * cell_matrix_C(i, j) * cell_matrix_D(k, j)); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_pressure_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_pressure_error</h4>
 * 

 * 
 * 这一部分是为了计算压力的 $L_2$ 误差。 我们定义一个向量，用来保存每个单元上的误差规范。接下来，我们使用 VectorTool::integrate_difference() 来计算每个单元上的 $L_2$ 准则的误差。然而，我们实际上只关心解向量的内部分量的误差（我们甚至不能评估正交点的界面压力，因为这些都位于单元格的内部），因此必须使用一个权重函数，确保解变量的界面分量被忽略。这是通过使用ComponentSelectFunction来实现的，其参数表明我们要选择哪个分量（零分量，即内部压力）以及总共有多少分量（两个）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::compute_pressure_error() 
 *   { 
 *     Vector<float> difference_per_cell(triangulation.n_active_cells()); 
 *     const ComponentSelectFunction<dim> select_interior_pressure(0, 2); 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       ExactPressure<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(fe.degree + 2), 
 *                                       VectorTools::L2_norm, 
 *                                       &select_interior_pressure); 
 * 
 *     const double L2_error = difference_per_cell.l2_norm(); 
 *     std::cout << "L2_error_pressure " << L2_error << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationdimcompute_velocity_error"></a> 
 * <h4>WGDarcyEquation<dim>::compute_velocity_error</h4>
 * 

 * 
 * 在这个函数中，我们评估每个单元的速度的 $L_2$ 误差，以及面的流量的 $L_2$ 误差。该函数依赖于之前计算过的`compute_postprocessed_velocity()`函数，该函数根据之前计算过的压力解来计算速度场。
 * 

 * 
 * 我们将评估每个单元的速度，并计算数值速度和精确速度之间的差异。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::compute_velocity_errors() 
 *   { 
 *     const QGauss<dim>     quadrature_formula(fe_dgrt.degree + 1); 
 *     const QGauss<dim - 1> face_quadrature_formula(fe_dgrt.degree + 1); 
 * 
 *     FEValues<dim> fe_values_dgrt(fe_dgrt, 
 *                                  quadrature_formula, 
 *                                  update_values | update_gradients | 
 *                                    update_quadrature_points | 
 *                                    update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_face_values_dgrt(fe_dgrt, 
 *                                           face_quadrature_formula, 
 *                                           update_values | 
 *                                             update_normal_vectors | 
 *                                             update_quadrature_points | 
 *                                             update_JxW_values); 
 * 
 *     const unsigned int n_q_points_dgrt = fe_values_dgrt.get_quadrature().size(); 
 *     const unsigned int n_face_q_points_dgrt = 
 *       fe_face_values_dgrt.get_quadrature().size(); 
 * 
 *     std::vector<Tensor<1, dim>> velocity_values(n_q_points_dgrt); 
 *     std::vector<Tensor<1, dim>> velocity_face_values(n_face_q_points_dgrt); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 * 
 *     const ExactVelocity<dim> exact_velocity; 
 * 
 *     double L2_err_velocity_cell_sqr_global = 0; 
 *     double L2_err_flux_sqr                 = 0; 
 * 
 * @endcode
 * 
 * 在之前计算了后处理的速度之后，我们在这里只需要提取每个单元和面的相应数值，并与精确的数值进行比较。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell_dgrt : dof_handler_dgrt.active_cell_iterators()) 
 *       { 
 *         fe_values_dgrt.reinit(cell_dgrt); 
 * 
 * @endcode
 * 
 * 首先计算后处理的速度场与精确速度场之间的 $L_2$ 误差。
 * 

 * 
 * 
 * @code
 *         fe_values_dgrt[velocities].get_function_values(darcy_velocity, 
 *                                                        velocity_values); 
 *         double L2_err_velocity_cell_sqr_local = 0; 
 *         for (unsigned int q = 0; q < n_q_points_dgrt; ++q) 
 *           { 
 *             const Tensor<1, dim> velocity = velocity_values[q]; 
 *             const Tensor<1, dim> true_velocity = 
 *               exact_velocity.value(fe_values_dgrt.quadrature_point(q)); 
 * 
 *             L2_err_velocity_cell_sqr_local += 
 *               ((velocity - true_velocity) * (velocity - true_velocity) * 
 *                fe_values_dgrt.JxW(q)); 
 *           } 
 *         L2_err_velocity_cell_sqr_global += L2_err_velocity_cell_sqr_local; 
 * 
 * @endcode
 * 
 * 为了重建通量，我们需要单元格和面的大小。由于通量是按面计算的，我们必须在每个单元的所有四个面上进行循环。为了计算面的速度，我们从之前计算的`darcy_velocity`中提取正交点的值。然后，我们计算法线方向的速度平方误差。最后，我们通过对面和单元面积的适当缩放来计算单元上的 $L_2$ 通量误差，并将其加入全局误差。
 * 

 * 
 * 
 * @code
 *         const double cell_area = cell_dgrt->measure(); 
 *         for (const auto &face_dgrt : cell_dgrt->face_iterators()) 
 *           { 
 *             const double face_length = face_dgrt->measure(); 
 *             fe_face_values_dgrt.reinit(cell_dgrt, face_dgrt); 
 *             fe_face_values_dgrt[velocities].get_function_values( 
 *               darcy_velocity, velocity_face_values); 
 * 
 *             double L2_err_flux_face_sqr_local = 0; 
 *             for (unsigned int q = 0; q < n_face_q_points_dgrt; ++q) 
 *               { 
 *                 const Tensor<1, dim> velocity = velocity_face_values[q]; 
 *                 const Tensor<1, dim> true_velocity = 
 *                   exact_velocity.value(fe_face_values_dgrt.quadrature_point(q)); 
 * 
 *                 const Tensor<1, dim> &normal = 
 *                   fe_face_values_dgrt.normal_vector(q); 
 * 
 *  
 *                   ((velocity * normal - true_velocity * normal) * 
 *                    (velocity * normal - true_velocity * normal) * 
 *                    fe_face_values_dgrt.JxW(q)); 
 *               } 
 *             const double err_flux_each_face = 
 *               L2_err_flux_face_sqr_local / face_length * cell_area; 
 *             L2_err_flux_sqr += err_flux_each_face; 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 将所有单元和面的误差相加后，我们进行平方根计算，得到速度和流量的 $L_2$ 误差。我们将这些数据输出到屏幕上。
 * 

 * 
 * 
 * @code
 *     const double L2_err_velocity_cell = 
 *       std::sqrt(L2_err_velocity_cell_sqr_global); 
 *     const double L2_err_flux_face = std::sqrt(L2_err_flux_sqr); 
 * 
 *     std::cout << "L2_error_vel:  " << L2_err_velocity_cell << std::endl 
 *               << "L2_error_flux: " << L2_err_flux_face << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationoutput_results"></a> 
 * <h4>WGDarcyEquation::output_results</h4>
 * 

 * 
 * 我们有两组结果要输出：内部解和骨架解。我们使用 <code>DataOut</code> 来显示内部结果。骨架结果的图形输出是通过使用DataOutFaces类完成的。
 * 

 * 
 * 在这两个输出文件中，内部和面的变量都被存储。对于界面输出，输出文件只是包含了内部压力对面的插值，但是因为没有确定从两个相邻的单元中得到的是哪一个内部压力变量，所以在界面输出文件中最好是忽略内部压力。相反，对于单元格内部输出文件，当然不可能显示任何界面压力 $p^\partial$ ，因为这些压力只适用于界面，而不是单元格内部。因此，你会看到它们被显示为一个无效的值（比如一个无穷大）。
 * 

 * 
 * 对于单元内部的输出，我们还想输出速度变量。这有点棘手，因为它生活在同一个网格上，但使用不同的DoFHandler对象（压力变量生活在`dof_handler`对象上，达西速度生活在`dof_handler_dgrt`对象上）。幸运的是， DataOut::add_data_vector() 函数有一些变化，允许指定一个矢量对应的DoFHandler，因此我们可以在同一个文件中对两个DoFHandler对象的数据进行可视化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::output_results() const 
 *   { 
 *     { 
 *       DataOut<dim> data_out; 
 * 
 * @endcode
 * 
 * 首先将压力解决方案附加到DataOut对象上。
 * 

 * 
 * 
 * @code
 *       const std::vector<std::string> solution_names = {"interior_pressure", 
 *                                                        "interface_pressure"}; 
 *       data_out.add_data_vector(dof_handler, solution, solution_names); 
 * 
 * @endcode
 * 
 * 然后对达西速度场做同样的处理，并继续将所有内容写进文件。
 * 

 * 
 * 
 * @code
 *       const std::vector<std::string> velocity_names(dim, "velocity"); 
 *       const std::vector< 
 *         DataComponentInterpretation::DataComponentInterpretation> 
 *         velocity_component_interpretation( 
 *           dim, DataComponentInterpretation::component_is_part_of_vector); 
 *       data_out.add_data_vector(dof_handler_dgrt, 
 *                                darcy_velocity, 
 *                                velocity_names, 
 *                                velocity_component_interpretation); 
 * 
 *       data_out.build_patches(fe.degree); 
 *       std::ofstream output("solution_interior.vtu"); 
 *       data_out.write_vtu(output); 
 *     } 
 * 
 *     { 
 *       DataOutFaces<dim> data_out_faces(false); 
 *       data_out_faces.attach_dof_handler(dof_handler); 
 *       data_out_faces.add_data_vector(solution, "Pressure_Face"); 
 *       data_out_faces.build_patches(fe.degree); 
 *       std::ofstream face_output("solution_interface.vtu"); 
 *       data_out_faces.write_vtu(face_output); 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="WGDarcyEquationrun"></a> 
 * <h4>WGDarcyEquation::run</h4>
 * 

 * 
 * 这是主类的最后一个函数。它调用我们类的其他函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WGDarcyEquation<dim>::run() 
 *   { 
 *     std::cout << "Solving problem in " << dim << " space dimensions." 
 *               << std::endl; 
 *     make_grid(); 
 *     setup_system(); 
 *     assemble_system(); 
 *     solve(); 
 *     compute_postprocessed_velocity(); 
 *     compute_pressure_error(); 
 *     compute_velocity_errors(); 
 *     output_results(); 
 *   } 
 * 
 * } // namespace Step61 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 这是主函数。我们可以在这里改变维度以在3D中运行。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step61::WGDarcyEquation<2> wg_darcy(0); 
 *       wg_darcy.run(); 
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
examples/step-61/doc/results.dox



<a name="Results"></a><h1>Results</h1>


我们在运行程序时，右手边会产生解  $p = \sin(\pi x) \sin(\pi y)$  ，并且在域  $\Omega = (0,1)^2$  中具有同质的迪里希特边界条件。此外，我们选择微分算子 $\mathbf{K}$ 中的系数矩阵作为身份矩阵。我们使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 、 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 和 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 元素组合测试这一设置，可以通过使用`main()`中`WGDarcyEquation`对象的适当构造参数来选择。然后我们将可视化单元内部和面上的压力值。随着网格的细化，压力、速度和流量的收敛率对于 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 应该是1，对于 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 是2，对于 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 是3。




<a name="TestresultsoniWGQsub0subQsub0subRTsub0subiiWGQsub0subQsub0subRTsub0subi"></a><h3>Test results on <i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i><i>WG(Q<sub>0</sub>,Q<sub>0</sub>;RT<sub>[0]</sub>)</i></h3> 。


下面的数字显示了使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 元素的内部压力和表面压力。网格分别细化了2倍（顶部）和4倍（底部）。(这个数字可以在`make_grid()`函数中调整)。当网格较粗时，可以看到面压 $p^\partial$ 整齐地位于两个相邻单元的内压 $p^\circ$ 的数值之间。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_2.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_2.png" alt=""></td>
  </tr>
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg000_3d_4.png" alt=""></td>
  </tr>
</table> 

从图中我们可以看出，随着网格的细化，最大和最小的压力值正在接近我们的预期值。由于网格是一个矩形网格，每个方向的单元数是偶数，所以我们有对称的解决方案。从右边的三维图中，我们可以看到在 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 上，压力在单元的内部是一个常数，正如预期的那样。

<a name="Convergencetableforik0iik0i"></a><h4>Convergence table for <i>k=0</i><i>k=0</i></h4> 。


我们用不同的细化网格（在 "make_grid() "函数中选择）运行代码，得到压力、速度和通量（如引言中定义的）的以下收敛率。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
   <td>   2                  </td><td>    1.587e-01        </td><td>        5.113e-01               </td><td>   7.062e-01 </td>
  </tr>
  <tr>
   <td>   3                  </td><td>    8.000e-02        </td><td>        2.529e-01               </td><td>   3.554e-01 </td>
  </tr>
  <tr>
   <td>   4                  </td><td>    4.006e-02        </td><td>        1.260e-01               </td><td>   1.780e-01 </td>
  </tr>
  <tr>
   <td>   5                  </td><td>    2.004e-02        </td><td>        6.297e-02               </td><td>   8.902e-02 </td>
  </tr>
  <tr>
   <th>Conv.rate             </th><th>      1.00           </th><th>          1.00                  </th><th>      1.00   </th>
  </tr>
</table> 

我们可以看到， $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 的收敛率在1左右。当然，这与我们的理论预期相符。




<a name="TestresultsoniWGQsub1subQsub1subRTsub1subiiWGQsub1subQsub1subRTsub1subi"></a><h3>Test results on <i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i><i>WG(Q<sub>1</sub>,Q<sub>1</sub>;RT<sub>[1]</sub>)</i></h3> 。


我们可以用下一个更高的多项式度数重复上面的实验。下面的数字是使用 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 实现的内部压力和表面压力。网格被细化了4次。  与之前使用 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 的数字相比，在每个单元上，解决方案不再是恒定的，因为我们现在使用双线性多项式来做近似。因此，在一个内部有4个压力值，在每个面上有2个压力值。

 <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_2d_4.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg111_3d_4.png" alt=""></td>
  </tr>
</table> 

与 $\mbox{WG}(Q_0,Q_0;RT_{[0]})$ 组合的相应图像相比，现在的解决方案大大增加了准确性，特别是在界面上如此接近于连续，以至于我们不再能够区分相邻单元上的界面压力 $p^\partial$ 和内部压力 $p^\circ$ 。

<a name="Convergencetableforik1iik1i"></a><h4>Convergence table for <i>k=1</i><i>k=1</i></h4> 。


以下是我们使用 $\mbox{WG}(Q_1,Q_1;RT_{[1]})$ 元素组合得到的压力、速度和通量的收敛率。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
    <td>  2           </td><td>           1.613e-02      </td><td>          5.093e-02     </td><td>             7.167e-02   </td>
  </tr>
  <tr>
    <td>  3           </td><td>           4.056e-03      </td><td>          1.276e-02     </td><td>             1.802e-02    </td>
  </tr>
  <tr>
    <td>  4           </td><td>           1.015e-03      </td><td>          3.191e-03     </td><td>             4.512e-03  </td>
  </tr>
  <tr>
    <td>  5           </td><td>           2.540e-04      </td><td>          7.979e-04     </td><td>             1.128e-03  </td>
  </tr>
  <tr>
    <th>Conv.rate     </th><th>              2.00        </th><th>             2.00       </th><th>                 2.00    </th>
  </tr>
</table> 

 $WG(Q_1,Q_1;RT_{[1]})$ 的收敛率在2左右，符合预期。




<a name="TestresultsoniWGQsub2subQsub2subRTsub2subiiWGQsub2subQsub2subRTsub2subi"></a><h3>Test results on <i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i><i>WG(Q<sub>2</sub>,Q<sub>2</sub>;RT<sub>[2]</sub>)</i></h3> 。


让我们再提高一个多项式等级。以下是使用 $WG(Q_2,Q_2;RT_{[2]})$ 实现的内部压力和表面压力，网格大小为 $h = 1/32$ （即5个全局网格细化步骤）。在程序中，我们在生成图形输出时使用`data_out_face.build_patches(fe.degree)`（参见 DataOut::build_patches()), 的文档，这里意味着我们将每个2d单元内部分成4个子单元，以便提供更好的二次多项式的可视化。   <table align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_2d_5.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-61.wg222_3d_5.png" alt=""></td>
  </tr>
</table> 




<a name="Convergencetableforik2iik2i"></a><h4>Convergence table for <i>k=2</i><i>k=2</i></h4> 。


和以前一样，我们可以使用 $L_2$ 组合生成压力、速度和流量的 $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 误差的收敛数据。

 <table align="center" class="doxtable">
  <tr>
   <th>number of refinements </th><th>  $\|p-p_h^\circ\|$  </th><th>  $\|\mathbf{u}-\mathbf{u}_h\|$ </th><th> $\|(\mathbf{u}-\mathbf{u}_h) \cdot \mathbf{n}\|$ </th>
  </tr>
  <tr>
     <td>  2               </td><td>       1.072e-03       </td><td>         3.375e-03       </td><td>           4.762e-03   </td>
  </tr>
  <tr>
    <td>   3               </td><td>       1.347e-04       </td><td>         4.233e-04       </td><td>           5.982e-04    </td>
  </tr>
  <tr>
    <td>   4               </td><td>       1.685e-05      </td><td>          5.295e-05       </td><td>           7.487e-05  </td>
  </tr>
  <tr>
    <td>   5               </td><td>       2.107e-06      </td><td>          6.620e-06       </td><td>           9.362e-06  </td>
  </tr>
  <tr>
    <th>Conv.rate          </th><th>         3.00         </th><th>            3.00          </th><th>              3.00    </th>
  </tr>
</table> 

再一次， $\mbox{WG}(Q_2,Q_2;RT_{[2]})$ 的收敛率符合预期，其数值在3左右。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-61.cc"
*/
