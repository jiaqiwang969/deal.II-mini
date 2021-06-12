/**
@page step_46 The step-46 tutorial program
This tutorial depends on step-8, step-22, step-27.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thegeneralidea">The general idea</a>
        <li><a href="#Implementation">Implementation</a>
        <li><a href="#Specificsoftheimplementation"> Specifics of the implementation </a>
      <ul>
        <li><a href="#Dealingwiththeinterfaceterms">Dealing with the interface terms</a>
        <li><a href="#Velocityboundaryconditionsontheinterface">Velocity boundary conditions on the interface</a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
      <ul>
        <li><a href="#Identifyingwhichsubdomainacellisin">Identifying which subdomain a cell is in</a>
        <li><a href="#Linearsolvers">Linear solvers</a>
        <li><a href="#Meshrefinement">Mesh refinement</a>
    </ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeFluidStructureProblemcodeclasstemplate">The <code>FluidStructureProblem</code> class template</a>
        <li><a href="#Boundaryvaluesandrighthandside">Boundary values and right hand side</a>
        <li><a href="#ThecodeFluidStructureProblemcodeimplementation">The <code>FluidStructureProblem</code> implementation</a>
      <ul>
        <li><a href="#Constructorsandhelperfunctions">Constructors and helper functions</a>
        <li><a href="#Meshesandassigningsubdomains">Meshes and assigning subdomains</a>
        <li><a href="#codeFluidStructureProblemsetup_dofscode"><code>FluidStructureProblem::setup_dofs</code></a>
        <li><a href="#codeFluidStructureProblemassemble_systemcode"><code>FluidStructureProblem::assemble_system</code></a>
        <li><a href="#codeFluidStructureProblemsolvecode"><code>FluidStructureProblem::solve</code></a>
        <li><a href="#codeFluidStructureProblemoutput_resultscode"><code>FluidStructureProblem::output_results</code></a>
        <li><a href="#codeFluidStructureProblemrefine_meshcode"><code>FluidStructureProblem::refine_mesh</code></a>
        <li><a href="#codeFluidStructureProblemruncode"><code>FluidStructureProblem::run</code></a>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#2dresults">2d results</a>
        <li><a href="#3dresults">3d results</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
        <li><a href="#Refinementindicators">Refinement indicators</a>
        <li><a href="#Verification">Verification</a>
        <li><a href="#Bettermodels">Better models</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-46/doc/intro.dox

 <br> 

<i>This program was contributed by Wolfgang Bangerth.
<br>
This material is based upon work partly supported by the National Science
Foundation under Award No. EAR-0949446 and The University of California
&ndash; Davis. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not necessarily
reflect the views of the National Science Foundation or of The University of
California &ndash; Davis.  </i>


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序处理的是在领域的不同部分耦合不同物理学的问题。具体来说，让我们考虑以下情况，即把斯托克斯流体与弹性固体耦合起来（这两个问题以前在步骤22和步骤8中分别讨论过，你可能想在那里阅读一下各个方程式）。

- 在 $\Omega$ 的 $\Omega_f$ 部分，我们有一个流动的流体，满足与时间无关的斯托克斯方程（以涉及应变张量的形式）。   @f{align*}


    -2\eta\nabla \cdot \varepsilon(\mathbf v) + \nabla p &= 0,
          \qquad \qquad && \text{in}\ \Omega_f\\


    -\nabla \cdot \mathbf v &= 0  && \text{in}\ \Omega_f.
  @f}

  这里， $\mathbf v, p$ 分别是流体的速度和压力。   我们规定了部分外部边界上的速度，@f{align*}
    \mathbf v = \mathbf v_0 \qquad\qquad
     \text{on}\ \Gamma_{f,1} \subset \partial\Omega \cap \partial\Omega_f
  @f} 。

  而我们假设外部边界的其余部分为自由流动条件，@f{align*}
    (2\eta \varepsilon(\mathbf v) - p \mathbf 1) \cdot \mathbf n = 0
     \qquad\qquad
     \text{on}\ \Gamma_{f,2} = \partial\Omega \cap \partial\Omega_f \backslash
     \Gamma_{f,1}.
  @f}



- 域的其余部分， $\Omega_s = \Omega \backslash \Omega_f$ 被一个固体占据，其变形场 $\mathbf u$ 满足弹性方程，@f{align*}


    -\nabla \cdot C \varepsilon(\mathbf u) = 0 \qquad\qquad
    & \text{in}\ \Omega_s,
  @f}。

  其中 $C$ 是等级4的弹性张量（我们将使用一个特别简单的形式，假设固体是各向异性的）。   它在对沿固体边界流动的流体所施加的力的反应中发生变形。我们假设这种变形非常小，以至于它对流体没有反馈作用，也就是说，这种耦合只是在一个方向。为了简单起见，我们将假设固体的外部边界是被夹紧的，即@f{align*}
    \mathbf u = \mathbf 0 \qquad\qquad
     \text{on}\ \Gamma_{s,1} = \partial\Omega \cap \partial\Omega_s
  @f} 。



- 作为小位移假设的结果，我们将对流体和固体之间的界面提出以下边界条件：首先，我们对流体没有滑移的边界条件，@f{align*}
    \mathbf v = \mathbf 0 \qquad\qquad
     \text{on}\ \Gamma_{i} = \partial\Omega_s \cap \partial\Omega_f.
  @f} 。

  其次，固体上的力（牵引力）等于来自流体的法向应力，@f{align*}
    (C \varepsilon(\mathbf u)) \mathbf n =
    (2 \eta \varepsilon(\mathbf v) - p \mathbf 1) \mathbf n \qquad\qquad
     \text{on}\ \Gamma_{i} = \partial\Omega_s \cap \partial\Omega_f,
  @f} 。

  其中 $\mathbf{n}$ 是 $\Gamma_{i}$ 上的法向量，从固体指向流体。

我们通过遵循我们通常的规则，即从左边乘以一个测试函数并在域上进行积分，得到这个问题的弱表述。那么它看起来像这样。找到 $y = \{\mathbf v, p,
\mathbf u\} \in Y \subset H^1(\Omega_f)^d \times L_2(\Omega_f) \times
H^1(\Omega_s)^d$ ，使得

@f{align*}
	2 \eta (\varepsilon(\mathbf a), \varepsilon(\mathbf v))_{\Omega_f}


	- (\nabla \cdot \mathbf a, p)_{\Omega_f}


	- (q, \nabla \cdot \mathbf v)_{\Omega_f} &
	\\
	+ (\varepsilon(\mathbf b), C \varepsilon(\mathbf u))_{\Omega_s} &
	\\


	- (\mathbf b,
           (2 \eta \varepsilon(\mathbf v) - p \mathbf 1) \mathbf n)_{\Gamma_i}
	&=
	0,


@f}

为所有测试函数 $\mathbf a, q, \mathbf b$ ；第一、第二和第三行分别对应于流体、固体和界面贡献。请注意， $Y$ 只是上述空间的一个子空间，以适应各种迪里切特边界条件。

当然，只要有两个Triangulation和两个DoFHandler对象，两个子域各一个，就可以实现这种耦合。另一方面，如果有一个知道整个问题离散化的单一DoFHandler对象，那么deal.II的使用就简单多了。

这个程序是关于如何实现这一点的。请注意，我们的目标并不是要提出一个特别有用的物理模型（一个现实的流固交互模型必须考虑到固体的有限变形和它对流体的影响）：这毕竟只是一个旨在演示技术的教程程序，而不是为了解决实际问题。此外，我们将假设子域之间的界面与粗略的网格单元面对齐。




<a name="Thegeneralidea"></a><h3>The general idea</h3>


在讨论更多细节之前，让我们先说明一下：这是一个有多个解变量的问题；为此，你可能想先阅读一下 @ref vector_valued 文档模块，它介绍了我们处理有多个解变量问题的基本哲学框架。但回到手头的问题上。

在deal.II中实现这类问题的基本思路如下：在问题表述中，速度和压力变量 $\mathbf v, p$ 只存在于流体子域 $\Omega_f$ 中。但我们假设将它们以零点扩展到整个域 $\Omega$ （在一般情况下，这意味着它们沿 $\Gamma_i$ 将是不连续的）。那么，什么是这些变量的适当函数空间呢？我们知道，在 $\Omega_f$ 上我们应该要求 $\mathbf v \in H^1(\Omega_f)^d, p \in L_2(\Omega_f)$ ，所以对于 $\tilde{\mathbf v}, \tilde p$ 到整个域的扩展，下面出现了一组有用的函数空间。

@f{align*}
  \tilde {\mathbf v} &\in V
   = \{\tilde {\mathbf v}|_{\Omega_f} \in H^1(\Omega_f)^d, \quad
       \tilde {\mathbf v}|_{\Omega_s} = 0 \}
  \\
  \tilde p &\in P
  = \{\tilde p|_{\Omega_f} \in L_2(\Omega_f), \quad
       \tilde p|_{\Omega_s} = 0 \}.


@f}

(由于这对目前的讨论并不重要，我们从函数空间的选择中省略了边界值的问题；这个问题也影响到我们是否可以为压力选择 $L_2$ 或者我们是否必须为压力选择空间 $L_{2,0}(\Omega_f)=\{q\in L_2(\Omega_f): \int_{\Omega_f} q
= 0\}$ 。不过，这些问题都与下面的讨论无关)。

请注意，这些确实是一个具有明显规范的线性函数空间。由于在实践中不可能发生混淆，因此我们今后将再次省略省略号，以表示一个函数对整个域的扩展，并简单地用 $\mathbf v, p$ 来指代原始函数和扩展函数。

对于离散化，我们需要 $V_h,P_h$ 的有限维子空间 $V, P$  。对于斯托克斯，我们从步骤22中知道，适当的选择是 $Q_{p+1}^d\times Q_P$ ，但这只适用于流体占据的那部分域。对于扩展场，让我们使用以下定义在三角形上的子空间  $\mathbb T$  。

@f{align*}
  V_h
   &= \{{\mathbf v}_h \quad | \quad
       \forall K \in {\mathbb T}:
       {\mathbf v}_h|_K \in Q_{p+1}^d\  \text{if}\ K\subset {\Omega_f}, \quad
       {\mathbf v}_h|_{\Omega_f}\ \text{is continuous}, \quad
       {\mathbf v}_h|_K = 0\ \text{if}\ K\subset {\Omega_s}\}
   && \subset V
  \\
  P_h
  &= \{ p_h \quad | \quad
       \forall K \in {\mathbb T}:
       p_h|_K \in Q_p\  \text{if}\ K\subset {\Omega_f}, \quad
       p_h|_{\Omega_f}\ \text{is continuous}, \quad
       p_h|_K = 0\ \text{if}\ K\subset {\Omega_s}\ \}
   && \subset P.


@f}

换句话说，在 $\Omega_f$ 上，我们选择了通常的离散空间，但我们保留了由零扩展的（不连续的）。要说明的是，我们现在需要一个描述在单元上为零的函数的有限元空间&mdash;而这正是FE_Nothing类的用武之地：它描述了一个常数为零的函数的有限维函数空间。这个奇特的线性向量空间的一个特殊属性是它没有自由度：它不仅仅是有限维度的，它实际上是零维的，因此对于这种类型的对象， FiniteElement::n_dofs_per_cell() 将返回零。为了下面的讨论，让我们给这个空间一个合适的符号。

@f[
  Z = \{ \varphi: \varphi(x)=0 \}.


@f]

符号 $Z$ 提醒了这个空间的函数为零的事实。很明显，我们选择 $Z_h=Z$  。

对于我们用来描述弹性方程的变量，上面的整个讨论都可以重复。在这里，对于扩展变量，我们有

@f{align*}
  \tilde {\mathbf u} &\in U
   = \{\tilde {\mathbf u}|_{\Omega_s} \in H^1(\Omega_f)^d, \quad
       \tilde {\mathbf u}|_{\Omega_f} \in Z(\Omega_s)^d \},


@f}

而我们通常会使用这样的一个有限元空间

@f{align*}
  U_h
   &= \{{\mathbf u}_h \quad | \quad
       \forall K \in {\mathbb T}:
       {\mathbf u}_h|_K \in Q_r^d\  \text{if}\ K\subset {\Omega_s}, \quad
       {\mathbf u}_h|_{\Omega_f}\ \text{is continuous}, \quad
       {\mathbf u}_h|_K \in Z^d\ \text{if}\ K\subset {\Omega_f}\}
   && \subset U


@f}

的多项式程度  $r$  。

因此，总结起来，我们要在以下空间寻找一个离散的矢量值解 $y_h = \{\mathbf v_h, p_h, \mathbf u_h\}$ 。

@f{align*}
  Y_h = \{
      & y_h = \{\mathbf v_h, p_h, \mathbf u_h\} : \\
      & y_h|_{\Omega_f} \in Q_{p+1}^d \times Q_p \times Z^d, \\
      & y_h|_{\Omega_s} \in Z^d \times Z \times Q_r^d \}.


@f}






<a name="Implementation"></a><h3>Implementation</h3>


那么，我们如何实现这种事情呢？首先，我们意识到，离散空间 $Y_h$ 本质上需要两个不同的有限元。首先，在流体子域上，我们需要元素 $Q_{p+1}^d \times Q_p \times Z^d$ ，这在deal.II中很容易通过以下方式实现

@code
  FESystem<dim> (FE_Q<dim>(p+1), dim,
		 FE_Q<dim>(p), 1,
		 FE_Nothing<dim>(), dim),
@endcode

其中 <code>FE_Nothing</code> 实现了永远为零的函数的空间。其次，在实体子域上，我们需要元素 $\in Z^d \times Z \times Q_r^d$ ，我们用以下方法得到它

@code
  FESystem<dim> (FE_Nothing<dim>(), dim,
		 FE_Nothing<dim>(), 1,
		 FE_Q<dim>(r), dim),
@endcode



下一步是，我们将这两个元素中的每一个都与占据两个子域的细胞联系起来。为此，我们认识到，从某种意义上说，这两个元素只是彼此的变化，因为它们具有相同数量的向量分量，但具有不同的多项式度数&mdash；这很像人们在 $hp$ 有限元方法中的做法，这也正是我们在这里要做的：我们将（ab）使用hp-namespace的类和设施，将不同的元素分配给不同的单元。换句话说，我们将使用 hp::FECollection, 中的两个有限元与一个适当的 hp::QCollection 集成，使用 hp::FEValues 对象，而我们的DoFHandler将处于<i>hp</i>模式下。你不妨看一下步骤27，了解所有这些概念的概况。

在继续描述测试案例之前，让我们先澄清一下<i>why</i>这种将函数以零为单位扩展到整个领域，然后将问题映射到hp-framework上的方法是有意义的。

- 它使事情变得统一。在所有单元格中，向量分量的数量是相同的（这里是 <code>2*dim+1</code> ）。这使得各种事情都成为可能，因为统一的描述允许代码的重复使用。例如，计算每个向量分量的自由度 (DoFTools::count_dofs_per_fe_component), 按分量对自由度进行排序 (DoFRenumbering::component_wise), ，随后将矩阵和向量分割成块，以及其他许多函数都能一如既往地工作，而不需要为它们添加特殊的逻辑来描述某些变量只存在于部分领域的情况。因此，在像现在这样的程序中，你已经有了各种工具，这些工具最初并不是为多物理场情况编写的，但在目前的背景下却能正常工作。

- 它可以方便地进行图形化输出。我们支持的所有图形输出格式都要求输出中的每个字段都定义在网格的所有节点上。但是考虑到现在所有的解决方案组件都存在于各个地方，我们现有的DataOut例程可以像以前一样工作，并产生适合于可视化的图形输出--这些字段将简单地被扩展为0，如果不需要，可视化程序可以很容易地过滤掉这个值。

- 基本上没有成本。FE_Nothing的技巧并没有给整个问题增加任何自由度，我们也不需要处理属于这些分量的形状函数&mdash；FE_Nothing没有自由度，也没有形状函数，它所做的只是占用了矢量分量。




<a name="Specificsoftheimplementation"></a><h3> Specifics of the implementation </h3>


更具体地说，在该方案中，我们必须解决以下几点。

- 实现双线性形式，特别是处理界面项，在矩阵和稀疏模式中都是如此。

- 在边界的外部和内部部分实施迪里希特边界条件  $\partial\Omega_f,\partial\Omega_s$  。




<a name="Dealingwiththeinterfaceterms"></a><h4>Dealing with the interface terms</h4>


让我们首先讨论实现双线性形式，在离散水平上，我们记得它是

@f{align*}
	2 \eta (\varepsilon(\mathbf a_h), \varepsilon(\mathbf v_h))_{\Omega_f}


	- (\nabla \cdot \mathbf a_h, p_h)_{\Omega_f}


	- (q_h, \nabla \cdot \mathbf v_h)_{\Omega_f} &
	\\
	+ (\varepsilon(\mathbf b_h), C \varepsilon(\mathbf u_h))_{\Omega_s} &
	\\


	- (\mathbf b_h,
           (2 \eta \varepsilon(\mathbf v_h) - p \mathbf 1) \mathbf n)_{\Gamma_i}
	&=
	0,


@f}

鉴于我们已经将场扩展为零，原则上我们可以将子域上的积分写成整个域 $\Omega$ ，尽管在决定对哪些项进行积分之前，首先询问一个单元是弹性区域还是流体区域的一部分，这没有什么额外的努力。实际上，对这些项进行积分并不十分困难；对于斯托克斯方程，相关步骤已在步骤22中显示，而对于弹性方程，我们基本上采取 @ref vector_valued 模块中的形式（而不是步骤8中的形式）。

更值得关注的是界面术语。

@f[


	-(\mathbf b_h,
           (2 \eta \varepsilon(\mathbf v_h) - p \mathbf 1) \mathbf n)_{\Gamma_i}.


@f]

基于我们假设界面 $\Gamma_i$ 与细胞边界重合，这实际上可以写成一组面积分。如果我们用提取器符号 $\psi_i\in Y_h$ 表示形状函数 $\psi_i[\mathbf v],\psi_i[p], \psi_i[\mathbf u]$ 的速度、压力和位移分量，那么上述项就会产生对全局矩阵项 $i,j$ 的如下贡献。

@f[


	-\sum_K (\psi_i[\mathbf u],
           (2 \eta \varepsilon(\psi_j[\mathbf v]) - \psi_j[p] \mathbf 1)
	   \mathbf n)_{\partial K \cap \Gamma_i}.


@f]

虽然不是很明显，但这个术语带来了一点复杂的问题：虽然 $\psi_i[\mathbf u]$ 和 $\mathbf n$ 是在界面的实体一侧评估的（它们分别是位移和对 $\Omega_s$ 的法向量的测试函数，我们需要在界面的流体一侧评估 $\psi_j[\mathbf v],\psi_j[p]$ ，因为它们对应于流体施加的应力/力。换句话说，在我们的实现中，我们将需要界面两边的FEFaceValue对象。让事情变得更糟糕的是，我们可能还必须处理这样一个事实，即一方或另一方可能被细化，使我们需要在一个面的部分区域进行整合。请看下面的实现，如何处理这个问题。

作为一个额外的复杂问题，由这个术语产生的矩阵条目需要以某种方式添加到矩阵的稀疏模式中。这是DoFTools命名空间中各种函数的领域，比如 DoFTools::make_sparsity_pattern 和 DoFTools::make_flux_sparsity_pattern. 本质上，这些函数所做的是模拟系统矩阵装配过程中发生的事情：每当装配将非零条目写入全局矩阵，DoFTools中的函数将添加一个条目到稀疏模式中。因此，我们可以这样做：让 DoFTools::make_sparsity_pattern 将所有由常规的逐个单元积分产生的条目添加到稀疏性模式中，然后用手做同样的事情，即由接口项产生的条目。如果你看一下下面的程序中界面积分的实现，那么如何做应该是显而易见的，最多只需要不超过100行的代码。

但我们是懒人：界面项是沿一个面耦合两个相邻单元的自由度，这正是人们在非连续Galerkin方案中要做的事情，函数 DoFTools::make_flux_sparsity_pattern 就是为此而写。与通常的 DoFTools::make_sparsity_pattern: 相比，这是一个矩阵条目的超集，它还将添加所有计算来自所有面的两侧自由度的耦合项的条目。不幸的是，对于这个函数的最简单版本，这是一个相当大的超集。例如，考虑下面这个有两个单元和一个 $Q_1$ 有限元的网格。

@code
  2---3---5
  |   |   |
  0---1---4
@endcode

这里，由 DoFTools::make_sparsity_pattern 产生的稀疏模式将只有在一个单元上耦合的自由度的条目。然而，它不会有稀疏模式条目 $(0,4),(0,5),(2,4),(2,5)$  。然而，由 DoFTools::make_flux_sparsity_pattern 生成的稀疏模式将有这些条目：它假定你想为一个双线性形式建立一个稀疏模式，该形式将<i>all</i>自由度从相邻的单元上耦合起来。这不是我们想要的：我们的界面项只作用于一小部分单元，我们当然不需要两个相邻流体单元或两个相邻固体单元之间的所有额外耦合。此外，我们使用高阶元素的事实意味着我们确实会产生比实际需要多得多的条目：在最粗的网格上，在2D中，44,207个非零条目而不是16,635个 DoFTools::make_sparsity_pattern, ，导致我们后来建立的矩阵中出现大量的零（当然，16,635是不够的，因为它们不包括界面条目）。这个比例在3D中会更糟糕。

所以极度懒惰是有代价的：矩阵中的条目太多。但我们可以适度偷懒：有一个 DoFTools::make_flux_sparsity_pattern 的变体，允许我们指定有限元的哪些矢量分量与哪些其他分量耦合，既可以用单元术语，也可以用面术语。对于实体子域中的单元，我们将所有位移相互耦合；对于流体单元，所有速度与所有速度和压力耦合，但压力不与自身耦合。由于没有一个单元同时拥有两组变量，因此没有必要区分这两种单元，所以我们可以这样写掩码。

@code
    Table<2,DoFTools::Coupling> cell_coupling (fe_collection.n_components(),
					       fe_collection.n_components());


    for (unsigned int c=0; c<fe_collection.n_components(); ++c)
      for (unsigned int d=0; d<fe_collection.n_components(); ++d)
	if (((c<dim+1) && (d<dim+1)
	     && !((c==dim) && (d==dim)))
	    ||
	    ((c>=dim+1) && (d>=dim+1)))
	  cell_coupling[c][d] = DoFTools::Coupling::always;
@endcode

在这里，我们使用了这样一个事实：有限元的第一个 <code>dim</code> 分量是速度，然后是压力，最后是 <code>dim</code> 位移。(我们也可以说，速度/压力也与位移耦合，因为没有一个单元同时拥有两组变量)。另一方面，界面条款需要一个像这样的掩码。

@code
    Table<2,DoFTools::Coupling> face_coupling (fe_collection.n_components(),
					       fe_collection.n_components());


    for (unsigned int c=0; c<fe_collection.n_components(); ++c)
      for (unsigned int d=0; d<fe_collection.n_components(); ++d)
	if ((c>=dim+1) && (d<dim+1))
	  face_coupling[c][d] = DoFTools::Coupling::always;
@endcode

换句话说，所有的位移测试函数（组件 <code>c@>=dim+1</code> ）与界面另一侧的所有速度和压力形状函数耦合。这并不完全正确，尽管很接近：事实上，界面的确切形式仅指那些在共同界面上确实为非零的压力位移形状函数，这对所有形状函数来说并不正确；另一方面，它确实耦合了所有速度（因为积分涉及速度形状函数的梯度，这些梯度在单元的所有面上均为非零）。然而，我们在上面建立的掩码，并不具备这些微妙的能力。尽管如此，通过这些掩码，我们设法将稀疏模式的条目数降低到21028个&mdash；目前来说已经足够了。




<a name="Velocityboundaryconditionsontheinterface"></a><h4>Velocity boundary conditions on the interface</h4>


第二个困难是，虽然我们知道如何在外部边界上强制执行速度或应力为零（使用 VectorTools::interpolate_boundary_values, 调用适当的分量掩码，并为固体和流体外部边界设置不同的边界指标），但现在我们还需要在内部界面上的速度为零，即 $\mathbf v|_{\Gamma_i}=0$  。在写这篇文章时，deal.II中没有处理这部分的函数，但用手实现并不特别困难：基本上，我们只需要在所有单元上循环，如果它是一个流体单元，而它的邻居是一个固体单元，然后添加约束，确保这个面上的速度自由度为零。在处理相邻的固体单元被细化的情况下，有必要进行一些处理，产生以下代码。

@code
std::vector<unsigned int> local_face_dof_indices (stokes_fe.dofs_per_face);
for (const auto &cell: dof_handler.active_cell_iterators())
  if (cell_is_in_fluid_domain (cell))
    for (const auto f : cell->face_indices())
      if (!cell->at_boundary(f))
        {
          bool face_is_on_interface = false;


          if ((cell->neighbor(f)->has_children() == false)
	          &&
	          (cell_is_in_solid_domain (cell->neighbor(f))))
	        face_is_on_interface = true;
          else if (cell->neighbor(f)->has_children() == true)
	        {
              // The neighbor does have children. See if any of the cells
              // on the other side are elastic
	          for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf)
	            if (cell_is_in_solid_domain (cell->neighbor_child_on_subface(f, sf)))
	              {
                   face_is_on_interface = true;
		            break;
	              }
	        }


          if (face_is_on_interface)
           {
             cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
             for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
             if (stokes_fe.face_system_to_component_index(i).first < dim)
               constraints.add_line (local_face_dof_indices[i]);
           }
        }
@endcode



调用 <code>constraints.add_line(t)</code> 告诉AffineConstraints为自由度 <code>t</code> 启动一个新的约束，其形式为 $x_t=\sum_{l=0}^{N-1} c_{tl} x_l +
b_t$  。通常情况下，我们会将单个系数 $c_{tl}$ 设置为非零值（使用 AffineConstraints::add_entry) 或将 $b_t$ 设置为非零值（使用 AffineConstraints::set_inhomogeneity); 像上面那样什么都不做，虽然看起来很有趣，但只是让约束成为 $x_t=0$ ，这正是我们在当前情况下需要的。对 FiniteElement::face_system_to_component_index 的调用确保了我们只将速度分量的边界值设置为零，而不是压力分量。

请注意，在有些情况下，这可能会产生不正确的结果：特别是，一旦我们找到当前流体单元的一个实体邻接子，我们就会假设共同面上的所有邻接子都在实体子域。但事实并非如此，例如，考虑以下的网格。

@code
+---------+----+----+
|         | f  |    |
|    f    +----+----+
|         | s  |    |
+---------+----+----+
@endcode



在这种情况下，我们将把左单元右面的所有速度自由度设置为零，这对该面的顶部自由度来说是不正确的。也就是说，只有当流体和固体子域不与一组完整的粗网格单元重合时才会发生这种情况&mdash;但这与本介绍第一节末尾所述的假设是矛盾的。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们将考虑以下情况作为一个测试案例。

 <img src="https://www.dealii.org/images/steps/developer/step-46.layout.png" alt=""> 

正如本文顶部所讨论的，我们需要在一些地方假设一个单元完全处于域的流体部分或固体部分，此外，一个不活动单元的所有子域也属于同一个子域。如果粗略网格已经将网格细分为实体和流体粗略网格单元，这一点肯定可以得到保证；考虑到上面概述的几何形状，我们可以通过使用 $8\times 8$ 粗略网格，方便地提供 GridGenerator::subdivided_hyper_rectangle 函数来实现。

底部的固定边界意味着 $\mathbf u=0$ ，我们也为顶部的流动规定了迪里希特条件，因此我们在左边得到流入，在右边得到流出。在左边和右边的边界，没有对流动施加明确的边界条件，产生隐性的无应力条件  $(2\eta
\varepsilon(\mathbf v) - p \mathbf 1) \cdot \mathbf n = 0$  。上面已经讨论了两个域之间的界面条件。

为了简单起见，我们选择材料参数为 $\eta=\lambda=\mu=1$  。在下面的结果部分，我们还将展示一个可以从同一程序中获得的三维模拟。边界条件和几何形状的定义几乎与上面的2d情况类似。




<a name="Identifyingwhichsubdomainacellisin"></a><h4>Identifying which subdomain a cell is in</h4>


在程序中，我们需要一种方法来识别一个细胞处于域的哪一部分。有许多不同的方法可以做到这一点。一个典型的方法是使用每个单元的 @ref GlossSubdomainId "subdomain_id "标签，尽管这个字段在%并行计算中具有特殊意义。另一种方法是 @ref GlossMaterialId "material_id "字段，也是每个单元格都有的。它有一个额外的优点，就是在网格细化时，它可以从母体继承到子体；换句话说，我们在创建网格时设置一次材料ID，即使经过几次细化循环，它对所有活动单元都是正确的。因此，我们采用这种方法：我们定义一个 <code>enum</code> ，用符号名称来表示材料ID的数字，并使用它们来识别单元在域的哪一部分。

其次，我们使用一个在<i>hp</i>模式下运行的DoFHandler类型的对象。该类需要知道哪些单元将使用斯托克斯有限元，哪些使用弹性有限元。因此，在每个细化周期的开始，我们必须走过所有的单元，并将（在hp-parlance中）活动FE索引设置为任何适合当前情况的索引。虽然我们可以使用符号名称来表示材料ID，但主动FE索引实际上是一个数字，经常用于索引对象的集合（例如 hp::FECollection 和 hp::QCollection); 类型），这意味着主动FE索引实际上对于领域的流体部分必须是0，对于弹性部分必须是1。




<a name="Linearsolvers"></a><h4>Linear solvers</h4>


这个程序主要是为了展示如何处理领域内不同部分的不同物理现象，以及如何在deal.II中实现这样的模型。因此，我们不会费力想出一个好的求解器：我们只是使用SparseDirectUMFPACK类，它总是有效的，即使不是最佳的复杂性。然而，我们将在<a href="#Results">results</a>部分对可能的其他求解器进行评论。




<a name="Meshrefinement"></a><h4>Mesh refinement</h4>


这个程序的一个比较棘手的方面是如何估计误差。因为它几乎适用于任何程序，所以我们想使用KellyErrorEstimator，在这里我们也可以用下面这样的代码相对容易地做到。

@code
  Vector<float> stokes_estimated_error_per_cell (triangulation.n_active_cells());
  Vector<float> elasticity_estimated_error_per_cell (triangulation.n_active_cells());


  std::vector<bool> stokes_component_mask (dim+1+dim, false);
  for (unsigned int d=0; d<dim; ++d)
    stokes_component_mask[d] = true;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      face_q_collection,
                                      std::map<types::boundary_id, const Function<dim>*>(),
                                      solution,
                                      stokes_estimated_error_per_cell,
                                      stokes_component_mask);


  std::vector<bool> elasticity_component_mask (dim+1+dim, false);
  for (unsigned int d=0; d<dim; ++d)
    elasticity_component_mask[dim+1+d] = true;
  KellyErrorEstimator<dim>::estimate (dof_handler,
                                      face_q_collection,
                                      std::map<types::boundary_id, const Function<dim>*>(),
                                      solution,
                                      elasticity_estimated_error_per_cell,
                                      elasticity_component_mask);
@endcode

这就为每个单元提供了两套误差指标。然后我们会以某种方式将它们合并成一个用于网格细化，例如使用类似下面的方法（注意，我们将两个向量中的平方误差指标归一化，因为误差量的物理单位在当前情况下并不匹配，导致两个子域之间的误差指标可能存在数量级的差异）。

@code
  stokes_estimated_error_per_cell /= stokes_estimated_error_per_cell.l2_norm();
  elasticity_estimated_error_per_cell /= elasticity_estimated_error_per_cell.l2_norm();


  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  estimated_error_per_cell += stokes_estimated_error_per_cell;
  estimated_error_per_cell += elasticity_estimated_error_per_cell;
@endcode

(在代码中，我们实际上以4:1的比例权衡误差指标，以支持在斯托克斯子域上计算的误差指标，因为细化在其他方面严重偏向弹性子域，但这只是一个技术问题。因素4已经被启发式地确定为相当好的工作。)

虽然这个原则是合理的，但它并不完全像预期的那样工作。原因是KellyErrorEstimator类是通过整合每个单元面周围的解的梯度跳跃来计算误差指标。这个跳跃在解不连续和扩展为零的地方可能非常大；它也不会随着网格的细化而变小。KellyErrorEstimator类不能忽视这个接口，因为它基本上只看到<i>hp</i>模式下的DoFHandler，其中元素类型从一个单元改变到另一个单元&mdash；正是<i>hp</i>模式所设计的东西，当前程序中的接口看起来与步骤27中的接口没有什么不同，例如，当然也没有更合理的。尽管如此，最终的结果是，在两个子域之间的界面两侧都有一层单元，其误差指标大得不合理。因此，大部分的网格细化工作都集中在界面上。

如果我们有一个真正理解问题的细化指标，并且在积分跳跃项时简单地忽略子域之间的界面，这显然就不会发生。另一方面，这个程序是关于展示如何表示我们在不同子域有不同物理学的问题，而不是关于KellyErrorEstimator的特殊性，因此我们诉诸于被称为 "启发式 "的大锤子：我们简单地把界面上的单元格的误差指标设置为零。这就切断了误差指标中的尖峰。乍看之下，人们也会认为它阻止了网格在界面上的细化，但是相邻的单元只能有一级细化的要求，仍然会导致一个合理的细化网格。

虽然这显然是一个次优的解决方案，但它目前是可行的，并为未来的改进留下了空间。


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
 * 这个程序的包含文件与之前许多其他程序的包含文件是一样的。唯一的新文件是在介绍中讨论的声明FE_Nothing的文件。hp目录下的文件已经在  step-27  中讨论过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_nothing.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/hp/fe_collection.h> 
 * #include <deal.II/hp/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * 
 * namespace Step46 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeFluidStructureProblemcodeclasstemplate"></a> 
 * <h3>The <code>FluidStructureProblem</code> class template</h3>
 * 

 * 
 * 这是主类。如果你想的话，它是 step-8 和 step-22 的组合，因为它的成员变量要么针对全局问题（Triangulation和DoFHandler对象，以及 hp::FECollection 和各种线性代数对象），要么与弹性或斯托克斯子问题有关。然而，该类的一般结构与其他大多数实现静止问题的程序一样。
 * 

 * 
 * 有几个不言自明的辅助函数（<code>cell_is_in_fluid_domain, cell_is_in_solid_domain</code>）（对两个子域的符号名称进行操作，这些名称将被用作属于子域的单元的 material_ids。正如介绍中所解释的那样）和几个函数（<code>make_grid, set_active_fe_indices, assemble_interface_terms</code>），这些函数已经从其他的函数中分离出来，可以在其他的教程程序中找到，我们将在实现它们的时候讨论。
 * 

 * 
 * 最后一组变量 (  <code>viscosity, lambda, eta</code>  ) 描述了用于两个物理模型的材料属性。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class FluidStructureProblem 
 *   { 
 *   public: 
 *     FluidStructureProblem(const unsigned int stokes_degree, 
 *                           const unsigned int elasticity_degree); 
 *     void run(); 
 * 
 *   private: 
 *     enum 
 *     { 
 *       fluid_domain_id, 
 *       solid_domain_id 
 *     }; 
 * 
 *     static bool cell_is_in_fluid_domain( 
 *       const typename DoFHandler<dim>::cell_iterator &cell); 
 * 
 *     static bool cell_is_in_solid_domain( 
 *       const typename DoFHandler<dim>::cell_iterator &cell); 
 * 
 *     void make_grid(); 
 *     void set_active_fe_indices(); 
 *     void setup_dofs(); 
 *     void assemble_system(); 
 *     void assemble_interface_term( 
 *       const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
 *       const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
 *       std::vector<Tensor<1, dim>> &         elasticity_phi, 
 *       std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
 *       std::vector<double> &                 stokes_phi_p, 
 *       FullMatrix<double> &                  local_interface_matrix) const; 
 *     void solve(); 
 *     void output_results(const unsigned int refinement_cycle) const; 
 *     void refine_mesh(); 
 * 
 *     const unsigned int stokes_degree; 
 *     const unsigned int elasticity_degree; 
 * 
 *     Triangulation<dim>    triangulation; 
 *     FESystem<dim>         stokes_fe; 
 *     FESystem<dim>         elasticity_fe; 
 *     hp::FECollection<dim> fe_collection; 
 *     DoFHandler<dim>       dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 * 
 *     const double viscosity; 
 *     const double lambda; 
 *     const double mu; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Boundaryvaluesandrighthandside"></a> 
 * <h3>Boundary values and right hand side</h3>
 * 

 * 
 * 下面这个类如其名。速度的边界值分别为2d的 
 * $\mathbf u=(0, \sin(\pi x))^T$ 和3d的 $\mathbf u=(0,
 * 0, \sin(\pi x)\sin(\pi y))^T$ 。
 * 这个问题的其余边界条件都是同质的，在介绍中已经讨论过。右边的强迫项对于流体和固体都是零，所以我们不需要为它设置额外的类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class StokesBoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     StokesBoundaryValues() 
 *       : Function<dim>(dim + 1 + dim) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  value) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double StokesBoundaryValues<dim>::value(const Point<dim> & p, 
 *                                           const unsigned int component) const 
 *   { 
 *     Assert(component < this->n_components, 
 *            ExcIndexRange(component, 0, this->n_components)); 
 * 
 *     if (component == dim - 1) 
 *       switch (dim) 
 *         { 
 *           case 2: 
 *             return std::sin(numbers::PI * p[0]); 
 *           case 3: 
 *             return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
 *           default: 
 *             Assert(false, ExcNotImplemented()); 
 *         } 
 * 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesBoundaryValues<dim>::vector_value(const Point<dim> &p, 
 *                                                Vector<double> &  values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = StokesBoundaryValues<dim>::value(p, c); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeFluidStructureProblemcodeimplementation"></a> 
 * <h3>The <code>FluidStructureProblem</code> implementation</h3>
 * 
 * <a name="Constructorsandhelperfunctions"></a> 
 * <h4>Constructors and helper functions</h4>
 * 

 * 
 * 现在我们来谈谈这个程序的主类的实现。最初的几个函数是构造函数和辅助函数，可以用来确定一个单元格在域的哪个部分。鉴于介绍中对这些主题的讨论，它们的实现是相当明显的。在构造函数中，注意我们必须从斯托克斯和弹性的基本元素中构造 hp::FECollection 对象；使用 hp::FECollection::push_back 函数在这个集合中为它们分配了0和1的位置，我们必须记住这个顺序，并在程序的其余部分一致使用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   FluidStructureProblem<dim>::FluidStructureProblem( 
 *     const unsigned int stokes_degree, 
 *     const unsigned int elasticity_degree) 
 *     : stokes_degree(stokes_degree) 
 *     , elasticity_degree(elasticity_degree) 
 *     , triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , stokes_fe(FE_Q<dim>(stokes_degree + 1), 
 *                 dim, 
 *                 FE_Q<dim>(stokes_degree), 
 *                 1, 
 *                 FE_Nothing<dim>(), 
 *                 dim) 
 *     , elasticity_fe(FE_Nothing<dim>(), 
 *                     dim, 
 *                     FE_Nothing<dim>(), 
 *                     1, 
 *                     FE_Q<dim>(elasticity_degree), 
 *                     dim) 
 *     , dof_handler(triangulation) 
 *     , viscosity(2) 
 *     , lambda(1) 
 *     , mu(1) 
 *   { 
 *     fe_collection.push_back(stokes_fe); 
 *     fe_collection.push_back(elasticity_fe); 
 *   } 
 * 
 *   template <int dim> 
 *   bool FluidStructureProblem<dim>::cell_is_in_fluid_domain( 
 *     const typename DoFHandler<dim>::cell_iterator &cell) 
 *   { 
 *     return (cell->material_id() == fluid_domain_id); 
 *   } 
 * 
 *   template <int dim> 
 *   bool FluidStructureProblem<dim>::cell_is_in_solid_domain( 
 *     const typename DoFHandler<dim>::cell_iterator &cell) 
 *   { 
 *     return (cell->material_id() == solid_domain_id); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Meshesandassigningsubdomains"></a> 
 * <h4>Meshes and assigning subdomains</h4>
 * 

 * 
 * 接下来的一对函数是处理生成网格，并确保所有表示子域的标志都是正确的。  <code>make_grid</code>  ，正如在介绍中所讨论的，生成一个 $8\times 8$ 的网格（或者一个 $8\times 8\times 8$ 的三维网格）以确保每个粗略的网格单元完全在一个子域内。生成这个网格后，我们在其边界上循环，并在顶部边界设置边界指标为1，这是我们设置非零迪里希特边界条件的唯一地方。在这之后，我们再次在所有单元上循环，设置材料指标&mdash;用来表示我们处于域的哪一部分，是流体指标还是固体指标。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::make_grid() 
 *   { 
 *     GridGenerator::subdivided_hyper_cube(triangulation, 8, -1, 1); 
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->at_boundary() && (face->center()[dim - 1] == 1)) 
 *           face->set_all_boundary_ids(1); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (((std::fabs(cell->center()[0]) < 0.25) && 
 *            (cell->center()[dim - 1] > 0.5)) || 
 *           ((std::fabs(cell->center()[0]) >= 0.25) && 
 *            (cell->center()[dim - 1] > -0.5))) 
 *         cell->set_material_id(fluid_domain_id); 
 *       else 
 *         cell->set_material_id(solid_domain_id); 
 *   } 
 * 
 * @endcode
 * 
 * 这对函数的第二部分决定在每个单元上使用哪个有限元。上面我们设置了每个粗略网格单元的材料指标，正如在介绍中提到的，这个信息在网格细化时将从母单元继承到子单元。
 * 

 * 
 * 换句话说，只要我们细化（或创建）了网格，我们就可以依靠材料指示器来正确描述一个单元所处的域的哪一部分。然后我们利用这一点将单元的活动FE索引设置为该类的 hp::FECollection 成员变量中的相应元素：流体单元为0，固体单元为1。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::set_active_fe_indices() 
 *   { 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         if (cell_is_in_fluid_domain(cell)) 
 *           cell->set_active_fe_index(0); 
 *         else if (cell_is_in_solid_domain(cell)) 
 *           cell->set_active_fe_index(1); 
 *         else 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemsetup_dofscode"></a> 
 * <h4><code>FluidStructureProblem::setup_dofs</code></h4>
 * 

 * 
 * 下一步是为线性系统设置数据结构。为此，我们首先要用上面的函数设置活动FE指数，然后分配自由度，再确定线性系统的约束。后者包括像往常一样的悬挂节点约束，但也包括顶部流体边界的不均匀边界值，以及沿固体子域周边的零边界值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::setup_dofs() 
 *   { 
 *     set_active_fe_indices(); 
 *     dof_handler.distribute_dofs(fe_collection); 
 * 
 *     { 
 *       constraints.clear(); 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 * 
 *       const FEValuesExtractors::Vector velocities(0); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                1, 
 *                                                StokesBoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe_collection.component_mask( 
 *                                                  velocities)); 
 * 
 *       const FEValuesExtractors::Vector displacements(dim + 1); 
 *       VectorTools::interpolate_boundary_values( 
 *         dof_handler, 
 *         0, 
 *         Functions::ZeroFunction<dim>(dim + 1 + dim), 
 *         constraints, 
 *         fe_collection.component_mask(displacements)); 
 *     } 
 * 
 * @endcode
 * 
 * 不过，我们还需要处理更多的约束条件：我们必须确保在流体和固体的界面上速度为零。下面这段代码已经在介绍中介绍过了。
 * 

 * 
 * 
 * @code
 *     { 
 *       std::vector<types::global_dof_index> local_face_dof_indices( 
 *         stokes_fe.n_dofs_per_face()); 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         if (cell_is_in_fluid_domain(cell)) 
 *           for (const auto face_no : cell->face_indices()) 
 *             if (cell->face(face_no)->at_boundary() == false) 
 *               { 
 *                 bool face_is_on_interface = false; 
 * 
 *                 if ((cell->neighbor(face_no)->has_children() == false) && 
 *                     (cell_is_in_solid_domain(cell->neighbor(face_no)))) 
 *                   face_is_on_interface = true; 
 *                 else if (cell->neighbor(face_no)->has_children() == true) 
 *                   { 
 *                     for (unsigned int sf = 0; 
 *                          sf < cell->face(face_no)->n_children(); 
 *                          ++sf) 
 *                       if (cell_is_in_solid_domain( 
 *                             cell->neighbor_child_on_subface(face_no, sf))) 
 *                         { 
 *                           face_is_on_interface = true; 
 *                           break; 
 *                         } 
 *                   } 
 * 
 *                 if (face_is_on_interface) 
 *                   { 
 *                     cell->face(face_no)->get_dof_indices(local_face_dof_indices, 
 *                                                          0); 
 *                     for (unsigned int i = 0; i < local_face_dof_indices.size(); 
 *                          ++i) 
 *                       if (stokes_fe.face_system_to_component_index(i).first < 
 *                           dim) 
 *                         constraints.add_line(local_face_dof_indices[i]); 
 *                   } 
 *               } 
 *     } 
 * 
 * @endcode
 * 
 * 在这一切结束后，我们可以向约束对象声明，我们现在已经准备好了所有的约束，并且该对象可以重建其内部数据结构以提高效率。
 * 

 * 
 * 
 * @code
 *     constraints.close(); 
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 * @endcode
 * 
 * 在这个函数的其余部分，我们创建了一个在介绍中广泛讨论的稀疏模式，并使用它来初始化矩阵；然后还将向量设置为正确的大小。
 * 

 * 
 * 
 * @code
 *     { 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 * 
 *       Table<2, DoFTools::Coupling> cell_coupling(fe_collection.n_components(), 
 *                                                  fe_collection.n_components()); 
 *       Table<2, DoFTools::Coupling> face_coupling(fe_collection.n_components(), 
 *                                                  fe_collection.n_components()); 
 * 
 *       for (unsigned int c = 0; c < fe_collection.n_components(); ++c) 
 *         for (unsigned int d = 0; d < fe_collection.n_components(); ++d) 
 *           { 
 *             if (((c < dim + 1) && (d < dim + 1) && 
 *                  !((c == dim) && (d == dim))) || 
 *                 ((c >= dim + 1) && (d >= dim + 1))) 
 *               cell_coupling[c][d] = DoFTools::always; 
 * 
 *             if ((c >= dim + 1) && (d < dim + 1)) 
 *               face_coupling[c][d] = DoFTools::always; 
 *           } 
 * 
 *       DoFTools::make_flux_sparsity_pattern(dof_handler, 
 *                                            dsp, 
 *                                            cell_coupling, 
 *                                            face_coupling); 
 *       constraints.condense(dsp); 
 *       sparsity_pattern.copy_from(dsp); 
 *     } 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemassemble_systemcode"></a> 
 * <h4><code>FluidStructureProblem::assemble_system</code></h4>
 * 

 * 
 * 下面是这个程序的中心函数：组装线性系统的函数。它在开始时有一长段设置辅助函数的内容：从创建正交公式到设置FEValues、FEFaceValues和FESubfaceValues对象，这些都是整合单元项以及界面项所必需的，以应对界面上的单元以相同大小或不同细化程度聚集在一起的情况...
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::assemble_system() 
 *   { 
 *     system_matrix = 0; 
 *     system_rhs    = 0; 
 * 
 *     const QGauss<dim> stokes_quadrature(stokes_degree + 2); 
 *     const QGauss<dim> elasticity_quadrature(elasticity_degree + 2); 
 * 
 *     hp::QCollection<dim> q_collection; 
 *     q_collection.push_back(stokes_quadrature); 
 *     q_collection.push_back(elasticity_quadrature); 
 * 
 *     hp::FEValues<dim> hp_fe_values(fe_collection, 
 *                                    q_collection, 
 *                                    update_values | update_quadrature_points | 
 *                                      update_JxW_values | update_gradients); 
 * 
 *     const QGauss<dim - 1> common_face_quadrature( 
 *       std::max(stokes_degree + 2, elasticity_degree + 2)); 
 * 
 *     FEFaceValues<dim>    stokes_fe_face_values(stokes_fe, 
 *                                             common_face_quadrature, 
 *                                             update_JxW_values | 
 *                                               update_gradients | update_values); 
 *     FEFaceValues<dim>    elasticity_fe_face_values(elasticity_fe, 
 *                                                 common_face_quadrature, 
 *                                                 update_normal_vectors | 
 *                                                   update_values); 
 *     FESubfaceValues<dim> stokes_fe_subface_values(stokes_fe, 
 *                                                   common_face_quadrature, 
 *                                                   update_JxW_values | 
 *                                                     update_gradients | 
 *                                                     update_values); 
 *     FESubfaceValues<dim> elasticity_fe_subface_values(elasticity_fe, 
 *                                                       common_face_quadrature, 
 *                                                       update_normal_vectors | 
 *                                                         update_values); 
 * 
 * @endcode
 * 
 * ...描述局部对全局线性系统贡献所需的对象...
 * 

 * 
 * 
 * @code
 *     const unsigned int stokes_dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
 *     const unsigned int elasticity_dofs_per_cell = 
 *       elasticity_fe.n_dofs_per_cell(); 
 * 
 *     FullMatrix<double> local_matrix; 
 *     FullMatrix<double> local_interface_matrix(elasticity_dofs_per_cell, 
 *                                               stokes_dofs_per_cell); 
 *     Vector<double>     local_rhs; 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 *     std::vector<types::global_dof_index> neighbor_dof_indices( 
 *       stokes_dofs_per_cell); 
 * 
 *     const Functions::ZeroFunction<dim> right_hand_side(dim + 1); 
 * 
 * @endcode
 * 
 * ...到变量，允许我们提取形状函数的某些成分并缓存它们的值，而不是在每个正交点重新计算它们。
 * 

 * 
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 *     const FEValuesExtractors::Vector displacements(dim + 1); 
 * 
 *     std::vector<SymmetricTensor<2, dim>> stokes_symgrad_phi_u( 
 *       stokes_dofs_per_cell); 
 *     std::vector<double> stokes_div_phi_u(stokes_dofs_per_cell); 
 *     std::vector<double> stokes_phi_p(stokes_dofs_per_cell); 
 * 
 *     std::vector<Tensor<2, dim>> elasticity_grad_phi(elasticity_dofs_per_cell); 
 *     std::vector<double>         elasticity_div_phi(elasticity_dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> elasticity_phi(elasticity_dofs_per_cell); 
 * 
 * @endcode
 * 
 * 然后是所有单元格的主循环，和 step-27 一样，初始化当前单元格的 hp::FEValues 对象，提取适合当前单元格的FEValues对象。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         hp_fe_values.reinit(cell); 
 * 
 *         const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values(); 
 * 
 *         local_matrix.reinit(cell->get_fe().n_dofs_per_cell(), 
 *                             cell->get_fe().n_dofs_per_cell()); 
 *         local_rhs.reinit(cell->get_fe().n_dofs_per_cell()); 
 * 
 * @endcode
 * 
 * 做完这些后，我们继续为属于斯托克斯和弹性区域的单元组装单元项。虽然我们原则上可以在一个公式中完成，实际上就是实现了介绍中所说的双线性形式，但我们意识到，我们的有限元空间的选择方式是，在每个单元上，有一组变量（速度和压力，或者位移）总是为零，因此，计算局部积分的更有效的方法是，根据测试我们处于域的哪一部分的 <code>if</code> 条款，只做必要的事情。
 * 

 * 
 * 局部矩阵的实际计算与 step-22 以及 @ref vector_valued 文件模块中给出的弹性方程的计算相同。
 * 

 * 
 * 
 * @code
 *         if (cell_is_in_fluid_domain(cell)) 
 *           { 
 *             const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
 *             Assert(dofs_per_cell == stokes_dofs_per_cell, ExcInternalError()); 
 * 
 *             for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
 *               { 
 *                 for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                   { 
 *                     stokes_symgrad_phi_u[k] = 
 *                       fe_values[velocities].symmetric_gradient(k, q); 
 *                     stokes_div_phi_u[k] = 
 *                       fe_values[velocities].divergence(k, q); 
 *                     stokes_phi_p[k] = fe_values[pressure].value(k, q); 
 *                   } 
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     local_matrix(i, j) += 
 *                       (2 * viscosity * stokes_symgrad_phi_u[i] * 
 *                          stokes_symgrad_phi_u[j] - 
 *                        stokes_div_phi_u[i] * stokes_phi_p[j] - 
 *                        stokes_phi_p[i] * stokes_div_phi_u[j]) * 
 *                       fe_values.JxW(q); 
 *               } 
 *           } 
 *         else 
 *           { 
 *             const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
 *             Assert(dofs_per_cell == elasticity_dofs_per_cell, 
 *                    ExcInternalError()); 
 * 
 *             for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
 *               { 
 *                 for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                   { 
 *                     elasticity_grad_phi[k] = 
 *                       fe_values[displacements].gradient(k, q); 
 *                     elasticity_div_phi[k] = 
 *                       fe_values[displacements].divergence(k, q); 
 *                   } 
 * 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     { 
 *                       local_matrix(i, j) += 
 *                         (lambda * elasticity_div_phi[i] * 
 *                            elasticity_div_phi[j] + 
 *                          mu * scalar_product(elasticity_grad_phi[i], 
 *                                              elasticity_grad_phi[j]) + 
 *                          mu * 
 *                            scalar_product(elasticity_grad_phi[i], 
 *                                           transpose(elasticity_grad_phi[j]))) * 
 *                         fe_values.JxW(q); 
 *                     } 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 一旦我们得到了单元积分的贡献，我们就把它们复制到全局矩阵中（通过 AffineConstraints::distribute_local_to_global 函数，立即处理约束）。请注意，我们没有向 <code>local_rhs</code> 变量中写入任何东西，尽管我们仍然需要传递它，因为消除非零边界值需要修改局部，因此也需要修改全局的右手值。
 * 

 * 
 * 
 * @code
 *         local_dof_indices.resize(cell->get_fe().n_dofs_per_cell()); 
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global(local_matrix, 
 *                                                local_rhs, 
 *                                                local_dof_indices, 
 *                                                system_matrix, 
 *                                                system_rhs); 
 * 
 * @endcode
 * 
 * 这个函数更有趣的部分是我们看到关于两个子域之间的界面上的脸部条款。为此，我们首先要确保我们只组装一次，即使在所有单元的所有面的循环中会遇到界面的每一部分两次。我们武断地决定，只有当当前单元是固体子域的一部分，并且因此一个面不在边界上，并且它后面的潜在邻居是流体域的一部分时，我们才会评估界面条款。让我们从这些条件开始。
 * 

 * 
 * 
 * @code
 *         if (cell_is_in_solid_domain(cell)) 
 *           for (const auto f : cell->face_indices()) 
 *             if (cell->face(f)->at_boundary() == false) 
 *               { 
 * 
 * @endcode
 * 
 * 在这一点上，我们知道当前的单元格是一个候选的整合对象，并且面 <code>f</code> 后面存在一个邻居。现在有三种可能性。           
 * 

 * 
 * - 邻居处于同一细化水平，并且没有孩子。     
 * 

 * 
 * - 邻居有子女。     
 * 

 * 
 * - 邻居比较粗糙。            在所有这三种情况下，我们只对它感兴趣，如果它是流体子域的一部分。因此，让我们从第一种最简单的情况开始：如果邻居处于同一层次，没有子女，并且是一个流体单元，那么这两个单元共享一个边界，这个边界是界面的一部分，我们想沿着这个边界整合界面项。我们所要做的就是用当前面和邻接单元的面初始化两个FEFaceValues对象（注意我们是如何找出邻接单元的哪个面与当前单元接壤的），然后把东西传给评估界面项的函数（这个函数的第三个到第五个参数为它提供了抓取数组）。然后，结果再次被复制到全局矩阵中，使用一个知道本地矩阵的行和列的DoF指数来自不同单元的函数。
 * 

 * 
 * 
 * @code
 *                 if ((cell->neighbor(f)->level() == cell->level()) && 
 *                     (cell->neighbor(f)->has_children() == false) && 
 *                     cell_is_in_fluid_domain(cell->neighbor(f))) 
 *                   { 
 *                     elasticity_fe_face_values.reinit(cell, f); 
 *                     stokes_fe_face_values.reinit(cell->neighbor(f), 
 *                                                  cell->neighbor_of_neighbor(f)); 
 * 
 *                     assemble_interface_term(elasticity_fe_face_values, 
 *                                             stokes_fe_face_values, 
 *                                             elasticity_phi, 
 *                                             stokes_symgrad_phi_u, 
 *                                             stokes_phi_p, 
 *                                             local_interface_matrix); 
 * 
 *                     cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
 *                     constraints.distribute_local_to_global( 
 *                       local_interface_matrix, 
 *                       local_dof_indices, 
 *                       neighbor_dof_indices, 
 *                       system_matrix); 
 *                   } 
 * 
 * @endcode
 * 
 * 第二种情况是，如果邻居还有更多的孩子。在这种情况下，我们必须在邻居的所有子女中进行循环，看他们是否属于流体子域的一部分。如果它们是，那么我们就在共同界面上进行整合，这个界面是邻居的一个面和当前单元的一个子面，要求我们对邻居使用FEFaceValues，对当前单元使用FESubfaceValues。
 * 

 * 
 * 
 * @code
 *                 else if ((cell->neighbor(f)->level() == cell->level()) && 
 *                          (cell->neighbor(f)->has_children() == true)) 
 *                   { 
 *                     for (unsigned int subface = 0; 
 *                          subface < cell->face(f)->n_children(); 
 *                          ++subface) 
 *                       if (cell_is_in_fluid_domain( 
 *                             cell->neighbor_child_on_subface(f, subface))) 
 *                         { 
 *                           elasticity_fe_subface_values.reinit(cell, f, subface); 
 *                           stokes_fe_face_values.reinit( 
 *                             cell->neighbor_child_on_subface(f, subface), 
 *                             cell->neighbor_of_neighbor(f)); 
 * 
 *                           assemble_interface_term(elasticity_fe_subface_values, 
 *                                                   stokes_fe_face_values, 
 *                                                   elasticity_phi, 
 *                                                   stokes_symgrad_phi_u, 
 *                                                   stokes_phi_p, 
 *                                                   local_interface_matrix); 
 * 
 *                           cell->neighbor_child_on_subface(f, subface) 
 *                             ->get_dof_indices(neighbor_dof_indices); 
 *                           constraints.distribute_local_to_global( 
 *                             local_interface_matrix, 
 *                             local_dof_indices, 
 *                             neighbor_dof_indices, 
 *                             system_matrix); 
 *                         } 
 *                   } 
 * 
 * @endcode
 * 
 * 最后一个选项是，邻居比较粗大。在这种情况下，我们必须为邻居使用一个FESubfaceValues对象，为当前单元使用一个FEFaceValues；其余部分与之前相同。
 * 

 * 
 * 
 * @code
 *                 else if (cell->neighbor_is_coarser(f) && 
 *                          cell_is_in_fluid_domain(cell->neighbor(f))) 
 *                   { 
 *                     elasticity_fe_face_values.reinit(cell, f); 
 *                     stokes_fe_subface_values.reinit( 
 *                       cell->neighbor(f), 
 *                       cell->neighbor_of_coarser_neighbor(f).first, 
 *                       cell->neighbor_of_coarser_neighbor(f).second); 
 * 
 *                     assemble_interface_term(elasticity_fe_face_values, 
 *                                             stokes_fe_subface_values, 
 *                                             elasticity_phi, 
 *                                             stokes_symgrad_phi_u, 
 *                                             stokes_phi_p, 
 *                                             local_interface_matrix); 
 * 
 *                     cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
 *                     constraints.distribute_local_to_global( 
 *                       local_interface_matrix, 
 *                       local_dof_indices, 
 *                       neighbor_dof_indices, 
 *                       system_matrix); 
 *                   } 
 *               } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 在组装全局系统的函数中，我们将计算接口条款传递给我们在此讨论的一个单独的函数。关键是，尽管我们无法预测FEFaceValues和FESubfaceValues对象的组合，但它们都是从FEFaceValuesBase类派生出来的，因此我们不必在意：该函数被简单地调用，有两个这样的对象表示面的两边的正交点上的形状函数值。然后我们做我们一直在做的事情：我们用形状函数的值和它们的导数来填充从头数组，然后循环计算矩阵的所有条目来计算局部积分。我们在这里评估的双线性形式的细节在介绍中给出。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::assemble_interface_term( 
 *     const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
 *     const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
 *     std::vector<Tensor<1, dim>> &         elasticity_phi, 
 *     std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
 *     std::vector<double> &                 stokes_phi_p, 
 *     FullMatrix<double> &                  local_interface_matrix) const 
 *   { 
 *     Assert(stokes_fe_face_values.n_quadrature_points == 
 *              elasticity_fe_face_values.n_quadrature_points, 
 *            ExcInternalError()); 
 *     const unsigned int n_face_quadrature_points = 
 *       elasticity_fe_face_values.n_quadrature_points; 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 *     const FEValuesExtractors::Vector displacements(dim + 1); 
 * 
 *     local_interface_matrix = 0; 
 *     for (unsigned int q = 0; q < n_face_quadrature_points; ++q) 
 *       { 
 *         const Tensor<1, dim> normal_vector = 
 *           elasticity_fe_face_values.normal_vector(q); 
 * 
 *         for (unsigned int k = 0; k < stokes_fe_face_values.dofs_per_cell; ++k) 
 *           { 
 *             stokes_symgrad_phi_u[k] = 
 *               stokes_fe_face_values[velocities].symmetric_gradient(k, q); 
 *             stokes_phi_p[k] = stokes_fe_face_values[pressure].value(k, q); 
 *           } 
 *         for (unsigned int k = 0; k < elasticity_fe_face_values.dofs_per_cell; 
 *              ++k) 
 *           elasticity_phi[k] = 
 *             elasticity_fe_face_values[displacements].value(k, q); 
 * 
 *         for (unsigned int i = 0; i < elasticity_fe_face_values.dofs_per_cell; 
 *              ++i) 
 *           for (unsigned int j = 0; j < stokes_fe_face_values.dofs_per_cell; ++j) 
 *             local_interface_matrix(i, j) += 
 *               -((2 * viscosity * (stokes_symgrad_phi_u[j] * normal_vector) - 
 *                  stokes_phi_p[j] * normal_vector) * 
 *                 elasticity_phi[i] * stokes_fe_face_values.JxW(q)); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemsolvecode"></a> 
 * <h4><code>FluidStructureProblem::solve</code></h4>
 * 

 * 
 * 正如介绍中所讨论的，我们在这里使用了一个相当琐碎的求解器：我们只是将线性系统传递给SparseDirectUMFPACK直接求解器（例如，见 step-29  ）。在求解之后，我们唯一要做的是确保悬挂的节点和边界值约束是正确的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::solve() 
 *   { 
 *     SparseDirectUMFPACK direct_solver; 
 *     direct_solver.initialize(system_matrix); 
 *     direct_solver.vmult(solution, system_rhs); 
 * 
 *     constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemoutput_resultscode"></a> 
 * <h4><code>FluidStructureProblem::output_results</code></h4>
 * 

 * 
 * 生成图形输出在这里相当简单：我们所要做的就是确定解向量的哪些成分属于标量和/或向量（例如，见 step-22 之前的例子），然后把它全部传递给DataOut类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::output_results( 
 *     const unsigned int refinement_cycle) const 
 *   { 
 *     std::vector<std::string> solution_names(dim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       solution_names.emplace_back("displacement"); 
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       data_component_interpretation.push_back( 
 *         DataComponentInterpretation::component_is_part_of_vector); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 *     data_out.add_data_vector(solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output( 
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemrefine_meshcode"></a> 
 * <h4><code>FluidStructureProblem::refine_mesh</code></h4>
 * 

 * 
 * 下一步是细化网格。正如在介绍中所讨论的，这有点棘手，主要是因为流体和固体子域使用的变量具有不同的物理尺寸，因此，误差估计的绝对大小不能直接比较。因此，我们将不得不对它们进行缩放。因此，在函数的顶部，我们首先分别计算不同变量的误差估计值（在流体域中使用速度而不是压力，在固体域中使用位移）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::refine_mesh() 
 *   { 
 *     Vector<float> stokes_estimated_error_per_cell( 
 *       triangulation.n_active_cells()); 
 *     Vector<float> elasticity_estimated_error_per_cell( 
 *       triangulation.n_active_cells()); 
 * 
 *     const QGauss<dim - 1> stokes_face_quadrature(stokes_degree + 2); 
 *     const QGauss<dim - 1> elasticity_face_quadrature(elasticity_degree + 2); 
 * 
 *     hp::QCollection<dim - 1> face_q_collection; 
 *     face_q_collection.push_back(stokes_face_quadrature); 
 *     face_q_collection.push_back(elasticity_face_quadrature); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       face_q_collection, 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       stokes_estimated_error_per_cell, 
 *       fe_collection.component_mask(velocities)); 
 * 
 *     const FEValuesExtractors::Vector displacements(dim + 1); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       face_q_collection, 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       elasticity_estimated_error_per_cell, 
 *       fe_collection.component_mask(displacements)); 
 * 
 * @endcode
 * 
 * 然后，我们通过除以误差估计值的法线对其进行归一化处理，并按照介绍中所讨论的那样，将流体误差指标按4的系数进行缩放。然后将这些结果加在一起，形成一个包含所有单元的误差指标的向量。
 * 

 * 
 * 
 * @code
 *     stokes_estimated_error_per_cell *= 
 *       4. / stokes_estimated_error_per_cell.l2_norm(); 
 *     elasticity_estimated_error_per_cell *= 
 *       1. / elasticity_estimated_error_per_cell.l2_norm(); 
 * 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     estimated_error_per_cell += stokes_estimated_error_per_cell; 
 *     estimated_error_per_cell += elasticity_estimated_error_per_cell; 
 * 
 * @endcode
 * 
 * 在实际细化网格之前，函数的倒数第二部分涉及到我们在介绍中已经提到的启发式方法：由于解是不连续的，KellyErrorEstimator类对位于子域之间边界的单元感到困惑：它认为那里的误差很大，因为梯度的跳跃很大，尽管这完全是预期的，事实上在精确解中也存在这一特征，因此不表明任何数值错误。
 * 

 * 
 * 因此，我们将界面上的所有单元的误差指标设置为零；决定影响哪些单元的条件略显尴尬，因为我们必须考虑到自适应细化网格的可能性，这意味着邻近的单元可能比当前的单元更粗，或者事实上可能被细化一些。这些嵌套条件的结构与我们在 <code>assemble_system</code> 中组装接口条款时遇到的情况基本相同。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       for (const auto f : cell->face_indices()) 
 *         if (cell_is_in_solid_domain(cell)) 
 *           { 
 *             if ((cell->at_boundary(f) == false) && 
 *                 (((cell->neighbor(f)->level() == cell->level()) && 
 *                   (cell->neighbor(f)->has_children() == false) && 
 *                   cell_is_in_fluid_domain(cell->neighbor(f))) || 
 *                  ((cell->neighbor(f)->level() == cell->level()) && 
 *                   (cell->neighbor(f)->has_children() == true) && 
 *                   (cell_is_in_fluid_domain( 
 *                     cell->neighbor_child_on_subface(f, 0)))) || 
 *                  (cell->neighbor_is_coarser(f) && 
 *                   cell_is_in_fluid_domain(cell->neighbor(f))))) 
 *               estimated_error_per_cell(cell->active_cell_index()) = 0; 
 *           } 
 *         else 
 *           { 
 *             if ((cell->at_boundary(f) == false) && 
 *                 (((cell->neighbor(f)->level() == cell->level()) && 
 *                   (cell->neighbor(f)->has_children() == false) && 
 *                   cell_is_in_solid_domain(cell->neighbor(f))) || 
 *                  ((cell->neighbor(f)->level() == cell->level()) && 
 *                   (cell->neighbor(f)->has_children() == true) && 
 *                   (cell_is_in_solid_domain( 
 *                     cell->neighbor_child_on_subface(f, 0)))) || 
 *                  (cell->neighbor_is_coarser(f) && 
 *                   cell_is_in_solid_domain(cell->neighbor(f))))) 
 *               estimated_error_per_cell(cell->active_cell_index()) = 0; 
 *           } 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.0); 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeFluidStructureProblemruncode"></a> 
 * <h4><code>FluidStructureProblem::run</code></h4>
 * 

 * 
 * 像往常一样，这是控制整个操作流程的函数。如果你读过教程程序  step-1  到  step-6  ，例如，那么你已经对以下结构相当熟悉。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void FluidStructureProblem<dim>::run() 
 *   { 
 *     make_grid(); 
 * 
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 10 - 2 * dim; 
 *          ++refinement_cycle) 
 *       { 
 *         std::cout << "Refinement cycle " << refinement_cycle << std::endl; 
 * 
 *         if (refinement_cycle > 0) 
 *           refine_mesh(); 
 * 
 *         setup_dofs(); 
 * 
 *         std::cout << "   Assembling..." << std::endl; 
 *         assemble_system(); 
 * 
 *         std::cout << "   Solving..." << std::endl; 
 *         solve(); 
 * 
 *         std::cout << "   Writing output..." << std::endl; 
 *         output_results(refinement_cycle); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step46 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h4>The <code>main()</code> function</h4>
 * 

 * 
 * 这个，最后的，函数所包含的内容几乎与其他大多数教程程序的内容完全一样。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step46; 
 * 
 *       FluidStructureProblem<2> flow_problem(1, 1); 
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
examples/step-46/doc/results.dox

<a name="Results"></a>

<a name="Results"></a><h1>Results</h1>


<a name="2dresults"></a><h3>2d results</h3>



当运行该程序时，你应该得到如下输出。

@code
Refinement cycle 0
   Number of active cells: 64
   Number of degrees of freedom: 531
   Assembling...
   Solving...
   Writing output...


Refinement cycle 1
   Number of active cells: 136
   Number of degrees of freedom: 1260
   Assembling...
   Solving...
   Writing output...


Refinement cycle 2
   Number of active cells: 436
   Number of degrees of freedom: 3723
   Assembling...
   Solving...
   Writing output...


Refinement cycle 3
   Number of active cells: 1072
   Number of degrees of freedom: 7493
   Assembling...
   Solving...
   Writing output...


Refinement cycle 4
   Number of active cells: 2632
   Number of degrees of freedom: 15005
   Assembling...
   Solving...
   Writing output...


Refinement cycle 5
   Number of active cells: 5944
   Number of degrees of freedom: 29437
   Assembling...
   Solving...
   Writing output...
@endcode



结果很容易被可视化。

 <table width="80%" align="center">
  <tr valign="top">
    <td valign="top" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-46.9.2.velocity-2d.png" alt="">
      <p align="center">
        Magnitude and vectors for the fluid velocity.
      </p>
    </td>
    <td valign="top" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-46.9.2.pressure-2d.png" alt="">
      <p align="center">
        Fluid pressure. The dynamic range has been truncated to cut off the
        pressure singularities at the top left and right corners of the domain
        as well as the top corners of the solid that forms re-entrant corners
        into the fluid domain.
      </p>
    </td>
    <td valign="top" align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-46.9.2.displacement-2d.png" alt="">
      <p align="center">
        Magnitude and vectors for the solid displacement.
      </p>
    </td>
  </tr>
</table> 

这些图很容易解释：当水流在固体直立部分的左边向下、右边向上时，它产生的压力在左边高，在右边低，这些力量使固体的垂直部分向右弯曲。




<a name="3dresults"></a><h3>3d results</h3>


通过将 <code>FluidStructureProblem</code> 类中的维度改为3，我们也可以运行同样的问题3D。你会得到如下的输出。

@code
Refinement cycle 0
   Number of active cells: 512
   Number of degrees of freedom: 11631
   Assembling...
   Solving...
   Writing output...


Refinement cycle 1
   Number of active cells: 1716
   Number of degrees of freedom: 48984
   Assembling...
   Solving...
   Writing output...


Refinement cycle 2
   Number of active cells: 8548
   Number of degrees of freedom: 245746
   Assembling...
   Solving...
@endcode

你会发现，最大的瓶颈是求解器。SparseDirectUmfpack在2016年的工作站上解决这个问题的最后一次迭代需要将近5个小时和大约80GB的内存（倒数第二次迭代只用了16分钟）。显然，这里需要一个更好的求解器，这个话题在下面讨论。

结果也可以被可视化，并产生良好的图片。这里有一张，显示了速度的矢量图（橙色），实体位移（蓝色），以及实体区域的阴影。

<p align="center">  <img src="https://www.dealii.org/images/steps/developer/step-46.9.2.3d.png" alt="">   </p>  。

除了缺乏一个好的求解器之外，网格也有点不平衡：网格细化严重偏向于流体子域（在2D中，情况恰恰相反，促使我们对流体误差指标的权重更高）。显然，如果想继续做更多的三维计算，对两个子域中误差指标的相对重要性进行一些调整是很重要的。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Linearsolversandpreconditioners"></a><h4>Linear solvers and preconditioners</h4>


改进程序的一个明显的地方是使用一个更复杂的求解器&mdash；特别是一个能很好地扩展并能解决现实的三维问题的求解器。这在这里其实不难实现，因为从流体到固体的单向耦合。为此，假设我们将自由度重新排序，首先是所有的速度和压力自由度，然后是所有的位移（使用 DoFRenumbering::component_wise). 很容易实现），那么系统矩阵可以被分成以下的块状形式。

@f[
  A_\text{global}
  =
  \begin{pmatrix}
    A_{\text{fluid}} & 0 \\
    B & A_{\text{solid}}
  \end{pmatrix}


@f]

其中 $A_{\text{fluid}}$ 是速度和压力的斯托克斯矩阵（它可以进一步细分为 $2\times 2$ 矩阵，如步骤22，尽管这对目前的目的不重要）， $A_{\text{solid}}$ 是位移的弹性方程的结果， $B$ 是来自界面条件的矩阵。现在注意到，这个矩阵

@f[
  A_\text{global}^{-1}
  =
  \begin{pmatrix}
    A_{\text{fluid}}^{-1} & 0 \\


    -A_\text{solid}^{-1} B
      A_\text{fluid}^{-1} & A_{\text{solid}}^{-1}
  \end{pmatrix}


@f]

是 $A_\text{global}$ 的逆数。应用这个矩阵只需要与 $A_\text{fluid}$ 和 $A_\text{solid}$ 各解一次，因为

@f[
  \begin{pmatrix}
    p_x \\ p_y
  \end{pmatrix}
  =
  \begin{pmatrix}
    A_{\text{fluid}}^{-1} & 0 \\


    -A_\text{solid}^{-1} B
      A_\text{fluid}^{-1} & A_{\text{solid}}^{-1}
  \end{pmatrix}
  \begin{pmatrix}
    x \\ y
  \end{pmatrix}


@f]

可以计算为 $p_x = A_{\text{fluid}}^{-1} x$ ，然后是 $p_y = A_{\text{solid}}^{-1} (y-Bp_x)$  。

因此，人们可以预期，

@f[
  \widetilde{A_\text{global}^{-1}}
  =
  \begin{pmatrix}
    \widetilde{A_{\text{fluid}}^{-1}} & 0 \\


    -\widetilde{A_\text{solid}^{-1}} B
      \widetilde{A_\text{fluid}^{-1}} & \widetilde{A_{\text{solid}}^{-1}}
  \end{pmatrix}


@f]

如果  $\widetilde{A_{\text{fluid}}^{-1}}
\approx A_{\text{fluid}}^{-1}, \widetilde{A_{\text{solid}}^{-1}}
\approx A_{\text{solid}}^{-1}$  ，将是一个好的预处理程序。

这意味着，我们只需要为斯托克斯和弹性方程分别提供良好的预处理。这些都是众所周知的：对于斯托克斯，我们可以使用step-22的结果部分所讨论的预处理程序；对于弹性，一个好的预处理程序将是一个几何或代数多重网格的单一V形周期。然而，还有更多的开放性问题。对于由两个子调节器构建的 "优化 "求解器块-三角调节器，经常出现的一点是，在为子调节器选择参数时，在单独解决两个问题时效果好的值，在组合成多物理学调节器时可能不是最佳值。  特别是，当单独解决固体或流体力学问题时，在收敛所需的迭代次数和每次迭代应用预调节器的成本之间进行平衡，可能导致人们为斯托克斯问题选择昂贵的预调节器，为弹性问题选择便宜的预调节器（或者反之）。  然而，当结合在一起时，还有一个额外的约束，即你希望两个子预处理程序以大致相同的速度收敛，否则便宜的预处理程序可能会增加全局的迭代次数，而昂贵的预处理程序则会增加每迭代的成本。例如，虽然单一的AMG V型循环本身是一个很好的弹性方法，但当结合到一个多物理问题时，可能会有动力使用一个完整的W型循环或多个循环来帮助降低总求解时间。




<a name="Refinementindicators"></a><h4>Refinement indicators</h4>


正如介绍中提到的，我们在这个程序中使用的细化指标是相当临时的。一个更好的会明白，解的梯度在界面上的跳跃并不是错误的指示，而是可以预期的，并且在积分跳跃项的时候忽略界面。然而，这并不是KellyErrorEstimator类所做的。另一个更大的问题是，这种估计器首先是否是一个好的策略：例如，如果我们想在位移的某个特定方面（例如实体右上角的位移）有最大的准确性，那么将流体和实体的误差指标放大到相同的程度是否合适？也许有必要用比固体更高的精度来解决流体问题，因为流体的解决方案直接影响到固体的解决方案？也许正好相反？

因此，改进该程序的一个明显的可能性是实施一个更好的细化标准。关于这个话题有一些文献；各种可能的起点之一是Thomas Wick的论文 "Adaptive finite elements for monolithic fluid-structure interaction on a prolongated domain:应用于心脏瓣膜模拟"，2011年机械学计算机方法会议论文集（CMM-2011），2011年5月9-12日，波兰华沙。




<a name="Verification"></a><h4>Verification</h4>


上面的结果纯粹是定性的，因为没有证据表明我们的方案实际上是收敛的。因此，一个显而易见的事情是增加一些量化的措施来检查该方案至少收敛到<i>something</i>。例如，我们可以为每个细化周期输出实体的右上角突出到流体子域的部分的偏移。或者我们可以计算出流体对实体施加的净力向量或扭矩。




<a name="Bettermodels"></a><h4>Better models</h4>


在现实中，大多数流体结构的相互作用问题是这样的：固体的运动确实影响了流体的流动。例如，空气在气箔周围的作用力导致它弯曲并改变其形状。同样地，一面旗帜在风中飘动，完全改变了它的形状。

这种双向耦合的问题通常在任意拉格朗日欧拉（ALE）框架下处理，其中固体的位移以某种平滑的方式扩展到流体领域，而不是像我们在这里所做的那样以零为单位。然后，扩展的位移场被用来对我们计算流体流动的网格进行变形。此外，界面上流体的边界条件不再是速度为零；相反，在一个时间相关的程序中，流体速度必须等于沿界面的位移的时间导数。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-46.cc"
*/
