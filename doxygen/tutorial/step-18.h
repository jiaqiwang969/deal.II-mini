/**
@page step_18 The step-18 tutorial program
This tutorial depends on step-17.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Quasistaticelasticdeformation">Quasistatic elastic deformation</a>
      <ul>
        <li><a href="#Motivationofthemodel">Motivation of the model</a>
        <li><a href="#Timediscretization">Time discretization</a>
        <li><a href="#Updatingthestressvariable">Updating the stress variable</a>
      </ul>
        <li><a href="#Parallelgraphicaloutput">Parallel graphical output</a>
        <li><a href="#Atriangulationwithautomaticpartitioning">A triangulation with automatic partitioning</a>
        <li><a href="#Overallstructureoftheprogram">Overall structure of the program</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ThecodePointHistorycodeclass">The <code>PointHistory</code> class</a>
        <li><a href="#Thestressstraintensor">The stress-strain tensor</a>
        <li><a href="#Auxiliaryfunctions">Auxiliary functions</a>
        <li><a href="#ThecodeTopLevelcodeclass">The <code>TopLevel</code> class</a>
        <li><a href="#ThecodeBodyForcecodeclass">The <code>BodyForce</code> class</a>
        <li><a href="#ThecodeIncrementalBoundaryValuecodeclass">The <code>IncrementalBoundaryValue</code> class</a>
        <li><a href="#ImplementationofthecodeTopLevelcodeclass">Implementation of the <code>TopLevel</code> class</a>
      <ul>
        <li><a href="#Thepublicinterface">The public interface</a>
        <li><a href="#TopLevelcreate_coarse_grid">TopLevel::create_coarse_grid</a>
        <li><a href="#TopLevelsetup_system">TopLevel::setup_system</a>
        <li><a href="#TopLevelassemble_system">TopLevel::assemble_system</a>
        <li><a href="#TopLevelsolve_timestep">TopLevel::solve_timestep</a>
        <li><a href="#TopLevelsolve_linear_problem">TopLevel::solve_linear_problem</a>
        <li><a href="#TopLeveloutput_results">TopLevel::output_results</a>
        <li><a href="#TopLeveldo_initial_timestep">TopLevel::do_initial_timestep</a>
        <li><a href="#TopLeveldo_timestep">TopLevel::do_timestep</a>
        <li><a href="#TopLevelrefine_initial_grid">TopLevel::refine_initial_grid</a>
        <li><a href="#TopLevelmove_mesh">TopLevel::move_mesh</a>
        <li><a href="#TopLevelsetup_quadrature_point_history">TopLevel::setup_quadrature_point_history</a>
        <li><a href="#TopLevelupdate_quadrature_point_history">TopLevel::update_quadrature_point_history</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
      <ul>
        <li><a href="#Plasticitymodels">Plasticity models</a>
        <li><a href="#Stabilizationissues">Stabilization issues</a>
        <li><a href="#Refinementduringtimesteps">Refinement during timesteps</a>
        <li><a href="#Ensuringmeshregularity">Ensuring mesh regularity</a>
    </ul>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-18/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这个教程程序是我们已经在步骤8和步骤17中开始的弹性问题系列中的另一个。它将其扩展到两个不同的方向：首先，它用拉格朗日网格运动方法解决了大变形的准静态但与时间相关的弹性问题。其次，它又展示了一些使用PETSc的线性代数的%并行处理来解决此类问题的技术。除此之外，我们还展示了如何解决step-17的两个主要瓶颈中的一个，即我们只从一个进程中产生图形输出，而这在更多的进程和大问题上的扩展性非常差。另一个瓶颈，即每个处理器都必须持有整个网格和DoFHandler，将在第40步中解决）。最后，我们还展示了许多以前的程序中未曾展示过的各种改进和技术。

如同前面的第17步，只要你安装了PETSc，程序在单机上的运行也是一样的。关于如何告诉deal.II你的系统上安装了PETSc的信息可以在deal.II的README文件中找到，该文件可以从你安装的deal.II的<a href="../../index.html">main
documentation page</a>中链接到，或者在<a href="http://www.dealii.org/">the
deal.II webpage</a>上。




<a name="Quasistaticelasticdeformation"></a><h3>Quasistatic elastic deformation</h3>


<a name="Motivationofthemodel"></a><h4>Motivation of the model</h4>


一般来说，随时间变化的小弹性变形是由弹性波方程描述的

@f[
  \rho \frac{\partial^2 \mathbf{u}}{\partial t^2}
  + c \frac{\partial \mathbf{u}}{\partial t}


  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) = \mathbf{f}
  \qquad
  \textrm{in}\ \Omega,


@f]

其中 $\mathbf{u}=\mathbf{u} (\mathbf{x},t)$ 是身体的变形， $\rho$ 和 $c$ 是密度和衰减系数，以及 $\mathbf{f}$ 外力。此外，初始条件

@f[
  \mathbf{u}(\cdot, 0) = \mathbf{u}_0(\cdot)
  \qquad
  \textrm{on}\ \Omega,


@f]

和Dirichlet（位移）或Neumann（牵引）边界条件，需要指定一个唯一的解决方案。

@f{eqnarray*}
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega\backslash\Gamma_D.


@f}

在上述公式中， $\varepsilon(\mathbf{u})= \frac 12 (\nabla \mathbf{u} + \nabla
\mathbf{u}^T)$ 是位移的对称梯度，也称为 <em> 应变 </em>  。   $C$ 是一个等级为4的张量，称为 <em> 应力-应变张量 </em> （<a
  href="https://en.wikipedia.org/wiki/Hooke%27s_law#Hooke's_law_for_continuous_media"><em>compliance
  tensor</em></a>的逆向）。]），它包含了材料弹性强度的知识；它的对称性特性确保它将秩为2的对称张量（&ldquo;矩阵&rdquo;的维数 $d$ ，其中 $d$ 是空间维数）映射到相同秩的对称张量上。我们将在下面更多地评论应变和应力张量的作用。现在只需要说，我们将术语 $\textrm{div}\  ( C \varepsilon(\mathbf{u}))$ 解释为具有分量 $\frac \partial{\partial x_j} C_{ijkl} \varepsilon(\mathbf{u})_{kl}$ 的向量，其中对指数 $j,k,l$ 的求和是隐含的。

这个方程的准静态极限的动机如下：身体的每个小扰动，例如边界条件或强迫函数的变化，将导致身体配置的相应变化。一般来说，这将是以波的形式从扰动的位置辐射出去。由于阻尼项的存在，这些波将在例如 $\tau$ 的时间尺度上被衰减。现在，假设所有外部强制力的变化发生在比 $\tau$ 大得多的时间尺度上。在这种情况下，变化的动态性质并不重要：我们可以认为身体总是处于静态平衡状态，也就是说，我们可以假设在任何时候，身体都满足于

@f{eqnarray*}


  - \textrm{div}\  ( C \varepsilon(\mathbf{u})) &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega,
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D,
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N.


@f}

请注意，微分方程不再包含任何时间导数 -- 所有的时间依赖性都是通过边界条件和可能的时间变化的力函数引入的  $\mathbf{f}(\mathbf{x},t)$  。因此，配置的变化可以被认为是瞬时静止的。对此的另一种看法是， $t$ 并不是真正的时间变量，而只是一个支配问题演变的类似时间的参数。

虽然这些方程足以描述小的变形，但计算大的变形就有点复杂了，一般来说，会导致非线性方程，如步骤-44中处理的那些。在下文中，让我们考虑在模拟变形成为<i>large</i>的问题时，人们会采用的一些工具。

 @note 我们下面要考虑的模型并不是建立在任何在数学上合理的基础上：我们将考虑一个模型，在这个模型中，我们产生一个小的变形，通过这个变形使身体的物理坐标变形，然后再考虑下一个加载步骤，作为一个线性问题。这并不一致，因为线性的假设意味着变形是无限小的，所以在解决下一个线性问题之前，在我们的网格顶点周围移动一个有限的量是不一致的做法。因此，我们应该注意到，在文献中找不到下面讨论的方程，这并不奇怪。<b>The model considered here has
little to do with reality!</b>另一方面，我们所考虑的实现技术正是人们在实现<i>real</i>模型时需要使用的，我们将在步骤-44中看到。


为了回到定义我们的 "人工 "模型，让我们首先引入一个张量的应力变量 $\sigma$ ，并以应力为基础写出微分方程。

@f{eqnarray*}


  - \textrm{div}\  \sigma &=& \mathbf{f}(\mathbf{x},t)
  \qquad
  \textrm{in}\ \Omega(t),
  \\
  \mathbf{u}(\mathbf{x},t) &=& \mathbf{d}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_D\subset\partial\Omega(t),
  \\
  \mathbf{n} \ C \varepsilon(\mathbf{u}(\mathbf{x},t)) &=& \mathbf{b}(\mathbf{x},t)
  \qquad
  \textrm{on}\ \Gamma_N=\partial\Omega(t)\backslash\Gamma_D.


@f}

注意这些方程是在一个随时间变化的域 $\Omega(t)$ 上提出的，边界根据边界上各点的位移 $\mathbf{u}(\mathbf{x},t)$ 而移动。为了完成这个系统，我们必须指定应力和应变之间的增量关系，如下所示。<a name="step_18.stress-strain"></a>

@f[
  \dot\sigma = C \varepsilon (\dot{\mathbf{u}}),
  \qquad
  \qquad
  \textrm{[stress-strain]}


@f]

其中点表示一个时间导数。应力 $\sigma$ 和应变 $\varepsilon(\mathbf{u})$ 都是等级2的对称张量。




<a name="Timediscretization"></a><h4>Time discretization</h4>


在数值上，该系统的求解方法如下：首先，我们使用后向欧拉方案对时间部分进行离散化。这导致了时间步长的离散平衡力  $n$  。

@f[


  -\textrm{div}\  \sigma^n = f^n,


@f]

其中

@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n),


@f]

和  $\Delta \mathbf{u}^n$  时间步长的增量位移  $n$  。此外，我们必须指定初始数据  $\mathbf{u}(\cdot,0)=\mathbf{u}_0$  。这样一来，如果我们想求解位移增量，我们必须求解以下系统。

@f{align*}


  - \textrm{div}\   C \varepsilon(\Delta\mathbf{u}^n) &= \mathbf{f} + \textrm{div}\  \sigma^{n-1}
  \qquad
  &&\textrm{in}\ \Omega(t_{n-1}),
  \\
  \Delta \mathbf{u}^n(\mathbf{x},t) &= \mathbf{d}(\mathbf{x},t_n) - \mathbf{d}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_D\subset\partial\Omega(t_{n-1}),
  \\
  \mathbf{n} \ C \varepsilon(\Delta \mathbf{u}^n(\mathbf{x},t)) &= \mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1})
  \qquad
  &&\textrm{on}\ \Gamma_N=\partial\Omega(t_{n-1})\backslash\Gamma_D.


@f}

这组方程的弱形式，像往常一样是有限元公式的基础，其内容如下：找到 $\Delta \mathbf{u}^n \in
\{v\in H^1(\Omega(t_{n-1}))^d: v|_{\Gamma_D}=\mathbf{d}(\cdot,t_n) - \mathbf{d}(\cdot,t_{n-1})\}$ ，使<a name="step_18.linear-system"></a>这样的方程。

@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}


  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  \\
  &\qquad
  +(\mathbf{b}(\mathbf{x},t_n)-\mathbf{b}(\mathbf{x},t_{n-1}), \varphi)_{\Gamma_N}
  +(\sigma^{n-1} \mathbf{n}, \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.


@f}

利用 $\sigma^{n-1} \mathbf{n}
            = [C \varepsilon(\mathbf{u}^{n-1})] \mathbf{n}
            = \mathbf{b}(\mathbf x, t_{n-1})$ ，这些方程可以简化为

@f{align*}
  (C \varepsilon(\Delta\mathbf{u}^n), \varepsilon(\varphi) )_{\Omega(t_{n-1})}
  &=
  (\mathbf{f}, \varphi)_{\Omega(t_{n-1})}


  -(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  +(\mathbf{b}(\mathbf{x},t_n),t_{n-1}), \varphi)_{\Gamma_N}
  \\
  &\qquad\qquad
  \forall \varphi \in \{\mathbf{v}\in H^1(\Omega(t_{n-1}))^d: \mathbf{v}|_{\Gamma_D}=0\}.
  \qquad
  \qquad
  \textrm{[linear-system]}


@f}



我们注意到，为了简单起见，在程序中我们总是假设没有边界力，即 $\mathbf{b} = 0$ ，并且身体的变形仅由身体力 $\mathbf{f}$ 和规定的边界位移 $\mathbf{d}$ 驱动。还值得注意的是，当通过部分积分时，我们会得到形式为 $(C \varepsilon(\Delta\mathbf{u}^n), \nabla \varphi
)_{\Omega(t_{n-1})}$ 的条款，但我们用涉及对称梯度的条款 $\varepsilon(\varphi)$ 而不是 $\nabla\varphi$ 来取代它们。由于 $C$ 的对称性，这两个项在数学上是等价的，但对称版本避免了可能出现的四舍五入错误，使得到的矩阵略显非对称性。

在时间步长 $n$ 的系统，要在旧域 $\Omega(t_{n-1})$ 上求解，其形式完全是一个静止的弹性问题，因此与我们在以前的例子程序中已经实现的类似。因此，除了说我们再次使用最低阶连续有限元之外，我们将不对空间离散化进行评论。

但也有不同之处。<ol>  <li>  我们必须在每个时间步骤之后移动（更新）网格，以便能够在新的领域上解决下一个时间步骤。

    <li> 我们需要知道 $\sigma^{n-1}$ 来计算下一个增量位移，也就是说，我们需要在时间步骤结束时计算它，以确保它可以用于下一个时间步骤。从本质上讲，应力变量是我们了解体的变形历史的窗口。   </ol>  这两个操作在程序中的 <code>move_mesh</code> 和 <code>update_quadrature_point_history</code> 函数中完成。移动网格只是一个技术问题，而更新应力则要复杂一些，将在下一节讨论。




<a name="Updatingthestressvariable"></a><h4>Updating the stress variable</h4>


如上所述，在计算时间步长 $n+1$ 时，我们需要有应力变量 $\sigma^n$ ，我们可以用<a name="step_18.stress-update"></a>来计算它。

@f[
  \sigma^n = \sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update]}


@f]

尽管这个方程表面上很简单，但有两个问题我们需要讨论。第一个问题是关于我们存储 $\sigma^n$ 的方式：即使我们使用最低阶有限元计算增量更新 $\Delta\mathbf{u}^n$ ，那么其对称梯度 $\varepsilon(\Delta\mathbf{u}^n)$ 一般来说仍然是一个不容易描述的函数。特别是，它不是一个片状常数函数，在一般的网格上（单元不是平行于坐标轴的矩形）或非恒定应力-应变张量 $C$ ，它甚至不是一个双线性或三线性函数。因此，如何在计算机程序中存储 $\sigma^n$ 是先验的。

要决定这一点，我们必须看它被用在什么地方。我们需要应力的唯一地方是在术语 $(\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}$ 中。在实践中，我们当然会用数值正交来代替这个项。

@f[
  (\sigma^{n-1},\varepsilon(\varphi))_{\Omega(t_{n-1})}
  =
  \sum_{K\subset {T}}
  (\sigma^{n-1},\varepsilon(\varphi))_K
  \approx
  \sum_{K\subset {T}}
  \sum_q
  w_q \ \sigma^{n-1}(\mathbf{x}_q) : \varepsilon(\varphi(\mathbf{x}_q),


@f]

其中 $w_q$ 是正交权重， $\mathbf{x}_q$ 是单元格 $K$ 上的正交点。这应该表明，我们真正需要的不是应力 $\sigma^{n-1}$ 本身，而只是所有单元上的正交点的应力值。然而，这是一个更简单的任务：我们只需要提供一个数据结构，能够为所有单元上的每个正交点（或者，由于我们是并行计算，目前的MPI进程拥有的所有单元的所有正交点&ldquo;rdquo;）容纳一个等级为2的对称张量。在每个时间步骤结束时，我们只需评估 $\varepsilon(\Delta \mathbf{u}^n(\mathbf{x}_q))$ ，将其乘以应力-应变张量 $C$ ，并使用该结果来更新正交点 $q$ 的应力 $\sigma^n(\mathbf{x}_q)$ 。

第二个复杂的问题在我们上面选择的符号中并不明显。这是由于我们在域 $\Omega(t_{n-1})$ 上计算 $\Delta u^n$ ，然后用这个位移增量来更新应力，同时移动网格节点，以达到 $\Omega(t_n)$ ，在此基础上计算下一个增量。在这种情况下，我们必须确定的是，移动网格不仅涉及到节点的移动，还涉及到应力变量的相应变化：更新的应力是一个相对于旧域中材料的坐标系而定义的变量，必须转移到新域中。其原因可以理解为：在局部，增量变形 $\Delta\mathbf{u}$ 可以分解为三个部分，线性平移（点附近的位移增量场的常数部分），扩张分量（位移场梯度中具有非零发散的那部分），以及旋转。材料的线性平移并不影响冻结在其中的应力--应力值只是沿着平移。扩张或压缩的变化产生相应的应力更新。然而，旋转分量不一定会引起非零的应力更新（想想，在2d中，例如 $\Delta\mathbf{u}=(y, -x)^T$  ，与 $\varepsilon(\Delta
\mathbf{u})=0$  的情况）。尽管如此，如果材料在某个方向上被预应力，那么这个方向将随着材料的旋转而旋转。  为此，我们必须定义一个旋转矩阵 $R(\Delta \mathbf{u}^n)$ ，描述在每一个点上由于位移增量而产生的旋转。不难看出， $R$ 对 $\Delta \mathbf{u}^n$ 的实际依赖只能是通过位移的卷曲，而不是位移本身或其全部梯度（如上所述，增量的常数分量描述平移，其发散描述扩张模式，而卷曲描述旋转模式）。由于 $R$ 的确切形式很麻烦，我们只在程序代码中说明，并注意到应力变量的正确更新公式是<a name="step_18.stress-update+rot"></a> 。

@f[
  \sigma^n
  =
  R(\Delta \mathbf{u}^n)^T
  [\sigma^{n-1} + C \varepsilon (\Delta \mathbf{u}^n)]
  R(\Delta \mathbf{u}^n).
  \qquad
  \qquad
  \textrm{[stress-update+rot]}


@f]



应力更新和旋转都是在示例程序的函数 <code>update_quadrature_point_history</code> 中实现的。




<a name="Parallelgraphicaloutput"></a><h3>Parallel graphical output</h3>


在步骤17中，就运行时间而言，平行计算的主要瓶颈是只有第一个处理器产生整个领域的输出。由于生成图形输出是很昂贵的，所以当涉及到更多数量的处理器时，这并不能很好地扩展。我们将在这里解决这个问题。关于程序 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"）。

基本上，我们需要做的是让每个进程为它所拥有的那个单元子集产生图形输出，将它们写进单独的文件，并有办法同时显示某个时间步长的所有文件。这样，代码在每个时间步长的每个进程产生一个 <code>.vtu</code> 文件。两个常见的VTK文件查看器ParaView和Viscit都支持一次打开一个以上的 <code>.vtu</code> 文件。为了简化挑选正确文件的过程，并允许在时间上移动，两者都支持记录文件，以引用特定时间步长的所有文件。遗憾的是，记录文件在VisIt和Paraview之间有不同的格式，所以我们把两种格式都写出来。

代码将生成文件 <code>solution-TTTT.NNN.vtu</code> ，其中 <code>TTTT</code> 是时间步数（从1开始）， <code>NNN</code> 是进程等级（从0开始）。这些文件包含时间段和处理器的本地所有单元。文件 <code>solution-TTTT.visit</code> 是时间段的访问记录 <code>TTTT</code>, while <code>solution-TTTT.pvtu</code> 对ParaView也是如此。(较新版本的VisIt实际上也可以读取 <code>.pvtu</code> 文件，但输出两种记录文件也无妨。)最后， <code>solution.pvd</code> 文件是只有ParaView支持的特殊记录，它引用所有的时间步骤。所以在ParaView中，只需要打开solution.pvd，而在VisIt中需要选择所有的.visit文件组，才能达到同样的效果。




<a name="Atriangulationwithautomaticpartitioning"></a><h3>A triangulation with automatic partitioning</h3>


在第17步中，我们使用了一个在每个处理器上简单复制的常规三角形，以及一个相应的DoFHandler。两者都不知道它们是在%并行环境下使用的--它们只是完整地存在于每个处理器上，我们认为这最终会成为一个主要的内存瓶颈。

我们在这里不解决这个问题（我们将在第40步中解决），但使情况稍微自动化一些。在第17步中，我们创建了三角形，然后手动 "分区"，也就是说，我们给每个单元分配了 @ref GlossSubdomainId "子域ID"，以表明哪个 @ref GlossMPIProcess "MPI进程""拥有 "该单元。在这里，我们使用了一个类 parallel::shared::Triangulation ，它至少自动完成了这一部分：每当你创建或完善这样一个三角图时，它都会自动在所有参与的进程之间进行划分（它知道这些进程，因为你必须告诉它在构建三角图时连接这些进程的 @ref GlossMPICommunicator "MPI通信器"）。否则， parallel::shared::Triangulation 看起来，就所有的实际目的而言，就像一个普通的Triangulation对象。

使用这个类的便利性不仅来自于能够避免手动调用 GridTools::partition(). ，相反，DoFHandler类现在也知道你想在并行环境下使用它，并且默认情况下会自动列举自由度，使进程0拥有的所有DoF先于进程1拥有的所有DoF，等等。换句话说，你也可以避免对 DoFRenumbering::subdomain_wise(). 的调用。

还有其他好处。例如，由于三角计算知道它生活在一个%parallel universe中，它也知道它 "拥有 "某些单元（即那些子域id等于其MPI等级的单元；以前，三角计算只存储这些子域id，但没有办法使它们有意义）。因此，在汇编函数中，你可以测试一个单元是否 "本地拥有"（即由当前进程拥有，见 @ref GlossLocallyOwnedCell ），当你在所有单元上循环时，使用以下语法

@code
  if (cell->is_locally_owned())
@endcode

这种知识延伸到建立在这种三角形上的DoFHandler对象，然后它可以通过 DoFHandler::compute_n_locally_owned_dofs_per_processor() 和 DoFTools::extract_locally_relevant_dofs(). 等调用来识别哪些自由度是本地拥有的（见 @ref GlossLocallyOwnedDof ）。最后，DataOut类也知道如何处理这种三角形，并将简单地跳过在非本地拥有的单元上生成图形输出。

当然，正如在第17步的讨论中多次指出的那样，在每个进程上保持整个三角形将无法扩展：大型问题可能根本无法再适合每个进程的内存，即使我们有足够多的进程在合理的时间内解决它们。在这种情况下， parallel::shared::Triangulation 不再是一个合理的计算基础，我们将在步骤40中展示如何使用 parallel::distributed::Triangulation 类来解决这个问题，即让每个进程只存储一个<i>part</i>的三角图。




<a name="Overallstructureoftheprogram"></a><h3>Overall structure of the program</h3>


程序的整体结构可以从 <code>run()</code> 函数中推断出来，该函数首先在第一个时间步骤中调用 <code>do_initial_timestep()</code> ，然后在所有后续时间步骤中调用 <code>do_timestep()</code> 。这些函数之间的区别仅仅在于，在第一个时间步骤中，我们从一个粗略的网格开始，在其上求解，自适应地细化网格，然后在新的网格上以干净的状态重新开始。这个过程给了我们一个更好的起始网格，尽管我们当然应该在迭代过程中不断调整网格--这个程序中没有这样做，但是下面会有评论。

这两个处理时间步骤的函数的共同部分是在本网格上的以下操作序列。   <ul>   <li>   <code>assemble_system ()</code> [via <code>solve_timestep ()</code>  ] 。   这第一个函数也是最有趣的一个。它组装了对应于方程<a href="#step_18.linear-system">[linear-system]</a>离散化版本的线性系统。这导致了一个系统矩阵 $A_{ij} = \sum_K
  A^K_{ij}$ ，由每个单元 $K$ 上的局部贡献组成，其条目为@f[
    A^K_{ij} = (C \varepsilon(\varphi_j), \varepsilon(\varphi_i))_K;
  @f] 。

  在实践中， $A^K$ 是根据公式@f[
    A^K_{ij} = \sum_q w_q [\varepsilon(\varphi_i(\mathbf{x}_q)) : C :
                           \varepsilon(\varphi_j(\mathbf{x}_q))],
  @f]使用数值正交计算出来的。

  与正交点 $\mathbf{x}_q$ 和权重 $w_q$  。我们之前在步骤8和步骤17中建立了这些贡献，但在这两种情况下，我们都是通过使用等级4张量 $C$ 的组成知识，以及考虑应变张量的单个元素 $\varepsilon(\varphi_i),\varepsilon(\varphi_j)$ ，相当笨拙地完成的。这其实并不方便，特别是如果我们想考虑比各向同性的情况更复杂的弹性模型，而 $C$ 有方便的形式  $C_{ijkl}  = \lambda \delta_{ij} \delta_{kl} + \mu (\delta_{ik} \delta_{jl}
  + \delta_{il} \delta_{jk})$  。虽然我们在本程序中没有使用比这更复杂的形式，但我们还是希望以一种容易实现的方式来编写它。因此，很自然地要引入代表等级为2（用于应变和应力）和4（用于应力-应变张量 $C$ ）的对称张量的类。幸运的是，deal.II提供了这些： <code>SymmetricTensor<rank,dim></code> 类模板提供了等级 <code>rank</code> （需要是偶数）和维度 <code>dim</code> 的这类张量的完整实现。

  然后我们需要的是两件事：一种创建应力-应变等级4张量 $C$ 的方法，以及从形状函数 $\varphi_i$ 的梯度在给定单元上的正交点 $\mathbf{x}_q$ 创建一个等级2的对称张量（应变张量）。在这个例子程序的执行顶部，你会发现这样的函数。第一个， <code>get_stress_strain_tensor</code>  ，需要两个参数，对应于Lam&eacute; 常数 $\lambda$ 和 $\mu$ ，并返回对应于这些常数的各向同性的应力应变张量（在程序中，我们将选择对应于钢的常数）；用一个计算各向异性的张量的函数来代替这个函数是很简单的，或者考虑到晶体对称性，比如说。第二个， <code>get_strain</code> takes an object of type <code>FEValues</code> 和指数 $i$ 和 $q$ ，返回对称梯度，即应变，对应于形状函数 $\varphi_i(\mathbf{x}_q)$ ，在 <code>FEValues</code> 对象最后被重新初始化的单元上评估。

  鉴于此， <code>assemble_system</code> 的最内部循环以下列优雅的方式计算对矩阵的局部贡献（变量 <code>stress_strain_tensor</code> ，对应于张量 $C$ ，之前已经用上述第一个函数的结果初始化了）。   @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  for (unsigned int j=0; j<dofs_per_cell; ++j)
    for (unsigned int q_point=0; q_point<n_q_points;
         ++q_point)
      {
        const SymmetricTensor<2,dim>
          eps_phi_i = get_strain (fe_values, i, q_point),
          eps_phi_j = get_strain (fe_values, j, q_point);


        cell_matrix(i,j)
          += (eps_phi_i * stress_strain_tensor * eps_phi_j *
              fe_values.JxW (q_point));
      }
  @endcode

  值得注意的是这段代码的表现力，并将其与我们在以前的例子中为弹性问题所经历的复杂情况进行比较。公平地说，在写这些以前的例子时，SymmetricTensor类模板还不存在）。为了简单起见， <code>operator*</code> 在这里规定了偶数等级的对称张量之间的（双重求和）积。

  组建本地捐款@f{eqnarray*}
      f^K_i &=&
      (\mathbf{f}, \varphi_i)_K -(\sigma^{n-1},\varepsilon(\varphi_i))_K
      \\
      &\approx&
      \sum_q
      w_q \left\{
        \mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q) -
        \sigma^{n-1}_q : \varepsilon(\varphi_i(\mathbf{x}_q))
      \right\}
  @f}

  到<a href="#step_18.linear-system">[linear-system]</a>的右手边同样是直接的（注意，我们在这里不考虑任何边界牵引 $\mathbf{b}$ ）。请记住，我们只需要在单元格的正交点上存储旧的应力。在程序中，我们将提供一个变量 <code>local_quadrature_points_data</code> ，允许访问每个正交点的应力 $\sigma^{n-1}_q$ 。有了这个，右手边的代码看起来就像这样，同样相当优雅。   @code
for (unsigned int i=0; i<dofs_per_cell; ++i)
  {
    const unsigned int
      component_i = fe.system_to_component_index(i).first;


    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
      {
        const SymmetricTensor<2,dim> &old_stress
          = local_quadrature_points_data[q_point].old_stress;


        cell_rhs(i) += (body_force_values[q_point](component_i) *
                        fe_values.shape_value (i,q_point)


                        -
                        old_stress *
                        get_strain (fe_values,i,q_point)) *
                       fe_values.JxW (q_point);
      }
  }
  @endcode

  请注意，在乘法 $\mathbf{f}(\mathbf{x}_q) \cdot \varphi_i(\mathbf{x}_q)$ 中，我们利用了这样一个事实：对于所选的有限元素， $\varphi_i$ 中只有一个向量分量（即 <code>component_i</code> ）是非零的，因此我们也只需要考虑 $\mathbf{f}(\mathbf{x}_q)$ 的一个分量。

  这基本上结束了我们在这个函数中提出的新材料。它后来必须处理边界条件以及悬挂节点约束，但这与我们以前在其他程序中已经要做的事情相类似。

 <li>   <code>solve_linear_problem ()</code> [via <code>solve_timestep ()</code>  ] 。   与前一个函数不同，这个函数其实并不有趣，因为它做的是以前所有教程程序中的类似函数--用CG方法求解线性系统，使用不完整的LU分解作为预处理程序（在%并行情况下，它分别使用每个处理器块的ILU）。它与第17步几乎没有变化。

 <li>   <code>update_quadrature_point_history ()</code>  [通过 <code>solve_timestep ()</code>  ] 。基于之前计算的位移场 $\Delta \mathbf{u}^n$ ，我们根据<a href="#step_18.stress-update">[stress-update]</a>和<a href="#step_18.stress-update+rot">[stress-update+rot]</a>更新所有正交点的应力值，包括坐标系的旋转。

 <li>   <code>move_mesh ()</code>  ：给定之前计算的解决方案，在这个函数中，我们通过移动每个顶点的位移矢量场来实现网格的变形。

 <li>   <code>output_results ()</code>  : 这个函数只是根据我们上面所说的输出解决方案，也就是说，每个处理器只对自己的那部分域计算输出。除了解决方案，我们还计算了每个单元上所有正交点平均的应力的规范。   </ul> 

有了这个代码的一般结构，我们只需要定义我们要解决的情况。在本程序中，我们选择模拟一个垂直圆柱体的准静态变形，其底部边界是固定的，顶部边界以规定的垂直速度被推倒。然而，顶层边界的水平速度没有被指定--我们可以把这种情况想象成一块油性良好的板从顶部推到圆柱体上，圆柱体顶层边界上的点被允许沿着板的表面水平滑动，但被板强迫向下移动。圆柱体的内部和外部边界是自由的，不受任何规定的偏转或牵引的影响。此外，重力作用于身体。

程序文本将揭示更多关于如何实现这种情况，而结果部分将显示这种模拟产生的位移模式。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 首先是通常的头文件列表，这些文件已经在以前的示例程序中使用过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/multithread_info.h> 
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/petsc_vector.h> 
 * #include <deal.II/lac/petsc_sparse_matrix.h> 
 * #include <deal.II/lac/petsc_solver.h> 
 * #include <deal.II/lac/petsc_precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/sparsity_tools.h> 
 * #include <deal.II/distributed/shared_tria.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 这里是头文件中仅有的三个新东西：一个包含文件，其中实现了等级为2和4的对称张量，正如介绍中所介绍的那样。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/symmetric_tensor.h> 
 * 
 * @endcode
 * 
 * 最后是一个包含一些函数的头文件，这些函数将帮助我们计算域中特定点的局部坐标系的旋转矩阵。
 * 

 * 
 * 
 * @code
 * #include <deal.II/physics/transformations.h> 
 * 
 * @endcode
 * 
 * 然后，这又是简单的C++。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * #include <iomanip> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step18 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodePointHistorycodeclass"></a> 
 * <h3>The <code>PointHistory</code> class</h3>
 * 

 * 
 * 正如介绍中提到的，我们必须在正交点存储旧的应力，这样我们就可以在下一个时间步骤中计算这一点的残余力。仅仅这一点还不能保证只有一个成员的结构，但在更复杂的应用中，我们还必须在正交点上存储更多的信息，比如塑性的历史变量等。从本质上讲，我们必须在这里存储所有影响材料当前状态的信息，在塑性中，这些信息是由变形历史变量决定的。
 * 

 * 
 * 除了能够存储数据之外，我们不会给这个类任何有意义的功能，也就是说，没有构造函数、析构函数或其他成员函数。在这种 "哑巴 "类的情况下，我们通常选择将其声明为  <code>struct</code> rather than <code>class</code>  ，以表明它们更接近于C语言风格的结构而不是C++风格的类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct PointHistory 
 *   { 
 *     SymmetricTensor<2, dim> old_stress; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Thestressstraintensor"></a> 
 * <h3>The stress-strain tensor</h3>
 * 

 * 
 * 接下来，我们定义弹性中的应力和应变的线性关系。它由一个等级为4的张量给出，通常被写成  $C_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) + \lambda \delta_{ij} \delta_{kl}$  的形式。这个张量将等级2的对称张量映射到等级2的对称张量。对于Lam&eacute;常数 $\lambda$ 和 $\mu$ 的给定值，一个实现其创建的函数是直接的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda, 
 *                                                    const double mu) 
 *   { 
 *     SymmetricTensor<4, dim> tmp; 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       for (unsigned int j = 0; j < dim; ++j) 
 *         for (unsigned int k = 0; k < dim; ++k) 
 *           for (unsigned int l = 0; l < dim; ++l) 
 *             tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) + 
 *                                ((i == l) && (j == k) ? mu : 0.0) + 
 *                                ((i == j) && (k == l) ? lambda : 0.0)); 
 *     return tmp; 
 *   } 
 * 
 * @endcode
 * 
 * 通过这个函数，我们将在下面的主类中定义一个静态成员变量，在整个程序中作为应力-应变张量使用。请注意，在更复杂的程序中，这可能是某个类的成员变量，或者是一个根据其他输入返回应力-应变关系的函数。例如，在损伤理论模型中，Lam&eacute;常数被认为是一个点的先前应力/应变历史的函数。相反，在塑性中，如果材料在某一点达到了屈服应力，那么应力-应变张量的形式就会被修改，而且可能还取决于其先前的历史。
 * 

 * 
 * 然而，在本程序中，我们假设材料是完全弹性和线性的，恒定的应力-应变张量对我们目前的目的来说是足够的。
 * 

 * 
 * 
 * <a name="Auxiliaryfunctions"></a> 
 * <h3>Auxiliary functions</h3>
 * 

 * 
 * 在程序的其他部分之前，这里有几个我们需要的函数作为工具。这些是在内循环中调用的小函数，所以我们把它们标记为  <code>inline</code>  。
 * 

 * 
 * 第一个是通过形成这个形状函数的对称梯度来计算形状函数 <code>shape_func</code> at quadrature point <code>q_point</code> 的对称应变张量。当我们想形成矩阵时，我们需要这样做，比如说。
 * 

 * 
 * 我们应该注意到，在以前处理矢量值问题的例子中，我们总是问有限元对象在哪个矢量分量中的形状函数实际上是不为零的，从而避免计算任何我们反正可以证明为零的项。为此，我们使用了 <code>fe.system_to_component_index</code> 函数来返回形状函数在哪个分量中为零，同时 <code>fe_values.shape_value</code> 和 <code>fe_values.shape_grad</code> 函数只返回形状函数的单个非零分量的值和梯度，如果这是一个矢量值元素。
 * 

 * 
 * 这是一个优化，如果不是非常关键的时间，我们可以用一个更简单的技术来解决：只需向 <code>fe_values</code> 询问一个给定形状函数的给定分量在给定正交点的值或梯度。这就是  <code>fe_values.shape_grad_component(shape_func,q_point,i)</code>  调用的作用：返回形状函数  <code>shape_func</code>  的第  <code>q_point</code>  个分量在正交点的全部梯度。如果某个形状函数的某个分量总是为零，那么这将简单地总是返回零。
 * 

 * 
 * 如前所述，使用 <code>fe_values.shape_grad_component</code> 而不是 <code>fe.system_to_component_index</code> 和 <code>fe_values.shape_grad</code> 的组合可能效率较低，但其实现已针对这种情况进行了优化，应该不会有很大的减慢。我们在这里演示这个技术，因为它是如此的简单和直接。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values, 
 *                                             const unsigned int   shape_func, 
 *                                             const unsigned int   q_point) 
 *   { 
 * 
 * @endcode
 * 
 * 声明一个将保存返回值的暂存器。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<2, dim> tmp; 
 * 
 * @endcode
 * 
 * 首先，填充对角线项，这只是矢量值形状函数的方向 <code>i</code> of the <code>i</code> 分量的导数。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i]; 
 * 
 * @endcode
 * 
 * 然后填充应变张量的其余部分。注意，由于张量是对称的，我们只需要计算一半（这里：右上角）的非对角线元素， <code>SymmetricTensor</code> 类的实现确保至少到外面的对称条目也被填充（实际上，这个类当然只存储一份）。在这里，我们选择了张量的右上半部分，但是左下半部分也一样好。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       for (unsigned int j = i + 1; j < dim; ++j) 
 *         tmp[i][j] = 
 *           (fe_values.shape_grad_component(shape_func, q_point, i)[j] + 
 *            fe_values.shape_grad_component(shape_func, q_point, j)[i]) / 
 *           2; 
 * 
 *     return tmp; 
 *   } 
 * 
 * @endcode
 * 
 * 第二个函数做了非常类似的事情（因此被赋予相同的名字）：从一个矢量值场的梯度计算对称应变张量。如果你已经有了一个解场， <code>fe_values.get_function_gradients</code> 函数允许你在一个正交点上提取解场的每个分量的梯度。它返回的是一个秩-1张量的矢量：解的每个矢量分量有一个秩-1张量（梯度）。由此，我们必须通过转换数据存储格式和对称化来重建（对称的）应变张量。我们用和上面一样的方法来做，也就是说，我们通过首先填充对角线，然后只填充对称张量的一半来避免一些计算（ <code>SymmetricTensor</code> 类确保只写两个对称分量中的一个就足够了）。
 * 

 * 
 * 不过在我们这样做之前，我们要确保输入有我们期望的那种结构：即有 <code>dim</code> 个矢量分量，即每个坐标方向有一个位移分量。我们用 <code>Assert</code> 宏来测试这一点，如果不符合条件，我们的程序就会被终止。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   inline SymmetricTensor<2, dim> 
 *   get_strain(const std::vector<Tensor<1, dim>> &grad) 
 *   { 
 *     Assert(grad.size() == dim, ExcInternalError()); 
 * 
 *     SymmetricTensor<2, dim> strain; 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       strain[i][i] = grad[i][i]; 
 * 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       for (unsigned int j = i + 1; j < dim; ++j) 
 *         strain[i][j] = (grad[i][j] + grad[j][i]) / 2; 
 * 
 *     return strain; 
 *   } 
 * 
 * @endcode
 * 
 * 最后，下面我们将需要一个函数来计算某一点的位移所引起的旋转矩阵。当然，事实上，单点的位移只有一个方向和一个幅度，诱发旋转的是方向和幅度的变化。实际上，旋转矩阵可以通过位移的梯度来计算，或者更具体地说，通过卷曲来计算。
 * 

 * 
 * 确定旋转矩阵的公式有点笨拙，特别是在三维中。对于2D来说，有一个更简单的方法，所以我们把这个函数实现了两次，一次用于2D，一次用于3D，这样我们就可以在两个空间维度上编译和使用这个程序，如果需要的话--毕竟，deal.II是关于独立维度编程和重复使用算法的，在2D的廉价计算中经过测试，在3D的更昂贵的计算中使用。下面是一种情况，我们必须为2D和3D实现不同的算法，但可以用独立于空间维度的方式来编写程序的其余部分。
 * 

 * 
 * 所以，不用再多说了，来看看2D的实现。
 * 

 * 
 * 
 * @code
 *   Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>> &grad_u) 
 *   { 
 * 
 * @endcode
 * 
 * 首先，根据梯度计算出速度场的卷曲。注意，我们是在2d中，所以旋转是一个标量。
 * 

 * 
 * 
 * @code
 *     const double curl = (grad_u[1][0] - grad_u[0][1]); 
 * 
 * @endcode
 * 
 * 由此计算出旋转的角度。
 * 

 * 
 * 
 * @code
 *     const double angle = std::atan(curl); 
 * 
 * @endcode
 * 
 * 由此，建立反对称的旋转矩阵。我们希望这个旋转矩阵能够代表本地坐标系相对于全局直角坐标系的旋转，所以我们用一个负的角度来构建它。因此，这个旋转矩阵代表了从本地坐标系移动到全局坐标系所需的旋转。
 * 

 * 
 * 
 * @code
 *     return Physics::Transformations::Rotations::rotation_matrix_2d(-angle); 
 *   } 
 * 
 * @endcode
 * 
 * 三维的情况就比较复杂了。
 * 

 * 
 * 
 * @code
 *   Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u) 
 *   { 
 * 
 * @endcode
 * 
 * 同样首先计算速度场的卷曲。这一次，它是一个实数向量。
 * 

 * 
 * 
 * @code
 *     const Point<3> curl(grad_u[2][1] - grad_u[1][2], 
 *                         grad_u[0][2] - grad_u[2][0], 
 *                         grad_u[1][0] - grad_u[0][1]); 
 * 
 * @endcode
 * 
 * 从这个矢量中，利用它的大小，计算出旋转角度的正切值，并由此计算出相对于直角坐标系的实际旋转角度。
 * 

 * 
 * 
 * @code
 *     const double tan_angle = std::sqrt(curl * curl); 
 *     const double angle     = std::atan(tan_angle); 
 * 
 * @endcode
 * 
 * 现在，这里有一个问题：如果旋转角度太小，那就意味着没有旋转发生（例如平移运动）。在这种情况下，旋转矩阵就是身份矩阵。
 * 

 * 
 * 我们强调这一点的原因是，在这种情况下，我们有  <code>tan_angle==0</code>  。再往下看，我们在计算旋转轴的时候需要除以这个数字，这样做除法的时候会遇到麻烦。因此，让我们走捷径，如果旋转角度真的很小，就简单地返回同一矩阵。
 * 

 * 
 * 
 * @code
 *     if (std::abs(angle) < 1e-9) 
 *       { 
 *         static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; 
 *         static const Tensor<2, 3> rot(rotation); 
 *         return rot; 
 *       } 
 * 
 * @endcode
 * 
 * 否则计算真实的旋转矩阵。为此，我们再次依靠一个预定义的函数来计算本地坐标系的旋转矩阵。
 * 

 * 
 * 
 * @code
 *     const Point<3> axis = curl / tan_angle; 
 *     return Physics::Transformations::Rotations::rotation_matrix_3d(axis, 
 *                                                                    -angle); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeTopLevelcodeclass"></a> 
 * <h3>The <code>TopLevel</code> class</h3>
 * 

 * 
 * 这就是程序的主类。由于命名空间已经表明了我们要解决的问题，让我们用它的作用来称呼它：它引导着程序的流程，也就是说，它是顶层驱动。
 * 

 * 
 * 这个类的成员变量基本上和以前一样，即它必须有一个三角形，一个DoF处理程序和相关的对象，如约束条件，描述线性系统的变量等。现在还有很多成员函数，我们将在下面解释。
 * 

 * 
 * 然而，该类的外部接口是不变的：它有一个公共的构造函数和析构函数，并且它有一个 <code>run</code> 函数来启动所有的工作。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class TopLevel 
 *   { 
 *   public: 
 *     TopLevel(); 
 *     ~TopLevel(); 
 *     void run(); 
 * 
 *   private: 
 * 
 * @endcode
 * 
 * 私有接口比  step-17  中的更加广泛。首先，我们显然需要创建初始网格的函数，设置描述当前网格上的线性系统的变量（即矩阵和向量），然后是实际组装系统的函数，指导每个时间步长中必须解决的问题，一个解决每个时间步长中出现的线性系统的函数（并返回它的迭代次数），最后在正确的网格上输出解向量。
 * 

 * 
 * 
 * @code
 *     void create_coarse_grid(); 
 * 
 *     void setup_system(); 
 * 
 *     void assemble_system(); 
 * 
 *     void solve_timestep(); 
 * 
 *     unsigned int solve_linear_problem(); 
 * 
 *     void output_results() const; 
 * 
 * @endcode
 * 
 * 除了前两个，所有这些函数都在每个时间步中被调用。由于第一个时间步骤有点特殊，我们有单独的函数来描述一个时间步骤中必须发生的事情：一个用于第一个时间步骤，一个用于所有后续时间步骤。
 * 

 * 
 * 
 * @code
 *     void do_initial_timestep(); 
 * 
 *     void do_timestep(); 
 * 
 * @endcode
 * 
 * 然后我们需要一大堆函数来做各种事情。第一个是细化初始网格：我们从原始状态的粗网格开始，解决这个问题，然后看一下，并相应地细化网格，然后重新开始同样的过程，再次以原始状态。因此，细化初始网格比在两个连续的时间步骤之间细化网格要简单一些，因为它不涉及将数据从旧的三角测量转移到新的三角测量，特别是存储在每个正交点的历史数据。
 * 

 * 
 * 
 * @code
 *     void refine_initial_grid(); 
 * 
 * @endcode
 * 
 * 在每个时间步骤结束时，我们要根据这个时间步骤计算的增量位移来移动网格顶点。这就是完成这个任务的函数。
 * 

 * 
 * 
 * @code
 *     void move_mesh(); 
 * 
 * @endcode
 * 
 * 接下来是两个处理存储在每个正交点的历史变量的函数。第一个函数在第一个时间步长之前被调用，为历史变量设置一个原始状态。它只对属于当前处理器的单元上的正交点起作用。
 * 

 * 
 * 
 * @code
 *     void setup_quadrature_point_history(); 
 * 
 * @endcode
 * 
 * 第二项是在每个时间段结束时更新历史变量。
 * 

 * 
 * 
 * @code
 *     void update_quadrature_point_history(); 
 * 
 * @endcode
 * 
 * 这是新的共享三角法。
 * 

 * 
 * 
 * @code
 *     parallel::shared::Triangulation<dim> triangulation; 
 * 
 *     FESystem<dim> fe; 
 * 
 *     DoFHandler<dim> dof_handler; 
 * 
 *     AffineConstraints<double> hanging_node_constraints; 
 * 
 * @endcode
 * 
 * 这个程序的一个不同之处在于，我们在类声明中声明了正交公式。原因是在所有其他程序中，如果我们在计算矩阵和右手边时使用不同的正交公式，并没有什么坏处，比如说。然而，在目前的情况下，它确实如此：我们在正交点中存储了信息，所以我们必须确保程序的所有部分都同意它们的位置以及每个单元格上有多少个。因此，让我们首先声明将在整个程序中使用的正交公式...。
 * 

 * 
 * 
 * @code
 *     const QGauss<dim> quadrature_formula; 
 * 
 * @endcode
 * 
 * ......然后也有一个历史对象的向量，在我们负责的那些单元格上的每个正交点都有一个（也就是说，我们不为其他处理器拥有的单元格上的正交点存储历史数据）。请注意，我们可以像在  step-44  中那样使用 CellDataStorage 类来代替我们自己存储和管理这些数据。然而，为了演示的目的，在这种情况下，我们手动管理存储。
 * 

 * 
 * 
 * @code
 *     std::vector<PointHistory<dim>> quadrature_point_history; 
 * 
 * @endcode
 * 
 * 这个对象的访问方式是通过每个单元格、面或边持有的 <code>user pointer</code> ：它是一个 <code>void*</code> 指针，可以被应用程序用来将任意的数据与单元格、面或边联系起来。程序对这些数据的实际操作属于自己的职责范围，库只是为这些指针分配了一些空间，而应用程序可以设置和读取这些对象中的每个指针。
 * 

 * 
 * 进一步说：我们需要待解的线性系统的对象，即矩阵、右手边的向量和解向量。由于我们预计要解决大问题，我们使用了与 step-17 中相同的类型，即建立在PETSc库之上的分布式%并行矩阵和向量。方便的是，它们也可以在只在一台机器上运行时使用，在这种情况下，这台机器正好是我们的%并行宇宙中唯一的机器。
 * 

 * 
 * 然而，与 step-17 不同的是，我们不以分布式方式存储解向量--这里是在每个时间步骤中计算的增量位移。也就是说，在计算时它当然必须是一个分布式矢量，但紧接着我们确保每个处理器都有一个完整的副本。原因是我们已经在 step-17 中看到，许多函数需要一个完整的副本。虽然得到它并不难，但这需要在网络上进行通信，因此很慢。此外，这些都是重复的相同操作，这当然是不可取的，除非不必总是存储整个向量的收益超过了它。在编写这个程序时，事实证明，我们在很多地方都需要一份完整的解决方案，以至于只在必要时才获得它似乎不值得。相反，我们选择一劳永逸地获得完整的副本，而立即摆脱分散的副本。因此，请注意， <code>incremental_displacement</code> 的声明并没有像中间命名空间 <code>MPI</code> 所表示的那样，表示一个分布式向量。
 * 

 * 
 * 
 * @code
 *     PETScWrappers::MPI::SparseMatrix system_matrix; 
 * 
 *     PETScWrappers::MPI::Vector system_rhs; 
 * 
 *     Vector<double> incremental_displacement; 
 * 
 * @endcode
 * 
 * 接下来的变量块与问题的时间依赖性有关：它们表示我们要模拟的时间间隔的长度，现在的时间和时间步数，以及现在时间步数的长度。
 * 

 * 
 * 
 * @code
 *     double       present_time; 
 *     double       present_timestep; 
 *     double       end_time; 
 *     unsigned int timestep_no; 
 * 
 * @endcode
 * 
 * 然后是几个与%并行处理有关的变量：首先，一个变量表示我们使用的MPI通信器，然后是两个数字，告诉我们有多少个参与的处理器，以及我们在这个世界上的位置。最后，一个流对象，确保只有一个处理器实际产生输出到控制台。这与  step-17  中的所有内容相同。
 * 

 * 
 * 
 * @code
 *     MPI_Comm mpi_communicator; 
 * 
 *     const unsigned int n_mpi_processes; 
 * 
 *     const unsigned int this_mpi_process; 
 * 
 *     ConditionalOStream pcout; 
 * 
 * @endcode
 * 
 * 我们正在存储本地拥有的和本地相关的索引。
 * 

 * 
 * 
 * @code
 *     IndexSet locally_owned_dofs; 
 *     IndexSet locally_relevant_dofs; 
 * 
 * @endcode
 * 
 * 最后，我们有一个静态变量，表示应力和应变之间的线性关系。由于它是一个不依赖任何输入的常量对象（至少在这个程序中不依赖），我们把它作为一个静态变量，并将在我们定义这个类的构造函数的同一个地方初始化它。
 * 

 * 
 * 
 * @code
 *     static const SymmetricTensor<4, dim> stress_strain_tensor; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ThecodeBodyForcecodeclass"></a> 
 * <h3>The <code>BodyForce</code> class</h3>
 * 

 * 
 * 在我们进入这个程序的主要功能之前，我们必须定义哪些力将作用在我们想要研究的变形的体上。这些力可以是体力，也可以是边界力。体力通常是由四种基本的物理力类型之一所介导的：重力、强弱相互作用和电磁力。除非人们想考虑亚原子物体（对于这些物体，无论如何准静态变形是不相关的，也是不合适的描述），否则只需要考虑引力和电磁力。为了简单起见，让我们假设我们的身体有一定的质量密度，但要么是非磁性的，不导电的，要么周围没有明显的电磁场。在这种情况下，身体的力只是 <code>rho g</code>, where <code>rho</code> 是材料密度， <code>g</code> 是一个负Z方向的矢量，大小为9.81米/秒^2。 密度和 <code>g</code> 都是在函数中定义的，我们把7700 kg/m^3作为密度，这是对钢材通常假定的值。
 * 

 * 
 * 为了更普遍一点，也为了能够在2d中进行计算，我们意识到体力总是一个返回 <code>dim</code> 维矢量的函数。我们假设重力沿着最后一个，即 <code>dim-1</code> 个坐标的负方向作用。考虑到以前的例子程序中的类似定义，这个函数的其余实现应该大部分是不言自明的。请注意，身体的力量与位置无关；为了避免编译器对未使用的函数参数发出警告，我们因此注释了 <code>vector_value</code> 函数的第一个参数的名称。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BodyForce : public Function<dim> 
 *   { 
 *   public: 
 *     BodyForce(); 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  values) const override; 
 * 
 *     virtual void 
 *     vector_value_list(const std::vector<Point<dim>> &points, 
 *                       std::vector<Vector<double>> &  value_list) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   BodyForce<dim>::BodyForce() 
 *     : Function<dim>(dim) 
 *   {} 
 * 
 *   template <int dim> 
 *   inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/, 
 *                                            Vector<double> &values) const 
 *   { 
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 
 * 
 *     const double g   = 9.81; 
 *     const double rho = 7700; 
 * 
 *     values          = 0; 
 *     values(dim - 1) = -rho * g; 
 *   } 
 * 
 *   template <int dim> 
 *   void BodyForce<dim>::vector_value_list( 
 *     const std::vector<Point<dim>> &points, 
 *     std::vector<Vector<double>> &  value_list) const 
 *   { 
 *     const unsigned int n_points = points.size(); 
 * 
 *     Assert(value_list.size() == n_points, 
 *            ExcDimensionMismatch(value_list.size(), n_points)); 
 * 
 *     for (unsigned int p = 0; p < n_points; ++p) 
 *       BodyForce<dim>::vector_value(points[p], value_list[p]); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeIncrementalBoundaryValuecodeclass"></a> 
 * <h3>The <code>IncrementalBoundaryValue</code> class</h3>
 * 

 * 
 * 除了身体的力之外，运动还可以由边界力和强制边界位移引起。后一种情况相当于以这样的方式选择力，使其诱发某种位移。
 * 

 * 
 * 对于准静态位移，典型的边界力是对一个体的压力，或者对另一个体的切向摩擦。我们在这里选择了一种更简单的情况：我们规定了边界（部分）的某种运动，或者至少是位移矢量的某些分量。我们用另一个矢量值函数来描述，对于边界上的某一点，返回规定的位移。
 * 

 * 
 * 由于我们有一个随时间变化的问题，边界的位移增量等于在时间段内累积的位移。因此，该类必须同时知道当前时间和当前时间步长，然后可以将位移增量近似为当前速度乘以当前时间步长。
 * 

 * 
 * 在本程序中，我们选择了一种简单的边界位移形式：我们以恒定的速度向下位移顶部的边界。边界的其余部分要么是固定的（然后用一个 <code>Functions::ZeroFunction</code> 类型的对象来描述），要么是自由的（Neumann类型，在这种情况下不需要做任何特殊的事情）。 利用我们在前面所有的例子程序中获得的知识，描述持续向下运动的类的实现应该是很明显的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class IncrementalBoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     IncrementalBoundaryValues(const double present_time, 
 *                               const double present_timestep); 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  values) const override; 
 * 
 *     virtual void 
 *     vector_value_list(const std::vector<Point<dim>> &points, 
 *                       std::vector<Vector<double>> &  value_list) const override; 
 * 
 *   private: 
 *     const double velocity; 
 *     const double present_time; 
 *     const double present_timestep; 
 *   }; 
 * 
 *   template <int dim> 
 *   IncrementalBoundaryValues<dim>::IncrementalBoundaryValues( 
 *     const double present_time, 
 *     const double present_timestep) 
 *     : Function<dim>(dim) 
 *     , velocity(.08) 
 *     , present_time(present_time) 
 *     , present_timestep(present_timestep) 
 *   {} 
 * 
 *   template <int dim> 
 *   void 
 *   IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/, 
 *                                                Vector<double> &values) const 
 *   { 
 *     Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 
 * 
 *     values    = 0; 
 *     values(2) = -present_timestep * velocity; 
 *   } 
 * 
 *   template <int dim> 
 *   void IncrementalBoundaryValues<dim>::vector_value_list( 
 *     const std::vector<Point<dim>> &points, 
 *     std::vector<Vector<double>> &  value_list) const 
 *   { 
 *     const unsigned int n_points = points.size(); 
 * 
 *     Assert(value_list.size() == n_points, 
 *            ExcDimensionMismatch(value_list.size(), n_points)); 
 * 
 *     for (unsigned int p = 0; p < n_points; ++p) 
 *       IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeTopLevelcodeclass"></a> 
 * <h3>Implementation of the <code>TopLevel</code> class</h3>
 * 

 * 
 * 现在是主类的实现。首先，我们初始化应力应变张量，我们将其声明为一个静态常量变量。我们选择了适合于钢铁的Lam&eacute;常数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   const SymmetricTensor<4, dim> TopLevel<dim>::stress_strain_tensor = 
 *     get_stress_strain_tensor<dim>(
 *       /*lambda = */ 9.695e10, 
 *       /*mu =  */  7.617e10)
 * 
 * @endcode
 * 
 * 
 * <a name="Thepublicinterface"></a> 
 * <h4>The public interface</h4>
 * 

 * 
 * 下一步是构造函数和析构函数的定义。这里没有什么惊喜：我们为解的每个 <code>dim</code> 矢量分量选择线性和连续的有限元，以及每个坐标方向上有2个点的高斯正交公式。解构器应该是显而易见的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   TopLevel<dim>::TopLevel() 
 *     : triangulation(MPI_COMM_WORLD) 
 *     , fe(FE_Q<dim>(1), dim) 
 *     , dof_handler(triangulation) 
 *     , quadrature_formula(fe.degree + 1) 
 *     , present_time(0.0) 
 *     , present_timestep(1.0) 
 *     , end_time(10.0) 
 *     , timestep_no(0) 
 *     , mpi_communicator(MPI_COMM_WORLD) 
 *     , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)) 
 *     , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)) 
 *     , pcout(std::cout, this_mpi_process == 0) 
 *   {} 
 * 
 *   template <int dim> 
 *   TopLevel<dim>::~TopLevel() 
 *   { 
 *     dof_handler.clear(); 
 *   } 
 * 
 * @endcode
 * 
 * 最后一个公共函数是指导所有工作的函数，  <code>run()</code>  。它初始化了描述我们目前所处时间位置的变量，然后运行第一个时间步骤，再循环所有其他时间步骤。请注意，为了简单起见，我们使用一个固定的时间步长，而一个更复杂的程序当然要以某种更合理的方式自适应地选择它。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::run() 
 *   { 
 *     do_initial_timestep(); 
 * 
 *     while (present_time < end_time) 
 *       do_timestep(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TopLevelcreate_coarse_grid"></a> 
 * <h4>TopLevel::create_coarse_grid</h4>
 * 

 * 
 * 按照上面声明的顺序，下一个函数是创建粗略网格的函数，我们从这里开始。在这个示例程序中，我们想计算一个圆柱体在轴向压缩下的变形。因此第一步是生成一个长度为3，内外半径分别为0.8和1的圆柱体的网格。幸运的是，有一个库函数可以生成这样的网格。
 * 

 * 
 * 在第二步中，我们必须在圆柱体的上表面和下表面关联边界条件。我们为边界面选择一个边界指示器0，这些边界面的中点的Z坐标为0（底面），Z=3的指示器为1（顶面）；最后，我们对圆柱体外壳内部的所有面使用边界指示器2，外部使用3。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::create_coarse_grid() 
 *   { 
 *     const double inner_radius = 0.8, outer_radius = 1; 
 *     GridGenerator::cylinder_shell(triangulation, 3, inner_radius, outer_radius); 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->at_boundary()) 
 *           { 
 *             const Point<dim> face_center = face->center(); 
 * 
 *             if (face_center[2] == 0) 
 *               face->set_boundary_id(0); 
 *             else if (face_center[2] == 3) 
 *               face->set_boundary_id(1); 
 *             else if (std::sqrt(face_center[0] * face_center[0] + 
 *                                face_center[1] * face_center[1]) < 
 *                      (inner_radius + outer_radius) / 2) 
 *               face->set_boundary_id(2); 
 *             else 
 *               face->set_boundary_id(3); 
 *           } 
 * 
 * @endcode
 * 
 * 一旦完成了这些，我们就可以对网格进行一次全面的细化。
 * 

 * 
 * 
 * @code
 *     triangulation.refine_global(1); 
 * 
 * @endcode
 * 
 * 作为最后一步，我们需要设置一个干净的数据状态，我们将这些数据存储在目前处理器上处理的所有单元的正交点中。
 * 

 * 
 * 
 * @code
 *     setup_quadrature_point_history(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_system"></a> 
 * <h4>TopLevel::setup_system</h4>
 * 

 * 
 * 下一个函数是为一个给定的网格设置数据结构。这与 step-17 中的方法基本相同：分配自由度，然后对这些自由度进行排序，使每个处理器得到一个连续的块。请注意，每个处理器的细分块是在创建或完善网格的函数中处理的，与之前的例子程序不同（发生这种情况的时间点主要是口味问题；在这里，我们选择在创建网格时进行，因为在 <code>do_initial_timestep</code> 和 <code>do_timestep</code> 函数中，我们想在还没有调用当前函数的时候输出每个处理器上的单元数量）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
 * 
 * @endcode
 * 
 * 下一步是设置由于悬挂节点而产生的约束。这在以前已经处理过很多次了。
 * 

 * 
 * 
 * @code
 *     hanging_node_constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 * @endcode
 * 
 * 然后我们要设置矩阵。这里我们偏离了  step-17  ，在那里我们简单地使用了PETSc的能力，即只知道矩阵的大小，随后分配那些被写入的非零元素。虽然从正确性的角度来看，这样做很好，但是效率却不高：如果我们不给PETSc提供关于哪些元素被写入的线索，那么当我们第一次设置矩阵中的元素时（即在第一个时间步中），它的速度会慢得令人难以忍受。后来，当元素被分配后，一切都快多了。在我们所做的实验中，如果我们指示PETSc哪些元素将被使用，哪些不被使用，那么第一个时间步骤可以加快近两个数量级。
 * 

 * 
 * 要做到这一点，我们首先要生成我们要处理的矩阵的稀疏模式，并确保浓缩的悬挂节点约束在稀疏模式中增加必要的额外条目。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs); 
 *     DoFTools::make_sparsity_pattern(dof_handler, 
 *                                     sparsity_pattern, 
 *                                     hanging_node_constraints, 
 *                                     /*保持约束性dofs  */ false)
 * 
 *     SparsityTools::distribute_sparsity_pattern(sparsity_pattern, 
 *                                                locally_owned_dofs, 
 *                                                mpi_communicator, 
 *                                                locally_relevant_dofs); 
 * 
 * @endcode
 * 
 * 注意，我们在这里使用了已经在 step-11 中介绍过的 <code>DynamicSparsityPattern</code> 类，而不是我们在所有其他情况下使用的 <code>SparsityPattern</code> 类。其原因是，为了使后一个类发挥作用，我们必须给每一行的条目数提供一个初始的上限，这项任务传统上是由 <code>DoFHandler::max_couplings_between_dofs()</code> 完成。然而，这个函数有一个严重的问题：它必须计算每一行中非零项的数量的上限，而这是一个相当复杂的任务，特别是在3D中。实际上，虽然它在2D中相当准确，但在3D中经常得出太大的数字，在这种情况下， <code>SparsityPattern</code> 一开始就分配了太多的内存，经常是几百MB。后来当 <code>DoFTools::make_sparsity_pattern</code> 被调用时，我们意识到我们不需要那么多的内存，但这时已经太晚了：对于大问题，临时分配太多的内存会导致内存不足的情况。
 * 

 * 
 * 为了避免这种情况，我们采用了 <code>DynamicSparsityPattern</code> 类，该类速度较慢，但不需要预先估计每行非零条目的数量。因此，它在任何时候都只分配它所需要的内存，而且我们甚至可以为大型的三维问题建立它。
 * 

 * 
 * 值得注意的是，由于 parallel::shared::Triangulation, 的特殊性，我们构建的稀疏模式是全局的，即包括所有的自由度，无论它们是属于我们所在的处理器还是另一个处理器（如果这个程序是通过MPI并行运行的）。这当然不是最好的--它限制了我们可以解决的问题的规模，因为在每个处理器上存储整个稀疏模式（即使只是短时间）的规模并不大。然而，在程序中还有几个地方我们是这样做的，例如，我们总是把全局三角测量和DoF处理对象保留在周围，即使我们只对它们的一部分进行工作。目前，deal.II没有必要的设施来完全分配这些对象（事实上，这项任务在自适应网格中很难实现，因为随着网格的自适应细化，领域的均衡分区往往会变得不均衡）。
 * 

 * 
 * 有了这个数据结构，我们就可以进入PETSc稀疏矩阵，告诉它预先分配所有我们以后要写入的条目。
 * 

 * 
 * 
 * @code
 *     system_matrix.reinit(locally_owned_dofs, 
 *                          locally_owned_dofs, 
 *                          sparsity_pattern, 
 *                          mpi_communicator); 
 * 
 * @endcode
 * 
 * 在这一点上，不再需要对稀疏模式有任何明确的了解，我们可以让 <code>sparsity_pattern</code> 这个变量离开范围，不会有任何问题。
 * 

 * 
 * 这个函数的最后一个任务是将右侧向量和求解向量重置为正确的大小；记住，求解向量是一个本地向量，不像右侧向量是一个分布式的%并行向量，因此需要知道MPI通信器，它应该通过这个通信器来传输消息。
 * 

 * 
 * 
 * @code
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
 *     incremental_displacement.reinit(dof_handler.n_dofs()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelassemble_system"></a> 
 * <h4>TopLevel::assemble_system</h4>
 * 

 * 
 * 同样，组装系统矩阵和右手边的结构与之前许多例子程序中的结构相同。特别是，它主要等同于 step-17 ，除了不同的右手边，现在只需要考虑到内部应力。此外，通过使用 <code>SymmetricTensor</code> 类，组装矩阵明显变得更加透明：请注意形成2级和4级对称张量的标量积的优雅性。这个实现也更加通用，因为它与我们可能使用或不使用各向同性的弹性张量这一事实无关。
 * 

 * 
 * 汇编程序的第一部分和以往一样。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::assemble_system() 
 *   { 
 *     system_rhs    = 0; 
 *     system_matrix = 0; 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     BodyForce<dim>              body_force; 
 *     std::vector<Vector<double>> body_force_values(n_q_points, 
 *                                                   Vector<double>(dim)); 
 * 
 * @endcode
 * 
 * 如同在  step-17  中一样，我们只需要在属于当前处理器的所有单元中进行循环。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           cell_matrix = 0; 
 *           cell_rhs    = 0; 
 * 
 *           fe_values.reinit(cell); 
 * 
 * @endcode
 * 
 * 然后在所有指数i,j和正交点上循环，并从这个单元中组合出系统矩阵的贡献。 注意我们如何从 <code>FEValues</code> 对象中提取给定正交点的形状函数的对称梯度（应变），以及我们如何优雅地形成三重收缩 <code>eps_phi_i : C : eps_phi_j</code> ；后者需要与 step-17 中需要的笨拙计算进行比较，无论是在介绍中还是在程序的相应位置。
 * 

 * 
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *                 { 
 *                   const SymmetricTensor<2, dim> 
 *                     eps_phi_i = get_strain(fe_values, i, q_point), 
 *                     eps_phi_j = get_strain(fe_values, j, q_point); 
 * 
 *                   cell_matrix(i, j) += (eps_phi_i *            // 
 *                                         stress_strain_tensor * // 
 *                                         eps_phi_j              // 
 *                                         ) *                    // 
 *                                        fe_values.JxW(q_point); // 
 *                 } 
 * 
 * @endcode
 * 
 * 然后也要组装本地的右手边贡献。为此，我们需要访问这个正交点的先验应力值。为了得到它，我们使用该单元的用户指针，该指针指向全局数组中与当前单元的第一个正交点相对应的正交点数据，然后添加一个与我们现在考虑的正交点的索引相对应的偏移量。
 * 

 * 
 * 
 * @code
 *           const PointHistory<dim> *local_quadrature_points_data = 
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 
 * 
 * @endcode
 * 
 * 此外，我们还需要这个单元上的正交点的外体力值。
 * 

 * 
 * 
 * @code
 *           body_force.vector_value_list(fe_values.get_quadrature_points(), 
 *                                        body_force_values); 
 * 
 * @endcode
 * 
 * 然后，我们可以循环计算这个单元上的所有自由度，并计算出对右侧的局部贡献。
 * 

 * 
 * 
 * @code
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const unsigned int component_i = 
 *                 fe.system_to_component_index(i).first; 
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *                 { 
 *                   const SymmetricTensor<2, dim> &old_stress = 
 *                     local_quadrature_points_data[q_point].old_stress; 
 * 
 *                   cell_rhs(i) += 
 *                     (body_force_values[q_point](component_i) * 
 *                        fe_values.shape_value(i, q_point) - 
 *                      old_stress * get_strain(fe_values, i, q_point)) * 
 *                     fe_values.JxW(q_point); 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 现在我们有了对线性系统的局部贡献，我们需要将其转移到全局对象中。这与  step-17  中的做法完全相同。
 * 

 * 
 * 
 * @code
 *           cell->get_dof_indices(local_dof_indices); 
 * 
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix, 
 *                                                               cell_rhs, 
 *                                                               local_dof_indices, 
 *                                                               system_matrix, 
 *                                                               system_rhs); 
 *         } 
 * 
 * @endcode
 * 
 * 现在压缩矢量和系统矩阵。
 * 

 * 
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 * 
 * @endcode
 * 
 * 最后一步是再次修复边界值，就像我们在以前的程序中已经做的那样。一个稍微复杂的问题是， <code>apply_boundary_values</code> 函数希望有一个与矩阵和右手边兼容的解向量（即这里是一个分布式的%并行向量，而不是我们在这个程序中使用的顺序向量），以便用正确的边界值预设解向量的条目。我们以临时向量的形式提供这样一个兼容向量，然后将其复制到顺序向量中。
 * 

 * 
 * 我们通过展示边界值的灵活使用来弥补这种复杂性：按照我们创建三角形的方式，有三个不同的边界指标用来描述领域，分别对应于底面和顶面，以及内/外表面。我们希望施加以下类型的边界条件。内外圆柱体表面没有外力，这一事实对应于自然（诺伊曼型）边界条件，我们不需要做任何事情。在底部，我们希望完全没有运动，对应于圆柱体在边界的这一部分被夹住或粘住。然而，在顶部，我们希望有一个规定的垂直向下的运动来压缩圆柱体；此外，我们只希望限制垂直运动，而不是水平运动--可以把这种情况看作是一块油性良好的板坐在圆柱体的顶部将其向下推：圆柱体的原子被迫向下移动，但它们可以自由地沿着板水平滑动。
 * 

 * 
 * 描述这种情况的方法如下：对于边界指标为零（底面）的边界，我们使用一个二维的零函数，代表在任何坐标方向都没有运动。对于指标1（顶面）的边界，我们使用 <code>IncrementalBoundaryValues</code> 类，但我们为 <code>VectorTools::interpolate_boundary_values</code> 函数指定一个额外的参数，表示它应该适用于哪些矢量分量；这是一个针对每个矢量分量的bools矢量，由于我们只想限制垂直运动，它只有最后一个分量的设置。
 * 

 * 
 * 
 * @code
 *     FEValuesExtractors::Scalar                z_component(dim - 1); 
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              Functions::ZeroFunction<dim>(dim), 
 *                                              boundary_values); 
 *     VectorTools::interpolate_boundary_values( 
 *       dof_handler, 
 *       1, 
 *       IncrementalBoundaryValues<dim>(present_time, present_timestep), 
 *       boundary_values, 
 *       fe.component_mask(z_component)); 
 * 
 *     PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 
 *     MatrixTools::apply_boundary_values( 
 *       boundary_values, system_matrix, tmp, system_rhs, false); 
 *     incremental_displacement = tmp; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_timestep"></a> 
 * <h4>TopLevel::solve_timestep</h4>
 * 

 * 
 * 下一个函数是控制一个时间段内必须发生的所有事情的函数。从函数名称上看，事情的顺序应该是相对不言自明的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::solve_timestep() 
 *   { 
 *     pcout << "    Assembling system..." << std::flush; 
 *     assemble_system(); 
 *     pcout << " norm of rhs is " << system_rhs.l2_norm() << std::endl; 
 * 
 *     const unsigned int n_iterations = solve_linear_problem(); 
 * 
 *     pcout << "    Solver converged in " << n_iterations << " iterations." 
 *           << std::endl; 
 * 
 *     pcout << "    Updating quadrature point data..." << std::flush; 
 *     update_quadrature_point_history(); 
 *     pcout << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelsolve_linear_problem"></a> 
 * <h4>TopLevel::solve_linear_problem</h4>
 * 

 * 
 * 再次求解线性系统的工作原理与之前基本相同。唯一不同的是，我们只想保留一份完整的本地解向量，而不是从PETSc的求解程序中得到的分布式向量。为此，我们为分布式向量声明一个本地临时变量，并用本地变量的内容对其进行初始化（记得 <code>apply_boundary_values</code> 中调用的 <code>assemble_system</code> 函数预设了该向量中边界节点的值），用它进行求解，并在函数结束时将其再次复制到我们声明为成员变量的完整本地向量中。然后，挂起的节点约束只分布在本地拷贝上，也就是说，在每个处理器上都是独立的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   unsigned int TopLevel<dim>::solve_linear_problem() 
 *   { 
 *     PETScWrappers::MPI::Vector distributed_incremental_displacement( 
 *       locally_owned_dofs, mpi_communicator); 
 *     distributed_incremental_displacement = incremental_displacement; 
 * 
 *     SolverControl solver_control(dof_handler.n_dofs(), 
 *                                  1e-16 * system_rhs.l2_norm()); 
 * 
 *     PETScWrappers::SolverCG cg(solver_control, mpi_communicator); 
 * 
 *     PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix); 
 * 
 *     cg.solve(system_matrix, 
 *              distributed_incremental_displacement, 
 *              system_rhs, 
 *              preconditioner); 
 * 
 *     incremental_displacement = distributed_incremental_displacement; 
 * 
 *     hanging_node_constraints.distribute(incremental_displacement); 
 * 
 *     return solver_control.last_step(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveloutput_results"></a> 
 * <h4>TopLevel::output_results</h4>
 * 

 * 
 * 这个函数生成.vtu格式的图形输出，正如介绍中所解释的。每个进程将只对其拥有的单元格进行工作，然后将结果写入自己的文件中。此外，处理器0将写下引用所有.vtu文件的记录文件。
 * 

 * 
 * 这个函数的关键部分是给 <code>DataOut</code> 类提供一种方法，使其只对当前进程拥有的单元格进行工作。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 * @endcode
 * 
 * 然后，
 * 就像在 step-17 中一样，定义求解变量的名称（这里是位移增量）并排队输出求解向量。请注意在下面的开关中，我们如何确保如果空间维度应该不被处理，我们抛出一个异常，说我们还没有实现这种情况（另一个防御性编程的案例）。
 * 

 * 
 * 
 * @code
 *     std::vector<std::string> solution_names; 
 *     switch (dim) 
 *       { 
 *         case 1: 
 *           solution_names.emplace_back("delta_x"); 
 *           break; 
 *         case 2: 
 *           solution_names.emplace_back("delta_x"); 
 *           solution_names.emplace_back("delta_y"); 
 *           break; 
 *         case 3: 
 *           solution_names.emplace_back("delta_x"); 
 *           solution_names.emplace_back("delta_y"); 
 *           solution_names.emplace_back("delta_z"); 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     data_out.add_data_vector(incremental_displacement, solution_names); 
 * 
 * @endcode
 * 
 * 接下来的事情是，我们想输出类似于我们在每个单元中存储的应力的平均规范。这看起来很复杂，因为在目前的处理器上，我们只在那些实际属于目前进程的单元格上存储正交点的应力。换句话说，我们似乎无法计算出所有单元的平均应力。然而，请记住，我们源自 <code>DataOut</code> 的类只迭代那些实际属于当前处理器的单元，也就是说，我们不必为所有其他单元计算任何东西，因为这些信息不会被触及。下面的小循环就是这样做的。我们将整个区块包围在一对大括号中，以确保迭代器变量不会在它们被使用的区块结束后仍然意外地可见。
 * 

 * 
 * 
 * @code
 *     Vector<double> norm_of_stress(triangulation.n_active_cells()); 
 *     { 
 * 
 * @endcode
 * 
 * 在所有的单元格上循环...
 * 

 * 
 * 
 * @code
 *       for (auto &cell : triangulation.active_cell_iterators()) 
 *         if (cell->is_locally_owned()) 
 *           { 
 * 
 * @endcode
 * 
 * 在这些单元上，将所有正交点的应力相加...
 * 

 * 
 * 
 * @code
 *             SymmetricTensor<2, dim> accumulated_stress; 
 *             for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
 *               accumulated_stress += 
 *                 reinterpret_cast<PointHistory<dim> *>(cell->user_pointer())[q] 
 *                   .old_stress; 
 * 
 * @endcode
 * 
 * ...然后把平均值的常数写到它们的目的地。
 * 

 * 
 * 
 * @code
 *             norm_of_stress(cell->active_cell_index()) = 
 *               (accumulated_stress / quadrature_formula.size()).norm(); 
 *           } 
 * 
 * @endcode
 * 
 * 在我们不感兴趣的单元格上，将向量中各自的值设置为一个假值（规范必须是正值，大的负值应该能吸引你的眼球），以确保如果我们的假设有误，即这些元素不会出现在输出文件中，我们会通过观察图形输出发现。
 * 

 * 
 * 
 * @code
 *         else 
 *           norm_of_stress(cell->active_cell_index()) = -1e+20; 
 *     } 
 * 
 * @endcode
 * 
 * 最后把这个向量也附在上面，以便进行输出处理。
 * 

 * 
 * 
 * @code
 *     data_out.add_data_vector(norm_of_stress, "norm_of_stress"); 
 * 
 * @endcode
 * 
 * 作为最后一个数据，如果这是一个并行作业，让我们也把域划分为与处理器相关的子域。这与 step-17 程序中的工作方式完全相同。
 * 

 * 
 * 
 * @code
 *     std::vector<types::subdomain_id> partition_int( 
 *       triangulation.n_active_cells()); 
 *     GridTools::get_subdomain_association(triangulation, partition_int); 
 *     const Vector<double> partitioning(partition_int.begin(), 
 *                                       partition_int.end()); 
 *     data_out.add_data_vector(partitioning, "partitioning"); 
 * 
 * @endcode
 * 
 * 最后，有了这些数据，我们可以指示deal.II对信息进行整合，并产生一些中间数据结构，其中包含所有这些解决方案和其他数据向量。
 * 

 * 
 * 
 * @code
 *     data_out.build_patches(); 
 * 
 * @endcode
 * 
 * 让我们调用一个函数，打开必要的输出文件，将我们生成的数据写入其中。该函数根据给定的目录名（第一个参数）和文件名基数（第二个参数）自动构建文件名。它通过由时间步数和 "片数 "产生的片断来增加所产生的字符串，"片数 "对应于整个域的一部分，可以由一个或多个子域组成。
 * 

 * 
 * 该函数还为Paraview写了一个记录文件（后缀为`.pvd`），描述了所有这些输出文件如何组合成这个单一时间步骤的数据。
 * 

 * 
 * 
 * @code
 *     const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
 *       "./", "solution", timestep_no, mpi_communicator, 4); 
 * 
 * @endcode
 * 
 * 记录文件必须只写一次，而不是由每个处理器来写，所以我们在0号处理器上做这个。
 * 

 * 
 * 
 * @code
 *     if (this_mpi_process == 0) 
 *       { 
 * 
 * @endcode
 * 
 * 最后，我们写入paraview记录，它引用了所有.pvtu文件和它们各自的时间。注意，变量times_and_names被声明为静态的，所以它将保留前几个时间段的条目。
 * 

 * 
 * 
 * @code
 *         static std::vector<std::pair<double, std::string>> times_and_names; 
 *         times_and_names.push_back( 
 *           std::pair<double, std::string>(present_time, pvtu_filename)); 
 *         std::ofstream pvd_output("solution.pvd"); 
 *         DataOutBase::write_pvd_record(pvd_output, times_and_names); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_initial_timestep"></a> 
 * <h4>TopLevel::do_initial_timestep</h4>
 * 

 * 
 * 这个函数和下一个函数分别处理第一个和下一个时间步骤的整体结构。第一个时间步骤的工作量稍大，因为我们要在连续细化的网格上多次计算，每次都从一个干净的状态开始。在这些计算的最后，我们每次都计算增量位移，我们使用最后得到的增量位移的结果来计算产生的应力更新并相应地移动网格。在这个新的网格上，我们再输出解决方案和任何我们认为重要的附加数据。
 * 

 * 
 * 所有这些都会穿插着产生输出到控制台，以更新屏幕上的人正在发生的事情。如同在 step-17 中一样，使用 <code>pcout</code> instead of <code>std::cout</code> 可以确保只有一个并行进程实际在向控制台写数据，而不需要在每个产生输出的地方明确地编码一个if语句。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::do_initial_timestep() 
 *   { 
 *     present_time += present_timestep; 
 *     ++timestep_no; 
 *     pcout << "Timestep " << timestep_no << " at time " << present_time 
 *           << std::endl; 
 * 
 *     for (unsigned int cycle = 0; cycle < 2; ++cycle) 
 *       { 
 *         pcout << "  Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           create_coarse_grid(); 
 *         else 
 *           refine_initial_grid(); 
 * 
 *         pcout << "    Number of active cells:       " 
 *               << triangulation.n_active_cells() << " (by partition:"; 
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p) 
 *           pcout << (p == 0 ? ' ' : '+') 
 *                 << (GridTools::count_cells_with_subdomain_association( 
 *                      triangulation, p)); 
 *         pcout << ")" << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (by partition:"; 
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p) 
 *           pcout << (p == 0 ? ' ' : '+') 
 *                 << (DoFTools::count_dofs_with_subdomain_association(dof_handler, 
 *                                                                     p)); 
 *         pcout << ")" << std::endl; 
 * 
 *         solve_timestep(); 
 *       } 
 * 
 *     move_mesh(); 
 *     output_results(); 
 * 
 *     pcout << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLeveldo_timestep"></a> 
 * <h4>TopLevel::do_timestep</h4>
 * 

 * 
 * 后续的时间步骤比较简单，鉴于上面对前一个函数的解释，可能不需要更多的文件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::do_timestep() 
 *   { 
 *     present_time += present_timestep; 
 *     ++timestep_no; 
 *     pcout << "Timestep " << timestep_no << " at time " << present_time 
 *           << std::endl; 
 *     if (present_time > end_time) 
 *       { 
 *         present_timestep -= (present_time - end_time); 
 *         present_time = end_time; 
 *       } 
 * 
 *     solve_timestep(); 
 * 
 *     move_mesh(); 
 *     output_results(); 
 * 
 *     pcout << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TopLevelrefine_initial_grid"></a> 
 * <h4>TopLevel::refine_initial_grid</h4>
 * 

 * 
 * 当在连续细化的网格上求解第一个时间步骤时，调用以下函数。每次迭代后，它都会计算一个细化准则，细化网格，并将每个正交点的历史变量再次设置为干净状态。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::refine_initial_grid() 
 *   { 
 * 
 * @endcode
 * 
 * 首先，让每个进程计算其拥有的单元格的误差指标。
 * 

 * 
 * 
 * @code
 *     Vector<float> error_per_cell(triangulation.n_active_cells()); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(fe.degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       incremental_displacement, 
 *       error_per_cell, 
 *       ComponentMask(), 
 *       nullptr, 
 *       MultithreadInfo::n_threads(), 
 *       this_mpi_process); 
 * 
 * @endcode
 * 
 * 然后建立一个全局向量，我们将来自每个%并行进程的局部指标合并到其中。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_local_cells = 
 *       triangulation.n_locally_owned_active_cells(); 
 * 
 *     PETScWrappers::MPI::Vector distributed_error_per_cell( 
 *       mpi_communicator, triangulation.n_active_cells(), n_local_cells); 
 * 
 *     for (unsigned int i = 0; i < error_per_cell.size(); ++i) 
 *       if (error_per_cell(i) != 0) 
 *         distributed_error_per_cell(i) = error_per_cell(i); 
 *     distributed_error_per_cell.compress(VectorOperation::insert); 
 * 
 * @endcode
 * 
 * 一旦我们有了这个，就把它复制回所有处理器上的本地副本，并相应地完善网格。
 * 

 * 
 * 
 * @code
 *     error_per_cell = distributed_error_per_cell; 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     error_per_cell, 
 *                                                     0.35, 
 *                                                     0.03); 
 *     triangulation.execute_coarsening_and_refinement(); 
 * 
 * @endcode
 * 
 * 最后，在新的网格上再次设置正交点数据，并且只在那些我们已经确定是我们的单元上设置。
 * 

 * 
 * 
 * @code
 *     setup_quadrature_point_history(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelmove_mesh"></a> 
 * <h4>TopLevel::move_mesh</h4>
 * 

 * 
 * 在每个时间步骤结束时，我们根据这个时间步骤计算的增量位移来移动网格的节点。为了做到这一点，我们保留一个标志的向量，为每个顶点指示我们是否已经移动过它，然后在所有单元中循环，移动那些尚未移动的单元顶点。值得注意的是，我们从某个顶点相邻的单元中移动这个顶点并不重要：因为我们使用连续有限元计算位移，位移场也是连续的，我们可以从每个相邻的单元中计算某个顶点的位移。我们只需要确保每个节点都精确地移动一次，这就是为什么我们要保留标志的矢量。
 * 

 * 
 * 在这个函数中，有两个值得注意的地方。首先，我们如何使用 <code>cell-@>vertex_dof_index(v,d)</code> 函数获得给定顶点的位移场，该函数返回给定单元的 <code>d</code>th degree of freedom at vertex <code>v</code> 的索引。在本例中，k-th坐标方向的位移对应于有限元的k-th分量。使用这样的函数有一定的风险，因为它使用了我们在 <code>FESystem</code> 元素中为这个程序共同采取的元素顺序的知识。如果我们决定增加一个额外的变量，例如用于稳定的压力变量，并碰巧将其作为元素的第一个变量插入，那么下面的计算将开始产生无意义的结果。此外，这种计算还依赖于其他假设：首先，我们使用的元素确实有与顶点相关的自由度。对于目前的Q1元素来说确实如此，对于所有多项式阶的Qp元素来说也是如此  <code>p</code>  。然而，这对不连续的元素或混合公式的元素来说是不成立的。其次，它还建立在这样的假设上：一个顶点的位移只由与这个顶点相关的自由度的值决定；换句话说，所有对应于其他自由度的形状函数在这个特定的顶点是零。同样，对于目前的元素来说是这样的，但对于目前在deal.II中的所有元素来说并非如此。尽管有风险，我们还是选择使用这种方式，以便提出一种查询与顶点相关的单个自由度的方法。
 * 

 * 
 * 在这种情况下，指出一种更普遍的方法是很有意义的。对于一般的有限元来说，应该采用正交公式，将正交点放在单元的顶点上。梯形规则的 <code>QTrapezoid</code> 公式正是这样做的。有了这个正交公式，我们就可以在每个单元格中初始化一个 <code>FEValues</code> 对象，并使用 <code>FEValues::get_function_values</code> 函数来获得正交点，即单元格顶点的解函数值。这些是我们真正需要的唯一数值，也就是说，我们对与这个特定正交公式相关的权重（或 <code>JxW</code> 值）完全不感兴趣，这可以作为 <code>FEValues</code> 构造器的最后一个参数来指定。这个方案中唯一的一点小麻烦是，我们必须弄清楚哪个正交点对应于我们目前考虑的顶点，因为它们可能是以相同的顺序排列，也可能不是。
 * 

 * 
 * 如果有限元在顶点上有支持点（这里的支持点是有的；关于支持点的概念，见 @ref GlossSupport "支持点"），这种不便就可以避免了。对于这种情况，我们可以使用 FiniteElement::get_unit_support_points(). 构建一个自定义的正交规则，然后第一个 <code>cell-&gt;n_vertices()*fe.dofs_per_vertex</code> 正交点将对应于单元格的顶点，其顺序与 <code>cell-@>vertex(i)</code> 一致，同时考虑到矢量元素的支持点将被重复 <code>fe.dofs_per_vertex</code> 次。
 * 

 * 
 * 关于这个短函数值得解释的另一点是三角形类输出其顶点信息的方式：通过 <code>Triangulation::n_vertices</code> 函数，它公布了三角形中有多少个顶点。并非所有的顶点都是一直在使用的--有些是之前被粗化的单元的遗留物，自从deal.II以来一直存在，一旦一个顶点出现，即使数量较少的顶点消失了，也不会改变它的编号。其次， <code>cell-@>vertex(v)</code> 返回的位置不仅是一个类型为 <code>Point@<dim@></code> 的只读对象，而且事实上是一个可以写入的引用。这允许相对容易地移动网格的节点，但值得指出的是，使用该功能的应用程序有责任确保所得到的单元仍然有用，即没有扭曲到单元退化的程度（例如，用负的雅各布系数表示）。请注意，我们在这个函数中没有任何规定来实际保证这一点，我们只是有信心。
 * 

 * 
 * 在这个冗长的介绍之后，下面是全部20行左右的代码。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::move_mesh() 
 *   { 
 *     pcout << "    Moving mesh..." << std::endl; 
 * 
 *     std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 
 *     for (auto &cell : dof_handler.active_cell_iterators()) 
 *       for (const auto v : cell->vertex_indices()) 
 *         if (vertex_touched[cell->vertex_index(v)] == false) 
 *           { 
 *             vertex_touched[cell->vertex_index(v)] = true; 
 * 
 *             Point<dim> vertex_displacement; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               vertex_displacement[d] = 
 *                 incremental_displacement(cell->vertex_dof_index(v, d)); 
 * 
 *             cell->vertex(v) += vertex_displacement; 
 *           } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TopLevelsetup_quadrature_point_history"></a> 
 * <h4>TopLevel::setup_quadrature_point_history</h4>
 * 

 * 
 * 在计算的开始，我们需要设置历史变量的初始值，例如材料中的现有应力，我们将其存储在每个正交点中。如上所述，我们使用每个单元中都有的 <code>user_pointer</code> 来做这个。
 * 

 * 
 * 为了从更大的角度看这个问题，我们注意到，如果我们的模型中有先前可用的应力（为了这个程序的目的，我们假定这些应力不存在），那么我们就需要将先前存在的应力场插值到正交点上。同样，如果我们要模拟具有硬化/软化的弹塑性材料，那么我们就必须在每个正交点存储额外的历史变量，如累积塑性应变的当前屈服应力。预先存在的硬化或弱化也将通过在当前函数中插值这些变量来实现。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::setup_quadrature_point_history() 
 *   { 
 * 
 * @endcode
 * 
 * 为了慎重起见，我们把所有单元格的用户指针，不管是不是我们的，都设置为空指针。这样，如果我们访问了不应该访问的单元格的用户指针，一个分段故障将让我们知道这不应该发生。
 * 

 * 
 * 
 * @code
 *     triangulation.clear_user_data(); 
 * 
 * @endcode
 * 
 * 接下来，分配属于这个处理器职责范围内的正交对象。当然，这等于属于这个处理器的单元格的数量乘以我们的正交公式在每个单元格上的正交点的数量。由于`resize()`函数在要求的新大小小于旧大小的情况下，实际上并没有缩小分配的内存量，所以我们采用了一个技巧，首先释放所有的内存，然后再重新分配：我们声明一个空向量作为临时变量，然后交换旧向量和这个临时变量的内容。这就确保了`正交点历史'现在确实是空的，我们可以让现在保存着以前的向量内容的临时变量超出范围并被销毁。在下一步中，我们可以根据需要重新分配尽可能多的元素，矢量默认初始化`PointHistory`对象，这包括将压力变量设置为零。
 * 

 * 
 * 
 * @code
 *     { 
 *       std::vector<PointHistory<dim>> tmp; 
 *       quadrature_point_history.swap(tmp); 
 *     } 
 *     quadrature_point_history.resize( 
 *       triangulation.n_locally_owned_active_cells() * quadrature_formula.size()); 
 * 
 * @endcode
 * 
 * 最后再次循环所有单元，并将属于本处理器的单元的用户指针设置为指向此类对象的向量中与本单元对应的第一个正交点对象。
 * 

 * 
 * 
 * @code
 *     unsigned int history_index = 0; 
 *     for (auto &cell : triangulation.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           cell->set_user_pointer(&quadrature_point_history[history_index]); 
 *           history_index += quadrature_formula.size(); 
 *         } 
 * 
 * @endcode
 * 
 * 最后，为了慎重起见，确保我们对元素的计数是正确的，而且我们已经用完了之前分配的所有对象，并且没有指向任何超出向量末端的对象。这样的防御性编程策略总是很好的检查，以避免意外的错误，并防止将来对这个函数的修改忘记同时更新一个变量的所有用途。回顾一下，使用 <code>Assert</code> 宏的构造在优化模式下被优化掉了，所以不影响优化运行的运行时间。
 * 

 * 
 * 
 * @code
 *     Assert(history_index == quadrature_point_history.size(), 
 *            ExcInternalError()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TopLevelupdate_quadrature_point_history"></a> 
 * <h4>TopLevel::update_quadrature_point_history</h4>
 * 

 * 
 * 在每个时间步骤结束时，我们应该计算出一个增量的位移更新，使材料在其新的配置中能够容纳这个时间步骤中施加的外部体和边界力减去通过预先存在的内部应力施加的力之间的差异。为了在下一个时间步骤中获得预先存在的应力，我们必须用本时间步骤中计算的增量位移引起的应力来更新预先存在的应力。理想情况下，所产生的内应力之和将完全抵消所有的外力。事实上，一个简单的实验可以确保这一点：如果我们选择边界条件和体力与时间无关，那么强迫项（外力和内应力之和）应该正好是零。如果你做了这个实验，你会从每个时间步长的右手边的规范输出中意识到这几乎是事实：它并不完全是零，因为在第一个时间步长中，增量位移和应力的更新是相对于未变形的网格计算的，然后再进行变形。在第二个时间步骤中，我们再次计算位移和应力的更新，但这次是在变形的网格中 -- 在那里，结果的更新非常小但不完全是零。这可以迭代，在每一次迭代中，残差，即右手边向量的法线，都会减少；如果做这个小实验，就会发现这个残差的法线会随着迭代次数的增加而呈指数下降，在最初的快速下降之后，每次迭代大约会减少3.5倍（对于我看的一个测试案例，其他测试案例和其他未知数都会改变这个系数，但不会改变指数下降的情况）。
 * 

 * 
 * 在某种意义上，这可以被认为是一个准时序方案，以解决在一个以拉格朗日方式移动的网格上解决大变形弹性的非线性问题。
 * 

 * 
 * 另一个复杂的问题是，现有的（旧的）应力是在旧的网格上定义的，我们将在更新应力后移动这个网格。如果这个网格的更新涉及到单元的旋转，那么我们也需要对更新的应力进行旋转，因为它是相对于旧单元的坐标系计算的。
 * 

 * 
 * 因此，我们需要的是：在当前处理器拥有的每个单元上，我们需要从每个正交点存储的数据中提取旧的应力，计算应力更新，将两者相加，然后将结果与从当前正交点的增量位移计算出来的增量旋转一起旋转。下面我们将详细介绍这些步骤。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TopLevel<dim>::update_quadrature_point_history() 
 *   { 
 * 
 * @endcode
 * 
 * 首先，建立一个 <code>FEValues</code> 对象，我们将通过它来评估正交点的增量位移及其梯度，还有一个保存这些信息的向量。
 * 

 * 
 * 
 * @code
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients); 
 * 
 *     std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads( 
 *       quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim)); 
 * 
 * @endcode
 * 
 * 然后在所有单元格上循环，在属于我们子域的单元格中进行工作。
 * 

 * 
 * 
 * @code
 *     for (auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 * 
 * @endcode
 * 
 * 接下来，获得一个指向当前单元本地正交点历史数据的指针，作为防御措施，确保这个指针在全局数组的范围内。
 * 

 * 
 * 
 * @code
 *           PointHistory<dim> *local_quadrature_points_history = 
 *             reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 
 *           Assert(local_quadrature_points_history >= 
 *                    &quadrature_point_history.front(), 
 *                  ExcInternalError()); 
 *           Assert(local_quadrature_points_history <= 
 *                    &quadrature_point_history.back(), 
 *                  ExcInternalError()); 
 * 
 * @endcode
 * 
 * 然后在本单元上初始化 <code>FEValues</code> 对象，并提取正交点上的位移梯度，以便以后计算应变。
 * 

 * 
 * 
 * @code
 *           fe_values.reinit(cell); 
 *           fe_values.get_function_gradients(incremental_displacement, 
 *                                            displacement_increment_grads); 
 * 
 * @endcode
 * 
 * 然后在这个单元的正交点上循环。
 * 

 * 
 * 
 * @code
 *           for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
 *             { 
 * 
 * @endcode
 * 
 * 在每个正交点上，从梯度中计算出应变增量，并将其乘以应力-应变张量，得到应力更新。然后将此更新添加到该点已有的应变中。
 * 

 * 
 * 
 * @code
 *               const SymmetricTensor<2, dim> new_stress = 
 *                 (local_quadrature_points_history[q].old_stress + 
 *                  (stress_strain_tensor * 
 *                   get_strain(displacement_increment_grads[q]))); 
 * 
 * @endcode
 * 
 * 最后，我们要对结果进行旋转。为此，我们首先要从增量位移中计算出目前正交点的旋转矩阵。事实上，它可以从梯度中计算出来，而且我们已经有一个函数用于这个目的。
 * 

 * 
 * 
 * @code
 *               const Tensor<2, dim> rotation = 
 *                 get_rotation_matrix(displacement_increment_grads[q]); 
 * 
 * @endcode
 * 
 * 注意这个结果，即旋转矩阵，一般来说是一个等级为2的反对称张量，所以我们必须把它作为一个完整的张量来存储。
 * 

 * 
 * 有了这个旋转矩阵，在我们将对称张量 <code>new_stress</code> 扩展为全张量之后，我们可以通过从左和右的收缩来计算旋转的张量。
 * 

 * 
 * 
 * @code
 *               const SymmetricTensor<2, dim> rotated_new_stress = 
 *                 symmetrize(transpose(rotation) * 
 *                            static_cast<Tensor<2, dim>>(new_stress) * rotation); 
 * 
 * @endcode
 * 
 * 注意，虽然这三个矩阵的乘法结果应该是对称的，但由于浮点舍入的原因，它并不是对称的：我们得到的结果的非对角线元素有1e-16的不对称性。当把结果赋给一个 <code>SymmetricTensor</code> 时，该类的构造函数会检查对称性并意识到它不是完全对称的；然后它会引发一个异常。为了避免这种情况，我们明确地对结果进行对称，使其完全对称。
 * 

 * 
 * 所有这些操作的结果会被写回到原来的地方。
 * 

 * 
 * 
 * @code
 *               local_quadrature_points_history[q].old_stress = 
 *                 rotated_new_stress; 
 *             } 
 *         } 
 *   } 
 * 
 * @endcode
 * 
 * 这就结束了项目特定的命名空间  <code>Step18</code>  。其余的和往常一样，并且在  step-17  中已经显示：一个  <code>main()</code>  函数初始化和终止 PETSc，调用做实际工作的类，并确保我们捕捉所有传播到这一点的异常。
 * 

 * 
 * 
 * @code
 * } // namespace Step18 
 * 
 * int main(int argc, char **argv) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step18; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       TopLevel<3> elastic_problem; 
 *       elastic_problem.run(); 
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
examples/step-18/doc/results.dox



<a name="Results"></a><h1>Results</h1>



如果使用调试模式，运行该程序需要很长时间；在我的i7台式机上需要大约11分钟。幸运的是，经过优化编译的版本要快得多；在同一台机器上用<tt>make release</tt>命令重新编译后，程序只需要大约1.5分钟，这个时间要合理得多。


如果运行，该程序会打印出以下输出，解释它在这段时间内做了什么。

@verbatim
\$ time make run
[ 66%] Built target \step-18
[100%] Run \step-18 with Release configuration
Timestep 1 at time 1
  Cycle 0:
    Number of active cells:       3712 (by partition: 3712)
    Number of degrees of freedom: 17226 (by partition: 17226)
    Assembling system... norm of rhs is 1.88062e+10
    Solver converged in 103 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       12812 (by partition: 12812)
    Number of degrees of freedom: 51738 (by partition: 51738)
    Assembling system... norm of rhs is 1.86145e+10
    Solver converged in 121 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 2 at time 2
    Assembling system... norm of rhs is 1.84169e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 3 at time 3
    Assembling system... norm of rhs is 1.82355e+10
    Solver converged in 122 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 4 at time 4
    Assembling system... norm of rhs is 1.80728e+10
    Solver converged in 117 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 5 at time 5
    Assembling system... norm of rhs is 1.79318e+10
    Solver converged in 116 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 6 at time 6
    Assembling system... norm of rhs is 1.78171e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 7 at time 7
    Assembling system... norm of rhs is 1.7737e+10
    Solver converged in 112 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 8 at time 8
    Assembling system... norm of rhs is 1.77127e+10
    Solver converged in 111 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 9 at time 9
    Assembling system... norm of rhs is 1.78207e+10
    Solver converged in 113 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 10 at time 10
    Assembling system... norm of rhs is 1.83544e+10
    Solver converged in 115 iterations.
    Updating quadrature point data...
    Moving mesh...


[100%] Built target run
make run  176.82s user 0.15s system 198% cpu 1:28.94 total
@endverbatim

换句话说，它是在12,000个单元和大约52,000个未知数的情况下进行计算。不是很多，但对于一个耦合的三维问题来说，足以让计算机忙上一阵子。在一天结束的时候，这就是我们的输出。

@verbatim
\$ ls -l *vtu *visit


-rw-r--r-- 1 drwells users 1706059 Feb 13 19:36 solution-0010.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0010.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0010.visit


-rw-r--r-- 1 drwells users 1707907 Feb 13 19:36 solution-0009.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:36 solution-0009.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:36 solution-0009.visit


-rw-r--r-- 1 drwells users 1703771 Feb 13 19:35 solution-0008.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0008.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0008.visit


-rw-r--r-- 1 drwells users 1693671 Feb 13 19:35 solution-0007.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0007.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0007.visit


-rw-r--r-- 1 drwells users 1681847 Feb 13 19:35 solution-0006.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0006.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0006.visit


-rw-r--r-- 1 drwells users 1670115 Feb 13 19:35 solution-0005.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0005.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0005.visit


-rw-r--r-- 1 drwells users 1658559 Feb 13 19:35 solution-0004.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0004.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0004.visit


-rw-r--r-- 1 drwells users 1639983 Feb 13 19:35 solution-0003.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0003.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0003.visit


-rw-r--r-- 1 drwells users 1625851 Feb 13 19:35 solution-0002.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:35 solution-0002.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:35 solution-0002.visit


-rw-r--r-- 1 drwells users 1616035 Feb 13 19:34 solution-0001.000.vtu


-rw-r--r-- 1 drwells users     761 Feb 13 19:34 solution-0001.pvtu


-rw-r--r-- 1 drwells users      33 Feb 13 19:34 solution-0001.visit
@endverbatim




如果我们用VisIt或Paraview将这些文件可视化，我们就能看到我们的强制压缩对圆柱体造成的灾难的全貌（图像中的颜色编码了材料中的应力规范）。


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0002.0000.png" alt="Time = 2" width="400"> </div> <div class="text" align="center"> Time = 2 </div> <div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0005.0000.png" alt="时间=5" width="400"> </div> <div class="text" align="center"> 时间=5 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0007.0000.png" alt="时间=7" width="400"> </div> <div class="text" align="center">时间=7 </div> </div> </div>


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0008.0000.png" alt="Time = 8" width="400"> </div> <div class="text" align="center"> 时间 = 8 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0009.0000.png" alt="时间=9" width="400"> </div> <div class="text" align="center"> 时间=9 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.sequential-0010.0000.png" alt="Time = 10" width="400"> </div> <div class="text" align="center"> Time = 10 </div> </div> </div>


可以清楚地看到，当我们不断压缩圆柱体时，它开始在完全约束的底面附近弯曲，并在大约8个时间单位后，以方位对称的方式弯曲。


虽然这个结果对于对称几何和加载来说似乎是合理的，但计算是否完全收敛还有待确定。为了确定是否收敛，我们再次运行程序，在开始时再进行一次全局细化，并将时间步长减半。这在单机上会花费很长的时间，所以我们使用了一个合适的工作站，在16个处理器上并行运行。现在输出的开头看起来像这样。

@verbatim
Timestep 1 at time 0.5
  Cycle 0:
    Number of active cells:       29696 (by partition: 1808+1802+1894+1881+1870+1840+1884+1810+1876+1818+1870+1884+1854+1903+1816+1886)
    Number of degrees of freedom: 113100 (by partition: 6936+6930+7305+7116+7326+6869+7331+6786+7193+6829+7093+7162+6920+7280+6843+7181)
    Assembling system... norm of rhs is 1.10765e+10
    Solver converged in 209 iterations.
    Updating quadrature point data...
  Cycle 1:
    Number of active cells:       102034 (by partition: 6387+6202+6421+6341+6408+6201+6428+6428+6385+6294+6506+6244+6417+6527+6299+6546)
    Number of degrees of freedom: 359337 (by partition: 23255+21308+24774+24019+22304+21415+22430+22184+22298+21796+22396+21592+22325+22553+21977+22711)
    Assembling system... norm of rhs is 1.35759e+10
    Solver converged in 268 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 2 at time 1
    Assembling system... norm of rhs is 1.34674e+10
    Solver converged in 267 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 3 at time 1.5
    Assembling system... norm of rhs is 1.33607e+10
    Solver converged in 265 iterations.
    Updating quadrature point data...
    Moving mesh...


Timestep 4 at time 2
    Assembling system... norm of rhs is 1.32558e+10
    Solver converged in 263 iterations.
    Updating quadrature point data...
    Moving mesh...


[...]


Timestep 20 at time 10
    Assembling system... norm of rhs is 1.47755e+10
    Solver converged in 425 iterations.
    Updating quadrature point data...
    Moving mesh...
@endverbatim

考虑到我们是在三维空间中，这是一个相当好的未知数的数量。这个程序的输出是每个时间步骤的16个文件。

@verbatim
\$ ls -l solution-0001*


-rw-r--r-- 1 wellsd2 user 761065 Feb 13 21:09 solution-0001.000.vtu


-rw-r--r-- 1 wellsd2 user 759277 Feb 13 21:09 solution-0001.001.vtu


-rw-r--r-- 1 wellsd2 user 761217 Feb 13 21:09 solution-0001.002.vtu


-rw-r--r-- 1 wellsd2 user 761605 Feb 13 21:09 solution-0001.003.vtu


-rw-r--r-- 1 wellsd2 user 756917 Feb 13 21:09 solution-0001.004.vtu


-rw-r--r-- 1 wellsd2 user 752669 Feb 13 21:09 solution-0001.005.vtu


-rw-r--r-- 1 wellsd2 user 735217 Feb 13 21:09 solution-0001.006.vtu


-rw-r--r-- 1 wellsd2 user 750065 Feb 13 21:09 solution-0001.007.vtu


-rw-r--r-- 1 wellsd2 user 760273 Feb 13 21:09 solution-0001.008.vtu


-rw-r--r-- 1 wellsd2 user 777265 Feb 13 21:09 solution-0001.009.vtu


-rw-r--r-- 1 wellsd2 user 772469 Feb 13 21:09 solution-0001.010.vtu


-rw-r--r-- 1 wellsd2 user 760833 Feb 13 21:09 solution-0001.011.vtu


-rw-r--r-- 1 wellsd2 user 782241 Feb 13 21:09 solution-0001.012.vtu


-rw-r--r-- 1 wellsd2 user 748905 Feb 13 21:09 solution-0001.013.vtu


-rw-r--r-- 1 wellsd2 user 738413 Feb 13 21:09 solution-0001.014.vtu


-rw-r--r-- 1 wellsd2 user 762133 Feb 13 21:09 solution-0001.015.vtu


-rw-r--r-- 1 wellsd2 user   1421 Feb 13 21:09 solution-0001.pvtu


-rw-r--r-- 1 wellsd2 user    364 Feb 13 21:09 solution-0001.visit
@endverbatim




这里首先是我们计算的网格，以及16个处理器的分区。


<div class="twocolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-000mesh.png" alt="Discretization" width="400"> </div> <div class="text" align="center"> Discretization </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.p.png" alt="Parallel partitioning" width="400"> </div> <div class="text" align="center"> Parallel partitioning</div> </div> </div>


最后，这里是与我们之前展示的更小的顺序情况相同的输出。

<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0002.s.png" alt="Time = 2" width="400"> </div> <div class="text" align="center"> 时间 = 2 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0005.s.png" alt="时间=5" width="400"> </div> <div class="text" align="center"> 时间=5 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0007.s.png" alt="Time = 7" width="400"> </div> <div class="text" align="center"> Time = 7 </div> </div> </div>


<div class="threecolumn" style="width: 80%"> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0008.s.png" alt="Time = 8" width="400"> </div> <div class="text" align="center"> 时间 = 8 </div> </div> <div class="parent"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0009.s.png" alt="Time = 9" width="400"> </div> <div class="text" align="center"> Time = 9 </div> </div> <div class="father"> <div class="img" align="center"> <img src="https://www.dealii.org/images/steps/developer/step-18.parallel-0010.s.png" alt="Time = 10" width="400"> </div> <div class="text" align="center"> Time = 10 </div> </div> </div>


和以前一样，我们观察到，在高轴向压缩时，圆柱体开始弯曲，但这一次最终是在自己身上塌陷。与我们的第一次运行相反，在模拟结束时，变形模式变得不对称（中心隆起向侧面偏转）。该模型显然没有规定这一点（我们所有的力和边界偏转都是对称的），但这种效果可能在物理上是正确的：在现实中，身体材料属性的小不均匀性会导致它向一侧弯曲以逃避强制力；在数值模拟中，小的扰动，如数值舍入或迭代求解器对线性系统的不精确求解，也会产生同样的效果。在自适应计算中，另一个典型的不对称来源是每一步只细化一定的单元，这可能导致不对称的网格，即使原来的粗网格是对称的。


如果将其与之前的运行相比较，结果在质和量上都有不同。因此，以前的计算肯定没有收敛，尽管我们不能肯定地说现在的计算有什么问题。我们需要一个更精细的计算来找出答案。然而，这一点可能是没有意义的：详细看一下最后一张图片，很明显，不仅我们选择的线性小变形模型是完全不够的，而且对于一个现实的模拟，我们还需要确保身体在变形过程中不相交（如果我们继续压缩圆柱体，我们会观察到一些自我相交）。如果没有这样的表述，我们就不能指望任何东西都有物理意义，即使它能产生漂亮的图片!




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这个程序并没有真正解决一个在实践中有很多应用的方程：基于纯弹性规律的准静态材料变形几乎是无聊的。然而，该程序可以作为更有趣的实验的起点，而这确实是编写该程序的最初动机。这里有一些建议，说明这个程序缺少什么，以及它可以在什么方向上进行扩展。

<a name="Plasticitymodels"></a><h5>Plasticity models</h5>


最明显的扩展是使用一个更现实的材料模型来处理大规模的静态变形。这方面的自然选择是塑性，其中应力和应变之间的非线性关系取代了方程<a href="#step_18.stress-strain">[stress-strain]</a>。塑性模型的编程通常相当复杂，因为应力-应变关系通常是非平滑的。可以认为材料只能承受一个最大的应力（屈服应力），之后它就会开始&ldquo;流动&rdquo;。这方面的数学描述可以以变分不等式的形式给出，也可以将其视为弹性能量的最小化

@f[
  E(\mathbf{u}) =
  (\varepsilon(\mathbf{u}), C\varepsilon(\mathbf{u}))_{\Omega}


  - (\mathbf{f}, \mathbf{u})_{\Omega} - (\mathbf{b}, \mathbf{u})_{\Gamma_N},


@f]

受制于约束条件

@f[
  f(\sigma(\mathbf{u})) \le 0


@f]

对应力的影响。这种扩展使得在每个时间步长中要解决的问题是非线性的，所以我们需要在每个时间步长中的另一个循环。

在不进一步了解这个模型的细节的情况下，我们可以参考Simo和Hughes关于&ldquo;计算非弹性&rdquo;的优秀书籍，以全面了解解决塑性模型的计算策略。另外，在S. Commend, A. Truty, and Th. Zimmermann的文章中，对塑性的算法做了简单而简洁的描述。Zimmermann;  @cite CTZ04  。




<a name="Stabilizationissues"></a><h5>Stabilization issues</h5>


我们选择的公式，即对位移矢量的所有分量使用分片（双，三）线性元素，并将应力视为依赖于位移的变量，对于大多数材料是合适的。然而，对于不可压缩或几乎不可压缩的材料，这种所谓的基于位移的公式变得不稳定，并表现出虚假的模式。虽然流体通常不是弹性的（在大多数情况下，应力取决于速度梯度，而不是位移梯度，但也有例外，如电流变流体），但也有少数固体是几乎不可压缩的，如橡胶。另一种情况是，许多塑性模型最终让材料变得不可压缩，尽管这不在本方案的范围之内。

不可压缩性是由泊松比来表征的

@f[
  \nu = \frac{\lambda}{2(\lambda+\mu)},


@f]

其中 $\lambda,\mu$ 是材料的Lam&eacute; 常数。物理约束表明 $-1\le \nu\le \frac 12$ （该条件也来自于数学稳定性考虑）。如果 $\nu$ 接近 $\frac 12$ ，则材料变得不可压缩。在这种情况下，纯粹的基于位移的公式不再适合于解决这类问题，必须采用稳定化技术以获得稳定和准确的解决方案。上面引用的书和论文给出了如何做到这一点的指示，但在这个问题上也有大量的文献；在H.-Y. Duan和Q. Lin的论文的参考文献中可以找到一个获得该主题概述的良好开端。H.-Y. Duan and Q. Lin;  @cite DL05  。




<a name="Refinementduringtimesteps"></a><h5>Refinement during timesteps</h5>


在目前的形式下，程序只对初始网格进行若干次细化，然后就不再进行细化。对于任何一种现实的模拟，我们都希望将其扩展到每隔几步就对网格进行细化和粗化。事实上，这并不难做到，但如果你愿意的话，可以留待将来的教程程序或作为练习。

我们必须克服的主要复杂问题是，我们必须将存储在旧网格单元的正交点中的数据转移到新网格中，最好是通过某种投影方案。这方面的一般方法是这样的。

- 开始时，数据只在各个单元的正交点上可用，而不是作为一个到处定义的有限元场。

- 所以让我们找到一个<i>is</i>处处定义的有限元场，这样我们以后就可以把它插到新网格的正交点上。一般来说，要找到一个与正交点中的数值完全匹配的连续有限元场是很困难的，因为这些场的自由度数与正交点的数量不匹配，这个全局场的节点值要么是过定的，要么是欠定的。但是找到一个与正交点数值相匹配的不连续场通常不是很困难；例如，如果你有一个QGauss(2)正交公式（即2d中每个单元4个点，3d中8个点），那么就可以使用FE_DGQ(1)类型的有限元，即双/三线性函数，因为这些函数在2d中每个单元有4个自由度，在3d中有8个自由度。

- 有一些函数可以使这种从单个点到全局场的转换更简单。如果你使用QGauss(2)正交公式，下面这段伪代码应该会有所帮助。请注意，下面的投影矩阵的乘法需要一个标量分量的向量，也就是说，我们一次只能将一组标量从正交点转换成自由度，反之亦然。所以我们需要分别存储每个应力分量，这需要 <code>dim*dim</code> 个向量。我们将把这组向量存储在一个二维数组中，以便于用读出应力张量的方式来读出分量。   因此，我们将对每个单元的应力分量进行循环，并将这些值存储在全局历史域中。(前缀 <code>history_</code> 表示我们的工作与正交点中定义的历史变量有关。)   @code
    FE_DGQ<dim>     history_fe (1);
    DoFHandler<dim> history_dof_handler (triangulation);
    history_dof_handler.distribute_dofs (history_fe);


    std::vector< std::vector< Vector<double> > >
                 history_field (dim, std::vector< Vector<double> >(dim)),
                 local_history_values_at_qpoints (dim, std::vector< Vector<double> >(dim)),
                 local_history_fe_values (dim, std::vector< Vector<double> >(dim));


    for (unsigned int i=0; i<dim; i++)
      for (unsigned int j=0; j<dim; j++)
      {
        history_field[i][j].reinit(history_dof_handler.n_dofs());
        local_history_values_at_qpoints[i][j].reinit(quadrature.size());
        local_history_fe_values[i][j].reinit(history_fe.n_dofs_per_cell());
      }


    FullMatrix<double> qpoint_to_dof_matrix (history_fe.dofs_per_cell,
                                             quadrature.size());
    FETools::compute_projection_from_quadrature_points_matrix
              (history_fe,
               quadrature, quadrature,
               qpoint_to_dof_matrix);


    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();


    for (; cell!=endc; ++cell, ++dg_cell)
      {


        PointHistory<dim> *local_quadrature_points_history
          = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());


        Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
                ExcInternalError());
        Assert (local_quadrature_points_history < &quadrature_point_history.back(),
                ExcInternalError());


        for (unsigned int i=0; i<dim; i++)
          for (unsigned int j=0; j<dim; j++)
          {
            for (unsigned int q=0; q<quadrature.size(); ++q)
              local_history_values_at_qpoints[i][j](q)
                = local_quadrature_points_history[q].old_stress[i][j];


            qpoint_to_dof_matrix.vmult (local_history_fe_values[i][j],
                                        local_history_values_at_qpoints[i][j]);


            dg_cell->set_dof_values (local_history_fe_values[i][j],
                                     history_field[i][j]);
          }
      }
  @endcode



- 现在我们有了一个全局场，我们可以像往常一样使用SolutionTransfer类来细化网格并转移history_field向量。这将把所有的东西从旧的网格插值到新的网格。

- 在最后一步，我们必须将数据从现在插值的全局场返回到新网格上的正交点。下面的代码将做到这一点。   @code
    FullMatrix<double> dof_to_qpoint_matrix (quadrature.size(),
                                             history_fe.n_dofs_per_cell());
    FETools::compute_interpolation_to_quadrature_points_matrix
              (history_fe,
               quadrature,
               dof_to_qpoint_matrix);


    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end(),
                                                   dg_cell = history_dof_handler.begin_active();


    for (; cell != endc; ++cell, ++dg_cell)
    {
      PointHistory<dim> *local_quadrature_points_history
       = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());


      Assert (local_quadrature_points_history >= &quadrature_point_history.front(),
              ExcInternalError());
      Assert (local_quadrature_points_history < &quadrature_point_history.back(),
              ExcInternalError());


      for (unsigned int i=0; i<dim; i++)
        for (unsigned int j=0; j<dim; j++)
        {
          dg_cell->get_dof_values (history_field[i][j],
                                   local_history_fe_values[i][j]);


          dof_to_qpoint_matrix.vmult (local_history_values_at_qpoints[i][j],
                                      local_history_fe_values[i][j]);


          for (unsigned int q=0; q<quadrature.size(); ++q)
            local_quadrature_points_history[q].old_stress[i][j]
              = local_history_values_at_qpoints[i][j](q);
      }
  @endcode



一旦我们并行运行程序，情况就变得有点复杂了，因为那时每个进程只为它在旧网格上拥有的单元存储这些数据。也就是说，如果你在正交点转移到全局向量之后，使用 <code>history_field</code> 的并行向量就可以做到这一点。




<a name="Ensuringmeshregularity"></a><h5>Ensuring mesh regularity</h5>


目前，程序没有尝试确保一个单元在时间步数结束时移动其顶点后，仍然具有有效的几何形状（即它的雅各布行列式是正的，并且在任何地方都远离零的界限）。事实上，设置边界值和强迫项并不难，这样就可以很快得到扭曲和倒置的单元。当然，在某些大变形的情况下，这在有限网格的情况下是不可避免的，但在其他一些情况下，通过适当的网格细化和/或减少时间步长，这应该是可以避免的。这个程序没有做到这一点，但是一个更复杂的版本肯定应该采用某种启发式方法来定义哪些单元的变形量是可以接受的，哪些是不可以的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-18.cc"
*/
