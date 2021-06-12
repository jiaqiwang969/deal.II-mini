/**
@page step_28 The step-28 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Theeigenvalueproblem">The eigenvalue problem</a>
        <li><a href="#Meshesandmeshrefinement">Meshes and mesh refinement</a>
      <ul>
        <li><a href="#Meshrefinement">Mesh refinement</a>
        <li><a href="#Assemblingtermsondifferentmeshes">Assembling terms on different meshes</a>
      </ul>
        <li><a href="#Descriptionofthetestcase">Description of the test case</a>
        <li><a href="#Whattheprogramdoesandhowitdoesthat">What the program does (and how it does that)</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Materialdata">Material data</a>
        <li><a href="#ThecodeEnergyGroupcodeclass">The <code>EnergyGroup</code> class</a>
      <ul>
      <ul>
        <li><a href="#codeEnergyGroupcodepublicmemberfunctions"><code>EnergyGroup</code> public member functions</a>
        <li><a href="#codeEnergyGroupcodepublicdatamembers"><code>EnergyGroup</code> public data members</a>
        <li><a href="#codeEnergyGroupcodeprivatedatamembers"><code>EnergyGroup</code> private data members</a>
        <li><a href="#codeEnergyGroupcodeprivatememberfunctions"><code>EnergyGroup</code> private member functions</a>
      </ul>
        <li><a href="#ImplementationofthecodeEnergyGroupcodeclass">Implementation of the <code>EnergyGroup</code> class</a>
      <ul>
        <li><a href="#codeEnergyGroupsetup_linear_systemcode"><code>EnergyGroup::setup_linear_system</code></a>
        <li><a href="#codeEnergyGroupassemble_system_matrixcode"><code>EnergyGroup::assemble_system_matrix</code></a>
        <li><a href="#codeEnergyGroupassemble_ingroup_rhscode"><code>EnergyGroup::assemble_ingroup_rhs</code></a>
        <li><a href="#codeEnergyGroupassemble_cross_group_rhscode"><code>EnergyGroup::assemble_cross_group_rhs</code></a>
        <li><a href="#codeEnergyGroupassemble_cross_group_rhs_recursivecode"><code>EnergyGroup::assemble_cross_group_rhs_recursive</code></a>
        <li><a href="#codeEnergyGroupget_fission_sourcecode"><code>EnergyGroup::get_fission_source</code></a>
        <li><a href="#codeEnergyGroupsolvecode"><code>EnergyGroup::solve</code></a>
        <li><a href="#codeEnergyGroupestimate_errorscode"><code>EnergyGroup::estimate_errors</code></a>
        <li><a href="#codeEnergyGrouprefine_gridcode"><code>EnergyGroup::refine_grid</code></a>
        <li><a href="#codeEnergyGroupoutput_resultscode"><code>EnergyGroup::output_results</code></a>
      </ul>
      </ul>
        <li><a href="#ThecodeNeutronDiffusionProblemcodeclasstemplate">The <code>NeutronDiffusionProblem</code> class template</a>
      <ul>
      <ul>
        <li><a href="#codeNeutronDiffusionProblemcodeprivatememberfunctions"><code>NeutronDiffusionProblem</code> private member functions</a>
        <li><a href="#codeNeutronDiffusionProblemcodeprivatemembervariables"><code>NeutronDiffusionProblem</code> private member variables</a>
      </ul>
        <li><a href="#ImplementationofthecodeParameterscodeclass">Implementation of the <code>Parameters</code> class</a>
        <li><a href="#ImplementationofthecodeNeutronDiffusionProblemcodeclass">Implementation of the <code>NeutronDiffusionProblem</code> class</a>
      <ul>
        <li><a href="#codeNeutronDiffusionProbleminitialize_problemcode"><code>NeutronDiffusionProblem::initialize_problem</code></a>
        <li><a href="#codeNeutronDiffusionProblemget_total_fission_sourcecode"><code>NeutronDiffusionProblem::get_total_fission_source</code></a>
        <li><a href="#codeNeutronDiffusionProblemrefine_gridcode"><code>NeutronDiffusionProblem::refine_grid</code></a>
        <li><a href="#codeNeutronDiffusionProblemruncode"><code>NeutronDiffusionProblem::run</code></a>
      </ul>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-28/doc/intro.dox

 <br> 

<i>This program was contributed by Yaqi Wang and Wolfgang
Bangerth. Results from this program are used and discussed in the publication
"Three-dimensional h-adaptivity for the multigroup neutron diffusion
equations" by Yaqi Wang, Wolfgang Bangerth and Jean Ragusa. The paper's full
bibliographic details are as follows:


@code
@Article{WBR09,
  author  = {Yaqi Wang and Wolfgang Bangerth and Jean Ragusa},
  title   = {Three-dimensional h-adaptivity for the multigroup
             neutron diffusion equations},
  journal = {Progr. Nucl. Energy},
  year    = 2009,
  volume  = 51,
  pages   = {543--555}
}
@endcode@code
@Article{WBR09,
  author  = {Yaqi Wang and Wolfgang Bangerth and Jean Ragusa},
  title   = {Three-dimensional h-adaptivity for the multigroup
             neutron diffusion equations},
  journal = {Progr. Nucl. Energy},
  year    = 2009,
  volume  = 51,
  pages   = {543--555}
}
@endcode


The paper is available <a target="_top"
href="https://www.semanticscholar.org/paper/Three-dimensional-h-adaptivity-for-the-multigroup-Wang-Bangerth/900592e8e891d9b888d59a69ec58bf2bbda56b4b">here</a><a target="_top"
href="https://www.semanticscholar.org/paper/Three-dimensional-h-adaptivity-for-the-multigroup-Wang-Bangerth/900592e8e891d9b888d59a69ec58bf2bbda56b4b">here</a>.
</i>

 <br> 


<a name="Introduction"></a><a name="Intro"></a> <h1>Introduction</h1> 在这个例子中，我们打算解决中子传输方程的多组扩散近似。基本上，看待这个问题的方式如下。在核反应堆中，中子以不同的能量飞驰，被吸收或散射，或开始一个新的裂变事件。如果从足够长的长度尺度来看，中子的运动可以被视为一个扩散过程。


对这一点的数学描述将把中子归入能级，并考虑每个能级或能量组中的中子通量的平衡方程。然后，散射、吸收和裂变事件将成为描述中子通量的扩散方程中的算子。假设我们有能量组 $g=1,\ldots,G$ ，按照惯例，我们假设能量最高的中子在1组，能量最低的中子在 $G$ 组。那么每组的中子通量满足以下公式。

@f{eqnarray*}
\frac 1{v_g}\frac{\partial \phi_g(x,t)}{\partial t}
&=&
\nabla \cdot(D_g(x) \nabla \phi_g(x,t))


-
\Sigma_{r,g}(x)\phi_g(x,t)
\\
&& \qquad
+
\chi_g\sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}(x,t)
+
\sum_{g'\ne g}\Sigma_{s,g'\to g}(x)\phi_{g'}(x,t)
+
s_{\mathrm{ext},g}(x,t)


@f}

通过适当的边界条件增强。这里， $v_g$ 是 $g$ 组内中子的速度。换句话说， $g$ 组中的中子通量的时间变化受以下过程支配。   <ul>   <li>  扩散  $\nabla \cdot(D_g(x) \nabla \phi_g(x,t))$  。这里， $D_g$ 是（空间上可变的）扩散系数。   <li>  吸收  $\Sigma_{r,g}(x)\phi_g(x,t)$  （注意是负号）。系数 $\Sigma_{r,g}$ 被称为<i>removal
  cross section</i>。   <li>  核裂变  $\chi_g\sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}(x,t)$  。   能量的中子 $g$ 的产生与能量的中子通量 $g'$ 乘以能量的中子 $\Sigma_{f,g'}$ 引起裂变事件的概率 $g'$ 乘以每个裂变事件中产生的中子数量 $\nu$ 乘以该事件中产生的中子具有能量的概率 $g$  。   $\nu\Sigma_{f,g'}$ 被称为<i>fission cross section</i>， $\chi_g$ 被称为<i>fission spectrum</i>。我们将把 $\chi_g\nu\Sigma_{f,g'}$ 这个词表示为程序中的<i>fission distribution cross
    section</i>。   <li> 能量为 $\sum_{g'\ne g}\Sigma_{s,g'\to g}(x)\phi_{g'}(x,t)$ 的中子的散射 $g'$ 产生能量为 $g$ 的中子 。   $\Sigma_{s,g'\to g}$ 被称为<i>scattering cross
    section</i>。弹性、群内散射 $g'=g$ 的情况也存在，但我们将其归入清除截面。 $g'<g$ 的情况被称为向下散射，因为中子在这样的事件中失去了能量。另一方面， $g'>g$ 对应于上散射：中子在散射事件中从其周围原子的热运动中获得能量；因此，上散射仅对动能已经与热动能处于同一等级（即在 $eV$ 以下范围）的中子是一个重要过程。   <li>  一个外源  $s_{\mathrm{ext},g}$  。   </ul> 

为了在反应堆分析中进行现实模拟，人们可能希望将中子能量的连续谱分成许多能量组，通常多达100个。然而，如果对某些类型的反应堆（例如压水反应堆，PWR）的中子能谱有足够的了解，那么只用2个能群就有可能获得满意的结果。

在本教程的程序中，我们提供的结构可以根据需要用尽可能多的能量组进行计算。然而，为了保持适度的计算时间，并避免将数百个系数列表，我们只提供上述两组模拟方程的系数，即  $g=1,2$  。然而，我们确实考虑了一个现实的情况，即假设系数不是恒定的，而是取决于以相当复杂的方式装配到反应堆燃料组件的材料（见下文）。




<a name="Theeigenvalueproblem"></a><h3>The eigenvalue problem</h3>


如果我们一次考虑所有的能量组，我们可以将上述方程写成以下算子形式。

@f{eqnarray*}
\frac 1v \frac{\partial \phi}{\partial t}
=


-L\phi
+
F\phi
+
X\phi
+
s_{\mathrm{ext}},


@f}

其中 $L,F,X$ 分别是沉降、裂变和散射算子。  这里的 $L$ 包括扩散和清除项。请注意， $L$ 是对称的，而 $F$ 和 $X$ 不是。

众所周知，如果算子 $-L+F+X$ 的所有特征值都是负的，这个方程就有一个稳定的解。这可以通过将方程乘以 $\phi$ 并在域上进行积分而很容易看出，从而导致

@f{eqnarray*}
  \frac 1{2v} \frac{\partial}{\partial t}  \|\phi\|^2 = ((-L+F+X)\phi,\phi).


@f}

稳定性意味着解决方案不会增长，也就是说，我们希望左手边的值小于零，如果右边算子的特征值都是负的，就会出现这种情况。由于显而易见的原因，如果核反应堆产生的中子通量呈指数增长，这不是很理想，所以特征值分析是核工程师的面包和主食。因此，该程序的要点是考虑特征值问题

@f{eqnarray*}
  (L-F-X) \phi = \lambda \phi,


@f}

其中我们要确保所有的特征值都是正的。请注意， $L$ ，即扩散算子加上吸收（去除），是正定的；因此，所有特征值为正的条件意味着我们要确保裂变和组间散射足够弱，不会使频谱转移到负值。

在核工程中，人们通常会关注特征值问题的一个稍微不同的表述。为此，我们不只是与 $\phi$ 相乘并整合，而是与 $\phi(L-X)^{-1}$ 相乘。然后我们得到以下演化方程。

@f{eqnarray*}
  \frac 1{2v} \frac{\partial}{\partial t}  \|\phi\|^2_{(L-X)^{-1}} = ((L-X)^{-1}(-L+F+X)\phi,\phi).


@f}

如果以下问题的特征值都是负的，那么稳定性就得到了保证。

@f{eqnarray*}
  (L-X)^{-1}(-L+F+X)\phi = \lambda_F \phi,


@f}

等同于特征值问题

@f{eqnarray*}
  (L-X)\phi = \frac 1{\lambda_F+1} F \phi.


@f}

核工程中的典型表述是将其写为

@f{eqnarray*}
  (L-X) \phi = \frac 1{k_{\mathrm{eff}}} F \phi,


@f}

其中 $k_{\mathrm{eff}}=\frac 1{\lambda^F+1}$  。直观地说， $k_{\mathrm{eff}}$ 是类似于每个典型时间尺度的中子的乘法系数，对于反应堆的稳定运行，它应该小于或等于1：如果它小于1，链式反应将消亡，而核弹等的 $k$ -特征值大于1。一个稳定的反应堆应该有 $k_{\mathrm{eff}}=1$  。

对于那些想知道在实践中如何才能做到这一点，而又不至于在不经意间变得稍大而引发核弹的人来说：首先，裂变过程发生在不同的时间尺度上。虽然大多数中子在裂变事件后很快就被释放出来，但少量的中子只有在裂变开始后经过几次进一步的衰变后，才会由子核释放出来，最长可达10-60秒。因此，如果一个人稍稍超过 $k_{\mathrm{eff}}=1$ ，就会有许多秒的反应时间，直到裂变中产生的所有中子重新进入裂变循环。然而，吸收中子的核反应堆中的控制棒--并因此减少 $k_{\mathrm{eff}}$ --被设计成最多在2秒内全部进入反应堆。

因此，如果 $k_{\mathrm{eff}}$ 在一段时间内大于1，正如不断增长的中子通量所显示的那样，人们有10-60秒的时间来调节核反应。调节可以通过持续监测中子通量来实现，必要时通过将吸收中子的控制棒移入或移出反应堆几毫米来增加或减少中子通量。在更大的范围内，冷却反应堆的水含有硼，一种良好的中子吸收剂。每隔几个小时，通过添加硼或稀释冷却剂来调整硼的浓度。

最后，一些吸收和散射反应有一些内置的稳定性；例如，较高的中子通量导致局部温度升高，这降低了水的密度，因此减少了散射体的数量，而这些散射体是在中子开始裂变事件本身之前将其从高能量调节到低能量所必需的。

在这个教程程序中，我们解决上述 $k$ -两个能量组的特征值问题，我们正在寻找最大的乘法因子 $k_{\mathrm{eff}}$  ，它与最小特征值的逆值加1成正比。为了解决特征值问题，我们一般使用<i>inverse power method</i>的修改版。该算法看起来像这样。

<ol>  <li>  用 $\phi_g^{(0)}$ 和 $k_{\mathrm{eff}}^{(0)}$ 初始化 $\phi_g$ 和 $k_{\mathrm{eff}}$ ，让 $n=1$  .

 <li>  通过@f{eqnarray*}
    s_f^{(n-1)}(x)
    =
    \frac{1}{k_{\mathrm{eff}}^{(n-1)}}
    \sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}^{(n-1)}(x).
  @f}定义所谓的<i>fission source</i>。



 <li>  利用@f{eqnarray*}


    -\nabla \cdot D_g\nabla \phi_g^{(n)}
    +
    \Sigma_{r,g}\phi_g^{(n)}
    =
    \chi_g s_f^{(n-1)}
    +
    \sum_{g'< g} \Sigma_{s,g'\to g} \phi_{g'}^{(n)}
    +
    \sum_{g'> g}\Sigma_{s,g'\to g}\phi_{g'}^{(n-1)}.
  @f}求解所有组通量 $\phi_g,g=1,\ldots,G$  。



 <li>  更新@f{eqnarray*}
    k_{\mathrm{eff}}^{(n)}
    =
    \sum_{g'=1}^G
    \int_{\Omega}\nu\Sigma_{f,g'}(x)
    \phi_{g'}^{(n)}(x)dx.
  @f} 。



 <li>  比较 $k_{\mathrm{eff}}^{(n)}$ 和 $k_{\mathrm{eff}}^{(n-1)}$  。   如果变化大于规定的公差，则设置  $n=n+1$  从步骤2开始重复迭代，否则结束迭代。   </ol> 

需要注意的是，在这个方案中，我们在每次幂迭代中并不完全求解群通量，而是考虑以前只计算  $\phi_{g'}^{(n)}$  的下散射事件  $g'<g$  。上散射仅通过使用旧的迭代器 $\phi_{g'}^{(n-1)}$ 来处理，实质上是假设散射算子是三角形的。这在物理上是有原因的，因为向上散射在中子散射中并不扮演太重要的角色。此外，实践表明，即使使用这种简化方法，反功率迭代也是稳定的。

还要注意的是，我们可以使用很多外推技术来加速上述的功率迭代。然而，这些都没有在这个例子中实现。




<a name="Meshesandmeshrefinement"></a><h3>Meshes and mesh refinement</h3>


人们可能会想，在同一网格上求解各个能量组方程的解是否合适。这个问题可以归结为： $\phi_g$ 和 $\phi_{g'}$ 是否会有类似的光滑度特性？如果是这样的话，那么对两者使用相同的网格是合适的；一个典型的应用可能是化学燃烧，通常所有或大多数化学物种的浓度在火焰前沿快速变化。事实证明，通过观察本教程程序结果部分显示的图形就会发现，然而这里的情况并非如此：由于不同能量组的扩散系数不同，快中子（在小群数 $g$ 的仓中）有一个非常平滑的通量函数，而慢中子（在大群数的仓中）受当地材料特性的影响更大，如果像我们这里计算的情况一样，系数粗糙，则有一个相应的粗糙解决方案。因此，我们要使用不同的网格来计算每个能量组。

这有两个影响，我们将不得不考虑。首先，我们需要找到一种方法来单独细化这些网格。第二，为逆功率迭代组装源项，我们必须将定义在网格 $g'$ 上的解 $\phi_{g'}^{(n)}$ 与定义在网格 $g$ 上的形状函数进行整合，这将成为一项更为复杂的任务。




<a name="Meshrefinement"></a><h4>Mesh refinement</h4>


我们使用通常的范式：在一个给定的网格上求解，然后为每个网格的每个单元评估一个误差指标。因为它非常方便，我们再次使用Kelly, Gago, Zienkiewicz和Babuska的事后误差估计器，它通过整合每个单元面的解的梯度跳跃来接近每个单元的误差。利用这一点，我们得到指标

@f{eqnarray*}
\eta_{g,K}, \qquad g=1,2,\ldots,G,\qquad K\in{\cal T}_g,


@f}

其中 ${\cal T}_g$ 是用于解决 $\phi_g$ 的三角法。问题是该如何处理这个问题。其一，很明显，只细化那些误差指标最高的单元可能会导致不好的结果。为了理解这一点，必须认识到 $\eta_{g,K}$ 与 $\phi_g$ 的二阶导数成比例。换句话说，如果我们有两个能量组 $g=1,2$ ，它们的解同样平滑，但其中一个大了一万倍，例如，那么只有该网格的单元被细化，而小幅度的解的网格将保持粗糙。这可能不是人们想要的，因为我们可以认为解的两个部分同样重要。

因此，从本质上讲，我们必须用一个重要系数 $z_g$ 来衡量 $\eta_{g,K}$ ，这个系数表示将 $\phi_g$ 解决到任何特定的精度有多重要。这样的重要系数可以用二元性技术来计算（例如，见step-14教程程序，以及那里引用的Bangerth和Rannacher的书的参考）。然而，我们不会去那里，而只是假设所有的能量组都是同等重要的，因此将以解 $\phi_g$ 的最大值来规范 $\eta_{g,K}$ 组的误差指标 $g$  。然后我们对误差满足以下条件的单元进行细化

@f{eqnarray*}
  \frac{\eta_{g,K}}{\|\phi_g\|_\infty}
  >
  \alpha_1
  \displaystyle{\max_{\begin{matrix}1\le g\le G \\ K\in {\cal T}_g\end{matrix}}
    \frac{\eta_{g,K}}{\|\phi_g\|_\infty}}


@f}

和粗化的细胞，其中

@f{eqnarray*}
  \frac{\eta_{g,K}}{\|\phi_g\|_\infty}
  <
  \alpha_2
  \displaystyle{\max_{\begin{matrix}1\le g\le G \\ K\in {\cal T}_g\end{matrix}}
    \frac{\eta_{g,K}}{\|\phi_g\|_\infty}}.


@f}

我们在代码中选择了［ $\alpha_1=0.3$ ］和［ $\alpha_2=0.01$ ］。注意，这当然会导致不同能量组的不同网格。

上面的策略实质上意味着以下几点。如果对于能量组 $g$ 来说，有许多单元 $K\in {\cal T}_g$ 的误差很大，例如因为解决方案在整体上非常粗糙，那么许多单元将高于阈值。另一方面，如果有几个单元的误差较大，而许多单元的误差较小，例如因为除少数地方外，解决方案总体上相当平滑，那么只有少数有较大误差的单元将被细化。因此，该策略允许网格很好地跟踪相应的解决方案的全局平滑性属性。




<a name="Assemblingtermsondifferentmeshes"></a><h4>Assembling terms on different meshes</h4>


如上所述，多组细化策略导致不同解的网格不同  $\phi_g$  。那么问题出在哪里呢？实质上是这样的：在特征值迭代的第3步中，我们要像往常一样通过与定义在网格上的测试函数 $\varphi_g^i$ 相乘来形成要计算的方程的弱形式 $\phi_g^{(n)}$ ；在这个过程中，我们要计算包含以下形式的项的右手向量。

@f{eqnarray*}
  F_i = \int_\Omega f(x) \varphi_g^i(x) \phi_{g'}(x) \ dx,


@f}

其中 $f(x)$ 是用于特征值方程右侧的系数函数 $\Sigma_{s,g'\to g}$ 或 $\nu\chi_g\Sigma_{f,g'}$ 中的一个。现在的困难是， $\phi_{g'}$ 是定义在能量组 $g'$ 的网格上，即它可以扩展为 $\phi_{g'}(x)=\sum_j\phi_{g'}^j \varphi_{g'}^j(x)$ ，基函数 $\varphi_{g'}^j(x)$ 定义在网格 $g'$ 。因此，对右边的贡献可以写成

@f{eqnarray*}
  F_i = \sum_j \left\{\int_\Omega f(x) \varphi_g^i(x) \varphi_{g'}^j(x)
  \ dx \right\} \phi_{g'}^j ,


@f}

另一方面，测试函数  $\varphi_g^i(x)$  是在网格  $g$  上定义的。这意味着我们不能将积分  $\Omega$  分割成网格  $g$  或  $g'$  上的积分，因为其他的基函数可能没有定义在这些单元上。

这个问题的解决方案在于， $g$ 和 $g'$ 的网格都是通过自适应细化从一个共同的粗略的网格中得到。因此，我们总能找到一组单元，我们用 ${\cal T}_g \cap
{\cal T}_{g'}$ 表示，它们满足以下条件。   <ul>   <li>  这些单元的联合覆盖了整个领域，并且  <li>  一个单元  $K \in {\cal T}_g \cap {\cal T}_{g'}$  在两个网格中至少有一个是活动的。   </ul>  构建这个集合的方法是，取粗略网格的每个单元，做以下步骤。(i) 如果该单元在 ${\cal T}_g$ 或 ${\cal T}_{g'}$ 上处于活动状态，则将该单元加入该集合；(ii) 否则，即如果该单元在两个网格上都有子节点，则对该单元的每个子节点进行步骤(i)。事实上，deal.II有一个函数 GridTools::get_finest_common_cells ，可以准确地计算出在两个网格中至少有一个处于活动状态的单元的集合。

有了这个，我们可以把上述积分写成如下。

@f{eqnarray*}
  F_i
  =
  \sum_{K \in {\cal T}_g \cap {\cal T}_{g'}}
  \sum_j \left\{\int_K f(x) \varphi_g^i(x) \varphi_{g'}^j(x)
  \ dx \right\} \phi_{g'}^j.


@f}

在代码中，我们在函数 <code>NeutronDiffusionProblem::assemble_rhs</code> 中计算右手边，其中（除其他外）我们在常见的最精炼的单元格集合上循环，对这些单元格的每一对调用函数 <code>NeutronDiffusionProblem::assemble_common_cell</code> 。

根据结构，现在有三种情况需要考虑。<ol>  <li>  单元 $K$ 在两个网格上都是有效的，也就是说，基函数 $\varphi_g^i$ 以及 $\varphi_{g'}^j$ 都是在 $K$ 上定义。   <li>  单元 $K$ 在网格 $g$ 上是有效的，但在 $g'$ 上不是，即 $\varphi_g^i$ 是在 $K$ 上定义的，而 $\varphi_{g'}^j$ 是在 $K$ 的子网格上定义的。   <li>  单元 $K$ 在网格 $g'$ 上是有效的，但在 $g$ 上不是，其结论与(ii)相反。   </ol> 

为了计算上面的右手边，我们就需要对这三种情况有不同的代码，如下所示。<ol>  <li>  如果单元 $K$ 在两个网格上都是活动的，那么我们可以直接评估积分。事实上，我们甚至不必理会基函数 $\varphi_{g'}$ ，因为我们所需要的只是 $\phi_{g'}$ 在正交点的值。我们可以使用 FEValues::get_function_values 函数来完成这个任务。这在 <code>NeutronDiffusionProblem::assemble_common_cell</code> 函数中直接完成。

 <li>  如果单元格  $K$  在网格  $g$  上是有效的，而不是  $g'$  ，那么基函数  $\varphi_{g'}^j$  只能定义在子单元  $K_c,0\le c<2^{\texttt{dim}}$  上，或者在这些子单元  $K$  上被精炼一次以上。

  让我们假设 $K$ 在网格 $g'$ 上只比在网格 $g$ 上多精炼一次。利用我们使用嵌入式有限元空间的事实，即一个网格上的每个基函数可以写成下一个细化网格上的基函数的线性组合，我们可以将 $\phi_g^i$ 对子单元 $K_c$ 的限制扩展为定义在该子单元上的基函数（即定义了基函数 $\varphi_{g'}^l$ 的单元上）。   @f{eqnarray*}
    \phi_g^i|_{K_c} = B_c^{il} \varphi_{g'}^l|_{K_c}.
  @f}

  在这里，以及在下文中，对出现两次的指数进行求和是隐含的。矩阵 $B_c$ 是将数据从一个单元格内插到其 $c$ 的子单元的矩阵。

  那么我们可以把单元格 $K$ 对右侧分量 $F_i$ 的贡献写成@f{eqnarray*}
    F_i|_K
    &=&
    \left\{ \int_K f(x) \varphi_g^i(x) \varphi_{g'}^j(x)
    \ dx \right\} \phi_{g'}^j
    \\
    &=&
    \left\{
    \sum_{0\le c<2^{\texttt{dim}}}
    B_c^{il} \int_{K_c} f(x) \varphi_{g'}^l(x) \varphi_{g'}^j(x)
    \ dx \right\} \phi_{g'}^j.
  @f} 。

  在矩阵符号中，这可以写成@f{eqnarray*}
    F_i|_K
    =
    \sum_{0\le c<2^{\texttt{dim}}}
    F_i|_{K_c},
    \qquad
    \qquad
    F_i|_{K_c} = B_c^{il} M_{K_c}^{lj}  \phi_{g'}^j
    = (B_c M_{K_c})^{ij} \phi_{g'}^j,
  @f}

  其中 $M_{K_c}^{lj}=\int_{K_c} f(x) \varphi_{g'}^l(x) \varphi_{g'}^j(x)$ 是单元格 $K$ 的子 $c$ 上的加权质量矩阵。

  下一个问题是，如果 $K$ 的一个子集 $K_c$ 没有被激活，会发生什么？然后，我们必须递归地应用这个过程，即我们必须将基础函数 $\varphi_g^i$ 插值到 $K$ 的子 $K_c$ 上，然后插值到该单元的子 $K_{cc'}$ 上，插值到该单元的子 $K_{cc'c''}$ 上，等等，直到我们找到一个活动单元。然后，我们必须将单元格 $K$ 的所有子代、孙代等的贡献相加，其贡献形式为@f{eqnarray*}
    F_i|_{K_{cc'}} = (B_cB_{c'} M_{K_{cc'}})^{ij}  \phi_{g'}^j,
  @f} 。

  或@f{eqnarray*}
    F_i|_{K_{cc'c''}} = (B_c B_{c'} B_{c''}M_{K_{cc'c''}})^{ij}
    \phi_{g'}^j,
  @f}

  等等。我们递归地做这个过程，即如果我们坐在单元格 $K$ 上，看到它在网格 $g'$ 上有孩子，那么我们就用一个身份矩阵调用一个函数 <code>assemble_case_2</code> ；该函数将把它的参数从左边乘以延长矩阵；如果该单元格还有孩子，它将用这个新矩阵调用自己，否则它将进行整合。

 <li>  最后一种情况是 $K$ 在网格 $g'$ 上是有效的，但在网格 $g$ 上不是。在这种情况下，我们必须用定义在单元格  $K$  上的基函数来表达基函数  $\varphi_{g'}^j$  ，而不是像以前那样用  $\varphi_g^i$  来表达。这当然是以完全相同的方式进行的。如果 $K$ 的子单元在网格 $g$ 上是活动的，那么就会导致表达式@f{eqnarray*}
    F_i|_K
    &=&
    \left\{ \int_K f(x) \varphi_g^i(x) \varphi_{g'}^j(x)
    \ dx \right\} \phi_{g'}^j
    \\
    &=&
    \left\{
    \sum_{0\le c<2^{\texttt{dim}}}
    \int_{K_c} f(x) \varphi_g^i(x) B_c^{jl} \varphi_{g}^l(x)
    \ dx \right\} \phi_{g'}^j.
  @f} 。

  在矩阵符号中，这个表达式现在读作@f{eqnarray*}
    F_i|_K
    =
    \sum_{0\le c<2^{\texttt{dim}}}
    F_i|_{K_c},
    \qquad
    \qquad
    F_i|_{K_c} = M_{K_c}^{il} B_c^{jl}  \phi_{g'}^j
    =
    (M_{K_c} B_c^T)^{ij} \phi_{g'}^j,
  @f} 。

  而相应地，对于单元格 $K$ 在网格 $g$ 上被精炼一次以上的情况：@f{eqnarray*}
    F_i|_{K_{cc'}} = (M_{K_{cc'}} B_{c'}^T B_c^T)^{ij}  \phi_{g'}^j,
  @f} 。

  或@f{eqnarray*}
    F_i|_{K_{cc'c''}} = (M_{K_{cc'c''}} B_{c''}^T B_{c'}^T B_c^T)^{ij}
    \phi_{g'}^j,
  @f}

  等。换句话说，这个过程与之前的工作方式完全相同，只是我们必须采取延长矩阵的转置，并需要从另一侧乘以质量矩阵。   </ol> 


情况（二）和（三）的表达式可以理解为将标量积 $(f \varphi_g^i, \varphi_{g'}^j)_K$ 中的左或右基函数反复插值到子单元上，然后在最后的单元上形成内积（质量矩阵）。为了使这些情况的对称性更加明显，我们可以这样写：对于情况（二），我们有

@f{eqnarray*}
  F_i|_{K_{cc'\cdots c^{(k)}}}
  = [B_c B_{c'} \cdots B_{c^{(k)}} M_{K_{cc'\cdots c^{(k)}}}]^{ij}
    \phi_{g'}^j,


@f}

而对于情况（三），我们得到

@f{eqnarray*}
  F_i|_{K_{cc'\cdots c^{(k)}}}
  = [(B_c B_{c'} \cdots B_{c^{(k)}} M_{K_{cc'\cdots c^{(k)}}})^T]^{ij}
    \phi_{g'}^j,


@f}






<a name="Descriptionofthetestcase"></a><h3>Description of the test case</h3>


一个核反应堆的堆芯是由不同类型的组件组成的。一个组件基本上是可以在反应堆内和外移动的最小单元，通常是矩形或方形。然而，组件并不是固定的单位，因为它们是由不同的燃料棒、控制棒和仪器元件组成的复杂晶格组装而成的，这些元件通过永久连接到燃料棒上的间隔件来保持彼此的位置。使事情更加复杂的是，在反应堆中同时使用不同种类的组件，这些组件在其组成的燃料棒的类型和排列上有所不同。

显然，组件的排列以及组件内燃料棒的排列都会影响反应堆内中子通量的分布（这一事实通过查看本方案结果部分中下面显示的解决方案就会很明显）。例如，燃料棒在铀235或钚239的富集程度上彼此不同。另一方面，控制棒具有零裂变，但散射和吸收截面不为零。

这整个安排将使描述或空间上的材料参数变得非常复杂。它不会变得更简单，但我们将做一个近似：我们将每个圆柱形棒和周围的水所居住的体积合并成二次横截面的体积，变成所谓的 "夹子单元"，用核数据库和中子光谱的知识获得这些单元的同质化材料数据。同质化使所有材料数据在带有新燃料的反应堆的求解域上成为片状常数。然后为一个点所在的四元组查询空间相关的材料参数，然后为这个四元组中的四元针单元查询。

在这个教程程序中，我们模拟了一个由 $4
\times 4$ 组件组成的反应堆的四分之一。我们使用对称性（Neumann）边界条件，将问题缩小到四分之一的领域，因此只模拟 $2\times 2$ 套组件。其中两个将是UO ${}_2$ 燃料，另外两个是MOX燃料。这些组件中的每一个都由不同成分的 $17\times 17$ 棒组成。因此，我们总共创造了一个 $34\times 34$ 棒子的晶格。为了使以后的事情更简单，我们通过创建一个 $34\times 34$ 单元的粗大网格来反映这一事实（尽管领域是一个正方形，我们通常会使用一个单元）。在deal.II中，每个单元都有一个 <code>material_id</code> ，可以用来将每个单元与一个特定的数字联系起来，识别这个单元的体积是由什么材料制成的；我们将使用这个材料ID来识别在这个测试案例中使用的8种不同的杆子中的哪一种组成了一个特定的单元。请注意，在网格细化后，单元格的子代会继承材料ID，这样即使在网格细化后也能简单地跟踪材料。

在结果部分显示的图像中，棒的排列将清晰可见。材料和两个能量组的横截面取自OECD/NEA的基准问题。详细的配置和材料数据在代码中给出。




<a name="Whattheprogramdoesandhowitdoesthat"></a><h3>What the program does (and how it does that)</h3>


作为对程序具体工作的粗略概述，这里是基本布局：从每个能量组相同的粗略网格开始，我们计算反特征值迭代以计算给定网格集上的 $k$ 特征值。当特征值的变化低于一定的容忍度时，我们停止这些迭代，然后写出每个能量组的网格和解，供图形程序检查。由于解决方案的网格是不同的，我们必须为每个能量组生成一个单独的输出文件，而不是能够将所有能量组的解决方案添加到同一个文件中。

在这之后，我们按照上面某一节的解释对每个网格的误差指标进行评估，并独立地对每个网格的单元进行细化和粗化。由于特征值迭代是相当昂贵的，我们不想在新的网格上重新开始；相反，我们使用SolutionTransfer类在网格细化时将前一个网格的解插到下一个网格。一个简单的实验会让你相信，这比我们省略这一步要省事得多。这样做之后，我们在下一组网格上继续进行特征值迭代。

该程序由一个参数文件控制，使用ParameterHandler类。我们将在本教程的结果部分展示一个参数文件。现在只需说它控制了所使用的有限元的多项式程度、能量组的数量（尽管目前实现的只是2组问题的系数）、停止反特征值迭代的容忍度以及我们要做的细化循环的数量。


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
 * 我们从一堆包含文件开始，这些文件已经在以前的教程程序中解释过了。一个新的文件是  <code>timer.h</code>  : 这是第一个使用Timer类的例子程序。Timer同时记录了经过的挂钟时间（即安装在墙上的时钟所测量的时间）和CPU时钟时间（当前进程在CPU上使用的时间）。我们将在下面使用一个Timer来测量每个网格细化周期所需的CPU时间。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/base/thread_management.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparsity_pattern.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 我们使用下一个include文件来访问块向量，它为我们提供了一种方便的方式来管理所有能量组的解和右手向量。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/block_vector.h> 
 * 
 * @endcode
 * 
 * 这个文件是用来将解从一个网格转移到另一个不同的网格。我们在每次网格迭代后初始化解法时使用它。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * @endcode
 * 
 * 当在一个网格上定义的函数与另一个网格上定义的形状函数进行整合时，我们需要一个函数  @p get_finest_common_cells  (在介绍中已经讨论过)，它定义在以下头文件中。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h> 
 * 
 * @endcode
 * 
 * 我们使用一个来自boost的小工具类来保存输出流的状态（见下面的 <code>run</code> 函数）。
 * 

 * 
 * 
 * @code
 * #include <boost/io/ios_state.hpp> 
 * 
 * @endcode
 * 
 * 这里还有两个C++标准头，我们用它们来定义列表数据类型，以及微调我们生成的输出。
 * 

 * 
 * 
 * @code
 * #include <list> 
 * #include <iomanip> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step28 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Materialdata"></a> 
 * <h3>Material data</h3>
 * 

 * 
 * 首先，我们需要定义一个类，为主类提供材料数据（包括扩散系数、清除截面、散射截面、裂变截面和裂变光谱）。
 * 

 * 
 * 构造函数的参数决定了我们为多少个能量组设置了相关的表格。目前，这个程序只包括2个能量组的数据，但是一个更复杂的程序可能也能为更多的能量组初始化数据结构，这取决于在参数文件中选择了多少个能量组。
 * 

 * 
 * 对于每个不同的系数类型，都有一个函数来返回该系数在特定能量组（或能量组的组合，如分布截面 $\chi_g\nu\Sigma_{f,g'}$ 或散射截面 $\Sigma_{s,g'\to g}$ ）的值。除了能量组之外，这些系数还取决于燃料或控制棒的类型，正如介绍中所解释的那样。因此，这些函数需要一个额外的参数， @p  material_id，以确定特定种类的棒。在这个程序中，我们使用 <code>n_materials=8</code> 不同种类的棒子。
 * 

 * 
 * 除了散射截面，每个系数都可以表示为一个二维浮点数组中的一个条目，该数组由能量组编号以及材料ID索引。表类模板是存储此类数据的理想方式。最后，散射系数取决于两个能量组的索引，因此需要存储在一个三维数组中，为此我们再次使用表类，这时第一个模板参数（表示数组的维度）当然需要是三。
 * 

 * 
 * 
 * @code
 *   class MaterialData 
 *   { 
 *   public: 
 *     MaterialData(const unsigned int n_groups); 
 * 
 *     double get_diffusion_coefficient(const unsigned int group, 
 *                                      const unsigned int material_id) const; 
 *     double get_removal_XS(const unsigned int group, 
 *                           const unsigned int material_id) const; 
 *     double get_fission_XS(const unsigned int group, 
 *                           const unsigned int material_id) const; 
 *     double get_fission_dist_XS(const unsigned int group_1, 
 *                                const unsigned int group_2, 
 *                                const unsigned int material_id) const; 
 *     double get_scattering_XS(const unsigned int group_1, 
 *                              const unsigned int group_2, 
 *                              const unsigned int material_id) const; 
 *     double get_fission_spectrum(const unsigned int group, 
 *                                 const unsigned int material_id) const; 
 * 
 *   private: 
 *     const unsigned int n_groups; 
 *     const unsigned int n_materials; 
 * 
 *     Table<2, double> diffusion; 
 *     Table<2, double> sigma_r; 
 *     Table<2, double> nu_sigma_f; 
 *     Table<3, double> sigma_s; 
 *     Table<2, double> chi; 
 *   }; 
 * 
 * @endcode
 * 
 * 该类的构造函数用于初始化所有材料数据数组。它需要能量组的数量作为参数（如果该值不等于2，就会抛出一个错误，因为目前只实现了两个能量组的数据；但是，使用这个参数，该函数仍然是灵活的，可以扩展到未来）。在开始的成员初始化部分，它也将数组的大小调整为正确的大小。
 * 

 * 
 * 目前，材料数据被存储为8种不同类型的材料。这一点在将来也可以很容易地被扩展。
 * 

 * 
 * 
 * @code
 *   MaterialData::MaterialData(const unsigned int n_groups) 
 *     : n_groups(n_groups) 
 *     , n_materials(8) 
 *     , diffusion(n_materials, n_groups) 
 *     , sigma_r(n_materials, n_groups) 
 *     , nu_sigma_f(n_materials, n_groups) 
 *     , sigma_s(n_materials, n_groups, n_groups) 
 *     , chi(n_materials, n_groups) 
 *   { 
 *     switch (this->n_groups) 
 *       { 
 *         case 2: 
 *           { 
 *             for (unsigned int m = 0; m < n_materials; ++m) 
 *               { 
 *                 diffusion[m][0] = 1.2; 
 *                 diffusion[m][1] = 0.4; 
 *                 chi[m][0]       = 1.0; 
 *                 chi[m][1]       = 0.0; 
 *                 sigma_r[m][0]   = 0.03; 
 *                 for (unsigned int group_1 = 0; group_1 < n_groups; ++group_1) 
 *                   for (unsigned int group_2 = 0; group_2 < n_groups; ++group_2) 
 *                     sigma_s[m][group_1][group_2] = 0.0; 
 *               } 
 * 
 *             diffusion[5][1] = 0.2; 
 * 
 *             sigma_r[4][0] = 0.026; 
 *             sigma_r[5][0] = 0.051; 
 *             sigma_r[6][0] = 0.026; 
 *             sigma_r[7][0] = 0.050; 
 * 
 *             sigma_r[0][1] = 0.100; 
 *             sigma_r[1][1] = 0.200; 
 *             sigma_r[2][1] = 0.250; 
 *             sigma_r[3][1] = 0.300; 
 *             sigma_r[4][1] = 0.020; 
 *             sigma_r[5][1] = 0.040; 
 *             sigma_r[6][1] = 0.020; 
 *             sigma_r[7][1] = 0.800; 
 * 
 *             nu_sigma_f[0][0] = 0.0050; 
 *             nu_sigma_f[1][0] = 0.0075; 
 *             nu_sigma_f[2][0] = 0.0075; 
 *             nu_sigma_f[3][0] = 0.0075; 
 *             nu_sigma_f[4][0] = 0.000; 
 *             nu_sigma_f[5][0] = 0.000; 
 *             nu_sigma_f[6][0] = 1e-7; 
 *             nu_sigma_f[7][0] = 0.00; 
 * 
 *             nu_sigma_f[0][1] = 0.125; 
 *             nu_sigma_f[1][1] = 0.300; 
 *             nu_sigma_f[2][1] = 0.375; 
 *             nu_sigma_f[3][1] = 0.450; 
 *             nu_sigma_f[4][1] = 0.000; 
 *             nu_sigma_f[5][1] = 0.000; 
 *             nu_sigma_f[6][1] = 3e-6; 
 *             nu_sigma_f[7][1] = 0.00; 
 * 
 *             sigma_s[0][0][1] = 0.020; 
 *             sigma_s[1][0][1] = 0.015; 
 *             sigma_s[2][0][1] = 0.015; 
 *             sigma_s[3][0][1] = 0.015; 
 *             sigma_s[4][0][1] = 0.025; 
 *             sigma_s[5][0][1] = 0.050; 
 *             sigma_s[6][0][1] = 0.025; 
 *             sigma_s[7][0][1] = 0.010; 
 * 
 *             break; 
 *           } 
 * 
 *         default: 
 *           Assert(false, 
 *                  ExcMessage( 
 *                    "Presently, only data for 2 groups is implemented")); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 接下来是返回给定材料和能量组的系数值的函数。它们所做的就是确保给定的参数在允许的范围内，然后在相应的表格中查找相应的值。
 * 

 * 
 * 
 * @code
 *   double 
 *   MaterialData::get_diffusion_coefficient(const unsigned int group, 
 *                                           const unsigned int material_id) const 
 *   { 
 *     Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
 *     Assert(material_id < n_materials, 
 *            ExcIndexRange(material_id, 0, n_materials)); 
 * 
 *     return diffusion[material_id][group]; 
 *   } 
 * 
 *   double MaterialData::get_removal_XS(const unsigned int group, 
 *                                       const unsigned int material_id) const 
 *   { 
 *     Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
 *     Assert(material_id < n_materials, 
 *            ExcIndexRange(material_id, 0, n_materials)); 
 * 
 *     return sigma_r[material_id][group]; 
 *   } 
 * 
 *   double MaterialData::get_fission_XS(const unsigned int group, 
 *                                       const unsigned int material_id) const 
 *   { 
 *     Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
 *     Assert(material_id < n_materials, 
 *            ExcIndexRange(material_id, 0, n_materials)); 
 * 
 *     return nu_sigma_f[material_id][group]; 
 *   } 
 * 
 *   double MaterialData::get_scattering_XS(const unsigned int group_1, 
 *                                          const unsigned int group_2, 
 *                                          const unsigned int material_id) const 
 *   { 
 *     Assert(group_1 < n_groups, ExcIndexRange(group_1, 0, n_groups)); 
 *     Assert(group_2 < n_groups, ExcIndexRange(group_2, 0, n_groups)); 
 *     Assert(material_id < n_materials, 
 *            ExcIndexRange(material_id, 0, n_materials)); 
 * 
 *     return sigma_s[material_id][group_1][group_2]; 
 *   } 
 * 
 *   double 
 *   MaterialData::get_fission_spectrum(const unsigned int group, 
 *                                      const unsigned int material_id) const 
 *   { 
 *     Assert(group < n_groups, ExcIndexRange(group, 0, n_groups)); 
 *     Assert(material_id < n_materials, 
 *            ExcIndexRange(material_id, 0, n_materials)); 
 * 
 *     return chi[material_id][group]; 
 *   } 
 * 
 * @endcode
 * 
 * 计算裂变分布截面的函数略有不同，因为它将其值计算为另外两个系数的乘积。我们不需要在这里检查参数，因为这在我们调用其他两个相关函数时已经发生了，尽管这样做可能也无妨。
 * 

 * 
 * 
 * @code
 *   double MaterialData::get_fission_dist_XS(const unsigned int group_1, 
 *                                            const unsigned int group_2, 
 *                                            const unsigned int material_id) const 
 *   { 
 *     return (get_fission_spectrum(group_1, material_id) * 
 *             get_fission_XS(group_2, material_id)); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeEnergyGroupcodeclass"></a> 
 * <h3>The <code>EnergyGroup</code> class</h3>
 * 

 * 
 * 第一个有趣的类是包含所有特定于单个能量组的东西。为了将那些属于同一个对象的东西分组，我们声明了一个结构，该结构包含了用于单个能量组的网格的Triangulation和DoFHandler对象，以及一些其他对象和成员函数，我们将在下面的章节中讨论。
 * 

 * 
 * 这个类的主要原因如下：对于正向问题（有指定的右手边）和特征值问题，人们通常要解决一连串的问题，而不是完全耦合的问题，每个能量组。一旦意识到单一能量组的系统矩阵是对称和正定的（它只是一个扩散算子），而完全耦合问题的矩阵通常是非对称和非定值的，这就可以理解了。如果涉及几个以上的能量组，它也是非常大和相当完整的。
 * 

 * 
 * 让我们先看看在有外部右手的情况下要解决的方程（对于时间无关的情况）。
 * @f{eqnarray*} -\nabla
 * \cdot(D_g(x) \nabla \phi_g(x)) + \Sigma_{r,g}(x)\phi_g(x) =
 * \chi_g\sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}(x) + \sum_{g'\ne
 * g}\Sigma_{s,g'\to g}(x)\phi_{g'}(x) + s_{\mathrm{ext},g}(x) 
 * @f}
 * 

 * 
 * 我们通常会通过将右手边的所有项与  $g'=g$  移到左手边来解决这个方程，并求出  $\phi_g$  。当然，我们还不知道 $\phi_{g'}$ ，因为这些变量的方程包括涉及 $\phi_g$ 的右侧项。在这种情况下，通常的做法是进行迭代：计算 
 * @f{eqnarray*} -\nabla \cdot(D_g(x) \nabla \phi^{(n)}_g(x)) &+&
 * \Sigma_{r,g}(x)\phi^{(n)}_g(x) \\ &=&
 * \chi_g\sum_{g'=1}^{g-1}\nu\Sigma_{f,g'}(x)\phi^{(n)}_{g'}(x) +
 * \chi_g\sum_{g'=g}^G\nu\Sigma_{f,g'}(x)\phi^{(n-1)}_{g'}(x) + \sum_{g'\ne
 * g, g'<g}\Sigma_{s,g'\to g}(x)\phi^{(n)}_{g'}(x) + \sum_{g'\ne g,
 * g'>g}\Sigma_{s,g'\to g}(x)\phi^{(n-1)}_{g'}(x) + s_{\mathrm{ext},g}(x)
 * @f} 。
 * 

 * 
 * 换句话说，我们一个一个地解方程，如果 $g'\ge g$ ，就用上一次迭代的 $\phi_{g'}$ 的值，如果 $g'<g$ ，就用本次迭代已经计算的 $\phi_{g'}$ 的值。
 * 

 * 
 * 在计算特征值时，我们做了一个非常类似的迭代，只是我们没有外部的右手边，而且每次迭代后的解都会被缩放，正如在介绍中所解释的。
 * 

 * 
 * 在任何一种情况下，如果我们所做的只是让下面这一类人具备这些能力，那么这两种情况就可以共同处理。(i) 形成左手边的矩阵，(ii) 形成组内右手边的贡献，即涉及不相干的来源，(iii) 形成源于组  $g'$  的对右手边的贡献。这个类正是做这些工作（以及一些簿记工作，如网格细化、设置矩阵和向量等）。另一方面，这个类本身并不知道有多少个能量组，特别是它们之间的相互作用，也就是说，外部迭代的样子（以及因此我们是解决一个特征值还是一个直接问题）的决定是留给本程序下面的NeutronDiffusionProblem类。
 * 

 * 
 * 所以让我们来看看这个类和它的接口。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class EnergyGroup 
 *   { 
 *   public: 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupcodepublicmemberfunctions"></a> 
 * <h5><code>EnergyGroup</code> public member functions</h5>
 * 

 * 
 * 该类有相当数量的公共成员函数，因为其操作方式是由外部控制的，因此所有做重要事情的函数都需要从另一个类中调用。让我们从记账开始：该类显然需要知道它所代表的能量组，使用哪些材料数据，以及从哪个粗略的网格开始。构造函数接收这些信息，并通过这些信息初始化相关的成员变量（见下文）。
 * 

 * 
 * 然后，我们还需要设置线性系统的函数，即在给定的有限元对象的情况下，正确地确定矩阵的大小和它的稀疏模式等。 <code>setup_linear_system</code> 函数就是这样做的。最后，对于这个初始块，有两个函数可以返回这个对象中使用的活动单元和自由度的数量--利用这一点，我们可以使三角形和DoF处理成员变量成为私有的，不必授予外部使用，增强了封装性。
 * 

 * 
 * 
 * @code
 *     EnergyGroup(const unsigned int        group, 
 *                 const MaterialData &      material_data, 
 *                 const Triangulation<dim> &coarse_grid, 
 *                 const FiniteElement<dim> &fe); 
 * 
 *     void setup_linear_system(); 
 * 
 *     unsigned int n_active_cells() const; 
 *     unsigned int n_dofs() const; 
 * 
 * @endcode
 * 
 * 然后是为每个迭代和当前能量组组装线性系统的函数。请注意，该矩阵与迭代次数无关，因此在每个细化周期只需计算一次。对于必须在每次逆功率迭代中更新的右手边来说，情况就有点复杂了，而且由于计算它可能涉及到几个不同的网格，正如介绍中所解释的那样，这就更复杂了。为了使事情在解决正向或特征值问题方面更加灵活，我们将右手边的计算分成一个函数，将无关的源和组内贡献（我们将其称为零函数，作为特征值问题的源项）和一个计算来自另一个能量组的右手边的贡献。
 * 

 * 
 * 
 * @code
 *     void assemble_system_matrix(); 
 *     void assemble_ingroup_rhs(const Function<dim> &extraneous_source); 
 *     void assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime); 
 * 
 * @endcode
 * 
 * 接下来我们需要一组函数来实际计算线性系统的解，并对其进行处理（比如计算介绍中提到的裂变源贡献，将图形信息写入输出文件，计算误差指标，或者根据这些标准和阈值实际细化和粗化网格）。所有这些函数以后都可以从驱动类 <code>NeutronDiffusionProblem</code> 中调用，或者你想实现的任何其他类来解决涉及中子通量方程的问题。
 * 

 * 
 * 
 * @code
 *     void solve(); 
 * 
 *     double get_fission_source() const; 
 * 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     void estimate_errors(Vector<float> &error_indicators) const; 
 * 
 *     void refine_grid(const Vector<float> &error_indicators, 
 *                      const double         refine_threshold, 
 *                      const double         coarsen_threshold); 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupcodepublicdatamembers"></a> 
 * <h5><code>EnergyGroup</code> public data members</h5>
 * 

 * 
 * 作为面向对象编程的良好实践，我们通过使它们成为私有的来隐藏大多数数据成员。然而，我们必须允许驱动进程的类访问解向量以及上一次迭代的解，因为在幂迭代中，解向量在每次迭代中都被我们正在寻找的特征值的当前猜测所缩放。
 * 

 * 
 * 
 * @code
 *   public: 
 *     Vector<double> solution; 
 *     Vector<double> solution_old; 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupcodeprivatedatamembers"></a> 
 * <h5><code>EnergyGroup</code> private data members</h5>
 * 

 * 
 * 其余的数据成员是私有的。与之前所有的教程程序相比，唯一的新数据成员是一个存储此对象所代表的能量组的整数，以及此对象的构造函数从驱动类中得到的材料数据对象的引用。同样地，构造函数得到了我们要使用的有限元对象的引用。
 * 

 * 
 * 最后，我们必须在每次迭代中对线性系统应用边界值，即相当频繁。我们不是每次都插值，而是在每个新的网格上插值一次，然后和这个类的所有其他数据一起存储。
 * 

 * 
 * 
 * @code
 *   private: 
 *     const unsigned int  group; 
 *     const MaterialData &material_data; 
 * 
 *     Triangulation<dim>        triangulation; 
 *     const FiniteElement<dim> &fe; 
 *     DoFHandler<dim>           dof_handler; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> system_rhs; 
 * 
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     AffineConstraints<double>                 hanging_node_constraints; 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupcodeprivatememberfunctions"></a> 
 * <h5><code>EnergyGroup</code> private member functions</h5>
 * 

 * 
 * 在这个类中有一个私有成员函数。它递归地走过两个网格的单元，以计算跨组的右手边项。这个算法在本程序的介绍中已经解释过了。这个函数的参数是对一个对象的引用，该对象代表了我们要整合的右手项的能量组，一个指向用于当前能量组的网格单元的迭代器，一个指向另一个网格上相应单元的迭代器，以及将自由度从两个单元中较粗的单元插补到较细的单元的矩阵。
 * 

 * 
 * 
 * @code
 *   private: 
 *     void assemble_cross_group_rhs_recursive( 
 *       const EnergyGroup<dim> &                       g_prime, 
 *       const typename DoFHandler<dim>::cell_iterator &cell_g, 
 *       const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
 *       const FullMatrix<double> &                     prolongation_matrix); 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeEnergyGroupcodeclass"></a> 
 * <h4>Implementation of the <code>EnergyGroup</code> class</h4>
 * 

 * 
 * 这个类的前几个函数大部分是不言自明的。构造函数只设置了几个数据成员，并创建了一个给定三角形的副本，作为该能量组使用的三角形的基础。接下来的两个函数只是从私有数据成员中返回数据，从而使我们能够使这些数据成员私有化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   EnergyGroup<dim>::EnergyGroup(const unsigned int        group, 
 *                                 const MaterialData &      material_data, 
 *                                 const Triangulation<dim> &coarse_grid, 
 *                                 const FiniteElement<dim> &fe) 
 *     : group(group) 
 *     , material_data(material_data) 
 *     , fe(fe) 
 *     , dof_handler(triangulation) 
 *   { 
 *     triangulation.copy_triangulation(coarse_grid); 
 *     dof_handler.distribute_dofs(fe); 
 *   } 
 * 
 *   template <int dim> 
 *   unsigned int EnergyGroup<dim>::n_active_cells() const 
 *   { 
 *     return triangulation.n_active_cells(); 
 *   } 
 * 
 *   template <int dim> 
 *   unsigned int EnergyGroup<dim>::n_dofs() const 
 *   { 
 *     return dof_handler.n_dofs(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupsetup_linear_systemcode"></a> 
 * <h5><code>EnergyGroup::setup_linear_system</code></h5>
 * 

 * 
 * 第一个 "真正的 "函数是在新的网格上或网格细化后设置网格、矩阵等。我们用这个函数来初始化稀疏系统矩阵，以及右手边的向量。如果求解向量之前从未被设置过（如用零大小表示），我们也会初始化它并将其设置为默认值。如果它已经有一个非零的大小，我们就不这么做了（也就是说，这个函数是在网格细化之后调用的），因为在这种情况下，我们希望在不同的网格细化中保留解决方案（这一点我们在 <code>EnergyGroup::refine_grid</code> 函数中做到了）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::setup_linear_system() 
 *   { 
 *     const unsigned int n_dofs = dof_handler.n_dofs(); 
 * 
 *     hanging_node_constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 *     system_matrix.clear(); 
 * 
 *     DynamicSparsityPattern dsp(n_dofs, n_dofs); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     hanging_node_constraints.condense(dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     system_rhs.reinit(n_dofs); 
 * 
 *     if (solution.size() == 0) 
 *       { 
 *         solution.reinit(n_dofs); 
 *         solution_old.reinit(n_dofs); 
 *         solution_old = 1.0; 
 *         solution     = solution_old; 
 *       } 
 * 
 * @endcode
 * 
 * 在这个函数的最后，我们更新边界节点列表和它们的数值，首先清除这个列表和重新插值的边界数值（记住，这个函数是在第一次设置网格后，以及每次网格细化后调用）。
 * 

 * 
 * 为了理解这段代码，有必要认识到我们使用 <code>GridGenerator::subdivided_hyper_rectangle</code> 函数来创建网格(在 <code>NeutronDiffusionProblem::initialize_problem</code> )，其中我们将最后一个参数设置为 <code>true</code>  。这意味着域的边界被 "着色"，也就是说，域的四个（或六个，在3D）边被分配了不同的边界指标。结果是，底部边界得到指标0，顶部一个边界指标1，左右边界分别得到指标2和3。
 * 

 * 
 * 在这个程序中，我们只模拟一个，即右上角的反应器的四分之一。也就是说，我们只想在顶部和右侧边界插值边界条件，而在底部和左侧边界不做任何事情（即施加自然的、无流量的诺伊曼边界条件）。这很容易被推广到任意维度，即我们想在指标为1、3、......的边界上插值，我们在下面的循环中这样做（注意，对 <code>VectorTools::interpolate_boundary_values</code> 的调用是加法的，即它们不首先清除边界值图）。
 * 

 * 
 * 
 * @code
 *     boundary_values.clear(); 
 * 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                2 * i + 1, 
 *                                                Functions::ZeroFunction<dim>(), 
 *                                                boundary_values); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupassemble_system_matrixcode"></a> 
 * <h5><code>EnergyGroup::assemble_system_matrix</code></h5>
 * 

 * 
 * 接下来我们需要函数来组装系统矩阵和右手边。考虑到介绍中列出的方程以及我们在以前的例子程序中看到的内容，组装矩阵是很简单的。注意使用 <code>cell->material_id()</code> 来获取一个单元的材料种类。还请注意我们是如何设置正交公式的顺序的，以便它总是适合于正在使用的有限元。
 * 

 * 
 * 最后，请注意，由于我们在这里只组装了系统矩阵，所以我们还不能消除边界值（我们需要右边的向量来实现）。我们将其推迟到 <code>EnergyGroup::solve</code> 函数中，这时所有的信息都可以得到。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::assemble_system_matrix() 
 *   { 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix = 0; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         const double diffusion_coefficient = 
 *           material_data.get_diffusion_coefficient(group, cell->material_id()); 
 *         const double removal_XS = 
 *           material_data.get_removal_XS(group, cell->material_id()); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               cell_matrix(i, j) += 
 *                 ((diffusion_coefficient * fe_values.shape_grad(i, q_point) * 
 *                     fe_values.shape_grad(j, q_point) + 
 *                   removal_XS * fe_values.shape_value(i, q_point) * 
 *                     fe_values.shape_value(j, q_point)) * 
 *                  fe_values.JxW(q_point)); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             system_matrix.add(local_dof_indices[i], 
 *                               local_dof_indices[j], 
 *                               cell_matrix(i, j)); 
 *       } 
 * 
 *     hanging_node_constraints.condense(system_matrix); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupassemble_ingroup_rhscode"></a> 
 * <h5><code>EnergyGroup::assemble_ingroup_rhs</code></h5>
 * 

 * 
 * 正如 <code>EnergyGroup</code> 类的文档中所解释的，我们将组装右手边分成两部分：组内耦合和跨组耦合。首先，我们需要一个函数来组装这里的一个特定组的右手边，即包括一个无关的源（我们将在特征值问题上设置为零）以及组内裂变贡献。 组内散射已经在清除截面的定义中得到了考虑）。这个函数的工作原理就组装右手边而言是非常标准的，因此不需要更多的评论，只是我们要提到在函数的开头将右手边的向量设置为零--这一点我们不打算为跨组项做，这些跨组项只是添加到右手边的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   EnergyGroup<dim>::assemble_ingroup_rhs(const Function<dim> &extraneous_source) 
 *   { 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values); 
 * 
 *     Vector<double>      cell_rhs(dofs_per_cell); 
 *     std::vector<double> extraneous_source_values(n_q_points); 
 *     std::vector<double> solution_old_values(n_q_points); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_rhs = 0; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         const double fission_dist_XS = 
 *           material_data.get_fission_dist_XS(group, group, cell->material_id()); 
 * 
 *         extraneous_source.value_list(fe_values.get_quadrature_points(), 
 *                                      extraneous_source_values); 
 * 
 *         fe_values.get_function_values(solution_old, solution_old_values); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             cell_rhs(i) += 
 *               ((extraneous_source_values[q_point] + 
 *                 fission_dist_XS * solution_old_values[q_point]) * 
 *                fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupassemble_cross_group_rhscode"></a> 
 * <h5><code>EnergyGroup::assemble_cross_group_rhs</code></h5>
 * 

 * 
 * 对于组装单一能量组方程的右手向量来说，更有趣的函数是将能量组  $g$  和  $g'$  耦合起来。正如介绍中所解释的，我们首先要找到两个能量组的网格所共有的单元集。首先我们调用 <code>get_finest_common_cells</code> 来获得这一对来自两个网格的共同单元的列表。一对单元格中的两个单元格可能都不活跃，但至少有一个是活跃的。然后我们将这些单元对中的每一个交给一个函数，该函数将递归地计算右手边的项。
 * 

 * 
 * 注意，组内耦合在之前已经处理过了，所以我们提前退出这个函数，如果  $g=g'$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   EnergyGroup<dim>::assemble_cross_group_rhs(const EnergyGroup<dim> &g_prime) 
 *   { 
 *     if (group == g_prime.group) 
 *       return; 
 * 
 *     const std::list<std::pair<typename DoFHandler<dim>::cell_iterator, 
 *                               typename DoFHandler<dim>::cell_iterator>> 
 *       cell_list = 
 *         GridTools::get_finest_common_cells(dof_handler, g_prime.dof_handler); 
 * 
 *     for (const auto &cell_pair : cell_list) 
 *       { 
 *         FullMatrix<double> unit_matrix(fe.n_dofs_per_cell()); 
 *         for (unsigned int i = 0; i < unit_matrix.m(); ++i) 
 *           unit_matrix(i, i) = 1; 
 *         assemble_cross_group_rhs_recursive(g_prime, 
 *                                            cell_pair.first, 
 *                                            cell_pair.second, 
 *                                            unit_matrix); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupassemble_cross_group_rhs_recursivecode"></a> 
 * <h5><code>EnergyGroup::assemble_cross_group_rhs_recursive</code></h5>
 * 

 * 
 * 这是最后一个处理在潜在的不同网格上递归组装右手边条款的函数，使用介绍中描述的算法。该函数接收一个代表能量组 $g'$ 的对象的引用，以及能量组 $g$ 和 $g'$ 的网格中相应单元的迭代器。起初，即从上面的函数中调用这个函数时，这两个单元将是两个网格上的匹配单元；然而，这两个单元中的一个可能被进一步细化，我们将递归地调用这个函数，两个迭代器中的一个被原始单元的一个子单元所取代。
 * 

 * 
 * 最后一个参数是介绍中的矩阵乘积矩阵 $B_{c^{(k)}}^T \cdots B_{c'}^T B_c^T$ ，它从两个单元中较粗的单元插值到较细的单元。如果这两个单元格匹配，那么这就是身份矩阵--正是我们最初传递给这个函数的。
 * 

 * 
 * 该函数必须考虑两种情况：两种单元格都没有进一步细化，即没有子代，在这种情况下，我们可以最终组装这对单元格的右侧贡献；以及两种单元格中的一种被进一步细化，在这种情况下，我们必须通过循环未被激活的单元格的子代来不断地进行循环。下面将讨论这两种情况。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::assemble_cross_group_rhs_recursive( 
 *     const EnergyGroup<dim> &                       g_prime, 
 *     const typename DoFHandler<dim>::cell_iterator &cell_g, 
 *     const typename DoFHandler<dim>::cell_iterator &cell_g_prime, 
 *     const FullMatrix<double> &                     prolongation_matrix) 
 *   { 
 * 
 * @endcode
 * 
 * 第一种情况是，两个单元格都没有进一步的细化。在这种情况下，我们可以组装相关条款（见介绍）。这涉及到在两个单元中较细的单元上组装质量矩阵（事实上有两个具有不同系数的质量矩阵，一个用于裂变分布截面 $\chi_g\nu\Sigma_{f,g'}$ ，一个用于散射截面 $\Sigma_{s,g'\to g}$ ）。这是直截了当的，但请注意我们如何通过查看两个单元的细化程度来确定哪个是更细的单元。
 * 

 * 
 * 
 * @code
 *     if (!cell_g->has_children() && !cell_g_prime->has_children()) 
 *       { 
 *         const QGauss<dim>  quadrature_formula(fe.degree + 1); 
 *         const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *         FEValues<dim> fe_values(fe, 
 *                                 quadrature_formula, 
 *                                 update_values | update_JxW_values); 
 * 
 *         if (cell_g->level() > cell_g_prime->level()) 
 *           fe_values.reinit(cell_g); 
 *         else 
 *           fe_values.reinit(cell_g_prime); 
 * 
 *         const double fission_dist_XS = 
 *           material_data.get_fission_dist_XS(group, 
 *                                             g_prime.group, 
 *                                             cell_g_prime->material_id()); 
 * 
 *         const double scattering_XS = 
 *           material_data.get_scattering_XS(g_prime.group, 
 *                                           group, 
 *                                           cell_g_prime->material_id()); 
 * 
 *         FullMatrix<double> local_mass_matrix_f(fe.n_dofs_per_cell(), 
 *                                                fe.n_dofs_per_cell()); 
 *         FullMatrix<double> local_mass_matrix_g(fe.n_dofs_per_cell(), 
 *                                                fe.n_dofs_per_cell()); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
 *             for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
 *               { 
 *                 local_mass_matrix_f(i, j) += 
 *                   (fission_dist_XS * fe_values.shape_value(i, q_point) * 
 *                    fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
 *                 local_mass_matrix_g(i, j) += 
 *                   (scattering_XS * fe_values.shape_value(i, q_point) * 
 *                    fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
 *               } 
 * 
 * @endcode
 * 
 * 现在我们有了所有的插值（延长）矩阵以及局部质量矩阵，所以我们只需要根据两个单元中哪一个更细，形成积
 * @f[
 * F_i|_{K_{cc'\cdots c^{(k)}}} = [B_c B_{c'} \cdots B_{c^{(k)}}
 * M_{K_{cc'\cdots c^{(k)}}}]^{ij} \phi_{g'}^j, 
 * @f]
 * 

 * 
 * 或
 * 

 * 
 * @f[
 * F_i|_{K_{cc'\cdots c^{(k)}}} = [(B_c B_{c'} \cdots B_{c^{(k)}}
 * M_{K_{cc'\cdots c^{(k)}}})^T]^{ij} \phi_{g'}^j, 
 * @f]
 * 。我们使用 <code>vmult</code> 函数提供的矩阵-向量乘积，或者使用 <code>Tvmult</code> 与转置矩阵的乘积来完成。这样做之后，我们将结果转移到能量组的全局右侧向量中  $g$  。
 * 

 * 
 * 

 * 
 * 
 * @code
 *         Vector<double> g_prime_new_values(fe.n_dofs_per_cell()); 
 *         Vector<double> g_prime_old_values(fe.n_dofs_per_cell()); 
 *         cell_g_prime->get_dof_values(g_prime.solution_old, g_prime_old_values); 
 *         cell_g_prime->get_dof_values(g_prime.solution, g_prime_new_values); 
 * 
 *         Vector<double> cell_rhs(fe.n_dofs_per_cell()); 
 *         Vector<double> tmp(fe.n_dofs_per_cell()); 
 * 
 *         if (cell_g->level() > cell_g_prime->level()) 
 *           { 
 *             prolongation_matrix.vmult(tmp, g_prime_old_values); 
 *             local_mass_matrix_f.vmult(cell_rhs, tmp); 
 * 
 *             prolongation_matrix.vmult(tmp, g_prime_new_values); 
 *             local_mass_matrix_g.vmult_add(cell_rhs, tmp); 
 *           } 
 *         else 
 *           { 
 *             local_mass_matrix_f.vmult(tmp, g_prime_old_values); 
 *             prolongation_matrix.Tvmult(cell_rhs, tmp); 
 * 
 *             local_mass_matrix_g.vmult(tmp, g_prime_new_values); 
 *             prolongation_matrix.Tvmult_add(cell_rhs, tmp); 
 *           } 
 * 
 *         std::vector<types::global_dof_index> local_dof_indices( 
 *           fe.n_dofs_per_cell()); 
 *         cell_g->get_dof_indices(local_dof_indices); 
 * 
 *         for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *       } 
 * 
 * @endcode
 * 
 * 另一种情况是，两个单元中的一个被进一步细化。在这种情况下，我们必须在所有的子单元上循环，将现有的矩阵的插值（延长）乘以从现在的单元到其子单元的插值（使用矩阵-矩阵乘法函数 <code>mmult</code> ），然后将结果再次交给这个非常相同的函数，但将有子单元替换为其子单元之一。
 * 

 * 
 * 
 * @code
 *     else 
 *       for (unsigned int child = 0; 
 *            child < GeometryInfo<dim>::max_children_per_cell; 
 *            ++child) 
 *         { 
 *           FullMatrix<double> new_matrix(fe.n_dofs_per_cell(), 
 *                                         fe.n_dofs_per_cell()); 
 *           fe.get_prolongation_matrix(child).mmult(new_matrix, 
 *                                                   prolongation_matrix); 
 * 
 *           if (cell_g->has_children()) 
 *             assemble_cross_group_rhs_recursive(g_prime, 
 *                                                cell_g->child(child), 
 *                                                cell_g_prime, 
 *                                                new_matrix); 
 *           else 
 *             assemble_cross_group_rhs_recursive(g_prime, 
 *                                                cell_g, 
 *                                                cell_g_prime->child(child), 
 *                                                new_matrix); 
 *         } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupget_fission_sourcecode"></a> 
 * <h5><code>EnergyGroup::get_fission_source</code></h5>
 * 

 * 
 * 在（反）功率迭代中，我们使用综合裂变源来更新 $k$  -特征值。鉴于其定义，以下函数基本上是不言自明的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double EnergyGroup<dim>::get_fission_source() const 
 *   { 
 *     const QGauss<dim>  quadrature_formula(fe.degree + 1); 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_JxW_values); 
 * 
 *     std::vector<double> solution_values(n_q_points); 
 * 
 *     double fission_source = 0; 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 * 
 *         const double fission_XS = 
 *           material_data.get_fission_XS(group, cell->material_id()); 
 * 
 *         fe_values.get_function_values(solution, solution_values); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           fission_source += 
 *             (fission_XS * solution_values[q_point] * fe_values.JxW(q_point)); 
 *       } 
 * 
 *     return fission_source; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupsolvecode"></a> 
 * <h5><code>EnergyGroup::solve</code></h5>
 * 

 * 
 * 接下来是一个解决之前组装的线性系统的函数。事情基本是标准的，只是我们把应用边界值的时间推迟到了这里，因为在之前的所有函数中，我们还是在为右边的向量做加法。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::solve() 
 *   { 
 *     hanging_node_constraints.condense(system_rhs); 
 *     MatrixTools::apply_boundary_values(boundary_values, 
 *                                        system_matrix, 
 *                                        solution, 
 *                                        system_rhs); 
 * 
 *     SolverControl            solver_control(system_matrix.m(), 
 *                                  1e-12 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *     hanging_node_constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupestimate_errorscode"></a> 
 * <h5><code>EnergyGroup::estimate_errors</code></h5>
 * 

 * 
 * 网格细化被分成两个函数。第一个函数估计每个单元的误差，用解的大小对其进行归一化处理，并将其返回到作为参数的向量中。调用函数收集所有能量组的所有误差指标，并计算出细化和粗化单元的阈值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::estimate_errors(Vector<float> &error_indicators) const 
 *   { 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(fe.degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       error_indicators); 
 *     error_indicators /= solution.linfty_norm(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGrouprefine_gridcode"></a> 
 * <h5><code>EnergyGroup::refine_grid</code></h5>
 * 

 * 
 * 第二部分是细化网格，给定前一个函数中计算的误差指标和误差阈值，超过这个阈值的单元应被细化，低于这个阈值的单元应被粗化。注意，我们在这里没有使用 <code>GridRefinement</code> 中的任何函数，而是自己设置细化标志。
 * 

 * 
 * 在设置完这些标志后，我们使用SolutionTransfer类将求解向量从旧网格转移到新网格。这里使用的程序在该类的文档中已有详细描述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::refine_grid(const Vector<float> &error_indicators, 
 *                                      const double         refine_threshold, 
 *                                      const double         coarsen_threshold) 
 *   { 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       if (error_indicators(cell->active_cell_index()) > refine_threshold) 
 *         cell->set_refine_flag(); 
 *       else if (error_indicators(cell->active_cell_index()) < coarsen_threshold) 
 *         cell->set_coarsen_flag(); 
 * 
 *     SolutionTransfer<dim> soltrans(dof_handler); 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 *     soltrans.prepare_for_coarsening_and_refinement(solution); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 *     dof_handler.distribute_dofs(fe); 
 *     setup_linear_system(); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     soltrans.interpolate(solution_old, solution); 
 * 
 * @endcode
 * 
 * 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。
 * 

 * 
 * 
 * @code
 *     hanging_node_constraints.distribute(solution); 
 * 
 *     solution_old.reinit(dof_handler.n_dofs()); 
 *     solution_old = solution; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeEnergyGroupoutput_resultscode"></a> 
 * <h5><code>EnergyGroup::output_results</code></h5>
 * 

 * 
 * 该类的最后一个函数在每次网格迭代后输出网格和解。这在以前已经展示过很多次了。唯一值得指出的是使用 <code>Utilities::int_to_string</code> 函数将一个整数转换为其字符串表示。该函数的第二个参数表示我们应使用多少个数字 -- 如果这个值大于1，那么数字将被填充前导零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EnergyGroup<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     const std::string filename = std::string("solution-") + 
 *                                  Utilities::int_to_string(group, 2) + "." + 
 *                                  Utilities::int_to_string(cycle, 2) + ".vtu"; 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeNeutronDiffusionProblemcodeclasstemplate"></a> 
 * <h3>The <code>NeutronDiffusionProblem</code> class template</h3>
 * 

 * 
 * 这是程序的主类，并不是因为它实现了所有的功能（事实上，大部分的功能都在 <code>EnergyGroup</code> 类中实现了），而是因为它包含了决定什么时候计算的驱动算法。它主要是像其他许多教程程序中显示的那样，它有一个公共的 <code>run</code> 函数和私有函数来做其他的事情。在一些地方，我们必须为所有能源组做一些事情，在这种情况下，如果deal.II被配置为多线程，我们将为每个组启动任务，让这些事情并行运行。 关于并行化的策略，请看一下 @ref threads 模块。
 * 

 * 
 * 与以前的例子程序最大的不同是，我们还声明了一个嵌套类，该类有成员变量，用于所有可在输入文件中传递给程序的运行时参数。现在，这些参数是能量组的数量、细化周期的数量、要使用的有限元的多项式程度以及用于确定反幂迭代何时收敛的公差。此外，我们有一个该类的构造函数，将所有这些值设置为默认值，还有一个函数 <code>declare_parameters</code> 向ParameterHandler类描述输入文件中接受哪些参数，还有一个函数 <code>get_parameters</code> 可以从ParameterHandler对象中提取这些参数的值。参见 step-29 ，了解另一个使用ParameterHandler的例子。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class NeutronDiffusionProblem 
 *   { 
 *   public: 
 *     class Parameters 
 *     { 
 *     public: 
 *       Parameters(); 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 *       void        get_parameters(ParameterHandler &prm); 
 * 
 *       unsigned int n_groups; 
 *       unsigned int n_refinement_cycles; 
 * 
 *       unsigned int fe_degree; 
 * 
 *       double convergence_tolerance; 
 *     }; 
 * 
 *     NeutronDiffusionProblem(const Parameters &parameters); 
 * 
 *     void run(); 
 * 
 *   private: 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProblemcodeprivatememberfunctions"></a> 
 * <h5><code>NeutronDiffusionProblem</code> private member functions</h5>
 * 

 * 
 * 这个类中没有那么多的成员函数，因为大部分的功能已经被移到了 <code>EnergyGroup</code> 类中，只是从这个类的 <code>run()</code> 成员函数中调用。留下来的成员函数有不言自明的名字。
 * 

 * 
 * 
 * @code
 *     void initialize_problem(); 
 * 
 *     void refine_grid(); 
 * 
 *     double get_total_fission_source() const; 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProblemcodeprivatemembervariables"></a> 
 * <h5><code>NeutronDiffusionProblem</code> private member variables</h5>
 * 

 * 
 * 接下来，我们有几个成员变量。特别是，这些是（i）对参数对象的引用（由本程序的主函数拥有，并传递给本类的构造函数），（ii）描述输入文件中要求的能量组数量的材料参数的对象，以及（iii）所有能量组将使用的有限元。
 * 

 * 
 * 
 * @code
 *     const Parameters & parameters; 
 *     const MaterialData material_data; 
 *     FE_Q<dim>          fe; 
 * 
 * @endcode
 * 
 * 此外，我们有(iv)目前迭代时计算的特征值的值。事实上，这是在所有能量组之间共享的解决方案的唯一部分--解决方案的所有其他部分，如中子通量是特定于一个或另一个能量组的，因此被存储在描述单一能量组的对象中。
 * 

 * 
 * 
 * @code
 *     double k_eff; 
 * 
 * @endcode
 * 
 * 最后一个计算对象（v）是一个指向能量组对象的数组。当然，这个数组的长度等于参数文件中指定的能量组的数量。
 * 

 * 
 * 
 * @code
 *     std::vector<std::unique_ptr<EnergyGroup<dim>>> energy_groups; 
 * 
 * @endcode
 * 
 * 最后(vi)我们有一个文件流，我们将把总结的输出保存到这个文件中。
 * 

 * 
 * 
 * @code
 *     std::ofstream convergence_table_stream; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeParameterscodeclass"></a> 
 * <h4>Implementation of the <code>Parameters</code> class</h4>
 * 

 * 
 * 在继续实现外层类之前，我们必须实现参数结构的功能。这是很直接的，事实上，对于所有使用ParameterHandler功能的这类参数类来说，看起来都是一样的。因此，我们将不再对此进行评论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   NeutronDiffusionProblem<dim>::Parameters::Parameters() 
 *     : n_groups(2) 
 *     , n_refinement_cycles(5) 
 *     , fe_degree(2) 
 *     , convergence_tolerance(1e-12) 
 *   {} 
 * 
 *   template <int dim> 
 *   void NeutronDiffusionProblem<dim>::Parameters::declare_parameters( 
 *     ParameterHandler &prm) 
 *   { 
 *     prm.declare_entry("Number of energy groups", 
 *                       "2", 
 *                       Patterns::Integer(), 
 *                       "The number of energy different groups considered"); 
 *     prm.declare_entry("Refinement cycles", 
 *                       "5", 
 *                       Patterns::Integer(), 
 *                       "Number of refinement cycles to be performed"); 
 *     prm.declare_entry("Finite element degree", 
 *                       "2", 
 *                       Patterns::Integer(), 
 *                       "Polynomial degree of the finite element to be used"); 
 *     prm.declare_entry( 
 *       "Power iteration tolerance", 
 *       "1e-12", 
 *       Patterns::Double(), 
 *       "Inner power iterations are stopped when the change in k_eff falls " 
 *       "below this tolerance"); 
 *   } 
 * 
 *   template <int dim> 
 *   void NeutronDiffusionProblem<dim>::Parameters::get_parameters( 
 *     ParameterHandler &prm) 
 *   { 
 *     n_groups              = prm.get_integer("Number of energy groups"); 
 *     n_refinement_cycles   = prm.get_integer("Refinement cycles"); 
 *     fe_degree             = prm.get_integer("Finite element degree"); 
 *     convergence_tolerance = prm.get_double("Power iteration tolerance"); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeNeutronDiffusionProblemcodeclass"></a> 
 * <h4>Implementation of the <code>NeutronDiffusionProblem</code> class</h4>
 * 

 * 
 * 现在是 <code>NeutronDiffusionProblem</code> 类。构造函数和析构函数没有什么值得注意的地方。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   NeutronDiffusionProblem<dim>::NeutronDiffusionProblem( 
 *     const Parameters &parameters) 
 *     : parameters(parameters) 
 *     , material_data(parameters.n_groups) 
 *     , fe(parameters.fe_degree) 
 *     , k_eff(std::numeric_limits<double>::quiet_NaN()) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProbleminitialize_problemcode"></a> 
 * <h5><code>NeutronDiffusionProblem::initialize_problem</code></h5>
 * 

 * 
 * 第一个感兴趣的函数是设置反应堆核心的几何形状的函数。这在介绍中会有更详细的描述。
 * 

 * 
 * 该函数的第一部分定义了几何数据，然后创建了一个粗略的网格，其单元数与我们模拟的那部分反应堆堆芯中的燃料棒（或针状单元）的数量相当。正如上面插值边界值时提到的， <code>GridGenerator::subdivided_hyper_rectangle</code> 函数的最后一个参数指定域的两侧应具有唯一的边界指标，这将使我们能够以简单的方式确定哪些边界具有诺伊曼条件，哪些边界具有迪里希特条件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void NeutronDiffusionProblem<dim>::initialize_problem() 
 *   { 
 *     const unsigned int rods_per_assembly_x = 17, rods_per_assembly_y = 17; 
 *     const double       pin_pitch_x = 1.26, pin_pitch_y = 1.26; 
 *     const double       assembly_height = 200; 
 * 
 *     const unsigned int assemblies_x = 2, assemblies_y = 2, assemblies_z = 1; 
 * 
 *     const Point<dim> bottom_left = Point<dim>(); 
 *     const Point<dim> upper_right = 
 *       (dim == 2 ? Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
 *                              assemblies_y * rods_per_assembly_y * pin_pitch_y) : 
 *                   Point<dim>(assemblies_x * rods_per_assembly_x * pin_pitch_x, 
 *                              assemblies_y * rods_per_assembly_y * pin_pitch_y, 
 *                              assemblies_z * assembly_height)); 
 * 
 *     std::vector<unsigned int> n_subdivisions; 
 *     n_subdivisions.push_back(assemblies_x * rods_per_assembly_x); 
 *     if (dim >= 2) 
 *       n_subdivisions.push_back(assemblies_y * rods_per_assembly_y); 
 *     if (dim >= 3) 
 *       n_subdivisions.push_back(assemblies_z); 
 * 
 *     Triangulation<dim> coarse_grid; 
 *     GridGenerator::subdivided_hyper_rectangle( 
 *       coarse_grid, n_subdivisions, bottom_left, upper_right, true); 
 * 
 * @endcode
 * 
 * 该函数的第二部分涉及每种类型的组件的引脚单元的材料数量。在这里，我们定义了四种不同类型的组件，对于这些组件，我们在以下表格中描述了燃料棒的排列。
 * 

 * 
 * 这里描述的装配体来自于介绍中提到的基准，它们是（按照这个顺序）。  <ol>  
 * <li>  'UX'组件。二氧化铀燃料组件，带有24个导向管和一个中央可移动裂变室  <li>  'UA' 组件。带有24个AIC的二氧化铀燃料组件和一个中央可移动裂变室  <li>  'PX'组件。MOX燃料组件，带有24个导向管和一个中央可移动裂变室  <li>  'R'组件：一个反射器。   </ol>  
 * 

 * 
 * 注意这里列出的数字和从基准描述中提取的数字，以良好的老Fortran方式，是基于一的。我们以后在给各个单元分配材料时将从每个数字中减去1，以便将事情转换为C语言风格的零基索引。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_assemblies = 4; 
 *     const unsigned int assembly_materials 
 *       [n_assemblies][rods_per_assembly_x][rods_per_assembly_y] = { 
 *         {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 5, 1, 1, 5, 1, 1, 7, 1, 1, 5, 1, 1, 5, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
 *         {{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 8, 1, 1, 8, 1, 1, 7, 1, 1, 8, 1, 1, 8, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, 
 *          {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}, 
 *         {{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}, 
 *          {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
 *          {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
 *          {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
 *          {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
 *          {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
 *          {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
 *          {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
 *          {2, 3, 5, 4, 4, 5, 4, 4, 7, 4, 4, 5, 4, 4, 5, 3, 2}, 
 *          {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
 *          {2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2}, 
 *          {2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2}, 
 *          {2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2}, 
 *          {2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2}, 
 *          {2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2}, 
 *          {2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2}, 
 *          {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}}, 
 *         {{6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}, 
 *          {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6}}}; 
 * 
 * @endcode
 * 
 * 在描述了组成装配体的材料之后，我们必须指定装配体在核心中的排列。我们使用一个对称的模式，实际上只使用'UX'和'PX'装配体。
 * 

 * 
 * 
 * @code
 *     const unsigned int core[assemblies_x][assemblies_y][assemblies_z] = { 
 *       {{0}, {2}}, {{2}, {0}}}; 
 * 
 * @endcode
 * 
 * 我们现在可以为每个单元实际设置材料ID。为此，我们在所有单元中循环，查看单元中心的位置，并确定它将在哪个组件和燃料棒中。我们增加了一些检查，以确保我们计算的位置在我们必须查找材料的数组的范围内）。在循环结束时，我们相应地设置材料标识符。
 * 

 * 
 * 
 * @code
 *     for (auto &cell : coarse_grid.active_cell_iterators()) 
 *       { 
 *         const Point<dim> cell_center = cell->center(); 
 * 
 *         const unsigned int tmp_x = int(cell_center[0] / pin_pitch_x); 
 *         const unsigned int ax    = tmp_x / rods_per_assembly_x; 
 *         const unsigned int cx    = tmp_x - ax * rods_per_assembly_x; 
 * 
 *         const unsigned     tmp_y = int(cell_center[1] / pin_pitch_y); 
 *         const unsigned int ay    = tmp_y / rods_per_assembly_y; 
 *         const unsigned int cy    = tmp_y - ay * rods_per_assembly_y; 
 * 
 *         const unsigned int az = 
 *           (dim == 2 ? 0 : int(cell_center[dim - 1] / assembly_height)); 
 * 
 *         Assert(ax < assemblies_x, ExcInternalError()); 
 *         Assert(ay < assemblies_y, ExcInternalError()); 
 *         Assert(az < assemblies_z, ExcInternalError()); 
 * 
 *         Assert(core[ax][ay][az] < n_assemblies, ExcInternalError()); 
 * 
 *         Assert(cx < rods_per_assembly_x, ExcInternalError()); 
 *         Assert(cy < rods_per_assembly_y, ExcInternalError()); 
 * 
 *         cell->set_material_id(assembly_materials[core[ax][ay][az]][cx][cy] - 1); 
 *       } 
 * 
 * @endcode
 * 
 * 有了这样初始化的粗网格，我们创建适当数量的能量组对象，让它们用上面生成的粗网格初始化各自的网格。
 * 

 * 
 * 
 * @code
 *     for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *       energy_groups.emplace_back(std::make_unique<EnergyGroup<dim>>( 
 *         group, material_data, coarse_grid, fe)); 
 *     convergence_table_stream.open("convergence_table"); 
 *     convergence_table_stream.precision(12); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProblemget_total_fission_sourcecode"></a> 
 * <h5><code>NeutronDiffusionProblem::get_total_fission_source</code></h5>
 * 

 * 
 * 在特征值计算中，我们需要在每次功率迭代后计算裂变中子源总量。然后用总功率来更新k-effective。
 * 

 * 
 * 由于总裂变源是所有能量组的总和，而且每个总和都可以独立计算，所以我们实际上是以并行方式进行的。其中一个问题是， <code>EnergyGroup</code> 类中计算裂变源的函数会返回一个值。我们想在循环本身中把这些值加在一起：理想的情况是，每个任务计算它的值，然后立即把它加到总数中。以这种方式对数值进行加总需要两个功能。  <ol>  
 * <li>  我们需要一种存储数值的方式，使多个线程能够以防止数据竞赛的方式并发地读写（即线程安全的读写）。 </li>  
 * <li>  我们需要一种方法来增加这样一个值，而且是线程安全的。 </li>  
 * </ol>  
 * 

 * 
 * 第一个特性可以通过模板类实现  <code>std::atomic</code>  。然而，第二个特性，由 <code>std::atomic<double>::fetch_add()</code> 实现，只在C++20及以后的版本中可用：由于deal.II支持旧版本的C++语言标准，我们还不能使用这个特性。因此，取而代之的是，我们简单地将每个组的值写成向量中的一个条目，并在函数的最后将这些值相加。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double NeutronDiffusionProblem<dim>::get_total_fission_source() const 
 *   { 
 *     std::vector<double>  fission_sources(parameters.n_groups); 
 *     Threads::TaskGroup<> tasks; 
 *     for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *       tasks += Threads::new_task<>([&, group]() { 
 *         fission_sources[group] = energy_groups[group]->get_fission_source(); 
 *       }); 
 *     tasks.join_all(); 
 * 
 *     return std::accumulate(fission_sources.begin(), fission_sources.end(), 0.0); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProblemrefine_gridcode"></a> 
 * <h5><code>NeutronDiffusionProblem::refine_grid</code></h5>
 * 

 * 
 * 下一个函数让各个能量组对象细化其网格。这其中的大部分，也是可以独立并行完成的任务：首先，让所有的能量组对象并行计算它们的误差指标，然后计算所有能量组的最大误差指标，并确定细化和粗化单元的阈值，然后要求所有的能量组相应地细化它们的网格，也是并行的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void NeutronDiffusionProblem<dim>::refine_grid() 
 *   { 
 *     std::vector<types::global_dof_index> n_cells(parameters.n_groups); 
 *     for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *       n_cells[group] = energy_groups[group]->n_active_cells(); 
 * 
 *     BlockVector<float> group_error_indicators(n_cells); 
 * 
 *     { 
 *       Threads::TaskGroup<> tasks; 
 *       for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *         tasks += Threads::new_task([&, group]() { 
 *           energy_groups[group]->estimate_errors( 
 *             group_error_indicators.block(group)); 
 *         }); 
 *     } 
 * 
 * @endcode
 * 
 * Threads::TaskGroup 的析构器连接所有线程，所以我们知道在我们退出范围时，计算已经完成。
 * 

 * 
 * 
 * @code
 *     const float max_error         = group_error_indicators.linfty_norm(); 
 *     const float refine_threshold  = 0.3 * max_error; 
 *     const float coarsen_threshold = 0.01 * max_error; 
 * 
 *     { 
 *       Threads::TaskGroup<void> tasks; 
 *       for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *         tasks += Threads::new_task([&, group]() { 
 *           energy_groups[group]->refine_grid(group_error_indicators.block(group), 
 *                                             refine_threshold, 
 *                                             coarsen_threshold); 
 *         }); 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeNeutronDiffusionProblemruncode"></a> 
 * <h5><code>NeutronDiffusionProblem::run</code></h5>
 * 

 * 
 * 最后，这就是肉的函数：在一连串的网格上进行迭代，并对每一个网格进行幂级迭代，以计算特征值。
 * 

 * 
 * 鉴于介绍中对算法的描述，实际上没有什么可评论的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void NeutronDiffusionProblem<dim>::run() 
 *   { 
 * 
 * @endcode
 * 
 * 我们希望只为这个函数改变输出精度，并在这个函数返回时恢复到 <code>std::cout</code> 的状态。因此，我们需要一种方法来撤销输出格式的改变。Boost提供了一种方便的方法来保存输出流的状态，并在当前块结束时（当调用 <code>restore_flags</code> 的析构器时）用 <code>ios_flags_saver</code> 类来恢复它，我们在这里使用了这种方法。
 * 

 * 
 * 
 * @code
 *     boost::io::ios_flags_saver restore_flags(std::cout); 
 *     std::cout << std::setprecision(12) << std::fixed; 
 * 
 * @endcode
 * 
 * 我们通过k_eff的变化来计算下面的误差（即k_eff_old的区别。
 * 

 * 
 * 
 * @code
 *     double k_eff_old = 0.0; 
 * 
 *     for (unsigned int cycle = 0; cycle < parameters.n_refinement_cycles; 
 *          ++cycle) 
 *       { 
 * 
 * @endcode
 * 
 * 我们将在下面测量每个周期所需的CPU时间。计时器的构造函数调用 Timer::start(), ，所以一旦我们创建了一个计时器，就可以查询它的信息。由于这个循环的许多部分是用任务并行化的，所以我们测量的CPU时间（如果我们用一个以上的线程运行）将大于墙的时间。
 * 

 * 
 * 
 * @code
 *         Timer timer; 
 * 
 *         std::cout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             initialize_problem(); 
 *             for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *               energy_groups[group]->setup_linear_system(); 
 *           } 
 * 
 *         else 
 *           { 
 *             refine_grid(); 
 *             for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *               energy_groups[group]->solution *= k_eff; 
 *           } 
 * 
 *         std::cout << "   Numbers of active cells:       "; 
 *         for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *           std::cout << energy_groups[group]->n_active_cells() << ' '; 
 *         std::cout << std::endl; 
 *         std::cout << "   Numbers of degrees of freedom: "; 
 *         for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *           std::cout << energy_groups[group]->n_dofs() << ' '; 
 *         std::cout << std::endl << std::endl; 
 * 
 *         Threads::TaskGroup<> tasks; 
 *         for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *           tasks += Threads::new_task( 
 *             [&, group]() { energy_groups[group]->assemble_system_matrix(); }); 
 *         tasks.join_all(); 
 * 
 *         double       error; 
 *         unsigned int iteration = 1; 
 *         do 
 *           { 
 *             for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *               { 
 *                 energy_groups[group]->assemble_ingroup_rhs( 
 *                   Functions::ZeroFunction<dim>()); 
 * 
 *                 for (unsigned int bgroup = 0; bgroup < parameters.n_groups; 
 *                      ++bgroup) 
 *                   energy_groups[group]->assemble_cross_group_rhs( 
 *                     *energy_groups[bgroup]); 
 * 
 *                 energy_groups[group]->solve(); 
 *               } 
 * 
 *             k_eff = get_total_fission_source(); 
 *             error = std::abs(k_eff - k_eff_old) / std::abs(k_eff); 
 *             const double flux_ratio = energy_groups[0]->solution.linfty_norm() / 
 *                                       energy_groups[1]->solution.linfty_norm(); 
 *             const double max_thermal = energy_groups[1]->solution.linfty_norm(); 
 *             std::cout << "Iter number:" << std::setw(2) << std::right 
 *                       << iteration << " k_eff=" << k_eff 
 *                       << " flux_ratio=" << flux_ratio 
 *                       << " max_thermal=" << max_thermal << std::endl; 
 *             k_eff_old = k_eff; 
 * 
 *             for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *               { 
 *                 energy_groups[group]->solution_old = 
 *                   energy_groups[group]->solution; 
 *                 energy_groups[group]->solution_old /= k_eff; 
 *               } 
 * 
 *             ++iteration; 
 *           } 
 *         while ((error > parameters.convergence_tolerance) && (iteration < 500)); 
 *         convergence_table_stream << cycle << " " << energy_groups[0]->n_dofs() 
 *                                  << " " << energy_groups[1]->n_dofs() << " " 
 *                                  << k_eff << " " 
 *                                  << energy_groups[0]->solution.linfty_norm() / 
 *                                       energy_groups[1]->solution.linfty_norm() 
 *                                  << '\n'; 
 * 
 *         for (unsigned int group = 0; group < parameters.n_groups; ++group) 
 *           energy_groups[group]->output_results(cycle); 
 * 
 * @endcode
 * 
 * 打印出关于模拟的信息以及耗费的CPU时间。我们可以不先调用 Timer::cpu_time() ，而直接调用 Timer::stop() ，以获得调用该函数时的已用CPU时间。
 * 

 * 
 * 
 * @code
 *         std::cout << std::endl; 
 *         std::cout << "   Cycle=" << cycle << ", n_dofs=" 
 *                   << energy_groups[0]->n_dofs() + energy_groups[1]->n_dofs() 
 *                   << ",  k_eff=" << k_eff << ", time=" << timer.cpu_time() 
 *                   << std::endl; 
 * 
 *         std::cout << std::endl << std::endl; 
 *       } 
 *   } 
 * } // namespace Step28 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * 程序中的最后一件事在 <code>main()</code> 函数中。其结构与其他大多数教程程序一样，唯一的例外是我们在这里处理一个参数文件。 为此，我们首先看一下传递给这个函数的命令行参数：如果在命令行上没有指定输入文件，那么就使用 "project.prm"，否则就取命令行上作为第一个参数给出的文件名。
 * 

 * 
 * 有了这个，我们创建一个ParameterHandler对象，让 <code>NeutronDiffusionProblem::Parameters</code> 类声明它想在输入文件中看到的所有参数（或者，采取默认值，如果参数文件中没有列出任何参数），然后读取输入文件，要求参数对象提取数值，最后把所有东西交给 <code>NeutronDiffusionProblem</code> 类型的对象来计算特征值。
 * 

 * 
 * 
 * @code
 * int main(int argc, char **argv) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step28; 
 * 
 *       std::string filename; 
 *       if (argc < 2) 
 *         filename = "project.prm"; 
 *       else 
 *         filename = argv[1]; 
 * 
 *       const unsigned int dim = 2; 
 * 
 *       ParameterHandler parameter_handler; 
 * 
 *       NeutronDiffusionProblem<dim>::Parameters parameters; 
 *       parameters.declare_parameters(parameter_handler); 
 * 
 *       parameter_handler.parse_input(filename); 
 * 
 *       parameters.get_parameters(parameter_handler); 
 * 
 *       NeutronDiffusionProblem<dim> neutron_diffusion_problem(parameters); 
 *       neutron_diffusion_problem.run(); 
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
 * @endcode
examples/step-28/doc/results.dox



<a name="Results"></a><h1>Results</h1>


我们可以用下面的输入文件运行该程序。

@code
# Listing of Parameters
# ---------------------
# Polynomial degree of the finite element to be used
set Finite element degree     = 2


# The number of energy different groups considered
set Number of energy groups   = 2


# Inner power iterations are stopped when the change in k_eff falls below this
# tolerance
set Power iteration tolerance = 1e-12


# Number of refinement cycles to be performed
set Refinement cycles         = 12
@endcode

这个程序的输出包括控制台输出，一个名为 "convergence_table "的文件，记录网格迭代的主要结果，以及vtu格式的图形输出。

控制台的输出看起来像这样。

@code
Cycle 0:
   Numbers of active cells:       1156 1156
   Numbers of degrees of freedom: 4761 4761


Iter number: 1 k_eff=319.375676634310 flux_ratio=6.836246075630 max_thermal=1.433899030144
Iter number: 2 k_eff=0.834072546055 flux_ratio=5.204601882144 max_thermal=0.004630925876
Iter number: 3 k_eff=0.862826188043 flux_ratio=4.645051765984 max_thermal=0.005380396338
...
Iter number:69 k_eff=0.906841960370 flux_ratio=4.384056022578 max_thermal=0.008466414246
Iter number:70 k_eff=0.906841960371 flux_ratio=4.384056022583 max_thermal=0.008466414246


   Cycle=0, n_dofs=9522,  k_eff=0.906841960371, time=7.623425000000



Cycle 1:
   Numbers of active cells:       1156 2380
   Numbers of degrees of freedom: 4761 10667


Iter number: 1 k_eff=0.906838267472 flux_ratio=4.385474405125 max_thermal=0.008463675976
...


Cycle 11:
   Numbers of active cells:       11749 47074
   Numbers of degrees of freedom: 50261 204523


Iter number: 1 k_eff=0.906798057750 flux_ratio=4.384878772166 max_thermal=0.008464822382
Iter number: 2 k_eff=0.906833008185 flux_ratio=4.384868138638 max_thermal=0.008465057191
...
Iter number:32 k_eff=0.906834736550 flux_ratio=4.384846081793 max_thermal=0.008465019607
Iter number:33 k_eff=0.906834736551 flux_ratio=4.384846081798 max_thermal=0.008465019607


   Cycle=11, n_dofs=254784,  k_eff=0.906834736551, time=238.593762000000
@endcode



我们看到动力迭代在第0周期后确实收敛得更快，这是因为用上一次网格迭代的解决方案进行了初始化。收敛表 "的内容是。

@code
0 4761 4761 0.906841960371 4.38405602258
1 4761 10667 0.906837901031 4.38548908776
2 4761 18805 0.906836075928 4.3854666475
3 6629 27301 0.90683550011 4.38540458087
4 12263 48095 0.906835001796 4.38538179873
5 17501 69297 0.906834858174 4.38485382341
6 19933 78605 0.90683482406 4.38485065879
7 23979 93275 0.906834787555 4.38484837926
8 30285 117017 0.906834761604 4.38484654495
9 40087 154355 0.906834746215 4.38484608319
10 45467 179469 0.906834740155 4.38484600505
11 50261 204523 0.906834736551 4.3848460818
@endcode

列的含义是：网格迭代次数，快速能量组的自由度数，热能组的自由度数，收敛的K效应和快速通量的最大值和热能的最大值之间的比率。

网格迭代#9时，快速和热能组的网格看起来如下。

 <img width="400" src="https://www.dealii.org/images/steps/developer/step-28.grid-0.9.order2.png" alt=""> &nbsp;  <img width="400" src="https://www.dealii.org/images/steps/developer/step-28.grid-1.9.order2.png" alt="">  。

我们看到，热组的网格比快组的网格要细得多。这些网格上的解决方案是：（注：通量被归一化，总裂变源等于1）。

 <img width="400" src="https://www.dealii.org/images/steps/developer/step-28.solution-0.9.order2.png" alt="">  &nbsp;  <img width="400" src="https://www.dealii.org/images/steps/developer/step-28.solution-1.9.order2.png" alt=""> 

然后我们绘制出多项式阶数等于1、2和3的收敛数据。

 <img src="https://www.dealii.org/images/steps/developer/step-28.convergence.png" alt=""> 

估计的 "精确的 "k-effective=0.906834721253，这只是从最后一次网格迭代的多项式阶数3减去2e-10。我们看到，h-adaptive计算提供了一个代数收敛性。多项式阶数越高，网格迭代收敛的速度越快。在我们的问题中，我们需要较少的DoFs数量来实现较高的多项式阶数下的相同精度。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-28.cc"
*/
