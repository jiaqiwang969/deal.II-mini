/**
@page step_9 The step-9 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Discretizingtheadvectionequation">Discretizing the advection equation</a>
        <li><a href="#Whyisthismethodcalledstreamlinediffusion">Why is this method called "streamline diffusion"?</a>
        <li><a href="#WhyisthismethodalsocalledPetrovGalerkin">Why is this method also called "Petrov-Galerkin"?</a>
        <li><a href="#Whyisthismethodalsocalledstreamlineupwind">Why is this method also called "streamline-upwind"?</a>
        <li><a href="#Solvingthelinearsystemthatcorrespondstotheadvectionequation">Solving the linear system that corresponds to the advection equation</a>
        <li><a href="#Thetestcase">The test case</a>
        <li><a href="#Asimplerefinementcriterion">A simple refinement criterion</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdatadeclaration">Equation data declaration</a>
        <li><a href="#AdvectionProblemclassdeclaration">AdvectionProblem class declaration</a>
        <li><a href="#GradientEstimationclassdeclaration">GradientEstimation class declaration</a>
        <li><a href="#AdvectionProblemclassimplementation">AdvectionProblem class implementation</a>
        <li><a href="#GradientEstimationclassimplementation">GradientEstimation class implementation</a>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-9/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



在这个例子中，我们的目的如下。<ol>  <li>  解决平流方程  $\beta \cdot \nabla u = f$  ；  <li>  显示如果我们有一台多处理器机器，我们如何使用多线程来快速获得结果；  <li>  开发一个简单的细化准则。   </ol> 虽然第二个目的在不参考代码的情况下很难进行一般性描述，但我们将在下文中讨论其他两个目的。然后将在程序中的相关地方详细介绍多线程的使用。然而，我们将遵循 @ref threads "多处理器访问共享内存的并行计算 "文件模块中详述的关于WorkStream方法的一般讨论。




<a name="Discretizingtheadvectionequation"></a><h3>Discretizing the advection equation</h3>


在本例程序中，我们要对平流方程的解进行数值近似计算

@f[
  \beta \cdot \nabla u = f,


@f]

其中 $\beta$ 是描述平流方向和速度的矢量场（如果 $\beta=\beta(\mathbf x)$ ，它可能取决于空间变量）， $f$ 是一个源函数， $u$ 是解。该方程描述的物理过程是一个给定的流场 $\beta$ ，另一种物质随其流动，其密度或浓度由 $u$ 给出。该方程不包含这第二种物质在其载体物质内的扩散，但有源项。

很明显，在流入地，上述方程需要用边界条件来增加。

@f[
  u = g \qquad\qquad \mathrm{on}\ \partial\Omega_-,


@f]

其中 $\partial\Omega_-$ 描述了边界的流入部分，正式定义为

@f[
  \partial\Omega_-
  =
  \{{\mathbf x}\in \partial\Omega: \beta\cdot{\mathbf n}({\mathbf x}) < 0\},


@f]

和 ${\mathbf n}({\mathbf x})$ 是点 ${\mathbf x}\in\partial\Omega$ 处的域的向外法线。这个定义非常直观，因为由于 ${\mathbf n}$ 指向外侧，如果传输方向 $\beta$ 指向内侧，即在流入边界，那么与 $\beta$ 的标量乘积只能是负数。数学理论规定，我们不能在边界的流出部分提出任何边界条件。

不幸的是，上述方程不能以稳定的方式用标准的有限元方法进行求解。问题是，这个方程的解在垂直于传输方向上具有不充分的规则性：虽然它们沿 "风场" $\beta$ 定义的流线是平滑的，但它们在垂直于这个方向上可能是不连续的。这很容易理解：方程 $\beta \cdot
\nabla u = f$ 的意思实质上就是<i>rate of change of $u$ in
direction $\beta$ equals $f$</i>。但该方程对垂直方向的导数没有影响，因此，如果 $u$ 在流入边界的某一点上不连续，那么这个不连续将简单地沿着从这个边界点开始的风场流线传输。这些不连续会导致数值不稳定，使标准的连续有限元离散化不可能获得稳定的解。

解决这一困难的标准方法是  <em>  "流线-上风Petrov-Galerkin"  </em>  (SUPG)方法，有时也称为流线扩散法。对该方法的良好解释可以在  @cite elman2005  中找到。从形式上看，这种方法取代了我们从强形式推导出微分方程的弱形式的步骤。我们不是用测试函数 $v$ 乘以方程并在域上积分，而是乘以 $v + \delta \beta\cdot\nabla v$  ，其中 $\delta$ 是在（局部）网格宽度范围内选择的参数 $h$ ；通过设置 $\delta=0.1h$ 通常可以获得良好的结果。 为什么这被称为 "流线扩散 "将在下面解释；目前，让我们简单地认为这是我们如何推导出稳定的离散公式。这里 $\delta$ 的值小到足以使我们不引入过度的扩散，但大到足以使所产生的问题得到良好的解决。

使用上面定义的测试函数，问题的初始弱形式将要求找到一个函数 $u_h$ ，以便对于所有测试函数 $v_h$ ，我们有

@f[
  (\beta \cdot \nabla u_h, v_h + \delta \beta\cdot\nabla v_h)_\Omega
  =
  (f, v_h + \delta \beta\cdot\nabla v_h)_\Omega.


@f]

然而，我们希望将流入的边界条件 $u=g$ 弱化到这个问题中，这可以通过要求除了上述方程之外，我们还有

@f[
  (u_h, w_h)_{\partial\Omega_-}
  =
  (g, w_h)_{\partial\Omega_-}


@f]

为所有住在边界上的测试函数 $w_h$ ，它们来自一个合适的测试空间。事实证明，一个合适的检验函数空间恰好是 $\beta\cdot {\mathbf n}$ 乘以我们已经用于域内微分方程的检验空间中的函数 $v_h$ 的踪迹。因此，我们要求对于所有测试函数 $v_h$ ，我们有

@f[
  (u_h, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}
  =
  (g, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}.


@f]

在不试图进行论证的情况下（请再次参阅关于一般的有限元方法，特别是流线扩散方法的文献），我们可以将微分方程和边界值的方程结合在我们的稳定化问题的以下弱表述中：找到一个离散函数 $u_h$ ，使得对于所有离散测试函数 $v_h$ ，存在着

@f[
  (\beta \cdot \nabla u_h, v_h + \delta \beta\cdot\nabla v_h)_\Omega


  -
  (u_h, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}
  =
  (f, v_h + \delta \beta\cdot\nabla v_h)_\Omega


  -
  (g, \beta\cdot {\mathbf n} v_h)_{\partial\Omega_-}.


@f]




人们会认为，这将导致一个系统矩阵被倒置，其形式为

@f[
  a_{ij} =
  (\beta \cdot \nabla \varphi_i,
   \varphi_j + \delta \beta\cdot\nabla \varphi_j)_\Omega


  -
  (\varphi_i, \beta\cdot {\mathbf n} \varphi_j)_{\partial\Omega_-},


@f]

与基函数  $\varphi_i,\varphi_j$  。  然而，这是每一个数值分析师至少会遇到一次的陷阱（包括作者）：我们在这里扩大了解决方案  $u_h = \sum_i U_i \varphi_i$  ，但如果我们这样做，我们将不得不解决这个问题

@f[
  U^T A = F^T,


@f]

其中 $U$ 是膨胀系数的向量，也就是说，我们必须解决我们可能天真的预期的转置问题。

这是我们在步骤3的介绍中提出的一个观点。在那里，我们认为为了避免这种问题，应该养成总是与检验函数<i>from the left</i>相乘的习惯，而不是从右边得到正确的矩阵。为了得到我们所需要的线性系统的形式，最好将弱式改写为

@f[
  (v_h + \delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega


  -
  (\beta\cdot {\mathbf n} v_h, u_h)_{\partial\Omega_-}
  =
  (v_h + \delta \beta\cdot\nabla v_h, f)_\Omega


  -
  (\beta\cdot {\mathbf n} v_h, g)_{\partial\Omega_-}


@f]

然后得到

@f[
  a_{ij} =
  (\varphi_i + \delta \beta \cdot \nabla \varphi_i,
   \beta\cdot\nabla \varphi_j)_\Omega


  -
  (\beta\cdot {\mathbf n} \varphi_i, \varphi_j)_{\partial\Omega_-},


@f]

作为系统矩阵。我们将在程序中组装这个矩阵。




<a name="Whyisthismethodcalledstreamlinediffusion"></a><h3>Why is this method called "streamline diffusion"?</h3>


看一下上面提到的双线性形式，我们看到离散解必须满足一个方程，其中弱形式的左手边有一个域项，即

@f[
  (v_h + \delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega,


@f]

或如果我们将其拆分，则形式为

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega
  +
  (\delta \beta\cdot\nabla v_h, \beta \cdot \nabla u_h)_\Omega.


@f]

如果我们想看看这将对应于什么强的方程形式，我们需要对第二项进行积分。这就产生了下面的公式，为了简单起见，我们暂时不考虑边界项。

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega


  -
  \left(v_h, \delta \nabla \cdot \left[\beta \left(\beta \cdot \nabla
  u_h\right)\right]\right)_\Omega
  +
  \text{boundary terms}.


@f]

让我们暂时假设风场 $\beta$ 是无发散的，即 $\nabla \cdot \beta = 0$  。然后将乘积法则应用于右边方括号内项的导数，并利用发散-绿色，我们将得到以下结果。

@f[
  (v_h, \beta \cdot \nabla u_h)_\Omega


  -
  \left(v_h, \delta \left[\beta \cdot \nabla\right] \left[\beta \cdot \nabla
  \right]u_h\right)_\Omega
  +
  \text{boundary terms}.


@f]

这意味着，方程的强势形式将是这样的

@f[
  \beta \cdot \nabla u_h


  -
  \delta
  \left[\beta \cdot \nabla\right] \left[\beta \cdot \nabla
  \right] u_h.


@f]

现在要认识到的是， $\beta\cdot\nabla$ 是 <em> 方向的导数 $\beta$   </em>  。因此，如果我们用 $\beta\cdot\nabla=\frac{\partial}{\partial \beta}$ 来表示（就像我们经常用 $\mathbf n\cdot\nabla=\frac{\partial}{\partial n}$ 来表示边界处法线方向的导数一样），那么方程的强形式是

@f[
  \beta \cdot \nabla u_h


  -
  \delta
  \frac{\partial^2}{\partial\beta^2} u_h.


@f]

换句话说，测试函数的不寻常选择相当于在强形式中增加了一个项，它对应于风场方向的二阶（即扩散）微分算子 $\beta$ ，即 "流线方向"。更全面的说明还必须探讨测试函数对边界值的影响，以及为什么有必要对右手边也使用相同的测试函数，但上面的讨论可能会使人明白该方法的 "流线扩散 "名称的由来。




<a name="WhyisthismethodalsocalledPetrovGalerkin"></a><h3>Why is this method also called "Petrov-Galerkin"?</h3>


"Galerkin方法 "是指通过将方程乘以测试函数 $v$ （然后在 $\Omega$ 上进行积分）来获得弱表述，其中函数 $v$ 与解 $u$ 来自同一空间（尽管可能具有不同的边界值）。但这并不是严格意义上的必要条件。我们也可以想象从不同的函数集中选择测试函数，只要这个不同的函数集具有与原始函数集 "同样多的维度"，这样我们最终就会有与自由度同样多的独立方程（在无限维的情况下，所有这些都需要适当地定义）。利用这种可能性的方法（即以不同的方式选择测试函数集和解决方案集）被称为 "Petrov-Galerkin "方法。在目前的情况下，测试函数的形式都是 $v+\beta\cdot\nabla v$ ，其中 $v$ 是来自解集。




<a name="Whyisthismethodalsocalledstreamlineupwind"></a><h3>Why is this method also called "streamline-upwind"?</h3>


[上风方法](https://en.wikipedia.org/wiki/Upwind_scheme)在推导平流方程的稳定方案方面有很长的历史。一般来说，这个想法是，我们不是在 "这里 "看一个函数，而是在 "上游 "或 "上风 "的一小段距离看它，也就是说，"这里 "的信息最初来自哪里。这可能意味着不考虑 $u(\mathbf x)$ ，而是像 $u(\mathbf x - \delta \beta)$ 这样的。或者，等于是在整合后，我们可以评估 $u(\mathbf x)$ ，而考虑 $v$ 的下游部分。   $v(\mathbf x+\delta \beta)$  .由于各种原因，这将是很麻烦的。首先，如果 $\mathbf x + \delta \beta$ 恰好在 $\Omega$ 之外，我们将不得不定义 $v$ 应该是什么；其次，数值计算积分将更加困难，因为我们不再在同一正交点评估 $u$ 和 $v$ 。但由于我们假设 $\delta$ 很小，我们可以做一个泰勒扩展。

@f[
  v(\mathbf x + \delta \beta)
  \approx
  v(\mathbf x) + \delta \beta \cdot \nabla v(\mathbf x).


@f]

这个测试函数的形式现在看起来应该很熟悉。




<a name="Solvingthelinearsystemthatcorrespondstotheadvectionequation"></a><h3>Solving the linear system that corresponds to the advection equation</h3>


由于产生的矩阵不再是对称正定的，我们不能使用通常的共轭梯度方法（在SolverCG类中实现）来解决这个系统。相反，我们使用GMRES（Generalized Minimum RESidual）方法（在SolverGMRES中实现），它适用于我们这里的那种问题。




<a name="Thetestcase"></a><h3>The test case</h3>


对于我们将在本教程程序中解决的问题，我们使用以下领域和函数（在 $d=2$ 空间维度）。

@f{eqnarray*}
  \Omega &=& [-1,1]^d \\
  \beta({\mathbf x})
  &=&
  \left(
    \begin{array}{c}2 \\ 1+\frac 45 \sin(8\pi x)\end{array}
  \right),
  \\
  s
  &=&
  0.1,
  \\
  f({\mathbf x})
  &=&
  \left\{
    \begin{array}{ll}
        \frac 1{10 s^d} &
        \mathrm{for}\ |{\mathbf x}-{\mathbf x}_0|<s, \\
        0 & \mathrm{else},
    \end{array}
  \right.
  \qquad\qquad
  {\mathbf x}_0
  =
  \left(
    \begin{array}{c} -\frac 34 \\ -\frac 34\end{array}
  \right),
  \\
  g
  &=&
  e^{5 (1 - |\mathbf x|^2)} \sin(16\pi|\mathbf x|^2).


@f}

对于 $d>2$ ，我们扩展了 $\beta$ 和 ${\mathbf x}_0$ ，只是将上面显示的最后一个组件再重复一次。

说了这么多，下面的评论是有必要的。<ol>  <li>  平流场 $\beta$ 大致以对角线方向从左下角向右上角传输解决方案，但叠加了一个摆动结构。   <li>  右手边在流入边界条件产生的场中加入了左下角的一个圆球，然后沿着这个圆球传输。   <li>  流入边界条件施加了一个加权的正弦结构，该结构与流场一起被传输。由于 $|{\mathbf x}|\ge 1$ 在边界上，加权项从未变得非常大。   </ol> 




<a name="Asimplerefinementcriterion"></a><h3>A simple refinement criterion</h3>


在以前所有的自适应细化的例子中，我们都使用了Kelly等人首先开发的误差估计器，它为每个单元 $K$ 分配了以下指标。

@f[
  \eta_K =
  \left(
    \frac {h_K}{24}
    \int_{\partial K}
      [\partial_n u_h]^2 \; d\sigma
  \right)^{1/2},


@f]

其中 $[\partial n u_h]$ 表示跨越单元格 $K$ 的一个面的法向导数的跳变。可以证明，这个误差指标使用的是二阶导数的离散类似物，由单元大小的一个幂加权，这个幂是根据这里假定使用的线性元素调整的。

@f[
  \eta_K \approx
  C h \| \nabla^2 u \|_K,


@f]

其本身与能量准则中的误差大小有关。

在目前的情况下，这个误差指标的问题是，它假定精确的解拥有二阶导数。在某些情况下，这对于拉普拉斯问题的解来说已经是个问题了，尽管那里大多数问题允许在 $H^2$ 中求解。如果解只在 $H^1$ 中，那么二阶导数在域的某些部分（低维）是奇异的，在网格细化的情况下，误差指标不会减少。因此，该算法将不断细化这些部分周围的单元，即细化为点或线（在2d中）。

然而，对于目前的情况，解通常不在 $H^1$ 中（而且这种缺失的规律性并不是像拉普拉斯方程那样的特殊情况），所以上述的误差指标并不真正适用。因此，我们将开发一个基于梯度的离散近似的指标。尽管梯度经常不存在，但这是我们唯一可用的标准，至少在我们使用连续元素时是如此。首先，我们注意到，给定两个单元 $K$ ， $K'$ ，其中心由矢量 ${\mathbf y}_{KK'}$ 连接，我们可以对一个函数 $u$ 的方向导数进行近似，如下所示。

@f[
  \frac{{\mathbf y}_{KK'}^T}{|{\mathbf y}_{KK'}|} \nabla u
  \approx
  \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|},


@f]

其中 $u(K)$ 和 $u(K')$ 表示 $u$ 在各自单元格中心的评价。现在我们将上述近似值乘以 ${\mathbf y}_{KK'}/|{\mathbf y}_{KK'}|$ ，并对 $K$ 的所有邻居 $K'$ 求和。

@f[
  \underbrace{
    \left(\sum_{K'} \frac{{\mathbf y}_{KK'} {\mathbf y}_{KK'}^T}
                         {|{\mathbf y}_{KK'}|^2}\right)}_{=:Y}
  \nabla u
  \approx
  \sum_{K'}
  \frac{{\mathbf y}_{KK'}}{|{\mathbf y}_{KK'}|}
  \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|}.


@f]

如果连接 ${\mathbf y}_{KK'}$ 与邻居的向量 $K$ 横跨整个空间（即大致为： $K$ 在所有方向都有邻居），那么左侧表达式中括号内的项形成一个正则矩阵，我们可以将其反转，得到 $u$ 对 $K$ 的梯度的近似。

@f[
  \nabla u
  \approx
  Y^{-1}
  \left(
    \sum_{K'}
    \frac{{\mathbf y}_{KK'}}{|{\mathbf y}_{KK'}|}
    \frac{u(K') - u(K)}{|{\mathbf y}_{KK'}|}
  \right).


@f]

我们将用 $\nabla_h u(K)$ 表示右手边的近似值，我们将使用以下数量作为细化标准。

@f[
  \eta_K = h^{1+d/2} |\nabla_h u_h(K)|,


@f]

这是受以下（不严谨的）论证的启发。

@f{eqnarray*}
  \|u-u_h\|^2_{L_2}
  &\le&
  C h^2 \|\nabla u\|^2_{L_2}
\\
  &\approx&
  C
  \sum_K
  h_K^2 \|\nabla u\|^2_{L_2(K)}
\\
  &\le&
  C
  \sum_K
  h_K^2 h_K^d \|\nabla u\|^2_{L_\infty(K)}
\\
  &\approx&
  C
  \sum_K
  h_K^{2+d} |\nabla_h u_h(K)|^2


@f}




 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 就像以前的例子一样，我们必须包括几个文件，其中的含义已经讨论过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/grid/grid_out.h> 
 * 
 * @endcode
 * 
 * 下面两个文件提供了多线程程序的类和信息。在第一个文件中，声明了我们需要做并行装配的类和函数（即 <code>WorkStream</code> 命名空间）。第二个文件有一个类MultithreadInfo，可以用来查询系统中的处理器数量，这在决定启动多少个并行线程时通常很有用。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/work_stream.h> 
 * #include <deal.II/base/multithread_info.h> 
 * 
 * @endcode
 * 
 * 下一个新的include文件声明了一个基类 <code>TensorFunction</code> ，与 <code>Function</code> 类不一样，但不同的是 TensorFunction::value 返回一个张量而不是一个标量。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/tensor_function.h> 
 * 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 这是C++，因为我们想把一些输出写入磁盘。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后一步和以前的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step9 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Equationdatadeclaration"></a> 
 * <h3>Equation data declaration</h3>
 * 

 * 
 * 接下来我们声明一个描述平流场的类。当然，这是一个矢量场，有多少分量就有多少空间维度。现在我们可以使用一个从 <code>Function</code> 基类派生出来的类，就像我们在前面的例子中对边界值和系数所做的那样，但是在库中还有另一种可能性，即一个描述张量值函数的基类。这比重写 Function::value() 知道多个函数成分的方法更方便：最后我们需要一个张量，所以我们不妨直接使用一个返回张量的类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class AdvectionField : public TensorFunction<1, dim> 
 *   { 
 *   public: 
 *     virtual Tensor<1, dim> value(const Point<dim> &p) const override; 
 * 
 * @endcode
 * 
 * 在前面的例子中，我们已经在多个地方使用了抛出异常的断言。但是，我们还没有看到如何声明这种异常。这可以这样做。
 * 

 * 
 * 
 * @code
 *     DeclException2(ExcDimensionMismatch, 
 *                    unsigned int, 
 *                    unsigned int, 
 *                    << "The vector has size " << arg1 << " but should have " 
 *                    << arg2 << " elements."); 
 * 
 * @endcode
 * 
 * 语法可能看起来有点奇怪，但很合理。其格式基本如下：使用其中一个宏的名称  <code>DeclExceptionN</code>, where <code>N</code>  表示异常对象应采取的附加参数的数量。在本例中，由于我们想在两个向量的大小不同时抛出异常，我们需要两个参数，所以我们使用  <code>DeclException2</code>  。第一个参数描述了异常的名称，而下面的参数则声明了参数的数据类型。最后一个参数是一连串的输出指令，这些指令将被输送到  <code>std::cerr</code>  对象中，因此出现了奇怪的格式，前面是  <code>@<@<</code>  操作符之类的。注意，我们可以通过使用名称  <code>arg1</code> through <code>argN</code>  来访问在构造时（即在  <code>Assert</code>  调用中）传递给异常的参数，其中  <code>N</code>  是通过使用各自的宏  <code>DeclExceptionN</code>  来定义的参数数。
 * 

 * 
 * 要了解预处理器如何将这个宏扩展为实际代码，请参考异常类的文档。简而言之，这个宏调用声明并定义了一个继承自 ExceptionBase 的类  <code>ExcDimensionMismatch</code>  ，它实现了所有必要的错误输出功能。
 * 

 * 
 * 
 * @code
 *   }; 
 * 
 * @endcode
 * 
 * 下面的两个函数实现了上述的接口。第一个简单地实现了介绍中所描述的函数，而第二个使用了同样的技巧来避免调用虚拟函数，在前面的例子程序中已经介绍过了。注意第二个函数中对参数的正确大小的检查，这种检查应该始终存在于这类函数中；根据我们的经验，许多甚至大多数编程错误都是由不正确的初始化数组、不兼容的函数参数等造成的；像本例中那样使用断言可以消除许多这样的问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   Tensor<1, dim> AdvectionField<dim>::value(const Point<dim> &p) const 
 *   { 
 *     Tensor<1, dim> value; 
 *     value[0] = 2; 
 *     for (unsigned int i = 1; i < dim; ++i) 
 *       value[i] = 1 + 0.8 * std::sin(8. * numbers::PI * p[0]); 
 * 
 *     return value; 
 *   } 
 * 
 * @endcode
 * 
 * 除了平流场，我们还需要两个描述源项（  <code>right hand side</code>  ）和边界值的函数。如介绍中所述，源是一个源点附近的常数函数，我们用常数静态变量  <code>center_point</code>  表示。我们使用与我们在 step-7 示例程序中所示相同的模板技巧来设置这个中心的值。剩下的就很简单了，之前已经展示过了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *   private: 
 *     static const Point<dim> center_point; 
 *   }; 
 * 
 *   template <> 
 *   const Point<1> RightHandSide<1>::center_point = Point<1>(-0.75); 
 * 
 *   template <> 
 *   const Point<2> RightHandSide<2>::center_point = Point<2>(-0.75, -0.75); 
 * 
 *   template <> 
 *   const Point<3> RightHandSide<3>::center_point = Point<3>(-0.75, -0.75, -0.75); 
 * 
 * @endcode
 * 
 * 这里唯一的新东西是我们检查 <code>component</code> 参数的值。由于这是一个标量函数，很明显，只有当所需分量的索引为0时才有意义，所以我们断言这确实是这样的。  <code>ExcIndexRange</code> 是一个全局预定义的异常（可能是最经常使用的异常，因此我们让它成为全局的，而不是某个类的局部），它需要三个参数：超出允许范围的索引，有效范围的第一个元素和超过最后一个的元素（即又是C++标准库中经常使用的半开放区间）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> & p, 
 *                                    const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 *     const double diameter = 0.1; 
 *     return ((p - center_point).norm_square() < diameter * diameter ? 
 *               0.1 / std::pow(diameter, dim) : 
 *               0.0); 
 *   } 
 * 
 * @endcode
 * 
 * 最后是边界值，这只是从 <code>Function</code> 基类派生的另一个类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double BoundaryValues<dim>::value(const Point<dim> & p, 
 *                                     const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 * 
 *     const double sine_term = std::sin(16. * numbers::PI * p.norm_square()); 
 *     const double weight    = std::exp(5. * (1. - p.norm_square())); 
 *     return weight * sine_term; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="AdvectionProblemclassdeclaration"></a> 
 * <h3>AdvectionProblem class declaration</h3>
 * 

 * 
 * 这里是这个程序的主类。它和前面的例子中的主类非常相似，所以我们再次只对其不同之处进行评论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class AdvectionProblem 
 *   { 
 *   public: 
 *     AdvectionProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 * 
 * @endcode
 * 
 * 下一组函数将被用来组装矩阵。然而，与前面的例子不同， <code>assemble_system()</code> 函数不会自己做这些工作，而是将实际的装配工作委托给辅助函数  <code>assemble_local_system()</code>  和  <code>copy_local_to_global()</code>  。其原理是，矩阵组装可以很好地并行化，因为每个单元的局部贡献的计算完全独立于其他单元，我们只需要在将一个单元的贡献添加到全局矩阵中时进行同步。
 * 

 * 
 * 我们在这里选择的并行化策略是文档中 @ref threads 模块中详细提及的可能性之一。具体来说，我们将使用那里讨论的WorkStream方法。由于这个模块有很多文档，我们不会在这里重复设计选择的理由（例如，如果你读完上面提到的模块，你会明白 <code>AssemblyScratchData</code> 和 <code>AssemblyCopyData</code> 结构的目的是什么）。相反，我们将只讨论具体的实现。
 * 

 * 
 * 如果你阅读了上面提到的页面，你会发现为了使汇编并行化，我们需要两个数据结构--一个对应于我们在局部集成过程中需要的数据（"scratch data"，即我们只需要作为临时存储的东西），另一个是将信息从局部集成携带到函数中，然后将局部贡献添加到全局矩阵的相应元素中。其中前者通常包含FEValues和FEFaceValues对象，而后者则有局部矩阵、局部右手边，以及关于哪些自由度生活在我们正在组装局部贡献的单元上的信息。有了这些信息，下面的内容应该是相对不言自明的。
 * 

 * 
 * 
 * @code
 *     struct AssemblyScratchData 
 *     { 
 *       AssemblyScratchData(const FiniteElement<dim> &fe); 
 *       AssemblyScratchData(const AssemblyScratchData &scratch_data); 
 * 
 * @endcode
 * 
 * FEValues和FEFaceValues是很昂贵的设置对象，所以我们把它们包含在scratch对象中，以便尽可能多的数据在单元格之间被重复使用。
 * 

 * 
 * 
 * @code
 *       FEValues<dim>     fe_values; 
 *       FEFaceValues<dim> fe_face_values; 
 * 
 * @endcode
 * 
 * 我们还存储了一些向量，我们将在每个单元格上填充数值。在通常情况下，设置这些对象是很便宜的；但是，它们需要内存分配，这在多线程应用程序中可能很昂贵。因此，我们把它们保存在这里，这样在一个单元格上的计算就不需要新的分配。
 * 

 * 
 * 
 * @code
 *       std::vector<double>         rhs_values; 
 *       std::vector<Tensor<1, dim>> advection_directions; 
 *       std::vector<double>         face_boundary_values; 
 *       std::vector<Tensor<1, dim>> face_advection_directions; 
 * 
 * @endcode
 * 
 * 最后，我们需要描述该问题数据的对象。
 * 

 * 
 * 
 * @code
 *       AdvectionField<dim> advection_field; 
 *       RightHandSide<dim>  right_hand_side; 
 *       BoundaryValues<dim> boundary_values; 
 *     }; 
 * 
 *     struct AssemblyCopyData 
 *     { 
 *       FullMatrix<double>                   cell_matrix; 
 *       Vector<double>                       cell_rhs; 
 *       std::vector<types::global_dof_index> local_dof_indices; 
 *     }; 
 * 
 *     void assemble_system(); 
 *     void local_assemble_system( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       AssemblyScratchData &                                 scratch, 
 *       AssemblyCopyData &                                    copy_data); 
 *     void copy_local_to_global(const AssemblyCopyData &copy_data); 
 * 
 * @endcode
 * 
 * 下面的函数又和前面的例子一样，后面的变量也是一样的。
 * 

 * 
 * 
 * @code
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     FE_Q<dim> fe; 
 * 
 *     AffineConstraints<double> hanging_node_constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="GradientEstimationclassdeclaration"></a> 
 * <h3>GradientEstimation class declaration</h3>
 * 

 * 
 * 现在，最后，这里有一个类，它将计算每个单元上梯度的差分近似值，并以网格大小的幂数进行权衡，如介绍中所述。这个类是库中 <code>DerivativeApproximation</code> 类的一个简单版本，它使用类似的技术来获得有限元场的梯度的有限差分近似值，或者更高导数。
 * 

 * 
 * 该类有一个公共静态函数 <code>estimate</code> ，被调用来计算误差指标的向量，还有一些私有函数，在所有活动单元上做实际工作。在库的其他部分，我们遵循一个非正式的惯例，使用浮点数向量作为误差指标，而不是常见的双数向量，因为对于估计值来说，额外的精度是没有必要的。
 * 

 * 
 * 除了这两个函数，该类还声明了两个异常，当一个单元在每个空间方向上都没有邻居时（在这种情况下，介绍中描述的矩阵将是奇异的，不能被倒置），而另一个异常用于更常见的函数参数无效的情况，即一个大小错误的向量。
 * 

 * 
 * 还有两点意见：首先，这个类没有非静态成员函数或变量，所以这不是一个真正的类，而是起到了C++中 <code>namespace</code> 的作用。我们选择类而不是命名空间的原因是，这种方式我们可以声明私有的函数。如果在命名空间的头文件中声明一些函数，并在实现文件中实现这些函数和其他函数，这也可以用命名空间来实现。没有在头文件中声明的函数仍然在名字空间中，但不能从外部调用。然而，由于我们这里只有一个文件，在目前的情况下不可能隐藏函数。
 * 

 * 
 * 第二个意见是，维度模板参数被附在函数上，而不是附在类本身。这样，你就不必像其他大多数情况下那样自己指定模板参数，而是编译器可以从作为第一个参数传递的DoFHandler对象的尺寸中自行计算出其值。
 * 

 * 
 * 在开始实施之前，让我们也来评论一下并行化策略。我们已经在上面这个程序的主类的声明中介绍了使用WorkStream概念的必要框架。我们将在这里再次使用它。在目前的情况下，这意味着我们必须定义 <ol>  。
 * <li> 类，用于抓取和复制对象， </li>  。
 * <li>  一个在一个单元上进行局部计算的函数，以及 </li>  
 * <li>  一个将本地结果复制到全局对象的函数。 </li>  
 * </ol>  鉴于这个总体框架，我们将稍微偏离它。特别是，WorkStream一般是为这样的情况而发明的，即每个单元上的局部计算<i>adds</i>到一个全局对象--例如，在组装线性系统时，我们将局部贡献添加到全局矩阵和右手边中。WorkStream的设计是为了处理多个线程试图同时进行这种添加的潜在冲突，因此必须提供一些方法来确保每次只有一个线程可以做这个。然而，这里的情况略有不同：我们单独计算每个单元的贡献，但随后我们需要做的是将它们放入每个单元独有的输出向量中的一个元素。因此，不存在来自两个单元的写操作可能发生冲突的风险，也没有必要使用WorkStream的复杂机制来避免冲突的写操作。因此，我们要做的就是这样。我们仍然需要一个持有例如 FEValues 对象的 scratch 对象。但是，我们只创建一个假的、空的拷贝数据结构。同样，我们确实需要计算本地贡献的函数，但由于它已经可以把结果放到最终位置，我们不需要一个从本地到全球的拷贝函数，而是给 WorkStream::run() 函数一个空函数对象--相当于一个NULL函数指针。
 * 

 * 
 * 
 * @code
 *   class GradientEstimation 
 *   { 
 *   public: 
 *     template <int dim> 
 *     static void estimate(const DoFHandler<dim> &dof, 
 *                          const Vector<double> & solution, 
 *                          Vector<float> &        error_per_cell); 
 * 
 *     DeclException2(ExcInvalidVectorLength, 
 *                    int, 
 *                    int, 
 *                    << "Vector has length " << arg1 << ", but should have " 
 *                    << arg2); 
 *     DeclException0(ExcInsufficientDirections); 
 * 
 *   private: 
 *     template <int dim> 
 *     struct EstimateScratchData 
 *     { 
 *       EstimateScratchData(const FiniteElement<dim> &fe, 
 *                           const Vector<double> &    solution, 
 *                           Vector<float> &           error_per_cell); 
 *       EstimateScratchData(const EstimateScratchData &data); 
 * 
 *       FEValues<dim> fe_midpoint_value; 
 *       std::vector<typename DoFHandler<dim>::active_cell_iterator> 
 *         active_neighbors; 
 * 
 *       const Vector<double> &solution; 
 *       Vector<float> &       error_per_cell; 
 * 
 *       std::vector<double> cell_midpoint_value; 
 *       std::vector<double> neighbor_midpoint_value; 
 *     }; 
 * 
 *     struct EstimateCopyData 
 *     {}; 
 * 
 *     template <int dim> 
 *     static void 
 *     estimate_cell(const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *                   EstimateScratchData<dim> &scratch_data, 
 *                   const EstimateCopyData &  copy_data); 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="AdvectionProblemclassimplementation"></a> 
 * <h3>AdvectionProblem class implementation</h3>
 * 

 * 
 * 现在是主类的实现。构造器、析构器和函数 <code>setup_system</code> 遵循之前使用的模式，所以我们不需要对这三个函数进行评论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   AdvectionProblem<dim>::AdvectionProblem() 
 *     : dof_handler(triangulation) 
 *     , fe(5) 
 *   {} 
 * 
 *   template <int dim> 
 *   void AdvectionProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     hanging_node_constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, 
 *                                     dsp, 
 *                                     hanging_node_constraints, 
 *                                     /*keep_constrained_dofs =  */ false);
 * 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * 
 * @endcode
 * 
 * 在下面的函数中，矩阵和右手被组装起来。正如上面main类的文档所述，它本身并不做这个，而是委托给接下来的函数，利用 @ref threads 中讨论的WorkStream概念。
 * 

 * 
 * 如果你看了 @ref threads 模块，你会发现并行装配并不需要大量的额外代码，只要你认真地描述什么是从头开始和复制数据对象，如果你为本地装配和从本地贡献到全局对象的复制操作定义了合适的函数。完成这些工作后，下面将完成所有繁重的工作，使这些操作在多个线程上完成，只要你的系统有多少个内核。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::assemble_system() 
 *   { 
 *     WorkStream::run(dof_handler.begin_active(), 
 *                     dof_handler.end(), 
 *                     *this, 
 *                     &AdvectionProblem::local_assemble_system, 
 *                     &AdvectionProblem::copy_local_to_global, 
 *                     AssemblyScratchData(fe), 
 *                     AssemblyCopyData()); 
 *   } 
 * 
 * @endcode
 * 
 * 正如上面已经提到的，我们需要有抓取对象来进行局部贡献的并行计算。这些对象包含FEValues和FEFaceValues对象（以及一些数组），因此我们需要有构造函数和复制构造函数，以便我们能够创建它们。对于单元项，我们需要形状函数的值和梯度、正交点以确定给定点的源密度和平流场，以及正交点的权重乘以这些点的雅各布系数的行列式。相反，对于边界积分，我们不需要梯度，而是需要单元的法向量。这决定了我们必须将哪些更新标志传递给类的成员的构造函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData( 
 *     const FiniteElement<dim> &fe) 
 *     : fe_values(fe, 
 *                 QGauss<dim>(fe.degree + 1), 
 *                 update_values | update_gradients | update_quadrature_points | 
 *                   update_JxW_values) 
 *     , fe_face_values(fe, 
 *                      QGauss<dim - 1>(fe.degree + 1), 
 *                      update_values | update_quadrature_points | 
 *                        update_JxW_values | update_normal_vectors) 
 *     , rhs_values(fe_values.get_quadrature().size()) 
 *     , advection_directions(fe_values.get_quadrature().size()) 
 *     , face_boundary_values(fe_face_values.get_quadrature().size()) 
 *     , face_advection_directions(fe_face_values.get_quadrature().size()) 
 *   {} 
 * 
 *   template <int dim> 
 *   AdvectionProblem<dim>::AssemblyScratchData::AssemblyScratchData( 
 *     const AssemblyScratchData &scratch_data) 
 *     : fe_values(scratch_data.fe_values.get_fe(), 
 *                 scratch_data.fe_values.get_quadrature(), 
 *                 update_values | update_gradients | update_quadrature_points | 
 *                   update_JxW_values) 
 *     , fe_face_values(scratch_data.fe_face_values.get_fe(), 
 *                      scratch_data.fe_face_values.get_quadrature(), 
 *                      update_values | update_quadrature_points | 
 *                        update_JxW_values | update_normal_vectors) 
 *     , rhs_values(scratch_data.rhs_values.size()) 
 *     , advection_directions(scratch_data.advection_directions.size()) 
 *     , face_boundary_values(scratch_data.face_boundary_values.size()) 
 *     , face_advection_directions(scratch_data.face_advection_directions.size()) 
 *   {} 
 * 
 * @endcode
 * 
 * 现在，这就是做实际工作的函数。它与前面例子程序中的 <code>assemble_system</code> 函数没有什么不同，所以我们将再次只对其不同之处进行评论。数学上的东西紧跟我们在介绍中所说的。
 * 

 * 
 * 不过，这里有一些值得一提的地方。首先，我们把FEValues和FEFaceValues对象移到了ScratchData对象中。我们这样做是因为我们每次进入这个函数时都要简单地创建一个，也就是在每个单元格上。现在发现，FEValues类的编写目标很明确，就是将所有从单元格到单元格保持不变的东西都移到对象的构造中，每当我们移到一个新单元格时，只在 FEValues::reinit() 做尽可能少的工作。这意味着在这个函数中创建一个这样的新对象是非常昂贵的，因为我们必须为每一个单元格都这样做--这正是我们想通过FEValues类来避免的事情。相反，我们所做的是在抓取对象中只创建一次（或少数几次），然后尽可能多地重复使用它。
 * 

 * 
 * 这就引出了一个问题：我们在这个函数中创建的其他对象，与它的使用相比，其创建成本很高。事实上，在函数的顶部，我们声明了各种各样的对象。 <code>AdvectionField</code>  ,  <code>RightHandSide</code> and <code>BoundaryValues</code> 的创建成本并不高，所以这里没有什么危害。然而，在创建 <code>rhs_values</code> 和下面类似的变量时，分配内存通常要花费大量的时间，而只是访问我们存储在其中的（临时）值。因此，这些将是移入 <code>AssemblyScratchData</code> 类的候选者。我们将把这作为一个练习。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::local_assemble_system( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     AssemblyScratchData &                                 scratch_data, 
 *     AssemblyCopyData &                                    copy_data) 
 *   { 
 * 
 * @endcode
 * 
 * 我们定义一些缩写，以避免不必要的长行。
 * 

 * 
 * 
 * @code
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points = 
 *       scratch_data.fe_values.get_quadrature().size(); 
 *     const unsigned int n_face_q_points = 
 *       scratch_data.fe_face_values.get_quadrature().size(); 
 * 
 * @endcode
 * 
 * 我们声明单元格矩阵和单元格右侧...
 * 

 * 
 * 
 * @code
 *     copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
 *     copy_data.cell_rhs.reinit(dofs_per_cell); 
 * 
 * @endcode
 * 
 * ...一个数组，用于保存我们目前正在处理的单元格的自由度的全局索引...
 * 

 * 
 * 
 * @code
 *     copy_data.local_dof_indices.resize(dofs_per_cell); 
 * 
 * @endcode
 * 
 * ...然后初始化 <code>FEValues</code> 对象...
 * 

 * 
 * 
 * @code
 *     scratch_data.fe_values.reinit(cell); 
 * 
 * @endcode
 * 
 * ... 获得正交点的右手边和平流方向的数值...
 * 

 * 
 * 
 * @code
 *     scratch_data.advection_field.value_list( 
 *       scratch_data.fe_values.get_quadrature_points(), 
 *       scratch_data.advection_directions); 
 *     scratch_data.right_hand_side.value_list( 
 *       scratch_data.fe_values.get_quadrature_points(), scratch_data.rhs_values); 
 * 
 * @endcode
 * 
 * ... 设置流线扩散参数的值，如介绍中所述...
 * 

 * 
 * 
 * @code
 *     const double delta = 0.1 * cell->diameter(); 
 * 
 * @endcode
 * 
 * ...... 并按照上面的讨论，集合对系统矩阵和右手边的局部贡献。
 * 

 * 
 * 
 * @code
 *     for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *         { 
 * 
 * @endcode
 * 
 * 别名AssemblyScratchData对象，以防止行数过长。
 * 

 * 
 * 
 * @code
 *           const auto &sd = scratch_data; 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             copy_data.cell_matrix(i, j) += 
 *               ((sd.fe_values.shape_value(i, q_point) +           // (phi_i + 
 *                 delta * (sd.advection_directions[q_point] *      // delta beta 
 *                          sd.fe_values.shape_grad(i, q_point))) * // grad phi_i) 
 *                sd.advection_directions[q_point] *                // beta 
 *                sd.fe_values.shape_grad(j, q_point)) *            // grad phi_j 
 *               sd.fe_values.JxW(q_point);                         // dx 
 * 
 *           copy_data.cell_rhs(i) += 
 *             (sd.fe_values.shape_value(i, q_point) +           // (phi_i + 
 *              delta * (sd.advection_directions[q_point] *      // delta beta 
 *                       sd.fe_values.shape_grad(i, q_point))) * // grad phi_i) 
 *             sd.rhs_values[q_point] *                          // f 
 *             sd.fe_values.JxW(q_point);                        // dx 
 *         } 
 * 
 * @endcode
 * 
 * 除了我们现在建立的单元项，本问题的双线性形式还包含域的边界上的项。因此，我们必须检查这个单元的任何一个面是否在域的边界上，如果是的话，也要把这个面的贡献集合起来。当然，双线性形式只包含来自边界 <code>inflow</code> 部分的贡献，但要找出本单元的某个面是否属于流入边界的一部分，我们必须有关于正交点的确切位置和该点的流动方向的信息；我们使用FEFaceValues对象获得这些信息，并只在主循环中决定某个正交点是否在流入边界上。
 * 

 * 
 * 
 * @code
 *     for (const auto &face : cell->face_iterators()) 
 *       if (face->at_boundary()) 
 *         { 
 * 
 * @endcode
 * 
 * 好的，当前单元格的这个面是在域的边界上。就像我们在前面的例子和上面的例子中使用的通常的FEValues对象一样，我们必须重新初始化当前面的FEFaceValues对象。
 * 

 * 
 * 
 * @code
 *           scratch_data.fe_face_values.reinit(cell, face); 
 * 
 * @endcode
 * 
 * 对于手头的正交点，我们要求提供流入函数的值和流动方向。
 * 

 * 
 * 
 * @code
 *           scratch_data.boundary_values.value_list( 
 *             scratch_data.fe_face_values.get_quadrature_points(), 
 *             scratch_data.face_boundary_values); 
 *           scratch_data.advection_field.value_list( 
 *             scratch_data.fe_face_values.get_quadrature_points(), 
 *             scratch_data.face_advection_directions); 
 * 
 * @endcode
 * 
 * 现在循环所有正交点，看看这个面是在边界的流入还是流出部分。法向量指向单元外：由于该面处于边界，法向量指向域外，所以如果平流方向指向域内，其与法向量的标量乘积一定是负的（要知道为什么会这样，请考虑使用余弦的标量乘积定义）。
 * 

 * 
 * 
 * @code
 *           for (unsigned int q_point = 0; q_point < n_face_q_points; ++q_point) 
 *             if (scratch_data.fe_face_values.normal_vector(q_point) * 
 *                   scratch_data.face_advection_directions[q_point] < 
 *                 0.) 
 * 
 * @endcode
 * 
 * 如果该面是流入边界的一部分，则使用从FEFaceValues对象中获得的值和介绍中讨论的公式，计算该面对全局矩阵和右侧的贡献。
 * 

 * 
 * 
 * @code
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     copy_data.cell_matrix(i, j) -= 
 *                       (scratch_data.face_advection_directions[q_point] * 
 *                        scratch_data.fe_face_values.normal_vector(q_point) * 
 *                        scratch_data.fe_face_values.shape_value(i, q_point) * 
 *                        scratch_data.fe_face_values.shape_value(j, q_point) * 
 *                        scratch_data.fe_face_values.JxW(q_point)); 
 * 
 *                   copy_data.cell_rhs(i) -= 
 *                     (scratch_data.face_advection_directions[q_point] * 
 *                      scratch_data.fe_face_values.normal_vector(q_point) * 
 *                      scratch_data.face_boundary_values[q_point] * 
 *                      scratch_data.fe_face_values.shape_value(i, q_point) * 
 *                      scratch_data.fe_face_values.JxW(q_point)); 
 *                 } 
 *         } 
 * 
 * @endcode
 * 
 * 复制程序需要的最后一条信息是这个单元上自由度的全局索引，所以我们最后把它们写到本地数组中。
 * 

 * 
 * 
 * @code
 *     cell->get_dof_indices(copy_data.local_dof_indices); 
 *   } 
 * 
 * @endcode
 * 
 * 我们需要写的第二个函数是将前一个函数计算出的本地贡献（并放入AssemblyCopyData对象）复制到全局矩阵和右侧向量对象。这基本上就是我们在每个单元上装配东西时，一直作为最后一块代码的内容。因此，下面的内容应该是很明显的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   AdvectionProblem<dim>::copy_local_to_global(const AssemblyCopyData &copy_data) 
 *   { 
 *     hanging_node_constraints.distribute_local_to_global( 
 *       copy_data.cell_matrix, 
 *       copy_data.cell_rhs, 
 *       copy_data.local_dof_indices, 
 *       system_matrix, 
 *       system_rhs); 
 *   } 
 * 
 * @endcode
 * 
 * 这里是线性求解程序。由于系统不再像以前的例子那样是对称正定的，我们不能再使用共轭梯度法。相反，我们使用一个更通用的，不依赖矩阵的任何特殊属性的求解器：GMRES方法。GMRES和共轭梯度法一样，需要一个合适的预处理程序：我们在这里使用一个雅可比预处理程序，它对这个问题来说足够好。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::solve() 
 *   { 
 *     SolverControl               solver_control(std::max<std::size_t>(1000, 
 *                                                        system_rhs.size() / 10), 
 *                                  1e-10 * system_rhs.l2_norm()); 
 *     SolverGMRES<Vector<double>> solver(solver_control); 
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.0); 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *     Vector<double> residual(dof_handler.n_dofs()); 
 * 
 *     system_matrix.vmult(residual, solution); 
 *     residual -= system_rhs; 
 *     std::cout << "   Iterations required for convergence: " 
 *               << solver_control.last_step() << '\n' 
 *               << "   Max norm of residual:                " 
 *               << residual.linfty_norm() << '\n'; 
 * 
 *     hanging_node_constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 下面的函数根据介绍中描述的数量来细化网格。各自的计算是在类  <code>GradientEstimation</code>  中进行的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::refine_grid() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     GradientEstimation::estimate(dof_handler, 
 *                                  solution, 
 *                                  estimated_error_per_cell); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.03); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 * @endcode
 * 
 * 这个函数与第6步中的函数类似，但由于我们使用的是高阶有限元，所以我们以不同的方式保存解决方案。像VisIt和Paraview这样的可视化程序通常只能理解与节点相关的数据：它们不能绘制五度基函数，这导致我们计算的解的图片非常不准确。为了解决这个问题，我们为每个单元保存了多个 <em> 补丁 </em> ：在二维中，我们为每个单元在VTU文件中保存64个双线性 "单元"，在三维中，我们保存512个。最终的结果是，可视化程序将使用立方体基础函数的片状线性插值：这捕捉到了解决方案的细节，并且在大多数屏幕分辨率下，看起来很平滑。我们在一个单独的步骤中保存网格，没有额外的补丁，这样我们就有了细胞面的视觉表现。
 * 

 * 
 * 9.1版本的deal.II获得了编写更高程度多项式（即为我们的片状二项式解决方案编写片状二项式可视化数据）VTK和VTU输出的能力：然而，并非所有最新版本的ParaView和Viscit（截至2018年）都能读取这种格式，所以我们在这里使用更古老、更通用（但效率较低）的方法。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     { 
 *       GridOut       grid_out; 
 *       std::ofstream output("grid-" + std::to_string(cycle) + ".vtu"); 
 *       grid_out.write_vtu(triangulation, output); 
 *     } 
 * 
 *     { 
 *       DataOut<dim> data_out; 
 *       data_out.attach_dof_handler(dof_handler); 
 *       data_out.add_data_vector(solution, "solution"); 
 *       data_out.build_patches(8); 
 * 
 * @endcode
 * 
 * VTU输出可能很昂贵，无论是计算还是写入磁盘。这里我们要求ZLib，一个压缩库，以最大限度地提高吞吐量的方式来压缩数据。
 * 

 * 
 * 
 * @code
 *       DataOutBase::VtkFlags vtk_flags; 
 *       vtk_flags.compression_level = 
 *         DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
 *       data_out.set_flags(vtk_flags); 
 * 
 *       std::ofstream output("solution-" + std::to_string(cycle) + ".vtu"); 
 *       data_out.write_vtu(output); 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * ... 如同主循环（设置-求解-细化）一样，除了循环次数和初始网格之外。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 10; ++cycle) 
 *       { 
 *         std::cout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation, -1, 1); 
 *             triangulation.refine_global(3); 
 *           } 
 *         else 
 *           { 
 *             refine_grid(); 
 *           } 
 * 
 *         std::cout << "   Number of active cells:              " 
 *                   << triangulation.n_active_cells() << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         std::cout << "   Number of degrees of freedom:        " 
 *                   << dof_handler.n_dofs() << std::endl; 
 * 
 *         assemble_system(); 
 *         solve(); 
 *         output_results(cycle); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="GradientEstimationclassimplementation"></a> 
 * <h3>GradientEstimation class implementation</h3>
 * 

 * 
 * 现在是 <code>GradientEstimation</code> 类的实现。让我们先为 <code>estimate_cell()</code> 函数所使用的 <code>EstimateScratchData</code> 类定义构造函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   GradientEstimation::EstimateScratchData<dim>::EstimateScratchData( 
 *     const FiniteElement<dim> &fe, 
 *     const Vector<double> &    solution, 
 *     Vector<float> &           error_per_cell) 
 *     : fe_midpoint_value(fe, 
 *                         QMidpoint<dim>(), 
 *                         update_values | update_quadrature_points) 
 *     , solution(solution) 
 *     , error_per_cell(error_per_cell) 
 *     , cell_midpoint_value(1)
 *     , neighbor_midpoint_value(1)
 *  
 *   {
 *  
 * @endcode
 * 
 * 我们分配一个向量来保存一个单元的所有活动邻居的迭代器。我们保留活动邻居的最大数量，以避免以后的重新分配。注意这个最大的活动邻居数是如何计算出来的。
 * 

 * 
 * 
 * @code
 *     active_neighbors.reserve(GeometryInfo<dim>::faces_per_cell * 
 *                              GeometryInfo<dim>::max_children_per_face); 
 *   } 
 * 
 *   template <int dim> 
 *   GradientEstimation::EstimateScratchData<dim>::EstimateScratchData( 
 *     const EstimateScratchData &scratch_data) 
 *     : fe_midpoint_value(scratch_data.fe_midpoint_value.get_fe(), 
 *                         scratch_data.fe_midpoint_value.get_quadrature(), 
 *                         update_values | update_quadrature_points) 
 *     , solution(scratch_data.solution) 
 *     , error_per_cell(scratch_data.error_per_cell) 
 *     , cell_midpoint_value(1) 
 *     , neighbor_midpoint_value(1) 
 *   {} 
 * 
 * @endcode
 * 
 * 接下来是对 <code>GradientEstimation</code> 类的实现。第一个函数除了将工作委托给另一个函数外，并没有做什么，但在顶部有一点设置。
 * 

 * 
 * 在开始工作之前，我们要检查写入结果的向量是否有正确的大小。在编程中，忘记在调用处正确确定参数大小的错误是很常见的。因为没有发现这种错误所造成的损失往往是微妙的（例如，内存中某个地方的数据损坏，或者是无法重现的结果），所以非常值得努力去检查这些东西。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void GradientEstimation::estimate(const DoFHandler<dim> &dof_handler, 
 *                                     const Vector<double> & solution, 
 *                                     Vector<float> &        error_per_cell) 
 *   { 
 *     Assert( 
 *       error_per_cell.size() == dof_handler.get_triangulation().n_active_cells(), 
 *       ExcInvalidVectorLength(error_per_cell.size(), 
 *                              dof_handler.get_triangulation().n_active_cells())); 
 * 
 *     WorkStream::run(dof_handler.begin_active(), 
 *                     dof_handler.end(), 
 *                     &GradientEstimation::template estimate_cell<dim>, 
 *                     std::function<void(const EstimateCopyData &)>(), 
 *                     EstimateScratchData<dim>(dof_handler.get_fe(), 
 *                                              solution, 
 *                                              error_per_cell), 
 *                     EstimateCopyData()); 
 *   } 
 * 
 * @endcode
 * 
 * 这里是通过计算梯度的有限差分近似值来估计局部误差的函数。该函数首先计算当前单元的活动邻居列表，然后为每个邻居计算介绍中描述的数量。之所以有这样的顺序，是因为在局部细化网格的情况下，要找到一个给定的邻居并不是一蹴而就的事情。原则上，一个优化的实现可以在一个步骤中找到邻域和取决于它们的量，而不是先建立一个邻域列表，然后在第二步中找到它们的贡献，但是我们很乐意将此作为一个练习。正如之前所讨论的，传递给 WorkStream::run 的工作者函数是在保留所有临时对象的 "scratch "对象上工作。这样，我们就不需要在每次为给定单元调用工作的函数内创建和初始化那些昂贵的对象了。这样的参数被作为第二个参数传递。第三个参数是一个 "copy-data "对象（更多信息见 @ref threads ），但我们在这里实际上没有使用这些对象。由于 WorkStream::run() 坚持传递三个参数，我们声明这个函数有三个参数，但简单地忽略了最后一个参数。
 * 

 * 
 * （从美学角度看，这是不令人满意的。它可以通过使用一个匿名（lambda）函数来避免。如果你允许的话，让我们在这里展示一下如何做。首先，假设我们已经声明这个函数只接受两个参数，省略了未使用的最后一个参数。现在， WorkStream::run 仍然想用三个参数来调用这个函数，所以我们需要找到一种方法来 "忘记 "调用中的第三个参数。简单地像上面那样把指针传给 WorkStream::run 这个函数是做不到的--编译器会抱怨一个声明为有两个参数的函数在调用时有三个参数。然而，我们可以通过将以下内容作为第三个参数传递给 WorkStream::run(): 来做到这一点 
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *  [](const typename DoFHandler<dim>::active_cell_iterator &cell,
 *     EstimateScratchData<dim> &                            scratch_data,
 *     EstimateCopyData &)
 *  {
 *    GradientEstimation::estimate_cell<dim>(cell, scratch_data);
 *  }
 *  @endcode 
 * </div>
 * 这并不比下面实现的解决方案好多少：要么例程本身必须带三个参数，要么它必须被带三个参数的东西包起来。我们不使用这种方法，因为在开始时添加未使用的参数更简单。
 * 

 * 
 * 现在来看看细节。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void GradientEstimation::estimate_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     EstimateScratchData<dim> &                            scratch_data, 
 *     const EstimateCopyData &) 
 *   { 
 * 
 * @endcode
 * 
 * 我们需要为张量 <code>Y</code> 提供空间，它是Y向量的外积之和。
 * 

 * 
 * 
 * @code
 *     Tensor<2, dim> Y; 
 * 
 * @endcode
 * 
 * 首先初始化  <code>FEValues</code>  对象，以及  <code>Y</code>  张量。
 * 

 * 
 * 
 * @code
 *     scratch_data.fe_midpoint_value.reinit(cell); 
 * 
 * @endcode
 * 
 * 现在，在我们继续之前，我们首先计算当前单元的所有活动邻居的列表。我们首先在所有面上进行循环，看那里的邻居是否处于活动状态，如果它与本单元在同一级别或更粗一级，就会出现这种情况（注意，一个邻居只能比本单元粗一次，因为我们在deal.II中只允许在一个面上有一个最大的细化差）。另外，邻居也可能在同一级别，并被进一步细化；那么我们必须找到它的哪些子单元与当前单元相邻，并选择这些子单元（注意，如果一个活动单元的邻居的一个子单元与这个活动单元相邻，那么它本身就必须是活动的，这是由于上面提到的一个细化规则）。
 * 

 * 
 * 在一个空间维度上，情况略有不同，因为在那里不存在单一细化规则：相邻的活动单元可以在任意多的细化级别上有所不同。在这种情况下，计算变得有点困难，但我们将在下面解释。
 * 

 * 
 * 在开始对当前单元的所有邻域进行循环之前，我们当然要清除存储活动邻域的迭代器的数组。
 * 

 * 
 * 
 * @code
 *     scratch_data.active_neighbors.clear(); 
 *     for (const auto face_n : cell->face_indices()) 
 *       if (!cell->at_boundary(face_n)) 
 *         { 
 * 
 * @endcode
 * 
 * 首先定义面的迭代器和邻居的缩写
 * 

 * 
 * 
 * @code
 *           const auto face     = cell->face(face_n); 
 *           const auto neighbor = cell->neighbor(face_n); 
 * 
 * @endcode
 * 
 * 然后检查邻居是否是活动的。如果是，那么它就在同一层或更粗的一层（如果我们不是在1D中），而且我们在任何情况下都会对它感兴趣。
 * 

 * 
 * 
 * @code
 *           if (neighbor->is_active()) 
 *             scratch_data.active_neighbors.push_back(neighbor); 
 *           else 
 *             { 
 * 
 * @endcode
 * 
 * 如果邻居没有活动，则检查其子女。
 * 

 * 
 * 
 * @code
 *               if (dim == 1) 
 *                 { 
 * 
 * @endcode
 * 
 * 要找到与本单元相邻的子单元，如果我们在本单元的左边（n==0），则依次去找其右边的子单元，如果我们在右边（n==1），则依次去找左边的子单元，直到找到一个活动单元。
 * 

 * 
 * 
 * @code
 *                   auto neighbor_child = neighbor; 
 *                   while (neighbor_child->has_children()) 
 *                     neighbor_child = neighbor_child->child(face_n == 0 ? 1 : 0); 
 * 
 * @endcode
 * 
 * 由于这使用了一些非微妙的几何直觉，我们可能想检查一下我们是否做对了，也就是说，检查我们找到的单元格的邻居是否确实是我们目前正在处理的单元。像这样的检查通常是有用的，并且经常发现像上面这一行的算法（不由自主地交换 <code>n==1</code> for <code>n==0</code> 或类似的算法是很简单的）和库中的错误（上面的算法所依据的假设可能是错误的，记录错误，或者由于库中的错误而被违反）。原则上，我们可以在程序运行一段时间后删除这样的检查，但是无论如何留下它来检查库中或上述算法中的变化可能是一件好事。                请注意，如果这个检查失败了，那么这肯定是一个无法恢复的错误，而且很可能被称为内部错误。因此我们在这里使用一个预定义的异常类来抛出。
 * 

 * 
 * 
 * @code
 *                   Assert(neighbor_child->neighbor(face_n == 0 ? 1 : 0) == cell, 
 *                          ExcInternalError()); 
 * 
 * @endcode
 * 
 * 如果检查成功，我们就把刚刚发现的活动邻居推到我们保留的堆栈中。
 * 

 * 
 * 
 * @code
 *                   scratch_data.active_neighbors.push_back(neighbor_child); 
 *                 } 
 *               else 
 * 
 * @endcode
 * 
 * 如果我们不在1d中，我们收集所有 "在 "当前面的子面后面的邻居孩子，然后继续前进。
 * 

 * 
 * 
 * @code
 *                 for (unsigned int subface_n = 0; subface_n < face->n_children(); 
 *                      ++subface_n) 
 *                   scratch_data.active_neighbors.push_back( 
 *                     cell->neighbor_child_on_subface(face_n, subface_n)); 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 好了，现在我们有了所有的邻居，让我们开始对他们每个人进行计算。首先，我们做一些预备工作：找出当前单元格的中心和该点的解决方案。后者是以正交点的函数值向量的形式得到的，当然，正交点只有一个。同样地，中心的位置是实空间中第一个（也是唯一的）正交点的位置。
 * 

 * 
 * 
 * @code
 *     const Point<dim> this_center = 
 *       scratch_data.fe_midpoint_value.quadrature_point(0); 
 * 
 *     scratch_data.fe_midpoint_value.get_function_values( 
 *       scratch_data.solution, scratch_data.cell_midpoint_value); 
 * 
 * @endcode
 * 
 * 现在在所有活动邻居上循环，收集我们需要的数据。
 * 

 * 
 * 
 * @code
 *     Tensor<1, dim> projected_gradient; 
 *     for (const auto &neighbor : scratch_data.active_neighbors) 
 *       { 
 * 
 * @endcode
 * 
 * 然后得到邻近单元的中心和该点的有限元函数值。注意，为了获得这些信息，我们必须重新初始化相邻单元的 <code>FEValues</code> 对象。
 * 

 * 
 * 
 * @code
 *         scratch_data.fe_midpoint_value.reinit(neighbor); 
 *         const Point<dim> neighbor_center = 
 *           scratch_data.fe_midpoint_value.quadrature_point(0); 
 * 
 *         scratch_data.fe_midpoint_value.get_function_values( 
 *           scratch_data.solution, scratch_data.neighbor_midpoint_value); 
 * 
 * @endcode
 * 
 * 计算连接两个单元格中心的向量 <code>y</code> 。注意，与介绍不同，我们用 <code>y</code> 表示归一化的差分向量，因为这是在计算中随处可见的数量。
 * 

 * 
 * 
 * @code
 *         Tensor<1, dim> y        = neighbor_center - this_center; 
 *         const double   distance = y.norm(); 
 *         y /= distance; 
 * 
 * @endcode
 * 
 * 然后把这个单元格对Y矩阵的贡献加起来...
 * 

 * 
 * 
 * @code
 *         for (unsigned int i = 0; i < dim; ++i) 
 *           for (unsigned int j = 0; j < dim; ++j) 
 *             Y[i][j] += y[i] * y[j]; 
 * 
 * @endcode
 * 
 * ...并更新差额商数之和。
 * 

 * 
 * 
 * @code
 *         projected_gradient += (scratch_data.neighbor_midpoint_value[0] - 
 *                                scratch_data.cell_midpoint_value[0]) / 
 *                               distance * y; 
 *       } 
 * 
 * @endcode
 * 
 * 如果现在，在收集了来自邻居的所有信息后，我们可以确定当前单元的梯度的近似值，那么我们需要经过跨越整个空间的向量 <code>y</code> ，否则我们就不会有梯度的所有成分。这可以通过矩阵的可逆性来说明。
 * 

 * 
 * 如果矩阵不可逆，那么当前单元的活动邻居数量不足。与之前所有的情况（我们提出了异常）相比，这不是一个编程错误：这是一个运行时错误，即使在调试模式下运行良好，也可能在优化模式下发生，所以在优化模式下尝试捕捉这个错误是合理的。对于这种情况，有一个 <code>AssertThrow</code> 宏：它像 <code>Assert</code> 宏一样检查条件，但不仅仅是在调试模式下；然后输出一个错误信息，但不是像 <code>Assert</code> 宏那样中止程序，而是使用C++的 <code>throw</code> 命令抛出异常。这样，人们就有可能捕捉到这个错误，并采取合理的应对措施。其中一个措施是在全局范围内细化网格，因为如果初始网格的每个单元都至少被细化过一次，就不会出现方向不足的情况。
 * 

 * 
 * 
 * @code
 *     AssertThrow(determinant(Y) != 0, ExcInsufficientDirections()); 
 * 
 * @endcode
 * 
 * 如果另一方面，矩阵是可反转的，那么就反转它，用它乘以其他数量，然后用这个数量和正确的网格宽度的幂来计算估计误差。
 * 

 * 
 * 
 * @code
 *     const Tensor<2, dim> Y_inverse = invert(Y); 
 * 
 *     const Tensor<1, dim> gradient = Y_inverse * projected_gradient; 
 * 
 * @endcode
 * 
 * 这个函数的最后一部分是将我们刚刚计算出来的内容写入输出向量的元素中。这个向量的地址已经存储在Scratch数据对象中，我们所要做的就是知道如何在这个向量中获得正确的元素--但我们可以问一下我们所在的单元格是第多少个活动单元。
 * 

 * 
 * 
 * @code
 *     scratch_data.error_per_cell(cell->active_cell_index()) = 
 *       (std::pow(cell->diameter(), 1 + 1.0 * dim / 2) * gradient.norm()); 
 *   } 
 * } // namespace Step9 
 * @endcode
 * 
 * 
 * <a name="Mainfunction"></a> 
 * <h3>Main function</h3>
 * 

 * 
 * <code>main</code> 函数与前面的例子类似。主要区别是我们使用MultithreadInfo来设置最大的线程数（更多信息请参见文档模块  @ref threads  "多处理器访问共享内存的并行计算"）。使用的线程数是环境变量DEAL_II_NUM_THREADS和  <code>set_thread_limit</code>  的参数的最小值。如果没有给  <code>set_thread_limit</code>  的值，则使用英特尔线程构建块（TBB）库的默认值。如果省略了对  <code>set_thread_limit</code>  的调用，线程的数量将由 TBB 选择，与 DEAL_II_NUM_THREADS无关。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   using namespace dealii; 
 *   try 
 *     {
 *       MultithreadInfo::set_thread_limit(); 
 * 
 *       Step9::AdvectionProblem<2> advection_problem_2d; 
 *       advection_problem_2d.run(); 
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
 * 
 * @endcode
examples/step-9/doc/results.dox



<a name="Results"></a><h1>Results</h1>



这个程序的结果并不特别引人注目。它们由控制台输出、一些网格文件和每个网格的解决方案组成。首先是控制台的输出。

@code
Cycle 0:
   Number of active cells:              64
   Number of degrees of freedom:        1681
   Iterations required for convergence: 298
   Max norm of residual:                3.60316e-12
Cycle 1:
   Number of active cells:              124
   Number of degrees of freedom:        3537
   Iterations required for convergence: 415
   Max norm of residual:                3.70682e-12
Cycle 2:
   Number of active cells:              247
   Number of degrees of freedom:        6734
   Iterations required for convergence: 543
   Max norm of residual:                7.19716e-13
Cycle 3:
   Number of active cells:              502
   Number of degrees of freedom:        14105
   Iterations required for convergence: 666
   Max norm of residual:                3.45628e-13
Cycle 4:
   Number of active cells:              1003
   Number of degrees of freedom:        27462
   Iterations required for convergence: 1064
   Max norm of residual:                1.86495e-13
Cycle 5:
   Number of active cells:              1993
   Number of degrees of freedom:        55044
   Iterations required for convergence: 1251
   Max norm of residual:                1.28765e-13
Cycle 6:
   Number of active cells:              3985
   Number of degrees of freedom:        108492
   Iterations required for convergence: 2035
   Max norm of residual:                6.78085e-14
Cycle 7:
   Number of active cells:              7747
   Number of degrees of freedom:        210612
   Iterations required for convergence: 2187
   Max norm of residual:                2.61457e-14
Cycle 8:
   Number of active cells:              15067
   Number of degrees of freedom:        406907
   Iterations required for convergence: 3079
   Max norm of residual:                2.9932e-14
Cycle 9:
   Number of active cells:              29341
   Number of degrees of freedom:        780591
   Iterations required for convergence: 3913
   Max norm of residual:                8.15689e-15
@endcode



相当多的单元格被用在最精细的层面上，以解决解决方案的特征。下面是第四和第十个网格。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-grid-3.png" alt="细化周期中的第四个网格，显示对特征的一些适应性。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-grid-9.png" alt="细化周期中的第十个网格，显示完全捕捉到波。" width="400" height="400"> </div> <div> 以及第四和第十的解决方案。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-3.png" alt="第四个解决方案，显示我们解决了大多数特征，但有些仍然没有解决，显得很模糊。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-9.png" alt="第十个解决方案，显示完全解决的流动。" width="400" height="400"> </div> </div> 以及网格和解决方案都被放大。<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-3-zoom.png" alt="第四个解决方案的细节，显示我们解决了大多数特征，但有些仍然没有解决，显得模糊不清。特别是，较大的单元格需要细化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-9-solution-9-zoom.png" alt="第十个解决方案的细节，显示我们需要比第四个解决方案中存在的更多单元格。" width="400" height="400"> </div> </div>

解决方案是由那部分沿摆动的平流场从左、下边界传送到右上方的部分，以及由左下角的源所产生的部分，其结果也是沿传送的。上面显示的网格很好地适应了解决这些特征。图中的比较表明，即使我们使用的是高阶近似，我们仍然需要自适应的网格细化来完全解决摆动。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-9.cc"
*/
