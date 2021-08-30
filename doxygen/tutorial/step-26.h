/**
@page step_26 The step-26 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Adaptingmeshesfortimedependentproblems"> Adapting meshes for time dependent problems </a>
        <li><a href="#WhatcouldpossiblygowrongVerifyingwhetherthecodeiscorrect"> What could possibly go wrong? Verifying whether the code is correct </a>
        <li><a href="#Thetestcase"> The testcase </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#ThecodeHeatEquationcodeclass">The <code>HeatEquation</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ThecodeHeatEquationcodeimplementation">The <code>HeatEquation</code> implementation</a>
      <ul>
        <li><a href="#codeHeatEquationsetup_systemcode"><code>HeatEquation::setup_system</code></a>
        <li><a href="#codeHeatEquationsolve_time_stepcode"><code>HeatEquation::solve_time_step</code></a>
        <li><a href="#codeHeatEquationoutput_resultscode"><code>HeatEquation::output_results</code></a>
        <li><a href="#codeHeatEquationrefine_meshcode"><code>HeatEquation::refine_mesh</code></a>
        <li><a href="#codeHeatEquationruncode"><code>HeatEquation::run</code></a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Adaptivetimestepping">Adaptive time stepping</a>
        <li><a href="#Bettertimesteppingmethods">Better time stepping methods</a>
        <li><a href="#Betterrefinementcriteria">Better refinement criteria</a>
        <li><a href="#Positivitypreservation">Positivity preservation</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-26/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{29,30}  （  @dealiiVideoLectureSeeAlso{31.7})  。


这个程序实现了热方程

@f{align*}
  \frac{\partial u(\mathbf x, t)}{\partial t}


  -
  \Delta u(\mathbf x, t)
  &=
  f(\mathbf x, t),
  \qquad\qquad &&
  \forall \mathbf x \in \Omega, t\in (0,T),
  \\
  u(\mathbf x, 0) &= u_0(\mathbf x) &&
  \forall \mathbf x \in \Omega, \\
  \\
  u(\mathbf x, t) &= g(\mathbf x,t) &&
  \forall \mathbf x \in \partial\Omega, t \in (0,T).


@f}

从某种意义上说，这个方程比我们在前面的程序步骤23、步骤24、步骤25中讨论的方程，即波浪方程要简单。这是由于热方程随着时间的推移使解变得平滑，因此在许多方面更宽容。例如，当使用隐式时间步长方法时，我们实际上可以采取大的时间步长，我们对通过每隔几步适应网格而引入的小干扰的麻烦较小，等等。

我们在这里的目标将是使用theta-scheme解决上述方程，该方程在时间上离散，使用以下方法，我们希望 $u^n(\mathbf x)$ 在某个时间 $t_n$ 近似 $u(\mathbf x, t_n)$  。

@f{align*}
  \frac{u^n(\mathbf x)-u^{n-1}(\mathbf x)}{k_n}


  -
  \left[
  (1-\theta)\Delta u^{n-1}(\mathbf x)
  +
  \theta\Delta u^n(\mathbf x)
  \right]
  &=
  \left[
  (1-\theta)f(\mathbf x, t_{n-1})
  +
  \theta f(\mathbf x, t_n)
  \right].


@f}

这里， $k_n=t_n-t_{n-1}$ 是时间步长。Theta-scheme概括了显式欧拉（  $\theta=0$  ）、隐式欧拉（  $\theta=1$  ）和Crank-Nicolson（  $\theta=\frac 12$  ）时间离散。由于后者具有最高的收敛顺序，我们将在下面的程序中选择 $\theta=\frac 12$ ，但要使这个参数的操作保持简单。如果你对玩更高阶的方法感兴趣，可以看一下步骤52）。

鉴于这种时间离散化，空间离散化会像往常一样发生，通过与测试函数相乘，通过部分积分，然后将一切限制在一个有限维的子空间。在与 $k_n$ 相乘之后，这就产生了以下一组完全离散的方程。

@f{align*}
  M U^n-MU^{n-1}
  +
  k_n \left[
  (1-\theta)A U^{n-1}
  +
  \theta A U^n
  \right]
  &=
  k_n
  \left[
  (1-\theta)F^{n-1}
  +
  \theta F^n
  \right],


@f}

其中 $M$ 是质量矩阵， $A$ 是将拉普拉斯离散化后的刚度矩阵。将所有已知的量带到右手边，就得到了我们在每一步都要解决的线性系统。

@f{align*}
  (M
  +
  k_n \theta A) U^n
  &=
  MU^{n-1}


  -
  k_n
  (1-\theta)A U^{n-1}
  +
  k_n
  \left[
  (1-\theta)F^{n-1}
  +
  \theta F^n
  \right].


@f}

左手边的线性系统是对称和正定的，所以我们用共轭梯度法解决它应该没有问题。

如果我们在初始时间有一组节点系数 $U^0$ ，我们可以开始上面的迭代。在这里，我们采取的是将初始值 $u_0(\mathbf x)$ 插值到用于第一个时间步长的网格上得到的那些。我们还需要选择一个时间步长；在这里我们只选择固定的时间步长，但显然先进的模拟器会希望自适应地选择它。我们将在<a href="#Results">results section
below</a>中简要地回到这个问题。




<a name="Adaptingmeshesfortimedependentproblems"></a><h3> Adapting meshes for time dependent problems </h3>


在前面几个程序中求解波浪方程及其变体时，我们保持了固定的网格。就像静止方程一样，我们可以很好地说明这不是最聪明的方法，通过调整网格可以节省大量资金。然而，与静止的情况相比，还有很大的困难。让我们依次来看看这些困难。

 <ul>   <li>  <i>Time step size and minimal mesh size</i> 。对于静止的问题，一般的做法是 "将网格做得越细越好"。对于有奇点的问题，这往往会导致我们在角落或界面上得到许多细化层次的情况。第一个使用自适应网格的教程，第6步，已经是一个案例的要点。

  然而，对于时间相关问题，我们通常需要选择与网格大小相关的时间步长。对于显式时间离散，这是显而易见的，因为我们需要遵守一个CFL条件，将时间步长与最小的网格尺寸联系起来。对于隐式时间离散，不存在这样的硬性限制，但在实践中，如果我们使网格尺寸变小，我们仍然希望使时间步长变小，因为我们通常有 $\|e\| \le {\cal O}(k^p + h^q)$ 形式的误差估计，其中 $p,q$ 分别是时间和空间离散的收敛阶。我们只有减少这两个项，才能使误差变小。理想情况下，这样的估计会建议选择 $k \propto h^{q/p}$  。因为，至少对于非光滑解的问题来说，误差通常集中在网格尺寸最小的单元中，所以我们必须确实选择 $k \propto h_{\text{min}}^{q/p}$  ，使用<i>smallest</i>网格尺寸。

  其结果是，在一个地方进一步细化网格不仅意味着略微增加自由度数量的适度额外努力，而且由于时间步长较小，必须更频繁地求解<i>global</i>线性系统的更大努力。

  在实践中，我们通常通过承认我们不能使时间步长任意地小，因此也不能使局部网格尺寸任意地小来处理这个问题。相反，我们设置了一个最大的细化水平，当我们标记单元进行细化时，我们只是不细化那些子单元会超过这个最大的细化水平。

  还有一个类似的问题是，我们将选择一个右手边，在不同的时间在领域的不同部分开启。为了避免在突然需要更细的网格的地方，被太粗的网格弄得措手不及，我们还将在程序中强制执行<i>minimal</i>的网格细化水平。

    <li>  <i>Test functions from different meshes</i>。让我们再次考虑我们在上面写下的半离散方程。   @f{align*}
    \frac{u^n(\mathbf x)-u^{n-1}(\mathbf x)}{k_n}


    -
    \left[
    (1-\theta)\Delta u^{n-1}(\mathbf x)
    +
    \theta\Delta u^n(\mathbf x)
    \right]
    &=
    \left[
    (1-\theta)f(\mathbf x, t_{n-1})
    +
    \theta f(\mathbf x, t_n)
    \right].
  @f}

  我们在这里可以把 $u^{n-1}$ 视为数据，因为它可能已经被计算过了。现在，让我们替换掉@f{align*}
    u^n(\mathbf x)\approx u_h^n(\mathbf x)
    =
    \sum_j U^n \varphi_j(\mathbf x),
  @f}。

  与测试函数 $\varphi_i(\mathbf x)$ 相乘，并在必要时按部分整合。在如上所述的过程中，这将产生@f{align*}
    \sum_j
    (M
    +
    k_n \theta A)_{ij} U^n_j
    &=
    (\varphi_i, u_h^{n-1})


    -
    k_n
    (1-\theta)(\nabla \varphi_i, \nabla u_h^{n-1})
    +
    k_n
    \left[
    (1-\theta)F^{n-1}
    +
    \theta F^n
    \right].
  @f}。

  现在想象一下，我们在时间步骤  $n-1$  和  $n$  之间改变了网格。那么问题来了，我们在 $u_h^n$ 和 $u^{n-1}$ 中使用的基函数是不同的！这与右边的项有关。这与右手边的项有关，其中第一个项我们可以更清楚地写成（第二个项遵循同样的模式）@f{align*}
    (\varphi_i, u_h^{n-1})
    =
    (\varphi_i^n, u_h^{n-1})
    =
    \sum_{j=1}^{N_{n-1}}
    (\varphi_i^n, \varphi_j^{n-1}) U^{n-1}_j,
    \qquad\qquad
    i=1\ldots N_n.
  @f}

  如果在这两个时间步骤中使用的网格是相同的，那么 $(\varphi_i^n, \varphi_j^{n-1})$ 就会形成一个方形的质量矩阵 $M_{ij}$  。然而，如果网格不一样，那么一般来说，矩阵是矩形的。更糟的是，甚至很难计算这些积分，因为如果我们在时间步长 $n$ 的网格单元上循环，那么我们需要在这些单元的正交点上评估 $\varphi_j^{n-1}$ ，但它们不一定对应于时间步长 $n-1$ 的网格单元，而且 $\varphi_j^{n-1}$ 不是通过这些单元定义的；如果我们想通过对网格 $n-1$ 的单元积分计算这些积分当然也是一样。

  在任何情况下，我们必须面对的情况是，我们需要整合定义在两个不同网格上的形状函数。这是可以做到的，事实上在步骤28中也有演示，但这个过程最多只能用 "笨拙 "一词来形容。

  在实践中，人们通常不希望这样做。相反，我们在每次调整网格的时候，通过从旧的网格插值到新的网格来避免整个情况。换句话说，我们不是解决上面的方程，而是解决@f{align*}
    \sum_j
    (M
    +
    k_n \theta A)_{ij} U^n_j
    &=
    (\varphi_i, I_h^n u_h^{n-1})


    -
    k_n
    (1-\theta)(\nabla \varphi_i, \nabla I_h^n u_h^{n-1})
    +
    k_n
    \left[
    (1-\theta)F^{n-1}
    +
    \theta F^n
    \right],
  @f}的问题

  其中 $I_h^n$ 是对时间步骤 $n$ 中使用的有限元空间的内插算子。这不是最佳的方法，因为它除了时间和空间离散化之外还引入了一个额外的误差，但这是一个务实的方法，使得做时间适应网格是可行的。   </ul> 




<a name="WhatcouldpossiblygowrongVerifyingwhetherthecodeiscorrect"></a><h3> What could possibly go wrong? Verifying whether the code is correct </h3>


在实现有限元代码时，通常有许多事情会出错。特别是对于时间相关问题，以下是常见的错误来源。

- 时间积分，例如把涉及当前和前一个时间步骤的条款前面的系数弄错了（例如，把一个系数 $\theta$ 混为 $1-\theta$ ）。

- 处理右侧，例如忘记了 $k_n$ 或 $\theta$ 的系数。

- 对边界值处理不当，例如再次忘记了 $k_n$ 或 $\theta$ 的系数，或忘记了不仅对右手边而且对系统矩阵应用非零边界值。

一个不太常见的问题是把初始条件弄错了，因为通常只要输出第一个时间步长就能看出它是错的。在任何情况下，为了验证代码的正确性，有一个测试协议是很有帮助的，它允许我们分别验证这些组件中的每一个。这意味着

- 用非零初始条件但零右手边和边界值测试代码，并验证时间演化是否正确。

- 然后用零初始条件和边界值但非零的右手边进行测试，再次确保正确性。

- 最后，用零初始条件和右手边但非零边界值进行测试。

这听起来很复杂，但幸运的是，对于像这里的无系数（或常数系数）的线性偏微分方程，有一个相当标准的协议，它基于以下观察：如果你选择一个正方形 $[0,1]^2$ 作为你的领域（或者，稍加修改，一个矩形），那么精确解可以写成

@f{align*}
  u(x,y,t) = a(t) \sin(n_x \pi x) \sin(n_y \pi y)


@f}

有整数常数 $n_x,n_y$ ），如果只有初始条件、右手边和边界值也都是 $\sin(n_x \pi x) \sin(n_y \pi y)$ 的形式。这是由于函数 $\sin(n_x \pi x) \sin(n_y \pi y)$ 是拉普拉斯算子的特征函数，允许我们以分析方式计算时间因子 $a(t)$ 等东西，并因此与我们的数值结果进行比较。

作为一个例子，让我们考虑我们有 $u_0(x,y)=\sin(n_x \pi x) \sin(n_x \pi y)$ 和 $f(x,y,t)=0$ 的情况。通过上面对 $u(x,y,t)$ 的形式的主张（ansatz），我们可以得到

@f{align*}
  \left(\frac{\partial}{\partial t} -\Delta\right)
  u(x,y,t)
  &=
  \left(\frac{\partial}{\partial t} -\Delta\right)
  a(t) \sin(n_x \pi x) \sin(n_y \pi y)
  \\
  &=
  \left(a'(t) + (n_x^2+n_y^2)\pi^2 a(t) \right) \sin(n_x \pi x) \sin(n_y \pi y).


@f}

为了使其等于 $f(x,y,t)=0$ ，我们需要： 1.

@f{align*}
  a'(t) + (n_x^2+n_y^2)\pi^2 a(t) = 0


@f}

并由于初始条件， $a(0)=1$  。这个微分方程可以被整合，得到

@f{align*}
  a(t) = - e^{-(n_x^2+n_y^2)\pi^2 t}.


@f}

换句话说，如果初始条件是正弦的乘积，那么解的形状与正弦的乘积完全相同，它以已知的时间依赖性衰减到零。如果你有足够细的网格和足够小的时间步长，这一点是很容易测试的。

如果你把时间积分方案弄错了（例如，在各条款前面有错误的 $\theta$ 或 $k$ 的因子），通常会发生的情况是你没有得到正确的解的时间行为。仔细检查各种因素，直到你得到正确的行为。你可能还想验证一下，时间衰减率（例如，通过绘制固定点的解的值来确定）在你每次将时间步长或网格大小增加一倍或减少一半的时候，都不会增加。你知道这不是对边界条件或右手边的处理，因为这些都是零。

如果你已经如此验证了时间积分器的正确性，那么就拿右手边非零但初始条件为零的情况来说。   $u_0(x,y)=0$  和  $f(x,y,t)=\sin(n_x \pi x) \sin(n_x \pi y)$  。再来看看。

@f{align*}
  \left(\frac{\partial}{\partial t} -\Delta\right)
  u(x,y,t)
  &=
  \left(\frac{\partial}{\partial t} -\Delta\right)
  a(t) \sin(n_x \pi x) \sin(n_y \pi y)
  \\
  &=
  \left(a'(t) + (n_x^2+n_y^2)\pi^2 a(t) \right) \sin(n_x \pi x) \sin(n_y \pi y),


@f}

而要使其等于 $f(x,y,t)$ ，我们需要的是

@f{align*}
  a'(t) + (n_x^2+n_y^2)\pi^2 a(t) = 1


@f}

并由于初始条件， $a(0)=0$  。对这个方程进行时间积分，可以得到

@f{align*}
  a(t) = \frac{1}{(n_x^2+n_y^2)\pi^2} \left[ 1 - e^{-(n_x^2+n_y^2)\pi^2 t} \right].


@f}



同样，如果你在右手边的条款前面有 $\theta$ 或 $k$ 的错误因子，你将不会得到正确的解的时间行为，或者它将收敛到 $\frac{1}{(n_x^2+n_y^2)\pi^2}$ 以外的最大值。

一旦我们验证了使用这个方案的时间积分和右手处理是正确的，我们就可以继续验证我们的边界值是否正确，使用非常类似的方法。




<a name="Thetestcase"></a><h3> The testcase </h3>


在一个简单的域上用一个简单的右手边求解热方程，几乎总是导致解非常无聊，因为它们很快就变得非常光滑，然后就不再有什么变化。相反，我们在这里用零迪里希特边界值和零初始条件在L形域上求解方程，但作为右手边我们选择

@f{align*}
  f(\mathbf x, t)
  =
  \left\{
  \begin{array}{ll}
    \chi_1(\mathbf x)
    & \text{if \(0\le t \le 0.2\tau\) or \(\tau\le t \le 1.2\tau\) or \(2\tau\le t
    \le 2.2\tau\), etc}
    \\
    \chi_2(\mathbf x)
    & \text{if \(0.5\le t \le 0.7\tau\) or \(1.5\tau\le t \le 1.7\tau\) or \(2.5\tau\le t
    \le 2.7\tau\), etc}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.


@f}

在这里。

@f{align*}
  \chi_1(\mathbf x) &=
  \left\{
  \begin{array}{ll}
    1
    & \text{if \(x>0.5\) and \(y>-0.5\)}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.
  \\
  \chi_2(\mathbf x) &=
  \left\{
  \begin{array}{ll}
    1
    & \text{if \(x>-0.5\) and \(y>0.5\)}
    \\
    0
    & \text{otherwise}
  \end{array}
  \right.


@f}

换句话说，在每一个长度为 $\tau$ 的周期中，右手边首先在域1中闪烁，然后完全关闭，然后在域2中打开，然后再次完全关闭。通过<a href="#Results">results
section</a>中显示的解决方案的小动画，这种模式可能是最好的观察。

如果你把热方程解释为寻找导电固体的空间和时间可变的温度分布，那么上面的测试案例相当于一个L形体，我们把边界保持在零温度，并在领域的两个部分交替加热。在加热的同时，这些地方的温度会上升，之后温度会扩散并再次减弱。这些初始条件的意义在于，它们为我们提供了一个在时间上（当源打开和关闭时）以及时间上（在再入角以及源作用区域的边缘和角落）都有奇点的解决方案。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 程序以通常的包含文件开始，所有这些文件你现在应该都见过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/solution_transfer.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 然后照例将这个程序的所有内容放入一个命名空间，并将deal.II命名空间导入到我们将要工作的命名空间中。
 * 

 * 
 * 
 * @code
 * namespace Step26 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeHeatEquationcodeclass"></a> 
 * <h3>The <code>HeatEquation</code> class</h3>
 * 

 * 
 * 下一个部分是这个程序的主类的声明。它沿用了以前的例子中公认的路径。如果你看过 step-6 ，例如，这里唯一值得注意的是，我们需要建立两个矩阵（质量和拉普拉斯矩阵），并保存当前和前一个时间步骤的解。然后，我们还需要存储当前时间、时间步长和当前时间步长的编号。最后一个成员变量表示介绍中讨论的theta参数，它允许我们在一个程序中处理显式和隐式欧拉方法，以及Crank-Nicolson方法和其他通用方法。
 * 

 * 
 * 就成员函数而言，唯一可能的惊喜是 <code>refine_mesh</code> 函数需要最小和最大的网格细化级别的参数。这样做的目的在介绍中已经讨论过了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class HeatEquation 
 *   { 
 *   public: 
 *     HeatEquation(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void solve_time_step(); 
 *     void output_results() const; 
 *     void refine_mesh(const unsigned int min_grid_level, 
 *                      const unsigned int max_grid_level); 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> mass_matrix; 
 *     SparseMatrix<double> laplace_matrix; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> old_solution; 
 *     Vector<double> system_rhs; 
 * 
 *     double       time; 
 *     double       time_step; 
 *     unsigned int timestep_number; 
 * 
 *     const double theta; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 在下面的类和函数中，我们实现了定义这个问题的各种数据（右手边和边界值），这些数据在这个程序中使用，我们需要函数对象。右手边的选择是在介绍的最后讨论的。对于边界值，我们选择零值，但这很容易在下面改变。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     RightHandSide() 
 *       : Function<dim>() 
 *       , period(0.2) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *   private: 
 *     const double period; 
 *   }; 
 * 
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> & p, 
 *                                    const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     AssertIndexRange(component, 1); 
 *     Assert(dim == 2, ExcNotImplemented()); 
 * 
 *     const double time = this->get_time(); 
 *     const double point_within_period = 
 *       (time / period - std::floor(time / period)); 
 * 
 *     if ((point_within_period >= 0.0) && (point_within_period <= 0.2)) 
 *       { 
 *         if ((p[0] > 0.5) && (p[1] > -0.5)) 
 *           return 1; 
 *         else 
 *           return 0; 
 *       } 
 *     else if ((point_within_period >= 0.5) && (point_within_period <= 0.7)) 
 *       { 
 *         if ((p[0] > -0.5) && (p[1] > 0.5)) 
 *           return 1; 
 *         else 
 *           return 0; 
 *       } 
 *     else 
 *       return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
 *                                     const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 *     return 0; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeHeatEquationcodeimplementation"></a> 
 * <h3>The <code>HeatEquation</code> implementation</h3>
 * 

 * 
 * 现在是实现主类的时候了。让我们从构造函数开始，它选择了一个线性元素，一个时间步长为1/500的常数（记得上面把右边的源的一个周期设置为0.2，所以我们用100个时间步长来解决每个周期），并通过设置  $\theta=1/2$  选择了Crank Nicolson方法.
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   HeatEquation<dim>::HeatEquation() 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *     , time_step(1. / 500) 
 *     , theta(0.5) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationsetup_systemcode"></a> 
 * <h4><code>HeatEquation::setup_system</code></h4>
 * 

 * 
 * 下一个函数是设置DoFHandler对象，计算约束，并将线性代数对象设置为正确的大小。我们还在这里通过简单地调用库中的两个函数来计算质量和拉普拉斯矩阵。
 * 

 * 
 * 注意我们在组装矩阵时不考虑悬挂节点的约束（两个函数都有一个AffineConstraints参数，默认为一个空对象）。这是因为我们要在结合当前时间步长的矩阵后，在run()中浓缩约束。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HeatEquation<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::cout << std::endl 
 *               << "===========================================" << std::endl 
 *               << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl 
 *               << std::endl; 
 * 
 *     constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *     constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, 
 *                                     dsp, 
 *                                     constraints, 
 *                                     /*keep_constrained_dofs =  */ true);
 * 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     mass_matrix.reinit(sparsity_pattern); 
 *     laplace_matrix.reinit(sparsity_pattern); 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler, 
 *                                       QGauss<dim>(fe.degree + 1), 
 *                                       mass_matrix); 
 *     MatrixCreator::create_laplace_matrix(dof_handler, 
 *                                          QGauss<dim>(fe.degree + 1), 
 *                                          laplace_matrix); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     old_solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationsolve_time_stepcode"></a> 
 * <h4><code>HeatEquation::solve_time_step</code></h4>
 * 

 * 
 * 下一个函数是解决单个时间步骤的实际线性系统的函数。这里没有什么值得惊讶的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HeatEquation<dim>::solve_time_step() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.0); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *     constraints.distribute(solution); 
 * 
 *     std::cout << "     " << solver_control.last_step() << " CG iterations." 
 *               << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationoutput_resultscode"></a> 
 * <h4><code>HeatEquation::output_results</code></h4>
 * 

 * 
 * 在生成图形输出方面也没有什么新东西，只是我们告诉DataOut对象当前的时间和时间步长是多少，以便将其写入输出文件中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HeatEquation<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "U"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number)); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtk"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationrefine_meshcode"></a> 
 * <h4><code>HeatEquation::refine_mesh</code></h4>
 * 

 * 
 * 这个函数是程序中最有趣的部分。它负责自适应网格细化的工作。这个函数执行的三个任务是：首先找出需要细化/粗化的单元，然后实际进行细化，最后在两个不同的网格之间传输解向量。第一个任务是通过使用成熟的凯利误差估计器来实现的。第二项任务是实际进行再细化。这也只涉及到基本的函数，例如 <code>refine_and_coarsen_fixed_fraction</code> ，它可以细化那些具有最大估计误差的单元，这些误差加起来占60%，并粗化那些具有最小误差的单元，这些单元加起来占40%的误差。请注意，对于像当前这样的问题，即有事发生的区域正在四处移动，我们希望积极地进行粗化，以便我们能够将单元格移动到有必要的地方。
 * 

 * 
 * 正如在介绍中已经讨论过的，太小的网格会导致太小的时间步长，而太大的网格会导致太小的分辨率。因此，在前两个步骤之后，我们有两个循环，将细化和粗化限制在一个允许的单元范围内。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HeatEquation<dim>::refine_mesh(const unsigned int min_grid_level, 
 *                                       const unsigned int max_grid_level) 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(fe.degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       estimated_error_per_cell); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, 
 *                                                       estimated_error_per_cell, 
 *                                                       0.6, 
 *                                                       0.4); 
 * 
 *     if (triangulation.n_levels() > max_grid_level) 
 *       for (const auto &cell : 
 *            triangulation.active_cell_iterators_on_level(max_grid_level)) 
 *         cell->clear_refine_flag(); 
 *     for (const auto &cell : 
 *          triangulation.active_cell_iterators_on_level(min_grid_level)) 
 *       cell->clear_coarsen_flag(); 
 * 
 * @endcode
 * 
 * 上面这两个循环略有不同，但这很容易解释。在第一个循环中，我们没有调用  <code>triangulation.end()</code>  ，而是调用  <code>triangulation.end_active(max_grid_level)</code>  。这两个调用应该产生相同的迭代器，因为迭代器是按级别排序的，不应该有任何级别高于 <code>max_grid_level</code> 的单元格。事实上，这段代码确保了这种情况的发生。
 * 

 * 
 * 作为网格细化的一部分，我们需要将旧的网格中的解向量转移到新的网格中。为此，我们使用了SolutionTransfer类，我们必须准备好需要转移到新网格的解向量（一旦完成细化，我们将失去旧的网格，所以转移必须与细化同时发生）。在我们调用这个函数的时候，我们将刚刚计算出解决方案，所以我们不再需要old_solution变量（它将在网格被细化后被解决方案覆盖，也就是在时间步长结束时；见下文）。换句话说，我们只需要一个求解向量，并将其复制到一个临时对象中，当我们进一步向下调用 <code>setup_system()</code> 时，它就不会被重置。
 * 

 * 
 * 因此，我们将一个SolutionTransfer对象附加到旧的DoF处理程序中，以初始化它。然后，我们准备好三角形和数据向量，以便进行细化（按照这个顺序）。
 * 

 * 
 * 
 * @code
 *     SolutionTransfer<dim> solution_trans(dof_handler); 
 * 
 *     Vector<double> previous_solution; 
 *     previous_solution = solution; 
 *     triangulation.prepare_coarsening_and_refinement(); 
 *     solution_trans.prepare_for_coarsening_and_refinement(previous_solution); 
 * 
 * @endcode
 * 
 * 现在一切都准备好了，所以进行细化并在新网格上重新创建DoF结构，最后在 <code>setup_system</code> 函数中初始化矩阵结构和新的向量。接下来，我们实际执行从旧网格到新网格的插值解。最后一步是对解向量应用悬空节点约束，即确保位于悬空节点上的自由度值，使解是连续的。这是必要的，因为SolutionTransfer只对单元格进行局部操作，不考虑邻域。
 * 

 * 
 * 
 * @code
 *     triangulation.execute_coarsening_and_refinement(); 
 *     setup_system(); 
 * 
 *     solution_trans.interpolate(previous_solution, solution); 
 *     constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeHeatEquationruncode"></a> 
 * <h4><code>HeatEquation::run</code></h4>
 * 

 * 
 * 这是程序的主要驱动，我们在这里循环所有的时间步骤。在函数的顶部，我们通过重复第一个时间步长，设置初始全局网格细化的数量和自适应网格细化的初始周期数量。然后，我们创建一个网格，初始化我们要处理的各种对象，设置一个标签，说明我们在重新运行第一个时间步长时应该从哪里开始，并将初始解插值到网格上（我们在这里选择了零函数，当然，我们可以用更简单的方法，直接将解向量设置为零）。我们还输出一次初始时间步长。
 * 

 * 
 * @note  如果你是一个有经验的程序员，你可能会对我们在这段代码中使用 <code>goto</code> 语句感到吃惊   <code>goto</code> 语句现在已经不是特别受人欢迎了，因为计算机科学界的大师之一Edsgar Dijkstra在1968年写了一封信，叫做 "Go To Statement considered harmful"（见<a href="http:en.wikipedia.org/wiki/Considered_harmful">here</a>）。这段代码的作者全心全意地赞同这一观念。  <code>goto</code> 是难以理解的。事实上，deal.II几乎不包含任何出现的情况：不包括基本上是从书本上转录的代码，也不计算重复的代码片断，在写这篇笔记时，大约60万行代码中有3个位置；我们还在4个教程程序中使用它，其背景与这里完全相同。与其在这里试图证明这种情况的出现，不如先看看代码，我们在函数的最后再来讨论这个问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void HeatEquation<dim>::run() 
 *   { 
 *     const unsigned int initial_global_refinement       = 2; 
 *     const unsigned int n_adaptive_pre_refinement_steps = 4; 
 * 
 *     GridGenerator::hyper_L(triangulation); 
 *     triangulation.refine_global(initial_global_refinement); 
 * 
 *     setup_system(); 
 * 
 *     unsigned int pre_refinement_step = 0; 
 * 
 *     Vector<double> tmp; 
 *     Vector<double> forcing_terms; 
 * 
 *   start_time_iteration: 
 * 
 *     time            = 0.0; 
 *     timestep_number = 0; 
 * 
 *     tmp.reinit(solution.size()); 
 *     forcing_terms.reinit(solution.size()); 
 * 
 *     VectorTools::interpolate(dof_handler, 
 *                              Functions::ZeroFunction<dim>(), 
 *                              old_solution); 
 *     solution = old_solution; 
 * 
 *     output_results(); 
 * 
 * @endcode
 * 
 * 然后我们开始主循环，直到计算的时间超过我们的结束时间0.5。第一个任务是建立我们需要在每个时间步骤中解决的线性系统的右手边。回顾一下，它包含项 $MU^{n-1}-(1-\theta)k_n AU^{n-1}$  。我们把这些项放到变量system_rhs中，借助于一个临时矢量。
 * 

 * 
 * 
 * @code
 *     while (time <= 0.5) 
 *       { 
 *         time += time_step; 
 *         ++timestep_number; 
 * 
 *         std::cout << "Time step " << timestep_number << " at t=" << time 
 *                   << std::endl; 
 * 
 *         mass_matrix.vmult(system_rhs, old_solution); 
 * 
 *         laplace_matrix.vmult(tmp, old_solution); 
 *         system_rhs.add(-(1 - theta) * time_step, tmp); 
 * 
 * @endcode
 * 
 * 第二块是计算源项的贡献。这与术语  $k_n \left[ (1-\theta)F^{n-1} + \theta F^n \right]$  相对应。下面的代码调用  VectorTools::create_right_hand_side  来计算向量  $F$  ，在这里我们在评估之前设置了右侧（源）函数的时间。这一切的结果最终都在forcing_terms变量中。
 * 

 * 
 * 
 * @code
 *         RightHandSide<dim> rhs_function; 
 *         rhs_function.set_time(time); 
 *         VectorTools::create_right_hand_side(dof_handler, 
 *                                             QGauss<dim>(fe.degree + 1), 
 *                                             rhs_function, 
 *                                             tmp); 
 *         forcing_terms = tmp; 
 *         forcing_terms *= time_step * theta; 
 * 
 *         rhs_function.set_time(time - time_step); 
 *         VectorTools::create_right_hand_side(dof_handler, 
 *                                             QGauss<dim>(fe.degree + 1), 
 *                                             rhs_function, 
 *                                             tmp); 
 * 
 *         forcing_terms.add(time_step * (1 - theta), tmp); 
 * 
 * @endcode
 * 
 * 接下来，我们将强迫项加入到来自时间步长的强迫项中，同时建立矩阵 $M+k_n\theta A$ ，我们必须在每个时间步长中进行反转。这些操作的最后一块是消除线性系统中悬挂的节点约束自由度。
 * 

 * 
 * 
 * @code
 *         system_rhs += forcing_terms; 
 * 
 *         system_matrix.copy_from(mass_matrix); 
 *         system_matrix.add(theta * time_step, laplace_matrix); 
 * 
 *         constraints.condense(system_matrix, system_rhs); 
 * 
 * @endcode
 * 
 * 在解决这个问题之前，我们还需要做一个操作：边界值。为此，我们创建一个边界值对象，将适当的时间设置为当前时间步长的时间，并像以前多次那样对其进行评估。其结果也被用来在线性系统中设置正确的边界值。
 * 

 * 
 * 
 * @code
 *         { 
 *           BoundaryValues<dim> boundary_values_function; 
 *           boundary_values_function.set_time(time); 
 * 
 *           std::map<types::global_dof_index, double> boundary_values; 
 *           VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                    0, 
 *                                                    boundary_values_function, 
 *                                                    boundary_values); 
 * 
 *           MatrixTools::apply_boundary_values(boundary_values, 
 *                                              system_matrix, 
 *                                              solution, 
 *                                              system_rhs); 
 *         } 
 * 
 * @endcode
 * 
 * 有了这些，我们要做的就是解决这个系统，生成图形数据，以及......
 * 

 * 
 * 
 * @code
 *         solve_time_step(); 
 * 
 *         output_results(); 
 * 
 * @endcode
 * 
 * ...负责网格的细化。在这里，我们要做的是：(i)在求解过程的最开始，细化所要求的次数，之后我们跳到顶部重新开始时间迭代，(ii)之后每隔五步细化一次。
 * 

 * 
 * 时间循环和程序的主要部分以开始进入下一个时间步骤结束，将old_solution设置为我们刚刚计算出的解决方案。
 * 

 * 
 * 
 * @code
 *         if ((timestep_number == 1) && 
 *             (pre_refinement_step < n_adaptive_pre_refinement_steps)) 
 *           { 
 *             refine_mesh(initial_global_refinement, 
 *                         initial_global_refinement + 
 *                           n_adaptive_pre_refinement_steps); 
 *             ++pre_refinement_step; 
 * 
 *             tmp.reinit(solution.size()); 
 *             forcing_terms.reinit(solution.size()); 
 * 
 *             std::cout << std::endl; 
 * 
 *             goto start_time_iteration; 
 *           } 
 *         else if ((timestep_number > 0) && (timestep_number % 5 == 0)) 
 *           { 
 *             refine_mesh(initial_global_refinement, 
 *                         initial_global_refinement + 
 *                           n_adaptive_pre_refinement_steps); 
 *             tmp.reinit(solution.size()); 
 *             forcing_terms.reinit(solution.size()); 
 *           } 
 * 
 *         old_solution = solution; 
 *       } 
 *   } 
 * } // namespace Step26 
 * 
 * @endcode
 * 
 * 现在你已经看到了这个函数的作用，让我们再来看看  <code>goto</code>  的问题。从本质上讲，代码所做的事情是这样的。
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *    void run ()
 *    {
 *      initialize;
 *    start_time_iteration:
 *      for (timestep=1...)
 *      {
 *         solve timestep;
 *         if (timestep==1 && not happy with the result)
 *         {
 *           adjust some data structures;
 *           goto start_time_iteration; simply try again
 *         }
 *         postprocess;
 *      }
 *    }
 *  @endcode 
 * </div>
 * 这里，"对结果满意 "的条件是我们想保留当前的网格，还是宁愿细化网格并在新网格上重新开始。我们当然可以用下面的方法来取代  <code>goto</code>  的使用。
 * <div class=CodeFragmentInTutorialComment>
 *  @code
 *    void run ()
 *    {
 *      initialize;
 *      while (true)
 *      {
 *         solve timestep;
 *         if (not happy with the result)
 *            adjust some data structures;
 *         else
 *            break;
 *      }
 *      postprocess;
 *

 *

 *      for (timestep=2...)
 *      {
 *         solve timestep;
 *         postprocess;
 *      }
 *    }
 *  @endcode 
 * </div>
 * 这样做的好处是摆脱了 <code>goto</code> ，但缺点是必须在两个不同的地方重复实现 "解算时间步长 "和 "后处理 "操作的代码。这可以通过将这些部分的代码（在上面的实际实现中是相当大的块）放到自己的函数中来解决，但是一个带有 <code>break</code> 语句的 <code>while(true)</code> 循环并不真的比 <code>goto</code> 容易阅读或理解。
 * 

 * 
 * 最后，人们可能会简单地同意，<i>in general</i> 。
 * <code>goto</code>  语句是个坏主意，但要务实地指出，在某些情况下，它们可以帮助避免代码重复和尴尬的控制流。这可能就是其中之一，它与Steve McConnell在他关于良好编程实践的优秀书籍 "Code Complete"  @cite CodeComplete 中采取的立场一致（见 step-1 的介绍中提到的这本书），该书花了惊人的10页来讨论一般的 <code>goto</code> 问题。
 * 

 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 走到这一步，这个程序的主函数又没有什么好讨论的了：它看起来就像自 step-6 以来的所有此类函数一样。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step26; 
 * 
 *       HeatEquation<2> heat_equation_solver; 
 *       heat_equation_solver.run(); 
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
examples/step-26/doc/results.dox



<a name="Results"></a><h1>Results</h1>


正如许多教程一样，程序的实际输出并不重要，重要的是我们是如何到达那里的。尽管如此，它还是来了。

@code
===========================================
Number of active cells: 48
Number of degrees of freedom: 65


Time step 1 at t=0.002
     7 CG iterations.


===========================================
Number of active cells: 60
Number of degrees of freedom: 81



Time step 1 at t=0.002
     7 CG iterations.


===========================================
Number of active cells: 105
Number of degrees of freedom: 136



Time step 1 at t=0.002
     7 CG iterations.


[...]


Time step 249 at t=0.498
     13 CG iterations.
Time step 250 at t=0.5
     14 CG iterations.


===========================================
Number of active cells: 1803
Number of degrees of freedom: 2109
@endcode



也许更有意义的是解决方案的可视化和计算的网格。

 <img src="https://www.dealii.org/images/steps/developer/step-26.movie.gif" alt="Animation of the solution of step 26."> 

这部电影显示了两个信号源是如何开启和关闭的，以及网格对此的反应。很明显，现在的网格可能不是我们能想出的最好的。我们将在下一节再来讨论这个问题。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


至少有两个方面可以大大改善这个程序：自适应时间步长和更好地选择网格。

<a name="Adaptivetimestepping"></a><h4>Adaptive time stepping</h4>


由于选择了隐式时间步进方案，我们不受任何类似CFL的时间步进条件的约束。此外，由于在热力方程中发生变化的时间尺度不受细胞直径的约束（与波浪方程的情况不同，在波浪方程中我们有一个固定的信息传输速度，将时间尺度和空间尺度联系起来），我们可以随意选择时间步长。或者，最好是按照我们认为必要的准确性来选择它。

看一下解决方案，很明显，行动并不是随时间均匀发生的：在我们打开一个源的时候，很多东西都在发生变化，一旦一个源开启了一段时间，事情就变得不那么戏剧化了，而当两个源都关闭的时候，我们进入了一个漫长的下降阶段。在这些时候，我们肯定可以在不牺牲太多精确度的情况下，采用比以前更大的时间步长。

文献中有许多关于如何自适应地选择时间步长的建议。例如，可以从ODE求解器选择其时间步长的方式中学到很多。我们也可以从后验误差估计中得到启发，理想情况下，后验误差估计可以写成由对整体误差的时间和空间贡献组成。如果时间上的贡献太大，我们应该选择一个较小的时间步长。例如，这个方向的想法可以在deal.II的前主要开发者Ralf Hartmann的博士论文中找到，该论文由德国海德堡大学在2002年出版。




<a name="Bettertimesteppingmethods"></a><h4>Better time stepping methods</h4>


我们在这里使用了一种比较简单的时间步进方法，即二阶时间的Crank-Nicolson方法。然而，更精确的方法如Runge-Kutta方法是可用的，并且应该使用，因为它们并不代表很多额外的努力。对于目前的程序来说，实现这一点并不困难，但在step-52中也给出了一个更系统的处理。




<a name="Betterrefinementcriteria"></a><h4>Better refinement criteria</h4>


如果你看一下上面电影中的网格，很明显，它们不是特别适合手头的任务。事实上，它们看起来相当随意。

有两个因素在起作用。首先，有一些岛屿，其中的单元已经被细化，但周围是未细化的单元（可能还有一些偶尔被粗化的岛屿）。这些并不可怕，因为它们大多数时候并不影响网格的近似质量，但是它们也没有帮助，因为它们的许多额外的自由度事实上受到悬挂节点约束的限制。也就是说，这很容易解决：Triangulation类在其构造函数中接受一个参数，表示 "网格平滑 "的程度。传递许多可能的标志之一，这将指示三角剖分细化一些额外的单元，或者不细化一些单元，这样得到的网格就不会有这些伪影。

第二个问题更为严重：网格似乎滞后于解。其根本原因是我们每隔五步才调整一次网格，而且在这些情况下只允许进行一次细化。每当一个源打开时，之前的解在这个区域是非常平滑的，因此网格也是相当粗糙的。这意味着在下一个时间步骤中，当我们细化网格时，我们会在这个区域多得到一个细化级别，五个时间步骤后再多一个级别，等等。但这还不够：首先，我们应该在一个源打开时立即进行细化（毕竟在当前情况下，我们至少知道右手边是什么），而且我们应该允许超过一个细化级别。当然，所有这些都可以用deal.II来完成，只是需要在如何使其发挥作用方面有一些算法上的思考!




<a name="Positivitypreservation"></a><h4>Positivity preservation</h4>


为了提高你的模拟在时间上的准确性和分辨率，通常会减少时间步长  $k_n$  。如果你在这个特定的例子中开始玩弄时间步长，你会注意到，如果 $k_n$ 低于某个阈值，解决方案会变成部分负数。这不是我们所期望发生的（在自然界）。

为了从数学上了解这种行为，让我们考虑一个一般的、完全离散的问题。

@f{align*}
  A u^{n} = B u^{n-1}.


@f}

那么 $i$ th方程的一般形式为：。

@f{align*}
  a_{ii} u^{n}_i &= b_{ii} u^{n-1}_i +
  \sum\limits_{j \in S_i} \left( b_{ij} u^{n-1}_j - a_{ij} u^{n}_j \right),


@f}

其中 $S_i$ 是与DoF $i$ 耦合的自由度集合（即矩阵 $A$ 或矩阵 $B$ 在位置 $(i,j)$ 有一个非零条目）。如果所有系数都满足以下条件。

@f{align*}
  a_{ii} &> 0, & b_{ii} &\geq 0, & a_{ij} &\leq 0, & b_{ij} &\geq 0,
  &
  \forall j &\in S_i,


@f}

所有的解决方案 $u^{n}$ 都保持其与之前的解决方案 $u^{n-1}$ 的符号，因此也保持其与初始值 $u^0$ 的符号。关于正性保留的更多信息，请参见例如 <a href="http://bookstore.siam.org/cs14/">Kuzmin, H&auml;m&auml;l&auml;inen</a> 。

根据要解决的PDE和使用的时间积分方案，人们能够推导出时间步长的条件  $k_n$  。对于采用Crank-Nicolson方案的热方程，<a href="https://doi.org/10.2478/cmam-2010-0025">Schatz et. al.</a>已将其转化为下列方案。

@f{align*}
  (1 - \theta) k a_{ii} &\leq m_{ii},\qquad \forall i,
  &
  \theta k \left| a_{ij} \right| &\geq m_{ij},\qquad j \neq i,


@f}

其中 $M = m_{ij}$ 表示质量矩阵， $A = a_{ij}$ 表示刚度矩阵， $a_{ij} \leq 0$ 分别表示 $j \neq i$ 。有了 $a_{ij} \leq 0$ ，我们可以制定全局时间步长 $k$ 的界限如下。

@f{align*}
  k_{\text{max}} &= \frac{ 1 }{ 1 - \theta }
  \min\left( \frac{ m_{ii} }{ a_{ii} } \right),~ \forall i,
  &
  k_{\text{min}} &= \frac{ 1 }{ \theta  }
  \max\left( \frac{ m_{ij} }{ \left|a_{ij}\right| } \right),~ j \neq i.


@f}

换句话说，在Crank-Nicolson方案的情况下，时间步长受到<i>both a lower
and upper bound</i>的限制。这些约束应与CFL条件一起考虑，以确保所进行的模拟的重要性。

无法使时间步长达到我们想要的小，以获得更多的精度，而又不失去积极性属性，这是令人讨厌的。这就提出了一个问题：在这个特定的教程中，我们是否至少可以<i>compute</i>选择最小的时间步长来确保正性的保留。事实上，我们可以使用通过MatrixCreator函数创建的质量和刚度的SparseMatrix对象。通过SparseMatrixIterators遍历每个条目，我们可以检查对角线和非对角线条目，从而动态地设置一个合适的时间步长。对于二次元矩阵，对角线元素被存储为一行的第一个成员（见SparseMatrix文档）。下面是一个关于如何从 <code>mass_matrix</code> 中抓取感兴趣的条目的示范性代码片断。

@code
Assert (mass_matrix.m() == mass_matrix.n(), ExcNotQuadratic());
const unsigned int num_rows = mass_matrix.m();
double mass_matrix_min_diag    = std::numeric_limits<double>::max(),
       mass_matrix_max_offdiag = 0.;


SparseMatrixIterators::Iterator<double,true> row_it (&mass_matrix, 0);


for(unsigned int m = 0; m<num_rows; ++m)
{
  // check the diagonal element
  row_it = mass_matrix.begin(m);
  mass_matrix_min_diag = std::min(row_it->value(), mass_matrix_min_diag);
  ++row_it;


  // check the off-diagonal elements
  for(; row_it != mass_matrix.end(m); ++row_it)
    mass_matrix_max_offdiag = std::max(row_it->value(), mass_matrix_max_offdiag);
}
@endcode



利用这样计算出来的信息，我们可以通过上述公式来约束时间步长。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-26.cc"
*/
