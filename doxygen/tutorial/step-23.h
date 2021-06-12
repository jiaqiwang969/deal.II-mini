/**
@page step_23 The step-23 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Timediscretization">Time discretization</a>
      <ul>
        <li><a href="#MethodoflinesorRothesmethod">Method of lines or Rothe's method?</a>
        <li><a href="#Rothesmethod">Rothe's method!</a>
      </ul>
        <li><a href="#Spacediscretization">Space discretization</a>
        <li><a href="#Energyconservation">Energy conservation</a>
        <li><a href="#WhoareCourantFriedrichsandLewy">Who are Courant, Friedrichs, and Lewy?</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeWaveEquationcodeclass">The <code>WaveEquation</code> class</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeWaveEquationcodeclass">Implementation of the <code>WaveEquation</code> class</a>
      <ul>
        <li><a href="#WaveEquationsetup_system">WaveEquation::setup_system</a>
        <li><a href="#WaveEquationsolve_uandWaveEquationsolve_v">WaveEquation::solve_u and WaveEquation::solve_v</a>
        <li><a href="#WaveEquationoutput_results">WaveEquation::output_results</a>
        <li><a href="#WaveEquationrun">WaveEquation::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-23/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{28} 

这是一系列教程程序中的第一个，它将最终涵盖 "真正的 "时间依赖问题，而不是步骤18中发现的略显奇怪的时间依赖形式或步骤21的DAE模型。特别是，这个程序介绍了有界域中的波浪方程。后来，第24步将考虑一个吸收边界条件的例子，以及 @ref
step_25 "第25步 "一种产生称为孤子的解的非线性波方程。

波浪方程的原型形式如下：找到满足以下条件的 $u(x,t), x\in\Omega, t\in[0,T]$ 。

@f{eqnarray*}
	\frac{\partial^2 u}{\partial t^2}


	-
	\Delta u &=& f
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	u(x,t) &=& g
	\qquad
	\textrm{on}\ \partial\Omega\times [0,T],
\\
	u(x,0) &=& u_0(x)
	\qquad
	\textrm{in}\ \Omega,
\\
	\frac{\partial u(x,0)}{\partial t} &=& u_1(x)
	\qquad
	\textrm{in}\ \Omega.


@f}

请注意，由于这是一个具有二阶时间导数的方程，我们需要提出两个初始条件，一个是值，一个是解的时间导数。

在物理上，该方程描述了弹性介质的运动。在二维空间中，我们可以考虑膜在受到力的作用下如何运动。上面的Dirichlet边界条件表明，膜被夹在边界的高度 $g(x,t)$ （这个高度可能也在移动&mdash；想想人们拿着毯子上下摇晃）。第一个初始条件等于膜的初始偏转，而第二个初始条件给出了其速度。例如，我们可以考虑用手指把膜推下去，然后在 $t=0$ 处让它离开（非零偏转但零初始速度），或者在 $t=0$ 处用锤子砸它（零偏转但非零速度）。这两种情况都会引起膜的运动。




<a name="Timediscretization"></a><h3>Time discretization</h3>


<a name="MethodoflinesorRothesmethod"></a><h4>Method of lines or Rothe's method?</h4>在数值分析界有一个长期的争论，即时间依赖方程的离散化是否应该首先离散时间变量，导致每个时间步长的静止PDE，然后用标准的有限元技术来解决（这被称为Rothe方法），或者是否应该首先离散空间变量，导致一个大型的常微分方程系统，然后用一个通常的ODE求解器来处理（这被称为线的方法）。


这两种方法都有优点和缺点。传统上，人们更倾向于线的方法，因为它允许使用非常发达的高阶ODE求解器，可用于由这种方法产生的相当刚性的ODE，包括步长控制和时间误差的估计。

另一方面，当使用高阶时间步长法时，罗特的方法变得很尴尬，因为这时人们必须写下一个PDE，将当前时间步长的解与前一个时间步长的解结合起来，而且可能还有更早的解，从而导致大量的条款。

由于这些原因，线条法在很长一段时间内都是人们的首选方法。然而，它有一个很大的缺点：如果我们先将空间变量离散化，导致一个大的ODE系统，我们必须一劳永逸地选择一个网格。如果我们愿意这样做，那么这就是一种合法的、可能是优越的方法。

另一方面，如果我们看的是波浪方程和其他许多与时间有关的问题，我们会发现，随着时间的推移，解的特征会发生变化。例如，对于波浪方程，我们可能有一个单一的波浪穿过域，在波浪的前后，解是平滑的，甚至是恒定的&mdash;自适应性对于这种情况确实很有用，但关键是我们需要细化网格的区域会随着时间步数的变化而变化!

如果我们打算这样做，即为每个时间步长（或一组时间步长）选择不同的网格，那么线段法就不再合适了：我们不是得到一个变量数等于有限元网格中未知数数量的ODE系统，而是未知数的数量一直在变化，这是标准ODE求解器肯定不准备处理的事实。另一方面，对于罗特方法，我们只是在每个时间步长得到一个PDE，我们可以选择独立于前一个时间步长所用的网格进行离散化；这种方法并非没有危险和困难，但至少是一个合理的、定义明确的程序。

由于所有这些原因，在本程序中，我们选择使用Rothe方法进行离散化，即我们首先在时间上进行离散化，然后在空间上离散化。我们实际上根本不会使用自适应网格，因为这涉及到大量的额外代码，但我们将在<a href="#Results">results section below</a>中对此作一些评论。




<a name="Rothesmethod"></a><h4>Rothe's method!</h4>


鉴于这些考虑，我们将这样做：让我们首先为这个二阶问题定义一个简单的时间步进方法，然后在第二步做空间离散化，即我们将遵循Rothe的方法。

对于第一步，让我们先绕一点路：为了离散化一个第二时间导数，我们可以直接离散化它，或者引入一个额外的变量，将系统转化为一阶系统。在许多情况下，这证明是等价的，但处理一阶系统往往更简单。为此，让我们引入

@f[
	v = \frac{\partial u}{\partial t},


@f]

并称这个变量为<i>velocity</i>，原因很明显。然后，我们可以将原来的波浪方程重新表述如下。

@f{eqnarray*}
	\frac{\partial u}{\partial t}


	-
	v
	&=& 0
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	\frac{\partial v}{\partial t}


	-
	\Delta u &=& f
	\qquad
	\textrm{in}\ \Omega\times [0,T],
\\
	u(x,t) &=& g
	\qquad
	\textrm{on}\ \partial\Omega\times [0,T],
\\
	u(x,0) &=& u_0(x)
	\qquad
	\textrm{in}\ \Omega,
\\
	v(x,0) &=& u_1(x)
	\qquad
	\textrm{in}\ \Omega.


@f}

这个公式的优点是它现在只包含两个变量的第一时间导数，对于它来说，写下时间步进方案很简单。请注意，我们起初没有 $v$ 的边界条件。然而，我们可以在边界上强制执行 $v=\frac{\partial
g}{\partial t}$ 。在数值例子中发现，这实际上是必要的：如果不这样做，解决方案看起来并不特别错误，但如果不执行这些边界条件，Crank-Nicolson方案并不保存能量。

有了这个公式，让我们引入以下时间离散化，其中上标 $n$ 表示一个时间步长， $k=t_n-t_{n-1}$ 是当前时间步长。

\f{eqnarray*}
  \frac{u^n - u^{n-1}}{k}


  - \left[\theta v^n + (1-\theta) v^{n-1}\right] &=& 0,
  \\
  \frac{v^n - v^{n-1}}{k}


  - \Delta\left[\theta u^n + (1-\theta) u^{n-1}\right]
  &=& \theta f^n + (1-\theta) f^{n-1}.
\f} 注意我们在这里引入了一个参数 $\theta$ 。例如，如果我们选择 $\theta=0$ ，第一个方程将简化为 $\frac{u^n - u^{n-1}}{k}  - v^{n-1} = 0$ ，这就是众所周知的正向或显式欧拉方法。另一方面，如果我们设定 $\theta=1$ ，那么我们将得到 $\frac{u^n - u^{n-1}}{k}  - v^n = 0$ ，这对应于后向或隐式欧拉方法。这两种方法都是一阶精确方法。它们实现起来很简单，但其实并不十分精确。

第三种情况是选择 $\theta=\frac 12$  。然后上面的第一个方程将变成 $\frac{u^n - u^{n-1}}{k}


- \frac 12 \left[v^n + v^{n-1}\right] = 0$  。这种方法被称为Crank-Nicolson方法，它的优点是二阶精确。此外，它还有一个很好的特性，即保留了溶液中的能量（从物理上讲，能量是膜中粒子的动能加上由于局部拉伸而存在的势能的总和；这个量在连续方程中是一个守恒量，但大多数时间步进方案在时间离散化后并不保留它）。由于 $v^n$ 也出现在 $u^n$ 的方程中，Crank-Nicolson方案也是隐式的。

在程序中，我们将把 $\theta$ 作为一个参数，这样就很容易发挥它的作用了。结果部分将显示一些比较不同方案的数字证据。

上面的方程（称为<i>semidiscretized</i>方程，因为我们只离散了时间，而没有离散空间），可以通过从第一个方程中消除 $v^n$ 和重新排列项来简化一下。然后我们得到

\f{eqnarray*}
  \left[ 1-k^2\theta^2\Delta \right] u^n &=&
  	 \left[ 1+k^2\theta(1-\theta)\Delta\right] u^{n-1} + k v^{n-1}
   	 + k^2\theta\left[\theta f^n + (1-\theta) f^{n-1}\right],\\
   v^n &=& v^{n-1} + k\Delta\left[ \theta u^n + (1-\theta) u^{n-1}\right]
   + k\left[\theta f^n + (1-\theta) f^{n-1}\right].
\f}在这种形式下，我们看到，如果我们得到了前一个时间段的解 $u^{n-1},v^{n-1}$ ，那么我们就可以分别求解变量 $u^n,v^n$ ，也就是一次一个。这很方便。此外，我们认识到第一个方程中的算子是正定的，而第二个方程看起来特别简单。




<a name="Spacediscretization"></a><h3>Space discretization</h3>


我们现在已经得出了将时间 $u^n(x)$ 的近似（半离散）解 $v^n(x)$ 及其时间导数 $t_n$ 与前一个时间步骤 $t_{n-1}$ 的解 $u^{n-1}(x),v^{n-1}(x)$ 相关的方程。下一步是使用通常的有限元方法将空间变量离散化。为此，我们将每个方程与一个测试函数相乘，在整个域上进行积分，并在必要时进行部分积分。这就导致了

\f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla u^{n-1},\nabla \varphi)
  +
  k(v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right],
  \\
  (v^n,\varphi)
   &=&
   (v^{n-1},\varphi)


    -
    k\left[ \theta (\nabla u^n,\nabla\varphi) +
    (1-\theta) (\nabla u^{n-1},\nabla \varphi)\right]
  + k
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right].
\f}

然后习惯于近似 $u^n(x) \approx u^n_h(x) = \sum_i
U_i^n\phi_i^n(x)$  ，其中 $\phi_i^n(x)$  是用于离散化 $n$  -个时间步长的形状函数， $U_i^n$  是解决方案的未知结点值。同样地， $v^n(x) \approx
v^n_h(x) = \sum_i V_i^n\phi_i^n(x)$  。最后，我们有前一个时间步骤的解，  $u^{n-1}(x) \approx u^{n-1}_h(x) = \sum_i
U_i^{n-1}\phi_i^{n-1}(x)$  和  $v^{n-1}(x) \approx v^{n-1}_h(x) = \sum_i
V_i^{n-1}\phi_i^{n-1}(x)$  。请注意，由于在我们到达时间步骤 $n$ 时，前一个时间步骤的解已经被计算出来了，所以 $U^{n-1},V^{n-1}$ 是已知的。此外，注意上一步的解可能是在不同的网格上计算的，所以我们必须使用形状函数  $\phi^{n-1}_i(x)$  。

如果我们将这些扩展插入上述方程，并用本网的测试函数进行测试，我们得到以下线性系统。

\f{eqnarray*}
  (M^n + k^2\theta^2 A^n)U^n &=&
  M^{n,n-1}U^{n-1} - k^2\theta(1-\theta) A^{n,n-1}U^{n-1}
  +
  kM^{n,n-1}V^{n-1}
  + k^2\theta
  \left[
  \theta F^n + (1-\theta) F^{n-1}
  \right],
  \\
  M^nV^n
   &=&
   M^{n,n-1}V^{n-1}


    -
    k\left[ \theta A^n U^n +
    (1-\theta) A^{n,n-1} U^{n-1}\right]
   + k
  \left[
  \theta F^n + (1-\theta) F^{n-1}
  \right],
\f} 其中

@f{eqnarray*}
	M^n_{ij} &=& (\phi_i^n, \phi_j^n),
	\\
	A^n_{ij} &=& (\nabla\phi_i^n, \nabla\phi_j^n),
	\\
	M^{n,n-1}_{ij} &=& (\phi_i^n, \phi_j^{n-1}),
	\\
	A^{n,n-1}_{ij} &=& (\nabla\phi_i^n, \nabla\phi_j^{n-1}),
	\\
	F^n_{i} &=& (f^n,\phi_i^n),
	\\
	F^{n-1}_{i} &=& (f^{n-1},\phi_i^n).


@f}



如果我们解决这两个方程，我们可以将解决方案向前推进一步，并进入下一个时间步骤。

值得注意的是，如果我们在每个时间步长选择相同的网格（事实上我们将在下面的程序中这样做），那么我们在时间步长 $n$ 和 $n-1$ 上有相同的形状函数，即 $\phi^n_i=\phi_i^{n-1}=\phi_i$  。因此，我们得到  $M^n=M^{n,n-1}=M$  和  $A^n=A^{n,n-1}=A$  。另一方面，如果我们使用了不同的形状函数，那么我们将不得不计算包含定义在两个网格上的形状函数的积分。这是一个有些混乱的过程，我们在此省略，但在步骤28中会有一些详细的处理。

在这些条件下（即网格不发生变化），我们可以通过基本消除第二个线性系统的解来优化求解过程。我们将在 @ref step_25 "step-25 "程序的介绍中讨论这个问题。

<a name="Energyconservation"></a><h3>Energy conservation</h3>


比较时间步进方案质量的一个方法是看数值近似是否保留了连续方程的守恒特性。对于波浪方程来说，自然要看的是能量。通过将波浪方程乘以 $u_t$ ，对 $\Omega$ 进行积分，并在必要时进行部分积分，我们发现

@f[
	\frac{d}{d t}
	\left[\frac 12 \int_\Omega \left(\frac{\partial u}{\partial
	t}\right)^2 + (\nabla u)^2 \; dx\right]
	=
	\int_\Omega f \frac{\partial u}{\partial t} \; dx
	+
	\int_{\partial\Omega} n\cdot\nabla u
	\frac{\partial g}{\partial t} \; dx.


@f]

因此，在没有体力和恒定边界值的情况下，我们得到的结果是

@f[
	E(t) = \frac 12 \int_\Omega \left(\frac{\partial u}{\partial
	t}\right)^2 + (\nabla u)^2 \; dx


@f]

是一个守恒量，即一个不随时间变化的量。我们将在每个时间步骤后计算这个量。很容易看出，如果我们用有限元近似值代替 $u$ ，用速度的有限元近似值代替 $\frac{\partial u}{\partial t}$ ，那么

@f[
	E(t_n) = \frac 12 \left<V^n, M^n V^n\right>
	+
	\frac 12 \left<U^n, A^n U^n\right>.


@f]

正如我们将在结果部分看到的，Crank-Nicolson方案确实保存了能量，而前向和后向Euler方案都没有。




<a name="WhoareCourantFriedrichsandLewy"></a><h3>Who are Courant, Friedrichs, and Lewy?</h3>


波浪方程的数值求解很麻烦，原因之一是显式时间离散化只有在时间步长足够小的情况下才稳定。特别是，它与空间网格宽度有耦合关系  $h$  。对于我们这里使用的最低阶离散化，其关系为

@f[
	k\le \frac hc


@f]

其中 $c$ 是波速，在我们对波浪方程的表述中，它已被归一。因此，除非我们使用带有 $\theta>0$ 的隐式方案，否则如果我们违反这一限制，我们的解在数值上是不稳定的。隐式方案在稳定性方面没有这个限制，但如果时间步长过大，它们就会变得不准确了。

这一条件是由库兰特、弗里德里希斯和卢伊首先认识到的；在1928年，远在计算机可用于数值计算之前！（这一结果出现在德语文章R.库兰特、K.弗里德里希斯和H.卢伊中。这个结果出现在德语文章R.Courant, K. Friedrichs and H. Lewy:<i>&Uuml;ber die partiellen
Differenzengleichungen der mathematischen Physik</i>, Mathematische Annalen, vol. 100, no. 1, pages 32-74, 1928.)这个关于时间步长的条件最常被称为<i>CFL</i>条件。直观地说，CFL条件说的是，时间步长不能大于一个波穿过一个单元的时间。

在程序中，我们将对正方形 $[-1,1]^2$ 均匀地细化七次，得到的网格尺寸为 $h=\frac 1{64}$  ，这就是我们设置的时间步长。我们在两个不同的地方分别设置时间步长和网格尺寸是很容易出错的：很容易再细化一次网格，却忘记同时调整时间步长。   @ref
step_24  "step-24 "显示了一个更好的方法来保持这些东西的同步。




<a name="Thetestcase"></a><h3>The test case</h3>


尽管该程序具有处理非零初始和边界条件以及体力的所有钩子，但我们采取一个简单的案例，即领域是一个正方形 $[-1,1]^2$ ，并且

@f{eqnarray*}
	f &=& 0,
	\\
	u_0 &=& 0,
	\\
	u_1 &=& 0,
	\\
	g &=& \left\{\begin{matrix}\sin (4\pi t)
	&\qquad& \text{for }\ t\le \frac 12, x=-1, -\frac 13<y<\frac 13
	\\
	 0
	&&\text{otherwise}
	\end{matrix}
	\right.


@f}

这相当于一个最初处于静止状态、四周被夹住的膜，有人将夹住的边界的一部分上下挥动一次，从而将波射入领域。


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
 * 我们从通常的各种各样的包含文件开始，我们在以前的许多测试中都看到过。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 这里是仅有的三个有一些新兴趣的包含文件。第一个文件已经被使用了，例如，用于 VectorTools::interpolate_boundary_values 和 MatrixTools::apply_boundary_values 函数。然而，我们在这里使用该类中的另一个函数， VectorTools::project 来计算我们的初始值，作为连续初始值的 $L^2$ 投影。此外，我们使用  VectorTools::create_right_hand_side  来生成积分  $(f^n,\phi^n_i)$  。这些以前总是由 <code>assemble_system</code> 或应用程序代码中的类似函数手工生成。然而，我们太懒了，不能在这里这么做，所以干脆使用库函数。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * @endcode
 * 
 * 与此非常相似，我们也懒得写代码来组装质量矩阵和拉普拉斯矩阵，尽管这只需要从以前的任何一个教程程序中复制相关代码。相反，我们想把重点放在这个程序中真正新的东西上，因此使用了 MatrixCreator::create_mass_matrix 和 MatrixCreator::create_laplace_matrix 函数。它们被声明在这里。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/matrix_tools.h> 
 * 
 * @endcode
 * 
 * 最后，这里有一个include文件，它包含了人们有时需要的各种工具函数。特别是，我们需要 Utilities::int_to_string 类，该类在给定一个整数参数后，返回它的字符串表示。它特别有用，因为它允许第二个参数，表明我们希望结果用前导零填充的数字数。我们将用它来写输出文件，其形式为 <code>solution-XXX.vtu</code> where <code>XXX</code> 表示时间步数，并且总是由三位数组成，即使我们仍然处于个位或两位数的时间步数中。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/utilities.h> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step23 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeWaveEquationcodeclass"></a> 
 * <h3>The <code>WaveEquation</code> class</h3>
 * 

 * 
 * 接下来是主类的声明。它的公共函数接口与其他大多数教程程序一样。值得一提的是，我们现在必须存储四个矩阵，而不是一个：质量矩阵  $M$  ，拉普拉斯矩阵  $A$  ，用于求解  $U^n$  的矩阵  $M+k^2\theta^2A$  ，以及用于求解  $V^n$  的带有边界条件的质量矩阵副本。请注意，在周围有一个额外的质量矩阵副本是有点浪费的。我们将在可能的改进部分讨论如何避免这种情况的策略。
 * 

 * 
 * 同样，我们需要 $U^n,V^n$ 的解向量，以及前一个时间步骤 $U^{n-1},V^{n-1}$ 的相应向量。 <code>system_rhs</code> 将用于我们在每个时间步骤中求解两个线性系统之一时的任何右手向量。这些将在两个函数  <code>solve_u</code>  和  <code>solve_v</code>  中解决。
 * 

 * 
 * 最后，变量 <code>theta</code> 用来表示参数 $\theta$ ，该参数用于定义使用哪种时间步进方案，这在介绍中已经说明。剩下的就不言而喻了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class WaveEquation 
 *   { 
 *   public: 
 *     WaveEquation(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void solve_u(); 
 *     void solve_v(); 
 *     void output_results() const; 
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
 *     SparseMatrix<double> matrix_u; 
 *     SparseMatrix<double> matrix_v; 
 * 
 *     Vector<double> solution_u, solution_v; 
 *     Vector<double> old_solution_u, old_solution_v; 
 *     Vector<double> system_rhs; 
 * 
 *     double       time_step; 
 *     double       time; 
 *     unsigned int timestep_number; 
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
 * 在我们继续填写主类的细节之前，让我们定义与问题相对应的方程数据，即解 $u$ 及其时间导数 $v$ 的初始值和边界值，以及一个右手类。我们使用从Function类模板派生出来的类来做这件事，这个模板以前已经用过很多次了，所以下面的内容不应该是一个惊喜。
 * 

 * 
 * 我们从初始值开始，对数值 $u$ 以及它的时间导数，即速度 $v$ 都选择零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class InitialValuesU : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & /*p*/, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 *       return 0; 
 *     } 
 *   }; 
 * 
 *   template <int dim> 
 *   class InitialValuesV : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & /*p*/, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 *       return 0; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 其次，我们有右手边的强制项。无聊的是，我们在这里也选择零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & /*p*/, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 *       return 0; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 最后，我们有  $u$  和  $v$  的边界值。它们与介绍中描述的一样，一个是另一个的时间导数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BoundaryValuesU : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 * 
 *       if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) && 
 *           (p[1] > -1. / 3)) 
 *         return std::sin(this->get_time() * 4 * numbers::PI); 
 *       else 
 *         return 0; 
 *     } 
 *   }; 
 * 
 *   template <int dim> 
 *   class BoundaryValuesV : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 * 
 *       if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) && 
 *           (p[1] > -1. / 3)) 
 *         return (std::cos(this->get_time() * 4 * numbers::PI) * 4 * numbers::PI); 
 *       else 
 *         return 0; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeWaveEquationcodeclass"></a> 
 * <h3>Implementation of the <code>WaveEquation</code> class</h3>
 * 

 * 
 * 实际逻辑的实现实际上是相当短的，因为我们把组装矩阵和右手边的向量等事情交给了库。其余的实际代码不超过130行，其中相当一部分是可以从以前的例子程序中获取的模板代码（例如，解决线性系统的函数，或生成输出的函数）。
 * 

 * 
 * 我们从构造函数开始（关于时间步长的选择的解释，请参见介绍中关于Courant, Friedrichs, and Lewy的部分）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   WaveEquation<dim>::WaveEquation() 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *     , time_step(1. / 64) 
 *     , time(time_step) 
 *     , timestep_number(1) 
 *     , theta(0.5) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="WaveEquationsetup_system"></a> 
 * <h4>WaveEquation::setup_system</h4>
 * 

 * 
 * 下一个函数是在程序开始时，也就是在第一个时间步骤之前，设置网格、DoFHandler以及矩阵和向量。如果你已经阅读了至少到 step-6 为止的教程程序，那么前几行是相当标准的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WaveEquation<dim>::setup_system() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -1, 1); 
 *     triangulation.refine_global(7); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl; 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl 
 *               << std::endl; 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 * @endcode
 * 
 * 然后，我们必须初始化程序过程中需要的3个矩阵：质量矩阵、拉普拉斯矩阵和在每个时间步长中求解 $M+k^2\theta^2A$ 时使用的矩阵 $U^n$ 。
 * 

 * 
 * 在设置这些矩阵时，请注意它们都是利用了相同的稀疏模式对象。最后，在deal.II中矩阵和稀疏模式是独立对象的原因（与其他许多有限元或线性代数类不同）变得很清楚：在相当一部分应用中，我们必须持有几个恰好具有相同稀疏模式的矩阵，它们没有理由不共享这一信息，而不是重新建立并多次浪费内存。
 * 

 * 
 * 在初始化所有这些矩阵后，我们调用库函数来建立拉普拉斯和质量矩阵。它们所需要的只是一个DoFHandler对象和一个将用于数值积分的正交公式对象。请注意，在许多方面，这些函数比我们通常在应用程序中做的要好，例如，如果一台机器有多个处理器，它们会自动并行构建矩阵：更多信息见WorkStream的文档或 @ref threads "多处理器并行计算 "模块。解决线性系统的矩阵将在run()方法中被填充，因为我们需要在每个时间步长中重新应用边界条件。
 * 

 * 
 * 
 * @code
 *     mass_matrix.reinit(sparsity_pattern); 
 *     laplace_matrix.reinit(sparsity_pattern); 
 *     matrix_u.reinit(sparsity_pattern); 
 *     matrix_v.reinit(sparsity_pattern); 
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler, 
 *                                       QGauss<dim>(fe.degree + 1), 
 *                                       mass_matrix); 
 *     MatrixCreator::create_laplace_matrix(dof_handler, 
 *                                          QGauss<dim>(fe.degree + 1), 
 *                                          laplace_matrix); 
 * 
 * @endcode
 * 
 * 该函数的其余部分用于将矢量大小设置为正确的值。最后一行关闭了悬挂的节点约束对象。由于我们在一个均匀细化的网格上工作，所以不存在或没有计算过约束条件（即没有必要像其他程序那样调用 DoFTools::make_hanging_node_constraints ），但无论如何，我们需要在下面的一个地方进一步设置一个约束对象。
 * 

 * 
 * 
 * @code
 *     solution_u.reinit(dof_handler.n_dofs()); 
 *     solution_v.reinit(dof_handler.n_dofs()); 
 *     old_solution_u.reinit(dof_handler.n_dofs()); 
 *     old_solution_v.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     constraints.close(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationsolve_uandWaveEquationsolve_v"></a> 
 * <h4>WaveEquation::solve_u and WaveEquation::solve_v</h4>
 * 

 * 
 * 接下来的两个函数是解决与  $U^n$  和  $V^n$  的方程有关的线性系统。这两个函数并不特别有趣，因为它们基本沿用了前面所有教程程序中的方案。
 * 

 * 
 * 我们可以对我们要反转的两个矩阵的预处理程序做一些小实验。然而，事实证明，对于这里的矩阵，使用雅可比或SSOR预处理器可以稍微减少解决线性系统所需的迭代次数，但由于应用预处理器的成本，在运行时间方面并不占优势。这也不是什么损失，但让我们保持简单，只做不做。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WaveEquation<dim>::solve_u() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity()); 
 * 
 *     std::cout << "   u-equation: " << solver_control.last_step() 
 *               << " CG iterations." << std::endl; 
 *   } 
 * 
 *   template <int dim> 
 *   void WaveEquation<dim>::solve_v() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity()); 
 * 
 *     std::cout << "   v-equation: " << solver_control.last_step() 
 *               << " CG iterations." << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationoutput_results"></a> 
 * <h4>WaveEquation::output_results</h4>
 * 

 * 
 * 同样地，下面的函数也和我们之前做的差不多。唯一值得一提的是，这里我们使用 Utilities::int_to_string 函数的第二个参数，生成了一个用前导零填充的时间步长的字符串表示，长度为3个字符。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WaveEquation<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution_u, "U"); 
 *     data_out.add_data_vector(solution_v, "V"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
 * 
 * @endcode
 * 
 * 像  step-15  一样，由于我们在每个时间步长写输出（而且我们要解决的系统相对简单），我们指示DataOut使用zlib压缩算法，该算法针对速度而不是磁盘使用进行了优化，因为否则绘制输出会成为一个瓶颈。
 * 

 * 
 * 
 * @code
 *     DataOutBase::VtkFlags vtk_flags; 
 *     vtk_flags.compression_level = 
 *       DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
 *     data_out.set_flags(vtk_flags); 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="WaveEquationrun"></a> 
 * <h4>WaveEquation::run</h4>
 * 

 * 
 * 下面是程序中唯一有趣的功能。它包含了所有时间步骤的循环，但在这之前我们必须设置网格、DoFHandler和矩阵。此外，我们必须以某种方式从初始值开始。为此，我们使用 VectorTools::project 函数，该函数接收一个描述连续函数的对象，并计算该函数在DoFHandler对象所描述的有限元空间的 $L^2$ 投影。没有比这更简单的了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void WaveEquation<dim>::run() 
 *   { 
 *     setup_system(); 
 * 
 *     VectorTools::project(dof_handler, 
 *                          constraints, 
 *                          QGauss<dim>(fe.degree + 1), 
 *                          InitialValuesU<dim>(), 
 *                          old_solution_u); 
 *     VectorTools::project(dof_handler, 
 *                          constraints, 
 *                          QGauss<dim>(fe.degree + 1), 
 *                          InitialValuesV<dim>(), 
 *                          old_solution_v); 
 * 
 * @endcode
 * 
 * 接下来是循环所有的时间步骤，直到我们到达结束时间（本例中为 $T=5$ ）。在每个时间步骤中，我们首先要解决 $U^n$ ，使用方程  $(M^n + k^2\theta^2 A^n)U^n =$  。
 * $(M^{n,n-1} - k^2\theta(1-\theta) A^{n,n-1})U^{n-1} + kM^{n,n-1}V^{n-1} +$  
 * $k\theta \left[k \theta F^n + k(1-\theta) F^{n-1} \right]$  . 请注意，我们在所有的时间步骤中使用相同的网格，因此， $M^n=M^{n,n-1}=M$  和  $A^n=A^{n,n-1}=A$  。因此，我们首先要做的是将 $MU^{n-1} - k^2\theta(1-\theta) AU^{n-1} + kMV^{n-1}$ 和强制项相加，并将结果放入 <code>system_rhs</code> 向量中。(对于这些加法，我们需要在循环之前声明一个临时向量，以避免在每个时间步骤中重复分配内存)。
 * 

 * 
 * 这里需要意识到的是我们如何将时间变量传达给描述右手边的对象：每个从函数类派生出来的对象都有一个时间字段，可以用 Function::set_time 来设置，用 Function::get_time. 来读取。 实质上，使用这种机制，所有空间和时间的函数因此被认为是在某个特定时间评估的空间的函数。这与我们在有限元程序中的典型需求非常吻合，在有限元程序中，我们几乎总是在一个时间步长上工作，而且从来没有发生过，例如，人们想在任何给定的空间位置上为所有时间评估一个时空函数。
 * 

 * 
 * 
 * @code
 *     Vector<double> tmp(solution_u.size()); 
 *     Vector<double> forcing_terms(solution_u.size()); 
 * 
 *     for (; time <= 5; time += time_step, ++timestep_number) 
 *       { 
 *         std::cout << "Time step " << timestep_number << " at t=" << time 
 *                   << std::endl; 
 * 
 *         mass_matrix.vmult(system_rhs, old_solution_u); 
 * 
 *         mass_matrix.vmult(tmp, old_solution_v); 
 *         system_rhs.add(time_step, tmp); 
 * 
 *         laplace_matrix.vmult(tmp, old_solution_u); 
 *         system_rhs.add(-theta * (1 - theta) * time_step * time_step, tmp); 
 * 
 *         RightHandSide<dim> rhs_function; 
 *         rhs_function.set_time(time); 
 *         VectorTools::create_right_hand_side(dof_handler, 
 *                                             QGauss<dim>(fe.degree + 1), 
 *                                             rhs_function, 
 *                                             tmp); 
 *         forcing_terms = tmp; 
 *         forcing_terms *= theta * time_step; 
 * 
 *         rhs_function.set_time(time - time_step); 
 *         VectorTools::create_right_hand_side(dof_handler, 
 *                                             QGauss<dim>(fe.degree + 1), 
 *                                             rhs_function, 
 *                                             tmp); 
 * 
 *         forcing_terms.add((1 - theta) * time_step, tmp); 
 * 
 *         system_rhs.add(theta * time_step, forcing_terms); 
 * 
 * @endcode
 * 
 * 如此构建了第一个方程的右手向量后，我们要做的就是应用正确的边界值。至于右手边，这是一个在特定时间评估的时空函数，我们在边界节点插值，然后像通常那样用结果来应用边界值。然后将结果交给solve_u()函数。
 * 

 * 
 * 
 * @code
 *         { 
 *           BoundaryValuesU<dim> boundary_values_u_function; 
 *           boundary_values_u_function.set_time(time); 
 * 
 *           std::map<types::global_dof_index, double> boundary_values; 
 *           VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                    0, 
 *                                                    boundary_values_u_function, 
 *                                                    boundary_values); 
 * 
 * @endcode
 * 
 * solve_u()的矩阵在每个时间步骤中都是相同的，所以人们可以认为只在模拟开始时做一次就足够了。然而，由于我们需要对线性系统应用边界值（消除了一些矩阵的行和列，并对右手边做出了贡献），在实际应用边界数据之前，我们必须在每个时间步骤中重新填充该矩阵。实际内容非常简单：它是质量矩阵和加权拉普拉斯矩阵的总和。
 * 

 * 
 * 
 * @code
 *           matrix_u.copy_from(mass_matrix); 
 *           matrix_u.add(theta * theta * time_step * time_step, laplace_matrix); 
 *           MatrixTools::apply_boundary_values(boundary_values, 
 *                                              matrix_u, 
 *                                              solution_u, 
 *                                              system_rhs); 
 *         } 
 *         solve_u(); 
 * 
 * @endcode
 * 
 * 第二步，即求解 $V^n$ ，工作原理类似，只是这次左边的矩阵是质量矩阵（我们再次复制，以便能够应用边界条件，而右边是 $MV^{n-1} - k\left[ \theta A U^n + (1-\theta) AU^{n-1}\right]$ 加上强制项。边界值的应用方式与之前相同，只是现在我们必须使用BoundaryValuesV类。
 * 

 * 
 * 
 * @code
 *         laplace_matrix.vmult(system_rhs, solution_u); 
 *         system_rhs *= -theta * time_step; 
 * 
 *         mass_matrix.vmult(tmp, old_solution_v); 
 *         system_rhs += tmp; 
 * 
 *         laplace_matrix.vmult(tmp, old_solution_u); 
 *         system_rhs.add(-time_step * (1 - theta), tmp); 
 * 
 *         system_rhs += forcing_terms; 
 * 
 *         { 
 *           BoundaryValuesV<dim> boundary_values_v_function; 
 *           boundary_values_v_function.set_time(time); 
 * 
 *           std::map<types::global_dof_index, double> boundary_values; 
 *           VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                    0, 
 *                                                    boundary_values_v_function, 
 *                                                    boundary_values); 
 *           matrix_v.copy_from(mass_matrix); 
 *           MatrixTools::apply_boundary_values(boundary_values, 
 *                                              matrix_v, 
 *                                              solution_v, 
 *                                              system_rhs); 
 *         } 
 *         solve_v(); 
 * 
 * @endcode
 * 
 * 最后，在计算完两个解的组成部分后，我们输出结果，计算解中的能量，并在将现在的解移入持有上一个时间步长的解的向量后，继续下一个时间步长。注意函数 SparseMatrix::matrix_norm_square 可以在一个步骤中计算 $\left<V^n,MV^n\right>$ 和 $\left<U^n,AU^n\right>$ ，为我们节省了一个临时向量和几行代码的费用。
 * 

 * 
 * 
 * @code
 *         output_results(); 
 * 
 *         std::cout << "   Total energy: " 
 *                   << (mass_matrix.matrix_norm_square(solution_v) + 
 *                       laplace_matrix.matrix_norm_square(solution_u)) / 
 *                        2 
 *                   << std::endl; 
 * 
 *         old_solution_u = solution_u; 
 *         old_solution_v = solution_v; 
 *       } 
 *   } 
 * } // namespace Step23 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 剩下的就是程序的主要功能了。这里没有什么是在前面几个程序中没有展示过的。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step23; 
 * 
 *       WaveEquation<2> wave_equation_solver; 
 *       wave_equation_solver.run(); 
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
examples/step-23/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当该程序运行时，它产生了以下输出。

@code
Number of active cells: 16384
Number of degrees of freedom: 16641


Time step 1 at t=0.015625
   u-equation: 8 CG iterations.
   v-equation: 22 CG iterations.
   Total energy: 1.17887
Time step 2 at t=0.03125
   u-equation: 8 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 2.9655
Time step 3 at t=0.046875
   u-equation: 8 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 4.33761
Time step 4 at t=0.0625
   u-equation: 7 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 5.35499
Time step 5 at t=0.078125
   u-equation: 7 CG iterations.
   v-equation: 21 CG iterations.
   Total energy: 6.18652
Time step 6 at t=0.09375
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 6.6799


...


Time step 31 at t=0.484375
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 21.9068
Time step 32 at t=0.5
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.3394
Time step 33 at t=0.515625
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019


...


Time step 319 at t=4.98438
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019
Time step 320 at t=5
   u-equation: 7 CG iterations.
   v-equation: 20 CG iterations.
   Total energy: 23.1019
@endcode



我们立即看到的是，至少在 $t=\frac 12$ 之后，能量是一个常数（在此之前，边界源项 $g$ 是非零的，向系统注入能量）。

除了屏幕输出外，程序还将每个时间步骤的解写到输出文件中。如果我们对其进行充分处理，并将其粘贴到电影中，我们会得到以下结果。

 <img src="https://www.dealii.org/images/steps/developer/step-23.movie.gif" alt="Animation of the solution of step 23."> 

影片显示了所产生的波在域中移动并返回，在夹持的边界处被反射。一些数值噪声跟在波的后面，这是由于网格尺寸过大造成的假象，可以通过减小网格宽度和时间步长来减少。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


如果你想探索一下，可以尝试以下一些东西。   <ul>   <li>  Varying  $\theta$  。这给出了不同的时间步进方案，其中一些是稳定的，而另一些则不是。看一看能量是如何演变的。

    <li>  不同的初始和边界条件，右手边。

    <li>  更复杂的域或更精细的网格。请记住，时间步长需要以网格宽度为界，所以改变网格也应该包括改变时间步长。我们将在第24步再来讨论这个问题。

    <li>  可变系数。在现实介质中，波速往往是可变的。特别是，现实介质中的 "真实 "波浪方程将读作@f[
     \rho(x) \frac{\partial^2 u}{\partial t^2}


     -
     \nabla \cdot
     a(x) \nabla u = f,
  @f]

  其中 $\rho(x)$ 是材料的密度， $a(x)$ 与刚度系数有关。那么波速就是  $c=\sqrt{a/\rho}$  。

  为了做出这样的改变，我们将不得不用可变系数来计算质量和拉普拉斯矩阵。幸运的是，这并不难：函数 MatrixCreator::create_laplace_matrix 和 MatrixCreator::create_mass_matrix 有额外的默认参数，可以用来向它们传递非恒定系数函数。因此，所需的变化相对较小。另一方面，必须再次注意确保时间步长在允许范围内。

    <li>  在代码内的注释中，我们讨论了这样一个事实：由于边界条件的原因，用于求解 $U^n$ 和 $V^n$ 的矩阵需要在每次都被重置，尽管实际内容没有变化。可以通过不消除线性系统中的列来避免复制，这可以通过在调用中附加一个 @p false 参数来实现。   @code
    MatrixTools::apply_boundary_values(boundary_values,
                                       matrix_u,
                                       solution_u,
                                       system_rhs,
                                       false);
  @endcode



    <li>  deal.II是一个支持自适应网格的库，如果这个程序支持每隔几步就改变网格，那当然很好。考虑到解决方案的结构&mdash; 一个穿越领域的波浪&mdash; 如果我们只在波浪目前所在的地方完善网格，而不是简单地在所有地方完善网格，这似乎是合适的。直观地看，我们应该能够通过这种方式节省大量的单元。虽然经过进一步的思考，我们意识到这只是在模拟的初始阶段。   一段时间后，对于波浪现象来说，域中充满了初始波的反射，向各个方向发展，充满了域中的每个角落。   在这一点上，一般来说，使用局部网格细化可以获得的好处不多）。)

  为了使自适应改变网格成为可能，基本上有两条路线。   "正确 "的方法是回到我们使用罗特方法得到的弱形式。例如，在每个时间步骤中要解决的两个方程中的第一个方程看起来是这样的。   \f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla u^{n-1},\nabla \varphi)
  +
  k(v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right].
  \f} 现在，注意我们在网格 ${\mathbb T}^n$ 上求解 $u^n$ ，因此测试函数 $\varphi$ 也必须来自空间 $V_h^n$ 。正如在介绍中所讨论的，像 $(u^{n-1},\varphi)$ 这样的条款要求我们将上一步的解决方案（可能是在不同的网格 ${\mathbb T}^{n-1}$ 上计算的）与当前网格的测试函数进行整合，从而得到一个矩阵 $M^{n,n-1}$ 。这个整合不同网格的形状函数的过程，充其量是尴尬的。它是可以做到的，但是因为很难保证 ${\mathbb T}^{n-1}$ 和 ${\mathbb T}^{n}$ 最多只相差一个细化级别，所以我们必须递归匹配两个网格的单元。这样做是可行的，但它会导致冗长的、不完全明显的代码。

  第二种方法如下：每当我们改变网格时，我们只需使用SolutionTransfer类将旧网格上的最后一个时间步长的解内插到新网格上。换句话说，我们将解决\f{eqnarray*}
  (u^n,\varphi) + k^2\theta^2(\nabla u^n,\nabla \varphi) &=&
  (I^n u^{n-1},\varphi) - k^2\theta(1-\theta)(\nabla I^n u^{n-1},\nabla \varphi)
  +
  k(I^n v^{n-1},\varphi)
  + k^2\theta
  \left[
  \theta (f^n,\varphi) + (1-\theta) (f^{n-1},\varphi)
  \right],
  \f}，其中 $I^n$ 将一个给定的函数插值到网格 ${\mathbb T}^n$ ，而不是上面的方程。   这是一个更简单的方法，因为在每个时间步长中，我们不再需要担心 $u^{n-1},v^{n-1}$ 是在我们现在使用的同一个网格上计算的，还是在不同的网格上计算的。因此，代码的唯一变化是增加了一个计算误差的函数，为细化标记单元，设置SolutionTransfer对象，将解转移到新的网格上，并在新的网格上重建矩阵和右手向量。建立矩阵和右手边的函数以及求解器都不需要改变。

  虽然严格来说，这第二种方法在罗特框架中是不太正确的（它引入了一个额外的误差源，即插值），然而这几乎是每个人在解决时间相关方程时所做的事情。我们将在步骤31中使用这种方法，例如。   </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-23.cc"
*/
