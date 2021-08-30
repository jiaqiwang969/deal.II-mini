/**
@page step_52 The step-52 tutorial program
This tutorial depends on step-26.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Problemstatement">Problem statement</a>
        <li><a href="#RungeKuttamethods">Runge-Kutta methods</a>
      <ul>
        <li><a href="#ExplicitRungeKuttamethods">Explicit Runge-Kutta methods</a>
        <li><a href="#EmbeddedRungeKuttamethods">Embedded Runge-Kutta methods</a>
        <li><a href="#ImplicitRungeKuttamethods">Implicit Runge-Kutta methods</a>
      </ul>
        <li><a href="#Spatiallydiscreteformulation">Spatially discrete formulation</a>
        <li><a href="#Notesonthetestcase">Notes on the testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeDiffusioncodeclass">The <code>Diffusion</code> class</a>
      <ul>
        <li><a href="#codeDiffusionassemble_systemcodeintDnablab_icdotnablab_jdboldsymbolrintSigma_ab_ib_jdboldsymbolrintb_ib_jdboldsymbolrcodeinverse_mass_matrixcodeM1"><code>Diffusion::assemble_system</code>}  在这个函数中，我们计算  $-\int D \nabla b_i \cdot \nabla b_j d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$  和质量矩阵  $\int b_i b_j d\boldsymbol{r}$  。然后使用直接求解器对质量矩阵进行反演；然后 <code>inverse_mass_matrix</code> 变量将存储质量矩阵的反值，这样 $M^{-1</a>
        <li><a href="#codeDiffusionget_sourcecode"><code>Diffusion::get_source</code></a>
        <li><a href="#codeDiffusionevaluate_diffusioncode"><code>Diffusion::evaluate_diffusion</code></a>
        <li><a href="#codeDiffusionid_minus_tau_J_inversecode"><code>Diffusion::id_minus_tau_J_inverse</code></a>
        <li><a href="#codeDiffusionoutput_resultscode"><code>Diffusion::output_results</code></a>
        <li><a href="#codeDiffusionexplicit_methodcode"><code>Diffusion::explicit_method</code></a>
        <li><a href="#codeDiffusionimplicit_methodcodecodeexplicit_methodcodeM1ftyleftItauM1fracpartialftypartialyright1"><code>Diffusion::implicit_method</code>}  这个函数等同于 <code>explicit_method</code> ，但用于隐式方法。当使用隐式方法时，我们需要评估 $M^{-1}(f(t,y))$ 和 $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1</a>
        <li><a href="#codeDiffusionembedded_explicit_methodcode"><code>Diffusion::embedded_explicit_method</code></a>
        <li><a href="#codeDiffusionruncode"><code>Diffusion::run</code></a>
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
examples/step-52/doc/intro.dox

 <br> 

<i>This program was contributed by Bruno Turcksin and Damien Lebrun-Grandie.</i>

 @note  为了运行这个程序，deal.II必须被配置为使用UMFPACK稀疏直接解算器。请参考<a
href="../../readme.html#umfpack">ReadMe</a>中的说明如何做到这一点。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序展示了如何使用Runge-Kutta方法来解决一个随时间变化的问题。它解决了首先在步骤26中讨论的热方程的一个小变化，但是由于这个程序的目的只是演示使用更高级的方法与deal.II的时间步进算法对接，所以只解决了一个均匀细化网格上的简单问题。




<a name="Problemstatement"></a><h3>Problem statement</h3>


在这个例子中，我们求解中子输运方程的单组时间依赖性扩散近似（关于时间依赖性多组扩散，见步骤28）。这是一个关于中子如何在高散射介质中移动的模型，因此它是时间依赖性扩散方程的一个变体--它只是步骤26中讨论的热方程的一个不同名称，加上一些额外的条款。我们假设介质是不可逆的，因此，中子通量满足以下方程。

@f{eqnarray*}
\frac{1}{v}\frac{\partial \phi(x,t)}{\partial t} = \nabla \cdot D(x) \nabla \phi(x,t)


- \Sigma_a(x) \phi(x,t) + S(x,t)


@f}

通过适当的边界条件增强。这里， $v$ 是中子的速度（为简单起见，我们假设它等于1，这可以通过简单地缩放时间变量来实现）， $D$ 是扩散系数， $\Sigma_a$ 是吸收截面， $S$ 是一个源。因为我们只对时间依赖性感兴趣，我们假设 $D$ 和 $\Sigma_a$ 是常数。

由于这个程序只打算演示如何使用先进的时间步进算法，我们将只寻找相对简单问题的解。具体来说，我们要在一个正方形域 $[0,b]\times[0,b]$ 上寻找一个解，其形式为

@f{eqnarray*}
\phi(x,t) = A\sin(\omega t)(bx-x^2).


@f}

通过使用二次有限元，我们可以在任何特定时间精确地表示这个函数，所有的误差都是由于时间离散化造成的。我们这样做是因为这样就很容易观察到我们将要考虑的各种时间步进方案的收敛顺序，而不需要将空间和时间误差分开。

我们施加以下边界条件：对 $x=0$ 和 $x=b$ 施加同质的迪里希特条件，对 $y=0$ 和 $y=b$ 施加同质的纽曼条件。我们选择源项，以便相应的解决方案实际上是上述的形式。

@f{eqnarray*}
S=A\left(\frac{1}{v}\omega \cos(\omega t)(bx -x^2) + \sin(\omega t)
\left(\Sigma_a (bx-x^2)+2D\right) \right).


@f}

因为解是时间上的正弦，我们知道精确解满足 $\phi\left(x,\frac{\pi}{\omega}\right) = 0$  。因此，时间 $t=\frac{\pi}{\omega}$ 的误差只是数值解的规范，即 $\|e(\cdot,t=\frac{\pi}{\omega})\|_{L_2} = \|\phi_h(\cdot,t=\frac{\pi}{\omega})\|_{L_2}$ ，而且特别容易评估。在代码中，我们评估 $l_2$ 的节点值的规范，而不是相关空间函数的 $L_2$ 规范，因为前者的计算更简单；然而，在均匀网格上，这两者只是由一个常数相关，因此我们可以用其中一个观察时间收敛顺序。




<a name="RungeKuttamethods"></a><h3>Runge-Kutta methods</h3>


在deal.II中实现的Runge-Kutta方法假定要解决的方程可以写成。

@f{eqnarray*}
\frac{dy}{dt} = g(t,y).


@f}

另一方面，当使用有限元时，离散化的时间导数总是导致左手边存在一个质量矩阵。这可以很容易地看出，如果上述方程中的解向量 $y(t)$ 实际上是节点系数的向量 $U(t)$ ，其形式为变量

@f{eqnarray*}
  u_h(x,t) = \sum_j U_j(t) \varphi_j(x)


@f}

用空间形状函数 $\varphi_j(x)$ ，然后乘以一个形式的方程

@f{eqnarray*}
  \frac{\partial u(x,t)}{\partial t} = q(t,u(x,t))


@f}

通过测试函数，对 $\Omega$ 进行积分，代入 $u\rightarrow u_h$ 并将测试函数限制在上面的 $\varphi_i(x)$ ，那么这个空间离散方程的形式为

@f{eqnarray*}
M\frac{dU}{dt} = f(t,U),


@f}

其中 $M$ 是质量矩阵， $f(t,U)$ 是 $q(t,u(x,t))$ 的空间离散版本（其中 $q$ 通常是出现空间导数的地方，但鉴于我们只考虑时间导数，这一点目前并不太关心）。换句话说，这种形式符合上面的一般方案，如果我们写成

@f{eqnarray*}
\frac{dy}{dt} = g(t,y) = M^{-1}f(t,y).


@f}



Runk-Kutta方法是一种时间步进方案，通过特定的一步法对 $y(t_n)\approx
y_{n}$ 进行近似。它们通常被写成以下形式

@f{eqnarray*}
y_{n+1} = y_n + \sum_{i=1}^s b_i k_i


@f}

其中对于上面的右手边的形式

@f{eqnarray*}
k_i = h M^{-1} f\left(t_n+c_ih,y_n+\sum_{j=1}^sa_{ij}k_j\right).


@f}

这里 $a_{ij}$ ,  $b_i$ , 和 $c_i$ 是已知的系数，确定你要使用的特定Runge-Kutta方案， $h=t_{n+1}-t_n$ 是使用的时间步长。Runge-Kutta类的不同时间步长方法在级数 $s$ 和系数 $a_{ij}$ 、 $b_i$ 和 $c_i$ 上有所不同，但由于可以查找这些系数的表格值，所以很容易实施。这些表格通常被称为Butcher tableaus）。

在编写本教程时，deal.II中实现的方法可分为三类。<ol>  <li>  显式Runge-Kutta；为了使一个方法成为显式，必须在上述定义 $k_i$ 的公式中， $k_i$ 不出现在右侧。换句话说，这些方法必须满足  $a_{ii}=0, i=1,\ldots,s$  。   <li>  嵌入式（或自适应）Runge-Kutta；我们将在下面讨论其特性。   <li>  隐式Runge-Kutta；这类方法需要解决可能是非线性系统的上述阶段 $k_i$ ，即它们至少有 $a_{ii}\neq 0$  个阶段 $i=1,\ldots,s$  。   </ol>  许多众所周知的时间步进方案，人们通常不会将其与Runge或Kutta的名字联系起来，事实上，它们也可以用这些类别来表达。它们往往代表这些系列的最低阶成员。




<a name="ExplicitRungeKuttamethods"></a><h4>Explicit Runge-Kutta methods</h4>


这些方法，只需要一个函数来评估 $M^{-1}f(t,y)$ ，但不需要（作为隐式方法）来解决涉及 $f(t,y)$ 的 $y$ 的方程。与所有显式时间步长方法一样，当选择的时间步长过大时，它们会变得不稳定。

这一类众所周知的方法包括正向欧拉、三阶Runge-Kutta和四阶Runge-Kutta（通常缩写为RK4）。




<a name="EmbeddedRungeKuttamethods"></a><h4>Embedded Runge-Kutta methods</h4>


这些方法同时使用低阶和高阶方法来估计误差，并决定是否需要缩短时间步长或可以增加。术语 "嵌入 "是指低阶方法不需要对函数 $M^{-1}f(\cdot,\cdot)$ 进行额外的评估，而是重复使用那些必须为高阶方法计算的数据。换句话说，它基本上是免费的，而我们得到的误差估计是使用高阶方法的副产品。

这类方法包括Heun-Euler、Bogacki-Shampine、Dormand-Prince（在Matlab中为ode45，通常缩写为RK45，表示这里使用的低阶和高阶方法分别为4阶和5阶Runge-Kutta方法），Fehlberg和Cash-Karp。

在撰写本文时，只有嵌入式的显式方法得到了实现。




<a name="ImplicitRungeKuttamethods"></a><h4>Implicit Runge-Kutta methods</h4>


隐式方法要求在每个（子）时间步中解决 $\alpha y = f(t,y)$ 形式的 $y$ 的（可能是非线性）系统。在内部，这是用牛顿式方法完成的，因此，它们要求用户提供能够评估 $M^{-1}f(t,y)$ 和 $\left(I-\tau M^{-1} \frac{\partial f}{\partial y}\right)^{-1}$ 或等价的 $\left(M - \tau \frac{\partial f}{\partial y}\right)^{-1} M$ 的函数。

这个算子的特殊形式来自于这样一个事实，即每一个牛顿步骤都需要解一个形式的方程

@f{align*}
  \left(M - \tau \frac{\partial f}{\partial y}\right) \Delta y
  = -M h(t,y)


@f}

对于一些（给定的） $h(t,y)$  。无论时间步长如何，隐式方法总是稳定的，但过大的时间步长当然会影响到解的<i>accuracy</i>，即使数值解仍然稳定且有界。

这类方法包括后退欧拉法、隐式中点法、Crank-Nicolson法和两阶段SDIRK法（"单对角隐式Runge-Kutta "的简称，这个术语是用来表示定义时间步进方法的对角线元素 $a_{ii}$ 都是相等的；这个特性使得牛顿矩阵 $I-\tau M^{-1}\frac{\partial f}{\partial y}$ 可以在各阶段之间重复使用，因为 $\tau$ 每次都是相同的）。




<a name="Spatiallydiscreteformulation"></a><h3>Spatially discrete formulation</h3>


通过扩大我们的模型问题的解决方案，一如既往地使用形状函数 $\psi_j$ 并写出

@f{eqnarray*}
\phi_h(x,t) = \sum_j U_j(t) \psi_j(x),


@f}

我们立即得到扩散方程的空间离散化版本为

@f{eqnarray*}
  M \frac{dU(t)}{dt}
  = -{\cal D} U(t) - {\cal A} U(t) + {\cal S}(t)


@f}

其中

@f{eqnarray*}
  M_{ij}  &=& (\psi_i,\psi_j), \\
  {\cal D}_{ij}  &=& (D\nabla\psi_i,\nabla\psi_j)_\Omega, \\
  {\cal A}_{ij}  &=& (\Sigma_a\psi_i,\psi_j)_\Omega, \\
  {\cal S}_{i}(t)  &=& (\psi_i,S(x,t))_\Omega.


@f}

参见第24步和第26步以了解我们如何到达这里。由于当前问题所选择的边界条件，边界项是没有必要的。为了使用Runge-Kutta方法，我们将其改写如下。

@f{eqnarray*}
f(y) = -{\cal D}y - {\cal A}y + {\cal S}.


@f}

在代码中，我们将需要能够评估这个函数 $f(U)$ 以及它的导数。

@f{eqnarray*}
\frac{\partial f}{\partial y} = -{\cal D} - {\cal A}.


@f}






<a name="Notesonthetestcase"></a><h3>Notes on the testcase</h3>


为了简化问题，域是二维的，网格是均匀细化的（不需要调整网格，因为我们使用的是二次有限元，而且精确解是二次的）。从二维域到三维域并不是很有挑战性。然而，如果你打算解决更复杂的问题，必须对网格进行调整（例如在步骤26中），那么就必须记住以下问题。

<ol>  <li>  在改变网格时，你需要将解投影到新的网格上。当然，从每个时间步长的开始到结束，所使用的网格应该是相同的，这个问题的出现是因为Runge-Kutta方法在每个时间步长内使用了多次方程求值。   <li>  每次改变网格时，你都需要更新质量矩阵和它的逆值。   </ol>  这些步骤的技术可以通过查看步骤26轻易获得。


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
 * 像往常一样，第一个任务是包括这些著名的deal.II库文件和一些C++头文件的功能。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/discrete_time.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_out.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * #include <cmath> 
 * #include <map> 
 * 
 * @endcode
 * 
 * 这是唯一一个新的包含文件：它包括所有的Runge-Kutta方法。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/time_stepping.h> 
 * 
 * @endcode
 * 
 * 接下来的步骤与之前所有的教程程序一样。我们把所有的东西放到一个自己的命名空间中，然后把deal.II的类和函数导入其中。
 * 

 * 
 * 
 * @code
 * namespace Step52 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeDiffusioncodeclass"></a> 
 * <h3>The <code>Diffusion</code> class</h3>
 * 

 * 
 * 下一块是主类的声明。这个类中的大部分函数并不新鲜，在以前的教程中已经解释过了。唯一有趣的函数是  <code>evaluate_diffusion()</code>  和  <code>id_minus_tau_J_inverse()</code>. <code>evaluate_diffusion()</code>  评估扩散方程，  $M^{-1}(f(t,y))$  ，在一个给定的时间和一个给定的  $y$  。  <code>id_minus_tau_J_inverse()</code>  在给定的时间和给定的 $\tau$ 和 $y$ 下，评估 $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ 或类似的 $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$ 。当使用隐式方法时，就需要这个函数。
 * 

 * 
 * 
 * @code
 *   class Diffusion 
 *   { 
 *   public: 
 *     Diffusion(); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 * 
 *     void assemble_system(); 
 * 
 *     double get_source(const double time, const Point<2> &point) const; 
 * 
 *     Vector<double> evaluate_diffusion(const double          time, 
 *                                       const Vector<double> &y) const; 
 * 
 *     Vector<double> id_minus_tau_J_inverse(const double          time, 
 *                                           const double          tau, 
 *                                           const Vector<double> &y); 
 * 
 *     void output_results(const double                     time, 
 *                         const unsigned int               time_step, 
 *                         TimeStepping::runge_kutta_method method) const; 
 * 
 * @endcode
 * 
 * 接下来的三个函数分别是显式方法、隐式方法和嵌入式显式方法的驱动。嵌入显式方法的驱动函数返回执行的步数，鉴于它只接受作为参数传递的时间步数作为提示，但内部计算了最佳时间步数本身。
 * 

 * 
 * 
 * @code
 *     void explicit_method(const TimeStepping::runge_kutta_method method, 
 *                          const unsigned int                     n_time_steps, 
 *                          const double                           initial_time, 
 *                          const double                           final_time); 
 * 
 *     void implicit_method(const TimeStepping::runge_kutta_method method, 
 *                          const unsigned int                     n_time_steps, 
 *                          const double                           initial_time, 
 *                          const double                           final_time); 
 * 
 *     unsigned int 
 *     embedded_explicit_method(const TimeStepping::runge_kutta_method method, 
 *                              const unsigned int n_time_steps, 
 *                              const double       initial_time, 
 *                              const double       final_time); 
 * 
 *     const unsigned int fe_degree; 
 * 
 *     const double diffusion_coefficient; 
 *     const double absorption_cross_section; 
 * 
 *     Triangulation<2> triangulation; 
 * 
 *     const FE_Q<2> fe; 
 * 
 *     DoFHandler<2> dof_handler; 
 * 
 *     AffineConstraints<double> constraint_matrix; 
 * 
 *     SparsityPattern sparsity_pattern; 
 * 
 *     SparseMatrix<double> system_matrix; 
 *     SparseMatrix<double> mass_matrix; 
 *     SparseMatrix<double> mass_minus_tau_Jacobian; 
 * 
 *     SparseDirectUMFPACK inverse_mass_matrix; 
 * 
 *     Vector<double> solution; 
 *   }; 
 * 
 * @endcode
 * 
 * 我们选择二次方有限元，并初始化参数。
 * 

 * 
 * 
 * @code
 *   Diffusion::Diffusion() 
 *     : fe_degree(2) 
 *     , diffusion_coefficient(1. / 30.) 
 *     , absorption_cross_section(1.) 
 *     , fe(fe_degree) 
 *     , dof_handler(triangulation) 
 *   {} 
 * 
 * @endcode
 * 
 * 现在，我们创建约束矩阵和稀疏模式。然后，我们初始化这些矩阵和求解向量。
 * 

 * 
 * 
 * @code
 *   void Diffusion::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              1, 
 *                                              Functions::ZeroFunction<2>(), 
 *                                              constraint_matrix); 
 *     constraint_matrix.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraint_matrix); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 *     mass_matrix.reinit(sparsity_pattern); 
 *     mass_minus_tau_Jacobian.reinit(sparsity_pattern); 
 *     solution.reinit(dof_handler.n_dofs()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{<code>Diffusion::assemble_system</code>}  在这个函数中，我们计算  $-\int D \nabla b_i \cdot \nabla b_j d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$  和质量矩阵  $\int b_i b_j d\boldsymbol{r}$  。然后使用直接求解器对质量矩阵进行反演；然后 <code>inverse_mass_matrix</code> 变量将存储质量矩阵的反值，这样 $M^{-1}$ 就可以使用该对象的 <code>vmult()</code> 函数应用于一个矢量。在内部，UMFPACK并没有真正存储矩阵的逆，而是存储它的LU因子；应用矩阵的逆相当于用这两个因子做一次正解和一次逆解，这与应用矩阵的显式逆具有相同的复杂性）。
 * 

 * 
 * 
 * @code
 *   void Diffusion::assemble_system() 
 *   { 
 *     system_matrix = 0.; 
 *     mass_matrix   = 0.; 
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1); 
 * 
 *     FEValues<2> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_gradients | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix      = 0.; 
 *         cell_mass_matrix = 0.; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               { 
 *                 cell_matrix(i, j) += 
 *                   ((-diffusion_coefficient *                // (-D 
 *                       fe_values.shape_grad(i, q_point) *    //  * grad phi_i 
 *                       fe_values.shape_grad(j, q_point)      //  * grad phi_j 
 *                     - absorption_cross_section *            //  -Sigma 
 *                         fe_values.shape_value(i, q_point) * //  * phi_i 
 *                         fe_values.shape_value(j, q_point))  //  * phi_j) 
 *                    * fe_values.JxW(q_point));               // * dx 
 *                 cell_mass_matrix(i, j) += fe_values.shape_value(i, q_point) * 
 *                                           fe_values.shape_value(j, q_point) * 
 *                                           fe_values.JxW(q_point); 
 *               } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         constraint_matrix.distribute_local_to_global(cell_matrix, 
 *                                                      local_dof_indices, 
 *                                                      system_matrix); 
 *         constraint_matrix.distribute_local_to_global(cell_mass_matrix, 
 *                                                      local_dof_indices, 
 *                                                      mass_matrix); 
 *       } 
 * 
 *     inverse_mass_matrix.initialize(mass_matrix); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionget_sourcecode"></a> 
 * <h4><code>Diffusion::get_source</code></h4>
 * 

 * 
 * 在这个函数中，计算出特定时间和特定点的方程的源项。
 * 

 * 
 * 
 * @code
 *   double Diffusion::get_source(const double time, const Point<2> &point) const 
 *   { 
 *     const double intensity = 10.; 
 *     const double frequency = numbers::PI / 10.; 
 *     const double b         = 5.; 
 *     const double x         = point(0); 
 * 
 *     return intensity * 
 *            (frequency * std::cos(frequency * time) * (b * x - x * x) + 
 *             std::sin(frequency * time) * 
 *               (absorption_cross_section * (b * x - x * x) + 
 *                2. * diffusion_coefficient)); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionevaluate_diffusioncode"></a> 
 * <h4><code>Diffusion::evaluate_diffusion</code></h4>
 * 

 * 
 * 接下来，我们在给定的时间  $t$  和给定的矢量  $y$  评价扩散方程的弱形式。换句话说，正如介绍中所述，我们评估  $M^{-1}(-{\cal D}y - {\cal A}y + {\cal S})$  。为此，我们必须将矩阵 $-{\cal D} - {\cal A}$ （之前计算并存储在变量 <code>system_matrix</code> 中）应用于 $y$ ，然后添加源项，我们像通常那样进行积分。(如果你想节省几行代码，或者想利用并行积分的优势，可以用 VectorTools::create_right_hand_side() 来进行积分。) 然后将结果乘以 $M^{-1}$  。
 * 

 * 
 * 
 * @code
 *   Vector<double> Diffusion::evaluate_diffusion(const double          time, 
 *                                                const Vector<double> &y) const 
 *   { 
 *     Vector<double> tmp(dof_handler.n_dofs()); 
 *     tmp = 0.; 
 *     system_matrix.vmult(tmp, y); 
 * 
 *     const QGauss<2> quadrature_formula(fe_degree + 1); 
 * 
 *     FEValues<2> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_quadrature_points | 
 *                             update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     Vector<double> cell_source(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_source = 0.; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           { 
 *             const double source = 
 *               get_source(time, fe_values.quadrature_point(q_point)); 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               cell_source(i) += fe_values.shape_value(i, q_point) * // phi_i(x) 
 *                                 source *                            // * S(x) 
 *                                 fe_values.JxW(q_point);             // * dx 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         constraint_matrix.distribute_local_to_global(cell_source, 
 *                                                      local_dof_indices, 
 *                                                      tmp); 
 *       } 
 * 
 *     Vector<double> value(dof_handler.n_dofs()); 
 *     inverse_mass_matrix.vmult(value, tmp); 
 * 
 *     return value; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionid_minus_tau_J_inversecode"></a> 
 * <h4><code>Diffusion::id_minus_tau_J_inverse</code></h4>
 * 

 * 
 * 我们计算  $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$  。这要分几个步骤进行。 
 * 

 * 
 * - 计算  $M-\tau \frac{\partial f}{\partial y}$  。  
 * 

 * 
 * - 反转矩阵，得到  $\left(M-\tau \frac{\partial f} {\partial y}\right)^{-1}$  。  
 * 

 * 
 * --计算 $tmp=My$ 。  
 * 

 * 
 * --计算 $z=\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} tmp =  \left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} My$ 。  
 * 

 * 
 * - 返回z。
 * 

 * 
 * 
 * @code
 *   Vector<double> Diffusion::id_minus_tau_J_inverse(const double /*time*/, 
 *                                                    const double          tau, 
 *                                                    const Vector<double> &y) 
 *   { 
 *     SparseDirectUMFPACK inverse_mass_minus_tau_Jacobian; 
 * 
 *     mass_minus_tau_Jacobian.copy_from(mass_matrix); 
 *     mass_minus_tau_Jacobian.add(-tau, system_matrix); 
 * 
 *     inverse_mass_minus_tau_Jacobian.initialize(mass_minus_tau_Jacobian); 
 * 
 *     Vector<double> tmp(dof_handler.n_dofs()); 
 *     mass_matrix.vmult(tmp, y); 
 * 
 *     Vector<double> result(y); 
 *     inverse_mass_minus_tau_Jacobian.vmult(result, tmp); 
 * 
 *     return result; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionoutput_resultscode"></a> 
 * <h4><code>Diffusion::output_results</code></h4>
 * 

 * 
 * 下面的函数将解决方案以vtu文件的形式输出，并以时间步长和时间步长方法的名称为索引。当然，所有时间步长方法的（精确）结果应该是一样的，但这里的输出至少可以让我们对它们进行比较。
 * 

 * 
 * 
 * @code
 *   void Diffusion::output_results(const double                     time, 
 *                                  const unsigned int               time_step, 
 *                                  TimeStepping::runge_kutta_method method) const 
 *   { 
 *     std::string method_name; 
 * 
 *     switch (method) 
 *       { 
 *         case TimeStepping::FORWARD_EULER: 
 *           { 
 *             method_name = "forward_euler"; 
 *             break; 
 *           } 
 *         case TimeStepping::RK_THIRD_ORDER: 
 *           { 
 *             method_name = "rk3"; 
 *             break; 
 *           } 
 *         case TimeStepping::RK_CLASSIC_FOURTH_ORDER: 
 *           { 
 *             method_name = "rk4"; 
 *             break; 
 *           } 
 *         case TimeStepping::BACKWARD_EULER: 
 *           { 
 *             method_name = "backward_euler"; 
 *             break; 
 *           } 
 *         case TimeStepping::IMPLICIT_MIDPOINT: 
 *           { 
 *             method_name = "implicit_midpoint"; 
 *             break; 
 *           } 
 *         case TimeStepping::SDIRK_TWO_STAGES: 
 *           { 
 *             method_name = "sdirk"; 
 *             break; 
 *           } 
 *         case TimeStepping::HEUN_EULER: 
 *           { 
 *             method_name = "heun_euler"; 
 *             break; 
 *           } 
 *         case TimeStepping::BOGACKI_SHAMPINE: 
 *           { 
 *             method_name = "bocacki_shampine"; 
 *             break; 
 *           } 
 *         case TimeStepping::DOPRI: 
 *           { 
 *             method_name = "dopri"; 
 *             break; 
 *           } 
 *         case TimeStepping::FEHLBERG: 
 *           { 
 *             method_name = "fehlberg"; 
 *             break; 
 *           } 
 *         case TimeStepping::CASH_KARP: 
 *           { 
 *             method_name = "cash_karp"; 
 *             break; 
 *           } 
 *         default: 
 *           { 
 *             break; 
 *           } 
 *       } 
 * 
 *     DataOut<2> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     data_out.set_flags(DataOutBase::VtkFlags(time, time_step)); 
 * 
 *     const std::string filename = "solution_" + method_name + "-" + 
 *                                  Utilities::int_to_string(time_step, 3) + 
 *                                  ".vtu"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 * 
 *     static std::vector<std::pair<double, std::string>> times_and_names; 
 * 
 *     static std::string method_name_prev = ""; 
 *     static std::string pvd_filename; 
 *     if (method_name_prev != method_name) 
 *       { 
 *         times_and_names.clear(); 
 *         method_name_prev = method_name; 
 *         pvd_filename     = "solution_" + method_name + ".pvd"; 
 *       } 
 *     times_and_names.emplace_back(time, filename); 
 *     std::ofstream pvd_output(pvd_filename); 
 *     DataOutBase::write_pvd_record(pvd_output, times_and_names); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionexplicit_methodcode"></a> 
 * <h4><code>Diffusion::explicit_method</code></h4>
 * 

 * 
 * 这个函数是所有显式方法的驱动。在顶部，它初始化了时间步长和解决方案（通过将其设置为零，然后确保边界值和悬挂节点约束得到尊重；当然，对于我们在这里使用的网格，悬挂节点约束实际上不是一个问题）。然后调用 <code>evolve_one_time_step</code> ，执行一个时间步骤。时间是通过DiscreteTime对象来存储和增加的。
 * 

 * 
 * 对于显式方法， <code>evolve_one_time_step</code> 需要评估 $M^{-1}(f(t,y))$ ，也就是说，它需要 <code>evaluate_diffusion</code>  。因为 <code>evaluate_diffusion</code> 是一个成员函数，它需要被绑定到 <code>this</code> 。在每个进化步骤之后，我们再次应用正确的边界值和悬挂节点约束。
 * 

 * 
 * 最后，每隔10个时间步骤就会输出解决方案。
 * 

 * 
 * 
 * @code
 *   void Diffusion::explicit_method(const TimeStepping::runge_kutta_method method, 
 *                                   const unsigned int n_time_steps, 
 *                                   const double       initial_time, 
 *                                   const double       final_time) 
 *   { 
 *     const double time_step = 
 *       (final_time - initial_time) / static_cast<double>(n_time_steps); 
 * 
 *     solution = 0.; 
 *     constraint_matrix.distribute(solution); 
 * 
 *     TimeStepping::ExplicitRungeKutta<Vector<double>> explicit_runge_kutta( 
 *       method); 
 *     output_results(initial_time, 0, method); 
 *     DiscreteTime time(initial_time, final_time, time_step); 
 *     while (time.is_at_end() == false) 
 *       { 
 *         explicit_runge_kutta.evolve_one_time_step( 
 *           [this](const double time, const Vector<double> &y) { 
 *             return this->evaluate_diffusion(time, y); 
 *           }, 
 *           time.get_current_time(), 
 *           time.get_next_step_size(), 
 *           solution); 
 *         time.advance_time(); 
 * 
 *         constraint_matrix.distribute(solution); 
 * 
 *         if (time.get_step_number() % 10 == 0) 
 *           output_results(time.get_current_time(), 
 *                          time.get_step_number(), 
 *                          method); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{<code>Diffusion::implicit_method</code>}  这个函数等同于 <code>explicit_method</code> ，但用于隐式方法。当使用隐式方法时，我们需要评估 $M^{-1}(f(t,y))$ 和 $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ ，为此我们使用之前介绍的两个成员函数。
 * 

 * 
 * 
 * @code
 *   void Diffusion::implicit_method(const TimeStepping::runge_kutta_method method, 
 *                                   const unsigned int n_time_steps, 
 *                                   const double       initial_time, 
 *                                   const double       final_time) 
 *   { 
 *     const double time_step = 
 *       (final_time - initial_time) / static_cast<double>(n_time_steps); 
 * 
 *     solution = 0.; 
 *     constraint_matrix.distribute(solution); 
 * 
 *     TimeStepping::ImplicitRungeKutta<Vector<double>> implicit_runge_kutta( 
 *       method); 
 *     output_results(initial_time, 0, method); 
 *     DiscreteTime time(initial_time, final_time, time_step); 
 *     while (time.is_at_end() == false) 
 *       { 
 *         implicit_runge_kutta.evolve_one_time_step( 
 *           [this](const double time, const Vector<double> &y) { 
 *             return this->evaluate_diffusion(time, y); 
 *           }, 
 *           [this](const double time, const double tau, const Vector<double> &y) { 
 *             return this->id_minus_tau_J_inverse(time, tau, y); 
 *           }, 
 *           time.get_current_time(), 
 *           time.get_next_step_size(), 
 *           solution); 
 *         time.advance_time(); 
 * 
 *         constraint_matrix.distribute(solution); 
 * 
 *         if (time.get_step_number() % 10 == 0) 
 *           output_results(time.get_current_time(), 
 *                          time.get_step_number(), 
 *                          method); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{<code>Diffusion::embedded_explicit_method</code>}  这个函数是嵌入式显式方法的驱动。它需要更多的参数。 
 * 

 * 
 * - coarsen_param：当误差低于阈值时，乘以当前时间步长的系数。 
 * 

 * 
 * - refine_param: 当误差高于阈值时，乘以当前时间步长的系数。 
 * 

 * 
 * - min_delta: 可接受的最小时间步长。 
 * 

 * 
 * - max_delta: 可接受的最大时间步长。 
 * 

 * 
 * - refine_tol：时间步长超过的阈值。 
 * 

 * 
 * - coarsen_tol：阈值，低于该阈值的时间步长将被粗化。
 * 

 * 
 * 嵌入方法使用一个猜测的时间步长。如果使用这个时间步长的误差太大，时间步长将被缩小。如果误差低于阈值，则在下一个时间步长时将尝试更大的时间步长。  <code>delta_t_guess</code> 是由嵌入式方法产生的猜测的时间步长。总之，时间步长有可能以三种方式修改。 
 * 

 * 
 * - 在 TimeStepping::EmbeddedExplicitRungeKutta::evolve_one_time_step(). 内减少或增加时间步长。  
 * 

 * 
 * - 使用计算出的  <code>delta_t_guess</code>  。 
 * 

 * 
 * - 自动调整最后一个时间步长，以确保模拟在  <code>final_time</code>  处精确结束。这种调整是在DiscreteTime实例中处理的。
 * 

 * 
 * 
 * @code
 *   unsigned int Diffusion::embedded_explicit_method( 
 *     const TimeStepping::runge_kutta_method method, 
 *     const unsigned int                     n_time_steps, 
 *     const double                           initial_time, 
 *     const double                           final_time) 
 *   { 
 *     const double time_step = 
 *       (final_time - initial_time) / static_cast<double>(n_time_steps); 
 *     const double coarsen_param = 1.2; 
 *     const double refine_param  = 0.8; 
 *     const double min_delta     = 1e-8; 
 *     const double max_delta     = 10 * time_step; 
 *     const double refine_tol    = 1e-1; 
 *     const double coarsen_tol   = 1e-5; 
 * 
 *     solution = 0.; 
 *     constraint_matrix.distribute(solution); 
 * 
 *     TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> 
 *       embedded_explicit_runge_kutta(method, 
 *                                     coarsen_param, 
 *                                     refine_param, 
 *                                     min_delta, 
 *                                     max_delta, 
 *                                     refine_tol, 
 *                                     coarsen_tol); 
 *     output_results(initial_time, 0, method); 
 *     DiscreteTime time(initial_time, final_time, time_step); 
 *     while (time.is_at_end() == false) 
 *       { 
 *         const double new_time = 
 *           embedded_explicit_runge_kutta.evolve_one_time_step( 
 *             [this](const double time, const Vector<double> &y) { 
 *               return this->evaluate_diffusion(time, y); 
 *             }, 
 *             time.get_current_time(), 
 *             time.get_next_step_size(), 
 *             solution); 
 *         time.set_next_step_size(new_time - time.get_current_time()); 
 *         time.advance_time(); 
 * 
 *         constraint_matrix.distribute(solution); 
 * 
 *         if (time.get_step_number() % 10 == 0) 
 *           output_results(time.get_current_time(), 
 *                          time.get_step_number(), 
 *                          method); 
 * 
 *         time.set_desired_next_step_size( 
 *           embedded_explicit_runge_kutta.get_status().delta_t_guess); 
 *       } 
 * 
 *     return time.get_step_number(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="codeDiffusionruncode"></a> 
 * <h4><code>Diffusion::run</code></h4>
 * 

 * 
 * 下面是该程序的主要功能。在顶部，我们创建网格（一个[0,5]x[0,5]的正方形）并对其进行四次细化，得到一个有16乘16单元的网格，共256个。 然后我们将边界指示器设置为1，用于边界中 $x=0$ 和 $x=5$ 的部分。
 * 

 * 
 * 
 * @code
 *   void Diffusion::run() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, 0., 5.); 
 *     triangulation.refine_global(4); 
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->at_boundary()) 
 *           { 
 *             if ((face->center()[0] == 0.) || (face->center()[0] == 5.)) 
 *               face->set_boundary_id(1); 
 *             else 
 *               face->set_boundary_id(0); 
 *           } 
 * 
 * @endcode
 * 
 * 接下来，我们设置线性系统并为其填充内容，以便在整个时间步进过程中使用它们。
 * 

 * 
 * 
 * @code
 *     setup_system(); 
 * 
 *     assemble_system(); 
 * 
 * @endcode
 * 
 * 最后，我们使用命名空间TimeStepping中实现的几种Runge-Kutta方法来解决扩散问题，每次都会在结束时输出误差。(正如介绍中所解释的，由于精确解在最后时间为零，所以误差等于数值解，只需取解向量的 $l_2$ 准则即可计算出来。)
 * 

 * 
 * 
 * @code
 *     unsigned int       n_steps      = 0; 
 *     const unsigned int n_time_steps = 200; 
 *     const double       initial_time = 0.; 
 *     const double       final_time   = 10.; 
 * 
 *     std::cout << "Explicit methods:" << std::endl; 
 *     explicit_method(TimeStepping::FORWARD_EULER, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Forward Euler:            error=" << solution.l2_norm() 
 *               << std::endl; 
 * 
 *     explicit_method(TimeStepping::RK_THIRD_ORDER, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Third order Runge-Kutta:  error=" << solution.l2_norm() 
 *               << std::endl; 
 * 
 *     explicit_method(TimeStepping::RK_CLASSIC_FOURTH_ORDER, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Fourth order Runge-Kutta: error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << std::endl; 
 * 
 *     std::cout << "Implicit methods:" << std::endl; 
 *     implicit_method(TimeStepping::BACKWARD_EULER, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Backward Euler:           error=" << solution.l2_norm() 
 *               << std::endl; 
 * 
 *     implicit_method(TimeStepping::IMPLICIT_MIDPOINT, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Implicit Midpoint:        error=" << solution.l2_norm() 
 *               << std::endl; 
 * 
 *     implicit_method(TimeStepping::CRANK_NICOLSON, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   Crank-Nicolson:           error=" << solution.l2_norm() 
 *               << std::endl; 
 * 
 *     implicit_method(TimeStepping::SDIRK_TWO_STAGES, 
 *                     n_time_steps, 
 *                     initial_time, 
 *                     final_time); 
 *     std::cout << "   SDIRK:                    error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << std::endl; 
 * 
 *     std::cout << "Embedded explicit methods:" << std::endl; 
 *     n_steps = embedded_explicit_method(TimeStepping::HEUN_EULER, 
 *                                        n_time_steps, 
 *                                        initial_time, 
 *                                        final_time); 
 *     std::cout << "   Heun-Euler:               error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << "                   steps performed=" << n_steps << std::endl; 
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::BOGACKI_SHAMPINE, 
 *                                        n_time_steps, 
 *                                        initial_time, 
 *                                        final_time); 
 *     std::cout << "   Bogacki-Shampine:         error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << "                   steps performed=" << n_steps << std::endl; 
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::DOPRI, 
 *                                        n_time_steps, 
 *                                        initial_time, 
 *                                        final_time); 
 *     std::cout << "   Dopri:                    error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << "                   steps performed=" << n_steps << std::endl; 
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::FEHLBERG, 
 *                                        n_time_steps, 
 *                                        initial_time, 
 *                                        final_time); 
 *     std::cout << "   Fehlberg:                 error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << "                   steps performed=" << n_steps << std::endl; 
 * 
 *     n_steps = embedded_explicit_method(TimeStepping::CASH_KARP, 
 *                                        n_time_steps, 
 *                                        initial_time, 
 *                                        final_time); 
 *     std::cout << "   Cash-Karp:                error=" << solution.l2_norm() 
 *               << std::endl; 
 *     std::cout << "                   steps performed=" << n_steps << std::endl; 
 *   } 
 * } // namespace Step52 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * 下面的 <code>main</code> 函数与前面的例子类似，不需要注释。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step52::Diffusion diffusion; 
 *       diffusion.run(); 
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
 *     }; 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-52/doc/results.dox



<a name="Results"></a><h1>Results</h1>


这个程序的重点不在于显示特定的结果，而在于显示它是如何做到的。这一点我们已经通过讨论上面的代码证明过了。因此，该程序的输出相对较少，只包括控制台输出和用于可视化的VTU格式的解决方案。

控制台输出既包含错误，也包含对某些方法所执行的步骤数量。

@code
Explicit methods:
   Forward Euler:            error=1.00883
   Third order Runge-Kutta:  error=0.000227982
   Fourth order Runge-Kutta: error=1.90541e-06


Implicit methods:
   Backward Euler:           error=1.03428
   Implicit Midpoint:        error=0.00862702
   Crank-Nicolson:           error=0.00862675
   SDIRK:                    error=0.0042349


Embedded explicit methods:
   Heun-Euler:               error=0.0073012
                   steps performed=284
   Bogacki-Shampine:         error=0.000408407
                   steps performed=181
   Dopri:                    error=0.000836695
                   steps performed=120
   Fehlberg:                 error=0.00248922
                   steps performed=106
   Cash-Karp:                error=0.0787735
                   steps performed=106
@endcode



正如预期的那样，高阶方法给出了（更）准确的解决方案。我们还看到，（相当不准确的）Heun-Euler方法增加了时间步数，以满足公差要求。另一方面，其他嵌入式方法使用的时间步数比规定的要少得多。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-52.cc"
*/
