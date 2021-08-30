/**
@page step_24 The step-24 tutorial program
This tutorial depends on step-23.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Theproblem">The problem</a>
        <li><a href="#Weakformanddiscretization">Weak form and discretization</a>
        <li><a href="#Whattheprogramdoes">What the program does</a>
        <li><a href="#AppendixPDEswithDiracdeltafunctionsasrighthandsideandtheirtransformationtoaninitialvalueproblem">Appendix: PDEs with Dirac delta functions as right hand side and their transformation to an initial value problem</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Theforwardproblemclasstemplate">The "forward problem" class template</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeTATForwardProblemcodeclass">Implementation of the <code>TATForwardProblem</code> class</a>
      <ul>
        <li><a href="#TATForwardProblemsetup_system">TATForwardProblem::setup_system</a>
        <li><a href="#TATForwardProblemsolve_pandTATForwardProblemsolve_v">TATForwardProblem::solve_p and TATForwardProblem::solve_v</a>
        <li><a href="#TATForwardProblemoutput_results">TATForwardProblem::output_results</a>
        <li><a href="#TATForwardProblemrun">TATForwardProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Oneabsorber"> One absorber </a>
        <li><a href="#Multipleabsorbers">Multiple absorbers</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-24/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个项目是由德克萨斯A&amp;M大学的Xing Jin的一个学生项目发展而来。本程序的大部分工作是由她完成的。这个教程程序的部分工作得到了美国国家科学基金会DMS-0604778号拨款的资助。

该计划是一个旨在模拟热声断层成像的项目的一部分。在热声断层成像中，脉冲电磁能量被送入生物问题。组织吸收一些这种能量，组织中吸收能量最多的那些部分通过热弹性膨胀产生热声波。对于成像来说，人们利用不同种类的组织，最重要的是健康和病变组织，吸收不同数量的能量，因此以不同的速度膨胀。实验装置是测量这些源在组织表面产生的压力波的振幅，并试图重建源的分布，这对吸收器的分布有指示作用，因此对不同种类的组织有指示作用。这个项目的一部分是将模拟数据与实际测量进行比较，因此必须解决 "正向问题"，即描述压力波在组织中传播的波浪方程。因此，这个程序是 @ref
step_23 "step-23 "的延续，其中首次介绍了波浪方程。




<a name="Theproblem"></a><h3>The problem</h3>


在忽略热扩散的情况下，某个位置的温度可以表示为

@f[
\rho C_p \frac{\partial}{\partial t}T(t,\mathbf r) = H(t,\mathbf r)


@f]



这里 $\rho (\mathbf r) $ 是密度； $C_p (\mathbf r) $ 是比热； $\frac{\partial T}{\partial t}(t,\mathbf r)$ 是由于传递的微波能量引起的温升； $H(t,\mathbf r)$ 是加热函数，定义为由沉积的微波能量转化的每一时间和体积的热能。

让我们假设组织具有异质的介电特性，但具有同质的声学特性。在声学均质介质中的基本声学生成方程可以描述如下：如果 $u$ 是矢量值的位移，那么组织肯定通过加速度对压力的变化做出反应。

@f[
\rho \frac{\partial^2}{\partial t^2}u(t,\mathbf r) =


-\nabla p(t,\mathbf r).


@f]

此外，它因压力过大而收缩，并根据温度的变化而膨胀。

@f[
\nabla \cdot u(t,\mathbf r) = -\frac{p(t,\mathbf r)}{\rho c_0^2}+\beta T(t,\mathbf r) .


@f]

这里， $\beta$ 是一个热膨胀系数。

现在让我们假设，加热只发生在比波在组织中传播短得多的时间尺度上（即加热组织的微波脉冲的时间长度远短于波穿过领域的时间）。在这种情况下，加热率 $H(t,\mathbf r)$ 可以写成 $H(t,\mathbf r) = a(\mathbf
r)\delta(t)$ （其中 $a(\mathbf r)$ 是微波能量的吸收强度图， $\delta(t)$ 是狄拉克三角函数），与上述第一个方程一起将产生温度 $T(\mathbf r)$ 在时间 $t=0$ 的瞬时跳跃。利用这一假设，并将所有方程放在一起，我们可以将上述内容重写并合并为以下内容。

@f[
\Delta p-\frac{1}{c_0^2} \frac{\partial^2 p}{\partial t^2} = \lambda
a(\mathbf r)\frac{d\delta(t)}{dt}


@f]

其中 $\lambda = - \frac{\beta}{C_p}$  。

这个有点奇怪的方程，右边是狄拉克三角函数的导数，可以重写为一个初值问题，如下所示。

@f{eqnarray*}
\Delta \bar{p}- \frac{1}{c_0^2} \frac{\partial^2 \bar{p}}{\partial t^2} & = &
0 \\
\bar{p}(0,\mathbf r) &=& c_0^2 \lambda a(\mathbf r) = b(\mathbf r)  \\
\frac{\partial\bar{p}(0,\mathbf r)}{\partial t} &=& 0.


@f}

(在本引言的最后，作为附录给出了这种转化为初值问题的推导)。

在逆向问题中，人们希望恢复的是初始条件 $b(\mathbf r) = c_0^2 \lambda a(\mathbf r)$ ，因为它是微波能量的吸收强度图，因此可能是分辨健康和病变组织的指标。

在实际应用中，热声源相对于介质来说是非常小的。  因此，热声波的传播路径可以被近似为从源头到无限远。此外，检测器离源头只有有限的距离。我们只需要评估热声波通过检测器时的数值，尽管它们确实继续超出。因此，这是一个我们只对无限介质的一小部分感兴趣的问题，我们不希望某个地方产生的波在我们认为有趣的领域的边界上被反射。相反，我们希望只模拟包含在感兴趣的领域内的那部分波场，而碰到该领域边界的波则不受干扰地通过边界。换句话说，我们希望边界能吸收撞击它的任何波。

一般来说，这是一个困难的问题：好的吸收边界条件是非线性的和/或数值上非常昂贵。因此，我们选择了一个简单的一阶近似吸收边界条件，其内容为

@f[
\frac{\partial\bar{p}}{\partial\mathbf n} =


-\frac{1}{c_0} \frac{\partial\bar{p}}{\partial t}


@f]

这里， $\frac{\partial\bar{p}}{\partial\mathbf n}$ 是边界处的法向导数。应该指出的是，这不是一个特别好的边界条件，但它是少数几个合理简单的实现条件之一。




<a name="Weakformanddiscretization"></a><h3>Weak form and discretization</h3>


如同步骤23，首先引入第二个变量，定义为压力势的导数。

@f[
v = \frac{\partial\bar{p}}{\partial t}


@f]



有了第二个变量，我们就可以将正向问题转化为两个独立的方程式。

@f{eqnarray*}
\bar{p}_{t} - v & = & 0 \\
\Delta\bar{p} - \frac{1}{c_0^2}\,v_{t} & = & f


@f}

具有初始条件。

@f{eqnarray*}
\bar{p}(0,\mathbf r) & = & b(r) \\
v(0,\mathbf r)=\bar{p}_t(0,\mathbf r) & = & 0.


@f}

注意，我们在这里引入了一个右手边 $f(t,\mathbf r)$ ，以显示如何在一般情况下推导这些公式，尽管在应用于热声问题时 $f=0$  。

然后，使用步骤23中介绍的一般 $\theta$ 方案，这个模型的半具体化、弱化版本是。

@f{eqnarray*}
\left(\frac{\bar{p}^n-\bar{p}^{n-1}}{k},\phi\right)_\Omega-
\left(\theta v^{n}+(1-\theta)v^{n-1},\phi\right)_\Omega & = & 0   \\


-\left(\nabla((\theta\bar{p}^n+(1-\theta)\bar{p}^{n-1})),\nabla\phi\right)_\Omega-
\frac{1}{c_0}\left(\frac{\bar{p}^n-\bar{p}^{n-1}}{k},\phi\right)_{\partial\Omega} -
\frac{1}{c_0^2}\left(\frac{v^n-v^{n-1}}{k},\phi\right)_\Omega & =
& \left(\theta f^{n}+(1-\theta)f^{n-1}, \phi\right)_\Omega,


@f}

其中 $\phi$ 是一个任意的测试函数，我们使用了吸收边界条件来进行部分积分：吸收边界条件通过使用以下方法被纳入到弱形式之中

@f[
\int_\Omega\varphi \, \Delta p\; dx =


-\int_\Omega\nabla \varphi \cdot \nabla p dx +
\int_{\partial\Omega}\varphi \frac{\partial p}{\partial {\mathbf n}}ds.


@f]



由此，我们通过引入有限数量的形状函数得到离散模型，并得到

@f{eqnarray*}
M\bar{p}^{n}-k \theta M v^n & = & M\bar{p}^{n-1}+k (1-\theta)Mv^{n-1},\\


(-c_0^2k \theta A-c_0 B)\bar{p}^n-Mv^{n} & = &
(c_0^2k(1-\theta)A-c_0B)\bar{p}^{n-1}-Mv^{n-1}+c_0^2k(\theta F^{n}+(1-\theta)F^{n-1}).


@f}

这里的矩阵 $M$ 和 $A$ 与步骤23相同，而边界质量矩阵

@f[
	B_{ij} = \left(\varphi_i,\varphi_j\right)_{\partial\Omega}


@f]

是使用吸收性边界条件的结果。

以上两个方程可以用矩阵形式重写，压力和它的导数是一个未知矢量。

@f[
\left(\begin{array}{cc}
 M         &       -k\theta M \\
c_0^2\,k\,\theta\,A+c_0\,B  &  M   \\
               \end{array} \right)\\
\left(\begin{array}{c}
 \bar{p}^{n}    \\
 \bar{v}^{n}
              \end{array}\right)=\\
\left(\begin{array}{l}
 G_1  \\
 G_2 -(\theta F^{n}+(1-\theta)F ^{n-1})c_{0}^{2}k \\
                \end{array}\right)


@f]



其中

@f[
\left(\begin{array}{c}
G_1 \\
G_2 \\
   \end{array} \right)=\\
\left(\begin{array}{l}
 M\bar{p}^{n-1}+k(1-\theta)Mv^{n-1}\\
 (-c_{0}^{2}k (1-\theta)A+c_0 B)\bar{p}^{n-1} +Mv^{n-1}
                \end{array}\right)


@f]



通过简单的转换，就可以得到压力势及其导数的两个方程，就像前面的教程程序一样。

@f{eqnarray*}
(M+(k\,\theta\,c_{0})^{2}A+c_0k\theta B)\bar{p}^{n} & = &
G_{1}+(k\, \theta)G_{2}-(c_0k)^2\theta (\theta F^{n}+(1-\theta)F^{n-1}) \\
Mv^n & = & -(c_0^2\,k\, \theta\, A+c_0B)\bar{p}^{n}+ G_2 -
c_0^2k(\theta F^{n}+(1-\theta)F^{n-1})


@f}






<a name="Whattheprogramdoes"></a><h3>What the program does</h3>


与Step-23相比，本程序增加了对简单吸收边界条件的处理。此外，它还处理了从实际实验测量得到的数据。为此，我们需要在实验也评估了真实压力场的点上评估解决方案。我们将看到如何使用 VectorTools::point_value 函数在下文中进一步做到这一点。




<a name="AppendixPDEswithDiracdeltafunctionsasrighthandsideandtheirtransformationtoaninitialvalueproblem"></a><h3>Appendix: PDEs with Dirac delta functions as right hand side and their transformation to an initial value problem</h3>


在推导波浪方程的初值问题时，我们最初发现该方程有一个狄拉克三角函数的导数作为右手边。

@f[
\Delta p-\frac{1}{c_0^2} \frac{\partial^2 p}{\partial t^2} = \lambda
a(\mathbf r)\frac{d\delta(t)}{dt}.


@f]

为了看看如何将这个单一的方程转化为具有初始条件的PDE的通常陈述，让我们假设物理上相当合理的介质最初处于静止状态，即 $p(t,\mathbf
r)=\frac{\partial p(t,\mathbf r)}{\partial t}=0$ 为 $t<0$  。接下来，让我们对两边的时间形成不确定的积分。

@f[
\int^t \Delta p\; dt -\int^t \frac{1}{c_0^2} \frac{\partial^2 p}{\partial t^2}
\; dt
=
\int^t \lambda a(\mathbf r)\frac{d\delta(t)}{dt} \;dt.


@f]

这立即引出了一个说法

@f[
P(t,\mathbf r) - \frac{1}{c_0^2} \frac{\partial p}{\partial t}
=
\lambda a(\mathbf r) \delta(t),


@f]

其中 $P(t,\mathbf r)$ 是这样的： $\frac{dP(t,\mathbf r)}{dt}=\Delta
p$  。接下来，我们对 $t=-\epsilon$ 到 $t=+\epsilon$ 的时间进行（定）积分，以求得

@f[
\int_{-\epsilon}^{\epsilon} P(t,\mathbf r)\; dt


- \frac{1}{c_0^2} \left[ p(\epsilon,\mathbf r) - p(-\epsilon,\mathbf r) \right]
=
\int_{-\epsilon}^{\epsilon} \lambda a(\mathbf r) \delta(t) \; dt.


@f]

如果我们利用三角洲函数的属性，即 $\int_{-\epsilon}^{\epsilon}
\delta(t)\; dt = 1$ ，并假设 $P$ 是一个时间上的连续函数，我们发现当我们让 $\epsilon$ 归零时，我们发现

@f[


- \lim_{\epsilon\rightarrow 0}\frac{1}{c_0^2} \left[ p(\epsilon,\mathbf r) - p(-\epsilon,\mathbf r) \right]
=
\lambda a(\mathbf r).


@f]

换句话说，利用 $p(-\epsilon,\mathbf r)=0$ ，我们找回了初始条件

@f[
  \frac{1}{c_0^2} p(0,\mathbf r)
  =
  \lambda a(\mathbf r).


@f]

同时，我们知道，对于每一个 $t>0$ ，三角洲函数都是零，所以对于 $0<t<T$ ，我们得到的方程式是

@f[
\Delta p-\frac{1}{c_0^2} \frac{\partial^2 p}{\partial t^2} = 0.


@f]

因此，我们从原来有些奇怪的方程中得到了一个波浪方程和一个初始条件的表示。

最后，由于我们这里有一个带有两个时间导数的方程，我们仍然需要第二个初始条件。为此，让我们回到方程中去

@f[
\Delta p-\frac{1}{c_0^2} \frac{\partial^2 p}{\partial t^2} = \lambda
a(\mathbf r)\frac{d\delta(t)}{dt}.


@f]

并从 $t=-\epsilon$ 到 $t=+\epsilon$ 进行时间整合。这就导致了

@f[
P(\epsilon)-P(-\epsilon)


-\frac{1}{c_0^2} \left[\frac{\partial p(\epsilon)}{\partial t} -
                       \frac{\partial p(-\epsilon)}{\partial t}\right]
 = \lambda a(\mathbf r) \int_{-\epsilon}^{\epsilon}\frac{d\delta(t)}{dt} \; dt.


@f]

使用部分整合的形式

@f[
  \int_{-\epsilon}^{\epsilon}\varphi(t)\frac{d\delta(t)}{dt} \; dt
  =


  -\int_{-\epsilon}^{\epsilon}\frac{d\varphi(t)}{dt} \delta(t)\; dt


@f]

在这里我们使用 $\delta(\pm \epsilon)=0$ 并插入 $\varphi(t)=1$ ，我们看到事实上

@f[
  \int_{-\epsilon}^{\epsilon}\frac{d\delta(t)}{dt} \; dt
  =
  0.


@f]



现在，让 $\epsilon\rightarrow 0$  。假设 $P$ 是一个时间上的连续函数，我们看到

@f[
  P(\epsilon)-P(-\epsilon) \rightarrow 0,


@f]

因此

@f[
  \frac{\partial p(\epsilon)}{\partial t} -
                       \frac{\partial p(-\epsilon)}{\partial t}
		       \rightarrow 0.


@f]

然而，我们已经假设 $\frac{\partial p(-\epsilon)}{\partial t}=0$  。因此，我们得到的第二个初始条件是

@f[
  \frac{\partial p(0)}{\partial t} = 0,


@f]

完成方程组。


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
 * 以下内容之前都已经介绍过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
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
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 这是唯一一个新的。我们将需要一个定义在GridTools命名空间的库函数，用来计算最小的单元格直径。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step24 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Theforwardproblemclasstemplate"></a> 
 * <h3>The "forward problem" class template</h3>
 * 

 * 
 * 主类的第一部分与 step-23 中的内容完全一致（除了名字）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class TATForwardProblem 
 *   { 
 *   public: 
 *     TATForwardProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void solve_p(); 
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
 *     SparseMatrix<double> system_matrix; 
 *     SparseMatrix<double> mass_matrix; 
 *     SparseMatrix<double> laplace_matrix; 
 * 
 *     Vector<double> solution_p, solution_v; 
 *     Vector<double> old_solution_p, old_solution_v; 
 *     Vector<double> system_rhs_p, system_rhs_v; 
 * 
 *     double       time_step, time; 
 *     unsigned int timestep_number; 
 *     const double theta; 
 * 
 * @endcode
 * 
 * 下面是新的内容：首先，我们需要从吸收边界条件出来的那个边界质量矩阵 $B$ 。同样，由于这次我们考虑的是一个现实的介质，我们必须有一个衡量波速的标准 $c_0$ ，它将进入所有与拉普拉斯矩阵（我们仍然定义为 $(\nabla \phi_i,\nabla \phi_j)$ ）有关的公式。
 * 

 * 
 * 
 * @code
 *     SparseMatrix<double> boundary_matrix; 
 *     const double         wave_speed; 
 * 
 * @endcode
 * 
 * 我们必须注意的最后一件事是，我们想在一定数量的检测器位置评估解决方案。我们需要一个数组来保存这些位置，在这里声明并在构造函数中填充。
 * 

 * 
 * 
 * @code
 *     std::vector<Point<dim>> detector_locations; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 像往常一样，我们必须定义我们的初始值、边界条件和右手边的函数。这次事情有点简单：我们考虑的是一个由初始条件驱动的问题，所以没有右手函数（尽管你可以在 step-23 中查找，看看如何做到这一点）。其次，没有边界条件：域的整个边界由吸收性边界条件组成。这就只剩下初始条件了，这里的事情也很简单，因为对于这个特殊的应用，只规定了压力的非零初始条件，而没有规定速度的非零初始条件（速度在初始时间为零）。
 * 

 * 
 * 所以这就是我们所需要的：一个指定压力初始条件的类。在本程序所考虑的物理环境中，这些是小的吸收器，我们将其建模为一系列的小圆圈，我们假设压力盈余为1，而其他地方没有吸收，因此没有压力盈余。我们是这样做的（注意，如果我们想把这个程序扩展到不仅可以编译，而且可以运行，我们将不得不用三维源的位置来初始化源）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class InitialValuesP : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/ = 0) const override 
 *     { 
 *       static const std::array<Source, 5> sources{ 
 *         {Source(Point<dim>(0, 0), 0.025), 
 *          Source(Point<dim>(-0.135, 0), 0.05), 
 *          Source(Point<dim>(0.17, 0), 0.03), 
 *          Source(Point<dim>(-0.25, 0), 0.02), 
 *          Source(Point<dim>(-0.05, -0.15), 0.015)}}; 
 * 
 *       for (const auto &source : sources) 
 *         if (p.distance(source.location) < source.radius) 
 *           return 1; 
 * 
 *       return 0; 
 *     } 
 * 
 *  
 *     struct Source 
 *     { 
 *       Source(const Point<dim> &l, const double r) 
 *         : location(l) 
 *         , radius(r) 
 *       {} 
 * 
 *       const Point<dim> location; 
 *       const double     radius; 
 *     }; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeTATForwardProblemcodeclass"></a> 
 * <h3>Implementation of the <code>TATForwardProblem</code> class</h3>
 * 

 * 
 * 让我们再从构造函数开始。设置成员变量是很直接的。我们使用矿物油的声波速度（单位为毫米/微秒，是实验性生物医学成像中的常用单位），因为我们想和输出的许多实验都是在这里进行的。再次使用Crank-Nicolson方案，即theta被设定为0.5。随后选择时间步长以满足 $k = \frac hc$ ：这里我们把它初始化为一个无效的数字。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   TATForwardProblem<dim>::TATForwardProblem() 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *     , time_step(std::numeric_limits<double>::quiet_NaN()) 
 *     , time(time_step) 
 *     , timestep_number(1) 
 *     , theta(0.5) 
 *     , wave_speed(1.437) 
 *   { 
 * 
 * @endcode
 * 
 * 构造函数中的第二个任务是初始化存放检测器位置的数组。这个程序的结果与实验进行了比较，其中检测器间距的步长为2.25度，对应160个检测器位置。扫描圆的半径被选为中心和边界之间的一半，以避免不完善的边界条件带来的剩余反射破坏我们的数值结果。
 * 

 * 
 * 然后按顺时针顺序计算探测器的位置。请注意，下面的内容当然只有在我们以2D计算时才有效，我们用一个断言来保护这个条件。如果我们以后想在三维中运行同样的程序，我们将不得不在这里添加代码来初始化三维中的探测器位置。由于断言的存在，我们不可能忘记这样做。
 * 

 * 
 * 
 * @code
 *     Assert(dim == 2, ExcNotImplemented()); 
 * 
 *     const double detector_step_angle = 2.25; 
 *     const double detector_radius     = 0.5; 
 * 
 *     for (double detector_angle = 2 * numbers::PI; detector_angle >= 0; 
 *          detector_angle -= detector_step_angle / 360 * 2 * numbers::PI) 
 *       detector_locations.push_back( 
 *         Point<dim>(std::cos(detector_angle), std::sin(detector_angle)) * 
 *         detector_radius); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TATForwardProblemsetup_system"></a> 
 * <h4>TATForwardProblem::setup_system</h4>
 * 

 * 
 * 下面的系统几乎就是我们在  step-23  中已经做过的，但有两个重要的区别。首先，我们必须在原点周围创建一个半径为1的圆形（或球形）网格。这并不新鲜：我们之前在 step-6 和 step-10 中已经这样做了，在那里我们还解释了PolarManifold或SphericalManifold对象如何在细化单元时将新点放在同心圆上，我们在这里也将使用它。
 * 

 * 
 * 我们必须确保的一点是，时间步长满足  step-23  的介绍中讨论的 CFL 条件。在那个程序中，我们通过设置一个与网格宽度相匹配的时间步长来确保这一点，但是这很容易出错，因为如果我们再细化一次网格，我们也必须确保时间步长有所改变。在这里，我们自动做到了这一点：我们向一个库函数询问任何单元的最小直径。然后我们设置 $k=\frac h{c_0}$  。唯一的问题是： $h$ 到底是什么？关键是，对于波浪方程来说，这个问题确实没有好的理论。众所周知，对于由矩形组成的均匀细化网格， $h$ 是最小边长。但对于一般四边形的网格，确切的关系似乎是未知的，也就是说，不知道单元格的什么属性与CFL条件有关。问题是，CFL条件来自于对拉普拉斯矩阵最小特征值的了解，而这只能对简单结构的网格进行分析计算。
 * 

 * 
 * 这一切的结果是，我们并不十分确定我们应该对 $h$ 采取什么措施。函数 GridTools::minimal_cell_diameter 计算了所有单元的最小直径。如果单元格都是正方形或立方体，那么最小边长就是最小直径除以 <code>std::sqrt(dim)</code>  。我们简单地将此概括为非均匀网格的情况，没有理论上的理由。
 * 

 * 
 * 唯一的其他重大变化是我们需要建立边界质量矩阵。我们将在下文中进一步评论这个问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TATForwardProblem<dim>::setup_system() 
 *   { 
 *     const Point<dim> center; 
 *     GridGenerator::hyper_ball(triangulation, center, 1.); 
 *     triangulation.refine_global(7); 
 * 
 *     time_step = GridTools::minimal_cell_diameter(triangulation) / wave_speed / 
 *                 std::sqrt(1. * dim); 
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
 *     system_matrix.reinit(sparsity_pattern); 
 *     mass_matrix.reinit(sparsity_pattern); 
 *     laplace_matrix.reinit(sparsity_pattern); 
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
 * 如前所述，与 step-23 的第二个区别是，我们需要建立从吸收性边界条件中生长出来的边界质量矩阵。
 * 

 * 
 * 第一个观察结果是，这个矩阵比常规质量矩阵要稀疏得多，因为没有一个具有纯内部支持的形状函数对这个矩阵有贡献。因此，我们可以根据这种情况优化存储模式，建立第二个稀疏模式，只包含我们需要的非零项。这里有一个权衡：首先，我们必须要有第二个稀疏模式对象，所以这需要花费内存。其次，与该稀疏性模式相连的矩阵将更小，因此需要更少的内存；用它进行矩阵-向量乘法也会更快。然而，最后一个论点是提示规模的论点：我们主要感兴趣的不是单独对边界矩阵进行矩阵-向量运算（尽管我们需要在每个时间步长对右侧向量进行一次运算），而是主要希望将其与两个方程中的第一个方程使用的其他矩阵相加，因为这是CG方法每个迭代都要与之相乘的，即明显更频繁。现在的情况是， SparseMatrix::add 类允许将一个矩阵添加到另一个矩阵中，但前提是它们使用相同的稀疏模式（原因是我们不能在稀疏模式创建后向矩阵添加非零条目，所以我们只是要求这两个矩阵具有相同的稀疏模式）。
 * 

 * 
 * 所以，我们就用这个方法吧。
 * 

 * 
 * 
 * @code
 *     boundary_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 第二件要做的事是实际建立矩阵。在这里，我们需要对单元格的面进行积分，所以首先我们需要一个能在 <code>dim-1</code> 维对象上工作的正交对象。其次，FEValues的变体FEFaceValues，正如它的名字所暗示的，它可以在面上工作。最后，其他的变量是组装机器的一部分。所有这些我们都放在大括号里，以便将这些变量的范围限制在我们真正需要它们的地方。
 * 然后
 * 组装矩阵的实际行为是相当直接的：我们在所有单元中循环，在每个单元的所有面中循环，然后只在特定的面位于域的边界时做一些事情。像这样。
 * 

 * 
 * 
 * @code
 *     { 
 *       const QGauss<dim - 1> quadrature_formula(fe.degree + 1); 
 *       FEFaceValues<dim>     fe_values(fe, 
 *                                   quadrature_formula, 
 *                                   update_values | update_JxW_values); 
 * 
 *       const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *       const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *       FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *       std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary()) 
 *             { 
 *               cell_matrix = 0; 
 * 
 *               fe_values.reinit(cell, face); 
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     cell_matrix(i, j) += (fe_values.shape_value(i, q_point) * 
 *                                           fe_values.shape_value(j, q_point) * 
 *                                           fe_values.JxW(q_point)); 
 * 
 *               cell->get_dof_indices(local_dof_indices); 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                   boundary_matrix.add(local_dof_indices[i], 
 *                                       local_dof_indices[j], 
 *                                       cell_matrix(i, j)); 
 *             } 
 *     } 
 * 
 *     system_matrix.copy_from(mass_matrix); 
 *     system_matrix.add(time_step * time_step * theta * theta * wave_speed * 
 *                         wave_speed, 
 *                       laplace_matrix); 
 *     system_matrix.add(wave_speed * theta * time_step, boundary_matrix); 
 * 
 *     solution_p.reinit(dof_handler.n_dofs()); 
 *     old_solution_p.reinit(dof_handler.n_dofs()); 
 *     system_rhs_p.reinit(dof_handler.n_dofs()); 
 * 
 *     solution_v.reinit(dof_handler.n_dofs()); 
 *     old_solution_v.reinit(dof_handler.n_dofs()); 
 *     system_rhs_v.reinit(dof_handler.n_dofs()); 
 * 
 *     constraints.close(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TATForwardProblemsolve_pandTATForwardProblemsolve_v"></a> 
 * <h4>TATForwardProblem::solve_p and TATForwardProblem::solve_v</h4>
 * 

 * 
 * 下面两个函数，解决压力和速度变量的线性系统，几乎是逐字逐句地从 step-23 中提取的（除了主变量的名字从 $u$ 改为 $p$ ）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TATForwardProblem<dim>::solve_p() 
 *   { 
 *     SolverControl solver_control(1000, 1e-8 * system_rhs_p.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     cg.solve(system_matrix, solution_p, system_rhs_p, PreconditionIdentity()); 
 * 
 *     std::cout << "   p-equation: " << solver_control.last_step() 
 *               << " CG iterations." << std::endl; 
 *   } 
 * 
 *   template <int dim> 
 *   void TATForwardProblem<dim>::solve_v() 
 *   { 
 *     SolverControl solver_control(1000, 1e-8 * system_rhs_v.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     cg.solve(mass_matrix, solution_v, system_rhs_v, PreconditionIdentity()); 
 * 
 *     std::cout << "   v-equation: " << solver_control.last_step() 
 *               << " CG iterations." << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TATForwardProblemoutput_results"></a> 
 * <h4>TATForwardProblem::output_results</h4>
 * 

 * 
 * 这里也是如此：该函数来自  step-23  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TATForwardProblem<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution_p, "P"); 
 *     data_out.add_data_vector(solution_v, "V"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
 *     DataOutBase::VtkFlags vtk_flags; 
 *     vtk_flags.compression_level = 
 *       DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TATForwardProblemrun"></a> 
 * <h4>TATForwardProblem::run</h4>
 * 

 * 
 * 这个做大部分工作的函数又和 step-23 中的差不多，尽管我们通过使用介绍中提到的向量G1和G2使事情变得更加清晰。与程序的整体内存消耗相比，引入几个临时向量并没有什么坏处。
 * 

 * 
 * 这个函数唯一的变化是：首先，我们不必为速度 $v$ 预测初始值，因为我们知道它是零。其次，我们在构造函数中计算的检测器位置上评估解决方案。这是用 VectorTools::point_value 函数完成的。然后，这些值被写入我们在函数开始时打开的一个文件中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TATForwardProblem<dim>::run() 
 *   { 
 *     setup_system(); 
 * 
 *     VectorTools::project(dof_handler, 
 *                          constraints, 
 *                          QGauss<dim>(fe.degree + 1), 
 *                          InitialValuesP<dim>(), 
 *                          old_solution_p); 
 *     old_solution_v = 0; 
 * 
 *     std::ofstream detector_data("detectors.dat"); 
 * 
 *     Vector<double> tmp(solution_p.size()); 
 *     Vector<double> G1(solution_p.size()); 
 *     Vector<double> G2(solution_v.size()); 
 * 
 *     const double end_time = 0.7; 
 *     for (time = time_step; time <= end_time; 
 *          time += time_step, ++timestep_number) 
 *       { 
 *         std::cout << std::endl; 
 *         std::cout << "time_step " << timestep_number << " @ t=" << time 
 *                   << std::endl; 
 * 
 *         mass_matrix.vmult(G1, old_solution_p); 
 *         mass_matrix.vmult(tmp, old_solution_v); 
 *         G1.add(time_step * (1 - theta), tmp); 
 * 
 *         mass_matrix.vmult(G2, old_solution_v); 
 *         laplace_matrix.vmult(tmp, old_solution_p); 
 *         G2.add(-wave_speed * wave_speed * time_step * (1 - theta), tmp); 
 * 
 *         boundary_matrix.vmult(tmp, old_solution_p); 
 *         G2.add(wave_speed, tmp); 
 * 
 *         system_rhs_p = G1; 
 *         system_rhs_p.add(time_step * theta, G2); 
 * 
 *         solve_p(); 
 * 
 *         system_rhs_v = G2; 
 *         laplace_matrix.vmult(tmp, solution_p); 
 *         system_rhs_v.add(-time_step * theta * wave_speed * wave_speed, tmp); 
 * 
 *         boundary_matrix.vmult(tmp, solution_p); 
 *         system_rhs_v.add(-wave_speed, tmp); 
 * 
 *         solve_v(); 
 * 
 *         output_results(); 
 * 
 *         detector_data << time; 
 *         for (unsigned int i = 0; i < detector_locations.size(); ++i) 
 *           detector_data << " " 
 *                         << VectorTools::point_value(dof_handler, 
 *                                                     solution_p, 
 *                                                     detector_locations[i]) 
 *                         << " "; 
 *         detector_data << std::endl; 
 * 
 *         old_solution_p = solution_p; 
 *         old_solution_v = solution_v; 
 *       } 
 *   } 
 * } // namespace Step24 
 * 
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
 *       using namespace Step24; 
 * 
 *       TATForwardProblem<2> forward_problem_solver; 
 *       forward_problem_solver.run(); 
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
examples/step-24/doc/results.dox



<a name="Results"></a><h1>Results</h1>


该程序将每个时间步骤的图形数据以及每个探测器位置的评估值都写入磁盘。然后我们将它们绘制成图。还收集了实验数据进行比较。目前，我们的实验只在二维空间通过圆形扫描单个探测器进行。这里的组织样本是 $X-Y$ 平面的薄片（ $Z=0$ ），我们假设其他 $Z$ 方向的信号不会对数据产生影响。因此，我们只需要将我们的实验数据与二维模拟数据进行比较。

<a name="Oneabsorber"></a><h3> One absorber </h3>


这部电影显示了由单个小吸收器产生的热声波在介质中传播（在我们的模拟中，我们假设介质是矿物油，其声速为1.437  $\frac{mm}{\mu s}$  ）。

 <img src="https://www.dealii.org/images/steps/developer/step-24.one_movie.gif" alt=""> 

对于单个吸收器，我们当然要相应地改变 <code>InitialValuesP</code> 类。

接下来，让我们比较一下实验和计算的结果。可视化使用了一种在地震学中长期使用的技术，即把每个探测器的数据全部绘制在一张图上。这样做的方法是将每个探测器的信号与前一个探测器相比偏移一点。例如，这里是前四个探测器的图（从下到上，时间从左到右为微秒），使用程序中使用的源设置，与目前只有一个源的情况相比，使事情更有趣。

 <img src="https://www.dealii.org/images/steps/developer/step-24.traces.png" alt=""> 

例如，可以看到的一点是，第二和第四个信号的到达时间在探测器数量较多的情况下（即最上面的探测器）会转移到较早的时间，但第一和第三信号则不然；这可以解释为这些信号的起源必须比前者更接近后一个探测器。

如果我们不仅将4个，而是将所有160个探测器堆叠在一张图中，单个线条就会变得模糊，但在它们运行在一起的地方，就会形成一种较深或较浅的灰度模式。  下面两张图显示了在以这种方式堆叠的探测器位置获得的结果。左图是由实验得到的，右图是模拟数据。在实验中，一个小的强吸收器被嵌入到较弱的吸收组织中。

 <table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.one.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.one_s.png" alt="">
</td>
</tr>
</table> 

很明显，在角度 $180^\circ$ 处，源的位置离探测器更近。在实验数据中可以看到的所有其他信号都是由于组织的其他部分也有弱的吸收体，这些吸收体环绕着中心的小强吸收体产生的信号。另一方面，在模拟数据中，我们只模拟了小的强吸收体。

在现实中，探测器的带宽有限。因此，通过探测器的热声波将被过滤掉。通过使用高通滤波器（在MATLAB中实现并针对本程序产生的数据文件运行），可以使模拟结果看起来更接近于实验数据。

 <img src="https://www.dealii.org/images/steps/developer/step-24.one_sf.png" alt=""> 

在我们的模拟中，我们看到主波后面的假信号是由数值伪影造成的。这个问题可以通过使用更细的网格来缓解，从而得到下面的图。

 <img src="https://www.dealii.org/images/steps/developer/step-24.one_s2.png" alt=""> 




<a name="Multipleabsorbers"></a><h3>Multiple absorbers</h3>


为了进一步验证该程序，我们还将展示多个吸收器的模拟结果。这与程序中实际实现的情况相对应。下面的影片显示了由多个吸收器产生的热声波在介质中的传播情况。

 <img src="https://www.dealii.org/images/steps/developer/step-24.multi_movie.gif" alt=""> 

实验数据和我们的模拟数据在以下两个图中进行了比较。   <table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.multi.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.multi_s.png" alt="">
</td>
</tr>
</table> 

请注意，在实验数据中，第一个信号（即最左边的暗线）来自于组织边界的吸收，因此首先到达检测器，比来自内部的任何信号都要早。这个信号在痕迹的末端也是微弱可见的，大约在30 $\mu s$ ，这表明信号穿过整个组织到达另一侧的探测器，在所有来自内部的信号到达它们之后。

和以前一样，通过应用符合探测器实际行为的带宽滤波器（左）和选择更细的网格（右），数值结果与实验结果更加匹配。

 <table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.multi_sf.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-24.multi_s2.png" alt="">
</td>
</tr>
</table> 

左图和右图的一个重要区别是，右图的曲线看起来没有那么多 "棱角"。角度来自于这样一个事实：虽然连续方程中的波在各个方向上的移动速度相同，但离散化后的情况并非如此：在那里，对角线上的波与平行于网格线的波的移动速度略有不同。这种各向异性导致波前不是完全的圆形（在堆积图中会产生正弦信号），而是在某些方向上凸出。更糟糕的是，我们使用的圆形网格（例如，见步骤6的粗略网格图）也不是各向同性的。最终的结果是，除非网格足够细，否则信号锋面不是正弦波的。右图在这方面要好得多，尽管仍然可以看到拖尾假波形式的伪影。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-24.cc"
*/
