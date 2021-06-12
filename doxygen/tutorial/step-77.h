/**
@page step_77 The step-77 tutorial program
This tutorial depends on step-15.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#HowdealIIinterfaceswithKINSOL"> How deal.II interfaces with KINSOL </a>
        <li><a href="#Detailsoftheimplementation"> Details of the implementation </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclasstemplate">The <code>MinimalSurfaceProblem</code> class template</a>
        <li><a href="#Boundarycondition">Boundary condition</a>
        <li><a href="#ThecodeMinimalSurfaceProblemcodeclassimplementation">The <code>MinimalSurfaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#Constructorandsetupfunctions">Constructor and set up functions</a>
        <li><a href="#AssemblingandfactorizingtheJacobianmatrix">Assembling and factorizing the Jacobian matrix</a>
        <li><a href="#Computingtheresidualvector">Computing the residual vector</a>
        <li><a href="#SolvinglinearsystemswiththeJacobianmatrix">Solving linear systems with the Jacobian matrix</a>
        <li><a href="#Refiningthemeshsettingboundaryvaluesandgeneratinggraphicaloutput">Refining the mesh, setting boundary values, and generating graphical output</a>
        <li><a href="#Therunfunctionandtheoveralllogicoftheprogram">The run() function and the overall logic of the program</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-77/doc/intro.dox

 <br> 

<i>
This program was contributed by Wolfgang Bangerth, Colorado State University.


This material is based upon work partially supported by National Science
Foundation grants OAC-1835673, DMS-1821210, and EAR-1925595;
and by the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-1550901 and The University of California-Davis.
</i> <br>  。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


第15步程序解决了以下描述最小表面问题的非线性方程。

@f{align*}{


    -\nabla \cdot \left( \frac{1}{\sqrt{1+|\nabla u|^{2}}}\nabla u \right) &= 0 \qquad
    \qquad &&\textrm{in} ~ \Omega
    \\
    u&=g \qquad\qquad &&\textrm{on} ~ \partial \Omega.


@f}

step-15使用的是牛顿方法，牛顿方法的工作原理是反复求解一个更新 $\delta u_k$ 的*线性化*问题--称为 "搜索方向"--计算 "步长" $\alpha_k$ ，然后将它们结合起来，通过以下方式计算出新的猜测解。

@f{align*}{
    u_{k+1} = u_k + \alpha_k \, \delta u_k.


@f}



在步骤15的讨论过程中，我们发现计算步长是很尴尬的，所以只是解决了简单的选择。总是选择 $\alpha_k=0.1$  。这当然是没有效率的。我们知道，只有当我们最终能够选择 $\alpha_k=1$ 时，我们才能实现牛顿的二次收敛率，尽管在最初的几次迭代中，我们可能不得不选择较小的步长，因为我们离使用这么长的步长还很遥远。

因此，本方案的目标之一是解决这一缺陷。由于行搜索算法的实现并不完全是微不足道的，因此人们无论如何都要做自己应该做的事。从外部库中导入复杂的功能。为此，我们将利用deal.II与大型非线性求解器包之一的接口，即[SUNDIALS](https://computing.llnl.gov/projects/sundials)套件的[KINSOL](https://computing.llnl.gov/projects/sundials/kinsol)子包。%SUNDIALS的核心是一个用于解决复杂的常微分方程（ODE）和微分代数方程（DAE）的软件包，deal.II接口允许通过SUNDIALS命名空间的类来实现。特别是 SUNDIALS::ARKode 和 SUNDIALS::IDA 类。但是，由于这是用隐式方法解决ODE和DAE的一个重要步骤，%SUNDIALS也有一个非线性问题的求解器，叫做KINSOL，deal.II有一个接口，以 SUNDIALS::KINSOL 类的形式与之连接。这就是我们将用于解决我们的问题的方法。

但是%SUNDIALS不仅仅是一个方便我们避免编写线搜索算法的方法。一般来说，非线性问题的解决是相当昂贵的，人们通常希望尽可能地节省计算时间。一个可以实现这个目标的方法是如下的。第15步中的算法将问题离散化，然后在每次迭代中求解形式为的线性系统

@f{align*}{
  J_k \, \delta U_k = -F_k


@f}

其中 $F_k$ 是使用当前节点值矢量 $U_k$ 计算的残差矢量， $J_k$ 是其导数（称为 "雅各布"），而 $\delta U_k$ 是对应于上述函数 $\delta u_k$ 的更新矢量。步骤15中已经彻底讨论了 $J_k,F_k$ 的构造，以及在每个牛顿迭代中解决线性系统的方法。因此，让我们关注一下非线性求解过程的另一个方面。计算 $F_k$ 是昂贵的，而组装矩阵 $J_k$ 更是如此。我们真的需要在每次迭代中都这样做吗？事实证明，在许多应用中，这实际上是没有必要的。即使我们用近似值 $\tilde J_k$ 代替 $J_k$ ，这些方法通常也能收敛，并解决了

@f{align*}{
  \tilde J_k \, \widetilde{\delta U}_k = -F_k


@f}

代替，然后更新

@f{align*}{
    U_{k+1} = U_k + \alpha_k \, \widetilde{\delta U}_k.


@f}

这可能需要多一两个迭代，因为我们的更新 $\widetilde{\delta U}_k$ 并不像 $\delta U_k$ 那样好，但它可能仍然是一个胜利，因为我们不必经常组装 $J_k$ 。

对于 $J_k$ ，我们希望得到什么样的近似值 $\tilde J_k$ ？理论上说，由于 $U_k$ 收敛于精确解 $U^\ast$ ，我们需要确保 $\tilde J_k$ 需要收敛于 $J^\ast = \nabla F(U^\ast)$  。特别是，由于  $J_k\rightarrow J^\ast$  ，有效的选择是  $\tilde J_k = J_k$  。但是每一次，比如说，第五次迭代选择 $\tilde J_k = J_k$ 也是如此，对于其他的迭代，我们选择 $\tilde J_k$ 等于最后计算的 $J_{k'}$  。这就是我们在这里要做的：我们将只是重新使用前一次迭代中的 $\tilde J_{k-1}$ ，这可能又是我们在之前的迭代中使用的， $\tilde J_{k-2}$  。

如果对于带有 $J_k$ 的线性系统的求解，我们不只是要组装一个矩阵，还要计算一个好的预处理程序，那么这个方案就变得更加有趣。例如，如果我们要通过SparseDirectUMFPACK类使用稀疏LU分解，或者使用几何或代数多重网格。在这些情况下，我们也不必更新预处理程序，因为预处理程序的计算时间可能和当初组装矩阵的时间一样长，甚至更长。事实上，在这种心态下，我们也许应该考虑使用我们能想到的*好的前置条件器，尽管它们的构造通常相当昂贵。我们希望通过将其应用于不止一个线性求解，来摊销计算这个预处理程序的成本。

当然，最大的问题是。我们根据什么标准来决定我们是否可以摆脱基于先前计算的雅各布矩阵 $J_{k-s}$ 的近似 $s$ 步，或者我们是否需要--至少在这个迭代中--实际重新计算雅各布 $J_k$ 和相应的前置条件器？这就像行搜索的问题一样，需要大量的代码来监控整个算法的收敛性。我们*可以*自己实现这些东西，但我们可能*不应该*。KINSOL已经为我们做了这些。它将告诉我们的代码何时要 "更新 "雅各布矩阵。

如果我们要使用迭代求解器而不是上面提到的稀疏直接求解器，还有最后一个考虑。在求解更新 $\delta U_k$ 时，不仅有可能用一些近似值 $J_k$ 代替 $\tilde J_k$ ，而且还可以问是否有必要求解线性系统

@f{align*}{
  \tilde J_k \widetilde{\delta U}_k = -F_k


@f}

准确度高。其思路是这样的。虽然我们目前的解决方案 $U_k$ 离 $U^\ast$ 还很远，但我们为什么要特别精确地解决这个线性系统？更新后的 $U_{k+1}=U_k + \widetilde{\delta U}_k$ 很可能仍然离精确的解决方案很远，那么为什么要花很多时间来解决这个线性系统的精确性？这就是 "Eisenstat-Walker技巧" @cite eiwa96 等算法的基础思维，在该算法中，人们被赋予一个公差，在迭代 $k$ 中必须解决上述线性系统，该公差取决于整个非线性求解器的进展。像以前一样，我们可以尝试自己实现，但是KINSOL已经为我们提供了这种信息--尽管我们不会在这个程序中使用它，因为我们使用的是直接求解器，不需要求解器的容忍度，只是精确求解线性系统到舍入。

作为对所有这些考虑的总结，我们可以说以下几点。没有必要重新发明轮子。就像deal.II提供了大量的有限元功能一样，%SUNDIALS的KINSOL软件包提供了大量的非线性求解器功能，我们最好使用它。




<a name="HowdealIIinterfaceswithKINSOL"></a><h3> How deal.II interfaces with KINSOL </h3>


KINSOL，像许多类似的软件包一样，以一种相当抽象的方式工作。在其核心部分，它看到了一个非线性问题，其形式为

@f{align*}{
    F(U) = 0


@f}

并构建一个迭代序列  $U_k$  ，一般来说，迭代序列是与函数返回的向量相同长度的向量  $F$  。要做到这一点，它需要从用户那里得到一些东西。

- 将一个给定的向量调整到正确大小的方法。

- 对于一个给定的向量 $U$ ，评估函数 $F(U)$ 的一种方法。这个函数通常被称为 "剩余 "操作，因为目标当然是找到一个点 $U^\ast$ ，对于这个点 $F(U^\ast)=0$ ；如果 $F(U)$ 返回一个非零向量，那么这就是<a href="https://en.wikipedia.org/wiki/Residual_(numerical_analysis)">"residual"</a>（即 "剩余"，或任何 "剩余"）。做到这一点的函数在本质上与步骤15中的右手边向量的计算相同，但有一个重要区别。   在那里，右手边表示的是残差的*负数，所以我们必须换一个符号。

- 计算矩阵 $J_k$ 的方法，如果这在当前迭代中是必要的，同时可能还有一个预处理程序或其他数据结构（例如，通过SparseDirectUMFPACK进行稀疏分解，如果我们选择用它来解决一个线性系统）。这个操作通常被称为 "设置 "操作。

- 用最后计算的任何矩阵 $\tilde J_k$ 来解决一个线性系统 $\tilde J_k x = b$ 的方法。这个操作一般被称为 "求解 "操作。

所有这些操作都需要由 [std::function](https://en.cppreference.com/w/cpp/utility/functional/function) 对象提供给KINSOL，这些对象接受适当的参数集，通常返回一个表示成功（返回值为零）或失败（返回值为非零）的整数。具体来说，我们要访问的对象是 SUNDIALS::KINSOL::reinit_vector,   SUNDIALS::KINSOL::residual,   SUNDIALS::KINSOL::setup_jacobian,  和  SUNDIALS::KINSOL::solve_jacobian_system  成员变量。(详见这些变量的文档。)在我们的实现中，我们将使用[lambda functions](https://en.cppreference.com/w/cpp/language/lambda)来实现这些 "回调"，反过来可以调用成员函数；然后KINSOL将在其内部算法认为有用时调用这些回调。




<a name="Detailsoftheimplementation"></a><h3> Details of the implementation </h3>


本教程程序的大部分代码与步骤15一样，我们将不作过多评论。实际上只有一个方面需要注意，即一方面给定一个向量 $U$ ，另一方面给定一个向量 $J(U)$ ，如何计算 $U$ 。起初，这似乎很简单：我们只需使用`assemble_system()`函数，在一种情况下抛出所有处理矩阵的代码，在另一种情况下抛出右手边的向量。就这样。问题解决了。

但它并不那么简单。这是因为如果我们有非零的Dirichlet边界值，这两者并不独立，就像我们在这里做的那样。我们要解决的线性系统包含内部和边界自由度，当从那些真正 "自由 "的自由度中消除这些自由度时，使用例如 AffineConstraints::distribute_local_to_global(), ，我们在组装右手边的向量时需要知道矩阵。

当然，这完全违背了原意。如果我们可以不组装矩阵，就不要*组装。我们解决这个问题的方法如下。

- 我们将解向量的起始猜测， $U_0$ ，设定为边界自由度已经有了正确的值。

- 这意味着所有的更新都可以有这些自由度的零更新，我们可以建立残差向量 $F(U_k)$ 和雅各布矩阵 $J_k$ ，对应于线性系统的解在这些向量分量中为零。对于这种特殊情况，矩阵和右手边向量的组装是独立的，可以分解成不同的函数。

这里有一个假设，即每当KINSOL要求用雅各布的（近似值）进行线性求解时，这将是为了更新 $\delta U$ （其边界值为零），其倍数将被添加到解决方案（其已经有正确的边界值）。  这可能不是真的，如果是的话，我们可能要重新考虑我们的方法。也就是说，事实证明，在实践中，这正是KINSOL在使用牛顿方法时的表现，因此我们的方法是成功的。


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
 * 这个程序开始时和其他大多数程序一样，有众所周知的包含文件。与 step-15 程序相比，我们在这里所做的大部分工作都是从该程序中复制的，唯一不同的是包括头文件，我们从该文件中导入了SparseDirectUMFPACK类和KINSOL的实际接口。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_accessor.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_q.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * #include <deal.II/sundials/kinsol.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step77 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class template</h3>
 * 

 * 
 * 同样地，这个程序的主类基本上是  step-15  中的一个副本。然而，该类确实将雅各布（系统）矩阵（以及使用直接求解器对其进行因式分解）和残差的计算分成了不同的函数，原因已在介绍中列出。出于同样的原因，该类也有一个指向雅各布矩阵因式分解的指针，该指针在我们每次更新雅各布矩阵时被重置。
 * 

 * 
 * （如果你想知道为什么程序对雅各布矩阵使用直接对象，而对因式分解使用指针。每次KINSOL要求更新雅各布矩阵时，我们可以简单地写`jacobian_matrix=0;`将其重置为一个空矩阵，然后我们可以再次填充。另一方面，SparseDirectUMFPACK类没有办法扔掉它的内容或用新的因式分解来替换它，所以我们使用一个指针。我们只是扔掉整个对象，并在我们有新的雅各布矩阵需要分解时创建一个新的对象。)
 * 

 * 
 * 最后，该类有一个定时器变量，我们将用它来评估程序的不同部分需要多长时间，这样我们就可以评估KINSOL的不重建矩阵及其因式分解的倾向是否合理。我们将在下面的 "结果 "部分讨论这个问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class MinimalSurfaceProblem 
 *   { 
 *   public: 
 *     MinimalSurfaceProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(const bool initial_step); 
 *     void solve(const Vector<double> &rhs, 
 *                Vector<double> &      solution, 
 *                const double          tolerance); 
 *     void refine_mesh(); 
 *     void output_results(const unsigned int refinement_cycle); 
 *     void set_boundary_values(); 
 *     void compute_and_factorize_jacobian(const Vector<double> &evaluation_point); 
 *     void compute_residual(const Vector<double> &evaluation_point, 
 *                           Vector<double> &      residual); 
 * 
 *     Triangulation<dim> triangulation; 
 * 
 *     DoFHandler<dim> dof_handler; 
 *     FE_Q<dim>       fe; 
 * 
 *     AffineConstraints<double> hanging_node_constraints; 
 * 
 *     SparsityPattern                      sparsity_pattern; 
 *     SparseMatrix<double>                 jacobian_matrix; 
 *     std::unique_ptr<SparseDirectUMFPACK> jacobian_matrix_factorization; 
 * 
 *     Vector<double> current_solution; 
 * 
 *     TimerOutput computing_timer; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="Boundarycondition"></a> 
 * <h3>Boundary condition</h3>
 * 

 * 
 * 实现边界值的类是对  step-15  的复制。
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
 *   double BoundaryValues<dim>::value(const Point<dim> &p, 
 *                                     const unsigned int /*component*/) const 
 *   { 
 *     return std::sin(2 * numbers::PI * (p[0] + p[1])); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ThecodeMinimalSurfaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>MinimalSurfaceProblem</code> class implementation</h3>
 * 
 * <a name="Constructorandsetupfunctions"></a> 
 * <h4>Constructor and set up functions</h4>
 * 

 * 
 * 下面的几个函数也基本上是复制了 step-15 已经做的事情，所以没有什么可讨论的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   MinimalSurfaceProblem<dim>::MinimalSurfaceProblem() 
 *     : dof_handler(triangulation) 
 *     , fe(1) 
 *     , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times) 
 *   {} 
 * 
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "set up"); 
 * 
 *     if (initial_step) 
 *       { 
 *         dof_handler.distribute_dofs(fe); 
 *         current_solution.reinit(dof_handler.n_dofs()); 
 * 
 *         hanging_node_constraints.clear(); 
 *         DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                                 hanging_node_constraints); 
 *         hanging_node_constraints.close(); 
 *       } 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 * 
 *     hanging_node_constraints.condense(dsp); 
 * 
 *     sparsity_pattern.copy_from(dsp); 
 *     jacobian_matrix.reinit(sparsity_pattern); 
 *     jacobian_matrix_factorization.reset(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="AssemblingandfactorizingtheJacobianmatrix"></a> 
 * <h4>Assembling and factorizing the Jacobian matrix</h4>
 * 

 * 
 * 然后，下面的函数负责对雅各布矩阵进行组装和因子化。该函数的前半部分实质上是 step-15 的`assemble_system()`函数，只是它没有处理同时形成右手边的向量（即残差），因为我们并不总是要同时做这些操作。
 * 

 * 
 * 我们把整个装配功能放在一个由大括号包围的代码块中，这样我们就可以用一个 TimerOutput::Scope 变量来衡量在这个代码块中花费了多少时间，不包括在这个函数中发生在匹配的闭合括号`}`之后的一切。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::compute_and_factorize_jacobian( 
 *     const Vector<double> &evaluation_point) 
 *   { 
 *     { 
 *       TimerOutput::Scope t(computing_timer, "assembling the Jacobian"); 
 * 
 *       std::cout << "  Computing Jacobian matrix" << std::endl; 
 * 
 *       const QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *       jacobian_matrix = 0; 
 * 
 *       FEValues<dim> fe_values(fe, 
 *                               quadrature_formula, 
 *                               update_gradients | update_quadrature_points | 
 *                                 update_JxW_values); 
 * 
 *       const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *       const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *       FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *       std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points); 
 * 
 *       std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         { 
 *           cell_matrix = 0; 
 * 
 *           fe_values.reinit(cell); 
 * 
 *           fe_values.get_function_gradients(evaluation_point, 
 *                                            evaluation_point_gradients); 
 * 
 *           for (unsigned int q = 0; q < n_q_points; ++q) 
 *             { 
 *               const double coeff = 
 *                 1.0 / std::sqrt(1 + evaluation_point_gradients[q] * 
 *                                       evaluation_point_gradients[q]); 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     cell_matrix(i, j) += 
 *                       (((fe_values.shape_grad(i, q)    // ((\nabla \phi_i 
 *                          * coeff                       //   * a_n 
 *                          * fe_values.shape_grad(j, q)) //   * \nabla \phi_j) 
 *                         -                              //  - 
 *                         (fe_values.shape_grad(i, q)    //  (\nabla \phi_i 
 *                          * coeff * coeff * coeff       //   * a_n^3 
 *                          * 
 *                          (fe_values.shape_grad(j, q)       //   * (\nabla \phi_j 
 *                           * evaluation_point_gradients[q]) //      * \nabla u_n) 
 *                          * evaluation_point_gradients[q])) //   * \nabla u_n))) 
 *                        * fe_values.JxW(q));                // * dx 
 *                 } 
 *             } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix, 
 *                                                               local_dof_indices, 
 *                                                               jacobian_matrix); 
 *         } 
 * 
 *       std::map<types::global_dof_index, double> boundary_values; 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                Functions::ZeroFunction<dim>(), 
 *                                                boundary_values); 
 *       Vector<double> dummy_solution(dof_handler.n_dofs()); 
 *       Vector<double> dummy_rhs(dof_handler.n_dofs()); 
 *       MatrixTools::apply_boundary_values(boundary_values, 
 *                                          jacobian_matrix, 
 *                                          dummy_solution, 
 *                                          dummy_rhs); 
 *     } 
 * 
 * @endcode
 * 
 * 该函数的后半部分是对计算出的矩阵进行因数分解。为此，我们首先创建一个新的SparseDirectUMFPACK对象，并将其分配给成员变量`jacobian_matrix_factorization`，同时销毁该指针之前指向的任何对象（如果有）。然后我们告诉该对象对雅各布系数进行分解。
 * 

 * 
 * 如上所述，我们把这段代码放在大括号里，用一个计时器来评估这部分程序所需的时间。
 * 

 * 
 * (严格来说，我们在这里完成后实际上不再需要矩阵了，我们可以把矩阵对象扔掉。一个旨在提高内存效率的代码会这样做，并且只在这个函数中创建矩阵对象，而不是作为周围类的成员变量。我们在这里省略了这一步，因为使用与以前的教程程序相同的编码风格可以培养对通用风格的熟悉，并有助于使这些教程程序更容易阅读)。
 * 

 * 
 * 
 * @code
 *     { 
 *       TimerOutput::Scope t(computing_timer, "factorizing the Jacobian"); 
 * 
 *       std::cout << "  Factorizing Jacobian matrix" << std::endl; 
 * 
 *       jacobian_matrix_factorization = std::make_unique<SparseDirectUMFPACK>(); 
 *       jacobian_matrix_factorization->factorize(jacobian_matrix); 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingtheresidualvector"></a> 
 * <h4>Computing the residual vector</h4>
 * 

 * 
 * `assemble_system()`在 step-15 中用来做的第二部分是计算残差向量，也就是牛顿线性系统的右手向量。我们把这一点从前面的函数中分解出来，但如果你理解了 step-15 中`assemble_system()`的作用，下面的函数就会很容易理解。然而，重要的是，我们需要计算的残差不是围绕当前解向量线性化的，而是我们从KINSOL得到的任何东西。这对于诸如直线搜索这样的操作是必要的，我们想知道在不同的 $\alpha_k$ 值下，残差 $F(U^k + \alpha_k \delta U^K)$ 是多少；在这些情况下，KINSOL只是给我们函数 $F$ 的参数，然后我们在这时计算残差 $F(\cdot)$ 。
 * 

 * 
 * 该函数在最后打印出如此计算的残差的规范，作为我们跟踪程序进展的一种方式。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::compute_residual( 
 *     const Vector<double> &evaluation_point, 
 *     Vector<double> &      residual) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "assembling the residual"); 
 * 
 *     std::cout << "  Computing residual vector..." << std::flush; 
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1); 
 *     FEValues<dim>     fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_gradients | update_quadrature_points | 
 *                               update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     Vector<double>              cell_residual(dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_residual = 0; 
 *         fe_values.reinit(cell); 
 * 
 *         fe_values.get_function_gradients(evaluation_point, 
 *                                          evaluation_point_gradients); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             const double coeff = 
 *               1.0 / std::sqrt(1 + evaluation_point_gradients[q] * 
 *                                     evaluation_point_gradients[q]); 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               cell_residual(i) = (fe_values.shape_grad(i, q) // \nabla \phi_i 
 *                                   * coeff                    // * a_n 
 *                                   * evaluation_point_gradients[q] // * u_n 
 *                                   * fe_values.JxW(q));            // * dx 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           residual(local_dof_indices[i]) += cell_residual(i); 
 *       } 
 * 
 *     hanging_node_constraints.condense(residual); 
 * 
 *     for (const types::global_dof_index i : 
 *          DoFTools::extract_boundary_dofs(dof_handler)) 
 *       residual(i) = 0; 
 * 
 *     for (const types::global_dof_index i : 
 *          DoFTools::extract_hanging_node_dofs(dof_handler)) 
 *       residual(i) = 0; 
 * 
 *     std::cout << " norm=" << residual.l2_norm() << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="SolvinglinearsystemswiththeJacobianmatrix"></a> 
 * <h4>Solving linear systems with the Jacobian matrix</h4>
 * 

 * 
 * 接下来是实现用雅各布矩阵解线性系统的函数。由于我们在建立矩阵时已经对矩阵进行了因式分解，所以解决线性系统的方法就是将逆矩阵应用于给定的右侧向量。这就是我们在这里使用的 SparseDirectUMFPACK::vmult() 函数的作用。在这之后，我们必须确保我们也能解决解向量中的悬空节点的值，而这是用 AffineConstraints::distribute(). 来完成的。
 * 

 * 
 * 该函数需要一个额外的，但未使用的参数`tolerance`，它表示我们必须解决线性系统的精确程度。这个参数的含义在介绍中结合 "Eisenstat Walker技巧 "进行了讨论，但由于我们使用的是直接求解器而不是迭代求解器，所以我们并没有利用这个机会只求解线性系统的不精确性。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::solve(const Vector<double> &rhs, 
 *                                          Vector<double> &      solution, 
 *                                          const double /*tolerance*/) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "linear system solve"); 
 * 
 *     std::cout << "  Solving linear system" << std::endl; 
 * 
 *     jacobian_matrix_factorization->vmult(solution, rhs); 
 * 
 *     hanging_node_constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Refiningthemeshsettingboundaryvaluesandgeneratinggraphicaloutput"></a> 
 * <h4>Refining the mesh, setting boundary values, and generating graphical output</h4>
 * 

 * 
 * 以下三个函数又是对  step-15  中的函数的简单复制。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::refine_mesh() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(fe.degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       current_solution, 
 *       estimated_error_per_cell); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.03); 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 * 
 *     SolutionTransfer<dim> solution_transfer(dof_handler); 
 *     solution_transfer.prepare_for_coarsening_and_refinement(current_solution); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     Vector<double> tmp(dof_handler.n_dofs()); 
 *     solution_transfer.interpolate(current_solution, tmp); 
 *     current_solution = std::move(tmp); 
 * 
 *     hanging_node_constraints.clear(); 
 * 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 *     hanging_node_constraints.distribute(current_solution); 
 * 
 *     set_boundary_values(); 
 * 
 *     setup_system(/*initial_step=*/false); 
 *   } 
 * 
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::set_boundary_values() 
 *   { 
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              BoundaryValues<dim>(), 
 *                                              boundary_values); 
 *     for (const auto &boundary_value : boundary_values) 
 *       current_solution(boundary_value.first) = boundary_value.second; 
 * 
 *     hanging_node_constraints.distribute(current_solution); 
 *   } 
 * 
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::output_results( 
 *     const unsigned int refinement_cycle) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "graphical output"); 
 * 
 *     DataOut<dim> data_out; 
 * 
 *  
 *     data_out.add_data_vector(current_solution, "solution"); 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Therunfunctionandtheoveralllogicoftheprogram"></a> 
 * <h4>The run() function and the overall logic of the program</h4>
 * 

 * 
 * 这个程序中唯一**有趣的函数是驱动整个算法的函数，即从一个粗大的网格开始，做一些网格细化循环，并在每个网格上使用KINSOL来寻找我们从这个网格上离散化得到的非线性代数方程的解。上面的`refine_mesh()`函数可以确保一个网格上的解被用作下一个网格的起始猜测。我们还使用一个TimerOutput对象来测量每个网格上的每一次操作所花费的时间，并在每个周期开始时重置该计时器。
 * 

 * 
 * 正如在介绍中所讨论的，没有必要特别精确地解决粗略网格上的问题，因为这些问题只能作为下一个网格的起始猜测来解决。因此，我们将在 $k$ 个网格细化周期中使用 $\tau=10^{-3} \frac{1}{10^k}$ 的目标公差。
 * 

 * 
 * 所有这些都在这个函数的第一部分进行了编码。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MinimalSurfaceProblem<dim>::run() 
 *   { 
 *     GridGenerator::hyper_ball(triangulation); 
 *     triangulation.refine_global(2); 
 * 
 *     setup_system(/*initial_step=*/true); 
 *     set_boundary_values(); 
 * 
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 6; 
 *          ++refinement_cycle) 
 *       { 
 *         computing_timer.reset(); 
 *         std::cout << "Mesh refinement step " << refinement_cycle << std::endl; 
 * 
 *         if (refinement_cycle != 0) 
 *           refine_mesh(); 
 * 
 *         const double target_tolerance = 1e-3 * std::pow(0.1, refinement_cycle); 
 *         std::cout << "  Target_tolerance: " << target_tolerance << std::endl 
 *                   << std::endl; 
 * 
 * @endcode
 * 
 * 这就是有趣的开始。在顶部，我们创建了KINSOL求解器对象，并给它提供了一个对象，该对象编码了一些额外的具体情况（其中我们只改变了我们想要达到的非线性容忍度；但你可能想看看 SUNDIALS::KINSOL::AdditionalData 类有哪些其他成员，并与它们一起玩）。
 * 

 * 
 * 
 * @code
 *         { 
 *           typename SUNDIALS::KINSOL<Vector<double>>::AdditionalData 
 *             additional_data; 
 *           additional_data.function_tolerance = target_tolerance; 
 * 
 *           SUNDIALS::KINSOL<Vector<double>> nonlinear_solver(additional_data); 
 * 
 * @endcode
 * 
 * 然后，我们必须描述在介绍中已经提到的操作。从本质上讲，我们必须教KINSOL如何(i)将一个向量调整到正确的大小，(ii)计算残差向量，(iii)计算雅各布矩阵（在这期间我们也计算其因式分解），以及(iv)用雅各布矩阵解一个线性系统。
 * 

 * 
 * 所有这四种操作都由 SUNDIALS::KINSOL 类的成员变量表示，这些成员变量的类型是 `std::function`, ，即它们是我们可以分配给一个函数的指针的对象，或者像我们在这里做的那样，一个 "lambda函数"，它接受相应的参数并返回相应的信息。按照惯例，KINSOL希望做一些不重要的事情的函数返回一个整数，其中0表示成功。事实证明，我们只需用25行代码就可以完成所有这些工作。
 * 

 * 
 * 如果你不知道什么是 "lambda函数"，可以看看 step-12 或[wikipedia页面](https:en.wikipedia.org/wiki/Anonymous_function)关于这个问题。lambda函数的想法是，人们想用一组参数来定义一个函数，但(i)不使它成为一个命名的函数，因为通常情况下，该函数只在一个地方使用，似乎没有必要给它一个全局名称；(ii)该函数可以访问存在于定义它的地方的一些变量，包括成员变量。lambda函数的语法很笨拙，但最终还是很有用的）。)
 * 

 * 
 * 在代码块的最后，我们告诉KINSOL去工作，解决我们的问题。从'residual'、'setup_jacobian'和'solve_jacobian_system'函数中调用的成员函数将向屏幕打印输出，使我们能够跟踪程序的进展情况。
 * 

 * 
 * 
 * @code
 *           nonlinear_solver.reinit_vector = [&](Vector<double> &x) { 
 *             x.reinit(dof_handler.n_dofs()); 
 *           }; 
 * 
 *           nonlinear_solver.residual = 
 *             [&](const Vector<double> &evaluation_point, 
 *                 Vector<double> &      residual) { 
 *               compute_residual(evaluation_point, residual); 
 * 
 *               return 0; 
 *             }; 
 * 
 *           nonlinear_solver.setup_jacobian = 
 *             [&](const Vector<double> &current_u, 
 *                 const Vector<double> & /*current_f*/) { 
 *               compute_and_factorize_jacobian(current_u); 
 * 
 *               return 0; 
 *             }; 
 * 
 *           nonlinear_solver.solve_with_jacobian = [&](const Vector<double> &rhs, 
 *                                                      Vector<double> &      dst, 
 *                                                      const double tolerance) { 
 *             this->solve(rhs, dst, tolerance); 
 * 
 *             return 0; 
 *           }; 
 * 
 *           nonlinear_solver.solve(current_solution); 
 *         } 
 * 
 * @endcode
 * 
 * 剩下的就只是内务整理了。将数据写入文件，以便进行可视化，并显示收集到的时间摘要，以便我们可以解释每个操作花了多长时间，执行的频率如何，等等。
 * 

 * 
 * 
 * @code
 *         output_results(refinement_cycle); 
 * 
 *         computing_timer.print_summary(); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step77 
 * 
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step77; 
 * 
 *       MinimalSurfaceProblem<2> laplace_problem_2d; 
 *       laplace_problem_2d.run(); 
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
 *   return 0; 
 * } 
 * 
 * @endcode
examples/step-77/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当运行该程序时，你得到的输出看起来像这样。

@code
Mesh refinement step 0
  Target_tolerance: 0.001


  Computing residual vector... norm=0.231202
  Computing Jacobian matrix
  Factorizing Jacobian matrix
  Solving linear system
  Computing residual vector... norm=0.231202
  Computing residual vector... norm=0.171585
  Solving linear system
  Computing residual vector... norm=0.171585
  Computing residual vector... norm=0.127245
  Computing residual vector... norm=0.0796471
  Solving linear system
  Computing residual vector... norm=0.0796471
  Computing residual vector... norm=0.0625301
  Solving linear system
  Computing residual vector... norm=0.0625301
  Computing residual vector... norm=0.0498864
  Solving linear system
  Computing residual vector... norm=0.0498864
  Computing residual vector... norm=0.0407765
  Solving linear system
  Computing residual vector... norm=0.0407765
  Computing residual vector... norm=0.0341589
  Solving linear system
  Computing residual vector... norm=0.0341589
  Computing residual vector... norm=0.0292867
  Solving linear system
  Computing residual vector... norm=0.0292867
  Computing residual vector... norm=0.0256309
  Computing residual vector... norm=0.0223448
  Solving linear system
  Computing residual vector... norm=0.0223448
  Computing residual vector... norm=0.0202797
  Computing residual vector... norm=0.0183817
  Solving linear system
  Computing residual vector... norm=0.0183817
  Computing residual vector... norm=0.0170464
  Computing residual vector... norm=0.0157967
  Computing Jacobian matrix
  Factorizing Jacobian matrix
  Solving linear system
  Computing residual vector... norm=0.0157967
  Computing residual vector... norm=0.0141572
  Computing residual vector... norm=0.012657
 Solving linear system
  Computing residual vector... norm=0.012657
  Computing residual vector... norm=0.0116863
  Computing residual vector... norm=0.0107696
  Solving linear system
  Computing residual vector... norm=0.0107696
  Computing residual vector... norm=0.0100986
  Computing residual vector... norm=0.00944829
  Computing residual vector... norm=0.00822576
  Solving linear system
  Computing residual vector... norm=0.00822576
  Computing residual vector... norm=0.00781983
  Computing residual vector... norm=0.00741619
  Computing residual vector... norm=0.00661792
  Solving linear system
  Computing residual vector... norm=0.00661792
  Computing residual vector... norm=0.00630571
  Computing residual vector... norm=0.00599457
  Computing residual vector... norm=0.00537663
  Solving linear system
  Computing residual vector... norm=0.00537663
  Computing residual vector... norm=0.00512813
  Computing residual vector... norm=0.00488033
  Computing residual vector... norm=0.00438751
  Computing residual vector... norm=0.00342052
  Solving linear system
  Computing residual vector... norm=0.00342052
  Computing residual vector... norm=0.00326581
  Computing residual vector... norm=0.00311176
  Computing residual vector... norm=0.00280617
  Computing residual vector... norm=0.00220992
  Solving linear system
  Computing residual vector... norm=0.00220992
  Computing residual vector... norm=0.00209976
  Computing residual vector... norm=0.00199943
  Solving linear system
  Computing residual vector... norm=0.00199942
  Computing residual vector... norm=0.00190953
  Computing residual vector... norm=0.00182005
  Computing residual vector... norm=0.00164259
  Computing residual vector... norm=0.00129652



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.192s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| assembling the Jacobian         |         2 |    0.0141s |       7.4% |
| assembling the residual         |        61 |     0.168s |        88% |
| factorizing the Jacobian        |         2 |    0.0016s |      0.83% |
| graphical output                |         1 |   0.00385s |         2% |
| linear system solve             |        19 |    0.0013s |      0.68% |
+---------------------------------+-----------+------------+------------+



Mesh refinement step 1
  Target_tolerance: 0.0001


  Computing residual vector... norm=0.0883422
  Computing Jacobian matrix
  Factorizing Jacobian matrix
  Solving linear system
  Computing residual vector... norm=0.0883422
  Computing residual vector... norm=0.0607066
  Solving linear system
  Computing residual vector... norm=0.0607066
  Computing residual vector... norm=0.0437266
  Solving linear system
  Computing residual vector... norm=0.0437266
  Computing residual vector... norm=0.0327999
  Solving linear system
  Computing residual vector... norm=0.0327999
  Computing residual vector... norm=0.0255418
  Solving linear system
  Computing residual vector... norm=0.0255417
  Computing residual vector... norm=0.0206042
  Solving linear system
  Computing residual vector... norm=0.0206042
  Computing residual vector... norm=0.0171602
  Solving linear system
  Computing residual vector... norm=0.0171602
  Computing residual vector... norm=0.014689
  Solving linear system


[...]
@endcode



通过查看第一个网格上的输出的前几行，应该最容易解释这种方式。

@code
Mesh refinement step 0
Mesh refinement step 0
  Target_tolerance: 0.001


  Computing residual vector... norm=0.231202
  Computing Jacobian matrix
  Factorizing Jacobian matrix
  Solving linear system
  Computing residual vector... norm=0.231202
  Computing residual vector... norm=0.171585
  Solving linear system
  Computing residual vector... norm=0.171585
  Computing residual vector... norm=0.127245
  Computing residual vector... norm=0.0796471
  Solving linear system
  Computing residual vector... norm=0.0796471
  ...
@endcode

现在的情况是这样的。

- 在第一次残差计算中，KINSOL计算残差以查看是否达到了所需的公差。答案是否定的，所以它要求用户程序计算雅各布矩阵（然后该函数还通过SparseDirectUMFPACK对矩阵进行因子化）。

- 然后KINSOL指示我们用这个矩阵和之前计算的残差向量来解决一个形式为 $J_k \, \delta U_k = -F_k$ 的线性系统。

- 然后就是确定我们要在这个方向上走多远，也就是做线搜索。为此，KINSOL要求我们计算不同步长 $\alpha_k$ 的残差向量 $F(U_k + \alpha_k \delta U_k)$  。对于上面的第一步，它在尝试了两次后找到了一个可接受的 $\alpha_k$ ，第二次则需要尝试三次。

- 在找到一个合适的更新解 $U_{k+1}$ 之后，这个过程被重复，只是现在KINSOL对当前的雅各布矩阵很满意，没有指示我们重新建立矩阵和它的因式分解，而是要求我们用同一个矩阵解决一个线性系统。

在每个网格细化周期结束时，程序也会将解决方案写入VTU文件，它看起来如下。   <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-77.solution.png" alt="">
    </td>
  </tr>
</table> 


该计划的主要收获信息如下。

- 这个解和我们在步骤15中计算的解是一样的，也就是说，%SUNDIALS的KINSOL包的接口确实做了它们应该做的事。这不应该是一个惊喜，但重要的一点是，我们不必自己花时间去实现高级非线性求解器所依据的复杂算法。

- KINSOL能够避免各种操作，比如在实际上没有必要的时候重建雅各布矩阵。将上述输出中的线性求解次数与我们重建雅各布矩阵和计算其因式分解的次数相比较，应该可以清楚地看到，这在计算时间上带来了非常可观的节省，而我们却不需要实现复杂的算法来确定何时需要重建这些信息。

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


除了我们在这里考虑的小问题之外，稀疏的直接求解器需要太多的时间和内存--我们需要一个迭代求解器，就像我们在许多其他程序中使用的那样。然而，在目前的情况下，构建一个昂贵的预处理程序（例如，一个几何或代数多重网格方法）的权衡是不同的。由于我们可以在许多线性求解中重复使用同一个矩阵，我们也可以对预处理程序做同样的处理，与我们只在单一线性求解中使用预处理程序相比，在构建一个好的预处理程序上投入更多的工作更容易被证明，就像在许多其他情况下一样。

但迭代求解器也提供了其他机会。例如（正如在介绍中简要讨论的那样），只要我们离实际的解还很远，我们可能不需要在早期的非线性迭代中解到非常高的精度（小公差）。这就是那里提到的Eisenstat-Walker技巧的基础。

KINSOL提供了做线性解的函数，有一个需要达到的目标公差。我们在上面的程序中忽略了它，因为我们使用的直接求解器不需要公差，而是精确地求解线性系统（当然是四舍五入），但是迭代求解器可以利用这种信息--事实上也应该如此。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-77.cc"
*/
