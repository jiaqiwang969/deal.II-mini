/**
@page step_55 The step-55 tutorial program
This tutorial depends on step-40, step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Optimalpreconditioners">Optimal preconditioners</a>
        <li><a href="#Thesolverandpreconditioner">The solver and preconditioner</a>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
        <li><a href="#Problemsetup">Problem setup</a>
        <li><a href="#Themainprogram">The main program</a>
        <li><a href="#SystemSetup">System Setup</a>
        <li><a href="#Assembly">Assembly</a>
        <li><a href="#Solving">Solving</a>
        <li><a href="#Therest">The rest</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#InvestigateTrilinositerations">Investigate Trilinos iterations</a>
        <li><a href="#SolvetheOseenprobleminsteadoftheStokessystem">Solve the Oseen problem instead of the Stokes system</a>
        <li><a href="#Adaptiverefinement">Adaptive refinement</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-55/doc/intro.dox

 <br> 

<i>This program was contributed by Timo Heister. Special thanks to Sander
Rhebergen for the inspiration to finally write this tutorial.


This material is based upon work partially supported by National Science
Foundation grant DMS1522191 and the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-0949446 and The University of California-Davis.


The authors would like to thank the Isaac Newton Institute for
Mathematical Sciences, Cambridge, for support and hospitality during
the programme Melt in the Mantle where work on this tutorial was
undertaken. This work was supported by EPSRC grant no EP/K032208/1.
</i>




 @note  作为这个程序的前提条件，你需要安装PETSc或Trilinos和p4est库。在<a href="../../readme.html"
target="body">README</a>文件中描述了deal.II与这些附加库的安装情况。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


在第40步的基础上，本教程展示了如何使用MPI与PETSc或Trilinos进行线性代数，并行解决具有多个组件的线性PDEs。为此，我们返回到步骤22中讨论的斯托克斯方程。编写本教程的动机是在第40步（并行拉普拉斯）和第32步（针对时间相关问题的并行耦合斯托克斯与布西尼斯克）之间提供一个中间步骤（双关）。

本教程的学习成果是。

- 你能够并行地解决有多个变量的PDEs，并能将其应用于不同的问题。

- 你了解最佳预处理程序的概念，并能对某一特定问题进行检查。

- 你能够使用免费的计算机algreba系统SymPy（https://sympy.org）来构建制造的解决方案。

- 你可以为并行程序实现各种其他任务：错误计算、编写图形输出等。

- 你可以将矢量场、流线和矢量的轮廓可视化。

我们要解决的是满足斯托克斯方程的速度 $\textbf{u}$ 和压力 $p$ ，其内容为

@f{eqnarray*}


  - \triangle \textbf{u} + \nabla p &=& \textbf{f}, \\


  -\textrm{div}\; \textbf{u} &=& 0.


@f}






<a name="Optimalpreconditioners"></a><h3>Optimal preconditioners</h3>


请确保你阅读（甚至更好：尝试）步骤22中 "可能的扩展 "部分的 "块舒尔补码预处理 "所描述的内容。就像那里描述的那样，我们将使用Krylov方法和块状预处理程序来解决块状系统。

我们的目标是为线性系统构造一个非常简单的（也许是最简单的）最优预处理。如果预处理系统的迭代次数与网格大小无关，则该预处理程序被称为 "最优 "或 "最优复杂性" $h$  。你可以把这个定义扩展到要求与使用的处理器数量无关（我们将在结果部分讨论这个问题），计算域和网格质量，测试案例本身，有限元空间的多项式程度，等等。

为什么恒定的迭代次数被认为是 "最佳 "的？假设离散化的PDE给出一个有N个未知数的线性系统。因为来自有限元离散化的矩阵是稀疏的，矩阵-向量乘积可以在O(N)时间内完成。先决条件的应用充其量也只能是O(N)（例如可以用多网格方法来做）。如果解决线性系统所需的迭代次数与 $h$ 无关（因此也与N无关），那么解决该系统的总成本将是O(N)。不可能战胜这个复杂度，因为即使是查看右手边的所有条目也已经需要O(N)的时间。更多信息见  @cite elman2005  ，第2.5章（多网格）。

这里描述的预处理程序甚至比步骤22中描述的更简单，通常需要更多的迭代，因此需要更多的时间来解决。在考虑预处理程序时，最优性并不是唯一重要的衡量标准。但是一个最优的、昂贵的预处理程序通常比一个更便宜的、非最优的预处理程序更可取。这是因为，最终，随着网格尺寸越来越小，线性问题越来越大，前者将最终击败后者。

<a name="Thesolverandpreconditioner"></a><h3>The solver and preconditioner</h3>


我们对线性系统进行预处理

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{c}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{c}
    F \\ 0
  \end{array}\right),


@f}



块状对角线预处理器

@f{eqnarray*}
  P^{-1}
  =
  \left(\begin{array}{cc}
    A & 0 \\ 0 & S
  \end{array}\right) ^{-1},
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ 0 & S^{-1}
  \end{array}\right),


@f}

其中 $S=-BA^{-1} B^T$ 是舒尔补。

对于 $P$ 的这种选择，假设我们准确地处理了 $A^{-1}$ 和 $S^{-1}$ （这是一种 "理想化 "的情况），预处理的线性系统有三个独立于 $h$ 的特征值，因此是 "最优 "的。  见  @cite elman2005  中的6.2.1节（特别是第292页）。作为比较，在第22步中使用理想版的上块三角预处理（也用于第56步）会使所有的特征值都等于1。

我们将使用 $P^{-1}$ 中的逆运算的近似值，它（几乎）独立于 $h$ 。在这种情况下，我们可以再次证明，特征值是独立于 $h$ 的。对于Krylov方法，我们选择MINRES，它对分析很有吸引力（迭代次数被证明与 $h$ 无关，见上述书中第6.2.1章的其余部分），从计算的角度看很好（例如比GMRES更简单、更便宜），而且适用（矩阵和预处理器是对称的）。

对于近似，我们将使用压力空间中的质量矩阵的CG解来近似 $S^{-1}$  的作用。请注意，质量矩阵在光谱上等同于 $S$  。我们可以预期CG迭代的数量与 $h$ 无关，即使使用ILU这样的简单预处理程序。

对于速度块 $A$ 的近似，我们将执行一个单一的AMG V-循环。在实践中，这种选择并不完全独立于 $h$ ，这可以解释迭代数的轻微增加。一个可能的解释是，最粗的层次将被精确解决，而最粗的矩阵的层次数和大小是不可预测的。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们将根据经典的Kovasznay问题构建一个制造的解决方案，见  @cite kovasznay1948laminar  。这里是一个由x速度着色的解决方案的图像，包括速度的流线。

   <img src="https://www.dealii.org/images/steps/developer/step-55.solution.png" alt=""> 

不过，我们在这里必须作弊，因为我们不是在解决非线性的纳维-斯托克斯方程，而是解决没有对流项的线性斯托克斯系统。因此，为了重现完全相同的解，我们用科瓦兹内问题的解来制造解的方法。这将有效地把对流项移到右手边  $f$  。

右手边是用脚本 "reference.py "计算的，我们使用精确的解决方案来计算边界条件和误差。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/timer.h> 
 * 
 * @endcode
 * 
 * 下面这块出场代码与 step-40 相同，可以在PETSc和Trilinos之间切换。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/generic_linear_algebra.h> 
 * 
 * /* #define FORCE_USE_OF_TRILINOS */ 
 * 
 * 
 * 
 * namespace LA 
 * { 
 * #if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \ 
 *   !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS)) 
 *   using namespace dealii::LinearAlgebraPETSc; 
 * #  define USE_PETSC_LA 
 * #elif defined(DEAL_II_WITH_TRILINOS) 
 *   using namespace dealii::LinearAlgebraTrilinos; 
 * #else 
 * #  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required 
 * #endif 
 * } // namespace LA 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/solver_minres.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * #include <deal.II/lac/petsc_sparse_matrix.h> 
 * #include <deal.II/lac/petsc_vector.h> 
 * #include <deal.II/lac/petsc_solver.h> 
 * #include <deal.II/lac/petsc_precondition.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/lac/sparsity_tools.h> 
 * #include <deal.II/distributed/tria.h> 
 * #include <deal.II/distributed/grid_refinement.h> 
 * 
 * #include <cmath> 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step55 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * 我们需要一些辅助类来表示我们在介绍中描述的求解器策略。
 * 

 * 
 * 
 * @code
 *   namespace LinearSolvers 
 *   { 
 * 
 * @endcode
 * 
 * 这个类暴露了通过函数 InverseMatrix::vmult(). 应用给定矩阵的逆的动作，在内部，逆不是显式形成的。相反，一个带有CG的线性求解器被执行。这个类扩展了 step-22 中的InverseMatrix类，增加了一个指定预处理程序的选项，并允许在vmult函数中使用不同的矢量类型。
 * 

 * 
 * 
 * @code
 *     template <class Matrix, class Preconditioner> 
 *     class InverseMatrix : public Subscriptor 
 *     { 
 *     public: 
 *       InverseMatrix(const Matrix &m, const Preconditioner &preconditioner); 
 * 
 *       template <typename VectorType> 
 *       void vmult(VectorType &dst, const VectorType &src) const; 
 * 
 *     private: 
 *       const SmartPointer<const Matrix> matrix; 
 *       const Preconditioner &           preconditioner; 
 *     }; 
 * 
 *     template <class Matrix, class Preconditioner> 
 *     InverseMatrix<Matrix, Preconditioner>::InverseMatrix( 
 *       const Matrix &        m, 
 *       const Preconditioner &preconditioner) 
 *       : matrix(&m) 
 *       , preconditioner(preconditioner) 
 *     {} 
 * 
 *     template <class Matrix, class Preconditioner> 
 *     template <typename VectorType> 
 *     void 
 *     InverseMatrix<Matrix, Preconditioner>::vmult(VectorType &      dst, 
 *                                                  const VectorType &src) const 
 *     { 
 *       SolverControl solver_control(src.size(), 1e-8 * src.l2_norm()); 
 *       SolverCG<LA::MPI::Vector> cg(solver_control); 
 *       dst = 0; 
 * 
 *       try 
 *         { 
 *           cg.solve(*matrix, dst, src, preconditioner); 
 *         } 
 *       catch (std::exception &e) 
 *         { 
 *           Assert(false, ExcMessage(e.what())); 
 *         } 
 *     } 
 * 
 * @endcode
 * 
 * 该类是一个简单的2x2矩阵的块状对角线预处理器的模板类。
 * 

 * 
 * 
 * @code
 *     template <class PreconditionerA, class PreconditionerS> 
 *     class BlockDiagonalPreconditioner : public Subscriptor 
 *     { 
 *     public: 
 *       BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A, 
 *                                   const PreconditionerS &preconditioner_S); 
 * 
 *       void vmult(LA::MPI::BlockVector &      dst, 
 *                  const LA::MPI::BlockVector &src) const; 
 * 
 *     private: 
 *       const PreconditionerA &preconditioner_A; 
 *       const PreconditionerS &preconditioner_S; 
 *     }; 
 * 
 *     template <class PreconditionerA, class PreconditionerS> 
 *     BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>:: 
 *       BlockDiagonalPreconditioner(const PreconditionerA &preconditioner_A, 
 *                                   const PreconditionerS &preconditioner_S) 
 *       : preconditioner_A(preconditioner_A) 
 *       , preconditioner_S(preconditioner_S) 
 *     {} 
 * 
 *     template <class PreconditionerA, class PreconditionerS> 
 *     void BlockDiagonalPreconditioner<PreconditionerA, PreconditionerS>::vmult( 
 *       LA::MPI::BlockVector &      dst, 
 *       const LA::MPI::BlockVector &src) const 
 *     { 
 *       preconditioner_A.vmult(dst.block(0), src.block(0)); 
 *       preconditioner_S.vmult(dst.block(1), src.block(1)); 
 *     } 
 * 
 *  
 * @endcode
 * 
 * 
 * <a name="Problemsetup"></a> 
 * <h3>Problem setup</h3>
 * 

 * 
 * 下面的类代表测试问题的右手边和精确解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     RightHandSide() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  value) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   void RightHandSide<dim>::vector_value(const Point<dim> &p, 
 *                                         Vector<double> &  values) const 
 *   { 
 *     const double R_x = p[0]; 
 *     const double R_y = p[1]; 
 * 
 *     const double pi  = numbers::PI; 
 *     const double pi2 = pi * pi; 
 *     values[0] = 
 *       -1.0L / 2.0L * (-2 * sqrt(25.0 + 4 * pi2) + 10.0) * 
 *         exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) - 
 *       0.4 * pi2 * exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 
 *       0.1 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 2) * 
 *         exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi); 
 *     values[1] = 0.2 * pi * (-sqrt(25.0 + 4 * pi2) + 5.0) * 
 *                   exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) - 
 *                 0.05 * pow(-sqrt(25.0 + 4 * pi2) + 5.0, 3) * 
 *                   exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) / 
 *                   pi; 
 *     values[2] = 0; 
 *   } 
 * 
 *   template <int dim> 
 *   class ExactSolution : public Function<dim> 
 *   { 
 *   public: 
 *     ExactSolution() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  value) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   void ExactSolution<dim>::vector_value(const Point<dim> &p, 
 *                                         Vector<double> &  values) const 
 *   { 
 *     const double R_x = p[0]; 
 *     const double R_y = p[1]; 
 * 
 *     const double pi  = numbers::PI; 
 *     const double pi2 = pi * pi; 
 *     values[0] = 
 *       -exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * cos(2 * R_y * pi) + 1; 
 *     values[1] = (1.0L / 2.0L) * (-sqrt(25.0 + 4 * pi2) + 5.0) * 
 *                 exp(R_x * (-sqrt(25.0 + 4 * pi2) + 5.0)) * sin(2 * R_y * pi) / 
 *                 pi; 
 *     values[2] = 
 *       -1.0L / 2.0L * exp(R_x * (-2 * sqrt(25.0 + 4 * pi2) + 10.0)) - 
 *       2.0 * 
 *         (-6538034.74494422 + 
 *          0.0134758939981709 * exp(4 * sqrt(25.0 + 4 * pi2))) / 
 *         (-80.0 * exp(3 * sqrt(25.0 + 4 * pi2)) + 
 *          16.0 * sqrt(25.0 + 4 * pi2) * exp(3 * sqrt(25.0 + 4 * pi2))) - 
 *       1634508.68623606 * exp(-3.0 * sqrt(25.0 + 4 * pi2)) / 
 *         (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)) + 
 *       (-0.00673794699908547 * exp(sqrt(25.0 + 4 * pi2)) + 
 *        3269017.37247211 * exp(-3 * sqrt(25.0 + 4 * pi2))) / 
 *         (-8 * sqrt(25.0 + 4 * pi2) + 40.0) + 
 *       0.00336897349954273 * exp(1.0 * sqrt(25.0 + 4 * pi2)) / 
 *         (-10.0 + 2.0 * sqrt(25.0 + 4 * pi2)); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainprogram"></a> 
 * <h3>The main program</h3>
 * 

 * 
 * 主类与  step-40  非常相似，只是矩阵和向量现在是块状的，而且我们为拥有的和相关的DoF存储一个  std::vector<IndexSet>  ，而不是一个IndexSet。我们正好有两个IndexSets，一个用于所有速度未知数，一个用于所有压力未知数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class StokesProblem 
 *   { 
 *   public: 
 *     StokesProblem(unsigned int velocity_degree); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     unsigned int velocity_degree; 
 *     double       viscosity; 
 *     MPI_Comm     mpi_communicator; 
 * 
 *     FESystem<dim>                             fe; 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 *     DoFHandler<dim>                           dof_handler; 
 * 
 *     std::vector<IndexSet> owned_partitioning; 
 *     std::vector<IndexSet> relevant_partitioning; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     LA::MPI::BlockSparseMatrix system_matrix; 
 *     LA::MPI::BlockSparseMatrix preconditioner_matrix; 
 *     LA::MPI::BlockVector       locally_relevant_solution; 
 *     LA::MPI::BlockVector       system_rhs; 
 * 
 *     ConditionalOStream pcout; 
 *     TimerOutput        computing_timer; 
 *   }; 
 * 
 *   template <int dim> 
 *   StokesProblem<dim>::StokesProblem(unsigned int velocity_degree) 
 *     : velocity_degree(velocity_degree) 
 *     , viscosity(0.1) 
 *     , mpi_communicator(MPI_COMM_WORLD) 
 *     , fe(FE_Q<dim>(velocity_degree), dim, FE_Q<dim>(velocity_degree - 1), 1) 
 *     , triangulation(mpi_communicator, 
 *                     typename Triangulation<dim>::MeshSmoothing( 
 *                       Triangulation<dim>::smoothing_on_refinement | 
 *                       Triangulation<dim>::smoothing_on_coarsening)) 
 *     , dof_handler(triangulation) 
 *     , pcout(std::cout, 
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
 *     , computing_timer(mpi_communicator, 
 *                       pcout, 
 *                       TimerOutput::summary, 
 *                       TimerOutput::wall_times) 
 *   {} 
 * 
 * @endcode
 * 
 * Kovasnay流定义在域[-0.5, 1.5]^2上，我们通过将最小和最大值传递给 GridGenerator::hyper_cube. 来创建这个域。
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::make_grid() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -0.5, 1.5); 
 *     triangulation.refine_global(3); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SystemSetup"></a> 
 * <h3>System Setup</h3>
 * 

 * 
 * 与 step-40 相比，块矩阵和向量的构造是新的，与 step-22 这样的串行代码相比也是不同的，因为我们需要提供属于我们处理器的行的集合。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::setup_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "setup"); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 将所有的昏暗速度放入0区块，压力放入1区块，然后按区块重新排列未知数。最后计算每块有多少个未知数。
 * 

 * 
 * 
 * @code
 *     std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0); 
 *     stokes_sub_blocks[dim] = 1; 
 *     DoFRenumbering::component_wise(dof_handler, stokes_sub_blocks); 
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, stokes_sub_blocks); 
 * 
 *     const unsigned int n_u = dofs_per_block[0]; 
 *     const unsigned int n_p = dofs_per_block[1]; 
 * 
 *     pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << " (" 
 *           << n_u << '+' << n_p << ')' << std::endl; 
 * 
 * @endcode
 * 
 * 我们根据我们想要创建块状矩阵和向量的方式，将本地拥有的和本地相关的DoF的IndexSet分割成两个IndexSets。
 * 

 * 
 * 
 * @code
 *     owned_partitioning.resize(2); 
 *     owned_partitioning[0] = dof_handler.locally_owned_dofs().get_view(0, n_u); 
 *     owned_partitioning[1] = 
 *       dof_handler.locally_owned_dofs().get_view(n_u, n_u + n_p); 
 * 
 *     IndexSet locally_relevant_dofs; 
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
 *     relevant_partitioning.resize(2); 
 *     relevant_partitioning[0] = locally_relevant_dofs.get_view(0, n_u); 
 *     relevant_partitioning[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p); 
 * 
 * @endcode
 * 
 * 设置边界条件和悬挂节点的约束与  step-40  相同。尽管我们没有任何悬空节点，因为我们只进行全局细化，但把这个函数调用放进去仍然是个好主意，以备以后引入自适应细化。
 * 

 * 
 * 
 * @code
 *     { 
 *       constraints.reinit(locally_relevant_dofs); 
 * 
 *       FEValuesExtractors::Vector velocities(0); 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                ExactSolution<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 *       constraints.close(); 
 *     } 
 * 
 * @endcode
 * 
 * 现在我们根据BlockDynamicSparsityPattern来创建系统矩阵。我们知道我们不会有不同速度分量之间的耦合（因为我们使用的是拉普拉斯而不是变形张量），也不会有压力与其测试函数之间的耦合，所以我们使用一个表来将这个耦合信息传达给  DoFTools::make_sparsity_pattern.  。
 * 
 * @code
 *     { 
 *       system_matrix.clear(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (c == dim && d == dim) 
 *             coupling[c][d] = DoFTools::none; 
 *           else if (c == dim || d == dim || c == d) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         dof_handler, coupling, dsp, constraints, false); 
 * 
 *       SparsityTools::distribute_sparsity_pattern( 
 *         dsp, 
 *         dof_handler.locally_owned_dofs(), 
 *         mpi_communicator, 
 *         locally_relevant_dofs); 
 * 
 *       system_matrix.reinit(owned_partitioning, dsp, mpi_communicator); 
 *     } 
 * 
 * @endcode
 * 
 * 先决条件矩阵有不同的耦合（我们只在1,1块中填入质量矩阵），否则这段代码与上面的system_matrix的构造是相同的。
 * 

 * 
 * 
 * @code
 *     { 
 *       preconditioner_matrix.clear(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (c == dim && d == dim) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         dof_handler, coupling, dsp, constraints, false); 
 *       SparsityTools::distribute_sparsity_pattern( 
 *         dsp, 
 *         Utilities::MPI::all_gather(mpi_communicator, 
 *                                    dof_handler.locally_owned_dofs()), 
 *         mpi_communicator, 
 *         locally_relevant_dofs); 
 *       preconditioner_matrix.reinit(owned_partitioning, 
 * 
 * @endcode
 * 
 * owned_partitioning。
 * 

 * 
 * 
 * @code
 *                                    dsp, 
 *                                    mpi_communicator); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，我们以正确的尺寸构建块状向量。带有两个 std::vector<IndexSet> 的函数调用将创建一个重影向量。
 * 

 * 
 * 
 * @code
 *     locally_relevant_solution.reinit(owned_partitioning, 
 *                                      relevant_partitioning, 
 *                                      mpi_communicator); 
 *     system_rhs.reinit(owned_partitioning, mpi_communicator); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Assembly"></a> 
 * <h3>Assembly</h3>
 * 

 * 
 * 这个函数将系统矩阵、预处理矩阵和右手边集合起来。其代码非常标准。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::assemble_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "assembly"); 
 * 
 *     system_matrix         = 0; 
 *     preconditioner_matrix = 0; 
 *     system_rhs            = 0; 
 * 
 *     const QGauss<dim> quadrature_formula(velocity_degree + 1); 
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
 *     FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     const RightHandSide<dim>    right_hand_side; 
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 
 * 
 *     std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
 *     std::vector<double>         div_phi_u(dofs_per_cell); 
 *     std::vector<double>         phi_p(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 *     const FEValuesExtractors::Vector     velocities(0); 
 *     const FEValuesExtractors::Scalar     pressure(dim); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           cell_matrix  = 0; 
 *           cell_matrix2 = 0; 
 *           cell_rhs     = 0; 
 * 
 *           fe_values.reinit(cell); 
 *           right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
 *                                             rhs_values); 
 *           for (unsigned int q = 0; q < n_q_points; ++q) 
 *             { 
 *               for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                 { 
 *                   grad_phi_u[k] = fe_values[velocities].gradient(k, q); 
 *                   div_phi_u[k]  = fe_values[velocities].divergence(k, q); 
 *                   phi_p[k]      = fe_values[pressure].value(k, q); 
 *                 } 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     { 
 *                       cell_matrix(i, j) += 
 *                         (viscosity * 
 *                            scalar_product(grad_phi_u[i], grad_phi_u[j]) - 
 *                          div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
 *                         fe_values.JxW(q); 
 * 
 *                       cell_matrix2(i, j) += 1.0 / viscosity * phi_p[i] * 
 *                                             phi_p[j] * fe_values.JxW(q); 
 *                     } 
 * 
 *                   const unsigned int component_i = 
 *                     fe.system_to_component_index(i).first; 
 *                   cell_rhs(i) += fe_values.shape_value(i, q) * 
 *                                  rhs_values[q](component_i) * fe_values.JxW(q); 
 *                 } 
 *             } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           constraints.distribute_local_to_global(cell_matrix, 
 *                                                  cell_rhs, 
 *                                                  local_dof_indices, 
 *                                                  system_matrix, 
 *                                                  system_rhs); 
 * 
 *           constraints.distribute_local_to_global(cell_matrix2, 
 *                                                  local_dof_indices, 
 *                                                  preconditioner_matrix); 
 *         } 
 * 
 *     system_matrix.compress(VectorOperation::add); 
 *     preconditioner_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Solving"></a> 
 * <h3>Solving</h3>
 * 

 * 
 * 这个函数用MINRES求解线性系统，如介绍中所述，对两个对角线块使用块状对角线预处理和AMG。预处理程序对0,0块应用v循环，对1,1块应用质量矩阵的CG（Schur补充）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::solve() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "solve"); 
 * 
 *     LA::MPI::PreconditionAMG prec_A; 
 *     { 
 *       LA::MPI::PreconditionAMG::AdditionalData data; 
 * 
 * #ifdef USE_PETSC_LA 
 *       data.symmetric_operator = true; 
 * #endif 
 *       prec_A.initialize(system_matrix.block(0, 0), data); 
 *     } 
 * 
 *     LA::MPI::PreconditionAMG prec_S; 
 *     { 
 *       LA::MPI::PreconditionAMG::AdditionalData data; 
 * 
 * #ifdef USE_PETSC_LA 
 *       data.symmetric_operator = true; 
 * #endif 
 *       prec_S.initialize(preconditioner_matrix.block(1, 1), data); 
 *     } 
 * 
 * @endcode
 * 
 * InverseMatrix用于解决质量矩阵的问题。
 * 

 * 
 * 
 * @code
 *     using mp_inverse_t = LinearSolvers::InverseMatrix<LA::MPI::SparseMatrix, 
 *                                                       LA::MPI::PreconditionAMG>; 
 *     const mp_inverse_t mp_inverse(preconditioner_matrix.block(1, 1), prec_S); 
 * 
 * @endcode
 * 
 * 这是在上面定义的各个块的预处理的基础上构造的块预处理。
 * 

 * 
 * 
 * @code
 *     const LinearSolvers::BlockDiagonalPreconditioner<LA::MPI::PreconditionAMG, 
 *                                                      mp_inverse_t> 
 *       preconditioner(prec_A, mp_inverse); 
 * 
 * @endcode
 * 
 * 有了这些，我们终于可以设置线性求解器并求解该系统。
 * 

 * 
 * 
 * @code
 *     SolverControl solver_control(system_matrix.m(), 
 *                                  1e-10 * system_rhs.l2_norm()); 
 * 
 *     SolverMinRes<LA::MPI::BlockVector> solver(solver_control); 
 * 
 *     LA::MPI::BlockVector distributed_solution(owned_partitioning, 
 *                                               mpi_communicator); 
 * 
 *     constraints.set_zero(distributed_solution); 
 * 
 *     solver.solve(system_matrix, 
 *                  distributed_solution, 
 *                  system_rhs, 
 *                  preconditioner); 
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations." 
 *           << std::endl; 
 * 
 *     constraints.distribute(distributed_solution); 
 * 
 * @endcode
 * 
 * 像在  step-56  中一样，我们减去平均压力，以便与我们的参考解决方案进行误差计算，该解决方案的平均值为零。
 * 

 * 
 * 
 * @code
 *     locally_relevant_solution = distributed_solution; 
 *     const double mean_pressure = 
 *       VectorTools::compute_mean_value(dof_handler, 
 *                                       QGauss<dim>(velocity_degree + 2), 
 *                                       locally_relevant_solution, 
 *                                       dim); 
 *     distributed_solution.block(1).add(-mean_pressure); 
 *     locally_relevant_solution.block(1) = distributed_solution.block(1); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Therest"></a> 
 * <h3>The rest</h3>
 * 

 * 
 * 其余处理网格细化、输出和主循环的代码非常标准。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::refine_grid() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "refine"); 
 * 
 *     triangulation.refine_global(); 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     { 
 *       const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1); 
 *       const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), 
 *                                                        dim + 1); 
 * 
 *       Vector<double> cellwise_errors(triangulation.n_active_cells()); 
 *       QGauss<dim>    quadrature(velocity_degree + 2); 
 * 
 *       VectorTools::integrate_difference(dof_handler, 
 *                                         locally_relevant_solution, 
 *                                         ExactSolution<dim>(), 
 *                                         cellwise_errors, 
 *                                         quadrature, 
 *                                         VectorTools::L2_norm, 
 *                                         &velocity_mask); 
 * 
 *       const double error_u_l2 = 
 *         VectorTools::compute_global_error(triangulation, 
 *                                           cellwise_errors, 
 *                                           VectorTools::L2_norm); 
 * 
 *       VectorTools::integrate_difference(dof_handler, 
 *                                         locally_relevant_solution, 
 *                                         ExactSolution<dim>(), 
 *                                         cellwise_errors, 
 *                                         quadrature, 
 *                                         VectorTools::L2_norm, 
 *                                         &pressure_mask); 
 * 
 *       const double error_p_l2 = 
 *         VectorTools::compute_global_error(triangulation, 
 *                                           cellwise_errors, 
 *                                           VectorTools::L2_norm); 
 * 
 *       pcout << "error: u_0: " << error_u_l2 << " p_0: " << error_p_l2 
 *             << std::endl; 
 *     } 
 * 
 *     std::vector<std::string> solution_names(dim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(locally_relevant_solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 * 
 *     LA::MPI::BlockVector interpolated; 
 *     interpolated.reinit(owned_partitioning, MPI_COMM_WORLD); 
 *     VectorTools::interpolate(dof_handler, ExactSolution<dim>(), interpolated); 
 * 
 *     LA::MPI::BlockVector interpolated_relevant(owned_partitioning, 
 *                                                relevant_partitioning, 
 *                                                MPI_COMM_WORLD); 
 *     interpolated_relevant = interpolated; 
 *     { 
 *       std::vector<std::string> solution_names(dim, "ref_u"); 
 *       solution_names.emplace_back("ref_p"); 
 *       data_out.add_data_vector(interpolated_relevant, 
 *                                solution_names, 
 *                                DataOut<dim>::type_dof_data, 
 *                                data_component_interpretation); 
 *     } 
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells()); 
 *     for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *       subdomain(i) = triangulation.locally_owned_subdomain(); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     data_out.write_vtu_with_pvtu_record( 
 *       "./", "solution", cycle, mpi_communicator, 2); 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::run() 
 *   { 
 * #ifdef USE_PETSC_LA 
 *     pcout << "Running using PETSc." << std::endl; 
 * #else 
 *     pcout << "Running using Trilinos." << std::endl; 
 * #endif 
 *     const unsigned int n_cycles = 5; 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
 *       { 
 *         pcout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           make_grid(); 
 *         else 
 *           refine_grid(); 
 * 
 *         setup_system(); 
 * 
 *         assemble_system(); 
 *         solve(); 
 * 
 *         if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32) 
 *           { 
 *             TimerOutput::Scope t(computing_timer, "output"); 
 *             output_results(cycle); 
 *           } 
 * 
 *         computing_timer.print_summary(); 
 *         computing_timer.reset(); 
 * 
 *         pcout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step55 
 * 
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step55; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       StokesProblem<2> problem(2); 
 *       problem.run(); 
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
examples/step-55/doc/results.dox



<a name="Results"></a><h1>Results</h1>


正如上面的讨论所预期的那样，迭代次数与处理器的数量无关，只与 $h$ 有非常小的关系。

 <table>
<tr>
  <th colspan="2">PETSc</th>
  <th colspan="8">number of processors</th>
</tr>
<tr>
  <th>cycle</th>
  <th>dofs</th>
  <th>1</th>
  <th>2</th>
  <th>4</th>
  <th>8</th>
  <th>16</th>
  <th>32</th>
  <th>64</th>
  <th>128</th>
</tr>
<tr>
  <td>0</td>
  <td>659</td>
  <td>49</td>
  <td>49</td>
  <td>49</td>
  <td>51</td>
  <td>51</td>
  <td>51</td>
  <td>49</td>
  <td>49</td>
</tr>
<tr>
  <td>1</td>
  <td>2467</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>52</td>
  <td>54</td>
  <td>54</td>
  <td>53</td>
</tr>
<tr>
  <td>2</td>
  <td>9539</td>
  <td>56</td>
  <td>56</td>
  <td>56</td>
  <td>54</td>
  <td>56</td>
  <td>56</td>
  <td>54</td>
  <td>56</td>
</tr>
<tr>
  <td>3</td>
  <td>37507</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>56</td>
  <td>57</td>
  <td>56</td>
</tr>
<tr>
  <td>4</td>
  <td>148739</td>
  <td>58</td>
  <td>59</td>
  <td>57</td>
  <td>59</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
  <td>57</td>
</tr>
<tr>
  <td>5</td>
  <td>592387</td>
  <td>60</td>
  <td>60</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
  <td>59</td>
</tr>
<tr>
  <td>6</td>
  <td>2364419</td>
  <td>62</td>
  <td>62</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
  <td>61</td>
</tr>
</table> 

 <table>
<tr>
  <th colspan="2">Trilinos</th>
  <th colspan="8">number of processors</th>
</tr>
<tr>
  <th>cycle</th>
  <th>dofs</th>
  <th>1</th>
  <th>2</th>
  <th>4</th>
  <th>8</th>
  <th>16</th>
  <th>32</th>
  <th>64</th>
  <th>128</th>
</tr>
<tr>
  <td>0</td>
  <td>659</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
  <td>37</td>
</tr>
<tr>
  <td>1</td>
  <td>2467</td>
  <td>92</td>
  <td>89</td>
  <td>89</td>
  <td>82</td>
  <td>86</td>
  <td>81</td>
  <td>78</td>
  <td>78</td>
</tr>
<tr>
  <td>2</td>
  <td>9539</td>
  <td>102</td>
  <td>99</td>
  <td>96</td>
  <td>95</td>
  <td>95</td>
  <td>88</td>
  <td>83</td>
  <td>95</td>
</tr>
<tr>
  <td>3</td>
  <td>37507</td>
  <td>107</td>
  <td>105</td>
  <td>104</td>
  <td>99</td>
  <td>100</td>
  <td>96</td>
  <td>96</td>
  <td>90</td>
</tr>
<tr>
  <td>4</td>
  <td>148739</td>
  <td>112</td>
  <td>112</td>
  <td>111</td>
  <td>111</td>
  <td>127</td>
  <td>126</td>
  <td>115</td>
  <td>117</td>
</tr>
<tr>
  <td>5</td>
  <td>592387</td>
  <td>116</td>
  <td>115</td>
  <td>114</td>
  <td>112</td>
  <td>118</td>
  <td>120</td>
  <td>131</td>
  <td>130</td>
</tr>
<tr>
  <td>6</td>
  <td>2364419</td>
  <td>130</td>
  <td>126</td>
  <td>120</td>
  <td>120</td>
  <td>121</td>
  <td>122</td>
  <td>121</td>
  <td>123</td>
</tr>
</table> 

虽然PETSc的结果显示迭代次数不变，但使用Trilinos时，迭代次数增加。这可能是由于AMG预处理程序的不同设置造成的。出于性能方面的考虑，我们不允许在几千个未知数以下进行粗化。由于粗解器是精确求解（我们默认使用LU），层数的变化将影响V型循环的质量。因此，对于较小的问题规模，V型循环更接近于精确求解器。

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="InvestigateTrilinositerations"></a><h4>Investigate Trilinos iterations</h4>


玩弄平滑器、平滑步骤和Trilinos AMG的其他属性，以实现最佳预处理。

<a name="SolvetheOseenprobleminsteadoftheStokessystem"></a><h4>Solve the Oseen problem instead of the Stokes system</h4>


这一变化需要将外部求解器改为GMRES或BiCGStab，因为系统不再是对称的了。

你可以在对流项 $b
\cdot \nabla u$ 中规定精确的流动解，即 $b$  。如果你把右手边设置为零，这应该可以得到与原问题相同的解。

<a name="Adaptiverefinement"></a><h4>Adaptive refinement</h4>


到目前为止，这个教程程序在每一步都会对网格进行全局细化。将 StokesProblem::refine_grid() 中的代码替换为如下内容

@code
Vector<float> estimated_error_per_cell(triangulation.n_active_cells());


FEValuesExtractors::Vector velocities(0);
KellyErrorEstimator<dim>::estimate(
  dof_handler,
  QGauss<dim - 1>(fe.degree + 1),
  std::map<types::boundary_id, const Function<dim> *>(),
  locally_relevant_solution,
  estimated_error_per_cell,
  fe.component_mask(velocities));
parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
  triangulation, estimated_error_per_cell, 0.3, 0.0);
triangulation.execute_coarsening_and_refinement();
@endcode

使得探索自适应网格细化变得简单。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-55.cc"
*/
