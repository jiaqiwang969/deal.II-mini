/**
@page step_50 The step-50 tutorial program
This tutorial depends on step-16, step-37.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The testcase</a>
        <li><a href="#Workloadimbalanceforgeometricmultigridmethods">Workload imbalance for geometric multigrid methods</a>
        <li><a href="#Workloadimbalanceforalgebraicmultigridmethods">Workload imbalance for algebraic multigrid methods</a>
        <li><a href="#Runningtheprogram">Running the program</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Coefficientsandhelperclasses">Coefficients and helper classes</a>
        <li><a href="#Runtimeparameters">Run time parameters</a>
        <li><a href="#LaplaceProblemclass">LaplaceProblem class</a>
      <ul>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system()</a>
        <li><a href="#LaplaceProblemsetup_multigrid">LaplaceProblem::setup_multigrid()</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system()</a>
        <li><a href="#LaplaceProblemassemble_multigrid">LaplaceProblem::assemble_multigrid()</a>
        <li><a href="#LaplaceProblemassemble_rhs">LaplaceProblem::assemble_rhs()</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve()</a>
      </ul>
        <li><a href="#Theerrorestimator">The error estimator</a>
      <ul>
        <li><a href="#LaplaceProblemrefine_grid">LaplaceProblem::refine_grid()</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results()</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run()</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Testingconvergenceandhigherorderelements"> Testing convergence and higher order elements </a>
        <li><a href="#Coarsesolver"> Coarse solver </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-50/doc/intro.dox

 <br> 

<i>
This program was contributed by Thomas C. Clevenger and Timo Heister.
<br>
This material is based upon work partly supported by the National
Science Foundation Award DMS-2028346, OAC-2015848, EAR-1925575, by the Computational
Infrastructure in Geodynamics initiative (CIG), through the NSF under Award
EAR-0949446 and EAR-1550901 and The University of California -- Davis.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.4004166,https://zenodo.org/badge/DOI/10.5281/zenodo.4004166.svg} 

 @note  作为这个程序的前提条件，你需要同时安装p4est和PETSc或Trilinos库。在<a href="../../readme.html" target="body">README</a>文件中描述了deal.II和这些附加库的安装情况。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这个例子显示了deal.II中的多级函数在并行、分布式网格上的应用，并给出了几何和代数多栅方法的比较。代数多网格(AMG)的前置条件与step-40中使用的相同。考虑了两种几何多网格（GMG）预处理方法：一种是类似于步骤16的基于矩阵的版本（但用于并行计算），另一种是步骤37中讨论的无矩阵版本。我们的目标是找出哪种方法能够为大型并行计算提供最佳解算器。

本教程是基于  @cite clevenger_par_gmg  中的一个数值例子。关于deal.II中多网格实现的详细背景，请参见该出版物。我们将在下面的文字中总结一些结果。

代数多网格方法显然是最容易用deal.II实现的，因为诸如 TrilinosWrappers::PreconditionAMG 和 PETScWrappers::PreconditionBoomerAMG 这样的类本质上是黑盒子预处理程序，即使是并行计算，也只需要几行就能设置好。另一方面，几何多网格方法需要对整个代码库进行修改 -- 不是很多，但必须知道自己在做什么。

这个程序的结果将显示，代数和几何多网格方法的性能大致相当<i>when using matrix-based formulations</i>，而无矩阵的几何多网格方法对于这里所考虑的问题要好很多。另一个结论是，当每个处理器的未知数小于20,000个时，基于矩阵的几何多网格方法真的不能很好地扩展。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们考虑变系数拉普拉斯的弱表述

@f{align*}
 (\epsilon \nabla u, \nabla v) = (f,v) \quad \forall v \in V_h


@f}

在域 $\Omega = [-1,1]^\text{dim} \setminus [0,1]^\text{dim}$ （二维的L形域和三维的Fichera角）上，如果 $\min(x,y,z)>-\frac{1}{2}$ ，则 $\epsilon = 100$ 。换句话说， $\epsilon$ 是沿着域的边缘或面跑到重入角的小，这在下图中会看到。

边界条件在整个边界上是 $u=0$ ，右手边是 $f=1$  。我们使用连续 $Q_2$ 元素来表示离散的有限元空间 $V_h$ ，并使用基于残差的、单元的后验误差估计器 $e(K) = e_{\text{cell}}(K) + e_{\text{face}}(K)$ ，来自 @cite karakashian2003posteriori 的

@f{align*}
 e_{\text{cell}}(K) &= h^2 \| f + \epsilon \triangle u \|_K^2, \\
 e_{\text{face}}(K) &= \sum_F h_F \| \jump{ \epsilon \nabla u \cdot n } \|_F^2,


@f}

来适应性地细化网格。(这是KellyErrorEstimator类中使用的Kelly误差估计器的概括，KellyErrorEstimator类驱动大多数其他教程程序中的网格细化。)下图显示了二维的求解和细化：  <img width="400px" src="https://www.dealii.org/images/steps/developer/step-50-2d-solution.png" alt="">  在三维中，求解看起来类似（见下文）。在左边你可以看到解决方案，在右边我们显示了靠近域中心的 $x$ 的切片，显示了自适应细化的网格。   <table width="60%" align="center">
  <tr>
    <td align="center">
      <img width="400px" src="https://www.dealii.org/images/steps/developer/step-50-3d-solution.png" alt="">
    </td>
    <td align="center">
      <img width="400px" src="https://www.dealii.org/images/steps/developer/step-50-refinement.png" alt="">
    </td>
  </tr>
</table> 在二维和三维中，你都可以看到自适应细化拾取了角部奇点和粘度跳跃的内部奇点，而沿分离两个粘度的线的界面（正确地）没有被细化，因为它被充分地解决。这是因为由系数跳跃导致的解决方案中的扭结与细胞界面对齐。




<a name="Workloadimbalanceforgeometricmultigridmethods"></a><h3>Workload imbalance for geometric multigrid methods</h3>


如上所述，这个程序的目的是展示代数和几何多网格方法在这个问题上的应用，并做到并行计算。使算法扩展到大型并行机器的一个重要组成部分是确保每个处理器都有相同的工作量。更准确地说，重要的是没有一小部分处理器比其他处理器有更多的工作，因为如果是这样的话，很大一部分处理器会闲置，等待小部分处理器完成。相反，一小部分处理器的工作大大超过<i>less</i>并不是问题，因为大多数处理器继续生产，只有一小部分处理器在完成工作后闲置。)

对于活跃的网格，我们使用 parallel::distributed::Triangulation 类，正如在步骤40中所做的那样，它使用外部库<a href="http://www.p4est.org/">p4est</a>中的功能在处理器之间分配活跃单元。对于多级层次结构中的非活动单元，deal.II实现了我们所说的 "第一子规则"，对于层次结构中的每个单元，我们递归地将一个单元的父级分配给第一个子单元的所有者。下面的数字给出了这样一个分布的例子。这里的左图表示使用空间填充曲线划分的二维网格样本的活动单元（这也是p4est用来划分单元的方法）；中间的图片给出了活动网格的树状表示；右图给出了单元的多级层次结构。颜色和数字代表不同的处理器。树上的圆形节点是非活动单元，使用 "长子规则 "进行分配。

 <img width="800px" src="https://www.dealii.org/images/steps/developer/step-50-workload-example.png" alt=""> 

在这个例子中，屏幕上的输出包括一个 "分区效率 "的值，这个值由 MGTools::workload_imbalance(). 给出，将用 $\mathbb{E}$ 表示，量化了多网格层次结构中每一层没有完美的工作平衡所产生的开销。这种不平衡在上面的例子中很明显：虽然 $\ell=2$ 层在三个处理器的四个单元中尽可能的平衡，但粗略的 $\ell=0$ 层只有一个处理器有工作，而 $\ell=1$ 层只有两个处理器有工作，其中一个处理器的工作是另一个的三倍。

对于定义 $\mathbb{E}$ ，需要注意的是，由于我们使用局部平滑来定义多网格层次（参见 @ref mg_paper "多网格论文 "中对局部平滑的描述），一个单元的细化水平对应于该单元的多网格水平。现在，让 $N_{\ell}$ 为 $\ell$ 层的单元数（包括活动和非活动单元）， $N_{\ell,p}$ 为进程 $p$ 所拥有的子集。我们还将用 $P$ 表示处理器的总数量。假设任何一个处理器的工作量与该处理器拥有的单元格数量成正比，每个处理器的最佳工作量为

@f{align*}
W_{\text{opt}} = \frac1{P}\sum_{\ell} N_{\ell} = \sum_{\ell}\left(\frac1{P}\sum_{p}N_{\ell,p}\right).


@f}

接下来，假设每一层的工作都是同步的（即在V型循环的每一层，在进入下一层之前，所有的处理器都必须完成工作），每一层的极限工作由以下公式给出

@f{align*}
W_\ell = \max_{p} N_{\ell,p},


@f}

和总的并行复杂性

@f{align*}
W = \sum_{\ell} W_\ell.


@f}

然后我们将 $\mathbb{E}$ 定义为最佳分区与当前分区的并行复杂度之比

@f{align*}
  \mathbb{E} = \frac{W_{\text{opt}}}{W}.


@f}

对于上面的例子分布，我们有

@f{align*}
W_{\text{opt}}&=\frac{1}{P}\sum_{\ell} N_{\ell} = \frac{1}{3} \left(1+4+4\right)= 3 \qquad
\\
W &= \sum_\ell W_\ell = 1 + 2 + 3 = 6
\\
\mathbb{E} &= \frac{W_{\text{opt}}}{W} = \frac12.


@f}

这个值 MGTools::workload_imbalance()  $= 1/\mathbb{E}$ 代表了我们对GMG方法（vmults、assembly等）所期望的时间增加的因素，因为与完全负载平衡的工作负载相比，网格分区的不平衡。我们将在下面的结果部分报告一连串的网格，并与观察到的减速进行比较，因为我们的处理器数量越来越大（通常，负载不平衡也会变大）。

这些考虑在 @cite clevenger_par_gmg 中得到了更详细的考虑，其中包含了对分区效率模型和不平衡对GMG V周期时间的影响的全面讨论。总之， $\mathbb{E}$ 的值高度依赖于所使用的局部网格细化程度，对于全局细化的网格有一个最佳值 $\mathbb{E} \approx 1$ 。通常对于自适应细化的网格，用于分配单个网格的处理器数量对 $\mathbb{E}$ 有负面影响，但只到一个平移点，即处理器数量增加时，不平衡度保持相对稳定，进一步细化对 $\mathbb{E}$ 的影响很小。最后， $1/\mathbb{E}$ 被证明可以准确地表示出对V型周期的计时所预期的并行扩展的减慢。

应该注意的是，在多级网格之间有可能存在一些异步工作，特别是纯粹的近邻MPI通信，而且可以构建一个自适应网格，由于异步工作 "掩盖 "了不平衡，效率模型将远远高估V-周期的减慢（假设各级同步）。然而，对于大多数现实的自适应网格来说，预期这种异步工作只会掩盖非常小的一部分不平衡，效率模型会很好地描述减速。




<a name="Workloadimbalanceforalgebraicmultigridmethods"></a><h3>Workload imbalance for algebraic multigrid methods</h3>


上面的考虑表明，我们必须期待在deal.II中实现的几何多网格算法的可扩展性有一定的限制，因为即使在网格的最细层是完全负载平衡的情况下，较粗层也可能不是。同时，较粗层的权重较小（ $W_\ell$ 对 $W$ 的贡献较小），因为较粗层的单元较少，因此，对整体运行时间的贡献不如较细层。换句话说，较粗层次的不平衡可能不会导致大局的影响。

代数多网格方法当然是基于一种完全不同的方法来创建层次结构的。特别是，他们纯粹是在分析系统矩阵的基础上创建这些层次，并且在作为 TrilinosWrappers::PreconditionAMG 和 PETScWrappers::PreconditionBoomerAMG 类基础的hypre和ML/MueLu包中都实现了非常复杂的算法，以确保问题在每个层次上都得到良好的负载平衡。在某种意义上，这些算法比几何多网格方法更简单，因为它们只处理矩阵本身，而不是所有的网格、邻居、父母和其他几何实体的内涵。同时，为了使代数多网格方法能够扩展到非常大的问题，人们也做了很多工作，包括将在某一层次上工作的处理器数量减少到所有处理器的一个子集，如果不这样的话，处理器花在计算上的时间会比花在通信上的时间少。(人们可能会注意到，在几何多网格算法中也有可能实现这些相同的想法，在这些算法中，人们有目的地将一些处理器闲置在较粗的层次上，以减少通信量。只是目前deal.II没有这样做。)

然而，这些并不是我们在这里通常需要担心的问题。在大多数情况下，我们使用代数多网格方法作为黑箱方法。




<a name="Runningtheprogram"></a><h3>Running the program</h3>


如上所述，这个程序可以使用三种不同的方式来求解线性系统：基于矩阵的几何多网格（"MB"），无矩阵几何多网格（"MF"）和代数多网格（"AMG"）。这个程序所在的目录有后缀为".prm "的输入文件，适用于所有这三种选项，以及2D和3D。

你可以按以下方式执行该程序

@code
  ./step-50 gmg_mb_2d.prm
@endcode

而这将从给定的输入文件（这里是`mg_mb_2d.prm`）中获取运行时参数。

该程序的目的是要并行运行，你可以使用诸如以下的命令来实现这一点

@code
  mpirun -np 4 ./step-50 gmg_mb_2d.prm
@endcode

如果你想，比如说，在四个处理器上运行。也就是说，如果你有多少个处理器，程序也可以在`-np 28672`下运行）。


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
 * 包含文件是  step-40  ,  step-16  , 和  step-37  的组合。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/data_out_base.h> 
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/distributed/grid_refinement.h> 
 * #include <deal.II/distributed/tria.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * 
 * @endcode
 * 
 * 我们使用与 step-40 相同的策略，在PETSc和Trilinos之间进行切换。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/generic_linear_algebra.h> 
 * 
 * @endcode
 * 
 * 如果你已经安装了PETSc和Trilinos，并且你喜欢在本例中使用PETSc，请将下面的预处理程序定义注释进去或退出。
 * 

 * 
 * 
 * @code
 * #define FORCE_USE_OF_TRILINOS 
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
 * #include <deal.II/matrix_free/matrix_free.h> 
 * #include <deal.II/matrix_free/operators.h> 
 * #include <deal.II/matrix_free/fe_evaluation.h> 
 * #include <deal.II/multigrid/mg_coarse.h> 
 * #include <deal.II/multigrid/mg_constrained_dofs.h> 
 * #include <deal.II/multigrid/mg_matrix.h> 
 * #include <deal.II/multigrid/mg_smoother.h> 
 * #include <deal.II/multigrid/mg_tools.h> 
 * #include <deal.II/multigrid/mg_transfer.h> 
 * #include <deal.II/multigrid/multigrid.h> 
 * #include <deal.II/multigrid/mg_transfer_matrix_free.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * @endcode
 * 
 * 以下文件用于组装误差估计器，如  step-12  。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_interface_values.h> 
 * #include <deal.II/meshworker/mesh_loop.h> 
 * 
 * using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Coefficientsandhelperclasses"></a> 
 * <h3>Coefficients and helper classes</h3>
 * 

 * 
 * MatrixFree运算符必须使用 dealii::LinearAlgebra::distributed::Vector 矢量类型。这里我们定义了复制到Trilinos向量的操作，以便与基于矩阵的代码兼容。请注意，目前PETSc矢量类型不存在这种功能，所以必须安装Trilinos来使用本教程中的MatrixFree求解器。
 * 

 * 
 * 
 * @code
 * namespace ChangeVectorTypes 
 * { 
 *   template <typename number> 
 *   void copy(LA::MPI::Vector &                                         out, 
 *             const dealii::LinearAlgebra::distributed::Vector<number> &in) 
 *   { 
 *     dealii::LinearAlgebra::ReadWriteVector<double> rwv( 
 *       out.locally_owned_elements()); 
 *     rwv.import(in, VectorOperation::insert); 
 * #ifdef USE_PETSC_LA 
 *     AssertThrow(false, 
 *                 ExcMessage("CopyVectorTypes::copy() not implemented for " 
 *                            "PETSc vector types.")); 
 * #else 
 *     out.import(rwv, VectorOperation::insert); 
 * #endif 
 *   } 
 * 
 *   template <typename number> 
 *   void copy(dealii::LinearAlgebra::distributed::Vector<number> &out, 
 *             const LA::MPI::Vector &                             in) 
 *   { 
 *     dealii::LinearAlgebra::ReadWriteVector<double> rwv; 
 * #ifdef USE_PETSC_LA 
 *     (void)in; 
 *     AssertThrow(false, 
 *                 ExcMessage("CopyVectorTypes::copy() not implemented for " 
 *                            "PETSc vector types.")); 
 * #else 
 *     rwv.reinit(in); 
 * #endif 
 *     out.import(rwv, VectorOperation::insert); 
 *   } 
 * } // namespace ChangeVectorTypes 
 * 
 * @endcode
 * 
 * 让我们继续描述我们要解决的问题。我们把右边的函数设置为1.0。 @p value 函数返回一个VectorizedArray，被无矩阵代码路径所使用。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * class RightHandSide : public Function<dim> 
 * { 
 * public: 
 *   virtual double value(const Point<dim> & /*p*/, 
 *                        const unsigned int /*component*/ = 0) const override 
 *   { 
 *     return 1.0; 
 *   } 
 * 
 *   template <typename number> 
 *   VectorizedArray<number> 
 *   value(const Point<dim, VectorizedArray<number>> & /*p*/, 
 *         const unsigned int /*component*/ = 0) const 
 *   { 
 *     return VectorizedArray<number>(1.0); 
 *   } 
 * }; 
 * 
 * @endcode
 * 
 * 接下来的这个类表示扩散系数。我们使用一个可变的系数，在任何一个至少有一个坐标小于-0.5的点上是100.0，在所有其他点上是1.0。如上所述，一个单独的value()返回一个VectorizedArray，用于无矩阵代码。一个 @p average()函数计算了一组点的算术平均。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * class Coefficient : public Function<dim> 
 * { 
 * public: 
 *   virtual double value(const Point<dim> &p, 
 *                        const unsigned int /*component*/ = 0) const override; 
 * 
 *   template <typename number> 
 *   VectorizedArray<number> value(const Point<dim, VectorizedArray<number>> &p, 
 *                                 const unsigned int /*component*/ = 0) const; 
 * 
 *   template <typename number> 
 *   number average_value(const std::vector<Point<dim, number>> &points) const; 
 * 
 * @endcode
 * 
 * 当在MatrixFree框架中使用一个系数时，我们还需要一个函数，为MatrixFree运算符参数提供的一组单元格创建一个系数表。
 * 

 * 
 * 
 * @code
 *   template <typename number> 
 *   std::shared_ptr<Table<2, VectorizedArray<number>>> make_coefficient_table( 
 *     const MatrixFree<dim, number, VectorizedArray<number>> &mf_storage) const; 
 * }; 
 * 
 * template <int dim> 
 * double Coefficient<dim>::value(const Point<dim> &p, const unsigned int) const 
 * { 
 *   for (int d = 0; d < dim; ++d) 
 *     { 
 *       if (p[d] < -0.5) 
 *         return 100.0; 
 *     } 
 *   return 1.0; 
 * } 
 * 
 * template <int dim> 
 * template <typename number> 
 * VectorizedArray<number> 
 * Coefficient<dim>::value(const Point<dim, VectorizedArray<number>> &p, 
 *                         const unsigned int) const 
 * { 
 *   VectorizedArray<number> return_value = VectorizedArray<number>(1.0); 
 *   for (unsigned int i = 0; i < VectorizedArray<number>::size(); ++i) 
 *     { 
 *       for (int d = 0; d < dim; ++d) 
 *         if (p[d][i] < -0.5) 
 *           { 
 *             return_value[i] = 100.0; 
 *             break; 
 *           } 
 *     } 
 * 
 *   return return_value; 
 * } 
 * 
 * template <int dim> 
 * template <typename number> 
 * number Coefficient<dim>::average_value( 
 *   const std::vector<Point<dim, number>> &points) const 
 * { 
 *   number average(0); 
 *   for (unsigned int i = 0; i < points.size(); ++i) 
 *     average += value(points[i]); 
 *   average /= points.size(); 
 * 
 *   return average; 
 * } 
 * 
 * template <int dim> 
 * template <typename number> 
 * std::shared_ptr<Table<2, VectorizedArray<number>>> 
 * Coefficient<dim>::make_coefficient_table( 
 *   const MatrixFree<dim, number, VectorizedArray<number>> &mf_storage) const 
 * { 
 *   auto coefficient_table = 
 *     std::make_shared<Table<2, VectorizedArray<number>>>(); 
 * 
 *   FEEvaluation<dim, -1, 0, 1, number> fe_eval(mf_storage); 
 * 
 *   const unsigned int n_cells    = mf_storage.n_cell_batches(); 
 *   const unsigned int n_q_points = fe_eval.n_q_points; 
 * 
 *   coefficient_table->reinit(n_cells, 1); 
 * 
 *   for (unsigned int cell = 0; cell < n_cells; ++cell) 
 *     { 
 *       fe_eval.reinit(cell); 
 * 
 *       VectorizedArray<number> average_value = 0.; 
 *       for (unsigned int q = 0; q < n_q_points; ++q) 
 *         average_value += value(fe_eval.quadrature_point(q)); 
 *       average_value /= n_q_points; 
 * 
 *       (*coefficient_table)(cell, 0) = average_value; 
 *     } 
 * 
 *   return coefficient_table; 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameters"></a> 
 * <h3>Run time parameters</h3>
 * 

 * 
 * 我们将使用ParameterHandler来在运行时传入参数。 该结构 @p Settings 解析并存储这些参数，以便在整个程序中进行查询。
 * 

 * 
 * 
 * @code
 * struct Settings 
 * { 
 *   bool try_parse(const std::string &prm_filename); 
 * 
 *   enum SolverType 
 *   { 
 *     gmg_mb, 
 *     gmg_mf, 
 *     amg 
 *   }; 
 * 
 *   SolverType solver; 
 * 
 *   int          dimension; 
 *   double       smoother_dampen; 
 *   unsigned int smoother_steps; 
 *   unsigned int n_steps; 
 *   bool         output; 
 * }; 
 * 
 * bool Settings::try_parse(const std::string &prm_filename) 
 * { 
 *   ParameterHandler prm; 
 *   prm.declare_entry("dim", "2", Patterns::Integer(), "The problem dimension."); 
 *   prm.declare_entry("n_steps", 
 *                     "10", 
 *                     Patterns::Integer(0), 
 *                     "Number of adaptive refinement steps."); 
 *   prm.declare_entry("smoother dampen", 
 *                     "1.0", 
 *                     Patterns::Double(0.0), 
 *                     "Dampen factor for the smoother."); 
 *   prm.declare_entry("smoother steps", 
 *                     "1", 
 *                     Patterns::Integer(1), 
 *                     "Number of smoother steps."); 
 *   prm.declare_entry("solver", 
 *                     "MF", 
 *                     Patterns::Selection("MF|MB|AMG"), 
 *                     "Switch between matrix-free GMG, " 
 *                     "matrix-based GMG, and AMG."); 
 *   prm.declare_entry("output", 
 *                     "false", 
 *                     Patterns::Bool(), 
 *                     "Output graphical results."); 
 * 
 *   if (prm_filename.size() == 0) 
 *     { 
 *       std::cout << "****  Error: No input file provided!\n" 
 *                 << "****  Error: Call this program as './step-50 input.prm\n" 
 *                 << "\n" 
 *                 << "****  You may want to use one of the input files in this\n" 
 *                 << "****  directory, or use the following default values\n" 
 *                 << "****  to create an input file:\n"; 
 *       if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
 *         prm.print_parameters(std::cout, ParameterHandler::Text); 
 *       return false; 
 *     } 
 * 
 *   try 
 *     { 
 *       prm.parse_input(prm_filename); 
 *     } 
 *   catch (std::exception &e) 
 *     { 
 *       if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
 *         std::cerr << e.what() << std::endl; 
 *       return false; 
 *     } 
 * 
 *   if (prm.get("solver") == "MF") 
 *     this->solver = gmg_mf; 
 *   else if (prm.get("solver") == "MB") 
 *     this->solver = gmg_mb; 
 *   else if (prm.get("solver") == "AMG") 
 *     this->solver = amg; 
 *   else 
 *     AssertThrow(false, ExcNotImplemented()); 
 * 
 *   this->dimension       = prm.get_integer("dim"); 
 *   this->n_steps         = prm.get_integer("n_steps"); 
 *   this->smoother_dampen = prm.get_double("smoother dampen"); 
 *   this->smoother_steps  = prm.get_integer("smoother steps"); 
 *   this->output          = prm.get_bool("output"); 
 * 
 *   return true; 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemclass"></a> 
 * <h3>LaplaceProblem class</h3>
 * 

 * 
 * 这是该程序的主类。它看起来与  step-16  ,  step-37  , 和  step-40  非常相似。对于MatrixFree的设置，我们使用 MatrixFreeOperators::LaplaceOperator 类，它在内部定义了`local_apply()`, `compute_diagonal()`, 和`set_coefficient()`函数。请注意，多项式的度数是这个类的一个模板参数。这对无矩阵代码来说是必要的。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * class LaplaceProblem 
 * { 
 * public: 
 *   LaplaceProblem(const Settings &settings); 
 *   void run(); 
 * 
 * private: 
 * 
 * @endcode
 * 
 * 我们将在整个程序中使用以下类型。首先是基于矩阵的类型，之后是无矩阵的类。对于无矩阵的实现，我们使用 @p float 作为水平运算符。
 * 

 * 
 * 
 * @code
 *   using MatrixType         = LA::MPI::SparseMatrix; 
 *   using VectorType         = LA::MPI::Vector; 
 *   using PreconditionAMG    = LA::MPI::PreconditionAMG; 
 *   using PreconditionJacobi = LA::MPI::PreconditionJacobi; 
 * 
 *   using MatrixFreeLevelMatrix = MatrixFreeOperators::LaplaceOperator< 
 *     dim, 
 *     degree, 
 *     degree + 1, 
 *     1, 
 *     LinearAlgebra::distributed::Vector<float>>; 
 *   using MatrixFreeActiveMatrix = MatrixFreeOperators::LaplaceOperator< 
 *     dim, 
 *     degree, 
 *     degree + 1, 
 *     1, 
 *     LinearAlgebra::distributed::Vector<double>>; 
 * 
 *   using MatrixFreeLevelVector  = LinearAlgebra::distributed::Vector<float>; 
 *   using MatrixFreeActiveVector = LinearAlgebra::distributed::Vector<double>; 
 * 
 *   void setup_system(); 
 *   void setup_multigrid(); 
 *   void assemble_system(); 
 *   void assemble_multigrid(); 
 *   void assemble_rhs(); 
 *   void solve(); 
 *   void estimate(); 
 *   void refine_grid(); 
 *   void output_results(const unsigned int cycle); 
 * 
 *   Settings settings; 
 * 
 *   MPI_Comm           mpi_communicator; 
 *   ConditionalOStream pcout; 
 * 
 *   parallel::distributed::Triangulation<dim> triangulation; 
 *   const MappingQ1<dim>                      mapping; 
 *   FE_Q<dim>                                 fe; 
 * 
 *   DoFHandler<dim> dof_handler; 
 * 
 *  
 *   IndexSet                  locally_relevant_dofs; 
 *   AffineConstraints<double> constraints; 
 * 
 *   MatrixType             system_matrix; 
 *   MatrixFreeActiveMatrix mf_system_matrix; 
 *   VectorType             solution; 
 *   VectorType             right_hand_side; 
 *   Vector<double>         estimated_error_square_per_cell; 
 * 
 *   MGLevelObject<MatrixType> mg_matrix; 
 *   MGLevelObject<MatrixType> mg_interface_in; 
 *   MGConstrainedDoFs         mg_constrained_dofs; 
 * 
 *   MGLevelObject<MatrixFreeLevelMatrix> mf_mg_matrix; 
 * 
 *   TimerOutput computing_timer; 
 * }; 
 * 
 * @endcode
 * 
 * 关于构造函数的唯一有趣的部分是，除非我们使用AMG，否则我们会构造多网格的层次结构。为此，我们需要在这个构造函数完成之前解析运行时参数。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * LaplaceProblem<dim, degree>::LaplaceProblem(const Settings &settings) 
 *   : settings(settings) 
 *   , mpi_communicator(MPI_COMM_WORLD) 
 *   , pcout(std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
 *   , triangulation(mpi_communicator, 
 *                   Triangulation<dim>::limit_level_difference_at_vertices, 
 *                   (settings.solver == Settings::amg) ? 
 *                     parallel::distributed::Triangulation<dim>::default_setting : 
 *                     parallel::distributed::Triangulation< 
 *                       dim>::construct_multigrid_hierarchy) 
 *   , mapping() 
 *   , fe(degree) 
 *   , dof_handler(triangulation) 
 *   , computing_timer(pcout, TimerOutput::never, TimerOutput::wall_times) 
 * { 
 *   GridGenerator::hyper_L(triangulation, -1., 1., /*colorize*/ false); 
 *   triangulation.refine_global(1); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system()</h4>
 * 

 * 
 * 与  step-16  和  step-37  不同，我们将设置分成两部分，setup_system() 和 setup_multigrid() 。下面是大多数教程中常见的主动网格的典型setup_system()函数。对于无矩阵，活动网格的设置类似于  step-37  ；对于基于矩阵（GMG和AMG求解器），设置类似于  step-40  。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::setup_system() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Setup"); 
 * 
 *   dof_handler.distribute_dofs(fe); 
 * 
 *   DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
 *   locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 * 
 *   solution.reinit(locally_owned_dofs, mpi_communicator); 
 *   right_hand_side.reinit(locally_owned_dofs, mpi_communicator); 
 *   constraints.reinit(locally_relevant_dofs); 
 *   DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 * 
 *   VectorTools::interpolate_boundary_values( 
 *     mapping, dof_handler, 0, Functions::ZeroFunction<dim>(), constraints); 
 *   constraints.close(); 
 * 
 *   switch (settings.solver) 
 *     { 
 *       case Settings::gmg_mf: 
 *         { 
 *           typename MatrixFree<dim, double>::AdditionalData additional_data; 
 *           additional_data.tasks_parallel_scheme = 
 *             MatrixFree<dim, double>::AdditionalData::none; 
 *           additional_data.mapping_update_flags = 
 *             (update_gradients | update_JxW_values | update_quadrature_points); 
 *           std::shared_ptr<MatrixFree<dim, double>> mf_storage = 
 *             std::make_shared<MatrixFree<dim, double>>(); 
 *           mf_storage->reinit(mapping, 
 *                              dof_handler, 
 *                              constraints, 
 *                              QGauss<1>(degree + 1), 
 *                              additional_data); 
 * 
 *           mf_system_matrix.initialize(mf_storage); 
 * 
 *           const Coefficient<dim> coefficient; 
 *           mf_system_matrix.set_coefficient( 
 *             coefficient.make_coefficient_table(*mf_storage)); 
 * 
 *           break; 
 *         } 
 * 
 *       case Settings::gmg_mb: 
 *       case Settings::amg: 
 *         { 
 * #ifdef USE_PETSC_LA 
 *           DynamicSparsityPattern dsp(locally_relevant_dofs); 
 *           DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
 * 
 *           SparsityTools::distribute_sparsity_pattern(dsp, 
 *                                                      locally_owned_dofs, 
 *                                                      mpi_communicator, 
 *                                                      locally_relevant_dofs); 
 * 
 *           system_matrix.reinit(locally_owned_dofs, 
 *                                locally_owned_dofs, 
 *                                dsp, 
 *                                mpi_communicator); 
 * #else 
 *           TrilinosWrappers::SparsityPattern dsp(locally_owned_dofs, 
 *                                                 locally_owned_dofs, 
 *                                                 locally_relevant_dofs, 
 *                                                 mpi_communicator); 
 *           DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
 *           dsp.compress(); 
 *           system_matrix.reinit(dsp); 
 * #endif 
 * 
 *           break; 
 *         } 
 * 
 *       default: 
 *         Assert(false, ExcNotImplemented()); 
 *     } 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_multigrid"></a> 
 * <h4>LaplaceProblem::setup_multigrid()</h4>
 * 

 * 
 * 该函数为无矩阵和基于矩阵的GMG进行多级设置。无矩阵的设置类似于 step-37 ，而基于矩阵的设置类似于 step-16 ，只是我们必须使用适当的分布式稀疏度模式。
 * 

 * 
 * 该函数没有被AMG方法调用，但为了安全起见，该函数的主`switch`语句还是确保了该函数只在已知的多网格设置下运行，如果该函数被调用到两种几何多网格方法以外的地方，则抛出一个断言。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::setup_multigrid() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Setup multigrid"); 
 * 
 *   dof_handler.distribute_mg_dofs(); 
 * 
 *   mg_constrained_dofs.clear(); 
 *   mg_constrained_dofs.initialize(dof_handler); 
 * 
 *   const std::set<types::boundary_id> boundary_ids = {types::boundary_id(0)}; 
 *   mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, boundary_ids); 
 * 
 *   const unsigned int n_levels = triangulation.n_global_levels(); 
 * 
 *   switch (settings.solver) 
 *     { 
 *       case Settings::gmg_mf: 
 *         { 
 *           mf_mg_matrix.resize(0, n_levels - 1); 
 * 
 *           for (unsigned int level = 0; level < n_levels; ++level) 
 *             { 
 *               IndexSet relevant_dofs; 
 *               DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
 *                                                             level, 
 *                                                             relevant_dofs); 
 *               AffineConstraints<double> level_constraints; 
 *               level_constraints.reinit(relevant_dofs); 
 *               level_constraints.add_lines( 
 *                 mg_constrained_dofs.get_boundary_indices(level)); 
 *               level_constraints.close(); 
 * 
 *               typename MatrixFree<dim, float>::AdditionalData additional_data; 
 *               additional_data.tasks_parallel_scheme = 
 *                 MatrixFree<dim, float>::AdditionalData::none; 
 *               additional_data.mapping_update_flags = 
 *                 (update_gradients | update_JxW_values | 
 *                  update_quadrature_points); 
 *               additional_data.mg_level = level; 
 *               std::shared_ptr<MatrixFree<dim, float>> mf_storage_level( 
 *                 new MatrixFree<dim, float>()); 
 *               mf_storage_level->reinit(mapping, 
 *                                        dof_handler, 
 *                                        level_constraints, 
 *                                        QGauss<1>(degree + 1), 
 *                                        additional_data); 
 * 
 *               mf_mg_matrix[level].initialize(mf_storage_level, 
 *                                              mg_constrained_dofs, 
 *                                              level); 
 * 
 *               const Coefficient<dim> coefficient; 
 *               mf_mg_matrix[level].set_coefficient( 
 *                 coefficient.make_coefficient_table(*mf_storage_level)); 
 * 
 *               mf_mg_matrix[level].compute_diagonal(); 
 *             } 
 * 
 *           break; 
 *         } 
 * 
 *       case Settings::gmg_mb: 
 *         { 
 *           mg_matrix.resize(0, n_levels - 1); 
 *           mg_matrix.clear_elements(); 
 *           mg_interface_in.resize(0, n_levels - 1); 
 *           mg_interface_in.clear_elements(); 
 * 
 *           for (unsigned int level = 0; level < n_levels; ++level) 
 *             { 
 *               IndexSet dof_set; 
 *               DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
 *                                                             level, 
 *                                                             dof_set); 
 * 
 *               { 
 * #ifdef USE_PETSC_LA 
 *                 DynamicSparsityPattern dsp(dof_set); 
 *                 MGTools::make_sparsity_pattern(dof_handler, dsp, level); 
 *                 dsp.compress(); 
 *                 SparsityTools::distribute_sparsity_pattern( 
 *                   dsp, 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   mpi_communicator, 
 *                   dof_set); 
 * 
 *                 mg_matrix[level].reinit( 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dsp, 
 *                   mpi_communicator); 
 * #else 
 *                 TrilinosWrappers::SparsityPattern dsp( 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_set, 
 *                   mpi_communicator); 
 *                 MGTools::make_sparsity_pattern(dof_handler, dsp, level); 
 * 
 *                 dsp.compress(); 
 *                 mg_matrix[level].reinit(dsp); 
 * #endif 
 *               } 
 * 
 *               { 
 * #ifdef USE_PETSC_LA 
 *                 DynamicSparsityPattern dsp(dof_set); 
 *                 MGTools::make_interface_sparsity_pattern(dof_handler, 
 *                                                          mg_constrained_dofs, 
 *                                                          dsp, 
 *                                                          level); 
 *                 dsp.compress(); 
 *                 SparsityTools::distribute_sparsity_pattern( 
 *                   dsp, 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   mpi_communicator, 
 *                   dof_set); 
 * 
 *                 mg_interface_in[level].reinit( 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dsp, 
 *                   mpi_communicator); 
 * #else 
 *                 TrilinosWrappers::SparsityPattern dsp( 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_handler.locally_owned_mg_dofs(level), 
 *                   dof_set, 
 *                   mpi_communicator); 
 * 
 *                 MGTools::make_interface_sparsity_pattern(dof_handler, 
 *                                                          mg_constrained_dofs, 
 *                                                          dsp, 
 *                                                          level); 
 *                 dsp.compress(); 
 *                 mg_interface_in[level].reinit(dsp); 
 * #endif 
 *               } 
 *             } 
 *           break; 
 *         } 
 * 
 *       default: 
 *         Assert(false, ExcNotImplemented()); 
 *     } 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system()</h4>
 * 

 * 
 * 汇编被分成三个部分：`assemble_system()`, `assemble_multigrid()`, 和`assemble_rhs()`。这里的`assemble_system()`函数组装并存储（全局）系统矩阵和基于矩阵的方法的右手边。它类似于  step-40  中的装配。
 * 

 * 
 * 注意，无矩阵方法不执行这个函数，因为它不需要组装矩阵，而是在assemble_rhs()中组装右手边。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::assemble_system() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Assemble"); 
 * 
 *   const QGauss<dim> quadrature_formula(degree + 1); 
 * 
 *   FEValues<dim> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_gradients | 
 *                             update_quadrature_points | update_JxW_values); 
 * 
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *   const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *   Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *   const Coefficient<dim> coefficient; 
 *   RightHandSide<dim>     rhs; 
 *   std::vector<double>    rhs_values(n_q_points); 
 * 
 *   for (const auto &cell : dof_handler.active_cell_iterators()) 
 *     if (cell->is_locally_owned()) 
 *       { 
 *         cell_matrix = 0; 
 *         cell_rhs    = 0; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         const double coefficient_value = 
 *           coefficient.average_value(fe_values.get_quadrature_points()); 
 *         rhs.value_list(fe_values.get_quadrature_points(), rhs_values); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 cell_matrix(i, j) += 
 *                   coefficient_value *                // epsilon(x) 
 *                   fe_values.shape_grad(i, q_point) * // * grad phi_i(x) 
 *                   fe_values.shape_grad(j, q_point) * // * grad phi_j(x) 
 *                   fe_values.JxW(q_point);            // * dx 
 * 
 *               cell_rhs(i) += 
 *                 fe_values.shape_value(i, q_point) * // grad phi_i(x) 
 *                 rhs_values[q_point] *               // * f(x) 
 *                 fe_values.JxW(q_point);             // * dx 
 *             } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global(cell_matrix, 
 *                                                cell_rhs, 
 *                                                local_dof_indices, 
 *                                                system_matrix, 
 *                                                right_hand_side); 
 *       } 
 * 
 *   system_matrix.compress(VectorOperation::add); 
 *   right_hand_side.compress(VectorOperation::add); 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_multigrid"></a> 
 * <h4>LaplaceProblem::assemble_multigrid()</h4>
 * 

 * 
 * 下面的函数为基于矩阵的GMG方法组装和存储多级矩阵。这个函数与 step-16 中的函数类似，只是在这里它适用于分布式网格。这个区别在于增加了一个条件，即我们只在本地拥有的水平单元上进行组装，并为每个被建立的矩阵调用压缩（）。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::assemble_multigrid() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Assemble multigrid"); 
 * 
 *   QGauss<dim> quadrature_formula(degree + 1); 
 * 
 *   FEValues<dim> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_gradients | 
 *                             update_quadrature_points | update_JxW_values); 
 * 
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *   const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *   const Coefficient<dim> coefficient; 
 * 
 *   std::vector<AffineConstraints<double>> boundary_constraints( 
 *     triangulation.n_global_levels()); 
 *   for (unsigned int level = 0; level < triangulation.n_global_levels(); ++level) 
 *     { 
 *       IndexSet dof_set; 
 *       DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
 *                                                     level, 
 *                                                     dof_set); 
 *       boundary_constraints[level].reinit(dof_set); 
 *       boundary_constraints[level].add_lines( 
 *         mg_constrained_dofs.get_refinement_edge_indices(level)); 
 *       boundary_constraints[level].add_lines( 
 *         mg_constrained_dofs.get_boundary_indices(level)); 
 * 
 *       boundary_constraints[level].close(); 
 *     } 
 * 
 *   for (const auto &cell : dof_handler.cell_iterators()) 
 *     if (cell->level_subdomain_id() == triangulation.locally_owned_subdomain()) 
 *       { 
 *         cell_matrix = 0; 
 *         fe_values.reinit(cell); 
 * 
 *         const double coefficient_value = 
 *           coefficient.average_value(fe_values.get_quadrature_points()); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               cell_matrix(i, j) += 
 *                 coefficient_value * fe_values.shape_grad(i, q_point) * 
 *                 fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point); 
 * 
 *         cell->get_mg_dof_indices(local_dof_indices); 
 * 
 *         boundary_constraints[cell->level()].distribute_local_to_global( 
 *           cell_matrix, local_dof_indices, mg_matrix[cell->level()]); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             if (mg_constrained_dofs.is_interface_matrix_entry( 
 *                   cell->level(), local_dof_indices[i], local_dof_indices[j])) 
 *               mg_interface_in[cell->level()].add(local_dof_indices[i], 
 *                                                  local_dof_indices[j], 
 *                                                  cell_matrix(i, j)); 
 *       } 
 * 
 *   for (unsigned int i = 0; i < triangulation.n_global_levels(); ++i) 
 *     { 
 *       mg_matrix[i].compress(VectorOperation::add); 
 *       mg_interface_in[i].compress(VectorOperation::add); 
 *     } 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_rhs"></a> 
 * <h4>LaplaceProblem::assemble_rhs()</h4>
 * 

 * 
 * 这个三要素中的最后一个函数为无矩阵方法组装右手边的向量--因为在无矩阵框架中，我们不需要组装矩阵，只需要组装右手边就可以了。我们可以通过从上面的`assemble_system()`函数中提取处理右手边的代码来做到这一点，但是我们决定完全采用无矩阵的方法，也用这种方法进行装配。
 * 

 * 
 * 结果是一个类似于 step-37 中 "使用 FEEvaluation::read_dof_values_plain() 来避免解决约束 "一节中的函数。
 * 

 * 
 * 这个函数的原因是MatrixFree运算符不考虑非同质的Dirichlet约束，而是将所有的Dirichlet约束视为同质的。为了说明这一点，这里的右手边被组装成残差 $r_0 = f-Au_0$ ，其中 $u_0$ 是一个零向量，除了在Dirichlet值中。然后在求解的时候，我们可以看到，解决方案是  $u = u_0 + A^{-1}r_0$  。这可以看作是对初始猜测为  $u_0$  的线性系统进行的牛顿迭代。下面`solve()`函数中的CG解计算了 $A^{-1}r_0$ ，调用`constraints.distribution()`（直接在后面）增加了 $u_0$  。
 * 

 * 
 * 显然，由于我们考虑的是一个零迪里希特边界的问题，我们可以采取类似于 step-37  `assemble_rhs()`的方法，但是这个额外的工作允许我们改变问题声明，如果我们选择的话。
 * 

 * 
 * 这个函数在积分循环中有两个部分：通过提交梯度的负值将矩阵  $A$  的负值应用于  $u_0$  ，并通过提交值  $f$  添加右手边的贡献。我们必须确保使用`read_dof_values_plain()`来评估 $u_0$ ，因为`read_dof_vaues()`会将所有Dirichlet值设置为0。
 * 

 * 
 * 最后，system_rhs向量的类型是 LA::MPI::Vector, ，但MatrixFree类只对 dealii::LinearAlgebra::distributed::Vector. 起作用，因此我们必须使用MatrixFree功能计算右手边，然后使用`ChangeVectorType`命名空间的函数将其复制到正确的类型。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::assemble_rhs() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Assemble right-hand side"); 
 * 
 *   MatrixFreeActiveVector solution_copy; 
 *   MatrixFreeActiveVector right_hand_side_copy; 
 *   mf_system_matrix.initialize_dof_vector(solution_copy); 
 *   mf_system_matrix.initialize_dof_vector(right_hand_side_copy); 
 * 
 *   solution_copy = 0.; 
 *   constraints.distribute(solution_copy); 
 *   solution_copy.update_ghost_values(); 
 *   right_hand_side_copy = 0; 
 *   const Table<2, VectorizedArray<double>> &coefficient = 
 *     *(mf_system_matrix.get_coefficient()); 
 * 
 *   RightHandSide<dim> right_hand_side_function; 
 * 
 *   FEEvaluation<dim, degree, degree + 1, 1, double> phi( 
 *     *mf_system_matrix.get_matrix_free()); 
 * 
 *   for (unsigned int cell = 0; 
 *        cell < mf_system_matrix.get_matrix_free()->n_cell_batches(); 
 *        ++cell) 
 *     { 
 *       phi.reinit(cell); 
 *       phi.read_dof_values_plain(solution_copy); 
 *       phi.evaluate(EvaluationFlags::gradients); 
 * 
 *       for (unsigned int q = 0; q < phi.n_q_points; ++q) 
 *         { 
 *           phi.submit_gradient(-1.0 * 
 *                                 (coefficient(cell, 0) * phi.get_gradient(q)), 
 *                               q); 
 *           phi.submit_value( 
 *             right_hand_side_function.value(phi.quadrature_point(q)), q); 
 *         } 
 * 
 *       phi.integrate_scatter(EvaluationFlags::values | 
 *                               EvaluationFlags::gradients, 
 *                             right_hand_side_copy); 
 *     } 
 * 
 *   right_hand_side_copy.compress(VectorOperation::add); 
 * 
 *   ChangeVectorTypes::copy(right_hand_side, right_hand_side_copy); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve()</h4>
 * 

 * 
 * 这里我们设置了多网格预处理程序，测试了单个V型周期的时间，并解决了线性系统。不出所料，这是三种方法差别最大的地方之一。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::solve() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Solve"); 
 * 
 *   SolverControl solver_control(1000, 1.e-10 * right_hand_side.l2_norm()); 
 *   solver_control.enable_history_data(); 
 * 
 *   solution = 0.; 
 * 
 * @endcode
 * 
 * 无矩阵GMG方法的求解器类似于  step-37  ，除了增加一些接口矩阵，完全类似于  step-16  。
 * 

 * 
 * 
 * @code
 *   switch (settings.solver) 
 *     { 
 *       case Settings::gmg_mf: 
 *         { 
 *           computing_timer.enter_subsection("Solve: Preconditioner setup"); 
 * 
 *           MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs); 
 *           mg_transfer.build(dof_handler); 
 * 
 *           SolverControl coarse_solver_control(1000, 1e-12, false, false); 
 *           SolverCG<MatrixFreeLevelVector> coarse_solver(coarse_solver_control); 
 *           PreconditionIdentity            identity; 
 *           MGCoarseGridIterativeSolver<MatrixFreeLevelVector, 
 *                                       SolverCG<MatrixFreeLevelVector>, 
 *                                       MatrixFreeLevelMatrix, 
 *                                       PreconditionIdentity> 
 *             coarse_grid_solver(coarse_solver, mf_mg_matrix[0], identity); 
 * 
 *           using Smoother = dealii::PreconditionJacobi<MatrixFreeLevelMatrix>; 
 *           MGSmootherPrecondition<MatrixFreeLevelMatrix, 
 *                                  Smoother, 
 *                                  MatrixFreeLevelVector> 
 *             smoother; 
 *           smoother.initialize(mf_mg_matrix, 
 *                               typename Smoother::AdditionalData( 
 *                                 settings.smoother_dampen)); 
 *           smoother.set_steps(settings.smoother_steps); 
 * 
 *           mg::Matrix<MatrixFreeLevelVector> mg_m(mf_mg_matrix); 
 * 
 *           MGLevelObject< 
 *             MatrixFreeOperators::MGInterfaceOperator<MatrixFreeLevelMatrix>> 
 *             mg_interface_matrices; 
 *           mg_interface_matrices.resize(0, triangulation.n_global_levels() - 1); 
 *           for (unsigned int level = 0; level < triangulation.n_global_levels(); 
 *                ++level) 
 *             mg_interface_matrices[level].initialize(mf_mg_matrix[level]); 
 *           mg::Matrix<MatrixFreeLevelVector> mg_interface(mg_interface_matrices); 
 * 
 *           Multigrid<MatrixFreeLevelVector> mg( 
 *             mg_m, coarse_grid_solver, mg_transfer, smoother, smoother); 
 *           mg.set_edge_matrices(mg_interface, mg_interface); 
 * 
 *           PreconditionMG<dim, 
 *                          MatrixFreeLevelVector, 
 *                          MGTransferMatrixFree<dim, float>> 
 *             preconditioner(dof_handler, mg, mg_transfer); 
 * 
 * @endcode
 * 
 * 将求解向量和右手边从 LA::MPI::Vector 复制到 dealii::LinearAlgebra::distributed::Vector ，这样我们就可以解决了。
 * 

 * 
 * 
 * @code
 *           MatrixFreeActiveVector solution_copy; 
 *           MatrixFreeActiveVector right_hand_side_copy; 
 *           mf_system_matrix.initialize_dof_vector(solution_copy); 
 *           mf_system_matrix.initialize_dof_vector(right_hand_side_copy); 
 * 
 *           ChangeVectorTypes::copy(solution_copy, solution); 
 *           ChangeVectorTypes::copy(right_hand_side_copy, right_hand_side); 
 *           computing_timer.leave_subsection("Solve: Preconditioner setup"); 
 * 
 * @endcode
 * 
 * 1个V型周期的时间安排。
 * 

 * 
 * 
 * @code
 *           { 
 *             TimerOutput::Scope timing(computing_timer, 
 *                                       "Solve: 1 multigrid V-cycle"); 
 *             preconditioner.vmult(solution_copy, right_hand_side_copy); 
 *           } 
 *           solution_copy = 0.; 
 * 
 * @endcode
 * 
 * 解出线性系统，更新解的鬼魂值，复制回 LA::MPI::Vector 并分配约束。
 * 

 * 
 * 
 * @code
 *           { 
 *             SolverCG<MatrixFreeActiveVector> solver(solver_control); 
 * 
 *             TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
 *             solver.solve(mf_system_matrix, 
 *                          solution_copy, 
 *                          right_hand_side_copy, 
 *                          preconditioner); 
 *           } 
 * 
 *           solution_copy.update_ghost_values(); 
 *           ChangeVectorTypes::copy(solution, solution_copy); 
 *           constraints.distribute(solution); 
 * 
 *           break; 
 *         } 
 * 
 * @endcode
 * 
 * 基于矩阵的GMG方法的求解器，类似于  step-16  ，只是使用了雅可比平滑器，而不是SOR平滑器（该平滑器没有并行实现）。
 * 

 * 
 * 
 * @code
 *       case Settings::gmg_mb: 
 *         { 
 *           computing_timer.enter_subsection("Solve: Preconditioner setup"); 
 * 
 *           MGTransferPrebuilt<VectorType> mg_transfer(mg_constrained_dofs); 
 *           mg_transfer.build(dof_handler); 
 * 
 *           SolverControl        coarse_solver_control(1000, 1e-12, false, false); 
 *           SolverCG<VectorType> coarse_solver(coarse_solver_control); 
 *           PreconditionIdentity identity; 
 *           MGCoarseGridIterativeSolver<VectorType, 
 *                                       SolverCG<VectorType>, 
 *                                       MatrixType, 
 *                                       PreconditionIdentity> 
 *             coarse_grid_solver(coarse_solver, mg_matrix[0], identity); 
 * 
 *           using Smoother = LA::MPI::PreconditionJacobi; 
 *           MGSmootherPrecondition<MatrixType, Smoother, VectorType> smoother; 
 * 
 * #ifdef USE_PETSC_LA 
 *           smoother.initialize(mg_matrix); 
 *           Assert( 
 *             settings.smoother_dampen == 1.0, 
 *             ExcNotImplemented( 
 *               "PETSc's PreconditionJacobi has no support for a damping parameter.")); 
 * #else 
 *           smoother.initialize(mg_matrix, settings.smoother_dampen); 
 * #endif 
 * 
 *           smoother.set_steps(settings.smoother_steps); 
 * 
 *           mg::Matrix<VectorType> mg_m(mg_matrix); 
 *           mg::Matrix<VectorType> mg_in(mg_interface_in); 
 *           mg::Matrix<VectorType> mg_out(mg_interface_in); 
 * 
 *           Multigrid<VectorType> mg( 
 *             mg_m, coarse_grid_solver, mg_transfer, smoother, smoother); 
 *           mg.set_edge_matrices(mg_out, mg_in); 
 * 
 *           PreconditionMG<dim, VectorType, MGTransferPrebuilt<VectorType>> 
 *             preconditioner(dof_handler, mg, mg_transfer); 
 * 
 *           computing_timer.leave_subsection("Solve: Preconditioner setup"); 
 * 
 * @endcode
 * 
 * 1个V型周期的计时。
 * 

 * 
 * 
 * @code
 *           { 
 *             TimerOutput::Scope timing(computing_timer, 
 *                                       "Solve: 1 multigrid V-cycle"); 
 *             preconditioner.vmult(solution, right_hand_side); 
 *           } 
 *           solution = 0.; 
 * 
 * @endcode
 * 
 * 解决线性系统和分配约束。
 * 

 * 
 * 
 * @code
 *           { 
 *             SolverCG<VectorType> solver(solver_control); 
 * 
 *             TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
 *             solver.solve(system_matrix, 
 *                          solution, 
 *                          right_hand_side, 
 *                          preconditioner); 
 *           } 
 * 
 *           constraints.distribute(solution); 
 * 
 *  
 *         } 
 * 
 * @endcode
 * 
 * AMG方法的求解器，类似于  step-40  。
 * 

 * 
 * 
 * @code
 *       case Settings::amg: 
 *         { 
 *           computing_timer.enter_subsection("Solve: Preconditioner setup"); 
 * 
 *           PreconditionAMG                 preconditioner; 
 *           PreconditionAMG::AdditionalData Amg_data; 
 * 
 * #ifdef USE_PETSC_LA 
 *           Amg_data.symmetric_operator = true; 
 * #else 
 *           Amg_data.elliptic              = true; 
 *           Amg_data.smoother_type         = "Jacobi"; 
 *           Amg_data.higher_order_elements = true; 
 *           Amg_data.smoother_sweeps       = settings.smoother_steps; 
 *           Amg_data.aggregation_threshold = 0.02; 
 * #endif 
 * 
 *           Amg_data.output_details = false; 
 * 
 *           preconditioner.initialize(system_matrix, Amg_data); 
 *           computing_timer.leave_subsection("Solve: Preconditioner setup"); 
 * 
 * @endcode
 * 
 * 1个V型周期的计时。
 * 

 * 
 * 
 * @code
 *           { 
 *             TimerOutput::Scope timing(computing_timer, 
 *                                       "Solve: 1 multigrid V-cycle"); 
 *             preconditioner.vmult(solution, right_hand_side); 
 *           } 
 *           solution = 0.; 
 * 
 * @endcode
 * 
 * 解决线性系统和分配约束。
 * 

 * 
 * 
 * @code
 *           { 
 *             SolverCG<VectorType> solver(solver_control); 
 * 
 *             TimerOutput::Scope timing(computing_timer, "Solve: CG"); 
 *             solver.solve(system_matrix, 
 *                          solution, 
 *                          right_hand_side, 
 *                          preconditioner); 
 *           } 
 *           constraints.distribute(solution); 
 * 
 *           break; 
 *         } 
 * 
 *       default: 
 *         Assert(false, ExcInternalError()); 
 *     } 
 * 
 *   pcout << "   Number of CG iterations:      " << solver_control.last_step() 
 *         << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="Theerrorestimator"></a> 
 * <h3>The error estimator</h3>
 * 

 * 
 * 我们使用FEInterfaceValues类来组装一个误差估计器，以决定哪些单元需要细化。请看介绍中对单元和面积分的确切定义。为了使用该方法，我们为 MeshWorker::mesh_loop() 定义了Scratch和Copy对象，下面的大部分代码本质上与 step-12 中已经设置的一样（或者至少精神上相似）。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * struct ScratchData 
 * { 
 *   ScratchData(const Mapping<dim> &      mapping, 
 *               const FiniteElement<dim> &fe, 
 *               const unsigned int        quadrature_degree, 
 *               const UpdateFlags         update_flags, 
 *               const UpdateFlags         interface_update_flags) 
 *     : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
 *     , fe_interface_values(mapping, 
 *                           fe, 
 *                           QGauss<dim - 1>(quadrature_degree), 
 *                           interface_update_flags) 
 *   {} 
 * 
 *   ScratchData(const ScratchData<dim> &scratch_data) 
 *     : fe_values(scratch_data.fe_values.get_mapping(), 
 *                 scratch_data.fe_values.get_fe(), 
 *                 scratch_data.fe_values.get_quadrature(), 
 *                 scratch_data.fe_values.get_update_flags()) 
 *     , fe_interface_values(scratch_data.fe_values.get_mapping(), 
 *                           scratch_data.fe_values.get_fe(), 
 *                           scratch_data.fe_interface_values.get_quadrature(), 
 *                           scratch_data.fe_interface_values.get_update_flags()) 
 *   {} 
 * 
 *   FEValues<dim>          fe_values; 
 *   FEInterfaceValues<dim> fe_interface_values; 
 * }; 
 * 
 * struct CopyData 
 * { 
 *   CopyData() 
 *     : cell_index(numbers::invalid_unsigned_int) 
 *     , value(0.) 
 *   {} 
 * 
 *   CopyData(const CopyData &) = default; 
 * 
 *   struct FaceData 
 *   { 
 *     unsigned int cell_indices[2]; 
 *     double       values[2]; 
 *   }; 
 * 
 *   unsigned int          cell_index; 
 *   double                value; 
 *   std::vector<FaceData> face_data; 
 * }; 
 * 
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::estimate() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Estimate"); 
 * 
 *   VectorType temp_solution; 
 *   temp_solution.reinit(locally_owned_dofs, 
 *                        locally_relevant_dofs, 
 *                        mpi_communicator); 
 *   temp_solution = solution; 
 * 
 *   const Coefficient<dim> coefficient; 
 * 
 *   estimated_error_square_per_cell.reinit(triangulation.n_active_cells()); 
 * 
 *   using Iterator = typename DoFHandler<dim>::active_cell_iterator; 
 * 
 * @endcode
 * 
 * 剩余单元的汇编程序  $h^2 \| f + \epsilon \triangle u \|_K^2$  。
 * 
 * @code
 *   auto cell_worker = [&](const Iterator &  cell, 
 *                          ScratchData<dim> &scratch_data, 
 *                          CopyData &        copy_data) { 
 *     FEValues<dim> &fe_values = scratch_data.fe_values; 
 *     fe_values.reinit(cell); 
 * 
 *     RightHandSide<dim> rhs; 
 *     const double       rhs_value = rhs.value(cell->center()); 
 * 
 *     const double nu = coefficient.value(cell->center()); 
 * 
 *  
 *     fe_values.get_function_hessians(temp_solution, hessians); 
 * 
 *     copy_data.cell_index = cell->active_cell_index(); 
 * 
 *     double residual_norm_square = 0.; 
 *     for (unsigned k = 0; k < fe_values.n_quadrature_points; ++k) 
 *       { 
 *         const double residual = (rhs_value + nu * trace(hessians[k])); 
 *         residual_norm_square += residual * residual * fe_values.JxW(k); 
 *       } 
 * 
 *     copy_data.value = 
 *       cell->diameter() * cell->diameter() * residual_norm_square; 
 *   }; 
 * 
 * @endcode
 * 
 * 脸部术语的汇编器  $\sum_F h_F \| \jump{\epsilon \nabla u \cdot n} \|_F^2$  。
 * 
 * @code
 *   auto face_worker = [&](const Iterator &    cell, 
 *                          const unsigned int &f, 
 *                          const unsigned int &sf, 
 *                          const Iterator &    ncell, 
 *                          const unsigned int &nf, 
 *                          const unsigned int &nsf, 
 *                          ScratchData<dim> &  scratch_data, 
 *                          CopyData &          copy_data) { 
 *     FEInterfaceValues<dim> &fe_interface_values = 
 *       scratch_data.fe_interface_values; 
 *     fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf); 
 * 
 *     copy_data.face_data.emplace_back(); 
 *     CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 
 * 
 *     copy_data_face.cell_indices[0] = cell->active_cell_index(); 
 *     copy_data_face.cell_indices[1] = ncell->active_cell_index(); 
 * 
 *     const double coeff1 = coefficient.value(cell->center()); 
 *     const double coeff2 = coefficient.value(ncell->center()); 
 * 
 *     std::vector<Tensor<1, dim>> grad_u[2]; 
 * 
 *     for (unsigned int i = 0; i < 2; ++i) 
 *       { 
 *         grad_u[i].resize(fe_interface_values.n_quadrature_points); 
 *         fe_interface_values.get_fe_face_values(i).get_function_gradients( 
 *           temp_solution, grad_u[i]); 
 *       } 
 * 
 *     double jump_norm_square = 0.; 
 * 
 *     for (unsigned int qpoint = 0; 
 *          qpoint < fe_interface_values.n_quadrature_points; 
 *          ++qpoint) 
 *       { 
 *         const double jump = 
 *           coeff1 * grad_u[0][qpoint] * fe_interface_values.normal(qpoint) - 
 *           coeff2 * grad_u[1][qpoint] * fe_interface_values.normal(qpoint); 
 * 
 *         jump_norm_square += jump * jump * fe_interface_values.JxW(qpoint); 
 *       } 
 * 
 *     const double h           = cell->face(f)->measure(); 
 *     copy_data_face.values[0] = 0.5 * h * jump_norm_square; 
 *     copy_data_face.values[1] = copy_data_face.values[0]; 
 *   }; 
 * 
 *   auto copier = [&](const CopyData &copy_data) { 
 *     if (copy_data.cell_index != numbers::invalid_unsigned_int) 
 *       estimated_error_square_per_cell[copy_data.cell_index] += copy_data.value; 
 * 
 *     for (auto &cdf : copy_data.face_data) 
 *       for (unsigned int j = 0; j < 2; ++j) 
 *         estimated_error_square_per_cell[cdf.cell_indices[j]] += cdf.values[j]; 
 *   }; 
 * 
 *   const unsigned int n_gauss_points = degree + 1; 
 *   ScratchData<dim>   scratch_data(mapping, 
 *                                 fe, 
 *                                 n_gauss_points, 
 *                                 update_hessians | update_quadrature_points | 
 *                                   update_JxW_values, 
 *                                 update_values | update_gradients | 
 *                                   update_JxW_values | update_normal_vectors); 
 *   CopyData           copy_data; 
 * 
 * @endcode
 * 
 * 我们需要对每个内部面进行一次装配，但我们需要确保两个进程都对本地拥有的单元和幽灵单元之间的面术语进行装配。这可以通过设置 MeshWorker::assemble_ghost_faces_both 标志来实现。我们需要这样做，因为我们不在这里交流误差估计器的贡献。
 * 

 * 
 * 
 * @code
 *   MeshWorker::mesh_loop(dof_handler.begin_active(), 
 *                         dof_handler.end(), 
 *                         cell_worker, 
 *                         copier, 
 *                         scratch_data, 
 *                         copy_data, 
 *                         MeshWorker::assemble_own_cells | 
 *                           MeshWorker::assemble_ghost_faces_both | 
 *                           MeshWorker::assemble_own_interior_faces_once, 
 *                           /*boundary_worker=*/nullptr, face_worker);
 * 
 * 
 *  
 * 
 *   const double global_error_estimate = 
 *     std::sqrt(Utilities::MPI::sum(estimated_error_square_per_cell.l1_norm(), 
 *                                   mpi_communicator)); 
 *   pcout << "   Global error estimate:        " << global_error_estimate 
 *         << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrefine_grid"></a> 
 * <h4>LaplaceProblem::refine_grid()</h4>
 * 

 * 
 * 我们使用存储在向量 @p estimate_vector 中的单元估计器，并细化固定数量的单元（这里选择的是每一步中大约两倍的DoFs数量）。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::refine_grid() 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Refine grid"); 
 * 
 *   const double refinement_fraction = 1. / (std::pow(2.0, dim) - 1.); 
 *   parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
 *     triangulation, estimated_error_square_per_cell, refinement_fraction, 0.0); 
 * 
 *   triangulation.execute_coarsening_and_refinement(); 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemoutput_results"></a> 
 * <h4>LaplaceProblem::output_results()</h4>
 * 

 * 
 * output_results()函数与许多教程中的函数类似（例如，见 step-40 ）。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::output_results(const unsigned int cycle) 
 * { 
 *   TimerOutput::Scope timing(computing_timer, "Output results"); 
 * 
 *   VectorType temp_solution; 
 *   temp_solution.reinit(locally_owned_dofs, 
 *                        locally_relevant_dofs, 
 *                        mpi_communicator); 
 *   temp_solution = solution; 
 * 
 *   DataOut<dim> data_out; 
 *   data_out.attach_dof_handler(dof_handler); 
 *   data_out.add_data_vector(temp_solution, "solution"); 
 * 
 *   Vector<float> subdomain(triangulation.n_active_cells()); 
 *   for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *     subdomain(i) = triangulation.locally_owned_subdomain(); 
 *   data_out.add_data_vector(subdomain, "subdomain"); 
 * 
 *   Vector<float> level(triangulation.n_active_cells()); 
 *   for (const auto &cell : triangulation.active_cell_iterators()) 
 *     level(cell->active_cell_index()) = cell->level(); 
 *   data_out.add_data_vector(level, "level"); 
 * 
 *   if (estimated_error_square_per_cell.size() > 0) 
 *     data_out.add_data_vector(estimated_error_square_per_cell, 
 *                              "estimated_error_square_per_cell"); 
 * 
 *   data_out.build_patches(); 
 * 
 *   const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
 *     "", "solution", cycle, mpi_communicator, 2 /*n_digits*/, 1 /*n_groups*/); 
 * 
 *   pcout << "   Wrote " << pvtu_filename << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run()</h4>
 * 

 * 
 * 和大多数教程一样，这个函数调用上面定义的各种函数来设置、组合、求解和输出结果。
 * 

 * 
 * 
 * @code
 * template <int dim, int degree> 
 * void LaplaceProblem<dim, degree>::run() 
 * { 
 *   for (unsigned int cycle = 0; cycle < settings.n_steps; ++cycle) 
 *     { 
 *       pcout << "Cycle " << cycle << ':' << std::endl; 
 *       if (cycle > 0) 
 *         refine_grid(); 
 * 
 *       pcout << "   Number of active cells:       " 
 *             << triangulation.n_global_active_cells(); 
 * 
 * @endcode
 * 
 * 我们只为GMG方法输出层次单元数据（与下面的DoF数据相同）。请注意，对于AMG来说，分区效率是不相关的，因为在计算过程中没有分布或使用层次结构。
 * 

 * 
 * 
 * @code
 *       if (settings.solver == Settings::gmg_mf || 
 *           settings.solver == Settings::gmg_mb) 
 *         pcout << " (" << triangulation.n_global_levels() << " global levels)" 
 *               << std::endl 
 *               << "   Partition efficiency:         " 
 *               << 1.0 / MGTools::workload_imbalance(triangulation); 
 *       pcout << std::endl; 
 * 
 *       setup_system(); 
 * 
 * @endcode
 * 
 * 只为GMG设置多级层次结构。
 * 

 * 
 * 
 * @code
 *       if (settings.solver == Settings::gmg_mf || 
 *           settings.solver == Settings::gmg_mb) 
 *         setup_multigrid(); 
 * 
 *       pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs(); 
 *       if (settings.solver == Settings::gmg_mf || 
 *           settings.solver == Settings::gmg_mb) 
 *         { 
 *           pcout << " (by level: "; 
 *           for (unsigned int level = 0; level < triangulation.n_global_levels(); 
 *                ++level) 
 *             pcout << dof_handler.n_dofs(level) 
 *                   << (level == triangulation.n_global_levels() - 1 ? ")" : 
 *                                                                      ", "); 
 *         } 
 *       pcout << std::endl; 
 * 
 * @endcode
 * 
 * 对于无矩阵的方法，我们只组装右手边。对于这两种基于矩阵的方法，我们同时装配主动矩阵和右手边，对于基于矩阵的GMG，我们只装配多网格矩阵。
 * 

 * 
 * 
 * @code
 *       if (settings.solver == Settings::gmg_mf) 
 *         assemble_rhs(); 
 *       else /*gmg_mb or amg*/ 
 *         { 
 *           assemble_system(); 
 *           if (settings.solver == Settings::gmg_mb) 
 *             assemble_multigrid(); 
 *         } 
 * 
 *       solve(); 
 *       estimate(); 
 * 
 *       if (settings.output) 
 *         output_results(cycle); 
 * 
 *       computing_timer.print_summary(); 
 *       computing_timer.reset(); 
 *     } 
 * } 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 这是一个类似于 step-40 的主函数，但我们要求用户传递一个.prm文件作为唯一的命令行参数（参见 step-29 和ParameterHandler类的文档，以了解关于参数文件的完整讨论）。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   using namespace dealii; 
 *   Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *   Settings settings; 
 *   if (!settings.try_parse((argc > 1) ? (argv[1]) : "")) 
 *     return 0; 
 * 
 *   try 
 *     { 
 *       constexpr unsigned int fe_degree = 2; 
 * 
 *       switch (settings.dimension) 
 *         { 
 *           case 2: 
 *             { 
 *               LaplaceProblem<2, fe_degree> test(settings); 
 *               test.run(); 
 * 
 *               break; 
 *             } 
 * 
 *           case 3: 
 *             { 
 *               LaplaceProblem<3, fe_degree> test(settings); 
 *               test.run(); 
 * 
 *               break; 
 *             } 
 * 
 *           default: 
 *             Assert(false, ExcMessage("This program only works in 2d and 3d.")); 
 *         } 
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
 *       MPI_Abort(MPI_COMM_WORLD, 1); 
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
 *       MPI_Abort(MPI_COMM_WORLD, 2); 
 *       return 1; 
 *     } 
 * 
 *   return 0; 
 * } 
 * 
 * @endcode
examples/step-50/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你使用以下命令运行该程序时

@code
mpirun -np 16 ./step-50  gmg_mf_2d.prm
@endcode

屏幕输出应该如下。

@code
Cycle 0:
   Number of active cells:       12 (2 global levels)
   Partition efficiency:         0.1875
   Number of degrees of freedom: 65 (by level: 21, 65)
   Number of CG iterations:      10
   Global error estimate:        0.355373
   Wrote solution_00.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |    0.0163s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble right-hand side        |         1 |  0.000374s |       2.3% |
| Estimate                        |         1 |  0.000724s |       4.4% |
| Output results                  |         1 |   0.00277s |        17% |
| Setup                           |         1 |   0.00225s |        14% |
| Setup multigrid                 |         1 |   0.00181s |        11% |
| Solve                           |         1 |   0.00364s |        22% |
| Solve: 1 multigrid V-cycle      |         1 |  0.000354s |       2.2% |
| Solve: CG                       |         1 |   0.00151s |       9.3% |
| Solve: Preconditioner setup     |         1 |   0.00125s |       7.7% |
+---------------------------------+-----------+------------+------------+


Cycle 1:
   Number of active cells:       24 (3 global levels)
   Partition efficiency:         0.276786
   Number of degrees of freedom: 139 (by level: 21, 65, 99)
   Number of CG iterations:      10
   Global error estimate:        0.216726
   Wrote solution_01.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |    0.0169s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble right-hand side        |         1 |  0.000309s |       1.8% |
| Estimate                        |         1 |   0.00156s |       9.2% |
| Output results                  |         1 |   0.00222s |        13% |
| Refine grid                     |         1 |   0.00278s |        16% |
| Setup                           |         1 |   0.00196s |        12% |
| Setup multigrid                 |         1 |    0.0023s |        14% |
| Solve                           |         1 |   0.00565s |        33% |
| Solve: 1 multigrid V-cycle      |         1 |  0.000349s |       2.1% |
| Solve: CG                       |         1 |   0.00285s |        17% |
| Solve: Preconditioner setup     |         1 |   0.00195s |        12% |
+---------------------------------+-----------+------------+------------+


Cycle 2:
   Number of active cells:       51 (4 global levels)
   Partition efficiency:         0.41875
   Number of degrees of freedom: 245 (by level: 21, 65, 225, 25)
   Number of CG iterations:      11
   Global error estimate:        0.112098
   Wrote solution_02.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |    0.0183s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble right-hand side        |         1 |  0.000274s |       1.5% |
| Estimate                        |         1 |   0.00127s |       6.9% |
| Output results                  |         1 |   0.00227s |        12% |
| Refine grid                     |         1 |    0.0024s |        13% |
| Setup                           |         1 |   0.00191s |        10% |
| Setup multigrid                 |         1 |   0.00295s |        16% |
| Solve                           |         1 |   0.00702s |        38% |
| Solve: 1 multigrid V-cycle      |         1 |  0.000398s |       2.2% |
| Solve: CG                       |         1 |   0.00376s |        21% |
| Solve: Preconditioner setup     |         1 |   0.00238s |        13% |
+---------------------------------+-----------+------------+------------+
.
.
.
@endcode

在这里，"solve() "函数的时间被分成三部分：设置多网格预处理程序，执行单一的多网格V型循环，以及CG求解器。被计时的V型循环对整个求解来说是不必要的，只是为了让人们了解AMG和GMG的不同成本。还应注意的是，当使用AMG求解器时，"工作量不平衡 "不包括在输出中，因为不需要粗网格的层次结构。

本节中的所有结果都是在英特尔至强铂金8280（Cascade Lake）节点上收集的，这些节点有56个内核，每个节点有192GB，支持AVX-512指令，允许对8个双倍数进行矢量化（矢量化仅用于无矩阵计算）。代码是用gcc 7.1.0和intel-mpi 17.0.3编译的。Trilinos 12.10.1被用于基于矩阵的GMG/AMG计算。

然后我们可以通过调用step-50目录下的输入文件来收集各种信息。使用这些文件，并调整网格细化步骤的数量，我们可以得出程序的扩展性如何的信息。

下表给出了该程序在高达256M自由度和7168个处理器上的弱比例计时。(回顾一下，在增加处理器数量的同时，弱缩放保持每个处理器的自由度数量不变；也就是说，它考虑的是越来越大的问题。)这里， $\mathbb{E}$ 是介绍中的分区效率（也等于1.0/工作量不平衡），"Setup "是设置、设置多重网格、装配和装配多重网格的组合，"Prec "是预处理程序的设置。理想情况下，在每个问题大小上，各个求解器的所有时间都保持不变，但由于分区效率从最大问题大小到最小问题大小从0.371下降到0.161，我们期望看到GMG的时间大约增加 $0.371/0.161=2.3$ 倍。事实上，这与我们实际得到的情况非常接近。

 <table align="center" class="doxtable">
<tr>
  <th colspan="4"></th>
  <th></th>
  <th colspan="4">MF-GMG</th>
  <th></th>
  <th colspan="4">MB-GMG</th>
  <th></th>
  <th colspan="4">AMG</th>
</tr>
<tr>
  <th align="right">Procs</th>
  <th align="right">Cycle</th>
  <th align="right">DoFs</th>
  <th align="right">$\mathbb{E}$</th>
  <th></th>
  <th align="right">Setup</th>
  <th align="right">Prec</th>
  <th align="right">Solve</th>
  <th align="right">Total</th>
  <th></th>
  <th align="right">Setup</th>
  <th align="right">Prec</th>
  <th align="right">Solve</th>
  <th align="right">Total</th>
  <th></th>
  <th align="right">Setup</th>
  <th align="right">Prec</th>
  <th align="right">Solve</th>
  <th align="right">Total</th>
</tr>
<tr>
  <td align="right">112</th>
  <td align="right">13</th>
  <td align="right">4M</th>
  <td align="right">0.37</th>
  <td></td>
  <td align="right">0.742</th>
  <td align="right">0.393</th>
  <td align="right">0.200</th>
  <td align="right">1.335</th>
  <td></td>
  <td align="right">1.714</th>
  <td align="right">2.934</th>
  <td align="right">0.716</th>
  <td align="right">5.364</th>
  <td></td>
  <td align="right">1.544</th>
  <td align="right">0.456</th>
  <td align="right">1.150</th>
  <td align="right">3.150</th>
</tr>
<tr>
  <td align="right">448</th>
  <td align="right">15</th>
  <td align="right">16M</th>
  <td align="right">0.29</th>
  <td></td>
  <td align="right">0.884</th>
  <td align="right">0.535</th>
  <td align="right">0.253</th>
  <td align="right">1.672</th>
  <td></td>
  <td align="right">1.927</th>
  <td align="right">3.776</th>
  <td align="right">1.190</th>
  <td align="right">6.893</th>
  <td></td>
  <td align="right">1.544</th>
  <td align="right">0.456</th>
  <td align="right">1.150</th>
  <td align="right">3.150</th>
</tr>
<tr>
  <td align="right">1,792</th>
  <td align="right">17</th>
  <td align="right">65M</th>
  <td align="right">0.22</th>
  <td></td>
  <td align="right">1.122</th>
  <td align="right">0.686</th>
  <td align="right">0.309</th>
  <td align="right">2.117</th>
  <td></td>
  <td align="right">2.171</th>
  <td align="right">4.862</th>
  <td align="right">1.660</th>
  <td align="right">8.693</th>
  <td></td>
  <td align="right">1.654</th>
  <td align="right">0.546</th>
  <td align="right">1.460</th>
  <td align="right">3.660</th>
</tr>
<tr>
  <td align="right">7,168</th>
  <td align="right">19</th>
  <td align="right">256M</th>
  <td align="right">0.16</th>
  <td></td>
  <td align="right">1.214</th>
  <td align="right">0.893</th>
  <td align="right">0.521</th>
  <td align="right">2.628</th>
  <td></td>
  <td align="right">2.386</th>
  <td align="right">7.260</th>
  <td align="right">2.560</th>
  <td align="right">12.206</th>
  <td></td>
  <td align="right">1.844</th>
  <td align="right">1.010</th>
  <td align="right">1.890</th>
  <td align="right">4.744</th>
</tr>
</table> 

另一方面，最后一组列中的代数多网格相对来说不受网格层次不平衡度增加的影响（因为它不使用网格层次），时间的增长反而是由文献中记载的其他因素驱动的（最明显的是，代数多网格方法的某些部分的算法复杂度似乎是 ${\cal O}(N
\log N)$ ，而不是几何多网格的 ${\cal O}(N)$ ）。

上表的短处是，无矩阵的几何多网格方法似乎是解决这个方程的最快方法，即使不是以很大的优势。另一方面，基于矩阵的方法始终是最差的。

下图提供了每种方法的强大扩展结果，也就是说，我们在越来越多的处理器上解决同一个问题。具体来说，我们考虑在56至28672个处理器上，经过16个网格细化周期（32M DoFs）和19个周期（256M DoFs）后的问题。

 <img width="600px" src="https://www.dealii.org/images/steps/developer/step-50-strong-scaling.png" alt=""> 

虽然基于矩阵的GMG求解器和AMG的规模相似，求解时间也相似（至少在每个处理器有大量未知数的情况下--比如说，几万个），但无矩阵的GMG求解器的规模要好得多，在较粗的网格上，只用八分之一的处理器就能解决较细的问题，与AMG求解器的时间大致相当。反之，在相同数量的处理器上，它可以用大约八分之一的时间解决同样的问题。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="Testingconvergenceandhigherorderelements"></a><h4> Testing convergence and higher order elements </h4>


有限元度目前是硬编码为2，见主类的模板参数。这很容易改变。为了测试，最好是切换到一个有参考解的测试问题。这样，你可以比较错误率。

<a name="Coarsesolver"></a><h4> Coarse solver </h4>


一个更有趣的例子是涉及到一个更复杂的粗大的网格（见步骤49的启发）。这种情况下的问题是，网格层次的最粗层实际上是相当大的，我们必须考虑如何有效地解决粗层问题。这对代数多网格方法来说不是一个问题，因为它们只是继续建立越来越粗的矩阵层次，而不管它们的几何来源）。

在这里的程序中，我们只是简单地用共轭梯度法解决粗级问题，没有任何预处理程序。如果粗略问题真的很小，这是可以接受的--例如，如果粗略网格有一个单元，那么粗略网格问题在2d中有一个 $9\times 9$ 矩阵，在3d中有一个 $27\times 27$ 矩阵；对于我们在当前程序的 $L$ 形域上使用的粗略网格，这些大小在2d中是 $21\times 21$ ，在3d中有 $117\times 117$ 。但如果粗略的网格由数百或数千个单元组成，这种方法将不再起作用，并可能开始主导每个V型单元的整体运行时间。一个常见的方法是使用代数多网格预处理程序来解决粗网格问题；然而，这将需要组装粗矩阵（即使是无矩阵版本）作为AMG实现的输入。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-50.cc"
*/
