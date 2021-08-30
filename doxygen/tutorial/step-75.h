/**
@page step_75 The step-75 tutorial program
This tutorial depends on step-27, step-37, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Loadbalancing">Load balancing</a>
        <li><a href="#hpdecisionindicators">hp-decision indicators</a>
        <li><a href="#Hybridgeometricmultigrid">Hybrid geometric multigrid</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeSolutioncodeclasstemplate">The <code>Solution</code> class template</a>
        <li><a href="#Parameters">Parameters</a>
        <li><a href="#MatrixfreeLaplaceoperator">Matrix-free Laplace operator</a>
        <li><a href="#Solverandpreconditioner">Solver and preconditioner</a>
      <ul>
        <li><a href="#Conjugategradientsolverwithmultigridpreconditioner">Conjugate-gradient solver with multigrid preconditioner</a>
        <li><a href="#Hybridpolynomialgeometricglobalcoarseningmultigridpreconditioner">Hybrid polynomial/geometric-global-coarsening multigrid preconditioner</a>
      </ul>
        <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The <code>LaplaceProblem</code> class template</a>
        <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The <code>LaplaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#LaplaceProbleminitialize_grid">LaplaceProblem::initialize_grid</a>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemprint_diagnostics">LaplaceProblem::print_diagnostics</a>
        <li><a href="#LaplaceProblemsolve_system">LaplaceProblem::solve_system</a>
        <li><a href="#LaplaceProblemcompute_indicators">LaplaceProblem::compute_indicators</a>
        <li><a href="#LaplaceProblemadapt_resolution">LaplaceProblem::adapt_resolution</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
        <li><a href="#main">main()</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Differenthpdecisionstrategies">Different hp-decision strategies</a>
        <li><a href="#Solvewithmatrixbasedmethods">Solve with matrix-based methods</a>
        <li><a href="#Multigridvariants">Multigrid variants</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-75/doc/intro.dox

 <br> 

<i>This program was contributed by Marc Fehling, Peter Munch and
Wolfgang Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. DMS-1821210, EAR-1550901, and
OAC-1835673. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the authors and do not
necessarily reflect the views of the National Science Foundation.
<br>
Peter Munch would like to thank Timo Heister, Martin Kronbichler, and
Laura Prieto Saavedra for many very interesting discussions.
</i>




 @note  作为这个程序的先决条件，你需要安装p4est库和Trilinos库。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II与这些附加库的安装情况。




<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


在有限元背景下，更多的自由度通常会产生一个更精确的解决方案，但也需要更多的计算工作。

在以前的整个教程中，我们找到了通过将网格分辨率与解的复杂性进行局部调整来有效分配自由度的方法（自适应网格细化，步骤6）。如果我们不仅单独调整网格，而且还局部调整每个单元上相关有限元的多项式程度，这种方法就特别有效（hp-adaptation，第27步）。

此外，分配更多的进程同时运行你的程序有助于在更短的时间内解决计算工作量。根据你的机器的硬件结构，你的程序必须为所有进程都能访问相同的内存（共享内存，第18步），或者进程被托管在几个独立的节点上（分布式内存，第40步）这种情况做好准备。

在高性能计算部分，内存访问变成了当前超级计算机的瓶颈。我们可以通过使用MatrixFree方法（第37步）来计算矩阵-向量乘积的效果，从而完全避免存储矩阵。它们可以用于几何多网格方法（步骤50），也可以用于多项式多网格方法，以极大地加快方程组的求解速度。

本教程结合所有这些特点，介绍了如何解决一个简单的拉普拉斯问题的最先进的方法：在具有分布式内存的机器上利用hp-适应和无矩阵混合多网格方法。




<a name="Loadbalancing"></a><h3>Load balancing</h3>


对于有限元的并行应用，我们将网格划分为子域（又称域分解），这些子域被分配给进程。这种划分发生在deal.II的活动单元上，如步骤40所示。在那里，每个单元都有相同的有限元和相同的自由度分配，以及大致相同的工作负荷。为了平衡所有进程的工作负荷，我们必须平衡所有参与进程上的单元数量。

在hp-adaptive方法中，情况不再如此：有限元类型可能因单元而异，因此自由度的数量也不同。匹配单元的数量并不能产生一个平衡的工作量。在无矩阵的情况下，可以假设工作量与每个过程的自由度数量成正比，因为在最好的情况下，只有源和目的向量需要被加载。

我们可以通过给每个单元分配权重来平衡工作量，这些权重与自由度的数量成正比，并平衡所有进程之间的所有权重之和。给每个单元分配单独的权重可以通过我们后面要使用的 parallel::CellWeights 类来实现。




<a name="hpdecisionindicators"></a><h3>hp-decision indicators</h3>


使用hp-adaptive方法，我们不仅要决定哪些单元需要细化或粗化，而且还可以选择如何做：要么调整网格分辨率，要么调整有限元的多项式程度。

我们将再次根据当前解决方案的（后验）计算误差估计值来决定哪些单元需要调整，例如，使用KellyErrorEstimator。我们将同样决定如何用（事后）计算的平滑度估计值进行调整：大的多项式度数对解决方案的平滑部分效果最好，而细的网格分辨率对不规则部分是有利的。在第27步中，我们提出了一种基于傅里叶系数衰减的平滑度估计的计算方法。让我们利用这个机会，提出一种遵循相同思路的替代方法，但采用Legendre系数。

我们将简要介绍这种新技术的思路，但为了简单起见，将其描述限制在一维。假设 $u_\text{hp}(x)$ 是一个有限元函数，在单元格 $K$ 上定义为

@f[
u_\text{hp}(x) = \sum c_i \varphi_i(x)


@f]

其中每个 $\varphi_i(x)$ 是一个形状函数。我们可以用Legendre多项式 $P_k$ 的基础等价表示 $u_\text{hp}(x)$ 为

@f[
u_\text{hp}(x) = \sum l_k P_k(x).


@f]

我们的目标是获得有限元系数 $c_i$ 和Legendre系数 $l_k$ 之间的映射。我们将通过把问题写成 $L^2$ 对 $u_\text{hp}(x)$ 在Legendre基础上的投影来实现这一目标。每个系数 $l_k$ 可以通过以下方式计算

@f[
l_k = \int_K u_\text{hp}(x) P_k(x) dx.


@f]

根据结构，Legendre多项式在 $L^2$ 上的内积下是正交的。此外，我们假设它们已经被归一化，所以它们的内积可以写成

@f[
\int_K P_i(x) P_j(x) dx = \det(J_K) \, \delta_{ij}


@f]

其中 $\delta_{ij}$ 是克朗克三角洲， $J_K$ 是 $\hat{K}$ 到 $K$ 的映射的雅各布，（在本教程中）假定它是常数（即，映射必须是仿射的）。

因此，结合所有这些假设，在Legendre基础上表达 $u_\text{hp}(x)$ 的投影矩阵只是 $\det(J_K) \,
\mathbb{I}$  -- 即 $\det(J_K)$ 乘以身份矩阵。让 $F_K$ 成为从 $K$ 到其参考单元 $\hat{K}$ 的映射。因此，投影系统中右侧的条目为：。

@f[
\int_K u_\text{hp}(x) P_k(x) dx
= \det(J_K) \int_\hat{K} u_\text{hp}(F_K(\hat{x})) P_k(F_K(\hat{x})) d\hat{x}.


@f]

回顾 $u_\text{hp}(x)$ 的形状函数表示，我们可以把它写成 $\det(J_K) \, \mathbf{C} \, \mathbf{c}$ ，其中 $\mathbf{C}$ 是改变基础的矩阵，条目是

@f[
\int_K P_i(x) \varphi_j(x) dx
= \det(J_K) \int_{\hat{K}} P_i(F_K(\hat{x})) \varphi_j(F_K(\hat{x})) d\hat{x}
= \det(J_K) \int_{\hat{K}} \hat{P}_i(\hat{x}) \hat{\varphi}_j(\hat{x}) d\hat{x}
\dealcoloneq \det(J_K) \, C_{ij}


@f]

所以 $\mathbf{C}$ 的值可以写成 <em> 独立于 </em> 的 $K$ ，在转换为参考坐标后，将 $\det(J_K)$ 从前面因式分解。因此，把这一切放在一起，投影问题可以写为

@f[
\det(J_K) \, \mathbb{I} \, \mathbf{l} = \det(J_K) \, \mathbf{C} \, \mathbf{c}


@f]

可以简单改写为

@f[
\mathbf{l} = \mathbf{C} \, \mathbf{c}.


@f]



在这一点上，我们需要强调的是，大多数有限元应用都使用非结构化网格，对于这些网格的映射几乎总是非affine的。换句话说： $J_K$ 在整个单元中是恒定的这一假设对于一般的网格来说是不正确的。因此， $l_k$ 的正确计算不仅要求我们为每一个单元计算相应的变换矩阵 $\mathbf{C}$ ，而且还要求我们在可能具有任意和非常复杂的几何形状的单元 $K$ 上定义一组类Legendre正交函数。特别是第二部分，在计算上非常昂贵。目前FESeries变换类的实现依赖于具有恒定雅各布系数所带来的简化，以提高性能，因此只对仿射映射产生正确结果。变换只用于平滑度估计的目的，以决定适应的类型，这不是有限元程序的一个关键组成部分。除此之外，这种情况对本教程不构成问题，因为我们只使用方形的单元。

Eibner和Melenk  @cite eibner2007hp  认为，当且仅当Legendre系数的绝对值随指数增加而衰减时，一个函数是解析的，即可以用幂级数表示  $k$  。

@f[
\exists C,\sigma > 0 : \quad \forall k \in \mathbb{N}_0 : \quad |l_k|
\leq C \exp\left( - \sigma k \right) .


@f]

衰减率 $\sigma$ 可以被解释为衡量该函数的平滑度。我们可以把它看成是转化系数的线性回归拟合的斜率。

@f[
\ln(|l_k|) \sim \ln(C) - \sigma k .


@f]



我们将对每个单元 $K$ 进行这种拟合，以获得对有限元近似的平滑度的局部估计。然后，衰减率 $\sigma_K$ 作为hp-adaptation的决策指标。对于单元上的有限元 $K$ 的多项式程度 $p$ ，计算 $k \leq (p+1)$ 的系数被证明是估计平稳性的合理选择。你可以在  @cite fehling2020  中找到更详细和独立于维度的描述。

以上所有内容已经在 FESeries::Legendre 类和 SmoothnessEstimator::Legendre 命名空间中实现。有了误差估计和平滑度指标，我们就可以对单元格进行实际细化和粗化了。来自 parallel::distributed::GridRefinement 和 hp::Refinement 命名空间的一些函数将在后面帮助我们完成这个任务。




<a name="Hybridgeometricmultigrid"></a><h3>Hybrid geometric multigrid</h3>


有限元矩阵通常是非常稀疏的。此外，hp-adaptive方法对应于每行非零项数量变化很大的矩阵。一些最先进的预处理程序，如Step-40中使用的代数多重网格（AMG），在这些情况下表现不佳。

因此，我们将依靠一个无矩阵的混合多网格预处理程序。Step-50已经证明了几何多网格方法与MatrixFree框架结合时的优越性。在hp-adaptive FEM上的应用需要一些额外的工作，因为一个单元的子代可能有不同的多项式程度。作为补救措施，我们首先对线性元素进行p松弛（类似于Mitchell @cite mitchell2010hpmg ），然后以常规方式进行h松弛。在最粗的层次上，我们应用代数多网格求解器。p-多栅、h-多栅和AMG的结合使求解器成为一个混合多栅求解器。

我们将通过使用MGTransferGlobalCoarsening，在现有的全局粗化基础设施的帮助下，创建一个具有上述特殊水平要求的自定义混合多网格预处理器。




<a name="Thetestcase"></a><h3>The test case</h3>


对于椭圆方程来说，每个再入角通常会引出一个奇点  @cite brenner2008  。我们可以利用这种情况对我们的HP决策算法进行测试：在所有要适应的单元上，我们倾向于在奇点附近采用精细的网格，而在其他情况下采用高的多项式程度。

作为在这些条件下要解决的最简单的椭圆问题，我们选择了L型域中的拉普拉斯方程，其再入角位于坐标系的原点。

为了能够确定实际的误差，我们制造一个有已知解的边界值问题。在上述领域，拉普拉斯方程的一个解是，在极坐标中， $(r, \varphi)$  。

@f[
u_\text{sol} = r^{2/3} \sin(2/3 \varphi).


@f]



参见  @cite brenner2008  或  @cite mitchell2014hp  。解决方案看起来如下。

<div style="text-align:center;"> <img src="https://www.dealii.org/images/steps/developer/step-75.solution.svg" alt="分析性解决方案。"> </div>

通过研究再入角附近的解决方案的梯度，即原点，奇异性变得很明显了。

@f[
\left\| \nabla u_\text{sol} \right\|_{2} = 2/3 r^{-1/3} , \quad
\lim\limits_{r \rightarrow 0} \left\| \nabla u_\text{sol} \right\|_{2} =
\infty .


@f]



由于我们知道奇点的位置，我们希望我们的hp-decision算法在这个特定的区域内决定采用精细的网格分辨率，而在其他地方采用高多项式程度。

因此，让我们看看情况是否真的如此，以及hp-adaptation与纯h-adaptation相比表现如何。但首先让我们详细看看实际的代码。


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
 * 在以前的教程程序中，特别是在 step-27 和 step-40 中，已经使用和讨论了以下包含文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/base/mpi.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/timer.h> 
 * 
 * #include <deal.II/distributed/grid_refinement.h> 
 * #include <deal.II/distributed/tria.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_series.h> 
 * 
 * #include <deal.II/hp/fe_collection.h> 
 * #include <deal.II/hp/refinement.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/vector.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/smoothness_estimator.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * #include <algorithm> 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 为了实现负载平衡，我们将在单元格上分配单独的权重，为此我们将使用类  parallel::CellWeights.  。
 * 
 * @code
 * #include <deal.II/distributed/cell_weights.h> 
 * 
 * @endcode
 * 
 * 求解函数需要从直角坐标到极坐标的转换。 GeometricUtilities::Coordinates 命名空间提供了必要的工具。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/geometric_utilities.h> 
 * 
 * @endcode
 * 
 * 以下包含的文件将启用MatrixFree功能。
 * 

 * 
 * 
 * @code
 * #include <deal.II/matrix_free/matrix_free.h> 
 * #include <deal.II/matrix_free/fe_evaluation.h> 
 * #include <deal.II/matrix_free/tools.h> 
 * 
 * @endcode
 * 
 * 我们将使用 LinearAlgebra::distributed::Vector 进行线性代数操作。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/la_parallel_vector.h> 
 * 
 * @endcode
 * 
 * 我们剩下的就是包含多网格求解器所需的文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/multigrid/mg_coarse.h> 
 * #include <deal.II/multigrid/mg_constrained_dofs.h> 
 * #include <deal.II/multigrid/mg_matrix.h> 
 * #include <deal.II/multigrid/mg_smoother.h> 
 * #include <deal.II/multigrid/mg_tools.h> 
 * #include <deal.II/multigrid/mg_transfer_global_coarsening.h> 
 * #include <deal.II/multigrid/multigrid.h> 
 * 
 * namespace Step75 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeSolutioncodeclasstemplate"></a> 
 * <h3>The <code>Solution</code> class template</h3>
 * 

 * 
 * 我们有一个分析性的方案可以使用。我们将用这个解来为问题的数值解施加边界条件。解决方案的表述需要转换为极坐标。为了从笛卡尔坐标转换到球面坐标，我们将使用 GeometricUtilities::Coordinates 命名空间的一个辅助函数。这个转换的前两个坐标对应于x-y面的极坐标。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solution : public Function<dim> 
 *   { 
 *   public: 
 *     Solution() 
 *       : Function<dim>() 
 *     {} 
 * 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/) const override 
 *     { 
 *       const std::array<double, dim> p_sphere = 
 *         GeometricUtilities::Coordinates::to_spherical(p); 
 * 
 *       constexpr const double alpha = 2. / 3.; 
 *       return std::pow(p_sphere[0], alpha) * std::sin(alpha * p_sphere[1]); 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="Parameters"></a> 
 * <h3>Parameters</h3>
 * 

 * 
 * 在本教程中，我们将使用一个简化的参数集。这里也可以使用ParameterHandler类，但为了使本教程简短，我们决定使用简单的结构。所有这些参数的实际意图将在接下来的类中描述，在它们各自使用的位置。
 * 

 * 
 * 下面的参数集控制着多网格机制的粗网格求解器、平滑器和网格间传输方案。我们用默认参数来填充它。
 * 

 * 
 * 
 * @code
 *   struct MultigridParameters 
 *   { 
 *     struct 
 *     { 
 *       std::string  type            = "cg_with_amg"; 
 *       unsigned int maxiter         = 10000; 
 *       double       abstol          = 1e-20; 
 *       double       reltol          = 1e-4; 
 *       unsigned int smoother_sweeps = 1; 
 *       unsigned int n_cycles        = 1; 
 *       std::string  smoother_type   = "ILU"; 
 *     } coarse_solver; 
 * 
 *     struct 
 *     { 
 *       std::string  type                = "chebyshev"; 
 *       double       smoothing_range     = 20; 
 *       unsigned int degree              = 5; 
 *       unsigned int eig_cg_n_iterations = 20; 
 *     } smoother; 
 * 
 *     struct 
 *     { 
 *       MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType 
 *         p_sequence = MGTransferGlobalCoarseningTools:: 
 *           PolynomialCoarseningSequenceType::decrease_by_one; 
 *       bool perform_h_transfer = true; 
 *     } transfer; 
 *   }; 
 * 
 * @endcode
 * 
 * 这是该问题类的一般参数结构。你会发现这个结构分为几个类别，包括一般的运行时参数、级别限制、细化和粗化分数，以及单元加权的参数。它还包含一个上述结构的实例，用于多网格参数，这些参数将被传递给多网格算法。
 * 

 * 
 * 
 * @code
 *   struct Parameters 
 *   { 
 *     unsigned int n_cycles         = 8; 
 *     double       tolerance_factor = 1e-12; 
 * 
 *     MultigridParameters mg_data; 
 * 
 *     unsigned int min_h_level            = 5; 
 *     unsigned int max_h_level            = 12; 
 *     unsigned int min_p_degree           = 2; 
 *     unsigned int max_p_degree           = 6; 
 *     unsigned int max_p_level_difference = 1; 
 * 
 *     double refine_fraction    = 0.3; 
 *     double coarsen_fraction   = 0.03; 
 *     double p_refine_fraction  = 0.9; 
 *     double p_coarsen_fraction = 0.9; 
 * 
 *     double weighting_factor   = 1e6; 
 *     double weighting_exponent = 1.; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="MatrixfreeLaplaceoperator"></a> 
 * <h3>Matrix-free Laplace operator</h3>
 * 

 * 
 * 这是一个无矩阵的拉普拉斯算子的实现，基本上将接管其他教程中的`assemble_system()`函数的部分。所有成员函数的含义将在后面的定义中解释。
 * 

 * 
 * 我们将使用FEEvaluation类来评估正交点的解向量并进行积分。与其他教程不同的是，模板参数`度数`被设置为  $-1$  ，`一维正交数`被设置为  $0$  。在这种情况下，FEEvaluation会动态地选择正确的多项式度数和正交点的数量。在这里，我们为FEEvaluation引入一个带有正确模板参数的别名，这样我们以后就不用担心这些参数了。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   class LaplaceOperator : public Subscriptor 
 *   { 
 *   public: 
 *     using VectorType = LinearAlgebra::distributed::Vector<number>; 
 * 
 *     using FECellIntegrator = FEEvaluation<dim, -1, 0, 1, number>; 
 * 
 *     LaplaceOperator() = default; 
 * 
 *     LaplaceOperator(const hp::MappingCollection<dim> &mapping, 
 *                     const DoFHandler<dim> &           dof_handler, 
 *                     const hp::QCollection<dim> &      quad, 
 *                     const AffineConstraints<number> & constraints, 
 *                     VectorType &                      system_rhs); 
 * 
 *     void reinit(const hp::MappingCollection<dim> &mapping, 
 *                 const DoFHandler<dim> &           dof_handler, 
 *                 const hp::QCollection<dim> &      quad, 
 *                 const AffineConstraints<number> & constraints, 
 *                 VectorType &                      system_rhs); 
 * 
 *     types::global_dof_index m() const; 
 * 
 *     number el(unsigned int, unsigned int) const; 
 * 
 *     void initialize_dof_vector(VectorType &vec) const; 
 * 
 *     void vmult(VectorType &dst, const VectorType &src) const; 
 * 
 *     void Tvmult(VectorType &dst, const VectorType &src) const; 
 * 
 *     const TrilinosWrappers::SparseMatrix &get_system_matrix() const; 
 * 
 *     void compute_inverse_diagonal(VectorType &diagonal) const; 
 * 
 *   private: 
 *     void do_cell_integral_local(FECellIntegrator &integrator) const; 
 * 
 *     void do_cell_integral_global(FECellIntegrator &integrator, 
 *                                  VectorType &      dst, 
 *                                  const VectorType &src) const; 
 * 
 *     void do_cell_integral_range( 
 *       const MatrixFree<dim, number> &              matrix_free, 
 *       VectorType &                                 dst, 
 *       const VectorType &                           src, 
 *       const std::pair<unsigned int, unsigned int> &range) const; 
 * 
 *     MatrixFree<dim, number> matrix_free; 
 * 
 * @endcode
 * 
 * 为了用AMG预处理程序解决最粗层次的方程系统，我们需要一个最粗层次的实际系统矩阵。为此，我们提供了一种机制，可以选择从无矩阵公式中计算出一个矩阵，为此我们引入了一个专门的SparseMatrix对象。在默认情况下，这个矩阵保持为空。一旦`get_system_matrix()`被调用，这个矩阵就会被填充（懒惰分配）。由于这是一个 "const "函数，我们需要在这里使用 "mutable "关键字。我们还需要一个约束对象来构建矩阵。
 * 

 * 
 * 
 * @code
 *     AffineConstraints<number>              constraints; 
 *     mutable TrilinosWrappers::SparseMatrix system_matrix; 
 *   }; 
 * 
 * @endcode
 * 
 * 下面的部分包含了初始化和重新初始化该类的函数。特别是，这些函数初始化了内部的MatrixFree实例。为了简单起见，我们还计算了系统右侧的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   LaplaceOperator<dim, number>::LaplaceOperator( 
 *     const hp::MappingCollection<dim> &mapping, 
 *     const DoFHandler<dim> &           dof_handler, 
 *     const hp::QCollection<dim> &      quad, 
 *     const AffineConstraints<number> & constraints, 
 *     VectorType &                      system_rhs) 
 *   { 
 *     this->reinit(mapping, dof_handler, quad, constraints, system_rhs); 
 *   } 
 * 
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::reinit( 
 *     const hp::MappingCollection<dim> &mapping, 
 *     const DoFHandler<dim> &           dof_handler, 
 *     const hp::QCollection<dim> &      quad, 
 *     const AffineConstraints<number> & constraints, 
 *     VectorType &                      system_rhs) 
 *   { 
 * 
 * @endcode
 * 
 * 清除内部数据结构（在操作者被重复使用的情况下）。
 * 

 * 
 * 
 * @code
 *     this->system_matrix.clear(); 
 * 
 * @endcode
 * 
 * 复制约束条件，因为以后在计算系统矩阵时可能需要它们。
 * 

 * 
 * 
 * @code
 *     this->constraints.copy_from(constraints); 
 * 
 * @endcode
 * 
 * 设置MatrixFree。在正交点，我们只需要评估解的梯度，并用形状函数的梯度进行测试，所以我们只需要设置标志`update_gradients`。
 * 

 * 
 * 
 * @code
 *     typename MatrixFree<dim, number>::AdditionalData data; 
 *     data.mapping_update_flags = update_gradients; 
 * 
 *     matrix_free.reinit(mapping, dof_handler, constraints, quad, data); 
 * 
 * @endcode
 * 
 * 计算右手边的向量。为此，我们设置了第二个MatrixFree实例，它使用一个修改过的AffineConstraints，不包含由于Dirichlet-边界条件的约束。这个修改过的算子被应用于一个只设置了迪里希特值的向量。其结果是负的右手边向量。
 * 

 * 
 * 
 * @code
 *     { 
 *       AffineConstraints<number> constraints_without_dbc; 
 * 
 *       IndexSet locally_relevant_dofs; 
 *       DoFTools::extract_locally_relevant_dofs(dof_handler, 
 *                                               locally_relevant_dofs); 
 *       constraints_without_dbc.reinit(locally_relevant_dofs); 
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                               constraints_without_dbc); 
 *       constraints_without_dbc.close(); 
 * 
 *       VectorType b, x; 
 * 
 *       this->initialize_dof_vector(system_rhs); 
 * 
 *       MatrixFree<dim, number> matrix_free; 
 *       matrix_free.reinit( 
 *         mapping, dof_handler, constraints_without_dbc, quad, data); 
 * 
 *       matrix_free.initialize_dof_vector(b); 
 *       matrix_free.initialize_dof_vector(x); 
 * 
 *       constraints.distribute(x); 
 * 
 *       matrix_free.cell_loop(&LaplaceOperator::do_cell_integral_range, 
 *                             this, 
 *                             b, 
 *                             x); 
 * 
 *       constraints.set_zero(b); 
 * 
 *       system_rhs -= b; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 以下函数是多网格算法隐含需要的，包括平滑器。
 * 

 * 
 * 由于我们没有矩阵，所以要向DoFHandler查询自由度的数量。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   types::global_dof_index LaplaceOperator<dim, number>::m() const 
 *   { 
 *     return matrix_free.get_dof_handler().n_dofs(); 
 *   } 
 * 
 * @endcode
 * 
 * 访问矩阵中的一个特定元素。这个函数既不需要也没有实现，但是，在编译程序时需要它。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   number LaplaceOperator<dim, number>::el(unsigned int, unsigned int) const 
 *   { 
 *     Assert(false, ExcNotImplemented()); 
 *     return 0; 
 *   } 
 * 
 * @endcode
 * 
 * 初始化给定的向量。我们只是把这个任务委托给同名的MatrixFree函数。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void 
 *   LaplaceOperator<dim, number>::initialize_dof_vector(VectorType &vec) const 
 *   { 
 *     matrix_free.initialize_dof_vector(vec); 
 *   } 
 * 
 * @endcode
 * 
 * 在MatrixFree的帮助下，通过在所有单元中循环进行运算评估，并评估单元积分的效果（参见。`do_cell_integral_local()`和`do_cell_integral_global()`）。)
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::vmult(VectorType &      dst, 
 *                                            const VectorType &src) const 
 *   { 
 *     this->matrix_free.cell_loop( 
 *       &LaplaceOperator::do_cell_integral_range, this, dst, src, true); 
 *   } 
 * 
 * @endcode
 * 
 * 执行转置的运算符评估。由于我们考虑的是对称的 "矩阵"，这个函数可以简单地将其任务委托给vmult()。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::Tvmult(VectorType &      dst, 
 *                                             const VectorType &src) const 
 *   { 
 *     this->vmult(dst, src); 
 *   } 
 * 
 * @endcode
 * 
 * 由于我们没有一个系统矩阵，我们不能循环计算矩阵的对角线项。相反，我们通过对单位基向量进行一连串的运算符评估来计算对角线。为此，我们使用了MatrixFreeTools命名空间中的一个优化函数。之后再手动进行反转。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::compute_inverse_diagonal( 
 *     VectorType &diagonal) const 
 *   { 
 *     MatrixFreeTools::compute_diagonal(matrix_free, 
 *                                       diagonal, 
 *                                       &LaplaceOperator::do_cell_integral_local, 
 *                                       this); 
 * 
 *     for (auto &i : diagonal) 
 *       i = (std::abs(i) > 1.0e-10) ? (1.0 / i) : 1.0; 
 *   } 
 * 
 * @endcode
 * 
 * 在无矩阵的情况下，在这个类的初始化过程中没有设置系统矩阵。因此，如果需要的话，它必须在这里被计算出来。由于矩阵在本教程中只对线性元素进行计算（在粗略的网格上），这一点是可以接受的。矩阵的条目是通过运算符的评估序列得到的。为此，使用了优化函数 MatrixFreeTools::compute_matrix() 。矩阵只有在尚未设置的情况下才会被计算（懒惰分配）。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   const TrilinosWrappers::SparseMatrix & 
 *   LaplaceOperator<dim, number>::get_system_matrix() const 
 *   { 
 *     if (system_matrix.m() == 0 && system_matrix.n() == 0) 
 *       { 
 *         const auto &dof_handler = this->matrix_free.get_dof_handler(); 
 * 
 *         TrilinosWrappers::SparsityPattern dsp( 
 *           dof_handler.locally_owned_dofs(), 
 *           dof_handler.get_triangulation().get_communicator()); 
 * 
 *         DoFTools::make_sparsity_pattern(dof_handler, dsp, this->constraints); 
 * 
 *         dsp.compress(); 
 *         system_matrix.reinit(dsp); 
 * 
 *         MatrixFreeTools::compute_matrix( 
 *           matrix_free, 
 *           constraints, 
 *           system_matrix, 
 *           &LaplaceOperator::do_cell_integral_local, 
 *           this); 
 *       } 
 * 
 *     return this->system_matrix; 
 *   } 
 * 
 * @endcode
 * 
 * 对一个单元格批处理进行单元格积分，不需要收集和分散数值。MatrixFreeTools函数需要这个函数，因为这些函数直接对FEEvaluation的缓冲区进行操作。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::do_cell_integral_local( 
 *     FECellIntegrator &integrator) const 
 *   { 
 *     integrator.evaluate(EvaluationFlags::gradients); 
 * 
 *     for (unsigned int q = 0; q < integrator.n_q_points; ++q) 
 *       integrator.submit_gradient(integrator.get_gradient(q), q); 
 * 
 *     integrator.integrate(EvaluationFlags::gradients); 
 *   } 
 * 
 * @endcode
 * 
 * 与上述相同，但可以访问全局向量。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::do_cell_integral_global( 
 *     FECellIntegrator &integrator, 
 *     VectorType &      dst, 
 *     const VectorType &src) const 
 *   { 
 *     integrator.gather_evaluate(src, EvaluationFlags::gradients); 
 * 
 *     for (unsigned int q = 0; q < integrator.n_q_points; ++q) 
 *       integrator.submit_gradient(integrator.get_gradient(q), q); 
 * 
 *     integrator.integrate_scatter(EvaluationFlags::gradients, dst); 
 *   } 
 * 
 * @endcode
 * 
 * 这个函数在一个单元格批次范围内的所有单元格批次上循环，并调用上述函数。
 * 

 * 
 * 
 * @code
 *   template <int dim, typename number> 
 *   void LaplaceOperator<dim, number>::do_cell_integral_range( 
 *     const MatrixFree<dim, number> &              matrix_free, 
 *     VectorType &                                 dst, 
 *     const VectorType &                           src, 
 *     const std::pair<unsigned int, unsigned int> &range) const 
 *   { 
 *     FECellIntegrator integrator(matrix_free, range); 
 * 
 *     for (unsigned cell = range.first; cell < range.second; ++cell) 
 *       { 
 *         integrator.reinit(cell); 
 * 
 *         do_cell_integral_global(integrator, dst, src); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Solverandpreconditioner"></a> 
 * <h3>Solver and preconditioner</h3>
 * 
 * <a name="Conjugategradientsolverwithmultigridpreconditioner"></a> 
 * <h4>Conjugate-gradient solver with multigrid preconditioner</h4>
 * 

 * 
 * 这个函数用一连串提供的多网格对象来解决方程组。它的目的是为了尽可能的通用，因此有许多模板参数。
 * 

 * 
 * 
 * @code
 *   template <typename VectorType, 
 *             int dim, 
 *             typename SystemMatrixType, 
 *             typename LevelMatrixType, 
 *             typename MGTransferType> 
 *   static void 
 *   mg_solve(SolverControl &            solver_control, 
 *            VectorType &               dst, 
 *            const VectorType &         src, 
 *            const MultigridParameters &mg_data, 
 *            const DoFHandler<dim> &    dof, 
 *            const SystemMatrixType &   fine_matrix, 
 *            const MGLevelObject<std::unique_ptr<LevelMatrixType>> &mg_matrices, 
 *            const MGTransferType &                                 mg_transfer) 
 *   { 
 *     AssertThrow(mg_data.coarse_solver.type == "cg_with_amg", 
 *                 ExcNotImplemented()); 
 *     AssertThrow(mg_data.smoother.type == "chebyshev", ExcNotImplemented()); 
 * 
 *     const unsigned int min_level = mg_matrices.min_level(); 
 *     const unsigned int max_level = mg_matrices.max_level(); 
 * 
 *     using SmootherPreconditionerType = DiagonalMatrix<VectorType>; 
 *     using SmootherType               = PreconditionChebyshev<LevelMatrixType, 
 *                                                VectorType, 
 *                                                SmootherPreconditionerType>; 
 *     using PreconditionerType = PreconditionMG<dim, VectorType, MGTransferType>; 
 * 
 * @endcode
 * 
 * 我们在这里初始化电平运算符和切比雪夫平滑器。
 * 

 * 
 * 
 * @code
 *     mg::Matrix<VectorType> mg_matrix(mg_matrices); 
 * 
 *     MGLevelObject<typename SmootherType::AdditionalData> smoother_data( 
 *       min_level, max_level); 
 * 
 *     for (unsigned int level = min_level; level <= max_level; level++) 
 *       { 
 *         smoother_data[level].preconditioner = 
 *           std::make_shared<SmootherPreconditionerType>(); 
 *         mg_matrices[level]->compute_inverse_diagonal( 
 *           smoother_data[level].preconditioner->get_vector()); 
 *         smoother_data[level].smoothing_range = mg_data.smoother.smoothing_range; 
 *         smoother_data[level].degree          = mg_data.smoother.degree; 
 *         smoother_data[level].eig_cg_n_iterations = 
 *           mg_data.smoother.eig_cg_n_iterations; 
 *       } 
 * 
 *     MGSmootherPrecondition<LevelMatrixType, SmootherType, VectorType> 
 *       mg_smoother; 
 *     mg_smoother.initialize(mg_matrices, smoother_data); 
 * 
 * @endcode
 * 
 * 接下来，我们初始化粗略网格求解器。我们使用共轭梯度法和AMG作为预处理程序。
 * 

 * 
 * 
 * @code
 *     ReductionControl coarse_grid_solver_control(mg_data.coarse_solver.maxiter, 
 *                                                 mg_data.coarse_solver.abstol, 
 *                                                 mg_data.coarse_solver.reltol, 
 *                                                 false, 
 *                                                 false); 
 *     SolverCG<VectorType> coarse_grid_solver(coarse_grid_solver_control); 
 * 
 *     std::unique_ptr<MGCoarseGridBase<VectorType>> mg_coarse; 
 * 
 *     TrilinosWrappers::PreconditionAMG                 precondition_amg; 
 *     TrilinosWrappers::PreconditionAMG::AdditionalData amg_data; 
 *     amg_data.smoother_sweeps = mg_data.coarse_solver.smoother_sweeps; 
 *     amg_data.n_cycles        = mg_data.coarse_solver.n_cycles; 
 *     amg_data.smoother_type   = mg_data.coarse_solver.smoother_type.c_str(); 
 * 
 *     precondition_amg.initialize(mg_matrices[min_level]->get_system_matrix(), 
 *                                 amg_data); 
 * 
 *     mg_coarse = 
 *       std::make_unique<MGCoarseGridIterativeSolver<VectorType, 
 *                                                    SolverCG<VectorType>, 
 *                                                    LevelMatrixType, 
 *                                                    decltype(precondition_amg)>>( 
 *         coarse_grid_solver, *mg_matrices[min_level], precondition_amg); 
 * 
 * @endcode
 * 
 * 最后，我们创建Multigrid对象，将其转换为预处理程序，并在共轭梯度求解器中使用它来解决线性方程组。
 * 

 * 
 * 
 * @code
 *     Multigrid<VectorType> mg( 
 *       mg_matrix, *mg_coarse, mg_transfer, mg_smoother, mg_smoother); 
 * 
 *     PreconditionerType preconditioner(dof, mg, mg_transfer); 
 * 
 *     SolverCG<VectorType>(solver_control) 
 *       .solve(fine_matrix, dst, src, preconditioner); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Hybridpolynomialgeometricglobalcoarseningmultigridpreconditioner"></a> 
 * <h4>Hybrid polynomial/geometric-global-coarsening multigrid preconditioner</h4>
 * 

 * 
 * 上述函数处理给定的多网格对象序列的实际解决方案。这个函数创建了实际的多重网格层次，特别是运算符，以及作为MGTransferGlobalCoarsening对象的转移运算符。
 * 

 * 
 * 
 * @code
 *   template <typename VectorType, typename OperatorType, int dim> 
 *   void solve_with_gmg(SolverControl &                  solver_control, 
 *                       const OperatorType &             system_matrix, 
 *                       VectorType &                     dst, 
 *                       const VectorType &               src, 
 *                       const MultigridParameters &      mg_data, 
 *                       const hp::MappingCollection<dim> mapping_collection, 
 *                       const DoFHandler<dim> &          dof_handler, 
 *                       const hp::QCollection<dim> &     quadrature_collection) 
 *   { 
 * 
 * @endcode
 * 
 * 为每个多网格层次创建一个DoFHandler和操作符，以及，创建转移操作符。为了能够设置运算符，我们需要一组DoFHandler，通过p或h的全局粗化来创建。
 * 

 * 
 * 如果没有要求h-transfer，我们为`emplace_back()`函数提供一个空的删除器，因为我们的DoFHandler的Triangulation是一个外部字段，其析构器在其他地方被调用。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<DoFHandler<dim>>                     dof_handlers; 
 *     MGLevelObject<std::unique_ptr<OperatorType>>       operators; 
 *     MGLevelObject<MGTwoLevelTransfer<dim, VectorType>> transfers; 
 * 
 *     std::vector<std::shared_ptr<const Triangulation<dim>>> 
 *       coarse_grid_triangulations; 
 *     if (mg_data.transfer.perform_h_transfer) 
 *       coarse_grid_triangulations = 
 *         MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence( 
 *           dof_handler.get_triangulation()); 
 *     else 
 *       coarse_grid_triangulations.emplace_back( 
 *         const_cast<Triangulation<dim> *>(&(dof_handler.get_triangulation())), 
 *         [](auto &) {}); 
 * 
 * @endcode
 * 
 * 确定多栅格操作的总层数，并为所有层数分配足够的内存。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_h_levels = coarse_grid_triangulations.size() - 1; 
 * 
 *     const auto get_max_active_fe_degree = [&](const auto &dof_handler) { 
 *       unsigned int max = 0; 
 * 
 *       for (auto &cell : dof_handler.active_cell_iterators()) 
 *         if (cell->is_locally_owned()) 
 *           max = 
 *             std::max(max, dof_handler.get_fe(cell->active_fe_index()).degree); 
 * 
 *       return Utilities::MPI::max(max, MPI_COMM_WORLD); 
 *     }; 
 * 
 *     const unsigned int n_p_levels = 
 *       MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence( 
 *         get_max_active_fe_degree(dof_handler), mg_data.transfer.p_sequence) 
 *         .size(); 
 * 
 *     std::map<unsigned int, unsigned int> fe_index_for_degree; 
 *     for (unsigned int i = 0; i < dof_handler.get_fe_collection().size(); ++i) 
 *       { 
 *         const unsigned int degree = dof_handler.get_fe(i).degree; 
 *         Assert(fe_index_for_degree.find(degree) == fe_index_for_degree.end(), 
 *                ExcMessage("FECollection does not contain unique degrees.")); 
 *         fe_index_for_degree[degree] = i; 
 *       } 
 * 
 *     unsigned int minlevel   = 0; 
 *     unsigned int minlevel_p = n_h_levels; 
 *     unsigned int maxlevel   = n_h_levels + n_p_levels - 1; 
 * 
 *     dof_handlers.resize(minlevel, maxlevel); 
 *     operators.resize(minlevel, maxlevel); 
 *     transfers.resize(minlevel, maxlevel); 
 * 
 * @endcode
 * 
 * 从最小（最粗）到最大（最细）级别的循环，并相应地设置DoFHandler。我们从h层开始，在这里我们分布在越来越细的网格上的线性元素。
 * 

 * 
 * 
 * @code
 *     for (unsigned int l = 0; l < n_h_levels; ++l) 
 *       { 
 *         dof_handlers[l].reinit(*coarse_grid_triangulations[l]); 
 *         dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection()); 
 *       } 
 * 
 * @endcode
 * 
 * 在我们达到最细的网格后，我们将调整每一层的多项式度数。我们反向迭代我们的数据结构，从包含所有活动FE指数信息的最细网格开始。然后我们逐级降低每个单元的多项式度数。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0, l = maxlevel; i < n_p_levels; ++i, --l) 
 *       { 
 *         dof_handlers[l].reinit(dof_handler.get_triangulation()); 
 * 
 *         if (l == maxlevel) // finest level 
 *           { 
 *             auto &dof_handler_mg = dof_handlers[l]; 
 * 
 *             auto cell_other = dof_handler.begin_active(); 
 *             for (auto &cell : dof_handler_mg.active_cell_iterators()) 
 *               { 
 *                 if (cell->is_locally_owned()) 
 *                   cell->set_active_fe_index(cell_other->active_fe_index()); 
 *                 cell_other++; 
 *               } 
 *           } 
 *         else // coarse level 
 *           { 
 *             auto &dof_handler_fine   = dof_handlers[l + 1]; 
 *             auto &dof_handler_coarse = dof_handlers[l + 0]; 
 * 
 *             auto cell_other = dof_handler_fine.begin_active(); 
 *             for (auto &cell : dof_handler_coarse.active_cell_iterators()) 
 *               { 
 *                 if (cell->is_locally_owned()) 
 *                   { 
 *                     const unsigned int next_degree = 
 *                       MGTransferGlobalCoarseningTools:: 
 *                         create_next_polynomial_coarsening_degree( 
 *                           cell_other->get_fe().degree, 
 *                           mg_data.transfer.p_sequence); 
 *                     Assert(fe_index_for_degree.find(next_degree) != 
 *                              fe_index_for_degree.end(), 
 *                            ExcMessage("Next polynomial degree in sequence " 
 *                                       "does not exist in FECollection.")); 
 * 
 *                     cell->set_active_fe_index(fe_index_for_degree[next_degree]); 
 *                   } 
 *                 cell_other++; 
 *               } 
 *           } 
 * 
 *         dof_handlers[l].distribute_dofs(dof_handler.get_fe_collection()); 
 *       } 
 * 
 * @endcode
 * 
 * 接下来，我们将在每个多重网格层面上创建所有额外需要的数据结构。这涉及到确定具有同质Dirichlet边界条件的约束，并像在活动层上一样建立运算器。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<AffineConstraints<typename VectorType::value_type>> 
 *       constraints(minlevel, maxlevel); 
 * 
 *     for (unsigned int level = minlevel; level <= maxlevel; ++level) 
 *       { 
 *         const auto &dof_handler = dof_handlers[level]; 
 *         auto &      constraint  = constraints[level]; 
 * 
 *         IndexSet locally_relevant_dofs; 
 *         DoFTools::extract_locally_relevant_dofs(dof_handler, 
 *                                                 locally_relevant_dofs); 
 *         constraint.reinit(locally_relevant_dofs); 
 * 
 *         DoFTools::make_hanging_node_constraints(dof_handler, constraint); 
 *         VectorTools::interpolate_boundary_values(mapping_collection, 
 *                                                  dof_handler, 
 *                                                  0, 
 *                                                  Functions::ZeroFunction<dim>(), 
 *                                                  constraint); 
 *         constraint.close(); 
 * 
 *         VectorType dummy; 
 * 
 *         operators[level] = std::make_unique<OperatorType>(mapping_collection, 
 *                                                           dof_handler, 
 *                                                           quadrature_collection, 
 *                                                           constraint, 
 *                                                           dummy); 
 *       } 
 * 
 * @endcode
 * 
 * 根据多网格求解器类的需要，在单个算子中设置网格间算子和收集转移算子。
 * 

 * 
 * 
 * @code
 *     for (unsigned int level = minlevel; level < minlevel_p; ++level) 
 *       transfers[level + 1].reinit_geometric_transfer(dof_handlers[level + 1], 
 *                                                      dof_handlers[level], 
 *                                                      constraints[level + 1], 
 *                                                      constraints[level]); 
 * 
 *     for (unsigned int level = minlevel_p; level < maxlevel; ++level) 
 *       transfers[level + 1].reinit_polynomial_transfer(dof_handlers[level + 1], 
 *                                                       dof_handlers[level], 
 *                                                       constraints[level + 1], 
 *                                                       constraints[level]); 
 * 
 *     MGTransferGlobalCoarsening<dim, VectorType> transfer( 
 *       transfers, [&](const auto l, auto &vec) { 
 *         operators[l]->initialize_dof_vector(vec); 
 *       }); 
 * 
 * @endcode
 * 
 * 最后，继续用多网格法解决问题。
 * 

 * 
 * 
 * @code
 *     mg_solve(solver_control, 
 *              dst, 
 *              src, 
 *              mg_data, 
 *              dof_handler, 
 *              system_matrix, 
 *              operators, 
 *              transfer); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceProblem</code> class template</h3>
 * 

 * 
 * 现在，我们将最后声明这个程序的主类，它在随后的精炼函数空间上求解拉普拉斯方程。它的结构看起来很熟悉，因为它与  step-27  和  step-40  的主类类似。基本上只增加了两个。
 * 

 * 
 * - 持有系统矩阵的SparseMatrix对象已经被MatrixFree公式中的LaplaceOperator类对象所取代。
 * 

 * 
 * - 加入了一个 parallel::CellWeights, 的对象，它将帮助我们实现负载平衡。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class LaplaceProblem 
 *   { 
 *   public: 
 *     LaplaceProblem(const Parameters &parameters); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void initialize_grid(); 
 *     void setup_system(); 
 *     void print_diagnostics(); 
 *     void solve_system(); 
 *     void compute_indicators(); 
 *     void adapt_resolution(); 
 *     void output_results(const unsigned int cycle); 
 * 
 *     MPI_Comm mpi_communicator; 
 * 
 *     const Parameters prm; 
 * 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 *     DoFHandler<dim>                           dof_handler; 
 * 
 *     hp::MappingCollection<dim> mapping_collection; 
 *     hp::FECollection<dim>      fe_collection; 
 *     hp::QCollection<dim>       quadrature_collection; 
 *     hp::QCollection<dim - 1>   face_quadrature_collection; 
 * 
 *     IndexSet locally_owned_dofs; 
 *     IndexSet locally_relevant_dofs; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     LaplaceOperator<dim, double>               laplace_operator; 
 *     LinearAlgebra::distributed::Vector<double> locally_relevant_solution; 
 *     LinearAlgebra::distributed::Vector<double> system_rhs; 
 * 
 *     std::unique_ptr<FESeries::Legendre<dim>>    legendre; 
 *     std::unique_ptr<parallel::CellWeights<dim>> cell_weights; 
 * 
 *     Vector<float> estimated_error_per_cell; 
 *     Vector<float> hp_decision_indicators; 
 * 
 *     ConditionalOStream pcout; 
 *     TimerOutput        computing_timer; 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>LaplaceProblem</code> class implementation</h3>
 * 
 * <a name="Constructor"></a> 
 * <h4>Constructor</h4>
 * 

 * 
 * 构造函数以一个初始化器列表开始，该列表看起来与  step-40  的列表相似。我们再次准备好ConditionalOStream对象，只允许第一个进程在控制台输出任何东西，并正确初始化计算计时器。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   LaplaceProblem<dim>::LaplaceProblem(const Parameters &parameters) 
 *     : mpi_communicator(MPI_COMM_WORLD) 
 *     , prm(parameters) 
 *     , triangulation(mpi_communicator) 
 *     , dof_handler(triangulation) 
 *     , pcout(std::cout, 
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
 *     , computing_timer(mpi_communicator, 
 *                       pcout, 
 *                       TimerOutput::summary, 
 *                       TimerOutput::wall_times) 
 *   { 
 *     Assert(prm.min_h_level <= prm.max_h_level, 
 *            ExcMessage( 
 *              "Triangulation level limits have been incorrectly set up.")); 
 *     Assert(prm.min_p_degree <= prm.max_p_degree, 
 *            ExcMessage("FECollection degrees have been incorrectly set up.")); 
 * 
 * @endcode
 * 
 * 我们需要在构造函数的实际主体中为hp-functionality准备数据结构，并在参数结构的指定范围内为每个度数创建相应的对象。由于我们只处理非扭曲的矩形单元，在这种情况下，一个线性映射对象就足够了。
 * 

 * 
 * 在参数结构中，我们为函数空间以合理的分辨率运行的层级提供范围。多网格算法需要在最粗的层次上使用线性元素。所以我们从最低的多项式度数开始，用连续的高度数填充集合，直到达到用户指定的最大值。
 * 

 * 
 * 
 * @code
 *     mapping_collection.push_back(MappingQ1<dim>()); 
 * 
 *     for (unsigned int degree = 1; degree <= prm.max_p_degree; ++degree) 
 *       { 
 *         fe_collection.push_back(FE_Q<dim>(degree)); 
 *         quadrature_collection.push_back(QGauss<dim>(degree + 1)); 
 *         face_quadrature_collection.push_back(QGauss<dim - 1>(degree + 1)); 
 *       } 
 * 
 * @endcode
 * 
 * 由于我们的FECollection包含的有限元比我们想用于求解的有限元近似值要多，我们想限制活动FE指数可以操作的范围。为此，FECollection类允许注册一个层次结构，在p-精简和p-粗化的情况下，分别决定后续的和前面的有限元。 hp::Refinement 命名空间中的所有函数都会参考这个层次结构来确定未来的FE指数。我们将注册这样一个层次结构，它只对建议范围内的多项式程度的有限元起作用  <code>[min_p_degree, max_p_degree]</code>  。
 * 

 * 
 * 
 * @code
 *     const unsigned int min_fe_index = prm.min_p_degree - 1; 
 *     fe_collection.set_hierarchy( 
 * 
 *    /*下一个_index=  */ 
 *       [](const typename hp::FECollection<dim> &fe_collection, 
 *          const unsigned int                    fe_index) -> unsigned int { 
 *         return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 : 
 *                                                          fe_index; 
 *       }, 
 *     /*上一页_index=  */ 
 *       [min_fe_index](const typename hp::FECollection<dim> &, 
 *                      const unsigned int fe_index) -> unsigned int { 
 *         Assert(fe_index >= min_fe_index, 
 *                ExcMessage("Finite element is not part of hierarchy!")); 
 *         return (fe_index > min_fe_index) ? fe_index - 1 : fe_index; 
 *       }); 
 * 
 * @endcode
 * 
 * 我们以默认配置初始化 FESeries::Legendre 对象，以便进行平滑度估计。
 * 

 * 
 * 
 * @code
 *     legendre = std::make_unique<FESeries::Legendre<dim>>( 
 *       SmoothnessEstimator::Legendre::default_fe_series(fe_collection)); 
 * 
 * @endcode
 * 
 * 接下来的部分会很棘手。在执行细化的过程中，有几个hp-算法需要干扰三角形对象上的实际细化过程。我们通过将几个函数连接到 Triangulation::Signals: 信号，在实际细化过程中的不同阶段被调用，并触发所有连接的函数来做到这一点。我们需要这个功能来实现负载平衡和限制相邻单元的多项式度数。
 * 

 * 
 * 对于前者，我们希望给每个单元分配一个权重，这个权重与它未来的有限元的自由度数成正比。该库提供了一个类 parallel::CellWeights ，允许在细化过程中的正确位置轻松地附加单个权重，即在所有细化和粗化标志被正确设置为hp-adaptation之后，以及在即将发生的负载平衡的重新划分之前。可以注册一些函数，这些函数将以  $a (n_\text{dofs})^b$  提供的一对参数的形式附加权重  $(a,b)$  。我们在下文中注册了这样一个函数。每个单元在创建时将被赋予一个恒定的权重，这个值是1000（见  Triangulation::Signals::cell_weight).  ）。
 * 

 * 
 * 为了实现负载平衡，像我们使用的高效求解器应该与拥有的自由度数量成线性比例。此外，为了增加我们想要附加的权重的影响，确保单个权重将超过这个基础权重的数量级。我们相应地设置单元加权的参数。大的加权系数为 $10^6$ ，指数为 $1$  。
 * 

 * 
 * 
 * @code
 *     cell_weights = std::make_unique<parallel::CellWeights<dim>>( 
 *       dof_handler, 
 *       parallel::CellWeights<dim>::ndofs_weighting( 
 *         {prm.weighting_factor, prm.weighting_exponent})); 
 * 
 * @endcode
 * 
 * 在h-adaptive应用中，我们通过限制相邻单元的细化水平的差异为1来确保2:1的网格平衡。通过下面代码片段中的第二个调用，我们将确保相邻单元的p级数也是如此：未来有限元的级数不允许相差超过指定的差值。函数 hp::Refinement::limit_p_level_difference 可以处理这个问题，但需要与并行环境中的一个非常特殊的信号相连。问题是，我们需要知道网格的实际细化情况，以便相应地设置未来的FE指数。由于我们要求p4est神谕进行细化，我们需要确保Triangulation已经先用神谕的适应标志进行了更新。 parallel::distributed::TemporarilyMatchRefineFlags 的实例化在其生命期内正是如此。因此，我们将在限制p级差之前创建这个类的对象，并将相应的lambda函数连接到信号 Triangulation::Signals::post_p4est_refinement, 上，该信号将在神谕被完善之后，但在三角法被完善之前被触发。此外，我们指定这个函数将被连接到信号的前面，以确保修改在连接到同一信号的任何其他函数之前进行。
 * 

 * 
 * 
 * @code
 *     triangulation.signals.post_p4est_refinement.connect( 
 *       [&, min_fe_index]() { 
 *         const parallel::distributed::TemporarilyMatchRefineFlags<dim> 
 *           refine_modifier(triangulation); 
 *         hp::Refinement::limit_p_level_difference(dof_handler, 
 *                                                  prm.max_p_level_difference, 
 *                                                  /*包含=  */ min_fe_index);
 *       }, 
 *       boost::signals2::at_front); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProbleminitialize_grid"></a> 
 * <h4>LaplaceProblem::initialize_grid</h4>
 * 

 * 
 * 对于L型域，我们可以使用 GridGenerator::hyper_L() 这个函数，如 step-50 中所演示的。然而在二维的情况下，该函数只去除第一象限，而在我们的方案中我们需要去除第四象限。因此，我们将使用一个不同的函数 GridGenerator::subdivided_hyper_L() ，它给我们更多的选择来创建网格。此外，我们在制定该函数时，也会生成一个三维网格：二维L型域基本上会在正Z方向上拉长1。
 * 

 * 
 * 我们首先假装建立一个  GridGenerator::subdivided_hyper_rectangle().  我们需要提供的参数是左下角和右上角的点对象，以及基本网格在每个方向的重复次数。我们为前两个维度提供这些参数，对更高的第三维度单独处理。
 * 

 * 
 * 为了创建一个L型域，我们需要去除多余的单元。为此，我们相应地指定 <code>cells_to_remove</code> 。我们希望从负方向的每一个单元格中移除一个单元格，但从正的x方向移除一个。
 * 

 * 
 * 最后，我们提供与所提供的最小网格细化水平相对应的初始细化数。此外，我们相应地设置初始活动FE指数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::initialize_grid() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "initialize grid"); 
 * 
 *     std::vector<unsigned int> repetitions(dim); 
 *     Point<dim>                bottom_left, top_right; 
 *     for (unsigned int d = 0; d < dim; ++d) 
 *       if (d < 2) 
 *         { 
 *           repetitions[d] = 2; 
 *           bottom_left[d] = -1.; 
 *           top_right[d]   = 1.; 
 *         } 
 *       else 
 *         { 
 *           repetitions[d] = 1; 
 *           bottom_left[d] = 0.; 
 *           top_right[d]   = 1.; 
 *         } 
 * 
 *     std::vector<int> cells_to_remove(dim, 1); 
 *     cells_to_remove[0] = -1; 
 * 
 *     GridGenerator::subdivided_hyper_L( 
 *       triangulation, repetitions, bottom_left, top_right, cells_to_remove); 
 * 
 *     triangulation.refine_global(prm.min_h_level); 
 * 
 *     const unsigned int min_fe_index = prm.min_p_degree - 1; 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         cell->set_active_fe_index(min_fe_index); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * 这个函数看起来和 step-40 的函数完全一样，但是你会注意到没有系统矩阵以及围绕它的脚手架。相反，我们将在这里初始化 <code>laplace_operator</code> 中的MatrixFree公式。对于边界条件，我们将使用本教程前面介绍的Solution类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::setup_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "setup system"); 
 * 
 *     dof_handler.distribute_dofs(fe_collection); 
 * 
 *     locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
 * 
 *     locally_relevant_solution.reinit(locally_owned_dofs, 
 *                                      locally_relevant_dofs, 
 *                                      mpi_communicator); 
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
 * 
 *     constraints.clear(); 
 *     constraints.reinit(locally_relevant_dofs); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *     VectorTools::interpolate_boundary_values( 
 *       mapping_collection, dof_handler, 0, Solution<dim>(), constraints); 
 *     constraints.close(); 
 * 
 *     laplace_operator.reinit(mapping_collection, 
 *                             dof_handler, 
 *                             quadrature_collection, 
 *                             constraints, 
 *                             system_rhs); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemprint_diagnostics"></a> 
 * <h4>LaplaceProblem::print_diagnostics</h4>
 * 

 * 
 * 这是一个打印关于方程组及其划分的额外诊断的函数。除了通常的全局活动单元数和自由度外，我们还输出它们的局部等价物。为了规范输出，我们将用 Utilities::MPI::gather 操作将局部数量传达给第一个进程，然后由该进程输出所有信息。本地量的输出只限于前8个进程，以避免终端的杂乱。
 * 

 * 
 * 此外，我们想打印数值离散化中的多项式度数的频率。由于这些信息只存储在本地，我们将计算本地拥有的单元上的有限元，随后通过 Utilities::MPI::sum. 进行交流。
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::print_diagnostics() 
 *   { 
 *     const unsigned int first_n_processes = 
 *       std::min<unsigned int>(8, 
 *                              Utilities::MPI::n_mpi_processes(mpi_communicator)); 
 *     const bool output_cropped = 
 *       first_n_processes < Utilities::MPI::n_mpi_processes(mpi_communicator); 
 * 
 *     { 
 *       pcout << "   Number of active cells:       " 
 *             << triangulation.n_global_active_cells() << std::endl 
 *             << "     by partition:              "; 
 * 
 *       std::vector<unsigned int> n_active_cells_per_subdomain = 
 *         Utilities::MPI::gather(mpi_communicator, 
 *                                triangulation.n_locally_owned_active_cells()); 
 *       for (unsigned int i = 0; i < first_n_processes; ++i) 
 *         pcout << ' ' << n_active_cells_per_subdomain[i]; 
 *       if (output_cropped) 
 *         pcout << " ..."; 
 *       pcout << std::endl; 
 *     } 
 * 
 *     { 
 *       pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *             << std::endl 
 *             << "     by partition:              "; 
 * 
 *       std::vector<types::global_dof_index> n_dofs_per_subdomain = 
 *         Utilities::MPI::gather(mpi_communicator, 
 *                                dof_handler.n_locally_owned_dofs()); 
 *       for (unsigned int i = 0; i < first_n_processes; ++i) 
 *         pcout << ' ' << n_dofs_per_subdomain[i]; 
 *       if (output_cropped) 
 *         pcout << " ..."; 
 *       pcout << std::endl; 
 *     } 
 * 
 *     { 
 *       std::vector<types::global_dof_index> n_constraints_per_subdomain = 
 *         Utilities::MPI::gather(mpi_communicator, constraints.n_constraints()); 
 * 
 *       pcout << "   Number of constraints:        " 
 *             << std::accumulate(n_constraints_per_subdomain.begin(), 
 *                                n_constraints_per_subdomain.end(), 
 *                                0) 
 *             << std::endl 
 *             << "     by partition:              "; 
 *       for (unsigned int i = 0; i < first_n_processes; ++i) 
 *         pcout << ' ' << n_constraints_per_subdomain[i]; 
 *       if (output_cropped) 
 *         pcout << " ..."; 
 *       pcout << std::endl; 
 *     } 
 * 
 *     { 
 *       std::vector<unsigned int> n_fe_indices(fe_collection.size(), 0); 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         if (cell->is_locally_owned()) 
 *           n_fe_indices[cell->active_fe_index()]++; 
 * 
 *       Utilities::MPI::sum(n_fe_indices, mpi_communicator, n_fe_indices); 
 * 
 *       pcout << "   Frequencies of poly. degrees:"; 
 *       for (unsigned int i = 0; i < fe_collection.size(); ++i) 
 *         if (n_fe_indices[i] > 0) 
 *           pcout << ' ' << fe_collection[i].degree << ":" << n_fe_indices[i]; 
 *       pcout << std::endl; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve_system"></a> 
 * <h4>LaplaceProblem::solve_system</h4>
 * 

 * 
 * 围绕解决方案的脚手架与  step-40  的类似。我们准备一个符合MatrixFree要求的向量，并收集本地相关的自由度，我们解决了方程系统。解决方法是通过前面介绍的函数进行的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::solve_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "solve system"); 
 * 
 *     LinearAlgebra::distributed::Vector<double> completely_distributed_solution; 
 *     laplace_operator.initialize_dof_vector(completely_distributed_solution); 
 * 
 *     SolverControl solver_control(system_rhs.size(), 
 *                                  prm.tolerance_factor * system_rhs.l2_norm()); 
 * 
 *     solve_with_gmg(solver_control, 
 *                    laplace_operator, 
 *                    completely_distributed_solution, 
 *                    system_rhs, 
 *                    prm.mg_data, 
 *                    mapping_collection, 
 *                    dof_handler, 
 *                    quadrature_collection); 
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations." 
 *           << std::endl; 
 * 
 *     constraints.distribute(completely_distributed_solution); 
 * 
 *     locally_relevant_solution.copy_locally_owned_data_from( 
 *       completely_distributed_solution); 
 *     locally_relevant_solution.update_ghost_values(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemcompute_indicators"></a> 
 * <h4>LaplaceProblem::compute_indicators</h4>
 * 

 * 
 * 这个函数只包含其他教程中典型的 <code>refine_grid</code> 函数的一部分，在这个意义上是新的。在这里，我们将只计算与实际细化网格相适应的所有指标。我们这样做的目的是将所有的指标写到文件系统中，以便为以后储存。
 * 

 * 
 * 由于我们处理的是一个椭圆问题，我们将再次利用KellyErrorEstimator，但有一点不同。修改底层面积分的缩放系数，使其取决于相邻元素的实际多项式程度，这对hp-adaptive应用是有利的  @cite davydov2017hp  。我们可以通过指定你所注意到的附加参数中的最后一个参数来做到这一点。其他的实际上只是默认的。
 * 

 * 
 * 为了hp-adaptation的目的，我们将用教程介绍中的策略来计算平滑度估计，并使用 SmoothnessEstimator::Legendre. 中的实现 在参数结构中，我们将最小多项式度数设置为2，因为似乎平滑度估计算法在处理线性元素时有问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::compute_indicators() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "compute indicators"); 
 * 
 *     estimated_error_per_cell.grow_or_shrink(triangulation.n_active_cells()); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       face_quadrature_collection, 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       locally_relevant_solution, 
 *       estimated_error_per_cell, 
 *       /*component_mask=  */ 
 *       ComponentMask(), 
 *       /*coefficients=  */ 
 *       nullptr, 
 *       /*n_threads=*/
 *       numbers::invalid_unsigned_int,  
 *       /*subdomain_id=*/ 
 *       numbers::invalid_subdomain_id,  
 *       /*material_id=*/ 
 *       numbers::invalid_material_id,  
 *       /*策略=  */ 
 *       KellyErrorEstimator<dim>::Strategy::face_diameter_over_twice_max_degree); 
 * 
 *     hp_decision_indicators.grow_or_shrink(triangulation.n_active_cells()); 
 *     SmoothnessEstimator::Legendre::coefficient_decay(*legendre, 
 *                                                      dof_handler, 
 *                                                      locally_relevant_solution, 
 *                                                      hp_decision_indicators); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemadapt_resolution"></a> 
 * <h4>LaplaceProblem::adapt_resolution</h4>
 * 

 * 
 * 有了之前计算出的指标，我们最终将标记所有单元进行适应，同时在这个函数中执行细化。和以前的教程一样，我们将使用 "固定数字 "策略，但现在是针对hp-adaptation。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::adapt_resolution() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "adapt resolution"); 
 * 
 * @endcode
 * 
 * 首先，我们将根据每个单元的误差估计值来设置细化和粗化标志。这里没有什么新东西。
 * 

 * 
 * 我们将使用在其他deal.II教程中阐述过的一般细化和粗化比例：使用固定数字策略，我们将标记所有单元中的30%进行细化，3%进行粗化，如参数结构中提供的。
 * 

 * 
 * 
 * @code
 *     parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
 *       triangulation, 
 *       estimated_error_per_cell, 
 *       prm.refine_fraction, 
 *       prm.coarsen_fraction); 
 * 
 * @endcode
 * 
 * 接下来，我们将对hp-adaptation进行所有调整。我们想细化和粗化那些在上一步中被标记的单元，但需要决定是通过调整网格分辨率还是调整多项式程度来实现。
 * 

 * 
 * 下一个函数调用根据之前计算的平滑度指标设置未来的FE指数，作为p-adaptation指标。这些指数将只设置在那些分配了细化或粗化标志的单元上。
 * 

 * 
 * 对于p-adaptation分数，我们将采取一个有根据的猜测。由于我们只期望在我们的方案中出现一个单一的奇点，即在域的原点，而在其他任何地方都有一个平滑的解决方案，所以我们希望强烈倾向于使用p-adaptation而不是h-adaptation。这反映在我们对p-精简和p-粗化都选择了90%的分数。
 * 

 * 
 * 
 * @code
 *     hp::Refinement::p_adaptivity_fixed_number(dof_handler, 
 *                                               hp_decision_indicators, 
 *                                               prm.p_refine_fraction, 
 *                                               prm.p_coarsen_fraction); 
 * 
 * @endcode
 * 
 * 在这个阶段，我们既有未来的FE指数，也有经典的细化和粗化标志，后者将由 Triangulation::execute_coarsening_and_refinement() 解释为h-适应性。我们希望只对细胞施加一种适应，这就是下一个函数将为我们解决的问题。简而言之，在分配有两种类型指标的单元格上，我们将倾向于p-适应的那一种，并删除h-适应的那一种。
 * 

 * 
 * 
 * @code
 *     hp::Refinement::choose_p_over_h(dof_handler); 
 * 
 * @endcode
 * 
 * 设置完所有指标后，我们将删除那些超过参数结构中提供的水平范围的指定限制的指标。由于提供的有限元数量有限，这种限制自然会出现在p-adaptation中。此外，我们在构造函数中为p-adaptation注册了一个自定义层次结构。现在，我们需要像  step-31  中那样，在h-adaptive的上下文中手动完成。
 * 

 * 
 * 我们将遍历指定的最小和最大层次上的所有单元格，并删除相应的标志。作为一种选择，我们也可以通过相应地设置未来的FE指数来标记这些单元的p适应性，而不是简单地清除细化和粗化的标志。
 * 

 * 
 * 
 * @code
 *     Assert(triangulation.n_levels() >= prm.min_h_level + 1 && 
 *              triangulation.n_levels() <= prm.max_h_level + 1, 
 *            ExcInternalError()); 
 * 
 *     if (triangulation.n_levels() > prm.max_h_level) 
 *       for (const auto &cell : 
 *            triangulation.active_cell_iterators_on_level(prm.max_h_level)) 
 *         cell->clear_refine_flag(); 
 * 
 *     for (const auto &cell : 
 *          triangulation.active_cell_iterators_on_level(prm.min_h_level)) 
 *       cell->clear_coarsen_flag(); 
 * 
 * @endcode
 * 
 * 最后，我们就剩下执行粗化和细化了。在这里，不仅网格会被更新，而且所有以前的未来FE指数也会变得活跃。
 * 

 * 
 * 记得我们在构造函数中为三角化信号附加了函数，将在这个函数调用中被触发。所以会有更多的事情发生：加权重新分区将被执行以确保负载平衡，以及我们将限制相邻单元之间的p级差。
 * 

 * 
 * 
 * @code
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemoutput_results"></a> 
 * <h4>LaplaceProblem::output_results</h4>
 * 

 * 
 * 在并行应用中向文件系统写入结果的工作方式与  step-40  中完全相同。除了我们在整个教程中准备的数据容器外，我们还想写出网格上每个有限元的多项式程度，以及每个单元所属的子域。我们在这个函数的范围内为此准备必要的容器。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "output results"); 
 * 
 *     Vector<float> fe_degrees(triangulation.n_active_cells()); 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         fe_degrees(cell->active_cell_index()) = cell->get_fe().degree; 
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells()); 
 *     for (auto &subd : subdomain) 
 *       subd = triangulation.locally_owned_subdomain(); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(locally_relevant_solution, "solution"); 
 *     data_out.add_data_vector(fe_degrees, "fe_degree"); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 *     data_out.add_data_vector(estimated_error_per_cell, "error"); 
 *     data_out.add_data_vector(hp_decision_indicators, "hp_indicator"); 
 *     data_out.build_patches(mapping_collection); 
 * 
 *     data_out.write_vtu_with_pvtu_record( 
 *       "./", "solution", cycle, mpi_communicator, 2, 1); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * 实际的运行函数看起来又和  step-40  非常相似。唯一增加的是实际循环之前的括号内的部分。在这里，我们将预先计算Legendre变换矩阵。一般来说，每当需要某个矩阵时，这些矩阵将通过懒惰分配的方式进行实时计算。然而，出于计时的目的，我们希望在实际的时间测量开始之前，一次性地计算它们。因此，我们将把它们的计算指定为自己的范围。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::run() 
 *   { 
 *     pcout << "Running with Trilinos on " 
 *           << Utilities::MPI::n_mpi_processes(mpi_communicator) 
 *           << " MPI rank(s)..." << std::endl; 
 * 
 *     { 
 *       pcout << "Calculating transformation matrices..." << std::endl; 
 *       TimerOutput::Scope t(computing_timer, "calculate transformation"); 
 *       legendre->precalculate_all_transformation_matrices(); 
 *     } 
 * 
 *     for (unsigned int cycle = 0; cycle < prm.n_cycles; ++cycle) 
 *       { 
 *         pcout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           initialize_grid(); 
 *         else 
 *           adapt_resolution(); 
 * 
 *         setup_system(); 
 * 
 *         print_diagnostics(); 
 * 
 *         solve_system(); 
 * 
 *         compute_indicators(); 
 * 
 *         if (Utilities::MPI::n_mpi_processes(mpi_communicator) <= 32) 
 *           output_results(cycle); 
 * 
 *         computing_timer.print_summary(); 
 *         computing_timer.reset(); 
 * 
 *         pcout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step75 
 * 
 * @endcode
 * 
 * 
 * <a name="main"></a> 
 * <h4>main()</h4>
 * 

 * 
 * 最后一个函数是 <code>main</code> 函数，它将最终创建并运行一个LaplaceOperator实例。它的结构与其他大多数教程程序相似。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step75; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       Parameters        prm; 
 *       LaplaceProblem<2> laplace_problem(prm); 
 *       laplace_problem.run(); 
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
examples/step-75/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你在释放模式下，在四个进程上用给定的参数运行该程序时，你的终端输出应该是这样的。

@code
Running with Trilinos on 4 MPI rank(s)...
Calculating transformation matrices...
Cycle 0:
   Number of active cells:       3072
     by partition:               768 768 768 768
   Number of degrees of freedom: 12545
     by partition:               3201 3104 3136 3104
   Number of constraints:        542
     by partition:               165 74 138 165
   Frequencies of poly. degrees: 2:3072
   Solved in 7 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.598s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| calculate transformation        |         1 |    0.0533s |       8.9% |
| compute indicators              |         1 |    0.0177s |         3% |
| initialize grid                 |         1 |    0.0397s |       6.6% |
| output results                  |         1 |    0.0844s |        14% |
| setup system                    |         1 |    0.0351s |       5.9% |
| solve system                    |         1 |     0.362s |        61% |
+---------------------------------+-----------+------------+------------+



Cycle 1:
   Number of active cells:       3351
     by partition:               875 761 843 872
   Number of degrees of freedom: 18223
     by partition:               4535 4735 4543 4410
   Number of constraints:        1202
     by partition:               303 290 326 283
   Frequencies of poly. degrees: 2:2523 3:828
   Solved in 7 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.442s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| adapt resolution                |         1 |    0.0189s |       4.3% |
| compute indicators              |         1 |    0.0135s |         3% |
| output results                  |         1 |     0.064s |        14% |
| setup system                    |         1 |    0.0232s |       5.2% |
| solve system                    |         1 |     0.322s |        73% |
+---------------------------------+-----------+------------+------------+



...



Cycle 7:
   Number of active cells:       5610
     by partition:               1324 1483 1482 1321
   Number of degrees of freedom: 82062
     by partition:               21116 19951 20113 20882
   Number of constraints:        14383
     by partition:               3825 3225 3557 3776
   Frequencies of poly. degrees: 2:1130 3:1283 4:2727 5:465 6:5
   Solved in 7 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.932s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| adapt resolution                |         1 |    0.0182s |       1.9% |
| compute indicators              |         1 |    0.0173s |       1.9% |
| output results                  |         1 |    0.0572s |       6.1% |
| setup system                    |         1 |    0.0252s |       2.7% |
| solve system                    |         1 |     0.813s |        87% |
+---------------------------------+-----------+------------+------------+
@endcode



当用更多的进程运行代码时，你会注意到活动单元和自由度的数量有轻微的差异。这是由于求解器和预处理程序取决于问题的分区，这可能会导致最后一位数的解决方案的微小差异，并最终产生不同的适应行为。

此外，尽管有hp-adaptation，求解器的迭代次数在所有周期中都保持不变，这表明所提出的算法的稳健性，并有望在更大的问题规模和更多的进程中具有良好的可扩展性。

让我们看一下程序的图形输出。在给定参数配置的所有细化循环之后，实际离散的函数空间看起来如下，左边是其在12个进程上的分区，右边是有限元的多项式程度。在左图中，每种颜色代表一个独特的子域。在右图中，最浅的颜色对应于多项式的2度，最深的对应于6度。

<div class="twocolumn" style="width: 80%; text-align: center;"> <div> <img src="https://www.dealii.org/images/steps/developer/step-75.subdomains-07.svg" alt="七次细化后的分区。"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-75.fedegrees-07.svg" alt="七次细化后的局部近似度。"> </div> <div>




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Differenthpdecisionstrategies"></a><h4>Different hp-decision strategies</h4>


deal.II库提供了多种策略来决定对单元格施加哪种类型的适应：要么调整网格分辨率，要么改变多项式程度。我们在本教程中只介绍了<i>Legendre
coefficient decay</i>策略，而Step-27则演示了相同想法的<i>Fourier</i>等值。

有关这些策略的概述，请参见步骤27的 "扩展的可能性 "部分，或相应的文件的详细描述。

在这里，提到了另一个迄今为止还没有在任何教程中展示过的策略：基于<i>refinement history</i>的策略。这种方法在并行分布式应用中的使用比其他方法更棘手，所以我们将强调随之而来的挑战。我们需要有关细化标志的最终状态的信息，并且我们需要在细化的网格之间转移解决方案。对于前者，我们需要将 hp::Refinement::predict_error() 函数附加到 Triangulation::Signals::post_p4est_refinement 信号上，其方式是将<i>after</i>的 hp::Refinement::limit_p_level_difference() 函数调用。在这个阶段，所有的细化标志和未来的FE指数都被终止设置，并且可以对误差进行可靠预测。然后，预测的误差需要借助于 parallel::distributed::CellDataTransfer. 在细化网格之间进行转移。

试着在本教程中实现这些策略之一，并观察结果的微妙变化。你会注意到，所有的策略都能够识别出重心角附近的奇点，并且会在这些区域进行 $h$ -精化，而在体域中更倾向于 $p$ -精化。这些策略的详细比较见于  @cite fehling2020  。




<a name="Solvewithmatrixbasedmethods"></a><h4>Solve with matrix-based methods</h4>


本教程只关注无矩阵策略。然而，所有的hp自适应算法在并行分布式背景下也可以使用基于矩阵的方法。

为了创建一个系统矩阵，你可以使用 LaplaceOperator::get_system_matrix() 函数，或者使用类似于步骤27的 <code>assemble_system()</code> 函数。然后你可以像往常一样将系统矩阵传递给求解器。

你可以对基于矩阵和无矩阵的实现结果进行计时，量化速度提升，并说服自己哪种变体更快。




<a name="Multigridvariants"></a><h4>Multigrid variants</h4>


为了简单起见，我们将自己限制在单一类型的粗网格求解器（带AMG的CG）、平滑器（带点Jacobi预处理的Chebyshev平滑器）以及多网格算法中的几何粗化方案（全局粗化）。请自由尝试替代方案并调查其性能和稳健性。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-75.cc"
*/
