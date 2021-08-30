/**
@page step_40 The step-40 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Thetestcase">The testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The <code>LaplaceProblem</code> class template</a>
        <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The <code>LaplaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#Constructor">Constructor</a>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#LaplaceProblemrefine_grid">LaplaceProblem::refine_grid</a>
        <li><a href="#LaplaceProblemoutput_results">LaplaceProblem::output_results</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
        <li><a href="#main">main()</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-40/doc/intro.dox

 <br> 

<i>This program was contributed by Timo Heister, Martin Kronbichler and Wolfgang
Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation or of The
California Institute of Technology.
</i>




 @note  作为这个程序的前提条件，你需要同时安装PETSc和p4est库。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II与这两个附加库的安装。还要注意的是，为了正常工作，本程序需要访问实现代数多网格的Hypre预处理程序包；它可以作为PETSc的一部分安装，但必须在配置PETSc时明确启用；参见PETSc安装说明中的链接页面。


<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{41.5,41.75} 

鉴于今天的计算机，大多数有限元计算可以在一台机器上完成。因此，以前的大多数教程程序只显示了这一点，可能是在一些处理器之间进行分工，但这些处理器都可以访问相同的共享内存空间。也就是说，有些问题对于单台机器来说实在是太大了，在这种情况下，必须以适当的方式将问题分割给多台机器，每台机器都为整体贡献自己的一部分。在第17步和第18步中展示了一个简单的方法，我们展示了一个程序如何使用<a
href="http://www.mpi-forum.org/" target="_top">MPI</a>来并行组装线性系统，存储它，解决它，并计算误差估计。所有这些操作的扩展都是相对微不足道的（关于操作 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"），但是有一个明显的缺点：为了使这个实现适度简单，每个MPI处理器都必须保留自己的整个Triangulation和DoFHandler对象的副本。因此，虽然我们可以怀疑（有充分的理由）上面列出的操作可以扩展到成千上万的计算机和数十亿个单元和数十亿个自由度的问题规模，但在每一个最后的处理器上为这成千上万的计算机所解决的整个问题建立一个大的网格显然是不能扩展的：这将需要永远，也许更重要的是没有一台机器会有足够的内存来存储一个有十亿个单元的网格（至少在写这篇文章时没有）。在现实中，像第17步和第18步这样的程序不可能在超过100或200个处理器上运行，即使在那里，存储Triangulation和DoFHandler对象也会消耗每台机器上的绝大部分内存。

因此，我们需要以不同的方式来处理这个问题：为了扩展到非常大的问题，每个处理器只能存储自己的一小块三角形和DoFHandler对象。deal.II在 parallel::distributed 命名空间和其中的类中实现了这样一个方案。它建立在一个外部库上，<a
href="http://www.p4est.org/">p4est</a>（对表达式<i>parallel forest</i>的发挥，描述了将分层构造的网格作为四叉树或八叉树的森林进行并行存储）。你需要<a
href="../../external-libs/p4est.html">install and configure p4est</a>，但除此之外，它的所有工作原理都隐藏在deal.II的表面之下。

本质上， parallel::distributed::Triangulation 类和DoFHandler类中的代码所做的是分割全局网格，使每个处理器只存储其 "拥有 "的一小部分，以及围绕其拥有的单元的一层 "幽灵 "单元。在我们想要解决偏微分方程的领域的其余部分发生了什么，对每个处理器来说都是未知的，如果需要这些信息，只能通过与其他机器的交流来推断。这意味着我们还必须以不同于例如第17步和第18步的方式来思考问题：例如，没有一个处理器可以拥有用于后处理的整个解矢量，程序的每一部分都必须被并行化，因为没有一个处理器拥有顺序操作所需的所有信息。

在 @ref distributed 文档模块中描述了这种并行化如何发生的一般概述。在阅读本程序的源代码之前，你应该先阅读它，以获得一个顶层的概述。在 @ref distributed_paper "分布式计算论文 "中也提供了关于我们将在程序中使用的许多术语的简明讨论。也许值得一读，以了解本程序内部如何工作的背景信息。




<a name="Thetestcase"></a><h3>The testcase</h3>


这个程序基本上重新解决了我们在步骤6中已经做的事情，即它解决了拉普拉斯方程

@f{align*}


  -\Delta u &= f \qquad &&\text{in}\ \Omega=[0,1]^2, \\
  u &= 0 \qquad &&\text{on}\ \partial\Omega.


@f}

当然不同的是，现在我们要在一个可能有十亿个单元，有十亿个左右自由度的网格上这样做。毫无疑问，对于这样一个简单的问题，这样做是完全愚蠢的，但毕竟一个教程程序的重点不是做一些有用的东西，而是展示如何使用deal.II来实现有用的程序。尽管如此，为了使事情至少有一点点有趣，我们选择右侧为一个不连续的函数。

@f{align*}
  f(x,y)
  =
  \left\{
  \begin{array}{ll}
    1 & \text{if}\ y > \frac 12 + \frac 14 \sin(4\pi x), \\


    -1 & \text{otherwise},
  \end{array}
  \right.


@f}

使得解沿着蜿蜒穿过域的正弦线有一个奇点。因此，网格的细化将集中在这条线上。你可以在下面结果部分的网格图中看到这一点。

与其在这里继续做冗长的介绍，不如让我们直接进入程序代码。如果你已经读完了步骤6和 @ref distributed 文档模块，大部分将要发生的事情你应该已经熟悉了。事实上，比较这两个程序，你会发现在%parallel中工作所需的额外努力几乎是微不足道的：这两个程序的代码行数差不多（尽管步骤6在处理系数和输出方面花费了更多的空间）。在任何情况下，下面的评论将只针对使step-40与step-6不同的事情，而且在 @ref distributed 文档模块中还没有涵盖。




 @note  这个程序将能够在你想扔给它的多少个处理器上进行计算，以及你有多少内存和耐心来解决多大的问题。然而，<i>is</i>有一个限制：未知数的数量不能超过可以用类型 types::global_dof_index. 的对象存储的最大数量。默认情况下，这是<code>unsigned int</code>的别名，在今天大多数机器上是一个32位的整数，限制了你大约40亿（实际上，由于这个程序使用PETSc，你将被限制在一半，因为PETSc使用有符号整数）。然而，这可以在配置过程中改变为使用64位整数，见ReadMe文件。这将使问题的大小在短期内不太可能超过。


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
 * 我们在这个程序中需要的大部分包含文件已经在以前的程序中讨论过了。特别是，以下所有的文件都应该已经是熟悉的朋友了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/timer.h> 
 * 
 * #include <deal.II/lac/generic_linear_algebra.h> 
 * 
 * @endcode
 * 
 * 这个程序可以使用PETSc或Trilinos来满足其并行代数的需要。默认情况下，如果deal.II已经被配置为PETSc，它将使用PETSc。否则，下面几行将检查deal.II是否已被配置为Trilinos，并采用它。
 * 

 * 
 * 但是在某些情况下，即使deal.II也被配置为PETSc，你还是想使用Trilinos，例如，比较这两个库的性能。要做到这一点，请在源代码中添加以下的\#define。
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *  #define FORCE_USE_OF_TRILINOS
 *  @endcode
 * </div>
 * 

 * 
 * 使用这个逻辑，下面几行将导入PETSc或Trilinos包装器到命名空间`LA`（代表 "线性代数"）。在前一种情况下，我们还要定义宏 `USE_PETSC_LA`，这样我们就可以检测到我们是否在使用PETSc（参见solve()中需要用到的例子）。
 * 

 * 
 * 
 * @code
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
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 然而，下面这些将是新的，或在新的角色中使用。让我们来看看它们。其中第一个将提供 Utilities::System 命名空间的工具，我们将用它来查询诸如与当前MPI宇宙相关的处理器数量，或者这个作业运行的处理器在这个宇宙中的编号。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/utilities.h> 
 * 
 * @endcode
 * 
 * 下一个提供了一个类，ConditionOStream，它允许我们编写代码，将东西输出到一个流中（例如在每个处理器上的 <code>std::cout</code> ，但在除了一个处理器以外的所有处理器上都将文本扔掉。我们可以通过简单地在每个可能产生输出的地方前面放一个 <code>if</code> 语句来实现同样的目的，但这并不能使代码更漂亮。此外，这个处理器是否应该向屏幕输出的条件每次都是一样的--因此，把它放在产生输出的语句中应该是很简单的。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * 
 * @endcode
 * 
 * 在这些预演之后，这里变得更加有趣。正如在 @ref distributed 模块中提到的，在大量处理器上解决问题的一个基本事实是，任何处理器都不可能存储所有的东西（例如，关于网格中所有单元的信息，所有的自由度，或者解向量中所有元素的值）。相反，每个处理器都会<i>own</i>其中的几个，如果有必要，还可能<i>know</i>另外几个，例如，位于与该处理器自己拥有的单元相邻的那些单元。我们通常称后者为<i>ghost cells</i>、<i>ghost nodes</i>或<i>ghost elements of a vector</i>。这里讨论的重点是，我们需要有一种方法来表明一个特定的处理器拥有或需要知道哪些元素。这就是IndexSet类的领域：如果总共有 $N$ 个单元、自由度或向量元素，与（非负）积分指数 $[0,N)$ 相关，那么当前处理器拥有的元素集以及它需要了解的（可能更大）指数集都是集合 $[0,N)$ 的子集。IndexSet是一个类，它以一种有效的格式存储这个集合的子集。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/index_set.h> 
 * 
 * @endcode
 * 
 * 下一个头文件是一个单一的函数所必需的，  SparsityTools::distribute_sparsity_pattern.  这个函数的作用将在下面解释。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparsity_tools.h> 
 * 
 * @endcode
 * 
 * 最后两个新的头文件提供了类 parallel::distributed::Triangulation ，它提供了分布在可能非常多的处理器上的网格，而第二个文件提供了命名空间 parallel::distributed::GridRefinement ，它提供了可以自适应细化这种分布式网格的函数。
 * 

 * 
 * 
 * @code
 * #include <deal.II/distributed/tria.h> 
 * #include <deal.II/distributed/grid_refinement.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step40 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceProblem</code> class template</h3>
 * 

 * 
 * 接下来我们来声明这个程序的主类。它的结构几乎与 step-6 的教程程序一模一样。唯一显著的区别是。
 * 

 * 
 * --  <code>mpi_communicator</code> 变量，它描述了我们希望这段代码运行在哪一组处理器上。在实践中，这将是MPI_COMM_WORLD，即批处理调度系统分配给这个特定作业的所有处理器。
 * 

 * 
 * - ConditionOStream类型的 <code>pcout</code> 变量的存在。
 * 

 * 
 * - 明显使用 parallel::distributed::Triangulation 而不是Triangulation。
 * 

 * 
 * - 两个IndexSet对象的存在，表示我们在当前处理器上拥有哪些自由度集（以及解和右手向量的相关元素），以及我们需要哪些（作为幽灵元素）来使本程序中的算法工作。
 * 

 * 
 * - 现在所有的矩阵和向量都是分布式的。我们使用PETSc或Trilinos包装类，这样我们就可以使用Hypre（使用PETSc）或ML（使用Trilinos）提供的复杂的预处理器之一。请注意，作为这个类的一部分，我们存储的解向量不仅包含当前处理器拥有的自由度，还包括（作为鬼魂元素）所有对应于 "本地相关 "自由度的向量元素（即所有生活在本地拥有的单元或围绕它的鬼魂单元层的自由度）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class LaplaceProblem 
 *   { 
 *   public: 
 *     LaplaceProblem(); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     MPI_Comm mpi_communicator; 
 * 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 * 
 *     FE_Q<dim>       fe; 
 *     DoFHandler<dim> dof_handler; 
 * 
 *     IndexSet locally_owned_dofs; 
 *     IndexSet locally_relevant_dofs; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     LA::MPI::SparseMatrix system_matrix; 
 *     LA::MPI::Vector       locally_relevant_solution; 
 *     LA::MPI::Vector       system_rhs; 
 * 
 *     ConditionalOStream pcout; 
 *     TimerOutput        computing_timer; 
 *   }; 
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
 * 构造函数和析构函数是相当微不足道的。除了我们在 step-6 中所做的，我们将我们想要工作的处理器集合设置为所有可用的机器（MPI_COMM_WORLD）；要求三角化以确保网格保持平滑并自由精炼岛屿，例如；并初始化 <code>pcout</code> 变量，只允许处理器0输出任何东西。最后一块是初始化一个定时器，我们用它来决定程序的不同部分需要多少计算时间。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   LaplaceProblem<dim>::LaplaceProblem() 
 *     : mpi_communicator(MPI_COMM_WORLD) 
 *     , triangulation(mpi_communicator, 
 *                     typename Triangulation<dim>::MeshSmoothing( 
 *                       Triangulation<dim>::smoothing_on_refinement | 
 *                       Triangulation<dim>::smoothing_on_coarsening)) 
 *     , fe(2) 
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
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * 下面这个函数可以说是整个程序中最有趣的一个，因为它涉及到了%并行  step-40  和顺序  step-6  的核心区别。
 * 

 * 
 * 在顶部我们做了我们一直在做的事情：告诉DoFHandler对象来分配自由度。由于我们在这里使用的三角测量是分布式的，DoFHandler对象足够聪明，它认识到在每个处理器上只能在它所拥有的单元上分配自由度；接下来是一个交换步骤，处理器互相告诉对方关于ghost单元的自由度。结果是DoFHandler知道本地拥有的单元和幽灵单元（即与本地拥有的单元相邻的单元）的自由度，但对更远的单元则一无所知，这与分布式计算的基本理念一致，即没有处理器可以知道所有的事情。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::setup_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "setup"); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 接下来的两行提取了一些我们以后需要的信息，即两个索引集，提供了关于哪些自由度为当前处理器所拥有的信息（这些信息将被用来初始化解和右手向量以及系统矩阵，表明哪些元素要存储在当前处理器上，哪些要期望存储在其他地方）；以及一个索引集，表明哪些自由度是本地相关的（即生活在当前处理器所拥有的单元上或本地所拥有的单元周围的鬼魂单元上；我们将把这些自由度存储在当前处理器上。 例如，生活在当前处理器拥有的单元上或本地拥有的单元周围的幽灵单元层上；例如，我们需要所有这些自由度来估计本地单元的误差）。)
 * 

 * 
 * 
 * @code
 *     locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 *     DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 
 * 
 * @endcode
 * 
 * 接下来，让我们初始化解和右手边的向量。如上所述，我们寻求的解向量不仅存储了我们自己的元素，还存储了幽灵条目；另一方面，右手向量只需要有当前处理器拥有的条目，因为我们所做的只是向其中写入，而不是从其中读取本地拥有的单元（当然，线性求解器会从其中读取，但它们并不关心自由度的几何位置）。
 * 

 * 
 * 
 * @code
 *     locally_relevant_solution.reinit(locally_owned_dofs, 
 *                                      locally_relevant_dofs, 
 *                                      mpi_communicator); 
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
 * 
 * @endcode
 * 
 * 下一步是计算悬挂节点和边界值约束，我们将其合并为一个存储所有约束的对象。
 * 

 * 
 * 就像在%parallel中的所有其他事情一样，口头禅必须是：没有一个处理器可以存储整个宇宙的所有信息。因此，我们需要告诉AffineConstraints对象哪些自由度可以存储约束条件，哪些可以不期望存储任何信息。在我们的例子中，正如 @ref distributed 模块所解释的，我们需要在每个处理器上关心的自由度是本地相关的自由度，所以我们把这个传递给 AffineConstraints::reinit 函数。顺便提一下，如果你忘记传递这个参数，AffineConstraints类将分配一个长度等于它目前看到的最大自由度索引的数组。对于MPI进程数很高的处理器来说，这可能是非常大的 -- 也许是数十亿的数量级。然后，程序将为这个单一的数组分配比其他所有操作加起来还要多的内存。
 * 

 * 
 * 
 * @code
 *     constraints.clear(); 
 *     constraints.reinit(locally_relevant_dofs); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              Functions::ZeroFunction<dim>(), 
 *                                              constraints); 
 *     constraints.close(); 
 * 
 * @endcode
 * 
 * 这个函数的最后一部分涉及到用伴随的稀疏模式初始化矩阵。和以前的教程程序一样，我们使用DynamicSparsityPattern作为一个中介，然后用它来初始化系统矩阵。为了做到这一点，我们必须告诉稀疏模式它的大小，但如上所述，所产生的对象不可能为每个全局自由度存储哪怕一个指针；我们最好的希望是它能存储每个局部相关自由度的信息，即所有我们在组装矩阵的过程中可能接触到的自由度（ @ref distributed_paper "分布式计算论文 "有很长的讨论，为什么我们真的需要局部相关自由度，而不是在此背景下的小的局部活动自由度集）。
 * 

 * 
 * 所以我们告诉稀疏模式它的大小和要存储什么自由度，然后要求 DoFTools::make_sparsity_pattern 来填充它（这个函数忽略了所有不属于本地的单元，模仿我们下面在装配过程中的做法）。在这之后，我们调用一个函数，在处理器之间交换这些稀疏模式的条目，以便最后每个处理器真正知道它将拥有的那部分有限元矩阵中的所有条目。最后一步是用稀疏模式初始化矩阵。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern dsp(locally_relevant_dofs); 
 * 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
 *     SparsityTools::distribute_sparsity_pattern(dsp, 
 *                                                dof_handler.locally_owned_dofs(), 
 *                                                mpi_communicator, 
 *                                                locally_relevant_dofs); 
 * 
 *     system_matrix.reinit(locally_owned_dofs, 
 *                          locally_owned_dofs, 
 *                          dsp, 
 *                          mpi_communicator); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system</h4>
 * 

 * 
 * 然后组装线性系统的函数相对来说比较无聊，几乎和我们之前看到的一模一样。需要注意的地方是。
 * 

 * 
 * - 装配必须只在本地拥有的单元上循环。有多种方法来测试；例如，我们可以将一个单元的subdomain_id与三角形的信息进行比较，如<code>cell->subdomain_id() == triangulation.local_owned_subdomain()</code>，或者跳过所有条件<code>cell->is_ghost() || cell->is_artificial()</code>为真的单元。然而，最简单的方法是简单地询问单元格是否为本地处理器所拥有。
 * 

 * 
 * - 将本地贡献复制到全局矩阵中必须包括分配约束和边界值。换句话说，我们不能（就像我们在 step-6 中所做的那样）首先将每个本地贡献复制到全局矩阵中，然后在后面的步骤中才处理悬挂节点的约束和边界值。原因是，正如在 step-17 中所讨论的那样，一旦矩阵中的任意元素被组装到矩阵中，并行矢量类就不能提供对这些元素的访问--部分原因是它们可能不再存在于当前的处理器中，而是被运到了不同的机器上。
 * 

 * 
 * - 我们计算右手边的方式（考虑到介绍中的公式）可能不是最优雅的，但对于重点在某个完全不同的地方的程序来说是可以的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::assemble_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "assembly"); 
 * 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1); 
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
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           cell_matrix = 0.; 
 *           cell_rhs    = 0.; 
 * 
 *           fe_values.reinit(cell); 
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             { 
 *               const double rhs_value = 
 *                 (fe_values.quadrature_point(q_point)[1] > 
 *                      0.5 + 
 *                        0.25 * std::sin(4.0 * numbers::PI * 
 *                                        fe_values.quadrature_point(q_point)[0]) ? 
 *                    1. : 
 *                    -1.); 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     cell_matrix(i, j) += fe_values.shape_grad(i, q_point) * 
 *                                          fe_values.shape_grad(j, q_point) * 
 *                                          fe_values.JxW(q_point); 
 * 
 *                   cell_rhs(i) += rhs_value *                         // 
 *                                  fe_values.shape_value(i, q_point) * // 
 *                                  fe_values.JxW(q_point); 
 *                 } 
 *             } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           constraints.distribute_local_to_global(cell_matrix, 
 *                                                  cell_rhs, 
 *                                                  local_dof_indices, 
 *                                                  system_matrix, 
 *                                                  system_rhs); 
 *         } 
 * 
 * @endcode
 * 
 * 注意，上面的装配只是一个局部操作。因此，为了形成 "全局 "线性系统，需要在所有处理器之间进行同步。这可以通过调用函数compress()来实现。参见 @ref GlossCompress "压缩分布式对象"，以了解更多关于compress()的设计目的的信息。
 * 

 * 
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * 尽管在可能是数以万计的处理器上求解线性系统到目前为止并不是一项微不足道的工作，但完成这项工作的函数--至少在外表上--相对简单。大部分的部分你都见过了。真正值得一提的只有两件事。
 * 

 * 
 * - 解算器和预处理器是建立在PETSc和Trilinos功能的deal.II包装上的。众所周知，大规模并行线性求解器的主要瓶颈实际上不是处理器之间的通信，而是很难产生能够很好地扩展到大量处理器的预处理程序。在21世纪前十年的后半段，代数多网格（AMG）方法在这种情况下显然是非常有效的，我们将使用其中的一种方法--要么是可以通过PETSc接口的Hypre软件包的BoomerAMG实现，要么是由ML提供的预处理程序，它是Trilinos的一部分--用于当前的程序。解算器本身的其余部分是模板，之前已经展示过了。由于线性系统是对称和正定的，我们可以使用CG方法作为外解器。
 * 

 * 
 * - 最终，我们想要一个向量，它不仅存储了当前处理器拥有的自由度的解的元素，而且还存储了所有其他本地相关的自由度。另一方面，求解器本身需要一个在处理器之间唯一分割的向量，没有任何重叠。因此，我们在这个函数的开头创建一个具有这些特性的向量，用它来求解线性系统，并在最后才把它分配给我们想要的向量。这最后一步确保所有的鬼魂元素也在必要时被复制。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::solve() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "solve"); 
 *     LA::MPI::Vector    completely_distributed_solution(locally_owned_dofs, 
 *                                                     mpi_communicator); 
 * 
 *     SolverControl solver_control(dof_handler.n_dofs(), 1e-12); 
 * 
 * #ifdef USE_PETSC_LA 
 *     LA::SolverCG solver(solver_control, mpi_communicator); 
 * #else 
 *     LA::SolverCG solver(solver_control); 
 * #endif 
 * 
 *     LA::MPI::PreconditionAMG preconditioner; 
 * 
 *     LA::MPI::PreconditionAMG::AdditionalData data; 
 * 
 * #ifdef USE_PETSC_LA 
 *     data.symmetric_operator = true; 
 * #else 
 * /* Trilinos的默认值是好的  */ 
 * #endif 
 *     preconditioner.initialize(system_matrix, data); 
 * 
 *     solver.solve(system_matrix, 
 *                  completely_distributed_solution, 
 *                  system_rhs, 
 *                  preconditioner); 
 * 
 *     pcout << "   Solved in " << solver_control.last_step() << " iterations." 
 *           << std::endl; 
 * 
 *     constraints.distribute(completely_distributed_solution); 
 * 
 *     locally_relevant_solution = completely_distributed_solution; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrefine_grid"></a> 
 * <h4>LaplaceProblem::refine_grid</h4>
 * 

 * 
 * 估计误差和细化网格的函数又与  step-6  中的函数几乎完全一样。唯一不同的是，标志着要细化的单元格的函数现在在命名空间  parallel::distributed::GridRefinement  中 -- 这个命名空间的函数可以在所有参与的处理器之间进行通信，并确定全局阈值，用于决定哪些单元格要细化，哪些要粗化。
 * 

 * 
 * 注意，我们不需要对KellyErrorEstimator类做任何特殊处理：我们只是给它一个向量，其元素数量与本地三角形的单元（本地拥有的单元、幽灵单元和人工单元）一样多，但它只填入那些对应于本地拥有的单元的条目。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::refine_grid() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "refine"); 
 * 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(fe.degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       locally_relevant_solution, 
 *       estimated_error_per_cell); 
 *     parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
 *       triangulation, estimated_error_per_cell, 0.3, 0.03); 
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
 * 与 step-6 中的相应函数相比，这里的函数要复杂一点。有两个原因：第一个原因是，我们不只是想输出解决方案，还想输出每个单元的处理器（即它在哪个 "子域"）。其次，正如在 step-17 和 step-18 中详细讨论的那样，生成图形数据可能是并行化的一个瓶颈。在 step-18 中，我们将这一步骤从实际计算中移出，而是将其转移到一个单独的程序中，随后将各个处理器的输出合并到一个文件中。但这并不具规模：如果处理器的数量很大，这可能意味着在单个处理器上合并数据的步骤后来成为程序中运行时间最长的部分，或者它可能产生一个大到无法再可视化的文件。我们在这里遵循一个更合理的方法，即为每个MPI进程创建单独的文件，并将其留给可视化程序来理解。
 * 

 * 
 * 首先，函数的顶部看起来和平时一样。除了附加解决方案向量（包含所有本地相关元素的条目，而不仅仅是本地拥有的元素）外，我们还附加一个数据向量，为每个单元存储该单元所属的子域。这稍微有点棘手，因为当然不是每个处理器都知道每个单元。因此，我们附加的向量有一个当前处理器在其网格中拥有的每个单元的条目（本地拥有的单元、幽灵单元和人造单元），但DataOut类将忽略所有对应于不属于当前处理器的单元的条目。因此，我们在这些向量条目中写入什么值实际上并不重要：我们只需用当前MPI进程的编号（即当前进程的子域_id）来填充整个向量；这就正确地设置了我们关心的值，即对应于本地拥有的单元的条目，而为所有其他元素提供了错误的值--但无论如何这些都会被忽略。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(locally_relevant_solution, "u"); 
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells()); 
 *     for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *       subdomain(i) = triangulation.locally_owned_subdomain(); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 * 
 *     data_out.build_patches(); 
 * 
 * @endcode
 * 
 * 下一步是把这些数据写到磁盘上。在MPI-IO的帮助下，我们最多可以并行写入8个VTU文件。此外，还产生了一个PVTU记录，它将写入的VTU文件分组。
 * 

 * 
 * 
 * @code
 *     data_out.write_vtu_with_pvtu_record( 
 *       "./", "solution", cycle, mpi_communicator, 2, 8); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * 控制程序整体行为的函数又和  step-6  中的一样。小的区别是使用 <code>pcout</code> instead of <code>std::cout</code> 来输出到控制台（也见 step-17 ），而且我们只在最多涉及32个处理器的情况下产生图形输出。如果没有这个限制，人们很容易在没有阅读这个程序的情况下粗心大意地运行这个程序，从而导致集群互连中断，并填满任何可用的文件系统 :-)
 * 

 * 
 * 与 step-6 的一个功能上的区别是使用了一个正方形域，并且我们从一个稍细的网格开始（5个全局细化周期）--在4个单元上开始显示一个大规模的%并行程序没有什么意义（尽管承认在1024单元上开始显示的意义只是稍强）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::run() 
 *   { 
 *     pcout << "Running with " 
 * #ifdef USE_PETSC_LA 
 *           << "PETSc" 
 * #else 
 *           << "Trilinos" 
 * #endif 
 *           << " on " << Utilities::MPI::n_mpi_processes(mpi_communicator) 
 *           << " MPI rank(s)..." << std::endl; 
 * 
 *     const unsigned int n_cycles = 8; 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
 *       { 
 *         pcout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation); 
 *             triangulation.refine_global(5); 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 *         setup_system(); 
 * 
 *         pcout << "   Number of active cells:       " 
 *               << triangulation.n_global_active_cells() << std::endl 
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
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
 * } // namespace Step40 
 * 
 * @endcode
 * 
 * 
 * <a name="main"></a> 
 * <h4>main()</h4>
 * 

 * 
 * 最后一个函数，  <code>main()</code>  ，同样具有与所有其他程序相同的结构，特别是  step-6  。像其他使用MPI的程序一样，我们必须初始化和最终确定MPI，这是用辅助对象  Utilities::MPI::MPI_InitFinalize.  完成的。该类的构造函数也初始化了依赖MPI的库，如p4est、PETSc、SLEPc和Zoltan（尽管最后两个在本教程中没有使用）。这里的顺序很重要：在这些库被初始化之前，我们不能使用它们，所以在创建  Utilities::MPI::MPI_InitFinalize.  的实例之前做任何事情都没有意义。
 * 

 * 
 * 在求解器完成后，LaplaceProblem解构器将运行，然后是 Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize().  这个顺序也很重要： Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize() 调用 <code>PetscFinalize</code> （以及其他库的最终确定函数），这将删除任何正在使用的PETSc对象。这必须在我们解构拉普拉斯求解器之后进行，以避免双重删除错误。幸运的是，由于C++的析构器调用顺序规则，我们不需要担心这些：一切都以正确的顺序发生（即，与构造顺序相反）。由 Utilities::MPI::MPI_InitFinalize::~MPI_InitFinalize() 调用的最后一个函数是 <code>MPI_Finalize</code> ：也就是说，一旦这个对象被析构，程序应该退出，因为MPI将不再可用。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step40; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       LaplaceProblem<2> laplace_problem_2d; 
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
 * 
 *   return 0; 
 * } 
 * 
 * @endcode
examples/step-40/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你在单个处理器上或在几个本地MPI安装上运行该程序时，你应该得到这样的输出。

@code
Cycle 0:
   Number of active cells:       1024
   Number of degrees of freedom: 4225
   Solved in 10 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.176s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| assembly                        |         1 |    0.0209s |        12% |
| output                          |         1 |    0.0189s |        11% |
| setup                           |         1 |    0.0299s |        17% |
| solve                           |         1 |    0.0419s |        24% |
+---------------------------------+-----------+------------+------------+



Cycle 1:
   Number of active cells:       1954
   Number of degrees of freedom: 8399
   Solved in 10 iterations.



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |     0.327s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| assembly                        |         1 |    0.0368s |        11% |
| output                          |         1 |    0.0208s |       6.4% |
| refine                          |         1 |     0.157s |        48% |
| setup                           |         1 |    0.0452s |        14% |
| solve                           |         1 |    0.0668s |        20% |
+---------------------------------+-----------+------------+------------+



Cycle 2:
   Number of active cells:       3664
   Number of degrees of freedom: 16183
   Solved in 11 iterations.


...
@endcode



确切的数字是不同的，这取决于我们使用多少个处理器；这是由于预处理程序取决于问题的分区，然后解决方案在最后几位上有所不同，因此，网格细化也略有不同。不过，这里最值得注意的是，迭代次数并不随问题的大小而增加。这保证了我们甚至可以有效地解决最大的问题。

当在足够多的机器上运行时（比如说几千台），这个程序可以相对容易地在不到一分钟的时间内解决有远超过10亿个未知数的问题。另一方面，这样的大问题已经不能被视觉化，所以我们也只在16个处理器上运行该程序。下面是一个网格，以及它在16个处理器上的划分，还有相应的解决方案。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.mesh.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.solution.png" alt="">
</td>
</tr>
</table> 

左边的网格仅有7,069个单元。当然，这个问题我们在单台处理器上使用step-6就已经很容易解决了，但是这个程序的重点是展示如何编写一个可以扩展到更多机器的程序。例如，这里有两张图，显示了如果我们采取越来越多的处理器，程序的大量部分的运行时间是如何在大约5200万和37500万自由度的问题上扩展的（这些和接下来的几张图取自 @ref distributed_paper "分布式计算论文 "的早期版本；显示在更大数量的处理器上运行数据的更新图，以及更多的解释可以在该论文的最终版本中找到）。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.strong2.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.strong.png" alt="">
</td>
</tr>
</table> 

可以清楚地看到，这个程序可以很好地扩展到非常多的处理器。关于我们认为的 "可扩展 "程序的讨论，见 @ref GlossParallelScaling "本词汇表条目"）。曲线，特别是线性求解器，在图形的右端变得有点摇摆不定，因为每个处理器要做的事情太少，无法抵消通信成本（在上面两个例子中，每个处理器要解决的整个问题的部分，在使用4,096个处理器时，只有13,000和90,000个自由度；一个好的经验法则是，如果每个处理器至少有100,000个未知数，并行程序就会运行良好）。

虽然上面的强扩展图显示，如果我们采取越来越多的处理器，我们可以越来越快地解决一个固定大小的问题，但更有趣的问题可能是，问题可以变得多大，以便在一个特定大小的机器上仍然可以在合理的时间内解决它们。我们在下面两张256和4096处理器的图中展示了这一点。

 <table width="100%">
<tr>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.256.png" alt="">
</td>
<td>
  <img src="https://www.dealii.org/images/steps/developer/step-40.4096.png" alt="">
</td>
</tr>
</table> 

这些图显示的是，程序的所有部分都随着自由度数的增加而线性扩展。这一次，由于局部问题的规模太小，线条在左边摇摆不定。关于这些结果的更多讨论，我们参考了 @ref distributed_paper "分布式计算论文"。

那么，一个人能够解决的最大问题是多大？在写这个问题的时候，限制因素是程序使用<a
href="http://acts.nersc.gov/hypre/" target="_top">Hypre package</a>中的BoomerAMG代数多网格方法作为预处理程序，不幸的是，它使用有符号的32位整数来索引%分布式矩阵的元素。这将问题的大小限制在 $2^{31}-1=2,147,483,647$ 个自由度。从上面的图中可以看出，可扩展性会超过这个数字，而且可以预期，给定超过上面显示的4096台机器也会进一步减少计算时间。也就是说，人们当然可以期待，这个限制最终会被hybre的开发者解除。

另一方面，这并不意味着deal.II不能解决更大的问题。事实上，step-37展示了如何解决不仅仅是一点点，而是大大超过我们在这里所展示的任何问题的问题。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


从某种意义上说，这个程序是拉普拉斯方程的终极解算器：只要你有足够的处理器，它基本上可以把方程解到你想要的精度。由于拉普拉斯方程本身在这种精度水平上并不十分有趣，因此，更有趣的扩展可能性不在于这个程序，而在于它之后的内容。例如，本教程中的其他几个程序都有相当长的运行时间，特别是在3D中。因此，使用这里解释的技术来扩展其他程序以支持并行的分布式计算将是有趣的。我们在step-32教程程序中对step-31做了这样的处理，但同样的做法也适用于，例如，用于双曲时间相关问题的step-23和step-25，用于气体动力学的step-33，或用于纳维-斯托克斯方程的step-35。

也许同样有趣的是后处理的问题。如上所述，我们只展示了16个处理器的解决方案和网格的图片，因为4,096个处理器解决10亿个未知数会产生几10G的图形输出。目前，除非在至少几百个处理器上运行，否则没有任何程序能够以任何合理的方式将如此大量的数据可视化。然而，有一些方法，可视化程序直接与每个处理器上的求解器进行通信，每个可视化进程渲染这个处理器上的求解器所计算的场景部分。实现这样的接口将允许快速可视化那些在其他方面不适合用图形显示的东西。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-40.cc"
*/
