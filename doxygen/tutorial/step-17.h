/**
@page step_17 The step-17 tutorial program
This tutorial depends on step-8.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Overview">Overview</a>
        <li><a href="#ParallelizingsoftwarewithMPI">Parallelizing software with MPI</a>
        <li><a href="#Whatthisprogramdoes">What this program does</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeElasticProblemcodeclasstemplate">The <code>ElasticProblem</code> class template</a>
        <li><a href="#Righthandsidevalues">Right hand side values</a>
        <li><a href="#ThecodeElasticProblemcodeclassimplementation">The <code>ElasticProblem</code> class implementation</a>
      <ul>
        <li><a href="#ElasticProblemElasticProblem">ElasticProblem::ElasticProblem</a>
        <li><a href="#ElasticProblemsetup_system">ElasticProblem::setup_system</a>
        <li><a href="#ElasticProblemassemble_system">ElasticProblem::assemble_system</a>
        <li><a href="#ElasticProblemsolve">ElasticProblem::solve</a>
        <li><a href="#ElasticProblemrefine_grid">ElasticProblem::refine_grid</a>
        <li><a href="#ElasticProblemoutput_results">ElasticProblem::output_results</a>
        <li><a href="#ElasticProblemrun">ElasticProblem::run</a>
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
examples/step-17/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Overview"></a><h3>Overview</h3>


这个程序没有引入任何新的数学思想；事实上，它所做的只是做与step-8已经做的完全相同的计算，但它以一种不同的方式来做：我们没有使用deal.II自己的线性代数类，而是在deal.II提供的类之上建立一切，这些类包裹着<a
href="http://www.mcs.anl.gov/petsc/" target="_top">PETSc</a>库的线性代数实现。由于PETSc允许将矩阵和向量分布在MPI网络中的几台计算机上，因此产生的代码甚至能够以%并行方式解决问题。如果你不知道PETSc是什么，那么这将是一个快速浏览其主页的好时机。

作为这个程序的先决条件，你需要安装PETSc，如果你想在一个集群上以%并行方式运行，你还需要<a
href="http://www-users.cs.umn.edu/~karypis/metis/index.html"
target="_top">METIS</a>来划分网格。在<a
href="../../readme.html" target="body">README</a>文件中描述了deal.II和这两个附加库的安装。

现在，关于细节：如前所述，该程序不计算任何新的东西，所以对有限元类等的使用与以前完全相同。与以前的程序不同的是，我们用几乎所有的类 <code>Vector</code> and <code>SparseMatrix</code> 代替了它们的近似值 <code>PETScWrappers::MPI::Vector</code> 和 <code>PETScWrappers::MPI::SparseMatrix</code> ，它们存储数据的方式使MPI网络中的每个处理器只存储矩阵或矢量的一部分。更具体地说，每个处理器将只存储与它 "拥有 "的自由度相对应的矩阵的那些行。对于向量，它们要么只存储与处理器拥有的自由度相对应的元素（这是右手边所需要的），要么也存储一些额外的元素，以确保每个处理器都能访问处理器拥有的单元（所谓 @ref GlossLocallyActiveDof "本地活动的自由度"）或邻近单元（所谓 @ref GlossLocallyRelevantDof "本地相关自由度"）上的解组件。

来自PETScWrapper命名空间的类所提供的接口与deal.II线性代数类的接口非常相似，但它们不是自己实现这一功能，而是简单地传递给它们相应的PETSc函数。因此，包装器只是用来给PETSc一个更现代的、面向对象的接口，并使PETSc和deal.II对象的使用尽可能地互换。使用PETSc的主要意义在于它可以在%并行状态下运行。我们将利用这一点，将域划分为与MPI网络中的进程一样多的块（"子域"）。同时，PETSc还提供了假的MPI存根，所以如果PETSc的配置中没有MPI，你可以在一台机器上运行这个程序。




<a name="ParallelizingsoftwarewithMPI"></a><h3>Parallelizing software with MPI</h3>


开发软件以通过MPI在%parallel中运行，需要改变一下思维方式，因为我们通常必须分割所有的数据结构，使每个处理器只存储整个问题的一部分。因此，你通常不能在每个处理器上访问一个解决方案向量的所有组成部分 -- 每个处理器可能根本没有足够的内存来容纳整个解决方案向量。由于数据被分割或 "分布 "在各个处理器上，我们把MPI使用的编程模型称为 "分布式内存计算"（与 "共享内存计算 "相反，后者意味着多个处理器可以访问一个内存空间中的所有数据，例如，当一台机器的多个核心在一个共同任务上工作时）。分布式内存计算的一些基本原理在 @ref distributed "使用分布式内存的多处理器并行计算 "文档模块中讨论，该模块本身是 @ref Parallel "并行计算 "模块的一个子模块。

一般来说，为了真正能够扩展到大量的处理器，我们需要在可用的处理器之间分割出<i>every</i>数据结构，其大小随着整个问题的大小而扩展。关于程序 "扩展 "的定义，见 @ref GlossParallelScaling "本词汇表条目"）。这包括，例如，三角形、矩阵和所有全局向量（解决方案，右手边）。如果不拆分所有这些对象，其中一个对象将被复制到所有的处理器上，如果问题大小（和可用的处理器数量）变得很大，最终会简单地变得太大。另一方面，在每个处理器上保留大小与整个问题大小无关的对象是完全可以的。例如，可执行文件的每个副本将创建自己的有限元对象，或者我们在汇编中使用的局部矩阵）。)

在当前的程序中（以及相关的第18步），我们不会走得这么远，而是对MPI的使用做一个比较温和的介绍。更具体地说，我们要并行化的数据结构只有矩阵和向量。然而，我们并没有拆分Triangulation和DoFHandler类：每个进程仍然拥有这些对象的完整副本，而且所有进程都拥有其他进程所拥有的确切副本。然后，我们只需在每个处理器上的三角形的每个副本中，标记哪个处理器拥有哪些单元。这个过程被称为将网格 "分割 "为 @ref GlossSubdomainId "子域"。

对于较大的问题，必须在每个处理器上存储<i>entire</i>网格，显然会产生一个瓶颈。分割网格是稍微的，虽然没有多复杂（从用户的角度来看，虽然它<i>much</i>下更复杂）来实现，我们将展示如何在step-40和其他一些程序中这样做。在讨论这个程序的某个功能如何工作的过程中，我们会多次评论它不会扩展到大型问题，以及为什么不会。所有这些问题都将在第18步，特别是第40步中得到解决，它可以扩展到非常多的进程。

从哲学上讲，MPI的运作方式如下。你通常通过以下方式运行一个程序

@code
  mpirun -np 32 ./step-17
@endcode

这意味着在（比如）32个处理器上运行它。如果你是在一个集群系统上，你通常需要<i>schedule</i>程序在32个处理器可用时运行；这将在你的集群的文档中描述。但是在系统内部，每当这些处理器可用时，通常会执行上述相同的调用）。)这样做的目的是，MPI系统将启动32个<i>copies</i>的 <code>step-17</code> 的可执行文件。(这些正在运行的可执行文件中的每一个的MPI术语是，你有32个 @ref GlossMPIProcess "MPI进程"。)这可能发生在不同的机器上，甚至不能从对方的内存空间中读取，也可能发生在同一台机器上，但最终的结果是一样的：这32个副本中的每一个都将以操作系统分配给它的一些内存运行，而且它不能直接读取其他31个副本的内存。为了在一个共同的任务中进行协作，这32个副本就必须<i>communicate</i>相互协作。MPI是<i>Message Passing Interface</i>的缩写，通过允许程序<i>send messages</i>来实现这一点。你可以把它看作是邮件服务：你可以把一封写给特定地址的信放入邮件，它将被送达。但这是你能控制事物的程度。如果你想让收信人对信的内容做些什么，例如把你想要的数据从那边返回给你，那么需要发生两件事。(i)接收方需要实际去检查他们的邮箱里是否有东西，(ii)如果有的话，做出适当的反应，比如说发送数据回来。如果你等待这个返回信息，但原来的接收者却心不在焉，没有注意，那么你就不走运了：你只需要等待，直到你在那边的请求将被解决。在某些情况下，错误会导致原始接收者永远不检查你的邮件，在这种情况下，你将永远等待--这被称为<i>deadlock</i>。(  @dealiiVideoLectureSeeAlso{39,41,41.25,41.5}) 

在实践中，人们通常不在发送和接收单个消息的层面上编程，而是使用更高层次的操作。例如，在程序中，我们将使用函数调用，从每个处理器获取一个数字，将它们全部相加，然后将总和返回给所有处理器。在内部，这是用单个消息实现的，但对用户来说，这是透明的。我们称这种操作为<i>collectives</i>，因为<i>all</i>处理器参与其中。集合体允许我们编写程序，其中不是每个可执行文件的副本都在做完全不同的事情（这将是难以置信的编程难度），但实质上所有副本都在为自己做同样的事情（尽管是在不同的数据上），通过相同的代码块运行；然后他们通过集合体进行数据通信；然后再回到为自己做事情，通过相同的数据块运行。这是能够编写程序的关键部分，也是确保程序能够在任何数量的处理器上运行的关键部分，因为我们不需要为每个参与的处理器编写不同的代码。

这并不是说程序从来都是以不同的处理器在其可执行文件的副本中运行不同的代码块的方式来编写的。程序内部也经常以其他方式而不是通过集合体进行通信。但是在实践中，%并行有限元代码几乎总是遵循这样的方案：程序的每个副本在同一时间运行相同的代码块，中间穿插着所有处理器相互交流的阶段）。)

在现实中，即使是调用MPI集体函数的水平也太低了。相反，下面的程序根本不会包含对MPI的任何直接调用，而只包含对deal.II的用户隐藏这种通信的函数。这样做的好处是，你不需要学习MPI的细节和相当复杂的函数调用。也就是说，你确实必须理解上文所述的MPI背后的一般哲学。




<a name="Whatthisprogramdoes"></a><h3>What this program does</h3>


然后，这个程序演示的技术是。

- 如何使用PETSc封装类；这在本程序的主类的声明中已经可以看到，  <code>ElasticProblem</code>  。

- 如何将网格划分为子域；这发生在 <code>ElasticProblem::setup_system()</code> 函数。

- 如何对运行在MPI网络上的作业进行并行化操作；在这里，这是在很多地方都要注意的，最明显的是在 <code>ElasticProblem::assemble_system()</code> 函数中。

- 如何处理只存储向量项子集的向量，对于这些向量，我们必须确保它们在当前处理器上存储我们需要的东西。例如见 <code>ElasticProblem::solve()</code> and <code>ElasticProblem::refine_grid()</code> 函数。

- 如何处理同时在多个处理器上运行的程序的状态输出。这是通过程序中的 <code>pcout</code> 变量完成的，在构造函数中初始化。

由于这一切只能用实际的代码来证明，让我们直接进入代码，不再多说。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 *  

 * 
 * 
 * @code
 *  * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004 
 *  *         Wolfgang Bangerth, Texas A&M University, 2016 
 *  */ 
 * 
 * 
 * @endcode
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * 首先是我们在以前的例子程序中已经使用过的常见的各种头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/multithread_info.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/sparsity_tools.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 这里是我们对这个例子程序特别需要的东西，而这些东西并不在  step-8  中。首先，我们替换掉标准输出 <code>std::cout</code> by a new stream <code>pcout</code> ，它在并行计算中只用于在其中一个MPI进程中生成输出。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * 
 * @endcode
 * 
 * 我们将通过调用  Utilities::MPI  名称空间中的相应函数来查询进程的数量和当前进程的数量。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/mpi.h> 
 * 
 * @endcode
 * 
 * 然后，我们要把所有涉及（全局）线性系统的线性代数组件替换成类，这些类围绕PETSc提供的接口与我们自己的线性代数类相似（PETSc是一个用C语言编写的库，而deal.II附带的包装类提供的PETSc功能的接口与我们自己的线性代数类已经有的接口相似）。特别是，我们需要在MPI程序中分布在几个 @ref GlossMPIProcess "进程 "中的向量和矩阵（如果只有一个进程，也就是说，如果你只在一台机器上运行，并且没有MPI支持，则简单映射为顺序的、本地的向量和矩阵）。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/petsc_vector.h> 
 * #include <deal.II/lac/petsc_sparse_matrix.h> 
 * 
 * @endcode
 * 
 * 然后，我们还需要PETSc提供的求解器和预处理器的接口。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/petsc_solver.h> 
 * #include <deal.II/lac/petsc_precondition.h> 
 * 
 * @endcode
 * 
 * 此外，我们还需要一些划分网格的算法，以便在MPI网络上有效地分布这些网格。分区算法在 <code>GridTools</code> 命名空间中实现，我们需要一个额外的包含文件，用于 <code>DoFRenumbering</code> 中的一个函数，该函数允许对与自由度相关的索引进行排序，以便根据它们所关联的子域进行编号。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * 
 * @endcode
 * 
 * 而这又是简单的C++。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step17 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclasstemplate"></a> 
 * <h3>The <code>ElasticProblem</code> class template</h3>
 * 

 * 
 * 该程序的第一个真正的部分是主类的声明。 正如在介绍中提到的，几乎所有的内容都是从 step-8 中逐字复制过来的，所以我们只对这两个教程之间的少数差异进行评论。 有一个（表面上的）变化是，我们让 <code>solve</code> 返回一个值，即收敛所需的迭代次数，这样我们就可以在适当的地方将其输出到屏幕上。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ElasticProblem 
 *   { 
 *   public: 
 *     ElasticProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void         setup_system(); 
 *     void         assemble_system(); 
 *     unsigned int solve(); 
 *     void         refine_grid(); 
 *     void         output_results(const unsigned int cycle) const; 
 * 
 * @endcode
 * 
 * 第一个变化是，我们必须声明一个变量，表明我们应该通过它来分配我们的计算的 @ref GlossMPICommunicator "MPI通信器"。
 * 

 * 
 * 
 * @code
 *     MPI_Comm mpi_communicator; 
 * 
 * @endcode
 * 
 * 然后我们有两个变量，告诉我们在并行世界中的位置。下面的第一个变量， <code>n_mpi_processes</code>  ，告诉我们总共有多少个MPI进程，而第二个变量， <code>this_mpi_process</code>  ，表示在这个进程空间中，目前进程的编号（在MPI语言中，这相当于进程的 @ref GlossMPIRank  "等级"）。后者对每个进程都有一个唯一的值，介于0和（小于） <code>n_mpi_processes</code>  之间。如果这个程序运行在没有MPI支持的单机上，那么它们的值分别为 <code>1</code> and <code>0</code>  ，。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_mpi_processes; 
 *     const unsigned int this_mpi_process; 
 * 
 * @endcode
 * 
 * 接下来是一个类似流的变量  <code>pcout</code>  。从本质上讲，它只是我们为了方便而使用的东西：在一个并行程序中，如果每个进程都输出状态信息，那么很快就会有很多杂乱的信息。相反，我们希望只让一个 @ref GlossMPIProcess "进程 "输出一次所有的信息，例如， @ref GlossMPIRank "等级 "为零的那个。同时，在我们创建输出的<i>every</i>地方加上 <code>if (my_rank==0)</code> 条件的前缀似乎很傻。
 * 

 * 
 * 为了使这个问题更简单，ConditionalOStream类正是这样做的：它就像一个流一样，但只有在一个标志被设置后才转发到一个真正的、底层的流。通过将这个条件设置为 <code>this_mpi_process==0</code> （其中 <code>this_mpi_process</code> 对应于MPI进程的等级），我们确保输出只从第一个进程中产生，并且我们不会在每个进程中重复得到同样的输出行。因此，我们可以在每一个地方和每一个进程中使用 <code>pcout</code> ，但是除了一个进程之外，所有的进程都不会发生通过 <code>operator&lt;&lt;</code> 输送到对象中的信息。
 * 

 * 
 * 
 * @code
 *     ConditionalOStream pcout; 
 * 
 * @endcode
 * 
 * 成员变量列表的其余部分与  step-8  中的内容基本相同。然而，我们改变了矩阵和矢量类型的声明，以使用并行的PETSc对象代替。请注意，我们没有使用单独的稀疏模式，因为PETSc将其作为矩阵数据结构的一部分进行内部管理。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     AffineConstraints<double> hanging_node_constraints; 
 * 
 *     PETScWrappers::MPI::SparseMatrix system_matrix; 
 * 
 *     PETScWrappers::MPI::Vector solution; 
 *     PETScWrappers::MPI::Vector system_rhs; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Righthandsidevalues"></a> 
 * <h3>Right hand side values</h3>
 * 

 * 
 * 以下内容取自 step-8 ，未作改动。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  values) const override 
 *     { 
 *       Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 
 *       Assert(dim >= 2, ExcInternalError()); 
 * 
 *       Point<dim> point_1, point_2; 
 *       point_1(0) = 0.5; 
 *       point_2(0) = -0.5; 
 * 
 *       if (((p - point_1).norm_square() < 0.2 * 0.2) || 
 *           ((p - point_2).norm_square() < 0.2 * 0.2)) 
 *         values(0) = 1; 
 *       else 
 *         values(0) = 0; 
 * 
 *       if (p.square() < 0.2 * 0.2) 
 *         values(1) = 1; 
 *       else 
 *         values(1) = 0; 
 *     } 
 * 
 *     virtual void 
 *     vector_value_list(const std::vector<Point<dim>> &points, 
 *                       std::vector<Vector<double>> &  value_list) const override 
 *     { 
 *       const unsigned int n_points = points.size(); 
 * 
 *       Assert(value_list.size() == n_points, 
 *              ExcDimensionMismatch(value_list.size(), n_points)); 
 * 
 *       for (unsigned int p = 0; p < n_points; ++p) 
 *         RightHandSide<dim>::vector_value(points[p], value_list[p]); 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclassimplementation"></a> 
 * <h3>The <code>ElasticProblem</code> class implementation</h3>
 * 
 * <a name="ElasticProblemElasticProblem"></a> 
 * <h4>ElasticProblem::ElasticProblem</h4>
 * 

 * 
 * 实际实现的第一步是主类的构造函数。除了初始化我们在 step-8 中已经有的相同成员变量外，我们在这里用连接所有进程的全局MPI通信器来初始化我们将使用的MPI通信器变量（在更复杂的应用中，可以在这里使用只连接所有进程的一个子集的通信器对象），并调用 Utilities::MPI 辅助函数来确定进程的数量以及当前进程在这个画面中的地位。此外，我们确保输出只由（全局）第一个进程产生。我们通过将我们想要输出的流传给 (<code>std::cout</code>) 和一个真/假标志作为参数，后者是通过测试当前执行构造函数调用的进程是否是MPI宇宙中的第一个来确定的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ElasticProblem<dim>::ElasticProblem() 
 *     : mpi_communicator(MPI_COMM_WORLD) 
 *     , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)) 
 *     , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)) 
 *     , pcout(std::cout, (this_mpi_process == 0)) 
 *     , fe(FE_Q<dim>(1), dim) 
 *     , dof_handler(triangulation) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsetup_system"></a> 
 * <h4>ElasticProblem::setup_system</h4>
 * 

 * 
 * 接下来，我们需要实现为要解决的全局线性系统设置各种变量的函数。
 * 

 * 
 * 然而，在我们进行这项工作之前，对于一个并行程序来说，有一件事要做：我们需要确定哪个MPI进程负责每个单元。在进程之间分割单元，通常称为 "划分网格"，是通过给每个单元分配一个 @ref GlossSubdomainId "子域id "来完成的。我们通过调用METIS库来完成这一工作，METIS库以一种非常有效的方式完成这一工作，试图将子域之间接口上的节点数量降到最低。我们没有尝试直接调用METIS，而是通过调用 GridTools::partition_triangulation() 函数来实现，该函数在更高的编程水平上实现了这一点。
 * 

 * 
 * @note  正如在介绍中提到的，如果我们使用 parallel::shared::Triangulation 类来代替三角形对象，我们就可以避免这个手动划分的步骤（正如我们在 step-18 中所做的）。  该类实质上做了所有常规三角形的工作，但它也在每次创建或细化网格操作后自动划分网格。
 * 

 * 
 * 在分割之后，我们需要像往常一样列举所有的自由度。 然而，我们希望列举自由度的方式是：所有与子域0（位于进程0）的单元相关的自由度都在与子域1的单元相关的自由度之前，在进程2的单元之前，以此类推。我们需要这样做，因为我们必须将全局向量的右手边和解决方案，以及矩阵分割成连续的行块，住在每个处理器上，而且我们希望以一种需要最小通信的方式来做。这个特殊的列举可以通过使用 DoFRenumbering::subdomain_wise(). 对自由度指数重新排序来获得。
 * 

 * 
 * 这个初始设置的最后一步是，我们为自己得到一个IndexSet，表示这个过程所负责的全局未知数的子集。(注意，一个自由度不一定是由拥有一个单元的进程所拥有，只是因为这个自由度生活在这个单元上：有些自由度生活在子域之间的接口上，因此只由这个接口附近的一个进程所拥有。)
 * 

 * 
 * 在我们继续之前，让我们回顾一下在介绍中已经讨论过的一个事实。我们在这里使用的三角形是在所有进程中复制的，每个进程都有整个三角形的完整副本，包括所有单元。分区只提供了一种方法来确定每个进程 "拥有 "哪些单元，但它知道所有单元的一切。同样，DoFHandler对象知道每个单元的一切，特别是每个单元上的自由度，无论它是否是当前进程拥有的单元。这不能扩展到大型问题，因为如果问题足够大，最终只是在每个进程中存储整个网格以及与之相关的所有内容将变得不可行。另一方面，如果我们将三角形分割成若干部分，使每个进程只存储它 "拥有 "的单元格，而不存储其他的单元格（或者，至少是其他单元格的一小部分），那么，只要我们将足够多的MPI进程扔给它们，我们就可以解决大问题。这就是我们在 step-40 中要做的，例如，使用 parallel::distributed::Triangulation 类。 另一方面，我们在当前程序中演示的其余大部分内容实际上将继续工作，无论我们有整个三角形的可用，还是只有其中的一部分。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::setup_system() 
 *   { 
 *     GridTools::partition_triangulation(n_mpi_processes, triangulation); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 *     DoFRenumbering::subdomain_wise(dof_handler); 
 * 
 * @endcode
 * 
 * 我们需要初始化表示当前网格的悬挂节点约束的对象。与三角形和DoFHandler对象一样，我们将简单地在每个进程上存储<i>all</i>约束；同样，这不会有规模，但我们在 step-40 中展示了如何通过在每个MPI进程上只存储对这个特定进程实际重要的自由度约束来解决这个问题。
 * 

 * 
 * 
 * @code
 *     hanging_node_constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                             hanging_node_constraints); 
 *     hanging_node_constraints.close(); 
 * 
 * @endcode
 * 
 * 现在我们为系统矩阵创建稀疏性模式。请注意，我们再次计算并存储所有条目，而不仅仅是与此相关的条目（参见 step-18 或 step-40 ，以获得更有效的处理方式）。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, 
 *                                     dsp, 
 *                                     hanging_node_constraints, 
 *                                     false); 
 * 
 * @endcode
 * 
 * 现在我们确定本地拥有的DoF的集合，并使用它来初始化并行向量和矩阵。由于矩阵和向量需要并行工作，我们必须向它们传递一个MPI通信对象，以及IndexSet  @p locally_owned_dofs. 中包含的分区信息。IndexSet包含关于全局大小（<i>total</i>自由度数）的信息，也包含要在本地存储哪些行的子集。 注意，系统矩阵需要该行和列的分区信息。对于正方形矩阵，就像这里的情况一样，列的划分方式应该与行的划分方式相同，但是对于矩形矩阵，我们必须按照与矩阵相乘的向量的划分方式来划分列，而行的划分方式必须与矩阵-向量乘法的目的向量相同。
 * 

 * 
 * 
 * @code
 *     const std::vector<IndexSet> locally_owned_dofs_per_proc = 
 *       DoFTools::locally_owned_dofs_per_subdomain(dof_handler); 
 *     const IndexSet locally_owned_dofs = 
 *       locally_owned_dofs_per_proc[this_mpi_process]; 
 * 
 *     system_matrix.reinit(locally_owned_dofs, 
 *                          locally_owned_dofs, 
 *                          dsp, 
 *                          mpi_communicator); 
 * 
 *     solution.reinit(locally_owned_dofs, mpi_communicator); 
 *     system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemassemble_system"></a> 
 * <h4>ElasticProblem::assemble_system</h4>
 * 

 * 
 * 我们现在组装矩阵和问题的右手边。在我们进行详细讨论之前，有一些事情值得一提。首先，我们将并行组装系统，也就是说，每个进程将负责在属于这个特定进程的单元上进行组装。请注意，自由度的分割方式是，单元内部和属于同一子域的单元之间的所有自由度都属于 <code>owns</code> 该单元的过程。然而，即使如此，我们有时也需要在一个单元上与属于不同过程的邻居集合，在这些情况下，当我们将局部贡献加到全局矩阵或右手向量中时，我们必须将这些条目转移到拥有这些元素的过程中。幸运的是，我们不需要用手去做这件事。PETSc为我们做了这一切，它在本地缓存了这些元素，当我们在这个函数的末尾对矩阵和向量调用 <code>compress()</code> 函数时，根据需要将它们发送给其他进程。
 * 

 * 
 * 第二点是，一旦我们把矩阵和向量的贡献交给了PETSc，那么，a）很难，b）要把它们拿回来进行修改，效率非常低。这不仅是PETSc的错，也是这个程序的分布式性质的结果：如果一个条目驻留在另一个处理器上，那么要得到它必然是很昂贵的。这样做的后果是，我们不应该试图首先组装矩阵和右手边，就像没有悬挂的节点约束和边界值一样，然后在第二步中消除这些约束（例如使用 AffineConstraints::condense()). ），相反，我们应该在将这些条目交给PETSc之前尝试消除悬挂的节点约束。这很容易：我们不需要手工复制元素到全局矩阵中（就像我们在 step-4 中做的那样），而是使用 AffineConstraints::distribute_local_to_global() 函数来同时处理悬空节点的问题。我们在  step-6  中也已经这样做了。第二步，消除边界节点，也可以这样做，把边界值放到与悬挂节点相同的AffineConstraints对象中（例如，见 step-6 中的方法）；但是，严格来说，在这里没有必要这样做，因为消除边界值可以只用每个进程本身存储的数据来完成，因此，我们使用之前在 step-4 中使用的方法，即通过 MatrixTools::apply_boundary_values().  
 * 

 * 
 * 说了这么多，下面是实际的实现，从辅助变量的一般设置开始。 请注意，我们仍然使用deal.II的全矩阵和向量类型的本地系统，因为这些类型很小，不需要在不同进程中共享）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::assemble_system() 
 *   { 
 *     QGauss<dim>   quadrature_formula(fe.degree + 1); 
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
 *     std::vector<double> lambda_values(n_q_points); 
 *     std::vector<double> mu_values(n_q_points); 
 * 
 *     Functions::ConstantFunction<dim> lambda(1.), mu(1.); 
 * 
 *     RightHandSide<dim>          right_hand_side; 
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim)); 
 * 
 * @endcode
 * 
 * 接下来是对所有元素的循环。请注意，我们不需要在每个进程上做<i>all</i>的工作：我们在这里的工作只是在实际属于这个MPI进程的单元上组装系统，所有其他的单元将由其他进程来处理。这就是紧随for-loop之后的if-clause所要处理的：它查询每个单元的子域标识符，这是一个与每个单元相关的数字，告诉我们所有者进程的情况。在更大的范围内，子域标识被用来将一个域分成几个部分（我们在上面 <code>setup_system()</code> 的开头就这样做了），并允许识别一个单元生活在哪个子域。在这个应用中，我们让每个进程恰好处理一个子域，所以我们确定了  <code>subdomain</code> and <code>MPI process</code>  的条款。
 * 

 * 
 * 除此以外，如果你已经了解了  step-8  中的组装方式，那么组装本地系统就相对不容易了。如上所述，将本地贡献分配到全局矩阵和右手边，也是以与  step-6  中相同的方式来处理悬挂节点约束。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->subdomain_id() == this_mpi_process) 
 *         { 
 *           cell_matrix = 0; 
 *           cell_rhs    = 0; 
 * 
 *           fe_values.reinit(cell); 
 * 
 *           lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
 *           mu.value_list(fe_values.get_quadrature_points(), mu_values); 
 * 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const unsigned int component_i = 
 *                 fe.system_to_component_index(i).first; 
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const unsigned int component_j = 
 *                     fe.system_to_component_index(j).first; 
 * 
 *                   for (unsigned int q_point = 0; q_point < n_q_points; 
 *                        ++q_point) 
 *                     { 
 *                       cell_matrix(i, j) += 
 *                         ((fe_values.shape_grad(i, q_point)[component_i] * 
 *                           fe_values.shape_grad(j, q_point)[component_j] * 
 *                           lambda_values[q_point]) + 
 *                          (fe_values.shape_grad(i, q_point)[component_j] * 
 *                           fe_values.shape_grad(j, q_point)[component_i] * 
 *                           mu_values[q_point]) + 
 *                          ((component_i == component_j) ? 
 *                             (fe_values.shape_grad(i, q_point) * 
 *                              fe_values.shape_grad(j, q_point) * 
 *                              mu_values[q_point]) : 
 *                             0)) * 
 *                         fe_values.JxW(q_point); 
 *                     } 
 *                 } 
 *             } 
 * 
 *           right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
 *                                             rhs_values); 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const unsigned int component_i = 
 *                 fe.system_to_component_index(i).first; 
 * 
 *               for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *                 cell_rhs(i) += fe_values.shape_value(i, q_point) * 
 *                                rhs_values[q_point](component_i) * 
 *                                fe_values.JxW(q_point); 
 *             } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           hanging_node_constraints.distribute_local_to_global(cell_matrix, 
 *                                                               cell_rhs, 
 *                                                               local_dof_indices, 
 *                                                               system_matrix, 
 *                                                               system_rhs); 
 *         } 
 * 
 * @endcode
 * 
 * 下一步是对向量和系统矩阵进行 "压缩"。这意味着每个进程将对矩阵和向量中那些自己不拥有的条目所做的添加发送给拥有这些条目的进程。在收到其他进程的这些加法后，每个进程再把它们加到它已经拥有的值上。这些加法是将生活在几个单元上的形状函数的积分贡献结合起来，就像在串行计算中一样，不同的是这些单元被分配给不同的进程。
 * 

 * 
 * 
 * @code
 *     system_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 * 
 * @endcode
 * 
 * 全局矩阵和右边的向量现在已经形成。我们仍然要应用边界值，方法与我们在 step-3 ,  step-4 , 和其他一些程序中的方法相同。
 * 

 * 
 * 下面调用 MatrixTools::apply_boundary_values() 的最后一个参数允许进行一些优化。它控制我们是否应该删除对应于边界节点的矩阵列中的条目（即，将其设置为零），或者保留它们（通过 <code>true</code> 意味着：是的，消除这些列）。如果我们消除了列，那么结果矩阵将再次成为对称的，如果我们不这样做，那么它将不会。不过，结果系统的解应该是一样的。我们想让系统重新成为对称的唯一原因是我们想使用CG方法，该方法只对对称矩阵有效。我们可能<i>not</i>想让矩阵对称的原因是，这将要求我们写进实际存在于其他进程中的列项，即涉及到数据的交流。这总是很昂贵的。
 * 

 * 
 * 经验告诉我们，如果我们不删除与边界节点相关的列，CG也可以工作（而且工作得几乎一样好），这可以用这种特殊的非对称性结构来解释。为了避免通信的费用，我们因此不消除受影响列中的条目。
 * 

 * 
 * 
 * @code
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              Functions::ZeroFunction<dim>(dim), 
 *                                              boundary_values); 
 *     MatrixTools::apply_boundary_values( 
 *       boundary_values, system_matrix, solution, system_rhs, false); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsolve"></a> 
 * <h4>ElasticProblem::solve</h4>
 * 

 * 
 * 组建了线性系统后，我们接下来需要解决它。PETSc提供了各种顺序和并行求解器，我们为这些求解器编写了包装器，其接口与之前所有示例程序中使用的deal.II求解器几乎相同。因此，下面的代码看起来应该相当熟悉。
 * 

 * 
 * 在该函数的顶部，我们设置了一个收敛监视器，并指定了我们希望解决线性系统的精度。接下来，我们使用PETSc的CG求解器创建一个实际的求解器对象，该求解器也可用于并行（分布式）矢量和矩阵。最后是一个预处理程序；我们选择使用一个块状雅可比预处理程序，它通过计算矩阵的每个对角线块的不完全LU分解来工作。 换句话说，每个MPI进程从其存储的行中计算出一个ILU，丢掉与本地未存储的行指数相对应的列；这就产生了一个方形的矩阵块，我们可以从中计算出一个ILU。这意味着如果你只用一个进程来运行程序，那么你将使用一个ILU(0)作为预处理程序，而如果它在许多进程上运行，那么我们将在对角线上有许多块，预处理程序是这些块中每个块的ILU(0)。在每个处理器只有一个自由度的极端情况下，这个预处理程序只是一个雅可比预处理程序，因为对角线矩阵块只由一个条目组成。这样的预处理程序相对容易计算，因为它不需要在处理器之间进行任何形式的通信，但一般来说，对于大量的处理器来说，它的效率并不高)。
 * 

 * 
 * 按照这样的设置，我们就可以解决这个线性系统。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   unsigned int ElasticProblem<dim>::solve() 
 *   { 
 *     SolverControl solver_control(solution.size(), 1e-8 * system_rhs.l2_norm()); 
 *     PETScWrappers::SolverCG cg(solver_control, mpi_communicator); 
 * 
 *     PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 * @endcode
 * 
 * 下一步是分配悬挂的节点约束。这有点麻烦，因为要填入一个约束节点的值，你需要访问它所约束的节点的值（例如，对于2d中的Q1元素，我们需要访问悬挂节点面的大边上的两个节点，以计算中间的约束节点的值）。
 * 

 * 
 * 问题是，我们已经建立了我们的向量（在 <code>setup_system()</code> 中），使每个进程只负责存储解向量中与该进程 "拥有 "的自由度相对应的那些元素。然而，在有些情况下，为了计算一个进程中受限自由度的向量项的值，我们需要访问存储在其他进程中的向量项。 PETSc（以及它所基于的MPI模型）不允许简单地查询存储在其他进程上的向量条目，所以我们在这里所做的是获得一个 "分布式 "向量的副本，我们将所有元素存储在本地。这很简单，因为deal.II包装器有一个针对deal.II Vector类的转换构造函数。这种转换当然需要通信，但实质上每个进程只需要将其数据批量发送给其他每个进程一次，而不需要对单个元素的查询做出回应）。
 * 

 * 
 * 
 * @code
 *     Vector<double> localized_solution(solution); 
 * 
 * @endcode
 * 
 * 当然，和以前的讨论一样，如果你想在大量进程上解决大问题，这样的步骤显然不能扩展得很远，因为现在每个进程都存储了<i>all elements</i>的解向量。(我们将在 step-40 中展示如何更好地做到这一点。) 另一方面，在这个本地副本上分配悬挂节点约束很简单，使用通常的函数 AffineConstraints::distributed().  特别是，我们可以计算<i>all</i>约束自由度的值，无论当前进程是否拥有它们。
 * 

 * 
 * 
 * @code
 *     hanging_node_constraints.distribute(localized_solution); 
 * 
 * @endcode
 * 
 * 然后把所有的东西都转回全局向量中。下面的操作是复制我们在分布式解决方案中本地存储的那些本地化解决方案的元素，而不碰其他的。由于我们在所有处理器上做同样的操作，我们最终得到一个分布式向量（即在每个进程上只存储与该进程拥有的自由度相对应的向量项），该向量的所有受限节点都被固定。
 * 

 * 
 * 我们通过返回收敛所需的迭代次数来结束这个函数，以允许一些输出。
 * 

 * 
 * 
 * @code
 *     solution = localized_solution; 
 * 
 *     return solver_control.last_step(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrefine_grid"></a> 
 * <h4>ElasticProblem::refine_grid</h4>
 * 

 * 
 * 使用某种细化指标，可以对网格进行细化。这个问题与分布悬挂节点约束基本相同：为了计算误差指标（即使我们只是对当前进程拥有的单元上的指标感兴趣），我们需要访问解向量的更多元素，而不仅仅是当前处理器存储的那些元素。为了实现这一点，我们基本上做了我们在 <code>solve()</code> 中已经做过的事情，即获取<i>complete</i>解向量的副本到每个进程中，并使用它来计算。如上所述，这本身就很昂贵，尤其是没有必要，因为我们刚刚在 <code>solve()</code> 中创建并销毁了这样一个向量，但效率并不是这个程序的重点，所以让我们选择一种设计，即每个函数都尽可能地独立。
 * 

 * 
 * 一旦我们有了这样一个包含<i>all</i>解向量元素的 "本地化 "向量，我们就可以计算属于当前过程的单元的指标。事实上，我们当然可以计算<i>all</i>细化指标，因为我们的Triangulation和DoFHandler对象存储了所有单元的信息，而且我们有一个完整的解向量副本。但是为了展示如何进行%并行操作，让我们演示一下，如果只计算<i>some</i>错误指标，然后与其他进程交换剩余的指标，会如何操作。(最终，每个进程都需要一套完整的细化指标，因为每个进程都需要细化他们的网格，并且需要以与其他进程完全相同的方式细化它。)
 * 

 * 
 * 所以，为了做到这一切，我们需要。
 * 

 * 
 * - 首先，获得分布式求解向量的本地拷贝。
 * 

 * 
 * - 第二，创建一个向量来存储细化指标。
 * 

 * 
 * - 第三，让KellyErrorEstimator计算属于当前子域/过程的所有单元的细化指标。调用的最后一个参数表明我们对哪个子域感兴趣。在它之前的三个参数是其他各种默认参数，通常不需要（也不说明数值，而是使用默认值），但我们必须在这里明确说明，因为我们要修改下面一个参数的值（即表示子域的参数）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::refine_grid() 
 *   { 
 *     const Vector<double> localized_solution(solution); 
 * 
 *     Vector<float> local_error_per_cell(triangulation.n_active_cells()); 
 *     KellyErrorEstimator<dim>::estimate(dof_handler, 
 *                                        QGauss<dim - 1>(fe.degree + 1), 
 *                                        {}, 
 *                                        localized_solution, 
 *                                        local_error_per_cell, 
 *                                        ComponentMask(), 
 *                                        nullptr, 
 *                                        MultithreadInfo::n_threads(), 
 *                                        this_mpi_process); 
 * 
 * @endcode
 * 
 * 现在所有进程都计算了自己单元格的错误指标，并将其存储在 <code>local_error_per_cell</code> 向量的相应元素中。这个向量中不属于本进程的单元格的元素为零。然而，由于所有进程都有整个三角形的副本，并需要保持这些副本的同步，他们需要三角形的所有单元的细化指标值。因此，我们需要分配我们的结果。我们通过创建一个分布式向量来做到这一点，每个进程都有自己的份额，并设置它所计算的元素。因此，当你把这个向量看作是一个存在于所有进程中的向量时，那么这个向量的每个元素都被设置过一次。然后，我们可以将这个并行向量分配给每个进程上的一个本地非并行向量，使<i>all</i>错误指示器在每个进程上都可用。
 * 因此，
 * 在第一步，我们需要设置一个并行向量。为了简单起见，每个进程都将拥有一个元素块，其数量与该进程拥有的单元格一样多，因此第一个元素块存储在进程0，下一个元素块存储在进程1，以此类推。然而，需要注意的是，这些元素不一定是我们要写入的元素。这是单元格排列顺序的结果，也就是说，向量中的元素对应单元格的顺序并不是根据这些单元格所属的子域来排序的。换句话说，如果在这个过程中，我们计算某个子域的单元的指标，我们可能会把结果写到分布式向量的或多或少的随机元素中；特别是，它们不一定位于我们在这个过程中拥有的向量块中。它们随后将不得不被复制到另一个进程的内存空间中，当我们调用 <code>compress()</code> 函数时，PETSc为我们做了这项操作。这种低效率可以通过更多的代码来避免，但我们不这样做，因为它不是程序总运行时间的一个主要因素。
 * 

 * 
 * 所以我们是这样做的：计算有多少个单元属于这个过程，建立一个有这么多元素的分布式向量存储在本地，将我们在本地计算的元素复制过去，最后将结果压缩。事实上，我们实际上只复制了非零的元素，所以我们可能会错过一些我们计算为零的元素，但这不会有什么影响，因为无论如何，向量的原始值是零。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_local_cells = 
 *       GridTools::count_cells_with_subdomain_association(triangulation, 
 *                                                         this_mpi_process); 
 *     PETScWrappers::MPI::Vector distributed_all_errors( 
 *       mpi_communicator, triangulation.n_active_cells(), n_local_cells); 
 * 
 *     for (unsigned int i = 0; i < local_error_per_cell.size(); ++i) 
 *       if (local_error_per_cell(i) != 0) 
 *         distributed_all_errors(i) = local_error_per_cell(i); 
 *     distributed_all_errors.compress(VectorOperation::insert); 
 * 
 * @endcode
 * 
 * 所以现在我们有了这个分布式向量，它包含了所有单元的细化指标。为了使用它，我们需要获得一个本地副本，然后用它来标记要细化或粗化的单元，并实际进行细化和粗化。重要的是要认识到，<i>every</i>过程对它自己的三角形副本做了这个工作，并且以完全相同的方式进行。
 * 

 * 
 * 
 * @code
 *     const Vector<float> localized_all_errors(distributed_all_errors); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     localized_all_errors, 
 *                                                     0.3, 
 *                                                     0.03); 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemoutput_results"></a> 
 * <h4>ElasticProblem::output_results</h4>
 * 

 * 
 * 最后一个有意义的函数是创建图形输出的函数。它的工作方式与 step-8 中的相同，但有两个小的区别。在讨论这些之前，让我们说明这个函数的一般工作原理：我们打算让所有的数据都在一个进程中产生，然后写入一个文件中。正如本程序的许多其他部分已经讨论过的那样，这不是一个可以扩展的东西。之前，我们认为我们会在三角计算、DoFHandlers和解决方案向量的副本方面遇到麻烦，每个进程都必须存储所有的数据，而且会出现一个点，即每个进程根本没有足够的内存来存储这么多数据。在这里，情况是不同的：不仅是内存，而且运行时间也是一个问题。如果一个进程负责处理<i>all</i>的数据，而其他所有的进程什么都不做，那么这一个函数最终会在程序的整个运行时间中占主导地位。 特别是，这个函数花费的时间将与问题的整体大小（以单元数或自由度数计算）成正比，与我们扔给它的进程数量无关。
 * 

 * 
 * 这种情况需要避免，我们将在 step-18 和 step-40 中展示如何解决这个问题。对于目前的问题，解决方案是让每个进程只为自己的本地单元产生输出数据，并将它们写入单独的文件，每个进程一个文件。这就是 step-18 的操作方式。另外，我们可以简单地把所有的东西放在一组独立的文件中，让可视化软件读取所有的文件（可能也使用多个处理器），并从所有的文件中创建一个单一的可视化；这就是 step-40 、 step-32 以及后来开发的所有其他并行程序的路径。
 * 

 * 
 * 更具体地说，对于当前的函数，所有的进程都调用这个函数，但不是所有的进程都需要做与生成输出相关的工作。事实上，它们不应该这样做，因为我们会试图一次多次地写到同一个文件。所以我们只让第一个进程做这件事，而其他所有的进程在这段时间内闲置（或者为下一次迭代开始工作，或者干脆把它们的CPU让给碰巧在同一时间运行的其他作业）。第二件事是，我们不仅要输出解决方案的向量，还要输出一个向量，表明每个单元属于哪个子域。这将使一些分区域的图片变得很好。
 * 

 * 
 * 为了实现这一点，过程0需要一个完整的本地向量中的解决方案组件。就像前面的函数一样，有效的方法是重新使用在 <code>solve()</code> 函数中已经创建的向量，但是为了使事情更加自洽，我们在这里简单地从分布式解决方案向量中重新创建一个向量。
 * 

 * 
 * 需要认识到的一个重要问题是，我们在所有的进程中都做了这个定位操作，而不是只有那个实际需要数据的进程。然而，这一点是无法避免的，在本教程程序中，我们对向量使用的MPI简化通信模型。MPI没有办法查询另一个进程的数据，双方必须在同一时间启动通信。因此，即使大多数进程不需要本地化的解决方案，我们也必须把将分布式转换为本地化向量的语句放在那里，以便所有进程都执行它。
 * 

 * 
 * （这项工作的一部分实际上可以避免。我们所做的是将所有进程的本地部分发送给所有其他进程。我们真正需要做的是在所有进程上发起一个操作，每个进程只需将其本地的数据块发送给进程0，因为只有这个进程才真正需要它，也就是说，我们需要类似于收集操作的东西。PETSc可以做到这一点，但是为了简单起见，我们在这里并不试图利用这一点。我们没有这样做，因为我们所做的事情在整个计划中并不昂贵：它是所有进程之间的一个矢量通信，这必须与我们在求解线性系统、为预处理程序设置块状ILU以及其他操作时必须进行的通信数量相比较。)
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     const Vector<double> localized_solution(solution); 
 * 
 * @endcode
 * 
 * 这样做后，零进程继续设置输出文件，如  step-8  ，并将（本地化的）解决方案矢量附加到输出对象上。
 * 

 * 
 * 
 * @code
 *     if (this_mpi_process == 0) 
 *       { 
 *         std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
 * 
 *         DataOut<dim> data_out; 
 *         data_out.attach_dof_handler(dof_handler); 
 * 
 *         std::vector<std::string> solution_names; 
 *         switch (dim) 
 *           { 
 *             case 1: 
 *               solution_names.emplace_back("displacement"); 
 *               break; 
 *             case 2: 
 *               solution_names.emplace_back("x_displacement"); 
 *               solution_names.emplace_back("y_displacement"); 
 *               break; 
 *             case 3: 
 *               solution_names.emplace_back("x_displacement"); 
 *               solution_names.emplace_back("y_displacement"); 
 *               solution_names.emplace_back("z_displacement"); 
 *               break; 
 *             default: 
 *               Assert(false, ExcInternalError()); 
 *           } 
 * 
 *         data_out.add_data_vector(localized_solution, solution_names); 
 * 
 * @endcode
 * 
 * 我们在这里做的唯一其他事情是，我们也为每个单元格输出一个值，表明它属于哪个子域（即MPI进程）。这需要一些转换工作，因为库提供给我们的数据不是输出类所期望的数据，但这并不困难。首先，设置一个整数向量，每个单元格一个，然后由每个单元格的子域id填充。
 * 

 * 
 * 这个向量的元素在第二步中被转换为浮点向量，这个向量被添加到DataOut对象中，然后它去创建VTK格式的输出。
 * 

 * 
 * 
 * @code
 *         std::vector<unsigned int> partition_int(triangulation.n_active_cells()); 
 *         GridTools::get_subdomain_association(triangulation, partition_int); 
 * 
 *         const Vector<double> partitioning(partition_int.begin(), 
 *                                           partition_int.end()); 
 * 
 *         data_out.add_data_vector(partitioning, "partitioning"); 
 * 
 *         data_out.build_patches(); 
 *         data_out.write_vtk(output); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrun"></a> 
 * <h4>ElasticProblem::run</h4>
 * 

 * 
 * 最后，这里是驱动程序的功能。它与 step-8 几乎完全没有变化，只是我们替换了 <code>std::cout</code> by the <code>pcout</code> 流。除此以外，唯一的表面变化是我们输出了每个进程有多少个自由度，以及线性求解器花了多少次收敛。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 10; ++cycle) 
 *       { 
 *         pcout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation, -1, 1); 
 *             triangulation.refine_global(3); 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 *         pcout << "   Number of active cells:       " 
 *               << triangulation.n_active_cells() << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         pcout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (by partition:"; 
 *         for (unsigned int p = 0; p < n_mpi_processes; ++p) 
 *           pcout << (p == 0 ? ' ' : '+') 
 *                 << (DoFTools::count_dofs_with_subdomain_association(dof_handler, 
 *                                                                     p)); 
 *         pcout << ")" << std::endl; 
 * 
 *         assemble_system(); 
 *         const unsigned int n_iterations = solve(); 
 * 
 *         pcout << "   Solver converged in " << n_iterations << " iterations." 
 *               << std::endl; 
 * 
 *         output_results(cycle); 
 *       } 
 *   } 
 * } // namespace Step17 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * <code>main()</code> 的工作方式与其他示例程序中的大多数主函数相同，即它将工作委托给管理对象的 <code>run</code> 函数，并且只将所有内容包装成一些代码来捕获异常。
 * 

 * 
 * 
 * @code
 * int main(int argc, char **argv) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step17; 
 * 
 * @endcode
 * 
 * 这里是唯一真正的区别。MPI和PETSc都要求我们在程序开始时初始化这些库，并在结束时解除初始化。MPI_InitFinalize类处理了所有这些。后面的参数`1`意味着我们确实想让每个MPI进程以单线程运行，这是PETSc并行线性代数的前提条件。
 * 

 * 
 * 
 * @code
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 * 
 *       ElasticProblem<2> elastic_problem; 
 *       elastic_problem.run(); 
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
examples/step-17/doc/results.dox



<a name="Results"></a><h1>Results</h1>



如果上述程序在单处理器机器上编译和运行，它产生的结果应该与我们已经通过步骤8得到的结果非常相似。然而，如果我们在多核机器或计算机集群上运行它，就会变得更加有趣。运行MPI程序的最基本方法是使用一个命令行，如

@code
  mpirun -np 32 ./step-17
@endcode

以32个处理器运行step-17可执行文件。

如果你在一个集群上工作，那么中间通常有一个步骤，你需要设置一个作业脚本，并将该脚本提交给调度器。只要调度器能够为你的工作分配32个未使用的处理器，它就会执行这个脚本。如何编写这样的作业脚本因集群而异，你应该找到你的集群的文档来看看如何做。在我的系统上，我必须使用 <code>qsub</code> 这个命令，加上一大堆的选项来并行运行一个作业）。)

无论是直接还是通过调度器，如果你在8个处理器上运行这个程序，你应该得到如下输出。

@code
Cycle 0:
   Number of active cells:       64
   Number of degrees of freedom: 162 (by partition: 22+22+20+20+18+16+20+24)
   Solver converged in 23 iterations.
Cycle 1:
   Number of active cells:       124
   Number of degrees of freedom: 302 (by partition: 38+42+36+34+44+44+36+28)
   Solver converged in 35 iterations.
Cycle 2:
   Number of active cells:       238
   Number of degrees of freedom: 570 (by partition: 68+80+66+74+58+68+78+78)
   Solver converged in 46 iterations.
Cycle 3:
   Number of active cells:       454
   Number of degrees of freedom: 1046 (by partition: 120+134+124+130+154+138+122+124)
   Solver converged in 55 iterations.
Cycle 4:
   Number of active cells:       868
   Number of degrees of freedom: 1926 (by partition: 232+276+214+248+230+224+234+268)
   Solver converged in 77 iterations.
Cycle 5:
   Number of active cells:       1654
   Number of degrees of freedom: 3550 (by partition: 418+466+432+470+442+474+424+424)
   Solver converged in 93 iterations.
Cycle 6:
   Number of active cells:       3136
   Number of degrees of freedom: 6702 (by partition: 838+796+828+892+866+798+878+806)
   Solver converged in 127 iterations.
Cycle 7:
   Number of active cells:       5962
   Number of degrees of freedom: 12446 (by partition: 1586+1484+1652+1552+1556+1576+1560+1480)
   Solver converged in 158 iterations.
Cycle 8:
   Number of active cells:       11320
   Number of degrees of freedom: 23586 (by partition: 2988+2924+2890+2868+2864+3042+2932+3078)
   Solver converged in 225 iterations.
Cycle 9:
   Number of active cells:       21424
   Number of degrees of freedom: 43986 (by partition: 5470+5376+5642+5450+5630+5470+5416+5532)
   Solver converged in 282 iterations.
Cycle 10:
   Number of active cells:       40696
   Number of degrees of freedom: 83754 (by partition: 10660+10606+10364+10258+10354+10322+10586+10604)
   Solver converged in 392 iterations.
Cycle 11:
   Number of active cells:       76978
   Number of degrees of freedom: 156490 (by partition: 19516+20148+19390+19390+19336+19450+19730+19530)
   Solver converged in 509 iterations.
Cycle 12:
   Number of active cells:       146206
   Number of degrees of freedom: 297994 (by partition: 37462+37780+37000+37060+37232+37328+36860+37272)
   Solver converged in 705 iterations.
Cycle 13:
   Number of active cells:       276184
   Number of degrees of freedom: 558766 (by partition: 69206+69404+69882+71266+70348+69616+69796+69248)
   Solver converged in 945 iterations.
Cycle 14:
   Number of active cells:       523000
   Number of degrees of freedom: 1060258 (by partition: 132928+132296+131626+132172+132170+133588+132252+133226)
   Solver converged in 1282 iterations.
Cycle 15:
   Number of active cells:       987394
   Number of degrees of freedom: 1994226 (by partition: 253276+249068+247430+248402+248496+251380+248272+247902)
   Solver converged in 1760 iterations.
Cycle 16:
   Number of active cells:       1867477
   Number of degrees of freedom: 3771884 (by partition: 468452+474204+470818+470884+469960+
471186+470686+475694)
   Solver converged in 2251 iterations.
@endcode

(这次运行比examples/目录中的代码多用了几个细化周期。该运行还使用了2004年的METIS版本，产生了不同的分区；因此，你今天得到的数字略有不同）。)

可以看出，我们可以很容易地达到近400万个未知数。事实上，这段代码在8个进程中的运行时间不到7分钟，直到（包括）第14个周期，14分钟包括倒数第二步。(这些数字与该代码最初编写的时间有关，即2004年。)虽然我失去了最后一步的时间信息，但你会明白的。所有这些都是在通过运行 <code>make release</code> 启用发布模式之后，并且由于上述程序注释中所述的原因，关闭了图形输出的生成。(  @dealiiVideoLectureSeeAlso{18})  我做的最大的2D计算大约有710万个未知数，在32个进程上完成。花了大约40分钟。毫不奇怪，一个人能够走多远的限制因素是他有多少内存，因为每个进程都必须持有整个网格和DoFHandler对象，尽管矩阵和向量被分割开来。对于7.1M的计算，每个未知数的内存消耗约为600字节，这并不坏，但我们必须考虑到这是针对每个未知数的，无论我们是否在本地存储矩阵和向量条目。




下面是在程序的第12个周期中产生的一些输出，即大约有30万个未知数。

 <table align="center" style="width:80%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-ux.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-uy.png" alt="" width="100%"></td>
  </tr>
</table> 

正如人们所希望的那样，这里显示的X位移（左）和Y位移（右）与我们在第8步中已经看到的密切相关。正如第22步所示，我们也可以制作一个位移场的矢量图，而不是把它绘制成两个独立的标量场。不过，可能更有趣的是，在这一步看一下网格和分区。

 <table align="center" width="80%">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-grid.png" alt="" width="100%"></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.12-partition.png" alt="" width="100%"></td>
  </tr>
</table> 

同样，网格（左边）显示了与之前看到的相同的细化模式。右图显示了8个过程中的域的划分，每个过程用不同的颜色表示。图片显示，在网格单元较小的地方，子域较小，考虑到分区算法试图平衡每个子域的单元数，这一事实是需要预期的；这种平衡在上图所示的输出中也很容易识别，每个子域的度数大致相同。




值得注意的是，如果我们用不同的进程数来运行同一个程序，我们可能会得到稍微不同的输出：不同的网格，不同的未知数和迭代收敛的次数。其原因是，虽然矩阵和右手边是相同的，与使用的进程数无关，但预处理程序不是：它对每个处理器的 <em>  矩阵块分别执行ILU(0)  </em>  。因此，随着进程数的增加，它作为预处理程序的有效性会降低，这使得迭代次数增加。由于不同的预处理程序会导致计算出的解有细微的变化，这将导致细化时标记的网格单元略有不同，在后续步骤中的差异也更大。不过，解决方案看起来总是非常相似的。




最后，这里是3D模拟的一些结果。你可以通过改变以下内容来重复这些结果

@code
        ElasticProblem<2> elastic_problem;
@endcode

至

@code
        ElasticProblem<3> elastic_problem;
@endcode

在主函数中。如果你再并行运行该程序，你会得到与此类似的东西（这是针对一个有16个进程的工作）。

@code
Cycle 0:
   Number of active cells:       512
   Number of degrees of freedom: 2187 (by partition: 114+156+150+114+114+210+105+102+120+120+96+123+141+183+156+183)
   Solver converged in 27 iterations.
Cycle 1:
   Number of active cells:       1604
   Number of degrees of freedom: 6549 (by partition: 393+291+342+354+414+417+570+366+444+288+543+525+345+387+489+381)
   Solver converged in 42 iterations.
Cycle 2:
   Number of active cells:       4992
   Number of degrees of freedom: 19167 (by partition: 1428+1266+1095+1005+1455+1257+1410+1041+1320+1380+1080+1050+963+1005+1188+1224)
   Solver converged in 65 iterations.
Cycle 3:
   Number of active cells:       15485
   Number of degrees of freedom: 56760 (by partition: 3099+3714+3384+3147+4332+3858+3615+3117+3027+3888+3942+3276+4149+3519+3030+3663)
   Solver converged in 96 iterations.
Cycle 4:
   Number of active cells:       48014
   Number of degrees of freedom: 168762 (by partition: 11043+10752+9846+10752+9918+10584+10545+11433+12393+11289+10488+9885+10056+9771+11031+8976)
   Solver converged in 132 iterations.
Cycle 5:
   Number of active cells:       148828
   Number of degrees of freedom: 492303 (by partition: 31359+30588+34638+32244+30984+28902+33297+31569+29778+29694+28482+28032+32283+30702+31491+28260)
   Solver converged in 179 iterations.
Cycle 6:
   Number of active cells:       461392
   Number of degrees of freedom: 1497951 (by partition: 103587+100827+97611+93726+93429+88074+95892+88296+96882+93000+87864+90915+92232+86931+98091+90594)
   Solver converged in 261 iterations.
@endcode






最后一步，达到150万个未知数，在8台双处理器机器（2003年可用的那种）上进行16个进程，需要大约55分钟。这个工作产生的图形输出相当大（第5周期已经打印了大约82MB的数据），所以我们要显示第4周期的输出。

 <table width="80%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.4-3d-partition.png" width="100%" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-17.4-3d-ux.png" alt="" width="100%"></td>
  </tr>
</table> 




左图显示的是将立方体划分为16个过程，而右图显示的是沿两个切面通过立方体的X位移。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


该程序在每个处理器上都保留一份三角形和DoFHandler对象的完整副本。它还创建了解决方案矢量的完整副本，并且只在一个处理器上创建输出。就并行化而言，所有这些显然是瓶颈。

在内部，在deal.II中，将分层和非结构化三角计算中使用的数据结构并行化是一个难点，我们又花了几年时间才实现了这一点。step-40教程程序和 @ref distributed 文档模块谈到了如何做这些步骤，以及从应用的角度来看需要什么。当前程序的一个明显的扩展是使用这个功能将计算完全分布到比这里使用的更多的处理器。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-17.cc"
*/
