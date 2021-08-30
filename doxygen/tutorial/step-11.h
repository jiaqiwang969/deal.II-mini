/**
@page step_11 The step-11 tutorial program
This tutorial depends on step-10.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-11/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


我们要考虑的问题是只带诺伊曼边界条件的拉普拉斯问题的解决方案。

@f{eqnarray*}


  -\Delta u &=& f \qquad \mathrm{in}\ \Omega,
  \\
  \partial_n u &=& g \qquad \mathrm{on}\ \partial\Omega.


@f}

众所周知，如果这个问题要有一个解决方案，那么力需要满足兼容性条件

@f[
  \int_\Omega f\; dx + \int_{\partial\Omega} g\; ds = 0.


@f]

我们将考虑这样的特殊情况： $\Omega$ 是围绕原点的半径为1的圆，而 $f=-2$  ， $g=1$  。这种选择满足了兼容性条件。

兼容性条件允许上述方程的解，但它仍然保留了一个模糊性：因为只有解的导数出现在方程中，解只确定到一个常数。出于这个原因，我们必须为数字解提出另一个条件，以固定这个常数。

对于这一点，有多种可能性。<ol>  <li>  将离散化的一个节点固定为零或任何其他固定值。   这相当于一个附加条件  $u_h(x_0)=0$  。虽然这是常见的做法，但不一定是个好主意，因为我们知道拉普拉斯方程的解只在  $H^1$  中，这不允许定义点值，因为它不是连续函数的子集。因此，即使固定一个节点对离散函数来说是允许的，但对连续函数来说是不允许的，在数值解的这一点上，人们常常可以看到由此产生的错误尖峰。

 <li>  将域上的均值固定为零或任何其他值。这在连续水平上是允许的，因为 $H^1(\Omega)\subset L^1(\Omega)$ 由Sobolev不等式决定，因此在离散水平上也是允许的，因为我们那里只考虑 $H^1$ 的子集。

 <li>  将域的边界上的均值固定为零或任何其他值。这在连续水平上也是允许的，因为 $H^{1/2}(\partial\Omega)\subset L^1(\partial\Omega)$  ，同样由Sobolev的不等式。   </ol>  我们将选择最后一种可能性，因为我们想用它来演示另一种技术。

虽然这描述了要解决的问题，但我们仍然要弄清楚如何实现它。基本上，除了额外的均值约束，我们已经多次解决了这个问题，使用的是迪里希特边界值，我们只需要放弃对迪里希特边界节点的处理。高阶映射的使用也是相当琐碎的，我们会在使用它的各个地方进行解释；在几乎所有可以想象的情况下，你只会把描述映射的对象视为一个黑盒子，你不需要担心，因为它们的唯一用途似乎是被传递到库的深处，在那里函数知道如何处理它们（即在 <code>FEValues</code> 类及其后代）。

这个程序中的棘手之处在于对均值约束的使用。幸运的是，库中有一个知道如何处理这种约束的类，我们已经经常使用它了，没有提到它的通用性。请注意，如果我们假设边界节点沿边界的间隔是相等的，那么均值约束

@f[
  \int_{\partial \Omega} u(x) \; ds = 0


@f]

可写为

@f[
  \sum_{i\in\partial\Omega_h} u_i = 0,


@f]

其中总和应贯穿位于计算域边界上的所有自由度指数。让我们用 $i_0$ 表示边界上数字最小的指数（或任何其他方便选择的指数），那么这个约束也可以用以下方式表示

@f[
  u_{i_0} = \sum_{i\in\partial\Omega_h\backslash i_0} -u_i.


@f]

幸运的是，这正是AffineConstraints类所设计的约束形式。请注意，我们在之前的几个例子中使用了这个类来表示悬空节点的约束，它也有这种形式：在这里，中间的顶点应具有相邻顶点的平均值。一般来说，AffineConstraints类被设计用来处理以下形式的仿生约束

@f[
  CU = b


@f]

其中 $C$ 表示一个矩阵， $b$ 表示一个向量， $U$ 是节点值的向量。在这种情况下，由于 $C$ 代表一个同质约束， $b$ 是零向量。

在这个例子中，沿边界的平均值允许这样的表示， $C$ 是一个只有一行的矩阵（即只有一个约束条件）。在实现中，我们将创建一个AffineConstraints对象，添加一个参考第一个边界节点 $i_0$ 的约束（即给矩阵添加另一行），并插入所有其他节点贡献的权重，在这个例子中刚好是 $-1$  。

稍后，我们将使用这个对象来消除线性方程组中的第一个边界节点，将其还原为一个没有常数偏移值的解。实施过程中的一个问题是，这个节点的明确消除会导致矩阵中出现一些额外的元素，我们事先不知道这些元素的位置，也不知道矩阵的每一行中会有多少额外的条目。我们将展示我们如何使用一个中间对象来解决这个问题。

但现在开始实施解决这个问题的方案......


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 像往常一样，程序以一个相当长的包含文件列表开始，你现在可能已经习惯了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/table_handler.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_q.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * @endcode
 * 
 * 只有这一条是新的：它声明了一个动态稀疏模式（DynamicSparsityPattern）类，我们将在下面进一步使用和解释。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * @endcode
 * 
 * 我们将使用C++标准库中的 std::find 算法，所以我们必须包括以下文件来声明它。
 * 

 * 
 * 
 * @code
 * #include <algorithm> 
 * #include <iostream> 
 * #include <iomanip> 
 * #include <cmath> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step11 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 然后我们声明一个表示拉普拉斯问题解决方案的类。由于这个例子程序是基于 step-5 ，这个类看起来相当相同，唯一的结构区别是函数 <code>assemble_system</code> now calls <code>solve</code> 本身，因此被称为 <code>assemble_and_solve</code> ，而且输出函数被删除，因为解函数非常无聊，不值得查看。
 * 

 * 
 * 其他唯一值得注意的变化是，构造函数取一个值，代表以后要使用的映射的多项式程度，而且它还有一个成员变量，正好代表这个映射。一般来说，这个变量在实际应用中会出现在声明或使用有限元的相同地方。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class LaplaceProblem 
 *   { 
 *   public: 
 *     LaplaceProblem(const unsigned int mapping_degree); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_and_solve(); 
 *     void solve(); 
 *     void write_high_order_mesh(const unsigned cycle); 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 *     MappingQ<dim>      mapping; 
 * 
 *     SparsityPattern           sparsity_pattern; 
 *     SparseMatrix<double>      system_matrix; 
 *     AffineConstraints<double> mean_value_constraints; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 * 
 *     TableHandler output_table; 
 *   }; 
 * 
 * @endcode
 * 
 * 构建这样一个对象，通过初始化变量。这里，我们使用线性有限元（ <code>fe</code> 变量的参数表示多项式的度数），以及给定阶数的映射。将我们要做的事情打印到屏幕上。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   LaplaceProblem<dim>::LaplaceProblem(const unsigned int mapping_degree) 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *     , mapping(mapping_degree) 
 *   { 
 *     std::cout << "Using mapping with degree " << mapping_degree << ":" 
 *               << std::endl 
 *               << "============================" << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 第一个任务是为这个问题设置变量。这包括生成一个有效的 <code>DoFHandler</code> 对象，以及矩阵的稀疏模式，和代表边界上自由度平均值为零的约束条件的对象。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::setup_system() 
 *   { 
 * 
 * @endcode
 * 
 * 第一个任务很简单：生成一个自由度的枚举，并将解和右手向量初始化为正确的大小。
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 * @endcode
 * 
 * 下一个任务是构建代表约束的对象，即边界上自由度的平均值应该是零。为此，我们首先需要一个实际在边界上的节点的列表。 <code>DoFTools</code> 命名空间有一个函数可以返回一个IndexSet对象，该对象包含所有在边界上的自由度的指数。
 * 

 * 
 * 一旦我们有了这个索引集，我们想知道哪个是对应于边界上的自由度的第一个索引。我们需要这个，因为我们想通过边界上所有其他自由度的值来约束边界上的一个节点。使用IndexSet类很容易得到这个 "第一个 "自由度的索引。
 * 

 * 
 * 
 * @code
 *     const IndexSet boundary_dofs = DoFTools::extract_boundary_dofs(dof_handler); 
 * 
 *     const types::global_dof_index first_boundary_dof = 
 *       boundary_dofs.nth_index_in_set(0); 
 * 
 * @endcode
 * 
 * 然后生成一个只有这一个约束的约束对象。首先清除所有以前的内容（这些内容可能来自以前在更粗的网格上的计算），然后添加这一行，将 <code>first_boundary_dof</code> 约束到其他边界DoF的总和，每一个权重为-1。最后，关闭约束对象，也就是说，对它做一些内部记录，以便更快地处理后面的内容。
 * 

 * 
 * 
 * @code
 *     mean_value_constraints.clear(); 
 *     mean_value_constraints.add_line(first_boundary_dof); 
 *     for (types::global_dof_index i : boundary_dofs) 
 *       if (i != first_boundary_dof) 
 *         mean_value_constraints.add_entry(first_boundary_dof, i, -1); 
 *     mean_value_constraints.close(); 
 * 
 * @endcode
 * 
 * 下一个任务是生成一个稀疏模式。这的确是一个棘手的任务。通常情况下，我们只需调用 <code>DoFTools::make_sparsity_pattern</code> 并使用悬挂节点约束来浓缩结果。我们在这里没有悬挂节点约束（因为我们在这个例子中只进行全局细化），但是我们在边界上有这个全局约束。在这种情况下，这带来了一个严重的问题： <code>SparsityPattern</code> 类希望我们事先说明每行的最大条目数，可以是所有行的，也可以是每行单独的。在库中有一些函数可以告诉你这个数字，如果你只有悬空的节点约束的话（即 DoFHandler::max_couplings_between_dofs), ，但这对现在的情况来说是怎样的？困难的出现是因为消除约束的自由度需要在矩阵中增加一些条目，而这些条目的位置并不那么容易确定。因此，如果我们在这里给出每行的最大条目数，我们就会有一个问题。
 * 

 * 
 * 由于这可能非常困难，以至于无法给出合理的答案，只能分配合理的内存量，所以有一个DynamicSparsityPattern类，它可以帮助我们解决这个问题。它不要求我们事先知道行可以有多少个条目，而是允许任何长度。因此，在你对行的长度没有很好的估计的情况下，它明显更灵活，但是代价是建立这样一个模式也比建立一个你事先有信息的模式要昂贵得多。尽管如此，由于我们在这里没有其他选择，我们将建立这样一个对象，用矩阵的尺寸初始化它，并调用另一个函数 <code>DoFTools::make_sparsity_pattern</code> 来获得由于微分算子引起的稀疏模式，然后用约束对象浓缩它，在稀疏模式中增加那些消除约束所需的位置。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     mean_value_constraints.condense(dsp); 
 * 
 * @endcode
 * 
 * 最后，一旦我们有了完整的模式，我们就可以从中初始化一个 <code>SparsityPattern</code> 类型的对象，并反过来用它初始化矩阵。请注意，这实际上是必要的，因为与 <code>SparsityPattern</code> 类相比，DynamicSparsityPattern的效率非常低，因为它必须使用更灵活的数据结构，所以我们不可能将稀疏矩阵类建立在它的基础上，而是需要一个 <code>SparsityPattern</code> 类型的对象，我们通过复制中间对象产生这个对象。
 * 

 * 
 * 作为进一步的附带说明，你会注意到我们在这里没有明确的  <code>compress</code>  稀疏模式。当然，这是由于 <code>copy_from</code> 函数从一开始就生成了一个压缩对象，你不能再向其添加新的条目。因此， <code>compress</code> 的调用是隐含在 <code>copy_from</code> 的调用中的。
 * 

 * 
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp); 
 *     system_matrix.reinit(sparsity_pattern); 
 *   } 
 * 
 * @endcode
 * 
 * 下一个函数接着组装线性方程组，对其进行求解，并对解进行评估。这样就有了三个动作，我们将把它们放到八个真实的语句中（不包括变量的声明，以及临时向量的处理）。因此，这个函数是为非常懒惰的人准备的。尽管如此，所调用的函数是相当强大的，通过它们，这个函数使用了整个库的大量内容。但让我们来看看每一个步骤。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::assemble_and_solve() 
 *   { 
 * 
 * @endcode
 * 
 * 首先，我们要把矩阵和右手边的内容组合起来。在之前的所有例子中，我们已经研究了如何手动完成这一工作的各种方法。然而，由于拉普拉斯矩阵和简单的右手边在应用中出现的频率很高，库中提供的函数实际上是为你做这件事的，也就是说，它们在所有单元格上进行循环，设置局部的矩阵和向量，并将它们放在一起，得到最终结果。
 * 

 * 
 * 以下是两个最常用的函数：创建拉普拉斯矩阵和创建来自体或边界力的右侧向量。它们需要映射对象、代表自由度和使用中的有限元的 <code>DoFHandler</code> 对象、要使用的正交公式以及输出对象。创建右手向量的函数还必须接受一个描述（连续）右手向量函数的函数对象。
 * 

 * 
 * 让我们来看看矩阵和体力的集成方式。
 * 

 * 
 * 
 * @code
 *     const unsigned int gauss_degree = 
 *       std::max(static_cast<unsigned int>( 
 *                  std::ceil(1. * (mapping.get_degree() + 1) / 2)), 
 *                2U); 
 *     MatrixTools::create_laplace_matrix(mapping, 
 *                                        dof_handler, 
 *                                        QGauss<dim>(gauss_degree), 
 *                                        system_matrix); 
 *     VectorTools::create_right_hand_side(mapping, 
 *                                         dof_handler, 
 *                                         QGauss<dim>(gauss_degree), 
 *                                         Functions::ConstantFunction<dim>(-2), 
 *                                         system_rhs); 
 * 
 * @endcode
 * 
 * 这很简单，对吗？
 * 

 * 
 * 不过，有两点需要注意。首先，这些函数在很多情况下都会用到。也许你想为一个矢量值有限元创建一个拉普拉斯或质量矩阵；或者你想使用默认的Q1映射；或者你想用拉普拉斯算子的一个系数来装配矩阵。由于这个原因，在 <code>MatrixCreator</code> 和 <code>MatrixTools</code> 命名空间中有相当多的这些函数的变种。每当你需要这些函数的一个与上面调用的略有不同的版本时，当然值得看一下文档，并检查一些东西是否适合你的需要。
 * 

 * 
 * 第二点是关于我们使用的正交公式：我们想对双线性形状函数进行积分，所以我们知道我们至少要使用二阶高斯正交公式。另一方面，我们希望正交规则至少有边界近似的阶数。因为有 $r$ 点的高斯规则的阶数是 $2r -1$  ，而使用 $p$ 度的多项式的边界近似的阶数是 $p+1$ ，我们知道 $2r \geq p$  。由于r必须是一个整数，并且（如上所述）必须至少是 $2$ ，这就弥补了上述公式计算 <code>gauss_degree</code> 。
 * 

 * 
 * 由于对右侧向量的体力贡献的生成是如此简单，我们对边界力也要重新做一遍：分配一个合适大小的向量并调用合适的函数。边界函数有常量值，所以我们可以从库中快速生成一个对象，我们使用与上面相同的正交公式，但这次的维度较低，因为我们现在是在面上而不是在单元上积分。
 * 

 * 
 * 
 * @code
 *     Vector<double> tmp(system_rhs.size()); 
 *     VectorTools::create_boundary_right_hand_side( 
 *       mapping, 
 *       dof_handler, 
 *       QGauss<dim - 1>(gauss_degree), 
 *       Functions::ConstantFunction<dim>(1), 
 *       tmp); 
 * 
 * @endcode
 * 
 * 然后将边界的贡献与域内部的贡献相加。
 * 

 * 
 * 
 * @code
 *     system_rhs += tmp; 
 * 
 * @endcode
 * 
 * 在组装右手边时，我们必须使用两个不同的矢量对象，然后将它们加在一起。我们不得不这样做的原因是， <code>VectorTools::create_right_hand_side</code> 和 <code>VectorTools::create_boundary_right_hand_side</code> 函数首先清除输出向量，而不是将它们的结果与之前的内容相加。这可以合理地称为库在起步阶段的设计缺陷，但不幸的是，事情现在已经是这样了，很难改变这种无声地破坏现有代码的事情，所以我们不得不接受。
 * 

 * 
 * 现在，线性系统已经建立起来了，所以我们可以从矩阵和右手向量中消除我们约束到边界上其他DoF的一个自由度的均值约束，并解决这个系统。之后，再次分配约束，在这种情况下，这意味着将被约束的自由度设置为适当的值
 * 

 * 
 * 
 * @code
 *     mean_value_constraints.condense(system_matrix); 
 *     mean_value_constraints.condense(system_rhs); 
 * 
 *     solve(); 
 *     mean_value_constraints.distribute(solution); 
 * 
 * @endcode
 * 
 * 最后，评估我们得到的解决方案。正如在介绍中所说，我们对解决方案的H1半正态感兴趣。在这里，我们在库中也有一个函数可以做到这一点，尽管是以一种稍微不明显的方式： <code>VectorTools::integrate_difference</code> 函数整合了一个有限元函数和一个连续函数之间的差值的规范。因此，如果我们想要一个有限元场的规范，我们只需将连续函数设为零。请注意，这个函数，就像库中的许多其他函数一样，至少有两个版本，一个是以映射为参数的（我们在这里使用），另一个是我们在以前的例子中使用的隐含的 <code>MappingQ1</code>  。 还要注意的是，我们采用的是高一级的正交公式，以避免出现超融合效应，即在某些点上的解特别接近精确解（我们不知道这里是否会出现这种情况，但有已知的案例，我们只是想确认一下）。
 * 

 * 
 * 
 * @code
 *     Vector<float> norm_per_cell(triangulation.n_active_cells()); 
 *     VectorTools::integrate_difference(mapping, 
 *                                       dof_handler, 
 *                                       solution, 
 *                                       Functions::ZeroFunction<dim>(), 
 *                                       norm_per_cell, 
 *                                       QGauss<dim>(gauss_degree + 1), 
 *                                       VectorTools::H1_seminorm); 
 * 
 * @endcode
 * 
 * 然后，刚刚调用的函数将其结果作为一个值的向量返回，每个值表示一个单元格上的法线。为了得到全局法线，我们要做以下工作。
 * 

 * 
 * 
 * @code
 *     const double norm = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         norm_per_cell, 
 *                                         VectorTools::H1_seminorm); 
 * 
 * @endcode
 * 
 * 最后一项任务--生成输出。
 * 

 * 
 * 
 * @code
 *     output_table.add_value("cells", triangulation.n_active_cells()); 
 *     output_table.add_value("|u|_1", norm); 
 *     output_table.add_value("error", 
 *                            std::fabs(norm - std::sqrt(3.14159265358 / 2))); 
 *   } 
 * 
 * @endcode
 * 
 * 下面这个解线性方程组的函数是从 step-5 中复制过来的，在那里有详细的解释。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::solve() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-12); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 *   } 
 * 
 * @endcode
 * 
 * 接下来，我们把解决方案以及材料ID写到一个VTU文件中。这与其他许多教程程序中的做法相似。这个教程程序中提出的新内容是，我们要确保写到文件中用于可视化的数据实际上是deal.II内部使用的数据的忠实代表。这是因为大多数可视化数据格式只用顶点坐标表示单元，但没有办法表示deal.II中使用高阶映射时的曲线边界--换句话说，你在可视化工具中看到的东西实际上不是你正在计算的东西。顺带一提，在使用高阶形状函数时也是如此。大多数可视化工具只呈现双线性/三线性的表示。这在 DataOut::build_patches().) 中有详细的讨论。
 * 

 * 
 * 所以我们需要确保高阶表示被写入文件中。我们需要考虑两个特别的话题。首先，我们通过 DataOutBase::VtkFlags 告诉DataOut对象，我们打算将元素的细分解释为高阶拉格朗日多项式，而不是双线性斑块的集合。最近的可视化程序，如ParaView 5.5版或更新版，然后可以呈现高阶解决方案（更多细节见<a
 * href="https:github.com/dealii/dealii/wiki/Notes-on-visualizing-high-order-output">wiki
 * page</a>）。其次，我们需要确保映射被传递给 DataOut::build_patches() 方法。最后，DataOut类默认只打印<i>boundary</i>单元的曲面，所以我们需要确保通过映射将内部单元也打印成曲面。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::write_high_order_mesh(const unsigned cycle) 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     DataOutBase::VtkFlags flags; 
 *     flags.write_higher_order_cells = true; 
 *     data_out.set_flags(flags); 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 * 
 *     data_out.build_patches(mapping, 
 *                            mapping.get_degree(), 
 *                            DataOut<dim>::curved_inner_cells); 
 * 
 *     std::ofstream file("solution-c=" + std::to_string(cycle) + 
 *                        ".p=" + std::to_string(mapping.get_degree()) + ".vtu"); 
 * 
 *     data_out.write_vtu(file); 
 *   } 
 * 
 * @endcode
 * 
 * 最后是控制要执行的不同步骤的主要函数。它的内容相当简单，生成一个圆的三角形，给它关联一个边界，然后在随后的更细的网格上做几个循环。请注意，我们将网格细化放到了循环头中；这对测试程序来说可能是件好事，但对实际应用来说，你应该考虑到这意味着网格是在循环最后一次执行后被细化的，因为增量子句（三部分循环头的最后一部分）是在比较部分（第二部分）之前执行的，如果网格已经相当细化了，这可能是相当昂贵的。在这种情况下，你应该安排代码，使网格在最后一次循环运行后不再被进一步细化（或者你应该在每次运行的开始就这样做，除了第一次）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::run() 
 *   { 
 *     GridGenerator::hyper_ball(triangulation); 
 * 
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle) 
 *       { 
 *         setup_system(); 
 *         assemble_and_solve(); 
 *         write_high_order_mesh(cycle); 
 * 
 *         triangulation.refine_global(); 
 *       } 
 * 
 * @endcode
 * 
 * 在所有的数据生成之后，将结果的表格写到屏幕上。
 * 

 * 
 * 
 * @code
 *     output_table.set_precision("|u|_1", 6); 
 *     output_table.set_precision("error", 6); 
 *     output_table.write_text(std::cout); 
 *     std::cout << std::endl; 
 *   } 
 * } // namespace Step11 
 * 
 * @endcode
 * 
 * 最后是主函数。它的结构与前面几个例子中使用的结构相同，所以可能不需要更多解释。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       std::cout.precision(5); 
 * 
 * @endcode
 * 
 * 这是主循环，用线性到立方的映射做计算。注意，由于我们只需要一次 <code>LaplaceProblem@<2@></code> 类型的对象，我们甚至不给它命名，而是创建一个未命名的这样的对象，并调用它的 <code>run</code> 函数，随后它又立即被销毁。
 * 

 * 
 * 
 * @code
 *       for (unsigned int mapping_degree = 1; mapping_degree <= 3; 
 *            ++mapping_degree) 
 *         Step11::LaplaceProblem<2>(mapping_degree).run(); 
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
examples/step-11/doc/results.dox



<a name="Results"></a><h1>Results</h1>


这是该程序的输出结果。

@code
Using mapping with degree 1:
============================
cells  |u|_1    error
    5 0.680402 0.572912
   20 1.088141 0.165173
   80 1.210142 0.043172
  320 1.242375 0.010939
 1280 1.250569 0.002745
 5120 1.252627 0.000687


Using mapping with degree 2:
============================
cells  |u|_1    error
    5 1.177062 0.076252
   20 1.228978 0.024336
   80 1.245175 0.008139
  320 1.250881 0.002433
 1280 1.252646 0.000668
 5120 1.253139 0.000175


Using mapping with degree 3:
============================
cells  |u|_1    error
    5 1.193493 0.059821
   20 1.229825 0.023489
   80 1.245221 0.008094
  320 1.250884 0.002430
 1280 1.252646 0.000668
 5120 1.253139 0.000175
@endcode

正如我们所期望的，每个不同的映射的收敛顺序显然是与网格大小成二次方的。不过 <em> 有趣的是，双线性映射（即1度）的误差比高阶映射的误差大三倍以上；因此在这种情况下，使用高阶映射显然是有利的，不是因为它提高了收敛顺序，而只是为了减少收敛顺序前的常数。另一方面，除了在非常粗的网格上，使用立方体映射只能进一步改善结果，幅度微乎其微。

我们还可以通过使用例如ParaView来可视化底层网格。下面的图片显示了不同映射度的初始网格。

 <img src="https://www.dealii.org/images/steps/developer/step-11.cycle_0.png" alt=""> 

显然，当我们从线性映射到二次映射时，这种影响是最明显的。这也反映在上表中给出的误差值中。从二次方到三次方的效果没有那么明显，但由于对圆形边界的描述更加准确，所以还是很明显的。

接下来，让我们看看三次全局细化后的网格

 <img src="https://www.dealii.org/images/steps/developer/step-11.cycle_3.png" alt=""> 

在这里，差异就不那么明显了，特别是对于高阶映射。事实上，在这个细化水平上，表格中报告的误差值在二度和三度的映射之间基本上是相同的。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-11.cc"
*/
