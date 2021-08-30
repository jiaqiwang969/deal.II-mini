/**
@page step_16 The step-16 tutorial program
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
        <li><a href="#TheScratchandCopyobjects">The Scratch and Copy objects</a>
        <li><a href="#ThecodeLaplaceProblemcodeclasstemplate">The <code>LaplaceProblem</code> class template</a>
        <li><a href="#ThecodeLaplaceProblemcodeclassimplementation">The <code>LaplaceProblem</code> class implementation</a>
      <ul>
        <li><a href="#LaplaceProblemsetup_system">LaplaceProblem::setup_system</a>
        <li><a href="#LaplaceProblemcell_worker">LaplaceProblem::cell_worker</a>
        <li><a href="#LaplaceProblemassemble_system">LaplaceProblem::assemble_system</a>
        <li><a href="#LaplaceProblemassemble_multigrid">LaplaceProblem::assemble_multigrid</a>
        <li><a href="#LaplaceProblemsolve">LaplaceProblem::solve</a>
        <li><a href="#Postprocessing">Postprocessing</a>
        <li><a href="#LaplaceProblemrun">LaplaceProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-16/doc/intro.dox

 <br> 

<i> Note: A variant called step-16b of this tutorial exists, that uses
MeshWorker and LocalIntegrators instead of assembling matrices manually as it
is done in this tutorial.
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这个例子展示了deal.II中多级函数的基本用法。它几乎解决了与步骤6中使用的相同的问题，但展示了使用多网格作为预处理程序时必须提供的东西。特别是，这要求我们定义一个层次结构，提供从一个层次到下一个层次以及返回的转移算子，并在每个层次上提供拉普拉斯算子的表示。

为了使微分方程系统和块状预处理程序具有足够的灵活性，在启动多级方法之前必须创建一些不同的对象，尽管大部分需要做的事情都是由deal.II本身提供。这些对象是

  - 网格之间的对象处理转移；我们使用MGTransferPrebuilt类来处理这个问题，它几乎完成了库内的所有工作。

  - 解算器的最粗层次；在这里，我们使用MGCoarseGridHouseholder。

  - 所有其他级别的平滑器，在我们的例子中，这将是使用SOR作为基本方法的 mg::SmootherRelaxation 类。

  - 和 mg::Matrix, 一个具有特殊水平乘法的类，也就是说，我们基本上每个网格水平存储一个矩阵并允许与之相乘。

这些对象中的大多数只需要在实际求解线性系统的函数中使用。在这里，这些对象被组合到一个多网格类型的对象中，其中包含V型循环的实现，它又被预设条件器PreconditionMG使用，准备插入LAC库的线性求解器中。

这里实现的自适应细化网格的多网格方法遵循 @ref mg_paper "多网格论文 "中的大纲，该论文在deal.II中描述了底层实现，也介绍了很多术语。首先，我们必须区分层次网格，即与粗网格有相同细化距离的单元，以及由层次中的活动单元组成的叶子网格（在较早的工作中，我们将其称为全局网格，但这个术语被过度使用）。最重要的是，叶子网格与最细层次上的层次网格不完全相同。下面的图片显示了我们认为的 "层次网"。

<p align="center">  @image html "multigrid.png" ""   </p>  。

这个网格中的精细层次只包括定义在精细单元上的自由度，但不延伸到领域中未被精细化的那部分。虽然这保证了整体的努力增长为 ${\cal O}(N)$ 的最佳多网格复杂度所必需的，但它导致了在定义平滑的地方和对定义在各个层次上的算子提出什么边界条件时的问题，如果层次边界不是外部边界。这些问题将在上面引用的文章中详细讨论。

<a name="Thetestcase"></a><h3>The testcase</h3>


我们在这里解决的问题与第6步类似，主要有两个不同点：第一，多网格预处理程序，显然。我们还改变了系数的不连续性，使局部装配器看起来不会比必要的更复杂。


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
 * 同样，前几个include文件已经知道了，所以我们不会对它们进行评论。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 现在，这些是多级方法所需的包括。第一个声明了如何处理多网格方法每个层次上的Dirichlet边界条件。对于自由度的实际描述，我们不需要任何新的包含文件，因为DoFHandler已经实现了所有必要的方法。我们只需要将自由度分配给更多的层次。
 * 

 * 
 * 其余的包含文件涉及到作为线性算子（求解器或预处理器）的多重网格的力学问题。
 * 

 * 
 * 
 * @code
 * #include <deal.II/multigrid/mg_constrained_dofs.h> 
 * #include <deal.II/multigrid/multigrid.h> 
 * #include <deal.II/multigrid/mg_transfer.h> 
 * #include <deal.II/multigrid/mg_tools.h> 
 * #include <deal.II/multigrid/mg_coarse.h> 
 * #include <deal.II/multigrid/mg_smoother.h> 
 * #include <deal.II/multigrid/mg_matrix.h> 
 * 
 * @endcode
 * 
 * 我们将使用 MeshWorker::mesh_loop 来对单元格进行循环，所以在这里包括它。
 * 

 * 
 * 
 * @code
 * #include <deal.II/meshworker/mesh_loop.h> 
 * 
 * @endcode
 * 
 * 这就是C++。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * 
 * using namespace dealii; 
 * 
 * namespace Step16 
 * { 
 * @endcode
 * 
 * 
 * <a name="TheScratchandCopyobjects"></a> 
 * <h3>The Scratch and Copy objects</h3>
 * 

 * 
 * 我们使用 MeshWorker::mesh_loop() 来组装我们的矩阵。为此，我们需要一个ScratchData对象来存储每个单元的临时数据（这只是FEValues对象）和一个CopyData对象，它将包含每个单元装配的输出。关于scratch和copy对象的用法的更多细节，请参见WorkStream命名空间。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct ScratchData 
 *   { 
 *     ScratchData(const Mapping<dim> &      mapping, 
 *                 const FiniteElement<dim> &fe, 
 *                 const unsigned int        quadrature_degree, 
 *                 const UpdateFlags         update_flags) 
 *       : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
 *     {} 
 * 
 *     ScratchData(const ScratchData<dim> &scratch_data) 
 *       : fe_values(scratch_data.fe_values.get_mapping(), 
 *                   scratch_data.fe_values.get_fe(), 
 *                   scratch_data.fe_values.get_quadrature(), 
 *                   scratch_data.fe_values.get_update_flags()) 
 *     {} 
 * 
 *     FEValues<dim> fe_values; 
 *   }; 
 * 
 *   struct CopyData 
 *   { 
 *     unsigned int                         level; 
 *     FullMatrix<double>                   cell_matrix; 
 *     Vector<double>                       cell_rhs; 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 * 
 *     template <class Iterator> 
 *     void reinit(const Iterator &cell, unsigned int dofs_per_cell) 
 *     { 
 *       cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
 *       cell_rhs.reinit(dofs_per_cell); 
 * 
 *       local_dof_indices.resize(dofs_per_cell); 
 *       cell->get_active_or_mg_dof_indices(local_dof_indices); 
 *       level = cell->level(); 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceProblem</code> class template</h3>
 * 

 * 
 * 这个主类与  step-6  中的同一类相似。就成员函数而言，唯一增加的是。
 * 

 * 
 * --  <code>assemble_multigrid</code> 的函数，该函数组装了对应于中间层离散运算符的矩阵。
 * 

 * 
 * -  <code>cell_worker</code> 函数，它将我们的PDE集合在一个单元上。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class LaplaceProblem 
 *   { 
 *   public: 
 *     LaplaceProblem(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     template <class Iterator> 
 *     void cell_worker(const Iterator &  cell, 
 *                      ScratchData<dim> &scratch_data, 
 *                      CopyData &        copy_data); 
 * 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void assemble_multigrid(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 * 
 *     const unsigned int degree; 
 * 
 * @endcode
 * 
 * 以下成员是多网格方法的基本数据结构。前四个表示稀疏模式和多级层次结构中各个层次的矩阵，非常类似于上面的全局网格的对象。
 * 

 * 
 * 然后，我们有两个新的矩阵，只需要在自适应网格上进行局部平滑的多网格方法。它们在细化区域的内部和细化边缘之间传递数据，在 @ref mg_paper "多网格论文 "中详细介绍过。
 * 

 * 
 * 最后一个对象存储了每个层次上的边界指数信息和位于两个不同细化层次之间的细化边缘上的指数信息。因此，它的作用与AffineConstraints相似，但在每个层次上。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<SparsityPattern> mg_sparsity_patterns; 
 *     MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns; 
 * 
 *     MGLevelObject<SparseMatrix<double>> mg_matrices; 
 *     MGLevelObject<SparseMatrix<double>> mg_interface_matrices; 
 *     MGConstrainedDoFs                   mg_constrained_dofs; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceProblemcodeclassimplementation"></a> 
 * <h3>The <code>LaplaceProblem</code> class implementation</h3>
 * 

 * 
 * 关于三角形的构造函数只有一个简短的评论：按照惯例，deal.II中所有自适应精化的三角形在单元格之间的面的变化不会超过一个级别。然而，对于我们的多网格算法，我们需要一个更严格的保证，即网格在连接两个单元的顶点上的变化也不超过细化级别。换句话说，我们必须防止出现以下情况。
 * 

 * 
 * @image html limit_level_difference_at_vertices.png ""  
 * 

 * 
 * 这可以通过向三角化类的构造函数传递 Triangulation::limit_level_difference_at_vertices 标志来实现。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree) 
 *     : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
 *     , fe(degree) 
 *     , dof_handler(triangulation) 
 *     , degree(degree) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsetup_system"></a> 
 * <h4>LaplaceProblem::setup_system</h4>
 * 

 * 
 * 除了只是在DoFHandler中分配自由度之外，我们在每一层都做同样的事情。然后，我们按照之前的程序，在叶子网格上设置系统。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     dof_handler.distribute_mg_dofs(); 
 * 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (by level: "; 
 *     for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
 *       std::cout << dof_handler.n_dofs(level) 
 *                 << (level == triangulation.n_levels() - 1 ? ")" : ", "); 
 *     std::cout << std::endl; 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 * 
 *     std::set<types::boundary_id> dirichlet_boundary_ids = {0}; 
 *     Functions::ZeroFunction<dim> homogeneous_dirichlet_bc; 
 *     const std::map<types::boundary_id, const Function<dim> *> 
 *       dirichlet_boundary_functions = { 
 *         {types::boundary_id(0), &homogeneous_dirichlet_bc}}; 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              dirichlet_boundary_functions, 
 *                                              constraints); 
 *     constraints.close(); 
 * 
 *     { 
 *       DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
 *       sparsity_pattern.copy_from(dsp); 
 *     } 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 多网格约束必须被初始化。他们需要知道在哪里规定了Dirichlet边界条件。
 * 

 * 
 * 
 * @code
 *     mg_constrained_dofs.clear(); 
 *     mg_constrained_dofs.initialize(dof_handler); 
 *     mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
 *                                                        dirichlet_boundary_ids); 
 * 
 * @endcode
 * 
 * 现在是关于多网格数据结构的事情。首先，我们调整多级对象的大小，以容纳每一级的矩阵和稀疏模式。粗略的级别是零（现在是强制性的，但在未来的修订中可能会改变）。注意，这些函数在这里采取的是一个完整的、包容的范围（而不是一个起始索引和大小），所以最细的级别是 <code>n_levels-1</code>  。我们首先要调整容纳SparseMatrix类的容器的大小，因为它们必须在调整大小时释放它们的SparsityPattern才能被销毁。
 * 

 * 
 * 
 * @code
 *     const unsigned int n_levels = triangulation.n_levels(); 
 * 
 *     mg_interface_matrices.resize(0, n_levels - 1); 
 *     mg_matrices.resize(0, n_levels - 1); 
 *     mg_sparsity_patterns.resize(0, n_levels - 1); 
 *     mg_interface_sparsity_patterns.resize(0, n_levels - 1); 
 * 
 * @endcode
 * 
 * 现在，我们必须在每个级别上提供一个矩阵。为此，我们首先使用 MGTools::make_sparsity_pattern 函数在每个层次上生成一个初步的压缩稀疏模式（关于这个主题的更多信息，请参见 @ref Sparsity 模块），然后将其复制到我们真正想要的那一个。下一步是用拟合的稀疏度模式初始化接口矩阵。
 * 

 * 
 * 值得指出的是，界面矩阵只包含位于较粗和较细的网格之间的自由度的条目。因此，它们甚至比我们的多网格层次结构中的各个层次的矩阵还要稀疏。因此，我们使用一个专门为此目的而建立的函数来生成它。
 * 

 * 
 * 
 * @code
 *     for (unsigned int level = 0; level < n_levels; ++level) 
 *       { 
 *         { 
 *           DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
 *                                      dof_handler.n_dofs(level)); 
 *           MGTools::make_sparsity_pattern(dof_handler, dsp, level); 
 * 
 *           mg_sparsity_patterns[level].copy_from(dsp); 
 *           mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
 *         } 
 *         { 
 *           DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
 *                                      dof_handler.n_dofs(level)); 
 *           MGTools::make_interface_sparsity_pattern(dof_handler, 
 *                                                    mg_constrained_dofs, 
 *                                                    dsp, 
 *                                                    level); 
 *           mg_interface_sparsity_patterns[level].copy_from(dsp); 
 *           mg_interface_matrices[level].reinit( 
 *             mg_interface_sparsity_patterns[level]); 
 *         } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemcell_worker"></a> 
 * <h4>LaplaceProblem::cell_worker</h4>
 * 

 * 
 * cell_worker函数用于在给定的单元上组装矩阵和右手边。这个函数用于活动单元生成system_matrix，并在每个层次上建立层次矩阵。
 * 

 * 
 * 注意，当从assemble_multigrid()调用时，我们也会组装一个右手边，尽管它没有被使用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   template <class Iterator> 
 *   void LaplaceProblem<dim>::cell_worker(const Iterator &  cell, 
 *                                         ScratchData<dim> &scratch_data, 
 *                                         CopyData &        copy_data) 
 *   { 
 *     FEValues<dim> &fe_values = scratch_data.fe_values; 
 *     fe_values.reinit(cell); 
 * 
 *     const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = fe_values.get_quadrature().size(); 
 * 
 *     copy_data.reinit(cell, dofs_per_cell); 
 * 
 *     const std::vector<double> &JxW = fe_values.get_JxW_values(); 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         const double coefficient = 
 *           (fe_values.get_quadrature_points()[q][0] < 0.0) ? 1.0 : 0.1; 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           { 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               { 
 *                 copy_data.cell_matrix(i, j) += 
 *                   coefficient * 
 *                   (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) * 
 *                   JxW[q]; 
 *               } 
 *             copy_data.cell_rhs(i) += 1.0 * fe_values.shape_value(i, q) * JxW[q]; 
 *           } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_system"></a> 
 * <h4>LaplaceProblem::assemble_system</h4>
 * 

 * 
 * 下面的函数将线性系统集合在网格的活动单元上。为此，我们向Mesh_loop()函数传递两个lambda函数。cell_worker函数重定向到同名的类成员函数，而copyer是这个函数特有的，它使用约束条件将本地矩阵和向量复制到相应的全局矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::assemble_system() 
 *   { 
 *     MappingQ1<dim> mapping; 
 * 
 *     auto cell_worker = 
 *       [&](const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *           ScratchData<dim> &                                    scratch_data, 
 *           CopyData &                                            copy_data) { 
 *         this->cell_worker(cell, scratch_data, copy_data); 
 *       }; 
 * 
 *     auto copier = [&](const CopyData &cd) { 
 *       this->constraints.distribute_local_to_global(cd.cell_matrix, 
 *                                                    cd.cell_rhs, 
 *                                                    cd.local_dof_indices, 
 *                                                    system_matrix, 
 *                                                    system_rhs); 
 *     }; 
 * 
 *     const unsigned int n_gauss_points = degree + 1; 
 * 
 *     ScratchData<dim> scratch_data(mapping, 
 *                                   fe, 
 *                                   n_gauss_points, 
 *                                   update_values | update_gradients | 
 *                                     update_JxW_values | 
 *                                     update_quadrature_points); 
 * 
 *     MeshWorker::mesh_loop(dof_handler.begin_active(), 
 *                           dof_handler.end(), 
 *                           cell_worker, 
 *                           copier, 
 *                           scratch_data, 
 *                           CopyData(), 
 *                           MeshWorker::assemble_own_cells); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemassemble_multigrid"></a> 
 * <h4>LaplaceProblem::assemble_multigrid</h4>
 * 

 * 
 * 下一个函数是建立矩阵，定义每一层网格上的多网格方法。集成的核心与上面的相同，但是下面的循环将遍历所有已存在的单元，而不仅仅是活动的单元，并且必须将结果输入正确的层矩阵。幸运的是，MeshWorker对我们隐藏了大部分的内容，因此这个函数和之前的函数的区别只在于汇编器的设置和循环中的不同迭代器。
 * 

 * 
 * 我们为每个层次生成一个AffineConstraints对象，其中包含边界和界面道夫作为约束条目。然后，相应的对象被用来生成层次矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::assemble_multigrid() 
 *   { 
 *     MappingQ1<dim>     mapping; 
 *     const unsigned int n_levels = triangulation.n_levels(); 
 * 
 *     std::vector<AffineConstraints<double>> boundary_constraints(n_levels); 
 *     for (unsigned int level = 0; level < n_levels; ++level) 
 *       { 
 *         IndexSet dofset; 
 *         DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
 *                                                       level, 
 *                                                       dofset); 
 *         boundary_constraints[level].reinit(dofset); 
 *         boundary_constraints[level].add_lines( 
 *           mg_constrained_dofs.get_refinement_edge_indices(level)); 
 *         boundary_constraints[level].add_lines( 
 *           mg_constrained_dofs.get_boundary_indices(level)); 
 *         boundary_constraints[level].close(); 
 *       } 
 * 
 *     auto cell_worker = 
 *       [&](const typename DoFHandler<dim>::level_cell_iterator &cell, 
 *           ScratchData<dim> &                                   scratch_data, 
 *           CopyData &                                           copy_data) { 
 *         this->cell_worker(cell, scratch_data, copy_data); 
 *       }; 
 * 
 *     auto copier = [&](const CopyData &cd) { 
 *       boundary_constraints[cd.level].distribute_local_to_global( 
 *         cd.cell_matrix, cd.local_dof_indices, mg_matrices[cd.level]); 
 * 
 *       const unsigned int dofs_per_cell = cd.local_dof_indices.size(); 
 * 
 * @endcode
 * 
 * 接口条目在填充mg_matrices[cd.level]时被上面的boundary_constraints对象所忽略。相反，我们手动将这些条目复制到当前级别的界面矩阵中。
 * 

 * 
 * 
 * @code
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *         for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *           if (mg_constrained_dofs.is_interface_matrix_entry( 
 *                 cd.level, cd.local_dof_indices[i], cd.local_dof_indices[j])) 
 *             { 
 *               mg_interface_matrices[cd.level].add(cd.local_dof_indices[i], 
 *                                                   cd.local_dof_indices[j], 
 *                                                   cd.cell_matrix(i, j)); 
 *             } 
 *     }; 
 * 
 *     const unsigned int n_gauss_points = degree + 1; 
 * 
 *     ScratchData<dim> scratch_data(mapping, 
 *                                   fe, 
 *                                   n_gauss_points, 
 *                                   update_values | update_gradients | 
 *                                     update_JxW_values | 
 *                                     update_quadrature_points); 
 * 
 *     MeshWorker::mesh_loop(dof_handler.begin_mg(), 
 *                           dof_handler.end_mg(), 
 *                           cell_worker, 
 *                           copier, 
 *                           scratch_data, 
 *                           CopyData(), 
 *                           MeshWorker::assemble_own_cells); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemsolve"></a> 
 * <h4>LaplaceProblem::solve</h4>
 * 

 * 
 * 这是另外一个在支持多栅求解器（或者说，事实上，我们使用多栅方法的前提条件）方面有明显不同的函数。
 * 

 * 
 * 让我们从建立多层次方法的两个组成部分开始：层次间的转移运算器和最粗层次上的求解器。在有限元方法中，转移算子来自所涉及的有限元函数空间，通常可以用独立于所考虑问题的通用方式计算。在这种情况下，我们可以使用MGTransferPrebuilt类，给定最终线性系统的约束和MGConstrainedDoFs对象，该对象知道每个层次的边界条件和不同细化层次之间接口的自由度，可以从具有层次自由度的DoFHandler对象中建立这些转移操作的矩阵。
 * 

 * 
 * 下面几行的第二部分是关于粗略网格求解器的。由于我们的粗网格确实非常粗，我们决定采用直接求解器（最粗层次矩阵的Householder分解），即使其实现不是特别复杂。如果我们的粗网格比这里的5个单元多得多，那么这里显然需要更合适的东西。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::solve() 
 *   { 
 *     MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs); 
 *     mg_transfer.build(dof_handler); 
 * 
 *     FullMatrix<double> coarse_matrix; 
 *     coarse_matrix.copy_from(mg_matrices[0]); 
 *     MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver; 
 *     coarse_grid_solver.initialize(coarse_matrix); 
 * 
 * @endcode
 * 
 * 多级求解器或预处理器的下一个组成部分是，我们需要在每一级上设置平滑器。这方面常见的选择是使用松弛方法的应用（如SOR、Jacobi或Richardson方法）或求解器方法的少量迭代（如CG或GMRES）。 mg::SmootherRelaxation 和MGSmootherPrecondition类为这两种平滑器提供支持。这里，我们选择应用单一的SOR迭代。为此，我们定义一个适当的别名，然后设置一个平滑器对象。
 * 

 * 
 * 最后一步是用我们的水平矩阵初始化平滑器对象，并设置一些平滑参数。 <code>initialize()</code> 函数可以有选择地接受额外的参数，这些参数将被传递给每一级的平滑器对象。在目前SOR平滑器的情况下，这可能包括一个松弛参数。然而，我们在这里将这些参数保留为默认值。对 <code>set_steps()</code> 的调用表明我们将在每个级别上使用两个前平滑步骤和两个后平滑步骤；为了在不同级别上使用可变数量的平滑器步骤，可以在对 <code>mg_smoother</code> 对象的构造函数调用中设置更多选项。
 * 

 * 
 * 最后一步的结果是我们使用SOR方法作为平滑器的事实
 * 

 * 
 * --这不是对称的
 * 

 * 
 * 但我们在下面使用共轭梯度迭代（需要对称的预处理），我们需要让多级预处理确保我们得到一个对称的算子，即使是非对称的平滑器。
 * 

 * 
 * 
 * @code
 *     using Smoother = PreconditionSOR<SparseMatrix<double>>; 
 *     mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother; 
 *     mg_smoother.initialize(mg_matrices); 
 *     mg_smoother.set_steps(2); 
 *     mg_smoother.set_symmetric(true); 
 * 
 * @endcode
 * 
 * 下一个准备步骤是，我们必须将我们的水平和接口矩阵包裹在一个具有所需乘法函数的对象中。我们将为从粗到细的接口对象创建两个对象，反之亦然；多网格算法将在后面的操作中使用转置运算器，允许我们用已经建立的矩阵初始化该运算器的上下版本。
 * 

 * 
 * 
 * @code
 *     mg::Matrix<Vector<double>> mg_matrix(mg_matrices); 
 *     mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices); 
 *     mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices); 
 * 
 * @endcode
 * 
 * 现在，我们准备设置V型循环算子和多级预处理程序。
 * 

 * 
 * 
 * @code
 *     Multigrid<Vector<double>> mg( 
 *       mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother); 
 *     mg.set_edge_matrices(mg_interface_down, mg_interface_up); 
 * 
 *     PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
 *       preconditioner(dof_handler, mg, mg_transfer); 
 * 
 * @endcode
 * 
 * 有了这一切，我们终于可以用通常的方法来解决这个线性系统了。
 * 

 * 
 * 
 * @code
 *     SolverControl            solver_control(1000, 1e-12); 
 *     SolverCG<Vector<double>> solver(solver_control); 
 * 
 *     solution = 0; 
 * 
 *     solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 *     std::cout << "   Number of CG iterations: " << solver_control.last_step() 
 *               << "\n" 
 *               << std::endl; 
 *     constraints.distribute(solution); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Postprocessing"></a> 
 * <h4>Postprocessing</h4>
 * 

 * 
 * 以下两个函数在计算出解决方案后对其进行后处理。特别是，第一个函数在每个周期开始时细化网格，第二个函数在每个周期结束时输出结果。这些函数与  step-6  中的函数几乎没有变化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::refine_grid() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(degree + 2), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       estimated_error_per_cell); 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.03); 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 *   template <int dim> 
 *   void LaplaceProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "solution"); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="LaplaceProblemrun"></a> 
 * <h4>LaplaceProblem::run</h4>
 * 

 * 
 * 和上面的几个函数一样，这几乎是对  step-6  中相应函数的复制。唯一的区别是对 <code>assemble_multigrid</code> 的调用，它负责形成我们在多网格方法中需要的每一层的矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void LaplaceProblem<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 8; ++cycle) 
 *       { 
 *         std::cout << "Cycle " << cycle << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_ball(triangulation); 
 *             triangulation.refine_global(2); 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 *         std::cout << "   Number of active cells:       " 
 *                   << triangulation.n_active_cells() << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         assemble_system(); 
 *         assemble_multigrid(); 
 * 
 *         solve(); 
 *         output_results(cycle); 
 *       } 
 *   } 
 * } // namespace Step16 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 这又是与 step-6 中相同的函数。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step16; 
 * 
 *       LaplaceProblem<2> laplace_problem(1); 
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
 * 
 * 
 * @endcode
examples/step-16/doc/results.dox



<a name="Results"></a><h1>Results</h1>


在最细的网格上，解决方案看起来像这样。

<p align="center">  <img src="https://www.dealii.org/images/steps/developer/step-16.solution.png" alt="">   </p>  。

更重要的是，我们想看看多网格方法是否真的改善了求解器的性能。因此，这里是文本输出。

<p>周期0 活动单元数：80 自由度数：89（按级别：8，25，89） CG迭代数。8

周期1 活动单元数：158 自由度数：183（按级别：8，25，89，138） CG迭代次数。9

周期2 活动单元数：302 自由度数：352（按级别：8、25、89、223、160） CG迭代次数。10

周期3 活动单元数：578 自由度数：649（按级别：8、25、89、231、494、66） CG迭代次数。10

第四周期 活跃单元数：1100 自由度数：1218（按级别：8、25、89、274、764、417、126） CG迭代次数。10

周期5 活动单元数：2096 自由度数：2317（按级别：8、25、89、304、779、1214、817） CG迭代次数。11

第6周期 活跃单元数：3986 自由度数：4366（按级别：8，25，89，337，836，2270，897，1617） CG迭代数。10

周期7 活动单元数：7574 自由度数：8350（按级别：8、25、89、337、1086、2835、2268、1789、3217） CG迭代次数。11 </pre>

这几乎是完美的多重网格性能：线性残差在10个迭代步骤中被减少了12个数量级，而且结果几乎与网格大小无关。这显然部分是由于所解决的问题的简单性质，但它显示了多梯度方法的力量。




<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>



我们鼓励你生成solve()调用的时间并与步骤6进行比较。你会看到多网格方法在粗网格上有相当大的开销，但由于其最佳的复杂性，它在细网格上总是胜过其他方法。

仔细检查这个程序的性能就会发现，它主要是由矩阵-向量操作主导的。step-37展示了一种通过使用无矩阵方法来避免这种情况的方法。

另一个途径是使用代数多网格方法。这里使用的几何多栅方法在实现上有时会有点尴尬，因为它需要所有这些额外的数据结构，如果程序要在%parallel的机器上通过MPI耦合运行，就会变得更加困难，例如。在这种情况下，如果能使用一个黑盒预处理程序，就会更简单，该程序使用某种多栅层次结构以获得良好的性能，但可以自己计算出水平矩阵和类似的东西。代数多栅方法正是这样做的，我们将在步骤31中使用它们来解决斯托克斯问题，在步骤32和步骤40中使用它们来解决平行变化。也就是说，在步骤50中可以找到这个例子程序的MPI并行版本。

最后，人们可能要考虑如何将几何多网格用于其他类型的问题，特别是 @ref vector_valued "矢量值问题"。这是step-56的主题，我们将这里的技术用于斯托克斯方程。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-16.cc"
*/
