CCTest_file/step-16.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2003 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Authors: Guido Kanschat, University of Heidelberg, 2003 
 *          Baerbel Janssen, University of Heidelberg, 2010 
 *          Wolfgang Bangerth, Texas A&M University, 2010 
 *          Timo Heister, Clemson University, 2018 
 */ 


// @sect3{Include files}  

// 同样，前几个include文件已经知道了，所以我们不会对它们进行评论。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 现在，这些是多级方法所需的包括。第一个声明了如何处理多网格方法每个层次上的Dirichlet边界条件。对于自由度的实际描述，我们不需要任何新的包含文件，因为DoFHandler已经实现了所有必要的方法。我们只需要将自由度分配给更多的层次。

// 其余的包含文件涉及到作为线性算子（求解器或预处理器）的多重网格的力学问题。

#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

// 我们将使用 MeshWorker::mesh_loop 来对单元格进行循环，所以在这里包括它。

#include <deal.II/meshworker/mesh_loop.h> 

// 这就是C++。

#include <iostream> 
#include <fstream> 

using namespace dealii; 

namespace Step16 
{ 
// @sect3{The Scratch and Copy objects}  

// 我们使用 MeshWorker::mesh_loop() 来组装我们的矩阵。为此，我们需要一个ScratchData对象来存储每个单元的临时数据（这只是FEValues对象）和一个CopyData对象，它将包含每个单元装配的输出。关于scratch和copy对象的用法的更多细节，请参见WorkStream命名空间。

  template <int dim> 
  struct ScratchData 
  { 
    ScratchData(const Mapping<dim> &      mapping, 
                const FiniteElement<dim> &fe, 
                const unsigned int        quadrature_degree, 
                const UpdateFlags         update_flags) 
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
    {} 

    ScratchData(const ScratchData<dim> &scratch_data) 
      : fe_values(scratch_data.fe_values.get_mapping(), 
                  scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  scratch_data.fe_values.get_update_flags()) 
    {} 

    FEValues<dim> fe_values; 
  }; 

  struct CopyData 
  { 
    unsigned int                         level; 
    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_rhs; 
    std::vector<types::global_dof_index> local_dof_indices; 

    template <class Iterator> 
    void reinit(const Iterator &cell, unsigned int dofs_per_cell) 
    { 
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
      cell_rhs.reinit(dofs_per_cell); 

      local_dof_indices.resize(dofs_per_cell); 
      cell->get_active_or_mg_dof_indices(local_dof_indices); 
      level = cell->level(); 
    } 
  }; 
// @sect3{The <code>LaplaceProblem</code> class template}  

// 这个主类与  step-6  中的同一类相似。就成员函数而言，唯一增加的是。

// --  <code>assemble_multigrid</code> 的函数，该函数组装了对应于中间层离散运算符的矩阵。

// -  <code>cell_worker</code> 函数，它将我们的PDE集合在一个单元上。

  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const unsigned int degree); 
    void run(); 

  private: 
    template <class Iterator> 
    void cell_worker(const Iterator &  cell, 
                     ScratchData<dim> &scratch_data, 
                     CopyData &        copy_data); 

    void setup_system(); 
    void assemble_system(); 
    void assemble_multigrid(); 
    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    AffineConstraints<double> constraints; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const unsigned int degree; 

// 以下成员是多网格方法的基本数据结构。前四个表示稀疏模式和多级层次结构中各个层次的矩阵，非常类似于上面的全局网格的对象。

// 然后，我们有两个新的矩阵，只需要在自适应网格上进行局部平滑的多网格方法。它们在细化区域的内部和细化边缘之间传递数据，在 @ref mg_paper "多网格论文 "中详细介绍过。

// 最后一个对象存储了每个层次上的边界指数信息和位于两个不同细化层次之间的细化边缘上的指数信息。因此，它的作用与AffineConstraints相似，但在每个层次上。

    MGLevelObject<SparsityPattern> mg_sparsity_patterns; 
    MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns; 

    MGLevelObject<SparseMatrix<double>> mg_matrices; 
    MGLevelObject<SparseMatrix<double>> mg_interface_matrices; 
    MGConstrainedDoFs                   mg_constrained_dofs; 
  }; 
// @sect3{The <code>LaplaceProblem</code> class implementation}  

// 关于三角形的构造函数只有一个简短的评论：按照惯例，deal.II中所有自适应精化的三角形在单元格之间的面的变化不会超过一个级别。然而，对于我们的多网格算法，我们需要一个更严格的保证，即网格在连接两个单元的顶点上的变化也不超过细化级别。换句话说，我们必须防止出现以下情况。

//  @image html limit_level_difference_at_vertices.png ""  

// 这可以通过向三角化类的构造函数传递 Triangulation::limit_level_difference_at_vertices 标志来实现。

  template <int dim> 
  LaplaceProblem<dim>::LaplaceProblem(const unsigned int degree) 
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , fe(degree) 
    , dof_handler(triangulation) 
    , degree(degree) 
  {} 

//  @sect4{LaplaceProblem::setup_system}  

// 除了只是在DoFHandler中分配自由度之外，我们在每一层都做同样的事情。然后，我们按照之前的程序，在叶子网格上设置系统。

  template <int dim> 
  void LaplaceProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (by level: "; 
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
      std::cout << dof_handler.n_dofs(level) 
                << (level == triangulation.n_levels() - 1 ? ")" : ", "); 
    std::cout << std::endl; 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    std::set<types::boundary_id> dirichlet_boundary_ids = {0}; 
    Functions::ZeroFunction<dim> homogeneous_dirichlet_bc; 
    const std::map<types::boundary_id, const Function<dim> *> 
      dirichlet_boundary_functions = { 
        {types::boundary_id(0), &homogeneous_dirichlet_bc}}; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             dirichlet_boundary_functions, 
                                             constraints); 
    constraints.close(); 

    { 
      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints); 
      sparsity_pattern.copy_from(dsp); 
    } 
    system_matrix.reinit(sparsity_pattern); 

// 多网格约束必须被初始化。他们需要知道在哪里规定了Dirichlet边界条件。

    mg_constrained_dofs.clear(); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary_ids); 

// 现在是关于多网格数据结构的事情。首先，我们调整多级对象的大小，以容纳每一级的矩阵和稀疏模式。粗略的级别是零（现在是强制性的，但在未来的修订中可能会改变）。注意，这些函数在这里采取的是一个完整的、包容的范围（而不是一个起始索引和大小），所以最细的级别是 <code>n_levels-1</code>  。我们首先要调整容纳SparseMatrix类的容器的大小，因为它们必须在调整大小时释放它们的SparsityPattern才能被销毁。

    const unsigned int n_levels = triangulation.n_levels(); 

    mg_interface_matrices.resize(0, n_levels - 1); 
    mg_matrices.resize(0, n_levels - 1); 
    mg_sparsity_patterns.resize(0, n_levels - 1); 
    mg_interface_sparsity_patterns.resize(0, n_levels - 1); 

// 现在，我们必须在每个级别上提供一个矩阵。为此，我们首先使用 MGTools::make_sparsity_pattern 函数在每个层次上生成一个初步的压缩稀疏模式（关于这个主题的更多信息，请参见 @ref Sparsity 模块），然后将其复制到我们真正想要的那一个。下一步是用拟合的稀疏度模式初始化接口矩阵。

// 值得指出的是，界面矩阵只包含位于较粗和较细的网格之间的自由度的条目。因此，它们甚至比我们的多网格层次结构中的各个层次的矩阵还要稀疏。因此，我们使用一个专门为此目的而建立的函数来生成它。

    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        { 
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                     dof_handler.n_dofs(level)); 
          MGTools::make_sparsity_pattern(dof_handler, dsp, level); 

          mg_sparsity_patterns[level].copy_from(dsp); 
          mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
        } 
        { 
          DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                     dof_handler.n_dofs(level)); 
          MGTools::make_interface_sparsity_pattern(dof_handler, 
                                                   mg_constrained_dofs, 
                                                   dsp, 
                                                   level); 
          mg_interface_sparsity_patterns[level].copy_from(dsp); 
          mg_interface_matrices[level].reinit( 
            mg_interface_sparsity_patterns[level]); 
        } 
      } 
  } 
// @sect4{LaplaceProblem::cell_worker}  

// cell_worker函数用于在给定的单元上组装矩阵和右手边。这个函数用于活动单元生成system_matrix，并在每个层次上建立层次矩阵。

// 注意，当从assemble_multigrid()调用时，我们也会组装一个右手边，尽管它没有被使用。

  template <int dim> 
  template <class Iterator> 
  void LaplaceProblem<dim>::cell_worker(const Iterator &  cell, 
                                        ScratchData<dim> &scratch_data, 
                                        CopyData &        copy_data) 
  { 
    FEValues<dim> &fe_values = scratch_data.fe_values; 
    fe_values.reinit(cell); 

    const unsigned int dofs_per_cell = fe_values.get_fe().n_dofs_per_cell(); 
    const unsigned int n_q_points    = fe_values.get_quadrature().size(); 

    copy_data.reinit(cell, dofs_per_cell); 

    const std::vector<double> &JxW = fe_values.get_JxW_values(); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const double coefficient = 
          (fe_values.get_quadrature_points()[q][0] < 0.0) ? 1.0 : 0.1; 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              { 
                copy_data.cell_matrix(i, j) += 
                  coefficient * 
                  (fe_values.shape_grad(i, q) * fe_values.shape_grad(j, q)) * 
                  JxW[q]; 
              } 
            copy_data.cell_rhs(i) += 1.0 * fe_values.shape_value(i, q) * JxW[q]; 
          } 
      } 
  } 

//  @sect4{LaplaceProblem::assemble_system}  

// 下面的函数将线性系统集合在网格的活动单元上。为此，我们向Mesh_loop()函数传递两个lambda函数。cell_worker函数重定向到同名的类成员函数，而copyer是这个函数特有的，它使用约束条件将本地矩阵和向量复制到相应的全局矩阵。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    MappingQ1<dim> mapping; 

    auto cell_worker = 
      [&](const typename DoFHandler<dim>::active_cell_iterator &cell, 
          ScratchData<dim> &                                    scratch_data, 
          CopyData &                                            copy_data) { 
        this->cell_worker(cell, scratch_data, copy_data); 
      }; 

    auto copier = [&](const CopyData &cd) { 
      this->constraints.distribute_local_to_global(cd.cell_matrix, 
                                                   cd.cell_rhs, 
                                                   cd.local_dof_indices, 
                                                   system_matrix, 
                                                   system_rhs); 
    }; 

    const unsigned int n_gauss_points = degree + 1; 

    ScratchData<dim> scratch_data(mapping, 
                                  fe, 
                                  n_gauss_points, 
                                  update_values | update_gradients | 
                                    update_JxW_values | 
                                    update_quadrature_points); 

    MeshWorker::mesh_loop(dof_handler.begin_active(), 
                          dof_handler.end(), 
                          cell_worker, 
                          copier, 
                          scratch_data, 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 
  } 
// @sect4{LaplaceProblem::assemble_multigrid}  

// 下一个函数是建立矩阵，定义每一层网格上的多网格方法。集成的核心与上面的相同，但是下面的循环将遍历所有已存在的单元，而不仅仅是活动的单元，并且必须将结果输入正确的层矩阵。幸运的是，MeshWorker对我们隐藏了大部分的内容，因此这个函数和之前的函数的区别只在于汇编器的设置和循环中的不同迭代器。

// 我们为每个层次生成一个AffineConstraints对象，其中包含边界和界面道夫作为约束条目。然后，相应的对象被用来生成层次矩阵。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_multigrid() 
  { 
    MappingQ1<dim>     mapping; 
    const unsigned int n_levels = triangulation.n_levels(); 

    std::vector<AffineConstraints<double>> boundary_constraints(n_levels); 
    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        IndexSet dofset; 
        DoFTools::extract_locally_relevant_level_dofs(dof_handler, 
                                                      level, 
                                                      dofset); 
        boundary_constraints[level].reinit(dofset); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_refinement_edge_indices(level)); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_boundary_indices(level)); 
        boundary_constraints[level].close(); 
      } 

    auto cell_worker = 
      [&](const typename DoFHandler<dim>::level_cell_iterator &cell, 
          ScratchData<dim> &                                   scratch_data, 
          CopyData &                                           copy_data) { 
        this->cell_worker(cell, scratch_data, copy_data); 
      }; 

    auto copier = [&](const CopyData &cd) { 
      boundary_constraints[cd.level].distribute_local_to_global( 
        cd.cell_matrix, cd.local_dof_indices, mg_matrices[cd.level]); 

      const unsigned int dofs_per_cell = cd.local_dof_indices.size(); 

// 接口条目在填充mg_matrices[cd.level]时被上面的boundary_constraints对象所忽略。相反，我们手动将这些条目复制到当前级别的界面矩阵中。

      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        for (unsigned int j = 0; j < dofs_per_cell; ++j) 
          if (mg_constrained_dofs.is_interface_matrix_entry( 
                cd.level, cd.local_dof_indices[i], cd.local_dof_indices[j])) 
            { 
              mg_interface_matrices[cd.level].add(cd.local_dof_indices[i], 
                                                  cd.local_dof_indices[j], 
                                                  cd.cell_matrix(i, j)); 
            } 
    }; 

    const unsigned int n_gauss_points = degree + 1; 

    ScratchData<dim> scratch_data(mapping, 
                                  fe, 
                                  n_gauss_points, 
                                  update_values | update_gradients | 
                                    update_JxW_values | 
                                    update_quadrature_points); 

    MeshWorker::mesh_loop(dof_handler.begin_mg(), 
                          dof_handler.end_mg(), 
                          cell_worker, 
                          copier, 
                          scratch_data, 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 
  } 

//  @sect4{LaplaceProblem::solve}  

// 这是另外一个在支持多栅求解器（或者说，事实上，我们使用多栅方法的前提条件）方面有明显不同的函数。

// 让我们从建立多层次方法的两个组成部分开始：层次间的转移运算器和最粗层次上的求解器。在有限元方法中，转移算子来自所涉及的有限元函数空间，通常可以用独立于所考虑问题的通用方式计算。在这种情况下，我们可以使用MGTransferPrebuilt类，给定最终线性系统的约束和MGConstrainedDoFs对象，该对象知道每个层次的边界条件和不同细化层次之间接口的自由度，可以从具有层次自由度的DoFHandler对象中建立这些转移操作的矩阵。

// 下面几行的第二部分是关于粗略网格求解器的。由于我们的粗网格确实非常粗，我们决定采用直接求解器（最粗层次矩阵的Householder分解），即使其实现不是特别复杂。如果我们的粗网格比这里的5个单元多得多，那么这里显然需要更合适的东西。

  template <int dim> 
  void LaplaceProblem<dim>::solve() 
  { 
    MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs); 
    mg_transfer.build(dof_handler); 

    FullMatrix<double> coarse_matrix; 
    coarse_matrix.copy_from(mg_matrices[0]); 
    MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver; 
    coarse_grid_solver.initialize(coarse_matrix); 

// 多级求解器或预处理器的下一个组成部分是，我们需要在每一级上设置平滑器。这方面常见的选择是使用松弛方法的应用（如SOR、Jacobi或Richardson方法）或求解器方法的少量迭代（如CG或GMRES）。 mg::SmootherRelaxation 和MGSmootherPrecondition类为这两种平滑器提供支持。这里，我们选择应用单一的SOR迭代。为此，我们定义一个适当的别名，然后设置一个平滑器对象。

// 最后一步是用我们的水平矩阵初始化平滑器对象，并设置一些平滑参数。 <code>initialize()</code> 函数可以有选择地接受额外的参数，这些参数将被传递给每一级的平滑器对象。在目前SOR平滑器的情况下，这可能包括一个松弛参数。然而，我们在这里将这些参数保留为默认值。对 <code>set_steps()</code> 的调用表明我们将在每个级别上使用两个前平滑步骤和两个后平滑步骤；为了在不同级别上使用可变数量的平滑器步骤，可以在对 <code>mg_smoother</code> 对象的构造函数调用中设置更多选项。

// 最后一步的结果是我们使用SOR方法作为平滑器的事实

// --这不是对称的

// 但我们在下面使用共轭梯度迭代（需要对称的预处理），我们需要让多级预处理确保我们得到一个对称的算子，即使是非对称的平滑器。

    using Smoother = PreconditionSOR<SparseMatrix<double>>; 
    mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother; 
    mg_smoother.initialize(mg_matrices); 
    mg_smoother.set_steps(2); 
    mg_smoother.set_symmetric(true); 

// 下一个准备步骤是，我们必须将我们的水平和接口矩阵包裹在一个具有所需乘法函数的对象中。我们将为从粗到细的接口对象创建两个对象，反之亦然；多网格算法将在后面的操作中使用转置运算器，允许我们用已经建立的矩阵初始化该运算器的上下版本。

    mg::Matrix<Vector<double>> mg_matrix(mg_matrices); 
    mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices); 
    mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices); 

// 现在，我们准备设置V型循环算子和多级预处理程序。

    Multigrid<Vector<double>> mg( 
      mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother); 
    mg.set_edge_matrices(mg_interface_down, mg_interface_up); 

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
      preconditioner(dof_handler, mg, mg_transfer); 

// 有了这一切，我们终于可以用通常的方法来解决这个线性系统了。

    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    solution = 0; 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 
    std::cout << "   Number of CG iterations: " << solver_control.last_step() 
              << "\n" 
              << std::endl; 
    constraints.distribute(solution); 
  } 

//  @sect4{Postprocessing}  

// 以下两个函数在计算出解决方案后对其进行后处理。特别是，第一个函数在每个周期开始时细化网格，第二个函数在每个周期结束时输出结果。这些函数与  step-6  中的函数几乎没有变化。

  template <int dim> 
  void LaplaceProblem<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(degree + 2), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell); 
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 
    triangulation.execute_coarsening_and_refinement(); 
  } 

  template <int dim> 
  void LaplaceProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
    data_out.write_vtk(output); 
  } 
// @sect4{LaplaceProblem::run}  

// 和上面的几个函数一样，这几乎是对  step-6  中相应函数的复制。唯一的区别是对 <code>assemble_multigrid</code> 的调用，它负责形成我们在多网格方法中需要的每一层的矩阵。

  template <int dim> 
  void LaplaceProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 8; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_ball(triangulation); 
            triangulation.refine_global(2); 
          } 
        else 
          refine_grid(); 

        std::cout << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        assemble_system(); 
        assemble_multigrid(); 

        solve(); 
        output_results(cycle); 
      } 
  } 
} // namespace Step16 
// @sect3{The main() function}  

// 这又是与 step-6 中相同的函数。

int main() 
{ 
  try 
    { 
      using namespace Step16; 

      LaplaceProblem<2> laplace_problem(1); 
      laplace_problem.run(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 

      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 



