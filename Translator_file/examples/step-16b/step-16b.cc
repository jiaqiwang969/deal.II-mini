

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

// 这些，现在，是多级方法所必需的包括。第一个声明了如何处理多网格方法每个层次上的Dirichlet边界条件。对于自由度的实际描述，我们不需要任何新的包含文件，因为DoFHandler已经实现了所有必要的方法。我们只需要将自由度分配给更多的层次。

// 其余的包含文件涉及到作为线性算子（求解器或预处理器）的多重网格的力学问题。

#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

// 最后我们包括MeshWorker框架。这个框架通过其函数loop()和integration_loop()，自动在单元格上进行循环，并将数据组装成向量、矩阵等。它自动服从约束。由于我们必须建立几个矩阵，并且必须注意几组约束，这将使我们省去很多麻烦。

#include <deal.II/meshworker/dof_info.h> 
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/simple.h> 
#include <deal.II/meshworker/output.h> 
#include <deal.II/meshworker/loop.h> 

// 为了节省精力，我们使用了在以下文件中找到的预先实现的拉普拉斯。

#include <deal.II/integrators/laplace.h> 
#include <deal.II/integrators/l2.h> 

// 这就是C++。

#include <iostream> 
#include <fstream> 

using namespace dealii; 

namespace Step16 
{ 
// @sect3{The integrator on each cell}  

//  MeshWorker::integration_loop() 希望有一个类能够提供在单元格和边界及内部面的积分功能。这是由下面的类来完成的。在构造函数中，我们告诉循环应该计算单元格积分（"真"），但不应该计算边界和内部面的积分（两个 "假"）。因此，我们只需要一个单元格函数，而不需要面的函数。

  template <int dim> 
  class LaplaceIntegrator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    LaplaceIntegrator(); 
    virtual void cell(MeshWorker::DoFInfo<dim> &        dinfo, 
                      MeshWorker::IntegrationInfo<dim> &info) const override; 
  }; 

  template <int dim> 
  LaplaceIntegrator<dim>::LaplaceIntegrator() 
    : MeshWorker::LocalIntegrator<dim>(true, false, false) 
  {} 

// 接下来是每个单元上的实际积分器。我们解决一个泊松问题，在右半平面上的系数为1，在左半平面上的系数为十分之一。

//  MeshWorker::LocalResults 的基类 MeshWorker::DoFInfo 包含可以在这个局部积分器中填充的对象。在MeshWorker框架内，有多少对象被创建是由装配器类决定的。在这里，我们举例测试一下，需要一个矩阵 (MeshWorker::LocalResults::n_matrices()).  矩阵是通过 MeshWorker::LocalResults::matrix(), 来访问的，它的第一个参数是矩阵的编号。第二个参数只用于面的积分，当每个测试函数使用两个矩阵时。那么，第二个指标为 "true "的矩阵将以相同的索引存在。

//  MeshWorker::IntegrationInfo 提供了一个或几个FEValues对象，下面这些对象被 LocalIntegrators::Laplace::cell_matrix() 或 LocalIntegrators::L2::L2(). 使用，因为我们只组装一个PDE，所以也只有一个索引为0的对象。

// 此外，我们注意到这个积分器的作用是计算多级预处理的矩阵，以及全局系统的矩阵和右手边。由于系统的汇编器需要一个额外的向量， MeshWorker::LocalResults::n_vectors() 要返回一个非零值。相应地，我们在这个函数的末尾填充了一个右边的向量。由于LocalResults可以处理多个BlockVector对象，但我们这里又是最简单的情况，所以我们将信息输入到零号向量的零号块中。

  template <int dim> 
  void 
  LaplaceIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &        dinfo, 
                               MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    AssertDimension(dinfo.n_matrices(), 1); 
    const double coefficient = (dinfo.cell->center()(0) > 0.) ? .1 : 1.; 

    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix, 
                                           info.fe_values(0), 
                                           coefficient); 

    if (dinfo.n_vectors() > 0) 
      { 
        std::vector<double> rhs(info.fe_values(0).n_quadrature_points, 1.); 
        LocalIntegrators::L2::L2(dinfo.vector(0).block(0), 
                                 info.fe_values(0), 
                                 rhs); 
      } 
  } 
// @sect3{The <code>LaplaceProblem</code> class template}  

// 这个主类与  step-6  中的类基本相同。就成员函数而言，唯一增加的是 <code>assemble_multigrid</code> 函数，它组装了对应于中间层离散运算符的矩阵。

  template <int dim> 
  class LaplaceProblem 
  { 
  public: 
    LaplaceProblem(const unsigned int degree); 
    void run(); 

  private: 
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

// 以下成员是多网格方法的基本数据结构。前两个表示稀疏模式和多级层次结构中各个层次的矩阵，非常类似于上面的全局网格的对象。

// 然后，我们有两个新的矩阵，只需要在自适应网格上进行局部平滑的多网格方法。它们在细化区域的内部和细化边缘之间传递数据，在 @ref mg_paper "多网格论文 "中详细介绍过。

// 最后一个对象存储了每个层次上的边界指数信息和位于两个不同细化层次之间的细化边缘上的指数信息。因此，它的作用与AffineConstraints类似，但在每个层次上。

    MGLevelObject<SparsityPattern>      mg_sparsity_patterns; 
    MGLevelObject<SparseMatrix<double>> mg_matrices; 
    MGLevelObject<SparseMatrix<double>> mg_interface_in; 
    MGLevelObject<SparseMatrix<double>> mg_interface_out; 
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

    deallog << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << " (by level: "; 
    for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
      deallog << dof_handler.n_dofs(level) 
              << (level == triangulation.n_levels() - 1 ? ")" : ", "); 
    deallog << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

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
    constraints.condense(dsp); 
    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 

// 多网格约束必须被初始化。他们也需要知道边界值，所以我们也在这里传递 <code>dirichlet_boundary</code> 。

    mg_constrained_dofs.clear(); 
    mg_constrained_dofs.initialize(dof_handler); 
    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, 
                                                       dirichlet_boundary_ids); 

// 现在是关于多网格数据结构的事情。首先，我们调整多级对象的大小，以容纳每一级的矩阵和稀疏模式。粗略的级别是零（现在是强制性的，但在未来的修订中可能会改变）。注意，这些函数在这里采取的是一个完整的、包容的范围（而不是一个起始索引和大小），所以最细的级别是 <code>n_levels-1</code>  。我们首先要调整容纳SparseMatrix类的容器的大小，因为它们必须在调整大小时释放它们的SparsityPattern才能被销毁。

    const unsigned int n_levels = triangulation.n_levels(); 

    mg_interface_in.resize(0, n_levels - 1); 
    mg_interface_in.clear_elements(); 
    mg_interface_out.resize(0, n_levels - 1); 
    mg_interface_out.clear_elements(); 
    mg_matrices.resize(0, n_levels - 1); 
    mg_matrices.clear_elements(); 
    mg_sparsity_patterns.resize(0, n_levels - 1); 

// 现在，我们必须在每个层面上提供一个矩阵。为此，我们首先使用 MGTools::make_sparsity_pattern 函数在每个层次上生成一个初步的压缩稀疏模式（关于这个主题的更多信息，请参见 @ref Sparsity 模块），然后把它复制到我们真正想要的那个层次上。下一步是用这些稀疏模式初始化两种层次矩阵。

// 值得指出的是，界面矩阵只有位于较粗的网格和较细的网格之间的界面上的自由度条目。因此，它们甚至比我们多网格层次结构中的各个层次的矩阵还要稀少。如果我们更关心内存的使用（可能还有我们使用这些矩阵的速度），我们应该对这两种矩阵使用不同的稀疏性模式。

    for (unsigned int level = 0; level < n_levels; ++level) 
      { 
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level), 
                                   dof_handler.n_dofs(level)); 
        MGTools::make_sparsity_pattern(dof_handler, dsp, level); 

        mg_sparsity_patterns[level].copy_from(dsp); 

        mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
        mg_interface_in[level].reinit(mg_sparsity_patterns[level]); 
        mg_interface_out[level].reinit(mg_sparsity_patterns[level]); 
      } 
  } 
// @sect4{LaplaceProblem::assemble_system}  

// 下面的函数将线性系统装配在网格的最细层上。由于我们想在下面的层次装配中重用这里的代码，我们使用本地积分器类LaplaceIntegrator，而将循环留给MeshWorker框架。因此，这个函数首先设置了这个框架所需的对象，即  

// - 一个 MeshWorker::IntegrationInfoBox 对象，它将提供单元格上正交点的所有需要的数据。这个对象可以看作是FEValues的扩展，提供更多的有用信息。 

// - 一个 MeshWorker::DoFInfo 对象，它一方面扩展了单元格迭代器的功能，另一方面也为其基类LocalResults的返回值提供了空间。 

// - 一个汇编器，在这里是指整个系统。这里的 "简单 "指的是全局系统没有一个块状结构。 

// - 本地集成器，它实现了实际的形式。

// 在循环将所有这些组合成一个矩阵和一个右手边之后，还有一件事要做：集合器对受限自由度的矩阵行和列不做任何处理。因此，我们在对角线上放一个一，使整个系统摆好。一的值或任何固定的值都有一个好处，即它对矩阵的频谱的影响很容易理解。由于相应的特征向量形成了一个不变的子空间，所选择的值不会影响Krylov空间求解器的收敛性。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_system() 
  { 
    MappingQ1<dim>                      mapping; 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags                         update_flags = 
      update_values | update_gradients | update_hessians; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>> 
      assembler; 
    assembler.initialize(constraints); 
    assembler.initialize(system_matrix, system_rhs); 

    LaplaceIntegrator<dim> matrix_integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           matrix_integrator, 
                                           assembler); 

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
      if (constraints.is_constrained(i)) 
        system_matrix.set(i, i, 1.); 
  } 
// @sect4{LaplaceProblem::assemble_multigrid}  

// 下一个函数是建立线性算子（矩阵），定义每一级网格上的多栅方法。积分的核心和上面的一样，但是下面的循环会遍历所有已有的单元，而不仅仅是活动的单元，而且结果必须输入正确的层次矩阵。幸运的是，MeshWorker对我们隐藏了大部分的内容，因此这个函数和之前的函数的区别只在于汇编器的设置和循环中不同的迭代器。另外，最后修复矩阵的过程也比较复杂。

  template <int dim> 
  void LaplaceProblem<dim>::assemble_multigrid() 
  { 
    MappingQ1<dim>                      mapping; 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags                         update_flags = 
      update_values | update_gradients | update_hessians; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler; 
    assembler.initialize(mg_constrained_dofs); 
    assembler.initialize(mg_matrices); 
    assembler.initialize_interfaces(mg_interface_in, mg_interface_out); 

    LaplaceIntegrator<dim> matrix_integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(), 
                                           dof_handler.end_mg(), 
                                           dof_info, 
                                           info_box, 
                                           matrix_integrator, 
                                           assembler); 

    const unsigned int nlevels = triangulation.n_levels(); 
    for (unsigned int level = 0; level < nlevels; ++level) 
      { 
        for (unsigned int i = 0; i < dof_handler.n_dofs(level); ++i) 
          if (mg_constrained_dofs.is_boundary_index(level, i) || 
              mg_constrained_dofs.at_refinement_edge(level, i)) 
            mg_matrices[level].set(i, i, 1.); 
      } 
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

// 多级求解器或预处理器的下一个组成部分是，我们需要在每一级上有一个平滑器。这方面常见的选择是使用松弛方法的应用（如SOR、Jacobi或Richardson方法）或求解器方法的少量迭代（如CG或GMRES）。 mg::SmootherRelaxation 和MGSmootherPrecondition类为这两种平滑器提供支持。这里，我们选择应用单一的SOR迭代。为此，我们定义一个适当的别名，然后设置一个平滑器对象。

// 最后一步是用我们的水平矩阵初始化平滑器对象，并设置一些平滑参数。 <code>initialize()</code> 函数可以有选择地接受额外的参数，这些参数将被传递给每一级的平滑器对象。在当前SOR平滑器的情况下，这可能包括一个松弛参数。然而，我们在这里将这些参数保留为默认值。对 <code>set_steps()</code> 的调用表明我们将在每个级别上使用两个前平滑步骤和两个后平滑步骤；为了在不同级别上使用可变数量的平滑器步骤，可以在对 <code>mg_smoother</code> 对象的构造函数调用中设置更多选项。

// 最后一步的结果是我们使用SOR方法作为平滑器的事实

// --这不是对称的

// 但我们在下面使用共轭梯度迭代（需要对称的预处理），我们需要让多级预处理确保我们得到一个对称的算子，即使是非对称的平滑器。

    using Smoother = PreconditionSOR<SparseMatrix<double>>; 
    mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother; 
    mg_smoother.initialize(mg_matrices); 
    mg_smoother.set_steps(2); 
    mg_smoother.set_symmetric(true); 

// 下一个准备步骤是，我们必须将我们的水平矩阵和接口矩阵包裹在一个具有所需乘法函数的对象中。我们将为从粗到细的接口对象创建两个对象，反之亦然；多网格算法将在以后的操作中使用转置运算器，允许我们用已经建立的矩阵初始化该运算器的上下版本。

    mg::Matrix<Vector<double>> mg_matrix(mg_matrices); 
    mg::Matrix<Vector<double>> mg_interface_up(mg_interface_in); 
    mg::Matrix<Vector<double>> mg_interface_down(mg_interface_out); 

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
    constraints.distribute(solution); 
  } 

//  @sect4{Postprocessing}  

// 下面两个函数在计算出解决方案后对其进行后处理。特别是，第一个函数在每个周期开始时细化网格，第二个函数在每个周期结束时输出结果。这些函数与 step-6 中的函数几乎没有变化，只有一个小的区别：我们以VTK格式生成输出，以使用当今更现代的可视化程序，而不是 step-6 编写时的那些。

  template <int dim> 
  void LaplaceProblem<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
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
        deallog << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_ball(triangulation); 
            triangulation.refine_global(1); 
          } 
        else 
          refine_grid(); 

        deallog << "   Number of active cells:       " 
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

      deallog.depth_console(2); 

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

