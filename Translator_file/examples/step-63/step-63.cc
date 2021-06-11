


/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
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
 * Authors: Thomas C. Clevenger, Clemson University 
 *          Timo Heister, Clemson University and University of Utah 
 */ 


// @sect3{Include files}  

// 标准deal.II需要的典型文件。

#include <deal.II/base/tensor_function.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/parameter_handler.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/relaxation_block.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_out.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/mapping_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

// 包括所有相关的多层次文件。

#include <deal.II/multigrid/mg_constrained_dofs.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 
#include <deal.II/multigrid/mg_matrix.h> 

// C++:

#include <algorithm> 
#include <fstream> 
#include <iostream> 
#include <random> 

// 我们将使用 MeshWorker::mesh_loop 功能来组装矩阵。

#include <deal.II/meshworker/mesh_loop.h> 
// @sect3{MeshWorker data}  

// 像往常一样，我们将把所有与这个程序有关的东西放到一个自己的命名空间中。

// 由于我们将使用MeshWorker框架，第一步是定义以下由 MeshWorker::mesh_loop(): 使用的assemble_cell()函数所需要的结构 `ScratchData`包含一个FEValues对象，这是组装一个单元的局部贡献所需要的，而`CopyData`包含一个单元的局部贡献的输出和复制到全局系统的必要信息。它们的目的在WorkStream类的文档中也有解释）。

namespace Step63 
{ 
  using namespace dealii; 

  template <int dim> 
  struct ScratchData 
  { 
    ScratchData(const FiniteElement<dim> &fe, 
                const unsigned int        quadrature_degree) 
      : fe_values(fe, 
                  QGauss<dim>(quadrature_degree), 
                  update_values | update_gradients | update_hessians | 
                    update_quadrature_points | update_JxW_values) 
    {} 

    ScratchData(const ScratchData<dim> &scratch_data) 
      : fe_values(scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  update_values | update_gradients | update_hessians | 
                    update_quadrature_points | update_JxW_values) 
    {} 

    FEValues<dim> fe_values; 
  }; 

  struct CopyData 
  { 
    CopyData() = default; 

    unsigned int level; 
    unsigned int dofs_per_cell; 

    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_rhs; 
    std::vector<types::global_dof_index> local_dof_indices; 
  }; 

//  @sect3{Problem parameters}  

// 第二步是定义处理要从输入文件中读取的运行时参数的类。

// 我们将使用ParameterHandler在运行时传入参数。结构`Settings`解析并存储整个程序要查询的参数。

  struct Settings 
  { 
    enum DoFRenumberingStrategy 
    { 
      none, 
      downstream, 
      upstream, 
      random 
    }; 

    void get_parameters(const std::string &prm_filename); 

    double                 epsilon; 
    unsigned int           fe_degree; 
    std::string            smoother_type; 
    unsigned int           smoothing_steps; 
    DoFRenumberingStrategy dof_renumbering; 
    bool                   with_streamline_diffusion; 
    bool                   output; 
  }; 

  void Settings::get_parameters(const std::string &prm_filename) 
  { 

/* 首先声明参数...   */ 

 
    ParameterHandler prm; 

    prm.declare_entry("Epsilon", 
                      "0.005", 
                      Patterns::Double(0), 
                      "Diffusion parameter"); 

    prm.declare_entry("Fe degree", 
                      "1", 
                      Patterns::Integer(1), 
                      "Finite Element degree"); 
    prm.declare_entry("Smoother type", 
                      "block SOR", 
                      Patterns::Selection("SOR|Jacobi|block SOR|block Jacobi"), 
                      "Select smoother: SOR|Jacobi|block SOR|block Jacobi"); 
    prm.declare_entry("Smoothing steps", 
                      "2", 
                      Patterns::Integer(1), 
                      "Number of smoothing steps"); 
    prm.declare_entry( 
      "DoF renumbering", 
      "downstream", 
      Patterns::Selection("none|downstream|upstream|random"), 
      "Select DoF renumbering: none|downstream|upstream|random"); 
    prm.declare_entry("With streamline diffusion", 
                      "true", 
                      Patterns::Bool(), 
                      "Enable streamline diffusion stabilization: true|false"); 
    prm.declare_entry("Output", 
                      "true", 
                      Patterns::Bool(), 
                      "Generate graphical output: true|false"); 
    /* ...然后尝试从输入文件中读取它们的值。  */ 
    if (prm_filename.empty()) 
      { 
        prm.print_parameters(std::cout, ParameterHandler::Text); 
        AssertThrow( 
          false, ExcMessage("Please pass a .prm file as the first argument!")); 
      } 

    prm.parse_input(prm_filename); 

    epsilon         = prm.get_double("Epsilon"); 
    fe_degree       = prm.get_integer("Fe degree"); 
    smoother_type   = prm.get("Smoother type"); 
    smoothing_steps = prm.get_integer("Smoothing steps"); 

    const std::string renumbering = prm.get("DoF renumbering"); 
    if (renumbering == "none") 
      dof_renumbering = DoFRenumberingStrategy::none; 
    else if (renumbering == "downstream") 
      dof_renumbering = DoFRenumberingStrategy::downstream; 
    else if (renumbering == "upstream") 
      dof_renumbering = DoFRenumberingStrategy::upstream; 
    else if (renumbering == "random") 
      dof_renumbering = DoFRenumberingStrategy::random; 
    else 
      AssertThrow(false, 
                  ExcMessage("The <DoF renumbering> parameter has " 
                             "an invalid value.")); 

    with_streamline_diffusion = prm.get_bool("With streamline diffusion"); 
    output                    = prm.get_bool("Output"); 
  } 
// @sect3{Cell permutations}  

// 遍历单元和自由度的顺序将对乘法的收敛速度起作用。在这里，我们定义了一些函数，这些函数返回单元格的特定顺序，供块平滑器使用。

// 对于每种类型的单元格排序，我们定义了一个用于活动网格的函数和一个用于水平网格的函数（即用于多网格层次结构中的某一层的单元格）。虽然求解系统所需的唯一重新排序是在水平网格上进行的，但为了可视化的目的，我们在output_results()中包含了主动网格的重新排序。

// 对于两个下游排序函数，我们首先创建一个包含所有相关单元的数组，然后使用一个 "比较器 "对象在下游方向进行排序。然后，函数的输出是一个简单的数组，包含了刚刚计算出来的单元格的索引。

  template <int dim> 
  std::vector<unsigned int> 
  create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler, 
                                  const Tensor<1, dim>   direction, 
                                  const unsigned int     level) 
  { 
    std::vector<typename DoFHandler<dim>::level_cell_iterator> ordered_cells; 
    ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level)); 
    for (const auto &cell : dof_handler.cell_iterators_on_level(level)) 
      ordered_cells.push_back(cell); 

    const DoFRenumbering:: 
      CompareDownstream<typename DoFHandler<dim>::level_cell_iterator, dim> 
        comparator(direction); 
    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator); 

    std::vector<unsigned> ordered_indices; 
    ordered_indices.reserve(dof_handler.get_triangulation().n_cells(level)); 

    for (const auto &cell : ordered_cells) 
      ordered_indices.push_back(cell->index()); 

    return ordered_indices; 
  } 

  template <int dim> 
  std::vector<unsigned int> 
  create_downstream_cell_ordering(const DoFHandler<dim> &dof_handler, 
                                  const Tensor<1, dim>   direction) 
  { 
    std::vector<typename DoFHandler<dim>::active_cell_iterator> ordered_cells; 
    ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells()); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      ordered_cells.push_back(cell); 

    const DoFRenumbering:: 
      CompareDownstream<typename DoFHandler<dim>::active_cell_iterator, dim> 
        comparator(direction); 
    std::sort(ordered_cells.begin(), ordered_cells.end(), comparator); 

    std::vector<unsigned int> ordered_indices; 
    ordered_indices.reserve(dof_handler.get_triangulation().n_active_cells()); 

    for (const auto &cell : ordered_cells) 
      ordered_indices.push_back(cell->index()); 

    return ordered_indices; 
  } 

// 产生随机排序的函数在精神上是相似的，它们首先将所有单元的信息放入一个数组。但是，它们不是对它们进行排序，而是利用C++提供的生成随机数的设施对元素进行随机洗牌。这样做的方式是在数组的所有元素上进行迭代，为之前的另一个元素抽取一个随机数，然后交换这些元素。其结果是对数组中的元素进行随机洗牌。

  template <int dim> 
  std::vector<unsigned int> 
  create_random_cell_ordering(const DoFHandler<dim> &dof_handler, 
                              const unsigned int     level) 
  { 
    std::vector<unsigned int> ordered_cells; 
    ordered_cells.reserve(dof_handler.get_triangulation().n_cells(level)); 
    for (const auto &cell : dof_handler.cell_iterators_on_level(level)) 
      ordered_cells.push_back(cell->index()); 

    std::mt19937 random_number_generator; 
    std::shuffle(ordered_cells.begin(), 
                 ordered_cells.end(), 
                 random_number_generator); 

    return ordered_cells; 
  } 

  template <int dim> 
  std::vector<unsigned int> 
  create_random_cell_ordering(const DoFHandler<dim> &dof_handler) 
  { 
    std::vector<unsigned int> ordered_cells; 
    ordered_cells.reserve(dof_handler.get_triangulation().n_active_cells()); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      ordered_cells.push_back(cell->index()); 

    std::mt19937 random_number_generator; 
    std::shuffle(ordered_cells.begin(), 
                 ordered_cells.end(), 
                 random_number_generator); 

    return ordered_cells; 
  } 
// @sect3{Right-hand side and boundary values}  

// 本教程中所解决的问题是对<a
//  href="https:global.oup.com/academic/product/finite-elements-and-fast-iterative-solvers-9780199678808">
//  Finite Elements and Fast Iterative Solvers: with Applications in
//  Incompressible Fluid Dynamics by Elman, Silvester, and Wathen</a>第118页上的例3.1.3的修改。主要的区别是我们在域的中心增加了一个洞，其边界条件为零的Dirichlet。

// 为了获得完整的描述，我们需要首先实现零右手边的类（当然，我们可以直接使用 Functions::ZeroFunction):  。
  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> &, 
                                   const unsigned int component) const 
  { 
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 
    (void)component; 

    return 0.0; 
  } 

  template <int dim> 
  void RightHandSide<dim>::value_list(const std::vector<Point<dim>> &points, 
                                      std::vector<double> &          values, 
                                      const unsigned int component) const 
  { 
    Assert(values.size() == points.size(), 
           ExcDimensionMismatch(values.size(), points.size())); 

    for (unsigned int i = 0; i < points.size(); ++i) 
      values[i] = RightHandSide<dim>::value(points[i], component); 
  } 

// 我们也有迪里希特的边界条件。在外部正方形边界的连接部分，我们将数值设置为1，其他地方（包括内部圆形边界）的数值设置为0。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & p, 
                                    const unsigned int component) const 
  { 
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 
    (void)component; 

// 如果  $x=1$  ，或如果  $x>0.5$  和  $y=-1$  ，则将边界设为 1。

    if (std::fabs(p[0] - 1) < 1e-8 || 
        (std::fabs(p[1] + 1) < 1e-8 && p[0] >= 0.5)) 
      { 
        return 1.0; 
      } 
    else 
      { 
        return 0.0; 
      } 
  } 

  template <int dim> 
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points, 
                                       std::vector<double> &          values, 
                                       const unsigned int component) const 
  { 
    Assert(values.size() == points.size(), 
           ExcDimensionMismatch(values.size(), points.size())); 

    for (unsigned int i = 0; i < points.size(); ++i) 
      values[i] = BoundaryValues<dim>::value(points[i], component); 
  } 

//  @sect3{Streamline diffusion implementation}  

// 流水线扩散方法有一个稳定常数，我们需要能够计算出来。这个参数的计算方式的选择取自于<a
//  href="https:link.springer.com/chapter/10.1007/978-3-540-34288-5_27">On
//  Discontinuity-Capturing Methods for Convection-Diffusion
//  Equations by Volker John and Petr Knobloch</a>。

  template <int dim> 
  double compute_stabilization_delta(const double         hk, 
                                     const double         eps, 
                                     const Tensor<1, dim> dir, 
                                     const double         pk) 
  { 
    const double Peclet = dir.norm() * hk / (2.0 * eps * pk); 
    const double coth = 
      (1.0 + std::exp(-2.0 * Peclet)) / (1.0 - std::exp(-2.0 * Peclet)); 

    return hk / (2.0 * dir.norm() * pk) * (coth - 1.0 / Peclet); 
  } 
// @sect3{<code>AdvectionProlem</code> class}  

// 这是程序的主类，看起来应该与  step-16  非常相似。主要的区别是，由于我们是在运行时定义我们的多网格平滑器，我们选择定义一个函数`create_smoother()`和一个类对象`mg_smoother`，这是一个  `std::unique_ptr`  派生于MGSmoother的平滑器。请注意，对于从RelaxationBlock派生的平滑器，我们必须为每个级别包括一个`smoother_data`对象。这将包含关于单元格排序和单元格矩阵倒置方法的信息。

  template <int dim> 
  class AdvectionProblem 
  { 
  public: 
    AdvectionProblem(const Settings &settings); 
    void run(); 

  private: 
    void setup_system(); 

    template <class IteratorType> 
    void assemble_cell(const IteratorType &cell, 
                       ScratchData<dim> &  scratch_data, 
                       CopyData &          copy_data); 
    void assemble_system_and_multigrid(); 

    void setup_smoother(); 

    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim> triangulation; 
    DoFHandler<dim>    dof_handler; 

    const FE_Q<dim>     fe; 
    const MappingQ<dim> mapping; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    MGLevelObject<SparsityPattern> mg_sparsity_patterns; 
    MGLevelObject<SparsityPattern> mg_interface_sparsity_patterns; 

    MGLevelObject<SparseMatrix<double>> mg_matrices; 
    MGLevelObject<SparseMatrix<double>> mg_interface_in; 
    MGLevelObject<SparseMatrix<double>> mg_interface_out; 

    mg::Matrix<Vector<double>> mg_matrix; 
    mg::Matrix<Vector<double>> mg_interface_matrix_in; 
    mg::Matrix<Vector<double>> mg_interface_matrix_out; 

    std::unique_ptr<MGSmoother<Vector<double>>> mg_smoother; 

    using SmootherType = 
      RelaxationBlock<SparseMatrix<double>, double, Vector<double>>; 
    using SmootherAdditionalDataType = SmootherType::AdditionalData; 
    MGLevelObject<SmootherAdditionalDataType> smoother_data; 

    MGConstrainedDoFs mg_constrained_dofs; 

    Tensor<1, dim> advection_direction; 

    const Settings settings; 
  }; 

  template <int dim> 
  AdvectionProblem<dim>::AdvectionProblem(const Settings &settings) 
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , dof_handler(triangulation) 
    , fe(settings.fe_degree) 
    , mapping(settings.fe_degree) 
    , settings(settings) 
  { 
    advection_direction[0] = -std::sin(numbers::PI / 6.0); 
    if (dim >= 2) 
      advection_direction[1] = std::cos(numbers::PI / 6.0); 
    if (dim >= 3) 
      AssertThrow(false, ExcNotImplemented()); 
  } 
// @sect4{<code>AdvectionProblem::setup_system()</code>}  

// 在这里，我们首先为活动和多网格级别的网格设置DoFHandler、AffineConstraints和SparsityPattern对象。

// 我们可以用DoFRenumbering类对活动DoF进行重新编号，但是平滑器只作用于多网格层，因此，这对计算并不重要。相反，我们将对每个多网格层的DoFs进行重新编号。

  template <int dim> 
  void AdvectionProblem<dim>::setup_system() 
  { 
    const unsigned int n_levels = triangulation.n_levels(); 

    dof_handler.distribute_dofs(fe); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    VectorTools::interpolate_boundary_values( 
      mapping, dof_handler, 0, BoundaryValues<dim>(), constraints); 
    VectorTools::interpolate_boundary_values( 
      mapping, dof_handler, 1, BoundaryValues<dim>(), constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ false);

    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 

    dof_handler.distribute_mg_dofs(); 

// 在列举了全局自由度以及（上面最后一行）水平自由度之后，让我们对水平自由度进行重新编号，以获得一个更好的平滑器，正如介绍中所解释的。 如果需要的话，下面的第一个区块会对下游或上游方向的每个层次的自由度进行重新编号。这只对点平滑器（SOR和Jacobi）有必要，因为块平滑器是在单元上操作的（见`create_smoother()`）。然后，下面的块也实现了随机编号。

    if (settings.smoother_type == "SOR" || settings.smoother_type == "Jacobi") 
      { 
        if (settings.dof_renumbering == 
              Settings::DoFRenumberingStrategy::downstream || 
            settings.dof_renumbering == 
              Settings::DoFRenumberingStrategy::upstream) 
          { 
            const Tensor<1, dim> direction = 
              (settings.dof_renumbering == 
                   Settings::DoFRenumberingStrategy::upstream ? 
                 -1.0 : 
                 1.0) * 
              advection_direction; 

            for (unsigned int level = 0; level < n_levels; ++level) 
              DoFRenumbering::downstream(dof_handler, 
                                         level, 
                                         direction, 
                                         /*dof_wise_renumbering =  */ true);

          } 
        else if (settings.dof_renumbering == 
                 Settings::DoFRenumberingStrategy::random) 
          { 
            for (unsigned int level = 0; level < n_levels; ++level) 
              DoFRenumbering::random(dof_handler, level); 
          } 
        else 
          Assert(false, ExcNotImplemented()); 
      } 

// 该函数的其余部分只是设置了数据结构。下面代码的最后几行与其他GMG教程不同，因为它同时设置了接口输入和输出矩阵。我们需要这样做，因为我们的问题是非对称性的。

    mg_constrained_dofs.clear(); 
    mg_constrained_dofs.initialize(dof_handler); 

    mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, {0, 1}); 

    mg_matrices.resize(0, n_levels - 1); 
    mg_matrices.clear_elements(); 
    mg_interface_in.resize(0, n_levels - 1); 
    mg_interface_in.clear_elements(); 
    mg_interface_out.resize(0, n_levels - 1); 
    mg_interface_out.clear_elements(); 
    mg_sparsity_patterns.resize(0, n_levels - 1); 
    mg_interface_sparsity_patterns.resize(0, n_levels - 1); 

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

          mg_interface_in[level].reinit(mg_interface_sparsity_patterns[level]); 
          mg_interface_out[level].reinit(mg_interface_sparsity_patterns[level]); 
        } 
      } 
  } 
// @sect4{<code>AdvectionProblem::assemble_cell()</code>}  

// 这里我们定义了每个单元上的线性系统的装配，以便被下面的Mesh_loop()函数使用。这个函数为活动单元或水平单元（不管它的第一个参数是什么）装配单元矩阵，并且只有在调用活动单元时才装配右手边。

  template <int dim> 
  template <class IteratorType> 
  void AdvectionProblem<dim>::assemble_cell(const IteratorType &cell, 
                                            ScratchData<dim> &  scratch_data, 
                                            CopyData &          copy_data) 
  { 
    copy_data.level = cell->level(); 

    const unsigned int dofs_per_cell = 
      scratch_data.fe_values.get_fe().n_dofs_per_cell(); 
    copy_data.dofs_per_cell = dofs_per_cell; 
    copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 

    const unsigned int n_q_points = 
      scratch_data.fe_values.get_quadrature().size(); 

    if (cell->is_level_cell() == false) 
      copy_data.cell_rhs.reinit(dofs_per_cell); 

    copy_data.local_dof_indices.resize(dofs_per_cell); 
    cell->get_active_or_mg_dof_indices(copy_data.local_dof_indices); 

    scratch_data.fe_values.reinit(cell); 

    RightHandSide<dim>  right_hand_side; 
    std::vector<double> rhs_values(n_q_points); 

    right_hand_side.value_list(scratch_data.fe_values.get_quadrature_points(), 
                               rhs_values); 

// 如果我们使用流线扩散，我们必须把它的贡献加到单元格矩阵和单元格的右手边。如果我们不使用流线扩散，设置 $\delta=0$ 就可以否定这个贡献，我们就可以使用标准的Galerkin有限元组合。

    const double delta = (settings.with_streamline_diffusion ? 
                            compute_stabilization_delta(cell->diameter(), 
                                                        settings.epsilon, 
                                                        advection_direction, 
                                                        settings.fe_degree) : 
                            0.0); 

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            { 

// 本地矩阵的组装有两个部分。首先是Galerkin贡献。

              copy_data.cell_matrix(i, j) += 
                (settings.epsilon * 
                 scratch_data.fe_values.shape_grad(i, q_point) * 
                 scratch_data.fe_values.shape_grad(j, q_point) * 
                 scratch_data.fe_values.JxW(q_point)) + 
                (scratch_data.fe_values.shape_value(i, q_point) * 
                 (advection_direction * 
                  scratch_data.fe_values.shape_grad(j, q_point)) * 
                 scratch_data.fe_values.JxW(q_point)) 

//然后是流线扩散贡献。

                + delta * 
                    (advection_direction * 
                     scratch_data.fe_values.shape_grad(j, q_point)) * 
                    (advection_direction * 
                     scratch_data.fe_values.shape_grad(i, q_point)) * 
                    scratch_data.fe_values.JxW(q_point) - 
                delta * settings.epsilon * 
                  trace(scratch_data.fe_values.shape_hessian(j, q_point)) * 
                  (advection_direction * 
                   scratch_data.fe_values.shape_grad(i, q_point)) * 
                  scratch_data.fe_values.JxW(q_point); 
            } 
          if (cell->is_level_cell() == false) 
            { 

// 同样的情况也适用于右手边。首先是Galerkin贡献。

              copy_data.cell_rhs(i) += 
                scratch_data.fe_values.shape_value(i, q_point) * 
                  rhs_values[q_point] * scratch_data.fe_values.JxW(q_point) 

// 然后是流线扩散贡献。

                + delta * rhs_values[q_point] * advection_direction * 
                    scratch_data.fe_values.shape_grad(i, q_point) * 
                    scratch_data.fe_values.JxW(q_point); 
            } 
        } 
  } 
// @sect4{<code>AdvectionProblem::assemble_system_and_multigrid()</code>}  

// 这里我们采用 MeshWorker::mesh_loop() 来翻阅单元格，为我们组装system_matrix、system_rhs和所有mg_matrices。

  template <int dim> 
  void AdvectionProblem<dim>::assemble_system_and_multigrid() 
  { 
    const auto cell_worker_active = 
      [&](const decltype(dof_handler.begin_active()) &cell, 
          ScratchData<dim> &                          scratch_data, 
          CopyData &                                  copy_data) { 
        this->assemble_cell(cell, scratch_data, copy_data); 
      }; 

    const auto copier_active = [&](const CopyData &copy_data) { 
      constraints.distribute_local_to_global(copy_data.cell_matrix, 
                                             copy_data.cell_rhs, 
                                             copy_data.local_dof_indices, 
                                             system_matrix, 
                                             system_rhs); 
    }; 

    MeshWorker::mesh_loop(dof_handler.begin_active(), 
                          dof_handler.end(), 
                          cell_worker_active, 
                          copier_active, 
                          ScratchData<dim>(fe, fe.degree + 1), 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 

// 与活动层的约束不同，我们选择在这个函数的本地为每个多网格层创建约束对象，因为它们在程序的其他地方从来不需要。

    std::vector<AffineConstraints<double>> boundary_constraints( 
      triangulation.n_global_levels()); 
    for (unsigned int level = 0; level < triangulation.n_global_levels(); 
         ++level) 
      { 
        IndexSet locally_owned_level_dof_indices; 
        DoFTools::extract_locally_relevant_level_dofs( 
          dof_handler, level, locally_owned_level_dof_indices); 
        boundary_constraints[level].reinit(locally_owned_level_dof_indices); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_refinement_edge_indices(level)); 
        boundary_constraints[level].add_lines( 
          mg_constrained_dofs.get_boundary_indices(level)); 
        boundary_constraints[level].close(); 
      } 

    const auto cell_worker_mg = 
      [&](const decltype(dof_handler.begin_mg()) &cell, 
          ScratchData<dim> &                      scratch_data, 
          CopyData &                              copy_data) { 
        this->assemble_cell(cell, scratch_data, copy_data); 
      }; 

    const auto copier_mg = [&](const CopyData &copy_data) { 
      boundary_constraints[copy_data.level].distribute_local_to_global( 
        copy_data.cell_matrix, 
        copy_data.local_dof_indices, 
        mg_matrices[copy_data.level]); 

// 如果 $(i,j)$ 是一个`interface_out` dof对，那么 $(j,i)$ 就是一个`interface_in` dof对。注意：对于 "interface_in"，我们加载接口条目的转置，即，dof对 $(j,i)$ 的条目被存储在 "interface_in(i,j)"。这是对对称情况的优化，允许在solve()中设置边缘矩阵时只使用一个矩阵。然而，在这里，由于我们的问题是非对称的，我们必须同时存储`interface_in`和`interface_out`矩阵。

      for (unsigned int i = 0; i < copy_data.dofs_per_cell; ++i) 
        for (unsigned int j = 0; j < copy_data.dofs_per_cell; ++j) 
          if (mg_constrained_dofs.is_interface_matrix_entry( 
                copy_data.level, 
                copy_data.local_dof_indices[i], 
                copy_data.local_dof_indices[j])) 
            { 
              mg_interface_out[copy_data.level].add( 
                copy_data.local_dof_indices[i], 
                copy_data.local_dof_indices[j], 
                copy_data.cell_matrix(i, j)); 
              mg_interface_in[copy_data.level].add( 
                copy_data.local_dof_indices[i], 
                copy_data.local_dof_indices[j], 
                copy_data.cell_matrix(j, i)); 
            } 
    }; 

    MeshWorker::mesh_loop(dof_handler.begin_mg(), 
                          dof_handler.end_mg(), 
                          cell_worker_mg, 
                          copier_mg, 
                          ScratchData<dim>(fe, fe.degree + 1), 
                          CopyData(), 
                          MeshWorker::assemble_own_cells); 
  } 
// @sect4{<code>AdvectionProblem::setup_smoother()</code>}  

// 接下来，我们根据`.prm`文件中的设置来设置平滑器。两个重要的选项是多网格v周期每一级的平滑前和平滑后步骤的数量以及松弛参数。

// 由于乘法往往比加法更强大，所以需要较少的平滑步骤来实现收敛，与网格大小无关。块平滑器比点平滑器也是如此。这反映在下面对每种平滑器的平滑步数的选择上。

// 点平滑器的松弛参数是在试验和错误的基础上选择的，它反映了在我们细化网格时保持GMRES求解的迭代次数不变（或尽可能接近）的必要值。在`.prm`文件中给 "Jacobi "和 "SOR "的两个值是针对1度和3度有限元的。如果用户想改成其他度数，他们可能需要调整这些数字。对于块平滑器，这个参数有一个更直接的解释，即对于二维的加法，一个DoF可以有多达4个单元的重复贡献，因此我们必须将这些方法放松0.25来补偿。对于乘法来说，这不是一个问题，因为每个单元的逆向应用都会给其所有的DoF带来新的信息。

// 最后，如上所述，点平滑器只对DoF进行操作，而块平滑器对单元进行操作，因此只有块平滑器需要被赋予有关单元排序的信息。点平滑器的DoF排序已经在`setup_system()`中得到了处理。

  template <int dim> 
  void AdvectionProblem<dim>::setup_smoother() 
  { 
    if (settings.smoother_type == "SOR") 
      { 
        using Smoother = PreconditionSOR<SparseMatrix<double>>; 

        auto smoother = 
          std::make_unique<MGSmootherPrecondition<SparseMatrix<double>, 
                                                  Smoother, 
                                                  Vector<double>>>(); 
        smoother->initialize(mg_matrices, 
                             Smoother::AdditionalData(fe.degree == 1 ? 1.0 : 
                                                                       0.62)); 
        smoother->set_steps(settings.smoothing_steps); 
        mg_smoother = std::move(smoother); 
      } 
    else if (settings.smoother_type == "Jacobi") 
      { 
        using Smoother = PreconditionJacobi<SparseMatrix<double>>; 
        auto smoother = 
          std::make_unique<MGSmootherPrecondition<SparseMatrix<double>, 
                                                  Smoother, 
                                                  Vector<double>>>(); 
        smoother->initialize(mg_matrices, 
                             Smoother::AdditionalData(fe.degree == 1 ? 0.6667 : 
                                                                       0.47)); 
        smoother->set_steps(settings.smoothing_steps); 
        mg_smoother = std::move(smoother); 
      } 
    else if (settings.smoother_type == "block SOR" || 
             settings.smoother_type == "block Jacobi") 
      { 
        smoother_data.resize(0, triangulation.n_levels() - 1); 

        for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
          { 
            DoFTools::make_cell_patches(smoother_data[level].block_list, 
                                        dof_handler, 
                                        level); 

            smoother_data[level].relaxation = 
              (settings.smoother_type == "block SOR" ? 1.0 : 0.25); 
            smoother_data[level].inversion = PreconditionBlockBase<double>::svd; 

            std::vector<unsigned int> ordered_indices; 
            switch (settings.dof_renumbering) 
              { 
                case Settings::DoFRenumberingStrategy::downstream: 
                  ordered_indices = 
                    create_downstream_cell_ordering(dof_handler, 
                                                    advection_direction, 
                                                    level); 
                  break; 

                case Settings::DoFRenumberingStrategy::upstream: 
                  ordered_indices = 
                    create_downstream_cell_ordering(dof_handler, 
                                                    -1.0 * advection_direction, 
                                                    level); 
                  break; 

                case Settings::DoFRenumberingStrategy::random: 
                  ordered_indices = 
                    create_random_cell_ordering(dof_handler, level); 
                  break; 

                case Settings::DoFRenumberingStrategy::none: 
                  break; 

                default: 
                  AssertThrow(false, ExcNotImplemented()); 
                  break; 
              } 

            smoother_data[level].order = 
              std::vector<std::vector<unsigned int>>(1, ordered_indices); 
          } 

        if (settings.smoother_type == "block SOR") 
          { 
            auto smoother = std::make_unique<MGSmootherPrecondition< 
              SparseMatrix<double>, 
              RelaxationBlockSOR<SparseMatrix<double>, double, Vector<double>>, 
              Vector<double>>>(); 
            smoother->initialize(mg_matrices, smoother_data); 
            smoother->set_steps(settings.smoothing_steps); 
            mg_smoother = std::move(smoother); 
          } 
        else if (settings.smoother_type == "block Jacobi") 
          { 
            auto smoother = std::make_unique< 
              MGSmootherPrecondition<SparseMatrix<double>, 
                                     RelaxationBlockJacobi<SparseMatrix<double>, 
                                                           double, 
                                                           Vector<double>>, 
                                     Vector<double>>>(); 
            smoother->initialize(mg_matrices, smoother_data); 
            smoother->set_steps(settings.smoothing_steps); 
            mg_smoother = std::move(smoother); 
          } 
      } 
    else 
      AssertThrow(false, ExcNotImplemented()); 
  } 
// @sect4{<code>AdvectionProblem::solve()</code>}  

// 在解决这个系统之前，我们必须首先设置多网格预处理程序。这需要设置各级之间的转换、粗略矩阵求解器和平滑器。这个设置几乎与 Step-16 相同，主要区别在于上面定义的各种平滑器，以及由于我们的问题是非对称的，我们需要不同的界面边缘矩阵。实际上，在本教程中，这些接口矩阵是空的，因为我们只使用全局细化，因此没有细化边。然而，我们在这里仍然包括了这两个矩阵，因为如果我们简单地切换到自适应细化方法，程序仍然可以正常运行）。)

// 最后要注意的是，由于我们的问题是非对称的，我们必须使用适当的Krylov子空间方法。我们在这里选择使用GMRES，因为它能保证在每次迭代中减少残差。GMRES的主要缺点是，每次迭代，存储的临时向量的数量都会增加一个，而且还需要计算与之前存储的所有向量的标量积。这是很昂贵的。通过使用重启的GMRES方法可以放松这一要求，该方法对我们在任何时候需要存储的向量数量设置了上限（这里我们在50个临时向量后重启，即48次迭代）。这样做的缺点是我们失去了在整个迭代过程中收集的信息，因此我们可以看到收敛速度较慢。因此，在哪里重启是一个平衡内存消耗、CPU工作量和收敛速度的问题。然而，本教程的目标是通过使用强大的GMG预处理程序来实现非常低的迭代次数，所以我们选择了重启长度，使下面显示的所有结果在重启发生之前就能收敛，因此我们有一个标准的GMRES方法。如果用户有兴趣，deal.II中提供的另一种合适的方法是BiCGStab。

  template <int dim> 
  void AdvectionProblem<dim>::solve() 
  { 
    const unsigned int max_iters       = 200; 
    const double       solve_tolerance = 1e-8 * system_rhs.l2_norm(); 
    SolverControl      solver_control(max_iters, solve_tolerance, true, true); 
    solver_control.enable_history_data(); 

    using Transfer = MGTransferPrebuilt<Vector<double>>; 
    Transfer mg_transfer(mg_constrained_dofs); 
    mg_transfer.build(dof_handler); 

    FullMatrix<double> coarse_matrix; 
    coarse_matrix.copy_from(mg_matrices[0]); 
    MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver; 
    coarse_grid_solver.initialize(coarse_matrix); 

    setup_smoother(); 

    mg_matrix.initialize(mg_matrices); 
    mg_interface_matrix_in.initialize(mg_interface_in); 
    mg_interface_matrix_out.initialize(mg_interface_out); 

    Multigrid<Vector<double>> mg( 
      mg_matrix, coarse_grid_solver, mg_transfer, *mg_smoother, *mg_smoother); 
    mg.set_edge_matrices(mg_interface_matrix_out, mg_interface_matrix_in); 

    PreconditionMG<dim, Vector<double>, Transfer> preconditioner(dof_handler, 
                                                                 mg, 
                                                                 mg_transfer); 

    std::cout << "     Solving with GMRES to tol " << solve_tolerance << "..." 
              << std::endl; 
    SolverGMRES<Vector<double>> solver( 
      solver_control, SolverGMRES<Vector<double>>::AdditionalData(50, true)); 

    Timer time; 
    time.start(); 
    solver.solve(system_matrix, solution, system_rhs, preconditioner); 
    time.stop(); 

    std::cout << "          converged in " << solver_control.last_step() 
              << " iterations" 
              << " in " << time.last_wall_time() << " seconds " << std::endl; 

    constraints.distribute(solution); 

    mg_smoother.release(); 
  } 
// @sect4{<code>AdvectionProblem::output_results()</code>}  

// 最后一个感兴趣的函数会生成图形输出。这里我们以.vtu格式输出解决方案和单元格排序。

// 在函数的顶部，我们为每个单元生成一个索引，以显示平滑器所使用的排序。请注意，我们只对活动单元而不是平滑器实际使用的层级做这个处理。对于点平滑器，我们对DoFs而不是单元进行重新编号，所以这只是对现实中发生的情况的一种近似。最后，这个随机排序不是我们实际使用的随机排序（见`create_smoother()`）。

// 然后，单元格的（整数）排序被复制到一个（浮点）矢量中，用于图形输出。

  template <int dim> 
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    const unsigned int n_active_cells = triangulation.n_active_cells(); 
    Vector<double>     cell_indices(n_active_cells); 
    { 
      std::vector<unsigned int> ordered_indices; 
      switch (settings.dof_renumbering) 
        { 
          case Settings::DoFRenumberingStrategy::downstream: 
            ordered_indices = 
              create_downstream_cell_ordering(dof_handler, advection_direction); 
            break; 

          case Settings::DoFRenumberingStrategy::upstream: 
            ordered_indices = 
              create_downstream_cell_ordering(dof_handler, 
                                              -1.0 * advection_direction); 
            break; 

          case Settings::DoFRenumberingStrategy::random: 
            ordered_indices = create_random_cell_ordering(dof_handler); 
            break; 

          case Settings::DoFRenumberingStrategy::none: 
            ordered_indices.resize(n_active_cells); 
            for (unsigned int i = 0; i < n_active_cells; ++i) 
              ordered_indices[i] = i; 
            break; 

          default: 
            AssertThrow(false, ExcNotImplemented()); 
            break; 
        } 

      for (unsigned int i = 0; i < n_active_cells; ++i) 
        cell_indices(ordered_indices[i]) = static_cast<double>(i); 
    } 

// 考虑到以前的教程程序，该函数的其余部分就很简单了。

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.add_data_vector(cell_indices, "cell_index"); 
    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(cycle) + ".vtu"; 
    std::ofstream output(filename.c_str()); 
    data_out.write_vtu(output); 
  } 
// @sect4{<code>AdvectionProblem::run()</code>}  

// 和大多数教程一样，这个函数创建/细化网格并调用上面定义的各种函数来设置、装配、求解和输出结果。

// 在第0个循环中，我们在正方形 <code>[-1,1]^dim</code> 上生成网格，半径为3/10个单位的孔以原点为中心。对于`manifold_id`等于1的对象（即与洞相邻的面），我们指定了一个球形流形。

  template <int dim> 
  void AdvectionProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < (settings.fe_degree == 1 ? 7 : 5); 
         ++cycle) 
      { 
        std::cout << "  Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube_with_cylindrical_hole(triangulation, 
                                                            0.3, 
                                                            1.0); 

            const SphericalManifold<dim> manifold_description(Point<dim>(0, 0)); 
            triangulation.set_manifold(1, manifold_description); 
          } 

        triangulation.refine_global(); 

        setup_system(); 

        std::cout << "     Number of active cells:       " 
                  << triangulation.n_active_cells() << " (" 
                  << triangulation.n_levels() << " levels)" << std::endl; 
        std::cout << "     Number of degrees of freedom: " 
                  << dof_handler.n_dofs() << std::endl; 

        assemble_system_and_multigrid(); 

        solve(); 

        if (settings.output) 
          output_results(cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step63 
// @sect3{The <code>main</code> function}  

// 最后，主函数和大多数教程一样。唯一有趣的一点是，我们要求用户传递一个`.prm`文件作为唯一的命令行参数。如果没有给出参数文件，程序将在屏幕上输出一个带有所有默认值的样本参数文件的内容，然后用户可以复制并粘贴到自己的`.prm`文件中。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      Step63::Settings settings; 
      settings.get_parameters((argc > 1) ? (argv[1]) : ""); 

      Step63::AdvectionProblem<2> advection_problem_2d(settings); 
      advection_problem_2d.run(); 
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


