

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2021 by the deal.II authors 
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
 * Author: Justin O'Connor, Colorado State University, 2021. 
 */ 


// @sect3{Preliminaries}  
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/tensor.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/signaling_nan.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/linear_operator.h> 
#include <deal.II/lac/packaged_operation.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <iostream> 
#include <fstream> 
#include <algorithm> 

// 以上是相当常见的包含文件。这些文件还包括稀疏直接类的文件 SparseDirectUMFPACK。这不是解决大型线性问题的最有效的方法，但现在可以了。

// 像往常一样，我们把所有的东西都放到一个共同的命名空间里。然后，我们开始声明一些常数的符号名称，这些常数将在本教程中使用。具体来说，我们在这个程序中有*多的变量（当然是密度和位移，但也有未过滤的密度和相当多的拉格朗日乘数）。我们很容易忘记这些变量在求解向量中的哪个位置，而且试图用数字来表示这些向量分量是一个错误的处方。相反，我们定义的静态变量可以在所有这些地方使用，而且只需初始化一次。在实践中，这将导致一些冗长的表达式，但它们更具可读性，而且不太可能出错。

// 一个类似的问题出现在系统矩阵和向量中块的排序上。矩阵中有 $9\times 9$ 块，而且很难记住哪个是哪个。对这些块也使用符号名称要容易得多。

// 最后，我们为我们将要使用的边界指标引入符号名称，与  step-19  中的精神相同。

// 在所有这些情况下，我们将这些变量声明为命名空间中的成员。在求解组件的情况下，这些变量的具体数值取决于空间维度，因此我们使用[模板变量](https:en.cppreference.com/w/cpp/language/variable_template)来使变量的数值取决于模板参数，就像我们经常使用模板函数一样。

namespace SAND 
{ 
  using namespace dealii; 

// 这个命名空间记录了我们的有限元系统中与每个变量相对应的第一个组件。

  namespace SolutionComponents 
  { 
    template <int dim> 
    constexpr unsigned int density = 0; 
    template <int dim> 
    constexpr unsigned int displacement = 1; 
    template <int dim> 
    constexpr unsigned int unfiltered_density = 1 + dim; 
    template <int dim> 
    constexpr unsigned int displacement_multiplier = 2 + dim; 
    template <int dim> 
    constexpr unsigned int unfiltered_density_multiplier = 2 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_lower_slack = 3 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_lower_slack_multiplier = 4 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_upper_slack = 5 + 2 * dim; 
    template <int dim> 
    constexpr unsigned int density_upper_slack_multiplier = 6 + 2 * dim; 
  } // namespace SolutionComponents 

// 这是一个命名空间，它记录了哪个区块对应于哪个变量。

  namespace SolutionBlocks 
  { 
    constexpr unsigned int density                        = 0; 
    constexpr unsigned int displacement                   = 1; 
    constexpr unsigned int unfiltered_density             = 2; 
    constexpr unsigned int displacement_multiplier        = 3; 
    constexpr unsigned int unfiltered_density_multiplier  = 4; 
    constexpr unsigned int density_lower_slack            = 5; 
    constexpr unsigned int density_lower_slack_multiplier = 6; 
    constexpr unsigned int density_upper_slack            = 7; 
    constexpr unsigned int density_upper_slack_multiplier = 8; 
  } // namespace SolutionBlocks 

  namespace BoundaryIds 
  { 
    constexpr types::boundary_id down_force = 101; 
    constexpr types::boundary_id no_force   = 102; 
  } // namespace BoundaryIds 

  namespace ValueExtractors 
  { 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      densities(SolutionComponents::density<dim>); 
    template <int dim> 
    const FEValuesExtractors::Vector 
      displacements(SolutionComponents::displacement<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      unfiltered_densities(SolutionComponents::unfiltered_density<dim>); 
    template <int dim> 
    const FEValuesExtractors::Vector displacement_multipliers( 
      SolutionComponents::displacement_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar unfiltered_density_multipliers( 
      SolutionComponents::unfiltered_density_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      density_lower_slacks(SolutionComponents::density_lower_slack<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar density_lower_slack_multipliers( 
      SolutionComponents::density_lower_slack_multiplier<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar 
      density_upper_slacks(SolutionComponents::density_upper_slack<dim>); 
    template <int dim> 
    const FEValuesExtractors::Scalar density_upper_slack_multipliers( 
      SolutionComponents::density_upper_slack_multiplier<dim>); 
  } // namespace ValueExtractors 
// @sect3{The SANDTopOpt main class}  

// 接下来是这个问题的主类。大多数函数都遵循教程程序的常规命名方式，不过有几个函数因为长度问题被从通常称为`setup_system()`的函数中分离出来，还有一些函数是处理优化算法的各个方面的。

// 作为额外的奖励，该程序将计算出的设计写成STL文件，例如，可以将其发送给3D打印机。

  template <int dim> 
  class SANDTopOpt 
  { 
  public: 
    SANDTopOpt(); 

    void run(); 

  private: 
    void create_triangulation(); 

    void setup_boundary_values(); 

    void setup_block_system(); 

    void setup_filter_matrix(); 

    void assemble_system(); 

    BlockVector<double> solve(); 

    std::pair<double, double> 
    calculate_max_step_size(const BlockVector<double> &state, 
                            const BlockVector<double> &step) const; 

    BlockVector<double> 
    calculate_test_rhs(const BlockVector<double> &test_solution) const; 

    double calculate_exact_merit(const BlockVector<double> &test_solution); 

    BlockVector<double> find_max_step(); 

    BlockVector<double> compute_scaled_step(const BlockVector<double> &state, 
                                            const BlockVector<double> &step, 
                                            const double descent_requirement); 

    bool check_convergence(const BlockVector<double> &state); 

    void output_results(const unsigned int j) const; 

    void write_as_stl(); 

    std::set<typename Triangulation<dim>::cell_iterator> 
    find_relevant_neighbors( 
      typename Triangulation<dim>::cell_iterator cell) const; 

// 大部分的成员变量也是标准的。但是，有一些变量是专门与优化算法有关的（比如下面的各种标量因子），以及过滤器矩阵，以确保设计保持平稳。

    Triangulation<dim>        triangulation; 
    FESystem<dim>             fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 

    std::map<types::global_dof_index, double> boundary_values; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    SparsityPattern      filter_sparsity_pattern; 
    SparseMatrix<double> filter_matrix; 

    BlockVector<double> system_rhs; 
    BlockVector<double> nonlinear_solution; 

    const double density_ratio; 
    const double density_penalty_exponent; 
    const double filter_r; 
    double       penalty_multiplier; 
    double       barrier_size; 

    TimerOutput timer; 
  }; 
// @sect3{Constructor and set-up functions}  

// 我们初始化一个由2  $\times$  dim `FE_Q(1)`元素组成的FES系统，用于位移变量及其拉格朗日乘数，以及7 `FE_DGQ(0)`元素。 这些片状常数函数用于与密度相关的变量：密度本身、未过滤的密度、用于未过滤的密度的下限和上限的松弛变量，然后是用于过滤和未过滤的密度之间的连接以及不等式约束的拉格朗日乘子。

// 这些元素出现的顺序在上面有记载。

  template <int dim> 
  SANDTopOpt<dim>::SANDTopOpt() 
    : fe(FE_DGQ<dim>(0), 
         1, 
         (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
         1, 
         FE_DGQ<dim>(0), 
         1, 
         (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
         1, 
         FE_DGQ<dim>(0), 
         5) 
    , dof_handler(triangulation) 
    , density_ratio(.5) 
    , density_penalty_exponent(3) 
    , filter_r(.251) 
    , penalty_multiplier(1) 
    , timer(std::cout, TimerOutput::summary, TimerOutput::wall_times) 
  { 
    Assert(dim > 1, ExcNotImplemented()); 
  } 

// 然后，第一步是创建与介绍中的问题描述相匹配的三角形--一个6乘1的矩形（或者一个6乘1乘1的3D盒子），在这个盒子的顶部中心将施加一个力。然后，这个三角形被均匀地细化若干次。

// 与本程序的其他部分相比，这个函数特别假定我们是在2D中，如果我们想转到3D模拟，就需要进行修改。我们通过函数顶部的断言来确保没有人试图不经修改就意外地在三维中运行。

  template <int dim> 
  void SANDTopOpt<dim>::create_triangulation() 
  { 
    Assert(dim == 2, ExcNotImplemented()); 
    GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                              {6, 1}, 
                                              Point<dim>(0, 0), 
                                              Point<dim>(6, 1)); 

    triangulation.refine_global(3); 

// 第二步是将边界指标应用于边界的一部分。下面的代码分别为盒子的底部、顶部、左侧和右侧的边界分配了边界指示器。顶部边界的中心区域被赋予一个单独的边界指示器。这就是我们要施加向下力的地方。

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 
                if (std::fabs(center(1) - 1) < 1e-12) 
                  { 
                    if ((std::fabs(center(0) - 3) < .3)) 
                      face->set_boundary_id(BoundaryIds::down_force); 
                    else 
                      face->set_boundary_id(BoundaryIds::no_force); 
                  } 
                else 
                  face->set_boundary_id(BoundaryIds::no_force); 
              } 
          } 
      } 
  } 

// 接下来，确定由于边界值而产生的约束。 域的底角在 $y$ 方向保持不变--左下角也在 $x$ 方向。deal.II通常认为边界值是附着在边界的片段上的，即面，而不是单个顶点。的确，从数学上讲，对于无穷大的偏微分方程，我们不能把边界值分配给单个点。但是，由于我们试图重现一个广泛使用的基准，我们还是要这样做，并牢记我们有一个有限维的问题，在单个节点上施加边界条件是有效的。

  template <int dim> 
  void SANDTopOpt<dim>::setup_boundary_values() 
  { 
    boundary_values.clear(); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 

// 检查当前面是否在底层边界上，如果是，则检查其顶点之一是否可能是左底层或右底层顶点。

                if (std::fabs(center(1) - 0) < 1e-12) 
                  { 
                    for (const auto vertex_number : cell->vertex_indices()) 
                      { 
                        const auto vert = cell->vertex(vertex_number); 

                        if (std::fabs(vert(0) - 0) < 1e-12 && 
                            std::fabs(vert(1) - 0) < 1e-12) 
                          { 
                            types::global_dof_index x_displacement = 
                              cell->vertex_dof_index(vertex_number, 0); 
                            types::global_dof_index y_displacement = 
                              cell->vertex_dof_index(vertex_number, 1); 
                            types::global_dof_index x_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 2); 
                            types::global_dof_index y_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 3); 

                            boundary_values[x_displacement]            = 0; 
                            boundary_values[y_displacement]            = 0; 
                            boundary_values[x_displacement_multiplier] = 0; 
                            boundary_values[y_displacement_multiplier] = 0; 
                          } 

                        else if (std::fabs(vert(0) - 6) < 1e-12 && 
                                 std::fabs(vert(1) - 0) < 1e-12) 
                          { 
                            types::global_dof_index y_displacement = 
                              cell->vertex_dof_index(vertex_number, 1); 
                            types::global_dof_index y_displacement_multiplier = 
                              cell->vertex_dof_index(vertex_number, 3); 

                            boundary_values[y_displacement]            = 0; 
                            boundary_values[y_displacement_multiplier] = 0; 
                          } 
                      } 
                  } 
              } 
          } 
      } 
  } 
// @sect3{Setting up block matrices and vectors}  

// 下一个函数制作了一个巨大的9乘9的块状矩阵，并且还设置了必要的块状向量。 这个矩阵的稀疏度模式包括滤波矩阵的稀疏度模式。它还初始化了我们将使用的任何块向量。

// 设置块本身并不复杂，并且遵循诸如  step-22  等程序中已经完成的工作，例如。

  template <int dim> 
  void SANDTopOpt<dim>::setup_block_system() 
  { 
    std::vector<unsigned int> block_component(9, 2); 
    block_component[0] = 0; 
    block_component[1] = 1; 
    const std::vector<types::global_dof_index> dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 

    const types::global_dof_index                     n_p = dofs_per_block[0]; 
    const types::global_dof_index                     n_u = dofs_per_block[1]; 
    const std::vector<BlockVector<double>::size_type> block_sizes = { 
      n_p, n_u, n_p, n_u, n_p, n_p, n_p, n_p, n_p}; 

    BlockDynamicSparsityPattern dsp(9, 9); 
    for (unsigned int k = 0; k < 9; ++k) 
      for (unsigned int j = 0; j < 9; ++j) 
        dsp.block(j, k).reinit(block_sizes[j], block_sizes[k]); 
    dsp.collect_sizes(); 

// 该函数的大部分内容是设置这些块中哪些将实际包含任何内容，即哪些变量与哪些其他变量相耦合。这很麻烦，但也是必要的，以确保我们不会为我们的矩阵分配大量的条目，而这些条目最终会变成零。

// 你在下面看到的具体模式可能需要在纸上画一次，但是从我们在每次非线性迭代中必须组装的双线性形式的许多项来看，它是相对直接的方式。

// 使用命名空间 "SolutionComponents "中定义的符号名称有助于理解下面每个项所对应的内容，但它也使表达式变得冗长而不流畅。像 `coupling[SolutionComponents::density_upper_slack_multiplier<dim>][SolutionComponents::density<dim>]` 这样的术语读起来就不太顺口，要么必须分成几行，要么几乎跑到每个屏幕的右边缘。因此，我们打开了一个大括号封闭的代码块，在这个代码块中，我们通过说 "使用命名空间SolutionComponents"，暂时使命名空间`SolutionComponents'中的名字可用，而不需要命名空间修饰语。

    Table<2, DoFTools::Coupling> coupling(2 * dim + 7, 2 * dim + 7); 
    { 
      using namespace SolutionComponents; 

      coupling[density<dim>][density<dim>] = DoFTools::always; 

      for (unsigned int i = 0; i < dim; ++i) 
        { 
          coupling[density<dim>][displacement<dim> + i] = DoFTools::always; 
          coupling[displacement<dim> + i][density<dim>] = DoFTools::always; 
        } 

      for (unsigned int i = 0; i < dim; ++i) 
        { 
          coupling[density<dim>][displacement_multiplier<dim> + i] = 
            DoFTools::always; 
          coupling[displacement_multiplier<dim> + i][density<dim>] = 
            DoFTools::always; 
        } 

      coupling[density<dim>][unfiltered_density_multiplier<dim>] = 
        DoFTools::always; 
      coupling[unfiltered_density_multiplier<dim>][density<dim>] = 
        DoFTools::always; 
      /*位移的联结  */ 
      for (unsigned int i = 0; i < dim; ++i) 
        { 
          for (unsigned int k = 0; k < dim; ++k) 
            { 
              coupling[displacement<dim> + i] 
                      [displacement_multiplier<dim> + k] = DoFTools::always; 
              coupling[displacement_multiplier<dim> + k] 
                      [displacement<dim> + i] = DoFTools::always; 
            } 
        } 
      /*松弛变量的耦合 */ 
      coupling[density_lower_slack<dim>][density_lower_slack<dim>] = 
        DoFTools::always; 
      coupling[density_lower_slack<dim>][density_upper_slack<dim>] = 
        DoFTools::always; 
      coupling[density_upper_slack<dim>][density_lower_slack<dim>] = 
        DoFTools::always; 

      coupling[density_lower_slack_multiplier<dim>] 
              [density_lower_slack_multiplier<dim>] = DoFTools::always; 
      coupling[density_lower_slack_multiplier<dim>] 
              [density_upper_slack_multiplier<dim>] = DoFTools::always; 
      coupling[density_upper_slack_multiplier<dim>] 
              [density_lower_slack_multiplier<dim>] = DoFTools::always; 
    } 

// 在创建稀疏模式之前，我们还必须设置约束。由于这个程序没有自适应地细化网格，我们唯一的约束是将所有的密度变量耦合在一起，强制执行体积约束。这将最终导致矩阵的密集子块，但我们对此没有什么办法。

    const ComponentMask density_mask = 
      fe.component_mask(ValueExtractors::densities<dim>); 
    const IndexSet density_dofs = 
      DoFTools::extract_dofs(dof_handler, density_mask); 

    types::global_dof_index last_density_dof = 
      density_dofs.nth_index_in_set(density_dofs.n_elements() - 1); 
    constraints.clear(); 
    constraints.add_line(last_density_dof); 
    for (unsigned int i = 0; i < density_dofs.n_elements() - 1; ++i) 
      constraints.add_entry(last_density_dof, 
                            density_dofs.nth_index_in_set(i), 
                            -1); 
    constraints.set_inhomogeneity(last_density_dof, 0); 

    constraints.close(); 

// 现在我们终于可以为矩阵创建稀疏模式了，考虑到哪些变量与哪些其他变量耦合，以及我们对密度的约束。

    DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp, constraints); 

// 矩阵中唯一没有处理的部分是过滤矩阵和它的转置。这些都是非局部（积分）运算符，目前deal.II还没有相关的函数。我们最终需要做的是遍历所有单元，并将此单元上的未过滤密度与小于阈值距离的相邻单元的所有过滤密度联系起来，反之亦然；目前，我们只关心建立与这种矩阵相对应的稀疏模式，所以我们执行等效循环，以后我们将写进矩阵的一个条目，现在我们只需向稀疏矩阵添加一个条目。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (const auto &check_cell : find_relevant_neighbors(cell)) 
          { 
            const double distance = 
              cell->center().distance(check_cell->center()); 
            if (distance < filter_r) 
              { 
                dsp 
                  .block(SolutionBlocks::unfiltered_density, 
                         SolutionBlocks::unfiltered_density_multiplier) 
                  .add(i, check_cell->active_cell_index()); 
                dsp 
                  .block(SolutionBlocks::unfiltered_density_multiplier, 
                         SolutionBlocks::unfiltered_density) 
                  .add(i, check_cell->active_cell_index()); 
              } 
          } 
      } 

// 在生成了 "动态 "稀疏度模式之后，我们终于可以将其复制到用于将矩阵与稀疏度模式联系起来的结构中。由于稀疏模式很大很复杂，我们还将其输出到一个自己的文件中，以达到可视化的目的--换句话说，是为了 "可视化调试"。

    sparsity_pattern.copy_from(dsp); 

    std::ofstream out("sparsity.plt"); 
    sparsity_pattern.print_gnuplot(out); 

    system_matrix.reinit(sparsity_pattern); 

// 剩下的就是正确确定各种向量及其块的大小，以及为（非线性）解向量的一些分量设置初始猜测。我们在这里使用解向量各个区块的符号分量名称，为了简洁起见，使用与上面的 "使用命名空间 "相同的技巧。

    nonlinear_solution.reinit(block_sizes); 
    system_rhs.reinit(block_sizes); 

    { 
      using namespace SolutionBlocks; 
      nonlinear_solution.block(density).add(density_ratio); 
      nonlinear_solution.block(unfiltered_density).add(density_ratio); 
      nonlinear_solution.block(unfiltered_density_multiplier) 
        .add(density_ratio); 
      nonlinear_solution.block(density_lower_slack).add(density_ratio); 
      nonlinear_solution.block(density_lower_slack_multiplier).add(50); 
      nonlinear_solution.block(density_upper_slack).add(1 - density_ratio); 
      nonlinear_solution.block(density_upper_slack_multiplier).add(50); 
    } 
  } 
// @sect3{Creating the filter matrix}  

// 接下来是一个在程序开始时使用一次的函数。它创建了一个矩阵 $H$ ，使过滤后的密度向量等于 $H$ 乘以未过滤的密度。 这个矩阵的创建是非同小可的，它在每次迭代中都会被使用，因此，与其像我们对牛顿矩阵那样对其进行改造，不如只做一次并单独存储。

// 这个矩阵的计算方式遵循上面已经使用过的大纲，以形成其稀疏模式。我们在这里对这个单独形成的矩阵的稀疏性模式重复这个过程，然后实际建立矩阵本身。你可能想看看本程序介绍中关于这个矩阵的定义。

  template <int dim> 
  void SANDTopOpt<dim>::setup_filter_matrix() 
  { 

// 滤波器的稀疏模式已经在setup_system()函数中确定并实现。我们从相应的块中复制该结构，并在这里再次使用它。

    filter_sparsity_pattern.copy_from( 
      sparsity_pattern.block(SolutionBlocks::unfiltered_density, 
                             SolutionBlocks::unfiltered_density_multiplier)); 
    filter_matrix.reinit(filter_sparsity_pattern); 

// 在建立了稀疏模式之后，现在我们重新做所有这些循环，以实际计算矩阵项的必要值。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (const auto &check_cell : find_relevant_neighbors(cell)) 
          { 
            const double distance = 
              cell->center().distance(check_cell->center()); 
            if (distance < filter_r) 
              { 
                filter_matrix.add(i, 
                                  check_cell->active_cell_index(), 
                                  filter_r - distance); 

//      

              } 
          } 
      } 

// 最后一步是对矩阵进行标准化处理，使每一行的条目之和等于1。

    for (unsigned int i = 0; i < filter_matrix.m(); ++i) 
      { 
        double denominator = 0; 
        for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             iter++) 
          denominator = denominator + iter->value(); 
        for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             iter++) 
          iter->value() = iter->value() / denominator; 
      } 
  } 

// 这个函数用于建立过滤矩阵。我们创建一个输入单元的一定半径内的所有单元迭代器的集合。这些是与过滤器有关的邻近单元。

  template <int dim> 
  std::set<typename Triangulation<dim>::cell_iterator> 
  SANDTopOpt<dim>::find_relevant_neighbors( 
    typename Triangulation<dim>::cell_iterator cell) const 
  { 
    std::set<unsigned int>                               neighbor_ids; 
    std::set<typename Triangulation<dim>::cell_iterator> cells_to_check; 

    neighbor_ids.insert(cell->active_cell_index()); 
    cells_to_check.insert(cell); 

    bool new_neighbors_found; 
    do 
      { 
        new_neighbors_found = false; 
        for (const auto &check_cell : 
             std::vector<typename Triangulation<dim>::cell_iterator>( 
               cells_to_check.begin(), cells_to_check.end())) 
          { 
            for (const auto n : check_cell->face_indices()) 
              { 
                if (!(check_cell->face(n)->at_boundary())) 
                  { 
                    const auto & neighbor = check_cell->neighbor(n); 
                    const double distance = 
                      cell->center().distance(neighbor->center()); 
                    if ((distance < filter_r) && 
                        !(neighbor_ids.count(neighbor->active_cell_index()))) 
                      { 
                        cells_to_check.insert(neighbor); 
                        neighbor_ids.insert(neighbor->active_cell_index()); 
                        new_neighbors_found = true; 
                      } 
                  } 
              } 
          } 
      } 
    while (new_neighbors_found); 
    return cells_to_check; 
  } 
// @sect3{Assembling the Newton matrix}  

// setup_filter_matrix函数建立了一个只要网格不改变就不变的矩阵（在这个程序中我们反正不改变），而下一个函数建立了每次迭代都要解决的矩阵。这就是奇迹发生的地方。描述牛顿求解KKT条件的方法的线性方程组的组成部分在这里实现。

// 这个函数的顶部与大多数此类函数一样，只是设置了实际装配所需的各种变量，包括一大堆提取器。如果你以前看过  step-22  ，整个设置应该看起来很熟悉，尽管有些冗长。

  template <int dim> 
  void SANDTopOpt<dim>::assemble_system() 
  { 
    TimerOutput::Scope t(timer, "assembly"); 

    system_matrix = 0; 
    system_rhs    = 0; 

    MappingQGeneric<dim> mapping(1); 
    QGauss<dim>          quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1>      face_quadrature_formula(fe.degree + 1); 
    FEValues<dim>        fe_values(mapping, 
                            fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim>    fe_face_values(mapping, 
                                     fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_normal_vectors | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     dummy_cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<double>                    lambda_values(n_q_points); 
    std::vector<double>                    mu_values(n_q_points); 
    const Functions::ConstantFunction<dim> lambda(1.); 
    const Functions::ConstantFunction<dim> mu(1.); 
    std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 

// 在这一点上，我们对未过滤的密度进行过滤，并对未过滤的密度乘法器进行邻接（转置）操作，都是对当前非线性解决方案的最佳猜测。后来我们用它来告诉我们，我们过滤的密度与应用于未过滤密度的过滤器有多大的偏差。这是因为在非线性问题的解中，我们有 $\rho=H\varrho$ ，但在中间迭代中，我们一般有 $\rho^k\neq H\varrho^k$ ，然后 "残差" $\rho^k-H\varrho^k$ 将出现在我们下面计算的牛顿更新方程中的右边。

    BlockVector<double> filtered_unfiltered_density_solution = 
      nonlinear_solution; 
    BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
      nonlinear_solution; 

    filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
                          SolutionBlocks::unfiltered_density), 
                        nonlinear_solution.block( 
                          SolutionBlocks::unfiltered_density)); 
    filter_matrix.Tvmult( 
      filter_adjoint_unfiltered_density_multiplier_solution.block( 
        SolutionBlocks::unfiltered_density_multiplier), 
      nonlinear_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 

    std::vector<double>                  old_density_values(n_q_points); 
    std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
    std::vector<double>                  old_displacement_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
    std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
    std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
      n_q_points); 
    std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
    std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
    std::vector<double> old_lower_slack_values(n_q_points); 
    std::vector<double> old_upper_slack_values(n_q_points); 
    std::vector<double> old_unfiltered_density_values(n_q_points); 
    std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
    std::vector<double> filtered_unfiltered_density_values(n_q_points); 
    std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
      n_q_points); 

    using namespace ValueExtractors; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 

        cell->get_dof_indices(local_dof_indices); 

        fe_values.reinit(cell); 

        lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
        mu.value_list(fe_values.get_quadrature_points(), mu_values); 

// 作为构建系统矩阵的一部分，我们需要从我们目前对解决方案的猜测中获取数值。以下几行代码将检索出所需的值。

        fe_values[densities<dim>].get_function_values(nonlinear_solution, 
                                                      old_density_values); 
        fe_values[displacements<dim>].get_function_values( 
          nonlinear_solution, old_displacement_values); 
        fe_values[displacements<dim>].get_function_divergences( 
          nonlinear_solution, old_displacement_divs); 
        fe_values[displacements<dim>].get_function_symmetric_gradients( 
          nonlinear_solution, old_displacement_symmgrads); 
        fe_values[displacement_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_displacement_multiplier_values); 
        fe_values[displacement_multipliers<dim>].get_function_divergences( 
          nonlinear_solution, old_displacement_multiplier_divs); 
        fe_values[displacement_multipliers<dim>] 
          .get_function_symmetric_gradients( 
            nonlinear_solution, old_displacement_multiplier_symmgrads); 
        fe_values[density_lower_slacks<dim>].get_function_values( 
          nonlinear_solution, old_lower_slack_values); 
        fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_lower_slack_multiplier_values); 
        fe_values[density_upper_slacks<dim>].get_function_values( 
          nonlinear_solution, old_upper_slack_values); 
        fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_upper_slack_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          nonlinear_solution, old_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          nonlinear_solution, old_unfiltered_density_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          filtered_unfiltered_density_solution, 
          filtered_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          filter_adjoint_unfiltered_density_multiplier_solution, 
          filter_adjoint_unfiltered_density_multiplier_values); 

        for (const auto q_point : fe_values.quadrature_point_indices()) 
          { 

// 我们还需要几个与来自拉格朗日的第一导数的测试函数相对应的数值，也就是 $d_{\bullet}$ 函数。这些都是在这里计算的。

            for (const auto i : fe_values.dof_indices()) 
              { 
                const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
                  fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
                const double displacement_phi_i_div = 
                  fe_values[displacements<dim>].divergence(i, q_point); 

                const SymmetricTensor<2, dim> 
                  displacement_multiplier_phi_i_symmgrad = 
                    fe_values[displacement_multipliers<dim>].symmetric_gradient( 
                      i, q_point); 
                const double displacement_multiplier_phi_i_div = 
                  fe_values[displacement_multipliers<dim>].divergence(i, 
                                                                      q_point); 

                const double density_phi_i = 
                  fe_values[densities<dim>].value(i, q_point); 
                const double unfiltered_density_phi_i = 
                  fe_values[unfiltered_densities<dim>].value(i, q_point); 
                const double unfiltered_density_multiplier_phi_i = 
                  fe_values[unfiltered_density_multipliers<dim>].value(i, 
                                                                       q_point); 

                const double lower_slack_multiplier_phi_i = 
                  fe_values[density_lower_slack_multipliers<dim>].value( 
                    i, q_point); 

                const double lower_slack_phi_i = 
                  fe_values[density_lower_slacks<dim>].value(i, q_point); 

                const double upper_slack_phi_i = 
                  fe_values[density_upper_slacks<dim>].value(i, q_point); 

                const double upper_slack_multiplier_phi_i = 
                  fe_values[density_upper_slack_multipliers<dim>].value( 
                    i, q_point); 

                for (const auto j : fe_values.dof_indices()) 
                  { 

// 最后，我们需要来自拉格朗日的第二轮导数的数值，即 $c_{\bullet}$ 函数。这些是在这里计算的。

                    const SymmetricTensor<2, dim> displacement_phi_j_symmgrad = 
                      fe_values[displacements<dim>].symmetric_gradient(j, 
                                                                       q_point); 
                    const double displacement_phi_j_div = 
                      fe_values[displacements<dim>].divergence(j, q_point); 

                    const SymmetricTensor<2, dim> 
                      displacement_multiplier_phi_j_symmgrad = 
                        fe_values[displacement_multipliers<dim>] 
                          .symmetric_gradient(j, q_point); 
                    const double displacement_multiplier_phi_j_div = 
                      fe_values[displacement_multipliers<dim>].divergence( 
                        j, q_point); 

                    const double density_phi_j = 
                      fe_values[densities<dim>].value(j, q_point); 

                    const double unfiltered_density_phi_j = 
                      fe_values[unfiltered_densities<dim>].value(j, q_point); 
                    const double unfiltered_density_multiplier_phi_j = 
                      fe_values[unfiltered_density_multipliers<dim>].value( 
                        j, q_point); 

                    const double lower_slack_phi_j = 
                      fe_values[density_lower_slacks<dim>].value(j, q_point); 

                    const double upper_slack_phi_j = 
                      fe_values[density_upper_slacks<dim>].value(j, q_point); 

                    const double lower_slack_multiplier_phi_j = 
                      fe_values[density_lower_slack_multipliers<dim>].value( 
                        j, q_point); 

                    const double upper_slack_multiplier_phi_j = 
                      fe_values[density_upper_slack_multipliers<dim>].value( 
                        j, q_point); 

// 这就是实际工作的开始。在下文中，我们将建立矩阵的所有项--它们数量众多，而且不完全是不言自明的，也取决于之前的解和它的导数（我们已经在上面评估了这些导数，并将其放入名为`old_*`的变量中）。为了理解这些条款的每一个对应的内容，你要看一下上面介绍中这些条款的明确形式。                    被驱动到0的方程的右边给出了寻找局部最小值的所有KKT条件--每个单独方程的描述都是随着右边的计算给出的。

                    /* 方程1  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      ( 
                        -density_phi_i * unfiltered_density_multiplier_phi_j 
                        + density_penalty_exponent * 
                            (density_penalty_exponent - 1) * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 2) * 
                            density_phi_i * density_phi_j * 
                            (old_displacement_multiplier_divs[q_point] * 
                               old_displacement_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_symmgrads[q_point] * 
                                old_displacement_multiplier_symmgrads[q_point])) 
                        + density_penalty_exponent * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 1) * 
                            density_phi_i * 
                            (displacement_multiplier_phi_j_div * 
                               old_displacement_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_symmgrads[q_point] * 
                                displacement_multiplier_phi_j_symmgrad)) 
                        + density_penalty_exponent * 
                            std::pow(old_density_values[q_point], 
                                     density_penalty_exponent - 1) * 
                            density_phi_i * 
                            (displacement_phi_j_div * 
                               old_displacement_multiplier_divs[q_point] * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (old_displacement_multiplier_symmgrads[q_point] * 
                                displacement_phi_j_symmgrad))); 
                   
                    /* 方程2  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (density_penalty_exponent * 
                         std::pow(old_density_values[q_point], 
                                  density_penalty_exponent - 1) * 
                         density_phi_j * 
                         (old_displacement_multiplier_divs[q_point] * 
                            displacement_phi_i_div * lambda_values[q_point] + 
                          2 * mu_values[q_point] * 
                            (old_displacement_multiplier_symmgrads[q_point] * 
                             displacement_phi_i_symmgrad)) 
                       + std::pow(old_density_values[q_point], 
                                  density_penalty_exponent) * 
                           (displacement_multiplier_phi_j_div * 
                              displacement_phi_i_div * lambda_values[q_point] + 
                            2 * mu_values[q_point] * 
                              (displacement_multiplier_phi_j_symmgrad * 
                               displacement_phi_i_symmgrad)) 
                      ); 

                   /*方程3，这与过滤器有关 */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (-1 * unfiltered_density_phi_i * 
                         lower_slack_multiplier_phi_j + 
                       unfiltered_density_phi_i * upper_slack_multiplier_phi_j); 

                     /* 方程4：原始可行性  */ 
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      ( 
                        density_penalty_exponent * 
                          std::pow(old_density_values[q_point], 
                                   density_penalty_exponent - 1) * 
                          density_phi_j * 
                          (old_displacement_divs[q_point] * 
                             displacement_multiplier_phi_i_div * 
                             lambda_values[q_point] + 
                           2 * mu_values[q_point] * 
                             (old_displacement_symmgrads[q_point] * 
                              displacement_multiplier_phi_i_symmgrad)) 

                        + std::pow(old_density_values[q_point], 
                                   density_penalty_exponent) * 
                            (displacement_phi_j_div * 
                               displacement_multiplier_phi_i_div * 
                               lambda_values[q_point] + 
                             2 * mu_values[q_point] * 
                               (displacement_phi_j_symmgrad * 
                                displacement_multiplier_phi_i_symmgrad))); 

                   /*等式5：原始可行性  */ 
                    cell_matrix(i, j) += 
                      -1 * fe_values.JxW(q_point) * 
                      lower_slack_multiplier_phi_i * 
                      (unfiltered_density_phi_j - lower_slack_phi_j); 
                  /* 等式6：原始可行性  */ 
                    cell_matrix(i, j) += 
                      -1 * fe_values.JxW(q_point) * 
                      upper_slack_multiplier_phi_i * 
                      (-1 * unfiltered_density_phi_j - upper_slack_phi_j); 
                    /* Equation 7: Primal feasibility - the part with the filter
                     * is added later */
                    cell_matrix(i, j) += -1 * fe_values.JxW(q_point) * 
                                         unfiltered_density_multiplier_phi_i * 
                                         (density_phi_j); 
                    /* Equation 8: Complementary slackness */
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (lower_slack_phi_i * lower_slack_multiplier_phi_j 

                       + lower_slack_phi_i * lower_slack_phi_j * 
                           old_lower_slack_multiplier_values[q_point] / 
                           old_lower_slack_values[q_point]); 
                    /* Equation 9: Complementary slackness */
                    cell_matrix(i, j) += 
                      fe_values.JxW(q_point) * 
                      (upper_slack_phi_i * upper_slack_multiplier_phi_j 

                       + upper_slack_phi_i * upper_slack_phi_j * 
                           old_upper_slack_multiplier_values[q_point] / 
                           old_upper_slack_values[q_point]); 
                  } 
              } 
          } 

// 现在我们已经把所有的东西都组装好了，我们要做的就是处理（Dirichlet）边界条件的影响和其他约束。我们将前者与当前单元的贡献结合在一起，然后让AffineConstraint类来处理后者，同时将当前单元的贡献复制到全局线性系统中。

        MatrixTools::local_apply_boundary_values(boundary_values, 
                                                 local_dof_indices, 
                                                 cell_matrix, 
                                                 dummy_cell_rhs, 
                                                 true); 

        constraints.distribute_local_to_global(cell_matrix, 
                                               local_dof_indices, 
                                               system_matrix); 
      } 

// 在积累了所有属于牛顿矩阵的项之后，我们现在还必须计算右手边的项（即负残差）。我们已经在另一个函数中做了这个工作，所以我们在这里调用它。

    system_rhs = calculate_test_rhs(nonlinear_solution); 

// 这里我们使用我们已经构建好的过滤器矩阵。我们只需要整合这个应用于测试函数的过滤器，它是片状常数，所以整合变成了简单的乘以单元格的度量。 遍历预制的过滤器矩阵可以让我们使用哪些单元格在过滤器中或不在过滤器中的信息，而不需要再次重复检查邻居单元格。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int i = cell->active_cell_index(); 
        for (typename SparseMatrix<double>::iterator iter = 
               filter_matrix.begin(i); 
             iter != filter_matrix.end(i); 
             ++iter) 
          { 
            const unsigned int j     = iter->column(); 
            const double       value = iter->value() * cell->measure(); 

            system_matrix 
              .block(SolutionBlocks::unfiltered_density_multiplier, 
                     SolutionBlocks::unfiltered_density) 
              .add(i, j, value); 
            system_matrix 
              .block(SolutionBlocks::unfiltered_density, 
                     SolutionBlocks::unfiltered_density_multiplier) 
              .add(j, i, value); 
          } 
      } 
  } 
// @sect3{Solving the Newton linear system}  

// 我们将需要在每次迭代中解决一个线性系统。我们暂时使用一个直接求解器--对于一个有这么多非零值的矩阵来说，这显然不是一个有效的选择，而且它不会扩展到任何有趣的地方。对于 "真正的 "应用，我们将需要一个迭代求解器，但系统的复杂性意味着一个迭代求解器的算法将需要大量的工作。因为这不是当前程序的重点，所以我们简单地坚持使用我们在这里的直接求解器--该函数遵循与 step-29 中使用的相同结构。

  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::solve() 
  { 
    TimerOutput::Scope t(timer, "solver"); 

    BlockVector<double> linear_solution; 
    linear_solution.reinit(nonlinear_solution); 

    SparseDirectUMFPACK A_direct; 
    A_direct.initialize(system_matrix); 
    A_direct.vmult(linear_solution, system_rhs); 

    constraints.distribute(linear_solution); 

    return linear_solution; 
  } 
// @sect3{Details of the optimization algorithm}  

// 接下来的几个函数处理优化算法的具体部分，最主要的是决定通过求解线性化（牛顿）系统计算出的方向是否可行，如果可行，我们要在这个方向上走多远。

//  @sect4{Computing step lengths}  

// 我们先用一个函数进行二进制搜索，找出符合对偶可行性的最大步骤--也就是说，我们能走多远，使  $s>0$  和  $z>0$  。该函数返回一对数值，分别代表 $s$ 和 $z$ 的松弛变量。

  template <int dim> 
  std::pair<double, double> SANDTopOpt<dim>::calculate_max_step_size( 
    const BlockVector<double> &state, 
    const BlockVector<double> &step) const 
  { 
    double       fraction_to_boundary; 
    const double min_fraction_to_boundary = .8; 
    const double max_fraction_to_boundary = 1. - 1e-5; 

    if (min_fraction_to_boundary < 1 - barrier_size) 
      { 
        if (1 - barrier_size < max_fraction_to_boundary) 
          fraction_to_boundary = 1 - barrier_size; 
        else 
          fraction_to_boundary = max_fraction_to_boundary; 
      } 
    else 
      fraction_to_boundary = min_fraction_to_boundary; 

    double step_size_s_low  = 0; 
    double step_size_z_low  = 0; 
    double step_size_s_high = 1; 
    double step_size_z_high = 1; 
    double step_size_s, step_size_z; 

    const int max_bisection_method_steps = 50; 
    for (unsigned int k = 0; k < max_bisection_method_steps; ++k) 
      { 
        step_size_s = (step_size_s_low + step_size_s_high) / 2; 
        step_size_z = (step_size_z_low + step_size_z_high) / 2; 

        const BlockVector<double> state_test_s = 
          (fraction_to_boundary * state) + (step_size_s * step); 

        const BlockVector<double> state_test_z = 
          (fraction_to_boundary * state) + (step_size_z * step); 

        const bool accept_s = 
          (state_test_s.block(SolutionBlocks::density_lower_slack) 
             .is_non_negative()) && 
          (state_test_s.block(SolutionBlocks::density_upper_slack) 
             .is_non_negative()); 
        const bool accept_z = 
          (state_test_z.block(SolutionBlocks::density_lower_slack_multiplier) 
             .is_non_negative()) && 
          (state_test_z.block(SolutionBlocks::density_upper_slack_multiplier) 
             .is_non_negative()); 

        if (accept_s) 
          step_size_s_low = step_size_s; 
        else 
          step_size_s_high = step_size_s; 

        if (accept_z) 
          step_size_z_low = step_size_z; 
        else 
          step_size_z_high = step_size_z; 
      } 

    return {step_size_s_low, step_size_z_low}; 
  } 
// @sect4{Computing residuals}  

// 下一个函数计算一个围绕 "测试解向量 "线性化的右手向量，我们可以用它来观察KKT条件的大小。 然后，这将用于在缩小障碍大小之前测试收敛性，以及计算 $l_1$ 的优点。

// 这个函数冗长而复杂，但它实际上只是复制了上面`assemble_system()`函数的右侧部分的内容。

  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::calculate_test_rhs( 
    const BlockVector<double> &test_solution) const 
  { 

// 我们首先创建一个零向量，其大小和阻塞为system_rhs

    BlockVector<double> test_rhs; 
    test_rhs.reinit(system_rhs); 

    MappingQGeneric<dim>  mapping(1); 
    const QGauss<dim>     quadrature_formula(fe.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
    FEValues<dim>         fe_values(mapping, 
                            fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim>     fe_face_values(mapping, 
                                     fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_normal_vectors | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>     cell_rhs(dofs_per_cell); 
    FullMatrix<double> dummy_cell_matrix(dofs_per_cell, dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<double> lambda_values(n_q_points); 
    std::vector<double> mu_values(n_q_points); 

    const Functions::ConstantFunction<dim> lambda(1.), mu(1.); 
    std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 

    BlockVector<double> filtered_unfiltered_density_solution = test_solution; 
    BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
      test_solution; 
    filtered_unfiltered_density_solution.block( 
      SolutionBlocks::unfiltered_density) = 0; 
    filter_adjoint_unfiltered_density_multiplier_solution.block( 
      SolutionBlocks::unfiltered_density_multiplier) = 0; 

    filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
                          SolutionBlocks::unfiltered_density), 
                        test_solution.block( 
                          SolutionBlocks::unfiltered_density)); 
    filter_matrix.Tvmult( 
      filter_adjoint_unfiltered_density_multiplier_solution.block( 
        SolutionBlocks::unfiltered_density_multiplier), 
      test_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 

    std::vector<double>                  old_density_values(n_q_points); 
    std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
    std::vector<double>                  old_displacement_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
    std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
    std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
    std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
      n_q_points); 
    std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
    std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
    std::vector<double> old_lower_slack_values(n_q_points); 
    std::vector<double> old_upper_slack_values(n_q_points); 
    std::vector<double> old_unfiltered_density_values(n_q_points); 
    std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
    std::vector<double> filtered_unfiltered_density_values(n_q_points); 
    std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
      n_q_points); 

    using namespace ValueExtractors; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_rhs = 0; 

        cell->get_dof_indices(local_dof_indices); 

        fe_values.reinit(cell); 

 
        mu.value_list(fe_values.get_quadrature_points(), mu_values); 

        fe_values[densities<dim>].get_function_values(test_solution, 
                                                      old_density_values); 
        fe_values[displacements<dim>].get_function_values( 
          test_solution, old_displacement_values); 
        fe_values[displacements<dim>].get_function_divergences( 
          test_solution, old_displacement_divs); 
        fe_values[displacements<dim>].get_function_symmetric_gradients( 
          test_solution, old_displacement_symmgrads); 
        fe_values[displacement_multipliers<dim>].get_function_values( 
          test_solution, old_displacement_multiplier_values); 
        fe_values[displacement_multipliers<dim>].get_function_divergences( 
          test_solution, old_displacement_multiplier_divs); 
        fe_values[displacement_multipliers<dim>] 
          .get_function_symmetric_gradients( 
            test_solution, old_displacement_multiplier_symmgrads); 
        fe_values[density_lower_slacks<dim>].get_function_values( 
          test_solution, old_lower_slack_values); 
        fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
          test_solution, old_lower_slack_multiplier_values); 
        fe_values[density_upper_slacks<dim>].get_function_values( 
          test_solution, old_upper_slack_values); 
        fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
          test_solution, old_upper_slack_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          test_solution, old_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          test_solution, old_unfiltered_density_multiplier_values); 
        fe_values[unfiltered_densities<dim>].get_function_values( 
          filtered_unfiltered_density_solution, 
          filtered_unfiltered_density_values); 
        fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
          filter_adjoint_unfiltered_density_multiplier_solution, 
          filter_adjoint_unfiltered_density_multiplier_values); 

        for (const auto q_point : fe_values.quadrature_point_indices()) 
          { 
            for (const auto i : fe_values.dof_indices()) 
              { 
                const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
                  fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
                const double displacement_phi_i_div = 
                  fe_values[displacements<dim>].divergence(i, q_point); 

                const SymmetricTensor<2, dim> 
                  displacement_multiplier_phi_i_symmgrad = 
                    fe_values[displacement_multipliers<dim>].symmetric_gradient( 
                      i, q_point); 
                const double displacement_multiplier_phi_i_div = 
                  fe_values[displacement_multipliers<dim>].divergence(i, 
                                                                      q_point); 

                const double density_phi_i = 
                  fe_values[densities<dim>].value(i, q_point); 
                const double unfiltered_density_phi_i = 
                  fe_values[unfiltered_densities<dim>].value(i, q_point); 
                const double unfiltered_density_multiplier_phi_i = 
                  fe_values[unfiltered_density_multipliers<dim>].value(i, 
                                                                       q_point); 

                const double lower_slack_multiplier_phi_i = 
                  fe_values[density_lower_slack_multipliers<dim>].value( 
                    i, q_point); 

                const double lower_slack_phi_i = 
                  fe_values[density_lower_slacks<dim>].value(i, q_point); 

                const double upper_slack_phi_i = 
                  fe_values[density_upper_slacks<dim>].value(i, q_point); 

                const double upper_slack_multiplier_phi_i = 
                  fe_values[density_upper_slack_multipliers<dim>].value( 
                    i, q_point); 

                /* 方程1：这个方程以及方程
                 * 2 and 3, are the variational derivatives of the 
                 * Lagrangian with respect to the decision 
                 * variables - the density, displacement, and 
                 * unfiltered density. */ 


                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (density_penalty_exponent * 
                     std::pow(old_density_values[q_point], 
                              density_penalty_exponent - 1) * 
                     density_phi_i * 
                     (old_displacement_multiplier_divs[q_point] * 
                        old_displacement_divs[q_point] * 
                        lambda_values[q_point] + 
                      2 * mu_values[q_point] * 
                        (old_displacement_symmgrads[q_point] * 
                         old_displacement_multiplier_symmgrads[q_point])) - 
                   density_phi_i * 
                     old_unfiltered_density_multiplier_values[q_point]); 

                /*方程2；边界项将被进一步添加。
                 * below. */ 


                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (std::pow(old_density_values[q_point], 
                            density_penalty_exponent) * 
                   (old_displacement_multiplier_divs[q_point] * 
                      displacement_phi_i_div * lambda_values[q_point] + 
                    2 * mu_values[q_point] * 
                      (old_displacement_multiplier_symmgrads[q_point] * 
                       displacement_phi_i_symmgrad))); 
//           
               /* 方程3  */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (unfiltered_density_phi_i * 
                     filter_adjoint_unfiltered_density_multiplier_values 
                       [q_point] + 
                   unfiltered_density_phi_i * 
                     old_upper_slack_multiplier_values[q_point] + 
                   -1 * unfiltered_density_phi_i * 
                     old_lower_slack_multiplier_values[q_point]); 

               /* 方程4；边界项将再次被处理。with below. 
                * This equation being driven to 0 ensures that the elasticity 
                * equation is met as a constraint. */ 
                cell_rhs(i) += -1 * fe_values.JxW(q_point) * 
                               (std::pow(old_density_values[q_point], 
                                         density_penalty_exponent) * 
                                (old_displacement_divs[q_point] * 
                                   displacement_multiplier_phi_i_div * 
                                   lambda_values[q_point] + 
                                 2 * mu_values[q_point] * 
                                   (displacement_multiplier_phi_i_symmgrad * 
                                    old_displacement_symmgrads[q_point]))); 

                /* 方程5：该方程设定了下限的松弛量， giving a minimum density of 0. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (lower_slack_multiplier_phi_i * 
                                (old_unfiltered_density_values[q_point] - 
                                 old_lower_slack_values[q_point])); 

                /* 方程6：该方程设定了上层松弛量variable equal to one minus the unfiltered density. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (upper_slack_multiplier_phi_i * 
                                (1 - old_unfiltered_density_values[q_point] - 
                                 old_upper_slack_values[q_point])); 

                /*等式7：这是在
                 * density and the filter applied to the 
                 * unfiltered density. This being driven to 0 by 
                 * the Newton steps ensures that the filter is 
                 * applied correctly. */ 
                cell_rhs(i) += fe_values.JxW(q_point) * 
                               (unfiltered_density_multiplier_phi_i * 
                                (old_density_values[q_point] - 
                                 filtered_unfiltered_density_values[q_point])); 

                /*方程8：这与方程9一起给出了
                 * requirement that $s*z = \alpha$ for the barrier 
                 * size alpha, and gives complementary slackness 
                 * from KKT conditions when $\alpha$ goes to 0. */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (lower_slack_phi_i * 
                   (old_lower_slack_multiplier_values[q_point] - 
                    barrier_size / old_lower_slack_values[q_point])); 

                /*方程9  */ 
                cell_rhs(i) += 
                  -1 * fe_values.JxW(q_point) * 
                  (upper_slack_phi_i * 
                   (old_upper_slack_multiplier_values[q_point] - 
                    barrier_size / old_upper_slack_values[q_point])); 
              } 
          } 

        for (const auto &face : cell->face_iterators()) 
          { 
            if (face->at_boundary() && 
                face->boundary_id() == BoundaryIds::down_force) 
              { 
                fe_face_values.reinit(cell, face); 

                for (const auto face_q_point : 
                     fe_face_values.quadrature_point_indices()) 
                  { 
                    for (const auto i : fe_face_values.dof_indices()) 
                      { 
                        Tensor<1, dim> traction; 
                        traction[1] = -1.; 

                        cell_rhs(i) += 
                          -1 * 
                          (traction * fe_face_values[displacements<dim>].value( 
                                        i, face_q_point)) * 
                          fe_face_values.JxW(face_q_point); 

                        cell_rhs(i) += 
                          (traction * 
                           fe_face_values[displacement_multipliers<dim>].value( 
                             i, face_q_point)) * 
                          fe_face_values.JxW(face_q_point); 
                      } 
                  } 
              } 
          } 

        MatrixTools::local_apply_boundary_values(boundary_values, 
                                                 local_dof_indices, 
                                                 dummy_cell_matrix, 
                                                 cell_rhs, 
                                                 true); 

        constraints.distribute_local_to_global(cell_rhs, 
                                               local_dof_indices, 
                                               test_rhs); 
      } 

    return test_rhs; 
  } 
  // @sect4{Computing the merit function}  

  // 我们在这里使用的算法使用一个 "看门狗 "策略来确定从当前迭代的位置和程度。 我们将看门狗策略建立在一个精确的 $l_1$ 功绩函数上。这个函数计算一个给定的、假定的、下一个迭代的精确 $l_1$ 功绩。

  //优点函数由目标函数的总和（简单来说就是外力的积分（在域的边界上）乘以测试解的位移值（通常是当前解加上牛顿更新的某个倍数），以及残差向量的拉格朗日乘数分量的 $l_1$ 准则组成。下面的代码依次计算这些部分。

  template <int dim> 
  double SANDTopOpt<dim>::calculate_exact_merit( 
    const BlockVector<double> &test_solution) 
  { 
    TimerOutput::Scope t(timer, "merit function"); 

    // 从计算目标函数开始。
    double objective_function_merit = 0; 
    { 
      MappingQGeneric<dim>  mapping(1); 
      const QGauss<dim>     quadrature_formula(fe.degree + 1); 
      const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
      FEValues<dim>         fe_values(mapping, 
                              fe, 
                              quadrature_formula, 
                              update_values | update_gradients | 
                                update_quadrature_points | update_JxW_values); 
      FEFaceValues<dim>     fe_face_values(mapping, 
                                       fe, 
                                       face_quadrature_formula, 
                                       update_values | 
                                         update_quadrature_points | 
                                         update_normal_vectors | 
                                         update_JxW_values); 

      const unsigned int n_face_q_points = face_quadrature_formula.size(); 

      std::vector<Tensor<1, dim>> displacement_face_values(n_face_q_points); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          for (const auto &face : cell->face_iterators()) 
            { 
              if (face->at_boundary() && 
                  face->boundary_id() == BoundaryIds::down_force) 
                { 
                  fe_face_values.reinit(cell, face); 
                  fe_face_values[ValueExtractors::displacements<dim>] 
                    .get_function_values(test_solution, 
                                         displacement_face_values); 
                  for (unsigned int face_q_point = 0; 
                       face_q_point < n_face_q_points; 
                       ++face_q_point) 
                    { 
                      Tensor<1, dim> traction; 
                      traction[1] = -1.; 

                      objective_function_merit += 
                        (traction * displacement_face_values[face_q_point]) * 
                        fe_face_values.JxW(face_q_point); 
                    } 
                } 
            } 
        } 
    } 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        objective_function_merit = 
          objective_function_merit - 
          barrier_size * cell->measure() * 
            std::log(test_solution.block( 
              SolutionBlocks::density_lower_slack)[cell->active_cell_index()]); 
        objective_function_merit = 
          objective_function_merit - 
          barrier_size * cell->measure() * 
            std::log(test_solution.block( 
              SolutionBlocks::density_upper_slack)[cell->active_cell_index()]); 
      } 
    //然后
    //计算残差，并取对应于拉格朗日多边形的组件的 $l_1$ 准则。我们把这些加到上面计算的目标函数中，并在底部返回总和。
    const BlockVector<double> test_rhs = calculate_test_rhs(test_solution); 

    const double elasticity_constraint_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::displacement_multiplier).l1_norm(); 
    const double filter_constraint_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::unfiltered_density_multiplier).l1_norm(); 
    const double lower_slack_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::density_lower_slack_multiplier).l1_norm(); 
    const double upper_slack_merit = 
      penalty_multiplier * 
      test_rhs.block(SolutionBlocks::density_upper_slack_multiplier).l1_norm(); 

    const double total_merit = 
      objective_function_merit + elasticity_constraint_merit + 
      filter_constraint_merit + lower_slack_merit + upper_slack_merit; 
    return total_merit; 
  } 

  //  @sect4{Finding a search direction}  

  // 接下来是实际计算从当前状态（作为第一个参数传递）开始的搜索方向并返回结果向量的函数。为此，该函数首先调用与牛顿系统相对应的线性系统的组合函数，并对其进行求解。

  // 这个函数还更新了优点函数中的惩罚乘数，然后返回最大比例的可行步骤。它使用`calculate_max_step_sizes()`函数来找到满足  $s>0$  和  $z>0$  的最大可行步骤。

  template <int dim> 
  BlockVector<double> SANDTopOpt<dim>::find_max_step() 
  { 
    assemble_system(); 
    BlockVector<double> step = solve(); 

    // 接下来我们要更新punice_multiplier。 从本质上讲，更大的惩罚乘数使我们更多考虑约束条件。 观察与我们的决策变量有关的Hessian和梯度，并将其与我们的约束误差的规范相比较，可以确保我们的优点函数是 "精确的"

    // 也就是说，它在与目标函数相同的位置有一个最小值。 由于我们的优点函数对任何超过某个最小值的惩罚乘数都是精确的，所以我们只保留计算值，如果它增加了惩罚乘数。

    const std::vector<unsigned int> decision_variables = { 
      SolutionBlocks::density, 
      SolutionBlocks::displacement, 
      SolutionBlocks::unfiltered_density, 
      SolutionBlocks::density_upper_slack, 
      SolutionBlocks::density_lower_slack}; 
    double hess_part = 0; 
    double grad_part = 0; 
    for (const unsigned int decision_variable_i : decision_variables) 
      { 
        for (const unsigned int decision_variable_j : decision_variables) 
          { 
            Vector<double> temp_vector(step.block(decision_variable_i).size()); 
            system_matrix.block(decision_variable_i, decision_variable_j) 
              .vmult(temp_vector, step.block(decision_variable_j)); 
            hess_part += step.block(decision_variable_i) * temp_vector; 
          } 
        grad_part -= system_rhs.block(decision_variable_i) * 
                     step.block(decision_variable_i); 
      } 

    const std::vector<unsigned int> equality_constraint_multipliers = { 
      SolutionBlocks::displacement_multiplier, 
      SolutionBlocks::unfiltered_density_multiplier, 
      SolutionBlocks::density_lower_slack_multiplier, 
      SolutionBlocks::density_upper_slack_multiplier}; 
    double constraint_norm = 0; 
    for (unsigned int multiplier_i : equality_constraint_multipliers) 
      constraint_norm += system_rhs.block(multiplier_i).linfty_norm(); 

    double test_penalty_multiplier; 
    if (hess_part > 0) 
      test_penalty_multiplier = 
        (grad_part + .5 * hess_part) / (.05 * constraint_norm); 
    else 
      test_penalty_multiplier = (grad_part) / (.05 * constraint_norm); 

    penalty_multiplier = std::max(penalty_multiplier, test_penalty_multiplier); 

    // 基于所有这些，我们现在可以计算出原始变量和对偶变量（拉格朗日乘数）的步长。一旦我们有了这些，我们就可以对解向量的分量进行缩放，这就是这个函数的回报。

    const std::pair<double, double> max_step_sizes = 
      calculate_max_step_size(nonlinear_solution, step); 
    const double step_size_s = max_step_sizes.first; 
    const double step_size_z = max_step_sizes.second; 

    step.block(SolutionBlocks::density) *= step_size_s; 
    step.block(SolutionBlocks::displacement) *= step_size_s; 
    step.block(SolutionBlocks::unfiltered_density) *= step_size_s; 
    step.block(SolutionBlocks::displacement_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::unfiltered_density_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::density_lower_slack) *= step_size_s; 
    step.block(SolutionBlocks::density_lower_slack_multiplier) *= step_size_z; 
    step.block(SolutionBlocks::density_upper_slack) *= step_size_s; 
    step.block(SolutionBlocks::density_upper_slack_multiplier) *= step_size_z; 

    return step; 
  } 

  //  @sect4{Computing a scaled step}  

  // 下一个函数接着实现了直线搜索的反向跟踪算法。它不断缩小步长，直到找到一个优点减少的步长，然后根据当前的状态向量，以及要进入的方向，乘以步长，返回新的位置。

  template <int dim> 
  BlockVector<double> 
  SANDTopOpt<dim>::compute_scaled_step(const BlockVector<double> &state, 
                                       const BlockVector<double> &max_step, 
                                       const double descent_requirement) 
  { 
    const double merit_derivative = 
      (calculate_exact_merit(state + 1e-4 * max_step) - 
       calculate_exact_merit(state)) / 
      1e-4; 
    double       step_size                 = 1; 
    unsigned int max_linesearch_iterations = 10; 
    for (unsigned int k = 0; k < max_linesearch_iterations; ++k) 
      { 
        if (calculate_exact_merit(state + step_size * max_step) < 
            calculate_exact_merit(state) + 
              step_size * descent_requirement * merit_derivative) 
          break; 
        else 
          step_size = step_size / 2; 
      } 
    return state + (step_size * max_step); 
  } 

  // @sect4{Checking for convergence}  

  // 本块中的最后一个辅助函数是检查是否充分满足KKT条件，以便整个算法可以降低障碍物的大小。它通过计算残差的 $l_1$ 准则来实现，这就是`calculate_test_rhs()`的计算。

  template <int dim> 
  bool SANDTopOpt<dim>::check_convergence(const BlockVector<double> &state) 
  { 
    const BlockVector<double> test_rhs      = calculate_test_rhs(state); 
    const double              test_rhs_norm = test_rhs.l1_norm(); 

    const double convergence_condition = 1e-2; 
    const double target_norm           = convergence_condition * barrier_size; 

    std::cout << "    Checking convergence. Current rhs norm is " 
              << test_rhs_norm << ", target is " << target_norm << std::endl; 

    return (test_rhs_norm < target_norm); 
  } 

  // @sect3{Postprocessing the solution}  

  // 后处理函数中的第一个函数在VTU文件中输出信息，用于可视化。它看起来很长，但实际上与  step-22  中所做的一样，例如，只是增加了（很多）解决方案的变量。

  template <int dim> 
  void SANDTopOpt<dim>::output_results(const unsigned int iteration) const 
  { 
    std::vector<std::string> solution_names(1, "density"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        1, DataComponentInterpretation::component_is_scalar); 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        solution_names.emplace_back("displacement"); 
        data_component_interpretation.push_back( 
          DataComponentInterpretation::component_is_part_of_vector); 
      } 
    solution_names.emplace_back("unfiltered_density"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        solution_names.emplace_back("displacement_multiplier"); 
        data_component_interpretation.push_back( 
          DataComponentInterpretation::component_is_part_of_vector); 
      } 
    solution_names.emplace_back("unfiltered_density_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("low_slack"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("low_slack_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("high_slack"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    solution_names.emplace_back("high_slack_multiplier"); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(nonlinear_solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output("solution" + std::to_string(iteration) + ".vtu"); 
    data_out.write_vtu(output); 
  } 

  // 其中第二个函数将解决方案输出为`.stl`文件，用于3D打印。STL](https:en.wikipedia.org/wiki/STL_(file_format))文件是由三角形和法线向量组成的，我们将用它来显示所有那些密度值大于0的单元，首先将网格从 $z$ 值挤出到 $z=0.25$  ，然后为密度值足够大的单元的每个面生成两个三角形。当从外面看时，三角形节点必须逆时针走，法向量必须是指向外部的单位向量，这需要进行一些检查。
  template <int dim> 
  void SANDTopOpt<dim>::write_as_stl() 
  { 
    static_assert(dim == 2, 
                  "This function is not implemented for anything " 
                  "other than the 2d case."); 

    std::ofstream stlfile; 
    stlfile.open("bridge.stl"); 

    stlfile << "solid bridge\n" << std::scientific; 
    double height = .25; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        if (nonlinear_solution.block( 
              SolutionBlocks::density)[cell->active_cell_index()] > 0.5) 
          { 
            // 我们现在已经找到了一个密度值大于0的单元。让我们先写出底部和顶部的面。由于上面提到的排序问题，我们必须确保了解一个单元的坐标系是右旋的还是左旋的。我们通过询问从顶点0开始的两条边的方向以及它们是否形成一个右手坐标系来做到这一点。
            const Tensor<1, dim> edge_directions[2] = {cell->vertex(1) - 
                                                         cell->vertex(0), 
                                                       cell->vertex(2) - 
                                                         cell->vertex(0)}; 
            const Tensor<2, dim> edge_tensor( 
              {{edge_directions[0][0], edge_directions[0][1]}, 
               {edge_directions[1][0], edge_directions[1][1]}}); 
            const bool is_right_handed_cell = (determinant(edge_tensor) > 0); 

            if (is_right_handed_cell) 
              { 

               /*在z=0处写出一个边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 

               /*在z=高度处写下一个边。  */  
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
              } 
            else /* The cell has a left-handed set up */ 
              { 
               /* 在z=0处写出一边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " "
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 

               /*在z=高度处写出一个边。  */ 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(0)[0] << " " 
                        << cell->vertex(0)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
                stlfile << "   facet normal " << 0.000000e+00 << " " 
                        << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
                stlfile << "      outer loop\n"; 
                stlfile << "         vertex " << cell->vertex(1)[0] << " " 
                        << cell->vertex(1)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(2)[0] << " " 
                        << cell->vertex(2)[1] << " " << height << "\n"; 
                stlfile << "         vertex " << cell->vertex(3)[0] << " " 
                        << cell->vertex(3)[1] << " " << height << "\n"; 
                stlfile << "      endloop\n"; 
                stlfile << "   endfacet\n"; 
              } 

            // 接下来我们需要处理单元格的四个面，扩展到 $z$ 方向。然而，我们只需要写这些面，如果该面在域的边界上，或者它是密度大于0.5的单元和密度小于0.5的单元之间的界面。

            for (unsigned int face_number = 0; 
                 face_number < GeometryInfo<dim>::faces_per_cell; 
                 ++face_number) 
              { 
                const typename DoFHandler<dim>::face_iterator face = 
                  cell->face(face_number); 

                if ((face->at_boundary()) || 
                    (!face->at_boundary() && 
                     (nonlinear_solution.block( 
                        0)[cell->neighbor(face_number)->active_cell_index()] < 
                      0.5))) 
                  { 
                    const Tensor<1, dim> normal_vector = 
                      (face->center() - cell->center()); 
                    const double normal_norm = normal_vector.norm(); 
                    if ((face->vertex(0)[0] - face->vertex(0)[0]) * 
                            (face->vertex(1)[1] - face->vertex(0)[1]) * 
                            0.000000e+00 + 
                          (face->vertex(0)[1] - face->vertex(0)[1]) * (0 - 0) * 
                            normal_vector[0] + 
                          (height - 0) * 
                            (face->vertex(1)[0] - face->vertex(0)[0]) * 
                            normal_vector[1] - 
                          (face->vertex(0)[0] - face->vertex(0)[0]) * (0 - 0) * 
                            normal_vector[1] - 
                          (face->vertex(0)[1] - face->vertex(0)[1]) * 
                            (face->vertex(1)[0] - face->vertex(0)[0]) * 
                            normal_vector[0] - 
                          (height - 0) * 
                            (face->vertex(1)[1] - face->vertex(0)[1]) * 0 > 
                        0) 
                      { 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                      } 
                    else 
                      { 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                        stlfile << "   facet normal " 
                                << normal_vector[0] / normal_norm << " " 
                                << normal_vector[1] / normal_norm << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "      outer loop\n"; 
                        stlfile << "         vertex " << face->vertex(0)[0] 
                                << " " << face->vertex(0)[1] << " " << height 
                                << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " 
                                << 0.000000e+00 << "\n"; 
                        stlfile << "         vertex " << face->vertex(1)[0] 
                                << " " << face->vertex(1)[1] << " " << height 
                                << "\n"; 
                        stlfile << "      endloop\n"; 
                        stlfile << "   endfacet\n"; 
                      } 
                  } 
              } 
          } 
      } 
    stlfile << "endsolid bridge"; 
  } 

  //  @sect3{The run() function driving the overall algorithm}  

  // 这个函数最终提供了整体的驱动逻辑。从总体上看，这是一个相当复杂的函数，主要是因为优化算法很困难：它不仅仅是像 step-15 中那样找到一个牛顿方向，然后在这个方向上再走一个固定的距离，而是要（i）确定当前步骤中的最佳对数障碍惩罚参数应该是什么，（ii）通过复杂的算法来确定我们要走多远，还有其他成分。让我们看看如何在下面的文件中把它分解成小块。

  // 该函数一开始就很简单，首先设置了网格、DoFHandler，然后是下面所需的各种线性代数对象。

  template <int dim> 
  void SANDTopOpt<dim>::run() 
  { 
    std::cout << "filter r is: " << filter_r << std::endl; 

    { 
      TimerOutput::Scope t(timer, "setup"); 

      create_triangulation(); 

      dof_handler.distribute_dofs(fe); 
      DoFRenumbering::component_wise(dof_handler); 

      setup_boundary_values(); 
      setup_block_system(); 
      setup_filter_matrix(); 
    } 

  // 然后，我们设置一些影响优化算法的对数屏障和直线搜索部分的参数。

    barrier_size                  = 25; 
    const double min_barrier_size = .0005; 

    const unsigned int max_uphill_steps    = 8; 
    const double       descent_requirement = .0001; 

    // 现在开始进行主迭代。整个算法通过使用一个外循环来工作，在这个外循环中，我们一直循环到（i）对数障碍参数变得足够小，或者（ii）我们已经达到收敛。在任何情况下，如果最终的迭代次数过多，我们就会终止。这个整体结构被编码为一个 "do{ ... } while (...)`循环，其中收敛条件在底部。

    unsigned int       iteration_number = 0; 
    const unsigned int max_iterations   = 10000; 

    do 
      { 
        std::cout << "Starting outer step in iteration " << iteration_number 
                  << " with barrier parameter " << barrier_size << std::endl; 

        // 在这个外循环中，我们有一个内循环，在这个内循环中，我们试图使用介绍中描述的看门狗算法找到一个更新方向。

        // 看门狗算法本身的总体思路是这样的。对于最大的`max_uphill_steps`（即上述 "内循环 "中的一个循环）的尝试，我们使用`find_max_step()`来计算牛顿更新步骤，并在`nonlinear_solution`向量中加上这些。 在每一次尝试中（从上一次尝试结束时到达的地方开始），我们检查我们是否已经达到了上述优点函数的目标值。目标值是根据本算法的起始位置（看门狗循环开始时的`nonlinear_solution'，保存为`看门狗_state'）和本循环第一个回合中`find_max_step()'提供的第一个建议方向（`k=0'情况）计算的。

        do 
          { 
            std::cout << "  Starting inner step in iteration " 
                      << iteration_number 
                      << " with merit function penalty multiplier " 
                      << penalty_multiplier << std::endl; 

            bool watchdog_step_found = false; 

            const BlockVector<double> watchdog_state = nonlinear_solution; 
            BlockVector<double>       first_step; 
            double target_merit     = numbers::signaling_nan<double>(); 
            double merit_derivative = numbers::signaling_nan<double>(); 

            for (unsigned int k = 0; k < max_uphill_steps; ++k) 
              { 
                ++iteration_number; 
                const BlockVector<double> update_step = find_max_step(); 

                if (k == 0) 
                  { 
                    first_step = update_step; 
                    merit_derivative = 
                      ((calculate_exact_merit(watchdog_state + 
                                              .0001 * first_step) - 
                        calculate_exact_merit(watchdog_state)) / 
                       .0001); 
                    target_merit = calculate_exact_merit(watchdog_state) + 
                                   descent_requirement * merit_derivative; 
                  } 

                nonlinear_solution += update_step; 
                const double current_merit = 
                  calculate_exact_merit(nonlinear_solution); 

                std::cout << "    current watchdog state merit is: " 
                          << current_merit << "; target merit is " 
                          << target_merit << std::endl; 

                if (current_merit < target_merit) 
                  { 
                    watchdog_step_found = true; 
                    std::cout << "    found workable step after " << k + 1 
                              << " iterations" << std::endl; 
                    break; 
                  } 
              } 
            //然后
            //算法的下一部分取决于上面的看门狗循环是否成功。如果成功了，那么我们就满意了，不需要进一步的行动。我们只是停留在原地。然而，如果我们在上面的循环中采取了最大数量的不成功的步骤，那么我们就需要做一些别的事情，这就是下面的代码块所做的。    具体来说，从上述循环的最后（不成功的）状态开始，我们再寻找一个更新方向，并采取所谓的 "伸展步骤"。如果该拉伸状态满足涉及优点函数的条件，那么我们就去那里。另一方面，如果拉伸状态也是不可接受的（就像上面所有的看门狗步骤一样），那么我们就放弃上面所有的看门狗步骤，在我们开始看门狗迭代的地方重新开始--那个地方被存储在上面的`看门狗_状态`变量中。更具体地说，下面的条件首先测试我们是否从`看门狗_state`方向的`first_step`走了一步，或者我们是否可以从拉伸状态再做一次更新来找到一个新的地方。有可能这两种情况实际上都不比我们在看门狗算法开始时的状态好，但即使是这样，那个地方显然是个困难的地方，离开后从另一个地方开始下一次迭代可能是一个有用的策略，最终收敛。    我们不断重复上面的看门狗步骤以及下面的逻辑，直到这个内部迭代最终收敛（或者如果我们遇到最大的迭代次数--在这里我们把线性求解的次数算作迭代次数，并在每次调用`find_max_step()`时增加计数器，因为这就是线性求解实际发生的地方）。在任何情况下，在这些内部迭代的每一次结束时，我们也会以适合可视化的形式输出解决方案。

            if (watchdog_step_found == false) 
              { 
                ++iteration_number; 
                const BlockVector<double> update_step = find_max_step(); 
                const BlockVector<double> stretch_state = 
                  compute_scaled_step(nonlinear_solution, 
                                      update_step, 
                                      descent_requirement); 

                // 如果我们没有得到一个成功的看门狗步骤，我们现在需要决定是回到我们开始的地方，还是使用最终状态。 我们比较这两个位置的优劣，然后从哪个位置取一个按比例的步长。 由于按比例的步长可以保证降低优点，所以我们最终会保留这两个位置中的一个。

                if ((calculate_exact_merit(nonlinear_solution) < 
                     calculate_exact_merit(watchdog_state)) || 
                    (calculate_exact_merit(stretch_state) < target_merit)) 
                  { 
                    std::cout << "    Taking scaled step from end of watchdog" 
                              << std::endl; 
                    nonlinear_solution = stretch_state; 
                  } 
                else 
                  { 
                    std::cout 
                      << "    Taking scaled step from beginning of watchdog" 
                      << std::endl; 
                    if (calculate_exact_merit(stretch_state) > 
                        calculate_exact_merit(watchdog_state)) 
                      { 
                        nonlinear_solution = 
                          compute_scaled_step(watchdog_state, 
                                              first_step, 
                                              descent_requirement); 
                      } 
                    else 
                      { 
                        ++iteration_number; 
                        nonlinear_solution = stretch_state; 
                        const BlockVector<double> stretch_step = 
                          find_max_step(); 
                        nonlinear_solution = 
                          compute_scaled_step(nonlinear_solution, 
                                              stretch_step, 
                                              descent_requirement); 
                      } 
                  } 
              } 

            output_results(iteration_number); 
          } 
        while ((iteration_number < max_iterations) && 
               (check_convergence(nonlinear_solution) == false)); 

        // 在外循环结束时，我们必须更新屏障参数，为此我们使用以下公式。该函数的其余部分只是检查外循环的收敛条件，如果我们决定终止计算，就把最终的 "设计 "写成STL文件，用于3D打印，并输出一些时间信息。

        const double barrier_size_multiplier = .8; 
        const double barrier_size_exponent   = 1.2; 

        barrier_size = 
          std::max(std::min(barrier_size * barrier_size_multiplier, 
                            std::pow(barrier_size, barrier_size_exponent)), 
                   min_barrier_size); 

        std::cout << std::endl; 
      } 
    while (((barrier_size > min_barrier_size) || 
            (check_convergence(nonlinear_solution) == false)) && 
           (iteration_number < max_iterations)); 

    write_as_stl(); 
    timer.print_summary(); 
  } 
} // namespace SAND 
// @sect3{The main function}  

// 余下的代码，即`main()`函数，和平常一样。

int main() 
{ 
  try 
    { 
      SAND::SANDTopOpt<2> elastic_problem_2d; 
      elastic_problem_2d.run(); 
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
                << "Aborting!" << std::endl;
      return 1;
    }
  return 0;
}
