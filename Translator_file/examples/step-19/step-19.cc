

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2020 - 2021 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Rene Gassmoeller, Peter Munch, 2020. 
 */ 


// @sect3{Include files}  

// 本程序中使用的大部分include文件都是 step-6 和类似程序中众所周知的。

#include <deal.II/base/quadrature_lib.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/fe/mapping_q.h> 
#include <deal.II/matrix_free/fe_point_evaluation.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/error_estimator.h> 

// 新的只有以下三个。第一个声明了DiscreteTime类，它帮助我们在时间相关的模拟中跟踪时间。后面两个提供了所有的粒子功能，即记录位于网格上的粒子的方法（ Particles::ParticleHandler 类）和为可视化目的输出这些粒子的位置及其属性的能力（ Particles::DataOut  类）。

#include <deal.II/base/discrete_time.h> 
#include <deal.II/particles/particle_handler.h> 
#include <deal.II/particles/data_out.h> 

#include <fstream> 

using namespace dealii; 
// @sect3{Global definitions}  

// 按照惯例，我们把所有与程序细节相对应的东西都放到一个自己的命名空间中。在顶部，我们定义了一些常量，我们宁愿使用符号名称而不是硬编码的数字。

// 具体来说，我们为几何学的各个部分定义了 @ref GlossBoundaryIndicator "边界指标 "的数字，以及电子的物理属性和我们在这里使用的其他具体设置。

// 对于边界指标，让我们从某个随机值101开始列举。这里的原则是要使用*不常见的数字。如果之前有`GridGenerator'函数设置的预定义边界指标，它们很可能是从0开始的小整数，但不是在这个相当随机的范围内。使用下面这样的数字可以避免冲突的可能性，同时也减少了在程序中直接拼出这些数字的诱惑（因为你可能永远不会记得哪个是哪个，而如果它们从0开始，你可能会受到诱惑）。

namespace Step19 
{ 
  namespace BoundaryIds 
  { 
    constexpr types::boundary_id open          = 101; 
    constexpr types::boundary_id cathode       = 102; 
    constexpr types::boundary_id focus_element = 103; 
    constexpr types::boundary_id anode         = 104; 
  } // namespace BoundaryIds 

  namespace Constants 
  { 
    constexpr double electron_mass   = 9.1093837015e-31; 
    constexpr double electron_charge = 1.602176634e-19; 

    constexpr double V0 = 1; 

    constexpr double E_threshold = 0.05; 

    constexpr double electrons_per_particle = 3e15; 
  } // namespace Constants 
// @sect3{The main class}  

// 然后，下面是这个程序的主类。从根本上说，它的结构与 step-6 和其他许多教程程序相同。这包括大部分的成员函数（其余部分的目的可能从它们的名字中不难看出），以及超出 step-6 的少量成员变量，所有这些都与处理粒子有关。

  template <int dim> 
  class CathodeRaySimulator 
  { 
  public: 
    CathodeRaySimulator(); 

    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve_field(); 
    void refine_grid(); 

    void create_particles(); 
    void move_particles(); 
    void track_lost_particle( 
      const typename Particles::ParticleIterator<dim> &        particle, 
      const typename Triangulation<dim>::active_cell_iterator &cell); 

    void update_timestep_size(); 
    void output_results() const; 

    Triangulation<dim>        triangulation; 
    MappingQGeneric<dim>      mapping; 
    FE_Q<dim>                 fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 

    SparseMatrix<double> system_matrix; 
    SparsityPattern      sparsity_pattern; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    Particles::ParticleHandler<dim> particle_handler; 
    types::particle_index           next_unused_particle_id; 
    types::particle_index           n_recently_lost_particles; 
    types::particle_index           n_total_lost_particles; 
    types::particle_index           n_particles_lost_through_anode; 

    DiscreteTime time; 
  }; 

//  @sect3{The <code>CathodeRaySimulator</code> class implementation}  
// @sect4{The <code>CathodeRaySimulator</code> constructor}  

// 那么，让我们开始执行。构造函数所做的实际上只是对顶部的所有成员变量进行简单的初始化。唯一值得一提的是`particle_handler'，它被交给了一个指向粒子所在的三角形的引用（目前当然还是空的，但是粒子处理程序存储了这个引用，一旦粒子被添加，就会使用它--这发生在三角形被构建之后）。它得到的另一个信息是每个粒子需要存储多少 "属性"。在这里，我们需要每个粒子记住的是它当前的速度，也就是一个带有`dim`分量的矢量。然而，每个粒子还有其他的内在属性， Particles::ParticleHandler 类会自动并始终确保这些属性是可用的；特别是，这些属性是粒子的当前位置、它所在的单元格、它在该单元格中的参考位置，以及粒子的ID。

// 唯一感兴趣的其他变量是 "时间"，一个DiscreteTime类型的对象。它记录了我们在一个随时间变化的模拟中的当前时间，并以开始时间（零）和结束时间（ $10^{-4}$ ）初始化。我们以后将在`update_timestep_size()`中设置时间步长。

// 构造函数的主体由我们在介绍中已经讨论过的一段代码组成。也就是说，我们要确保每次有粒子离开域时，`track_lost_particle()`函数都会被`particle_handler`对象调用。

  template <int dim> 
  CathodeRaySimulator<dim>::CathodeRaySimulator() 
    : mapping(1) 
    , fe(2) 
    , dof_handler(triangulation) 
    , particle_handler(triangulation, mapping, /*n_properties=*/dim) 
    , next_unused_particle_id(0) 
    , n_recently_lost_particles(0) 
    , n_total_lost_particles(0) 
    , n_particles_lost_through_anode(0) 
    , time(0, 1e-4) 
  { 
    particle_handler.signals.particle_lost.connect( 
      [this](const typename Particles::ParticleIterator<dim> &        particle, 
             const typename Triangulation<dim>::active_cell_iterator &cell) { 
        this->track_lost_particle(particle, cell); 
      }); 
  } 

//  @sect4{The <code>CathodeRaySimulator::make_grid</code> function}  

// 下一个函数是负责生成我们要解决的网格。回顾一下域的样子。    
// <p align="center">
//      <img
//      src="https:www.dealii.org/images/steps/developer/step-19.geometry.png"
//           alt="The geometry used in this program"
//           width="600">
//    </p>  我们把这个几何体细分为 $4\times 2$ 个单元的网格，看起来像这样。
//    @code
//    *---*---*---*---*
//    \   |   |   |   |
//     *--*---*---*---*
//    /   |   |   |   |
//    *---*---*---*---*
//  @endcode 
//  这样做的方法是首先定义 $15=5\times 3$ 顶点的位置--在这里，我们说它们在整数点上，左边的中间点向右移动了`delta=0.5`的值。

// 在下文中，我们必须说明哪些顶点共同组成了8个单元。下面的代码就完全等同于我们在 step-14 中的做法。

  template <int dim> 
  void CathodeRaySimulator<dim>::make_grid() 
  { 
    static_assert(dim == 2, 
                  "This function is currently only implemented for 2d."); 

    const double       delta = 0.5; 
    const unsigned int nx    = 5; 
    const unsigned int ny    = 3; 

    const std::vector<Point<dim>> vertices // 
      = {{0, 0}, 
         {1, 0}, 
         {2, 0}, 
         {3, 0}, 
         {4, 0}, 
         {delta, 1}, 
         {1, 1}, 
         {2, 1}, 
         {3, 1}, 
         {4, 1}, 
         {0, 2}, 
         {1, 2}, 
         {2, 2}, 
         {3, 2}, 
         {4, 2}}; 
    AssertDimension(vertices.size(), nx * ny); 

    const std::vector<unsigned int> cell_vertices[(nx - 1) * (ny - 1)] = { 
      {0, 1, nx + 0, nx + 1}, 
      {1, 2, nx + 1, nx + 2}, 
      {2, 3, nx + 2, nx + 3}, 
      {3, 4, nx + 3, nx + 4}, 

      {5, nx + 1, 2 * nx + 0, 2 * nx + 1}, 
      {nx + 1, nx + 2, 2 * nx + 1, 2 * nx + 2}, 
      {nx + 2, nx + 3, 2 * nx + 2, 2 * nx + 3}, 
      {nx + 3, nx + 4, 2 * nx + 3, 2 * nx + 4}}; 

// 有了这些数组，我们可以转向稍高的高层数据结构。我们创建一个CellData对象的向量，为每个要创建的单元存储相关的顶点以及 @ref GlossMaterialId "材料ID"（我们在这里将其简单地设置为0，因为我们在程序中不使用它）。

// 然后，这些信息将被传递给 Triangulation::create_triangulation() 函数，并对网格进行两次全局细化。

    std::vector<CellData<dim>> cells((nx - 1) * (ny - 1), CellData<dim>()); 
    for (unsigned int i = 0; i < cells.size(); ++i) 
      { 
        cells[i].vertices    = cell_vertices[i]; 
        cells[i].material_id = 0; 
      } 

    triangulation.create_triangulation( 
      vertices, 
      cells, 
      SubCellData()); // No boundary information 

    triangulation.refine_global(2); 

// 该函数的其余部分循环所有的单元格和它们的面，如果一个面在边界上，则决定哪个边界指标应该应用于它。如果你将代码与上面的几何图形相比较，各种条件应该是有意义的。

// 一旦完成了这一步，我们再全局地细化一下网格。

    for (auto &cell : triangulation.active_cell_iterators()) 
      for (auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          { 
            if ((face->center()[0] > 0) && (face->center()[0] < 0.5) && 
                (face->center()[1] > 0) && (face->center()[1] < 2)) 
              face->set_boundary_id(BoundaryIds::cathode); 
            else if ((face->center()[0] > 0) && (face->center()[0] < 2)) 
              face->set_boundary_id(BoundaryIds::focus_element); 
            else if ((face->center()[0] > 4 - 1e-12) && 
                     ((face->center()[1] > 1.5) || (face->center()[1] < 0.5))) 
              face->set_boundary_id(BoundaryIds::anode); 
            else 
              face->set_boundary_id(BoundaryIds::open); 
          } 

    triangulation.refine_global(1); 
  } 
// @sect4{The <code>CathodeRaySimulator::setup_system</code> function}  

// 本程序中的下一个函数是处理与解决偏微分方程有关的各种对象的设置。它本质上是对 step-6 中相应函数的复制，不需要进一步讨论。

  template <int dim> 
  void CathodeRaySimulator<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::cathode, 
                                             Functions::ConstantFunction<dim>( 
                                               -Constants::V0), 
                                             constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::focus_element, 
                                             Functions::ConstantFunction<dim>( 
                                               -Constants::V0), 
                                             constraints); 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             BoundaryIds::anode, 
                                             Functions::ConstantFunction<dim>( 
                                               +Constants::V0), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    dsp, 
                                    constraints, 
                                    /*keep_constrained_dofs =  */ false);

    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
  } 
// @sect4{The <code>CathodeRaySimulator::assemble_system</code> function}  

// 计算矩阵项的函数实质上还是复制了  step-6  中的相应函数。

  template <int dim> 
  void CathodeRaySimulator<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.dofs_per_cell; 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        fe_values.reinit(cell); 

        for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
          for (const unsigned int i : fe_values.dof_indices()) 
            { 
              for (const unsigned int j : fe_values.dof_indices()) 
                cell_matrix(i, j) += 
                  (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
                   fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
                   fe_values.JxW(q_index));           // dx 
            } 

// 这个函数唯一有趣的部分是它是如何形成线性系统的右手边的。回顾一下，PDE的右边是
// @f[
//    \sum_p (N e)\delta(\mathbf x-\mathbf x_p),
//  @f]
//  ，在这里我们用 $p$ 来索引粒子，以避免与形状函数 $\varphi_i$ 混淆； $\mathbf x_p$ 是第 $p$ 个粒子的位置。

// 当与测试函数 $\varphi_i$ 相乘并在域上积分时，会得到一个右手边的向量
// @f{align*}{
//    F_i &= \int_\Omega \varphi_i (\mathbf x)\left[
//                 \sum_p (N e)\delta(\mathbf x-\mathbf x_p) \right] dx
//    \\  &=  \sum_p (N e) \varphi_i(\mathbf x_p).
//  @f} 
//  注意最后一行不再包含一个积分，因此也没有出现 $dx$ ，这需要在我们的代码中出现`JxW`符号。
// 
// 对于一个给定的单元 $K$ ，这个单元对右边的贡献是
//  @f{align*}{
//    F_i^K &= \sum_{p, \mathbf x_p\in K} (N e) \varphi_i(\mathbf x_p),
//  @f}，
//  也就是说，我们只需要担心那些实际位于当前单元 $K$ 上的粒子。

// 在实践中，我们在这里所做的是以下几点。如果当前单元格上有任何粒子，那么我们首先获得一个迭代器范围，指向该单元格的第一个粒子以及该单元格上最后一个粒子之后的粒子（或结束迭代器）--即C++函数中常见的半开放范围。现在知道了粒子的列表，我们查询它们的参考位置（相对于参考单元），评估这些参考位置的形状函数，并根据上面的公式计算力（没有任何  FEValues::JxW).  ）。
// @note  值得指出的是，调用 Particles::ParticleHandler::particles_in_cell() 和 Particles::ParticleHandler::n_particles_in_cell() 函数在有大量粒子的问题上不是很有效。但是它说明了写这个算法的最简单的方法，所以我们愿意为了说明问题而暂时承担这个代价。  我们在下面的<a href="#extensions">"possibilities for extensions" section</a>中更详细地讨论了这个问题，并在 step-70 中使用了一个更好的方法，例如：。

        if (particle_handler.n_particles_in_cell(cell) > 0) 
          for (const auto &particle : particle_handler.particles_in_cell(cell)) 
            { 
              const Point<dim> &reference_location = 
                particle.get_reference_location(); 
              for (const unsigned int i : fe_values.dof_indices()) 
                cell_rhs(i) += 
                  (fe.shape_value(i, reference_location) * // phi_i(x_p) 
                   (-Constants::electrons_per_particle *   // N 
                    Constants::electron_charge));          // e 
            } 

// 最后，我们可以把这个单元格的贡献复制到全局矩阵和右边的向量中。

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global( 
          cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
      } 
  } 
// @sect4{CathodeRaySimulator::solve}  

// 解决线性系统的函数又与 step-6 中的完全一样。

  template <int dim> 
  void CathodeRaySimulator<dim>::solve_field() 
  { 
    SolverControl            solver_control(1000, 1e-12); 
    SolverCG<Vector<double>> solver(solver_control); 

    PreconditionSSOR<SparseMatrix<double>> preconditioner; 
    preconditioner.initialize(system_matrix, 1.2); 

    solver.solve(system_matrix, solution, system_rhs, preconditioner); 

    constraints.distribute(solution); 
  } 
// @sect4{CathodeRaySimulator::refine_grid}  

// 最后一个与场相关的函数是细化网格的函数。我们将在第一个时间步骤中多次调用它，以获得一个能很好地适应解的结构的网格，特别是解决解中由于重心角和边界条件类型变化的地方而产生的各种奇异现象。你可能想再参考一下 step-6 以了解更多的细节。

  template <int dim> 
  void CathodeRaySimulator<dim>::refine_grid() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate(dof_handler, 
                                       QGauss<dim - 1>(fe.degree + 1), 
                                       {}, 
                                       solution, 
                                       estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.1, 
                                                    0.03); 

    triangulation.execute_coarsening_and_refinement(); 
  } 
// @sect4{CathodeRaySimulator::create_particles}  

// 现在让我们来看看处理粒子的函数。第一个是关于粒子的创建。正如介绍中提到的，如果电场 $\mathbf E=\nabla V$ 超过某个阈值，即如果 $|\mathbf E| \ge E_\text{threshold}$ ，并且如果电场进一步指向域内（即如果 $\mathbf E \cdot \mathbf n < 0$ ），我们希望在阴极的各点创建一个粒子。正如有限元方法中常见的那样，我们在特定的评估点评估场（及其导数）；通常，这些是 "正交点"，因此我们创建了一个 "正交公式"，我们将用它来指定我们要评估解决方案的点。在这里，我们将简单地采用QMidpoint，意味着我们将只在面的中点检查阈值条件。然后我们用它来初始化一个FEFaceValues类型的对象来评估这些点的解。

// 然后，所有这些将被用于所有单元格、它们的面，特别是那些位于边界的面，而且是边界的阴极部分的循环中。

  template <int dim> 
  void CathodeRaySimulator<dim>::create_particles() 
  { 
    FEFaceValues<dim> fe_face_values(fe, 
                                     QMidpoint<dim - 1>(), 
                                     update_quadrature_points | 
                                       update_gradients | 
                                       update_normal_vectors); 

    std::vector<Tensor<1, dim>> solution_gradients( 
      fe_face_values.n_quadrature_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary() && 
            (face->boundary_id() == BoundaryIds::cathode)) 
          { 
            fe_face_values.reinit(cell, face); 

// 所以我们已经找到了阴极上的一个面。接下来，我们让FEFaceValues对象计算每个 "正交 "点的解的梯度，并通过 @ref vector_valued "矢量值问题 "文件模块中讨论的方法，以张量变量的形式从梯度中提取电场向量。

            const FEValuesExtractors::Scalar electric_potential(0); 
            fe_face_values[electric_potential].get_function_gradients( 
              solution, solution_gradients); 
            for (const unsigned int q_point : 
                 fe_face_values.quadrature_point_indices()) 
              { 
                const Tensor<1, dim> E = solution_gradients[q_point]; 

// 只有当电场强度超过阈值时，电子才能逃离阴极，而且关键是，如果电场指向*域内，电子才能逃离阴极。      一旦我们检查了这一点，我们就在这个位置创建一个新的 Particles::Particle 对象，并将其插入到 Particles::ParticleHandler 对象中，并设置一个唯一的ID。            这里唯一不明显的是，我们还将这个粒子与我们当前所在的单元格的参考坐标中的位置联系起来。这样做是因为我们将在下游函数中计算诸如粒子位置的电场等量（例如，在每个时间步长中更新其位置时计算作用于它的力）。在任意坐标上评估有限元场是一个相当昂贵的操作，因为形状函数实际上只定义在参考单元上，所以当要求一个任意点的电场时，我们首先要确定这个点的参考坐标是什么。为了避免反复操作，我们一次性地确定这些坐标，然后将这些参考坐标直接存储在粒子上。

                if ((E * fe_face_values.normal_vector(q_point) < 0) && 
                    (E.norm() > Constants::E_threshold)) 
                  { 
                    const Point<dim> &location = 
                      fe_face_values.quadrature_point(q_point); 

                    Particles::Particle<dim> new_particle; 
                    new_particle.set_location(location); 
                    new_particle.set_reference_location( 
                      mapping.transform_real_to_unit_cell(cell, location)); 
                    new_particle.set_id(next_unused_particle_id); 
                    particle_handler.insert_particle(new_particle, cell); 

                    ++next_unused_particle_id; 
                  } 
              } 
          } 

// 在所有这些插入结束时，我们让`particle_handler`更新它所存储的粒子的一些内部统计数据。

    particle_handler.update_cached_numbers(); 
  } 
// @sect4{CathodeRaySimulator::move_particles}  

// 第二个与粒子有关的函数是在每个时间步骤中移动粒子的函数。要做到这一点，我们必须在所有的单元格、每个单元格中的粒子上循环，并评估每个粒子位置的电场。

// 这里使用的方法在概念上与`assemble_system()`函数中使用的相同。我们在所有单元中循环，找到位于那里的粒子（同样要注意这里用来寻找这些粒子的算法的低效率），并使用FEPointEvaluation对象来评估这些位置的梯度。

  template <int dim> 
  void CathodeRaySimulator<dim>::move_particles() 
  { 
    const double dt = time.get_next_step_size(); 

    Vector<double>            solution_values(fe.n_dofs_per_cell()); 
    FEPointEvaluation<1, dim> evaluator(mapping, fe, update_gradients); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (particle_handler.n_particles_in_cell(cell) > 0) 
        { 
          const typename Particles::ParticleHandler< 
            dim>::particle_iterator_range particles_in_cell = 
            particle_handler.particles_in_cell(cell); 

          std::vector<Point<dim>> particle_positions; 
          for (const auto &particle : particles_in_cell) 
            particle_positions.push_back(particle.get_reference_location()); 

          cell->get_dof_values(solution, solution_values); 

// 然后，我们可以向FEPointEvaluation对象询问这些位置的解决方案的梯度（即电场 $\mathbf E$ ），并在各个粒子上循环。

          evaluator.reinit(cell, particle_positions); 
          evaluator.evaluate(make_array_view(solution_values), 
                             EvaluationFlags::gradients); 

          { 
            typename Particles::ParticleHandler<dim>::particle_iterator 
              particle = particles_in_cell.begin(); 
            for (unsigned int particle_index = 0; 
                 particle != particles_in_cell.end(); 
                 ++particle, ++particle_index) 
              { 
                const Tensor<1, dim> &E = 
                  evaluator.get_gradient(particle_index); 

// 现在我们已经得到了其中一个粒子位置的电场，我们首先用它来更新速度，然后更新位置。为此，我们首先从粒子的属性中获取旧的速度，计算加速度，更新速度，并将这个新的速度再次存储在粒子的属性中。回顾一下，这对应于介绍中所讨论的以下一组更新方程中的第一个。      
      //  @f{align*}{
      //      \frac{{\mathbf v}_i^{(n)}
      //            -{\mathbf v}_i^{(n-1)}}{\Delta t}
      //      &= \frac{e\nabla V^{(n)}}{m}
      //   \\ \frac{{\mathbf x}_i^{(n)}-{\mathbf x}_i^{(n-1)}}
      //           {\Delta t} &= {\mathbf v}_i^{(n)}.
      //  @f}

                const Tensor<1, dim> old_velocity(particle->get_properties()); 

                const Tensor<1, dim> acceleration = 
                  Constants::electron_charge / Constants::electron_mass * E; 

                const Tensor<1, dim> new_velocity = 
                  old_velocity + acceleration * dt; 

                particle->set_properties(make_array_view(new_velocity)); 

// 有了新的速度，我们也就可以更新粒子的位置，并告诉粒子这个位置。

                const Point<dim> new_location = 
                  particle->get_location() + dt * new_velocity; 
                particle->set_location(new_location); 
              } 
          } 
        } 

// 在更新了所有粒子的位置和属性（即速度）之后，我们需要确保`particle_handler`再次知道它们在哪个单元中，以及它们在参考单元坐标系中的位置。下面的函数就是这样做的。(它还确保在并行计算中，如果粒子从一个处理器拥有的子域移动到另一个处理器拥有的子域，那么粒子会从一个处理器移动到另一个处理器。)

    particle_handler.sort_particles_into_subdomains_and_cells(); 
  } 
// @sect4{CathodeRaySimulator::track_lost_particle}  

// 最后一个与粒子相关的函数是当一个粒子从模拟中丢失时被调用的函数。这通常发生在它离开域的时候。如果发生这种情况，这个函数会同时调用单元（我们可以询问它的新位置）和它之前所在的单元。然后，该函数不断跟踪更新这个时间步骤中丢失的粒子数，丢失的粒子总数，然后估计该粒子是否通过阳极中间的孔离开。我们这样做，首先检查它最后所在的单元是否有一个 $x$ 坐标在右边边界的左边（位于 $x=4$ ），而粒子现在的位置在右边边界的右边。如果是这样的话，我们就计算出它的运动方向矢量，这个方向矢量被归一化了，所以方向矢量的 $x$ 分量等于 $1$  。有了这个方向矢量，我们可以计算出它与直线 $x=4$ 的相交位置。如果这个相交点在 $0.5$ 和 $1.5$ 之间，那么我们就声称粒子从孔中离开，并增加一个计数器。

  template <int dim> 
  void CathodeRaySimulator<dim>::track_lost_particle( 
    const typename Particles::ParticleIterator<dim> &        particle, 
    const typename Triangulation<dim>::active_cell_iterator &cell) 
  { 
    ++n_recently_lost_particles; 
    ++n_total_lost_particles; 

    const Point<dim> current_location              = particle->get_location(); 
    const Point<dim> approximate_previous_location = cell->center(); 

    if ((approximate_previous_location[0] < 4) && (current_location[0] > 4)) 
      { 
        const Tensor<1, dim> direction = 
          (current_location - approximate_previous_location) / 
          (current_location[0] - approximate_previous_location[0]); 

        const double right_boundary_intercept = 
          approximate_previous_location[1] + 
          (4 - approximate_previous_location[0]) * direction[1]; 
        if ((right_boundary_intercept > 0.5) && 
            (right_boundary_intercept < 1.5)) 
          ++n_particles_lost_through_anode; 
      } 
  } 

//  @sect4{CathodeRaySimulator::update_timestep_size}  

// 正如在介绍中详细讨论的那样，我们需要尊重一个时间步长条件，即颗粒在一个时间步长中不能移动超过一个单元。为了确保这一点，我们首先计算每个单元上所有粒子的最大速度，然后用该速度除以单元大小。然后，我们使用介绍中讨论的安全系数，将下一个时间步长计算为所有单元上这个量的最小值，并使用 DiscreteTime::set_desired_time_step_size() 函数将其设定为所需的时间步长。

  template <int dim> 
  void CathodeRaySimulator<dim>::update_timestep_size() 
  { 
    if (time.get_step_number() > 0) 
      { 
        double min_cell_size_over_velocity = std::numeric_limits<double>::max(); 

        for (const auto &cell : dof_handler.active_cell_iterators()) 
          if (particle_handler.n_particles_in_cell(cell) > 0) 
            { 
              const double cell_size = cell->minimum_vertex_distance(); 

              double max_particle_velocity(0.0); 

              for (const auto &particle : 
                   particle_handler.particles_in_cell(cell)) 
                { 
                  const Tensor<1, dim> velocity(particle.get_properties()); 
                  max_particle_velocity = 
                    std::max(max_particle_velocity, velocity.norm()); 
                } 

              if (max_particle_velocity > 0) 
                min_cell_size_over_velocity = 
                  std::min(min_cell_size_over_velocity, 
                           cell_size / max_particle_velocity); 
            } 

        constexpr double c_safety = 0.5; 
        time.set_desired_next_step_size(c_safety * 0.5 * 
                                        min_cell_size_over_velocity); 
      } 

// 正如在介绍中提到的，我们必须以不同的方式对待第一个时间步长，因为在那里，粒子还没有出现，或者还没有我们计算合理步长所需的相关信息。下面的公式遵循介绍中的讨论。

    else 
      { 
        const QTrapezoid<dim> vertex_quadrature; 
        FEValues<dim> fe_values(fe, vertex_quadrature, update_gradients); 

        std::vector<Tensor<1, dim>> field_gradients(vertex_quadrature.size()); 

        double min_timestep = std::numeric_limits<double>::max(); 

        for (const auto &cell : dof_handler.active_cell_iterators()) 
          if (particle_handler.n_particles_in_cell(cell) > 0) 
            { 
              const double cell_size = cell->minimum_vertex_distance(); 

              fe_values.reinit(cell); 
              fe_values.get_function_gradients(solution, field_gradients); 

              double max_E = 0; 
              for (const auto q_point : fe_values.quadrature_point_indices()) 
                max_E = std::max(max_E, field_gradients[q_point].norm()); 

              if (max_E > 0) 
                min_timestep = 
                  std::min(min_timestep, 
                           std::sqrt(0.5 * cell_size * 
                                     Constants::electron_mass / 
                                     Constants::electron_charge / max_E)); 
            } 

        time.set_desired_next_step_size(min_timestep); 
      } 
  } 

//  @sect4{The <code>CathodeRaySimulator::output_results()</code> function}  

// 实现整个算法的最后一个函数是生成图形输出的函数。在目前的情况下，我们想同时输出电势场以及粒子的位置和速度。但我们也想输出电场，即解决方案的梯度。

// deal.II有一个一般的方法，可以从解决方案中计算出派生量，并输出这些量。在这里，这是电场，但也可以是其他的量--比如说，电场的法线，或者事实上任何其他人们想从解 $V_h(\mathbf x)$ 或其导数中计算的量。这个一般的解决方案使用了DataPostprocessor类，在像这里的情况下，我们想输出一个代表矢量场的量，则使用DataPostprocessorVector类。

// 与其尝试解释这个类是如何工作的，不如让我们简单地参考一下DataPostprocessorVector类的文档，这个案例基本上是一个有据可查的例子。

  template <int dim> 
  class ElectricFieldPostprocessor : public DataPostprocessorVector<dim> 
  { 
  public: 
    ElectricFieldPostprocessor() 
      : DataPostprocessorVector<dim>("electric_field", update_gradients) 
    {} 

    virtual void evaluate_scalar_field( 
      const DataPostprocessorInputs::Scalar<dim> &input_data, 
      std::vector<Vector<double>> &computed_quantities) const override 
    { 
      AssertDimension(input_data.solution_gradients.size(), 
                      computed_quantities.size()); 

      for (unsigned int p = 0; p < input_data.solution_gradients.size(); ++p) 
        { 
          AssertDimension(computed_quantities[p].size(), dim); 
          for (unsigned int d = 0; d < dim; ++d) 
            computed_quantities[p][d] = input_data.solution_gradients[p][d]; 
        } 
    } 
  }; 

// 有了这个，`output_results()`函数就变得相对简单了。我们使用DataOut类，就像我们在以前几乎所有的教程程序中使用的那样，来输出解决方案（"电动势"），我们使用上面定义的后处理程序来输出其梯度（"电场"）。这些都被写入一个VTU格式的文件中，同时将当前时间和时间步长与该文件联系起来。

  template <int dim> 
  void CathodeRaySimulator<dim>::output_results() const 
  { 
    { 
      ElectricFieldPostprocessor<dim> electric_field; 
      DataOut<dim>                    data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "electric_potential"); 
      data_out.add_data_vector(solution, electric_field); 
      data_out.build_patches(); 

      data_out.set_flags( 
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number())); 

      std::ofstream output("solution-" + 
                           Utilities::int_to_string(time.get_step_number(), 4) + 
                           ".vtu"); 
      data_out.write_vtu(output); 
    } 

// 输出粒子的位置和属性并不复杂。 Particles::DataOut 类扮演了粒子的DataOut类的角色，我们所要做的就是告诉该类从哪里获取粒子，以及如何解释属性中的`dim`分量--即作为表示速度的单一矢量，而不是作为`dim`标量属性。剩下的就和上面一样了。

    { 
      Particles::DataOut<dim, dim> particle_out; 
      particle_out.build_patches( 
        particle_handler, 
        std::vector<std::string>(dim, "velocity"), 
        std::vector<DataComponentInterpretation::DataComponentInterpretation>( 
          dim, DataComponentInterpretation::component_is_part_of_vector)); 

      particle_out.set_flags( 
        DataOutBase::VtkFlags(time.get_current_time(), time.get_step_number())); 

      std::ofstream output("particles-" + 
                           Utilities::int_to_string(time.get_step_number(), 4) + 
                           ".vtu"); 
      particle_out.write_vtu(output); 
    } 
  } 
// @sect4{CathodeRaySimulator::run}  

// 这个程序的主类的最后一个成员函数是驱动。在顶层，它通过在一连串越来越细的网格上求解问题（尚未创建粒子），对网格进行多次细化。

  template <int dim> 
  void CathodeRaySimulator<dim>::run() 
  { 
    make_grid(); 

//在前面做几个细化循环

    const unsigned int n_pre_refinement_cycles = 3; 
    for (unsigned int refinement_cycle = 0; 
         refinement_cycle < n_pre_refinement_cycles; 
         ++refinement_cycle) 
      { 
        setup_system(); 
        assemble_system(); 
        solve_field(); 
        refine_grid(); 
      } 

// 现在进行时间上的循环。这个步骤的顺序紧跟介绍中讨论的算法大纲。正如在DiscreteTime类的文档中详细讨论的那样，虽然我们将场和粒子信息向前移动了一个时间步长，但存储在`time`变量中的时间与这些量的（部分）位置不一致（在DiscreteTime的字典中，这就是 "更新阶段"）。对`time.advance_time()`的调用通过将`time`变量设置为场和粒子已经处于的时间而使一切重新保持一致，一旦我们处于这个 "一致阶段"，我们就可以生成图形输出并将模拟的当前状态的信息写入屏幕。

    setup_system(); 
    do 
      { 
        std::cout << "Timestep " << time.get_step_number() + 1 << std::endl; 
        std::cout << "  Field degrees of freedom:                 " 
                  << dof_handler.n_dofs() << std::endl; 

        assemble_system(); 
        solve_field(); 

        create_particles(); 
        std::cout << "  Total number of particles in simulation:  " 
                  << particle_handler.n_global_particles() << std::endl; 

        n_recently_lost_particles = 0; 
        update_timestep_size(); 
        move_particles(); 

        time.advance_time(); 

        output_results(); 

        std::cout << "  Number of particles lost this time step:  " 
                  << n_recently_lost_particles << std::endl; 
        if (n_total_lost_particles > 0) 
          std::cout << "  Fraction of particles lost through anode: " 
                    << 1. * n_particles_lost_through_anode / 
                         n_total_lost_particles 
                    << std::endl; 

        std::cout << std::endl 
                  << "  Now at t=" << time.get_current_time() 
                  << ", dt=" << time.get_previous_step_size() << '.' 
                  << std::endl 
                  << std::endl; 
      } 
    while (time.is_at_end() == false); 
  } 
} // namespace Step19 

//  @sect3{The <code>main</code> function}  

// 程序的最后一个函数又是`main()`函数。自 step-6 以来，它在所有的教程程序中都没有变化，因此没有什么新的内容需要讨论。

int main() 
{ 
  try 
    { 
      Step19::CathodeRaySimulator<2> cathode_ray_simulator_2d; 
      cathode_ray_simulator_2d.run(); 
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


