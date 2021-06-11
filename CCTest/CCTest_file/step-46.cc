CCTest_file/step-46.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2011 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2011 
 */ 


// @sect3{Include files}  

// 这个程序的包含文件与之前许多其他程序的包含文件是一样的。唯一的新文件是在介绍中讨论的声明FE_Nothing的文件。hp目录下的文件已经在  step-27  中讨论过了。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_nothing.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/hp/fe_collection.h> 
#include <deal.II/hp/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <iostream> 
#include <fstream> 

namespace Step46 
{ 
  using namespace dealii; 
// @sect3{The <code>FluidStructureProblem</code> class template}  

// 这是主类。如果你想的话，它是 step-8 和 step-22 的组合，因为它的成员变量要么针对全局问题（Triangulation和DoFHandler对象，以及 hp::FECollection 和各种线性代数对象），要么与弹性或斯托克斯子问题有关。然而，该类的一般结构与其他大多数实现静止问题的程序一样。

// 有几个不言自明的辅助函数（<code>cell_is_in_fluid_domain, cell_is_in_solid_domain</code>）（对两个子域的符号名称进行操作，这些名称将被用作属于子域的单元的 material_ids。正如介绍中所解释的那样）和几个函数（<code>make_grid, set_active_fe_indices, assemble_interface_terms</code>），这些函数已经从其他的函数中分离出来，可以在其他的教程程序中找到，我们将在实现它们的时候讨论。

// 最后一组变量 (  <code>viscosity, lambda, eta</code>  ) 描述了用于两个物理模型的材料属性。

  template <int dim> 
  class FluidStructureProblem 
  { 
  public: 
    FluidStructureProblem(const unsigned int stokes_degree, 
                          const unsigned int elasticity_degree); 
    void run(); 

  private: 
    enum 
    { 
      fluid_domain_id, 
      solid_domain_id 
    }; 

    static bool cell_is_in_fluid_domain( 
      const typename DoFHandler<dim>::cell_iterator &cell); 

    static bool cell_is_in_solid_domain( 
      const typename DoFHandler<dim>::cell_iterator &cell); 

    void make_grid(); 
    void set_active_fe_indices(); 
    void setup_dofs(); 
    void assemble_system(); 
    void assemble_interface_term( 
      const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
      const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
      std::vector<Tensor<1, dim>> &         elasticity_phi, 
      std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
      std::vector<double> &                 stokes_phi_p, 
      FullMatrix<double> &                  local_interface_matrix) const; 
    void solve(); 
    void output_results(const unsigned int refinement_cycle) const; 
    void refine_mesh(); 

    const unsigned int stokes_degree; 
    const unsigned int elasticity_degree; 

    Triangulation<dim>    triangulation; 
    FESystem<dim>         stokes_fe; 
    FESystem<dim>         elasticity_fe; 
    hp::FECollection<dim> fe_collection; 
    DoFHandler<dim>       dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 

    const double viscosity; 
    const double lambda; 
    const double mu; 
  }; 
// @sect3{Boundary values and right hand side}  

// 下面这个类如其名。速度的边界值分别为2d的 
//  $\mathbf u=(0, \sin(\pi x))^T$ 和3d的 $\mathbf u=(0,
//  0, \sin(\pi x)\sin(\pi y))^T$ 。
//  这个问题的其余边界条件都是同质的，在介绍中已经讨论过。右边的强迫项对于流体和固体都是零，所以我们不需要为它设置额外的类。

  template <int dim> 
  class StokesBoundaryValues : public Function<dim> 
  { 
  public: 
    StokesBoundaryValues() 
      : Function<dim>(dim + 1 + dim) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  double StokesBoundaryValues<dim>::value(const Point<dim> & p, 
                                          const unsigned int component) const 
  { 
    Assert(component < this->n_components, 
           ExcIndexRange(component, 0, this->n_components)); 

    if (component == dim - 1) 
      switch (dim) 
        { 
          case 2: 
            return std::sin(numbers::PI * p[0]); 
          case 3: 
            return std::sin(numbers::PI * p[0]) * std::sin(numbers::PI * p[1]); 
          default: 
            Assert(false, ExcNotImplemented()); 
        } 

    return 0; 
  } 

  template <int dim> 
  void StokesBoundaryValues<dim>::vector_value(const Point<dim> &p, 
                                               Vector<double> &  values) const 
  { 
    for (unsigned int c = 0; c < this->n_components; ++c) 
      values(c) = StokesBoundaryValues<dim>::value(p, c); 
  } 

//  @sect3{The <code>FluidStructureProblem</code> implementation}  
// @sect4{Constructors and helper functions}  

// 现在我们来谈谈这个程序的主类的实现。最初的几个函数是构造函数和辅助函数，可以用来确定一个单元格在域的哪个部分。鉴于介绍中对这些主题的讨论，它们的实现是相当明显的。在构造函数中，注意我们必须从斯托克斯和弹性的基本元素中构造 hp::FECollection 对象；使用 hp::FECollection::push_back 函数在这个集合中为它们分配了0和1的位置，我们必须记住这个顺序，并在程序的其余部分一致使用。

  template <int dim> 
  FluidStructureProblem<dim>::FluidStructureProblem( 
    const unsigned int stokes_degree, 
    const unsigned int elasticity_degree) 
    : stokes_degree(stokes_degree) 
    , elasticity_degree(elasticity_degree) 
    , triangulation(Triangulation<dim>::maximum_smoothing) 
    , stokes_fe(FE_Q<dim>(stokes_degree + 1), 
                dim, 
                FE_Q<dim>(stokes_degree), 
                1, 
                FE_Nothing<dim>(), 
                dim) 
    , elasticity_fe(FE_Nothing<dim>(), 
                    dim, 
                    FE_Nothing<dim>(), 
                    1, 
                    FE_Q<dim>(elasticity_degree), 
                    dim) 
    , dof_handler(triangulation) 
    , viscosity(2) 
    , lambda(1) 
    , mu(1) 
  { 
    fe_collection.push_back(stokes_fe); 
    fe_collection.push_back(elasticity_fe); 
  } 

  template <int dim> 
  bool FluidStructureProblem<dim>::cell_is_in_fluid_domain( 
    const typename DoFHandler<dim>::cell_iterator &cell) 
  { 
    return (cell->material_id() == fluid_domain_id); 
  } 

  template <int dim> 
  bool FluidStructureProblem<dim>::cell_is_in_solid_domain( 
    const typename DoFHandler<dim>::cell_iterator &cell) 
  { 
    return (cell->material_id() == solid_domain_id); 
  } 
// @sect4{Meshes and assigning subdomains}  

// 接下来的一对函数是处理生成网格，并确保所有表示子域的标志都是正确的。  <code>make_grid</code>  ，正如在介绍中所讨论的，生成一个 $8\times 8$ 的网格（或者一个 $8\times 8\times 8$ 的三维网格）以确保每个粗略的网格单元完全在一个子域内。生成这个网格后，我们在其边界上循环，并在顶部边界设置边界指标为1，这是我们设置非零迪里希特边界条件的唯一地方。在这之后，我们再次在所有单元上循环，设置材料指标&mdash;用来表示我们处于域的哪一部分，是流体指标还是固体指标。

  template <int dim> 
  void FluidStructureProblem<dim>::make_grid() 
  { 
    GridGenerator::subdivided_hyper_cube(triangulation, 8, -1, 1); 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary() && (face->center()[dim - 1] == 1)) 
          face->set_all_boundary_ids(1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (((std::fabs(cell->center()[0]) < 0.25) && 
           (cell->center()[dim - 1] > 0.5)) || 
          ((std::fabs(cell->center()[0]) >= 0.25) && 
           (cell->center()[dim - 1] > -0.5))) 
        cell->set_material_id(fluid_domain_id); 
      else 
        cell->set_material_id(solid_domain_id); 
  } 

// 这对函数的第二部分决定在每个单元上使用哪个有限元。上面我们设置了每个粗略网格单元的材料指标，正如在介绍中提到的，这个信息在网格细化时将从母单元继承到子单元。

// 换句话说，只要我们细化（或创建）了网格，我们就可以依靠材料指示器来正确描述一个单元所处的域的哪一部分。然后我们利用这一点将单元的活动FE索引设置为该类的 hp::FECollection 成员变量中的相应元素：流体单元为0，固体单元为1。

  template <int dim> 
  void FluidStructureProblem<dim>::set_active_fe_indices() 
  { 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        if (cell_is_in_fluid_domain(cell)) 
          cell->set_active_fe_index(0); 
        else if (cell_is_in_solid_domain(cell)) 
          cell->set_active_fe_index(1); 
        else 
          Assert(false, ExcNotImplemented()); 
      } 
  } 
// @sect4{<code>FluidStructureProblem::setup_dofs</code>}  

// 下一步是为线性系统设置数据结构。为此，我们首先要用上面的函数设置活动FE指数，然后分配自由度，再确定线性系统的约束。后者包括像往常一样的悬挂节点约束，但也包括顶部流体边界的不均匀边界值，以及沿固体子域周边的零边界值。

  template <int dim> 
  void FluidStructureProblem<dim>::setup_dofs() 
  { 
    set_active_fe_indices(); 
    dof_handler.distribute_dofs(fe_collection); 

    { 
      constraints.clear(); 
      DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

      const FEValuesExtractors::Vector velocities(0); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               1, 
                                               StokesBoundaryValues<dim>(), 
                                               constraints, 
                                               fe_collection.component_mask( 
                                                 velocities)); 

      const FEValuesExtractors::Vector displacements(dim + 1); 
      VectorTools::interpolate_boundary_values( 
        dof_handler, 
        0, 
        Functions::ZeroFunction<dim>(dim + 1 + dim), 
        constraints, 
        fe_collection.component_mask(displacements)); 
    } 

// 不过，我们还需要处理更多的约束条件：我们必须确保在流体和固体的界面上速度为零。下面这段代码已经在介绍中介绍过了。

    { 
      std::vector<types::global_dof_index> local_face_dof_indices( 
        stokes_fe.n_dofs_per_face()); 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (cell_is_in_fluid_domain(cell)) 
          for (const auto face_no : cell->face_indices()) 
            if (cell->face(face_no)->at_boundary() == false) 
              { 
                bool face_is_on_interface = false; 

                if ((cell->neighbor(face_no)->has_children() == false) && 
                    (cell_is_in_solid_domain(cell->neighbor(face_no)))) 
                  face_is_on_interface = true; 
                else if (cell->neighbor(face_no)->has_children() == true) 
                  { 
                    for (unsigned int sf = 0; 
                         sf < cell->face(face_no)->n_children(); 
                         ++sf) 
                      if (cell_is_in_solid_domain( 
                            cell->neighbor_child_on_subface(face_no, sf))) 
                        { 
                          face_is_on_interface = true; 
                          break; 
                        } 
                  } 

                if (face_is_on_interface) 
                  { 
                    cell->face(face_no)->get_dof_indices(local_face_dof_indices, 
                                                         0); 
                    for (unsigned int i = 0; i < local_face_dof_indices.size(); 
                         ++i) 
                      if (stokes_fe.face_system_to_component_index(i).first < 
                          dim) 
                        constraints.add_line(local_face_dof_indices[i]); 
                  } 
              } 
    } 

// 在这一切结束后，我们可以向约束对象声明，我们现在已经准备好了所有的约束，并且该对象可以重建其内部数据结构以提高效率。

    constraints.close(); 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

// 在这个函数的其余部分，我们创建了一个在介绍中广泛讨论的稀疏模式，并使用它来初始化矩阵；然后还将向量设置为正确的大小。

    { 
      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 

      Table<2, DoFTools::Coupling> cell_coupling(fe_collection.n_components(), 
                                                 fe_collection.n_components()); 
      Table<2, DoFTools::Coupling> face_coupling(fe_collection.n_components(), 
                                                 fe_collection.n_components()); 

      for (unsigned int c = 0; c < fe_collection.n_components(); ++c) 
        for (unsigned int d = 0; d < fe_collection.n_components(); ++d) 
          { 
            if (((c < dim + 1) && (d < dim + 1) && 
                 !((c == dim) && (d == dim))) || 
                ((c >= dim + 1) && (d >= dim + 1))) 
              cell_coupling[c][d] = DoFTools::always; 

            if ((c >= dim + 1) && (d < dim + 1)) 
              face_coupling[c][d] = DoFTools::always; 
          } 

      DoFTools::make_flux_sparsity_pattern(dof_handler, 
                                           dsp, 
                                           cell_coupling, 
                                           face_coupling); 
      constraints.condense(dsp); 
      sparsity_pattern.copy_from(dsp); 
    } 

    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 

//  @sect4{<code>FluidStructureProblem::assemble_system</code>}  

// 下面是这个程序的中心函数：组装线性系统的函数。它在开始时有一长段设置辅助函数的内容：从创建正交公式到设置FEValues、FEFaceValues和FESubfaceValues对象，这些都是整合单元项以及界面项所必需的，以应对界面上的单元以相同大小或不同细化程度聚集在一起的情况...

  template <int dim> 
  void FluidStructureProblem<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim> stokes_quadrature(stokes_degree + 2); 
    const QGauss<dim> elasticity_quadrature(elasticity_degree + 2); 

    hp::QCollection<dim> q_collection; 
    q_collection.push_back(stokes_quadrature); 
    q_collection.push_back(elasticity_quadrature); 

    hp::FEValues<dim> hp_fe_values(fe_collection, 
                                   q_collection, 
                                   update_values | update_quadrature_points | 
                                     update_JxW_values | update_gradients); 

    const QGauss<dim - 1> common_face_quadrature( 
      std::max(stokes_degree + 2, elasticity_degree + 2)); 

    FEFaceValues<dim>    stokes_fe_face_values(stokes_fe, 
                                            common_face_quadrature, 
                                            update_JxW_values | 
                                              update_gradients | update_values); 
    FEFaceValues<dim>    elasticity_fe_face_values(elasticity_fe, 
                                                common_face_quadrature, 
                                                update_normal_vectors | 
                                                  update_values); 
    FESubfaceValues<dim> stokes_fe_subface_values(stokes_fe, 
                                                  common_face_quadrature, 
                                                  update_JxW_values | 
                                                    update_gradients | 
                                                    update_values); 
    FESubfaceValues<dim> elasticity_fe_subface_values(elasticity_fe, 
                                                      common_face_quadrature, 
                                                      update_normal_vectors | 
                                                        update_values); 

// ...描述局部对全局线性系统贡献所需的对象...

    const unsigned int stokes_dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
    const unsigned int elasticity_dofs_per_cell = 
      elasticity_fe.n_dofs_per_cell(); 

    FullMatrix<double> local_matrix; 
    FullMatrix<double> local_interface_matrix(elasticity_dofs_per_cell, 
                                              stokes_dofs_per_cell); 
    Vector<double>     local_rhs; 

    std::vector<types::global_dof_index> local_dof_indices; 
    std::vector<types::global_dof_index> neighbor_dof_indices( 
      stokes_dofs_per_cell); 

    const Functions::ZeroFunction<dim> right_hand_side(dim + 1); 

// ...到变量，允许我们提取形状函数的某些成分并缓存它们的值，而不是在每个正交点重新计算它们。

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Vector displacements(dim + 1); 

    std::vector<SymmetricTensor<2, dim>> stokes_symgrad_phi_u( 
      stokes_dofs_per_cell); 
    std::vector<double> stokes_div_phi_u(stokes_dofs_per_cell); 
    std::vector<double> stokes_phi_p(stokes_dofs_per_cell); 

    std::vector<Tensor<2, dim>> elasticity_grad_phi(elasticity_dofs_per_cell); 
    std::vector<double>         elasticity_div_phi(elasticity_dofs_per_cell); 
    std::vector<Tensor<1, dim>> elasticity_phi(elasticity_dofs_per_cell); 

// 然后是所有单元格的主循环，和 step-27 一样，初始化当前单元格的 hp::FEValues 对象，提取适合当前单元格的FEValues对象。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        hp_fe_values.reinit(cell); 

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values(); 

        local_matrix.reinit(cell->get_fe().n_dofs_per_cell(), 
                            cell->get_fe().n_dofs_per_cell()); 
        local_rhs.reinit(cell->get_fe().n_dofs_per_cell()); 

// 做完这些后，我们继续为属于斯托克斯和弹性区域的单元组装单元项。虽然我们原则上可以在一个公式中完成，实际上就是实现了介绍中所说的双线性形式，但我们意识到，我们的有限元空间的选择方式是，在每个单元上，有一组变量（速度和压力，或者位移）总是为零，因此，计算局部积分的更有效的方法是，根据测试我们处于域的哪一部分的 <code>if</code> 条款，只做必要的事情。

// 局部矩阵的实际计算与 step-22 以及 @ref vector_valued 文件模块中给出的弹性方程的计算相同。

        if (cell_is_in_fluid_domain(cell)) 
          { 
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
            Assert(dofs_per_cell == stokes_dofs_per_cell, ExcInternalError()); 

            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
              { 
                for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                  { 
                    stokes_symgrad_phi_u[k] = 
                      fe_values[velocities].symmetric_gradient(k, q); 
                    stokes_div_phi_u[k] = 
                      fe_values[velocities].divergence(k, q); 
                    stokes_phi_p[k] = fe_values[pressure].value(k, q); 
                  } 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    local_matrix(i, j) += 
                      (2 * viscosity * stokes_symgrad_phi_u[i] * 
                         stokes_symgrad_phi_u[j] - 
                       stokes_div_phi_u[i] * stokes_phi_p[j] - 
                       stokes_phi_p[i] * stokes_div_phi_u[j]) * 
                      fe_values.JxW(q); 
              } 
          } 
        else 
          { 
            const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell(); 
            Assert(dofs_per_cell == elasticity_dofs_per_cell, 
                   ExcInternalError()); 

            for (unsigned int q = 0; q < fe_values.n_quadrature_points; ++q) 
              { 
                for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                  { 
                    elasticity_grad_phi[k] = 
                      fe_values[displacements].gradient(k, q); 
                    elasticity_div_phi[k] = 
                      fe_values[displacements].divergence(k, q); 
                  } 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    { 
                      local_matrix(i, j) += 
                        (lambda * elasticity_div_phi[i] * 
                           elasticity_div_phi[j] + 
                         mu * scalar_product(elasticity_grad_phi[i], 
                                             elasticity_grad_phi[j]) + 
                         mu * 
                           scalar_product(elasticity_grad_phi[i], 
                                          transpose(elasticity_grad_phi[j]))) * 
                        fe_values.JxW(q); 
                    } 
              } 
          } 

// 一旦我们得到了单元积分的贡献，我们就把它们复制到全局矩阵中（通过 AffineConstraints::distribute_local_to_global 函数，立即处理约束）。请注意，我们没有向 <code>local_rhs</code> 变量中写入任何东西，尽管我们仍然需要传递它，因为消除非零边界值需要修改局部，因此也需要修改全局的右手值。

        local_dof_indices.resize(cell->get_fe().n_dofs_per_cell()); 
        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global(local_matrix, 
                                               local_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               system_rhs); 

// 这个函数更有趣的部分是我们看到关于两个子域之间的界面上的脸部条款。为此，我们首先要确保我们只组装一次，即使在所有单元的所有面的循环中会遇到界面的每一部分两次。我们武断地决定，只有当当前单元是固体子域的一部分，并且因此一个面不在边界上，并且它后面的潜在邻居是流体域的一部分时，我们才会评估界面条款。让我们从这些条件开始。

        if (cell_is_in_solid_domain(cell)) 
          for (const auto f : cell->face_indices()) 
            if (cell->face(f)->at_boundary() == false) 
              { 

// 在这一点上，我们知道当前的单元格是一个候选的整合对象，并且面 <code>f</code> 后面存在一个邻居。现在有三种可能性。           

// - 邻居处于同一细化水平，并且没有孩子。     

// - 邻居有子女。     

// - 邻居比较粗糙。            在所有这三种情况下，我们只对它感兴趣，如果它是流体子域的一部分。因此，让我们从第一种最简单的情况开始：如果邻居处于同一层次，没有子女，并且是一个流体单元，那么这两个单元共享一个边界，这个边界是界面的一部分，我们想沿着这个边界整合界面项。我们所要做的就是用当前面和邻接单元的面初始化两个FEFaceValues对象（注意我们是如何找出邻接单元的哪个面与当前单元接壤的），然后把东西传给评估界面项的函数（这个函数的第三个到第五个参数为它提供了抓取数组）。然后，结果再次被复制到全局矩阵中，使用一个知道本地矩阵的行和列的DoF指数来自不同单元的函数。

                if ((cell->neighbor(f)->level() == cell->level()) && 
                    (cell->neighbor(f)->has_children() == false) && 
                    cell_is_in_fluid_domain(cell->neighbor(f))) 
                  { 
                    elasticity_fe_face_values.reinit(cell, f); 
                    stokes_fe_face_values.reinit(cell->neighbor(f), 
                                                 cell->neighbor_of_neighbor(f)); 

                    assemble_interface_term(elasticity_fe_face_values, 
                                            stokes_fe_face_values, 
                                            elasticity_phi, 
                                            stokes_symgrad_phi_u, 
                                            stokes_phi_p, 
                                            local_interface_matrix); 

                    cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
                    constraints.distribute_local_to_global( 
                      local_interface_matrix, 
                      local_dof_indices, 
                      neighbor_dof_indices, 
                      system_matrix); 
                  } 

// 第二种情况是，如果邻居还有更多的孩子。在这种情况下，我们必须在邻居的所有子女中进行循环，看他们是否属于流体子域的一部分。如果它们是，那么我们就在共同界面上进行整合，这个界面是邻居的一个面和当前单元的一个子面，要求我们对邻居使用FEFaceValues，对当前单元使用FESubfaceValues。

                else if ((cell->neighbor(f)->level() == cell->level()) && 
                         (cell->neighbor(f)->has_children() == true)) 
                  { 
                    for (unsigned int subface = 0; 
                         subface < cell->face(f)->n_children(); 
                         ++subface) 
                      if (cell_is_in_fluid_domain( 
                            cell->neighbor_child_on_subface(f, subface))) 
                        { 
                          elasticity_fe_subface_values.reinit(cell, f, subface); 
                          stokes_fe_face_values.reinit( 
                            cell->neighbor_child_on_subface(f, subface), 
                            cell->neighbor_of_neighbor(f)); 

                          assemble_interface_term(elasticity_fe_subface_values, 
                                                  stokes_fe_face_values, 
                                                  elasticity_phi, 
                                                  stokes_symgrad_phi_u, 
                                                  stokes_phi_p, 
                                                  local_interface_matrix); 

                          cell->neighbor_child_on_subface(f, subface) 
                            ->get_dof_indices(neighbor_dof_indices); 
                          constraints.distribute_local_to_global( 
                            local_interface_matrix, 
                            local_dof_indices, 
                            neighbor_dof_indices, 
                            system_matrix); 
                        } 
                  } 

// 最后一个选项是，邻居比较粗大。在这种情况下，我们必须为邻居使用一个FESubfaceValues对象，为当前单元使用一个FEFaceValues；其余部分与之前相同。

                else if (cell->neighbor_is_coarser(f) && 
                         cell_is_in_fluid_domain(cell->neighbor(f))) 
                  { 
                    elasticity_fe_face_values.reinit(cell, f); 
                    stokes_fe_subface_values.reinit( 
                      cell->neighbor(f), 
                      cell->neighbor_of_coarser_neighbor(f).first, 
                      cell->neighbor_of_coarser_neighbor(f).second); 

                    assemble_interface_term(elasticity_fe_face_values, 
                                            stokes_fe_subface_values, 
                                            elasticity_phi, 
                                            stokes_symgrad_phi_u, 
                                            stokes_phi_p, 
                                            local_interface_matrix); 

                    cell->neighbor(f)->get_dof_indices(neighbor_dof_indices); 
                    constraints.distribute_local_to_global( 
                      local_interface_matrix, 
                      local_dof_indices, 
                      neighbor_dof_indices, 
                      system_matrix); 
                  } 
              } 
      } 
  } 

// 在组装全局系统的函数中，我们将计算接口条款传递给我们在此讨论的一个单独的函数。关键是，尽管我们无法预测FEFaceValues和FESubfaceValues对象的组合，但它们都是从FEFaceValuesBase类派生出来的，因此我们不必在意：该函数被简单地调用，有两个这样的对象表示面的两边的正交点上的形状函数值。然后我们做我们一直在做的事情：我们用形状函数的值和它们的导数来填充从头数组，然后循环计算矩阵的所有条目来计算局部积分。我们在这里评估的双线性形式的细节在介绍中给出。

  template <int dim> 
  void FluidStructureProblem<dim>::assemble_interface_term( 
    const FEFaceValuesBase<dim> &         elasticity_fe_face_values, 
    const FEFaceValuesBase<dim> &         stokes_fe_face_values, 
    std::vector<Tensor<1, dim>> &         elasticity_phi, 
    std::vector<SymmetricTensor<2, dim>> &stokes_symgrad_phi_u, 
    std::vector<double> &                 stokes_phi_p, 
    FullMatrix<double> &                  local_interface_matrix) const 
  { 
    Assert(stokes_fe_face_values.n_quadrature_points == 
             elasticity_fe_face_values.n_quadrature_points, 
           ExcInternalError()); 
    const unsigned int n_face_quadrature_points = 
      elasticity_fe_face_values.n_quadrature_points; 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Vector displacements(dim + 1); 

    local_interface_matrix = 0; 
    for (unsigned int q = 0; q < n_face_quadrature_points; ++q) 
      { 
        const Tensor<1, dim> normal_vector = 
          elasticity_fe_face_values.normal_vector(q); 

        for (unsigned int k = 0; k < stokes_fe_face_values.dofs_per_cell; ++k) 
          { 
            stokes_symgrad_phi_u[k] = 
              stokes_fe_face_values[velocities].symmetric_gradient(k, q); 
            stokes_phi_p[k] = stokes_fe_face_values[pressure].value(k, q); 
          } 
        for (unsigned int k = 0; k < elasticity_fe_face_values.dofs_per_cell; 
             ++k) 
          elasticity_phi[k] = 
            elasticity_fe_face_values[displacements].value(k, q); 

        for (unsigned int i = 0; i < elasticity_fe_face_values.dofs_per_cell; 
             ++i) 
          for (unsigned int j = 0; j < stokes_fe_face_values.dofs_per_cell; ++j) 
            local_interface_matrix(i, j) += 
              -((2 * viscosity * (stokes_symgrad_phi_u[j] * normal_vector) - 
                 stokes_phi_p[j] * normal_vector) * 
                elasticity_phi[i] * stokes_fe_face_values.JxW(q)); 
      } 
  } 
// @sect4{<code>FluidStructureProblem::solve</code>}  

// 正如介绍中所讨论的，我们在这里使用了一个相当琐碎的求解器：我们只是将线性系统传递给SparseDirectUMFPACK直接求解器（例如，见 step-29  ）。在求解之后，我们唯一要做的是确保悬挂的节点和边界值约束是正确的。

  template <int dim> 
  void FluidStructureProblem<dim>::solve() 
  { 
    SparseDirectUMFPACK direct_solver; 
    direct_solver.initialize(system_matrix); 
    direct_solver.vmult(solution, system_rhs); 

    constraints.distribute(solution); 
  } 

//  @sect4{<code>FluidStructureProblem::output_results</code>}  

// 生成图形输出在这里相当简单：我们所要做的就是确定解向量的哪些成分属于标量和/或向量（例如，见 step-22 之前的例子），然后把它全部传递给DataOut类。

  template <int dim> 
  void FluidStructureProblem<dim>::output_results( 
    const unsigned int refinement_cycle) const 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("pressure"); 
    for (unsigned int d = 0; d < dim; ++d) 
      solution_names.emplace_back("displacement"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    for (unsigned int d = 0; d < dim; ++d) 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_part_of_vector); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    data_out.add_data_vector(solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output( 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk"); 
    data_out.write_vtk(output); 
  } 
// @sect4{<code>FluidStructureProblem::refine_mesh</code>}  

// 下一步是细化网格。正如在介绍中所讨论的，这有点棘手，主要是因为流体和固体子域使用的变量具有不同的物理尺寸，因此，误差估计的绝对大小不能直接比较。因此，我们将不得不对它们进行缩放。因此，在函数的顶部，我们首先分别计算不同变量的误差估计值（在流体域中使用速度而不是压力，在固体域中使用位移）。

  template <int dim> 
  void FluidStructureProblem<dim>::refine_mesh() 
  { 
    Vector<float> stokes_estimated_error_per_cell( 
      triangulation.n_active_cells()); 
    Vector<float> elasticity_estimated_error_per_cell( 
      triangulation.n_active_cells()); 

    const QGauss<dim - 1> stokes_face_quadrature(stokes_degree + 2); 
    const QGauss<dim - 1> elasticity_face_quadrature(elasticity_degree + 2); 

    hp::QCollection<dim - 1> face_q_collection; 
    face_q_collection.push_back(stokes_face_quadrature); 
    face_q_collection.push_back(elasticity_face_quadrature); 

    const FEValuesExtractors::Vector velocities(0); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_q_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      stokes_estimated_error_per_cell, 
      fe_collection.component_mask(velocities)); 

    const FEValuesExtractors::Vector displacements(dim + 1); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      face_q_collection, 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      elasticity_estimated_error_per_cell, 
      fe_collection.component_mask(displacements)); 

// 然后，我们通过除以误差估计值的法线对其进行归一化处理，并按照介绍中所讨论的那样，将流体误差指标按4的系数进行缩放。然后将这些结果加在一起，形成一个包含所有单元的误差指标的向量。

    stokes_estimated_error_per_cell *= 
      4. / stokes_estimated_error_per_cell.l2_norm(); 
    elasticity_estimated_error_per_cell *= 
      1. / elasticity_estimated_error_per_cell.l2_norm(); 

    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    estimated_error_per_cell += stokes_estimated_error_per_cell; 
    estimated_error_per_cell += elasticity_estimated_error_per_cell; 

// 在实际细化网格之前，函数的倒数第二部分涉及到我们在介绍中已经提到的启发式方法：由于解是不连续的，KellyErrorEstimator类对位于子域之间边界的单元感到困惑：它认为那里的误差很大，因为梯度的跳跃很大，尽管这完全是预期的，事实上在精确解中也存在这一特征，因此不表明任何数值错误。

// 因此，我们将界面上的所有单元的误差指标设置为零；决定影响哪些单元的条件略显尴尬，因为我们必须考虑到自适应细化网格的可能性，这意味着邻近的单元可能比当前的单元更粗，或者事实上可能被细化一些。这些嵌套条件的结构与我们在 <code>assemble_system</code> 中组装接口条款时遇到的情况基本相同。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto f : cell->face_indices()) 
        if (cell_is_in_solid_domain(cell)) 
          { 
            if ((cell->at_boundary(f) == false) && 
                (((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == false) && 
                  cell_is_in_fluid_domain(cell->neighbor(f))) || 
                 ((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == true) && 
                  (cell_is_in_fluid_domain( 
                    cell->neighbor_child_on_subface(f, 0)))) || 
                 (cell->neighbor_is_coarser(f) && 
                  cell_is_in_fluid_domain(cell->neighbor(f))))) 
              estimated_error_per_cell(cell->active_cell_index()) = 0; 
          } 
        else 
          { 
            if ((cell->at_boundary(f) == false) && 
                (((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == false) && 
                  cell_is_in_solid_domain(cell->neighbor(f))) || 
                 ((cell->neighbor(f)->level() == cell->level()) && 
                  (cell->neighbor(f)->has_children() == true) && 
                  (cell_is_in_solid_domain( 
                    cell->neighbor_child_on_subface(f, 0)))) || 
                 (cell->neighbor_is_coarser(f) && 
                  cell_is_in_solid_domain(cell->neighbor(f))))) 
              estimated_error_per_cell(cell->active_cell_index()) = 0; 
          } 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.0); 
    triangulation.execute_coarsening_and_refinement(); 
  } 

//  @sect4{<code>FluidStructureProblem::run</code>}  

// 像往常一样，这是控制整个操作流程的函数。如果你读过教程程序  step-1  到  step-6  ，例如，那么你已经对以下结构相当熟悉。

  template <int dim> 
  void FluidStructureProblem<dim>::run() 
  { 
    make_grid(); 

    for (unsigned int refinement_cycle = 0; refinement_cycle < 10 - 2 * dim; 
         ++refinement_cycle) 
      { 
        std::cout << "Refinement cycle " << refinement_cycle << std::endl; 

        if (refinement_cycle > 0) 
          refine_mesh(); 

        setup_dofs(); 

        std::cout << "   Assembling..." << std::endl; 
        assemble_system(); 

        std::cout << "   Solving..." << std::endl; 
        solve(); 

        std::cout << "   Writing output..." << std::endl; 
        output_results(refinement_cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step46 

//  @sect4{The <code>main()</code> function}  

// 这个，最后的，函数所包含的内容几乎与其他大多数教程程序的内容完全一样。

int main() 
{ 
  try 
    { 
      using namespace Step46; 

      FluidStructureProblem<2> flow_problem(1, 1); 
      flow_problem.run(); 
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


