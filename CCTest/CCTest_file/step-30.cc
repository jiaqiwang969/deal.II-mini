

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2007 - 2021 by the deal.II authors 
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
 * Author: Tobias Leicht, 2007 
 */ 



// deal.II包括的文件已经在前面的例子中介绍过了，因此不再做进一步的评论。

#include <deal.II/base/function.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/lac/precondition_block.h> 
#include <deal.II/lac/solver_richardson.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q1.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/derivative_approximation.h> 

// 而这又是C++。

#include <array> 
#include <iostream> 
#include <fstream> 

// 最后一步和以前所有的程序一样。

namespace Step30 
{ 
  using namespace dealii; 
// @sect3{Equation data}  

// 描述方程数据的类和单个术语的实际装配几乎完全照搬自  step-12  。我们将对差异进行评论。

  template <int dim> 
  class RHS : public Function<dim> 
  { 
  public: 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int /*component*/ = 0) const override 
    { 
      (void)points; 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      std::fill(values.begin(), values.end(), 0.); 
    } 
  }; 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int /*component*/ = 0) const override 
    { 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      for (unsigned int i = 0; i < values.size(); ++i) 
        { 
          if (points[i](0) < 0.5) 
            values[i] = 1.; 
          else 
            values[i] = 0.; 
        } 
    } 
  }; 

  template <int dim> 
  class Beta 
  { 
  public: 

//流场选择为逆时针方向的四分之一圆，原点为域的右半部分的中点，数值为正 $x$ ，而在域的左边部分，流速只是向左走，与从右边进来的流速一致。在圆形部分，流速的大小与离原点的距离成正比。这与 step-12 不同，在该定义中，到处都是1。新定义导致 $\beta$ 沿单元的每个给定面的线性变化。另一方面， $u(x,y)$ 的解决方案与之前完全相同。

    void value_list(const std::vector<Point<dim>> &points, 
                    std::vector<Point<dim>> &      values) const 
    { 
      Assert(values.size() == points.size(), 
             ExcDimensionMismatch(values.size(), points.size())); 

      for (unsigned int i = 0; i < points.size(); ++i) 
        { 
          if (points[i](0) > 0) 
            { 
              values[i](0) = -points[i](1); 
              values[i](1) = points[i](0); 
            } 
          else 
            { 
              values[i]    = Point<dim>(); 
              values[i](0) = -points[i](1); 
            } 
        } 
    } 
  }; 

//  @sect3{Class: DGTransportEquation}  

// 这个类的声明完全不受我们目前的变化影响。

  template <int dim> 
  class DGTransportEquation 
  { 
  public: 
    DGTransportEquation(); 

    void assemble_cell_term(const FEValues<dim> &fe_v, 
                            FullMatrix<double> & ui_vi_matrix, 
                            Vector<double> &     cell_vector) const; 

    void assemble_boundary_term(const FEFaceValues<dim> &fe_v, 
                                FullMatrix<double> &     ui_vi_matrix, 
                                Vector<double> &         cell_vector) const; 

    void assemble_face_term(const FEFaceValuesBase<dim> &fe_v, 
                            const FEFaceValuesBase<dim> &fe_v_neighbor, 
                            FullMatrix<double> &         ui_vi_matrix, 
                            FullMatrix<double> &         ue_vi_matrix, 
                            FullMatrix<double> &         ui_ve_matrix, 
                            FullMatrix<double> &         ue_ve_matrix) const; 

  private: 
    const Beta<dim>           beta_function; 
    const RHS<dim>            rhs_function; 
    const BoundaryValues<dim> boundary_function; 
  }; 

// 同样地，该类的构造函数以及组装对应于单元格内部和边界面的术语的函数与之前没有变化。装配单元间面术语的函数也没有改变，因为它所做的只是对两个FEFaceValuesBase类型的对象进行操作（它是FEFaceValues和FESubfaceValues的基类）。这些对象从何而来，即它们是如何被初始化的，对这个函数来说并不重要：它只是假设这两个对象所代表的面或子面上的正交点与物理空间中的相同点相对应。

  template <int dim> 
  DGTransportEquation<dim>::DGTransportEquation() 
    : beta_function() 
    , rhs_function() 
    , boundary_function() 
  {} 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_cell_term( 
    const FEValues<dim> &fe_v, 
    FullMatrix<double> & ui_vi_matrix, 
    Vector<double> &     cell_vector) const 
  { 
    const std::vector<double> &JxW = fe_v.get_JxW_values(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 
    std::vector<double>     rhs(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 
    rhs_function.value_list(fe_v.get_quadrature_points(), rhs); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
        { 
          for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
            ui_vi_matrix(i, j) -= beta[point] * fe_v.shape_grad(i, point) * 
                                  fe_v.shape_value(j, point) * JxW[point]; 

          cell_vector(i) += 
            rhs[point] * fe_v.shape_value(i, point) * JxW[point]; 
        } 
  } 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_boundary_term( 
    const FEFaceValues<dim> &fe_v, 
    FullMatrix<double> &     ui_vi_matrix, 
    Vector<double> &         cell_vector) const 
  { 
    const std::vector<double> &        JxW     = fe_v.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 
    std::vector<double>     g(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 
    boundary_function.value_list(fe_v.get_quadrature_points(), g); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      { 
        const double beta_n = beta[point] * normals[point]; 
        if (beta_n > 0) 
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
              ui_vi_matrix(i, j) += beta_n * fe_v.shape_value(j, point) * 
                                    fe_v.shape_value(i, point) * JxW[point]; 
        else 
          for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
            cell_vector(i) -= 
              beta_n * g[point] * fe_v.shape_value(i, point) * JxW[point]; 
      } 
  } 

  template <int dim> 
  void DGTransportEquation<dim>::assemble_face_term( 
    const FEFaceValuesBase<dim> &fe_v, 
    const FEFaceValuesBase<dim> &fe_v_neighbor, 
    FullMatrix<double> &         ui_vi_matrix, 
    FullMatrix<double> &         ue_vi_matrix, 
    FullMatrix<double> &         ui_ve_matrix, 
    FullMatrix<double> &         ue_ve_matrix) const 
  { 
    const std::vector<double> &        JxW     = fe_v.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

    std::vector<Point<dim>> beta(fe_v.n_quadrature_points); 

    beta_function.value_list(fe_v.get_quadrature_points(), beta); 

    for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
      { 
        const double beta_n = beta[point] * normals[point]; 
        if (beta_n > 0) 
          { 
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
                ui_vi_matrix(i, j) += beta_n * fe_v.shape_value(j, point) * 
                                      fe_v.shape_value(i, point) * JxW[point]; 

            for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k) 
              for (unsigned int j = 0; j < fe_v.dofs_per_cell; ++j) 
                ui_ve_matrix(k, j) -= beta_n * fe_v.shape_value(j, point) * 
                                      fe_v_neighbor.shape_value(k, point) * 
                                      JxW[point]; 
          } 
        else 
          { 
            for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
              for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l) 
                ue_vi_matrix(i, l) += beta_n * 
                                      fe_v_neighbor.shape_value(l, point) * 
                                      fe_v.shape_value(i, point) * JxW[point]; 

            for (unsigned int k = 0; k < fe_v_neighbor.dofs_per_cell; ++k) 
              for (unsigned int l = 0; l < fe_v_neighbor.dofs_per_cell; ++l) 
                ue_ve_matrix(k, l) -= 
                  beta_n * fe_v_neighbor.shape_value(l, point) * 
                  fe_v_neighbor.shape_value(k, point) * JxW[point]; 
          } 
      } 
  } 
// @sect3{Class: DGMethod}  

// 这个声明很像  step-12  的声明。然而，我们引入了一个新的例程（set_anisotropic_flags）并修改了另一个例程（refine_grid）。

  template <int dim> 
  class DGMethod 
  { 
  public: 
    DGMethod(const bool anisotropic); 

    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(Vector<double> &solution); 
    void refine_grid(); 
    void set_anisotropic_flags(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 

// 我们再次希望使用程度为1的DG元素（但这只在构造函数中指定）。如果你想使用不同程度的DG方法，请在构造函数中用新的程度替换1。

    const unsigned int degree; 
    FE_DGQ<dim>        fe; 
    DoFHandler<dim>    dof_handler; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

// 这是新的，在介绍中解释的各向异性跳跃指标的评估中使用的阈值。它的值在构造函数中被设置为3.0，但它可以很容易地被改变为一个大于1的不同值。

    const double anisotropic_threshold_ratio; 

// 这是一个指示是否使用各向异性细化的bool标志。它由构造函数设置，构造函数需要一个同名的参数。

    const bool anisotropic; 

    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 

    Vector<double> solution2; 
    Vector<double> right_hand_side; 

    const DGTransportEquation<dim> dg; 
  }; 

  template <int dim> 
  DGMethod<dim>::DGMethod(const bool anisotropic) 
    : mapping() 
    , 

// 对于不同程度的DG方法，在这里进行修改。

    degree(1) 
    , fe(degree) 
    , dof_handler(triangulation) 
    , anisotropic_threshold_ratio(3.) 
    , anisotropic(anisotropic) 
    , 

// 由于β是一个线性函数，我们可以选择正交的度数，对于这个度数，所得的积分是正确的。因此，我们选择使用  <code>degree+1</code>  高斯点，这使我们能够准确地积分度数为  <code>2*degree+1</code>  的多项式，足以满足我们在本程序中要进行的所有积分。

    quadrature(degree + 1) 
    , face_quadrature(degree + 1) 
    , dg() 
  {} 

  template <int dim> 
  void DGMethod<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    sparsity_pattern.reinit(dof_handler.n_dofs(), 
                            dof_handler.n_dofs(), 
                            (GeometryInfo<dim>::faces_per_cell * 
                               GeometryInfo<dim>::max_children_per_face + 
                             1) * 
                              fe.n_dofs_per_cell()); 

    DoFTools::make_flux_sparsity_pattern(dof_handler, sparsity_pattern); 

    sparsity_pattern.compress(); 

    system_matrix.reinit(sparsity_pattern); 

    solution2.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{Function: assemble_system}  

// 我们继续使用 <code>assemble_system</code> 函数来实现DG离散化。这个函数与 step-12 中的 <code>assemble_system</code> 函数的作用相同（但没有MeshWorker）。 一个单元的邻居关系所考虑的四种情况与各向同性的情况相同，即a)单元在边界上，b)有更细的邻居单元，c)邻居既不粗也不细，d)邻居更粗。 然而，我们决定哪种情况的方式是按照介绍中描述的方式进行修改的。

  template <int dim> 
  void DGMethod<dim>::assemble_system() 
  { 
    const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
    std::vector<types::global_dof_index> dofs(dofs_per_cell); 
    std::vector<types::global_dof_index> dofs_neighbor(dofs_per_cell); 

    const UpdateFlags update_flags = update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values; 

    const UpdateFlags face_update_flags = 
      update_values | update_quadrature_points | update_JxW_values | 
      update_normal_vectors; 

    const UpdateFlags neighbor_face_update_flags = update_values; 

    FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags); 
    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         neighbor_face_update_flags); 

    FullMatrix<double> ui_vi_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> ue_vi_matrix(dofs_per_cell, dofs_per_cell); 

    FullMatrix<double> ui_ve_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> ue_ve_matrix(dofs_per_cell, dofs_per_cell); 

    Vector<double> cell_vector(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        ui_vi_matrix = 0; 
        cell_vector  = 0; 

        fe_v.reinit(cell); 

        dg.assemble_cell_term(fe_v, ui_vi_matrix, cell_vector); 

        cell->get_dof_indices(dofs); 

        for (const auto face_no : cell->face_indices()) 
          { 
            const auto face = cell->face(face_no); 

// 情况(a)。该面在边界上。

            if (face->at_boundary()) 
              { 
                fe_v_face.reinit(cell, face_no); 

                dg.assemble_boundary_term(fe_v_face, ui_vi_matrix, cell_vector); 
              } 
            else 
              { 
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, 
                       ExcInternalError()); 
                const auto neighbor = cell->neighbor(face_no); 

// 情况(b)。这是一个内部面，邻居是精炼的（我们可以通过询问当前单元格的面是否有孩子来测试）。在这种情况下，我们需要对 "子面 "进行整合，即当前单元格的面的子女。            (有一个稍微令人困惑的角落案例。如果我们是在1d中--诚然，当前的程序和它对各向异性细化的演示并不特别相关--那么单元间的面总是相同的：它们只是顶点。换句话说，在1d中，我们不希望对不同层次的单元之间的面进行不同的处理。我们在这里检查的条件`face->has_children()`确保了这一点：在1d中，这个函数总是返回`false`，因此在1d中我们不会进入这个`if`分支。但我们将不得不在下面的情况（c）中回到这个角落。

                if (face->has_children()) 
                  { 

// 我们需要知道，哪个邻居的面朝向我们单元格的方向。使用  @p  neighbor_face_no 函数，我们可以得到粗邻和非粗邻的这些信息。

                    const unsigned int neighbor2 = 
                      cell->neighbor_face_no(face_no); 

// 现在我们对所有的子脸进行循环，也就是当前脸的子脸和可能的孙子脸。

                    for (unsigned int subface_no = 0; 
                         subface_no < face->n_active_descendants(); 
                         ++subface_no) 
                      { 

// 为了得到当前子面后面的单元，我们可以使用 @p neighbor_child_on_subface 函数。它照顾到了所有各向异性细化和非标准面的复杂情况。

                        const auto neighbor_child = 
                          cell->neighbor_child_on_subface(face_no, subface_no); 
                        Assert(!neighbor_child->has_children(), 
                               ExcInternalError()); 

// 这个案例的其余部分没有变化。

                        ue_vi_matrix = 0; 
                        ui_ve_matrix = 0; 
                        ue_ve_matrix = 0; 

                        fe_v_subface.reinit(cell, face_no, subface_no); 
                        fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 

                        dg.assemble_face_term(fe_v_subface, 
                                              fe_v_face_neighbor, 
                                              ui_vi_matrix, 
                                              ue_vi_matrix, 
                                              ui_ve_matrix, 
                                              ue_ve_matrix); 

                        neighbor_child->get_dof_indices(dofs_neighbor); 

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                            { 
                              system_matrix.add(dofs[i], 
                                                dofs_neighbor[j], 
                                                ue_vi_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs[j], 
                                                ui_ve_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs_neighbor[j], 
                                                ue_ve_matrix(i, j)); 
                            } 
                      } 
                  } 
                else 
                  { 

//情况(c)。如果这是一个内部面，并且邻居没有进一步细化，我们就会得到这里（或者，如上所述，我们是在1d中，在这种情况下，我们对每个内部面都会得到这里）。然后我们需要决定是否要对当前面进行整合。如果邻居实际上更粗，那么我们就忽略这个面，而是在访问邻居单元并查看当前面的时候处理它（除了在1d中，如上所述，这不会发生）。

                    if (dim > 1 && cell->neighbor_is_coarser(face_no)) 
                      continue; 

// 另一方面，如果邻居是更精细的，那么我们已经处理了上面(b)情况下的脸（1d除外）。所以对于2d和3d，我们只需要决定是要处理来自当前一侧的同一层次的单元格之间的面，还是来自邻接一侧的面。 我们通过引入一个平局来做到这一点。          我们只取索引较小的单元格（在当前细化级别内）。在1d中，我们取较粗的单元，或者如果它们在同一层次，则取该层次中指数较小的单元。这就导致了一个复杂的条件，希望在上面的描述中可以理解。

                    if (((dim > 1) && (cell->index() < neighbor->index())) || 
                        ((dim == 1) && ((cell->level() < neighbor->level()) || 
                                        ((cell->level() == neighbor->level()) && 
                                         (cell->index() < neighbor->index()))))) 
                      { 

// 这里我们知道，邻居不是更粗的，所以我们可以使用通常的  @p neighbor_of_neighbor  函数。然而，我们也可以使用更通用的 @p neighbor_face_no 函数。

                        const unsigned int neighbor2 = 
                          cell->neighbor_of_neighbor(face_no); 

                        ue_vi_matrix = 0; 
                        ui_ve_matrix = 0; 
                        ue_ve_matrix = 0; 

                        fe_v_face.reinit(cell, face_no); 
                        fe_v_face_neighbor.reinit(neighbor, neighbor2); 

                        dg.assemble_face_term(fe_v_face, 
                                              fe_v_face_neighbor, 
                                              ui_vi_matrix, 
                                              ue_vi_matrix, 
                                              ui_ve_matrix, 
                                              ue_ve_matrix); 

                        neighbor->get_dof_indices(dofs_neighbor); 

                        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                            { 
                              system_matrix.add(dofs[i], 
                                                dofs_neighbor[j], 
                                                ue_vi_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs[j], 
                                                ui_ve_matrix(i, j)); 
                              system_matrix.add(dofs_neighbor[i], 
                                                dofs_neighbor[j], 
                                                ue_ve_matrix(i, j)); 
                            } 
                      } 

// 我们不需要考虑情况(d)，因为这些面在情况(b)中被 "从另一侧 "处理。

                  } 
              } 
          } 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(dofs[i], dofs[j], ui_vi_matrix(i, j)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          right_hand_side(dofs[i]) += cell_vector(i); 
      } 
  } 
// @sect3{Solver}  

// 对于这个简单的问题，我们再次使用简单的Richardson迭代法。该求解器完全不受我们各向异性变化的影响。

  template <int dim> 
  void DGMethod<dim>::solve(Vector<double> &solution) 
  { 
    SolverControl                    solver_control(1000, 1e-12, false, false); 
    SolverRichardson<Vector<double>> solver(solver_control); 

    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 

    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 

    solver.solve(system_matrix, solution, right_hand_side, preconditioner); 
  } 
// @sect3{Refinement}  

// 我们根据 step-12 中使用的相同的简单细化标准来细化网格，即对解的梯度的近似。

  template <int dim> 
  void DGMethod<dim>::refine_grid() 
  { 
    Vector<float> gradient_indicator(triangulation.n_active_cells()); 

// 我们对梯度进行近似计算。

    DerivativeApproximation::approximate_gradient(mapping, 
                                                  dof_handler, 
                                                  solution2, 
                                                  gradient_indicator); 

//并对其进行缩放，以获得一个误差指标。

    for (const auto &cell : triangulation.active_cell_iterators()) 
      gradient_indicator[cell->active_cell_index()] *= 
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2); 

// 然后我们用这个指标来标记误差指标最高的30%的单元格来进行精炼。

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    gradient_indicator, 
                                                    0.3, 
                                                    0.1); 

// 现在，细化标志被设置为那些具有大误差指标的单元。如果不做任何改变，这些单元将被等向细化。如果给这个函数的 @p anisotropic 标志被设置，我们现在调用set_anisotropic_flags()函数，该函数使用跳转指标将一些细化标志重置为各向异性细化。

    if (anisotropic) 
      set_anisotropic_flags(); 

// 现在执行考虑各向异性以及各向同性的细化标志的细化。

    triangulation.execute_coarsening_and_refinement(); 
  } 

// 一旦错误指标被评估，误差最大的单元被标记为细化，我们要再次循环这些被标记的单元，以决定它们是否需要各向同性的细化或各向异性的细化更为合适。这就是在介绍中解释的各向异性跳跃指标。

  template <int dim> 
  void DGMethod<dim>::set_anisotropic_flags() 
  { 

// 我们想在被标记的单元格的面上评估跳跃，所以我们需要一些对象来评估面上的解决方案的值。

    UpdateFlags face_update_flags = 
      UpdateFlags(update_values | update_JxW_values); 

    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         update_values); 

// 现在我们需要对所有活动单元进行循环。

    for (const auto &cell : dof_handler.active_cell_iterators()) 

// 我们只需要考虑那些被标记为细化的单元。

      if (cell->refine_flag_set()) 
        { 
          Point<dim> jump; 
          Point<dim> area; 

          for (const auto face_no : cell->face_indices()) 
            { 
              const auto face = cell->face(face_no); 

              if (!face->at_boundary()) 
                { 
                  Assert(cell->neighbor(face_no).state() == 
                           IteratorState::valid, 
                         ExcInternalError()); 
                  const auto neighbor = cell->neighbor(face_no); 

                  std::vector<double> u(fe_v_face.n_quadrature_points); 
                  std::vector<double> u_neighbor(fe_v_face.n_quadrature_points); 

// 在汇编例程中看到的四种不同的邻居关系的情况在这里以同样的方式重复。

                  if (face->has_children()) 
                    { 

// 邻居被完善。 首先，我们存储信息，即邻居的哪个面指向我们当前单元的方向。这个属性将被继承给子代。

                      unsigned int neighbor2 = cell->neighbor_face_no(face_no); 

// 现在我们对所有的子面进行循环。

                      for (unsigned int subface_no = 0; 
                           subface_no < face->n_active_descendants(); 
                           ++subface_no) 
                        { 

//得到一个迭代器，指向当前子面后面的单元格...

                          const auto neighbor_child = 
                            cell->neighbor_child_on_subface(face_no, 
                                                            subface_no); 
                          Assert(!neighbor_child->has_children(), 
                                 ExcInternalError()); 

// ...并重新启动各自的FEFaceValues和FESSubFaceValues对象。

                          fe_v_subface.reinit(cell, face_no, subface_no); 
                          fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 

// 我们获得了函数值

                          fe_v_subface.get_function_values(solution2, u); 
                          fe_v_face_neighbor.get_function_values(solution2, 
                                                                 u_neighbor); 

//以及正交权重，乘以雅各布行列式。

                          const std::vector<double> &JxW = 
                            fe_v_subface.get_JxW_values(); 

// 现在我们在所有的正交点上循环。

                          for (unsigned int x = 0; 
                               x < fe_v_subface.n_quadrature_points; 
                               ++x) 
                            { 

//并整合解决方案的跳跃的绝对值，即分别从当前单元和邻近单元看到的函数值的绝对值。我们知道，前两个面与单元格上的第一个坐标方向正交，后两个面与第二个坐标方向正交，以此类推，所以我们将这些值累积成具有 <code>dim</code> 成分的向量。

                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 

// 我们还将缩放后的权重相加，以获得脸部的量度。

                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                    } 
                  else 
                    { 
                      if (!cell->neighbor_is_coarser(face_no)) 
                        { 

// 我们的当前单元和邻居在所考虑的面有相同的细化。除此以外，我们的做法与上述情况下的一个子单元的做法基本相同。

                          unsigned int neighbor2 = 
                            cell->neighbor_of_neighbor(face_no); 

                          fe_v_face.reinit(cell, face_no); 
                          fe_v_face_neighbor.reinit(neighbor, neighbor2); 

                          fe_v_face.get_function_values(solution2, u); 
                          fe_v_face_neighbor.get_function_values(solution2, 
                                                                 u_neighbor); 

                          const std::vector<double> &JxW = 
                            fe_v_face.get_JxW_values(); 

                          for (unsigned int x = 0; 
                               x < fe_v_face.n_quadrature_points; 
                               ++x) 
                            { 
                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 
                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                      else // i.e. neighbor is coarser than cell 
                        { 

// 现在邻居实际上更粗了。这种情况是新的，因为它没有出现在汇编程序中。在这里，我们必须考虑它，但这并不太复杂。我们只需使用  @p  neighbor_of_coarser_neighbor 函数，它再次自行处理各向异性的细化和非标准面的方向。

                          std::pair<unsigned int, unsigned int> 
                            neighbor_face_subface = 
                              cell->neighbor_of_coarser_neighbor(face_no); 
                          Assert(neighbor_face_subface.first < cell->n_faces(), 
                                 ExcInternalError()); 
                          Assert(neighbor_face_subface.second < 
                                   neighbor->face(neighbor_face_subface.first) 
                                     ->n_active_descendants(), 
                                 ExcInternalError()); 
                          Assert(neighbor->neighbor_child_on_subface( 
                                   neighbor_face_subface.first, 
                                   neighbor_face_subface.second) == cell, 
                                 ExcInternalError()); 

                          fe_v_face.reinit(cell, face_no); 
                          fe_v_subface.reinit(neighbor, 
                                              neighbor_face_subface.first, 
                                              neighbor_face_subface.second); 

                          fe_v_face.get_function_values(solution2, u); 
                          fe_v_subface.get_function_values(solution2, 
                                                           u_neighbor); 

                          const std::vector<double> &JxW = 
                            fe_v_face.get_JxW_values(); 

                          for (unsigned int x = 0; 
                               x < fe_v_face.n_quadrature_points; 
                               ++x) 
                            { 
                              jump[face_no / 2] += 
                                std::abs(u[x] - u_neighbor[x]) * JxW[x]; 
                              area[face_no / 2] += JxW[x]; 
                            } 
                        } 
                    } 
                } 
            } 

// 现在我们分析一下平均跳动的大小，我们用跳动除以各面的度量得到。

          std::array<double, dim> average_jumps; 
          double                  sum_of_average_jumps = 0.; 
          for (unsigned int i = 0; i < dim; ++i) 
            { 
              average_jumps[i] = jump(i) / area(i); 
              sum_of_average_jumps += average_jumps[i]; 
            } 

// 现在我们在单元格的 <code>dim</code> 坐标方向上进行循环，并比较与该方向正交的面的平均跳跃和与其余方向正交的面的平均跳跃。如果前者比后者大一个给定的系数，我们只沿帽轴进行细化。否则，我们不改变细化标志，导致各向同性的细化。

          for (unsigned int i = 0; i < dim; ++i) 
            if (average_jumps[i] > anisotropic_threshold_ratio * 
                                     (sum_of_average_jumps - average_jumps[i])) 
              cell->set_refine_flag(RefinementCase<dim>::cut_axis(i)); 
        } 
  } 
// @sect3{The Rest}  

// 程序的其余部分非常遵循之前教程程序的方案。我们以VTU格式输出网格（就像我们在 step-1 中所做的那样，例如），并以VTU格式输出可视化，我们几乎总是这样做。

  template <int dim> 
  void DGMethod<dim>::output_results(const unsigned int cycle) const 
  { 
    std::string refine_type; 
    if (anisotropic) 
      refine_type = ".aniso"; 
    else 
      refine_type = ".iso"; 

    { 
      const std::string filename = 
        "grid-" + std::to_string(cycle) + refine_type + ".svg"; 
      std::cout << "   Writing grid to <" << filename << ">..." << std::endl; 
      std::ofstream svg_output(filename); 

      GridOut grid_out; 
      grid_out.write_svg(triangulation, svg_output); 
    } 

    { 
      const std::string filename = 
        "sol-" + std::to_string(cycle) + refine_type + ".vtu"; 
      std::cout << "   Writing solution to <" << filename << ">..." 
                << std::endl; 
      std::ofstream gnuplot_output(filename); 

      DataOut<dim> data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution2, "u"); 

      data_out.build_patches(degree); 

      data_out.write_vtu(gnuplot_output); 
    } 
  } 

  template <int dim> 
  void DGMethod<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          { 

// 创建矩形域。

            Point<dim> p1, p2; 
            p1(0) = 0; 
            p1(0) = -1; 
            for (unsigned int i = 0; i < dim; ++i) 
              p2(i) = 1.; 

// 调整不同方向的单元数，以获得原始网格的完全各向同性的单元。

            std::vector<unsigned int> repetitions(dim, 1); 
            repetitions[0] = 2; 
            GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                                      repetitions, 
                                                      p1, 
                                                      p2); 

            triangulation.refine_global(5 - dim); 
          } 
        else 
          refine_grid(); 

        std::cout << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl; 

        Timer assemble_timer; 
        assemble_system(); 
        std::cout << "   Time of assemble_system: " << assemble_timer.cpu_time() 
                  << std::endl; 
        solve(solution2); 

        output_results(cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step30 

int main() 
{ 
  try 
    { 
      using namespace Step30; 

// 如果你想以3D方式运行程序，只需将下面一行改为 <code>const unsigned int dim = 3;</code>  。

      const unsigned int dim = 2; 

      { 

// 首先，我们用各向同性的细化方法进行一次运行。

        std::cout << "Performing a " << dim 
                  << "D run with isotropic refinement..." << std::endl 
                  << "------------------------------------------------" 
                  << std::endl; 
        DGMethod<dim> dgmethod_iso(false); 
        dgmethod_iso.run(); 
      } 

      { 

// 现在我们进行第二次运行，这次是各向异性的细化。

        std::cout << std::endl 
                  << "Performing a " << dim 
                  << "D run with anisotropic refinement..." << std::endl 
                  << "--------------------------------------------------" 
                  << std::endl; 
        DGMethod<dim> dgmethod_aniso(true); 
        dgmethod_aniso.run(); 
      } 
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
    }; 

  return 0; 
} 

