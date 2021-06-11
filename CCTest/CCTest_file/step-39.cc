CCTest_file/step-39.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2010 - 2020 by the deal.II authors 
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
 * Author: Guido Kanschat, Texas A&M University, 2009 
 */ 



// 线性代数的包含文件。一个普通的SparseMatrix，它又将包括SparsityPattern和Vector类的必要文件。

#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/precondition_block.h> 
#include <deal.II/lac/block_vector.h> 

// 包括用于设置网格的文件

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

// FiniteElement类和DoFHandler的包含文件。

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_dgp.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/dofs/dof_tools.h> 

// 使用MeshWorker框架的包含文件

#include <deal.II/meshworker/dof_info.h> 
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/assembler.h> 
#include <deal.II/meshworker/loop.h> 

// 与拉普拉斯相关的局部积分器的包含文件

#include <deal.II/integrators/laplace.h> 

// 支持多网格方法

#include <deal.II/multigrid/mg_tools.h> 
#include <deal.II/multigrid/multigrid.h> 
#include <deal.II/multigrid/mg_matrix.h> 
#include <deal.II/multigrid/mg_transfer.h> 
#include <deal.II/multigrid/mg_coarse.h> 
#include <deal.II/multigrid/mg_smoother.h> 

// 最后，我们从库中取出我们的精确解，以及正交和附加工具。

#include <deal.II/base/function_lib.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <iostream> 
#include <fstream> 

// deal.II库的所有类都在dealii命名空间中。为了节省打字，我们告诉编译器也要在其中搜索名字。

namespace Step39 
{ 
  using namespace dealii; 

// 这是我们用来设置边界值的函数，也是我们比较的精确解。

  Functions::SlitSingularityFunction<2> exact_solution; 
// @sect3{The local integrators}  

// MeshWorker将局部积分与单元格和面的循环分离开来。因此，我们必须编写局部积分类来生成矩阵、右手边和误差估计器。

// 所有这些类都有相同的三个函数，分别用于对单元、边界面和内部面的积分。局部积分所需的所有信息都由 MeshWorker::IntegrationInfo<dim>. 提供。请注意，函数的签名不能改变，因为它是由 MeshWorker::integration_loop(). 所期望的。

// 第一个定义局部积分器的类负责计算单元和面矩阵。它被用来组装全局矩阵以及水平矩阵。

  template <int dim> 
  class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
              typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void 
         boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
                  typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
              MeshWorker::DoFInfo<dim> &                 dinfo2, 
              typename MeshWorker::IntegrationInfo<dim> &info1, 
              typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
  }; 

// 在每个单元上，我们对Dirichlet形式进行积分。我们使用LocalIntegrators中的现成积分库来避免自己编写这些循环。同样地，我们实现了Nitsche边界条件和单元间的内部惩罚通量。

// 边界和通量项需要一个惩罚参数，这个参数应该根据单元的大小和多项式的度数来调整。在 LocalIntegrators::Laplace::compute_penalty() 中可以找到关于这个参数的安全选择，我们在下面使用这个参数。

  template <int dim> 
  void MatrixIntegrator<dim>::cell( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix, 
                                           info.fe_values()); 
  } 

  template <int dim> 
  void MatrixIntegrator<dim>::boundary( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const unsigned int degree = info.fe_values(0).get_fe().tensor_degree(); 
    LocalIntegrators::Laplace::nitsche_matrix( 
      dinfo.matrix(0, false).matrix, 
      info.fe_values(0), 
      LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, degree, degree)); 
  } 

// 内部面使用内部惩罚方法

  template <int dim> 
  void MatrixIntegrator<dim>::face( 
    MeshWorker::DoFInfo<dim> &                 dinfo1, 
    MeshWorker::DoFInfo<dim> &                 dinfo2, 
    typename MeshWorker::IntegrationInfo<dim> &info1, 
    typename MeshWorker::IntegrationInfo<dim> &info2) const 
  { 
    const unsigned int degree = info1.fe_values(0).get_fe().tensor_degree(); 
    LocalIntegrators::Laplace::ip_matrix( 
      dinfo1.matrix(0, false).matrix, 
      dinfo1.matrix(0, true).matrix, 
      dinfo2.matrix(0, true).matrix, 
      dinfo2.matrix(0, false).matrix, 
      info1.fe_values(0), 
      info2.fe_values(0), 
      LocalIntegrators::Laplace::compute_penalty( 
        dinfo1, dinfo2, degree, degree)); 
  } 

// 第二个局部积分器建立了右手边。在我们的例子中，右手边的函数为零，这样，这里只设置了弱形式的边界条件。

  template <int dim> 
  class RHSIntegrator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
              typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void 
         boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
                  typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
              MeshWorker::DoFInfo<dim> &                 dinfo2, 
              typename MeshWorker::IntegrationInfo<dim> &info1, 
              typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
  }; 

  template <int dim> 
  void 
  RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &, 
                           typename MeshWorker::IntegrationInfo<dim> &) const 
  {} 

  template <int dim> 
  void RHSIntegrator<dim>::boundary( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const FEValuesBase<dim> &fe           = info.fe_values(); 
    Vector<double> &         local_vector = dinfo.vector(0).block(0); 

    std::vector<double> boundary_values(fe.n_quadrature_points); 
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values); 

    const unsigned int degree = fe.get_fe().tensor_degree(); 
    const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
                           dinfo.cell->measure(); 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) 
        local_vector(i) += 
          (-penalty * fe.shape_value(i, k)              // (-sigma * v_i(x_k) 
           + fe.normal_vector(k) * fe.shape_grad(i, k)) // + n * grad v_i(x_k)) 
          * boundary_values[k] * fe.JxW(k);             // u^D(x_k) * dx 
  } 

  template <int dim> 
  void 
  RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &, 
                           MeshWorker::DoFInfo<dim> &, 
                           typename MeshWorker::IntegrationInfo<dim> &, 
                           typename MeshWorker::IntegrationInfo<dim> &) const 
  {} 

//第三个局部积分器负责对误差估计的贡献。这是由Karakashian和Pascal（2003）提出的标准能量估计器。

  template <int dim> 
  class Estimator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
              typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void 
         boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
                  typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
              MeshWorker::DoFInfo<dim> &                 dinfo2, 
              typename MeshWorker::IntegrationInfo<dim> &info1, 
              typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
  }; 

// 单元的贡献是离散解的拉普拉斯，因为右手边是零。

  template <int dim> 
  void 
  Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
                       typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const FEValuesBase<dim> &fe = info.fe_values(); 

    const std::vector<Tensor<2, dim>> &DDuh = info.hessians[0][0]; 
    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        const double t = dinfo.cell->diameter() * trace(DDuh[k]); 
        dinfo.value(0) += t * t * fe.JxW(k); 
      } 
    dinfo.value(0) = std::sqrt(dinfo.value(0)); 
  } 

// 在边界，我们简单地使用边界残差的加权形式，即有限元解和正确边界条件之间的差值的规范。

  template <int dim> 
  void Estimator<dim>::boundary( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const FEValuesBase<dim> &fe = info.fe_values(); 

    std::vector<double> boundary_values(fe.n_quadrature_points); 
    exact_solution.value_list(fe.get_quadrature_points(), boundary_values); 

    const std::vector<double> &uh = info.values[0][0]; 

    const unsigned int degree = fe.get_fe().tensor_degree(); 
    const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
                           dinfo.cell->measure(); 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        const double diff = boundary_values[k] - uh[k]; 
        dinfo.value(0) += penalty * diff * diff * fe.JxW(k); 
      } 
    dinfo.value(0) = std::sqrt(dinfo.value(0)); 
  } 

// 最后，在内部面，估计器由解的跳跃和它的法向导数组成，并进行适当的加权。

  template <int dim> 
  void 
  Estimator<dim>::face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
                       MeshWorker::DoFInfo<dim> &                 dinfo2, 
                       typename MeshWorker::IntegrationInfo<dim> &info1, 
                       typename MeshWorker::IntegrationInfo<dim> &info2) const 
  { 
    const FEValuesBase<dim> &          fe   = info1.fe_values(); 
    const std::vector<double> &        uh1  = info1.values[0][0]; 
    const std::vector<double> &        uh2  = info2.values[0][0]; 
    const std::vector<Tensor<1, dim>> &Duh1 = info1.gradients[0][0]; 
    const std::vector<Tensor<1, dim>> &Duh2 = info2.gradients[0][0]; 

    const unsigned int degree = fe.get_fe().tensor_degree(); 
    const double       penalty1 = 
      degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure(); 
    const double penalty2 = 
      degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure(); 
    const double penalty = penalty1 + penalty2; 
    const double h       = dinfo1.face->measure(); 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        const double diff1 = uh1[k] - uh2[k]; 
        const double diff2 = 
          fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k]; 
        dinfo1.value(0) += 
          (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k); 
      } 
    dinfo1.value(0) = std::sqrt(dinfo1.value(0)); 
    dinfo2.value(0) = dinfo1.value(0); 
  } 

// 最后我们有一个误差的积分器。由于不连续Galerkin问题的能量准则不仅涉及到单元内部的梯度差，还涉及到跨面和边界的跳跃项，所以我们不能仅仅使用  VectorTools::integrate_difference().  而是使用MeshWorker接口来自己计算误差。

//有几种不同的方法来定义这个能量准则，但是所有的方法都是随着网格大小的变化而等价的（有些不是随着多项式程度的变化而等价）。这里，我们选择
// @f[ \|u\|_{1,h} =
//  \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
//  4\sigma_F\|\average{ u \mathbf n}\|^2_F + \sum_{F \in F_h^b}
//  2\sigma_F\|u\|^2_F 
//  @f]

  template <int dim> 
  class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim> 
  { 
  public: 
    void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
              typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void 
         boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
                  typename MeshWorker::IntegrationInfo<dim> &info) const override; 
    void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
              MeshWorker::DoFInfo<dim> &                 dinfo2, 
              typename MeshWorker::IntegrationInfo<dim> &info1, 
              typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
  }; 

// 这里我们有关于单元格的集成。目前MeshWorker中还没有很好的接口可以让我们访问正交点中的正则函数值。因此，我们必须在单元格积分器中创建精确函数值和梯度的向量。之后，一切照旧，我们只需将差值的平方加起来。

// 除了计算能量准则的误差，我们还利用网格工作者的能力同时计算两个函数并在同一个循环中计算<i>L<sup>2</sup></i>的误差。很明显，这个函数没有任何跳跃项，只出现在单元格的积分中。

  template <int dim> 
  void ErrorIntegrator<dim>::cell( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const FEValuesBase<dim> &   fe = info.fe_values(); 
    std::vector<Tensor<1, dim>> exact_gradients(fe.n_quadrature_points); 
    std::vector<double>         exact_values(fe.n_quadrature_points); 

    exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients); 
    exact_solution.value_list(fe.get_quadrature_points(), exact_values); 

    const std::vector<Tensor<1, dim>> &Duh = info.gradients[0][0]; 
    const std::vector<double> &        uh  = info.values[0][0]; 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        double sum = 0; 
        for (unsigned int d = 0; d < dim; ++d) 
          { 
            const double diff = exact_gradients[k][d] - Duh[k][d]; 
            sum += diff * diff; 
          } 
        const double diff = exact_values[k] - uh[k]; 
        dinfo.value(0) += sum * fe.JxW(k); 
        dinfo.value(1) += diff * diff * fe.JxW(k); 
      } 
    dinfo.value(0) = std::sqrt(dinfo.value(0)); 
    dinfo.value(1) = std::sqrt(dinfo.value(1)); 
  } 

  template <int dim> 
  void ErrorIntegrator<dim>::boundary( 
    MeshWorker::DoFInfo<dim> &                 dinfo, 
    typename MeshWorker::IntegrationInfo<dim> &info) const 
  { 
    const FEValuesBase<dim> &fe = info.fe_values(); 

    std::vector<double> exact_values(fe.n_quadrature_points); 
    exact_solution.value_list(fe.get_quadrature_points(), exact_values); 

    const std::vector<double> &uh = info.values[0][0]; 

    const unsigned int degree = fe.get_fe().tensor_degree(); 
    const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
                           dinfo.cell->measure(); 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        const double diff = exact_values[k] - uh[k]; 
        dinfo.value(0) += penalty * diff * diff * fe.JxW(k); 
      } 
    dinfo.value(0) = std::sqrt(dinfo.value(0)); 
  } 

  template <int dim> 
  void ErrorIntegrator<dim>::face( 
    MeshWorker::DoFInfo<dim> &                 dinfo1, 
    MeshWorker::DoFInfo<dim> &                 dinfo2, 
    typename MeshWorker::IntegrationInfo<dim> &info1, 
    typename MeshWorker::IntegrationInfo<dim> &info2) const 
  { 
    const FEValuesBase<dim> &  fe  = info1.fe_values(); 
    const std::vector<double> &uh1 = info1.values[0][0]; 
    const std::vector<double> &uh2 = info2.values[0][0]; 

    const unsigned int degree = fe.get_fe().tensor_degree(); 
    const double       penalty1 = 
      degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure(); 
    const double penalty2 = 
      degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure(); 
    const double penalty = penalty1 + penalty2; 

    for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
      { 
        const double diff = uh1[k] - uh2[k]; 
        dinfo1.value(0) += (penalty * diff * diff) * fe.JxW(k); 
      } 
    dinfo1.value(0) = std::sqrt(dinfo1.value(0)); 
    dinfo2.value(0) = dinfo1.value(0); 
  } 

//  @sect3{The main class}  

// 这个类做主要的工作，就像前面的例子一样。关于这里声明的函数的描述，请参考下面的实现。

  template <int dim> 
  class InteriorPenaltyProblem 
  { 
  public: 
    using CellInfo = MeshWorker::IntegrationInfo<dim>; 

    InteriorPenaltyProblem(const FiniteElement<dim> &fe); 

    void run(unsigned int n_steps); 

  private: 
    void   setup_system(); 
    void   assemble_matrix(); 
    void   assemble_mg_matrix(); 
    void   assemble_right_hand_side(); 
    void   error(); 
    double estimate(); 
    void   solve(); 
    void   output_results(const unsigned int cycle) const; 

// 与离散化有关的成员对象在这里。

    Triangulation<dim>        triangulation; 
    const MappingQ1<dim>      mapping; 
    const FiniteElement<dim> &fe; 
    DoFHandler<dim>           dof_handler; 

// 然后，我们有与全局离散系统相关的矩阵和向量。

    SparsityPattern      sparsity; 
    SparseMatrix<double> matrix; 
    Vector<double>       solution; 
    Vector<double>       right_hand_side; 
    BlockVector<double>  estimates; 

// 最后，我们有一组与多级预处理程序相关的稀疏模式和稀疏矩阵。 首先，我们有一个水平矩阵和它的稀疏性模式。

    MGLevelObject<SparsityPattern>      mg_sparsity; 
    MGLevelObject<SparseMatrix<double>> mg_matrix; 

// 当我们在局部细化的网格上进行局部平滑的多重网格时，需要额外的矩阵；见Kanschat（2004）。这里是这些边缘矩阵的稀疏性模式。我们只需要一个，因为上矩阵的模式是下矩阵的转置。实际上，我们并不太关心这些细节，因为MeshWorker正在填充这些矩阵。

    MGLevelObject<SparsityPattern> mg_sparsity_dg_interface; 

// 精细化边缘的通量矩阵，将精细级自由度与粗略级自由度相耦合。

    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_down; 

// 精细化边缘的通量矩阵的转置，将粗级自由度耦合到精细级。

    MGLevelObject<SparseMatrix<double>> mg_matrix_dg_up; 
  }; 

// 构造函数简单地设置了粗略的网格和DoFHandler。FiniteElement作为一个参数被提供，以实现灵活性。

  template <int dim> 
  InteriorPenaltyProblem<dim>::InteriorPenaltyProblem( 
    const FiniteElement<dim> &fe) 
    : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
    , mapping() 
    , fe(fe) 
    , dof_handler(triangulation) 
    , estimates(1) 
  { 
    GridGenerator::hyper_cube_slit(triangulation, -1, 1); 
  } 

// 在这个函数中，我们设置了线性系统的维度和全局矩阵以及水平矩阵的稀疏性模式。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::setup_system() 
  { 

// 首先，我们用有限元将自由度分布在网格上并对其进行编号。

    dof_handler.distribute_dofs(fe); 
    dof_handler.distribute_mg_dofs(); 
    unsigned int n_dofs = dof_handler.n_dofs(); 

// 然后，我们已经知道代表有限元函数的向量的大小。

    solution.reinit(n_dofs); 
    right_hand_side.reinit(n_dofs); 

// 接下来，我们为全局矩阵设置稀疏性模式。由于我们事先不知道行的大小，所以我们首先填充一个临时的DynamicSparsityPattern对象，一旦完成，就将其复制到常规的SparsityPattern中。

    DynamicSparsityPattern dsp(n_dofs); 
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
    sparsity.copy_from(dsp); 
    matrix.reinit(sparsity); 

    const unsigned int n_levels = triangulation.n_levels(); 

// 全局系统已经设置好了，现在我们来关注一下级别矩阵。我们调整所有矩阵对象的大小，以便每一级都有一个矩阵。

    mg_matrix.resize(0, n_levels - 1); 
    mg_matrix.clear_elements(); 
    mg_matrix_dg_up.resize(0, n_levels - 1); 
    mg_matrix_dg_up.clear_elements(); 
    mg_matrix_dg_down.resize(0, n_levels - 1); 
    mg_matrix_dg_down.clear_elements(); 

// 在为水平矩阵调用<tt>clear()</tt>之后更新稀疏模式很重要，因为矩阵通过SmartPointer和Subscriptor机制锁定了稀疏模式。

    mg_sparsity.resize(0, n_levels - 1); 
    mg_sparsity_dg_interface.resize(0, n_levels - 1); 

// 现在，所有的对象都准备好了，可以在每一层容纳一个稀疏模式或矩阵。剩下的就是在每一层设置稀疏模式了。

    for (unsigned int level = mg_sparsity.min_level(); 
         level <= mg_sparsity.max_level(); 
         ++level) 
      { 

// 这些与上面的全局矩阵的行数大致相同，现在是每个级别的。

        DynamicSparsityPattern dsp(dof_handler.n_dofs(level)); 
        MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level); 
        mg_sparsity[level].copy_from(dsp); 
        mg_matrix[level].reinit(mg_sparsity[level]); 

// 另外，我们需要初始化各层之间细化边缘的转移矩阵。它们被存储在两个索引中较细的索引处，因此在0层没有这样的对象。

        if (level > 0) 
          { 
            DynamicSparsityPattern dsp; 
            dsp.reinit(dof_handler.n_dofs(level - 1), 
                       dof_handler.n_dofs(level)); 
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, dsp, level); 
            mg_sparsity_dg_interface[level].copy_from(dsp); 
            mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]); 
            mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]); 
          } 
      } 
  } 

// 在这个函数中，我们组装全局系统矩阵，这里的全局是指我们解决的离散系统的矩阵，它覆盖了整个网格。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::assemble_matrix() 
  { 

// 首先，我们需要设置提供我们集成值的对象。这个对象包含了所有需要的FEValues和FEFaceValues对象，并且自动维护它们，使它们总是指向当前单元。为此，我们首先需要告诉它，在哪里计算，计算什么。由于我们没有做任何花哨的事情，我们可以依靠他们对正交规则的标准选择。

// 由于他们的默认更新标志是最小的，我们另外添加我们需要的东西，即所有对象（单元格、边界和内部面）上的形状函数的值和梯度。之后，我们准备初始化容器，它将创建所有必要的FEValuesBase对象进行整合。

    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags update_flags = update_values | update_gradients; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

// 这就是我们整合本地数据的对象。它由MatrixIntegrator中的局部整合例程填充，然后由汇编器用来将信息分配到全局矩阵中。

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

// 此外，我们还需要一个将局部矩阵装配到全局矩阵的对象。这些装配器对象拥有目标对象结构的所有知识，在这里是一个稀疏矩阵，可能的约束和网格结构。

    MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler; 
    assembler.initialize(matrix); 

// 现在是我们自己编码的部分，局部积分器。这是唯一与问题有关的部分。

    MatrixIntegrator<dim> integrator; 

// 现在，我们把所有的东西都扔到 MeshWorker::loop(), 中，在这里遍历网格的所有活动单元，计算单元和面的矩阵，并把它们集合到全局矩阵中。我们在这里使用变量<tt>dof_handler</tt>，以便使用全局自由度的编号。

    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           integrator, 
                                           assembler); 
  } 

// 现在，我们对水平矩阵做同样的处理。不太令人惊讶的是，这个函数看起来像前一个函数的孪生兄弟。事实上，只有两个小的区别。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::assemble_mg_matrix() 
  { 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags update_flags = update_values | update_gradients; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

// 很明显，需要用一个填充水平矩阵的汇编器来代替。请注意，它也会自动填充边缘矩阵。

    MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler; 
    assembler.initialize(mg_matrix); 
    assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down); 

    MatrixIntegrator<dim> integrator; 

// 这里是与前一个函数的另一个不同之处：我们在所有单元上运行，而不仅仅是活动单元。而且我们使用以 <code>_mg</code> 结尾的函数，因为我们需要每一层的自由度，而不是全局的编号。

    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(), 
                                           dof_handler.end_mg(), 
                                           dof_info, 
                                           info_box, 
                                           integrator, 
                                           assembler); 
  } 

// 这里我们有另一个assemble函数的克隆。与组装系统矩阵的区别在于，我们在这里组装了一个向量。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::assemble_right_hand_side() 
  { 
    MeshWorker::IntegrationInfoBox<dim> info_box; 
    UpdateFlags                         update_flags = 
      update_quadrature_points | update_values | update_gradients; 
    info_box.add_update_flags_all(update_flags); 
    info_box.initialize(fe, mapping); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

// 因为这个汇编器允许我们填充多个向量，所以接口要比上面复杂一些。向量的指针必须存储在一个AnyData对象中。虽然这在这里似乎造成了两行额外的代码，但实际上在更复杂的应用中它是很方便的。

    MeshWorker::Assembler::ResidualSimple<Vector<double>> assembler; 
    AnyData                                               data; 
    data.add<Vector<double> *>(&right_hand_side, "RHS"); 
    assembler.initialize(data); 

    RHSIntegrator<dim> integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           integrator, 
                                           assembler); 

    right_hand_side *= -1.; 
  } 

// 现在，我们已经对构建离散线性系统的所有函数进行了编码，现在是我们实际解决它的时候了。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::solve() 
  { 

// 选择的求解器是共轭梯度。

    SolverControl            control(1000, 1.e-12); 
    SolverCG<Vector<double>> solver(control); 

// 现在我们正在设置多级预处理程序的组件。首先，我们需要在网格层之间进行转移。我们在这里使用的对象为这些转移生成了稀疏矩阵。

    MGTransferPrebuilt<Vector<double>> mg_transfer; 
    mg_transfer.build(dof_handler); 

// 然后，我们需要一个精确的解算器来解算最粗层次上的矩阵。

    FullMatrix<double> coarse_matrix; 
    coarse_matrix.copy_from(mg_matrix[0]); 
    MGCoarseGridHouseholder<double, Vector<double>> mg_coarse; 
    mg_coarse.initialize(coarse_matrix); 

// 虽然转移和粗略网格求解器几乎是通用的，但为平滑器提供了更多的灵活性。首先，我们选择Gauss-Seidel作为我们的平滑方法。

    GrowingVectorMemory<Vector<double>> mem; 
    using RELAXATION = PreconditionSOR<SparseMatrix<double>>; 
    mg::SmootherRelaxation<RELAXATION, Vector<double>> mg_smoother; 
    RELAXATION::AdditionalData                         smoother_data(1.); 
    mg_smoother.initialize(mg_matrix, smoother_data); 

// 在每个级别上做两个平滑步骤。

    mg_smoother.set_steps(2); 

// 由于SOR方法不是对称的，但我们在下面使用共轭梯度迭代，这里有一个技巧，使多级预处理器成为对称算子，即使是对非对称平滑器。

    mg_smoother.set_symmetric(true); 

// 平滑器类可以选择实现变量V型循环，我们在这里不需要。

    mg_smoother.set_variable(false); 

// 最后，我们必须将我们的矩阵包裹在一个具有所需乘法函数的对象中。

    mg::Matrix<Vector<double>> mgmatrix(mg_matrix); 
    mg::Matrix<Vector<double>> mgdown(mg_matrix_dg_down); 
    mg::Matrix<Vector<double>> mgup(mg_matrix_dg_up); 

// 现在，我们准备设置V型循环算子和多级预处理程序。

    Multigrid<Vector<double>> mg( 
      mgmatrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 

// 让我们不要忘记因为自适应细化而需要的边缘矩阵。

    mg.set_edge_flux_matrices(mgdown, mgup); 

// 在所有的准备工作完成后，将Multigrid对象包装成另一个对象，它可以作为一个普通的预处理程序使用。

    PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
      preconditioner(dof_handler, mg, mg_transfer); 

// 并用它来解决这个系统。

    solver.solve(matrix, solution, right_hand_side, preconditioner); 
  } 

// 另一个克隆的集合函数。与之前的最大区别是，这里我们也有一个输入向量。

  template <int dim> 
  double InteriorPenaltyProblem<dim>::estimate() 
  { 

// 估算器的结果存储在一个每个单元格有一个条目的向量中。由于deal.II中的单元格没有编号，我们必须建立自己的编号，以便使用这个向量。对于下面使用的汇编器来说，结果存储在向量的哪个分量中的信息是由每个单元的user_index变量传送的。我们需要在这里设置这个编号。

// 另一方面，有人可能已经使用了用户指数。所以，让我们做个好公民，在篡改它们之前保存它们。

    std::vector<unsigned int> old_user_indices; 
    triangulation.save_user_indices(old_user_indices); 

    estimates.block(0).reinit(triangulation.n_active_cells()); 
    unsigned int i = 0; 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      cell->set_user_index(i++); 

// 这就像以前一样开始。

    MeshWorker::IntegrationInfoBox<dim> info_box; 
    const unsigned int                  n_gauss_points = 
      dof_handler.get_fe().tensor_degree() + 1; 
    info_box.initialize_gauss_quadrature(n_gauss_points, 
                                         n_gauss_points + 1, 
                                         n_gauss_points); 

// 但现在我们需要通知信息框我们要在正交点上评估的有限元函数。首先，我们用这个向量创建一个AnyData对象，这个向量就是我们刚刚计算的解。

    AnyData solution_data; 
    solution_data.add<const Vector<double> *>(&solution, "solution"); 

// 然后，我们告诉单元格的 Meshworker::VectorSelector ，我们需要这个解决方案的二次导数（用来计算拉普拉斯）。因此，选择函数值和第一导数的布尔参数是假的，只有选择第二导数的最后一个参数是真的。

    info_box.cell_selector.add("solution", false, false, true); 

// 在内部和边界面，我们需要函数值和第一导数，但不需要第二导数。

    info_box.boundary_selector.add("solution", true, true, false); 
    info_box.face_selector.add("solution", true, true, false); 

// 我们继续像以前一样，除了默认的更新标志已经被调整为我们上面要求的值和导数之外。

    info_box.add_update_flags_boundary(update_quadrature_points); 
    info_box.initialize(fe, mapping, solution_data, solution); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

// 汇编器在每个单元格中存储一个数字，否则这与右侧的计算是一样的。

    MeshWorker::Assembler::CellsAndFaces<double> assembler; 
    AnyData                                      out_data; 
    out_data.add<BlockVector<double> *>(&estimates, "cells"); 
    assembler.initialize(out_data, false); 

    Estimator<dim> integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           integrator, 
                                           assembler); 

// 就在我们返回错误估计的结果之前，我们恢复旧的用户索引。

    triangulation.load_user_indices(old_user_indices); 
    return estimates.block(0).l2_norm(); 
  } 

// 这里我们把我们的有限元解和（已知的）精确解进行比较，计算梯度和函数本身的平均二次误差。这个函数是上面那个估计函数的克隆。

// 由于我们分别计算能量和<i>L<sup>2</sup></i>-norm的误差，我们的块向量在这里需要两个块。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::error() 
  { 
    BlockVector<double> errors(2); 
    errors.block(0).reinit(triangulation.n_active_cells()); 
    errors.block(1).reinit(triangulation.n_active_cells()); 

    std::vector<unsigned int> old_user_indices; 
    triangulation.save_user_indices(old_user_indices); 
    unsigned int i = 0; 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      cell->set_user_index(i++); 

    MeshWorker::IntegrationInfoBox<dim> info_box; 
    const unsigned int                  n_gauss_points = 
      dof_handler.get_fe().tensor_degree() + 1; 
    info_box.initialize_gauss_quadrature(n_gauss_points, 
                                         n_gauss_points + 1, 
                                         n_gauss_points); 

    AnyData solution_data; 
    solution_data.add<Vector<double> *>(&solution, "solution"); 

    info_box.cell_selector.add("solution", true, true, false); 
    info_box.boundary_selector.add("solution", true, false, false); 
    info_box.face_selector.add("solution", true, false, false); 

    info_box.add_update_flags_cell(update_quadrature_points); 
    info_box.add_update_flags_boundary(update_quadrature_points); 
    info_box.initialize(fe, mapping, solution_data, solution); 

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

    MeshWorker::Assembler::CellsAndFaces<double> assembler; 
    AnyData                                      out_data; 
    out_data.add<BlockVector<double> *>(&errors, "cells"); 
    assembler.initialize(out_data, false); 

    ErrorIntegrator<dim> integrator; 
    MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
                                           dof_handler.end(), 
                                           dof_info, 
                                           info_box, 
                                           integrator, 
                                           assembler); 
    triangulation.load_user_indices(old_user_indices); 

    deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl; 
    deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl; 
  } 

// 创建图形输出。我们通过整理其各个组成部分的名称来产生文件名，包括我们用两个数字输出的细化周期。

  template <int dim> 
  void 
  InteriorPenaltyProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    const std::string filename = 
      "sol-" + Utilities::int_to_string(cycle, 2) + ".gnuplot"; 

    deallog << "Writing solution to <" << filename << ">..." << std::endl 
            << std::endl; 
    std::ofstream gnuplot_output(filename); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "u"); 
    data_out.add_data_vector(estimates.block(0), "est"); 

    data_out.build_patches(); 

    data_out.write_gnuplot(gnuplot_output); 
  } 

// 最后是自适应循环，或多或少和前面的例子一样。

  template <int dim> 
  void InteriorPenaltyProblem<dim>::run(unsigned int n_steps) 
  { 
    deallog << "Element: " << fe.get_name() << std::endl; 
    for (unsigned int s = 0; s < n_steps; ++s) 
      { 
        deallog << "Step " << s << std::endl; 
        if (estimates.block(0).size() == 0) 
          triangulation.refine_global(1); 
        else 
          { 
            GridRefinement::refine_and_coarsen_fixed_fraction( 
              triangulation, estimates.block(0), 0.5, 0.0); 
            triangulation.execute_coarsening_and_refinement(); 
          } 

        deallog << "Triangulation " << triangulation.n_active_cells() 
                << " cells, " << triangulation.n_levels() << " levels" 
                << std::endl; 

        setup_system(); 
        deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs"; 
        for (unsigned int l = 0; l < triangulation.n_levels(); ++l) 
          deallog << ' ' << dof_handler.n_dofs(l); 
        deallog << std::endl; 

        deallog << "Assemble matrix" << std::endl; 
        assemble_matrix(); 
        deallog << "Assemble multilevel matrix" << std::endl; 
        assemble_mg_matrix(); 
        deallog << "Assemble right hand side" << std::endl; 
        assemble_right_hand_side(); 
        deallog << "Solve" << std::endl; 
        solve(); 
        error(); 
        deallog << "Estimate " << estimate() << std::endl; 
        output_results(s); 
      } 
  } 
} // namespace Step39 

int main() 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step39; 

      deallog.depth_console(2); 
      std::ofstream logfile("deallog"); 
      deallog.attach(logfile); 
      FE_DGQ<2>                 fe1(3); 
      InteriorPenaltyProblem<2> test1(fe1); 
      test1.run(12); 
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

