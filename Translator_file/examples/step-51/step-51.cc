


/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2013 - 2021 by the deal.II authors 
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
 * Author: Martin Kronbichler, Technische Universität München, 
 *         Scott T. Miller, The Pennsylvania State University, 2013 
 */ 


// @sect3{Include files}  

// 大多数deal.II的include文件已经在前面的例子中涉及到了，没有注释。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/tensor_function.h> 
#include <deal.II/base/exceptions.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/work_stream.h> 
#include <deal.II/base/convergence_table.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_bicgstab.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

// 然而，我们确实有一些新的包括在这个例子中。第一个定义了三角形面的有限元空间，我们把它称为 "骨架"。这些有限元在元素内部没有任何支持，它们代表的是在每个模数一的表面上有一个单一的值的多项式，但在模数二的表面上允许有不连续。

#include <deal.II/fe/fe_face.h> 

// 我们包含的第二个新文件定义了一种新的稀疏矩阵类型。 常规的 <code>SparseMatrix</code> 类型存储了所有非零条目的索引。  <code>ChunkSparseMatrix</code> 则是利用了DG解的耦合性。 它存储了一个指定大小的矩阵子块的索引。 在HDG背景下，这个子块大小实际上是由骨架解场定义的每个面的自由度数量。这使得矩阵的内存消耗减少了三分之一，并且在求解器中使用矩阵时也会有类似的速度提升。

#include <deal.II/lac/chunk_sparse_matrix.h> 

// 这个例子的最后一个新的包括涉及到数据输出。 由于我们在网格的骨架上定义了一个有限元场，我们希望能够直观地看到这个解决方案的实际情况。DataOutFaces正是这样做的；它的接口与我们熟悉的DataOut几乎一样，但输出的数据只有模拟的二维1数据。

#include <deal.II/numerics/data_out_faces.h> 

#include <iostream> 

// 我们首先将所有的类放入自己的命名空间。

namespace Step51 
{ 
  using namespace dealii; 
// @sect3{Equation data}  

//分析解的结构与 step-7 中相同。有两个例外情况。首先，我们也为3D情况创建了一个解决方案，其次，我们对解决方案进行了缩放，使其在解决方案的所有宽度值上的规范是统一的。

  template <int dim> 
  class SolutionBase 
  { 
  protected: 
    static const unsigned int n_source_centers = 3; 
    static const Point<dim>   source_centers[n_source_centers]; 
    static const double       width; 
  }; 

  template <> 
  const Point<1> 
    SolutionBase<1>::source_centers[SolutionBase<1>::n_source_centers] = 
      {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}; 

  template <> 
  const Point<2> 
    SolutionBase<2>::source_centers[SolutionBase<2>::n_source_centers] = 
      {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}; 

  template <> 
  const Point<3> 
    SolutionBase<3>::source_centers[SolutionBase<3>::n_source_centers] = { 
      Point<3>(-0.5, +0.5, 0.25), 
      Point<3>(-0.6, -0.5, -0.125), 
      Point<3>(+0.5, -0.5, 0.5)}; 

  template <int dim> 
  const double SolutionBase<dim>::width = 1. / 5.; 

  template <int dim> 
  class Solution : public Function<dim>, protected SolutionBase<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      double sum = 0; 
      for (unsigned int i = 0; i < this->n_source_centers; ++i) 
        { 
          const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 
          sum += 
            std::exp(-x_minus_xi.norm_square() / (this->width * this->width)); 
        } 

      return sum / 
             std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
    } 

    virtual Tensor<1, dim> 
    gradient(const Point<dim> &p, 
             const unsigned int /*component*/ = 0) const override 
    { 
      Tensor<1, dim> sum; 
      for (unsigned int i = 0; i < this->n_source_centers; ++i) 
        { 
          const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 

          sum += 
            (-2 / (this->width * this->width) * 
             std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) * 
             x_minus_xi); 
        } 

      return sum / 
             std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
    } 
  }; 

// 这个类实现了一个函数，标量解和它的负梯度被收集在一起。这个函数在计算HDG近似的误差时使用，它的实现是简单地调用Solution类的值和梯度函数。

  template <int dim> 
  class SolutionAndGradient : public Function<dim>, protected SolutionBase<dim> 
  { 
  public: 
    SolutionAndGradient() 
      : Function<dim>(dim + 1) 
    {} 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  v) const override 
    { 
      AssertDimension(v.size(), dim + 1); 
      Solution<dim>  solution; 
      Tensor<1, dim> grad = solution.gradient(p); 
      for (unsigned int d = 0; d < dim; ++d) 
        v[d] = -grad[d]; 
      v[dim] = solution.value(p); 
    } 
  }; 

// 接下来是对流速度的实现。如介绍中所述，我们选择的速度场在二维是 $(y, -x)$ ，在三维是 $(y, -x, 1)$ 。这就得到了一个无发散的速度场。

  template <int dim> 
  class ConvectionVelocity : public TensorFunction<1, dim> 
  { 
  public: 
    ConvectionVelocity() 
      : TensorFunction<1, dim>() 
    {} 

    virtual Tensor<1, dim> value(const Point<dim> &p) const override 
    { 
      Tensor<1, dim> convection; 
      switch (dim) 
        { 
          case 1: 
            convection[0] = 1; 
            break; 
          case 2: 
            convection[0] = p[1]; 
            convection[1] = -p[0]; 
            break; 
          case 3: 
            convection[0] = p[1]; 
            convection[1] = -p[0]; 
            convection[2] = 1; 
            break; 
          default: 
            Assert(false, ExcNotImplemented()); 
        } 
      return convection; 
    } 
  }; 

// 我们实现的最后一个函数是用于制造解决方案的右手边。它与 step-7 非常相似，不同的是我们现在有一个对流项而不是反应项。由于速度场是不可压缩的，即 $\nabla \cdot \mathbf{c} =0$ ，对流项简单读作  $\mathbf{c} \nabla u$  。

  template <int dim> 
  class RightHandSide : public Function<dim>, protected SolutionBase<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      ConvectionVelocity<dim> convection_velocity; 
      Tensor<1, dim>          convection = convection_velocity.value(p); 
      double                  sum        = 0; 
      for (unsigned int i = 0; i < this->n_source_centers; ++i) 
        { 
          const Tensor<1, dim> x_minus_xi = p - this->source_centers[i]; 

          sum += 
            ((2 * dim - 2 * convection * x_minus_xi - 
              4 * x_minus_xi.norm_square() / (this->width * this->width)) / 
             (this->width * this->width) * 
             std::exp(-x_minus_xi.norm_square() / (this->width * this->width))); 
        } 

      return sum / 
             std::pow(2. * numbers::PI * this->width * this->width, dim / 2.); 
    } 
  }; 

//  @sect3{The HDG solver class}  

// HDG的求解过程与  step-7  的求解过程非常相似。主要区别在于使用了三套不同的DoFHandler和FE对象，以及ChunkSparseMatrix和相应的解决方案向量。我们还使用WorkStream来实现多线程的本地求解过程，该过程利用了本地求解器的尴尬的并行性质。对于WorkStream，我们定义了对单元格的本地操作和复制到全局矩阵和向量的函数。我们这样做既是为了装配（装配要运行两次，一次是在我们生成系统矩阵时，另一次是在我们从骨架值计算元素内部解时），也是为了后处理，在后处理中我们提取一个在高阶收敛的解。

  template <int dim> 
  class HDG 
  { 
  public: 
    enum RefinementMode 
    { 
      global_refinement, 
      adaptive_refinement 
    }; 

    HDG(const unsigned int degree, const RefinementMode refinement_mode); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(const bool reconstruct_trace = false); 
    void solve(); 
    void postprocess(); 
    void refine_grid(const unsigned int cycle); 
    void output_results(const unsigned int cycle); 

// 用于组装和解决原始变量的数据。

    struct PerTaskData; 
    struct ScratchData; 

// 对解决方案进行后处理以获得  $u^*$  是一个逐个元素的过程；因此，我们不需要组装任何全局数据，也不需要声明任何 "任务数据 "供WorkStream使用。

    struct PostProcessScratchData; 

// 以下三个函数被 WorkStream 用来完成程序的实际工作。

    void assemble_system_one_cell( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      ScratchData &                                         scratch, 
      PerTaskData &                                         task_data); 

    void copy_local_to_global(const PerTaskData &data); 

    void postprocess_one_cell( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      PostProcessScratchData &                              scratch, 
      unsigned int &                                        empty_data); 

    Triangulation<dim> triangulation; 

// "局部 "解是每个元素的内部。 这些代表了原始解场  $u$  以及辅助场  $\mathbf{q}$  。

    FESystem<dim>   fe_local; 
    DoFHandler<dim> dof_handler_local; 
    Vector<double>  solution_local; 

// 新的有限元类型和相应的 <code>DoFHandler</code> 被用于耦合元素级局部解的全局骨架解。

    FE_FaceQ<dim>   fe; 
    DoFHandler<dim> dof_handler; 
    Vector<double>  solution; 
    Vector<double>  system_rhs; 

// 如介绍中所述，HDG解可以通过后处理达到  $\mathcal{O}(h^{p+2})$  的超收敛率。 后处理的解是一个不连续的有限元解，代表每个单元内部的原始变量。 我们定义了一个程度为 $p+1$ 的FE类型来表示这个后处理的解，我们只在构造后用于输出。

    FE_DGQ<dim>     fe_u_post; 
    DoFHandler<dim> dof_handler_u_post; 
    Vector<double>  solution_u_post; 

// 与骨架相对应的自由度强烈地执行Dirichlet边界条件，就像在连续Galerkin有限元方法中一样。我们可以通过AffineConstraints对象以类似的方式强制执行边界条件。此外，悬挂节点的处理方式与连续有限元的处理方式相同。对于只在面定义自由度的面元素，这个过程将精炼面的解设置为与粗略面的表示相吻合。

// 请注意，对于HDG来说，消除悬空节点并不是唯一的可能性，就HDG理论而言，我们也可以使用精炼侧的未知数，通过精炼侧的跟踪值来表达粗略侧的局部解。然而，这样的设置在deal.II循环方面并不容易实现，因此没有进一步分析。

    AffineConstraints<double> constraints; 

// ChunkSparseMatrix类的用法与通常的稀疏矩阵类似。你需要一个ChunkSparsityPattern类型的稀疏模式和实际的矩阵对象。在创建稀疏模式时，我们只需要额外传递局部块的大小。

    ChunkSparsityPattern      sparsity_pattern; 
    ChunkSparseMatrix<double> system_matrix; 

// 与  step-7  相同。

    const RefinementMode refinement_mode; 
    ConvergenceTable     convergence_table; 
  }; 
// @sect3{The HDG class implementation}  
// @sect4{Constructor}  该构造函数与其他例子中的构造函数类似，除了处理多个DoFHandler和FiniteElement对象。请注意，我们为局部DG部分创建了一个有限元系统，包括梯度/通量部分和标量部分。

  template <int dim> 
  HDG<dim>::HDG(const unsigned int degree, const RefinementMode refinement_mode) 
    : fe_local(FE_DGQ<dim>(degree), dim, FE_DGQ<dim>(degree), 1) 
    , dof_handler_local(triangulation) 
    , fe(degree) 
    , dof_handler(triangulation) 
    , fe_u_post(degree + 1) 
    , dof_handler_u_post(triangulation) 
    , refinement_mode(refinement_mode) 
  {} 

//  @sect4{HDG::setup_system}  HDG解决方案的系统是以类似于其他大多数教程程序的方式设置的。 我们小心翼翼地用我们所有的DoFHandler对象来分配道夫。  @p solution 和 @p system_matrix 对象与全局骨架解决方案一起。

  template <int dim> 
  void HDG<dim>::setup_system() 
  { 
    dof_handler_local.distribute_dofs(fe_local); 
    dof_handler.distribute_dofs(fe); 
    dof_handler_u_post.distribute_dofs(fe_u_post); 

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    solution_local.reinit(dof_handler_local.n_dofs()); 
    solution_u_post.reinit(dof_handler_u_post.n_dofs()); 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
    std::map<types::boundary_id, const Function<dim> *> boundary_functions; 
    Solution<dim>                                       solution_function; 
    boundary_functions[0] = &solution_function; 
    VectorTools::project_boundary_values(dof_handler, 
                                         boundary_functions, 
                                         QGauss<dim - 1>(fe.degree + 1), 
                                         constraints); 
    constraints.close(); 

// 在创建块状稀疏模式时，我们首先创建通常的动态稀疏模式，然后设置块状大小，该大小等于一个面的道夫数，当把它复制到最终的稀疏模式时。

    { 
      DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
      sparsity_pattern.copy_from(dsp, fe.n_dofs_per_face()); 
    } 
    system_matrix.reinit(sparsity_pattern); 
  } 

//  @sect4{HDG::PerTaskData}  接下来是定义并行装配的本地数据结构。第一个结构 @p PerTaskData 包含了被写入全局矩阵的本地向量和矩阵，而ScratchData包含了我们在本地装配中需要的所有数据。这里有一个变量值得注意，即布尔变量 @p  trace_reconstruct。正如介绍中提到的，我们分两步解决HDG系统。首先，我们为骨架系统创建一个线性系统，通过舒尔补码 $D-CA^{-1}B$  将局部部分浓缩到其中。然后，我们用骨架的解来解决局部部分。对于这两个步骤，我们需要两次元素上的相同矩阵，我们希望通过两个装配步骤来计算。由于大部分的代码是相似的，我们用相同的函数来做这件事，但只是根据我们在开始装配时设置的一个标志在两者之间切换。因为我们需要把这个信息传递给本地的工作程序，所以我们把它存储在任务数据中一次。

  template <int dim> 
  struct HDG<dim>::PerTaskData 
  { 
    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_vector; 
    std::vector<types::global_dof_index> dof_indices; 

    bool trace_reconstruct; 

    PerTaskData(const unsigned int n_dofs, const bool trace_reconstruct) 
      : cell_matrix(n_dofs, n_dofs) 
      , cell_vector(n_dofs) 
      , dof_indices(n_dofs) 
      , trace_reconstruct(trace_reconstruct) 
    {} 
  }; 

//  @sect4{HDG::ScratchData}  
// @p ScratchData  包含WorkStream中每个线程的持久化数据。 FEValues、矩阵和矢量对象现在应该很熟悉了。 有两个对象需要讨论。  `std::vector<std::vector<unsigned  int> > fe_local_support_on_face` 和  `std::vector<std::vector<unsigned  int> > fe_support_on_face`。 这些用于指示所选择的有限元是否在与 @p fe_local 相关的局部部分和骨架部分 @p fe. 的参考单元的特定面上有支持（非零值）。 我们在构造函数中提取这一信息，并为我们工作的所有单元存储一次。 如果我们不存储这一信息，我们将被迫在每个单元上装配大量的零项，这将大大降低程序的速度。

  template <int dim> 
  struct HDG<dim>::ScratchData 
  { 
    FEValues<dim>     fe_values_local; 
    FEFaceValues<dim> fe_face_values_local; 
    FEFaceValues<dim> fe_face_values; 

    FullMatrix<double> ll_matrix; 
    FullMatrix<double> lf_matrix; 
    FullMatrix<double> fl_matrix; 
    FullMatrix<double> tmp_matrix; 
    Vector<double>     l_rhs; 
    Vector<double>     tmp_rhs; 

    std::vector<Tensor<1, dim>> q_phi; 
    std::vector<double>         q_phi_div; 
    std::vector<double>         u_phi; 
    std::vector<Tensor<1, dim>> u_phi_grad; 
    std::vector<double>         tr_phi; 
    std::vector<double>         trace_values; 

    std::vector<std::vector<unsigned int>> fe_local_support_on_face; 
    std::vector<std::vector<unsigned int>> fe_support_on_face; 

    ConvectionVelocity<dim> convection_velocity; 
    RightHandSide<dim>      right_hand_side; 
    const Solution<dim>     exact_solution; 

    ScratchData(const FiniteElement<dim> &fe, 
                const FiniteElement<dim> &fe_local, 
                const QGauss<dim> &       quadrature_formula, 
                const QGauss<dim - 1> &   face_quadrature_formula, 
                const UpdateFlags         local_flags, 
                const UpdateFlags         local_face_flags, 
                const UpdateFlags         flags) 
      : fe_values_local(fe_local, quadrature_formula, local_flags) 
      , fe_face_values_local(fe_local, 
                             face_quadrature_formula, 
                             local_face_flags) 
      , fe_face_values(fe, face_quadrature_formula, flags) 
      , ll_matrix(fe_local.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
      , lf_matrix(fe_local.n_dofs_per_cell(), fe.n_dofs_per_cell()) 
      , fl_matrix(fe.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
      , tmp_matrix(fe.n_dofs_per_cell(), fe_local.n_dofs_per_cell()) 
      , l_rhs(fe_local.n_dofs_per_cell()) 
      , tmp_rhs(fe_local.n_dofs_per_cell()) 
      , q_phi(fe_local.n_dofs_per_cell()) 
      , q_phi_div(fe_local.n_dofs_per_cell()) 
      , u_phi(fe_local.n_dofs_per_cell()) 
      , u_phi_grad(fe_local.n_dofs_per_cell()) 
      , tr_phi(fe.n_dofs_per_cell()) 
      , trace_values(face_quadrature_formula.size()) 
      , fe_local_support_on_face(GeometryInfo<dim>::faces_per_cell) 
      , fe_support_on_face(GeometryInfo<dim>::faces_per_cell) 
      , exact_solution() 
    { 
      for (unsigned int face_no : GeometryInfo<dim>::face_indices()) 
        for (unsigned int i = 0; i < fe_local.n_dofs_per_cell(); ++i) 
          { 
            if (fe_local.has_support_on_face(i, face_no)) 
              fe_local_support_on_face[face_no].push_back(i); 
          } 

      for (unsigned int face_no : GeometryInfo<dim>::face_indices()) 
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
          { 
            if (fe.has_support_on_face(i, face_no)) 
              fe_support_on_face[face_no].push_back(i); 
          } 
    } 

    ScratchData(const ScratchData &sd) 
      : fe_values_local(sd.fe_values_local.get_fe(), 
                        sd.fe_values_local.get_quadrature(), 
                        sd.fe_values_local.get_update_flags()) 
      , fe_face_values_local(sd.fe_face_values_local.get_fe(), 
                             sd.fe_face_values_local.get_quadrature(), 
                             sd.fe_face_values_local.get_update_flags()) 
      , fe_face_values(sd.fe_face_values.get_fe(), 
                       sd.fe_face_values.get_quadrature(), 
                       sd.fe_face_values.get_update_flags()) 
      , ll_matrix(sd.ll_matrix) 
      , lf_matrix(sd.lf_matrix) 
      , fl_matrix(sd.fl_matrix) 
      , tmp_matrix(sd.tmp_matrix) 
      , l_rhs(sd.l_rhs) 
      , tmp_rhs(sd.tmp_rhs) 
      , q_phi(sd.q_phi) 
      , q_phi_div(sd.q_phi_div) 
      , u_phi(sd.u_phi) 
      , u_phi_grad(sd.u_phi_grad) 
      , tr_phi(sd.tr_phi) 
      , trace_values(sd.trace_values) 
      , fe_local_support_on_face(sd.fe_local_support_on_face) 
      , fe_support_on_face(sd.fe_support_on_face) 
      , exact_solution() 
    {} 
  }; 

//  @sect4{HDG::PostProcessScratchData}  
// @p PostProcessScratchData  包含WorkStream在对本地解决方案进行后处理时使用的数据  $u^*$  。 它与  @p ScratchData.  类似，但要简单得多。
  template <int dim> 
  struct HDG<dim>::PostProcessScratchData 
  { 
    FEValues<dim> fe_values_local; 
    FEValues<dim> fe_values; 

    std::vector<double>         u_values; 
    std::vector<Tensor<1, dim>> u_gradients; 
    FullMatrix<double>          cell_matrix; 

    Vector<double> cell_rhs; 
    Vector<double> cell_sol; 

    PostProcessScratchData(const FiniteElement<dim> &fe, 
                           const FiniteElement<dim> &fe_local, 
                           const QGauss<dim> &       quadrature_formula, 
                           const UpdateFlags         local_flags, 
                           const UpdateFlags         flags) 
      : fe_values_local(fe_local, quadrature_formula, local_flags) 
      , fe_values(fe, quadrature_formula, flags) 
      , u_values(quadrature_formula.size()) 
      , u_gradients(quadrature_formula.size()) 
      , cell_matrix(fe.n_dofs_per_cell(), fe.n_dofs_per_cell()) 
      , cell_rhs(fe.n_dofs_per_cell()) 
      , cell_sol(fe.n_dofs_per_cell()) 
    {} 

    PostProcessScratchData(const PostProcessScratchData &sd) 
      : fe_values_local(sd.fe_values_local.get_fe(), 
                        sd.fe_values_local.get_quadrature(), 
                        sd.fe_values_local.get_update_flags()) 
      , fe_values(sd.fe_values.get_fe(), 
                  sd.fe_values.get_quadrature(), 
                  sd.fe_values.get_update_flags()) 
      , u_values(sd.u_values) 
      , u_gradients(sd.u_gradients) 
      , cell_matrix(sd.cell_matrix) 
      , cell_rhs(sd.cell_rhs) 
      , cell_sol(sd.cell_sol) 
    {} 
  }; 

//  @sect4{HDG::assemble_system}   @p assemble_system 函数与 Step-32 上的函数类似，其中正交公式和更新标志被设置，然后 <code>WorkStream</code> 被用来以多线程的方式进行工作。  @p trace_reconstruct  输入参数用于决定我们是求全局骨架解（false）还是局部解（true）。

// 对于汇编的多线程执行，有一点值得注意的是，`assemble_system_one_cell()`中的局部计算会调用BLAS和LAPACK函数，如果这些函数在deal.II中可用。因此，底层的BLAS/LAPACK库必须支持同时来自多个线程的调用。大多数实现都支持这一点，但有些库需要以特定方式构建以避免问题。例如，在BLAS/LAPACK调用内部没有多线程的情况下编译的OpenBLAS需要在构建时将一个名为`USE_LOCKING'的标志设置为true。

  template <int dim> 
  void HDG<dim>::assemble_system(const bool trace_reconstruct) 
  { 
    const QGauss<dim>     quadrature_formula(fe.degree + 1); 
    const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    const UpdateFlags local_flags(update_values | update_gradients | 
                                  update_JxW_values | update_quadrature_points); 

    const UpdateFlags local_face_flags(update_values); 

    const UpdateFlags flags(update_values | update_normal_vectors | 
                            update_quadrature_points | update_JxW_values); 

    PerTaskData task_data(fe.n_dofs_per_cell(), trace_reconstruct); 
    ScratchData scratch(fe, 
                        fe_local, 
                        quadrature_formula, 
                        face_quadrature_formula, 
                        local_flags, 
                        local_face_flags, 
                        flags); 

    WorkStream::run(dof_handler.begin_active(), 
                    dof_handler.end(), 
                    *this, 
                    &HDG<dim>::assemble_system_one_cell, 
                    &HDG<dim>::copy_local_to_global, 
                    scratch, 
                    task_data); 
  } 

//  @sect4{HDG::assemble_system_one_cell}  HDG程序的实际工作由  @p assemble_system_one_cell.  组装局部矩阵  $A, B, C$  在这里完成，同时还有全局矩阵的局部贡献  $D$  。

  template <int dim> 
  void HDG<dim>::assemble_system_one_cell( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    ScratchData &                                         scratch, 
    PerTaskData &                                         task_data) 
  { 

//为Dof_handler_local构建迭代器，用于FEValues的reinit函数。

    typename DoFHandler<dim>::active_cell_iterator loc_cell(&triangulation, 
                                                            cell->level(), 
                                                            cell->index(), 
                                                            &dof_handler_local); 

    const unsigned int n_q_points = 
      scratch.fe_values_local.get_quadrature().size(); 
    const unsigned int n_face_q_points = 
      scratch.fe_face_values_local.get_quadrature().size(); 

    const unsigned int loc_dofs_per_cell = 
      scratch.fe_values_local.get_fe().n_dofs_per_cell(); 

    const FEValuesExtractors::Vector fluxes(0); 
    const FEValuesExtractors::Scalar scalar(dim); 

    scratch.ll_matrix = 0; 
    scratch.l_rhs     = 0; 
    if (!task_data.trace_reconstruct) 
      { 
        scratch.lf_matrix     = 0; 
        scratch.fl_matrix     = 0; 
        task_data.cell_matrix = 0; 
        task_data.cell_vector = 0; 
      } 
    scratch.fe_values_local.reinit(loc_cell); 

// 我们首先计算对应于局部-局部耦合的 @p ll_matrix 矩阵（在介绍中称为矩阵 $A$ ）的单元内部贡献，以及局部右手向量。 我们在每个正交点存储基函数、右手边值和对流速度的值，以便快速访问这些场。

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const double rhs_value = scratch.right_hand_side.value( 
          scratch.fe_values_local.quadrature_point(q)); 
        const Tensor<1, dim> convection = scratch.convection_velocity.value( 
          scratch.fe_values_local.quadrature_point(q)); 
        const double JxW = scratch.fe_values_local.JxW(q); 
        for (unsigned int k = 0; k < loc_dofs_per_cell; ++k) 
          { 
            scratch.q_phi[k] = scratch.fe_values_local[fluxes].value(k, q); 
            scratch.q_phi_div[k] = 
              scratch.fe_values_local[fluxes].divergence(k, q); 
            scratch.u_phi[k] = scratch.fe_values_local[scalar].value(k, q); 
            scratch.u_phi_grad[k] = 
              scratch.fe_values_local[scalar].gradient(k, q); 
          } 
        for (unsigned int i = 0; i < loc_dofs_per_cell; ++i) 
          { 
            for (unsigned int j = 0; j < loc_dofs_per_cell; ++j) 
              scratch.ll_matrix(i, j) += 
                (scratch.q_phi[i] * scratch.q_phi[j] - 
                 scratch.q_phi_div[i] * scratch.u_phi[j] + 
                 scratch.u_phi[i] * scratch.q_phi_div[j] - 
                 (scratch.u_phi_grad[i] * convection) * scratch.u_phi[j]) * 
                JxW; 
            scratch.l_rhs(i) += scratch.u_phi[i] * rhs_value * JxW; 
          } 
      } 

// 脸部条款是在所有元素的所有面上集合起来的。这与更传统的DG方法相反，在组装过程中，每个面只被访问一次。

    for (const auto face_no : cell->face_indices()) 
      { 
        scratch.fe_face_values_local.reinit(loc_cell, face_no); 
        scratch.fe_face_values.reinit(cell, face_no); 

// 在求解局部变量时需要已经得到的  $\hat{u}$  值。

        if (task_data.trace_reconstruct) 
          scratch.fe_face_values.get_function_values(solution, 
                                                     scratch.trace_values); 

        for (unsigned int q = 0; q < n_face_q_points; ++q) 
          { 
            const double     JxW = scratch.fe_face_values.JxW(q); 
            const Point<dim> quadrature_point = 
              scratch.fe_face_values.quadrature_point(q); 
            const Tensor<1, dim> normal = 
              scratch.fe_face_values.normal_vector(q); 
            const Tensor<1, dim> convection = 
              scratch.convection_velocity.value(quadrature_point); 

// 这里我们计算介绍中讨论的稳定参数：由于扩散是1，并且扩散长度尺度被设定为1/5，它只是导致扩散部分的贡献为5，而对流部分的贡献是通过元素边界的居中方案中的对流大小。

            const double tau_stab = (5. + std::abs(convection * normal)); 

// 我们存储非零通量和标量值，利用我们在 @p ScratchData. 中创建的 support_on_face 信息。
            for (unsigned int k = 0; 
                 k < scratch.fe_local_support_on_face[face_no].size(); 
                 ++k) 
              { 
                const unsigned int kk = 
                  scratch.fe_local_support_on_face[face_no][k]; 
                scratch.q_phi[k] = 
                  scratch.fe_face_values_local[fluxes].value(kk, q); 
                scratch.u_phi[k] = 
                  scratch.fe_face_values_local[scalar].value(kk, q); 
              } 

// 当  @p trace_reconstruct=false,  我们准备为骨架变量  $\hat{u}$  组装系统。如果是这种情况，我们必须组装所有与问题相关的局部矩阵：局部-局部、局部-面部、面部-局部和面部-面部。 面-面矩阵被存储为 @p TaskData::cell_matrix, ，这样就可以通过 @p copy_local_to_global将其组装到全局系统中。

            if (!task_data.trace_reconstruct) 
              { 
                for (unsigned int k = 0; 
                     k < scratch.fe_support_on_face[face_no].size(); 
                     ++k) 
                  scratch.tr_phi[k] = scratch.fe_face_values.shape_value( 
                    scratch.fe_support_on_face[face_no][k], q); 
                for (unsigned int i = 0; 
                     i < scratch.fe_local_support_on_face[face_no].size(); 
                     ++i) 
                  for (unsigned int j = 0; 
                       j < scratch.fe_support_on_face[face_no].size(); 
                       ++j) 
                    { 
                      const unsigned int ii = 
                        scratch.fe_local_support_on_face[face_no][i]; 
                      const unsigned int jj = 
                        scratch.fe_support_on_face[face_no][j]; 
                      scratch.lf_matrix(ii, jj) += 
                        ((scratch.q_phi[i] * normal + 
                          (convection * normal - tau_stab) * scratch.u_phi[i]) * 
                         scratch.tr_phi[j]) * 
                        JxW; 

// 注意face_no-local矩阵的符号。 我们在组装时否定了这个符号，这样我们就可以在计算舒尔补时使用 FullMatrix::mmult 的加法。

                      scratch.fl_matrix(jj, ii) -= 
                        ((scratch.q_phi[i] * normal + 
                          tau_stab * scratch.u_phi[i]) * 
                         scratch.tr_phi[j]) * 
                        JxW; 
                    } 

                for (unsigned int i = 0; 
                     i < scratch.fe_support_on_face[face_no].size(); 
                     ++i) 
                  for (unsigned int j = 0; 
                       j < scratch.fe_support_on_face[face_no].size(); 
                       ++j) 
                    { 
                      const unsigned int ii = 
                        scratch.fe_support_on_face[face_no][i]; 
                      const unsigned int jj = 
                        scratch.fe_support_on_face[face_no][j]; 
                      task_data.cell_matrix(ii, jj) += 
                        ((convection * normal - tau_stab) * scratch.tr_phi[i] * 
                         scratch.tr_phi[j]) * 
                        JxW; 
                    } 

                if (cell->face(face_no)->at_boundary() && 
                    (cell->face(face_no)->boundary_id() == 1)) 
                  { 
                    const double neumann_value = 
                      -scratch.exact_solution.gradient(quadrature_point) * 
                        normal + 
                      convection * normal * 
                        scratch.exact_solution.value(quadrature_point); 
                    for (unsigned int i = 0; 
                         i < scratch.fe_support_on_face[face_no].size(); 
                         ++i) 
                      { 
                        const unsigned int ii = 
                          scratch.fe_support_on_face[face_no][i]; 
                        task_data.cell_vector(ii) += 
                          scratch.tr_phi[i] * neumann_value * JxW; 
                      } 
                  } 
              } 

// 这最后一个项将 $\left<w,\tau u_h\right>_{\partial \mathcal T}$ 项的贡献加入到本地矩阵中。相对于上面的脸部矩阵，我们在两个装配阶段都需要它。

            for (unsigned int i = 0; 
                 i < scratch.fe_local_support_on_face[face_no].size(); 
                 ++i) 
              for (unsigned int j = 0; 
                   j < scratch.fe_local_support_on_face[face_no].size(); 
                   ++j) 
                { 
                  const unsigned int ii = 
                    scratch.fe_local_support_on_face[face_no][i]; 
                  const unsigned int jj = 
                    scratch.fe_local_support_on_face[face_no][j]; 
                  scratch.ll_matrix(ii, jj) += 
                    tau_stab * scratch.u_phi[i] * scratch.u_phi[j] * JxW; 
                } 

// 当 @p trace_reconstruct=true, 时，我们在逐个元素的基础上求解局部解。 局部右手边的计算是通过用计算值 @p trace_values替换 @p 计算中的基函数 @p tr_phi。 当然，现在矩阵的符号是减号，因为我们已经把所有的东西移到了方程的另一边。

            if (task_data.trace_reconstruct) 
              for (unsigned int i = 0; 
                   i < scratch.fe_local_support_on_face[face_no].size(); 
                   ++i) 
                { 
                  const unsigned int ii = 
                    scratch.fe_local_support_on_face[face_no][i]; 
                  scratch.l_rhs(ii) -= 
                    (scratch.q_phi[i] * normal + 
                     scratch.u_phi[i] * (convection * normal - tau_stab)) * 
                    scratch.trace_values[q] * JxW; 
                } 
          } 
      } 

// 一旦完成所有局部贡献的组装，我们必须：（1）组装全局系统；（2）计算局部贡献。(1)组装全局系统，或者(2)计算局部解值并保存。无论哪种情况，第一步都是对局部-局部矩阵进行反转。

    scratch.ll_matrix.gauss_jordan(); 

// 对于(1)，我们计算舒尔补码，并将其添加到 @p  cell_matrix，介绍中的矩阵 $D$ 。

    if (task_data.trace_reconstruct == false) 
      { 
        scratch.fl_matrix.mmult(scratch.tmp_matrix, scratch.ll_matrix); 
        scratch.tmp_matrix.vmult_add(task_data.cell_vector, scratch.l_rhs); 
        scratch.tmp_matrix.mmult(task_data.cell_matrix, 
                                 scratch.lf_matrix, 
                                 true); 
        cell->get_dof_indices(task_data.dof_indices); 
      } 

// 对于(2)，我们只是求解(ll_matrix). (solution_local) = (l_rhs)。因此，我们用 @p l_rhs 乘以我们已经倒置的局部-局部矩阵，并用 <code>set_dof_values</code> 函数来存储结果。

    else 
      { 
        scratch.ll_matrix.vmult(scratch.tmp_rhs, scratch.l_rhs); 
        loc_cell->set_dof_values(scratch.tmp_rhs, solution_local); 
      } 
  } 

// 如果我们处于解题的第一步，即 @sect4{HDG::copy_local_to_global} ，那么我们就把局部矩阵组装到全局系统中。

  template <int dim> 
  void HDG<dim>::copy_local_to_global(const PerTaskData &data) 
  { 
    if (data.trace_reconstruct == false) 
      constraints.distribute_local_to_global(data.cell_matrix, 
                                             data.cell_vector, 
                                             data.dof_indices, 
                                             system_matrix, 
                                             system_rhs); 
  } 

//  @sect4{HDG::solve}  骨架解是通过使用带有身份预处理程序的BiCGStab求解器来解决的。

  template <int dim> 
  void HDG<dim>::solve() 
  { 
    SolverControl                  solver_control(system_matrix.m() * 10, 
                                 1e-11 * system_rhs.l2_norm()); 
    SolverBicgstab<Vector<double>> solver(solver_control); 
    solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity()); 

    std::cout << "   Number of BiCGStab iterations: " 
              << solver_control.last_step() << std::endl; 

    system_matrix.clear(); 
    sparsity_pattern.reinit(0, 0, 0, 1); 

    constraints.distribute(solution); 

// 一旦我们求出了骨架解，我们就可以以逐个元素的方式求出局部解。 我们通过重新使用相同的 @p assemble_system 函数来做到这一点，但将 @p trace_reconstruct 切换为真。

    assemble_system(true); 
  } 

//  @sect4{HDG::postprocess}  

// 后处理方法有两个目的。首先，我们要在度数为 $p+1$ 的元素空间中构造一个后处理的标量变量，我们希望它能在阶 $p+2$ 上收敛。这也是一个逐个元素的过程，只涉及标量解以及局部单元上的梯度。为了做到这一点，我们引入了已经定义好的从头开始的数据以及一些更新标志，并运行工作流来并行地完成这一工作。

// 第二，我们要计算离散化误差，就像我们在  step-7  中做的那样。整个过程与调用 VectorTools::integrate_difference. 相似，区别在于我们如何计算标量变量和梯度变量的误差。在 step-7 中，我们通过计算 @p L2_norm 或 @p H1_seminorm 的贡献来做到这一点。在这里，我们有一个DoFHandler，计算了这两个贡献，并按其矢量分量排序， <code>[0, dim)</code> 为梯度， @p dim 为标量。为了计算它们的值，我们用一个ComponentSelectFunction来计算它们中的任何一个，再加上上面介绍的 @p SolutionAndGradient类，它包含了它们中任何一个的分析部分。最终，我们还计算了后处理的解决方案的L2-误差，并将结果添加到收敛表中。

  template <int dim> 
  void HDG<dim>::postprocess() 
  { 
    { 
      const QGauss<dim> quadrature_formula(fe_u_post.degree + 1); 
      const UpdateFlags local_flags(update_values); 
      const UpdateFlags flags(update_values | update_gradients | 
                              update_JxW_values); 

      PostProcessScratchData scratch( 
        fe_u_post, fe_local, quadrature_formula, local_flags, flags); 

      WorkStream::run( 
        dof_handler_u_post.begin_active(), 
        dof_handler_u_post.end(), 
        [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
               PostProcessScratchData &                              scratch, 
               unsigned int &                                        data) { 
          this->postprocess_one_cell(cell, scratch, data); 
        }, 
        std::function<void(const unsigned int &)>(), 
        scratch, 
        0U); 
    } 

    Vector<float> difference_per_cell(triangulation.n_active_cells()); 

    ComponentSelectFunction<dim> value_select(dim, dim + 1); 
    VectorTools::integrate_difference(dof_handler_local, 
                                      solution_local, 
                                      SolutionAndGradient<dim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm, 
                                      &value_select); 
    const double L2_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 

    ComponentSelectFunction<dim> gradient_select( 
      std::pair<unsigned int, unsigned int>(0, dim), dim + 1); 
    VectorTools::integrate_difference(dof_handler_local, 
                                      solution_local, 
                                      SolutionAndGradient<dim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 2), 
                                      VectorTools::L2_norm, 
                                      &gradient_select); 
    const double grad_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 

    VectorTools::integrate_difference(dof_handler_u_post, 
                                      solution_u_post, 
                                      Solution<dim>(), 
                                      difference_per_cell, 
                                      QGauss<dim>(fe.degree + 3), 
                                      VectorTools::L2_norm); 
    const double post_error = 
      VectorTools::compute_global_error(triangulation, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 

    convergence_table.add_value("cells", triangulation.n_active_cells()); 
    convergence_table.add_value("dofs", dof_handler.n_dofs()); 

    convergence_table.add_value("val L2", L2_error); 
    convergence_table.set_scientific("val L2", true); 
    convergence_table.set_precision("val L2", 3); 

    convergence_table.add_value("grad L2", grad_error); 
    convergence_table.set_scientific("grad L2", true); 
    convergence_table.set_precision("grad L2", 3); 

    convergence_table.add_value("val L2-post", post_error); 
    convergence_table.set_scientific("val L2-post", true); 
    convergence_table.set_precision("val L2-post", 3); 
  } 

//  @sect4{HDG::postprocess_one_cell}  

// 这是为后处理所做的实际工作。根据介绍中的讨论，我们需要建立一个系统，将DG解的梯度部分投影到后处理变量的梯度上。此外，我们还需要将新的后处理变量的平均值设置为等于标量DG解在单元上的平均值。

// 从技术上讲，梯度的投影是一个有可能填满我们的 @p dofs_per_cell 乘以 @p dofs_per_cell 矩阵的系统，但它是单数（所有行的总和为零，因为常数函数的梯度为零）。因此，我们拿掉一行，用它来强加标量值的平均值。我们为标量部分挑选第一行，尽管我们可以为 $\mathcal Q_{-p}$ 元素挑选任何一行。然而，如果我们使用FE_DGP元素，第一行将对应常数部分，删除例如最后一行将得到一个奇异系统。这样一来，我们的程序也可以用于这些元素。

  template <int dim> 
  void HDG<dim>::postprocess_one_cell( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    PostProcessScratchData &                              scratch, 
    unsigned int &) 
  { 
    typename DoFHandler<dim>::active_cell_iterator loc_cell(&triangulation, 
                                                            cell->level(), 
                                                            cell->index(), 
                                                            &dof_handler_local); 

    scratch.fe_values_local.reinit(loc_cell); 
    scratch.fe_values.reinit(cell); 

    FEValuesExtractors::Vector fluxes(0); 
    FEValuesExtractors::Scalar scalar(dim); 

    const unsigned int n_q_points = scratch.fe_values.get_quadrature().size(); 
    const unsigned int dofs_per_cell = scratch.fe_values.dofs_per_cell; 

    scratch.fe_values_local[scalar].get_function_values(solution_local, 
                                                        scratch.u_values); 
    scratch.fe_values_local[fluxes].get_function_values(solution_local, 
                                                        scratch.u_gradients); 

    double sum = 0; 
    for (unsigned int i = 1; i < dofs_per_cell; ++i) 
      { 
        for (unsigned int j = 0; j < dofs_per_cell; ++j) 
          { 
            sum = 0; 
            for (unsigned int q = 0; q < n_q_points; ++q) 
              sum += (scratch.fe_values.shape_grad(i, q) * 
                      scratch.fe_values.shape_grad(j, q)) * 
                     scratch.fe_values.JxW(q); 
            scratch.cell_matrix(i, j) = sum; 
          } 

        sum = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          sum -= (scratch.fe_values.shape_grad(i, q) * scratch.u_gradients[q]) * 
                 scratch.fe_values.JxW(q); 
        scratch.cell_rhs(i) = sum; 
      } 
    for (unsigned int j = 0; j < dofs_per_cell; ++j) 
      { 
        sum = 0; 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          sum += scratch.fe_values.shape_value(j, q) * scratch.fe_values.JxW(q); 
        scratch.cell_matrix(0, j) = sum; 
      } 
    { 
      sum = 0; 
      for (unsigned int q = 0; q < n_q_points; ++q) 
        sum += scratch.u_values[q] * scratch.fe_values.JxW(q); 
      scratch.cell_rhs(0) = sum; 
    } 

// 集合了所有条款后，我们又可以继续解决这个线性系统。我们对矩阵进行反转，然后将反转结果乘以右手边。另一种方法（数字上更稳定）是只对矩阵进行因式分解，然后应用因式分解。

    scratch.cell_matrix.gauss_jordan(); 
    scratch.cell_matrix.vmult(scratch.cell_sol, scratch.cell_rhs); 
    cell->distribute_local_to_global(scratch.cell_sol, solution_u_post); 
  } 

//  @sect4{HDG::output_results}  我们有三组我们想输出的结果：局部解决方案，后处理的局部解决方案，以及骨架解决方案。前两个结果都 "活 "在元素体积上，而后者则活在三角形的一维表面上。 我们的 @p output_results 函数将所有的局部解决方案写入同一个vtk文件，尽管它们对应于不同的DoFHandler对象。 骨架变量的图形输出是通过使用DataOutFaces类完成的。

  template <int dim> 
  void HDG<dim>::output_results(const unsigned int cycle) 
  { 
    std::string filename; 
    switch (refinement_mode) 
      { 
        case global_refinement: 
          filename = "solution-global"; 
          break; 
        case adaptive_refinement: 
          filename = "solution-adaptive"; 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    std::string face_out(filename); 
    face_out += "-face"; 

    filename += "-q" + Utilities::int_to_string(fe.degree, 1); 
    filename += "-" + Utilities::int_to_string(cycle, 2); 
    filename += ".vtk"; 
    std::ofstream output(filename); 

    DataOut<dim> data_out; 

// 我们首先定义本地解决方案的名称和类型，并将数据添加到  @p data_out.  中。
    std::vector<std::string> names(dim, "gradient"); 
    names.emplace_back("solution"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      component_interpretation( 
        dim + 1, DataComponentInterpretation::component_is_part_of_vector); 
    component_interpretation[dim] = 
      DataComponentInterpretation::component_is_scalar; 
    data_out.add_data_vector(dof_handler_local, 
                             solution_local, 
                             names, 
                             component_interpretation); 

// 我们添加的第二个数据项是后处理的解决方案。在这种情况下，它是一个属于不同DoFHandler的单一标量变量。

    std::vector<std::string> post_name(1, "u_post"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      post_comp_type(1, DataComponentInterpretation::component_is_scalar); 
    data_out.add_data_vector(dof_handler_u_post, 
                             solution_u_post, 
                             post_name, 
                             post_comp_type); 

    data_out.build_patches(fe.degree); 
    data_out.write_vtk(output); 

    face_out += "-q" + Utilities::int_to_string(fe.degree, 1); 
    face_out += "-" + Utilities::int_to_string(cycle, 2); 
    face_out += ".vtk"; 
    std::ofstream face_output(face_out); 

//  <code>DataOutFaces</code> 类的工作原理与 <code>DataOut</code> class when we have a <code>DoFHandler</code> 类似，后者定义了三角形骨架上的解决方案。 我们在此将其视为如此，代码与上面类似。

    DataOutFaces<dim>        data_out_face(false); 
    std::vector<std::string> face_name(1, "u_hat"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      face_component_type(1, DataComponentInterpretation::component_is_scalar); 

    data_out_face.add_data_vector(dof_handler, 
                                  solution, 
                                  face_name, 
                                  face_component_type); 

    data_out_face.build_patches(fe.degree); 
    data_out_face.write_vtk(face_output); 
  } 
// @sect4{HDG::refine_grid}  

// 我们为HDG实现了两种不同的细化情况，就像在 <code>Step-7</code> 中一样：adaptive_refinement和global_refinement。 global_refinement选项每次都会重新创建整个三角形。这是因为我们想使用比一个细化步骤更细的网格序列，即每个方向2、3、4、6、8、12、16...个元素。

// adaptive_refinement模式使用 <code>KellyErrorEstimator</code> 对标量局部解中的非规则区域给出一个体面的指示。

  template <int dim> 
  void HDG<dim>::refine_grid(const unsigned int cycle) 
  { 
    if (cycle == 0) 
      { 
        GridGenerator::subdivided_hyper_cube(triangulation, 2, -1, 1); 
        triangulation.refine_global(3 - dim); 
      } 
    else 
      switch (refinement_mode) 
        { 
          case global_refinement: 
            { 
              triangulation.clear(); 
              GridGenerator::subdivided_hyper_cube(triangulation, 
                                                   2 + (cycle % 2), 
                                                   -1, 
                                                   1); 
              triangulation.refine_global(3 - dim + cycle / 2); 
              break; 
            } 

          case adaptive_refinement: 
            { 
              Vector<float> estimated_error_per_cell( 
                triangulation.n_active_cells()); 

              FEValuesExtractors::Scalar scalar(dim); 
              std::map<types::boundary_id, const Function<dim> *> 
                neumann_boundary; 
              KellyErrorEstimator<dim>::estimate(dof_handler_local, 
                                                 QGauss<dim - 1>(fe.degree + 1), 
                                                 neumann_boundary, 
                                                 solution_local, 
                                                 estimated_error_per_cell, 
                                                 fe_local.component_mask( 
                                                   scalar)); 

              GridRefinement::refine_and_coarsen_fixed_number( 
                triangulation, estimated_error_per_cell, 0.3, 0.); 

              triangulation.execute_coarsening_and_refinement(); 

              break; 
            } 

          default: 
            { 
              Assert(false, ExcNotImplemented()); 
            } 
        } 

// 就像在 step-7 中一样，我们将其中两个面的边界指标设置为1，在这里我们要指定诺伊曼边界条件而不是迪里希特条件。由于我们每次都会为全局细化重新创建三角形，所以在每个细化步骤中都会设置标志，而不仅仅是在开始时。

    for (const auto &cell : triangulation.cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          if ((std::fabs(face->center()(0) - (-1)) < 1e-12) || 
              (std::fabs(face->center()(1) - (-1)) < 1e-12)) 
            face->set_boundary_id(1); 
  } 
// @sect4{HDG::run}  这里的功能与 <code>Step-7</code>  基本相同。我们在10个周期中循环，在每个周期中细化网格。 在最后，收敛表被创建。

  template <int dim> 
  void HDG<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 10; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << ':' << std::endl; 

        refine_grid(cycle); 
        setup_system(); 
        assemble_system(false); 
        solve(); 
        postprocess(); 
        output_results(cycle); 
      } 

// 与 step-7 相比，收敛表有一个微小的变化：由于我们没有在每个周期内以2的系数细化我们的网格（而是使用2，3，4，6，8，12，...的序列），我们需要告诉收敛率评估这一点。我们通过设置单元格数量作为参考列，并额外指定问题的维度来实现这一目的，这为单元格数量和网格大小之间的关系提供了必要的信息。

    if (refinement_mode == global_refinement) 
      { 
        convergence_table.evaluate_convergence_rates( 
          "val L2", "cells", ConvergenceTable::reduction_rate_log2, dim); 
        convergence_table.evaluate_convergence_rates( 
          "grad L2", "cells", ConvergenceTable::reduction_rate_log2, dim); 
        convergence_table.evaluate_convergence_rates( 
          "val L2-post", "cells", ConvergenceTable::reduction_rate_log2, dim); 
      } 
    convergence_table.write_text(std::cout); 
  } 

} // end of namespace Step51 

int main() 
{ 
  const unsigned int dim = 2; 

  try 
    { 

// 现在是对主类的三次调用，完全类似于  step-7  。

      { 
        std::cout << "Solving with Q1 elements, adaptive refinement" 
                  << std::endl 
                  << "=============================================" 
                  << std::endl 
                  << std::endl; 

        Step51::HDG<dim> hdg_problem(1, Step51::HDG<dim>::adaptive_refinement); 
        hdg_problem.run(); 

        std::cout << std::endl; 
      } 

      { 
        std::cout << "Solving with Q1 elements, global refinement" << std::endl 
                  << "===========================================" << std::endl 
                  << std::endl; 

        Step51::HDG<dim> hdg_problem(1, Step51::HDG<dim>::global_refinement); 
        hdg_problem.run(); 

        std::cout << std::endl; 
      } 

      { 
        std::cout << "Solving with Q3 elements, global refinement" << std::endl 
                  << "===========================================" << std::endl 
                  << std::endl; 

        Step51::HDG<dim> hdg_problem(3, Step51::HDG<dim>::global_refinement); 
        hdg_problem.run(); 

        std::cout << std::endl; 
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
    } 

  return 0; 
} 

