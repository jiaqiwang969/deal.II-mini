

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2010 - 2021 by the deal.II authors 
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
 * Authors: Chih-Che Chueh, University of Victoria, 2010 
 *          Wolfgang Bangerth, Texas A&M University, 2010 
 */ 


// @sect3{Include files}  

// 像往常一样，第一步是包括一些deal.II和C++头文件的功能。

// 列表中包括一些提供向量、矩阵和预处理类的头文件，这些头文件实现了各自Trilinos类的接口；关于这些的一些更多信息可以在  step-31  中找到。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/tensor_function.h> 
#include <deal.II/base/index_set.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/block_sparsity_pattern.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/solution_transfer.h> 

#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_block_sparse_matrix.h> 
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_parallel_block_vector.h> 
#include <deal.II/lac/trilinos_precondition.h> 

#include <iostream> 
#include <fstream> 
#include <memory> 

// 在这个顶层设计的最后，我们为当前项目开辟一个命名空间，下面的所有材料都将进入这个命名空间，然后将所有deal.II名称导入这个命名空间。

namespace Step43 
{ 
  using namespace dealii; 
// @sect3{Boundary and initial value classes}  

// 下面的部分直接取自 step-21 ，所以没有必要重复那里的描述。

  template <int dim> 
  class PressureBoundaryValues : public Function<dim> 
  { 
  public: 
    PressureBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double 
  PressureBoundaryValues<dim>::value(const Point<dim> &p, 
                                     const unsigned int /*component*/) const 
  { 
    return 1 - p[0]; 
  } 

  template <int dim> 
  class SaturationBoundaryValues : public Function<dim> 
  { 
  public: 
    SaturationBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double 
  SaturationBoundaryValues<dim>::value(const Point<dim> &p, 
                                       const unsigned int /*component*/) const 
  { 
    if (p[0] == 0) 
      return 1; 
    else 
      return 0; 
  } 

  template <int dim> 
  class SaturationInitialValues : public Function<dim> 
  { 
  public: 
    SaturationInitialValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  double 
  SaturationInitialValues<dim>::value(const Point<dim> & /*p*/, 
                                      const unsigned int /*component*/) const 
  { 
    return 0.2; 
  } 

  template <int dim> 
  void SaturationInitialValues<dim>::vector_value(const Point<dim> &p, 
                                                  Vector<double> &values) const 
  { 
    for (unsigned int c = 0; c < this->n_components; ++c) 
      values(c) = SaturationInitialValues<dim>::value(p, c); 
  } 
// @sect3{Permeability models}  

// 在本教程中，我们仍然使用之前在 step-21 中使用的两个渗透率模型，所以我们再次避免对它们进行详细评论。

  namespace SingleCurvingCrack 
  { 
    template <int dim> 
    class KInverse : public TensorFunction<2, dim> 
    { 
    public: 
      KInverse() 
        : TensorFunction<2, dim>() 
      {} 

      virtual void 
      value_list(const std::vector<Point<dim>> &points, 
                 std::vector<Tensor<2, dim>> &  values) const override; 
    }; 

    template <int dim> 
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
                                   std::vector<Tensor<2, dim>> &  values) const 
    { 
      Assert(points.size() == values.size(), 
             ExcDimensionMismatch(points.size(), values.size())); 

      for (unsigned int p = 0; p < points.size(); ++p) 
        { 
          values[p].clear(); 

          const double distance_to_flowline = 
            std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0])); 

          const double permeability = 
            std::max(std::exp(-(distance_to_flowline * distance_to_flowline) / 
                              (0.1 * 0.1)), 
                     0.01); 

          for (unsigned int d = 0; d < dim; ++d) 
            values[p][d][d] = 1. / permeability; 
        } 
    } 
  } // namespace SingleCurvingCrack 

  namespace RandomMedium 
  { 
    template <int dim> 
    class KInverse : public TensorFunction<2, dim> 
    { 
    public: 
      KInverse() 
        : TensorFunction<2, dim>() 
      {} 

      virtual void 
      value_list(const std::vector<Point<dim>> &points, 
                 std::vector<Tensor<2, dim>> &  values) const override; 

    private: 
      static std::vector<Point<dim>> centers; 
    }; 

    template <int dim> 
    std::vector<Point<dim>> KInverse<dim>::centers = []() { 
      const unsigned int N = 
        (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented())); 

      std::vector<Point<dim>> centers_list(N); 
      for (unsigned int i = 0; i < N; ++i) 
        for (unsigned int d = 0; d < dim; ++d) 
          centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX; 

      return centers_list; 
    }(); 

    template <int dim> 
    void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
                                   std::vector<Tensor<2, dim>> &  values) const 
    { 
      AssertDimension(points.size(), values.size()); 

      for (unsigned int p = 0; p < points.size(); ++p) 
        { 
          values[p].clear(); 

          double permeability = 0; 
          for (unsigned int i = 0; i < centers.size(); ++i) 
            permeability += 
              std::exp(-(points[p] - centers[i]).norm_square() / (0.05 * 0.05)); 

          const double normalized_permeability = 
            std::min(std::max(permeability, 0.01), 4.); 

          for (unsigned int d = 0; d < dim; ++d) 
            values[p][d][d] = 1. / normalized_permeability; 
        } 
    } 
  } // namespace RandomMedium 
// @sect3{Physical quantities}  

// 所有物理量的实现，如总流动性 $\lambda_t$ 和水的部分流量 $F$ 都来自 step-21 ，所以我们也没有对它们做任何评论。与 step-21 相比，我们增加了检查，即传递给这些函数的饱和度实际上是在物理上有效的范围内。此外，鉴于润湿相以速度 $\mathbf u F'(S)$ 移动，很明显 $F'(S)$ 必须大于或等于零，所以我们也断言，以确保我们的计算得到的导数公式是合理的。

  double mobility_inverse(const double S, const double viscosity) 
  { 
    return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S)); 
  } 

  double fractional_flow(const double S, const double viscosity) 
  { 
    Assert((S >= 0) && (S <= 1), 
           ExcMessage("Saturation is outside its physically valid range.")); 

    return S * S / (S * S + viscosity * (1 - S) * (1 - S)); 
  } 

  double fractional_flow_derivative(const double S, const double viscosity) 
  { 
    Assert((S >= 0) && (S <= 1), 
           ExcMessage("Saturation is outside its physically valid range.")); 

    const double temp = (S * S + viscosity * (1 - S) * (1 - S)); 

    const double numerator = 
      2.0 * S * temp - S * S * (2.0 * S - 2.0 * viscosity * (1 - S)); 
    const double denominator = std::pow(temp, 2.0); 

    const double F_prime = numerator / denominator; 

    Assert(F_prime >= 0, ExcInternalError()); 

    return F_prime; 
  } 
// @sect3{Helper classes for solvers and preconditioners}  

// 在这第一部分中，我们定义了一些我们在构建线性求解器和预处理器时需要的类。这一部分与  step-31  中使用的基本相同。唯一不同的是，原来的变量名称stokes_matrix被另一个名称darcy_matrix取代，以配合我们的问题。

  namespace LinearSolvers 
  { 
    template <class MatrixType, class PreconditionerType> 
    class InverseMatrix : public Subscriptor 
    { 
    public: 
      InverseMatrix(const MatrixType &        m, 
                    const PreconditionerType &preconditioner); 

      template <typename VectorType> 
      void vmult(VectorType &dst, const VectorType &src) const; 

    private: 
      const SmartPointer<const MatrixType> matrix; 
      const PreconditionerType &           preconditioner; 
    }; 

    template <class MatrixType, class PreconditionerType> 
    InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
      const MatrixType &        m, 
      const PreconditionerType &preconditioner) 
      : matrix(&m) 
      , preconditioner(preconditioner) 
    {} 

    template <class MatrixType, class PreconditionerType> 
    template <typename VectorType> 
    void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
      VectorType &      dst, 
      const VectorType &src) const 
    { 
      SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm()); 
      SolverCG<VectorType> cg(solver_control); 

      dst = 0; 

      try 
        { 
          cg.solve(*matrix, dst, src, preconditioner); 
        } 
      catch (std::exception &e) 
        { 
          Assert(false, ExcMessage(e.what())); 
        } 
    } 

    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    class BlockSchurPreconditioner : public Subscriptor 
    { 
    public: 
      BlockSchurPreconditioner( 
        const TrilinosWrappers::BlockSparseMatrix &S, 
        const InverseMatrix<TrilinosWrappers::SparseMatrix, 
                            PreconditionerTypeMp> &Mpinv, 
        const PreconditionerTypeA &                Apreconditioner); 

      void vmult(TrilinosWrappers::MPI::BlockVector &      dst, 
                 const TrilinosWrappers::MPI::BlockVector &src) const; 

    private: 
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> 
        darcy_matrix; 
      const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix, 
                                             PreconditionerTypeMp>> 
                                 m_inverse; 
      const PreconditionerTypeA &a_preconditioner; 

      mutable TrilinosWrappers::MPI::Vector tmp; 
    }; 

    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>:: 
      BlockSchurPreconditioner( 
        const TrilinosWrappers::BlockSparseMatrix &S, 
        const InverseMatrix<TrilinosWrappers::SparseMatrix, 
                            PreconditionerTypeMp> &Mpinv, 
        const PreconditionerTypeA &                Apreconditioner) 
      : darcy_matrix(&S) 
      , m_inverse(&Mpinv) 
      , a_preconditioner(Apreconditioner) 
      , tmp(complete_index_set(darcy_matrix->block(1, 1).m())) 
    {} 

    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    void 
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult( 
      TrilinosWrappers::MPI::BlockVector &      dst, 
      const TrilinosWrappers::MPI::BlockVector &src) const 
    { 
      a_preconditioner.vmult(dst.block(0), src.block(0)); 
      darcy_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1)); 
      tmp *= -1; 
      m_inverse->vmult(dst.block(1), tmp); 
    } 
  } // namespace LinearSolvers 
// @sect3{The TwoPhaseFlowProblem class}  

// 定义解决随时间变化的平流主导的两相流问题（或Buckley-Leverett问题[Buckley 1942]）的顶层逻辑的类的定义主要基于教程程序 step-21 和 step-33 ，特别是 step-31 ，我们在这里使用的一般结构基本相同。与 step-31 一样，在下面的实现中需要寻找的关键例程是 <code>run()</code> and <code>solve()</code> 函数。

// 与 step-31 的主要区别是，由于考虑了自适应算子拆分，我们需要多几个成员变量来保存最近两次计算的达西（速度/压力）解，以及当前的达西（直接计算，或从前两次计算中推断），我们需要记住最近两次计算的达西解。我们还需要一个辅助函数来确定我们是否真的需要重新计算达西解。

// 与 step-31 不同，这一步多用了一个AffineConstraints对象，叫做darcy_preconditioner_constraints。这个约束对象只用于为Darcy预处理程序组装矩阵，包括悬挂节点约束以及压力变量的Dirichlet边界值约束。我们需要这个，因为我们正在为压力建立一个拉普拉斯矩阵，作为舒尔补码的近似值），如果应用边界条件，这个矩阵是正定的。

// 这样在这个类中声明的成员函数和变量的集合与  step-31  中的相当相似。

  template <int dim> 
  class TwoPhaseFlowProblem 
  { 
  public: 
    TwoPhaseFlowProblem(const unsigned int degree); 
    void run(); 

  private: 
    void setup_dofs(); 
    void assemble_darcy_preconditioner(); 
    void build_darcy_preconditioner(); 
    void assemble_darcy_system(); 
    void assemble_saturation_system(); 
    void assemble_saturation_matrix(); 
    void assemble_saturation_rhs(); 
    void assemble_saturation_rhs_cell_term( 
      const FEValues<dim> &                       saturation_fe_values, 
      const FEValues<dim> &                       darcy_fe_values, 
      const double                                global_max_u_F_prime, 
      const double                                global_S_variation, 
      const std::vector<types::global_dof_index> &local_dof_indices); 
    void assemble_saturation_rhs_boundary_term( 
      const FEFaceValues<dim> &                   saturation_fe_face_values, 
      const FEFaceValues<dim> &                   darcy_fe_face_values, 
      const std::vector<types::global_dof_index> &local_dof_indices); 
    void solve(); 
    void refine_mesh(const unsigned int min_grid_level, 
                     const unsigned int max_grid_level); 
    void output_results() const; 

// 我们接下来会有一些辅助函数，这些函数在整个程序中的不同地方都会用到。

    double                    get_max_u_F_prime() const; 
    std::pair<double, double> get_extrapolated_saturation_range() const; 
    bool   determine_whether_to_solve_for_pressure_and_velocity() const; 
    void   project_back_saturation(); 
    double compute_viscosity( 
      const std::vector<double> &        old_saturation, 
      const std::vector<double> &        old_old_saturation, 
      const std::vector<Tensor<1, dim>> &old_saturation_grads, 
      const std::vector<Tensor<1, dim>> &old_old_saturation_grads, 
      const std::vector<Vector<double>> &present_darcy_values, 
      const double                       global_max_u_F_prime, 
      const double                       global_S_variation, 
      const double                       cell_diameter) const; 

// 接下来是成员变量，其中大部分与 step-31 中的变量类似，但与速度/压力系统的宏观时间步长有关的变量除外。

    Triangulation<dim> triangulation; 
    double             global_Omega_diameter; 

    const unsigned int degree; 

    const unsigned int        darcy_degree; 
    FESystem<dim>             darcy_fe; 
    DoFHandler<dim>           darcy_dof_handler; 
    AffineConstraints<double> darcy_constraints; 

    AffineConstraints<double> darcy_preconditioner_constraints; 

    TrilinosWrappers::BlockSparseMatrix darcy_matrix; 
    TrilinosWrappers::BlockSparseMatrix darcy_preconditioner_matrix; 

    TrilinosWrappers::MPI::BlockVector darcy_solution; 
    TrilinosWrappers::MPI::BlockVector darcy_rhs; 

    TrilinosWrappers::MPI::BlockVector last_computed_darcy_solution; 
    TrilinosWrappers::MPI::BlockVector second_last_computed_darcy_solution; 

    const unsigned int        saturation_degree; 
    FE_Q<dim>                 saturation_fe; 
    DoFHandler<dim>           saturation_dof_handler; 
    AffineConstraints<double> saturation_constraints; 

    TrilinosWrappers::SparseMatrix saturation_matrix; 

    TrilinosWrappers::MPI::Vector saturation_solution; 
    TrilinosWrappers::MPI::Vector old_saturation_solution; 
    TrilinosWrappers::MPI::Vector old_old_saturation_solution; 
    TrilinosWrappers::MPI::Vector saturation_rhs; 

    TrilinosWrappers::MPI::Vector 
      saturation_matching_last_computed_darcy_solution; 

    const double saturation_refinement_threshold; 

    double       time; 
    const double end_time; 

    double current_macro_time_step; 
    double old_macro_time_step; 

    double       time_step; 
    double       old_time_step; 
    unsigned int timestep_number; 

    const double viscosity; 
    const double porosity; 
    const double AOS_threshold; 

    std::shared_ptr<TrilinosWrappers::PreconditionIC> Amg_preconditioner; 
    std::shared_ptr<TrilinosWrappers::PreconditionIC> Mp_preconditioner; 

    bool rebuild_saturation_matrix; 

// 在最后，我们声明一个变量，表示材料模型。与 step-21 相比，我们在这里把它作为一个成员变量，因为我们想在不同的地方使用它，所以有一个声明这样一个变量的中心位置，将使我们更容易用另一个类来替换 RandomMedium::KInverse （例如，用 SingleCurvingCrack::KInverse). 替换 RandomMedium::KInverse ）。
    const RandomMedium::KInverse<dim> k_inverse; 
  }; 
// @sect3{TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem}  

// 这个类的构造函数是对  step-21  和  step-31  中的构造函数的扩展。我们需要添加涉及饱和度的各种变量。正如介绍中所讨论的，我们将再次使用 $Q_2 \times Q_1$ （Taylor-Hood）元素来处理Darcy系统，这是一个满足Ladyzhenskaya-Babuska-Brezzi（LBB）条件的元素组合[Brezzi and Fortin 1991, Chen 2005]，并使用 $Q_1$ 元素处理饱和度。然而，通过使用存储Darcy和温度有限元的多项式程度的变量，可以很容易地持续修改这些元素的程度以及在其上使用的所有正交公式的下游。此外，我们还初始化了与算子分割有关的时间步进变量，以及矩阵装配和预处理的选项。

  template <int dim> 
  TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree) 
    : triangulation(Triangulation<dim>::maximum_smoothing) 
    , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()) 
    , degree(degree) 
    , darcy_degree(degree) 
    , darcy_fe(FE_Q<dim>(darcy_degree + 1), dim, FE_Q<dim>(darcy_degree), 1) 
    , darcy_dof_handler(triangulation) 
    , 

    saturation_degree(degree + 1) 
    , saturation_fe(saturation_degree) 
    , saturation_dof_handler(triangulation) 
    , 

    saturation_refinement_threshold(0.5) 
    , 

    time(0) 
    , end_time(10) 
    , 

    current_macro_time_step(0) 
    , old_macro_time_step(0) 
    , 

    time_step(0) 
    , old_time_step(0) 
    , timestep_number(0) 
    , viscosity(0.2) 
    , porosity(1.0) 
    , AOS_threshold(3.0) 
    , 

    rebuild_saturation_matrix(true) 
  {} 
// @sect3{TwoPhaseFlowProblem<dim>::setup_dofs}  

// 这个函数设置了我们这里的DoFHandler对象（一个用于Darcy部分，一个用于饱和部分），以及将本程序中线性代数所需的各种对象设置为合适的尺寸。其基本操作与 step-31 所做的类似。

// 该函数的主体首先列举了达西和饱和系统的所有自由度。对于Darcy部分，自由度会被排序，以确保速度优先于压力DoF，这样我们就可以将Darcy矩阵划分为一个 $2 \times 2$ 矩阵。
//然后，
//我们需要将悬挂节点约束和Dirichlet边界值约束纳入 darcy_preconditioner_constraints。 边界条件约束只设置在压力分量上，因为对应于非混合形式的多孔介质流算子的Schur complement预处理程序 $-\nabla \cdot [\mathbf K \lambda_t(S)]\nabla$  ，只作用于压力变量。因此，我们使用一个过滤掉速度分量的分量掩码，这样就可以只对压力自由度进行缩减。

// 做完这些后，我们计算各个块中的自由度数量。然后，这些信息被用来创建达西和饱和系统矩阵的稀疏模式，以及用于建立达西预处理的预处理矩阵。如同 step-31 ，我们选择使用DynamicSparsityPattern的封锁版本来创建模式。因此，对于这一点，我们遵循与 step-31 相同的方式，对于成员函数的其他部分，我们不必再重复描述。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::setup_dofs() 
  { 
    std::vector<unsigned int> darcy_block_component(dim + 1, 0); 
    darcy_block_component[dim] = 1; 
    { 
      darcy_dof_handler.distribute_dofs(darcy_fe); 
      DoFRenumbering::Cuthill_McKee(darcy_dof_handler); 
      DoFRenumbering::component_wise(darcy_dof_handler, darcy_block_component); 

      darcy_constraints.clear(); 
      DoFTools::make_hanging_node_constraints(darcy_dof_handler, 
                                              darcy_constraints); 
      darcy_constraints.close(); 
    } 
    { 
      saturation_dof_handler.distribute_dofs(saturation_fe); 

      saturation_constraints.clear(); 
      DoFTools::make_hanging_node_constraints(saturation_dof_handler, 
                                              saturation_constraints); 
      saturation_constraints.close(); 
    } 
    { 
      darcy_preconditioner_constraints.clear(); 

      FEValuesExtractors::Scalar pressure(dim); 

      DoFTools::make_hanging_node_constraints(darcy_dof_handler, 
                                              darcy_preconditioner_constraints); 
      DoFTools::make_zero_boundary_constraints(darcy_dof_handler, 
                                               darcy_preconditioner_constraints, 
                                               darcy_fe.component_mask( 
                                                 pressure)); 

      darcy_preconditioner_constraints.close(); 
    } 

    const std::vector<types::global_dof_index> darcy_dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(darcy_dof_handler, 
                                        darcy_block_component); 
    const unsigned int n_u = darcy_dofs_per_block[0], 
                       n_p = darcy_dofs_per_block[1], 
                       n_s = saturation_dof_handler.n_dofs(); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << " (on " << triangulation.n_levels() << " levels)" << std::endl 
              << "Number of degrees of freedom: " << n_u + n_p + n_s << " (" 
              << n_u << '+' << n_p << '+' << n_s << ')' << std::endl 
              << std::endl; 

    { 
      darcy_matrix.clear(); 

      BlockDynamicSparsityPattern dsp(2, 2); 

      dsp.block(0, 0).reinit(n_u, n_u); 
      dsp.block(0, 1).reinit(n_u, n_p); 
      dsp.block(1, 0).reinit(n_p, n_u); 
      dsp.block(1, 1).reinit(n_p, n_p); 

      dsp.collect_sizes(); 

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 

      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (!((c == dim) && (d == dim))) 
            coupling[c][d] = DoFTools::always; 
          else 
            coupling[c][d] = DoFTools::none; 

      DoFTools::make_sparsity_pattern( 
        darcy_dof_handler, coupling, dsp, darcy_constraints, false); 

      darcy_matrix.reinit(dsp); 
    } 

    { 
      Amg_preconditioner.reset(); 
      Mp_preconditioner.reset(); 
      darcy_preconditioner_matrix.clear(); 

      BlockDynamicSparsityPattern dsp(2, 2); 

      dsp.block(0, 0).reinit(n_u, n_u); 
      dsp.block(0, 1).reinit(n_u, n_p); 
      dsp.block(1, 0).reinit(n_p, n_u); 
      dsp.block(1, 1).reinit(n_p, n_p); 

      dsp.collect_sizes(); 

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (c == d) 
            coupling[c][d] = DoFTools::always; 
          else 
            coupling[c][d] = DoFTools::none; 

      DoFTools::make_sparsity_pattern( 
        darcy_dof_handler, coupling, dsp, darcy_constraints, false); 

      darcy_preconditioner_matrix.reinit(dsp); 
    } 

    { 
      saturation_matrix.clear(); 

      DynamicSparsityPattern dsp(n_s, n_s); 

      DoFTools::make_sparsity_pattern(saturation_dof_handler, 
                                      dsp, 
                                      saturation_constraints, 
                                      false); 

      saturation_matrix.reinit(dsp); 
    } 

    std::vector<IndexSet> darcy_partitioning(2); 
    darcy_partitioning[0] = complete_index_set(n_u); 
    darcy_partitioning[1] = complete_index_set(n_p); 
    darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD); 
    darcy_solution.collect_sizes(); 

    last_computed_darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD); 
    last_computed_darcy_solution.collect_sizes(); 

    second_last_computed_darcy_solution.reinit(darcy_partitioning, 
                                               MPI_COMM_WORLD); 
    second_last_computed_darcy_solution.collect_sizes(); 

    darcy_rhs.reinit(darcy_partitioning, MPI_COMM_WORLD); 
    darcy_rhs.collect_sizes(); 

    IndexSet saturation_partitioning = complete_index_set(n_s); 
    saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 
    old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 
    old_old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 

    saturation_matching_last_computed_darcy_solution.reinit( 
      saturation_partitioning, MPI_COMM_WORLD); 

    saturation_rhs.reinit(saturation_partitioning, MPI_COMM_WORLD); 
  } 
// @sect3{Assembling matrices and preconditioners}  

// 接下来的几个函数专门用来设置我们在这个程序中必须处理的各种系统和预处理矩阵及右手边。

//  @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner}  

// 这个函数组装我们用于预处理达西系统的矩阵。我们需要的是在速度分量上用 $\left(\mathbf{K} \lambda_t\right)^{-1}$ 加权的向量质量矩阵和在压力分量上用 $\left(\mathbf{K} \lambda_t\right)$ 加权的质量矩阵。我们首先生成一个适当阶数的正交对象，即FEValues对象，可以给出正交点的数值和梯度（连同正交权重）。接下来我们为单元格矩阵和局部与全局DoF之间的关系创建数据结构。向量phi_u和grad_phi_p将保存基函数的值，以便更快地建立局部矩阵，正如在  step-22  中已经做的。在我们开始对所有活动单元进行循环之前，我们必须指定哪些成分是压力，哪些是速度。

// 局部矩阵的创建是相当简单的。只有一个由 $\left(\mathbf{K} \lambda_t\right)^{-1}$ 加权的项（关于速度）和一个由 $\left(\mathbf{K} \lambda_t\right)$ 加权的拉普拉斯矩阵需要生成，所以局部矩阵的创建基本上只需要两行就可以完成。由于该文件顶部的材料模型函数只提供了渗透率和迁移率的倒数，我们必须根据给定的数值手工计算 $\mathbf K$ 和 $\lambda_t$ ，每个正交点一次。

// 一旦本地矩阵准备好了（在每个正交点上对本地矩阵的行和列进行循环），我们就可以得到本地的DoF指数，并将本地信息写入全局矩阵中。我们通过直接应用约束条件（即darcy_preconditioner_constraints）来做到这一点，该约束条件负责处理悬挂节点和零Dirichlet边界条件约束。这样做，我们就不必事后再做，以后也不必使用 AffineConstraints::condense 和 MatrixTools::apply_boundary_values, 这两个需要修改矩阵和向量项的函数，因此对于我们不能立即访问单个内存位置的特里诺斯类来说，很难编写。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner() 
  { 
    std::cout << "   Rebuilding darcy preconditioner..." << std::endl; 

    darcy_preconditioner_matrix = 0; 

    const QGauss<dim> quadrature_formula(darcy_degree + 2); 
    FEValues<dim>     darcy_fe_values(darcy_fe, 
                                  quadrature_formula, 
                                  update_JxW_values | update_values | 
                                    update_gradients | 
                                    update_quadrature_points); 
    FEValues<dim>     saturation_fe_values(saturation_fe, 
                                       quadrature_formula, 
                                       update_values); 

    const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

    std::vector<double> old_saturation_values(n_q_points); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
    std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

    auto       cell            = darcy_dof_handler.begin_active(); 
    const auto endc            = darcy_dof_handler.end(); 
    auto       saturation_cell = saturation_dof_handler.begin_active(); 

    for (; cell != endc; ++cell, ++saturation_cell) 
      { 
        darcy_fe_values.reinit(cell); 
        saturation_fe_values.reinit(saturation_cell); 

        local_matrix = 0; 

        saturation_fe_values.get_function_values(old_saturation_solution, 
                                                 old_saturation_values); 

        k_inverse.value_list(darcy_fe_values.get_quadrature_points(), 
                             k_inverse_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double old_s = old_saturation_values[q]; 

            const double inverse_mobility = mobility_inverse(old_s, viscosity); 
            const double mobility         = 1.0 / inverse_mobility; 
            const Tensor<2, dim> permeability = invert(k_inverse_values[q]); 

            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                phi_u[k]      = darcy_fe_values[velocities].value(k, q); 
                grad_phi_p[k] = darcy_fe_values[pressure].gradient(k, q); 
              } 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  local_matrix(i, j) += 
                    (k_inverse_values[q] * inverse_mobility * phi_u[i] * 
                       phi_u[j] + 
                     permeability * mobility * grad_phi_p[i] * grad_phi_p[j]) * 
                    darcy_fe_values.JxW(q); 
                } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        darcy_preconditioner_constraints.distribute_local_to_global( 
          local_matrix, local_dof_indices, darcy_preconditioner_matrix); 
      } 
  } 
// @sect4{TwoPhaseFlowProblem<dim>::build_darcy_preconditioner}  

// 在调用上述函数组装预处理矩阵后，该函数生成将用于舒尔补块预处理的内部预处理器。前置条件需要在每个饱和时间步长时重新生成，因为它们取决于随时间变化的饱和度  $S$  。

// 在这里，我们为速度-速度矩阵  $\mathbf{M}^{\mathbf{u}}$  和Schur补码  $\mathbf{S}$  设置了预处理器。正如介绍中所解释的，我们将使用一个基于矢量矩阵 $\mathbf{M}^{\mathbf{u}}$ 的IC预处理器和另一个基于标量拉普拉斯矩阵 $\tilde{\mathbf{S}}^p$ 的IC预处理器（它在频谱上与达西矩阵的舒尔补码接近）。通常， TrilinosWrappers::PreconditionIC 类可以被看作是一个很好的黑盒预处理程序，不需要对矩阵结构和/或背后的算子有任何特殊的了解。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::build_darcy_preconditioner() 
  { 
    assemble_darcy_preconditioner(); 

    Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
    Amg_preconditioner->initialize(darcy_preconditioner_matrix.block(0, 0)); 

    Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
    Mp_preconditioner->initialize(darcy_preconditioner_matrix.block(1, 1)); 
  } 
// @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_system}  

// 这是为达西系统组装线性系统的函数。

// 关于执行的技术细节，其程序与  step-22  和  step-31  中的程序相似。我们重置矩阵和向量，在单元格上创建正交公式，然后创建相应的FEValues对象。

// 有一件事需要评论：由于我们有一个单独的有限元和DoFHandler来处理饱和问题，我们需要生成第二个FEValues对象来正确评估饱和解。要实现这一点并不复杂：只需使用饱和结构，并为基函数值设置一个更新标志，我们需要对饱和解进行评估。这里需要记住的唯一重要部分是，两个FEValues对象使用相同的正交公式，以确保我们在循环计算两个对象的正交点时获得匹配的信息。

// 声明的过程中，对数组的大小、本地矩阵的创建、右手边以及与全局系统相比较的本地道夫指数的向量都有一些快捷方式。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_darcy_system() 
  { 
    darcy_matrix = 0; 
    darcy_rhs    = 0; 

    QGauss<dim>     quadrature_formula(darcy_degree + 2); 
    QGauss<dim - 1> face_quadrature_formula(darcy_degree + 2); 

    FEValues<dim> darcy_fe_values(darcy_fe, 
                                  quadrature_formula, 
                                  update_values | update_gradients | 
                                    update_quadrature_points | 
                                    update_JxW_values); 

    FEValues<dim> saturation_fe_values(saturation_fe, 
                                       quadrature_formula, 
                                       update_values); 

    FEFaceValues<dim> darcy_fe_face_values(darcy_fe, 
                                           face_quadrature_formula, 
                                           update_values | 
                                             update_normal_vectors | 
                                             update_quadrature_points | 
                                             update_JxW_values); 

    const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell(); 

    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const Functions::ZeroFunction<dim> pressure_right_hand_side; 
    const PressureBoundaryValues<dim>  pressure_boundary_values; 

    std::vector<double>         pressure_rhs_values(n_q_points); 
    std::vector<double>         boundary_values(n_face_q_points); 
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

// 接下来我们需要一个向量，该向量将包含前一时间层在正交点的饱和解的值，以组装达西方程中的饱和相关系数。

// 我们接下来创建的向量集包含了基函数的评价以及它们的梯度，将用于创建矩阵。把这些放到自己的数组中，而不是每次都向FEValues对象索取这些信息，是为了加速装配过程的优化，详情请见 step-22 。

// 最后两个声明是用来从整个FE系统中提取各个块（速度、压力、饱和度）的。

    std::vector<double> old_saturation_values(n_q_points); 

    std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
    std::vector<double>         div_phi_u(dofs_per_cell); 
    std::vector<double>         phi_p(dofs_per_cell); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

// 现在开始对问题中的所有单元格进行循环。我们在这个装配例程中使用了两个不同的DoFHandlers，所以我们必须为使用中的两个对象设置两个不同的单元格迭代器。这可能看起来有点奇怪，但是由于达西系统和饱和系统都使用相同的网格，我们可以假设这两个迭代器在两个DoFHandler对象的单元格中同步运行。

// 循环中的第一条语句又是非常熟悉的，按照更新标志的规定对有限元数据进行更新，将局部数组清零，并得到正交点上的旧解的值。 在这一点上，我们还必须在正交点上获得前一个时间步长的饱和函数的值。为此，我们可以使用 FEValues::get_function_values （之前已经在 step-9 、 step-14 和 step-15 中使用），这个函数接收一个解向量，并返回当前单元的正交点的函数值列表。事实上，它返回每个正交点的完整矢量值解，即不仅是饱和度，还有速度和压力。

// 然后，我们就可以在单元格上的正交点上进行循环，以进行积分。这方面的公式直接来自介绍中所讨论的内容。

// 一旦这样做了，我们就开始在局部矩阵的行和列上进行循环，并将相关的乘积输入矩阵中。

// 循环所有单元的最后一步是将本地贡献输入到全局矩阵和向量结构中，并在local_dof_indices中指定位置。同样，我们让AffineConstraints类将单元格矩阵元素插入到全局矩阵中，全局矩阵已经浓缩了悬挂节点的约束。

    auto       cell            = darcy_dof_handler.begin_active(); 
    const auto endc            = darcy_dof_handler.end(); 
    auto       saturation_cell = saturation_dof_handler.begin_active(); 

    for (; cell != endc; ++cell, ++saturation_cell) 
      { 
        darcy_fe_values.reinit(cell); 
        saturation_fe_values.reinit(saturation_cell); 

        local_matrix = 0; 
        local_rhs    = 0; 

        saturation_fe_values.get_function_values(old_saturation_solution, 
                                                 old_saturation_values); 

        pressure_right_hand_side.value_list( 
          darcy_fe_values.get_quadrature_points(), pressure_rhs_values); 
        k_inverse.value_list(darcy_fe_values.get_quadrature_points(), 
                             k_inverse_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                phi_u[k]     = darcy_fe_values[velocities].value(k, q); 
                div_phi_u[k] = darcy_fe_values[velocities].divergence(k, q); 
                phi_p[k]     = darcy_fe_values[pressure].value(k, q); 
              } 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                const double old_s = old_saturation_values[q]; 
                for (unsigned int j = 0; j <= i; ++j) 
                  { 
                    local_matrix(i, j) += 
                      (phi_u[i] * k_inverse_values[q] * 
                         mobility_inverse(old_s, viscosity) * phi_u[j] - 
                       div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
                      darcy_fe_values.JxW(q); 
                  } 

                local_rhs(i) += 
                  (-phi_p[i] * pressure_rhs_values[q]) * darcy_fe_values.JxW(q); 
              } 
          } 

        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              darcy_fe_face_values.reinit(cell, face); 

              pressure_boundary_values.value_list( 
                darcy_fe_face_values.get_quadrature_points(), boundary_values); 

              for (unsigned int q = 0; q < n_face_q_points; ++q) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  { 
                    const Tensor<1, dim> phi_i_u = 
                      darcy_fe_face_values[velocities].value(i, q); 

                    local_rhs(i) += 
                      -(phi_i_u * darcy_fe_face_values.normal_vector(q) * 
                        boundary_values[q] * darcy_fe_face_values.JxW(q)); 
                  } 
            } 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
            local_matrix(i, j) = local_matrix(j, i); 

        cell->get_dof_indices(local_dof_indices); 

        darcy_constraints.distribute_local_to_global( 
          local_matrix, local_rhs, local_dof_indices, darcy_matrix, darcy_rhs); 
      } 
  } 
// @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_system}  

// 这个函数是为了组装饱和传输方程的线性系统。如果有必要，它会调用另外两个成员函数：assemble_saturation_matrix()和assemble_saturation_rhs()。前一个函数然后组装饱和度矩阵，只需要偶尔改变。另一方面，后一个组装右手边的函数必须在每个饱和时间步骤中调用。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_saturation_system() 
  { 
    if (rebuild_saturation_matrix == true) 
      { 
        saturation_matrix = 0; 
        assemble_saturation_matrix(); 
      } 

    saturation_rhs = 0; 
    assemble_saturation_rhs(); 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_matrix}  

// 这个函数很容易理解，因为它只是通过基函数phi_i_s和phi_j_s为饱和线性系统的左侧形成一个简单的质量矩阵。最后，像往常一样，我们通过在local_dof_indices中指定位置将局部贡献输入全局矩阵。这是通过让AffineConstraints类将单元矩阵元素插入全局矩阵来完成的，全局矩阵已经浓缩了悬挂节点约束。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix() 
  { 
    QGauss<dim> quadrature_formula(saturation_degree + 2); 

    FEValues<dim> saturation_fe_values(saturation_fe, 
                                       quadrature_formula, 
                                       update_values | update_JxW_values); 

    const unsigned int dofs_per_cell = saturation_fe.n_dofs_per_cell(); 

    const unsigned int n_q_points = quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
      { 
        saturation_fe_values.reinit(cell); 
        local_matrix = 0; 
        local_rhs    = 0; 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const double phi_i_s = saturation_fe_values.shape_value(i, q); 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const double phi_j_s = saturation_fe_values.shape_value(j, q); 
                  local_matrix(i, j) += 
                    porosity * phi_i_s * phi_j_s * saturation_fe_values.JxW(q); 
                } 
            } 
        cell->get_dof_indices(local_dof_indices); 

        saturation_constraints.distribute_local_to_global(local_matrix, 
                                                          local_dof_indices, 
                                                          saturation_matrix); 
      } 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs}  

// 这个函数是用来组装饱和传输方程的右边。在进行这项工作之前，我们必须为达西系统和饱和系统分别创建两个FEValues对象，此外，还必须为这两个系统创建两个FEFaceValues对象，因为我们在饱和方程的弱形式中存在一个边界积分项。对于饱和系统的FEFaceValues对象，我们还需要法向量，我们使用update_normal_vectors标志来申请。

// 接下来，在对所有单元进行循环之前，我们必须计算一些参数（例如global_u_infty、global_S_variation和global_Omega_diameter），这是人工黏度 $\nu$ 需要的。这与 step-31 中的做法基本相同，所以你可以在那里看到更多的信息。

// 真正的工作是从循环所有的饱和和Darcy单元开始的，以便将局部贡献放到全局矢量中。在这个循环中，为了简化实现，我们把一些工作分成两个辅助函数：assemble_saturation_rhs_cell_term和assemble_saturation_rhs_boundary_term。 我们注意到，我们在这两个函数中把细胞或边界贡献插入全局向量，而不是在本函数中。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs() 
  { 
    QGauss<dim>     quadrature_formula(saturation_degree + 2); 
    QGauss<dim - 1> face_quadrature_formula(saturation_degree + 2); 

    FEValues<dim> saturation_fe_values(saturation_fe, 
                                       quadrature_formula, 
                                       update_values | update_gradients | 
                                         update_quadrature_points | 
                                         update_JxW_values); 
    FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values); 
    FEFaceValues<dim> saturation_fe_face_values(saturation_fe, 
                                                face_quadrature_formula, 
                                                update_values | 
                                                  update_normal_vectors | 
                                                  update_quadrature_points | 
                                                  update_JxW_values); 
    FEFaceValues<dim> darcy_fe_face_values(darcy_fe, 
                                           face_quadrature_formula, 
                                           update_values); 
    FEFaceValues<dim> saturation_fe_face_values_neighbor( 
      saturation_fe, face_quadrature_formula, update_values); 

    const unsigned int dofs_per_cell = 
      saturation_dof_handler.get_fe().n_dofs_per_cell(); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const double                    global_max_u_F_prime = get_max_u_F_prime(); 
    const std::pair<double, double> global_S_range = 
      get_extrapolated_saturation_range(); 
    const double global_S_variation = 
      global_S_range.second - global_S_range.first; 

    auto       cell       = saturation_dof_handler.begin_active(); 
    const auto endc       = saturation_dof_handler.end(); 
    auto       darcy_cell = darcy_dof_handler.begin_active(); 
    for (; cell != endc; ++cell, ++darcy_cell) 
      { 
        saturation_fe_values.reinit(cell); 
        darcy_fe_values.reinit(darcy_cell); 

        cell->get_dof_indices(local_dof_indices); 

        assemble_saturation_rhs_cell_term(saturation_fe_values, 
                                          darcy_fe_values, 
                                          global_max_u_F_prime, 
                                          global_S_variation, 
                                          local_dof_indices); 

        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              darcy_fe_face_values.reinit(darcy_cell, face); 
              saturation_fe_face_values.reinit(cell, face); 
              assemble_saturation_rhs_boundary_term(saturation_fe_face_values, 
                                                    darcy_fe_face_values, 
                                                    local_dof_indices); 
            } 
      } 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term}  

// 这个函数负责整合饱和度方程右边的单元项，然后将其组装成全局右边的矢量。鉴于介绍中的讨论，这些贡献的形式很清楚。唯一棘手的部分是获得人工黏度和计算它所需的一切。该函数的前半部分专门用于这项任务。

// 该函数的最后一部分是将局部贡献复制到全局向量中，其位置由local_dof_indices指定。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term( 
    const FEValues<dim> &                       saturation_fe_values, 
    const FEValues<dim> &                       darcy_fe_values, 
    const double                                global_max_u_F_prime, 
    const double                                global_S_variation, 
    const std::vector<types::global_dof_index> &local_dof_indices) 
  { 
    const unsigned int dofs_per_cell = saturation_fe_values.dofs_per_cell; 
    const unsigned int n_q_points    = saturation_fe_values.n_quadrature_points; 

    std::vector<double>         old_saturation_solution_values(n_q_points); 
    std::vector<double>         old_old_saturation_solution_values(n_q_points); 
    std::vector<Tensor<1, dim>> old_grad_saturation_solution_values(n_q_points); 
    std::vector<Tensor<1, dim>> old_old_grad_saturation_solution_values( 
      n_q_points); 
    std::vector<Vector<double>> present_darcy_solution_values( 
      n_q_points, Vector<double>(dim + 1)); 

    saturation_fe_values.get_function_values(old_saturation_solution, 
                                             old_saturation_solution_values); 
    saturation_fe_values.get_function_values( 
      old_old_saturation_solution, old_old_saturation_solution_values); 
    saturation_fe_values.get_function_gradients( 
      old_saturation_solution, old_grad_saturation_solution_values); 
    saturation_fe_values.get_function_gradients( 
      old_old_saturation_solution, old_old_grad_saturation_solution_values); 
    darcy_fe_values.get_function_values(darcy_solution, 
                                        present_darcy_solution_values); 

    const double nu = 
      compute_viscosity(old_saturation_solution_values, 
                        old_old_saturation_solution_values, 
                        old_grad_saturation_solution_values, 
                        old_old_grad_saturation_solution_values, 
                        present_darcy_solution_values, 
                        global_max_u_F_prime, 
                        global_S_variation, 
                        saturation_fe_values.get_cell()->diameter()); 

    Vector<double> local_rhs(dofs_per_cell); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          const double   old_s = old_saturation_solution_values[q]; 
          Tensor<1, dim> present_u; 
          for (unsigned int d = 0; d < dim; ++d) 
            present_u[d] = present_darcy_solution_values[q](d); 

          const double         phi_i_s = saturation_fe_values.shape_value(i, q); 
          const Tensor<1, dim> grad_phi_i_s = 
            saturation_fe_values.shape_grad(i, q); 

          local_rhs(i) += 
            (time_step * fractional_flow(old_s, viscosity) * present_u * 
               grad_phi_i_s - 
             time_step * nu * old_grad_saturation_solution_values[q] * 
               grad_phi_i_s + 
             porosity * old_s * phi_i_s) * 
            saturation_fe_values.JxW(q); 
        } 

    saturation_constraints.distribute_local_to_global(local_rhs, 
                                                      local_dof_indices, 
                                                      saturation_rhs); 
  } 
// @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term}  

// 下一个函数负责饱和方程右侧形式中的边界积分项。 对于这些，我们必须计算全局边界面上的上行通量，也就是说，我们只对全局边界的流入部分弱加迪里切特边界条件。如前所述，这在 step-21 中已经描述过了，所以我们不对其进行更多的描述。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term( 
    const FEFaceValues<dim> &                   saturation_fe_face_values, 
    const FEFaceValues<dim> &                   darcy_fe_face_values, 
    const std::vector<types::global_dof_index> &local_dof_indices) 
  { 
    const unsigned int dofs_per_cell = saturation_fe_face_values.dofs_per_cell; 
    const unsigned int n_face_q_points = 
      saturation_fe_face_values.n_quadrature_points; 

    Vector<double> local_rhs(dofs_per_cell); 

 
    std::vector<Vector<double>> present_darcy_solution_values_face( 
      n_face_q_points, Vector<double>(dim + 1)); 
    std::vector<double> neighbor_saturation(n_face_q_points); 

    saturation_fe_face_values.get_function_values( 
      old_saturation_solution, old_saturation_solution_values_face); 
    darcy_fe_face_values.get_function_values( 
      darcy_solution, present_darcy_solution_values_face); 

    SaturationBoundaryValues<dim> saturation_boundary_values; 
    saturation_boundary_values.value_list( 
      saturation_fe_face_values.get_quadrature_points(), neighbor_saturation); 

    for (unsigned int q = 0; q < n_face_q_points; ++q) 
      { 
        Tensor<1, dim> present_u_face; 
        for (unsigned int d = 0; d < dim; ++d) 
          present_u_face[d] = present_darcy_solution_values_face[q](d); 

 
          present_u_face * saturation_fe_face_values.normal_vector(q); 

        const bool is_outflow_q_point = (normal_flux >= 0); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          local_rhs(i) -= 
            time_step * normal_flux * 
            fractional_flow((is_outflow_q_point == true ? 
                               old_saturation_solution_values_face[q] : 
                               neighbor_saturation[q]), 
                            viscosity) * 
            saturation_fe_face_values.shape_value(i, q) * 
            saturation_fe_face_values.JxW(q); 
      } 
    saturation_constraints.distribute_local_to_global(local_rhs, 
                                                      local_dof_indices, 
                                                      saturation_rhs); 
  } 
// @sect3{TwoPhaseFlowProblem<dim>::solve}  

// 该函数实现了算子分割算法，即在每个时间步长中，它要么重新计算达西系统的解，要么从以前的时间步长中推算出速度/压力，然后确定时间步长的大小，然后更新饱和度变量。其实现主要遵循  step-31  中的类似代码。除了run()函数外，它是本程序中的核心函数。

// 在函数的开始，我们询问是否要通过评估后验准则来解决压力-速度部分（见下面的函数）。如果有必要，我们将使用GMRES求解器和Schur补充块预处理来求解压力-速度部分，如介绍中所述。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::solve() 
  { 
    const bool solve_for_pressure_and_velocity = 
      determine_whether_to_solve_for_pressure_and_velocity(); 

    if (solve_for_pressure_and_velocity == true) 
      { 
        std::cout << "   Solving Darcy (pressure-velocity) system..." 
                  << std::endl; 

        assemble_darcy_system(); 
        build_darcy_preconditioner(); 

        { 
          const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix, 
                                             TrilinosWrappers::PreconditionIC> 
            mp_inverse(darcy_preconditioner_matrix.block(1, 1), 
                       *Mp_preconditioner); 

          const LinearSolvers::BlockSchurPreconditioner< 
            TrilinosWrappers::PreconditionIC, 
            TrilinosWrappers::PreconditionIC> 
            preconditioner(darcy_matrix, mp_inverse, *Amg_preconditioner); 

          SolverControl solver_control(darcy_matrix.m(), 
                                       1e-16 * darcy_rhs.l2_norm()); 

          SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres( 
            solver_control, 
            SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData( 
              100)); 

          for (unsigned int i = 0; i < darcy_solution.size(); ++i) 
            if (darcy_constraints.is_constrained(i)) 
              darcy_solution(i) = 0; 

          gmres.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner); 

          darcy_constraints.distribute(darcy_solution); 

          std::cout << "        ..." << solver_control.last_step() 
                    << " GMRES iterations." << std::endl; 
        } 

        { 
  }; 
          last_computed_darcy_solution        = darcy_solution; 

          saturation_matching_last_computed_darcy_solution = 
            saturation_solution; 
        } 
      } 

// 另一方面，如果我们决定不计算当前时间步长的达西系统的解，那么我们需要简单地将前两个达西解外推到与我们计算速度/压力的时间相同。我们做一个简单的线性外推，即给定从上次计算达西解到现在的宏观时间步长 $dt$ （由 <code>current_macro_time_step</code> 给出），以及 $DT$ 上一个宏观时间步长（由 <code>old_macro_time_step</code> 给出），然后得到 $u^\ast = u_p + dt \frac{u_p-u_{pp}}{DT} = (1+dt/DT)u_p - dt/DT u_{pp}$  ，其中 $u_p$ 和 $u_{pp}$ 是最近两个计算的达西解。我们只需用两行代码就可以实现这个公式。

// 请注意，这里的算法只有在我们至少有两个先前计算的Darcy解，我们可以从中推断出当前的时间，这一点通过要求重新计算前两个时间步骤的Darcy解来保证。

    else 
      { 
        darcy_solution = last_computed_darcy_solution; 
        darcy_solution.sadd(1 + current_macro_time_step / old_macro_time_step, 
                            -current_macro_time_step / old_macro_time_step, 
                            second_last_computed_darcy_solution); 
      } 

// 用这样计算出来的速度矢量，根据介绍中讨论的CFL标准计算出最佳时间步长......

    { 
      old_time_step = time_step; 

      const double max_u_F_prime = get_max_u_F_prime(); 
      if (max_u_F_prime > 0) 
        time_step = porosity * GridTools::minimal_cell_diameter(triangulation) / 
                    saturation_degree / max_u_F_prime / 50; 
      else 
        time_step = end_time - time; 
    } 

// ......然后在我们处理时间步长的时候，还要更新我们使用的宏观时间步长。具体而言，这涉及到。(i) 如果我们刚刚重新计算了达西解，那么之前的宏观时间步长现在是固定的，当前的宏观时间步长，到现在为止，只是当前（微观）时间步长。(ii) 如果我们没有重新计算达西解，那么当前的宏观时间步长刚刚增长了 <code>time_step</code>  。

    if (solve_for_pressure_and_velocity == true) 
      { 
        old_macro_time_step     = current_macro_time_step; 
        current_macro_time_step = time_step; 
      } 
    else 
      current_macro_time_step += time_step; 

// 这个函数的最后一步是根据我们刚刚得到的速度场重新计算饱和解。这自然发生在每一个时间步骤中，我们不会跳过这些计算。在计算饱和度的最后，我们投射回允许的区间 $[0,1]$ ，以确保我们的解保持物理状态。

    { 
      std::cout << "   Solving saturation transport equation..." << std::endl; 

      assemble_saturation_system(); 

      SolverControl solver_control(saturation_matrix.m(), 
                                   1e-16 * saturation_rhs.l2_norm()); 
      SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 

      TrilinosWrappers::PreconditionIC preconditioner; 
      preconditioner.initialize(saturation_matrix); 

      cg.solve(saturation_matrix, 
               saturation_solution, 
               saturation_rhs, 
               preconditioner); 

      saturation_constraints.distribute(saturation_solution); 
      project_back_saturation(); 

      std::cout << "        ..." << solver_control.last_step() 
                << " CG iterations." << std::endl; 
    } 
  } 
// @sect3{TwoPhaseFlowProblem<dim>::refine_mesh}  

// 下一个函数是对网格进行细化和粗化。它的工作分三块进行。(i) 计算细化指标，方法是通过使用各自的时间步长（如果这是第一个时间步长，则取唯一的解决方案），从前两个时间步长中线性推断出的解决方案向量的梯度。(ii) 在梯度大于或小于某一阈值的单元中标记出细化和粗化的单元，保留网格细化的最小和最大水平。(iii) 将解决方案从旧网格转移到新网格。这些都不是特别困难。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::refine_mesh(const unsigned int min_grid_level, 
                                             const unsigned int max_grid_level) 
  { 
    Vector<double> refinement_indicators(triangulation.n_active_cells()); 
    { 
      const QMidpoint<dim>        quadrature_formula; 
      FEValues<dim>               fe_values(saturation_fe, 
                              quadrature_formula, 
                              update_gradients); 
      std::vector<Tensor<1, dim>> grad_saturation(1); 

      TrilinosWrappers::MPI::Vector extrapolated_saturation_solution( 
        saturation_solution); 
      if (timestep_number != 0) 
        extrapolated_saturation_solution.sadd((1. + time_step / old_time_step), 
                                              time_step / old_time_step, 
                                              old_saturation_solution); 

      for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
        { 
          const unsigned int cell_no = cell->active_cell_index(); 
          fe_values.reinit(cell); 
          fe_values.get_function_gradients(extrapolated_saturation_solution, 
                                           grad_saturation); 

          refinement_indicators(cell_no) = grad_saturation[0].norm(); 
        } 
    } 

    { 
      for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
        { 
          const unsigned int cell_no = cell->active_cell_index(); 
          cell->clear_coarsen_flag(); 
          cell->clear_refine_flag(); 

          if ((static_cast<unsigned int>(cell->level()) < max_grid_level) && 
              (std::fabs(refinement_indicators(cell_no)) > 
               saturation_refinement_threshold)) 
            cell->set_refine_flag(); 
          else if ((static_cast<unsigned int>(cell->level()) > 
                    min_grid_level) && 
                   (std::fabs(refinement_indicators(cell_no)) < 
                    0.5 * saturation_refinement_threshold)) 
            cell->set_coarsen_flag(); 
        } 
    } 

    triangulation.prepare_coarsening_and_refinement(); 

    { 
      std::vector<TrilinosWrappers::MPI::Vector> x_saturation(3); 
      x_saturation[0] = saturation_solution; 
      x_saturation[1] = old_saturation_solution; 
      x_saturation[2] = saturation_matching_last_computed_darcy_solution; 

      std::vector<TrilinosWrappers::MPI::BlockVector> x_darcy(2); 
      x_darcy[0] = last_computed_darcy_solution; 
      x_darcy[1] = second_last_computed_darcy_solution; 

      SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> saturation_soltrans( 
        saturation_dof_handler); 

      SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> darcy_soltrans( 
        darcy_dof_handler); 

      triangulation.prepare_coarsening_and_refinement(); 
      saturation_soltrans.prepare_for_coarsening_and_refinement(x_saturation); 

      darcy_soltrans.prepare_for_coarsening_and_refinement(x_darcy); 

      triangulation.execute_coarsening_and_refinement(); 
      setup_dofs(); 

      std::vector<TrilinosWrappers::MPI::Vector> tmp_saturation(3); 
      tmp_saturation[0].reinit(saturation_solution); 
      tmp_saturation[1].reinit(saturation_solution); 
      tmp_saturation[2].reinit(saturation_solution); 
      saturation_soltrans.interpolate(x_saturation, tmp_saturation); 

      saturation_solution                              = tmp_saturation[0]; 
      old_saturation_solution                          = tmp_saturation[1]; 
      saturation_matching_last_computed_darcy_solution = tmp_saturation[2]; 

      saturation_constraints.distribute(saturation_solution); 
      saturation_constraints.distribute(old_saturation_solution); 
      saturation_constraints.distribute( 
        saturation_matching_last_computed_darcy_solution); 

      std::vector<TrilinosWrappers::MPI::BlockVector> tmp_darcy(2); 
      tmp_darcy[0].reinit(darcy_solution); 
      tmp_darcy[1].reinit(darcy_solution); 
      darcy_soltrans.interpolate(x_darcy, tmp_darcy); 

      last_computed_darcy_solution        = tmp_darcy[0]; 
      second_last_computed_darcy_solution = tmp_darcy[1]; 

      darcy_constraints.distribute(last_computed_darcy_solution); 
      darcy_constraints.distribute(second_last_computed_darcy_solution); 

      rebuild_saturation_matrix = true; 
    } 
  } 

//  @sect3{TwoPhaseFlowProblem<dim>::output_results}  

// 这个函数生成图形输出。它实质上是对  step-31  中实现的复制。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::output_results() const 
  { 
    const FESystem<dim> joint_fe(darcy_fe, 1, saturation_fe, 1); 
    DoFHandler<dim>     joint_dof_handler(triangulation); 
    joint_dof_handler.distribute_dofs(joint_fe); 
    Assert(joint_dof_handler.n_dofs() == 
             darcy_dof_handler.n_dofs() + saturation_dof_handler.n_dofs(), 
           ExcInternalError()); 

    Vector<double> joint_solution(joint_dof_handler.n_dofs()); 

    { 
      std::vector<types::global_dof_index> local_joint_dof_indices( 
        joint_fe.n_dofs_per_cell()); 
      std::vector<types::global_dof_index> local_darcy_dof_indices( 
        darcy_fe.n_dofs_per_cell()); 
      std::vector<types::global_dof_index> local_saturation_dof_indices( 
        saturation_fe.n_dofs_per_cell()); 

      auto       joint_cell      = joint_dof_handler.begin_active(); 
      const auto joint_endc      = joint_dof_handler.end(); 
      auto       darcy_cell      = darcy_dof_handler.begin_active(); 
      auto       saturation_cell = saturation_dof_handler.begin_active(); 

      for (; joint_cell != joint_endc; 
           ++joint_cell, ++darcy_cell, ++saturation_cell) 
        { 
          joint_cell->get_dof_indices(local_joint_dof_indices); 
          darcy_cell->get_dof_indices(local_darcy_dof_indices); 
          saturation_cell->get_dof_indices(local_saturation_dof_indices); 

          for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i) 
            if (joint_fe.system_to_base_index(i).first.first == 0) 
              { 
                Assert(joint_fe.system_to_base_index(i).second < 
                         local_darcy_dof_indices.size(), 
                       ExcInternalError()); 
                joint_solution(local_joint_dof_indices[i]) = darcy_solution( 
                  local_darcy_dof_indices[joint_fe.system_to_base_index(i) 
                                            .second]); 
              } 
            else 
              { 
                Assert(joint_fe.system_to_base_index(i).first.first == 1, 
                       ExcInternalError()); 
                Assert(joint_fe.system_to_base_index(i).second < 
                         local_darcy_dof_indices.size(), 
                       ExcInternalError()); 
                joint_solution(local_joint_dof_indices[i]) = 
                  saturation_solution( 
                    local_saturation_dof_indices 
                      [joint_fe.system_to_base_index(i).second]); 
              } 
        } 
    } 
    std::vector<std::string> joint_solution_names(dim, "velocity"); 
    joint_solution_names.emplace_back("pressure"); 
    joint_solution_names.emplace_back("saturation"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(joint_dof_handler); 
    data_out.add_data_vector(joint_solution, 
                             joint_solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 

    data_out.build_patches(); 

    std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 5) + ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect3{Tool functions}  
// @sect4{TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity}  

// 这个函数实现了自适应运算符拆分的后验标准。考虑到我们在上面实现其他函数的方式，并考虑到论文中得出的准则公式，该函数是相对简单的。

// 如果我们决定要采用原始的IMPES方法，即在每个时间步长中求解Darcy方程，那么可以通过将阈值 <code>AOS_threshold</code> （默认为 $5.0$ ）设置为0来实现，从而迫使该函数总是返回true。

// 最后，请注意，该函数在前两个时间步骤中无条件地返回真，以确保我们在跳过达西系统的解时总是至少解了两次，从而允许我们从 <code>solve()</code> 中的最后两次解中推算出速度。

  template <int dim> 
  bool TwoPhaseFlowProblem< 
    dim>::determine_whether_to_solve_for_pressure_and_velocity() const 
  { 
    if (timestep_number <= 2) 
      return true; 

    const QGauss<dim>  quadrature_formula(saturation_degree + 2); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(saturation_fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points); 

    std::vector<double> old_saturation_after_solving_pressure(n_q_points); 
    std::vector<double> present_saturation(n_q_points); 

    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

    double max_global_aop_indicator = 0.0; 

    for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
      { 
        double max_local_mobility_reciprocal_difference = 0.0; 
        double max_local_permeability_inverse_l1_norm   = 0.0; 

        fe_values.reinit(cell); 
        fe_values.get_function_values( 
          saturation_matching_last_computed_darcy_solution, 
          old_saturation_after_solving_pressure); 
        fe_values.get_function_values(saturation_solution, present_saturation); 

        k_inverse.value_list(fe_values.get_quadrature_points(), 
                             k_inverse_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double mobility_reciprocal_difference = std::fabs( 
              mobility_inverse(present_saturation[q], viscosity) - 
              mobility_inverse(old_saturation_after_solving_pressure[q], 
                               viscosity)); 

            max_local_mobility_reciprocal_difference = 
              std::max(max_local_mobility_reciprocal_difference, 
                       mobility_reciprocal_difference); 

            max_local_permeability_inverse_l1_norm = 
              std::max(max_local_permeability_inverse_l1_norm, 
                       l1_norm(k_inverse_values[q])); 
          } 

        max_global_aop_indicator = 
          std::max(max_global_aop_indicator, 
                   (max_local_mobility_reciprocal_difference * 
                    max_local_permeability_inverse_l1_norm)); 
      } 

    return (max_global_aop_indicator > AOS_threshold); 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::project_back_saturation}  

// 下一个函数只是确保饱和度值始终保持在  $[0,1]$  的物理合理范围内。虽然连续方程保证了这一点，但离散方程并没有。然而，如果我们允许离散解逃脱这个范围，我们就会遇到麻烦，因为像 $F(S)$ 和 $F'(S)$ 这样的项会产生不合理的结果（例如 $F'(S)<0$ 为 $S<0$ ，这将意味着润湿液相的流动方向为<i>against</i>的散流体速度））。因此，在每个时间步骤结束时，我们只需将饱和场投射回物理上合理的区域。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::project_back_saturation() 
  { 
    for (unsigned int i = 0; i < saturation_solution.size(); ++i) 
      if (saturation_solution(i) < 0.2) 
        saturation_solution(i) = 0.2; 
      else if (saturation_solution(i) > 1) 
        saturation_solution(i) = 1; 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::get_max_u_F_prime}  

// 另一个比较简单的辅助函数。计算总速度乘以分数流函数的导数的最大值，即计算  $\|\mathbf{u} F'(S)\|_{L_\infty(\Omega)}$  。这个项既用于时间步长的计算，也用于人工黏度中熵留项的正常化。

  template <int dim> 
  double TwoPhaseFlowProblem<dim>::get_max_u_F_prime() const 
  { 
    const QGauss<dim>  quadrature_formula(darcy_degree + 2); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values); 
    FEValues<dim> saturation_fe_values(saturation_fe, 
                                       quadrature_formula, 
                                       update_values); 

    std::vector<Vector<double>> darcy_solution_values(n_q_points, 
                                                      Vector<double>(dim + 1)); 
    std::vector<double>         saturation_values(n_q_points); 

    double max_velocity_times_dF_dS = 0; 

    auto       cell            = darcy_dof_handler.begin_active(); 
    const auto endc            = darcy_dof_handler.end(); 
    auto       saturation_cell = saturation_dof_handler.begin_active(); 
    for (; cell != endc; ++cell, ++saturation_cell) 
      { 
        darcy_fe_values.reinit(cell); 
        saturation_fe_values.reinit(saturation_cell); 

        darcy_fe_values.get_function_values(darcy_solution, 
                                            darcy_solution_values); 
        saturation_fe_values.get_function_values(old_saturation_solution, 
                                                 saturation_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            Tensor<1, dim> velocity; 
            for (unsigned int i = 0; i < dim; ++i) 
              velocity[i] = darcy_solution_values[q](i); 

            const double dF_dS = 
              fractional_flow_derivative(saturation_values[q], viscosity); 

            max_velocity_times_dF_dS = 
              std::max(max_velocity_times_dF_dS, velocity.norm() * dF_dS); 
          } 
      } 

    return max_velocity_times_dF_dS; 
  } 
// @sect4{TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range}  

// 为了计算稳定化项，我们需要知道饱和变量的范围。与 step-31 不同，这个范围很容易被区间 $[0,1]$ 所约束，但是我们可以通过在正交点的集合上循环，看看那里的值是多少，从而做得更好。如果可以的话，也就是说，如果周围至少有两个时间步长，我们甚至可以把这些值推算到下一个时间步长。

// 和以前一样，这个函数是在对  step-31  进行最小修改后取的。

  template <int dim> 
  std::pair<double, double> 
  TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range() const 
  { 
    const QGauss<dim>  quadrature_formula(saturation_degree + 2); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(saturation_fe, quadrature_formula, update_values); 
    std::vector<double> old_saturation_values(n_q_points); 
    std::vector<double> old_old_saturation_values(n_q_points); 

    if (timestep_number != 0) 
      { 
        double min_saturation = std::numeric_limits<double>::max(), 
               max_saturation = -std::numeric_limits<double>::max(); 

        for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
          { 
            fe_values.reinit(cell); 
            fe_values.get_function_values(old_saturation_solution, 
                                          old_saturation_values); 
            fe_values.get_function_values(old_old_saturation_solution, 
                                          old_old_saturation_values); 

            for (unsigned int q = 0; q < n_q_points; ++q) 
              { 
                const double saturation = 
                  (1. + time_step / old_time_step) * old_saturation_values[q] - 
                  time_step / old_time_step * old_old_saturation_values[q]; 

                min_saturation = std::min(min_saturation, saturation); 
                max_saturation = std::max(max_saturation, saturation); 
              } 
          } 

        return std::make_pair(min_saturation, max_saturation); 
      } 
    else 
      { 
        double min_saturation = std::numeric_limits<double>::max(), 
               max_saturation = -std::numeric_limits<double>::max(); 

        for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
          { 
            fe_values.reinit(cell); 
            fe_values.get_function_values(old_saturation_solution, 
                                          old_saturation_values); 

            for (unsigned int q = 0; q < n_q_points; ++q) 
              { 
                const double saturation = old_saturation_values[q]; 

                min_saturation = std::min(min_saturation, saturation); 
                max_saturation = std::max(max_saturation, saturation); 
              } 
          } 

        return std::make_pair(min_saturation, max_saturation); 
      } 
  } 

//  @sect4{TwoPhaseFlowProblem<dim>::compute_viscosity}  

// 最后一个工具函数是用来计算给定单元上的人工粘度的。如果你面前有它的公式，这并不特别复杂，看一下  step-31  中的实现。与那个教程程序的主要区别是，这里的速度不是简单的 $\mathbf u$ ，而是 $\mathbf u F'(S)$ ，一些公式需要做相应的调整。

  template <int dim> 
  double TwoPhaseFlowProblem<dim>::compute_viscosity( 
    const std::vector<double> &        old_saturation, 
    const std::vector<double> &        old_old_saturation, 
    const std::vector<Tensor<1, dim>> &old_saturation_grads, 
    const std::vector<Tensor<1, dim>> &old_old_saturation_grads, 
    const std::vector<Vector<double>> &present_darcy_values, 
    const double                       global_max_u_F_prime, 
    const double                       global_S_variation, 
    const double                       cell_diameter) const 
  { 
    const double beta  = .4 * dim; 
    const double alpha = 1; 

    if (global_max_u_F_prime == 0) 
      return 5e-3 * cell_diameter; 

    const unsigned int n_q_points = old_saturation.size(); 

    double max_residual             = 0; 
    double max_velocity_times_dF_dS = 0; 

    const bool use_dF_dS = true; 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        Tensor<1, dim> u; 
        for (unsigned int d = 0; d < dim; ++d) 
          u[d] = present_darcy_values[q](d); 

        const double dS_dt = porosity * 
                             (old_saturation[q] - old_old_saturation[q]) / 
                             old_time_step; 

        const double dF_dS = fractional_flow_derivative( 
          (old_saturation[q] + old_old_saturation[q]) / 2.0, viscosity); 

        const double u_grad_S = 
          u * dF_dS * (old_saturation_grads[q] + old_old_saturation_grads[q]) / 
          2.0; 

        const double residual = 
          std::abs((dS_dt + u_grad_S) * 
                   std::pow((old_saturation[q] + old_old_saturation[q]) / 2, 
                            alpha - 1.)); 

        max_residual = std::max(residual, max_residual); 
        max_velocity_times_dF_dS = 
          std::max(std::sqrt(u * u) * (use_dF_dS ? std::max(dF_dS, 1.) : 1), 
                   max_velocity_times_dF_dS); 
      } 

    const double c_R            = 1.0; 
    const double global_scaling = c_R * porosity * 
                                  (global_max_u_F_prime)*global_S_variation / 
                                  std::pow(global_Omega_diameter, alpha - 2.); 

    return (beta * 
            (max_velocity_times_dF_dS)*std::min(cell_diameter, 
                                                std::pow(cell_diameter, alpha) * 
                                                  max_residual / 
                                                  global_scaling)); 
  } 
// @sect3{TwoPhaseFlowProblem<dim>::run}  

// 除了 <code>solve()</code> 之外，这个函数是这个程序的主要功能，因为它控制了迭代的时间，以及何时将解决方案写入输出文件，何时进行网格细化。

// 除了启动代码通过 <code>goto start_time_iteration</code> 标签循环回到函数的开头外，一切都应该是相对简单的。无论如何，它模仿了  step-31  中的相应函数。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::run() 
  { 
    const unsigned int initial_refinement     = (dim == 2 ? 5 : 2); 
    const unsigned int n_pre_refinement_steps = (dim == 2 ? 3 : 2); 

    GridGenerator::hyper_cube(triangulation, 0, 1); 
    triangulation.refine_global(initial_refinement); 
    global_Omega_diameter = GridTools::diameter(triangulation); 

    setup_dofs(); 

    unsigned int pre_refinement_step = 0; 

  start_time_iteration: 

    VectorTools::project(saturation_dof_handler, 
                         saturation_constraints, 
                         QGauss<dim>(saturation_degree + 2), 
                         SaturationInitialValues<dim>(), 
                         old_saturation_solution); 

    time_step = old_time_step = 0; 
    current_macro_time_step = old_macro_time_step = 0; 

    time = 0; 

    do 
      { 
        std::cout << "Timestep " << timestep_number << ":  t=" << time 
                  << ", dt=" << time_step << std::endl; 

        solve(); 

        std::cout << std::endl; 

        if (timestep_number % 200 == 0) 
          output_results(); 

        if (timestep_number % 25 == 0) 
          refine_mesh(initial_refinement, 
                      initial_refinement + n_pre_refinement_steps); 

        if ((timestep_number == 0) && 
            (pre_refinement_step < n_pre_refinement_steps)) 
          { 
            ++pre_refinement_step; 
            goto start_time_iteration; 
          } 

        time += time_step; 
        ++timestep_number; 

        old_old_saturation_solution = old_saturation_solution; 
        old_saturation_solution     = saturation_solution; 
      } 
    while (time <= end_time); 
  } 
} // namespace Step43 

//  @sect3{The <code>main()</code> function}  

// 主函数看起来与所有其他程序几乎一样。对于使用Trilinos的程序来说，需要初始化MPI子系统--即使是那些实际上没有并行运行的程序--在  step-31  中有解释。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step43; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 

// 这个程序只能在串行中运行。否则，将抛出一个异常。

      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
                  ExcMessage( 
                    "This program can only be run in serial, use ./step-43")); 

      TwoPhaseFlowProblem<2> two_phase_flow_problem(1); 
      two_phase_flow_problem.run(); 
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


