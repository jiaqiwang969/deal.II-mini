

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2006 - 2021 by the deal.II authors 
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
 * Authors: Yan Li, Wolfgang Bangerth, Texas A&M University, 2006 
 */ 



// 这个程序是对  step-20  的改编，包括一些来自  step-12  的DG方法的技术。因此，该程序的很大一部分与  step-20  非常相似，我们将不再对这些部分进行评论。只有新的东西才会被详细讨论。

//  @sect3{Include files}  

// 这些include文件以前都用过了。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_raviart_thomas.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <iostream> 
#include <fstream> 

// 在这个程序中，我们使用一个张量值的系数。由于它可能具有空间依赖性，我们认为它是一个张量值的函数。下面的include文件提供了提供这种功能的 <code>TensorFunction</code> 类。

#include <deal.II/base/tensor_function.h> 

// 此外，我们使用 <code>DiscreteTime</code> 类来执行与时间递增有关的操作。

#include <deal.II/base/discrete_time.h> 

// 最后一步和以前所有的程序一样。

namespace Step21 
{ 
  using namespace dealii; 
// @sect3{The <code>TwoPhaseFlowProblem</code> class}  

// 这是该程序的主类。它与 step-20 中的类很接近，但增加了一些功能。

//  <ul>  
// <li>  
// <code>assemble_rhs_S</code> 集合了饱和度方程的右侧。正如介绍中所解释的，这不能被集成到 <code>assemble_rhs</code> 中，因为它取决于在时间步长的第一部分计算的速度。

//  <li>  
// <code>get_maximal_velocity</code> 的作用正如其名称所示。这个函数用于计算时间步长。

//  <li>  
// <code>project_back_saturation</code>  将所有饱和度小于0的自由度重置为0，所有饱和度大于1的自由度重置为1。   </ul>  

// 该类的其余部分应该是非常明显的。变量 <code>viscosity</code> 存储粘度 $\mu$ ，它进入了非线性方程中的几个公式。变量 <code>time</code> 记录了模拟过程中的时间信息。

  template <int dim> 
  class TwoPhaseFlowProblem 
  { 
  public: 
    TwoPhaseFlowProblem(const unsigned int degree); 
    void run(); 

  private: 
    void   make_grid_and_dofs(); 
    void   assemble_system(); 
    void   assemble_rhs_S(); 
    double get_maximal_velocity() const; 
    void   solve(); 
    void   project_back_saturation(); 
    void   output_results() const; 

    const unsigned int degree; 

    Triangulation<dim> triangulation; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    const unsigned int n_refinement_steps; 

    DiscreteTime time; 
    double       viscosity; 

    BlockVector<double> solution; 
    BlockVector<double> old_solution; 
    BlockVector<double> system_rhs; 
  }; 
// @sect3{Equation data}  
// @sect4{Pressure right hand side}  

// 目前，压力方程的右侧仅仅是零函数。但是，如果需要的话，程序的其余部分完全可以处理其他的东西。

  template <int dim> 
  class PressureRightHandSide : public Function<dim> 
  { 
  public: 
    PressureRightHandSide() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      return 0; 
    } 
  }; 

//  @sect4{Pressure boundary values}  

// 接下来是压力边界值。正如介绍中提到的，我们选择一个线性压力场。

  template <int dim> 
  class PressureBoundaryValues : public Function<dim> 
  { 
  public: 
    PressureBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      return 1 - p[0]; 
    } 
  }; 

//  @sect4{Saturation boundary values}  

// 然后，我们还需要边界的流入部分的边界值。某物是否为流入部分的问题是在组装右手边时决定的，我们只需要提供边界值的功能描述。这正如介绍中所解释的。

  template <int dim> 
  class SaturationBoundaryValues : public Function<dim> 
  { 
  public: 
    SaturationBoundaryValues() 
      : Function<dim>(1) 
    {} 

    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      if (p[0] == 0) 
        return 1; 
      else 
        return 0; 
    } 
  }; 

//  @sect4{Initial data}  

// 最后，我们需要初始数据。实际上，我们只需要饱和度的初始数据，但我们很懒，所以以后在第一个时间步骤之前，我们会简单地从一个包含所有矢量分量的函数中插值出前一个时间步骤的整个解决方案。
//因此，
//我们简单地创建一个所有分量都返回0的函数。我们通过简单地将每个函数转发到 Functions::ZeroFunction 类来做到这一点。为什么不在这个程序中我们目前使用 <code>InitialValues</code> 类的地方立即使用呢？因为这样，以后再回去选择不同的函数来做初始值就更简单了。

  template <int dim> 
  class InitialValues : public Function<dim> 
  { 
  public: 
    InitialValues() 
      : Function<dim>(dim + 2) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      return Functions::ZeroFunction<dim>(dim + 2).value(p, component); 
    } 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override 
    { 
      Functions::ZeroFunction<dim>(dim + 2).vector_value(p, values); 
    } 
  }; 

//  @sect3{The inverse permeability tensor}  

// 正如介绍中所宣布的，我们实现了两个不同的渗透率张量场。我们把它们各自放入一个命名空间，这样以后就可以很容易地在代码中用另一个来代替一个。

//  @sect4{Single curving crack permeability}  

// 渗透率的第一个函数是模拟单个弯曲裂缝的函数。它在 step-20 的结尾已经使用过了，它的函数形式在本教程程序的介绍中给出。和以前的一些程序一样，我们必须声明KInverse类的一个（似乎是不必要的）默认构造函数，以避免某些编译器的警告。

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
                 std::vector<Tensor<2, dim>> &  values) const override 
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
    }; 
  } // namespace SingleCurvingCrack 
// @sect4{Random medium permeability}  

// 这个函数的作用与介绍中公布的一样，即在随机的地方创建一个叠加的指数。对于这个类，有一件事值得考虑。这个问题的核心是，这个类使用随机函数创建指数的中心。如果我们因此在每次创建本类型的对象时都创建中心，我们每次都会得到一个不同的中心列表。这不是我们对这种类型的类的期望：它们应该可靠地表示同一个函数。

// 解决这个问题的方法是使中心列表成为这个类的静态成员变量，也就是说，在整个程序中只存在一个这样的变量，而不是为这个类型的每个对象。这正是我们所要做的。

// 然而，接下来的问题是，我们需要一种方法来初始化这个变量。由于这个变量是在程序开始时初始化的，我们不能使用普通的成员函数来实现，因为当时身边可能没有这个类型的对象。因此C++标准规定，只有非成员函数和静态成员函数可以用来初始化静态变量。我们通过定义一个函数 <code>get_centers</code> 来使用后一种可能性，该函数在调用时计算中心点的列表。

// 注意，这个类在2D和3D中都能正常工作，唯一的区别是我们在3D中使用了更多的点：通过实验我们发现，我们在3D中比2D中需要更多的指数（毕竟我们有更多的地方需要覆盖，如果我们想保持中心之间的距离大致相等），所以我们在2D中选择40，在3D中选择100。对于任何其他维度，该函数目前不知道该怎么做，所以只是抛出一个异常，表明这一点。

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
                 std::vector<Tensor<2, dim>> &  values) const override 
      { 
        Assert(points.size() == values.size(), 
               ExcDimensionMismatch(points.size(), values.size())); 

        for (unsigned int p = 0; p < points.size(); ++p) 
          { 
            values[p].clear(); 

            double permeability = 0; 
            for (unsigned int i = 0; i < centers.size(); ++i) 
              permeability += std::exp(-(points[p] - centers[i]).norm_square() / 
                                       (0.05 * 0.05)); 

            const double normalized_permeability = 
              std::min(std::max(permeability, 0.01), 4.); 

            for (unsigned int d = 0; d < dim; ++d) 
              values[p][d][d] = 1. / normalized_permeability; 
          } 
      } 

    private: 
      static std::vector<Point<dim>> centers; 

      static std::vector<Point<dim>> get_centers() 
      { 
        const unsigned int N = 
          (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented())); 

        std::vector<Point<dim>> centers_list(N); 
        for (unsigned int i = 0; i < N; ++i) 
          for (unsigned int d = 0; d < dim; ++d) 
            centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX; 

        return centers_list; 
      } 
    }; 

    template <int dim> 
    std::vector<Point<dim>> 
      KInverse<dim>::centers = KInverse<dim>::get_centers(); 
  } // namespace RandomMedium 

//  @sect3{The inverse mobility and saturation functions}  

// 还有两个数据我们需要描述，即反流动性函数和饱和度曲线。它们的形式也在介绍中给出。

  double mobility_inverse(const double S, const double viscosity) 
  { 
    return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S)); 
  } 

  double fractional_flow(const double S, const double viscosity) 
  { 
    return S * S / (S * S + viscosity * (1 - S) * (1 - S)); 
  } 

//  @sect3{Linear solvers and preconditioners}  

// 我们使用的线性求解器也完全类似于  step-20  中使用的。因此，下面的类是逐字逐句从那里复制过来的。请注意，这里的类不仅是从 step-20 中复制的，而且在deal.II中也有重复的类。在这个例子的未来版本中，它们应该被一个有效的方法所取代，不过。有一个变化：如果线性系统的尺寸很小，即当网格很粗时，那么在 <code>src.size()</code> 函数中的求解器收敛之前，设置 <code>vmult()</code> CG迭代的最大值有时是不够的。(当然，这是数值取舍的结果，因为我们知道在纸面上，CG方法最多在 <code>src.size()</code> 步内收敛)。因此，我们将最大的迭代次数设定为等于线性系统的最大规模和200。

  template <class MatrixType> 
  class InverseMatrix : public Subscriptor 
  { 
  public: 
    InverseMatrix(const MatrixType &m) 
      : matrix(&m) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      SolverControl solver_control(std::max<unsigned int>(src.size(), 200), 
                                   1e-8 * src.l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

      dst = 0; 

      cg.solve(*matrix, dst, src, PreconditionIdentity()); 
    } 

  private: 
    const SmartPointer<const MatrixType> matrix; 
  }; 

  class SchurComplement : public Subscriptor 
  { 
  public: 
    SchurComplement(const BlockSparseMatrix<double> &          A, 
                    const InverseMatrix<SparseMatrix<double>> &Minv) 
      : system_matrix(&A) 
      , m_inverse(&Minv) 
      , tmp1(A.block(0, 0).m()) 
      , tmp2(A.block(0, 0).m()) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      system_matrix->block(0, 1).vmult(tmp1, src); 
      m_inverse->vmult(tmp2, tmp1); 
      system_matrix->block(1, 0).vmult(dst, tmp2); 
    } 

  private: 
    const SmartPointer<const BlockSparseMatrix<double>>           system_matrix; 
    const SmartPointer<const InverseMatrix<SparseMatrix<double>>> m_inverse; 

    mutable Vector<double> tmp1, tmp2; 
  }; 

  class ApproximateSchurComplement : public Subscriptor 
  { 
  public: 
    ApproximateSchurComplement(const BlockSparseMatrix<double> &A) 
      : system_matrix(&A) 
      , tmp1(A.block(0, 0).m()) 
      , tmp2(A.block(0, 0).m()) 
    {} 

    void vmult(Vector<double> &dst, const Vector<double> &src) const 
    { 
      system_matrix->block(0, 1).vmult(tmp1, src); 
      system_matrix->block(0, 0).precondition_Jacobi(tmp2, tmp1); 
      system_matrix->block(1, 0).vmult(dst, tmp2); 
    } 

  private: 
    const SmartPointer<const BlockSparseMatrix<double>> system_matrix; 

    mutable Vector<double> tmp1, tmp2; 
  }; 

//  @sect3{<code>TwoPhaseFlowProblem</code> class implementation}  

// 现在是主类的实现。它的大部分内容实际上是从  step-20  中复制过来的，所以我们不会对它进行详细的评论。你应该试着先熟悉一下那个程序，然后这里发生的大部分事情就应该很清楚了。

//  @sect4{TwoPhaseFlowProblem::TwoPhaseFlowProblem}  

// 首先是构造函数。我们使用 $RT_k \times DQ_k \times DQ_k$ 空间。对于初始化DiscreteTime对象，我们不在构造函数中设置时间步长，因为我们还没有它的值。时间步长最初被设置为零，但在需要增量时间之前，它将被计算出来，正如介绍的一个小节中所描述的。时间对象在内部阻止自己在 $dt = 0$ 时被递增，迫使我们在推进时间之前为 $dt$ 设置一个非零的期望大小。

  template <int dim> 
  TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree) 
    : degree(degree) 
    , fe(FE_RaviartThomas<dim>(degree), 
         1, 
         FE_DGQ<dim>(degree), 
         1, 
         FE_DGQ<dim>(degree), 
         1) 
    , dof_handler(triangulation) 
    , n_refinement_steps(5) 
    , time(/*start time*/ 0., /*end time*/ 1.) 
    , viscosity(0.2) 
  {} 

//  @sect4{TwoPhaseFlowProblem::make_grid_and_dofs}  

// 下一个函数从众所周知的函数调用开始，创建和细化一个网格，然后将自由度与之关联。它所做的事情与 step-20 中的相同，只是现在是三个组件而不是两个。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, 0, 1); 
    triangulation.refine_global(n_refinement_steps); 

    dof_handler.distribute_dofs(fe); 
    DoFRenumbering::component_wise(dof_handler); 

    const std::vector<types::global_dof_index> dofs_per_component = 
      DoFTools::count_dofs_per_fe_component(dof_handler); 
    const unsigned int n_u = dofs_per_component[0], 
                       n_p = dofs_per_component[dim], 
                       n_s = dofs_per_component[dim + 1]; 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (" << n_u << '+' << n_p << '+' << n_s << ')' << std::endl 
              << std::endl; 

    const unsigned int n_couplings = dof_handler.max_couplings_between_dofs(); 

    sparsity_pattern.reinit(3, 3); 
    sparsity_pattern.block(0, 0).reinit(n_u, n_u, n_couplings); 
    sparsity_pattern.block(1, 0).reinit(n_p, n_u, n_couplings); 
    sparsity_pattern.block(2, 0).reinit(n_s, n_u, n_couplings); 
    sparsity_pattern.block(0, 1).reinit(n_u, n_p, n_couplings); 
    sparsity_pattern.block(1, 1).reinit(n_p, n_p, n_couplings); 
    sparsity_pattern.block(2, 1).reinit(n_s, n_p, n_couplings); 
    sparsity_pattern.block(0, 2).reinit(n_u, n_s, n_couplings); 
    sparsity_pattern.block(1, 2).reinit(n_p, n_s, n_couplings); 
    sparsity_pattern.block(2, 2).reinit(n_s, n_s, n_couplings); 

    sparsity_pattern.collect_sizes(); 

    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern); 
    sparsity_pattern.compress(); 

    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(3); 
    solution.block(0).reinit(n_u); 
    solution.block(1).reinit(n_p); 
    solution.block(2).reinit(n_s); 
    solution.collect_sizes(); 

    old_solution.reinit(3); 
    old_solution.block(0).reinit(n_u); 
    old_solution.block(1).reinit(n_p); 
    old_solution.block(2).reinit(n_s); 
    old_solution.collect_sizes(); 

    system_rhs.reinit(3); 
    system_rhs.block(0).reinit(n_u); 
    system_rhs.block(1).reinit(n_p); 
    system_rhs.block(2).reinit(n_s); 
    system_rhs.collect_sizes(); 
  } 
// @sect4{TwoPhaseFlowProblem::assemble_system}  

// 这是组装线性系统的函数，或者至少是除了(1,3)块之外的所有东西，它取决于在这个时间步长中计算的仍然未知的速度（我们在 <code>assemble_rhs_S</code> 中处理这个问题）。它的大部分内容与 step-20 一样，但这次我们必须处理一些非线性的问题。 然而，该函数的顶部与往常一样（注意我们在开始时将矩阵和右手边设置为零&mdash; 对于静止问题我们不必这样做，因为在那里我们只使用一次矩阵对象，而且在开始时它是空的）。

// 注意，在目前的形式下，该函数使用 RandomMedium::KInverse 类中实现的渗透率。切换到单曲裂缝渗透率函数就像改变命名空间名称一样简单。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_system() 
  { 
    system_matrix = 0; 
    system_rhs    = 0; 

    QGauss<dim>     quadrature_formula(degree + 2); 
    QGauss<dim - 1> face_quadrature_formula(degree + 2); 

    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const PressureRightHandSide<dim>  pressure_right_hand_side; 
    const PressureBoundaryValues<dim> pressure_boundary_values; 
    const RandomMedium::KInverse<dim> k_inverse; 

    std::vector<double>         pressure_rhs_values(n_q_points); 
    std::vector<double>         boundary_values(n_face_q_points); 
    std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 

    std::vector<Vector<double>>              old_solution_values(n_q_points, 
                                                                 Vector<double>(dim + 2)); 
    std::vector<std::vector<Tensor<1, dim>>> old_solution_grads( 
      n_q_points, std::vector<Tensor<1, dim>>(dim + 2)); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 
    const FEValuesExtractors::Scalar saturation(dim + 1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        local_matrix = 0; 
        local_rhs    = 0; 

// 这里是第一个重要的区别。我们必须在正交点上获得前一个时间步骤的饱和函数值。为此，我们可以使用 FEValues::get_function_values （之前已经在 step-9 、 step-14 和 step-15 中使用），这个函数接收一个解向量并返回当前单元的正交点的函数值列表。事实上，它返回每个正交点的完整矢量值解，即不仅是饱和度，还有速度和压力。

        fe_values.get_function_values(old_solution, old_solution_values); 

// 然后，我们还必须得到压力的右手边和反渗透性张量在正交点的数值。

        pressure_right_hand_side.value_list(fe_values.get_quadrature_points(), 
                                            pressure_rhs_values); 
        k_inverse.value_list(fe_values.get_quadrature_points(), 
                             k_inverse_values); 

// 有了这些，我们现在可以在这个单元格上的所有正交点和形状函数上进行循环，并将我们在这个函数中处理的矩阵和右手边的那些部分组合起来。考虑到引言中所述的双线性形式的明确形式，贡献中的各个条款应该是不言自明的。

        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const double old_s = old_solution_values[q](dim + 1); 

              const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q); 
              const double div_phi_i_u = fe_values[velocities].divergence(i, q); 
              const double phi_i_p     = fe_values[pressure].value(i, q); 
              const double phi_i_s     = fe_values[saturation].value(i, q); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const Tensor<1, dim> phi_j_u = 
                    fe_values[velocities].value(j, q); 
                  const double div_phi_j_u = 
                    fe_values[velocities].divergence(j, q); 
                  const double phi_j_p = fe_values[pressure].value(j, q); 
                  const double phi_j_s = fe_values[saturation].value(j, q); 

                  local_matrix(i, j) += 
                    (phi_i_u * k_inverse_values[q] * 
                       mobility_inverse(old_s, viscosity) * phi_j_u - 
                     div_phi_i_u * phi_j_p - phi_i_p * div_phi_j_u + 
                     phi_i_s * phi_j_s) * 
                    fe_values.JxW(q); 
                } 

              local_rhs(i) += 
                (-phi_i_p * pressure_rhs_values[q]) * fe_values.JxW(q); 
            } 

// 接下来，我们还必须处理压力边界值。这一点，还是和 step-20 中一样。

        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              fe_face_values.reinit(cell, face); 

              pressure_boundary_values.value_list( 
                fe_face_values.get_quadrature_points(), boundary_values); 

              for (unsigned int q = 0; q < n_face_q_points; ++q) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  { 
                    const Tensor<1, dim> phi_i_u = 
                      fe_face_values[velocities].value(i, q); 

                    local_rhs(i) += 
                      -(phi_i_u * fe_face_values.normal_vector(q) * 
                        boundary_values[q] * fe_face_values.JxW(q)); 
                  } 
            } 

// 在所有单元的循环中，最后一步是将局部贡献转移到全局矩阵和右侧向量中。

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            system_matrix.add(local_dof_indices[i], 
                              local_dof_indices[j], 
                              local_matrix(i, j)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += local_rhs(i); 
      } 
  } 

// 矩阵和右手边的组装就这么多了。请注意，我们不需要插值和应用边界值，因为它们都已经在弱式中被处理过了。

//  @sect4{TwoPhaseFlowProblem::assemble_rhs_S}  

// 正如在介绍中所解释的，我们只有在计算出速度后才能评估饱和方程的右边。因此，我们有这个单独的函数来实现这个目的。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::assemble_rhs_S() 
  { 
    QGauss<dim>       quadrature_formula(degree + 2); 
    QGauss<dim - 1>   face_quadrature_formula(degree + 2); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    FEFaceValues<dim> fe_face_values(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_normal_vectors | 
                                       update_quadrature_points | 
                                       update_JxW_values); 
    FEFaceValues<dim> fe_face_values_neighbor(fe, 
                                              face_quadrature_formula, 
                                              update_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    Vector<double> local_rhs(dofs_per_cell); 

    std::vector<Vector<double>> old_solution_values(n_q_points, 
                                                    Vector<double>(dim + 2)); 
    std::vector<Vector<double>> old_solution_values_face(n_face_q_points, 
                                                         Vector<double>(dim + 
                                                                        2)); 
    std::vector<Vector<double>> old_solution_values_face_neighbor( 
      n_face_q_points, Vector<double>(dim + 2)); 
    std::vector<Vector<double>> present_solution_values(n_q_points, 
                                                        Vector<double>(dim + 
                                                                       2)); 
    std::vector<Vector<double>> present_solution_values_face( 
      n_face_q_points, Vector<double>(dim + 2)); 

    std::vector<double>                  neighbor_saturation(n_face_q_points); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    SaturationBoundaryValues<dim> saturation_boundary_values; 

    const FEValuesExtractors::Scalar saturation(dim + 1); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        local_rhs = 0; 
        fe_values.reinit(cell); 

        fe_values.get_function_values(old_solution, old_solution_values); 
        fe_values.get_function_values(solution, present_solution_values); 

// 首先是单元格条款。按照介绍中的公式，这些是  $(S^n,\sigma)-(F(S^n) \mathbf{v}^{n+1},\nabla \sigma)$  ，其中  $\sigma$  是测试函数的饱和成分。

        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const double   old_s = old_solution_values[q](dim + 1); 
              Tensor<1, dim> present_u; 
              for (unsigned int d = 0; d < dim; ++d) 
                present_u[d] = present_solution_values[q](d); 

              const double         phi_i_s = fe_values[saturation].value(i, q); 
              const Tensor<1, dim> grad_phi_i_s = 
                fe_values[saturation].gradient(i, q); 

              local_rhs(i) += 
                (time.get_next_step_size() * fractional_flow(old_s, viscosity) * 
                   present_u * grad_phi_i_s + 
                 old_s * phi_i_s) * 
                fe_values.JxW(q); 
            } 

// 其次，我们必须处理面的边界上的通量部分。这就有点麻烦了，因为我们首先要确定哪些是细胞边界的流入和流出部分。如果我们有一个流入的边界，我们需要评估面的另一边的饱和度（或者边界值，如果我们在域的边界上）。

// 所有这些都有点棘手，但在  step-9  中已经有了一些详细的解释。请看这里，这应该是如何工作的!

        for (const auto face_no : cell->face_indices()) 
          { 
            fe_face_values.reinit(cell, face_no); 

            fe_face_values.get_function_values(old_solution, 
                                               old_solution_values_face); 
            fe_face_values.get_function_values(solution, 
                                               present_solution_values_face); 

            if (cell->at_boundary(face_no)) 
              saturation_boundary_values.value_list( 
                fe_face_values.get_quadrature_points(), neighbor_saturation); 
            else 
              { 
                const auto         neighbor = cell->neighbor(face_no); 
                const unsigned int neighbor_face = 
                  cell->neighbor_of_neighbor(face_no); 

                fe_face_values_neighbor.reinit(neighbor, neighbor_face); 

                fe_face_values_neighbor.get_function_values( 
                  old_solution, old_solution_values_face_neighbor); 

                for (unsigned int q = 0; q < n_face_q_points; ++q) 
                  neighbor_saturation[q] = 
                    old_solution_values_face_neighbor[q](dim + 1); 
              } 

            for (unsigned int q = 0; q < n_face_q_points; ++q) 
              { 
                Tensor<1, dim> present_u_face; 
                for (unsigned int d = 0; d < dim; ++d) 
                  present_u_face[d] = present_solution_values_face[q](d); 

                const double normal_flux = 
                  present_u_face * fe_face_values.normal_vector(q); 

                const bool is_outflow_q_point = (normal_flux >= 0); 

                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  local_rhs(i) -= 
                    time.get_next_step_size() * normal_flux * 
                    fractional_flow((is_outflow_q_point == true ? 
                                       old_solution_values_face[q](dim + 1) : 
                                       neighbor_saturation[q]), 
                                    viscosity) * 
                    fe_face_values[saturation].value(i, q) * 
                    fe_face_values.JxW(q); 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          system_rhs(local_dof_indices[i]) += local_rhs(i); 
      } 
  } 

//  @sect4{TwoPhaseFlowProblem::solve}  

// 在所有这些准备工作之后，我们最终以与  step-20  相同的方式解决速度和压力的线性系统。在这之后，我们必须处理饱和方程（见下文）。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::solve() 
  { 
    const InverseMatrix<SparseMatrix<double>> m_inverse( 
      system_matrix.block(0, 0)); 
    Vector<double> tmp(solution.block(0).size()); 
    Vector<double> schur_rhs(solution.block(1).size()); 
    Vector<double> tmp2(solution.block(2).size()); 

// 首先是压力，使用前两个方程的压力舒尔补。

    { 
      m_inverse.vmult(tmp, system_rhs.block(0)); 
      system_matrix.block(1, 0).vmult(schur_rhs, tmp); 
      schur_rhs -= system_rhs.block(1); 

      SchurComplement schur_complement(system_matrix, m_inverse); 

      ApproximateSchurComplement approximate_schur_complement(system_matrix); 

      InverseMatrix<ApproximateSchurComplement> preconditioner( 
        approximate_schur_complement); 

      SolverControl            solver_control(solution.block(1).size(), 
                                   1e-12 * schur_rhs.l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

      cg.solve(schur_complement, solution.block(1), schur_rhs, preconditioner); 

      std::cout << "   " << solver_control.last_step() 
                << " CG Schur complement iterations for pressure." << std::endl; 
    } 

// 现在是速度。

    { 
      system_matrix.block(0, 1).vmult(tmp, solution.block(1)); 
      tmp *= -1; 
      tmp += system_rhs.block(0); 

      m_inverse.vmult(solution.block(0), tmp); 
    } 

// 最后，我们必须处理好饱和度方程。在这里，我们要做的第一件事是使用介绍中的公式来确定时间步长。知道了我们领域的形状，以及我们通过有规律地划分单元来创建网格，我们可以很容易地计算出每个单元的直径（事实上我们使用的是单元坐标方向上的线性扩展，而不是直径）。请注意，我们将在 step-24 中学习一种更通用的方法，在那里我们使用 GridTools::minimal_cell_diameter 函数。

// 我们使用一个辅助函数来计算下面定义的最大速度，有了这些，我们就可以评估我们新的时间步长了。我们使用方法 DiscreteTime::set_desired_next_time_step() 来向DiscreteTime对象建议新的时间步长的计算值。在大多数情况下，时间对象使用精确提供的值来增加时间。在某些情况下，时间对象可以进一步修改步骤大小。例如，如果计算出的时间增量超过了结束时间，它将被相应地截断。

    time.set_desired_next_step_size(std::pow(0.5, double(n_refinement_steps)) / 
                                    get_maximal_velocity()); 

// 下一步是组装右手边，然后把所有的东西都传给解。最后，我们把饱和度投射回物理上合理的范围。

    assemble_rhs_S(); 
    { 
      SolverControl            solver_control(system_matrix.block(2, 2).m(), 
                                   1e-8 * system_rhs.block(2).l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 
      cg.solve(system_matrix.block(2, 2), 
               solution.block(2), 
               system_rhs.block(2), 
               PreconditionIdentity()); 

      project_back_saturation(); 

      std::cout << "   " << solver_control.last_step() 
                << " CG iterations for saturation." << std::endl; 
    } 

    old_solution = solution; 
  } 
// @sect4{TwoPhaseFlowProblem::output_results}  

// 这里没有什么值得惊讶的。由于程序会做大量的时间步骤，我们只在每第五个时间步骤创建一个输出文件，并在文件的顶部已经跳过所有其他时间步骤。

// 在为接近函数底部的输出创建文件名时，我们将时间步长的数字转换为字符串表示，用前导零填充到四位数。我们这样做是因为这样所有的输出文件名都有相同的长度，因此在创建目录列表时可以很好地排序。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::output_results() const 
  { 
    if (time.get_step_number() % 5 != 0) 
      return; 

    std::vector<std::string> solution_names; 
    switch (dim) 
      { 
        case 2: 
          solution_names = {"u", "v", "p", "S"}; 
          break; 

        case 3: 
          solution_names = {"u", "v", "w", "p", "S"}; 
          break; 

        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, solution_names); 

    data_out.build_patches(degree + 1); 

    std::ofstream output("solution-" + 
                         Utilities::int_to_string(time.get_step_number(), 4) + 
                         ".vtk"); 
    data_out.write_vtk(output); 
  } 

//  @sect4{TwoPhaseFlowProblem::project_back_saturation}  

// 在这个函数中，我们简单地遍历所有的饱和自由度，并确保如果它们离开了物理上的合理范围，它们将被重置到区间  $[0,1]$  。要做到这一点，我们只需要循环解决向量的所有饱和分量；这些分量存储在块2中（块0是速度，块1是压力）。

// 值得注意的是，当时间步长选择如介绍中提到的那样时，这个函数几乎从未触发过，这一点可能很有启发。然而，如果我们只选择稍大的时间步长，我们会得到大量超出适当范围的数值。严格来说，如果我们选择的时间步长足够小，这个函数因此是不必要的。从某种意义上说，这个函数只是一个安全装置，以避免由于个别自由度在几个时间步长之前变得不符合物理条件而导致我们的整个解决方案变得不符合物理条件的情况。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::project_back_saturation() 
  { 
    for (unsigned int i = 0; i < solution.block(2).size(); ++i) 
      if (solution.block(2)(i) < 0) 
        solution.block(2)(i) = 0; 
      else if (solution.block(2)(i) > 1) 
        solution.block(2)(i) = 1; 
  } 
// @sect4{TwoPhaseFlowProblem::get_maximal_velocity}  

// 下面的函数用于确定允许的最大时间步长。它的作用是在域中的所有正交点上循环，找出速度的最大幅度。

  template <int dim> 
  double TwoPhaseFlowProblem<dim>::get_maximal_velocity() const 
  { 
    QGauss<dim>        quadrature_formula(degree + 2); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(fe, quadrature_formula, update_values); 
    std::vector<Vector<double>> solution_values(n_q_points, 
                                                Vector<double>(dim + 2)); 
    double                      max_velocity = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        fe_values.get_function_values(solution, solution_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            Tensor<1, dim> velocity; 
            for (unsigned int i = 0; i < dim; ++i) 
              velocity[i] = solution_values[q](i); 

            max_velocity = std::max(max_velocity, velocity.norm()); 
          } 
      } 

    return max_velocity; 
  } 
// @sect4{TwoPhaseFlowProblem::run}  

// 这是我们主类的最后一个函数。它的简洁不言自明。只有两点是值得注意的。首先，该函数在开始时将初始值投射到有限元空间上； VectorTools::project 函数这样做需要一个表明悬挂节点约束的参数。我们在这个程序中没有（我们在一个均匀细化的网格上计算），但是这个函数当然需要这个参数。所以我们必须创建一个约束对象。在原始状态下，约束对象是没有排序的，在使用前必须进行排序（使用 AffineConstraints::close 函数）。这就是我们在这里所做的，这也是为什么我们不能简单地用一个匿名的临时对象 <code>AffineConstraints<double>()</code> 作为第二个参数来调用 VectorTools::project 函数。

// 值得一提的第二点是，我们只在求解每个时间步长对应的线性系统的过程中计算当前时间步长。因此，我们只有在时间步长结束时才能输出一个时间步长的当前时间。我们通过调用循环内的方法 DiscreteTime::advance_time() 来增加时间。由于我们在增量后报告时间和dt，我们必须调用方法 DiscreteTime::get_previous_step_size() ，而不是 DiscreteTime::get_next_step_size(). 。 经过许多步，当模拟到达结束时间时，最后的dt由DiscreteTime类选择，其方式是最后一步正好在结束时间完成。

  template <int dim> 
  void TwoPhaseFlowProblem<dim>::run() 
  { 
    make_grid_and_dofs(); 

    { 
      AffineConstraints<double> constraints; 
      constraints.close(); 

      VectorTools::project(dof_handler, 
                           constraints, 
                           QGauss<dim>(degree + 2), 
                           InitialValues<dim>(), 
                           old_solution); 
    } 

    do 
      { 
        std::cout << "Timestep " << time.get_step_number() + 1 << std::endl; 

        assemble_system(); 

        solve(); 

        output_results(); 

        time.advance_time(); 
        std::cout << "   Now at t=" << time.get_current_time() 
                  << ", dt=" << time.get_previous_step_size() << '.' 
                  << std::endl 
                  << std::endl; 
      } 
    while (time.is_at_end() == false); 
  } 
} // namespace Step21 
// @sect3{The <code>main</code> function}  

// 这就是了。在主函数中，我们将有限元空间的度数传递给TwoPhaseFlowProblem对象的构造函数。 这里，我们使用零度元素，即 $RT_0\times DQ_0 \times DQ_0$  。其余部分与其他所有程序一样。

int main() 
{ 
  try 
    { 
      using namespace Step21; 

      TwoPhaseFlowProblem<2> two_phase_flow_problem(0); 
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


