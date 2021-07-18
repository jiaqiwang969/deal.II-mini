

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
 * Authors: Martin Kronbichler, Uppsala University, 
 *          Wolfgang Bangerth, Texas A&M University 2007, 2008 
 */ 


// @sect3{Include files}  

// 像往常一样，第一步是包括这些著名的deal.II库文件和一些C++头文件的功能。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/block_sparsity_pattern.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/solution_transfer.h> 

// 然后我们需要包括一些头文件，这些文件提供了矢量、矩阵和预处理类，这些类实现了各自Trilinos类的接口。特别是，我们将需要基于Trilinos的矩阵和向量类以及Trilinos预处理程序的接口。

#include <deal.II/base/index_set.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_block_sparse_matrix.h> 
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_parallel_block_vector.h> 
#include <deal.II/lac/trilinos_precondition.h> 

// 最后，这里有几个C++头文件还没有被上述头文件中的某个文件所包含。

#include <iostream> 
#include <fstream> 
#include <memory> 
#include <limits> 

// 在这个顶层事项的最后，我们将所有deal.II的名字导入到全局命名空间。

namespace Step31 
{ 
  using namespace dealii; 
// @sect3{Equation data}  

// 同样，程序的下一阶段是定义方程数据，即各种边界条件、右手边和初始条件（记住，我们要解决的是一个时间依赖型系统）。这个定义的基本策略与  step-22  中的相同。不过关于细节，还是有一些区别。

// 首先，我们没有在速度上设置任何不均匀的边界条件，因为正如介绍中所解释的，我们将使用无流条件  $\mathbf{n}\cdot\mathbf{u}=0$  。所以剩下的是应力张量法线分量的切向部分的条件 <code>dim-1</code> ， $\textbf{n} \cdot [p \textbf{1} - \eta\varepsilon(\textbf{u})]$ ；我们假定这些分量的值是同质的，也就是说，一个自然的边界条件，不需要具体的动作（它作为零项出现在弱形式的右边）。

// 对于温度  $T$  ，我们假设没有热能通量，即  $\mathbf{n} \cdot \kappa \nabla T=0$  。这也是一个边界条件，不需要我们做任何特别的事情。

// 第二，我们必须设定温度的初始条件（速度和压力不需要初始条件，因为我们在这里考虑的准稳态情况下的斯托克斯方程没有速度或压力的时间导数）。在这里，我们选择一个非常简单的测试案例，即初始温度为零，所有的动力学都由温度的右手边驱动。

// 第三，我们需要定义温度方程的右边。我们选择它在域的底部某处的三个圆（或三维球）内为常数，如介绍中所解释的那样，而在域外为零。

// 最后，或者说首先，在这个命名空间的顶部，我们定义我们需要的各种材料常数（ $\eta,\kappa$ ，密度 $\rho$ 和热膨胀系数 $\beta$ ）。

  namespace EquationData 
  { 
    constexpr double eta     = 1; 
    constexpr double kappa   = 1e-6; 
    constexpr double beta    = 10; 
    constexpr double density = 1; 

    template <int dim> 
    class TemperatureInitialValues : public Function<dim> 
    { 
    public: 
      TemperatureInitialValues() 
        : Function<dim>(1) 
      {} 

      virtual double value(const Point<dim> & /*p*/, 
                           const unsigned int /*component*/ = 0) const override 
      { 
        return 0; 
      } 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  value) const override 
      { 
        for (unsigned int c = 0; c < this->n_components; ++c) 
          value(c) = TemperatureInitialValues<dim>::value(p, c); 
      } 
    }; 

    template <int dim> 
    class TemperatureRightHandSide : public Function<dim> 
    { 
    public: 
      TemperatureRightHandSide() 
        : Function<dim>(1) 
      {} 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override 
      { 
        (void)component; 
        Assert(component == 0, 
               ExcMessage("Invalid operation for a scalar function.")); 

        Assert((dim == 2) || (dim == 3), ExcNotImplemented()); 

        static const Point<dim> source_centers[3] = { 
          (dim == 2 ? Point<dim>(.3, .1) : Point<dim>(.3, .5, .1)), 
          (dim == 2 ? Point<dim>(.45, .1) : Point<dim>(.45, .5, .1)), 
          (dim == 2 ? Point<dim>(.75, .1) : Point<dim>(.75, .5, .1))}; 
        static const double source_radius = (dim == 2 ? 1. / 32 : 1. / 8); 

        return ((source_centers[0].distance(p) < source_radius) || 
                    (source_centers[1].distance(p) < source_radius) || 
                    (source_centers[2].distance(p) < source_radius) ? 
                  1 : 
                  0); 
      } 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  value) const override 
      { 
        for (unsigned int c = 0; c < this->n_components; ++c) 
          value(c) = TemperatureRightHandSide<dim>::value(p, c); 
      } 
    }; 
  } // namespace EquationData 

//  @sect3{Linear solvers and preconditioners}  

// 本节介绍了一些用于求解斯托克斯系统线性方程的对象，我们需要在每个时间步长中求解。这里使用的许多想法与 step-20 相同，其中介绍了基于Schur补的预处理程序和求解器，实际接口来自 step-22 （特别是 step-22 中 "结果 "部分的讨论，其中我们介绍了直接Schur补方法的替代品）。但是请注意，在这里我们不使用Schur补数来解决Stokes方程，尽管预处理程序中出现了一个近似的Schur补数（压力空间的质量矩阵）。

  namespace LinearSolvers 
  { 
// @sect4{The <code>InverseMatrix</code> class template}  

// 这个类是一个接口，用于计算 "倒置 "矩阵对向量的作用（使用 <code>vmult</code> 操作），其方式与 step-22 中的相应类相同：当请求这个类的对象的乘积时，我们使用CG方法解决与该矩阵有关的线性方程组，通过（模板化） <code>PreconditionerType</code> 类的预处理器加速。

// 与 step-22 中同一类别的实现略有不同，我们让 <code>vmult</code> 函数接受任何类型的向量类型（但是，如果矩阵不允许与这种向量进行矩阵-向量乘积，它将产生编译器错误）。

// 第二，我们捕捉解算器可能抛出的任何异常。原因如下。在调试这样的程序时，偶尔会犯一个错误，即把一个不确定或不对称的矩阵或预处理程序传递给当前的类。在这种情况下，求解器将不能收敛并抛出一个运行时异常。如果在这里没有被捕捉到，它就会在调用堆栈中传播，最后可能会在 <code>main()</code> 中出现，在那里我们会输出一个错误信息，说CG求解器失败。那么问题来了。哪个CG求解器？倒置质量矩阵的那个？用拉普拉斯算子反转左上角块的那个？还是在当前代码中我们使用线性求解器的其他几个嵌套位置中的一个CG求解器？在运行时异常中没有这方面的指示，因为它没有存储我们到达产生异常的地方的调用栈。
//所以
//与其让异常自由传播到 <code>main()</code> ，不如意识到如果内部求解器失败，外部函数能做的很少，不如将运行时异常转化为一个断言，该断言失败后会触发对 <code>abort()</code> 的调用，允许我们在调试器中追溯我们如何到达当前位置。

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
// @sect4{Schur complement preconditioner}  

// 这是在介绍中详细描述的舒尔补码预处理程序的实现。与 step-20 和 step-22 相反，我们使用GMRES一次性解决块系统，并使用块结构矩阵的Schur补码来建立一个良好的预处理程序。

// 让我们看看介绍中描述的理想预处理矩阵  $P=\left(\begin{array}{cc} A & 0 \\ B & -S \end{array}\right)$  。如果我们在线性系统的求解中应用这个矩阵，迭代式GMRES求解器的收敛性将受矩阵
// @f{eqnarray*} P^{-1}\left(\begin{array}{cc} A &
//  B^T \\ B & 0 \end{array}\right) = \left(\begin{array}{cc} I & A^{-1}
//  B^T \\ 0 & I \end{array}\right), 
//  @f}
//  的制约，这确实非常简单。基于精确矩阵的GMRES求解器将在一次迭代中收敛，因为所有的特征值都是相等的（任何Krylov方法最多需要多少次迭代就有多少个不同的特征值）。Silvester和Wathen提出了这样一个用于受阻斯托克斯系统的预处理程序（"稳定的斯托克斯系统的快速迭代解第二部分。 Using general block preconditioners", SIAM J. Numer. Anal., 31 (1994), pp.1352-1367）。)

//用 $\tilde{P}$ 代替 $P$ 可以保持这种精神：乘积 $P^{-1} A$ 仍将接近于特征值为1的矩阵，其分布不取决于问题大小。这让我们希望能够得到一个与问题规模无关的GMRES迭代次数。

// 已经通过 step-20 和 step-22 教程的deal.II用户当然可以想象我们将如何实现这一点。 我们用一些由InverseMatrix类构建的近似逆矩阵取代 $P^{-1}$ 中的精确逆矩阵，逆舒尔补码将由压力质量矩阵 $M_p$ 近似（如介绍中提到的由 $\eta^{-1}$ 加权）。正如在 step-22 的结果部分所指出的，我们可以通过应用一个预处理程序来取代 $A$ 的精确逆，在这种情况下，如介绍中所解释的那样，在一个矢量拉普拉斯矩阵上。这确实增加了（外部）GMRES的迭代次数，但仍然比精确的逆运算便宜得多，因为 <em> 的每个 </em> 外部求解器步骤（使用AMG预处理程序）需要20到35次CG迭代。

// 考虑到上述解释，我们定义了一个具有 <code>vmult</code> 功能的预处理类，这就是我们在程序代码中进一步与通常的求解器函数交互所需要的。

// 首先是声明。这与 step-20 中Schur补码的定义相似，不同的是我们在构造函数中需要更多的预处理程序，而且我们在这里使用的矩阵是建立在Trilinos之上的。

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
        stokes_matrix; 
      const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix, 
                                             PreconditionerTypeMp>> 
                                 m_inverse; 
      const PreconditionerTypeA &a_preconditioner; 

      mutable TrilinosWrappers::MPI::Vector tmp; 
    }; 

// 当使用 TrilinosWrappers::MPI::Vector 或 TrilinosWrappers::MPI::BlockVector, 时，Vector被使用IndexSet初始化。IndexSet不仅用于调整 TrilinosWrappers::MPI::Vector 的大小，而且还将 TrilinosWrappers::MPI::Vector 中的一个索引与一个自由度联系起来（更详细的解释见 step-40 ）。函数complete_index_set()创建了一个IndexSet，每个有效的索引都是这个集合的一部分。请注意，这个程序只能按顺序运行，如果并行使用，将抛出一个异常。

    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>:: 
      BlockSchurPreconditioner( 
        const TrilinosWrappers::BlockSparseMatrix &S, 
        const InverseMatrix<TrilinosWrappers::SparseMatrix, 
                            PreconditionerTypeMp> &Mpinv, 
        const PreconditionerTypeA &                Apreconditioner) 
      : stokes_matrix(&S) 
      , m_inverse(&Mpinv) 
      , a_preconditioner(Apreconditioner) 
      , tmp(complete_index_set(stokes_matrix->block(1, 1).m())) 
    {} 

// 接下来是 <code>vmult</code> 函数。我们以三个连续的步骤实现上述 $P^{-1}$ 的动作。 在公式中，我们要计算 $Y=P^{-1}X$ ，其中 $X,Y$ 都是有两个块成分的向量。

// 第一步用矩阵 $A$ 的预处理乘以矢量的速度部分，即计算 $Y_0={\tilde A}^{-1}X_0$  。 然后将得到的速度矢量乘以 $B$ 并减去压力，即我们要计算 $X_1-BY_0$  。这第二步只作用于压力向量，由我们矩阵类的残差函数完成，只是符号不对。因此，我们改变临时压力向量中的符号，最后乘以反压力质量矩阵，得到最终的压力向量，完成我们对斯托克斯预处理的工作。

    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    void 
    BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult( 
      TrilinosWrappers::MPI::BlockVector &      dst, 
      const TrilinosWrappers::MPI::BlockVector &src) const 
    { 
      a_preconditioner.vmult(dst.block(0), src.block(0)); 
      stokes_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1)); 
      tmp *= -1; 
      m_inverse->vmult(dst.block(1), tmp); 
    } 
  } // namespace LinearSolvers 

//  @sect3{The <code>BoussinesqFlowProblem</code> class template}  

// 定义了解决随时间变化的Boussinesq问题的顶层逻辑的类的定义主要是基于 step-22 的教程程序。主要的区别在于，现在我们还必须求解温度方程，这迫使我们为温度变量准备第二个DoFHandler对象，以及当前和之前的时间步骤的矩阵、右手边和求解向量。正如介绍中提到的，所有的线性代数对象都将使用相应的Trilinos功能的包装器。

// 这个类的成员函数让人想起 step-21 ，在那里我们也使用了一个交错的方案，首先解决流动方程（这里是斯托克斯方程， step-21 是达西流），然后更新平流量（这里是温度，那里是饱和度）。新的函数主要涉及到确定时间步长，以及人工粘性稳定的适当大小。

// 最后三个变量表示在下次调用相应的建立函数时，是否需要重建各种矩阵或预处理程序。这使得我们可以将相应的 <code>if</code> 移到相应的函数中，从而使我们的主 <code>run()</code> 函数保持干净，易于阅读。

  template <int dim> 
  class BoussinesqFlowProblem 
  { 
  public: 
    BoussinesqFlowProblem(); 
    void run(); 

  private: 
    void   setup_dofs(); 
    void   assemble_stokes_preconditioner(); 
    void   build_stokes_preconditioner(); 
    void   assemble_stokes_system(); 
    void   assemble_temperature_system(const double maximal_velocity); 
    void   assemble_temperature_matrix(); 
    double get_maximal_velocity() const; 
    std::pair<double, double> get_extrapolated_temperature_range() const; 
    void                      solve(); 
    void                      output_results() const; 
    void                      refine_mesh(const unsigned int max_grid_level); 

    double compute_viscosity( 
      const std::vector<double> &        old_temperature, 
      const std::vector<double> &        old_old_temperature, 
      const std::vector<Tensor<1, dim>> &old_temperature_grads, 
      const std::vector<Tensor<1, dim>> &old_old_temperature_grads, 
      const std::vector<double> &        old_temperature_laplacians, 
      const std::vector<double> &        old_old_temperature_laplacians, 
      const std::vector<Tensor<1, dim>> &old_velocity_values, 
      const std::vector<Tensor<1, dim>> &old_old_velocity_values, 
      const std::vector<double> &        gamma_values, 
      const double                       global_u_infty, 
      const double                       global_T_variation, 
      const double                       cell_diameter) const; 

    Triangulation<dim> triangulation; 
    double             global_Omega_diameter; 

    const unsigned int        stokes_degree; 
    FESystem<dim>             stokes_fe; 
    DoFHandler<dim>           stokes_dof_handler; 
    AffineConstraints<double> stokes_constraints; 

    std::vector<IndexSet>               stokes_partitioning; 
    TrilinosWrappers::BlockSparseMatrix stokes_matrix; 
    TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix; 

    TrilinosWrappers::MPI::BlockVector stokes_solution; 
    TrilinosWrappers::MPI::BlockVector old_stokes_solution; 
    TrilinosWrappers::MPI::BlockVector stokes_rhs; 

    const unsigned int        temperature_degree; 
    FE_Q<dim>                 temperature_fe; 
    DoFHandler<dim>           temperature_dof_handler; 
    AffineConstraints<double> temperature_constraints; 

    TrilinosWrappers::SparseMatrix temperature_mass_matrix; 
    TrilinosWrappers::SparseMatrix temperature_stiffness_matrix; 
    TrilinosWrappers::SparseMatrix temperature_matrix; 

    TrilinosWrappers::MPI::Vector temperature_solution; 
    TrilinosWrappers::MPI::Vector old_temperature_solution; 
    TrilinosWrappers::MPI::Vector old_old_temperature_solution; 
    TrilinosWrappers::MPI::Vector temperature_rhs; 

    double       time_step; 
    double       old_time_step; 
    unsigned int timestep_number; 

    std::shared_ptr<TrilinosWrappers::PreconditionAMG> Amg_preconditioner; 
    std::shared_ptr<TrilinosWrappers::PreconditionIC>  Mp_preconditioner; 

    bool rebuild_stokes_matrix; 
    bool rebuild_temperature_matrices; 
    bool rebuild_stokes_preconditioner; 
  }; 
// @sect3{BoussinesqFlowProblem class implementation}  
// @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}  

// 这个类的构造函数是对  step-22  中的构造函数的扩展。我们需要添加涉及温度的各种变量。正如介绍中所讨论的，我们将再次使用 $Q_2\times Q_1$ （Taylor-Hood）元素来表示斯托克斯部分，并使用 $Q_2$ 元素表示温度。然而，通过使用存储斯托克斯和温度有限元的多项式程度的变量，可以很容易地持续修改这些元素的程度以及下游使用的所有正交公式。此外，我们还初始化了时间步长以及矩阵组合和预处理的选项。

  template <int dim> 
  BoussinesqFlowProblem<dim>::BoussinesqFlowProblem() 
    : triangulation(Triangulation<dim>::maximum_smoothing) 
    , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()) 
    , stokes_degree(1) 
    , stokes_fe(FE_Q<dim>(stokes_degree + 1), dim, FE_Q<dim>(stokes_degree), 1) 
    , stokes_dof_handler(triangulation) 
    , 

    temperature_degree(2) 
    , temperature_fe(temperature_degree) 
    , temperature_dof_handler(triangulation) 
    , 

    time_step(0) 
    , old_time_step(0) 
    , timestep_number(0) 
    , rebuild_stokes_matrix(true) 
    , rebuild_temperature_matrices(true) 
    , rebuild_stokes_preconditioner(true) 
  {} 

//  @sect4{BoussinesqFlowProblem::get_maximal_velocity}  

// 开始这个类的真正功能是一个辅助函数，确定域内（事实上是正交点）的最大（ $L_\infty$  ）速度。它是如何工作的，对所有已经达到本教程这一点的人来说应该是比较明显的。请注意，由于我们只对速度感兴趣，我们不使用 <code>stokes_fe_values.get_function_values</code> 来获取整个斯托克斯解的值（速度和压力），而是使用 <code>stokes_fe_values[velocities].get_function_values</code> 来提取速度部分。这样做的额外好处是，我们得到的是张量<1,dim>，而不是向量<double>中的一些分量，这样我们就可以马上用 <code>norm()</code> 函数来处理它，得到速度的大小。

// 唯一值得思考的一点是如何选择我们在这里使用的正交点。由于这个函数的目标是通过查看每个单元格上的正交点来寻找域内的最大速度。所以我们应该问，我们应该如何最好地选择每个单元上的这些正交点。为此，回顾一下，如果我们有一个单一的 $Q_1$ 场（而不是高阶的矢量值场），那么最大值将在网格的一个顶点达到。换句话说，我们应该使用QTrapezoid类，它的正交点只在单元的顶点。

// 对于高阶形状函数，情况更为复杂：最大值和最小值可能在形状函数的支持点之间达到（对于通常的 $Q_p$ 元素，支持点是等距的Lagrange插值点）；此外，由于我们正在寻找一个矢量值的最大幅值，我们更不能肯定地说潜在的最大点集合在哪里。然而，从直觉上讲，即使不能证明，拉格朗日插值点似乎也是比高斯点更好的选择。

// 现在有不同的方法来产生一个正交公式，其正交点等于有限元的插值点。一种选择是使用 FiniteElement::get_unit_support_points() 函数，将输出减少到一组唯一的点以避免重复的函数评估，并使用这些点创建一个正交对象。另一个选择，这里选择的是使用QTrapezoid类，并将其与QIterated类相结合，该类在每个坐标方向的若干子单元上重复QTrapezoid公式。为了覆盖所有的支持点，我们需要对其进行 <code>stokes_degree+1</code> 次迭代，因为这是使用中的斯托克斯元素的多项式程度。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::get_maximal_velocity() const 
  { 
    const QIterated<dim> quadrature_formula(QTrapezoid<1>(), stokes_degree + 1); 
    const unsigned int   n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(stokes_fe, quadrature_formula, update_values); 
    std::vector<Tensor<1, dim>> velocity_values(n_q_points); 
    double                      max_velocity = 0; 

    const FEValuesExtractors::Vector velocities(0); 

    for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        fe_values[velocities].get_function_values(stokes_solution, 
                                                  velocity_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          max_velocity = std::max(max_velocity, velocity_values[q].norm()); 
      } 

    return max_velocity; 
  } 

//  @sect4{BoussinesqFlowProblem::get_extrapolated_temperature_range}  

// 接下来是一个函数，确定从前两个时间步长推算到当前步长时， $\Omega$ 内正交点的最低和最高温度。我们在计算人工粘性参数 $\nu$ 时需要这个信息，正如在介绍中所讨论的那样。

// 外推温度的公式是  $\left(1+\frac{k_n}{k_{n-1}} \right)T^{n-1} + \frac{k_n}{k_{n-1}} T^{n-2}$  。计算的方法是在所有正交点上循环，如果当前值比前一个值大/小，则更新最大和最小值。在对所有正交点进行循环之前，我们将存储最大和最小值的变量初始化为可表示为双数的最小和最大的数字。这样我们就知道它比最小/最大值大/小，并且所有正交点的循环最终会用正确的值来更新初始值。

// 这里唯一值得一提的复杂情况是，在第一个时间步骤中， $T^{k-2}$ 当然还不能使用。在这种情况下，我们只能使用 $T^{k-1}$ ，这是我们从初始温度得到的。作为正交点，我们使用与前一个函数相同的选择，但不同的是，现在重复的数量由温度场的多项式程度决定。

  template <int dim> 
  std::pair<double, double> 
  BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range() const 
  { 
    const QIterated<dim> quadrature_formula(QTrapezoid<1>(), 
                                            temperature_degree); 
    const unsigned int   n_q_points = quadrature_formula.size(); 

    FEValues<dim> fe_values(temperature_fe, quadrature_formula, update_values); 
    std::vector<double> old_temperature_values(n_q_points); 
    std::vector<double> old_old_temperature_values(n_q_points); 

    if (timestep_number != 0) 
      { 
        double min_temperature = std::numeric_limits<double>::max(), 
               max_temperature = -std::numeric_limits<double>::max(); 

        for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
          { 
            fe_values.reinit(cell); 
            fe_values.get_function_values(old_temperature_solution, 
                                          old_temperature_values); 
            fe_values.get_function_values(old_old_temperature_solution, 
                                          old_old_temperature_values); 

            for (unsigned int q = 0; q < n_q_points; ++q) 
              { 
                const double temperature = 
                  (1. + time_step / old_time_step) * old_temperature_values[q] - 
                  time_step / old_time_step * old_old_temperature_values[q]; 

                min_temperature = std::min(min_temperature, temperature); 
                max_temperature = std::max(max_temperature, temperature); 
              } 
          } 

        return std::make_pair(min_temperature, max_temperature); 
      } 
    else 
      { 
        double min_temperature = std::numeric_limits<double>::max(), 
               max_temperature = -std::numeric_limits<double>::max(); 

        for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
          { 
            fe_values.reinit(cell); 
            fe_values.get_function_values(old_temperature_solution, 
                                          old_temperature_values); 

            for (unsigned int q = 0; q < n_q_points; ++q) 
              { 
                const double temperature = old_temperature_values[q]; 

                min_temperature = std::min(min_temperature, temperature); 
                max_temperature = std::max(max_temperature, temperature); 
              } 
          } 

        return std::make_pair(min_temperature, max_temperature); 
      } 
  } 

//  @sect4{BoussinesqFlowProblem::compute_viscosity}  

// 最后一个工具函数计算单元 $\nu|_K$ 上的人工粘度参数 $K$ ，作为外推温度、其梯度和Hessian（二阶导数）、速度、当前单元正交点上的所有右手 $\gamma$ 和其他各种参数的函数，在介绍中已详细说明。

// 这里有一些值得一提的通用常数。首先，我们需要固定 $\beta$ ；我们选择 $\beta=0.017\cdot dim$ ，这个选择在本教程程序的结果部分有详细讨论。其次是指数 $\alpha$ ； $\alpha=1$ 对于目前的程序似乎很好用，尽管选择 $\alpha = 2$ 可能会有一些额外的好处。最后，有一件事需要特别说明。在第一个时间步骤中，速度等于零， $\nu|_K$ 的公式没有定义。在这种情况下，我们返回 $\nu|_K=5\cdot 10^3 \cdot h_K$ ，这个选择无疑更多的是出于启发式的考虑（不过，它与第二个时间步骤中大多数单元的返回值处于同一数量级）。

// 根据介绍中讨论的材料，该函数的其余部分应该是显而易见的。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::compute_viscosity( 
    const std::vector<double> &        old_temperature, 
    const std::vector<double> &        old_old_temperature, 
    const std::vector<Tensor<1, dim>> &old_temperature_grads, 
    const std::vector<Tensor<1, dim>> &old_old_temperature_grads, 
    const std::vector<double> &        old_temperature_laplacians, 
    const std::vector<double> &        old_old_temperature_laplacians, 
    const std::vector<Tensor<1, dim>> &old_velocity_values, 
    const std::vector<Tensor<1, dim>> &old_old_velocity_values, 
    const std::vector<double> &        gamma_values, 
    const double                       global_u_infty, 
    const double                       global_T_variation, 
    const double                       cell_diameter) const 
  { 
    constexpr double beta  = 0.017 * dim; 
    constexpr double alpha = 1.0; 

    if (global_u_infty == 0) 
      return 5e-3 * cell_diameter; 

    const unsigned int n_q_points = old_temperature.size(); 

    double max_residual = 0; 
    double max_velocity = 0; 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const Tensor<1, dim> u = 
          (old_velocity_values[q] + old_old_velocity_values[q]) / 2; 

        const double dT_dt = 
          (old_temperature[q] - old_old_temperature[q]) / old_time_step; 
        const double u_grad_T = 
          u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2; 

        const double kappa_Delta_T = 
          EquationData::kappa * 
          (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) / 
          2; 

        const double residual = 
          std::abs((dT_dt + u_grad_T - kappa_Delta_T - gamma_values[q]) * 
                   std::pow((old_temperature[q] + old_old_temperature[q]) / 2, 
                            alpha - 1.)); 

        max_residual = std::max(residual, max_residual); 
        max_velocity = std::max(std::sqrt(u * u), max_velocity); 
      } 

    const double c_R            = std::pow(2., (4. - 2 * alpha) / dim); 
    const double global_scaling = c_R * global_u_infty * global_T_variation * 
                                  std::pow(global_Omega_diameter, alpha - 2.); 

    return ( 
      beta * max_velocity * 
      std::min(cell_diameter, 
               std::pow(cell_diameter, alpha) * max_residual / global_scaling)); 
  } 

//  @sect4{BoussinesqFlowProblem::setup_dofs}  

// 这是一个函数，用于设置我们这里的DoFHandler对象（一个用于斯托克斯部分，一个用于温度部分），以及将本程序中线性代数所需的各种对象设置为合适的尺寸。它的基本操作与我们在  step-22  中的操作类似。

// 该函数的主体首先列举了斯托克斯和温度系统的所有自由度。对于斯托克斯部分，自由度被排序，以确保速度优先于压力自由度，这样我们就可以将斯托克斯矩阵划分为一个 $2\times 2$ 矩阵。作为与 step-22 的区别，我们不进行任何额外的DoF重新编号。在那个程序中，它得到了回报，因为我们的求解器严重依赖ILU，而我们在这里使用AMG，它对DoF编号不敏感。用于压力质量矩阵反演的IC预处理程序当然会利用类似Cuthill-McKee的重新编号，但是与速度部分相比，其成本很低，所以额外的工作并没有得到回报。

// 然后，我们继续生成悬挂的节点约束，这些约束来自两个DoFHandler对象的自适应网格细化。对于速度，我们通过向已经存储了悬挂节点约束矩阵的对象添加约束来施加无流边界条件 $\mathbf{u}\cdot \mathbf{n}=0$ 。函数中的第二个参数描述了总dof向量中的第一个速度分量，这里是零。变量 <code>no_normal_flux_boundaries</code> 表示要设置无通量边界条件的边界指标；这里是边界指标0。

// 做完这些后，我们计算各块中的自由度数量。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::setup_dofs() 
  { 
    std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0); 
    stokes_sub_blocks[dim] = 1; 

    { 
      stokes_dof_handler.distribute_dofs(stokes_fe); 
      DoFRenumbering::component_wise(stokes_dof_handler, stokes_sub_blocks); 

      stokes_constraints.clear(); 
      DoFTools::make_hanging_node_constraints(stokes_dof_handler, 
                                              stokes_constraints); 
      std::set<types::boundary_id> no_normal_flux_boundaries; 
      no_normal_flux_boundaries.insert(0); 
      VectorTools::compute_no_normal_flux_constraints(stokes_dof_handler, 
                                                      0, 
                                                      no_normal_flux_boundaries, 
                                                      stokes_constraints); 
      stokes_constraints.close(); 
    } 
    { 
      temperature_dof_handler.distribute_dofs(temperature_fe); 

      temperature_constraints.clear(); 
      DoFTools::make_hanging_node_constraints(temperature_dof_handler, 
                                              temperature_constraints); 
      temperature_constraints.close(); 
    } 

    const std::vector<types::global_dof_index> stokes_dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(stokes_dof_handler, stokes_sub_blocks); 

    const unsigned int n_u = stokes_dofs_per_block[0], 
                       n_p = stokes_dofs_per_block[1], 
                       n_T = temperature_dof_handler.n_dofs(); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << " (on " << triangulation.n_levels() << " levels)" << std::endl 
              << "Number of degrees of freedom: " << n_u + n_p + n_T << " (" 
              << n_u << '+' << n_p << '+' << n_T << ')' << std::endl 
              << std::endl; 

// 下一步是创建斯托克斯和温度系统矩阵的稀疏模式，以及建立斯托克斯预处理矩阵的预处理。如同在 step-22 中一样，我们选择使用DynamicSparsityPattern的封锁版本来创建模式。

// 因此，我们首先释放存储在矩阵中的内存，然后建立一个BlockDynamicSparsityPattern类型的对象，该对象由 $2\times 2$ 块（用于斯托克斯系统矩阵和预处理器）或DynamicSparsityPattern（用于温度部分）组成。然后我们用非零模式填充这些对象，考虑到对于斯托克斯系统矩阵，在压力-压力块中没有条目（但所有速度矢量分量相互耦合并与压力耦合）。同样，在斯托克斯预处理矩阵中，只有对角线块是非零的，因为我们使用了介绍中讨论的矢量拉普拉斯。这个算子只把拉普拉斯的每个矢量分量与它自己联系起来，而不是与其他矢量分量联系起来。然而，应用无流量边界条件产生的约束条件将在边界处再次耦合向量分量）。

// 在生成稀疏模式时，我们直接应用悬挂节点和无流边界条件的约束。这种方法在 step-27 中已经使用过了，但与早期教程中的方法不同，在早期教程中我们先建立原始的稀疏模式，然后才加入约束条件产生的条目。这样做的原因是，在以后的装配过程中，我们要在将本地道夫转移到全局道夫时立即分配约束。因此，在受限自由度的位置不会有数据写入，所以我们可以通过将最后一个布尔标志设置为 <code>false</code> ，让 DoFTools::make_sparsity_pattern 函数省略这些条目。一旦稀疏性模式准备好了，我们就可以用它来初始化特里诺斯矩阵。由于Trilinos矩阵在内部存储了稀疏模式，所以在初始化矩阵之后，没有必要再保留稀疏模式。

    stokes_partitioning.resize(2); 
    stokes_partitioning[0] = complete_index_set(n_u); 
    stokes_partitioning[1] = complete_index_set(n_p); 
    { 
      stokes_matrix.clear(); 

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
        stokes_dof_handler, coupling, dsp, stokes_constraints, false); 

      stokes_matrix.reinit(dsp); 
    } 

    { 
      Amg_preconditioner.reset(); 
      Mp_preconditioner.reset(); 
      stokes_preconditioner_matrix.clear(); 

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
        stokes_dof_handler, coupling, dsp, stokes_constraints, false); 

      stokes_preconditioner_matrix.reinit(dsp); 
    } 

// 温度矩阵（或者说是矩阵，因为我们提供了一个温度质量矩阵和一个温度刚度矩阵，它们将在时间离散化中被加在一起）的创建与斯托克斯矩阵的生成相同；只是在这里要简单得多，因为我们不需要照顾任何块或组件之间的耦合。注意我们是如何初始化三个温度矩阵的。我们只使用稀疏模式对第一个矩阵进行再初始化，而对其余两个再初始化则使用先前生成的矩阵。这样做的原因是，从一个已经生成的矩阵进行重新初始化，可以让Trilinos重新使用稀疏模式，而不是为每个副本生成一个新的模式。这样可以节省一些时间和内存。

    { 
      temperature_mass_matrix.clear(); 
      temperature_stiffness_matrix.clear(); 
      temperature_matrix.clear(); 

      DynamicSparsityPattern dsp(n_T, n_T); 
      DoFTools::make_sparsity_pattern(temperature_dof_handler, 
                                      dsp, 
                                      temperature_constraints, 
                                      false); 

      temperature_matrix.reinit(dsp); 
      temperature_mass_matrix.reinit(temperature_matrix); 
      temperature_stiffness_matrix.reinit(temperature_matrix); 
    } 

// 最后，我们将斯托克斯解的向量 $\mathbf u^{n-1}$ 和 $\mathbf u^{n-2}$ ，以及温度 $T^{n}$ 、 $T^{n-1}$ 和 $T^{n-2}$ （时间步进所需）和所有系统的右手边设置为正确的大小和块结构。

    IndexSet temperature_partitioning = complete_index_set(n_T); 
    stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD); 
    old_stokes_solution.reinit(stokes_partitioning, MPI_COMM_WORLD); 
    stokes_rhs.reinit(stokes_partitioning, MPI_COMM_WORLD); 

    temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD); 
    old_temperature_solution.reinit(temperature_partitioning, MPI_COMM_WORLD); 
    old_old_temperature_solution.reinit(temperature_partitioning, 
                                        MPI_COMM_WORLD); 

    temperature_rhs.reinit(temperature_partitioning, MPI_COMM_WORLD); 
  } 

//  @sect4{BoussinesqFlowProblem::assemble_stokes_preconditioner}  

// 这个函数组装了我们用于预处理斯托克斯系统的矩阵。我们需要的是速度分量上的矢量拉普拉斯矩阵和压力分量上的质量矩阵，并以 $\eta^{-1}$ 加权。我们首先生成一个适当阶数的正交对象，即FEValues对象，它可以给出正交点的值和梯度（连同正交权重）。接下来我们为单元格矩阵和局部与全局DoF之间的关系创建数据结构。向量 <code>grad_phi_u</code> and <code>phi_p</code> 将保存基函数的值，以便更快地建立局部矩阵，正如在 step-22 中已经完成的那样。在我们开始对所有活动单元进行循环之前，我们必须指定哪些成分是压力，哪些是速度。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner() 
  { 
    stokes_preconditioner_matrix = 0; 

    const QGauss<dim> quadrature_formula(stokes_degree + 2); 
    FEValues<dim>     stokes_fe_values(stokes_fe, 
                                   quadrature_formula, 
                                   update_JxW_values | update_values | 
                                     update_gradients); 

    const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
    std::vector<double>         phi_p(dofs_per_cell); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

    for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
      { 
        stokes_fe_values.reinit(cell); 
        local_matrix = 0; 

// 本地矩阵的创建相当简单。只有一个拉普拉斯项（关于速度）和一个由 $\eta^{-1}$ 加权的质量矩阵需要生成，所以本地矩阵的创建在两行中完成。一旦本地矩阵准备好了（在每个正交点上循环查看本地矩阵的行和列），我们就可以得到本地的DoF指数，并将本地信息写入全局矩阵中。我们像在 step-27 中那样做，也就是说，我们直接应用本地悬挂节点的约束。这样做，我们就不必事后再做，而且我们也不会在消除约束时将矩阵的条目写成实际上将再次设置为零。

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                grad_phi_u[k] = stokes_fe_values[velocities].gradient(k, q); 
                phi_p[k]      = stokes_fe_values[pressure].value(k, q); 
              } 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                local_matrix(i, j) += 
                  (EquationData::eta * 
                     scalar_product(grad_phi_u[i], grad_phi_u[j]) + 
                   (1. / EquationData::eta) * phi_p[i] * phi_p[j]) * 
                  stokes_fe_values.JxW(q); 
          } 

        cell->get_dof_indices(local_dof_indices); 
        stokes_constraints.distribute_local_to_global( 
          local_matrix, local_dof_indices, stokes_preconditioner_matrix); 
      } 
  } 

//  @sect4{BoussinesqFlowProblem::build_stokes_preconditioner}  

// 这个函数生成将用于Schur互补块预处理的内部预处理。由于只有当矩阵发生变化时才需要重新生成预处理程序，因此在矩阵没有变化的情况下，该函数不需要做任何事情（即标志 <code>rebuild_stokes_preconditioner</code> 的值为 <code>false</code> ）。否则，它的第一个任务是调用 <code>assemble_stokes_preconditioner</code> 来生成预处理矩阵。

// 接下来，我们为速度-速度矩阵  $A$  设置预处理程序。正如介绍中所解释的，我们将使用基于矢量拉普拉斯矩阵 $\hat{A}$ 的AMG预处理器（它在频谱上与斯托克斯矩阵 $A$ 接近）。通常， TrilinosWrappers::PreconditionAMG 类可以被看作是一个好的黑箱预处理程序，不需要任何特殊的知识。然而，在这种情况下，我们必须小心：因为我们为一个矢量问题建立了一个AMG，我们必须告诉预处理程序设置哪个道夫属于哪个矢量成分。我们使用 DoFTools::extract_constant_modes, 函数来做这件事，该函数生成一组 <code>dim</code> 向量，其中每个向量在向量问题的相应分量中为1，在其他地方为0。因此，这些是每个分量上的常数模式，这解释了变量的名称。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::build_stokes_preconditioner() 
  { 
    if (rebuild_stokes_preconditioner == false) 
      return; 

    std::cout << "   Rebuilding Stokes preconditioner..." << std::flush; 

    assemble_stokes_preconditioner(); 

    Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>(); 

    std::vector<std::vector<bool>> constant_modes; 
    FEValuesExtractors::Vector     velocity_components(0); 
    DoFTools::extract_constant_modes(stokes_dof_handler, 
                                     stokes_fe.component_mask( 
                                       velocity_components), 
                                     constant_modes); 
    TrilinosWrappers::PreconditionAMG::AdditionalData amg_data; 
    amg_data.constant_modes = constant_modes; 

// 接下来，我们再设置一些AMG预处理程序的选项。特别是，我们需要告诉AMG设置，我们对速度矩阵使用二次基函数（这意味着矩阵中有更多的非零元素，因此需要在内部选择一种更稳健的算法）。此外，我们希望能够控制粗化结构的建立方式。Trilinos平滑聚合AMG的方法是寻找哪些矩阵条目与对角线条目大小相似，以便代数式地建立一个粗网格结构。通过将参数 <code>aggregation_threshold</code> 设置为0.02，我们指定所有尺寸超过该行中一些对角线枢轴的百分之二的条目应该形成一个粗网格点。这个参数是比较特别的，对它进行一些微调会影响预处理程序的性能。根据经验，较大的 <code>aggregation_threshold</code> 值会减少迭代次数，但增加每次迭代的成本。看一下Trilinos的文档会提供更多关于这些参数的信息。有了这个数据集，我们就用我们想要的矩阵来初始化预处理程序。

// 最后，我们也初始化预处理程序以反转压力质量矩阵。这个矩阵是对称的，表现良好，所以我们可以选择一个简单的预处理程序。我们坚持使用不完全Cholesky（IC）因子化预处理器，它是为对称矩阵设计的。我们也可以选择SSOR预处理器，其松弛系数约为1.2，但IC对我们的例子来说更便宜。我们把预处理程序包成一个 <code>std::shared_ptr</code> 指针，这使得下次重新创建预处理程序更加容易，因为我们不必关心破坏以前使用的对象。

    amg_data.elliptic              = true; 
    amg_data.higher_order_elements = true; 
    amg_data.smoother_sweeps       = 2; 
    amg_data.aggregation_threshold = 0.02; 
    Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0, 0), 
                                   amg_data); 

    Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
    Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1, 1)); 

    std::cout << std::endl; 

    rebuild_stokes_preconditioner = false; 
  } 

//  @sect4{BoussinesqFlowProblem::assemble_stokes_system}  

// 我们用于推进耦合的斯托克斯-温度系统的时滞方案迫使我们将装配（以及线性系统的解）分成两步。第一步是创建斯托克斯系统的矩阵和右手边，第二步是创建温度道夫的矩阵和右手边，这取决于速度的线性系统的结果。

// 该函数在每个时间步长的开始时被调用。在第一个时间步骤中，或者如果网格已经改变，由 <code>rebuild_stokes_matrix</code> 表示，我们需要组装斯托克斯矩阵；另一方面，如果网格没有改变，矩阵已经有了，这就没有必要了，我们需要做的就是组装右手边的向量，它在每个时间步骤中都会改变。

// 关于实现的技术细节，与  step-22  相比没有太大变化。我们重置矩阵和向量，在单元格上创建正交公式，然后创建相应的FEValues对象。对于更新标志，我们只在完全装配的情况下需要基函数导数，因为右手边不需要它们；像往常一样，根据当前需要选择最小的标志集，使程序中进一步调用  FEValues::reinit  的效率更高。

// 有一件事需要评论&ndash；因为我们有一个单独的有限元和DoFHandler来处理温度问题，所以我们需要生成第二个FEValues对象来正确评估温度解决方案。要实现这一点并不复杂：只需使用温度结构，并为我们需要用于评估温度解决方案的基函数值设置一个更新标志。这里需要记住的唯一重要部分是，两个FEValues对象使用相同的正交公式，以确保我们在循环计算两个对象的正交点时得到匹配的信息。

// 声明的过程中，有一些关于数组大小的快捷方式，本地矩阵和右手的创建，以及与全局系统相比，本地道夫的索引的向量。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_stokes_system() 
  { 
    std::cout << "   Assembling..." << std::flush; 

    if (rebuild_stokes_matrix == true) 
      stokes_matrix = 0; 

    stokes_rhs = 0; 

    const QGauss<dim> quadrature_formula(stokes_degree + 2); 
    FEValues<dim>     stokes_fe_values( 
      stokes_fe, 
      quadrature_formula, 
      update_values | update_quadrature_points | update_JxW_values | 
        (rebuild_stokes_matrix == true ? update_gradients : UpdateFlags(0))); 

    FEValues<dim> temperature_fe_values(temperature_fe, 
                                        quadrature_formula, 
                                        update_values); 

    const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 接下来我们需要一个向量，它将包含前一个时间层的温度解在正交点的值，以组装动量方程右侧的源项。让我们把这个向量称为  <code>old_solution_values</code>  。

// 我们接下来创建的向量集包含了基函数的评估以及它们的梯度和对称梯度，将用于创建矩阵。将这些放到自己的数组中，而不是每次都向FEValues对象索取这些信息，是为了加速装配过程的优化，详情请参见 step-22 。

// 最后两个声明是用来从整个FE系统中提取各个块（速度、压力、温度）的。

    std::vector<double> old_temperature_values(n_q_points); 

    std::vector<Tensor<1, dim>>          phi_u(dofs_per_cell); 
    std::vector<SymmetricTensor<2, dim>> grads_phi_u(dofs_per_cell); 
    std::vector<double>                  div_phi_u(dofs_per_cell); 
    std::vector<double>                  phi_p(dofs_per_cell); 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

// 现在开始对问题中的所有单元格进行循环。我们正在为这个装配例程处理两个不同的DoFHandlers，所以我们必须为使用中的两个对象设置两个不同的单元格迭代器。这可能看起来有点奇怪，因为斯托克斯系统和温度系统都使用相同的网格，但这是保持自由度同步的唯一方法。循环中的第一条语句也是非常熟悉的，按照更新标志的规定对有限元数据进行更新，将局部数组清零，并在正交点处获得旧解的值。然后我们准备在单元格上的正交点上循环。

    auto       cell             = stokes_dof_handler.begin_active(); 
    const auto endc             = stokes_dof_handler.end(); 
    auto       temperature_cell = temperature_dof_handler.begin_active(); 

    for (; cell != endc; ++cell, ++temperature_cell) 
      { 
        stokes_fe_values.reinit(cell); 
        temperature_fe_values.reinit(temperature_cell); 

        local_matrix = 0; 
        local_rhs    = 0; 

        temperature_fe_values.get_function_values(old_temperature_solution, 
                                                  old_temperature_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double old_temperature = old_temperature_values[q]; 

// 接下来我们提取与内积中的条款相关的基础函数的值和梯度。如 step-22 所示，这有助于加速装配。    一旦完成，我们开始在本地矩阵的行和列上进行循环，并将相关的乘积送入矩阵。右手边是由温度驱动的重力方向（在我们的例子中是垂直方向）的强迫项。 请注意，右手边的项总是生成的，而矩阵的贡献只有在 <code>rebuild_matrices</code> 标志要求时才会更新。

            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                phi_u[k] = stokes_fe_values[velocities].value(k, q); 
                if (rebuild_stokes_matrix) 
                  { 
                    grads_phi_u[k] = 
                      stokes_fe_values[velocities].symmetric_gradient(k, q); 
                    div_phi_u[k] = 
                      stokes_fe_values[velocities].divergence(k, q); 
                    phi_p[k] = stokes_fe_values[pressure].value(k, q); 
                  } 
              } 

            if (rebuild_stokes_matrix) 
              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  local_matrix(i, j) += 
                    (EquationData::eta * 2 * (grads_phi_u[i] * grads_phi_u[j]) - 
                     div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
                    stokes_fe_values.JxW(q); 

            const Point<dim> gravity = 
              -((dim == 2) ? (Point<dim>(0, 1)) : (Point<dim>(0, 0, 1))); 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              local_rhs(i) += (-EquationData::density * EquationData::beta * 
                               gravity * phi_u[i] * old_temperature) * 
                              stokes_fe_values.JxW(q); 
          } 

// 循环所有单元的最后一步是将局部贡献输入到全局矩阵和向量结构中，并将其输入到  <code>local_dof_indices</code>  指定的位置。 同样，我们让AffineConstraints类来完成将单元格矩阵元素插入全局矩阵的工作，这已经浓缩了悬挂的节点约束。

        cell->get_dof_indices(local_dof_indices); 

        if (rebuild_stokes_matrix == true) 
          stokes_constraints.distribute_local_to_global(local_matrix, 
                                                        local_rhs, 
                                                        local_dof_indices, 
                                                        stokes_matrix, 
                                                        stokes_rhs); 
        else 
          stokes_constraints.distribute_local_to_global(local_rhs, 
                                                        local_dof_indices, 
                                                        stokes_rhs); 
      } 

    rebuild_stokes_matrix = false; 

    std::cout << std::endl; 
  } 

//  @sect4{BoussinesqFlowProblem::assemble_temperature_matrix}  

// 这个函数组装温度方程中的矩阵。温度矩阵由两部分组成，质量矩阵和时间步长乘以刚度矩阵，由拉普拉斯项乘以扩散量给出。由于该矩阵取决于时间步长（从一个步长到另一个步长），温度矩阵需要在每个时间步长进行更新。我们可以简单地在每个时间步长中重新生成矩阵，但这并不真正有效，因为质量和拉普拉斯矩阵只有在我们改变网格时才会改变。因此，我们通过在这个函数中生成两个单独的矩阵，一个是质量矩阵，一个是刚度（扩散）矩阵，这样做更有效率。一旦我们知道了实际的时间步长，我们将把这个矩阵加上刚度矩阵乘以时间步长的总和。

// 所以这第一步的细节非常简单。为了防止我们需要重建矩阵（即网格发生了变化），我们将数据结构归零，得到一个正交公式和一个FEValues对象，并为基函数创建局部矩阵、局部dof指数和评估结构。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_temperature_matrix() 
  { 
    if (rebuild_temperature_matrices == false) 
      return; 

    temperature_mass_matrix      = 0; 
    temperature_stiffness_matrix = 0; 

    QGauss<dim>   quadrature_formula(temperature_degree + 2); 
    FEValues<dim> temperature_fe_values(temperature_fe, 
                                        quadrature_formula, 
                                        update_values | update_gradients | 
                                          update_JxW_values); 

    const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> local_mass_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> local_stiffness_matrix(dofs_per_cell, dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    std::vector<double>         phi_T(dofs_per_cell); 
    std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell); 

// 现在，让我们开始在三角结构中的所有单元上进行循环。我们需要将局部矩阵清零，更新有限元评估，然后在每个正交点上循环矩阵的行和列，然后我们创建质量矩阵和刚度矩阵（拉普拉斯项乘以扩散  <code>EquationData::kappa</code>  。最后，我们让约束对象将这些值插入全局矩阵中，并直接将约束条件浓缩到矩阵中。

    for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
      { 
        local_mass_matrix      = 0; 
        local_stiffness_matrix = 0; 

        temperature_fe_values.reinit(cell); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                grad_phi_T[k] = temperature_fe_values.shape_grad(k, q); 
                phi_T[k]      = temperature_fe_values.shape_value(k, q); 
              } 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  local_mass_matrix(i, j) += 
                    (phi_T[i] * phi_T[j] * temperature_fe_values.JxW(q)); 
                  local_stiffness_matrix(i, j) += 
                    (EquationData::kappa * grad_phi_T[i] * grad_phi_T[j] * 
                     temperature_fe_values.JxW(q)); 
                } 
          } 

        cell->get_dof_indices(local_dof_indices); 

        temperature_constraints.distribute_local_to_global( 
          local_mass_matrix, local_dof_indices, temperature_mass_matrix); 
        temperature_constraints.distribute_local_to_global( 
          local_stiffness_matrix, 
          local_dof_indices, 
          temperature_stiffness_matrix); 
      } 

    rebuild_temperature_matrices = false; 
  } 

//  @sect4{BoussinesqFlowProblem::assemble_temperature_system}  

// 这个函数对温度矩阵进行第二部分的装配工作，实际添加压力质量和刚度矩阵（时间步长在这里起作用），以及创建依赖于速度的右手边。这个函数中的右侧装配的声明与其他装配例程中使用的声明基本相同，只是这次我们把自己限制在矢量上。我们将计算温度系统的残差，这意味着我们必须评估二阶导数，由更新标志 <code>update_hessians</code> 指定。

// 温度方程通过流体速度与斯托克斯系统相耦合。解决方案的这两部分与不同的DoFHandlers相关联，因此我们需要再次创建第二个FEValues对象来评估正交点的速度。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_temperature_system( 
    const double maximal_velocity) 
  { 
    const bool use_bdf2_scheme = (timestep_number != 0); 

    if (use_bdf2_scheme == true) 
      { 
        temperature_matrix.copy_from(temperature_mass_matrix); 
        temperature_matrix *= 
          (2 * time_step + old_time_step) / (time_step + old_time_step); 
        temperature_matrix.add(time_step, temperature_stiffness_matrix); 
      } 
    else 
      { 
        temperature_matrix.copy_from(temperature_mass_matrix); 
        temperature_matrix.add(time_step, temperature_stiffness_matrix); 
      } 

    temperature_rhs = 0; 

    const QGauss<dim> quadrature_formula(temperature_degree + 2); 
    FEValues<dim>     temperature_fe_values(temperature_fe, 
                                        quadrature_formula, 
                                        update_values | update_gradients | 
                                          update_hessians | 
                                          update_quadrature_points | 
                                          update_JxW_values); 
    FEValues<dim>     stokes_fe_values(stokes_fe, 
                                   quadrature_formula, 
                                   update_values); 

    const unsigned int dofs_per_cell = temperature_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double> local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

// 接下来是向量的声明，用来保存旧的和更早的解决方案的值（分别作为时间级别 $n-1$ 和 $n-2$ 的符号）和当前单元的正交点的梯度。我们还声明了一个对象来保存温度的右侧值（ <code>gamma_values</code> ），并且我们再次使用温度基函数的快捷方式。最终，我们需要找到温度极值和计算域的直径，这将用于稳定参数的定义（我们得到了最大速度作为这个函数的输入）。

    std::vector<Tensor<1, dim>> old_velocity_values(n_q_points); 
    std::vector<Tensor<1, dim>> old_old_velocity_values(n_q_points); 
    std::vector<double>         old_temperature_values(n_q_points); 
    std::vector<double>         old_old_temperature_values(n_q_points); 
    std::vector<Tensor<1, dim>> old_temperature_grads(n_q_points); 
    std::vector<Tensor<1, dim>> old_old_temperature_grads(n_q_points); 
    std::vector<double>         old_temperature_laplacians(n_q_points); 
    std::vector<double>         old_old_temperature_laplacians(n_q_points); 

    EquationData::TemperatureRightHandSide<dim> temperature_right_hand_side; 
    std::vector<double>                         gamma_values(n_q_points); 

    std::vector<double>         phi_T(dofs_per_cell); 
    std::vector<Tensor<1, dim>> grad_phi_T(dofs_per_cell); 

    const std::pair<double, double> global_T_range = 
      get_extrapolated_temperature_range(); 

    const FEValuesExtractors::Vector velocities(0); 

// 现在，让我们开始在三角结构中的所有单元格上进行循环。同样，我们需要两个单元格迭代器，平行走过两个参与的DoFHandler对象的单元格，用于斯托克斯和温度部分。在这个循环中，我们首先将局部rhs设置为零，然后在正交点上获得旧的解函数的值和导数，因为它们将被用于稳定参数的定义和作为方程中的系数，分别需要。请注意，由于温度有自己的DoFHandler和FEValues对象，我们在正交点得到整个解（反正只有标量温度场），而对于斯托克斯部分，我们仅限于通过使用 <code>stokes_fe_values[velocities].get_function_values</code> 提取速度部分（而忽略压力部分）。

    auto       cell        = temperature_dof_handler.begin_active(); 
    const auto endc        = temperature_dof_handler.end(); 
    auto       stokes_cell = stokes_dof_handler.begin_active(); 

    for (; cell != endc; ++cell, ++stokes_cell) 
      { 
        local_rhs = 0; 

        temperature_fe_values.reinit(cell); 
        stokes_fe_values.reinit(stokes_cell); 

        temperature_fe_values.get_function_values(old_temperature_solution, 
                                                  old_temperature_values); 
        temperature_fe_values.get_function_values(old_old_temperature_solution, 
                                                  old_old_temperature_values); 

        temperature_fe_values.get_function_gradients(old_temperature_solution, 
                                                     old_temperature_grads); 
        temperature_fe_values.get_function_gradients( 
          old_old_temperature_solution, old_old_temperature_grads); 

        temperature_fe_values.get_function_laplacians( 
          old_temperature_solution, old_temperature_laplacians); 
        temperature_fe_values.get_function_laplacians( 
          old_old_temperature_solution, old_old_temperature_laplacians); 

        temperature_right_hand_side.value_list( 
          temperature_fe_values.get_quadrature_points(), gamma_values); 

        stokes_fe_values[velocities].get_function_values(stokes_solution, 
                                                         old_velocity_values); 
        stokes_fe_values[velocities].get_function_values( 
          old_stokes_solution, old_old_velocity_values); 

// 接下来，我们根据介绍中的讨论，使用专用函数计算用于稳定的人工粘性。有了这个，我们就可以进入正交点和局部rhs矢量分量的循环了。这里的术语相当冗长，但其定义遵循本方案介绍中开发的时间-离散系统。BDF-2方案比用于第一时间步的后向欧拉方案多需要一个旧时间步的术语（并且涉及更复杂的因素）。当所有这些都完成后，我们将局部向量分配到全局向量中（包括悬挂节点约束）。

        const double nu = 
          compute_viscosity(old_temperature_values, 
                            old_old_temperature_values, 
                            old_temperature_grads, 
                            old_old_temperature_grads, 
                            old_temperature_laplacians, 
                            old_old_temperature_laplacians, 
                            old_velocity_values, 
                            old_old_velocity_values, 
                            gamma_values, 
                            maximal_velocity, 
                            global_T_range.second - global_T_range.first, 
                            cell->diameter()); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                grad_phi_T[k] = temperature_fe_values.shape_grad(k, q); 
                phi_T[k]      = temperature_fe_values.shape_value(k, q); 
              } 

            const double T_term_for_rhs = 
              (use_bdf2_scheme ? 
                 (old_temperature_values[q] * (1 + time_step / old_time_step) - 
                  old_old_temperature_values[q] * (time_step * time_step) / 
                    (old_time_step * (time_step + old_time_step))) : 
                 old_temperature_values[q]); 

            const Tensor<1, dim> ext_grad_T = 
              (use_bdf2_scheme ? 
                 (old_temperature_grads[q] * (1 + time_step / old_time_step) - 
                  old_old_temperature_grads[q] * time_step / old_time_step) : 
                 old_temperature_grads[q]); 

            const Tensor<1, dim> extrapolated_u = 
              (use_bdf2_scheme ? 
                 (old_velocity_values[q] * (1 + time_step / old_time_step) - 
                  old_old_velocity_values[q] * time_step / old_time_step) : 
                 old_velocity_values[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              local_rhs(i) += 
                (T_term_for_rhs * phi_T[i] - 
                 time_step * extrapolated_u * ext_grad_T * phi_T[i] - 
                 time_step * nu * ext_grad_T * grad_phi_T[i] + 
                 time_step * gamma_values[q] * phi_T[i]) * 
                temperature_fe_values.JxW(q); 
          } 

        cell->get_dof_indices(local_dof_indices); 
        temperature_constraints.distribute_local_to_global(local_rhs, 
                                                           local_dof_indices, 
                                                           temperature_rhs); 
      } 
  } 

//  @sect4{BoussinesqFlowProblem::solve}  

// 这个函数可以解决线性方程组的问题。在介绍之后，我们从斯托克斯系统开始，在这里我们需要生成我们的块状舒尔预处理器。由于所有相关的动作都在类 <code>BlockSchurPreconditioner</code> 中实现，我们所要做的就是适当地初始化这个类。我们需要传递的是一个用于压力质量矩阵的 <code>InverseMatrix</code> 对象，我们使用相应的类和我们已经生成的IC预处理器以及用于速度-速度矩阵的AMG预处理器一起设置。注意， <code>Mp_preconditioner</code> 和 <code>Amg_preconditioner</code> 都只是指针，所以我们用 <code>*</code> 来传递实际的预处理对象。

// 一旦预处理程序准备好了，我们就为该块系统创建一个GMRES求解器。由于我们使用的是Trilinos数据结构，我们必须在求解器中设置相应的模板参数。GMRES需要在内部存储每次迭代的临时向量（见 step-22 的结果部分的讨论）&ndash；它可以使用的向量越多，一般来说性能越好。为了控制内存需求，我们将向量的数量设置为100。这意味着在求解器的100次迭代中，每个临时向量都可以被存储。如果求解器需要更频繁地迭代以获得指定的容忍度，它将通过每100次迭代重新开始，在一个减少的向量集上工作。

// 有了这些设置，我们求解系统并在斯托克斯系统中分配约束条件，即悬挂节点和无流体边界条件，以便即使在受约束的道夫下也有适当的解值。最后，我们把迭代次数写到屏幕上。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::solve() 
  { 
    std::cout << "   Solving..." << std::endl; 

    { 
      const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix, 
                                         TrilinosWrappers::PreconditionIC> 
        mp_inverse(stokes_preconditioner_matrix.block(1, 1), 
                   *Mp_preconditioner); 

      const LinearSolvers::BlockSchurPreconditioner< 
        TrilinosWrappers::PreconditionAMG, 
        TrilinosWrappers::PreconditionIC> 
        preconditioner(stokes_matrix, mp_inverse, *Amg_preconditioner); 

      SolverControl solver_control(stokes_matrix.m(), 
                                   1e-6 * stokes_rhs.l2_norm()); 

      SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres( 
        solver_control, 
        SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData(100)); 

      for (unsigned int i = 0; i < stokes_solution.size(); ++i) 
        if (stokes_constraints.is_constrained(i)) 
          stokes_solution(i) = 0; 

      gmres.solve(stokes_matrix, stokes_solution, stokes_rhs, preconditioner); 

      stokes_constraints.distribute(stokes_solution); 

      std::cout << "   " << solver_control.last_step() 
                << " GMRES iterations for Stokes subsystem." << std::endl; 
    } 

// 一旦我们知道了斯托克斯解，我们就可以根据最大速度确定新的时间步长。我们必须这样做以满足CFL条件，因为对流项在温度方程中得到了明确的处理，正如在介绍中所讨论的那样。这里使用的时间步长公式的确切形式将在本程序的结果部分讨论。

// 这里有一个插曲。该公式包含了对速度最大值的除法。然而，在计算开始时，我们有一个恒定的温度场（我们以恒定的温度开始，只有在源作用的第一个时间步长后，它才会变成非恒定的）。恒定温度意味着没有浮力作用，所以速度为零。除以它不可能得到什么好结果。

// 为了避免产生无限的时间步长，我们问最大速度是否非常小（特别是小于我们在接下来的任何时间步长中遇到的值），如果是，我们就不除以零，而是除以一个小值，从而产生一个大的但有限的时间步长。

    old_time_step                 = time_step; 
    const double maximal_velocity = get_maximal_velocity(); 

    if (maximal_velocity >= 0.01) 
      time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree * 
                  GridTools::minimal_cell_diameter(triangulation) / 
                  maximal_velocity; 
    else 
      time_step = 1. / (1.7 * dim * std::sqrt(1. * dim)) / temperature_degree * 
                  GridTools::minimal_cell_diameter(triangulation) / .01; 

    std::cout << "   " 
              << "Time step: " << time_step << std::endl; 

    temperature_solution = old_temperature_solution; 

// 接下来我们用函数  <code>assemble_temperature_system()</code>  设置温度系统和右手边。 知道了温度方程的矩阵和右手边，我们设置了一个预处理程序和一个求解器。温度矩阵是一个质量矩阵（特征值在1左右）加上一个拉普拉斯矩阵（特征值在0和 $ch^{-2}$ 之间）乘以一个与时间步长成正比的小数字  $k_n$  。因此，产生的对称和正定矩阵的特征值在 $[1,1+k_nh^{-2}]$ 范围内（至于常数）。这个矩阵即使对于小的网格尺寸也只是适度的条件不良，我们通过简单的方法得到一个相当好的预处理，例如用一个不完全的Cholesky分解预处理（IC），我们也用它来预处理压力质量矩阵求解器。作为一个求解器，我们选择共轭梯度法CG。和以前一样，我们通过模板参数 <code>TrilinosWrappers::MPI::Vector</code> 告诉求解器使用Trilinos向量。最后，我们求解，分配悬挂节点约束，并写出迭代次数。

    assemble_temperature_system(maximal_velocity); 
    { 
      SolverControl solver_control(temperature_matrix.m(), 
                                   1e-8 * temperature_rhs.l2_norm()); 
      SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 

      TrilinosWrappers::PreconditionIC preconditioner; 
      preconditioner.initialize(temperature_matrix); 

      cg.solve(temperature_matrix, 
               temperature_solution, 
               temperature_rhs, 
               preconditioner); 

      temperature_constraints.distribute(temperature_solution); 

      std::cout << "   " << solver_control.last_step() 
                << " CG iterations for temperature." << std::endl; 

// 在这个函数的结尾，我们在向量中步进并读出最大和最小的温度值，我们也想输出这些值。在本程序的结果部分讨论的确定时间步长的正确常数时，这将非常有用。

      double min_temperature = temperature_solution(0), 
             max_temperature = temperature_solution(0); 
      for (unsigned int i = 0; i < temperature_solution.size(); ++i) 
        { 
          min_temperature = 
            std::min<double>(min_temperature, temperature_solution(i)); 
          max_temperature = 
            std::max<double>(max_temperature, temperature_solution(i)); 
        } 

      std::cout << "   Temperature range: " << min_temperature << ' ' 
                << max_temperature << std::endl; 
    } 
  } 

//  @sect4{BoussinesqFlowProblem::output_results}  

// 该函数将解决方案写入VTK输出文件，用于可视化，每隔10个时间步长就会完成。这通常是一个相当简单的任务，因为deal.II库提供的函数几乎为我们完成了所有的工作。与以前的例子相比，有一个新的函数。我们想把斯托克斯解和温度都看作一个数据集，但是我们已经根据两个不同的DoFHandler对象完成了所有的计算。幸运的是，DataOut类已经准备好处理这个问题。我们所要做的就是不要在一开始就附加一个单一的DoFHandler，然后将其用于所有添加的向量，而是为每个向量分别指定DoFHandler。剩下的就像  step-22  中所做的那样。我们创建解决方案的名称（这些名称将出现在各个组件的可视化程序中）。第一个 <code>dim</code> 分量是矢量速度，然后我们有斯托克斯部分的压力，而温度是标量。这些信息是用DataComponentInterpretation辅助类读出来的。接下来，我们将数据向量与它们的DoFHandler对象连接起来，根据自由度建立补丁，这些补丁是描述可视化程序数据的（子）元素。最后，我们打开一个文件（包括时间步数）并将vtk数据写入其中。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::output_results() const 
  { 
    if (timestep_number % 10 != 0) 
      return; 

    std::vector<std::string> stokes_names(dim, "velocity"); 
    stokes_names.emplace_back("p"); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      stokes_component_interpretation( 
        dim + 1, DataComponentInterpretation::component_is_scalar); 
    for (unsigned int i = 0; i < dim; ++i) 
      stokes_component_interpretation[i] = 
        DataComponentInterpretation::component_is_part_of_vector; 

    DataOut<dim> data_out; 
    data_out.add_data_vector(stokes_dof_handler, 
                             stokes_solution, 
                             stokes_names, 
                             stokes_component_interpretation); 
    data_out.add_data_vector(temperature_dof_handler, 
                             temperature_solution, 
                             "T"); 
    data_out.build_patches(std::min(stokes_degree, temperature_degree)); 

    std::ofstream output("solution-" + 
                         Utilities::int_to_string(timestep_number, 4) + ".vtk"); 
    data_out.write_vtk(output); 
  } 

//  @sect4{BoussinesqFlowProblem::refine_mesh}  

// 这个函数负责处理自适应网格细化。这个函数执行的三个任务是：首先找出需要细化/粗化的单元，然后实际进行细化，并最终在两个不同的网格之间传输解向量。第一个任务是通过对温度使用成熟的凯利误差估计器来实现的（对于这个程序，我们主要关注的是温度，我们需要在高温度梯度的区域保持精确，同时也要避免有太多的数值扩散）。第二项任务是实际进行再塑形。这也只涉及到基本函数，例如 <code>refine_and_coarsen_fixed_fraction</code> ，它可以细化那些具有最大估计误差的单元，这些误差合计占80%，并粗化那些具有最小误差的单元，这些误差合计占10%。

// 如果像这样实施，我们会得到一个不会有太大进展的程序。请记住，我们期望的温度场几乎是不连续的（扩散率 $\kappa$ 毕竟非常小），因此我们可以预期，一个自由适应的网格会越来越细化到大梯度的区域。网格大小的减少将伴随着时间步长的减少，需要大量的时间步长来解决给定的最终时间。这也会导致在几个网格细化周期后，网格的不连续性解决得比开始时好得多。

// 特别是为了防止时间步长的减少和相应的大量时间步长，我们限制了网格的最大细化深度。为此，在细化指标应用于单元格后，我们简单地在最细层的所有单元格上循环，如果它们会导致网格层次过高，则取消对它们的细化选择。

  template <int dim> 
  void 
  BoussinesqFlowProblem<dim>::refine_mesh(const unsigned int max_grid_level) 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate(temperature_dof_handler, 
                                       QGauss<dim - 1>(temperature_degree + 1), 
                                       {}, 
                                       temperature_solution, 
                                       estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.8, 
                                                      0.1); 
    if (triangulation.n_levels() > max_grid_level) 
      for (auto &cell : 
           triangulation.active_cell_iterators_on_level(max_grid_level)) 
        cell->clear_refine_flag(); 

// 作为网格细化的一部分，我们需要将旧的网格中的解决方案向量转移到新的网格中。为此，我们使用SolutionTransfer类，我们必须准备好需要转移到新网格的解向量（一旦完成细化，我们将失去旧的网格，所以转移必须与细化同时发生）。我们肯定需要的是当前温度和旧温度（BDF-2时间步长需要两个旧的解决方案）。由于SolutionTransfer对象只支持在每个dof处理程序中传输一个对象，我们需要在一个数据结构中收集两个温度解决方案。此外，我们也选择转移斯托克斯解，因为我们需要前两个时间步长的速度，其中只有一个是在飞行中计算的。

// 因此，我们为斯托克斯和温度的DoFHandler对象初始化了两个SolutionTransfer对象，将它们附加到旧的dof处理程序中。有了这个，我们就可以准备三角测量和数据向量的细化了（按这个顺序）。

    std::vector<TrilinosWrappers::MPI::Vector> x_temperature(2); 
    x_temperature[0]                            = temperature_solution; 
    x_temperature[1]                            = old_temperature_solution; 
    TrilinosWrappers::MPI::BlockVector x_stokes = stokes_solution; 

    SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> temperature_trans( 
      temperature_dof_handler); 
    SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> stokes_trans( 
      stokes_dof_handler); 

    triangulation.prepare_coarsening_and_refinement(); 
    temperature_trans.prepare_for_coarsening_and_refinement(x_temperature); 
    stokes_trans.prepare_for_coarsening_and_refinement(x_stokes); 

// 现在一切都准备好了，所以进行细化，在新的网格上重新创建dof结构，并初始化矩阵结构和 <code>setup_dofs</code> 函数中的新向量。接下来，我们实际执行网格之间的插值解。我们为温度创建另一份临时向量（现在与新网格相对应），并让插值函数完成这项工作。然后，产生的向量数组被写入各自的向量成员变量中。

// 记住，约束集将在setup_dofs()调用中为新的三角结构进行更新。

    triangulation.execute_coarsening_and_refinement(); 
    setup_dofs(); 

    std::vector<TrilinosWrappers::MPI::Vector> tmp(2); 
    tmp[0].reinit(temperature_solution); 
    tmp[1].reinit(temperature_solution); 
    temperature_trans.interpolate(x_temperature, tmp); 

    temperature_solution     = tmp[0]; 
    old_temperature_solution = tmp[1]; 

// 在解决方案被转移后，我们再对被转移的解决方案实施约束。

    temperature_constraints.distribute(temperature_solution); 
    temperature_constraints.distribute(old_temperature_solution); 

// 对于斯托克斯矢量，一切都一样&ndash;除了我们不需要另一个临时矢量，因为我们只是插值了一个矢量。最后，我们必须告诉程序，矩阵和预处理程序需要重新生成，因为网格已经改变。

    stokes_trans.interpolate(x_stokes, stokes_solution); 

    stokes_constraints.distribute(stokes_solution); 

    rebuild_stokes_matrix         = true; 
    rebuild_temperature_matrices  = true; 
    rebuild_stokes_preconditioner = true; 
  } 

//  @sect4{BoussinesqFlowProblem::run}  

// 这个函数执行Boussinesq程序中的所有基本步骤。它首先设置一个网格（根据空间维度，我们选择一些不同级别的初始细化和额外的自适应细化步骤，然后在 <code>dim</code> 维度上创建一个立方体，并首次设置了道夫。由于我们想用一个自适应细化的网格开始时间步进，我们执行一些预细化步骤，包括所有的装配、求解和细化，但实际上没有在时间上推进。相反，我们使用被人诟病的 <code>goto</code> 语句，在网格细化后立即跳出时间循环，从 <code>start_time_iteration</code> 标签开始的新网格上重新开始。( <code>goto</code> 的使用将在 step-26 中讨论) 。

// 在我们开始之前，我们将初始值投影到网格上，并获得 <code>old_temperature_solution</code> 矢量的第一个数据。然后，我们初始化时间步数和时间步长，开始时间循环。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::run() 
  { 
    const unsigned int initial_refinement     = (dim == 2 ? 4 : 2); 
    const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3); 

    GridGenerator::hyper_cube(triangulation); 
    global_Omega_diameter = GridTools::diameter(triangulation); 

    triangulation.refine_global(initial_refinement); 

    setup_dofs(); 

    unsigned int pre_refinement_step = 0; 

  start_time_iteration: 

    VectorTools::project(temperature_dof_handler, 
                         temperature_constraints, 
                         QGauss<dim>(temperature_degree + 2), 
                         EquationData::TemperatureInitialValues<dim>(), 
                         old_temperature_solution); 

    timestep_number = 0; 
    time_step = old_time_step = 0; 

    double time = 0; 

    do 
      { 
        std::cout << "Timestep " << timestep_number << ":  t=" << time 
                  << std::endl; 

// 时间循环的第一步都是显而易见的；我们组装斯托克斯系统、预处理程序、温度矩阵（矩阵和预处理程序实际上只在我们之前重新处理的情况下发生变化），然后进行求解。在继续下一个时间步骤之前，我们必须检查我们是否应该首先完成预精炼步骤，或者是否应该重新啮合（每五个时间步骤），精炼到一个与初始精炼和预精炼步骤一致的水平。循环的最后一个步骤是推进解，即把解复制到下一个 "较早 "的时间层。

        assemble_stokes_system(); 
        build_stokes_preconditioner(); 
        assemble_temperature_matrix(); 

        solve(); 

        output_results(); 

        std::cout << std::endl; 

        if ((timestep_number == 0) && 
            (pre_refinement_step < n_pre_refinement_steps)) 
          { 
            refine_mesh(initial_refinement + n_pre_refinement_steps); 
            ++pre_refinement_step; 
            goto start_time_iteration; 
          } 
        else if ((timestep_number > 0) && (timestep_number % 5 == 0)) 
          refine_mesh(initial_refinement + n_pre_refinement_steps); 

        time += time_step; 
        ++timestep_number; 

        old_stokes_solution          = stokes_solution; 
        old_old_temperature_solution = old_temperature_solution; 
        old_temperature_solution     = temperature_solution; 
      } 

// 做以上所有的工作，直到我们到达时间100。

    while (time <= 100); 
  } 
} // namespace Step31 

//  @sect3{The <code>main</code> function}  

// 主函数看起来与所有其他程序几乎一样。

// 有一个区别是我们必须要注意的。这个程序使用了Trilinos，而通常情况下，Trilinos被配置为可以使用MPI在%parallel中运行。这并不意味着它<i>has</i>可以在%parallel中运行，事实上这个程序（不像 step-32 ）根本没有尝试使用MPI在%parallel中做任何事情。然而，Trilinos希望MPI系统被初始化。我们通过创建一个类型为 Utilities::MPI::MPI_InitFinalize 的对象来做到这一点，该对象使用给main()的参数（即 <code>argc</code> 和 <code>argv</code> ）初始化MPI（如果可用的话），并在对象超出范围时再次去初始化它。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step31; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 

// 这个程序只能在串行中运行。否则，将抛出一个异常。

      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
                  ExcMessage( 
                    "This program can only be run in serial, use ./step-31")); 

      BoussinesqFlowProblem<2> flow_problem; 
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


