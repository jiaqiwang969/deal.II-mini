

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2008 - 2021 by the deal.II authors 
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
 *          Wolfgang Bangerth, Texas A&M University, 
 *          Timo Heister, University of Goettingen, 2008-2011 
 */ 


// @sect3{Include files}  

// 像往常一样，第一个任务是包括这些著名的deal.II库文件和一些C++头文件的功能。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/work_stream.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/parameter_handler.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/solver_bicgstab.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/block_sparsity_pattern.h> 
#include <deal.II/lac/trilinos_parallel_block_vector.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_block_sparse_matrix.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_solver.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/filtered_iterator.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_dgq.h> 
#include <deal.II/fe/fe_dgp.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/solution_transfer.h> 

#include <fstream> 
#include <iostream> 
#include <limits> 
#include <locale> 
#include <string> 

// 这是唯一一个新的包含文件：它引入了相当于 parallel::distributed::SolutionTransfer 的 dealii::SolutionTransfer 类，用于在网格细化时将解决方案从一个网格带到下一个网格，但在并行分布式三角形计算的情况下。

#include <deal.II/distributed/solution_transfer.h> 

// 以下是用于并行分布式计算的类，在 step-40 中已经全部介绍过。

#include <deal.II/base/index_set.h> 
#include <deal.II/distributed/tria.h> 
#include <deal.II/distributed/grid_refinement.h> 

// 接下来的步骤与之前所有的教程程序一样。我们把所有东西放到一个自己的命名空间中，然后把deal.II的类和函数导入其中。

namespace Step32 
{ 
  using namespace dealii; 
// @sect3{Equation data}  

// 在以下命名空间中，我们定义了描述问题的各种方程数据。这对应于使问题至少有一点现实性的各个方面，并且在介绍中对测试案例的描述中已经详尽地讨论了这些方面。

// 我们从一些具有常数的系数开始（数值后面的注释表示其物理单位）。

  namespace EquationData 
  { 
    constexpr double eta                   = 1e21;    /* Pa s       */ 


    constexpr double kappa                 = 1e-6;    /* m^2 / s    */ 


    constexpr double reference_density     = 3300;    /* kg / m^3   */ 


    constexpr double reference_temperature = 293;     /* K          */ 


    constexpr double expansion_coefficient = 2e-5;    /* 1/K        */ 


    constexpr double specific_heat         = 1250;    /* J / K / kg */ 


    constexpr double radiogenic_heating    = 7.4e-12; /* W / kg     */ 



    constexpr double R0 = 6371000. - 2890000.; /* m          */ 


    constexpr double R1 = 6371000. - 35000.;   /* m          */ 



    constexpr double T0 = 4000 + 273; /* K          */ 


    constexpr double T1 = 700 + 273;  /* K          */ 



// 下一组定义是用于编码密度与温度的函数、重力矢量和温度的初始值的函数。同样，所有这些（以及它们所计算的值）都在介绍中讨论过。

    double density(const double temperature) 
    { 
      return ( 
        reference_density * 
        (1 - expansion_coefficient * (temperature - reference_temperature))); 
    } 

    template <int dim> 
    Tensor<1, dim> gravity_vector(const Point<dim> &p) 
    { 
      const double r = p.norm(); 
      return -(1.245e-6 * r + 7.714e13 / r / r) * p / r; 
    } 

    template <int dim> 
    class TemperatureInitialValues : public Function<dim> 
    { 
    public: 
      TemperatureInitialValues() 
        : Function<dim>(1) 
      {} 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  value) const override; 
    }; 

    template <int dim> 
    double TemperatureInitialValues<dim>::value(const Point<dim> &p, 
                                                const unsigned int) const 
    { 
      const double r = p.norm(); 
      const double h = R1 - R0; 

      const double s = (r - R0) / h; 
      const double q = 
        (dim == 3) ? std::max(0.0, cos(numbers::PI * abs(p(2) / R1))) : 1.0; 
      const double phi = std::atan2(p(0), p(1)); 
      const double tau = s + 0.2 * s * (1 - s) * std::sin(6 * phi) * q; 

      return T0 * (1.0 - tau) + T1 * tau; 
    } 

    template <int dim> 
    void 
    TemperatureInitialValues<dim>::vector_value(const Point<dim> &p, 
                                                Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = TemperatureInitialValues<dim>::value(p, c); 
    } 

// 正如介绍中所提到的，我们需要重新调整压力的比例，以避免动量和质量守恒方程的相对条件不良。比例系数为 $\frac{\eta}{L}$ ，其中 $L$ 是一个典型的长度尺度。通过实验发现，一个好的长度尺度是烟羽的直径，大约是10公里。

    constexpr double pressure_scaling = eta / 10000; 

// 这个命名空间的最后一个数字是一个常数，表示每（平均，热带）年的秒数。我们只在生成屏幕输出时使用它：在内部，这个程序的所有计算都是以SI单位（公斤、米、秒）进行的，但是用秒来写地质学时间产生的数字无法与现实联系起来，所以我们用这里定义的系数转换为年。

    const double year_in_seconds = 60 * 60 * 24 * 365.2425; 

  } // namespace EquationData 

//  @sect3{Preconditioning the Stokes system}  

// 这个命名空间实现了预处理程序。正如介绍中所讨论的，这个预处理程序在一些关键部分与  step-31  中使用的预处理程序不同。具体来说，它是一个右预处理程序，实现了矩阵
//  @f{align*}
//    \left(\begin{array}{cc}A^{-1} & B^T
//                         \\0 & S^{-1}
//  \end{array}\right)
//  @f}
//  中的两个逆矩阵操作由线性求解器近似，或者，如果给这个类的构造函数加上右标志，则由速度块的单个AMG V-循环实现。 <code>vmult</code> 函数的三个代码块实现了与该预处理矩阵的三个块的乘法运算，如果你读过 step-31 或 step-20 中关于组成求解器的讨论，应该是不言自明的。

  namespace LinearSolvers 
  { 
    template <class PreconditionerTypeA, class PreconditionerTypeMp> 
    class BlockSchurPreconditioner : public Subscriptor 
    { 
    public: 
      BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix &S, 
                               const TrilinosWrappers::BlockSparseMatrix &Spre, 
                               const PreconditionerTypeMp &Mppreconditioner, 
                               const PreconditionerTypeA & Apreconditioner, 
                               const bool                  do_solve_A) 
        : stokes_matrix(&S) 
        , stokes_preconditioner_matrix(&Spre) 
        , mp_preconditioner(Mppreconditioner) 
        , a_preconditioner(Apreconditioner) 
        , do_solve_A(do_solve_A) 
      {} 

      void vmult(TrilinosWrappers::MPI::BlockVector &      dst, 
                 const TrilinosWrappers::MPI::BlockVector &src) const 
      { 
        TrilinosWrappers::MPI::Vector utmp(src.block(0)); 

        { 
          SolverControl solver_control(5000, 1e-6 * src.block(1).l2_norm()); 

          SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control); 

          solver.solve(stokes_preconditioner_matrix->block(1, 1), 
                       dst.block(1), 
                       src.block(1), 
                       mp_preconditioner); 

          dst.block(1) *= -1.0; 
        } 

        { 
          stokes_matrix->block(0, 1).vmult(utmp, dst.block(1)); 
          utmp *= -1.0; 
          utmp.add(src.block(0)); 
        } 

        if (do_solve_A == true) 
          { 
            SolverControl solver_control(5000, utmp.l2_norm() * 1e-2); 
            TrilinosWrappers::SolverCG solver(solver_control); 
            solver.solve(stokes_matrix->block(0, 0), 
                         dst.block(0), 
                         utmp, 
                         a_preconditioner); 
          } 
        else 
          a_preconditioner.vmult(dst.block(0), utmp); 
      } 

    private: 
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> 
        stokes_matrix; 
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> 
                                  stokes_preconditioner_matrix; 
      const PreconditionerTypeMp &mp_preconditioner; 
      const PreconditionerTypeA & a_preconditioner; 
      const bool                  do_solve_A; 
    }; 
  } // namespace LinearSolvers 

//  @sect3{Definition of assembly data structures}  

// 如介绍中所述，我们将使用 @ref threads 模块中讨论的WorkStream机制来实现单台机器的处理器之间的并行操作。WorkStream类要求数据在两种数据结构中传递，一种是用于抓取数据，一种是将数据从装配函数传递到将本地贡献复制到全局对象的函数。

// 下面的命名空间（以及两个子命名空间）包含了服务于这一目的的数据结构的集合，介绍中讨论的四种操作中的每一种都有一对，我们将想把它们并行化。每个装配例程都会得到两组数据：一个是Scratch数组，收集所有用于计算单元格贡献的类和数组，另一个是CopyData数组，保存将被写入全局矩阵的本地矩阵和向量。而CopyData是一个容器，用来存放最终写入全局矩阵和向量的数据（因此是绝对必要的），Scratch数组只是出于性能的考虑而存在；在每个单元上设置一个FEValues对象，要比只创建一次并更新一些导数数据要昂贵得多。

//  Step-31 有四个汇编程序。一个用于斯托克斯系统的预处理矩阵，一个用于斯托克斯矩阵和右手边，一个用于温度矩阵，一个用于温度方程的右手边。我们在这里使用 <code>struct</code> 环境为这四个汇编组件中的每一个组织从头数组和CopyData对象（因为我们认为这些是我们传递的临时对象，而不是实现自己功能的类，尽管这是区分 <code>struct</code>s and <code>class</code> es的一个比较主观的观点）。

// 关于Scratch对象，每个结构都配备了一个构造函数，可以使用 @ref FiniteElement 、正交、 @ref Mapping （描述弯曲边界的插值）和 @ref UpdateFlags 实例创建一个 @ref FEValues 对象。此外，我们手动实现了一个复制构造函数（因为FEValues类本身是不可复制的），并提供了一些额外的矢量字段，用于在计算局部贡献时保存中间数据。

// 让我们从抓取数组开始，特别是用于组装斯托克斯预处理程序的数组。

  namespace Assembly 
  { 
    namespace Scratch 
    { 
      template <int dim> 
      struct StokesPreconditioner 
      { 
        StokesPreconditioner(const FiniteElement<dim> &stokes_fe, 
                             const Quadrature<dim> &   stokes_quadrature, 
                             const Mapping<dim> &      mapping, 
                             const UpdateFlags         update_flags); 

        StokesPreconditioner(const StokesPreconditioner &data); 

        FEValues<dim> stokes_fe_values; 

        std::vector<Tensor<2, dim>> grad_phi_u; 
        std::vector<double>         phi_p; 
      }; 

      template <int dim> 
      StokesPreconditioner<dim>::StokesPreconditioner( 
        const FiniteElement<dim> &stokes_fe, 
        const Quadrature<dim> &   stokes_quadrature, 
        const Mapping<dim> &      mapping, 
        const UpdateFlags         update_flags) 
        : stokes_fe_values(mapping, stokes_fe, stokes_quadrature, update_flags) 
        , grad_phi_u(stokes_fe.n_dofs_per_cell()) 
        , phi_p(stokes_fe.n_dofs_per_cell()) 
      {} 

      template <int dim> 
      StokesPreconditioner<dim>::StokesPreconditioner( 
        const StokesPreconditioner &scratch) 
        : stokes_fe_values(scratch.stokes_fe_values.get_mapping(), 
                           scratch.stokes_fe_values.get_fe(), 
                           scratch.stokes_fe_values.get_quadrature(), 
                           scratch.stokes_fe_values.get_update_flags()) 
        , grad_phi_u(scratch.grad_phi_u) 
        , phi_p(scratch.phi_p) 
      {} 

// 下一个是用于组装完整的斯托克斯系统的从头对象。请注意，我们从上面的StokesPreconditioner类派生出StokesSystem scratch类。我们这样做是因为所有用于组装预处理程序的对象也需要用于实际的矩阵系统和右手边，还有一些额外的数据。这使得程序更加紧凑。还需要注意的是，斯托克斯系统的装配和进一步的温度右手边分别需要温度和速度的数据，所以我们实际上需要两个FEValues对象来处理这两种情况。

      template <int dim> 
      struct StokesSystem : public StokesPreconditioner<dim> 
      { 
        StokesSystem(const FiniteElement<dim> &stokes_fe, 
                     const Mapping<dim> &      mapping, 
                     const Quadrature<dim> &   stokes_quadrature, 
                     const UpdateFlags         stokes_update_flags, 
                     const FiniteElement<dim> &temperature_fe, 
                     const UpdateFlags         temperature_update_flags); 

        StokesSystem(const StokesSystem<dim> &data); 

        FEValues<dim> temperature_fe_values; 

        std::vector<Tensor<1, dim>>          phi_u; 
        std::vector<SymmetricTensor<2, dim>> grads_phi_u; 
        std::vector<double>                  div_phi_u; 

        std::vector<double> old_temperature_values; 
      }; 

      template <int dim> 
      StokesSystem<dim>::StokesSystem( 
        const FiniteElement<dim> &stokes_fe, 
        const Mapping<dim> &      mapping, 
        const Quadrature<dim> &   stokes_quadrature, 
        const UpdateFlags         stokes_update_flags, 
        const FiniteElement<dim> &temperature_fe, 
        const UpdateFlags         temperature_update_flags) 
        : StokesPreconditioner<dim>(stokes_fe, 
                                    stokes_quadrature, 
                                    mapping, 
                                    stokes_update_flags) 
        , temperature_fe_values(mapping, 
                                temperature_fe, 
                                stokes_quadrature, 
                                temperature_update_flags) 
        , phi_u(stokes_fe.n_dofs_per_cell()) 
        , grads_phi_u(stokes_fe.n_dofs_per_cell()) 
        , div_phi_u(stokes_fe.n_dofs_per_cell()) 
        , old_temperature_values(stokes_quadrature.size()) 
      {} 

      template <int dim> 
      StokesSystem<dim>::StokesSystem(const StokesSystem<dim> &scratch) 
        : StokesPreconditioner<dim>(scratch) 
        , temperature_fe_values( 
            scratch.temperature_fe_values.get_mapping(), 
            scratch.temperature_fe_values.get_fe(), 
            scratch.temperature_fe_values.get_quadrature(), 
            scratch.temperature_fe_values.get_update_flags()) 
        , phi_u(scratch.phi_u) 
        , grads_phi_u(scratch.grads_phi_u) 
        , div_phi_u(scratch.div_phi_u) 
        , old_temperature_values(scratch.old_temperature_values) 
      {} 

// 在定义了用于组装斯托克斯系统的对象之后，我们对温度系统所需的矩阵的组装也做了同样的工作。一般的结构是非常相似的。

      template <int dim> 
      struct TemperatureMatrix 
      { 
        TemperatureMatrix(const FiniteElement<dim> &temperature_fe, 
                          const Mapping<dim> &      mapping, 
                          const Quadrature<dim> &   temperature_quadrature); 

        TemperatureMatrix(const TemperatureMatrix &data); 

        FEValues<dim> temperature_fe_values; 

        std::vector<double>         phi_T; 
        std::vector<Tensor<1, dim>> grad_phi_T; 
      }; 

      template <int dim> 
      TemperatureMatrix<dim>::TemperatureMatrix( 
        const FiniteElement<dim> &temperature_fe, 
        const Mapping<dim> &      mapping, 
        const Quadrature<dim> &   temperature_quadrature) 
        : temperature_fe_values(mapping, 
                                temperature_fe, 
                                temperature_quadrature, 
                                update_values | update_gradients | 
                                  update_JxW_values) 
        , phi_T(temperature_fe.n_dofs_per_cell()) 
        , grad_phi_T(temperature_fe.n_dofs_per_cell()) 
      {} 

      template <int dim> 
      TemperatureMatrix<dim>::TemperatureMatrix( 
        const TemperatureMatrix &scratch) 
        : temperature_fe_values( 
            scratch.temperature_fe_values.get_mapping(), 
            scratch.temperature_fe_values.get_fe(), 
            scratch.temperature_fe_values.get_quadrature(), 
            scratch.temperature_fe_values.get_update_flags()) 
        , phi_T(scratch.phi_T) 
        , grad_phi_T(scratch.grad_phi_T) 
      {} 

// 最后的划痕对象被用于温度系统右侧的装配。这个对象比上面的对象要大得多，因为有更多的量进入温度方程右边的计算中。特别是，前两个时间步骤的温度值和梯度需要在正交点评估，还有速度和应变率（即速度的对称梯度），它们作为摩擦加热项进入右侧。尽管有很多条款，但以下内容应该是不言自明的。

      template <int dim> 
      struct TemperatureRHS 
      { 
        TemperatureRHS(const FiniteElement<dim> &temperature_fe, 
                       const FiniteElement<dim> &stokes_fe, 
                       const Mapping<dim> &      mapping, 
                       const Quadrature<dim> &   quadrature); 

        TemperatureRHS(const TemperatureRHS &data); 

        FEValues<dim> temperature_fe_values; 
        FEValues<dim> stokes_fe_values; 

        std::vector<double>         phi_T; 
        std::vector<Tensor<1, dim>> grad_phi_T; 

        std::vector<Tensor<1, dim>> old_velocity_values; 
        std::vector<Tensor<1, dim>> old_old_velocity_values; 

        std::vector<SymmetricTensor<2, dim>> old_strain_rates; 
        std::vector<SymmetricTensor<2, dim>> old_old_strain_rates; 

        std::vector<double>         old_temperature_values; 
        std::vector<double>         old_old_temperature_values; 
        std::vector<Tensor<1, dim>> old_temperature_grads; 
        std::vector<Tensor<1, dim>> old_old_temperature_grads; 
        std::vector<double>         old_temperature_laplacians; 
        std::vector<double>         old_old_temperature_laplacians; 
      }; 

      template <int dim> 
      TemperatureRHS<dim>::TemperatureRHS( 
        const FiniteElement<dim> &temperature_fe, 
        const FiniteElement<dim> &stokes_fe, 
        const Mapping<dim> &      mapping, 
        const Quadrature<dim> &   quadrature) 
        : temperature_fe_values(mapping, 
                                temperature_fe, 
                                quadrature, 
                                update_values | update_gradients | 
                                  update_hessians | update_quadrature_points | 
                                  update_JxW_values) 
        , stokes_fe_values(mapping, 
                           stokes_fe, 
                           quadrature, 
                           update_values | update_gradients) 
        , phi_T(temperature_fe.n_dofs_per_cell()) 
        , grad_phi_T(temperature_fe.n_dofs_per_cell()) 
        , 

        old_velocity_values(quadrature.size()) 
        , old_old_velocity_values(quadrature.size()) 
        , old_strain_rates(quadrature.size()) 
        , old_old_strain_rates(quadrature.size()) 
        , 

        old_temperature_values(quadrature.size()) 
        , old_old_temperature_values(quadrature.size()) 
        , old_temperature_grads(quadrature.size()) 
        , old_old_temperature_grads(quadrature.size()) 
        , old_temperature_laplacians(quadrature.size()) 
        , old_old_temperature_laplacians(quadrature.size()) 
      {} 

      template <int dim> 
      TemperatureRHS<dim>::TemperatureRHS(const TemperatureRHS &scratch) 
        : temperature_fe_values( 
            scratch.temperature_fe_values.get_mapping(), 
            scratch.temperature_fe_values.get_fe(), 
            scratch.temperature_fe_values.get_quadrature(), 
            scratch.temperature_fe_values.get_update_flags()) 
        , stokes_fe_values(scratch.stokes_fe_values.get_mapping(), 
                           scratch.stokes_fe_values.get_fe(), 
                           scratch.stokes_fe_values.get_quadrature(), 
                           scratch.stokes_fe_values.get_update_flags()) 
        , phi_T(scratch.phi_T) 
        , grad_phi_T(scratch.grad_phi_T) 
        , 

        old_velocity_values(scratch.old_velocity_values) 
        , old_old_velocity_values(scratch.old_old_velocity_values) 
        , old_strain_rates(scratch.old_strain_rates) 
        , old_old_strain_rates(scratch.old_old_strain_rates) 
        , 

        old_temperature_values(scratch.old_temperature_values) 
        , old_old_temperature_values(scratch.old_old_temperature_values) 
        , old_temperature_grads(scratch.old_temperature_grads) 
        , old_old_temperature_grads(scratch.old_old_temperature_grads) 
        , old_temperature_laplacians(scratch.old_temperature_laplacians) 
        , old_old_temperature_laplacians(scratch.old_old_temperature_laplacians) 
      {} 
    } // namespace Scratch 

// CopyData对象比Scratch对象更简单，因为它们所要做的就是存储本地计算的结果，直到它们可以被复制到全局矩阵或向量对象中。因此，这些结构只需要提供一个构造函数，一个复制操作，以及一些用于本地矩阵、本地向量和本地与全局自由度之间关系的数组（又称 <code>local_dof_indices</code> ）。同样，我们为我们将使用WorkStream类并行化的四个操作中的每一个都有一个这样的结构。

    namespace CopyData 
    { 
      template <int dim> 
      struct StokesPreconditioner 
      { 
        StokesPreconditioner(const FiniteElement<dim> &stokes_fe); 
        StokesPreconditioner(const StokesPreconditioner &data); 
        StokesPreconditioner &operator=(const StokesPreconditioner &) = default; 

        FullMatrix<double>                   local_matrix; 
        std::vector<types::global_dof_index> local_dof_indices; 
      }; 

      template <int dim> 
      StokesPreconditioner<dim>::StokesPreconditioner( 
        const FiniteElement<dim> &stokes_fe) 
        : local_matrix(stokes_fe.n_dofs_per_cell(), stokes_fe.n_dofs_per_cell()) 
        , local_dof_indices(stokes_fe.n_dofs_per_cell()) 
      {} 

      template <int dim> 
      StokesPreconditioner<dim>::StokesPreconditioner( 
        const StokesPreconditioner &data) 
        : local_matrix(data.local_matrix) 
        , local_dof_indices(data.local_dof_indices) 
      {} 

      template <int dim> 
      struct StokesSystem : public StokesPreconditioner<dim> 
      { 
        StokesSystem(const FiniteElement<dim> &stokes_fe); 

        Vector<double> local_rhs; 
      }; 

      template <int dim> 
      StokesSystem<dim>::StokesSystem(const FiniteElement<dim> &stokes_fe) 
        : StokesPreconditioner<dim>(stokes_fe) 
        , local_rhs(stokes_fe.n_dofs_per_cell()) 
      {} 

      template <int dim> 
      struct TemperatureMatrix 
      { 
        TemperatureMatrix(const FiniteElement<dim> &temperature_fe); 

        FullMatrix<double>                   local_mass_matrix; 
        FullMatrix<double>                   local_stiffness_matrix; 
        std::vector<types::global_dof_index> local_dof_indices; 
      }; 

      template <int dim> 
      TemperatureMatrix<dim>::TemperatureMatrix( 
        const FiniteElement<dim> &temperature_fe) 
        : local_mass_matrix(temperature_fe.n_dofs_per_cell(), 
                            temperature_fe.n_dofs_per_cell()) 
        , local_stiffness_matrix(temperature_fe.n_dofs_per_cell(), 
                                 temperature_fe.n_dofs_per_cell()) 
        , local_dof_indices(temperature_fe.n_dofs_per_cell()) 
      {} 

      template <int dim> 
      struct TemperatureRHS 
      { 
        TemperatureRHS(const FiniteElement<dim> &temperature_fe); 

        Vector<double>                       local_rhs; 
        std::vector<types::global_dof_index> local_dof_indices; 
        FullMatrix<double>                   matrix_for_bc; 
      }; 

      template <int dim> 
      TemperatureRHS<dim>::TemperatureRHS( 
        const FiniteElement<dim> &temperature_fe) 
        : local_rhs(temperature_fe.n_dofs_per_cell()) 
        , local_dof_indices(temperature_fe.n_dofs_per_cell()) 
        , matrix_for_bc(temperature_fe.n_dofs_per_cell(), 
                        temperature_fe.n_dofs_per_cell()) 
      {} 
    } // namespace CopyData 
  }   // namespace Assembly 

//  @sect3{The <code>BoussinesqFlowProblem</code> class template}  

// 这是主类的声明。它与 step-31 非常相似，但有一些区别我们将在下面评论。

// 该类的顶部与 step-31 中的内容基本相同，列出了公共方法和一组做重活的私有函数。与 step-31 相比，这部分只增加了两个：计算所有单元的最大CFL数的函数 <code>get_cfl_number()</code> ，然后我们根据它计算全局时间步长；以及用于计算熵值稳定的函数 <code>get_entropy_variation()</code> 。它类似于我们在 step-31 中用于此目的的 <code>get_extrapolated_temperature_range()</code> ，但它的工作对象是熵而不是温度。

  template <int dim> 
  class BoussinesqFlowProblem 
  { 
  public: 
    struct Parameters; 
    BoussinesqFlowProblem(Parameters &parameters); 
    void run(); 

  private: 
    void   setup_dofs(); 
    void   assemble_stokes_preconditioner(); 
    void   build_stokes_preconditioner(); 
    void   assemble_stokes_system(); 
    void   assemble_temperature_matrix(); 
    void   assemble_temperature_system(const double maximal_velocity); 
    double get_maximal_velocity() const; 
    double get_cfl_number() const; 
    double get_entropy_variation(const double average_temperature) const; 
    std::pair<double, double> get_extrapolated_temperature_range() const; 
    void                      solve(); 
    void                      output_results(); 
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
      const std::vector<SymmetricTensor<2, dim>> &old_strain_rates, 
      const std::vector<SymmetricTensor<2, dim>> &old_old_strain_rates, 
      const double                                global_u_infty, 
      const double                                global_T_variation, 
      const double                                average_temperature, 
      const double                                global_entropy_variation, 
      const double                                cell_diameter) const; 

  public: 

// 第一个重要的新组件是根据介绍中的讨论为参数定义了一个结构。这个结构是在构建这个对象的过程中通过读取参数文件来初始化的。

    struct Parameters 
    { 
      Parameters(const std::string &parameter_filename); 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 

      double end_time; 

      unsigned int initial_global_refinement; 
      unsigned int initial_adaptive_refinement; 

      bool         generate_graphical_output; 
      unsigned int graphical_output_interval; 

      unsigned int adaptive_refinement_interval; 

      double stabilization_alpha; 
      double stabilization_c_R; 
      double stabilization_beta; 

      unsigned int stokes_velocity_degree; 
      bool         use_locally_conservative_discretization; 

      unsigned int temperature_degree; 
    }; 

  private: 
    Parameters &parameters; 

//  <code>pcout</code> （用于<i>%parallel <code>std::cout</code></i>）对象被用来简化输出的书写：每个MPI进程都可以像往常一样使用它来产生输出，但由于这些进程中的每一个都会（希望）产生相同的输出，它只是被重复了许多次；使用ConditionalOStream类，只有一个MPI进程产生的输出会真正被打印到屏幕上，而所有其他线程的输出将只是被遗忘。

    ConditionalOStream pcout; 

// 下面的成员变量将再次与 step-31 中的成员变量相似（也与其他教程程序相似）。正如介绍中提到的，我们完全分布计算，所以我们将不得不使用 parallel::distributed::Triangulation 类（见 step-40  ），但这些变量的其余部分相当标准，有两个例外。



// --  <code>mapping</code> 这个变量是用来表示高阶多项式映射的。正如在介绍中提到的，我们在通过正交形成积分时使用这个映射，用于所有与我们域的内边界或外边界相邻的单元，其中边界是弯曲的。



// - 在命名混乱的情况下，你会注意到下面一些来自命名空间TrilinosWrappers的变量取自命名空间 TrilinosWrappers::MPI （比如右手边的向量），而其他变量则不是（比如各种矩阵）。这是由于遗留的原因。我们经常需要查询任意正交点的速度和温度；因此，每当我们需要访问与本地相关但属于另一个处理器的自由度时，我们不是导入矢量的幽灵信息，而是以%并行方式求解线性系统，但随后立即初始化一个矢量，包括求解的幽灵条目，以便进一步处理。因此，各种 <code>*_solution</code> 向量在以%parallel求解各自的线性系统后立即被填充，并且总是包含所有 @ref GlossLocallyRelevantDof "本地相关自由度 "的值；我们从求解过程中获得的完全分布的向量，只包含 @ref GlossLocallyOwnedDof "本地拥有的自由度"，在求解过程后，在我们将相关值复制到成员变量向量后立即销毁。

    parallel::distributed::Triangulation<dim> triangulation; 
    double                                    global_Omega_diameter; 

    const MappingQ<dim> mapping; 

    const FESystem<dim>       stokes_fe; 
    DoFHandler<dim>           stokes_dof_handler; 
    AffineConstraints<double> stokes_constraints; 

    TrilinosWrappers::BlockSparseMatrix stokes_matrix; 
    TrilinosWrappers::BlockSparseMatrix stokes_preconditioner_matrix; 

    TrilinosWrappers::MPI::BlockVector stokes_solution; 
    TrilinosWrappers::MPI::BlockVector old_stokes_solution; 
    TrilinosWrappers::MPI::BlockVector stokes_rhs; 

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

    std::shared_ptr<TrilinosWrappers::PreconditionAMG>    Amg_preconditioner; 
    std::shared_ptr<TrilinosWrappers::PreconditionJacobi> Mp_preconditioner; 
    std::shared_ptr<TrilinosWrappers::PreconditionJacobi> T_preconditioner; 

    bool rebuild_stokes_matrix; 
    bool rebuild_stokes_preconditioner; 
    bool rebuild_temperature_matrices; 
    bool rebuild_temperature_preconditioner; 

// 下一个成员变量， <code>computing_timer</code> 是用来方便地计算在某些重复输入的代码 "部分 "所花费的计算时间。例如，我们将进入（和离开）斯托克斯矩阵装配的部分，并希望在所有的时间步骤中累积在这部分花费的运行时间。每隔一段时间，以及在程序结束时（通过TimerOutput类的析构器），我们将产生一个很好的总结，即在不同部分花费的时间，我们把这个程序的运行时间归类为不同部分。

    TimerOutput computing_timer; 

// 在这些成员变量之后，我们有一些辅助函数，这些函数已经从上面列出的那些函数中分解出来。具体来说，首先有三个我们从 <code>setup_dofs</code> 中调用的函数，然后是做线性系统组装的函数。

    void setup_stokes_matrix( 
      const std::vector<IndexSet> &stokes_partitioning, 
      const std::vector<IndexSet> &stokes_relevant_partitioning); 
    void setup_stokes_preconditioner( 
      const std::vector<IndexSet> &stokes_partitioning, 
      const std::vector<IndexSet> &stokes_relevant_partitioning); 
    void setup_temperature_matrices( 
      const IndexSet &temperature_partitioning, 
      const IndexSet &temperature_relevant_partitioning); 

// 遵循 @ref MTWorkStream "基于任务的并行化 "范式，我们将所有的汇编例程分成两部分：第一部分可以在某个单元上做所有的计算，而不需要照顾其他线程；第二部分（就是将本地数据写入全局矩阵和向量中），每次只能由一个线程进入。为了实现这一点，我们为这一程序中使用的所有四个汇编例程的这两个步骤分别提供了函数。下面的八个函数正是这样做的。

    void local_assemble_stokes_preconditioner( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      Assembly::Scratch::StokesPreconditioner<dim> &        scratch, 
      Assembly::CopyData::StokesPreconditioner<dim> &       data); 

    void copy_local_to_global_stokes_preconditioner( 
      const Assembly::CopyData::StokesPreconditioner<dim> &data); 

    void local_assemble_stokes_system( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      Assembly::Scratch::StokesSystem<dim> &                scratch, 
      Assembly::CopyData::StokesSystem<dim> &               data); 

    void copy_local_to_global_stokes_system( 
      const Assembly::CopyData::StokesSystem<dim> &data); 

    void local_assemble_temperature_matrix( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      Assembly::Scratch::TemperatureMatrix<dim> &           scratch, 
      Assembly::CopyData::TemperatureMatrix<dim> &          data); 

    void copy_local_to_global_temperature_matrix( 
      const Assembly::CopyData::TemperatureMatrix<dim> &data); 

    void local_assemble_temperature_rhs( 
      const std::pair<double, double> global_T_range, 
      const double                    global_max_velocity, 
      const double                    global_entropy_variation, 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      Assembly::Scratch::TemperatureRHS<dim> &              scratch, 
      Assembly::CopyData::TemperatureRHS<dim> &             data); 

    void copy_local_to_global_temperature_rhs( 
      const Assembly::CopyData::TemperatureRHS<dim> &data); 

// 最后，我们向前声明一个成员类，我们将在以后定义这个成员类，它将被用来从我们的解决方案向量中计算一些数量，我们希望将这些数量放入输出文件中，以便进行可视化。

    class Postprocessor; 
  }; 
// @sect3{BoussinesqFlowProblem class implementation}  
// @sect4{BoussinesqFlowProblem::Parameters}  

// 这里是对斯托克斯问题的参数的定义。我们允许设置模拟的结束时间、细化水平（包括全局细化和自适应细化，总的来说就是允许单元的最大细化水平），以及细化的时间间隔。

// 然后，我们让用户指定稳定参数的常数（如介绍中所讨论的）、斯托克斯速度空间的多项式程度、是否对压力使用基于FE_DGP元素的局部保守离散化（对压力使用FE_Q元素）、以及温度插值的多项式程度。

// 构造函数检查是否有有效的输入文件（如果没有，将写一个带有默认参数的文件），并最终解析参数。

  template <int dim> 
  BoussinesqFlowProblem<dim>::Parameters::Parameters( 
    const std::string &parameter_filename) 
    : end_time(1e8) 
    , initial_global_refinement(2) 
    , initial_adaptive_refinement(2) 
    , adaptive_refinement_interval(10) 
    , stabilization_alpha(2) 
    , stabilization_c_R(0.11) 
    , stabilization_beta(0.078) 
    , stokes_velocity_degree(2) 
    , use_locally_conservative_discretization(true) 
    , temperature_degree(2) 
  { 
    ParameterHandler prm; 
    BoussinesqFlowProblem<dim>::Parameters::declare_parameters(prm); 

    std::ifstream parameter_file(parameter_filename); 

    if (!parameter_file) 
      { 
        parameter_file.close(); 

        std::ofstream parameter_out(parameter_filename); 
        prm.print_parameters(parameter_out, ParameterHandler::Text); 

        AssertThrow( 
          false, 
          ExcMessage( 
            "Input parameter file <" + parameter_filename + 
            "> not found. Creating a template file of the same name.")); 
      } 

    prm.parse_input(parameter_file); 
    parse_parameters(prm); 
  } 

// 接下来我们有一个函数，声明我们在输入文件中期望的参数，以及它们的数据类型、默认值和描述。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::Parameters::declare_parameters( 
    ParameterHandler &prm) 
  { 
    prm.declare_entry("End time", 
                      "1e8", 
                      Patterns::Double(0), 
                      "The end time of the simulation in years."); 
    prm.declare_entry("Initial global refinement", 
                      "2", 
                      Patterns::Integer(0), 
                      "The number of global refinement steps performed on " 
                      "the initial coarse mesh, before the problem is first " 
                      "solved there."); 
    prm.declare_entry("Initial adaptive refinement", 
                      "2", 
                      Patterns::Integer(0), 
                      "The number of adaptive refinement steps performed after " 
                      "initial global refinement."); 
    prm.declare_entry("Time steps between mesh refinement", 
                      "10", 
                      Patterns::Integer(1), 
                      "The number of time steps after which the mesh is to be " 
                      "adapted based on computed error indicators."); 
    prm.declare_entry("Generate graphical output", 
                      "false", 
                      Patterns::Bool(), 
                      "Whether graphical output is to be generated or not. " 
                      "You may not want to get graphical output if the number " 
                      "of processors is large."); 
    prm.declare_entry("Time steps between graphical output", 
                      "50", 
                      Patterns::Integer(1), 
                      "The number of time steps between each generation of " 
                      "graphical output files."); 

    prm.enter_subsection("Stabilization parameters"); 
    { 
      prm.declare_entry("alpha", 
                        "2", 
                        Patterns::Double(1, 2), 
                        "The exponent in the entropy viscosity stabilization."); 
      prm.declare_entry("c_R", 
                        "0.11", 
                        Patterns::Double(0), 
                        "The c_R factor in the entropy viscosity " 
                        "stabilization."); 
      prm.declare_entry("beta", 
                        "0.078", 
                        Patterns::Double(0), 
                        "The beta factor in the artificial viscosity " 
                        "stabilization. An appropriate value for 2d is 0.052 " 
                        "and 0.078 for 3d."); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Discretization"); 
    { 
      prm.declare_entry( 
        "Stokes velocity polynomial degree", 
        "2", 
        Patterns::Integer(1), 
        "The polynomial degree to use for the velocity variables " 
        "in the Stokes system."); 
      prm.declare_entry( 
        "Temperature polynomial degree", 
        "2", 
        Patterns::Integer(1), 
        "The polynomial degree to use for the temperature variable."); 
      prm.declare_entry( 
        "Use locally conservative discretization", 
        "true", 
        Patterns::Bool(), 
        "Whether to use a Stokes discretization that is locally " 
        "conservative at the expense of a larger number of degrees " 
        "of freedom, or to go with a cheaper discretization " 
        "that does not locally conserve mass (although it is " 
        "globally conservative."); 
    } 
    prm.leave_subsection(); 
  } 

// 然后，我们需要一个函数来读取我们通过读取输入文件得到的ParameterHandler对象的内容，并将结果放入储存我们之前声明的参数值的变量。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::Parameters::parse_parameters( 
    ParameterHandler &prm) 
  { 
    end_time                  = prm.get_double("End time"); 
    initial_global_refinement = prm.get_integer("Initial global refinement"); 
    initial_adaptive_refinement = 
      prm.get_integer("Initial adaptive refinement"); 

    adaptive_refinement_interval = 
      prm.get_integer("Time steps between mesh refinement"); 

    generate_graphical_output = prm.get_bool("Generate graphical output"); 
    graphical_output_interval = 
      prm.get_integer("Time steps between graphical output"); 

    prm.enter_subsection("Stabilization parameters"); 
    { 
      stabilization_alpha = prm.get_double("alpha"); 
      stabilization_c_R   = prm.get_double("c_R"); 
      stabilization_beta  = prm.get_double("beta"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Discretization"); 
    { 
      stokes_velocity_degree = 
        prm.get_integer("Stokes velocity polynomial degree"); 
      temperature_degree = prm.get_integer("Temperature polynomial degree"); 
      use_locally_conservative_discretization = 
        prm.get_bool("Use locally conservative discretization"); 
    } 
    prm.leave_subsection(); 
  } 

//  @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}  

// 该问题的构造函数与  step-31  中的构造函数非常相似。不同的是%并行通信。Trilinos使用消息传递接口（MPI）进行数据分配。当进入BoussinesqFlowProblem类时，我们必须决定如何进行并行化。我们选择一个相当简单的策略，让所有正在运行程序的处理器一起工作，由通信器  <code>MPI_COMM_WORLD</code>  指定。接下来，我们创建输出流（就像我们在 step-18 中已经做的那样），它只在第一个MPI进程上产生输出，而在其他所有进程上则完全不考虑。这个想法的实现是在 <code>pcout</code> 得到一个真实参数时检查进程号，它使用 <code>std::cout</code> 流进行输出。例如，如果我们是一个处理器五，那么我们将给出一个  <code>false</code> argument to <code>pcout</code>  ，这意味着该处理器的输出将不会被打印。除了映射对象（我们对其使用4度的多项式），除了最后的成员变量外，其他都与  step-31  中的完全相同。

// 这个最后的对象，TimerOutput对象，然后被告知限制输出到 <code>pcout</code> 流（处理器0），然后我们指定要在程序结束时得到一个汇总表，该表显示我们的壁挂时钟时间（而不是CPU时间）。我们还将在下面的 <code>run()</code> 函数中手动请求每隔这么多时间步的中间总结。

  template <int dim> 
  BoussinesqFlowProblem<dim>::BoussinesqFlowProblem(Parameters &parameters_) 
    : parameters(parameters_) 
    , pcout(std::cout, (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)) 
    , 

    triangulation(MPI_COMM_WORLD, 
                  typename Triangulation<dim>::MeshSmoothing( 
                    Triangulation<dim>::smoothing_on_refinement | 
                    Triangulation<dim>::smoothing_on_coarsening)) 
    , 

    global_Omega_diameter(0.) 
    , 

    mapping(4) 
    , 

    stokes_fe(FE_Q<dim>(parameters.stokes_velocity_degree), 
              dim, 
              (parameters.use_locally_conservative_discretization ? 
                 static_cast<const FiniteElement<dim> &>( 
                   FE_DGP<dim>(parameters.stokes_velocity_degree - 1)) : 
                 static_cast<const FiniteElement<dim> &>( 
                   FE_Q<dim>(parameters.stokes_velocity_degree - 1))), 
              1) 
    , 

    stokes_dof_handler(triangulation) 
    , 

    temperature_fe(parameters.temperature_degree) 
    , temperature_dof_handler(triangulation) 
    , 

    time_step(0) 
    , old_time_step(0) 
    , timestep_number(0) 
    , rebuild_stokes_matrix(true) 
    , rebuild_stokes_preconditioner(true) 
    , rebuild_temperature_matrices(true) 
    , rebuild_temperature_preconditioner(true) 
    , 

    computing_timer(MPI_COMM_WORLD, 
                    pcout, 
                    TimerOutput::summary, 
                    TimerOutput::wall_times) 
  {} 

//  @sect4{The BoussinesqFlowProblem helper functions}  
// @sect5{BoussinesqFlowProblem::get_maximal_velocity}  

// 除了两个小细节外，计算速度全局最大值的函数与 step-31 中的相同。第一个细节实际上是所有在三角形的所有单元上实现循环的函数所共有的。当以%并行方式操作时，每个处理器只能处理一大块单元，因为每个处理器只拥有整个三角结构的某一部分。我们要处理的这块单元是通过所谓的 <code>subdomain_id</code> 来确定的，正如我们在 step-18 中做的那样。因此，我们需要改变的是只对当前进程所拥有的单元格（相对于幽灵或人造单元格）进行与单元格相关的操作，即对子域id等于进程ID的数字。由于这是一个常用的操作，所以这个操作有一个快捷方式：我们可以用 <code>cell-@>is_locally_owned()</code> 询问单元格是否为当前处理器所拥有。

// 第二个区别是我们计算最大值的方式。以前，我们可以简单地有一个 <code>double</code> 变量，在每个单元的每个正交点上进行检查。现在，我们必须更加小心，因为每个处理器只对单元格的一个子集进行操作。我们要做的是，首先让每个处理器计算其单元中的最大值，然后做一个全局通信操作 <code>Utilities::MPI::max</code> ，计算各个处理器所有最大值中的最大值。MPI提供了这样的调用，但更简单的是使用MPI通信器对象在命名空间 Utilities::MPI 中使用相应的函数，因为即使我们没有MPI并且只在一台机器上工作，这也会做正确的事情。对 <code>Utilities::MPI::max</code> 的调用需要两个参数，即本地最大值（input）和MPI通信器，在这个例子中是MPI_COMM_WORLD。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::get_maximal_velocity() const 
  { 
    const QIterated<dim> quadrature_formula(QTrapezoid<1>(), 
                                            parameters.stokes_velocity_degree); 
    const unsigned int   n_q_points = quadrature_formula.size(); 

    FEValues<dim>               fe_values(mapping, 
                            stokes_fe, 
                            quadrature_formula, 
                            update_values); 
    std::vector<Tensor<1, dim>> velocity_values(n_q_points); 

    const FEValuesExtractors::Vector velocities(0); 

    double max_local_velocity = 0; 

    for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          fe_values[velocities].get_function_values(stokes_solution, 
                                                    velocity_values); 

          for (unsigned int q = 0; q < n_q_points; ++q) 
            max_local_velocity = 
              std::max(max_local_velocity, velocity_values[q].norm()); 
        } 

    return Utilities::MPI::max(max_local_velocity, MPI_COMM_WORLD); 
  } 
// @sect5{BoussinesqFlowProblem::get_cfl_number}  

// 下一个函数做了类似的事情，但我们现在计算CFL数，即一个单元上的最大速度除以单元直径。这个数字对于确定时间步长是必要的，因为我们对温度方程使用半显式的时间步长方案（讨论见 step-31 ）。我们用上述同样的方法计算它。在所有本地拥有的单元上计算本地最大值，然后通过MPI交换，找到全球最大值。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::get_cfl_number() const 
  { 
    const QIterated<dim> quadrature_formula(QTrapezoid<1>(), 
                                            parameters.stokes_velocity_degree); 
    const unsigned int   n_q_points = quadrature_formula.size(); 

    FEValues<dim>               fe_values(mapping, 
                            stokes_fe, 
                            quadrature_formula, 
                            update_values); 
    std::vector<Tensor<1, dim>> velocity_values(n_q_points); 

    const FEValuesExtractors::Vector velocities(0); 

    double max_local_cfl = 0; 

    for (const auto &cell : stokes_dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          fe_values[velocities].get_function_values(stokes_solution, 
                                                    velocity_values); 

          double max_local_velocity = 1e-10; 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            max_local_velocity = 
              std::max(max_local_velocity, velocity_values[q].norm()); 
          max_local_cfl = 
            std::max(max_local_cfl, max_local_velocity / cell->diameter()); 
        } 

    return Utilities::MPI::max(max_local_cfl, MPI_COMM_WORLD); 
  } 
// @sect5{BoussinesqFlowProblem::get_entropy_variation}  

// 接下来是计算全局熵的变化 $\|E(T)-\bar{E}(T)\|_\infty$ ，其中熵 $E$ 的定义如介绍中所讨论的。 这对于评估温度方程中的稳定度是必要的，正如介绍中所解释的。实际上，只有当我们在残差计算中使用 $\alpha=2$ 作为幂时，才需要熵的变化。无限准则是由正交点上的最大值计算出来的，就像离散计算中通常的那样。

// 为了计算这个量，我们首先要找到空间平均数 $\bar{E}(T)$ ，然后评估最大值。然而，这意味着我们需要执行两个循环。我们可以通过注意到 $\|E(T)-\bar{E}(T)\|_\infty =
//  \max\big(E_{\textrm{max}}(T)-\bar{E}(T),
//  \bar{E}(T)-E_{\textrm{min}}(T)\big)$ ，即正负方向上与平均熵的偏差的最大值来避免开销。我们在后一个公式中需要的四个量（最大熵、最小熵、平均熵、面积）都可以在所有单元格的同一个循环中进行评估，所以我们选择这个更简单的变体。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::get_entropy_variation( 
    const double average_temperature) const 
  { 
    if (parameters.stabilization_alpha != 2) 
      return 1.; 

    const QGauss<dim>  quadrature_formula(parameters.temperature_degree + 1); 
    const unsigned int n_q_points = quadrature_formula.size(); 

    FEValues<dim>       fe_values(temperature_fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values); 
    std::vector<double> old_temperature_values(n_q_points); 
    std::vector<double> old_old_temperature_values(n_q_points); 

// 在上面的两个函数中，我们计算了所有非负数的最大值，所以我们知道0肯定是一个下限。另一方面，在这里我们需要找到与平均值的最大偏差，也就是说，我们需要知道熵的最大和最小值，而这些值的符号我们并不事先知道。

// 为了计算它，我们可以从我们可以存储在一个双精度数字中的最大和最小的可能值开始。最小值被初始化为一个更大的数字，最大值被初始化为一个比将要出现的任何一个数字都小的数字。然后，我们保证这些数字将在第一个单元的循环中被覆盖，或者，如果这个处理器不拥有任何单元，最迟在通信步骤中被覆盖。下面的循环将计算最小和最大的局部熵，并跟踪我们局部拥有的域的面积/体积，以及对其熵的积分。

    double min_entropy = std::numeric_limits<double>::max(), 
           max_entropy = -std::numeric_limits<double>::max(), area = 0, 
           entropy_integrated = 0; 

    for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          fe_values.get_function_values(old_temperature_solution, 
                                        old_temperature_values); 
          fe_values.get_function_values(old_old_temperature_solution, 
                                        old_old_temperature_values); 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 
              const double T = 
                (old_temperature_values[q] + old_old_temperature_values[q]) / 2; 
              const double entropy = 
                ((T - average_temperature) * (T - average_temperature)); 

              min_entropy = std::min(min_entropy, entropy); 
              max_entropy = std::max(max_entropy, entropy); 
              area += fe_values.JxW(q); 
              entropy_integrated += fe_values.JxW(q) * entropy; 
            } 
        } 

// 现在我们只需要在处理器之间交换数据：我们需要将两个积分相加（  <code>area</code>, <code>entropy_integrated</code>  ），并得到最大和最小的极值。我们可以通过四个不同的数据交换来完成这个任务，但我们只需要两个就可以了。  Utilities::MPI::sum 也有一个变体，它接受一个数组的值，这些值都是要加起来的。我们还可以利用 Utilities::MPI::max 函数，认识到在最小熵上形成最小值等于在最小熵的负值上形成最大值的负值；然后这个最大值可以与在最大熵上形成最大值结合起来。

    const double local_sums[2]   = {entropy_integrated, area}, 
                 local_maxima[2] = {-min_entropy, max_entropy}; 
    double global_sums[2], global_maxima[2]; 

    Utilities::MPI::sum(local_sums, MPI_COMM_WORLD, global_sums); 
    Utilities::MPI::max(local_maxima, MPI_COMM_WORLD, global_maxima); 

// 以这种方式计算了所有的东西之后，我们就可以计算平均熵，并通过取最大值或最小值与平均值的偏差中的较大值来找到 $L^\infty$ 准则。

    const double average_entropy = global_sums[0] / global_sums[1]; 
    const double entropy_diff    = std::max(global_maxima[1] - average_entropy, 
                                         average_entropy - (-global_maxima[0])); 
    return entropy_diff; 
  } 

//  @sect5{BoussinesqFlowProblem::get_extrapolated_temperature_range}  

// 下一个函数是计算整个领域内外推温度的最小值和最大值。同样，这只是  step-31  中相应函数的一个略微修改的版本。和上面的函数一样，我们收集局部最小值和最大值，然后用上面的技巧计算全局极值。

// 正如在 step-31 中已经讨论过的，该函数需要区分第一个时间步长和所有后续时间步长，因为当至少有两个以前的时间步长时，它使用了一个高阶温度外推方案。

  template <int dim> 
  std::pair<double, double> 
  BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range() const 
  { 
    const QIterated<dim> quadrature_formula(QTrapezoid<1>(), 
                                            parameters.temperature_degree); 
    const unsigned int   n_q_points = quadrature_formula.size(); 

    FEValues<dim>       fe_values(mapping, 
                            temperature_fe, 
                            quadrature_formula, 
                            update_values); 
    std::vector<double> old_temperature_values(n_q_points); 
    std::vector<double> old_old_temperature_values(n_q_points); 

    double min_local_temperature = std::numeric_limits<double>::max(), 
           max_local_temperature = -std::numeric_limits<double>::max(); 

    if (timestep_number != 0) 
      { 
        for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
          if (cell->is_locally_owned()) 
            { 
              fe_values.reinit(cell); 
              fe_values.get_function_values(old_temperature_solution, 
                                            old_temperature_values); 
              fe_values.get_function_values(old_old_temperature_solution, 
                                            old_old_temperature_values); 

              for (unsigned int q = 0; q < n_q_points; ++q) 
                { 
                  const double temperature = 
                    (1. + time_step / old_time_step) * 
                      old_temperature_values[q] - 
                    time_step / old_time_step * old_old_temperature_values[q]; 

                  min_local_temperature = 
                    std::min(min_local_temperature, temperature); 
                  max_local_temperature = 
                    std::max(max_local_temperature, temperature); 
                } 
            } 
      } 
    else 
      { 
        for (const auto &cell : temperature_dof_handler.active_cell_iterators()) 
          if (cell->is_locally_owned()) 
            { 
              fe_values.reinit(cell); 
              fe_values.get_function_values(old_temperature_solution, 
                                            old_temperature_values); 

              for (unsigned int q = 0; q < n_q_points; ++q) 
                { 
                  const double temperature = old_temperature_values[q]; 

                  min_local_temperature = 
                    std::min(min_local_temperature, temperature); 
                  max_local_temperature = 
                    std::max(max_local_temperature, temperature); 
                } 
            } 
      } 

    double local_extrema[2] = {-min_local_temperature, max_local_temperature}; 
    double global_extrema[2]; 
    Utilities::MPI::max(local_extrema, MPI_COMM_WORLD, global_extrema); 

    return std::make_pair(-global_extrema[0], global_extrema[1]); 
  } 
// @sect5{BoussinesqFlowProblem::compute_viscosity}  

// 计算粘度的函数是纯粹的本地函数，所以根本不需要通信。它与 step-31 中的内容基本相同，但如果选择 $\alpha=2$ ，则会有一个最新的粘度表述。

  template <int dim> 
  double BoussinesqFlowProblem<dim>::compute_viscosity( 
    const std::vector<double> &                 old_temperature, 
    const std::vector<double> &                 old_old_temperature, 
    const std::vector<Tensor<1, dim>> &         old_temperature_grads, 
    const std::vector<Tensor<1, dim>> &         old_old_temperature_grads, 
    const std::vector<double> &                 old_temperature_laplacians, 
    const std::vector<double> &                 old_old_temperature_laplacians, 
    const std::vector<Tensor<1, dim>> &         old_velocity_values, 
    const std::vector<Tensor<1, dim>> &         old_old_velocity_values, 
    const std::vector<SymmetricTensor<2, dim>> &old_strain_rates, 
    const std::vector<SymmetricTensor<2, dim>> &old_old_strain_rates, 
    const double                                global_u_infty, 
    const double                                global_T_variation, 
    const double                                average_temperature, 
    const double                                global_entropy_variation, 
    const double                                cell_diameter) const 
  { 
    if (global_u_infty == 0) 
      return 5e-3 * cell_diameter; 

    const unsigned int n_q_points = old_temperature.size(); 

    double max_residual = 0; 
    double max_velocity = 0; 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const Tensor<1, dim> u = 
          (old_velocity_values[q] + old_old_velocity_values[q]) / 2; 

        const SymmetricTensor<2, dim> strain_rate = 
          (old_strain_rates[q] + old_old_strain_rates[q]) / 2; 

        const double T = (old_temperature[q] + old_old_temperature[q]) / 2; 
        const double dT_dt = 
          (old_temperature[q] - old_old_temperature[q]) / old_time_step; 
        const double u_grad_T = 
          u * (old_temperature_grads[q] + old_old_temperature_grads[q]) / 2; 

        const double kappa_Delta_T = 
          EquationData::kappa * 
          (old_temperature_laplacians[q] + old_old_temperature_laplacians[q]) / 
          2; 
        const double gamma = 
          ((EquationData::radiogenic_heating * EquationData::density(T) + 
            2 * EquationData::eta * strain_rate * strain_rate) / 
           (EquationData::density(T) * EquationData::specific_heat)); 

        double residual = std::abs(dT_dt + u_grad_T - kappa_Delta_T - gamma); 
        if (parameters.stabilization_alpha == 2) 
          residual *= std::abs(T - average_temperature); 

        max_residual = std::max(residual, max_residual); 
        max_velocity = std::max(std::sqrt(u * u), max_velocity); 
      } 

    const double max_viscosity = 
      (parameters.stabilization_beta * max_velocity * cell_diameter); 
    if (timestep_number == 0) 
      return max_viscosity; 
    else 
      { 
        Assert(old_time_step > 0, ExcInternalError()); 

        double entropy_viscosity; 
        if (parameters.stabilization_alpha == 2) 
          entropy_viscosity = 
            (parameters.stabilization_c_R * cell_diameter * cell_diameter * 
             max_residual / global_entropy_variation); 
        else 
          entropy_viscosity = 
            (parameters.stabilization_c_R * cell_diameter * 
             global_Omega_diameter * max_velocity * max_residual / 
             (global_u_infty * global_T_variation)); 

        return std::min(max_viscosity, entropy_viscosity); 
      } 
  } 

//  @sect4{The BoussinesqFlowProblem setup functions}  

// 以下三个函数设置了斯托克斯矩阵、用于斯托克斯预调节器的矩阵和温度矩阵。这些代码与 step-31 中的代码基本相同，但为了简单起见，将其分成了三个自己的函数。

// 这里的代码与 step-31 中的代码在功能上的主要区别是，我们要建立的矩阵是分布在多个处理器上的。由于我们仍然希望出于效率的原因先建立起稀疏性模式，我们可以继续将<i>entire</i>的稀疏性模式作为BlockDynamicSparsityPattern来建立，正如我们在 step-31 中所做的那样。然而，这将是低效的：每个处理器将建立相同的稀疏性模式，但只用它初始化矩阵的一小部分。这也违反了一个原则，即每个处理器应该只对它所拥有的单元格（如果有必要的话，还有它周围的幽灵单元格层）工作。

// 相反，我们使用一个类型为 TrilinosWrappers::BlockSparsityPattern, 的对象，它（显然）是对Trilinos提供的稀疏模式对象的一个封装。这样做的好处是Trilinos稀疏模式类可以在多个处理器之间进行通信：如果这个处理器填入它所拥有的单元格产生的所有非零条目，并且其他每个处理器也这样做，那么在由 <code>compress()</code> 调用发起的MPI通信结束后，我们将有全局组装的稀疏模式可用，全局矩阵可以被初始化。

// 在并行初始化Trilinos稀疏度模式时，有一个重要的方面。除了通过 @p stokes_partitioning 索引集指定矩阵的本地拥有的行和列之外，我们还提供了在某个处理器上装配时可能要写进的所有行的信息。本地相关行的集合包含了所有这样的行（可能还有一些不必要的行，但在实际获得所有单元格的索引和解决约束之前，很难找到确切的行索引）。这种额外的信息可以准确地确定在装配过程中发现的非处理器数据的结构。虽然Trilinos矩阵也能在飞行中收集这些信息（当从其他一些reinit方法初始化它们时），但效率较低，在用多线程组装矩阵时，会导致问题。在这个程序中，我们悲观地假设每次只有一个处理器可以在组装时写入矩阵（而计算是并行的），这对特里诺斯矩阵是没有问题的。在实践中，可以通过在不共享顶点的单元中提示WorkStream来做得更好，允许这些单元之间的并行性（参见图形着色算法和带有彩色迭代器的WorkStream参数）。然而，这只在只有一个MPI处理器的情况下有效，因为Trilinos的内部数据结构在飞行中积累非处理器的数据，不是线程安全。有了这里介绍的初始化，就不存在这样的问题，人们可以安全地为这个算法引入图形着色。

// 我们唯一需要做的改变是告诉 DoFTools::make_sparsity_pattern() 函数，它只应该在一个单元格子集上工作，即那些 <code>subdomain_id</code> 等于当前处理器数量的单元格，而忽略所有其他单元格。

// 这个策略被复制到以下三个函数中。

// 注意，Trilinos 矩阵存储的信息包含在稀疏模式中，所以一旦矩阵被赋予稀疏结构，我们就可以安全地释放  <code>sp</code>  变量。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::setup_stokes_matrix( 
    const std::vector<IndexSet> &stokes_partitioning, 
    const std::vector<IndexSet> &stokes_relevant_partitioning) 
  { 
    stokes_matrix.clear(); 

    TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, 
                                              stokes_partitioning, 
                                              stokes_relevant_partitioning, 
                                              MPI_COMM_WORLD); 

    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
    for (unsigned int c = 0; c < dim + 1; ++c) 
      for (unsigned int d = 0; d < dim + 1; ++d) 
        if (!((c == dim) && (d == dim))) 
          coupling[c][d] = DoFTools::always; 
        else 
          coupling[c][d] = DoFTools::none; 

    DoFTools::make_sparsity_pattern(stokes_dof_handler, 
                                    coupling, 
                                    sp, 
                                    stokes_constraints, 
                                    false, 
                                    Utilities::MPI::this_mpi_process( 
                                      MPI_COMM_WORLD)); 
    sp.compress(); 

    stokes_matrix.reinit(sp); 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::setup_stokes_preconditioner( 
    const std::vector<IndexSet> &stokes_partitioning, 
    const std::vector<IndexSet> &stokes_relevant_partitioning) 
  { 
    Amg_preconditioner.reset(); 
    Mp_preconditioner.reset(); 

    stokes_preconditioner_matrix.clear(); 

    TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, 
                                              stokes_partitioning, 
                                              stokes_relevant_partitioning, 
                                              MPI_COMM_WORLD); 

    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
    for (unsigned int c = 0; c < dim + 1; ++c) 
      for (unsigned int d = 0; d < dim + 1; ++d) 
        if (c == d) 
          coupling[c][d] = DoFTools::always; 
        else 
          coupling[c][d] = DoFTools::none; 

    DoFTools::make_sparsity_pattern(stokes_dof_handler, 
                                    coupling, 
                                    sp, 
                                    stokes_constraints, 
                                    false, 
                                    Utilities::MPI::this_mpi_process( 
                                      MPI_COMM_WORLD)); 
    sp.compress(); 

    stokes_preconditioner_matrix.reinit(sp); 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::setup_temperature_matrices( 
    const IndexSet &temperature_partitioner, 
    const IndexSet &temperature_relevant_partitioner) 
  { 
    T_preconditioner.reset(); 
    temperature_mass_matrix.clear(); 
    temperature_stiffness_matrix.clear(); 
    temperature_matrix.clear(); 

    TrilinosWrappers::SparsityPattern sp(temperature_partitioner, 
                                         temperature_partitioner, 
                                         temperature_relevant_partitioner, 
                                         MPI_COMM_WORLD); 
    DoFTools::make_sparsity_pattern(temperature_dof_handler, 
                                    sp, 
                                    temperature_constraints, 
                                    false, 
                                    Utilities::MPI::this_mpi_process( 
                                      MPI_COMM_WORLD)); 
    sp.compress(); 

    temperature_matrix.reinit(sp); 
    temperature_mass_matrix.reinit(sp); 
    temperature_stiffness_matrix.reinit(sp); 
  } 

// 设置函数的其余部分（在拆分出上面的三个函数后）主要是处理我们需要做的跨处理器并行化的事情。因为设置所有这些都是程序的一个重要的计算时间支出，所以我们把在这里做的所有事情都放到一个定时器组中，这样我们就可以在程序结束时得到关于这部分时间的总结信息。

// 像往常一样，我们在顶部列举自由度，并按照组件/块进行排序，然后从零号处理器开始将它们的数字写到屏幕上。当 DoFHandler::distributed_dofs() 函数应用于 parallel::distributed::Triangulation 对象时，对自由度的排序是这样的：所有与子域0相关的自由度排在所有与子域1相关的自由度之前，等等。对于斯托克斯部分，这意味着速度和压力会混在一起，但这可以通过再次按块排序来解决；值得注意的是，后一种操作只保留了所有速度和压力的相对顺序，即在速度块内，我们仍然会将所有与子域零相关的速度放在与子域一相关的速度之前，等等。这一点很重要，因为我们把这个矩阵的每一个块都分布在所有的处理器上，并且希望这样做的方式是，每个处理器存储的矩阵部分与它将实际工作的单元上的自由度大致相等。

// 在打印自由度的数字时，注意如果我们使用许多处理器，这些数字将会很大。因此，我们让流在每三个数字之间放一个逗号分隔符。流的状态，使用locale，从这个操作之前保存到之后。虽然有点不透明，但这段代码是有效的，因为默认的locale（我们使用构造函数调用 <code>std::locale("")</code> 得到的）意味着打印数字时，每三位数字都有一个逗号分隔符（即千、百万、亿）。

// 在这个函数以及下面的许多函数中，我们测量了我们在这里花费的时间，并将其收集在一个叫做 "设置dof系统 "的部分，跨函数调用。这是用一个 TimerOutput::Scope 对象完成的，该对象在构建本地变量时，在上述名称为`computing_timer`的部分启动一个定时器；当`timing_section`变量的析构器被调用时，该定时器再次停止。 当然，这要么发生在函数的末尾，要么我们通过`return`语句离开函数，或者在某处抛出异常时--换句话说，只要我们以任何方式离开这个函数。因此，使用这种 "范围 "对象可以确保我们不必手动添加代码，告诉定时器在每个可能离开这个函数的地方停止。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::setup_dofs() 
  { 
    TimerOutput::Scope timing_section(computing_timer, "Setup dof systems"); 

    stokes_dof_handler.distribute_dofs(stokes_fe); 

    std::vector<unsigned int> stokes_sub_blocks(dim + 1, 0); 
    stokes_sub_blocks[dim] = 1; 
    DoFRenumbering::component_wise(stokes_dof_handler, stokes_sub_blocks); 

    temperature_dof_handler.distribute_dofs(temperature_fe); 

    const std::vector<types::global_dof_index> stokes_dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(stokes_dof_handler, stokes_sub_blocks); 

    const unsigned int n_u = stokes_dofs_per_block[0], 
                       n_p = stokes_dofs_per_block[1], 
                       n_T = temperature_dof_handler.n_dofs(); 

    std::locale s = pcout.get_stream().getloc(); 
    pcout.get_stream().imbue(std::locale("")); 
    pcout << "Number of active cells: " << triangulation.n_global_active_cells() 
          << " (on " << triangulation.n_levels() << " levels)" << std::endl 
          << "Number of degrees of freedom: " << n_u + n_p + n_T << " (" << n_u 
          << '+' << n_p << '+' << n_T << ')' << std::endl 
          << std::endl; 
    pcout.get_stream().imbue(s); 

// 在这之后，我们必须设置各种分区器（类型为 <code>IndexSet</code>  ，见介绍），描述每个矩阵或向量的哪些部分将被存储在哪里，然后调用实际设置矩阵的函数，在最后还要调整我们在这个程序中保留的各种向量的大小。

    std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning; 
    IndexSet              temperature_partitioning(n_T), 
      temperature_relevant_partitioning(n_T); 
    IndexSet stokes_relevant_set; 
    { 
      IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs(); 
      stokes_partitioning.push_back(stokes_index_set.get_view(0, n_u)); 
      stokes_partitioning.push_back(stokes_index_set.get_view(n_u, n_u + n_p)); 

      DoFTools::extract_locally_relevant_dofs(stokes_dof_handler, 
                                              stokes_relevant_set); 
      stokes_relevant_partitioning.push_back( 
        stokes_relevant_set.get_view(0, n_u)); 
      stokes_relevant_partitioning.push_back( 
        stokes_relevant_set.get_view(n_u, n_u + n_p)); 

      temperature_partitioning = temperature_dof_handler.locally_owned_dofs(); 
      DoFTools::extract_locally_relevant_dofs( 
        temperature_dof_handler, temperature_relevant_partitioning); 
    } 

// 在这之后，我们可以计算求解向量的约束，包括悬挂节点约束和斯托克斯和温度场的同质和非同质边界值。请注意，和其他一切一样，约束对象不能在每个处理器上都持有<i>all</i>约束。相反，鉴于每个处理器只在其拥有的单元上组装线性系统，因此每个处理器只需要存储那些对正确性实际必要的约束。正如在 @ref distributed_paper "本文 "中所讨论的，我们需要了解的约束集正是所有本地相关自由度的约束集，所以这就是我们用来初始化约束对象的。

    { 
      stokes_constraints.clear(); 
      stokes_constraints.reinit(stokes_relevant_set); 

      DoFTools::make_hanging_node_constraints(stokes_dof_handler, 
                                              stokes_constraints); 

      FEValuesExtractors::Vector velocity_components(0); 
      VectorTools::interpolate_boundary_values( 
        stokes_dof_handler, 
        0, 
        Functions::ZeroFunction<dim>(dim + 1), 
        stokes_constraints, 
        stokes_fe.component_mask(velocity_components)); 

      std::set<types::boundary_id> no_normal_flux_boundaries; 
      no_normal_flux_boundaries.insert(1); 
      VectorTools::compute_no_normal_flux_constraints(stokes_dof_handler, 
                                                      0, 
                                                      no_normal_flux_boundaries, 
                                                      stokes_constraints, 
                                                      mapping); 
      stokes_constraints.close(); 
    } 
    { 
      temperature_constraints.clear(); 
      temperature_constraints.reinit(temperature_relevant_partitioning); 

      DoFTools::make_hanging_node_constraints(temperature_dof_handler, 
                                              temperature_constraints); 
      VectorTools::interpolate_boundary_values( 
        temperature_dof_handler, 
        0, 
        EquationData::TemperatureInitialValues<dim>(), 
        temperature_constraints); 
      VectorTools::interpolate_boundary_values( 
        temperature_dof_handler, 
        1, 
        EquationData::TemperatureInitialValues<dim>(), 
        temperature_constraints); 
      temperature_constraints.close(); 
    } 

// 做完这些，我们就可以将各种矩阵和向量对象初始化到合适的大小。在最后，我们还记录了所有的矩阵和前置条件器必须在下一个时间步长开始时重新计算。注意我们是如何初始化斯托克斯和温度右侧的向量的。这些是可写的向量（最后一个布尔参数设置为 @p true) ），具有正确的一对一的本地拥有元素的分区，但仍被赋予相关的分区，以弄清要立即设置的向量条目。至于矩阵，这允许用多个线程将本地贡献写入向量（总是假设同一向量条目不被多个线程同时访问）。其他向量只允许对单个元素的读取访问，包括鬼魂，但不适合求解器。

    setup_stokes_matrix(stokes_partitioning, stokes_relevant_partitioning); 
    setup_stokes_preconditioner(stokes_partitioning, 
                                stokes_relevant_partitioning); 
    setup_temperature_matrices(temperature_partitioning, 
                               temperature_relevant_partitioning); 

    stokes_rhs.reinit(stokes_partitioning, 
                      stokes_relevant_partitioning, 
                      MPI_COMM_WORLD, 
                      true); 
    stokes_solution.reinit(stokes_relevant_partitioning, MPI_COMM_WORLD); 
    old_stokes_solution.reinit(stokes_solution); 

    temperature_rhs.reinit(temperature_partitioning, 
                           temperature_relevant_partitioning, 
                           MPI_COMM_WORLD, 
                           true); 
    temperature_solution.reinit(temperature_relevant_partitioning, 
                                MPI_COMM_WORLD); 
    old_temperature_solution.reinit(temperature_solution); 
    old_old_temperature_solution.reinit(temperature_solution); 

    rebuild_stokes_matrix              = true; 
    rebuild_stokes_preconditioner      = true; 
    rebuild_temperature_matrices       = true; 
    rebuild_temperature_preconditioner = true; 
  } 

//  @sect4{The BoussinesqFlowProblem assembly functions}  

// 按照介绍和 @ref threads 模块中的讨论，我们将装配功能分成不同的部分。

//  <ul>  
//  <li>  矩阵和右手边的局部计算，给定某个单元作为输入（这些函数被命名为下面的 <code>local_assemble_*</code> ）。换句话说，得出的函数基本上是 step-31 中所有单元格的循环体。然而，请注意，这些函数将本地计算的结果存储在CopyData命名空间的类的变量中。

//  <li>  然后这些对象被交给第二步，将本地数据写入全局数据结构中（这些函数被命名为下面的 <code>copy_local_to_global_*</code> ）。这些函数是相当琐碎的。

//  <li>  然后这两个子函数被用于各自的汇编例程（下面称为 <code>assemble_*</code> ），在那里，一个WorkStream对象被设置并在属于处理器子域的所有单元中运行。   </ul>  
// @sect5{Stokes preconditioner assembly}  

// 让我们从构建斯托克斯预处理的函数开始。考虑到上面的讨论，其中的前两个是非常微不足道的。请特别注意，使用scratch数据对象的主要意义在于，我们希望避免每次访问新单元时在自由空间上分配任何对象。因此，下面的汇编函数只有自动的局部变量，其他的都是通过从头开始的数据对象访问的，在我们开始对所有单元进行循环之前，只分配了一次。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::local_assemble_stokes_preconditioner( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    Assembly::Scratch::StokesPreconditioner<dim> &        scratch, 
    Assembly::CopyData::StokesPreconditioner<dim> &       data) 
  { 
    const unsigned int dofs_per_cell = stokes_fe.n_dofs_per_cell(); 
    const unsigned int n_q_points = 
      scratch.stokes_fe_values.n_quadrature_points; 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

    scratch.stokes_fe_values.reinit(cell); 
    cell->get_dof_indices(data.local_dof_indices); 

    data.local_matrix = 0; 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          { 
            scratch.grad_phi_u[k] = 
              scratch.stokes_fe_values[velocities].gradient(k, q); 
            scratch.phi_p[k] = scratch.stokes_fe_values[pressure].value(k, q); 
          } 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            data.local_matrix(i, j) += 
              (EquationData::eta * 
                 scalar_product(scratch.grad_phi_u[i], scratch.grad_phi_u[j]) + 
               (1. / EquationData::eta) * EquationData::pressure_scaling * 
                 EquationData::pressure_scaling * 
                 (scratch.phi_p[i] * scratch.phi_p[j])) * 
              scratch.stokes_fe_values.JxW(q); 
      } 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::copy_local_to_global_stokes_preconditioner( 
    const Assembly::CopyData::StokesPreconditioner<dim> &data) 
  { 
    stokes_constraints.distribute_local_to_global(data.local_matrix, 
                                                  data.local_dof_indices, 
                                                  stokes_preconditioner_matrix); 
  } 

// 现在是真正把事情放在一起的函数，使用WorkStream函数。   WorkStream::run 需要一个开始和结束迭代器来列举它应该工作的单元格。通常情况下，我们会使用 DoFHandler::begin_active() 和 DoFHandler::end() 来实现这一点，但在这里，我们实际上只想获得事实上由当前处理器拥有的单元格子集。这就是FilteredIterator类发挥作用的地方：你给它一个单元格范围，它提供一个迭代器，只迭代满足某个谓词的单元格子集（谓词是一个参数的函数，要么返回真，要么返回假）。我们在这里使用的谓词是 IteratorFilters::LocallyOwnedCell, ，也就是说，如果单元格为当前处理器所拥有，它就会准确返回真。这样得到的迭代器范围正是我们需要的。

// 有了这个障碍，我们用这组单元格、scratch和copy对象以及两个函数的指针来调用 WorkStream::run 函数：本地装配和copy-local-to-global函数。这些函数需要有非常具体的签名：前者有三个参数，后者有一个参数（关于这些参数的含义，请参见 WorkStream::run 函数的文档）。注意我们是如何使用lambda函数来创建一个满足这一要求的函数对象的。它使用了指定单元格、抓取数据和复制数据的本地装配函数的函数参数，以及期望将数据写入全局矩阵的复制函数的函数参数（也可参见 step-13 的 <code>assemble_linear_system()</code> 函数中的讨论）。另一方面，成员函数的隐含的第2个参数（即该成员函数要操作的对象的 <code>this</code> 指针）是<i>bound</i>到当前函数的 <code>this</code> 指针的，并被捕获。因此， WorkStream::run 函数不需要知道这些函数所操作的对象的任何信息。

// 当WorkStream被执行时，它将为几个单元创建几个第一类的本地装配例程，并让一些可用的处理器对其工作。然而，需要同步的函数，即写进全局矩阵的操作，每次只由一个线程按照规定的顺序执行。当然，这只适用于单个MPI进程上的并行化。不同的MPI进程将有自己的WorkStream对象，并完全独立地进行这项工作（并且在不同的内存空间）。在分布式计算中，一些数据将积累在不属于各自处理器的自由度上。如果每次遇到这样的自由度就把数据送来送去，那就没有效率了。取而代之的是，Trilinos稀疏矩阵将保留这些数据，并在装配结束时通过调用 <code>compress()</code> 命令将其发送给所有者。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner() 
  { 
    stokes_preconditioner_matrix = 0; 

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree + 1); 

    using CellFilter = 
      FilteredIterator<typename DoFHandler<2>::active_cell_iterator>; 

    auto worker = 
      [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
             Assembly::Scratch::StokesPreconditioner<dim> &        scratch, 
             Assembly::CopyData::StokesPreconditioner<dim> &       data) { 
        this->local_assemble_stokes_preconditioner(cell, scratch, data); 
      }; 

    auto copier = 
      [this](const Assembly::CopyData::StokesPreconditioner<dim> &data) { 
        this->copy_local_to_global_stokes_preconditioner(data); 
      }; 

    WorkStream::run(CellFilter(IteratorFilters::LocallyOwnedCell(), 
                               stokes_dof_handler.begin_active()), 
                    CellFilter(IteratorFilters::LocallyOwnedCell(), 
                               stokes_dof_handler.end()), 
                    worker, 
                    copier, 
                    Assembly::Scratch::StokesPreconditioner<dim>( 
                      stokes_fe, 
                      quadrature_formula, 
                      mapping, 
                      update_JxW_values | update_values | update_gradients), 
                    Assembly::CopyData::StokesPreconditioner<dim>(stokes_fe)); 

    stokes_preconditioner_matrix.compress(VectorOperation::add); 
  } 

// 这个模块的最后一个函数启动了斯托克斯预处理矩阵的装配，然后实际上是建立了斯托克斯预处理程序。它与串行情况下的功能基本相同。与 step-31 唯一不同的是，我们对压力质量矩阵使用雅可比预处理，而不是IC，这一点在介绍中已经讨论过。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::build_stokes_preconditioner() 
  { 
    if (rebuild_stokes_preconditioner == false) 
      return; 

    TimerOutput::Scope timer_section(computing_timer, 
                                     "   Build Stokes preconditioner"); 
    pcout << "   Rebuilding Stokes preconditioner..." << std::flush; 

    assemble_stokes_preconditioner(); 

    std::vector<std::vector<bool>> constant_modes; 
    FEValuesExtractors::Vector     velocity_components(0); 
    DoFTools::extract_constant_modes(stokes_dof_handler, 
                                     stokes_fe.component_mask( 
                                       velocity_components), 
                                     constant_modes); 

    Mp_preconditioner = 
      std::make_shared<TrilinosWrappers::PreconditionJacobi>(); 
    Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionAMG>(); 

    TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data; 
    Amg_data.constant_modes        = constant_modes; 
    Amg_data.elliptic              = true; 
    Amg_data.higher_order_elements = true; 
    Amg_data.smoother_sweeps       = 2; 
    Amg_data.aggregation_threshold = 0.02; 

    Mp_preconditioner->initialize(stokes_preconditioner_matrix.block(1, 1)); 
    Amg_preconditioner->initialize(stokes_preconditioner_matrix.block(0, 0), 
                                   Amg_data); 

    rebuild_stokes_preconditioner = false; 

    pcout << std::endl; 
  } 
// @sect5{Stokes system assembly}  

// 接下来的三个函数实现了斯托克斯系统的装配，同样分为执行局部计算的部分，将局部数据写入全局矩阵和向量的部分，以及在WorkStream类的帮助下实际运行所有单元的循环。请注意，只有在我们改变了网格的情况下才需要进行斯托克斯矩阵的组装。否则，这里只需要计算（与温度有关的）右手边。由于我们正在处理分布式矩阵和向量，我们必须在装配结束时调用相应的 <code>compress()</code> 函数，以便将非本地数据发送到所有者进程。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::local_assemble_stokes_system( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    Assembly::Scratch::StokesSystem<dim> &                scratch, 
    Assembly::CopyData::StokesSystem<dim> &               data) 
  { 
    const unsigned int dofs_per_cell = 
      scratch.stokes_fe_values.get_fe().n_dofs_per_cell(); 
    const unsigned int n_q_points = 
      scratch.stokes_fe_values.n_quadrature_points; 

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

    scratch.stokes_fe_values.reinit(cell); 

    typename DoFHandler<dim>::active_cell_iterator temperature_cell( 
      &triangulation, cell->level(), cell->index(), &temperature_dof_handler); 
    scratch.temperature_fe_values.reinit(temperature_cell); 

    if (rebuild_stokes_matrix) 
      data.local_matrix = 0; 
    data.local_rhs = 0; 

    scratch.temperature_fe_values.get_function_values( 
      old_temperature_solution, scratch.old_temperature_values); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        const double old_temperature = scratch.old_temperature_values[q]; 

        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          { 
            scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value(k, q); 
            if (rebuild_stokes_matrix) 
              { 
                scratch.grads_phi_u[k] = 
                  scratch.stokes_fe_values[velocities].symmetric_gradient(k, q); 
                scratch.div_phi_u[k] = 
                  scratch.stokes_fe_values[velocities].divergence(k, q); 
                scratch.phi_p[k] = 
                  scratch.stokes_fe_values[pressure].value(k, q); 
              } 
          } 

        if (rebuild_stokes_matrix == true) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              data.local_matrix(i, j) += 
                (EquationData::eta * 2 * 
                   (scratch.grads_phi_u[i] * scratch.grads_phi_u[j]) - 
                 (EquationData::pressure_scaling * scratch.div_phi_u[i] * 
                  scratch.phi_p[j]) - 
                 (EquationData::pressure_scaling * scratch.phi_p[i] * 
                  scratch.div_phi_u[j])) * 
                scratch.stokes_fe_values.JxW(q); 

        const Tensor<1, dim> gravity = EquationData::gravity_vector( 
          scratch.stokes_fe_values.quadrature_point(q)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          data.local_rhs(i) += (EquationData::density(old_temperature) * 
                                gravity * scratch.phi_u[i]) * 
                               scratch.stokes_fe_values.JxW(q); 
      } 

    cell->get_dof_indices(data.local_dof_indices); 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::copy_local_to_global_stokes_system( 
    const Assembly::CopyData::StokesSystem<dim> &data) 
  { 
    if (rebuild_stokes_matrix == true) 
      stokes_constraints.distribute_local_to_global(data.local_matrix, 
                                                    data.local_rhs, 
                                                    data.local_dof_indices, 
                                                    stokes_matrix, 
                                                    stokes_rhs); 
    else 
      stokes_constraints.distribute_local_to_global(data.local_rhs, 
                                                    data.local_dof_indices, 
                                                    stokes_rhs); 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_stokes_system() 
  { 
    TimerOutput::Scope timer_section(computing_timer, 
                                     "   Assemble Stokes system"); 

    if (rebuild_stokes_matrix == true) 
      stokes_matrix = 0; 

    stokes_rhs = 0; 

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree + 1); 

    using CellFilter = 
      FilteredIterator<typename DoFHandler<2>::active_cell_iterator>; 

    WorkStream::run( 
      CellFilter(IteratorFilters::LocallyOwnedCell(), 
                 stokes_dof_handler.begin_active()), 
      CellFilter(IteratorFilters::LocallyOwnedCell(), stokes_dof_handler.end()), 
      [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
             Assembly::Scratch::StokesSystem<dim> &                scratch, 
             Assembly::CopyData::StokesSystem<dim> &               data) { 
        this->local_assemble_stokes_system(cell, scratch, data); 
      }, 
      [this](const Assembly::CopyData::StokesSystem<dim> &data) { 
        this->copy_local_to_global_stokes_system(data); 
      }, 
      Assembly::Scratch::StokesSystem<dim>( 
        stokes_fe, 
        mapping, 
        quadrature_formula, 
        (update_values | update_quadrature_points | update_JxW_values | 
         (rebuild_stokes_matrix == true ? update_gradients : UpdateFlags(0))), 
        temperature_fe, 
        update_values), 
      Assembly::CopyData::StokesSystem<dim>(stokes_fe)); 

    if (rebuild_stokes_matrix == true) 
      stokes_matrix.compress(VectorOperation::add); 
    stokes_rhs.compress(VectorOperation::add); 

    rebuild_stokes_matrix = false; 

    pcout << std::endl; 
  } 
// @sect5{Temperature matrix assembly}  

// 下面三个函数要完成的任务是计算温度系统的质量矩阵和拉普拉斯矩阵。这些将被结合起来，以产生半隐式时间步进矩阵，该矩阵由质量矩阵加上一个与时间 step- 相关的权重系数乘以拉普拉斯矩阵组成。这个函数本质上还是从 step-31 开始的所有单元的循环主体。

// 下面两个函数的功能与上面的类似。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::local_assemble_temperature_matrix( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    Assembly::Scratch::TemperatureMatrix<dim> &           scratch, 
    Assembly::CopyData::TemperatureMatrix<dim> &          data) 
  { 
    const unsigned int dofs_per_cell = 
      scratch.temperature_fe_values.get_fe().n_dofs_per_cell(); 
    const unsigned int n_q_points = 
      scratch.temperature_fe_values.n_quadrature_points; 

    scratch.temperature_fe_values.reinit(cell); 
    cell->get_dof_indices(data.local_dof_indices); 

    data.local_mass_matrix      = 0; 
    data.local_stiffness_matrix = 0; 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          { 
            scratch.grad_phi_T[k] = 
              scratch.temperature_fe_values.shape_grad(k, q); 
            scratch.phi_T[k] = scratch.temperature_fe_values.shape_value(k, q); 
          } 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            { 
              data.local_mass_matrix(i, j) += 
                (scratch.phi_T[i] * scratch.phi_T[j] * 
                 scratch.temperature_fe_values.JxW(q)); 
              data.local_stiffness_matrix(i, j) += 
                (EquationData::kappa * scratch.grad_phi_T[i] * 
                 scratch.grad_phi_T[j] * scratch.temperature_fe_values.JxW(q)); 
            } 
      } 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::copy_local_to_global_temperature_matrix( 
    const Assembly::CopyData::TemperatureMatrix<dim> &data) 
  { 
    temperature_constraints.distribute_local_to_global(data.local_mass_matrix, 
                                                       data.local_dof_indices, 
                                                       temperature_mass_matrix); 
    temperature_constraints.distribute_local_to_global( 
      data.local_stiffness_matrix, 
      data.local_dof_indices, 
      temperature_stiffness_matrix); 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::assemble_temperature_matrix() 
  { 
    if (rebuild_temperature_matrices == false) 
      return; 

    TimerOutput::Scope timer_section(computing_timer, 
                                     "   Assemble temperature matrices"); 
    temperature_mass_matrix      = 0; 
    temperature_stiffness_matrix = 0; 

    const QGauss<dim> quadrature_formula(parameters.temperature_degree + 2); 

    using CellFilter = 
      FilteredIterator<typename DoFHandler<2>::active_cell_iterator>; 

    WorkStream::run( 
      CellFilter(IteratorFilters::LocallyOwnedCell(), 
                 temperature_dof_handler.begin_active()), 
      CellFilter(IteratorFilters::LocallyOwnedCell(), 
                 temperature_dof_handler.end()), 
      [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
             Assembly::Scratch::TemperatureMatrix<dim> &           scratch, 
             Assembly::CopyData::TemperatureMatrix<dim> &          data) { 
        this->local_assemble_temperature_matrix(cell, scratch, data); 
      }, 
      [this](const Assembly::CopyData::TemperatureMatrix<dim> &data) { 
        this->copy_local_to_global_temperature_matrix(data); 
      }, 
      Assembly::Scratch::TemperatureMatrix<dim>(temperature_fe, 
                                                mapping, 
                                                quadrature_formula), 
      Assembly::CopyData::TemperatureMatrix<dim>(temperature_fe)); 

    temperature_mass_matrix.compress(VectorOperation::add); 
    temperature_stiffness_matrix.compress(VectorOperation::add); 

    rebuild_temperature_matrices       = false; 
    rebuild_temperature_preconditioner = true; 
  } 
// @sect5{Temperature right hand side assembly}  

// 这是最后一个装配函数。它计算温度系统的右侧，其中包括对流和稳定项。它包括对正交点上的旧解的大量评估（这对于计算稳定化的人工粘性是必要的），但在其他方面与其他装配函数类似。请注意，我们再次解决了具有不均匀边界条件的困境，只是在这一点上做了一个右手边（比较上面对 <code>project()</code> 函数的评论）。我们创建一些矩阵列，其值正好是为温度刚度矩阵输入的值，如果我们有不均匀约束的DFS的话。这将说明右边的向量与温度矩阵系统的正确平衡。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::local_assemble_temperature_rhs( 
    const std::pair<double, double> global_T_range, 
    const double                    global_max_velocity, 
    const double                    global_entropy_variation, 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    Assembly::Scratch::TemperatureRHS<dim> &              scratch, 
    Assembly::CopyData::TemperatureRHS<dim> &             data) 
  { 
    const bool use_bdf2_scheme = (timestep_number != 0); 

    const unsigned int dofs_per_cell = 
      scratch.temperature_fe_values.get_fe().n_dofs_per_cell(); 
    const unsigned int n_q_points = 
      scratch.temperature_fe_values.n_quadrature_points; 

    const FEValuesExtractors::Vector velocities(0); 

    data.local_rhs     = 0; 
    data.matrix_for_bc = 0; 
    cell->get_dof_indices(data.local_dof_indices); 

    scratch.temperature_fe_values.reinit(cell); 

    typename DoFHandler<dim>::active_cell_iterator stokes_cell( 
      &triangulation, cell->level(), cell->index(), &stokes_dof_handler); 
    scratch.stokes_fe_values.reinit(stokes_cell); 

    scratch.temperature_fe_values.get_function_values( 
      old_temperature_solution, scratch.old_temperature_values); 
    scratch.temperature_fe_values.get_function_values( 
      old_old_temperature_solution, scratch.old_old_temperature_values); 

    scratch.temperature_fe_values.get_function_gradients( 
      old_temperature_solution, scratch.old_temperature_grads); 
    scratch.temperature_fe_values.get_function_gradients( 
      old_old_temperature_solution, scratch.old_old_temperature_grads); 

    scratch.temperature_fe_values.get_function_laplacians( 
      old_temperature_solution, scratch.old_temperature_laplacians); 
    scratch.temperature_fe_values.get_function_laplacians( 
      old_old_temperature_solution, scratch.old_old_temperature_laplacians); 

    scratch.stokes_fe_values[velocities].get_function_values( 
      stokes_solution, scratch.old_velocity_values); 
    scratch.stokes_fe_values[velocities].get_function_values( 
      old_stokes_solution, scratch.old_old_velocity_values); 
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients( 
      stokes_solution, scratch.old_strain_rates); 
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients( 
      old_stokes_solution, scratch.old_old_strain_rates); 

    const double nu = 
      compute_viscosity(scratch.old_temperature_values, 
                        scratch.old_old_temperature_values, 
                        scratch.old_temperature_grads, 
                        scratch.old_old_temperature_grads, 
                        scratch.old_temperature_laplacians, 
                        scratch.old_old_temperature_laplacians, 
                        scratch.old_velocity_values, 
                        scratch.old_old_velocity_values, 
                        scratch.old_strain_rates, 
                        scratch.old_old_strain_rates, 
                        global_max_velocity, 
                        global_T_range.second - global_T_range.first, 
                        0.5 * (global_T_range.second + global_T_range.first), 
                        global_entropy_variation, 
                        cell->diameter()); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          { 
            scratch.phi_T[k] = scratch.temperature_fe_values.shape_value(k, q); 
            scratch.grad_phi_T[k] = 
              scratch.temperature_fe_values.shape_grad(k, q); 
          } 

        const double T_term_for_rhs = 
          (use_bdf2_scheme ? 
             (scratch.old_temperature_values[q] * 
                (1 + time_step / old_time_step) - 
              scratch.old_old_temperature_values[q] * (time_step * time_step) / 
                (old_time_step * (time_step + old_time_step))) : 
             scratch.old_temperature_values[q]); 

        const double ext_T = 
          (use_bdf2_scheme ? (scratch.old_temperature_values[q] * 
                                (1 + time_step / old_time_step) - 
                              scratch.old_old_temperature_values[q] * 
                                time_step / old_time_step) : 
                             scratch.old_temperature_values[q]); 

        const Tensor<1, dim> ext_grad_T = 
          (use_bdf2_scheme ? (scratch.old_temperature_grads[q] * 
                                (1 + time_step / old_time_step) - 
                              scratch.old_old_temperature_grads[q] * time_step / 
                                old_time_step) : 
                             scratch.old_temperature_grads[q]); 

        const Tensor<1, dim> extrapolated_u = 
          (use_bdf2_scheme ? 
             (scratch.old_velocity_values[q] * (1 + time_step / old_time_step) - 
              scratch.old_old_velocity_values[q] * time_step / old_time_step) : 
             scratch.old_velocity_values[q]); 

        const SymmetricTensor<2, dim> extrapolated_strain_rate = 
          (use_bdf2_scheme ? 
             (scratch.old_strain_rates[q] * (1 + time_step / old_time_step) - 
              scratch.old_old_strain_rates[q] * time_step / old_time_step) : 
             scratch.old_strain_rates[q]); 

        const double gamma = 
          ((EquationData::radiogenic_heating * EquationData::density(ext_T) + 
            2 * EquationData::eta * extrapolated_strain_rate * 
              extrapolated_strain_rate) / 
           (EquationData::density(ext_T) * EquationData::specific_heat)); 

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          { 
            data.local_rhs(i) += 
              (T_term_for_rhs * scratch.phi_T[i] - 
               time_step * extrapolated_u * ext_grad_T * scratch.phi_T[i] - 
               time_step * nu * ext_grad_T * scratch.grad_phi_T[i] + 
               time_step * gamma * scratch.phi_T[i]) * 
              scratch.temperature_fe_values.JxW(q); 

            if (temperature_constraints.is_inhomogeneously_constrained( 
                  data.local_dof_indices[i])) 
              { 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  data.matrix_for_bc(j, i) += 
                    (scratch.phi_T[i] * scratch.phi_T[j] * 
                       (use_bdf2_scheme ? ((2 * time_step + old_time_step) / 
                                           (time_step + old_time_step)) : 
                                          1.) + 
                     scratch.grad_phi_T[i] * scratch.grad_phi_T[j] * 
                       EquationData::kappa * time_step) * 
                    scratch.temperature_fe_values.JxW(q); 
              } 
          } 
      } 
  } 

  template <int dim> 
  void BoussinesqFlowProblem<dim>::copy_local_to_global_temperature_rhs( 
    const Assembly::CopyData::TemperatureRHS<dim> &data) 
  { 
    temperature_constraints.distribute_local_to_global(data.local_rhs, 
                                                       data.local_dof_indices, 
                                                       temperature_rhs, 
                                                       data.matrix_for_bc); 
  } 

// 在运行实际计算右手边的WorkStream的函数中，我们也生成了最终矩阵。如上所述，它是质量矩阵和拉普拉斯矩阵的总和，再加上一些与时间 step- 相关的权重。这个权重是由BDF-2时间积分方案指定的，见  step-31  中的介绍。这个教程程序的新内容（除了使用MPI并行化和WorkStream类），是我们现在也预先计算了温度预处理程序。原因是与求解器相比，设置雅可比预处理器需要明显的时间，因为我们通常只需要10到20次迭代来求解温度系统（这听起来很奇怪，因为雅可比实际上只包括对角线，但在特里诺斯，它是从更普遍的点松弛预处理器框架中衍生出来的，效率有点低）。因此，尽管由于时间步长可能会发生变化，矩阵条目可能会略有变化，但预先计算预处理程序的效率更高。这不是太大的问题，因为我们每隔几步就重新网格化（然后重新生成预处理程序）。

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

    if (rebuild_temperature_preconditioner == true) 
      { 
        T_preconditioner = 
          std::make_shared<TrilinosWrappers::PreconditionJacobi>(); 
        T_preconditioner->initialize(temperature_matrix); 
        rebuild_temperature_preconditioner = false; 
      } 

// 接下来的部分是计算右手边的向量。 为此，我们首先计算平均温度  $T_m$  ，我们通过残差  $E(T) =
//  (T-T_m)^2$  来评估人工黏度的稳定。我们通过在熵粘度的定义中把最高和最低温度之间的中点定义为平均温度来做到这一点。另一种方法是使用积分平均，但结果对这种选择不是很敏感。那么剩下的就只需要再次调用 WorkStream::run ，将每次调用都相同的 <code>local_assemble_temperature_rhs</code> 函数的参数绑定到正确的值中。

    temperature_rhs = 0; 

    const QGauss<dim> quadrature_formula(parameters.temperature_degree + 2); 
    const std::pair<double, double> global_T_range = 
      get_extrapolated_temperature_range(); 

    const double average_temperature = 
      0.5 * (global_T_range.first + global_T_range.second); 
    const double global_entropy_variation = 
      get_entropy_variation(average_temperature); 

    using CellFilter = 
      FilteredIterator<typename DoFHandler<2>::active_cell_iterator>; 

    auto worker = 
      [this, global_T_range, maximal_velocity, global_entropy_variation]( 
        const typename DoFHandler<dim>::active_cell_iterator &cell, 
        Assembly::Scratch::TemperatureRHS<dim> &              scratch, 
        Assembly::CopyData::TemperatureRHS<dim> &             data) { 
        this->local_assemble_temperature_rhs(global_T_range, 
                                             maximal_velocity, 
                                             global_entropy_variation, 
                                             cell, 
                                             scratch, 
                                             data); 
      }; 

    auto copier = [this](const Assembly::CopyData::TemperatureRHS<dim> &data) { 
      this->copy_local_to_global_temperature_rhs(data); 
    }; 

    WorkStream::run(CellFilter(IteratorFilters::LocallyOwnedCell(), 
                               temperature_dof_handler.begin_active()), 
                    CellFilter(IteratorFilters::LocallyOwnedCell(), 
                               temperature_dof_handler.end()), 
                    worker, 
                    copier, 
                    Assembly::Scratch::TemperatureRHS<dim>( 
                      temperature_fe, stokes_fe, mapping, quadrature_formula), 
                    Assembly::CopyData::TemperatureRHS<dim>(temperature_fe)); 

    temperature_rhs.compress(VectorOperation::add); 
  } 

//  @sect4{BoussinesqFlowProblem::solve}  

// 这个函数在Boussinesq问题的每个时间步长中求解线性系统。首先，我们在斯托克斯系统上工作，然后在温度系统上工作。从本质上讲，它与  step-31  中的相应函数做了同样的事情。然而，这里有一些变化。

// 第一个变化与我们存储解决方案的方式有关：我们在每个MPI节点上保留具有本地拥有的自由度的向量加鬼魂节点。当我们进入一个应该用分布式矩阵进行矩阵-向量乘积的求解器时，这不是合适的形式，虽然。在那里，我们希望求解向量的分布方式与矩阵的分布方式相同，即没有任何重影。所以我们首先要做的是生成一个名为 <code>distributed_stokes_solution</code> 的分布式向量，并只将本地拥有的dof放入其中，这可以通过特里诺向量的 <code>operator=</code> 整齐地完成。

// 接下来，我们为求解器缩放压力解（或者说，初始猜测），使其与矩阵中的长度尺度相匹配，正如在介绍中讨论的那样。在求解完成后，我们也会立即将压力值缩回到正确的单位。 我们还需要将悬挂节点的压力值设置为零。这一点我们在 step-31 中也做过，以避免一些在求解阶段实际上无关紧要的向量项干扰舒尔补数。与 step-31 不同的是，这里我们只对局部拥有的压力道夫进行了处理。在对斯托克斯解进行求解后，每个处理器将分布式解复制到解向量中，其中也包括鬼元素。

// 第三个也是最明显的变化是，我们有两种斯托克斯求解器的变体。一种是有时会崩溃的快速求解器，另一种是速度较慢的稳健求解器。这就是我们在介绍中已经讨论过的。以下是我们如何实现它的。首先，我们用快速求解器进行30次迭代，该求解器是基于AMG V型循环的简单预处理，而不是近似求解（这由 <code>false</code> 对象的 <code>LinearSolvers::BlockSchurPreconditioner</code> 参数表示）。如果我们收敛了，一切都很好。如果我们没有收敛，求解器控制对象将抛出一个异常 SolverControl::NoConvergence. 通常，这将中止程序，因为我们在通常的 <code>solve()</code> 函数中没有捕捉它们。这当然不是我们想在这里发生的。相反，我们希望切换到强求解器，并继续用我们目前得到的任何矢量进行求解。因此，我们用C++的try/catch机制来捕获这个异常。然后我们在 <code>catch</code> 子句中简单地再次经历相同的求解器序列，这次我们将 @p true 标志传递给强求解器的预处理程序，表示近似CG求解。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::solve() 
  { 
    { 
      TimerOutput::Scope timer_section(computing_timer, 
                                       "   Solve Stokes system"); 

      pcout << "   Solving Stokes system... " << std::flush; 

      TrilinosWrappers::MPI::BlockVector distributed_stokes_solution( 
        stokes_rhs); 
      distributed_stokes_solution = stokes_solution; 

      distributed_stokes_solution.block(1) /= EquationData::pressure_scaling; 

      const unsigned int 
        start = (distributed_stokes_solution.block(0).size() + 
                 distributed_stokes_solution.block(1).local_range().first), 
        end   = (distributed_stokes_solution.block(0).size() + 
               distributed_stokes_solution.block(1).local_range().second); 
      for (unsigned int i = start; i < end; ++i) 
        if (stokes_constraints.is_constrained(i)) 
          distributed_stokes_solution(i) = 0; 

      PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem; 

      unsigned int  n_iterations     = 0; 
      const double  solver_tolerance = 1e-8 * stokes_rhs.l2_norm(); 
      SolverControl solver_control(30, solver_tolerance); 

      try 
        { 
          const LinearSolvers::BlockSchurPreconditioner< 
            TrilinosWrappers::PreconditionAMG, 
            TrilinosWrappers::PreconditionJacobi> 
            preconditioner(stokes_matrix, 
                           stokes_preconditioner_matrix, 
                           *Mp_preconditioner, 
                           *Amg_preconditioner, 
                           false); 

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver( 
            solver_control, 
            mem, 
            SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData( 
              30)); 
          solver.solve(stokes_matrix, 
                       distributed_stokes_solution, 
                       stokes_rhs, 
                       preconditioner); 

          n_iterations = solver_control.last_step(); 
        } 

      catch (SolverControl::NoConvergence &) 
        { 
          const LinearSolvers::BlockSchurPreconditioner< 
            TrilinosWrappers::PreconditionAMG, 
            TrilinosWrappers::PreconditionJacobi> 
            preconditioner(stokes_matrix, 
                           stokes_preconditioner_matrix, 
                           *Mp_preconditioner, 
                           *Amg_preconditioner, 
                           true); 

          SolverControl solver_control_refined(stokes_matrix.m(), 
                                               solver_tolerance); 
          SolverFGMRES<TrilinosWrappers::MPI::BlockVector> solver( 
            solver_control_refined, 
            mem, 
            SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData( 
              50)); 
          solver.solve(stokes_matrix, 
                       distributed_stokes_solution, 
                       stokes_rhs, 
                       preconditioner); 

          n_iterations = 
            (solver_control.last_step() + solver_control_refined.last_step()); 
        } 

      stokes_constraints.distribute(distributed_stokes_solution); 

      distributed_stokes_solution.block(1) *= EquationData::pressure_scaling; 

      stokes_solution = distributed_stokes_solution; 
      pcout << n_iterations << " iterations." << std::endl; 
    } 

// 现在让我们转到温度部分。首先，我们计算时间步长。我们发现，对于壳的几何形状，我们需要三维的时间步长比二维的小。这是因为在这种情况下，单元格的变形更大（决定CFL数值的是最小的边长）。我们不是像 step-31 中那样从最大速度和最小网格尺寸计算时间步长，而是计算局部的CFL数，即在每个单元上计算最大速度乘以网格尺寸，并计算它们的最大值。因此，我们需要将时间步长前面的因子选择得稍小一些。

// 在温度的右手边装配后，我们解决温度的线性系统（有完全分布的向量，没有任何鬼魂），应用约束条件，并将向量复制回有鬼魂的向量。

// 最后，我们提取与 step-31 类似的温度范围，以产生一些输出（例如为了帮助我们选择稳定常数，如介绍中所讨论的）。唯一的区别是，我们需要在所有处理器上交换最大值。

    { 
      TimerOutput::Scope timer_section(computing_timer, 
                                       "   Assemble temperature rhs"); 

      old_time_step = time_step; 

      const double scaling = (dim == 3 ? 0.25 : 1.0); 
      time_step            = (scaling / (2.1 * dim * std::sqrt(1. * dim)) / 
                   (parameters.temperature_degree * get_cfl_number())); 

      const double maximal_velocity = get_maximal_velocity(); 
      pcout << "   Maximal velocity: " 
            << maximal_velocity * EquationData::year_in_seconds * 100 
            << " cm/year" << std::endl; 
      pcout << "   " 
            << "Time step: " << time_step / EquationData::year_in_seconds 
            << " years" << std::endl; 

      temperature_solution = old_temperature_solution; 
      assemble_temperature_system(maximal_velocity); 
    } 

    { 
      TimerOutput::Scope timer_section(computing_timer, 
                                       "   Solve temperature system"); 

      SolverControl solver_control(temperature_matrix.m(), 
                                   1e-12 * temperature_rhs.l2_norm()); 
      SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 

      TrilinosWrappers::MPI::Vector distributed_temperature_solution( 
        temperature_rhs); 
      distributed_temperature_solution = temperature_solution; 

      cg.solve(temperature_matrix, 
               distributed_temperature_solution, 
               temperature_rhs, 
               *T_preconditioner); 

      temperature_constraints.distribute(distributed_temperature_solution); 
      temperature_solution = distributed_temperature_solution; 

      pcout << "   " << solver_control.last_step() 
            << " CG iterations for temperature" << std::endl; 

      double temperature[2] = {std::numeric_limits<double>::max(), 
                               -std::numeric_limits<double>::max()}; 
      double global_temperature[2]; 

      for (unsigned int i = 
             distributed_temperature_solution.local_range().first; 
           i < distributed_temperature_solution.local_range().second; 
           ++i) 
        { 
          temperature[0] = 
            std::min<double>(temperature[0], 
                             distributed_temperature_solution(i)); 
          temperature[1] = 
            std::max<double>(temperature[1], 
                             distributed_temperature_solution(i)); 
        } 

      temperature[0] *= -1.0; 
      Utilities::MPI::max(temperature, MPI_COMM_WORLD, global_temperature); 
      global_temperature[0] *= -1.0; 

      pcout << "   Temperature range: " << global_temperature[0] << ' ' 
            << global_temperature[1] << std::endl; 
    } 
  } 
// @sect4{BoussinesqFlowProblem::output_results}  

// 接下来是生成输出的函数。输出的数量可以像我们在  step-31  中那样手动引入。另一种方法是把这个任务交给一个继承自DataPostprocessor类的PostProcessor，它可以被附加到DataOut。这允许我们从解决方案中输出派生量，比如本例中包含的摩擦热。它重载了虚拟函数 DataPostprocessor::evaluate_vector_field(), ，然后从 DataOut::build_patches(). 内部调用。 我们必须给它数值解、它的导数、单元的法线、实际评估点和任何额外的数量。这与 step-29 和其他程序中讨论的程序相同。

  template <int dim> 
  class BoussinesqFlowProblem<dim>::Postprocessor 
    : public DataPostprocessor<dim> 
  { 
  public: 
    Postprocessor(const unsigned int partition, const double minimal_pressure); 

    virtual void evaluate_vector_field( 
      const DataPostprocessorInputs::Vector<dim> &inputs, 
      std::vector<Vector<double>> &computed_quantities) const override; 

    virtual std::vector<std::string> get_names() const override; 

    virtual std::vector< 
      DataComponentInterpretation::DataComponentInterpretation> 
    get_data_component_interpretation() const override; 

    virtual UpdateFlags get_needed_update_flags() const override; 

  private: 
    const unsigned int partition; 
    const double       minimal_pressure; 
  }; 

  template <int dim> 
  BoussinesqFlowProblem<dim>::Postprocessor::Postprocessor( 
    const unsigned int partition, 
    const double       minimal_pressure) 
    : partition(partition) 
    , minimal_pressure(minimal_pressure) 
  {} 

// 这里我们定义了要输出的变量的名称。这些是速度、压力和温度的实际求解值，以及摩擦热和对每个单元拥有的处理器的编号。这使我们能够直观地看到处理器之间的领域划分。除了速度是矢量值的，其他的量都是标量。

  template <int dim> 
  std::vector<std::string> 
  BoussinesqFlowProblem<dim>::Postprocessor::get_names() const 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("p"); 
    solution_names.emplace_back("T"); 
    solution_names.emplace_back("friction_heating"); 
    solution_names.emplace_back("partition"); 

    return solution_names; 
  } 

  template <int dim> 
  std::vector<DataComponentInterpretation::DataComponentInterpretation> 
  BoussinesqFlowProblem<dim>::Postprocessor::get_data_component_interpretation() 
    const 
  { 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      interpretation(dim, 
                     DataComponentInterpretation::component_is_part_of_vector); 

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 

    return interpretation; 
  } 

  template <int dim> 
  UpdateFlags 
  BoussinesqFlowProblem<dim>::Postprocessor::get_needed_update_flags() const 
  { 
    return update_values | update_gradients | update_quadrature_points; 
  } 

// 现在我们实现计算派生量的函数。正如我们对输出所做的那样，我们将速度从其SI单位重新调整为更容易阅读的单位，即厘米/年。接下来，压力被缩放为0和最大压力之间。这使得它更容易比较--本质上是使所有的压力变量变成正数或零。温度按原样计算，摩擦热按  $2 \eta \varepsilon(\mathbf{u}) \cdot \varepsilon(\mathbf{u})$  计算。

// 我们在这里输出的数量更多的是为了说明问题，而不是为了实际的科学价值。我们在本程序的结果部分简要地回到这一点，并解释人们实际上可能感兴趣的是什么。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::Postprocessor::evaluate_vector_field( 
    const DataPostprocessorInputs::Vector<dim> &inputs, 
    std::vector<Vector<double>> &               computed_quantities) const 
  { 
    const unsigned int n_quadrature_points = inputs.solution_values.size(); 
    Assert(inputs.solution_gradients.size() == n_quadrature_points, 
           ExcInternalError()); 
    Assert(computed_quantities.size() == n_quadrature_points, 
           ExcInternalError()); 
    Assert(inputs.solution_values[0].size() == dim + 2, ExcInternalError()); 

    for (unsigned int q = 0; q < n_quadrature_points; ++q) 
      { 
        for (unsigned int d = 0; d < dim; ++d) 
          computed_quantities[q](d) = (inputs.solution_values[q](d) * 
                                       EquationData::year_in_seconds * 100); 

        const double pressure = 
          (inputs.solution_values[q](dim) - minimal_pressure); 
        computed_quantities[q](dim) = pressure; 

        const double temperature        = inputs.solution_values[q](dim + 1); 
        computed_quantities[q](dim + 1) = temperature; 

        Tensor<2, dim> grad_u; 
        for (unsigned int d = 0; d < dim; ++d) 
          grad_u[d] = inputs.solution_gradients[q][d]; 
        const SymmetricTensor<2, dim> strain_rate = symmetrize(grad_u); 
        computed_quantities[q](dim + 2) = 
          2 * EquationData::eta * strain_rate * strain_rate; 

        computed_quantities[q](dim + 3) = partition; 
      } 
  } 

//  <code>output_results()</code> 函数的任务与 step-31 中的类似。然而，在这里我们将演示一种不同的技术，即如何合并来自不同DoFHandler对象的输出。我们要实现这种重组的方法是创建一个联合的DoFHandler，收集两个部分，斯托克斯解和温度解。这可以通过将两个系统的有限元结合起来形成一个FES系统来很好地完成，并让这个集体系统定义一个新的DoFHandler对象。为了确保一切都做得很正确，我们进行了一次理智的检查，确保我们从斯托克斯和温度两个系统中得到了所有的道夫，甚至是在组合系统中。然后我们将数据向量合并。不幸的是，没有直接的关系告诉我们如何将斯托克斯和温度矢量分类到联合矢量中。我们可以绕过这个麻烦的方法是依靠FES系统中收集的信息。对于一个单元上的每个dof，联合有限元知道它属于哪个方程分量（速度分量、压力或温度）--这就是我们所需要的信息！这就是我们所需要的。因此，我们通过所有单元（迭代器进入所有三个DoFHandlers同步移动），对于每个联合单元dof，我们使用 FiniteElement::system_to_base_index 函数读出该分量（关于其返回值的各个部分的描述见那里）。我们还需要跟踪我们是在斯托克斯道次还是温度道次，这包含在joint_fe.system_to_base_index(i).first.first中。最终，三个系统中的任何一个系统的dof_indices数据结构都会告诉我们全局矢量和局部dof之间的关系在当前单元上是怎样的，这就结束了这项繁琐的工作。我们确保每个处理器在建立联合求解向量时，只在其本地拥有的子域上工作（而不是在幽灵或人工单元上）。然后在 DataOut::build_patches(), 中也要这样做，但该函数会自动这样做。

// 我们最终得到的是一组补丁，我们可以使用DataOutBase中的函数以各种输出格式编写补丁。在这里，我们必须注意，每个处理器所写的实际上只是它自己领域的一部分，也就是说，我们要把每个处理器的贡献写进一个单独的文件。我们通过在写解决方案时给文件名添加一个额外的数字来做到这一点。这其实并不新鲜，我们在  step-40  中也是这样做的。注意，我们用压缩格式 @p .vtu 而不是普通的vtk文件来写，这样可以节省不少存储空间。

// 所有其余的工作都在后处理程序类中完成。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::output_results() 
  { 
    TimerOutput::Scope timer_section(computing_timer, "Postprocessing"); 

    const FESystem<dim> joint_fe(stokes_fe, 1, temperature_fe, 1); 

    DoFHandler<dim> joint_dof_handler(triangulation); 
    joint_dof_handler.distribute_dofs(joint_fe); 
    Assert(joint_dof_handler.n_dofs() == 
             stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(), 
           ExcInternalError()); 

    TrilinosWrappers::MPI::Vector joint_solution; 
    joint_solution.reinit(joint_dof_handler.locally_owned_dofs(), 
                          MPI_COMM_WORLD); 

    { 
      std::vector<types::global_dof_index> local_joint_dof_indices( 
        joint_fe.n_dofs_per_cell()); 
      std::vector<types::global_dof_index> local_stokes_dof_indices( 
        stokes_fe.n_dofs_per_cell()); 
      std::vector<types::global_dof_index> local_temperature_dof_indices( 
        temperature_fe.n_dofs_per_cell()); 

      typename DoFHandler<dim>::active_cell_iterator 
        joint_cell       = joint_dof_handler.begin_active(), 
        joint_endc       = joint_dof_handler.end(), 
        stokes_cell      = stokes_dof_handler.begin_active(), 
        temperature_cell = temperature_dof_handler.begin_active(); 
      for (; joint_cell != joint_endc; 
           ++joint_cell, ++stokes_cell, ++temperature_cell) 
        if (joint_cell->is_locally_owned()) 
          { 
            joint_cell->get_dof_indices(local_joint_dof_indices); 
            stokes_cell->get_dof_indices(local_stokes_dof_indices); 
            temperature_cell->get_dof_indices(local_temperature_dof_indices); 

            for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i) 
              if (joint_fe.system_to_base_index(i).first.first == 0) 
                { 
                  Assert(joint_fe.system_to_base_index(i).second < 
                           local_stokes_dof_indices.size(), 
                         ExcInternalError()); 

                  joint_solution(local_joint_dof_indices[i]) = stokes_solution( 
                    local_stokes_dof_indices[joint_fe.system_to_base_index(i) 
                                               .second]); 
                } 
              else 
                { 
                  Assert(joint_fe.system_to_base_index(i).first.first == 1, 
                         ExcInternalError()); 
                  Assert(joint_fe.system_to_base_index(i).second < 
                           local_temperature_dof_indices.size(), 
                         ExcInternalError()); 
                  joint_solution(local_joint_dof_indices[i]) = 
                    temperature_solution( 
                      local_temperature_dof_indices 
                        [joint_fe.system_to_base_index(i).second]); 
                } 
          } 
    } 

    joint_solution.compress(VectorOperation::insert); 

    IndexSet locally_relevant_joint_dofs(joint_dof_handler.n_dofs()); 
    DoFTools::extract_locally_relevant_dofs(joint_dof_handler, 
                                            locally_relevant_joint_dofs); 
    TrilinosWrappers::MPI::Vector locally_relevant_joint_solution; 
    locally_relevant_joint_solution.reinit(locally_relevant_joint_dofs, 
                                           MPI_COMM_WORLD); 
    locally_relevant_joint_solution = joint_solution; 

    Postprocessor postprocessor(Utilities::MPI::this_mpi_process( 
                                  MPI_COMM_WORLD), 
                                stokes_solution.block(1).min()); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(joint_dof_handler); 
    data_out.add_data_vector(locally_relevant_joint_solution, postprocessor); 
    data_out.build_patches(); 

    static int out_index = 0; 
    data_out.write_vtu_with_pvtu_record( 
      "./", "solution", out_index, MPI_COMM_WORLD, 5); 

    out_index++; 
  } 

//  @sect4{BoussinesqFlowProblem::refine_mesh}  

// 这个函数也不是真正的新函数。因为我们在中间调用的 <code>setup_dofs</code> 函数有自己的定时器部分，所以我们把这个函数的定时分成两部分。这也可以让我们很容易地识别出这两个中哪个更昂贵。
//但是，
//有一点需要注意的是，我们只想在本地拥有的子域上计算错误指标。为了达到这个目的，我们向 KellyErrorEstimator::estimate 函数传递一个额外的参数。请注意，用于误差估计的向量被调整为当前进程上存在的活动单元的数量，它小于所有处理器上活动单元的总数（但大于本地拥有的活动单元的数量）；每个处理器只有本地拥有的单元周围有一些粗略的单元，这在  step-40  中也有解释。

// 本地误差估计值然后被交给GridRefinement的%并行版本（在命名空间 parallel::distributed::GridRefinement, 中，也见 step-40 ），它查看误差并通过比较各处理器的误差值找到需要细化的单元。正如在 step-31 中，我们希望限制最大的网格级别。因此，万一有些单元格已经被标记为最精细的级别，我们只需清除细化标志。

  template <int dim> 
  void 
  BoussinesqFlowProblem<dim>::refine_mesh(const unsigned int max_grid_level) 
  { 
    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> 
      temperature_trans(temperature_dof_handler); 
    parallel::distributed::SolutionTransfer<dim, 
                                            TrilinosWrappers::MPI::BlockVector> 
      stokes_trans(stokes_dof_handler); 

    { 
      TimerOutput::Scope timer_section(computing_timer, 
                                       "Refine mesh structure, part 1"); 

      Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

      KellyErrorEstimator<dim>::estimate( 
        temperature_dof_handler, 
        QGauss<dim - 1>(parameters.temperature_degree + 1), 
        std::map<types::boundary_id, const Function<dim> *>(), 
        temperature_solution, 
        estimated_error_per_cell, 
        ComponentMask(), 
        nullptr, 
        0, 
        triangulation.locally_owned_subdomain()); 

      parallel::distributed::GridRefinement::refine_and_coarsen_fixed_fraction( 
        triangulation, estimated_error_per_cell, 0.3, 0.1); 

      if (triangulation.n_levels() > max_grid_level) 
        for (typename Triangulation<dim>::active_cell_iterator cell = 
               triangulation.begin_active(max_grid_level); 
             cell != triangulation.end(); 
             ++cell) 
          cell->clear_refine_flag(); 

// 有了所有的标记，我们就可以告诉 parallel::distributed::SolutionTransfer 对象准备将数据从一个网格转移到下一个网格，当Triangulation作为 @p execute_coarsening_and_refinement() 调用的一部分通知他们时，他们就会这样做。语法类似于非%并行解决方案的传输（例外的是这里有一个指向向量项的指针就足够了）。下面函数的其余部分是在网格细化后再次设置数据结构，并在新的网格上恢复求解向量。

      std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature(2); 
      x_temperature[0] = &temperature_solution; 
      x_temperature[1] = &old_temperature_solution; 
      std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes(2); 
      x_stokes[0] = &stokes_solution; 
      x_stokes[1] = &old_stokes_solution; 

      triangulation.prepare_coarsening_and_refinement(); 

      temperature_trans.prepare_for_coarsening_and_refinement(x_temperature); 
      stokes_trans.prepare_for_coarsening_and_refinement(x_stokes); 

      triangulation.execute_coarsening_and_refinement(); 
    } 

    setup_dofs(); 

    { 
      TimerOutput::Scope timer_section(computing_timer, 
                                       "Refine mesh structure, part 2"); 

      { 
        TrilinosWrappers::MPI::Vector distributed_temp1(temperature_rhs); 
        TrilinosWrappers::MPI::Vector distributed_temp2(temperature_rhs); 

        std::vector<TrilinosWrappers::MPI::Vector *> tmp(2); 
        tmp[0] = &(distributed_temp1); 
        tmp[1] = &(distributed_temp2); 
        temperature_trans.interpolate(tmp); 

// 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。

        temperature_constraints.distribute(distributed_temp1); 
        temperature_constraints.distribute(distributed_temp2); 

        temperature_solution     = distributed_temp1; 
        old_temperature_solution = distributed_temp2; 
      } 

      { 
        TrilinosWrappers::MPI::BlockVector distributed_stokes(stokes_rhs); 
        TrilinosWrappers::MPI::BlockVector old_distributed_stokes(stokes_rhs); 

        std::vector<TrilinosWrappers::MPI::BlockVector *> stokes_tmp(2); 
        stokes_tmp[0] = &(distributed_stokes); 
        stokes_tmp[1] = &(old_distributed_stokes); 

        stokes_trans.interpolate(stokes_tmp); 

// 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。

        stokes_constraints.distribute(distributed_stokes); 
        stokes_constraints.distribute(old_distributed_stokes); 

        stokes_solution     = distributed_stokes; 
        old_stokes_solution = old_distributed_stokes; 
      } 
    } 
  } 

//  @sect4{BoussinesqFlowProblem::run}  

// 这是这个类中的最后一个控制函数。事实上，它运行了整个程序的其余部分，并且再次与  step-31  非常相似。唯一的实质性区别是我们现在使用了一个不同的网格（一个 GridGenerator::hyper_shell 而不是一个简单的立方体几何）。

  template <int dim> 
  void BoussinesqFlowProblem<dim>::run() 
  { 
    GridGenerator::hyper_shell(triangulation, 
                               Point<dim>(), 
                               EquationData::R0, 
                               EquationData::R1, 
                               (dim == 3) ? 96 : 12, 
                               true); 

    global_Omega_diameter = GridTools::diameter(triangulation); 

    triangulation.refine_global(parameters.initial_global_refinement); 

    setup_dofs(); 

    unsigned int pre_refinement_step = 0; 

  start_time_iteration: 

    { 
      TrilinosWrappers::MPI::Vector solution( 
        temperature_dof_handler.locally_owned_dofs()); 
// VectorTools::project 通过deal.II自己的本地MatrixFree框架支持具有大多数标准有限元素的并行矢量类：由于我们使用中等阶数的标准拉格朗日元素，这个函数在这里工作得很好。

      VectorTools::project(temperature_dof_handler, 
                           temperature_constraints, 
                           QGauss<dim>(parameters.temperature_degree + 2), 
                           EquationData::TemperatureInitialValues<dim>(), 
                           solution); 

// 在如此计算了当前的温度字段之后，让我们设置保存温度节点的成员变量。严格来说，我们真的只需要设置 <code>old_temperature_solution</code> ，因为我们要做的第一件事是计算斯托克斯解，它只需要前一个时间步长的温度场。尽管如此，如果我们想扩展我们的数值方法或物理模型，不初始化其他的向量也不会有什么好处（特别是这是一个相对便宜的操作，我们只需要在程序开始时做一次），所以我们也初始化 <code>old_temperature_solution</code> 和 <code>old_old_temperature_solution</code> 。这个赋值确保了左边的向量（初始化后也包含鬼魂元素）也得到了正确的鬼魂元素。换句话说，这里的赋值需要处理器之间的通信。

      temperature_solution         = solution; 
      old_temperature_solution     = solution; 
      old_old_temperature_solution = solution; 
    } 

    timestep_number = 0; 
    time_step = old_time_step = 0; 

    double time = 0; 

    do 
      { 
        pcout << "Timestep " << timestep_number 
              << ":  t=" << time / EquationData::year_in_seconds << " years" 
              << std::endl; 

        assemble_stokes_system(); 
        build_stokes_preconditioner(); 
        assemble_temperature_matrix(); 

        solve(); 

        pcout << std::endl; 

        if ((timestep_number == 0) && 
            (pre_refinement_step < parameters.initial_adaptive_refinement)) 
          { 
            refine_mesh(parameters.initial_global_refinement + 
                        parameters.initial_adaptive_refinement); 
            ++pre_refinement_step; 
            goto start_time_iteration; 
          } 
        else if ((timestep_number > 0) && 
                 (timestep_number % parameters.adaptive_refinement_interval == 
                  0)) 
          refine_mesh(parameters.initial_global_refinement + 
                      parameters.initial_adaptive_refinement); 

        if ((parameters.generate_graphical_output == true) && 
            (timestep_number % parameters.graphical_output_interval == 0)) 
          output_results(); 

// 为了加快线性求解器的速度，我们从旧的时间水平上推断出新的解决方案。这可以提供一个非常好的初始猜测，使求解器所需的迭代次数减少一半以上。我们不需要在最后一次迭代中进行推断，所以如果我们达到了最后的时间，我们就在这里停止。

// 作为一个时间步长的最后一件事（在实际提高时间步长之前），我们检查当前的时间步长是否被100整除，如果是的话，我们让计算计时器打印一个到目前为止所花费的CPU时间的总结。

        if (time > parameters.end_time * EquationData::year_in_seconds) 
          break; 

        TrilinosWrappers::MPI::BlockVector old_old_stokes_solution; 
        old_old_stokes_solution      = old_stokes_solution; 
        old_stokes_solution          = stokes_solution; 
        old_old_temperature_solution = old_temperature_solution; 
        old_temperature_solution     = temperature_solution; 
        if (old_time_step > 0) 
          { 

// Trilinos sadd不喜欢鬼魂向量，即使作为输入。暂时复制到分布式向量中。

            { 
              TrilinosWrappers::MPI::BlockVector distr_solution(stokes_rhs); 
              distr_solution = stokes_solution; 
              TrilinosWrappers::MPI::BlockVector distr_old_solution(stokes_rhs); 
              distr_old_solution = old_old_stokes_solution; 
              distr_solution.sadd(1. + time_step / old_time_step, 
                                  -time_step / old_time_step, 
                                  distr_old_solution); 
              stokes_solution = distr_solution; 
            } 
            { 
              TrilinosWrappers::MPI::Vector distr_solution(temperature_rhs); 
              distr_solution = temperature_solution; 
              TrilinosWrappers::MPI::Vector distr_old_solution(temperature_rhs); 
              distr_old_solution = old_old_temperature_solution; 
              distr_solution.sadd(1. + time_step / old_time_step, 
                                  -time_step / old_time_step, 
                                  distr_old_solution); 
              temperature_solution = distr_solution; 
            } 
          } 

        if ((timestep_number > 0) && (timestep_number % 100 == 0)) 
          computing_timer.print_summary(); 

        time += time_step; 
        ++timestep_number; 
      } 
    while (true); 

// 如果我们要生成图形输出，也要对最后一个时间步骤这样做，除非我们在离开do-while循环之前刚刚这样做。

    if ((parameters.generate_graphical_output == true) && 
        !((timestep_number - 1) % parameters.graphical_output_interval == 0)) 
      output_results(); 
  } 
} // namespace Step32 

//  @sect3{The <code>main</code> function}  

// 主函数像往常一样简短，与  step-31  中的函数非常相似。由于我们使用了一个参数文件，该文件在命令行中被指定为参数，所以我们必须在这里读取它，并将其传递给参数类进行解析。如果命令行中没有给出文件名，我们就简单地使用与程序一起分发的  <code>\step-32.prm</code>  文件。

// 由于三维计算非常缓慢，除非你投入大量的处理器，程序默认为二维。你可以通过把下面的常数维度改为3来获得三维版本。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace Step32; 
      using namespace dealii; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 

      std::string parameter_filename; 
      if (argc >= 2) 
        parameter_filename = argv[1]; 
      else 
        parameter_filename = "step-32.prm"; 

      const int                              dim = 2; 
      BoussinesqFlowProblem<dim>::Parameters parameters(parameter_filename); 
      BoussinesqFlowProblem<dim>             flow_problem(parameters); 
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


