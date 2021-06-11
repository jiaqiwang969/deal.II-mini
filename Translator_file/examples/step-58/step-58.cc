CCTest_file/step-58.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE at 
 * the top level of the deal.II distribution. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Author: Wolfgang Bangerth, Colorado State University 
 *         Yong-Yong Cai, Beijing Computational Science Research Center 
 */ 


// @sect3{Include files}  程序以通常的包含文件开始，所有这些文件你现在应该都见过了。

#include <deal.II/base/logstream.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/matrix_tools.h> 

#include <fstream> 
#include <iostream> 

// 然后按照惯例将这个程序的所有内容放入一个命名空间，并将deal.II命名空间导入到我们将要工作的命名空间中。

namespace Step58 
{ 
  using namespace dealii; 
// @sect3{The <code>NonlinearSchroedingerEquation</code> class}  

// 然后是主类。它看起来非常像  step-4  或  step-6  中的相应类，唯一的例外是，矩阵和向量以及其他所有与线性系统相关的元素现在都存储为  `std::complex<double>`  类型，而不仅仅是 `double`。

  template <int dim> 
  class NonlinearSchroedingerEquation 
  { 
  public: 
    NonlinearSchroedingerEquation(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_matrices(); 
    void do_half_phase_step(); 
    void do_full_spatial_step(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<std::complex<double>> constraints; 

    SparsityPattern                    sparsity_pattern; 
    SparseMatrix<std::complex<double>> system_matrix; 
    SparseMatrix<std::complex<double>> rhs_matrix; 

    Vector<std::complex<double>> solution; 
    Vector<std::complex<double>> system_rhs; 

    double       time; 
    double       time_step; 
    unsigned int timestep_number; 

    double kappa; 
  }; 

//  @sect3{Equation data}  

// 在我们继续填写主类的细节之前，让我们定义与问题相对应的方程数据，即初始值，以及一个右手类。我们将把初始条件也用于边界值，我们只是保持边界值不变）。我们使用派生自Function类模板的类来做这件事，这个模板之前已经用过很多次了，所以下面的内容看起来并不令人惊讶。唯一值得注意的是，我们这里有一个复值问题，所以我们必须提供Function类的第二个模板参数（否则会默认为`double`）。此外，`value()`函数的返回类型当然也是复数。

// 这些函数精确地返回什么，在介绍部分的最后已经讨论过了。

  template <int dim> 
  class InitialValues : public Function<dim, std::complex<double>> 
  { 
  public: 
    InitialValues() 
      : Function<dim, std::complex<double>>(1) 
    {} 

    virtual std::complex<double> 
    value(const Point<dim> &p, const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  std::complex<double> 
  InitialValues<dim>::value(const Point<dim> & p, 
                            const unsigned int component) const 
  { 
    static_assert(dim == 2, "This initial condition only works in 2d."); 

    (void)component; 
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 

    const std::vector<Point<dim>> vortex_centers = {{0, -0.3}, 
                                                    {0, +0.3}, 
                                                    {+0.3, 0}, 
                                                    {-0.3, 0}}; 

    const double R = 0.1; 
    const double alpha = 
      1. / (std::pow(R, dim) * std::pow(numbers::PI, dim / 2.)); 

    double sum = 0; 
    for (const auto &vortex_center : vortex_centers) 
      { 
        const Tensor<1, dim> distance = p - vortex_center; 
        const double         r        = distance.norm(); 

        sum += alpha * std::exp(-(r * r) / (R * R)); 
      } 

    return {std::sqrt(sum), 0.}; 
  } 

  template <int dim> 
  class Potential : public Function<dim> 
  { 
  public: 
    Potential() = default; 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double Potential<dim>::value(const Point<dim> & p, 
                               const unsigned int component) const 
  { 
    (void)component; 
    Assert(component == 0, ExcIndexRange(component, 0, 1)); 

    return (Point<dim>().distance(p) > 0.7 ? 1000 : 0); 
  } 

//  @sect3{Implementation of the <code>NonlinearSchroedingerEquation</code> class}  

// 我们首先指定了类的构造函数的实现。

  template <int dim> 
  NonlinearSchroedingerEquation<dim>::NonlinearSchroedingerEquation() 
    : fe(2) 
    , dof_handler(triangulation) 
    , time(0) 
    , time_step(1. / 128) 
    , timestep_number(0) 
    , kappa(1) 
  {} 
// @sect4{Setting up data structures and assembling matrices}  

// 下一个函数是在程序开始时，也就是在第一个时间步骤之前，设置网格、DoFHandler以及矩阵和向量。如果你已经阅读了至少到 step-6 为止的教程程序，那么前几行是相当标准的。

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::setup_system() 
  { 
    GridGenerator::hyper_cube(triangulation, -1, 1); 
    triangulation.refine_global(6); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl; 

    dof_handler.distribute_dofs(fe); 

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl 
              << std::endl; 

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
    rhs_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 

    constraints.close(); 
  } 

// 接下来，我们组装相关的矩阵。按照我们对斯特朗分裂的空间步骤（即每个时间步骤中三个部分步骤中的第二个步骤）的Crank-Nicolson离散化的写法，我们被引导到线性系统  
// $\left[ -iM  +  \frac 14 k_{n+1} A + \frac 12 k_{n+1} W \right]
//    \Psi^{(n,2)}
//   =
//   \left[ -iM  -  \frac 14 k_{n+1} A - \frac 12 k_{n+1} W \right]
//    \Psi^{(n,1)}$ 
   
//     换句话说，这里有两个矩阵在起作用--一个用于左手边，一个用于右手边。我们分别建立这些矩阵。我们可以避免建立右手边的矩阵，而只是在每个时间步长中形成矩阵的*作用* $\Psi^{(n,1)}$ 。这可能更有效，也可能不有效，但是对于这个程序来说，效率并不是最重要的）。)

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::assemble_matrices() 
  { 
    const QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<std::complex<double>> cell_matrix_lhs(dofs_per_cell, 
                                                     dofs_per_cell); 
    FullMatrix<std::complex<double>> cell_matrix_rhs(dofs_per_cell, 
                                                     dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
    std::vector<double>                  potential_values(n_q_points); 
    const Potential<dim>                 potential; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix_lhs = std::complex<double>(0.); 
        cell_matrix_rhs = std::complex<double>(0.); 

        fe_values.reinit(cell); 

        potential.value_list(fe_values.get_quadrature_points(), 
                             potential_values); 

        for (unsigned int q_index = 0; q_index < n_q_points; ++q_index) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                for (unsigned int l = 0; l < dofs_per_cell; ++l) 
                  { 
                    const std::complex<double> i = {0, 1}; 

                    cell_matrix_lhs(k, l) += 
                      (-i * fe_values.shape_value(k, q_index) * 
                         fe_values.shape_value(l, q_index) + 
                       time_step / 4 * fe_values.shape_grad(k, q_index) * 
                         fe_values.shape_grad(l, q_index) + 
                       time_step / 2 * potential_values[q_index] * 
                         fe_values.shape_value(k, q_index) * 
                         fe_values.shape_value(l, q_index)) * 
                      fe_values.JxW(q_index); 

                    cell_matrix_rhs(k, l) += 
                      (-i * fe_values.shape_value(k, q_index) * 
                         fe_values.shape_value(l, q_index) - 
                       time_step / 4 * fe_values.shape_grad(k, q_index) * 
                         fe_values.shape_grad(l, q_index) - 
                       time_step / 2 * potential_values[q_index] * 
                         fe_values.shape_value(k, q_index) * 
                         fe_values.shape_value(l, q_index)) * 
                      fe_values.JxW(q_index); 
                  } 
              } 
          } 

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global(cell_matrix_lhs, 
                                               local_dof_indices, 
                                               system_matrix); 
        constraints.distribute_local_to_global(cell_matrix_rhs, 
                                               local_dof_indices, 
                                               rhs_matrix); 
      } 
  } 
// @sect4{Implementing the Strang splitting steps}  

// 在建立了上述所有数据结构后，我们现在可以实现构成斯特朗分裂方案的部分步骤。我们从推进阶段的半步开始，这被用作每个时间步骤的第一和最后部分。

// 为此，回顾一下，对于第一个半步，我们需要计算  $\psi^{(n,1)} = e^{-i\kappa|\psi^{(n,0)}|^2 \tfrac 12\Delta t} \; \psi^{(n,0)}$  。这里， $\psi^{(n,0)}=\psi^{(n)}$ 和 $\psi^{(n,1)}$ 是空间的函数，分别对应于前一个完整时间步骤的输出和三个部分步骤中第一个步骤的结果。必须为第三个部分步骤计算相应的解决方案，即  $\psi^{(n,3)} = e^{-i\kappa|\psi^{(n,2)}|^2 \tfrac 12\Delta t} \; \psi^{(n,2)}$  ，其中  $\psi^{(n,3)}=\psi^{(n+1)}$  是整个时间步骤的结果，其输入  $\psi^{(n,2)}$  是斯特朗分割的空间步骤的结果。

// 一个重要的认识是，虽然 $\psi^{(n,0)}(\mathbf x)$ 可能是一个有限元函数（即，是片状多项式），但对于我们使用指数因子更新相位的 "旋转 "函数来说，不一定是这样的（回顾一下，该函数的振幅在该步骤中保持不变）。换句话说，我们可以在每一个点 $\psi^{(n,1)}(\mathbf x)$ *计算 $\mathbf x\in\Omega$ ，但我们不能在网格上表示它，因为它不是一个片状多项式函数。在一个离散的环境中，我们能做的最好的事情就是计算一个投影或内插。换句话说，我们可以计算 $\psi_h^{(n,1)}(\mathbf x) = \Pi_h \left(e^{-i\kappa|\psi_h^{(n,0)}(\mathbf x)|^2 \tfrac 12\Delta t} \; \psi_h^{(n,0)}(\mathbf x) \right)$ ，其中 $\Pi_h$ 是一个投影或内插算子。如果我们选择插值，情况就特别简单。那么，我们需要计算的就是*在节点点上的右手边的值，并将这些作为自由度向量 $\Psi^{(n,1)}$ 的节点值。这很容易做到，因为在这里使用的拉格朗日有限元的节点点上评估右手边，需要我们只看节点向量的一个（复值）条目。换句话说，我们需要做的是计算 $\Psi^{(n,1)}_j = e^{-i\kappa|\Psi^{(n,0)}_j|^2 \tfrac 12\Delta t} \; \Psi^{(n,0)}_j$ ，其中 $j$ 在我们的解向量的所有条目上循环。这就是下面的函数所做的--事实上，它甚至没有为 $\Psi^{(n,0)}$ 和 $\Psi^{(n,1)}$ 使用单独的向量，而只是适当地更新同一个向量。

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::do_half_phase_step() 
  { 
    for (auto &value : solution) 
      { 
        const std::complex<double> i         = {0, 1}; 
        const double               magnitude = std::abs(value); 

        value = std::exp(-i * kappa * magnitude * magnitude * (time_step / 2)) * 
                value; 
      } 
  } 

// 下一步是求解每个时间步骤中的线性系统，即我们使用的Strang分割的后半步。记得它的形式是 $C\Psi^{(n,2)} = R\Psi^{(n,1)}$ ，其中 $C$ 和 $R$ 是我们之前组装的矩阵。

// 我们在这里解决这个问题的方法是使用直接求解器。我们首先使用 $r=R\Psi^{(n,1)}$ 函数形成右边的 SparseMatrix::vmult() ，并将结果放入`system_rhs`变量。然后我们调用 SparseDirectUMFPACK::solver() ，该函数以矩阵 $C$ 和右手边的向量为参数，并在同一向量`system_rhs`中返回解。最后一步是将计算出的解放回`solution`变量中。

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::do_full_spatial_step() 
  { 
    rhs_matrix.vmult(system_rhs, solution); 

    SparseDirectUMFPACK direct_solver; 
    direct_solver.solve(system_matrix, system_rhs); 

    solution = system_rhs; 
  } 

//  @sect4{Creating graphical output}  

// 我们应该讨论的最后一个辅助函数和类是那些创建图形输出的函数。对斯特朗分裂的局部和空间部分运行半步和全步的结果是，我们在每个时间步数结束时将`solution`向量 $\Psi^n$ 更新为正确的值。它的条目包含有限元网格节点上的解的复数。

// 复数不容易被视觉化。我们可以输出它们的实部和虚部，即字段 $\text{Re}(\psi_h^{(n)}(\mathbf x))$ 和 $\text{Im}(\psi_h^{(n)}(\mathbf x))$ ，这正是DataOut类在通过 DataOut::add_data_vector() 附加复数向量，然后调用 DataOut::build_patches(). 时所做的事情，这确实是我们下面要做的。

// 但很多时候，我们对解向量的实部和虚部并不特别感兴趣，而是对解的幅度 $|\psi|$ 和相位角 $\text{arg}(\psi)$ 等衍生量感兴趣。在这里这样的量子系统的背景下，幅度本身并不那么有趣，相反，"振幅"， $|\psi|^2$ 才是一个物理属性：它对应于在一个特定的状态场所找到一个粒子的概率密度。将计算出的量放入输出文件以实现可视化的方法--正如在以前的许多教程程序中使用的那样--是使用数据后处理程序和派生类的设施。具体来说，一个复数的振幅和它的相位角都是标量，因此DataPostprocessorScalar类是我们要做的正确工具。

// 因此，我们在这里要做的是实现两个类`ComplexAmplitude`和`ComplexPhase`，为DataOut决定生成输出的每个点计算解决方案的振幅 $|\psi_h|^2$ 和相位 $\text{arg}(\psi_h)$ ，以便进行可视化。下面有大量的模板代码，这两个类中的第一个唯一有趣的部分是它的`evaluate_vector_field()`函数如何计算`computed_quantities`对象。

//（还有一个相当尴尬的事实是，<a
//  href="https:en.cppreference.com/w/cpp/numeric/complex/norm">std::norm()</a>函数并没有计算人们天真的想象，即 $|\psi|$  ，而是返回 $|\psi|^2$ 。一个标准函数以这样的方式被错误地命名，这当然是相当令人困惑的......)

  namespace DataPostprocessors 
  { 
    template <int dim> 
    class ComplexAmplitude : public DataPostprocessorScalar<dim> 
    { 
    public: 
      ComplexAmplitude(); 

      virtual void evaluate_vector_field( 
        const DataPostprocessorInputs::Vector<dim> &inputs, 
        std::vector<Vector<double>> &computed_quantities) const override; 
    }; 

    template <int dim> 
    ComplexAmplitude<dim>::ComplexAmplitude() 
      : DataPostprocessorScalar<dim>("Amplitude", update_values) 
    {} 

    template <int dim> 
    void ComplexAmplitude<dim>::evaluate_vector_field( 
      const DataPostprocessorInputs::Vector<dim> &inputs, 
      std::vector<Vector<double>> &               computed_quantities) const 
    { 
      Assert(computed_quantities.size() == inputs.solution_values.size(), 
             ExcDimensionMismatch(computed_quantities.size(), 
                                  inputs.solution_values.size())); 

      for (unsigned int q = 0; q < computed_quantities.size(); ++q) 
        { 
          Assert(computed_quantities[q].size() == 1, 
                 ExcDimensionMismatch(computed_quantities[q].size(), 1)); 
          Assert(inputs.solution_values[q].size() == 2, 
                 ExcDimensionMismatch(inputs.solution_values[q].size(), 2)); 

          const std::complex<double> psi(inputs.solution_values[q](0), 
                                         inputs.solution_values[q](1)); 
          computed_quantities[q](0) = std::norm(psi); 
        } 
    } 

// 这些后处理程序类中的第二个是计算每一个点的复值解决方案的相位角。换句话说，如果我们表示  $\psi(\mathbf x,t)=r(\mathbf x,t) e^{i\varphi(\mathbf x,t)}$  ，那么这个类就会计算  $\varphi(\mathbf x,t)$  。函数 <a href="https:en.cppreference.com/w/cpp/numeric/complex/arg">std::arg</a> 为我们做这个，并将角度作为实数返回  $-\pi$  和  $+\pi$  之间。

// 由于我们将在结果部分详细解释的原因，我们实际上没有在产生输出的每个位置输出这个值。相反，我们取相位所有评估点的最大值，然后用这个最大值填充每个评估点的输出字段--实质上，我们将相位角作为一个片状常数字段输出，其中每个单元都有自己的常数值。一旦你读完下面的讨论就会明白其中的原因。

    template <int dim> 
    class ComplexPhase : public DataPostprocessorScalar<dim> 
    { 
    public: 
      ComplexPhase(); 

      virtual void evaluate_vector_field( 
        const DataPostprocessorInputs::Vector<dim> &inputs, 
        std::vector<Vector<double>> &computed_quantities) const override; 
    }; 

    template <int dim> 
    ComplexPhase<dim>::ComplexPhase() 
      : DataPostprocessorScalar<dim>("Phase", update_values) 
    {} 

    template <int dim> 
    void ComplexPhase<dim>::evaluate_vector_field( 
      const DataPostprocessorInputs::Vector<dim> &inputs, 
      std::vector<Vector<double>> &               computed_quantities) const 
    { 
      Assert(computed_quantities.size() == inputs.solution_values.size(), 
             ExcDimensionMismatch(computed_quantities.size(), 
                                  inputs.solution_values.size())); 

      double max_phase = -numbers::PI; 
      for (unsigned int q = 0; q < computed_quantities.size(); ++q) 
        { 
          Assert(computed_quantities[q].size() == 1, 
                 ExcDimensionMismatch(computed_quantities[q].size(), 1)); 
          Assert(inputs.solution_values[q].size() == 2, 
                 ExcDimensionMismatch(inputs.solution_values[q].size(), 2)); 

          max_phase = 
            std::max(max_phase, 
                     std::arg( 
                       std::complex<double>(inputs.solution_values[q](0), 
                                            inputs.solution_values[q](1)))); 
        } 

      for (auto &output : computed_quantities) 
        output(0) = max_phase; 
    } 

  } // namespace DataPostprocessors 

// 在这样实现了这些后处理程序后，我们像往常一样创建输出。与其他许多时间相关的教程程序一样，我们给DataOut附加标志，表示时间步数和当前模拟时间。

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::output_results() const 
  { 
    const DataPostprocessors::ComplexAmplitude<dim> complex_magnitude; 
    const DataPostprocessors::ComplexPhase<dim>     complex_phase; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "Psi"); 
    data_out.add_data_vector(solution, complex_magnitude); 
    data_out.add_data_vector(solution, complex_phase); 
    data_out.build_patches(); 

    data_out.set_flags(DataOutBase::VtkFlags(time, timestep_number)); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect4{Running the simulation}  

// 剩下的步骤是我们如何设置这个程序的整体逻辑。这其实是比较简单的。设置数据结构；将初始条件插值到有限元空间；然后迭代所有时间步长，在每个时间步长上执行斯特朗分割法的三个部分。每隔10个时间步长，我们就生成图形输出。这就是了。

  template <int dim> 
  void NonlinearSchroedingerEquation<dim>::run() 
  { 
    setup_system(); 
    assemble_matrices(); 

    time = 0; 
    VectorTools::interpolate(dof_handler, InitialValues<dim>(), solution); 
    output_results(); 

    const double end_time = 1; 
    for (; time <= end_time; time += time_step) 
      { 
        ++timestep_number; 

        std::cout << "Time step " << timestep_number << " at t=" << time 
                  << std::endl; 

        do_half_phase_step(); 
        do_full_spatial_step(); 
        do_half_phase_step(); 

        if (timestep_number % 1 == 0) 
          output_results(); 
      } 
  } 
} // namespace Step58 

//  @sect4{The main() function}  

// 其余的又是锅炉板，和以前几乎所有的教程程序完全一样。

int main() 
{ 
  try 
    { 
      using namespace Step58; 

      NonlinearSchroedingerEquation<2> nse; 
      nse.run(); 
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

