

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2014 - 2021 by the deal.II authors 
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
 * Authors: Damien Lebrun-Grandie, Bruno Turcksin, 2014 
 */ 


// @sect3{Include files}  

// 像往常一样，第一个任务是包括这些著名的deal.II库文件和一些C++头文件的功能。

#include <deal.II/base/discrete_time.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/quadrature_lib.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_out.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparse_direct.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iostream> 
#include <cmath> 
#include <map> 

// 这是唯一一个新的包含文件：它包括所有的Runge-Kutta方法。

#include <deal.II/base/time_stepping.h> 

// 接下来的步骤与之前所有的教程程序一样。我们把所有的东西放到一个自己的命名空间中，然后把deal.II的类和函数导入其中。

namespace Step52 
{ 
  using namespace dealii; 
// @sect3{The <code>Diffusion</code> class}  

// 下一块是主类的声明。这个类中的大部分函数并不新鲜，在以前的教程中已经解释过了。唯一有趣的函数是  <code>evaluate_diffusion()</code>  和  <code>id_minus_tau_J_inverse()</code>. <code>evaluate_diffusion()</code>  评估扩散方程，  $M^{-1}(f(t,y))$  ，在一个给定的时间和一个给定的  $y$  。  <code>id_minus_tau_J_inverse()</code>  在给定的时间和给定的 $\tau$ 和 $y$ 下，评估 $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ 或类似的 $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$ 。当使用隐式方法时，就需要这个函数。

  class Diffusion 
  { 
  public: 
    Diffusion(); 

    void run(); 

  private: 
    void setup_system(); 

    void assemble_system(); 

    double get_source(const double time, const Point<2> &point) const; 

    Vector<double> evaluate_diffusion(const double          time, 
                                      const Vector<double> &y) const; 

    Vector<double> id_minus_tau_J_inverse(const double          time, 
                                          const double          tau, 
                                          const Vector<double> &y); 

    void output_results(const double                     time, 
                        const unsigned int               time_step, 
                        TimeStepping::runge_kutta_method method) const; 

// 接下来的三个函数分别是显式方法、隐式方法和嵌入式显式方法的驱动。嵌入显式方法的驱动函数返回执行的步数，鉴于它只接受作为参数传递的时间步数作为提示，但内部计算了最佳时间步数本身。

    void explicit_method(const TimeStepping::runge_kutta_method method, 
                         const unsigned int                     n_time_steps, 
                         const double                           initial_time, 
                         const double                           final_time); 

    void implicit_method(const TimeStepping::runge_kutta_method method, 
                         const unsigned int                     n_time_steps, 
                         const double                           initial_time, 
                         const double                           final_time); 

    unsigned int 
    embedded_explicit_method(const TimeStepping::runge_kutta_method method, 
                             const unsigned int n_time_steps, 
                             const double       initial_time, 
                             const double       final_time); 

    const unsigned int fe_degree; 

    const double diffusion_coefficient; 
    const double absorption_cross_section; 

    Triangulation<2> triangulation; 

    const FE_Q<2> fe; 

    DoFHandler<2> dof_handler; 

    AffineConstraints<double> constraint_matrix; 

    SparsityPattern sparsity_pattern; 

    SparseMatrix<double> system_matrix; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> mass_minus_tau_Jacobian; 

    SparseDirectUMFPACK inverse_mass_matrix; 

    Vector<double> solution; 
  }; 

// 我们选择二次方有限元，并初始化参数。

  Diffusion::Diffusion() 
    : fe_degree(2) 
    , diffusion_coefficient(1. / 30.) 
    , absorption_cross_section(1.) 
    , fe(fe_degree) 
    , dof_handler(triangulation) 
  {} 

// 现在，我们创建约束矩阵和稀疏模式。然后，我们初始化这些矩阵和求解向量。

  void Diffusion::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             1, 
                                             Functions::ZeroFunction<2>(), 
                                             constraint_matrix); 
    constraint_matrix.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraint_matrix); 
    sparsity_pattern.copy_from(dsp); 

    system_matrix.reinit(sparsity_pattern); 
    mass_matrix.reinit(sparsity_pattern); 
    mass_minus_tau_Jacobian.reinit(sparsity_pattern); 
    solution.reinit(dof_handler.n_dofs()); 
  } 

//  @sect4{<code>Diffusion::assemble_system</code>}  在这个函数中，我们计算  $-\int D \nabla b_i \cdot \nabla b_j d\boldsymbol{r} - \int \Sigma_a b_i b_j d\boldsymbol{r}$  和质量矩阵  $\int b_i b_j d\boldsymbol{r}$  。然后使用直接求解器对质量矩阵进行反演；然后 <code>inverse_mass_matrix</code> 变量将存储质量矩阵的反值，这样 $M^{-1}$ 就可以使用该对象的 <code>vmult()</code> 函数应用于一个矢量。在内部，UMFPACK并没有真正存储矩阵的逆，而是存储它的LU因子；应用矩阵的逆相当于用这两个因子做一次正解和一次逆解，这与应用矩阵的显式逆具有相同的复杂性）。

  void Diffusion::assemble_system() 
  { 
    system_matrix = 0.; 
    mass_matrix   = 0.; 

    const QGauss<2> quadrature_formula(fe_degree + 1); 

    FEValues<2> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_gradients | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_matrix      = 0.; 
        cell_mass_matrix = 0.; 

        fe_values.reinit(cell); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              { 
                cell_matrix(i, j) += 
                  ((-diffusion_coefficient *                // (-D 
                      fe_values.shape_grad(i, q_point) *    //  * grad phi_i 
                      fe_values.shape_grad(j, q_point)      //  * grad phi_j 
                    - absorption_cross_section *            //  -Sigma 
                        fe_values.shape_value(i, q_point) * //  * phi_i 
                        fe_values.shape_value(j, q_point))  //  * phi_j) 
                   * fe_values.JxW(q_point));               // * dx 
                cell_mass_matrix(i, j) += fe_values.shape_value(i, q_point) * 
                                          fe_values.shape_value(j, q_point) * 
                                          fe_values.JxW(q_point); 
              } 

        cell->get_dof_indices(local_dof_indices); 

        constraint_matrix.distribute_local_to_global(cell_matrix, 
                                                     local_dof_indices, 
                                                     system_matrix); 
        constraint_matrix.distribute_local_to_global(cell_mass_matrix, 
                                                     local_dof_indices, 
                                                     mass_matrix); 
      } 

    inverse_mass_matrix.initialize(mass_matrix); 
  } 

//  @sect4{<code>Diffusion::get_source</code>}  

// 在这个函数中，计算出特定时间和特定点的方程的源项。

  double Diffusion::get_source(const double time, const Point<2> &point) const 
  { 
    const double intensity = 10.; 
    const double frequency = numbers::PI / 10.; 
    const double b         = 5.; 
    const double x         = point(0); 

    return intensity * 
           (frequency * std::cos(frequency * time) * (b * x - x * x) + 
            std::sin(frequency * time) * 
              (absorption_cross_section * (b * x - x * x) + 
               2. * diffusion_coefficient)); 
  } 

//  @sect4{<code>Diffusion::evaluate_diffusion</code>}  

// 接下来，我们在给定的时间  $t$  和给定的矢量  $y$  评价扩散方程的弱形式。换句话说，正如介绍中所述，我们评估  $M^{-1}(-{\cal D}y - {\cal A}y + {\cal S})$  。为此，我们必须将矩阵 $-{\cal D} - {\cal A}$ （之前计算并存储在变量 <code>system_matrix</code> 中）应用于 $y$ ，然后添加源项，我们像通常那样进行积分。(如果你想节省几行代码，或者想利用并行积分的优势，可以用 VectorTools::create_right_hand_side() 来进行积分。) 然后将结果乘以 $M^{-1}$  。

  Vector<double> Diffusion::evaluate_diffusion(const double          time, 
                                               const Vector<double> &y) const 
  { 
    Vector<double> tmp(dof_handler.n_dofs()); 
    tmp = 0.; 
    system_matrix.vmult(tmp, y); 

    const QGauss<2> quadrature_formula(fe_degree + 1); 

    FEValues<2> fe_values(fe, 
                          quadrature_formula, 
                          update_values | update_quadrature_points | 
                            update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double> cell_source(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_source = 0.; 

        fe_values.reinit(cell); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          { 
            const double source = 
              get_source(time, fe_values.quadrature_point(q_point)); 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_source(i) += fe_values.shape_value(i, q_point) * // phi_i(x) 
                                source *                            // * S(x) 
                                fe_values.JxW(q_point);             // * dx 
          } 

        cell->get_dof_indices(local_dof_indices); 

        constraint_matrix.distribute_local_to_global(cell_source, 
                                                     local_dof_indices, 
                                                     tmp); 
      } 

    Vector<double> value(dof_handler.n_dofs()); 
    inverse_mass_matrix.vmult(value, tmp); 

    return value; 
  } 
// @sect4{<code>Diffusion::id_minus_tau_J_inverse</code>}  

// 我们计算  $\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} M$  。这要分几个步骤进行。 

// - 计算  $M-\tau \frac{\partial f}{\partial y}$  。  

// - 反转矩阵，得到  $\left(M-\tau \frac{\partial f} {\partial y}\right)^{-1}$  。  

// --计算 $tmp=My$ 。  

// --计算 $z=\left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} tmp =  \left(M-\tau \frac{\partial f}{\partial y}\right)^{-1} My$ 。  

// - 返回z。

  Vector<double> Diffusion::id_minus_tau_J_inverse(const double /*time*/, 
                                                   const double          tau, 
                                                   const Vector<double> &y) 
  { 
    SparseDirectUMFPACK inverse_mass_minus_tau_Jacobian; 

    mass_minus_tau_Jacobian.copy_from(mass_matrix); 
    mass_minus_tau_Jacobian.add(-tau, system_matrix); 

    inverse_mass_minus_tau_Jacobian.initialize(mass_minus_tau_Jacobian); 

    Vector<double> tmp(dof_handler.n_dofs()); 
    mass_matrix.vmult(tmp, y); 

    Vector<double> result(y); 
    inverse_mass_minus_tau_Jacobian.vmult(result, tmp); 

    return result; 
  } 

//  @sect4{<code>Diffusion::output_results</code>}  

// 下面的函数将解决方案以vtu文件的形式输出，并以时间步长和时间步长方法的名称为索引。当然，所有时间步长方法的（精确）结果应该是一样的，但这里的输出至少可以让我们对它们进行比较。

  void Diffusion::output_results(const double                     time, 
                                 const unsigned int               time_step, 
                                 TimeStepping::runge_kutta_method method) const 
  { 
    std::string method_name; 

    switch (method) 
      { 
        case TimeStepping::FORWARD_EULER: 
          { 
            method_name = "forward_euler"; 
            break; 
          } 
        case TimeStepping::RK_THIRD_ORDER: 
          { 
            method_name = "rk3"; 
            break; 
          } 
        case TimeStepping::RK_CLASSIC_FOURTH_ORDER: 
          { 
            method_name = "rk4"; 
            break; 
          } 
        case TimeStepping::BACKWARD_EULER: 
          { 
            method_name = "backward_euler"; 
            break; 
          } 
        case TimeStepping::IMPLICIT_MIDPOINT: 
          { 
            method_name = "implicit_midpoint"; 
            break; 
          } 
        case TimeStepping::SDIRK_TWO_STAGES: 
          { 
            method_name = "sdirk"; 
            break; 
          } 
        case TimeStepping::HEUN_EULER: 
          { 
            method_name = "heun_euler"; 
            break; 
          } 
        case TimeStepping::BOGACKI_SHAMPINE: 
          { 
            method_name = "bocacki_shampine"; 
            break; 
          } 
        case TimeStepping::DOPRI: 
          { 
            method_name = "dopri"; 
            break; 
          } 
        case TimeStepping::FEHLBERG: 
          { 
            method_name = "fehlberg"; 
            break; 
          } 
        case TimeStepping::CASH_KARP: 
          { 
            method_name = "cash_karp"; 
            break; 
          } 
        default: 
          { 
            break; 
          } 
      } 

    DataOut<2> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 

    data_out.build_patches(); 

    data_out.set_flags(DataOutBase::VtkFlags(time, time_step)); 

    const std::string filename = "solution_" + method_name + "-" + 
                                 Utilities::int_to_string(time_step, 3) + 
                                 ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 

    static std::vector<std::pair<double, std::string>> times_and_names; 

    static std::string method_name_prev = ""; 
    static std::string pvd_filename; 
    if (method_name_prev != method_name) 
      { 
        times_and_names.clear(); 
        method_name_prev = method_name; 
        pvd_filename     = "solution_" + method_name + ".pvd"; 
      } 
    times_and_names.emplace_back(time, filename); 
    std::ofstream pvd_output(pvd_filename); 
    DataOutBase::write_pvd_record(pvd_output, times_and_names); 
  } 
// @sect4{<code>Diffusion::explicit_method</code>}  

// 这个函数是所有显式方法的驱动。在顶部，它初始化了时间步长和解决方案（通过将其设置为零，然后确保边界值和悬挂节点约束得到尊重；当然，对于我们在这里使用的网格，悬挂节点约束实际上不是一个问题）。然后调用 <code>evolve_one_time_step</code> ，执行一个时间步骤。时间是通过DiscreteTime对象来存储和增加的。

// 对于显式方法， <code>evolve_one_time_step</code> 需要评估 $M^{-1}(f(t,y))$ ，也就是说，它需要 <code>evaluate_diffusion</code>  。因为 <code>evaluate_diffusion</code> 是一个成员函数，它需要被绑定到 <code>this</code> 。在每个进化步骤之后，我们再次应用正确的边界值和悬挂节点约束。

// 最后，每隔10个时间步骤就会输出解决方案。

  void Diffusion::explicit_method(const TimeStepping::runge_kutta_method method, 
                                  const unsigned int n_time_steps, 
                                  const double       initial_time, 
                                  const double       final_time) 
  { 
    const double time_step = 
      (final_time - initial_time) / static_cast<double>(n_time_steps); 

    solution = 0.; 
    constraint_matrix.distribute(solution); 

    TimeStepping::ExplicitRungeKutta<Vector<double>> explicit_runge_kutta( 
      method); 
    output_results(initial_time, 0, method); 
    DiscreteTime time(initial_time, final_time, time_step); 
    while (time.is_at_end() == false) 
      { 
        explicit_runge_kutta.evolve_one_time_step( 
          [this](const double time, const Vector<double> &y) { 
            return this->evaluate_diffusion(time, y); 
          }, 
          time.get_current_time(), 
          time.get_next_step_size(), 
          solution); 
        time.advance_time(); 

        constraint_matrix.distribute(solution); 

        if (time.get_step_number() % 10 == 0) 
          output_results(time.get_current_time(), 
                         time.get_step_number(), 
                         method); 
      } 
  } 

//  @sect4{<code>Diffusion::implicit_method</code>}  这个函数等同于 <code>explicit_method</code> ，但用于隐式方法。当使用隐式方法时，我们需要评估 $M^{-1}(f(t,y))$ 和 $\left(I-\tau M^{-1} \frac{\partial f(t,y)}{\partial y}\right)^{-1}$ ，为此我们使用之前介绍的两个成员函数。

  void Diffusion::implicit_method(const TimeStepping::runge_kutta_method method, 
                                  const unsigned int n_time_steps, 
                                  const double       initial_time, 
                                  const double       final_time) 
  { 
    const double time_step = 
      (final_time - initial_time) / static_cast<double>(n_time_steps); 

    solution = 0.; 
    constraint_matrix.distribute(solution); 

    TimeStepping::ImplicitRungeKutta<Vector<double>> implicit_runge_kutta( 
      method); 
    output_results(initial_time, 0, method); 
    DiscreteTime time(initial_time, final_time, time_step); 
    while (time.is_at_end() == false) 
      { 
        implicit_runge_kutta.evolve_one_time_step( 
          [this](const double time, const Vector<double> &y) { 
            return this->evaluate_diffusion(time, y); 
          }, 
          [this](const double time, const double tau, const Vector<double> &y) { 
            return this->id_minus_tau_J_inverse(time, tau, y); 
          }, 
          time.get_current_time(), 
          time.get_next_step_size(), 
          solution); 
        time.advance_time(); 

        constraint_matrix.distribute(solution); 

        if (time.get_step_number() % 10 == 0) 
          output_results(time.get_current_time(), 
                         time.get_step_number(), 
                         method); 
      } 
  } 

//  @sect4{<code>Diffusion::embedded_explicit_method</code>}  这个函数是嵌入式显式方法的驱动。它需要更多的参数。 

// - coarsen_param：当误差低于阈值时，乘以当前时间步长的系数。 

// - refine_param: 当误差高于阈值时，乘以当前时间步长的系数。 

// - min_delta: 可接受的最小时间步长。 

// - max_delta: 可接受的最大时间步长。 

// - refine_tol：时间步长超过的阈值。 

// - coarsen_tol：阈值，低于该阈值的时间步长将被粗化。

// 嵌入方法使用一个猜测的时间步长。如果使用这个时间步长的误差太大，时间步长将被缩小。如果误差低于阈值，则在下一个时间步长时将尝试更大的时间步长。  <code>delta_t_guess</code> 是由嵌入式方法产生的猜测的时间步长。总之，时间步长有可能以三种方式修改。 

// - 在 TimeStepping::EmbeddedExplicitRungeKutta::evolve_one_time_step(). 内减少或增加时间步长。  

// - 使用计算出的  <code>delta_t_guess</code>  。 

// - 自动调整最后一个时间步长，以确保模拟在  <code>final_time</code>  处精确结束。这种调整是在DiscreteTime实例中处理的。

  unsigned int Diffusion::embedded_explicit_method( 
    const TimeStepping::runge_kutta_method method, 
    const unsigned int                     n_time_steps, 
    const double                           initial_time, 
    const double                           final_time) 
  { 
    const double time_step = 
      (final_time - initial_time) / static_cast<double>(n_time_steps); 
    const double coarsen_param = 1.2; 
    const double refine_param  = 0.8; 
    const double min_delta     = 1e-8; 
    const double max_delta     = 10 * time_step; 
    const double refine_tol    = 1e-1; 
    const double coarsen_tol   = 1e-5; 

    solution = 0.; 
    constraint_matrix.distribute(solution); 

    TimeStepping::EmbeddedExplicitRungeKutta<Vector<double>> 
      embedded_explicit_runge_kutta(method, 
                                    coarsen_param, 
                                    refine_param, 
                                    min_delta, 
                                    max_delta, 
                                    refine_tol, 
                                    coarsen_tol); 
    output_results(initial_time, 0, method); 
    DiscreteTime time(initial_time, final_time, time_step); 
    while (time.is_at_end() == false) 
      { 
        const double new_time = 
          embedded_explicit_runge_kutta.evolve_one_time_step( 
            [this](const double time, const Vector<double> &y) { 
              return this->evaluate_diffusion(time, y); 
            }, 
            time.get_current_time(), 
            time.get_next_step_size(), 
            solution); 
        time.set_next_step_size(new_time - time.get_current_time()); 
        time.advance_time(); 

        constraint_matrix.distribute(solution); 

        if (time.get_step_number() % 10 == 0) 
          output_results(time.get_current_time(), 
                         time.get_step_number(), 
                         method); 

        time.set_desired_next_step_size( 
          embedded_explicit_runge_kutta.get_status().delta_t_guess); 
      } 

    return time.get_step_number(); 
  } 

//  @sect4{<code>Diffusion::run</code>}  

// 下面是该程序的主要功能。在顶部，我们创建网格（一个[0,5]x[0,5]的正方形）并对其进行四次细化，得到一个有16乘16单元的网格，共256个。 然后我们将边界指示器设置为1，用于边界中 $x=0$ 和 $x=5$ 的部分。

  void Diffusion::run() 
  { 
    GridGenerator::hyper_cube(triangulation, 0., 5.); 
    triangulation.refine_global(4); 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          { 
            if ((face->center()[0] == 0.) || (face->center()[0] == 5.)) 
              face->set_boundary_id(1); 
            else 
              face->set_boundary_id(0); 
          } 

// 接下来，我们设置线性系统并为其填充内容，以便在整个时间步进过程中使用它们。

    setup_system(); 

    assemble_system(); 

// 最后，我们使用命名空间TimeStepping中实现的几种Runge-Kutta方法来解决扩散问题，每次都会在结束时输出误差。(正如介绍中所解释的，由于精确解在最后时间为零，所以误差等于数值解，只需取解向量的 $l_2$ 准则即可计算出来。)

    unsigned int       n_steps      = 0; 
    const unsigned int n_time_steps = 200; 
    const double       initial_time = 0.; 
    const double       final_time   = 10.; 

    std::cout << "Explicit methods:" << std::endl; 
    explicit_method(TimeStepping::FORWARD_EULER, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Forward Euler:            error=" << solution.l2_norm() 
              << std::endl; 

    explicit_method(TimeStepping::RK_THIRD_ORDER, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Third order Runge-Kutta:  error=" << solution.l2_norm() 
              << std::endl; 

    explicit_method(TimeStepping::RK_CLASSIC_FOURTH_ORDER, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Fourth order Runge-Kutta: error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << std::endl; 

    std::cout << "Implicit methods:" << std::endl; 
    implicit_method(TimeStepping::BACKWARD_EULER, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Backward Euler:           error=" << solution.l2_norm() 
              << std::endl; 

    implicit_method(TimeStepping::IMPLICIT_MIDPOINT, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Implicit Midpoint:        error=" << solution.l2_norm() 
              << std::endl; 

    implicit_method(TimeStepping::CRANK_NICOLSON, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   Crank-Nicolson:           error=" << solution.l2_norm() 
              << std::endl; 

    implicit_method(TimeStepping::SDIRK_TWO_STAGES, 
                    n_time_steps, 
                    initial_time, 
                    final_time); 
    std::cout << "   SDIRK:                    error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << std::endl; 

    std::cout << "Embedded explicit methods:" << std::endl; 
    n_steps = embedded_explicit_method(TimeStepping::HEUN_EULER, 
                                       n_time_steps, 
                                       initial_time, 
                                       final_time); 
    std::cout << "   Heun-Euler:               error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << "                   steps performed=" << n_steps << std::endl; 

    n_steps = embedded_explicit_method(TimeStepping::BOGACKI_SHAMPINE, 
                                       n_time_steps, 
                                       initial_time, 
                                       final_time); 
    std::cout << "   Bogacki-Shampine:         error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << "                   steps performed=" << n_steps << std::endl; 

    n_steps = embedded_explicit_method(TimeStepping::DOPRI, 
                                       n_time_steps, 
                                       initial_time, 
                                       final_time); 
    std::cout << "   Dopri:                    error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << "                   steps performed=" << n_steps << std::endl; 

    n_steps = embedded_explicit_method(TimeStepping::FEHLBERG, 
                                       n_time_steps, 
                                       initial_time, 
                                       final_time); 
    std::cout << "   Fehlberg:                 error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << "                   steps performed=" << n_steps << std::endl; 

    n_steps = embedded_explicit_method(TimeStepping::CASH_KARP, 
                                       n_time_steps, 
                                       initial_time, 
                                       final_time); 
    std::cout << "   Cash-Karp:                error=" << solution.l2_norm() 
              << std::endl; 
    std::cout << "                   steps performed=" << n_steps << std::endl; 
  } 
} // namespace Step52 

//  @sect3{The <code>main()</code> function}  

// 下面的 <code>main</code> 函数与前面的例子类似，不需要注释。

int main() 
{ 
  try 
    { 
      Step52::Diffusion diffusion; 
      diffusion.run(); 
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


