

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Colorado State University, 2021. 
 * Based on step-15 by Sven Wetterauer, University of Heidelberg, 2012. 
 */ 


// @sect3{Include files}  

// 这个程序开始时和其他大多数程序一样，有众所周知的包含文件。与 step-15 程序相比，我们在这里所做的大部分工作都是从该程序中复制的，唯一不同的是包括头文件，我们从该文件中导入了SparseDirectUMFPACK类和KINSOL的实际接口。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparse_direct.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_accessor.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/solution_transfer.h> 

#include <deal.II/sundials/kinsol.h> 

#include <fstream> 
#include <iostream> 

namespace Step77 
{ 
  using namespace dealii; 
// @sect3{The <code>MinimalSurfaceProblem</code> class template}  

// 同样地，这个程序的主类基本上是  step-15  中的一个副本。然而，该类确实将雅各布（系统）矩阵（以及使用直接求解器对其进行因式分解）和残差的计算分成了不同的函数，原因已在介绍中列出。出于同样的原因，该类也有一个指向雅各布矩阵因式分解的指针，该指针在我们每次更新雅各布矩阵时被重置。

// （如果你想知道为什么程序对雅各布矩阵使用直接对象，而对因式分解使用指针。每次KINSOL要求更新雅各布矩阵时，我们可以简单地写`jacobian_matrix=0;`将其重置为一个空矩阵，然后我们可以再次填充。另一方面，SparseDirectUMFPACK类没有办法扔掉它的内容或用新的因式分解来替换它，所以我们使用一个指针。我们只是扔掉整个对象，并在我们有新的雅各布矩阵需要分解时创建一个新的对象。)

// 最后，该类有一个定时器变量，我们将用它来评估程序的不同部分需要多长时间，这样我们就可以评估KINSOL的不重建矩阵及其因式分解的倾向是否合理。我们将在下面的 "结果 "部分讨论这个问题。

  template <int dim> 
  class MinimalSurfaceProblem 
  { 
  public: 
    MinimalSurfaceProblem(); 
    void run(); 

  private: 
    void setup_system(const bool initial_step); 
    void solve(const Vector<double> &rhs, 
               Vector<double> &      solution, 
               const double          tolerance); 
    void refine_mesh(); 
    void output_results(const unsigned int refinement_cycle); 
    void set_boundary_values(); 
    void compute_and_factorize_jacobian(const Vector<double> &evaluation_point); 
    void compute_residual(const Vector<double> &evaluation_point, 
                          Vector<double> &      residual); 

    Triangulation<dim> triangulation; 

    DoFHandler<dim> dof_handler; 
    FE_Q<dim>       fe; 

    AffineConstraints<double> hanging_node_constraints; 

    SparsityPattern                      sparsity_pattern; 
    SparseMatrix<double>                 jacobian_matrix; 
    std::unique_ptr<SparseDirectUMFPACK> jacobian_matrix_factorization; 

    Vector<double> current_solution; 

    TimerOutput computing_timer; 
  }; 

//  @sect3{Boundary condition}  

// 实现边界值的类是对  step-15  的复制。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> &p, 
                                    const unsigned int /*component*/) const 
  { 
    return std::sin(2 * numbers::PI * (p[0] + p[1])); 
  } 
// @sect3{The <code>MinimalSurfaceProblem</code> class implementation}  
// @sect4{Constructor and set up functions}  

// 下面的几个函数也基本上是复制了 step-15 已经做的事情，所以没有什么可讨论的。

  template <int dim> 
  MinimalSurfaceProblem<dim>::MinimalSurfaceProblem() 
    : dof_handler(triangulation) 
    , fe(1) 
    , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times) 
  {} 

  template <int dim> 
  void MinimalSurfaceProblem<dim>::setup_system(const bool initial_step) 
  { 
    TimerOutput::Scope t(computing_timer, "set up"); 

    if (initial_step) 
      { 
        dof_handler.distribute_dofs(fe); 
        current_solution.reinit(dof_handler.n_dofs()); 

        hanging_node_constraints.clear(); 
        DoFTools::make_hanging_node_constraints(dof_handler, 
                                                hanging_node_constraints); 
        hanging_node_constraints.close(); 
      } 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    hanging_node_constraints.condense(dsp); 

    sparsity_pattern.copy_from(dsp); 
    jacobian_matrix.reinit(sparsity_pattern); 
    jacobian_matrix_factorization.reset(); 
  } 

//  @sect4{Assembling and factorizing the Jacobian matrix}  

// 然后，下面的函数负责对雅各布矩阵进行组装和因子化。该函数的前半部分实质上是 step-15 的`assemble_system()`函数，只是它没有处理同时形成右手边的向量（即残差），因为我们并不总是要同时做这些操作。

// 我们把整个装配功能放在一个由大括号包围的代码块中，这样我们就可以用一个 TimerOutput::Scope 变量来衡量在这个代码块中花费了多少时间，不包括在这个函数中发生在匹配的闭合括号`}`之后的一切。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::compute_and_factorize_jacobian( 
    const Vector<double> &evaluation_point) 
  { 
    { 
      TimerOutput::Scope t(computing_timer, "assembling the Jacobian"); 

      std::cout << "  Computing Jacobian matrix" << std::endl; 

      const QGauss<dim> quadrature_formula(fe.degree + 1); 

      jacobian_matrix = 0; 

      FEValues<dim> fe_values(fe, 
                              quadrature_formula, 
                              update_gradients | update_quadrature_points | 
                                update_JxW_values); 

      const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
      const unsigned int n_q_points    = quadrature_formula.size(); 

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 

      std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points); 

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          cell_matrix = 0; 

          fe_values.reinit(cell); 

          fe_values.get_function_gradients(evaluation_point, 
                                           evaluation_point_gradients); 

          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 
              const double coeff = 
                1.0 / std::sqrt(1 + evaluation_point_gradients[q] * 
                                      evaluation_point_gradients[q]); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    cell_matrix(i, j) += 
                      (((fe_values.shape_grad(i, q)    // ((\nabla \phi_i 
                         * coeff                       //   * a_n 
                         * fe_values.shape_grad(j, q)) //   * \nabla \phi_j) 
                        -                              //  - 
                        (fe_values.shape_grad(i, q)    //  (\nabla \phi_i 
                         * coeff * coeff * coeff       //   * a_n^3 
                         * 
                         (fe_values.shape_grad(j, q)       //   * (\nabla \phi_j 
                          * evaluation_point_gradients[q]) //      * \nabla u_n) 
                         * evaluation_point_gradients[q])) //   * \nabla u_n))) 
                       * fe_values.JxW(q));                // * dx 
                } 
            } 

          cell->get_dof_indices(local_dof_indices); 
          hanging_node_constraints.distribute_local_to_global(cell_matrix, 
                                                              local_dof_indices, 
                                                              jacobian_matrix); 
        } 

      std::map<types::global_dof_index, double> boundary_values; 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               Functions::ZeroFunction<dim>(), 
                                               boundary_values); 
      Vector<double> dummy_solution(dof_handler.n_dofs()); 
      Vector<double> dummy_rhs(dof_handler.n_dofs()); 
      MatrixTools::apply_boundary_values(boundary_values, 
                                         jacobian_matrix, 
                                         dummy_solution, 
                                         dummy_rhs); 
    } 

// 该函数的后半部分是对计算出的矩阵进行因数分解。为此，我们首先创建一个新的SparseDirectUMFPACK对象，并将其分配给成员变量`jacobian_matrix_factorization`，同时销毁该指针之前指向的任何对象（如果有）。然后我们告诉该对象对雅各布系数进行分解。

// 如上所述，我们把这段代码放在大括号里，用一个计时器来评估这部分程序所需的时间。

// (严格来说，我们在这里完成后实际上不再需要矩阵了，我们可以把矩阵对象扔掉。一个旨在提高内存效率的代码会这样做，并且只在这个函数中创建矩阵对象，而不是作为周围类的成员变量。我们在这里省略了这一步，因为使用与以前的教程程序相同的编码风格可以培养对通用风格的熟悉，并有助于使这些教程程序更容易阅读)。

    { 
      TimerOutput::Scope t(computing_timer, "factorizing the Jacobian"); 

      std::cout << "  Factorizing Jacobian matrix" << std::endl; 

      jacobian_matrix_factorization = std::make_unique<SparseDirectUMFPACK>(); 
      jacobian_matrix_factorization->factorize(jacobian_matrix); 
    } 
  } 

//  @sect4{Computing the residual vector}  

// `assemble_system()`在 step-15 中用来做的第二部分是计算残差向量，也就是牛顿线性系统的右手向量。我们把这一点从前面的函数中分解出来，但如果你理解了 step-15 中`assemble_system()`的作用，下面的函数就会很容易理解。然而，重要的是，我们需要计算的残差不是围绕当前解向量线性化的，而是我们从KINSOL得到的任何东西。这对于诸如直线搜索这样的操作是必要的，我们想知道在不同的 $\alpha_k$ 值下，残差 $F(U^k + \alpha_k \delta U^K)$ 是多少；在这些情况下，KINSOL只是给我们函数 $F$ 的参数，然后我们在这时计算残差 $F(\cdot)$ 。

// 该函数在最后打印出如此计算的残差的规范，作为我们跟踪程序进展的一种方式。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::compute_residual( 
    const Vector<double> &evaluation_point, 
    Vector<double> &      residual) 
  { 
    TimerOutput::Scope t(computing_timer, "assembling the residual"); 

    std::cout << "  Computing residual vector..." << std::flush; 

    const QGauss<dim> quadrature_formula(fe.degree + 1); 
    FEValues<dim>     fe_values(fe, 
                            quadrature_formula, 
                            update_gradients | update_quadrature_points | 
                              update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    Vector<double>              cell_residual(dofs_per_cell); 
    std::vector<Tensor<1, dim>> evaluation_point_gradients(n_q_points); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        cell_residual = 0; 
        fe_values.reinit(cell); 

        fe_values.get_function_gradients(evaluation_point, 
                                         evaluation_point_gradients); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            const double coeff = 
              1.0 / std::sqrt(1 + evaluation_point_gradients[q] * 
                                    evaluation_point_gradients[q]); 

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_residual(i) = (fe_values.shape_grad(i, q) // \nabla \phi_i 
                                  * coeff                    // * a_n 
                                  * evaluation_point_gradients[q] // * u_n 
                                  * fe_values.JxW(q));            // * dx 
          } 

        cell->get_dof_indices(local_dof_indices); 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          residual(local_dof_indices[i]) += cell_residual(i); 
      } 

    hanging_node_constraints.condense(residual); 

    for (const types::global_dof_index i : 
         DoFTools::extract_boundary_dofs(dof_handler)) 
      residual(i) = 0; 

    for (const types::global_dof_index i : 
         DoFTools::extract_hanging_node_dofs(dof_handler)) 
      residual(i) = 0; 

    std::cout << " norm=" << residual.l2_norm() << std::endl; 
  } 

//  @sect4{Solving linear systems with the Jacobian matrix}  

// 接下来是实现用雅各布矩阵解线性系统的函数。由于我们在建立矩阵时已经对矩阵进行了因式分解，所以解决线性系统的方法就是将逆矩阵应用于给定的右侧向量。这就是我们在这里使用的 SparseDirectUMFPACK::vmult() 函数的作用。在这之后，我们必须确保我们也能解决解向量中的悬空节点的值，而这是用 AffineConstraints::distribute(). 来完成的。

// 该函数需要一个额外的，但未使用的参数`tolerance`，它表示我们必须解决线性系统的精确程度。这个参数的含义在介绍中结合 "Eisenstat Walker技巧 "进行了讨论，但由于我们使用的是直接求解器而不是迭代求解器，所以我们并没有利用这个机会只求解线性系统的不精确性。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::solve(const Vector<double> &rhs, 
                                         Vector<double> &      solution, 
                                         const double /*tolerance*/) 
  { 
    TimerOutput::Scope t(computing_timer, "linear system solve"); 

    std::cout << "  Solving linear system" << std::endl; 

    jacobian_matrix_factorization->vmult(solution, rhs); 

    hanging_node_constraints.distribute(solution); 
  } 

//  @sect4{Refining the mesh, setting boundary values, and generating graphical output}  

// 以下三个函数又是对  step-15  中的函数的简单复制。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::refine_mesh() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      current_solution, 
      estimated_error_per_cell); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.03); 

    triangulation.prepare_coarsening_and_refinement(); 

    SolutionTransfer<dim> solution_transfer(dof_handler); 
    solution_transfer.prepare_for_coarsening_and_refinement(current_solution); 

    triangulation.execute_coarsening_and_refinement(); 

    dof_handler.distribute_dofs(fe); 

    Vector<double> tmp(dof_handler.n_dofs()); 
    solution_transfer.interpolate(current_solution, tmp); 
    current_solution = std::move(tmp); 

    hanging_node_constraints.clear(); 

    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 

    hanging_node_constraints.distribute(current_solution); 

    set_boundary_values(); 

    setup_system(/*initial_step=*/false); 
  } 

  template <int dim> 
  void MinimalSurfaceProblem<dim>::set_boundary_values() 
  { 
    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             BoundaryValues<dim>(), 
                                             boundary_values); 
    for (const auto &boundary_value : boundary_values) 
      current_solution(boundary_value.first) = boundary_value.second; 

    hanging_node_constraints.distribute(current_solution); 
  } 

  template <int dim> 
  void MinimalSurfaceProblem<dim>::output_results( 
    const unsigned int refinement_cycle) 
  { 
    TimerOutput::Scope t(computing_timer, "graphical output"); 

    DataOut<dim> data_out; 

 
    data_out.add_data_vector(current_solution, "solution"); 
    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtu"; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect4{The run() function and the overall logic of the program}  

// 这个程序中唯一**有趣的函数是驱动整个算法的函数，即从一个粗大的网格开始，做一些网格细化循环，并在每个网格上使用KINSOL来寻找我们从这个网格上离散化得到的非线性代数方程的解。上面的`refine_mesh()`函数可以确保一个网格上的解被用作下一个网格的起始猜测。我们还使用一个TimerOutput对象来测量每个网格上的每一次操作所花费的时间，并在每个周期开始时重置该计时器。

// 正如在介绍中所讨论的，没有必要特别精确地解决粗略网格上的问题，因为这些问题只能作为下一个网格的起始猜测来解决。因此，我们将在 $k$ 个网格细化周期中使用 $\tau=10^{-3} \frac{1}{10^k}$ 的目标公差。

// 所有这些都在这个函数的第一部分进行了编码。

  template <int dim> 
  void MinimalSurfaceProblem<dim>::run() 
  { 
    GridGenerator::hyper_ball(triangulation); 
    triangulation.refine_global(2); 

    setup_system(/*initial_step=*/true); 
    set_boundary_values(); 

    for (unsigned int refinement_cycle = 0; refinement_cycle < 6; 
         ++refinement_cycle) 
      { 
        computing_timer.reset(); 
        std::cout << "Mesh refinement step " << refinement_cycle << std::endl; 

        if (refinement_cycle != 0) 
          refine_mesh(); 

        const double target_tolerance = 1e-3 * std::pow(0.1, refinement_cycle); 
        std::cout << "  Target_tolerance: " << target_tolerance << std::endl 
                  << std::endl; 

// 这就是有趣的开始。在顶部，我们创建了KINSOL求解器对象，并给它提供了一个对象，该对象编码了一些额外的具体情况（其中我们只改变了我们想要达到的非线性容忍度；但你可能想看看 SUNDIALS::KINSOL::AdditionalData 类有哪些其他成员，并与它们一起玩）。

        { 
          typename SUNDIALS::KINSOL<Vector<double>>::AdditionalData 
            additional_data; 
          additional_data.function_tolerance = target_tolerance; 

          SUNDIALS::KINSOL<Vector<double>> nonlinear_solver(additional_data); 

// 然后，我们必须描述在介绍中已经提到的操作。从本质上讲，我们必须教KINSOL如何(i)将一个向量调整到正确的大小，(ii)计算残差向量，(iii)计算雅各布矩阵（在这期间我们也计算其因式分解），以及(iv)用雅各布矩阵解一个线性系统。

// 所有这四种操作都由 SUNDIALS::KINSOL 类的成员变量表示，这些成员变量的类型是 `std::function`, ，即它们是我们可以分配给一个函数的指针的对象，或者像我们在这里做的那样，一个 "lambda函数"，它接受相应的参数并返回相应的信息。按照惯例，KINSOL希望做一些不重要的事情的函数返回一个整数，其中0表示成功。事实证明，我们只需用25行代码就可以完成所有这些工作。

// 如果你不知道什么是 "lambda函数"，可以看看 step-12 或[wikipedia页面](https:en.wikipedia.org/wiki/Anonymous_function)关于这个问题。lambda函数的想法是，人们想用一组参数来定义一个函数，但(i)不使它成为一个命名的函数，因为通常情况下，该函数只在一个地方使用，似乎没有必要给它一个全局名称；(ii)该函数可以访问存在于定义它的地方的一些变量，包括成员变量。lambda函数的语法很笨拙，但最终还是很有用的）。)

// 在代码块的最后，我们告诉KINSOL去工作，解决我们的问题。从'residual'、'setup_jacobian'和'solve_jacobian_system'函数中调用的成员函数将向屏幕打印输出，使我们能够跟踪程序的进展情况。

          nonlinear_solver.reinit_vector = [&](Vector<double> &x) { 
            x.reinit(dof_handler.n_dofs()); 
          }; 

          nonlinear_solver.residual = 
            [&](const Vector<double> &evaluation_point, 
                Vector<double> &      residual) { 
              compute_residual(evaluation_point, residual); 

              return 0; 
            }; 

          nonlinear_solver.setup_jacobian = 
            [&](const Vector<double> &current_u, 
                const Vector<double> & /*current_f*/) { 
              compute_and_factorize_jacobian(current_u); 

              return 0; 
            }; 

          nonlinear_solver.solve_with_jacobian = [&](const Vector<double> &rhs, 
                                                     Vector<double> &      dst, 
                                                     const double tolerance) { 
            this->solve(rhs, dst, tolerance); 

            return 0; 
          }; 

          nonlinear_solver.solve(current_solution); 
        } 

// 剩下的就只是内务整理了。将数据写入文件，以便进行可视化，并显示收集到的时间摘要，以便我们可以解释每个操作花了多长时间，执行的频率如何，等等。

        output_results(refinement_cycle); 

        computing_timer.print_summary(); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step77 

int main() 
{ 
  try 
    { 
      using namespace Step77; 

      MinimalSurfaceProblem<2> laplace_problem_2d; 
      laplace_problem_2d.run(); 
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

