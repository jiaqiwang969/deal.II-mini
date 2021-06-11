CCTest_file/step-41.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2011 - 2021 by the deal.II authors 
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
 * Authors: Joerg Frohne, Texas A&M University and 
 *                        University of Siegen, 2011, 2012 
 *          Wolfgang Bangerth, Texas A&M University, 2012 
 */ 


// @sect3{Include files}  

// 像往常一样，在开始的时候，我们把所有我们需要的头文件都包含在这里。除了为Trilinos库提供接口的各种文件外，没有什么意外。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/index_set.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_precondition.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <iostream> 

namespace Step41 
{ 
  using namespace dealii; 
// @sect3{The <code>ObstacleProblem</code> class template}  

// 该类提供了描述障碍问题所需的所有函数和变量。它与我们在 step-4 中要做的事情很接近，所以相对简单。唯一真正的新组件是计算主动集合的update_solution_and_constraints函数和一些描述线性系统原始（无约束）形式所需的变量（ <code>complete_system_matrix</code> 和 <code>complete_system_rhs</code> ），以及主动集合本身和主动集合公式中用于缩放拉格朗日乘数的质量矩阵 $B$ 的对角线。其余的内容与 step-4 相同。

  template <int dim> 
  class ObstacleProblem 
  { 
  public: 
    ObstacleProblem(); 
    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void 
         assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix); 
    void update_solution_and_constraints(); 
    void solve(); 
    void output_results(const unsigned int iteration) const; 

    Triangulation<dim>        triangulation; 
    FE_Q<dim>                 fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 
    IndexSet                  active_set; 

    TrilinosWrappers::SparseMatrix system_matrix; 
    TrilinosWrappers::SparseMatrix complete_system_matrix; 

    TrilinosWrappers::MPI::Vector solution; 
    TrilinosWrappers::MPI::Vector system_rhs; 
    TrilinosWrappers::MPI::Vector complete_system_rhs; 
    TrilinosWrappers::MPI::Vector diagonal_of_mass_matrix; 
    TrilinosWrappers::MPI::Vector contact_force; 
  }; 
// @sect3{Right hand side, boundary values, and the obstacle}  

// 在下文中，我们定义了描述右侧函数、Dirichlet边界值以及作为 $\mathbf x$ 函数的障碍物高度的类。在这三种情况下，我们都从函数 @<dim@>, 派生出这些类，尽管在 <code>RightHandSide</code> 和 <code>Obstacle</code> 的情况下，这更多的是出于惯例而非必要，因为我们从未将此类对象传递给库。在任何情况下，鉴于我们选择了 $f=-10$  ,  $u|_{\partial\Omega}=0$  ...，右手和边界值类的定义是显而易见的。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      AssertIndexRange(component, 1); 

      return -10; 
    } 
  }; 

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & /*p*/, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      AssertIndexRange(component, 1); 

      return 0; 
    } 
  }; 

// 我们用一个级联的障碍物来描述障碍物的功能（想想看：楼梯的阶梯）。

  template <int dim> 
  class Obstacle : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override 
    { 
      (void)component; 
      Assert(component == 0, ExcIndexRange(component, 0, 1)); 

      if (p(0) < -0.5) 
        return -0.2; 
      else if (p(0) >= -0.5 && p(0) < 0.0) 
        return -0.4; 
      else if (p(0) >= 0.0 && p(0) < 0.5) 
        return -0.6; 
      else 
        return -0.8; 
    } 
  }; 

//  @sect3{Implementation of the <code>ObstacleProblem</code> class}  
// @sect4{ObstacleProblem::ObstacleProblem}  

// 对每个看过前几个教程程序的人来说，构造函数是完全显而易见的。

  template <int dim> 
  ObstacleProblem<dim>::ObstacleProblem() 
    : fe(1) 
    , dof_handler(triangulation) 
  {} 
// @sect4{ObstacleProblem::make_grid}  

// 我们在二维的正方形 $[-1,1]\times [-1,1]$ 上解决我们的障碍物问题。因此这个函数只是设置了一个最简单的网格。

  template <int dim> 
  void ObstacleProblem<dim>::make_grid() 
  { 
    GridGenerator::hyper_cube(triangulation, -1, 1); 
    triangulation.refine_global(7); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Total number of cells: " << triangulation.n_cells() 
              << std::endl; 
  } 
// @sect4{ObstacleProblem::setup_system}  

// 在这个值得注意的第一个函数中，我们设置了自由度处理程序，调整了向量和矩阵的大小，并处理了约束。最初，约束条件当然只是由边界值给出的，所以我们在函数的顶部对它们进行插值。

  template <int dim> 
  void ObstacleProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    active_set.set_size(dof_handler.n_dofs()); 

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl 
              << std::endl; 

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             BoundaryValues<dim>(), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 

    system_matrix.reinit(dsp); 
    complete_system_matrix.reinit(dsp); 

    IndexSet solution_index_set = dof_handler.locally_owned_dofs(); 
    solution.reinit(solution_index_set, MPI_COMM_WORLD); 
    system_rhs.reinit(solution_index_set, MPI_COMM_WORLD); 
    complete_system_rhs.reinit(solution_index_set, MPI_COMM_WORLD); 
    contact_force.reinit(solution_index_set, MPI_COMM_WORLD); 

// 这里唯一要做的事情是计算 $B$ 矩阵中的因子，该矩阵用于缩放残差。正如在介绍中所讨论的，我们将使用一个小技巧来使这个质量矩阵成为对角线，在下文中，首先将所有这些计算成一个矩阵，然后提取对角线元素供以后使用。

    TrilinosWrappers::SparseMatrix mass_matrix; 
    mass_matrix.reinit(dsp); 
    assemble_mass_matrix_diagonal(mass_matrix); 
    diagonal_of_mass_matrix.reinit(solution_index_set); 
    for (unsigned int j = 0; j < solution.size(); j++) 
      diagonal_of_mass_matrix(j) = mass_matrix.diag_element(j); 
  } 
// @sect4{ObstacleProblem::assemble_system}  

// 这个函数一次就把系统矩阵和右手边集合起来，并把约束条件（由于活动集以及来自边界值）应用到我们的系统中。否则，它在功能上等同于例如  step-4  中的相应函数。

  template <int dim> 
  void ObstacleProblem<dim>::assemble_system() 
  { 
    std::cout << "   Assembling system..." << std::endl; 

    system_matrix = 0; 
    system_rhs    = 0; 

    const QGauss<dim>  quadrature_formula(fe.degree + 1); 
    RightHandSide<dim> right_hand_side; 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        cell_matrix = 0; 
        cell_rhs    = 0; 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                cell_matrix(i, j) += 
                  (fe_values.shape_grad(i, q_point) * 
                   fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point)); 

              cell_rhs(i) += 
                (fe_values.shape_value(i, q_point) * 
                 right_hand_side.value(fe_values.quadrature_point(q_point)) * 
                 fe_values.JxW(q_point)); 
            } 

        cell->get_dof_indices(local_dof_indices); 

        constraints.distribute_local_to_global(cell_matrix, 
                                               cell_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               system_rhs, 
                                               true); 
      } 
  } 

//  @sect4{ObstacleProblem::assemble_mass_matrix_diagonal}  

// 下一个函数用于计算对角线质量矩阵 $B$ ，用于在主动集方法中缩放变量。正如介绍中所讨论的，我们通过选择正交的梯形规则来获得质量矩阵的对角线。这样一来，我们就不再需要在正交点、指数 $i$ 和指数 $j$ 上进行三重循环，而是可以直接使用双重循环。考虑到我们在以前的许多教程程序中讨论过的内容，该函数的其余部分是显而易见的。

// 注意在调用这个函数的时候，约束对象只包含边界值约束；因此我们在最后的复制-本地-全局步骤中不必注意保留矩阵项的值，这些项以后可能会受到活动集的约束。

// 还需要注意的是，只有在我们拥有 $Q_1$ 元素的情况下，使用梯形规则的技巧才有效。对于更高阶的元素，我们需要使用一个正交公式，在有限元的所有支持点都有正交点。构建这样一个正交公式其实并不难，但不是这里的重点，所以我们只是在函数的顶部断言我们对有限元的隐含假设实际上得到了满足。

  template <int dim> 
  void ObstacleProblem<dim>::assemble_mass_matrix_diagonal( 
    TrilinosWrappers::SparseMatrix &mass_matrix) 
  { 
    Assert(fe.degree == 1, ExcNotImplemented()); 

    const QTrapezoid<dim> quadrature_formula; 
    FEValues<dim>         fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

 
      { 
        fe_values.reinit(cell); 
        cell_matrix = 0; 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            cell_matrix(i, i) += 
              (fe_values.shape_value(i, q_point) * 
               fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 

        cell->get_dof_indices(local_dof_indices); 

        constraints.distribute_local_to_global(cell_matrix, 
                                               local_dof_indices, 
                                               mass_matrix); 
      } 
  } 
// @sect4{ObstacleProblem::update_solution_and_constraints}  

// 在某种意义上，这是本程序的核心功能。 它更新了介绍中所讨论的受限自由度的活动集，并从中计算出一个AffineConstraints对象，然后可以用来在下一次迭代的解中消除受限自由度。同时，我们将解决方案的受限自由度设置为正确的值，即障碍物的高度。

// 从根本上说，这个函数是相当简单的。我们必须在所有自由度上循环，并检查函数 $\Lambda^k_i + c([BU^k]_i - G_i) = \Lambda^k_i + cB_i(U^k_i - [g_h]_i)$ 的符号，因为在我们的例子中 $G_i = B_i[g_h]_i$  。为此，我们使用介绍中给出的公式，通过该公式我们可以计算出拉格朗日乘数，作为原始线性系统的残差（通过变量 <code>complete_system_matrix</code> and <code>complete_system_rhs</code> 给出。在这个函数的顶部，我们使用一个属于矩阵类的函数来计算这个残差。

  template <int dim> 
  void ObstacleProblem<dim>::update_solution_and_constraints() 
  { 
    std::cout << "   Updating active set..." << std::endl; 

    const double penalty_parameter = 100.0; 

    TrilinosWrappers::MPI::Vector lambda( 
      complete_index_set(dof_handler.n_dofs())); 
    complete_system_matrix.residual(lambda, solution, complete_system_rhs); 

// 计算 contact_force[i] =

// - lambda[i] * diagonal_of_mass_matrix[i]。

    contact_force = lambda; 
    contact_force.scale(diagonal_of_mass_matrix); 
    contact_force *= -1; 

// 下一步是重置活动集和约束对象，并在所有自由度上开始循环。由于我们不能只是在解向量的所有元素上循环，所以这变得稍微复杂了一些，因为我们没有办法找出一个自由度与哪个位置相关；但是，我们需要这个位置来测试一个自由度的位移是大于还是小于这个位置的障碍物高度。

// 我们通过在所有单元和定义在每个单元上的DoF上循环来解决这个问题。我们在这里使用一个 $Q_1$ 函数来描述位移，对于该函数，自由度总是位于单元格的顶点上；因此，我们可以通过询问顶点来获得每个自由度的索引及其位置。另一方面，这显然对高阶元素不起作用，因此我们添加了一个断言，确保我们只处理所有自由度都位于顶点的元素，以避免万一有人想玩增加解的多项式程度时用非功能性代码绊倒自己。

// 循环单元格而不是自由度的代价是我们可能会多次遇到一些自由度，即每次我们访问与给定顶点相邻的一个单元格时。因此，我们必须跟踪我们已经接触过的顶点和尚未接触的顶点。我们通过使用一个标志数组来做到这一点  <code>dof_touched</code>  。

    constraints.clear(); 
    active_set.clear(); 

    const Obstacle<dim> obstacle; 
    std::vector<bool>   dof_touched(dof_handler.n_dofs(), false); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto v : cell->vertex_indices()) 
        { 
          Assert(dof_handler.get_fe().n_dofs_per_cell() == cell->n_vertices(), 
                 ExcNotImplemented()); 

          const unsigned int dof_index = cell->vertex_dof_index(v, 0); 

          if (dof_touched[dof_index] == false) 
            dof_touched[dof_index] = true; 
          else 
            continue; 

// 现在我们知道我们还没有触及这个DoF，让我们得到那里的位移函数的值以及障碍函数的值，并使用这个来决定当前DoF是否属于活动集。为此，我们使用上面和介绍中给出的函数。

// 如果我们决定该DoF应该是活动集的一部分，我们将其索引添加到活动集中，在AffineConstraints对象中引入一个不均匀的平等约束，并将解的值重置为障碍物的高度。最后，系统的非接触部分的残差作为一个额外的控制（残差等于剩余的、未计算的力，在接触区之外应该为零），所以我们把残差向量的分量（即拉格朗日乘数lambda）清零，这些分量对应于身体接触的区域；在所有单元的循环结束时，残差将因此只包括非接触区的残差。我们在循环结束后输出这个残差的准则和活动集的大小。

          const double obstacle_value = obstacle.value(cell->vertex(v)); 
          const double solution_value = solution(dof_index); 

          if (lambda(dof_index) + penalty_parameter * 
                                    diagonal_of_mass_matrix(dof_index) * 
                                    (solution_value - obstacle_value) < 
              0) 
            { 
              active_set.add_index(dof_index); 
              constraints.add_line(dof_index); 
              constraints.set_inhomogeneity(dof_index, obstacle_value); 

              solution(dof_index) = obstacle_value; 

              lambda(dof_index) = 0; 
            } 
        } 
    std::cout << "      Size of active set: " << active_set.n_elements() 
              << std::endl; 

    std::cout << "   Residual of the non-contact part of the system: " 
              << lambda.l2_norm() << std::endl; 

// 在最后一步中，我们将迄今为止从活动集合中得到的对DoF的约束加入到那些由Dirichlet边界值产生的约束中，并关闭约束对象。

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             BoundaryValues<dim>(), 
                                             constraints); 
    constraints.close(); 
  } 
// @sect4{ObstacleProblem::solve}  

// 关于求解函数，其实没有什么可说的。在牛顿方法的背景下，我们通常对非常高的精度不感兴趣（为什么要求一个高度精确的线性问题的解，而我们知道它只能给我们一个非线性问题的近似解），所以我们使用ReductionControl类，当达到一个绝对公差（为此我们选择 $10^{-12}$ ）或者当残差减少一定的系数（这里是 $10^{-3}$ ）时停止反复运算。

  template <int dim> 
  void ObstacleProblem<dim>::solve() 
  { 
    std::cout << "   Solving system..." << std::endl; 

    ReductionControl                        reduction_control(100, 1e-12, 1e-3); 
    SolverCG<TrilinosWrappers::MPI::Vector> solver(reduction_control); 
    TrilinosWrappers::PreconditionAMG       precondition; 
    precondition.initialize(system_matrix); 

    solver.solve(system_matrix, solution, system_rhs, precondition); 
    constraints.distribute(solution); 

    std::cout << "      Error: " << reduction_control.initial_value() << " -> " 
              << reduction_control.last_value() << " in " 
              << reduction_control.last_step() << " CG iterations." 
              << std::endl; 
  } 
// @sect4{ObstacleProblem::output_results}  

// 我们使用vtk-format进行输出。 该文件包含位移和活动集的数字表示。

  template <int dim> 
  void ObstacleProblem<dim>::output_results(const unsigned int iteration) const 
  { 
    std::cout << "   Writing graphical output..." << std::endl; 

    TrilinosWrappers::MPI::Vector active_set_vector( 
      dof_handler.locally_owned_dofs(), MPI_COMM_WORLD); 
    for (const auto index : active_set) 
      active_set_vector[index] = 1.; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "displacement"); 
    data_out.add_data_vector(active_set_vector, "active_set"); 
    data_out.add_data_vector(contact_force, "lambda"); 

    data_out.build_patches(); 

    std::ofstream output_vtk("output_" + 
                             Utilities::int_to_string(iteration, 3) + ".vtk"); 
    data_out.write_vtk(output_vtk); 
  } 

//  @sect4{ObstacleProblem::run}  

// 这是一个对所有事情都有最高级别控制的函数。 它并不长，而且事实上相当直接：在主动集方法的每一次迭代中，我们都要组装线性系统，求解它，更新主动集并将解投射回可行集，然后输出结果。只要主动集在前一次迭代中没有变化，迭代就会终止。

// 唯一比较棘手的部分是，我们必须在第一次迭代组装好线性系统（即矩阵和右手边）后保存它。原因是这是唯一一个我们可以在没有任何接触约束的情况下访问线性系统的步骤。我们需要这个来计算其他迭代中解决方案的残差，但是在其他迭代中，我们形成的线性系统中对应于约束自由度的行和列都被消除了，因此我们不能再访问原始方程的全部残差。

  template <int dim> 
  void ObstacleProblem<dim>::run() 
  { 
    make_grid(); 
    setup_system(); 

    IndexSet active_set_old(active_set); 
    for (unsigned int iteration = 0; iteration <= solution.size(); ++iteration) 
      { 
        std::cout << "Newton iteration " << iteration << std::endl; 

        assemble_system(); 

        if (iteration == 0) 
          { 
            complete_system_matrix.copy_from(system_matrix); 
            complete_system_rhs = system_rhs; 
          } 

        solve(); 
        update_solution_and_constraints(); 
        output_results(iteration); 

        if (active_set == active_set_old) 
          break; 

        active_set_old = active_set; 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step41 
// @sect3{The <code>main</code> function}  

// 这就是主函数。它遵循所有其他主函数的模式。调用初始化MPI是因为我们在这个程序中建立线性求解器的Trilinos库需要它。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step41; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 

// 这个程序只能在串行中运行。否则，将抛出一个异常。

      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
                  ExcMessage( 
                    "This program can only be run in serial, use ./step-41")); 

      ObstacleProblem<2> obstacle_problem; 
      obstacle_problem.run(); 
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


