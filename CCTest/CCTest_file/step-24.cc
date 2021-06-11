CCTest_file/step-24.cc

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
 * Author: Xing Jin, Wolfgang Bangerth, Texas A&M University, 2006 
 */ 


// @sect3{Include files}  

// 以下内容之前都已经介绍过了。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/vector_tools.h> 

#include <fstream> 
#include <iostream> 

// 这是唯一一个新的。我们将需要一个定义在GridTools命名空间的库函数，用来计算最小的单元格直径。

#include <deal.II/grid/grid_tools.h> 

// 最后一步和以前所有的程序一样。

namespace Step24 
{ 
  using namespace dealii; 
// @sect3{The "forward problem" class template}  

// 主类的第一部分与 step-23 中的内容完全一致（除了名字）。

  template <int dim> 
  class TATForwardProblem 
  { 
  public: 
    TATForwardProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void solve_p(); 
    void solve_v(); 
    void output_results() const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 
    SparseMatrix<double> mass_matrix; 
    SparseMatrix<double> laplace_matrix; 

    Vector<double> solution_p, solution_v; 
    Vector<double> old_solution_p, old_solution_v; 
    Vector<double> system_rhs_p, system_rhs_v; 

    double       time_step, time; 
    unsigned int timestep_number; 
    const double theta; 

// 下面是新的内容：首先，我们需要从吸收边界条件出来的那个边界质量矩阵 $B$ 。同样，由于这次我们考虑的是一个现实的介质，我们必须有一个衡量波速的标准 $c_0$ ，它将进入所有与拉普拉斯矩阵（我们仍然定义为 $(\nabla \phi_i,\nabla \phi_j)$ ）有关的公式。

    SparseMatrix<double> boundary_matrix; 
    const double         wave_speed; 

// 我们必须注意的最后一件事是，我们想在一定数量的检测器位置评估解决方案。我们需要一个数组来保存这些位置，在这里声明并在构造函数中填充。

    std::vector<Point<dim>> detector_locations; 
  }; 
// @sect3{Equation data}  

// 像往常一样，我们必须定义我们的初始值、边界条件和右手边的函数。这次事情有点简单：我们考虑的是一个由初始条件驱动的问题，所以没有右手函数（尽管你可以在 step-23 中查找，看看如何做到这一点）。其次，没有边界条件：域的整个边界由吸收性边界条件组成。这就只剩下初始条件了，这里的事情也很简单，因为对于这个特殊的应用，只规定了压力的非零初始条件，而没有规定速度的非零初始条件（速度在初始时间为零）。

// 所以这就是我们所需要的：一个指定压力初始条件的类。在本程序所考虑的物理环境中，这些是小的吸收器，我们将其建模为一系列的小圆圈，我们假设压力盈余为1，而其他地方没有吸收，因此没有压力盈余。我们是这样做的（注意，如果我们想把这个程序扩展到不仅可以编译，而且可以运行，我们将不得不用三维源的位置来初始化源）。

  template <int dim> 
  class InitialValuesP : public Function<dim> 
  { 
  public: 
    virtual double value(const Point<dim> &p, 
                         const unsigned int /*component*/ = 0) const override 
    { 
      static const std::array<Source, 5> sources{ 
        {Source(Point<dim>(0, 0), 0.025), 
         Source(Point<dim>(-0.135, 0), 0.05), 
         Source(Point<dim>(0.17, 0), 0.03), 
         Source(Point<dim>(-0.25, 0), 0.02), 
         Source(Point<dim>(-0.05, -0.15), 0.015)}}; 

      for (const auto &source : sources) 
        if (p.distance(source.location) < source.radius) 
          return 1; 

      return 0; 
    } 

 
    struct Source 
    { 
      Source(const Point<dim> &l, const double r) 
        : location(l) 
        , radius(r) 
      {} 

      const Point<dim> location; 
      const double     radius; 
    }; 
  }; 
// @sect3{Implementation of the <code>TATForwardProblem</code> class}  

// 让我们再从构造函数开始。设置成员变量是很直接的。我们使用矿物油的声波速度（单位为毫米/微秒，是实验性生物医学成像中的常用单位），因为我们想和输出的许多实验都是在这里进行的。再次使用Crank-Nicolson方案，即theta被设定为0.5。随后选择时间步长以满足 $k = \frac hc$ ：这里我们把它初始化为一个无效的数字。

  template <int dim> 
  TATForwardProblem<dim>::TATForwardProblem() 
    : fe(1) 
    , dof_handler(triangulation) 
    , time_step(std::numeric_limits<double>::quiet_NaN()) 
    , time(time_step) 
    , timestep_number(1) 
    , theta(0.5) 
    , wave_speed(1.437) 
  { 

// 构造函数中的第二个任务是初始化存放检测器位置的数组。这个程序的结果与实验进行了比较，其中检测器间距的步长为2.25度，对应160个检测器位置。扫描圆的半径被选为中心和边界之间的一半，以避免不完善的边界条件带来的剩余反射破坏我们的数值结果。

// 然后按顺时针顺序计算探测器的位置。请注意，下面的内容当然只有在我们以2D计算时才有效，我们用一个断言来保护这个条件。如果我们以后想在三维中运行同样的程序，我们将不得不在这里添加代码来初始化三维中的探测器位置。由于断言的存在，我们不可能忘记这样做。

    Assert(dim == 2, ExcNotImplemented()); 

    const double detector_step_angle = 2.25; 
    const double detector_radius     = 0.5; 

    for (double detector_angle = 2 * numbers::PI; detector_angle >= 0; 
         detector_angle -= detector_step_angle / 360 * 2 * numbers::PI) 
      detector_locations.push_back( 
        Point<dim>(std::cos(detector_angle), std::sin(detector_angle)) * 
        detector_radius); 
  } 

//  @sect4{TATForwardProblem::setup_system}  

// 下面的系统几乎就是我们在  step-23  中已经做过的，但有两个重要的区别。首先，我们必须在原点周围创建一个半径为1的圆形（或球形）网格。这并不新鲜：我们之前在 step-6 和 step-10 中已经这样做了，在那里我们还解释了PolarManifold或SphericalManifold对象如何在细化单元时将新点放在同心圆上，我们在这里也将使用它。

// 我们必须确保的一点是，时间步长满足  step-23  的介绍中讨论的 CFL 条件。在那个程序中，我们通过设置一个与网格宽度相匹配的时间步长来确保这一点，但是这很容易出错，因为如果我们再细化一次网格，我们也必须确保时间步长有所改变。在这里，我们自动做到了这一点：我们向一个库函数询问任何单元的最小直径。然后我们设置 $k=\frac h{c_0}$  。唯一的问题是： $h$ 到底是什么？关键是，对于波浪方程来说，这个问题确实没有好的理论。众所周知，对于由矩形组成的均匀细化网格， $h$ 是最小边长。但对于一般四边形的网格，确切的关系似乎是未知的，也就是说，不知道单元格的什么属性与CFL条件有关。问题是，CFL条件来自于对拉普拉斯矩阵最小特征值的了解，而这只能对简单结构的网格进行分析计算。

// 这一切的结果是，我们并不十分确定我们应该对 $h$ 采取什么措施。函数 GridTools::minimal_cell_diameter 计算了所有单元的最小直径。如果单元格都是正方形或立方体，那么最小边长就是最小直径除以 <code>std::sqrt(dim)</code>  。我们简单地将此概括为非均匀网格的情况，没有理论上的理由。

// 唯一的其他重大变化是我们需要建立边界质量矩阵。我们将在下文中进一步评论这个问题。

  template <int dim> 
  void TATForwardProblem<dim>::setup_system() 
  { 
    const Point<dim> center; 
    GridGenerator::hyper_ball(triangulation, center, 1.); 
    triangulation.refine_global(7); 

    time_step = GridTools::minimal_cell_diameter(triangulation) / wave_speed / 
                std::sqrt(1. * dim); 

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
    mass_matrix.reinit(sparsity_pattern); 
    laplace_matrix.reinit(sparsity_pattern); 

    MatrixCreator::create_mass_matrix(dof_handler, 
                                      QGauss<dim>(fe.degree + 1), 
                                      mass_matrix); 
    MatrixCreator::create_laplace_matrix(dof_handler, 
                                         QGauss<dim>(fe.degree + 1), 
                                         laplace_matrix); 

// 如前所述，与 step-23 的第二个区别是，我们需要建立从吸收性边界条件中生长出来的边界质量矩阵。

// 第一个观察结果是，这个矩阵比常规质量矩阵要稀疏得多，因为没有一个具有纯内部支持的形状函数对这个矩阵有贡献。因此，我们可以根据这种情况优化存储模式，建立第二个稀疏模式，只包含我们需要的非零项。这里有一个权衡：首先，我们必须要有第二个稀疏模式对象，所以这需要花费内存。其次，与该稀疏性模式相连的矩阵将更小，因此需要更少的内存；用它进行矩阵-向量乘法也会更快。然而，最后一个论点是提示规模的论点：我们主要感兴趣的不是单独对边界矩阵进行矩阵-向量运算（尽管我们需要在每个时间步长对右侧向量进行一次运算），而是主要希望将其与两个方程中的第一个方程使用的其他矩阵相加，因为这是CG方法每个迭代都要与之相乘的，即明显更频繁。现在的情况是， SparseMatrix::add 类允许将一个矩阵添加到另一个矩阵中，但前提是它们使用相同的稀疏模式（原因是我们不能在稀疏模式创建后向矩阵添加非零条目，所以我们只是要求这两个矩阵具有相同的稀疏模式）。

// 所以，我们就用这个方法吧。

    boundary_matrix.reinit(sparsity_pattern); 

// 第二件要做的事是实际建立矩阵。在这里，我们需要对单元格的面进行积分，所以首先我们需要一个能在 <code>dim-1</code> 维对象上工作的正交对象。其次，FEValues的变体FEFaceValues，正如它的名字所暗示的，它可以在面上工作。最后，其他的变量是组装机器的一部分。所有这些我们都放在大括号里，以便将这些变量的范围限制在我们真正需要它们的地方。
//然后
//组装矩阵的实际行为是相当直接的：我们在所有单元中循环，在每个单元的所有面中循环，然后只在特定的面位于域的边界时做一些事情。像这样。

    { 
      const QGauss<dim - 1> quadrature_formula(fe.degree + 1); 
      FEFaceValues<dim>     fe_values(fe, 
                                  quadrature_formula, 
                                  update_values | update_JxW_values); 

      const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
      const unsigned int n_q_points    = quadrature_formula.size(); 

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 

      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary()) 
            { 
              cell_matrix = 0; 

              fe_values.reinit(cell, face); 

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    cell_matrix(i, j) += (fe_values.shape_value(i, q_point) * 
                                          fe_values.shape_value(j, q_point) * 
                                          fe_values.JxW(q_point)); 

              cell->get_dof_indices(local_dof_indices); 
              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                  boundary_matrix.add(local_dof_indices[i], 
                                      local_dof_indices[j], 
                                      cell_matrix(i, j)); 
            } 
    } 

    system_matrix.copy_from(mass_matrix); 
    system_matrix.add(time_step * time_step * theta * theta * wave_speed * 
                        wave_speed, 
                      laplace_matrix); 
    system_matrix.add(wave_speed * theta * time_step, boundary_matrix); 

    solution_p.reinit(dof_handler.n_dofs()); 
    old_solution_p.reinit(dof_handler.n_dofs()); 
    system_rhs_p.reinit(dof_handler.n_dofs()); 

    solution_v.reinit(dof_handler.n_dofs()); 
    old_solution_v.reinit(dof_handler.n_dofs()); 
    system_rhs_v.reinit(dof_handler.n_dofs()); 

    constraints.close(); 
  } 
// @sect4{TATForwardProblem::solve_p and TATForwardProblem::solve_v}  

// 下面两个函数，解决压力和速度变量的线性系统，几乎是逐字逐句地从 step-23 中提取的（除了主变量的名字从 $u$ 改为 $p$ ）。

  template <int dim> 
  void TATForwardProblem<dim>::solve_p() 
  { 
    SolverControl solver_control(1000, 1e-8 * system_rhs_p.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(system_matrix, solution_p, system_rhs_p, PreconditionIdentity()); 

    std::cout << "   p-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 

  template <int dim> 
  void TATForwardProblem<dim>::solve_v() 
  { 
    SolverControl solver_control(1000, 1e-8 * system_rhs_v.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    cg.solve(mass_matrix, solution_v, system_rhs_v, PreconditionIdentity()); 

    std::cout << "   v-equation: " << solver_control.last_step() 
              << " CG iterations." << std::endl; 
  } 

//  @sect4{TATForwardProblem::output_results}  

// 这里也是如此：该函数来自  step-23  。

  template <int dim> 
  void TATForwardProblem<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution_p, "P"); 
    data_out.add_data_vector(solution_v, "V"); 

    data_out.build_patches(); 

    const std::string filename = 
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
    DataOutBase::VtkFlags vtk_flags; 
    vtk_flags.compression_level = 
      DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
    std::ofstream output(filename); 
    data_out.write_vtu(output); 
  } 

//  @sect4{TATForwardProblem::run}  

// 这个做大部分工作的函数又和 step-23 中的差不多，尽管我们通过使用介绍中提到的向量G1和G2使事情变得更加清晰。与程序的整体内存消耗相比，引入几个临时向量并没有什么坏处。

// 这个函数唯一的变化是：首先，我们不必为速度 $v$ 预测初始值，因为我们知道它是零。其次，我们在构造函数中计算的检测器位置上评估解决方案。这是用 VectorTools::point_value 函数完成的。然后，这些值被写入我们在函数开始时打开的一个文件中。

  template <int dim> 
  void TATForwardProblem<dim>::run() 
  { 
    setup_system(); 

    VectorTools::project(dof_handler, 
                         constraints, 
                         QGauss<dim>(fe.degree + 1), 
                         InitialValuesP<dim>(), 
                         old_solution_p); 
    old_solution_v = 0; 

    std::ofstream detector_data("detectors.dat"); 

    Vector<double> tmp(solution_p.size()); 
    Vector<double> G1(solution_p.size()); 
    Vector<double> G2(solution_v.size()); 

    const double end_time = 0.7; 
    for (time = time_step; time <= end_time; 
         time += time_step, ++timestep_number) 
      { 
        std::cout << std::endl; 
        std::cout << "time_step " << timestep_number << " @ t=" << time 
                  << std::endl; 

        mass_matrix.vmult(G1, old_solution_p); 
        mass_matrix.vmult(tmp, old_solution_v); 
        G1.add(time_step * (1 - theta), tmp); 

        mass_matrix.vmult(G2, old_solution_v); 
        laplace_matrix.vmult(tmp, old_solution_p); 
        G2.add(-wave_speed * wave_speed * time_step * (1 - theta), tmp); 

        boundary_matrix.vmult(tmp, old_solution_p); 
        G2.add(wave_speed, tmp); 

        system_rhs_p = G1; 
        system_rhs_p.add(time_step * theta, G2); 

        solve_p(); 

        system_rhs_v = G2; 
        laplace_matrix.vmult(tmp, solution_p); 
        system_rhs_v.add(-time_step * theta * wave_speed * wave_speed, tmp); 

        boundary_matrix.vmult(tmp, solution_p); 
        system_rhs_v.add(-wave_speed, tmp); 

        solve_v(); 

        output_results(); 

        detector_data << time; 
        for (unsigned int i = 0; i < detector_locations.size(); ++i) 
          detector_data << " " 
                        << VectorTools::point_value(dof_handler, 
                                                    solution_p, 
                                                    detector_locations[i]) 
                        << " "; 
        detector_data << std::endl; 

        old_solution_p = solution_p; 
        old_solution_v = solution_v; 
      } 
  } 
} // namespace Step24 

//  @sect3{The <code>main</code> function}  

// 剩下的就是程序的主要功能了。这里没有什么是在前面几个程序中没有展示过的。

int main() 
{ 
  try 
    { 
      using namespace Step24; 

      TATForwardProblem<2> forward_problem_solver; 
      forward_problem_solver.run(); 
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


