CCTest_file/step-36.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 * Authors: Toby D. Young, Polish Academy of Sciences, 
 *          Wolfgang Bangerth, Texas A&M University 
 */ 


// @sect3{Include files}  

// 正如介绍中提到的，本程序基本上只是  step-4  的一个小修改版本。因此，以下大部分的include文件都是在那里使用的，或者至少是在以前的教程程序中已经使用的。

#include <deal.II/base/logstream.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/function_parser.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/full_matrix.h> 

// IndexSet用于设置每个  PETScWrappers::MPI::Vector:  的大小。
#include <deal.II/base/index_set.h> 

// PETSc出现在这里是因为SLEPc依赖于这个库。

#include <deal.II/lac/petsc_sparse_matrix.h> 
#include <deal.II/lac/petsc_vector.h> 

// 然后我们需要实际导入SLEPc提供的求解器接口。

#include <deal.II/lac/slepc_solver.h> 

// 我们还需要一些标准的C++。

#include <fstream> 
#include <iostream> 

// 最后，和以前的程序一样，我们将所有的deal.II类和函数名导入到本程序中所有的名字空间中。

namespace Step36 
{ 
  using namespace dealii; 
// @sect3{The <code>EigenvalueProblem</code> class template}  

// 下面是主类模板的类声明。它看起来和在  step-4  中已经展示过的差不多了。

  template <int dim> 
  class EigenvalueProblem 
  { 
  public: 
    EigenvalueProblem(const std::string &prm_file); 
    void run(); 

  private: 
    void         make_grid_and_dofs(); 
    void         assemble_system(); 
    unsigned int solve(); 
    void         output_results() const; 

    Triangulation<dim> triangulation; 
    FE_Q<dim>          fe; 
    DoFHandler<dim>    dof_handler; 

// 有了这些例外情况。对于我们的特征值问题，我们既需要左手边的刚度矩阵，也需要右手边的质量矩阵。我们还需要的不仅仅是一个解函数，而是一整套我们想要计算的特征函数，以及相应的特征值。

    PETScWrappers::SparseMatrix             stiffness_matrix, mass_matrix; 
    std::vector<PETScWrappers::MPI::Vector> eigenfunctions; 
    std::vector<double>                     eigenvalues; 

// 然后，我们需要一个对象来存储几个运行时参数，我们将在输入文件中指定。

    ParameterHandler parameters; 

// 最后，我们将有一个对象，包含对我们自由度的 "约束"。如果我们有自适应细化的网格（目前的程序中没有），这可能包括悬挂节点约束。这里，我们将存储边界节点的约束  $U_i=0$  。

    AffineConstraints<double> constraints; 
  }; 
// @sect3{Implementation of the <code>EigenvalueProblem</code> class}  
// @sect4{EigenvalueProblem::EigenvalueProblem}  

// 首先是构造函数。主要的新部分是处理运行时的输入参数。我们需要首先声明它们的存在，然后从输入文件中读取它们的值，该文件的名称被指定为该函数的参数。

  template <int dim> 
  EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file) 
    : fe(1) 
    , dof_handler(triangulation) 
  { 

// TODO研究为什么获得正确的特征值退化所需的最小细化步骤数为6

    parameters.declare_entry( 
      "Global mesh refinement steps", 
      "5", 
      Patterns::Integer(0, 20), 
      "The number of times the 1-cell coarse mesh should " 
      "be refined globally for our computations."); 
    parameters.declare_entry("Number of eigenvalues/eigenfunctions", 
                             "5", 
                             Patterns::Integer(0, 100), 
                             "The number of eigenvalues/eigenfunctions " 
                             "to be computed."); 
    parameters.declare_entry("Potential", 
                             "0", 
                             Patterns::Anything(), 
                             "A functional description of the potential."); 

    parameters.parse_input(prm_file); 
  } 
// @sect4{EigenvalueProblem::make_grid_and_dofs}  

// 下一个函数在域 $[-1,1]^d$ 上创建一个网格，根据输入文件的要求对其进行多次细化，然后给它附加一个DoFHandler，将矩阵和向量初始化为正确的大小。我们还建立了对应于边界值的约束  $u|_{\partial\Omega}=0$  。

// 对于矩阵，我们使用PETSc包装器。这些包装器能够在非零条目被添加时分配必要的内存。这看起来效率很低：我们可以先计算稀疏模式，用它来初始化矩阵，然后在我们插入条目时，我们可以确定我们不需要重新分配内存和释放之前使用的内存。一种方法是使用这样的代码。用
// @code
//    DynamicSparsityPattern
//       dsp (dof_handler.n_dofs(),
//            dof_handler.n_dofs());
//    DoFTools::make_sparsity_pattern (dof_handler, dsp);
//    dsp.compress ();
//    stiffness_matrix.reinit (dsp);
//    mass_matrix.reinit (dsp);
//  @endcode
//  代替下面两个 <code>reinit()</code> 的刚度和质量矩阵的调用。

// 不幸的是，这并不完全可行。上面的代码可能会导致在非零模式下的一些条目，我们只写零条目；最值得注意的是，对于那些属于边界节点的行和列的非对角线条目，这一点是成立的。这不应该是一个问题，但是不管什么原因，PETSc的ILU预处理程序（我们用来解决特征值求解器中的线性系统）不喜欢这些额外的条目，并以错误信息中止。

// 在没有任何明显的方法来避免这种情况的情况下，我们干脆选择第二种最好的方法，即让PETSc在必要时分配内存。也就是说，由于这不是一个时间上的关键部分，这整个事件就不再重要了。

  template <int dim> 
  void EigenvalueProblem<dim>::make_grid_and_dofs() 
  { 
    GridGenerator::hyper_cube(triangulation, -1, 1); 
    triangulation.refine_global( 
      parameters.get_integer("Global mesh refinement steps")); 
    dof_handler.distribute_dofs(fe); 

    DoFTools::make_zero_boundary_constraints(dof_handler, constraints); 
    constraints.close(); 

    stiffness_matrix.reinit(dof_handler.n_dofs(), 
                            dof_handler.n_dofs(), 
                            dof_handler.max_couplings_between_dofs()); 
    mass_matrix.reinit(dof_handler.n_dofs(), 
                       dof_handler.n_dofs(), 
                       dof_handler.max_couplings_between_dofs()); 

// 下一步是处理特征谱的问题。在这种情况下，输出是特征值和特征函数，所以我们将特征函数和特征值列表的大小设置为与我们在输入文件中要求的一样大。当使用 PETScWrappers::MPI::Vector, 时，Vector是使用IndexSet初始化的。IndexSet不仅用于调整 PETScWrappers::MPI::Vector 的大小，而且还将 PETScWrappers::MPI::Vector 中的一个索引与一个自由度联系起来（更详细的解释见 step-40 ）。函数complete_index_set()创建了一个IndexSet，每个有效的索引都是这个集合的一部分。请注意，这个程序只能按顺序运行，如果并行使用，将抛出一个异常。

    IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs(); 
    eigenfunctions.resize( 
      parameters.get_integer("Number of eigenvalues/eigenfunctions")); 
    for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
      eigenfunctions[i].reinit(eigenfunction_index_set, MPI_COMM_WORLD); 

    eigenvalues.resize(eigenfunctions.size()); 
  } 
// @sect4{EigenvalueProblem::assemble_system}  

// 在这里，我们从局部贡献 $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf x) \cdot \nabla\varphi_j(\mathbf x) + V(\mathbf x)\varphi_i(\mathbf x)\varphi_j(\mathbf x)$ 和 $M^K_{ij} = \int_K \varphi_i(\mathbf x)\varphi_j(\mathbf x)$ 中分别组合出全局刚度和质量矩阵。如果你看过以前的教程程序，这个函数应该会很熟悉。唯一新的东西是使用我们从输入文件中得到的表达式，设置一个描述势 $V(\mathbf x)$ 的对象。然后我们需要在每个单元的正交点上评估这个对象。如果你见过如何评估函数对象（例如，见 step-5 中的系数），这里的代码也会显得相当熟悉。

  template <int dim> 
  void EigenvalueProblem<dim>::assemble_system() 
  { 
    QGauss<dim> quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    FunctionParser<dim> potential; 
    potential.initialize(FunctionParser<dim>::default_variable_names(), 
                         parameters.get("Potential"), 
                         typename FunctionParser<dim>::ConstMap()); 

    std::vector<double> potential_values(n_q_points); 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        cell_stiffness_matrix = 0; 
        cell_mass_matrix      = 0; 

        potential.value_list(fe_values.get_quadrature_points(), 
                             potential_values); 

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              { 
                cell_stiffness_matrix(i, j) +=           // 
                  (fe_values.shape_grad(i, q_point) *    // 
                     fe_values.shape_grad(j, q_point)    // 
                   +                                     // 
                   potential_values[q_point] *           // 
                     fe_values.shape_value(i, q_point) * // 
                     fe_values.shape_value(j, q_point)   // 
                   ) *                                   // 
                  fe_values.JxW(q_point);                // 

                cell_mass_matrix(i, j) +=              // 
                  (fe_values.shape_value(i, q_point) * // 
                   fe_values.shape_value(j, q_point)   // 
                   ) *                                 // 
                  fe_values.JxW(q_point);              // 
              } 

// 现在我们有了本地矩阵的贡献，我们把它们转移到全局对象中，并处理好零边界约束。

        cell->get_dof_indices(local_dof_indices); 

        constraints.distribute_local_to_global(cell_stiffness_matrix, 
                                               local_dof_indices, 
                                               stiffness_matrix); 
        constraints.distribute_local_to_global(cell_mass_matrix, 
                                               local_dof_indices, 
                                               mass_matrix); 
      } 

// 在函数的最后，我们告诉PETSc，矩阵现在已经完全组装好了，稀疏矩阵表示法现在可以被压缩了，因为不会再添加任何条目。

    stiffness_matrix.compress(VectorOperation::add); 
    mass_matrix.compress(VectorOperation::add); 

// 在离开函数之前，我们计算虚假的特征值，这些特征值是由零Dirichlet约束引入到系统中的。正如介绍中所讨论的，使用Dirichlet边界条件，加上位于域的边界的自由度仍然是我们所求解的线性系统的一部分，引入了一些虚假的特征值。下面，我们输出它们所处的区间，以确保我们在计算中出现时可以忽略它们。

    double min_spurious_eigenvalue = std::numeric_limits<double>::max(), 
           max_spurious_eigenvalue = -std::numeric_limits<double>::max(); 

    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
      if (constraints.is_constrained(i)) 
        { 
          const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i); 
          min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev); 
          max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev); 
        } 

    std::cout << "   Spurious eigenvalues are all in the interval " 
              << "[" << min_spurious_eigenvalue << "," 
              << max_spurious_eigenvalue << "]" << std::endl; 
  } 
// @sect4{EigenvalueProblem::solve}  

// 这是该程序的关键新功能。现在系统已经设置好了，现在是实际解决问题的好时机：和其他例子一样，这是使用 "解决 "程序来完成的。从本质上讲，它的工作原理与其他程序一样：你设置一个SolverControl对象，描述我们要解决的线性系统的精度，然后我们选择我们想要的解算器类型。这里我们选择了SLEPc的Krylov-Schur求解器，对于这类问题来说，这是一个相当快速和强大的选择。

  template <int dim> 
  unsigned int EigenvalueProblem<dim>::solve() 
  { 

// 我们从这里开始，就像我们通常做的那样，指定我们想要的收敛控制。

    SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9); 
    SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control); 

// 在我们实际求解特征函数和-值之前，我们还必须选择哪一组特征值来求解。让我们选择那些实部最小的特征值和相应的特征函数（事实上，我们在这里解决的问题是对称的，所以特征值是纯实部的）。之后，我们就可以真正让SLEPc做它的工作了。

    eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL); 

    eigensolver.set_problem_type(EPS_GHEP); 

    eigensolver.solve(stiffness_matrix, 
                      mass_matrix, 
                      eigenvalues, 
                      eigenfunctions, 
                      eigenfunctions.size()); 

// 上述调用的输出是一组向量和数值。在特征值问题中，特征函数只确定到一个常数，这个常数可以很随意地固定。由于对特征值问题的原点一无所知，SLEPc除了将特征向量归一到 $l_2$ （向量）准则外，没有其他选择。不幸的是，这个规范与我们从特征函数角度可能感兴趣的任何规范没有什么关系： $L_2(\Omega)$ 规范，或者也许是 $L_\infty(\Omega)$ 规范。

//让我们选择后者，重新划分特征函数的尺度，使其具有 $\|\phi_i(\mathbf x)\|_{L^\infty(\Omega)}=1$ 而不是 $\|\Phi\|_{l_2}=1$ （其中 $\phi_i$ 是 $i$ 第三个特征<i>function</i>， $\Phi_i$ 是相应的结点值矢量）。对于这里选择的 $Q_1$ 元素，我们知道函数 $\phi_i(\mathbf x)$ 的最大值是在其中一个节点达到的，所以 $\max_{\mathbf x}\phi_i(\mathbf x)=\max_j (\Phi_i)_j$ ，使得在 $L_\infty$ 准则下的归一化是微不足道的。请注意，如果我们选择 $Q_k$ 元素与 $k>1$ ，这就不容易了：在那里，一个函数的最大值不一定要在一个节点上达到，所以 $\max_{\mathbf x}\phi_i(\mathbf x)\ge\max_j (\Phi_i)_j$ （尽管平等通常几乎是真的）。

    for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
      eigenfunctions[i] /= eigenfunctions[i].linfty_norm(); 

// 最后返回收敛所需的迭代次数。

    return solver_control.last_step(); 
  } 
// @sect4{EigenvalueProblem::output_results}  

// 这是本程序的最后一个重要功能。它使用DataOut类来生成特征函数的图形输出，以便以后进行可视化。它的工作原理与其他许多教程中的程序一样。

// 整个函数的集合被输出为一个单一的VTK文件。

  template <int dim> 
  void EigenvalueProblem<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 

    for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
      data_out.add_data_vector(eigenfunctions[i], 
                               std::string("eigenfunction_") + 
                                 Utilities::int_to_string(i)); 

// 唯一值得讨论的可能是，由于势在输入文件中被指定为函数表达式，因此最好能将其与特征函数一起以图形形式表示。实现这一目的的过程相对简单：我们建立一个代表 $V(\mathbf x)$ 的对象，然后将这个连续函数插值到有限元空间。我们还将结果附加到DataOut对象上，以便进行可视化。

    Vector<double> projected_potential(dof_handler.n_dofs()); 
    { 
      FunctionParser<dim> potential; 
      potential.initialize(FunctionParser<dim>::default_variable_names(), 
                           parameters.get("Potential"), 
                           typename FunctionParser<dim>::ConstMap()); 
      VectorTools::interpolate(dof_handler, potential, projected_potential); 
    } 
    data_out.add_data_vector(projected_potential, "interpolated_potential"); 

    data_out.build_patches(); 

    std::ofstream output("eigenvectors.vtk"); 
    data_out.write_vtk(output); 
  } 
// @sect4{EigenvalueProblem::run}  

// 这是一个对一切都有顶层控制的函数。它几乎与  step-4  中的内容完全相同。

  template <int dim> 
  void EigenvalueProblem<dim>::run() 
  { 
    make_grid_and_dofs(); 

    std::cout << "   Number of active cells:       " 
              << triangulation.n_active_cells() << std::endl 
              << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

    assemble_system(); 

    const unsigned int n_iterations = solve(); 
    std::cout << "   Solver converged in " << n_iterations << " iterations." 
              << std::endl; 

    output_results(); 

    std::cout << std::endl; 
    for (unsigned int i = 0; i < eigenvalues.size(); ++i) 
      std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i] 
                << std::endl; 
  } 
} // namespace Step36 
// @sect3{The <code>main</code> function}  
int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step36; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

// 这个程序只能在串行中运行。否则，将抛出一个异常。

      AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
                  ExcMessage( 
                    "This program can only be run in serial, use ./step-36")); 

      EigenvalueProblem<2> problem("step-36.prm"); 
      problem.run(); 
    } 

// 在这期间，我们一直在注意是否有任何异常应该被生成。如果是这样的话，我们就会惊慌失措...

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

// 如果没有抛出异常，我们就告诉程序不要再胡闹了，乖乖地退出。

  std::cout << std::endl << "   Job done." << std::endl; 

  return 0; 
} 


