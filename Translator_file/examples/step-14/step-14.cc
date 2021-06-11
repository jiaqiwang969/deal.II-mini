CCTest_file/step-14.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2002 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, ETH Zurich, 2002 
 */ 



// 从众所周知的事情开始......

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/thread_management.h> 
#include <deal.II/base/work_stream.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_tools.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <algorithm> 
#include <fstream> 
#include <iostream> 
#include <list> 
#include <memory> 
#include <numeric> 

// 最后一步和以前所有的程序一样。

namespace Step14 
{ 
  using namespace dealii; 
// @sect3{Evaluating the solution}  

// 正如介绍中提到的，该程序的重要部分只是从 step-13 的例子程序中拿过来的。因此，我们只对那些新的东西进行评论。

// 首先，评估解决方案的框架没有改变，即基类是相同的，评估网格点上的解决方案的类也没有改变。

  namespace Evaluation 
  { 
// @sect4{The EvaluationBase class}  
    template <int dim> 
    class EvaluationBase 
    { 
    public: 
      virtual ~EvaluationBase() = default; 

      void set_refinement_cycle(const unsigned int refinement_cycle); 

      virtual void operator()(const DoFHandler<dim> &dof_handler, 
                              const Vector<double> & solution) const = 0; 

    protected: 
      unsigned int refinement_cycle; 
    }; 

    template <int dim> 
    void EvaluationBase<dim>::set_refinement_cycle(const unsigned int step) 
    { 
      refinement_cycle = step; 
    } 
// @sect4{The PointValueEvaluation class}  
    template <int dim> 
    class PointValueEvaluation : public EvaluationBase<dim> 
    { 
    public: 
      PointValueEvaluation(const Point<dim> &evaluation_point); 

      virtual void operator()(const DoFHandler<dim> &dof_handler, 
                              const Vector<double> & solution) const override; 

      DeclException1( 
        ExcEvaluationPointNotFound, 
        Point<dim>, 
        << "The evaluation point " << arg1 
        << " was not found among the vertices of the present grid."); 

    private: 
      const Point<dim> evaluation_point; 
    }; 

    template <int dim> 
    PointValueEvaluation<dim>::PointValueEvaluation( 
      const Point<dim> &evaluation_point) 
      : evaluation_point(evaluation_point) 
    {} 

    template <int dim> 
    void PointValueEvaluation<dim>:: 
         operator()(const DoFHandler<dim> &dof_handler, 
               const Vector<double> & solution) const 
    { 
      double point_value = 1e20; 

      bool evaluation_point_found = false; 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (!evaluation_point_found) 
          for (const auto vertex : cell->vertex_indices()) 
            if (cell->vertex(vertex).distance(evaluation_point) < 
                cell->diameter() * 1e-8) 
              { 
                point_value = solution(cell->vertex_dof_index(vertex, 0)); 

                evaluation_point_found = true; 
                break; 
              } 

      AssertThrow(evaluation_point_found, 
                  ExcEvaluationPointNotFound(evaluation_point)); 

      std::cout << "   Point value=" << point_value << std::endl; 
    } 
// @sect4{The PointXDerivativeEvaluation class}  

// 除了实现在一个点上求解的类，我们在这里提供一个在网格点上求梯度的类。由于一般情况下，有限元函数的梯度在一个顶点上是不连续的，所以我们在这里要稍微小心一点。我们要做的是在所有单元中循环，即使我们已经在一个单元中找到了点，也要使用所有相邻单元中的顶点梯度的平均值。

// 鉴于 <code>PointValueEvaluation</code> 类的接口，这个类的声明没有提供什么惊喜，构造函数也没有。

    template <int dim> 
    class PointXDerivativeEvaluation : public EvaluationBase<dim> 
    { 
    public: 
      PointXDerivativeEvaluation(const Point<dim> &evaluation_point); 

      virtual void operator()(const DoFHandler<dim> &dof_handler, 
                              const Vector<double> & solution) const; 

      DeclException1( 
        ExcEvaluationPointNotFound, 
        Point<dim>, 
        << "The evaluation point " << arg1 
        << " was not found among the vertices of the present grid."); 

    private: 
      const Point<dim> evaluation_point; 
    }; 

    template <int dim> 
    PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation( 
      const Point<dim> &evaluation_point) 
      : evaluation_point(evaluation_point) 
    {} 

// 更有趣的事情发生在进行实际评估的函数中。

    template <int dim> 
    void PointXDerivativeEvaluation<dim>:: 
         operator()(const DoFHandler<dim> &dof_handler, 
               const Vector<double> & solution) const 
    { 

// 这次用一些有用的东西来初始化返回值，因为我们要把一些贡献加起来，之后再取平均值。

      double point_derivative = 0; 

// ...然后有一些对象，其含义将在下面变得清晰...

      QTrapezoid<dim>             vertex_quadrature; 
      FEValues<dim>               fe_values(dof_handler.get_fe(), 
                              vertex_quadrature, 
                              update_gradients | update_quadrature_points); 
      std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size()); 

// ...接下来循环所有单元格及其顶点，并计算顶点被发现的频率。

      unsigned int evaluation_point_hits = 0; 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto vertex : cell->vertex_indices()) 
          if (cell->vertex(vertex) == evaluation_point) 
            { 

// 现在事情不再那么简单了，因为我们不能像以前那样得到有限元场的梯度，我们只需要在一个顶点上选择一个自由度。        相反，我们必须在这个单元上，在某一点上评估有限元场。如你所知，在某一点上评估有限元场是通过 <code>FEValues</code> 类完成的，所以我们使用它。问题是： <code>FEValues</code> 对象需要给定一个正交公式，然后可以计算正交点的有限元量值。在这里，我们并不想做正交，我们只是想指定一些点!         尽管如此，还是选择同样的方式：使用一个特殊的正交规则，点在顶点，因为这些是我们感兴趣的。适当的规则是梯形规则，所以这就是我们上面使用该规则的原因。        因此：在这个单元格上初始化 <code>FEValues</code> 对象。

              fe_values.reinit(cell); 

// 并在顶点提取解向量的梯度。

              fe_values.get_function_gradients(solution, solution_gradients); 

// 现在我们有了所有顶点的梯度，所以选出属于评估点的那一个（注意顶点的顺序不一定和正交点的顺序一样）。

              unsigned int q_point = 0; 
              for (; q_point < solution_gradients.size(); ++q_point) 
                if (fe_values.quadrature_point(q_point) == evaluation_point) 
                  break; 

// 检查是否确实找到了评估点。

              Assert(q_point < solution_gradients.size(), ExcInternalError()); 

// 如果是这样，就把那里的梯度的X导数作为我们感兴趣的值，并增加计数器，表示我们向该变量添加了多少次。

              point_derivative += solution_gradients[q_point][0]; 
              ++evaluation_point_hits; 

// 最后跳出最内层的循环，遍历当前单元格的顶点，因为如果我们在一个顶点找到了评估点，就不可能在后面的顶点也找到。

              break; 
            } 

// 现在我们已经循环了所有的单元和顶点，所以检查是否找到了这个点。

      AssertThrow(evaluation_point_hits > 0, 
                  ExcEvaluationPointNotFound(evaluation_point)); 

// 我们已经简单地将所有相邻的单元格的贡献相加，所以我们仍然要计算出平均值。一旦完成，报告状态。

      point_derivative /= evaluation_point_hits; 
      std::cout << "   Point x-derivative=" << point_derivative << std::endl; 
    } 

//  @sect4{The GridOutput class}  

// 由于这个程序有一个更困难的结构（它除了计算一个原始解之外，还计算了一个对偶解），所以写出解不再由一个评估对象来完成，因为我们想把两个解同时写进一个文件，这需要一些比评估类可用的信息。

// 然而，我们也想看看生成的网格。这也可以通过一个这样的类来完成。它的结构类似于前面例子程序中的 <code>SolutionOutput</code> 类，所以我们在这里不做更详细的讨论。此外，这里所使用的一切都已经在前面的例子程序中使用过了。

    template <int dim> 
    class GridOutput : public EvaluationBase<dim> 
    { 
    public: 
      GridOutput(const std::string &output_name_base); 

      virtual void operator()(const DoFHandler<dim> &dof_handler, 
                              const Vector<double> & solution) const override; 

    private: 
      const std::string output_name_base; 
    }; 

    template <int dim> 
    GridOutput<dim>::GridOutput(const std::string &output_name_base) 
      : output_name_base(output_name_base) 
    {} 

    template <int dim> 
    void GridOutput<dim>::operator()(const DoFHandler<dim> &dof_handler, 
                                     const Vector<double> & /*solution*/) const 
    { 
      std::ofstream out(output_name_base + "-" + 
                        std::to_string(this->refinement_cycle) + ".svg"); 
      GridOut().write_svg(dof_handler.get_triangulation(), out); 
    } 
  } // namespace Evaluation 
// @sect3{The Laplace solver classes}  

// 接下来是实际的求解器类。同样，我们只讨论与之前程序的不同之处。

  namespace LaplaceSolver 
  { 
// @sect4{The Laplace solver base class}  

// 这个类几乎没有变化，只是多声明了两个函数。  <code>output_solution</code> 将用于从派生类计算的实际解决方案中生成输出文件，以及 <code>set_refinement_cycle</code> 函数，测试框架通过该函数将细化周期的编号设置为该类中的一个局部变量；该编号随后将用于生成解决方案输出的文件名。

    template <int dim> 
    class Base 
    { 
    public: 
      Base(Triangulation<dim> &coarse_grid); 
      virtual ~Base() = default; 

      virtual void solve_problem() = 0; 
      virtual void postprocess( 
        const Evaluation::EvaluationBase<dim> &postprocessor) const = 0; 
      virtual void         refine_grid()                            = 0; 
      virtual unsigned int n_dofs() const                           = 0; 

      virtual void set_refinement_cycle(const unsigned int cycle); 

      virtual void output_solution() const = 0; 

    protected: 
      const SmartPointer<Triangulation<dim>> triangulation; 

      unsigned int refinement_cycle; 
    }; 

    template <int dim> 
    Base<dim>::Base(Triangulation<dim> &coarse_grid) 
      : triangulation(&coarse_grid) 
      , refinement_cycle(numbers::invalid_unsigned_int) 
    {} 

    template <int dim> 
    void Base<dim>::set_refinement_cycle(const unsigned int cycle) 
    { 
      refinement_cycle = cycle; 
    } 
// @sect4{The Laplace Solver class}  

// 同样地， <code>Solver</code> 类完全没有变化，因此将不进行讨论。

    template <int dim> 
    class Solver : public virtual Base<dim> 
    { 
    public: 
      Solver(Triangulation<dim> &       triangulation, 
             const FiniteElement<dim> & fe, 
             const Quadrature<dim> &    quadrature, 
             const Quadrature<dim - 1> &face_quadrature, 
             const Function<dim> &      boundary_values); 
      virtual ~Solver() override; 

      virtual void solve_problem() override; 

      virtual void postprocess( 
        const Evaluation::EvaluationBase<dim> &postprocessor) const override; 

      virtual unsigned int n_dofs() const override; 

    protected: 
      const SmartPointer<const FiniteElement<dim>>  fe; 
      const SmartPointer<const Quadrature<dim>>     quadrature; 
      const SmartPointer<const Quadrature<dim - 1>> face_quadrature; 
      DoFHandler<dim>                               dof_handler; 
      Vector<double>                                solution; 
      const SmartPointer<const Function<dim>>       boundary_values; 

      virtual void assemble_rhs(Vector<double> &rhs) const = 0; 

    private: 
      struct LinearSystem 
      { 
        LinearSystem(const DoFHandler<dim> &dof_handler); 

        void solve(Vector<double> &solution) const; 

        AffineConstraints<double> hanging_node_constraints; 
        SparsityPattern           sparsity_pattern; 
        SparseMatrix<double>      matrix; 
        Vector<double>            rhs; 
      }; 

// 该类的其余部分基本上也是 step-13 的副本，包括使用WorkStream框架并行计算线性系统所需的数据结构和函数。

      struct AssemblyScratchData 
      { 
        AssemblyScratchData(const FiniteElement<dim> &fe, 
                            const Quadrature<dim> &   quadrature); 
        AssemblyScratchData(const AssemblyScratchData &scratch_data); 

        FEValues<dim> fe_values; 
      }; 

      struct AssemblyCopyData 
      { 
        FullMatrix<double>                   cell_matrix; 
        std::vector<types::global_dof_index> local_dof_indices; 
      }; 

      void assemble_linear_system(LinearSystem &linear_system); 

      void local_assemble_matrix( 
        const typename DoFHandler<dim>::active_cell_iterator &cell, 
        AssemblyScratchData &                                 scratch_data, 
        AssemblyCopyData &                                    copy_data) const; 

      void copy_local_to_global(const AssemblyCopyData &copy_data, 
                                LinearSystem &          linear_system) const; 
    }; 

    template <int dim> 
    Solver<dim>::Solver(Triangulation<dim> &       triangulation, 
                        const FiniteElement<dim> & fe, 
                        const Quadrature<dim> &    quadrature, 
                        const Quadrature<dim - 1> &face_quadrature, 
                        const Function<dim> &      boundary_values) 
      : Base<dim>(triangulation) 
      , fe(&fe) 
      , quadrature(&quadrature) 
      , face_quadrature(&face_quadrature) 
      , dof_handler(triangulation) 
      , boundary_values(&boundary_values) 
    {} 

    template <int dim> 
    Solver<dim>::~Solver() 
    { 
      dof_handler.clear(); 
    } 

    template <int dim> 
    void Solver<dim>::solve_problem() 
    { 
      dof_handler.distribute_dofs(*fe); 
      solution.reinit(dof_handler.n_dofs()); 

      LinearSystem linear_system(dof_handler); 
      assemble_linear_system(linear_system); 
      linear_system.solve(solution); 
    } 

    template <int dim> 
    void Solver<dim>::postprocess( 
      const Evaluation::EvaluationBase<dim> &postprocessor) const 
    { 
      postprocessor(dof_handler, solution); 
    } 

    template <int dim> 
    unsigned int Solver<dim>::n_dofs() const 
    { 
      return dof_handler.n_dofs(); 
    } 

// 以下几个函数和构造函数是逐字复制的，来自  step-13  。

    template <int dim> 
    void Solver<dim>::assemble_linear_system(LinearSystem &linear_system) 
    { 
      Threads::Task<void> rhs_task = 
        Threads::new_task(&Solver<dim>::assemble_rhs, *this, linear_system.rhs); 

      auto worker = 
        [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
               AssemblyScratchData &scratch_data, 
               AssemblyCopyData &   copy_data) { 
          this->local_assemble_matrix(cell, scratch_data, copy_data); 
        }; 

      auto copier = [this, &linear_system](const AssemblyCopyData &copy_data) { 
        this->copy_local_to_global(copy_data, linear_system); 
      }; 

      WorkStream::run(dof_handler.begin_active(), 
                      dof_handler.end(), 
                      worker, 
                      copier, 
                      AssemblyScratchData(*fe, *quadrature), 
                      AssemblyCopyData()); 
      linear_system.hanging_node_constraints.condense(linear_system.matrix); 

      std::map<types::global_dof_index, double> boundary_value_map; 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               0, 
                                               *boundary_values, 
                                               boundary_value_map); 

      rhs_task.join(); 
      linear_system.hanging_node_constraints.condense(linear_system.rhs); 

      MatrixTools::apply_boundary_values(boundary_value_map, 
                                         linear_system.matrix, 
                                         solution, 
                                         linear_system.rhs); 
    } 

    template <int dim> 
    Solver<dim>::AssemblyScratchData::AssemblyScratchData( 
      const FiniteElement<dim> &fe, 
      const Quadrature<dim> &   quadrature) 
      : fe_values(fe, quadrature, update_gradients | update_JxW_values) 
    {} 

    template <int dim> 
    Solver<dim>::AssemblyScratchData::AssemblyScratchData( 
      const AssemblyScratchData &scratch_data) 
      : fe_values(scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  update_gradients | update_JxW_values) 
    {} 

    template <int dim> 
    void Solver<dim>::local_assemble_matrix( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      AssemblyScratchData &                                 scratch_data, 
      AssemblyCopyData &                                    copy_data) const 
    { 
      const unsigned int dofs_per_cell = fe->n_dofs_per_cell(); 
      const unsigned int n_q_points    = quadrature->size(); 

      copy_data.cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 

      copy_data.local_dof_indices.resize(dofs_per_cell); 

      scratch_data.fe_values.reinit(cell); 

      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < dofs_per_cell; ++j) 
            copy_data.cell_matrix(i, j) += 
              (scratch_data.fe_values.shape_grad(i, q_point) * 
               scratch_data.fe_values.shape_grad(j, q_point) * 
               scratch_data.fe_values.JxW(q_point)); 

      cell->get_dof_indices(copy_data.local_dof_indices); 
    } 

    template <int dim> 
    void Solver<dim>::copy_local_to_global(const AssemblyCopyData &copy_data, 
                                           LinearSystem &linear_system) const 
    { 
      for (unsigned int i = 0; i < copy_data.local_dof_indices.size(); ++i) 
        for (unsigned int j = 0; j < copy_data.local_dof_indices.size(); ++j) 
          linear_system.matrix.add(copy_data.local_dof_indices[i], 
                                   copy_data.local_dof_indices[j], 
                                   copy_data.cell_matrix(i, j)); 
    } 

// 现在是实现线性系统类中动作的函数。首先，构造函数将所有数据元素初始化为正确的大小，并设置了一些额外的数据结构，例如由于悬挂节点而产生的约束。由于设置悬空节点和找出矩阵的非零元素是独立的，所以我们以并行方式进行（如果库被配置为使用并发，至少是这样；否则，这些动作是按顺序执行的）。注意，我们只启动一个线程，并在主线程中做第二个动作。由于只生成一个线程，我们在这里不使用 <code>Threads::TaskGroup</code> 类，而是直接使用创建的一个任务对象来等待这个特定任务的退出。这个方法与我们在上面 <code>Solver::assemble_linear_system()</code> 中使用的方法大致相同。

// 注意，获取 <code>DoFTools::make_hanging_node_constraints</code> 函数的地址有点麻烦，因为实际上有三个这个名字的函数，每个支持的空间维度都有一个。在C++中，获取重载函数的地址有些复杂，因为在这种情况下，操作符 <code>&</code> 会返回一组值（所有具有该名称的函数的地址），然后选择正确的函数是下一步的事情。如果上下文决定采取哪一个（例如通过分配给一个已知类型的函数指针），那么编译器可以自己做，但是如果这组指针应作为一个采取模板的函数的参数，编译器可以选择所有的，而不偏向于一个。因此，我们必须向编译器说明我们想要哪一个；为此，我们可以使用cast，但为了更清楚，我们把它分配给一个具有正确类型的临时 <code>mhnc_p</code> （简称<code>pointer to make_hanging_node_constraints</code>），并使用这个指针代替。

    template <int dim> 
    Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler) 
    { 
      hanging_node_constraints.clear(); 

      void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) = 
        &DoFTools::make_hanging_node_constraints; 

// 启动一个辅助任务，然后在主线程上继续进行

      Threads::Task<void> side_task = 
        Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints); 

      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp); 

// 等到边上的任务完成后再继续前进

      side_task.join(); 

      hanging_node_constraints.close(); 
      hanging_node_constraints.condense(dsp); 
      sparsity_pattern.copy_from(dsp); 

      matrix.reinit(sparsity_pattern); 
      rhs.reinit(dof_handler.n_dofs()); 
    } 

    template <int dim> 
    void Solver<dim>::LinearSystem::solve(Vector<double> &solution) const 
    { 
      SolverControl            solver_control(5000, 1e-12); 
      SolverCG<Vector<double>> cg(solver_control); 

      PreconditionSSOR<SparseMatrix<double>> preconditioner; 
      preconditioner.initialize(matrix, 1.2); 

      cg.solve(matrix, solution, rhs, preconditioner); 

      hanging_node_constraints.distribute(solution); 
    } 

//  @sect4{The PrimalSolver class}  

//  <code>PrimalSolver</code> 类除了实现 <code>output_solution</code> 函数外，也基本没有变化。我们在这个程序中保留了 <code>GlobalRefinement</code> and <code>RefinementKelly</code> 类，它们就可以依赖这个函数的默认实现，这个函数只是输出原始解。实现双重加权误差估计的类将自行重载这个函数，以同时输出双重解。

    template <int dim> 
    class PrimalSolver : public Solver<dim> 
    { 
    public: 
      PrimalSolver(Triangulation<dim> &       triangulation, 
                   const FiniteElement<dim> & fe, 
                   const Quadrature<dim> &    quadrature, 
                   const Quadrature<dim - 1> &face_quadrature, 
                   const Function<dim> &      rhs_function, 
                   const Function<dim> &      boundary_values); 

      virtual void output_solution() const override; 

    protected: 
      const SmartPointer<const Function<dim>> rhs_function; 
      virtual void assemble_rhs(Vector<double> &rhs) const override; 
    }; 

    template <int dim> 
    PrimalSolver<dim>::PrimalSolver(Triangulation<dim> &       triangulation, 
                                    const FiniteElement<dim> & fe, 
                                    const Quadrature<dim> &    quadrature, 
                                    const Quadrature<dim - 1> &face_quadrature, 
                                    const Function<dim> &      rhs_function, 
                                    const Function<dim> &      boundary_values) 
      : Base<dim>(triangulation) 
      , Solver<dim>(triangulation, 
                    fe, 
                    quadrature, 
                    face_quadrature, 
                    boundary_values) 
      , rhs_function(&rhs_function) 
    {} 

    template <int dim> 
    void PrimalSolver<dim>::output_solution() const 
    { 
      DataOut<dim> data_out; 
      data_out.attach_dof_handler(this->dof_handler); 
      data_out.add_data_vector(this->solution, "solution"); 
      data_out.build_patches(); 

      std::ofstream out("solution-" + std::to_string(this->refinement_cycle) + 
                        ".vtu"); 
      data_out.write(out, DataOutBase::vtu); 
    } 

    template <int dim> 
    void PrimalSolver<dim>::assemble_rhs(Vector<double> &rhs) const 
    { 
      FEValues<dim> fe_values(*this->fe, 
                              *this->quadrature, 
                              update_values | update_quadrature_points | 
                                update_JxW_values); 

      const unsigned int dofs_per_cell = this->fe->n_dofs_per_cell(); 
      const unsigned int n_q_points    = this->quadrature->size(); 

      Vector<double>                       cell_rhs(dofs_per_cell); 
      std::vector<double>                  rhs_values(n_q_points); 
      std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

      for (const auto &cell : this->dof_handler.active_cell_iterators()) 
        { 
          cell_rhs = 0; 

          fe_values.reinit(cell); 

          rhs_function->value_list(fe_values.get_quadrature_points(), 
                                   rhs_values); 

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q) 
                              rhs_values[q_point] *               // f((x_q) 
                              fe_values.JxW(q_point));            // dx 

          cell->get_dof_indices(local_dof_indices); 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            rhs(local_dof_indices[i]) += cell_rhs(i); 
        } 
    } 
// @sect4{The RefinementGlobal and RefinementKelly classes}  

// 对于下面的两个类，与上面的大多数情况相同：类是按原样取自前面的例子。

    template <int dim> 
    class RefinementGlobal : public PrimalSolver<dim> 
    { 
    public: 
      RefinementGlobal(Triangulation<dim> &       coarse_grid, 
                       const FiniteElement<dim> & fe, 
                       const Quadrature<dim> &    quadrature, 
                       const Quadrature<dim - 1> &face_quadrature, 
                       const Function<dim> &      rhs_function, 
                       const Function<dim> &      boundary_values); 

      virtual void refine_grid() override; 
    }; 

    template <int dim> 
    RefinementGlobal<dim>::RefinementGlobal( 
      Triangulation<dim> &       coarse_grid, 
      const FiniteElement<dim> & fe, 
      const Quadrature<dim> &    quadrature, 
      const Quadrature<dim - 1> &face_quadrature, 
      const Function<dim> &      rhs_function, 
      const Function<dim> &      boundary_values) 
      : Base<dim>(coarse_grid) 
      , PrimalSolver<dim>(coarse_grid, 
                          fe, 
                          quadrature, 
                          face_quadrature, 
                          rhs_function, 
                          boundary_values) 
    {} 

    template <int dim> 
    void RefinementGlobal<dim>::refine_grid() 
    { 
      this->triangulation->refine_global(1); 
    } 

    template <int dim> 
    class RefinementKelly : public PrimalSolver<dim> 
    { 
    public: 
      RefinementKelly(Triangulation<dim> &       coarse_grid, 
                      const FiniteElement<dim> & fe, 
                      const Quadrature<dim> &    quadrature, 
                      const Quadrature<dim - 1> &face_quadrature, 
                      const Function<dim> &      rhs_function, 
                      const Function<dim> &      boundary_values); 

      virtual void refine_grid() override; 
    }; 

    template <int dim> 
    RefinementKelly<dim>::RefinementKelly( 
      Triangulation<dim> &       coarse_grid, 
      const FiniteElement<dim> & fe, 
      const Quadrature<dim> &    quadrature, 
      const Quadrature<dim - 1> &face_quadrature, 
      const Function<dim> &      rhs_function, 
      const Function<dim> &      boundary_values) 
      : Base<dim>(coarse_grid) 
      , PrimalSolver<dim>(coarse_grid, 
                          fe, 
                          quadrature, 
                          face_quadrature, 
                          rhs_function, 
                          boundary_values) 
    {} 

    template <int dim> 
    void RefinementKelly<dim>::refine_grid() 
    { 
      Vector<float> estimated_error_per_cell( 
        this->triangulation->n_active_cells()); 
      KellyErrorEstimator<dim>::estimate( 
        this->dof_handler, 
        QGauss<dim - 1>(this->fe->degree + 1), 
        std::map<types::boundary_id, const Function<dim> *>(), 
        this->solution, 
        estimated_error_per_cell); 
      GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.3, 
                                                      0.03); 
      this->triangulation->execute_coarsening_and_refinement(); 
    } 

//  @sect4{The RefinementWeightedKelly class}  

// 这个类是前一个类的变种，它允许通过一些函数来加权我们从库的凯利指标中得到的细化指标。我们包括这个类，因为这个例子程序的目标是展示自动细化标准，即使是复杂的输出量，如点值或应力。如果我们不解决一个对偶问题并计算其中的权重，我们很可能会想给指标一个手工制作的权重，以说明我们要评估这些数量的事实。这个类接受这样一个加权函数作为其构造函数的参数。

    template <int dim> 
    class RefinementWeightedKelly : public PrimalSolver<dim> 
    { 
    public: 
      RefinementWeightedKelly(Triangulation<dim> &       coarse_grid, 
                              const FiniteElement<dim> & fe, 
                              const Quadrature<dim> &    quadrature, 
                              const Quadrature<dim - 1> &face_quadrature, 
                              const Function<dim> &      rhs_function, 
                              const Function<dim> &      boundary_values, 
                              const Function<dim> &      weighting_function); 

      virtual void refine_grid() override; 

    private: 
      const SmartPointer<const Function<dim>> weighting_function; 
    }; 

    template <int dim> 
    RefinementWeightedKelly<dim>::RefinementWeightedKelly( 
      Triangulation<dim> &       coarse_grid, 
      const FiniteElement<dim> & fe, 
      const Quadrature<dim> &    quadrature, 
      const Quadrature<dim - 1> &face_quadrature, 
      const Function<dim> &      rhs_function, 
      const Function<dim> &      boundary_values, 
      const Function<dim> &      weighting_function) 
      : Base<dim>(coarse_grid) 
      , PrimalSolver<dim>(coarse_grid, 
                          fe, 
                          quadrature, 
                          face_quadrature, 
                          rhs_function, 
                          boundary_values) 
      , weighting_function(&weighting_function) 
    {} 

// 现在，这里是主函数，包括加权。

    template <int dim> 
    void RefinementWeightedKelly<dim>::refine_grid() 
    { 

// 首先通过库中已经实现的方法为所有单元计算一些基于残差的误差指标。我们在这里计算的具体内容在该类的文档中会有更详细的描述。

      Vector<float> estimated_error_per_cell( 
        this->triangulation->n_active_cells()); 
      std::map<types::boundary_id, const Function<dim> *> dummy_function_map; 
      KellyErrorEstimator<dim>::estimate(this->dof_handler, 
                                         *this->face_quadrature, 
                                         dummy_function_map, 
                                         this->solution, 
                                         estimated_error_per_cell); 

// 接下来用给与构造函数的值来衡量指标向量中的每个条目，在单元格中心进行评估。我们需要将结果写入对应于当前单元的向量条目中，我们可以通过使用 CellAccessor::active_cell_index(). 询问该单元在所有活动单元中的索引来获得（实际上，对于我们在循环中处理的第一个单元，该索引为0，第二个单元为1，等等，我们也可以使用一个整数计数器来跟踪该索引；但是使用 CellAccessor::active_cell_index() 使其更加明确。

      for (const auto &cell : this->dof_handler.active_cell_iterators()) 
        estimated_error_per_cell(cell->active_cell_index()) *= 
          weighting_function->value(cell->center()); 

      GridRefinement::refine_and_coarsen_fixed_number(*this->triangulation, 
                                                      estimated_error_per_cell, 
                                                      0.3, 
                                                      0.03); 
      this->triangulation->execute_coarsening_and_refinement(); 
    } 

  } // namespace LaplaceSolver 
// @sect3{Equation data}  

// 在这个例子中，我们使用的数据集和前面的一样，但由于可能有人想用不同的边界值和右手函数来运行程序，或者在不同的网格上运行，我们展示了一个简单的技术来做到这一点。为了更加清晰，我们进一步将所有与方程数据有关的东西都打包到一个自己的命名空间中。

// 我们的基本假设是，这是一个研究项目，我们经常有一些测试案例，包括一个域，一个右手边，边界值，可能还有一个指定的系数，以及一些其他参数。当从一个例子转移到另一个例子时，它们常常同时变化。为了使处理这样的问题描述参数集变得简单，是以下的目标。

// 基本上，这个想法是这样的：让我们为每一组数据都有一个结构，在这个结构中，我们把描述一个测试案例的所有东西都打包：这里，这些是两个子类，一个叫 <code>BoundaryValues</code> ，用于精确解的边界值，一个叫 <code>RightHandSide</code>  ，然后是生成粗略网格的方法。由于前面的例子程序的解看起来像弯曲的山脊，所以我们在这里用这个名字来表示包围的类。请注意，两个内层类的名称对于所有包围的测试案例类必须是相同的，同时我们将维度模板参数附加到包围类而不是内层类，以使进一步的处理更简单。 从语言的角度来看，用命名空间来封装这些内部类会比用结构来封装更好。然而，命名空间不能作为模板参数给出，所以我们使用一个结构来允许第二个对象从其给定的参数中选择。当然，这个封闭的结构除了它所声明的类之外，没有任何成员变量，还有一个静态函数来生成粗略的网格；一般来说，它永远不会被实例化）。)

//然后
//想法是如下的（这是正确的时间，也可以简单看看下面的代码）：我们可以为边界值和右手边生成对象，只需将外层类的名字作为模板参数给一个类，我们在这里称之为 <code>Data::SetUp</code> ，然后它为内部类创建对象。在这种情况下，为了得到所有描述弧形山脊解决方案的特征，我们将简单地生成一个 <code>Data::SetUp@<Data::CurvedRidge@></code> 的实例，而我们需要知道的关于该解决方案的一切都将是该对象的静态成员变量和函数。

// 在这种情况下，这种方法可能显得有些多余，但是一旦某种设定不仅有迪里希特边界值和右手函数的特征，而且还有材料属性、诺伊曼值、不同的边界描述符等，就会变得非常方便。在这种情况下， <code>SetUp</code> 类可能由十几个对象组成，而每个描述符类（如下面的 <code>CurvedRidges</code> 类）都必须提供这些对象。然后，你会很高兴，只需在一个地方改变 <code>SetUp</code> 类的模板参数，而不是在很多地方改变，就能从一组数据改变到另一组。

// 有了这个不同测试用例的框架，我们就快完成了，但还有一件事：现在我们可以通过改变一个模板参数，静态地选择要选择的数据集。为了能够动态地做到这一点，即在运行时，我们需要一个基类。我们以明显的方式提供这个基类，见下文，用虚拟抽象函数。这迫使我们引入第二个模板参数 <code>dim</code> ，我们需要这个基类（这可以通过一些模板魔法来避免，但我们省略），但这就是全部。

// 添加新的测试用例现在很简单，你不需要接触框架类，只需要一个类似于  <code>CurvedRidges</code>  的结构。

  namespace Data 
  { 
// @sect4{The SetUpBase and SetUp classes}  

// 基于上述描述， <code>SetUpBase</code> 类就看起来如下。为了允许用这个类来使用 <code>SmartPointer</code> 类，我们从 <code>Subscriptor</code> 类派生出来。

    template <int dim> 
    struct SetUpBase : public Subscriptor 
    { 
      virtual const Function<dim> &get_boundary_values() const = 0; 

      virtual const Function<dim> &get_right_hand_side() const = 0; 

      virtual void 
      create_coarse_grid(Triangulation<dim> &coarse_grid) const = 0; 
    }; 

// 现在是接受模板参数的派生类，如上所述。

// 在这里，我们把数据元素打包成私有变量，并允许通过基类的方法来访问它们。

    template <class Traits, int dim> 
    struct SetUp : public SetUpBase<dim> 
    { 
      virtual const Function<dim> &get_boundary_values() const override; 

      virtual const Function<dim> &get_right_hand_side() const override; 

      virtual void 
      create_coarse_grid(Triangulation<dim> &coarse_grid) const override; 

    private: 
      static const typename Traits::BoundaryValues boundary_values; 
      static const typename Traits::RightHandSide  right_hand_side; 
    }; 

// 我们必须为上述类的静态成员变量提供定义。

    template <class Traits, int dim> 
    const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values; 
    template <class Traits, int dim> 
    const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side; 

// 还有成员函数的定义。

    template <class Traits, int dim> 
    const Function<dim> &SetUp<Traits, dim>::get_boundary_values() const 
    { 
      return boundary_values; 
    } 

    template <class Traits, int dim> 
    const Function<dim> &SetUp<Traits, dim>::get_right_hand_side() const 
    { 
      return right_hand_side; 
    } 

    template <class Traits, int dim> 
    void SetUp<Traits, dim>::create_coarse_grid( 
      Triangulation<dim> &coarse_grid) const 
    { 
      Traits::create_coarse_grid(coarse_grid); 
    } 
// @sect4{The CurvedRidges class}  

// 用于描述 <code>curved ridge</code> 问题的边界值和右手边的类已经在 step-13 示例程序中使用，那么就像这样。

    template <int dim> 
    struct CurvedRidges 
    { 
      class BoundaryValues : public Function<dim> 
      { 
      public: 
        virtual double value(const Point<dim> & p, 
                             const unsigned int component) const;
      }; 

      class RightHandSide : public Function<dim> 
      { 
      public: 
        virtual double value(const Point<dim> & p, 
                             const unsigned int component) const; 
      }; 

      static void create_coarse_grid(Triangulation<dim> &coarse_grid); 
    }; 

    template <int dim> 
    double CurvedRidges<dim>::BoundaryValues::value( 
      const Point<dim> &p, 
      const unsigned int /*component*/) const 
    { 
      double q = p(0); 
      for (unsigned int i = 1; i < dim; ++i) 
        q += std::sin(10 * p(i) + 5 * p(0) * p(0)); 
      const double exponential = std::exp(q); 
      return exponential; 
    } 

    template <int dim> 
    double CurvedRidges<dim>::RightHandSide::value( 
      const Point<dim> &p, 
      const unsigned int /*component*/) const 
    { 
      double q = p(0); 
      for (unsigned int i = 1; i < dim; ++i) 
        q += std::sin(10 * p(i) + 5 * p(0) * p(0)); 
      const double u  = std::exp(q); 
      double       t1 = 1, t2 = 0, t3 = 0; 
      for (unsigned int i = 1; i < dim; ++i) 
        { 
          t1 += std::cos(10 * p(i) + 5 * p(0) * p(0)) * 10 * p(0); 
          t2 += 10 * std::cos(10 * p(i) + 5 * p(0) * p(0)) - 
                100 * std::sin(10 * p(i) + 5 * p(0) * p(0)) * p(0) * p(0); 
          t3 += 100 * std::cos(10 * p(i) + 5 * p(0) * p(0)) * 
                  std::cos(10 * p(i) + 5 * p(0) * p(0)) - 
                100 * std::sin(10 * p(i) + 5 * p(0) * p(0)); 
        } 
      t1 = t1 * t1; 

      return -u * (t1 + t2 + t3); 
    } 

    template <int dim> 
    void CurvedRidges<dim>::create_coarse_grid(Triangulation<dim> &coarse_grid) 
    { 
      GridGenerator::hyper_cube(coarse_grid, -1, 1); 
      coarse_grid.refine_global(2); 
    } 
// @sect4{The Exercise_2_3 class}  

// 这个例子程序是在为自适应有限元方法和基于对偶性的误差估计的讲座提供实践课程时写的。在这些课程中，我们有一个练习，要求在一个中心有方孔的正方形域上求解右方恒定的拉普拉斯方程，并且边界值为零。由于这个问题的属性在这里的实现特别简单，所以让我们来做。由于练习的编号是2.3，所以我们也擅自为这个类保留这个名称。

    template <int dim> 
    struct Exercise_2_3 
    { 

// 我们需要一个类来表示问题的边界值。在这种情况下，这很简单：它是零函数，所以甚至不需要声明一个类，只需要一个别名。

      using BoundaryValues = Functions::ZeroFunction<dim>; 

// 第二，一个表示右边的类。因为它们是常数，所以只要把库中相应的类子类化就可以了。

      class RightHandSide : public Functions::ConstantFunction<dim> 
      { 
      public: 
        RightHandSide() 
          : Functions::ConstantFunction<dim>(1.) 
        {} 
      }; 

// 最后是一个生成粗略网格的函数。这里的情况有些复杂，请看下面的内容。

      static void create_coarse_grid(Triangulation<dim> &coarse_grid); 
    }; 

// 如上所述，本例的网格是正方形[-1,1]^2，其中的正方形[-1/2,1/2]^2为孔。我们将粗略的网格创建为4乘以4的单元，中间的四个单元缺失。要了解网格的具体样子，最简单的方法可能是先看一下本教程程序的 "结果 "部分。一般来说，如果你想了解更多关于创建网格的信息，无论是像我们在这里所做的那样从头开始，还是使用其他技术，你都应该看一下  step-49  。

// 当然，这个例子可以扩展到3d，但是由于这个函数不能以独立于维度的方式编写，我们选择不在这里实现，而是只对dim=2的模板进行专业化处理。如果你编译3d的程序，你会从链接器中得到一个信息，即这个函数在3d中没有实现，需要提供。

// 对于这个几何体的创建，库中没有预定义的方法。在这种情况下，几何体还是很简单的，可以用手来创建，而不是用网格发生器。

    template <> 
    void Exercise_2_3<2>::create_coarse_grid(Triangulation<2> &coarse_grid) 
    { 

// 我们首先定义空间维度，以便让函数中那些实际上与维度无关的部分使用这个变量。这样，如果你以后把它作为一个起点来实现这个网格的三维版本，就会更简单。下一步是要有一个顶点的列表。这里，它们是24个（5乘以5，中间的省略）。最好的办法是在这里画一个草图。

      const unsigned int dim = 2; 

      const std::vector<Point<2>> vertices = { 
        {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, // 
        {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, // 
        {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               // 
        {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, // 
        {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}}; 

// 接下来，我们要定义单元格和它们所包含的顶点。

      const std::vector<std::array<int, GeometryInfo<dim>::vertices_per_cell>> 
        cell_vertices = {{{0, 1, 5, 6}}, 
                         {{1, 2, 6, 7}}, 
                         {{2, 3, 7, 8}}, 
                         {{3, 4, 8, 9}}, 
                         {{5, 6, 10, 11}}, 
                         {{8, 9, 12, 13}}, 
                         {{10, 11, 14, 15}}, 
                         {{12, 13, 17, 18}}, 
                         {{14, 15, 19, 20}}, 
                         {{15, 16, 20, 21}}, 
                         {{16, 17, 21, 22}}, 
                         {{17, 18, 22, 23}}}; 

      const unsigned int n_cells = cell_vertices.size(); 

// 我们再次从中生成一个C++向量类型，但这次是通过在单元格上循环来实现的（是的，这很无聊）。此外，我们将所有单元格的材料指标设置为零。

      std::vector<CellData<dim>> cells(n_cells, CellData<dim>()); 
      for (unsigned int i = 0; i < n_cells; ++i) 
        { 
          for (unsigned int j = 0; j < cell_vertices[i].size(); ++j) 
            cells[i].vertices[j] = cell_vertices[i][j]; 
          cells[i].material_id = 0; 
        } 

// 最后将所有这些信息传递给库，以生成一个三角图。最后一个参数可以用来将三角形的某些面的非零边界指标的信息传递给库，但我们在这里不希望这样，所以我们给出一个空对象。

      coarse_grid.create_triangulation(vertices, cells, SubCellData()); 

// 因为我们希望本例中的评估点（3/4,3/4）是一个网格点，所以我们在全局范围内细化一次。

      coarse_grid.refine_global(1); 
    } 
  } // namespace Data 
// @sect4{Discussion}  

// 你现在已经读完了这个框架，你可能会想，为什么我们没有选择直接把实现某种设置的类（比如 <code>CurvedRidges</code> 类）作为派生自 <code>Data::SetUpBase</code> 的类来实现。事实上，我们可以很好地这样做。唯一的原因是，这样我们就必须在 <code>CurvedRidges</code> 类中为解决方案和右手边的类设置成员变量，以及重载基类的抽象函数来访问这些成员变量的成员函数。 <code>SetUp</code> 类的唯一原因是让我们不必再重申这些成员变量和函数，而这些成员变量和函数在所有这些类中都是必要的。在某种程度上，这里的模板机制只是提供了一种方法，为一些依赖于外部量的函数提供默认的实现，因此不能使用正常的虚拟函数来提供，至少在没有模板的帮助下不能。

// 然而，可能有很好的理由来实际实现从 <code>Data::SetUpBase</code> 派生的类，例如，如果解或右手边的类需要带参数的构造函数，而 <code>Data::SetUpBase</code> 类无法提供。在这种情况下，子类化是一个值得考虑的策略。对于特殊情况的其他可能性是，从 <code>Data::SetUp@<SomeSetUp@></code> 派生，其中 <code>SomeSetUp</code> 表示一个类，或者甚至明确地专门化 <code>Data::SetUp@<SomeSetUp@></code>  。后者允许透明地使用 <code>SetUp</code> 类用于其他设置的方式，但对特殊参数采取特殊行动。

// 赞成这里采取的方法的最后一个意见是：我们无数次发现，当开始一个项目时，参数的数量（通常是边界值，右侧，粗略的网格，就像这里）很小，测试案例的数量也很小。然后，人们一开始就把它们手工编码成一些 <code>switch</code> 的语句。随着时间的推移，项目的增长，测试用例的数量也在增长。 <code>switch</code> 语句的数量也随之增长，它们的长度也是如此，人们开始想办法考虑不可能的例子，其中域、边界值和右手边不再适合在一起，并且开始失去对整个结构的概述。事实证明，把属于某个测试用例的所有东西都封装到一个自己的结构中是值得的，因为它把属于一个测试用例的所有东西都放在一个地方。此外，它允许把这些东西都放在一个或多个文件中，这些文件只用于测试用例和它们的数据，而不需要把它们的实际实现与程序的其他部分联系起来。

//  @sect3{Dual functionals}  

// 和程序的其他部分一样，我们把所有需要描述对偶函数的东西放到一个自己的命名空间中，并定义一个抽象的基类，提供解决对偶问题的类在工作中需要的接口。

// 然后，我们将实现两个这样的类，用于评估一个点的值和该点的解的导数。对于这些函数，我们已经有了相应的评估对象，所以它们是互补的。

  namespace DualFunctional 
  { 
// @sect4{The DualFunctionalBase class}  

// 首先从对偶函数的基类开始。因为对于线性问题来说，对偶问题的特征只在右手边起作用，所以我们只需要提供一个函数来组装给定离散化的右手边。

    template <int dim> 
    class DualFunctionalBase : public Subscriptor 
    { 
    public: 
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                Vector<double> &       rhs) const = 0; 
    }; 
// @sect4{The dual functional PointValueEvaluation class}  

// 作为第一个应用，我们考虑对应于在一个给定的点上评估解决方案的值的函数，我们再次假设该点为一个顶点。除了接受和存储评估点的构造函数之外，这个类只包括实现组装右手边的函数。

    template <int dim> 
    class PointValueEvaluation : public DualFunctionalBase<dim> 
    { 
    public: 
      PointValueEvaluation(const Point<dim> &evaluation_point); 

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                Vector<double> &       rhs) const override; 

      DeclException1( 
        ExcEvaluationPointNotFound, 
        Point<dim>, 
        << "The evaluation point " << arg1 
        << " was not found among the vertices of the present grid."); 

    protected: 
      const Point<dim> evaluation_point; 
    }; 

    template <int dim> 
    PointValueEvaluation<dim>::PointValueEvaluation( 
      const Point<dim> &evaluation_point) 
      : evaluation_point(evaluation_point) 
    {} 

// 至于做这门课的主要目的，组装右手边，让我们首先考虑什么是必要的。对偶问题的右手边是一个值的向量J(phi_i)，其中J是误差函数，phi_i是第i个形状函数。这里，J是在点x0处的评价，即J(phi_i)=phi_i(x0)。

// 现在，我们已经假定评价点是一个顶点。因此，对于我们在这个程序中可能使用的通常的有限元来说，我们可以想当然地认为在这样一个点上正好有一个形状函数是不为零的，特别是其值为1。因此，我们将右手边的向量设置为全零，然后寻找与该点相关的形状函数，并将右手边向量的相应值设置为1。

    template <int dim> 
    void 
    PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                            Vector<double> &       rhs) const 
    { 

// 所以，首先把所有东西都设为零......

      rhs.reinit(dof_handler.n_dofs()); 

// ...然后在单元格上循环，在顶点中找到评估点（或者非常接近顶点，由于浮点舍入，可能会出现这种情况）。

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto vertex : cell->vertex_indices()) 
          if (cell->vertex(vertex).distance(evaluation_point) < 
              cell->diameter() * 1e-8) 
            { 

// 好的，找到了，所以设置相应的条目，然后离开函数，因为我们已经完成了。

              rhs(cell->vertex_dof_index(vertex, 0)) = 1; 
              return; 
            } 

// 最后，一个理智的检查：如果我们以某种方式来到这里，那么我们一定是错过了评估点，所以无条件地引发一个异常。

      AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point)); 
    } 
// @sect4{The dual functional PointXDerivativeEvaluation class}  

// 作为第二个应用，我们再次考虑在一个点上对解决方案的x导数进行评估。同样，这个类的声明和它的构造函数的实现也不是太有趣。

    template <int dim> 
    class PointXDerivativeEvaluation : public DualFunctionalBase<dim> 
    { 
    public: 
      PointXDerivativeEvaluation(const Point<dim> &evaluation_point); 

      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                Vector<double> &       rhs) const; 

      DeclException1( 
        ExcEvaluationPointNotFound, 
        Point<dim>, 
        << "The evaluation point " << arg1 
        << " was not found among the vertices of the present grid."); 

    protected: 
      const Point<dim> evaluation_point; 
    }; 

    template <int dim> 
    PointXDerivativeEvaluation<dim>::PointXDerivativeEvaluation( 
      const Point<dim> &evaluation_point) 
      : evaluation_point(evaluation_point) 
    {} 

//有趣的是这个函数的实现：这里，J(phi_i)=d/dx phi_i(x0)。

// 我们可以像实现各自的评价对象那样，在这个评价点上取每个形状函数phi_i的梯度的平均值。然而，我们采取了一个略微不同的方法：我们简单地取该点周围所有单元格的平均值。哪些单元 <code>surrounds</code> 是评估点，这个问题取决于网格宽度，包括那些单元的中点到评估点的距离小于单元的直径的单元。

// 在这些单元的面积/体积上取梯度的平均值，可以得到一个与梯度的点评估非常接近的对偶解。从理论上讲，这并没有明显改变方法，这一点很简单。

    template <int dim> 
    void PointXDerivativeEvaluation<dim>::assemble_rhs( 
      const DoFHandler<dim> &dof_handler, 
      Vector<double> &       rhs) const 
    { 

// 同样，首先将所有条目设置为零。

      rhs.reinit(dof_handler.n_dofs()); 

// 用正交公式初始化一个 <code>FEValues</code> 对象，有正交点数量和形状函数的缩写......

      QGauss<dim>        quadrature(dof_handler.get_fe().degree + 1); 
      FEValues<dim>      fe_values(dof_handler.get_fe(), 
                              quadrature, 
                              update_gradients | update_quadrature_points | 
                                update_JxW_values); 
      const unsigned int n_q_points    = fe_values.n_quadrature_points; 
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell; 

// ...并有两个对象用于存储单元上自由度的全局指数，以及正交点上形状函数的梯度值。

      Vector<double>            cell_rhs(dofs_per_cell); 
      std::vector<unsigned int> local_dof_indices(dofs_per_cell); 

// 最后有一个变量，我们将通过对这些单元上的单位函数进行积分，总结出这些单元的面积/体积。

      double total_volume = 0; 

// 然后在所有单元上开始循环，并选择那些与评估点足够接近的单元。

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (cell->center().distance(evaluation_point) <= cell->diameter()) 
          { 

// 如果我们找到了这样的单元，那么就初始化 <code>FEValues</code> 对象，并整合每个形状函数梯度的x分量，以及总面积/体积的单位函数。

            fe_values.reinit(cell); 
            cell_rhs = 0; 

            for (unsigned int q = 0; q < n_q_points; ++q) 
              { 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  cell_rhs(i) += 
                    fe_values.shape_grad(i, q)[0] // (d/dx phi_i(x_q)) 
                    * fe_values.JxW(q);           // * dx 
                total_volume += fe_values.JxW(q); 
              } 

// 如果我们有本地的贡献，把它们分配到全局矢量。

            cell->get_dof_indices(local_dof_indices); 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              rhs(local_dof_indices[i]) += cell_rhs(i); 
          } 

// 在我们循环了所有的单元格之后，检查我们是否找到了任何单元格，确保其体积不为零。如果不是，那么结果将是错误的，因为这时的右边应该仍然是零，所以抛出一个异常。

      AssertThrow(total_volume > 0, 
                  ExcEvaluationPointNotFound(evaluation_point)); 

// 最后，我们现在只整合了形状函数的梯度，而没有取其平均值。我们通过除以我们所积分的体积的大小来解决这个问题。

      rhs /= total_volume; 
    } 

  } // namespace DualFunctional 
// @sect3{Extending the LaplaceSolver namespace}  
  namespace LaplaceSolver 
  { 
// @sect4{The DualSolver class}  

// 与上面的  <code>PrimalSolver</code>  类相同，我们现在实现一个  <code>DualSolver</code>  。它具有所有相同的特征，唯一的区别是它不接受一个表示右侧对象的函数对象，而现在接受一个 <code>DualFunctionalBase</code> 对象，它将集合对偶问题的右侧向量。这个类的其余部分是相当琐碎的。

// 由于原始求解器和对偶求解器将使用相同的三角形，但是不同的离散化，现在很清楚为什么我们将 <code>Base</code> 类变成了虚拟类：因为最终类将从 <code>PrimalSolver</code> 以及 <code>DualSolver</code>, it would have two <code>Base</code> 实例中派生，我们是不是应该将继承标记为虚拟。因为在许多应用中，基类会存储更多的信息，而不仅仅是需要在原始求解器和对偶求解器之间共享的三角形，所以我们通常不希望使用两个这样的基类。

    template <int dim> 
    class DualSolver : public Solver<dim> 
    { 
    public: 
      DualSolver( 
        Triangulation<dim> &                           triangulation, 
        const FiniteElement<dim> &                     fe, 
        const Quadrature<dim> &                        quadrature, 
        const Quadrature<dim - 1> &                    face_quadrature, 
        const DualFunctional::DualFunctionalBase<dim> &dual_functional); 

    protected: 
      const SmartPointer<const DualFunctional::DualFunctionalBase<dim>> 
                   dual_functional; 
      virtual void assemble_rhs(Vector<double> &rhs) const override; 

      static const Functions::ZeroFunction<dim> boundary_values; 
    }; 

    template <int dim> 
    const Functions::ZeroFunction<dim> DualSolver<dim>::boundary_values; 

    template <int dim> 
    DualSolver<dim>::DualSolver( 
      Triangulation<dim> &                           triangulation, 
      const FiniteElement<dim> &                     fe, 
      const Quadrature<dim> &                        quadrature, 
      const Quadrature<dim - 1> &                    face_quadrature, 
      const DualFunctional::DualFunctionalBase<dim> &dual_functional) 
      : Base<dim>(triangulation) 
      , Solver<dim>(triangulation, 
                    fe, 
                    quadrature, 
                    face_quadrature, 
                    boundary_values) 
      , dual_functional(&dual_functional) 
    {} 

    template <int dim> 
    void DualSolver<dim>::assemble_rhs(Vector<double> &rhs) const 
    { 
      dual_functional->assemble_rhs(this->dof_handler, rhs); 
    } 
// @sect4{The WeightedResidual class}  

// 这里终于出现了这个程序的主类，也就是实现双重加权残差估计器的类。它连接了原始和对偶求解器类，用于计算原始和对偶解，并实现了误差表示公式，用于误差估计和网格细化。

// 该类的前几个函数大多是对基类各自函数的覆盖。

    template <int dim> 
    class WeightedResidual : public PrimalSolver<dim>, public DualSolver<dim> 
    { 
    public: 
      WeightedResidual( 
        Triangulation<dim> &                           coarse_grid, 
        const FiniteElement<dim> &                     primal_fe, 
        const FiniteElement<dim> &                     dual_fe, 
        const Quadrature<dim> &                        quadrature, 
        const Quadrature<dim - 1> &                    face_quadrature, 
        const Function<dim> &                          rhs_function, 
        const Function<dim> &                          boundary_values, 
        const DualFunctional::DualFunctionalBase<dim> &dual_functional); 

      virtual void solve_problem() override; 

      virtual void postprocess( 
        const Evaluation::EvaluationBase<dim> &postprocessor) const override; 

      virtual unsigned int n_dofs() const override; 

      virtual void refine_grid() override; 

      virtual void output_solution() const override; 

    private: 

// 在私有部分，我们有两个函数，用来调用原始类和双基类的 <code>solve_problem</code> 函数。这两个函数将被本类的 <code>solve_problem</code> 函数所平行调用。

      void solve_primal_problem(); 
      void solve_dual_problem(); 

// 然后声明活动单元迭代器的缩写，以避免我们不得不重复写这个冗长的名字。

      using active_cell_iterator = 
        typename DoFHandler<dim>::active_cell_iterator; 

// 接下来，声明一个数据类型，我们将用它来存储面对误差估计器的贡献。我们的想法是，我们可以计算从两个单元格中的每一个到这个面的脸部条款，因为从两边看时它们是一样的。我们要做的是，根据下面解释的一些规则，只计算一次，由相邻的两个单元负责计算。然后，我们将每个面的贡献存储在一个映射面与它们的值的地图中，并通过第二次在单元格上循环并从地图上抓取值来收集每个单元格的贡献。

// 这个地图的数据类型在此声明。

      using FaceIntegrals = 
        typename std::map<typename DoFHandler<dim>::face_iterator, double>; 

// 在计算单元和面的误差估计时，我们需要一些辅助对象，例如 <code>FEValues</code> 和 <code>FEFaceValues</code> 函数，但也需要一些临时对象来存储原始和对偶解的值和梯度，例如。这些字段在三个函数中都是需要的，这些函数分别在单元格、规则面和不规则面上做积分。

// 有三种合理的方式来提供这些字段：第一，作为需要它们的函数中的局部变量；第二，作为本类的成员变量；第三，作为参数传递给该函数。

// 这三种方式都有缺点：第三种是它们的数量不可忽略，会使调用这些函数成为一项漫长的事业。第二种方法的缺点是不允许并行化，因为计算错误估计值的线程必须各自拥有这些变量的副本，所以包围类的成员变量将不起作用。第一种方法虽然直接，但有一个微妙但重要的缺点：我们会反复调用这些函数，也许是成千上万次；现在证明，从堆中分配向量和其他需要内存的对象在运行时间上是很昂贵的，因为当涉及到几个线程时，内存分配很昂贵。因此，只分配一次内存，并尽可能频繁地回收这些对象是明显更好的做法。

// 该怎么做呢？我们的答案是使用第三种策略的一个变种。事实上，这正是WorkStream概念所要做的（我们已经在上面介绍了它，但也可以参见 @ref threads ）。为了避免我们必须给这些函数十几个参数，我们将所有这些变量打包成两个结构，一个用于单元格的计算，另一个用于面的计算。然后，这两个结构被加入到WeightedResidualScratchData类中，该类将作为WorkStream概念的 "划痕数据 "类。

      struct CellData 
      { 
        FEValues<dim>                           fe_values; 
        const SmartPointer<const Function<dim>> right_hand_side; 

        std::vector<double> cell_residual; 
        std::vector<double> rhs_values; 
        std::vector<double> dual_weights; 
        std::vector<double> cell_laplacians; 
        CellData(const FiniteElement<dim> &fe, 
                 const Quadrature<dim> &   quadrature, 
                 const Function<dim> &     right_hand_side); 
        CellData(const CellData &cell_data); 
      }; 

      struct FaceData 
      { 
        FEFaceValues<dim>    fe_face_values_cell; 
        FEFaceValues<dim>    fe_face_values_neighbor; 
        FESubfaceValues<dim> fe_subface_values_cell; 

        std::vector<double>                  jump_residual; 
        std::vector<double>                  dual_weights; 
        typename std::vector<Tensor<1, dim>> cell_grads; 
        typename std::vector<Tensor<1, dim>> neighbor_grads;
        FaceData(const FiniteElement<dim> & fe, 
                 const Quadrature<dim - 1> &face_quadrature); 
        FaceData(const FaceData &face_data); 
      }; 

      struct WeightedResidualScratchData 
      { 
        WeightedResidualScratchData(
          const FiniteElement<dim> & primal_fe,
          const Quadrature<dim> &    primal_quadrature, 
          const Quadrature<dim - 1> &primal_face_quadrature, 
          const Function<dim> &      rhs_function, 
          const Vector<double> &     primal_solution, 
          const Vector<double> &     dual_weights); 

        WeightedResidualScratchData( 
          const WeightedResidualScratchData &scratch_data); 

        CellData       cell_data; 
        FaceData       face_data; 
        Vector<double> primal_solution;
        Vector<double> dual_weights;
      }; 
// WorkStream::run 一般要有一个从头开始的对象和一个拷贝对象。在这里，由于与我们在 step-9 中讨论梯度的近似计算时类似的原因，我们实际上不需要一个 "拷贝数据 "结构。既然WorkStream坚持要有一个这样的结构，我们就声明一个空的结构，除了存在之外什么都不做。

      struct WeightedResidualCopyData 
      {}; 

// 关于误差估计器的评估，我们有一个驱动函数，使用  WorkStream::run()  在每个单元上调用第二个函数。

      void estimate_error(Vector<float> &error_indicators) const; 

      void estimate_on_one_cell(const active_cell_iterator & cell, 
                                WeightedResidualScratchData &scratch_data, 
                                WeightedResidualCopyData &   copy_data, 
                                Vector<float> &              error_indicators, 
                                FaceIntegrals &face_integrals) const; 

// 然后我们有函数对误差表示公式进行实际积分。它们将分别处理单元格内部、没有悬挂节点的面和有悬挂节点的面的条款。

      void integrate_over_cell(const active_cell_iterator &cell, 
                               const Vector<double> &      primal_solution, 
                               const Vector<double> &      dual_weights, 
                               CellData &                  cell_data, 
                               Vector<float> &error_indicators) const; 

      void integrate_over_regular_face(const active_cell_iterator &cell, 
                                       const unsigned int          face_no, 
                                       const Vector<double> &primal_solution, 
                                       const Vector<double> &dual_weights, 
                                       FaceData &            face_data, 
                                       FaceIntegrals &face_integrals) const;  
      void integrate_over_irregular_face(const active_cell_iterator &cell, 
                                         const unsigned int          face_no, 
                                         const Vector<double> &primal_solution, 
                                         const Vector<double> &dual_weights, 
                                         FaceData &            face_data,  
                                         FaceIntegrals &face_integrals) const; 
    }; 

// 在这个类的实现中，我们首先有 <code>CellData</code> and <code>FaceData</code> 成员类的构造函数，以及 <code>WeightedResidual</code> 构造函数。它们只将字段初始化为正确的长度，所以我们不必过多地讨论它们。

    template <int dim> 
    WeightedResidual<dim>::CellData::CellData( 
      const FiniteElement<dim> &fe, 
      const Quadrature<dim> &   quadrature, 
      const Function<dim> &     right_hand_side) 
      : fe_values(fe, 
                  quadrature, 
                  update_values | update_hessians | update_quadrature_points | 
                    update_JxW_values)  
      , right_hand_side(&right_hand_side) 
      , cell_residual(quadrature.size()) 
      , rhs_values(quadrature.size()) 
      , dual_weights(quadrature.size()) 
      , cell_laplacians(quadrature.size()) 
    {} 

    template <int dim> 
    WeightedResidual<dim>::CellData::CellData(const CellData &cell_data) 
      : fe_values(cell_data.fe_values.get_fe(), 
                  cell_data.fe_values.get_quadrature(), 
                  update_values | update_hessians | update_quadrature_points | 
                    update_JxW_values) 
      , right_hand_side(cell_data.right_hand_side) 
      , cell_residual(cell_data.cell_residual) 
      , rhs_values(cell_data.rhs_values) 
      , dual_weights(cell_data.dual_weights) 
      , cell_laplacians(cell_data.cell_laplacians) 
    {} 

    template <int dim> 
    WeightedResidual<dim>::FaceData::FaceData( 
      const FiniteElement<dim> & fe, 
      const Quadrature<dim - 1> &face_quadrature) 
      : fe_face_values_cell(fe, 
                            face_quadrature, 
                            update_values | update_gradients | 
                              update_JxW_values | update_normal_vectors) 
      , fe_face_values_neighbor(fe, 
                                face_quadrature, 
                                update_values | update_gradients | 
                                  update_JxW_values | update_normal_vectors) 
      , fe_subface_values_cell(fe, face_quadrature, update_gradients) 
    { 
      const unsigned int n_face_q_points = face_quadrature.size(); 

      jump_residual.resize(n_face_q_points); 
      dual_weights.resize(n_face_q_points); 
      cell_grads.resize(n_face_q_points); 
      neighbor_grads.resize(n_face_q_points); 
    } 

    template <int dim> 
    WeightedResidual<dim>::FaceData::FaceData(const FaceData &face_data) 
      : fe_face_values_cell(face_data.fe_face_values_cell.get_fe(), 
                            face_data.fe_face_values_cell.get_quadrature(), 
                            update_values | update_gradients | 
                              update_JxW_values | update_normal_vectors) 
      , fe_face_values_neighbor( 
          face_data.fe_face_values_neighbor.get_fe(), 
          face_data.fe_face_values_neighbor.get_quadrature(), 
          update_values | update_gradients | update_JxW_values | 
            update_normal_vectors) 
      , fe_subface_values_cell( 
          face_data.fe_subface_values_cell.get_fe(), 
          face_data.fe_subface_values_cell.get_quadrature(), 
          update_gradients)  
      , jump_residual(face_data.jump_residual) 
      , dual_weights(face_data.dual_weights) 
      , cell_grads(face_data.cell_grads) 
      , neighbor_grads(face_data.neighbor_grads) 
    {} 

    template <int dim> 
    WeightedResidual<dim>::WeightedResidualScratchData::  
      WeightedResidualScratchData( 
        const FiniteElement<dim> & primal_fe, 
        const Quadrature<dim> &    primal_quadrature, 
        const Quadrature<dim - 1> &primal_face_quadrature, 
        const Function<dim> &      rhs_function, 
        const Vector<double> &     primal_solution, 
        const Vector<double> &     dual_weights) 
      : cell_data(primal_fe, primal_quadrature, rhs_function) 
      , face_data(primal_fe, primal_face_quadrature) 
      , primal_solution(primal_solution) 
      , dual_weights(dual_weights) 
    {} 

    template <int dim> 
    WeightedResidual<dim>::WeightedResidualScratchData:: 
      WeightedResidualScratchData( 
        const WeightedResidualScratchData &scratch_data) 
      : cell_data(scratch_data.cell_data) 
      , face_data(scratch_data.face_data) 
      , primal_solution(scratch_data.primal_solution) 
      , dual_weights(scratch_data.dual_weights) 
    {} 

    template <int dim> 
    WeightedResidual<dim>::WeightedResidual( 
      Triangulation<dim> &                           coarse_grid, 
      const FiniteElement<dim> &                     primal_fe, 
      const FiniteElement<dim> &                     dual_fe, 
      const Quadrature<dim> &                        quadrature, 
      const Quadrature<dim - 1> &                    face_quadrature,  
      const Function<dim> &                          rhs_function, 
      const Function<dim> &                          bv, 
      const DualFunctional::DualFunctionalBase<dim> &dual_functional) 
      : Base<dim>(coarse_grid) 
      , PrimalSolver<dim>(coarse_grid, 
                          primal_fe, 
                          quadrature, 
                          face_quadrature, 
                          rhs_function, 
                          bv) 
      , DualSolver<dim>(coarse_grid, 
                        dual_fe, 
                        quadrature, 
                        face_quadrature, 
                        dual_functional) 
    {} 

// 接下来的五个函数很无聊，因为它们只是简单地将它们的工作传递给基类。第一个函数并行地调用原始和对偶求解器，而解的后处理和检索自由度的数量则由原始类完成。

    template <int dim> 
    void WeightedResidual<dim>::solve_problem() 
    { 
      Threads::TaskGroup<void> tasks; 
      tasks += 
        Threads::new_task(&WeightedResidual<dim>::solve_primal_problem, *this); 
      tasks += 
        Threads::new_task(&WeightedResidual<dim>::solve_dual_problem, *this); 
      tasks.join_all(); 
    } 

    template <int dim> 
    void WeightedResidual<dim>::solve_primal_problem() 
    { 
      PrimalSolver<dim>::solve_problem(); 
    } 

    template <int dim> 
    void WeightedResidual<dim>::solve_dual_problem() 
    { 
      DualSolver<dim>::solve_problem(); 
    } 

    template <int dim> 
    void WeightedResidual<dim>::postprocess( 
      const Evaluation::EvaluationBase<dim> &postprocessor) const 
    { 
      PrimalSolver<dim>::postprocess(postprocessor); 
    } 

    template <int dim> 
    unsigned int WeightedResidual<dim>::n_dofs() const 
    { 
      return PrimalSolver<dim>::n_dofs(); 
    } 

// 现在，变得更加有趣了： <code>refine_grid()</code> 函数要求误差估计器计算单元格的误差指标，然后使用其绝对值进行网格细化。

    template <int dim> 
    void WeightedResidual<dim>::refine_grid() 
    { 

// 首先调用计算单元格和全局误差的函数。

      Vector<float> error_indicators(this->triangulation->n_active_cells()); 
      estimate_error(error_indicators); 

//然后
//注意，只有当所有的指标都是正数时，标记单元的细化或粗化才会起作用，以便于它们的比较。因此，去掉所有这些指标上的符号。

      for (float &error_indicator : error_indicators) 
        error_indicator = std::fabs(error_indicator); 

// 最后，我们可以选择不同的细化策略。这里默认的是细化那些误差指标最大、占总误差80%的单元格，而我们粗化那些指标最小、占总误差2%的单元格。

      GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation, 
                                                        error_indicators, 
                                                        0.8, 
                                                        0.02); 
      this->triangulation->execute_coarsening_and_refinement(); 
    } 

// 由于我们想同时输出原始解和对偶解，我们重载了 <code>output_solution</code> 函数。这个函数唯一有趣的特点是，原始解和对偶解是在不同的有限元空间上定义的，这不是 <code>DataOut</code> 类所期望的格式。因此，我们必须将它们转移到一个共同的有限元空间。由于我们只想从质量上看到这些解，所以我们要争夺将对偶解内插到（较小的）原始空间。对于插值，有一个库函数，它接收一个AffineConstraints对象，包括悬挂节点约束。其余的都是标准的。

    template <int dim> 
    void WeightedResidual<dim>::output_solution() const 
    { 
      AffineConstraints<double> primal_hanging_node_constraints; 
      DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler, 
                                              primal_hanging_node_constraints); 
      primal_hanging_node_constraints.close(); 
      Vector<double> dual_solution(PrimalSolver<dim>::dof_handler.n_dofs()); 
      FETools::interpolate(DualSolver<dim>::dof_handler, 
                           DualSolver<dim>::solution, 
                           PrimalSolver<dim>::dof_handler, 
                           primal_hanging_node_constraints, 
                           dual_solution); 

      DataOut<dim> data_out; 
      data_out.attach_dof_handler(PrimalSolver<dim>::dof_handler); 

// 添加我们想要输出的数据向量。两个都加上， <code>DataOut</code> 函数可以处理你想写到输出的多少个数据向量。

      data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution"); 
      data_out.add_data_vector(dual_solution, "dual_solution"); 

      data_out.build_patches(); 

      std::ofstream out("solution-" + std::to_string(this->refinement_cycle) + 
                        ".vtu"); 
      data_out.write(out, DataOutBase::vtu); 
    } 
// @sect3{Estimating errors}  
// @sect4{Error estimation driver functions}  

// 至于误差估计的实际计算，让我们从驱动这一切的函数开始，即调用那些真正做工作的函数，并最终收集结果。

    template <int dim> 
    void 
    WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const 
    { 

// 计算误差的第一个任务是建立向量，表示原始解，以及权重(z-z_h)=(z-I_hz)，两者都在我们已经计算出对偶解的有限元空间。为此，我们必须将原始解内插到对偶有限元空间，并将计算出的对偶解内插到原始有限元空间。幸运的是，库中提供了插值到更大或更小的有限元空间的函数，所以这一点是很明显的。

// 首先，让我们对原始解进行插值：它被逐格插值到我们已经解决了对偶问题的有限元空间中：但是，还是和 <code>WeightedResidual::output_solution</code> 函数一样，我们首先需要创建一个包括悬挂节点约束的AffineConstraints对象，但这次是对偶有限元空间的。

      AffineConstraints<double> dual_hanging_node_constraints; 
      DoFTools::make_hanging_node_constraints(DualSolver<dim>::dof_handler, 
                                              dual_hanging_node_constraints); 
      dual_hanging_node_constraints.close(); 
      Vector<double> primal_solution(DualSolver<dim>::dof_handler.n_dofs()); 
      FETools::interpolate(PrimalSolver<dim>::dof_handler, 
                           PrimalSolver<dim>::solution, 
                           DualSolver<dim>::dof_handler, 
                           dual_hanging_node_constraints, 
                           primal_solution); 
//然后
//为了计算数值逼近的对偶解z插值到原始解的有限元空间并从z中减去：使用 <code>interpolate_difference</code> 函数，在对偶解的元素空间中得到(z-I_hz)。

      AffineConstraints<double> primal_hanging_node_constraints; 
      DoFTools::make_hanging_node_constraints(PrimalSolver<dim>::dof_handler, 
                                              primal_hanging_node_constraints); 
      primal_hanging_node_constraints.close(); 
      Vector<double> dual_weights(DualSolver<dim>::dof_handler.n_dofs()); 
      FETools::interpolation_difference(DualSolver<dim>::dof_handler, 
                                        dual_hanging_node_constraints, 
                                        DualSolver<dim>::solution, 
                                        PrimalSolver<dim>::dof_handler, 
                                        primal_hanging_node_constraints, 
                                        dual_weights); 

// 请注意，这可能会更有效，因为这些约束条件在之前为原始问题组装矩阵和右手边以及写出对偶解时已经使用过了。我们把这方面的程序优化作为一个练习。

// 在计算了对偶权重之后，我们现在开始计算原始解的单元和面的残差。首先，我们在面的迭代器和它们的面的跳跃项对误差估计器的贡献之间建立一个映射。原因是我们只计算了一次跳跃项，从面的一侧开始，并且希望在第二次循环所有单元时才收集它们。

// 我们已经用一个-1e20的值初始化了这个地图，因为如果出现问题，我们因为某些原因无法计算某个面的值，这个值就会在结果中显示出来。其次，这个初始化已经让 std::map 对象分配了它可能需要的所有对象。这一点很重要，因为我们将从并行线程写进这个结构，如果地图需要分配内存，从而重塑其数据结构，那么这样做就不是线程安全的。换句话说，初始化使我们不必在线程每次写入（和修改）该地图的结构时通过互斥来同步。

      FaceIntegrals face_integrals; 
      for (const auto &cell : 
           DualSolver<dim>::dof_handler.active_cell_iterators()) 
        for (const auto &face : cell->face_iterators()) 
          face_integrals[face] = -1e20; 

      auto worker = [this, 
                     &error_indicators, 
                     &face_integrals](const active_cell_iterator & cell, 
                                      WeightedResidualScratchData &scratch_data, 
                                      WeightedResidualCopyData &   copy_data) { 
        this->estimate_on_one_cell( 
          cell, scratch_data, copy_data, error_indicators, face_integrals); 
      }; 

      auto do_nothing_copier = 
        std::function<void(const WeightedResidualCopyData &)>(); 

// 然后将其全部交给 WorkStream::run() ，以并行计算所有单元的估计器。

      WorkStream::run( 
        DualSolver<dim>::dof_handler.begin_active(), 
        DualSolver<dim>::dof_handler.end(), 
        worker, 
        do_nothing_copier, 
        WeightedResidualScratchData(*DualSolver<dim>::fe, 
                                    *DualSolver<dim>::quadrature, 
                                    *DualSolver<dim>::face_quadrature, 
                                    *this->rhs_function, 
                                    primal_solution, 
                                    dual_weights), 
        WeightedResidualCopyData()); 

// 一旦计算出误差贡献，就把它们加起来。为此，请注意，单元格项已经设置好了，只有边缘项需要收集。因此，在所有单元格和它们的面中循环，确保每个面的贡献都在那里，然后把它们加起来。只需要减去一半的跳跃项，因为另一半将被邻近的单元格拿走。

      unsigned int present_cell = 0; 
      for (const auto &cell : 
           DualSolver<dim>::dof_handler.active_cell_iterators()) 
        { 
          for (const auto &face : cell->face_iterators()) 
            { 
              Assert(face_integrals.find(face) != face_integrals.end(), 
                     ExcInternalError()); 
              error_indicators(present_cell) -= 0.5 * face_integrals[face]; 
            } 
          ++present_cell; 
        } 
      std::cout << "   Estimated error=" 
                << std::accumulate(error_indicators.begin(), 
                                   error_indicators.end(), 
                                   0.) 
                << std::endl; 
    } 
// @sect4{Estimating on a single cell}  

// 接下来我们有一个函数，它被调用来估计单个单元的误差。如果库被配置为使用多线程，该函数可能被多次调用。下面是它的内容。

    template <int dim> 
    void WeightedResidual<dim>::estimate_on_one_cell( 
      const active_cell_iterator & cell, 
      WeightedResidualScratchData &scratch_data, 
      WeightedResidualCopyData &   copy_data, 
      Vector<float> &              error_indicators, 
      FaceIntegrals &              face_integrals) const 
    { 

// 由于WorkStream的原因， estimate_on_one_cell需要一个CopyData对象，即使它没有被使用。下一行将对这个未使用的变量发出警告。

      (void)copy_data; 

// 每个单元的第一个任务是计算这个单元的剩余贡献，并把它们放入  <code>error_indicators</code>  变量中。

      integrate_over_cell(cell, 
                          scratch_data.primal_solution, 
                          scratch_data.dual_weights, 
                          scratch_data.cell_data, 
                          error_indicators); 

// 计算完单元格条款后，转向面条款。为此，在当前单元格的所有面中进行循环，看看是否需要对其进行计算。

      for (const auto face_no : cell->face_indices()) 
        { 

// 首先，如果这个面是边界的一部分，那么就没有什么可做的。然而，为了在汇总单元格的面的贡献时使事情变得简单，我们把这个面输入对误差贡献为零的面的列表中。

          if (cell->face(face_no)->at_boundary()) 
            { 
              face_integrals[cell->face(face_no)] = 0; 
              continue; 
            } 

// 接下来，请注意，由于我们想在每个面上只计算一次跳跃项，尽管我们访问它两次（如果它不在边界），我们必须定义一些规则，由谁负责在一个面上计算。

// 首先，如果相邻的单元格与这个单元格处于同一层次，也就是说，既不进一步细化，也不进一步粗化，那么这个层次中索引较低的单元格就负责计算。换句话说：如果另一个单元的索引更低，那么就跳过这个面的工作。

          if ((cell->neighbor(face_no)->has_children() == false) && 
              (cell->neighbor(face_no)->level() == cell->level()) && 
              (cell->neighbor(face_no)->index() < cell->index())) 
            continue; 

// 同样地，如果这个单元和它的邻居在细化程度上有差异，我们总是从较粗的单元开始工作。因此，如果相邻的单元比现在的单元细化程度低，那么就什么都不做，因为我们在访问粗略的单元时对子面进行整合。

          if (cell->at_boundary(face_no) == false) 
            if (cell->neighbor(face_no)->level() < cell->level()) 
              continue; 

// 现在我们知道，我们在这里负责，所以实际上是在计算面的跳跃项。如果这个面是一个规则的面，即另一边的单元格既不比这个单元格粗也不比这个单元格细，那么就调用一个函数，如果另一边的单元格进一步细化，那么就用另一个函数。请注意，另一边的单元格更粗的情况不可能发生，因为我们上面已经决定，当我们传递到另一个单元格时，我们会处理这种情况。

          if (cell->face(face_no)->has_children() == false) 
            integrate_over_regular_face(cell, 
                                        face_no, 
                                        scratch_data.primal_solution, 
                                        scratch_data.dual_weights, 
                                        scratch_data.face_data, 
                                        face_integrals); 
          else 
            integrate_over_irregular_face(cell, 
                                          face_no, 
                                          scratch_data.primal_solution, 
                                          scratch_data.dual_weights, 
                                          scratch_data.face_data, 
                                          face_integrals); 
        } 
    } 
// @sect4{Computing cell term error contributions}  

// 关于误差贡献的实际计算，首先转向单元条款。

    template <int dim> 
    void WeightedResidual<dim>::integrate_over_cell( 
      const active_cell_iterator &cell, 
      const Vector<double> &      primal_solution, 
      const Vector<double> &      dual_weights, 
      CellData &                  cell_data, 
      Vector<float> &             error_indicators) const 
    { 

// 需要完成的任务是通过观察误差估计公式看起来很自然的事情：首先在正交点得到单元残差的数值解的右手边和拉普拉斯。

      cell_data.fe_values.reinit(cell); 
      cell_data.right_hand_side->value_list( 
        cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values); 
      cell_data.fe_values.get_function_laplacians(primal_solution, 
                                                  cell_data.cell_laplacians); 

// ...然后得到双重权重...

      cell_data.fe_values.get_function_values(dual_weights, 
                                              cell_data.dual_weights); 

// ...最后建立所有正交点的总和，并将其存储在当前单元格中。

      double sum = 0; 
      for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p) 
        sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) * 
                cell_data.dual_weights[p] * cell_data.fe_values.JxW(p)); 
      error_indicators(cell->active_cell_index()) += sum; 
    } 
// @sect4{Computing edge term error contributions -- 1}  

// 另一方面，误差估计的边缘项的计算并不那么简单。首先，我们必须区分有悬挂节点和无悬挂节点的面。因为这是一种简单的情况，我们首先考虑一个面上没有悬挂节点的情况（我们称之为 "常规 "情况）。

    template <int dim> 
    void WeightedResidual<dim>::integrate_over_regular_face( 
      const active_cell_iterator &cell, 
      const unsigned int          face_no, 
      const Vector<double> &      primal_solution, 
      const Vector<double> &      dual_weights, 
      FaceData &                  face_data, 
      FaceIntegrals &             face_integrals) const 
    { 
      const unsigned int n_q_points = 
        face_data.fe_face_values_cell.n_quadrature_points; 

// 第一步是获取本单元上有限元场的正交点的梯度值。为此，初始化 <code>FEFaceValues</code> 对象，对应于面的这一侧，并使用该对象提取梯度。

      face_data.fe_face_values_cell.reinit(cell, face_no); 
      face_data.fe_face_values_cell.get_function_gradients( 
        primal_solution, face_data.cell_grads); 

// 第二步是提取面的另一侧正交点上的有限元解的梯度，即从相邻单元提取。

// 为此，在之前做一个理智的检查：确保邻居确实存在（是的，如果邻居不存在，我们就不应该来这里，但是在复杂的软件中会有bug，所以最好检查一下），如果不是这样，就扔一个错误。

      Assert(cell->neighbor(face_no).state() == IteratorState::valid, 
             ExcInternalError()); 

// 如果我们有了这个，那么我们需要找出相邻单元格的哪个面，也就是说， <code>how-many'th</code> 这个单元格是这个面后面的单元格的相邻面。为此，有一个函数，我们将结果放入一个变量，名称为 <code>neighbor_neighbor</code>  。

      const unsigned int neighbor_neighbor = 
        cell->neighbor_of_neighbor(face_no); 

// 然后定义一个邻近单元的缩写，在该单元上初始化 <code>FEFaceValues</code> 对象，并提取该单元上的梯度。

      const active_cell_iterator neighbor = cell->neighbor(face_no); 
      face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor); 
      face_data.fe_face_values_neighbor.get_function_gradients( 
        primal_solution, face_data.neighbor_grads); 

// 现在我们有了这个单元和邻近单元的梯度，通过将梯度的跳跃与法向量相乘来计算跳跃残差。

      for (unsigned int p = 0; p < n_q_points; ++p) 
        face_data.jump_residual[p] = 
          ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) * 
           face_data.fe_face_values_cell.normal_vector(p)); 

// 接下来得到这个面的双重权重。

      face_data.fe_face_values_cell.get_function_values(dual_weights, 
                                                        face_data.dual_weights); 

// 最后，我们要计算跳跃残差、对偶权重和正交权重的总和，以得到这个面的结果。

      double face_integral = 0; 
      for (unsigned int p = 0; p < n_q_points; ++p) 
        face_integral += 
          (face_data.jump_residual[p] * face_data.dual_weights[p] * 
           face_data.fe_face_values_cell.JxW(p)); 

// 仔细检查该元素是否已经存在，是否已经被写入...

      Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(), 
             ExcInternalError()); 
      Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError()); 

// ...然后在指定的位置存储计算值。注意，存储的值不包含错误表示中出现的因子1/2。原因是，如果我们在三角形的所有面上进行循环，这个项实际上没有这个因子，但只有当我们把它写成所有单元和每个单元的所有面的总和时才会出现；因此我们两次访问同一个面。我们稍后在对每个单元的贡献进行单独求和时，会使用这个因子-1/2来考虑这个问题。

      face_integrals[cell->face(face_no)] = face_integral; 
    } 
// @sect4{Computing edge term error contributions -- 2}  

// 我们仍然缺少有悬挂节点的面的情况。这就是这个函数中所涉及的内容。

    template <int dim> 
    void WeightedResidual<dim>::integrate_over_irregular_face( 
      const active_cell_iterator &cell, 
      const unsigned int          face_no, 
      const Vector<double> &      primal_solution, 
      const Vector<double> &      dual_weights, 
      FaceData &                  face_data, 
      FaceIntegrals &             face_integrals) const 
    { 

// 首先还是两个缩写，以及一些一致性检查，以确定该函数是否只在它应该被调用的面上被调用。

      const unsigned int n_q_points = 
        face_data.fe_face_values_cell.n_quadrature_points; 

      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no); 
      const typename DoFHandler<dim>::cell_iterator neighbor = 
        cell->neighbor(face_no); 
      Assert(neighbor.state() == IteratorState::valid, ExcInternalError()); 
      Assert(neighbor->has_children(), ExcInternalError()); 
      (void)neighbor; 

// 然后找出当前单元格是相邻单元格的哪个邻居。请注意，我们将对这个相邻单元的子女进行操作，但他们的方向与他们的母亲相同，也就是说，邻居的方向是一样的。

      const unsigned int neighbor_neighbor = 
        cell->neighbor_of_neighbor(face_no); 

// 然后简单地对所有的子面做我们在前面的函数中对一个面所做的一切。

      for (unsigned int subface_no = 0; subface_no < face->n_children(); 
           ++subface_no) 
        { 

// 再从一些检查开始：得到一个指向当前子面后面的单元格的迭代器，并检查其面是否是我们正在考虑的子面。如果不是这样，那么要么是上面调用的 <code>neighbor_neighbor</code> 函数存在错误，要么--更糟糕的是--库中的某些函数没有遵守关于单元格、它们的子面的一些基本假设。在任何情况下，即使这个断言不应该被触发，谨慎一点也无妨，而且在优化模式的计算中，这个断言无论如何都会被删除。

          const active_cell_iterator neighbor_child = 
            cell->neighbor_child_on_subface(face_no, subface_no); 
          Assert(neighbor_child->face(neighbor_neighbor) == 
                   cell->face(face_no)->child(subface_no), 
                 ExcInternalError()); 

// 现在开始工作，首先在界面的这一侧再次得到解决方案的梯度。

          face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no); 
          face_data.fe_subface_values_cell.get_function_gradients( 
            primal_solution, face_data.cell_grads); 

// 然后在另一边。

          face_data.fe_face_values_neighbor.reinit(neighbor_child, 
                                                   neighbor_neighbor); 
          face_data.fe_face_values_neighbor.get_function_gradients( 
            primal_solution, face_data.neighbor_grads); 

//最后建立跳跃残差。因为这次我们从另一个单元格中取法向量，所以与其他函数相比，将第一项的符号还原。

          for (unsigned int p = 0; p < n_q_points; ++p) 
            face_data.jump_residual[p] = 
              ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) * 
               face_data.fe_face_values_neighbor.normal_vector(p)); 

// 然后得到双重权重。

          face_data.fe_face_values_neighbor.get_function_values( 
            dual_weights, face_data.dual_weights); 

// 最后，总结这个子面的贡献，并将其设置在全局图中。

          double face_integral = 0; 
          for (unsigned int p = 0; p < n_q_points; ++p) 
            face_integral += 
              (face_data.jump_residual[p] * face_data.dual_weights[p] * 
               face_data.fe_face_values_neighbor.JxW(p)); 
          face_integrals[neighbor_child->face(neighbor_neighbor)] = 
            face_integral; 
        } 

// 一旦所有子面的贡献被计算出来，循环收集所有子面，并将其与母面一起存储，以便以后收集单元格的误差项时简单使用。再次进行安全检查，确保子面的条目已经被计算出来，并且不带有无效的值。

      double sum = 0; 
      for (unsigned int subface_no = 0; subface_no < face->n_children(); 
           ++subface_no) 
        { 
          Assert(face_integrals.find(face->child(subface_no)) != 
                   face_integrals.end(), 
                 ExcInternalError());  
          Assert(face_integrals[face->child(subface_no)] != -1e20, 
                 ExcInternalError()); 

          sum += face_integrals[face->child(subface_no)]; 
        } 

// 最后将该值存储在父脸。

      face_integrals[face] = sum; 
    } 

  } // namespace LaplaceSolver 
// @sect3{A simulation framework}  

// 在前面的例子程序中，我们有两个函数，用来驱动在随后的更细的网格上求解的过程。我们在这里进行了扩展，允许向这些函数传递一些参数，并将这些参数全部放入框架类中。

// 你会注意到这个程序是由许多小部分组成的（评估函数、实现各种细化方法的求解器类、不同的对偶函数、不同的问题和数据描述），这使得程序的扩展相对简单，但也允许通过用一个部分替换另一个部分来解决大量的不同问题。我们通过在下面的框架类中声明一个结构来体现这种灵活性，该结构持有一些参数，可以设置这些参数来测试这个程序的各个部分的组合，可以用简单的方式在各种问题和离散度上进行测试。

  template <int dim> 
  struct Framework 
  { 
  public: 

// 首先，我们声明两个缩写，以便简单使用各自的数据类型。

    using Evaluator     = Evaluation::EvaluationBase<dim>; 
    using EvaluatorList = std::list<Evaluator *>; 

// 然后我们有一个结构，它声明了所有可能被设置的参数。在该结构的默认构造函数中，这些值都被设置为默认值，以便简单使用。

    struct ProblemDescription 
    { 

// 首先允许原始和对偶问题离散化的片状多项式的度数。对于原始问题，它们默认为（双，三）线性分解函数，对于对偶问题，默认为（双，三）二次函数。如果选择了一个不需要解决对偶问题的细化准则，对偶有限元度的值当然会被忽略。

      unsigned int primal_fe_degree; 
      unsigned int dual_fe_degree; 

// 然后有一个描述问题类型的对象，即右手边、领域、边界值等。这里需要的指针默认为Null指针，也就是说，你必须在这个对象的实际实例中设置它，才能使它发挥作用。

      std::unique_ptr<const Data::SetUpBase<dim>> data; 

// 由于我们允许使用不同的细化标准（全局细化、通过凯利误差指标细化，可能有一个权重，以及使用双重估计器），定义一些枚举值，并随后定义一个该类型的变量。它将默认为 <code>dual_weighted_error_estimator</code>  。

      enum RefinementCriterion 
      { 
        dual_weighted_error_estimator, 
        global_refinement, 
        kelly_indicator, 
        weighted_kelly_indicator 
      }; 

      RefinementCriterion refinement_criterion; 

// 接下来是一个描述双重函数的对象。只有在选择双重加权残差细化时才需要这个对象，并且默认为一个空指针。

      std::unique_ptr<const DualFunctional::DualFunctionalBase<dim>> 
        dual_functional; 

// 然后是一个评估对象的列表。其默认值为空，即没有评价对象。

      EvaluatorList evaluator_list; 

// 接下来是一个函数，作为 <code>RefinementWeightedKelly</code> 类的权重。这个指针的默认值是零，但是如果你想使用 <code>weighted_kelly_indicator</code> 的细化标准，你必须把它设置成其他的值。

      std::unique_ptr<const Function<dim>> kelly_weight; 

// 最后，我们有一个变量，表示我们允许的（原始）离散化的最大自由度数。如果超过这个数值，我们将停止解算和间歇性网格细化的过程。其默认值为20,000。

      unsigned int max_degrees_of_freedom; 

// 最后是这个类的默认构造函数。

      ProblemDescription(); 
    }; 

// 驱动程序框架类只有一个方法，它断断续续地调用求解器和网格细化，并在中间做一些其他的小任务。由于它除了给它的参数外不需要其他数据，我们把它变成静态的。

    static void run(const ProblemDescription &descriptor); 
  }; 

// 至于实现，首先是参数对象的构造函数，将所有的值设置为默认值。

  template <int dim> 
  Framework<dim>::ProblemDescription::ProblemDescription() 
    : primal_fe_degree(1) 
    , dual_fe_degree(2) 
    , refinement_criterion(dual_weighted_error_estimator) 
    , max_degrees_of_freedom(20000) 
  {} 

// 然后是驱动整个过程的函数。

  template <int dim> 
  void Framework<dim>::run(const ProblemDescription &descriptor) 
  { 

// 首先从给定的数据对象中创建一个三角图。

    Triangulation<dim> triangulation( 
      Triangulation<dim>::smoothing_on_refinement); 
    descriptor.data->create_coarse_grid(triangulation); 

// 然后是一组有限元和适当的正交公式。

    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree); 
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree); 
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1); 
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1); 

// 接下来，从实现不同细化标准的类中选择一个。

    std::unique_ptr<LaplaceSolver::Base<dim>> solver; 
    switch (descriptor.refinement_criterion) 
      { 
        case ProblemDescription::dual_weighted_error_estimator: 
          { 
            solver = std::make_unique<LaplaceSolver::WeightedResidual<dim>>( 
              triangulation, 
              primal_fe, 
              dual_fe, 
              quadrature, 
              face_quadrature, 
              descriptor.data->get_right_hand_side(), 
              descriptor.data->get_boundary_values(), 
              *descriptor.dual_functional); 
            break; 
          } 

        case ProblemDescription::global_refinement: 
          { 
            solver = std::make_unique<LaplaceSolver::RefinementGlobal<dim>>( 
              triangulation, 
              primal_fe, 
              quadrature, 
              face_quadrature, 
              descriptor.data->get_right_hand_side(), 
              descriptor.data->get_boundary_values()); 
            break; 
          } 

        case ProblemDescription::kelly_indicator: 
          { 
            solver = std::make_unique<LaplaceSolver::RefinementKelly<dim>>( 
              triangulation, 
              primal_fe, 
              quadrature, 
              face_quadrature, 
              descriptor.data->get_right_hand_side(), 
              descriptor.data->get_boundary_values()); 
            break; 
          } 

        case ProblemDescription::weighted_kelly_indicator: 
          { 
            solver = 
              std::make_unique<LaplaceSolver::RefinementWeightedKelly<dim>>( 
                triangulation, 
                primal_fe, 
                quadrature, 
                face_quadrature, 
                descriptor.data->get_right_hand_side(), 
                descriptor.data->get_boundary_values(), 
                *descriptor.kelly_weight); 
            break; 
          } 

        default: 
          AssertThrow(false, ExcInternalError()); 
      } 

// 现在所有对象都到位了，运行主循环。停止的标准在循环的底部实现。

// 在循环中，首先设置新的循环数，然后解决问题，输出它的解，对它应用评估对象，然后决定我们是否要进一步细化网格并在这个网格上再次解决问题，或者跳出循环。

    for (unsigned int step = 0; true; ++step) 
      { 
        std::cout << "Refinement cycle: " << step << std::endl; 

        solver->set_refinement_cycle(step); 
        solver->solve_problem(); 
        solver->output_solution(); 

        std::cout << "   Number of degrees of freedom=" << solver->n_dofs() 
                  << std::endl; 

        for (const auto &evaluator : descriptor.evaluator_list) 
          { 
            evaluator->set_refinement_cycle(step); 
            solver->postprocess(*evaluator); 
          } 

        if (solver->n_dofs() < descriptor.max_degrees_of_freedom) 
          solver->refine_grid(); 
        else 
          break; 
      } 

// 循环运行后清理屏幕。

    std::cout << std::endl; 
  } 

} // namespace Step14 

//  @sect3{The main function}  

// 这里最后是主函数。它通过指定一组用于模拟的参数（多项式度数、评估和对偶函数等）来驱动整个过程，并将它们打包成一个结构传给上面的框架工作类。

int main() 
{ 
  try 
    { 
      using namespace Step14; 

//描述我们要在这里解决的问题，将一个描述符对象传递给做其他工作的函数。

      const unsigned int                 dim = 2; 
      Framework<dim>::ProblemDescription descriptor; 

// 首先设置我们希望使用的细化标准。

      descriptor.refinement_criterion = 
        Framework<dim>::ProblemDescription::dual_weighted_error_estimator; 

// 在这里，我们也可以使用  <code>global_refinement</code>  或  <code>weighted_kelly_indicator</code>  。请注意，所给出的关于对偶有限元、对偶函数等信息只对给定的细化准则选择很重要，否则就会被忽略。

// 然后设置原始问题和对偶问题的多项式程度。我们在这里选择双线性和双二次方的问题。

      descriptor.primal_fe_degree = 1; 
      descriptor.dual_fe_degree   = 2; 

// 然后设置测试案例的描述，即域、边界值和右手边。这些都是预先打包在类中的。我们在这里采用  <code>Exercise_2_3</code>  的描述，但你也可以使用  <code>CurvedRidges@<dim@></code>  。

      descriptor.data = 
        std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>(); 

// 接下来先设置一个二元函数，然后设置一个评价对象的列表。我们默认选择在一个评价点对值进行评价，由评价和二元函数类命名空间中的 <code>PointValueEvaluation</code> 类代表。你也可以设置 <code>PointXDerivativeEvaluation</code> 类来代替评价点上的值的x-derivative。

// 注意，双功能和评价对象应该匹配。然而，你可以给你想要的评价函数，所以你可以在每一步之后让点值和导数都得到评价。 一个这样的附加评价是在每一步中输出网格。

      const Point<dim> evaluation_point(0.75, 0.75); 
      descriptor.dual_functional = 
        std::make_unique<DualFunctional::PointValueEvaluation<dim>>( 
          evaluation_point); 

      Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point); 
      Evaluation::GridOutput<dim>           postprocessor2("grid"); 

      descriptor.evaluator_list.push_back(&postprocessor1); 
      descriptor.evaluator_list.push_back(&postprocessor2); 

// 设置最大的自由度数，在这个自由度数之后，我们希望程序停止进一步细化网格。

      descriptor.max_degrees_of_freedom = 20000; 

// 最后将描述符对象传递给一个函数，用它来运行整个解决方案。

      Framework<dim>::run(descriptor); 
    } 

// 捕获异常以提供有关失败的信息。

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



