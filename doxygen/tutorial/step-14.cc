

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


namespace Step14 
{ 
  using namespace dealii; 



  namespace Evaluation 
  { 
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


    template <int dim> 
    void PointXDerivativeEvaluation<dim>:: 
         operator()(const DoFHandler<dim> &dof_handler, 
               const Vector<double> & solution) const 
    { 


      double point_derivative = 0; 


      QTrapezoid<dim>             vertex_quadrature; 
      FEValues<dim>               fe_values(dof_handler.get_fe(), 
                              vertex_quadrature, 
                              update_gradients | update_quadrature_points); 
      std::vector<Tensor<1, dim>> solution_gradients(vertex_quadrature.size()); 


      unsigned int evaluation_point_hits = 0; 
      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto vertex : cell->vertex_indices()) 
          if (cell->vertex(vertex) == evaluation_point) 
            { 


              fe_values.reinit(cell); 


              fe_values.get_function_gradients(solution, solution_gradients); 


              unsigned int q_point = 0; 
              for (; q_point < solution_gradients.size(); ++q_point) 
                if (fe_values.quadrature_point(q_point) == evaluation_point) 
                  break; 


              Assert(q_point < solution_gradients.size(), ExcInternalError()); 


              point_derivative += solution_gradients[q_point][0]; 
              ++evaluation_point_hits; 


              break; 
            } 


      AssertThrow(evaluation_point_hits > 0, 
                  ExcEvaluationPointNotFound(evaluation_point)); 


      point_derivative /= evaluation_point_hits; 
      std::cout << "   Point x-derivative=" << point_derivative << std::endl; 
    } 




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


  namespace LaplaceSolver 
  { 


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



    template <int dim> 
    Solver<dim>::LinearSystem::LinearSystem(const DoFHandler<dim> &dof_handler) 
    { 
      hanging_node_constraints.clear(); 

      void (*mhnc_p)(const DoFHandler<dim> &, AffineConstraints<double> &) = 
        &DoFTools::make_hanging_node_constraints; 


      Threads::Task<void> side_task = 
        Threads::new_task(mhnc_p, dof_handler, hanging_node_constraints); 

      DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler, dsp); 


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


    template <int dim> 
    void RefinementWeightedKelly<dim>::refine_grid() 
    { 


      Vector<float> estimated_error_per_cell( 
        this->triangulation->n_active_cells()); 
      std::map<types::boundary_id, const Function<dim> *> dummy_function_map; 
      KellyErrorEstimator<dim>::estimate(this->dof_handler, 
                                         *this->face_quadrature, 
                                         dummy_function_map, 
                                         this->solution, 
                                         estimated_error_per_cell); 


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








  namespace Data 
  { 


    template <int dim> 
    struct SetUpBase : public Subscriptor 
    { 
      virtual const Function<dim> &get_boundary_values() const = 0; 

      virtual const Function<dim> &get_right_hand_side() const = 0; 

      virtual void 
      create_coarse_grid(Triangulation<dim> &coarse_grid) const = 0; 
    }; 



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


    template <class Traits, int dim> 
    const typename Traits::BoundaryValues SetUp<Traits, dim>::boundary_values; 
    template <class Traits, int dim> 
    const typename Traits::RightHandSide SetUp<Traits, dim>::right_hand_side; 


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


    template <int dim> 
    struct Exercise_2_3 
    { 


      using BoundaryValues = Functions::ZeroFunction<dim>; 


      class RightHandSide : public Functions::ConstantFunction<dim> 
      { 
      public: 
        RightHandSide() 
          : Functions::ConstantFunction<dim>(1.) 
        {} 
      }; 


      static void create_coarse_grid(Triangulation<dim> &coarse_grid); 
    }; 




    template <> 
    void Exercise_2_3<2>::create_coarse_grid(Triangulation<2> &coarse_grid) 
    { 


      const unsigned int dim = 2; 

      const std::vector<Point<2>> vertices = { 
        {-1.0, -1.0}, {-0.5, -1.0}, {+0.0, -1.0}, {+0.5, -1.0}, {+1.0, -1.0}, // 
        {-1.0, -0.5}, {-0.5, -0.5}, {+0.0, -0.5}, {+0.5, -0.5}, {+1.0, -0.5}, // 
        {-1.0, +0.0}, {-0.5, +0.0}, {+0.5, +0.0}, {+1.0, +0.0},               // 
        {-1.0, +0.5}, {-0.5, +0.5}, {+0.0, +0.5}, {+0.5, +0.5}, {+1.0, +0.5}, // 
        {-1.0, +1.0}, {-0.5, +1.0}, {+0.0, +1.0}, {+0.5, +1.0}, {+1.0, +1.0}}; 


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


      std::vector<CellData<dim>> cells(n_cells, CellData<dim>()); 
      for (unsigned int i = 0; i < n_cells; ++i) 
        { 
          for (unsigned int j = 0; j < cell_vertices[i].size(); ++j) 
            cells[i].vertices[j] = cell_vertices[i][j]; 
          cells[i].material_id = 0; 
        } 


      coarse_grid.create_triangulation(vertices, cells, SubCellData()); 


      coarse_grid.refine_global(1); 
    } 
  } // namespace Data 







  namespace DualFunctional 
  { 


    template <int dim> 
    class DualFunctionalBase : public Subscriptor 
    { 
    public: 
      virtual void assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                Vector<double> &       rhs) const = 0; 
    }; 


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



    template <int dim> 
    void 
    PointValueEvaluation<dim>::assemble_rhs(const DoFHandler<dim> &dof_handler, 
                                            Vector<double> &       rhs) const 
    { 


      rhs.reinit(dof_handler.n_dofs()); 


      for (const auto &cell : dof_handler.active_cell_iterators()) 
        for (const auto vertex : cell->vertex_indices()) 
          if (cell->vertex(vertex).distance(evaluation_point) < 
              cell->diameter() * 1e-8) 
            { 


              rhs(cell->vertex_dof_index(vertex, 0)) = 1; 
              return; 
            } 


      AssertThrow(false, ExcEvaluationPointNotFound(evaluation_point)); 
    } 


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




    template <int dim> 
    void PointXDerivativeEvaluation<dim>::assemble_rhs( 
      const DoFHandler<dim> &dof_handler, 
      Vector<double> &       rhs) const 
    { 


      rhs.reinit(dof_handler.n_dofs()); 


      QGauss<dim>        quadrature(dof_handler.get_fe().degree + 1); 
      FEValues<dim>      fe_values(dof_handler.get_fe(), 
                              quadrature, 
                              update_gradients | update_quadrature_points | 
                                update_JxW_values); 
      const unsigned int n_q_points    = fe_values.n_quadrature_points; 
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell; 


      Vector<double>            cell_rhs(dofs_per_cell); 
      std::vector<unsigned int> local_dof_indices(dofs_per_cell); 


      double total_volume = 0; 


      for (const auto &cell : dof_handler.active_cell_iterators()) 
        if (cell->center().distance(evaluation_point) <= cell->diameter()) 
          { 


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


            cell->get_dof_indices(local_dof_indices); 
            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              rhs(local_dof_indices[i]) += cell_rhs(i); 
          } 


      AssertThrow(total_volume > 0, 
                  ExcEvaluationPointNotFound(evaluation_point)); 


      rhs /= total_volume; 
    } 

  } // namespace DualFunctional 
  namespace LaplaceSolver 
  { 



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


      void solve_primal_problem(); 
      void solve_dual_problem(); 


      using active_cell_iterator = 
        typename DoFHandler<dim>::active_cell_iterator; 



      using FaceIntegrals = 
        typename std::map<typename DoFHandler<dim>::face_iterator, double>; 





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

      struct WeightedResidualCopyData 
      {}; 


      void estimate_error(Vector<float> &error_indicators) const; 

      void estimate_on_one_cell(const active_cell_iterator & cell, 
                                WeightedResidualScratchData &scratch_data, 
                                WeightedResidualCopyData &   copy_data, 
                                Vector<float> &              error_indicators, 
                                FaceIntegrals &face_integrals) const; 


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


    template <int dim> 
    void WeightedResidual<dim>::refine_grid() 
    { 


      Vector<float> error_indicators(this->triangulation->n_active_cells()); 
      estimate_error(error_indicators); 


      for (float &error_indicator : error_indicators) 
        error_indicator = std::fabs(error_indicator); 


      GridRefinement::refine_and_coarsen_fixed_fraction(*this->triangulation, 
                                                        error_indicators, 
                                                        0.8, 
                                                        0.02); 
      this->triangulation->execute_coarsening_and_refinement(); 
    } 


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


      data_out.add_data_vector(PrimalSolver<dim>::solution, "primal_solution"); 
      data_out.add_data_vector(dual_solution, "dual_solution"); 

      data_out.build_patches(); 

      std::ofstream out("solution-" + std::to_string(this->refinement_cycle) + 
                        ".vtu"); 
      data_out.write(out, DataOutBase::vtu); 
    } 


    template <int dim> 
    void 
    WeightedResidual<dim>::estimate_error(Vector<float> &error_indicators) const 
    { 



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


    template <int dim> 
    void WeightedResidual<dim>::estimate_on_one_cell( 
      const active_cell_iterator & cell, 
      WeightedResidualScratchData &scratch_data, 
      WeightedResidualCopyData &   copy_data, 
      Vector<float> &              error_indicators, 
      FaceIntegrals &              face_integrals) const 
    { 


      (void)copy_data; 


      integrate_over_cell(cell, 
                          scratch_data.primal_solution, 
                          scratch_data.dual_weights, 
                          scratch_data.cell_data, 
                          error_indicators); 


      for (const auto face_no : cell->face_indices()) 
        { 


          if (cell->face(face_no)->at_boundary()) 
            { 
              face_integrals[cell->face(face_no)] = 0; 
              continue; 
            } 



          if ((cell->neighbor(face_no)->has_children() == false) && 
              (cell->neighbor(face_no)->level() == cell->level()) && 
              (cell->neighbor(face_no)->index() < cell->index())) 
            continue; 


          if (cell->at_boundary(face_no) == false) 
            if (cell->neighbor(face_no)->level() < cell->level()) 
              continue; 


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


    template <int dim> 
    void WeightedResidual<dim>::integrate_over_cell( 
      const active_cell_iterator &cell, 
      const Vector<double> &      primal_solution, 
      const Vector<double> &      dual_weights, 
      CellData &                  cell_data, 
      Vector<float> &             error_indicators) const 
    { 


      cell_data.fe_values.reinit(cell); 
      cell_data.right_hand_side->value_list( 
        cell_data.fe_values.get_quadrature_points(), cell_data.rhs_values); 
      cell_data.fe_values.get_function_laplacians(primal_solution, 
                                                  cell_data.cell_laplacians); 


      cell_data.fe_values.get_function_values(dual_weights, 
                                              cell_data.dual_weights); 


      double sum = 0; 
      for (unsigned int p = 0; p < cell_data.fe_values.n_quadrature_points; ++p) 
        sum += ((cell_data.rhs_values[p] + cell_data.cell_laplacians[p]) * 
                cell_data.dual_weights[p] * cell_data.fe_values.JxW(p)); 
      error_indicators(cell->active_cell_index()) += sum; 
    } 


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


      face_data.fe_face_values_cell.reinit(cell, face_no); 
      face_data.fe_face_values_cell.get_function_gradients( 
        primal_solution, face_data.cell_grads); 



      Assert(cell->neighbor(face_no).state() == IteratorState::valid, 
             ExcInternalError()); 


      const unsigned int neighbor_neighbor = 
        cell->neighbor_of_neighbor(face_no); 


      const active_cell_iterator neighbor = cell->neighbor(face_no); 
      face_data.fe_face_values_neighbor.reinit(neighbor, neighbor_neighbor); 
      face_data.fe_face_values_neighbor.get_function_gradients( 
        primal_solution, face_data.neighbor_grads); 


      for (unsigned int p = 0; p < n_q_points; ++p) 
        face_data.jump_residual[p] = 
          ((face_data.cell_grads[p] - face_data.neighbor_grads[p]) * 
           face_data.fe_face_values_cell.normal_vector(p)); 


      face_data.fe_face_values_cell.get_function_values(dual_weights, 
                                                        face_data.dual_weights); 


      double face_integral = 0; 
      for (unsigned int p = 0; p < n_q_points; ++p) 
        face_integral += 
          (face_data.jump_residual[p] * face_data.dual_weights[p] * 
           face_data.fe_face_values_cell.JxW(p)); 


      Assert(face_integrals.find(cell->face(face_no)) != face_integrals.end(), 
             ExcInternalError()); 
      Assert(face_integrals[cell->face(face_no)] == -1e20, ExcInternalError()); 


      face_integrals[cell->face(face_no)] = face_integral; 
    } 


    template <int dim> 
    void WeightedResidual<dim>::integrate_over_irregular_face( 
      const active_cell_iterator &cell, 
      const unsigned int          face_no, 
      const Vector<double> &      primal_solution, 
      const Vector<double> &      dual_weights, 
      FaceData &                  face_data, 
      FaceIntegrals &             face_integrals) const 
    { 


      const unsigned int n_q_points = 
        face_data.fe_face_values_cell.n_quadrature_points; 

      const typename DoFHandler<dim>::face_iterator face = cell->face(face_no); 
      const typename DoFHandler<dim>::cell_iterator neighbor = 
        cell->neighbor(face_no); 
      Assert(neighbor.state() == IteratorState::valid, ExcInternalError()); 
      Assert(neighbor->has_children(), ExcInternalError()); 
      (void)neighbor; 


      const unsigned int neighbor_neighbor = 
        cell->neighbor_of_neighbor(face_no); 


      for (unsigned int subface_no = 0; subface_no < face->n_children(); 
           ++subface_no) 
        { 


          const active_cell_iterator neighbor_child = 
            cell->neighbor_child_on_subface(face_no, subface_no); 
          Assert(neighbor_child->face(neighbor_neighbor) == 
                   cell->face(face_no)->child(subface_no), 
                 ExcInternalError()); 


          face_data.fe_subface_values_cell.reinit(cell, face_no, subface_no); 
          face_data.fe_subface_values_cell.get_function_gradients( 
            primal_solution, face_data.cell_grads); 


          face_data.fe_face_values_neighbor.reinit(neighbor_child, 
                                                   neighbor_neighbor); 
          face_data.fe_face_values_neighbor.get_function_gradients( 
            primal_solution, face_data.neighbor_grads); 


          for (unsigned int p = 0; p < n_q_points; ++p) 
            face_data.jump_residual[p] = 
              ((face_data.neighbor_grads[p] - face_data.cell_grads[p]) * 
               face_data.fe_face_values_neighbor.normal_vector(p)); 


          face_data.fe_face_values_neighbor.get_function_values( 
            dual_weights, face_data.dual_weights); 


          double face_integral = 0; 
          for (unsigned int p = 0; p < n_q_points; ++p) 
            face_integral += 
              (face_data.jump_residual[p] * face_data.dual_weights[p] * 
               face_data.fe_face_values_neighbor.JxW(p)); 
          face_integrals[neighbor_child->face(neighbor_neighbor)] = 
            face_integral; 
        } 


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


      face_integrals[face] = sum; 
    } 

  } // namespace LaplaceSolver 



  template <int dim> 
  struct Framework 
  { 
  public: 


    using Evaluator     = Evaluation::EvaluationBase<dim>; 
    using EvaluatorList = std::list<Evaluator *>; 


    struct ProblemDescription 
    { 


      unsigned int primal_fe_degree; 
      unsigned int dual_fe_degree; 


      std::unique_ptr<const Data::SetUpBase<dim>> data; 


      enum RefinementCriterion 
      { 
        dual_weighted_error_estimator, 
        global_refinement, 
        kelly_indicator, 
        weighted_kelly_indicator 
      }; 

      RefinementCriterion refinement_criterion; 


      std::unique_ptr<const DualFunctional::DualFunctionalBase<dim>> 
        dual_functional; 


      EvaluatorList evaluator_list; 


      std::unique_ptr<const Function<dim>> kelly_weight; 


      unsigned int max_degrees_of_freedom; 


      ProblemDescription(); 
    }; 


    static void run(const ProblemDescription &descriptor); 
  }; 


  template <int dim> 
  Framework<dim>::ProblemDescription::ProblemDescription() 
    : primal_fe_degree(1) 
    , dual_fe_degree(2) 
    , refinement_criterion(dual_weighted_error_estimator) 
    , max_degrees_of_freedom(20000) 
  {} 


  template <int dim> 
  void Framework<dim>::run(const ProblemDescription &descriptor) 
  { 


    Triangulation<dim> triangulation( 
      Triangulation<dim>::smoothing_on_refinement); 
    descriptor.data->create_coarse_grid(triangulation); 


    const FE_Q<dim>       primal_fe(descriptor.primal_fe_degree); 
    const FE_Q<dim>       dual_fe(descriptor.dual_fe_degree); 
    const QGauss<dim>     quadrature(descriptor.dual_fe_degree + 1); 
    const QGauss<dim - 1> face_quadrature(descriptor.dual_fe_degree + 1); 


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


    std::cout << std::endl; 
  } 

} // namespace Step14 



int main() 
{ 
  try 
    { 
      using namespace Step14; 


      const unsigned int                 dim = 2; 
      Framework<dim>::ProblemDescription descriptor; 


      descriptor.refinement_criterion = 
        Framework<dim>::ProblemDescription::dual_weighted_error_estimator; 



      descriptor.primal_fe_degree = 1; 
      descriptor.dual_fe_degree   = 2; 


      descriptor.data = 
        std::make_unique<Data::SetUp<Data::Exercise_2_3<dim>, dim>>(); 



      const Point<dim> evaluation_point(0.75, 0.75); 
      descriptor.dual_functional = 
        std::make_unique<DualFunctional::PointValueEvaluation<dim>>( 
          evaluation_point); 

      Evaluation::PointValueEvaluation<dim> postprocessor1(evaluation_point); 
      Evaluation::GridOutput<dim>           postprocessor2("grid"); 

      descriptor.evaluator_list.push_back(&postprocessor1); 
      descriptor.evaluator_list.push_back(&postprocessor2); 


      descriptor.max_degrees_of_freedom = 20000; 


      Framework<dim>::run(descriptor); 
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



