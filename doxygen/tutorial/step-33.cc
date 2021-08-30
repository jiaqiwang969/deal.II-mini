

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
 * Author: David Neckels, Boulder, Colorado, 2007, 2008 
 */ 




#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/function_parser.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/conditional_ostream.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_in.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/mapping_q1.h> 
#include <deal.II/fe/fe_q.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/solution_transfer.h> 



#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_solver.h> 


#include <Sacado.hpp> 


#include <iostream> 
#include <fstream> 
#include <vector> 
#include <memory> 
#include <array> 


namespace Step33 
{ 
  using namespace dealii; 


  template <int dim> 
  struct EulerEquations 
  { 


    static const unsigned int n_components             = dim + 2; 
    static const unsigned int first_momentum_component = 0; 
    static const unsigned int density_component        = dim; 
    static const unsigned int energy_component         = dim + 1; 


    static std::vector<std::string> component_names() 
    { 
      std::vector<std::string> names(dim, "momentum"); 
      names.emplace_back("density"); 
      names.emplace_back("energy_density"); 

      return names; 
    } 

    static std::vector<DataComponentInterpretation::DataComponentInterpretation> 
    component_interpretation() 
    { 
      std::vector<DataComponentInterpretation::DataComponentInterpretation> 
        data_component_interpretation( 
          dim, DataComponentInterpretation::component_is_part_of_vector); 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 
      data_component_interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 

      return data_component_interpretation; 
    } 


    static const double gas_gamma; 


    template <typename InputVector> 
    static typename InputVector::value_type 
    compute_kinetic_energy(const InputVector &W) 
    { 
      typename InputVector::value_type kinetic_energy = 0; 
      for (unsigned int d = 0; d < dim; ++d) 
        kinetic_energy += 
          W[first_momentum_component + d] * W[first_momentum_component + d]; 
      kinetic_energy *= 1. / (2 * W[density_component]); 

      return kinetic_energy; 
    } 

    template <typename InputVector> 
    static typename InputVector::value_type 
    compute_pressure(const InputVector &W) 
    { 
      return ((gas_gamma - 1.0) * 
              (W[energy_component] - compute_kinetic_energy(W))); 
    } 



    template <typename InputVector> 
    static void compute_flux_matrix(const InputVector &W, 
                                    ndarray<typename InputVector::value_type, 
                                            EulerEquations<dim>::n_components, 
                                            dim> &     flux) 
    { 


      const typename InputVector::value_type pressure = compute_pressure(W); 

      for (unsigned int d = 0; d < dim; ++d) 
        { 
          for (unsigned int e = 0; e < dim; ++e) 
            flux[first_momentum_component + d][e] = 
              W[first_momentum_component + d] * 
              W[first_momentum_component + e] / W[density_component]; 

          flux[first_momentum_component + d][d] += pressure; 
        } 


      for (unsigned int d = 0; d < dim; ++d) 
        flux[density_component][d] = W[first_momentum_component + d]; 

      for (unsigned int d = 0; d < dim; ++d) 
        flux[energy_component][d] = W[first_momentum_component + d] / 
                                    W[density_component] * 
                                    (W[energy_component] + pressure); 
    } 


    template <typename InputVector> 
    static void numerical_normal_flux( 
      const Tensor<1, dim> &                                      normal, 
      const InputVector &                                         Wplus, 
      const InputVector &                                         Wminus, 
      const double                                                alpha, 
      std::array<typename InputVector::value_type, n_components> &normal_flux) 
    { 
      ndarray<typename InputVector::value_type, 
              EulerEquations<dim>::n_components, 
              dim> 
        iflux, oflux; 

      compute_flux_matrix(Wplus, iflux); 
      compute_flux_matrix(Wminus, oflux); 

      for (unsigned int di = 0; di < n_components; ++di) 
        { 
          normal_flux[di] = 0; 
          for (unsigned int d = 0; d < dim; ++d) 
            normal_flux[di] += 0.5 * (iflux[di][d] + oflux[di][d]) * normal[d]; 

          normal_flux[di] += 0.5 * alpha * (Wplus[di] - Wminus[di]); 
        } 
    } 


    template <typename InputVector> 
    static void compute_forcing_vector( 
      const InputVector &                                         W, 
      std::array<typename InputVector::value_type, n_components> &forcing) 
    { 
      const double gravity = -1.0; 

      for (unsigned int c = 0; c < n_components; ++c) 
        switch (c) 
          { 
            case first_momentum_component + dim - 1: 
              forcing[c] = gravity * W[density_component]; 
              break; 
            case energy_component: 
              forcing[c] = gravity * W[first_momentum_component + dim - 1]; 
              break; 
            default: 
              forcing[c] = 0; 
          } 
    } 


    enum BoundaryKind 
    { 
      inflow_boundary, 
      outflow_boundary, 
      no_penetration_boundary, 
      pressure_boundary 
    }; 




    template <typename DataVector> 
    static void 
    compute_Wminus(const std::array<BoundaryKind, n_components> &boundary_kind, 
                   const Tensor<1, dim> &                        normal_vector, 
                   const DataVector &                            Wplus, 
                   const Vector<double> &boundary_values, 
                   const DataVector &    Wminus) 
    { 
      for (unsigned int c = 0; c < n_components; c++) 
        switch (boundary_kind[c]) 
          { 
            case inflow_boundary: 
              { 
                Wminus[c] = boundary_values(c); 
                break; 
              } 

            case outflow_boundary: 
              { 
                Wminus[c] = Wplus[c]; 
                break; 
              } 


            case pressure_boundary: 
              { 
                const typename DataVector::value_type density = 
                  (boundary_kind[density_component] == inflow_boundary ? 
                     boundary_values(density_component) : 
                     Wplus[density_component]); 

                typename DataVector::value_type kinetic_energy = 0; 
                for (unsigned int d = 0; d < dim; ++d) 
                  if (boundary_kind[d] == inflow_boundary) 
                    kinetic_energy += boundary_values(d) * boundary_values(d); 
                  else 
                    kinetic_energy += Wplus[d] * Wplus[d]; 
                kinetic_energy *= 1. / 2. / density; 

                Wminus[c] = 
                  boundary_values(c) / (gas_gamma - 1.0) + kinetic_energy; 

                break; 
              } 

            case no_penetration_boundary: 
              { 


                typename DataVector::value_type vdotn = 0; 
                for (unsigned int d = 0; d < dim; d++) 
                  { 
                    vdotn += Wplus[d] * normal_vector[d]; 
                  } 

                Wminus[c] = Wplus[c] - 2.0 * vdotn * normal_vector[c]; 
                break; 
              } 

            default: 
              Assert(false, ExcNotImplemented()); 
          } 
    } 



    static void 
    compute_refinement_indicators(const DoFHandler<dim> &dof_handler, 
                                  const Mapping<dim> &   mapping, 
                                  const Vector<double> & solution, 
                                  Vector<double> &       refinement_indicators) 
    { 
      const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 
      std::vector<unsigned int> dofs(dofs_per_cell); 

      const QMidpoint<dim> quadrature_formula; 
      const UpdateFlags    update_flags = update_gradients; 
      FEValues<dim>        fe_v(mapping, 
                         dof_handler.get_fe(), 
                         quadrature_formula, 
                         update_flags); 

      std::vector<std::vector<Tensor<1, dim>>> dU( 
        1, std::vector<Tensor<1, dim>>(n_components)); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          const unsigned int cell_no = cell->active_cell_index(); 
          fe_v.reinit(cell); 
          fe_v.get_function_gradients(solution, dU); 

          refinement_indicators(cell_no) = std::log( 
            1 + std::sqrt(dU[0][density_component] * dU[0][density_component])); 
        } 
    } 






    class Postprocessor : public DataPostprocessor<dim> 
    { 
    public: 
      Postprocessor(const bool do_schlieren_plot); 

      virtual void evaluate_vector_field( 
        const DataPostprocessorInputs::Vector<dim> &inputs, 
        std::vector<Vector<double>> &computed_quantities) const override; 

      virtual std::vector<std::string> get_names() const override; 

      virtual std::vector< 
        DataComponentInterpretation::DataComponentInterpretation> 
      get_data_component_interpretation() const override; 

      virtual UpdateFlags get_needed_update_flags() const override; 

    private: 
      const bool do_schlieren_plot; 
    }; 
  }; 

  template <int dim> 
  const double EulerEquations<dim>::gas_gamma = 1.4; 

  template <int dim> 
  EulerEquations<dim>::Postprocessor::Postprocessor( 
    const bool do_schlieren_plot) 
    : do_schlieren_plot(do_schlieren_plot) 
  {} 


  template <int dim> 
  void EulerEquations<dim>::Postprocessor::evaluate_vector_field( 
    const DataPostprocessorInputs::Vector<dim> &inputs, 
    std::vector<Vector<double>> &               computed_quantities) const 
  { 


    const unsigned int n_quadrature_points = inputs.solution_values.size(); 

    if (do_schlieren_plot == true) 
      Assert(inputs.solution_gradients.size() == n_quadrature_points, 
             ExcInternalError()); 

    Assert(computed_quantities.size() == n_quadrature_points, 
           ExcInternalError()); 

    Assert(inputs.solution_values[0].size() == n_components, 
           ExcInternalError()); 

    if (do_schlieren_plot == true) 
      { 
        Assert(computed_quantities[0].size() == dim + 2, ExcInternalError()); 
      } 
    else 
      { 
        Assert(computed_quantities[0].size() == dim + 1, ExcInternalError()); 
      } 


    for (unsigned int q = 0; q < n_quadrature_points; ++q) 
      { 
        const double density = inputs.solution_values[q](density_component); 

        for (unsigned int d = 0; d < dim; ++d) 
          computed_quantities[q](d) = 
            inputs.solution_values[q](first_momentum_component + d) / density; 

        computed_quantities[q](dim) = 
          compute_pressure(inputs.solution_values[q]); 

        if (do_schlieren_plot == true) 
          computed_quantities[q](dim + 1) = 
            inputs.solution_gradients[q][density_component] * 
            inputs.solution_gradients[q][density_component]; 
      } 
  } 

  template <int dim> 
  std::vector<std::string> EulerEquations<dim>::Postprocessor::get_names() const 
  { 
    std::vector<std::string> names; 
    for (unsigned int d = 0; d < dim; ++d) 
      names.emplace_back("velocity"); 
    names.emplace_back("pressure"); 

    if (do_schlieren_plot == true) 
      names.emplace_back("schlieren_plot"); 

    return names; 
  } 

  template <int dim> 
  std::vector<DataComponentInterpretation::DataComponentInterpretation> 
  EulerEquations<dim>::Postprocessor::get_data_component_interpretation() const 
  { 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      interpretation(dim, 
                     DataComponentInterpretation::component_is_part_of_vector); 

    interpretation.push_back(DataComponentInterpretation::component_is_scalar); 

    if (do_schlieren_plot == true) 
      interpretation.push_back( 
        DataComponentInterpretation::component_is_scalar); 

    return interpretation; 
  } 

  template <int dim> 
  UpdateFlags 
  EulerEquations<dim>::Postprocessor::get_needed_update_flags() const 
  { 
    if (do_schlieren_plot == true) 
      return update_values | update_gradients; 
    else 
      return update_values; 
  } 





  namespace Parameters 
  { 







    struct Solver 
    { 
      enum SolverType 
      { 
        gmres, 
        direct 
      }; 
      SolverType solver; 

      enum OutputType 
      { 
        quiet, 
        verbose 
      }; 
      OutputType output; 

      double linear_residual; 
      int    max_iterations; 

      double ilut_fill; 
      double ilut_atol; 
      double ilut_rtol; 
      double ilut_drop; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Solver::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("linear solver"); 
      { 
        prm.declare_entry( 
          "output", 
          "quiet", 
          Patterns::Selection("quiet|verbose"), 
          "State whether output from solver runs should be printed. " 
          "Choices are <quiet|verbose>."); 
        prm.declare_entry("method", 
                          "gmres", 
                          Patterns::Selection("gmres|direct"), 
                          "The kind of solver for the linear system. " 
                          "Choices are <gmres|direct>."); 
        prm.declare_entry("residual", 
                          "1e-10", 
                          Patterns::Double(), 
                          "Linear solver residual"); 
        prm.declare_entry("max iters", 
                          "300", 
                          Patterns::Integer(), 
                          "Maximum solver iterations"); 
        prm.declare_entry("ilut fill", 
                          "2", 
                          Patterns::Double(), 
                          "Ilut preconditioner fill"); 
        prm.declare_entry("ilut absolute tolerance", 
                          "1e-9", 
                          Patterns::Double(), 
                          "Ilut preconditioner tolerance"); 
        prm.declare_entry("ilut relative tolerance", 
                          "1.1", 
                          Patterns::Double(), 
                          "Ilut relative tolerance"); 
        prm.declare_entry("ilut drop tolerance", 
                          "1e-10", 
                          Patterns::Double(), 
                          "Ilut drop tolerance"); 
      } 
      prm.leave_subsection(); 
    } 

    void Solver::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("linear solver"); 
      { 
        const std::string op = prm.get("output"); 
        if (op == "verbose") 
          output = verbose; 
        if (op == "quiet") 
          output = quiet; 

        const std::string sv = prm.get("method"); 
        if (sv == "direct") 
          solver = direct; 
        else if (sv == "gmres") 
          solver = gmres; 

        linear_residual = prm.get_double("residual"); 
        max_iterations  = prm.get_integer("max iters"); 
        ilut_fill       = prm.get_double("ilut fill"); 
        ilut_atol       = prm.get_double("ilut absolute tolerance"); 
        ilut_rtol       = prm.get_double("ilut relative tolerance"); 
        ilut_drop       = prm.get_double("ilut drop tolerance"); 
      } 
      prm.leave_subsection(); 
    } 



    struct Refinement 
    { 
      bool   do_refine; 
      double shock_val; 
      double shock_levels; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Refinement::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("refinement"); 
      { 
        prm.declare_entry("refinement", 
                          "true", 
                          Patterns::Bool(), 
                          "Whether to perform mesh refinement or not"); 
        prm.declare_entry("refinement fraction", 
                          "0.1", 
                          Patterns::Double(), 
                          "Fraction of high refinement"); 
        prm.declare_entry("unrefinement fraction", 
                          "0.1", 
                          Patterns::Double(), 
                          "Fraction of low unrefinement"); 
        prm.declare_entry("max elements", 
                          "1000000", 
                          Patterns::Double(), 
                          "maximum number of elements"); 
        prm.declare_entry("shock value", 
                          "4.0", 
                          Patterns::Double(), 
                          "value for shock indicator"); 
        prm.declare_entry("shock levels", 
                          "3.0", 
                          Patterns::Double(), 
                          "number of shock refinement levels"); 
      } 
      prm.leave_subsection(); 
    } 

    void Refinement::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("refinement"); 
      { 
        do_refine    = prm.get_bool("refinement"); 
        shock_val    = prm.get_double("shock value"); 
        shock_levels = prm.get_double("shock levels"); 
      } 
      prm.leave_subsection(); 
    } 



    struct Flux 
    { 
      enum StabilizationKind 
      { 
        constant, 
        mesh_dependent 
      }; 
      StabilizationKind stabilization_kind; 

      double stabilization_value; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Flux::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("flux"); 
      { 
        prm.declare_entry( 
          "stab", 
          "mesh", 
          Patterns::Selection("constant|mesh"), 
          "Whether to use a constant stabilization parameter or " 
          "a mesh-dependent one"); 
        prm.declare_entry("stab value", 
                          "1", 
                          Patterns::Double(), 
                          "alpha stabilization"); 
      } 
      prm.leave_subsection(); 
    } 

    void Flux::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("flux"); 
      { 
        const std::string stab = prm.get("stab"); 
        if (stab == "constant") 
          stabilization_kind = constant; 
        else if (stab == "mesh") 
          stabilization_kind = mesh_dependent; 
        else 
          AssertThrow(false, ExcNotImplemented()); 

        stabilization_value = prm.get_double("stab value"); 
      } 
      prm.leave_subsection(); 
    } 



    struct Output 
    { 
      bool   schlieren_plot; 
      double output_step; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    void Output::declare_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("output"); 
      { 
        prm.declare_entry("schlieren plot", 
                          "true", 
                          Patterns::Bool(), 
                          "Whether or not to produce schlieren plots"); 
        prm.declare_entry("step", 
                          "-1", 
                          Patterns::Double(), 
                          "Output once per this period"); 
      } 
      prm.leave_subsection(); 
    } 

    void Output::parse_parameters(ParameterHandler &prm) 
    { 
      prm.enter_subsection("output"); 
      { 
        schlieren_plot = prm.get_bool("schlieren plot"); 
        output_step    = prm.get_double("step"); 
      } 
      prm.leave_subsection(); 
    } 








    template <int dim> 
    struct AllParameters : public Solver, 
                           public Refinement, 
                           public Flux, 
                           public Output 
    { 
      static const unsigned int max_n_boundaries = 10; 

      struct BoundaryConditions 
      { 
        std::array<typename EulerEquations<dim>::BoundaryKind, 
                   EulerEquations<dim>::n_components> 
          kind; 

        FunctionParser<dim> values; 

        BoundaryConditions(); 
      }; 

      AllParameters(); 

      double diffusion_power; 

      double time_step, final_time; 
      double theta; 
      bool   is_stationary; 

      std::string mesh_filename; 

 
      BoundaryConditions  boundary_conditions[max_n_boundaries]; 

      static void declare_parameters(ParameterHandler &prm); 
      void        parse_parameters(ParameterHandler &prm); 
    }; 

    template <int dim> 
    AllParameters<dim>::BoundaryConditions::BoundaryConditions() 
      : values(EulerEquations<dim>::n_components) 
    { 
      std::fill(kind.begin(), 
                kind.end(), 
                EulerEquations<dim>::no_penetration_boundary); 
    } 

    template <int dim> 
    AllParameters<dim>::AllParameters() 
      : diffusion_power(0.) 
      , time_step(1.) 
      , final_time(1.) 
      , theta(.5) 
      , is_stationary(true) 
      , initial_conditions(EulerEquations<dim>::n_components) 
    {} 

    template <int dim> 
    void AllParameters<dim>::declare_parameters(ParameterHandler &prm) 
    { 
      prm.declare_entry("mesh", 
                        "grid.inp", 
                        Patterns::Anything(), 
                        "input file name"); 

      prm.declare_entry("diffusion power", 
                        "2.0", 
                        Patterns::Double(), 
                        "power of mesh size for diffusion"); 

      prm.enter_subsection("time stepping"); 
      { 
        prm.declare_entry("time step", 
                          "0.1", 
                          Patterns::Double(0), 
                          "simulation time step"); 
        prm.declare_entry("final time", 
                          "10.0", 
                          Patterns::Double(0), 
                          "simulation end time"); 
        prm.declare_entry("theta scheme value", 
                          "0.5", 
                          Patterns::Double(0, 1), 
                          "value for theta that interpolated between explicit " 
                          "Euler (theta=0), Crank-Nicolson (theta=0.5), and " 
                          "implicit Euler (theta=1)."); 
      } 
      prm.leave_subsection(); 

      for (unsigned int b = 0; b < max_n_boundaries; ++b) 
        { 
          prm.enter_subsection("boundary_" + Utilities::int_to_string(b)); 
          { 
            prm.declare_entry("no penetration", 
                              "false", 
                              Patterns::Bool(), 
                              "whether the named boundary allows gas to " 
                              "penetrate or is a rigid wall"); 

            for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
                 ++di) 
              { 
                prm.declare_entry("w_" + Utilities::int_to_string(di), 
                                  "outflow", 
                                  Patterns::Selection( 
                                    "inflow|outflow|pressure"), 
                                  "<inflow|outflow|pressure>"); 

                prm.declare_entry("w_" + Utilities::int_to_string(di) + 
                                    " value", 
                                  "0.0", 
                                  Patterns::Anything(), 
                                  "expression in x,y,z"); 
              } 
          } 
          prm.leave_subsection(); 
        } 

      prm.enter_subsection("initial condition"); 
      { 
        for (unsigned int di = 0; di < EulerEquations<dim>::n_components; ++di) 
          prm.declare_entry("w_" + Utilities::int_to_string(di) + " value", 
                            "0.0", 
                            Patterns::Anything(), 
                            "expression in x,y,z"); 
      } 
      prm.leave_subsection(); 

      Parameters::Solver::declare_parameters(prm); 
      Parameters::Refinement::declare_parameters(prm); 
      Parameters::Flux::declare_parameters(prm); 
      Parameters::Output::declare_parameters(prm); 
    } 

    template <int dim> 
    void AllParameters<dim>::parse_parameters(ParameterHandler &prm) 
    { 
      mesh_filename   = prm.get("mesh"); 
      diffusion_power = prm.get_double("diffusion power"); 

      prm.enter_subsection("time stepping"); 
      { 
        time_step = prm.get_double("time step"); 
        if (time_step == 0) 
          { 
            is_stationary = true; 
            time_step     = 1.0; 
            final_time    = 1.0; 
          } 
        else 
          is_stationary = false; 

        final_time = prm.get_double("final time"); 
        theta      = prm.get_double("theta scheme value"); 
      } 
      prm.leave_subsection(); 

      for (unsigned int boundary_id = 0; boundary_id < max_n_boundaries; 
           ++boundary_id) 
        { 
          prm.enter_subsection("boundary_" + 
                               Utilities::int_to_string(boundary_id)); 
          { 
            std::vector<std::string> expressions( 
              EulerEquations<dim>::n_components, "0.0"); 

            const bool no_penetration = prm.get_bool("no penetration"); 

            for (unsigned int di = 0; di < EulerEquations<dim>::n_components; 
                 ++di) 
              { 
                const std::string boundary_type = 
                  prm.get("w_" + Utilities::int_to_string(di)); 

                if ((di < dim) && (no_penetration == true)) 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::no_penetration_boundary; 
                else if (boundary_type == "inflow") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::inflow_boundary; 
                else if (boundary_type == "pressure") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::pressure_boundary; 
                else if (boundary_type == "outflow") 
                  boundary_conditions[boundary_id].kind[di] = 
                    EulerEquations<dim>::outflow_boundary; 
                else 
                  AssertThrow(false, ExcNotImplemented()); 

                expressions[di] = 
                  prm.get("w_" + Utilities::int_to_string(di) + " value"); 
              } 

            boundary_conditions[boundary_id].values.initialize( 
              FunctionParser<dim>::default_variable_names(), 
              expressions, 
              std::map<std::string, double>()); 
          } 
          prm.leave_subsection(); 
        } 

      prm.enter_subsection("initial condition"); 
      { 
        std::vector<std::string> expressions(EulerEquations<dim>::n_components, 
                                             "0.0"); 
        for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++) 
          expressions[di] = 
            prm.get("w_" + Utilities::int_to_string(di) + " value"); 
        initial_conditions.initialize( 
          FunctionParser<dim>::default_variable_names(), 
          expressions, 
          std::map<std::string, double>()); 
      } 
      prm.leave_subsection(); 

      Parameters::Solver::parse_parameters(prm); 
      Parameters::Refinement::parse_parameters(prm); 
      Parameters::Flux::parse_parameters(prm); 
      Parameters::Output::parse_parameters(prm); 
    } 
  } // namespace Parameters 



  template <int dim> 
  class ConservationLaw 
  { 
  public: 
    ConservationLaw(const char *input_filename); 
    void run(); 

  private: 
    void setup_system(); 

    void assemble_system(); 
    void assemble_cell_term(const FEValues<dim> &                       fe_v, 
                            const std::vector<types::global_dof_index> &dofs); 
    void assemble_face_term( 
      const unsigned int                          face_no, 
      const FEFaceValuesBase<dim> &               fe_v, 
      const FEFaceValuesBase<dim> &               fe_v_neighbor, 
      const std::vector<types::global_dof_index> &dofs, 
      const std::vector<types::global_dof_index> &dofs_neighbor, 
      const bool                                  external_face, 
      const unsigned int                          boundary_id, 
      const double                                face_diameter); 

    std::pair<unsigned int, double> solve(Vector<double> &solution); 

    void compute_refinement_indicators(Vector<double> &indicator) const; 
    void refine_grid(const Vector<double> &indicator); 

    void output_results() const; 


    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 

    const FESystem<dim> fe; 
    DoFHandler<dim>     dof_handler; 

    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 


    Vector<double> old_solution; 
    Vector<double> current_solution; 
    Vector<double> predictor; 

    Vector<double> right_hand_side; 


    TrilinosWrappers::SparseMatrix system_matrix; 

    Parameters::AllParameters<dim> parameters; 
    ConditionalOStream             verbose_cout; 
  }; 


  template <int dim> 
  ConservationLaw<dim>::ConservationLaw(const char *input_filename) 
    : mapping() 
    , fe(FE_Q<dim>(1), EulerEquations<dim>::n_components) 
    , dof_handler(triangulation) 
    , quadrature(fe.degree + 1) 
    , face_quadrature(fe.degree + 1) 
    , verbose_cout(std::cout, false) 
  { 
    ParameterHandler prm; 
    Parameters::AllParameters<dim>::declare_parameters(prm); 

    prm.parse_input(input_filename); 
    parameters.parse_parameters(prm); 

    verbose_cout.set_condition(parameters.output == 
                               Parameters::Solver::verbose); 
  } 



  template <int dim> 
  void ConservationLaw<dim>::setup_system() 
  { 
    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
    DoFTools::make_sparsity_pattern(dof_handler, dsp); 

    system_matrix.reinit(dsp); 
  } 




  template <int dim> 
  void ConservationLaw<dim>::assemble_system() 
  { 
    const unsigned int dofs_per_cell = dof_handler.get_fe().n_dofs_per_cell(); 

    std::vector<types::global_dof_index> dof_indices(dofs_per_cell); 
    std::vector<types::global_dof_index> dof_indices_neighbor(dofs_per_cell); 

    const UpdateFlags update_flags = update_values | update_gradients | 
                                     update_quadrature_points | 
                                     update_JxW_values, 
                      face_update_flags = 
                        update_values | update_quadrature_points | 
                        update_JxW_values | update_normal_vectors, 
                      neighbor_face_update_flags = update_values; 

    FEValues<dim>        fe_v(mapping, fe, quadrature, update_flags); 
    FEFaceValues<dim>    fe_v_face(mapping, 
                                fe, 
                                face_quadrature, 
                                face_update_flags); 
    FESubfaceValues<dim> fe_v_subface(mapping, 
                                      fe, 
                                      face_quadrature, 
                                      face_update_flags); 
    FEFaceValues<dim>    fe_v_face_neighbor(mapping, 
                                         fe, 
                                         face_quadrature, 
                                         neighbor_face_update_flags); 
    FESubfaceValues<dim> fe_v_subface_neighbor(mapping, 
                                               fe, 
                                               face_quadrature, 
                                               neighbor_face_update_flags); 


    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 
        cell->get_dof_indices(dof_indices); 

        assemble_cell_term(fe_v, dof_indices); 


        for (const auto face_no : cell->face_indices()) 
          if (cell->at_boundary(face_no)) 
            { 
              fe_v_face.reinit(cell, face_no); 
              assemble_face_term(face_no, 
                                 fe_v_face, 
                                 fe_v_face, 
                                 dof_indices, 
                                 std::vector<types::global_dof_index>(), 
                                 true, 
                                 cell->face(face_no)->boundary_id(), 
                                 cell->face(face_no)->diameter()); 
            } 






          else 
            { 
              if (cell->neighbor(face_no)->has_children()) 
                { 
                  const unsigned int neighbor2 = 
                    cell->neighbor_of_neighbor(face_no); 

                  for (unsigned int subface_no = 0; 
                       subface_no < cell->face(face_no)->n_children(); 
                       ++subface_no) 
                    { 
                      const typename DoFHandler<dim>::active_cell_iterator 
                        neighbor_child = 
                          cell->neighbor_child_on_subface(face_no, subface_no); 

                      Assert(neighbor_child->face(neighbor2) == 
                               cell->face(face_no)->child(subface_no), 
                             ExcInternalError()); 
                      Assert(neighbor_child->is_active(), ExcInternalError()); 

                      fe_v_subface.reinit(cell, face_no, subface_no); 
                      fe_v_face_neighbor.reinit(neighbor_child, neighbor2); 

                      neighbor_child->get_dof_indices(dof_indices_neighbor); 

                      assemble_face_term( 
                        face_no, 
                        fe_v_subface, 
                        fe_v_face_neighbor, 
                        dof_indices, 
                        dof_indices_neighbor, 
                        false, 
                        numbers::invalid_unsigned_int, 
                        neighbor_child->face(neighbor2)->diameter()); 
                    } 
                } 


              else if (cell->neighbor(face_no)->level() != cell->level()) 
                { 
                  const typename DoFHandler<dim>::cell_iterator neighbor = 
                    cell->neighbor(face_no); 
                  Assert(neighbor->level() == cell->level() - 1, 
                         ExcInternalError()); 

                  neighbor->get_dof_indices(dof_indices_neighbor); 

                  const std::pair<unsigned int, unsigned int> faceno_subfaceno = 
                    cell->neighbor_of_coarser_neighbor(face_no); 
                  const unsigned int neighbor_face_no = faceno_subfaceno.first, 
                                     neighbor_subface_no = 
                                       faceno_subfaceno.second; 

                  Assert(neighbor->neighbor_child_on_subface( 
                           neighbor_face_no, neighbor_subface_no) == cell, 
                         ExcInternalError()); 

                  fe_v_face.reinit(cell, face_no); 
                  fe_v_subface_neighbor.reinit(neighbor, 
                                               neighbor_face_no, 
                                               neighbor_subface_no); 

                  assemble_face_term(face_no, 
                                     fe_v_face, 
                                     fe_v_subface_neighbor, 
                                     dof_indices, 
                                     dof_indices_neighbor, 
                                     false, 
                                     numbers::invalid_unsigned_int, 
                                     cell->face(face_no)->diameter()); 
                } 
            } 
      } 
  } 






  template <int dim> 
  void ConservationLaw<dim>::assemble_cell_term( 
    const FEValues<dim> &                       fe_v, 
    const std::vector<types::global_dof_index> &dof_indices) 
  { 
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 
    const unsigned int n_q_points    = fe_v.n_quadrature_points; 

    Table<2, Sacado::Fad::DFad<double>> W(n_q_points, 
                                          EulerEquations<dim>::n_components); 

    Table<2, double> W_old(n_q_points, EulerEquations<dim>::n_components); 

    Table<3, Sacado::Fad::DFad<double>> grad_W( 
      n_q_points, EulerEquations<dim>::n_components, dim); 

    Table<3, double> grad_W_old(n_q_points, 
                                EulerEquations<dim>::n_components, 
                                dim); 

    std::vector<double> residual_derivatives(dofs_per_cell); 


    std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
      dofs_per_cell); 
    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
      independent_local_dof_values[i] = current_solution(dof_indices[i]); 



    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
      independent_local_dof_values[i].diff(i, dofs_per_cell); 



    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int c = 0; c < EulerEquations<dim>::n_components; ++c) 
        { 
          W[q][c]     = 0; 
          W_old[q][c] = 0; 
          for (unsigned int d = 0; d < dim; ++d) 
            { 
              grad_W[q][c][d]     = 0; 
              grad_W_old[q][c][d] = 0; 
            } 
        } 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          const unsigned int c = 
            fe_v.get_fe().system_to_component_index(i).first; 

          W[q][c] += independent_local_dof_values[i] * 
                     fe_v.shape_value_component(i, q, c); 
          W_old[q][c] += 
            old_solution(dof_indices[i]) * fe_v.shape_value_component(i, q, c); 

          for (unsigned int d = 0; d < dim; d++) 
            { 
              grad_W[q][c][d] += independent_local_dof_values[i] * 
                                 fe_v.shape_grad_component(i, q, c)[d]; 
              grad_W_old[q][c][d] += old_solution(dof_indices[i]) * 
                                     fe_v.shape_grad_component(i, q, c)[d]; 
            } 
        } 


    std::vector<ndarray<Sacado::Fad::DFad<double>, 
                        EulerEquations<dim>::n_components, 
                        dim>> 
      flux(n_q_points); 

    std::vector<ndarray<double, EulerEquations<dim>::n_components, dim>> 
      flux_old(n_q_points); 

 
      std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
      forcing(n_q_points); 

    std::vector<std::array<double, EulerEquations<dim>::n_components>> 
      forcing_old(n_q_points); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        EulerEquations<dim>::compute_flux_matrix(W_old[q], flux_old[q]); 
        EulerEquations<dim>::compute_forcing_vector(W_old[q], forcing_old[q]); 
        EulerEquations<dim>::compute_flux_matrix(W[q], flux[q]); 
        EulerEquations<dim>::compute_forcing_vector(W[q], forcing[q]); 
      } 



    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
      { 
        Sacado::Fad::DFad<double> R_i = 0; 

        const unsigned int component_i = 
          fe_v.get_fe().system_to_component_index(i).first; 


        for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
          { 
            if (parameters.is_stationary == false) 
              R_i += 1.0 / parameters.time_step * 
                     (W[point][component_i] - W_old[point][component_i]) * 
                     fe_v.shape_value_component(i, point, component_i) * 
                     fe_v.JxW(point); 

            for (unsigned int d = 0; d < dim; d++) 
              R_i -= 
                (parameters.theta * flux[point][component_i][d] + 
                 (1.0 - parameters.theta) * flux_old[point][component_i][d]) * 
                fe_v.shape_grad_component(i, point, component_i)[d] * 
                fe_v.JxW(point); 

            for (unsigned int d = 0; d < dim; d++) 
              R_i += 
                1.0 * 
                std::pow(fe_v.get_cell()->diameter(), 
                         parameters.diffusion_power) * 
                (parameters.theta * grad_W[point][component_i][d] + 
                 (1.0 - parameters.theta) * grad_W_old[point][component_i][d]) * 
                fe_v.shape_grad_component(i, point, component_i)[d] * 
                fe_v.JxW(point); 

            R_i -= 
              (parameters.theta * forcing[point][component_i] + 
               (1.0 - parameters.theta) * forcing_old[point][component_i]) * 
              fe_v.shape_value_component(i, point, component_i) * 
              fe_v.JxW(point); 
          } 


        for (unsigned int k = 0; k < dofs_per_cell; ++k) 
          residual_derivatives[k] = R_i.fastAccessDx(k); 
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 
        right_hand_side(dof_indices[i]) -= R_i.val(); 
      } 
  } 


  template <int dim> 
  void ConservationLaw<dim>::assemble_face_term( 
    const unsigned int                          face_no, 
    const FEFaceValuesBase<dim> &               fe_v, 
    const FEFaceValuesBase<dim> &               fe_v_neighbor, 
    const std::vector<types::global_dof_index> &dof_indices, 
    const std::vector<types::global_dof_index> &dof_indices_neighbor, 
    const bool                                  external_face, 
    const unsigned int                          boundary_id, 
    const double                                face_diameter) 
  { 
    const unsigned int n_q_points    = fe_v.n_quadrature_points; 
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell; 

    std::vector<Sacado::Fad::DFad<double>> independent_local_dof_values( 
      dofs_per_cell), 
      independent_neighbor_dof_values(external_face == false ? dofs_per_cell : 
                                                               0); 

    const unsigned int n_independent_variables = 
      (external_face == false ? 2 * dofs_per_cell : dofs_per_cell); 

    for (unsigned int i = 0; i < dofs_per_cell; i++) 
      { 
        independent_local_dof_values[i] = current_solution(dof_indices[i]); 
        independent_local_dof_values[i].diff(i, n_independent_variables); 
      } 

    if (external_face == false) 
      for (unsigned int i = 0; i < dofs_per_cell; i++) 
        { 
          independent_neighbor_dof_values[i] = 
            current_solution(dof_indices_neighbor[i]); 
          independent_neighbor_dof_values[i].diff(i + dofs_per_cell, 
                                                  n_independent_variables); 
        } 


    Table<2, Sacado::Fad::DFad<double>> Wplus( 
      n_q_points, EulerEquations<dim>::n_components), 
      Wminus(n_q_points, EulerEquations<dim>::n_components); 
    Table<2, double> Wplus_old(n_q_points, EulerEquations<dim>::n_components), 
      Wminus_old(n_q_points, EulerEquations<dim>::n_components); 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      for (unsigned int i = 0; i < dofs_per_cell; ++i) 
        { 
          const unsigned int component_i = 
            fe_v.get_fe().system_to_component_index(i).first; 
          Wplus[q][component_i] += 
            independent_local_dof_values[i] * 
            fe_v.shape_value_component(i, q, component_i); 
          Wplus_old[q][component_i] += 
            old_solution(dof_indices[i]) * 
            fe_v.shape_value_component(i, q, component_i); 
        } 


    if (external_face == false) 
      { 
        for (unsigned int q = 0; q < n_q_points; ++q) 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const unsigned int component_i = 
                fe_v_neighbor.get_fe().system_to_component_index(i).first; 
              Wminus[q][component_i] += 
                independent_neighbor_dof_values[i] * 
                fe_v_neighbor.shape_value_component(i, q, component_i); 
              Wminus_old[q][component_i] += 
                old_solution(dof_indices_neighbor[i]) * 
                fe_v_neighbor.shape_value_component(i, q, component_i); 
            } 
      } 




    else 
      { 
        Assert(boundary_id < Parameters::AllParameters<dim>::max_n_boundaries, 
               ExcIndexRange(boundary_id, 
                             0, 
                             Parameters::AllParameters<dim>::max_n_boundaries)); 

        std::vector<Vector<double>> boundary_values( 
          n_q_points, Vector<double>(EulerEquations<dim>::n_components)); 
        parameters.boundary_conditions[boundary_id].values.vector_value_list( 
          fe_v.get_quadrature_points(), boundary_values); 

        for (unsigned int q = 0; q < n_q_points; q++) 
          { 
            EulerEquations<dim>::compute_Wminus( 
              parameters.boundary_conditions[boundary_id].kind, 
              fe_v.normal_vector(q), 
              Wplus[q], 
              boundary_values[q], 
              Wminus[q]); 


            EulerEquations<dim>::compute_Wminus( 
              parameters.boundary_conditions[boundary_id].kind, 
              fe_v.normal_vector(q), 
              Wplus_old[q], 
              boundary_values[q], 
              Wminus_old[q]); 
          } 
      } 


    std::vector< 
      std::array<Sacado::Fad::DFad<double>, EulerEquations<dim>::n_components>> 
      normal_fluxes(n_q_points); 
    std::vector<std::array<double, EulerEquations<dim>::n_components>> 
      normal_fluxes_old(n_q_points); 

    double alpha; 

 
      { 
        case Parameters::Flux::constant: 
          alpha = parameters.stabilization_value; 
          break; 
        case Parameters::Flux::mesh_dependent: 
          alpha = face_diameter / (2.0 * parameters.time_step); 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
          alpha = 1; 
      } 

    for (unsigned int q = 0; q < n_q_points; ++q) 
      { 
        EulerEquations<dim>::numerical_normal_flux( 
          fe_v.normal_vector(q), Wplus[q], Wminus[q], alpha, normal_fluxes[q]); 
        EulerEquations<dim>::numerical_normal_flux(fe_v.normal_vector(q), 
                                                   Wplus_old[q], 
                                                   Wminus_old[q], 
                                                   alpha, 
                                                   normal_fluxes_old[q]); 
      } 


    std::vector<double> residual_derivatives(dofs_per_cell); 
    for (unsigned int i = 0; i < fe_v.dofs_per_cell; ++i) 
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true) 
        { 
          Sacado::Fad::DFad<double> R_i = 0; 

          for (unsigned int point = 0; point < n_q_points; ++point) 
            { 
              const unsigned int component_i = 
                fe_v.get_fe().system_to_component_index(i).first; 

              R_i += (parameters.theta * normal_fluxes[point][component_i] + 
                      (1.0 - parameters.theta) * 
                        normal_fluxes_old[point][component_i]) * 
                     fe_v.shape_value_component(i, point, component_i) * 
                     fe_v.JxW(point); 
            } 

          for (unsigned int k = 0; k < dofs_per_cell; ++k) 
            residual_derivatives[k] = R_i.fastAccessDx(k); 
          system_matrix.add(dof_indices[i], dof_indices, residual_derivatives); 

          if (external_face == false) 
            { 
              for (unsigned int k = 0; k < dofs_per_cell; ++k) 
                residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell + k); 
              system_matrix.add(dof_indices[i], 
                                dof_indices_neighbor, 
                                residual_derivatives); 
            } 

          right_hand_side(dof_indices[i]) -= R_i.val(); 
        } 
  } 


  template <int dim> 
  std::pair<unsigned int, double> 
  ConservationLaw<dim>::solve(Vector<double> &newton_update) 
  { 
    switch (parameters.solver) 
      { 


        case Parameters::Solver::direct: 
          { 
            SolverControl                                  solver_control(1, 0); 
            TrilinosWrappers::SolverDirect::AdditionalData data( 
              parameters.output == Parameters::Solver::verbose); 
            TrilinosWrappers::SolverDirect direct(solver_control, data); 

            direct.solve(system_matrix, newton_update, right_hand_side); 

            return {solver_control.last_step(), solver_control.last_value()}; 
          } 




        case Parameters::Solver::gmres: 
          { 
            Epetra_Vector x(View, 
                            system_matrix.trilinos_matrix().DomainMap(), 
                            newton_update.begin()); 
            Epetra_Vector b(View, 
                            system_matrix.trilinos_matrix().RangeMap(), 
                            right_hand_side.begin()); 

            AztecOO solver; 
            solver.SetAztecOption( 
              AZ_output, 
              (parameters.output == Parameters::Solver::quiet ? AZ_none : 
                                                                AZ_all)); 
            solver.SetAztecOption(AZ_solver, AZ_gmres); 
            solver.SetRHS(&b); 
            solver.SetLHS(&x); 

            solver.SetAztecOption(AZ_precond, AZ_dom_decomp); 
            solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut); 
            solver.SetAztecOption(AZ_overlap, 0); 
            solver.SetAztecOption(AZ_reorder, 0); 

            solver.SetAztecParam(AZ_drop, parameters.ilut_drop); 
            solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill); 
            solver.SetAztecParam(AZ_athresh, parameters.ilut_atol); 
            solver.SetAztecParam(AZ_rthresh, parameters.ilut_rtol); 

            solver.SetUserMatrix( 
              const_cast<Epetra_CrsMatrix *>(&system_matrix.trilinos_matrix())); 

            solver.Iterate(parameters.max_iterations, 
                           parameters.linear_residual); 

            return {solver.NumIters(), solver.TrueResidual()}; 
          } 
      } 

    Assert(false, ExcNotImplemented()); 
    return {0, 0}; 
  } 


  template <int dim> 
  void ConservationLaw<dim>::compute_refinement_indicators( 
    Vector<double> &refinement_indicators) const 
  { 
    EulerEquations<dim>::compute_refinement_indicators(dof_handler, 
                                                       mapping, 
                                                       predictor, 
                                                       refinement_indicators); 
  } 



  template <int dim> 
  void 
  ConservationLaw<dim>::refine_grid(const Vector<double> &refinement_indicators) 
  { 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        const unsigned int cell_no = cell->active_cell_index(); 
        cell->clear_coarsen_flag(); 
        cell->clear_refine_flag(); 

        if ((cell->level() < parameters.shock_levels) && 
            (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val)) 
          cell->set_refine_flag(); 
        else if ((cell->level() > 0) && 
                 (std::fabs(refinement_indicators(cell_no)) < 
                  0.75 * parameters.shock_val)) 
          cell->set_coarsen_flag(); 
      } 


    std::vector<Vector<double>> transfer_in; 
    std::vector<Vector<double>> transfer_out; 

    transfer_in.push_back(old_solution); 
    transfer_in.push_back(predictor); 

    triangulation.prepare_coarsening_and_refinement(); 

    SolutionTransfer<dim> soltrans(dof_handler); 
    soltrans.prepare_for_coarsening_and_refinement(transfer_in); 

    triangulation.execute_coarsening_and_refinement(); 

    dof_handler.clear(); 
    dof_handler.distribute_dofs(fe); 

    { 
      Vector<double> new_old_solution(1); 
      Vector<double> new_predictor(1); 

      transfer_out.push_back(new_old_solution); 
      transfer_out.push_back(new_predictor); 
      transfer_out[0].reinit(dof_handler.n_dofs()); 
      transfer_out[1].reinit(dof_handler.n_dofs()); 
    } 

    soltrans.interpolate(transfer_in, transfer_out); 

 
    old_solution = transfer_out[0]; 

    predictor.reinit(transfer_out[1].size()); 
    predictor = transfer_out[1]; 

    current_solution.reinit(dof_handler.n_dofs()); 
    current_solution = old_solution; 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 



  template <int dim> 
  void ConservationLaw<dim>::output_results() const 
  { 
    typename EulerEquations<dim>::Postprocessor postprocessor( 
      parameters.schlieren_plot); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 

    data_out.add_data_vector(current_solution, 
                             EulerEquations<dim>::component_names(), 
                             DataOut<dim>::type_dof_data, 
                             EulerEquations<dim>::component_interpretation()); 

    data_out.add_data_vector(current_solution, postprocessor); 

    data_out.build_patches(); 

    static unsigned int output_file_number = 0; 
    std::string         filename = 
      "solution-" + Utilities::int_to_string(output_file_number, 3) + ".vtk"; 
    std::ofstream output(filename); 
    data_out.write_vtk(output); 

    ++output_file_number; 
  } 




  template <int dim> 
  void ConservationLaw<dim>::run() 
  { 
    { 
      GridIn<dim> grid_in; 
      grid_in.attach_triangulation(triangulation); 

      std::ifstream input_file(parameters.mesh_filename); 
      Assert(input_file, ExcFileNotOpen(parameters.mesh_filename.c_str())); 

      grid_in.read_ucd(input_file); 
    } 

    dof_handler.clear(); 
    dof_handler.distribute_dofs(fe); 


    old_solution.reinit(dof_handler.n_dofs()); 
    current_solution.reinit(dof_handler.n_dofs()); 
    predictor.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 

    setup_system(); 

    VectorTools::interpolate(dof_handler, 
                             parameters.initial_conditions, 
                             old_solution); 
    current_solution = old_solution; 
    predictor        = old_solution; 

    if (parameters.do_refine == true) 
      for (unsigned int i = 0; i < parameters.shock_levels; ++i) 
        { 
          Vector<double> refinement_indicators(triangulation.n_active_cells()); 

          compute_refinement_indicators(refinement_indicators); 
          refine_grid(refinement_indicators); 

          setup_system(); 

          VectorTools::interpolate(dof_handler, 
                                   parameters.initial_conditions, 
                                   old_solution); 
          current_solution = old_solution; 
          predictor        = old_solution; 
        } 

    output_results(); 


    Vector<double> newton_update(dof_handler.n_dofs()); 

    double time        = 0; 
    double next_output = time + parameters.output_step; 

    predictor = old_solution; 
    while (time < parameters.final_time) 
      { 
        std::cout << "T=" << time << std::endl 
                  << "   Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl 
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl 
                  << std::endl; 

        std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl 
                  << "   _____________________________________" << std::endl; 



        unsigned int nonlin_iter = 0; 
        current_solution         = predictor; 
        while (true) 
          { 
            system_matrix = 0; 

            right_hand_side = 0; 
            assemble_system(); 

            const double res_norm = right_hand_side.l2_norm(); 
            if (std::fabs(res_norm) < 1e-10) 
              { 
                std::printf("   %-16.3e (converged)\n\n", res_norm); 
                break; 
              } 
            else 
              { 
                newton_update = 0; 

                std::pair<unsigned int, double> convergence = 
                  solve(newton_update); 

                current_solution += newton_update; 

                std::printf("   %-16.3e %04d        %-5.2e\n", 
                            res_norm, 
                            convergence.first, 
                            convergence.second); 
              } 

            ++nonlin_iter; 
            AssertThrow(nonlin_iter <= 10, 
                        ExcMessage("No convergence in nonlinear solver")); 
          } 



        time += parameters.time_step; 

        if (parameters.output_step < 0) 
          output_results(); 
        else if (time >= next_output) 
          { 
            output_results(); 
            next_output += parameters.output_step; 
          } 

        predictor = current_solution; 
        predictor.sadd(2.0, -1.0, old_solution); 

        old_solution = current_solution; 

        if (parameters.do_refine == true) 
          { 
            Vector<double> refinement_indicators( 
              triangulation.n_active_cells()); 
            compute_refinement_indicators(refinement_indicators); 

            refine_grid(refinement_indicators); 
            setup_system(); 

            newton_update.reinit(dof_handler.n_dofs()); 
          } 
      } 
  } 
} // namespace Step33 


int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step33; 

      if (argc != 2) 
        { 
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl; 
          std::exit(1); 
        } 

      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, dealii::numbers::invalid_unsigned_int); 

      ConservationLaw<2> cons(argv[1]); 
      cons.run(); 
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


