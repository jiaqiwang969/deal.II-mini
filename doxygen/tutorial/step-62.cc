

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
 * Author: Daniel Garcia-Sanchez, CNRS, 2019 
 */ 




#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/function.h> 

#include <deal.II/base/index_set.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/generic_linear_algebra.h> 
#include <deal.II/lac/petsc_solver.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

#include <fstream> 
#include <iostream> 


#include <deal.II/base/tensor.h> 


#include <deal.II/base/hdf5.h> 


#include <deal.II/numerics/vector_tools.h> 

#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/grid_tools_cache.h> 

namespace step62 
{ 
  using namespace dealii; 


  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide(HDF5::Group &data); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component) const override; 

  private: 


    HDF5::Group data; 


    const double     max_force_amplitude; 
    const double     force_sigma_x; 
    const double     force_sigma_y; 
    const double     max_force_width_x; 
    const double     max_force_width_y; 
    const Point<dim> force_center; 

  public: 


    const unsigned int force_component = 0; 
  }; 

  template <int dim> 
  class PML : public Function<dim, std::complex<double>> 
  { 
  public: 
    PML(HDF5::Group &data); 

    virtual std::complex<double> 
    value(const Point<dim> &p, const unsigned int component) const override; 

  private: 

    HDF5::Group data; 


    const double pml_coeff; 
    const int    pml_coeff_degree; 
    const double dimension_x; 
    const double dimension_y; 
    const bool   pml_x; 
    const bool   pml_y; 
    const double pml_width_x; 
    const double pml_width_y; 
    const double a_coeff_x; 
    const double a_coeff_y; 
  }; 


  template <int dim> 
  class Rho : public Function<dim> 
  { 
  public: 
    Rho(HDF5::Group &data); 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

  private: 

    HDF5::Group data; 


    const double       lambda; 
    const double       mu; 
    const double       material_a_rho; 
    const double       material_b_rho; 
    const double       cavity_resonance_frequency; 
    const unsigned int nb_mirror_pairs; 
    const double       dimension_y; 
    const unsigned int grid_level; 
    double             average_rho_width; 
  }; 


  template <int dim> 
  class Parameters 
  { 
  public: 
    Parameters(HDF5::Group &data); 

    HDF5::Group data; 


    const std::string        simulation_name; 
    const bool               save_vtu_files; 
    const double             start_frequency; 
    const double             stop_frequency; 
    const unsigned int       nb_frequency_points; 
    const double             lambda; 
    const double             mu; 
    const double             dimension_x; 
    const double             dimension_y; 
    const unsigned int       nb_probe_points; 
    const unsigned int       grid_level; 
    const Point<dim>         probe_start_point; 
    const Point<dim>         probe_stop_point; 
    const RightHandSide<dim> right_hand_side; 
    const PML<dim>           pml; 
    const Rho<dim>           rho; 

  private: 
    const double comparison_float_constant = 1e-12; 
  }; 


  template <int dim> 
  class QuadratureCache 
  { 
  public: 
    QuadratureCache(const unsigned int dofs_per_cell); 

  private: 
    unsigned int dofs_per_cell; 

  public: 


    FullMatrix<std::complex<double>>  mass_coefficient; 
    FullMatrix<std::complex<double>>  stiffness_coefficient; 
    std::vector<std::complex<double>> right_hand_side; 
    double                            JxW; 
  }; 



  template <int dim> 
  SymmetricTensor<4, dim> get_stiffness_tensor(const double lambda, 
                                               const double mu) 
  { 
    SymmetricTensor<4, dim> stiffness_tensor; 
    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = 0; j < dim; ++j) 
        for (unsigned int k = 0; k < dim; ++k) 
          for (unsigned int l = 0; l < dim; ++l) 
            stiffness_tensor[i][j][k][l] = 
              (((i == k) && (j == l) ? mu : 0.0) + 
               ((i == l) && (j == k) ? mu : 0.0) + 
               ((i == j) && (k == l) ? lambda : 0.0)); 
    return stiffness_tensor; 
  } 






  template <int dim> 
  class ElasticWave 
  { 
  public: 
    ElasticWave(const Parameters<dim> &parameters); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(const double omega, 
                         const bool   calculate_quadrature_data); 
    void solve(); 
    void initialize_probe_positions_vector(); 
    void store_frequency_step_data(const unsigned int frequency_idx); 
    void output_results(); 


    void setup_quadrature_cache(); 


    void frequency_sweep(); 


    Parameters<dim> parameters; 

    MPI_Comm mpi_communicator; 

    parallel::distributed::Triangulation<dim> triangulation; 

    QGauss<dim> quadrature_formula; 


    std::vector<QuadratureCache<dim>> quadrature_cache; 

    FESystem<dim>   fe; 
    DoFHandler<dim> dof_handler; 

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

    AffineConstraints<std::complex<double>> constraints; 

    LinearAlgebraPETSc::MPI::SparseMatrix system_matrix; 
    LinearAlgebraPETSc::MPI::Vector       locally_relevant_solution; 
    LinearAlgebraPETSc::MPI::Vector       system_rhs; 


    std::vector<double> frequency; 


    FullMatrix<double> probe_positions; 


    HDF5::DataSet frequency_dataset; 
    HDF5::DataSet probe_positions_dataset; 


    HDF5::DataSet displacement; 

    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 
  }; 



  template <int dim> 
  RightHandSide<dim>::RightHandSide(HDF5::Group &data) 
    : Function<dim>(dim) 
    , data(data) 
    , max_force_amplitude(data.get_attribute<double>("max_force_amplitude")) 
    , force_sigma_x(data.get_attribute<double>("force_sigma_x")) 
    , force_sigma_y(data.get_attribute<double>("force_sigma_y")) 
    , max_force_width_x(data.get_attribute<double>("max_force_width_x")) 
    , max_force_width_y(data.get_attribute<double>("max_force_width_y")) 
    , force_center(Point<dim>(data.get_attribute<double>("force_x_pos"), 
                              data.get_attribute<double>("force_y_pos"))) 
  {} 


  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> & p, 
                                   const unsigned int component) const 
  { 
    if (component == force_component) 
      { 
        if (std::abs(p[0] - force_center[0]) < max_force_width_x / 2 && 
            std::abs(p[1] - force_center[1]) < max_force_width_y / 2) 
          { 
            return max_force_amplitude * 
                   std::exp(-(std::pow(p[0] - force_center[0], 2) / 
                                (2 * std::pow(force_sigma_x, 2)) + 
                              std::pow(p[1] - force_center[1], 2) / 
                                (2 * std::pow(force_sigma_y, 2)))); 
          } 
        else 
          { 
            return 0; 
          } 
      } 
    else 
      { 
        return 0; 
      } 
  } 



  template <int dim> 
  PML<dim>::PML(HDF5::Group &data) 
    : Function<dim, std::complex<double>>(dim) 
    , data(data) 
    , pml_coeff(data.get_attribute<double>("pml_coeff")) 
    , pml_coeff_degree(data.get_attribute<int>("pml_coeff_degree")) 
    , dimension_x(data.get_attribute<double>("dimension_x")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , pml_x(data.get_attribute<bool>("pml_x")) 
    , pml_y(data.get_attribute<bool>("pml_y")) 
    , pml_width_x(data.get_attribute<double>("pml_width_x")) 
    , pml_width_y(data.get_attribute<double>("pml_width_y")) 
    , a_coeff_x(pml_coeff / std::pow(pml_width_x, pml_coeff_degree)) 
    , a_coeff_y(pml_coeff / std::pow(pml_width_y, pml_coeff_degree)) 
  {} 

  template <int dim> 
  std::complex<double> PML<dim>::value(const Point<dim> & p, 
                                       const unsigned int component) const 
  { 
    double calculated_pml_x_coeff = 0; 
    double calculated_pml_y_coeff = 0; 

    if ((component == 0) && pml_x) 
      { 
        const double pml_x_start_position = dimension_x / 2 - pml_width_x; 
        if (std::abs(p[0]) > pml_x_start_position) 
          { 
            const double x_prime = std::abs(p[0]) - pml_x_start_position; 
            calculated_pml_x_coeff = 
              a_coeff_x * std::pow(x_prime, pml_coeff_degree); 
          } 
      } 

    if ((component == 1) && pml_y) 
      { 
        const double pml_y_start_position = dimension_y / 2 - pml_width_y; 
        if (std::abs(p[1]) > pml_y_start_position) 
          { 
            const double y_prime = std::abs(p[1]) - pml_y_start_position; 
            calculated_pml_y_coeff = 
              a_coeff_y * std::pow(y_prime, pml_coeff_degree); 
          } 
      } 

    return 1. + std::max(calculated_pml_x_coeff, calculated_pml_y_coeff) * 
                  std::complex<double>(0., 1.); 
  } 



  template <int dim> 
  Rho<dim>::Rho(HDF5::Group &data) 
    : Function<dim>(1) 
    , data(data) 
    , lambda(data.get_attribute<double>("lambda")) 
    , mu(data.get_attribute<double>("mu")) 
    , material_a_rho(data.get_attribute<double>("material_a_rho")) 
    , material_b_rho(data.get_attribute<double>("material_b_rho")) 
    , cavity_resonance_frequency( 
        data.get_attribute<double>("cavity_resonance_frequency")) 
    , nb_mirror_pairs(data.get_attribute<int>("nb_mirror_pairs")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , grid_level(data.get_attribute<int>("grid_level")) 
  { 


    average_rho_width = dimension_y / (std::pow(2.0, grid_level)); 
    data.set_attribute("average_rho_width", average_rho_width); 
  } 

  template <int dim> 
  double Rho<dim>::value(const Point<dim> &p, 
                         const unsigned int /*component*/) const 
  { 


    double elastic_constant; 
    if (dim == 2) 
      { 
        elastic_constant = 4 * mu * (lambda + mu) / (lambda + 2 * mu); 
      } 
    else if (dim == 3) 
      { 
        elastic_constant = mu * (3 * lambda + 2 * mu) / (lambda + mu); 
      } 
    else 
      { 
        Assert(false, ExcInternalError()); 
      } 
    const double material_a_speed_of_sound = 
      std::sqrt(elastic_constant / material_a_rho); 
    const double material_a_wavelength = 
      material_a_speed_of_sound / cavity_resonance_frequency; 
    const double material_b_speed_of_sound = 
      std::sqrt(elastic_constant / material_b_rho); 
    const double material_b_wavelength = 
      material_b_speed_of_sound / cavity_resonance_frequency; 


    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_transition_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4); 
        if (std::abs(p[0]) >= 
              (layer_transition_center - average_rho_width / 2) && 
            std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2)) 
          { 
            const double coefficient = 
              (std::abs(p[0]) - 
               (layer_transition_center - average_rho_width / 2)) / 
              average_rho_width; 
            return (1 - coefficient) * material_a_rho + 
                   coefficient * material_b_rho; 
          } 
      } 


    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_transition_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 4; 
        if (std::abs(p[0]) >= 
              (layer_transition_center - average_rho_width / 2) && 
            std::abs(p[0]) <= (layer_transition_center + average_rho_width / 2)) 
          { 
            const double coefficient = 
              (std::abs(p[0]) - 
               (layer_transition_center - average_rho_width / 2)) / 
              average_rho_width; 
            return (1 - coefficient) * material_b_rho + 
                   coefficient * material_a_rho; 
          } 
      } 


    if (std::abs(p[0]) <= material_a_wavelength / 2) 
      { 
        return material_a_rho; 
      } 


    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 4 + material_a_wavelength / 8; 
        const double layer_width = material_a_wavelength / 4; 
        if (std::abs(p[0]) >= (layer_center - layer_width / 2) && 
            std::abs(p[0]) <= (layer_center + layer_width / 2)) 
          { 
            return material_a_rho; 
          } 
      } 


    for (unsigned int idx = 0; idx < nb_mirror_pairs; idx++) 
      { 
        const double layer_center = 
          material_a_wavelength / 2 + 
          idx * (material_b_wavelength / 4 + material_a_wavelength / 4) + 
          material_b_wavelength / 8; 
        const double layer_width = material_b_wavelength / 4; 
        if (std::abs(p[0]) >= (layer_center - layer_width / 2) && 
            std::abs(p[0]) <= (layer_center + layer_width / 2)) 
          { 
            return material_b_rho; 
          } 
      } 


    return material_a_rho; 
  } 



  template <int dim> 
  Parameters<dim>::Parameters(HDF5::Group &data) 
    : data(data) 
    , simulation_name(data.get_attribute<std::string>("simulation_name")) 
    , save_vtu_files(data.get_attribute<bool>("save_vtu_files")) 
    , start_frequency(data.get_attribute<double>("start_frequency")) 
    , stop_frequency(data.get_attribute<double>("stop_frequency")) 
    , nb_frequency_points(data.get_attribute<int>("nb_frequency_points")) 
    , lambda(data.get_attribute<double>("lambda")) 
    , mu(data.get_attribute<double>("mu")) 
    , dimension_x(data.get_attribute<double>("dimension_x")) 
    , dimension_y(data.get_attribute<double>("dimension_y")) 
    , nb_probe_points(data.get_attribute<int>("nb_probe_points")) 
    , grid_level(data.get_attribute<int>("grid_level")) 
    , probe_start_point(data.get_attribute<double>("probe_pos_x"), 
                        data.get_attribute<double>("probe_pos_y") - 
                          data.get_attribute<double>("probe_width_y") / 2) 
    , probe_stop_point(data.get_attribute<double>("probe_pos_x"), 
                       data.get_attribute<double>("probe_pos_y") + 
                         data.get_attribute<double>("probe_width_y") / 2) 
    , right_hand_side(data) 
    , pml(data) 
    , rho(data) 
  {} 



  template <int dim> 
  QuadratureCache<dim>::QuadratureCache(const unsigned int dofs_per_cell) 
    : dofs_per_cell(dofs_per_cell) 
    , mass_coefficient(dofs_per_cell, dofs_per_cell) 
    , stiffness_coefficient(dofs_per_cell, dofs_per_cell) 
    , right_hand_side(dofs_per_cell) 
  {} 



  template <int dim> 
  ElasticWave<dim>::ElasticWave(const Parameters<dim> &parameters) 
    : parameters(parameters) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , triangulation(mpi_communicator, 
                    typename Triangulation<dim>::MeshSmoothing( 
                      Triangulation<dim>::smoothing_on_refinement | 
                      Triangulation<dim>::smoothing_on_coarsening)) 
    , quadrature_formula(2) 
    , fe(FE_Q<dim>(1), dim) 
    , dof_handler(triangulation) 
    , frequency(parameters.nb_frequency_points) 
    , probe_positions(parameters.nb_probe_points, dim) 
    , frequency_dataset(parameters.data.template create_dataset<double>( 
        "frequency", 
        std::vector<hsize_t>{parameters.nb_frequency_points})) 
    , probe_positions_dataset(parameters.data.template create_dataset<double>( 
        "position", 
        std::vector<hsize_t>{parameters.nb_probe_points, dim})) 
    , displacement( 
        parameters.data.template create_dataset<std::complex<double>>( 
          "displacement", 
          std::vector<hsize_t>{parameters.nb_probe_points, 
                               parameters.nb_frequency_points})) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(mpi_communicator, 
                      pcout, 
                      TimerOutput::summary, 
                      TimerOutput::wall_times) 
  {} 



  template <int dim> 
  void ElasticWave<dim>::setup_system() 
  { 
    TimerOutput::Scope t(computing_timer, "setup"); 

    dof_handler.distribute_dofs(fe); 

    locally_owned_dofs = dof_handler.locally_owned_dofs(); 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

    locally_relevant_solution.reinit(locally_owned_dofs, 
                                     locally_relevant_dofs, 
                                     mpi_communicator); 

    system_rhs.reinit(locally_owned_dofs, mpi_communicator); 

    constraints.clear(); 
    constraints.reinit(locally_relevant_dofs); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    constraints.close(); 

    DynamicSparsityPattern dsp(locally_relevant_dofs); 

    DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
    SparsityTools::distribute_sparsity_pattern(dsp, 
                                               locally_owned_dofs, 
                                               mpi_communicator, 
                                               locally_relevant_dofs); 

    system_matrix.reinit(locally_owned_dofs, 
                         locally_owned_dofs, 
                         dsp, 
                         mpi_communicator); 
  } 



  template <int dim> 
  void ElasticWave<dim>::assemble_system(const double omega, 
                                         const bool   calculate_quadrature_data) 
  { 
    TimerOutput::Scope t(computing_timer, "assembly"); 

    FEValues<dim>      fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 
    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<std::complex<double>> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<std::complex<double>>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 


    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim)); 
    std::vector<double>         rho_values(n_q_points); 
    std::vector<Vector<std::complex<double>>> pml_values( 
      n_q_points, Vector<std::complex<double>>(dim)); 


    const SymmetricTensor<4, dim> stiffness_tensor = 
      get_stiffness_tensor<dim>(parameters.lambda, parameters.mu); 


    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell_matrix = 0; 
          cell_rhs    = 0; 


          if (calculate_quadrature_data) 
            { 
              fe_values.reinit(cell); 

              parameters.right_hand_side.vector_value_list( 
                fe_values.get_quadrature_points(), rhs_values); 
              parameters.rho.value_list(fe_values.get_quadrature_points(), 
                                        rho_values); 
              parameters.pml.vector_value_list( 
                fe_values.get_quadrature_points(), pml_values); 
            } 


          QuadratureCache<dim> *local_quadrature_points_data = 
            reinterpret_cast<QuadratureCache<dim> *>(cell->user_pointer()); 
          Assert(local_quadrature_points_data >= &quadrature_cache.front(), 
                 ExcInternalError()); 
          Assert(local_quadrature_points_data <= &quadrature_cache.back(), 
                 ExcInternalError()); 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 


              QuadratureCache<dim> &quadrature_data = 
                local_quadrature_points_data[q]; 


              Tensor<1, dim>                       force; 
              Tensor<1, dim, std::complex<double>> s; 
              std::complex<double>                 xi(1, 0); 


              if (calculate_quadrature_data) 
                { 


                  quadrature_data.JxW = fe_values.JxW(q); 

                  for (unsigned int component = 0; component < dim; ++component) 
                    { 


                      force[component] = rhs_values[q][component]; 
                      s[component]     = pml_values[q][component]; 
                      xi *= s[component]; 
                    } 


                  Tensor<4, dim, std::complex<double>> alpha; 
                  Tensor<4, dim, std::complex<double>> beta; 
                  for (unsigned int m = 0; m < dim; ++m) 
                    for (unsigned int n = 0; n < dim; ++n) 
                      for (unsigned int k = 0; k < dim; ++k) 
                        for (unsigned int l = 0; l < dim; ++l) 
                          { 
                            alpha[m][n][k][l] = xi * 
                                                stiffness_tensor[m][n][k][l] / 
                                                (2.0 * s[n] * s[k]); 
                            beta[m][n][k][l] = xi * 
                                               stiffness_tensor[m][n][k][l] / 
                                               (2.0 * s[n] * s[l]); 
                          } 

                  for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                    { 
                      const Tensor<1, dim> phi_i = 
                        fe_values[displacement].value(i, q); 
                      const Tensor<2, dim> grad_phi_i = 
                        fe_values[displacement].gradient(i, q); 

                      for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                        { 
                          const Tensor<1, dim> phi_j = 
                            fe_values[displacement].value(j, q); 
                          const Tensor<2, dim> grad_phi_j = 
                            fe_values[displacement].gradient(j, q); 


                          quadrature_data.mass_coefficient[i][j] = 
                            rho_values[q] * xi * phi_i * phi_j; 


                          std::complex<double> stiffness_coefficient = 0; 
                          for (unsigned int m = 0; m < dim; ++m) 
                            for (unsigned int n = 0; n < dim; ++n) 
                              for (unsigned int k = 0; k < dim; ++k) 
                                for (unsigned int l = 0; l < dim; ++l) 
                                  { 


                                    stiffness_coefficient += 
                                      grad_phi_i[m][n] * 
                                      (alpha[m][n][k][l] * grad_phi_j[l][k] + 
                                       beta[m][n][k][l] * grad_phi_j[k][l]); 
                                  } 


                          quadrature_data.stiffness_coefficient[i][j] = 
                            stiffness_coefficient; 
                        } 


 
                        phi_i * force * fe_values.JxW(q); 
                    } 
                } 


              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    { 
                      std::complex<double> matrix_sum = 0; 
                      matrix_sum += -std::pow(omega, 2) * 
                                    quadrature_data.mass_coefficient[i][j]; 
                      matrix_sum += quadrature_data.stiffness_coefficient[i][j]; 
                      cell_matrix(i, j) += matrix_sum * quadrature_data.JxW; 
                    } 
                  cell_rhs(i) += quadrature_data.right_hand_side[i]; 
                } 
            } 
          cell->get_dof_indices(local_dof_indices); 
          constraints.distribute_local_to_global(cell_matrix, 
                                                 cell_rhs, 
                                                 local_dof_indices, 
                                                 system_matrix, 
                                                 system_rhs); 
        } 

    system_matrix.compress(VectorOperation::add); 
    system_rhs.compress(VectorOperation::add); 
  } 


  template <int dim> 
  void ElasticWave<dim>::solve() 
  { 
    TimerOutput::Scope              t(computing_timer, "solve"); 
    LinearAlgebraPETSc::MPI::Vector completely_distributed_solution( 
      locally_owned_dofs, mpi_communicator); 

    SolverControl                    solver_control; 
    PETScWrappers::SparseDirectMUMPS solver(solver_control, mpi_communicator); 
    solver.solve(system_matrix, completely_distributed_solution, system_rhs); 

    pcout << "   Solved in " << solver_control.last_step() << " iterations." 
          << std::endl; 
    constraints.distribute(completely_distributed_solution); 
    locally_relevant_solution = completely_distributed_solution; 
  } 


  template <int dim> 
  void ElasticWave<dim>::initialize_probe_positions_vector() 
  { 
    for (unsigned int position_idx = 0; 
         position_idx < parameters.nb_probe_points; 
         ++position_idx) 
      { 



        const Point<dim> p = 
          (position_idx / ((double)(parameters.nb_probe_points - 1))) * 
            (parameters.probe_stop_point + (-parameters.probe_start_point)) + 
          parameters.probe_start_point; 
        probe_positions[position_idx][0] = p[0]; 
        probe_positions[position_idx][1] = p[1]; 
        if (dim == 3) 
          { 
            probe_positions[position_idx][2] = p[2]; 
          } 
      } 
  } 


  template <int dim> 
  void 
  ElasticWave<dim>::store_frequency_step_data(const unsigned int frequency_idx) 
  { 
    TimerOutput::Scope t(computing_timer, "store_frequency_step_data"); 


    const unsigned int probe_displacement_component = 0; 


    std::vector<hsize_t>              coordinates; 
    std::vector<std::complex<double>> displacement_data; 

    const auto &mapping = get_default_linear_mapping(triangulation); 
    GridTools::Cache<dim, dim> cache(triangulation, mapping); 
    typename Triangulation<dim, dim>::active_cell_iterator cell_hint{}; 
    std::vector<bool>                                      marked_vertices = {}; 
    const double                                           tolerance = 1.e-10; 

    for (unsigned int position_idx = 0; 
         position_idx < parameters.nb_probe_points; 
         ++position_idx) 
      { 
        Point<dim> point; 
        for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
          { 
            point[dim_idx] = probe_positions[position_idx][dim_idx]; 
          } 
        bool point_in_locally_owned_cell = false; 
        { 
          auto cell_and_ref_point = GridTools::find_active_cell_around_point( 
            cache, point, cell_hint, marked_vertices, tolerance); 
          if (cell_and_ref_point.first.state() == IteratorState::valid) 
            { 
              cell_hint = cell_and_ref_point.first; 
              point_in_locally_owned_cell = 
                cell_and_ref_point.first->is_locally_owned(); 
            } 
        } 
        if (point_in_locally_owned_cell) 
          { 


            Vector<std::complex<double>> tmp_vector(dim); 
            VectorTools::point_value(dof_handler, 
                                     locally_relevant_solution, 
                                     point, 
                                     tmp_vector); 
            coordinates.emplace_back(position_idx); 
            coordinates.emplace_back(frequency_idx); 
            displacement_data.emplace_back( 
              tmp_vector(probe_displacement_component)); 
          } 
      } 


    if (coordinates.size() > 0) 
      { 
        displacement.write_selection(displacement_data, coordinates); 
      } 

    else 
      { 
        displacement.write_none<std::complex<double>>(); 
      } 


    if (parameters.save_vtu_files) 
      { 
        std::vector<std::string> solution_names(dim, "displacement"); 
        std::vector<DataComponentInterpretation::DataComponentInterpretation> 
          interpretation( 
            dim, DataComponentInterpretation::component_is_part_of_vector); 

        DataOut<dim> data_out; 
        data_out.add_data_vector(dof_handler, 
                                 locally_relevant_solution, 
                                 solution_names, 
                                 interpretation); 
        Vector<float> subdomain(triangulation.n_active_cells()); 
        for (unsigned int i = 0; i < subdomain.size(); ++i) 
          subdomain(i) = triangulation.locally_owned_subdomain(); 
        data_out.add_data_vector(subdomain, "subdomain"); 

        std::vector<Vector<double>> force( 
          dim, Vector<double>(triangulation.n_active_cells())); 
        std::vector<Vector<double>> pml( 
          dim, Vector<double>(triangulation.n_active_cells())); 
        Vector<double> rho(triangulation.n_active_cells()); 

        for (auto &cell : triangulation.active_cell_iterators()) 
          { 
            if (cell->is_locally_owned()) 
              { 
                for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
                  { 
                    force[dim_idx](cell->active_cell_index()) = 
                      parameters.right_hand_side.value(cell->center(), dim_idx); 
                    pml[dim_idx](cell->active_cell_index()) = 
                      parameters.pml.value(cell->center(), dim_idx).imag(); 
                  } 
                rho(cell->active_cell_index()) = 
                  parameters.rho.value(cell->center()); 
              } 


            else 
              { 
                for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
                  { 
                    force[dim_idx](cell->active_cell_index()) = -1e+20; 
                    pml[dim_idx](cell->active_cell_index())   = -1e+20; 
                  } 
                rho(cell->active_cell_index()) = -1e+20; 
              } 
          } 

        for (unsigned int dim_idx = 0; dim_idx < dim; ++dim_idx) 
          { 
            data_out.add_data_vector(force[dim_idx], 
                                     "force_" + std::to_string(dim_idx)); 
            data_out.add_data_vector(pml[dim_idx], 
                                     "pml_" + std::to_string(dim_idx)); 
          } 
        data_out.add_data_vector(rho, "rho"); 

        data_out.build_patches(); 

        std::stringstream  frequency_idx_stream; 
        const unsigned int nb_number_positions = 
          ((unsigned int)std::log10(parameters.nb_frequency_points)) + 1; 
        frequency_idx_stream << std::setw(nb_number_positions) 
                             << std::setfill('0') << frequency_idx; 
        std::string filename = (parameters.simulation_name + "_" + 
                                frequency_idx_stream.str() + ".vtu"); 
        data_out.write_vtu_in_parallel(filename.c_str(), mpi_communicator); 
      } 
  } 



  template <int dim> 
  void ElasticWave<dim>::output_results() 
  { 

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
      { 
        frequency_dataset.write(frequency); 
        probe_positions_dataset.write(probe_positions); 
      } 
    else 
      { 
        frequency_dataset.write_none<double>(); 
        probe_positions_dataset.write_none<double>(); 
      } 
  } 



  template <int dim> 
  void ElasticWave<dim>::setup_quadrature_cache() 
  { 
    triangulation.clear_user_data(); 

    { 
      std::vector<QuadratureCache<dim>> tmp; 
      quadrature_cache.swap(tmp); 
    } 

    quadrature_cache.resize(triangulation.n_locally_owned_active_cells() * 
                              quadrature_formula.size(), 
                            QuadratureCache<dim>(fe.n_dofs_per_cell())); 
    unsigned int cache_index = 0; 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell->set_user_pointer(&quadrature_cache[cache_index]); 
          cache_index += quadrature_formula.size(); 
        } 
    Assert(cache_index == quadrature_cache.size(), ExcInternalError()); 
  } 



  template <int dim> 
  void ElasticWave<dim>::frequency_sweep() 
  { 
    for (unsigned int frequency_idx = 0; 
         frequency_idx < parameters.nb_frequency_points; 
         ++frequency_idx) 
      { 
        pcout << parameters.simulation_name + " frequency idx: " 
              << frequency_idx << '/' << parameters.nb_frequency_points - 1 
              << std::endl; 

        setup_system(); 
        if (frequency_idx == 0) 
          { 
            pcout << "   Number of active cells :       " 
                  << triangulation.n_active_cells() << std::endl; 
            pcout << "   Number of degrees of freedom : " 
                  << dof_handler.n_dofs() << std::endl; 
          } 

        if (frequency_idx == 0) 
          { 


            parameters.data.set_attribute("active_cells", 
                                          triangulation.n_active_cells()); 
            parameters.data.set_attribute("degrees_of_freedom", 
                                          dof_handler.n_dofs()); 
          } 


        const double current_loop_frequency = 
          (parameters.start_frequency + 
           frequency_idx * 
             (parameters.stop_frequency - parameters.start_frequency) / 
             (parameters.nb_frequency_points - 1)); 
        const double current_loop_omega = 
          2 * numbers::PI * current_loop_frequency; 


        assemble_system(current_loop_omega, 
                        (frequency_idx == 0) ? true : false); 
        solve(); 

        frequency[frequency_idx] = current_loop_frequency; 
        store_frequency_step_data(frequency_idx); 

        computing_timer.print_summary(); 
        computing_timer.reset(); 
        pcout << std::endl; 
      } 
  } 



  template <int dim> 
  void ElasticWave<dim>::run() 
  { 
#ifdef DEBUG 
    pcout << "Debug mode" << std::endl; 
#else 
    pcout << "Release mode" << std::endl; 
#endif 

    { 
      Point<dim> p1; 
      p1(0) = -parameters.dimension_x / 2; 
      p1(1) = -parameters.dimension_y / 2; 
      if (dim == 3) 
        { 
          p1(2) = -parameters.dimension_y / 2; 
        } 
      Point<dim> p2; 
      p2(0) = parameters.dimension_x / 2; 
      p2(1) = parameters.dimension_y / 2; 
      if (dim == 3) 
        { 
          p2(2) = parameters.dimension_y / 2; 
        } 
      std::vector<unsigned int> divisions(dim); 
      divisions[0] = int(parameters.dimension_x / parameters.dimension_y); 
      divisions[1] = 1; 
      if (dim == 3) 
        { 
          divisions[2] = 1; 
        } 
      GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                                divisions, 
                                                p1, 
                                                p2); 
    } 

    triangulation.refine_global(parameters.grid_level); 

    setup_quadrature_cache(); 

    initialize_probe_positions_vector(); 

    frequency_sweep(); 

    output_results(); 
  } 
} // namespace step62 



int main(int argc, char *argv[]) 
{ 
  try 
    { 
      using namespace dealii; 
      const unsigned int dim = 2; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      HDF5::File data_file("results.h5", 
                           HDF5::File::FileAccessMode::create, 
                           MPI_COMM_WORLD); 
      auto       data = data_file.create_group("data"); 


      const std::vector<std::string> group_names = {"displacement", 
                                                    "calibration"}; 
      for (auto group_name : group_names) 
        { 













          auto group = data.create_group(group_name); 

          group.set_attribute<double>("dimension_x", 2e-5); 
          group.set_attribute<double>("dimension_y", 2e-8); 
          group.set_attribute<double>("probe_pos_x", 8e-6); 
          group.set_attribute<double>("probe_pos_y", 0); 
          group.set_attribute<double>("probe_width_y", 2e-08); 
          group.set_attribute<unsigned int>("nb_probe_points", 5); 
          group.set_attribute<unsigned int>("grid_level", 1); 
          group.set_attribute<double>("cavity_resonance_frequency", 20e9); 
          group.set_attribute<unsigned int>("nb_mirror_pairs", 15); 

          group.set_attribute<double>("poissons_ratio", 0.27); 
          group.set_attribute<double>("youngs_modulus", 270000000000.0); 
          group.set_attribute<double>("material_a_rho", 3200); 

          if (group_name == std::string("displacement")) 
            group.set_attribute<double>("material_b_rho", 2000); 
          else 
            group.set_attribute<double>("material_b_rho", 3200); 

          group.set_attribute( 
            "lambda", 
            group.get_attribute<double>("youngs_modulus") * 
              group.get_attribute<double>("poissons_ratio") / 
              ((1 + group.get_attribute<double>("poissons_ratio")) * 
               (1 - 2 * group.get_attribute<double>("poissons_ratio")))); 
          group.set_attribute("mu", 
                              group.get_attribute<double>("youngs_modulus") / 
                                (2 * (1 + group.get_attribute<double>( 
                                            "poissons_ratio")))); 

          group.set_attribute<double>("max_force_amplitude", 1e26); 
          group.set_attribute<double>("force_sigma_x", 1e-7); 
          group.set_attribute<double>("force_sigma_y", 1); 
          group.set_attribute<double>("max_force_width_x", 3e-7); 
          group.set_attribute<double>("max_force_width_y", 2e-8); 
          group.set_attribute<double>("force_x_pos", -8e-6); 
          group.set_attribute<double>("force_y_pos", 0); 

          group.set_attribute<bool>("pml_x", true); 
          group.set_attribute<bool>("pml_y", false); 
          group.set_attribute<double>("pml_width_x", 1.8e-6); 
          group.set_attribute<double>("pml_width_y", 5e-7); 
          group.set_attribute<double>("pml_coeff", 1.6); 
          group.set_attribute<unsigned int>("pml_coeff_degree", 2); 

          group.set_attribute<double>("center_frequency", 20e9); 
          group.set_attribute<double>("frequency_range", 0.5e9); 
          group.set_attribute<double>( 
            "start_frequency", 
            group.get_attribute<double>("center_frequency") - 
              group.get_attribute<double>("frequency_range") / 2); 
          group.set_attribute<double>( 
            "stop_frequency", 
            group.get_attribute<double>("center_frequency") + 
              group.get_attribute<double>("frequency_range") / 2); 
          group.set_attribute<unsigned int>("nb_frequency_points", 400); 

          if (group_name == std::string("displacement")) 
            group.set_attribute<std::string>( 
              "simulation_name", std::string("phononic_cavity_displacement")); 
          else 
            group.set_attribute<std::string>( 
              "simulation_name", std::string("phononic_cavity_calibration")); 

          group.set_attribute<bool>("save_vtu_files", false); 
        } 

      { 


        auto                    displacement = data.open_group("displacement"); 
        step62::Parameters<dim> parameters(displacement); 

        step62::ElasticWave<dim> elastic_problem(parameters); 
        elastic_problem.run(); 
      } 

      { 


        auto                    calibration = data.open_group("calibration"); 
        step62::Parameters<dim> parameters(calibration); 

        step62::ElasticWave<dim> elastic_problem(parameters); 
        elastic_problem.run(); 
      } 
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

