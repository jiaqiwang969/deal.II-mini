

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2012 - 2021 by the deal.II authors 
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
 *                        University of Siegen, 2012, 2013 
 *          Wolfgang Bangerth, Texas A&M University, 2012, 2013 
 *          Timo Heister, Texas A&M University, 2013 
 */ 



#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/index_set.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/timer.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/sparsity_tools.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/block_sparsity_pattern.h> 
#include <deal.II/lac/solver_bicgstab.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/trilinos_sparse_matrix.h> 
#include <deal.II/lac/trilinos_block_sparse_matrix.h> 
#include <deal.II/lac/trilinos_vector.h> 
#include <deal.II/lac/trilinos_parallel_block_vector.h> 
#include <deal.II/lac/trilinos_precondition.h> 
#include <deal.II/lac/trilinos_solver.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/distributed/tria.h> 
#include <deal.II/distributed/grid_refinement.h> 
#include <deal.II/distributed/solution_transfer.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 
#include <deal.II/numerics/fe_field_function.h> 

#include <fstream> 
#include <iostream> 


#include <sys/stat.h> 
#include <cerrno> 

namespace Step42 
{ 
  using namespace dealii; 



  template <int dim> 
  class ConstitutiveLaw 
  { 
  public: 
    ConstitutiveLaw(const double E, 
                    const double nu, 
                    const double sigma_0, 
                    const double gamma); 

    void set_sigma_0(double sigma_zero); 

    bool get_stress_strain_tensor( 
      const SymmetricTensor<2, dim> &strain_tensor, 
      SymmetricTensor<4, dim> &      stress_strain_tensor) const; 

    void get_linearized_stress_strain_tensors( 
      const SymmetricTensor<2, dim> &strain_tensor, 
      SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
      SymmetricTensor<4, dim> &      stress_strain_tensor) const; 

  private: 
    const double kappa; 
    const double mu; 
    double       sigma_0; 
    const double gamma; 

    const SymmetricTensor<4, dim> stress_strain_tensor_kappa; 
    const SymmetricTensor<4, dim> stress_strain_tensor_mu; 
  }; 


  template <int dim> 
  ConstitutiveLaw<dim>::ConstitutiveLaw(double E, 
                                        double nu, 
                                        double sigma_0, 
                                        double gamma) 
    : kappa(E / (3 * (1 - 2 * nu))) 
    , mu(E / (2 * (1 + nu))) 
    , sigma_0(sigma_0) 
    , gamma(gamma) 
    , stress_strain_tensor_kappa(kappa * 
                                 outer_product(unit_symmetric_tensor<dim>(), 
                                               unit_symmetric_tensor<dim>())) 
    , stress_strain_tensor_mu( 
        2 * mu * 
        (identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(), 
                                                unit_symmetric_tensor<dim>()) / 
                                    3.0)) 
  {} 

  template <int dim> 
  void ConstitutiveLaw<dim>::set_sigma_0(double sigma_zero) 
  { 
    sigma_0 = sigma_zero; 
  } 



  template <int dim> 
  bool ConstitutiveLaw<dim>::get_stress_strain_tensor( 
    const SymmetricTensor<2, dim> &strain_tensor, 
    SymmetricTensor<4, dim> &      stress_strain_tensor) const 
  { 
    Assert(dim == 3, ExcNotImplemented()); 

    SymmetricTensor<2, dim> stress_tensor; 
    stress_tensor = 
      (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 

    const SymmetricTensor<2, dim> deviator_stress_tensor = 
      deviator(stress_tensor); 
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 

    stress_strain_tensor = stress_strain_tensor_mu; 
    if (deviator_stress_tensor_norm > sigma_0) 
      { 
        const double beta = sigma_0 / deviator_stress_tensor_norm; 
        stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
      } 

    stress_strain_tensor += stress_strain_tensor_kappa; 

    return (deviator_stress_tensor_norm > sigma_0); 
  } 


  template <int dim> 
  void ConstitutiveLaw<dim>::get_linearized_stress_strain_tensors( 
    const SymmetricTensor<2, dim> &strain_tensor, 
    SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
    SymmetricTensor<4, dim> &      stress_strain_tensor) const 
  { 
    Assert(dim == 3, ExcNotImplemented()); 

    SymmetricTensor<2, dim> stress_tensor; 
    stress_tensor = 
      (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 

    stress_strain_tensor            = stress_strain_tensor_mu; 
    stress_strain_tensor_linearized = stress_strain_tensor_mu; 

    SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor); 
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 

    if (deviator_stress_tensor_norm > sigma_0) 
      { 
        const double beta = sigma_0 / deviator_stress_tensor_norm; 
        stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
        stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta); 
        deviator_stress_tensor /= deviator_stress_tensor_norm; 
        stress_strain_tensor_linearized -= 
          (1 - gamma) * beta * 2 * mu * 
          outer_product(deviator_stress_tensor, deviator_stress_tensor); 
      } 

    stress_strain_tensor += stress_strain_tensor_kappa; 
    stress_strain_tensor_linearized += stress_strain_tensor_kappa; 
  } 


  namespace EquationData 
  { 
    template <int dim> 
    class BoundaryForce : public Function<dim> 
    { 
    public: 
      BoundaryForce(); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 
    }; 

    template <int dim> 
    BoundaryForce<dim>::BoundaryForce() 
      : Function<dim>(dim) 
    {} 

    template <int dim> 
    double BoundaryForce<dim>::value(const Point<dim> &, 
                                     const unsigned int) const 
    { 
      return 0.; 
    } 

    template <int dim> 
    void BoundaryForce<dim>::vector_value(const Point<dim> &p, 
                                          Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = BoundaryForce<dim>::value(p, c); 
    } 

    template <int dim> 
    class BoundaryValues : public Function<dim> 
    { 
    public: 
      BoundaryValues(); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    BoundaryValues<dim>::BoundaryValues() 
      : Function<dim>(dim) 
    {} 

    template <int dim> 
    double BoundaryValues<dim>::value(const Point<dim> &, 
                                      const unsigned int) const 
    { 
      return 0.; 
    } 



    template <int dim> 
    class SphereObstacle : public Function<dim> 
    { 
    public: 
      SphereObstacle(const double z_surface); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 

    private: 
      const double z_surface; 
    }; 

    template <int dim> 
    SphereObstacle<dim>::SphereObstacle(const double z_surface) 
      : Function<dim>(dim) 
      , z_surface(z_surface) 
    {} 

    template <int dim> 
    double SphereObstacle<dim>::value(const Point<dim> & p, 
                                      const unsigned int component) const 
    { 
      if (component == 0) 
        return p(0); 
      else if (component == 1) 
        return p(1); 
      else if (component == 2) 
        { 
          if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5) < 0.36) 
            return (-std::sqrt(0.36 - (p(0) - 0.5) * (p(0) - 0.5) - 
                               (p(1) - 0.5) * (p(1) - 0.5)) + 
                    z_surface + 0.59); 
          else 
            return 1000; 
        } 

      Assert(false, ExcNotImplemented()); 
      return 1e9; // an unreasonable value; ignored in debug mode because of the 


    } 

    template <int dim> 
    void SphereObstacle<dim>::vector_value(const Point<dim> &p, 
                                           Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = SphereObstacle<dim>::value(p, c); 
    } 




    template <int dim> 
    class BitmapFile 
    { 
    public: 
      BitmapFile(const std::string &name); 

      double get_value(const double x, const double y) const; 

    private: 
      std::vector<double> obstacle_data; 
      double              hx, hy; 
      int                 nx, ny; 

      double get_pixel_value(const int i, const int j) const; 
    }; 


    template <int dim> 
    BitmapFile<dim>::BitmapFile(const std::string &name) 
      : obstacle_data(0) 
      , hx(0) 
      , hy(0) 
      , nx(0) 
      , ny(0) 
    { 
      std::ifstream f(name); 
      AssertThrow(f, 
                  ExcMessage(std::string("Can't read from file <") + name + 
                             ">!")); 

      std::string temp; 
      f >> temp >> nx >> ny; 

      AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format.")); 

      for (int k = 0; k < nx * ny; ++k) 
        { 
          double val; 
          f >> val; 
          obstacle_data.push_back(val); 
        } 

      hx = 1.0 / (nx - 1); 
      hy = 1.0 / (ny - 1); 

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
        std::cout << "Read obstacle from file <" << name << ">" << std::endl 
                  << "Resolution of the scanned obstacle picture: " << nx 
                  << " x " << ny << std::endl; 
    } 


    template <int dim> 
    double BitmapFile<dim>::get_pixel_value(const int i, const int j) const 
    { 
      assert(i >= 0 && i < nx); 
      assert(j >= 0 && j < ny); 
      return obstacle_data[nx * (ny - 1 - j) + i]; 
    } 

    template <int dim> 
    double BitmapFile<dim>::get_value(const double x, const double y) const 
    { 
      const int ix = std::min(std::max(static_cast<int>(x / hx), 0), nx - 2); 
      const int iy = std::min(std::max(static_cast<int>(y / hy), 0), ny - 2); 

      const double xi  = std::min(std::max((x - ix * hx) / hx, 1.), 0.); 
      const double eta = std::min(std::max((y - iy * hy) / hy, 1.), 0.); 

      return ((1 - xi) * (1 - eta) * get_pixel_value(ix, iy) + 
              xi * (1 - eta) * get_pixel_value(ix + 1, iy) + 
              (1 - xi) * eta * get_pixel_value(ix, iy + 1) + 
              xi * eta * get_pixel_value(ix + 1, iy + 1)); 
    } 


    template <int dim> 
    class ChineseObstacle : public Function<dim> 
    { 
    public: 
      ChineseObstacle(const std::string &filename, const double z_surface); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void vector_value(const Point<dim> &p, 
                                Vector<double> &  values) const override; 

    private: 
      const BitmapFile<dim> input_obstacle; 
      double                z_surface; 
    }; 

    template <int dim> 
    ChineseObstacle<dim>::ChineseObstacle(const std::string &filename, 
                                          const double       z_surface) 
      : Function<dim>(dim) 
      , input_obstacle(filename) 
      , z_surface(z_surface) 
    {} 

    template <int dim> 
    double ChineseObstacle<dim>::value(const Point<dim> & p, 
                                       const unsigned int component) const 
    { 
      if (component == 0) 
        return p(0); 
      if (component == 1) 
        return p(1); 
      else if (component == 2) 
        { 
          if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0) 
            return z_surface + 0.999 - input_obstacle.get_value(p(0), p(1)); 
        } 

      Assert(false, ExcNotImplemented()); 
      return 1e9; // an unreasonable value; ignored in debug mode because of the 


    } 

    template <int dim> 
    void ChineseObstacle<dim>::vector_value(const Point<dim> &p, 
                                            Vector<double> &  values) const 
    { 
      for (unsigned int c = 0; c < this->n_components; ++c) 
        values(c) = ChineseObstacle<dim>::value(p, c); 
    } 
  } // namespace EquationData 




  template <int dim> 
  class PlasticityContactProblem 
  { 
  public: 
    PlasticityContactProblem(const ParameterHandler &prm); 

    void run(); 

    static void declare_parameters(ParameterHandler &prm); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void compute_dirichlet_constraints(); 
    void update_solution_and_constraints(); 
    void 
         assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix); 
    void assemble_newton_system( 
      const TrilinosWrappers::MPI::Vector &linearization_point); 
    void compute_nonlinear_residual( 
      const TrilinosWrappers::MPI::Vector &linearization_point); 
    void solve_newton_system(); 
    void solve_newton(); 
    void refine_grid(); 
    void move_mesh(const TrilinosWrappers::MPI::Vector &displacement) const; 
    void output_results(const unsigned int current_refinement_cycle); 

    void output_contact_force() const; 


    MPI_Comm           mpi_communicator; 
    ConditionalOStream pcout; 
    TimerOutput        computing_timer; 



    const unsigned int                        n_initial_global_refinements; 
    parallel::distributed::Triangulation<dim> triangulation; 

    const unsigned int fe_degree; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

    AffineConstraints<double> constraints_hanging_nodes; 
    AffineConstraints<double> constraints_dirichlet_and_hanging_nodes; 
    AffineConstraints<double> all_constraints; 

    IndexSet      active_set; 
    Vector<float> fraction_of_plastic_q_points_per_cell; 


    TrilinosWrappers::SparseMatrix newton_matrix; 

    TrilinosWrappers::MPI::Vector solution; 
    TrilinosWrappers::MPI::Vector newton_rhs; 
    TrilinosWrappers::MPI::Vector newton_rhs_uncondensed; 
    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector; 


    const double         e_modulus, nu, gamma, sigma_0; 
    ConstitutiveLaw<dim> constitutive_law; 


    const std::string                          base_mesh; 
    const std::shared_ptr<const Function<dim>> obstacle; 

    struct RefinementStrategy 
    { 
      enum value 
      { 
        refine_global, 
        refine_percentage, 
        refine_fix_dofs 
      }; 
    }; 
    typename RefinementStrategy::value refinement_strategy; 

    const bool         transfer_solution; 
    std::string        output_dir; 
    const unsigned int n_refinement_cycles; 
    unsigned int       current_refinement_cycle; 
  }; 


  template <int dim> 
  void PlasticityContactProblem<dim>::declare_parameters(ParameterHandler &prm) 
  { 
    prm.declare_entry( 
      "polynomial degree", 
      "1", 
      Patterns::Integer(), 
      "Polynomial degree of the FE_Q finite element space, typically 1 or 2."); 
    prm.declare_entry("number of initial refinements", 
                      "2", 
                      Patterns::Integer(), 
                      "Number of initial global mesh refinement steps before " 
                      "the first computation."); 
    prm.declare_entry( 
      "refinement strategy", 
      "percentage", 
      Patterns::Selection("global|percentage"), 
      "Mesh refinement strategy:\n" 
      " global: one global refinement\n" 
      " percentage: a fixed percentage of cells gets refined using the Kelly estimator."); 
    prm.declare_entry("number of cycles", 
                      "5", 
                      Patterns::Integer(), 
                      "Number of adaptive mesh refinement cycles to run."); 
    prm.declare_entry( 
      "obstacle", 
      "sphere", 
      Patterns::Selection("sphere|read from file"), 
      "The name of the obstacle to use. This may either be 'sphere' if we should " 
      "use a spherical obstacle, or 'read from file' in which case the obstacle " 
      "will be read from a file named 'obstacle.pbm' that is supposed to be in " 
      "ASCII PBM format."); 
    prm.declare_entry( 
      "output directory", 
      "", 
      Patterns::Anything(), 
      "Directory for output files (graphical output and benchmark " 
      "statistics). If empty, use the current directory."); 
    prm.declare_entry( 
      "transfer solution", 
      "false", 
      Patterns::Bool(), 
      "Whether the solution should be used as a starting guess " 
      "for the next finer mesh. If false, then the iteration starts at " 
      "zero on every mesh."); 
    prm.declare_entry("base mesh", 
                      "box", 
                      Patterns::Selection("box|half sphere"), 
                      "Select the shape of the domain: 'box' or 'half sphere'"); 
  } 


  template <int dim> 
  PlasticityContactProblem<dim>::PlasticityContactProblem( 
    const ParameterHandler &prm) 
    : mpi_communicator(MPI_COMM_WORLD) 
    , pcout(std::cout, 
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
    , computing_timer(MPI_COMM_WORLD, 
                      pcout, 
                      TimerOutput::never, 
                      TimerOutput::wall_times) 

    , n_initial_global_refinements( 
        prm.get_integer("number of initial refinements")) 
    , triangulation(mpi_communicator) 
    , fe_degree(prm.get_integer("polynomial degree")) 
    , fe(FE_Q<dim>(QGaussLobatto<1>(fe_degree + 1)), dim) 
    , dof_handler(triangulation) 

    , e_modulus(200000) 
    , nu(0.3) 
    , gamma(0.01) 
    , sigma_0(400.0) 
    , constitutive_law(e_modulus, nu, sigma_0, gamma) 

    , base_mesh(prm.get("base mesh")) 
    , obstacle(prm.get("obstacle") == "read from file" ? 
                 static_cast<const Function<dim> *>( 
                   new EquationData::ChineseObstacle<dim>( 
                     "obstacle.pbm", 
                     (base_mesh == "box" ? 1.0 : 0.5))) : 
                 static_cast<const Function<dim> *>( 
                   new EquationData::SphereObstacle<dim>( 
                     base_mesh == "box" ? 1.0 : 0.5))) 

    , transfer_solution(prm.get_bool("transfer solution")) 
    , n_refinement_cycles(prm.get_integer("number of cycles")) 
    , current_refinement_cycle(0) 

  { 
    std::string strat = prm.get("refinement strategy"); 
    if (strat == "global") 
      refinement_strategy = RefinementStrategy::refine_global; 
    else if (strat == "percentage") 
      refinement_strategy = RefinementStrategy::refine_percentage; 
    else 
      AssertThrow(false, ExcNotImplemented()); 

    output_dir = prm.get("output directory"); 
    if (output_dir != "" && *(output_dir.rbegin()) != '/') 
      output_dir += "/"; 


    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
      { 
        const int ierr = mkdir(output_dir.c_str(), 0777); 
        AssertThrow(ierr == 0 || errno == EEXIST, ExcIO()); 
      } 

    pcout << "    Using output directory '" << output_dir << "'" << std::endl; 
    pcout << "    FE degree " << fe_degree << std::endl; 
    pcout << "    transfer solution " << (transfer_solution ? "true" : "false") 
          << std::endl; 
  } 




  Point<3> rotate_half_sphere(const Point<3> &in) 
  { 
    return {in(2), in(1), -in(0)}; 
  } 

  template <int dim> 
  void PlasticityContactProblem<dim>::make_grid() 
  { 
    if (base_mesh == "half sphere") 
      { 
        const Point<dim> center(0, 0, 0); 
        const double     radius = 0.8; 
        GridGenerator::half_hyper_ball(triangulation, center, radius); 


        triangulation.reset_all_manifolds(); 

        GridTools::transform(&rotate_half_sphere, triangulation); 
        GridTools::shift(Point<dim>(0.5, 0.5, 0.5), triangulation); 

        SphericalManifold<dim> manifold_description(Point<dim>(0.5, 0.5, 0.5)); 
        GridTools::copy_boundary_to_manifold_id(triangulation); 
        triangulation.set_manifold(0, manifold_description); 
      } 


    else 
      { 
        const Point<dim> p1(0, 0, 0); 
        const Point<dim> p2(1.0, 1.0, 1.0); 

        GridGenerator::hyper_rectangle(triangulation, p1, p2); 

        for (const auto &cell : triangulation.active_cell_iterators()) 
          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary()) 
              { 
                if (std::fabs(face->center()[2] - p2[2]) < 1e-12) 
                  face->set_boundary_id(1); 
                if (std::fabs(face->center()[0] - p1[0]) < 1e-12 || 
                    std::fabs(face->center()[0] - p2[0]) < 1e-12 || 
                    std::fabs(face->center()[1] - p1[1]) < 1e-12 || 
                    std::fabs(face->center()[1] - p2[1]) < 1e-12) 
                  face->set_boundary_id(8); 
                if (std::fabs(face->center()[2] - p1[2]) < 1e-12) 
                  face->set_boundary_id(6); 
              } 
      } 

    triangulation.refine_global(n_initial_global_refinements); 
  } 




  template <int dim> 
  void PlasticityContactProblem<dim>::setup_system() 
  { 

/* 设置dofs，并为本地拥有的相关dofs获取索引集  */ 

    { 
      TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs"); 
      dof_handler.distribute_dofs(fe); 

      locally_owned_dofs = dof_handler.locally_owned_dofs(); 
      locally_relevant_dofs.clear(); 
      DoFTools::extract_locally_relevant_dofs(dof_handler, 
                                              locally_relevant_dofs); 
    } 

/*设置悬挂节点和Dirichlet约束 */ 

 
    { 
      TimerOutput::Scope t(computing_timer, "Setup: constraints"); 
      constraints_hanging_nodes.reinit(locally_relevant_dofs); 
      DoFTools::make_hanging_node_constraints(dof_handler, 
                                              constraints_hanging_nodes); 
      constraints_hanging_nodes.close(); 

      pcout << "   Number of active cells: " 
            << triangulation.n_global_active_cells() << std::endl 
            << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
            << std::endl; 

      compute_dirichlet_constraints(); 
    } 

/* 初始化向量和活动集  */ 

    { 
      TimerOutput::Scope t(computing_timer, "Setup: vectors"); 
      solution.reinit(locally_relevant_dofs, mpi_communicator); 
      newton_rhs.reinit(locally_owned_dofs, mpi_communicator); 
      newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator); 
      diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator); 
      fraction_of_plastic_q_points_per_cell.reinit( 
        triangulation.n_active_cells()); 

      active_set.clear(); 
      active_set.set_size(dof_handler.n_dofs()); 
    } 


    { 
      TimerOutput::Scope                t(computing_timer, "Setup: matrix"); 
      TrilinosWrappers::SparsityPattern sp(locally_owned_dofs, 
                                           mpi_communicator); 

      DoFTools::make_sparsity_pattern(dof_handler, 
                                      sp, 
                                      constraints_dirichlet_and_hanging_nodes, 
                                      false, 
                                      Utilities::MPI::this_mpi_process( 
                                        mpi_communicator)); 
      sp.compress(); 
      newton_matrix.reinit(sp); 

      TrilinosWrappers::SparseMatrix &mass_matrix = newton_matrix; 

      assemble_mass_matrix_diagonal(mass_matrix); 

      const unsigned int start = (newton_rhs.local_range().first), 
                         end   = (newton_rhs.local_range().second); 
      for (unsigned int j = start; j < end; ++j) 
        diag_mass_matrix_vector(j) = mass_matrix.diag_element(j); 
      diag_mass_matrix_vector.compress(VectorOperation::insert); 

      mass_matrix = 0; 
    } 
  } 





  template <int dim> 
  void PlasticityContactProblem<dim>::compute_dirichlet_constraints() 
  { 
    constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs); 
    constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes); 

    if (base_mesh == "box") 
      { 


        VectorTools::interpolate_boundary_values( 
          dof_handler, 
          6, 
          EquationData::BoundaryValues<dim>(), 
          constraints_dirichlet_and_hanging_nodes, 
          ComponentMask()); 


        const FEValuesExtractors::Scalar x_displacement(0); 
        const FEValuesExtractors::Scalar y_displacement(1); 
        VectorTools::interpolate_boundary_values( 
          dof_handler, 
          8, 
          EquationData::BoundaryValues<dim>(), 
          constraints_dirichlet_and_hanging_nodes, 
          (fe.component_mask(x_displacement) | 
           fe.component_mask(y_displacement))); 
      } 
    else 
      VectorTools::interpolate_boundary_values( 
        dof_handler, 
        0, 
        EquationData::BoundaryValues<dim>(), 
        constraints_dirichlet_and_hanging_nodes, 
        ComponentMask()); 

    constraints_dirichlet_and_hanging_nodes.close(); 
  } 



  template <int dim> 
  void PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal( 
    TrilinosWrappers::SparseMatrix &mass_matrix) 
  { 
    QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 
              cell_matrix = 0; 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                  cell_matrix(i, i) += 
                    (fe_values_face[displacement].value(i, q_point) * 
                     fe_values_face[displacement].value(i, q_point) * 
                     fe_values_face.JxW(q_point)); 

              cell->get_dof_indices(local_dof_indices); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                mass_matrix.add(local_dof_indices[i], 
                                local_dof_indices[i], 
                                cell_matrix(i, i)); 
            } 
    mass_matrix.compress(VectorOperation::add); 
  } 



  template <int dim> 
  void PlasticityContactProblem<dim>::update_solution_and_constraints() 
  { 
    std::vector<bool> dof_touched(dof_handler.n_dofs(), false); 

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 
    distributed_solution = solution; 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = newton_rhs_uncondensed; 

    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant( 
      locally_relevant_dofs, mpi_communicator); 
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector; 

    all_constraints.reinit(locally_relevant_dofs); 
    active_set.clear(); 


    Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points()); 
    FEFaceValues<dim>   fe_values_face(fe, 
                                     face_quadrature, 
                                     update_quadrature_points); 

    const unsigned int dofs_per_face   = fe.n_dofs_per_face(); 
    const unsigned int n_face_q_points = face_quadrature.size(); 

    std::vector<types::global_dof_index> dof_indices(dofs_per_face); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (!cell->is_artificial()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 
              face->get_dof_indices(dof_indices); 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                { 


                  const unsigned int component = 
                    fe.face_system_to_component_index(q_point).first; 

                  const unsigned int index_z = dof_indices[q_point]; 

                  if ((component == 2) && (dof_touched[index_z] == false)) 
                    { 
                      dof_touched[index_z] = true; 

                      const Point<dim> this_support_point = 
                        fe_values_face.quadrature_point(q_point); 

                      const double obstacle_value = 
                        obstacle->value(this_support_point, 2); 
                      const double solution_here = solution(index_z); 
                      const double undeformed_gap = 
                        obstacle_value - this_support_point(2); 

                      const double c = 100.0 * e_modulus; 
                      if ((lambda(index_z) / 
                               diag_mass_matrix_vector_relevant(index_z) + 
                             c * (solution_here - undeformed_gap) > 
                           0) && 
                          !constraints_hanging_nodes.is_constrained(index_z)) 
                        { 
                          all_constraints.add_line(index_z); 
                          all_constraints.set_inhomogeneity(index_z, 
                                                            undeformed_gap); 
                          distributed_solution(index_z) = undeformed_gap; 

                          active_set.add_index(index_z); 
                        } 
                    } 
                } 
            } 


    distributed_solution.compress(VectorOperation::insert); 
    solution = distributed_solution; 

    all_constraints.close(); 
    all_constraints.merge(constraints_dirichlet_and_hanging_nodes); 

 
          << Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(), 
                                 mpi_communicator) 
          << std::endl; 
  } 


  template <int dim> 
  void PlasticityContactProblem<dim>::assemble_newton_system( 
    const TrilinosWrappers::MPI::Vector &linearization_point) 
  { 
    TimerOutput::Scope t(computing_timer, "Assembling"); 

    QGauss<dim>     quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const EquationData::BoundaryForce<dim> boundary_force; 
    std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
                                                      Vector<double>(dim)); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          cell_matrix = 0; 
          cell_rhs    = 0; 

          std::vector<SymmetricTensor<2, dim>> strain_tensor(n_q_points); 
          fe_values[displacement].get_function_symmetric_gradients( 
            linearization_point, strain_tensor); 

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              SymmetricTensor<4, dim> stress_strain_tensor_linearized; 
              SymmetricTensor<4, dim> stress_strain_tensor; 
              constitutive_law.get_linearized_stress_strain_tensors( 
                strain_tensor[q_point], 
                stress_strain_tensor_linearized, 
                stress_strain_tensor); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 


                  const SymmetricTensor<2, dim> stress_phi_i = 
                    stress_strain_tensor_linearized * 
                    fe_values[displacement].symmetric_gradient(i, q_point); 

                  for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                    cell_matrix(i, j) += 
                      (stress_phi_i * 
                       fe_values[displacement].symmetric_gradient(j, q_point) * 
                       fe_values.JxW(q_point)); 

                  cell_rhs(i) += 
                    ((stress_phi_i - 
                      stress_strain_tensor * 
                        fe_values[displacement].symmetric_gradient(i, 
                                                                   q_point)) * 
                     strain_tensor[q_point] * fe_values.JxW(q_point)); 
                } 
            } 

          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary() && face->boundary_id() == 1) 
              { 
                fe_values_face.reinit(cell, face); 

                boundary_force.vector_value_list( 
                  fe_values_face.get_quadrature_points(), 
                  boundary_force_values); 

                for (unsigned int q_point = 0; q_point < n_face_q_points; 
                     ++q_point) 
                  { 
                    Tensor<1, dim> rhs_values; 
                    rhs_values[2] = boundary_force_values[q_point][2]; 
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      cell_rhs(i) += 
                        (fe_values_face[displacement].value(i, q_point) * 
                         rhs_values * fe_values_face.JxW(q_point)); 
                  } 
              } 

          cell->get_dof_indices(local_dof_indices); 
          all_constraints.distribute_local_to_global(cell_matrix, 
                                                     cell_rhs, 
                                                     local_dof_indices, 
                                                     newton_matrix, 
                                                     newton_rhs, 
                                                     true); 
        } 

    newton_matrix.compress(VectorOperation::add); 
    newton_rhs.compress(VectorOperation::add); 
  } 





  template <int dim> 
  void PlasticityContactProblem<dim>::compute_nonlinear_residual( 
    const TrilinosWrappers::MPI::Vector &linearization_point) 
  { 
    QGauss<dim>     quadrature_formula(fe.degree + 1); 
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_JxW_values); 

    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_quadrature_points | 
                                       update_JxW_values); 

    const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points      = quadrature_formula.size(); 
    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const EquationData::BoundaryForce<dim> boundary_force; 
    std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
                                                      Vector<double>(dim)); 

    Vector<double> cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const FEValuesExtractors::Vector displacement(0); 

    newton_rhs             = 0; 
    newton_rhs_uncondensed = 0; 

    fraction_of_plastic_q_points_per_cell = 0; 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          fe_values.reinit(cell); 
          cell_rhs = 0; 

          std::vector<SymmetricTensor<2, dim>> strain_tensors(n_q_points); 
          fe_values[displacement].get_function_symmetric_gradients( 
            linearization_point, strain_tensors); 

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              SymmetricTensor<4, dim> stress_strain_tensor; 
              const bool              q_point_is_plastic = 
                constitutive_law.get_stress_strain_tensor( 
                  strain_tensors[q_point], stress_strain_tensor); 
              if (q_point_is_plastic) 
                ++fraction_of_plastic_q_points_per_cell( 
                  cell->active_cell_index()); 

              for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                { 
                  cell_rhs(i) -= 
                    (strain_tensors[q_point] * stress_strain_tensor * 
                     fe_values[displacement].symmetric_gradient(i, q_point) * 
                     fe_values.JxW(q_point)); 

                  Tensor<1, dim> rhs_values; 
                  rhs_values = 0; 
                  cell_rhs(i) += (fe_values[displacement].value(i, q_point) * 
                                  rhs_values * fe_values.JxW(q_point)); 
                } 
            } 

          for (const auto &face : cell->face_iterators()) 
            if (face->at_boundary() && face->boundary_id() == 1) 
              { 
                fe_values_face.reinit(cell, face); 

                boundary_force.vector_value_list( 
                  fe_values_face.get_quadrature_points(), 
                  boundary_force_values); 

                for (unsigned int q_point = 0; q_point < n_face_q_points; 
                     ++q_point) 
                  { 
                    Tensor<1, dim> rhs_values; 
                    rhs_values[2] = boundary_force_values[q_point][2]; 
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                      cell_rhs(i) += 
                        (fe_values_face[displacement].value(i, q_point) * 
                         rhs_values * fe_values_face.JxW(q_point)); 
                  } 
              } 

          cell->get_dof_indices(local_dof_indices); 
          constraints_dirichlet_and_hanging_nodes.distribute_local_to_global( 
            cell_rhs, local_dof_indices, newton_rhs); 

          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            newton_rhs_uncondensed(local_dof_indices[i]) += cell_rhs(i); 
        } 

    fraction_of_plastic_q_points_per_cell /= quadrature_formula.size(); 
    newton_rhs.compress(VectorOperation::add); 
    newton_rhs_uncondensed.compress(VectorOperation::add); 
  } 








  template <int dim> 
  void PlasticityContactProblem<dim>::solve_newton_system() 
  { 
    TimerOutput::Scope t(computing_timer, "Solve"); 

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 
    distributed_solution = solution; 

    constraints_hanging_nodes.set_zero(distributed_solution); 
    constraints_hanging_nodes.set_zero(newton_rhs); 

 
    { 
      TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner"); 

      std::vector<std::vector<bool>> constant_modes; 
      DoFTools::extract_constant_modes(dof_handler, 
                                       ComponentMask(), 
                                       constant_modes); 

      TrilinosWrappers::PreconditionAMG::AdditionalData additional_data; 
      additional_data.constant_modes        = constant_modes; 
      additional_data.elliptic              = true; 
      additional_data.n_cycles              = 1; 
      additional_data.w_cycle               = false; 
      additional_data.output_details        = false; 
      additional_data.smoother_sweeps       = 2; 
      additional_data.aggregation_threshold = 1e-2; 

      preconditioner.initialize(newton_matrix, additional_data); 
    } 

    { 
      TimerOutput::Scope t(computing_timer, "Solve: iterate"); 

      TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 

      const double relative_accuracy = 1e-8; 
      const double solver_tolerance = 
        relative_accuracy * 
        newton_matrix.residual(tmp, distributed_solution, newton_rhs); 

      SolverControl solver_control(newton_matrix.m(), solver_tolerance); 
      SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control); 
      solver.solve(newton_matrix, 
                   distributed_solution, 
                   newton_rhs, 
                   preconditioner); 

      pcout << "         Error: " << solver_control.initial_value() << " -> " 
            << solver_control.last_value() << " in " 
            << solver_control.last_step() << " Bicgstab iterations." 
            << std::endl; 
    } 

    all_constraints.distribute(distributed_solution); 

    solution = distributed_solution; 
  } 



  template <int dim> 
  void PlasticityContactProblem<dim>::solve_newton() 
  { 
    TrilinosWrappers::MPI::Vector old_solution(locally_owned_dofs, 
                                               mpi_communicator); 
    TrilinosWrappers::MPI::Vector residual(locally_owned_dofs, 
                                           mpi_communicator); 
    TrilinosWrappers::MPI::Vector tmp_vector(locally_owned_dofs, 
                                             mpi_communicator); 
    TrilinosWrappers::MPI::Vector locally_relevant_tmp_vector( 
      locally_relevant_dofs, mpi_communicator); 
    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                       mpi_communicator); 

    double residual_norm; 
    double previous_residual_norm = -std::numeric_limits<double>::max(); 

    const double correct_sigma = sigma_0; 

    IndexSet old_active_set(active_set); 

    for (unsigned int newton_step = 1; newton_step <= 100; ++newton_step) 
      { 
        if (newton_step == 1 && 
            ((transfer_solution && current_refinement_cycle == 0) || 
             !transfer_solution)) 
          constitutive_law.set_sigma_0(1e+10); 
        else if (newton_step == 2 || current_refinement_cycle > 0 || 
                 !transfer_solution) 
          constitutive_law.set_sigma_0(correct_sigma); 

        pcout << " " << std::endl; 
        pcout << "   Newton iteration " << newton_step << std::endl; 
        pcout << "      Updating active set..." << std::endl; 

        { 
          TimerOutput::Scope t(computing_timer, "update active set"); 
          update_solution_and_constraints(); 
        } 

        pcout << "      Assembling system... " << std::endl; 
        newton_matrix = 0; 
        newton_rhs    = 0; 
        assemble_newton_system(solution); 

        pcout << "      Solving system... " << std::endl; 
        solve_newton_system(); 




        if ((newton_step == 1) || 
            (transfer_solution && newton_step == 2 && 
             current_refinement_cycle == 0) || 
            (!transfer_solution && newton_step == 2)) 
          { 
            compute_nonlinear_residual(solution); 
            old_solution = solution; 

            residual                     = newton_rhs; 
            const unsigned int start_res = (residual.local_range().first), 
                               end_res   = (residual.local_range().second); 
            for (unsigned int n = start_res; n < end_res; ++n) 
              if (all_constraints.is_inhomogeneously_constrained(n)) 
                residual(n) = 0; 

            residual.compress(VectorOperation::insert); 

            residual_norm = residual.l2_norm(); 

            pcout << "      Accepting Newton solution with residual: " 
                  << residual_norm << std::endl; 
          } 
        else 
          { 
            for (unsigned int i = 0; i < 5; ++i) 
              { 
                distributed_solution = solution; 

                const double alpha = std::pow(0.5, static_cast<double>(i)); 
                tmp_vector         = old_solution; 
                tmp_vector.sadd(1 - alpha, alpha, distributed_solution); 

                TimerOutput::Scope t(computing_timer, "Residual and lambda"); 

                locally_relevant_tmp_vector = tmp_vector; 
                compute_nonlinear_residual(locally_relevant_tmp_vector); 
                residual = newton_rhs; 

                const unsigned int start_res = (residual.local_range().first), 
                                   end_res   = (residual.local_range().second); 
                for (unsigned int n = start_res; n < end_res; ++n) 
                  if (all_constraints.is_inhomogeneously_constrained(n)) 
                    residual(n) = 0; 

                residual.compress(VectorOperation::insert); 

                residual_norm = residual.l2_norm(); 

 
                  << "      Residual of the non-contact part of the system: " 
                  << residual_norm << std::endl 
                  << "         with a damping parameter alpha = " << alpha 
                  << std::endl; 

                if (residual_norm < previous_residual_norm) 
                  break; 
              } 

            solution     = tmp_vector; 
            old_solution = solution; 
          } 

        previous_residual_norm = residual_norm; 


        if (Utilities::MPI::sum((active_set == old_active_set) ? 0 : 1, 
                                mpi_communicator) == 0) 
          { 
            pcout << "      Active set did not change!" << std::endl; 
            if (residual_norm < 1e-10) 
              break; 
          } 

        old_active_set = active_set; 
      } 
  } 


  template <int dim> 
  void PlasticityContactProblem<dim>::refine_grid() 
  { 
    if (refinement_strategy == RefinementStrategy::refine_global) 
      { 
        for (typename Triangulation<dim>::active_cell_iterator cell = 
               triangulation.begin_active(); 
             cell != triangulation.end(); 
             ++cell) 
          if (cell->is_locally_owned()) 
            cell->set_refine_flag(); 
      } 
    else 
      { 
        Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
        KellyErrorEstimator<dim>::estimate( 
          dof_handler, 
          QGauss<dim - 1>(fe.degree + 2), 
          std::map<types::boundary_id, const Function<dim> *>(), 
          solution, 
          estimated_error_per_cell); 

        parallel::distributed::GridRefinement ::refine_and_coarsen_fixed_number( 
          triangulation, estimated_error_per_cell, 0.3, 0.03); 
      } 

    triangulation.prepare_coarsening_and_refinement(); 

    parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> 
      solution_transfer(dof_handler); 
    if (transfer_solution) 
      solution_transfer.prepare_for_coarsening_and_refinement(solution); 

    triangulation.execute_coarsening_and_refinement(); 

    setup_system(); 

    if (transfer_solution) 
      { 
        TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
                                                           mpi_communicator); 
        solution_transfer.interpolate(distributed_solution); 


        constraints_hanging_nodes.distribute(distributed_solution); 

        solution = distributed_solution; 
        compute_nonlinear_residual(solution); 
      } 
  } 



  template <int dim> 
  void PlasticityContactProblem<dim>::move_mesh( 
    const TrilinosWrappers::MPI::Vector &displacement) const 
  { 
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto v : cell->vertex_indices()) 
          if (vertex_touched[cell->vertex_index(v)] == false) 
            { 
              vertex_touched[cell->vertex_index(v)] = true; 

              Point<dim> vertex_displacement; 
              for (unsigned int d = 0; d < dim; ++d) 
                vertex_displacement[d] = 
                  displacement(cell->vertex_dof_index(v, d)); 

              cell->vertex(v) += vertex_displacement; 
            } 
  } 



  template <int dim> 
  void PlasticityContactProblem<dim>::output_results( 
    const unsigned int current_refinement_cycle) 
  { 
    TimerOutput::Scope t(computing_timer, "Graphical output"); 

    pcout << "      Writing graphical output... " << std::flush; 

    move_mesh(solution); 


    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
                                                     mpi_communicator); 
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
                       end_res = (newton_rhs_uncondensed.local_range().second); 
    for (unsigned int n = start_res; n < end_res; ++n) 
      if (all_constraints.is_inhomogeneously_constrained(n)) 
        distributed_lambda(n) = 
          newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
    distributed_lambda.compress(VectorOperation::insert); 
    constraints_hanging_nodes.distribute(distributed_lambda); 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = distributed_lambda; 

    TrilinosWrappers::MPI::Vector distributed_active_set_vector( 
      locally_owned_dofs, mpi_communicator); 
    distributed_active_set_vector = 0.; 
    for (const auto index : active_set) 
      distributed_active_set_vector[index] = 1.; 
    distributed_lambda.compress(VectorOperation::insert); 

    TrilinosWrappers::MPI::Vector active_set_vector(locally_relevant_dofs, 
                                                    mpi_communicator); 
    active_set_vector = distributed_active_set_vector; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 

    const std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_out.add_data_vector(solution, 
                             std::vector<std::string>(dim, "displacement"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.add_data_vector(lambda, 
                             std::vector<std::string>(dim, "contact_force"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.add_data_vector(active_set_vector, 
                             std::vector<std::string>(dim, "active_set"), 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 

    Vector<float> subdomain(triangulation.n_active_cells()); 
    for (unsigned int i = 0; i < subdomain.size(); ++i) 
      subdomain(i) = triangulation.locally_owned_subdomain(); 
    data_out.add_data_vector(subdomain, "subdomain"); 

    data_out.add_data_vector(fraction_of_plastic_q_points_per_cell, 
                             "fraction_of_plastic_q_points"); 

    data_out.build_patches(); 


    const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
      output_dir, "solution", current_refinement_cycle, mpi_communicator, 2); 
    pcout << pvtu_filename << std::endl; 

    TrilinosWrappers::MPI::Vector tmp(solution); 
    tmp *= -1; 
    move_mesh(tmp); 
  } 


  template <int dim> 
  void PlasticityContactProblem<dim>::output_contact_force() const 
  { 
    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
                                                     mpi_communicator); 
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
                       end_res = (newton_rhs_uncondensed.local_range().second); 
    for (unsigned int n = start_res; n < end_res; ++n) 
      if (all_constraints.is_inhomogeneously_constrained(n)) 
        distributed_lambda(n) = 
          newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
      else 
        distributed_lambda(n) = 0; 
    distributed_lambda.compress(VectorOperation::insert); 
    constraints_hanging_nodes.distribute(distributed_lambda); 

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
                                         mpi_communicator); 
    lambda = distributed_lambda; 

    double contact_force = 0.0; 

    QGauss<dim - 1>   face_quadrature_formula(fe.degree + 1); 
    FEFaceValues<dim> fe_values_face(fe, 
                                     face_quadrature_formula, 
                                     update_values | update_JxW_values); 

    const unsigned int n_face_q_points = face_quadrature_formula.size(); 

    const FEValuesExtractors::Vector displacement(0); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        for (const auto &face : cell->face_iterators()) 
          if (face->at_boundary() && face->boundary_id() == 1) 
            { 
              fe_values_face.reinit(cell, face); 

              std::vector<Tensor<1, dim>> lambda_values(n_face_q_points); 
              fe_values_face[displacement].get_function_values(lambda, 
                                                               lambda_values); 

              for (unsigned int q_point = 0; q_point < n_face_q_points; 
                   ++q_point) 
                contact_force += 
                  lambda_values[q_point][2] * fe_values_face.JxW(q_point); 
            } 
    contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD); 

    pcout << "Contact force = " << contact_force << std::endl; 
  } 


  template <int dim> 
  void PlasticityContactProblem<dim>::run() 
  { 
    computing_timer.reset(); 
    for (; current_refinement_cycle < n_refinement_cycles; 
         ++current_refinement_cycle) 
      { 
        { 
          TimerOutput::Scope t(computing_timer, "Setup"); 

          pcout << std::endl; 
          pcout << "Cycle " << current_refinement_cycle << ':' << std::endl; 

          if (current_refinement_cycle == 0) 
            { 
              make_grid(); 
              setup_system(); 
            } 
          else 
            { 
              TimerOutput::Scope t(computing_timer, "Setup: refine mesh"); 
              refine_grid(); 
            } 
        } 

        solve_newton(); 

        output_results(current_refinement_cycle); 

        computing_timer.print_summary(); 
        computing_timer.reset(); 

        Utilities::System::MemoryStats stats; 
        Utilities::System::get_memory_stats(stats); 
        pcout << "Peak virtual memory used, resident in kB: " << stats.VmSize 
              << " " << stats.VmRSS << std::endl; 

        if (base_mesh == "box") 
          output_contact_force(); 
      } 
  } 
} // namespace Step42 


int main(int argc, char *argv[]) 
{ 
  using namespace dealii; 
  using namespace Step42; 

  try 
    { 
      ParameterHandler prm; 
      PlasticityContactProblem<3>::declare_parameters(prm); 
      if (argc != 2) 
        { 
          std::cerr << "*** Call this program as <./step-42 input.prm>" 
                    << std::endl; 
          return 1; 
        } 

      prm.parse_input(argv[1]); 
      Utilities::MPI::MPI_InitFinalize mpi_initialization( 
        argc, argv, numbers::invalid_unsigned_int); 
      { 
        PlasticityContactProblem<3> problem(prm); 
        problem.run(); 
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


