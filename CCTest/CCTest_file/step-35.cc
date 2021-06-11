

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
 * Author: Abner Salgado, Texas A&M University 2009 
 */ 


// @sect3{Include files}  

// 我们首先包括所有必要的deal.II头文件和一些C++相关的文件。它们中的每一个都已经在以前的教程程序中讨论过了，所以我们在这里就不做详细介绍了。

#include <deal.II/base/parameter_handler.h> 
#include <deal.II/base/point.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/multithread_info.h> 
#include <deal.II/base/thread_management.h> 
#include <deal.II/base/work_stream.h> 
#include <deal.II/base/parallel.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/conditional_ostream.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/sparse_ilu.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/grid_in.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/dofs/dof_renumbering.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_tools.h> 
#include <deal.II/fe/fe_system.h> 

#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 

#include <fstream> 
#include <cmath> 
#include <iostream> 

// 最后这和以前的所有程序一样。

namespace Step35 
{ 
  using namespace dealii; 
// @sect3{Run time parameters}  

// 由于我们的方法有几个可以微调的参数，我们把它们放到一个外部文件中，这样就可以在运行时确定它们。

// 这尤其包括辅助变量  $\phi$  的方程表述，为此我们声明一个  <code>enum</code>  。接下来，我们声明一个类，它将读取和存储我们的程序运行所需的所有参数。

  namespace RunTimeParameters 
  { 
    enum class Method 
    { 
      standard, 
      rotational 
    }; 

    class Data_Storage 
    { 
    public: 
      Data_Storage(); 

      void read_data(const std::string &filename); 

      Method form; 

      double dt; 
      double initial_time; 
      double final_time; 

      double Reynolds; 

      unsigned int n_global_refines; 

      unsigned int pressure_degree; 

      unsigned int vel_max_iterations; 
      unsigned int vel_Krylov_size; 
      unsigned int vel_off_diagonals; 
      unsigned int vel_update_prec; 
      double       vel_eps; 
      double       vel_diag_strength; 

      bool         verbose; 
      unsigned int output_interval; 

    protected: 
      ParameterHandler prm; 
    }; 

// 在这个类的构造函数中，我们声明所有的参数。这方面的细节已经在其他地方讨论过了，例如在  step-29  。

    Data_Storage::Data_Storage() 
      : form(Method::rotational) 
      , dt(5e-4) 
      , initial_time(0.) 
      , final_time(1.) 
      , Reynolds(1.) 
      , n_global_refines(0) 
      , pressure_degree(1) 
      , vel_max_iterations(1000) 
      , vel_Krylov_size(30) 
      , vel_off_diagonals(60) 
      , vel_update_prec(15) 
      , vel_eps(1e-12) 
      , vel_diag_strength(0.01) 
      , verbose(true) 
      , output_interval(15) 
    { 
      prm.declare_entry("Method_Form", 
                        "rotational", 
                        Patterns::Selection("rotational|standard"), 
                        " Used to select the type of method that we are going " 
                        "to use. "); 
      prm.enter_subsection("Physical data"); 
      { 
        prm.declare_entry("initial_time", 
                          "0.", 
                          Patterns::Double(0.), 
                          " The initial time of the simulation. "); 
        prm.declare_entry("final_time", 
                          "1.", 
                          Patterns::Double(0.), 
                          " The final time of the simulation. "); 
        prm.declare_entry("Reynolds", 
                          "1.", 
                          Patterns::Double(0.), 
                          " The Reynolds number. "); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Time step data"); 
      { 
        prm.declare_entry("dt", 
                          "5e-4", 
                          Patterns::Double(0.), 
                          " The time step size. "); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Space discretization"); 
      { 
        prm.declare_entry("n_of_refines", 
                          "0", 
                          Patterns::Integer(0, 15), 
                          " The number of global refines we do on the mesh. "); 
        prm.declare_entry("pressure_fe_degree", 
                          "1", 
                          Patterns::Integer(1, 5), 
                          " The polynomial degree for the pressure space. "); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Data solve velocity"); 
      { 
        prm.declare_entry( 
          "max_iterations", 
          "1000", 
          Patterns::Integer(1, 1000), 
          " The maximal number of iterations GMRES must make. "); 
        prm.declare_entry("eps", 
                          "1e-12", 
                          Patterns::Double(0.), 
                          " The stopping criterion. "); 
        prm.declare_entry("Krylov_size", 
                          "30", 
                          Patterns::Integer(1), 
                          " The size of the Krylov subspace to be used. "); 
        prm.declare_entry("off_diagonals", 
                          "60", 
                          Patterns::Integer(0), 
                          " The number of off-diagonal elements ILU must " 
                          "compute. "); 
        prm.declare_entry("diag_strength", 
                          "0.01", 
                          Patterns::Double(0.), 
                          " Diagonal strengthening coefficient. "); 
        prm.declare_entry("update_prec", 
                          "15", 
                          Patterns::Integer(1), 
                          " This number indicates how often we need to " 
                          "update the preconditioner"); 
      } 
      prm.leave_subsection(); 

      prm.declare_entry("verbose", 
                        "true", 
                        Patterns::Bool(), 
                        " This indicates whether the output of the solution " 
                        "process should be verbose. "); 

      prm.declare_entry("output_interval", 
                        "1", 
                        Patterns::Integer(1), 
                        " This indicates between how many time steps we print " 
                        "the solution. "); 
    } 

    void Data_Storage::read_data(const std::string &filename) 
    { 
      std::ifstream file(filename); 
      AssertThrow(file, ExcFileNotOpen(filename)); 

      prm.parse_input(file); 

      if (prm.get("Method_Form") == std::string("rotational")) 
        form = Method::rotational; 
      else 
        form = Method::standard; 

      prm.enter_subsection("Physical data"); 
      { 
        initial_time = prm.get_double("initial_time"); 
        final_time   = prm.get_double("final_time"); 
        Reynolds     = prm.get_double("Reynolds"); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Time step data"); 
      { 
        dt = prm.get_double("dt"); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Space discretization"); 
      { 
        n_global_refines = prm.get_integer("n_of_refines"); 
        pressure_degree  = prm.get_integer("pressure_fe_degree"); 
      } 
      prm.leave_subsection(); 

      prm.enter_subsection("Data solve velocity"); 
      { 
        vel_max_iterations = prm.get_integer("max_iterations"); 
        vel_eps            = prm.get_double("eps"); 
        vel_Krylov_size    = prm.get_integer("Krylov_size"); 
        vel_off_diagonals  = prm.get_integer("off_diagonals"); 
        vel_diag_strength  = prm.get_double("diag_strength"); 
        vel_update_prec    = prm.get_integer("update_prec"); 
      } 
      prm.leave_subsection(); 

      verbose = prm.get_bool("verbose"); 

      output_interval = prm.get_integer("output_interval"); 
    } 
  } // namespace RunTimeParameters 

//  @sect3{Equation data}  

// 在下一个命名空间中，我们声明初始和边界条件。

  namespace EquationData 
  { 

// 由于我们选择了一个完全解耦的公式，我们将不利用deal.II处理矢量值问题的能力。然而，我们确实希望为方程数据使用一个独立于维度的接口。为了做到这一点，我们的函数应该能够知道我们目前在哪个空间分量上工作，而且我们应该能够有一个通用的接口来做到这一点。下面的类是在这个方向上的一个尝试。

    template <int dim> 
    class MultiComponentFunction : public Function<dim> 
    { 
    public: 
      MultiComponentFunction(const double initial_time = 0.); 
      void set_component(const unsigned int d); 

    protected: 
      unsigned int comp; 
    }; 

    template <int dim> 
    MultiComponentFunction<dim>::MultiComponentFunction( 
      const double initial_time) 
      : Function<dim>(1, initial_time) 
      , comp(0) 
    {} 

    template <int dim> 
    void MultiComponentFunction<dim>::set_component(const unsigned int d) 
    { 
      Assert(d < dim, ExcIndexRange(d, 0, dim)); 
      comp = d; 
    } 

// 有了这个类的定义，我们声明描述速度和压力的边界条件的类。

    template <int dim> 
    class Velocity : public MultiComponentFunction<dim> 
    { 
    public: 
      Velocity(const double initial_time = 0.0); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void value_list(const std::vector<Point<dim>> &points, 
                              std::vector<double> &          values, 
                              const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    Velocity<dim>::Velocity(const double initial_time) 
      : MultiComponentFunction<dim>(initial_time) 
    {} 

    template <int dim> 
    void Velocity<dim>::value_list(const std::vector<Point<dim>> &points, 
                                   std::vector<double> &          values, 
                                   const unsigned int) const 
    { 
      const unsigned int n_points = points.size(); 
      Assert(values.size() == n_points, 
             ExcDimensionMismatch(values.size(), n_points)); 
      for (unsigned int i = 0; i < n_points; ++i) 
        values[i] = Velocity<dim>::value(points[i]); 
    } 

    template <int dim> 
    double Velocity<dim>::value(const Point<dim> &p, const unsigned int) const 
    { 
      if (this->comp == 0) 
        { 
          const double Um = 1.5; 
          const double H  = 4.1; 
          return 4. * Um * p(1) * (H - p(1)) / (H * H); 
        } 
      else 
        return 0.; 
    } 

    template <int dim> 
    class Pressure : public Function<dim> 
    { 
    public: 
      Pressure(const double initial_time = 0.0); 

      virtual double value(const Point<dim> & p, 
                           const unsigned int component = 0) const override; 

      virtual void value_list(const std::vector<Point<dim>> &points, 
                              std::vector<double> &          values, 
                              const unsigned int component = 0) const override; 
    }; 

    template <int dim> 
    Pressure<dim>::Pressure(const double initial_time) 
      : Function<dim>(1, initial_time) 
    {} 

    template <int dim> 
    double Pressure<dim>::value(const Point<dim> & p, 
                                const unsigned int component) const 
    { 
      (void)component; 
      AssertIndexRange(component, 1); 
      return 25. - p(0); 
    } 

    template <int dim> 
    void Pressure<dim>::value_list(const std::vector<Point<dim>> &points, 
                                   std::vector<double> &          values, 
                                   const unsigned int component) const 
    { 
      (void)component; 
      AssertIndexRange(component, 1); 
      const unsigned int n_points = points.size(); 
      Assert(values.size() == n_points, 
             ExcDimensionMismatch(values.size(), n_points)); 
      for (unsigned int i = 0; i < n_points; ++i) 
        values[i] = Pressure<dim>::value(points[i]); 
    } 
  } // namespace EquationData 

//  @sect3{The <code>NavierStokesProjection</code> class}  

// 现在是程序的主类。它实现了纳维-斯托克斯方程的各种版本的投影方法。考虑到介绍中给出的实现细节，所有方法和成员变量的名称应该是不言自明的。

  template <int dim> 
  class NavierStokesProjection 
  { 
  public: 
    NavierStokesProjection(const RunTimeParameters::Data_Storage &data); 

    void run(const bool verbose = false, const unsigned int n_plots = 10); 

  protected: 
    RunTimeParameters::Method type; 

    const unsigned int deg; 
    const double       dt; 
    const double       t_0; 
    const double       T; 
    const double       Re; 

    EquationData::Velocity<dim>               vel_exact; 
    std::map<types::global_dof_index, double> boundary_values; 
    std::vector<types::boundary_id>           boundary_ids; 

    Triangulation<dim> triangulation; 

    FE_Q<dim> fe_velocity; 
    FE_Q<dim> fe_pressure; 

    DoFHandler<dim> dof_handler_velocity; 
    DoFHandler<dim> dof_handler_pressure; 

    QGauss<dim> quadrature_pressure; 
    QGauss<dim> quadrature_velocity; 

    SparsityPattern sparsity_pattern_velocity; 
    SparsityPattern sparsity_pattern_pressure; 
    SparsityPattern sparsity_pattern_pres_vel; 

    SparseMatrix<double> vel_Laplace_plus_Mass; 
    SparseMatrix<double> vel_it_matrix[dim]; 
    SparseMatrix<double> vel_Mass; 
    SparseMatrix<double> vel_Laplace; 
    SparseMatrix<double> vel_Advection; 
    SparseMatrix<double> pres_Laplace; 
    SparseMatrix<double> pres_Mass; 
    SparseMatrix<double> pres_Diff[dim]; 
    SparseMatrix<double> pres_iterative; 

    Vector<double> pres_n; 
    Vector<double> pres_n_minus_1; 
    Vector<double> phi_n; 
    Vector<double> phi_n_minus_1; 
    Vector<double> u_n[dim]; 
    Vector<double> u_n_minus_1[dim]; 
    Vector<double> u_star[dim]; 
    Vector<double> force[dim]; 
    Vector<double> v_tmp; 
    Vector<double> pres_tmp; 
    Vector<double> rot_u; 

    SparseILU<double>   prec_velocity[dim]; 
    SparseILU<double>   prec_pres_Laplace; 
    SparseDirectUMFPACK prec_mass; 
    SparseDirectUMFPACK prec_vel_mass; 

    DeclException2(ExcInvalidTimeStep, 
                   double, 
                   double, 
                   << " The time step " << arg1 << " is out of range." 
                   << std::endl 
                   << " The permitted range is (0," << arg2 << "]"); 

    void create_triangulation_and_dofs(const unsigned int n_refines); 

    void initialize(); 

    void interpolate_velocity(); 

    void diffusion_step(const bool reinit_prec); 

    void projection_step(const bool reinit_prec); 

    void update_pressure(const bool reinit_prec); 

  private: 
    unsigned int vel_max_its; 
    unsigned int vel_Krylov_size; 
    unsigned int vel_off_diagonals; 
    unsigned int vel_update_prec; 
    double       vel_eps; 
    double       vel_diag_strength; 

    void initialize_velocity_matrices(); 

    void initialize_pressure_matrices(); 

// 接下来的几个结构和函数是用来做各种并行的事情。它们遵循 @ref threads 中规定的方案，使用WorkStream类。正如那里所解释的，这需要我们为每个装配器声明两个结构，一个是每个任务的数据，一个是scratch数据结构。然后，这些结构被移交给组装本地贡献的函数，并将这些本地贡献复制到全局对象上。

// 这个程序的一个特殊之处在于，我们并不是只有一个DoFHandler对象来代表速度和压力，而是为这两种变量使用单独的DoFHandler对象。当我们想把涉及这两个变量的条款，如速度的发散和压力的梯度，乘以各自的测试函数时，我们要为这种优化付费。在这样做的时候，我们不能只使用一个FEValues对象，而是需要两个，而且需要用单元格迭代器来初始化它们，这些单元格迭代器指向三角形中的同一个单元格，但不同的DoFHandlers。

// 为了在实践中做到这一点，我们声明一个 "同步 "迭代器--一个内部由几个（在我们的例子中是两个）迭代器组成的对象，每次同步迭代器向前移动一步，内部存储的每个迭代器也向前移动一步，从而始终保持同步。碰巧的是，有一个deal.II类可以促进这种事情。这里重要的是要知道，建立在同一个三角形上的两个DoFHandler对象将以相同的顺序走过三角形的单元。

    using IteratorTuple = 
      std::tuple<typename DoFHandler<dim>::active_cell_iterator, 
                 typename DoFHandler<dim>::active_cell_iterator>; 

    using IteratorPair = SynchronousIterators<IteratorTuple>; 

    void initialize_gradient_operator(); 

    struct InitGradPerTaskData 
    { 
      unsigned int                         d; 
      unsigned int                         vel_dpc; 
      unsigned int                         pres_dpc; 
      FullMatrix<double>                   local_grad; 
      std::vector<types::global_dof_index> vel_local_dof_indices; 
      std::vector<types::global_dof_index> pres_local_dof_indices; 

      InitGradPerTaskData(const unsigned int dd, 
                          const unsigned int vdpc, 
                          const unsigned int pdpc) 
        : d(dd) 
        , vel_dpc(vdpc) 
        , pres_dpc(pdpc) 
        , local_grad(vdpc, pdpc) 
        , vel_local_dof_indices(vdpc) 
        , pres_local_dof_indices(pdpc) 
      {} 
    }; 

    struct InitGradScratchData 
    { 
      unsigned int  nqp; 
      FEValues<dim> fe_val_vel; 
      FEValues<dim> fe_val_pres; 
      InitGradScratchData(const FE_Q<dim> &  fe_v, 
                          const FE_Q<dim> &  fe_p, 
                          const QGauss<dim> &quad, 
                          const UpdateFlags  flags_v, 
                          const UpdateFlags  flags_p) 
        : nqp(quad.size()) 
        , fe_val_vel(fe_v, quad, flags_v) 
        , fe_val_pres(fe_p, quad, flags_p) 
      {} 
      InitGradScratchData(const InitGradScratchData &data) 
        : nqp(data.nqp) 
        , fe_val_vel(data.fe_val_vel.get_fe(), 
                     data.fe_val_vel.get_quadrature(), 
                     data.fe_val_vel.get_update_flags()) 
        , fe_val_pres(data.fe_val_pres.get_fe(), 
                      data.fe_val_pres.get_quadrature(), 
                      data.fe_val_pres.get_update_flags()) 
      {} 
    }; 

    void assemble_one_cell_of_gradient(const IteratorPair & SI, 
                                       InitGradScratchData &scratch, 
                                       InitGradPerTaskData &data); 

    void copy_gradient_local_to_global(const InitGradPerTaskData &data); 

// 同样的一般布局也适用于以下实现平流项组装的类和函数。

    void assemble_advection_term(); 

    struct AdvectionPerTaskData 
    { 
      FullMatrix<double>                   local_advection; 
      std::vector<types::global_dof_index> local_dof_indices; 
      AdvectionPerTaskData(const unsigned int dpc) 
        : local_advection(dpc, dpc) 
        , local_dof_indices(dpc) 
      {} 
    }; 

    struct AdvectionScratchData 
    { 
      unsigned int                nqp; 
      unsigned int                dpc; 
      std::vector<Point<dim>>     u_star_local; 
      std::vector<Tensor<1, dim>> grad_u_star; 
      std::vector<double>         u_star_tmp; 
      FEValues<dim>               fe_val; 
      AdvectionScratchData(const FE_Q<dim> &  fe, 
                           const QGauss<dim> &quad, 
                           const UpdateFlags  flags) 
        : nqp(quad.size()) 
        , dpc(fe.n_dofs_per_cell()) 
        , u_star_local(nqp) 
        , grad_u_star(nqp) 
        , u_star_tmp(nqp) 
        , fe_val(fe, quad, flags) 
      {} 

      AdvectionScratchData(const AdvectionScratchData &data) 
        : nqp(data.nqp) 
        , dpc(data.dpc) 
        , u_star_local(nqp) 
        , grad_u_star(nqp) 
        , u_star_tmp(nqp) 
        , fe_val(data.fe_val.get_fe(), 
                 data.fe_val.get_quadrature(), 
                 data.fe_val.get_update_flags()) 
      {} 
    }; 

    void assemble_one_cell_of_advection( 
      const typename DoFHandler<dim>::active_cell_iterator &cell, 
      AdvectionScratchData &                                scratch, 
      AdvectionPerTaskData &                                data); 

    void copy_advection_local_to_global(const AdvectionPerTaskData &data); 

// 最后几个函数实现了扩散解以及输出的后处理，包括计算速度的曲线。

    void diffusion_component_solve(const unsigned int d); 

    void output_results(const unsigned int step); 

    void assemble_vorticity(const bool reinit_prec); 
  }; 

//  @sect4{ <code>NavierStokesProjection::NavierStokesProjection</code> }  

// 在构造函数中，我们只是从作为参数传递的 <code>Data_Storage</code> 对象中读取所有数据，验证我们读取的数据是否合理，最后，创建三角形并加载初始数据。

  template <int dim> 
  NavierStokesProjection<dim>::NavierStokesProjection( 
    const RunTimeParameters::Data_Storage &data) 
    : type(data.form) 
    , deg(data.pressure_degree) 
    , dt(data.dt) 
    , t_0(data.initial_time) 
    , T(data.final_time) 
    , Re(data.Reynolds) 
    , vel_exact(data.initial_time) 
    , fe_velocity(deg + 1) 
    , fe_pressure(deg) 
    , dof_handler_velocity(triangulation) 
    , dof_handler_pressure(triangulation) 
    , quadrature_pressure(deg + 1) 
    , quadrature_velocity(deg + 2) 
    , vel_max_its(data.vel_max_iterations) 
    , vel_Krylov_size(data.vel_Krylov_size) 
    , vel_off_diagonals(data.vel_off_diagonals) 
    , vel_update_prec(data.vel_update_prec) 
    , vel_eps(data.vel_eps) 
    , vel_diag_strength(data.vel_diag_strength) 
  { 
    if (deg < 1) 
      std::cout 
        << " WARNING: The chosen pair of finite element spaces is not stable." 
        << std::endl 
        << " The obtained results will be nonsense" << std::endl; 

    AssertThrow(!((dt <= 0.) || (dt > .5 * T)), ExcInvalidTimeStep(dt, .5 * T)); 

    create_triangulation_and_dofs(data.n_global_refines); 
    initialize(); 
  } 
// @sect4{<code>NavierStokesProjection::create_triangulation_and_dofs</code>}  

// 创建三角形的方法，并将其细化到所需的次数。在创建三角形之后，它创建了与网格相关的数据，即分配自由度和重新编号，并初始化我们将使用的矩阵和向量。

  template <int dim> 
  void NavierStokesProjection<dim>::create_triangulation_and_dofs( 
    const unsigned int n_refines) 
  { 
    GridIn<dim> grid_in; 
    grid_in.attach_triangulation(triangulation); 

    { 
      std::string   filename = "nsbench2.inp"; 
      std::ifstream file(filename); 
      Assert(file, ExcFileNotOpen(filename.c_str())); 
      grid_in.read_ucd(file); 
    } 

    std::cout << "Number of refines = " << n_refines << std::endl; 
    triangulation.refine_global(n_refines); 
    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl; 

    boundary_ids = triangulation.get_boundary_ids(); 

 
    DoFRenumbering::boost::Cuthill_McKee(dof_handler_velocity); 
    dof_handler_pressure.distribute_dofs(fe_pressure); 
    DoFRenumbering::boost::Cuthill_McKee(dof_handler_pressure); 

    initialize_velocity_matrices(); 
    initialize_pressure_matrices(); 
    initialize_gradient_operator(); 

    pres_n.reinit(dof_handler_pressure.n_dofs()); 
    pres_n_minus_1.reinit(dof_handler_pressure.n_dofs()); 
    phi_n.reinit(dof_handler_pressure.n_dofs()); 
    phi_n_minus_1.reinit(dof_handler_pressure.n_dofs()); 
    pres_tmp.reinit(dof_handler_pressure.n_dofs()); 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        u_n[d].reinit(dof_handler_velocity.n_dofs()); 
        u_n_minus_1[d].reinit(dof_handler_velocity.n_dofs()); 
        u_star[d].reinit(dof_handler_velocity.n_dofs()); 
        force[d].reinit(dof_handler_velocity.n_dofs()); 
      } 
    v_tmp.reinit(dof_handler_velocity.n_dofs()); 
    rot_u.reinit(dof_handler_velocity.n_dofs()); 

    std::cout << "dim (X_h) = " << (dof_handler_velocity.n_dofs() * dim) // 
              << std::endl                                               // 
              << "dim (M_h) = " << dof_handler_pressure.n_dofs()         // 
              << std::endl                                               // 
              << "Re        = " << Re << std::endl                       // 
              << std::endl; 
  } 
// @sect4{ <code>NavierStokesProjection::initialize</code> }  

// 该方法创建常数矩阵并加载初始数据。

  template <int dim> 
  void NavierStokesProjection<dim>::initialize() 
  { 
    vel_Laplace_plus_Mass = 0.; 
    vel_Laplace_plus_Mass.add(1. / Re, vel_Laplace); 
    vel_Laplace_plus_Mass.add(1.5 / dt, vel_Mass); 

    EquationData::Pressure<dim> pres(t_0); 
    VectorTools::interpolate(dof_handler_pressure, pres, pres_n_minus_1); 
    pres.advance_time(dt); 
    VectorTools::interpolate(dof_handler_pressure, pres, pres_n); 
    phi_n         = 0.; 
    phi_n_minus_1 = 0.; 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        vel_exact.set_time(t_0); 
        vel_exact.set_component(d); 
        VectorTools::interpolate(dof_handler_velocity, 
                                 vel_exact, 
                                 u_n_minus_1[d]); 
        vel_exact.advance_time(dt); 
        VectorTools::interpolate(dof_handler_velocity, vel_exact, u_n[d]); 
      } 
  } 
// @sect4{ <code>NavierStokesProjection::initialize_*_matrices</code> }  

// 在这组方法中，我们初始化了稀疏模式、约束条件（如果有的话）并组装了不依赖于时间步长的矩阵  <code>dt</code>  。请注意，对于拉普拉斯矩阵和质量矩阵，我们可以使用库中的函数来做这件事。因为这个函数的昂贵操作--创建两个矩阵--是完全独立的，我们原则上可以把它们标记为可以使用 Threads::new_task 函数进行%并行工作的任务。我们不会在这里这样做，因为这些函数在内部已经被并行化了，特别是由于当前的函数在每个程序运行中只被调用一次，所以在每个时间步长中不会产生费用。然而，必要的修改将是非常直接的。

  template <int dim> 
  void NavierStokesProjection<dim>::initialize_velocity_matrices() 
  { 
    { 
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(), 
                                 dof_handler_velocity.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler_velocity, dsp); 
      sparsity_pattern_velocity.copy_from(dsp); 
    } 
    vel_Laplace_plus_Mass.reinit(sparsity_pattern_velocity); 
    for (unsigned int d = 0; d < dim; ++d) 
      vel_it_matrix[d].reinit(sparsity_pattern_velocity); 
    vel_Mass.reinit(sparsity_pattern_velocity); 
    vel_Laplace.reinit(sparsity_pattern_velocity); 
    vel_Advection.reinit(sparsity_pattern_velocity); 

    MatrixCreator::create_mass_matrix(dof_handler_velocity, 
                                      quadrature_velocity, 
                                      vel_Mass); 
    MatrixCreator::create_laplace_matrix(dof_handler_velocity, 
                                         quadrature_velocity, 
                                         vel_Laplace); 
  } 

//作用于压力空间的矩阵的初始化与作用于速度空间的矩阵相似。

  template <int dim> 
  void NavierStokesProjection<dim>::initialize_pressure_matrices() 
  { 
    { 
      DynamicSparsityPattern dsp(dof_handler_pressure.n_dofs(), 
                                 dof_handler_pressure.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler_pressure, dsp); 
      sparsity_pattern_pressure.copy_from(dsp); 
    } 

    pres_Laplace.reinit(sparsity_pattern_pressure); 
    pres_iterative.reinit(sparsity_pattern_pressure); 
    pres_Mass.reinit(sparsity_pattern_pressure); 

    MatrixCreator::create_laplace_matrix(dof_handler_pressure, 
                                         quadrature_pressure, 
                                         pres_Laplace); 
    MatrixCreator::create_mass_matrix(dof_handler_pressure, 
                                      quadrature_pressure, 
                                      pres_Mass); 
  } 

// 对于梯度算子，我们从初始化稀疏模式和压缩它开始。这里需要注意的是，梯度算子从压力空间作用到速度空间，所以我们必须处理两个不同的有限元空间。为了保持循环的同步，我们使用之前定义的别名，即 <code>PairedIterators</code> and <code>IteratorPair</code>  。

  template <int dim> 
  void NavierStokesProjection<dim>::initialize_gradient_operator() 
  { 
    { 
      DynamicSparsityPattern dsp(dof_handler_velocity.n_dofs(), 
                                 dof_handler_pressure.n_dofs()); 
      DoFTools::make_sparsity_pattern(dof_handler_velocity, 
                                      dof_handler_pressure, 
                                      dsp); 
      sparsity_pattern_pres_vel.copy_from(dsp); 
    } 

    InitGradPerTaskData per_task_data(0, 
                                      fe_velocity.n_dofs_per_cell(), 
                                      fe_pressure.n_dofs_per_cell()); 
    InitGradScratchData scratch_data(fe_velocity, 
                                     fe_pressure, 
                                     quadrature_velocity, 
                                     update_gradients | update_JxW_values, 
                                     update_values); 

    for (unsigned int d = 0; d < dim; ++d) 
      { 
        pres_Diff[d].reinit(sparsity_pattern_pres_vel); 
        per_task_data.d = d; 
        WorkStream::run( 
          IteratorPair(IteratorTuple(dof_handler_velocity.begin_active(), 
                                     dof_handler_pressure.begin_active())), 
          IteratorPair(IteratorTuple(dof_handler_velocity.end(), 
                                     dof_handler_pressure.end())), 
          *this, 
          &NavierStokesProjection<dim>::assemble_one_cell_of_gradient, 
          &NavierStokesProjection<dim>::copy_gradient_local_to_global, 
          scratch_data, 
          per_task_data); 
      } 
  } 

  template <int dim> 
  void NavierStokesProjection<dim>::assemble_one_cell_of_gradient( 
    const IteratorPair & SI, 
    InitGradScratchData &scratch, 
    InitGradPerTaskData &data) 
  { 
    scratch.fe_val_vel.reinit(std::get<0>(*SI)); 
    scratch.fe_val_pres.reinit(std::get<1>(*SI)); 

    std::get<0>(*SI)->get_dof_indices(data.vel_local_dof_indices); 
    std::get<1>(*SI)->get_dof_indices(data.pres_local_dof_indices); 

    data.local_grad = 0.; 
    for (unsigned int q = 0; q < scratch.nqp; ++q) 
      { 
        for (unsigned int i = 0; i < data.vel_dpc; ++i) 
          for (unsigned int j = 0; j < data.pres_dpc; ++j) 
            data.local_grad(i, j) += 
              -scratch.fe_val_vel.JxW(q) * 
              scratch.fe_val_vel.shape_grad(i, q)[data.d] * 
              scratch.fe_val_pres.shape_value(j, q); 
      } 
  } 

  template <int dim> 
  void NavierStokesProjection<dim>::copy_gradient_local_to_global( 
    const InitGradPerTaskData &data) 
  { 
    for (unsigned int i = 0; i < data.vel_dpc; ++i) 
      for (unsigned int j = 0; j < data.pres_dpc; ++j) 
        pres_Diff[data.d].add(data.vel_local_dof_indices[i], 
                              data.pres_local_dof_indices[j], 
                              data.local_grad(i, j)); 
  } 
// @sect4{ <code>NavierStokesProjection::run</code> }  

// 这是时间行进函数，从 <code>t_0</code> 开始，使用时间步长 <code>dt</code> 的投影法在时间上前进，直到 <code>T</code>  。

// 它的第二个参数 <code>verbose</code> 表示该函数是否应该输出它在任何特定时刻正在做什么的信息：例如，它将说明我们是否正在进行扩散、投影子步骤；更新前置条件器等等。我们没有使用像
// @code
//    if (verbose) std::cout << "something";
//  @endcode
//  那样的代码来实现这种输出，而是使用ConditionalOStream类来为我们做这个。该类接受一个输出流和一个条件，该条件表明你传递给它的东西是否应该被传递到给定的输出流，或者应该被忽略。这样，上面的代码就变成了
//  @code
//    verbose_cout << "something";
//  @endcode，并且在任何情况下都会做正确的事情。

  template <int dim> 
  void NavierStokesProjection<dim>::run(const bool         verbose, 
                                        const unsigned int output_interval) 
  { 
    ConditionalOStream verbose_cout(std::cout, verbose); 

    const auto n_steps = static_cast<unsigned int>((T - t_0) / dt); 
    vel_exact.set_time(2. * dt); 
    output_results(1); 
    for (unsigned int n = 2; n <= n_steps; ++n) 
      { 
        if (n % output_interval == 0) 
          { 
            verbose_cout << "Plotting Solution" << std::endl; 
            output_results(n); 
          } 
        std::cout << "Step = " << n << " Time = " << (n * dt) << std::endl; 
        verbose_cout << "  Interpolating the velocity " << std::endl; 

        interpolate_velocity(); 
        verbose_cout << "  Diffusion Step" << std::endl; 
        if (n % vel_update_prec == 0) 
          verbose_cout << "    With reinitialization of the preconditioner" 
                       << std::endl; 
        diffusion_step((n % vel_update_prec == 0) || (n == 2)); 
        verbose_cout << "  Projection Step" << std::endl; 
        projection_step((n == 2)); 
        verbose_cout << "  Updating the Pressure" << std::endl; 
        update_pressure((n == 2)); 
        vel_exact.advance_time(dt); 
      } 
    output_results(n_steps); 
  } 

  template <int dim> 
  void NavierStokesProjection<dim>::interpolate_velocity() 
  { 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        u_star[d].equ(2., u_n[d]); 
        u_star[d] -= u_n_minus_1[d]; 
      } 
  } 
// @sect4{<code>NavierStokesProjection::diffusion_step</code>}  

// 扩散步骤的实现。请注意，昂贵的操作是函数末尾的扩散解，我们必须为每个速度分量做一次。为了加快进度，我们允许以%并行方式进行，使用 Threads::new_task 函数，确保 <code>dim</code> 的求解都得到处理，并被安排到可用的处理器上：如果你的机器有一个以上的处理器核心，并且这个程序的其他部分目前没有使用资源，那么扩散求解将以%并行方式运行。另一方面，如果你的系统只有一个处理器核心，那么以%并行方式运行将是低效的（因为它导致了，例如，缓存拥堵），事情将被顺序地执行。

  template <int dim> 
  void NavierStokesProjection<dim>::diffusion_step(const bool reinit_prec) 
  { 
    pres_tmp.equ(-1., pres_n); 
    pres_tmp.add(-4. / 3., phi_n, 1. / 3., phi_n_minus_1); 

    assemble_advection_term(); 

    for (unsigned int d = 0; d < dim; ++d) 
      { 
        force[d] = 0.; 
        v_tmp.equ(2. / dt, u_n[d]); 
        v_tmp.add(-.5 / dt, u_n_minus_1[d]); 
        vel_Mass.vmult_add(force[d], v_tmp); 

        pres_Diff[d].vmult_add(force[d], pres_tmp); 
        u_n_minus_1[d] = u_n[d]; 

        vel_it_matrix[d].copy_from(vel_Laplace_plus_Mass); 
        vel_it_matrix[d].add(1., vel_Advection); 

        vel_exact.set_component(d); 
        boundary_values.clear(); 
        for (const auto &boundary_id : boundary_ids) 
          { 
            switch (boundary_id) 
              { 
                case 1: 
                  VectorTools::interpolate_boundary_values( 
                    dof_handler_velocity, 
                    boundary_id, 
                    Functions::ZeroFunction<dim>(), 
                    boundary_values); 
                  break; 
                case 2: 
                  VectorTools::interpolate_boundary_values(dof_handler_velocity, 
                                                           boundary_id, 
                                                           vel_exact, 
                                                           boundary_values); 
                  break; 
                case 3: 
                  if (d != 0) 
                    VectorTools::interpolate_boundary_values( 
                      dof_handler_velocity, 
                      boundary_id, 
                      Functions::ZeroFunction<dim>(), 
                      boundary_values); 
                  break; 
                case 4: 
                  VectorTools::interpolate_boundary_values( 
                    dof_handler_velocity, 
                    boundary_id, 
                    Functions::ZeroFunction<dim>(), 
                    boundary_values); 
                  break; 
                default: 
                  Assert(false, ExcNotImplemented()); 
              } 
          } 
        MatrixTools::apply_boundary_values(boundary_values, 
                                           vel_it_matrix[d], 
                                           u_n[d], 
                                           force[d]); 
      } 

    Threads::TaskGroup<void> tasks; 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        if (reinit_prec) 
          prec_velocity[d].initialize(vel_it_matrix[d], 
                                      SparseILU<double>::AdditionalData( 
                                        vel_diag_strength, vel_off_diagonals)); 
        tasks += Threads::new_task( 
          &NavierStokesProjection<dim>::diffusion_component_solve, *this, d); 
      } 
    tasks.join_all(); 
  } 

  template <int dim> 
  void 
  NavierStokesProjection<dim>::diffusion_component_solve(const unsigned int d) 
  { 
    SolverControl solver_control(vel_max_its, vel_eps * force[d].l2_norm()); 
    SolverGMRES<Vector<double>> gmres( 
      solver_control, 
      SolverGMRES<Vector<double>>::AdditionalData(vel_Krylov_size)); 
    gmres.solve(vel_it_matrix[d], u_n[d], force[d], prec_velocity[d]); 
  } 
// @sect4{ <code>NavierStokesProjection::assemble_advection_term</code> }  

// 下面的几个函数是关于集合平流项的，平流项是扩散步骤的系统矩阵的一部分，在每个时间步骤中都会发生变化。如上所述，我们将使用WorkStream类和文件模块 @ref threads 中描述的其他设施，在所有单元上平行运行装配循环。

  template <int dim> 
  void NavierStokesProjection<dim>::assemble_advection_term() 
  { 
    vel_Advection = 0.; 
    AdvectionPerTaskData data(fe_velocity.n_dofs_per_cell()); 
    AdvectionScratchData scratch(fe_velocity, 
                                 quadrature_velocity, 
                                 update_values | update_JxW_values | 
                                   update_gradients); 
    WorkStream::run( 
      dof_handler_velocity.begin_active(), 
      dof_handler_velocity.end(), 
      *this, 
      &NavierStokesProjection<dim>::assemble_one_cell_of_advection, 
      &NavierStokesProjection<dim>::copy_advection_local_to_global, 
      scratch, 
      data); 
  } 

  template <int dim> 
  void NavierStokesProjection<dim>::assemble_one_cell_of_advection( 
    const typename DoFHandler<dim>::active_cell_iterator &cell, 
    AdvectionScratchData &                                scratch, 
    AdvectionPerTaskData &                                data) 
  { 
    scratch.fe_val.reinit(cell); 
    cell->get_dof_indices(data.local_dof_indices); 
    for (unsigned int d = 0; d < dim; ++d) 
      { 
        scratch.fe_val.get_function_values(u_star[d], scratch.u_star_tmp); 
        for (unsigned int q = 0; q < scratch.nqp; ++q) 
          scratch.u_star_local[q](d) = scratch.u_star_tmp[q]; 
      } 

    for (unsigned int d = 0; d < dim; ++d) 
      { 
        scratch.fe_val.get_function_gradients(u_star[d], scratch.grad_u_star); 
        for (unsigned int q = 0; q < scratch.nqp; ++q) 
          { 
            if (d == 0) 
              scratch.u_star_tmp[q] = 0.; 
            scratch.u_star_tmp[q] += scratch.grad_u_star[q][d]; 
          } 
      } 

    data.local_advection = 0.; 
    for (unsigned int q = 0; q < scratch.nqp; ++q) 
      for (unsigned int i = 0; i < scratch.dpc; ++i) 
        for (unsigned int j = 0; j < scratch.dpc; ++j) 
          data.local_advection(i, j) += (scratch.u_star_local[q] *            // 
                                           scratch.fe_val.shape_grad(j, q) *  // 
                                           scratch.fe_val.shape_value(i, q)   // 
                                         +                                    // 
                                         0.5 *                                // 
                                           scratch.u_star_tmp[q] *            // 
                                           scratch.fe_val.shape_value(i, q) * // 
                                           scratch.fe_val.shape_value(j, q))  // 
                                        * scratch.fe_val.JxW(q); 
  } 

  template <int dim> 
  void NavierStokesProjection<dim>::copy_advection_local_to_global( 
    const AdvectionPerTaskData &data) 
  { 
    for (unsigned int i = 0; i < fe_velocity.n_dofs_per_cell(); ++i) 
      for (unsigned int j = 0; j < fe_velocity.n_dofs_per_cell(); ++j) 
        vel_Advection.add(data.local_dof_indices[i], 
                          data.local_dof_indices[j], 
                          data.local_advection(i, j)); 
  } 

//  @sect4{<code>NavierStokesProjection::projection_step</code>}  

// 这实现了投影的步骤。

  template <int dim> 
  void NavierStokesProjection<dim>::projection_step(const bool reinit_prec) 
  { 
    pres_iterative.copy_from(pres_Laplace); 

    pres_tmp = 0.; 
    for (unsigned d = 0; d < dim; ++d) 
      pres_Diff[d].Tvmult_add(pres_tmp, u_n[d]); 

    phi_n_minus_1 = phi_n; 

    static std::map<types::global_dof_index, double> bval; 
    if (reinit_prec) 
      VectorTools::interpolate_boundary_values(dof_handler_pressure, 
                                               3, 
                                               Functions::ZeroFunction<dim>(), 
                                               bval); 

    MatrixTools::apply_boundary_values(bval, pres_iterative, phi_n, pres_tmp); 

    if (reinit_prec) 
      prec_pres_Laplace.initialize(pres_iterative, 
                                   SparseILU<double>::AdditionalData( 
                                     vel_diag_strength, vel_off_diagonals)); 

    SolverControl solvercontrol(vel_max_its, vel_eps * pres_tmp.l2_norm()); 
    SolverCG<Vector<double>> cg(solvercontrol); 
    cg.solve(pres_iterative, phi_n, pres_tmp, prec_pres_Laplace); 

    phi_n *= 1.5 / dt; 
  } 
// @sect4{ <code>NavierStokesProjection::update_pressure</code> }  

// 这是投影法的压力更新步骤。它实现了该方法的标准表述，即
//  @f[ p^{n+1} = p^n +
//  \phi^{n+1}, 
//  @f]
//  或旋转形式，即
//  @f[ p^{n+1} = p^n +
//  \phi^{n+1} - \frac{1}{Re} \nabla\cdot u^{n+1}. 
//  @f] 

  template <int dim> 
  void NavierStokesProjection<dim>::update_pressure(const bool reinit_prec) 
  { 
    pres_n_minus_1 = pres_n; 
    switch (type) 
      { 
        case RunTimeParameters::Method::standard: 
          pres_n += phi_n; 
          break; 
        case RunTimeParameters::Method::rotational: 
          if (reinit_prec) 
            prec_mass.initialize(pres_Mass); 
          pres_n = pres_tmp; 
          prec_mass.solve(pres_n); 
          pres_n.sadd(1. / Re, 1., pres_n_minus_1); 
          pres_n += phi_n; 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      }; 
  } 
// @sect4{ <code>NavierStokesProjection::output_results</code> }  

// 该方法绘制了当前的解决方案。主要的困难是，我们想创建一个单一的输出文件，其中包含所有的速度分量、压力以及流动的涡度的数据。另一方面，速度和压力存在于不同的DoFHandler对象中，因此不能用一个DataOut对象写入同一个文件。因此，我们必须更努力地把各种数据放到一个DoFHandler对象中，然后用它来驱动图形输出。

// 我们不会在这里详细说明这个过程，而是参考  step-32  ，那里使用了一个类似的程序（并有记录），为所有变量创建一个联合的 DoFHandler 对象。

// 我们还注意到，我们在这里将涡度作为一个单独的函数中的标量来计算，使用 $L^2$ 量的投影 $\text{curl} u$ 到用于速度成分的有限元空间。但原则上，我们也可以从速度中计算出一个点状量，并通过 step-29 和 step-33 中讨论的DataPostprocessor机制实现。

  template <int dim> 
  void NavierStokesProjection<dim>::output_results(const unsigned int step) 
  { 
    assemble_vorticity((step == 1)); 
    const FESystem<dim> joint_fe( 
      fe_velocity, dim, fe_pressure, 1, fe_velocity, 1); 
    DoFHandler<dim> joint_dof_handler(triangulation); 
    joint_dof_handler.distribute_dofs(joint_fe); 
    Assert(joint_dof_handler.n_dofs() == 
             ((dim + 1) * dof_handler_velocity.n_dofs() + 
              dof_handler_pressure.n_dofs()), 
           ExcInternalError()); 
    Vector<double> joint_solution(joint_dof_handler.n_dofs()); 
    std::vector<types::global_dof_index> loc_joint_dof_indices( 
      joint_fe.n_dofs_per_cell()), 
      loc_vel_dof_indices(fe_velocity.n_dofs_per_cell()), 
      loc_pres_dof_indices(fe_pressure.n_dofs_per_cell()); 
    typename DoFHandler<dim>::active_cell_iterator 
      joint_cell = joint_dof_handler.begin_active(), 
      joint_endc = joint_dof_handler.end(), 
      vel_cell   = dof_handler_velocity.begin_active(), 
      pres_cell  = dof_handler_pressure.begin_active(); 
    for (; joint_cell != joint_endc; ++joint_cell, ++vel_cell, ++pres_cell) 
      { 
        joint_cell->get_dof_indices(loc_joint_dof_indices); 
        vel_cell->get_dof_indices(loc_vel_dof_indices); 
        pres_cell->get_dof_indices(loc_pres_dof_indices); 
        for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i) 
          switch (joint_fe.system_to_base_index(i).first.first) 
            { 
              case 0: 
                Assert(joint_fe.system_to_base_index(i).first.second < dim, 
                       ExcInternalError()); 
                joint_solution(loc_joint_dof_indices[i]) = 
                  u_n[joint_fe.system_to_base_index(i).first.second]( 
                    loc_vel_dof_indices[joint_fe.system_to_base_index(i) 
                                          .second]); 
                break; 
              case 1: 
                Assert(joint_fe.system_to_base_index(i).first.second == 0, 
                       ExcInternalError()); 
                joint_solution(loc_joint_dof_indices[i]) = 
                  pres_n(loc_pres_dof_indices[joint_fe.system_to_base_index(i) 
                                                .second]); 
                break; 
              case 2: 
                Assert(joint_fe.system_to_base_index(i).first.second == 0, 
                       ExcInternalError()); 
                joint_solution(loc_joint_dof_indices[i]) = rot_u( 
                  loc_vel_dof_indices[joint_fe.system_to_base_index(i).second]); 
                break; 
              default: 
                Assert(false, ExcInternalError()); 
            } 
      } 
    std::vector<std::string> joint_solution_names(dim, "v"); 
    joint_solution_names.emplace_back("p"); 
    joint_solution_names.emplace_back("rot_u"); 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(joint_dof_handler); 
    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      component_interpretation( 
        dim + 2, DataComponentInterpretation::component_is_part_of_vector); 
    component_interpretation[dim] = 
      DataComponentInterpretation::component_is_scalar; 
    component_interpretation[dim + 1] = 
      DataComponentInterpretation::component_is_scalar; 
    data_out.add_data_vector(joint_solution, 
                             joint_solution_names, 
                             DataOut<dim>::type_dof_data, 
                             component_interpretation); 
    data_out.build_patches(deg + 1); 
    std::ofstream output("solution-" + Utilities::int_to_string(step, 5) + 
                         ".vtk"); 
    data_out.write_vtk(output); 
  } 

// 下面是一个辅助函数，通过将 $\text{curl} u$ 项投影到用于速度分量的有限元空间来计算涡度。这个函数只有在我们生成图形输出时才会被调用，所以不是很频繁，因此我们没有像对待其他装配函数那样，麻烦地使用WorkStream概念来并行化它。不过，如果需要的话，这应该不会太复杂。此外，我们在这里的实现只适用于2D，所以如果不是这种情况，我们就放弃了。

  template <int dim> 
  void NavierStokesProjection<dim>::assemble_vorticity(const bool reinit_prec) 
  { 
    Assert(dim == 2, ExcNotImplemented()); 
    if (reinit_prec) 
      prec_vel_mass.initialize(vel_Mass); 

    FEValues<dim>      fe_val_vel(fe_velocity, 
                             quadrature_velocity, 
                             update_gradients | update_JxW_values | 
                               update_values); 
    const unsigned int dpc = fe_velocity.n_dofs_per_cell(), 
                       nqp = quadrature_velocity.size(); 
    std::vector<types::global_dof_index> ldi(dpc); 
    Vector<double>                       loc_rot(dpc); 

    std::vector<Tensor<1, dim>> grad_u1(nqp), grad_u2(nqp); 
    rot_u = 0.; 

    for (const auto &cell : dof_handler_velocity.active_cell_iterators()) 
      { 
        fe_val_vel.reinit(cell); 
        cell->get_dof_indices(ldi); 
        fe_val_vel.get_function_gradients(u_n[0], grad_u1); 
        fe_val_vel.get_function_gradients(u_n[1], grad_u2); 
        loc_rot = 0.; 
        for (unsigned int q = 0; q < nqp; ++q) 
          for (unsigned int i = 0; i < dpc; ++i) 
            loc_rot(i) += (grad_u2[q][0] - grad_u1[q][1]) * // 
                          fe_val_vel.shape_value(i, q) *    // 
                          fe_val_vel.JxW(q); 

        for (unsigned int i = 0; i < dpc; ++i) 
          rot_u(ldi[i]) += loc_rot(i); 
      } 

    prec_vel_mass.solve(rot_u); 
  } 
} // namespace Step35 
// @sect3{ The main function }  

// 主函数看起来和其他所有的教程程序非常相似，所以这里没有什么可评论的。

int main() 
{ 
  try 
    { 
      using namespace Step35; 

      RunTimeParameters::Data_Storage data; 
      data.read_data("parameter-file.prm"); 

      deallog.depth_console(data.verbose ? 2 : 0); 

      NavierStokesProjection<2> test(data); 
      test.run(data.verbose, data.output_interval); 
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
  std::cout << "----------------------------------------------------" 
            << std::endl 
            << "Apparently everything went fine!" << std::endl 
            << "Don't forget to brush your teeth :-)" << std::endl 
            << std::endl; 
  return 0; 
} 

