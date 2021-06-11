

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 - 2021 by the deal.II authors
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
 * Authors: Luca Heltai, Bruno Blais, Rene Gassmoeller, 2020
 */


// @sect3{Include files}
// 其中大部分已经在其他地方介绍过了，我们只对新的部分进行评论。靠近顶部的开关允许在
// PETSc 和 Trilinos 线性代数功能之间进行选择，这与  step-40  和  step-50
// 中的开关类似。

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/generic_linear_algebra.h>
#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/linear_operator_tools.h>

#define FORCE_USE_OF_TRILINOS

namespace LA
{
#if defined(DEAL_II_WITH_PETSC) && !defined(DEAL_II_PETSC_WITH_COMPLEX) && \ 
  !(defined(DEAL_II_WITH_TRILINOS) && defined(FORCE_USE_OF_TRILINOS))
  using namespace dealii::LinearAlgebraPETSc;
#  define USE_PETSC_LA
#elif defined(DEAL_II_WITH_TRILINOS)
  using namespace dealii::LinearAlgebraTrilinos;
#else
#  error DEAL_II_WITH_PETSC or DEAL_II_WITH_TRILINOS required
#endif
} // namespace LA

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe_field.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_minres.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/vector_tools.h>

// 这些是关于  step-60
// 的唯一新的包含文件。在本教程中，实体和流体之间的非匹配耦合是通过一个中间数据结构来计算的，该结构记录了实体的正交点在流体网格中的位置如何演变。这个数据结构需要跟踪描述实体域的每个单元上的正交点的位置、正交权重，如果实体域是同维度的，还需要跟踪每个点的法向量。

// Deal.II通过ParticleHandler类在Particles命名空间中提供这些设施。ParticleHandler是一个允许你管理粒子集合的类（类型为
// Particles::Particle),
// 的对象，代表具有一些附加属性（如id）的点的集合，漂浮在一个
// parallel::distributed::Triangulation.
// 命名空间中的方法和类允许人们轻松实现Particle-In-Cell方法和在分布式三角形上的粒子追踪。

// 我们 "滥用
// "这个数据结构来存储嵌入周围流体网格中的实体正交点的位置信息，包括积分权重，以及可能的表面法线。我们之所以使用这个额外的数据结构，是因为实体网格和流体网格可能是不重叠的，如果我们使用两个独立的三角计算对象，那么它们将独立地分布在并行进程中。

// 为了耦合这两个问题，我们依靠ParticleHandler类，在每个粒子中存储一个实体正交点的位置（一般来说，它不与任何流体正交点对齐），它的权重，以及耦合这两个问题可能需要的任何其他信息。这些位置然后与固体叶轮的（规定）速度一起传播。

// 固体正交点的所有权最初是从固体网格本身的MPI分区中继承的。这样产生的粒子后来通过ParticleHandler类的方法分配到流体网格中。这允许MPI进程之间透明地交换关于流体单元和实体正交点之间的重叠模式的信息。

#include <deal.II/particles/data_out.h>
#include <deal.II/particles/generators.h>
#include <deal.II/particles/particle_handler.h>
#include <deal.II/particles/utilities.h>

// 在生成网格时，我们允许从文件中读取它，如果deal.II已经建立了OpenCASCADE支持，我们也允许读取CAD文件，并将它们作为网格的流形描述符（参见
// step-54 对OpenCASCADE命名空间中的各种流形描述符的详细描述）。

#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>
#ifdef DEAL_II_WITH_OPENCASCADE
#  include <TopoDS.hxx>
#endif

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>

namespace Step70
{
  using namespace dealii;
  // @sect3{Run-time parameter handling}

  // 与我们在 step-60
  // 中所做的类似，我们建立了一个持有我们问题的所有参数的类，并从ParameterAcceptor类中派生出来以简化参数文件的管理和创建。

  // ParameterAcceptor范式要求所有参数都可以被ParameterAcceptor方法写入。为了避免出现很难追踪的错误（比如写成`time
  // = 0`而不是`time ==
  // 0`），我们在一个外部类中声明所有的参数，该类在实际的`StokesImmersedProblem`类之前被初始化，并将其作为`const`引用传递给主类。

  // 该类的构造函数负责该类的成员与ParameterHandler中的相应条目之间的连接。由于使用了
  // ParameterHandler::add_parameter()
  // 方法，这种连接是微不足道的，但要求这个类的所有成员都是可写的。

  template <int dim, int spacedim = dim>
  class StokesImmersedProblemParameters : public ParameterAcceptor
  {
  public:
    StokesImmersedProblemParameters();

    // 然而，由于这个类将作为一个`const`引用传递给StokesImmersedProblem类，我们必须确保我们仍然可以在这里定义的Function类派生的对象中正确设置时间。为了做到这一点，我们声明
    // `StokesImmersedProblemParameters::rhs` 和
    // `StokesImmersedProblemParameters::angular_velocity` 成员都是 "可变
    // "的，并定义以下的小辅助方法，将它们的时间设置为正确的值。

    void
    set_time(const double &time) const
    {
      rhs.set_time(time);
      angular_velocity.set_time(time);
    }

    // 该类的其余部分主要由描述模拟及其离散化细节的成员变量组成。下面的参数是关于输出的位置、空间和时间离散化（默认是
    // $Q_2\times Q_1$
    // Taylor-Hood离散化，它使用2度的多项式来计算速度），以及在我们再次生成图形输出之前应该经过多少时间步长。

    std::string output_directory = ".";

    unsigned int velocity_degree = 2;

    unsigned int number_of_time_steps = 501;
    double       final_time           = 1.0;

    unsigned int output_frequency = 1;

    // 我们允许每个网格独立地被细化。在本教程中，固体网格上没有解决物理问题，其速度被作为基准点给出。然而，在本教程中加入一些弹性模型，并将其转化为一个完全成熟的FSI求解器是相对简单的。

    unsigned int initial_fluid_refinement      = 5;
    unsigned int initial_solid_refinement      = 5;
    unsigned int particle_insertion_refinement = 3;

    // 为了提供对流体领域的粗略描述，我们使用extract_rtree_level()方法，该方法适用于流体三角结构中每个局部拥有的单元的边界盒树。树的级别越高，提取的边界盒数量就越多，对流体领域的描述也就越准确。然而，大量的边界盒也意味着巨大的通信成本，因为边界盒的收集是由所有进程收集的。

    unsigned int fluid_rtree_extraction_level = 1;

    // 方程中使用的唯一两个数值参数是流体的粘度，以及Nitsche公式中使用的惩罚项
    // $\beta$ 。

    double viscosity    = 1.0;
    double penalty_term = 100;

    // 默认情况下，我们创建一个没有着色的hyper_cube，并且我们使用同质的Dirichlet边界条件。在这个集合中，我们存储了设置边界条件时要使用的边界ID。

    std::list<types::boundary_id> homogeneous_dirichlet_ids{0};

    // 我们在此说明另一种从参数文件中创建三角形的方法，使用
    // GridGenerator::generate_from_name_and_arguments(),
    // ，该方法接收GridGenerator命名空间中的函数名称，其参数为一个字符串，代表参数的元组。

    // 在 Patterns::Tools::Convert
    // 类中详细解释了将参数从字符串解析成字符串的机制，该类用于将字符串翻译成大多数基本STL类型（向量、映射、图元）和基本deal.II类型（点、张量、BoundingBox等）。

    // 一般来说，可以用等级1的统一元素表示的对象（即 std::vector<double>,
    // Point<dim>,  std::set<int>,
    // 等）是用逗号分开的。额外的等级采取分号，允许你将字符串解析为
    // `std::vector<std::vector<double>>`, 或例如 `std::vector<Point<dim>>`,
    // 类型的对象，如`0.0, 0.1; 0.1,
    // 0.2`。这个字符串可以被解释为两个Point对象的向量，或者一个双数向量的向量。

    // 当条目不统一时，比如在元组的情况下，我们用冒号来分隔各个条目。例如，像`5:
    // 0.1, 0.2`这样的字符串可以用来解析一个类型为 `std::pair<int,
    // Point<2>>的对象或者一个 `std::tuple<int,  的对象。 std::vector<double>>`.

    // 在我们的例子中，大多数参数是点对象（代表中心、角、细分元素等）、整数值（细分数量）、双倍值（半径、长度等）或布尔选项（如许多GridGenerator函数采取的`colorize`选项）。

    // 在下面的例子中，我们设置了合理的默认值，但这些值可以在运行时通过选择GridGenerator命名空间的任何其他支持的函数来改变。如果GridGenerator函数失败，本程序将把网格的名称解释为vtk网格文件名，把参数解释为从manifold_id到描述域的几何形状的CAD文件的映射。每个CAD文件都将被分析，并根据CAD文件本身的内容生成OpenCASCADE命名空间的Manifold。

    // 为了尽可能的通用，我们对每个生成的网格都这样做：流体网格、固体网格，但也包括使用三角法生成的示踪粒子。

    std::string name_of_fluid_grid       = "hyper_cube";
    std::string arguments_for_fluid_grid = "-1: 1: false";
    std::string name_of_solid_grid       = "hyper_rectangle";
    std::string arguments_for_solid_grid = spacedim == 2 ?
                                             "-.5, -.1: .5, .1: false" :
                                             "-.5, -.1, -.1: .5, .1, .1: false";
    std::string name_of_particle_grid    = "hyper_ball";
    std::string arguments_for_particle_grid =
      spacedim == 2 ? "0.3, 0.3: 0.1: false" : "0.3, 0.3, 0.3 : 0.1: false";

    // 同样地，我们允许不同的局部细化策略。特别是，我们限制了细化水平的最大数量，以控制流体网格的最小尺寸，并保证它与实体网格兼容。细化级数的最小值也得到了控制，以确保在流动的大部分地区有足够的精度。此外，我们根据流体速度场的标准误差估计器进行局部细化。

    // 我们允许用户选择两种最常见的细化策略，即 "fixed_number "或
    // "fixed_fraction"，这两种策略参考了
    // GridRefinement::refine_and_coarsen_fixed_fraction() 和
    // GridRefinement::refine_and_coarsen_fixed_number(). 方法。

    // 细化可以每隔几个时间步骤进行一次，而不是连续进行，我们通过`细化_频率`参数来控制这个值。

    int          max_level_refinement = 8;
    int          min_level_refinement = 5;
    std::string  refinement_strategy  = "fixed_fraction";
    double       coarsening_fraction  = 0.3;
    double       refinement_fraction  = 0.3;
    unsigned int max_cells            = 20000;
    int          refinement_frequency = 5;

    // 最后，以下两个函数对象被用来控制斯托克斯流的源项和我们移动固体体的角速度。在一个更现实的模拟中，实体速度或其变形将来自于实体域上的辅助问题的解决。在这个例子中，我们把这部分放在一边，只是在浸没的固体上沿Z轴施加一个固定的旋转速度场，由一个可以在参数文件中指定的函数来控制。

    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> rhs;
    mutable ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>>
      angular_velocity;
  };

  // 还有一个任务就是声明我们在输入文件中可以接受哪些运行时参数。我们将这些参数分成不同的类别，把它们放在ParameterHandler类的不同部分。我们首先在全局范围内声明StokesImmersedProblem使用的所有全局参数。

  template <int dim, int spacedim>
  StokesImmersedProblemParameters<dim,
                                  spacedim>::StokesImmersedProblemParameters()
    : ParameterAcceptor("Stokes Immersed Problem/")
    , rhs("Right hand side", spacedim + 1)
    , angular_velocity("Angular velocity")
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", this->prm, Patterns::Integer(1));

    add_parameter("Number of time steps", number_of_time_steps);
    add_parameter("Output frequency", output_frequency);

    add_parameter("Output directory", output_directory);

    add_parameter("Final time", final_time);

    add_parameter("Viscosity", viscosity);

    add_parameter("Nitsche penalty term", penalty_term);

    add_parameter("Initial fluid refinement",
                  initial_fluid_refinement,
                  "Initial mesh refinement used for the fluid domain Omega");

    add_parameter("Initial solid refinement",
                  initial_solid_refinement,
                  "Initial mesh refinement used for the solid domain Gamma");

    add_parameter("Fluid bounding boxes extraction level",
                  fluid_rtree_extraction_level,
                  "Extraction level of the rtree used to construct global "
                  "bounding boxes");

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles");

    add_parameter(
      "Homogeneous Dirichlet boundary ids",
      homogeneous_dirichlet_ids,
      "Boundary Ids over which homogeneous Dirichlet boundary conditions are applied");

    // 下一节专门介绍用于创建各种网格的参数。我们将需要三种不同的三角形。流体网格
    // "用于定义流体领域，"固体网格 "用于定义固体领域，"粒子网格
    // "用于分布一些示踪粒子，这些粒子随速度漂移，只作为被动示踪物使用。

    enter_subsection("Grid generation");
    {
      add_parameter("Fluid grid generator", name_of_fluid_grid);
      add_parameter("Fluid grid generator arguments", arguments_for_fluid_grid);

      add_parameter("Solid grid generator", name_of_solid_grid);
      add_parameter("Solid grid generator arguments", arguments_for_solid_grid);


      add_parameter("Particle grid generator arguments",
                    arguments_for_particle_grid);
    }
    leave_subsection();

    enter_subsection("Refinement and remeshing");
    {
      add_parameter("Refinement step frequency", refinement_frequency);
      add_parameter("Refinement maximal level", max_level_refinement);
      add_parameter("Refinement minimal level", min_level_refinement);
      add_parameter("Refinement strategy",
                    refinement_strategy,
                    "",
                    this->prm,
                    Patterns::Selection("fixed_fraction|fixed_number"));
      add_parameter("Refinement coarsening fraction", coarsening_fraction);
      add_parameter("Refinement fraction", refinement_fraction);
      add_parameter("Maximum number of cells", max_cells);
    }
    leave_subsection();

    // 最后的任务是修正右侧函数的默认尺寸，并定义一个有意义的默认角速度而不是零。

    rhs.declare_parameters_call_back.connect([&]() {
      Functions::ParsedFunction<spacedim>::declare_parameters(this->prm,
                                                              spacedim + 1);
    });
    angular_velocity.declare_parameters_call_back.connect([&]() {
      this->prm.set("Function expression",
                    "t < .500001 ? 6.283185 : -6.283185");
    });
  }

  // 一旦角速度被提供为一个函数对象，我们就通过下面这个派生自函数类的类来重建点状实体速度。它通过假设实体以给定的角速度绕原点（或3D中的
  // $z$ 轴）旋转，提供实体在给定位置的速度值。

  template <int spacedim>
  class SolidVelocity : public Function<spacedim>
  {
  public:
    static_assert(spacedim > 1,
                  "Cannot instantiate SolidVelocity for spacedim == 1");

    SolidVelocity(const Functions::ParsedFunction<spacedim> &angular_velocity)
      : angular_velocity(angular_velocity)
    {}

    virtual double
    value(const Point<spacedim> &p, unsigned int component = 0) const override
    {
      Tensor<1, spacedim> velocity;

      // 我们假设角速度是沿Z轴方向的，也就是说，我们把实际的角速度模拟成二维旋转，而不考虑`spacedim`的实际值。

      const double omega = angular_velocity.value(p);
      velocity[0]        = -omega * p[1];
      velocity[1]        = omega * p[0];

      return velocity[component];
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
  };

  // 同样地，我们假设固体的位置可以在每个时间步长明确地计算出来，利用角速度的知识。我们计算固体粒子的确切位置，假定固体的旋转量等于时间步长乘以在`p`点计算的角速度。

  template <int spacedim>
  class SolidPosition : public Function<spacedim>
  {
  public:
    static_assert(spacedim > 1,
                  "Cannot instantiate SolidPosition for spacedim == 1");

    SolidPosition(const Functions::ParsedFunction<spacedim> &angular_velocity,
                  const double                               time_step)
      : Function<spacedim>(spacedim)
      , angular_velocity(angular_velocity)
      , time_step(time_step)
    {}

    virtual double
    value(const Point<spacedim> &p, unsigned int component = 0) const override
    {
      Point<spacedim> new_position = p;

      double dtheta = angular_velocity.value(p) * time_step;

      new_position[0] = std::cos(dtheta) * p[0] - std::sin(dtheta) * p[1];
      new_position[1] = std::sin(dtheta) * p[0] + std::cos(dtheta) * p[1];

      return new_position[component];
    }

    void
    set_time_step(const double new_time_step)
    {
      time_step = new_time_step;
    }

  private:
    const Functions::ParsedFunction<spacedim> &angular_velocity;
    double                                     time_step;
  };
  // @sect3{The StokesImmersedProblem class declaration}

  // 我们现在准备介绍我们的教程程序的主类。像往常一样，除了构造函数外，我们只留下一个公共入口：`run()`方法。其他的都是
  // "私有 "的，并通过run方法本身进行访问。

  template <int dim, int spacedim = dim>
  class StokesImmersedProblem
  {
  public:
    StokesImmersedProblem(
      const StokesImmersedProblemParameters<dim, spacedim> &par);

    void
    run();

    // 接下来的部分包含了该类的`private`成员。第一个方法类似于前一个例子中的方法。然而，它不仅负责生成流体的网格，而且还负责生成固体的网格。第二个方法是计算最大的时间步长，保证每个粒子最多移动一个单元。这对于确保
    // Particles::ParticleHandler
    // 能够找到粒子最终所在的单元是非常重要的，因为它只能从一个单元看向它的近邻（因为在并行设置中，每个MPI进程只知道它拥有的单元以及它们的近邻）。

  private:
    void
    make_grid();

    double
    compute_time_step() const;

    // 接下来的两个函数将初始化这个类中使用的 Particles::ParticleHandler
    // 对象。我们有两个这样的对象。一个代表被动追踪器，用于绘制流体粒子的轨迹，而另一个代表固体的材料粒子，它们被放置在固体网格的正交点上。

    void
    setup_tracer_particles();
    void
    setup_solid_particles();

    // 剩下的设置分为两部分。以下两个函数中的第一个创建了每次模拟需要的所有对象，而另一个则设置了所有需要在每个细化步骤中重新初始化的对象。

    void
    initial_setup();
    void
    setup_dofs();

    // 装配例程与其他斯托克斯装配例程非常相似，但Nitsche限制部分除外，它利用其中一个粒子处理程序在流体域的非匹配部分进行积分，对应于固体的位置。我们将这两部分分成两个独立的函数。

    void
    assemble_stokes_system();
    void
    assemble_nitsche_restriction();

    // 其余的函数求解线性系统（看起来与 step-60
    // 中的线性系统几乎相同），然后对解进行后处理。refine_and_transfer()方法仅在每一个`refinement_frequency`步骤中被调用，以适应网格，并确保所有在细化前的时间步骤中计算的场都正确地转移到新的网格中。这包括矢量场，以及粒子信息。同样地，我们每隔`output_frequency`步就会调用两个输出方法。

    void
    solve();

    void
    refine_and_transfer();

    void
    output_results(const unsigned int cycle, const double time) const;
    void
    output_particles(const Particles::ParticleHandler<spacedim> &particles,
                     std::string                                 fprefix,
                     const unsigned int                          iter,
                     const double                                time) const;

    //接下来让我们来看看这个类的成员函数。第一个是处理从参数文件中读取的运行时参数。如前所述，我们通过使其成为一个`const`引用来确保我们不能从这个类中修改这个对象。

    const StokesImmersedProblemParameters<dim, spacedim> &par;

    // 然后还有MPI通信器对象，如果程序是并行运行的，我们将用它来让进程在网络上发送信息，还有`pcout`对象和定时器信息，也被
    // step-40 采用，例如。

    MPI_Comm mpi_communicator;

    ConditionalOStream pcout;

    mutable TimerOutput computing_timer;

    // 接下来是关于  step-60
    // 的主要创新点之一。这里我们假设固体和流体都是完全分布的三角形。这使得问题可以扩展到非常大的自由度，代价是要沟通所有非匹配三角形之间的重叠区域。这一点特别棘手，因为我们没有对两个三角形的各个子域的相对位置或分布做出假设。特别是，我们假设每个进程只拥有
    // "solid_tria "的一部分，以及 "fluid_tria
    // "的一部分，不一定在同一个物理区域，也不一定重叠。

    // 我们原则上可以尝试创建初始分区，使每个过程的子域在固体和流体区域之间重叠。然而，这种重叠在模拟过程中会被破坏，我们将不得不一次又一次地重新分配DoF。我们在本教程中采用的方法更加灵活，而且成本也不高。我们在模拟开始时进行两次全对全的通信，以交换每个处理器的几何占用信息（近似）（通过包围盒的集合完成）。

    // 这个信息被 Particles::ParticleHandler
    // 类用来交换（使用某对某的通信模式）所有的粒子，因此每个进程都知道生活在它所拥有的流体子域所占区域的粒子。

    // 为了把重叠的区域连接起来，我们利用了ParticleHandler类中实现的设施。

    parallel::distributed::Triangulation<spacedim>      fluid_tria;
    parallel::distributed::Triangulation<dim, spacedim> solid_tria;

    // 接下来是对所使用的有限元的描述，以及适当的正交公式和相应的DoFHandler对象。在目前的实现中，只有`fluid_fe`是真正必要的。为了完整起见，并便于扩展，我们还保留了`solid_fe`，但它被初始化为一个FE_Nothing有限元空间，即没有自由度的空间。

    // 我们将这两个有限元空间声明为 `std::unique_ptr`
    // 对象，而不是普通的成员变量，以便在`StokesImmersedProblemParameters'被初始化后生成它们。特别是，它们将在`initial_setup()`方法中被初始化。

    std::unique_ptr<FiniteElement<spacedim>>      fluid_fe;
    std::unique_ptr<FiniteElement<dim, spacedim>> solid_fe;

    std::unique_ptr<Quadrature<spacedim>> fluid_quadrature_formula;
    std::unique_ptr<Quadrature<dim>>      solid_quadrature_formula;

    DoFHandler<spacedim>      fluid_dh;
    DoFHandler<dim, spacedim> solid_dh;

    std::unique_ptr<MappingFEField<dim, spacedim>> solid_mapping;

    // 与 step-22
    // 中的做法类似，我们使用一个块状系统来处理问题的斯托克斯部分，并非常密切地遵循那里的做法。

    std::vector<IndexSet> fluid_owned_dofs;
    std::vector<IndexSet> solid_owned_dofs;

    std::vector<IndexSet> fluid_relevant_dofs;
    std::vector<IndexSet> solid_relevant_dofs;

    // 利用这种自由度的划分，我们就可以定义所有必要的对象来描述有关的线性系统。

    AffineConstraints<double> constraints;

    LA::MPI::BlockSparseMatrix system_matrix;
    LA::MPI::BlockSparseMatrix preconditioner_matrix;

    LA::MPI::BlockVector solution;
    LA::MPI::BlockVector locally_relevant_solution;
    LA::MPI::BlockVector system_rhs;

    // 让我们转到这个程序的粒子方面。有两个 Particles::ParticleHandler
    // 对象用于耦合固体和流体，以及描述被动追踪器。在许多方面，这些对象的作用类似于离散化中使用的DoFHandler类，也就是说，它们提供了粒子的枚举，并允许查询每个粒子的信息。

    Particles::ParticleHandler<spacedim> tracer_particle_handler;
    Particles::ParticleHandler<spacedim> solid_particle_handler;

    // 对于每个追踪器粒子，我们需要计算其当前位置的速度场，并使用离散时间步进方案更新其位置。我们使用分布式线性代数对象来做这件事，这些对象存储了每个粒子的位置或速度的坐标。也就是说，这些向量有`tracer_particle_handler.n_global_particles()
    // * spacedim`项，我们将以一种方式来存储这些向量的一部分，以便在所有进程中进行划分。(隐含地，我们在此假设每个粒子的`spacedim'坐标被存储在向量的连续条目中)。因此，我们需要确定每个向量条目的所有者是谁。我们将这个所有者设定为等于在时间 $t=0$ 产生该粒子的进程。这个信息对每一个进程都存储在`locally_owned_tracer_particle_coordinates`索引集里。

    // 一旦粒子被分配到与拥有粒子所在区域的进程相匹配，我们将需要从该进程读取相应的速度场。我们通过填充一个只读的速度矢量场来实现这一目标，该矢量场包含了Ghost条目中的相关信息。这是通过`locally_relevant_tracer_particle_coordinates`索引集实现的，该索引集记录了模拟过程中的变化情况，也就是说，它记录了当前进程拥有的粒子最终出现在哪里，以及谁拥有最终出现在我的子域的粒子。

    // 虽然这不是最有效的策略，但我们保持这种方式是为了说明事情在真实的流固耦合（FSI）问题中是如何运作的。如果一个粒子与一个特定的固体自由度相联系，我们就不能自由选择谁拥有它，我们必须把这个信息传达给周围的人。我们在这里说明了这一点，并表明通信模式是点对点的，就算法的总成本而言可以忽略不计。

    // 然后，基于这些细分定义的向量被用来存储粒子的速度（只读，有幽灵条目）和它们的位移（读/写，没有幽灵条目）。

    IndexSet locally_owned_tracer_particle_coordinates;
    IndexSet locally_relevant_tracer_particle_coordinates;

    LA::MPI::Vector tracer_particle_velocities;
    LA::MPI::Vector relevant_tracer_particle_displacements;

    // 本教程程序的关键点之一是两个独立的 parallel::distributed::Triangulation
    // 对象之间的耦合，其中一个对象可能相对于另一个对象移动和变形（可能有较大的变形）。当流体和实体的三角形都是
    // parallel::distributed::Triangulation,
    // 类型时，每个进程只能访问这两个三角形中每个单元的局部拥有的部分。如上所述，一般情况下，本地拥有的域是不重叠的。

    // 为了允许在不重叠的 parallel::distributed::Triangulation
    // 对象之间有效地交换信息，该库的一些算法要求用户提供三角形的本地拥有部分所占区域的粗略描述，其形式是每个进程的轴对齐的边界盒集合，这些边界盒提供了域的本地拥有部分的完整覆盖。这种信息就可以用于这样的情况：人们需要向已知位置周围的单元格的所有者发送信息，而不知道这个所有者实际上是谁。但是，如果我们知道每个进程拥有的几何区域或体积的边界盒集合，那么我们就可以确定可能拥有该位置所在单元的所有进程的一个子集：即其边界盒包含该点的所有进程。与其向所有进程发送与该位置相关的信息，不如只向具有点对点通信基元的一小部分进程发送信息。你会注意到，这也允许典型的时间与内存的权衡：我们愿意存储的关于每个进程拥有的区域的数据越多--以更精细的边界框信息的形式--我们必须执行的通信就越少）。

    // 我们通过收集一个向量（长度为 Utilities::MPI::n_mpi_processes())
    // 的BoundingBox对象的向量）来构建这些信息。我们用extract_rtree_level()函数填充这个向量，并允许用户选择要提取的树的哪一级。这个
    // "级别 "对应的是与边界框重叠的区域应该有多粗/多细。

    // 作为一个例子，这是由extract_rtree_level()函数对一个分布在三个过程中的二维超球所提取的结果。每张图片中，绿色显示的是与每个进程上的三角形的本地所有单元相关的边界框，紫色显示的是从rtree中提取的边界框。

    //  @image html rtree-process-0.png
    // @image html rtree-process-1.png
    // @image html rtree-process-2.png

    // 我们将这些盒子存储在一个全局成员变量中，在每个细化步骤中都会更新。

    std::vector<std::vector<BoundingBox<spacedim>>> global_fluid_bounding_boxes;
  };

  //  @sect3{The StokesImmersedProblem class implementation}
  // @sect4{Object construction and mesh initialization functions}

  // 在构造函数中，我们创建了mpi_communicator，以及流体和实体的三角计算和dof_handler。通过使用mpi_communicator，我们构建了ConditionalOStream和TimerOutput对象。

  template <int dim, int spacedim>
  StokesImmersedProblem<dim, spacedim>::StokesImmersedProblem(
    const StokesImmersedProblemParameters<dim, spacedim> &par)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , pcout(std::cout,
            (Utilities::MPI::this_mpi_process(mpi_communicator) == 0))
    , computing_timer(mpi_communicator,
                      pcout,
                      TimerOutput::summary,
                      TimerOutput::wall_times)
    , fluid_tria(mpi_communicator,
                 typename Triangulation<spacedim>::MeshSmoothing(
                   Triangulation<spacedim>::smoothing_on_refinement |
                   Triangulation<spacedim>::smoothing_on_coarsening))
    , solid_tria(mpi_communicator,
                 typename Triangulation<dim, spacedim>::MeshSmoothing(
                   Triangulation<dim, spacedim>::smoothing_on_refinement |
                   Triangulation<dim, spacedim>::smoothing_on_coarsening))
    , fluid_dh(fluid_tria)
    , solid_dh(solid_tria)
  {}

  // 为了生成网格，我们首先尝试使用deal.II
  // GridGenerator命名空间中的函数，通过利用
  // GridGenerator::generate_from_name_and_argument().
  // 如果这个函数失败，那么我们使用以下方法，名称被解释为文件名，参数被解释为从流形ID到CAD文件的映射，并使用OpenCASCADE命名空间设施转换为流形描述符。在顶部，我们把文件读成一个三角图。

  template <int dim, int spacedim>
  void
  read_grid_and_cad_files(const std::string &           grid_file_name,
                          const std::string &           ids_and_cad_file_names,
                          Triangulation<dim, spacedim> &tria)
  {
    GridIn<dim, spacedim> grid_in;
    grid_in.attach_triangulation(tria);
    grid_in.read(grid_file_name);

    // 如果我们走到这一步，那么三角图已经被读取，我们已经准备好将正确的流形描述附加到它上面。只有在deal.II支持OpenCASCADE的情况下，我们才会执行接下来的几行代码。对于地图中的每个条目，我们尝试打开相应的CAD文件，分析它，并根据其内容，选择一个
    // OpenCASCADE::ArcLengthProjectionLineManifold
    // （如果CAD文件包含一个`TopoDS_Edge'或一个`TopoDS_Wire'）或一个
    // OpenCASCADE::NURBSPatchManifold,
    // ，如果文件包含一个面。请注意，如果CAD文件不包含单一的线、边或面，在生成Manifold时将会抛出一个断言。

    // 我们使用 Patterns::Tools::Convert
    // 类来完成从字符串到歧管ID和文件名之间的映射的转换。

#ifdef DEAL_II_WITH_OPENCASCADE
    using map_type  = std::map<types::manifold_id, std::string>;
    using Converter = Patterns::Tools::Convert<map_type>;

    for (const auto &pair : Converter::to_value(ids_and_cad_file_names))
      {
        const auto &manifold_id   = pair.first;
        const auto &cad_file_name = pair.second;

        const auto extension = boost::algorithm::to_lower_copy(
          cad_file_name.substr(cad_file_name.find_last_of('.') + 1));

        TopoDS_Shape shape;
        if (extension == "iges" || extension == "igs")
          shape = OpenCASCADE::read_IGES(cad_file_name);
        else if (extension == "step" || extension == "stp")
          shape = OpenCASCADE::read_STEP(cad_file_name);
        else
          AssertThrow(false,
                      ExcNotImplemented("We found an extension that we "
                                        "do not recognize as a CAD file "
                                        "extension. Bailing out."));

        // 现在我们检查一下这个 "形状
        // "中包含了多少个面。OpenCASCADE本质上是三维的，所以如果这个数字是零，我们就把它解释为线状流形，否则就解释为`spacedim`=3中的
        // OpenCASCADE::NormalToMeshProjectionManifold ，或者`spacedim`=2中的
        // OpenCASCADE::NURBSPatchManifold 。

        const auto n_elements = OpenCASCADE::count_elements(shape);
        if ((std::get<0>(n_elements) == 0))
          tria.set_manifold(
            manifold_id,
            OpenCASCADE::ArclengthProjectionLineManifold<dim, spacedim>(shape));
        else if (spacedim == 3)
          {
            // 我们使用这个技巧，因为
            // OpenCASCADE::NormalToMeshProjectionManifold 只在spacedim =
            // 3的情况下实现。上面的检查保证了事情的实际运作是正确的。

            const auto t = reinterpret_cast<Triangulation<dim, 3> *>(&tria);
            t->set_manifold(manifold_id,
                            OpenCASCADE::NormalToMeshProjectionManifold<dim, 3>(
                              shape));
          }
        else

          // 我们也允许基于单个NURBS补丁的二维空间的曲面描述。要做到这一点，CAD文件必须包含一个单一的`TopoDS_Face`。

          tria.set_manifold(manifold_id,
                            OpenCASCADE::NURBSPatchManifold<dim, spacedim>(
                              TopoDS::Face(shape)));
      }
#else
    (void)ids_and_cad_file_names;
    AssertThrow(false, ExcNotImplemented("Generation of the grid failed."));
#endif
  }

  // 现在让我们把东西放在一起，并制作所有必要的网格。如上所述，我们首先尝试在内部生成网格，如果我们失败了（即如果我们最终进入了`catch'子句），那么我们就继续执行上述函数。

  // 我们对流体和固体网格都重复这个模式。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::make_grid()
  {
    try
      {
        GridGenerator::generate_from_name_and_arguments(
          fluid_tria, par.name_of_fluid_grid, par.arguments_for_fluid_grid);
      }
    catch (...)
      {
        pcout << "Generating from name and argument failed." << std::endl
              << "Trying to read from file name." << std::endl;
        read_grid_and_cad_files(par.name_of_fluid_grid,
                                par.arguments_for_fluid_grid,
                                fluid_tria);
      }
    fluid_tria.refine_global(par.initial_fluid_refinement);

    try
      {
        GridGenerator::generate_from_name_and_arguments(
          solid_tria, par.name_of_solid_grid, par.arguments_for_solid_grid);
      }
    catch (...)
      {
        read_grid_and_cad_files(par.name_of_solid_grid,
                                par.arguments_for_solid_grid,
                                solid_tria);
      }

    solid_tria.refine_global(par.initial_solid_refinement);
  }
  // @sect4{Particle initialization functions}

  // 一旦固体和流体网格被创建，我们就开始填充 Particles::ParticleHandler
  // 对象。我们要处理的第一个对象是用来跟踪流体中的被动追踪器的对象。这些东西只是沿途传送，从某种意义上说，它们的位置并不重要：我们只是想用它们来观察流体被传送的位置。我们可以使用任何我们选择的方式来确定它们的初始位置。一个方便的方法是将初始位置创建为我们所选择的形状的网格顶点，这个选择由参数文件中的一个运行时参数决定。

  // 在这个实现中，我们使用FE_Q有限元空间的支持点来创建追踪器，这些支持点定义在一个临时网格上，然后被丢弃。在这个网格中，我们只保留与支撑点相关的
  // Particles::Particle 对象（存储在 Particles::ParticleHandler 类中）。

  //  Particles::ParticleHandler
  //  类提供了插入一组粒子的可能性，这些粒子实际生活在活动过程所拥有的域的一部分。然而，在这种情况下，这个功能是不够的。作为任意网格（与流体网格不匹配）上的FE_Q对象的本地拥有的支持点所产生的粒子没有理由位于流体网格的本地拥有的子域的同一物理区域内。事实上，这种情况几乎不会发生，尤其是我们想要跟踪粒子本身发生了什么。

  // 在粒子入室方法（PIC）中，人们通常习惯于将粒子的所有权分配给粒子所在的过程。在本教程中，我们说明了一种不同的方法，如果想跟踪与粒子有关的信息，这种方法是很有用的（例如，如果一个粒子与一个特定的自由度有关，而这个自由度是由一个特定的过程所拥有的，不一定是在任何特定时间拥有该粒子所在的流体单元的同一个过程）。在这里使用的方法中，粒子的所有权在开始时被分配一次，每当原始所有者需要从拥有粒子所在单元的进程中获得信息时，就会发生一对一的通信。我们确保使用初始粒子分布来设置粒子的所有权，并在程序的整个执行过程中保持相同的所有权。

  // 有了这个概述，让我们看看这个函数做什么。在顶部，我们创建了一个临时的三角形和DoFHandler对象，我们将从中获取初始粒子位置的节点位置。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::setup_tracer_particles()
  {
    parallel::distributed::Triangulation<spacedim> particle_insert_tria(
      mpi_communicator);
    GridGenerator::generate_from_name_and_arguments(
      particle_insert_tria,
      par.name_of_particle_grid,
      par.arguments_for_particle_grid);
    particle_insert_tria.refine_global(par.particle_insertion_refinement);

    FE_Q<spacedim>       particles_fe(1);
    DoFHandler<spacedim> particles_dof_handler(particle_insert_tria);
    particles_dof_handler.distribute_dofs(particles_fe);

    // 这就是事情开始变得复杂的地方。由于我们可能会在并行环境中运行这个程序，每个并行进程现在都会创建这些临时三角形和DoFHandlers。但是，在完全分布式三角形中，活动进程只知道本地拥有的单元，而不知道其他进程是如何分布自己的单元的。这对于上面创建的临时三角形以及我们想嵌入粒子的流体三角形都是如此。另一方面，一般来说，这两个三角形的局部已知部分不会重合。也就是说，我们将从临时网格的节点位置创建的粒子的位置是任意的，并且可能落在当前进程无法访问的流体三角结构的区域内（即流体领域中细胞是人工的区域）。为了了解将这些粒子发送给谁，我们需要对流体网格在处理器中的分布有一个（粗略的）概念。

    // 我们通过以下方式来构建这一信息：首先建立一个以本地拥有的单元为边界的盒子的索引树，然后提取该树的第一层中的一个。

    std::vector<BoundingBox<spacedim>> all_boxes;
    all_boxes.reserve(fluid_tria.n_locally_owned_active_cells());
    for (const auto &cell : fluid_tria.active_cell_iterators())
      if (cell->is_locally_owned())
        all_boxes.emplace_back(cell->bounding_box());

    const auto tree = pack_rtree(all_boxes);
    const auto local_boxes =
      extract_rtree_level(tree, par.fluid_rtree_extraction_level);

    // 每个进程现在都有一个完全包围所有本地拥有的进程的边界盒集合（但可能与其他进程的边界盒相重叠）。然后我们在所有参与的进程之间交换这些信息，这样每个进程都知道所有其他进程的边界盒。

    // 有了这些信息，我们就可以将`tracer_particle_handler`初始化到流体网格，并从（临时）tracer
    // particles
    // triangulation的支持点生成粒子。这个函数调用使用了我们刚刚构建的`global_bounding_boxes`对象，以确定将位置来自`particles_dof_handler`的本地拥有部分的粒子发送到何处。在这个调用结束时，每个粒子将被分配到正确的进程（即拥有粒子所在的流体单元的进程）。在这一点上，我们也将他们的编号输出到屏幕上。

    global_fluid_bounding_boxes =
      Utilities::MPI::all_gather(mpi_communicator, local_boxes);

    tracer_particle_handler.initialize(fluid_tria,
                                       StaticMappingQ1<spacedim>::mapping);

    Particles::Generators::dof_support_points(particles_dof_handler,
                                              global_fluid_bounding_boxes,
                                              tracer_particle_handler);

    pcout << "Tracer particles: "
          << tracer_particle_handler.n_global_particles() << std::endl;

    // 这样创建的每个粒子都有一个唯一的ID。在下面的算法中的某个时刻，我们将需要包含每个粒子的位置和速度信息的向量。这个向量的大小为`n_particles
    // * // spacedim`，我们需要为每个粒子提供位置和速度信息。
    // spacedim`，我们将不得不以一种方式来存储这个向量的元素，以便每个并行进程
    // "拥有
    // "与它拥有的粒子的坐标相对应的那些元素。换句话说，我们必须在所有进程中划分0和`n_particles
    // * spacedim`之间的索引空间。我们可以通过查询`tracer_particle_handler`的本地相关粒子的ID来做到这一点，并构建需要的索引，将所有粒子的位置和速度存储在一个（平行分布的）矢量中，其中我们隐含地假设我们将每个位置或速度的坐标存储在`spacedim`连续的矢量元素中（这就是 IndexSet::tensor_priduct() 函数的作用）。

    locally_owned_tracer_particle_coordinates =
      tracer_particle_handler.locally_owned_particle_ids().tensor_product(
        complete_index_set(spacedim));

    // 在模拟开始时，所有粒子都在它们的原始位置。当粒子移动时，它们可能会穿越到另一个进程所拥有的领域的某个部分。如果发生这种情况，当前进程会正式保持对粒子的
    // "所有权"，但可能需要从粒子落地的进程中读取访问。我们将这一信息保存在另一个索引集中，该索引集存储了当前进程的子域中的所有粒子的索引，不管它们是否一直在这里。

    // 保留这个索引集使我们能够利用线性代数类来进行有关粒子位置和速度的所有通信。这模拟了在固体域中解决另一个问题的情况下会发生的情况（如在流体-结构相互作用中。在后一种情况下，实体域上的额外DOFs将被耦合到流体域中发生的情况。

    locally_relevant_tracer_particle_coordinates =
      locally_owned_tracer_particle_coordinates;

    // 最后，我们要确保在细化时，粒子被正确转移。在进行局部细化或粗化时，粒子会落在另一个单元中。原则上，我们可以在细化后重新分配所有的粒子，但是这将是非常昂贵的。

    //  Particles::ParticleHandler
    //  类有一种方法可以在细化时将信息从一个单元转移到它的子单元或它的父单元，而不需要重构整个数据结构。这是通过向三角结构注册两个回调函数来实现的。这些函数将在细化即将发生和刚刚发生时收到一个信号，并将以最小的计算成本将所有信息转移到新的细化网格中。

    fluid_tria.signals.pre_distributed_refinement.connect(
      [&]() { tracer_particle_handler.register_store_callback_function(); });

    fluid_tria.signals.post_distributed_refinement.connect([&]() {
      tracer_particle_handler.register_load_callback_function(false);
    });
  }

  // 与我们对被动追踪器所做的类似，我们接下来设置追踪实体网格的正交点的粒子。这里的主要区别是，我们还想给每个粒子附加一个权重值（正交点的
  // "JxW "值），这样我们就可以在不直接访问原始实体网格的情况下计算积分。

  // 这是通过利用 Particles::Particle 类的 "属性
  // "概念实现的。它可以（以一种有效的内存方式）在一个
  // Particles::ParticleHandler 对象内为每个 Particles::Particle
  // 对象存储任意数量的`双`数字。我们利用这种可能性来存储实体网格的正交点的JxW值。

  // 在我们的例子中，我们只需要为每个粒子存储一个属性：实体网格上的积分的JxW值。这将在构造时作为最后一个参数传递给solid_particle_handler对象。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::setup_solid_particles()
  {
    QGauss<dim> quadrature(fluid_fe->degree + 1);

    const unsigned int n_properties = 1;
    solid_particle_handler.initialize(fluid_tria,
                                      StaticMappingQ1<spacedim>::mapping,
                                      n_properties);

    // 我们在本地生成的粒子数等于本地拥有的单元总数乘以每个单元中使用的正交点的数量。我们将所有这些点存储在一个向量中，并将其相应的属性存储在一个向量的向量中。

    std::vector<Point<spacedim>> quadrature_points_vec;
    quadrature_points_vec.reserve(quadrature.size() *
                                  solid_tria.n_locally_owned_active_cells());

    std::vector<std::vector<double>> properties;
    properties.reserve(quadrature.size() *
                       solid_tria.n_locally_owned_active_cells());

    FEValues<dim, spacedim> fe_v(*solid_fe,
                                 quadrature,
                                 update_JxW_values | update_quadrature_points);
    for (const auto &cell : solid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          fe_v.reinit(cell);
          const auto &points = fe_v.get_quadrature_points();
          const auto &JxW    = fe_v.get_JxW_values();

          for (unsigned int q = 0; q < points.size(); ++q)
            {
              quadrature_points_vec.emplace_back(points[q]);
              properties.emplace_back(
                std::vector<double>(n_properties, JxW[q]));
            }
        }

    // 我们以处理示踪粒子的同样方式进行，重新使用计算出的边界盒。然而，我们首先检查`global_fluid_bounding_boxes`对象是否已经被填充。这里当然应该是这样的，因为这个方法是在初始化示踪粒子的方法之后调用的。然而，我们要确保，如果将来有人决定（无论出于什么原因）先初始化固体粒子处理程序，或者只复制教程的这一部分，当事情没有按照预期进行时，会抛出一个有意义的异常。

    // 由于我们已经存储了正交点的位置，我们可以使用这些位置来直接使用`solid_particle_handler`插入粒子，而不必通过
    // Particles::Generators 函数。

    Assert(!global_fluid_bounding_boxes.empty(),
           ExcInternalError(
             "I was expecting the "
             "global_fluid_bounding_boxes to be filled at this stage. "
             "Make sure you fill this vector before trying to use it "
             "here. Bailing out."));

    solid_particle_handler.insert_global_particles(quadrature_points_vec,
                                                   global_fluid_bounding_boxes,
                                                   properties);

    // 和前面的函数一样，我们最后要确保在细化时，粒子被正确转移。

    fluid_tria.signals.pre_distributed_refinement.connect(
      [&]() { solid_particle_handler.register_store_callback_function(); });

    fluid_tria.signals.post_distributed_refinement.connect(
      [&]() { solid_particle_handler.register_load_callback_function(false); });

    pcout << "Solid particles: " << solid_particle_handler.n_global_particles()
          << std::endl;
  }

  //  @sect4{DoF initialization functions}

  // 我们设置了有限元空间和整个步骤中使用的正交公式。对于流体，我们使用Taylor-Hood元素（例如
  // $Q_k \times Q_{k-1}$
  // ）。由于我们没有解决固体领域的任何方程，所以产生了一个空的有限元空间。这个程序的一个自然扩展是解决流体结构的相互作用问题，这就要求`solid_fe`使用更有用的FiniteElement类。

  // 和其他许多函数一样，我们在这里存储了进行操作所需的时间。当前的函数把它的时间信息放到一个标签为
  // "初始设置
  // "的部分。在不同的函数中对这个定时器进行了许多其他的调用。它们允许监测每个单独函数的绝对和相对成本，以确定瓶颈。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::initial_setup()
  {
    TimerOutput::Scope t(computing_timer, "Initial setup");

    fluid_fe =
      std::make_unique<FESystem<spacedim>>(FE_Q<spacedim>(par.velocity_degree),
                                           spacedim,
                                           FE_Q<spacedim>(par.velocity_degree -
                                                          1),
                                           1);

    solid_fe = std::make_unique<FE_Nothing<dim, spacedim>>();
    solid_dh.distribute_dofs(*solid_fe);

    fluid_quadrature_formula =
      std::make_unique<QGauss<spacedim>>(par.velocity_degree + 1);
    solid_quadrature_formula =
      std::make_unique<QGauss<dim>>(par.velocity_degree + 1);
  }

  // 接下来我们构建分布式块状矩阵和向量，用于解决问题中出现的线性方程。这个函数改编自
  // step-55 ，我们参考这个步骤进行全面解释。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::setup_dofs()
  {
    TimerOutput::Scope t(computing_timer, "Setup dofs");

    fluid_dh.distribute_dofs(*fluid_fe);

    std::vector<unsigned int> stokes_sub_blocks(spacedim + 1, 0);
    stokes_sub_blocks[spacedim] = 1;
    DoFRenumbering::component_wise(fluid_dh, stokes_sub_blocks);

    auto dofs_per_block =
      DoFTools::count_dofs_per_fe_block(fluid_dh, stokes_sub_blocks);

    const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

    pcout << "   Number of degrees of freedom: " << fluid_dh.n_dofs() << " ("
          << n_u << '+' << n_p << " -- "
          << solid_particle_handler.n_global_particles() << '+'
          << tracer_particle_handler.n_global_particles() << ')' << std::endl;

    fluid_owned_dofs.resize(2);
    fluid_owned_dofs[0] = fluid_dh.locally_owned_dofs().get_view(0, n_u);
    fluid_owned_dofs[1] =
      fluid_dh.locally_owned_dofs().get_view(n_u, n_u + n_p);

    IndexSet locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);
    fluid_relevant_dofs.resize(2);
    fluid_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    fluid_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    {
      constraints.reinit(locally_relevant_dofs);

      FEValuesExtractors::Vector velocities(0);
      DoFTools::make_hanging_node_constraints(fluid_dh, constraints);
      VectorTools::interpolate_boundary_values(
        fluid_dh,
        0,
        Functions::ZeroFunction<spacedim>(spacedim + 1),
        constraints,
        fluid_fe->component_mask(velocities));
      constraints.close();
    }

    auto locally_owned_dofs_per_processor =
      Utilities::MPI::all_gather(mpi_communicator,
                                 fluid_dh.locally_owned_dofs());
    {
      system_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::none;
          else if (c == spacedim || d == spacedim || c == d)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);

      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);

      system_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    {
      preconditioner_matrix.clear();

      Table<2, DoFTools::Coupling> coupling(spacedim + 1, spacedim + 1);
      for (unsigned int c = 0; c < spacedim + 1; ++c)
        for (unsigned int d = 0; d < spacedim + 1; ++d)
          if (c == spacedim && d == spacedim)
            coupling[c][d] = DoFTools::always;
          else
            coupling[c][d] = DoFTools::none;

      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

      DoFTools::make_sparsity_pattern(
        fluid_dh, coupling, dsp, constraints, false);
      SparsityTools::distribute_sparsity_pattern(
        dsp,
        locally_owned_dofs_per_processor,
        mpi_communicator,
        locally_relevant_dofs);
      preconditioner_matrix.reinit(fluid_owned_dofs, dsp, mpi_communicator);
    }

    locally_relevant_solution.reinit(fluid_owned_dofs,
                                     fluid_relevant_dofs,
                                     mpi_communicator);
    system_rhs.reinit(fluid_owned_dofs, mpi_communicator);
    solution.reinit(fluid_owned_dofs, mpi_communicator);
  }
  // @sect4{Assembly functions}

  // 我们将系统矩阵、预处理矩阵和右手边组合起来。这段代码改编自 step-55
  // ，基本上是 step-27
  // 的内容，如果你知道斯托克斯方程是什么样子的，就会觉得很标准。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::assemble_stokes_system()
  {
    system_matrix         = 0;
    preconditioner_matrix = 0;
    system_rhs            = 0;

    TimerOutput::Scope t(computing_timer, "Assemble Stokes terms");

    FEValues<spacedim> fe_values(*fluid_fe,
                                 *fluid_quadrature_formula,
                                 update_values | update_gradients |
                                   update_quadrature_points |
                                   update_JxW_values);

    const unsigned int dofs_per_cell = fluid_fe->n_dofs_per_cell();
    const unsigned int n_q_points    = fluid_quadrature_formula->size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_matrix2(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(spacedim + 1));

    std::vector<Tensor<2, spacedim>> grad_phi_u(dofs_per_cell);
    std::vector<double>              div_phi_u(dofs_per_cell);
    std::vector<double>              phi_p(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    const FEValuesExtractors::Vector     velocities(0);
    const FEValuesExtractors::Scalar     pressure(spacedim);

    for (const auto &cell : fluid_dh.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          cell_matrix  = 0;
          cell_matrix2 = 0;
          cell_rhs     = 0;

          fe_values.reinit(cell);
          par.rhs.vector_value_list(fe_values.get_quadrature_points(),
                                    rhs_values);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                  grad_phi_u[k] = fe_values[velocities].gradient(k, q);
                  div_phi_u[k]  = fe_values[velocities].divergence(k, q);
                  phi_p[k]      = fe_values[pressure].value(k, q);
                }

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    {
                      cell_matrix(i, j) +=
                        (par.viscosity *
                           scalar_product(grad_phi_u[i], grad_phi_u[j]) -
                         div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) *
                        fe_values.JxW(q);

                      cell_matrix2(i, j) += 1.0 / par.viscosity * phi_p[i] *
                                            phi_p[j] * fe_values.JxW(q);
                    }

                  const unsigned int component_i =
                    fluid_fe->system_to_component_index(i).first;
                  cell_rhs(i) += fe_values.shape_value(i, q) *
                                 rhs_values[q](component_i) * fe_values.JxW(q);
                }
            }

          cell->get_dof_indices(local_dof_indices);
          constraints.distribute_local_to_global(cell_matrix, 
                                                 cell_rhs, 
                                                 local_dof_indices, 
                                                 system_matrix, 
 

          constraints.distribute_local_to_global(cell_matrix2, 
                                                 local_dof_indices, 
                                                 preconditioner_matrix);
        }

    system_matrix.compress(VectorOperation::add);
    preconditioner_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }

  // 下面的方法是处理因对叶轮施加速度而产生的惩罚项。从某种意义上说，它是本教程的核心，但它相对简单。这里我们利用`solid_particle_handler`来计算Nitsche限制或嵌入域中的惩罚。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::assemble_nitsche_restriction()
  {
    TimerOutput::Scope t(computing_timer, "Assemble Nitsche terms");

    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(spacedim);

    SolidVelocity<spacedim> solid_velocity(par.angular_velocity);

    std::vector<types::global_dof_index> fluid_dof_indices(
      fluid_fe->n_dofs_per_cell());

    FullMatrix<double>     local_matrix(fluid_fe->n_dofs_per_cell(),
                                    fluid_fe->n_dofs_per_cell());
    dealii::Vector<double> local_rhs(fluid_fe->n_dofs_per_cell());

    const auto penalty_parameter =
      1.0 / GridTools::minimal_cell_diameter(fluid_tria);

    // 我们在所有的本地粒子上循环。虽然这可以直接通过循环所有的单元格来实现，但这将迫使我们循环许多不包含粒子的单元格。因此，我们在所有的粒子上循环，但是，我们得到粒子所在的单元格的参考，然后在该单元格中循环所有的粒子。这使得我们能够跳过不包含粒子的单元格，但又能集合每个单元格的局部矩阵和rhs来应用Nitsche的限制。一旦我们完成了一个单元格上的所有粒子，我们就将`粒子`迭代器推进到当前单元格上的粒子的末端（这是`while`循环体的最后一行）。

    auto particle = solid_particle_handler.begin();
    while (particle != solid_particle_handler.end())
      {
        local_matrix = 0;
        local_rhs    = 0;

        // 我们从粒子本身得到一个通往粒子所在单元的迭代器。然后，我们就可以像通常那样在系统矩阵和右手边组装附加项了。

        const auto &cell = particle->get_surrounding_cell(fluid_tria);
        const auto &dh_cell =
          typename DoFHandler<spacedim>::cell_iterator(*cell, &fluid_dh);
        dh_cell->get_dof_indices(fluid_dof_indices);
        //所以
        //然后让我们得到位于这个单元格上的单元格集合，并对它们进行迭代。从每个粒子中，我们收集该粒子的位置和参考位置，以及附加在该粒子上的额外信息。在本例中，这些信息是用于生成粒子的正交点的
        //"JxW"。

        // 利用这些信息，我们可以将正交点的贡献加入到local_matrix和local_rhs中。我们可以利用每个粒子的参考位置，轻松地评估其位置上的形状函数值。

        const auto pic = solid_particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (const auto &p : pic)
          {
            const auto &ref_q  = p.get_reference_location();
            const auto &real_q = p.get_location();
            const auto &JxW    = p.get_properties()[0];

            for (unsigned int i = 0; i < fluid_fe->n_dofs_per_cell(); ++i)
              {
                const auto comp_i =
                  fluid_fe->system_to_component_index(i).first;
                if (comp_i < spacedim)
                  {
                    for (unsigned int j = 0; j < fluid_fe->n_dofs_per_cell();
                         ++j)
                      {
                        const auto comp_j =
                          fluid_fe->system_to_component_index(j).first;
                        if (comp_i == comp_j)
                          local_matrix(i, j) +=
                            penalty_parameter * par.penalty_term *
                            fluid_fe->shape_value(i, ref_q) *
                            fluid_fe->shape_value(j, ref_q) * JxW;
                      }
                    local_rhs(i) += penalty_parameter * par.penalty_term *
                                    solid_velocity.value(real_q, comp_i) *
                                    fluid_fe->shape_value(i, ref_q) * JxW;
                  }
              }
          }

        constraints.distribute_local_to_global(local_matrix,
                                               local_rhs,
                                               fluid_dof_indices,
                                               system_matrix,
                                               system_rhs);
        particle = pic.end();
      }

    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);
  }
  // @sect4{Solving the linear system}

  // 这个函数用FGMRES求解线性系统，有一个对角线块的预处理和一个对角线块的代数多重网格（AMG）方法。该预处理程序对
  // $(0,0)$ （即速度-速度）块应用V循环，对 $(1,1)$
  // 块应用质量矩阵的CG（这是我们对舒尔补码的近似值：上面组装的压力质量矩阵）。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::solve()
  {
    TimerOutput::Scope t(computing_timer, "Solve");

    LA::MPI::PreconditionAMG prec_A;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_A.initialize(system_matrix.block(0, 0), data);
    }

    LA::MPI::PreconditionAMG prec_S;
    {
      LA::MPI::PreconditionAMG::AdditionalData data;

#ifdef USE_PETSC_LA
      data.symmetric_operator = true;
#endif
      prec_S.initialize(preconditioner_matrix.block(1, 1), data);
    }

    const auto A = linear_operator<LA::MPI::Vector>(system_matrix.block(0, 0));
    const auto amgA = linear_operator(A, prec_A);

    const auto S =
      linear_operator<LA::MPI::Vector>(preconditioner_matrix.block(1, 1));
    const auto amgS = linear_operator(S, prec_S);

    ReductionControl          inner_solver_control(100,
                                          1e-8 * system_rhs.l2_norm(),
                                          1.e-2);
    SolverCG<LA::MPI::Vector> cg(inner_solver_control);

    const auto invS = inverse_operator(S, cg, amgS);

    const auto P = block_diagonal_operator<2, LA::MPI::BlockVector>(
      std::array<
        dealii::LinearOperator<typename LA::MPI::BlockVector::BlockType>,
        2>{{amgA, amgS}});

    SolverControl solver_control(system_matrix.m(),
                                 1e-10 * system_rhs.l2_norm());

    SolverFGMRES<LA::MPI::BlockVector> solver(solver_control);

    constraints.set_zero(solution);



    pcout << "   Solved in " << solver_control.last_step() << " iterations."
          << std::endl;

    constraints.distribute(solution);

    locally_relevant_solution = solution;
    const double mean_pressure =
      VectorTools::compute_mean_value(fluid_dh,
                                      QGauss<spacedim>(par.velocity_degree + 2),
                                      locally_relevant_solution,
                                      spacedim);
    solution.block(1).add(-mean_pressure);
    locally_relevant_solution.block(1) = solution.block(1);
  }

  //  @sect4{Mesh refinement}

  // 我们以一种完全标准的方式处理网格细化问题。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::refine_and_transfer()
  {
    TimerOutput::Scope               t(computing_timer, "Refine");
    const FEValuesExtractors::Vector velocity(0);

    Vector<float> error_per_cell(fluid_tria.n_active_cells());
    KellyErrorEstimator<spacedim>::estimate(fluid_dh,
                                            QGauss<spacedim - 1>(
                                              par.velocity_degree + 1),
                                            {},
                                            locally_relevant_solution,
                                            error_per_cell,
                                            fluid_fe->component_mask(velocity));

    if (par.refinement_strategy == "fixed_fraction")
      {
        parallel::distributed::GridRefinement::
          refine_and_coarsen_fixed_fraction(fluid_tria,
                                            error_per_cell,
                                            par.refinement_fraction,
                                            par.coarsening_fraction);
      }
    else if (par.refinement_strategy == "fixed_number")
      {
        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
          fluid_tria,
          error_per_cell,
          par.refinement_fraction,
          par.coarsening_fraction,
          par.max_cells);
      }

    for (const auto &cell : fluid_tria.active_cell_iterators())
      {
        if (cell->refine_flag_set() &&
            cell->level() == par.max_level_refinement)
          cell->clear_refine_flag();
        if (cell->coarsen_flag_set() &&
            cell->level() == par.min_level_refinement)
          cell->clear_coarsen_flag();
      }

    parallel::distributed::SolutionTransfer<spacedim, LA::MPI::BlockVector>
      transfer(fluid_dh);
    fluid_tria.prepare_coarsening_and_refinement();
    transfer.prepare_for_coarsening_and_refinement(locally_relevant_solution);
    fluid_tria.execute_coarsening_and_refinement();

    setup_dofs();

    transfer.interpolate(solution);
    constraints.distribute(solution);
    locally_relevant_solution = solution;
  }
  // @sect4{Creating output for visualization}

  // 我们使用deal.II的标准并行功能在流体域上输出结果（速度和压力）。编写一个压缩的vtu文件，将所有处理器的信息聚集在一起。另外写一个`.pvd`记录，将物理时间与vtu文件联系起来。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::output_results(const unsigned int cycle,
                                                       double time) const
  {
    TimerOutput::Scope t(computing_timer, "Output fluid");

    std::vector<std::string> solution_names(spacedim, "velocity");
    solution_names.emplace_back("pressure");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        spacedim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    DataOut<spacedim> data_out;
    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(locally_relevant_solution,
                             solution_names,
                             DataOut<spacedim>::type_dof_data,
                             data_component_interpretation);

    Vector<float> subdomain(fluid_tria.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = fluid_tria.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches();

    const std::string filename =
      "solution-" + Utilities::int_to_string(cycle) + ".vtu";
    data_out.write_vtu_in_parallel(par.output_directory + "/" + filename,
                                   mpi_communicator);

    static std::vector<std::pair<double, std::string>> times_and_names;
    times_and_names.push_back(std::make_pair(time, filename));
    std::ofstream ofile(par.output_directory + "/" + "solution.pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names);
  }

  // 同样地，我们通过 Particles::DataOut
  // 对象将粒子（无论是来自实体还是追踪器）写成一个单一的压缩vtu文件。这个简单的对象并不写作为
  // "属性
  // "附加到粒子上的额外信息，而只写它们的id--但是，无论如何，我们并不关心这些粒子位置的
  // "JxW "值，所以我们可能想要可视化的信息并没有丢失。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::output_particles(
    const Particles::ParticleHandler<spacedim> &particles,
    std::string                                 fprefix,
    const unsigned int                          iter,
    const double                                time) const
  {
    Particles::DataOut<spacedim> particles_out;
    particles_out.build_patches(particles);
    const std::string filename =
      (fprefix + "-" + Utilities::int_to_string(iter) + ".vtu");
    particles_out.write_vtu_in_parallel(par.output_directory + "/" + filename,
                                        mpi_communicator);

    static std::map<std::string, std::vector<std::pair<double, std::string>>>
      times_and_names;
    if (times_and_names.find(fprefix) != times_and_names.end())
      times_and_names[fprefix].push_back(std::make_pair(time, filename));
    else
      times_and_names[fprefix] = {std::make_pair(time, filename)};
    std::ofstream ofile(par.output_directory + "/" + fprefix + ".pvd");
    DataOutBase::write_pvd_record(ofile, times_and_names[fprefix]);
  }
  // @sect4{The "run" function}

  // 这个函数现在负责协调整个模拟过程。它与其他时间相关的教程程序非常相似--以
  // step-21 或 step-26
  // 为例。在开始的时候，我们会输出一些状态信息，同时将所有的当前参数保存到输出目录下的文件中，以利于重现。

  template <int dim, int spacedim>
  void
  StokesImmersedProblem<dim, spacedim>::run()
  {
#ifdef USE_PETSC_LA
    pcout << "Running StokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using PETSc."
          << std::endl;
#else
    pcout << "Running StokesImmersedProblem<"
          << Utilities::dim_string(dim, spacedim) << "> using Trilinos."
          << std::endl;
#endif
    par.prm.print_parameters(par.output_directory + "/" + "used_parameters_" +
                               std::to_string(dim) + std::to_string(spacedim) +
                               ".prm",
                             ParameterHandler::Short);

    // 然后我们开始时间循环。我们在第一个循环中初始化模拟的所有元素

    const double time_step    = par.final_time / (par.number_of_time_steps - 1);
    double       time         = 0;
    unsigned int output_cycle = 0;

    for (unsigned int cycle = 0; cycle < par.number_of_time_steps;
         ++cycle, time += time_step)
      {
        par.set_time(time);
        pcout << "Cycle " << cycle << ':' << std::endl
              << "Time : " << time << ", time step: " << time_step << std::endl;

        if (cycle == 0)
          {
            make_grid();
            initial_setup();
            setup_dofs();
            setup_tracer_particles();
            setup_solid_particles();
            tracer_particle_velocities.reinit(
              locally_owned_tracer_particle_coordinates, mpi_communicator);
            output_results(output_cycle, time);
            {
              TimerOutput::Scope t(computing_timer, "Output tracer particles");
              output_particles(tracer_particle_handler,
                               "tracer",
                               output_cycle,
                               time);
            }
            {
              TimerOutput::Scope t(computing_timer, "Output solid particles");
              output_particles(solid_particle_handler,
                               "solid",
                               output_cycle,
                               time);
            }
          }

        // 在第一个时间步长之后，我们在每个时间步长的开始时对实体进行位移，以考虑到它已经移动的事实。

        else
          {
            TimerOutput::Scope t(computing_timer,
                                 "Set solid particle position");

            SolidPosition<spacedim> solid_position(par.angular_velocity,
                                                   time_step);
            solid_particle_handler.set_particle_positions(solid_position,
                                                          false);
          }

        // 为了更新系统的状态，我们首先对示踪粒子位置的流体速度进行插值，并采用天真的显式欧拉方案对无质量示踪粒子进行漂移。

        {
          TimerOutput::Scope t(computing_timer, "Set tracer particle motion");
          Particles::Utilities::interpolate_field_on_particles(
            fluid_dh,
            tracer_particle_handler,
            locally_relevant_solution,
            tracer_particle_velocities,
            fluid_fe->component_mask(FEValuesExtractors::Vector(0)));

          tracer_particle_velocities *= time_step;

          locally_relevant_tracer_particle_coordinates =
            tracer_particle_handler.locally_owned_particle_ids().tensor_product(
              complete_index_set(spacedim));

          relevant_tracer_particle_displacements.reinit(
            locally_owned_tracer_particle_coordinates,
            locally_relevant_tracer_particle_coordinates,
            mpi_communicator);

          relevant_tracer_particle_displacements = tracer_particle_velocities;

          tracer_particle_handler.set_particle_positions(
            relevant_tracer_particle_displacements);
        }

        // 利用这些新的位置，我们就可以组装斯托克斯系统，并解决它。

        assemble_stokes_system();
        assemble_nitsche_restriction();
        solve();

        // 在适当的频率下，我们再将固体粒子、示踪粒子和流体领域的信息写入文件，以便进行可视化，并通过适应网格来结束时间步骤。

        if (cycle % par.output_frequency == 0)
          {
            output_results(output_cycle, time);
            {
              TimerOutput::Scope t(computing_timer, "Output tracer particles");
              output_particles(tracer_particle_handler,
                               "tracer",
                               output_cycle,
                               time);
            }
            {
              TimerOutput::Scope t(computing_timer, "Output solid particles");
              output_particles(solid_particle_handler,
                               "solid",
                               output_cycle,
                               time);
            }
            ++output_cycle;
          }
        if (cycle % par.refinement_frequency == 0 &&
            cycle != par.number_of_time_steps - 1)
          refine_and_transfer();
      }
  }

} // namespace Step70
// @sect3{The main() function}

// 代码的其余部分，即`main()`函数，是标准的，除了对输入参数文件的处理。我们允许用户指定一个可选的参数文件作为程序的参数。如果没有指定，我们就使用默认文件
// "parameters.prm"，如果不存在，我们就创建这个文件。文件名首先被扫描为字符串
// "23"，然后是 "3"。如果文件名包含字符串
// "23"，问题类将分别以模板参数2和3进行实例化。如果只找到 "3
// "这个字符串，那么两个模板参数都被设置为3，否则都被设置为2。

// 如果程序被调用时没有任何命令行参数（即`argc==1`），那么我们就默认使用
// "参数.prm"。

int
main(int argc, char *argv[])
{
  using namespace Step70;
  using namespace dealii;
  deallog.depth_console(1);
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

      std::string prm_file;
      if (argc > 1)
        prm_file = argv[1];
      else
        prm_file = "parameters.prm";

      if (prm_file.find("23") != std::string::npos)
        {
          StokesImmersedProblemParameters<2, 3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2, 3> problem(par);
          problem.run();
        }
      else if (prm_file.find("3") != std::string::npos)
        {
          StokesImmersedProblemParameters<3> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<3> problem(par);
          problem.run();
        }
      else
        {
          StokesImmersedProblemParameters<2> par;
          ParameterAcceptor::initialize(prm_file);

          StokesImmersedProblem<2> problem(par);
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
