

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
 * Authors: Bruno Blais, Toni El Geitani Nehme, Rene Gassmoeller, Peter Munch
 */


// @sect3{Include files}
#include <deal.II/base/bounding_box.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/discrete_time.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parameter_acceptor.h>
#include <deal.II/base/timer.h>

#include <deal.II/distributed/cell_weights.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// 从下面的include文件中，我们导入了ParticleHandler类，该类允许你管理漂浮在
// Particles::Particle),
// 类型的粒子集合（代表具有一些附加属性（例如，一个id）的点集合的对象）。
// Particles命名空间中的方法和类允许人们轻松实现Particle-In-Cell方法和分布式三角形上的粒子追踪。

#include <deal.II/particles/particle_handler.h>

// 我们导入粒子发生器，使我们能够插入粒子。在本步骤中，粒子是通过非匹配的超壳三角形全局插入的。

#include <deal.II/particles/generators.h>

// 由于粒子没有形成三角形，它们有自己特定的DataOut类，这将使我们能够把它们写成常用的并行vtu格式（或其他任何数量的文件格式）。

#include <deal.II/particles/data_out.h>

#include <cmath>
#include <iostream>

namespace Step68
{
  using namespace dealii;
  // @sect3{Run-time parameter handling}

  // 与 step-60
  // 中的做法类似，我们建立了一个持有我们问题的所有参数的类，并从ParameterAcceptor类中派生出来以简化参数文件的管理和创建。

  // ParameterAcceptor范式要求所有的参数都可以被ParameterAcceptor方法写入。为了避免出现很难追踪的bug（比如写成`if
  // (time = 0)`而不是`if(time ==
  // 0)`），我们在一个外部类中声明所有的参数，该类在实际的`ParticleTracking`类之前被初始化，并将其作为`const`引用传递给主类。

  // 该类的构造函数负责该类的成员与ParameterHandler中的相应条目之间的连接。由于使用了
  // ParameterHandler::add_parameter()
  // 方法，这种连接是微不足道的，但要求这个类的所有成员都是可写的。

  class ParticleTrackingParameters : public ParameterAcceptor
  {
  public:
    ParticleTrackingParameters();

    // 该类主要由成员变量组成，描述了粒子跟踪模拟及其离散化的细节。下面的参数是关于输出应该写到哪里，速度的空间离散化（默认是
    // $Q_1$
    // ），时间步长和输出频率（在我们再次生成图形输出之前应该经过多少时间步长）。

    std::string output_directory = "./";

    unsigned int velocity_degree       = 1;
    double       time_step             = 0.002;
    double       final_time            = 4.0;
    unsigned int output_frequency      = 10;
    unsigned int repartition_frequency = 5;

    // 我们允许每个网格独立地被细化。在本教程中，流体网格上没有解决物理问题，其速度是通过分析计算得出的。

    unsigned int fluid_refinement              = 4;
    unsigned int particle_insertion_refinement = 3;
  };

  // 还有一个任务就是声明我们在输入文件中可以接受哪些运行时参数。由于我们的参数数量非常有限，所有的参数都在同一章节中声明。

  ParticleTrackingParameters::ParticleTrackingParameters()
    : ParameterAcceptor("Particle Tracking Problem/")
  {
    add_parameter(
      "Velocity degree", velocity_degree, "", prm, Patterns::Integer(1));

    add_parameter("Output frequency",
                  output_frequency,
                  "Iteration frequency at which output results are written",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Repartition frequency",
                  repartition_frequency,
                  "Iteration frequency at which the mesh is load balanced",
                  prm,
                  Patterns::Integer(1));

    add_parameter("Output directory", output_directory);

    add_parameter("Time step", time_step, "", prm, Patterns::Double());

    add_parameter("Final time",
                  final_time,
                  "End time of the simulation",
                  prm,
                  Patterns::Double());

    add_parameter("Fluid refinement",
                  fluid_refinement,
                  "Refinement level of the fluid domain",
                  prm,
                  Patterns::Integer(0));

    add_parameter(
      "Particle insertion refinement",
      particle_insertion_refinement,
      "Refinement of the volumetric mesh used to insert the particles",
      prm,
      Patterns::Integer(0));
  }

  //  @sect3{Velocity profile}

  // 速度曲线是作为一个函数对象提供的。这个函数在例子中是硬编码的。

  template <int dim>
  class Vortex : public Function<dim>
  {
  public:
    Vortex()
      : Function<dim>(dim)
    {}

    virtual void
    vector_value(const Point<dim> &point,
                 Vector<double> &  values) const override;
  };

  // Rayleigh-Kothe顶点的速度曲线是随时间变化的。因此，必须从函数对象中收集模拟的当前时间（t）。

  template <int dim>
  void
  Vortex<dim>::vector_value(const Point<dim> &point,
                            Vector<double> &  values) const
  {
    const double T = 4;
    const double t = this->get_time();

    const double px = numbers::PI * point(0);
    const double py = numbers::PI * point(1);
    const double pt = numbers::PI / T * t;

    values[0] = -2 * cos(pt) * pow(sin(px), 2) * sin(py) * cos(py);
    values[1] = 2 * cos(pt) * pow(sin(py), 2) * sin(px) * cos(px);
    if (dim == 3)
      {
        values[2] = 0;
      }
  }

  //  @sect3{The <code>ParticleTracking</code> class declaration}

  // 我们现在准备介绍我们的教程程序的主类。

  template <int dim>
  class ParticleTracking
  {
  public:
    ParticleTracking(const ParticleTrackingParameters &par,
                     const bool                        interpolated_velocity);
    void
    run();

  private:
    // 这个函数负责在背景网格之上初始生成粒子。

    void
    generate_particles();

    // 当速度曲线被内插到粒子的位置时，必须首先使用自由度来存储。因此，和其他并行情况一样（例如
    // step-40 ），我们在背景网格上初始化自由度。

    void
    setup_background_dofs();

    // 在其中一个测试案例中，该函数被映射到背景网格上，并使用有限元插值来计算粒子位置的速度。这个函数计算三角形的支持点处的函数值。

    void
    interpolate_function_to_field();

    // 下面两个函数分别负责对速度场在粒子位置插值或分析计算的情况下进行显式欧拉时间积分的步骤。

    void
    euler_step_interpolated(const double dt);
    void
    euler_step_analytical(const double dt);

    // `cell_weight()`函数向三角计算表明在这个单元上预计会发生多少计算工作，因此需要对域进行划分，以使每个MPI等级得到大致相等的工作量（可能不是相等的单元数量）。虽然该函数是从外部调用的，但它与该类内部的相应信号相连，因此它可以是
    // "私有 "的。

    unsigned int
    cell_weight(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const typename parallel::distributed::Triangulation<dim>::CellStatus
        status) const;

    // 以下两个函数分别负责输出粒子的模拟结果和背景网格上的速度曲线。

    void
    output_particles(const unsigned int it);
    void
    output_background(const unsigned int it);

    // 该类的私有成员与其他并行deal.II例子相似。参数被存储为`const`成员。值得注意的是，我们保留了`Vortex`类的成员，因为它的时间必须随着模拟的进行而被修改。

    const ParticleTrackingParameters &par;

    MPI_Comm                                  mpi_communicator;
    parallel::distributed::Triangulation<dim> background_triangulation;
    Particles::ParticleHandler<dim>           particle_handler;

    DoFHandler<dim>                            fluid_dh;
    FESystem<dim>                              fluid_fe;
    MappingQ1<dim>                             mapping;
    LinearAlgebra::distributed::Vector<double> velocity_field;

    Vortex<dim> velocity;

    ConditionalOStream pcout;

    bool interpolated_velocity;
  };

  //  @sect3{The <code>PatricleTracking</code> class implementation}
  // @sect4{Constructor}

  // 构造函数和析构函数是相当微不足道的。它们与  step-40
  // 中的做法非常相似。我们将我们要工作的处理器设置为所有可用的机器（`MPI_COMM_WORLD`），并初始化
  // <code>pcout</code>  变量，只允许处理器0输出任何东西到标准输出。

  template <int dim>
  ParticleTracking<dim>::ParticleTracking(const ParticleTrackingParameters &par,
                                          const bool interpolated_velocity)
    : par(par)
    , mpi_communicator(MPI_COMM_WORLD)
    , background_triangulation(mpi_communicator)
    , fluid_dh(background_triangulation)
    , fluid_fe(FE_Q<dim>(par.velocity_degree), dim)
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    , interpolated_velocity(interpolated_velocity)

  {}

  //  @sect4{Cell weight}

  // 这个函数是让我们动态平衡本例中计算负载的关键部分。该函数为每个单元赋予一个权重，代表该单元的计算工作。在这里，大部分的工作预计会发生在粒子上，因此这个函数的返回值（代表
  // "这个单元的工作"）是根据当前单元中的粒子数量来计算。该函数与三角形内部的cell_weight()信号相连，每一个单元将被调用一次，每当三角形在等级之间重新划分领域时（该连接是在该类的generate_particles()函数中创建的）。

  template <int dim>
  unsigned int
  ParticleTracking<dim>::cell_weight(
    const typename parallel::distributed::Triangulation<dim>::cell_iterator
      &                                                                  cell,
    const typename parallel::distributed::Triangulation<dim>::CellStatus status)
    const
  {
    // 我们不给我们不拥有的细胞分配任何权重（即人工或幽灵细胞）。

    if (!cell->is_locally_owned())
      return 0;

    // 这决定了粒子工作与细胞工作相比有多重要（默认情况下每个细胞的权重为1000）。我们将每个粒子的权重设置得更高，以表明在这个例子中，粒子的负载是唯一对分配单元很重要的。这个数字的最佳值取决于应用，可以从0（廉价的粒子操作，昂贵的单元操作）到远远大于1000（昂贵的粒子操作，廉价的单元操作，像本例中假定的那样）。

    const unsigned int particle_weight = 10000;

    // 这个例子没有使用自适应细化，因此每个单元都应该有`CELL_PERSIST`的状态。然而这个函数也可以用来在细化过程中分配负载，因此我们也考虑细化或粗化的单元。

    if (status == parallel::distributed::Triangulation<dim>::CELL_PERSIST ||
        status == parallel::distributed::Triangulation<dim>::CELL_REFINE)
      {
        const unsigned int n_particles_in_cell =
          particle_handler.n_particles_in_cell(cell);
        return n_particles_in_cell * particle_weight;
      }
    else if (status == parallel::distributed::Triangulation<dim>::CELL_COARSEN)
      {
        unsigned int n_particles_in_cell = 0;

        for (unsigned int child_index = 0; child_index < cell->n_children();
             ++child_index)
          n_particles_in_cell +=
            particle_handler.n_particles_in_cell(cell->child(child_index));

        return n_particles_in_cell * particle_weight;
      }

    Assert(false, ExcInternalError());
    return 0;
  }

  //  @sect4{Particles generation}

  // 这个函数生成示踪粒子和这些粒子演化的背景三角图。

  template <int dim>
  void
  ParticleTracking<dim>::generate_particles()
  {
    // 我们创建一个超立方体三角形，并对其进行全局细化。这个三角形覆盖了粒子的全部运动轨迹。

    GridGenerator::hyper_cube(background_triangulation, 0, 1);
    background_triangulation.refine_global(par.fluid_refinement);

    // 为了在重新划分三角形时考虑粒子，该算法需要知道三件事。

    // 1.给每个单元分配多少权重（里面有多少粒子）；2.在运送数据之前如何包装粒子；3.在重新分区之后如何拆开粒子。

    // 我们将正确的函数附加到信号里面  parallel::distributed::Triangulation.
    // 这些信号将在每次调用repartition()函数时被调用。这些连接只需要创建一次，所以我们不妨在这个类的构造函数中设置它们，但为了这个例子，我们要把粒子相关的指令分组。

    background_triangulation.signals.cell_weight.connect(
      [&](
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
          &cell,
        const typename parallel::distributed::Triangulation<dim>::CellStatus
          status) -> unsigned int { return this->cell_weight(cell, status); });

    background_triangulation.signals.pre_distributed_repartition.connect(
      [this]() { this->particle_handler.register_store_callback_function(); });

    background_triangulation.signals.post_distributed_repartition.connect(
      [&]() { this->particle_handler.register_load_callback_function(false); });

    // 这将初始化粒子所处的背景三角，以及粒子的属性数量。

    particle_handler.initialize(background_triangulation, mapping, 1 + dim);

    // 我们创建了一个粒子三角图，这个三角图只用来生成将用于插入粒子的点。这个三角形是一个偏离模拟域中心的超壳。这将被用来生成一个充满粒子的圆盘，这将使我们能够很容易地监测由于涡流而产生的运动。

    Point<dim> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (dim == 3)
      center[2] = 0.5;

    const double outer_radius = 0.15;
    const double inner_radius = 0.01;

    parallel::distributed::Triangulation<dim> particle_triangulation(
      MPI_COMM_WORLD);

    GridGenerator::hyper_shell(
      particle_triangulation, center, inner_radius, outer_radius, 6);
    particle_triangulation.refine_global(par.particle_insertion_refinement);

    // 我们为粒子发生器生成必要的边界盒。这些边界框是快速识别插入的粒子位于哪个进程的子域中，以及哪个单元拥有它的必要条件。

    const auto my_bounding_box = GridTools::compute_mesh_predicate_bounding_box(
      background_triangulation, IteratorFilters::LocallyOwnedCell());
    const auto global_bounding_boxes =
      Utilities::MPI::all_gather(MPI_COMM_WORLD, my_bounding_box);

    // 我们生成一个空的属性向量。一旦粒子生成，我们将把这些属性赋予它们。

    std::vector<std::vector<double>> properties(
      particle_triangulation.n_locally_owned_active_cells(),
      std::vector<double>(dim + 1, 0.));

    // 我们在单点正交的位置生成粒子。因此，在每个单元的中心点将生成一个粒子。

    Particles::Generators::quadrature_points(particle_triangulation,
                                             QMidpoint<dim>(),
                                             global_bounding_boxes,
                                             particle_handler,
                                             mapping,
                                             properties);

    pcout << "Number of particles inserted: "
          << particle_handler.n_global_particles() << std::endl;
  }

  //  @sect4{Background DOFs and interpolation}

  // 这个函数设置了用于速度插值的背景自由度，并分配了存储整个速度场解决方案的场向量。

  template <int dim>
  void
  ParticleTracking<dim>::setup_background_dofs()
  {
    fluid_dh.distribute_dofs(fluid_fe);
    const IndexSet locally_owned_dofs = fluid_dh.locally_owned_dofs();
    IndexSet       locally_relevant_dofs;
    DoFTools::extract_locally_relevant_dofs(fluid_dh, locally_relevant_dofs);

    velocity_field.reinit(locally_owned_dofs,
                          locally_relevant_dofs,
                          mpi_communicator);
  }

  // 这个函数负责将涡流速度场插值到场矢量上。这可以通过使用
  // VectorTools::interpolate() 函数相当容易地实现。

  template <int dim>
  void
  ParticleTracking<dim>::interpolate_function_to_field()
  {
    velocity_field.zero_out_ghost_values();
    VectorTools::interpolate(mapping, fluid_dh, velocity, velocity_field);
    velocity_field.update_ghost_values();
  }

  //  @sect4{Time integration of the trajectories}

  // 我们使用分析定义的速度场来整合粒子的轨迹。这展示了粒子的一个相对微不足道的用法。

  template <int dim>
  void
  ParticleTracking<dim>::euler_step_analytical(const double dt)
  {
    const unsigned int this_mpi_rank =
      Utilities::MPI::this_mpi_process(mpi_communicator);
    Vector<double> particle_velocity(dim);

    // 使用粒子迭代器在域中的所有粒子上进行循环操作

    for (auto &particle : particle_handler)
      {
        // 我们使用粒子的当前位置来计算它们的速度。

        Point<dim> particle_location = particle.get_location();
        velocity.vector_value(particle_location, particle_velocity);

        // 这就更新了粒子的位置，并将旧的位置设定为等于粒子的新位置。

        for (int d = 0; d < dim; ++d)
          particle_location[d] += particle_velocity[d] * dt;

        particle.set_location(particle_location);

        // 我们在粒子属性中存储处理器ID（标量）和粒子速度（矢量）。在这个例子中，这样做纯粹是为了可视化的目的。

        ArrayView<double> properties = particle.get_properties();
        for (int d = 0; d < dim; ++d)
          properties[d] = particle_velocity[d];
        properties[dim] = this_mpi_rank;
      }
  }

  // 与前面的函数不同，在这个函数中，我们通过将自由度处的速度场值插值到粒子的位置来积分粒子的轨迹。

  template <int dim>
  void
  ParticleTracking<dim>::euler_step_interpolated(const double dt)
  {
    Vector<double> local_dof_values(fluid_fe.dofs_per_cell);

    // 我们在所有的本地粒子上循环。虽然这可以直接通过循环所有的单元格来实现，但这将迫使我们循环许多不包含粒子的单元格。相反，我们在所有的粒子上循环，但是，我们得到粒子所在的单元格的引用，然后在该单元格中循环所有的粒子。这使我们能够从
    // "velocity_field "向量中收集一次速度值，并将其用于该单元中的所有粒子。

    auto particle = particle_handler.begin();
    while (particle != particle_handler.end())
      {
        const auto cell =
          particle->get_surrounding_cell(background_triangulation);
        const auto dh_cell =
          typename DoFHandler<dim>::cell_iterator(*cell, &fluid_dh);

        dh_cell->get_dof_values(velocity_field, local_dof_values);

        // 接下来，通过评估粒子位置的有限元解来计算粒子位置的速度。这基本上是第19步中粒子平流功能的优化版本，但我们不是为每个单元创建正交对象和FEValues对象，而是用手进行评估，这在一定程度上更有效率，而且只对本教程有意义，因为粒子工作是整个程序的主要成本。

        const auto pic = particle_handler.particles_in_cell(cell);
        Assert(pic.begin() == particle, ExcInternalError());
        for (auto &p : pic)
          {
            const Point<dim> reference_location = p.get_reference_location();
            Tensor<1, dim>   particle_velocity;
            for (unsigned int j = 0; j < fluid_fe.dofs_per_cell; ++j)
              {
                const auto comp_j = fluid_fe.system_to_component_index(j);

                particle_velocity[comp_j.first] +=
                  fluid_fe.shape_value(j, reference_location) *
                  local_dof_values[j];
              }

            Point<dim> particle_location = particle->get_location();
            for (int d = 0; d < dim; ++d)
              particle_location[d] += particle_velocity[d] * dt;
            p.set_location(particle_location);

            // 同样，我们在粒子属性中存储了粒子速度和处理器ID，以便于可视化。

            ArrayView<double> properties = p.get_properties();
            for (int d = 0; d < dim; ++d)
              properties[d] = particle_velocity[d];

            properties[dim] =
              Utilities::MPI::this_mpi_process(mpi_communicator);

            ++particle;
          }
      }
  }

  //  @sect4{Data output}

  // 接下来的两个函数负责将粒子和背景网格用pvtu记录写入vtu中。这可以确保在并行启动仿真时，仿真结果可以被可视化。

  template <int dim>
  void
  ParticleTracking<dim>::output_particles(const unsigned int it)
  {
    Particles::DataOut<dim, dim> particle_output;

    std::vector<std::string> solution_names(dim, "velocity");
    solution_names.push_back("process_id");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(
      DataComponentInterpretation::component_is_scalar);

    particle_output.build_patches(particle_handler,
                                  solution_names,
                                  data_component_interpretation);
    const std::string output_folder(par.output_directory);
    const std::string file_name(interpolated_velocity ?
                                  "interpolated-particles" :
                                  "analytical-particles");

    pcout << "Writing particle output file: " << file_name << "-" << it
          << std::endl;

    particle_output.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  template <int dim>
  void
  ParticleTracking<dim>::output_background(const unsigned int it)
  {
    std::vector<std::string> solution_names(dim, "velocity");
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation(
        dim, DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim> data_out;

    // 将解决方案的数据附加到data_out对象上

    data_out.attach_dof_handler(fluid_dh);
    data_out.add_data_vector(velocity_field,
                             solution_names,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);
    Vector<float> subdomain(background_triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = background_triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.build_patches(mapping);

    const std::string output_folder(par.output_directory);
    const std::string file_name("background");

    pcout << "Writing background field file: " << file_name << "-" << it
          << std::endl;

    data_out.write_vtu_with_pvtu_record(
      output_folder, file_name, it, mpi_communicator, 6);
  }

  //  @sect4{Running the simulation}
  //  这个函数协调了整个模拟过程。它与其他时间相关的教程程序非常相似--以 step-21
  //  或 step-26 为例。注意，我们使用DiscreteTime类来监控时间、时间步长和 step-
  //  号。这个函数相对来说是比较简单的。

  template <int dim>
  void
  ParticleTracking<dim>::run()
  {
    DiscreteTime discrete_time(0, par.final_time, par.time_step);

    generate_particles();

    pcout << "Repartitioning triangulation after particle generation"
          << std::endl;
    background_triangulation.repartition();

    // 我们通过在分析法和插值法的情况下进行时间步长为0的显式欧拉迭代来设置粒子的初始属性。

    if (interpolated_velocity)
      {
        setup_background_dofs();
        interpolate_function_to_field();
        euler_step_interpolated(0.);
      }
    else
      euler_step_analytical(0.);

    output_particles(discrete_time.get_step_number());
    if (interpolated_velocity)
      output_background(discrete_time.get_step_number());

    // 粒子通过循环的方式随时间推移而平移。

    while (!discrete_time.is_at_end())
      {
        discrete_time.advance_time();
        velocity.set_time(discrete_time.get_previous_time());

        if ((discrete_time.get_step_number() % par.repartition_frequency) == 0)
          {
            background_triangulation.repartition();
            if (interpolated_velocity)
              setup_background_dofs();
          }

        if (interpolated_velocity)
          {
            interpolate_function_to_field();
            euler_step_interpolated(discrete_time.get_previous_step_size());
          }
        else
          euler_step_analytical(discrete_time.get_previous_step_size());

        // 在粒子被移动之后，有必要确定它们现在所在的单元。这可以通过调用
        // <code>sort_particles_into_subdomains_and_cells</code> 来实现。
        particle_handler.sort_particles_into_subdomains_and_cells();

        if ((discrete_time.get_step_number() % par.output_frequency) == 0)
          {
            output_particles(discrete_time.get_step_number());
            if (interpolated_velocity)
              output_background(discrete_time.get_step_number());
          }
      }
  }

} // namespace Step68

//  @sect3{The main() function}

// 代码的其余部分，即`main()`函数，是标准的。我们注意到，我们用分析速度和插值速度运行粒子跟踪，并产生两种结果

int
main(int argc, char *argv[])
{
  using namespace Step68;
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

      ParticleTrackingParameters par;
      ParameterAcceptor::initialize(prm_file);
      {
        Step68::ParticleTracking<2> particle_tracking(par, false);
        particle_tracking.run();
      }
      {
        Step68::ParticleTracking<2> particle_tracking(par, true);
        particle_tracking.run();
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
