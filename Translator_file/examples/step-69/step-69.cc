

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2019 - 2021 by the deal.II authors 
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
 * Authors: Matthias Maier, Texas A&M University; 
 *          Ignacio Tomas, Texas A&M University, Sandia National Laboratories 
 * 
 * Sandia National Laboratories is a multimission laboratory managed and 
 * operated by National Technology & Engineering Solutions of Sandia, LLC, a 
 * wholly owned subsidiary of Honeywell International Inc., for the U.S. 
 * Department of Energy's National Nuclear Security Administration under 
 * contract DE-NA0003525. This document describes objective technical results 
 * and analysis. Any subjective views or opinions that might be expressed in 
 * the paper do not necessarily represent the views of the U.S. Department of 
 * Energy or the United States Government. 
 */ 


// @sect3{Include files}  

// 包含文件的集合是相当标准的。最耐人寻味的部分是，我们将完全依靠deal.II数据结构进行MPI并行化，特别是通过 parallel::distributed::Triangulation 和 LinearAlgebra::distributed::Vector 包含的 <code>distributed/tria.h</code> 和 <code>lac/la_parallel_vector.h</code>  。我们将使用非分布式的  dealii::SparseMatrix  (  <code>lac/sparse_matrix.h</code>  ) 来存储  $\mathbf{c}_{ij}$  、  $\mathbf{n}_{ij}$  和  $d_{ij}$  矩阵的本地部分，而不是 Trilinos 或 PETSc 特定的矩阵类。

#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/parallel.h> 
#include <deal.II/base/parameter_acceptor.h> 
#include <deal.II/base/partitioner.h> 
#include <deal.II/base/quadrature.h> 
#include <deal.II/base/timer.h> 
#include <deal.II/base/work_stream.h> 

#include <deal.II/distributed/tria.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/la_parallel_vector.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/sparse_matrix.templates.h> 
#include <deal.II/lac/vector.h> 

#include <deal.II/meshworker/scratch_data.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

// 除了上述deal.II的具体内容外，我们还包括四个提升头文件。前两个是二进制文件，我们将用它来实现检查点和重启机制。

#include <boost/archive/binary_iarchive.hpp> 
#include <boost/archive/binary_oarchive.hpp> 

// 最后两个头文件是用来在整数间隔上创建自定义迭代器范围。

#include <deal.II/base/std_cxx20/iota_view.h> 
#include <boost/range/iterator_range.hpp> 

//用于  std::isnan,  。
// std::isinf,  
// std::ifstream,  
// std::async, 和 std::future  。
#include <cmath> 
#include <fstream> 
#include <future> 
// @sect3{Class template declarations}  

// 我们开始实际的实现，先声明所有的类及其数据结构和方法。与之前的例子步骤相比，我们使用了更精细的概念、数据结构和参数封装到各个类中。因此，一个单一的类通常围绕着一个单一的数据结构（如 <code>Discretization</code> 类中的Triangulation），或者一个单一的方法（如 <code>%TimeStepping</code> 类的 <code>make_one_step()</code> 函数）。我们通常声明参数变量和从头开始的数据对象为`private'，而使其他类使用的方法和数据结构为`public'。

//  @note  一个更简洁的方法是通过<a
//  href="https:en.wikipedia.org/wiki/Mutator_method">getter/setter functions</a>来保护对所有数据结构的访问。为了简洁起见，我们不采用这种方法。

// 我们还注意到，绝大多数的类都是从ParameterAcceptor派生的。这有利于将所有的全局参数归入一个（全局）ParameterHandler。关于从ParameterAcceptor继承作为全局订阅机制的更多解释可以在  step-60  中找到。

namespace Step69 
{ 
  using namespace dealii; 

// 我们首先定义一些 types::boundary_id 常量，用于整个教程步骤。这使得我们可以用一个助记符（如 <code>do_nothing</code>  ）而不是一个数值来指代边界类型。

  namespace Boundaries 
  { 
    constexpr types::boundary_id do_nothing = 0; 
    constexpr types::boundary_id free_slip  = 1; 
    constexpr types::boundary_id dirichlet  = 2; 
  } // namespace Boundaries 
// @sect4{The <code>Discretization</code> class}  

//  <code>Discretization</code> 类包含所有关于问题的网格（三角形）和离散化（映射、有限元、正交）的数据结构。如前所述，我们使用ParameterAcceptor类来自动填充特定问题的参数，如几何信息（ <code>length</code>  等）或来自参数文件的细化水平（ <code>refinement</code>  ）。这就要求我们把数据结构的初始化分成两个函数。我们在构造函数中初始化所有不依赖参数的东西，并将网格的创建推迟到 <code>setup()</code> 方法中，一旦所有参数通过 ParameterAcceptor::initialize(). 读入，就可以调用该方法。
  template <int dim> 
  class Discretization : public ParameterAcceptor 
  { 
  public: 
    Discretization(const MPI_Comm     mpi_communicator, 
                   TimerOutput &      computing_timer, 
                   const std::string &subsection = "Discretization"); 

    void setup(); 

    const MPI_Comm mpi_communicator; 

    parallel::distributed::Triangulation<dim> triangulation; 

    const MappingQ<dim>   mapping; 
    const FE_Q<dim>       finite_element; 
    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> face_quadrature; 

  private: 
    TimerOutput &computing_timer; 

    double length; 
    double height; 
    double disk_position; 
    double disk_diameter; 

    unsigned int refinement; 
  }; 
// @sect4{The <code>OfflineData</code> class}  

//  <code>OfflineData</code> 类包含了离散化中几乎所有不随时间演变的组件，特别是DoFHandler、SparsityPattern、边界图、块状质量矩阵、 $\mathbf{c}_{ij}$ 和 $\mathbf{n}_{ij}$ 矩阵。这里，术语<i>offline</i>指的是 <code>OfflineData</code> 的所有类成员都有明确定义的值，与当前时间步长无关。这意味着它们可以提前初始化（在<i>time step zero</i>），并且不意味着在任何后来的时间步长中被修改。例如，稀疏模式不应该随着时间的推进而改变（我们在空间上不做任何形式的适应性）。同样地，包络质量矩阵的条目也不应该随着时间的推进而被修改。

// 我们还计算并存储一个 <code>boundary_normal_map</code> ，它包含一个从边界自由度的 types::global_dof_index 类型的全局索引到一个由法向量、边界ID和与自由度相关的位置组成的元组的映射。我们必须在这个类中计算和存储这些几何信息，因为我们在后面的代数循环中无法获得几何（或基于单元）的信息。

// 尽管这个类目前没有任何可以从参数文件中读入的参数，但我们还是从ParameterAcceptor派生出来，并遵循与Discretization类相同的习惯，提供一个 <code>setup()</code> (and <code>assemble()</code>  )方法。

  template <int dim> 
  class OfflineData : public ParameterAcceptor 
  { 
  public: 
    using BoundaryNormalMap = 
      std::map<types::global_dof_index, 
               std::tuple<Tensor<1, dim>, types::boundary_id, Point<dim>>>; 

    OfflineData(const MPI_Comm             mpi_communicator, 
                TimerOutput &              computing_timer, 
                const Discretization<dim> &discretization, 
                const std::string &        subsection = "OfflineData"); 

    void setup(); 
    void assemble(); 

    DoFHandler<dim> dof_handler; 

    std::shared_ptr<const Utilities::MPI::Partitioner> partitioner; 

    unsigned int n_locally_owned; 
    unsigned int n_locally_relevant; 

    SparsityPattern sparsity_pattern; 

    BoundaryNormalMap boundary_normal_map; 

    SparseMatrix<double>                  lumped_mass_matrix; 
    std::array<SparseMatrix<double>, dim> cij_matrix; 
    std::array<SparseMatrix<double>, dim> nij_matrix; 
    SparseMatrix<double>                  norm_matrix; 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const Discretization<dim>> discretization; 
  }; 
// @sect4{The <code>ProblemDescription</code> class}  

// 该类的成员函数是欧拉方程特有的实用函数和数据结构。

// - 类型别名 <code>state_type</code> 用于状态 $\mathbf{U}_i^n$  。

// - 类型别名  <code>flux_type</code>  用来表示通量  $\mathbb{f}(\mathbf{U}_j^n)$  。

// -  <code>momentum</code> 函数从状态向量 $[\rho,\textbf{m},E]$ 中提取 $\textbf{m}$ 并存储在一个 <code>Tensor<1, dim></code> 中。

// -  <code>internal_energy</code> 函数从给定的状态向量 $[\rho,\textbf{m},E]$ 中计算 $E - \frac{|\textbf{m}|^2}{2\rho}$  。

// 类成员  <code>component_names</code>  ,  <code>pressure</code>, and <code>speed_of_sound</code>  的目的从它们的名字中就可以看出。我们还提供了一个函数  <code>compute_lambda_max()</code>  ，用于计算上面提到的波速估计，  $\lambda_{max}(\mathbf{U},\mathbf{V},\mathbf{n})$  ，用于计算  $d_{ij}$  矩阵。

//  @note   <code>DEAL_II_ALWAYS_INLINE</code> 宏扩展为一个（编译器特定的）pragma，确保这个类中定义的相应函数总是内联的，也就是说，每次调用该函数时，函数体都被放在原位，而不会产生调用（和代码转接）。这比 <code>inline</code> 关键字要强，后者或多或少是对编译器的一个（温和的）建议，即程序员认为内联函数是有益的。  <code>DEAL_II_ALWAYS_INLINE</code> 只应在很少的情况下谨慎使用，比如在这种情况下，我们实际上知道（由于基准测试）内联有关的函数可以提高性能。

// 最后，我们注意到这是本教程步骤中唯一一个与特定的 "物理学 "或 "双曲守恒定律"（本例中为欧拉方程）相关的类。所有其他的类主要是 "离散化 "类，与所求解的特定物理学无关。

  template <int dim> 
  class ProblemDescription 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 2 + dim; 

    using state_type = Tensor<1, problem_dimension>; 
    using flux_type  = Tensor<1, problem_dimension, Tensor<1, dim>>; 

    const static std::array<std::string, problem_dimension> component_names; 

    static constexpr double gamma = 7. / 5.; 

    static DEAL_II_ALWAYS_INLINE inline Tensor<1, dim> 
    momentum(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    internal_energy(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double pressure(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    speed_of_sound(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline flux_type flux(const state_type &U); 

    static DEAL_II_ALWAYS_INLINE inline double 
    compute_lambda_max(const state_type &    U_i, 
                       const state_type &    U_j, 
                       const Tensor<1, dim> &n_ij); 
  }; 
// @sect4{The <code>InitialValues</code> class}  

//  <code>InitialValues</code> 类的唯一公共数据属性是一个 std::function 。
// <code>initial_state</code> ，用于计算给定的点和时间的初始状态。这个函数用于填充初始流场，以及在每个时间步长中明确设置迪里切特边界条件（在流入边界）。

// 在这个例子的步骤中，我们简单地实现了一个均匀的流场，其方向和一维原始状态（密度、速度、压力）从参数文件中读取。

// 最好是一次性初始化这个类：初始化/设置参数并定义依赖于这些默认参数的类成员。然而，由于我们不知道参数的实际值，这在一般情况下是毫无意义和不安全的（我们希望有机制来检查输入参数的一致性）。我们没有定义另一个 <code>setup()</code> 方法在调用 ParameterAcceptor::initialize() 后被调用（手动），而是为类成员 <code>parse_parameters_call_back()</code> 提供了一个 "实现"，当调用 ParameterAcceptor::initialize() 时，每个继承自ParameterAceptor的类都会自动调用。

  template <int dim> 
  class InitialValues : public ParameterAcceptor 
  { 
  public: 
    using state_type = typename ProblemDescription<dim>::state_type; 

    InitialValues(const std::string &subsection = "InitialValues"); 

    std::function<state_type(const Point<dim> &point, double t)> initial_state; 

  private: 

// 我们声明一个私有的回调函数，它将与 ParameterAcceptor::parse_parameters_call_back 信号相连接。

    void parse_parameters_callback(); 

    Tensor<1, dim> initial_direction; 
    Tensor<1, 3>   initial_1d_state; 
  }; 
// @sect4{The <code>%TimeStepping</code> class}  

// 有了 <code>OfflineData</code> and <code>ProblemDescription</code> 类在手，我们现在可以实现上面讨论中介绍的显式时间步进方案。 <code>%TimeStepping</code> 类的主要方法是<code>make_one_step(vector_type &U, double t)</code>，它接受对状态向量 <code>U</code> and a time point <code>t</code> 的引用（作为输入参数）计算更新的解决方案，将其存储在向量 <code>temp</code>, swaps its contents with the vector <code>U</code> 中，并返回选择的 step- 大小 $\tau$  。

// 另一个重要的方法是  <code>prepare()</code>  ，主要是为临时向量  <code>temp</code> and the matrix <code>dij_matrix</code>  分别设置适当的分区和稀疏模式。

  template <int dim> 
  class TimeStepping : public ParameterAcceptor 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

    using state_type = typename ProblemDescription<dim>::state_type; 
    using flux_type  = typename ProblemDescription<dim>::flux_type; 

    using vector_type = 
      std::array<LinearAlgebra::distributed::Vector<double>, problem_dimension>; 

    TimeStepping(const MPI_Comm            mpi_communicator, 
                 TimerOutput &             computing_timer, 
                 const OfflineData<dim> &  offline_data, 
                 const InitialValues<dim> &initial_values, 
                 const std::string &       subsection = "TimeStepping"); 

    void prepare(); 

    double make_one_step(vector_type &U, double t); 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const OfflineData<dim>>   offline_data; 
    SmartPointer<const InitialValues<dim>> initial_values; 

    SparseMatrix<double> dij_matrix; 

    vector_type temporary_vector; 

    double cfl_update; 
  }; 
// @sect4{The <code>SchlierenPostprocessor</code> class}  

// 在其核心中，Schlieren类实现了类成员  <code>compute_schlieren()</code>  。这个类成员的主要目的是计算一个辅助的有限元场 <code>schlieren</code>  ，它在每个节点上由\f[ \text{schlieren}[i] = e^{\beta \frac{ |\nabla r_i| - \min_j |\nabla r_j| }{\max_j |\nabla r_j| - \min_j |\nabla r_j| } }, \f]定义，其中 $r$ 原则上可以是任何标量。但在实践中，密度是一个自然的候选量，即 $r \dealcoloneq \rho$  。<a href="https:en.wikipedia.org/wiki/Schlieren">Schlieren</a>后处理是一种标准的方法，用于增强可视化的对比度，其灵感来自实际的实验X射线和可视化的阴影技术。(参见  step-67  另一个例子，我们创建了一个Schlieren图。)

  template <int dim> 
  class SchlierenPostprocessor : public ParameterAcceptor 
  { 
  public: 
    static constexpr unsigned int problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

    using state_type = typename ProblemDescription<dim>::state_type; 

    using vector_type = 
      std::array<LinearAlgebra::distributed::Vector<double>, problem_dimension>; 

    SchlierenPostprocessor( 
      const MPI_Comm          mpi_communicator, 
      TimerOutput &           computing_timer, 
      const OfflineData<dim> &offline_data, 
      const std::string &     subsection = "SchlierenPostprocessor"); 

    void prepare(); 

    void compute_schlieren(const vector_type &U); 

    LinearAlgebra::distributed::Vector<double> schlieren; 

  private: 
    const MPI_Comm mpi_communicator; 
    TimerOutput &  computing_timer; 

    SmartPointer<const OfflineData<dim>> offline_data; 

    Vector<double> r; 

    unsigned int schlieren_index; 
    double       schlieren_beta; 
  }; 
// @sect4{The <code>MainLoop</code> class}  

// 现在，剩下的就是把 <code>%TimeStepping</code>, <code>InitialValues</code>  , 和  <code>SchlierenPostprocessor</code>  类中实现的方法连在一起。我们在一个单独的类 <code>MainLoop</code> 中做到这一点，该类包含每个类的一个对象，并在ParameterAcceptor类的帮助下再次读入一些参数。

  template <int dim> 
  class MainLoop : public ParameterAcceptor 
  { 
  public: 
    using vector_type = typename TimeStepping<dim>::vector_type; 

    MainLoop(const MPI_Comm mpi_communnicator); 

    void run(); 

  private: 
    vector_type interpolate_initial_values(const double t = 0); 

    void output(const vector_type &U, 
                const std::string &name, 
                double             t, 
                unsigned int       cycle, 
                bool               checkpoint = false); 

    const MPI_Comm     mpi_communicator; 
    std::ostringstream timer_output; 
    TimerOutput        computing_timer; 

    ConditionalOStream pcout; 

    std::string base_name; 
    double      t_final; 
    double      output_granularity; 

    bool asynchronous_writeback; 

    bool resume; 

    Discretization<dim>         discretization; 
    OfflineData<dim>            offline_data; 
    InitialValues<dim>          initial_values; 
    TimeStepping<dim>           time_stepping; 
    SchlierenPostprocessor<dim> schlieren_postprocessor; 

    vector_type output_vector; 

    std::future<void> background_thread_state; 
  }; 
// @sect3{Implementation}  
// @sect4{Grid generation, setup of data structures}  

// 手头的第一个主要任务是典型的网格生成、数据结构的设置和装配这三者。在这个例子的步骤中，一个值得注意的创新是使用ParameterAcceptor类，我们用它来填充参数值：我们首先初始化ParameterAcceptor类，用一个字符串 <code>subsection</code> 表示参数文件中的正确分节，调用它的构造器。然后，在构造函数中，每个参数值都被初始化为一个合理的默认值，并通过调用 ParameterAcceptor::add_parameter(). 向ParameterAcceptor类注册。
  template <int dim> 
  Discretization<dim>::Discretization(const MPI_Comm     mpi_communicator, 
                                      TimerOutput &      computing_timer, 
                                      const std::string &subsection) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , triangulation(mpi_communicator) 
    , mapping(1) 
    , finite_element(1) 
    , quadrature(3) 
    , face_quadrature(3) 
    , computing_timer(computing_timer) 
  { 
    length = 4.; 
    add_parameter("length", length, "Length of computational domain"); 

    height = 2.; 
    add_parameter("height", height, "Height of computational domain"); 

    disk_position = 0.6; 
    add_parameter("object position", 
                  disk_position, 
                  "x position of immersed disk center point"); 

    disk_diameter = 0.5; 
    add_parameter("object diameter", 
                  disk_diameter, 
                  "Diameter of immersed disk"); 

    refinement = 5; 
    add_parameter("refinement", 
                  refinement, 
                  "Number of refinement steps of the geometry"); 
  } 

// 注意在前面的构造函数中，我们只把MPI通信器传给了 <code>triangulation</code> ，但我们仍然没有初始化底层几何体/网格。如前所述，我们必须将这项任务推迟到 <code>setup()</code> 函数，在 ParameterAcceptor::initialize() 函数用从参数文件中读取的最终值填充所有参数变量后，再调用该函数。

//  <code>setup()</code> 函数是最后一个必须实现的类成员。它创建了实际的三角结构，这是一个基准配置，由一个带有盘状障碍物的通道组成，见  @cite GuermondEtAl2018  。我们通过修改 GridGenerator::hyper_cube_with_cylindrical_hole(). 生成的网格来构建几何体。我们参考 step-49 、 step-53 和 step-54 来了解如何创建高级网格。我们首先创建4个临时的（非分布式的）粗略三角形，用 GridGenerator::merge_triangulation() 函数将其缝合起来。我们在 $(0,0)$ 处将圆盘居中，直径为 <code>disk_diameter</code>  。通道的左下角有坐标（  <code>-disk_position</code>, <code>-height/2</code>  ），右上角有（  <code>length-disk_position</code>  ,  <code>height/2</code>  ）。

  template <int dim> 
  void Discretization<dim>::setup() 
  { 
    TimerOutput::Scope scope(computing_timer, "discretization - setup"); 

    triangulation.clear(); 

    Triangulation<dim> tria1, tria2, tria3, tria4, tria5, tria6; 

    GridGenerator::hyper_cube_with_cylindrical_hole( 
      tria1, disk_diameter / 2., disk_diameter, 0.5, 1, false); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria2, 
      {2, 1}, 
      Point<2>(-disk_diameter, disk_diameter), 
      Point<2>(disk_diameter, height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria3, 
      {2, 1}, 
      Point<2>(-disk_diameter, -disk_diameter), 
      Point<2>(disk_diameter, -height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria4, 
      {6, 2}, 
      Point<2>(disk_diameter, -disk_diameter), 
      Point<2>(length - disk_position, disk_diameter)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria5, 
      {6, 1}, 
      Point<2>(disk_diameter, disk_diameter), 
      Point<2>(length - disk_position, height / 2.)); 

    GridGenerator::subdivided_hyper_rectangle( 
      tria6, 
      {6, 1}, 
      Point<2>(disk_diameter, -height / 2.), 
      Point<2>(length - disk_position, -disk_diameter)); 

    GridGenerator::merge_triangulations( 
      {&tria1, &tria2, &tria3, &tria4, &tria5, &tria6}, 
      triangulation, 
      1.e-12, 
      true); 

    triangulation.set_manifold(0, PolarManifold<2>(Point<2>())); 

// 我们必须修复目前位于 $x=-$ 的左边缘。
// <code>disk_diameter</code> ，必须移到 $x=-$ 。
// <code>disk_position</code>  . 作为最后一步，边界必须被着色，右边是 <code>Boundaries::do_nothing</code> ， <code>dirichlet</code> on the left and <code>free_slip</code> 是上、下外边界和障碍物。

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto v : cell->vertex_indices()) 
          { 
            if (cell->vertex(v)[0] <= -disk_diameter + 1.e-6) 
              cell->vertex(v)[0] = -disk_position; 
          } 
      } 

    for (const auto &cell : triangulation.active_cell_iterators()) 
      { 
        for (const auto f : cell->face_indices()) 
          { 
            const auto face = cell->face(f); 

            if (face->at_boundary()) 
              { 
                const auto center = face->center(); 

                if (center[0] > length - disk_position - 1.e-6) 
                  face->set_boundary_id(Boundaries::do_nothing); 
                else if (center[0] < -disk_position + 1.e-6) 
                  face->set_boundary_id(Boundaries::dirichlet); 
                else 
                  face->set_boundary_id(Boundaries::free_slip); 
              } 
          } 
      } 

    triangulation.refine_global(refinement); 
  } 
// @sect4{Assembly of offline matrices}  

// 在 <code>OfflineData</code> 的构造函数中，除了在初始化列表中初始化相应的类成员外，没有做太多的工作。

  template <int dim> 
  OfflineData<dim>::OfflineData(const MPI_Comm             mpi_communicator, 
                                TimerOutput &              computing_timer, 
                                const Discretization<dim> &discretization, 
                                const std::string &        subsection) 
    : ParameterAcceptor(subsection) 
    , dof_handler(discretization.triangulation) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , discretization(&discretization) 
  {} 

// 现在我们可以初始化DoFHandler，为本地拥有的和本地相关的DOF提取IndexSet对象，并初始化一个 Utilities::MPI::Partitioner 对象，这是分布式向量需要的。

  template <int dim> 
  void OfflineData<dim>::setup() 
  { 
    IndexSet locally_owned; 
    IndexSet locally_relevant; 

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "offline_data - distribute dofs"); 

      dof_handler.distribute_dofs(discretization->finite_element); 

      locally_owned   = dof_handler.locally_owned_dofs(); 
      n_locally_owned = locally_owned.n_elements(); 

      DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant); 
      n_locally_relevant = locally_relevant.n_elements(); 

      partitioner = 
        std::make_shared<Utilities::MPI::Partitioner>(locally_owned, 
                                                      locally_relevant, 
                                                      mpi_communicator); 
    } 
// @sect4{Translation to local index ranges}  

// 我们现在可以为我们的矩阵创建稀疏模式了。有相当多的特殊性需要详细解释。我们避免使用分布式矩阵类（例如由Trilinos或PETSc提供的），而是依靠deal.II自己的SparseMatrix对象来存储所有矩阵的局部部分。这一设计决定的动机是：(a)我们实际上从未进行过矩阵-向量乘法，(b)我们总是可以在一个给定的MPI等级上专门组装矩阵的局部部分。相反，我们将计算非线性更新，同时迭代连通性模版的（局部）部分；这是deal.II自己的SparsityPattern专门为之优化的任务。

// 不过，这种设计考虑有一个注意事项。让deal.II SparseMatrix类变得快速的是SparsityPattern中使用的<a
//  href="https:en.wikipedia.org/wiki/Sparse_matrix">compressed row
//  storage (CSR)</a>（见 @ref Sparsity  ）。不幸的是，这与全局分布式索引范围不相称，因为具有CSR的稀疏模式不能在索引范围内包含 "洞"。deal.II提供的分布式矩阵通过将全局索引范围转化为连续的局部索引范围来避免这一点。但这正是我们在迭代模版时想要避免的索引操作类型，因为它产生了可衡量的开销。

//  Utilities::MPI::Partitioner 类已经实现了从全局索引范围到连续的局部（每个MPI等级）索引范围的转换：我们不需要重新发明轮子。我们只需要使用这种转换能力（一次，而且只有一次），以便为连续的索引范围创建一个 "本地 "稀疏模式  $[0,$  。
// <code>n_locally_relevant</code>  
// $)$  . 这种能力可以通过 Utilities::MPI::Partitioner::global_to_local() 函数来调用。一旦使用本地索引创建了稀疏模式，剩下要做的就是确保（在实现我们的scatter和gather辅助函数时）我们总是通过调用 LinearAlgebra::distributed::Vector::local_element(). 来访问分布式向量的元素，这样我们就完全避免了索引转换，并完全使用本地索引进行操作。

    { 
      TimerOutput::Scope scope( 
        computing_timer, 
        "offline_data - create sparsity pattern and set up matrices"); 

// 我们必须手工创建 "本地 "稀疏模式。因此，我们在所有本地拥有的和重影的单元上循环（见  @ref  GlossArtificialCell），并提取与单元DOF相关的（全局）  <code>dof_indices</code>  ，并使用  <code>partitioner->global_to_local(index)</code>  重新编号。

// 在本地拥有的DOF的情况下，这种重新编号包括应用一个移位（即我们减去一个偏移量），这样，现在它们将成为整数区间 $[0,$ 中的一个数字。
// <code>n_locally_owned</code>   $)$  .
//然而，在重影道次的情况下（即不是本地拥有的），情况就完全不同了，因为与重影道次相关的全局指数（一般来说）不会是一个连续的整数集。

      DynamicSparsityPattern dsp(n_locally_relevant, n_locally_relevant); 

      const auto dofs_per_cell = 
        discretization->finite_element.n_dofs_per_cell(); 
      std::vector<types::global_dof_index> dof_indices(dofs_per_cell); 

      for (const auto &cell : dof_handler.active_cell_iterators()) 
        { 
          if (cell->is_artificial()) 
            continue; 

          /* We transform the set of global dof indices on the cell to the
           * corresponding "local" index range on the MPI process: */
          cell->get_dof_indices(dof_indices); 
          std::transform(dof_indices.begin(), 
                         dof_indices.end(), 
                         dof_indices.begin(), 
                         [&](types::global_dof_index index) { 
                           return partitioner->global_to_local(index); 
                         }); 

/* 为每个dof简单地添加一个与所有其他 "本地 "的联接。 */

           /* dofs on the cell: */ 


          for (const auto dof : dof_indices) 
            dsp.add_entries(dof, dof_indices.begin(), dof_indices.end()); 
        } 

      sparsity_pattern.copy_from(dsp); 

      lumped_mass_matrix.reinit(sparsity_pattern); 
      norm_matrix.reinit(sparsity_pattern); 
      for (auto &matrix : cij_matrix) 
        matrix.reinit(sparsity_pattern); 
      for (auto &matrix : nij_matrix) 
        matrix.reinit(sparsity_pattern); 
    } 
  } 

// DoFHandler和SparseMatrix对象的设置到此结束。接下来，我们要组装各种矩阵。我们在一个匿名命名空间中定义了一些辅助函数和数据结构。

  namespace 
  { 
// <code>CopyData</code> 类，将用于使用WorkStream组装离线数据矩阵。它作为一个容器：它只是一个结构，WorkStream在其中存储本地单元的贡献。请注意，它还包含一个类成员 <code>local_boundary_normal_map</code> ，用于存储计算边界法线所需的局部贡献。

    template <int dim> 
    struct CopyData 
    { 
      bool                                         is_artificial; 
      std::vector<types::global_dof_index>         local_dof_indices; 
      typename OfflineData<dim>::BoundaryNormalMap local_boundary_normal_map; 
      FullMatrix<double>                           cell_lumped_mass_matrix; 
      std::array<FullMatrix<double>, dim>          cell_cij_matrix; 
    }; 

// 接下来我们介绍一些辅助函数，它们都是关于读写矩阵和向量条目的。它们的主要动机是提供稍微有效的代码和<a href="https:en.wikipedia.org/wiki/Syntactic_sugar"> syntactic sugar</a>的代码，否则就有些乏味了。

// 我们介绍的第一个函数  <code>get_entry()</code>  ，将用于读取SparsityPattern迭代器  <code>it</code> of <code>matrix</code>  指向的条目所存储的值。该函数绕过了SparseMatrix接口中的一个小缺陷。SparsityPattern关注的是以CRS格式存储的稀疏矩阵的所有索引操作。因此，迭代器已经知道存储在SparseMatrix对象中的低级向量中相应矩阵条目的全局索引。由于SparseMatrix中缺乏直接用SparsityPattern迭代器访问该元素的接口，不幸的是我们必须创建一个临时的SparseMatrix迭代器。我们只需将其隐藏在 <code>get_entry()</code> 函数中。

    template <typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline SparseMatrix<double>::value_type 
    get_entry(const SparseMatrix<double> &matrix, const IteratorType &it) 
    { 
      const SparseMatrix<double>::const_iterator matrix_iterator( 
        &matrix, it->global_index()); 
      return matrix_iterator->value(); 
    } 

//  <code>set_entry()</code> 帮助器是 <code>get_value()</code> 的逆运算：给定一个迭代器和一个值，它在矩阵中设置迭代器所指向的条目。

    template <typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline void 
    set_entry(SparseMatrix<double> &           matrix, 
              const IteratorType &             it, 
              SparseMatrix<double>::value_type value) 
    { 
      SparseMatrix<double>::iterator matrix_iterator(&matrix, 
                                                     it->global_index()); 
      matrix_iterator->value() = value; 
    } 
// <code>gather_get_entry()</code>  : 我们注意到 $\mathbf{c}_{ij} \in \mathbb{R}^d$  。如果 $d=2$ ，那么 $\mathbf{c}_{ij} = [\mathbf{c}_{ij}^1,\mathbf{c}_{ij}^2]^\top$  。这基本上意味着我们需要每个空间维度的一个矩阵来存储 $\mathbf{c}_{ij}$ 向量。对于矩阵 $\mathbf{n}_{ij}$ 也有类似的观察。 <code>gather_get_entry()</code> 的目的是检索这些条目并将其存储到 <code>Tensor<1, dim></code> 中，以方便我们使用。

    template <std::size_t k, typename IteratorType> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather_get_entry(const std::array<SparseMatrix<double>, k> &c_ij, 
                     const IteratorType                         it) 
    { 
      Tensor<1, k> result; 
      for (unsigned int j = 0; j < k; ++j) 
        result[j] = get_entry(c_ij[j], it); 
      return result; 
    } 
// <code>gather()</code> （第一个接口）：这个第一个函数签名，有三个输入参数，将被用来检索矩阵的各个组成部分 <code>(i,l)</code> 。 <code>gather_get_entry()</code> 和 <code>gather()</code> 的功能非常相同，但它们的背景不同：函数 <code>gather()</code> 不依赖迭代器（实际上知道指向的值），而是依赖条目的索引 <code>(i,l)</code> ，以便检索其实际值。我们应该期望  <code>gather()</code>  比  <code>gather_get_entry()</code>  稍微昂贵一些。 <code>gather()</code> 的使用将限于计算代数粘度 $d_{ij}$ 的任务，在特殊情况下，当 $i$ 和 $j$ 都位于边界时。

//  @note  读者应该知道，访问一个矩阵的任意 <code>(i,l)</code> 条目（例如Trilinos或PETSc矩阵）一般来说是昂贵得不可接受的。在这里，我们可能要注意复杂度：我们希望这个操作有恒定的复杂度，这就是目前使用deal.II矩阵的实现的情况。

    template <std::size_t k> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather(const std::array<SparseMatrix<double>, k> &n_ij, 
           const unsigned int                         i, 
           const unsigned int                         j) 
    { 
      Tensor<1, k> result; 
      for (unsigned int l = 0; l < k; ++l) 
        result[l] = n_ij[l](i, j); 
      return result; 
    } 
// <code>gather()</code> （第二个接口）：这个有两个输入参数的第二个函数签名将被用来收集节点 <code>i</code> 的状态，并作为 <code>Tensor<1,problem_dimension></code> 返回，以方便我们使用。

    template <std::size_t k> 
    DEAL_II_ALWAYS_INLINE inline Tensor<1, k> 
    gather(const std::array<LinearAlgebra::distributed::Vector<double>, k> &U, 
           const unsigned int                                               i) 
    { 
      Tensor<1, k> result; 
      for (unsigned int j = 0; j < k; ++j) 
        result[j] = U[j].local_element(i); 
      return result; 
    } 
// <code>scatter()</code>  ：这个函数有三个输入参数，第一个是指一个 "全局对象"（比如一个本地拥有的或本地相关的矢量），第二个参数可以是一个 <code>Tensor<1,problem_dimension></code>  ，最后一个参数代表全局对象的索引。这个函数主要用于将更新的节点值（存储为 <code>Tensor<1,problem_dimension></code> ）写入全局对象中。

    template <std::size_t k, int k2> 
    DEAL_II_ALWAYS_INLINE inline void 
    scatter(std::array<LinearAlgebra::distributed::Vector<double>, k> &U, 
            const Tensor<1, k2> &                                      tensor, 
            const unsigned int                                         i) 
    { 
      static_assert(k == k2, 
                    "The dimensions of the input arguments must agree"); 
      for (unsigned int j = 0; j < k; ++j) 
        U[j].local_element(i) = tensor[j]; 
    } 
  } // namespace 

// 我们现在可以将储存在 <code>OfflineData</code> 中的所有矩阵集合起来：质量分录 $m_i$ ，矢量值矩阵 $\mathbf{c}_{ij}$ 和 $\mathbf{n}_{ij} = \frac{\mathbf{c}_{ij}}{|\mathbf{c}_{ij}|}$ ，以及边界法线 $\boldsymbol{\nu}_i$  。

// 为了利用线程并行化，我们使用了 @ref threads "多处理器的并行计算 "中详述的WorkStream方法来访问共享内存。按照惯例，这需要定义 

// - 抓取数据（即进行计算所需的输入信息）：在这种情况下，它是  <code>scratch_data</code>  。

// - 工作者：在我们的例子中，这是一个 <code>local_assemble_system()</code> 函数，它实际上是从抓取数据中计算出本地（即当前单元）贡献。

// - 拷贝数据：一个包含所有本地装配贡献的结构，在这里是  <code>CopyData<dim>()</code>  。

// - 一个拷贝数据的程序：在这种情况下，它是 <code>copy_local_to_global()</code> ，负责将这些局部贡献实际复制到全局对象（矩阵和/或矢量）中。

// 下面的大部分行是用来定义工作者  <code>local_assemble_system()</code>  和复制数据例程  <code>copy_local_to_global()</code>  的。关于WorkStream框架没有太多可说的，因为绝大多数的想法在  step-9  、  step-13  和  step-32  等文件中都有合理的记载。

// 最后，假设 $\mathbf{x}_i$ 是边界上的一个支持点，（节点）法线定义为。

// 
// @f{align*}
//  \widehat{\boldsymbol{\nu}}_i \dealcoloneq
//   \frac{\int_{\partial\Omega} \phi_i \widehat{\boldsymbol{\nu}} \,
//   \, \mathrm{d}\mathbf{s}}{\big|\int_{\partial\Omega} \phi_i
//   \widehat{\boldsymbol{\nu}} \, \mathrm{d}\mathbf{s}\big|}
//  @f}

// 我们将首先计算这个表达式的分子，并将其存储在  <code>OfflineData<dim>::BoundaryNormalMap</code>  中。我们将在一个后置循环中对这些向量进行归一化处理。

  template <int dim> 
  void OfflineData<dim>::assemble() 
  { 
    lumped_mass_matrix = 0.; 
    norm_matrix        = 0.; 
    for (auto &matrix : cij_matrix) 
      matrix = 0.; 
    for (auto &matrix : nij_matrix) 
      matrix = 0.; 

    unsigned int dofs_per_cell = 
      discretization->finite_element.n_dofs_per_cell(); 
    unsigned int n_q_points = discretization->quadrature.size(); 

// 下面是WorkStream所需的从动数据的初始化过程

    MeshWorker::ScratchData<dim> scratch_data( 
      discretization->mapping, 
      discretization->finite_element, 
      discretization->quadrature, 
      update_values | update_gradients | update_quadrature_points | 
        update_JxW_values, 
      discretization->face_quadrature, 
      update_normal_vectors | update_values | update_JxW_values); 

    { 
      TimerOutput::Scope scope( 
        computing_timer, 
        "offline_data - assemble lumped mass matrix, and c_ij"); 

      const auto local_assemble_system = // 
        [&](const typename DoFHandler<dim>::cell_iterator &cell, 
            MeshWorker::ScratchData<dim> &                 scratch, 
            CopyData<dim> &                                copy) { 
          copy.is_artificial = cell->is_artificial(); 
          if (copy.is_artificial) 
            return; 

          copy.local_boundary_normal_map.clear(); 
          copy.cell_lumped_mass_matrix.reinit(dofs_per_cell, dofs_per_cell); 
          for (auto &matrix : copy.cell_cij_matrix) 
            matrix.reinit(dofs_per_cell, dofs_per_cell); 

          const auto &fe_values = scratch.reinit(cell); 

          copy.local_dof_indices.resize(dofs_per_cell); 
          cell->get_dof_indices(copy.local_dof_indices); 

          std::transform(copy.local_dof_indices.begin(), 
                         copy.local_dof_indices.end(), 
                         copy.local_dof_indices.begin(), 
                         [&](types::global_dof_index index) { 
                           return partitioner->global_to_local(index); 
                         }); 

// 我们以通常的方式计算凑合质量矩阵项 $m_i$ 和向量 $c_{ij}$ 的局部贡献。

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              const auto JxW = fe_values.JxW(q_point); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const auto value_JxW = 
                    fe_values.shape_value(j, q_point) * JxW; 
                  const auto grad_JxW = fe_values.shape_grad(j, q_point) * JxW; 

                  copy.cell_lumped_mass_matrix(j, j) += value_JxW; 

                  for (unsigned int i = 0; i < dofs_per_cell; ++i) 
                    { 
                      const auto value = fe_values.shape_value(i, q_point); 
                      for (unsigned int d = 0; d < dim; ++d) 
                        copy.cell_cij_matrix[d](i, j) += value * grad_JxW[d]; 

                    } /* i */ 


                }     /* j */ 


            }         /* q */ 



// 现在我们要计算边界法线。请注意，除非该元素在域的边界上有面，否则下面的循环不会有什么作用。

          for (const auto f : cell->face_indices()) 
            { 
              const auto face = cell->face(f); 
              const auto id   = face->boundary_id(); 

              if (!face->at_boundary()) 
                continue; 

              const auto &fe_face_values = scratch.reinit(cell, f); 

              const unsigned int n_face_q_points = 
                fe_face_values.get_quadrature().size(); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  if (!discretization->finite_element.has_support_on_face(j, f)) 
                    continue; 

// 注意，"normal "只代表形状函数phi_j支持下的一个面的贡献。所以我们不能在这里对这个局部贡献进行归一化处理，我们必须 "原封不动 "地接受它，存储它并将它传递给复制数据例程。正确的归一化需要在节点上增加一个循环。这在下面的复制函数中完成。

                  Tensor<1, dim> normal; 
                  if (id == Boundaries::free_slip) 
                    { 
                      for (unsigned int q = 0; q < n_face_q_points; ++q) 
                        normal += fe_face_values.normal_vector(q) * 
                                  fe_face_values.shape_value(j, q); 
                    } 

                  const auto index = copy.local_dof_indices[j]; 

                  Point<dim> position; 
                  for (const auto v : cell->vertex_indices()) 
                    if (cell->vertex_dof_index(v, 0) == 
                        partitioner->local_to_global(index)) 
                      { 
                        position = cell->vertex(v); 
                        break; 
                      } 

                  const auto old_id = 
                    std::get<1>(copy.local_boundary_normal_map[index]); 
                  copy.local_boundary_normal_map[index] = 
                    std::make_tuple(normal, std::max(old_id, id), position); 
                } 
            } 
        }; 

// 最后，我们根据WorkStream的要求，提供一个copy_local_to_global函数

      const auto copy_local_to_global = [&](const CopyData<dim> &copy) { 
        if (copy.is_artificial) 
          return; 

        for (const auto &it : copy.local_boundary_normal_map) 
          { 
            std::get<0>(boundary_normal_map[it.first]) += 
              std::get<0>(it.second); 
            std::get<1>(boundary_normal_map[it.first]) = 
              std::max(std::get<1>(boundary_normal_map[it.first]), 
                       std::get<1>(it.second)); 
            std::get<2>(boundary_normal_map[it.first]) = std::get<2>(it.second); 
          } 

        lumped_mass_matrix.add(copy.local_dof_indices, 
                               copy.cell_lumped_mass_matrix); 

        for (int k = 0; k < dim; ++k) 
          { 
            cij_matrix[k].add(copy.local_dof_indices, copy.cell_cij_matrix[k]); 
            nij_matrix[k].add(copy.local_dof_indices, copy.cell_cij_matrix[k]); 
          } 
      }; 

      WorkStream::run(dof_handler.begin_active(), 
                      dof_handler.end(), 
                      local_assemble_system, 
                      copy_local_to_global, 
                      scratch_data, 
                      CopyData<dim>()); 
    } 

// 此时我们已经完成了 $m_i$ 和 $\mathbf{c}_{ij}$ 的计算，但到目前为止，矩阵 <code>nij_matrix</code> 只包含矩阵 <code>cij_matrix</code> 的一个副本。这不是我们真正想要的：我们必须对其条目进行标准化处理。此外，我们还没有填充矩阵 <code>norm_matrix</code> 的条目，存储在映射 <code>OfflineData<dim>::BoundaryNormalMap</code> 中的向量没有被归一化。

// 原则上，这只是离线数据，过度优化它们的计算并没有什么意义，因为它们的成本会在我们将要使用的许多时间步骤中得到摊销。然而，计算/存储矩阵 <code>norm_matrix</code> and the normalization of <code>nij_matrix</code> 的条目是说明线程并行节点循环的最佳方式。

// 我们要访问网格/稀疏图中的每个节点 $i$ 。

// - 对于每一个这样的节点，我们要访问每一个 $j$ ，以便 $\mathbf{c}_{ij} \not \equiv 0$  。

// 从代数的角度来看，这相当于：访问矩阵中的每一行，并对这些行中的每一行在列上执行循环。节点循环是本教程步骤的一个核心主题（见介绍中的伪代码），会反复出现。这就是为什么现在是介绍它们的恰当时机。

// 我们有线程并行化能力 parallel::apply_to_subranges() ，在某种程度上比WorkStream框架更通用。特别是， parallel::apply_to_subranges() 可以用于我们的节点循环。这个功能需要四个输入参数，我们详细解释一下（针对我们的线程并行节点循环的具体案例）。

// - 迭代器  <code>indices.begin()</code>  指向一个行索引。

// - 迭代器 <code>indices.end()</code> 指向一个数字上更高的行索引。

// - 函数 <code>on_subranges(i1,i2)</code> (where <code>i1</code> 和 <code>i2</code> 在前面两个子弹中定义的end和begin迭代器所跨越的范围内定义了一个子范围）对这个子范围内的每个迭代器应用一个操作。我们也可以把 <code>on_subranges</code> 称为 "工作者"。

// - Grainsize：每个线程处理的最小迭代器（在本例中代表行）的数量。我们决定最小为4096行。

// 一个小的注意事项是，提供给 parallel::apply_to_subranges() 的迭代器 <code>indices.begin()</code> 和 <code>indices.end()</code> 必须是随机访问的迭代器：在内部， parallel::apply_to_subranges() 将把 <code>indices.begin()</code> 和 <code>indices.end()</code> 迭代器定义的范围分成子范围（我们希望能够以恒定的复杂性读取这些子范围的任何条目）。为了提供这样的迭代器，我们求助于 std_cxx20::ranges::iota_view. 。

// 下面这段代码的大部分是用来定义 "工作者" <code>on_subranges</code> ：即在子范围的每一行应用的操作。给定一个固定的 <code>row_index</code> ，我们要访问这一行的每一列/每一个条目。为了执行这样的列-循环，我们使用标准库中的<a href="http:www.cplusplus.com/reference/algorithm/for_each/"> std::for_each</a>，其中。

// -  <code>sparsity_pattern.begin(row_index)</code> 给我们一个迭代器，从该行的第一列开始。

// -  <code>sparsity_pattern.end(row_index)</code> 是一个指向该行最后一列的迭代器。

// -  `std::for_each` 所要求的最后一个参数是应用于该行的每个非零条目（本例中为lambda表达式）的操作。

// 我们注意到， parallel::apply_to_subranges() 将对不相交的行集（子行）进行操作，我们的目标是写入这些行中。由于我们要进行的操作的简单性质（法线的计算和存储，以及条目 $\mathbf{c}_{ij}$ 的规范化），线程在试图写同一个条目时不会发生冲突（我们不需要一个调度器）。

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "offline_data - compute |c_ij|, and n_ij"); 

      const std_cxx20::ranges::iota_view<unsigned int, unsigned int> indices( 
        0, n_locally_relevant); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto row_index : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 

// 第一个列循环：我们计算并存储矩阵norm_matrix的条目，并将归一化的条目写入矩阵nij_matrix中。

              std::for_each( 
                sparsity_pattern.begin(row_index), 
                sparsity_pattern.end(row_index), 
                [&](const dealii::SparsityPatternIterators::Accessor &jt) { 
                  const auto   c_ij = gather_get_entry(cij_matrix, &jt); 
                  const double norm = c_ij.norm(); 

                  set_entry(norm_matrix, &jt, norm); 
                  for (unsigned int j = 0; j < dim; ++j) 
                    set_entry(nij_matrix[j], &jt, c_ij[j] / norm); 
                }); 
            } 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 

// 最后，我们对存储在  <code>OfflineData<dim>::BoundaryNormalMap</code>  中的向量进行规范化。这个操作没有被线程并行化，因为它既不能说明任何重要的概念，也不能带来任何明显的速度提升。

      for (auto &it : boundary_normal_map) 
        { 
          auto &normal = std::get<0>(it.second); 
          normal /= (normal.norm() + std::numeric_limits<double>::epsilon()); 
        } 
    } 
  } 

// 在这一点上，我们已经很好地完成了与离线数据有关的事情。

//  @sect4{Equation of state and approximate Riemann solver}  

// 在这一节中，我们描述了 <code>ProblemDescription</code> 类的成员的实现。这里的大部分代码都是针对具有理想气体定律的可压缩欧拉方程的。如果我们想把 step-69 重新用于不同的守恒定律（例如：浅水方程），那么这个类的大部分实现就必须改变。但是其他大部分的类（尤其是那些定义循环结构的类）将保持不变。

// 我们首先实现一些小的成员函数来计算 <code>momentum</code>, <code>internal_energy</code> 、 <code>pressure</code>, <code>speed_of_sound</code> 和系统的通量 <code>f</code> 。这些函数中的每一个的功能都可以从它们的名字中不难看出。

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline Tensor<1, dim> 
  ProblemDescription<dim>::momentum(const state_type &U) 
  { 
    Tensor<1, dim> result; 
    std::copy_n(&U[1], dim, &result[0]); 
    return result; 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::internal_energy(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const auto    m   = momentum(U); 
    const double &E   = U[dim + 1]; 
    return E - 0.5 * m.norm_square() / rho; 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::pressure(const state_type &U) 
  { 
    return (gamma - 1.) * internal_energy(U); 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::speed_of_sound(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const double  p   = pressure(U); 

    return std::sqrt(gamma * p / rho); 
  } 

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline typename ProblemDescription<dim>::flux_type 
  ProblemDescription<dim>::flux(const state_type &U) 
  { 
    const double &rho = U[0]; 
    const auto    m   = momentum(U); 
    const auto    p   = pressure(U); 
    const double &E   = U[dim + 1]; 

    flux_type result; 

    result[0] = m; 
    for (unsigned int i = 0; i < dim; ++i) 
      { 
        result[1 + i] = m * m[i] / rho; 
        result[1 + i][i] += p; 
      } 
    result[dim + 1] = m / rho * (E + p); 

    return result; 
  } 

// 现在我们讨论  $\lambda_{\text{max}} (\mathbf{U}_i^{n},\mathbf{U}_j^{n}, \textbf{n}_{ij})$  的计算。黎曼问题的最大波速的尖锐上界的分析和推导是一个非常技术性的工作，我们不能在本教程中对其进行高级讨论。在这部分文档中，我们将仅限于简述我们实现函数的主要功能，并指出具体的学术参考文献，以帮助（感兴趣的）读者追溯这些想法的来源（和适当的数学证明）。

// 一般来说，要获得最大波速的尖锐保证上界需要解决一个相当昂贵的标量非线性问题。这通常是通过一个迭代求解器来完成的。为了简化本例中的表述，我们决定不包括这样的迭代方案。相反，我们将只是使用一个初始猜测作为最大波速的上限猜测。更确切地说， @cite GuermondPopov2016b 的方程（2.11）（3.7）、（3.8）和（4.3）足以定义最大波速的保证上限。这个估计值通过调用函数  <code>lambda_max_two_rarefaction()</code>  来返回。在其核心部分，这样一个上界的构造使用了所谓的中间压力的二赖式近似  $p^*$  ，例如，见公式（4.46），在  @cite Toro2009  第128页。

// 由 <code>lambda_max_two_rarefaction()</code> 返回的估计值保证是一个上界，它在一般情况下是相当尖锐的，而且对我们的目的来说总体上是足够的。然而，对于一些特定的情况（特别是当其中一个状态接近真空条件时），这样的估计会过于悲观。这就是为什么我们使用第二个估计来避免这种退化，它将通过调用函数  <code>lambda_max_expansion()</code>  来调用。这里最重要的函数是  <code>compute_lambda_max()</code>  ，它取的是  <code>lambda_max_two_rarefaction()</code>  和  <code>lambda_max_expansion()</code>  所返回的估计值之间的最小值。

// 我们再次开始定义几个辅助函数。

// 第一个函数接收一个状态 <code>U</code> 和一个单位向量 <code>n_ij</code> ，并按照单位向量的方向计算<i>projected</i>一维状态。

  namespace 
  { 
    template <int dim> 
    DEAL_II_ALWAYS_INLINE inline std::array<double, 4> riemann_data_from_state( 
      const typename ProblemDescription<dim>::state_type U, 
      const Tensor<1, dim> &                             n_ij) 
    { 
      Tensor<1, 3> projected_U; 
      projected_U[0] = U[0]; 

// 为此，我们必须将动量改为 $\textbf{m}\cdot n_{ij}$ ，并且必须从总能量中减去垂直部分的动能。

      const auto m   = ProblemDescription<dim>::momentum(U); 
      projected_U[1] = n_ij * m; 

      const auto perpendicular_m = m - projected_U[1] * n_ij; 
      projected_U[2] = U[1 + dim] - 0.5 * perpendicular_m.norm_square() / U[0]; 

// 我们以<i>primitive</i>变量而不是守恒量来返回一维状态。返回数组包括密度  $\rho$  、速度  $u$  、压力  $p$  和局部声速  $a$  。

      return {{projected_U[0], 
               projected_U[1] / projected_U[0], 
               ProblemDescription<1>::pressure(projected_U), 
               ProblemDescription<1>::speed_of_sound(projected_U)}}; 
    } 

// 在这一点上，我们还定义了两个小函数，用来返回一个双数的正负部分。

    DEAL_II_ALWAYS_INLINE inline double positive_part(const double number) 
    { 
      return std::max(number, 0.); 
    } 

    DEAL_II_ALWAYS_INLINE inline double negative_part(const double number) 
    { 
      return -std::min(number, 0.); 
    } 

// 接下来，我们需要两个本地文数，它们是以原始状态 $[\rho, u, p, a]$ 和给定压力 $p^\ast$ 为条件定义的。
// @cite GuermondPopov2016  公式（3.7）。
// @f{align*}
//    \lambda^- = u - a\,\sqrt{1 + \frac{\gamma+1}{2\gamma}
//    \left(\frac{p^\ast-p}{p}\right)_+}
//  @f} 
//  这里， $(\cdot)_{+}$  表示给定参数的正数部分。

    DEAL_II_ALWAYS_INLINE inline double 
    lambda1_minus(const std::array<double, 4> &riemann_data, 
                  const double                 p_star) 
    { 

      /* Implements formula (3.7) in Guermond-Popov-2016 */


      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u     = riemann_data[1]; 
      const auto       p     = riemann_data[2]; 
      const auto       a     = riemann_data[3]; 

      const double factor = (gamma + 1.0) / 2.0 / gamma; 
      const double tmp    = positive_part((p_star - p) / p); 
      return u - a * std::sqrt(1.0 + factor * tmp); 
    } 

// Analougously  @cite GuermondPopov2016  方程（3.8）。
// @f{align*}
//    \lambda^+ = u + a\,\sqrt{1 + \frac{\gamma+1}{2\gamma}
//    \left(\frac{p^\ast-p}{p}\right)_+}
//  @f}

    DEAL_II_ALWAYS_INLINE inline double 
    lambda3_plus(const std::array<double, 4> &riemann_data, const double p_star) 
    { 

      /* Implements formula (3.8) in Guermond-Popov-2016 */

      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u     = riemann_data[1]; 
      const auto       p     = riemann_data[2]; 
      const auto       a     = riemann_data[3]; 

      const double factor = (gamma + 1.0) / 2.0 / gamma; 
      const double tmp    = positive_part((p_star - p) / p); 
      return u + a * std::sqrt(1.0 + factor * tmp); 
    } 

// 剩下的就是计算从左和右原始状态计算出来的 $\lambda^-$ 和 $\lambda^+$ 的最大值（  @cite GuermondPopov2016  公式（2.11）），其中 $p^\ast$ 由 @cite GuermondPopov2016  公式（4.3）给出。

    DEAL_II_ALWAYS_INLINE inline double 
    lambda_max_two_rarefaction(const std::array<double, 4> &riemann_data_i, 
                               const std::array<double, 4> &riemann_data_j) 
    { 
      constexpr double gamma = ProblemDescription<1>::gamma; 
      const auto       u_i   = riemann_data_i[1]; 
      const auto       p_i   = riemann_data_i[2]; 
      const auto       a_i   = riemann_data_i[3]; 
      const auto       u_j   = riemann_data_j[1]; 
      const auto       p_j   = riemann_data_j[2]; 
      const auto       a_j   = riemann_data_j[3]; 

      const double numerator = a_i + a_j - (gamma - 1.) / 2. * (u_j - u_i); 

      const double denominator = 
        a_i * std::pow(p_i / p_j, -1. * (gamma - 1.) / 2. / gamma) + a_j * 1.; 

/* Guermond-Popov-2016 */ 

中的公式（4.3）。
      const double p_star = 
        p_j * std::pow(numerator / denominator, 2. * gamma / (gamma - 1)); 

      const double lambda1 = lambda1_minus(riemann_data_i, p_star); 
      const double lambda3 = lambda3_plus(riemann_data_j, p_star); 

/* Guermond-Popov-2016中的公式（2.11）  */ 

      return std::max(positive_part(lambda3), negative_part(lambda1)); 
    } 

// 我们计算出最大波速的第二个上界，一般来说，它不像二重化估计那样尖锐。但在接近真空的条件下，当二赖子近似值可能达到极端值时，它将挽救一切。
// @f{align*}
//    \lambda_{\text{exp}} = \max(u_i,u_j) + 5. \max(a_i, a_j).
//  @f} 
// @note  常数5.0乘以声速的最大值是<i>neither</i>一个临时的常数，<i>nor</i>一个调整参数。它为任何  $\gamma \in [0,5/3]$  定义了一个上限。请不要玩弄它!

    DEAL_II_ALWAYS_INLINE inline double 
    lambda_max_expansion(const std::array<double, 4> &riemann_data_i, 
                         const std::array<double, 4> &riemann_data_j) 
    { 
      const auto u_i = riemann_data_i[1]; 
      const auto a_i = riemann_data_i[3]; 
      const auto u_j = riemann_data_j[1]; 
      const auto a_j = riemann_data_j[3]; 

      return std::max(std::abs(u_i), std::abs(u_j)) + 5. * std::max(a_i, a_j); 
    } 
  } // namespace 

// 下面是我们要调用的主函数，以计算  $\lambda_{\text{max}} (\mathbf{U}_i^{n},\mathbf{U}_j^{n}, \textbf{n}_{ij})$  。我们简单地计算两个最大的波速估计值并返回最小值。

  template <int dim> 
  DEAL_II_ALWAYS_INLINE inline double 
  ProblemDescription<dim>::compute_lambda_max(const state_type &    U_i, 
                                              const state_type &    U_j, 
                                              const Tensor<1, dim> &n_ij) 
  { 
    const auto riemann_data_i = riemann_data_from_state(U_i, n_ij); 
    const auto riemann_data_j = riemann_data_from_state(U_j, n_ij); 

    const double lambda_1 = 
      lambda_max_two_rarefaction(riemann_data_i, riemann_data_j); 

    const double lambda_2 = 
      lambda_max_expansion(riemann_data_i, riemann_data_j); 

    return std::min(lambda_1, lambda_2); 
  } 

// 我们通过定义静态数组 <code>component_names</code> 来结束本节，这些静态数组包含描述我们的状态向量的组件名称的字符串。我们对维度一、二和三进行了模板特化，这在后面的DataOut中被用来命名相应的组件。

  template <> 
  const std::array<std::string, 3> ProblemDescription<1>::component_names{ 
    {"rho", "m", "E"}}; 

  template <> 
  const std::array<std::string, 4> ProblemDescription<2>::component_names{ 
    {"rho", "m_1", "m_2", "E"}}; 

  template <> 
  const std::array<std::string, 5> ProblemDescription<3>::component_names{ 
    {"rho", "m_1", "m_2", "m_3", "E"}}; 
// @sect4{Initial values}  

// 在我们讨论正向欧拉方案的实现之前，最后一个准备步骤是简单地实现`InitialValues`类。

// 在构造函数中，我们用默认值初始化所有参数，为`参数接受器`类声明所有参数，并将 <code>parse_parameters_call_back</code> 槽连接到相应的信号。

//  <code>parse_parameters_call_back</code> 槽将在调用 ParameterAcceptor::initialize(). 后从ParameterAceptor中调用。 在这方面，它的使用适合于参数必须被后处理（在某种意义上）或必须检查参数之间的某些一致性条件的情况。

  template <int dim> 
  InitialValues<dim>::InitialValues(const std::string &subsection) 
    : ParameterAcceptor(subsection) 
  { 

    /* We wire up the slot InitialValues<dim>::parse_parameters_callback to
       the ParameterAcceptor::parse_parameters_call_back signal: */


    ParameterAcceptor::parse_parameters_call_back.connect( 
      std::bind(&InitialValues<dim>::parse_parameters_callback, this)); 

    initial_direction[0] = 1.; 
    add_parameter("initial direction", 
                  initial_direction, 
                  "Initial direction of the uniform flow field"); 

    initial_1d_state[0] = ProblemDescription<dim>::gamma; 
    initial_1d_state[1] = 3.; 
    initial_1d_state[2] = 1.; 
    add_parameter("initial 1d state", 
                  initial_1d_state, 
                  "Initial 1d state (rho, u, p) of the uniform flow field"); 
  } 

// 到目前为止， <code>InitialValues</code> 的构造函数已经为两个私有成员 <code>initial_direction</code> and <code>initial_1d_state</code> 定义了默认值，并将它们添加到参数列表中。但是我们还没有定义我们真正关心的唯一公共成员的实现，也就是 <code>initial_state()</code> （我们将调用这个函数来实际评估网格节点的初始解）。在该函数的顶部，我们必须确保提供的初始方向不是零矢量。

//  @note  正如所评论的，我们可以避免使用方法  <code>parse_parameters_call_back </code>  并定义一个类成员  <code>setup()</code>  以便定义  <code>initial_state()</code>  的实现。但为了说明问题，我们想在这里记录一种不同的方式，并使用ParameterAcceptor的回调信号。

  template <int dim> 
  void InitialValues<dim>::parse_parameters_callback() 
  { 
    AssertThrow(initial_direction.norm() != 0., 
                ExcMessage( 
                  "Initial shock front direction is set to the zero vector.")); 
    initial_direction /= initial_direction.norm(); 

// 接下来，我们用一个计算均匀流场的lambda函数来实现 <code>initial_state</code> 函数对象。为此，我们必须将给定的原始1d状态（密度 $\rho$ 、速度 $u$ 和压力 $p$ ）转换为保守的n维状态（密度 $\rho$ 、动量 $\mathbf{m}$ 和总能量 $E$  ）。

    initial_state = [this](const Point<dim> & /*point*/, double /*t*/) { 
      const double            rho   = initial_1d_state[0]; 
      const double            u     = initial_1d_state[1]; 
      const double            p     = initial_1d_state[2]; 
      static constexpr double gamma = ProblemDescription<dim>::gamma; 

      state_type state; 

      state[0] = rho; 
      for (unsigned int i = 0; i < dim; ++i) 
        state[1 + i] = rho * u * initial_direction[i]; 

      state[dim + 1] = p / (gamma - 1.) + 0.5 * rho * u * u; 

      return state; 
    }; 
  } 
// @sect4{The Forward Euler step}  

//  <code>%TimeStepping</code> 类的构造函数不包含任何令人惊讶的代码。

  template <int dim> 
  TimeStepping<dim>::TimeStepping( 
    const MPI_Comm            mpi_communicator, 
    TimerOutput &             computing_timer, 
    const OfflineData<dim> &  offline_data, 
    const InitialValues<dim> &initial_values, 
    const std::string &       subsection /*= "TimeStepping"*/) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , offline_data(&offline_data) 
    , initial_values(&initial_values) 
  { 
    cfl_update = 0.80; 
    add_parameter("cfl update", 
                  cfl_update, 
                  "Relative CFL constant used for update"); 
  } 

// 在类成员  <code>prepare()</code>  中我们初始化了临时向量  <code>temp</code> and the matrix <code>dij_matrix</code>  。该向量 <code>temp</code> 将在其内容与旧向量交换之前用于临时存储解决方案的更新。

  template <int dim> 
  void TimeStepping<dim>::prepare() 
  { 
    TimerOutput::Scope scope(computing_timer, 
                             "time_stepping - prepare scratch space"); 

    for (auto &it : temporary_vector) 
      it.reinit(offline_data->partitioner); 

    dij_matrix.reinit(offline_data->sparsity_pattern); 
  } 

// 现在是实现正向欧拉步骤的时候了。给出一个在时间 $t$ 的旧状态 <code>U</code> 的（可写引用），我们就地更新状态 <code>U</code> ，并返回所选择的时间步长。我们首先声明一些对各种不同变量和数据结构的只读引用。我们这样做主要是为了有更短的变量名称（例如， <code>sparsity</code> 而不是 <code>offline_data->sparsity_pattern</code> ）。

  template <int dim> 
  double TimeStepping<dim>::make_one_step(vector_type &U, double t) 
  { 
    const auto &n_locally_owned    = offline_data->n_locally_owned; 
    const auto &n_locally_relevant = offline_data->n_locally_relevant; 

    const std_cxx20::ranges::iota_view<unsigned int, unsigned int> 
      indices_owned(0, n_locally_owned); 
    const std_cxx20::ranges::iota_view<unsigned int, unsigned int> 
      indices_relevant(0, n_locally_relevant); 

    const auto &sparsity = offline_data->sparsity_pattern; 

    const auto &lumped_mass_matrix = offline_data->lumped_mass_matrix; 
    const auto &norm_matrix        = offline_data->norm_matrix; 
    const auto &nij_matrix         = offline_data->nij_matrix; 
    const auto &cij_matrix         = offline_data->cij_matrix; 

    const auto &boundary_normal_map = offline_data->boundary_normal_map; 
//<b>Step 1</b>: 计算 $d_{ij}$ 图的粘性矩阵。

// 需要强调的是，粘度矩阵必须是对称的，即  $d_{ij} = d_{ji}$  。在这方面我们注意到， $\int_{\Omega} \nabla \phi_j \phi_i \, \mathrm{d}\mathbf{x}= -
//  \int_{\Omega} \nabla \phi_i \phi_j \, \mathrm{d}\mathbf{x}$ （或等同于 $\mathbf{c}_{ij} = - \mathbf{c}_{ji}$ ）提供了 $\mathbf{x}_i$ 或 $\mathbf{x}_j$ 是一个位于远离边界的支持点。在这种情况下，我们可以通过构造检查出 $\lambda_{\text{max}} (\mathbf{U}_i^{n}, \mathbf{U}_j^{n},
//  \textbf{n}_{ij}) = \lambda_{\text{max}} (\mathbf{U}_j^{n},
//  \mathbf{U}_i^{n},\textbf{n}_{ji})$ ，这保证了 $d_{ij} = d_{ji}$ 的属性。

// 然而，如果两个支持点 $\mathbf{x}_i$ 或 $\mathbf{x}_j$ 恰好都位于边界上，那么，等式 $\mathbf{c}_{ij} =
// - \mathbf{c}_{ji}$ 和 $\lambda_{\text{max}} (\mathbf{U}_i^{n},
//  \mathbf{U}_j^{n}, \textbf{n}_{ij}) = \lambda_{\text{max}}
//  (\mathbf{U}_j^{n}, \mathbf{U}_i^{n}, \textbf{n}_{ji})$ 就不一定成立。对于这个难题，数学上唯一安全的解决方案是计算 $d_{ij}$ 和 $d_{ji}$ ，并取其最大值。

// 总体而言， $d_{ij}$ 的计算是相当昂贵的。为了节省一些计算时间，我们利用了粘度矩阵必须是对称的这一事实（如上所述）：我们只计算 $d_{ij}$ 的上三角条目，并将相应的条目复制到下三角的对应项上。

// 我们再次使用 parallel::apply_to_subranges() 来实现线程并行的for loops。我们在讨论矩阵的组装 <code>norm_matrix</code> 和上面 <code>nij_matrix</code> 的归一化时介绍的几乎所有并行遍历的想法都在这里得到了应用。

// 我们再次定义了一个 "工作者 "函数 <code>on_subranges</code> ，计算列索引子范围[i1, i2]的黏度 $d_{ij}$ 。

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 1 compute d_ij"); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 
              const auto U_i = gather(U, i); 

// 对于一个给定的列索引i，我们遍历从 <code>sparsity.begin(i)</code> 到 <code>sparsity.end(i)</code> 的稀疏模式的列。

              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 

// 我们只计算 $d_{ij}$ ，如果 $j < i$ （上三角条目），随后将数值复制到 $d_{ji}$  。

                  if (j >= i) 
                    continue; 

                  const auto U_j = gather(U, j); 

                  const auto   n_ij = gather_get_entry(nij_matrix, jt); 
                  const double norm = get_entry(norm_matrix, jt); 

                  const auto lambda_max = 
                    ProblemDescription<dim>::compute_lambda_max(U_i, U_j, n_ij); 

                  double d = norm * lambda_max; 

// 如果两个支持点刚好都在边界上，我们也要计算 $d_{ji}$ ，然后再取 $\max(d_{ij},d_{ji})$  。在这之后，我们可以最终设定上三角和下三角的条目。

                  if (boundary_normal_map.count(i) != 0 && 
                      boundary_normal_map.count(j) != 0) 
                    { 
                      const auto n_ji = gather(nij_matrix, j, i); 
                      const auto lambda_max_2 = 
                        ProblemDescription<dim>::compute_lambda_max(U_j, 
                                                                    U_i, 
                                                                    n_ji); 
                      const double norm_2 = norm_matrix(j, i); 

                      d = std::max(d, norm_2 * lambda_max_2); 
                    } 

                  set_entry(dij_matrix, jt, d); 
                  dij_matrix(j, i) = d; 
                } 
            } 
        }; 

      parallel::apply_to_subranges(indices_relevant.begin(), 
                                   indices_relevant.end(), 
                                   on_subranges, 
                                   4096); 
    } 
//<b>Step 2</b>: 计算对角线项  $d_{ii}$  和  $\tau_{\text{max}}$  。

// 到目前为止，我们已经计算了矩阵 <code>dij_matrix</code> 的所有非对角线项。我们仍然需要填补其对角线项，定义为  $d_{ii}^n = - \sum_{j \in \mathcal{I}(i)\backslash \{i\}} d_{ij}^n$  。我们再次使用 parallel::apply_to_subranges() 来实现这一目的。在计算 $d_{ii}$ s的同时，我们也确定了最大的可接受的时间步长，定义为
// \f[
//    \tau_n \dealcoloneq c_{\text{cfl}}\,\min_{i\in\mathcal{V}}
//    \left(\frac{m_i}{-2\,d_{ii}^{n}}\right) \, .
//  \f] 
//  注意， $\min_{i \in \mathcal{V}}$ 的操作本质上是全局的，它在所有节点上操作：首先我们必须在所有线程（特定节点的）上取最小值，然后我们必须在所有MPI进程上取最小值。在目前的实现中。

// - 我们将 <code>tau_max</code> （每个节点）存储为<a href="http:www.cplusplus.com/reference/atomic/atomic/"><code>std::atomic<double></code></a>。   <code>std::atomic</code> 的内部实现将在一个以上的线程试图同时读取和/或写入 <code>tau_max</code> 时，负责保护任何可能的竞赛条件。

// - 为了取所有MPI进程的最小值，我们使用实用函数  <code>Utilities::MPI::min</code>  。

    std::atomic<double> tau_max{std::numeric_limits<double>::infinity()}; 

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 2 compute d_ii, and tau_max"); 

// on_subranges()将在每个线程上单独执行。因此，变量 <code>tau_max_on_subrange</code> 被存储在线程本地。

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          double tau_max_on_subrange = std::numeric_limits<double>::infinity(); 

          for (const auto i : 
               std_cxx20::ranges::iota_view<unsigned int, unsigned int>(*i1, 
                                                                        *i2)) 
            { 
              double d_sum = 0.; 

              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 

                  if (j == i) 
                    continue; 

                  d_sum -= get_entry(dij_matrix, jt); 
                } 

// 我们将d_ij项的负数之和存储在对角线的位置上。

              dij_matrix.diag_element(i) = d_sum; 

// 并计算出最大的局部时间步长  <code>tau</code>  。

              const double mass   = lumped_mass_matrix.diag_element(i); 
              const double tau    = cfl_update * mass / (-2. * d_sum); 
              tau_max_on_subrange = std::min(tau_max_on_subrange, tau); 
            } 
// <code>tau_max_on_subrange</code>  包含为（线程局部）子范围计算的最大可能的时间步长。在这一点上，我们必须在所有线程上同步该值。这就是我们使用<a
//  href="http:www.cplusplus.com/reference/atomic/atomic/"><code>std::atomic<double></code></a> 的原因。
//<i>compare exchange</i> 更新机制。

          double current_tau_max = tau_max.load(); 
          while (current_tau_max > tau_max_on_subrange && 
                 !tau_max.compare_exchange_weak(current_tau_max, 
                                                tau_max_on_subrange)) 
            ; 
        }; 

      parallel::apply_to_subranges(indices_relevant.begin(), 
                                   indices_relevant.end(), 
                                   on_subranges, 
                                   4096); 

// 在所有线程完成后，我们可以简单地在所有MPI进程中同步该值。

      tau_max.store(Utilities::MPI::min(tau_max.load(), mpi_communicator)); 

// 这是一个验证计算出的 <code>tau_max</code> 确实是一个有效浮点数的好时机。

      AssertThrow( 
        !std::isnan(tau_max.load()) && !std::isinf(tau_max.load()) && 
          tau_max.load() > 0., 
        ExcMessage( 
          "I'm sorry, Dave. I'm afraid I can't do that. - We crashed.")); 
    } 
//<b>Step 3</b>: 执行更新。

// 在这一点上，我们已经计算了所有的粘性系数  $d_{ij}$  并且我们知道最大的可接受的时间步长  $\tau_{\text{max}}$  。这意味着我们现在可以计算更新了。

// \f[
   \mathbf{U}_i^{n+1} = \mathbf{U}_i^{n} - \frac{\tau_{\text{max}} }{m_i}
   \sum_{j \in \mathcal{I}(i)} (\mathbb{f}(\mathbf{U}_j^{n}) -
   \mathbb{f}(\mathbf{U}_i^{n})) \cdot \mathbf{c}_{ij} - d_{ij}
   (\mathbf{U}_j^{n} - \mathbf{U}_i^{n})
 \f]

// 这个更新公式与介绍中讨论的略有不同（在伪代码中）。然而，可以证明这两个公式在代数上是等价的（它们将产生相同的数值）。我们更倾向于第二个公式，因为它具有自然的取消属性，可能有助于避免数字上的伪影。

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 3 perform update"); 

      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              const auto U_i = gather(U, i); 

              const auto   f_i = ProblemDescription<dim>::flux(U_i); 
              const double m_i = lumped_mass_matrix.diag_element(i); 

              auto U_i_new = U_i; 

              for (auto jt = sparsity.begin(i); jt != sparsity.end(i); ++jt) 
                { 
                  const auto j = jt->column(); 

                  const auto U_j = gather(U, j); 
                  const auto f_j = ProblemDescription<dim>::flux(U_j); 

                  const auto c_ij = gather_get_entry(cij_matrix, jt); 
                  const auto d_ij = get_entry(dij_matrix, jt); 

                  for (unsigned int k = 0; k < problem_dimension; ++k) 
                    { 
                      U_i_new[k] += 
                        tau_max / m_i * 
                        (-(f_j[k] - f_i[k]) * c_ij + d_ij * (U_j[k] - U_i[k])); 
                    } 
                } 

              scatter(temporary_vector, U_i_new, i); 
            } 
        }; 

      parallel::apply_to_subranges(indices_owned.begin(), 
                                   indices_owned.end(), 
                                   on_subranges, 
                                   4096); 
    } 
//<b>Step 4</b>: 修复了边界状态。

// 作为正向欧拉方法的最后一步，我们必须修复所有的边界状态。正如在介绍中所讨论的，我们

// 在完全不满足边界条件的情况下进行时间推进。

// -- 在时间步长结束时，在后处理步骤中强力执行边界条件。

// 在这里，我们计算修正\f[
//    \mathbf{m}_i \dealcoloneq \mathbf{m}_i - (\boldsymbol{\nu}_i \cdot
//    \mathbf{m}_i) \boldsymbol{\nu}_i,
//  \f]，它消除了 $\mathbf{m}$ 的法线成分。

    { 
      TimerOutput::Scope scope(computing_timer, 
                               "time_stepping - 4 fix boundary states"); 

      for (auto it : boundary_normal_map) 
        { 
          const auto i = it.first; 

// 我们只对本地拥有的子集进行迭代。

          if (i >= n_locally_owned) 
            continue; 

          const auto &normal   = std::get<0>(it.second); 
          const auto &id       = std::get<1>(it.second); 
          const auto &position = std::get<2>(it.second); 

          auto U_i = gather(temporary_vector, i); 

// 在自由滑移的边界上，我们去除动量的法向分量。

          if (id == Boundaries::free_slip) 
            { 
              auto m = ProblemDescription<dim>::momentum(U_i); 
              m -= (m * normal) * normal; 
              for (unsigned int k = 0; k < dim; ++k) 
                U_i[k + 1] = m[k]; 
            } 

// 在Dirichlet边界上，我们强行执行初始条件。

          else if (id == Boundaries::dirichlet) 
            { 
              U_i = initial_values->initial_state(position, t + tau_max); 
            } 

          scatter(temporary_vector, U_i, i); 
        } 
    } 
//<b>Step 5</b>: 我们现在在所有MPI行列上更新幽灵层，将临时向量与解决方案向量交换  <code>U</code>  （将通过引用返回），并返回选择的时间步长  $\tau_{\text{max}}$  。

    for (auto &it : temporary_vector) 
      it.update_ghost_values(); 

    U.swap(temporary_vector); 

    return tau_max; 
  } 
// @sect4{Schlieren postprocessing}  

// 在不同的时间间隔内，我们将输出解决方案的当前状态 <code>U</code> 以及所谓的Schlieren图。 <code>SchlierenPostprocessor</code> 类的构造函数同样不包含任何惊喜。我们只是提供默认值并注册两个参数。

// - schlieren_beta: 是一个临时的正向放大系数，以增强可视化中的对比度。它的实际值是一个品味问题。

// - schlieren_index: 是一个整数，表示我们将使用状态 $[\rho, \mathbf{m},E]$ 中的哪个组件来生成可视化。

  template <int dim> 
  SchlierenPostprocessor<dim>::SchlierenPostprocessor( 
    const MPI_Comm          mpi_communicator, 
    TimerOutput &           computing_timer, 
    const OfflineData<dim> &offline_data, 
    const std::string &     subsection /*= "SchlierenPostprocessor"*/) 
    : ParameterAcceptor(subsection) 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(computing_timer) 
    , offline_data(&offline_data) 
  { 
    schlieren_beta = 10.; 
    add_parameter("schlieren beta", 
                  schlieren_beta, 
                  "Beta factor used in Schlieren-type postprocessor"); 

    schlieren_index = 0; 
    add_parameter("schlieren index", 
                  schlieren_index, 
                  "Use the corresponding component of the state vector for the " 
                  "schlieren plot"); 
  } 

// 同样， <code>prepare()</code> 函数初始化了两个临时向量（  <code>r</code> and <code>schlieren</code>  ）。

  template <int dim> 
  void SchlierenPostprocessor<dim>::prepare() 
  { 
    TimerOutput::Scope scope(computing_timer, 
                             "schlieren_postprocessor - prepare scratch space"); 

    r.reinit(offline_data->n_locally_relevant); 
    schlieren.reinit(offline_data->partitioner); 
  } 

// 我们现在讨论类成员 <code>SchlierenPostprocessor<dim>::compute_schlieren()</code> 的实现，它基本上是取状态向量 <code>U</code> 的一个分量并计算该分量的Schlieren指标（Schlieren指标的公式可以在类的声明 <code>SchlierenPostprocessor</code> 之前找到）。我们首先注意到这个公式需要 "结点梯度"  $\nabla r_j$  。然而，对于  $\mathcal{C}^0$  有限元函数来说，梯度的节点值并没有定义。更为普遍的是，梯度的点值对于 $W^{1,p}(\Omega)$ 函数没有定义。我们可以用最简单的技术来恢复节点的梯度，即加权平均法。

// \f[ \nabla r_j \dealcoloneq \frac{1}{\int_{S_i} \omega_i(\mathbf{x}) \,
//  \mathrm{d}\mathbf{x}}
//   \int_{S_i} r_h(\mathbf{x}) \omega_i(\mathbf{x}) \, \mathrm{d}\mathbf{x}
//  \ \ \ \ \ \mathbf{(*)} \f]

// 其中 $S_i$ 是形状函数 $\phi_i$ 的支持，而 $\omega_i(\mathbf{x})$ 是权重。权重可以是任何正函数，如 $\omega_i(\mathbf{x}) \equiv 1$ （这将使我们恢复通常的均值概念）。但是像往常一样，我们的目标是尽可能多地重复使用离线数据。在这个意义上，最自然的权重选择是 $\omega_i = \phi_i$  。将这种权重的选择和扩展 $r_h(\mathbf{x}) = \sum_{j \in \mathcal{V}} r_j \phi_j(\mathbf{x})$ 插入 $\mathbf{(*)}$ 中，我们得到:

//  \f[
//  \nabla r_j \dealcoloneq \frac{1}{m_i} \sum_{j \in \mathcal{I}(i)} r_j
// \mathbf{c}_{ij} \ \ \ \ \ \mathbf{(**)} \, . 
//  \f]

// 使用这最后一个公式，我们可以恢复平均的节点梯度，而不需要借助任何形式的正交。这个想法与基于边缘的方案（或代数方案）的整体精神非常吻合，我们希望尽可能直接对矩阵和向量进行操作，以避免使用双线性形式、单元环、正交，或在输入参数（上一时间步的状态）和计算更新所需的实际矩阵和向量之间的任何其他中间结构/操作。

// 第二件要注意的事情是，我们必须计算全局最小和最大  $\max_j |\nabla r_j|$  和  $\min_j |\nabla r_j|$  。按照在类成员 <code>%TimeStepping\<dim>::%step()</code> 中用于计算时间步长的相同思路，我们将 $\max_j |\nabla r_j|$ 和 $\min_j |\nabla r_j|$ 定义为原子双数，以解决线程之间的任何冲突。像往常一样，我们使用 <code>Utilities::MPI::max()</code> 和 <code>Utilities::MPI::min()</code> 来寻找所有MPI进程中的全局最大/最小值。

// 最后，不可能在所有节点上单次循环计算Schlieren指标。整个操作需要在节点上进行两次循环。



// - 第一个循环对网格中所有的 $|\nabla r_i|$ 进行计算，并计算边界 $\max_j |\nabla r_j|$ 和 $\min_j |\nabla r_j|$  。

// - 第二个循环最后用公式计算Schlieren指标

// \f[ \text{schlieren}[i] = e^{\beta \frac{ |\nabla r_i|
//  - \min_j |\nabla r_j| }{\max_j |\nabla r_j| - \min_j |\nabla r_j| } }
//  \, . 
//  \f]

// 这意味着我们将不得不为每一个阶段定义两个工作者 <code>on_subranges</code> 。

  template <int dim> 
  void SchlierenPostprocessor<dim>::compute_schlieren(const vector_type &U) 
  { 
    TimerOutput::Scope scope( 
      computing_timer, "schlieren_postprocessor - compute schlieren plot"); 

    const auto &sparsity            = offline_data->sparsity_pattern; 
    const auto &lumped_mass_matrix  = offline_data->lumped_mass_matrix; 
    const auto &cij_matrix          = offline_data->cij_matrix; 
    const auto &boundary_normal_map = offline_data->boundary_normal_map; 
    const auto &n_locally_owned     = offline_data->n_locally_owned; 

    const auto indices = 
      std_cxx20::ranges::iota_view<unsigned int, unsigned int>(0, 
                                                               n_locally_owned); 

// 我们将当前MPI进程中的r_i_max和r_i_min定义为原子倍数，以避免线程之间的竞赛条件。

    std::atomic<double> r_i_max{0.}; 
    std::atomic<double> r_i_min{std::numeric_limits<double>::infinity()}; 

// 第一个循环：计算每个节点的平均梯度以及梯度的全局最大值和最小值。

    { 
      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          double r_i_max_on_subrange = 0.; 
          double r_i_min_on_subrange = std::numeric_limits<double>::infinity(); 

          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              Tensor<1, dim> r_i; 

 
                { 
                  const auto j = jt->column(); 

                  if (i == j) 
                    continue; 

                  const auto U_js = U[schlieren_index].local_element(j); 
                  const auto c_ij = gather_get_entry(cij_matrix, jt); 
                  r_i += c_ij * U_js; 
                } 

// 我们在自由滑移边界固定梯度r_i，类似于我们在正向欧拉步骤中固定边界状态的方式。    这样可以避免在自由滑移边界的Schlieren图中出现尖锐的、人为的梯度，这纯粹是一种外观上的选择。

              const auto bnm_it = boundary_normal_map.find(i); 
              if (bnm_it != boundary_normal_map.end()) 
                { 
                  const auto &normal = std::get<0>(bnm_it->second); 
                  const auto &id     = std::get<1>(bnm_it->second); 

                  if (id == Boundaries::free_slip) 
                    r_i -= 1. * (r_i * normal) * normal; 
                  else 
                    r_i = 0.; 
                } 

// 我们提醒读者，我们对结点梯度本身并不感兴趣。我们只想得到它们的规范，以便计算Schlieren指标（用块状质量矩阵 $m_i$  加权）。

              const double m_i    = lumped_mass_matrix.diag_element(i); 
              r[i]                = r_i.norm() / m_i; 
              r_i_max_on_subrange = std::max(r_i_max_on_subrange, r[i]); 
              r_i_min_on_subrange = std::min(r_i_min_on_subrange, r[i]); 
            } 

// 我们将current_r_i_max和current_r_i_min（在当前子范围内）与r_i_max和r_i_min（对于当前MPI进程）进行比较，并在必要时进行更新。

          double current_r_i_max = r_i_max.load(); 
          while (current_r_i_max < r_i_max_on_subrange && 
                 !r_i_max.compare_exchange_weak(current_r_i_max, 
                                                r_i_max_on_subrange)) 
            ; 

          double current_r_i_min = r_i_min.load(); 
          while (current_r_i_min > r_i_min_on_subrange && 
                 !r_i_min.compare_exchange_weak(current_r_i_min, 
                                                r_i_min_on_subrange)) 
            ; 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 
    } 

// 在所有MPI进程中同步 <code>r_i_max</code> and <code>r_i_min</code> 。

    r_i_max.store(Utilities::MPI::max(r_i_max.load(), mpi_communicator)); 
    r_i_min.store(Utilities::MPI::min(r_i_min.load(), mpi_communicator)); 

// 第二个循环：我们现在有了矢量 <code>r</code> 和标量 <code>r_i_max</code> and <code>r_i_min</code> 可以使用。这样我们就可以实际计算Schlieren指标了。

    { 
      const auto on_subranges = // 
        [&](const auto i1, const auto i2) { 
          for (const auto i : boost::make_iterator_range(i1, i2)) 
            { 
              Assert(i < n_locally_owned, ExcInternalError()); 

              schlieren.local_element(i) = 
                1. - std::exp(-schlieren_beta * (r[i] - r_i_min) / 
                              (r_i_max - r_i_min)); 
            } 
        }; 

      parallel::apply_to_subranges(indices.begin(), 
                                   indices.end(), 
                                   on_subranges, 
                                   4096); 
    } 

// 最后，交换幽灵元素。

    schlieren.update_ghost_values(); 
  } 
// @sect4{The main loop}  

// 在实现了所有的类之后，是时候创建一个 <code>Discretization<dim></code>, <code>OfflineData<dim></code> 、 <code>InitialValues<dim></code>, <code>%TimeStepping\<dim></code> 和 <code>SchlierenPostprocessor<dim></code> 的实例，并在一个循环中运行欧拉正步。

// 在 <code>MainLoop<dim></code> 的构造函数中，我们现在初始化所有类的实例，并声明一些控制输出的参数。最值得注意的是，我们声明了一个布尔参数 <code>resume</code> ，它将控制程序是否试图从中断的计算中重新启动。

  template <int dim> 
  MainLoop<dim>::MainLoop(const MPI_Comm mpi_communicator) 
    : ParameterAcceptor("A - MainLoop") 
    , mpi_communicator(mpi_communicator) 
    , computing_timer(mpi_communicator, 
                      timer_output, 
                      TimerOutput::never, 
                      TimerOutput::cpu_and_wall_times) 
    , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
    , discretization(mpi_communicator, computing_timer, "B - Discretization") 
    , offline_data(mpi_communicator, 
                   computing_timer, 
                   discretization, 
                   "C - OfflineData") 
    , initial_values("D - InitialValues") 
    , time_stepping(mpi_communicator, 
                    computing_timer, 
                    offline_data, 
                    initial_values, 
                    "E - TimeStepping") 
    , schlieren_postprocessor(mpi_communicator, 
                              computing_timer, 
                              offline_data, 
                              "F - SchlierenPostprocessor") 
  { 
    base_name = "test"; 
    add_parameter("basename", base_name, "Base name for all output files"); 

    t_final = 4.; 
    add_parameter("final time", t_final, "Final time"); 

    output_granularity = 0.02; 
    add_parameter("output granularity", 
                  output_granularity, 
                  "time interval for output"); 

    asynchronous_writeback = true; 
    add_parameter("asynchronous writeback", 
                  asynchronous_writeback, 
                  "Write out solution in a background thread performing IO"); 

    resume = false; 
    add_parameter("resume", resume, "Resume an interrupted computation."); 
  } 

// 我们首先在匿名命名空间中实现一个辅助函数 <code>print_head()</code> ，用来在终端输出带有一些漂亮格式的信息。

  namespace 
  { 
    void print_head(ConditionalOStream &pcout, 
                    const std::string & header, 
                    const std::string & secondary = "") 
    { 
      const auto header_size   = header.size(); 
      const auto padded_header = std::string((34 - header_size) / 2, ' ') + 
                                 header + 
                                 std::string((35 - header_size) / 2, ' '); 

      const auto secondary_size = secondary.size(); 
      const auto padded_secondary = 
        std::string((34 - secondary_size) / 2, ' ') + secondary + 
        std::string((35 - secondary_size) / 2, ' '); 

 /* 关闭clang-format  */ 

      pcout << std::endl; 
      pcout << "    ####################################################" << std::endl; 
      pcout << "    #########                                  #########" << std::endl; 
      pcout << "    #########"     <<  padded_header   <<     "#########" << std::endl; 
      pcout << "    #########"     << padded_secondary <<     "#########" << std::endl; 
      pcout << "    #########                                  #########" << std::endl; 
      pcout << "    ####################################################" << std::endl; 
      pcout << std::endl; 
    /* clang-format on  */ 
    } 
  } // namespace 

// 有了 <code>print_head</code> ，现在是时候实现 <code>MainLoop<dim>::run()</code> 了，它包含了我们程序的主循环。

  template <int dim> 
  void MainLoop<dim>::run() 
  { 

// 我们开始读入参数并初始化所有对象。我们在这里注意到，对 ParameterAcceptor::initialize 的调用是从参数文件（其名称作为一个字符串参数给出）中读入所有参数。ParameterAcceptor处理一个全局的ParameterHandler，它被初始化为所有从ParameterAceptor派生的类实例的子节和参数声明。调用initialize进入每个每个派生类的分节，并设置所有使用 ParameterAcceptor::add_parameter() 添加的变量。
    pcout << "Reading parameters and allocating objects... " << std::flush; 

    ParameterAcceptor::initialize("step-69.prm"); 
    pcout << "done" << std::endl; 

// 接下来我们创建三角形，集合所有的矩阵，设置划痕空间，并初始化DataOut<dim>对象。

    { 
      print_head(pcout, "create triangulation"); 
      discretization.setup(); 

      pcout << "Number of active cells:       " 
            << discretization.triangulation.n_global_active_cells() 
            << std::endl; 

      print_head(pcout, "compute offline data"); 
      offline_data.setup(); 
      offline_data.assemble(); 

      pcout << "Number of degrees of freedom: " 
            << offline_data.dof_handler.n_dofs() << std::endl; 

      print_head(pcout, "set up time step"); 
      time_stepping.prepare(); 
      schlieren_postprocessor.prepare(); 
    } 

// 我们将在变量  <code>t</code> and vector <code>U</code>  中存储当前的时间和状态。

    double       t            = 0.; 
    unsigned int output_cycle = 0; 

    print_head(pcout, "interpolate initial values"); 
    vector_type U = interpolate_initial_values(); 
// @sect5{Resume}  

// 默认情况下，布尔值 <code>resume</code> 被设置为false，也就是说，下面的代码段不会被运行。然而，如果 <code>resume==true</code> ，我们表明我们确实有一个中断的计算，程序应重新启动，从检查点文件中读入由 <code>t</code> 、 <code>output_cycle</code>, and <code>U</code> 组成的旧状态。这些检查点文件将在下面讨论的 <code>output()</code> 程序中创建。

    if (resume) 
      { 
        print_head(pcout, "restore interrupted computation"); 

        const unsigned int i = 
          discretization.triangulation.locally_owned_subdomain(); 

        const std::string name = base_name + "-checkpoint-" + 
                                 Utilities::int_to_string(i, 4) + ".archive"; 
        std::ifstream file(name, std::ios::binary); 

// 我们使用一个 <code>boost::archive</code> 来存储和读入检查点状态的内容。

        boost::archive::binary_iarchive ia(file); 
        ia >> t >> output_cycle; 

        for (auto &it1 : U) 
          { 
// <code>it1</code>  遍历状态向量的所有组件  <code>U</code>  。我们依次读入分量的每一个条目，之后更新ghost层。

            for (auto &it2 : it1) 
              ia >> it2; 
            it1.update_ghost_values(); 
          } 
      } 

// 随着初始状态的建立，或中断状态的恢复，是时候进入主循环了。

    output(U, base_name, t, output_cycle++); 

    print_head(pcout, "enter main loop"); 

    for (unsigned int cycle = 1; t < t_final; ++cycle) 
      { 

// 我们首先打印一个信息性的状态信息

        std::ostringstream head; 
        std::ostringstream secondary; 

        head << "Cycle  " << Utilities::int_to_string(cycle, 6) << "  (" // 
             << std::fixed << std::setprecision(1) << t / t_final * 100  // 
             << "%)"; 
        secondary << "at time t = " << std::setprecision(8) << std::fixed << t; 

        print_head(pcout, head.str(), secondary.str()); 

// 然后执行一个单一的前向欧拉步骤。请注意，状态向量 <code>U</code> 被就地更新， <code>time_stepping.make_one_step()</code> 返回选择的步长。

        t += time_stepping.make_one_step(U, t); 

// 后期处理、生成输出和写出当前状态是一个CPU和IO密集型的任务，我们不能在每个时间步长进行处理

// -- 特别是在显式时间步进中。因此，我们只在超过 <code>output_granularity</code> 设定的阈值时，通过调用 <code>output()</code> 函数安排输出。

        if (t > output_cycle * output_granularity) 
          { 
            output(U, base_name, t, output_cycle, true); 
            ++output_cycle; 
          } 
      } 

// 我们等待任何剩余的后台输出线程完成，然后打印一个摘要并退出。

    if (background_thread_state.valid()) 
      background_thread_state.wait(); 

    computing_timer.print_summary(); 
    pcout << timer_output.str() << std::endl; 
  } 

//  <code>interpolate_initial_values</code> 将初始时间 "t "作为输入参数，并在 <code>InitialValues<dim>::initial_state</code> 对象的帮助下填充状态向量 <code>U</code> 。

  template <int dim> 
  typename MainLoop<dim>::vector_type 
  MainLoop<dim>::interpolate_initial_values(const double t) 
  { 
    pcout << "MainLoop<dim>::interpolate_initial_values(t = " << t << ")" 
          << std::endl; 
    TimerOutput::Scope scope(computing_timer, 
                             "main_loop - setup scratch space"); 

    vector_type U; 

    for (auto &it : U) 
      it.reinit(offline_data.partitioner); 

    constexpr auto problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

//  <code>InitialValues<dim>::initial_state</code> 的函数签名对于 VectorTools::interpolate(). 来说不太合适。我们通过以下方式来解决这个问题：首先，创建一个lambda函数，对于给定的位置 <code>x</code> 只返回 <code>i</code> 的第三部分的值。在ScalarFunctionFromFunctionObject包装器的帮助下，这个lambda又被转换为一个 dealii::Function 。

    for (unsigned int i = 0; i < problem_dimension; ++i) 
      VectorTools::interpolate(offline_data.dof_handler, 
                               ScalarFunctionFromFunctionObject<dim, double>( 
                                 [&](const Point<dim> &x) { 
                                   return initial_values.initial_state(x, t)[i]; 
                                 }), 
                               U[i]); 

    for (auto &it : U) 
      it.update_ghost_values(); 

    return U; 
  } 
// @sect5{Output and checkpointing}  

// 写出最终的 vtk 文件是一项相当密集的 IO 任务，会让主循环停滞一段时间。为了避免这种情况，我们使用了<a
//  href="https:en.wikipedia.org/wiki/Asynchronous_I/O">asynchronous
//  IO</a>的策略，即创建一个后台线程，在主循环被允许继续的情况下执行IO。为了使其发挥作用，我们必须注意两件事。

// - 在运行  <code>output_worker</code>  线程之前，我们必须创建一个状态向量  <code>U</code>  的副本。我们把它存储在向量  <code>output_vector</code>  中。

// - 我们必须避免在后台线程中进行任何MPI通信，否则程序可能会出现死锁。这意味着我们必须在工作线程之外运行后处理程序。

  template <int dim> 
  void MainLoop<dim>::output(const typename MainLoop<dim>::vector_type &U, 
                             const std::string &                        name, 
                             const double                               t, 
                             const unsigned int                         cycle, 
                             const bool checkpoint) 
  { 
    pcout << "MainLoop<dim>::output(t = " << t 
          << ", checkpoint = " << checkpoint << ")" << std::endl; 

// 如果设置了异步回写选项，我们会启动一个后台线程，执行所有的慢速IO到磁盘。在这种情况下，我们必须确保后台线程确实完成了运行。如果没有，我们必须等待它完成。我们用<a
//  href="https:en.cppreference.com/w/cpp/thread/async"><code>std::async()</code></a>启动上述背景线程，该线程返回<a
//  href="https:en.cppreference.com/w/cpp/thread/future"><code>std::future</code></a>对象。这个 <code>std::future</code> 对象包含了函数的返回值，在我们的例子中就是 <code>void</code>  。

    if (background_thread_state.valid()) 
      { 
        TimerOutput::Scope timer(computing_timer, "main_loop - stalled output"); 
        background_thread_state.wait(); 
      } 

    constexpr auto problem_dimension = 
      ProblemDescription<dim>::problem_dimension; 

// 在这一点上，我们制作一份状态向量的副本，运行schlieren后处理器，并运行 DataOut<dim>::build_patches()  实际输出代码是标准的。我们创建一个DataOut实例，附加所有我们想要输出的数据向量，并调用 DataOut<dim>::build_patches(). ，但是有一个转折。为了在后台线程上执行异步IO，我们将DataOut<dim>对象创建为一个共享指针，传递给工作线程，以确保一旦我们退出这个函数，工作线程完成后，DataOut<dim>对象再次被销毁。

    for (unsigned int i = 0; i < problem_dimension; ++i) 
      { 
        output_vector[i] = U[i]; 
        output_vector[i].update_ghost_values(); 
      } 

    schlieren_postprocessor.compute_schlieren(output_vector); 

    auto data_out = std::make_shared<DataOut<dim>>(); 

    data_out->attach_dof_handler(offline_data.dof_handler); 

    const auto &component_names = ProblemDescription<dim>::component_names; 

    for (unsigned int i = 0; i < problem_dimension; ++i) 
      data_out->add_data_vector(output_vector[i], component_names[i]); 

    data_out->add_data_vector(schlieren_postprocessor.schlieren, 
                              "schlieren_plot"); 

    data_out->build_patches(discretization.mapping, 
                            discretization.finite_element.degree - 1); 

// 接下来我们为后台线程创建一个lambda函数。我们 <a href="https:en.cppreference.com/w/cpp/language/lambda">capture</a>  <code>this</code>  指针以及输出函数的大部分参数的值，这样我们就可以在lambda函数中访问它们。

    const auto output_worker = [this, name, t, cycle, checkpoint, data_out]() { 
      if (checkpoint) 
        { 

// 我们通过对<a href="Resume">resume logic</a>的精确反向操作来检查当前状态。

          const unsigned int i = 
            discretization.triangulation.locally_owned_subdomain(); 
          std::string filename = 
            name + "-checkpoint-" + Utilities::int_to_string(i, 4) + ".archive"; 

          std::ofstream file(filename, std::ios::binary | std::ios::trunc); 

          boost::archive::binary_oarchive oa(file); 
          oa << t << cycle; 
          for (const auto &it1 : output_vector) 
            for (const auto &it2 : it1) 
              oa << it2; 
        } 

      DataOutBase::VtkFlags flags(t, 
                                  cycle, 
                                  true, 
                                  DataOutBase::VtkFlags::best_speed); 
      data_out->set_flags(flags); 

      data_out->write_vtu_with_pvtu_record( 
        "", name + "-solution", cycle, mpi_communicator, 6); 
    }; 

// 如果设置了异步回写选项，我们在<a
//  href="https:en.cppreference.com/w/cpp/thread/async"><code>std::async</code></a>函数的帮助下启动一个新的后台线程。该函数返回一个<a
//  href="https:en.cppreference.com/w/cpp/thread/future"><code>std::future</code></a>对象，我们可以用它来查询后台线程的状态。在这一点上，我们可以从 <code>output()</code> 函数中返回，继续在主循环中进行时间步进

// - 该线程将在后台运行。

    if (asynchronous_writeback) 
      { 
        background_thread_state = std::async(std::launch::async, output_worker); 
      } 
    else 
      { 
        output_worker(); 
      } 
  } 

} // namespace Step69 

// 最后是主函数。

int main(int argc, char *argv[]) 
{ 
  try 
    { 
      constexpr int dim = 2; 

      using namespace dealii; 
      using namespace Step69; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv); 

      MPI_Comm      mpi_communicator(MPI_COMM_WORLD); 
      MainLoop<dim> main_loop(mpi_communicator); 

      main_loop.run(); 
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
} 


