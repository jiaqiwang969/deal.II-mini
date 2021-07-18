

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2000 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004, 2005, 
 * Timo Heister, 2013 
 */ 



// 首先是通常的头文件列表，这些文件已经在以前的示例程序中使用过了。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/multithread_info.h> 
#include <deal.II/base/conditional_ostream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/petsc_vector.h> 
#include <deal.II/lac/petsc_sparse_matrix.h> 
#include <deal.II/lac/petsc_solver.h> 
#include <deal.II/lac/petsc_precondition.h> 
#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparsity_tools.h> 
#include <deal.II/distributed/shared_tria.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 这里是头文件中仅有的三个新东西：一个包含文件，其中实现了等级为2和4的对称张量，正如介绍中所介绍的那样。

#include <deal.II/base/symmetric_tensor.h> 

// 最后是一个包含一些函数的头文件，这些函数将帮助我们计算域中特定点的局部坐标系的旋转矩阵。

#include <deal.II/physics/transformations.h> 

// 然后，这又是简单的C++。

#include <fstream> 
#include <iostream> 
#include <iomanip> 

// 最后一步和以前所有的程序一样。

namespace Step18 
{ 
  using namespace dealii; 
// @sect3{The <code>PointHistory</code> class}  

// 正如介绍中提到的，我们必须在正交点存储旧的应力，这样我们就可以在下一个时间步骤中计算这一点的残余力。仅仅这一点还不能保证只有一个成员的结构，但在更复杂的应用中，我们还必须在正交点上存储更多的信息，比如塑性的历史变量等。从本质上讲，我们必须在这里存储所有影响材料当前状态的信息，在塑性中，这些信息是由变形历史变量决定的。

// 除了能够存储数据之外，我们不会给这个类任何有意义的功能，也就是说，没有构造函数、析构函数或其他成员函数。在这种 "哑巴 "类的情况下，我们通常选择将其声明为  <code>struct</code> rather than <code>class</code>  ，以表明它们更接近于C语言风格的结构而不是C++风格的类。

  template <int dim> 
  struct PointHistory 
  { 
    SymmetricTensor<2, dim> old_stress; 
  }; 
// @sect3{The stress-strain tensor}  

// 接下来，我们定义弹性中的应力和应变的线性关系。它由一个等级为4的张量给出，通常被写成  $C_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}) + \lambda \delta_{ij} \delta_{kl}$  的形式。这个张量将等级2的对称张量映射到等级2的对称张量。对于Lam&eacute;常数 $\lambda$ 和 $\mu$ 的给定值，一个实现其创建的函数是直接的。

  template <int dim> 
  SymmetricTensor<4, dim> get_stress_strain_tensor(const double lambda, 
                                                   const double mu) 
  { 
    SymmetricTensor<4, dim> tmp; 
    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = 0; j < dim; ++j) 
        for (unsigned int k = 0; k < dim; ++k) 
          for (unsigned int l = 0; l < dim; ++l) 
            tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) + 
                               ((i == l) && (j == k) ? mu : 0.0) + 
                               ((i == j) && (k == l) ? lambda : 0.0)); 
    return tmp; 
  } 

// 通过这个函数，我们将在下面的主类中定义一个静态成员变量，在整个程序中作为应力-应变张量使用。请注意，在更复杂的程序中，这可能是某个类的成员变量，或者是一个根据其他输入返回应力-应变关系的函数。例如，在损伤理论模型中，Lam&eacute;常数被认为是一个点的先前应力/应变历史的函数。相反，在塑性中，如果材料在某一点达到了屈服应力，那么应力-应变张量的形式就会被修改，而且可能还取决于其先前的历史。

// 然而，在本程序中，我们假设材料是完全弹性和线性的，恒定的应力-应变张量对我们目前的目的来说是足够的。

//  @sect3{Auxiliary functions}  

// 在程序的其他部分之前，这里有几个我们需要的函数作为工具。这些是在内循环中调用的小函数，所以我们把它们标记为  <code>inline</code>  。

// 第一个是通过形成这个形状函数的对称梯度来计算形状函数 <code>shape_func</code> at quadrature point <code>q_point</code> 的对称应变张量。当我们想形成矩阵时，我们需要这样做，比如说。

// 我们应该注意到，在以前处理矢量值问题的例子中，我们总是问有限元对象在哪个矢量分量中的形状函数实际上是不为零的，从而避免计算任何我们反正可以证明为零的项。为此，我们使用了 <code>fe.system_to_component_index</code> 函数来返回形状函数在哪个分量中为零，同时 <code>fe_values.shape_value</code> 和 <code>fe_values.shape_grad</code> 函数只返回形状函数的单个非零分量的值和梯度，如果这是一个矢量值元素。

// 这是一个优化，如果不是非常关键的时间，我们可以用一个更简单的技术来解决：只需向 <code>fe_values</code> 询问一个给定形状函数的给定分量在给定正交点的值或梯度。这就是  <code>fe_values.shape_grad_component(shape_func,q_point,i)</code>  调用的作用：返回形状函数  <code>shape_func</code>  的第  <code>q_point</code>  个分量在正交点的全部梯度。如果某个形状函数的某个分量总是为零，那么这将简单地总是返回零。

// 如前所述，使用 <code>fe_values.shape_grad_component</code> 而不是 <code>fe.system_to_component_index</code> 和 <code>fe_values.shape_grad</code> 的组合可能效率较低，但其实现已针对这种情况进行了优化，应该不会有很大的减慢。我们在这里演示这个技术，因为它是如此的简单和直接。

  template <int dim> 
  inline SymmetricTensor<2, dim> get_strain(const FEValues<dim> &fe_values, 
                                            const unsigned int   shape_func, 
                                            const unsigned int   q_point) 
  { 

// 声明一个将保存返回值的暂存器。

    SymmetricTensor<2, dim> tmp; 

// 首先，填充对角线项，这只是矢量值形状函数的方向 <code>i</code> of the <code>i</code> 分量的导数。

    for (unsigned int i = 0; i < dim; ++i) 
      tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i]; 

// 然后填充应变张量的其余部分。注意，由于张量是对称的，我们只需要计算一半（这里：右上角）的非对角线元素， <code>SymmetricTensor</code> 类的实现确保至少到外面的对称条目也被填充（实际上，这个类当然只存储一份）。在这里，我们选择了张量的右上半部分，但是左下半部分也一样好。

    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = i + 1; j < dim; ++j) 
        tmp[i][j] = 
          (fe_values.shape_grad_component(shape_func, q_point, i)[j] + 
           fe_values.shape_grad_component(shape_func, q_point, j)[i]) / 
          2; 

    return tmp; 
  } 

// 第二个函数做了非常类似的事情（因此被赋予相同的名字）：从一个矢量值场的梯度计算对称应变张量。如果你已经有了一个解场， <code>fe_values.get_function_gradients</code> 函数允许你在一个正交点上提取解场的每个分量的梯度。它返回的是一个秩-1张量的矢量：解的每个矢量分量有一个秩-1张量（梯度）。由此，我们必须通过转换数据存储格式和对称化来重建（对称的）应变张量。我们用和上面一样的方法来做，也就是说，我们通过首先填充对角线，然后只填充对称张量的一半来避免一些计算（ <code>SymmetricTensor</code> 类确保只写两个对称分量中的一个就足够了）。

// 不过在我们这样做之前，我们要确保输入有我们期望的那种结构：即有 <code>dim</code> 个矢量分量，即每个坐标方向有一个位移分量。我们用 <code>Assert</code> 宏来测试这一点，如果不符合条件，我们的程序就会被终止。

  template <int dim> 
  inline SymmetricTensor<2, dim> 
  get_strain(const std::vector<Tensor<1, dim>> &grad) 
  { 
    Assert(grad.size() == dim, ExcInternalError()); 

    SymmetricTensor<2, dim> strain; 
    for (unsigned int i = 0; i < dim; ++i) 
      strain[i][i] = grad[i][i]; 

    for (unsigned int i = 0; i < dim; ++i) 
      for (unsigned int j = i + 1; j < dim; ++j) 
        strain[i][j] = (grad[i][j] + grad[j][i]) / 2; 

    return strain; 
  } 

// 最后，下面我们将需要一个函数来计算某一点的位移所引起的旋转矩阵。当然，事实上，单点的位移只有一个方向和一个幅度，诱发旋转的是方向和幅度的变化。实际上，旋转矩阵可以通过位移的梯度来计算，或者更具体地说，通过卷曲来计算。

// 确定旋转矩阵的公式有点笨拙，特别是在三维中。对于2D来说，有一个更简单的方法，所以我们把这个函数实现了两次，一次用于2D，一次用于3D，这样我们就可以在两个空间维度上编译和使用这个程序，如果需要的话--毕竟，deal.II是关于独立维度编程和重复使用算法的，在2D的廉价计算中经过测试，在3D的更昂贵的计算中使用。下面是一种情况，我们必须为2D和3D实现不同的算法，但可以用独立于空间维度的方式来编写程序的其余部分。

// 所以，不用再多说了，来看看2D的实现。

  Tensor<2, 2> get_rotation_matrix(const std::vector<Tensor<1, 2>> &grad_u) 
  { 

// 首先，根据梯度计算出速度场的卷曲。注意，我们是在2d中，所以旋转是一个标量。

    const double curl = (grad_u[1][0] - grad_u[0][1]); 

// 由此计算出旋转的角度。

    const double angle = std::atan(curl); 

// 由此，建立反对称的旋转矩阵。我们希望这个旋转矩阵能够代表本地坐标系相对于全局直角坐标系的旋转，所以我们用一个负的角度来构建它。因此，这个旋转矩阵代表了从本地坐标系移动到全局坐标系所需的旋转。

    return Physics::Transformations::Rotations::rotation_matrix_2d(-angle); 
  } 

// 三维的情况就比较复杂了。

  Tensor<2, 3> get_rotation_matrix(const std::vector<Tensor<1, 3>> &grad_u) 
  { 

// 同样首先计算速度场的卷曲。这一次，它是一个实数向量。

    const Point<3> curl(grad_u[2][1] - grad_u[1][2], 
                        grad_u[0][2] - grad_u[2][0], 
                        grad_u[1][0] - grad_u[0][1]); 

// 从这个矢量中，利用它的大小，计算出旋转角度的正切值，并由此计算出相对于直角坐标系的实际旋转角度。

    const double tan_angle = std::sqrt(curl * curl); 
    const double angle     = std::atan(tan_angle); 

// 现在，这里有一个问题：如果旋转角度太小，那就意味着没有旋转发生（例如平移运动）。在这种情况下，旋转矩阵就是身份矩阵。

// 我们强调这一点的原因是，在这种情况下，我们有  <code>tan_angle==0</code>  。再往下看，我们在计算旋转轴的时候需要除以这个数字，这样做除法的时候会遇到麻烦。因此，让我们走捷径，如果旋转角度真的很小，就简单地返回同一矩阵。

    if (std::abs(angle) < 1e-9) 
      { 
        static const double rotation[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; 
        static const Tensor<2, 3> rot(rotation); 
        return rot; 
      } 

// 否则计算真实的旋转矩阵。为此，我们再次依靠一个预定义的函数来计算本地坐标系的旋转矩阵。

    const Point<3> axis = curl / tan_angle; 
    return Physics::Transformations::Rotations::rotation_matrix_3d(axis, 
                                                                   -angle); 
  } 

//  @sect3{The <code>TopLevel</code> class}  

// 这就是程序的主类。由于命名空间已经表明了我们要解决的问题，让我们用它的作用来称呼它：它引导着程序的流程，也就是说，它是顶层驱动。

// 这个类的成员变量基本上和以前一样，即它必须有一个三角形，一个DoF处理程序和相关的对象，如约束条件，描述线性系统的变量等。现在还有很多成员函数，我们将在下面解释。

// 然而，该类的外部接口是不变的：它有一个公共的构造函数和析构函数，并且它有一个 <code>run</code> 函数来启动所有的工作。

  template <int dim> 
  class TopLevel 
  { 
  public: 
    TopLevel(); 
    ~TopLevel(); 
    void run(); 

  private: 

// 私有接口比  step-17  中的更加广泛。首先，我们显然需要创建初始网格的函数，设置描述当前网格上的线性系统的变量（即矩阵和向量），然后是实际组装系统的函数，指导每个时间步长中必须解决的问题，一个解决每个时间步长中出现的线性系统的函数（并返回它的迭代次数），最后在正确的网格上输出解向量。

    void create_coarse_grid(); 

    void setup_system(); 

    void assemble_system(); 

    void solve_timestep(); 

    unsigned int solve_linear_problem(); 

    void output_results() const; 

// 除了前两个，所有这些函数都在每个时间步中被调用。由于第一个时间步骤有点特殊，我们有单独的函数来描述一个时间步骤中必须发生的事情：一个用于第一个时间步骤，一个用于所有后续时间步骤。

    void do_initial_timestep(); 

    void do_timestep(); 

// 然后我们需要一大堆函数来做各种事情。第一个是细化初始网格：我们从原始状态的粗网格开始，解决这个问题，然后看一下，并相应地细化网格，然后重新开始同样的过程，再次以原始状态。因此，细化初始网格比在两个连续的时间步骤之间细化网格要简单一些，因为它不涉及将数据从旧的三角测量转移到新的三角测量，特别是存储在每个正交点的历史数据。

    void refine_initial_grid(); 

// 在每个时间步骤结束时，我们要根据这个时间步骤计算的增量位移来移动网格顶点。这就是完成这个任务的函数。

    void move_mesh(); 

// 接下来是两个处理存储在每个正交点的历史变量的函数。第一个函数在第一个时间步长之前被调用，为历史变量设置一个原始状态。它只对属于当前处理器的单元上的正交点起作用。

    void setup_quadrature_point_history(); 

// 第二项是在每个时间段结束时更新历史变量。

    void update_quadrature_point_history(); 

// 这是新的共享三角法。

    parallel::shared::Triangulation<dim> triangulation; 

    FESystem<dim> fe; 

    DoFHandler<dim> dof_handler; 

    AffineConstraints<double> hanging_node_constraints; 

// 这个程序的一个不同之处在于，我们在类声明中声明了正交公式。原因是在所有其他程序中，如果我们在计算矩阵和右手边时使用不同的正交公式，并没有什么坏处，比如说。然而，在目前的情况下，它确实如此：我们在正交点中存储了信息，所以我们必须确保程序的所有部分都同意它们的位置以及每个单元格上有多少个。因此，让我们首先声明将在整个程序中使用的正交公式...。

    const QGauss<dim> quadrature_formula; 

// ......然后也有一个历史对象的向量，在我们负责的那些单元格上的每个正交点都有一个（也就是说，我们不为其他处理器拥有的单元格上的正交点存储历史数据）。请注意，我们可以像在  step-44  中那样使用 CellDataStorage 类来代替我们自己存储和管理这些数据。然而，为了演示的目的，在这种情况下，我们手动管理存储。

    std::vector<PointHistory<dim>> quadrature_point_history; 

// 这个对象的访问方式是通过每个单元格、面或边持有的 <code>user pointer</code> ：它是一个 <code>void*</code> 指针，可以被应用程序用来将任意的数据与单元格、面或边联系起来。程序对这些数据的实际操作属于自己的职责范围，库只是为这些指针分配了一些空间，而应用程序可以设置和读取这些对象中的每个指针。

// 进一步说：我们需要待解的线性系统的对象，即矩阵、右手边的向量和解向量。由于我们预计要解决大问题，我们使用了与 step-17 中相同的类型，即建立在PETSc库之上的分布式%并行矩阵和向量。方便的是，它们也可以在只在一台机器上运行时使用，在这种情况下，这台机器正好是我们的%并行宇宙中唯一的机器。

// 然而，与 step-17 不同的是，我们不以分布式方式存储解向量--这里是在每个时间步骤中计算的增量位移。也就是说，在计算时它当然必须是一个分布式矢量，但紧接着我们确保每个处理器都有一个完整的副本。原因是我们已经在 step-17 中看到，许多函数需要一个完整的副本。虽然得到它并不难，但这需要在网络上进行通信，因此很慢。此外，这些都是重复的相同操作，这当然是不可取的，除非不必总是存储整个向量的收益超过了它。在编写这个程序时，事实证明，我们在很多地方都需要一份完整的解决方案，以至于只在必要时才获得它似乎不值得。相反，我们选择一劳永逸地获得完整的副本，而立即摆脱分散的副本。因此，请注意， <code>incremental_displacement</code> 的声明并没有像中间命名空间 <code>MPI</code> 所表示的那样，表示一个分布式向量。

    PETScWrappers::MPI::SparseMatrix system_matrix; 

    PETScWrappers::MPI::Vector system_rhs; 

    Vector<double> incremental_displacement; 

// 接下来的变量块与问题的时间依赖性有关：它们表示我们要模拟的时间间隔的长度，现在的时间和时间步数，以及现在时间步数的长度。

    double       present_time; 
    double       present_timestep; 
    double       end_time; 
    unsigned int timestep_no; 

// 然后是几个与%并行处理有关的变量：首先，一个变量表示我们使用的MPI通信器，然后是两个数字，告诉我们有多少个参与的处理器，以及我们在这个世界上的位置。最后，一个流对象，确保只有一个处理器实际产生输出到控制台。这与  step-17  中的所有内容相同。

    MPI_Comm mpi_communicator; 

    const unsigned int n_mpi_processes; 

    const unsigned int this_mpi_process; 

    ConditionalOStream pcout; 

// 我们正在存储本地拥有的和本地相关的索引。

    IndexSet locally_owned_dofs; 
    IndexSet locally_relevant_dofs; 

// 最后，我们有一个静态变量，表示应力和应变之间的线性关系。由于它是一个不依赖任何输入的常量对象（至少在这个程序中不依赖），我们把它作为一个静态变量，并将在我们定义这个类的构造函数的同一个地方初始化它。

    static const SymmetricTensor<4, dim> stress_strain_tensor; 
  }; 
// @sect3{The <code>BodyForce</code> class}  

// 在我们进入这个程序的主要功能之前，我们必须定义哪些力将作用在我们想要研究的变形的体上。这些力可以是体力，也可以是边界力。体力通常是由四种基本的物理力类型之一所介导的：重力、强弱相互作用和电磁力。除非人们想考虑亚原子物体（对于这些物体，无论如何准静态变形是不相关的，也是不合适的描述），否则只需要考虑引力和电磁力。为了简单起见，让我们假设我们的身体有一定的质量密度，但要么是非磁性的，不导电的，要么周围没有明显的电磁场。在这种情况下，身体的力只是 <code>rho g</code>, where <code>rho</code> 是材料密度， <code>g</code> 是一个负Z方向的矢量，大小为9.81米/秒^2。 密度和 <code>g</code> 都是在函数中定义的，我们把7700 kg/m^3作为密度，这是对钢材通常假定的值。

// 为了更普遍一点，也为了能够在2d中进行计算，我们意识到体力总是一个返回 <code>dim</code> 维矢量的函数。我们假设重力沿着最后一个，即 <code>dim-1</code> 个坐标的负方向作用。考虑到以前的例子程序中的类似定义，这个函数的其余实现应该大部分是不言自明的。请注意，身体的力量与位置无关；为了避免编译器对未使用的函数参数发出警告，我们因此注释了 <code>vector_value</code> 函数的第一个参数的名称。

  template <int dim> 
  class BodyForce : public Function<dim> 
  { 
  public: 
    BodyForce(); 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override; 

    virtual void 
    vector_value_list(const std::vector<Point<dim>> &points, 
                      std::vector<Vector<double>> &  value_list) const override; 
  }; 

  template <int dim> 
  BodyForce<dim>::BodyForce() 
    : Function<dim>(dim) 
  {} 

  template <int dim> 
  inline void BodyForce<dim>::vector_value(const Point<dim> & /*p*/, 
                                           Vector<double> &values) const 
  { 
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 

    const double g   = 9.81; 
    const double rho = 7700; 

    values          = 0; 
    values(dim - 1) = -rho * g; 
  } 

  template <int dim> 
  void BodyForce<dim>::vector_value_list( 
    const std::vector<Point<dim>> &points, 
    std::vector<Vector<double>> &  value_list) const 
  { 
    const unsigned int n_points = points.size(); 

    Assert(value_list.size() == n_points, 
           ExcDimensionMismatch(value_list.size(), n_points)); 

    for (unsigned int p = 0; p < n_points; ++p) 
      BodyForce<dim>::vector_value(points[p], value_list[p]); 
  } 

//  @sect3{The <code>IncrementalBoundaryValue</code> class}  

// 除了身体的力之外，运动还可以由边界力和强制边界位移引起。后一种情况相当于以这样的方式选择力，使其诱发某种位移。

// 对于准静态位移，典型的边界力是对一个体的压力，或者对另一个体的切向摩擦。我们在这里选择了一种更简单的情况：我们规定了边界（部分）的某种运动，或者至少是位移矢量的某些分量。我们用另一个矢量值函数来描述，对于边界上的某一点，返回规定的位移。

// 由于我们有一个随时间变化的问题，边界的位移增量等于在时间段内累积的位移。因此，该类必须同时知道当前时间和当前时间步长，然后可以将位移增量近似为当前速度乘以当前时间步长。

// 在本程序中，我们选择了一种简单的边界位移形式：我们以恒定的速度向下位移顶部的边界。边界的其余部分要么是固定的（然后用一个 <code>Functions::ZeroFunction</code> 类型的对象来描述），要么是自由的（Neumann类型，在这种情况下不需要做任何特殊的事情）。 利用我们在前面所有的例子程序中获得的知识，描述持续向下运动的类的实现应该是很明显的。

  template <int dim> 
  class IncrementalBoundaryValues : public Function<dim> 
  { 
  public: 
    IncrementalBoundaryValues(const double present_time, 
                              const double present_timestep); 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  values) const override; 

    virtual void 
    vector_value_list(const std::vector<Point<dim>> &points, 
                      std::vector<Vector<double>> &  value_list) const override; 

  private: 
    const double velocity; 
    const double present_time; 
    const double present_timestep; 
  }; 

  template <int dim> 
  IncrementalBoundaryValues<dim>::IncrementalBoundaryValues( 
    const double present_time, 
    const double present_timestep) 
    : Function<dim>(dim) 
    , velocity(.08) 
    , present_time(present_time) 
    , present_timestep(present_timestep) 
  {} 

  template <int dim> 
  void 
  IncrementalBoundaryValues<dim>::vector_value(const Point<dim> & /*p*/, 
                                               Vector<double> &values) const 
  { 
    Assert(values.size() == dim, ExcDimensionMismatch(values.size(), dim)); 

    values    = 0; 
    values(2) = -present_timestep * velocity; 
  } 

  template <int dim> 
  void IncrementalBoundaryValues<dim>::vector_value_list( 
    const std::vector<Point<dim>> &points, 
    std::vector<Vector<double>> &  value_list) const 
  { 
    const unsigned int n_points = points.size(); 

    Assert(value_list.size() == n_points, 
           ExcDimensionMismatch(value_list.size(), n_points)); 

    for (unsigned int p = 0; p < n_points; ++p) 
      IncrementalBoundaryValues<dim>::vector_value(points[p], value_list[p]); 
  } 

//  @sect3{Implementation of the <code>TopLevel</code> class}  

// 现在是主类的实现。首先，我们初始化应力应变张量，我们将其声明为一个静态常量变量。我们选择了适合于钢铁的Lam&eacute;常数。

  template <int dim> 
  const SymmetricTensor<4, dim> TopLevel<dim>::stress_strain_tensor = 
    get_stress_strain_tensor<dim>(
      /*lambda = */ 9.695e10, 
      /*mu =  */  7.617e10)

//  @sect4{The public interface}  

// 下一步是构造函数和析构函数的定义。这里没有什么惊喜：我们为解的每个 <code>dim</code> 矢量分量选择线性和连续的有限元，以及每个坐标方向上有2个点的高斯正交公式。解构器应该是显而易见的。

  template <int dim> 
  TopLevel<dim>::TopLevel() 
    : triangulation(MPI_COMM_WORLD) 
    , fe(FE_Q<dim>(1), dim) 
    , dof_handler(triangulation) 
    , quadrature_formula(fe.degree + 1) 
    , present_time(0.0) 
    , present_timestep(1.0) 
    , end_time(10.0) 
    , timestep_no(0) 
    , mpi_communicator(MPI_COMM_WORLD) 
    , n_mpi_processes(Utilities::MPI::n_mpi_processes(mpi_communicator)) 
    , this_mpi_process(Utilities::MPI::this_mpi_process(mpi_communicator)) 
    , pcout(std::cout, this_mpi_process == 0) 
  {} 

  template <int dim> 
  TopLevel<dim>::~TopLevel() 
  { 
    dof_handler.clear(); 
  } 

// 最后一个公共函数是指导所有工作的函数，  <code>run()</code>  。它初始化了描述我们目前所处时间位置的变量，然后运行第一个时间步骤，再循环所有其他时间步骤。请注意，为了简单起见，我们使用一个固定的时间步长，而一个更复杂的程序当然要以某种更合理的方式自适应地选择它。

  template <int dim> 
  void TopLevel<dim>::run() 
  { 
    do_initial_timestep(); 

    while (present_time < end_time) 
      do_timestep(); 
  } 
// @sect4{TopLevel::create_coarse_grid}  

// 按照上面声明的顺序，下一个函数是创建粗略网格的函数，我们从这里开始。在这个示例程序中，我们想计算一个圆柱体在轴向压缩下的变形。因此第一步是生成一个长度为3，内外半径分别为0.8和1的圆柱体的网格。幸运的是，有一个库函数可以生成这样的网格。

// 在第二步中，我们必须在圆柱体的上表面和下表面关联边界条件。我们为边界面选择一个边界指示器0，这些边界面的中点的Z坐标为0（底面），Z=3的指示器为1（顶面）；最后，我们对圆柱体外壳内部的所有面使用边界指示器2，外部使用3。

  template <int dim> 
  void TopLevel<dim>::create_coarse_grid() 
  { 
    const double inner_radius = 0.8, outer_radius = 1; 
    GridGenerator::cylinder_shell(triangulation, 3, inner_radius, outer_radius); 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->at_boundary()) 
          { 
            const Point<dim> face_center = face->center(); 

            if (face_center[2] == 0) 
              face->set_boundary_id(0); 
            else if (face_center[2] == 3) 
              face->set_boundary_id(1); 
            else if (std::sqrt(face_center[0] * face_center[0] + 
                               face_center[1] * face_center[1]) < 
                     (inner_radius + outer_radius) / 2) 
              face->set_boundary_id(2); 
            else 
              face->set_boundary_id(3); 
          } 

// 一旦完成了这些，我们就可以对网格进行一次全面的细化。

    triangulation.refine_global(1); 

// 作为最后一步，我们需要设置一个干净的数据状态，我们将这些数据存储在目前处理器上处理的所有单元的正交点中。

    setup_quadrature_point_history(); 
  } 

//  @sect4{TopLevel::setup_system}  

// 下一个函数是为一个给定的网格设置数据结构。这与 step-17 中的方法基本相同：分配自由度，然后对这些自由度进行排序，使每个处理器得到一个连续的块。请注意，每个处理器的细分块是在创建或完善网格的函数中处理的，与之前的例子程序不同（发生这种情况的时间点主要是口味问题；在这里，我们选择在创建网格时进行，因为在 <code>do_initial_timestep</code> 和 <code>do_timestep</code> 函数中，我们想在还没有调用当前函数的时候输出每个处理器上的单元数量）。

  template <int dim> 
  void TopLevel<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 
    locally_owned_dofs = dof_handler.locally_owned_dofs(); 
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs); 

// 下一步是设置由于悬挂节点而产生的约束。这在以前已经处理过很多次了。

    hanging_node_constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, 
                                            hanging_node_constraints); 
    hanging_node_constraints.close(); 

// 然后我们要设置矩阵。这里我们偏离了  step-17  ，在那里我们简单地使用了PETSc的能力，即只知道矩阵的大小，随后分配那些被写入的非零元素。虽然从正确性的角度来看，这样做很好，但是效率却不高：如果我们不给PETSc提供关于哪些元素被写入的线索，那么当我们第一次设置矩阵中的元素时（即在第一个时间步中），它的速度会慢得令人难以忍受。后来，当元素被分配后，一切都快多了。在我们所做的实验中，如果我们指示PETSc哪些元素将被使用，哪些不被使用，那么第一个时间步骤可以加快近两个数量级。

// 要做到这一点，我们首先要生成我们要处理的矩阵的稀疏模式，并确保浓缩的悬挂节点约束在稀疏模式中增加必要的额外条目。

    DynamicSparsityPattern sparsity_pattern(locally_relevant_dofs); 
    DoFTools::make_sparsity_pattern(dof_handler, 
                                    sparsity_pattern, 
                                    hanging_node_constraints, 
                                    /*保持约束性dofs  */ false)

    SparsityTools::distribute_sparsity_pattern(sparsity_pattern, 
                                               locally_owned_dofs, 
                                               mpi_communicator, 
                                               locally_relevant_dofs); 

// 注意，我们在这里使用了已经在 step-11 中介绍过的 <code>DynamicSparsityPattern</code> 类，而不是我们在所有其他情况下使用的 <code>SparsityPattern</code> 类。其原因是，为了使后一个类发挥作用，我们必须给每一行的条目数提供一个初始的上限，这项任务传统上是由 <code>DoFHandler::max_couplings_between_dofs()</code> 完成。然而，这个函数有一个严重的问题：它必须计算每一行中非零项的数量的上限，而这是一个相当复杂的任务，特别是在3D中。实际上，虽然它在2D中相当准确，但在3D中经常得出太大的数字，在这种情况下， <code>SparsityPattern</code> 一开始就分配了太多的内存，经常是几百MB。后来当 <code>DoFTools::make_sparsity_pattern</code> 被调用时，我们意识到我们不需要那么多的内存，但这时已经太晚了：对于大问题，临时分配太多的内存会导致内存不足的情况。

// 为了避免这种情况，我们采用了 <code>DynamicSparsityPattern</code> 类，该类速度较慢，但不需要预先估计每行非零条目的数量。因此，它在任何时候都只分配它所需要的内存，而且我们甚至可以为大型的三维问题建立它。

// 值得注意的是，由于 parallel::shared::Triangulation, 的特殊性，我们构建的稀疏模式是全局的，即包括所有的自由度，无论它们是属于我们所在的处理器还是另一个处理器（如果这个程序是通过MPI并行运行的）。这当然不是最好的--它限制了我们可以解决的问题的规模，因为在每个处理器上存储整个稀疏模式（即使只是短时间）的规模并不大。然而，在程序中还有几个地方我们是这样做的，例如，我们总是把全局三角测量和DoF处理对象保留在周围，即使我们只对它们的一部分进行工作。目前，deal.II没有必要的设施来完全分配这些对象（事实上，这项任务在自适应网格中很难实现，因为随着网格的自适应细化，领域的均衡分区往往会变得不均衡）。

// 有了这个数据结构，我们就可以进入PETSc稀疏矩阵，告诉它预先分配所有我们以后要写入的条目。

    system_matrix.reinit(locally_owned_dofs, 
                         locally_owned_dofs, 
                         sparsity_pattern, 
                         mpi_communicator); 

// 在这一点上，不再需要对稀疏模式有任何明确的了解，我们可以让 <code>sparsity_pattern</code> 这个变量离开范围，不会有任何问题。

// 这个函数的最后一个任务是将右侧向量和求解向量重置为正确的大小；记住，求解向量是一个本地向量，不像右侧向量是一个分布式的%并行向量，因此需要知道MPI通信器，它应该通过这个通信器来传输消息。

    system_rhs.reinit(locally_owned_dofs, mpi_communicator); 
    incremental_displacement.reinit(dof_handler.n_dofs()); 
  } 

//  @sect4{TopLevel::assemble_system}  

// 同样，组装系统矩阵和右手边的结构与之前许多例子程序中的结构相同。特别是，它主要等同于 step-17 ，除了不同的右手边，现在只需要考虑到内部应力。此外，通过使用 <code>SymmetricTensor</code> 类，组装矩阵明显变得更加透明：请注意形成2级和4级对称张量的标量积的优雅性。这个实现也更加通用，因为它与我们可能使用或不使用各向同性的弹性张量这一事实无关。

// 汇编程序的第一部分和以往一样。

  template <int dim> 
  void TopLevel<dim>::assemble_system() 
  { 
    system_rhs    = 0; 
    system_matrix = 0; 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients | 
                              update_quadrature_points | update_JxW_values); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
    const unsigned int n_q_points    = quadrature_formula.size(); 

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
    Vector<double>     cell_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    BodyForce<dim>              body_force; 
    std::vector<Vector<double>> body_force_values(n_q_points, 
                                                  Vector<double>(dim)); 

// 如同在  step-17  中一样，我们只需要在属于当前处理器的所有单元中进行循环。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell_matrix = 0; 
          cell_rhs    = 0; 

          fe_values.reinit(cell); 

// 然后在所有指数i,j和正交点上循环，并从这个单元中组合出系统矩阵的贡献。 注意我们如何从 <code>FEValues</code> 对象中提取给定正交点的形状函数的对称梯度（应变），以及我们如何优雅地形成三重收缩 <code>eps_phi_i : C : eps_phi_j</code> ；后者需要与 step-17 中需要的笨拙计算进行比较，无论是在介绍中还是在程序的相应位置。

          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < dofs_per_cell; ++j) 
              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                { 
                  const SymmetricTensor<2, dim> 
                    eps_phi_i = get_strain(fe_values, i, q_point), 
                    eps_phi_j = get_strain(fe_values, j, q_point); 

                  cell_matrix(i, j) += (eps_phi_i *            // 
                                        stress_strain_tensor * // 
                                        eps_phi_j              // 
                                        ) *                    // 
                                       fe_values.JxW(q_point); // 
                } 

// 然后也要组装本地的右手边贡献。为此，我们需要访问这个正交点的先验应力值。为了得到它，我们使用该单元的用户指针，该指针指向全局数组中与当前单元的第一个正交点相对应的正交点数据，然后添加一个与我们现在考虑的正交点的索引相对应的偏移量。

          const PointHistory<dim> *local_quadrature_points_data = 
            reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 

// 此外，我们还需要这个单元上的正交点的外体力值。

          body_force.vector_value_list(fe_values.get_quadrature_points(), 
                                       body_force_values); 

// 然后，我们可以循环计算这个单元上的所有自由度，并计算出对右侧的局部贡献。

          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const unsigned int component_i = 
                fe.system_to_component_index(i).first; 

              for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
                { 
                  const SymmetricTensor<2, dim> &old_stress = 
                    local_quadrature_points_data[q_point].old_stress; 

                  cell_rhs(i) += 
                    (body_force_values[q_point](component_i) * 
                       fe_values.shape_value(i, q_point) - 
                     old_stress * get_strain(fe_values, i, q_point)) * 
                    fe_values.JxW(q_point); 
                } 
            } 

// 现在我们有了对线性系统的局部贡献，我们需要将其转移到全局对象中。这与  step-17  中的做法完全相同。

          cell->get_dof_indices(local_dof_indices); 

          hanging_node_constraints.distribute_local_to_global(cell_matrix, 
                                                              cell_rhs, 
                                                              local_dof_indices, 
                                                              system_matrix, 
                                                              system_rhs); 
        } 

// 现在压缩矢量和系统矩阵。

    system_matrix.compress(VectorOperation::add); 
    system_rhs.compress(VectorOperation::add); 

// 最后一步是再次修复边界值，就像我们在以前的程序中已经做的那样。一个稍微复杂的问题是， <code>apply_boundary_values</code> 函数希望有一个与矩阵和右手边兼容的解向量（即这里是一个分布式的%并行向量，而不是我们在这个程序中使用的顺序向量），以便用正确的边界值预设解向量的条目。我们以临时向量的形式提供这样一个兼容向量，然后将其复制到顺序向量中。

// 我们通过展示边界值的灵活使用来弥补这种复杂性：按照我们创建三角形的方式，有三个不同的边界指标用来描述领域，分别对应于底面和顶面，以及内/外表面。我们希望施加以下类型的边界条件。内外圆柱体表面没有外力，这一事实对应于自然（诺伊曼型）边界条件，我们不需要做任何事情。在底部，我们希望完全没有运动，对应于圆柱体在边界的这一部分被夹住或粘住。然而，在顶部，我们希望有一个规定的垂直向下的运动来压缩圆柱体；此外，我们只希望限制垂直运动，而不是水平运动--可以把这种情况看作是一块油性良好的板坐在圆柱体的顶部将其向下推：圆柱体的原子被迫向下移动，但它们可以自由地沿着板水平滑动。

//描述这种情况的方法如下：对于边界指标为零（底面）的边界，我们使用一个二维的零函数，代表在任何坐标方向都没有运动。对于指标1（顶面）的边界，我们使用 <code>IncrementalBoundaryValues</code> 类，但我们为 <code>VectorTools::interpolate_boundary_values</code> 函数指定一个额外的参数，表示它应该适用于哪些矢量分量；这是一个针对每个矢量分量的bools矢量，由于我们只想限制垂直运动，它只有最后一个分量的设置。

    FEValuesExtractors::Scalar                z_component(dim - 1); 
    std::map<types::global_dof_index, double> boundary_values; 
    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             Functions::ZeroFunction<dim>(dim), 
                                             boundary_values); 
    VectorTools::interpolate_boundary_values( 
      dof_handler, 
      1, 
      IncrementalBoundaryValues<dim>(present_time, present_timestep), 
      boundary_values, 
      fe.component_mask(z_component)); 

    PETScWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 
    MatrixTools::apply_boundary_values( 
      boundary_values, system_matrix, tmp, system_rhs, false); 
    incremental_displacement = tmp; 
  } 

//  @sect4{TopLevel::solve_timestep}  

// 下一个函数是控制一个时间段内必须发生的所有事情的函数。从函数名称上看，事情的顺序应该是相对不言自明的。

  template <int dim> 
  void TopLevel<dim>::solve_timestep() 
  { 
    pcout << "    Assembling system..." << std::flush; 
    assemble_system(); 
    pcout << " norm of rhs is " << system_rhs.l2_norm() << std::endl; 

    const unsigned int n_iterations = solve_linear_problem(); 

    pcout << "    Solver converged in " << n_iterations << " iterations." 
          << std::endl; 

    pcout << "    Updating quadrature point data..." << std::flush; 
    update_quadrature_point_history(); 
    pcout << std::endl; 
  } 

//  @sect4{TopLevel::solve_linear_problem}  

// 再次求解线性系统的工作原理与之前基本相同。唯一不同的是，我们只想保留一份完整的本地解向量，而不是从PETSc的求解程序中得到的分布式向量。为此，我们为分布式向量声明一个本地临时变量，并用本地变量的内容对其进行初始化（记得 <code>apply_boundary_values</code> 中调用的 <code>assemble_system</code> 函数预设了该向量中边界节点的值），用它进行求解，并在函数结束时将其再次复制到我们声明为成员变量的完整本地向量中。然后，挂起的节点约束只分布在本地拷贝上，也就是说，在每个处理器上都是独立的。

  template <int dim> 
  unsigned int TopLevel<dim>::solve_linear_problem() 
  { 
    PETScWrappers::MPI::Vector distributed_incremental_displacement( 
      locally_owned_dofs, mpi_communicator); 
    distributed_incremental_displacement = incremental_displacement; 

    SolverControl solver_control(dof_handler.n_dofs(), 
                                 1e-16 * system_rhs.l2_norm()); 

    PETScWrappers::SolverCG cg(solver_control, mpi_communicator); 

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix); 

    cg.solve(system_matrix, 
             distributed_incremental_displacement, 
             system_rhs, 
             preconditioner); 

    incremental_displacement = distributed_incremental_displacement; 

    hanging_node_constraints.distribute(incremental_displacement); 

    return solver_control.last_step(); 
  } 

//  @sect4{TopLevel::output_results}  

// 这个函数生成.vtu格式的图形输出，正如介绍中所解释的。每个进程将只对其拥有的单元格进行工作，然后将结果写入自己的文件中。此外，处理器0将写下引用所有.vtu文件的记录文件。

// 这个函数的关键部分是给 <code>DataOut</code> 类提供一种方法，使其只对当前进程拥有的单元格进行工作。

  template <int dim> 
  void TopLevel<dim>::output_results() const 
  { 
    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
//然后，
//就像在 step-17 中一样，定义求解变量的名称（这里是位移增量）并排队输出求解向量。请注意在下面的开关中，我们如何确保如果空间维度应该不被处理，我们抛出一个异常，说我们还没有实现这种情况（另一个防御性编程的案例）。

    std::vector<std::string> solution_names; 
    switch (dim) 
      { 
        case 1: 
          solution_names.emplace_back("delta_x"); 
          break; 
        case 2: 
          solution_names.emplace_back("delta_x"); 
          solution_names.emplace_back("delta_y"); 
          break; 
        case 3: 
          solution_names.emplace_back("delta_x"); 
          solution_names.emplace_back("delta_y"); 
          solution_names.emplace_back("delta_z"); 
          break; 
        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    data_out.add_data_vector(incremental_displacement, solution_names); 

// 接下来的事情是，我们想输出类似于我们在每个单元中存储的应力的平均规范。这看起来很复杂，因为在目前的处理器上，我们只在那些实际属于目前进程的单元格上存储正交点的应力。换句话说，我们似乎无法计算出所有单元的平均应力。然而，请记住，我们源自 <code>DataOut</code> 的类只迭代那些实际属于当前处理器的单元，也就是说，我们不必为所有其他单元计算任何东西，因为这些信息不会被触及。下面的小循环就是这样做的。我们将整个区块包围在一对大括号中，以确保迭代器变量不会在它们被使用的区块结束后仍然意外地可见。

    Vector<double> norm_of_stress(triangulation.n_active_cells()); 
    { 

// 在所有的单元格上循环...

      for (auto &cell : triangulation.active_cell_iterators()) 
        if (cell->is_locally_owned()) 
          { 

// 在这些单元上，将所有正交点的应力相加...

            SymmetricTensor<2, dim> accumulated_stress; 
            for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
              accumulated_stress += 
                reinterpret_cast<PointHistory<dim> *>(cell->user_pointer())[q] 
                  .old_stress; 

// ...然后把平均值的常数写到它们的目的地。

            norm_of_stress(cell->active_cell_index()) = 
              (accumulated_stress / quadrature_formula.size()).norm(); 
          } 

// 在我们不感兴趣的单元格上，将向量中各自的值设置为一个假值（规范必须是正值，大的负值应该能吸引你的眼球），以确保如果我们的假设有误，即这些元素不会出现在输出文件中，我们会通过观察图形输出发现。

        else 
          norm_of_stress(cell->active_cell_index()) = -1e+20; 
    } 

// 最后把这个向量也附在上面，以便进行输出处理。

    data_out.add_data_vector(norm_of_stress, "norm_of_stress"); 

// 作为最后一个数据，如果这是一个并行作业，让我们也把域划分为与处理器相关的子域。这与 step-17 程序中的工作方式完全相同。

    std::vector<types::subdomain_id> partition_int( 
      triangulation.n_active_cells()); 
    GridTools::get_subdomain_association(triangulation, partition_int); 
    const Vector<double> partitioning(partition_int.begin(), 
                                      partition_int.end()); 
    data_out.add_data_vector(partitioning, "partitioning"); 

// 最后，有了这些数据，我们可以指示deal.II对信息进行整合，并产生一些中间数据结构，其中包含所有这些解决方案和其他数据向量。

    data_out.build_patches(); 

// 让我们调用一个函数，打开必要的输出文件，将我们生成的数据写入其中。该函数根据给定的目录名（第一个参数）和文件名基数（第二个参数）自动构建文件名。它通过由时间步数和 "片数 "产生的片断来增加所产生的字符串，"片数 "对应于整个域的一部分，可以由一个或多个子域组成。

// 该函数还为Paraview写了一个记录文件（后缀为`.pvd`），描述了所有这些输出文件如何组合成这个单一时间步骤的数据。

    const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
      "./", "solution", timestep_no, mpi_communicator, 4); 

// 记录文件必须只写一次，而不是由每个处理器来写，所以我们在0号处理器上做这个。

    if (this_mpi_process == 0) 
      { 

// 最后，我们写入paraview记录，它引用了所有.pvtu文件和它们各自的时间。注意，变量times_and_names被声明为静态的，所以它将保留前几个时间段的条目。

        static std::vector<std::pair<double, std::string>> times_and_names; 
        times_and_names.push_back( 
          std::pair<double, std::string>(present_time, pvtu_filename)); 
        std::ofstream pvd_output("solution.pvd"); 
        DataOutBase::write_pvd_record(pvd_output, times_and_names); 
      } 
  } 

//  @sect4{TopLevel::do_initial_timestep}  

// 这个函数和下一个函数分别处理第一个和下一个时间步骤的整体结构。第一个时间步骤的工作量稍大，因为我们要在连续细化的网格上多次计算，每次都从一个干净的状态开始。在这些计算的最后，我们每次都计算增量位移，我们使用最后得到的增量位移的结果来计算产生的应力更新并相应地移动网格。在这个新的网格上，我们再输出解决方案和任何我们认为重要的附加数据。

// 所有这些都会穿插着产生输出到控制台，以更新屏幕上的人正在发生的事情。如同在 step-17 中一样，使用 <code>pcout</code> instead of <code>std::cout</code> 可以确保只有一个并行进程实际在向控制台写数据，而不需要在每个产生输出的地方明确地编码一个if语句。

  template <int dim> 
  void TopLevel<dim>::do_initial_timestep() 
  { 
    present_time += present_timestep; 
    ++timestep_no; 
    pcout << "Timestep " << timestep_no << " at time " << present_time 
          << std::endl; 

    for (unsigned int cycle = 0; cycle < 2; ++cycle) 
      { 
        pcout << "  Cycle " << cycle << ':' << std::endl; 

        if (cycle == 0) 
          create_coarse_grid(); 
        else 
          refine_initial_grid(); 

        pcout << "    Number of active cells:       " 
              << triangulation.n_active_cells() << " (by partition:"; 
        for (unsigned int p = 0; p < n_mpi_processes; ++p) 
          pcout << (p == 0 ? ' ' : '+') 
                << (GridTools::count_cells_with_subdomain_association( 
                     triangulation, p)); 
        pcout << ")" << std::endl; 

        setup_system(); 

        pcout << "    Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (by partition:"; 
        for (unsigned int p = 0; p < n_mpi_processes; ++p) 
          pcout << (p == 0 ? ' ' : '+') 
                << (DoFTools::count_dofs_with_subdomain_association(dof_handler, 
                                                                    p)); 
        pcout << ")" << std::endl; 

        solve_timestep(); 
      } 

    move_mesh(); 
    output_results(); 

    pcout << std::endl; 
  } 

//  @sect4{TopLevel::do_timestep}  

// 后续的时间步骤比较简单，鉴于上面对前一个函数的解释，可能不需要更多的文件。

  template <int dim> 
  void TopLevel<dim>::do_timestep() 
  { 
    present_time += present_timestep; 
    ++timestep_no; 
    pcout << "Timestep " << timestep_no << " at time " << present_time 
          << std::endl; 
    if (present_time > end_time) 
      { 
        present_timestep -= (present_time - end_time); 
        present_time = end_time; 
      } 

    solve_timestep(); 

    move_mesh(); 
    output_results(); 

    pcout << std::endl; 
  } 
// @sect4{TopLevel::refine_initial_grid}  

// 当在连续细化的网格上求解第一个时间步骤时，调用以下函数。每次迭代后，它都会计算一个细化准则，细化网格，并将每个正交点的历史变量再次设置为干净状态。

  template <int dim> 
  void TopLevel<dim>::refine_initial_grid() 
  { 

// 首先，让每个进程计算其拥有的单元格的误差指标。

    Vector<float> error_per_cell(triangulation.n_active_cells()); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(fe.degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      incremental_displacement, 
      error_per_cell, 
      ComponentMask(), 
      nullptr, 
      MultithreadInfo::n_threads(), 
      this_mpi_process); 

// 然后建立一个全局向量，我们将来自每个%并行进程的局部指标合并到其中。

    const unsigned int n_local_cells = 
      triangulation.n_locally_owned_active_cells(); 

    PETScWrappers::MPI::Vector distributed_error_per_cell( 
      mpi_communicator, triangulation.n_active_cells(), n_local_cells); 

    for (unsigned int i = 0; i < error_per_cell.size(); ++i) 
      if (error_per_cell(i) != 0) 
        distributed_error_per_cell(i) = error_per_cell(i); 
    distributed_error_per_cell.compress(VectorOperation::insert); 

// 一旦我们有了这个，就把它复制回所有处理器上的本地副本，并相应地完善网格。

    error_per_cell = distributed_error_per_cell; 
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    error_per_cell, 
                                                    0.35, 
                                                    0.03); 
    triangulation.execute_coarsening_and_refinement(); 

// 最后，在新的网格上再次设置正交点数据，并且只在那些我们已经确定是我们的单元上设置。

    setup_quadrature_point_history(); 
  } 

//  @sect4{TopLevel::move_mesh}  

// 在每个时间步骤结束时，我们根据这个时间步骤计算的增量位移来移动网格的节点。为了做到这一点，我们保留一个标志的向量，为每个顶点指示我们是否已经移动过它，然后在所有单元中循环，移动那些尚未移动的单元顶点。值得注意的是，我们从某个顶点相邻的单元中移动这个顶点并不重要：因为我们使用连续有限元计算位移，位移场也是连续的，我们可以从每个相邻的单元中计算某个顶点的位移。我们只需要确保每个节点都精确地移动一次，这就是为什么我们要保留标志的矢量。

// 在这个函数中，有两个值得注意的地方。首先，我们如何使用 <code>cell-@>vertex_dof_index(v,d)</code> 函数获得给定顶点的位移场，该函数返回给定单元的 <code>d</code>th degree of freedom at vertex <code>v</code> 的索引。在本例中，k-th坐标方向的位移对应于有限元的k-th分量。使用这样的函数有一定的风险，因为它使用了我们在 <code>FESystem</code> 元素中为这个程序共同采取的元素顺序的知识。如果我们决定增加一个额外的变量，例如用于稳定的压力变量，并碰巧将其作为元素的第一个变量插入，那么下面的计算将开始产生无意义的结果。此外，这种计算还依赖于其他假设：首先，我们使用的元素确实有与顶点相关的自由度。对于目前的Q1元素来说确实如此，对于所有多项式阶的Qp元素来说也是如此  <code>p</code>  。然而，这对不连续的元素或混合公式的元素来说是不成立的。其次，它还建立在这样的假设上：一个顶点的位移只由与这个顶点相关的自由度的值决定；换句话说，所有对应于其他自由度的形状函数在这个特定的顶点是零。同样，对于目前的元素来说是这样的，但对于目前在deal.II中的所有元素来说并非如此。尽管有风险，我们还是选择使用这种方式，以便提出一种查询与顶点相关的单个自由度的方法。

// 在这种情况下，指出一种更普遍的方法是很有意义的。对于一般的有限元来说，应该采用正交公式，将正交点放在单元的顶点上。梯形规则的 <code>QTrapezoid</code> 公式正是这样做的。有了这个正交公式，我们就可以在每个单元格中初始化一个 <code>FEValues</code> 对象，并使用 <code>FEValues::get_function_values</code> 函数来获得正交点，即单元格顶点的解函数值。这些是我们真正需要的唯一数值，也就是说，我们对与这个特定正交公式相关的权重（或 <code>JxW</code> 值）完全不感兴趣，这可以作为 <code>FEValues</code> 构造器的最后一个参数来指定。这个方案中唯一的一点小麻烦是，我们必须弄清楚哪个正交点对应于我们目前考虑的顶点，因为它们可能是以相同的顺序排列，也可能不是。

// 如果有限元在顶点上有支持点（这里的支持点是有的；关于支持点的概念，见 @ref GlossSupport "支持点"），这种不便就可以避免了。对于这种情况，我们可以使用 FiniteElement::get_unit_support_points(). 构建一个自定义的正交规则，然后第一个 <code>cell-&gt;n_vertices()*fe.dofs_per_vertex</code> 正交点将对应于单元格的顶点，其顺序与 <code>cell-@>vertex(i)</code> 一致，同时考虑到矢量元素的支持点将被重复 <code>fe.dofs_per_vertex</code> 次。

// 关于这个短函数值得解释的另一点是三角形类输出其顶点信息的方式：通过 <code>Triangulation::n_vertices</code> 函数，它公布了三角形中有多少个顶点。并非所有的顶点都是一直在使用的--有些是之前被粗化的单元的遗留物，自从deal.II以来一直存在，一旦一个顶点出现，即使数量较少的顶点消失了，也不会改变它的编号。其次， <code>cell-@>vertex(v)</code> 返回的位置不仅是一个类型为 <code>Point@<dim@></code> 的只读对象，而且事实上是一个可以写入的引用。这允许相对容易地移动网格的节点，但值得指出的是，使用该功能的应用程序有责任确保所得到的单元仍然有用，即没有扭曲到单元退化的程度（例如，用负的雅各布系数表示）。请注意，我们在这个函数中没有任何规定来实际保证这一点，我们只是有信心。

// 在这个冗长的介绍之后，下面是全部20行左右的代码。

  template <int dim> 
  void TopLevel<dim>::move_mesh() 
  { 
    pcout << "    Moving mesh..." << std::endl; 

    std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 
    for (auto &cell : dof_handler.active_cell_iterators()) 
      for (const auto v : cell->vertex_indices()) 
        if (vertex_touched[cell->vertex_index(v)] == false) 
          { 
            vertex_touched[cell->vertex_index(v)] = true; 

            Point<dim> vertex_displacement; 
            for (unsigned int d = 0; d < dim; ++d) 
              vertex_displacement[d] = 
                incremental_displacement(cell->vertex_dof_index(v, d)); 

            cell->vertex(v) += vertex_displacement; 
          } 
  } 
// @sect4{TopLevel::setup_quadrature_point_history}  

// 在计算的开始，我们需要设置历史变量的初始值，例如材料中的现有应力，我们将其存储在每个正交点中。如上所述，我们使用每个单元中都有的 <code>user_pointer</code> 来做这个。

// 为了从更大的角度看这个问题，我们注意到，如果我们的模型中有先前可用的应力（为了这个程序的目的，我们假定这些应力不存在），那么我们就需要将先前存在的应力场插值到正交点上。同样，如果我们要模拟具有硬化/软化的弹塑性材料，那么我们就必须在每个正交点存储额外的历史变量，如累积塑性应变的当前屈服应力。预先存在的硬化或弱化也将通过在当前函数中插值这些变量来实现。

  template <int dim> 
  void TopLevel<dim>::setup_quadrature_point_history() 
  { 

// 为了慎重起见，我们把所有单元格的用户指针，不管是不是我们的，都设置为空指针。这样，如果我们访问了不应该访问的单元格的用户指针，一个分段故障将让我们知道这不应该发生。

    triangulation.clear_user_data(); 

// 接下来，分配属于这个处理器职责范围内的正交对象。当然，这等于属于这个处理器的单元格的数量乘以我们的正交公式在每个单元格上的正交点的数量。由于`resize()`函数在要求的新大小小于旧大小的情况下，实际上并没有缩小分配的内存量，所以我们采用了一个技巧，首先释放所有的内存，然后再重新分配：我们声明一个空向量作为临时变量，然后交换旧向量和这个临时变量的内容。这就确保了`正交点历史'现在确实是空的，我们可以让现在保存着以前的向量内容的临时变量超出范围并被销毁。在下一步中，我们可以根据需要重新分配尽可能多的元素，矢量默认初始化`PointHistory`对象，这包括将压力变量设置为零。

    { 
      std::vector<PointHistory<dim>> tmp; 
      quadrature_point_history.swap(tmp); 
    } 
    quadrature_point_history.resize( 
      triangulation.n_locally_owned_active_cells() * quadrature_formula.size()); 

// 最后再次循环所有单元，并将属于本处理器的单元的用户指针设置为指向此类对象的向量中与本单元对应的第一个正交点对象。

    unsigned int history_index = 0; 
    for (auto &cell : triangulation.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 
          cell->set_user_pointer(&quadrature_point_history[history_index]); 
          history_index += quadrature_formula.size(); 
        } 

// 最后，为了慎重起见，确保我们对元素的计数是正确的，而且我们已经用完了之前分配的所有对象，并且没有指向任何超出向量末端的对象。这样的防御性编程策略总是很好的检查，以避免意外的错误，并防止将来对这个函数的修改忘记同时更新一个变量的所有用途。回顾一下，使用 <code>Assert</code> 宏的构造在优化模式下被优化掉了，所以不影响优化运行的运行时间。

    Assert(history_index == quadrature_point_history.size(), 
           ExcInternalError()); 
  } 

//  @sect4{TopLevel::update_quadrature_point_history}  

// 在每个时间步骤结束时，我们应该计算出一个增量的位移更新，使材料在其新的配置中能够容纳这个时间步骤中施加的外部体和边界力减去通过预先存在的内部应力施加的力之间的差异。为了在下一个时间步骤中获得预先存在的应力，我们必须用本时间步骤中计算的增量位移引起的应力来更新预先存在的应力。理想情况下，所产生的内应力之和将完全抵消所有的外力。事实上，一个简单的实验可以确保这一点：如果我们选择边界条件和体力与时间无关，那么强迫项（外力和内应力之和）应该正好是零。如果你做了这个实验，你会从每个时间步长的右手边的规范输出中意识到这几乎是事实：它并不完全是零，因为在第一个时间步长中，增量位移和应力的更新是相对于未变形的网格计算的，然后再进行变形。在第二个时间步骤中，我们再次计算位移和应力的更新，但这次是在变形的网格中 -- 在那里，结果的更新非常小但不完全是零。这可以迭代，在每一次迭代中，残差，即右手边向量的法线，都会减少；如果做这个小实验，就会发现这个残差的法线会随着迭代次数的增加而呈指数下降，在最初的快速下降之后，每次迭代大约会减少3.5倍（对于我看的一个测试案例，其他测试案例和其他未知数都会改变这个系数，但不会改变指数下降的情况）。

// 在某种意义上，这可以被认为是一个准时序方案，以解决在一个以拉格朗日方式移动的网格上解决大变形弹性的非线性问题。

// 另一个复杂的问题是，现有的（旧的）应力是在旧的网格上定义的，我们将在更新应力后移动这个网格。如果这个网格的更新涉及到单元的旋转，那么我们也需要对更新的应力进行旋转，因为它是相对于旧单元的坐标系计算的。

// 因此，我们需要的是：在当前处理器拥有的每个单元上，我们需要从每个正交点存储的数据中提取旧的应力，计算应力更新，将两者相加，然后将结果与从当前正交点的增量位移计算出来的增量旋转一起旋转。下面我们将详细介绍这些步骤。

  template <int dim> 
  void TopLevel<dim>::update_quadrature_point_history() 
  { 

// 首先，建立一个 <code>FEValues</code> 对象，我们将通过它来评估正交点的增量位移及其梯度，还有一个保存这些信息的向量。

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_gradients); 

    std::vector<std::vector<Tensor<1, dim>>> displacement_increment_grads( 
      quadrature_formula.size(), std::vector<Tensor<1, dim>>(dim)); 

// 然后在所有单元格上循环，在属于我们子域的单元格中进行工作。

    for (auto &cell : dof_handler.active_cell_iterators()) 
      if (cell->is_locally_owned()) 
        { 

// 接下来，获得一个指向当前单元本地正交点历史数据的指针，作为防御措施，确保这个指针在全局数组的范围内。

          PointHistory<dim> *local_quadrature_points_history = 
            reinterpret_cast<PointHistory<dim> *>(cell->user_pointer()); 
          Assert(local_quadrature_points_history >= 
                   &quadrature_point_history.front(), 
                 ExcInternalError()); 
          Assert(local_quadrature_points_history <= 
                   &quadrature_point_history.back(), 
                 ExcInternalError()); 

// 然后在本单元上初始化 <code>FEValues</code> 对象，并提取正交点上的位移梯度，以便以后计算应变。

          fe_values.reinit(cell); 
          fe_values.get_function_gradients(incremental_displacement, 
                                           displacement_increment_grads); 

// 然后在这个单元的正交点上循环。

          for (unsigned int q = 0; q < quadrature_formula.size(); ++q) 
            { 

// 在每个正交点上，从梯度中计算出应变增量，并将其乘以应力-应变张量，得到应力更新。然后将此更新添加到该点已有的应变中。

              const SymmetricTensor<2, dim> new_stress = 
                (local_quadrature_points_history[q].old_stress + 
                 (stress_strain_tensor * 
                  get_strain(displacement_increment_grads[q]))); 

// 最后，我们要对结果进行旋转。为此，我们首先要从增量位移中计算出目前正交点的旋转矩阵。事实上，它可以从梯度中计算出来，而且我们已经有一个函数用于这个目的。

              const Tensor<2, dim> rotation = 
                get_rotation_matrix(displacement_increment_grads[q]); 

// 注意这个结果，即旋转矩阵，一般来说是一个等级为2的反对称张量，所以我们必须把它作为一个完整的张量来存储。

// 有了这个旋转矩阵，在我们将对称张量 <code>new_stress</code> 扩展为全张量之后，我们可以通过从左和右的收缩来计算旋转的张量。

              const SymmetricTensor<2, dim> rotated_new_stress = 
                symmetrize(transpose(rotation) * 
                           static_cast<Tensor<2, dim>>(new_stress) * rotation); 

// 注意，虽然这三个矩阵的乘法结果应该是对称的，但由于浮点舍入的原因，它并不是对称的：我们得到的结果的非对角线元素有1e-16的不对称性。当把结果赋给一个 <code>SymmetricTensor</code> 时，该类的构造函数会检查对称性并意识到它不是完全对称的；然后它会引发一个异常。为了避免这种情况，我们明确地对结果进行对称，使其完全对称。

// 所有这些操作的结果会被写回到原来的地方。

              local_quadrature_points_history[q].old_stress = 
                rotated_new_stress; 
            } 
        } 
  } 

// 这就结束了项目特定的命名空间  <code>Step18</code>  。其余的和往常一样，并且在  step-17  中已经显示：一个  <code>main()</code>  函数初始化和终止 PETSc，调用做实际工作的类，并确保我们捕捉所有传播到这一点的异常。

} // namespace Step18 

int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step18; 

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 

      TopLevel<3> elastic_problem; 
      elastic_problem.run(); 
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


