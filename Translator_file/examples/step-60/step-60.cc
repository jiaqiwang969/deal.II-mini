CCTest_file/step-60.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2018 - 2020 by the deal.II authors 
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
 * Authors: Luca Heltai, Giovanni Alzetta, 
 * International School for Advanced Studies, Trieste, 2018 
 */ 


// @sect3{Include files}  其中大部分已经在其他地方介绍过了，我们只对新的部分进行评论。

#include <deal.II/base/logstream.h> 
#include <deal.II/base/utilities.h> 
#include <deal.II/base/timer.h> 

// 参数接受器类是本教程程序的第一个新颖之处：一般来说，参数文件是用来在运行时引导程序的执行。虽然即使是简单的方法也能节省编译时间，因为同一个可执行文件可以用不同的参数设置来运行，但要同时处理数百个参数，同时保持不同程序之间的兼容性，会变得很困难。这就是ParameterAcceptor类证明有用的地方。

// 该类用于定义一个公共接口，供那些希望使用一个全局ParameterHandler来处理参数的类使用。该类提供了一个静态的ParameterHandler成员，即 ParameterAcceptor::prm, ，并实现了 "命令设计模式"（例如，见E. Gamma, R. Helm, R. Johnson, J. Vlissides, 设计模式。Elements of Reusable Object-Oriented Software, Addison-Wesley Professional, 1994. https:goo.gl/NYByc）。)

// ParameterAcceptor提供了一个全局订阅机制。每当一个从ParameterAcceptor派生出来的类的对象被构造出来时，一个指向该派生型对象的指针就会被注册，同时在参数文件中还有一个部分条目。这种注册表在调用单一函数 ParameterAcceptor::initialize("file.prm") 时被遍历，该函数反过来确保所有存储在全局注册表中的类都声明它们将使用的参数，在声明了这些参数后，它读取`file.prm`的内容来解析实际参数。

// 如果你为你想在代码中使用的每个参数调用方法 ParameterHandler::add_parameter ，你就不需要做其他事情了。如果你使用的是一个已经存在的类，它提供了`declare_parameters`和`parse_parameters`这两个函数，你仍然可以使用ParameterAcceptor，方法是将现有的类封装成ParameterAcceptorProxy类。

// 在这个例子中，我们将同时使用这两种策略，为deal.II类使用ParameterAcceptorProxy，并直接从ParameterAcceptor派生出我们自己的参数类。

#include <deal.II/base/parameter_acceptor.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 

// 另一个新的包含文件是包含 GridTools::Cache 类的文件。deal.II的结构和许多现代数值库一样，是按照有向无环图（DAG）组织的。DAG是一个具有拓扑顺序的有向图：每个节点在结构上代表一个对象，并通过一条（或多条）有向边与非根节点相连，从父节点到子节点。这种结构最重要的例子是Triangulation及其 Triangulation::cell_iterator 结构。从Triangulation（主节点），我们可以访问每个单元（三角形的子节点）。从单元本身，我们可以访问该单元的所有顶点。在这个简单的例子中，DAG结构可以表示为三种节点类型（三角结构、单元格迭代器和顶点），由从三角结构到单元格迭代器，以及从单元格迭代器到顶点的定向边连接。这有几个优点，但它本质上造成了 "不对称"，使某些操作变得很快，而它们的逆向操作却很慢：寻找一个单元的顶点的计算成本很低，可以通过简单地遍历DAG来完成，而寻找所有共享一个顶点的单元则需要进行非难事的计算，除非添加一个新的DAG数据结构来表示逆向搜索。

// 由于在有限元代码中通常不需要逆向操作，所以在GridTools中实现了这些操作，而没有使用与Triangulation相关的额外数据结构，这将使它们的速度更快。例如，这样的数据结构是一个从三角结构的顶点到所有共享这些顶点的单元的映射，这将减少回答前面问题所需的计算。

// 有些方法，例如 GridTools::find_active_cell_around_point, 大量使用了这些非标准操作。如果你需要多次调用这些方法，那么将这些数据结构存储在某个地方就变得很方便。  GridTools::Cache 正是这样做的，它让你可以访问以前计算过的对象，或者在飞行中计算它们（然后将它们存储在类中供以后使用），并确保每当三角测量被更新时，相关的数据结构也被重新计算。

#include <deal.II/grid/grid_tools_cache.h> 

#include <deal.II/fe/fe.h> 
#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 

// 在这个例子中，我们将使用一个参考域来描述一个嵌入的三角形，通过一个有限元矢量场进行变形。

// 接下来的两个包含文件包含了在这些情况下可以使用的两个类的定义。MappingQEulerian允许人们通过一个*位移*场来描述一个域，基于FESystem[FE_Q(p)^spacedim]有限元空间。第二种是比较通用的，允许你使用任意的矢量有限元空间，只要它们能提供*连续*的
//对你的领域的描述。在这种情况下，这种描述是通过实际的*变形*场，而不是*位移*场来完成的。

// 使用哪一个取决于用户想如何指定参考域，和/或实际配置。我们将提供这两种选择，并在本教程程序的结果部分做一些实验。

#include <deal.II/fe/mapping_q_eulerian.h> 
#include <deal.II/fe/mapping_fe_field.h> 

#include <deal.II/dofs/dof_tools.h> 

// 被解析的函数类是另一个新条目。它允许人们创建一个Function对象，从参数文件中的一个字符串开始，它被解析成一个对象，你可以在deal.II接受Function的任何地方使用（例如，用于插值、边界条件等）。

#include <deal.II/base/parsed_function.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 

// 这是本教程程序的最后一个新条目。NonMatching命名空间包含了一些方法，这些方法在对非匹配网格或与底层网格不对齐的曲线进行计算时很有用。

// 我们将在后面的 "setup_coupling "方法中详细讨论其用途。

#include <deal.II/non_matching/coupling.h> 

#include <deal.II/lac/affine_constraints.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/sparse_direct.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/linear_operator.h> 
#include <deal.II/lac/linear_operator_tools.h> 

#include <iostream> 
#include <fstream> 

namespace Step60 
{ 
  using namespace dealii; 
// @sect3{DistributedLagrangeProblem}  

// 在DistributedLagrangeProblem中，我们需要两个参数来描述域 $\Gamma$ （`dim`）和域 $\Omega$ （`spacedim`）的尺寸。

// 这些参数将用于初始化一个Triangulation<dim,spacedim>（用于 $\Gamma$  ）和一个Triangulation<spacedim,spacedim>（用于 $\Omega$  ）。

// 与其他教程程序相比，一个新奇之处在于大量使用了  std::unique_ptr.  这些指针的行为与经典指针一样，其优点是可以进行自动的内部管理：一旦unique_ptr超出范围，所包含的对象就会被自动销毁，即使它是在一个容器内或者有一个异常。此外，它不允许有重复的指针，这可以防止所有权问题。我们这样做，是因为我们希望能够 i) 构建问题，ii) 读取参数，iii) 根据参数文件中指定的内容初始化所有对象。

// 我们在派生于ParameterAcceptor的内部类`Parameters`中构建我们问题的参数。`DistributedLagrangeProblem`类需要一个对`Parameters`对象的常量引用，因此不可能从DistributedLagrangeProblem类本身修改参数。

// 我们可以先初始化参数，然后将参数传递给DistributedLagrangeProblem，假设所有条目都设置为所需的值，但这有两个缺点。



// - 我们不应该对用户如何初始化一个不受我们直接控制的类做出假设。如果用户未能初始化该类，我们应该注意到并抛出一个异常。



// - 当我们构造参数时，并不是所有需要从参数文件中读取参数的对象都可以使用；对于复杂的程序，有多种物理现象，或者我们在一些外部类中重复使用现有的代码，往往是这种情况。我们通过将一些 "复杂 "的对象，如ParsedFunction对象，保留在`DistributedLagrangeProblem`内而不是`Parameters`内来模拟这种情况。

// 这里我们假设在构建时，构建我们问题的类还不能使用。解析参数文件是确保我们有所有的成分来建立我们的类，我们的设计是，如果解析失败，或者没有被执行，运行就会被中止。

  template <int dim, int spacedim = dim> 
  class DistributedLagrangeProblem 
  { 
  public: 

// `Parameters`类是由ParameterAcceptor派生的。这使得我们可以在其构造函数中使用 ParameterAcceptor::add_parameter() 方法。

// 这个函数的成员都是非常量的，但是`DistributedLagrangeProblem`类需要一个对`Parameters`对象的常量引用：这确保参数不会从`DistributedLagrangeProblem`类中被修改。

    class Parameters : public ParameterAcceptor 
    { 
    public: 
      Parameters(); 

// 现在描述的参数都可以通过参数文件在外部设置：如果运行可执行文件时没有参数文件，程序将创建一个 "参数.prm "文件，并在这里定义默认值，然后中止，让用户有机会修改参数.prm文件。

// 嵌入网格的初始细化，对应于域 $\Omega$  。

      unsigned int initial_refinement = 4; 

// 嵌入网格 $\Omega$ 和嵌入网格 $\Gamma$ 之间的交互是通过 $C$ 的计算来处理的，这涉及到 $\Omega$ 的所有单元与 $\Gamma$ 的部分重叠：这些单元的更高细化可能提高我们计算的质量。出于这个原因，我们定义了`delta_refinement'：如果它大于零，那么我们将空间网格中包含嵌入网格顶点的每个单元及其邻居标记出来，执行细化，并重复这个过程`delta_refinement'的次数。

      unsigned int delta_refinement = 3; 

// 开始细化嵌入网格，对应于域  $\Gamma$  。

      unsigned int initial_embedded_refinement = 8; 

// 边界id列表，我们在其中施加同质Dirichlet边界条件。在其余的边界id上（如果有的话），我们施加同质的诺伊曼边界条件。作为一个默认问题，我们在 $\partial \Omega$ 上设置了零迪里希特边界条件。
      std::list<types::boundary_id> homogeneous_dirichlet_ids{0, 1, 2, 3}; 

// 嵌入空间的有限元素程度。  $V_h(\Omega)$  
      unsigned int embedding_space_finite_element_degree = 1; 

// 嵌入空间的有限元素度。  $Q_h(\Gamma)$  
      unsigned int embedded_space_finite_element_degree = 1; 

// 用于描述嵌入域变形的空间的有限元度

      unsigned int embedded_configuration_finite_element_degree = 1; 

// 用于积分耦合的正交公式的阶数

      unsigned int coupling_quadrature_order = 3; 

// 如果设置为 "true"，则嵌入式配置函数被解释为位移函数。

      bool use_displacement = false; 

// 在输出中使用的粗略程度

      unsigned int verbosity_level = 10; 

// 用来记录我们是否被初始化的一个标志。

      bool initialized = false; 
    }; 

    DistributedLagrangeProblem(const Parameters &parameters); 

// 分布式拉格朗日问题的入口点

    void run(); 

  private: 

// 包含实际参数的对象

    const Parameters &parameters; 

// 下面的函数与其他所有的教程程序相似，不同的是，我们现在需要为两个不同系列的对象设置东西，即与*嵌入*网格有关的对象，以及与*嵌入*有关的对象。

    void setup_grids_and_dofs(); 

    void setup_embedding_dofs(); 

    void setup_embedded_dofs(); 

// 这里唯一的非常规函数是`setup_coupling()`方法，用于生成耦合矩阵的稀疏性模式  $C$  。

    void setup_coupling(); 

    void assemble_system(); 

    void solve(); 

    void output_results(); 

// 首先我们收集所有与嵌入空间几何学有关的对象

    std::unique_ptr<Triangulation<spacedim>> space_grid; 
    std::unique_ptr<GridTools::Cache<spacedim, spacedim>> 
                                             space_grid_tools_cache; 
    std::unique_ptr<FiniteElement<spacedim>> space_fe; 
    std::unique_ptr<DoFHandler<spacedim>>    space_dh; 

// 然后是与嵌入网格有关的，与拉格朗日乘数`lambda`有关的DoFHandler

    std::unique_ptr<Triangulation<dim, spacedim>> embedded_grid; 
    std::unique_ptr<FiniteElement<dim, spacedim>> embedded_fe; 
    std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_dh; 

// 最后，需要对嵌入的三角形进行*变形的所有内容

    std::unique_ptr<FiniteElement<dim, spacedim>> embedded_configuration_fe; 
    std::unique_ptr<DoFHandler<dim, spacedim>>    embedded_configuration_dh; 
    Vector<double>                                embedded_configuration; 

// ParameterAcceptorProxy类是一个 "透明 "的包装器，它同时来自ParameterAcceptor和作为其模板参数传递的类型。在构造时，参数被分成两部分：第一个参数是一个 std::string, ，转发给ParameterAcceptor类，并包含应该用于该类的部分名称，而所有其余的参数都转发给模板类型的构造器，在本例中，转发给 Functions::ParsedFunction 构造器。

// 这个类允许你结合ParameterAcceptor注册机制使用现有的类，只要这些类有`declare_parameters()`和`parse_parameters()`成员。

// 这里就是这种情况，使得利用 Functions::ParsedFunction 类变得相当容易：而不是要求用户在代码中为RHS、边界函数等创建新的Function对象。在这里，我们允许用户使用deal.II接口到muParser（http:muparser.beltoforion.de），在这里，函数的指定不是在编译时完成的，而是在运行时完成的，使用一个字符串被解析成一个实际的Function对象。

// 在这种情况下，`embedded_configuration_function`是一个向量值的Function，根据`parameters.use_displacement`的布尔值，可以被解释为*变形*或*位移*。分量的数量将在后面的结构中指定。

    ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> 
      embedded_configuration_function; 

    std::unique_ptr<Mapping<dim, spacedim>> embedded_mapping; 

// 我们做同样的事情来指定函数的值  $g$  ，这就是我们希望我们的解决方案在嵌入空间中的值。在这种情况下，该函数是一个标量函数。

    ParameterAcceptorProxy<Functions::ParsedFunction<spacedim>> 
      embedded_value_function; 

// 与我们对 Functions::ParsedFunction 类所做的类似，我们对ReductionControl类重复同样的做法，允许我们为我们以后使用的Schur补码迭代求解器指定所有可能的停止标准。

    ParameterAcceptorProxy<ReductionControl> schur_solver_control; 

// 接下来我们收集所有我们需要的SparsityPattern, SparseMatrix, 和Vector对象

    SparsityPattern stiffness_sparsity; 
    SparsityPattern coupling_sparsity; 

    SparseMatrix<double> stiffness_matrix; 
    SparseMatrix<double> coupling_matrix; 

    AffineConstraints<double> constraints; 

    Vector<double> solution; 
    Vector<double> rhs; 

    Vector<double> lambda; 
    Vector<double> embedded_rhs; 
    Vector<double> embedded_value; 

// TimerOutput类用来提供一些关于我们程序性能的统计数据。

    TimerOutput monitor; 
  }; 
// @sect3{DistributedLagrangeProblem::Parameters}  

// 在构造时，我们也初始化ParameterAcceptor类，在解析参数文件时，我们希望问题使用的部分名称。

// 参数文件可以被组织成section/subsection/etc.：这样做的好处是，当共享同一section/subsection/etc.时，定义的对象可以共享参数。ParameterAcceptor允许使用Unix约定的路径来指定部分名称。如果部分名称以斜线（"/"）开头，那么该部分将被解释为*绝对路径*，ParameterAcceptor为路径中的每个目录输入一个小节，使用它遇到的最后一个名称作为当前类的着陆小节。

// 例如，如果你使用`ParameterAcceptor("/first/second/third/My Class")`构建你的类，参数将被组织成如下样子。

// 
// @code
//  # Example parameter file
//  subsection first
//    subsection second
//      subsection third
//        subsection My Class
//         ... # all the parameters
//        end
//      end
//    end
//  end
//  @endcode

// 在内部，存储在ParameterAcceptor中的*当前路径现在被认为是"/first/second/third/"，也就是说，当你指定一个绝对路径时，ParameterAcceptor将*当前的部分改为当前路径，即改为直到*最后一个"/"的部分名称的路径。

// 你现在可以使用相对路径（例如：`ParameterAcceptor("My Other Class")`）而不是绝对路径（例如：`ParameterAcceptor("/first/second/third/My Other Class")`）来构造从ParameterAcceptor派生的另一个类，获得。
// @code
//  # Example parameter file
//  subsection first
//    subsection second
//      subsection third
//        subsection My Class
//          ... # all the parameters
//        end
//        subsection My Other Class
//          ... # all the parameters of MyOtherClass
//        end
//      end
//    end
//  end
//  @endcode

// 如果部分名称*以斜线结尾，那么后续的类将把它解释为完整的路径：例如，与上面类似，如果我们有两个类，一个用`ParameterAcceptor("/first/second/third/My Class/")`初始化，另一个用`ParameterAcceptor("My Other Class")`，那么得到的参数文件将看起来像。

// 
// @code
//  # Example parameter file
//  subsection first
//    subsection second
//      subsection third
//        subsection My Class
//          ... # all the parameters of MyClass
//          ... # notice My Class subsection does not end here
//          subsection My Other Class
//            ... # all the parameters of MyOtherClass
//          end # of subsection My Other Class
//        end # of subsection My Class
//      end
//    end
//  end
//  @endcode

// 我们将利用这一点，使我们的`Parameters`成为所有后续构建的类的*父类。由于大多数其他的类都是`DistributedLagrangeProblem`的成员，这就允许，例如，为两个不同的维度构造两个`DistributedLagrangeProblem`，而不会在这两个问题的参数上产生冲突。

  template <int dim, int spacedim> 
  DistributedLagrangeProblem<dim, spacedim>::Parameters::Parameters() 
    : ParameterAcceptor("/Distributed Lagrange<" + 
                        Utilities::int_to_string(dim) + "," + 
                        Utilities::int_to_string(spacedim) + ">/") 
  { 

//  ParameterAcceptor::add_parameter() 函数做了几件事。



// - 在构造时向ParameterAcceptor输入指定的分段



// - 调用 ParameterAcceptor::prm.add_parameter() 函数



// - 调用你可能附加到 ParameterAcceptor::declare_parameters_call_back 的任何信号。



// --离开分节

// 反过来， ParameterAcceptor::prm.add_parameter  。



// - 在参数处理程序中为给定的变量声明一个条目。



// - 取得该变量的当前值



// - 将其转换为一个字符串，作为参数文件的默认值使用



// - 在 ParameterAcceptor::prm 中附加一个*动作，当文件被解析时，或者当一个条目被设置时，它会更新传递给`add_parameter()`的变量的值，将其设置为输入文件中指定的值（当然，是在输入文件被解析并将文本表示转换为变量的类型后）。

    add_parameter("Initial embedding space refinement", initial_refinement); 

    add_parameter("Initial embedded space refinement", 
                  initial_embedded_refinement); 

    add_parameter("Local refinements steps near embedded domain", 
                  delta_refinement); 

    add_parameter("Homogeneous Dirichlet boundary ids", 
                  homogeneous_dirichlet_ids); 

    add_parameter("Use displacement in embedded interface", use_displacement); 

    add_parameter("Embedding space finite element degree", 
                  embedding_space_finite_element_degree); 

 
                  embedded_space_finite_element_degree); 

    add_parameter("Embedded configuration finite element degree", 
                  embedded_configuration_finite_element_degree); 

    add_parameter("Coupling quadrature order", coupling_quadrature_order); 

    add_parameter("Verbosity level", verbosity_level); 

// 一旦参数文件被解析，那么参数就可以使用了。将内部变量`initialized`设置为true。

    parse_parameters_call_back.connect([&]() -> void { initialized = true; }); 
  } 

// 构造函数是非常标准的，除了前面解释的`ParameterAcceptorProxy`对象之外。

  template <int dim, int spacedim> 
  DistributedLagrangeProblem<dim, spacedim>::DistributedLagrangeProblem( 
    const Parameters &parameters) 
    : parameters(parameters) 
    , embedded_configuration_function("Embedded configuration", spacedim) 
    , embedded_value_function("Embedded value") 
    , schur_solver_control("Schur solver control") 
    , monitor(std::cout, TimerOutput::summary, TimerOutput::cpu_and_wall_times) 
  { 

// 下面是一个为使用ParameterAcceptorProxy构建的ParameterAcceptor类设置默认值的方法。

// 在这种情况下，我们将嵌入式网格的默认变形设置为半径为 $R$ 、中心为 $(Cx, Cy)$ 的圆，我们将embedded_value_function的默认值设置为常数，并为SolverControl对象指定一些合理的值。

// 嵌入 $\Gamma$ 是最基本的：从 $C_{\alpha j}$ 的定义可以看出，如果 $\Gamma \not\subseteq \Omega$ ，矩阵 $C$ 的某些行将是零。这将是一个问题，因为舒尔补码法要求 $C$ 具有全列秩。

    embedded_configuration_function.declare_parameters_call_back.connect( 
      []() -> void { 
        ParameterAcceptor::prm.set("Function constants", "R=.3, Cx=.4, Cy=.4"); 

        ParameterAcceptor::prm.set("Function expression", 
                                   "R*cos(2*pi*x)+Cx; R*sin(2*pi*x)+Cy"); 
      }); 

    embedded_value_function.declare_parameters_call_back.connect( 
      []() -> void { ParameterAcceptor::prm.set("Function expression", "1"); }); 

    schur_solver_control.declare_parameters_call_back.connect([]() -> void { 
      ParameterAcceptor::prm.set("Max steps", "1000"); 
      ParameterAcceptor::prm.set("Reduction", "1.e-12"); 
      ParameterAcceptor::prm.set("Tolerance", "1.e-12"); 
    }); 
  } 
// @sect3{Set up}  

//函数 `DistributedLagrangeProblem::setup_grids_and_dofs()` 是用来设置有限元空间的。注意 `std::make_unique` 是如何用来创建包裹在 `std::unique_ptr` 对象里面的对象的。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::setup_grids_and_dofs() 
  { 
    TimerOutput::Scope timer_section(monitor, "Setup grids and dofs"); 

// 初始化  $\Omega$  : 构建三角形，并将其包装成一个  `std::unique_ptr`  对象

    space_grid = std::make_unique<Triangulation<spacedim>>(); 

// 接下来，我们使用 GridGenerator::hyper_cube(). 实际创建三角形，最后一个参数设置为true：这激活了着色（即为边界的不同部分分配不同的边界指标），我们用它来分配Dirichlet和Neumann条件。

    GridGenerator::hyper_cube(*space_grid, 0, 1, true); 

// 一旦我们构建了一个三角形，我们就根据参数文件中的规格对其进行全局细化，并用它构建一个 GridTools::Cache 。

    space_grid->refine_global(parameters.initial_refinement); 
    space_grid_tools_cache = 
      std::make_unique<GridTools::Cache<spacedim, spacedim>>(*space_grid); 

// 对嵌入式网格也是这样做的。由于嵌入式网格是变形的，我们首先需要设置变形映射。我们在下面几行中这样做。

    embedded_grid = std::make_unique<Triangulation<dim, spacedim>>(); 
    GridGenerator::hyper_cube(*embedded_grid); 
    embedded_grid->refine_global(parameters.initial_embedded_refinement); 

    embedded_configuration_fe = std::make_unique<FESystem<dim, spacedim>>( 
      FE_Q<dim, spacedim>( 
        parameters.embedded_configuration_finite_element_degree), 
      spacedim); 

    embedded_configuration_dh = 
      std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid); 

    embedded_configuration_dh->distribute_dofs(*embedded_configuration_fe); 
    embedded_configuration.reinit(embedded_configuration_dh->n_dofs()); 

// 一旦我们为变形定义了一个有限维度的空间，我们就对参数文件中定义的 "嵌入_配置_函数 "进行插值。

    VectorTools::interpolate(*embedded_configuration_dh, 
                             embedded_configuration_function, 
                             embedded_configuration); 

// 现在我们可以根据用户在参数文件中指定的内容来解释它：作为位移，在这种情况下，我们构建一个映射，将我们的配置有限元空间的每个支持点的位置在相应的配置矢量上*移开指定的数量，或者作为赦免位置。

// 在第一种情况下，MappingQEulerian类提供其服务，而在第二种情况下，我们将使用MappingFEField类。它们实际上是非常相似的。MappingQEulerian只适用于FE_Q有限元空间系统，其中位移矢量存储在FESystem的第一个`spacedim`分量中，并且在构造时作为参数给出的度数必须与第一个`spacedim`分量的度数一致。

// MappingFEField类稍显一般，因为它允许你在构造近似值时选择任意的FiniteElement类型。当然，根据你选择的FiniteElement的类型，一些选择可能（也可能没有）意义。MappingFEField实现了纯粹的等参量概念，例如，可以通过与FE_Bernstein有限元类结合，在deal.II中实现等参量分析代码。在这个例子中，我们将考虑到一个配置将是一个 "位移"，而另一个将是一个绝对的 "变形 "场，从而将这两者互换使用。

    if (parameters.use_displacement == true) 
      embedded_mapping = 
        std::make_unique<MappingQEulerian<dim, Vector<double>, spacedim>>( 
          parameters.embedded_configuration_finite_element_degree, 
          *embedded_configuration_dh, 
          embedded_configuration); 
    else 
      embedded_mapping = 
        std::make_unique<MappingFEField<dim, spacedim, Vector<double>>>( 
          *embedded_configuration_dh, embedded_configuration); 

    setup_embedded_dofs(); 

// 在本教程中，我们不仅对 $\Omega$ 进行全局细化，还允许根据 $\Gamma$ 的位置进行局部细化，根据`parameters.delta_refinement`的值，我们用来决定对 $\Omega$ 进行多少轮局部细化，对应于 $\Gamma$  的位置。

// 有了这个映射，现在就可以通过调用方法 DoFTools::map_dofs_to_support_points. 来查询与`embedded_dh`相关的所有支持点的位置了。

// 这个方法有两个变种。一种是不接受Mapping，另一种是接受Mapping。如果你使用第二种类型，就像我们在这种情况下所做的那样，支持点是通过指定的映射来计算的，它可以对它们进行相应操作。

// 这正是`embedded_mapping`的作用。

    std::vector<Point<spacedim>> support_points(embedded_dh->n_dofs()); 
    if (parameters.delta_refinement != 0) 
      DoFTools::map_dofs_to_support_points(*embedded_mapping, 
                                           *embedded_dh, 
                                           support_points); 

// 一旦我们有了嵌入有限元空间的支持点，我们就想确定嵌入空间的哪些单元包含哪些支持点，以便在有必要的地方，也就是嵌入网格的地方，获得完善嵌入网格的机会。这可以手动完成，在每个支持点上循环，然后为嵌入空间的每个单元调用方法 Mapping::transform_real_to_unit_cell ，直到我们找到一个返回单位参考单元中的点，或者可以用更智能的方式完成。

//  GridTools::find_active_cell_around_point 是一个可能的选择，它以更便宜的方式执行上述任务，首先确定嵌入三角的最接近目标点的顶点，然后只对那些共享找到的顶点的单元格调用 Mapping::transform_real_to_unit_cell 。

// 事实上，在GridTools命名空间中，有一些算法利用 GridTools::Cache 对象，可能还有KDTree对象来尽可能地加快这些操作。

// 利用最大速度的最简单方法是调用一个专门的方法， GridTools::compute_point_locations, ，该方法将在第一个点的搜索过程中存储大量有用的信息和数据结构，然后在随后的点中重复使用所有这些信息。

//  GridTools::compute_point_locations 返回一个元组，其中第一个元素是一个包含输入点的单元格向量，在这里是support_points。对于细化来说，这是我们唯一需要的信息，而这正是现在所发生的。

// 然而，当我们需要组装一个耦合矩阵时，我们还需要每个点的参考位置来评估嵌入空间的基础函数。由 GridTools::compute_point_locations 返回的元组的其他元素允许你重建，对于每个点，什么单元包含它，以及什么是给定点的参考单元的位置。由于这些信息最好被分组到单元格中，那么这就是算法返回的内容：一个元组，包含所有单元格中至少有一个点的向量，以及所有参考点的列表和它们在原始向量中的相应索引。

// 在下面的循环中，我们将忽略所有返回的对象，除了第一个，确定所有单元格至少包含一个嵌入空间的支持点。这允许一个简单的自适应细化策略：细化这些单元和它们的邻居。

// 注意，我们需要做一些理智的检查，在这个意义上，我们希望有一个嵌入网格，它在嵌入网格周围被很好地细化，但其中两个连续的支持点要么位于同一个单元，要么位于邻近的嵌入单元。

// 这只有在我们确保嵌入网格的最小单元尺寸仍然大于嵌入网格的最大单元尺寸时才有可能。由于用户可以修改细化水平，以及他们希望在嵌入网格周围进行的局部细化的数量，我们要确保所得到的网格满足我们的要求，如果不是这样，我们就用一个例外来放弃。

    for (unsigned int i = 0; i < parameters.delta_refinement; ++i) 
      { 
        const auto point_locations = 
          GridTools::compute_point_locations(*space_grid_tools_cache, 
                                             support_points); 
        const auto &cells = std::get<0>(point_locations); 
        for (auto &cell : cells) 
          { 
            cell->set_refine_flag(); 
            for (const auto face_no : cell->face_indices()) 
              if (!cell->at_boundary(face_no)) 
                cell->neighbor(face_no)->set_refine_flag(); 
          } 
        space_grid->execute_coarsening_and_refinement(); 
      } 

// 为了构造一个良好的耦合插值算子 $C$ ，对嵌入域和被嵌入域之间的网格的相对尺寸有一些限制。耦合算子 $C$ 和空间 $V$ 和 $Q$ 必须满足一个inf-sup条件，以使问题有一个解决方案。事实证明，只要空间 $V$ 和 $Q$ 之间相互兼容（例如，只要它们被选为引言中所述的空间），非匹配 $L^2$ 投影就满足这种inf-sup条件。

// 然而，*离散*的inf-sup条件也必须成立。这里没有出现复杂的情况，但事实证明，当非匹配网格的局部直径相距太远时，离散的inf-sup常数会恶化。特别是，事实证明，如果你选择一个相对于嵌入网格*细的嵌入网格，inf-sup常数的恶化程度要比你让嵌入网格更细的情况下大得多。

// 为了避免问题，在本教程中，如果用户选择的参数使嵌入网格的最大直径大于嵌入网格的最小直径，我们将抛出一个异常。

// 这个选择保证了嵌入网格的几乎每一个单元都不超过嵌入网格的两个单元，只有一些罕见的例外，这些例外在结果的inf-sup方面可以忽略不计。

    const double embedded_space_maximal_diameter = 
      GridTools::maximal_cell_diameter(*embedded_grid, *embedded_mapping); 
    double embedding_space_minimal_diameter = 
      GridTools::minimal_cell_diameter(*space_grid); 

    deallog << "Embedding minimal diameter: " 
            << embedding_space_minimal_diameter 
            << ", embedded maximal diameter: " 
            << embedded_space_maximal_diameter << ", ratio: " 
            << embedded_space_maximal_diameter / 
                 embedding_space_minimal_diameter 
            << std::endl; 

    AssertThrow(embedded_space_maximal_diameter < 
                  embedding_space_minimal_diameter, 
                ExcMessage( 
                  "The embedding grid is too refined (or the embedded grid " 
                  "is too coarse). Adjust the parameters so that the minimal " 
                  "grid size of the embedding grid is larger " 
                  "than the maximal grid size of the embedded grid.")); 
// $\Omega$ 已经被完善，我们现在可以设置它的DoF了

    setup_embedding_dofs(); 
  } 

// 我们现在设置 $\Omega$ 和 $\Gamma$ 的DoF：因为它们基本上是独立的（除了 $\Omega$ 的网格在 $\Gamma$ 周围更加精细），这个过程是标准的。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::setup_embedding_dofs() 
  { 
    space_dh = std::make_unique<DoFHandler<spacedim>>(*space_grid); 
    space_fe = std::make_unique<FE_Q<spacedim>>( 
      parameters.embedding_space_finite_element_degree); 
    space_dh->distribute_dofs(*space_fe); 

    DoFTools::make_hanging_node_constraints(*space_dh, constraints); 
    for (auto id : parameters.homogeneous_dirichlet_ids) 
      { 
        VectorTools::interpolate_boundary_values( 
          *space_dh, id, Functions::ZeroFunction<spacedim>(), constraints); 
      } 
    constraints.close(); 

// 根据定义，刚度矩阵只涉及 $\Omega$ 的DoF。

    DynamicSparsityPattern dsp(space_dh->n_dofs(), space_dh->n_dofs()); 
    DoFTools::make_sparsity_pattern(*space_dh, dsp, constraints); 
    stiffness_sparsity.copy_from(dsp); 
    stiffness_matrix.reinit(stiffness_sparsity); 
    solution.reinit(space_dh->n_dofs()); 
    rhs.reinit(space_dh->n_dofs()); 

    deallog << "Embedding dofs: " << space_dh->n_dofs() << std::endl; 
  } 

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::setup_embedded_dofs() 
  { 
    embedded_dh = std::make_unique<DoFHandler<dim, spacedim>>(*embedded_grid); 
    embedded_fe = std::make_unique<FE_Q<dim, spacedim>>( 
      parameters.embedded_space_finite_element_degree); 
    embedded_dh->distribute_dofs(*embedded_fe); 

// 根据定义，我们要解决的系统的rhs只涉及一个零向量和 $G$ ，它只用 $\Gamma$ 的DoF计算。

    lambda.reinit(embedded_dh->n_dofs()); 
    embedded_rhs.reinit(embedded_dh->n_dofs()); 
    embedded_value.reinit(embedded_dh->n_dofs()); 

    deallog << "Embedded dofs: " << embedded_dh->n_dofs() << std::endl; 
  } 

// 创建耦合稀疏模式是一个复杂的操作，但可以使用 NonMatching::create_coupling_sparsity_pattern, 轻松完成，它需要两个DoFHandler对象、耦合的正交点、一个DynamicSparsityPattern（然后需要像往常一样复制到稀疏模式中）、嵌入和嵌入三角形的组件掩码（我们留空）以及嵌入和嵌入三角形的映射关系。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::setup_coupling() 
  { 
    TimerOutput::Scope timer_section(monitor, "Setup coupling"); 

    QGauss<dim> quad(parameters.coupling_quadrature_order); 

    DynamicSparsityPattern dsp(space_dh->n_dofs(), embedded_dh->n_dofs()); 

    NonMatching::create_coupling_sparsity_pattern(*space_grid_tools_cache, 
                                                  *space_dh, 
                                                  *embedded_dh, 
                                                  quad, 
                                                  dsp, 
                                                  AffineConstraints<double>(), 
                                                  ComponentMask(), 
                                                  ComponentMask(), 
                                                  *embedded_mapping); 
    coupling_sparsity.copy_from(dsp); 
    coupling_matrix.reinit(coupling_sparsity); 
  } 
// @sect3{Assembly}  

// 以下是创建矩阵的函数：如前所述，计算刚度矩阵和rhs是一个标准程序。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::assemble_system() 
  { 
    { 
      TimerOutput::Scope timer_section(monitor, "Assemble system"); 

//嵌入刚度矩阵  $K$  ，以及右手边  $G$  。

      MatrixTools::create_laplace_matrix( 
        *space_dh, 
        QGauss<spacedim>(2 * space_fe->degree + 1), 
        stiffness_matrix, 
        static_cast<const Function<spacedim> *>(nullptr), 
        constraints); 

      VectorTools::create_right_hand_side(*embedded_mapping, 
                                          *embedded_dh, 
                                          QGauss<dim>(2 * embedded_fe->degree + 
                                                      1), 
                                          embedded_value_function, 
                                          embedded_rhs); 
    } 
    { 
      TimerOutput::Scope timer_section(monitor, "Assemble coupling system"); 

// 为了计算耦合矩阵，我们使用 NonMatching::create_coupling_mass_matrix 工具，其工作原理与 NonMatching::create_coupling_sparsity_pattern. 类似。
      QGauss<dim> quad(parameters.coupling_quadrature_order); 
      NonMatching::create_coupling_mass_matrix(*space_grid_tools_cache, 
                                               *space_dh, 
                                               *embedded_dh, 
                                               quad, 
                                               coupling_matrix, 
                                               AffineConstraints<double>(), 
                                               ComponentMask(), 
                                               ComponentMask(), 
                                               *embedded_mapping); 

      VectorTools::interpolate(*embedded_mapping, 
                               *embedded_dh, 
                               embedded_value_function, 
                               embedded_value); 
    } 
  } 
// @sect3{Solve}  

// 所有的部分都已经组装好了：我们用舒尔补数法解决这个系统

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::solve() 
  { 
    TimerOutput::Scope timer_section(monitor, "Solve system"); 

// 从创建反刚度矩阵开始

    SparseDirectUMFPACK K_inv_umfpack; 
    K_inv_umfpack.initialize(stiffness_matrix); 

//初始化运算符，如介绍中所述

    auto K  = linear_operator(stiffness_matrix); 
    auto Ct = linear_operator(coupling_matrix); 
    auto C  = transpose_operator(Ct); 

    auto K_inv = linear_operator(K, K_inv_umfpack); 

// 使用舒尔补数法

    auto                     S = C * K_inv * Ct; 
    SolverCG<Vector<double>> solver_cg(schur_solver_control); 
    auto S_inv = inverse_operator(S, solver_cg, PreconditionIdentity()); 

    lambda = S_inv * embedded_rhs; 

    solution = K_inv * Ct * lambda; 

    constraints.distribute(solution); 
  } 

// 下面的函数只是在两个独立的文件上生成标准结果输出，每个网格一个。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::output_results() 
  { 
    TimerOutput::Scope timer_section(monitor, "Output results"); 

    DataOut<spacedim> embedding_out; 

    std::ofstream embedding_out_file("embedding.vtu"); 

    embedding_out.attach_dof_handler(*space_dh); 
    embedding_out.add_data_vector(solution, "solution"); 
    embedding_out.build_patches( 
      parameters.embedding_space_finite_element_degree); 
    embedding_out.write_vtu(embedding_out_file); 

// 这两个输出例程之间的唯一区别是，在第二种情况下，我们想在当前配置上输出数据，而不是在参考配置上。这可以通过将实际的embedded_mapping传递给 DataOut::build_patches 函数来实现。该映射将负责在实际变形的配置上输出结果。

    DataOut<dim, spacedim> embedded_out; 

    std::ofstream embedded_out_file("embedded.vtu"); 

    embedded_out.attach_dof_handler(*embedded_dh); 
    embedded_out.add_data_vector(lambda, "lambda"); 
    embedded_out.add_data_vector(embedded_value, "g"); 
    embedded_out.build_patches(*embedded_mapping, 
                               parameters.embedded_space_finite_element_degree); 
    embedded_out.write_vtu(embedded_out_file); 
  } 

// 与所有其他教程程序类似，`run()`函数只是按照正确的顺序调用所有其他方法。没有什么特别需要注意的，只是在我们实际尝试运行我们的程序之前，我们检查是否完成了解析。

  template <int dim, int spacedim> 
  void DistributedLagrangeProblem<dim, spacedim>::run() 
  { 
    AssertThrow(parameters.initialized, ExcNotInitialized()); 
    deallog.depth_console(parameters.verbosity_level); 

    setup_grids_and_dofs(); 
    setup_coupling(); 
    assemble_system(); 
    solve(); 
    output_results(); 
  } 
} // namespace Step60 

int main(int argc, char **argv) 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step60; 

      const unsigned int dim = 1, spacedim = 2; 

// 与其他教程程序中的情况不同，这里我们使用ParameterAcceptor风格的初始化，即首先构建所有对象，然后对静态方法 ParameterAcceptor::initialize 发出一次调用，以填充从ParameterAcceptor派生的类的所有参数。

// 我们检查用户是否在程序启动时指定了一个要使用的参数文件名。如果是，就尝试读取该参数文件，否则就尝试读取文件 "parameters.prm"。

// 如果指定的参数文件（隐式或显式）不存在， ParameterAcceptor::initialize 将为你创建一个，并退出程序。

      DistributedLagrangeProblem<dim, spacedim>::Parameters parameters; 
      DistributedLagrangeProblem<dim, spacedim>             problem(parameters); 

      std::string parameter_file; 
      if (argc > 1) 
        parameter_file = argv[1]; 
      else 
        parameter_file = "parameters.prm"; 

      ParameterAcceptor::initialize(parameter_file, "used_parameters.prm"); 
      problem.run(); 
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


