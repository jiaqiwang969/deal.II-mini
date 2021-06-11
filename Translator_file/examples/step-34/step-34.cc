CCTest_file/step-34.cc

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
 * Author: Luca Heltai, Cataldo Manigrasso, 2009 
 */ 


// @sect3{Include files}  

// 程序一开始就包括了一堆include文件，我们将在程序的各个部分使用这些文件。其中大部分在以前的教程中已经讨论过了。

#include <deal.II/base/smartpointer.h> 
#include <deal.II/base/convergence_table.h> 
#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/quadrature_selector.h> 
#include <deal.II/base/parsed_function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/solver_control.h> 
#include <deal.II/lac/solver_gmres.h> 
#include <deal.II/lac/precondition.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_in.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/manifold_lib.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

// 这里有一些我们需要的C++标准头文件。

#include <cmath> 
#include <iostream> 
#include <fstream> 
#include <string> 

// 这个序言的最后部分是将dealii命名空间中的所有内容导入到这个程序中的所有内容中。

namespace Step34 
{ 
  using namespace dealii; 
// @sect3{Single and double layer operator kernels}  

// 首先，让我们定义一下边界积分方程的机制。

// 以下两个函数是单层和双层势能核的实际计算，即  $G$  和  $\nabla G$  。只有当矢量 $R = \mathbf{y}-\mathbf{x}$ 不同于零时，它们才是定义良好的。

  namespace LaplaceKernel 
  { 
    template <int dim> 
    double single_layer(const Tensor<1, dim> &R) 
    { 
      switch (dim) 
        { 
          case 2: 
            return (-std::log(R.norm()) / (2 * numbers::PI)); 

          case 3: 
            return (1. / (R.norm() * 4 * numbers::PI)); 

          default: 
            Assert(false, ExcInternalError()); 
            return 0.; 
        } 
    } 

    template <int dim> 
    Tensor<1, dim> double_layer(const Tensor<1, dim> &R) 
    { 
      switch (dim) 
        { 
          case 2: 
            return R / (-2 * numbers::PI * R.norm_square()); 
          case 3: 
            return R / (-4 * numbers::PI * R.norm_square() * R.norm()); 

          default: 
            Assert(false, ExcInternalError()); 
            return Tensor<1, dim>(); 
        } 
    } 
  } // namespace LaplaceKernel 
// @sect3{The BEMProblem class}  

// 边界元素方法代码的结构与有限元素代码的结构非常相似，所以这个类的成员函数与其他大多数教程程序的成员函数一样。特别是，现在你应该熟悉从外部文件中读取参数，以及将不同的任务分割成不同的模块。这同样适用于边界元素方法，我们不会对其进行过多的评论，只是对其中的差异进行评论。

  template <int dim> 
  class BEMProblem 
  { 
  public: 
    BEMProblem(const unsigned int fe_degree      = 1, 
               const unsigned int mapping_degree = 1); 

    void run(); 

  private: 
    void read_parameters(const std::string &filename); 

    void read_domain(); 

    void refine_and_resize(); 

// 我们在这里发现的唯一真正不同的函数是装配程序。我们以最通用的方式编写了这个函数，以便能够方便地推广到高阶方法和不同的基本解（例如斯托克斯或麦克斯韦）。

// 最明显的区别是，最终的矩阵是完整的，而且我们在通常的单元格循环内有一个嵌套的循环，访问所有自由度的支持点。 此外，当支持点位于我们所访问的单元内时，我们所执行的积分就会变成单数。

// 实际的结果是，我们有两套正交公式、有限元值和临时存储，一套用于标准积分，另一套用于奇异积分，在必要时使用。

    void assemble_system(); 

// 对于这个问题的解决有两种选择。第一个是使用直接求解器，第二个是使用迭代求解器。我们选择了第二种方案。

// 我们组装的矩阵不是对称的，我们选择使用GMRES方法；然而为边界元素方法构建一个有效的预处理程序并不是一个简单的问题。这里我们使用一个非预处理的GMRES求解器。迭代求解器的选项，如公差、最大迭代次数等，都是通过参数文件选择的。

    void solve_system(); 

// 一旦我们得到了解决方案，我们将计算计算出的势的 $L^2$ 误差，以及实体角的近似值的 $L^\infty$ 误差。我们使用的网格是平滑曲线的近似值，因此计算出的角的分量或实体角的对角线矩阵  $\alpha(\mathbf{x})$  应该一直等于  $\frac 12$  。在这个例程中，我们输出势的误差和计算角度的近似值的误差。注意，后者的误差实际上不是计算角度的误差，而是衡量我们对球体和圆的近似程度。

// 对角度的计算做一些实验，对于较简单的几何形状，可以得到非常准确的结果。为了验证这一点，你可以在read_domain()方法中注释掉tria.set_manifold(1, manifold)一行，并检查程序生成的alpha。通过删除这个调用，每当细化网格时，新的节点将沿着构成粗略网格的直线放置，而不是被拉到我们真正想要近似的表面。在三维案例中，球体的粗网格是从一个立方体开始得到的，得到的字母值正好是面的节点上的 $\frac 12$ ，边的节点上的 $\frac 34$ 和顶点的8个节点上的 $\frac 78$ 。

    void compute_errors(const unsigned int cycle); 

// 一旦我们在一维领域得到了一个解决方案，我们就想把它插值到空间的其他部分。这可以通过在compute_exterior_solution()函数中再次进行解与核的卷积来实现。

// 我们想绘制速度变量，也就是势解的梯度。势解只在边界上是已知的，但我们使用与基本解的卷积在标准的二维连续有限元空间上进行插值。外推解的梯度图将给我们提供我们想要的速度。

// 除了外域上的解，我们还在output_results()函数中输出域的边界上的解，当然了。

    void compute_exterior_solution(); 

    void output_results(const unsigned int cycle); 

// 为了实现不受维度限制的编程，我们对这个单一的函数进行了专业化处理，以提取整合单元内部的奇异核所需的奇异正交公式。

    const Quadrature<dim - 1> &get_singular_quadrature( 
      const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, 
      const unsigned int index) const; 

// 通常的deal.II类可以通过指定问题的 "二维 "来用于边界元素方法。这是通过将Triangulation, FiniteElement和DoFHandler的可选第二模板参数设置为嵌入空间的维度来实现的。在我们的例子中，我们生成了嵌入在二维或三维空间的一维或二维网格。

// 可选参数默认等于第一个参数，并产生我们在之前所有例子中看到的通常的有限元类。

// 该类的构造方式是允许任意的域（通过高阶映射）和有限元空间的逼近顺序。有限元空间和映射的顺序可以在该类的构造函数中选择。

    Triangulation<dim - 1, dim> tria; 
    FE_Q<dim - 1, dim>          fe; 
    DoFHandler<dim - 1, dim>    dof_handler; 
    MappingQ<dim - 1, dim>      mapping; 

// 在BEM方法中，生成的矩阵是密集的。根据问题的大小，最终的系统可能通过直接的LU分解来解决，或者通过迭代方法来解决。在这个例子中，我们使用了一个无条件的GMRES方法。为BEM方法建立一个预处理程序是不容易的，我们在此不做处理。

    FullMatrix<double> system_matrix; 
    Vector<double>     system_rhs; 

// 接下来的两个变量将表示解决方案 $\phi$ 以及一个向量，它将保存 $\alpha(\mathbf x)$ 的值（从一个点 $\mathbf x$ 可见的 $\Omega$ 的部分）在我们形状函数的支持点。

    Vector<double> phi; 
    Vector<double> alpha; 

// 收敛表是用来输出精确解和计算出的字母的误差的。

    ConvergenceTable convergence_table; 

// 下面的变量是我们通过参数文件来填充的。 本例中我们使用的新对象是 Functions::ParsedFunction 对象和QuadratureSelector对象。

//  Functions::ParsedFunction 类允许我们通过参数文件方便快捷地定义新的函数对象，自定义的定义可以非常复杂（关于所有可用的选项，见该类的文档）。

// 我们将使用QuadratureSelector类来分配正交对象，该类允许我们根据一个识别字符串和公式本身的可能程度来生成正交公式。我们用它来允许自定义选择标准积分的正交公式，并定义奇异正交规则的顺序。

// 我们还定义了几个参数，这些参数是在我们想把解决方案扩展到整个领域的情况下使用的。

    Functions::ParsedFunction<dim> wind; 
    Functions::ParsedFunction<dim> exact_solution; 

    unsigned int                         singular_quadrature_order; 
    std::shared_ptr<Quadrature<dim - 1>> quadrature; 

    SolverControl solver_control; 

    unsigned int n_cycles; 
    unsigned int external_refinement; 

    bool run_in_this_dimension; 
    bool extend_solution; 
  }; 
// @sect4{BEMProblem::BEMProblem and BEMProblem::read_parameters}  

//构造函数初始化各种对象的方式与有限元程序（如  step-4  或  step-6  ）中的方式基本相同。这里唯一的新成分是ParsedFunction对象，它在构造时需要说明组件的数量。

// 对于精确解来说，向量分量的数量是1，而且不需要任何操作，因为1是ParsedFunction对象的默认值。然而，风需要指定dim组件。注意，在为 Functions::ParsedFunction, 的表达式声明参数文件中的条目时，我们需要明确指定分量的数量，因为函数 Functions::ParsedFunction::declare_parameters 是静态的，对分量的数量没有了解。

  template <int dim> 
  BEMProblem<dim>::BEMProblem(const unsigned int fe_degree, 
                              const unsigned int mapping_degree) 
    : fe(fe_degree) 
    , dof_handler(tria) 
    , mapping(mapping_degree, true) 
    , wind(dim) 
    , singular_quadrature_order(5) 
    , n_cycles(4) 
    , external_refinement(5) 
    , run_in_this_dimension(true) 
    , extend_solution(true) 
  {} 

  template <int dim> 
  void BEMProblem<dim>::read_parameters(const std::string &filename) 
  { 
    deallog << std::endl 
            << "Parsing parameter file " << filename << std::endl 
            << "for a " << dim << " dimensional simulation. " << std::endl; 

    ParameterHandler prm; 

    prm.declare_entry("Number of cycles", "4", Patterns::Integer()); 
    prm.declare_entry("External refinement", "5", Patterns::Integer()); 
    prm.declare_entry("Extend solution on the -2,2 box", 
                      "true", 
                      Patterns::Bool()); 
    prm.declare_entry("Run 2d simulation", "true", Patterns::Bool()); 
    prm.declare_entry("Run 3d simulation", "true", Patterns::Bool()); 

    prm.enter_subsection("Quadrature rules"); 
    { 
      prm.declare_entry( 
        "Quadrature type", 
        "gauss", 
        Patterns::Selection( 
          QuadratureSelector<(dim - 1)>::get_quadrature_names())); 
      prm.declare_entry("Quadrature order", "4", Patterns::Integer()); 
      prm.declare_entry("Singular quadrature order", "5", Patterns::Integer()); 
    } 
    prm.leave_subsection(); 

// 对于二维和三维，我们将默认的输入数据设置为：解为  $x+y$  或  $x+y+z$  。实际计算出的解在无穷大时的数值为零。在这种情况下，这与精确解相吻合，不需要额外的修正，但是你应该注意，我们任意设置了 $\phi_\infty$ ，而我们传递给程序的精确解需要在无穷远处有相同的值，才能正确计算出误差。

//  Functions::ParsedFunction 对象的使用是非常直接的。 Functions::ParsedFunction::declare_parameters 函数需要一个额外的整数参数，指定给定函数的分量数量。它的默认值是1。当相应的 Functions::ParsedFunction::parse_parameters 方法被调用时，调用对象必须有与这里定义的相同数量的组件，否则会产生异常。

// 在声明条目时，我们同时声明了二维和三维的函数。然而只有二维的最终被解析。这使得我们对二维和三维问题都只需要一个参数文件。

// 注意，从数学的角度来看，边界上的风函数应该满足条件 $\int_{\partial\Omega} \mathbf{v}\cdot \mathbf{n} d \Gamma = 0$  ，这样问题才会有解。如果不满足这个条件，那么就找不到解，求解器也就不会收敛。

    prm.enter_subsection("Wind function 2d"); 
    { 
      Functions::ParsedFunction<2>::declare_parameters(prm, 2); 
      prm.set("Function expression", "1; 1"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Wind function 3d"); 
    { 
      Functions::ParsedFunction<3>::declare_parameters(prm, 3); 
      prm.set("Function expression", "1; 1; 1"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution 2d"); 
    { 
      Functions::ParsedFunction<2>::declare_parameters(prm); 
      prm.set("Function expression", "x+y"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution 3d"); 
    { 
      Functions::ParsedFunction<3>::declare_parameters(prm); 
      prm.set("Function expression", "x+y+z"); 
    } 
    prm.leave_subsection(); 

// 在求解器部分，我们设置所有的SolverControl参数。然后，该对象将在solve_system()函数中被送入GMRES求解器。

    prm.enter_subsection("Solver"); 
    SolverControl::declare_parameters(prm); 
    prm.leave_subsection(); 

// 在向ParameterHandler对象声明了所有这些参数后，让我们读取一个输入文件，该文件将为这些参数提供其值。然后我们继续从ParameterHandler对象中提取这些值。

    prm.parse_input(filename); 

    n_cycles            = prm.get_integer("Number of cycles"); 
    external_refinement = prm.get_integer("External refinement"); 
    extend_solution     = prm.get_bool("Extend solution on the -2,2 box"); 

    prm.enter_subsection("Quadrature rules"); 
    { 
      quadrature = std::shared_ptr<Quadrature<dim - 1>>( 
        new QuadratureSelector<dim - 1>(prm.get("Quadrature type"), 
                                        prm.get_integer("Quadrature order"))); 
      singular_quadrature_order = prm.get_integer("Singular quadrature order"); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Wind function " + std::to_string(dim) + "d"); 
    { 
      wind.parse_parameters(prm); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Exact solution " + std::to_string(dim) + "d"); 
    { 
      exact_solution.parse_parameters(prm); 
    } 
    prm.leave_subsection(); 

    prm.enter_subsection("Solver"); 
    solver_control.parse_parameters(prm); 
    prm.leave_subsection(); 

// 最后，这里是另一个如何在独立维度编程中使用参数文件的例子。 如果我们想关闭两个模拟中的一个，我们可以通过设置相应的 "运行2D模拟 "或 "运行3D模拟 "标志为假来实现。

    run_in_this_dimension = 
      prm.get_bool("Run " + std::to_string(dim) + "d simulation"); 
  } 
// @sect4{BEMProblem::read_domain}  

// 边界元素法三角剖分与（dim-1）维三角剖分基本相同，不同之处在于顶点属于（dim）维空间。

// deal.II中支持的一些网格格式默认使用三维点来描述网格。这些格式与deal.II的边界元素方法功能兼容。特别是我们可以使用UCD或GMSH格式。在这两种情况下，我们必须特别注意网格的方向，因为与标准有限元的情况不同，这里没有进行重新排序或兼容性检查。 所有的网格都被认为是有方向性的，因为它们被嵌入到一个高维空间中。参见GridIn和Triangulation的文档，以进一步了解三角结构中单元的方向。在我们的例子中，网格的法线是外在于2D的圆或3D的球体。

// 对边界元素网格进行适当细化所需要的另一个细节是对网格所逼近的流形的准确描述。对于标准有限元网格的边界，我们已经多次看到了这一点（例如在 step-5 和 step-6 中），这里的原理和用法是一样的，只是SphericalManifold类需要一个额外的模板参数来指定嵌入空间维度。

  template <int dim> 
  void BEMProblem<dim>::read_domain() 
  { 
    const Point<dim>                      center = Point<dim>(); 
    const SphericalManifold<dim - 1, dim> manifold(center); 

    std::ifstream in; 
    switch (dim) 
      { 
        case 2: 
          in.open("coarse_circle.inp"); 
          break; 

        case 3: 
          in.open("coarse_sphere.inp"); 
          break; 

        default: 
          Assert(false, ExcNotImplemented()); 
      } 

    GridIn<dim - 1, dim> gi; 
    gi.attach_triangulation(tria); 
    gi.read_ucd(in); 

    tria.set_all_manifold_ids(1); 

// 对  Triangulation::set_manifold  的调用复制了流形（通过  Manifold::clone()),  所以我们不需要担心对  <code>manifold</code>  的无效指针。

    tria.set_manifold(1, manifold); 
  } 
// @sect4{BEMProblem::refine_and_resize}  

// 这个函数对网格进行全局细化，分配自由度，并调整矩阵和向量的大小。

  template <int dim> 
  void BEMProblem<dim>::refine_and_resize() 
  { 
    tria.refine_global(1); 

    dof_handler.distribute_dofs(fe); 

    const unsigned int n_dofs = dof_handler.n_dofs(); 

    system_matrix.reinit(n_dofs, n_dofs); 

    system_rhs.reinit(n_dofs); 
    phi.reinit(n_dofs); 
    alpha.reinit(n_dofs); 
  } 
// @sect4{BEMProblem::assemble_system}  

// 下面是这个程序的主要功能，组装与边界积分方程相对应的矩阵。

  template <int dim> 
  void BEMProblem<dim>::assemble_system() 
  { 

// 首先我们用正交公式初始化一个FEValues对象，用于在非奇异单元中进行内核积分。这个正交公式是通过参数文件选择的，并且需要相当精确，因为我们要积分的函数不是多项式函数。

    FEValues<dim - 1, dim> fe_v(mapping, 
                                fe, 
                                *quadrature, 
                                update_values | update_normal_vectors | 
                                  update_quadrature_points | update_JxW_values); 

    const unsigned int n_q_points = fe_v.n_quadrature_points; 

    std::vector<types::global_dof_index> local_dof_indices( 
      fe.n_dofs_per_cell()); 

    std::vector<Vector<double>> cell_wind(n_q_points, Vector<double>(dim)); 
    double                      normal_wind; 

// 与有限元方法不同的是，如果我们使用拼合边界元方法，那么在每个装配循环中，我们只装配与一个自由度（与支撑点 $i$ 相关的自由度）和当前单元之间的耦合信息。这是用fe.dofs_per_cell元素的向量完成的，然后将其分配到全局行的矩阵中  $i$  。以下对象将持有这些信息。

    Vector<double> local_matrix_row_i(fe.n_dofs_per_cell()); 

// 索引  $i$  运行在拼合点上，这是  $i$  第三个基函数的支持点，而  $j$  运行在内部积分点上。

// 我们构建一个支持点的向量，它将用于局部积分。

    std::vector<Point<dim>> support_points(dof_handler.n_dofs()); 
    DoFTools::map_dofs_to_support_points<dim - 1, dim>(mapping, 
                                                       dof_handler, 
                                                       support_points); 

// 这样做之后，我们就可以开始对所有单元进行积分循环，首先初始化FEValues对象，得到正交点的 $\mathbf{\tilde v}$ 的值（这个向量场应该是常数，但更通用也无妨）。

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 
        cell->get_dof_indices(local_dof_indices); 

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
        const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 
        wind.vector_value_list(q_points, cell_wind); 

// 然后我们在当前单元上形成所有自由度的积分（注意，这包括不在当前单元上的自由度，这与通常的有限元积分有偏差）。如果其中一个局部自由度与支持点 $i$ 相同，我们需要执行的积分是单数。因此，在循环的开始，我们检查是否是这种情况，并存储哪一个是奇异指数。

        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
          { 
            local_matrix_row_i = 0; 

            bool         is_singular    = false; 
            unsigned int singular_index = numbers::invalid_unsigned_int; 

            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              if (local_dof_indices[j] == i) 
                { 
                  singular_index = j; 
                  is_singular    = true; 
                  break; 
                } 

// 然后我们进行积分。如果指数 $i$ 不是局部自由度之一，我们只需将单层项加到右边，将双层项加到矩阵中。

            if (is_singular == false) 
              { 
                for (unsigned int q = 0; q < n_q_points; ++q) 
                  { 
                    normal_wind = 0; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal_wind += normals[q][d] * cell_wind[q](d); 

                    const Tensor<1, dim> R = q_points[q] - support_points[i]; 

                    system_rhs(i) += (LaplaceKernel::single_layer(R) * 
                                      normal_wind * fe_v.JxW(q)); 

                    for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 

                      local_matrix_row_i(j) -= 
                        ((LaplaceKernel::double_layer(R) * normals[q]) * 
                         fe_v.shape_value(j, q) * fe_v.JxW(q)); 
                  } 
              } 
            else 
              { 

// 现在我们处理更微妙的情况。如果我们在这里，这意味着在 $j$ 索引上运行的单元包含support_point[i]。在这种情况下，单层和双层势都是单数，它们需要特殊处理。            
//每当在给定单元内进行积分时，就会使用一个特殊的正交公式，允许人们对参考单元上的奇异权重进行任意函数的积分。            
//正确的正交公式由get_singular_quadrature函数选择，下面将详细说明。

                Assert(singular_index != numbers::invalid_unsigned_int, 
                       ExcInternalError()); 

                const Quadrature<dim - 1> &singular_quadrature = 
                  get_singular_quadrature(cell, singular_index); 

                FEValues<dim - 1, dim> fe_v_singular( 
                  mapping, 
                  fe, 
                  singular_quadrature, 
                  update_jacobians | update_values | update_normal_vectors | 
                    update_quadrature_points); 

                fe_v_singular.reinit(cell); 

                std::vector<Vector<double>> singular_cell_wind( 
                  singular_quadrature.size(), Vector<double>(dim)); 

                const std::vector<Tensor<1, dim>> &singular_normals = 
                  fe_v_singular.get_normal_vectors(); 
                const std::vector<Point<dim>> &singular_q_points = 
                  fe_v_singular.get_quadrature_points(); 

                wind.vector_value_list(singular_q_points, singular_cell_wind); 

                for (unsigned int q = 0; q < singular_quadrature.size(); ++q) 
                  { 
                    const Tensor<1, dim> R = 
                      singular_q_points[q] - support_points[i]; 
                    double normal_wind = 0; 
                    for (unsigned int d = 0; d < dim; ++d) 
                      normal_wind += 
                        (singular_cell_wind[q](d) * singular_normals[q][d]); 

                    system_rhs(i) += (LaplaceKernel::single_layer(R) * 
                                      normal_wind * fe_v_singular.JxW(q)); 

                    for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
                      { 
                        local_matrix_row_i(j) -= 
                          ((LaplaceKernel::double_layer(R) * 
                            singular_normals[q]) * 
                           fe_v_singular.shape_value(j, q) * 
                           fe_v_singular.JxW(q)); 
                      } 
                  } 
              } 

// 最后，我们需要将当前单元格的贡献添加到全局矩阵中。

            for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
              system_matrix(i, local_dof_indices[j]) += local_matrix_row_i(j); 
          } 
      } 

// 积分运算符的第二部分是术语  $\alpha(\mathbf{x}_i) \phi_j(\mathbf{x}_i)$  。由于我们使用的是配位方案， $\phi_j(\mathbf{x}_i)=\delta_{ij}$  而相应的矩阵是一个对角线的矩阵，其条目等于 $\alpha(\mathbf{x}_i)$  。

// 计算这个实体角的对角矩阵的一个快速方法是使用诺伊曼矩阵本身。只需将该矩阵与一个元素都等于-1的向量相乘，就可以得到阿尔法角或实体角的对角线矩阵（见介绍中的公式）。然后将这个结果加回到系统矩阵对象上，得到矩阵的最终形式。

    Vector<double> ones(dof_handler.n_dofs()); 
    ones.add(-1.); 

    system_matrix.vmult(alpha, ones); 
    alpha.add(1); 
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
      system_matrix(i, i) += alpha(i); 
  } 
// @sect4{BEMProblem::solve_system}  

// 下一个函数简单地解决了线性系统。

  template <int dim> 
  void BEMProblem<dim>::solve_system() 
  { 
    SolverGMRES<Vector<double>> solver(solver_control); 
    solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity()); 
  } 
// @sect4{BEMProblem::compute_errors}  

// 误差的计算在其他所有的例子程序中都是完全一样的，我们就不做过多的评论。请注意，在有限元方法中使用的方法在这里也可以使用。

  template <int dim> 
  void BEMProblem<dim>::compute_errors(const unsigned int cycle) 
  { 
    Vector<float> difference_per_cell(tria.n_active_cells()); 
    VectorTools::integrate_difference(mapping, 
                                      dof_handler, 
                                      phi, 
                                      exact_solution, 
                                      difference_per_cell, 
                                      QGauss<(dim - 1)>(2 * fe.degree + 1), 
                                      VectorTools::L2_norm); 
    const double L2_error = 
      VectorTools::compute_global_error(tria, 
                                        difference_per_cell, 
                                        VectorTools::L2_norm); 

//可以直接使用 Vector::linfty_norm() 函数来计算α向量的误差，因为在每个节点上，该值应该是 $\frac 12$  。然后，所有的误差都会被输出并附加到我们的ConvergenceTable对象中，以便以后计算收敛率。

    Vector<double> difference_per_node(alpha); 
    difference_per_node.add(-.5); 

    const double       alpha_error    = difference_per_node.linfty_norm(); 
    const unsigned int n_active_cells = tria.n_active_cells(); 
    const unsigned int n_dofs         = dof_handler.n_dofs(); 

    deallog << "Cycle " << cycle << ':' << std::endl 
            << "   Number of active cells:       " << n_active_cells 
            << std::endl 
            << "   Number of degrees of freedom: " << n_dofs << std::endl; 

    convergence_table.add_value("cycle", cycle); 
    convergence_table.add_value("cells", n_active_cells); 
    convergence_table.add_value("dofs", n_dofs); 
    convergence_table.add_value("L2(phi)", L2_error); 
    convergence_table.add_value("Linfty(alpha)", alpha_error); 
  } 

// 奇异积分需要仔细选择正交规则。特别是deal.II库提供了为对数奇异性（QGaussLog, QGaussLogR）以及1/R奇异性（QGaussOneOverR）量身定制的正交规则。

// 奇异积分通常是通过构建具有奇异权重的加权正交公式得到的，因此可以写成

// \f[ \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i) \f]

// 其中 $s(x)$ 是一个给定的奇点，权重和正交点 $w_i,q_i$ 是精心选择的，以使上述公式对某类函数 $f(x)$ 是一个等式。

// 在我们目前看到的所有有限元例子中，正交点本身的权重（即函数  $s(x)$  ），总是不断地等于1。 对于奇异积分，我们有两个选择：我们可以使用上面的定义，从积分中剔除奇异性（即用特殊的正交规则对 $f(x)$ 进行积分），或者我们可以要求正交规则用 $s(q_i)$ 对权重 $w_i$ 进行 "标准化"。

// \f[ \int_K f(x) s(x) dx = \int_K g(x) dx = \sum_{i=1}^N \frac{w_i}{s(q_i)} g(q_i) \f]

// 我们通过QGaussLogR和QGaussOneOverR的 @p factor_out_singularity 参数来使用这第二种选择。

// 这些积分有些微妙，特别是在二维空间，由于从实数到参考单元的转换，积分的变量是以转换的行列式为尺度的。

// 在二维空间中，这个过程不仅会导致一个因子作为常数出现在整个积分上，而且还会导致一个需要评估的额外积分。

// \f[ \int_0^1 f(x)\ln(x/\alpha) dx = \int_0^1 f(x)\ln(x) dx - \int_0^1  f(x) \ln(\alpha) dx.  \f]

// 这个过程由QGaussLogR类的构造函数来处理，它增加了额外的正交点和权重，以考虑到积分的第二部分。

// 类似的推理应该在三维情况下进行，因为奇异正交是在参考单元的半径 $r$ 的逆上定制的，而我们的奇异函数生活在实空间，然而在三维情况下一切都更简单，因为奇异性与变换的行列式成线性比例。这使我们可以只建立一次奇异的二维正交规则，并在所有单元格中重复使用。

// 在一维的奇异积分中，这是不可能的，因为我们需要知道正交的缩放参数，而这个参数并不是先验的。这里，正交规则本身也取决于当前单元格的大小。出于这个原因，有必要为每个单数积分创建一个新的正交。

// 不同的正交规则是在get_singular_quadrature中建立的，它专门用于dim=2和dim=3，它们在assemble_system函数中被检索。作为参数给出的索引是奇异点所在的单位支持点的索引。

  template <> 
  const Quadrature<2> &BEMProblem<3>::get_singular_quadrature( 
    const DoFHandler<2, 3>::active_cell_iterator &, 
    const unsigned int index) const 
  { 
    Assert(index < fe.n_dofs_per_cell(), 
           ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 

    static std::vector<QGaussOneOverR<2>> quadratures; 
    if (quadratures.size() == 0) 
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
        quadratures.emplace_back(singular_quadrature_order, 
                                 fe.get_unit_support_points()[i], 
                                 true); 
    return quadratures[index]; 
  } 

  template <> 
  const Quadrature<1> &BEMProblem<2>::get_singular_quadrature( 
    const DoFHandler<1, 2>::active_cell_iterator &cell, 
    const unsigned int                            index) const 
  { 
    Assert(index < fe.n_dofs_per_cell(), 
           ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 

    static Quadrature<1> *q_pointer = nullptr; 
    if (q_pointer) 
      delete q_pointer; 

    q_pointer = new QGaussLogR<1>(singular_quadrature_order, 
                                  fe.get_unit_support_points()[index], 
                                  1. / cell->measure(), 
                                  true); 
    return (*q_pointer); 
  } 

//  @sect4{BEMProblem::compute_exterior_solution}  

// 我们还想知道一些关于外域中电势 $\phi$ 的值：毕竟我们考虑边界积分问题的动机是想知道外域中的速度!

// 为此，我们在此假设边界元素域包含在盒子 $[-2,2]^{\text{dim}}$ 中，我们用与基本解的卷积来推算这个盒子内的实际解。这方面的公式在介绍中已经给出。

// 整个空间的解的重建是在一个连续的、尺寸为dim的有限元网格上完成的。这些都是常用的，我们不做进一步评论。在函数的最后，我们再次以通常的方式输出这个外部解。

  template <int dim> 
  void BEMProblem<dim>::compute_exterior_solution() 
  { 
    Triangulation<dim> external_tria; 
    GridGenerator::hyper_cube(external_tria, -2, 2); 

    FE_Q<dim>       external_fe(1); 
    DoFHandler<dim> external_dh(external_tria); 
    Vector<double>  external_phi; 

    external_tria.refine_global(external_refinement); 
    external_dh.distribute_dofs(external_fe); 
    external_phi.reinit(external_dh.n_dofs()); 

    FEValues<dim - 1, dim> fe_v(mapping, 
                                fe, 
                                *quadrature, 
                                update_values | update_normal_vectors | 
                                  update_quadrature_points | update_JxW_values); 

    const unsigned int n_q_points = fe_v.n_quadrature_points; 

    std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell()); 

    std::vector<double>         local_phi(n_q_points); 
    std::vector<double>         normal_wind(n_q_points); 
    std::vector<Vector<double>> local_wind(n_q_points, Vector<double>(dim)); 

    std::vector<Point<dim>> external_support_points(external_dh.n_dofs()); 
    DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping, 
                                              external_dh, 
                                              external_support_points); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_v.reinit(cell); 

        const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
        const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 

        cell->get_dof_indices(dofs); 
        fe_v.get_function_values(phi, local_phi); 

        wind.vector_value_list(q_points, local_wind); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            normal_wind[q] = 0; 
            for (unsigned int d = 0; d < dim; ++d) 
              normal_wind[q] += normals[q][d] * local_wind[q](d); 
          } 

        for (unsigned int i = 0; i < external_dh.n_dofs(); ++i) 
          for (unsigned int q = 0; q < n_q_points; ++q) 
            { 
              const Tensor<1, dim> R = q_points[q] - external_support_points[i]; 

              external_phi(i) += 
                ((LaplaceKernel::single_layer(R) * normal_wind[q] + 
                  (LaplaceKernel::double_layer(R) * normals[q]) * 
                    local_phi[q]) * 
                 fe_v.JxW(q)); 
            } 
      } 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(external_dh); 
    data_out.add_data_vector(external_phi, "external_phi"); 
    data_out.build_patches(); 

    const std::string filename = std::to_string(dim) + "d_external.vtk"; 
    std::ofstream     file(filename); 

 
  } 
// @sect4{BEMProblem::output_results}  

// 输出我们的计算结果是一个相当机械的任务。这个函数的所有组成部分之前已经讨论过了。

  template <int dim> 
  void BEMProblem<dim>::output_results(const unsigned int cycle) 
  { 
    DataOut<dim - 1, dim> dataout; 

    dataout.attach_dof_handler(dof_handler); 
    dataout.add_data_vector(phi, "phi", DataOut<dim - 1, dim>::type_dof_data); 
    dataout.add_data_vector(alpha, 
                            "alpha", 
                            DataOut<dim - 1, dim>::type_dof_data); 
    dataout.build_patches(mapping, 
                          mapping.get_degree(), 
                          DataOut<dim - 1, dim>::curved_inner_cells); 

    const std::string filename = std::to_string(dim) + "d_boundary_solution_" + 
                                 std::to_string(cycle) + ".vtk"; 
    std::ofstream file(filename); 

    dataout.write_vtk(file); 

    if (cycle == n_cycles - 1) 
      { 
        convergence_table.set_precision("L2(phi)", 3); 
        convergence_table.set_precision("Linfty(alpha)", 3); 

        convergence_table.set_scientific("L2(phi)", true); 
        convergence_table.set_scientific("Linfty(alpha)", true); 

 
          "L2(phi)", ConvergenceTable::reduction_rate_log2); 
        convergence_table.evaluate_convergence_rates( 
          "Linfty(alpha)", ConvergenceTable::reduction_rate_log2); 
        deallog << std::endl; 
        convergence_table.write_text(std::cout); 
      } 
  } 
// @sect4{BEMProblem::run}  

// 这是最主要的功能。它应该是不言自明的。

  template <int dim> 
  void BEMProblem<dim>::run() 
  { 
    read_parameters("parameters.prm"); 

    if (run_in_this_dimension == false) 
      { 
        deallog << "Run in dimension " << dim 
                << " explicitly disabled in parameter file. " << std::endl; 
        return; 
      } 

    read_domain(); 

    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        refine_and_resize(); 
        assemble_system(); 
        solve_system(); 
        compute_errors(cycle); 
        output_results(cycle); 
      } 

    if (extend_solution == true) 
      compute_exterior_solution(); 
  } 
} // namespace Step34 
// @sect3{The main() function}  

// 这是本程序的主要功能。它与以前所有的教程程序完全一样。

int main() 
{ 
  try 
    { 
      using namespace Step34; 

      const unsigned int degree         = 1; 
      const unsigned int mapping_degree = 1; 

      deallog.depth_console(3); 
      BEMProblem<2> laplace_problem_2d(degree, mapping_degree); 
      laplace_problem_2d.run(); 

      BEMProblem<3> laplace_problem_3d(degree, mapping_degree); 
      laplace_problem_3d.run(); 
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


