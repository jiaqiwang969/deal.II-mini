

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
 * Author: Guido Kanschat, Texas A&M University, 2009 
 */ 



// 前面几个文件已经在前面的例子中讲过了，因此不再做进一步的评论。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 
#include <deal.II/lac/vector.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_refinement.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/fe/mapping_q1.h> 

// 这里定义了不连续的有限元。它们的使用方式与所有其他有限元相同，不过--正如你在以前的教程程序中所看到的--用户与有限元类的交互并不多：它们被传递给 <code>DoFHandler</code> 和 <code>FEValues</code> 对象，就这样了。

#include <deal.II/fe/fe_dgq.h> 

// 我们将使用最简单的求解器，称为Richardson迭代，它代表了一个简单的缺陷修正。这与一个块状SSOR预处理器（定义在precondition_block.h中）相结合，该预处理器使用DG离散化产生的系统矩阵的特殊块状结构。

#include <deal.II/lac/solver_richardson.h> 
#include <deal.II/lac/precondition_block.h> 

// 我们将使用梯度作为细化指标。

#include <deal.II/numerics/derivative_approximation.h> 

// 这里是使用MeshWorker框架的新的包含文件。第一个文件包含了 MeshWorker::DoFInfo, 类，它为局部积分器提供了局部与全局自由度之间的映射。在第二个文件中，我们发现一个类型为 MeshWorker::IntegrationInfo, 的对象，它主要是对一组FEValues对象的封装。文件<tt>meshworker/simple.h</tt>包含了将局部集成数据组装成只包含一个矩阵的全局系统的类。最后，我们将需要在所有的网格单元和面中运行循环的文件。

#include <deal.II/meshworker/dof_info.h> 
#include <deal.II/meshworker/integration_info.h> 
#include <deal.II/meshworker/simple.h> 
#include <deal.II/meshworker/loop.h> 

// 像所有的程序一样，我们在完成这一部分时要包括所需的C++头文件，并声明我们要使用dealii命名空间中的对象，不加前缀。

#include <iostream> 
#include <fstream> 

namespace Step12 
{ 
  using namespace dealii; 
// @sect3{Equation data}  

// 首先，我们定义一个描述不均匀边界数据的类。由于只使用它的值，我们实现value_list()，但不定义Function的所有其他函数。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() = default; 
    virtual void value_list(const std::vector<Point<dim>> &points, 
                            std::vector<double> &          values, 
                            const unsigned int component = 0) const override; 
  }; 

// 考虑到流动方向，单位方块 $[0,1]^2$ 的流入边界为右边界和下边界。我们在x轴上规定了不连续的边界值1和0，在右边界上规定了值0。该函数在流出边界上的值将不会在DG方案中使用。

  template <int dim> 
  void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points, 
                                       std::vector<double> &          values, 
                                       const unsigned int component) const 
  { 
    (void)component; 
    AssertIndexRange(component, 1); 
    Assert(values.size() == points.size(), 
           ExcDimensionMismatch(values.size(), points.size())); 

    for (unsigned int i = 0; i < values.size(); ++i) 
      { 
        if (points[i](0) < 0.5) 
          values[i] = 1.; 
        else 
          values[i] = 0.; 
      } 
  } 

// 最后，一个计算并返回风场的函数  $\beta=\beta(\mathbf x)$  。正如在介绍中所解释的，在2D中我们将使用一个围绕原点的旋转场。在3D中，我们只需不设置 $z$ 分量（即为零），而这个函数在目前的实现中不能用于1D。

  template <int dim> 
  Tensor<1, dim> beta(const Point<dim> &p) 
  { 
    Assert(dim >= 2, ExcNotImplemented()); 

    Tensor<1, dim> wind_field; 
    wind_field[0] = -p[1]; 
    wind_field[1] = p[0]; 
    wind_field /= wind_field.norm(); 

    return wind_field; 
  } 
// @sect3{The AdvectionProblem class}  

// 在这个准备工作之后，我们继续进行这个程序的主类，叫做AdvectionProblem。它基本上是  step-6  的主类。我们没有AffineConstraints对象，因为在DG离散中没有悬挂节点约束。

// 主要的区别只出现在集合函数的实现上，因为在这里，我们不仅需要覆盖面上的通量积分，我们还使用MeshWorker接口来简化涉及的循环。

  template <int dim> 
  class AdvectionProblem 
  { 
  public: 
    AdvectionProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(Vector<double> &solution); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 

// 此外，我们想使用程度为1的DG元素（但这只在构造函数中指定）。如果你想使用不同度数的DG方法，整个程序保持不变，只需在构造函数中用所需的多项式度数替换1。

    FE_DGQ<dim>     fe; 
    DoFHandler<dim> dof_handler; 

// 接下来的四个成员代表要解决的线性系统。  <code>system_matrix</code> and <code>right_hand_side</code> 是由 <code>assemble_system()</code>, the <code>solution</code> 产生的， <code>solve()</code>. The <code>sparsity_pattern</code> 是用来确定 <code>system_matrix</code> 中非零元素的位置。

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> right_hand_side; 

// 最后，我们必须提供集合单元、边界和内表面条款的函数。在MeshWorker框架中，所有单元的循环和大部分操作的设置都将在这个类之外完成，所以我们所要提供的只是这三个操作。他们将在中间对象上工作，首先，我们在这里定义了交给本地集成函数的信息对象的别名，以使我们的生活更轻松。

    using DoFInfo  = MeshWorker::DoFInfo<dim>; 
    using CellInfo = MeshWorker::IntegrationInfo<dim>; 

// 下面的三个函数是在所有单元和面的通用循环中被调用的。它们是进行实际整合的函数。

// 在我们下面的代码中，这些函数并不访问当前类的成员变量，所以我们可以将它们标记为 <code>static</code> ，并简单地将这些函数的指针传递给MeshWorker框架。然而，如果这些函数想要访问成员变量（或者需要额外的参数，而不是下面指定的参数），我们可以使用lambda函数的设施来为MeshWorker框架提供对象，这些对象就像它们拥有所需的参数数量和类型一样，但实际上已经绑定了其他参数。

    static void integrate_cell_term(DoFInfo &dinfo, CellInfo &info); 
    static void integrate_boundary_term(DoFInfo &dinfo, CellInfo &info); 
    static void integrate_face_term(DoFInfo & dinfo1, 
                                    DoFInfo & dinfo2, 
                                    CellInfo &info1, 
                                    CellInfo &info2); 
  }; 

// 我们从构造函数开始。 <code>fe</code> 的构造器调用中的1是多项式的度数。

  template <int dim> 
  AdvectionProblem<dim>::AdvectionProblem() 
    : mapping() 
    , fe(1) 
    , dof_handler(triangulation) 
  {} 

  template <int dim> 
  void AdvectionProblem<dim>::setup_system() 
  { 

// 在设置通常的有限元数据结构的函数中，我们首先需要分配DoF。

    dof_handler.distribute_dofs(fe); 

// 我们从生成稀疏模式开始。为此，我们首先用系统中出现的耦合物填充一个动态稀疏模式（DynamicSparsityPattern）类型的中间对象。在建立模式之后，这个对象被复制到 <code>sparsity_pattern</code> 并可以被丢弃。

// 为了建立DG离散的稀疏模式，我们可以调用类似于 DoFTools::make_sparsity_pattern, 的函数，它被称为 DoFTools::make_flux_sparsity_pattern:  
    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

// 最后，我们设置了线性系统的所有组成部分的结构。

    system_matrix.reinit(sparsity_pattern); 
    solution.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{The assemble_system function}  

// 这里我们看到了与手工组装的主要区别。我们不需要在单元格和面上写循环，而是将这一切交给MeshWorker框架。为了做到这一点，我们只需要定义局部的集成函数，并使用命名空间 MeshWorker::Assembler 中的一个类来构建全局系统。

  template <int dim> 
  void AdvectionProblem<dim>::assemble_system() 
  { 

// 这是一个神奇的对象，它知道关于数据结构和局部集成的一切。 这是在函数 MeshWorker::loop(), 中做工作的对象，它被下面的 MeshWorker::integration_loop() 隐式调用。在我们提供指针的函数完成局部积分后， MeshWorker::Assembler::SystemSimple 对象将这些数据分配到全局稀疏矩阵和右手边的向量。

    MeshWorker::IntegrationInfoBox<dim> info_box; 

// 首先，我们在工作者基类中初始化正交公式和更新标志。对于正交，我们采取安全措施，使用QGauss公式，其点数比使用的多项式度数高一个。由于单元格、边界和内部面的正交率可以独立选择，我们必须把这个值交给三次。

    const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1; 
    info_box.initialize_gauss_quadrature(n_gauss_points, 
                                         n_gauss_points, 
                                         n_gauss_points); 

// 这些是我们整合系统时需要的数值类型。它们被添加到单元格、边界和内部面以及内部邻居面所使用的标志中，这是由四个 @p true 值强制执行的。

    info_box.initialize_update_flags(); 
    UpdateFlags update_flags = 
      update_quadrature_points | update_values | update_gradients; 
    info_box.add_update_flags(update_flags, true, true, true, true); 

// 在准备好<tt>info_box</tt>中的所有数据后，我们初始化其中的FEValues对象。

    info_box.initialize(fe, mapping); 

// 到目前为止创建的对象帮助我们在每个单元和面进行局部积分。现在，我们需要一个对象来接收整合后的（本地）数据，并将它们转发给装配程序。

    MeshWorker::DoFInfo<dim> dof_info(dof_handler); 

// 现在，我们必须创建装配器对象，并告诉它将本地数据放在哪里。这些将是我们的系统矩阵和右手边的数据。

    MeshWorker::Assembler::SystemSimple<SparseMatrix<double>, Vector<double>> 
      assembler; 
    assembler.initialize(system_matrix, right_hand_side); 

// 最后，在所有活动单元上进行积分循环（由第一个参数决定，它是一个活动迭代器）。

// 正如在类声明中声明局部积分函数时的讨论中所指出的，装配积分器类所期望的参数实际上不是函数指针。相反，它们是可以像函数一样被调用的对象，有一定数量的参数。因此，我们也可以在这里传递具有适当的operator()实现的对象，或者如果本地集成器是，例如，非静态成员函数，则可以传递lambda函数。

    MeshWorker::loop<dim, 
                     dim, 
                     MeshWorker::DoFInfo<dim>, 
                     MeshWorker::IntegrationInfoBox<dim>>( 
      dof_handler.begin_active(), 
      dof_handler.end(), 
      dof_info, 
      info_box, 
      &AdvectionProblem<dim>::integrate_cell_term, 
      &AdvectionProblem<dim>::integrate_boundary_term, 
      &AdvectionProblem<dim>::integrate_face_term, 
      assembler); 
  } 
// @sect4{The local integrators}  

// 这些是给上面调用的 MeshWorker::integration_loop() 的函数。它们计算单元格和面中对系统矩阵和右手边的局部贡献。

  template <int dim> 
  void AdvectionProblem<dim>::integrate_cell_term(DoFInfo & dinfo, 
                                                  CellInfo &info) 
  { 

// 首先，让我们从 @p info. 中检索这里使用的一些对象。注意，这些对象可以处理更复杂的结构，因此这里的访问看起来比看起来更复杂。

    const FEValuesBase<dim> &  fe_values    = info.fe_values(); 
    FullMatrix<double> &       local_matrix = dinfo.matrix(0).matrix; 
    const std::vector<double> &JxW          = fe_values.get_JxW_values(); 

// 有了这些对象，我们像往常一样继续进行局部积分。首先，我们在正交点上循环，计算当前点的平流矢量。

    for (unsigned int point = 0; point < fe_values.n_quadrature_points; ++point) 
      { 
        const Tensor<1, dim> beta_at_q_point = 
          beta(fe_values.quadrature_point(point)); 

// 我们求解的是一个同质方程，因此在单元项中没有显示出右手。 剩下的就是对矩阵项的积分。

        for (unsigned int i = 0; i < fe_values.dofs_per_cell; ++i) 
          for (unsigned int j = 0; j < fe_values.dofs_per_cell; ++j) 
            local_matrix(i, j) += -beta_at_q_point *                // 
                                  fe_values.shape_grad(i, point) *  // 
                                  fe_values.shape_value(j, point) * // 
                                  JxW[point]; 
      } 
  } 

// 现在对边界条款也是如此。注意，现在我们使用FEValuesBase，即FEFaceValues和FESubfaceValues的基类，以便获得法向量。

  template <int dim> 
  void AdvectionProblem<dim>::integrate_boundary_term(DoFInfo & dinfo, 
                                                      CellInfo &info) 
  { 
    const FEValuesBase<dim> &fe_face_values = info.fe_values(); 
    FullMatrix<double> &     local_matrix   = dinfo.matrix(0).matrix; 
    Vector<double> &         local_vector   = dinfo.vector(0).block(0); 

    const std::vector<double> &        JxW = fe_face_values.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = 
      fe_face_values.get_normal_vectors(); 

    std::vector<double> g(fe_face_values.n_quadrature_points); 

    static BoundaryValues<dim> boundary_function; 
    boundary_function.value_list(fe_face_values.get_quadrature_points(), g); 

    for (unsigned int point = 0; point < fe_face_values.n_quadrature_points; 
         ++point) 
      { 
        const double beta_dot_n = 
          beta(fe_face_values.quadrature_point(point)) * normals[point]; 
        if (beta_dot_n > 0) 
          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i) 
            for (unsigned int j = 0; j < fe_face_values.dofs_per_cell; ++j) 
              local_matrix(i, j) += beta_dot_n *                           // 
                                    fe_face_values.shape_value(j, point) * // 
                                    fe_face_values.shape_value(i, point) * // 
                                    JxW[point]; 
        else 
          for (unsigned int i = 0; i < fe_face_values.dofs_per_cell; ++i) 
            local_vector(i) += -beta_dot_n *                          // 
                               g[point] *                             // 
                               fe_face_values.shape_value(i, point) * // 
                               JxW[point]; 
      } 
  } 

// 最后是内部面的条款。这里的区别是，我们收到了两个信息对象，相邻面的每个单元都有一个，我们组装了四个矩阵，每个单元一个，两个用于来回耦合。

  template <int dim> 
  void AdvectionProblem<dim>::integrate_face_term(DoFInfo & dinfo1, 
                                                  DoFInfo & dinfo2, 
                                                  CellInfo &info1, 
                                                  CellInfo &info2) 
  { 

// 对于正交点、权重等，我们使用第一个参数的FEValuesBase对象。

    const FEValuesBase<dim> &fe_face_values = info1.fe_values(); 
    const unsigned int       dofs_per_cell  = fe_face_values.dofs_per_cell; 

// 对于额外的形状函数，我们必须询问邻居的FEValuesBase。

    const FEValuesBase<dim> &fe_face_values_neighbor = info2.fe_values(); 
    const unsigned int       neighbor_dofs_per_cell = 
      fe_face_values_neighbor.dofs_per_cell; 

// 然后我们得到对四个局部矩阵的引用。字母u和v分别指的是试验和测试函数。%的数字表示由info1和info2提供的单元。按照惯例，每个信息对象中的两个矩阵指的是各自单元上的试验函数。第一个矩阵包含该单元的内部耦合，而第二个矩阵包含单元之间的耦合。

    FullMatrix<double> &u1_v1_matrix = dinfo1.matrix(0, false).matrix; 
    FullMatrix<double> &u2_v1_matrix = dinfo1.matrix(0, true).matrix; 
    FullMatrix<double> &u1_v2_matrix = dinfo2.matrix(0, true).matrix; 
    FullMatrix<double> &u2_v2_matrix = dinfo2.matrix(0, false).matrix; 

// 在这里，按照前面的函数，我们会有本地的右手边向量。幸运的是，界面条款只涉及到解决方案，右手边没有收到任何贡献。

    const std::vector<double> &        JxW = fe_face_values.get_JxW_values(); 
    const std::vector<Tensor<1, dim>> &normals = 
      fe_face_values.get_normal_vectors(); 

    for (unsigned int point = 0; point < fe_face_values.n_quadrature_points; 
         ++point) 
      { 
        const double beta_dot_n = 
          beta(fe_face_values.quadrature_point(point)) * normals[point]; 
        if (beta_dot_n > 0) 
          { 

// 这个词我们已经看过了。

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                u1_v1_matrix(i, j) += beta_dot_n *                           // 
                                      fe_face_values.shape_value(j, point) * // 
                                      fe_face_values.shape_value(i, point) * // 
                                      JxW[point]; 

// 我们另外组装术语  $(\beta\cdot n u,\hat v)_{\partial \kappa_+}$  。

            for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k) 
              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                u1_v2_matrix(k, j) += 
                  -beta_dot_n *                                   // 
                  fe_face_values.shape_value(j, point) *          // 
                  fe_face_values_neighbor.shape_value(k, point) * // 
                  JxW[point]; 
          } 
        else 
          { 

// 这个我们也已经看过了。

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l) 
                u2_v1_matrix(i, l) += 
                  beta_dot_n *                                    // 
                  fe_face_values_neighbor.shape_value(l, point) * // 
                  fe_face_values.shape_value(i, point) *          // 
                  JxW[point]; 

// 而这是另一个新的。 $(\beta\cdot n \hat u,\hat v)_{\partial \kappa_-}$  :

            for (unsigned int k = 0; k < neighbor_dofs_per_cell; ++k) 
              for (unsigned int l = 0; l < neighbor_dofs_per_cell; ++l) 
                u2_v2_matrix(k, l) += 
                  -beta_dot_n *                                   // 
                  fe_face_values_neighbor.shape_value(l, point) * // 
                  fe_face_values_neighbor.shape_value(k, point) * // 
                  JxW[point]; 
          } 
      } 
  } 
// @sect3{All the rest}  

// 对于这个简单的问题，我们使用了最简单的求解器，称为Richardson迭代，它代表了简单的缺陷修正。这与一个块状SSOR预处理相结合，该预处理使用DG离散化产生的系统矩阵的特殊块状结构。这些块的大小是每个单元的DoF数量。在这里，我们使用SSOR预处理，因为我们没有根据流场对DoFs进行重新编号。如果在流的下游方向对DoFs进行重新编号，那么块状的Gauss-Seidel预处理（见PreconditionBlockSOR类，放松=1）会做得更好。

  template <int dim> 
  void AdvectionProblem<dim>::solve(Vector<double> &solution) 
  { 
    SolverControl                    solver_control(1000, 1e-12); 
    SolverRichardson<Vector<double>> solver(solver_control); 

// 这里我们创建了预处理程序。

    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 

// 然后将矩阵分配给它，并设置正确的块大小。

    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 

// 做完这些准备工作后，我们就可以启动线性求解器了。

    solver.solve(system_matrix, solution, right_hand_side, preconditioner); 
  } 

// 我们根据一个非常简单的细化标准来细化网格，即对解的梯度的近似。由于这里我们考虑的是DG(1)方法（即我们使用片状双线性形状函数），我们可以简单地计算每个单元的梯度。但是我们并不希望我们的细化指标只建立在每个单元的梯度上，而是希望同时建立在相邻单元之间的不连续解函数的跳跃上。最简单的方法是通过差分商计算近似梯度，包括考虑中的单元和其相邻的单元。这是由 <code>DerivativeApproximation</code> 类完成的，它计算近似梯度的方式类似于本教程 step-9 中描述的 <code>GradientEstimation</code> 。事实上， <code>DerivativeApproximation</code> 类是在 step-9 的 <code>GradientEstimation</code> 类之后开发的。与  step-9  中的讨论相关，这里我们考虑  $h^{1+d/2}|\nabla_h u_h|$  。此外，我们注意到，我们不考虑近似的二次导数，因为线性平流方程的解一般不在 $H^2$ 中，而只在 $H^1$ 中（或者，更准确地说：在 $H^1_\beta$ 中，即在方向 $\beta$ 中的导数是可平方整除的函数空间）。

  template <int dim> 
  void AdvectionProblem<dim>::refine_grid() 
  { 

//  <code>DerivativeApproximation</code> 类将梯度计算为浮点精度。这已经足够了，因为它们是近似的，只作为细化指标。

    Vector<float> gradient_indicator(triangulation.n_active_cells()); 

// 现在，近似梯度被计算出来了

    DerivativeApproximation::approximate_gradient(mapping, 
                                                  dof_handler, 
                                                  solution, 
                                                  gradient_indicator); 

//并且它们被单元格按比例放大，系数为 $h^{1+d/2}$  。
    unsigned int cell_no = 0; 
    for (const auto &cell : dof_handler.active_cell_iterators()) 
      gradient_indicator(cell_no++) *= 
        std::pow(cell->diameter(), 1 + 1.0 * dim / 2); 

// 最后它们作为细化指标。

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    gradient_indicator, 
                                                    0.3, 
                                                    0.1); 

    triangulation.execute_coarsening_and_refinement(); 
  } 

// 这个程序的输出包括自适应细化网格的eps文件和gnuplot格式的数值解。

  template <int dim> 
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
  { 

// 首先将网格写成eps格式。

    { 
      const std::string filename = "grid-" + std::to_string(cycle) + ".eps"; 
      deallog << "Writing grid to <" << filename << ">" << std::endl; 
      std::ofstream eps_output(filename); 

      GridOut grid_out; 
      grid_out.write_eps(triangulation, eps_output); 
    } 

// 然后以gnuplot格式输出解决方案。

    { 
      const std::string filename = "sol-" + std::to_string(cycle) + ".gnuplot"; 
      deallog << "Writing solution to <" << filename << ">" << std::endl; 
      std::ofstream gnuplot_output(filename); 

      DataOut<dim> data_out; 
      data_out.attach_dof_handler(dof_handler); 
      data_out.add_data_vector(solution, "u"); 

      data_out.build_patches(); 

      data_out.write_gnuplot(gnuplot_output); 
    } 
  } 

// 下面的 <code>run</code> 函数与前面的例子类似。

  template <int dim> 
  void AdvectionProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      { 
        deallog << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube(triangulation); 

            triangulation.refine_global(3); 
          } 
        else 
          refine_grid(); 

        deallog << "Number of active cells:       " 
                << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        deallog << "Number of degrees of freedom: " << dof_handler.n_dofs() 
                << std::endl; 

        assemble_system(); 
        solve(solution); 

        output_results(cycle); 
      } 
  } 
} // namespace Step12 

// 下面的 <code>main</code> 函数与前面的例子也类似，不需要注释。

int main() 
{ 
  try 
    { 
      dealii::deallog.depth_console(5); 

      Step12::AdvectionProblem<2> dgmethod; 
      dgmethod.run(); 
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


