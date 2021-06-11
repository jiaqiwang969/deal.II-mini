CCTest_file/step-47.cc

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
 * Authors: Natasha Sharma, University of Texas at El Paso, 
 *          Guido Kanschat, University of Heidelberg 
 *          Timo Heister, Clemson University 
 *          Wolfgang Bangerth, Colorado State University 
 *          Zhuroan Wang, Colorado State University 
 */ 


// @sect3{Include files}  

// 前面的几个include文件已经在前面的例子中使用过了，所以我们在这里不再解释它们的含义。该程序的主要结构与例如 step-4 的结构非常相似，因此我们包含了许多相同的头文件。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/function.h> 

#include <deal.II/lac/vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/dynamic_sparsity_pattern.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/sparse_matrix.h> 
#include <deal.II/lac/sparse_direct.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_values.h> 
#include <deal.II/fe/mapping_q.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/data_out.h> 

// 最有趣的两个头文件将是这两个。

#include <deal.II/fe/fe_interface_values.h> 
#include <deal.II/meshworker/mesh_loop.h> 

// 其中第一个文件负责提供FEInterfaceValues类，该类可用于评估单元间界面的形状函数（或其梯度）的跳跃或平均值等数量。这个类在评估C0IP公式中出现的惩罚项时将相当有用。

#include <fstream> 
#include <iostream> 
#include <cmath> 

namespace Step47 
{ 
  using namespace dealii; 

// 在下面的命名空间中，让我们定义精确解，我们将与数值计算的解进行比较。它的形式是 $u(x,y) = \sin(\pi x) \sin(\pi y)$ （只实现了2d的情况），该命名空间还包含一个对应于产生该解的右手边的类。

  namespace ExactSolution 
  { 
    using numbers::PI; 

    template <int dim> 
    class Solution : public Function<dim> 
    { 
    public: 
      static_assert(dim == 2, "Only dim==2 is implemented."); 

      virtual double value(const Point<dim> &p, 
                           const unsigned int /*component*/ = 0) const override 
      { 
        return std::sin(PI * p[0]) * std::sin(PI * p[1]); 
      } 

      virtual Tensor<1, dim> 
      gradient(const Point<dim> &p, 
               const unsigned int /*component*/ = 0) const override 
      { 
        Tensor<1, dim> r; 
        r[0] = PI * std::cos(PI * p[0]) * std::sin(PI * p[1]); 
        r[1] = PI * std::cos(PI * p[1]) * std::sin(PI * p[0]); 
        return r; 
      } 

      virtual void 
      hessian_list(const std::vector<Point<dim>> &       points, 
                   std::vector<SymmetricTensor<2, dim>> &hessians, 
                   const unsigned int /*component*/ = 0) const override 
      { 
        for (unsigned i = 0; i < points.size(); ++i) 
          { 
            const double x = points[i][0]; 
            const double y = points[i][1]; 

            hessians[i][0][0] = -PI * PI * std::sin(PI * x) * std::sin(PI * y); 
            hessians[i][0][1] = PI * PI * std::cos(PI * x) * std::cos(PI * y); 
            hessians[i][1][1] = -PI * PI * std::sin(PI * x) * std::sin(PI * y); 
          } 
      } 
    }; 

    template <int dim> 
    class RightHandSide : public Function<dim> 
    { 
    public: 
      static_assert(dim == 2, "Only dim==2 is implemented"); 

      virtual double value(const Point<dim> &p, 
                           const unsigned int /*component*/ = 0) const override 

      { 
        return 4 * std::pow(PI, 4.0) * std::sin(PI * p[0]) * 
               std::sin(PI * p[1]); 
      } 
    }; 
  } // namespace ExactSolution 

//  @sect3{The main class}  

// 以下是本教程程序的主类。它具有许多其他教程程序的结构，其内容和后面的构造函数应该没有什么特别令人惊讶的地方。

  template <int dim> 
  class BiharmonicProblem 
  { 
  public: 
    BiharmonicProblem(const unsigned int fe_degree); 

    void run(); 

  private: 
    void make_grid(); 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void compute_errors(); 
    void output_results(const unsigned int iteration) const; 

    Triangulation<dim> triangulation; 

    MappingQ<dim> mapping; 

    FE_Q<dim>                 fe; 
    DoFHandler<dim>           dof_handler; 
    AffineConstraints<double> constraints; 

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> system_rhs; 
  }; 

  template <int dim> 
  BiharmonicProblem<dim>::BiharmonicProblem(const unsigned int fe_degree) 
    : mapping(1) 
    , fe(fe_degree) 
    , dof_handler(triangulation) 
  {} 

// 接下来是创建初始网格（一次精炼的单元格）和设置每个网格的约束、向量和矩阵的函数。同样，这两个函数与之前的许多教程程序基本没有变化。

  template <int dim> 
  void BiharmonicProblem<dim>::make_grid() 
  { 
    GridGenerator::hyper_cube(triangulation, 0., 1.); 
    triangulation.refine_global(1); 

    std::cout << "Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "Total number of cells: " << triangulation.n_cells() 
              << std::endl; 
  } 

  template <int dim> 
  void BiharmonicProblem<dim>::setup_system() 
  { 
    dof_handler.distribute_dofs(fe); 

    std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << std::endl; 

    constraints.clear(); 
    DoFTools::make_hanging_node_constraints(dof_handler, constraints); 

    VectorTools::interpolate_boundary_values(dof_handler, 
                                             0, 
                                             ExactSolution::Solution<dim>(), 
                                             constraints); 
    constraints.close(); 

    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, true); 
    sparsity_pattern.copy_from(dsp); 
    system_matrix.reinit(sparsity_pattern); 

    solution.reinit(dof_handler.n_dofs()); 
    system_rhs.reinit(dof_handler.n_dofs()); 
  } 

//  @sect4{Assembling the linear system}  

// 下面的几段代码更有意思。它们都与线性系统的组装有关。虽然组装单元格内部项的难度不大--这在本质上就像组装拉普拉斯方程的相应项一样，你已经在 step-4 或 step-6 中看到了这是如何工作的，例如，困难在于公式中的惩罚项。这需要在单元格的界面上对形状函数的梯度进行评估。因此，至少需要使用两个FEFaceValues对象，但如果其中一个面是自适应细化的，那么实际上需要一个FEFaceValues和一个FESubfaceValues对象；我们还需要跟踪哪些形状函数在哪里，最后我们需要确保每个面只被访问一次。所有这些对于我们真正想要实现的逻辑（即双线性形式中的惩罚项）来说都是一笔不小的开销。因此，我们将使用FEInterfaceValues类--这是deal.II中的一个辅助类，它允许我们抽象出两个FEFaceValues或FESubfaceValues对象，直接访问我们真正关心的东西：跳跃、平均等。

// 但这还没有解决我们的问题，即当我们在所有单元格和它们的所有面中循环时，必须跟踪我们已经访问过哪些面。为了使这个过程更简单，我们使用了 MeshWorker::mesh_loop() 函数，它为这个任务提供了一个简单的接口：基于WorkStream命名空间文档中概述的想法， MeshWorker::mesh_loop() 需要三个函数对单元、内部面和边界面进行工作。这些函数在抓取对象上工作，以获得中间结果，然后将其计算结果复制到复制数据对象中，由一个复制器函数将其复制到全局矩阵和右侧对象中。

// 然后，下面的结构提供了这种方法所需的从头开始和复制对象。你可以查阅WorkStream命名空间以及 @ref threads "多处理器并行计算 "模块，了解更多关于它们通常如何工作的信息。

  template <int dim> 
  struct ScratchData 
  { 
    ScratchData(const Mapping<dim> &      mapping, 
                const FiniteElement<dim> &fe, 
                const unsigned int        quadrature_degree, 
                const UpdateFlags         update_flags, 
                const UpdateFlags         interface_update_flags) 
      : fe_values(mapping, fe, QGauss<dim>(quadrature_degree), update_flags) 
      , fe_interface_values(mapping, 
                            fe, 
                            QGauss<dim - 1>(quadrature_degree), 
                            interface_update_flags) 
    {} 

    ScratchData(const ScratchData<dim> &scratch_data) 
      : fe_values(scratch_data.fe_values.get_mapping(), 
                  scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  scratch_data.fe_values.get_update_flags()) 
      , fe_interface_values(scratch_data.fe_values.get_mapping(), 
                            scratch_data.fe_values.get_fe(), 
                            scratch_data.fe_interface_values.get_quadrature(), 
                            scratch_data.fe_interface_values.get_update_flags()) 
    {} 

    FEValues<dim>          fe_values; 
    FEInterfaceValues<dim> fe_interface_values; 
  }; 

  struct CopyData 
  { 
    CopyData(const unsigned int dofs_per_cell) 
      : cell_matrix(dofs_per_cell, dofs_per_cell) 
      , cell_rhs(dofs_per_cell) 
      , local_dof_indices(dofs_per_cell) 
    {} 

    CopyData(const CopyData &) = default; 

    CopyData(CopyData &&) = default; 

    ~CopyData() = default; 

    CopyData &operator=(const CopyData &) = default; 

    CopyData &operator=(CopyData &&) = default; 

    struct FaceData 
    { 
      FullMatrix<double>                   cell_matrix; 
      std::vector<types::global_dof_index> joint_dof_indices; 
    }; 

    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_rhs; 
    std::vector<types::global_dof_index> local_dof_indices; 
    std::vector<FaceData>                face_data; 
  }; 

// 更有趣的部分是我们实际组装线性系统的地方。从根本上说，这个函数有五个部分。

// - `cell_worker`λ函数的定义，这是一个定义在`assemble_system()`函数中的小函数，它将负责计算单个单元上的局部积分。它将在`ScratchData`类的副本上工作，并将其结果放入相应的`CopyData`对象。

// - `face_worker` lambda函数的定义，它将对单元格之间的界面上的所有项进行积分。

// - 定义了`boundary_worker`函数，对位于域的边界上的单元面做同样的工作。

// - `copier`函数的定义，该函数负责将前面三个函数中的所有数据复制到单个单元的复制对象中，并复制到全局矩阵和右侧。

// 第五部分是我们把所有这些都集中在一起。

// 让我们轮流浏览一下这些组装所需的每一块。

  template <int dim> 
  void BiharmonicProblem<dim>::assemble_system() 
  { 
    using Iterator = typename DoFHandler<dim>::active_cell_iterator; 

// 第一部分是`cell_worker'，它在细胞内部进行组装。它是一个（lambda）函数，以一个单元格（输入）、一个抓取对象和一个复制对象（输出）为参数。它看起来像许多其他教程程序的装配函数，或者至少是所有单元格的循环主体。

// 我们在这里整合的条款是单元格对全局矩阵的贡献
// @f{align*}{
//     A^K_{ij} = \int_K \nabla^2\varphi_i(x) : \nabla^2\varphi_j(x) dx
//  @f} ，
//  以及对右侧向量的贡献
//  @f{align*}{
//     f^K_i = \int_K \varphi_i(x) f(x) dx
//  @f}

// 我们使用与组装 step-22 相同的技术来加速该函数。我们不在最里面的循环中调用`fe_values.shape_hessian(i, qpoint)`，而是创建一个变量`hessian_i`，在循环中对`i`进行一次评估，在循环中对`j`重新使用如此评估的值。为了对称，我们对变量`hessian_j`也做了同样的处理，尽管它确实只用了一次，而且我们可以在计算两个项之间标量乘积的指令中留下对`fe_values.shape_hessian(j,qpoint)`的调用。

    auto cell_worker = [&](const Iterator &  cell, 
                           ScratchData<dim> &scratch_data, 
                           CopyData &        copy_data) { 
      copy_data.cell_matrix = 0; 
      copy_data.cell_rhs    = 0; 

      FEValues<dim> &fe_values = scratch_data.fe_values; 
      fe_values.reinit(cell); 

      cell->get_dof_indices(copy_data.local_dof_indices); 

      const ExactSolution::RightHandSide<dim> right_hand_side; 

      const unsigned int dofs_per_cell = 
        scratch_data.fe_values.get_fe().n_dofs_per_cell(); 

      for (unsigned int qpoint = 0; qpoint < fe_values.n_quadrature_points; 
           ++qpoint) 
        { 
          for (unsigned int i = 0; i < dofs_per_cell; ++i) 
            { 
              const Tensor<2, dim> &hessian_i = 
                fe_values.shape_hessian(i, qpoint); 

              for (unsigned int j = 0; j < dofs_per_cell; ++j) 
                { 
                  const Tensor<2, dim> &hessian_j = 
                    fe_values.shape_hessian(j, qpoint); 

                  copy_data.cell_matrix(i, j) += 
                    scalar_product(hessian_i,   // nabla^2 phi_i(x) 
                                   hessian_j) * // nabla^2 phi_j(x) 
                    fe_values.JxW(qpoint);      // dx 
                } 

              copy_data.cell_rhs(i) += 
                fe_values.shape_value(i, qpoint) * // phi_i(x) 
                right_hand_side.value( 
                  fe_values.quadrature_point(qpoint)) * // f(x) 
                fe_values.JxW(qpoint);                  // dx 
            } 
        } 
    }; 

// 下一个构建模块是在网格的每个内部面组装惩罚项。正如 MeshWorker::mesh_loop(), 文档中所描述的，这个函数接收到的参数表示一个单元和它的相邻单元，以及（对于这两个单元中的每一个）我们必须整合的面（以及潜在的子面）。同样地，我们也得到了一个从头开始的对象，以及一个用于放置结果的拷贝对象。

// 这个函数本身有三个部分。在顶部，我们初始化FEInterfaceValues对象，并创建一个新的 `CopyData::FaceData` 对象来存储我们的输入。这将被推到`copy_data.face_data`变量的末尾。我们需要这样做，因为我们对一个给定单元进行积分的面（或子面）的数量因单元而异，而且这些矩阵的大小也不同，取决于面或子面相邻的自由度。正如 MeshWorker::mesh_loop(), 文档中所讨论的，每次访问一个新的单元时，复制对象都会被重置，所以我们推到`copy_data.face_data()`末尾的内容实际上就是后来的`copier`函数在复制每个单元的贡献到全局矩阵和右侧对象时所能看到的。

    auto face_worker = [&](const Iterator &    cell, 
                           const unsigned int &f, 
                           const unsigned int &sf, 
                           const Iterator &    ncell, 
                           const unsigned int &nf, 
                           const unsigned int &nsf, 
                           ScratchData<dim> &  scratch_data, 
                           CopyData &          copy_data) { 
      FEInterfaceValues<dim> &fe_interface_values = 
        scratch_data.fe_interface_values; 
      fe_interface_values.reinit(cell, f, sf, ncell, nf, nsf); 

      copy_data.face_data.emplace_back(); 
      CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 

      copy_data_face.joint_dof_indices = 
        fe_interface_values.get_interface_dof_indices(); 

      const unsigned int n_interface_dofs = 
        fe_interface_values.n_current_interface_dofs(); 
      copy_data_face.cell_matrix.reinit(n_interface_dofs, n_interface_dofs); 

// 第二部分涉及到确定惩罚参数应该是什么。通过观察双线性形式中各种项的单位，很明显，惩罚必须具有 $\frac{\gamma}{h_K}$ 的形式（即，超过长度尺度的一个），但如何选择无维数 $\gamma$ 并不是先验的。从拉普拉斯方程的不连续Galerkin理论来看，人们可能猜想正确的选择是 $\gamma=p(p+1)$ 是正确的选择，其中 $p$ 是所用有限元的多项式程度。我们将在本程序的结果部分更详细地讨论这个选择。

// 在上面的公式中， $h_K$  是单元格  $K$  的大小。但这也不是很简单的事情。如果使用高度拉伸的单元格，那么一个更复杂的理论说， $h$ 应该被单元格 $K$ 的直径取代，该直径是有关边缘方向的法线。 事实证明，在deal.II中有一个函数用于此。其次，当从一个面的两个不同侧面看时， $h_K$ 可能是不同的。

// 为了安全起见，我们取这两个值的最大值。我们将注意到，如果使用自适应网格细化所产生的悬空节点，有可能需要进一步调整这一计算方法。

      const unsigned int p = fe.degree; 
      const double       gamma_over_h = 
        std::max((1.0 * p * (p + 1) / 
                  cell->extent_in_direction( 
                    GeometryInfo<dim>::unit_normal_direction[f])), 
                 (1.0 * p * (p + 1) / 
                  ncell->extent_in_direction( 
                    GeometryInfo<dim>::unit_normal_direction[nf]))); 

// 最后，像往常一样，我们在正交点和指数`i`和`j`上循环，把这个面或子面的贡献加起来。然后将这些数据存储在上面创建的`copy_data.face_data`对象中。至于单元格工作者，如果可能的话，我们将平均数和跳跃的评估从循环中拉出来，引入局部变量来存储这些结果。然后组件只需要在最里面的循环中使用这些局部变量。关于这段代码实现的具体公式，回顾一下，双线性形式的接口项如下。
// @f{align*}{
//   -\sum_{e \in \mathbb{F}} \int_{e}
//   \jump{ \frac{\partial v_h}{\partial \mathbf n}}
//   \average{\frac{\partial^2 u_h}{\partial \mathbf n^2}} \ ds
//  -\sum_{e \in \mathbb{F}} \int_{e}
//  \average{\frac{\partial^2 v_h}{\partial \mathbf n^2}}
//  \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds
//  + \sum_{e \in \mathbb{F}}
//  \frac{\gamma}{h_e}
//  \int_e
//  \jump{\frac{\partial v_h}{\partial \mathbf n}}
//  \jump{\frac{\partial u_h}{\partial \mathbf n}} \ ds.
//  @f}

      for (unsigned int qpoint = 0; 
           qpoint < fe_interface_values.n_quadrature_points; 
           ++qpoint) 
        { 
          const auto &n = fe_interface_values.normal(qpoint); 

          for (unsigned int i = 0; i < n_interface_dofs; ++i) 
            { 
              const double av_hessian_i_dot_n_dot_n = 
                (fe_interface_values.average_hessian(i, qpoint) * n * n); 
              const double jump_grad_i_dot_n = 
                (fe_interface_values.jump_gradient(i, qpoint) * n); 

              for (unsigned int j = 0; j < n_interface_dofs; ++j) 
                { 
                  const double av_hessian_j_dot_n_dot_n = 
                    (fe_interface_values.average_hessian(j, qpoint) * n * n); 
                  const double jump_grad_j_dot_n = 
                    (fe_interface_values.jump_gradient(j, qpoint) * n); 

                  copy_data_face.cell_matrix(i, j) += 
                    (-av_hessian_i_dot_n_dot_n       // - {grad^2 v n n } 
                       * jump_grad_j_dot_n           // [grad u n] 
                     - av_hessian_j_dot_n_dot_n      // - {grad^2 u n n } 
                         * jump_grad_i_dot_n         // [grad v n] 
                     +                               // + 
                     gamma_over_h *                  // gamma/h 
                       jump_grad_i_dot_n *           // [grad v n] 
                       jump_grad_j_dot_n) *          // [grad u n] 
                    fe_interface_values.JxW(qpoint); // dx 
                } 
            } 
        } 
    }; 

// 第三块是对处于边界的面做同样的装配。当然，想法和上面一样，唯一不同的是，现在有惩罚条款也进入了右手边。

// 和以前一样，这个函数的第一部分只是设置了一些辅助对象。

    auto boundary_worker = [&](const Iterator &    cell, 
                               const unsigned int &face_no, 
                               ScratchData<dim> &  scratch_data, 
                               CopyData &          copy_data) { 
      FEInterfaceValues<dim> &fe_interface_values = 
        scratch_data.fe_interface_values; 
      fe_interface_values.reinit(cell, face_no); 
      const auto &q_points = fe_interface_values.get_quadrature_points(); 

      copy_data.face_data.emplace_back(); 
      CopyData::FaceData &copy_data_face = copy_data.face_data.back(); 

      const unsigned int n_dofs = 
        fe_interface_values.n_current_interface_dofs(); 
      copy_data_face.joint_dof_indices = 
        fe_interface_values.get_interface_dof_indices(); 

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs); 

      const std::vector<double> &JxW = fe_interface_values.get_JxW_values(); 
      const std::vector<Tensor<1, dim>> &normals = 
        fe_interface_values.get_normal_vectors(); 

      const ExactSolution::Solution<dim> exact_solution; 
      std::vector<Tensor<1, dim>>        exact_gradients(q_points.size()); 
      exact_solution.gradient_list(q_points, exact_gradients); 

// 从正面看，由于我们现在只处理与面相邻的一个单元（因为我们在边界上），惩罚因子 $\gamma$ 的计算大大简化了。

      const unsigned int p = fe.degree; 
      const double       gamma_over_h = 
        (1.0 * p * (p + 1) / 
         cell->extent_in_direction( 
           GeometryInfo<dim>::unit_normal_direction[face_no])); 

// 第三块是术语的组合。由于这些条款包含了矩阵的条款和右手边的条款，所以现在稍微有些麻烦。前者与上面所说的内部面完全相同，如果我们只是适当地定义了跳跃和平均（这就是FEInterfaceValues类所做的）。后者需要我们评估边界条件 $j(\mathbf x)$ ，在当前情况下（我们知道确切的解决方案），我们从 $j(\mathbf x) = \frac{\partial u(\mathbf x)}{\partial {\mathbf n}}$ 中计算出来。然后，要添加到右侧向量的项是  $\frac{\gamma}{h_e}\int_e \jump{\frac{\partial v_h}{\partial \mathbf n}} j \ ds$  。

      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint) 
        { 
          const auto &n = normals[qpoint]; 

          for (unsigned int i = 0; i < n_dofs; ++i) 
            { 
              const double av_hessian_i_dot_n_dot_n = 
                (fe_interface_values.average_hessian(i, qpoint) * n * n); 
              const double jump_grad_i_dot_n = 
                (fe_interface_values.jump_gradient(i, qpoint) * n); 

              for (unsigned int j = 0; j < n_dofs; ++j) 
                { 
                  const double av_hessian_j_dot_n_dot_n = 
                    (fe_interface_values.average_hessian(j, qpoint) * n * n); 
                  const double jump_grad_j_dot_n = 
                    (fe_interface_values.jump_gradient(j, qpoint) * n); 

                  copy_data_face.cell_matrix(i, j) += 
                    (-av_hessian_i_dot_n_dot_n  // - {grad^2 v n n} 
                       * jump_grad_j_dot_n      //   [grad u n] 

//                                      

                     - av_hessian_j_dot_n_dot_n // - {grad^2 u n n} 
                         * jump_grad_i_dot_n    //   [grad v n] 

//                                      

                     + gamma_over_h             //  gamma/h 
                         * jump_grad_i_dot_n    // [grad v n] 
                         * jump_grad_j_dot_n    // [grad u n] 
                     ) * 
                    JxW[qpoint]; // dx 
                } 

              copy_data.cell_rhs(i) += 
                (-av_hessian_i_dot_n_dot_n *       // - {grad^2 v n n } 
                   (exact_gradients[qpoint] * n)   //   (grad u_exact . n) 
                 +                                 // + 
                 gamma_over_h                      //  gamma/h 
                   * jump_grad_i_dot_n             // [grad v n] 
                   * (exact_gradients[qpoint] * n) // (grad u_exact . n) 
                 ) * 
                JxW[qpoint]; // dx 
            } 
        } 
    }; 

// 第四部分是一个小函数，它将上面的单元格、内部和边界面装配程序产生的数据复制到全局矩阵和右手向量中。这里真的没有什么可做的。我们分配单元格矩阵和右侧贡献，就像我们在其他几乎所有的教程程序中使用约束对象那样。然后，我们还必须对面矩阵的贡献做同样的处理，这些贡献已经获得了面（内部和边界）的内容，并且`面_工作`和`边界_工作`已经添加到`copy_data.face_data`阵列中。

    auto copier = [&](const CopyData &copy_data) { 
      constraints.distribute_local_to_global(copy_data.cell_matrix, 
                                             copy_data.cell_rhs, 
                                             copy_data.local_dof_indices, 
                                             system_matrix, 
                                             system_rhs); 

      for (auto &cdf : copy_data.face_data) 
        { 
          constraints.distribute_local_to_global(cdf.cell_matrix, 
                                                 cdf.joint_dof_indices, 
                                                 system_matrix); 
        } 
    }; 

// 在设置了所有这些之后，剩下的就是创建一个从头开始和复制数据的对象，并调用 MeshWorker::mesh_loop() 函数，然后遍历所有的单元格和面，调用它们各自的工作器，然后是复制器函数，将东西放入全局矩阵和右侧。作为一个额外的好处， MeshWorker::mesh_loop() 以并行方式完成所有这些工作，使用你的机器恰好有多少个处理器核心。

    const unsigned int n_gauss_points = dof_handler.get_fe().degree + 1; 
    ScratchData<dim>   scratch_data(mapping, 
                                  fe, 
                                  n_gauss_points, 
                                  update_values | update_gradients | 
                                    update_hessians | update_quadrature_points | 
                                    update_JxW_values, 
                                  update_values | update_gradients | 
                                    update_hessians | update_quadrature_points | 
                                    update_JxW_values | update_normal_vectors); 
    CopyData           copy_data(dof_handler.get_fe().n_dofs_per_cell()); 
    MeshWorker::mesh_loop(dof_handler.begin_active(), 
                          dof_handler.end(), 
                          cell_worker, 
                          copier, 
                          scratch_data, 
                          copy_data, 
                          MeshWorker::assemble_own_cells | 
                            MeshWorker::assemble_boundary_faces | 
                            MeshWorker::assemble_own_interior_faces_once, 
                          boundary_worker, 
                          face_worker); 
  } 

//  @sect4{Solving the linear system and postprocessing}  

// 到此为止，节目基本上结束了。其余的函数并不太有趣或新颖。第一个函数只是用一个直接求解器来求解线性系统（也见 step-29  ）。

  template <int dim> 
  void BiharmonicProblem<dim>::solve() 
  { 
    std::cout << "   Solving system..." << std::endl; 

    SparseDirectUMFPACK A_direct; 
    A_direct.initialize(system_matrix); 
    A_direct.vmult(solution, system_rhs); 

    constraints.distribute(solution); 
  } 

// 下一个函数评估了计算出的解和精确解之间的误差（在这里是已知的，因为我们选择了右手边和边界值的方式，所以我们知道相应的解）。在下面的前两个代码块中，我们计算了 $L_2$ 准则和 $H^1$ 半准则下的误差。

  template <int dim> 
  void BiharmonicProblem<dim>::compute_errors() 
  { 
    { 
      Vector<float> norm_per_cell(triangulation.n_active_cells()); 
      VectorTools::integrate_difference(mapping, 
                                        dof_handler, 
                                        solution, 
                                        ExactSolution::Solution<dim>(), 
                                        norm_per_cell, 
                                        QGauss<dim>(fe.degree + 2), 
                                        VectorTools::L2_norm); 
      const double error_norm = 
        VectorTools::compute_global_error(triangulation, 
                                          norm_per_cell, 
                                          VectorTools::L2_norm); 
      std::cout << "   Error in the L2 norm           :     " << error_norm 
                << std::endl; 
    } 

    { 
      Vector<float> norm_per_cell(triangulation.n_active_cells()); 
      VectorTools::integrate_difference(mapping, 
                                        dof_handler, 
                                        solution, 
                                        ExactSolution::Solution<dim>(), 
                                        norm_per_cell, 
                                        QGauss<dim>(fe.degree + 2), 
                                        VectorTools::H1_seminorm); 
      const double error_norm = 
        VectorTools::compute_global_error(triangulation, 
                                          norm_per_cell, 
                                          VectorTools::H1_seminorm); 
      std::cout << "   Error in the H1 seminorm       : " << error_norm 
                << std::endl; 
    } 

// 现在也计算一下 $H^2$ 半正态误差的近似值。实际的 $H^2$ 半规范要求我们对解决方案 $u_h$ 的二阶导数进行积分，但是考虑到我们使用的拉格朗日形状函数， $u_h$ 当然在单元间的界面上有结点，因此二阶导数在界面是奇异的。因此，我们实际上只对单元的内部进行积分，而忽略了界面的贡献。这不是*等同于问题的能量准则，但是仍然可以让我们了解误差收敛的速度。

// 我们注意到，我们可以通过定义一个等同于能量准则的准则来解决这个问题。这将涉及到不仅要像我们下面做的那样将细胞内部的积分相加，而且还要为 $u_h$ 的导数在界面上的跳跃添加惩罚项，并对这两种项进行适当的缩放。我们将把这个问题留给以后的工作。

    { 
      const QGauss<dim>            quadrature_formula(fe.degree + 2); 
      ExactSolution::Solution<dim> exact_solution; 
      Vector<double> error_per_cell(triangulation.n_active_cells()); 

      FEValues<dim> fe_values(mapping, 
                              fe, 
                              quadrature_formula, 
                              update_values | update_hessians | 
                                update_quadrature_points | update_JxW_values); 

      FEValuesExtractors::Scalar scalar(0); 
      const unsigned int         n_q_points = quadrature_formula.size(); 

      std::vector<SymmetricTensor<2, dim>> exact_hessians(n_q_points); 
      std::vector<Tensor<2, dim>>          hessians(n_q_points); 
      for (auto &cell : dof_handler.active_cell_iterators()) 
        { 
          fe_values.reinit(cell); 
          fe_values[scalar].get_function_hessians(solution, hessians); 
          exact_solution.hessian_list(fe_values.get_quadrature_points(), 
                                      exact_hessians); 

          double local_error = 0; 
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
            { 
              local_error += 
                ((exact_hessians[q_point] - hessians[q_point]).norm_square() * 
                 fe_values.JxW(q_point)); 
            } 
          error_per_cell[cell->active_cell_index()] = std::sqrt(local_error); 
        } 

      const double error_norm = error_per_cell.l2_norm(); 
      std::cout << "   Error in the broken H2 seminorm: " << error_norm 
                << std::endl; 
    } 
  } 

// 同样无趣的是生成图形输出的函数。它看起来和  step-6  中的一模一样，比如说。

  template <int dim> 
  void 
  BiharmonicProblem<dim>::output_results(const unsigned int iteration) const 
  { 
    std::cout << "   Writing graphical output..." << std::endl; 

    DataOut<dim> data_out; 

    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, "solution"); 
    data_out.build_patches(); 

    const std::string filename = 
      ("output_" + Utilities::int_to_string(iteration, 6) + ".vtu"); 
    std::ofstream output_vtu(filename); 
    data_out.write_vtu(output_vtu); 
  } 

// `run()`函数的情况也是如此。就像在以前的程序中一样。

  template <int dim> 
  void BiharmonicProblem<dim>::run() 
  { 
    make_grid(); 

    const unsigned int n_cycles = 4; 
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        std::cout << "Cycle " << cycle << " of " << n_cycles << std::endl; 

        triangulation.refine_global(1); 
        setup_system(); 

        assemble_system(); 
        solve(); 

        output_results(cycle); 

        compute_errors(); 
        std::cout << std::endl; 
      } 
  } 
} // namespace Step47 

//  @sect3{The main() function}  

// 最后是 "main() "函数。同样，这里没有什么可看的。它看起来和以前的教程程序中的一样。有一个变量，可以选择我们要用来解方程的元素的多项式程度。因为我们使用的C0IP公式要求元素的度数至少为2，所以我们用一个断言来检查，无论为多项式度数设置什么都是有意义的。

int main() 
{ 
  try 
    { 
      using namespace dealii; 
      using namespace Step47; 

      const unsigned int fe_degree = 2; 
      Assert(fe_degree >= 2, 
             ExcMessage("The C0IP formulation for the biharmonic problem " 
                        "only works if one uses elements of polynomial " 
                        "degree at least 2.")); 

      BiharmonicProblem<2> biharmonic_problem(fe_degree); 
      biharmonic_problem.run(); 
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


