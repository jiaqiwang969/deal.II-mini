CCTest_file/step-12.cc

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
 *         Timo Heister, Clemson University, 2019 
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
#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/dofs/dof_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/fe/mapping_q1.h> 

// 这里定义了不连续的有限元。它们的使用方式与所有其他有限元相同，不过--正如你在以前的教程程序中所看到的--用户与有限元类的交互根本不多：它们被传递给 <code>DoFHandler</code> 和 <code>FEValues</code> 对象，仅此而已。

#include <deal.II/fe/fe_dgq.h> 

// FEInterfaceValues需要这个头来计算界面上的积分。

#include <deal.II/fe/fe_interface_values.h> 

// 我们将使用最简单的求解器，称为Richardson迭代，它代表了一个简单的缺陷修正。这与一个块状SSOR预处理器（定义在precondition_block.h中）相结合，该预处理器使用DG离散产生的系统矩阵的特殊块状结构。

#include <deal.II/lac/solver_richardson.h> 
#include <deal.II/lac/precondition_block.h> 

// 我们将使用梯度作为细化指标。

#include <deal.II/numerics/derivative_approximation.h> 

// 最后，新的包含文件用于使用MeshWorker框架中的Mesh_loop。

#include <deal.II/meshworker/mesh_loop.h> 

// 像所有的程序一样，我们在完成这一部分时，要包括所需的C++头文件，并声明我们要使用dealii命名空间中的对象，不含前缀。

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

// 考虑到流动方向，单位方块 $[0,1]^2$ 的流入边界是右边界和下边界。我们在x轴上规定了不连续的边界值1和0，在右边界上规定了值0。该函数在流出边界上的值将不会在DG方案中使用。

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

    if (wind_field.norm() > 1e-10) 
      wind_field /= wind_field.norm(); 

    return wind_field; 
  } 
// @sect3{The ScratchData and CopyData classes}  

// 以下对象是我们在调用 MeshWorker::mesh_loop(). 时使用的抓取和复制对象 新对象是FEInterfaceValues对象，它的工作原理类似于FEValues或FEFacesValues，只是它作用于两个单元格之间的接口，并允许我们以我们的弱形式组装接口条款。

  template <int dim> 
  struct ScratchData 
  { 
    ScratchData(const Mapping<dim> &       mapping, 
                const FiniteElement<dim> & fe, 
                const Quadrature<dim> &    quadrature, 
                const Quadrature<dim - 1> &quadrature_face, 
                const UpdateFlags          update_flags = update_values | 
                                                 update_gradients | 
                                                 update_quadrature_points | 
                                                 update_JxW_values, 
                const UpdateFlags interface_update_flags = 
                  update_values | update_gradients | update_quadrature_points | 
                  update_JxW_values | update_normal_vectors) 
      : fe_values(mapping, fe, quadrature, update_flags) 
      , fe_interface_values(mapping, 
                            fe, 
                            quadrature_face, 
                            interface_update_flags) 
    {} 

    ScratchData(const ScratchData<dim> &scratch_data) 
      : fe_values(scratch_data.fe_values.get_mapping(), 
                  scratch_data.fe_values.get_fe(), 
                  scratch_data.fe_values.get_quadrature(), 
                  scratch_data.fe_values.get_update_flags()) 
      , fe_interface_values(scratch_data.fe_interface_values.get_mapping(), 
                            scratch_data.fe_interface_values.get_fe(), 
                            scratch_data.fe_interface_values.get_quadrature(), 
                            scratch_data.fe_interface_values.get_update_flags()) 
    {} 

    FEValues<dim>          fe_values; 
    FEInterfaceValues<dim> fe_interface_values; 
  }; 

  struct CopyDataFace 
  { 
    FullMatrix<double>                   cell_matrix; 
    std::vector<types::global_dof_index> joint_dof_indices; 
  }; 

  struct CopyData 
  { 
    FullMatrix<double>                   cell_matrix; 
    Vector<double>                       cell_rhs; 
    std::vector<types::global_dof_index> local_dof_indices; 
    std::vector<CopyDataFace>            face_data; 

    template <class Iterator> 
    void reinit(const Iterator &cell, unsigned int dofs_per_cell) 
    { 
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
      cell_rhs.reinit(dofs_per_cell); 

      local_dof_indices.resize(dofs_per_cell); 
      cell->get_dof_indices(local_dof_indices); 
    } 
  }; 
// @sect3{The AdvectionProblem class}  

// 在这个准备工作之后，我们继续进行这个程序的主类，称为AdvectionProblem。

// 这对你来说应该是非常熟悉的。有趣的细节只有在实现集合函数的时候才会出现。

  template <int dim> 
  class AdvectionProblem 
  { 
  public: 
    AdvectionProblem(); 
    void run(); 

  private: 
    void setup_system(); 
    void assemble_system(); 
    void solve(); 
    void refine_grid(); 
    void output_results(const unsigned int cycle) const; 

    Triangulation<dim>   triangulation; 
    const MappingQ1<dim> mapping; 

// 此外，我们要使用DG元素。

    const FE_DGQ<dim> fe; 
    DoFHandler<dim>   dof_handler; 

    const QGauss<dim>     quadrature; 
    const QGauss<dim - 1> quadrature_face; 

// 接下来的四个成员代表要解决的线性系统。  <code>system_matrix</code> and <code>right_hand_side</code> 是由 <code>assemble_system()</code>, the <code>solution</code> 产生的，在 <code>solve()</code>. The <code>sparsity_pattern</code> 中计算，用于确定 <code>system_matrix</code> 中非零元素的位置。

    SparsityPattern      sparsity_pattern; 
    SparseMatrix<double> system_matrix; 

    Vector<double> solution; 
    Vector<double> right_hand_side; 
  }; 

// 我们从构造函数开始。 <code>fe</code> 的构造器调用中的1是多项式的度数。

  template <int dim> 
  AdvectionProblem<dim>::AdvectionProblem() 
    : mapping() 
    , fe(1) 
    , dof_handler(triangulation) 
    , quadrature(fe.tensor_degree() + 1) 
    , quadrature_face(fe.tensor_degree() + 1) 
  {} 

  template <int dim> 
  void AdvectionProblem<dim>::setup_system() 
  { 

// 在设置通常的有限元数据结构的函数中，我们首先需要分配DoF。

    dof_handler.distribute_dofs(fe); 

// 我们从生成稀疏模式开始。为此，我们首先用系统中出现的耦合物填充一个动态稀疏模式（DynamicSparsityPattern）类型的中间对象。在建立模式之后，这个对象被复制到 <code>sparsity_pattern</code> 并可以被丢弃。

// 为了建立DG离散的稀疏模式，我们可以调用类似于 DoFTools::make_sparsity_pattern, 的函数，该函数被称为 DoFTools::make_flux_sparsity_pattern:  。
    DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
    DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
    sparsity_pattern.copy_from(dsp); 

// 最后，我们设置了线性系统的所有组成部分的结构。

    system_matrix.reinit(sparsity_pattern); 
    solution.reinit(dof_handler.n_dofs()); 
    right_hand_side.reinit(dof_handler.n_dofs()); 
  } 
// @sect4{The assemble_system function}  

// 这里我们看到了与手工组装的主要区别。我们不需要在单元格和面上写循环，而是在调用 MeshWorker::mesh_loop() 时包含逻辑，我们只需要指定在每个单元格、每个边界面和每个内部面应该发生什么。这三个任务是由下面的函数里面的lambda函数处理的。

  template <int dim> 
  void AdvectionProblem<dim>::assemble_system() 
  { 
    using Iterator = typename DoFHandler<dim>::active_cell_iterator; 
    const BoundaryValues<dim> boundary_function; 

// 这是将对每个单元格执行的函数。

    const auto cell_worker = [&](const Iterator &  cell, 
                                 ScratchData<dim> &scratch_data, 
                                 CopyData &        copy_data) { 
      const unsigned int n_dofs = 
        scratch_data.fe_values.get_fe().n_dofs_per_cell(); 
      copy_data.reinit(cell, n_dofs); 
      scratch_data.fe_values.reinit(cell); 

      const auto &q_points = scratch_data.fe_values.get_quadrature_points(); 

      const FEValues<dim> &      fe_v = scratch_data.fe_values; 
      const std::vector<double> &JxW  = fe_v.get_JxW_values(); 

// 我们解决的是一个同质方程，因此在单元项中没有显示出右手。 剩下的就是整合矩阵条目。

      for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
        { 
          auto beta_q = beta(q_points[point]); 
          for (unsigned int i = 0; i < n_dofs; ++i) 
            for (unsigned int j = 0; j < n_dofs; ++j) 
              { 
                copy_data.cell_matrix(i, j) += 
                  -beta_q                      // -\beta 
                  * fe_v.shape_grad(i, point)  // \nabla \phi_i 
                  * fe_v.shape_value(j, point) // \phi_j 
                  * JxW[point];                // dx 
              } 
        } 
    }; 

// 这是为边界面调用的函数，包括使用FEFaceValues的正常积分。新的逻辑是决定该术语是进入系统矩阵（流出）还是进入右手边（流入）。

    const auto boundary_worker = [&](const Iterator &    cell, 
                                     const unsigned int &face_no, 
                                     ScratchData<dim> &  scratch_data, 
                                     CopyData &          copy_data) { 
      scratch_data.fe_interface_values.reinit(cell, face_no); 
      const FEFaceValuesBase<dim> &fe_face = 
        scratch_data.fe_interface_values.get_fe_face_values(0); 

      const auto &q_points = fe_face.get_quadrature_points(); 

      const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell(); 
      const std::vector<double> &        JxW     = fe_face.get_JxW_values(); 
      const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors(); 

      std::vector<double> g(q_points.size()); 
      boundary_function.value_list(q_points, g); 

      for (unsigned int point = 0; point < q_points.size(); ++point) 
        { 
          const double beta_dot_n = beta(q_points[point]) * normals[point]; 

          if (beta_dot_n > 0) 
            { 
              for (unsigned int i = 0; i < n_facet_dofs; ++i) 
                for (unsigned int j = 0; j < n_facet_dofs; ++j) 
                  copy_data.cell_matrix(i, j) += 
                    fe_face.shape_value(i, point)   // \phi_i 
                    * fe_face.shape_value(j, point) // \phi_j 
                    * beta_dot_n                    // \beta . n 
                    * JxW[point];                   // dx 
            } 
          else 
            for (unsigned int i = 0; i < n_facet_dofs; ++i) 
              copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i 
                                       * g[point]                     // g 
                                       * beta_dot_n  // \beta . n 
                                       * JxW[point]; // dx 
        } 
    }; 

// 这是在内部面调用的函数。参数指定了单元格、面和子面的指数（用于自适应细化）。我们只是将它们传递给FEInterfaceValues的reinit()函数。

    const auto face_worker = [&](const Iterator &    cell, 
                                 const unsigned int &f, 
                                 const unsigned int &sf, 
                                 const Iterator &    ncell, 
                                 const unsigned int &nf, 
                                 const unsigned int &nsf, 
                                 ScratchData<dim> &  scratch_data, 
                                 CopyData &          copy_data) { 
      FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values; 
      fe_iv.reinit(cell, f, sf, ncell, nf, nsf); 
      const auto &q_points = fe_iv.get_quadrature_points(); 

      copy_data.face_data.emplace_back(); 
      CopyDataFace &copy_data_face = copy_data.face_data.back(); 

      const unsigned int n_dofs        = fe_iv.n_current_interface_dofs(); 
      copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices(); 

      copy_data_face.cell_matrix.reinit(n_dofs, n_dofs); 

      const std::vector<double> &        JxW     = fe_iv.get_JxW_values(); 
      const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors(); 

      for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint) 
        { 
          const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint]; 
          for (unsigned int i = 0; i < n_dofs; ++i) 
            for (unsigned int j = 0; j < n_dofs; ++j) 
              copy_data_face.cell_matrix(i, j) +=  
                fe_iv.jump(i, qpoint) // [\phi_i] 
                * 
                fe_iv.shape_value((beta_dot_n > 0), j, qpoint) // phi_j^{upwind} 
                * beta_dot_n                                   // (\beta . n) 
                * JxW[qpoint];                                 // dx 
        }  
    }; 

// 下面的lambda函数将处理从单元格和面组件中复制数据到全局矩阵和右侧的问题。

// 虽然我们不需要AffineConstraints对象，因为在DG离散中没有悬空节点约束，但我们在这里使用一个空对象，因为这允许我们使用其`copy_local_to_global`功能。

    const AffineConstraints<double> constraints; 

    const auto copier = [&](const CopyData &c) { 
      constraints.distribute_local_to_global(c.cell_matrix, 
                                             c.cell_rhs, 
                                             c.local_dof_indices, 
                                             system_matrix, 
                                             right_hand_side); 

      for (auto &cdf : c.face_data) 
        { 
          constraints.distribute_local_to_global(cdf.cell_matrix, 
                                                 cdf.joint_dof_indices, 
                                                 system_matrix); 
        } 
    }; 

    ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face); 
    CopyData         copy_data; 

// 在这里，我们最终处理了装配问题。我们传入ScratchData和CopyData对象，以及上面的lambda函数，并指定我们要对内部面进行一次装配。

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
// @sect3{All the rest}  

// 对于这个简单的问题，我们使用了最简单的求解器，称为Richardson迭代，它代表了简单的缺陷修正。这与一个块状SSOR预处理相结合，该预处理使用DG离散化产生的系统矩阵的特殊块状结构。这些块的大小是每个单元的DoF数量。这里，我们使用SSOR预处理，因为我们没有根据流场对DoFs进行重新编号。如果在流的下游方向对DoFs进行重新编号，那么块状的Gauss-Seidel预处理（见PreconditionBlockSOR类，放松=1）会做得更好。

  template <int dim> 
  void AdvectionProblem<dim>::solve() 
  { 
    SolverControl                    solver_control(1000, 1e-12); 
    SolverRichardson<Vector<double>> solver(solver_control); 

// 这里我们创建了预处理程序。

    PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 

// 然后将矩阵分配给它，并设置正确的块大小。

    preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 

// 做完这些准备工作后，我们就可以启动线性求解器了。

    solver.solve(system_matrix, solution, right_hand_side, preconditioner); 

    std::cout << "  Solver converged in " << solver_control.last_step() 
              << " iterations." << std::endl; 
  } 

// 我们根据一个非常简单的细化标准来细化网格，即对解的梯度的近似。由于这里我们考虑的是DG(1)方法（即我们使用片状双线性形状函数），我们可以简单地计算每个单元的梯度。但是我们并不希望我们的细化指标只建立在每个单元的梯度上，而是希望同时建立在相邻单元之间的不连续解函数的跳跃上。最简单的方法是通过差分商计算近似梯度，包括考虑中的单元和其相邻的单元。这是由 <code>DerivativeApproximation</code> 类完成的，它计算近似梯度的方式类似于本教程 step-9 中描述的 <code>GradientEstimation</code> 。事实上， <code>DerivativeApproximation</code> 类是在 step-9 的 <code>GradientEstimation</code> 类之后开发的。与  step-9  中的讨论相关，这里我们考虑  $h^{1+d/2}|\nabla_h u_h|$  。此外，我们注意到，我们不考虑近似的二次导数，因为线性平流方程的解一般不在 $H^2$ 中，而只在 $H^1$ 中（或者，更准确地说：在 $H^1_\beta$ 中，即在方向 $\beta$ 上的导数是可平方整除的函数空间）。

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

//并且它们的单元格按系数 $h^{1+d/2}$ 进行缩放。
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

// 这个程序的输出包括一个自适应细化网格的vtk文件和数值解。最后，我们还用 VectorTools::integrate_difference(). 计算了解的L-无穷大规范。
  template <int dim> 
  void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
  { 
    const std::string filename = "solution-" + std::to_string(cycle) + ".vtk"; 
    std::cout << "  Writing solution to <" << filename << ">" << std::endl; 
    std::ofstream output(filename); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler);  
    data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data); 

    data_out.build_patches(mapping); 

    data_out.write_vtk(output); 

    { 
      Vector<float> values(triangulation.n_active_cells()); 
      VectorTools::integrate_difference(mapping, 
                                        dof_handler, 
                                        solution, 
                                        Functions::ZeroFunction<dim>(), 
                                        values, 
                                        quadrature, 
                                        VectorTools::Linfty_norm); 
      const double l_infty = 
        VectorTools::compute_global_error(triangulation, 
                                          values,  
                                          VectorTools::Linfty_norm); 
      std::cout << "  L-infinity norm: " << l_infty << std::endl; 
    } 
  } 

// 下面的 <code>run</code> 函数与前面的例子类似。

  template <int dim> 
  void AdvectionProblem<dim>::run() 
  { 
    for (unsigned int cycle = 0; cycle < 6; ++cycle) 
      {  
        std::cout << "Cycle " << cycle << std::endl; 

        if (cycle == 0) 
          { 
            GridGenerator::hyper_cube(triangulation); 
            triangulation.refine_global(3); 
          } 
        else 
          refine_grid(); 

        std::cout << "  Number of active cells:       " 
                  << triangulation.n_active_cells() << std::endl; 

        setup_system(); 

        std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs() 
                  << std::endl; 

        assemble_system(); 
        solve(); 

        output_results(cycle); 
      } 
  } 
} // namespace Step12 

// 下面的 <code>main</code> 函数与前面的例子也类似，不需要注释。

int main() 
{ 
  try 
    { 
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


