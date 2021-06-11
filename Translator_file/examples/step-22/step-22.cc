

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2008 - 2021 by the deal.II authors 
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2008 
 */ 


// @sect3{Include files}  

// 像往常一样，我们从包括一些著名的文件开始。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/logstream.h> 
#include <deal.II/base/function.h> 
#include <deal.II/base/utilities.h> 

#include <deal.II/lac/block_vector.h> 
#include <deal.II/lac/full_matrix.h> 
#include <deal.II/lac/block_sparse_matrix.h> 
#include <deal.II/lac/solver_cg.h> 
#include <deal.II/lac/precondition.h> 
#include <deal.II/lac/affine_constraints.h> 

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/grid_refinement.h> 

#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/dofs/dof_renumbering.h> 
#include <deal.II/dofs/dof_tools.h> 

#include <deal.II/fe/fe_q.h> 
#include <deal.II/fe/fe_system.h> 
#include <deal.II/fe/fe_values.h> 

#include <deal.II/numerics/vector_tools.h> 
#include <deal.II/numerics/matrix_tools.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/error_estimator.h> 

// 然后我们需要包括稀疏直接求解器UMFPACK的头文件。

#include <deal.II/lac/sparse_direct.h> 

// 这包括不完全LU因子化的库，它将被用作3D的预处理程序。

#include <deal.II/lac/sparse_ilu.h> 

// 这是C++语言。

#include <iostream> 
#include <fstream> 
#include <memory> 

// 和所有的程序一样，名字空间dealii被包括在内。

namespace Step22 
{ 
  using namespace dealii; 
// @sect3{Defining the inner preconditioner type}  

// 正如介绍中所解释的，我们将分别对两个和三个空间维度使用不同的预处理程序。我们通过使用空间维度作为模板参数来区分它们。关于模板的细节，请参见 step-4 。我们不打算在这里创建任何预处理对象，我们所做的只是创建一个持有确定预处理类的本地别名的类，这样我们就可以以独立于维度的方式编写我们的程序。

  template <int dim> 
  struct InnerPreconditioner; 

// 在二维中，我们将使用一个稀疏的直接求解器作为预处理程序。

  template <> 
  struct InnerPreconditioner<2> 
  { 
    using type = SparseDirectUMFPACK; 
  }; 

// 还有三维的ILU预处理，由SparseILU调用。

  template <> 
  struct InnerPreconditioner<3> 
  { 
    using type = SparseILU<double>; 
  }; 
// @sect3{The <code>StokesProblem</code> class template}  

// 这是对 step-20 的改编，所以主类和数据类型与那里使用的几乎相同。唯一不同的是，我们有一个额外的成员  <code>preconditioner_matrix</code>  ，用于预处理Schur补码，以及一个相应的稀疏模式  <code>preconditioner_sparsity_pattern</code>  。此外，我们没有依赖LinearOperator，而是实现了我们自己的InverseMatrix类。

// 在这个例子中，我们还使用了自适应网格细化，其处理方式与  step-6  类似。根据介绍中的讨论，我们也将使用AffineConstraints对象来实现Dirichlet边界条件。因此，我们改变名称  <code>hanging_node_constraints</code> into <code>constraints</code>  。

  template <int dim> 
  class StokesProblem 
  { 
  public: 
    StokesProblem(const unsigned int degree); 
    void run(); 

  private: 
    void setup_dofs(); 
    void assemble_system(); 
    void solve(); 
    void output_results(const unsigned int refinement_cycle) const; 
    void refine_mesh(); 

    const unsigned int degree; 

    Triangulation<dim> triangulation; 
    FESystem<dim>      fe; 
    DoFHandler<dim>    dof_handler; 

    AffineConstraints<double> constraints; 

    BlockSparsityPattern      sparsity_pattern; 
    BlockSparseMatrix<double> system_matrix; 

    BlockSparsityPattern      preconditioner_sparsity_pattern; 
    BlockSparseMatrix<double> preconditioner_matrix; 

    BlockVector<double> solution; 
    BlockVector<double> system_rhs; 

// 这一条是新的：我们将使用一个所谓的共享指针结构来访问预处理程序。共享指针本质上只是指针的一种方便形式。几个共享指针可以指向同一个对象（就像普通的指针一样），但是当最后一个指向前提器对象的共享指针对象被删除时（例如共享指针对象超出了范围，它所在的类被销毁，或者指针被分配给了不同的前提器对象），那么指向的前提器对象也被销毁。这确保了我们不必手动跟踪有多少地方仍在引用一个前置条件器对象，它永远不会产生内存泄漏，也不会产生一个指向已被销毁对象的悬空指针。

    std::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner; 
  }; 
// @sect3{Boundary values and right hand side}  

// 与 step-20 和其他大多数例子程序一样，下一个任务是定义PDE的数据：对于斯托克斯问题，我们将在部分边界上使用自然边界值（即同质诺伊曼型），对于这些边界，我们不必做任何特殊处理（同质性意味着弱形式中的相应项只是零），而在边界的其余部分使用速度的边界条件（迪里希勒型），如介绍中所述。

// 为了强制执行速度上的Dirichlet边界值，我们将像往常一样使用 VectorTools::interpolate_boundary_values 函数，这要求我们写一个具有与有限元一样多分量的函数对象。换句话说，我们必须在 $(u,p)$ -空间上定义函数，但在插值边界值时，我们要过滤掉压力分量。

// 下面的函数对象是介绍中描述的边界值的表示。

  template <int dim> 
  class BoundaryValues : public Function<dim> 
  { 
  public: 
    BoundaryValues() 
      : Function<dim>(dim + 1) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  double BoundaryValues<dim>::value(const Point<dim> & p, 
                                    const unsigned int component) const 
  { 
    Assert(component < this->n_components, 
           ExcIndexRange(component, 0, this->n_components)); 

    if (component == 0) 
      return (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0)); 
    return 0; 
  } 

  template <int dim> 
  void BoundaryValues<dim>::vector_value(const Point<dim> &p, 
                                         Vector<double> &  values) const 
  { 
    for (unsigned int c = 0; c < this->n_components; ++c) 
      values(c) = BoundaryValues<dim>::value(p, c); 
  } 

// 我们为右手边实现类似的函数，在目前的例子中，右手边只是零。

  template <int dim> 
  class RightHandSide : public Function<dim> 
  { 
  public: 
    RightHandSide() 
      : Function<dim>(dim + 1) 
    {} 

    virtual double value(const Point<dim> & p, 
                         const unsigned int component = 0) const override; 

    virtual void vector_value(const Point<dim> &p, 
                              Vector<double> &  value) const override; 
  }; 

  template <int dim> 
  double RightHandSide<dim>::value(const Point<dim> & /*p*/, 
                                   const unsigned int /*component*/) const 
  { 
    return 0; 
  } 

  template <int dim> 
  void RightHandSide<dim>::vector_value(const Point<dim> &p, 
                                        Vector<double> &  values) const 
  { 
    for (unsigned int c = 0; c < this->n_components; ++c) 
      values(c) = RightHandSide<dim>::value(p, c); 
  } 
// @sect3{Linear solvers and preconditioners}  

// 在介绍中广泛讨论了线性求解器和预处理器。在这里，我们创建将被使用的各自对象。

//  @sect4{The <code>InverseMatrix</code> class template}   <code>InverseMatrix</code> 类表示逆矩阵的数据结构。与 step-20 不同，我们用一个类来实现，而不是用辅助函数inverse_linear_operator()，我们将把这个类应用于不同种类的矩阵，这些矩阵需要不同的预处理程序（在 step-20 中，我们只对质量矩阵使用非同一性预处理程序）。矩阵和预处理器的类型通过模板参数传递给这个类，当创建 <code>InverseMatrix</code> 对象时，这些类型的矩阵和预处理器对象将被传递给构造器。成员函数 <code>vmult</code> 是通过解决一个线性系统得到的。

  template <class MatrixType, class PreconditionerType> 
  class InverseMatrix : public Subscriptor 
  { 
  public: 
    InverseMatrix(const MatrixType &        m, 
                  const PreconditionerType &preconditioner); 

    void vmult(Vector<double> &dst, const Vector<double> &src) const; 

  private: 
    const SmartPointer<const MatrixType>         matrix; 
    const SmartPointer<const PreconditionerType> preconditioner; 
  }; 

  template <class MatrixType, class PreconditionerType> 
  InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
    const MatrixType &        m, 
    const PreconditionerType &preconditioner) 
    : matrix(&m) 
    , preconditioner(&preconditioner) 
  {} 

// 这就是 <code>vmult</code> 函数的实现。

// 在这个类中，我们对解算器控制使用了一个相当大的容忍度。这样做的原因是，该函数被频繁使用，因此，任何使CG求解中的残差变小的额外努力都会使求解更加昂贵。请注意，我们不仅将该类作为Schur补码的预处理程序，而且在形成拉普拉斯矩阵的逆时也使用该类；因此，该类直接对解本身的精度负责，所以我们也不能选择太大的公差。

  template <class MatrixType, class PreconditionerType> 
  void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
    Vector<double> &      dst, 
    const Vector<double> &src) const 
  { 
    SolverControl            solver_control(src.size(), 1e-6 * src.l2_norm()); 
    SolverCG<Vector<double>> cg(solver_control); 

    dst = 0; 

    cg.solve(*matrix, dst, src, *preconditioner); 
  } 
// @sect4{The <code>SchurComplement</code> class template}  

// 这个类实现了介绍中讨论的Schur补码。它与  step-20  相类似。 不过，我们现在用一个模板参数 <code>PreconditionerType</code> 来调用它，以便在指定逆矩阵类的各自类型时访问它。作为上述定义的结果，声明  <code>InverseMatrix</code>  现在包含了上述预处理类的第二个模板参数，这也影响到  <code>SmartPointer</code> object <code>m_inverse</code>  。

  template <class PreconditionerType> 
  class SchurComplement : public Subscriptor 
  { 
  public: 
    SchurComplement( 
      const BlockSparseMatrix<double> &system_matrix, 
      const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse); 

    void vmult(Vector<double> &dst, const Vector<double> &src) const; 

  private: 
    const SmartPointer<const BlockSparseMatrix<double>> system_matrix; 
    const SmartPointer< 
      const InverseMatrix<SparseMatrix<double>, PreconditionerType>> 
      A_inverse; 

    mutable Vector<double> tmp1, tmp2; 
  }; 

  template <class PreconditionerType> 
  SchurComplement<PreconditionerType>::SchurComplement( 
    const BlockSparseMatrix<double> &system_matrix, 
    const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse) 
    : system_matrix(&system_matrix) 
    , A_inverse(&A_inverse) 
    , tmp1(system_matrix.block(0, 0).m()) 
    , tmp2(system_matrix.block(0, 0).m()) 
  {} 

  template <class PreconditionerType> 
  void 
  SchurComplement<PreconditionerType>::vmult(Vector<double> &      dst, 
                                             const Vector<double> &src) const 
  { 
    system_matrix->block(0, 1).vmult(tmp1, src); 
    A_inverse->vmult(tmp2, tmp1); 
    system_matrix->block(1, 0).vmult(dst, tmp2); 
  } 
// @sect3{StokesProblem class implementation}  
// @sect4{StokesProblem::StokesProblem}  

// 这个类的构造函数看起来与  step-20  的构造函数非常相似。构造函数初始化了多项式程度、三角形、有限元系统和dof处理器的变量。矢量速度分量的基础多项式函数的阶数为 <code>degree+1</code> ，压力的阶数为 <code>degree</code> 。 这就得到了LBB稳定元对 $Q_{degree+1}^d\times Q_{degree}$ ，通常被称为泰勒-霍德元。

// 请注意，我们用MeshSmoothing参数初始化三角形，这可以确保单元的细化是以PDE解的近似保持良好的方式进行的（如果网格过于非结构化就会出现问题），详情请参见 <code>Triangulation::MeshSmoothing</code> 的文档。

  template <int dim> 
  StokesProblem<dim>::StokesProblem(const unsigned int degree) 
    : degree(degree) 
    , triangulation(Triangulation<dim>::maximum_smoothing) 
    , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1) 
    , dof_handler(triangulation) 
  {} 
// @sect4{StokesProblem::setup_dofs}  

// 给定一个网格，该函数将自由度与之关联，并创建相应的矩阵和向量。在开始的时候，它还释放了指向预处理对象的指针（如果共享指针在此时指向任何东西的话），因为在这之后肯定不会再需要它了，在组装矩阵后必须重新计算，并将稀疏矩阵从其稀疏模式对象中解开。

// 然后，我们继续分配自由度并重新编号。为了使ILU预处理程序（在3D中）有效地工作，重要的是以这样的方式列举自由度，以减少矩阵的带宽，或者也许更重要的是：以这样的方式使ILU尽可能地接近于真正的LU分解。另一方面，我们需要保留在  step-20  和  step-21  中已经看到的速度和压力的块状结构。这将分两步完成。首先，对所有的道次进行重新编号，以改善ILU，然后我们再一次按组件重新编号。由于 <code>DoFRenumbering::component_wise</code> 没有触及单个块内的重新编号，所以第一步的基本重新编号仍然存在。至于如何对自由度进行重新编号以提高ILU：deal.II有许多算法试图找到排序以提高ILU，或减少矩阵的带宽，或优化其他方面。DoFRenumbering命名空间显示了我们在本教程程序中基于这里讨论的测试案例而获得的几种算法的结果比较。在这里，我们将使用传统的Cuthill-McKee算法，该算法已经在之前的一些教程程序中使用。 在<a href="#improved-ilu">section on improved ILU</a>中我们将更详细地讨论这个问题。
//与以前的教程程序相比，
//还有一个变化。没有理由对 <code>dim</code> 的速度成分进行单独排序。事实上，与其先列举所有 $x$ -velocities，再列举所有 $y$ -velocities，等等，我们希望将所有速度放在一起，只在速度（所有分量）和压力之间分开。默认情况下， DoFRenumbering::component_wise 函数不是这样做的：它把每个矢量分量分开处理；我们要做的是把几个分量分成 "块"，并把这个块结构传递给该函数。因此，我们分配一个矢量 <code>block_component</code> ，有多少个元素就有多少个分量，描述所有的速度分量对应于块0，而压力分量将形成块1。

  template <int dim> 
  void StokesProblem<dim>::setup_dofs() 
  { 
    A_preconditioner.reset(); 
    system_matrix.clear(); 
    preconditioner_matrix.clear(); 

    dof_handler.distribute_dofs(fe); 
    DoFRenumbering::Cuthill_McKee(dof_handler); 

    std::vector<unsigned int> block_component(dim + 1, 0); 
    block_component[dim] = 1; 
    DoFRenumbering::component_wise(dof_handler, block_component); 

// 现在是对Dirichlet边界条件的实现，在介绍中的讨论之后，这应该是很明显的。所有的变化是，这个函数已经出现在设置函数中，而我们习惯于在一些汇编例程中看到它。在我们设置网格的下面，我们将把施加Dirichlet边界条件的顶部边界与边界指标1联系起来。 我们必须将这个边界指标作为第二个参数传递给下面的插值函数。 不过，还有一件事。 描述Dirichlet条件的函数是为所有分量定义的，包括速度和压力。然而，Dirichlet条件只为速度而设置。 为此，我们使用一个只选择速度分量的ComponentMask。通过指定我们想要的特定分量，从有限元中获得该分量掩码。由于我们使用自适应细化网格，仿生约束对象需要首先填充由DoF处理程序生成的悬挂节点约束。注意这两个函数的顺序；我们首先计算悬挂节点约束，然后将边界值插入约束对象。这确保了我们在有悬挂节点的边界上尊重H<sup>1</sup>一致性（在三个空间维度上），悬挂节点需要支配Dirichlet边界值。

    { 
      constraints.clear(); 

      FEValuesExtractors::Vector velocities(0); 
      DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
      VectorTools::interpolate_boundary_values(dof_handler, 
                                               1, 
                                               BoundaryValues<dim>(), 
                                               constraints, 
                                               fe.component_mask(velocities)); 
    } 

    constraints.close(); 

// 与 step-20 相类似，我们计算各个组件中的道夫。我们可以用与那里相同的方式来做，但我们想在我们已经用于重新编号的块结构上进行操作。函数  <code>DoFTools::count_dofs_per_fe_block</code>  的作用与  <code>DoFTools::count_dofs_per_fe_component</code>  相同，但现在通过  <code>block_component</code>  将速度和压力块分组。

    const std::vector<types::global_dof_index> dofs_per_block = 
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
    const unsigned int n_u = dofs_per_block[0]; 
    const unsigned int n_p = dofs_per_block[1]; 

    std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
              << std::endl 
              << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
              << " (" << n_u << '+' << n_p << ')' << std::endl; 

// 下一个任务是为我们将创建的系统矩阵分配一个稀疏模式，为预处理矩阵分配一个稀疏模式。我们可以用与 step-20 相同的方式来做这件事，即通过 DoFTools::make_sparsity_pattern. 直接建立一个SparsityPattern类型的对象，但是，有一个重要的理由不这样做。在3D中，函数 DoFTools::max_couplings_between_dofs 对各个道夫之间的耦合产生了一个保守但相当大的数字，因此，最初为创建矩阵的稀疏模式提供的内存太多--实际上，对于中等大小的3D问题，初始稀疏模式甚至无法放入大多数系统的物理内存中，也请参见 step-18  中的讨论。相反，我们首先建立临时对象，使用不同的数据结构，不需要分配更多的内存，但不适合作为SparseMatrix或BlockSparseMatrix对象的基础；在第二步，我们将这些对象复制到BlockSparsityPattern类型的对象中。这完全类似于我们在  step-11  和  step-18  中已经做过的事情。特别是，我们利用了这样一个事实，即我们永远不会写入系统矩阵的 $(1,1)$ 块中，而且这是唯一需要填充的预处理矩阵块。

// 所有这些都是在新范围内完成的，这意味着一旦信息被复制到  <code>sparsity_pattern</code>  ，  <code>dsp</code>  的内存将被释放。

    { 
      BlockDynamicSparsityPattern dsp(2, 2); 

      dsp.block(0, 0).reinit(n_u, n_u); 
      dsp.block(1, 0).reinit(n_p, n_u); 
      dsp.block(0, 1).reinit(n_u, n_p); 
      dsp.block(1, 1).reinit(n_p, n_p); 

      dsp.collect_sizes(); 

      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 

      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (!((c == dim) && (d == dim))) 
            coupling[c][d] = DoFTools::always; 
          else 
            coupling[c][d] = DoFTools::none; 

      DoFTools::make_sparsity_pattern( 
        dof_handler, coupling, dsp, constraints, false); 

 
    } 

    { 
      BlockDynamicSparsityPattern preconditioner_dsp(2, 2); 

      preconditioner_dsp.block(0, 0).reinit(n_u, n_u); 
      preconditioner_dsp.block(1, 0).reinit(n_p, n_u); 
      preconditioner_dsp.block(0, 1).reinit(n_u, n_p); 
      preconditioner_dsp.block(1, 1).reinit(n_p, n_p); 

      preconditioner_dsp.collect_sizes(); 

      Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1); 

      for (unsigned int c = 0; c < dim + 1; ++c) 
        for (unsigned int d = 0; d < dim + 1; ++d) 
          if (((c == dim) && (d == dim))) 
            preconditioner_coupling[c][d] = DoFTools::always; 
          else 
            preconditioner_coupling[c][d] = DoFTools::none; 

      DoFTools::make_sparsity_pattern(dof_handler, 
                                      preconditioner_coupling, 
                                      preconditioner_dsp, 
                                      constraints, 
                                      false); 

      preconditioner_sparsity_pattern.copy_from(preconditioner_dsp); 
    } 

// 最后，与  step-20  中的方法类似，从块状结构中创建系统矩阵、前导矩阵、解决方案和右侧向量。

    system_matrix.reinit(sparsity_pattern); 
    preconditioner_matrix.reinit(preconditioner_sparsity_pattern); 

    solution.reinit(2); 
    solution.block(0).reinit(n_u); 
    solution.block(1).reinit(n_p); 
    solution.collect_sizes(); 

    system_rhs.reinit(2); 
    system_rhs.block(0).reinit(n_u); 
    system_rhs.block(1).reinit(n_p); 
    system_rhs.collect_sizes(); 
  } 
// @sect4{StokesProblem::assemble_system}  

// 汇编过程遵循 step-20 和介绍中的讨论。我们使用众所周知的缩写来表示保存本单元自由度的局部矩阵、右手边和全局编号的数据结构。

  template <int dim> 
  void StokesProblem<dim>::assemble_system() 
  { 
    system_matrix         = 0; 
    system_rhs            = 0; 
    preconditioner_matrix = 0; 

    QGauss<dim> quadrature_formula(degree + 2); 

    FEValues<dim> fe_values(fe, 
                            quadrature_formula, 
                            update_values | update_quadrature_points | 
                              update_JxW_values | update_gradients); 

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 

    const unsigned int n_q_points = quadrature_formula.size(); 

    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
    FullMatrix<double> local_preconditioner_matrix(dofs_per_cell, 
                                                   dofs_per_cell); 
    Vector<double>     local_rhs(dofs_per_cell); 

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 

    const RightHandSide<dim>    right_hand_side; 
    std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 

// 接下来，我们需要两个对象，作为FEValues对象的提取器。它们的用途在  @ref  vector_valued 的报告中详细解释。

    const FEValuesExtractors::Vector velocities(0); 
    const FEValuesExtractors::Scalar pressure(dim); 

// 作为对 step-20 和 step-21 的扩展，我们包括了一些优化，使这个特定问题的装配速度大大加快。这些改进是基于这样的观察：当我们像 step-20 那样做时，我们做了太多次的计算：对称梯度实际上在每个正交点有 <code>dofs_per_cell</code> 个不同的值，但是我们从FEValues对象中提取了 <code>dofs_per_cell*dofs_per_cell</code> 次。

// - 在 <code>i</code> 的循环和 <code>j</code> 的内循环中。在3D中，这意味着评估它 $89^2=7921$ 次而不是 $89$ 次，这是一个不小的差别。

// 所以我们在这里要做的是，在开始对单元上的道夫进行循环之前，在正交点得到一个秩-2张量的向量（类似的还有压力上的发散和基函数值）来避免这种重复计算。首先，我们创建各自的对象来保存这些值。然后，我们开始在所有单元上进行循环，并在正交点上进行循环，在那里我们首先提取这些值。我们在这里还实现了一个优化：本地矩阵（以及全局矩阵）将是对称的，因为所有涉及的操作都是相对于 $i$ 和 $j$ 对称的。这可以通过简单地运行内循环而不是 <code>dofs_per_cell</code>, but only up to <code>i</code> 来实现，即外循环的索引。

    std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell); 
    std::vector<double>                  div_phi_u(dofs_per_cell); 
    std::vector<double>                  phi_p(dofs_per_cell); 

    for (const auto &cell : dof_handler.active_cell_iterators()) 
      { 
        fe_values.reinit(cell); 
        local_matrix                = 0; 
        local_preconditioner_matrix = 0; 
        local_rhs                   = 0; 

        right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
                                          rhs_values); 

        for (unsigned int q = 0; q < n_q_points; ++q) 
          { 
            for (unsigned int k = 0; k < dofs_per_cell; ++k) 
              { 
                symgrad_phi_u[k] = 
                  fe_values[velocities].symmetric_gradient(k, q); 
                div_phi_u[k] = fe_values[velocities].divergence(k, q); 
                phi_p[k]     = fe_values[pressure].value(k, q); 
              } 

// 最后是系统矩阵和我们用于预处理程序的矩阵的双线性形式。回顾一下，这两个的公式分别是
  //  @f{align*}{
  //    A_{ij} &= a(\varphi_i,\varphi_j)
  //    \\     &= \underbrace{2(\varepsilon(\varphi_{i,\textbf{u}}),
  //                            \varepsilon(\varphi_{j,\textbf{u}}))_{\Omega}}
  //                         _{(1)}
  //            \;
  //              \underbrace{- (\textrm{div}\; \varphi_{i,\textbf{u}},
  //                             \varphi_{j,p})_{\Omega}}
  //                         _{(2)}
  //            \;
  //              \underbrace{- (\varphi_{i,p},
  //                             \textrm{div}\;
  //                             \varphi_{j,\textbf{u}})_{\Omega}}
  //                         _{(3)}
  //  @f}
  //  和
  //  @f{align*}{
  //    M_{ij} &= \underbrace{(\varphi_{i,p},
  //                           \varphi_{j,p})_{\Omega}}
  //                         _{(4)},
  //  @f} ， 
  //  其中 $\varphi_{i,\textbf{u}}$ 和 $\varphi_{i,p}$ 是 $i$ th形状函数的速度和压力成分。然后，上述各种术语在下面的实现中很容易识别。

            for (unsigned int i = 0; i < dofs_per_cell; ++i) 
              { 
                for (unsigned int j = 0; j <= i; ++j) 
                  { 
                    local_matrix(i, j) += 
                      (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) // (1) 
                       - div_phi_u[i] * phi_p[j]                 // (2) 
                       - phi_p[i] * div_phi_u[j])                // (3) 
                      * fe_values.JxW(q);                        // * dx 

                    local_preconditioner_matrix(i, j) += 
                      (phi_p[i] * phi_p[j]) // (4) 
                      * fe_values.JxW(q);   // * dx 
                  } 

// 注意在上述（1）的实现中，`operator*`被重载用于对称张量，产生两个张量之间的标量乘积。            对于右手边，我们利用形状函数只在一个分量中不为零的事实（因为我们的元素是原始的）。 我们不是将代表形状函数i的dim+1值的张量与整个右手边的向量相乘，而是只看唯一的非零分量。函数 FiniteElement::system_to_component_index 将返回这个形状函数所处的分量（0=x速度，1=y速度，2=2d中的压力），我们用它来挑选出右手边向量的正确分量来相乘。

                const unsigned int component_i = 
                  fe.system_to_component_index(i).first; 
                local_rhs(i) += (fe_values.shape_value(i, q)   // (phi_u_i(x_q) 
                                 * rhs_values[q](component_i)) // * f(x_q)) 
                                * fe_values.JxW(q);            // * dx 
              } 
          } 

// 在我们将局部数据写入全局矩阵之前（同时使用AffineConstraints对象来应用Dirichlet边界条件并消除悬挂的节点约束，正如我们在介绍中讨论的那样），我们必须注意一件事。由于对称性，我们只建立了一半的局部矩阵，但我们要保存完整的矩阵，以便使用标准函数进行解算。这是通过翻转指数来实现的，以防我们指向本地矩阵的空部分。

        for (unsigned int i = 0; i < dofs_per_cell; ++i) 
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
            { 
              local_matrix(i, j) = local_matrix(j, i); 
              local_preconditioner_matrix(i, j) = 
                local_preconditioner_matrix(j, i); 
            } 

        cell->get_dof_indices(local_dof_indices); 
        constraints.distribute_local_to_global(local_matrix, 
                                               local_rhs, 
                                               local_dof_indices, 
                                               system_matrix, 
                                               system_rhs); 
        constraints.distribute_local_to_global(local_preconditioner_matrix, 
                                               local_dof_indices, 
                                               preconditioner_matrix); 
      } 

// 在我们要解决这个线性系统之前，我们为速度-速度矩阵生成一个预处理程序，即系统矩阵中的 <code>block(0,0)</code> 。如上所述，这取决于空间维度。由于 <code>InnerPreconditioner::type</code> 别名所描述的两个类具有相同的接口，因此无论我们想使用稀疏直接求解器还是ILU，都不需要做任何不同的事情。

    std::cout << "   Computing preconditioner..." << std::endl << std::flush; 

    A_preconditioner = 
      std::make_shared<typename InnerPreconditioner<dim>::type>(); 
    A_preconditioner->initialize( 
      system_matrix.block(0, 0), 
      typename InnerPreconditioner<dim>::type::AdditionalData()); 
  } 

//  @sect4{StokesProblem::solve}  

// 经过前面介绍中的讨论和各自类的定义， <code>solve</code> 函数的实现是相当直接的，其方式与 step-20 类似。首先，我们需要一个 <code>InverseMatrix</code> 类的对象，代表矩阵A的逆。正如在介绍中所描述的，在  <code>InnerPreconditioner::type</code>  类型的内部预处理器的帮助下，生成了逆。

  template <int dim> 
  void StokesProblem<dim>::solve() 
  { 
    const InverseMatrix<SparseMatrix<double>, 
                        typename InnerPreconditioner<dim>::type> 
                   A_inverse(system_matrix.block(0, 0), *A_preconditioner); 
    Vector<double> tmp(solution.block(0).size()); 

// 这与  step-20  中的情况一样。我们生成 Schur 补数的右手边  $B A^{-1} F - G$  和一个代表各自线性运算的对象  $B A^{-1} B^T$  ，现在有一个模板参数表示预处理器

// - 按照类的定义。

    { 
      Vector<double> schur_rhs(solution.block(1).size()); 
      A_inverse.vmult(tmp, system_rhs.block(0)); 
      system_matrix.block(1, 0).vmult(schur_rhs, tmp); 
      schur_rhs -= system_rhs.block(1); 

      SchurComplement<typename InnerPreconditioner<dim>::type> schur_complement( 
        system_matrix, A_inverse); 

// 解算器调用的常规控制结构被创建...

      SolverControl            solver_control(solution.block(1).size(), 
                                   1e-6 * schur_rhs.l2_norm()); 
      SolverCG<Vector<double>> cg(solver_control); 

// 现在是对舒尔补码的预处理。正如介绍中所解释的，预处理是由压力变量的质量矩阵来完成的。

// 实际上，求解器需要有 $P^{-1}$ 形式的预处理，所以我们需要创建一个逆运算。我们再次使用一个 <code>InverseMatrix</code> 类的对象，它实现了求解器需要的 <code>vmult</code> 操作。 在这种情况下，我们必须对压力质量矩阵进行反转。正如在早期的教程程序中已经证明的那样，质量矩阵的反转是一个相当便宜和简单的操作（与拉普拉斯矩阵等相比）。带有ILU预处理的CG方法在5-10步内收敛，与网格大小无关。 这正是我们在这里所做的。我们选择另一个ILU预处理，并通过相应的模板参数将其带入InverseMatrix对象。 然后在逆矩阵的vmult操作中调用一个CG求解器。

// 另一种方法是选择因子为1.2的SSOR预处理器，这种方法构建成本较低，但之后需要更多的迭代。它需要大约两倍的迭代次数，但其生成的成本几乎可以忽略不计。

      SparseILU<double> preconditioner; 
      preconditioner.initialize(preconditioner_matrix.block(1, 1), 
                                SparseILU<double>::AdditionalData()); 

      InverseMatrix<SparseMatrix<double>, SparseILU<double>> m_inverse( 
        preconditioner_matrix.block(1, 1), preconditioner); 

// 有了舒尔补码和高效的预处理程序，我们可以用通常的方法解决压力的相关方程（即解向量中的0块）。

      cg.solve(schur_complement, solution.block(1), schur_rhs, m_inverse); 

// 在这第一个求解步骤之后，必须将悬挂的节点约束分布到求解中，以实现一致的压力场。

      constraints.distribute(solution); 

      std::cout << "  " << solver_control.last_step() 
                << " outer CG Schur complement iterations for pressure" 
                << std::endl; 
    } 

// 和 step-20 一样，我们最后需要解速度方程，在这里我们插入压力方程的解。这只涉及我们已经知道的对象

// 所以我们只需用 $p$ 乘以 $B^T$ ，减去右边的部分，再乘以 $A$ 的逆数。最后，我们需要分配悬挂节点的约束，以获得一个一致的流场。

    { 
      system_matrix.block(0, 1).vmult(tmp, solution.block(1)); 
      tmp *= -1; 
      tmp += system_rhs.block(0); 

      A_inverse.vmult(solution.block(0), tmp); 

      constraints.distribute(solution); 
    } 
  } 
// @sect4{StokesProblem::output_results}  

// 下一个函数生成图形输出。在这个例子中，我们将使用VTK文件格式。 我们给问题中的各个变量附上名字： <code>velocity</code> to the <code>dim</code> 速度的组成部分和 <code>pressure</code> 压力的组成部分。

// 并非所有的可视化程序都有能力将各个矢量分量组合成一个矢量来提供矢量图；特别是对于一些基于VTK的可视化程序来说，这一点是成立的。在这种情况下，在包含数据的文件中应该已经描述了组件的逻辑分组为矢量的情况。换句话说，我们需要做的是为我们的输出编写者提供一种方法，让他们知道有限元的哪些分量在逻辑上形成一个矢量（在 $d$ 空间维度上有 $d$ 分量），而不是让他们假设我们只是有一堆标量场。 这是用 <code>DataComponentInterpretation</code> 命名空间的成员实现的：和文件名一样，我们创建一个矢量，其中第一个 <code>dim</code> 分量指的是速度，并被赋予 DataComponentInterpretation::component_is_part_of_vector; 标签，我们最后推一个标签 DataComponentInterpretation::component_is_scalar 来描述压力变量的分组。

// 然后函数的其余部分与  step-20  中的相同。

  template <int dim> 
  void 
  StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const 
  { 
    std::vector<std::string> solution_names(dim, "velocity"); 
    solution_names.emplace_back("pressure"); 

    std::vector<DataComponentInterpretation::DataComponentInterpretation> 
      data_component_interpretation( 
        dim, DataComponentInterpretation::component_is_part_of_vector); 
    data_component_interpretation.push_back( 
      DataComponentInterpretation::component_is_scalar); 

    DataOut<dim> data_out; 
    data_out.attach_dof_handler(dof_handler); 
    data_out.add_data_vector(solution, 
                             solution_names, 
                             DataOut<dim>::type_dof_data, 
                             data_component_interpretation); 
    data_out.build_patches(); 

    std::ofstream output( 
      "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk"); 
    data_out.write_vtk(output); 
  } 
// @sect4{StokesProblem::refine_mesh}  

// 这是 <code>StokesProblem</code> 类中最后一个有趣的函数。 正如它的名字所示，它获取问题的解决方案，并在需要时细化网格。其过程与 step-6 中的相应步骤相同，不同的是我们只根据压力的变化进行细化，也就是说，我们用ComponentMask类型的掩码对象调用Kelly误差估计器，选择我们感兴趣的压力的单一标量分量（我们通过指定我们想要的分量从有限元类中得到这样一个掩码）。此外，我们没有再次粗化网格。

  template <int dim> 
  void StokesProblem<dim>::refine_mesh() 
  { 
    Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 

    FEValuesExtractors::Scalar pressure(dim); 
    KellyErrorEstimator<dim>::estimate( 
      dof_handler, 
      QGauss<dim - 1>(degree + 1), 
      std::map<types::boundary_id, const Function<dim> *>(), 
      solution, 
      estimated_error_per_cell, 
      fe.component_mask(pressure)); 

    GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
                                                    estimated_error_per_cell, 
                                                    0.3, 
                                                    0.0); 
    triangulation.execute_coarsening_and_refinement(); 
  } 
// @sect4{StokesProblem::run}  

// 在斯托克斯类中的最后一步，像往常一样，是生成初始网格的函数，并按各自的顺序调用其他函数。

// 我们从一个大小为 $4 \times 1$ （2D）或 $4 \times 1 \times 1$ （3D）的矩形开始，在 $R^2/R^3$ 中分别放置为 $(-2,2)\times(-1,0)$ 或 $(-2,2)\times(0,1)\times(-1,0)$  。在每个方向上以相等的网格大小开始是很自然的，所以我们在第一个坐标方向上将初始矩形细分四次。为了将创建网格所涉及的变量的范围限制在我们实际需要的范围内，我们将整个块放在一对大括号之间。

  template <int dim> 
  void StokesProblem<dim>::run() 
  { 
    { 
      std::vector<unsigned int> subdivisions(dim, 1); 
      subdivisions[0] = 4; 

      const Point<dim> bottom_left = (dim == 2 ?                // 
                                        Point<dim>(-2, -1) :    // 2d case 
                                        Point<dim>(-2, 0, -1)); // 3d case 

      const Point<dim> top_right = (dim == 2 ?              // 
                                      Point<dim>(2, 0) :    // 2d case 
                                      Point<dim>(2, 1, 0)); // 3d case 

      GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                                subdivisions, 
                                                bottom_left, 
                                                top_right); 
    } 

// 边界指标1被设置为所有受Dirichlet边界条件约束的边界，即位于最后一个坐标方向上的0的面。详见上面的例子描述。

    for (const auto &cell : triangulation.active_cell_iterators()) 
      for (const auto &face : cell->face_iterators()) 
        if (face->center()[dim - 1] == 0) 
          face->set_all_boundary_ids(1); 

// 然后，在第一次求解之前，我们应用一个初始细化。在3D中，会有更多的自由度，所以我们在那里细化得更少。

    triangulation.refine_global(4 - dim); 

// 正如在 step-6 中第一次看到的那样，我们在不同的细化级别上循环细化（除了第一个循环），设置自由度和矩阵，组装，求解和创建输出。

    for (unsigned int refinement_cycle = 0; refinement_cycle < 6; 
         ++refinement_cycle) 
      { 
        std::cout << "Refinement cycle " << refinement_cycle << std::endl; 

        if (refinement_cycle > 0) 
          refine_mesh(); 

        setup_dofs(); 

        std::cout << "   Assembling..." << std::endl << std::flush; 
        assemble_system(); 

        std::cout << "   Solving..." << std::flush; 
        solve(); 

        output_results(refinement_cycle); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step22 
// @sect3{The <code>main</code> function}  

// 主函数与  step-20  中的相同。我们将元素度数作为参数传递，并在众所周知的模板槽中选择空间尺寸。

int main() 
{ 
  try 
    { 
      using namespace Step22; 

      StokesProblem<2> flow_problem(1); 
      flow_problem.run(); 
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


