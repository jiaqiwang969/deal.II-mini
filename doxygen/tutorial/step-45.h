/**
@page step_45 The step-45 tutorial program
This tutorial depends on step-6.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Settingupperiodicityconstraintsondistributedtriangulations">Setting up periodicity constraints on distributed triangulations</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-45/doc/intro.dox

 <br> 

［<i>This program was contributed by Daniel Arndt and Matthias Maier.</i>］ ［<a name="Intro"></a>］

<a name="Introduction"></a><h1>Introduction</h1>


在这个例子中，我们介绍了如何在deal.II中使用周期性边界条件。周期性边界条件是代数约束，通常出现在一个大域的代表性区域的计算中，这些区域在一个或多个方向上重复。

一个例子是模拟光子晶体的电子结构，因为它们有一个类似格子的结构，因此，往往只需要在格子的一个盒子上进行实际计算。为了能够以这种方式进行，我们必须假设该模型可以周期性地扩展到其他盒子上；这要求解决方案具有周期性结构。

<a name="Procedure"></a>

<a name="Procedure"></a><h1>Procedure</h1>


deal.II提供了一些高水平的入口来施加周期性边界条件。应用周期性边界条件的一般方法包括三个步骤（也见 @ref GlossPeriodicConstraints "关于周期性边界条件的词汇表条目"）。

-# 创建一个网格

-# 识别边界不同部分的那些面对，在这些面对上的解应该是对称的，使用 GridTools::collect_periodic_faces()  。

-# 使用 parallel::distributed::Triangulation::add_periodicity() 将周期性信息添加到网格中。

-# 使用 DoFTools::make_periodicity_constraints() 添加周期性约束。

第二和第三步对于使用 parallel::distributed::Triangulation 类的平行网格是必要的，以确保位于域的对面但由周期性面连接的单元是鬼层的一部分，如果它们中的一个存储在本地处理器上。如果三角结构不是 parallel::distributed::Triangulation, 类，这些步骤就没有必要。

第一步包括收集匹配的周期性面，并将它们存储在 <code>std::vector</code> 的 GridTools::PeriodicFacePair. 中，这是通过函数 GridTools::collect_periodic_faces() 来完成的，例如可以这样调用。

@code
GridTools::collect_periodic_faces(dof_handler,
                                  b_id1,
                                  b_id2,
                                  direction,
                                  matched_pairs,
                                  offset = <default value>,
                                  matrix = <default value>,
                                  first_vector_components = <default value>);
@endcode



这个调用在周期性边界上的容器dof_handler的所有面进行循环，其边界指标分别为 @p b_id1  和  @p b_id2, 。(你可以在创建粗略网格后手工分配这些边界指标，见 @ref GlossBoundaryIndicator  "边界指标"。另外，如果你指定了 "着色 "标志，你也可以让命名空间GridGenerator中的许多函数来做这件事；在这种情况下，这些函数会给边界的不同部分分配不同的边界指标，细节通常在这些函数的文档中详细说明）。)

具体来说，如果 $\text{vertices}_{1/2}$ 是两个面 $\text{face}_{1/2}$ 的顶点，那么上面的函数调用将匹配成对的面（和道夫），使得 $\text{vertices}_2$ 和 $matrix\cdot \text{vertices}_1+\text{offset}$ 之间的差异在除方向之外的每个分量中都消失，并将产生的对与相关数据存储在 @p matched_pairs. 中（关于匹配过程的详细信息，见 GridTools::orthogonal_equality() ）。

例如，考虑彩色单位方块 $\Omega=[0,1]^2$ ，其边界指标0在左边，1在右边，2在下面，3在上面的面。参见 GridGenerator::hyper_cube() 的文件，了解这个关于如何分配边界指标的惯例）。然后。

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs);
@endcode

将产生周期性约束，使 $u(0,y)=u(1,y)$ 为所有 $y\in[0,1]$ 。

如果我们转而考虑由 $(0,0)$ ,  $(1,1)$ ,  $(1,2)$ ,  $(0,1)$ 的凸壳给出的平行四边形，我们可以通过指定一个 @p offset: 来实现约束 $u(0,y)=u(1,y+1)$  。

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(0.,1.));
@endcode

或

@code
GridTools::collect_periodic_faces(dof_handler,
                                  /*b_id1*/ 0,
                                  /*b_id2*/ 1,
                                  /*arbitrary direction*/ 0,
                                  matched_pairs,
                                  Tensor<1, 2>(1.,1.));
@endcode

在这里，边界指标0和1的分配也是源于 GridGenerator::parallelogram() 的文件。

由此产生的 @p matched_pairs 可以在 DoFTools::make_periodicity_constraints 中使用，用于为AffineConstraints对象填充周期性约束。

@code
DoFTools::make_periodicity_constraints(matched_pairs, constraints);
@endcode



除了这个高水平的接口外，还有结合这两个步骤的 DoFTools::make_periodicity_constraints 的变体（见 DofTools::make_periodicity_constraints). 的变体）可用。

如果需要更多的灵活性，也有一个与 DoFTools::make_periodicity_constraints 的低级接口。该低级变体允许直接指定两个应被约束的面。

@code
using namespace DoFTools;
make_periodicity_constraints(face_1,
                             face_2,
                             affine_constraints,
                             component_mask = <default value>;
                             face_orientation = <default value>,
                             face_flip = <default value>,
                             face_rotation = <default value>,
                             matrix = <default value>);
@endcode

在这里，我们需要使用 @p face_orientation,  @p face_flip 和 @p face_orientation. 来指定两个面的方向。 要想了解更详细的描述，请看 DoFTools::make_periodicity_constraints. 的文档。除了自我解释的 @p component_mask 和 @p affine_constraints. 之外，其余的参数与高级接口相同。


<a name="problem"></a>

<a name="Apracticalexample"></a><h1>A practical example</h1>


在下文中，我们将展示如何在一个更复杂的例子中使用上述函数。任务是对斯托克斯流的速度分量实施旋转的周期性约束。

在由 $\Omega=\{{\bf x}\in(0,1)^2:\|{\bf x}\|\in (0.5,1)\}$ 定义的四分之一圆上，我们要解决斯托克斯问题

@f{eqnarray*}


  -\Delta \; \textbf{u} + \nabla p &=& (\exp(-100\|{\bf x}-(.75,0.1)^T\|^2),0)^T, \\


  -\textrm{div}\;  \textbf{u}&=&0,\\
  \textbf{u}|_{\Gamma_1}&=&{\bf 0},


@f}

其中边界 $\Gamma_1$ 被定义为 $\Gamma_1 \dealcoloneq \{x\in \partial\Omega: \|x\|\in\{0.5,1\}\}$  。对于边界的其余部分，我们将使用周期性边界条件，即

@f{align*}
  u_x(0,\nu)&=-u_y(\nu,0)&\nu&\in[0,1]\\
  u_y(0,\nu)&=u_x(\nu,0)&\nu&\in[0,1].


@f}



网格将由 GridGenerator::quarter_hyper_shell(), 生成，该文件还记录了如果其`colorize'参数设置为`true'，它是如何为其各种边界分配边界指标的。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 这个例子程序是对 step-22 的轻微修改，使用Trilinos并行运行，以演示交易.II中周期性边界条件的使用。因此我们不讨论大部分的源代码，只对处理周期性约束的部分进行评论。其余的请看 step-22 和底部的完整源代码。
 * 

 * 
 * 为了实现周期性边界条件，只有两个函数需要修改。
 * 

 * 
 * -  <code>StokesProblem<dim>::setup_dofs()</code>  : 用周期性约束来填充AffineConstraints对象
 * 

 * 
 * -  <code>StokesProblem<dim>::create_mesh()</code>  : 为分布式三角形提供周期性信息。
 * 

 * 
 * 程序的其余部分与 step-22 相同，所以让我们跳过这一部分，只在下面显示这两个函数。完整的程序可以在下面的 "普通程序 "部分找到）。
 * 

 * 
 * @cond  跳过
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * 
 * #include <deal.II/distributed/grid_refinement.h> 
 * 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/lac/trilinos_solver.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_parallel_block_vector.h> 
 * #include <deal.II/lac/block_sparsity_pattern.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/mapping_q.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * namespace Step45 
 * { 
 *   using namespace dealii; 
 * 
 *   template <int dim> 
 *   class StokesProblem 
 *   { 
 *   public: 
 *     StokesProblem(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     void create_mesh(); 
 *     void setup_dofs(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void output_results(const unsigned int refinement_cycle) const; 
 *     void refine_mesh(); 
 * 
 *     const unsigned int degree; 
 * 
 *     MPI_Comm mpi_communicator; 
 * 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 *     FESystem<dim>                             fe; 
 *     DoFHandler<dim>                           dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 *     std::vector<IndexSet>     owned_partitioning; 
 *     std::vector<IndexSet>     relevant_partitioning; 
 * 
 *     TrilinosWrappers::BlockSparseMatrix system_matrix; 
 * 
 *     TrilinosWrappers::BlockSparseMatrix preconditioner_matrix; 
 * 
 *     TrilinosWrappers::MPI::BlockVector solution; 
 *     TrilinosWrappers::MPI::BlockVector system_rhs; 
 * 
 *     ConditionalOStream pcout; 
 * 
 *     MappingQ<dim> mapping; 
 *   }; 
 * 
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     BoundaryValues() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  value) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double BoundaryValues<dim>::value(const Point<dim> & /*p*/, 
 *                                     const unsigned int component) const 
 *   { 
 *     (void)component; 
 *     Assert(component < this->n_components, 
 *            ExcIndexRange(component, 0, this->n_components)); 
 * 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   void BoundaryValues<dim>::vector_value(const Point<dim> &p, 
 *                                          Vector<double> &  values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = BoundaryValues<dim>::value(p, c); 
 *   } 
 * 
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     RightHandSide() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *     virtual void vector_value(const Point<dim> &p, 
 *                               Vector<double> &  value) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double RightHandSide<dim>::value(const Point<dim> & p, 
 *                                    const unsigned int component) const 
 *   { 
 *     const Point<dim> center(0.75, 0.1); 
 *     const double     r = (p - center).norm(); 
 * 
 *     if (component == 0) 
 *       return std::exp(-100. * r * r); 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   void RightHandSide<dim>::vector_value(const Point<dim> &p, 
 *                                         Vector<double> &  values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = RightHandSide<dim>::value(p, c); 
 *   } 
 * 
 *   template <class MatrixType, class PreconditionerType> 
 *   class InverseMatrix : public Subscriptor 
 *   { 
 *   public: 
 *     InverseMatrix(const MatrixType &        m, 
 *                   const PreconditionerType &preconditioner, 
 *                   const IndexSet &          locally_owned, 
 *                   const MPI_Comm &          mpi_communicator); 
 * 
 *     void vmult(TrilinosWrappers::MPI::Vector &      dst, 
 *                const TrilinosWrappers::MPI::Vector &src) const; 
 * 
 *   private: 
 *     const SmartPointer<const MatrixType>         matrix; 
 *     const SmartPointer<const PreconditionerType> preconditioner; 
 * 
 *     const MPI_Comm *                      mpi_communicator; 
 *     mutable TrilinosWrappers::MPI::Vector tmp; 
 *   }; 
 * 
 *   template <class MatrixType, class PreconditionerType> 
 *   InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
 *     const MatrixType &        m, 
 *     const PreconditionerType &preconditioner, 
 *     const IndexSet &          locally_owned, 
 *     const MPI_Comm &          mpi_communicator) 
 *     : matrix(&m) 
 *     , preconditioner(&preconditioner) 
 *     , mpi_communicator(&mpi_communicator) 
 *     , tmp(locally_owned, mpi_communicator) 
 *   {} 
 * 
 *   template <class MatrixType, class PreconditionerType> 
 *   void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
 *     TrilinosWrappers::MPI::Vector &      dst, 
 *     const TrilinosWrappers::MPI::Vector &src) const 
 *   { 
 *     SolverControl              solver_control(src.size(), 1e-6 * src.l2_norm()); 
 *     TrilinosWrappers::SolverCG cg(solver_control, 
 *                                   TrilinosWrappers::SolverCG::AdditionalData()); 
 * 
 *     tmp = 0.; 
 *     cg.solve(*matrix, tmp, src, *preconditioner); 
 *     dst = tmp; 
 *   } 
 * 
 *   template <class PreconditionerType> 
 *   class SchurComplement : public TrilinosWrappers::SparseMatrix 
 *   { 
 *   public: 
 *     SchurComplement(const TrilinosWrappers::BlockSparseMatrix &system_matrix, 
 *                     const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                                         PreconditionerType> &  A_inverse, 
 *                     const IndexSet &                           owned_pres, 
 *                     const MPI_Comm &mpi_communicator); 
 * 
 *     void vmult(TrilinosWrappers::MPI::Vector &      dst, 
 *                const TrilinosWrappers::MPI::Vector &src) const; 
 * 
 *   private: 
 *     const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> system_matrix; 
 *     const SmartPointer< 
 *       const InverseMatrix<TrilinosWrappers::SparseMatrix, PreconditionerType>> 
 *                                           A_inverse; 
 *     mutable TrilinosWrappers::MPI::Vector tmp1, tmp2; 
 *   }; 
 * 
 *   template <class PreconditionerType> 
 *   SchurComplement<PreconditionerType>::SchurComplement( 
 *     const TrilinosWrappers::BlockSparseMatrix &system_matrix, 
 *     const InverseMatrix<TrilinosWrappers::SparseMatrix, PreconditionerType> 
 *       &             A_inverse, 
 *     const IndexSet &owned_vel, 
 *     const MPI_Comm &mpi_communicator) 
 *     : system_matrix(&system_matrix) 
 *     , A_inverse(&A_inverse) 
 *     , tmp1(owned_vel, mpi_communicator) 
 *     , tmp2(tmp1) 
 *   {} 
 * 
 *   template <class PreconditionerType> 
 *   void SchurComplement<PreconditionerType>::vmult( 
 *     TrilinosWrappers::MPI::Vector &      dst, 
 *     const TrilinosWrappers::MPI::Vector &src) const 
 *   { 
 *     system_matrix->block(0, 1).vmult(tmp1, src); 
 *     A_inverse->vmult(tmp2, tmp1); 
 *     system_matrix->block(1, 0).vmult(dst, tmp2); 
 *   } 
 * 
 *   template <int dim> 
 *   StokesProblem<dim>::StokesProblem(const unsigned int degree) 
 *     : degree(degree) 
 *     , mpi_communicator(MPI_COMM_WORLD) 
 *     , triangulation(mpi_communicator) 
 *     , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1) 
 *     , dof_handler(triangulation) 
 *     , pcout(std::cout, Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
 *     , mapping(degree + 1) 
 *   {} 
 * @endcode
 * 
 * @endcond  
 * 
 * <a name="Settingupperiodicityconstraintsondistributedtriangulations"></a> 
 * <h3>Setting up periodicity constraints on distributed triangulations</h3>
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::create_mesh() 
 *   { 
 *     Point<dim>   center; 
 *     const double inner_radius = .5; 
 *     const double outer_radius = 1.; 
 * 
 *     GridGenerator::quarter_hyper_shell( 
 *       triangulation, center, inner_radius, outer_radius, 0, true); 
 * 
 * @endcode
 * 
 * 在我们可以规定周期性约束之前，我们需要确保位于域的对面但由周期性面连接的单元是幽灵层的一部分，如果其中一个单元存储在本地处理器上。在这一点上，我们需要考虑我们要如何规定周期性。左边边界上的面的顶点 $\text{vertices}_2$ 应该与下面边界上的面的顶点 $\text{vertices}_1$ 相匹配，由 $\text{vertices}_2=R\cdot \text{vertices}_1+b$ 给出，其中旋转矩阵 $R$ 和偏移量 $b$ 由
 * @f{align*}
 * R=\begin{pmatrix}
 * 0&1\\-1&0
 * \end{pmatrix},
 * \quad
 * b=\begin{pmatrix}0&0\end{pmatrix}.
 * @f}
 * 给出。 我们将所得信息保存到这里的数据结构是基于三角结构的。
 * 

 * 
 * 
 * @code
 *     std::vector<GridTools::PeriodicFacePair< 
 *       typename parallel::distributed::Triangulation<dim>::cell_iterator>> 
 *       periodicity_vector; 
 * 
 *     FullMatrix<double> rotation_matrix(dim); 
 *     rotation_matrix[0][1] = 1.; 
 *     rotation_matrix[1][0] = -1.; 
 * 
 *     GridTools::collect_periodic_faces(triangulation, 
 *                                       2, 
 *                                       3, 
 *                                       1, 
 *                                       periodicity_vector, 
 *                                       Tensor<1, dim>(), 
 *                                       rotation_matrix); 
 * 
 * @endcode
 * 
 * 现在，只要调用 parallel::distributed::Triangulation::add_periodicity. 就可以告诉三角函数所需的周期性，特别容易。
 * 
 * @code
 *     triangulation.add_periodicity(periodicity_vector); 
 * 
 *     triangulation.refine_global(4 - dim); 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::setup_dofs() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::vector<unsigned int> block_component(dim + 1, 0); 
 *     block_component[dim] = 1; 
 *     DoFRenumbering::component_wise(dof_handler, block_component); 
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 *     const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1]; 
 * 
 *     { 
 *       owned_partitioning.clear(); 
 *       IndexSet locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 *       owned_partitioning.push_back(locally_owned_dofs.get_view(0, n_u)); 
 *       owned_partitioning.push_back(locally_owned_dofs.get_view(n_u, n_u + n_p)); 
 * 
 *       relevant_partitioning.clear(); 
 *       IndexSet locally_relevant_dofs; 
 *       DoFTools::extract_locally_relevant_dofs(dof_handler, 
 *                                               locally_relevant_dofs); 
 *       relevant_partitioning.push_back(locally_relevant_dofs.get_view(0, n_u)); 
 *       relevant_partitioning.push_back( 
 *         locally_relevant_dofs.get_view(n_u, n_u + n_p)); 
 * 
 *       constraints.clear(); 
 *       constraints.reinit(locally_relevant_dofs); 
 * 
 *       FEValuesExtractors::Vector velocities(0); 
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *       VectorTools::interpolate_boundary_values(mapping, 
 *                                                dof_handler, 
 *                                                0, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 *       VectorTools::interpolate_boundary_values(mapping, 
 *                                                dof_handler, 
 *                                                1, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 * 
 * @endcode
 * 
 * 在我们为网格提供了周期性约束的必要信息后，我们现在可以实际创建它们。对于描述匹配，我们使用与之前相同的方法，也就是说，左边边界上的一个面的 $\text{vertices}_2$ 应该与下面边界上的一个面的顶点 $\text{vertices}_1$ 匹配，由 $\text{vertices}_2=R\cdot \text{vertices}_1+b$  ]，其中旋转矩阵 $R$ 和偏移量 $b$ 由
 * @f{align*}
 * R=\begin{pmatrix}
 * 0&1\\-1&0
 * \end{pmatrix},
 * \quad
 * b=\begin{pmatrix}0&0\end{pmatrix}.
 * @f}
 * 给出。 这两个对象不仅描述了应该如何匹配面，而且还描述了解决方案应该从 $\text{face}_2$ 转换到 $\text{face}_1$ 的意义。
 * 

 * 
 * 
 * @code
 *       FullMatrix<double> rotation_matrix(dim); 
 *       rotation_matrix[0][1] = 1.; 
 *       rotation_matrix[1][0] = -1.; 
 * 
 *       Tensor<1, dim> offset; 
 * 
 * @endcode
 * 
 * 为了设置约束，我们首先将周期性信息存储在一个类型为 <code>std::vector@<GridTools::PeriodicFacePair<typename 的辅助对象中。
 * DoFHandler@<dim@>::%cell_iterator@>  </code>。周期性边界的边界指标为2（x=0）和3（y=0）。所有其他的参数我们之前已经设置好了。在这种情况下，方向并不重要。由于 $\text{vertices}_2=R\cdot \text{vertices}_1+b$ 这正是我们想要的。
 * 

 * 
 * 
 * @code
 *       std::vector< 
 *         GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>> 
 *         periodicity_vector; 
 * 
 *       const unsigned int direction = 1; 
 * 
 *       GridTools::collect_periodic_faces(dof_handler, 
 *                                         2, 
 *                                         3, 
 *                                         direction, 
 *                                         periodicity_vector, 
 *                                         offset, 
 *                                         rotation_matrix); 
 * 
 * @endcode
 * 
 * 接下来，我们需要提供关于解决方案中哪些矢量值分量应该被旋转的信息。由于我们在这里选择只约束速度，并且从解决方案矢量的第一个分量开始，我们只需插入一个0。
 * 

 * 
 * 
 * @code
 *       std::vector<unsigned int> first_vector_components; 
 *       first_vector_components.push_back(0); 
 * 
 * @endcode
 * 
 * 在设置了周期性_vector中的所有信息之后，我们要做的就是告诉make_periodicity_constraints来创建所需的约束。
 * 

 * 
 * 
 * @code
 *       DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector, 
 *                                                        constraints, 
 *                                                        fe.component_mask( 
 *                                                          velocities), 
 *                                                        first_vector_components); 
 * 
 *       VectorTools::interpolate_boundary_values(mapping, 
 *                                                dof_handler, 
 *                                                0, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 *       VectorTools::interpolate_boundary_values(mapping, 
 *                                                dof_handler, 
 *                                                1, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 *     } 
 * 
 *     constraints.close(); 
 * 
 *     { 
 *       TrilinosWrappers::BlockSparsityPattern bsp(owned_partitioning, 
 *                                                  owned_partitioning, 
 *                                                  relevant_partitioning, 
 *                                                  mpi_communicator); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (!((c == dim) && (d == dim))) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler, 
 *                                       coupling, 
 *                                       bsp, 
 *                                       constraints, 
 *                                       false, 
 *                                       Utilities::MPI::this_mpi_process( 
 *                                         mpi_communicator)); 
 * 
 *       bsp.compress(); 
 * 
 *       system_matrix.reinit(bsp); 
 *     } 
 * 
 *     { 
 *       TrilinosWrappers::BlockSparsityPattern preconditioner_bsp( 
 *         owned_partitioning, 
 *         owned_partitioning, 
 *         relevant_partitioning, 
 *         mpi_communicator); 
 * 
 *       Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if ((c == dim) && (d == dim)) 
 *             preconditioner_coupling[c][d] = DoFTools::always; 
 *           else 
 *             preconditioner_coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler, 
 *                                       preconditioner_coupling, 
 *                                       preconditioner_bsp, 
 *                                       constraints, 
 *                                       false, 
 *                                       Utilities::MPI::this_mpi_process( 
 *                                         mpi_communicator)); 
 * 
 *       preconditioner_bsp.compress(); 
 * 
 *       preconditioner_matrix.reinit(preconditioner_bsp); 
 *     } 
 * 
 *     system_rhs.reinit(owned_partitioning, mpi_communicator); 
 *     solution.reinit(owned_partitioning, 
 *                     relevant_partitioning, 
 *                     mpi_communicator); 
 *   } 
 * 
 * @endcode
 * 
 * 然后程序的其余部分又与  step-22  相同。我们现在省略它，但和以前一样，你可以在下面的 "普通程序 "部分找到这些部分。
 * 

 * 
 * @cond  SKIP
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::assemble_system() 
 *   { 
 *     system_matrix         = 0.; 
 *     system_rhs            = 0.; 
 *     preconditioner_matrix = 0.; 
 * 
 *     QGauss<dim> quadrature_formula(degree + 2); 
 * 
 *     FEValues<dim> fe_values(mapping, 
 *                             fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values | update_gradients); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> local_preconditioner_matrix(dofs_per_cell, 
 *                                                    dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const RightHandSide<dim>    right_hand_side; 
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell); 
 *     std::vector<double>                  div_phi_u(dofs_per_cell); 
 *     std::vector<double>                  phi_p(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           fe_values.reinit(cell); 
 *           local_matrix                = 0; 
 *           local_preconditioner_matrix = 0; 
 *           local_rhs                   = 0; 
 * 
 *           right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
 *                                             rhs_values); 
 * 
 *           for (unsigned int q = 0; q < n_q_points; ++q) 
 *             { 
 *               for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *                 { 
 *                   symgrad_phi_u[k] = 
 *                     fe_values[velocities].symmetric_gradient(k, q); 
 *                   div_phi_u[k] = fe_values[velocities].divergence(k, q); 
 *                   phi_p[k]     = fe_values[pressure].value(k, q); 
 *                 } 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   for (unsigned int j = 0; j <= i; ++j) 
 *                     { 
 *                       local_matrix(i, j) += 
 *                         (symgrad_phi_u[i] * symgrad_phi_u[j] // diffusion 
 *                          - div_phi_u[i] * phi_p[j]           // pressure force 
 *                          - phi_p[i] * div_phi_u[j])          // divergence 
 *                         * fe_values.JxW(q); 
 * 
 *                       local_preconditioner_matrix(i, j) += 
 *                         (phi_p[i] * phi_p[j]) * fe_values.JxW(q); 
 *                     } 
 * 
 *                   const unsigned int component_i = 
 *                     fe.system_to_component_index(i).first; 
 *                   local_rhs(i) += fe_values.shape_value(i, q)  // 
 *                                   * rhs_values[q](component_i) // 
 *                                   * fe_values.JxW(q); 
 *                 } 
 *             } 
 * 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
 *               { 
 *                 local_matrix(i, j) = local_matrix(j, i); 
 *                 local_preconditioner_matrix(i, j) = 
 *                   local_preconditioner_matrix(j, i); 
 *               } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           constraints.distribute_local_to_global(local_matrix, 
 *                                                  local_rhs, 
 *                                                  local_dof_indices, 
 *                                                  system_matrix, 
 *                                                  system_rhs); 
 *           constraints.distribute_local_to_global(local_preconditioner_matrix, 
 *                                                  local_dof_indices, 
 *                                                  preconditioner_matrix); 
 *         } 
 * 
 *     system_matrix.compress(VectorOperation::add); 
 *     system_rhs.compress(VectorOperation::add); 
 * 
 *     pcout << "   Computing preconditioner..." << std::endl << std::flush; 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::solve() 
 *   { 
 *     TrilinosWrappers::PreconditionJacobi A_preconditioner; 
 *     A_preconditioner.initialize(system_matrix.block(0, 0)); 
 * 
 *     const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                         TrilinosWrappers::PreconditionJacobi> 
 *       A_inverse(system_matrix.block(0, 0), 
 *                 A_preconditioner, 
 *                 owned_partitioning[0], 
 *                 mpi_communicator); 
 * 
 *     TrilinosWrappers::MPI::BlockVector tmp(owned_partitioning, 
 *                                            mpi_communicator); 
 * 
 *  
 *       TrilinosWrappers::MPI::Vector schur_rhs(owned_partitioning[1], 
 *                                               mpi_communicator); 
 *       A_inverse.vmult(tmp.block(0), system_rhs.block(0)); 
 *       system_matrix.block(1, 0).vmult(schur_rhs, tmp.block(0)); 
 *       schur_rhs -= system_rhs.block(1); 
 * 
 *       SchurComplement<TrilinosWrappers::PreconditionJacobi> schur_complement( 
 *         system_matrix, A_inverse, owned_partitioning[0], mpi_communicator); 
 * 
 *       SolverControl solver_control(solution.block(1).size(), 
 *                                    1e-6 * schur_rhs.l2_norm()); 
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 
 * 
 *       TrilinosWrappers::PreconditionAMG preconditioner; 
 *       preconditioner.initialize(preconditioner_matrix.block(1, 1)); 
 * 
 *       InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                     TrilinosWrappers::PreconditionAMG> 
 *         m_inverse(preconditioner_matrix.block(1, 1), 
 *                   preconditioner, 
 *                   owned_partitioning[1], 
 *                   mpi_communicator); 
 * 
 *       cg.solve(schur_complement, tmp.block(1), schur_rhs, preconditioner); 
 * 
 *       constraints.distribute(tmp); 
 *       solution.block(1) = tmp.block(1); 
 *     } 
 * 
 *     { 
 *       system_matrix.block(0, 1).vmult(tmp.block(0), tmp.block(1)); 
 *       tmp.block(0) *= -1; 
 *       tmp.block(0) += system_rhs.block(0); 
 * 
 *       A_inverse.vmult(tmp.block(0), tmp.block(0)); 
 * 
 *       constraints.distribute(tmp); 
 *       solution.block(0) = tmp.block(0); 
 *     } 
 *   } 
 * 
 *   template <int dim> 
 *   void 
 *   StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const 
 *   { 
 *     std::vector<std::string> solution_names(dim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     Vector<float> subdomain(triangulation.n_active_cells()); 
 *     for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *       subdomain(i) = triangulation.locally_owned_subdomain(); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 *     data_out.build_patches(mapping, degree + 1); 
 * 
 *     data_out.write_vtu_with_pvtu_record( 
 *       "./", "solution", refinement_cycle, MPI_COMM_WORLD, 2); 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::refine_mesh() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     FEValuesExtractors::Scalar pressure(dim); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       estimated_error_per_cell, 
 *       fe.component_mask(pressure)); 
 * 
 *     parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number( 
 *       triangulation, estimated_error_per_cell, 0.3, 0.0); 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 *   template <int dim> 
 *   void StokesProblem<dim>::run() 
 *   { 
 *     create_mesh(); 
 * 
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 9; 
 *          ++refinement_cycle) 
 *       { 
 *         pcout << "Refinement cycle " << refinement_cycle << std::endl; 
 * 
 *         if (refinement_cycle > 0) 
 *           refine_mesh(); 
 * 
 *         setup_dofs(); 
 * 
 *         pcout << "   Assembling..." << std::endl << std::flush; 
 *         assemble_system(); 
 * 
 *         pcout << "   Solving..." << std::flush; 
 *         solve(); 
 * 
 *         output_results(refinement_cycle); 
 * 
 *         pcout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step45 
 * 
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step45; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
 *       StokesProblem<2>                 flow_problem(1); 
 *       flow_problem.run(); 
 *     } 
 *   catch (std::exception &exc) 
 *     { 
 *       std::cerr << std::endl 
 *                 << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       std::cerr << "Exception on processing: " << std::endl 
 *                 << exc.what() << std::endl 
 *                 << "Aborting!" << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 * 
 *       return 1; 
 *     } 
 *   catch (...) 
 *     { 
 *       std::cerr << std::endl 
 *                 << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       std::cerr << "Unknown exception!" << std::endl 
 *                 << "Aborting!" << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       return 1; 
 *     } 
 * 
 *   return 0; 
 * } 
 * @endcode
 * 
 * @endcond  
 * 

 * 
 * 

 * 
examples/step-45/doc/results.dox



<a name="Results"></a><h1>Results</h1>


创建的输出结果并不十分令人惊讶。我们只是看到，相对于左、下边界而言，解决方案是周期性的。

 <img src="https://www.dealii.org/images/steps/developer/step-45.periodic.png" alt=""> 

如果没有周期性约束，我们最终会得到以下解决方案。

 <img src="https://www.dealii.org/images/steps/developer/step-45.non_periodic.png" alt=""> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-45.cc"
*/
