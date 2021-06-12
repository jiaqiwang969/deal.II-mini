/**
@page step_12 The step-12 tutorial program
This tutorial depends on step-7.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Overview">Overview</a>
        <li><a href="#Theequation">The equation</a>
        <li><a href="#Thetestproblem">The test problem</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#TheScratchDataandCopyDataclasses">The ScratchData and CopyData classes</a>
        <li><a href="#TheAdvectionProblemclass">The AdvectionProblem class</a>
      <ul>
        <li><a href="#Theassemble_systemfunction">The assemble_system function</a>
      </ul>
        <li><a href="#Alltherest">All the rest</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Whyusediscontinuouselements">Why use discontinuous elements</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-12/doc/intro.dox

 <br> 

<i> Note: A variant called step-12b of this tutorial exists, using
MeshWorker and LocalIntegrators instead of assembling matrices using
FEInterfaceValues as is done in this tutorial.
</i>

<a name="Intro"></a>

<a name="AnexampleofanadvectionproblemusingtheDiscountinuousGalerkinmethod"></a><h1>An example of an advection problem using the Discountinuous Galerkin method</h1>


<a name="Overview"></a><h3>Overview</h3>


本例专门介绍了 <em> 非连续Galerkin方法 </em> ，简称为DG方法。它包括以下内容。<ol>  <li>  用DG方法对线性平流方程进行离散化。     <li>  使用FEInterfaceValues组装跳跃项和单元间界面上的其他表达式。     <li>  使用 MeshWorker::mesh_loop().   </ol>  组合系统矩阵。

这个项目特别关注的是DG方法的循环。这些被证明是特别复杂的，主要是因为对于面的条件，我们必须分别区分边界、规则的内部面和有悬挂节点的内部面的情况。 MeshWorker::mesh_loop() 处理了单元和面迭代的复杂性，并允许为不同的单元和面项指定 "工作者"。面条款本身的整合，包括对自适应细化面的整合，是通过FEInterfaceValues类完成的。

<a name="Theequation"></a><h3>The equation</h3>


本例中解决的模型问题是线性平流方程

@f[
  \nabla\cdot \left({\mathbf \beta} u\right)=0 \qquad\mbox{in }\Omega,


@f]

受制于边界条件

@f[
u=g\quad\mbox{on }\Gamma_-,


@f]

在域的边界 $\Gamma=\partial\Omega$ 的流入部分 $\Gamma_-$ 。  这里， ${\mathbf \beta}={\mathbf \beta}({\bf x})$ 表示一个矢量场， $u$ 是（标量）解函数， $g$ 是边界值函数。

@f[
\Gamma_- \dealcoloneq \{{\bf x}\in\Gamma, {\mathbf \beta}({\bf x})\cdot{\bf n}({\bf x})<0\}


@f]

表示域边界的流入部分， ${\bf n}$ 表示边界的单位外向法线 $\Gamma$  。这个方程是本教程第9步中已经考虑过的平流方程的保守版本。


在每个单元格 $T$ 上，我们从左边乘以一个测试函数 $v_h$ ，并通过部分整合得到。

@f[
  \left( v_h, \nabla \cdot (\beta u_h) \right)_T
= -(\nabla v_h, \beta u_h) + \int_\Gamma v_h u_h \beta \cdot n


@f]

当对所有单元 $T$ 求和时，边界积分是在所有内部和外部面进行的，因此，有三种情况。<ol>  <li>  流入的外部边界（我们用给定的 $g$ 代替 $u_h$ ）。     $\int_{\Gamma_-} v_h g \beta \cdot n$   <li>  流出的外部边界。     $\int_{\Gamma_+} v_h u_h \beta \cdot n$   <li> 内面（两边的积分变成了跳跃，我们使用上风速度）。     $\int_F [v_h] u_h^{\text{upwind}} \beta \cdot n$   </ol> 。

这里，跳跃被定义为 $[v] = v^+ - v^-$ ，其中上标指的是面的左（'+'）和右（'-'）值。如果 $\beta \cdot n>0$ ，上风值 $u^{\text{upwind}}$ 被定义为 $u^+$ ，否则为 $u^-$ 。

因此，依赖网格的弱形式为：。

@f[
\sum_{T\in \mathbb T_h} -\bigl(\nabla \phi_i,{\mathbf \beta}\cdot \phi_j \bigr)_T +
\sum_{F\in\mathbb F_h^i} \bigl< [\phi_i], \phi_j^{upwind} \beta\cdot \mathbf n\bigr>_{F} +
\bigl<\phi_i, \phi_j \beta\cdot \mathbf n\bigr>_{\Gamma_+}
= -\bigl<\phi_i, g \beta\cdot\mathbf n\bigr>_{\Gamma_-}.


@f]

这里， $\mathbb T_h$ 是三角形的所有活动单元的集合， $\mathbb F_h^i$ 是所有活动内部面的集合。这种公式被称为上风非连续Galerkin方法。

为了实现这种双线性形式，我们需要用通常的方法计算单元项（第一个和）来实现单元上的积分，用FEInterfaceValues计算界面项（第二个和），以及边界项（另外两个项）。所有这些的求和是通过 MeshWorker::mesh_loop(). 完成的。




<a name="Thetestproblem"></a><h3>The test problem</h3>


我们在 $\Omega=[0,1]^2$ 上求解平流方程， ${\mathbf \beta}=\frac{1}{|x|}(-x_2, x_1)$ 代表一个圆形的逆时针流场， $g=1$ 代表 ${\bf x}\in\Gamma_-^1 := [0,0.5]\times\{0\}$ ， $g=0$ 代表 ${\bf x}\in
\Gamma_-\setminus \Gamma_-^1$  。

我们通过估计每个单元的梯度规范，自适应地细化网格，在一连串的网格上求解。在每个网格上求解后，我们以vtk格式输出解决方案，并计算解决方案的 $L^\infty$ 准则。由于精确解是0或1，我们可以用它来衡量数值解的过冲程度。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 前面几个文件已经在前面的例子中讲过了，因此不再做进一步的评论。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/fe/mapping_q1.h> 
 * 
 * @endcode
 * 
 * 这里定义了不连续的有限元。它们的使用方式与所有其他有限元相同，不过--正如你在以前的教程程序中所看到的--用户与有限元类的交互根本不多：它们被传递给 <code>DoFHandler</code> 和 <code>FEValues</code> 对象，仅此而已。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_dgq.h> 
 * 
 * @endcode
 * 
 * FEInterfaceValues需要这个头来计算界面上的积分。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_interface_values.h> 
 * 
 * @endcode
 * 
 * 我们将使用最简单的求解器，称为Richardson迭代，它代表了一个简单的缺陷修正。这与一个块状SSOR预处理器（定义在precondition_block.h中）相结合，该预处理器使用DG离散产生的系统矩阵的特殊块状结构。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/solver_richardson.h> 
 * #include <deal.II/lac/precondition_block.h> 
 * 
 * @endcode
 * 
 * 我们将使用梯度作为细化指标。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/derivative_approximation.h> 
 * 
 * @endcode
 * 
 * 最后，新的包含文件用于使用MeshWorker框架中的Mesh_loop。
 * 

 * 
 * 
 * @code
 * #include <deal.II/meshworker/mesh_loop.h> 
 * 
 * @endcode
 * 
 * 像所有的程序一样，我们在完成这一部分时，要包括所需的C++头文件，并声明我们要使用dealii命名空间中的对象，不含前缀。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * 
 * namespace Step12 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 首先，我们定义一个描述不均匀边界数据的类。由于只使用它的值，我们实现value_list()，但不定义Function的所有其他函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     BoundaryValues() = default; 
 *     virtual void value_list(const std::vector<Point<dim>> &points, 
 *                             std::vector<double> &          values, 
 *                             const unsigned int component = 0) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 考虑到流动方向，单位方块 $[0,1]^2$ 的流入边界是右边界和下边界。我们在x轴上规定了不连续的边界值1和0，在右边界上规定了值0。该函数在流出边界上的值将不会在DG方案中使用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BoundaryValues<dim>::value_list(const std::vector<Point<dim>> &points, 
 *                                        std::vector<double> &          values, 
 *                                        const unsigned int component) const 
 *   {  
 *     (void)component; 
 *     AssertIndexRange(component, 1); 
 *     Assert(values.size() == points.size(), 
 *            ExcDimensionMismatch(values.size(), points.size())); 
 * 
 *     for (unsigned int i = 0; i < values.size(); ++i) 
 *       { 
 *         if (points[i](0) < 0.5)  
 *           values[i] = 1.; 
 *         else 
 *           values[i] = 0.; 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 最后，一个计算并返回风场的函数  $\beta=\beta(\mathbf x)$  。正如在介绍中所解释的，在2D中我们将使用一个围绕原点的旋转场。在3D中，我们只需不设置 $z$ 分量（即为零），而这个函数在目前的实现中不能用于1D。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   Tensor<1, dim> beta(const Point<dim> &p) 
 *   { 
 *     Assert(dim >= 2, ExcNotImplemented()); 
 * 
 *     Tensor<1, dim> wind_field; 
 *     wind_field[0] = -p[1]; 
 *     wind_field[1] = p[0]; 
 * 
 *     if (wind_field.norm() > 1e-10) 
 *       wind_field /= wind_field.norm(); 
 * 
 *     return wind_field; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TheScratchDataandCopyDataclasses"></a> 
 * <h3>The ScratchData and CopyData classes</h3>
 * 

 * 
 * 以下对象是我们在调用 MeshWorker::mesh_loop(). 时使用的抓取和复制对象 新对象是FEInterfaceValues对象，它的工作原理类似于FEValues或FEFacesValues，只是它作用于两个单元格之间的接口，并允许我们以我们的弱形式组装接口条款。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct ScratchData 
 *   { 
 *     ScratchData(const Mapping<dim> &       mapping, 
 *                 const FiniteElement<dim> & fe, 
 *                 const Quadrature<dim> &    quadrature, 
 *                 const Quadrature<dim - 1> &quadrature_face, 
 *                 const UpdateFlags          update_flags = update_values | 
 *                                                  update_gradients | 
 *                                                  update_quadrature_points | 
 *                                                  update_JxW_values, 
 *                 const UpdateFlags interface_update_flags = 
 *                   update_values | update_gradients | update_quadrature_points | 
 *                   update_JxW_values | update_normal_vectors) 
 *       : fe_values(mapping, fe, quadrature, update_flags) 
 *       , fe_interface_values(mapping, 
 *                             fe, 
 *                             quadrature_face, 
 *                             interface_update_flags) 
 *     {} 
 * 
 *     ScratchData(const ScratchData<dim> &scratch_data) 
 *       : fe_values(scratch_data.fe_values.get_mapping(), 
 *                   scratch_data.fe_values.get_fe(), 
 *                   scratch_data.fe_values.get_quadrature(), 
 *                   scratch_data.fe_values.get_update_flags()) 
 *       , fe_interface_values(scratch_data.fe_interface_values.get_mapping(), 
 *                             scratch_data.fe_interface_values.get_fe(), 
 *                             scratch_data.fe_interface_values.get_quadrature(), 
 *                             scratch_data.fe_interface_values.get_update_flags()) 
 *     {} 
 * 
 *     FEValues<dim>          fe_values; 
 *     FEInterfaceValues<dim> fe_interface_values; 
 *   }; 
 * 
 *   struct CopyDataFace 
 *   { 
 *     FullMatrix<double>                   cell_matrix; 
 *     std::vector<types::global_dof_index> joint_dof_indices; 
 *   }; 
 * 
 *   struct CopyData 
 *   { 
 *     FullMatrix<double>                   cell_matrix; 
 *     Vector<double>                       cell_rhs; 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 *     std::vector<CopyDataFace>            face_data; 
 * 
 *     template <class Iterator> 
 *     void reinit(const Iterator &cell, unsigned int dofs_per_cell) 
 *     { 
 *       cell_matrix.reinit(dofs_per_cell, dofs_per_cell); 
 *       cell_rhs.reinit(dofs_per_cell); 
 * 
 *       local_dof_indices.resize(dofs_per_cell); 
 *       cell->get_dof_indices(local_dof_indices); 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="TheAdvectionProblemclass"></a> 
 * <h3>The AdvectionProblem class</h3>
 * 

 * 
 * 在这个准备工作之后，我们继续进行这个程序的主类，称为AdvectionProblem。
 * 

 * 
 * 这对你来说应该是非常熟悉的。有趣的细节只有在实现集合函数的时候才会出现。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class AdvectionProblem 
 *   { 
 *   public: 
 *     AdvectionProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     Triangulation<dim>   triangulation; 
 *     const MappingQ1<dim> mapping; 
 * 
 * @endcode
 * 
 * 此外，我们要使用DG元素。
 * 

 * 
 * 
 * @code
 *     const FE_DGQ<dim> fe; 
 *     DoFHandler<dim>   dof_handler; 
 * 
 *     const QGauss<dim>     quadrature; 
 *     const QGauss<dim - 1> quadrature_face; 
 * 
 * @endcode
 * 
 * 接下来的四个成员代表要解决的线性系统。  <code>system_matrix</code> and <code>right_hand_side</code> 是由 <code>assemble_system()</code>, the <code>solution</code> 产生的，在 <code>solve()</code>. The <code>sparsity_pattern</code> 中计算，用于确定 <code>system_matrix</code> 中非零元素的位置。
 * 

 * 
 * 
 * @code
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> right_hand_side; 
 *   }; 
 * 
 * @endcode
 * 
 * 我们从构造函数开始。 <code>fe</code> 的构造器调用中的1是多项式的度数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   AdvectionProblem<dim>::AdvectionProblem() 
 *     : mapping() 
 *     , fe(1) 
 *     , dof_handler(triangulation) 
 *     , quadrature(fe.tensor_degree() + 1) 
 *     , quadrature_face(fe.tensor_degree() + 1) 
 *   {} 
 * 
 *   template <int dim> 
 *   void AdvectionProblem<dim>::setup_system() 
 *   { 
 * 
 * @endcode
 * 
 * 在设置通常的有限元数据结构的函数中，我们首先需要分配DoF。
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 我们从生成稀疏模式开始。为此，我们首先用系统中出现的耦合物填充一个动态稀疏模式（DynamicSparsityPattern）类型的中间对象。在建立模式之后，这个对象被复制到 <code>sparsity_pattern</code> 并可以被丢弃。
 * 

 * 
 * 为了建立DG离散的稀疏模式，我们可以调用类似于 DoFTools::make_sparsity_pattern, 的函数，该函数被称为 DoFTools::make_flux_sparsity_pattern:  。
 * 
 * @code
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 * @endcode
 * 
 * 最后，我们设置了线性系统的所有组成部分的结构。
 * 

 * 
 * 
 * @code
 *     system_matrix.reinit(sparsity_pattern); 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     right_hand_side.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Theassemble_systemfunction"></a> 
 * <h4>The assemble_system function</h4>
 * 

 * 
 * 这里我们看到了与手工组装的主要区别。我们不需要在单元格和面上写循环，而是在调用 MeshWorker::mesh_loop() 时包含逻辑，我们只需要指定在每个单元格、每个边界面和每个内部面应该发生什么。这三个任务是由下面的函数里面的lambda函数处理的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::assemble_system() 
 *   { 
 *     using Iterator = typename DoFHandler<dim>::active_cell_iterator; 
 *     const BoundaryValues<dim> boundary_function; 
 * 
 * @endcode
 * 
 * 这是将对每个单元格执行的函数。
 * 

 * 
 * 
 * @code
 *     const auto cell_worker = [&](const Iterator &  cell, 
 *                                  ScratchData<dim> &scratch_data, 
 *                                  CopyData &        copy_data) { 
 *       const unsigned int n_dofs = 
 *         scratch_data.fe_values.get_fe().n_dofs_per_cell(); 
 *       copy_data.reinit(cell, n_dofs); 
 *       scratch_data.fe_values.reinit(cell); 
 * 
 *       const auto &q_points = scratch_data.fe_values.get_quadrature_points(); 
 * 
 *       const FEValues<dim> &      fe_v = scratch_data.fe_values; 
 *       const std::vector<double> &JxW  = fe_v.get_JxW_values(); 
 * 
 * @endcode
 * 
 * 我们解决的是一个同质方程，因此在单元项中没有显示出右手。 剩下的就是整合矩阵条目。
 * 

 * 
 * 
 * @code
 *       for (unsigned int point = 0; point < fe_v.n_quadrature_points; ++point) 
 *         { 
 *           auto beta_q = beta(q_points[point]); 
 *           for (unsigned int i = 0; i < n_dofs; ++i) 
 *             for (unsigned int j = 0; j < n_dofs; ++j) 
 *               { 
 *                 copy_data.cell_matrix(i, j) += 
 *                   -beta_q                      // -\beta 
 *                   * fe_v.shape_grad(i, point)  // \nabla \phi_i 
 *                   * fe_v.shape_value(j, point) // \phi_j 
 *                   * JxW[point];                // dx 
 *               } 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 这是为边界面调用的函数，包括使用FEFaceValues的正常积分。新的逻辑是决定该术语是进入系统矩阵（流出）还是进入右手边（流入）。
 * 

 * 
 * 
 * @code
 *     const auto boundary_worker = [&](const Iterator &    cell, 
 *                                      const unsigned int &face_no, 
 *                                      ScratchData<dim> &  scratch_data, 
 *                                      CopyData &          copy_data) { 
 *       scratch_data.fe_interface_values.reinit(cell, face_no); 
 *       const FEFaceValuesBase<dim> &fe_face = 
 *         scratch_data.fe_interface_values.get_fe_face_values(0); 
 * 
 *       const auto &q_points = fe_face.get_quadrature_points(); 
 * 
 *       const unsigned int n_facet_dofs = fe_face.get_fe().n_dofs_per_cell(); 
 *       const std::vector<double> &        JxW     = fe_face.get_JxW_values(); 
 *       const std::vector<Tensor<1, dim>> &normals = fe_face.get_normal_vectors(); 
 * 
 *       std::vector<double> g(q_points.size()); 
 *       boundary_function.value_list(q_points, g); 
 * 
 *       for (unsigned int point = 0; point < q_points.size(); ++point) 
 *         { 
 *           const double beta_dot_n = beta(q_points[point]) * normals[point]; 
 * 
 *           if (beta_dot_n > 0) 
 *             { 
 *               for (unsigned int i = 0; i < n_facet_dofs; ++i) 
 *                 for (unsigned int j = 0; j < n_facet_dofs; ++j) 
 *                   copy_data.cell_matrix(i, j) += 
 *                     fe_face.shape_value(i, point)   // \phi_i 
 *                     * fe_face.shape_value(j, point) // \phi_j 
 *                     * beta_dot_n                    // \beta . n 
 *                     * JxW[point];                   // dx 
 *             } 
 *           else 
 *             for (unsigned int i = 0; i < n_facet_dofs; ++i) 
 *               copy_data.cell_rhs(i) += -fe_face.shape_value(i, point) // \phi_i 
 *                                        * g[point]                     // g 
 *                                        * beta_dot_n  // \beta . n 
 *                                        * JxW[point]; // dx 
 *         } 
 *     }; 
 * 
 * @endcode
 * 
 * 这是在内部面调用的函数。参数指定了单元格、面和子面的指数（用于自适应细化）。我们只是将它们传递给FEInterfaceValues的reinit()函数。
 * 

 * 
 * 
 * @code
 *     const auto face_worker = [&](const Iterator &    cell, 
 *                                  const unsigned int &f, 
 *                                  const unsigned int &sf, 
 *                                  const Iterator &    ncell, 
 *                                  const unsigned int &nf, 
 *                                  const unsigned int &nsf, 
 *                                  ScratchData<dim> &  scratch_data, 
 *                                  CopyData &          copy_data) { 
 *       FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values; 
 *       fe_iv.reinit(cell, f, sf, ncell, nf, nsf); 
 *       const auto &q_points = fe_iv.get_quadrature_points(); 
 * 
 *       copy_data.face_data.emplace_back(); 
 *       CopyDataFace &copy_data_face = copy_data.face_data.back(); 
 * 
 *       const unsigned int n_dofs        = fe_iv.n_current_interface_dofs(); 
 *       copy_data_face.joint_dof_indices = fe_iv.get_interface_dof_indices(); 
 * 
 *       copy_data_face.cell_matrix.reinit(n_dofs, n_dofs); 
 * 
 *       const std::vector<double> &        JxW     = fe_iv.get_JxW_values(); 
 *       const std::vector<Tensor<1, dim>> &normals = fe_iv.get_normal_vectors(); 
 * 
 *       for (unsigned int qpoint = 0; qpoint < q_points.size(); ++qpoint) 
 *         { 
 *           const double beta_dot_n = beta(q_points[qpoint]) * normals[qpoint]; 
 *           for (unsigned int i = 0; i < n_dofs; ++i) 
 *             for (unsigned int j = 0; j < n_dofs; ++j) 
 *               copy_data_face.cell_matrix(i, j) +=  
 *                 fe_iv.jump(i, qpoint) // [\phi_i] 
 *                 * 
 *                 fe_iv.shape_value((beta_dot_n > 0), j, qpoint) // phi_j^{upwind} 
 *                 * beta_dot_n                                   // (\beta . n) 
 *                 * JxW[qpoint];                                 // dx 
 *         }  
 *     }; 
 * 
 * @endcode
 * 
 * 下面的lambda函数将处理从单元格和面组件中复制数据到全局矩阵和右侧的问题。
 * 

 * 
 * 虽然我们不需要AffineConstraints对象，因为在DG离散中没有悬空节点约束，但我们在这里使用一个空对象，因为这允许我们使用其`copy_local_to_global`功能。
 * 

 * 
 * 
 * @code
 *     const AffineConstraints<double> constraints; 
 * 
 *     const auto copier = [&](const CopyData &c) { 
 *       constraints.distribute_local_to_global(c.cell_matrix, 
 *                                              c.cell_rhs, 
 *                                              c.local_dof_indices, 
 *                                              system_matrix, 
 *                                              right_hand_side); 
 * 
 *       for (auto &cdf : c.face_data) 
 *         { 
 *           constraints.distribute_local_to_global(cdf.cell_matrix, 
 *                                                  cdf.joint_dof_indices, 
 *                                                  system_matrix); 
 *         } 
 *     }; 
 * 
 *     ScratchData<dim> scratch_data(mapping, fe, quadrature, quadrature_face); 
 *     CopyData         copy_data; 
 * 
 * @endcode
 * 
 * 在这里，我们最终处理了装配问题。我们传入ScratchData和CopyData对象，以及上面的lambda函数，并指定我们要对内部面进行一次装配。
 * 

 * 
 * 
 * @code
 *     MeshWorker::mesh_loop(dof_handler.begin_active(), 
 *                           dof_handler.end(), 
 *                           cell_worker, 
 *                           copier, 
 *                           scratch_data, 
 *                           copy_data, 
 *                           MeshWorker::assemble_own_cells | 
 *                             MeshWorker::assemble_boundary_faces | 
 *                             MeshWorker::assemble_own_interior_faces_once, 
 *                           boundary_worker, 
 *                           face_worker); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Alltherest"></a> 
 * <h3>All the rest</h3>
 * 

 * 
 * 对于这个简单的问题，我们使用了最简单的求解器，称为Richardson迭代，它代表了简单的缺陷修正。这与一个块状SSOR预处理相结合，该预处理使用DG离散化产生的系统矩阵的特殊块状结构。这些块的大小是每个单元的DoF数量。这里，我们使用SSOR预处理，因为我们没有根据流场对DoFs进行重新编号。如果在流的下游方向对DoFs进行重新编号，那么块状的Gauss-Seidel预处理（见PreconditionBlockSOR类，放松=1）会做得更好。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::solve() 
 *   { 
 *     SolverControl                    solver_control(1000, 1e-12); 
 *     SolverRichardson<Vector<double>> solver(solver_control); 
 * 
 * @endcode
 * 
 * 这里我们创建了预处理程序。
 * 

 * 
 * 
 * @code
 *     PreconditionBlockSSOR<SparseMatrix<double>> preconditioner; 
 * 
 * @endcode
 * 
 * 然后将矩阵分配给它，并设置正确的块大小。
 * 

 * 
 * 
 * @code
 *     preconditioner.initialize(system_matrix, fe.n_dofs_per_cell()); 
 * 
 * @endcode
 * 
 * 做完这些准备工作后，我们就可以启动线性求解器了。
 * 

 * 
 * 
 * @code
 *     solver.solve(system_matrix, solution, right_hand_side, preconditioner); 
 * 
 *     std::cout << "  Solver converged in " << solver_control.last_step() 
 *               << " iterations." << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 我们根据一个非常简单的细化标准来细化网格，即对解的梯度的近似。由于这里我们考虑的是DG(1)方法（即我们使用片状双线性形状函数），我们可以简单地计算每个单元的梯度。但是我们并不希望我们的细化指标只建立在每个单元的梯度上，而是希望同时建立在相邻单元之间的不连续解函数的跳跃上。最简单的方法是通过差分商计算近似梯度，包括考虑中的单元和其相邻的单元。这是由 <code>DerivativeApproximation</code> 类完成的，它计算近似梯度的方式类似于本教程 step-9 中描述的 <code>GradientEstimation</code> 。事实上， <code>DerivativeApproximation</code> 类是在 step-9 的 <code>GradientEstimation</code> 类之后开发的。与  step-9  中的讨论相关，这里我们考虑  $h^{1+d/2}|\nabla_h u_h|$  。此外，我们注意到，我们不考虑近似的二次导数，因为线性平流方程的解一般不在 $H^2$ 中，而只在 $H^1$ 中（或者，更准确地说：在 $H^1_\beta$ 中，即在方向 $\beta$ 上的导数是可平方整除的函数空间）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::refine_grid() 
 *   { 
 * 
 * @endcode
 * 
 * <code>DerivativeApproximation</code> 类将梯度计算为浮点精度。这已经足够了，因为它们是近似的，只作为细化指标。
 * 

 * 
 * 
 * @code
 *     Vector<float> gradient_indicator(triangulation.n_active_cells()); 
 * 
 * @endcode
 * 
 * 现在，近似梯度被计算出来了
 * 

 * 
 * 
 * @code
 *     DerivativeApproximation::approximate_gradient(mapping, 
 *                                                   dof_handler, 
 *                                                   solution, 
 *                                                   gradient_indicator); 
 * 
 * @endcode
 * 
 * 并且它们的单元格按系数 $h^{1+d/2}$ 进行缩放。
 * 
 * @code
 *     unsigned int cell_no = 0; 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       gradient_indicator(cell_no++) *= 
 *         std::pow(cell->diameter(), 1 + 1.0 * dim / 2); 
 * 
 * @endcode
 * 
 * 最后它们作为细化指标。
 * 

 * 
 * 
 * @code
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     gradient_indicator, 
 *                                                     0.3, 
 *                                                     0.1); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * 
 * @endcode
 * 
 * 这个程序的输出包括一个自适应细化网格的vtk文件和数值解。最后，我们还用 VectorTools::integrate_difference(). 计算了解的L-无穷大规范。
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     const std::string filename = "solution-" + std::to_string(cycle) + ".vtk"; 
 *     std::cout << "  Writing solution to <" << filename << ">" << std::endl; 
 *     std::ofstream output(filename); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler);  
 *     data_out.add_data_vector(solution, "u", DataOut<dim>::type_dof_data); 
 * 
 *     data_out.build_patches(mapping); 
 * 
 *     data_out.write_vtk(output); 
 * 
 *     { 
 *       Vector<float> values(triangulation.n_active_cells()); 
 *       VectorTools::integrate_difference(mapping, 
 *                                         dof_handler, 
 *                                         solution, 
 *                                         Functions::ZeroFunction<dim>(), 
 *                                         values, 
 *                                         quadrature, 
 *                                         VectorTools::Linfty_norm); 
 *       const double l_infty = 
 *         VectorTools::compute_global_error(triangulation, 
 *                                           values,  
 *                                           VectorTools::Linfty_norm); 
 *       std::cout << "  L-infinity norm: " << l_infty << std::endl; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 下面的 <code>run</code> 函数与前面的例子类似。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void AdvectionProblem<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 6; ++cycle) 
 *       {  
 *         std::cout << "Cycle " << cycle << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation); 
 *             triangulation.refine_global(3); 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 *         std::cout << "  Number of active cells:       " 
 *                   << triangulation.n_active_cells() << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         std::cout << "  Number of degrees of freedom: " << dof_handler.n_dofs() 
 *                   << std::endl; 
 * 
 *         assemble_system(); 
 *         solve(); 
 * 
 *         output_results(cycle); 
 *       } 
 *   } 
 * } // namespace Step12 
 * 
 * @endcode
 * 
 * 下面的 <code>main</code> 函数与前面的例子也类似，不需要注释。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step12::AdvectionProblem<2> dgmethod; 
 *       dgmethod.run(); 
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
 * 
 * 
 * @endcode
examples/step-12/doc/results.dox



<a name="Results"></a><h1>Results</h1>



这个程序的输出包括控制台输出和vtk格式的解决方案。

@code
Cycle 0
  Number of active cells:       64
  Number of degrees of freedom: 256
  Solver converged in 4 iterations.
  Writing solution to <solution-0.vtk>
  L-infinity norm: 1.09057
Cycle 1
  Number of active cells:       112
  Number of degrees of freedom: 448
  Solver converged in 9 iterations.
  Writing solution to <solution-1.vtk>
  L-infinity norm: 1.10402
Cycle 2
  Number of active cells:       214
  Number of degrees of freedom: 856
  Solver converged in 16 iterations.
  Writing solution to <solution-2.vtk>
  L-infinity norm: 1.09813
Cycle 3
  Number of active cells:       415
  Number of degrees of freedom: 1660
  Solver converged in 26 iterations.
  Writing solution to <solution-3.vtk>
  L-infinity norm: 1.09579
Cycle 4
  Number of active cells:       796
  Number of degrees of freedom: 3184
  Solver converged in 44 iterations.
  Writing solution to <solution-4.vtk>
  L-infinity norm: 1.09612
Cycle 5
  Number of active cells:       1561
  Number of degrees of freedom: 6244
  Solver converged in 81 iterations.
  Writing solution to <solution-5.vtk>
@endcode



我们展示了初始网格的解决方案，以及经过两个和五个自适应细化步骤后的网格。

 <img src="https://www.dealii.org/images/steps/developer/step-12.sol-0.png" alt="">   <img src="https://www.dealii.org/images/steps/developer/step-12.sol-2.png" alt="">   <img src="https://www.dealii.org/images/steps/developer/step-12.sol-5.png" alt=""> 。

最后我们展示一个3D计算的图。

 <img src="https://www.dealii.org/images/steps/developer/step-12.sol-5-3d.png" alt=""> 


<a name="dg-vs-cg"></a>

<a name="Whyusediscontinuouselements"></a><h3>Why use discontinuous elements</h3>


在这个程序中，我们使用了不连续的元素。这是一个合理的问题，为什么不简单地使用正常的、连续的元素呢？当然，对于每个有数值方法背景的人来说，答案是显而易见的：连续Galerkin（cG）方法对于输运方程是不稳定的，除非特别添加稳定项。然而，DG方法<i>is</i>则是稳定的。用目前的程序来说明这一点并不十分困难；事实上，只需要做以下的小修改就可以了。

- 将该元素改为FE_Q，而不是FE_DGQ。

- 以与步骤6完全相同的方式增加对悬挂节点约束的处理。

- 我们需要一个不同的求解器；步骤29中的直接求解器是一个方便的选择。一个有经验的deal.II用户将能够在10分钟内完成这一工作。

虽然上面显示了二维解决方案，在界面上含有一些小的尖峰，但是在网格细化的情况下，这些尖峰的高度是稳定的，当使用连续元素时，结果看起来有很大不同。

 <table align="center">
  <tr>
    <td valign="top">
      0 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-0.png" alt="">
    </td>
    <td valign="top">
      1 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-1.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      2 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-2.png" alt="">
    </td>
    <td valign="top">
      3 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-3.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      4 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-4.png" alt="">
    </td>
    <td valign="top">
      5 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-12.cg.sol-5.png" alt="">
    </td>
  </tr>
</table> 

在细化迭代5中，图像不能再以合理的方式绘制成三维图。因此我们展示了一个范围为 $[-1,2]$ 的彩色图（当然，精确解的解值位于 $[0,1]$ ）。在任何情况下，很明显，连续Galerkin解表现出振荡行为，随着网格的细化越来越差。

如果出于某种原因想使用连续元素，有很多策略可以稳定cG方法。讨论这些方法超出了本教程程序的范围；例如，感兴趣的读者可以看看步骤31。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


鉴于在这种情况下确切的解是已知的，进一步扩展的一个有趣的途径是确认这个程序的收敛顺序。在目前的情况下，解是非光滑的，因此我们不能期望得到特别高的收敛阶，即使我们使用高阶元素。但是，即使解<i>is</i>光滑，方程也不是椭圆的，因此，我们应该得到等于最优插值估计的收敛阶，这一点并不明显（例如，我们通过使用二次元会得到 $h^3$ 在 $L^2$ 准则下的收敛）。

事实上，对于双曲方程来说，理论预测常常表明，我们所能希望的最好结果是比插值估计值低二分之一的阶。例如，对于流线扩散法（此处用于稳定传输方程解的DG法的替代方法），可以证明对于度数为 $p$ 的元素，在任意网格上的收敛阶为 $p+\frac 12$ 。虽然在均匀细化的网格上观察到的顺序经常是 $p+1$ ，但人们可以构建所谓的彼得森网格，在这些网格上实际上达到了更差的理论约束。这应该是比较容易验证的，例如使用 VectorTools::integrate_difference 函数。

一个不同的方向是观察运输问题的解决经常有不连续性，因此，我们<i>bisect</i>在每个坐标方向上的每个单元的网格可能不是最佳的。相反，一个更好的策略是只在平行于不连续的方向上切割单元。这被称为<i>anisotropic mesh refinement</i>，是步骤30的主题。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-12.cc"
*/
