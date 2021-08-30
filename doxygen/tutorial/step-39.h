/**
@page step_39 The step-39 tutorial program
This tutorial depends on step-12b.

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
        <li><a href="#Thelocalintegrators">The local integrators</a>
        <li><a href="#Themainclass">The main class</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Logfileoutput">Logfile output</a>
        <li><a href="#Postprocessingofthelogfile">Postprocessing of the logfile</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-39/doc/intro.dox

<a name="Intro"></a>

在这个程序中，我们使用内部惩罚方法和Nitsche的弱边界条件来解决Poisson方程。我们在局部细化的网格上使用多网格方法，这些网格是用一个体块准则和一个基于单元和面残差的标准误差估计器生成的。所有的运算符都是用MeshWorker接口实现的。

像步骤12一样，离散化依赖于有限元空间，它在网格单元 $K\in \mathbb T_h$ 内是多项式的，但在单元之间没有连续性。由于这种函数在每个内部面 $F\in \mathbb F_h^i$ 上有两个值，每边一个，我们定义均值和跳跃算子如下：让<i>K</i><sub>1</sub>和<i>K</i><sub>2</sub>是共享一个面的两个单元，让函数的轨迹<i>u<sub>i</sub></i>和外法向量<b>n</b><i><sub>i</sub></i>相应地被标记。然后，在这个面上，我们让

@f[
	\average{ u } = \frac{u_1 + u_2}2


@f]



注意，如果这样的表达式包含一个法向量，那么平均运算符就会变成一个跳跃。该问题的内部惩罚方法

@f[


  -\Delta u = f \text{ in }\Omega \qquad u = u^D \text{ on } \partial\Omega


@f]

成为

@f{multline*}
  \sum_{K\in \mathbb T_h} (\nabla u, \nabla v)_K
  \\
  + \sum_{F \in F_h^i} \biggl\{4\sigma_F (\average{ u \mathbf n}, \average{ v \mathbf n })_F


  - 2 (\average{ \nabla u },\average{ v\mathbf n })_F


  - 2 (\average{ \nabla v },\average{ u\mathbf n })_F
  \biggr\}
  \\
  + \sum_{F \in F_h^b} \biggl\{2\sigma_F (u, v)_F


  - (\partial_n u,v)_F


  - (\partial_n v,u)_F
  \biggr\}
  \\
  = (f, v)_\Omega + \sum_{F \in F_h^b} \biggl\{
  2\sigma_F (u^D, v)_F - (\partial_n v,u^D)_F
  \biggr\}.


@f}



这里， $\sigma_F$ 是惩罚参数，其选择如下：对于<i>F</i>单元格<i>K</i>的一个面，计算数值

@f[
\sigma_{F,K} = p(p+1) \frac{|F|_{d-1}}{|K|_d},


@f]

其中<i>p</i>是有限元函数的多项式程度， $|\cdot|_d$ 和 $|\cdot|_{d-1}$ 表示相应对象的 $d$ 和 $d-1$ 维度的Hausdorff度量。如果面在边界上，选择 $\sigma_F = \sigma_{F,K}$  。对于一个内部的面，我们取这个面的两个值的平均值。

在我们的有限元程序中，我们区分了三种不同的积分，分别对应于上面的单元、内部面和边界面的总和。由于 MeshWorker::loop 为我们组织了这些和，我们只需要实现对每个网格元素的积分。下面的MatrixIntegrator类有这三个函数用于公式的左边，RHSIntegrator类用于右边。

正如我们将在下面看到的，甚至误差估计也是相同的结构，因为它可以写成

@f{align*}
  \eta^2 &= \eta_K^2 + \eta_F^2 + \eta_B^2
  \\
  \eta_K^2 &= \sum_{K\in \mathbb T_h} h^2 \|f + \Delta u_h\|^2
  \\
  \eta_F^2 &= \sum_{F \in F_h^i} \biggl\{
    4 \sigma_F \| \average{u_h\mathbf n} \|^2 + h \|\average{\partial_n u_h}\|^2 \biggr\}
  \\
  \eta_B^2 &= \sum_{F \in F_h^b} 2\sigma_F \| u_h-u^D \|^2.


@f}



因此，下面用于组装矩阵、右手和误差估计的函数显示，这些循环都是通用的，可以用同样的方式进行编程。

这个程序与步骤12b有关，因为它使用MeshWorker和非连续Galerkin方法。在那里，我们解决的是一个平流问题，而这里是一个扩散问题。在这里，我们还使用了多网格预处理和一个理论上合理的误差估计器，见Karakashian和Pascal（2003）。Kanschat (2004)详细讨论了多层次方案。Hoppe, Kanschat, and Warburton (2009)讨论了自适应迭代及其收敛性（对于三角形网格）。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 线性代数的包含文件。一个普通的SparseMatrix，它又将包括SparsityPattern和Vector类的必要文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/precondition_block.h> 
 * #include <deal.II/lac/block_vector.h> 
 * 
 * @endcode
 * 
 * 包括用于设置网格的文件
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * @endcode
 * 
 * FiniteElement类和DoFHandler的包含文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_dgp.h> 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * @endcode
 * 
 * 使用MeshWorker框架的包含文件
 * 

 * 
 * 
 * @code
 * #include <deal.II/meshworker/dof_info.h> 
 * #include <deal.II/meshworker/integration_info.h> 
 * #include <deal.II/meshworker/assembler.h> 
 * #include <deal.II/meshworker/loop.h> 
 * 
 * @endcode
 * 
 * 与拉普拉斯相关的局部积分器的包含文件
 * 

 * 
 * 
 * @code
 * #include <deal.II/integrators/laplace.h> 
 * 
 * @endcode
 * 
 * 支持多网格方法
 * 

 * 
 * 
 * @code
 * #include <deal.II/multigrid/mg_tools.h> 
 * #include <deal.II/multigrid/multigrid.h> 
 * #include <deal.II/multigrid/mg_matrix.h> 
 * #include <deal.II/multigrid/mg_transfer.h> 
 * #include <deal.II/multigrid/mg_coarse.h> 
 * #include <deal.II/multigrid/mg_smoother.h> 
 * 
 * @endcode
 * 
 * 最后，我们从库中取出我们的精确解，以及正交和附加工具。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function_lib.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * 
 * @endcode
 * 
 * deal.II库的所有类都在dealii命名空间中。为了节省打字，我们告诉编译器也要在其中搜索名字。
 * 

 * 
 * 
 * @code
 * namespace Step39 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 这是我们用来设置边界值的函数，也是我们比较的精确解。
 * 

 * 
 * 
 * @code
 *   Functions::SlitSingularityFunction<2> exact_solution; 
 * @endcode
 * 
 * 
 * <a name="Thelocalintegrators"></a> 
 * <h3>The local integrators</h3>
 * 

 * 
 * MeshWorker将局部积分与单元格和面的循环分离开来。因此，我们必须编写局部积分类来生成矩阵、右手边和误差估计器。
 * 

 * 
 * 所有这些类都有相同的三个函数，分别用于对单元、边界面和内部面的积分。局部积分所需的所有信息都由 MeshWorker::IntegrationInfo<dim>. 提供。请注意，函数的签名不能改变，因为它是由 MeshWorker::integration_loop(). 所期望的。
 * 

 * 
 * 第一个定义局部积分器的类负责计算单元和面矩阵。它被用来组装全局矩阵以及水平矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class MatrixIntegrator : public MeshWorker::LocalIntegrator<dim> 
 *   { 
 *   public: 
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void 
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *               MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *               typename MeshWorker::IntegrationInfo<dim> &info1, 
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 在每个单元上，我们对Dirichlet形式进行积分。我们使用LocalIntegrators中的现成积分库来避免自己编写这些循环。同样地，我们实现了Nitsche边界条件和单元间的内部惩罚通量。
 * 

 * 
 * 边界和通量项需要一个惩罚参数，这个参数应该根据单元的大小和多项式的度数来调整。在 LocalIntegrators::Laplace::compute_penalty() 中可以找到关于这个参数的安全选择，我们在下面使用这个参数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MatrixIntegrator<dim>::cell( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0, false).matrix, 
 *                                            info.fe_values()); 
 *   } 
 * 
 *   template <int dim> 
 *   void MatrixIntegrator<dim>::boundary( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const unsigned int degree = info.fe_values(0).get_fe().tensor_degree(); 
 *     LocalIntegrators::Laplace::nitsche_matrix( 
 *       dinfo.matrix(0, false).matrix, 
 *       info.fe_values(0), 
 *       LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, degree, degree)); 
 *   } 
 * 
 * @endcode
 * 
 * 内部面使用内部惩罚方法
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MatrixIntegrator<dim>::face( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *     MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *     typename MeshWorker::IntegrationInfo<dim> &info1, 
 *     typename MeshWorker::IntegrationInfo<dim> &info2) const 
 *   { 
 *     const unsigned int degree = info1.fe_values(0).get_fe().tensor_degree(); 
 *     LocalIntegrators::Laplace::ip_matrix( 
 *       dinfo1.matrix(0, false).matrix, 
 *       dinfo1.matrix(0, true).matrix, 
 *       dinfo2.matrix(0, true).matrix, 
 *       dinfo2.matrix(0, false).matrix, 
 *       info1.fe_values(0), 
 *       info2.fe_values(0), 
 *       LocalIntegrators::Laplace::compute_penalty( 
 *         dinfo1, dinfo2, degree, degree)); 
 *   } 
 * 
 * @endcode
 * 
 * 第二个局部积分器建立了右手边。在我们的例子中，右手边的函数为零，这样，这里只设置了弱形式的边界条件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RHSIntegrator : public MeshWorker::LocalIntegrator<dim> 
 *   { 
 *   public: 
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void 
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *               MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *               typename MeshWorker::IntegrationInfo<dim> &info1, 
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   void 
 *   RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &, 
 *                            typename MeshWorker::IntegrationInfo<dim> &) const 
 *   {} 
 * 
 *   template <int dim> 
 *   void RHSIntegrator<dim>::boundary( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const FEValuesBase<dim> &fe           = info.fe_values(); 
 *     Vector<double> &         local_vector = dinfo.vector(0).block(0); 
 * 
 *     std::vector<double> boundary_values(fe.n_quadrature_points); 
 *     exact_solution.value_list(fe.get_quadrature_points(), boundary_values); 
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree(); 
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
 *                            dinfo.cell->measure(); 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       for (unsigned int i = 0; i < fe.dofs_per_cell; ++i) 
 *         local_vector(i) += 
 *           (-penalty * fe.shape_value(i, k)              // (-sigma * v_i(x_k) 
 *            + fe.normal_vector(k) * fe.shape_grad(i, k)) // + n * grad v_i(x_k)) 
 *           * boundary_values[k] * fe.JxW(k);             // u^D(x_k) * dx 
 *   } 
 * 
 *   template <int dim> 
 *   void 
 *   RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &, 
 *                            MeshWorker::DoFInfo<dim> &, 
 *                            typename MeshWorker::IntegrationInfo<dim> &, 
 *                            typename MeshWorker::IntegrationInfo<dim> &) const 
 *   {} 
 * 
 * @endcode
 * 
 * 第三个局部积分器负责对误差估计的贡献。这是由Karakashian和Pascal（2003）提出的标准能量估计器。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Estimator : public MeshWorker::LocalIntegrator<dim> 
 *   { 
 *   public: 
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void 
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *               MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *               typename MeshWorker::IntegrationInfo<dim> &info1, 
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 单元的贡献是离散解的拉普拉斯，因为右手边是零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *                        typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const FEValuesBase<dim> &fe = info.fe_values(); 
 * 
 *     const std::vector<Tensor<2, dim>> &DDuh = info.hessians[0][0]; 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         const double t = dinfo.cell->diameter() * trace(DDuh[k]); 
 *         dinfo.value(0) += t * t * fe.JxW(k); 
 *       } 
 *     dinfo.value(0) = std::sqrt(dinfo.value(0)); 
 *   } 
 * 
 * @endcode
 * 
 * 在边界，我们简单地使用边界残差的加权形式，即有限元解和正确边界条件之间的差值的规范。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Estimator<dim>::boundary( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const FEValuesBase<dim> &fe = info.fe_values(); 
 * 
 *     std::vector<double> boundary_values(fe.n_quadrature_points); 
 *     exact_solution.value_list(fe.get_quadrature_points(), boundary_values); 
 * 
 *     const std::vector<double> &uh = info.values[0][0]; 
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree(); 
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
 *                            dinfo.cell->measure(); 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         const double diff = boundary_values[k] - uh[k]; 
 *         dinfo.value(0) += penalty * diff * diff * fe.JxW(k); 
 *       } 
 *     dinfo.value(0) = std::sqrt(dinfo.value(0)); 
 *   } 
 * 
 * @endcode
 * 
 * 最后，在内部面，估计器由解的跳跃和它的法向导数组成，并进行适当的加权。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   Estimator<dim>::face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *                        MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *                        typename MeshWorker::IntegrationInfo<dim> &info1, 
 *                        typename MeshWorker::IntegrationInfo<dim> &info2) const 
 *   { 
 *     const FEValuesBase<dim> &          fe   = info1.fe_values(); 
 *     const std::vector<double> &        uh1  = info1.values[0][0]; 
 *     const std::vector<double> &        uh2  = info2.values[0][0]; 
 *     const std::vector<Tensor<1, dim>> &Duh1 = info1.gradients[0][0]; 
 *     const std::vector<Tensor<1, dim>> &Duh2 = info2.gradients[0][0]; 
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree(); 
 *     const double       penalty1 = 
 *       degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure(); 
 *     const double penalty2 = 
 *       degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure(); 
 *     const double penalty = penalty1 + penalty2; 
 *     const double h       = dinfo1.face->measure(); 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         const double diff1 = uh1[k] - uh2[k]; 
 *         const double diff2 = 
 *           fe.normal_vector(k) * Duh1[k] - fe.normal_vector(k) * Duh2[k]; 
 *         dinfo1.value(0) += 
 *           (penalty * diff1 * diff1 + h * diff2 * diff2) * fe.JxW(k); 
 *       } 
 *     dinfo1.value(0) = std::sqrt(dinfo1.value(0)); 
 *     dinfo2.value(0) = dinfo1.value(0); 
 *   } 
 * 
 * @endcode
 * 
 * 最后我们有一个误差的积分器。由于不连续Galerkin问题的能量准则不仅涉及到单元内部的梯度差，还涉及到跨面和边界的跳跃项，所以我们不能仅仅使用  VectorTools::integrate_difference().  而是使用MeshWorker接口来自己计算误差。
 * 

 * 
 * 有几种不同的方法来定义这个能量准则，但是所有的方法都是随着网格大小的变化而等价的（有些不是随着多项式程度的变化而等价）。这里，我们选择
 * @f[ \|u\|_{1,h} =
 * \sum_{K\in \mathbb T_h} \|\nabla u\|_K^2 + \sum_{F \in F_h^i}
 * 4\sigma_F\|\average{ u \mathbf n}\|^2_F + \sum_{F \in F_h^b}
 * 2\sigma_F\|u\|^2_F 
 * @f]
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim> 
 *   { 
 *   public: 
 *     void cell(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *               typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void 
 *          boundary(MeshWorker::DoFInfo<dim> &                 dinfo, 
 *                   typename MeshWorker::IntegrationInfo<dim> &info) const override; 
 *     void face(MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *               MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *               typename MeshWorker::IntegrationInfo<dim> &info1, 
 *               typename MeshWorker::IntegrationInfo<dim> &info2) const override; 
 *   }; 
 * 
 * @endcode
 * 
 * 这里我们有关于单元格的集成。目前MeshWorker中还没有很好的接口可以让我们访问正交点中的正则函数值。因此，我们必须在单元格积分器中创建精确函数值和梯度的向量。之后，一切照旧，我们只需将差值的平方加起来。
 * 

 * 
 * 除了计算能量准则的误差，我们还利用网格工作者的能力同时计算两个函数并在同一个循环中计算<i>L<sup>2</sup></i>的误差。很明显，这个函数没有任何跳跃项，只出现在单元格的积分中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ErrorIntegrator<dim>::cell( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const FEValuesBase<dim> &   fe = info.fe_values(); 
 *     std::vector<Tensor<1, dim>> exact_gradients(fe.n_quadrature_points); 
 *     std::vector<double>         exact_values(fe.n_quadrature_points); 
 * 
 *     exact_solution.gradient_list(fe.get_quadrature_points(), exact_gradients); 
 *     exact_solution.value_list(fe.get_quadrature_points(), exact_values); 
 * 
 *     const std::vector<Tensor<1, dim>> &Duh = info.gradients[0][0]; 
 *     const std::vector<double> &        uh  = info.values[0][0]; 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         double sum = 0; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           { 
 *             const double diff = exact_gradients[k][d] - Duh[k][d]; 
 *             sum += diff * diff; 
 *           } 
 *         const double diff = exact_values[k] - uh[k]; 
 *         dinfo.value(0) += sum * fe.JxW(k); 
 *         dinfo.value(1) += diff * diff * fe.JxW(k); 
 *       } 
 *     dinfo.value(0) = std::sqrt(dinfo.value(0)); 
 *     dinfo.value(1) = std::sqrt(dinfo.value(1)); 
 *   } 
 * 
 *   template <int dim> 
 *   void ErrorIntegrator<dim>::boundary( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo, 
 *     typename MeshWorker::IntegrationInfo<dim> &info) const 
 *   { 
 *     const FEValuesBase<dim> &fe = info.fe_values(); 
 * 
 *     std::vector<double> exact_values(fe.n_quadrature_points); 
 *     exact_solution.value_list(fe.get_quadrature_points(), exact_values); 
 * 
 *     const std::vector<double> &uh = info.values[0][0]; 
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree(); 
 *     const double penalty = 2. * degree * (degree + 1) * dinfo.face->measure() / 
 *                            dinfo.cell->measure(); 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         const double diff = exact_values[k] - uh[k]; 
 *         dinfo.value(0) += penalty * diff * diff * fe.JxW(k); 
 *       } 
 *     dinfo.value(0) = std::sqrt(dinfo.value(0)); 
 *   } 
 * 
 *   template <int dim> 
 *   void ErrorIntegrator<dim>::face( 
 *     MeshWorker::DoFInfo<dim> &                 dinfo1, 
 *     MeshWorker::DoFInfo<dim> &                 dinfo2, 
 *     typename MeshWorker::IntegrationInfo<dim> &info1, 
 *     typename MeshWorker::IntegrationInfo<dim> &info2) const 
 *   { 
 *     const FEValuesBase<dim> &  fe  = info1.fe_values(); 
 *     const std::vector<double> &uh1 = info1.values[0][0]; 
 *     const std::vector<double> &uh2 = info2.values[0][0]; 
 * 
 *     const unsigned int degree = fe.get_fe().tensor_degree(); 
 *     const double       penalty1 = 
 *       degree * (degree + 1) * dinfo1.face->measure() / dinfo1.cell->measure(); 
 *     const double penalty2 = 
 *       degree * (degree + 1) * dinfo2.face->measure() / dinfo2.cell->measure(); 
 *     const double penalty = penalty1 + penalty2; 
 * 
 *     for (unsigned k = 0; k < fe.n_quadrature_points; ++k) 
 *       { 
 *         const double diff = uh1[k] - uh2[k]; 
 *         dinfo1.value(0) += (penalty * diff * diff) * fe.JxW(k); 
 *       } 
 *     dinfo1.value(0) = std::sqrt(dinfo1.value(0)); 
 *     dinfo2.value(0) = dinfo1.value(0); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Themainclass"></a> 
 * <h3>The main class</h3>
 * 

 * 
 * 这个类做主要的工作，就像前面的例子一样。关于这里声明的函数的描述，请参考下面的实现。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class InteriorPenaltyProblem 
 *   { 
 *   public: 
 *     using CellInfo = MeshWorker::IntegrationInfo<dim>; 
 * 
 *     InteriorPenaltyProblem(const FiniteElement<dim> &fe); 
 * 
 *     void run(unsigned int n_steps); 
 * 
 *   private: 
 *     void   setup_system(); 
 *     void   assemble_matrix(); 
 *     void   assemble_mg_matrix(); 
 *     void   assemble_right_hand_side(); 
 *     void   error(); 
 *     double estimate(); 
 *     void   solve(); 
 *     void   output_results(const unsigned int cycle) const; 
 * 
 * @endcode
 * 
 * 与离散化有关的成员对象在这里。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim>        triangulation; 
 *     const MappingQ1<dim>      mapping; 
 *     const FiniteElement<dim> &fe; 
 *     DoFHandler<dim>           dof_handler; 
 * 
 * @endcode
 * 
 * 然后，我们有与全局离散系统相关的矩阵和向量。
 * 

 * 
 * 
 * @code
 *     SparsityPattern      sparsity; 
 *     SparseMatrix<double> matrix; 
 *     Vector<double>       solution; 
 *     Vector<double>       right_hand_side; 
 *     BlockVector<double>  estimates; 
 * 
 * @endcode
 * 
 * 最后，我们有一组与多级预处理程序相关的稀疏模式和稀疏矩阵。 首先，我们有一个水平矩阵和它的稀疏性模式。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<SparsityPattern>      mg_sparsity; 
 *     MGLevelObject<SparseMatrix<double>> mg_matrix; 
 * 
 * @endcode
 * 
 * 当我们在局部细化的网格上进行局部平滑的多重网格时，需要额外的矩阵；见Kanschat（2004）。这里是这些边缘矩阵的稀疏性模式。我们只需要一个，因为上矩阵的模式是下矩阵的转置。实际上，我们并不太关心这些细节，因为MeshWorker正在填充这些矩阵。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<SparsityPattern> mg_sparsity_dg_interface; 
 * 
 * @endcode
 * 
 * 精细化边缘的通量矩阵，将精细级自由度与粗略级自由度相耦合。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<SparseMatrix<double>> mg_matrix_dg_down; 
 * 
 * @endcode
 * 
 * 精细化边缘的通量矩阵的转置，将粗级自由度耦合到精细级。
 * 

 * 
 * 
 * @code
 *     MGLevelObject<SparseMatrix<double>> mg_matrix_dg_up; 
 *   }; 
 * 
 * @endcode
 * 
 * 构造函数简单地设置了粗略的网格和DoFHandler。FiniteElement作为一个参数被提供，以实现灵活性。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   InteriorPenaltyProblem<dim>::InteriorPenaltyProblem( 
 *     const FiniteElement<dim> &fe) 
 *     : triangulation(Triangulation<dim>::limit_level_difference_at_vertices) 
 *     , mapping() 
 *     , fe(fe) 
 *     , dof_handler(triangulation) 
 *     , estimates(1) 
 *   { 
 *     GridGenerator::hyper_cube_slit(triangulation, -1, 1); 
 *   } 
 * 
 * @endcode
 * 
 * 在这个函数中，我们设置了线性系统的维度和全局矩阵以及水平矩阵的稀疏性模式。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::setup_system() 
 *   { 
 * 
 * @endcode
 * 
 * 首先，我们用有限元将自由度分布在网格上并对其进行编号。
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 *     dof_handler.distribute_mg_dofs(); 
 *     unsigned int n_dofs = dof_handler.n_dofs(); 
 * 
 * @endcode
 * 
 * 然后，我们已经知道代表有限元函数的向量的大小。
 * 

 * 
 * 
 * @code
 *     solution.reinit(n_dofs); 
 *     right_hand_side.reinit(n_dofs); 
 * 
 * @endcode
 * 
 * 接下来，我们为全局矩阵设置稀疏性模式。由于我们事先不知道行的大小，所以我们首先填充一个临时的DynamicSparsityPattern对象，一旦完成，就将其复制到常规的SparsityPattern中。
 * 

 * 
 * 
 * @code
 *     DynamicSparsityPattern dsp(n_dofs); 
 *     DoFTools::make_flux_sparsity_pattern(dof_handler, dsp); 
 *     sparsity.copy_from(dsp); 
 *     matrix.reinit(sparsity); 
 * 
 *     const unsigned int n_levels = triangulation.n_levels(); 
 * 
 * @endcode
 * 
 * 全局系统已经设置好了，现在我们来关注一下级别矩阵。我们调整所有矩阵对象的大小，以便每一级都有一个矩阵。
 * 

 * 
 * 
 * @code
 *     mg_matrix.resize(0, n_levels - 1); 
 *     mg_matrix.clear_elements(); 
 *     mg_matrix_dg_up.resize(0, n_levels - 1); 
 *     mg_matrix_dg_up.clear_elements(); 
 *     mg_matrix_dg_down.resize(0, n_levels - 1); 
 *     mg_matrix_dg_down.clear_elements(); 
 * 
 * @endcode
 * 
 * 在为水平矩阵调用<tt>clear()</tt>之后更新稀疏模式很重要，因为矩阵通过SmartPointer和Subscriptor机制锁定了稀疏模式。
 * 

 * 
 * 
 * @code
 *     mg_sparsity.resize(0, n_levels - 1); 
 *     mg_sparsity_dg_interface.resize(0, n_levels - 1); 
 * 
 * @endcode
 * 
 * 现在，所有的对象都准备好了，可以在每一层容纳一个稀疏模式或矩阵。剩下的就是在每一层设置稀疏模式了。
 * 

 * 
 * 
 * @code
 *     for (unsigned int level = mg_sparsity.min_level(); 
 *          level <= mg_sparsity.max_level(); 
 *          ++level) 
 *       { 
 * 
 * @endcode
 * 
 * 这些与上面的全局矩阵的行数大致相同，现在是每个级别的。
 * 

 * 
 * 
 * @code
 *         DynamicSparsityPattern dsp(dof_handler.n_dofs(level)); 
 *         MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level); 
 *         mg_sparsity[level].copy_from(dsp); 
 *         mg_matrix[level].reinit(mg_sparsity[level]); 
 * 
 * @endcode
 * 
 * 另外，我们需要初始化各层之间细化边缘的转移矩阵。它们被存储在两个索引中较细的索引处，因此在0层没有这样的对象。
 * 

 * 
 * 
 * @code
 *         if (level > 0) 
 *           { 
 *             DynamicSparsityPattern dsp; 
 *             dsp.reinit(dof_handler.n_dofs(level - 1), 
 *                        dof_handler.n_dofs(level)); 
 *             MGTools::make_flux_sparsity_pattern_edge(dof_handler, dsp, level); 
 *             mg_sparsity_dg_interface[level].copy_from(dsp); 
 *             mg_matrix_dg_up[level].reinit(mg_sparsity_dg_interface[level]); 
 *             mg_matrix_dg_down[level].reinit(mg_sparsity_dg_interface[level]); 
 *           } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 在这个函数中，我们组装全局系统矩阵，这里的全局是指我们解决的离散系统的矩阵，它覆盖了整个网格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::assemble_matrix() 
 *   { 
 * 
 * @endcode
 * 
 * 首先，我们需要设置提供我们集成值的对象。这个对象包含了所有需要的FEValues和FEFaceValues对象，并且自动维护它们，使它们总是指向当前单元。为此，我们首先需要告诉它，在哪里计算，计算什么。由于我们没有做任何花哨的事情，我们可以依靠他们对正交规则的标准选择。
 * 

 * 
 * 由于他们的默认更新标志是最小的，我们另外添加我们需要的东西，即所有对象（单元格、边界和内部面）上的形状函数的值和梯度。之后，我们准备初始化容器，它将创建所有必要的FEValuesBase对象进行整合。
 * 

 * 
 * 
 * @code
 *     MeshWorker::IntegrationInfoBox<dim> info_box; 
 *     UpdateFlags update_flags = update_values | update_gradients; 
 *     info_box.add_update_flags_all(update_flags); 
 *     info_box.initialize(fe, mapping); 
 * 
 * @endcode
 * 
 * 这就是我们整合本地数据的对象。它由MatrixIntegrator中的局部整合例程填充，然后由汇编器用来将信息分配到全局矩阵中。
 * 

 * 
 * 
 * @code
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler); 
 * 
 * @endcode
 * 
 * 此外，我们还需要一个将局部矩阵装配到全局矩阵的对象。这些装配器对象拥有目标对象结构的所有知识，在这里是一个稀疏矩阵，可能的约束和网格结构。
 * 

 * 
 * 
 * @code
 *     MeshWorker::Assembler::MatrixSimple<SparseMatrix<double>> assembler; 
 *     assembler.initialize(matrix); 
 * 
 * @endcode
 * 
 * 现在是我们自己编码的部分，局部积分器。这是唯一与问题有关的部分。
 * 

 * 
 * 
 * @code
 *     MatrixIntegrator<dim> integrator; 
 * 
 * @endcode
 * 
 * 现在，我们把所有的东西都扔到 MeshWorker::loop(), 中，在这里遍历网格的所有活动单元，计算单元和面的矩阵，并把它们集合到全局矩阵中。我们在这里使用变量<tt>dof_handler</tt>，以便使用全局自由度的编号。
 * 

 * 
 * 
 * @code
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
 *                                            dof_handler.end(), 
 *                                            dof_info, 
 *                                            info_box, 
 *                                            integrator, 
 *                                            assembler); 
 *   } 
 * 
 * @endcode
 * 
 * 现在，我们对水平矩阵做同样的处理。不太令人惊讶的是，这个函数看起来像前一个函数的孪生兄弟。事实上，只有两个小的区别。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::assemble_mg_matrix() 
 *   { 
 *     MeshWorker::IntegrationInfoBox<dim> info_box; 
 *     UpdateFlags update_flags = update_values | update_gradients; 
 *     info_box.add_update_flags_all(update_flags); 
 *     info_box.initialize(fe, mapping); 
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler); 
 * 
 * @endcode
 * 
 * 很明显，需要用一个填充水平矩阵的汇编器来代替。请注意，它也会自动填充边缘矩阵。
 * 

 * 
 * 
 * @code
 *     MeshWorker::Assembler::MGMatrixSimple<SparseMatrix<double>> assembler; 
 *     assembler.initialize(mg_matrix); 
 *     assembler.initialize_fluxes(mg_matrix_dg_up, mg_matrix_dg_down); 
 * 
 *     MatrixIntegrator<dim> integrator; 
 * 
 * @endcode
 * 
 * 这里是与前一个函数的另一个不同之处：我们在所有单元上运行，而不仅仅是活动单元。而且我们使用以 <code>_mg</code> 结尾的函数，因为我们需要每一层的自由度，而不是全局的编号。
 * 

 * 
 * 
 * @code
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_mg(), 
 *                                            dof_handler.end_mg(), 
 *                                            dof_info, 
 *                                            info_box, 
 *                                            integrator, 
 *                                            assembler); 
 *   } 
 * 
 * @endcode
 * 
 * 这里我们有另一个assemble函数的克隆。与组装系统矩阵的区别在于，我们在这里组装了一个向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::assemble_right_hand_side() 
 *   { 
 *     MeshWorker::IntegrationInfoBox<dim> info_box; 
 *     UpdateFlags                         update_flags = 
 *       update_quadrature_points | update_values | update_gradients; 
 *     info_box.add_update_flags_all(update_flags); 
 *     info_box.initialize(fe, mapping); 
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler); 
 * 
 * @endcode
 * 
 * 因为这个汇编器允许我们填充多个向量，所以接口要比上面复杂一些。向量的指针必须存储在一个AnyData对象中。虽然这在这里似乎造成了两行额外的代码，但实际上在更复杂的应用中它是很方便的。
 * 

 * 
 * 
 * @code
 *     MeshWorker::Assembler::ResidualSimple<Vector<double>> assembler; 
 *     AnyData                                               data; 
 *     data.add<Vector<double> *>(&right_hand_side, "RHS"); 
 *     assembler.initialize(data); 
 * 
 *     RHSIntegrator<dim> integrator; 
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
 *                                            dof_handler.end(), 
 *                                            dof_info, 
 *                                            info_box, 
 *                                            integrator, 
 *                                            assembler); 
 * 
 *     right_hand_side *= -1.; 
 *   } 
 * 
 * @endcode
 * 
 * 现在，我们已经对构建离散线性系统的所有函数进行了编码，现在是我们实际解决它的时候了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::solve() 
 *   { 
 * 
 * @endcode
 * 
 * 选择的求解器是共轭梯度。
 * 

 * 
 * 
 * @code
 *     SolverControl            control(1000, 1.e-12); 
 *     SolverCG<Vector<double>> solver(control); 
 * 
 * @endcode
 * 
 * 现在我们正在设置多级预处理程序的组件。首先，我们需要在网格层之间进行转移。我们在这里使用的对象为这些转移生成了稀疏矩阵。
 * 

 * 
 * 
 * @code
 *     MGTransferPrebuilt<Vector<double>> mg_transfer; 
 *     mg_transfer.build(dof_handler); 
 * 
 * @endcode
 * 
 * 然后，我们需要一个精确的解算器来解算最粗层次上的矩阵。
 * 

 * 
 * 
 * @code
 *     FullMatrix<double> coarse_matrix; 
 *     coarse_matrix.copy_from(mg_matrix[0]); 
 *     MGCoarseGridHouseholder<double, Vector<double>> mg_coarse; 
 *     mg_coarse.initialize(coarse_matrix); 
 * 
 * @endcode
 * 
 * 虽然转移和粗略网格求解器几乎是通用的，但为平滑器提供了更多的灵活性。首先，我们选择Gauss-Seidel作为我们的平滑方法。
 * 

 * 
 * 
 * @code
 *     GrowingVectorMemory<Vector<double>> mem; 
 *     using RELAXATION = PreconditionSOR<SparseMatrix<double>>; 
 *     mg::SmootherRelaxation<RELAXATION, Vector<double>> mg_smoother; 
 *     RELAXATION::AdditionalData                         smoother_data(1.); 
 *     mg_smoother.initialize(mg_matrix, smoother_data); 
 * 
 * @endcode
 * 
 * 在每个级别上做两个平滑步骤。
 * 

 * 
 * 
 * @code
 *     mg_smoother.set_steps(2); 
 * 
 * @endcode
 * 
 * 由于SOR方法不是对称的，但我们在下面使用共轭梯度迭代，这里有一个技巧，使多级预处理器成为对称算子，即使是对非对称平滑器。
 * 

 * 
 * 
 * @code
 *     mg_smoother.set_symmetric(true); 
 * 
 * @endcode
 * 
 * 平滑器类可以选择实现变量V型循环，我们在这里不需要。
 * 

 * 
 * 
 * @code
 *     mg_smoother.set_variable(false); 
 * 
 * @endcode
 * 
 * 最后，我们必须将我们的矩阵包裹在一个具有所需乘法函数的对象中。
 * 

 * 
 * 
 * @code
 *     mg::Matrix<Vector<double>> mgmatrix(mg_matrix); 
 *     mg::Matrix<Vector<double>> mgdown(mg_matrix_dg_down); 
 *     mg::Matrix<Vector<double>> mgup(mg_matrix_dg_up); 
 * 
 * @endcode
 * 
 * 现在，我们准备设置V型循环算子和多级预处理程序。
 * 

 * 
 * 
 * @code
 *     Multigrid<Vector<double>> mg( 
 *       mgmatrix, mg_coarse, mg_transfer, mg_smoother, mg_smoother); 
 * 
 * @endcode
 * 
 * 让我们不要忘记因为自适应细化而需要的边缘矩阵。
 * 

 * 
 * 
 * @code
 *     mg.set_edge_flux_matrices(mgdown, mgup); 
 * 
 * @endcode
 * 
 * 在所有的准备工作完成后，将Multigrid对象包装成另一个对象，它可以作为一个普通的预处理程序使用。
 * 

 * 
 * 
 * @code
 *     PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
 *       preconditioner(dof_handler, mg, mg_transfer); 
 * 
 * @endcode
 * 
 * 并用它来解决这个系统。
 * 

 * 
 * 
 * @code
 *     solver.solve(matrix, solution, right_hand_side, preconditioner); 
 *   } 
 * 
 * @endcode
 * 
 * 另一个克隆的集合函数。与之前的最大区别是，这里我们也有一个输入向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double InteriorPenaltyProblem<dim>::estimate() 
 *   { 
 * 
 * @endcode
 * 
 * 估算器的结果存储在一个每个单元格有一个条目的向量中。由于deal.II中的单元格没有编号，我们必须建立自己的编号，以便使用这个向量。对于下面使用的汇编器来说，结果存储在向量的哪个分量中的信息是由每个单元的user_index变量传送的。我们需要在这里设置这个编号。
 * 

 * 
 * 另一方面，有人可能已经使用了用户指数。所以，让我们做个好公民，在篡改它们之前保存它们。
 * 

 * 
 * 
 * @code
 *     std::vector<unsigned int> old_user_indices; 
 *     triangulation.save_user_indices(old_user_indices); 
 * 
 *     estimates.block(0).reinit(triangulation.n_active_cells()); 
 *     unsigned int i = 0; 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       cell->set_user_index(i++); 
 * 
 * @endcode
 * 
 * 这就像以前一样开始。
 * 

 * 
 * 
 * @code
 *     MeshWorker::IntegrationInfoBox<dim> info_box; 
 *     const unsigned int                  n_gauss_points = 
 *       dof_handler.get_fe().tensor_degree() + 1; 
 *     info_box.initialize_gauss_quadrature(n_gauss_points, 
 *                                          n_gauss_points + 1, 
 *                                          n_gauss_points); 
 * 
 * @endcode
 * 
 * 但现在我们需要通知信息框我们要在正交点上评估的有限元函数。首先，我们用这个向量创建一个AnyData对象，这个向量就是我们刚刚计算的解。
 * 

 * 
 * 
 * @code
 *     AnyData solution_data; 
 *     solution_data.add<const Vector<double> *>(&solution, "solution"); 
 * 
 * @endcode
 * 
 * 然后，我们告诉单元格的 Meshworker::VectorSelector ，我们需要这个解决方案的二次导数（用来计算拉普拉斯）。因此，选择函数值和第一导数的布尔参数是假的，只有选择第二导数的最后一个参数是真的。
 * 

 * 
 * 
 * @code
 *     info_box.cell_selector.add("solution", false, false, true); 
 * 
 * @endcode
 * 
 * 在内部和边界面，我们需要函数值和第一导数，但不需要第二导数。
 * 

 * 
 * 
 * @code
 *     info_box.boundary_selector.add("solution", true, true, false); 
 *     info_box.face_selector.add("solution", true, true, false); 
 * 
 * @endcode
 * 
 * 我们继续像以前一样，除了默认的更新标志已经被调整为我们上面要求的值和导数之外。
 * 

 * 
 * 
 * @code
 *     info_box.add_update_flags_boundary(update_quadrature_points); 
 *     info_box.initialize(fe, mapping, solution_data, solution); 
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler); 
 * 
 * @endcode
 * 
 * 汇编器在每个单元格中存储一个数字，否则这与右侧的计算是一样的。
 * 

 * 
 * 
 * @code
 *     MeshWorker::Assembler::CellsAndFaces<double> assembler; 
 *     AnyData                                      out_data; 
 *     out_data.add<BlockVector<double> *>(&estimates, "cells"); 
 *     assembler.initialize(out_data, false); 
 * 
 *     Estimator<dim> integrator; 
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
 *                                            dof_handler.end(), 
 *                                            dof_info, 
 *                                            info_box, 
 *                                            integrator, 
 *                                            assembler); 
 * 
 * @endcode
 * 
 * 就在我们返回错误估计的结果之前，我们恢复旧的用户索引。
 * 

 * 
 * 
 * @code
 *     triangulation.load_user_indices(old_user_indices); 
 *     return estimates.block(0).l2_norm(); 
 *   } 
 * 
 * @endcode
 * 
 * 这里我们把我们的有限元解和（已知的）精确解进行比较，计算梯度和函数本身的平均二次误差。这个函数是上面那个估计函数的克隆。
 * 

 * 
 * 由于我们分别计算能量和<i>L<sup>2</sup></i>-norm的误差，我们的块向量在这里需要两个块。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::error() 
 *   { 
 *     BlockVector<double> errors(2); 
 *     errors.block(0).reinit(triangulation.n_active_cells()); 
 *     errors.block(1).reinit(triangulation.n_active_cells()); 
 * 
 *     std::vector<unsigned int> old_user_indices; 
 *     triangulation.save_user_indices(old_user_indices); 
 *     unsigned int i = 0; 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       cell->set_user_index(i++); 
 * 
 *     MeshWorker::IntegrationInfoBox<dim> info_box; 
 *     const unsigned int                  n_gauss_points = 
 *       dof_handler.get_fe().tensor_degree() + 1; 
 *     info_box.initialize_gauss_quadrature(n_gauss_points, 
 *                                          n_gauss_points + 1, 
 *                                          n_gauss_points); 
 * 
 *     AnyData solution_data; 
 *     solution_data.add<Vector<double> *>(&solution, "solution"); 
 * 
 *     info_box.cell_selector.add("solution", true, true, false); 
 *     info_box.boundary_selector.add("solution", true, false, false); 
 *     info_box.face_selector.add("solution", true, false, false); 
 * 
 *     info_box.add_update_flags_cell(update_quadrature_points); 
 *     info_box.add_update_flags_boundary(update_quadrature_points); 
 *     info_box.initialize(fe, mapping, solution_data, solution); 
 * 
 *     MeshWorker::DoFInfo<dim> dof_info(dof_handler); 
 * 
 *     MeshWorker::Assembler::CellsAndFaces<double> assembler; 
 *     AnyData                                      out_data; 
 *     out_data.add<BlockVector<double> *>(&errors, "cells"); 
 *     assembler.initialize(out_data, false); 
 * 
 *     ErrorIntegrator<dim> integrator; 
 *     MeshWorker::integration_loop<dim, dim>(dof_handler.begin_active(), 
 *                                            dof_handler.end(), 
 *                                            dof_info, 
 *                                            info_box, 
 *                                            integrator, 
 *                                            assembler); 
 *     triangulation.load_user_indices(old_user_indices); 
 * 
 *     deallog << "energy-error: " << errors.block(0).l2_norm() << std::endl; 
 *     deallog << "L2-error:     " << errors.block(1).l2_norm() << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 创建图形输出。我们通过整理其各个组成部分的名称来产生文件名，包括我们用两个数字输出的细化周期。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   InteriorPenaltyProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     const std::string filename = 
 *       "sol-" + Utilities::int_to_string(cycle, 2) + ".gnuplot"; 
 * 
 *     deallog << "Writing solution to <" << filename << ">..." << std::endl 
 *             << std::endl; 
 *     std::ofstream gnuplot_output(filename); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "u"); 
 *     data_out.add_data_vector(estimates.block(0), "est"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     data_out.write_gnuplot(gnuplot_output); 
 *   } 
 * 
 * @endcode
 * 
 * 最后是自适应循环，或多或少和前面的例子一样。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void InteriorPenaltyProblem<dim>::run(unsigned int n_steps) 
 *   { 
 *     deallog << "Element: " << fe.get_name() << std::endl; 
 *     for (unsigned int s = 0; s < n_steps; ++s) 
 *       { 
 *         deallog << "Step " << s << std::endl; 
 *         if (estimates.block(0).size() == 0) 
 *           triangulation.refine_global(1); 
 *         else 
 *           { 
 *             GridRefinement::refine_and_coarsen_fixed_fraction( 
 *               triangulation, estimates.block(0), 0.5, 0.0); 
 *             triangulation.execute_coarsening_and_refinement(); 
 *           } 
 * 
 *         deallog << "Triangulation " << triangulation.n_active_cells() 
 *                 << " cells, " << triangulation.n_levels() << " levels" 
 *                 << std::endl; 
 * 
 *         setup_system(); 
 *         deallog << "DoFHandler " << dof_handler.n_dofs() << " dofs, level dofs"; 
 *         for (unsigned int l = 0; l < triangulation.n_levels(); ++l) 
 *           deallog << ' ' << dof_handler.n_dofs(l); 
 *         deallog << std::endl; 
 * 
 *         deallog << "Assemble matrix" << std::endl; 
 *         assemble_matrix(); 
 *         deallog << "Assemble multilevel matrix" << std::endl; 
 *         assemble_mg_matrix(); 
 *         deallog << "Assemble right hand side" << std::endl; 
 *         assemble_right_hand_side(); 
 *         deallog << "Solve" << std::endl; 
 *         solve(); 
 *         error(); 
 *         deallog << "Estimate " << estimate() << std::endl; 
 *         output_results(s); 
 *       } 
 *   } 
 * } // namespace Step39 
 * 
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step39; 
 * 
 *       deallog.depth_console(2); 
 *       std::ofstream logfile("deallog"); 
 *       deallog.attach(logfile); 
 *       FE_DGQ<2>                 fe1(3); 
 *       InteriorPenaltyProblem<2> test1(fe1); 
 *       test1.run(12); 
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
 * @endcode
examples/step-39/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Logfileoutput"></a><h3>Logfile output</h3> 首先，该程序产生通常的日志文件，在这里存储在<tt>deallog</tt>。它的内容是（省略了中间的步骤


@code
DEAL::Element: FE_DGQ<2>(3)
DEAL::Step 0
DEAL::Triangulation 16 cells, 2 levels
DEAL::DoFHandler 256 dofs, level dofs 64 256
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 13 value 1.64974e-13
DEAL::energy-error: 0.297419
DEAL::L2-error:     0.00452447
DEAL::Estimate 0.990460
DEAL::Writing solution to <sol-00.gnuplot>...
DEAL::
DEAL::Step 1
DEAL::Triangulation 25 cells, 3 levels
DEAL::DoFHandler 400 dofs, level dofs 64 256 192
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 14 value 3.72262e-13
DEAL::energy-error: 0.258559
DEAL::L2-error:     0.00288510
DEAL::Estimate 0.738624
DEAL::Writing solution to <sol-01.gnuplot>...
DEAL::
DEAL::Step 2
DEAL::Triangulation 34 cells, 4 levels
DEAL::DoFHandler 544 dofs, level dofs 64 256 256 128
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 37.4071
DEAL:cg::Convergence step 15 value 1.91610e-13
DEAL::energy-error: 0.189234
DEAL::L2-error:     0.00147954
DEAL::Estimate 0.657507
DEAL::Writing solution to <sol-02.gnuplot>...


...


DEAL::Step 10
DEAL::Triangulation 232 cells, 11 levels
DEAL::DoFHandler 3712 dofs, level dofs 64 256 896 768 768 640 512 256 256 256 256
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 51.1571
DEAL:cg::Convergence step 15 value 7.19599e-13
DEAL::energy-error: 0.0132475
DEAL::L2-error:     1.00423e-05
DEAL::Estimate 0.0470724
DEAL::Writing solution to <sol-10.gnuplot>...
DEAL::
DEAL::Step 11
DEAL::Triangulation 322 cells, 12 levels
DEAL::DoFHandler 5152 dofs, level dofs 64 256 1024 1024 896 768 768 640 448 320 320 320
DEAL::Assemble matrix
DEAL::Assemble multilevel matrix
DEAL::Assemble right hand side
DEAL::Solve
DEAL:cg::Starting value 52.2226
DEAL:cg::Convergence step 15 value 8.15195e-13
DEAL::energy-error: 0.00934891
DEAL::L2-error:     5.41095e-06
DEAL::Estimate 0.0329102
DEAL::Writing solution to <sol-11.gnuplot>...
DEAL::
@endcode



例如，该日志显示共轭梯度迭代步骤的数量恒定在大约15个。

<a name="Postprocessingofthelogfile"></a><h3>Postprocessing of the logfile</h3>


 <img src="https://www.dealii.org/images/steps/developer/step-39-convergence.svg" alt="">  使用perl脚本<tt>postprocess.pl</tt>，我们提取相关数据到<tt>output.dat</tt>，可以用<tt>gnuplot</tt>绘制图形。例如，上面的图是用gnuplot脚本<tt>plot_errors.gpl</tt>制作的，通过

@code
perl postprocess.pl deallog &> output.dat
gnuplot plot_errors.gpl
@endcode



参考数据可以在<tt>output.reference.dat</tt>中找到。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-39.cc"
*/
