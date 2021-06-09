//include/deal.II-translator/A-headers/vector_valued_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



/**
 *    @defgroup vector_valued Handling vector valued problems
 * 矢量值问题是偏微分方程的系统。这些问题的解变量不是一个标量函数，而是一个矢量值函数或一组函数。这包括，例如。   <ul>   <li>   step-8 、 step-17 和 step-18 中讨论的弹性方程，其解是每一点的向量值位移。     <li>  在 step-20 和 step-21 中讨论的混合拉普拉斯方程及其扩展，其解是每一点的标量压力和矢量速度。     <li>  在 step-22 中讨论的斯托克斯方程及其扩展，以及 step-31 ，其中的解也是每一点的标量压力和矢量速度。     <li>  由实部和虚部组成的复值解，例如在  step-29  中讨论。   </ul>
 * 本页概述了如何在deal.II中轻松实现此类矢量值问题。特别是，它解释了FESystem类的用法，它允许我们为偏微分系统编写代码，就像我们为单个方程编写代码一样。
 * @dealiiVideoLecture{19,20,21} <table class="tutorial" width="50%">
 * <tr><th><b>%Table of contents</b></th></tr> <tr><td width="100%"
 * valign="top">
 * <ol>
 * <li> @ref VVExamples "Examples of vector-valued problems"
 * <li> @ref VVFEs "Describing finite element spaces"
 * <li> @ref VVAssembling "Assembling linear systems"
 * <li> @ref VVAlternative "An alternative approach"
 * <li> @ref VVBlockSolvers "Block solvers"
 * <li> @ref VVExtracting "Extracting data from solutions"
 * <li> @ref VVOutput "Generating graphical output"
 * </ol> </td> </tr> </table>
 *
 *
 * @anchor  VVExamples <h3>Examples of vector-valued problems</h3>。
 * 系统地处理向量值问题的方式与标量问题没有根本的不同：首先，我们需要一个考虑到所有解变量的问题的弱（变）式。在我们这样做之后，生成系统矩阵和求解线性系统遵循我们已经习惯的大纲。
 * <h4>Linear elasticity</h4> 让我们以 step-8
 * 中的弹性问题为例，甚至通过选择 $\lambda = 0$ 和 $\mu = 1$
 * 来简化它以突出重要的概念。因此，让我们考虑下面的弱表述：找到
 * $\mathbf u \in \mathbf V = H^1_0(\Omega; \mathbb R^3)$ ，使所有
 * $\mathbf
 * v\in V$ 都持有@f[
 * a(u,v) \equiv 2\int_{\Omega} \mathbf D\mathbf u : \mathbf D\mathbf
 * v\,dx = \int_\Omega \mathbf f\cdot \mathbf v \,dx.
 * @f] 这里，<b>D</b>表示由 $\mathbf Du = \tfrac12 (\nabla \mathbf u + (\nabla \mathbf u)^T)$ 定义的对称梯度，冒号表示两个等级为2的张量的双重收缩（Frobenius内部积）。这种双线性形式看起来确实非常像  step-3  中泊松问题的双线性形式。唯一的区别是  <ol>   <li>  我们用对称梯度替换了梯度算子；这实际上不是一个重大的区别，如果你用  $\nabla$  替换  $\mathbf D$  ，这里说的一切都是真的。事实上，让我们这样做来简化讨论。@f[
 * a(u,v) \equiv \int_{\Omega} \nabla\mathbf u : \nabla\mathbf v\,dx =
 * \int_\Omega \mathbf f\cdot \mathbf v \,dx.
 * @f]但请注意，这个系统并不十分令人兴奋，因为我们可以分别解决<b>u</b>的三个组成部分。
 * <li>  现在的试验和测试函数来自空间 $H^1_0(\Omega; \mathbb R^3)$ ，它可以被视为标量空间 $H^1_0(\Omega)$ 的三份副本。而这正是我们将在下面使用FESystem来实现这个空间的方法。   </ol>
 * 但现在，让我们再仔细看看这个系统。首先，让我们利用<b>u</b>=(<i>u</i><sub>1</sub>, <i>u</i><sub>2</sub>, <i>u</i><sub>3</sub>)<sup>T</sup>和<b>v</b>相应。然后，我们可以把坐标中的简化方程写成@f[
 * a(u,v) = \int_\Omega \bigl(\nabla u_1\cdot \nabla v_1
 * +\nabla u_2\cdot \nabla v_2+\nabla u_3\cdot \nabla v_3\bigr)\,dx
 * = \int_\Omega \bigl(f_1v_1 + f_2 v_2 + f_3 v_3\bigr)\,dx.
 * @f]。 我们看到，这只是拉普拉斯的双线性形式的三份拷贝，一份应用于每个分量（这是与 $\mathbf D$ 的表述更令人兴奋的地方，我们想推导出一个也适用于该表述的框架）。我们可以通过选择特殊的测试函数使这个弱形式再次成为微分方程组：首先，选择<b>v</b>=（<i>v</i><sub>1</sub>,0,0）<sup>T</sup>，然后<b>v</b>=（0，<i>v</i><sub>2</sub>,0)<sup>T</sup>, 最后<b>v</b>=(0,0,<i>v</i><sub>3</sub>)<sup>T</sup>. 将这些结果写在彼此下面，我们得到系统@f[
 * \begin{matrix} (\nabla u_1,\nabla v_1) &&& = (f_1, v_1) \\ & (\nabla
 * u_2,\nabla v_2) && = (f_2, v_2) \\ && (\nabla u_3,\nabla v_3) & = (f_3,
 * v_3) \end{matrix} @f]，我们使用标准内积符号  $(\mathbf f,\mathbf
 * g) =
 * \int_\Omega \mathbf f \cdot \mathbf g \,dx$  。对我们的理解很重要的是，我们要记住，后一种形式作为PDE系统完全等同于双线性形式的原始定义<i>a</i>(<i>u</i>,<i>v</i>)，它并不立即表现出这种系统结构。最后，让我们写出具有对称梯度的弹性方程<b>D</b>的完整系统。@f[
 * \begin{matrix}
 * (\nabla u_1,\nabla v_1) + (\partial_1 u_1,\partial_1 v_1)
 * & (\partial_1 u_2,\partial_2 v_1)
 * & (\partial_1 u_3,\partial_3 v_1)
 * & = (f_1, v_1)
 * \\
 * (\partial_2 u_1,\partial_1 v_2)
 * & (\nabla u_2,\nabla v_2) + (\partial_2 u_2,\partial_2 v_2)
 * & (\partial_2 u_3,\partial_3 v_2)
 * & = (f_2, v_2)
 * \\
 * (\partial_3 u_1,\partial_1 v_3)
 * & (\partial_3 u_2,\partial_2 v_3)
 * & (\nabla u_3,\nabla v_3) + (\partial_3 u_3,\partial_3 v_3)
 * & = (f_3, v_3)
 * \end{matrix}.
 * @f] 非常正式地，如果我们相信算子值矩阵，我们可以将其改写为<b>v</b><sup>T</sup><b>Au</b> = <b>v</b><sup>T</sup><b>f</b> 或 @f[
 * \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}^T \begin{pmatrix} (\nabla
 * \cdot,\nabla \cdot) + (\partial_1 \cdot,\partial_1 \cdot) & (\partial_1
 * \cdot,\partial_2 \cdot) & (\partial_1 \cdot,\partial_3 \cdot) \\
 * (\partial_2 \cdot,\partial_1 \cdot) & (\nabla \cdot,\nabla \cdot) +
 * (\partial_2 \cdot,\partial_2 \cdot) & (\partial_2 \cdot,\partial_3 \cdot)
 * \\ (\partial_3 \cdot,\partial_1 \cdot) & (\partial_3 \cdot,\partial_2
 * \cdot) & (\nabla \cdot,\nabla \cdot) + (\partial_3 \cdot,\partial_3 \cdot)
 * \end{pmatrix} \begin{pmatrix} u_1 \\ u_2 \\ u_3 \end{pmatrix} =
 * \begin{pmatrix} v_1 \\ v_2 \\ v_3 \end{pmatrix}^T \begin{pmatrix} f_1 \\
 * f_2 \\ f_3\end{pmatrix} @f] 的形式。
 * *<h4>Mixed elliptic problems</h4> 现在，让我们考虑一个更复杂的例子，即 step-20 中讨论的三维混合拉普拉斯方程：@f{eqnarray*}
 * \textbf{u} + \nabla p &=& 0,
 * \\
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -\textrm{div}\; \textbf{u} &=& f,
 * @f}
 *
 * 这里，我们有四个解分量：标量压力 $p \in L^2(\Omega)$
 * 和有三个矢量分量的矢量值速度 $\mathbf u \in \mathbf V =
 * H^{\text{div}}_0(\Omega)$
 * 。请注意，作为与前一个例子的重要区别，向量空间<b>V</b>并不仅仅是三个相同空间的简单复制/。
 * 对于这个问题和其他矢量问题，获得弱式或变式的系统方法是，首先将其视为一个问题，其中运算符和解变量以矢量和矩阵形式书写。对于这个例子，这将读作：@f{eqnarray*}
 * \left(
 * \begin{array}{cc} \mathbf 1 & \nabla \\
 *
 * -\nabla^T & 0 \end{array}
 * \right)
 * \left(
 * \begin{array}{c} \mathbf u \\ p \end{array}
 * \right)
 * =
 * \left(
 * \begin{array}{c} \mathbf 0 \\ f \end{array}
 * \right)
 * @f} 。
 *
 * 这就清楚地表明，解@f{eqnarray*}
 * U =
 * \left(
 * \begin{array}{c} \mathbf u \\ p \end{array}
 * \right)
 * @f}有四个分量。
 * 确实有四个成分。我们注意到，如果我们同时改变矩阵运算符的列，我们可以改变解的成分
 * $\textbf u$ 和 $p$ 在 $U$ 里面的顺序。
 * 接下来，我们需要考虑测试函数  $V$  。我们要将方程的两边都与它们相乘，然后对 $\Omega$ 进行积分。结果应该是一个标量的相等。我们可以通过选择 $V$ 来实现这一点，因为@f{eqnarray*}
 * V =
 * \left(
 * \begin{array}{c} \mathbf v \\ q \end{array}
 * \right).
 * @f}也是矢量值的。
 *
 * 将矩阵-向量方程与测试函数从左边相乘是很方便的，因为这样我们以后会自动得到正确的矩阵（在线性系统中，矩阵也是从右边与解变量相乘的，而不是从左边），而如果我们从右边相乘，那么这样汇集的矩阵就是我们真正想要的矩阵的转置。
 * 考虑到这一点，让我们乘以 $V$ 并进行积分，得到以下方程，该方程对所有测试函数 $V$ 必须成立：@f{eqnarray*}
 * \int_\Omega
 * \left(
 * \begin{array}{c} \mathbf v \\ q \end{array}
 * \right)^T
 * \left(
 * \begin{array}{cc} \mathbf 1 & \nabla \\
 *
 * -\nabla^T & 0 \end{array}
 * \right)
 * \left(
 * \begin{array}{c} \mathbf u \\ p \end{array}
 * \right)
 * \ dx
 * =
 * \int_\Omega
 * \left(
 * \begin{array}{c} \mathbf v \\ q \end{array}
 * \right)^T
 * \left(
 * \begin{array}{c} \mathbf 0 \\ f \end{array}
 * \right)
 * \ dx,
 * @f}
 * 或者等同于：@f{eqnarray*}
 * (\mathbf v, \mathbf u)
 * +
 * (\mathbf v, \nabla p)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 * (q, \mathrm{div}\ \mathbf u)
 * =
 * (q,f),
 * @f} 。
 *
 *
 * 我们通过对第二项的部分积分得到最终形式：@f{eqnarray*}
 * (\mathbf v, \mathbf u)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 * (\mathrm{div}\ \mathbf v, p)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 * (q, \mathrm{div}\ \mathbf u)
 * =
 * (q,f)
 *
 * - (\mathbf n\cdot\mathbf v, p)_{\partial\Omega}.
 * @f} 。
 *
 * 正是这种形式，我们以后将用于将离散的弱形式组装成一个矩阵和一个右手边的向量：在这种形式中，我们有解和测试函数
 * $U,V$ ，每个函数都由一些向量分量组成，我们可以提取。
 *
 *  @anchor  VVFEs<h3>Describing finite element spaces</h3> 。
 * 一旦我们确定了双线性形式和函数设置，我们就需要找到一种方法来描述矢量值有限元空间，并从中提取解和检验函数。这就是FESystem类的作用：它由较简单的空间组成矢量值的有限元空间。在弹性问题的例子中，我们需要同一元素的
 * <code>dim</code> 份，例如
 * @code
 * FESystem<dim> elasticity_element (FE_Q<dim>(1), dim);
 * @endcode
 * 这将产生一个维度为 <code>dim</code>
 * 的向量估值空间，其中每个分量都是FE_Q类型的连续双线性元素。它将有
 * <code>dim</code>
 * 倍于相应的FE_Q的基函数，这些基函数中的每一个都是FE_Q的基函数，被提升到矢量的一个分量中。
 * 对于混合拉普拉斯，情况更为复杂。首先，我们必须确定一对离散空间
 * $\mathbf V_h \times Q_h \subset H^{\text{div}}_0(\Omega) \times
 * L^2_0(\Omega)$  。一种选择是稳定的Raviart-Thomas对
 * @code
 * FESystem<dim> rt_element (FE_RaviartThomas<dim>(1), 1,
 *                           FE_DGQ<dim>(1),          1);
 * @endcode
 * 这个系统中的第一个元素已经是一个维度为 <code>dim</code>
 * 的矢量值元素，而第二个元素是一个普通的标量元素。
 * 除了使用稳定的Raviart-Thomas对之外，我们还可以考虑混合拉普拉斯的稳定公式，例如LDG方法。在这里，我们可以选择使用相同的空间来计算速度分量和压力，即
 * @code
 * FESystem<dim> ldg_convoluted_element_1 (FE_DGQ<dim>(1), dim+1);
 * @endcode
 * 这个系统只是有 <code>dim+1</code>
 * 个相同的不连续元素的相等拷贝，这并没有真正反映系统的结构。因此，我们倾向于
 * @code
 * FESystem<dim> ldg_equal_element (FESystem<dim>(FE_DGQ<dim>(1), dim), 1,
 *                                  FE_DGQ<dim>(1),                     1);
 * @endcode
 * 这里，我们有一个由两个元素组成的系统，一个是矢量值，一个是标量值，很像与
 * <code>rt_element</code>
 * 。事实上，在许多代码中，这两个可以互换。这个元素也允许我们很容易地切换到速度的低阶近似的LDG方法，即
 * @code
 * FESystem<dim> ldg_unequal_element (FESystem<dim>(FE_DGQ<dim>(1), dim), 1,
 *                                    FE_DGQ<dim>(2),                     1);
 * @endcode
 * 必须指出的是，这个元素不同于
 * @code
 * FESystem<dim> ldg_convoluted_element_2 (FE_DGQ<dim>(1), dim,
 *                                         FE_DGQ<dim>(2), 1);
 * @endcode
 * 虽然构造函数的调用与 <code>rt_element</code>
 * 非常相似，但结果实际上更像
 * <code>ldg_convoluted_element_1</code> ，因为这个元素产生
 * <code>dim+1</code>
 * 的独立组件。下面是对产生的FESystem对象的更详细的比较。
 * <h4>Internal structure of FESystem</h4>
 * FESystem有一些内部变量，反映了构造函数所设置的内部结构。然后这些也可以被应用程序用来给矩阵组装和线性代数提供结构。我们在下表中给出了上述例子中这些变量的名称和值。<table
 * border="1"> <tr><th>系统元素</th>
 * <th>FiniteElementData::n_blocks()</th>
 * <th>FiniteElementData::n_components()</th>
 * <th>FiniteElement::n_base_elements()</th>  </tr> <tr><td>
 * <code>elasticity_element</code></td><td><code>dim</code></td><td><code>dim</code>
 * </td><td>1</td> </tr> <tr><td>
 * <code>rt_element</code></td><td>2</td><td><code>dim+1</code>
 * </td><td>2</td> </tr> <tr><td>
 * <code>ldg_equal_element</code></td><td>2</td><td><code>dim+1</code>
 * </td><td>2</td> </tr> <tr><td>
 * <code>ldg_convoluted_element_1</code></td><td><code>dim+1</code></td><td><code>dim+1</code>
 * </td><td>1</td> </tr> <tr><td>
 * <code>ldg_convoluted_element_2</code></td><td><code>dim+1</code></td><td><code>dim+1</code>
 * </td><td>2</tr> </table> 从这个表中可以看出，FES系统反映了
 * <code>rt_element</code> 和 <code>ldg_equal_element</code>
 * 情况下微分方程组的很多结构，因为我们有一个矢量值和一个标量变量。另一方面，卷积元素没有这种结构，我们必须在组装系统时以某种方式重构它，如下所述。
 * 在这一点上，需要注意的是，两个FES系统对象的嵌套可以给整个FES系统带来更丰富的结构，而不仅仅是将它们串联起来。这种结构可以被应用程序所利用，但不是自动的。
 * @anchor  VVAssembling <h3>Assembling linear systems</h3>
 * 下一步是对线性系统进行组装。对于标量问题的简单情况，如何做到这一点已经在许多教程程序中显示出来，首先是
 * step-3
 * 。在这里，我们将展示如何对矢量问题进行处理。对应于上述弱式的不同特征和创建的不同系统元素，我们有几个选择，概述如下。
 * 整个概念可能最好的解释是通过展示一个例子，说明如何组装一个单元对上述混合拉普拉斯方程的弱形式的局部贡献。
 * <h4>A single FEValues and FEValuesExtractors</h4> 这基本上是 step-20
 * 的做法。
 * @code
 * const FEValuesExtractors::Vector velocities (0);
 * const FEValuesExtractors::Scalar pressure (dim);
 *
 * ...
 *
 * typename DoFHandler<dim>::active_cell_iterator
 *  cell = dof_handler.begin_active(),
 *  endc = dof_handler.end();
 * for (; cell!=endc; ++cell)
 *  {
 *    fe_values.reinit (cell);
 *    local_matrix = 0;
 *    local_rhs = 0;
 *
 *    right_hand_side.value_list (fe_values.get_quadrature_points(),
 *                                rhs_values);
 *
 *    for (unsigned int q=0; q<n_q_points; ++q)
 *      for (unsigned int i=0; i<dofs_per_cell; ++i)
 *        {
 *          for (unsigned int j=0; j<dofs_per_cell; ++j)
 *            local_matrix(i,j) += (fe_values[velocities].value (i, q)
 *                                  fe_values[velocities].value (j, q)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *                                  fe_values[velocities].divergence (i, q)
 *                                  fe_values[pressure].value (j, q)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *                                  fe_values[pressure].value (i, q)
 *                                  fe_values[velocities].divergence (j, q))
 *                                  fe_values.JxW(q);
 *
 *          local_rhs(i) +=
 *
 * - fe_values[pressure].value (i, q)
 *                          rhs_values[q]
 *                          fe_values.JxW(q);
 *        }
 * @endcode
 *
 * 所以这里是发生了什么。   <ul>   <li>  我们做的第一件事是声明 "抽取器"（见FEValuesExtractors命名空间）。这些对象除了存储矢量值有限元的哪些分量构成单一的标量分量或秩1的张量（即我们所说的 "物理矢量"，总是由 <code>dim</code> 分量组成）外，没有什么作用。在这里，我们声明一个对象，表示由 <code>dim</code> 分量组成的速度，从零分量开始，以及压力的提取器，它是位置 <code>dim</code> 的标量分量。
 * <li>
 * 然后我们对所有单元、形状函数和正交点进行常规循环。在最内部的循环中，我们计算一对形状函数对全局矩阵和右手向量的局部贡献。回顾一下，根据形状函数
 * $V_i=\left(\begin{array}{c}\mathbf v_i \\ q_i\end{array}\right),
 *       V_j=\left(\begin{array}{c}\mathbf v_j \\ q_j\end{array}\right)$ ：@f{eqnarray*}
 *          (\mathbf v_i, \mathbf v_j)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *          (\mathrm{div}\ \mathbf v_i, q_j)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *          (q_i, \mathrm{div}\ \mathbf v_j)
 *        @f}，单元格对双线性形式的贡献（即忽略边界条款）看起来如下
 * 而实施起来则是这样的。
 * @code
 *            local_matrix(i,j) += (fe_values[velocities].value (i, q)
 *                                  fe_values[velocities].value (j, q)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *                                  fe_values[velocities].divergence (i, q)
 *                                  fe_values[pressure].value (j, q)
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -
 *                                  fe_values[pressure].value (i, q)
 *                                  fe_values[velocities].divergence (j, q)
 *                                 )
 *                                 fe_values.JxW(q);
 * @endcode
 * 相似之处是相当明显的。
 * <li>  基本上，上述代码中发生的情况是这样的：当你执行
 * <code>fe_values[pressure]</code>  时，会创建一个所谓的
 * "视图"，即一个与完整的FEValues对象不同的对象，它不代表有限元的所有组件，而只代表提取器对象
 * <code>pressure</code>  或  <code>velocities</code>
 * 所代表的那（些）组件。
 * <li>
 * 然后可以向这些视图询问关于这些单独组件的信息。例如，当你写
 * <code>fe_values[pressure].value(i,q)</code> 时，你会得到 $i$
 * 第1个形状函数 $V_i$ 在 $q$
 * 第1个正交点的压力分量的值。因为提取器
 * <code>pressure</code> 代表一个标量分量，运算器
 * <code>fe_values[pressure].value(i,q)</code>
 * 的结果是一个标量数。另一方面，调用
 * <code>fe_values[velocities].value(i,q)</code> 将产生整组
 * <code>dim</code> 分量的值，其类型为 <code>Tensor@<1,dim@></code>
 * 。
 * <li>
 * 其他可以用视图做的事情是要求提取器所描述的特定形状函数的分量的梯度。例如，
 * <code>fe_values[pressure].gradient(i,q)</code>
 * 表示标量压力分量的梯度，其类型为
 * <code>Tensor@<1,dim@></code> ，而速度分量的梯度，
 * <code>fe_values[velocities].gradient(i,q)</code> 为
 * <code>Tensor@<2,dim@></code> ，即一个由条目
 * $G_{ij}=\frac{\partial\phi_i}{\partial x_j}$ 组成的矩阵 $G_{ij}$
 * 。最后，标量和矢量视图都可以询问二阶导数（"Hessians"），矢量视图可以询问对称梯度，定义为
 * $S_{ij}=\frac 12 \left[\frac{\partial\phi_i}{\partial x_j}
 * + \frac{\partial\phi_j}{\partial x_i}\right]$  以及分歧  $\sum_{d=0}^{dim-1} \frac{\partial\phi_d}{\partial x_d}$  。   </ul>  其他使用提取器和视图的例子见教程程序  step-21  ,  step-22  ,  step-31  和其他几个程序。
 *
 * @note
 * 在目前的背景下，当我们谈论一个矢量时（例如在提取上面的速度分量时），我们指的是物理学上使用的这个词：它有
 * <code>spacedim</code>
 * 分量，在坐标系变换下以特定方式表现出来。例子包括速度场或位移场。这与数学中使用
 * "矢量
 * "一词的方式相反（以及我们在库中的其他上下文中使用这个词的方式，例如在矢量类中），在那里它真正代表了一个数字的集合。后者的一个例子是火焰中化学物种浓度的集合；然而，这些实际上只是标量变量的集合，因为如果坐标系被旋转，它们不会改变，不像速度矢量的分量，因此，这个
 * FEValuesExtractors::Vector 类不应该被用于这种情况。
 *
 *  @anchor  VVAlternative <h3>An alternative approach</h3>。
 * 在有些情况下，我们可以利用所使用的有限元的知识，对矩阵或右手边向量的装配进行一些优化。例如，考虑我们在
 * step-8  中首先关注的弹性方程的双线性形式。
 * @f[
 * a({\mathbf u}, {\mathbf v}) =
 * \left(
 *  \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
 * \right)_\Omega
 * +
 * \sum_{i,j}
 * \left(
 *  \mu \partial_i u_j, \partial_i v_j
 * \right)_\Omega,
 * +
 * \sum_{i,j}
 * \left(
 *  \mu \partial_i u_j, \partial_j v_i
 * \right)_\Omega,
 * @f]
 * 这里， $\mathbf u$ 是一个具有 <code>dim</code>
 * 分量的向量函数， $\mathbf v$ 是相应的测试函数，
 * $\lambda,\mu$
 * 是材料参数。鉴于我们上面的讨论，实现这种双线性形式的明显方法如下，使用一个提取器对象，将有限元的所有
 * <code>dim</code>
 * 分量解释为单一矢量，而不是不相交的标量分量。
 *
 * @code
 *    const FEValuesExtractors::Vector displacements (0);
 *
 *    ...
 *
 *    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
 *      for (unsigned int i=0; i<dofs_per_cell; ++i)
 *        {
 *          const Tensor<2,dim> phi_i_grad
 *            = fe_values[displacements].gradient (i,q_point);
 *          const double phi_i_div
 *            = fe_values[displacements].divergence (i,q_point);
 *
 *          for (unsigned int j=0; j<dofs_per_cell; ++j)
 *            {
 *              const Tensor<2,dim> phi_j_grad
 *                = fe_values[displacements].gradient (j,q_point);
 *              const double phi_j_div
 *                = fe_values[displacements].divergence (j,q_point);
 *
 *              cell_matrix(i,j)
 *                +=  (lambda_values[q_point]
 *                     phi_i_div phi_j_div
 *                     +
 *                     mu_values[q_point]
 *                     double_contract(phi_i_grad, phi_j_grad)
 *                     +
 *                     mu_values[q_point]
 *                     double_contract(phi_i_grad, transpose(phi_j_grad))
 *                    )
 *                    fe_values.JxW(q_point);
 *            }
 *        }
 * @endcode
 * 现在，这不是 step-8 中使用的代码。事实上，如果在该程序中实现的代码之上使用上述代码，它的运行速度会慢8%左右。通过仔细研究双线性形式，它可以得到改善（将惩罚降低到大约4%）。事实上，我们可以将其转换为：@f{eqnarray*}
 * a({\mathbf u}, {\mathbf v})
 * &=&
 * \left(
 *  \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
 * \right)_\Omega
 * +
 * \sum_{i,j}
 * \left(
 *  \mu \partial_i u_j, \partial_i v_j
 * \right)_\Omega
 * +
 * \sum_{i,j}
 * \left(
 *  \mu \partial_i u_j, \partial_j v_i
 * \right)_\Omega
 * \\
 * &=&
 * \left(
 *  \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
 * \right)_\Omega
 * +
 * 2
 * \sum_{i,j}
 * \left(
 *  \mu \partial_i u_j, \frac 12[\partial_i v_j + \partial_j v_i]
 * \right)_\Omega
 * \\
 * &=&
 * \left(
 *  \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
 * \right)_\Omega
 * +
 * 2
 * \sum_{i,j}
 * \left(
 *  \mu \frac 12[\partial_i u_j + \partial_j u_i], \frac 12[\partial_i v_j + \partial_j v_i]
 * \right)_\Omega
 * \\
 * &=&
 * \left(
 *  \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
 * \right)_\Omega
 * +
 * 2
 * \sum_{i,j}
 * \left(
 * \mu \varepsilon(\mathbf u), \varepsilon(\mathbf v)
 * \right)_\Omega,
 * @f} 。
 * 其中 $\varepsilon(\mathbf u) = \frac 12 \left([\nabla\mathbf u] +
 * [\nabla\mathbf u]^T\right)$
 * 是对称的梯度。在第二至最后一步，我们用任意张量
 * $\nabla\mathbf u$ 和对称张量 $\frac 12[\partial_i v_j + \partial_j
 * v_i]$
 * 之间的标量积等于前者的对称部分与第二张量的标量积。使用上面讨论的技术，实现这一点的明显方法是这样的。
 *
 * @code
 *    for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
 *      for (unsigned int i=0; i<dofs_per_cell; ++i)
 *        {
 *          const SymmetricTensor<2,dim> phi_i_symmgrad
 *            = fe_values[displacements].symmetric_gradient (i,q_point);
 *          const double phi_i_div
 *            = fe_values[displacements].divergence (i,q_point);
 *
 *          for (unsigned int j=0; j<dofs_per_cell; ++j)
 *            {
 *              const SymmetricTensor<2,dim> phi_j_symmgrad
 *                = fe_values[displacements].symmetric_gradient (j,q_point);
 *              const double phi_j_div
 *                = fe_values[displacements].divergence (j,q_point);
 *
 *              cell_matrix(i,j)
 *                +=  (phi_i_div phi_j_div
 *                     lambda_values[q_point]
 *                     +
 *                     2
 *                     (phi_i_symmgrad phi_j_symmgrad)
 *                     mu_values[q_point])
 *                    fe_values.JxW(q_point);
 *            }
 *        }
 * @endcode
 *
 * 那么，如果同样，这不是我们在 step-8
 * 中使用的代码，我们在那里做什么？答案取决于我们使用的有限元。在
 * step-8 中，我们使用以下元素。
 * @code
 * FESystem<dim> finite_element (FE_Q<dim>(1), dim);
 * @endcode
 * 换句话说，我们使用的有限元由同一标量元素的 <code>dim</code> 份组成。这就是我们所说的 @ref GlossPrimitive "原始 "
 * 元素：一个可能是矢量值的元素，但每个形状函数正好有一个非零成分。换句话说：如果一个位移形状函数的
 * $x$ -分量是非零的，那么 $y$ -和 $z$
 * -分量必须是零，其他分量也是如此。这意味着基于形状函数的衍生量也继承了这种稀疏性。例如：矢量值形状函数
 * $\Phi(x,y,z)=(\varphi_x(x,y,z), \varphi_y(x,y,z), \varphi_z(x,y,z))^T$
 * 的发散 $\mathrm{div}\ \Phi(x,y,z)=\partial_x\varphi_x(x,y,z) +
 * \partial_y\varphi_y(x,y,z) + \partial_z\varphi_z(x,y,z)$ 在本例中是
 * $\mathrm{div}\ \Phi(x,y,z)=\partial_x\varphi_x(x,y,z)$ 、 $\mathrm{div}\
 * \Phi(x,y,z)=\partial_y\varphi_y(x,y,z)$ 或 $\mathrm{div}\
 * \Phi(x,y,z)=\partial_z\varphi_z(x,y,z)$ ，因为 $\varphi_\ast$
 * 中正好有一个是非零。知道这一点意味着我们可以节省一些计算，如果我们要做这些计算，只会产生零的加法。
 * 类似地，如果一个形状函数只有一个分量是非零的，那么它的梯度
 * $\nabla\Phi$ 就只有一行是非零。这对于像 $(\mu
 * \nabla\Phi_i,\nabla\Phi_j)$
 * 这样的术语意味着什么，其中两个张量之间的标量积被定义为
 * $(\tau, \gamma)_\Omega=\int_\Omega \sum_{i,j=1}^d \tau_{ij} \gamma_{ij}$
 * ，即只有当两个张量的非零项在同一行时，该术语才是非零的，这意味着两个形状函数必须在同一位置有其单一非零分量。
 * 如果我们使用这种知识，那么我们可以在第一步避免计算梯度张量，如果我们可以预先确定它们的标量乘积将是非零的，在第二步避免建立整个张量，只得到它的非零分量，在最后一步简化标量乘积，只考虑一个非零行的索引
 * $i$ ，而不是乘以和增加零。
 * 这一切的载体是确定哪个向量分量将是非零的能力。这个信息是由
 * FiniteElement::system_to_component_index 函数提供的。在  step-8
 * 中详细解释了用它可以做什么，使用上面的例子。
 *
 *  @anchor  VVBlockSolvers <h3>Block solvers</h3>。
 * 使用如上所示的技术，对于一个矢量值问题，组装线性系统，即矩阵和右手边，并不特别复杂。然而，然后它还必须被解决。这就比较复杂了。直观地说，我们可以只把矩阵作为一个整体来考虑。对于大多数问题，这个矩阵不会是确定的（除了特殊情况，如
 * step-8 和 step-17
 * 中涉及的弹性方程）。它通常也不是对称的。这类相当普遍的矩阵给迭代求解器带来了问题：由于缺乏结构特性，无法使用最有效的方法和预处理器。虽然可以做到这一点，但求解过程往往比必要的要慢。
 * *这个问题的答案是利用问题的结构。例如，对于上面讨论的混合拉普拉斯方程，算子的形式是@f{eqnarray*}
 * \left(
 * \begin{array}{cc} \mathbf 1 & \nabla \\
 *
 * -\nabla^T & 0 \end{array}
 * \right)
 * @f} 。
 *
 * 如果这种结构也能在线性系统中恢复，那就更好了。例如，在离散化之后，我们希望有一个具有以下块状结构的矩阵：@f{eqnarray*}
 * \left(
 * \begin{array}{cc} M & B \\ B^T & 0 \end{array}
 * \right),
 * @f}
 * 其中 $M$ 代表离散化身份算子 $\mathbf 1$ 产生的质量矩阵，
 * $B$ 是梯度算子的等价物。
 * 然而，在默认情况下，这并不是发生的情况。相反，deal.II以一种相当随机的方式给自由度分配%的数字。因此，如果你用自由度的值组成一个向量，就不会像@f{eqnarray*}
 * \left(
 * \begin{array}{c} U \\ P \end{array}
 * \right).
 * @f}那样整齐地排列在一个向量中。
 * 相反，它将是一个排列组合，与速度和压力相对应的自由度的%数混合在一起。因此，系统矩阵也不会有上面提到的漂亮结构，而是有相同的排列组合或行和列。
 * 但之后我们仍然要利用它，也就是说，我们必须拿出一个使用该结构的求解器。例如，在 step-20 中，我们对线性系统@f{eqnarray*}
 * \left(
 * \begin{array}{cc} M & B \\ B^T & 0 \end{array}
 * \right)
 * \left(
 * \begin{array}{c} U \\ P \end{array}
 * \right)
 * =
 * \left(
 * \begin{array}{c} F \\ G \end{array}
 * \right).
 * @f}做了一个块消除。
 * 当然，这个系统的含义是@f{eqnarray*}
 * MU + BP &=& F,\\
 * B^TU  &=& G.
 * @f} 。
 *
 * 因此，如果我们用第一个方程乘以 $B^TM^{-1}$ ，再从结果中减去第二个方程，我们就得到@f{eqnarray*}
 * B^TM^{-1}BP &=& B^TM^{-1}F-G.
 * @f} 。
 *
 * 这是一个现在只包含压力变量的方程。如果我们能解决这个问题，我们可以在第二步用@f{eqnarray*}
 * MU = F-BP.
 * @f}来解决速度问题。
 *
 * 这样做的好处是，我们要解决的矩阵 $B^TM^{-1}B$ 和 $M$
 * 都是对称的和正定的，而不是我们之前的大整数矩阵。
 * 像这样的求解器是如何实现的，在 @ref step_20 "  step-20 "
 * 、 step-31
 * 和其他一些教程程序中有更详细的解释。我们在这里想指出的是，我们现在需要一种方法来提取矩阵或向量的某些部分：如果我们要将，比如说，解向量的
 * $U$ 部分与全局矩阵的 $M$
 * 部分相乘，那么我们需要有一种方法来访问整体的这些部分。
 * 这就是BlockVector、BlockSparseMatrix和类似的类的用处。为了所有的实际目的，那么可以作为常规的向量或稀疏矩阵使用，也就是说，它们提供元素访问，提供常规的向量操作，并实现例如矩阵-向量的乘法。换句话说，组装矩阵和右手边的工作方式与非块版本完全相同。也就是说，在内部，它们以
 * "块
 * "的形式存储向量和矩阵的元素；例如，BlockVector类不是使用一个大数组，而是将其存储为一组数组，每个数组我们称之为一个块。这样做的好处是，虽然整个东西可以作为一个向量使用，但人们也可以访问一个单独的块，然后，它又是一个具有所有向量操作的向量。
 * 为了说明如何做到这一点，让我们考虑上面要解决的第二个方程
 * $MU=F-BP$ 。这可以通过以下类似于我们在 step-20
 * 中的序列来实现。
 * @code
 *  Vector<double> tmp (solution.block(0).size());
 *  system_matrix.block(0,1).vmult (tmp, solution.block(1));
 *  tmp=
 *
 * -1;
 *  tmp += system_rhs.block(0);
 *
 *
 *
 *
 *  SolverControl solver_control (solution.block(0).size(),
 *                                1e-8*tmp.l2_norm());
 *  SolverCG<> cg (solver_control, vector_memory);
 *
 *  cg.solve (system_matrix.block(0,0),
 *            solution.block(0),
 *            tmp,
 *            PreconditionIdentity());
 * @endcode
 *
 * 这里发生的事情是，我们分配了一个临时向量，其元素数与解向量的第一块，即速度分量
 * $U$ 相同。然后我们将这个临时向量设置为等于矩阵的
 * $(0,1)$ 块，即 $B$
 * ，乘以解决方案的分量1，即之前计算的压力 $P$
 * 。结果乘以 $-1$ ，右手边的0分量 $F$
 * 被添加到其中。现在的临时向量包含  $F-BP$
 * 。剩下的代码片段只是解决了一个线性系统， $F-BP$
 * 为右手边，全局矩阵的 $(0,0)$ 块，即 $M$
 * 。因此，以这种方式使用块状向量和矩阵，我们可以很容易地编写相当复杂的求解器，利用线性系统的块状结构。
 *
 *   @anchor  VVExtracting <h3>Extracting data from solutions</h3>。
 * 一旦计算出一个解决方案，往往需要在正交点上进行评估，例如为下一次牛顿迭代评估非线性残差，为误差估计器评估有限元残差，或者为时间相关问题的下一个时间步骤计算右手边。
 * 这样做的方法是再次使用FEValues对象来评估正交点的形状函数，并且用这些来评估有限元函数的值。对于上面的混合拉普拉斯问题的例子，请考虑解算后的以下代码。
 * @code
 * std::vector<Vector<double> > local_solution_values (n_q_points,
 *                                                    Vector<double> (dim+1));
 *
 * typename DoFHandler<dim>::active_cell_iterator
 *  cell = dof_handler.begin_active(),
 *  endc = dof_handler.end();
 * for (; cell!=endc; ++cell)
 *  {
 *    fe_values.reinit (cell);
 *
 *    fe_values.get_function_values (solution,
 *                                   local_solution_values);
 * @endcode
 *
 * 在这之后，变量 <code>local_solution_values</code>
 * 是一个长度等于我们初始化FEValues对象的正交点数量的向量列表；每个向量都有
 * <code>dim+1</code> 元素，包含 <code>dim</code>
 * 速度的值和正交点的一个压力。
 * 我们可以用这些值来构建其他的东西，如残差。然而，这个构造有点尴尬。首先，我们有一个
 * <code>std::vector</code> of <code>dealii::Vector</code>
 * s，这看起来总是很奇怪。它也是低效的，因为它意味着为外向量以及所有内向量分配动态内存。其次，也许我们只对速度感兴趣，例如，在第二阶段解决一个平流问题（例如，在
 * step-21 或 step-31
 * 中）。在这种情况下，我们必须像这样手工提取这些值。
 * @code
 * for (unsigned int q=0; q<n_q_points; ++q)
 *   {
 *     Tensor<1,dim> velocity;
 *     for (unsigned int d=0; d<dim; ++d)
 *       velocity[d] = local_solution_values[q](d);
 *
 *     ... do something with this velocity ...
 * @endcode
 * 注意我们如何从 dealii::Vector
 * （它只是一个矢量元素的集合）转换为
 * <code>Tensor@<1,dim@></code> ，因为速度是一个由 <code>dim</code>
 * 元素表征的量，在坐标系的旋转下具有某些变换特性。
 * 这段代码可以用下面这样的代码写得更优雅、更有效。
 * @code
 * std::vector<Tensor<1,dim> > local_velocity_values (n_q_points);
 *
 * const FEValuesExtractors::Vector velocities (0);
 *
 * typename DoFHandler<dim>::active_cell_iterator
 *  cell = dof_handler.begin_active(),
 *  endc = dof_handler.end();
 * for (; cell!=endc; ++cell)
 *  {
 *    fe_values.reinit (cell);
 *
 *    fe_values[velocities].get_function_values (solution,
 *                                               local_velocity_values);
 * @endcode
 *
 * 结果，我们在这里马上得到了速度，而且是正确的数据类型（因为我们已经用提取器描述了有限元的第一个
 * <code>dim</code>
 * 分量属于一起，形成一个张量）。这段代码也更有效率：它需要更少的动态内存分配，因为张量类将其成分作为成员变量而不是在堆上分配，而且我们节省了周期，因为我们甚至不需要费力计算正交点上的压力变量值。另一方面，如果我们只对压力而不是速度感兴趣，那么下面提取标量值的代码就可以了。
 * @code
 * std::vector<double> local_pressure_values (n_q_points);
 *
 * const FEValuesExtractors::Scalar pressure (dim);
 *
 * typename DoFHandler<dim>::active_cell_iterator
 *  cell = dof_handler.begin_active(),
 *  endc = dof_handler.end();
 * for (; cell!=endc; ++cell)
 *  {
 *    fe_values.reinit (cell);
 *
 *    fe_values[pressure].get_function_values (solution,
 *                                             local_pressure_values);
 * @endcode
 *
 * 在类似情况下，有时需要解的梯度或二阶导数，或者个别标量或矢量分量的梯度或二阶导数。为了获得解的所有分量的梯度，函数
 * FEValuesBase::get_function_gradients 和
 * FEValuesBase::get_function_hessians 相当于上面使用的函数
 * FEValuesBase::get_function_values 。
 * 同样，要提取标量分量的梯度，
 * FEValuesViews::Scalar::get_function_gradients 和
 * FEValuesViews::Scalar::get_function_hessians
 * 就可以完成这项工作。对于矢量（张量）值的量，有函数
 * FEValuesViews::Vector::get_function_gradients 和
 * FEValuesViews::Vector::get_function_hessians, ，此外还有
 * FEValuesViews::Vector::get_function_symmetric_gradients 和
 * FEValuesViews::Vector::get_function_divergences.  。
 * 此外，在只需要解的拉普拉斯（即豫备的轨迹）的情况下，还有一个捷径，可用于标量和矢量值的问题，如
 * FEValuesViews::Scalar::get_function_laplacians 和
 * FEValuesViews::Vector::get_function_laplacians.  。
 *
 *  @anchor  VVOutput <h3>Generating graphical output</h3>。
 * 如上所述，一个FESystem对象可能持有多个向量组件，但它并不清楚这些组件的实际含义。作为一个例子，以对象
 * @code
 * FESystem<dim> finite_element (FE_Q<dim>(1), dim+1);
 * @endcode
 * 它有 <code>dim+1</code>
 * 个矢量分量，但它们是什么意思？它们是速度矢量的
 * <code>dim</code> 分量加上一个压力吗？它们是压力加上
 * <code>dim</code> 速度分量吗？还是它们是一个标量的集合？
 * 关键是，FESystem类并不关心。元素的<i>interpretation</i>含义是由后来使用该元素的人决定的，例如在组装线性表格时，或者在下一个牛顿步骤中为线性化系统提取数据解决方案组件时。几乎在所有情况下，这种解释都发生在需要它的地方。
 * 然而，有一种情况是必须明确的，那就是在生成图形输出时。原因是许多用于可视化的文件格式希望表示矢量的数据（如速度、位移等）与标量（压力、密度等）分开存储，而且通常没有办法在可视化程序中把一堆标量分组为一个矢量场。
 * 为了实现这一点，我们需要让DataOut类和朋友们知道FESystem的哪些成分形成了向量（有
 * <code>dim</code> 成分），哪些是标量。例如，这在 step-22
 * 中显示，我们产生的输出如下。
 * @code
 * std::vector<std::string> solution_names (dim, "velocity");
 * solution_names.push_back ("pressure");
 *
 * std::vector<DataComponentInterpretation::DataComponentInterpretation>
 *  data_component_interpretation
 *  (dim, DataComponentInterpretation::component_is_part_of_vector);
 * data_component_interpretation
 *  .push_back (DataComponentInterpretation::component_is_scalar);
 *
 * DataOut<dim> data_out;
 * data_out.attach_dof_handler (dof_handler);
 * data_out.add_data_vector (solution, solution_names,
 *                          DataOut<dim>::type_dof_data,
 *                          data_component_interpretation);
 * data_out.build_patches ();
 * @endcode
 * 换句话说，我们在这里创建了一个 <code>dim+1</code>
 * 元素的数组，在这个数组中，我们存储了有限元中哪些元素是向量，哪些是标量；这个数组被填充了
 * <code>dim</code> 的副本和
 * DataComponentInterpretation::component_is_scalar
 * 的一个尾部元素。然后，该数组被作为一个额外的参数给到
 * DataOut::add_data_vector
 * ，以解释如何解释给定解向量中的数据。像VisIt和Paraview这样的可视化程序将提供显示这些
 * <code>dim</code> 组件的矢量场，而不是单个标量场。
 *
 * @ingroup feall feaccess
 *
 */


