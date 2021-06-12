/**
@page step_8 The step-8 tutorial program
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
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeElasticProblemcodeclasstemplate">The <code>ElasticProblem</code> class template</a>
        <li><a href="#Righthandsidevalues">Right hand side values</a>
        <li><a href="#ThecodeElasticProblemcodeclassimplementation">The <code>ElasticProblem</code> class implementation</a>
      <ul>
        <li><a href="#ElasticProblemElasticProblemconstructor">ElasticProblem::ElasticProblem constructor</a>
        <li><a href="#ElasticProblemsetup_system">ElasticProblem::setup_system</a>
        <li><a href="#ElasticProblemassemble_system">ElasticProblem::assemble_system</a>
        <li><a href="#ElasticProblemsolve">ElasticProblem::solve</a>
        <li><a href="#ElasticProblemrefine_grid">ElasticProblem::refine_grid</a>
        <li><a href="#ElasticProblemoutput_results">ElasticProblem::output_results</a>
        <li><a href="#ElasticProblemrun">ElasticProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-8/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



在现实生活中，大多数偏微分方程实际上是方程组。相应地，解通常是矢量值的。deal.II库支持这样的问题（见 @ref vector_valued 模块中的大量文档），我们将表明这大多是相当简单的。唯一比较复杂的问题是在组装矩阵和右手边，但这些也很容易理解。

 @dealiiVideoLecture{19} 

在这个教程程序中，我们将想解决<a href="https://en.wikipedia.org/wiki/Linear_elasticity">elastic equations</a>。它们是对拉普拉斯方程的扩展，有一个矢量值的解，描述了受力的刚体在每个空间方向的位移。当然，力也是矢量值的，意味着在每一个点上它都有一个方向和一个绝对值。

人们可以用多种方式来写弹性方程。以最明显的方式显示与拉普拉斯方程的对称性的是将其写成

@f[


  -
  \text{div}\,
  ({\mathbf C} \nabla \mathbf{u})
  =
  \mathbf f,


@f]

其中 $\mathbf u$ 是每一点的矢量值位移， $\mathbf f$ 是力， ${\mathbf C}$ 是一个等级4的张量（即它有四个指数），编码应力-应变关系--本质上，它代表胡克斯定律中的<a href="https://en.wikipedia.org/wiki/Hooke%27s_law">"spring constant"</a>，将位移与力联系起来。  在许多情况下，如果我们想要模拟的物体的变形是由不同的材料组成的，那么 ${\mathbf C}$ 将取决于 $\mathbf x$ 。

虽然上述方程的形式是正确的，但这并不是它们通常的推导方式。事实上，位移的梯度 $\nabla\mathbf u$ （一个矩阵）没有物理意义，而其对称版本。

@f[
\varepsilon(\mathbf u)_{kl} =\frac{1}{2}(\partial_k u_l + \partial_l u_k),


@f]

做，通常被称为 "应变"。(在这里和下文中， $\partial_k=\frac{\partial}{\partial x_k}$  。我们还将使用<a href="https://en.wikipedia.org/wiki/Einstein_notation">Einstein summation
convention</a>，即只要同一指数在方程式中出现两次，就意味着对该指数进行求和；但是，我们将不区分上下指数）。)有了这个应变的定义，弹性方程就读作

@f[


  -
  \text{div}\,
  ({\mathbf C} \varepsilon(\mathbf u))
  =
  \mathbf f,


@f]

你可以把它看作是拉普拉斯方程对矢量值问题的更自然的概括。(首先显示的形式等同于这种形式，因为张量 ${\mathbf C}$ 具有某些对称性，即 $C_{ijkl}=C_{ijlk}$  ，因此 ${\mathbf C} \varepsilon(\mathbf u)_{kl}
= {\mathbf C} \nabla\mathbf u$  。)

当然，我们也可以把这些方程写成组件形式。

@f[


  -
  \partial_j (c_{ijkl} \varepsilon_{kl})
  =
  f_i,
  \qquad
  i=1\ldots d.


@f]



在许多情况下，我们知道所考虑的材料是各向同性的，在这种情况下，通过引入两个系数 $\lambda$ 和 $\mu$ ，系数张量减少为

@f[
  c_{ijkl}
  =
  \lambda \delta_{ij} \delta_{kl} +
  \mu (\delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk}).


@f]



然后，弹性方程可以用更简单的形式重写。

@f[


   -
   \nabla \lambda (\nabla\cdot {\mathbf u})


   -
   (\nabla \cdot \mu \nabla) {\mathbf u}


   -
   \nabla\cdot \mu (\nabla {\mathbf u})^T
   =
   {\mathbf f},


@f]

而各自的双线性形式则是

@f[
  a({\mathbf u}, {\mathbf v}) =
  \left(
    \lambda \nabla\cdot {\mathbf u}, \nabla\cdot {\mathbf v}
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega,


@f]

或将第一项写成成分之和。

@f[
  a({\mathbf u}, {\mathbf v}) =
  \sum_{k,l}
  \left(
    \lambda \partial_l u_l, \partial_k v_k
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_k v_l
  \right)_\Omega
  +
  \sum_{k,l}
  \left(
    \mu \partial_k u_l, \partial_l v_k
  \right)_\Omega.


@f]



 @note 按照写法，如果位移很小，我们可以假设<a
href="http://en.wikipedia.org/wiki/Hookes_law">Hooke's law</a>是有效的，上面的方程一般被认为是对三维物体位移的正确描述。在这种情况下，上面的指数 $i,j,k,l$ 都是在集合 $\{1,2,3\}$ 上运行的（或者，在C++源中，在 $\{0,1,2\}$ 上运行）。然而，按照目前的情况，程序是在2d中运行的，虽然上面的方程在这种情况下也有数学意义，但它们只能描述一个真正的二维实体。特别是，它们不是对 $x-y$ 方向上无限大的体的横截面的适当描述；这与其他许多二维方程相反，这些方程可以通过假设体在 $z$ -方向上具有无限大的范围和解函数不依赖于 $z$ 坐标来获得。另一方面，也有二维弹性模型的方程；例如，见维基百科上的<a
href="http://en.wikipedia.org/wiki/Infinitesimal_strain_theory#Special_cases">plane
strain</a>、<a
href="http://en.wikipedia.org/wiki/Antiplane_shear">antiplane shear</a>和<a
href="http://en.wikipedia.org/wiki/Plane_stress#Plane_stress">plan stress</a>文章。

但让我们回到最初的问题上。我们如何为这样一个方程组装矩阵？在 @ref vector_valued 模块的文档中给出了一个很长的答案，其中有许多不同的选择。从历史上看，下面所示的解决方案是该库早期唯一可用的解决方案。事实证明，它也是最快的。另一方面，如果百分之几的计算时间并不重要，还有比下面讨论的更简单、更直观的方法来组装线性系统，但这些方法直到本教程首次编写后的几年才可用；如果你对它们感兴趣，可以看看  @ref vector_valued  模块。

让我们回到如何组装线性系统的问题上来。首先我们需要一些关于形状函数在矢量值有限元情况下如何工作的知识。基本上，这归结为以下几点：让 $n$ 为我们建立矢量元素的标量有限元素的形状函数的数量（例如，我们将对矢量值有限元素的每个分量使用双线性函数，所以标量有限元素是我们在以前的例子中已经使用过的 <code>FE_Q(1)</code> 元素，以及两个空间维度的 $n=4$ ）。此外，让 $N$ 为矢量元素的形状函数数量；在两个空间维度中，我们需要为矢量的每个分量提供 $n$ 个形状函数，因此 $N=2n$  。那么，矢量元素的 $i$ 个形状函数的形式为

@f[
  \Phi_i({\mathbf x}) = \varphi_{\text{base}(i)}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)},


@f]

其中 $e_l$ 是第 $l$ 个单位向量， $\text{comp}(i)$ 是告诉我们 $\Phi_i$ 的哪个分量是不为零的函数（对于每个向量形状函数，只有一个分量是不为零的，其他都是零）。   $\varphi_{\text{base}(i)}(x)$ 描述了形状函数的空间依赖性，它被认为是标量元素的第 $\text{base}(i)$ 个形状函数。当然，虽然 $i$ 的范围是 $0,\ldots,N-1$ ，但函数 $\text{comp}(i)$ 和 $\text{base}(i)$ 的范围分别为 $0,1$ （在二维）和 $0,\ldots,n-1$ 。

例如（尽管这种形状函数的顺序不被保证，你也不应该依赖它），下面的布局可以被库使用。

@f{eqnarray*}
  \Phi_0({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_0({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_1({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_0({\mathbf x})
  \end{array}\right),
  \\
  \Phi_2({\mathbf x}) &=&
  \left(\begin{array}{c}
    \varphi_1({\mathbf x}) \\ 0
  \end{array}\right),
  \\
  \Phi_3({\mathbf x}) &=&
  \left(\begin{array}{c}
    0 \\ \varphi_1({\mathbf x})
  \end{array}\right),
  \ldots


@f}

在这里

@f[
  \text{comp}(0)=0, \quad  \text{comp}(1)=1, \quad  \text{comp}(2)=0, \quad  \text{comp}(3)=1, \quad  \ldots


@f]



@f[
  \text{base}(0)=0, \quad  \text{base}(1)=0, \quad  \text{base}(2)=1, \quad  \text{base}(3)=1, \quad  \ldots


@f]



除了非常罕见的情况，你不需要知道标量元素的哪个形状函数 $\varphi_{\text{base}(i)}$ 属于矢量元素的一个形状函数 $\Phi_i$ 。因此，让我们定义

@f[
  \phi_i = \varphi_{\text{base}(i)}


@f]

据此，我们可以将矢量形状函数写为

@f[
  \Phi_i({\mathbf x}) = \phi_{i}({\mathbf x})\ {\mathbf e}_{\text{comp}(i)}.


@f]

现在你可以安全地忘记函数 $\text{base}(i)$ 了，至少在这个例子程序的其余部分。

现在使用这个矢量形状函数，我们可以将离散的有限元解写为

@f[
  {\mathbf u}_h({\mathbf x}) =
  \sum_i \Phi_i({\mathbf x})\ U_i


@f]

具有标量系数  $U_i$  。如果我们定义一个模拟函数 ${\mathbf v}_h$ 作为测试函数，我们可以将离散问题写成如下。找出系数 $U_i$ ，使得

@f[
  a({\mathbf u}_h, {\mathbf v}_h) = ({\mathbf f}, {\mathbf v}_h)
  \qquad
  \forall {\mathbf v}_h.


@f]



如果我们把双线性形式的定义和 ${\mathbf u}_h$ 和 ${\mathbf v}_h$ 的表示插入这个公式。

@f{eqnarray*}
  \sum_{i,j}
    U_i V_j
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_\Omega
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_\Omega
  \right\}
\\
=
  \sum_j V_j
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_\Omega.


@f}

我们注意到，在这里和下文中，指数 $k,l$ 在空间方向上运行，即 $0\le k,l < d$  ，而指数 $i,j$ 在自由度上运行。

因此，单元 $K$ 上的局部刚度矩阵有以下条目。

@f[
  A^K_{ij}
  =
  \sum_{k,l}
  \left\{
  \left(
    \lambda \partial_l (\Phi_i)_l, \partial_k (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_l (\Phi_j)_k
  \right)_K
  +
  \left(
    \mu \partial_l (\Phi_i)_k, \partial_k (\Phi_j)_l
  \right)_K
  \right\},


@f]

其中 $i,j$ 现在是局部自由度，因此 $0\le i,j < N$  。在这些公式中，我们总是取矢量形状函数 $\Phi_i$ 的一些分量，当然，这些分量是如下给出的（见其定义）。

@f[
  (\Phi_i)_l = \phi_i \delta_{l,\text{comp}(i)},


@f]

与克朗克符号  $\delta_{nm}$  。由于这一点，我们可以删除一些对  $k$  和  $l$  的和。

@f{eqnarray*}
  A^K_{ij}
  &=&
  \sum_{k,l}
  \Bigl\{
  \left(
    \lambda \partial_l \phi_i\ \delta_{l,\text{comp}(i)},
            \partial_k \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
\\
  &\qquad\qquad& +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_l \phi_j\ \delta_{k,\text{comp}(j)}
  \right)_K
  +
  \left(
    \mu \partial_l \phi_i\ \delta_{k,\text{comp}(i)},
        \partial_k \phi_j\ \delta_{l,\text{comp}(j)}
  \right)_K
  \Bigr\}
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \sum_l
  \left(
    \mu \partial_l \phi_i,
        \partial_l \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K
\\
  &=&
  \left(
    \lambda \partial_{\text{comp}(i)} \phi_i,
            \partial_{\text{comp}(j)} \phi_j
  \right)_K
  +
  \left(
    \mu \nabla \phi_i,
        \nabla \phi_j
  \right)_K
  \ \delta_{\text{comp}(i),\text{comp}(j)}
  +
  \left(
    \mu \partial_{\text{comp}(j)} \phi_i,
        \partial_{\text{comp}(i)} \phi_j
  \right)_K.


@f}



同样地，单元格 $K$ 对右侧向量的贡献是

@f{eqnarray*}
  f^K_j
  &=&
  \sum_l
  \left(
    f_l,
    (\Phi_j)_l
  \right)_K
\\
  &=&
  \sum_l
  \left(
    f_l,
    \phi_j \delta_{l,\text{comp}(j)}
  \right)_K
\\
  &=&
  \left(
    f_{\text{comp}(j)},
    \phi_j
  \right)_K.


@f}



这就是我们要实现局部刚度矩阵和右手边向量的形式。

作为最后的说明：在第17步的例子程序中，我们将重新审视这里提出的弹性问题，并将展示如何在一个计算机集群上以%并行的方式解决这个问题。因此，所产生的程序将能够以更高的精度解决这个问题，而且如果需要的话，效率更高。此外，在第20步， @ref step_21 "第21步"，以及其他一些后来的教程程序中，我们将重新审视一些矢量值问题，并展示一些技术，这些技术可能使实际通过上面显示的所有东西更简单，与 FiniteElement::system_to_component_index 等。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
 * 

 * 
 * 像往常一样，前几个include文件已经知道了，所以我们将不再评论它们。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/tensor.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 在这个例子中，我们需要矢量值的有限元。对这些的支持可以在下面的include文件中找到。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_system.h> 
 * 
 * @endcode
 * 
 * 我们将用常规的Q1元素组成矢量值的有限元素，这些元素可以在这里找到，像往常一样。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_q.h> 
 * 
 * @endcode
 * 
 * 这又是C++语言。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后一步和以前的程序一样。特别是，就像在 step-7 中一样，我们把这个程序所特有的一切都打包到一个自己的命名空间中。
 * 

 * 
 * 
 * @code
 * namespace Step8 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclasstemplate"></a> 
 * <h3>The <code>ElasticProblem</code> class template</h3>
 * 

 * 
 * 主类除了名称外，与 step-6 的例子相比几乎没有变化。
 * 

 * 
 * 唯一的变化是为 <code>fe</code> 变量使用了一个不同的类。我们现在使用的不是FE_Q这样具体的有限元类，而是一个更通用的类，FESystem。事实上，FESystem本身并不是一个真正的有限元，因为它没有实现自己的形状函数。相反，它是一个可以用来将其他几个元素堆叠在一起形成一个矢量值的有限元的类。在我们的例子中，我们将组成 <code>FE_Q(1)</code> 对象的矢量值元素，如下所示，在这个类的构造函数中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ElasticProblem 
 *   { 
 *   public: 
 *     ElasticProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void refine_grid(); 
 *     void output_results(const unsigned int cycle) const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     FESystem<dim> fe; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 * 
 *     Vector<double> solution; 
 *     Vector<double> system_rhs; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Righthandsidevalues"></a> 
 * <h3>Right hand side values</h3>
 * 

 * 
 * 在进入主类的实现之前，我们声明并定义描述右手边的函数。这一次，右手边是向量值，解决方案也是如此，所以我们将更详细地描述为此所需的变化。
 * 

 * 
 * 为了防止出现返回向量没有被设置成正确大小的情况，我们对这种情况进行了测试，否则将在函数的开始部分抛出一个异常。请注意，强制输出参数已经具有正确的大小是deal.II中的一个惯例，并且几乎在所有地方都强制执行。原因是，否则我们将不得不在函数开始时检查，并可能改变输出向量的大小。这很昂贵，而且几乎总是不必要的（对函数的第一次调用会将向量设置为正确的大小，随后的调用只需要做多余的检查）。此外，如果我们不能依赖向量已经具有正确大小的假设，那么检查和可能调整向量大小的操作是不能被删除的；这与Assert调用是一个契约，如果程序在优化模式下编译，Assert调用将被完全删除。
 * 

 * 
 * 同样，如果由于某种意外，有人试图在只有一个空间维度的情况下编译和运行程序（在这种情况下，弹性方程没有什么意义，因为它们还原为普通的拉普拉斯方程），我们在第二个断言中终止程序。然而，该程序在三维空间中也能正常工作。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void right_hand_side(const std::vector<Point<dim>> &points, 
 *                        std::vector<Tensor<1, dim>> &  values) 
 *   { 
 *     Assert(values.size() == points.size(), 
 *            ExcDimensionMismatch(values.size(), points.size())); 
 *     Assert(dim >= 2, ExcNotImplemented()); 
 * 
 * @endcode
 * 
 * 该函数的其余部分实现了计算力值。我们将使用一个位于(0.5,0)和(-0.5,0)点周围的两个小圆圈（或球体，在3D中）的X方向的恒定（单位）力，以及位于原点周围的Y方向的力；在3D中，这些中心的Z分量也是零。
 * 

 * 
 * 为此，让我们首先定义两个对象，表示这些区域的中心。请注意，在构建点对象时，所有的分量都被设置为零。
 * 

 * 
 * 
 * @code
 *     Point<dim> point_1, point_2; 
 *     point_1(0) = 0.5; 
 *     point_2(0) = -0.5; 
 * 
 *     for (unsigned int point_n = 0; point_n < points.size(); ++point_n) 
 *       { 
 * 
 * @endcode
 * 
 * 如果 <code>points[point_n]</code> 处于围绕这些点之一的半径为0.2的圆（球）中，那么将X方向的力设置为1，否则为0。
 * 

 * 
 * 
 * @code
 *         if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) || 
 *             ((points[point_n] - point_2).norm_square() < 0.2 * 0.2)) 
 *           values[point_n][0] = 1.0; 
 *         else 
 *           values[point_n][0] = 0.0; 
 * 
 * @endcode
 * 
 * 同样地，如果 <code>points[point_n]</code> 在原点附近，那么将y力设置为1，否则为0。
 * 

 * 
 * 
 * @code
 *         if (points[point_n].norm_square() < 0.2 * 0.2) 
 *           values[point_n][1] = 1.0; 
 *         else 
 *           values[point_n][1] = 0.0; 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeElasticProblemcodeclassimplementation"></a> 
 * <h3>The <code>ElasticProblem</code> class implementation</h3>
 * 
 * <a name="ElasticProblemElasticProblemconstructor"></a> 
 * <h4>ElasticProblem::ElasticProblem constructor</h4>
 * 

 * 
 * 下面是主类的构造函数。如前所述，我们想构造一个由多个标量有限元组成的矢量值有限元（即，我们想构造矢量值元素，使其每个矢量成分都由一个标量元素的形状函数组成）。当然，我们想堆叠在一起的标量有限元的数量等于解函数的分量数量，由于我们考虑每个空间方向上的位移，所以是 <code>dim</code> 。FESystem类可以处理这个问题：我们传递给它我们想组成系统的有限元，以及它的重复频率。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ElasticProblem<dim>::ElasticProblem() 
 *     : dof_handler(triangulation) 
 *     , fe(FE_Q<dim>(1), dim) 
 *   {} 
 * 
 * @endcode
 * 
 * 事实上，FESystem类还有几个构造函数，可以进行更复杂的操作，而不仅仅是将几个相同类型的标量有限元堆叠在一起；我们将在后面的例子中了解这些可能性。
 * 

 * 
 * 
 * <a name="ElasticProblemsetup_system"></a> 
 * <h4>ElasticProblem::setup_system</h4>
 * 

 * 
 * 设置方程组与 step-6 例子中使用的函数相同。DoFHandler类和这里使用的所有其他类都完全知道我们要使用的有限元是矢量值的，并且照顾到了有限元本身的矢量值。(事实上，它们不知道，但这不需要困扰你：因为它们只需要知道每个顶点、直线和单元有多少个自由度，它们不问它们代表什么，也就是说，考虑的有限元是矢量值的，还是例如在每个顶点上有几个自由度的标量Hermite元)。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 * 
 *     constraints.clear(); 
 *     DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              Functions::ZeroFunction<dim>(dim), 
 *                                              constraints); 
 *     constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, 
 *                                     dsp, 
 *                                     constraints, 
 *                                     /*keep_constrained_dofs =  */ false);
 * 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemassemble_system"></a> 
 * <h4>ElasticProblem::assemble_system</h4>
 * 

 * 
 * 这个程序中最大的变化是创建矩阵和右手边，因为它们是取决于问题的。我们将一步一步地完成这个过程  step-  ，因为它比以前的例子要复杂一些。
 * 

 * 
 * 然而，这个函数的前几部分和以前一样：设置一个合适的正交公式，为我们使用的（矢量值）有限元以及正交对象初始化一个FEValues对象，并声明了一些辅助数组。此外，我们还声明了永远相同的两个缩写。  <code>n_q_points</code>  和  <code>dofs_per_cell</code>  。每个单元的自由度数量，我们现在显然是从组成的有限元中询问，而不是从底层的标量Q1元中询问。在这里，它是 <code>dim</code> 乘以Q1元素的每个单元的自由度数，尽管这不是我们需要关心的明确知识。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::assemble_system() 
 *   { 
 *     QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 正如前面的例子所示，我们需要一个地方来存储单元格上所有正交点的系数值。在目前的情况下，我们有两个系数，lambda和mu。
 * 

 * 
 * 
 * @code
 *     std::vector<double> lambda_values(n_q_points); 
 *     std::vector<double> mu_values(n_q_points); 
 * 
 * @endcode
 * 
 * 好吧，我们也可以省略上面的两个数组，因为我们将对lambda和mu使用常数系数，可以这样声明。它们都代表函数总是返回常量值1.0。尽管我们可以在矩阵的组合中省略各自的系数，但为了演示，我们在这里使用它们。
 * 

 * 
 * 
 * @code
 *     Functions::ConstantFunction<dim> lambda(1.), mu(1.); 
 * 
 * @endcode
 * 
 * 和上面的两个常量函数一样，我们将在每个单元格中只调用一次函数right_hand_side，以使事情更简单。
 * 

 * 
 * 
 * @code
 *     std::vector<Tensor<1, dim>> rhs_values(n_q_points); 
 * 
 * @endcode
 * 
 * 现在我们可以开始对所有单元格进行循环。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix = 0; 
 *         cell_rhs    = 0; 
 * 
 *         fe_values.reinit(cell); 
 * 
 * @endcode
 * 
 * 接下来我们得到正交点的系数值。同样，对于右手边也是如此。
 * 

 * 
 * 
 * @code
 *         lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
 *         mu.value_list(fe_values.get_quadrature_points(), mu_values); 
 *         right_hand_side(fe_values.get_quadrature_points(), rhs_values); 
 * 
 * @endcode
 * 
 * 然后将局部刚度矩阵的条目和右手边的向量组合起来。这几乎是一对一地遵循本例介绍中描述的模式。 在位的几个评论之一是，我们可以计算数字  <code>comp(i)</code>  ，即使用下面的  <code>fe.system_to_component_index(i).first</code>  函数调用形状函数  <code>i</code>  的唯一非零向量成分的索引。
 * 

 * 
 * （通过访问 <code>system_to_component_index</code> 函数返回值的 <code>first</code> 变量，你可能已经猜到其中还有更多的内容。事实上，该函数返回一个 <code>std::pair@<unsigned int，无符号int @></code>, ，其中第一个元素是 <code>comp(i)</code> ，第二个元素是介绍中也指出的值 <code>base(i)</code> ，即这个形状函数在这个组件中所有非零的形状函数中的索引，即介绍中的字典 <code>base(i)</code> 。不过，这不是我们通常感兴趣的数字）。)
 * 

 * 
 * 有了这些知识，我们就可以把局部矩阵的贡献集合起来。
 * 

 * 
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices()) 
 *           { 
 *             const unsigned int component_i = 
 *               fe.system_to_component_index(i).first; 
 * 
 *             for (const unsigned int j : fe_values.dof_indices()) 
 *               { 
 *                 const unsigned int component_j = 
 *                   fe.system_to_component_index(j).first; 
 * 
 *                 for (const unsigned int q_point : 
 *                      fe_values.quadrature_point_indices()) 
 *                   { 
 *                     cell_matrix(i, j) += 
 * 
 * @endcode
 * 
 * 第一个项是  $\lambda \partial_i u_i, \partial_j v_j) + (\mu \partial_i u_j, \partial_j v_i)$  。注意， <code>shape_grad(i,q_point)</code> 返回正交点q_point处第i个形状函数的唯一非零分量的梯度。梯度的分量 <code>comp(i)</code> 是第i个形状函数的唯一非零矢量分量相对于comp(i)th坐标的导数，由附加的括号访问。
 * 

 * 
 * 
 * @code
 *                       (                                                  // 
 *                         (fe_values.shape_grad(i, q_point)[component_i] * // 
 *                          fe_values.shape_grad(j, q_point)[component_j] * // 
 *                          lambda_values[q_point])                         // 
 *                         +                                                // 
 *                         (fe_values.shape_grad(i, q_point)[component_j] * // 
 *                          fe_values.shape_grad(j, q_point)[component_i] * // 
 *                          mu_values[q_point])                             // 
 *                         +                                                // 
 * 
 * @endcode
 * 
 * 第二个项是  $(\mu \nabla u_i, \nabla v_j)$  。我们不需要访问梯度的具体分量，因为我们只需要计算两个梯度的标量乘积，这个问题由<tt>operator*</tt>的重载版本来负责，就像前面的例子一样。                            注意，通过使用<tt>?:</tt>操作符，我们只在<tt>component_i</tt>等于<tt>component_j</tt>时才这样做，否则会加上一个零（编译器会将其优化掉）。
 * 

 * 
 * 
 * @code
 *                         ((component_i == component_j) ?        // 
 *                            (fe_values.shape_grad(i, q_point) * // 
 *                             fe_values.shape_grad(j, q_point) * // 
 *                             mu_values[q_point]) :              // 
 *                            0)                                  // 
 *                         ) *                                    // 
 *                       fe_values.JxW(q_point);                  // 
 *                   } 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 组装右手边也和介绍中讨论的一样。
 * 

 * 
 * 
 * @code
 *         for (const unsigned int i : fe_values.dof_indices()) 
 *           { 
 *             const unsigned int component_i = 
 *               fe.system_to_component_index(i).first; 
 * 
 *             for (const unsigned int q_point : 
 *                  fe_values.quadrature_point_indices()) 
 *               cell_rhs(i) += fe_values.shape_value(i, q_point) * 
 *                              rhs_values[q_point][component_i] * 
 *                              fe_values.JxW(q_point); 
 *           } 
 * 
 * @endcode
 * 
 * 从局部自由度到全局矩阵和右手向量的转移不取决于所考虑的方程，因此与之前所有的例子相同。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global( 
 *           cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemsolve"></a> 
 * <h4>ElasticProblem::solve</h4>
 * 

 * 
 * 解算器并不关心方程组的来源，只要它保持正定和对称（这是使用CG解算器的要求），而这个方程组确实是这样。因此，我们不需要改变任何东西。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::solve() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-12); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *     constraints.distribute(solution); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrefine_grid"></a> 
 * <h4>ElasticProblem::refine_grid</h4>
 * 

 * 
 * 对网格进行细化的函数与 step-6 的例子相同。正交公式再次适应了线性元素。请注意，误差估计器默认情况下是将从有限元解的所有分量中得到的估计值相加，也就是说，它使用所有方向的位移，权重相同。如果我们希望网格只适应x方向的位移，我们可以给函数传递一个额外的参数，告诉它这样做，而不考虑其他所有方向的位移作为误差指标。然而，对于目前的问题，似乎应该考虑所有的位移分量，而且权重相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::refine_grid() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     KellyErrorEstimator<dim>::estimate(dof_handler, 
 *                                        QGauss<dim - 1>(fe.degree + 1), 
 *                                        {}, 
 *                                        solution, 
 *                                        estimated_error_per_cell); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.03); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemoutput_results"></a> 
 * <h4>ElasticProblem::output_results</h4>
 * 

 * 
 * 输出的情况与之前的例子中已经显示过的差不多了。唯一的区别是，求解函数是矢量值的。DataOut类会自动处理这个问题，但我们必须给求解向量的每个分量一个不同的名字。
 * 

 * 
 * 为了做到这一点， DataOut::add_vector() 函数想要一个字符串的向量。由于分量的数量与我们工作的维数相同，我们使用下面的 <code>switch</code> 语句。
 * 

 * 
 * 我们注意到，一些图形程序对变量名称中允许的字符有限制。因此，deal.II只支持所有程序都支持的这些字符的最小子集。基本上，这些字符是字母、数字、下划线和其他一些字符，但特别是没有空格和减号/横线。否则该库将抛出一个异常，至少在调试模式下是这样。
 * 

 * 
 * 在列出了1d、2d和3d的情况后，如果我们遇到一个我们没有考虑到的情况，让程序死亡是一种很好的风格。请记住，如果第一个参数中的条件没有得到满足，Assert宏会产生一个异常。当然，条件 <code>false</code> 永远不可能被满足，所以只要程序运行到默认语句，就会中止。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::output_results(const unsigned int cycle) const 
 *   { 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 *     std::vector<std::string> solution_names; 
 *     switch (dim) 
 *       { 
 *         case 1: 
 *           solution_names.emplace_back("displacement"); 
 *           break; 
 *         case 2: 
 *           solution_names.emplace_back("x_displacement"); 
 *           solution_names.emplace_back("y_displacement"); 
 *           break; 
 *         case 3: 
 *           solution_names.emplace_back("x_displacement"); 
 *           solution_names.emplace_back("y_displacement"); 
 *           solution_names.emplace_back("z_displacement"); 
 *           break; 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 * @endcode
 * 
 * 在为解向量的不同组成部分设置了名称之后，我们可以将解向量添加到计划输出的数据向量列表中。请注意，下面的函数需要一个字符串向量作为第二个参数，而我们在以前所有例子中使用的函数在那里接受一个字符串。(事实上，我们之前使用的函数会将单个字符串转换成只有一个元素的向量，并将其转发给另一个函数)。
 * 

 * 
 * 
 * @code
 *     data_out.add_data_vector(solution, solution_names); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output("solution-" + std::to_string(cycle) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ElasticProblemrun"></a> 
 * <h4>ElasticProblem::run</h4>
 * 

 * 
 * <code>run</code> 函数所做的事情与 step-6 中的相同，比如说。这一次，我们使用平方[-1,1]^d作为域，在开始第一次迭代之前，我们在全局上对其进行了四次细化。
 * 

 * 
 * 细化的原因有点意外：我们使用QGauss正交公式，在每个方向上有两个点用于整合右手边；这意味着每个单元上有四个正交点（在二维）。如果我们只对初始网格进行一次全局细化，那么在域上每个方向上就只有四个正交点。然而，右侧函数被选择为相当局部的，在这种情况下，纯属偶然，恰好所有的正交点都位于右侧函数为零的点上（用数学术语来说，正交点恰好在右侧函数的<i>support</i>之外的点上）。这样一来，用正交计算的右手向量将只包含零（尽管如果我们完全用积分计算右手向量的话，它当然会是非零的），方程组的解就是零向量，也就是一个处处为零的有限元函数。从某种意义上说，我们不应该对这种情况的发生感到惊讶，因为我们选择了一个完全不适合手头问题的初始网格。
 * 

 * 
 * 不幸的是，如果离散解是常数，那么KellyErrorEstimator类计算的误差指标对每个单元来说也是零，对 Triangulation::refine_and_coarsen_fixed_number() 的调用将不会标记任何单元进行细化（如果每个单元的指示误差是零，为什么要这样做？因此，下一次迭代中的网格也将只由四个单元组成，同样的问题再次发生。
 * 

 * 
 * 结论是：虽然我们当然不会把初始网格选择得非常适合问题的精确解决，但我们至少必须选择它，使它有机会捕捉到解决方案的重要特征。在这种情况下，它需要能够看到右手边的情况。因此，我们进行了四次全局细化。(任何更大的全局细化步骤当然也可以。)
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ElasticProblem<dim>::run() 
 *   { 
 *     for (unsigned int cycle = 0; cycle < 8; ++cycle) 
 *       { 
 *         std::cout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *         if (cycle == 0) 
 *           { 
 *             GridGenerator::hyper_cube(triangulation, -1, 1); 
 *             triangulation.refine_global(4); 
 *           } 
 *         else 
 *           refine_grid(); 
 * 
 *         std::cout << "   Number of active cells:       " 
 *                   << triangulation.n_active_cells() << std::endl; 
 * 
 *         setup_system(); 
 * 
 *         std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *                   << std::endl; 
 * 
 *         assemble_system(); 
 *         solve(); 
 *         output_results(cycle); 
 *       } 
 *   } 
 * } // namespace Step8 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 在上面最后一行关闭了 <code>Step8</code> 命名空间后，下面是程序的主要功能，又和 step-6 中一模一样（当然，除了改变了类名）。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       Step8::ElasticProblem<2> elastic_problem_2d; 
 *       elastic_problem_2d.run(); 
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
 * 
 * @endcode
examples/step-8/doc/results.dox



<a name="Results"></a><h1>Results</h1>



关于这个程序的结果，除了它们看起来很好之外，没有什么可说的。所有图片都是用VisIt从程序写入磁盘的输出文件中制作的。前两张图片显示了 $x$ -和 $y$ -位移的标量分量。

 <table width="100%">
<tr>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.x.png" alt="">
</td>
<td>
<img src="https://www.dealii.org/images/steps/developer/step-8.y.png" alt="">
</td>
</tr>
</table> 


你可以清楚地看到 $x$ 周围的位移 $x=0.5$ 和 $x=-0.5$ 的来源，以及 $y$ 在原点的位移。

人们经常想做的是将位移显示为一个矢量场，也就是说，每一个点的矢量都说明了位移的方向和大小。不幸的是，这就有点麻烦了。为了理解为什么会这样，请记住，我们刚刚将我们的有限元定义为两个分量的集合（在 <code>dim=2</code> 维度）。我们没有说过这不仅仅是一个压力和一个浓度（两个标量），而是说这两个分量实际上是一个矢量值量的一部分，即位移。如果没有这方面的知识，DataOut类就会假定我们打印的所有单个变量都是独立的标量，然后VisIt和Paraview就会忠实地假定这确实是这样的。换句话说，一旦我们把数据写成标量，这些程序中就没有任何东西可以让我们把这两个标量字段粘贴到一起作为一个矢量字段。我们必须从根本上解决这个问题，即在  <code>ElasticProblem::output_results()</code>  。我们不会在这里这样做，而是让读者参考step-22程序，在那里我们展示了如何在一个更普遍的情况下这样做。话虽如此，我们还是忍不住要生成数据，以显示如果按照步骤22中讨论的方式实施，这将是什么样子。矢量场看起来是这样的（VisIt和Paraview随机选择几百个顶点来绘制矢量；从每个顶点绘制矢量会使图片无法阅读）。

 <img src="https://www.dealii.org/images/steps/developer/step-8.vectors.png" alt=""> 


我们注意到，由于 $x$ -和 $y$ -力相对于这些轴是对称的，人们可能直观地期望解是关于 $x$ -和 $y$ -轴的对称。然而，作为矢量的力是不对称的，因此解决方案也不对称。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-8.cc"
*/
