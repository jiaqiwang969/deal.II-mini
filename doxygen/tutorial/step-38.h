/**
@page step_38 The step-38 tutorial program
This tutorial depends on step-34.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Testcase">Testcase</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeLaplaceBeltramiProblemcodeclasstemplate">The <code>LaplaceBeltramiProblem</code> class template</a>
        <li><a href="#Equationdata">Equation data</a>
        <li><a href="#ImplementationofthecodeLaplaceBeltramiProblemcodeclass">Implementation of the <code>LaplaceBeltramiProblem</code> class</a>
      <ul>
        <li><a href="#LaplaceBeltramiProblemmake_grid_and_dofs">LaplaceBeltramiProblem::make_grid_and_dofs</a>
        <li><a href="#LaplaceBeltramiProblemassemble_system">LaplaceBeltramiProblem::assemble_system</a>
        <li><a href="#LaplaceBeltramiProblemsolve">LaplaceBeltramiProblem::solve</a>
        <li><a href="#LaplaceBeltramiProblemoutput_result">LaplaceBeltramiProblem::output_result</a>
        <li><a href="#LaplaceBeltramiProblemcompute_error">LaplaceBeltramiProblem::compute_error</a>
        <li><a href="#LaplaceBeltramiProblemrun">LaplaceBeltramiProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main() function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-38/doc/intro.dox

 <br> 

<i>This program was contributed by Andrea Bonito and M. Sebastian Pauletti,
with editing and writing by Wolfgang Bangerth.
<br>
This material is based upon work supported by the National Science
Foundation under Grant No. DMS-0914977. Any opinions, findings and conclusions
or recommendations expressed in this material are those of the author(s) and
do not necessarily reflect the views of the National Science Foundation
(NSF).
</i>

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


在这个例子中，我们展示了如何解决由四边形组成的一维曲面 $\Gamma \subset \mathbb R^3$ 上的偏微分方程（PDE），即在三维的曲面或二维的直线上。我们重点讨论以下的椭圆二阶PDE

@f{align*}


-\Delta_\Gamma u &= f \qquad \text{on } \qquad \Gamma,\\
u  &= g \qquad \text{on} \qquad \partial \Gamma,


@f}

它概括了我们以前在几个早期教程程序中解决的拉普拉斯方程。我们的实现是基于step-4的。step-34也可以解决低维曲面上的问题；但是，在那里我们只考虑不涉及解变量导数的积分方程，而在这里我们实际上要研究只在一个（可能是弯曲的）曲面上定义的函数的导数是什么意思。

为了定义上述算子，我们首先要介绍一些符号。让 $\mathbf x_S:\hat S \rightarrow S$ 是一个由参考元素 $\hat S \subset \mathbb R^2$ 构成的曲面 $S$ 的参数化，即每个点 $\hat{\mathbf x}\in\hat S$ 诱导出一个点 ${\mathbf
  x}_S(\hat{\mathbf x}) \in S$  。那么让

@f[
G_S\dealcoloneq (D \mathbf{x}_S)^T \ D \mathbf{x}_S


@f]

表示相应的第一基本形式，其中 $D
\mathbf{x}_S=\left(\frac{\partial x_{S,i}(\hat{\mathbf x})}{\partial \hat x_j}\right)_{ij}$ 是映射的导数（雅各布）。在下文中， $S$ 将是整个表面 $\Gamma$ ，或者对有限元方法更方便的是任何面 $S \in
{\mathbb T}$ ，其中 ${\mathbb T}$ 是由四边形构成的 $\Gamma$ 的分区（三角化）。我们现在可以定义一个函数 $v : S \rightarrow \mathbb
R$ 的切向梯度，即

@f[
(\nabla_S v)\circ \mathbf x_S \dealcoloneq  D \mathbf x_S \ G_S^{-1} \ \nabla (v \circ \mathbf x_S).


@f]

表面拉普拉斯(也叫拉普拉斯-贝特拉米算子)的定义是  $\Delta_S \dealcoloneq \nabla_S \cdot \nabla_S$  。请注意，在光滑表面上计算表面梯度的另一种方法  $\Gamma$  是

@f[
\nabla_S v = \nabla \tilde v - \mathbf n (\mathbf n \cdot \nabla \tilde v),


@f]

其中 $\tilde v$ 是 $v$ 在 $\Gamma$ 的管状邻域的 "平滑 "扩展， $\mathbf n$ 是 $\Gamma$ 的法线。由于 $\Delta_S = \nabla_S \cdot \nabla_S$ ，我们推导出

@f[
\Delta_S v = \Delta \tilde v - \mathbf n^T \ D^2 \tilde v \ \mathbf n - (\mathbf n \cdot \nabla \tilde v) (\nabla \cdot \mathbf n - \mathbf n^T \ D \mathbf n \ \mathbf n ).


@f]

值得一提的是，上述表达式中出现的术语 $\nabla \cdot \mathbf n - \mathbf n \ D \mathbf n \ \mathbf n$ 是曲面的总曲率（主曲率之和）。

像往常一样，我们只对弱解感兴趣，为此我们可以使用 $C^0$ 有限元（而不是像强解那样要求 $C^1$ 的连续性）。因此，我们求助于弱的表述

@f[
\int_\Gamma \nabla_\Gamma u \cdot
\nabla_\Gamma v = \int_\Gamma f \ v  \qquad \forall v \in H^1_0(\Gamma)


@f]

并利用分区 ${\mathbb T}$ 的优势，进一步编写

@f[
\sum_{K\in  {\mathbb T}}\int_K \nabla_{K} u \cdot \nabla_{K} v = \sum_{K\in
  {\mathbb T}} \int_K f \ v  \qquad \forall v \in H^1_0(\Gamma).


@f]

此外，上述表达式中的每个积分都是在参考元素 $\hat K \dealcoloneq [0,1]^2$ 中计算的，因此

@f{align*}
\int_{K} \nabla_{K} u \cdot \nabla_{K} v
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} (D \mathbf
  x_K)^T D \mathbf x_K G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}
\\
&=
\int_{\hat K} \nabla (u \circ \mathbf x_K)^T G_K^{-1} \nabla (v \circ \mathbf x_K) \sqrt{\det
    (G_K)}


@f}

和

@f[
\int_{K} f \ v = \int_{\hat K} (f \circ \mathbf x_K) (v \circ \mathbf
x_K)  \sqrt{\det
    (G_K)}.


@f]

最后，我们使用由点 $\{p_l\}_{l=1}^N\subset
\hat K$ 和权重 $\{w_l\}_{l=1}^N \subset \mathbb R^+_*$ 定义的正交公式来评估上述积分，得到

@f[\int_{K} \nabla_{K} u \cdot \nabla_{K} v \approx \sum_{l=1}^N
 (\nabla (u \circ \mathbf x_K)(p_l))^T G^{-1}(p_l)  \nabla (v \circ \mathbf x_K)
(p_l) \sqrt{\det (G(p_l))} \ w_l


@f]

和

@f[
\int_{K} f \ v \approx \sum_{l=1}^N (f \circ \mathbf x_K)(p_l) \ (v \circ \mathbf x_K)(p_l) \sqrt{\det (G(p_l))} \ w_l.


@f]




幸运的是，deal.II已经有了所有的工具来计算上述表达式。事实上，它们与我们求解通常的拉普拉斯的方法几乎没有区别，只需要在FEValues类的构造函数中提供表面坐标映射。这个曲面描述给定，在二维曲面的情况下，两个例程 FEValues::shape_grad 和 FEValues::JxW 会返回

@f{align*}
\text{FEValues::shape\_grad}(i,l)&=D \mathbf x_K(p_l) G^{-1}(p_l)D(\varphi_i \circ \mathbf x_K)
  (p_l)
\\
\text{FEValues::JxW}(l) &=  \sqrt{\det (G(p_l))} \ w_l.


@f}

这正好提供了我们的计算所需的术语。

在更广泛的意义上，表面有限元逼近的细节可以在[Dziuk, in Partial differential equations and calculus of variations 1357, Lecture Notes in Math., 1988], [Demlow, SIAM J. Numer. Anal. 47(2), 2009] 和 [Bonito, Nochetto, and Pauletti, SIAM J. Numer. Anal. 48(5), 2010] 中找到。




<a name="Testcase"></a><h3>Testcase</h3>


一般来说，当你想在数值上测试一个算法的准确性和/或收敛性，你需要提供一个精确的解决方案。通常的技巧是选择一个我们希望成为解决方案的函数，然后对其应用微分算子，为右侧定义一个强制项。这就是我们在这个例子中所做的。在当前情况下，域的形式显然也是至关重要的。

我们为二维问题制作一个测试案例，为三维问题制作另一个测试案例。

 <ul>   <li>  在2d中，让我们选择一个半圆作为域。在这个域上，我们选择函数 $u(\mathbf x)=-2x_1x_2$ 作为解决方案。为了计算右手边，我们必须计算解函数的表面拉普拉斯。有（至少）两种方法可以做到这一点。第一种是使用 $u(\mathbf x)$ 的自然延伸（仍然用 $u$ 表示）在 $\mathbb R^d$ 上投影掉上面描述的法向导数，即计算@f[


    -\Delta_\Gamma u =  \Delta u - \mathbf n^T \ D^2 u \ \mathbf n - (\mathbf n \cdot \nabla u)\ \kappa,
  @f] 。

  其中  $\kappa$  是  $\Gamma$  的总曲率。   由于我们在单位圆上， $\mathbf n=\mathbf x$ 和 $\kappa = 1$ 所以@f[


    -\Delta_\Gamma u = -8 x_1x_2.
  @f]



  一个更简单的方法，至少对于目前二维空间的曲线的情况，是注意到我们可以用变换 $t \in
  [0,\pi]$ 将区间 $\Omega$ 映射到域 $\mathbf x(t)= \left(\begin{array}{c} \cos t \\ \sin t \end{array}\right)$ 。   在位置  $\mathbf x=\mathbf x(t)$  上，解的值是  $u(\mathbf x(t)) = -2\cos t \sin t$  。   考虑到转换是保长的，即长度为 $dt$ 的线段被映射到完全相同长度的曲线上，那么切向拉普拉斯就满足@f{align*}
    \Delta_\Gamma u
    &= \frac{d^2}{dt^2}(-2\cos t \sin t)
    = -2 \frac{d}{dt}(-\sin^2 t + \cos^2 t)
    = -2 (-2 \sin t \cos t - 2 \cos t \sin t)
    \\
    &= 8 \sin t \cos t
    \\
    &= 8 x_1x_2,
  @f} 。

  这当然和我们上面的结果是一样的。   </li>   <li>  在三维中，域又是单位球表面的一半，即半球或圆顶。我们选择 $u(\mathbf x)=-2\sin(\pi x_1)\cos(\pi x_2)e^z$ 作为解决方案。我们可以用上面的方法计算方程的右边， $f=-\Delta_\Gamma u$ ，（用 $\kappa = 2$ ），得到一个笨拙而冗长的表达。你可以在源代码中找到完整的表达式。   </li>   </ul> 。

在程序中，我们还将计算出解的 $H^1$ 半规范误差。由于解函数及其数值近似只在流形上定义，这个误差函数的明显定义是 $| e |_{H^1(\Gamma)}
  = | \nabla_\Gamma e |_{L_2(\Gamma)}
  = \left( \int_\Gamma | \nabla_\Gamma (u-u_h) |^2 \right)^{1/2}$  。这就要求我们为函数 VectorTools::integrate_difference （在步骤7中首次引入）提供<i>tangential</i>梯度 $\nabla_\Gamma u$ ，我们将通过在下面的程序中实现函数 <code>Solution::gradient</code> 来实现。




<a name="Implementation"></a><h3>Implementation</h3>


如果你已经读完了第4步，并且理解了上面关于解和右边如何对应的讨论，你也会立即熟悉这个程序。事实上，只有两件事是有意义的。

- 我们生成三角计算域的网格的方式。

- 我们使用映射对象的方式来描述，我们解决偏微分方程的领域不是平面的，实际上是弯曲的。

在第10步和第11步中已经介绍了映射对象，正如那里所解释的，只要你对边界的样子有一个有效的描述，你通常不需要知道它们是如何工作的。从本质上讲，我们将简单地声明一个适当的MappingQ类型的对象，它将自动从三角图中获得边界描述。然后，该映射对象将被传递给适当的函数，我们将得到库中预定义的半圆或半球形的边界描述。

该程序的其余部分紧跟步骤4，至于计算误差，则是步骤7。这个程序的某些方面，特别是在Triangulation、DoFHandler和类似的类上使用两个模板参数，已经在步骤34中作了详细描述；你不妨也读一读这个教程程序。


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
 * 如果你读过 step-4 和 step-7 ，你会认识到我们已经在那里使用了以下所有的包含文件。因此，我们不会在这里再次解释它们的含义。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/solver_control.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_q.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step38 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeLaplaceBeltramiProblemcodeclasstemplate"></a> 
 * <h3>The <code>LaplaceBeltramiProblem</code> class template</h3>
 * 

 * 
 * 这个类几乎与 step-4 中的 <code>LaplaceProblem</code> 类完全相似。
 * 

 * 
 * 本质上的区别是这样的。
 * 

 * 
 * 

 * 
 * 

 * 
 * - 模板参数现在表示嵌入空间的维度，它不再与域和我们计算的三角形的维度相同。我们通过调用参数 @p spacedim, 并引入一个等于域的维度的常数 @p dim 来表明这一点--这里等于 <code>spacedim-1</code>  。
 * 

 * 
 * - 所有具有几何特征的成员变量现在都需要知道它们自己的维度以及嵌入空间的维度。因此，我们需要指定它们的模板参数，一个是网格的维度 @p dim, ，另一个是嵌入空间的维度， @p spacedim.  这正是我们在 step-34 中所做的，请看那里有更深的解释。
 * 

 * 
 * - 我们需要一个对象来描述从参考单元到三角形组成的单元所使用的哪种映射。从Mapping基类派生出来的类正是这样做的。在deal.II的大部分时间里，如果你不做任何事情，图书馆会假定你想要一个使用（双，三）线性映射的MappingQ1对象。在许多情况下，这就足够了，这就是为什么这些对象的使用大多是可选的：例如，如果你有一个二维空间中的多边形二维域，参考单元到三角形单元的双线性映射会产生该域的精确表示。如果你有一个弯曲的域，你可能想对那些位于域的边界的单元使用一个高阶映射--例如，这就是我们在 step-11 中所做的。然而，在这里我们有一个弯曲的域，而不仅仅是一个弯曲的边界，虽然我们可以用双线性映射的单元来近似它，但对所有单元使用高阶映射才是真正谨慎的。因此，这个类有一个MappingQ类型的成员变量；我们将选择映射的多项式程度等于计算中使用的有限元的多项式程度，以确保最佳近似，尽管这种等参数性不是必须的。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   class LaplaceBeltramiProblem 
 *   { 
 *   public: 
 *     LaplaceBeltramiProblem(const unsigned degree = 2); 
 *     void run(); 
 * 
 *   private: 
 *     static constexpr unsigned int dim = spacedim - 1; 
 * 
 *     void make_grid_and_dofs(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void output_results() const; 
 *     void compute_error() const; 
 * 
 *     Triangulation<dim, spacedim> triangulation; 
 *     FE_Q<dim, spacedim>          fe; 
 *     DoFHandler<dim, spacedim>    dof_handler; 
 *     MappingQ<dim, spacedim>      mapping; 
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
 * <a name="Equationdata"></a> 
 * <h3>Equation data</h3>
 * 

 * 
 * 接下来，让我们定义描述问题的精确解和右手边的类。这与 step-4 和 step-7 相类似，在那里我们也定义了此类对象。鉴于介绍中的讨论，实际的公式应该是不言自明的。值得关注的一点是，我们是如何使用一般模板的明确特化，分别定义2D和3D情况下的值和梯度函数的。另一种方法是定义通用模板，并为空间维度的每个可能的值设置一个 <code>switch</code> 语句（或一串 <code>if</code> s）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solution : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 * 
 *     virtual Tensor<1, dim> 
 *     gradient(const Point<dim> & p, 
 *              const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <> 
 *   double Solution<2>::value(const Point<2> &p, const unsigned int) const 
 *   { 
 *     return (-2. * p(0) * p(1)); 
 *   } 
 * 
 *   template <> 
 *   Tensor<1, 2> Solution<2>::gradient(const Point<2> &p, 
 *                                      const unsigned int) const 
 *   { 
 *     Tensor<1, 2> return_value; 
 *     return_value[0] = -2. * p(1) * (1 - 2. * p(0) * p(0)); 
 *     return_value[1] = -2. * p(0) * (1 - 2. * p(1) * p(1)); 
 * 
 *     return return_value; 
 *   } 
 * 
 *   template <> 
 *   double Solution<3>::value(const Point<3> &p, const unsigned int) const 
 *   { 
 *     return (std::sin(numbers::PI * p(0)) * std::cos(numbers::PI * p(1)) * 
 *             exp(p(2))); 
 *   } 
 * 
 *   template <> 
 *   Tensor<1, 3> Solution<3>::gradient(const Point<3> &p, 
 *                                      const unsigned int) const 
 *   { 
 *     using numbers::PI; 
 * 
 *     Tensor<1, 3> return_value; 
 * 
 *     return_value[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 *     return_value[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 *     return_value[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 * 
 *     return return_value; 
 *   } 
 * 
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <> 
 *   double RightHandSide<2>::value(const Point<2> &p, 
 *                                  const unsigned int /*component*/) const 
 *   { 
 *     return (-8. * p(0) * p(1)); 
 *   } 
 * 
 *   template <> 
 *   double RightHandSide<3>::value(const Point<3> &p, 
 *                                  const unsigned int /*component*/) const 
 *   { 
 *     using numbers::PI; 
 * 
 *     Tensor<2, 3> hessian; 
 * 
 *     hessian[0][0] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 *     hessian[1][1] = -PI * PI * sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 *     hessian[2][2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 * 
 *     hessian[0][1] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 *     hessian[1][0] = -PI * PI * cos(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 * 
 *     hessian[0][2] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 *     hessian[2][0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 * 
 *     hessian[1][2] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 *     hessian[2][1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 * 
 *     Tensor<1, 3> gradient; 
 *     gradient[0] = PI * cos(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 *     gradient[1] = -PI * sin(PI * p(0)) * sin(PI * p(1)) * exp(p(2)); 
 *     gradient[2] = sin(PI * p(0)) * cos(PI * p(1)) * exp(p(2)); 
 * 
 *     Point<3> normal = p; 
 *     normal /= p.norm(); 
 * 
 *     return (-trace(hessian) + 2 * (gradient * normal) + 
 *             (hessian * normal) * normal); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeLaplaceBeltramiProblemcodeclass"></a> 
 * <h3>Implementation of the <code>LaplaceBeltramiProblem</code> class</h3>
 * 

 * 
 * 如果你知道  step-4  ，程序的其余部分实际上是很不引人注目的。我们的第一步是定义构造函数，设置有限元和映射的多项式程度，并将DoF处理程序与三角形关联。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   LaplaceBeltramiProblem<spacedim>::LaplaceBeltramiProblem( 
 *     const unsigned degree) 
 *     : fe(degree) 
 *     , dof_handler(triangulation) 
 *     , mapping(degree) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemmake_grid_and_dofs"></a> 
 * <h4>LaplaceBeltramiProblem::make_grid_and_dofs</h4>
 * 

 * 
 * 下一步是创建网格，分配自由度，并设置描述线性系统的各种变量。所有这些步骤都是标准的，只有如何创建一个描述曲面的网格除外。我们可以为我们感兴趣的领域生成一个网格，用一个网格生成器生成一个三角形，然后用GridIn类将其读入。或者，就像我们在这里做的那样，我们使用GridGenerator命名空间的设施来生成网格。
 * 

 * 
 * 具体来说，我们要做的是这样的（在下面的大括号中）：我们使用 <code>spacedim</code> 函数为半圆盘（2D）或半球（3D）生成一个 GridGenerator::half_hyper_ball 维度的网格。这个函数将位于圆盘/球周边的所有面的边界指标设置为零，而在将整个圆盘/球分成两半的直线部分设置为零。下一步是主要的一点。 GridGenerator::extract_boundary_mesh 函数创建的网格是由那些作为前一个网格的面的单元组成的，也就是说，它描述了原始（体积）网格的<i>surface</i>单元。然而，我们不需要所有的面：只需要那些在圆盘或球的周边，边界指示器为零的面；我们可以使用一组边界指示器来选择这些单元，并传递给 GridGenerator::extract_boundary_mesh. 。
 * 

 * 
 * 有一点需要提及。为了在流形是弯曲的情况下适当地细化表面网格（类似于细化与弯曲边界相邻的单元面），三角形必须要有一个对象附加在上面，描述新顶点应该位于何处。如果你不附加这样的边界对象，它们将位于现有顶点之间的中间位置；如果你有一个具有直线边界的域（例如多边形），这是很合适的，但如果像这里一样，流形具有曲率，则不合适。因此，为了让事情正常进行，我们需要将流形对象附加到我们的（表面）三角形上，其方式与我们在1d中为边界所做的大致相同。我们创建这样一个对象，并将其附加到三角剖面上。
 * 

 * 
 * 创建网格的最后一步是对其进行多次细化。该函数的其余部分与之前的教程程序中相同。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::make_grid_and_dofs() 
 *   { 
 *     { 
 *       Triangulation<spacedim> volume_mesh; 
 *       GridGenerator::half_hyper_ball(volume_mesh); 
 * 
 *       std::set<types::boundary_id> boundary_ids; 
 *       boundary_ids.insert(0); 
 * 
 *       GridGenerator::extract_boundary_mesh(volume_mesh, 
 *                                            triangulation, 
 *                                            boundary_ids); 
 *     } 
 *     triangulation.set_all_manifold_ids(0); 
 *     triangulation.set_manifold(0, SphericalManifold<dim, spacedim>()); 
 * 
 *     triangulation.refine_global(4); 
 * 
 *     std::cout << "Surface mesh has " << triangulation.n_active_cells() 
 *               << " cells." << std::endl; 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::cout << "Surface mesh has " << dof_handler.n_dofs() 
 *               << " degrees of freedom." << std::endl; 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemassemble_system"></a> 
 * <h4>LaplaceBeltramiProblem::assemble_system</h4>
 * 

 * 
 * 下面是这个程序的中心函数，即组装与表面拉普拉斯（Laplace-Beltrami算子）相对应的矩阵。也许令人惊讶的是，它实际上与例如在  step-4  中讨论的普通拉普拉斯算子看起来完全一样。关键是 FEValues::shape_grad() 函数发挥了魔力：它返回 $i$ 第1个形状函数在 $q$ 第1个正交点的表面梯度 $\nabla_K \phi_i(x_q)$ 。其余的也不需要任何改变。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::assemble_system() 
 *   { 
 *     system_matrix = 0; 
 *     system_rhs    = 0; 
 * 
 *     const QGauss<dim>       quadrature_formula(2 * fe.degree); 
 *     FEValues<dim, spacedim> fe_values(mapping, 
 *                                       fe, 
 *                                       quadrature_formula, 
 *                                       update_values | update_gradients | 
 *                                         update_quadrature_points | 
 *                                         update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<double>                  rhs_values(n_q_points); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     RightHandSide<spacedim> rhs; 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix = 0; 
 *         cell_rhs    = 0; 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         rhs.value_list(fe_values.get_quadrature_points(), rhs_values); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *               cell_matrix(i, j) += fe_values.shape_grad(i, q_point) * 
 *                                    fe_values.shape_grad(j, q_point) * 
 *                                    fe_values.JxW(q_point); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             cell_rhs(i) += fe_values.shape_value(i, q_point) * 
 *                            rhs_values[q_point] * fe_values.JxW(q_point); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           { 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               system_matrix.add(local_dof_indices[i], 
 *                                 local_dof_indices[j], 
 *                                 cell_matrix(i, j)); 
 * 
 *             system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *           } 
 *       } 
 * 
 *     std::map<types::global_dof_index, double> boundary_values; 
 *     VectorTools::interpolate_boundary_values( 
 *       mapping, dof_handler, 0, Solution<spacedim>(), boundary_values); 
 * 
 *     MatrixTools::apply_boundary_values( 
 *       boundary_values, system_matrix, solution, system_rhs, false); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemsolve"></a> 
 * <h4>LaplaceBeltramiProblem::solve</h4>
 * 

 * 
 * 下一个函数是解决线性系统的函数。在这里，也不需要做任何改变。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::solve() 
 *   { 
 *     SolverControl solver_control(solution.size(), 1e-7 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution, system_rhs, preconditioner); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemoutput_result"></a> 
 * <h4>LaplaceBeltramiProblem::output_result</h4>
 * 

 * 
 * 这是一个从解决方案中生成图形输出的函数。它的大部分都是模板代码，但有两点值得指出。
 * 

 * 
 * 

 * 
 * 

 * 
 * -  DataOut::add_data_vector() 函数可以接受两种向量。  一种是之前通过 DataOut::attach_dof_handler(); 连接的DoFHandler对象定义的每个自由度有一个值的向量，另一种是三角测量的每个单元有一个值的向量，例如，输出每个单元的估计误差。通常，DataOut类知道如何区分这两种向量：自由度几乎总是比单元格多，所以我们可以通过两种向量的长度来区分。我们在这里也可以这样做，但只是因为我们很幸运：我们使用了一个半球体。如果我们用整个球体作为域和 $Q_1$ 元素，我们将有相同数量的单元格作为顶点，因此这两种向量将有相同数量的元素。为了避免由此产生的混乱，我们必须告诉 DataOut::add_data_vector() 函数我们有哪种矢量。DoF数据。这就是该函数的第三个参数的作用。
 * 

 * 
 * -  DataOut::build_patches() 函数可以生成细分每个单元的输出，这样可视化程序可以更好地解决弯曲流形或更高的多项式程度的形状函数。在这里，我们在每个坐标方向上对每个元素进行细分，细分的次数与使用的有限元的多项式程度相同。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::output_results() const 
 *   { 
 *     DataOut<dim, spacedim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, 
 *                              "solution", 
 *                              DataOut<dim, spacedim>::type_dof_data); 
 *     data_out.build_patches(mapping, mapping.get_degree()); 
 * 
 *     const std::string filename = 
 *       "solution-" + std::to_string(spacedim) + "d.vtk"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtk(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemcompute_error"></a> 
 * <h4>LaplaceBeltramiProblem::compute_error</h4>
 * 

 * 
 * 这是最后一块功能：我们要计算数值解的误差。它是之前在  step-7  中展示和讨论的代码的逐字复制。正如介绍中提到的， <code>Solution</code> 类提供了解决方案的（切向）梯度。为了避免只评估超收敛点的误差，我们选择一个足够高阶的正交规则。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::compute_error() const 
 *   { 
 *     Vector<float> difference_per_cell(triangulation.n_active_cells()); 
 *     VectorTools::integrate_difference(mapping, 
 *                                       dof_handler, 
 *                                       solution, 
 *                                       Solution<spacedim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(2 * fe.degree + 1), 
 *                                       VectorTools::H1_norm); 
 * 
 *     double h1_error = VectorTools::compute_global_error(triangulation, 
 *                                                         difference_per_cell, 
 *                                                         VectorTools::H1_norm); 
 *     std::cout << "H1 error = " << h1_error << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="LaplaceBeltramiProblemrun"></a> 
 * <h4>LaplaceBeltramiProblem::run</h4>
 * 

 * 
 * 最后一个函数提供了顶层的逻辑。它的内容是不言自明的。
 * 

 * 
 * 
 * @code
 *   template <int spacedim> 
 *   void LaplaceBeltramiProblem<spacedim>::run() 
 *   { 
 *     make_grid_and_dofs(); 
 *     assemble_system(); 
 *     solve(); 
 *     output_results(); 
 *     compute_error(); 
 *   } 
 * } // namespace Step38 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 该程序的其余部分由 <code>main()</code> 函数占据。它完全遵循首次在 step-6 中介绍的一般布局，并在随后的所有教程程序中使用。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step38; 
 * 
 *       LaplaceBeltramiProblem<3> laplace_beltrami; 
 *       laplace_beltrami.run(); 
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
examples/step-38/doc/results.dox



<a name="Results"></a><h1>Results</h1>


当你运行该程序时，应在屏幕上打印出以下输出。

@verbatim
Surface mesh has 1280 cells.
Surface mesh has 5185 degrees of freedom.
H1 error = 0.0217136
@endverbatim




通过在 <code>LaplaceBeltrami::make_grid_and_dofs</code> 函数中玩弄全局细化的数量，可以增加或减少网格细化。例如，多做一次细化，只运行三维曲面问题，得到的输出结果如下。

@verbatim
Surface mesh has 5120 cells.
Surface mesh has 20609 degrees of freedom.
H1 error = 0.00543481
@endverbatim



这就是我们所期望的：将网格尺寸缩小2倍，误差下降4倍（记住我们使用的是双二次元）。从一到五次细化的全部误差序列看起来是这样的，整齐地遵循理论上预测的模式。

@verbatim
0.339438
0.0864385
0.0217136
0.00543481
0.00135913
@endverbatim



最后，该程序产生图形输出，我们可以将其可视化。下面是一个结果图。

 <img src="https://www.dealii.org/images/steps/developer/step-38.solution-3d.png" alt=""> 

该程序也适用于2D中的1D曲线，而不仅仅是3D中的2D曲面。你可以通过改变 <code>main()</code> 中的模板参数来测试这一点，像这样。

@code
      LaplaceBeltramiProblem<2> laplace_beltrami;
@endcode

域是一条2D的曲线，我们可以通过使用第三维（和颜色）来表示函数 $u(x)$ 的值来可视化解决方案。这样看起来就像这样（白色的曲线是域，彩色的曲线是被挤压到第三维的解决方案，清楚地显示了当曲线从域的一个象限移动到相邻的象限时符号的变化）。

 <img src="https://www.dealii.org/images/steps/developer/step-38.solution-2d.png" alt=""> 


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


只有当表面比仅仅是一个半球体更有趣时，在表面上的计算才会变得有趣。为了达到这个目的，deal.II可以通过通常的GridIn类读取描述曲面的网格。或者，万一你有一个分析性的描述，一个简单的网格有时可以被拉伸和弯曲成我们感兴趣的形状。

让我们考虑一个相对简单的例子：我们把之前用过的半球体，在Z方向上拉伸10倍，然后把X和Y坐标拼一下。让我们先展示一下计算域和解决方案，然后再讨论下面的实现细节。

 <img src="https://www.dealii.org/images/steps/developer/step-38.warp-1.png" alt=""> 

 <img src="https://www.dealii.org/images/steps/developer/step-38.warp-2.png" alt=""> 

产生这种网格的方法是使用 GridTools::transform() 函数。它需要一个方法来转换每个单独的网格点到不同的位置。让我们在这里使用下面这个相当简单的函数（记住：在一个方向上拉伸，在另外两个方向上拼凑）。

@code
template <int spacedim>
Point<spacedim> warp(const Point<spacedim> &p)
{
  Point<spacedim> q = p;
  q[spacedim-1] *= 10;


  if (spacedim >= 2)
    q[0] += 2*std::sin(q[spacedim-1]);
  if (spacedim >= 3)
    q[1] += 2*std::cos(q[spacedim-1]);


  return q;
}
@endcode



如果我们遵循 <code>LaplaceBeltrami::make_grid_and_dofs</code> 函数，我们会像以前一样提取半球形表面网格，将其扭曲成我们想要的形状，并根据需要经常进行细化。但这并不像我们所希望的那样简单：细化需要我们有一个适当的流形对象附加到三角形上，描述细化时网格的新顶点应该位于何处。我相信可以通过简单地撤销上面的变换（重新得到球面），找到球面上新的点的位置，然后重新扭曲结果，以一种不太复杂的方式描述这个流形。但我是个懒人，既然这样做并不是真正的重点，我们还是让我们的生活变得简单一点：我们将提取半球体，根据需要对其进行细化，摆脱描述流形的对象，因为我们现在不再需要它，然后最后对网格进行扭曲。使用上面的函数，这将看起来如下。

@code
template <int spacedim>
void LaplaceBeltrami<spacedim>::make_grid_and_dofs()
{
  {
    Triangulation<spacedim> volume_mesh;
    GridGenerator::half_hyper_ball(volume_mesh);


    volume_mesh.refine_global(4);


    std::set<types::boundary_id> boundary_ids;
    boundary_ids.insert(0);


    GridGenerator::extract_boundary_mesh(volume_mesh, triangulation,
                                         boundary_ids);
    GridTools::transform(&warp<spacedim>, triangulation);       /* ** */
    std::ofstream x("x"), y("y");
    GridOut().write_gnuplot(volume_mesh, x);
    GridOut().write_gnuplot(triangulation, y);
  }


  std::cout << "Surface mesh has " << triangulation.n_active_cells()
            << " cells."
            << std::endl;
  ...
}
@endcode



请注意，唯一必要的补充是标有星号的那一行。不过值得指出的是：由于我们将流形描述从表面网格中分离出来，所以当我们在程序的其余部分使用映射对象时，它不再有曲线边界描述可言。相反，它将不得不使用隐含的FlatManifold类，该类用于域的所有未明确指定不同流形对象的部分。因此，无论我们使用MappingQ(2)、MappingQ(15)还是MappingQ1，我们的网格的每个单元都将使用双线性近似进行映射。

撇开所有这些缺点不谈，得到的图片还是很好看的。与步骤38中的内容唯一不同的是，我们把右手边改为 $f(\mathbf x)=\sin x_3$ ，把边界值（通过 <code>Solution</code> 类）改为 $u(\mathbf x)|_{\partial\Omega}=\cos x_3$  。当然，我们现在已经不知道确切的解决方案，所以在 <code>LaplaceBeltrami::run</code> 末尾的误差计算将得到一个毫无意义的数字。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-38.cc"
*/
