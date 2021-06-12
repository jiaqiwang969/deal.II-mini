/**
@page step_20 The step-20 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Theequations">The equations</a>
        <li><a href="#Formulationweakformanddiscreteproblem">Formulation, weak form, and discrete problem</a>
        <li><a href="#Assemblingthelinearsystem">Assembling the linear system</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
      <ul>
        <li><a href="#SolvingusingtheSchurcomplement">Solving using the Schur complement</a>
        <li><a href="#TheLinearOperatorframeworkindealII">The LinearOperator framework in deal.II</a>
        <li><a href="#ApreconditionerfortheSchurcomplement">A preconditioner for the Schur complement</a>
      </ul>
        <li><a href="#Definitionofthetestcase">Definition of the test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeMixedLaplaceProblemcodeclasstemplate">The <code>MixedLaplaceProblem</code> class template</a>
        <li><a href="#Righthandsideboundaryvaluesandexactsolution">Right hand side, boundary values, and exact solution</a>
        <li><a href="#Theinversepermeabilitytensor">The inverse permeability tensor</a>
        <li><a href="#MixedLaplaceProblemclassimplementation">MixedLaplaceProblem class implementation</a>
      <ul>
        <li><a href="#MixedLaplaceProblemMixedLaplaceProblem">MixedLaplaceProblem::MixedLaplaceProblem</a>
        <li><a href="#MixedLaplaceProblemmake_grid_and_dofs">MixedLaplaceProblem::make_grid_and_dofs</a>
        <li><a href="#MixedLaplaceProblemassemble_system">MixedLaplaceProblem::assemble_system</a>
      </ul>
        <li><a href="#Implementationoflinearsolversandpreconditioners">Implementation of linear solvers and preconditioners</a>
      <ul>
        <li><a href="#MixedLaplacesolve">MixedLaplace::solve</a>
      </ul>
        <li><a href="#MixedLaplaceProblemclassimplementationcontinued">MixedLaplaceProblem class implementation (continued)</a>
      <ul>
        <li><a href="#MixedLaplacecompute_errors">MixedLaplace::compute_errors</a>
        <li><a href="#MixedLaplaceoutput_results">MixedLaplace::output_results</a>
        <li><a href="#MixedLaplacerun">MixedLaplace::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Outputoftheprogramandgraphicalvisualization">Output of the program and graphical visualization</a>
        <li><a href="#Convergence">Convergence</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Morerealisticpermeabilityfields">More realistic permeability fields</a>
        <li><a href="#Betterlinearsolvers">Better linear solvers</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-20/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{19,20,21} 

这个程序致力于两个方面：使用混合有限元--特别是Raviart-Thomas元--以及使用块状矩阵来定义求解器、预处理器和使用系统矩阵的子结构的嵌套版本。我们要解决的方程仍然是泊松方程，虽然有一个矩阵值的系数。

@f{eqnarray*}


  -\nabla \cdot K({\mathbf x}) \nabla p &=& f \qquad {\textrm{in}\ } \Omega, \\
  p &=& g \qquad {\textrm{on}\ }\partial\Omega.


@f}

 $K({\mathbf x})$ 被假定为均匀正定，即有 $\alpha>0$ ，使得 $K(x)$ 的特征值 $\lambda_i({\mathbf x})$ 满足 $\lambda_i({\mathbf x})\ge \alpha$  。使用符号 $p$ 而不是通常的 $u$ 作为解变量将在下一节中变得清晰。

在讨论了方程和我们要用来解决它的公式之后，这个介绍将包括块状矩阵和向量的使用，求解器和预处理器的定义，以及最后我们要解决的实际测试案例。

我们将在第21步中扩展这个教程程序，不仅要解决混合拉普拉斯方程，还要增加另一个描述两种流体混合物运输的方程。

这里所涉及的方程属于矢量值问题的范畴。这个主题的顶层概述可以在 @ref vector_valued 模块中找到。




<a name="Theequations"></a><h3>The equations</h3>


在上述形式中，泊松方程（即具有非零右手边的拉普拉斯方程）通常被认为是流体在多孔介质中流动的良好模型方程。当然，人们通常通过<a href="https://en.wikipedia.org/wiki/Navier%E2%80%93Stokes_equations">Navier-Stokes
equations</a>来模拟流体流动，或者，如果流体速度很慢或粘度很大，则通过<a href="https://en.wikipedia.org/wiki/Stokes_flow">Stokes
equations</a>（我们在步骤22中涉及）。在这两个模型中，第一个模型的作用力是惯性和粘性摩擦力，而在第二个模型中，只有粘性摩擦力--即一个流体粒子对附近的粒子施加的力。如果你在一个大的领域里有自由流动，例如管道、河流或空气中，这是很合适的。另一方面，如果流体被限制在孔隙中，那么孔壁对流体施加的摩擦力变得越来越重要，而内部粘性摩擦力变得越来越不重要。如果这两种效应都很重要，那么建立模型首先会导致<a href="https://en.wikipedia.org/wiki/Darcy%27s_law#Brinkman_form_of_Darcy's_law">Brinkman
model</a>，而在非常小的孔隙的限制下会导致<a href="https://en.wikipedia.org/wiki/Darcy%27s_law">Darcy equations</a>。后者只是泊松方程或拉普拉斯方程的不同名称，其内涵是人们想应用它的领域：多孔介质中的缓慢流动。本质上，它说速度与驱动流体通过多孔介质的负压梯度成正比。

达西方程对驱动流动的这种压力进行建模。由于解变量是压力，我们在这里使用 $p$ 这个名字，而不是更常用于偏微分方程解的 $u$ 这个名字）。这种拉普拉斯方程观点的典型应用是为地下水流建模，或者为油藏中的碳氢化合物流动建模。在这些应用中， $K$ 是渗透性张量，即衡量土壤或岩石基质对流体流动的阻力大小。

在上述应用中，数值方案的一个理想特征是它应该是局部保守的，也就是说，无论什么东西流入一个单元，也会从该单元流出（或者如果源不为零，则差值等于每个单元的源项的积分）。然而，事实证明，拉普拉斯方程的通常离散化（如步骤3、步骤4或步骤6中使用的那些）并不满足这一特性。但是，人们可以通过选择问题的不同表述和有限元空间的特定组合来实现这一点。




<a name="Formulationweakformanddiscreteproblem"></a><h3>Formulation, weak form, and discrete problem</h3>


为此，我们首先引入了第二个变量，称为速度，  ${\mathbf u}=-K\nabla p$  。根据其定义，速度是压力梯度的负方向的一个矢量，乘以渗透性张量。如果渗透率张量与单位矩阵成正比，这个方程就很容易理解和直观：渗透率越高，速度越高；速度与压力梯度成正比，从高压区到低压区（因此是负号）。

有了这第二个变量，就可以找到拉普拉斯方程的另一个版本，称为<i>mixed formulation</i>。

@f{eqnarray*}
  K^{-1} {\mathbf u} + \nabla p &=& 0 \qquad {\textrm{in}\ } \Omega, \\


  -{\textrm{div}}\ {\mathbf u} &=& -f \qquad {\textrm{in}\ }\Omega, \\
  p &=& g \qquad {\textrm{on}\ } \partial\Omega.


@f}

这里，我们将定义速度的方程 ${\mathbf
u}$ 乘以 $K^{-1}$ ，因为这使得方程组是对称的：其中一个方程有梯度，第二个方程有负发散，这两个当然是彼此相邻的，结果是一个对称的双线性形式，因此在 $K$ 是对称张量的共同假设下，是一个对称的系统矩阵。

这个问题的弱表述是通过将两个方程与测试函数相乘，并对一些项进行分项积分来找到的。

@f{eqnarray*}
  A(\{{\mathbf u},p\},\{{\mathbf v},q\}) = F(\{{\mathbf v},q\}),


@f}

其中

@f{eqnarray*}
  A(\{{\mathbf u},p\},\{{\mathbf v},q\})
  &=&
  ({\mathbf v}, K^{-1}{\mathbf u})_\Omega - ({\textrm{div}}\ {\mathbf v}, p)_\Omega


  - (q,{\textrm{div}}\ {\mathbf u})_\Omega
  \\
  F(\{{\mathbf v},q\}) &=& -(g,{\mathbf v}\cdot {\mathbf n})_{\partial\Omega} - (f,q)_\Omega.


@f}

这里， ${\mathbf n}$ 是边界处的外向法向量。请注意，在这个公式中，原问题的迪里希特边界值被纳入到弱形式中。

为了得到良好的解决，我们必须在空间 $H({\textrm{div}})=\{{\mathbf w}\in L^2(\Omega)^d:\ {\textrm{div}}\ {\mathbf w}\in L^2\}$ 中寻找 $\mathbf u$ 、 $\mathbf v$ 和 $L^2$ 的解和检验函数。几乎每一本关于有限元理论的书都说过一个众所周知的事实，如果选择离散的有限元空间来逼近 ${\mathbf u},p$ 是不恰当的，那么产生的离散问题是不稳定的，离散的解将不会收敛到精确的解。这里考虑的问题的一些细节--属于 "鞍点问题 "的范畴

--可以在维基百科上找到<a
href="https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition">Ladyzhenskaya-Babuska-Brezzi
(LBB) condition</a>的页面）。)

为了克服这个问题，已经为 ${\mathbf u},p$ 开发了一些不同的有限元对，导致了稳定的离散问题。其中一个对子是对速度 ${\mathbf u}$ 使用Raviart-Thomas空间，对压力 $p$ 使用类 $DQ(k)$ 的不连续元素。关于这些空间的细节，我们特别参考Brezzi和Fortin的关于混合有限元方法的书，但许多其他关于有限元理论的书，例如Brenner和Scott的经典书，也说明了相关结果。在任何情况下，在适当选择函数空间的情况下，离散的表述如下。找到 ${\mathbf
u}_h,p_h$ ，以便

@f{eqnarray*}
  A(\{{\mathbf u}_h,p_h\},\{{\mathbf v}_h,q_h\}) = F(\{{\mathbf v}_h,q_h\})
  \qquad\qquad \forall {\mathbf v}_h,q_h.


@f}




在继续之前，让我们简单地停顿一下，说明上面的函数空间的选择为我们提供了所需的局部守恒特性。特别是，由于压力空间由不连续的片断多项式组成，我们可以选择测试函数 $q$ 作为在任何给定单元 $K$ 上等于1，其他地方为0的函数。如果我们也到处选择 ${\mathbf v}=0$ （记住，上面的弱式对<i>all</i>离散测试函数 $q,v$ 必须成立），那么把这些测试函数的选择放入上面的弱式中，特别意味着

@f{eqnarray*}


  - (1,{\textrm{div}}\ {\mathbf u}_h)_K
  =


  -(1,f)_K,


@f}

当然，我们可以将其以更明确的形式写为

@f{eqnarray*}
  \int_K {\textrm{div}}\ {\mathbf u}_h
 =
  \int_K f.


@f}

应用发散定理的结果是，对于每个单元 ${\mathbf
u}_h$ 的选择， $K$ 必须满足以下关系

@f{eqnarray*}
  \int_{\partial K} {\mathbf u}_h\cdot{\mathbf n}
  =
  \int_K f.


@f}

如果你现在记得 ${\mathbf u}$ 是速度，那么左边的积分正好是穿过单元格边界的（离散）通量 $K$  。那么声明是，通量必须等于对 $K$ 内的源的积分。特别是，如果没有源（即 $f=0$ 在 $K$ 中），那么声明是<i>total</i>通量为零，也就是说，任何流入一个单元的东西都必须通过单元边界的其他部分流出。这就是我们所说的<i>local conservation</i>，因为它对每个细胞都是成立的。

另一方面，通常的连续 $Q_k$ 元素在用于压力时不会产生这种性质（例如，我们在步骤-43中所做的），因为我们不能选择一个离散的测试函数 $q_h$ ，在单元 $K$ 上为1，其他地方为0：它将是不连续的，因此不在有限元空间内。严格来说，我们只能说上面的证明对连续元素不起作用。这些元素是否仍然可能导致局部守恒是一个不同的问题，因为人们可以认为不同的证明可能仍然有效；然而，在现实中，这个属性确实不成立）。)




<a name="Assemblingthelinearsystem"></a><h3>Assembling the linear system</h3>


deal.II库（当然）实现了任意阶的Raviart-Thomas元素 $RT(k)$ ，以及不连续元素 $DG(k)$  。如果我们暂时忘记它们的特殊属性，那么我们就必须解决一个离散的问题

@f{eqnarray*}
  A(x_h,w_h) = F(w_h),


@f}

的双线性形式和右手边，以及 $x_h=\{{\mathbf u}_h,p_h\}$  ,  $w_h=\{{\mathbf v}_h,q_h\}$  。 $x_h$ 和 $w_h$ 都来自空间 $X_h=RT(k)\times DQ(k)$ ，其中 $RT(k)$ 本身就是一个 $dim$ 维函数的空间，以适应流速为矢量值的事实。那么必要的问题是：我们如何在程序中做到这一点？

矢量值元素已经在以前的教程程序中讨论过了，第一次是在步骤8中详细讨论。那里的主要区别是，矢量值空间 $V_h$ 的所有分量都是统一的：位移矢量的 $dim$ 分量都是相等的，来自同一个函数空间。因此，我们可以做的是将 $V_h$ 建立为 $dim$ 乘以通常的 $Q(1)$ 有限元空间的外积，并以此确保我们所有的形状函数只有一个非零矢量分量。因此，我们在步骤8中所做的不是处理矢量值的形状函数，而是查看（标量）唯一的非零分量，并使用 <code>fe.system_to_component_index(i).first</code> 调用来计算这实际上是哪个分量。

这对Raviart-Thomas元素不起作用：由于它们的构造满足空间 $H({\textrm{div}})$ 的某些规则性属性， $RT(k)$ 的形状函数通常在其所有矢量分量中都是不为零。由于这个原因，如果应用 <code>fe.system_to_component_index(i).first</code> 来确定形状函数 $i$ 的唯一非零分量，就会产生一个例外。我们真正需要做的是在 <em> 中获得一个形状函数的所有 </em> 向量分量。在deal.II的字典中，我们称这样的有限元为 <em> 非原始 </em> ，而那些标量的有限元或者每个矢量值的形状函数只在一个矢量分量中不为零的有限元被称为 <em>  原始 </em> 。

那么，对于非原始元素，我们要怎么做呢？为了弄清楚这个问题，让我们回到教程程序中，几乎是最开始的时候。在那里，我们了解到我们使用 <code>FEValues</code> 类来确定正交点的形状函数的值和梯度。例如，我们会调用  <code>fe_values.shape_value(i,q_point)</code>  来获得  <code>i</code>  第三个形状函数在编号为  <code>q_point</code>  的正交点的值。后来，在step-8和其他教程程序中，我们了解到这个函数调用也适用于矢量值的形状函数（原始有限元），它返回形状函数  <code>i</code>  在正交点  <code>q_point</code>  的唯一非零分量的值。

对于非原始形状函数，这显然是行不通的：形状函数 <code>i</code> 没有单一的非零向量分量，因此调用 <code>fe_values.shape_value(i,q_point)</code> 也就没有什么意义。然而，deal.II提供了第二个函数调用， <code>fe_values.shape_value_component(i,q_point,comp)</code> ，返回正交点 <code>comp</code>th vector component of shape function  <code>i</code> 的值 <code>q_point</code>, where <code>comp</code> 是一个介于零和当前有限元的矢量分量数量之间的索引；例如，我们将用于描述速度和压力的元素将有 $dim+1$ 分量。值得注意的是，这个函数调用也可用于原始形状函数：它将简单地对除一个分量外的所有分量返回零；对于非原始形状函数，它一般会对不止一个分量返回非零值。

我们现在可以尝试用矢量分量来重写上面的双线性形式。例如，在2d中，第一项可以这样改写（注意， $u_0=x_0, u_1=x_1, p=x_2$ ）。

@f{eqnarray*}
  ({\mathbf u}_h^i, K^{-1}{\mathbf u}_h^j)
  =
  &\left((x_h^i)_0, K^{-1}_{00} (x_h^j)_0\right) +
   \left((x_h^i)_0, K^{-1}_{01} (x_h^j)_1\right) + \\
  &\left((x_h^i)_1, K^{-1}_{10} (x_h^j)_0\right) +
   \left((x_h^i)_1, K^{-1}_{11} (x_h^j)_1\right).


@f}

如果我们实现了这一点，我们会得到这样的代码。

@code
  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        local_matrix(i,j) += (k_inverse_values[q][0][0] *
                              fe_values.shape_value_component(i,q,0) *
                              fe_values.shape_value_component(j,q,0)
                              +
                              k_inverse_values[q][0][1] *
                              fe_values.shape_value_component(i,q,0) *
                              fe_values.shape_value_component(j,q,1)
                              +
                              k_inverse_values[q][1][0] *
                              fe_values.shape_value_component(i,q,1) *
                              fe_values.shape_value_component(j,q,0)
                              +
                              k_inverse_values[q][1][1] *
                              fe_values.shape_value_component(i,q,1) *
                              fe_values.shape_value_component(j,q,1)
                             ) *
                             fe_values.JxW(q);
@endcode



这充其量是繁琐的，容易出错的，而且不是独立的维度。有一些明显的方法可以使事情与维度无关，但最终，代码根本不漂亮。如果我们能够简单地提取形状函数 ${\mathbf u}$ 和 $p$ 的分量，那就更好了。在程序中，我们以如下方式进行。

@code
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);


  ...


  for (unsigned int q=0; q<n_q_points; ++q)
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        local_matrix(i,j) += (fe_values[velocities].value (i, q) *
                              k_inverse_values[q] *
                              fe_values[velocities].value (j, q)


                              -
                              fe_values[velocities].divergence (i, q) *
                              fe_values[pressure].value (j, q)


                              -
                              fe_values[pressure].value (i, q) *
                              fe_values[velocities].divergence (j, q)) *
                              fe_values.JxW(q);
@endcode



事实上，这不仅是双线性形式的第一项，而且是整个事情（不包括边界贡献）。

这段代码的作用是，给定一个 <code>fe_values</code> 对象，提取形状函数 <code>i</code> 在正交点 <code>q</code> 的第一个 $dim$ 分量的值，也就是该形状函数的速度部分。换句话说，如果我们把形状函数 $x_h^i$ 写成元组 $\{{\mathbf u}_h^i,p_h^i\}$ ，那么该函数返回这个元组的速度部分。请注意，速度当然是一个 <code>dim</code> 维的张量，函数返回一个相应的对象。同样地，在我们用压力提取器下标的地方，我们提取标量压力分量。整个机制在 @ref vector_valued 模块中有更详细的描述。

在实践中，如果我们在每个最外层的循环中只评估一次形状函数、它们的梯度和发散，并存储结果，我们可以做得更好一些，因为这样可以节省一些重复的计算（通过提前计算所有相关的量，然后只在实际的循环中插入结果，甚至可以节省更多的重复操作，关于这种方法的实现见步骤22）。最后的结果是这样的，在每个空间维度上都是如此。

@code
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;


      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      k_inverse.value_list (fe_values.get_quadrature_points(),
                            k_inverse_values);


      for (unsigned int q=0; q<n_q_points; ++q)
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim> phi_i_u     = fe_values[velocities].value (i, q);
            const double        div_phi_i_u = fe_values[velocities].divergence (i, q);
            const double        phi_i_p     = fe_values[pressure].value (i, q);


            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u     = fe_values[velocities].value (j, q);
                const double        div_phi_j_u = fe_values[velocities].divergence (j, q);
                const double        phi_j_p     = fe_values[pressure].value (j, q);


                local_matrix(i,j) += (phi_i_u * k_inverse_values[q] * phi_j_u


                                      - div_phi_i_u * phi_j_p


                                      - phi_i_p * div_phi_j_u) *
                                     fe_values.JxW(q);
              }


            local_rhs(i) += -phi_i_p *
                            rhs_values[q] *
                            fe_values.JxW(q);
          }
@endcode



这非常类似于我们最初写下的双线性形式和右手边的形式。

有一个最后的项我们必须注意：右手边包含项 $(g,{\mathbf v}\cdot {\mathbf n})_{\partial\Omega}$ ，构成压力边界条件的弱执行。我们已经在步骤7中看到了如何处理面积分：本质上与域积分完全相同，只是我们必须使用FEFaceValues类而不是 <code>FEValues</code> 。为了计算边界项，我们只需在所有的边界面上进行循环并在那里进行积分。该机制的工作方式与上述相同，即提取器类也对FEFaceValues对象工作。

@code
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary())
            {
              fe_face_values.reinit(cell, face);


              pressure_boundary_values.value_list(
                fe_face_values.get_quadrature_points(), boundary_values);


              for (unsigned int q = 0; q < n_face_q_points; ++q)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  local_rhs(i) += -(fe_face_values[velocities].value(i, q) *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
@endcode



在本程序的源代码中，你会发现与上面完全相同的代码。因此，我们在下文中不做过多评论。




<a name="Linearsolversandpreconditioners"></a><h3>Linear solvers and preconditioners</h3>


在组装好线性系统后，我们就面临着解决它的任务。这里的问题是，矩阵拥有两个不理想的特性。

- 它是<a href="https://en.wikipedia.org/wiki/Definiteness_of_a_matrix">indefinite</a>，也就是说，它有正负两个特征值。   我们不想在这里证明这个属性，但要注意，对于所有形式为 $\left(\begin{array}{cc} M & B \\ B^T & 0 \end{array}\right)$ 的矩阵，如这里的 $M$ 是正定的矩阵，这都是真的。

- 矩阵的右下方有一个零块（在双线性形式中没有将压力 $p$ 与压力测试函数 $q$ 耦合的项）。

至少它是对称的，但是上面的第一个问题仍然意味着共轭梯度法是行不通的，因为它只适用于矩阵是对称和正定的问题。我们将不得不求助于其他迭代求解器，如MinRes、SymmLQ或GMRES，它们可以处理不确定的系统。然而，下一个问题立即浮现。由于零块，对角线上有零，通常的 "简单 "预处理程序（Jacobi、SSOR）都不能工作，因为它们需要除以对角线上的元素。

对于我们期望用这个程序运行的矩阵大小来说，迄今为止最简单的方法是直接使用一个直接求解器（特别是与deal.II捆绑的SparseDirectUMFPACK类）。Step-29走的就是这条路线，它表明只需3、4行代码就可以完成<i>any</i>线性系统的求解。

但是，这是一个教程。我们教的是如何做事情。因此，在下文中，我们将介绍一些可用于类似这种情况的技术。也就是说，我们将考虑线性系统不是由一个大矩阵和向量组成的，而是要将矩阵分解为<i>blocks</i>，这些矩阵对应于系统中出现的各个运算符。我们注意到，所产生的求解器并不是最优的--有更好的方法来有效地计算该系统，例如在步骤22的结果部分所解释的方法，或者我们在步骤43中用于类似于当前问题的方法。在这里，我们的目标只是介绍新的求解技术以及如何在交易中实现它们。




<a name="SolvingusingtheSchurcomplement"></a><h4>Solving using the Schur complement</h4>


鉴于使用上述标准求解器和预处理器的困难，让我们再看一下矩阵。如果我们对自由度进行排序，使所有的速度变量排在所有的压力变量之前，那么我们可以将线性系统 $Ax=b$ 细分为以下几个块。

@f{eqnarray*}
  \left(\begin{array}{cc}
    M & B \\ B^T & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),


@f}

其中 $U,P$ 分别是速度和压力自由度的值， $M$ 是速度空间上的质量矩阵， $B^T$ 对应于负发散算子， $B$ 是其转置，对应于梯度。

通过区块消除法，我们就可以按以下方式对这个系统重新排序（用系统的第一行乘以 $B^TM^{-1}$ ，然后用第二行减去）。

@f{eqnarray*}
  B^TM^{-1}B P &=& B^TM^{-1} F - G, \\
  MU &=& F - BP.


@f}

这里，矩阵 $S=B^TM^{-1}B$ （称为 $A$ 的<a href="https://en.wikipedia.org/wiki/Schur_complement">Schur complement</a>）显然是对称的，由于 $M$ 的正定性和 $B$ 具有全列秩， $S$ 也是正定的。

因此，如果我们能计算出 $S$ ，我们就可以对其应用共轭梯度法。然而，计算 $S$ 是昂贵的，因为它需要我们计算（可能很大的）矩阵 $M$ 的逆；而且 $S$ 实际上也是一个完整的矩阵，因为即使 $M$ 是稀疏的，它的逆 $M^{-1}$ 通常是一个密集的矩阵。另一方面，CG算法并不要求我们真正拥有 $S$ 的表示：只需与它形成矩阵-向量乘积即可。我们可以利用矩阵乘积是关联的这一事实，分步进行（也就是说，我们可以通过设置括号来使乘积的计算更加方便）。为了计算 $Sv=(B^TM^{-1}B)v=B^T(M^{-1}(Bv))$ ，我们<ol>  <li> 计算 $w = B v$ ； <li> 解 $My = w$ 为 $y=M^{-1}w$ ，使用CG方法应用于正定和对称质量矩阵 $M$ ； <li> 计算 $z=B^Ty$ ，得到 $z=Sv$  。   </ol>  注意我们如何从右到左评估表达式 $B^TM^{-1}Bv$ 以避免矩阵-矩阵乘积；这样，我们所要做的就是评估矩阵-向量乘积。

在下文中，我们将不得不想出表示矩阵 $S$ 的方法，以便它可以用于共轭梯度求解器，以及定义我们可以预设涉及 $S$ 的线性系统解决方案的方法，并处理与矩阵 $M$ 的线性系统求解（上述第二步骤）。

 @note  这个考虑的关键点是要认识到，为了实现CG或GMRES这样的迭代求解器，我们实际上从来不需要矩阵的实际<i>elements</i>!所需要的只是我们能够形成矩阵-向量乘积。对于预处理程序也是如此。在deal.II中，我们对这一要求进行了编码，只要求给予求解器类的矩阵和预处理器有一个 <code>vmult()</code> 成员函数来做矩阵-向量乘积。一个类如何选择实现这个函数对求解器来说并不重要。因此，类可以通过，例如，做一连串的乘积和线性求解来实现它，正如上面所讨论的。




<a name="TheLinearOperatorframeworkindealII"></a><h4>The LinearOperator framework in deal.II</h4>


deal.II包括支持以一种非常普遍的方式来描述这种线性操作。这是由LinearOperator类完成的，与 @ref ConceptMatrixType "MatrixType概念 "一样，它定义了<i>applying</i>对矢量进行线性操作的最小接口。

@code
    std::function<void(Range &, const Domain &)> vmult;
    std::function<void(Range &, const Domain &)> vmult_add;
    std::function<void(Domain &, const Range &)> Tvmult;
    std::function<void(Domain &, const Range &)> Tvmult_add;
@endcode

然而，LinearOperator和普通矩阵的关键区别在于，LinearOperator不允许对底层对象进行任何进一步的访问。你能用LinearOperator做的就是把它的 "动作 "应用到一个向量上。我们借此机会介绍一下LinearOperator的概念，因为它是一个非常有用的工具，可以让你以非常直观的方式构造复杂的求解器和预处理器。

作为第一个例子，让我们构建一个代表  $M^{-1}$  的 LinearOperator 对象。这意味着每当这个运算符的 <code>vmult()</code> 函数被调用时，它必须解决一个线性系统。这就要求我们指定一个解算器（和相应的）前置条件。假设 <code>M</code> 是对系统矩阵左上块的引用，我们可以写出。

@code
    const auto op_M = linear_operator(M);


    PreconditionJacobi<SparseMatrix<double>> preconditioner_M;
    preconditioner_M.initialize(M);


    ReductionControl reduction_control_M(2000, 1.0e-18, 1.0e-10);
    SolverCG<Vector<double>>    solver_M(reduction_control_M);


    const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
@endcode

我们没有使用SolverControl，而是在这里使用ReductionControl类，当达到绝对公差（我们选择 $10^{-18}$ ）或者当残差减少了某个系数（这里是 $10^{-10}$ ）时，它就会停止迭代。相反，SolverControl类只检查绝对公差。在我们的案例中，我们必须使用ReductionControl来解决一个小问题。我们将送入 <code>op_M_inv</code> 的右手边基本上是由残差形成的，随着外部迭代的进行，残差的规范自然会大大减少。这使得绝对公差的控制非常容易出错。

我们现在有一个LinearOperator  <code>op_M_inv</code> ，我们可以用它来构造更复杂的运算符，如Schur补码  $S$  。假设 <code>B</code> 是对右上角区块的引用，构造一个LinearOperator  <code>op_S</code> 只需两行即可。

@code
    const auto op_B = linear_operator(B);
    const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
@endcode

在这里，三个LinearOperator对象的乘法产生了一个复合对象 <code>op_S</code> whose <code>vmult()</code> ，该函数首先应用 $B$ ，然后是 $M^{-1}$ （即用 $M$ 解方程），最后是 $B^T$ 到任何指定的输入矢量。在这个意义上， <code>op_S.vmult()</code> 类似于以下代码。

@code
    B.vmult (tmp1, src); // multiply with the top right block: B
    solver_M(M, tmp2, tmp1, preconditioner_M); // multiply with M^-1
    B.Tvmult (dst, tmp2); // multiply with the bottom left block: B^T
@endcode

(  <code>tmp1</code> and <code>tmp2</code> 是两个临时向量)。这种方法背后的关键点是，我们实际上从未创建一个矩阵的内积。相反，每当我们要用 <code>op_S</code> 进行矩阵向量乘法时，我们只需按上述顺序运行所有单独的 <code>vmult</code> 操作。

 @note 我们可以通过实现一个专门的类 <code>SchurComplement</code> ，提供一个合适的 <code>vmult()</code> 函数，来实现创建一个 "类似矩阵 "的对象的相同目标。跳过一些细节，这可能看起来像下面这样。

@code
class SchurComplement
{
  public:


  // ...


  void SchurComplement::vmult (Vector<double>       &dst,
                               const Vector<double> &src) const
  {
    B.vmult (tmp1, src);
    solver_M(M, tmp2, tmp1, preconditioner_M);
    B.Tvmult (dst, tmp2);
  }
};
@endcode

尽管这两种方法完全等同，但LinearOperator类比这种手工方法有很大的优势。它提供了所谓的［<i><a href="https://en.wikipedia.org/wiki/Syntactic_sugar">syntactic sugar</a>］［<a href="https://en.wikipedia.org/wiki/Syntactic_sugar">syntactic sugar</a></i>］。在数学上，我们认为 $S$ 是复合矩阵 $S=B^TM^{-1}B$ ，LinearOperator类允许你或多或少地逐字写出这一点。

@code
const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M);
const auto op_S = transpose_operator(op_B) * op_M_inv * op_B;
@endcode

另一方面，人工方法掩盖了这一事实。

现在我们要做的就是形成定义 $P$ 和 $U$ 的两个方程的右手边，然后分别用舒尔补码矩阵和质量矩阵来解决它们。例如，第一个方程的右手边为 $B^TM^{-1}F-G$  。这可以通过以下方式实现。

@code
    Vector<double> schur_rhs (P.size());
    Vector<double> tmp (U.size());
    op_M_inv.vmult (tmp, F);
    transpose_operator(op_B).vmult (schur_rhs, tmp);
    schur_rhs -= G;
@endcode

同样，这是一个完全有效的方法，但是deal.II要求我们手动调整最终向量和临时向量的大小，而且每一个操作都要占用一个新的行，这让我们很难阅读。这就是线性运算符框架中的第二个类可以将帮助我们的地方。与LinearOperator的精神类似，一个PackagedOperation存储了一个 "计算"。

@code
    std::function<void(Range &)> apply;
    std::function<void(Range &)> apply_add;
@endcode

该类允许<a href="https://en.wikipedia.org/wiki/Lazy_evaluation">lazy evaluation</a>涉及向量和线性运算符的表达式。这是通过存储计算表达式来实现的，只有当对象被转换为向量对象，或者 PackagedOperation::apply() （或 PackagedOperation::apply_add()) 被手动调用时才执行计算。假设 <code>F</code> and <code>G</code> 是右手边的两个向量，我们可以简单地写。

@code
    const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G;
@endcode

这里， <code>schur_rhs</code> 是一个打包操作，<i>records</i>是我们指定的计算。它不会立即创建一个带有实际结果的向量。

有了这些先决条件，解决 $P$ 和 $U$ 的问题就是创建另一个解算器和反演。

@code
    SolverControl solver_control_S(2000, 1.e-12);
    SolverCG<Vector<double>>    solver_S(solver_control_S);
    PreconditionIdentity preconditioner_S;


    const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S);


    P = op_S_inv * schur_rhs;
    U = op_M_inv * (F - op_B * P);
@endcode



 @note  我们在这个例子中手工开发的功能在库中已经是现成的。看看schur_complement(), condense_schur_rhs(), and postprocess_schur_solution()。




<a name="ApreconditionerfortheSchurcomplement"></a><h4>A preconditioner for the Schur complement</h4>


有人可能会问，如果我们有一个Schur补数的预处理程序，是否会有帮助  $S=B^TM^{-1}B$  。一般来说，答案是：当然。问题是，我们对这个舒尔补码矩阵一无所知。我们不知道它的条目，我们所知道的只是它的作用。另一方面，我们必须认识到，我们的求解器是昂贵的，因为在每次迭代中，我们必须与舒尔补矩阵做一次矩阵-向量乘积，这意味着我们必须在每次迭代中对质量矩阵做一次反转。

对这样一个矩阵的预处理有不同的方法。一个极端是使用便宜的东西，因此对每次迭代的工作没有实际影响。另一个极端是使用本身非常昂贵的预处理程序，但作为回报，它确实降低了用 $S$ 求解所需的迭代次数。

我们将按照第二种方法进行尝试，既是为了提高程序的性能，也是为了展示一些技术。为此，让我们回顾一下，理想的预处理程序当然是 $S^{-1}$ ，但这是无法实现的。然而，如何

@f{eqnarray*}
  \tilde S^{-1} = [B^T ({\textrm{diag}\ }M)^{-1}B]^{-1}


@f}

作为一个预处理程序？这就意味着，每次我们要做一个预处理步骤时，实际上都要用 $\tilde S$ 来解。起初，这看起来几乎和立即用 $S$ 求解一样昂贵。然而，请注意，在内迭代中，我们不必计算 $M^{-1}$ ，而只需计算其对角线的逆值，这很便宜。

值得庆幸的是，LinearOperator框架使得这一点非常容易写出来。我们之前已经对 <code>preconditioner_M</code> 矩阵使用了雅可比预处理程序（ $M$ ）。所以剩下的就是写出近似的舒尔补码应该是什么样子。

@code
    const auto op_aS =
      transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B;
@endcode

注意这个运算符的不同之处在于，它只是做了一次雅可比扫频（即与对角线的逆数相乘），而不是与整个 $M^{-1}$ 相乘。]的定义：它是与 $M$ 的对角线的倒数相乘；换句话说，对向量 $x$ 的 $({\textrm{diag}\ }M)^{-1}x$ 操作正是PreconditionJacobi所做的）。)

有了这些，我们几乎已经完成了预处理程序：它应该是近似Schur补码的逆。我们再次通过使用inverse_operator()函数创建一个线性算子来实现这一点。然而这一次我们想为CG求解器选择一个相对较小的容忍度（即反转 <code>op_aS</code> ）。理由是 <code>op_aS</code> is only coarse approximation to <code>op_S</code> ，所以我们实际上不需要完全反转它。然而，这产生了一个微妙的问题： <code>preconditioner_S</code> 将被用于最后的外层CG迭代以创建一个正交基础。但为了使其发挥作用，每次调用都必须是精确的线性操作。我们通过使用IterationNumberControl来确保这一点，它允许我们将执行的CG迭代次数固定为一个固定的小数字（在我们的例子中为30）。

@code
    IterationNumberControl iteration_number_control_aS(30, 1.e-18);
    SolverCG<Vector<double>>           solver_aS(iteration_number_control_aS);
    PreconditionIdentity preconditioner_aS;
    const auto preconditioner_S =
      inverse_operator(op_aS, solver_aS, preconditioner_aS);
@endcode



就这样吧!

很明显，应用这个近似Schur补数的逆运算是一个非常昂贵的预处理程序，几乎和反转Schur补数本身一样昂贵。我们可以期望它能大大减少Schur补数所需的外部迭代次数。事实上，它确实如此：在使用0阶元素的7次细化网格的典型运行中，外部迭代次数从592次下降到39次。另一方面，我们现在必须应用一个非常昂贵的预处理程序25次。因此，更好的衡量标准只是程序的运行时间：在目前的笔记本电脑上（截至2019年1月），对于这个测试案例，它从3.57秒下降到2.05秒。这似乎并不令人印象深刻，但在更细的网格和更高阶的元素上，节省的时间变得更加明显了。例如，一个7倍细化的网格和使用2阶元素（相当于约40万个自由度）产生了1134次到83次的外部迭代的改进，运行时间从168秒到40秒。虽然不是惊天动地，但意义重大。




<a name="Definitionofthetestcase"></a><h3>Definition of the test case</h3>


在这个教程程序中，我们将求解上述混合公式中的拉普拉斯方程。由于我们想在程序中监测解的收敛性，我们选择右手边、边界条件和系数，以便恢复我们已知的解函数。特别是，我们选择压力解

@f{eqnarray*}
  p = -\left(\frac \alpha 2 xy^2 + \beta x - \frac \alpha 6 x^3\right),


@f}

而对于系数，为了简单起见，我们选择单位矩阵 $K_{ij}=\delta_{ij}$ 。因此，确切的速度满足于

@f{eqnarray*}
  {\mathbf u} =
  \left(\begin{array}{cc}
    \frac \alpha 2 y^2 + \beta - \frac \alpha 2 x^2 \\
    \alpha xy
  \end{array}\right).


@f}

选择这个解决方案是因为它完全没有发散，使得它成为不可压缩流体流动的一个现实的测试案例。因此，右手边等于 $f=0$ ，作为边界值，我们必须选择 $g=p|_{\partial\Omega}$ 。

对于本程序中的计算，我们选择  $\alpha=0.3,\beta=1$  。你可以在<a name="#Results">results section
below</a>中找到结果的解决方案，在注释的程序之后。


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
 * 由于这个程序只是对 step-4 的改编，所以在头文件方面没有太多的新东西。在deal.II中，我们通常按照base-lac-grid-dofs-fe-numerics的顺序列出包含文件，然后是C++标准包含文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/function.h> 
 * 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * 
 * @endcode
 * 
 * 唯一值得关注的两个新头文件是LinearOperator和PackagedOperation类的文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/linear_operator.h> 
 * #include <deal.II/lac/packaged_operation.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 这是唯一重要的新标题，即声明Raviart-Thomas有限元的标题。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_raviart_thomas.h> 
 * 
 * @endcode
 * 
 * 最后，作为本程序中的一项奖励，我们将使用一个张量系数。由于它可能具有空间依赖性，我们认为它是一个张量值的函数。下面的include文件提供了 <code>TensorFunction</code> 类，提供了这样的功能。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/tensor_function.h> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。我们把所有与这个程序相关的代码放到一个命名空间中。(这个想法在  step-7  中首次提出) 。
 * 

 * 
 * 
 * @code
 * namespace Step20 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeMixedLaplaceProblemcodeclasstemplate"></a> 
 * <h3>The <code>MixedLaplaceProblem</code> class template</h3>
 * 

 * 
 * 同样，由于这是对 step-6 的改编，主类与该教程程序中的主类几乎相同。就成员函数而言，主要区别在于构造函数将Raviart-Thomas元素的度数作为参数（并且有一个相应的成员变量来存储这个值），并且增加了 <code>compute_error</code> 函数，在这个函数中，不出意外，我们将计算精确解和数值解之间的差异，以确定我们计算的收敛性。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class MixedLaplaceProblem 
 *   { 
 *   public: 
 *     MixedLaplaceProblem(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid_and_dofs(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void compute_errors() const; 
 *     void output_results() const; 
 * 
 *     const unsigned int degree; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 * @endcode
 * 
 * 第二个区别是疏散模式、系统矩阵、解和右手向量现在被封锁了。这意味着什么，人们可以用这些对象做什么，在本程序的介绍中已经解释过了，下面我们在解释这个问题的线性求解器和预处理器时也会进一步解释。
 * 

 * 
 * 
 * @code
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> system_matrix; 
 * 
 *     BlockVector<double> solution; 
 *     BlockVector<double> system_rhs; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandexactsolution"></a> 
 * <h3>Right hand side, boundary values, and exact solution</h3>
 * 

 * 
 * 我们的下一个任务是定义我们问题的右手边（即原始拉普拉斯方程中压力的标量右手边），压力的边界值，以及一个描述压力和精确解的速度的函数，以便以后计算误差。请注意，这些函数分别有一个、一个和 <code>dim+1</code> 个分量，我们将分量的数量传递给 <code>Function@<dim@></code> 基类。对于精确解，我们只声明实际一次性返回整个解向量（即其中的所有成分）的函数。下面是各自的声明。
 * 

 * 
 * 
 * @code
 *   namespace PrescribedSolution 
 *   { 
 *     constexpr double alpha = 0.3; 
 *     constexpr double beta  = 1; 
 * 
 *     template <int dim> 
 *     class RightHandSide : public Function<dim> 
 *     { 
 *     public: 
 *       RightHandSide() 
 *         : Function<dim>(1) 
 *       {} 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 *     }; 
 * 
 *     template <int dim> 
 *     class PressureBoundaryValues : public Function<dim> 
 *     { 
 *     public: 
 *       PressureBoundaryValues() 
 *         : Function<dim>(1) 
 *       {} 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 *     }; 
 * 
 *     template <int dim> 
 *     class ExactSolution : public Function<dim> 
 *     { 
 *     public: 
 *       ExactSolution() 
 *         : Function<dim>(dim + 1) 
 *       {} 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  value) const override; 
 *     }; 
 * 
 * @endcode
 * 
 * 然后我们还必须定义这些各自的函数，当然了。鉴于我们在介绍中讨论了解决方案应该是怎样的，下面的计算应该是很简单的。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     double RightHandSide<dim>::value(const Point<dim> & /*p*/, 
 *                                      const unsigned int /*component*/) const 
 *     { 
 *       return 0; 
 *     } 
 * 
 *     template <int dim> 
 *     double 
 *     PressureBoundaryValues<dim>::value(const Point<dim> &p, 
 *                                        const unsigned int /*component*/) const 
 *     { 
 *       return -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] - 
 *                alpha * p[0] * p[0] * p[0] / 6); 
 *     } 
 * 
 *     template <int dim> 
 *     void ExactSolution<dim>::vector_value(const Point<dim> &p, 
 *                                           Vector<double> &  values) const 
 *     { 
 *       Assert(values.size() == dim + 1, 
 *              ExcDimensionMismatch(values.size(), dim + 1)); 
 * 
 *       values(0) = alpha * p[1] * p[1] / 2 + beta - alpha * p[0] * p[0] / 2; 
 *       values(1) = alpha * p[0] * p[1]; 
 *       values(2) = -(alpha * p[0] * p[1] * p[1] / 2 + beta * p[0] - 
 *                     alpha * p[0] * p[0] * p[0] / 6); 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="Theinversepermeabilitytensor"></a> 
 * <h3>The inverse permeability tensor</h3>
 * 

 * 
 * 除了其他方程数据外，我们还想使用渗透性张量，或者更好的是--因为这是在弱形式中出现的全部内容--渗透性张量的逆，  <code>KInverse</code>  。对于验证解的精确性和确定收敛顺序的目的来说，这个张量的作用大于帮助。因此，我们将简单地把它设置为同一矩阵。
 * 

 * 
 * 然而，在现实生活中的多孔介质流动模拟中，空间变化的渗透率张量是不可缺少的，我们想利用这个机会来展示使用张量值函数的技术。
 * 

 * 
 * 可能不足为奇，deal.II也有一个基类，不仅适用于标量和一般的矢量值函数（ <code>Function</code> 基类），也适用于返回固定维度和等级的张量的函数， <code>TensorFunction</code> 模板。在这里，所考虑的函数返回一个dim-by-dim矩阵，即一个等级为2、维度为 <code>dim</code> 的张量。然后我们适当地选择基类的模板参数。
 * 

 * 
 * <code>TensorFunction</code> 类提供的接口本质上等同于 <code>Function</code> 类。特别是，存在一个 <code>value_list</code> 函数，它接收一个评估函数的点的列表，并在第二个参数中返回函数的值，一个张量的列表。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class KInverse : public TensorFunction<2, dim> 
 *     { 
 *     public: 
 *       KInverse() 
 *         : TensorFunction<2, dim>() 
 *       {} 
 * 
 *       virtual void 
 *       value_list(const std::vector<Point<dim>> &points, 
 *                  std::vector<Tensor<2, dim>> &  values) const override; 
 *     }; 
 * 
 * @endcode
 * 
 * 实现起来就不那么有趣了。和以前的例子一样，我们在类的开头添加一个检查，以确保输入和输出参数的大小是相同的（关于这个技术的讨论见 step-5 ）。然后我们在所有的评估点上循环，对于每一个评估点，将输出张量设置为身份矩阵。
 * 

 * 
 * 在函数的顶部有一个奇怪的地方（`(void)point;`语句），值得讨论。我们放到输出`values`数组中的值实际上并不取决于函数被评估的坐标`points`数组。换句话说，`points'参数实际上是不用的，如果我们想的话，可以不给它起名字。但是我们想用`points`对象来检查`values`对象是否有正确的大小。问题是，在发布模式下，`AssertDimension`被定义为一个宏，扩展为空；然后编译器会抱怨`points`对象没有使用。消除这个警告的习惯方法是有一个评估（读取）变量的语句，但实际上不做任何事情：这就是`(void)points;`所做的：它从`points`中读取，然后将读取的结果转换为`void`，也就是什么都没有。换句话说，这句话是完全没有意义的，除了向编译器解释是的，这个变量事实上是被使用的，即使是在发布模式下。(在调试模式下，`AssertDimension`宏会扩展为从变量中读出的东西，所以在调试模式下，这个有趣的语句是没有必要的)。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
 *                                    std::vector<Tensor<2, dim>> &  values) const 
 *     { 
 *       (void)points; 
 *       AssertDimension(points.size(), values.size()); 
 * 
 *       for (auto &value : values) 
 *         value = unit_symmetric_tensor<dim>(); 
 *     } 
 *   } // namespace PrescribedSolution 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemclassimplementation"></a> 
 * <h3>MixedLaplaceProblem class implementation</h3>
 * 
 * <a name="MixedLaplaceProblemMixedLaplaceProblem"></a> 
 * <h4>MixedLaplaceProblem::MixedLaplaceProblem</h4>
 * 

 * 
 * 在这个类的构造函数中，我们首先存储传入的关于我们将使用的有限元的度数的值（例如，度数为0，意味着使用RT(0)和DG(0)），然后构造属于介绍中描述的空间 $X_h$ 的向量值的元素。构造函数的其余部分与早期的教程程序一样。
 * 

 * 
 * 这里唯一值得描述的是，这个变量所属的 <code>fe</code> variable. The <code>FESystem</code> 类的构造函数调用有很多不同的构造函数，它们都是指将较简单的元素绑定在一起，成为一个较大的元素。在目前的情况下，我们想把一个RT(度)元素与一个DQ(度)元素结合起来。这样做的 <code>FESystem</code> 构造函数要求我们首先指定第一个基本元素（给定程度的 <code>FE_RaviartThomas</code> 对象），然后指定这个基本元素的副本数量，然后类似地指定 <code>FE_DGQ</code> 元素的种类和数量。注意Raviart-Thomas元素已经有 <code>dim</code> 个矢量分量，所以耦合元素将有 <code>dim+1</code> 个矢量分量，其中第一个 <code>dim</code> 个对应于速度变量，最后一个对应于压力。
 * 

 * 
 * 我们从基本元素中构建这个元素的方式与我们在 step-8 中的方式也值得比较：在那里，我们将其构建为 <code>fe (FE_Q@<dim@>(1), dim)</code> ，即我们简单地使用 <code>dim</code> copies of the <code>FE_Q(1)</code> 元素，每个坐标方向上的位移都有一份。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   MixedLaplaceProblem<dim>::MixedLaplaceProblem(const unsigned int degree) 
 *     : degree(degree) 
 *     , fe(FE_RaviartThomas<dim>(degree), 1, FE_DGQ<dim>(degree), 1) 
 *     , dof_handler(triangulation) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemmake_grid_and_dofs"></a> 
 * <h4>MixedLaplaceProblem::make_grid_and_dofs</h4>
 * 

 * 
 * 接下来的函数开始于众所周知的函数调用，创建和细化一个网格，然后将自由度与之关联。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::make_grid_and_dofs() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -1, 1); 
 *     triangulation.refine_global(5); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 然而，接下来事情就变得不同了。正如介绍中提到的，我们要将矩阵细分为对应于速度和压力这两种不同的变量的块。为此，我们首先要确保与速度和压力相对应的指数不会混在一起。首先是所有速度自由度，然后是所有压力自由度。这样一来，全局矩阵就很好地分离成一个 $2 \times 2$ 系统。为了达到这个目的，我们必须根据自由度的矢量分量对其重新编号，这个操作已经很方便地实现了。
 * 

 * 
 * 
 * @code
 *     DoFRenumbering::component_wise(dof_handler); 
 * 
 * @endcode
 * 
 * 接下来，我们要弄清楚这些块的大小，以便我们可以分配适当的空间量。为此，我们调用了 DoFTools::count_dofs_per_fe_component() 函数，该函数计算了某个向量分量的形状函数非零的数量。我们有 <code>dim+1</code> 个向量分量， DoFTools::count_dofs_per_fe_component() 将计算有多少个形状函数属于这些分量中的每个。
 * 

 * 
 * 这里有一个问题。正如该函数的文档所描述的，它 <i>wants</i> 将  $x$  -速度形状函数的数量放入  <code>dofs_per_component[0]</code>  中，将  $y$  -速度形状函数的数量放入  <code>dofs_per_component[1]</code>  中（以及类似的3d），并将压力形状函数的数量放入  <code>dofs_per_component[dim]</code>  中 。但是，Raviart-Thomas元素的特殊性在于它是非 @ref GlossPrimitive "原始 "的，也就是说，对于Raviart-Thomas元素，所有的速度形状函数在所有分量中都是非零。换句话说，该函数不能区分 $x$ 和 $y$ 速度函数，因为<i>is</i>没有这种区分。因此，它将速度的总体数量放入 <code>dofs_per_component[c]</code>  ,  $0\le c\le \text{dim}$ 中的每一个。另一方面，压力变量的数量等于在dim-th分量中不为零的形状函数的数量。
 * 

 * 
 * 利用这些知识，我们可以从 <code>dofs_per_component</code> 的第一个 <code>dim</code> 元素中的任何一个得到速度形状函数的数量，然后用下面这个来初始化向量和矩阵块的大小，以及创建输出。
 * 

 * 
 * @note  如果你觉得这个概念难以理解，你可以考虑用函数  DoFTools::count_dofs_per_fe_block()  来代替，就像我们在  step-22  的相应代码中做的那样。你可能还想阅读一下术语表中 @ref GlossBlock "块 "和 @ref GlossComponent "组件 "的区别。
 * 

 * 
 * 
 * @code
 *     const std::vector<types::global_dof_index> dofs_per_component = 
 *       DoFTools::count_dofs_per_fe_component(dof_handler); 
 *     const unsigned int n_u = dofs_per_component[0], 
 *                        n_p = dofs_per_component[dim]; 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "Total number of cells: " << triangulation.n_cells() 
 *               << std::endl 
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (" << n_u << '+' << n_p << ')' << std::endl; 
 * 
 * @endcode
 * 
 * 下一个任务是为我们将要创建的矩阵分配一个稀疏模式。我们使用与前面步骤一样的压缩稀疏模式，但是由于 <code>system_matrix</code> 是一个块状矩阵，我们使用 <code>BlockDynamicSparsityPattern</code> 类，而不仅仅是 <code>DynamicSparsityPattern</code>  。这种块状稀疏模式在 $2 \times 2$ 模式下有四个块。块的大小取决于 <code>n_u</code> and <code>n_p</code> ，它持有速度和压力变量的数量。在第二步中，我们必须指示块系统更新它所管理的块的大小的知识；这发生在 <code>dsp.collect_sizes ()</code> 的调用中。
 * 

 * 
 * 
 * @code
 *     BlockDynamicSparsityPattern dsp(2, 2); 
 *     dsp.block(0, 0).reinit(n_u, n_u); 
 *     dsp.block(1, 0).reinit(n_p, n_u); 
 *     dsp.block(0, 1).reinit(n_u, n_p); 
 *     dsp.block(1, 1).reinit(n_p, n_p); 
 *     dsp.collect_sizes(); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 * 
 * @endcode
 * 
 * 我们以与非区块版本相同的方式使用压缩的区块稀疏模式，以创建稀疏模式，然后创建系统矩阵。
 * 

 * 
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp); 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 然后，我们必须以与块压缩稀疏度模式完全相同的方式调整解决方案和右侧向量的大小。
 * 

 * 
 * 
 * @code
 *     solution.reinit(2); 
 *     solution.block(0).reinit(n_u); 
 *     solution.block(1).reinit(n_p); 
 *     solution.collect_sizes(); 
 * 
 *     system_rhs.reinit(2); 
 *     system_rhs.block(0).reinit(n_u); 
 *     system_rhs.block(1).reinit(n_p); 
 *     system_rhs.collect_sizes(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemassemble_system"></a> 
 * <h4>MixedLaplaceProblem::assemble_system</h4>
 * 

 * 
 * 同样地，组装线性系统的函数在这个例子的介绍中已经讨论过很多了。在它的顶部，发生的是所有常见的步骤，此外，我们不仅为单元项分配正交和 <code>FEValues</code> 对象，而且还为面项分配。之后，我们为变量定义通常的缩写，并为本地矩阵和右手贡献分配空间，以及保存当前单元的全局自由度数的数组。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::assemble_system() 
 *   { 
 *     QGauss<dim>     quadrature_formula(degree + 2); 
 *     QGauss<dim - 1> face_quadrature_formula(degree + 2); 
 * 
 *     FEValues<dim>     fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 *     FEFaceValues<dim> fe_face_values(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_normal_vectors | 
 *                                        update_quadrature_points | 
 *                                        update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points      = quadrature_formula.size(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 下一步是声明代表方程中源项、压力边界值和系数的对象。除了这些代表连续函数的对象外，我们还需要数组来保存它们在各个单元格（或面，对于边界值）的正交点的值。请注意，在系数的情况下，数组必须是矩阵的一种。
 * 

 * 
 * 
 * @code
 *     const PrescribedSolution::RightHandSide<dim> right_hand_side; 
 *     const PrescribedSolution::PressureBoundaryValues<dim> 
 *                                             pressure_boundary_values; 
 *     const PrescribedSolution::KInverse<dim> k_inverse; 
 * 
 *     std::vector<double>         rhs_values(n_q_points); 
 *     std::vector<double>         boundary_values(n_face_q_points); 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 
 * 
 * @endcode
 * 
 * 最后，我们需要几个提取器，用来获取矢量值形状函数的速度和压力成分。它们的功能和使用在 @ref vector_valued报告中有详细描述。基本上，我们将把它们作为下面FEValues对象的下标：FEValues对象描述了形状函数的所有矢量分量，而在订阅后，它将只指速度（一组从零分量开始的 <code>dim</code> 分量）或压力（位于 <code>dim</code> 位置的标量分量）。
 * 

 * 
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 * @endcode
 * 
 * 有了这些，我们就可以继续对所有单元进行循环。这个循环的主体已经在介绍中讨论过了，这里就不再做任何评论了。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 *         right_hand_side.value_list(fe_values.get_quadrature_points(), 
 *                                    rhs_values); 
 *         k_inverse.value_list(fe_values.get_quadrature_points(), 
 *                              k_inverse_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const Tensor<1, dim> phi_i_u = fe_values[velocities].value(i, q); 
 *               const double div_phi_i_u = fe_values[velocities].divergence(i, q); 
 *               const double phi_i_p     = fe_values[pressure].value(i, q); 
 * 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const Tensor<1, dim> phi_j_u = 
 *                     fe_values[velocities].value(j, q); 
 *                   const double div_phi_j_u = 
 *                     fe_values[velocities].divergence(j, q); 
 *                   const double phi_j_p = fe_values[pressure].value(j, q); 
 * 
 *                   local_matrix(i, j) += 
 *                     (phi_i_u * k_inverse_values[q] * phi_j_u // 
 *                      - phi_i_p * div_phi_j_u                 // 
 *                      - div_phi_i_u * phi_j_p)                // 
 *                     * fe_values.JxW(q); 
 *                 } 
 * 
 *               local_rhs(i) += -phi_i_p * rhs_values[q] * fe_values.JxW(q); 
 *             } 
 * 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary()) 
 *             { 
 *               fe_face_values.reinit(cell, face); 
 * 
 *               pressure_boundary_values.value_list( 
 *                 fe_face_values.get_quadrature_points(), boundary_values); 
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   local_rhs(i) += -(fe_face_values[velocities].value(i, q) * // 
 *                                     fe_face_values.normal_vector(q) *        // 
 *                                     boundary_values[q] *                     // 
 *                                     fe_face_values.JxW(q)); 
 *             } 
 * 
 * @endcode
 * 
 * 循环所有单元的最后一步是将局部贡献转移到全局矩阵和右手向量中。请注意，我们使用的接口与之前的例子完全相同，尽管我们现在使用的是块状矩阵和向量，而不是常规的。换句话说，对于外界来说，块对象具有与矩阵和向量相同的接口，但它们还允许访问单个块。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             system_matrix.add(local_dof_indices[i], 
 *                               local_dof_indices[j], 
 *                               local_matrix(i, j)); 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           system_rhs(local_dof_indices[i]) += local_rhs(i); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Implementationoflinearsolversandpreconditioners"></a> 
 * <h3>Implementation of linear solvers and preconditioners</h3>
 * 

 * 
 * 我们在这个例子中使用的线性求解器和预处理器已经在介绍中进行了详细的讨论。因此，我们在这里不再讨论我们的方法的原理，而只是对剩下的一些实现方面进行评论。
 * 

 * 
 * 
 * <a name="MixedLaplacesolve"></a> 
 * <h4>MixedLaplace::solve</h4>
 * 

 * 
 * 正如在介绍中所概述的那样，求解函数基本上由两个步骤组成。首先，我们必须形成涉及舒尔补数的第一个方程，并求解压力（解决方案的第一部分）。然后，我们可以从第二个方程（解的第0部分）中重构速度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::solve() 
 *   { 
 * 
 * @endcode
 * 
 * 作为第一步，我们声明对矩阵的所有块状成分、右手边和我们将需要的解向量的引用。
 * 

 * 
 * 
 * @code
 *     const auto &M = system_matrix.block(0, 0); 
 *     const auto &B = system_matrix.block(0, 1); 
 * 
 *     const auto &F = system_rhs.block(0); 
 *     const auto &G = system_rhs.block(1); 
 * 
 *     auto &U = solution.block(0); 
 *     auto &P = solution.block(1); 
 * 
 * @endcode
 * 
 * 然后，我们将创建相应的LinearOperator对象并创建 <code>op_M_inv</code> 运算器。
 * 

 * 
 * 
 * @code
 *     const auto op_M = linear_operator(M); 
 *     const auto op_B = linear_operator(B); 
 * 
 *     ReductionControl         reduction_control_M(2000, 1.0e-18, 1.0e-10); 
 *     SolverCG<Vector<double>> solver_M(reduction_control_M); 
 *     PreconditionJacobi<SparseMatrix<double>> preconditioner_M; 
 * 
 *     preconditioner_M.initialize(M); 
 * 
 *     const auto op_M_inv = inverse_operator(op_M, solver_M, preconditioner_M); 
 * 
 * @endcode
 * 
 * 这样我们就可以声明舒尔补数  <code>op_S</code>  和近似舒尔补数  <code>op_aS</code>  。
 * 

 * 
 * 
 * @code
 *     const auto op_S = transpose_operator(op_B) * op_M_inv * op_B; 
 *     const auto op_aS = 
 *       transpose_operator(op_B) * linear_operator(preconditioner_M) * op_B; 
 * 
 * @endcode
 * 
 * 我们现在从 <code>op_aS</code> 中创建一个预处理程序，应用固定数量的30次（便宜的）CG迭代。
 * 

 * 
 * 
 * @code
 *     IterationNumberControl   iteration_number_control_aS(30, 1.e-18); 
 *     SolverCG<Vector<double>> solver_aS(iteration_number_control_aS); 
 * 
 *     const auto preconditioner_S = 
 *       inverse_operator(op_aS, solver_aS, PreconditionIdentity()); 
 * 
 * @endcode
 * 
 * 现在来看看第一个方程。它的右边是 $B^TM^{-1}F-G$  ，这就是我们在前几行计算的结果。然后我们用CG求解器和我们刚刚声明的预处理程序来解决第一个方程。
 * 

 * 
 * 
 * @code
 *     const auto schur_rhs = transpose_operator(op_B) * op_M_inv * F - G; 
 * 
 *     SolverControl            solver_control_S(2000, 1.e-12); 
 *     SolverCG<Vector<double>> solver_S(solver_control_S); 
 * 
 *     const auto op_S_inv = inverse_operator(op_S, solver_S, preconditioner_S); 
 * 
 *     P = op_S_inv * schur_rhs; 
 * 
 *     std::cout << solver_control_S.last_step() 
 *               << " CG Schur complement iterations to obtain convergence." 
 *               << std::endl; 
 * 
 * @endcode
 * 
 * 得到压力后，我们可以计算速度。方程为 $MU=-BP+F$  ，我们通过首先计算右手边，然后与代表质量矩阵逆的对象相乘来解决这个问题。
 * 

 * 
 * 
 * @code
 *     U = op_M_inv * (F - op_B * P); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceProblemclassimplementationcontinued"></a> 
 * <h3>MixedLaplaceProblem class implementation (continued)</h3>
 * 
 * <a name="MixedLaplacecompute_errors"></a> 
 * <h4>MixedLaplace::compute_errors</h4>
 * 

 * 
 * 在我们处理完线性求解器和预处理器之后，我们继续实现我们的主类。特别是，下一个任务是计算我们数值解的误差，包括压力和速度。
 * 

 * 
 * 为了计算解的误差，我们已经在  step-7  和  step-11  中介绍了  <code>VectorTools::integrate_difference</code>  函数。然而，在那里我们只处理了标量解，而在这里我们有一个矢量值的解，其组成部分甚至表示不同的量，并且可能有不同的收敛阶数（由于所使用的有限元的选择，这里不是这种情况，但在混合有限元应用中经常出现这种情况）。因此，我们要做的是 "掩盖 "我们感兴趣的成分。这很容易做到： <code>VectorTools::integrate_difference</code> 函数将一个指向权重函数的指针作为其参数之一（该参数默认为空指针，意味着单位权重）。我们要做的是传递一个函数对象，在我们感兴趣的成分中等于1，而在其他成分中等于0。例如，为了计算压力误差，我们应该传入一个函数，该函数在分量 <code>dim</code> 中代表单位值的常数向量，而对于速度，常数向量在第一个 <code>dim</code> 分量中应该是1，而在压力的位置是0。
 * 

 * 
 * 在deal.II中， <code>ComponentSelectFunction</code> 正是这样做的：它想知道它要表示的函数应该有多少个向量分量（在我们的例子中，这将是 <code>dim+1</code> ，用于联合速度-压力空间），哪个个体或范围的分量应该等于1。因此，我们在函数的开头定义了两个这样的掩码，接下来是一个代表精确解的对象和一个向量，我们将在其中存储由 <code>integrate_difference</code> 计算的单元误差。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::compute_errors() const 
 *   { 
 *     const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1); 
 *     const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), 
 *                                                      dim + 1); 
 * 
 *     PrescribedSolution::ExactSolution<dim> exact_solution; 
 *     Vector<double> cellwise_errors(triangulation.n_active_cells()); 
 * 
 * @endcode
 * 
 * 正如在 step-7 中已经讨论过的那样，我们必须认识到，不可能精确地整合误差。我们所能做的就是用正交法对这个积分进行近似。这实际上在这里提出了一个小小的转折：如果我们像人们可能倾向于做的那样天真地选择一个 <code>QGauss@<dim@>(degree+1)</code> 类型的对象（这就是我们用于积分线性系统的对象），就会发现误差非常小，根本不遵循预期的收敛曲线。现在的情况是，对于这里使用的混合有限元，高斯点恰好是超收敛点，其中的点误差要比其他地方小得多（而且收敛的阶数更高）。因此，这些点不是特别好的积分点。为了避免这个问题，我们只需使用梯形法则，并在每个坐标方向上迭代 <code>degree+2</code> 次（同样如 step-7 中的解释）。
 * 

 * 
 * 
 * @code
 *     QTrapezoid<1>  q_trapez; 
 *     QIterated<dim> quadrature(q_trapez, degree + 2); 
 * 
 * @endcode
 * 
 * 有了这个，我们就可以让库计算出误差并将其输出到屏幕上。
 * 

 * 
 * 
 * @code
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       exact_solution, 
 *                                       cellwise_errors, 
 *                                       quadrature, 
 *                                       VectorTools::L2_norm, 
 *                                       &pressure_mask); 
 *     const double p_l2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         cellwise_errors, 
 *                                         VectorTools::L2_norm); 
 * 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       exact_solution, 
 *                                       cellwise_errors, 
 *                                       quadrature, 
 *                                       VectorTools::L2_norm, 
 *                                       &velocity_mask); 
 *     const double u_l2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         cellwise_errors, 
 *                                         VectorTools::L2_norm); 
 * 
 *     std::cout << "Errors: ||e_p||_L2 = " << p_l2_error 
 *               << ",   ||e_u||_L2 = " << u_l2_error << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="MixedLaplaceoutput_results"></a> 
 * <h4>MixedLaplace::output_results</h4>
 * 

 * 
 * 最后一个有趣的函数是我们生成图形输出的函数。请注意，所有的速度分量都得到相同的解名 "u"。再加上使用 DataComponentInterpretation::component_is_part_of_vector ，这将导致 DataOut<dim>::write_vtu() 生成各个速度分量的矢量表示，更多信息请参见 step-22 或 @ref VVOutput 模块中的 "生成图形输出 "部分。最后，对于高阶元素来说，在图形输出中每个单元只显示一个双线性四边形似乎不合适。因此，我们生成大小为(度数+1)x(度数+1)的斑块来捕捉解决方案的全部信息内容。有关这方面的更多信息，请参见 step-7 的教程程序。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::output_results() const 
 *   { 
 *     std::vector<std::string> solution_names(dim, "u"); 
 *     solution_names.emplace_back("p"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       interpretation(dim, 
 *                      DataComponentInterpretation::component_is_part_of_vector); 
 *     interpretation.push_back(DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.add_data_vector(dof_handler, 
 *                              solution, 
 *                              solution_names, 
 *                              interpretation); 
 * 
 *     data_out.build_patches(degree + 1); 
 * 
 *     std::ofstream output("solution.vtu"); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="MixedLaplacerun"></a> 
 * <h4>MixedLaplace::run</h4>
 * 

 * 
 * 这是我们主类的最后一个函数。它唯一的工作是按照自然顺序调用其他函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void MixedLaplaceProblem<dim>::run() 
 *   { 
 *     make_grid_and_dofs(); 
 *     assemble_system(); 
 *     solve(); 
 *     compute_errors(); 
 *     output_results(); 
 *   } 
 * } // namespace Step20 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 我们从  step-6  而不是  step-4  那里偷来的主函数。它几乎等同于 step-6 中的函数（当然，除了改变的类名），唯一的例外是我们将有限元空间的度数传递给混合拉普拉斯问题的构造函数（这里，我们使用零阶元素）。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step20; 
 * 
 *       const unsigned int     fe_degree = 0; 
 *       MixedLaplaceProblem<2> mixed_laplace_problem(fe_degree); 
 *       mixed_laplace_problem.run(); 
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
 * 
 * 
 * @endcode
examples/step-20/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Outputoftheprogramandgraphicalvisualization"></a><h3>Output of the program and graphical visualization</h3>



如果我们按原样运行程序，对于我们使用的 $32\times 32$ 网格，我们得到这样的输出（因为我们使用片状常数，所以总共有1024个单元，1024个压力自由度，2112个速度，因为Raviart-Thomas元素定义了每个面的一个自由度，有 $1024 + 32 = 1056$ 个面与 $x$  -轴平行，有同样数量的面与 $y$  -轴平行）。

@verbatim
\$ make run
[ 66%] Built target \step-20
Scanning dependencies of target run
[100%] Run \step-20 with Release configuration
Number of active cells: 1024
Total number of cells: 1365
Number of degrees of freedom: 3136 (2112+1024)
24 CG Schur complement iterations to obtain convergence.
Errors: ||e_p||_L2 = 0.0445032,   ||e_u||_L2 = 0.010826
[100%] Built target run
@endverbatim



当然，迭代次数如此之少的事实是由于我们开发的良好（但昂贵！）的预处理程序。为了获得对解决方案的信心，让我们看一下它。下面三张图片显示了（从左到右）X-速度、Y-速度和压力。

 <table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u_new.jpg" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.v_new.jpg" width="400" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.p_new.jpg" width="400" alt=""></td>
  </tr>
</table> 




让我们从压力开始：它在左边是最高的，在右边是最低的，所以流动将是从左到右。此外，虽然在图中几乎看不出来，但我们选择的压力场是这样的：从左到右的流动首先是向中心流动，然后再向外流动。因此，X-速度必须增加以使流动通过狭窄的部分，这一点在左图中很容易看到。中间的图像表示在域的左端有Y方向的内流，而在域的右端有Y方向的外流。




作为补充说明，请注意左图中的x-速度在x方向上是连续的，而y-速度在y方向上是连续的。其他方向上的流场是不连续的。这非常明显地反映了Raviart-Thomas元素的连续性特性，事实上，它只在空间H(div)而不是在空间 $H^1$ 。最后，压力场是完全不连续的，但鉴于我们选择了 <code>FE_DGQ(0)</code> 作为该解分量的有限元，这不应该令人惊讶。




<a name="Convergence"></a><h3>Convergence</h3>



该程序提供了两个明显的玩耍和观察收敛的地方：使用的有限元的程度（传递给 <code>MixedLaplaceProblem</code> class from <code>main()</code> 的构造器），和细化水平（在 <code>MixedLaplaceProblem::make_grid_and_dofs</code> 中确定）。人们可以做的是改变这些值，观察后来在程序运行过程中计算出的误差。




如果这样做，就会发现压力变量中 $L_2$ 的错误有如下模式。   <table align="center" class="doxtable">
  <tr>
    <th></th>
    <th colspan="3" align="center">Finite element order</th>
  </tr>
  <tr>
    <th>Refinement level</th>
    <th>0</th>
    <th>1</th>
    <th>2</th>
  </tr>
  <tr>
    <th>0</th>  <td>1.45344</td>  <td>0.0831743</td>  <td>0.0235186</td>
  </tr>
  <tr>
    <th>1</th>  <td>0.715099</td>  <td>0.0245341</td>  <td>0.00293983</td>
  </tr>
  <tr>
    <th>2</th>  <td>0.356383</td>  <td>0.0063458</td>  <td>0.000367478</td>
  </tr>
  <tr>
    <th>3</th>  <td>0.178055</td>  <td>0.00159944</td>  <td>4.59349e-05</td>
  </tr>
  <tr>
    <th>4</th>  <td>0.0890105</td>  <td>0.000400669</td>  <td>5.74184e-06</td>
  </tr>
  <tr>
    <th>5</th>  <td>0.0445032</td>  <td>0.000100218</td>  <td>7.17799e-07</td>
  </tr>
  <tr>
    <th>6</th>  <td>0.0222513</td>  <td>2.50576e-05</td>  <td>9.0164e-08</td>
  </tr>
  <tr>
    <th></th>  <th>$O(h)$</th>  <th>$O(h^2)$</th>  <th>$O(h^3)$</th>
  </tr>
</table> 

理论上预期的收敛顺序很好地反映在表中最后一行所显示的实验观察结果中。




我们可以用速度变量的 $L_2$ 误差做同样的实验。   <table align="center" class="doxtable">
  <tr>
    <th></th>
    <th colspan="3" align="center">Finite element order</th>
  </tr>
  <tr>
    <th>Refinement level</th>
    <th>0</th>
    <th>1</th>
    <th>2</th>
  </tr>
  <tr>
    <th>0</th> <td>0.367423</td> <td>0.127657</td> <td>5.10388e-14</td>
  </tr>
  <tr>
    <th>1</th> <td>0.175891</td> <td>0.0319142</td> <td>9.04414e-15</td>
  </tr>
  <tr>
    <th>2</th> <td>0.0869402</td> <td>0.00797856</td> <td>1.23723e-14</td>
  </tr>
  <tr>
    <th>3</th> <td>0.0433435</td> <td>0.00199464</td> <td>1.86345e-07</td>
  </tr>
  <tr>
    <th>4</th> <td>0.0216559</td> <td>0.00049866</td> <td>2.72566e-07</td>
  </tr>
  <tr>
    <th>5</th> <td>0.010826</td> <td>0.000124664</td> <td>3.57141e-07</td>
  </tr>
  <tr>
    <th>6</th> <td>0.00541274</td> <td>3.1166e-05</td> <td>4.46124e-07</td>
  </tr>
  <tr>
    <th></th>  <td>$O(h)$</td>  <td>$O(h^2)$</td>  <td>$O(h^3)$</td>
  </tr>
</table>  这里关于收敛顺序的结果是一样的。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Morerealisticpermeabilityfields"></a><h4>More realistic permeability fields</h4>


用于地下水或油藏模拟的现实流动计算不会使用恒定的渗透率。下面是改变这种情况的第一个相当简单的方法：我们使用一个渗透率，它在远离中心流线的地方迅速衰减，直到达到一个背景值0.001。这是为了模仿流体在砂岩中的行为：在大部分领域中，砂岩是均匀的，虽然对流体有渗透性，但不是过分的渗透；在另一块石头上，石头沿着一条线裂开了，或者说断层了，流体沿着这条大裂缝更容易流动。下面是我们如何实现这样的东西。

@code
template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));


  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();


      const double distance_to_flowline
        = std::fabs(points[p][1]-0.2*std::sin(10*points[p][0]));


      const double permeability = std::max(std::exp(-(distance_to_flowline*
                                                      distance_to_flowline)
                                                    / (0.1 * 0.1)),
                                           0.001);


      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 1./permeability;
    }
}
@endcode

记住，该函数返回渗透率张量的逆值。




通过一个明显更高的网格分辨率，我们可以直观地看到这一点，这里有x-和y-速度。

 <table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u-wiggle.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.v-wiggle.png" alt=""></td>
  </tr>
</table> 

很明显，流体基本上只沿着中线流动，而不是其他地方。




另一种可能性是使用一个随机的渗透率场。实现这一点的一个简单方法是在领域周围分散一些中心，然后使用一个渗透率场，这个渗透率场是这些中心的（负）指数之和。然后，流动将试图从一个高渗透率的中心跳到下一个中心。这是描述随机介质的一种完全不科学的尝试，但实现这种行为的一种可能性是这样的。

@code
template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    KInverse ();


    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<Tensor<2,dim> >    &values) const;


  private:
    std::vector<Point<dim> > centers;
};



template <int dim>
KInverse<dim>::KInverse ()
{
  const unsigned int N = 40;
  centers.resize (N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int d=0; d<dim; ++d)
      centers[i][d] = 2.*rand()/RAND_MAX-1;
}



template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));


  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();


      double permeability = 0;
      for (unsigned int i=0; i<centers.size(); ++i)
        permeability += std::exp(-(points[p] - centers[i]).norm_square() / (0.1 * 0.1));


      const double normalized_permeability
        = std::max(permeability, 0.005);


      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 1./normalized_permeability;
    }
}
@endcode



这个张量的对角线元素的片状常数插值（即 <code>normalized_permeability</code> ）看起来如下。

 <img src="https://www.dealii.org/images/steps/developer/step-20.k-random.png" alt=""> 


有了这样一个渗透率场，我们将得到如下的x-velocities和压力。

 <table style="width:60%" align="center">
  <tr>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.u-random.png" alt=""></td>
    <td><img src="https://www.dealii.org/images/steps/developer/step-20.p-random.png" alt=""></td>
  </tr>
</table> 

我们将在步骤21和步骤43中再次使用这些渗透率场。




<a name="Betterlinearsolvers"></a><h4>Better linear solvers</h4>


正如介绍中提到的，这里使用的Schur补码求解器并不是可以想象的最好的（也不打算成为一个特别好的）。更好的解算器可以在文献中找到，并且可以使用这里介绍的相同的块矩阵技术来构建。我们将在第22步中再次讨论这个主题，在这里我们首先为斯托克斯方程建立一个Schur补数求解器，然后在<a
href="step_22.html#improved-solver">Improved Solvers</a>部分讨论基于整体求解系统但基于单个块的预处理的更好方法。我们还将在第43步中再来讨论这个问题。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-20.cc"
*/
