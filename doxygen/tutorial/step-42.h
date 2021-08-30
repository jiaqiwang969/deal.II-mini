/**
@page step_42 The step-42 tutorial program
This tutorial depends on step-41, step-40.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Introduction">Introduction</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Reformulationasavariationalinequality">Reformulation as a variational inequality</a>
        <li><a href="#ANewtonmethodfortheplasticnonlinearity">A Newton method for the plastic nonlinearity</a>
        <li><a href="#ActiveSetmethodstosolvethesaddlepointproblem">Active Set methods to solve the saddle point problem</a>
        <li><a href="#Overallalgorithm">Overall algorithm</a>
        <li><a href="#Adaptivemeshrefinement">Adaptive mesh refinement</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeConstitutiveLawcodeclasstemplate">The <code>ConstitutiveLaw</code> class template</a>
      <ul>
        <li><a href="#ConstitutiveLawget_stress_strain_tensor">ConstitutiveLaw::get_stress_strain_tensor</a>
        <li><a href="#ConstitutiveLawget_linearized_stress_strain_tensors">ConstitutiveLaw::get_linearized_stress_strain_tensors</a>
        <li><a href="#ThecodeSphereObstaclecodeclass">The <code>SphereObstacle</code> class</a>
        <li><a href="#ThecodeBitmapFilecodeandcodeChineseObstaclecodeclasses">The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes</a>
      </ul>
        <li><a href="#ThecodePlasticityContactProblemcodeclasstemplate">The <code>PlasticityContactProblem</code> class template</a>
        <li><a href="#ImplementationofthecodePlasticityContactProblemcodeclass">Implementation of the <code>PlasticityContactProblem</code> class</a>
      <ul>
        <li><a href="#PlasticityContactProblemdeclare_parameters">PlasticityContactProblem::declare_parameters</a>
        <li><a href="#ThecodePlasticityContactProblemcodeconstructor">The <code>PlasticityContactProblem</code> constructor</a>
        <li><a href="#PlasticityContactProblemmake_grid">PlasticityContactProblem::make_grid</a>
        <li><a href="#PlasticityContactProblemsetup_system">PlasticityContactProblem::setup_system</a>
        <li><a href="#PlasticityContactProblemcompute_dirichlet_constraints">PlasticityContactProblem::compute_dirichlet_constraints</a>
        <li><a href="#PlasticityContactProblemassemble_mass_matrix_diagonal">PlasticityContactProblem::assemble_mass_matrix_diagonal</a>
        <li><a href="#PlasticityContactProblemupdate_solution_and_constraints">PlasticityContactProblem::update_solution_and_constraints</a>
        <li><a href="#PlasticityContactProblemassemble_newton_system">PlasticityContactProblem::assemble_newton_system</a>
        <li><a href="#PlasticityContactProblemcompute_nonlinear_residual">PlasticityContactProblem::compute_nonlinear_residual</a>
        <li><a href="#PlasticityContactProblemsolve_newton_system">PlasticityContactProblem::solve_newton_system</a>
        <li><a href="#PlasticityContactProblemsolve_newton">PlasticityContactProblem::solve_newton</a>
        <li><a href="#PlasticityContactProblemrefine_grid">PlasticityContactProblem::refine_grid</a>
        <li><a href="#PlasticityContactProblemmove_mesh">PlasticityContactProblem::move_mesh</a>
        <li><a href="#PlasticityContactProblemoutput_results">PlasticityContactProblem::output_results</a>
        <li><a href="#PlasticityContactProblemoutput_contact_force">PlasticityContactProblem::output_contact_force</a>
        <li><a href="#PlasticityContactProblemrun">PlasticityContactProblem::run</a>
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
examples/step-42/doc/intro.dox

 <br> 

<i>This program was contributed by Jörg Frohne (University of Siegen,
Germany) while on a long-term visit to Texas A&amp;M University, with significant
contributions by Timo Heister and Wolfgang Bangerth.
<br>
<br>
The code described here provides the basis for the numerical experiments shown
in the following paper:
<br>
  J. Frohne, T. Heister, W. Bangerth: <b>Efficient numerical methods for the large-scale, parallel
                  solution of elastoplastic contact problems</b><b>Efficient numerical methods for the large-scale, parallel
                  solution of elastoplastic contact problems</b>.
  Accepted for publication in International Journal for Numerical Methods in Engineering, 2015.
</i> 。




<a name="Intro"></a>

<a name="Introduction"></a><h3>Introduction</h3>


这个例子是第41步的延伸，考虑的是三维接触问题，具有各向同性硬化的弹塑性材料行为。换句话说，它考虑的是，如果把一个刚性的障碍物推到一个三维体上，它是如何变形的（接触问题），其中的变形受弹塑性材料法则（一种只能容纳一定最大应力的材料）的制约，随着变形的累积，该材料会硬化。为了说明我们打算做什么，在讨论太多细节之前，让我们只展示一张解决方案的图片（可变形体是一个立方体--实际上只显示了一半--，障碍物对应于一个汉字，将在下面讨论）。

 <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionLi2.png" alt=""> 


这个问题的描述意味着，与第41步相比，我们必须照顾到一个额外的非线性因素：材料行为。由于我们在这里考虑的是一个三维问题，我们还必须考虑到一个事实，即现在接触区是在可变形体的边界，而不是在内部。最后，与第41步相比，我们还必须在处理线性系统和不等式约束时处理悬空节点，因为我们希望使用自适应网格；在后一种情况下，我们将不得不处理优先考虑悬空节点的约束还是不等式的约束更重要。

由于在三维空间中很容易达到几百万个自由度，即使使用自适应网格细化，我们决定使用Trilinos和p4est来并行运行我们的代码，在步骤40的框架上进行并行化。并行化的其他指针可以在步骤32中找到。




<a name="Classicalformulation"></a><h3>Classical formulation</h3>


该问题的经典表述具有以下形式。

@f{align*}
 \varepsilon(\mathbf u) &= A\sigma + \varepsilon^p & &\quad\text{in } \Omega,\\


  -\textrm{div}\ \sigma &= \mathbf f & &\quad\text{in } \Omega,\\
  \varepsilon^p:(\tau - \sigma) &\geq 0\quad\forall\tau\text{ with
  }\mathcal{F}(\tau)\leq 0 & &\quad\text{in } \Omega,\\
  \mathbf u &= 0 & &\quad\text{on }\Gamma_D,\\
  \sigma \mathbf n - [\mathbf n \cdot(\sigma \mathbf n)]\mathbf n &= 0,
  \quad \mathbf n \cdot (\sigma
  \mathbf n) \leq 0 & &\quad\text{on }\Gamma_C,\\
  (\mathbf n \cdot (\sigma
  \mathbf n))(\mathbf n \cdot \mathbf u - g) &= 0,\quad \mathbf n
  \cdot \mathbf u - g \leq 0 & &\quad\text{on } \Gamma_C.


@f}

这里，这些方程的第一个定义了应变 $\varepsilon(\mathbf u)=\frac{1}{2}\left(\nabla \mathbf u
  + \nabla \mathbf u^T\right)$ 和应力 $\sigma$ 之间的关系，通过四阶顺应性张量 $A$ ； $\varepsilon^p$ 提供了应变的塑性成分，确保应力不超过屈服应力。我们将只考虑各向同性的材料，对于这些材料， $A$ 可以用Lam&eacute;模量 $\lambda$ 和 $\mu$ 表示，或者用体模量 $\kappa$ 和 $\mu$ 表示。第二个方程是力的平衡；我们在此不考虑任何体力，并假定 $\mathbf f=0$  。第三行的互补条件意味着，如果 $\mathcal{F}(\sigma)< 0$ ，则 $\varepsilon^p=0$ ，但当且仅当 $\mathcal{F}(\sigma) = 0$ ， $\varepsilon^p$ 可能是一个非零张量，特别是在这种情况下， $\varepsilon^p$ 必须指向 $\partial
\mathcal{F}(\sigma)/\partial \sigma$ 的方向。不等式 $\mathcal{F}(\sigma)\le 0$ 是塑性材料只能支持有限的应力；换句话说，如果外力会导致 $\sigma$ 的应力，那么它们就会产生塑性变形 $\varepsilon^p$ 的反应。这种<i>yield function</i>的典型形式是 $\mathcal{F}(\sigma)=|\sigma^D|-\sigma_{\text{yield}}$ ，其中 $\tau^D
= \tau - \dfrac{1}{3}tr(\tau)I$ 是张量的偏离部分， $|\cdot|$ 表示弗罗本尼斯规范。

进一步的方程描述了 $\Gamma_D$ 上固定的零位移，在可能出现接触的表面 $\Gamma_C=\partial\Omega\backslash\Gamma_D$ 上，障碍物施加的法向力 $\sigma_n=\mathbf n \cdot (\sigma(\mathbf u)
  \mathbf n)$ 是向内的（障碍物对我们的身体没有 "拉力"），切向分量为零 $\mathbf \sigma_t= \sigma \mathbf n - \mathbf \sigma_n \mathbf n
= \sigma \mathbf n - [\mathbf n \cdot(\sigma \mathbf n)]\mathbf n$  。最后一个条件又是一个互补条件，意味着在 $\Gamma_C$ 上，只有当身体与障碍物接触时，法向力才能非零；第二部分描述了障碍物和身体的不可穿透性。最后两个方程通常被称为Signorini接触条件。

大多数材料--尤其是金属--都有这样的特性，即它们在变形时表现出一定的硬化。换句话说， $\sigma_{\text{yield}}$ 随着变形而增加。在实践中，导致硬化的不是弹性变形，而是塑性成分。有不同的构成法则来描述这些材料行为。最简单的称为线性各向同性硬化，由流动函数  $\mathcal{F}(\sigma,\varepsilon^p) = \vert\sigma^D\vert - (\sigma_0 +
\gamma^{\text{iso}}|\varepsilon^p|)$  描述。




<a name="Reformulationasavariationalinequality"></a><h3>Reformulation as a variational inequality</h3>


一般来说，处理不等式是相当笨拙的。在这里，我们必须处理两个问题：塑性和接触问题。正如本页顶部提到的论文中详细描述的那样，我们至少可以重新表述塑性，使其看起来像一个非线性，然后我们可以用牛顿方法处理。这在数学上略显棘手，因为非线性不只是一些平滑的函数，而是在应力达到屈服应力的地方有结点；然而，对于这样的<i>semismooth</i>函数，可以证明牛顿方法仍然收敛。

在不涉及细节的情况下，我们也将摆脱作为独立变量的应力，而完全用位移来工作  $\mathbf u$  。最终，这种重构的目标是，我们希望最终得到一个对称的、正定的问题--比如一个线性化的弹性问题，其空间变量系数由塑性行为产生--需要在每个牛顿步骤中解决。我们希望如此，因为有高效和可扩展的方法来解决这样的线性系统，如用代数多重网格的CG预处理。这与我们继续使用包含位移和应力的混合公式所得到的类似于混合拉普拉斯的鞍点问题（见第20步）是相反的，第20步已经提示了构建良好的求解器和预处理器是多么困难。

说到这里，让我们简单陈述一下我们在重构后得到的问题（同样，细节可以在论文中找到）。找到一个位移 $\mathbf u \in
V^+$ ，以便

@f{align*}
\left(P_{\Pi}(C\varepsilon(\mathbf u)),\varepsilon(\varphi) - \varepsilon(\mathbf u)\right) \geq 0,\quad \forall \varphi\in V^+.


@f}

其中投影仪 $P_\Pi$ 被定义为

@f{align*}
 P_{\Pi}(\tau) \dealcoloneq \begin{cases}
    \tau, & \text{if }\vert\tau^D\vert \leq \sigma_0,\\
    \left[
      \dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}} +
      \left(1-\dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}\right)\dfrac{\sigma_0}{\vert\tau^D\vert}
    \right]\tau^D
    + \dfrac{1}{3}\text{trace}(\tau) I, & \text{if }\vert\tau^D\vert >
    \sigma_0,
  \end{cases}


@f}

和空间 $V^+$ 是满足接触条件的所有位移的空间。

@f{align*}
  V
  &=
  \left\{ \mathbf u\in \left[H^1(\Omega)\right]^{d}:
    \mathbf u = 0 \text{ on } \Gamma_D\right\},
  \\
  V^+
  &=
  \left\{ \mathbf u\in V: \mathbf n \cdot \mathbf u\leq g \text{ on } \Gamma_C \right\}.


@f}



在实际代码中，我们将使用缩写  $\gamma=\dfrac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}$  。

鉴于这种表述，我们将应用两种技术。

- 运行牛顿方法来迭代出投影仪的非线性。

- 为接触条件运行一个主动设置方法，方法与我们在步骤41中所做的基本相同。

一个严格的方法是在我们迭代牛顿方法到收敛时保持活动集的固定（或者也许反过来：在进入下一个牛顿迭代之前找到最终的活动集）。在实践中，事实证明，每个活动集迭代只做一个牛顿步骤就足够了，所以我们将同时迭代它们。我们还将每隔一段时间细化一下网格。




<a name="ANewtonmethodfortheplasticnonlinearity"></a><h3>A Newton method for the plastic nonlinearity</h3>


如前所述，我们将通过应用牛顿方法来处理算子 $P_\Pi$ 的非线性，尽管该算子在严格意义上是不可微的。然而，它满足了<i>slant</i>的可微条件，这就足以使牛顿方法发挥作用。由此产生的方法被称为<i>semi-smooth Newton method</i>，听起来令人印象深刻，但实际上只是一个牛顿方法应用于一个具有适当选择的 "导数 "的半光滑函数。

在目前的情况下，我们将通过在每个迭代 $i$ 中求解以下方程来运行我们的迭代（仍然是不等式，但是线性化）。

@f{align*}
  \label{eq:linearization}
  \left(I_{\Pi}\varepsilon(\tilde {\mathbf u}^{i}),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right) \geq
  \left(\left(I_{\Pi}\varepsilon({\mathbf u}^{i-1}),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right) -
  \left(P_{\Pi}(C\varepsilon({\mathbf u}^{i-1})),
    \varepsilon(\varphi) - \varepsilon(\tilde {\mathbf u}^{i})\right)\right),
  \quad \forall \varphi\in V^+,


@f}

其中，等级4张量 $I_\Pi=I_\Pi(\varepsilon^D(\mathbf u^{i-1}))$ 由以下公式给出

@f{align}
  I_\Pi = \begin{cases}
    C_{\mu} + C_{\kappa}, & \hspace{-8em} \text{if } \vert C\varepsilon^D(\mathbf u^{i-1}) \vert \leq \sigma_0,
    \\
    \frac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}} C_{\mu} + \frac{\left(1-\frac{\gamma^{\text{iso}}}{2\mu + \gamma^{\text{iso}}}\right)\sigma_0}{\vert C\varepsilon^D(\mathbf u^{i-1}) \vert}\left(C_{\mu} -
      2\mu\dfrac{C\varepsilon^D(\mathbf u^{i-1})\otimes C\varepsilon^D(\mathbf
        u^{i-1})}{\vert C\varepsilon^D(\mathbf u^{i-1})\vert^2}\right) + C_{\kappa}, & \text{ else.}
\end{cases}


@f}

这个张量是 $P_\Pi(C\cdot)$ 围绕 $\varepsilon^D(\mathbf u^{i-1})$ 的（形式）线性化。对于我们这里考虑的线性各向同性材料，投影仪的体积和剪切分量由以下公式给出

@f{gather*}
  C_{\kappa} = \kappa I\otimes I,
  \qquad\qquad\qquad\qquad
  C_{\mu} = 2\mu\left(\mathbb{I}  - \dfrac{1}{3} I\otimes
    I\right),


@f}

其中 $I$ 和 $\mathbb{I}$ 分别是等级为2和4的认同张量。

请注意，这个问题对应于线性弹性接触问题，其中 $I_\Pi$ 扮演弹性张量的角色  $C=A^{-1}$  。事实上，如果材料在某一点上没有塑性，那么 $I_\Pi=C$  。然而，在材料具有塑性的地方， $I_\Pi$ 是一个空间变化的函数。在任何情况下，我们必须解决牛顿迭代的系统 $\tilde {\mathbf u}^{i}$ 使我们更接近重写我们问题的目标，使我们能够使用众所周知的椭圆系统的求解器和预处理器。

作为对牛顿方法的最后说明，让我们提一下，正如牛顿方法常见的那样，我们需要通过控制步长来使其全球化。换句话说，虽然上面的系统求解的是 $\tilde {\mathbf u}^{i}$ ，但最后的迭代结果将是

@f{align*}
  {\mathbf u}^{i} = {\mathbf u}^{i-1} + \alpha_i (\tilde {\mathbf u}^{i} - {\mathbf u}^{i-1})


@f}

其中右边括号中的差值扮演了传统牛顿方向的角色， $\delta {\mathbf u}^{i}$  。我们将用标准的直线搜索来确定 $\alpha^i$ 。




<a name="ActiveSetmethodstosolvethesaddlepointproblem"></a><h3>Active Set methods to solve the saddle point problem</h3>


这个要在每个牛顿步骤中解决的线性化问题基本上与步骤41一样。唯一的区别在于接触区是在边界而不是在域中。但这没有进一步的后果，所以我们参考步骤41的文件，唯一的提示是 $\mathcal{S}$ 这次包含了接触边界的所有顶点 $\Gamma_C$ 。和那里一样，我们需要做的是保持一个自由度子集的固定，导致额外的约束，可以写成一个鞍点问题。然而，正如论文中所讨论的，通过以适当的方式写这些约束，消除自由度之间的耦合，我们最终会得到一组节点，这些节点基本上只是附加了Dirichlet值。




<a name="Overallalgorithm"></a><h3>Overall algorithm</h3>


上述算法结合了阻尼半光滑牛顿法（我们用于非线性构成法）和半光滑牛顿法用于接触。它的工作原理如下。<ol>  <li>  初始化活动和非活动集 $\mathcal{A}_i$ 和 $\mathcal{F}_i$ ，使 $\mathcal{S} = \mathcal{A}_i \cup \mathcal{F}_i$ 和 $\mathcal{A}_i \cap
 \mathcal{F}_i = \emptyset$ 和集 $i = 1$  。这里， $\mathcal{S}$ 是位于可能发生接触的域的表面的所有自由度的集合。  起始值 $\hat U^0 \dealcoloneq
 P_{\mathcal{A}_k}(0)$ 满足我们的障碍条件，也就是说，我们将初始零位移投射到可行位移集合上。

   <li>  组装牛顿矩阵  $A_{pq} \dealcoloneq a'(
 U^{i-1};\varphi_p,\varphi_q)$  和右侧  $F(\hat U^{i-1})$  。  这些对应于线性化的牛顿步骤，暂时忽略了接触不等式。

   <li>  找到满足@f{align*}
 A\tilde U^i + B\Lambda^i & = F, &\\
 \left[B^T\tilde U^i\right]_p & = G_p & \forall p\in\mathcal{A}_i,\\
 \Lambda^i_p & = 0 & \forall p\in\mathcal{F}_i.
 @f}的原始-双数对 $(\tilde U^i,\Lambda^i)$  。

如同步骤-41，我们可以通过消除第一个方程中 ${\cal A}_i$ 的那些自由度来获得这个问题的解决方案，并获得一个线性系统 $\hat {\hat A}(U^{i-1}) \tilde U^i = \hat {\hat H}(U^{i-1})$  。




   <li>  通过应用直线搜索和计算 $U^{i-1}$ 和 $\tilde U^i$ 的线性组合来减弱 $i>2$ 的牛顿迭代。这需要找到一个 $\alpha^i_l \dealcoloneq 2^{-l},(l=0,\ldots,10)$ ，以便@f{gather*}U^i \dealcoloneq \alpha^i_l\bar U^i +
 (1-\alpha^i_l)U^{i-1}@f}。

满足@f{gather*}
   \vert {\hat R}\left({\mathbf u}^{i}\right) \vert < \vert {\hat R}\left({\mathbf u}^{i-1}\right) \vert.
 \f}与 ${\hat R}\left({\mathbf u}\right)=\left(P_{Pi}(C\varepsilon(u)),\varepsilon(\varphi^{i}_p\right)$ ，除了(i)元素 $p\in\mathcal{A}_i$ ，我们设置 ${\hat R}\left({\mathbf u}\right)=0$ ，和(ii)对应于悬挂节点的元素，我们以通常方式消除。

   <li>  通过@f{gather*}\mathcal{A}_{i+1} \dealcoloneq \lbrace p\in\mathcal{S}:\Lambda^i_p +
 c\left(\left[B^TU^i\right]_p - G_p\right) > 0\rbrace,@f}定义新的活动和非活动集。

@f{gather*}\mathcal{F}_{i+1} \dealcoloneq \lbrace p\in\mathcal{S}:\Lambda^i_p +
 c\left(\left[B^TU^i\right]_p - G_p\right) \leq 0\rbrace.@f}



   <li> 项目 $U^i$ ，使其满足接触不等式，@f{gather*}\hat U^i \dealcoloneq P_{\mathcal{A}_{i+1}}(U^i).@f} 。

这里， $P_{\mathcal{A}}(U)$ 是 $\mathcal{A}$ 中的活性成分对间隙@f{gather*}P_{\mathcal{A}}(U)_p \dealcoloneq \begin{cases}
 U_p, & \textrm{if}\quad p\notin\mathcal{A}\\
 g_{h,p}, & \textrm{if}\quad
 p\in\mathcal{A},
 \end{cases}@f}的投影。

其中 $g_{h,p}$ 是<i>gap</i>，表示障碍物与身体未位移配置的距离。

   <li>  如果 $\mathcal{A}_{i+1} = \mathcal{A}_k$ 和 $\left\|
 {\hat R}\left({\mathbf u}^{i}\right) \right\|_{\ell_2} < \delta$ 则停止，否则设置 $i=i+1$ 并转到步骤（1）。这一步确保我们只有在找到正确的活动集和塑性已经迭代到足够的精度时才停止迭代。   </ol> 

在这个算法的第3步中，矩阵 $B\in\mathbb{R}^{n\times m}$ ,  $n>m$ 描述了位移和拉格朗日乘数（接触力）的基数的耦合，在我们的情况下它不是二次的，因为 $\Lambda^k$ 只定义在 $\Gamma_C$ ，即可能发生接触的面。如文中所示，我们可以选择 $B$ 是一个每行只有一个条目的矩阵，（另见H&uuml;eber, Wohlmuth:A primal-dual active set strategy for non-linear multibody contact problems, Comput.Method Appl. Mech.Engrg.194, 2005, pp.3147-3166）。)矢量 $G$ 是由间隙 $g_h$ 的合适近似值定义的。

@f{gather*}G_p = \begin{cases}
g_{h,p}, & \text{if}\quad p\in\mathcal{S}\\
0, & \text{if}\quad p\notin\mathcal{S}.
\end{cases}@f}






<a name="Adaptivemeshrefinement"></a><h3>Adaptive mesh refinement</h3>


由于我们的程序是在三维空间中运行的，所以程序执行的计算很昂贵。因此，使用自适应网格细化是在可接受的运行时间内的一个重要步骤。为了使我们的生活更轻松，我们简单地选择已经在deal.II中实现的KellyErrorEstimator。我们把包含位移 $u$ 的解向量交给它。正如我们将在结果中看到的，它产生了一个相当合理的接触区和塑性的自适应网格。




<a name="Implementation"></a><h3>Implementation</h3>


本教程实质上是步骤40和步骤41的混合体，但我们没有使用PETSc，而是让Trilinos库来处理线性代数的并行化问题（就像步骤32一样）。由于我们试图解决一个类似于步骤41的问题，我们将使用同样的方法，但现在是并行的。

一个困难是处理来自Dirichlet条件的约束，悬挂节点和由接触产生的不平等条件。为此，我们创建了三个AffineConstraints类型的对象，它们描述了各种约束条件，我们将在每次迭代中适当地组合它们。

与第41步相比，该计划有一些新的课程。

 <ul>   <li>   <code>ConstitutiveLaw</code>  描述材料的塑性行为。

 <li>   <code>SphereObstacle</code> 描述一个球体，作为被推入可变形弹性体的障碍物。   是用这个还是下一个类来描述障碍物，由输入参数文件决定。

 <li>   <code>ChineseObstacle</code> （和一个辅助类）是一个允许我们从一个文件中读入障碍物的类。在我们将在结果部分展示的例子中，这个文件将是 <code>'obstacle_file.dat'</code> ，并对应于显示力或力量的中文、日文或韩文符号的数据（见http://www.orientaloutpost.com/："这个词可用于激励--它也可以指力量/运动/推进/力。它可以是任何使你继续前进的内部或外部事物。这是用中文表达动机的最安全方式。如果你的听众是日本人，请看另一个关于动机的条目。这是日语和韩语中的一个词，但它的意思是 "动力 "或 "动能"（没有你可能正在寻找的动机的意思）"）。实质上，我们将假装有一个印章（即对应于平底障碍物的面具，没有中间高度的碎片），我们把它压在身体里。有关的符号看起来如下（也可参见本节顶部的图片，了解最终结果是怎样的）。

    <img src="https://www.dealii.org/images/steps/developer/step-42.character.png" alt="" width="25%">   </ul> 。

除此以外，让我们只对以下方面进行评论。   <ul>   <li>  程序允许你通过参数文件从两个不同的粗略网格中进行选择。这些是立方体 $[0,1]^3$ 或半球体，其开放面朝向正 $z$ 方向。

 <li> 在这两种情况下，我们将假设可能与障碍物接触的边界部分具有边界指标一的惯例。对于这两种网格，我们假定这是一个自由表面，即身体要么在那里接触，要么没有力作用在它身上。对于半球体，弯曲部分的边界指标为零，我们在那里施加零位移。对于盒子，我们沿底部施加零位移，但允许沿边的垂直位移（尽管没有水平位移）。   </ul> 


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name=""></a> 
 * @sect3{Include files}  这组包含文件在这个时候已经没有什么惊喜了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/conditional_ostream.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/index_set.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/timer.h> 
 * 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparsity_tools.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/block_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_bicgstab.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_vector.h> 
 * #include <deal.II/lac/trilinos_parallel_block_vector.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * #include <deal.II/lac/trilinos_solver.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * 
 * #include <deal.II/distributed/tria.h> 
 * #include <deal.II/distributed/grid_refinement.h> 
 * #include <deal.II/distributed/solution_transfer.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * #include <deal.II/numerics/fe_field_function.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后，我们包括两个系统头文件，让我们为输出文件创建一个目录。第一个头文件提供了 <code>mkdir</code> 的功能，第二个头文件让我们确定在 <code>mkdir</code> 失败时发生了什么。
 * 

 * 
 * 
 * @code
 * #include <sys/stat.h> 
 * #include <cerrno> 
 * 
 * namespace Step42 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeConstitutiveLawcodeclasstemplate"></a> 
 * <h3>The <code>ConstitutiveLaw</code> class template</h3>
 * 

 * 
 * 该类提供了一个构成法的接口，即应变  $\varepsilon(\mathbf u)$  和应力  $\sigma$  之间的关系。在这个例子中，我们使用的是具有线性、各向同性硬化的弹塑性材料行为。这种材料的特点是杨氏模量  $E$  ，泊松比  $\nu$  ，初始屈服应力  $\sigma_0$  和各向同性硬化参数  $\gamma$  。 对于 $\gamma = 0$ ，我们得到完美的弹塑性行为。
 * 

 * 
 * 正如描述这个程序的论文所解释的那样，第一个牛顿步骤是用一个完全弹性材料模型来解决的，以避免同时处理两种非线性（塑性和接触）。为此，这个类有一个函数 <code>set_sigma_0()</code> ，我们在后面使用这个函数，简单地将 $\sigma_0$ 设置为一个非常大的值--基本上保证了实际应力不会超过它，从而产生一个弹性材料。当我们准备使用塑性模型时，我们使用相同的函数将 $\sigma_0$ 设置回其适当的值。 由于这种方法，我们需要将 <code>sigma_0</code> 作为这个类的唯一非静态成员变量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ConstitutiveLaw 
 *   { 
 *   public: 
 *     ConstitutiveLaw(const double E, 
 *                     const double nu, 
 *                     const double sigma_0, 
 *                     const double gamma); 
 * 
 *     void set_sigma_0(double sigma_zero); 
 * 
 *     bool get_stress_strain_tensor( 
 *       const SymmetricTensor<2, dim> &strain_tensor, 
 *       SymmetricTensor<4, dim> &      stress_strain_tensor) const; 
 * 
 *     void get_linearized_stress_strain_tensors( 
 *       const SymmetricTensor<2, dim> &strain_tensor, 
 *       SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
 *       SymmetricTensor<4, dim> &      stress_strain_tensor) const; 
 * 
 *   private: 
 *     const double kappa; 
 *     const double mu; 
 *     double       sigma_0; 
 *     const double gamma; 
 * 
 *     const SymmetricTensor<4, dim> stress_strain_tensor_kappa; 
 *     const SymmetricTensor<4, dim> stress_strain_tensor_mu; 
 *   }; 
 * 
 * @endcode
 * 
 * ConstitutiveLaw类的构造函数为我们的可变形体设置所需的材料参数。弹性各向同性介质的材料参数可以用多种方式定义，如一对 $E, \nu$ （弹性模量和泊松数），使用Lam&eacute;参数 $\lambda,mu$ 或其他几种常用的约定。在这里，构造器采用 $E,\nu$ 形式的材料参数描述，但由于这证明这些不是出现在塑性投影仪方程中的系数，我们立即将它们转换为更合适的体模和剪模集合 $\kappa,\mu$ 。 此外，构造器以 $\sigma_0$ （无任何塑性应变的屈服应力）和 $\gamma$ （硬化参数）作为参数。在这个构造函数中，我们还计算了应力-应变关系的两个主成分及其线性化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ConstitutiveLaw<dim>::ConstitutiveLaw(double E, 
 *                                         double nu, 
 *                                         double sigma_0, 
 *                                         double gamma) 
 *     : kappa(E / (3 * (1 - 2 * nu))) 
 *     , mu(E / (2 * (1 + nu))) 
 *     , sigma_0(sigma_0) 
 *     , gamma(gamma) 
 *     , stress_strain_tensor_kappa(kappa * 
 *                                  outer_product(unit_symmetric_tensor<dim>(), 
 *                                                unit_symmetric_tensor<dim>())) 
 *     , stress_strain_tensor_mu( 
 *         2 * mu * 
 *         (identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(), 
 *                                                 unit_symmetric_tensor<dim>()) / 
 *                                     3.0)) 
 *   {} 
 * 
 *   template <int dim> 
 *   void ConstitutiveLaw<dim>::set_sigma_0(double sigma_zero) 
 *   { 
 *     sigma_0 = sigma_zero; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConstitutiveLawget_stress_strain_tensor"></a> 
 * <h4>ConstitutiveLaw::get_stress_strain_tensor</h4>
 * 

 * 
 * 这是构成法则的主成分。它计算的是四阶对称张量，根据上面给出的投影，当在一个特定的应变点上评估时，该张量将应变与应力联系起来。我们需要这个函数来计算 <code>PlasticityContactProblem::residual_nl_system()</code> 中的非线性残差，我们将这个张量与正交点的应变相乘。计算遵循介绍中列出的公式。在比较那里的公式和下面的实现时，记得 $C_\mu : \varepsilon = \tau_D$ 和 $C_\kappa : \varepsilon = \kappa \text{trace}(\varepsilon) I = \frac 13 \text{trace}(\tau) I$  。
 * 

 * 
 * 该函数返回正交点是否是塑性的，以便在下游对有多少正交点是塑性的，有多少是弹性的进行一些统计。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   bool ConstitutiveLaw<dim>::get_stress_strain_tensor( 
 *     const SymmetricTensor<2, dim> &strain_tensor, 
 *     SymmetricTensor<4, dim> &      stress_strain_tensor) const 
 *   { 
 *     Assert(dim == 3, ExcNotImplemented()); 
 * 
 *     SymmetricTensor<2, dim> stress_tensor; 
 *     stress_tensor = 
 *       (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 
 * 
 *     const SymmetricTensor<2, dim> deviator_stress_tensor = 
 *       deviator(stress_tensor); 
 *     const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 
 * 
 *     stress_strain_tensor = stress_strain_tensor_mu; 
 *     if (deviator_stress_tensor_norm > sigma_0) 
 *       { 
 *         const double beta = sigma_0 / deviator_stress_tensor_norm; 
 *         stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
 *       } 
 * 
 *     stress_strain_tensor += stress_strain_tensor_kappa; 
 * 
 *     return (deviator_stress_tensor_norm > sigma_0); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ConstitutiveLawget_linearized_stress_strain_tensors"></a> 
 * <h4>ConstitutiveLaw::get_linearized_stress_strain_tensors</h4>
 * 

 * 
 * 该函数返回线性化的应力应变张量，围绕前一个牛顿步骤 $u^{i-1}$ 的解进行线性化  $i-1$  。 参数 <code>strain_tensor</code> （通常表示为 $\varepsilon(u^{i-1})$ ）必须作为参数传递，并作为线性化点。该函数在变量stress_strain_tensor中返回非线性构成法的导数，在stress_strain_tensor_linearized中返回线性化问题的应力-应变张量。 参见 PlasticityContactProblem::assemble_nl_system ，其中使用了这个函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ConstitutiveLaw<dim>::get_linearized_stress_strain_tensors( 
 *     const SymmetricTensor<2, dim> &strain_tensor, 
 *     SymmetricTensor<4, dim> &      stress_strain_tensor_linearized, 
 *     SymmetricTensor<4, dim> &      stress_strain_tensor) const 
 *   { 
 *     Assert(dim == 3, ExcNotImplemented()); 
 * 
 *     SymmetricTensor<2, dim> stress_tensor; 
 *     stress_tensor = 
 *       (stress_strain_tensor_kappa + stress_strain_tensor_mu) * strain_tensor; 
 * 
 *     stress_strain_tensor            = stress_strain_tensor_mu; 
 *     stress_strain_tensor_linearized = stress_strain_tensor_mu; 
 * 
 *     SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor); 
 *     const double deviator_stress_tensor_norm = deviator_stress_tensor.norm(); 
 * 
 *     if (deviator_stress_tensor_norm > sigma_0) 
 *       { 
 *         const double beta = sigma_0 / deviator_stress_tensor_norm; 
 *         stress_strain_tensor *= (gamma + (1 - gamma) * beta); 
 *         stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta); 
 *         deviator_stress_tensor /= deviator_stress_tensor_norm; 
 *         stress_strain_tensor_linearized -= 
 *           (1 - gamma) * beta * 2 * mu * 
 *           outer_product(deviator_stress_tensor, deviator_stress_tensor); 
 *       } 
 * 
 *     stress_strain_tensor += stress_strain_tensor_kappa; 
 *     stress_strain_tensor_linearized += stress_strain_tensor_kappa; 
 *   } 
 * @endcode
 * 
 * <h3>Equation data: boundary forces, boundary values, obstacles</h3>
 * 

 * 
 * 下面的内容应该是比较标准的。我们需要边界强迫项（我们在此选择为零）和不属于接触面的边界部分的边界值（在此也选择为零）的类。
 * 

 * 
 * 
 * @code
 *   namespace EquationData 
 *   { 
 *     template <int dim> 
 *     class BoundaryForce : public Function<dim> 
 *     { 
 *     public: 
 *       BoundaryForce(); 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  values) const override; 
 *     }; 
 * 
 *     template <int dim> 
 *     BoundaryForce<dim>::BoundaryForce() 
 *       : Function<dim>(dim) 
 *     {} 
 * 
 *     template <int dim> 
 *     double BoundaryForce<dim>::value(const Point<dim> &, 
 *                                      const unsigned int) const 
 *     { 
 *       return 0.; 
 *     } 
 * 
 *     template <int dim> 
 *     void BoundaryForce<dim>::vector_value(const Point<dim> &p, 
 *                                           Vector<double> &  values) const 
 *     { 
 *       for (unsigned int c = 0; c < this->n_components; ++c) 
 *         values(c) = BoundaryForce<dim>::value(p, c); 
 *     } 
 * 
 *     template <int dim> 
 *     class BoundaryValues : public Function<dim> 
 *     { 
 *     public: 
 *       BoundaryValues(); 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 *     }; 
 * 
 *     template <int dim> 
 *     BoundaryValues<dim>::BoundaryValues() 
 *       : Function<dim>(dim) 
 *     {} 
 * 
 *     template <int dim> 
 *     double BoundaryValues<dim>::value(const Point<dim> &, 
 *                                       const unsigned int) const 
 *     { 
 *       return 0.; 
 *     } 
 * 
 * @endcode
 * 
 * 
 * <a name="ThecodeSphereObstaclecodeclass"></a> 
 * <h4>The <code>SphereObstacle</code> class</h4>
 * 

 * 
 * 下面这个类是可以从输入文件中选择的两个障碍物中的第一个。它描述了一个以位置 $x=y=0.5, z=z_{\text{surface}}+0.59$ 和半径 $r=0.6$ 为中心的球体，其中 $z_{\text{surface}}$ 是可变形体的（平）表面的垂直位置。该函数的 <code>value</code> 返回给定 $x,y$ 值的障碍物位置，如果该点实际位于球体下方，则返回一个不可能干扰变形的大正值，如果它位于球体的 "阴影 "之外。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class SphereObstacle : public Function<dim> 
 *     { 
 *     public: 
 *       SphereObstacle(const double z_surface); 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  values) const override; 
 * 
 *     private: 
 *       const double z_surface; 
 *     }; 
 * 
 *     template <int dim> 
 *     SphereObstacle<dim>::SphereObstacle(const double z_surface) 
 *       : Function<dim>(dim) 
 *       , z_surface(z_surface) 
 *     {} 
 * 
 *     template <int dim> 
 *     double SphereObstacle<dim>::value(const Point<dim> & p, 
 *                                       const unsigned int component) const 
 *     { 
 *       if (component == 0) 
 *         return p(0); 
 *       else if (component == 1) 
 *         return p(1); 
 *       else if (component == 2) 
 *         { 
 *           if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5) < 0.36) 
 *             return (-std::sqrt(0.36 - (p(0) - 0.5) * (p(0) - 0.5) - 
 *                                (p(1) - 0.5) * (p(1) - 0.5)) + 
 *                     z_surface + 0.59); 
 *           else 
 *             return 1000; 
 *         } 
 * 
 *       Assert(false, ExcNotImplemented()); 
 *       return 1e9; // an unreasonable value; ignored in debug mode because of the 
 * 
 * @endcode
 * 
 * 前面的断言
 * 

 * 
 * 
 * @code
 *     } 
 * 
 *     template <int dim> 
 *     void SphereObstacle<dim>::vector_value(const Point<dim> &p, 
 *                                            Vector<double> &  values) const 
 *     { 
 *       for (unsigned int c = 0; c < this->n_components; ++c) 
 *         values(c) = SphereObstacle<dim>::value(p, c); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="ThecodeBitmapFilecodeandcodeChineseObstaclecodeclasses"></a> 
 * <h4>The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes</h4>
 * 

 * 
 * 下面两个类描述了介绍中概述的障碍物，即汉字。两个中的第一个， <code>BitmapFile</code> 负责从一个以pbm ascii格式存储的图片文件中读入数据。这个数据将被双线性插值，从而提供一个描述障碍物的函数。(下面的代码显示了如何通过在给定的数据点之间进行内插来构造一个函数。人们可以使用在这个教程程序写完后引入的 Functions::InterpolatedUniformGridData, ，它正是我们在这里想要的，但看看如何手工操作是有启发的）。)
 * 

 * 
 * 我们从文件中读取的数据将被存储在一个名为 obstacle_data 的双 std::vector 中。 这个向量构成了计算单片双线性函数的基础，作为一个多项式插值。我们将从文件中读取的数据由零（白色）和一（黑色）组成。
 * 

 * 
 * <code>hx,hy</code> 变量表示 $x$ 和 $y$ 方向的像素之间的间距。  <code>nx,ny</code> 是这些方向上的像素的数量。   <code>get_value()</code> 返回图像在给定位置的值，由相邻像素值插值而成。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class BitmapFile 
 *     { 
 *     public: 
 *       BitmapFile(const std::string &name); 
 * 
 *       double get_value(const double x, const double y) const; 
 * 
 *     private: 
 *       std::vector<double> obstacle_data; 
 *       double              hx, hy; 
 *       int                 nx, ny; 
 * 
 *       double get_pixel_value(const int i, const int j) const; 
 *     }; 
 * 
 * @endcode
 * 
 * 该类的构造函数从给定的文件名中读入描述障碍物的数据。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     BitmapFile<dim>::BitmapFile(const std::string &name) 
 *       : obstacle_data(0) 
 *       , hx(0) 
 *       , hy(0) 
 *       , nx(0) 
 *       , ny(0) 
 *     { 
 *       std::ifstream f(name); 
 *       AssertThrow(f, 
 *                   ExcMessage(std::string("Can't read from file <") + name + 
 *                              ">!")); 
 * 
 *       std::string temp; 
 *       f >> temp >> nx >> ny; 
 * 
 *       AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format.")); 
 * 
 *       for (int k = 0; k < nx * ny; ++k) 
 *         { 
 *           double val; 
 *           f >> val; 
 *           obstacle_data.push_back(val); 
 *         } 
 * 
 *       hx = 1.0 / (nx - 1); 
 *       hy = 1.0 / (ny - 1); 
 * 
 *       if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0) 
 *         std::cout << "Read obstacle from file <" << name << ">" << std::endl 
 *                   << "Resolution of the scanned obstacle picture: " << nx 
 *                   << " x " << ny << std::endl; 
 *     } 
 * 
 * @endcode
 * 
 * 下面两个函数返回坐标为 $i,j$ 的给定像素的值，我们将其与定义在位置 <code>i*hx, j*hy</code> 的函数值和任意坐标 $x,y$ 的函数值相识别，在这里我们对两个函数中第一个函数返回的点值进行双线性内插。在第二个函数中，对于每个 $x,y$ ，我们首先计算离 $x,y$ 左下方最近的像素坐标的（整数）位置，然后计算这个像素内的坐标 $\xi,\eta$ 。我们从下方和上方截断这两种变量，以避免在评估函数时超出其定义的范围而可能发生的舍入误差问题。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     double BitmapFile<dim>::get_pixel_value(const int i, const int j) const 
 *     { 
 *       assert(i >= 0 && i < nx); 
 *       assert(j >= 0 && j < ny); 
 *       return obstacle_data[nx * (ny - 1 - j) + i]; 
 *     } 
 * 
 *     template <int dim> 
 *     double BitmapFile<dim>::get_value(const double x, const double y) const 
 *     { 
 *       const int ix = std::min(std::max(static_cast<int>(x / hx), 0), nx - 2); 
 *       const int iy = std::min(std::max(static_cast<int>(y / hy), 0), ny - 2); 
 * 
 *       const double xi  = std::min(std::max((x - ix * hx) / hx, 1.), 0.); 
 *       const double eta = std::min(std::max((y - iy * hy) / hy, 1.), 0.); 
 * 
 *       return ((1 - xi) * (1 - eta) * get_pixel_value(ix, iy) + 
 *               xi * (1 - eta) * get_pixel_value(ix + 1, iy) + 
 *               (1 - xi) * eta * get_pixel_value(ix, iy + 1) + 
 *               xi * eta * get_pixel_value(ix + 1, iy + 1)); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，这是一个实际使用上面的类的类。它有一个BitmapFile对象作为成员，描述障碍物的高度。如上所述，BitmapFile类将为我们提供一个掩码，即要么是0，要么是1的值（如果你要求的是像素之间的位置，则是在0和1之间插值的值）。这个类将其转化为高度，即低于可变形体表面的0.001（如果BitmapFile类在此位置报告为1）或高于障碍物的0.999（如果BitmapFile类报告为0）。那么下面的函数应该是不言自明的。
 * 

 * 
 * 
 * @code
 *     template <int dim> 
 *     class ChineseObstacle : public Function<dim> 
 *     { 
 *     public: 
 *       ChineseObstacle(const std::string &filename, const double z_surface); 
 * 
 *       virtual double value(const Point<dim> & p, 
 *                            const unsigned int component = 0) const override; 
 * 
 *       virtual void vector_value(const Point<dim> &p, 
 *                                 Vector<double> &  values) const override; 
 * 
 *     private: 
 *       const BitmapFile<dim> input_obstacle; 
 *       double                z_surface; 
 *     }; 
 * 
 *     template <int dim> 
 *     ChineseObstacle<dim>::ChineseObstacle(const std::string &filename, 
 *                                           const double       z_surface) 
 *       : Function<dim>(dim) 
 *       , input_obstacle(filename) 
 *       , z_surface(z_surface) 
 *     {} 
 * 
 *     template <int dim> 
 *     double ChineseObstacle<dim>::value(const Point<dim> & p, 
 *                                        const unsigned int component) const 
 *     { 
 *       if (component == 0) 
 *         return p(0); 
 *       if (component == 1) 
 *         return p(1); 
 *       else if (component == 2) 
 *         { 
 *           if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0) 
 *             return z_surface + 0.999 - input_obstacle.get_value(p(0), p(1)); 
 *         } 
 * 
 *       Assert(false, ExcNotImplemented()); 
 *       return 1e9; // an unreasonable value; ignored in debug mode because of the 
 * 
 * @endcode
 * 
 * 前面的断言
 * 

 * 
 * 
 * @code
 *     } 
 * 
 *     template <int dim> 
 *     void ChineseObstacle<dim>::vector_value(const Point<dim> &p, 
 *                                             Vector<double> &  values) const 
 *     { 
 *       for (unsigned int c = 0; c < this->n_components; ++c) 
 *         values(c) = ChineseObstacle<dim>::value(p, c); 
 *     } 
 *   } // namespace EquationData 
 * @endcode
 * 
 * 
 * <a name="ThecodePlasticityContactProblemcodeclasstemplate"></a> 
 * <h3>The <code>PlasticityContactProblem</code> class template</h3>
 * 

 * 
 * 这是本程序的主类，提供了描述非线性接触问题所需的所有函数和变量。它接近于 step-41 ，但有一些额外的功能，如处理悬挂节点，牛顿方法，使用Trilinos和p4est进行并行分布式计算。处理悬空节点使生活变得有点复杂，因为我们现在需要另一个AffineConstraints对象。我们为接触情况下的主动集合方法创建一个牛顿方法，并处理构成法的非线性算子。
 * 

 * 
 * 这个类的总体布局与其他大多数教程程序非常相似。为了使我们的生活更容易一些，这个类从输入文件中读取一组输入参数。这些参数，使用ParameterHandler类，在 <code>declare_parameters</code> 函数中声明（该函数是静态的，因此它可以在我们创建当前类型的对象之前被调用），然后一个已经用于读取输入文件的ParameterHandler对象将被传递给该类的构造函数。
 * 

 * 
 * 其余的成员函数大体上与我们在其他几个教程程序中看到的一样，虽然为当前的非线性系统增加了一些内容。我们将在下文中对它们的用途进行评论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class PlasticityContactProblem 
 *   { 
 *   public: 
 *     PlasticityContactProblem(const ParameterHandler &prm); 
 * 
 *     void run(); 
 * 
 *     static void declare_parameters(ParameterHandler &prm); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void compute_dirichlet_constraints(); 
 *     void update_solution_and_constraints(); 
 *     void 
 *          assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix); 
 *     void assemble_newton_system( 
 *       const TrilinosWrappers::MPI::Vector &linearization_point); 
 *     void compute_nonlinear_residual( 
 *       const TrilinosWrappers::MPI::Vector &linearization_point); 
 *     void solve_newton_system(); 
 *     void solve_newton(); 
 *     void refine_grid(); 
 *     void move_mesh(const TrilinosWrappers::MPI::Vector &displacement) const; 
 *     void output_results(const unsigned int current_refinement_cycle); 
 * 
 *     void output_contact_force() const; 
 * 
 * @endcode
 * 
 * 就成员变量而言，我们先用一个变量来表示这个程序运行的MPI宇宙，一个我们用来让确切的一个处理器产生输出到控制台的流（见 step-17  ）和一个用来为程序的各个部分计时的变量。
 * 

 * 
 * 
 * @code
 *     MPI_Comm           mpi_communicator; 
 *     ConditionalOStream pcout; 
 *     TimerOutput        computing_timer; 
 * 
 * @endcode
 * 
 * 下一组描述网格和有限元空间。特别是，对于这个并行程序，有限元空间有与之相关的变量，表明哪些自由度存在于当前的处理器上（索引集，也见 step-40 和 @ref distributed 文档模块），以及各种约束：那些由悬挂节点，由Dirichlet边界条件，以及由接触节点的活动集施加的约束。在这里定义的三个AffineConstraints变量中，第一个变量只包含悬挂节点的约束，第二个变量也包含与Dirichlet边界条件相关的约束，第三个变量包含这些约束和接触约束。
 * 

 * 
 * 变量 <code>active_set</code> 包括那些由接触约束的自由度，我们用 <code>fraction_of_plastic_q_points_per_cell</code> 来跟踪每个单元上应力等于屈服应力的正交点的分数。后者仅用于创建显示塑性区的图形输出，但不用于任何进一步的计算；该变量是该类的成员变量，因为该信息是作为计算残差的副产品计算的，但仅在很晚的时候使用。(注意，该向量是一个长度等于<i>local mesh</i>上活动单元数量的向量；它从未被用来在处理器之间交换信息，因此可以是一个普通的deal.II向量)。
 * 

 * 
 * 
 * @code
 *     const unsigned int                        n_initial_global_refinements; 
 *     parallel::distributed::Triangulation<dim> triangulation; 
 * 
 *     const unsigned int fe_degree; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     IndexSet locally_owned_dofs; 
 *     IndexSet locally_relevant_dofs; 
 * 
 *     AffineConstraints<double> constraints_hanging_nodes; 
 *     AffineConstraints<double> constraints_dirichlet_and_hanging_nodes; 
 *     AffineConstraints<double> all_constraints; 
 * 
 *     IndexSet      active_set; 
 *     Vector<float> fraction_of_plastic_q_points_per_cell; 
 * 
 * @endcode
 * 
 * 下一个变量块对应的是解决方案和我们需要形成的线性系统。特别是，这包括牛顿矩阵和右手边；与残差（即牛顿右手边）相对应的向量，但我们没有消除其中的各种约束，该向量用于确定在下一次迭代中需要约束哪些自由度；以及一个与介绍中简要提到的 $B$ 矩阵的对角线相对应的向量，并在随文中讨论。
 * 

 * 
 * 
 * @code
 *     TrilinosWrappers::SparseMatrix newton_matrix; 
 * 
 *     TrilinosWrappers::MPI::Vector solution; 
 *     TrilinosWrappers::MPI::Vector newton_rhs; 
 *     TrilinosWrappers::MPI::Vector newton_rhs_uncondensed; 
 *     TrilinosWrappers::MPI::Vector diag_mass_matrix_vector; 
 * 
 * @endcode
 * 
 * 下一个块包含描述材料响应的变量。
 * 

 * 
 * 
 * @code
 *     const double         e_modulus, nu, gamma, sigma_0; 
 *     ConstitutiveLaw<dim> constitutive_law; 
 * 
 * @endcode
 * 
 * 然后是各种各样的其他变量，用于识别参数文件所选择的要求我们建立的网格，被推入可变形体的障碍物，网格细化策略，是否将解决方案从一个网格转移到下一个网格，以及要执行多少个网格细化循环。在可能的情况下，我们将这些类型的变量标记为 <code>const</code> ，以帮助读者识别哪些变量以后可能会被修改，哪些可能不会被修改（输出目录是一个例外--它在构造函数之外从不被修改，但在构造函数中冒号后面的成员初始化列表中初始化是很尴尬的，因为在那里我们只有一次机会设置它；网格细化准则也是如此）。
 * 

 * 
 * 
 * @code
 *     const std::string                          base_mesh; 
 *     const std::shared_ptr<const Function<dim>> obstacle; 
 * 
 *     struct RefinementStrategy 
 *     { 
 *       enum value 
 *       { 
 *         refine_global, 
 *         refine_percentage, 
 *         refine_fix_dofs 
 *       }; 
 *     }; 
 *     typename RefinementStrategy::value refinement_strategy; 
 * 
 *     const bool         transfer_solution; 
 *     std::string        output_dir; 
 *     const unsigned int n_refinement_cycles; 
 *     unsigned int       current_refinement_cycle; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodePlasticityContactProblemcodeclass"></a> 
 * <h3>Implementation of the <code>PlasticityContactProblem</code> class</h3>
 * 
 * <a name="PlasticityContactProblemdeclare_parameters"></a> 
 * <h4>PlasticityContactProblem::declare_parameters</h4>
 * 

 * 
 * 让我们从声明可在输入文件中选择的运行时参数开始。这些值将在本类的构造函数中读回，以初始化本类的成员变量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::declare_parameters(ParameterHandler &prm) 
 *   { 
 *     prm.declare_entry( 
 *       "polynomial degree", 
 *       "1", 
 *       Patterns::Integer(), 
 *       "Polynomial degree of the FE_Q finite element space, typically 1 or 2."); 
 *     prm.declare_entry("number of initial refinements", 
 *                       "2", 
 *                       Patterns::Integer(), 
 *                       "Number of initial global mesh refinement steps before " 
 *                       "the first computation."); 
 *     prm.declare_entry( 
 *       "refinement strategy", 
 *       "percentage", 
 *       Patterns::Selection("global|percentage"), 
 *       "Mesh refinement strategy:\n" 
 *       " global: one global refinement\n" 
 *       " percentage: a fixed percentage of cells gets refined using the Kelly estimator."); 
 *     prm.declare_entry("number of cycles", 
 *                       "5", 
 *                       Patterns::Integer(), 
 *                       "Number of adaptive mesh refinement cycles to run."); 
 *     prm.declare_entry( 
 *       "obstacle", 
 *       "sphere", 
 *       Patterns::Selection("sphere|read from file"), 
 *       "The name of the obstacle to use. This may either be 'sphere' if we should " 
 *       "use a spherical obstacle, or 'read from file' in which case the obstacle " 
 *       "will be read from a file named 'obstacle.pbm' that is supposed to be in " 
 *       "ASCII PBM format."); 
 *     prm.declare_entry( 
 *       "output directory", 
 *       "", 
 *       Patterns::Anything(), 
 *       "Directory for output files (graphical output and benchmark " 
 *       "statistics). If empty, use the current directory."); 
 *     prm.declare_entry( 
 *       "transfer solution", 
 *       "false", 
 *       Patterns::Bool(), 
 *       "Whether the solution should be used as a starting guess " 
 *       "for the next finer mesh. If false, then the iteration starts at " 
 *       "zero on every mesh."); 
 *     prm.declare_entry("base mesh", 
 *                       "box", 
 *                       Patterns::Selection("box|half sphere"), 
 *                       "Select the shape of the domain: 'box' or 'half sphere'"); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ThecodePlasticityContactProblemcodeconstructor"></a> 
 * <h4>The <code>PlasticityContactProblem</code> constructor</h4>
 * 

 * 
 * 鉴于成员变量的声明以及从输入文件中读取的运行时参数的声明，在这个构造函数中没有任何令人惊讶的地方。在正文中，我们初始化了网格细化策略和输出目录，必要时创建这样一个目录。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   PlasticityContactProblem<dim>::PlasticityContactProblem( 
 *     const ParameterHandler &prm) 
 *     : mpi_communicator(MPI_COMM_WORLD) 
 *     , pcout(std::cout, 
 *             (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)) 
 *     , computing_timer(MPI_COMM_WORLD, 
 *                       pcout, 
 *                       TimerOutput::never, 
 *                       TimerOutput::wall_times) 
 * 
 *     , n_initial_global_refinements( 
 *         prm.get_integer("number of initial refinements")) 
 *     , triangulation(mpi_communicator) 
 *     , fe_degree(prm.get_integer("polynomial degree")) 
 *     , fe(FE_Q<dim>(QGaussLobatto<1>(fe_degree + 1)), dim) 
 *     , dof_handler(triangulation) 
 * 
 *     , e_modulus(200000) 
 *     , nu(0.3) 
 *     , gamma(0.01) 
 *     , sigma_0(400.0) 
 *     , constitutive_law(e_modulus, nu, sigma_0, gamma) 
 * 
 *     , base_mesh(prm.get("base mesh")) 
 *     , obstacle(prm.get("obstacle") == "read from file" ? 
 *                  static_cast<const Function<dim> *>( 
 *                    new EquationData::ChineseObstacle<dim>( 
 *                      "obstacle.pbm", 
 *                      (base_mesh == "box" ? 1.0 : 0.5))) : 
 *                  static_cast<const Function<dim> *>( 
 *                    new EquationData::SphereObstacle<dim>( 
 *                      base_mesh == "box" ? 1.0 : 0.5))) 
 * 
 *     , transfer_solution(prm.get_bool("transfer solution")) 
 *     , n_refinement_cycles(prm.get_integer("number of cycles")) 
 *     , current_refinement_cycle(0) 
 * 
 *   { 
 *     std::string strat = prm.get("refinement strategy"); 
 *     if (strat == "global") 
 *       refinement_strategy = RefinementStrategy::refine_global; 
 *     else if (strat == "percentage") 
 *       refinement_strategy = RefinementStrategy::refine_percentage; 
 *     else 
 *       AssertThrow(false, ExcNotImplemented()); 
 * 
 *     output_dir = prm.get("output directory"); 
 *     if (output_dir != "" && *(output_dir.rbegin()) != '/') 
 *       output_dir += "/"; 
 * 
 * @endcode
 * 
 * 如果有必要，为输出创建一个新的目录。
 * 

 * 
 * 
 * @code
 *     if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) 
 *       { 
 *         const int ierr = mkdir(output_dir.c_str(), 0777); 
 *         AssertThrow(ierr == 0 || errno == EEXIST, ExcIO()); 
 *       } 
 * 
 *     pcout << "    Using output directory '" << output_dir << "'" << std::endl; 
 *     pcout << "    FE degree " << fe_degree << std::endl; 
 *     pcout << "    transfer solution " << (transfer_solution ? "true" : "false") 
 *           << std::endl; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemmake_grid"></a> 
 * <h4>PlasticityContactProblem::make_grid</h4>
 * 

 * 
 * 下一个区块是关于构建起始网格的。我们将使用下面的辅助函数和 <code>make_grid()</code> 的第一个块来构造一个对应于半球形的网格。deal.II有一个函数可以创建这样的网格，但是它的位置和方向都是错误的，所以我们需要在使用它之前对它进行一些位移和旋转。
 * 

 * 
 * 供以后参考，如 GridGenerator::half_hyper_ball(), 文件中所述，半球体的平坦表面的边界指标为零，而其余部分的边界指标为一。
 * 

 * 
 * 
 * @code
 *   Point<3> rotate_half_sphere(const Point<3> &in) 
 *   { 
 *     return {in(2), in(1), -in(0)}; 
 *   } 
 * 
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::make_grid() 
 *   { 
 *     if (base_mesh == "half sphere") 
 *       { 
 *         const Point<dim> center(0, 0, 0); 
 *         const double     radius = 0.8; 
 *         GridGenerator::half_hyper_ball(triangulation, center, radius); 
 * 
 * @endcode
 * 
 * 由于我们将在下面附加一个不同的流形，我们立即清除默认的流形描述。
 * 

 * 
 * 
 * @code
 *         triangulation.reset_all_manifolds(); 
 * 
 *         GridTools::transform(&rotate_half_sphere, triangulation); 
 *         GridTools::shift(Point<dim>(0.5, 0.5, 0.5), triangulation); 
 * 
 *         SphericalManifold<dim> manifold_description(Point<dim>(0.5, 0.5, 0.5)); 
 *         GridTools::copy_boundary_to_manifold_id(triangulation); 
 *         triangulation.set_manifold(0, manifold_description); 
 *       } 
 * 
 * @endcode
 * 
 * 或者，创建一个超立方体网格。创建后，按如下方式分配边界指标。
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *  >     _______
 *  >    /  1    /|
 *  >   /______ / |
 *  >  |       | 8|
 *  >  |   8   | /
 *  >  |_______|/
 *  >      6
 *  @endcode
 * </div>
 * 换句话说，立方体的边的边界指标是8。底部的边界指标是6，顶部的指标是1。我们通过循环所有面的所有单元并查看单元中心的坐标值来设置这些指标，并在以后评估哪个边界将携带迪里希特边界条件或将受到潜在接触时使用这些指标。(在目前的情况下，网格只包含一个单元，它的所有面都在边界上，所以严格来说，所有单元的循环和查询一个面是否在边界上都是不必要的；我们保留它们只是出于习惯：这种代码可以在许多程序中找到，基本上都是这种形式。)
 * 

 * 
 * 
 * @code
 *     else 
 *       { 
 *         const Point<dim> p1(0, 0, 0); 
 *         const Point<dim> p2(1.0, 1.0, 1.0); 
 * 
 *         GridGenerator::hyper_rectangle(triangulation, p1, p2); 
 * 
 *         for (const auto &cell : triangulation.active_cell_iterators()) 
 *           for (const auto &face : cell->face_iterators()) 
 *             if (face->at_boundary()) 
 *               { 
 *                 if (std::fabs(face->center()[2] - p2[2]) < 1e-12) 
 *                   face->set_boundary_id(1); 
 *                 if (std::fabs(face->center()[0] - p1[0]) < 1e-12 || 
 *                     std::fabs(face->center()[0] - p2[0]) < 1e-12 || 
 *                     std::fabs(face->center()[1] - p1[1]) < 1e-12 || 
 *                     std::fabs(face->center()[1] - p2[1]) < 1e-12) 
 *                   face->set_boundary_id(8); 
 *                 if (std::fabs(face->center()[2] - p1[2]) < 1e-12) 
 *                   face->set_boundary_id(6); 
 *               } 
 *       } 
 * 
 *     triangulation.refine_global(n_initial_global_refinements); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsetup_system"></a> 
 * <h4>PlasticityContactProblem::setup_system</h4>
 * 

 * 
 * 谜题的下一块是设置DoFHandler，调整向量大小，并处理其他各种状态变量，如索引集和约束矩阵。
 * 

 * 
 * 在下面的内容中，每一组操作都被放入一个大括号封闭的块中，该块的顶部声明的变量正在进行计时（ TimerOutput::Scope 变量的构造器开始计时部分，在块的末端调用的析构器再次停止计时）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::setup_system() 
 *   { 
 * 
 * /* 设置dofs，并为本地拥有的相关dofs获取索引集  */ 
 * 
 *     { 
 *       TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs"); 
 *       dof_handler.distribute_dofs(fe); 
 * 
 *       locally_owned_dofs = dof_handler.locally_owned_dofs(); 
 *       locally_relevant_dofs.clear(); 
 *       DoFTools::extract_locally_relevant_dofs(dof_handler, 
 *                                               locally_relevant_dofs); 
 *     } 
 * 
 * /*设置悬挂节点和Dirichlet约束 */ 
 * 
 *  
 *     { 
 *       TimerOutput::Scope t(computing_timer, "Setup: constraints"); 
 *       constraints_hanging_nodes.reinit(locally_relevant_dofs); 
 *       DoFTools::make_hanging_node_constraints(dof_handler, 
 *                                               constraints_hanging_nodes); 
 *       constraints_hanging_nodes.close(); 
 * 
 *       pcout << "   Number of active cells: " 
 *             << triangulation.n_global_active_cells() << std::endl 
 *             << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *             << std::endl; 
 * 
 *       compute_dirichlet_constraints(); 
 *     } 
 * 
 * /* 初始化向量和活动集  */ 
 * 
 *     { 
 *       TimerOutput::Scope t(computing_timer, "Setup: vectors"); 
 *       solution.reinit(locally_relevant_dofs, mpi_communicator); 
 *       newton_rhs.reinit(locally_owned_dofs, mpi_communicator); 
 *       newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator); 
 *       diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator); 
 *       fraction_of_plastic_q_points_per_cell.reinit( 
 *         triangulation.n_active_cells()); 
 * 
 *       active_set.clear(); 
 *       active_set.set_size(dof_handler.n_dofs()); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，我们设置了稀疏模式和矩阵。我们暂时（ab）用系统矩阵来同时建立（对角线）矩阵，用于消除与障碍物接触的自由度，但我们随后立即将牛顿矩阵设回零。
 * 

 * 
 * 
 * @code
 *     { 
 *       TimerOutput::Scope                t(computing_timer, "Setup: matrix"); 
 *       TrilinosWrappers::SparsityPattern sp(locally_owned_dofs, 
 *                                            mpi_communicator); 
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler, 
 *                                       sp, 
 *                                       constraints_dirichlet_and_hanging_nodes, 
 *                                       false, 
 *                                       Utilities::MPI::this_mpi_process( 
 *                                         mpi_communicator)); 
 *       sp.compress(); 
 *       newton_matrix.reinit(sp); 
 * 
 *       TrilinosWrappers::SparseMatrix &mass_matrix = newton_matrix; 
 * 
 *       assemble_mass_matrix_diagonal(mass_matrix); 
 * 
 *       const unsigned int start = (newton_rhs.local_range().first), 
 *                          end   = (newton_rhs.local_range().second); 
 *       for (unsigned int j = start; j < end; ++j) 
 *         diag_mass_matrix_vector(j) = mass_matrix.diag_element(j); 
 *       diag_mass_matrix_vector.compress(VectorOperation::insert); 
 * 
 *       mass_matrix = 0; 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemcompute_dirichlet_constraints"></a> 
 * <h4>PlasticityContactProblem::compute_dirichlet_constraints</h4>
 * 

 * 
 * 这个函数从前面的函数中分离出来，计算与迪里切特型边界条件相关的约束，并通过与来自悬挂节点的约束合并，将其放入 <code>constraints_dirichlet_and_hanging_nodes</code> 变量。
 * 

 * 
 * 正如在介绍中所阐述的，我们需要区分两种情况。
 * 

 * 
 * - 如果域是一个盒子，我们将底部的位移设置为零，并允许沿侧面的Z方向的垂直运动。如 <code>make_grid()</code> 函数所示，前者对应于边界指标6，后者对应于8。
 * 

 * 
 * - 如果域是一个半球形，那么我们沿边界的弯曲部分施加零位移，与边界指标0相关。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::compute_dirichlet_constraints() 
 *   { 
 *     constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs); 
 *     constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes); 
 * 
 *     if (base_mesh == "box") 
 *       { 
 * 
 * @endcode
 * 
 * 插值解决方案的所有组成部分
 * 

 * 
 * 
 * @code
 *         VectorTools::interpolate_boundary_values( 
 *           dof_handler, 
 *           6, 
 *           EquationData::BoundaryValues<dim>(), 
 *           constraints_dirichlet_and_hanging_nodes, 
 *           ComponentMask()); 
 * 
 * @endcode
 * 
 * 对解决方案的X和Y分量进行插值（这是一个位掩码，所以应用运算器|）。
 * 

 * 
 * 
 * @code
 *         const FEValuesExtractors::Scalar x_displacement(0); 
 *         const FEValuesExtractors::Scalar y_displacement(1); 
 *         VectorTools::interpolate_boundary_values( 
 *           dof_handler, 
 *           8, 
 *           EquationData::BoundaryValues<dim>(), 
 *           constraints_dirichlet_and_hanging_nodes, 
 *           (fe.component_mask(x_displacement) | 
 *            fe.component_mask(y_displacement))); 
 *       } 
 *     else 
 *       VectorTools::interpolate_boundary_values( 
 *         dof_handler, 
 *         0, 
 *         EquationData::BoundaryValues<dim>(), 
 *         constraints_dirichlet_and_hanging_nodes, 
 *         ComponentMask()); 
 * 
 *     constraints_dirichlet_and_hanging_nodes.close(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemassemble_mass_matrix_diagonal"></a> 
 * <h4>PlasticityContactProblem::assemble_mass_matrix_diagonal</h4>
 * 

 * 
 * 下一个辅助函数计算（对角线）质量矩阵，用于确定我们在接触算法中使用的主动集合方法的主动集合。这个矩阵是质量矩阵类型的，但与标准质量矩阵不同，我们可以通过使用正交公式使其成为对角线（即使在高阶元素的情况下），该公式的正交点与有限元插值点的位置完全相同。我们通过使用QGaussLobatto正交公式来实现这一点，同时用一组从同一正交公式得出的插值点初始化有限元。该函数的其余部分相对简单：我们将得到的矩阵放入给定的参数中；因为我们知道矩阵是对角线的，所以只需在 $i$ 而不是 $j$ 上有一个循环即可。严格来说，我们甚至可以避免在正交点 <code>q_point</code> 处将形状函数的值与自身相乘，因为我们知道形状值是一个恰好有一个的向量，当与自身相点时产生1。由于这个函数不是时间关键，为了清楚起见，我们添加了这个术语。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal( 
 *     TrilinosWrappers::SparseMatrix &mass_matrix) 
 *   { 
 *     QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1); 
 * 
 *     FEFaceValues<dim> fe_values_face(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector displacement(0); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary() && face->boundary_id() == 1) 
 *             { 
 *               fe_values_face.reinit(cell, face); 
 *               cell_matrix = 0; 
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                    ++q_point) 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   cell_matrix(i, i) += 
 *                     (fe_values_face[displacement].value(i, q_point) * 
 *                      fe_values_face[displacement].value(i, q_point) * 
 *                      fe_values_face.JxW(q_point)); 
 * 
 *               cell->get_dof_indices(local_dof_indices); 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 mass_matrix.add(local_dof_indices[i], 
 *                                 local_dof_indices[i], 
 *                                 cell_matrix(i, i)); 
 *             } 
 *     mass_matrix.compress(VectorOperation::add); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemupdate_solution_and_constraints"></a> 
 * <h4>PlasticityContactProblem::update_solution_and_constraints</h4>
 * 

 * 
 * 下面的函数是我们在 <code>solve_newton()</code> 函数中每次牛顿迭代时调用的第一个函数。它的作用是将解决方案投射到可行集上，并更新接触或穿透障碍物的自由度的活动集。
 * 

 * 
 * 为了实现这个功能，我们首先需要做一些记账工作。我们需要写入解决方案向量（我们只能用没有鬼魂元素的完全分布的向量来做），我们需要从各自的向量中读取拉格朗日乘数和对角线质量矩阵的元素（我们只能用有鬼魂元素的向量来做），所以我们创建各自的向量。然后我们还要初始化约束对象，该对象将包含来自接触和所有其他来源的约束，以及一个包含所有属于接触的本地自由度的索引集的对象。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::update_solution_and_constraints() 
 *   { 
 *     std::vector<bool> dof_touched(dof_handler.n_dofs(), false); 
 * 
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
 *                                                        mpi_communicator); 
 *     distributed_solution = solution; 
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
 *                                          mpi_communicator); 
 *     lambda = newton_rhs_uncondensed; 
 * 
 *     TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant( 
 *       locally_relevant_dofs, mpi_communicator); 
 *     diag_mass_matrix_vector_relevant = diag_mass_matrix_vector; 
 * 
 *     all_constraints.reinit(locally_relevant_dofs); 
 *     active_set.clear(); 
 * 
 * @endcode
 * 
 * 第二部分是在所有单元格上的循环，在这个循环中，我们看每一个自由度被定义的点的活动集条件是否为真，我们需要把这个自由度加入到接触节点的活动集中。正如我们一直所做的，如果我们想在单个点上评估函数，我们用一个FEValues对象（或者，这里是FEFaceValues对象，因为我们需要检查表面的接触）和一个适当选择的正交对象来做。我们通过选择定义在单元格面上的形状函数的 "支持点 "来创建这个面的正交对象（关于支持点的更多信息，请参见这个 @ref GlossSupport "词汇表条目"）。因此，我们有多少个正交点，就有多少个面的形状函数，在正交点上循环就相当于在面的形状函数上循环。有了这个，代码看起来如下。
 * 

 * 
 * 
 * @code
 *     Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points()); 
 *     FEFaceValues<dim>   fe_values_face(fe, 
 *                                      face_quadrature, 
 *                                      update_quadrature_points); 
 * 
 *     const unsigned int dofs_per_face   = fe.n_dofs_per_face(); 
 *     const unsigned int n_face_q_points = face_quadrature.size(); 
 * 
 *     std::vector<types::global_dof_index> dof_indices(dofs_per_face); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (!cell->is_artificial()) 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary() && face->boundary_id() == 1) 
 *             { 
 *               fe_values_face.reinit(cell, face); 
 *               face->get_dof_indices(dof_indices); 
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                    ++q_point) 
 *                 { 
 * 
 * @endcode
 * 
 * 在每个正交点（即位于接触边界上的自由度的每个支持点），我们再询问它是否是z-位移自由度的一部分，如果我们还没有遇到这个自由度（对于那些位于面之间的边缘的自由度可能发生），我们需要评估变形物体与障碍物之间的间隙。如果活动集条件为真，那么我们在AffineConstraints对象中添加一个约束，下一次牛顿更新需要满足这个约束，将求解向量的相应元素设置为正确的值，并将索引添加到IndexSet对象中，该索引存储哪个自由度是接触的一部分。
 * 

 * 
 * 
 * @code
 *                   const unsigned int component = 
 *                     fe.face_system_to_component_index(q_point).first; 
 * 
 *                   const unsigned int index_z = dof_indices[q_point]; 
 * 
 *                   if ((component == 2) && (dof_touched[index_z] == false)) 
 *                     { 
 *                       dof_touched[index_z] = true; 
 * 
 *                       const Point<dim> this_support_point = 
 *                         fe_values_face.quadrature_point(q_point); 
 * 
 *                       const double obstacle_value = 
 *                         obstacle->value(this_support_point, 2); 
 *                       const double solution_here = solution(index_z); 
 *                       const double undeformed_gap = 
 *                         obstacle_value - this_support_point(2); 
 * 
 *                       const double c = 100.0 * e_modulus; 
 *                       if ((lambda(index_z) / 
 *                                diag_mass_matrix_vector_relevant(index_z) + 
 *                              c * (solution_here - undeformed_gap) > 
 *                            0) && 
 *                           !constraints_hanging_nodes.is_constrained(index_z)) 
 *                         { 
 *                           all_constraints.add_line(index_z); 
 *                           all_constraints.set_inhomogeneity(index_z, 
 *                                                             undeformed_gap); 
 *                           distributed_solution(index_z) = undeformed_gap; 
 * 
 *                           active_set.add_index(index_z); 
 *                         } 
 *                     } 
 *                 } 
 *             } 
 * 
 * @endcode
 * 
 * 在这个函数的最后，我们在处理器之间交换数据，更新 <code>solution</code> 变量中那些已经被其他处理器写入的幽灵元素。然后我们将Dirichlet约束和那些来自悬挂节点的约束合并到已经包含活动集的AffineConstraints对象中。我们通过输出主动约束自由度的总数来结束这个函数，对于这个自由度，我们对每个处理器拥有的主动约束自由度的数量进行加总。这个本地拥有的受限自由度的数量当然是活动集和本地拥有的自由度集的交集的元素数量，我们可以通过在两个IndexSets上使用 <code>operator&</code> 得到。
 * 

 * 
 * 
 * @code
 *     distributed_solution.compress(VectorOperation::insert); 
 *     solution = distributed_solution; 
 * 
 *     all_constraints.close(); 
 *     all_constraints.merge(constraints_dirichlet_and_hanging_nodes); 
 * 
 *  
 *           << Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(), 
 *                                  mpi_communicator) 
 *           << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemassemble_newton_system"></a> 
 * <h4>PlasticityContactProblem::assemble_newton_system</h4>
 * 

 * 
 * 鉴于问题的复杂性，可能会让人感到惊讶的是，在每次牛顿迭代中组装我们要解决的线性系统实际上是相当简单的。下面的函数建立了牛顿的右手边和牛顿矩阵。它看起来相当简单，因为繁重的工作发生在对 <code>ConstitutiveLaw::get_linearized_stress_strain_tensors()</code> 的调用中，特别是在 AffineConstraints::distribute_local_to_global(), 中使用我们之前计算的约束。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::assemble_newton_system( 
 *     const TrilinosWrappers::MPI::Vector &linearization_point) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Assembling"); 
 * 
 *     QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *     QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_values_face(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_quadrature_points | 
 *                                        update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points      = quadrature_formula.size(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     const EquationData::BoundaryForce<dim> boundary_force; 
 *     std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
 *                                                       Vector<double>(dim)); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector displacement(0); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           fe_values.reinit(cell); 
 *           cell_matrix = 0; 
 *           cell_rhs    = 0; 
 * 
 *           std::vector<SymmetricTensor<2, dim>> strain_tensor(n_q_points); 
 *           fe_values[displacement].get_function_symmetric_gradients( 
 *             linearization_point, strain_tensor); 
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             { 
 *               SymmetricTensor<4, dim> stress_strain_tensor_linearized; 
 *               SymmetricTensor<4, dim> stress_strain_tensor; 
 *               constitutive_law.get_linearized_stress_strain_tensors( 
 *                 strain_tensor[q_point], 
 *                 stress_strain_tensor_linearized, 
 *                 stress_strain_tensor); 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 * 
 * @endcode
 * 
 * 在计算了应力-应变张量及其线性化之后，我们现在可以把矩阵和右手边的部分放在一起。在这两部分中，我们需要线性化的应力-应变张量乘以 $\varphi_i$ 的对称梯度，即 $I_\Pi\varepsilon(\varphi_i)$ 项，因此我们引入这个项的缩写。回顾一下，该矩阵对应于随附出版物的符号中的双线性形式 $A_{ij}=(I_\Pi\varepsilon(\varphi_i),\varepsilon(\varphi_j))$ ，而右手边是 $F_i=([I_\Pi-P_\Pi C]\varepsilon(\varphi_i),\varepsilon(\mathbf u))$ ，其中 $u$ 是当前的线性化点（通常是最后的解）。这可能表明，如果材料是完全弹性的（其中 $I_\Pi=P_\Pi$ ），右手边将为零，但这忽略了一个事实，即右手边还将包含由于接触而产生的非均质约束的贡献。                
 * 接下来的代码块增加了由于边界力的贡献，如果有的话。
 * 

 * 
 * 
 * @code
 *                   const SymmetricTensor<2, dim> stress_phi_i = 
 *                     stress_strain_tensor_linearized * 
 *                     fe_values[displacement].symmetric_gradient(i, q_point); 
 * 
 *                   for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                     cell_matrix(i, j) += 
 *                       (stress_phi_i * 
 *                        fe_values[displacement].symmetric_gradient(j, q_point) * 
 *                        fe_values.JxW(q_point)); 
 * 
 *                   cell_rhs(i) += 
 *                     ((stress_phi_i - 
 *                       stress_strain_tensor * 
 *                         fe_values[displacement].symmetric_gradient(i, 
 *                                                                    q_point)) * 
 *                      strain_tensor[q_point] * fe_values.JxW(q_point)); 
 *                 } 
 *             } 
 * 
 *           for (const auto &face : cell->face_iterators()) 
 *             if (face->at_boundary() && face->boundary_id() == 1) 
 *               { 
 *                 fe_values_face.reinit(cell, face); 
 * 
 *                 boundary_force.vector_value_list( 
 *                   fe_values_face.get_quadrature_points(), 
 *                   boundary_force_values); 
 * 
 *                 for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                      ++q_point) 
 *                   { 
 *                     Tensor<1, dim> rhs_values; 
 *                     rhs_values[2] = boundary_force_values[q_point][2]; 
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                       cell_rhs(i) += 
 *                         (fe_values_face[displacement].value(i, q_point) * 
 *                          rhs_values * fe_values_face.JxW(q_point)); 
 *                   } 
 *               } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           all_constraints.distribute_local_to_global(cell_matrix, 
 *                                                      cell_rhs, 
 *                                                      local_dof_indices, 
 *                                                      newton_matrix, 
 *                                                      newton_rhs, 
 *                                                      true); 
 *         } 
 * 
 *     newton_matrix.compress(VectorOperation::add); 
 *     newton_rhs.compress(VectorOperation::add); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemcompute_nonlinear_residual"></a> 
 * <h4>PlasticityContactProblem::compute_nonlinear_residual</h4>
 * 

 * 
 * 下面的函数计算给定当前解（或任何其他线性化点）的方程的非线性残差。这在线性搜索算法中是需要的，我们需要尝试之前和当前（试验）解的各种线性组合来计算当前牛顿步骤的（真实的、全局化的）解。
 * 

 * 
 * 说到这里，在稍微滥用函数名称的情况下，它实际上做了很多事情。例如，它还计算出与牛顿残差相对应的矢量，但没有消除受限自由度。我们需要这个向量来计算接触力，并最终计算出下一个活动集。同样，通过跟踪我们在每个单元上遇到的显示塑性屈服的正交点的数量，我们也可以计算出 <code>fraction_of_plastic_q_points_per_cell</code> 矢量，随后我们可以输出这个矢量来可视化塑性区。在这两种情况下，作为线条搜索的一部分，这些结果是不必要的，因此我们可能会浪费少量的时间来计算它们。同时，无论如何，这些信息是我们在这里需要做的事情的自然副产品，而且我们想在每个牛顿步骤结束时收集一次，所以我们不妨在这里做。
 * 

 * 
 * 这个函数的实际实现应该是相当明显的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::compute_nonlinear_residual( 
 *     const TrilinosWrappers::MPI::Vector &linearization_point) 
 *   { 
 *     QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *     QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_JxW_values); 
 * 
 *     FEFaceValues<dim> fe_values_face(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_quadrature_points | 
 *                                        update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell   = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points      = quadrature_formula.size(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     const EquationData::BoundaryForce<dim> boundary_force; 
 *     std::vector<Vector<double>> boundary_force_values(n_face_q_points, 
 *                                                       Vector<double>(dim)); 
 * 
 *     Vector<double> cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector displacement(0); 
 * 
 *     newton_rhs             = 0; 
 *     newton_rhs_uncondensed = 0; 
 * 
 *     fraction_of_plastic_q_points_per_cell = 0; 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         { 
 *           fe_values.reinit(cell); 
 *           cell_rhs = 0; 
 * 
 *           std::vector<SymmetricTensor<2, dim>> strain_tensors(n_q_points); 
 *           fe_values[displacement].get_function_symmetric_gradients( 
 *             linearization_point, strain_tensors); 
 * 
 *           for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *             { 
 *               SymmetricTensor<4, dim> stress_strain_tensor; 
 *               const bool              q_point_is_plastic = 
 *                 constitutive_law.get_stress_strain_tensor( 
 *                   strain_tensors[q_point], stress_strain_tensor); 
 *               if (q_point_is_plastic) 
 *                 ++fraction_of_plastic_q_points_per_cell( 
 *                   cell->active_cell_index()); 
 * 
 *               for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                 { 
 *                   cell_rhs(i) -= 
 *                     (strain_tensors[q_point] * stress_strain_tensor * 
 *                      fe_values[displacement].symmetric_gradient(i, q_point) * 
 *                      fe_values.JxW(q_point)); 
 * 
 *                   Tensor<1, dim> rhs_values; 
 *                   rhs_values = 0; 
 *                   cell_rhs(i) += (fe_values[displacement].value(i, q_point) * 
 *                                   rhs_values * fe_values.JxW(q_point)); 
 *                 } 
 *             } 
 * 
 *           for (const auto &face : cell->face_iterators()) 
 *             if (face->at_boundary() && face->boundary_id() == 1) 
 *               { 
 *                 fe_values_face.reinit(cell, face); 
 * 
 *                 boundary_force.vector_value_list( 
 *                   fe_values_face.get_quadrature_points(), 
 *                   boundary_force_values); 
 * 
 *                 for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                      ++q_point) 
 *                   { 
 *                     Tensor<1, dim> rhs_values; 
 *                     rhs_values[2] = boundary_force_values[q_point][2]; 
 *                     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                       cell_rhs(i) += 
 *                         (fe_values_face[displacement].value(i, q_point) * 
 *                          rhs_values * fe_values_face.JxW(q_point)); 
 *                   } 
 *               } 
 * 
 *           cell->get_dof_indices(local_dof_indices); 
 *           constraints_dirichlet_and_hanging_nodes.distribute_local_to_global( 
 *             cell_rhs, local_dof_indices, newton_rhs); 
 * 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             newton_rhs_uncondensed(local_dof_indices[i]) += cell_rhs(i); 
 *         } 
 * 
 *     fraction_of_plastic_q_points_per_cell /= quadrature_formula.size(); 
 *     newton_rhs.compress(VectorOperation::add); 
 *     newton_rhs_uncondensed.compress(VectorOperation::add); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsolve_newton_system"></a> 
 * <h4>PlasticityContactProblem::solve_newton_system</h4>
 * 

 * 
 * 在我们讨论单个网格上的实际牛顿迭代之前的最后一块是线性系统的求解器。有几个复杂的问题使代码略显模糊，但大多数情况下，它只是设置然后求解。在这些复杂的问题中，包括。
 * 

 * 
 * 

 * 
 * 

 * 
 * 对于悬空节点，我们必须将 AffineConstraints::set_zero 函数应用于newton_rhs。  如果一个求解值为 $x_0$ 的悬空节点有一个与障碍物接触的数值为 $x_1$ 的邻居和一个没有接触的邻居 $x_2$ ，这就有必要。因为前者的更新将是规定的，所以悬挂的节点约束将有一个不均匀性，看起来像  $x_0 = x_1/2 +   \text{gap}/2$  。所以右侧的相应条目是无意义的非零值。这些值我们必须设置为零。
 * 

 * 
 * - 就像在  step-40  中一样，在求解或使用解决方案时，我们需要在有和没有鬼魂元素的向量之间进行洗牌。
 * 

 * 
 * 该函数的其余部分与 step-40 和 step-41 类似，只是我们使用BiCGStab求解器而不是CG。这是由于对于非常小的硬化参数 $\gamma$ ，线性系统变得几乎是半无限的，尽管仍然是对称的。BiCGStab似乎更容易处理这种线性系统。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::solve_newton_system() 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Solve"); 
 * 
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
 *                                                        mpi_communicator); 
 *     distributed_solution = solution; 
 * 
 *     constraints_hanging_nodes.set_zero(distributed_solution); 
 *     constraints_hanging_nodes.set_zero(newton_rhs); 
 * 
 *  
 *     { 
 *       TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner"); 
 * 
 *       std::vector<std::vector<bool>> constant_modes; 
 *       DoFTools::extract_constant_modes(dof_handler, 
 *                                        ComponentMask(), 
 *                                        constant_modes); 
 * 
 *       TrilinosWrappers::PreconditionAMG::AdditionalData additional_data; 
 *       additional_data.constant_modes        = constant_modes; 
 *       additional_data.elliptic              = true; 
 *       additional_data.n_cycles              = 1; 
 *       additional_data.w_cycle               = false; 
 *       additional_data.output_details        = false; 
 *       additional_data.smoother_sweeps       = 2; 
 *       additional_data.aggregation_threshold = 1e-2; 
 * 
 *       preconditioner.initialize(newton_matrix, additional_data); 
 *     } 
 * 
 *     { 
 *       TimerOutput::Scope t(computing_timer, "Solve: iterate"); 
 * 
 *       TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator); 
 * 
 *       const double relative_accuracy = 1e-8; 
 *       const double solver_tolerance = 
 *         relative_accuracy * 
 *         newton_matrix.residual(tmp, distributed_solution, newton_rhs); 
 * 
 *       SolverControl solver_control(newton_matrix.m(), solver_tolerance); 
 *       SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control); 
 *       solver.solve(newton_matrix, 
 *                    distributed_solution, 
 *                    newton_rhs, 
 *                    preconditioner); 
 * 
 *       pcout << "         Error: " << solver_control.initial_value() << " -> " 
 *             << solver_control.last_value() << " in " 
 *             << solver_control.last_step() << " Bicgstab iterations." 
 *             << std::endl; 
 *     } 
 * 
 *     all_constraints.distribute(distributed_solution); 
 * 
 *     solution = distributed_solution; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemsolve_newton"></a> 
 * <h4>PlasticityContactProblem::solve_newton</h4>
 * 

 * 
 * 最后，这是在当前网格上实现阻尼牛顿方法的函数。这里有两个嵌套的循环：外循环用于牛顿迭代，内循环用于直线搜索，只有在必要时才会使用。为了获得一个好的和合理的起始值，我们在每个网格上的第一个牛顿步骤中解决一个弹性问题（如果我们在网格之间转移解决方案，则只在第一个网格上解决）。我们通过在这些迭代中将屈服应力设置为一个不合理的大值，然后在随后的迭代中将其设置为正确值。
 * 

 * 
 * 除此以外，这个函数的顶部部分应该是相当明显的。我们将变量 <code>previous_residual_norm</code> 初始化为可以用双精度数字表示的最大负值，以便在第一步中比较当前残差是否小于前一步的残差时总是失败。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::solve_newton() 
 *   { 
 *     TrilinosWrappers::MPI::Vector old_solution(locally_owned_dofs, 
 *                                                mpi_communicator); 
 *     TrilinosWrappers::MPI::Vector residual(locally_owned_dofs, 
 *                                            mpi_communicator); 
 *     TrilinosWrappers::MPI::Vector tmp_vector(locally_owned_dofs, 
 *                                              mpi_communicator); 
 *     TrilinosWrappers::MPI::Vector locally_relevant_tmp_vector( 
 *       locally_relevant_dofs, mpi_communicator); 
 *     TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
 *                                                        mpi_communicator); 
 * 
 *     double residual_norm; 
 *     double previous_residual_norm = -std::numeric_limits<double>::max(); 
 * 
 *     const double correct_sigma = sigma_0; 
 * 
 *     IndexSet old_active_set(active_set); 
 * 
 *     for (unsigned int newton_step = 1; newton_step <= 100; ++newton_step) 
 *       { 
 *         if (newton_step == 1 && 
 *             ((transfer_solution && current_refinement_cycle == 0) || 
 *              !transfer_solution)) 
 *           constitutive_law.set_sigma_0(1e+10); 
 *         else if (newton_step == 2 || current_refinement_cycle > 0 || 
 *                  !transfer_solution) 
 *           constitutive_law.set_sigma_0(correct_sigma); 
 * 
 *         pcout << " " << std::endl; 
 *         pcout << "   Newton iteration " << newton_step << std::endl; 
 *         pcout << "      Updating active set..." << std::endl; 
 * 
 *         { 
 *           TimerOutput::Scope t(computing_timer, "update active set"); 
 *           update_solution_and_constraints(); 
 *         } 
 * 
 *         pcout << "      Assembling system... " << std::endl; 
 *         newton_matrix = 0; 
 *         newton_rhs    = 0; 
 *         assemble_newton_system(solution); 
 * 
 *         pcout << "      Solving system... " << std::endl; 
 *         solve_newton_system(); 
 * 
 * @endcode
 * 
 * 在我们计算了当前牛顿步骤的试解 $\tilde{\mathbf u}$ 之后，情况就变得有点棘手了。我们处理的是一个高度非线性的问题，所以我们必须用直线搜索的方式来抑制牛顿方法。为了理解我们如何做到这一点，请回顾一下，在我们的表述中，我们在每一个牛顿步骤中计算一个试解，而不是在新旧解之间进行更新。由于解集是一个凸集，我们将使用直线搜索，尝试以前的解和试验解的线性组合，以保证阻尼解再次出现在我们的解集中。我们最多应用5个阻尼步骤。
 * 

 * 
 * 在我们使用直线搜索的时候有一些例外情况。首先，如果这是任何网格上的第一个牛顿步骤，那么我们就没有任何点来比较残差，所以我们总是接受一个完整的步骤。同样地，如果这是第一个网格上的第二个牛顿步骤（如果我们不在网格之间转移解决方案，则是任何网格上的第二个牛顿步骤），则我们只用弹性模型计算了其中的第一个步骤（见上文我们如何将屈服应力σ设置为一个不合理的大值）。在这种情况下，第一个牛顿解是一个纯粹的弹性解，第二个牛顿解是一个塑性解，任何线性组合都不一定会位于可行的集合中--所以我们只是接受我们刚刚得到的解。
 * 

 * 
 * 在这两种情况下，我们绕过直线搜索，只是在必要时更新残差和其他向量。
 * 

 * 
 * 
 * @code
 *         if ((newton_step == 1) || 
 *             (transfer_solution && newton_step == 2 && 
 *              current_refinement_cycle == 0) || 
 *             (!transfer_solution && newton_step == 2)) 
 *           { 
 *             compute_nonlinear_residual(solution); 
 *             old_solution = solution; 
 * 
 *             residual                     = newton_rhs; 
 *             const unsigned int start_res = (residual.local_range().first), 
 *                                end_res   = (residual.local_range().second); 
 *             for (unsigned int n = start_res; n < end_res; ++n) 
 *               if (all_constraints.is_inhomogeneously_constrained(n)) 
 *                 residual(n) = 0; 
 * 
 *             residual.compress(VectorOperation::insert); 
 * 
 *             residual_norm = residual.l2_norm(); 
 * 
 *             pcout << "      Accepting Newton solution with residual: " 
 *                   << residual_norm << std::endl; 
 *           } 
 *         else 
 *           { 
 *             for (unsigned int i = 0; i < 5; ++i) 
 *               { 
 *                 distributed_solution = solution; 
 * 
 *                 const double alpha = std::pow(0.5, static_cast<double>(i)); 
 *                 tmp_vector         = old_solution; 
 *                 tmp_vector.sadd(1 - alpha, alpha, distributed_solution); 
 * 
 *                 TimerOutput::Scope t(computing_timer, "Residual and lambda"); 
 * 
 *                 locally_relevant_tmp_vector = tmp_vector; 
 *                 compute_nonlinear_residual(locally_relevant_tmp_vector); 
 *                 residual = newton_rhs; 
 * 
 *                 const unsigned int start_res = (residual.local_range().first), 
 *                                    end_res   = (residual.local_range().second); 
 *                 for (unsigned int n = start_res; n < end_res; ++n) 
 *                   if (all_constraints.is_inhomogeneously_constrained(n)) 
 *                     residual(n) = 0; 
 * 
 *                 residual.compress(VectorOperation::insert); 
 * 
 *                 residual_norm = residual.l2_norm(); 
 * 
 *  
 *                   << "      Residual of the non-contact part of the system: " 
 *                   << residual_norm << std::endl 
 *                   << "         with a damping parameter alpha = " << alpha 
 *                   << std::endl; 
 * 
 *                 if (residual_norm < previous_residual_norm) 
 *                   break; 
 *               } 
 * 
 *             solution     = tmp_vector; 
 *             old_solution = solution; 
 *           } 
 * 
 *         previous_residual_norm = residual_norm; 
 * 
 * @endcode
 * 
 * 最后一步是检查收敛情况。如果活动集在所有处理器中都没有变化，并且残差小于阈值 $10^{-10}$  ，那么我们就终止对当前网格的迭代。
 * 

 * 
 * 
 * @code
 *         if (Utilities::MPI::sum((active_set == old_active_set) ? 0 : 1, 
 *                                 mpi_communicator) == 0) 
 *           { 
 *             pcout << "      Active set did not change!" << std::endl; 
 *             if (residual_norm < 1e-10) 
 *               break; 
 *           } 
 * 
 *         old_active_set = active_set; 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemrefine_grid"></a> 
 * <h4>PlasticityContactProblem::refine_grid</h4>
 * 

 * 
 * 如果你已经在deal.II教程中做到了这一点，下面这个细化网格的函数应该不会再对你构成任何挑战。它对网格进行细化，可以是全局的，也可以是使用Kelly误差估计器的，如果这样要求的话，还可以将上一个网格的解转移到下一个网格。在后一种情况下，我们还需要再次计算活动集和其他数量，为此我们需要由  <code>compute_nonlinear_residual()</code>  计算的信息。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::refine_grid() 
 *   { 
 *     if (refinement_strategy == RefinementStrategy::refine_global) 
 *       { 
 *         for (typename Triangulation<dim>::active_cell_iterator cell = 
 *                triangulation.begin_active(); 
 *              cell != triangulation.end(); 
 *              ++cell) 
 *           if (cell->is_locally_owned()) 
 *             cell->set_refine_flag(); 
 *       } 
 *     else 
 *       { 
 *         Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 *         KellyErrorEstimator<dim>::estimate( 
 *           dof_handler, 
 *           QGauss<dim - 1>(fe.degree + 2), 
 *           std::map<types::boundary_id, const Function<dim> *>(), 
 *           solution, 
 *           estimated_error_per_cell); 
 * 
 *         parallel::distributed::GridRefinement ::refine_and_coarsen_fixed_number( 
 *           triangulation, estimated_error_per_cell, 0.3, 0.03); 
 *       } 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 * 
 *     parallel::distributed::SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> 
 *       solution_transfer(dof_handler); 
 *     if (transfer_solution) 
 *       solution_transfer.prepare_for_coarsening_and_refinement(solution); 
 * 
 *     triangulation.execute_coarsening_and_refinement(); 
 * 
 *     setup_system(); 
 * 
 *     if (transfer_solution) 
 *       { 
 *         TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, 
 *                                                            mpi_communicator); 
 *         solution_transfer.interpolate(distributed_solution); 
 * 
 * @endcode
 * 
 * 强制执行约束条件，使插值后的解决方案在新的网格上符合要求。
 * 

 * 
 * 
 * @code
 *         constraints_hanging_nodes.distribute(distributed_solution); 
 * 
 *         solution = distributed_solution; 
 *         compute_nonlinear_residual(solution); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemmove_mesh"></a> 
 * <h4>PlasticityContactProblem::move_mesh</h4>
 * 

 * 
 * 在我们到达 <code>run()</code> 之前的其余三个函数都与生成输出有关。下面一个是尝试显示变形体的变形构造。为此，这个函数接收一个位移矢量场，通过先前计算的位移来移动网格（局部）的每个顶点。在生成图形输出之前，我们将以当前的位移场调用该函数，在生成图形输出之后，我们将以负的位移场再次调用该函数，以撤销对网格所做的修改。
 * 

 * 
 * 这个函数本身是非常简单的。我们所要做的就是跟踪我们已经接触过的顶点，因为我们在单元格上循环时多次遇到相同的顶点。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::move_mesh( 
 *     const TrilinosWrappers::MPI::Vector &displacement) const 
 *   { 
 *     std::vector<bool> vertex_touched(triangulation.n_vertices(), false); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         for (const auto v : cell->vertex_indices()) 
 *           if (vertex_touched[cell->vertex_index(v)] == false) 
 *             { 
 *               vertex_touched[cell->vertex_index(v)] = true; 
 * 
 *               Point<dim> vertex_displacement; 
 *               for (unsigned int d = 0; d < dim; ++d) 
 *                 vertex_displacement[d] = 
 *                   displacement(cell->vertex_dof_index(v, d)); 
 * 
 *               cell->vertex(v) += vertex_displacement; 
 *             } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemoutput_results"></a> 
 * <h4>PlasticityContactProblem::output_results</h4>
 * 

 * 
 * 接下来是我们用来实际生成图形输出的函数。这个函数有点繁琐，但实际上并不特别复杂。它在顶部移动网格（最后再把它移回来），然后计算沿接触面的接触力。我们可以通过取未处理的残差向量，并通过询问它们是否有与之相关的不均匀约束来确定哪些自由度对应于有接触的自由度（如随文所示）。一如既往，我们需要注意的是，我们只能写进完全分布的向量（即没有鬼魂元素的向量），但当我们想产生输出时，我们需要的向量确实对所有局部相关的自由度都有鬼魂项。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::output_results( 
 *     const unsigned int current_refinement_cycle) 
 *   { 
 *     TimerOutput::Scope t(computing_timer, "Graphical output"); 
 * 
 *     pcout << "      Writing graphical output... " << std::flush; 
 * 
 *     move_mesh(solution); 
 * 
 * @endcode
 * 
 * 接触力的计算
 * 

 * 
 * 
 * @code
 *     TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
 *                                                      mpi_communicator); 
 *     const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
 *                        end_res = (newton_rhs_uncondensed.local_range().second); 
 *     for (unsigned int n = start_res; n < end_res; ++n) 
 *       if (all_constraints.is_inhomogeneously_constrained(n)) 
 *         distributed_lambda(n) = 
 *           newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
 *     distributed_lambda.compress(VectorOperation::insert); 
 *     constraints_hanging_nodes.distribute(distributed_lambda); 
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
 *                                          mpi_communicator); 
 *     lambda = distributed_lambda; 
 * 
 *     TrilinosWrappers::MPI::Vector distributed_active_set_vector( 
 *       locally_owned_dofs, mpi_communicator); 
 *     distributed_active_set_vector = 0.; 
 *     for (const auto index : active_set) 
 *       distributed_active_set_vector[index] = 1.; 
 *     distributed_lambda.compress(VectorOperation::insert); 
 * 
 *     TrilinosWrappers::MPI::Vector active_set_vector(locally_relevant_dofs, 
 *                                                     mpi_communicator); 
 *     active_set_vector = distributed_active_set_vector; 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 *     const std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_out.add_data_vector(solution, 
 *                              std::vector<std::string>(dim, "displacement"), 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.add_data_vector(lambda, 
 *                              std::vector<std::string>(dim, "contact_force"), 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.add_data_vector(active_set_vector, 
 *                              std::vector<std::string>(dim, "active_set"), 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 * 
 *     Vector<float> subdomain(triangulation.n_active_cells()); 
 *     for (unsigned int i = 0; i < subdomain.size(); ++i) 
 *       subdomain(i) = triangulation.locally_owned_subdomain(); 
 *     data_out.add_data_vector(subdomain, "subdomain"); 
 * 
 *     data_out.add_data_vector(fraction_of_plastic_q_points_per_cell, 
 *                              "fraction_of_plastic_q_points"); 
 * 
 *     data_out.build_patches(); 
 * 
 * @endcode
 * 
 * 在函数的其余部分，我们在每个处理器上生成一个VTU文件，以这个处理器的子域ID为索引。在第一个处理器上，我们随后还创建了一个 <code>.pvtu</code> 文件，对VTU文件的<i>all</i>进行索引，这样就可以一次性读取整个输出文件集。这些 <code>.pvtu</code> 被Paraview用来描述整个并行计算的输出文件。然后我们再为Paraview的竞争者--VisIt可视化程序做同样的事情，创建一个匹配的 <code>.visit</code> 文件。
 * 

 * 
 * 
 * @code
 *     const std::string pvtu_filename = data_out.write_vtu_with_pvtu_record( 
 *       output_dir, "solution", current_refinement_cycle, mpi_communicator, 2); 
 *     pcout << pvtu_filename << std::endl; 
 * 
 *     TrilinosWrappers::MPI::Vector tmp(solution); 
 *     tmp *= -1; 
 *     move_mesh(tmp); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemoutput_contact_force"></a> 
 * <h4>PlasticityContactProblem::output_contact_force</h4>
 * 

 * 
 * 这最后一个辅助函数通过计算接触面积上Z方向的接触压力的积分来计算接触力。为此，我们将所有非活动因子的接触压力lambda设置为0（一个自由度是否是接触的一部分，就像我们在前一个函数中做的那样）。对于所有活动的自由度，lambda包含非线性残差（newton_rhs_uncondensed）和质量矩阵（diag_mass_matrix_vector）的相应对角线条目的商数。因为悬空节点出现在接触区的可能性不小，所以对分布式_lambda向量应用constraints_hanging_nodes.distribution是很重要的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::output_contact_force() const 
 *   { 
 *     TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, 
 *                                                      mpi_communicator); 
 *     const unsigned int start_res = (newton_rhs_uncondensed.local_range().first), 
 *                        end_res = (newton_rhs_uncondensed.local_range().second); 
 *     for (unsigned int n = start_res; n < end_res; ++n) 
 *       if (all_constraints.is_inhomogeneously_constrained(n)) 
 *         distributed_lambda(n) = 
 *           newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n); 
 *       else 
 *         distributed_lambda(n) = 0; 
 *     distributed_lambda.compress(VectorOperation::insert); 
 *     constraints_hanging_nodes.distribute(distributed_lambda); 
 * 
 *     TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, 
 *                                          mpi_communicator); 
 *     lambda = distributed_lambda; 
 * 
 *     double contact_force = 0.0; 
 * 
 *     QGauss<dim - 1>   face_quadrature_formula(fe.degree + 1); 
 *     FEFaceValues<dim> fe_values_face(fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_JxW_values); 
 * 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     const FEValuesExtractors::Vector displacement(0); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       if (cell->is_locally_owned()) 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary() && face->boundary_id() == 1) 
 *             { 
 *               fe_values_face.reinit(cell, face); 
 * 
 *               std::vector<Tensor<1, dim>> lambda_values(n_face_q_points); 
 *               fe_values_face[displacement].get_function_values(lambda, 
 *                                                                lambda_values); 
 * 
 *               for (unsigned int q_point = 0; q_point < n_face_q_points; 
 *                    ++q_point) 
 *                 contact_force += 
 *                   lambda_values[q_point][2] * fe_values_face.JxW(q_point); 
 *             } 
 *     contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD); 
 * 
 *     pcout << "Contact force = " << contact_force << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="PlasticityContactProblemrun"></a> 
 * <h4>PlasticityContactProblem::run</h4>
 * 

 * 
 * 和其他所有的教程程序一样， <code>run()</code> 函数包含了整体逻辑。这里没有太多的内容：本质上，它在所有的网格细化循环中执行循环，并在每个循环中，将事情交给 <code>solve_newton()</code> 中的牛顿求解器，并调用函数来创建如此计算的解决方案的图形输出。然后输出一些关于运行时间和内存消耗的统计数据，这些数据是在这个网格的计算过程中收集的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void PlasticityContactProblem<dim>::run() 
 *   { 
 *     computing_timer.reset(); 
 *     for (; current_refinement_cycle < n_refinement_cycles; 
 *          ++current_refinement_cycle) 
 *       { 
 *         { 
 *           TimerOutput::Scope t(computing_timer, "Setup"); 
 * 
 *           pcout << std::endl; 
 *           pcout << "Cycle " << current_refinement_cycle << ':' << std::endl; 
 * 
 *           if (current_refinement_cycle == 0) 
 *             { 
 *               make_grid(); 
 *               setup_system(); 
 *             } 
 *           else 
 *             { 
 *               TimerOutput::Scope t(computing_timer, "Setup: refine mesh"); 
 *               refine_grid(); 
 *             } 
 *         } 
 * 
 *         solve_newton(); 
 * 
 *         output_results(current_refinement_cycle); 
 * 
 *         computing_timer.print_summary(); 
 *         computing_timer.reset(); 
 * 
 *         Utilities::System::MemoryStats stats; 
 *         Utilities::System::get_memory_stats(stats); 
 *         pcout << "Peak virtual memory used, resident in kB: " << stats.VmSize 
 *               << " " << stats.VmRSS << std::endl; 
 * 
 *         if (base_mesh == "box") 
 *           output_contact_force(); 
 *       } 
 *   } 
 * } // namespace Step42 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * <code>main()</code> 函数真的没有什么内容。看起来他们总是这样做。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   using namespace dealii; 
 *   using namespace Step42; 
 * 
 *   try 
 *     { 
 *       ParameterHandler prm; 
 *       PlasticityContactProblem<3>::declare_parameters(prm); 
 *       if (argc != 2) 
 *         { 
 *           std::cerr << "*** Call this program as <./step-42 input.prm>" 
 *                     << std::endl; 
 *           return 1; 
 *         } 
 * 
 *       prm.parse_input(argv[1]); 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization( 
 *         argc, argv, numbers::invalid_unsigned_int); 
 *       { 
 *         PlasticityContactProblem<3> problem(prm); 
 *         problem.run(); 
 *       } 
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
 * @endcode
examples/step-42/doc/results.dox



<a name="Results"></a><h1>Results</h1>


包含这个程序的目录还包含一些输入参数文件，可以用来创建各种不同的模拟。例如，用 <code>p1_adaptive.prm</code> 参数文件（用球作为障碍物，用盒子作为领域）在16个核心上运行该程序会产生这样的输出。

@code
    Using output directory 'p1adaptive/'
    FE degree 1
    transfer solution false


Cycle 0:
   Number of active cells: 512
   Number of degrees of freedom: 2187


  Newton iteration 1
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 173.076 -> 1.64265e-06 in 7 Bicgstab iterations.
      Accepting Newton solution with residual: 1.64265e-06


   Newton iteration 2
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 57.3622 -> 3.23721e-07 in 8 Bicgstab iterations.
      Accepting Newton solution with residual: 24.9028
      Active set did not change!


   Newton iteration 3
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 24.9028 -> 9.94326e-08 in 7 Bicgstab iterations.
      Residual of the non-contact part of the system: 1.63333
         with a damping parameter alpha = 1
      Active set did not change!


...


  Newton iteration 6
      Updating active set...
         Size of active set: 1
      Assembling system...
      Solving system...
         Error: 1.43188e-07 -> 3.56218e-16 in 8 Bicgstab iterations.
      Residual of the non-contact part of the system: 4.298e-14
         with a damping parameter alpha = 1
      Active set did not change!
      Writing graphical output... p1_adaptive/solution-00.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      1.13s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assembling                      |         6 |     0.463s |        41% |
| Graphical output                |         1 |    0.0257s |       2.3% |
| Residual and lambda             |         4 |    0.0754s |       6.7% |
| Setup                           |         1 |     0.227s |        20% |
| Setup: constraints              |         1 |    0.0347s |       3.1% |
| Setup: distribute DoFs          |         1 |    0.0441s |       3.9% |
| Setup: matrix                   |         1 |    0.0119s |       1.1% |
| Setup: vectors                  |         1 |   0.00155s |      0.14% |
| Solve                           |         6 |     0.246s |        22% |
| Solve: iterate                  |         6 |    0.0631s |       5.6% |
| Solve: setup preconditioner     |         6 |     0.167s |        15% |
| update active set               |         6 |    0.0401s |       3.6% |
+---------------------------------+-----------+------------+------------+


Peak virtual memory used, resident in kB: 541884 77464
Contact force = 37.3058


...


Cycle 3:
   Number of active cells: 14652
   Number of degrees of freedom: 52497


   Newton iteration 1
      Updating active set...
         Size of active set: 145
      Assembling system...
      Solving system...
         Error: 296.309 -> 2.72484e-06 in 10 Bicgstab iterations.
      Accepting Newton solution with residual: 2.72484e-06


...


   Newton iteration 10
      Updating active set...
         Size of active set: 145
      Assembling system...
      Solving system...
         Error: 2.71541e-07 -> 1.5428e-15 in 27 Bicgstab iterations.
      Residual of the non-contact part of the system: 1.89261e-13
         with a damping parameter alpha = 1
      Active set did not change!
      Writing graphical output... p1_adaptive/solution-03.pvtu



+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    |      38.4s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assembling                      |        10 |      22.5s |        58% |
| Graphical output                |         1 |     0.327s |      0.85% |
| Residual and lambda             |         9 |      3.75s |       9.8% |
| Setup                           |         1 |      4.83s |        13% |
| Setup: constraints              |         1 |     0.578s |       1.5% |
| Setup: distribute DoFs          |         1 |      0.71s |       1.8% |
| Setup: matrix                   |         1 |     0.111s |      0.29% |
| Setup: refine mesh              |         1 |      4.83s |        13% |
| Setup: vectors                  |         1 |   0.00548s |     0.014% |
| Solve                           |        10 |      5.49s |        14% |
| Solve: iterate                  |        10 |       3.5s |       9.1% |
| Solve: setup preconditioner     |        10 |      1.84s |       4.8% |
| update active set               |        10 |     0.662s |       1.7% |
+---------------------------------+-----------+------------+------------+


Peak virtual memory used, resident in kB: 566052 105788
Contact force = 56.794


...
@endcode



每个周期结束时的表格显示了最近一次网格细化周期的计算时间（这些数字当然是针对产生该输出的机器而言的）和程序不同部分的调用次数，如装配或计算残差。上面的一些数字可以通过将解决方案从一个网格转移到下一个网格来改善，我们在这里没有行使这个选项。当然，你也可以通过使用更多的处理器来使程序运行得更快，特别是在后期的细化周期中：附带的论文显示，至少有1000个内核的良好扩展性。

在一个典型的运行中，你可以看到，对于每一个细化步骤，活动集--接触点--首先被迭代出来。之后，牛顿方法只需要解决塑性问题。对于更细的网格，在最后4或5次牛顿迭代中可以看到二次收敛。

我们不会在这里详细讨论每个输入文件的情况。相反，让我们只展示解决方案的图片（如果单元格的正交点为零，塑性不等式处于活动状态，则域的左半部分被省略）。

 <table align="center">
  <tr>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionColorbar.png">
    </td>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionBall2.png" alt="" width="70%">
    </td>
    <td valign="top">
      &nbsp;
    </td>
    <td>
    <img src="https://www.dealii.org/images/steps/developer/step-42.CellConstitutionLi2.png" alt="" alt="" width="70%">
    </td>
  </tr>
</table> 

图中显示了适应性细化以及细胞在与球接触过程中的塑化程度。请记住，我们考虑每个正交点的应力偏差部分的规范，以查看是否有弹性或塑性行为。蓝色意味着这个单元只包含弹性正交点，与所有正交点都被塑化的红色单元相反。在顶面的中间--网格最细的地方--非常仔细地看可以看到由障碍物引起的凹陷。这是 <code>move_mesh()</code> 函数的结果。然而，由于我们在这里考虑的障碍物的压痕非常小，所以很难辨别这种效果；我们可以玩玩将网格的顶点按计算出的位移的倍数进行位移。

关于使用该程序可以获得的结果的进一步讨论，见本页面最上方提到的出版物。


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h1>Possibilities for extensions</h1>


像往常一样，有多种可能性来扩展这个程序。从算法的角度来看，这个程序在写作时已经达到了我们所能达到的程度，使用了接触不等式、塑性非线性和线性求解器的最佳可用算法。然而，就更现实的情况而言，人们希望用这个程序做一些事情。   <ul>   <li>  将程序从静态扩展到准静态情况，也许可以通过选择后向欧拉模式来实现时间离散化。一些理论结果可以在Jörg Frohne的博士论文中找到，<i>FEM-Simulation
der Umformtechnik metallischer Oberfl&auml;chen im Mikrokosmos</i>，德国锡根大学，2011。

 <li> 考虑有摩擦力的接触问题也将是一个有趣的进步。在几乎每个机械过程中，摩擦都有很大的影响。  为了模拟这种情况，我们必须考虑到接触面的切向应力。摩擦也给我们的问题增加了另一个不等式，因为只要切向应力不超过某个极限，身体和障碍物通常会粘在一起，超过这个极限，两个身体就会互相滑过。

 <li>  如果我们已经模拟了摩擦性接触，下一步要考虑的是接触区的发热。由两个物体之间的摩擦引起的热量会提高可变形物体的温度，并导致一些材料参数的变化。

 <li>  对于接触以及塑性，实施更精确的、与问题相适应的误差估计器可能是有意义的。   </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-42.cc"
*/
