/**
@page step_43 The step-43 tutorial program
This tutorial depends on step-31.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Advectiondominatedtwophaseflowmathematicalmodel">Advection-dominated two-phase flow mathematical model.</a>
        <li><a href="#Adaptiveoperatorsplittingandtimestepping">Adaptive operator splitting and time stepping.</a>
        <li><a href="#Timediscretization">Time discretization.</a>
        <li><a href="#Weakformspacediscretizationforthepressurevelocitypart">Weak form, space discretization for the pressure-velocity part.</a>
        <li><a href="#Stabilizationweakformandspacediscretizationforthesaturationtransportequation">Stabilization, weak form and space discretization for the saturation transport equation.</a>
        <li><a href="#Adaptivemeshrefinement">Adaptive mesh refinement.</a>
        <li><a href="#Linearsystemanditspreconditioning">Linear system and its preconditioning.</a>
        <li><a href="#Thetestcases">The test cases.</a>
        <li><a href="#Listofreferences">List of references</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Boundaryandinitialvalueclasses">Boundary and initial value classes</a>
        <li><a href="#Permeabilitymodels">Permeability models</a>
        <li><a href="#Physicalquantities">Physical quantities</a>
        <li><a href="#Helperclassesforsolversandpreconditioners">Helper classes for solvers and preconditioners</a>
        <li><a href="#TheTwoPhaseFlowProblemclass">The TwoPhaseFlowProblem class</a>
        <li><a href="#TwoPhaseFlowProblemdimTwoPhaseFlowProblem">TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem</a>
        <li><a href="#TwoPhaseFlowProblemdimsetup_dofs">TwoPhaseFlowProblem<dim>::setup_dofs</a>
        <li><a href="#Assemblingmatricesandpreconditioners">Assembling matrices and preconditioners</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemdimassemble_darcy_preconditioner">TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner</a>
        <li><a href="#TwoPhaseFlowProblemdimbuild_darcy_preconditioner">TwoPhaseFlowProblem<dim>::build_darcy_preconditioner</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_darcy_system">TwoPhaseFlowProblem<dim>::assemble_darcy_system</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_system">TwoPhaseFlowProblem<dim>::assemble_saturation_system</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_matrix">TwoPhaseFlowProblem<dim>::assemble_saturation_matrix</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs_cell_term">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term</a>
        <li><a href="#TwoPhaseFlowProblemdimassemble_saturation_rhs_boundary_term">TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term</a>
      </ul>
        <li><a href="#TwoPhaseFlowProblemdimsolve">TwoPhaseFlowProblem<dim>::solve</a>
        <li><a href="#TwoPhaseFlowProblemdimrefine_mesh">TwoPhaseFlowProblem<dim>::refine_mesh</a>
        <li><a href="#TwoPhaseFlowProblemdimoutput_results">TwoPhaseFlowProblem<dim>::output_results</a>
        <li><a href="#Toolfunctions">Tool functions</a>
      <ul>
        <li><a href="#TwoPhaseFlowProblemdimdetermine_whether_to_solve_for_pressure_and_velocity">TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity</a>
        <li><a href="#TwoPhaseFlowProblemdimproject_back_saturation">TwoPhaseFlowProblem<dim>::project_back_saturation</a>
        <li><a href="#TwoPhaseFlowProblemdimget_max_u_F_prime">TwoPhaseFlowProblem<dim>::get_max_u_F_prime</a>
        <li><a href="#TwoPhaseFlowProblemdimget_extrapolated_saturation_range">TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range</a>
        <li><a href="#TwoPhaseFlowProblemdimcompute_viscosity">TwoPhaseFlowProblem<dim>::compute_viscosity</a>
      </ul>
        <li><a href="#TwoPhaseFlowProblemdimrun">TwoPhaseFlowProblem<dim>::run</a>
        <li><a href="#Thecodemaincodefunction">The <code>main()</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-43/doc/intro.dox

 <br> 

<i>
This program was contributed by Chih-Che Chueh (University of Victoria) and
Wolfgang Bangerth. Results from this program are used and discussed in the
following publications (in particular in the second one):


- Chih-Che Chueh, Marc Secanell, Wolfgang Bangerth, Ned Djilali. Multi-level
  adaptive simulation of transient two-phase flow in heterogeneous porous
  media. Computers &amp; Fluids, 39:1585-1596, 2010


- Chih-Che Chueh, Ned Djilali, Wolfgang Bangerth. An h-adaptive operator
  splitting method for two-phase flow in 3D heterogeneous porous
  media. SIAM Journal on Scientific Computing, 35:B149-B175, 2013.


The implementation discussed here uses and extends
parts of the step-21 and step-31 tutorial programs.


The work of the Chih-Che Chueh was funded through the Canada Research Chairs
Program and the MITACS Network of Centres of Excellence. Parts of the work by
Wolfgang Bangerth were funded through Award No. KUS-C1-016-04, made by the King
Abdullah University of Science and Technology, and through an Alfred P. Sloan
Research Fellowship.
This material is also in parts based upon work supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology; and in a continuation by the National Science
Foundation under Award No. EAR-0949446 and The University of California
&ndash; Davis. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation, The
California Institute of Technology, or of The University of California
&ndash; Davis.
</i>


<a name="Introduction"></a><a name="Intro"></a><h1>Introduction</h1> 。


多孔介质中的多相流模拟是一个无处不在的问题，我们以前在步骤20和步骤21中已经以某种形式解决了这个问题。然而，正如在那里很容易看到的那样，它面临两个主要困难：数值精度和效率。第一个问题在第20步的静止求解器中很容易看到：使用最低阶的Raviart-Thomas元素不可能产生高度精确的解决方案。我们需要更精确的方法。第二个原因从时间相关的步骤-21中可以看出：该程序慢得令人发指，没有希望在合理的时间范围内得到高度准确的三维解。

在这个项目中，为了克服这两个问题，有五个方面我们正在努力改进，以实现高性能的模拟器。

 <ul>   <li>  高阶空间离散  <li>  自适应网格细化  <li>  自适应时间步进  <li>  运算器分割  <li>  高效求解器和预处理  </ul> 

这个计划的大部分灵感来自第31步，但这里讨论的几个技术是原创的。




<a name="Advectiondominatedtwophaseflowmathematicalmodel"></a><h3>Advection-dominated two-phase flow mathematical model.</h3>


我们考虑的是两相不相溶的不可压缩流体的流动。毛细管和重力效应被忽略了，粘性效应被假定为主导。这种流动的管理方程与步骤21中使用的方程相同，为

@f{align*}
  \mathbf{u}_t &= - \mathbf{K} \lambda_t \left(S\right) \nabla p, \\
  \nabla \cdot \mathbf{u}_t &= q, \\
  \epsilon \frac{\partial S}{\partial t} + \nabla \cdot \left( \mathbf{u}_t  F\left( S \right) \right)&=0,


@f}

其中 $S$ 是第二（润湿）相的饱和度（体积分数在零和一之间）， $p$ 是压力， $\mathbf{K}$ 是渗透率张量， $\lambda_t$ 是总流动性， $\epsilon$ 是孔隙度， $F$ 是湿润相的分流量， $q$ 是源项， $\mathbf{u}_t$ 是总速度。总流动性、润湿相的部分流量和总速度分别由以下公式给出

@f{align*}
   \lambda_t(S)&= \lambda_w + \lambda_{nw} = \frac{k_{rw}(S)}{\mu_w} + \frac{k_{rnw}(S)}{\mu_{nw}}, \\
   F(S) &= \frac{\lambda_w}{\lambda_t} = \frac{\lambda_w}{\lambda_w + \lambda_{nw}} = \frac{k_{rw}(S)/\mu_w}{k_{rw}(S)/\mu_w + k_{rnw}(S)/\mu_{nw}}, \\
   \mathbf{u}_t &= \mathbf{u}_w + \mathbf{u}_{nw} = -\lambda_t(S)\mathbf{K} \cdot \nabla p,


@f}

其中下标 $w, nw$ 分别代表湿润和非湿润阶段。

为方便起见，饱和度方程中的孔隙度 $\epsilon$ 可被视为时间变量的比例系数，被设定为1。根据相对渗透率 $k_{rw}$ 和 $k_{rnw}$ 对饱和度的依赖性的常用规定，我们用

@f{align*}
   k_{rw}  &= S^2, \qquad&\qquad
   k_{rnw} &= \left( 1-S \right)^2.


@f}



上面的多孔介质方程由饱和度的初始条件和压力的边界条件来补充。由于饱和度和压力梯度唯一地决定了速度，所以速度的边界条件是没有必要的。由于流动方程不包含时间导数，因此不需要速度和压力变量的初始条件。流场将边界分为流入或流出部分。具体来说。

@f[
   \mathbf{\Gamma}_{in}(t) = \left\{\mathbf{x} \in \partial \Omega:\mathbf{n} \cdot \mathbf{u}_t<0\right\},


@f]

我们通过在流入边界上施加饱和变量的边界值，得出一个完整的模型  $\mathbf{\Gamma}_{in}$  。




<a name="Adaptiveoperatorsplittingandtimestepping"></a><h3>Adaptive operator splitting and time stepping.</h3>


从第21步可以看出，一旦我们知道了流量变量，求解速度和压力的流量方程是程序中花费时间远大于饱和度变量的（明确）更新步骤的部分。另一方面，压力和速度对饱和度的依赖性很弱，因此可以考虑每隔几步只求解压力和速度，而每步更新饱和度。如果我们能找到一个关于何时需要更新流量变量的标准，我们把这种拆分称为 "自适应算子拆分 "方案。

在这里，我们使用以下后验标准来决定何时重新计算压力和速度变量（详细的推导和描述可以在[Chueh, Djilali and Bangerth 2011]中找到）。

@f{align*}
  \theta(n,n_p)
  =
    \max_{\kappa\in{\mathbb T}}
    \left(
    \left\|
      \frac 1{\lambda_t\left(S^{(n-1)}\right)}


      - \frac 1{\lambda_t\left(S^{(n_p)}\right)} \right\|_{L^\infty(\kappa)}
    \left\|\|\mathbf{K}^{-1}\|_1\right\|_{L^\infty(\kappa)}
    \right).


@f}

其中括号内的上标表示定义任何数量的饱和时间步数， $n_p<n$ 代表我们实际计算压力和速度的最后一步。如果 $\theta(n,n_p)$ 超过某个阈值，我们就重新计算流量变量；否则，我们在时间步骤 $n$ 中跳过这个计算，只将饱和变量向前移动一个时间步骤。

简而言之，该算法允许我们执行若干长度为 $\Delta t_c^{(n)}=t^{(n)}_c-t^{(n-1)}_c$ 的饱和时间步长，直到上述标准告诉我们重新计算速度和压力变量，导致一个长度为

@f[
   \Delta t_p^{(n)} = \sum_{i=n_p+1}^{n} \Delta t_c^{(i)}.


@f]

我们根据Courant-Friedrichs-Lewy（CFL）限制来选择（微型）步骤的长度，标准是

@f[
  \Delta t_c = \frac{\textrm{min}_{K}h_{K}}{7 \|\mathbf{u}_t\|_{L^{\infty}\left(\Omega\right)}},


@f]

我们已经证实，对于下面讨论的饱和方程的有限元和时间步长方案的选择是稳定的（ $h_K$ 表示单元 $K$ 的直径）。其结果是一个方案，微观和宏观的时间步长都不统一，两者都是自适应选择。

<a name="Timediscretization"></a><h3>Time discretization.</h3> 利用这种时间离散化，我们从IMPES方法中得到每个时间步骤的以下方程组（见步骤21）。


@f{align*}
   \mathbf{u}^{(n)}_t + \lambda_t\left(S^{(n-1)}\right) \mathbf{K} \nabla p^{(n)} =0, \\
   \nabla \cdot \mathbf{u}^{(n)}_t = q, \\
   \epsilon \left( \frac{S^{(n-1)}-S^{(n)}}{\Delta t^{(n)}_c} \right) + \mathbf{u}^{(n)}_t \cdot \nabla F\left(S^{(n-1)}\right) + F\left(S^{(n-1)}\right) \nabla \cdot \mathbf{u}^{(n)}_t =0.


@f}




利用 $\nabla \cdot \mathbf{u}_t = q$ 这一事实，时间离散的饱和度方程变为

@f{align*}
  &\epsilon \left( \frac{S^{(n)}-S^{(n-1)}}{\Delta t^{(n)}_c} \right) + \mathbf{u}^{(n)}_t \cdot \nabla F\left(S^{(n-1)}\right) + F\left(S^{(n-1)}\right)q=0.


@f}



<a name="Weakformspacediscretizationforthepressurevelocitypart"></a><h3>Weak form, space discretization for the pressure-velocity part.</h3>


通过将定义总速度的方程 $\mathbf u_t^{(n)}$ 和用源项表示其发散的方程分别与测试函数 $\mathbf{v}$ 和 $w$ 相乘，然后根据需要进行分项积分，问题的弱形式为。找出 $\mathbf u, p$ ，以便对所有测试函数 $\mathbf{v}, w$ 而言，存在

@f{gather*}
   \left( \left( \mathbf{K} \lambda_t\left(S^{(n-1)}\right) \right)^{-1} \mathbf{u}^{(n)}_t, \mathbf{v}\right)_{\Omega} - \left(p^{(n)}, \nabla \cdot \mathbf{v}\right)_{\Omega} = -\left(p^{(n)}, \mathbf{n} \cdot \mathbf{v} \right)_{\partial \Omega}, \\


   - \left( \nabla \cdot \mathbf{u}^{(n)}_t,w\right)_{\Omega} = - \big(q,w\big)_{\Omega}.


@f}

这里， $\mathbf{n}$ 代表 $\partial
\Omega$ 的单位外向法向量，压力 $p^{(n)}$ 可以在边界 $\partial \Omega$ 的开放部分弱化规定，而在那些规定了速度的部分（例如具有 $\mathbf n \cdot \mathbf
u=0$ 的不渗透边界，该术语完全消失了，因为 $\mathbf n \cdot \mathbf
v=0$  。

我们使用连续有限元来离散速度和压力方程。具体来说，我们使用混合有限元来确保同时对矢量变量（如流体速度）和标量变量（如压力）进行高阶逼近。对于鞍点问题，公认的是需要满足所谓的Babuska-Brezzi或Ladyzhenskaya-Babuska-Brezzi（LBB）条件[Brezzi 1991, Chen 2005]以确保压力-速度系统的稳定性。在本工作中，通过使用比压力高一阶的速度元素，即 $u_h \in Q^d_{p+1}$ 和 $p_h \in Q_p$ 来满足这些稳定性条件，其中 $p=1$ ， $d$ 是空间维度， $Q_s$ 表示每个变量的张量积Lagrange多项式的空间 $s$ 。

<a name="Stabilizationweakformandspacediscretizationforthesaturationtransportequation"></a><h3>Stabilization, weak form and space discretization for the saturation transport equation.</h3>为饱和方程选择的 $Q_1$ 元素在没有上卷或其他类型的稳定化的情况下不会导致稳定的离散化，并且在数值解中会出现虚假的震荡。添加一个人工扩散项是消除这些振荡的一种方法[Chen 2005]。另一方面，添加过多的扩散项会在解中涂抹出尖锐的锋面，并且会出现网格定向困难[Chen 2005]。为了避免这些影响，我们使用了由[Guermond和Pasquetti 2008]提出并在[Chueh, Djilali, Bangerth 2011]和[Kronbichler, Heister and Bangerth, 2011]以及步骤31中验证的人工扩散项。


这种方法修改了饱和度方程的（离散）弱形式，改为

@f{align*}
  \left(\epsilon \frac{\partial S_h}{\partial t},\sigma_h\right)


  -
  \left(\mathbf{u}_t  F\left( S_h \right),
    \nabla \sigma_h\right)
  +
  \left(\mathbf n \cdot \mathbf{u}_t  \hat F\left( S_h \right),
    \sigma_h\right)_{\partial\Omega}
  +
  (\nu(S_h) \nabla S_h, \nabla \sigma_h)
  &=0
  \qquad
  \forall \sigma_h,


@f}

其中 $\nu$ 是人工扩散参数， $\hat F$ 是域的边界上适当选择的数值通量（我们为此选择明显的全上风通量）。

根据[Guermond and Pasquetti 2008]（以及[Chueh, Djilali and Bangerth 2011]中的详细说明），我们将参数作为一个片状常数函数，设置在直径为 $K$ 的每个单元上，为

@f[
   \nu(S_h)|_{K} = \beta \| \mathbf{u}_t \max\{F'(S_h),1\} \|_{L^{\infty}(K)} \textrm{min} \left\{ h_{K},h^{\alpha}_{K} \frac{\|\textrm{Res}(S_h)\|_{L^{\infty}(K)}}{c(\mathbf{u}_t,S)} \right\}


@f]

其中 $\alpha$ 为稳定化指数， $\beta$ 为用户定义的无量纲稳定化常数。按照[Guermond和Pasquetti 2008]以及步骤31的实现，速度和饱和度全局归一化常数 $c(\mathbf{u}_t,S)$ 和残差 $\textrm{Res}(S)$ 分别为

@f[
   c(\mathbf{u}_t,S) = c_R \|\mathbf{u}_t \max\{F'(S),1\}\|_{L^{\infty}(\Omega)} \textrm{var}(S)^\alpha | \textrm{diam} (\Omega) |^{\alpha - 2}


@f]

和

@f[
   \textrm{Res}(S) = \left( \epsilon \frac{\partial S}{\partial t} + \mathbf{u}_t \cdot \nabla F(S) + F(S)q \right) \cdot S^{\alpha - 1}


@f]

其中 $c_R$ 是用户定义的第二个无维常数， $\textrm{diam}(\Omega)$ 是域的直径， $\textrm{var}(S) =
\textrm{max}_{\Omega} S - \textrm{min}_{\Omega} S$ 是整个计算域中目前饱和值的范围 $\Omega$  。

这种稳定方案与更简单的方案，如有限体积（或不连续Galerkin）方法或流线型上风Petrov Galerkin（SUPG）离散法相比有很多优点。特别是，人工扩散项主要作用于不连续点附近，因为在饱和度平稳的地区，残差很小。因此，它提供了一个更高的精度。另一方面，它是非线性的，因为  $\nu$  取决于饱和度  $S$  。我们通过明确处理所有的非线性项来避免这一困难，这导致了以下时间步长的完全离散问题  $n$  。

@f{align*}
   &\left( \epsilon S_h^{(n)},\sigma_h\right)_{\Omega} - \Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\mathbf{u}^{*}_t,\nabla\sigma_h\Big)_{\Omega} + \Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{*}_t\right),\sigma_h\Big)_{\partial\Omega} \nonumber \\
   & \quad = \left( \epsilon S_h^{(n-1)},\sigma_h\right)_{\Omega} - \Delta t^{(n)}_c \bigg(\nu\left(S_h^{(n-1)}\right)\nabla S_h^{(n-1)},\nabla\sigma_h\bigg)_{\Omega} \nonumber \\
   & \qquad + \Delta t^{(n)}_c \bigg(\mathbf{n}\cdot\nu\left(S_h^{(n-1)}\right)\nabla S^{(n-1)},\sigma_h\bigg)_{\partial\Omega}


@f}

其中 $\mathbf{u}_t^{*}$ 是从 $\mathbf{u}^{(n_p)}_t$ 和 $\mathbf{u}^{(n_{pp})}_t$ 线性外推到当前时间 $t^{(n)}$ 的速度，如果 $\theta<\theta^*$ ，而 $\mathbf{u}_t^{*}$ 是 $\mathbf{u}^{(n_p)}_t$ ，如果 $\theta>\theta^*$  。因此，该方程在 $S_h^{(n)}$ 中是线性的，所需要的是用饱和空间上的质量矩阵来解决。

由于饱和度的Dirichlet边界条件只施加在流入边界上，所以上述方程左边的第三个项需要进一步分成两部分。

@f{align*}
  &\Delta t^{(n)}_c \Big(F\left(S_h^{(n-1)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_t\right),\sigma_h\Big)_{\partial\Omega} \nonumber \\
  &\qquad= \Delta t^{(n)}_c \Big(F\left(S^{(n-1)}_{(+)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_{t(+)}\right),\sigma_h\Big)_{\partial\Omega_{(+)}} + \Delta t^{(n)}_c \Big(F\left(S^{(n-1)}_{(-)}\right)\left(\mathbf{n}\cdot\mathbf{u}^{(n)}_{t(-)}\right),\sigma_h\Big)_{\partial\Omega_{(-)}}


@f}

其中 $\partial\Omega_{(-)} = \left\{\mathbf{x} \in \partial\Omega : \mathbf{n}
  \cdot \mathbf{u}_t<0\right\}$ 和 $\partial\Omega_{(+)} = \left\{\mathbf{x} \in \partial\Omega : \mathbf{n} \cdot
  \mathbf{u}_t>0\right\}$ 分别代表流入和流出的边界。我们使用上风公式选择数值，即 $S^{(n-1)}_{(+)}$ 和 $\mathbf{u}^{(n)}_{t(+)}$ 对应于从当前单元中提取的数值，而 $S^{(n-1)}_{(-)}$ 和 $\mathbf{u}^{(n)}_{t(-)}$ 的数值是来自邻近的边界 $\partial\Omega_{(-)}$ 。




<a name="Adaptivemeshrefinement"></a><h3>Adaptive mesh refinement.</h3>


适应性地选择网格以解决尖锐的饱和前沿是我们算法中实现效率的一个基本要素。在这里，我们使用[Chueh, Djilali and Bangerth 2011]中使用的相同的冲击型细化方法来选择那些应该被细化或粗化的单元。三角形的每个单元 $K$ 的细化指标是通过以下方式计算的

@f[
   \eta_{K} = |\nabla S_h(\mathbf x_K)|


@f]

其中 $\nabla S_h(\mathbf x_K)$ 是在 $\mathbf x_K$ 单元的中心评价的离散饱和变量的梯度。这种方法类似于可压缩流动问题中经常使用的方法，即用密度梯度来表示细化。也就是说，正如我们将在<a href="#Results">results section</a>的结尾处讨论的那样，这被证明不是一个非常有用的标准，因为它基本上到处都导致细化。我们在这里只是为了说明问题而展示它。




<a name="Linearsystemanditspreconditioning"></a><h3>Linear system and its preconditioning.</h3>


按照上面讨论的治理方程的离散化，我们得到一个时间步长为 $(n)$ 的线性方程组，形式如下。

@f[
 \left(
  \begin{array}{ccc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T} & \mathbf{0}  \\
   \mathbf{B}           & \mathbf{0}     & \mathbf{0}   \\
   \mathbf{H}           & \mathbf{0}     & \mathbf{M}^{S}
  \end{array}
 \right)
 \left(
  \begin{array}{c}
   \mathbf{U}^{(n)} \\
   \mathbf{P}^{(n)} \\
   \mathbf{S}^{(n)}
  \end{array}
 \right)
 =
 \left(
  \begin{array}{c}
   0 \\
   \mathbf{F}_{2} \\
   \mathbf{F}_{3}
  \end{array}
 \right)


@f]

其中各个矩阵和向量的定义如下，使用形状函数 $\mathbf{v}_i$ 表示速度， $\phi_i$ 表示压力和饱和度。

@f{align*}
  \mathbf{M}^{\mathbf{u}}_{ij}
  &= \left( \left( \mathbf{K} \lambda_t\left(S^{(n-1)}\right) \right)^{-1}
  \mathbf{v}_{i},\mathbf{v}_{j}\right)_{\Omega},
  &
  \mathbf{M}^{S}_{ij}           &= \left(\epsilon \phi_i,\phi_j\right)_{\Omega}
  \\
  \mathbf{B}_{ij}
  &= - \left( \nabla \cdot \mathbf{v}_{j},\phi_{i}\right)_{\Omega},
  &
  \mathbf{H}_{ij}
  &= - \Delta t^{(n)}_c \Big( F\left(S^{(n-1)}\right) \mathbf{v}_i,\nabla\phi_j\Big)_{\Omega}
  \\
  \left(\mathbf{F}_{2}\right)_i
  &= - \big(F\left(S^{(n-1)}\right)q,\phi_i\big)_{\Omega},


@f}

和 $\mathbf{F}_{3}$ 在稳定传输方程的定义中给出。

如果我们把左上角的 $2\times 2$ 板块的矩阵视为一个板块，那么上面的线性系统是块状三角形形式。因此，我们可以首先求解速度和压力（除非我们决定用 $\mathbf U^{(n_p)}$ 来代替速度），然后再求解饱和度变量。其中第一个步骤要求我们解决

@f[
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T}  \\
   \mathbf{B}           & \mathbf{0}
  \end{array}
 \right)
 \left(
  \begin{array}{c}
   \mathbf{U}^{(n)} \\
   \mathbf{P}^{(n)}
  \end{array}
 \right)
 =
 \left(
  \begin{array}{c}
   0 \\
   \mathbf{F}_{2}
  \end{array}
 \right)


@f]

我们对这个线性系统采用广义最小残差（GMRES）方法[Saad和Schultz 1986]。速度-压力系统的理想预处理方法是

@f{align*}
\mathbf{P} =
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} &  \mathbf{0}  \\
   \mathbf{B}           & -\mathbf{S}
  \end{array}
 \right),
 & \qquad
 \mathbf{P}^{-1} =
 \left(
  \begin{array}{cc}
   \left(\mathbf{M}^{\mathbf{u}}\right)^{-1}                              &  \mathbf{0}  \\
   \mathbf{S}^{-1} \mathbf{B} \left(\mathbf{M}^{\mathbf{u}}\right)^{-1}   & -\mathbf{S}^{-1}
  \end{array}
 \right)
 @f}

其中 $\mathbf{S}=\mathbf{B}\left(\mathbf{M}^{\mathbf{u}}\right)^{-1}\mathbf{B}^T$ 是系统的Schur补充[Zhang 2005]。这个预处理程序是最优的，因为

@f{align*}
 \mathbf{P}^{-1}
 \left(
  \begin{array}{cc}
   \mathbf{M}^{\mathbf{u}} & \mathbf{B}^{T}  \\
   \mathbf{B}           & \mathbf{0}
  \end{array}
 \right)
 =
  \left(
  \begin{array}{cc}
   \mathbf{I}         &  \left(\mathbf{M}^{\mathbf{u}}\right)^{-1} \mathbf{B}^{T}  \\
   \mathbf{0}         &  \mathbf{I}
  \end{array}
 \right),


@f}

对其而言，可以证明GMRES在两次迭代中收敛。

然而，我们当然不能指望使用速度质量矩阵和Schur补数的精确求逆。因此，我们采用[Silvester and Wathen 1994]最初为斯托克斯系统提出的方法。将其适用于当前的方程组，得到预处理程序

@f{align*}
 \mathbf{\tilde{P}}^{-1} =
 \left(
  \begin{array}{cc}
   \widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}
                              &  \mathbf{0}  \\
   \widetilde{\mathbf{{S}}^{-1}} \mathbf{B} \widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}   & -\widetilde{\mathbf{{S}}^{-1}}
  \end{array}
 \right)


@f}

其中蒂尔德表示精确逆矩阵的近似值。特别是，由于 $\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}=\left( \left(
    \mathbf{K} \lambda_t \right)^{-1}
  \mathbf{v}_{i},\mathbf{v}_{j}\right)_{\Omega}$ 是一个稀疏的对称和正定矩阵，我们为 $\widetilde{\left(\mathbf{{M}}^{\mathbf{u}}\right)^{-1}}$ 选择了这个矩阵的稀疏不完全Cholesky分解的单一应用[Golub和Van Loan 1996]。我们注意到，对应于非混合形式的多孔介质流动算子的舒尔补， $-\nabla \cdot [\mathbf K
\lambda_t(S)]\nabla$ 和 $\mathbf{\tilde {S}} = \left( \left( \mathbf{K} \lambda_t \right) \nabla \phi_{i},\nabla \phi_{j}\right)_{\Omega}$ 应该是实际舒尔补矩阵 $\mathbf
S$ 的良好近似。由于这两个矩阵又都是对称和正定的，所以我们用 $\mathbf{\tilde S}$ 的不完全Cholesky分解来表示 $\widetilde
{\mathbf{{S}}^{-1}}$ 。需要注意的是， $\mathbf{\tilde S}$ 需要用Dirichlet边界条件建立，以确保其可逆性。

一旦有了速度 $\mathbf{U}^{(n)} \equiv \mathbf{u}^*_t$ ，我们就可以把 $\mathbf{H}$ 和 $\mathbf{F}_{3}$ 组合起来，用以下方法解决饱和度的问题

@f{align*}
  \mathbf{M}^{S} \mathbf{S}^{(n)} = \mathbf{F}_{3} - \mathbf{H} \mathbf{U}^{(n)}.


@f}

其中质量矩阵 $\mathbf{M}^{S}$ 用共轭梯度法求解，再一次使用不完全的Cholesky分解作为预处理。

<a name="Thetestcases"></a><h3>The test cases.</h3>


 @note  这里讨论的实现使用并扩展了这个库的步骤21、步骤31和步骤33教程的部分程序。特别是，如果你想了解它是如何工作的，请参考step-21关于数学问题的讨论，以及step-31，大部分的实现都来自于此。我们将不讨论在步骤31中已经讨论过的实现的各个方面。

我们展示了一些两相流方程的数值结果，这些方程通过适当的初始和边界条件，结合两种不同的渗透率模型的选择而得到增强。在所考虑的问题中，没有内部源项（ $q=0$ ）。如上所述，定量的数值结果在[Chueh, Djilali and Bangerth 2011]中提出。

为了简单起见，我们选择了 $\Omega=[0,1]^d,d=2,3$ ，尽管所有的方法（以及我们的实现）在一般的非结构化网格上都应该同样工作。

初始条件只需要饱和变量，我们选择 $S(\mathbf{x},0)=0.2$ ，即多孔介质最初是由非湿润（80%）和湿润（20%）相的混合物填充。这与步骤21中的初始条件不同，在该步骤中我们采用了 $S(\mathbf{x},0)=0$ ，但由于复杂的数学原因，在那里的长篇评论中提到，目前使用基于熵的人工扩散项的方法在不对方法进行额外修改的情况下不能收敛到这个初始条件的粘度解。因此，我们在目前的计划中选择了这个修改过的版本。

此外，我们在边界上规定了一个线性压力。

@f[
   p(\mathbf{x},t) = 1 - x \qquad
   \textrm{on} \quad \partial \Omega \times [0,T].


@f]

压力和饱和度唯一地决定了速度，而速度决定了一个边界段是流入还是流出的边界。在边界的流入部分， $\mathbf{\Gamma}_{in}(t)$ ，我们规定

@f{align*}
   S(\mathbf{x},t) = 1 \qquad & \textrm{on} \quad \mathbf{\Gamma}_{in}(t) \cap \left\{x = 0\right\}, \\
   S(\mathbf{x},t) = 0 \qquad & \textrm{on} \quad \mathbf{\Gamma}_{in}(t) \backslash \left\{x = 0\right\}.


@f}

换句话说，该领域被来自左边的湿润相淹没。对于边界的流出部分，不需要饱和的边界条件。

所有用于二维/三维案例的数值和物理参数都列在下表中。

 <table align="center" class="tutorial" width="50%">
<tr>
    <th>Parameter                           </th><th>Symbol          </th><th>Value               </th><th>units     </th></tr><tr>
    <td>Porosity                            </td><td>$\epsilon$      </td><td>1.0                 </td><td>-                   </td></tr><tr>
    <td>Viscosity (wetting)                 </td><td>$\mu_w$         </td><td>0.2                 </td><td>$kg \cdot m^{-1} \cdot sec^{-1}$   </td></tr><tr>
    <td>Viscosity (nonwetting)              </td><td>$\mu_{nw}$      </td><td>1.0                 </td><td>$kg \cdot m^{-1} \cdot sec^{-1}$      </td></tr><tr>
    <td>Stabilization exponent              </td><td>$\alpha$        </td><td>1.0                 </td><td>-     </td></tr><tr>
    <td>Stabilization constant              </td><td>$\beta$         </td><td>2D: 0.3; 3D: 0.27   </td><td>- </td></tr><tr>
    <td>Normalization constant              </td><td>$c_R$           </td><td>1.0                 </td><td>- </td></tr><tr>
    <td>Number of high-permeability regions </td><td>$N$             </td><td>50; 200             </td><td>- </td></tr><tr>
    <td>Operator splitting threshold        </td><td>$\theta^\ast$   </td><td>5.0              </td><td>- </td></tr>
</table> 




<a name="Listofreferences"></a><h3>List of references</h3>



<ol>  <li>  CC Chueh, N Djilali and W Bangerth.   <br>  三维异质多孔介质中两相流的h-适应性算子分割方法。   <br>  SIAM科学计算杂志，第35卷（2013），第B149-B175页

 <li>  M. Kronbichler, T. Heister, and W. Bangerth  <br>  通过现代数值方法进行高精度地幔对流模拟。   <br>  Geophysics Journal International, vol. 191 (2012), pp.

 <li>  F Brezzi和M Fortin。   <br>  <i>Mixed and Hybrid Finite Element Methods</i>.   <br>  Springer-Verlag, 1991.

 <li>  Z陈。   <br>  <i>Finite Element Methods and Their Applications</i>.   <br>  Springer, 2005.

 <li>  JL Guermond和R Pasquetti.   <br>  基于熵的非线性粘度的守恒定律的傅里叶近似。   <br>  <i>Comptes Rendus Mathematique</i>, 346(13-14): 801-806, 2008.

 <li>  CC Chueh, M Secanell, W Bangerth, and N Djilali.   <br>  异质多孔介质中瞬态两相流的多级自适应模拟。   <br>  <i>Computers and Fluids</i>, 39:1585-1596, 2010.

 <li>  Y Saad和MH Schultz。   <br>  Gmres:用于解决非对称线性系统的广义最小残差算法。   <br>  <i>SIAM Journal on Scientific and Statistical Computing</i>, 7(3):856-869, 1986.

 <li>  F张。   <br>  <i>The Schur Complement and its Applications</i>.   <br>  Springer, 2005.

 <li>  D Silvester和A Wathen。   <br>  稳定的斯托克斯系统的快速迭代解第二部分：使用一般的块状先决条件。   <br>  <i>SIAM Journal on Numerical Analysis</i>, 31(5):1352-1367, 1994.

 <li>  GH Golub和CF van Loan。   <br>  <i>Matrix Computations</i>.   <br>  第三版，约翰霍普金斯大学，1996年。

 <li>  SE Buckley和MC Leverett。   <br>  沙子中流体位移的机制。   <br>  <i>AIME Trans.</i>, 146:107-116, 1942.

 </ol> 


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
 * 像往常一样，第一步是包括一些deal.II和C++头文件的功能。
 * 

 * 
 * 列表中包括一些提供向量、矩阵和预处理类的头文件，这些头文件实现了各自Trilinos类的接口；关于这些的一些更多信息可以在  step-31  中找到。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/tensor_function.h> 
 * #include <deal.II/base/index_set.h> 
 * 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/block_sparsity_pattern.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
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
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_block_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_vector.h> 
 * #include <deal.II/lac/trilinos_parallel_block_vector.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * #include <memory> 
 * 
 * @endcode
 * 
 * 在这个顶层设计的最后，我们为当前项目开辟一个命名空间，下面的所有材料都将进入这个命名空间，然后将所有deal.II名称导入这个命名空间。
 * 

 * 
 * 
 * @code
 * namespace Step43 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Boundaryandinitialvalueclasses"></a> 
 * <h3>Boundary and initial value classes</h3>
 * 

 * 
 * 下面的部分直接取自 step-21 ，所以没有必要重复那里的描述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class PressureBoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     PressureBoundaryValues() 
 *       : Function<dim>(1) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double 
 *   PressureBoundaryValues<dim>::value(const Point<dim> &p, 
 *                                      const unsigned int /*component*/) const 
 *   { 
 *     return 1 - p[0]; 
 *   } 
 * 
 *   template <int dim> 
 *   class SaturationBoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     SaturationBoundaryValues() 
 *       : Function<dim>(1) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double 
 *   SaturationBoundaryValues<dim>::value(const Point<dim> &p, 
 *                                        const unsigned int /*component*/) const 
 *   { 
 *     if (p[0] == 0) 
 *       return 1; 
 *     else 
 *       return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   class SaturationInitialValues : public Function<dim> 
 *   { 
 *   public: 
 *     SaturationInitialValues() 
 *       : Function<dim>(1) 
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
 *   double 
 *   SaturationInitialValues<dim>::value(const Point<dim> & /*p*/, 
 *                                       const unsigned int /*component*/) const 
 *   { 
 *     return 0.2; 
 *   } 
 * 
 *   template <int dim> 
 *   void SaturationInitialValues<dim>::vector_value(const Point<dim> &p, 
 *                                                   Vector<double> &values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = SaturationInitialValues<dim>::value(p, c); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Permeabilitymodels"></a> 
 * <h3>Permeability models</h3>
 * 

 * 
 * 在本教程中，我们仍然使用之前在 step-21 中使用的两个渗透率模型，所以我们再次避免对它们进行详细评论。
 * 

 * 
 * 
 * @code
 *   namespace SingleCurvingCrack 
 *   { 
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
 *     template <int dim> 
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
 *                                    std::vector<Tensor<2, dim>> &  values) const 
 *     { 
 *       Assert(points.size() == values.size(), 
 *              ExcDimensionMismatch(points.size(), values.size())); 
 * 
 *       for (unsigned int p = 0; p < points.size(); ++p) 
 *         { 
 *           values[p].clear(); 
 * 
 *           const double distance_to_flowline = 
 *             std::fabs(points[p][1] - 0.5 - 0.1 * std::sin(10 * points[p][0])); 
 * 
 *           const double permeability = 
 *             std::max(std::exp(-(distance_to_flowline * distance_to_flowline) / 
 *                               (0.1 * 0.1)), 
 *                      0.01); 
 * 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             values[p][d][d] = 1. / permeability; 
 *         } 
 *     } 
 *   } // namespace SingleCurvingCrack 
 * 
 *   namespace RandomMedium 
 *   { 
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
 * 
 *     private: 
 *       static std::vector<Point<dim>> centers; 
 *     }; 
 * 
 *     template <int dim> 
 *     std::vector<Point<dim>> KInverse<dim>::centers = []() { 
 *       const unsigned int N = 
 *         (dim == 2 ? 40 : (dim == 3 ? 100 : throw ExcNotImplemented())); 
 * 
 *       std::vector<Point<dim>> centers_list(N); 
 *       for (unsigned int i = 0; i < N; ++i) 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           centers_list[i][d] = static_cast<double>(rand()) / RAND_MAX; 
 * 
 *       return centers_list; 
 *     }(); 
 * 
 *     template <int dim> 
 *     void KInverse<dim>::value_list(const std::vector<Point<dim>> &points, 
 *                                    std::vector<Tensor<2, dim>> &  values) const 
 *     { 
 *       AssertDimension(points.size(), values.size()); 
 * 
 *       for (unsigned int p = 0; p < points.size(); ++p) 
 *         { 
 *           values[p].clear(); 
 * 
 *           double permeability = 0; 
 *           for (unsigned int i = 0; i < centers.size(); ++i) 
 *             permeability += 
 *               std::exp(-(points[p] - centers[i]).norm_square() / (0.05 * 0.05)); 
 * 
 *           const double normalized_permeability = 
 *             std::min(std::max(permeability, 0.01), 4.); 
 * 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             values[p][d][d] = 1. / normalized_permeability; 
 *         } 
 *     } 
 *   } // namespace RandomMedium 
 * @endcode
 * 
 * 
 * <a name="Physicalquantities"></a> 
 * <h3>Physical quantities</h3>
 * 

 * 
 * 所有物理量的实现，如总流动性 $\lambda_t$ 和水的部分流量 $F$ 都来自 step-21 ，所以我们也没有对它们做任何评论。与 step-21 相比，我们增加了检查，即传递给这些函数的饱和度实际上是在物理上有效的范围内。此外，鉴于润湿相以速度 $\mathbf u F'(S)$ 移动，很明显 $F'(S)$ 必须大于或等于零，所以我们也断言，以确保我们的计算得到的导数公式是合理的。
 * 

 * 
 * 
 * @code
 *   double mobility_inverse(const double S, const double viscosity) 
 *   { 
 *     return 1.0 / (1.0 / viscosity * S * S + (1 - S) * (1 - S)); 
 *   } 
 * 
 *   double fractional_flow(const double S, const double viscosity) 
 *   { 
 *     Assert((S >= 0) && (S <= 1), 
 *            ExcMessage("Saturation is outside its physically valid range.")); 
 * 
 *     return S * S / (S * S + viscosity * (1 - S) * (1 - S)); 
 *   } 
 * 
 *   double fractional_flow_derivative(const double S, const double viscosity) 
 *   { 
 *     Assert((S >= 0) && (S <= 1), 
 *            ExcMessage("Saturation is outside its physically valid range.")); 
 * 
 *     const double temp = (S * S + viscosity * (1 - S) * (1 - S)); 
 * 
 *     const double numerator = 
 *       2.0 * S * temp - S * S * (2.0 * S - 2.0 * viscosity * (1 - S)); 
 *     const double denominator = std::pow(temp, 2.0); 
 * 
 *     const double F_prime = numerator / denominator; 
 * 
 *     Assert(F_prime >= 0, ExcInternalError()); 
 * 
 *     return F_prime; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Helperclassesforsolversandpreconditioners"></a> 
 * <h3>Helper classes for solvers and preconditioners</h3>
 * 

 * 
 * 在这第一部分中，我们定义了一些我们在构建线性求解器和预处理器时需要的类。这一部分与  step-31  中使用的基本相同。唯一不同的是，原来的变量名称stokes_matrix被另一个名称darcy_matrix取代，以配合我们的问题。
 * 

 * 
 * 
 * @code
 *   namespace LinearSolvers 
 *   { 
 *     template <class MatrixType, class PreconditionerType> 
 *     class InverseMatrix : public Subscriptor 
 *     { 
 *     public: 
 *       InverseMatrix(const MatrixType &        m, 
 *                     const PreconditionerType &preconditioner); 
 * 
 *       template <typename VectorType> 
 *       void vmult(VectorType &dst, const VectorType &src) const; 
 * 
 *     private: 
 *       const SmartPointer<const MatrixType> matrix; 
 *       const PreconditionerType &           preconditioner; 
 *     }; 
 * 
 *     template <class MatrixType, class PreconditionerType> 
 *     InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
 *       const MatrixType &        m, 
 *       const PreconditionerType &preconditioner) 
 *       : matrix(&m) 
 *       , preconditioner(preconditioner) 
 *     {} 
 * 
 *     template <class MatrixType, class PreconditionerType> 
 *     template <typename VectorType> 
 *     void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
 *       VectorType &      dst, 
 *       const VectorType &src) const 
 *     { 
 *       SolverControl        solver_control(src.size(), 1e-7 * src.l2_norm()); 
 *       SolverCG<VectorType> cg(solver_control); 
 * 
 *       dst = 0; 
 * 
 *       try 
 *         { 
 *           cg.solve(*matrix, dst, src, preconditioner); 
 *         } 
 *       catch (std::exception &e) 
 *         { 
 *           Assert(false, ExcMessage(e.what())); 
 *         } 
 *     } 
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     class BlockSchurPreconditioner : public Subscriptor 
 *     { 
 *     public: 
 *       BlockSchurPreconditioner( 
 *         const TrilinosWrappers::BlockSparseMatrix &S, 
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                             PreconditionerTypeMp> &Mpinv, 
 *         const PreconditionerTypeA &                Apreconditioner); 
 * 
 *       void vmult(TrilinosWrappers::MPI::BlockVector &      dst, 
 *                  const TrilinosWrappers::MPI::BlockVector &src) const; 
 * 
 *     private: 
 *       const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> 
 *         darcy_matrix; 
 *       const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                                              PreconditionerTypeMp>> 
 *                                  m_inverse; 
 *       const PreconditionerTypeA &a_preconditioner; 
 * 
 *       mutable TrilinosWrappers::MPI::Vector tmp; 
 *     }; 
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>:: 
 *       BlockSchurPreconditioner( 
 *         const TrilinosWrappers::BlockSparseMatrix &S, 
 *         const InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                             PreconditionerTypeMp> &Mpinv, 
 *         const PreconditionerTypeA &                Apreconditioner) 
 *       : darcy_matrix(&S) 
 *       , m_inverse(&Mpinv) 
 *       , a_preconditioner(Apreconditioner) 
 *       , tmp(complete_index_set(darcy_matrix->block(1, 1).m())) 
 *     {} 
 * 
 *     template <class PreconditionerTypeA, class PreconditionerTypeMp> 
 *     void 
 *     BlockSchurPreconditioner<PreconditionerTypeA, PreconditionerTypeMp>::vmult( 
 *       TrilinosWrappers::MPI::BlockVector &      dst, 
 *       const TrilinosWrappers::MPI::BlockVector &src) const 
 *     { 
 *       a_preconditioner.vmult(dst.block(0), src.block(0)); 
 *       darcy_matrix->block(1, 0).residual(tmp, dst.block(0), src.block(1)); 
 *       tmp *= -1; 
 *       m_inverse->vmult(dst.block(1), tmp); 
 *     } 
 *   } // namespace LinearSolvers 
 * @endcode
 * 
 * 
 * <a name="TheTwoPhaseFlowProblemclass"></a> 
 * <h3>The TwoPhaseFlowProblem class</h3>
 * 

 * 
 * 定义解决随时间变化的平流主导的两相流问题（或Buckley-Leverett问题[Buckley 1942]）的顶层逻辑的类的定义主要基于教程程序 step-21 和 step-33 ，特别是 step-31 ，我们在这里使用的一般结构基本相同。与 step-31 一样，在下面的实现中需要寻找的关键例程是 <code>run()</code> and <code>solve()</code> 函数。
 * 

 * 
 * 与 step-31 的主要区别是，由于考虑了自适应算子拆分，我们需要多几个成员变量来保存最近两次计算的达西（速度/压力）解，以及当前的达西（直接计算，或从前两次计算中推断），我们需要记住最近两次计算的达西解。我们还需要一个辅助函数来确定我们是否真的需要重新计算达西解。
 * 

 * 
 * 与 step-31 不同，这一步多用了一个AffineConstraints对象，叫做darcy_preconditioner_constraints。这个约束对象只用于为Darcy预处理程序组装矩阵，包括悬挂节点约束以及压力变量的Dirichlet边界值约束。我们需要这个，因为我们正在为压力建立一个拉普拉斯矩阵，作为舒尔补码的近似值），如果应用边界条件，这个矩阵是正定的。
 * 

 * 
 * 这样在这个类中声明的成员函数和变量的集合与  step-31  中的相当相似。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class TwoPhaseFlowProblem 
 *   { 
 *   public: 
 *     TwoPhaseFlowProblem(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_dofs(); 
 *     void assemble_darcy_preconditioner(); 
 *     void build_darcy_preconditioner(); 
 *     void assemble_darcy_system(); 
 *     void assemble_saturation_system(); 
 *     void assemble_saturation_matrix(); 
 *     void assemble_saturation_rhs(); 
 *     void assemble_saturation_rhs_cell_term( 
 *       const FEValues<dim> &                       saturation_fe_values, 
 *       const FEValues<dim> &                       darcy_fe_values, 
 *       const double                                global_max_u_F_prime, 
 *       const double                                global_S_variation, 
 *       const std::vector<types::global_dof_index> &local_dof_indices); 
 *     void assemble_saturation_rhs_boundary_term( 
 *       const FEFaceValues<dim> &                   saturation_fe_face_values, 
 *       const FEFaceValues<dim> &                   darcy_fe_face_values, 
 *       const std::vector<types::global_dof_index> &local_dof_indices); 
 *     void solve(); 
 *     void refine_mesh(const unsigned int min_grid_level, 
 *                      const unsigned int max_grid_level); 
 *     void output_results() const; 
 * 
 * @endcode
 * 
 * 我们接下来会有一些辅助函数，这些函数在整个程序中的不同地方都会用到。
 * 

 * 
 * 
 * @code
 *     double                    get_max_u_F_prime() const; 
 *     std::pair<double, double> get_extrapolated_saturation_range() const; 
 *     bool   determine_whether_to_solve_for_pressure_and_velocity() const; 
 *     void   project_back_saturation(); 
 *     double compute_viscosity( 
 *       const std::vector<double> &        old_saturation, 
 *       const std::vector<double> &        old_old_saturation, 
 *       const std::vector<Tensor<1, dim>> &old_saturation_grads, 
 *       const std::vector<Tensor<1, dim>> &old_old_saturation_grads, 
 *       const std::vector<Vector<double>> &present_darcy_values, 
 *       const double                       global_max_u_F_prime, 
 *       const double                       global_S_variation, 
 *       const double                       cell_diameter) const; 
 * 
 * @endcode
 * 
 * 接下来是成员变量，其中大部分与 step-31 中的变量类似，但与速度/压力系统的宏观时间步长有关的变量除外。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 *     double             global_Omega_diameter; 
 * 
 *     const unsigned int degree; 
 * 
 *     const unsigned int        darcy_degree; 
 *     FESystem<dim>             darcy_fe; 
 *     DoFHandler<dim>           darcy_dof_handler; 
 *     AffineConstraints<double> darcy_constraints; 
 * 
 *     AffineConstraints<double> darcy_preconditioner_constraints; 
 * 
 *     TrilinosWrappers::BlockSparseMatrix darcy_matrix; 
 *     TrilinosWrappers::BlockSparseMatrix darcy_preconditioner_matrix; 
 * 
 *     TrilinosWrappers::MPI::BlockVector darcy_solution; 
 *     TrilinosWrappers::MPI::BlockVector darcy_rhs; 
 * 
 *     TrilinosWrappers::MPI::BlockVector last_computed_darcy_solution; 
 *     TrilinosWrappers::MPI::BlockVector second_last_computed_darcy_solution; 
 * 
 *     const unsigned int        saturation_degree; 
 *     FE_Q<dim>                 saturation_fe; 
 *     DoFHandler<dim>           saturation_dof_handler; 
 *     AffineConstraints<double> saturation_constraints; 
 * 
 *     TrilinosWrappers::SparseMatrix saturation_matrix; 
 * 
 *     TrilinosWrappers::MPI::Vector saturation_solution; 
 *     TrilinosWrappers::MPI::Vector old_saturation_solution; 
 *     TrilinosWrappers::MPI::Vector old_old_saturation_solution; 
 *     TrilinosWrappers::MPI::Vector saturation_rhs; 
 * 
 *     TrilinosWrappers::MPI::Vector 
 *       saturation_matching_last_computed_darcy_solution; 
 * 
 *     const double saturation_refinement_threshold; 
 * 
 *     double       time; 
 *     const double end_time; 
 * 
 *     double current_macro_time_step; 
 *     double old_macro_time_step; 
 * 
 *     double       time_step; 
 *     double       old_time_step; 
 *     unsigned int timestep_number; 
 * 
 *     const double viscosity; 
 *     const double porosity; 
 *     const double AOS_threshold; 
 * 
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC> Amg_preconditioner; 
 *     std::shared_ptr<TrilinosWrappers::PreconditionIC> Mp_preconditioner; 
 * 
 *     bool rebuild_saturation_matrix; 
 * 
 * @endcode
 * 
 * 在最后，我们声明一个变量，表示材料模型。与 step-21 相比，我们在这里把它作为一个成员变量，因为我们想在不同的地方使用它，所以有一个声明这样一个变量的中心位置，将使我们更容易用另一个类来替换 RandomMedium::KInverse （例如，用 SingleCurvingCrack::KInverse). 替换 RandomMedium::KInverse ）。
 * 
 * @code
 *     const RandomMedium::KInverse<dim> k_inverse; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimTwoPhaseFlowProblem"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem</h3>
 * 

 * 
 * 这个类的构造函数是对  step-21  和  step-31  中的构造函数的扩展。我们需要添加涉及饱和度的各种变量。正如介绍中所讨论的，我们将再次使用 $Q_2 \times Q_1$ （Taylor-Hood）元素来处理Darcy系统，这是一个满足Ladyzhenskaya-Babuska-Brezzi（LBB）条件的元素组合[Brezzi and Fortin 1991, Chen 2005]，并使用 $Q_1$ 元素处理饱和度。然而，通过使用存储Darcy和温度有限元的多项式程度的变量，可以很容易地持续修改这些元素的程度以及在其上使用的所有正交公式的下游。此外，我们还初始化了与算子分割有关的时间步进变量，以及矩阵装配和预处理的选项。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem(const unsigned int degree) 
 *     : triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , global_Omega_diameter(std::numeric_limits<double>::quiet_NaN()) 
 *     , degree(degree) 
 *     , darcy_degree(degree) 
 *     , darcy_fe(FE_Q<dim>(darcy_degree + 1), dim, FE_Q<dim>(darcy_degree), 1) 
 *     , darcy_dof_handler(triangulation) 
 *     , 
 * 
 *     saturation_degree(degree + 1) 
 *     , saturation_fe(saturation_degree) 
 *     , saturation_dof_handler(triangulation) 
 *     , 
 * 
 *     saturation_refinement_threshold(0.5) 
 *     , 
 * 
 *     time(0) 
 *     , end_time(10) 
 *     , 
 * 
 *     current_macro_time_step(0) 
 *     , old_macro_time_step(0) 
 *     , 
 * 
 *     time_step(0) 
 *     , old_time_step(0) 
 *     , timestep_number(0) 
 *     , viscosity(0.2) 
 *     , porosity(1.0) 
 *     , AOS_threshold(3.0) 
 *     , 
 * 
 *     rebuild_saturation_matrix(true) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimsetup_dofs"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::setup_dofs</h3>
 * 

 * 
 * 这个函数设置了我们这里的DoFHandler对象（一个用于Darcy部分，一个用于饱和部分），以及将本程序中线性代数所需的各种对象设置为合适的尺寸。其基本操作与 step-31 所做的类似。
 * 

 * 
 * 该函数的主体首先列举了达西和饱和系统的所有自由度。对于Darcy部分，自由度会被排序，以确保速度优先于压力DoF，这样我们就可以将Darcy矩阵划分为一个 $2 \times 2$ 矩阵。
 * 然后，
 * 我们需要将悬挂节点约束和Dirichlet边界值约束纳入 darcy_preconditioner_constraints。 边界条件约束只设置在压力分量上，因为对应于非混合形式的多孔介质流算子的Schur complement预处理程序 $-\nabla \cdot [\mathbf K \lambda_t(S)]\nabla$  ，只作用于压力变量。因此，我们使用一个过滤掉速度分量的分量掩码，这样就可以只对压力自由度进行缩减。
 * 

 * 
 * 做完这些后，我们计算各个块中的自由度数量。然后，这些信息被用来创建达西和饱和系统矩阵的稀疏模式，以及用于建立达西预处理的预处理矩阵。如同 step-31 ，我们选择使用DynamicSparsityPattern的封锁版本来创建模式。因此，对于这一点，我们遵循与 step-31 相同的方式，对于成员函数的其他部分，我们不必再重复描述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::setup_dofs() 
 *   { 
 *     std::vector<unsigned int> darcy_block_component(dim + 1, 0); 
 *     darcy_block_component[dim] = 1; 
 *     { 
 *       darcy_dof_handler.distribute_dofs(darcy_fe); 
 *       DoFRenumbering::Cuthill_McKee(darcy_dof_handler); 
 *       DoFRenumbering::component_wise(darcy_dof_handler, darcy_block_component); 
 * 
 *       darcy_constraints.clear(); 
 *       DoFTools::make_hanging_node_constraints(darcy_dof_handler, 
 *                                               darcy_constraints); 
 *       darcy_constraints.close(); 
 *     } 
 *     { 
 *       saturation_dof_handler.distribute_dofs(saturation_fe); 
 * 
 *       saturation_constraints.clear(); 
 *       DoFTools::make_hanging_node_constraints(saturation_dof_handler, 
 *                                               saturation_constraints); 
 *       saturation_constraints.close(); 
 *     } 
 *     { 
 *       darcy_preconditioner_constraints.clear(); 
 * 
 *       FEValuesExtractors::Scalar pressure(dim); 
 * 
 *       DoFTools::make_hanging_node_constraints(darcy_dof_handler, 
 *                                               darcy_preconditioner_constraints); 
 *       DoFTools::make_zero_boundary_constraints(darcy_dof_handler, 
 *                                                darcy_preconditioner_constraints, 
 *                                                darcy_fe.component_mask( 
 *                                                  pressure)); 
 * 
 *       darcy_preconditioner_constraints.close(); 
 *     } 
 * 
 *     const std::vector<types::global_dof_index> darcy_dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(darcy_dof_handler, 
 *                                         darcy_block_component); 
 *     const unsigned int n_u = darcy_dofs_per_block[0], 
 *                        n_p = darcy_dofs_per_block[1], 
 *                        n_s = saturation_dof_handler.n_dofs(); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << " (on " << triangulation.n_levels() << " levels)" << std::endl 
 *               << "Number of degrees of freedom: " << n_u + n_p + n_s << " (" 
 *               << n_u << '+' << n_p << '+' << n_s << ')' << std::endl 
 *               << std::endl; 
 * 
 *     { 
 *       darcy_matrix.clear(); 
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2); 
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u); 
 *       dsp.block(0, 1).reinit(n_u, n_p); 
 *       dsp.block(1, 0).reinit(n_p, n_u); 
 *       dsp.block(1, 1).reinit(n_p, n_p); 
 * 
 *       dsp.collect_sizes(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 * 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (!((c == dim) && (d == dim))) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         darcy_dof_handler, coupling, dsp, darcy_constraints, false); 
 * 
 *       darcy_matrix.reinit(dsp); 
 *     } 
 * 
 *     { 
 *       Amg_preconditioner.reset(); 
 *       Mp_preconditioner.reset(); 
 *       darcy_preconditioner_matrix.clear(); 
 * 
 *       BlockDynamicSparsityPattern dsp(2, 2); 
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u); 
 *       dsp.block(0, 1).reinit(n_u, n_p); 
 *       dsp.block(1, 0).reinit(n_p, n_u); 
 *       dsp.block(1, 1).reinit(n_p, n_p); 
 * 
 *       dsp.collect_sizes(); 
 * 
 *       Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (c == d) 
 *             coupling[c][d] = DoFTools::always; 
 *           else 
 *             coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern( 
 *         darcy_dof_handler, coupling, dsp, darcy_constraints, false); 
 * 
 *       darcy_preconditioner_matrix.reinit(dsp); 
 *     } 
 * 
 *     { 
 *       saturation_matrix.clear(); 
 * 
 *       DynamicSparsityPattern dsp(n_s, n_s); 
 * 
 *       DoFTools::make_sparsity_pattern(saturation_dof_handler, 
 *                                       dsp, 
 *                                       saturation_constraints, 
 *                                       false); 
 * 
 *       saturation_matrix.reinit(dsp); 
 *     } 
 * 
 *     std::vector<IndexSet> darcy_partitioning(2); 
 *     darcy_partitioning[0] = complete_index_set(n_u); 
 *     darcy_partitioning[1] = complete_index_set(n_p); 
 *     darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD); 
 *     darcy_solution.collect_sizes(); 
 * 
 *     last_computed_darcy_solution.reinit(darcy_partitioning, MPI_COMM_WORLD); 
 *     last_computed_darcy_solution.collect_sizes(); 
 * 
 *     second_last_computed_darcy_solution.reinit(darcy_partitioning, 
 *                                                MPI_COMM_WORLD); 
 *     second_last_computed_darcy_solution.collect_sizes(); 
 * 
 *     darcy_rhs.reinit(darcy_partitioning, MPI_COMM_WORLD); 
 *     darcy_rhs.collect_sizes(); 
 * 
 *     IndexSet saturation_partitioning = complete_index_set(n_s); 
 *     saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 
 *     old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 
 *     old_old_saturation_solution.reinit(saturation_partitioning, MPI_COMM_WORLD); 
 * 
 *     saturation_matching_last_computed_darcy_solution.reinit( 
 *       saturation_partitioning, MPI_COMM_WORLD); 
 * 
 *     saturation_rhs.reinit(saturation_partitioning, MPI_COMM_WORLD); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Assemblingmatricesandpreconditioners"></a> 
 * <h3>Assembling matrices and preconditioners</h3>
 * 

 * 
 * 接下来的几个函数专门用来设置我们在这个程序中必须处理的各种系统和预处理矩阵及右手边。
 * 

 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_darcy_preconditioner"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner</h4>
 * 

 * 
 * 这个函数组装我们用于预处理达西系统的矩阵。我们需要的是在速度分量上用 $\left(\mathbf{K} \lambda_t\right)^{-1}$ 加权的向量质量矩阵和在压力分量上用 $\left(\mathbf{K} \lambda_t\right)$ 加权的质量矩阵。我们首先生成一个适当阶数的正交对象，即FEValues对象，可以给出正交点的数值和梯度（连同正交权重）。接下来我们为单元格矩阵和局部与全局DoF之间的关系创建数据结构。向量phi_u和grad_phi_p将保存基函数的值，以便更快地建立局部矩阵，正如在  step-22  中已经做的。在我们开始对所有活动单元进行循环之前，我们必须指定哪些成分是压力，哪些是速度。
 * 

 * 
 * 局部矩阵的创建是相当简单的。只有一个由 $\left(\mathbf{K} \lambda_t\right)^{-1}$ 加权的项（关于速度）和一个由 $\left(\mathbf{K} \lambda_t\right)$ 加权的拉普拉斯矩阵需要生成，所以局部矩阵的创建基本上只需要两行就可以完成。由于该文件顶部的材料模型函数只提供了渗透率和迁移率的倒数，我们必须根据给定的数值手工计算 $\mathbf K$ 和 $\lambda_t$ ，每个正交点一次。
 * 

 * 
 * 一旦本地矩阵准备好了（在每个正交点上对本地矩阵的行和列进行循环），我们就可以得到本地的DoF指数，并将本地信息写入全局矩阵中。我们通过直接应用约束条件（即darcy_preconditioner_constraints）来做到这一点，该约束条件负责处理悬挂节点和零Dirichlet边界条件约束。这样做，我们就不必事后再做，以后也不必使用 AffineConstraints::condense 和 MatrixTools::apply_boundary_values, 这两个需要修改矩阵和向量项的函数，因此对于我们不能立即访问单个内存位置的特里诺斯类来说，很难编写。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner() 
 *   { 
 *     std::cout << "   Rebuilding darcy preconditioner..." << std::endl; 
 * 
 *     darcy_preconditioner_matrix = 0; 
 * 
 *     const QGauss<dim> quadrature_formula(darcy_degree + 2); 
 *     FEValues<dim>     darcy_fe_values(darcy_fe, 
 *                                   quadrature_formula, 
 *                                   update_JxW_values | update_values | 
 *                                     update_gradients | 
 *                                     update_quadrature_points); 
 *     FEValues<dim>     saturation_fe_values(saturation_fe, 
 *                                        quadrature_formula, 
 *                                        update_values); 
 * 
 *     const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 
 * 
 *     std::vector<double> old_saturation_values(n_q_points); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> grad_phi_p(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 *     auto       cell            = darcy_dof_handler.begin_active(); 
 *     const auto endc            = darcy_dof_handler.end(); 
 *     auto       saturation_cell = saturation_dof_handler.begin_active(); 
 * 
 *     for (; cell != endc; ++cell, ++saturation_cell) 
 *       { 
 *         darcy_fe_values.reinit(cell); 
 *         saturation_fe_values.reinit(saturation_cell); 
 * 
 *         local_matrix = 0; 
 * 
 *         saturation_fe_values.get_function_values(old_saturation_solution, 
 *                                                  old_saturation_values); 
 * 
 *         k_inverse.value_list(darcy_fe_values.get_quadrature_points(), 
 *                              k_inverse_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             const double old_s = old_saturation_values[q]; 
 * 
 *             const double inverse_mobility = mobility_inverse(old_s, viscosity); 
 *             const double mobility         = 1.0 / inverse_mobility; 
 *             const Tensor<2, dim> permeability = invert(k_inverse_values[q]); 
 * 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 phi_u[k]      = darcy_fe_values[velocities].value(k, q); 
 *                 grad_phi_p[k] = darcy_fe_values[pressure].gradient(k, q); 
 *               } 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   local_matrix(i, j) += 
 *                     (k_inverse_values[q] * inverse_mobility * phi_u[i] * 
 *                        phi_u[j] + 
 *                      permeability * mobility * grad_phi_p[i] * grad_phi_p[j]) * 
 *                     darcy_fe_values.JxW(q); 
 *                 } 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         darcy_preconditioner_constraints.distribute_local_to_global( 
 *           local_matrix, local_dof_indices, darcy_preconditioner_matrix); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimbuild_darcy_preconditioner"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::build_darcy_preconditioner</h4>
 * 

 * 
 * 在调用上述函数组装预处理矩阵后，该函数生成将用于舒尔补块预处理的内部预处理器。前置条件需要在每个饱和时间步长时重新生成，因为它们取决于随时间变化的饱和度  $S$  。
 * 

 * 
 * 在这里，我们为速度-速度矩阵  $\mathbf{M}^{\mathbf{u}}$  和Schur补码  $\mathbf{S}$  设置了预处理器。正如介绍中所解释的，我们将使用一个基于矢量矩阵 $\mathbf{M}^{\mathbf{u}}$ 的IC预处理器和另一个基于标量拉普拉斯矩阵 $\tilde{\mathbf{S}}^p$ 的IC预处理器（它在频谱上与达西矩阵的舒尔补码接近）。通常， TrilinosWrappers::PreconditionIC 类可以被看作是一个很好的黑盒预处理程序，不需要对矩阵结构和/或背后的算子有任何特殊的了解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::build_darcy_preconditioner() 
 *   { 
 *     assemble_darcy_preconditioner(); 
 * 
 *     Amg_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
 *     Amg_preconditioner->initialize(darcy_preconditioner_matrix.block(0, 0)); 
 * 
 *     Mp_preconditioner = std::make_shared<TrilinosWrappers::PreconditionIC>(); 
 *     Mp_preconditioner->initialize(darcy_preconditioner_matrix.block(1, 1)); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_darcy_system"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_darcy_system</h4>
 * 

 * 
 * 这是为达西系统组装线性系统的函数。
 * 

 * 
 * 关于执行的技术细节，其程序与  step-22  和  step-31  中的程序相似。我们重置矩阵和向量，在单元格上创建正交公式，然后创建相应的FEValues对象。
 * 

 * 
 * 有一件事需要评论：由于我们有一个单独的有限元和DoFHandler来处理饱和问题，我们需要生成第二个FEValues对象来正确评估饱和解。要实现这一点并不复杂：只需使用饱和结构，并为基函数值设置一个更新标志，我们需要对饱和解进行评估。这里需要记住的唯一重要部分是，两个FEValues对象使用相同的正交公式，以确保我们在循环计算两个对象的正交点时获得匹配的信息。
 * 

 * 
 * 声明的过程中，对数组的大小、本地矩阵的创建、右手边以及与全局系统相比较的本地道夫指数的向量都有一些快捷方式。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_darcy_system() 
 *   { 
 *     darcy_matrix = 0; 
 *     darcy_rhs    = 0; 
 * 
 *     QGauss<dim>     quadrature_formula(darcy_degree + 2); 
 *     QGauss<dim - 1> face_quadrature_formula(darcy_degree + 2); 
 * 
 *     FEValues<dim> darcy_fe_values(darcy_fe, 
 *                                   quadrature_formula, 
 *                                   update_values | update_gradients | 
 *                                     update_quadrature_points | 
 *                                     update_JxW_values); 
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe, 
 *                                        quadrature_formula, 
 *                                        update_values); 
 * 
 *     FEFaceValues<dim> darcy_fe_face_values(darcy_fe, 
 *                                            face_quadrature_formula, 
 *                                            update_values | 
 *                                              update_normal_vectors | 
 *                                              update_quadrature_points | 
 *                                              update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = darcy_fe.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points      = quadrature_formula.size(); 
 *     const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const Functions::ZeroFunction<dim> pressure_right_hand_side; 
 *     const PressureBoundaryValues<dim>  pressure_boundary_values; 
 * 
 *     std::vector<double>         pressure_rhs_values(n_q_points); 
 *     std::vector<double>         boundary_values(n_face_q_points); 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 
 * 
 * @endcode
 * 
 * 接下来我们需要一个向量，该向量将包含前一时间层在正交点的饱和解的值，以组装达西方程中的饱和相关系数。
 * 

 * 
 * 我们接下来创建的向量集包含了基函数的评价以及它们的梯度，将用于创建矩阵。把这些放到自己的数组中，而不是每次都向FEValues对象索取这些信息，是为了加速装配过程的优化，详情请见 step-22 。
 * 

 * 
 * 最后两个声明是用来从整个FE系统中提取各个块（速度、压力、饱和度）的。
 * 

 * 
 * 
 * @code
 *     std::vector<double> old_saturation_values(n_q_points); 
 * 
 *     std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
 *     std::vector<double>         div_phi_u(dofs_per_cell); 
 *     std::vector<double>         phi_p(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 * @endcode
 * 
 * 现在开始对问题中的所有单元格进行循环。我们在这个装配例程中使用了两个不同的DoFHandlers，所以我们必须为使用中的两个对象设置两个不同的单元格迭代器。这可能看起来有点奇怪，但是由于达西系统和饱和系统都使用相同的网格，我们可以假设这两个迭代器在两个DoFHandler对象的单元格中同步运行。
 * 

 * 
 * 循环中的第一条语句又是非常熟悉的，按照更新标志的规定对有限元数据进行更新，将局部数组清零，并得到正交点上的旧解的值。 在这一点上，我们还必须在正交点上获得前一个时间步长的饱和函数的值。为此，我们可以使用 FEValues::get_function_values （之前已经在 step-9 、 step-14 和 step-15 中使用），这个函数接收一个解向量，并返回当前单元的正交点的函数值列表。事实上，它返回每个正交点的完整矢量值解，即不仅是饱和度，还有速度和压力。
 * 

 * 
 * 然后，我们就可以在单元格上的正交点上进行循环，以进行积分。这方面的公式直接来自介绍中所讨论的内容。
 * 

 * 
 * 一旦这样做了，我们就开始在局部矩阵的行和列上进行循环，并将相关的乘积输入矩阵中。
 * 

 * 
 * 循环所有单元的最后一步是将本地贡献输入到全局矩阵和向量结构中，并在local_dof_indices中指定位置。同样，我们让AffineConstraints类将单元格矩阵元素插入到全局矩阵中，全局矩阵已经浓缩了悬挂节点的约束。
 * 

 * 
 * 
 * @code
 *     auto       cell            = darcy_dof_handler.begin_active(); 
 *     const auto endc            = darcy_dof_handler.end(); 
 *     auto       saturation_cell = saturation_dof_handler.begin_active(); 
 * 
 *     for (; cell != endc; ++cell, ++saturation_cell) 
 *       { 
 *         darcy_fe_values.reinit(cell); 
 *         saturation_fe_values.reinit(saturation_cell); 
 * 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 *         saturation_fe_values.get_function_values(old_saturation_solution, 
 *                                                  old_saturation_values); 
 * 
 *         pressure_right_hand_side.value_list( 
 *           darcy_fe_values.get_quadrature_points(), pressure_rhs_values); 
 *         k_inverse.value_list(darcy_fe_values.get_quadrature_points(), 
 *                              k_inverse_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 phi_u[k]     = darcy_fe_values[velocities].value(k, q); 
 *                 div_phi_u[k] = darcy_fe_values[velocities].divergence(k, q); 
 *                 phi_p[k]     = darcy_fe_values[pressure].value(k, q); 
 *               } 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               { 
 *                 const double old_s = old_saturation_values[q]; 
 *                 for (unsigned int j = 0; j <= i; ++j) 
 *                   { 
 *                     local_matrix(i, j) += 
 *                       (phi_u[i] * k_inverse_values[q] * 
 *                          mobility_inverse(old_s, viscosity) * phi_u[j] - 
 *                        div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * 
 *                       darcy_fe_values.JxW(q); 
 *                   } 
 * 
 *                 local_rhs(i) += 
 *                   (-phi_p[i] * pressure_rhs_values[q]) * darcy_fe_values.JxW(q); 
 *               } 
 *           } 
 * 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary()) 
 *             { 
 *               darcy_fe_face_values.reinit(cell, face); 
 * 
 *               pressure_boundary_values.value_list( 
 *                 darcy_fe_face_values.get_quadrature_points(), boundary_values); 
 * 
 *               for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *                 for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *                   { 
 *                     const Tensor<1, dim> phi_i_u = 
 *                       darcy_fe_face_values[velocities].value(i, q); 
 * 
 *                     local_rhs(i) += 
 *                       -(phi_i_u * darcy_fe_face_values.normal_vector(q) * 
 *                         boundary_values[q] * darcy_fe_face_values.JxW(q)); 
 *                   } 
 *             } 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
 *             local_matrix(i, j) = local_matrix(j, i); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         darcy_constraints.distribute_local_to_global( 
 *           local_matrix, local_rhs, local_dof_indices, darcy_matrix, darcy_rhs); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_system"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_system</h4>
 * 

 * 
 * 这个函数是为了组装饱和传输方程的线性系统。如果有必要，它会调用另外两个成员函数：assemble_saturation_matrix()和assemble_saturation_rhs()。前一个函数然后组装饱和度矩阵，只需要偶尔改变。另一方面，后一个组装右手边的函数必须在每个饱和时间步骤中调用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_system() 
 *   { 
 *     if (rebuild_saturation_matrix == true) 
 *       { 
 *         saturation_matrix = 0; 
 *         assemble_saturation_matrix(); 
 *       } 
 * 
 *     saturation_rhs = 0; 
 *     assemble_saturation_rhs(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_matrix"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_matrix</h4>
 * 

 * 
 * 这个函数很容易理解，因为它只是通过基函数phi_i_s和phi_j_s为饱和线性系统的左侧形成一个简单的质量矩阵。最后，像往常一样，我们通过在local_dof_indices中指定位置将局部贡献输入全局矩阵。这是通过让AffineConstraints类将单元矩阵元素插入全局矩阵来完成的，全局矩阵已经浓缩了悬挂节点约束。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix() 
 *   { 
 *     QGauss<dim> quadrature_formula(saturation_degree + 2); 
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe, 
 *                                        quadrature_formula, 
 *                                        update_values | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = saturation_fe.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *       { 
 *         saturation_fe_values.reinit(cell); 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               const double phi_i_s = saturation_fe_values.shape_value(i, q); 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 { 
 *                   const double phi_j_s = saturation_fe_values.shape_value(j, q); 
 *                   local_matrix(i, j) += 
 *                     porosity * phi_i_s * phi_j_s * saturation_fe_values.JxW(q); 
 *                 } 
 *             } 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         saturation_constraints.distribute_local_to_global(local_matrix, 
 *                                                           local_dof_indices, 
 *                                                           saturation_matrix); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs</h4>
 * 

 * 
 * 这个函数是用来组装饱和传输方程的右边。在进行这项工作之前，我们必须为达西系统和饱和系统分别创建两个FEValues对象，此外，还必须为这两个系统创建两个FEFaceValues对象，因为我们在饱和方程的弱形式中存在一个边界积分项。对于饱和系统的FEFaceValues对象，我们还需要法向量，我们使用update_normal_vectors标志来申请。
 * 

 * 
 * 接下来，在对所有单元进行循环之前，我们必须计算一些参数（例如global_u_infty、global_S_variation和global_Omega_diameter），这是人工黏度 $\nu$ 需要的。这与 step-31 中的做法基本相同，所以你可以在那里看到更多的信息。
 * 

 * 
 * 真正的工作是从循环所有的饱和和Darcy单元开始的，以便将局部贡献放到全局矢量中。在这个循环中，为了简化实现，我们把一些工作分成两个辅助函数：assemble_saturation_rhs_cell_term和assemble_saturation_rhs_boundary_term。 我们注意到，我们在这两个函数中把细胞或边界贡献插入全局向量，而不是在本函数中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs() 
 *   { 
 *     QGauss<dim>     quadrature_formula(saturation_degree + 2); 
 *     QGauss<dim - 1> face_quadrature_formula(saturation_degree + 2); 
 * 
 *     FEValues<dim> saturation_fe_values(saturation_fe, 
 *                                        quadrature_formula, 
 *                                        update_values | update_gradients | 
 *                                          update_quadrature_points | 
 *                                          update_JxW_values); 
 *     FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values); 
 *     FEFaceValues<dim> saturation_fe_face_values(saturation_fe, 
 *                                                 face_quadrature_formula, 
 *                                                 update_values | 
 *                                                   update_normal_vectors | 
 *                                                   update_quadrature_points | 
 *                                                   update_JxW_values); 
 *     FEFaceValues<dim> darcy_fe_face_values(darcy_fe, 
 *                                            face_quadrature_formula, 
 *                                            update_values); 
 *     FEFaceValues<dim> saturation_fe_face_values_neighbor( 
 *       saturation_fe, face_quadrature_formula, update_values); 
 * 
 *     const unsigned int dofs_per_cell = 
 *       saturation_dof_handler.get_fe().n_dofs_per_cell(); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const double                    global_max_u_F_prime = get_max_u_F_prime(); 
 *     const std::pair<double, double> global_S_range = 
 *       get_extrapolated_saturation_range(); 
 *     const double global_S_variation = 
 *       global_S_range.second - global_S_range.first; 
 * 
 *     auto       cell       = saturation_dof_handler.begin_active(); 
 *     const auto endc       = saturation_dof_handler.end(); 
 *     auto       darcy_cell = darcy_dof_handler.begin_active(); 
 *     for (; cell != endc; ++cell, ++darcy_cell) 
 *       { 
 *         saturation_fe_values.reinit(cell); 
 *         darcy_fe_values.reinit(darcy_cell); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         assemble_saturation_rhs_cell_term(saturation_fe_values, 
 *                                           darcy_fe_values, 
 *                                           global_max_u_F_prime, 
 *                                           global_S_variation, 
 *                                           local_dof_indices); 
 * 
 *         for (const auto &face : cell->face_iterators()) 
 *           if (face->at_boundary()) 
 *             { 
 *               darcy_fe_face_values.reinit(darcy_cell, face); 
 *               saturation_fe_face_values.reinit(cell, face); 
 *               assemble_saturation_rhs_boundary_term(saturation_fe_face_values, 
 *                                                     darcy_fe_face_values, 
 *                                                     local_dof_indices); 
 *             } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs_cell_term"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term</h4>
 * 

 * 
 * 这个函数负责整合饱和度方程右边的单元项，然后将其组装成全局右边的矢量。鉴于介绍中的讨论，这些贡献的形式很清楚。唯一棘手的部分是获得人工黏度和计算它所需的一切。该函数的前半部分专门用于这项任务。
 * 

 * 
 * 该函数的最后一部分是将局部贡献复制到全局向量中，其位置由local_dof_indices指定。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term( 
 *     const FEValues<dim> &                       saturation_fe_values, 
 *     const FEValues<dim> &                       darcy_fe_values, 
 *     const double                                global_max_u_F_prime, 
 *     const double                                global_S_variation, 
 *     const std::vector<types::global_dof_index> &local_dof_indices) 
 *   { 
 *     const unsigned int dofs_per_cell = saturation_fe_values.dofs_per_cell; 
 *     const unsigned int n_q_points    = saturation_fe_values.n_quadrature_points; 
 * 
 *     std::vector<double>         old_saturation_solution_values(n_q_points); 
 *     std::vector<double>         old_old_saturation_solution_values(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_grad_saturation_solution_values(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_old_grad_saturation_solution_values( 
 *       n_q_points); 
 *     std::vector<Vector<double>> present_darcy_solution_values( 
 *       n_q_points, Vector<double>(dim + 1)); 
 * 
 *     saturation_fe_values.get_function_values(old_saturation_solution, 
 *                                              old_saturation_solution_values); 
 *     saturation_fe_values.get_function_values( 
 *       old_old_saturation_solution, old_old_saturation_solution_values); 
 *     saturation_fe_values.get_function_gradients( 
 *       old_saturation_solution, old_grad_saturation_solution_values); 
 *     saturation_fe_values.get_function_gradients( 
 *       old_old_saturation_solution, old_old_grad_saturation_solution_values); 
 *     darcy_fe_values.get_function_values(darcy_solution, 
 *                                         present_darcy_solution_values); 
 * 
 *     const double nu = 
 *       compute_viscosity(old_saturation_solution_values, 
 *                         old_old_saturation_solution_values, 
 *                         old_grad_saturation_solution_values, 
 *                         old_old_grad_saturation_solution_values, 
 *                         present_darcy_solution_values, 
 *                         global_max_u_F_prime, 
 *                         global_S_variation, 
 *                         saturation_fe_values.get_cell()->diameter()); 
 * 
 *     Vector<double> local_rhs(dofs_per_cell); 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *         { 
 *           const double   old_s = old_saturation_solution_values[q]; 
 *           Tensor<1, dim> present_u; 
 *           for (unsigned int d = 0; d < dim; ++d) 
 *             present_u[d] = present_darcy_solution_values[q](d); 
 * 
 *           const double         phi_i_s = saturation_fe_values.shape_value(i, q); 
 *           const Tensor<1, dim> grad_phi_i_s = 
 *             saturation_fe_values.shape_grad(i, q); 
 * 
 *           local_rhs(i) += 
 *             (time_step * fractional_flow(old_s, viscosity) * present_u * 
 *                grad_phi_i_s - 
 *              time_step * nu * old_grad_saturation_solution_values[q] * 
 *                grad_phi_i_s + 
 *              porosity * old_s * phi_i_s) * 
 *             saturation_fe_values.JxW(q); 
 *         } 
 * 
 *     saturation_constraints.distribute_local_to_global(local_rhs, 
 *                                                       local_dof_indices, 
 *                                                       saturation_rhs); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimassemble_saturation_rhs_boundary_term"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term</h4>
 * 

 * 
 * 下一个函数负责饱和方程右侧形式中的边界积分项。 对于这些，我们必须计算全局边界面上的上行通量，也就是说，我们只对全局边界的流入部分弱加迪里切特边界条件。如前所述，这在 step-21 中已经描述过了，所以我们不对其进行更多的描述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term( 
 *     const FEFaceValues<dim> &                   saturation_fe_face_values, 
 *     const FEFaceValues<dim> &                   darcy_fe_face_values, 
 *     const std::vector<types::global_dof_index> &local_dof_indices) 
 *   { 
 *     const unsigned int dofs_per_cell = saturation_fe_face_values.dofs_per_cell; 
 *     const unsigned int n_face_q_points = 
 *       saturation_fe_face_values.n_quadrature_points; 
 * 
 *     Vector<double> local_rhs(dofs_per_cell); 
 * 
 *  
 *     std::vector<Vector<double>> present_darcy_solution_values_face( 
 *       n_face_q_points, Vector<double>(dim + 1)); 
 *     std::vector<double> neighbor_saturation(n_face_q_points); 
 * 
 *     saturation_fe_face_values.get_function_values( 
 *       old_saturation_solution, old_saturation_solution_values_face); 
 *     darcy_fe_face_values.get_function_values( 
 *       darcy_solution, present_darcy_solution_values_face); 
 * 
 *     SaturationBoundaryValues<dim> saturation_boundary_values; 
 *     saturation_boundary_values.value_list( 
 *       saturation_fe_face_values.get_quadrature_points(), neighbor_saturation); 
 * 
 *     for (unsigned int q = 0; q < n_face_q_points; ++q) 
 *       { 
 *         Tensor<1, dim> present_u_face; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           present_u_face[d] = present_darcy_solution_values_face[q](d); 
 * 
 *  
 *           present_u_face * saturation_fe_face_values.normal_vector(q); 
 * 
 *         const bool is_outflow_q_point = (normal_flux >= 0); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           local_rhs(i) -= 
 *             time_step * normal_flux * 
 *             fractional_flow((is_outflow_q_point == true ? 
 *                                old_saturation_solution_values_face[q] : 
 *                                neighbor_saturation[q]), 
 *                             viscosity) * 
 *             saturation_fe_face_values.shape_value(i, q) * 
 *             saturation_fe_face_values.JxW(q); 
 *       } 
 *     saturation_constraints.distribute_local_to_global(local_rhs, 
 *                                                       local_dof_indices, 
 *                                                       saturation_rhs); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimsolve"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::solve</h3>
 * 

 * 
 * 该函数实现了算子分割算法，即在每个时间步长中，它要么重新计算达西系统的解，要么从以前的时间步长中推算出速度/压力，然后确定时间步长的大小，然后更新饱和度变量。其实现主要遵循  step-31  中的类似代码。除了run()函数外，它是本程序中的核心函数。
 * 

 * 
 * 在函数的开始，我们询问是否要通过评估后验准则来解决压力-速度部分（见下面的函数）。如果有必要，我们将使用GMRES求解器和Schur补充块预处理来求解压力-速度部分，如介绍中所述。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::solve() 
 *   { 
 *     const bool solve_for_pressure_and_velocity = 
 *       determine_whether_to_solve_for_pressure_and_velocity(); 
 * 
 *     if (solve_for_pressure_and_velocity == true) 
 *       { 
 *         std::cout << "   Solving Darcy (pressure-velocity) system..." 
 *                   << std::endl; 
 * 
 *         assemble_darcy_system(); 
 *         build_darcy_preconditioner(); 
 * 
 *         { 
 *           const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix, 
 *                                              TrilinosWrappers::PreconditionIC> 
 *             mp_inverse(darcy_preconditioner_matrix.block(1, 1), 
 *                        *Mp_preconditioner); 
 * 
 *           const LinearSolvers::BlockSchurPreconditioner< 
 *             TrilinosWrappers::PreconditionIC, 
 *             TrilinosWrappers::PreconditionIC> 
 *             preconditioner(darcy_matrix, mp_inverse, *Amg_preconditioner); 
 * 
 *           SolverControl solver_control(darcy_matrix.m(), 
 *                                        1e-16 * darcy_rhs.l2_norm()); 
 * 
 *           SolverGMRES<TrilinosWrappers::MPI::BlockVector> gmres( 
 *             solver_control, 
 *             SolverGMRES<TrilinosWrappers::MPI::BlockVector>::AdditionalData( 
 *               100)); 
 * 
 *           for (unsigned int i = 0; i < darcy_solution.size(); ++i) 
 *             if (darcy_constraints.is_constrained(i)) 
 *               darcy_solution(i) = 0; 
 * 
 *           gmres.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner); 
 * 
 *           darcy_constraints.distribute(darcy_solution); 
 * 
 *           std::cout << "        ..." << solver_control.last_step() 
 *                     << " GMRES iterations." << std::endl; 
 *         } 
 * 
 *         { 
 *   }; 
 *           last_computed_darcy_solution        = darcy_solution; 
 * 
 *           saturation_matching_last_computed_darcy_solution = 
 *             saturation_solution; 
 *         } 
 *       } 
 * 
 * @endcode
 * 
 * 另一方面，如果我们决定不计算当前时间步长的达西系统的解，那么我们需要简单地将前两个达西解外推到与我们计算速度/压力的时间相同。我们做一个简单的线性外推，即给定从上次计算达西解到现在的宏观时间步长 $dt$ （由 <code>current_macro_time_step</code> 给出），以及 $DT$ 上一个宏观时间步长（由 <code>old_macro_time_step</code> 给出），然后得到 $u^\ast = u_p + dt \frac{u_p-u_{pp}}{DT} = (1+dt/DT)u_p - dt/DT u_{pp}$  ，其中 $u_p$ 和 $u_{pp}$ 是最近两个计算的达西解。我们只需用两行代码就可以实现这个公式。
 * 

 * 
 * 请注意，这里的算法只有在我们至少有两个先前计算的Darcy解，我们可以从中推断出当前的时间，这一点通过要求重新计算前两个时间步骤的Darcy解来保证。
 * 

 * 
 * 
 * @code
 *     else 
 *       { 
 *         darcy_solution = last_computed_darcy_solution; 
 *         darcy_solution.sadd(1 + current_macro_time_step / old_macro_time_step, 
 *                             -current_macro_time_step / old_macro_time_step, 
 *                             second_last_computed_darcy_solution); 
 *       } 
 * 
 * @endcode
 * 
 * 用这样计算出来的速度矢量，根据介绍中讨论的CFL标准计算出最佳时间步长......
 * 

 * 
 * 
 * @code
 *     { 
 *       old_time_step = time_step; 
 * 
 *       const double max_u_F_prime = get_max_u_F_prime(); 
 *       if (max_u_F_prime > 0) 
 *         time_step = porosity * GridTools::minimal_cell_diameter(triangulation) / 
 *                     saturation_degree / max_u_F_prime / 50; 
 *       else 
 *         time_step = end_time - time; 
 *     } 
 * 
 * @endcode
 * 
 * ......然后在我们处理时间步长的时候，还要更新我们使用的宏观时间步长。具体而言，这涉及到。(i) 如果我们刚刚重新计算了达西解，那么之前的宏观时间步长现在是固定的，当前的宏观时间步长，到现在为止，只是当前（微观）时间步长。(ii) 如果我们没有重新计算达西解，那么当前的宏观时间步长刚刚增长了 <code>time_step</code>  。
 * 

 * 
 * 
 * @code
 *     if (solve_for_pressure_and_velocity == true) 
 *       { 
 *         old_macro_time_step     = current_macro_time_step; 
 *         current_macro_time_step = time_step; 
 *       } 
 *     else 
 *       current_macro_time_step += time_step; 
 * 
 * @endcode
 * 
 * 这个函数的最后一步是根据我们刚刚得到的速度场重新计算饱和解。这自然发生在每一个时间步骤中，我们不会跳过这些计算。在计算饱和度的最后，我们投射回允许的区间 $[0,1]$ ，以确保我们的解保持物理状态。
 * 

 * 
 * 
 * @code
 *     { 
 *       std::cout << "   Solving saturation transport equation..." << std::endl; 
 * 
 *       assemble_saturation_system(); 
 * 
 *       SolverControl solver_control(saturation_matrix.m(), 
 *                                    1e-16 * saturation_rhs.l2_norm()); 
 *       SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control); 
 * 
 *       TrilinosWrappers::PreconditionIC preconditioner; 
 *       preconditioner.initialize(saturation_matrix); 
 * 
 *       cg.solve(saturation_matrix, 
 *                saturation_solution, 
 *                saturation_rhs, 
 *                preconditioner); 
 * 
 *       saturation_constraints.distribute(saturation_solution); 
 *       project_back_saturation(); 
 * 
 *       std::cout << "        ..." << solver_control.last_step() 
 *                 << " CG iterations." << std::endl; 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimrefine_mesh"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::refine_mesh</h3>
 * 

 * 
 * 下一个函数是对网格进行细化和粗化。它的工作分三块进行。(i) 计算细化指标，方法是通过使用各自的时间步长（如果这是第一个时间步长，则取唯一的解决方案），从前两个时间步长中线性推断出的解决方案向量的梯度。(ii) 在梯度大于或小于某一阈值的单元中标记出细化和粗化的单元，保留网格细化的最小和最大水平。(iii) 将解决方案从旧网格转移到新网格。这些都不是特别困难。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::refine_mesh(const unsigned int min_grid_level, 
 *                                              const unsigned int max_grid_level) 
 *   { 
 *     Vector<double> refinement_indicators(triangulation.n_active_cells()); 
 *     { 
 *       const QMidpoint<dim>        quadrature_formula; 
 *       FEValues<dim>               fe_values(saturation_fe, 
 *                               quadrature_formula, 
 *                               update_gradients); 
 *       std::vector<Tensor<1, dim>> grad_saturation(1); 
 * 
 *       TrilinosWrappers::MPI::Vector extrapolated_saturation_solution( 
 *         saturation_solution); 
 *       if (timestep_number != 0) 
 *         extrapolated_saturation_solution.sadd((1. + time_step / old_time_step), 
 *                                               time_step / old_time_step, 
 *                                               old_saturation_solution); 
 * 
 *       for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *         { 
 *           const unsigned int cell_no = cell->active_cell_index(); 
 *           fe_values.reinit(cell); 
 *           fe_values.get_function_gradients(extrapolated_saturation_solution, 
 *                                            grad_saturation); 
 * 
 *           refinement_indicators(cell_no) = grad_saturation[0].norm(); 
 *         } 
 *     } 
 * 
 *     { 
 *       for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *         { 
 *           const unsigned int cell_no = cell->active_cell_index(); 
 *           cell->clear_coarsen_flag(); 
 *           cell->clear_refine_flag(); 
 * 
 *           if ((static_cast<unsigned int>(cell->level()) < max_grid_level) && 
 *               (std::fabs(refinement_indicators(cell_no)) > 
 *                saturation_refinement_threshold)) 
 *             cell->set_refine_flag(); 
 *           else if ((static_cast<unsigned int>(cell->level()) > 
 *                     min_grid_level) && 
 *                    (std::fabs(refinement_indicators(cell_no)) < 
 *                     0.5 * saturation_refinement_threshold)) 
 *             cell->set_coarsen_flag(); 
 *         } 
 *     } 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 * 
 *     { 
 *       std::vector<TrilinosWrappers::MPI::Vector> x_saturation(3); 
 *       x_saturation[0] = saturation_solution; 
 *       x_saturation[1] = old_saturation_solution; 
 *       x_saturation[2] = saturation_matching_last_computed_darcy_solution; 
 * 
 *       std::vector<TrilinosWrappers::MPI::BlockVector> x_darcy(2); 
 *       x_darcy[0] = last_computed_darcy_solution; 
 *       x_darcy[1] = second_last_computed_darcy_solution; 
 * 
 *       SolutionTransfer<dim, TrilinosWrappers::MPI::Vector> saturation_soltrans( 
 *         saturation_dof_handler); 
 * 
 *       SolutionTransfer<dim, TrilinosWrappers::MPI::BlockVector> darcy_soltrans( 
 *         darcy_dof_handler); 
 * 
 *       triangulation.prepare_coarsening_and_refinement(); 
 *       saturation_soltrans.prepare_for_coarsening_and_refinement(x_saturation); 
 * 
 *       darcy_soltrans.prepare_for_coarsening_and_refinement(x_darcy); 
 * 
 *       triangulation.execute_coarsening_and_refinement(); 
 *       setup_dofs(); 
 * 
 *       std::vector<TrilinosWrappers::MPI::Vector> tmp_saturation(3); 
 *       tmp_saturation[0].reinit(saturation_solution); 
 *       tmp_saturation[1].reinit(saturation_solution); 
 *       tmp_saturation[2].reinit(saturation_solution); 
 *       saturation_soltrans.interpolate(x_saturation, tmp_saturation); 
 * 
 *       saturation_solution                              = tmp_saturation[0]; 
 *       old_saturation_solution                          = tmp_saturation[1]; 
 *       saturation_matching_last_computed_darcy_solution = tmp_saturation[2]; 
 * 
 *       saturation_constraints.distribute(saturation_solution); 
 *       saturation_constraints.distribute(old_saturation_solution); 
 *       saturation_constraints.distribute( 
 *         saturation_matching_last_computed_darcy_solution); 
 * 
 *       std::vector<TrilinosWrappers::MPI::BlockVector> tmp_darcy(2); 
 *       tmp_darcy[0].reinit(darcy_solution); 
 *       tmp_darcy[1].reinit(darcy_solution); 
 *       darcy_soltrans.interpolate(x_darcy, tmp_darcy); 
 * 
 *       last_computed_darcy_solution        = tmp_darcy[0]; 
 *       second_last_computed_darcy_solution = tmp_darcy[1]; 
 * 
 *       darcy_constraints.distribute(last_computed_darcy_solution); 
 *       darcy_constraints.distribute(second_last_computed_darcy_solution); 
 * 
 *       rebuild_saturation_matrix = true; 
 *     } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimoutput_results"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::output_results</h3>
 * 

 * 
 * 这个函数生成图形输出。它实质上是对  step-31  中实现的复制。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::output_results() const 
 *   { 
 *     const FESystem<dim> joint_fe(darcy_fe, 1, saturation_fe, 1); 
 *     DoFHandler<dim>     joint_dof_handler(triangulation); 
 *     joint_dof_handler.distribute_dofs(joint_fe); 
 *     Assert(joint_dof_handler.n_dofs() == 
 *              darcy_dof_handler.n_dofs() + saturation_dof_handler.n_dofs(), 
 *            ExcInternalError()); 
 * 
 *     Vector<double> joint_solution(joint_dof_handler.n_dofs()); 
 * 
 *     { 
 *       std::vector<types::global_dof_index> local_joint_dof_indices( 
 *         joint_fe.n_dofs_per_cell()); 
 *       std::vector<types::global_dof_index> local_darcy_dof_indices( 
 *         darcy_fe.n_dofs_per_cell()); 
 *       std::vector<types::global_dof_index> local_saturation_dof_indices( 
 *         saturation_fe.n_dofs_per_cell()); 
 * 
 *       auto       joint_cell      = joint_dof_handler.begin_active(); 
 *       const auto joint_endc      = joint_dof_handler.end(); 
 *       auto       darcy_cell      = darcy_dof_handler.begin_active(); 
 *       auto       saturation_cell = saturation_dof_handler.begin_active(); 
 * 
 *       for (; joint_cell != joint_endc; 
 *            ++joint_cell, ++darcy_cell, ++saturation_cell) 
 *         { 
 *           joint_cell->get_dof_indices(local_joint_dof_indices); 
 *           darcy_cell->get_dof_indices(local_darcy_dof_indices); 
 *           saturation_cell->get_dof_indices(local_saturation_dof_indices); 
 * 
 *           for (unsigned int i = 0; i < joint_fe.n_dofs_per_cell(); ++i) 
 *             if (joint_fe.system_to_base_index(i).first.first == 0) 
 *               { 
 *                 Assert(joint_fe.system_to_base_index(i).second < 
 *                          local_darcy_dof_indices.size(), 
 *                        ExcInternalError()); 
 *                 joint_solution(local_joint_dof_indices[i]) = darcy_solution( 
 *                   local_darcy_dof_indices[joint_fe.system_to_base_index(i) 
 *                                             .second]); 
 *               } 
 *             else 
 *               { 
 *                 Assert(joint_fe.system_to_base_index(i).first.first == 1, 
 *                        ExcInternalError()); 
 *                 Assert(joint_fe.system_to_base_index(i).second < 
 *                          local_darcy_dof_indices.size(), 
 *                        ExcInternalError()); 
 *                 joint_solution(local_joint_dof_indices[i]) = 
 *                   saturation_solution( 
 *                     local_saturation_dof_indices 
 *                       [joint_fe.system_to_base_index(i).second]); 
 *               } 
 *         } 
 *     } 
 *     std::vector<std::string> joint_solution_names(dim, "velocity"); 
 *     joint_solution_names.emplace_back("pressure"); 
 *     joint_solution_names.emplace_back("saturation"); 
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(joint_dof_handler); 
 *     data_out.add_data_vector(joint_solution, 
 *                              joint_solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 * 
 *     data_out.build_patches(); 
 * 
 *     std::string filename = 
 *       "solution-" + Utilities::int_to_string(timestep_number, 5) + ".vtu"; 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Toolfunctions"></a> 
 * <h3>Tool functions</h3>
 * 
 * <a name="TwoPhaseFlowProblemdimdetermine_whether_to_solve_for_pressure_and_velocity"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity</h4>
 * 

 * 
 * 这个函数实现了自适应运算符拆分的后验标准。考虑到我们在上面实现其他函数的方式，并考虑到论文中得出的准则公式，该函数是相对简单的。
 * 

 * 
 * 如果我们决定要采用原始的IMPES方法，即在每个时间步长中求解Darcy方程，那么可以通过将阈值 <code>AOS_threshold</code> （默认为 $5.0$ ）设置为0来实现，从而迫使该函数总是返回true。
 * 

 * 
 * 最后，请注意，该函数在前两个时间步骤中无条件地返回真，以确保我们在跳过达西系统的解时总是至少解了两次，从而允许我们从 <code>solve()</code> 中的最后两次解中推算出速度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   bool TwoPhaseFlowProblem< 
 *     dim>::determine_whether_to_solve_for_pressure_and_velocity() const 
 *   { 
 *     if (timestep_number <= 2) 
 *       return true; 
 * 
 *     const QGauss<dim>  quadrature_formula(saturation_degree + 2); 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(saturation_fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points); 
 * 
 *     std::vector<double> old_saturation_after_solving_pressure(n_q_points); 
 *     std::vector<double> present_saturation(n_q_points); 
 * 
 *     std::vector<Tensor<2, dim>> k_inverse_values(n_q_points); 
 * 
 *     double max_global_aop_indicator = 0.0; 
 * 
 *     for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *       { 
 *         double max_local_mobility_reciprocal_difference = 0.0; 
 *         double max_local_permeability_inverse_l1_norm   = 0.0; 
 * 
 *         fe_values.reinit(cell); 
 *         fe_values.get_function_values( 
 *           saturation_matching_last_computed_darcy_solution, 
 *           old_saturation_after_solving_pressure); 
 *         fe_values.get_function_values(saturation_solution, present_saturation); 
 * 
 *         k_inverse.value_list(fe_values.get_quadrature_points(), 
 *                              k_inverse_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             const double mobility_reciprocal_difference = std::fabs( 
 *               mobility_inverse(present_saturation[q], viscosity) - 
 *               mobility_inverse(old_saturation_after_solving_pressure[q], 
 *                                viscosity)); 
 * 
 *             max_local_mobility_reciprocal_difference = 
 *               std::max(max_local_mobility_reciprocal_difference, 
 *                        mobility_reciprocal_difference); 
 * 
 *             max_local_permeability_inverse_l1_norm = 
 *               std::max(max_local_permeability_inverse_l1_norm, 
 *                        l1_norm(k_inverse_values[q])); 
 *           } 
 * 
 *         max_global_aop_indicator = 
 *           std::max(max_global_aop_indicator, 
 *                    (max_local_mobility_reciprocal_difference * 
 *                     max_local_permeability_inverse_l1_norm)); 
 *       } 
 * 
 *     return (max_global_aop_indicator > AOS_threshold); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimproject_back_saturation"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::project_back_saturation</h4>
 * 

 * 
 * 下一个函数只是确保饱和度值始终保持在  $[0,1]$  的物理合理范围内。虽然连续方程保证了这一点，但离散方程并没有。然而，如果我们允许离散解逃脱这个范围，我们就会遇到麻烦，因为像 $F(S)$ 和 $F'(S)$ 这样的项会产生不合理的结果（例如 $F'(S)<0$ 为 $S<0$ ，这将意味着润湿液相的流动方向为<i>against</i>的散流体速度））。因此，在每个时间步骤结束时，我们只需将饱和场投射回物理上合理的区域。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::project_back_saturation() 
 *   { 
 *     for (unsigned int i = 0; i < saturation_solution.size(); ++i) 
 *       if (saturation_solution(i) < 0.2) 
 *         saturation_solution(i) = 0.2; 
 *       else if (saturation_solution(i) > 1) 
 *         saturation_solution(i) = 1; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimget_max_u_F_prime"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::get_max_u_F_prime</h4>
 * 

 * 
 * 另一个比较简单的辅助函数。计算总速度乘以分数流函数的导数的最大值，即计算  $\|\mathbf{u} F'(S)\|_{L_\infty(\Omega)}$  。这个项既用于时间步长的计算，也用于人工黏度中熵留项的正常化。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double TwoPhaseFlowProblem<dim>::get_max_u_F_prime() const 
 *   { 
 *     const QGauss<dim>  quadrature_formula(darcy_degree + 2); 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> darcy_fe_values(darcy_fe, quadrature_formula, update_values); 
 *     FEValues<dim> saturation_fe_values(saturation_fe, 
 *                                        quadrature_formula, 
 *                                        update_values); 
 * 
 *     std::vector<Vector<double>> darcy_solution_values(n_q_points, 
 *                                                       Vector<double>(dim + 1)); 
 *     std::vector<double>         saturation_values(n_q_points); 
 * 
 *     double max_velocity_times_dF_dS = 0; 
 * 
 *     auto       cell            = darcy_dof_handler.begin_active(); 
 *     const auto endc            = darcy_dof_handler.end(); 
 *     auto       saturation_cell = saturation_dof_handler.begin_active(); 
 *     for (; cell != endc; ++cell, ++saturation_cell) 
 *       { 
 *         darcy_fe_values.reinit(cell); 
 *         saturation_fe_values.reinit(saturation_cell); 
 * 
 *         darcy_fe_values.get_function_values(darcy_solution, 
 *                                             darcy_solution_values); 
 *         saturation_fe_values.get_function_values(old_saturation_solution, 
 *                                                  saturation_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             Tensor<1, dim> velocity; 
 *             for (unsigned int i = 0; i < dim; ++i) 
 *               velocity[i] = darcy_solution_values[q](i); 
 * 
 *             const double dF_dS = 
 *               fractional_flow_derivative(saturation_values[q], viscosity); 
 * 
 *             max_velocity_times_dF_dS = 
 *               std::max(max_velocity_times_dF_dS, velocity.norm() * dF_dS); 
 *           } 
 *       } 
 * 
 *     return max_velocity_times_dF_dS; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimget_extrapolated_saturation_range"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range</h4>
 * 

 * 
 * 为了计算稳定化项，我们需要知道饱和变量的范围。与 step-31 不同，这个范围很容易被区间 $[0,1]$ 所约束，但是我们可以通过在正交点的集合上循环，看看那里的值是多少，从而做得更好。如果可以的话，也就是说，如果周围至少有两个时间步长，我们甚至可以把这些值推算到下一个时间步长。
 * 

 * 
 * 和以前一样，这个函数是在对  step-31  进行最小修改后取的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<double, double> 
 *   TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range() const 
 *   { 
 *     const QGauss<dim>  quadrature_formula(saturation_degree + 2); 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FEValues<dim> fe_values(saturation_fe, quadrature_formula, update_values); 
 *     std::vector<double> old_saturation_values(n_q_points); 
 *     std::vector<double> old_old_saturation_values(n_q_points); 
 * 
 *     if (timestep_number != 0) 
 *       { 
 *         double min_saturation = std::numeric_limits<double>::max(), 
 *                max_saturation = -std::numeric_limits<double>::max(); 
 * 
 *         for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *           { 
 *             fe_values.reinit(cell); 
 *             fe_values.get_function_values(old_saturation_solution, 
 *                                           old_saturation_values); 
 *             fe_values.get_function_values(old_old_saturation_solution, 
 *                                           old_old_saturation_values); 
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q) 
 *               { 
 *                 const double saturation = 
 *                   (1. + time_step / old_time_step) * old_saturation_values[q] - 
 *                   time_step / old_time_step * old_old_saturation_values[q]; 
 * 
 *                 min_saturation = std::min(min_saturation, saturation); 
 *                 max_saturation = std::max(max_saturation, saturation); 
 *               } 
 *           } 
 * 
 *         return std::make_pair(min_saturation, max_saturation); 
 *       } 
 *     else 
 *       { 
 *         double min_saturation = std::numeric_limits<double>::max(), 
 *                max_saturation = -std::numeric_limits<double>::max(); 
 * 
 *         for (const auto &cell : saturation_dof_handler.active_cell_iterators()) 
 *           { 
 *             fe_values.reinit(cell); 
 *             fe_values.get_function_values(old_saturation_solution, 
 *                                           old_saturation_values); 
 * 
 *             for (unsigned int q = 0; q < n_q_points; ++q) 
 *               { 
 *                 const double saturation = old_saturation_values[q]; 
 * 
 *                 min_saturation = std::min(min_saturation, saturation); 
 *                 max_saturation = std::max(max_saturation, saturation); 
 *               } 
 *           } 
 * 
 *         return std::make_pair(min_saturation, max_saturation); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimcompute_viscosity"></a> 
 * <h4>TwoPhaseFlowProblem<dim>::compute_viscosity</h4>
 * 

 * 
 * 最后一个工具函数是用来计算给定单元上的人工粘度的。如果你面前有它的公式，这并不特别复杂，看一下  step-31  中的实现。与那个教程程序的主要区别是，这里的速度不是简单的 $\mathbf u$ ，而是 $\mathbf u F'(S)$ ，一些公式需要做相应的调整。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double TwoPhaseFlowProblem<dim>::compute_viscosity( 
 *     const std::vector<double> &        old_saturation, 
 *     const std::vector<double> &        old_old_saturation, 
 *     const std::vector<Tensor<1, dim>> &old_saturation_grads, 
 *     const std::vector<Tensor<1, dim>> &old_old_saturation_grads, 
 *     const std::vector<Vector<double>> &present_darcy_values, 
 *     const double                       global_max_u_F_prime, 
 *     const double                       global_S_variation, 
 *     const double                       cell_diameter) const 
 *   { 
 *     const double beta  = .4 * dim; 
 *     const double alpha = 1; 
 * 
 *     if (global_max_u_F_prime == 0) 
 *       return 5e-3 * cell_diameter; 
 * 
 *     const unsigned int n_q_points = old_saturation.size(); 
 * 
 *     double max_residual             = 0; 
 *     double max_velocity_times_dF_dS = 0; 
 * 
 *     const bool use_dF_dS = true; 
 * 
 *     for (unsigned int q = 0; q < n_q_points; ++q) 
 *       { 
 *         Tensor<1, dim> u; 
 *         for (unsigned int d = 0; d < dim; ++d) 
 *           u[d] = present_darcy_values[q](d); 
 * 
 *         const double dS_dt = porosity * 
 *                              (old_saturation[q] - old_old_saturation[q]) / 
 *                              old_time_step; 
 * 
 *         const double dF_dS = fractional_flow_derivative( 
 *           (old_saturation[q] + old_old_saturation[q]) / 2.0, viscosity); 
 * 
 *         const double u_grad_S = 
 *           u * dF_dS * (old_saturation_grads[q] + old_old_saturation_grads[q]) / 
 *           2.0; 
 * 
 *         const double residual = 
 *           std::abs((dS_dt + u_grad_S) * 
 *                    std::pow((old_saturation[q] + old_old_saturation[q]) / 2, 
 *                             alpha - 1.)); 
 * 
 *         max_residual = std::max(residual, max_residual); 
 *         max_velocity_times_dF_dS = 
 *           std::max(std::sqrt(u * u) * (use_dF_dS ? std::max(dF_dS, 1.) : 1), 
 *                    max_velocity_times_dF_dS); 
 *       } 
 * 
 *     const double c_R            = 1.0; 
 *     const double global_scaling = c_R * porosity * 
 *                                   (global_max_u_F_prime)*global_S_variation / 
 *                                   std::pow(global_Omega_diameter, alpha - 2.); 
 * 
 *     return (beta * 
 *             (max_velocity_times_dF_dS)*std::min(cell_diameter, 
 *                                                 std::pow(cell_diameter, alpha) * 
 *                                                   max_residual / 
 *                                                   global_scaling)); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TwoPhaseFlowProblemdimrun"></a> 
 * <h3>TwoPhaseFlowProblem<dim>::run</h3>
 * 

 * 
 * 除了 <code>solve()</code> 之外，这个函数是这个程序的主要功能，因为它控制了迭代的时间，以及何时将解决方案写入输出文件，何时进行网格细化。
 * 

 * 
 * 除了启动代码通过 <code>goto start_time_iteration</code> 标签循环回到函数的开头外，一切都应该是相对简单的。无论如何，它模仿了  step-31  中的相应函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void TwoPhaseFlowProblem<dim>::run() 
 *   { 
 *     const unsigned int initial_refinement     = (dim == 2 ? 5 : 2); 
 *     const unsigned int n_pre_refinement_steps = (dim == 2 ? 3 : 2); 
 * 
 *     GridGenerator::hyper_cube(triangulation, 0, 1); 
 *     triangulation.refine_global(initial_refinement); 
 *     global_Omega_diameter = GridTools::diameter(triangulation); 
 * 
 *     setup_dofs(); 
 * 
 *     unsigned int pre_refinement_step = 0; 
 * 
 *   start_time_iteration: 
 * 
 *     VectorTools::project(saturation_dof_handler, 
 *                          saturation_constraints, 
 *                          QGauss<dim>(saturation_degree + 2), 
 *                          SaturationInitialValues<dim>(), 
 *                          old_saturation_solution); 
 * 
 *     time_step = old_time_step = 0; 
 *     current_macro_time_step = old_macro_time_step = 0; 
 * 
 *     time = 0; 
 * 
 *     do 
 *       { 
 *         std::cout << "Timestep " << timestep_number << ":  t=" << time 
 *                   << ", dt=" << time_step << std::endl; 
 * 
 *         solve(); 
 * 
 *         std::cout << std::endl; 
 * 
 *         if (timestep_number % 200 == 0) 
 *           output_results(); 
 * 
 *         if (timestep_number % 25 == 0) 
 *           refine_mesh(initial_refinement, 
 *                       initial_refinement + n_pre_refinement_steps); 
 * 
 *         if ((timestep_number == 0) && 
 *             (pre_refinement_step < n_pre_refinement_steps)) 
 *           { 
 *             ++pre_refinement_step; 
 *             goto start_time_iteration; 
 *           } 
 * 
 *         time += time_step; 
 *         ++timestep_number; 
 * 
 *         old_old_saturation_solution = old_saturation_solution; 
 *         old_saturation_solution     = saturation_solution; 
 *       } 
 *     while (time <= end_time); 
 *   } 
 * } // namespace Step43 
 * 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main()</code> function</h3>
 * 

 * 
 * 主函数看起来与所有其他程序几乎一样。对于使用Trilinos的程序来说，需要初始化MPI子系统--即使是那些实际上没有并行运行的程序--在  step-31  中有解释。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step43; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization( 
 *         argc, argv, numbers::invalid_unsigned_int); 
 * 
 * @endcode
 * 
 * 这个程序只能在串行中运行。否则，将抛出一个异常。
 * 

 * 
 * 
 * @code
 *       AssertThrow(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD) == 1, 
 *                   ExcMessage( 
 *                     "This program can only be run in serial, use ./step-43")); 
 * 
 *       TwoPhaseFlowProblem<2> two_phase_flow_problem(1); 
 *       two_phase_flow_problem.run(); 
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
examples/step-43/doc/results.dox



<a name="Results"></a><h1>Results</h1>



这个程序的输出与第21步的输出其实没有什么不同：毕竟它解决的是同一个问题。更重要的是定量指标，如解决方案的准确性以及计算所需的时间。这些在本页顶部列出的两份出版物中都有详细记载，我们在此不再重复。

也就是说，如果没有几张好的照片，任何教程程序都是不完整的，所以这里有一些三维运行的输出。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.velocity.png" alt="">
	<p align="center">
        Velocity vectors of flow through the porous medium with random
        permeability model. Streaming paths of high permeability and resulting
        high velocity are clearly visible.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.streamlines.png" alt="">
	<p align="center">
        Streamlines colored by the saturation along the streamline path. Blue
        streamlines indicate low saturations, i.e., the flow along these
	streamlines must be slow or else more fluid would have been
        transported along them. On the other hand, green paths indicate high
        velocities since the fluid front has already reached further into the
        domain.
	</p>
    </td>
  </tr>
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.saturation.png" alt="">
	<p align="center">
        Streamlines with a volume rendering of the saturation, showing how far
        the fluid front has advanced at this time.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-43.3d.mesh.png" alt="">
	<p align="center">
	Surface of the mesh showing the adaptive refinement along the front.
	</p>
    </td>
  </tr>
</table> 


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


人们对这个程序的主要反对意见是它仍然太慢了：在合理的细网格上的三维计算实在是太昂贵了，无法以合理的快速周转来进行常规计算。这与我们写step-31时的情况相似，这个程序从它那里得到了很多灵感。解决办法也是类似的，因为它也在那里。我们需要以类似于从第31步衍生出第32步的方式来并行化这个程序。事实上，步骤32中使用的所有技术也可以转移到这个程序中，使程序立即在几十或几百个处理器上运行。

一个不同的方向是使该程序与许多其他多孔介质的应用更加相关。具体来说，一个途径是去找多孔介质流动模拟器的主要用户，即石油工业。在那里，该领域的应用以多相流（即超过我们这里的两相）为主，以及它们之间可能发生的反应（或任何其他相的质量交换方式，如通过溶解和从油相中冒出的气体）。此外，气体的存在往往会导致流体的可压缩性效应。这些效应通常共同组成了广泛使用的 "黑油模型"。在考虑储层中石油的控制性燃烧以提高压力和温度时，多相之间的真正反应也在油藏模型中发挥作用。不过，这些问题要复杂得多，留待今后的项目研究。

最后，从数学的角度来看，我们得出了在某一时间步长重新计算速度/压力解的标准，其前提是我们要把在当前时间步长会得到的解与上次实际解这个系统时计算的解进行比较。然而，在程序中，每当我们没有重新计算解决方案时，我们并不只是使用之前计算的解决方案，而是从之前两次求解系统的结果中推算出来。因此，该标准被悲观地表述为：我们真正应该比较的是在当前时间步长得到的解与外推的解。在这方面重述该定理是一个练习。

也有其他方法可以扩展这个程序的数学基础；例如，人们可以说，我们关心的不是速度，而实际上是饱和度。因此，人们可能会问，我们在这里用来决定 $\mathbf u$ 是否需要重新计算的标准是否合适；例如，人们可能会提出，决定一个错误的速度场事实上是否会影响饱和方程的解（以及影响的程度）也很重要。这自然会导致敏感性分析。

从算法的角度来看，我们在这里使用了一个工程中经常使用的细化标准，即通过查看解的梯度。然而，如果你检查解决方案，你会发现它几乎在所有地方都迅速导致细化，甚至在明显没有必要的区域：因此经常使用并不需要暗示它是一个有用的标准开始。另一方面，用一个不同的、更好的标准来取代这个标准应该不是很困难。例如，许多其他程序中使用的KellyErrorEstimator类当然也应该适用于当前的问题。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-43.cc"
*/
