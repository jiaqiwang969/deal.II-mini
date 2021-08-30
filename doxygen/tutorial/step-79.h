/**
@page step_79 The step-79 tutorial program
This tutorial depends on step-8, step-15.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#SolidIsotropicMaterialwithPenalization">Solid Isotropic Material with Penalization</a>
        <li><a href="#ElasticityEquation">Elasticity Equation</a>
        <li><a href="#Makingthesolutionmeshindependent">Making the solution mesh-independent</a>
        <li><a href="#CompleteProblemFormulation">Complete Problem Formulation</a>
        <li><a href="#Solutionprocedure">Solution procedure</a>
        <li><a href="#Discretization">Discretization</a>
        <li><a href="#NonlinearAlgorithm">Nonlinear Algorithm</a>
        <li><a href="#MeritFunction">Merit Function</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Preliminaries">Preliminaries</a>
        <li><a href="#TheSANDTopOptmainclass">The SANDTopOpt main class</a>
        <li><a href="#Constructorandsetupfunctions">Constructor and set-up functions</a>
        <li><a href="#Settingupblockmatricesandvectors">Setting up block matrices and vectors</a>
        <li><a href="#Creatingthefiltermatrix">Creating the filter matrix</a>
        <li><a href="#AssemblingtheNewtonmatrix">Assembling the Newton matrix</a>
        <li><a href="#SolvingtheNewtonlinearsystem">Solving the Newton linear system</a>
        <li><a href="#Detailsoftheoptimizationalgorithm">Details of the optimization algorithm</a>
      <ul>
        <li><a href="#Computingsteplengths">Computing step lengths</a>
        <li><a href="#Computingresiduals">Computing residuals</a>
        <li><a href="#Computingthemeritfunction">Computing the merit function</a>
        <li><a href="#Findingasearchdirection">Finding a search direction</a>
        <li><a href="#Computingascaledstep">Computing a scaled step</a>
        <li><a href="#Checkingforconvergence">Checking for convergence</a>
      </ul>
        <li><a href="#Postprocessingthesolution">Postprocessing the solution</a>
        <li><a href="#Therunfunctiondrivingtheoverallalgorithm">The run() function driving the overall algorithm</a>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#TestProblem">Test Problem</a>
      <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-79/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


弹性介质的拓扑优化是一种用于优化承受某种载荷的结构的技术。理想情况下，我们希望通过选择一个放置材料的区域 $E$ ，使置于结构上的最大应力最小化。换句话说。

@f[
  \text{minimize}\| \boldsymbol{\sigma} (\mathbf{u}) \|_\infty


@f]



@f[
  \text{subject to } |E|\leq V_{\max},


@f]



@f[
  \text{and } \nabla \cdot \boldsymbol{\sigma} + \mathbf{F} = \mathbf{0}.


@f]



这里， $\boldsymbol{\sigma} = \mathbf{C} : \boldsymbol{\varepsilon}(\mathbf{u})$ 是由外力 $\mathbf F$ 引起的体内应力，为了简单起见，我们假设材料是线性弹性的，因此 $\mathbf{C}$ 是应力-应变张量， $\boldsymbol{\varepsilon}(\mathbf{u})=\frac{1}{2} (\nabla \mathbf{u} + (\nabla\mathbf{u})^T)$ 是作为位移 $\mathbf{u}$ 函数的小变形应变--关于线性弹性的详情，见步骤8 和步骤17。在上面的表述中， $V_\text{max}$ 是我们愿意为构建物体提供的最大材料量。最后一个约束条件是与应力 $\boldsymbol{\sigma}$ 和力 $\mathbf F$ 有关的偏微分方程，它只是稳态力平衡。

也就是说，上面的无穷大准则产生了一个问题：作为材料位置的函数，这个目标函数必然是不可微分的，使优化的前景相当暗淡。因此，取而代之的是，拓扑优化的一个常见方法是通过优化一个相关的问题来找到一个近似的解决方案：我们希望最小化应变能量。这是对物体因变形而储存的势能的衡量，同时也是对结构总变形的衡量。

@f[
  \text{minimize  } \int_E \frac{1}{2}\boldsymbol{\sigma} : \boldsymbol{\varepsilon} dV


@f]



@f[
  \text{subject to } \|E\| \leq V_{\max}


@f]



@f[
  \text{and } \nabla \cdot \boldsymbol{\sigma} + \mathbf{F} = \mathbf{0}


@f]



目标函数的值是用有限元方法计算的，其中的解决方案是位移。这被放置在一个非线性求解器的循环中，求解一个表示材料放置的向量。

<a name="SolidIsotropicMaterialwithPenalization"></a><h3>Solid Isotropic Material with Penalization</h3>


在实际操作中，我们只能建造材料在任何给定的点上要么存在，要么不存在的物体--也就是说，我们会有一个描述材料填充区域的指标函数 $\rho_E(\mathbf{x})\in \{0,1\}$ ，并且我们想通过优化问题找到这个指标。在这种情况下，优化问题变成了组合性的，而且解决起来非常昂贵。取而代之的是，我们使用一种叫做各向同性的固体材料与惩罚的方法，或SIMP。   @cite Bendse2004 

SIMP方法是基于一个想法，即允许材料存在于密度 $\rho$ 在0和1之间的位置。密度为0表明材料不存在，它不是结构的一部分，而密度为1表明材料存在。0和1之间的值并不反映我们在现实世界中可以创造的设计，但允许我们将组合问题变成一个连续问题。然后我们看一下密度值  $\rho$  ，约束条件是  $0 < \rho_{\min} \leq \rho \leq 1$  。最小值 $\rho_{\min}$ ，通常选择在 $10^{-3}$ 左右，避免了出现无限应变能量的可能性，但小到足以提供准确的结果。

这种 "密度 "对介质弹性的影响的直接应用是简单地将介质的刚度张量 $\mathbf{C}_0$ 乘以给定的密度，即 $\mathbf{C} = \rho \mathbf{C}_0$  。然而，这种方法经常给出密度值离0和1都很远的最佳解决方案。由于人们希望找到一个现实世界的解决方案，即材料要么存在，要么不存在，因此对这些介于两者之间的值进行惩罚。一个简单有效的方法是将刚度张量乘以密度，并将其提高到某个整数功率的惩罚参数  $p$  ，因此  $\mathbf{C} = \rho^p \mathbf{C}_0$  。这使得远离0或1的密度值变得不那么有效。已经证明，使用 $p=3$ 足够高，可以产生'黑白'的解决方案：也就是说，可以得到最佳的解决方案，其中材料在所有点上要么存在，要么不存在。

更多的材料应该总是提供一个具有较低应变能量的结构，因此不等式约束可以被看作是一个等式，其中使用的总体积是最大体积。

使用这种密度思想也使我们能够重新构建优化问题的体积约束。使用SIMP后，优化问题就变成了以下内容。

@f[
  \text{minimize  } \int_\Omega \frac{1}{2}\boldsymbol{\sigma}(\rho) : \boldsymbol{\varepsilon}(\rho) d\Omega


@f]



@f[
  \text{subject to } \int_\Omega \rho(x) d\Omega= V_{\max},


@f]



@f[
  0<\rho_{\min}\leq \rho(x) \leq 1,


@f]



@f[


  \nabla \cdot \boldsymbol{\sigma}(\rho) + \mathbf{F} = 0 \quad \text{on } \Omega


@f]

最后一个约束，即线性动量的平衡（我们将称之为弹性方程），给出了一种在给定密度 $\boldsymbol{\sigma}$ 和 $\boldsymbol{\varepsilon}$ 的情况下寻找 $\rho$  的方法。

<a name="ElasticityEquation"></a><h3>Elasticity Equation</h3> 在与时间无关的极限中，弹性方程为


@f[
  \nabla \cdot \boldsymbol{\sigma} + \mathbf{F} = \mathbf{0} .


@f]

在我们将关注的情况下，我们将假设介质具有线性材料响应，在这种情况下，我们有

@f[
  \boldsymbol{\sigma} = \mathbf{C} : \boldsymbol{\varepsilon} = \rho^p \mathbf{C}_0 : \boldsymbol{\varepsilon}(\mathbf{u})
   = \rho^p \mathbf{C}_0 : \left[\frac{1}{2} (\nabla \mathbf{u} + (\nabla \mathbf{u})^T) \right] .


@f]

在我们下面要做的一切中，我们将始终把位移场 $\mathbf{u}$ 视为唯一的解变量，而不是把 $\mathbf{u}$ 和 $\boldsymbol{\sigma}$ 视为解变量（像在混合公式中那样）。

此外，我们将假设材料是线性各向同性的，在这种情况下，应力-应变张量可以用Lam&eacute;参数 $\lambda,\mu$ 来表示，例如

@f{align}
  \boldsymbol{\sigma} &= \rho^p (\lambda \text{tr}(\boldsymbol{\varepsilon}) \mathbf{I} + 2 \mu \boldsymbol{\varepsilon}) , \\
  \sigma_{i,j} &= \rho^p (\lambda \varepsilon_{k,k} \delta_{i,j} + 2 \mu \varepsilon_{i,j}) .


@f}

参见步骤8，了解这种转变的原理。

对目标函数进行分项积分，得到

@f[
  \int_\Omega \boldsymbol{\sigma}(\rho) : (\nabla \mathbf{u} + (\nabla \mathbf{u}))^T  d\Omega+
  \int_\Omega (\nabla \cdot \boldsymbol{\sigma}(\rho)) \cdot \mathbf{u}  d\Omega=
  \int_{\partial \Omega} \mathbf{t} \cdot \mathbf{u} d\partial\Omega ,


@f]

然后将线性弹性方程代入其中，可以得到

@f[
  \int_\Omega \boldsymbol{\sigma}(\rho) : (\nabla \mathbf{u} + (\nabla \mathbf{u})^T) d\Omega =
  \int_\Omega \mathbf{F}\cdot \mathbf{u} d\Omega+
  \int_{\partial \Omega} \mathbf{t} \cdot \mathbf{u} d\partial\Omega .


@f]

因为我们假设没有身体的力量，这进一步简化为

@f[
  \int_\Omega \boldsymbol{\sigma}(\rho) : (\nabla \mathbf{u} + (\nabla \mathbf{u})^T) d\Omega
  = \int_{\partial \Omega} \mathbf{t} \cdot \mathbf{u} d\partial\Omega,


@f]

这就是我们从现在开始要考虑的治理方程的最终形式。

<a name="Makingthesolutionmeshindependent"></a><h3>Making the solution mesh-independent</h3>


通常情况下，拓扑优化问题的解决方案是依赖于网格的，因此问题是不成立的。这是因为随着网格的进一步细化，往往会形成分形结构。随着网格分辨率的提高，最优解通常会获得越来越小的结构。对于这个问题，有一些相互竞争的解决方法，但对于一阶优化来说，最流行的是灵敏度滤波器，而二阶优化方法则倾向于使用密度滤波器。

由于滤波器会影响应变能量的梯度和Hessian（即目标函数），所以滤波器的选择会对问题的解决产生影响。作为二阶方法的一部分，密度滤波器的工作原理是引入一个未经过滤的密度，我们称之为 $\varrho$  ，然后要求密度是未经过滤的密度的卷积。

@f[
  \rho = H(\varrho).


@f]

这里， $H$ 是一个运算符，因此 $\rho(\mathbf{x})$ 是 $\varrho$ 在 $\mathbf{x}$ 周围区域的某种平均值 -- 即，它是 $\varrho$ 的平滑版本。

这可以防止棋盘效应；滤波器的半径允许用户为我们寻求的最佳结构定义一个有效的最小光束宽度。

<div style="text-align:center;"> <img src="https://www.dealii.org/images/steps/developer/step-79.checkerboard.png" alt="Checkerboarding occurring in an MBB Beam"> </div>

<a name="CompleteProblemFormulation"></a><h3>Complete Problem Formulation</h3>


现在的最小化问题是

@f[
  \min_{\rho,\varrho,\mathbf{u}} \int_{\partial\Omega} \mathbf{u} \cdot \mathbf{t} d\partial\Omega


@f]



@f[
  \text{subject to   } \rho = H(\varrho)


@f]



@f[
  \int_\Omega \rho^p \left(\frac{\mu}{2}\left(\boldsymbol{\varepsilon}(\mathbf{v}):
  \boldsymbol{\varepsilon}(\mathbf{u})) \right) + \lambda \left( \nabla \cdot \mathbf{u} \nabla
  \cdot \mathbf{v} \right)  \right) d\Omega = \int_{\partial \Omega} \mathbf{v} \cdot
  \mathbf{t} d\partial\Omega


@f]



@f[
  \int_\Omega \rho d\Omega= V


@f]



@f[
  0\leq \varrho \leq 1


@f]



处理不等式约束的方法是，首先引入松弛变量，其次使用对数障碍来确保我们得到一个内点方法。惩罚参数将是 $\alpha$  ，下面的松弛变量是<ol>  <li>   $s_1$  -对应于下限的松弛变量  </li>   <li>   $s_2$  -对应于上限的松弛变量。 </li>   </ol>  现在得出以下问题。

@f[
  \min_{\rho,\varrho,\mathbf{u}, s_1, s_2} \int_{\partial\Omega} \mathbf{u} \cdot
  \mathbf{t} d\partial\Omega- \alpha \int_\Omega \left(\log(s_1) + \log(s_2)\right) d\Omega


@f]



@f[
  \text{subject to   } \rho = H(\varrho)


@f]



@f[
  \int_\Omega \rho^p \left(\frac{\mu}{2}\left(\boldsymbol{\varepsilon}(\mathbf{v}):
  \boldsymbol{\varepsilon}(\mathbf{u})) \right) + \lambda \left( \nabla \cdot \mathbf{u} \nabla
  \cdot \mathbf{v} \right)  \right) d\Omega = \int_{\partial \Omega} \mathbf{v} \cdot
  \mathbf{t} d\partial\Omega


@f]



@f[
  \int_\Omega \rho d\Omega = V


@f]



@f[
  \varrho = s_1


@f]



@f[
  1-\varrho = s_2


@f]



有了这些变量，我们就可以按照通常的方法来解决限制性优化问题。我们引入一个拉格朗日，通过将约束条件乘以拉格朗日乘数，将目标函数和约束条件结合起来。具体来说，我们将使用以下符号表示各种约束条件的拉格朗日乘数。<ol>  <li>   $\mathbf{y}_1 $  ：对应于弹性约束的拉格朗日乘数，  </li>   <li>   $y_2$  ：对应于卷积过滤器约束的拉格朗日乘数，  </li>   <li>   $z_1$  ：对应于下层松弛变量的拉格朗日乘数，以及  </li>   <li>   $z_2$  ：对应于上限松弛变量的拉格朗日乘数。   </li>   </ol>  有了这些变量，拉格朗日函数的内容如下。

@f{align}{
  \mathcal{L} =& \int_{\partial\Omega} \mathbf{u} \cdot \mathbf{t} d\partial\Omega


   - \alpha \int_\Omega \left(\log(s_1) + \log(s_2)\right) d\Omega-  \int_\Omega
   \rho^p \left(\frac{\mu}{2}\left(\boldsymbol{\varepsilon}(\mathbf{y}_1):\boldsymbol{\varepsilon}(\mathbf{u}))
   \right) + \lambda \left( \nabla \cdot \mathbf{u} \nabla \cdot \mathbf{y}_1
   \right)\right) d\Omega - \int_{\partial \Omega} \mathbf{y}_1 \cdot \mathbf{t} d\partial\Omega  \\
   & -\int_\Omega y_2 (\rho - H(\varrho)) d\Omega - \int_\Omega z_1 (\varrho-s_1) d\Omega


   - \int_\Omega z_2 (1 - s_2 -\varrho) d\Omega


@f}



然后优化问题的解决方案需要满足所谓的[Karush-Kuhn-Tucker（KKT）条件]（https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions）。拉格朗日相对于其所有参数的导数需要等于零，而且由于我们有不等式约束，我们也有 "互补性 "条件。由于我们这里有一个无穷大的问题，这些条件都涉及到拉格朗日相对于某些测试函数的方向性导数--换句话说，所有这些条件都必须以弱形式表述，因为这通常是有限元方法的基础。

障碍法允许我们最初削弱典型的KKT条件所要求的 "补充松弛"。通常情况下，我们会要求 $s_i z_i = 0$  ，但屏障公式给出的KKT条件是 $s_i z_i = \alpha$  ，其中 $\alpha$  是我们的屏障参数。作为障碍法的一部分，这个参数必须被驱动到接近0，以便对原始问题有一个良好的近似。

在下文中，让我们陈述所有这些条件，其中 $d_{\{\bullet\}}$ 是一个测试函数，它与拉格朗日相对于 $\{\bullet\}$ 函数的变异导数自然成对。为了简单起见，我们引入 $\Gamma$ 来表示边界上受力的部分，并使用诺伊曼边界条件。

<ol>  <li>  静止性。

@f[
  \int_\Omega  - d_\rho y_2 + p\rho^{p-1}d_\rho \left[\lambda
  (\nabla \cdot \mathbf{y}_1) (\nabla \cdot \mathbf{u}) +
  \mu \boldsymbol{\varepsilon}(\mathbf{u}):\boldsymbol{\varepsilon}(\mathbf{y}_1)\right] d\Omega=0\;\;
  \forall d_\rho


@f]



@f[
  \int_\Gamma \mathbf d_\mathbf{u} \cdot \mathbf{t} d\partial\Omega+ \int_\Omega p\rho^{p}
  \left[\lambda (\nabla \cdot \mathbf d_\mathbf{u})( \nabla \cdot \mathbf{y}_1)
  + \mu \boldsymbol{\varepsilon}(\mathbf d_\mathbf{u}):\boldsymbol{\varepsilon}(\mathbf{y}_1)\right] d\Omega=0\;\;
  \forall \mathbf{d}_\mathbf{u}


@f]



@f[
  \int_\Omega -d_\varrho z_1 + d_\varrho z_2 + H(d_\varrho)y_2 d\Omega= 0\;\;\forall
  d_\varrho


@f]

 </li>   <li>  原始的可行性。

@f[
  \int_\Omega \rho^{p}\lambda (\nabla \cdot \mathbf d_{\mathbf{y}_1})
  (\nabla \cdot \mathbf{u}) +  \rho^{p}\mu  \boldsymbol{\varepsilon}(\mathbf
  d_{\mathbf{y}_1}) : \boldsymbol{\varepsilon}(\mathbf{u}) d\Omega - \int_\Gamma \mathbf
  d_{\mathbf{y}_1} \cdot \mathbf{t} d\partial\Omega =0 \;\;\forall \mathbf{d}_{\mathbf{y}_1}


@f]



@f[
  \int_\Omega d_{z_1}(\varrho - s_1) d\Omega = 0\;\;\forall d_{z_1}


@f]



@f[
  \int_\Omega d_{z_z}(1-\varrho-s_2) d\Omega = 0\;\;\forall d_{z_2}


@f]



@f[
  \int_\Omega d_{y_2}(\rho - H(\varrho)) d\Omega = 0\;\;\forall d_{y_2}


@f]

 </li>   <li>  互补性松弛。

@f[
  \int_\Omega d_{s_1}(s_1z_1 - \alpha) d\Omega = 0 \;\;\forall d_{s_1} ,\;\;\;
  \alpha \to 0


@f]



@f[
  \int_\Omega d_{s_2}(s_2z_2 - \alpha) d\Omega = 0  \;\;\forall d_{s_2} ,\;\;\;
  \alpha \to 0


@f]

 </li>   <li>  双重可行性。

@f[
  s_{1,i},s_{2,i},z_{1,i},z_{2,i} \geq 0 \;\;\;\; \forall i


@f]

 </li>  </ol>  。

<a name="Solutionprocedure"></a><h3>Solution procedure</h3>


上面的优化条件除了复杂之外，还属于不容易解决的类型。它们通常是非线性的，而且有些关系也是不等式的。我们将使用牛顿方法计算搜索方向来解决非线性问题，并在下面讨论步长程序时再来讨论如何处理不等式问题。

牛顿方法应用于上述方程的结果是下面列出的方程组。其中，关于 $\{\bullet\}$ 变量的变异导数在 $c_{\{\bullet\}}$ 方向取值。

<ol>  <li>  静止性。这些方程确保我们在受约束时处于目标函数的临界点。

方程式1

@f{align}{
  &\int_\Omega-d_\rho c_{y_2} + p(p-1) \rho^{p-2} d_\rho c_\rho [\lambda \nabla
  \cdot \mathbf{y}_1 \nabla \cdot \mathbf{u} +  \mu  \boldsymbol{\varepsilon}(\mathbf{u})
  \boldsymbol{\varepsilon}(\mathbf{y}_1)]
  + p \rho^{p-1} d_\rho[\lambda \nabla \cdot
  \mathbf{c}_{\mathbf{y}_1} \nabla \cdot \mathbf{u} +   \mu  \boldsymbol{\varepsilon}
  (\mathbf{u}) \boldsymbol{\varepsilon}(\mathbf{c}_{\mathbf{y}_1})]  +  p \rho^{p-1} d_\rho
  [\lambda \nabla \cdot {\mathbf{y}_1} \nabla \cdot \mathbf{c}_\mathbf{u} +
  \mu  \boldsymbol{\varepsilon}(\mathbf{c}_\mathbf{u}) \boldsymbol{\varepsilon}(\mathbf{y}_1)] d\Omega \\
  &= -\int_\Omega -d_\rho z_1 + d_\rho z_2 - d_\rho y_2 + p\rho^{p-1}d_\rho
[\lambda \nabla \cdot \mathbf{y}_1 \nabla \cdot \mathbf{u} + \mu \boldsymbol{\varepsilon}
(\mathbf{u})\boldsymbol{\varepsilon}(\mathbf{y}_1)] d\Omega


@f}



方程式2

@f{align}{
  &\int_\Omega p \rho^{p-1} c_\rho [\lambda \nabla \cdot {\mathbf{y}_1} \nabla
  \cdot \mathbf{d}_\mathbf{u} +  \mu  \boldsymbol{\varepsilon}(\mathbf{d}_\mathbf{u})
  \boldsymbol{\varepsilon}(\mathbf{y})] + \rho^{p} [\lambda \nabla \cdot
  \mathbf{c}_{\mathbf{y}_1} \nabla \cdot \mathbf{d}_\mathbf{u} +  \mu
  \boldsymbol{\varepsilon}(\mathbf{d}_\mathbf{u})\boldsymbol{\varepsilon}(\mathbf{c}_{\mathbf{y}_1})] d\Omega \\
  &= -\int_\Gamma \mathbf{d}_\mathbf{u} \cdot \mathbf{t} -\int_\Omega \rho^{p}
  [\lambda \nabla \cdot \mathbf{y} \nabla \cdot \mathbf{d}_\mathbf{u} + \mu
  \boldsymbol{\varepsilon}(d_\mathbf{u})\boldsymbol{\varepsilon}(\mathbf{y}_1)] d\Omega


@f}



方程3

@f[
  \int_\Omega  - d_\varrho c_{z_1} +d_\varrho c_{z_2}  + H(d_\varrho)c_{y_2}  d\Omega =


  -\int_\Omega -d_\varrho z_1 + d_\varrho z_2 + H(d_\varrho)y_2 d\Omega


@f]

 </li> 

 <li>  原始可行性。这些方程保证了平等约束的满足。

方程4

@f{align}{
  &\int_\Omega p \rho^{p-1} c_p[\lambda \nabla \cdot
  \mathbf{d}_{\mathbf{y}_1} \nabla \cdot \mathbf{u} +  \mu
  \boldsymbol{\varepsilon}(\mathbf{u}) \boldsymbol{\varepsilon}(\mathbf{d}_{\mathbf{y}_1})] +
  \rho^{p}[\lambda \nabla \cdot \mathbf{d}_{\mathbf{y}_1} \nabla \cdot
  \mathbf{c}_\mathbf{u} +  \mu  \boldsymbol{\varepsilon}(\mathbf{c}_\mathbf{u})
  \boldsymbol{\varepsilon}(\mathbf{d}_{\mathbf{y}_1})] d\Omega \\
  &= -\int_\Omega \rho^{p}[\lambda \nabla \cdot \mathbf{d}_{\mathbf{y}_1} \nabla
  \cdot \mathbf{u} + \mu  \boldsymbol{\varepsilon}(\mathbf{u}) \boldsymbol{\varepsilon}
  (\mathbf{d}_{\mathbf{y}_1})]  + \int_\Gamma  \mathbf{d}_{\mathbf{y}_1}
  \cdot \mathbf{t} d\partial\Omega


@f}



方程5

@f[


  -\int_\Omega d_{z_1}(c_\varrho - c_{s_1}) d\Omega=\int_\Omega d_{z_1} (\varrho - s_1) d\Omega


@f]



方程6

@f[


  -\int_\Omega d_{z_2}(-c_\varrho-c_{s_2}) d\Omega= \int_\Omega d_{z_2} (1-\varrho-s_2) d\Omega


@f]



方程7

@f[


  -\int_\Omega   d_{y_2}(c_\rho - H(c_\varrho)) d\Omega=\int_\Omega d_{y_2}
  (\rho - H(\varrho)) d\Omega


@f]

 </li> 

 <li>  互补松弛性。这些方程基本上确保了障碍的满足--在最终的解决方案中，我们需要  $s^T z = 0$  。

方程8

@f[
  \int_\Omega d_{s_1}(c_{s_1}z_1/s_1 +  c_{z_1} ) d\Omega=-\int_\Omega d_{s_1}
  (z_1 - \alpha/s_1) d\Omega ,\;\;\; \alpha \to 0


@f]



方程9

@f[
  \int_\Omega d_{s_2} (c_{s_2}z_2/s_2 + c_{z_2} ) d\Omega=-\int_\Omega d_{s_2}
  (z_2 - \alpha/s_2)  d\Omega,\;\;\; \alpha \to 0


@f]

 </li> 

 <li>  双重可行性。松弛和松弛变量的拉格朗日乘数必须保持大于0。（这是唯一没有在 `SANDTopOpt::assemble_system()` 函数中实现的部分）。

@f[
  s,z \geq 0


@f]

 </li>   </ol> 




<a name="Discretization"></a><h3>Discretization</h3>我们使用带有 $Q_1$ 元素的四边形网格来离散位移和位移Lagrange乘数。分片常数 $DGQ_0$ 元素被用来离散密度、未过滤密度、密度松弛变量以及松弛变量和过滤约束的乘数。


<a name="NonlinearAlgorithm"></a><h3>Nonlinear Algorithm</h3>


虽然上面的大部分讨论都是按照传统的和众所周知的方法来解决非线性优化问题，但事实证明，这个问题在实践中其实是相当难解决的。特别是，它是相当非线性的，一个重要的问题不仅仅是像上面讨论的基于牛顿方法的搜索方向 $c_{\{\bullet\}}$ ，而是人们需要花相当多的注意力在这个方向上要走多远。这通常被称为 "线搜索"，归结为如何选择步长 $\alpha_k \in (0,1]$ 的问题，以便我们以尽可能有效的方式从当前迭代 $\mathbf{x}_k$ 移动到下一个迭代 $\mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k \mathbf{x}_k$ 。众所周知，我们最终需要选择 $\alpha_k=1$ 来实现牛顿方法的二次收敛；然而，在早期迭代中，采取如此长的步长实际上可能会使事情变得更糟，要么导致一个目标函数更差的点，要么在这个点上的约束条件的满足程度不如在 $\mathbf{x}_k$ 时。

已经提出了非常复杂的算法来处理这个问题  @cite Nocedal2009   @cite Waechter2005  。在这里，我们实现了一个看门狗搜索算法  @cite Nocedal2006  。在讨论这个算法时，我们将使用向量 $\mathbf{x}$ 来表示所有的原始变量--过滤和未过滤的密度、松弛变量和位移，并使用向量 $\mathbf{y}$ 来表示所有的对偶向量。上述非线性方程组的（增量）解决方案现在将被称为 $\Delta \mathbf{x}$ 和 $\Delta
\mathbf{y}$ ，而不是 $c_{\{\bullet\}}$  。一个优点函数（后面有详细解释）在这里被称为 $\phi(\mathbf{x,\mathbf{y}})$  。

应用于具有给定障碍参数的子问题的看门狗算法以如下方式工作。首先，当前迭代被保存为 "看门狗 "状态，并记录看门狗状态的优点。然后采取一个最大的可行的牛顿步骤。如果功绩比第一步充分减少，则接受这个新步骤。如果不是，则采取另一个最大可行的牛顿步骤，并再次将功绩与看门狗的功绩进行比较。如果经过一定数量（通常在5到8之间）的牛顿步骤后，功绩没有充分减少，算法从看门狗状态或最后一次迭代中选择一个缩放的牛顿步骤，以保证功绩充分减少，该步骤被接受。一旦一个步骤被接受，就会测量KKT误差的规范，如果它足够小，就会减少障碍值。如果不够小，则将最后接受的步骤作为新的看门狗步骤，并重复这一过程。


以上，"最大可行步长 "是对牛顿步长在原始变量和对偶变量中的一个缩放，其公式为

@f[
  \beta^\mathbf{y} = \min\{1,\max \beta \text{ such that }\left(\mathbf{z}_{k+i}
   + \beta^\mathbf{z}_{k+i} \Delta \mathbf{z}_{k+i}\right)_j \geq \zeta
   \mathbf{z}_{k+i,j} \forall j\}


@f]



@f[
  \beta^\mathbf{x} = \min\{1,\max \beta \text{ such that }\left(\mathbf{s}_{k+i}
   + \beta^\mathbf{s}_{k+i} \Delta \mathbf{s}_{k+i}\right)_j \geq \zeta
   \mathbf{s}_{k+i,j} \forall j\}


@f]



以上， $\zeta$ 是任何步骤上允许的 "到边界的分数"。由于导数在边界附近变得条件不良，这种技术代表了[信任区域](https://en.wikipedia.org/wiki/Trust_region)，对于确保未来的良好近似是必要的。   $\zeta$ 被认为是 $\max\{0.8, 1-\alpha\}$ ，这允许随着障碍物变小而向边界靠近。未来，在实施减少障碍物的LOQO算法时，必须将其保持在0.8，因为障碍物参数可能变化很大。

另外，我们需要处理我们用来强制执行松弛变量的正性约束的对数障碍 $s_1,s_2$ ：在我们解决的最终优化问题的声明中，我们添加了术语

@f[


  -\alpha \int_\Omega (\log(s_1) + \log(s_2)) d\Omega.


@f]

问题是我们应该如何选择惩罚因子  $\alpha$  。与所有的惩罚方法一样，我们实际上只对极限 $\alpha\to 0$ 感兴趣，因为这才是我们真正想要解决的问题，受松弛变量的正性约束。另一方面，我们需要选择足够大的 $\alpha$ 来使问题在实践中可以解决。因此，实际的实现从较大的 $\alpha$ 值开始，并随着外迭代的进行而逐渐减小它。

在这里实现的单调方法中，每当在当前的障碍参数下达到某种程度的收敛时，就会更新障碍参数。我们使用KKT条件的 $l_\infty$ 准则来检查每个障碍大小的收敛情况。要求是 $\|KKT\|_{l_\infty} < c \cdot \alpha$ ，其中 $c$ 是任何障碍大小的常数， $\alpha$ 是障碍参数。这迫使在以后的迭代中更好地收敛，这与[IPOPT](https://coin-or.github.io/Ipopt/)（一个用于大规模非线性优化的开源软件包）中的要求相同。

在这里，障碍值在较大的数值下是线性减少的，在较小的数值下是超线性的。在较大的数值下，它被乘以一个常数（大约0.6），而在较低的数值下，障碍值被提高到某个指数（大约1.2）的障碍值所取代。事实证明，这种方法能够有效地保持大障碍值下子问题的可解性，同时在较小的障碍值下仍然允许超线性收敛。在实践中，这看起来像以下情况。

@f[
  \alpha_{k+1} = \min\{\alpha_k^{1.2},0.6\alpha_k\}


@f]



虽然在达到收敛时大步减少障碍物的大小被广泛使用，但最近的研究表明，通常使用每次迭代自适应更新障碍物的算法会更快，也就是说，我们在每次迭代结束时使用具体的标准来决定下一次迭代中的惩罚参数应该是什么，而不是使用独立于当前解决方案的减少因素。也就是说，这样的方法也比较复杂，我们在此不做介绍。

<a name="MeritFunction"></a><h3>Merit %Function</h3>


上面概述的算法利用了 "优点函数"。功绩函数用于确定从 $x_k$ 到建议点 $x_{k+1}$ 的一步是否有利。在无约束的优化问题中，人们可以简单地用我们试图最小化的目标函数来检查，通常使用[沃尔夫和戈尔茨坦条件]（https://en.wikipedia.org/wiki/Wolfe_conditions）等条件。

在有约束的优化问题中，问题是如何平衡目标函数的减少和可能增加的对约束的违反。一个建议的步骤可能会使目标函数变小，但离满足约束条件的点集更远，或者相反。这种权衡通常通过使用结合这两个标准的优点函数来解决。

在这里，我们使用一个精确的 $l_1$ 功绩函数来测试步骤。

@f{align}{
  \phi(\mathbf{x},\mathbf{y}) =& \int_{\partial \Omega} \mathbf{u}\cdot
  \mathbf{t} d\partial\Omega- \alpha \int_\Omega (\log(s_1) + \log(s_2)) + p \sum_i\left|
  \int_\Omega y_{2,i}(H(\varrho) - \rho) d\Omega \right| \\
  & + p \sum_i\left| \int_{\partial \Omega} \mathbf{y}_{1,i}\cdot \mathbf{t}  d\partial\Omega


  - \int_\Omega \rho^p[\lambda \nabla \cdot \mathbf{u} \nabla \cdot \mathbf{y}_{1,i}
  + \mu \boldsymbol{\varepsilon}{\mathbf{u}}\boldsymbol{\varepsilon}{\mathbf{y}_{1,i}}] d\Omega \right|
  + p \sum_i\left| \int_\Omega z_{1,i}(s_1 - \varrho) d\Omega\right|
  + p \sum_i\left| \int_\Omega z_{2,i}(1-\varrho - s_2) d\Omega\right|


@f}



这里， $p$ 是一个惩罚参数。这个优点函数是精确的，意味着存在一些 $p_0$ ，以便对于任何 $p > p_0$ ，优点函数的最小值与原始问题的位置相同。这个惩罚参数被更新（根据Nocedal和Wright @cite Benson2002 的建议），如下。

@f[
  p > \frac{\frac{1}{2} \mathbf{x}^T \cdot \mathbf{H} \cdot \mathbf{x} - \mathbf{x}^T \cdot \nabla f}{\|c_i\|_{l_\infty}}
  \quad , i \in \mathcal{E},


@f]

其中 $\mathbf{H}$ 是目标函数的Hessian， $\mathbf{x}$ 是我们的决策（原始）变量的矢量， $f$ 是目标函数， $c_i$ 是当前平等约束的误差。

我们使用这种方法的部分原因是在寻找右手边时已经计算了大部分必要的部分，而且使用精确的优点函数可以确保它在与整个问题相同的位置被最小化。最近的研究表明，人们可以用所谓的 "滤波方法 "代替优点函数，人们应该考虑使用这些方法，因为它们被证明是更有效的。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Preliminaries"></a> 
 * <h3>Preliminaries</h3>
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/tensor.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/signaling_nan.h> 
 * 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/linear_operator.h> 
 * #include <deal.II/lac/packaged_operation.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/fe_dgq.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_q.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * #include <algorithm> 
 * 
 * @endcode
 * 
 * 以上是相当常见的包含文件。这些文件还包括稀疏直接类的文件 SparseDirectUMFPACK。这不是解决大型线性问题的最有效的方法，但现在可以了。
 * 

 * 
 * 像往常一样，我们把所有的东西都放到一个共同的命名空间里。然后，我们开始声明一些常数的符号名称，这些常数将在本教程中使用。具体来说，我们在这个程序中有*多的变量（当然是密度和位移，但也有未过滤的密度和相当多的拉格朗日乘数）。我们很容易忘记这些变量在求解向量中的哪个位置，而且试图用数字来表示这些向量分量是一个错误的处方。相反，我们定义的静态变量可以在所有这些地方使用，而且只需初始化一次。在实践中，这将导致一些冗长的表达式，但它们更具可读性，而且不太可能出错。
 * 

 * 
 * 一个类似的问题出现在系统矩阵和向量中块的排序上。矩阵中有 $9\times 9$ 块，而且很难记住哪个是哪个。对这些块也使用符号名称要容易得多。
 * 

 * 
 * 最后，我们为我们将要使用的边界指标引入符号名称，与  step-19  中的精神相同。
 * 

 * 
 * 在所有这些情况下，我们将这些变量声明为命名空间中的成员。在求解组件的情况下，这些变量的具体数值取决于空间维度，因此我们使用[模板变量](https:en.cppreference.com/w/cpp/language/variable_template)来使变量的数值取决于模板参数，就像我们经常使用模板函数一样。
 * 

 * 
 * 
 * @code
 * namespace SAND 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 这个命名空间记录了我们的有限元系统中与每个变量相对应的第一个组件。
 * 

 * 
 * 
 * @code
 *   namespace SolutionComponents 
 *   { 
 *     template <int dim> 
 *     constexpr unsigned int density = 0; 
 *     template <int dim> 
 *     constexpr unsigned int displacement = 1; 
 *     template <int dim> 
 *     constexpr unsigned int unfiltered_density = 1 + dim; 
 *     template <int dim> 
 *     constexpr unsigned int displacement_multiplier = 2 + dim; 
 *     template <int dim> 
 *     constexpr unsigned int unfiltered_density_multiplier = 2 + 2 * dim; 
 *     template <int dim> 
 *     constexpr unsigned int density_lower_slack = 3 + 2 * dim; 
 *     template <int dim> 
 *     constexpr unsigned int density_lower_slack_multiplier = 4 + 2 * dim; 
 *     template <int dim> 
 *     constexpr unsigned int density_upper_slack = 5 + 2 * dim; 
 *     template <int dim> 
 *     constexpr unsigned int density_upper_slack_multiplier = 6 + 2 * dim; 
 *   } // namespace SolutionComponents 
 * 
 * @endcode
 * 
 * 这是一个命名空间，它记录了哪个区块对应于哪个变量。
 * 

 * 
 * 
 * @code
 *   namespace SolutionBlocks 
 *   { 
 *     constexpr unsigned int density                        = 0; 
 *     constexpr unsigned int displacement                   = 1; 
 *     constexpr unsigned int unfiltered_density             = 2; 
 *     constexpr unsigned int displacement_multiplier        = 3; 
 *     constexpr unsigned int unfiltered_density_multiplier  = 4; 
 *     constexpr unsigned int density_lower_slack            = 5; 
 *     constexpr unsigned int density_lower_slack_multiplier = 6; 
 *     constexpr unsigned int density_upper_slack            = 7; 
 *     constexpr unsigned int density_upper_slack_multiplier = 8; 
 *   } // namespace SolutionBlocks 
 * 
 *   namespace BoundaryIds 
 *   { 
 *     constexpr types::boundary_id down_force = 101; 
 *     constexpr types::boundary_id no_force   = 102; 
 *   } // namespace BoundaryIds 
 * 
 *   namespace ValueExtractors 
 *   { 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar 
 *       densities(SolutionComponents::density<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Vector 
 *       displacements(SolutionComponents::displacement<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar 
 *       unfiltered_densities(SolutionComponents::unfiltered_density<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Vector displacement_multipliers( 
 *       SolutionComponents::displacement_multiplier<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar unfiltered_density_multipliers( 
 *       SolutionComponents::unfiltered_density_multiplier<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar 
 *       density_lower_slacks(SolutionComponents::density_lower_slack<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar density_lower_slack_multipliers( 
 *       SolutionComponents::density_lower_slack_multiplier<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar 
 *       density_upper_slacks(SolutionComponents::density_upper_slack<dim>); 
 *     template <int dim> 
 *     const FEValuesExtractors::Scalar density_upper_slack_multipliers( 
 *       SolutionComponents::density_upper_slack_multiplier<dim>); 
 *   } // namespace ValueExtractors 
 * @endcode
 * 
 * 
 * <a name="TheSANDTopOptmainclass"></a> 
 * <h3>The SANDTopOpt main class</h3>
 * 

 * 
 * 接下来是这个问题的主类。大多数函数都遵循教程程序的常规命名方式，不过有几个函数因为长度问题被从通常称为`setup_system()`的函数中分离出来，还有一些函数是处理优化算法的各个方面的。
 * 

 * 
 * 作为额外的奖励，该程序将计算出的设计写成STL文件，例如，可以将其发送给3D打印机。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class SANDTopOpt 
 *   { 
 *   public: 
 *     SANDTopOpt(); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void create_triangulation(); 
 * 
 *     void setup_boundary_values(); 
 * 
 *     void setup_block_system(); 
 * 
 *     void setup_filter_matrix(); 
 * 
 *     void assemble_system(); 
 * 
 *     BlockVector<double> solve(); 
 * 
 *     std::pair<double, double> 
 *     calculate_max_step_size(const BlockVector<double> &state, 
 *                             const BlockVector<double> &step) const; 
 * 
 *     BlockVector<double> 
 *     calculate_test_rhs(const BlockVector<double> &test_solution) const; 
 * 
 *     double calculate_exact_merit(const BlockVector<double> &test_solution); 
 * 
 *     BlockVector<double> find_max_step(); 
 * 
 *     BlockVector<double> compute_scaled_step(const BlockVector<double> &state, 
 *                                             const BlockVector<double> &step, 
 *                                             const double descent_requirement); 
 * 
 *     bool check_convergence(const BlockVector<double> &state); 
 * 
 *     void output_results(const unsigned int j) const; 
 * 
 *     void write_as_stl(); 
 * 
 *     std::set<typename Triangulation<dim>::cell_iterator> 
 *     find_relevant_neighbors( 
 *       typename Triangulation<dim>::cell_iterator cell) const; 
 * 
 * @endcode
 * 
 * 大部分的成员变量也是标准的。但是，有一些变量是专门与优化算法有关的（比如下面的各种标量因子），以及过滤器矩阵，以确保设计保持平稳。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim>        triangulation; 
 *     FESystem<dim>             fe; 
 *     DoFHandler<dim>           dof_handler; 
 *     AffineConstraints<double> constraints; 
 * 
 *     std::map<types::global_dof_index, double> boundary_values; 
 * 
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> system_matrix; 
 * 
 *     SparsityPattern      filter_sparsity_pattern; 
 *     SparseMatrix<double> filter_matrix; 
 * 
 *     BlockVector<double> system_rhs; 
 *     BlockVector<double> nonlinear_solution; 
 * 
 *     const double density_ratio; 
 *     const double density_penalty_exponent; 
 *     const double filter_r; 
 *     double       penalty_multiplier; 
 *     double       barrier_size; 
 * 
 *     TimerOutput timer; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Constructorandsetupfunctions"></a> 
 * <h3>Constructor and set-up functions</h3>
 * 

 * 
 * 我们初始化一个由2  $\times$  dim `FE_Q(1)`元素组成的FES系统，用于位移变量及其拉格朗日乘数，以及7 `FE_DGQ(0)`元素。 这些片状常数函数用于与密度相关的变量：密度本身、未过滤的密度、用于未过滤的密度的下限和上限的松弛变量，然后是用于过滤和未过滤的密度之间的连接以及不等式约束的拉格朗日乘子。
 * 

 * 
 * 这些元素出现的顺序在上面有记载。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   SANDTopOpt<dim>::SANDTopOpt() 
 *     : fe(FE_DGQ<dim>(0), 
 *          1, 
 *          (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
 *          1, 
 *          FE_DGQ<dim>(0), 
 *          1, 
 *          (FESystem<dim>(FE_Q<dim>(1) ^ dim)), 
 *          1, 
 *          FE_DGQ<dim>(0), 
 *          5) 
 *     , dof_handler(triangulation) 
 *     , density_ratio(.5) 
 *     , density_penalty_exponent(3) 
 *     , filter_r(.251) 
 *     , penalty_multiplier(1) 
 *     , timer(std::cout, TimerOutput::summary, TimerOutput::wall_times) 
 *   { 
 *     Assert(dim > 1, ExcNotImplemented()); 
 *   } 
 * 
 * @endcode
 * 
 * 然后，第一步是创建与介绍中的问题描述相匹配的三角形--一个6乘1的矩形（或者一个6乘1乘1的3D盒子），在这个盒子的顶部中心将施加一个力。然后，这个三角形被均匀地细化若干次。
 * 

 * 
 * 与本程序的其他部分相比，这个函数特别假定我们是在2D中，如果我们想转到3D模拟，就需要进行修改。我们通过函数顶部的断言来确保没有人试图不经修改就意外地在三维中运行。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::create_triangulation() 
 *   { 
 *     Assert(dim == 2, ExcNotImplemented()); 
 *     GridGenerator::subdivided_hyper_rectangle(triangulation, 
 *                                               {6, 1}, 
 *                                               Point<dim>(0, 0), 
 *                                               Point<dim>(6, 1)); 
 * 
 *     triangulation.refine_global(3); 
 * 
 * @endcode
 * 
 * 第二步是将边界指标应用于边界的一部分。下面的代码分别为盒子的底部、顶部、左侧和右侧的边界分配了边界指示器。顶部边界的中心区域被赋予一个单独的边界指示器。这就是我们要施加向下力的地方。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       { 
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             if (face->at_boundary()) 
 *               { 
 *                 const auto center = face->center(); 
 *                 if (std::fabs(center(1) - 1) < 1e-12) 
 *                   { 
 *                     if ((std::fabs(center(0) - 3) < .3)) 
 *                       face->set_boundary_id(BoundaryIds::down_force); 
 *                     else 
 *                       face->set_boundary_id(BoundaryIds::no_force); 
 *                   } 
 *                 else 
 *                   face->set_boundary_id(BoundaryIds::no_force); 
 *               } 
 *           } 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 接下来，确定由于边界值而产生的约束。 域的底角在 $y$ 方向保持不变--左下角也在 $x$ 方向。deal.II通常认为边界值是附着在边界的片段上的，即面，而不是单个顶点。的确，从数学上讲，对于无穷大的偏微分方程，我们不能把边界值分配给单个点。但是，由于我们试图重现一个广泛使用的基准，我们还是要这样做，并牢记我们有一个有限维的问题，在单个节点上施加边界条件是有效的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::setup_boundary_values() 
 *   { 
 *     boundary_values.clear(); 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             if (face->at_boundary()) 
 *               { 
 *                 const auto center = face->center(); 
 * 
 * @endcode
 * 
 * 检查当前面是否在底层边界上，如果是，则检查其顶点之一是否可能是左底层或右底层顶点。
 * 

 * 
 * 
 * @code
 *                 if (std::fabs(center(1) - 0) < 1e-12) 
 *                   { 
 *                     for (const auto vertex_number : cell->vertex_indices()) 
 *                       { 
 *                         const auto vert = cell->vertex(vertex_number); 
 * 
 *                         if (std::fabs(vert(0) - 0) < 1e-12 && 
 *                             std::fabs(vert(1) - 0) < 1e-12) 
 *                           { 
 *                             types::global_dof_index x_displacement = 
 *                               cell->vertex_dof_index(vertex_number, 0); 
 *                             types::global_dof_index y_displacement = 
 *                               cell->vertex_dof_index(vertex_number, 1); 
 *                             types::global_dof_index x_displacement_multiplier = 
 *                               cell->vertex_dof_index(vertex_number, 2); 
 *                             types::global_dof_index y_displacement_multiplier = 
 *                               cell->vertex_dof_index(vertex_number, 3); 
 * 
 *                             boundary_values[x_displacement]            = 0; 
 *                             boundary_values[y_displacement]            = 0; 
 *                             boundary_values[x_displacement_multiplier] = 0; 
 *                             boundary_values[y_displacement_multiplier] = 0; 
 *                           } 
 * 
 *                         else if (std::fabs(vert(0) - 6) < 1e-12 && 
 *                                  std::fabs(vert(1) - 0) < 1e-12) 
 *                           { 
 *                             types::global_dof_index y_displacement = 
 *                               cell->vertex_dof_index(vertex_number, 1); 
 *                             types::global_dof_index y_displacement_multiplier = 
 *                               cell->vertex_dof_index(vertex_number, 3); 
 * 
 *                             boundary_values[y_displacement]            = 0; 
 *                             boundary_values[y_displacement_multiplier] = 0; 
 *                           } 
 *                       } 
 *                   } 
 *               } 
 *           } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Settingupblockmatricesandvectors"></a> 
 * <h3>Setting up block matrices and vectors</h3>
 * 

 * 
 * 下一个函数制作了一个巨大的9乘9的块状矩阵，并且还设置了必要的块状向量。 这个矩阵的稀疏度模式包括滤波矩阵的稀疏度模式。它还初始化了我们将使用的任何块向量。
 * 

 * 
 * 设置块本身并不复杂，并且遵循诸如  step-22  等程序中已经完成的工作，例如。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::setup_block_system() 
 *   { 
 *     std::vector<unsigned int> block_component(9, 2); 
 *     block_component[0] = 0; 
 *     block_component[1] = 1; 
 *     const std::vector<types::global_dof_index> dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 * 
 *     const types::global_dof_index                     n_p = dofs_per_block[0]; 
 *     const types::global_dof_index                     n_u = dofs_per_block[1]; 
 *     const std::vector<BlockVector<double>::size_type> block_sizes = { 
 *       n_p, n_u, n_p, n_u, n_p, n_p, n_p, n_p, n_p}; 
 * 
 *     BlockDynamicSparsityPattern dsp(9, 9); 
 *     for (unsigned int k = 0; k < 9; ++k) 
 *       for (unsigned int j = 0; j < 9; ++j) 
 *         dsp.block(j, k).reinit(block_sizes[j], block_sizes[k]); 
 *     dsp.collect_sizes(); 
 * 
 * @endcode
 * 
 * 该函数的大部分内容是设置这些块中哪些将实际包含任何内容，即哪些变量与哪些其他变量相耦合。这很麻烦，但也是必要的，以确保我们不会为我们的矩阵分配大量的条目，而这些条目最终会变成零。
 * 

 * 
 * 你在下面看到的具体模式可能需要在纸上画一次，但是从我们在每次非线性迭代中必须组装的双线性形式的许多项来看，它是相对直接的方式。
 * 

 * 
 * 使用命名空间 "SolutionComponents "中定义的符号名称有助于理解下面每个项所对应的内容，但它也使表达式变得冗长而不流畅。像 `coupling[SolutionComponents::density_upper_slack_multiplier<dim>][SolutionComponents::density<dim>]` 这样的术语读起来就不太顺口，要么必须分成几行，要么几乎跑到每个屏幕的右边缘。因此，我们打开了一个大括号封闭的代码块，在这个代码块中，我们通过说 "使用命名空间SolutionComponents"，暂时使命名空间`SolutionComponents'中的名字可用，而不需要命名空间修饰语。
 * 

 * 
 * 
 * @code
 *     Table<2, DoFTools::Coupling> coupling(2 * dim + 7, 2 * dim + 7); 
 *     { 
 *       using namespace SolutionComponents; 
 * 
 *       coupling[density<dim>][density<dim>] = DoFTools::always; 
 * 
 *       for (unsigned int i = 0; i < dim; ++i) 
 *         { 
 *           coupling[density<dim>][displacement<dim> + i] = DoFTools::always; 
 *           coupling[displacement<dim> + i][density<dim>] = DoFTools::always; 
 *         } 
 * 
 *       for (unsigned int i = 0; i < dim; ++i) 
 *         { 
 *           coupling[density<dim>][displacement_multiplier<dim> + i] = 
 *             DoFTools::always; 
 *           coupling[displacement_multiplier<dim> + i][density<dim>] = 
 *             DoFTools::always; 
 *         } 
 * 
 *       coupling[density<dim>][unfiltered_density_multiplier<dim>] = 
 *         DoFTools::always; 
 *       coupling[unfiltered_density_multiplier<dim>][density<dim>] = 
 *         DoFTools::always; 
 *       /*位移的联结  */ 
 *       for (unsigned int i = 0; i < dim; ++i) 
 *         { 
 *           for (unsigned int k = 0; k < dim; ++k) 
 *             { 
 *               coupling[displacement<dim> + i] 
 *                       [displacement_multiplier<dim> + k] = DoFTools::always; 
 *               coupling[displacement_multiplier<dim> + k] 
 *                       [displacement<dim> + i] = DoFTools::always; 
 *             } 
 *         } 
 *       /*松弛变量的耦合 */ 
 *       coupling[density_lower_slack<dim>][density_lower_slack<dim>] = 
 *         DoFTools::always; 
 *       coupling[density_lower_slack<dim>][density_upper_slack<dim>] = 
 *         DoFTools::always; 
 *       coupling[density_upper_slack<dim>][density_lower_slack<dim>] = 
 *         DoFTools::always; 
 * 
 *       coupling[density_lower_slack_multiplier<dim>] 
 *               [density_lower_slack_multiplier<dim>] = DoFTools::always; 
 *       coupling[density_lower_slack_multiplier<dim>] 
 *               [density_upper_slack_multiplier<dim>] = DoFTools::always; 
 *       coupling[density_upper_slack_multiplier<dim>] 
 *               [density_lower_slack_multiplier<dim>] = DoFTools::always; 
 *     } 
 * 
 * @endcode
 * 
 * 在创建稀疏模式之前，我们还必须设置约束。由于这个程序没有自适应地细化网格，我们唯一的约束是将所有的密度变量耦合在一起，强制执行体积约束。这将最终导致矩阵的密集子块，但我们对此没有什么办法。
 * 

 * 
 * 
 * @code
 *     const ComponentMask density_mask = 
 *       fe.component_mask(ValueExtractors::densities<dim>); 
 *     const IndexSet density_dofs = 
 *       DoFTools::extract_dofs(dof_handler, density_mask); 
 * 
 *     types::global_dof_index last_density_dof = 
 *       density_dofs.nth_index_in_set(density_dofs.n_elements() - 1); 
 *     constraints.clear(); 
 *     constraints.add_line(last_density_dof); 
 *     for (unsigned int i = 0; i < density_dofs.n_elements() - 1; ++i) 
 *       constraints.add_entry(last_density_dof, 
 *                             density_dofs.nth_index_in_set(i), 
 *                             -1); 
 *     constraints.set_inhomogeneity(last_density_dof, 0); 
 * 
 *     constraints.close(); 
 * 
 * @endcode
 * 
 * 现在我们终于可以为矩阵创建稀疏模式了，考虑到哪些变量与哪些其他变量耦合，以及我们对密度的约束。
 * 

 * 
 * 
 * @code
 *     DoFTools::make_sparsity_pattern(dof_handler, coupling, dsp, constraints); 
 * 
 * @endcode
 * 
 * 矩阵中唯一没有处理的部分是过滤矩阵和它的转置。这些都是非局部（积分）运算符，目前deal.II还没有相关的函数。我们最终需要做的是遍历所有单元，并将此单元上的未过滤密度与小于阈值距离的相邻单元的所有过滤密度联系起来，反之亦然；目前，我们只关心建立与这种矩阵相对应的稀疏模式，所以我们执行等效循环，以后我们将写进矩阵的一个条目，现在我们只需向稀疏矩阵添加一个条目。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         const unsigned int i = cell->active_cell_index(); 
 *         for (const auto &check_cell : find_relevant_neighbors(cell)) 
 *           { 
 *             const double distance = 
 *               cell->center().distance(check_cell->center()); 
 *             if (distance < filter_r) 
 *               { 
 *                 dsp 
 *                   .block(SolutionBlocks::unfiltered_density, 
 *                          SolutionBlocks::unfiltered_density_multiplier) 
 *                   .add(i, check_cell->active_cell_index()); 
 *                 dsp 
 *                   .block(SolutionBlocks::unfiltered_density_multiplier, 
 *                          SolutionBlocks::unfiltered_density) 
 *                   .add(i, check_cell->active_cell_index()); 
 *               } 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 在生成了 "动态 "稀疏度模式之后，我们终于可以将其复制到用于将矩阵与稀疏度模式联系起来的结构中。由于稀疏模式很大很复杂，我们还将其输出到一个自己的文件中，以达到可视化的目的--换句话说，是为了 "可视化调试"。
 * 

 * 
 * 
 * @code
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     std::ofstream out("sparsity.plt"); 
 *     sparsity_pattern.print_gnuplot(out); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 剩下的就是正确确定各种向量及其块的大小，以及为（非线性）解向量的一些分量设置初始猜测。我们在这里使用解向量各个区块的符号分量名称，为了简洁起见，使用与上面的 "使用命名空间 "相同的技巧。
 * 

 * 
 * 
 * @code
 *     nonlinear_solution.reinit(block_sizes); 
 *     system_rhs.reinit(block_sizes); 
 * 
 *     { 
 *       using namespace SolutionBlocks; 
 *       nonlinear_solution.block(density).add(density_ratio); 
 *       nonlinear_solution.block(unfiltered_density).add(density_ratio); 
 *       nonlinear_solution.block(unfiltered_density_multiplier) 
 *         .add(density_ratio); 
 *       nonlinear_solution.block(density_lower_slack).add(density_ratio); 
 *       nonlinear_solution.block(density_lower_slack_multiplier).add(50); 
 *       nonlinear_solution.block(density_upper_slack).add(1 - density_ratio); 
 *       nonlinear_solution.block(density_upper_slack_multiplier).add(50); 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Creatingthefiltermatrix"></a> 
 * <h3>Creating the filter matrix</h3>
 * 

 * 
 * 接下来是一个在程序开始时使用一次的函数。它创建了一个矩阵 $H$ ，使过滤后的密度向量等于 $H$ 乘以未过滤的密度。 这个矩阵的创建是非同小可的，它在每次迭代中都会被使用，因此，与其像我们对牛顿矩阵那样对其进行改造，不如只做一次并单独存储。
 * 

 * 
 * 这个矩阵的计算方式遵循上面已经使用过的大纲，以形成其稀疏模式。我们在这里对这个单独形成的矩阵的稀疏性模式重复这个过程，然后实际建立矩阵本身。你可能想看看本程序介绍中关于这个矩阵的定义。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::setup_filter_matrix() 
 *   { 
 * 
 * @endcode
 * 
 * 滤波器的稀疏模式已经在setup_system()函数中确定并实现。我们从相应的块中复制该结构，并在这里再次使用它。
 * 

 * 
 * 
 * @code
 *     filter_sparsity_pattern.copy_from( 
 *       sparsity_pattern.block(SolutionBlocks::unfiltered_density, 
 *                              SolutionBlocks::unfiltered_density_multiplier)); 
 *     filter_matrix.reinit(filter_sparsity_pattern); 
 * 
 * @endcode
 * 
 * 在建立了稀疏模式之后，现在我们重新做所有这些循环，以实际计算矩阵项的必要值。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         const unsigned int i = cell->active_cell_index(); 
 *         for (const auto &check_cell : find_relevant_neighbors(cell)) 
 *           { 
 *             const double distance = 
 *               cell->center().distance(check_cell->center()); 
 *             if (distance < filter_r) 
 *               { 
 *                 filter_matrix.add(i, 
 *                                   check_cell->active_cell_index(), 
 *                                   filter_r - distance); 
 * 
 * @endcode
 * 
 * 
 * 
 * 

 * 
 * 
 * @code
 *               } 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 最后一步是对矩阵进行标准化处理，使每一行的条目之和等于1。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < filter_matrix.m(); ++i) 
 *       { 
 *         double denominator = 0; 
 *         for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
 *              iter != filter_matrix.end(i); 
 *              iter++) 
 *           denominator = denominator + iter->value(); 
 *         for (SparseMatrix<double>::iterator iter = filter_matrix.begin(i); 
 *              iter != filter_matrix.end(i); 
 *              iter++) 
 *           iter->value() = iter->value() / denominator; 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 这个函数用于建立过滤矩阵。我们创建一个输入单元的一定半径内的所有单元迭代器的集合。这些是与过滤器有关的邻近单元。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::set<typename Triangulation<dim>::cell_iterator> 
 *   SANDTopOpt<dim>::find_relevant_neighbors( 
 *     typename Triangulation<dim>::cell_iterator cell) const 
 *   { 
 *     std::set<unsigned int>                               neighbor_ids; 
 *     std::set<typename Triangulation<dim>::cell_iterator> cells_to_check; 
 * 
 *     neighbor_ids.insert(cell->active_cell_index()); 
 *     cells_to_check.insert(cell); 
 * 
 *     bool new_neighbors_found; 
 *     do 
 *       { 
 *         new_neighbors_found = false; 
 *         for (const auto &check_cell : 
 *              std::vector<typename Triangulation<dim>::cell_iterator>( 
 *                cells_to_check.begin(), cells_to_check.end())) 
 *           { 
 *             for (const auto n : check_cell->face_indices()) 
 *               { 
 *                 if (!(check_cell->face(n)->at_boundary())) 
 *                   { 
 *                     const auto & neighbor = check_cell->neighbor(n); 
 *                     const double distance = 
 *                       cell->center().distance(neighbor->center()); 
 *                     if ((distance < filter_r) && 
 *                         !(neighbor_ids.count(neighbor->active_cell_index()))) 
 *                       { 
 *                         cells_to_check.insert(neighbor); 
 *                         neighbor_ids.insert(neighbor->active_cell_index()); 
 *                         new_neighbors_found = true; 
 *                       } 
 *                   } 
 *               } 
 *           } 
 *       } 
 *     while (new_neighbors_found); 
 *     return cells_to_check; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="AssemblingtheNewtonmatrix"></a> 
 * <h3>Assembling the Newton matrix</h3>
 * 

 * 
 * setup_filter_matrix函数建立了一个只要网格不改变就不变的矩阵（在这个程序中我们反正不改变），而下一个函数建立了每次迭代都要解决的矩阵。这就是奇迹发生的地方。描述牛顿求解KKT条件的方法的线性方程组的组成部分在这里实现。
 * 

 * 
 * 这个函数的顶部与大多数此类函数一样，只是设置了实际装配所需的各种变量，包括一大堆提取器。如果你以前看过  step-22  ，整个设置应该看起来很熟悉，尽管有些冗长。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::assemble_system() 
 *   { 
 *     TimerOutput::Scope t(timer, "assembly"); 
 * 
 *     system_matrix = 0; 
 *     system_rhs    = 0; 
 * 
 *     MappingQGeneric<dim> mapping(1); 
 *     QGauss<dim>          quadrature_formula(fe.degree + 1); 
 *     QGauss<dim - 1>      face_quadrature_formula(fe.degree + 1); 
 *     FEValues<dim>        fe_values(mapping, 
 *                             fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 *     FEFaceValues<dim>    fe_face_values(mapping, 
 *                                      fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_quadrature_points | 
 *                                        update_normal_vectors | 
 *                                        update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.dofs_per_cell; 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     dummy_cell_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     std::vector<double>                    lambda_values(n_q_points); 
 *     std::vector<double>                    mu_values(n_q_points); 
 *     const Functions::ConstantFunction<dim> lambda(1.); 
 *     const Functions::ConstantFunction<dim> mu(1.); 
 *     std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 
 * 
 * @endcode
 * 
 * 在这一点上，我们对未过滤的密度进行过滤，并对未过滤的密度乘法器进行邻接（转置）操作，都是对当前非线性解决方案的最佳猜测。后来我们用它来告诉我们，我们过滤的密度与应用于未过滤密度的过滤器有多大的偏差。这是因为在非线性问题的解中，我们有 $\rho=H\varrho$ ，但在中间迭代中，我们一般有 $\rho^k\neq H\varrho^k$ ，然后 "残差" $\rho^k-H\varrho^k$ 将出现在我们下面计算的牛顿更新方程中的右边。
 * 

 * 
 * 
 * @code
 *     BlockVector<double> filtered_unfiltered_density_solution = 
 *       nonlinear_solution; 
 *     BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
 *       nonlinear_solution; 
 * 
 *     filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
 *                           SolutionBlocks::unfiltered_density), 
 *                         nonlinear_solution.block( 
 *                           SolutionBlocks::unfiltered_density)); 
 *     filter_matrix.Tvmult( 
 *       filter_adjoint_unfiltered_density_multiplier_solution.block( 
 *         SolutionBlocks::unfiltered_density_multiplier), 
 *       nonlinear_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 
 * 
 *     std::vector<double>                  old_density_values(n_q_points); 
 *     std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
 *     std::vector<double>                  old_displacement_divs(n_q_points); 
 *     std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
 *     std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
 *     std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
 *       n_q_points); 
 *     std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
 *     std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
 *     std::vector<double> old_lower_slack_values(n_q_points); 
 *     std::vector<double> old_upper_slack_values(n_q_points); 
 *     std::vector<double> old_unfiltered_density_values(n_q_points); 
 *     std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
 *     std::vector<double> filtered_unfiltered_density_values(n_q_points); 
 *     std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
 *       n_q_points); 
 * 
 *     using namespace ValueExtractors; 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_matrix = 0; 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         fe_values.reinit(cell); 
 * 
 *         lambda.value_list(fe_values.get_quadrature_points(), lambda_values); 
 *         mu.value_list(fe_values.get_quadrature_points(), mu_values); 
 * 
 * @endcode
 * 
 * 作为构建系统矩阵的一部分，我们需要从我们目前对解决方案的猜测中获取数值。以下几行代码将检索出所需的值。
 * 

 * 
 * 
 * @code
 *         fe_values[densities<dim>].get_function_values(nonlinear_solution, 
 *                                                       old_density_values); 
 *         fe_values[displacements<dim>].get_function_values( 
 *           nonlinear_solution, old_displacement_values); 
 *         fe_values[displacements<dim>].get_function_divergences( 
 *           nonlinear_solution, old_displacement_divs); 
 *         fe_values[displacements<dim>].get_function_symmetric_gradients( 
 *           nonlinear_solution, old_displacement_symmgrads); 
 *         fe_values[displacement_multipliers<dim>].get_function_values( 
 *           nonlinear_solution, old_displacement_multiplier_values); 
 *         fe_values[displacement_multipliers<dim>].get_function_divergences( 
 *           nonlinear_solution, old_displacement_multiplier_divs); 
 *         fe_values[displacement_multipliers<dim>] 
 *           .get_function_symmetric_gradients( 
 *             nonlinear_solution, old_displacement_multiplier_symmgrads); 
 *         fe_values[density_lower_slacks<dim>].get_function_values( 
 *           nonlinear_solution, old_lower_slack_values); 
 *         fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
 *           nonlinear_solution, old_lower_slack_multiplier_values); 
 *         fe_values[density_upper_slacks<dim>].get_function_values( 
 *           nonlinear_solution, old_upper_slack_values); 
 *         fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
 *           nonlinear_solution, old_upper_slack_multiplier_values); 
 *         fe_values[unfiltered_densities<dim>].get_function_values( 
 *           nonlinear_solution, old_unfiltered_density_values); 
 *         fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
 *           nonlinear_solution, old_unfiltered_density_multiplier_values); 
 *         fe_values[unfiltered_densities<dim>].get_function_values( 
 *           filtered_unfiltered_density_solution, 
 *           filtered_unfiltered_density_values); 
 *         fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
 *           filter_adjoint_unfiltered_density_multiplier_solution, 
 *           filter_adjoint_unfiltered_density_multiplier_values); 
 * 
 *         for (const auto q_point : fe_values.quadrature_point_indices()) 
 *           { 
 * 
 * @endcode
 * 
 * 我们还需要几个与来自拉格朗日的第一导数的测试函数相对应的数值，也就是 $d_{\bullet}$ 函数。这些都是在这里计算的。
 * 

 * 
 * 
 * @code
 *             for (const auto i : fe_values.dof_indices()) 
 *               { 
 *                 const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
 *                   fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
 *                 const double displacement_phi_i_div = 
 *                   fe_values[displacements<dim>].divergence(i, q_point); 
 * 
 *                 const SymmetricTensor<2, dim> 
 *                   displacement_multiplier_phi_i_symmgrad = 
 *                     fe_values[displacement_multipliers<dim>].symmetric_gradient( 
 *                       i, q_point); 
 *                 const double displacement_multiplier_phi_i_div = 
 *                   fe_values[displacement_multipliers<dim>].divergence(i, 
 *                                                                       q_point); 
 * 
 *                 const double density_phi_i = 
 *                   fe_values[densities<dim>].value(i, q_point); 
 *                 const double unfiltered_density_phi_i = 
 *                   fe_values[unfiltered_densities<dim>].value(i, q_point); 
 *                 const double unfiltered_density_multiplier_phi_i = 
 *                   fe_values[unfiltered_density_multipliers<dim>].value(i, 
 *                                                                        q_point); 
 * 
 *                 const double lower_slack_multiplier_phi_i = 
 *                   fe_values[density_lower_slack_multipliers<dim>].value( 
 *                     i, q_point); 
 * 
 *                 const double lower_slack_phi_i = 
 *                   fe_values[density_lower_slacks<dim>].value(i, q_point); 
 * 
 *                 const double upper_slack_phi_i = 
 *                   fe_values[density_upper_slacks<dim>].value(i, q_point); 
 * 
 *                 const double upper_slack_multiplier_phi_i = 
 *                   fe_values[density_upper_slack_multipliers<dim>].value( 
 *                     i, q_point); 
 * 
 *                 for (const auto j : fe_values.dof_indices()) 
 *                   { 
 * 
 * @endcode
 * 
 * 最后，我们需要来自拉格朗日的第二轮导数的数值，即 $c_{\bullet}$ 函数。这些是在这里计算的。
 * 

 * 
 * 
 * @code
 *                     const SymmetricTensor<2, dim> displacement_phi_j_symmgrad = 
 *                       fe_values[displacements<dim>].symmetric_gradient(j, 
 *                                                                        q_point); 
 *                     const double displacement_phi_j_div = 
 *                       fe_values[displacements<dim>].divergence(j, q_point); 
 * 
 *                     const SymmetricTensor<2, dim> 
 *                       displacement_multiplier_phi_j_symmgrad = 
 *                         fe_values[displacement_multipliers<dim>] 
 *                           .symmetric_gradient(j, q_point); 
 *                     const double displacement_multiplier_phi_j_div = 
 *                       fe_values[displacement_multipliers<dim>].divergence( 
 *                         j, q_point); 
 * 
 *                     const double density_phi_j = 
 *                       fe_values[densities<dim>].value(j, q_point); 
 * 
 *                     const double unfiltered_density_phi_j = 
 *                       fe_values[unfiltered_densities<dim>].value(j, q_point); 
 *                     const double unfiltered_density_multiplier_phi_j = 
 *                       fe_values[unfiltered_density_multipliers<dim>].value( 
 *                         j, q_point); 
 * 
 *                     const double lower_slack_phi_j = 
 *                       fe_values[density_lower_slacks<dim>].value(j, q_point); 
 * 
 *                     const double upper_slack_phi_j = 
 *                       fe_values[density_upper_slacks<dim>].value(j, q_point); 
 * 
 *                     const double lower_slack_multiplier_phi_j = 
 *                       fe_values[density_lower_slack_multipliers<dim>].value( 
 *                         j, q_point); 
 * 
 *                     const double upper_slack_multiplier_phi_j = 
 *                       fe_values[density_upper_slack_multipliers<dim>].value( 
 *                         j, q_point); 
 * 
 * @endcode
 * 
 * 这就是实际工作的开始。在下文中，我们将建立矩阵的所有项--它们数量众多，而且不完全是不言自明的，也取决于之前的解和它的导数（我们已经在上面评估了这些导数，并将其放入名为`old_*`的变量中）。为了理解这些条款的每一个对应的内容，你要看一下上面介绍中这些条款的明确形式。                    被驱动到0的方程的右边给出了寻找局部最小值的所有KKT条件--每个单独方程的描述都是随着右边的计算给出的。
 * 

 * 
 * 
 * @code
 *                     /* 方程1  */ 
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       ( 
 *                         -density_phi_i * unfiltered_density_multiplier_phi_j 
 *                         + density_penalty_exponent * 
 *                             (density_penalty_exponent - 1) * 
 *                             std::pow(old_density_values[q_point], 
 *                                      density_penalty_exponent - 2) * 
 *                             density_phi_i * density_phi_j * 
 *                             (old_displacement_multiplier_divs[q_point] * 
 *                                old_displacement_divs[q_point] * 
 *                                lambda_values[q_point] + 
 *                              2 * mu_values[q_point] * 
 *                                (old_displacement_symmgrads[q_point] * 
 *                                 old_displacement_multiplier_symmgrads[q_point])) 
 *                         + density_penalty_exponent * 
 *                             std::pow(old_density_values[q_point], 
 *                                      density_penalty_exponent - 1) * 
 *                             density_phi_i * 
 *                             (displacement_multiplier_phi_j_div * 
 *                                old_displacement_divs[q_point] * 
 *                                lambda_values[q_point] + 
 *                              2 * mu_values[q_point] * 
 *                                (old_displacement_symmgrads[q_point] * 
 *                                 displacement_multiplier_phi_j_symmgrad)) 
 *                         + density_penalty_exponent * 
 *                             std::pow(old_density_values[q_point], 
 *                                      density_penalty_exponent - 1) * 
 *                             density_phi_i * 
 *                             (displacement_phi_j_div * 
 *                                old_displacement_multiplier_divs[q_point] * 
 *                                lambda_values[q_point] + 
 *                              2 * mu_values[q_point] * 
 *                                (old_displacement_multiplier_symmgrads[q_point] * 
 *                                 displacement_phi_j_symmgrad))); 
 *                    
 *                     /* 方程2  */ 
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       (density_penalty_exponent * 
 *                          std::pow(old_density_values[q_point], 
 *                                   density_penalty_exponent - 1) * 
 *                          density_phi_j * 
 *                          (old_displacement_multiplier_divs[q_point] * 
 *                             displacement_phi_i_div * lambda_values[q_point] + 
 *                           2 * mu_values[q_point] * 
 *                             (old_displacement_multiplier_symmgrads[q_point] * 
 *                              displacement_phi_i_symmgrad)) 
 *                        + std::pow(old_density_values[q_point], 
 *                                   density_penalty_exponent) * 
 *                            (displacement_multiplier_phi_j_div * 
 *                               displacement_phi_i_div * lambda_values[q_point] + 
 *                             2 * mu_values[q_point] * 
 *                               (displacement_multiplier_phi_j_symmgrad * 
 *                                displacement_phi_i_symmgrad)) 
 *                       ); 
 * 
 *                    /*方程3，这与过滤器有关 */ 
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       (-1 * unfiltered_density_phi_i * 
 *                          lower_slack_multiplier_phi_j + 
 *                        unfiltered_density_phi_i * upper_slack_multiplier_phi_j); 
 * 
 *                      /* 方程4：原始可行性  */ 
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       ( 
 *                         density_penalty_exponent * 
 *                           std::pow(old_density_values[q_point], 
 *                                    density_penalty_exponent - 1) * 
 *                           density_phi_j * 
 *                           (old_displacement_divs[q_point] * 
 *                              displacement_multiplier_phi_i_div * 
 *                              lambda_values[q_point] + 
 *                            2 * mu_values[q_point] * 
 *                              (old_displacement_symmgrads[q_point] * 
 *                               displacement_multiplier_phi_i_symmgrad)) 
 * 
 *                         + std::pow(old_density_values[q_point], 
 *                                    density_penalty_exponent) * 
 *                             (displacement_phi_j_div * 
 *                                displacement_multiplier_phi_i_div * 
 *                                lambda_values[q_point] + 
 *                              2 * mu_values[q_point] * 
 *                                (displacement_phi_j_symmgrad * 
 *                                 displacement_multiplier_phi_i_symmgrad))); 
 * 
 *                    /*等式5：原始可行性  */ 
 *                     cell_matrix(i, j) += 
 *                       -1 * fe_values.JxW(q_point) * 
 *                       lower_slack_multiplier_phi_i * 
 *                       (unfiltered_density_phi_j - lower_slack_phi_j); 
 *                   /* 等式6：原始可行性  */ 
 *                     cell_matrix(i, j) += 
 *                       -1 * fe_values.JxW(q_point) * 
 *                       upper_slack_multiplier_phi_i * 
 *                       (-1 * unfiltered_density_phi_j - upper_slack_phi_j); 
 *                     /* Equation 7: Primal feasibility - the part with the filter
 *                      * is added later */
 *                     cell_matrix(i, j) += -1 * fe_values.JxW(q_point) * 
 *                                          unfiltered_density_multiplier_phi_i * 
 *                                          (density_phi_j); 
 *                     /* Equation 8: Complementary slackness */
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       (lower_slack_phi_i * lower_slack_multiplier_phi_j 
 * 
 *                        + lower_slack_phi_i * lower_slack_phi_j * 
 *                            old_lower_slack_multiplier_values[q_point] / 
 *                            old_lower_slack_values[q_point]); 
 *                     /* Equation 9: Complementary slackness */
 *                     cell_matrix(i, j) += 
 *                       fe_values.JxW(q_point) * 
 *                       (upper_slack_phi_i * upper_slack_multiplier_phi_j 
 * 
 *                        + upper_slack_phi_i * upper_slack_phi_j * 
 *                            old_upper_slack_multiplier_values[q_point] / 
 *                            old_upper_slack_values[q_point]); 
 *                   } 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 现在我们已经把所有的东西都组装好了，我们要做的就是处理（Dirichlet）边界条件的影响和其他约束。我们将前者与当前单元的贡献结合在一起，然后让AffineConstraint类来处理后者，同时将当前单元的贡献复制到全局线性系统中。
 * 

 * 
 * 
 * @code
 *         MatrixTools::local_apply_boundary_values(boundary_values, 
 *                                                  local_dof_indices, 
 *                                                  cell_matrix, 
 *                                                  dummy_cell_rhs, 
 *                                                  true); 
 * 
 *         constraints.distribute_local_to_global(cell_matrix, 
 *                                                local_dof_indices, 
 *                                                system_matrix); 
 *       } 
 * 
 * @endcode
 * 
 * 在积累了所有属于牛顿矩阵的项之后，我们现在还必须计算右手边的项（即负残差）。我们已经在另一个函数中做了这个工作，所以我们在这里调用它。
 * 

 * 
 * 
 * @code
 *     system_rhs = calculate_test_rhs(nonlinear_solution); 
 * 
 * @endcode
 * 
 * 这里我们使用我们已经构建好的过滤器矩阵。我们只需要整合这个应用于测试函数的过滤器，它是片状常数，所以整合变成了简单的乘以单元格的度量。 遍历预制的过滤器矩阵可以让我们使用哪些单元格在过滤器中或不在过滤器中的信息，而不需要再次重复检查邻居单元格。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         const unsigned int i = cell->active_cell_index(); 
 *         for (typename SparseMatrix<double>::iterator iter = 
 *                filter_matrix.begin(i); 
 *              iter != filter_matrix.end(i); 
 *              ++iter) 
 *           { 
 *             const unsigned int j     = iter->column(); 
 *             const double       value = iter->value() * cell->measure(); 
 * 
 *             system_matrix 
 *               .block(SolutionBlocks::unfiltered_density_multiplier, 
 *                      SolutionBlocks::unfiltered_density) 
 *               .add(i, j, value); 
 *             system_matrix 
 *               .block(SolutionBlocks::unfiltered_density, 
 *                      SolutionBlocks::unfiltered_density_multiplier) 
 *               .add(j, i, value); 
 *           } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SolvingtheNewtonlinearsystem"></a> 
 * <h3>Solving the Newton linear system</h3>
 * 

 * 
 * 我们将需要在每次迭代中解决一个线性系统。我们暂时使用一个直接求解器--对于一个有这么多非零值的矩阵来说，这显然不是一个有效的选择，而且它不会扩展到任何有趣的地方。对于 "真正的 "应用，我们将需要一个迭代求解器，但系统的复杂性意味着一个迭代求解器的算法将需要大量的工作。因为这不是当前程序的重点，所以我们简单地坚持使用我们在这里的直接求解器--该函数遵循与 step-29 中使用的相同结构。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BlockVector<double> SANDTopOpt<dim>::solve() 
 *   { 
 *     TimerOutput::Scope t(timer, "solver"); 
 * 
 *     BlockVector<double> linear_solution; 
 *     linear_solution.reinit(nonlinear_solution); 
 * 
 *     SparseDirectUMFPACK A_direct; 
 *     A_direct.initialize(system_matrix); 
 *     A_direct.vmult(linear_solution, system_rhs); 
 * 
 *     constraints.distribute(linear_solution); 
 * 
 *     return linear_solution; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Detailsoftheoptimizationalgorithm"></a> 
 * <h3>Details of the optimization algorithm</h3>
 * 

 * 
 * 接下来的几个函数处理优化算法的具体部分，最主要的是决定通过求解线性化（牛顿）系统计算出的方向是否可行，如果可行，我们要在这个方向上走多远。
 * 

 * 
 * 
 * <a name="Computingsteplengths"></a> 
 * <h4>Computing step lengths</h4>
 * 

 * 
 * 我们先用一个函数进行二进制搜索，找出符合对偶可行性的最大步骤--也就是说，我们能走多远，使  $s>0$  和  $z>0$  。该函数返回一对数值，分别代表 $s$ 和 $z$ 的松弛变量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<double, double> SANDTopOpt<dim>::calculate_max_step_size( 
 *     const BlockVector<double> &state, 
 *     const BlockVector<double> &step) const 
 *   { 
 *     double       fraction_to_boundary; 
 *     const double min_fraction_to_boundary = .8; 
 *     const double max_fraction_to_boundary = 1. - 1e-5; 
 * 
 *     if (min_fraction_to_boundary < 1 - barrier_size) 
 *       { 
 *         if (1 - barrier_size < max_fraction_to_boundary) 
 *           fraction_to_boundary = 1 - barrier_size; 
 *         else 
 *           fraction_to_boundary = max_fraction_to_boundary; 
 *       } 
 *     else 
 *       fraction_to_boundary = min_fraction_to_boundary; 
 * 
 *     double step_size_s_low  = 0; 
 *     double step_size_z_low  = 0; 
 *     double step_size_s_high = 1; 
 *     double step_size_z_high = 1; 
 *     double step_size_s, step_size_z; 
 * 
 *     const int max_bisection_method_steps = 50; 
 *     for (unsigned int k = 0; k < max_bisection_method_steps; ++k) 
 *       { 
 *         step_size_s = (step_size_s_low + step_size_s_high) / 2; 
 *         step_size_z = (step_size_z_low + step_size_z_high) / 2; 
 * 
 *         const BlockVector<double> state_test_s = 
 *           (fraction_to_boundary * state) + (step_size_s * step); 
 * 
 *         const BlockVector<double> state_test_z = 
 *           (fraction_to_boundary * state) + (step_size_z * step); 
 * 
 *         const bool accept_s = 
 *           (state_test_s.block(SolutionBlocks::density_lower_slack) 
 *              .is_non_negative()) && 
 *           (state_test_s.block(SolutionBlocks::density_upper_slack) 
 *              .is_non_negative()); 
 *         const bool accept_z = 
 *           (state_test_z.block(SolutionBlocks::density_lower_slack_multiplier) 
 *              .is_non_negative()) && 
 *           (state_test_z.block(SolutionBlocks::density_upper_slack_multiplier) 
 *              .is_non_negative()); 
 * 
 *         if (accept_s) 
 *           step_size_s_low = step_size_s; 
 *         else 
 *           step_size_s_high = step_size_s; 
 * 
 *         if (accept_z) 
 *           step_size_z_low = step_size_z; 
 *         else 
 *           step_size_z_high = step_size_z; 
 *       } 
 * 
 *     return {step_size_s_low, step_size_z_low}; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Computingresiduals"></a> 
 * <h4>Computing residuals</h4>
 * 

 * 
 * 下一个函数计算一个围绕 "测试解向量 "线性化的右手向量，我们可以用它来观察KKT条件的大小。 然后，这将用于在缩小障碍大小之前测试收敛性，以及计算 $l_1$ 的优点。
 * 

 * 
 * 这个函数冗长而复杂，但它实际上只是复制了上面`assemble_system()`函数的右侧部分的内容。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BlockVector<double> SANDTopOpt<dim>::calculate_test_rhs( 
 *     const BlockVector<double> &test_solution) const 
 *   { 
 * 
 * @endcode
 * 
 * 我们首先创建一个零向量，其大小和阻塞为system_rhs
 * 

 * 
 * 
 * @code
 *     BlockVector<double> test_rhs; 
 *     test_rhs.reinit(system_rhs); 
 * 
 *     MappingQGeneric<dim>  mapping(1); 
 *     const QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *     const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 *     FEValues<dim>         fe_values(mapping, 
 *                             fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 *     FEFaceValues<dim>     fe_face_values(mapping, 
 *                                      fe, 
 *                                      face_quadrature_formula, 
 *                                      update_values | update_quadrature_points | 
 *                                        update_normal_vectors | 
 *                                        update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.dofs_per_cell; 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     Vector<double>     cell_rhs(dofs_per_cell); 
 *     FullMatrix<double> dummy_cell_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     std::vector<double> lambda_values(n_q_points); 
 *     std::vector<double> mu_values(n_q_points); 
 * 
 *     const Functions::ConstantFunction<dim> lambda(1.), mu(1.); 
 *     std::vector<Tensor<1, dim>>            rhs_values(n_q_points); 
 * 
 *     BlockVector<double> filtered_unfiltered_density_solution = test_solution; 
 *     BlockVector<double> filter_adjoint_unfiltered_density_multiplier_solution = 
 *       test_solution; 
 *     filtered_unfiltered_density_solution.block( 
 *       SolutionBlocks::unfiltered_density) = 0; 
 *     filter_adjoint_unfiltered_density_multiplier_solution.block( 
 *       SolutionBlocks::unfiltered_density_multiplier) = 0; 
 * 
 *     filter_matrix.vmult(filtered_unfiltered_density_solution.block( 
 *                           SolutionBlocks::unfiltered_density), 
 *                         test_solution.block( 
 *                           SolutionBlocks::unfiltered_density)); 
 *     filter_matrix.Tvmult( 
 *       filter_adjoint_unfiltered_density_multiplier_solution.block( 
 *         SolutionBlocks::unfiltered_density_multiplier), 
 *       test_solution.block(SolutionBlocks::unfiltered_density_multiplier)); 
 * 
 *     std::vector<double>                  old_density_values(n_q_points); 
 *     std::vector<Tensor<1, dim>>          old_displacement_values(n_q_points); 
 *     std::vector<double>                  old_displacement_divs(n_q_points); 
 *     std::vector<SymmetricTensor<2, dim>> old_displacement_symmgrads(n_q_points); 
 *     std::vector<Tensor<1, dim>> old_displacement_multiplier_values(n_q_points); 
 *     std::vector<double>         old_displacement_multiplier_divs(n_q_points); 
 *     std::vector<SymmetricTensor<2, dim>> old_displacement_multiplier_symmgrads( 
 *       n_q_points); 
 *     std::vector<double> old_lower_slack_multiplier_values(n_q_points); 
 *     std::vector<double> old_upper_slack_multiplier_values(n_q_points); 
 *     std::vector<double> old_lower_slack_values(n_q_points); 
 *     std::vector<double> old_upper_slack_values(n_q_points); 
 *     std::vector<double> old_unfiltered_density_values(n_q_points); 
 *     std::vector<double> old_unfiltered_density_multiplier_values(n_q_points); 
 *     std::vector<double> filtered_unfiltered_density_values(n_q_points); 
 *     std::vector<double> filter_adjoint_unfiltered_density_multiplier_values( 
 *       n_q_points); 
 * 
 *     using namespace ValueExtractors; 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         cell_rhs = 0; 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         fe_values.reinit(cell); 
 * 
 *  
 *         mu.value_list(fe_values.get_quadrature_points(), mu_values); 
 * 
 *         fe_values[densities<dim>].get_function_values(test_solution, 
 *                                                       old_density_values); 
 *         fe_values[displacements<dim>].get_function_values( 
 *           test_solution, old_displacement_values); 
 *         fe_values[displacements<dim>].get_function_divergences( 
 *           test_solution, old_displacement_divs); 
 *         fe_values[displacements<dim>].get_function_symmetric_gradients( 
 *           test_solution, old_displacement_symmgrads); 
 *         fe_values[displacement_multipliers<dim>].get_function_values( 
 *           test_solution, old_displacement_multiplier_values); 
 *         fe_values[displacement_multipliers<dim>].get_function_divergences( 
 *           test_solution, old_displacement_multiplier_divs); 
 *         fe_values[displacement_multipliers<dim>] 
 *           .get_function_symmetric_gradients( 
 *             test_solution, old_displacement_multiplier_symmgrads); 
 *         fe_values[density_lower_slacks<dim>].get_function_values( 
 *           test_solution, old_lower_slack_values); 
 *         fe_values[density_lower_slack_multipliers<dim>].get_function_values( 
 *           test_solution, old_lower_slack_multiplier_values); 
 *         fe_values[density_upper_slacks<dim>].get_function_values( 
 *           test_solution, old_upper_slack_values); 
 *         fe_values[density_upper_slack_multipliers<dim>].get_function_values( 
 *           test_solution, old_upper_slack_multiplier_values); 
 *         fe_values[unfiltered_densities<dim>].get_function_values( 
 *           test_solution, old_unfiltered_density_values); 
 *         fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
 *           test_solution, old_unfiltered_density_multiplier_values); 
 *         fe_values[unfiltered_densities<dim>].get_function_values( 
 *           filtered_unfiltered_density_solution, 
 *           filtered_unfiltered_density_values); 
 *         fe_values[unfiltered_density_multipliers<dim>].get_function_values( 
 *           filter_adjoint_unfiltered_density_multiplier_solution, 
 *           filter_adjoint_unfiltered_density_multiplier_values); 
 * 
 *         for (const auto q_point : fe_values.quadrature_point_indices()) 
 *           { 
 *             for (const auto i : fe_values.dof_indices()) 
 *               { 
 *                 const SymmetricTensor<2, dim> displacement_phi_i_symmgrad = 
 *                   fe_values[displacements<dim>].symmetric_gradient(i, q_point); 
 *                 const double displacement_phi_i_div = 
 *                   fe_values[displacements<dim>].divergence(i, q_point); 
 * 
 *                 const SymmetricTensor<2, dim> 
 *                   displacement_multiplier_phi_i_symmgrad = 
 *                     fe_values[displacement_multipliers<dim>].symmetric_gradient( 
 *                       i, q_point); 
 *                 const double displacement_multiplier_phi_i_div = 
 *                   fe_values[displacement_multipliers<dim>].divergence(i, 
 *                                                                       q_point); 
 * 
 *                 const double density_phi_i = 
 *                   fe_values[densities<dim>].value(i, q_point); 
 *                 const double unfiltered_density_phi_i = 
 *                   fe_values[unfiltered_densities<dim>].value(i, q_point); 
 *                 const double unfiltered_density_multiplier_phi_i = 
 *                   fe_values[unfiltered_density_multipliers<dim>].value(i, 
 *                                                                        q_point); 
 * 
 *                 const double lower_slack_multiplier_phi_i = 
 *                   fe_values[density_lower_slack_multipliers<dim>].value( 
 *                     i, q_point); 
 * 
 *                 const double lower_slack_phi_i = 
 *                   fe_values[density_lower_slacks<dim>].value(i, q_point); 
 * 
 *                 const double upper_slack_phi_i = 
 *                   fe_values[density_upper_slacks<dim>].value(i, q_point); 
 * 
 *                 const double upper_slack_multiplier_phi_i = 
 *                   fe_values[density_upper_slack_multipliers<dim>].value( 
 *                     i, q_point); 
 * 
 *                 /* 方程1：这个方程以及方程
 *                  * 2 and 3, are the variational derivatives of the 
 *                  * Lagrangian with respect to the decision 
 *                  * variables - the density, displacement, and 
 *                  * unfiltered density. */ 
 * 
 * 
 *                 cell_rhs(i) += 
 *                   -1 * fe_values.JxW(q_point) * 
 *                   (density_penalty_exponent * 
 *                      std::pow(old_density_values[q_point], 
 *                               density_penalty_exponent - 1) * 
 *                      density_phi_i * 
 *                      (old_displacement_multiplier_divs[q_point] * 
 *                         old_displacement_divs[q_point] * 
 *                         lambda_values[q_point] + 
 *                       2 * mu_values[q_point] * 
 *                         (old_displacement_symmgrads[q_point] * 
 *                          old_displacement_multiplier_symmgrads[q_point])) - 
 *                    density_phi_i * 
 *                      old_unfiltered_density_multiplier_values[q_point]); 
 * 
 *                 /*方程2；边界项将被进一步添加。
 *                  * below. */ 
 * 
 * 
 *                 cell_rhs(i) += 
 *                   -1 * fe_values.JxW(q_point) * 
 *                   (std::pow(old_density_values[q_point], 
 *                             density_penalty_exponent) * 
 *                    (old_displacement_multiplier_divs[q_point] * 
 *                       displacement_phi_i_div * lambda_values[q_point] + 
 *                     2 * mu_values[q_point] * 
 *                       (old_displacement_multiplier_symmgrads[q_point] * 
 *                        displacement_phi_i_symmgrad))); 
 * @endcode
 * 
 * 
 * 
 * 
 * @code
 *                /* 方程3  */ 
 *                 cell_rhs(i) += 
 *                   -1 * fe_values.JxW(q_point) * 
 *                   (unfiltered_density_phi_i * 
 *                      filter_adjoint_unfiltered_density_multiplier_values 
 *                        [q_point] + 
 *                    unfiltered_density_phi_i * 
 *                      old_upper_slack_multiplier_values[q_point] + 
 *                    -1 * unfiltered_density_phi_i * 
 *                      old_lower_slack_multiplier_values[q_point]); 
 * 
 *                /* 方程4；边界项将再次被处理。with below. 
 *                 * This equation being driven to 0 ensures that the elasticity 
 *                 * equation is met as a constraint. */ 
 *                 cell_rhs(i) += -1 * fe_values.JxW(q_point) * 
 *                                (std::pow(old_density_values[q_point], 
 *                                          density_penalty_exponent) * 
 *                                 (old_displacement_divs[q_point] * 
 *                                    displacement_multiplier_phi_i_div * 
 *                                    lambda_values[q_point] + 
 *                                  2 * mu_values[q_point] * 
 *                                    (displacement_multiplier_phi_i_symmgrad * 
 *                                     old_displacement_symmgrads[q_point]))); 
 * 
 *                 /* 方程5：该方程设定了下限的松弛量， giving a minimum density of 0. */ 
 *                 cell_rhs(i) += fe_values.JxW(q_point) * 
 *                                (lower_slack_multiplier_phi_i * 
 *                                 (old_unfiltered_density_values[q_point] - 
 *                                  old_lower_slack_values[q_point])); 
 * 
 *                 /* 方程6：该方程设定了上层松弛量variable equal to one minus the unfiltered density. */ 
 *                 cell_rhs(i) += fe_values.JxW(q_point) * 
 *                                (upper_slack_multiplier_phi_i * 
 *                                 (1 - old_unfiltered_density_values[q_point] - 
 *                                  old_upper_slack_values[q_point])); 
 * 
 *                 /*等式7：这是在
 *                  * density and the filter applied to the 
 *                  * unfiltered density. This being driven to 0 by 
 *                  * the Newton steps ensures that the filter is 
 *                  * applied correctly. */ 
 *                 cell_rhs(i) += fe_values.JxW(q_point) * 
 *                                (unfiltered_density_multiplier_phi_i * 
 *                                 (old_density_values[q_point] - 
 *                                  filtered_unfiltered_density_values[q_point])); 
 * 
 *                 /*方程8：这与方程9一起给出了
 *                  * requirement that s*z = \alpha for the barrier 
 *                  * size alpha, and gives complementary slackness 
 *                  * from KKT conditions when \alpha goes to 0. */ 
 *                 cell_rhs(i) += 
 *                   -1 * fe_values.JxW(q_point) * 
 *                   (lower_slack_phi_i * 
 *                    (old_lower_slack_multiplier_values[q_point] - 
 *                     barrier_size / old_lower_slack_values[q_point])); 
 * 
 *                 /*方程9  */ 
 *                 cell_rhs(i) += 
 *                   -1 * fe_values.JxW(q_point) * 
 *                   (upper_slack_phi_i * 
 *                    (old_upper_slack_multiplier_values[q_point] - 
 *                     barrier_size / old_upper_slack_values[q_point])); 
 *               } 
 *           } 
 * 
 *         for (const auto &face : cell->face_iterators()) 
 *           { 
 *             if (face->at_boundary() && 
 *                 face->boundary_id() == BoundaryIds::down_force) 
 *               { 
 *                 fe_face_values.reinit(cell, face); 
 * 
 *                 for (const auto face_q_point : 
 *                      fe_face_values.quadrature_point_indices()) 
 *                   { 
 *                     for (const auto i : fe_face_values.dof_indices()) 
 *                       { 
 *                         Tensor<1, dim> traction; 
 *                         traction[1] = -1.; 
 * 
 *                         cell_rhs(i) += 
 *                           -1 * 
 *                           (traction * fe_face_values[displacements<dim>].value( 
 *                                         i, face_q_point)) * 
 *                           fe_face_values.JxW(face_q_point); 
 * 
 *                         cell_rhs(i) += 
 *                           (traction * 
 *                            fe_face_values[displacement_multipliers<dim>].value( 
 *                              i, face_q_point)) * 
 *                           fe_face_values.JxW(face_q_point); 
 *                       } 
 *                   } 
 *               } 
 *           } 
 * 
 *         MatrixTools::local_apply_boundary_values(boundary_values, 
 *                                                  local_dof_indices, 
 *                                                  dummy_cell_matrix, 
 *                                                  cell_rhs, 
 *                                                  true); 
 * 
 *         constraints.distribute_local_to_global(cell_rhs, 
 *                                                local_dof_indices, 
 *                                                test_rhs); 
 *       } 
 * 
 *     return test_rhs; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Computingthemeritfunction"></a> 
 * <h4>Computing the merit function</h4>
 * 

 * 
 * 我们在这里使用的算法使用一个 "看门狗 "策略来确定从当前迭代的位置和程度。 我们将看门狗策略建立在一个精确的 $l_1$ 功绩函数上。这个函数计算一个给定的、假定的、下一个迭代的精确 $l_1$ 功绩。
 * 

 * 
 * 优点函数由目标函数的总和（简单来说就是外力的积分（在域的边界上）乘以测试解的位移值（通常是当前解加上牛顿更新的某个倍数），以及残差向量的拉格朗日乘数分量的 $l_1$ 准则组成。下面的代码依次计算这些部分。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double SANDTopOpt<dim>::calculate_exact_merit( 
 *     const BlockVector<double> &test_solution) 
 *   { 
 *     TimerOutput::Scope t(timer, "merit function"); 
 * 
 * @endcode
 * 
 * 从计算目标函数开始。
 * 
 * @code
 *     double objective_function_merit = 0; 
 *     { 
 *       MappingQGeneric<dim>  mapping(1); 
 *       const QGauss<dim>     quadrature_formula(fe.degree + 1); 
 *       const QGauss<dim - 1> face_quadrature_formula(fe.degree + 1); 
 *       FEValues<dim>         fe_values(mapping, 
 *                               fe, 
 *                               quadrature_formula, 
 *                               update_values | update_gradients | 
 *                                 update_quadrature_points | update_JxW_values); 
 *       FEFaceValues<dim>     fe_face_values(mapping, 
 *                                        fe, 
 *                                        face_quadrature_formula, 
 *                                        update_values | 
 *                                          update_quadrature_points | 
 *                                          update_normal_vectors | 
 *                                          update_JxW_values); 
 * 
 *       const unsigned int n_face_q_points = face_quadrature_formula.size(); 
 * 
 *       std::vector<Tensor<1, dim>> displacement_face_values(n_face_q_points); 
 * 
 *       for (const auto &cell : dof_handler.active_cell_iterators()) 
 *         { 
 *           for (const auto &face : cell->face_iterators()) 
 *             { 
 *               if (face->at_boundary() && 
 *                   face->boundary_id() == BoundaryIds::down_force) 
 *                 { 
 *                   fe_face_values.reinit(cell, face); 
 *                   fe_face_values[ValueExtractors::displacements<dim>] 
 *                     .get_function_values(test_solution, 
 *                                          displacement_face_values); 
 *                   for (unsigned int face_q_point = 0; 
 *                        face_q_point < n_face_q_points; 
 *                        ++face_q_point) 
 *                     { 
 *                       Tensor<1, dim> traction; 
 *                       traction[1] = -1.; 
 * 
 *                       objective_function_merit += 
 *                         (traction * displacement_face_values[face_q_point]) * 
 *                         fe_face_values.JxW(face_q_point); 
 *                     } 
 *                 } 
 *             } 
 *         } 
 *     } 
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       { 
 *         objective_function_merit = 
 *           objective_function_merit - 
 *           barrier_size * cell->measure() * 
 *             std::log(test_solution.block( 
 *               SolutionBlocks::density_lower_slack)[cell->active_cell_index()]); 
 *         objective_function_merit = 
 *           objective_function_merit - 
 *           barrier_size * cell->measure() * 
 *             std::log(test_solution.block( 
 *               SolutionBlocks::density_upper_slack)[cell->active_cell_index()]); 
 *       } 
 * @endcode
 * 
 * 然后
 * 计算残差，并取对应于拉格朗日多边形的组件的 $l_1$ 准则。我们把这些加到上面计算的目标函数中，并在底部返回总和。
 * 
 * @code
 *     const BlockVector<double> test_rhs = calculate_test_rhs(test_solution); 
 * 
 *     const double elasticity_constraint_merit = 
 *       penalty_multiplier * 
 *       test_rhs.block(SolutionBlocks::displacement_multiplier).l1_norm(); 
 *     const double filter_constraint_merit = 
 *       penalty_multiplier * 
 *       test_rhs.block(SolutionBlocks::unfiltered_density_multiplier).l1_norm(); 
 *     const double lower_slack_merit = 
 *       penalty_multiplier * 
 *       test_rhs.block(SolutionBlocks::density_lower_slack_multiplier).l1_norm(); 
 *     const double upper_slack_merit = 
 *       penalty_multiplier * 
 *       test_rhs.block(SolutionBlocks::density_upper_slack_multiplier).l1_norm(); 
 * 
 *     const double total_merit = 
 *       objective_function_merit + elasticity_constraint_merit + 
 *       filter_constraint_merit + lower_slack_merit + upper_slack_merit; 
 *     return total_merit; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Findingasearchdirection"></a> 
 * <h4>Finding a search direction</h4>
 * 

 * 
 * 接下来是实际计算从当前状态（作为第一个参数传递）开始的搜索方向并返回结果向量的函数。为此，该函数首先调用与牛顿系统相对应的线性系统的组合函数，并对其进行求解。
 * 

 * 
 * 这个函数还更新了优点函数中的惩罚乘数，然后返回最大比例的可行步骤。它使用`calculate_max_step_sizes()`函数来找到满足  $s>0$  和  $z>0$  的最大可行步骤。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BlockVector<double> SANDTopOpt<dim>::find_max_step() 
 *   { 
 *     assemble_system(); 
 *     BlockVector<double> step = solve(); 
 * 
 * @endcode
 * 
 * 接下来我们要更新punice_multiplier。 从本质上讲，更大的惩罚乘数使我们更多考虑约束条件。 观察与我们的决策变量有关的Hessian和梯度，并将其与我们的约束误差的规范相比较，可以确保我们的优点函数是 "精确的"
 * 

 * 
 * 也就是说，它在与目标函数相同的位置有一个最小值。 由于我们的优点函数对任何超过某个最小值的惩罚乘数都是精确的，所以我们只保留计算值，如果它增加了惩罚乘数。
 * 

 * 
 * 
 * @code
 *     const std::vector<unsigned int> decision_variables = { 
 *       SolutionBlocks::density, 
 *       SolutionBlocks::displacement, 
 *       SolutionBlocks::unfiltered_density, 
 *       SolutionBlocks::density_upper_slack, 
 *       SolutionBlocks::density_lower_slack}; 
 *     double hess_part = 0; 
 *     double grad_part = 0; 
 *     for (const unsigned int decision_variable_i : decision_variables) 
 *       { 
 *         for (const unsigned int decision_variable_j : decision_variables) 
 *           { 
 *             Vector<double> temp_vector(step.block(decision_variable_i).size()); 
 *             system_matrix.block(decision_variable_i, decision_variable_j) 
 *               .vmult(temp_vector, step.block(decision_variable_j)); 
 *             hess_part += step.block(decision_variable_i) * temp_vector; 
 *           } 
 *         grad_part -= system_rhs.block(decision_variable_i) * 
 *                      step.block(decision_variable_i); 
 *       } 
 * 
 *     const std::vector<unsigned int> equality_constraint_multipliers = { 
 *       SolutionBlocks::displacement_multiplier, 
 *       SolutionBlocks::unfiltered_density_multiplier, 
 *       SolutionBlocks::density_lower_slack_multiplier, 
 *       SolutionBlocks::density_upper_slack_multiplier}; 
 *     double constraint_norm = 0; 
 *     for (unsigned int multiplier_i : equality_constraint_multipliers) 
 *       constraint_norm += system_rhs.block(multiplier_i).linfty_norm(); 
 * 
 *     double test_penalty_multiplier; 
 *     if (hess_part > 0) 
 *       test_penalty_multiplier = 
 *         (grad_part + .5 * hess_part) / (.05 * constraint_norm); 
 *     else 
 *       test_penalty_multiplier = (grad_part) / (.05 * constraint_norm); 
 * 
 *     penalty_multiplier = std::max(penalty_multiplier, test_penalty_multiplier); 
 * 
 * @endcode
 * 
 * 基于所有这些，我们现在可以计算出原始变量和对偶变量（拉格朗日乘数）的步长。一旦我们有了这些，我们就可以对解向量的分量进行缩放，这就是这个函数的回报。
 * 

 * 
 * 
 * @code
 *     const std::pair<double, double> max_step_sizes = 
 *       calculate_max_step_size(nonlinear_solution, step); 
 *     const double step_size_s = max_step_sizes.first; 
 *     const double step_size_z = max_step_sizes.second; 
 * 
 *     step.block(SolutionBlocks::density) *= step_size_s; 
 *     step.block(SolutionBlocks::displacement) *= step_size_s; 
 *     step.block(SolutionBlocks::unfiltered_density) *= step_size_s; 
 *     step.block(SolutionBlocks::displacement_multiplier) *= step_size_z; 
 *     step.block(SolutionBlocks::unfiltered_density_multiplier) *= step_size_z; 
 *     step.block(SolutionBlocks::density_lower_slack) *= step_size_s; 
 *     step.block(SolutionBlocks::density_lower_slack_multiplier) *= step_size_z; 
 *     step.block(SolutionBlocks::density_upper_slack) *= step_size_s; 
 *     step.block(SolutionBlocks::density_upper_slack_multiplier) *= step_size_z; 
 * 
 *     return step; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Computingascaledstep"></a> 
 * <h4>Computing a scaled step</h4>
 * 

 * 
 * 下一个函数接着实现了直线搜索的反向跟踪算法。它不断缩小步长，直到找到一个优点减少的步长，然后根据当前的状态向量，以及要进入的方向，乘以步长，返回新的位置。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BlockVector<double> 
 *   SANDTopOpt<dim>::compute_scaled_step(const BlockVector<double> &state, 
 *                                        const BlockVector<double> &max_step, 
 *                                        const double descent_requirement) 
 *   { 
 *     const double merit_derivative = 
 *       (calculate_exact_merit(state + 1e-4 * max_step) - 
 *        calculate_exact_merit(state)) / 
 *       1e-4; 
 *     double       step_size                 = 1; 
 *     unsigned int max_linesearch_iterations = 10; 
 *     for (unsigned int k = 0; k < max_linesearch_iterations; ++k) 
 *       { 
 *         if (calculate_exact_merit(state + step_size * max_step) < 
 *             calculate_exact_merit(state) + 
 *               step_size * descent_requirement * merit_derivative) 
 *           break; 
 *         else 
 *           step_size = step_size / 2; 
 *       } 
 *     return state + (step_size * max_step); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Checkingforconvergence"></a> 
 * <h4>Checking for convergence</h4>
 * 

 * 
 * 本块中的最后一个辅助函数是检查是否充分满足KKT条件，以便整个算法可以降低障碍物的大小。它通过计算残差的 $l_1$ 准则来实现，这就是`calculate_test_rhs()`的计算。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   bool SANDTopOpt<dim>::check_convergence(const BlockVector<double> &state) 
 *   { 
 *     const BlockVector<double> test_rhs      = calculate_test_rhs(state); 
 *     const double              test_rhs_norm = test_rhs.l1_norm(); 
 * 
 *     const double convergence_condition = 1e-2; 
 *     const double target_norm           = convergence_condition * barrier_size; 
 * 
 *     std::cout << "    Checking convergence. Current rhs norm is " 
 *               << test_rhs_norm << ", target is " << target_norm << std::endl; 
 * 
 *     return (test_rhs_norm < target_norm); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Postprocessingthesolution"></a> 
 * <h3>Postprocessing the solution</h3>
 * 

 * 
 * 后处理函数中的第一个函数在VTU文件中输出信息，用于可视化。它看起来很长，但实际上与  step-22  中所做的一样，例如，只是增加了（很多）解决方案的变量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::output_results(const unsigned int iteration) const 
 *   { 
 *     std::vector<std::string> solution_names(1, "density"); 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         1, DataComponentInterpretation::component_is_scalar); 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       { 
 *         solution_names.emplace_back("displacement"); 
 *         data_component_interpretation.push_back( 
 *           DataComponentInterpretation::component_is_part_of_vector); 
 *       } 
 *     solution_names.emplace_back("unfiltered_density"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     for (unsigned int i = 0; i < dim; ++i) 
 *       { 
 *         solution_names.emplace_back("displacement_multiplier"); 
 *         data_component_interpretation.push_back( 
 *           DataComponentInterpretation::component_is_part_of_vector); 
 *       } 
 *     solution_names.emplace_back("unfiltered_density_multiplier"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     solution_names.emplace_back("low_slack"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     solution_names.emplace_back("low_slack_multiplier"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     solution_names.emplace_back("high_slack"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     solution_names.emplace_back("high_slack_multiplier"); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(nonlinear_solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output("solution" + std::to_string(iteration) + ".vtu"); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * @endcode
 * 
 * 其中第二个函数将解决方案输出为`.stl`文件，用于3D打印。STL](https:en.wikipedia.org/wiki/STL_(file_format))文件是由三角形和法线向量组成的，我们将用它来显示所有那些密度值大于0的单元，首先将网格从 $z$ 值挤出到 $z=0.25$  ，然后为密度值足够大的单元的每个面生成两个三角形。当从外面看时，三角形节点必须逆时针走，法向量必须是指向外部的单位向量，这需要进行一些检查。
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::write_as_stl() 
 *   { 
 *     static_assert(dim == 2, 
 *                   "This function is not implemented for anything " 
 *                   "other than the 2d case."); 
 * 
 *     std::ofstream stlfile; 
 *     stlfile.open("bridge.stl"); 
 * 
 *     stlfile << "solid bridge\n" << std::scientific; 
 *     double height = .25; 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         if (nonlinear_solution.block( 
 *               SolutionBlocks::density)[cell->active_cell_index()] > 0.5) 
 *           { 
 * @endcode
 * 
 * 我们现在已经找到了一个密度值大于0的单元。让我们先写出底部和顶部的面。由于上面提到的排序问题，我们必须确保了解一个单元的坐标系是右旋的还是左旋的。我们通过询问从顶点0开始的两条边的方向以及它们是否形成一个右手坐标系来做到这一点。
 * 
 * @code
 *             const Tensor<1, dim> edge_directions[2] = {cell->vertex(1) - 
 *                                                          cell->vertex(0), 
 *                                                        cell->vertex(2) - 
 *                                                          cell->vertex(0)}; 
 *             const Tensor<2, dim> edge_tensor( 
 *               {{edge_directions[0][0], edge_directions[0][1]}, 
 *                {edge_directions[1][0], edge_directions[1][1]}}); 
 *             const bool is_right_handed_cell = (determinant(edge_tensor) > 0); 
 * 
 *             if (is_right_handed_cell) 
 *               { 
 * 
 *                /*在z=0处写出一个边。  */ 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(0)[0] << " " 
 *                         << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(3)[0] << " " 
 *                         << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 * 
 *                /*在z=高度处写下一个边。  */  
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(0)[0] << " " 
 *                         << cell->vertex(0)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << height << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(3)[0] << " " 
 *                         << cell->vertex(3)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << height << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *               } 
 *             else /* The cell has a left-handed set up */ 
 *               { 
 *                /* 在z=0处写出一边。  */ 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " "
 *                         << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(0)[0] << " " 
 *                         << cell->vertex(0)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << -1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(3)[0] << " " 
 *                         << cell->vertex(3)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << 0.000000e+00 << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 * 
 *                /*在z=高度处写出一个边。  */ 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(0)[0] << " " 
 *                         << cell->vertex(0)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << height << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *                 stlfile << "   facet normal " << 0.000000e+00 << " " 
 *                         << 0.000000e+00 << " " << 1.000000e+00 << "\n"; 
 *                 stlfile << "      outer loop\n"; 
 *                 stlfile << "         vertex " << cell->vertex(1)[0] << " " 
 *                         << cell->vertex(1)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(2)[0] << " " 
 *                         << cell->vertex(2)[1] << " " << height << "\n"; 
 *                 stlfile << "         vertex " << cell->vertex(3)[0] << " " 
 *                         << cell->vertex(3)[1] << " " << height << "\n"; 
 *                 stlfile << "      endloop\n"; 
 *                 stlfile << "   endfacet\n"; 
 *               } 
 * 
 * @endcode
 * 
 * 接下来我们需要处理单元格的四个面，扩展到 $z$ 方向。然而，我们只需要写这些面，如果该面在域的边界上，或者它是密度大于0.5的单元和密度小于0.5的单元之间的界面。
 * 

 * 
 * 
 * @code
 *             for (unsigned int face_number = 0; 
 *                  face_number < GeometryInfo<dim>::faces_per_cell; 
 *                  ++face_number) 
 *               { 
 *                 const typename DoFHandler<dim>::face_iterator face = 
 *                   cell->face(face_number); 
 * 
 *                 if ((face->at_boundary()) || 
 *                     (!face->at_boundary() && 
 *                      (nonlinear_solution.block( 
 *                         0)[cell->neighbor(face_number)->active_cell_index()] < 
 *                       0.5))) 
 *                   { 
 *                     const Tensor<1, dim> normal_vector = 
 *                       (face->center() - cell->center()); 
 *                     const double normal_norm = normal_vector.norm(); 
 *                     if ((face->vertex(0)[0] - face->vertex(0)[0]) * 
 *                             (face->vertex(1)[1] - face->vertex(0)[1]) * 
 *                             0.000000e+00 + 
 *                           (face->vertex(0)[1] - face->vertex(0)[1]) * (0 - 0) * 
 *                             normal_vector[0] + 
 *                           (height - 0) * 
 *                             (face->vertex(1)[0] - face->vertex(0)[0]) * 
 *                             normal_vector[1] - 
 *                           (face->vertex(0)[0] - face->vertex(0)[0]) * (0 - 0) * 
 *                             normal_vector[1] - 
 *                           (face->vertex(0)[1] - face->vertex(0)[1]) * 
 *                             (face->vertex(1)[0] - face->vertex(0)[0]) * 
 *                             normal_vector[0] - 
 *                           (height - 0) * 
 *                             (face->vertex(1)[1] - face->vertex(0)[1]) * 0 > 
 *                         0) 
 *                       { 
 *                         stlfile << "   facet normal " 
 *                                 << normal_vector[0] / normal_norm << " " 
 *                                 << normal_vector[1] / normal_norm << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      outer loop\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      endloop\n"; 
 *                         stlfile << "   endfacet\n"; 
 *                         stlfile << "   facet normal " 
 *                                 << normal_vector[0] / normal_norm << " " 
 *                                 << normal_vector[1] / normal_norm << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      outer loop\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      endloop\n"; 
 *                         stlfile << "   endfacet\n"; 
 *                       } 
 *                     else 
 *                       { 
 *                         stlfile << "   facet normal " 
 *                                 << normal_vector[0] / normal_norm << " " 
 *                                 << normal_vector[1] / normal_norm << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      outer loop\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "      endloop\n"; 
 *                         stlfile << "   endfacet\n"; 
 *                         stlfile << "   facet normal " 
 *                                 << normal_vector[0] / normal_norm << " " 
 *                                 << normal_vector[1] / normal_norm << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "      outer loop\n"; 
 *                         stlfile << "         vertex " << face->vertex(0)[0] 
 *                                 << " " << face->vertex(0)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " 
 *                                 << 0.000000e+00 << "\n"; 
 *                         stlfile << "         vertex " << face->vertex(1)[0] 
 *                                 << " " << face->vertex(1)[1] << " " << height 
 *                                 << "\n"; 
 *                         stlfile << "      endloop\n"; 
 *                         stlfile << "   endfacet\n"; 
 *                       } 
 *                   } 
 *               } 
 *           } 
 *       } 
 *     stlfile << "endsolid bridge"; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Therunfunctiondrivingtheoverallalgorithm"></a> 
 * <h3>The run() function driving the overall algorithm</h3>
 * 

 * 
 * 这个函数最终提供了整体的驱动逻辑。从总体上看，这是一个相当复杂的函数，主要是因为优化算法很困难：它不仅仅是像 step-15 中那样找到一个牛顿方向，然后在这个方向上再走一个固定的距离，而是要（i）确定当前步骤中的最佳对数障碍惩罚参数应该是什么，（ii）通过复杂的算法来确定我们要走多远，还有其他成分。让我们看看如何在下面的文件中把它分解成小块。
 * 

 * 
 * 该函数一开始就很简单，首先设置了网格、DoFHandler，然后是下面所需的各种线性代数对象。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SANDTopOpt<dim>::run() 
 *   { 
 *     std::cout << "filter r is: " << filter_r << std::endl; 
 * 
 *     { 
 *       TimerOutput::Scope t(timer, "setup"); 
 * 
 *       create_triangulation(); 
 * 
 *       dof_handler.distribute_dofs(fe); 
 *       DoFRenumbering::component_wise(dof_handler); 
 * 
 *       setup_boundary_values(); 
 *       setup_block_system(); 
 *       setup_filter_matrix(); 
 *     } 
 * 
 * @endcode
 * 
 * 然后，我们设置一些影响优化算法的对数屏障和直线搜索部分的参数。
 * 

 * 
 * 
 * @code
 *     barrier_size                  = 25; 
 *     const double min_barrier_size = .0005; 
 * 
 *     const unsigned int max_uphill_steps    = 8; 
 *     const double       descent_requirement = .0001; 
 * 
 * @endcode
 * 
 * 现在开始进行主迭代。整个算法通过使用一个外循环来工作，在这个外循环中，我们一直循环到（i）对数障碍参数变得足够小，或者（ii）我们已经达到收敛。在任何情况下，如果最终的迭代次数过多，我们就会终止。这个整体结构被编码为一个 "do{ ... } while (...)`循环，其中收敛条件在底部。
 * 

 * 
 * 
 * @code
 *     unsigned int       iteration_number = 0; 
 *     const unsigned int max_iterations   = 10000; 
 * 
 *     do 
 *       { 
 *         std::cout << "Starting outer step in iteration " << iteration_number 
 *                   << " with barrier parameter " << barrier_size << std::endl; 
 * 
 * @endcode
 * 
 * 在这个外循环中，我们有一个内循环，在这个内循环中，我们试图使用介绍中描述的看门狗算法找到一个更新方向。
 * 

 * 
 * 看门狗算法本身的总体思路是这样的。对于最大的`max_uphill_steps`（即上述 "内循环 "中的一个循环）的尝试，我们使用`find_max_step()`来计算牛顿更新步骤，并在`nonlinear_solution`向量中加上这些。 在每一次尝试中（从上一次尝试结束时到达的地方开始），我们检查我们是否已经达到了上述优点函数的目标值。目标值是根据本算法的起始位置（看门狗循环开始时的`nonlinear_solution'，保存为`看门狗_state'）和本循环第一个回合中`find_max_step()'提供的第一个建议方向（`k=0'情况）计算的。
 * 

 * 
 * 
 * @code
 *         do 
 *           { 
 *             std::cout << "  Starting inner step in iteration " 
 *                       << iteration_number 
 *                       << " with merit function penalty multiplier " 
 *                       << penalty_multiplier << std::endl; 
 * 
 *             bool watchdog_step_found = false; 
 * 
 *             const BlockVector<double> watchdog_state = nonlinear_solution; 
 *             BlockVector<double>       first_step; 
 *             double target_merit     = numbers::signaling_nan<double>(); 
 *             double merit_derivative = numbers::signaling_nan<double>(); 
 * 
 *             for (unsigned int k = 0; k < max_uphill_steps; ++k) 
 *               { 
 *                 ++iteration_number; 
 *                 const BlockVector<double> update_step = find_max_step(); 
 * 
 *                 if (k == 0) 
 *                   { 
 *                     first_step = update_step; 
 *                     merit_derivative = 
 *                       ((calculate_exact_merit(watchdog_state + 
 *                                               .0001 * first_step) - 
 *                         calculate_exact_merit(watchdog_state)) / 
 *                        .0001); 
 *                     target_merit = calculate_exact_merit(watchdog_state) + 
 *                                    descent_requirement * merit_derivative; 
 *                   } 
 * 
 *                 nonlinear_solution += update_step; 
 *                 const double current_merit = 
 *                   calculate_exact_merit(nonlinear_solution); 
 * 
 *                 std::cout << "    current watchdog state merit is: " 
 *                           << current_merit << "; target merit is " 
 *                           << target_merit << std::endl; 
 * 
 *                 if (current_merit < target_merit) 
 *                   { 
 *                     watchdog_step_found = true; 
 *                     std::cout << "    found workable step after " << k + 1 
 *                               << " iterations" << std::endl; 
 *                     break; 
 *                   } 
 *               } 
 * @endcode
 * 
 * 然后
 * 算法的下一部分取决于上面的看门狗循环是否成功。如果成功了，那么我们就满意了，不需要进一步的行动。我们只是停留在原地。然而，如果我们在上面的循环中采取了最大数量的不成功的步骤，那么我们就需要做一些别的事情，这就是下面的代码块所做的。    具体来说，从上述循环的最后（不成功的）状态开始，我们再寻找一个更新方向，并采取所谓的 "伸展步骤"。如果该拉伸状态满足涉及优点函数的条件，那么我们就去那里。另一方面，如果拉伸状态也是不可接受的（就像上面所有的看门狗步骤一样），那么我们就放弃上面所有的看门狗步骤，在我们开始看门狗迭代的地方重新开始--那个地方被存储在上面的`看门狗_状态`变量中。更具体地说，下面的条件首先测试我们是否从`看门狗_state`方向的`first_step`走了一步，或者我们是否可以从拉伸状态再做一次更新来找到一个新的地方。有可能这两种情况实际上都不比我们在看门狗算法开始时的状态好，但即使是这样，那个地方显然是个困难的地方，离开后从另一个地方开始下一次迭代可能是一个有用的策略，最终收敛。    我们不断重复上面的看门狗步骤以及下面的逻辑，直到这个内部迭代最终收敛（或者如果我们遇到最大的迭代次数--在这里我们把线性求解的次数算作迭代次数，并在每次调用`find_max_step()`时增加计数器，因为这就是线性求解实际发生的地方）。在任何情况下，在这些内部迭代的每一次结束时，我们也会以适合可视化的形式输出解决方案。
 * 

 * 
 * 
 * @code
 *             if (watchdog_step_found == false) 
 *               { 
 *                 ++iteration_number; 
 *                 const BlockVector<double> update_step = find_max_step(); 
 *                 const BlockVector<double> stretch_state = 
 *                   compute_scaled_step(nonlinear_solution, 
 *                                       update_step, 
 *                                       descent_requirement); 
 * 
 * @endcode
 * 
 * 如果我们没有得到一个成功的看门狗步骤，我们现在需要决定是回到我们开始的地方，还是使用最终状态。 我们比较这两个位置的优劣，然后从哪个位置取一个按比例的步长。 由于按比例的步长可以保证降低优点，所以我们最终会保留这两个位置中的一个。
 * 

 * 
 * 
 * @code
 *                 if ((calculate_exact_merit(nonlinear_solution) < 
 *                      calculate_exact_merit(watchdog_state)) || 
 *                     (calculate_exact_merit(stretch_state) < target_merit)) 
 *                   { 
 *                     std::cout << "    Taking scaled step from end of watchdog" 
 *                               << std::endl; 
 *                     nonlinear_solution = stretch_state; 
 *                   } 
 *                 else 
 *                   { 
 *                     std::cout 
 *                       << "    Taking scaled step from beginning of watchdog" 
 *                       << std::endl; 
 *                     if (calculate_exact_merit(stretch_state) > 
 *                         calculate_exact_merit(watchdog_state)) 
 *                       { 
 *                         nonlinear_solution = 
 *                           compute_scaled_step(watchdog_state, 
 *                                               first_step, 
 *                                               descent_requirement); 
 *                       } 
 *                     else 
 *                       { 
 *                         ++iteration_number; 
 *                         nonlinear_solution = stretch_state; 
 *                         const BlockVector<double> stretch_step = 
 *                           find_max_step(); 
 *                         nonlinear_solution = 
 *                           compute_scaled_step(nonlinear_solution, 
 *                                               stretch_step, 
 *                                               descent_requirement); 
 *                       } 
 *                   } 
 *               } 
 * 
 *             output_results(iteration_number); 
 *           } 
 *         while ((iteration_number < max_iterations) && 
 *                (check_convergence(nonlinear_solution) == false)); 
 * 
 * @endcode
 * 
 * 在外循环结束时，我们必须更新屏障参数，为此我们使用以下公式。该函数的其余部分只是检查外循环的收敛条件，如果我们决定终止计算，就把最终的 "设计 "写成STL文件，用于3D打印，并输出一些时间信息。
 * 

 * 
 * 
 * @code
 *         const double barrier_size_multiplier = .8; 
 *         const double barrier_size_exponent   = 1.2; 
 * 
 *         barrier_size = 
 *           std::max(std::min(barrier_size * barrier_size_multiplier, 
 *                             std::pow(barrier_size, barrier_size_exponent)), 
 *                    min_barrier_size); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *     while (((barrier_size > min_barrier_size) || 
 *             (check_convergence(nonlinear_solution) == false)) && 
 *            (iteration_number < max_iterations)); 
 * 
 *     write_as_stl(); 
 *     timer.print_summary(); 
 *   } 
 * } // namespace SAND 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 

 * 
 * 余下的代码，即`main()`函数，和平常一样。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       SAND::SANDTopOpt<2> elastic_problem_2d; 
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
 *                 << "Aborting!" << std::endl;
 *       return 1;
 *     }
 *   return 0;
 * }
 * @endcode
examples/step-79/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="TestProblem"></a><h3>Test Problem</h3>上面使用的算法是针对一个传统的拓扑优化问题进行测试的，这个问题叫做Messerschmitt-Bolkow-Blohm Beam（MBB Beam）。


这个问题考虑的是在一个6个单位宽、1个单位高的矩形上可以建立的最佳二维结构。底部的角在 $y$ 方向用零Dirichlet边界条件固定住，通过强制执行Neumann边界条件在梁的顶部中心施加一个向下的力。边界的其余部分被允许移动，并且没有施加任何外力，这采取了零诺伊曼边界条件的形式。从本质上讲，我们提出了以下问题。我们应该如何设计一座桥，使桥的左下角和右下角的点在滚轮上，允许这些点在水平方向上移动，但不允许在垂直方向上移动，从而使响应于中心的垂直力的位移最小。

虽然领域的总体积是6个单位，但结构允许有3个单位的材料。由于问题的对称性，可以在一个宽为3、高为1的矩形上提出，方法是将原域切成两半，并沿切边在 $x$ 方向使用零迪里希特边界条件。也就是说，解决方案的对称性是一个很好的指标，表明程序正在按预期工作，所以我们在整个领域上解决问题，如下图所示。   @cite Bendse2004 

<div style="text-align:center;"> <img src="https://www.dealii.org/images/steps/developer/step-79.mbbgeometry.png" alt="MBB问题域和边界条件"> </div>


使用上面讨论的程序，我们找到了MBB梁的最小体积，解决方案的各个组成部分看起来如下。

<div class="onecolumn" style="width: 80%; text-align: center;"> <div> <img src="https://www.dealii.org/images/steps/developer/step-79.filtereddensity.png" alt="过滤的密度溶液"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step-79.unfiltereddensity.png" alt="未过滤的密度溶液"> </div> </div>


这些图片表明，我们在这里发现的情况与人们通常在关于该主题的其他出版物中看到的情况相一致  @cite Bendse2004  。也许更有趣的是，结果看起来像一座桁架桥（除了我们在桁架的顶部施加负载，而不是像真正的桁架桥那样在底部施加负载，类似于 "桥面桁架 "桥），这表明几个世纪以来一直用于桥梁建设的设计确实是基于我们现在可以证明在某种意义上是最佳的想法。




<a name="Possibilitiesforextensions"></a><h4>Possibilities for extensions</h4>


上面显示的结果花了大约75次迭代才找到，考虑到在每次迭代中解决大型线性系统的费用，这相当令人担忧。看一下演化过程，收敛确实有快速发生和缓慢发生的时候。我们认为这是由于在何时和如何减少边界值方面缺乏精确性，以及我们对优点函数的选择不够理想。在未来，用LOQO障碍更新代替单调还原，以及用马尔科夫滤波器代替优点函数，将大大减少必要的迭代次数。

障碍物的减少在收敛的中间阶段最为敏感，这是有问题的，因为我们似乎需要它快速减少，然后缓慢减少，然后又快速减少。

其次，这里使用的线性求解器只是基于SparseDirectUMFPACK类的稀疏直接求解器。这在小问题上效果还不错，但是上面详述的优化问题的表述有相当多的变量，因此线性问题不仅大，而且在许多行中有很多非零项，即使在总体上仍然比较粗糙的网格上。因此，解算器的时间在计算中占主导地位，需要采用更复杂的方法来解决线性系统。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-79.cc"
*/
