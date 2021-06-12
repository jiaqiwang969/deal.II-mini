/**
@page step_34 The step-34 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Irrotationalflow"> Irrotational flow </a>
        <li><a href="#Thenumericalapproximation">The numerical approximation</a>
        <li><a href="#Collocationboundaryelementmethod"> Collocation boundary element method </a>
        <li><a href="#Treatingthesingularintegrals"> Treating the singular integrals. </a>
        <li><a href="#Implementation">Implementation</a>
        <li><a href="#Testcase">Testcase</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Singleanddoublelayeroperatorkernels">Single and double layer operator kernels</a>
        <li><a href="#TheBEMProblemclass">The BEMProblem class</a>
      <ul>
        <li><a href="#BEMProblemBEMProblemandBEMProblemread_parameters">BEMProblem::BEMProblem and BEMProblem::read_parameters</a>
        <li><a href="#BEMProblemread_domain">BEMProblem::read_domain</a>
        <li><a href="#BEMProblemrefine_and_resize">BEMProblem::refine_and_resize</a>
        <li><a href="#BEMProblemassemble_system">BEMProblem::assemble_system</a>
        <li><a href="#BEMProblemsolve_system">BEMProblem::solve_system</a>
        <li><a href="#BEMProblemcompute_errors">BEMProblem::compute_errors</a>
        <li><a href="#BEMProblemcompute_exterior_solution">BEMProblem::compute_exterior_solution</a>
        <li><a href="#BEMProblemoutput_results">BEMProblem::output_results</a>
        <li><a href="#BEMProblemrun">BEMProblem::run</a>
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
examples/step-34/doc/intro.dox

 <br> 

<i>This program was contributed by Luca Heltai (thanks to Michael
Gratton for pointing out what the exact solution should have been in
the three dimensional case).  </i>

 @dealiiTutorialDOI{10.5281/zenodo.495473,https://zenodo.org/badge/DOI/10.5281/zenodo.495473.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="Irrotationalflow"></a><h3> Irrotational flow </h3> 无粘性流体经过一个物体（例如空气经过飞机机翼，或空气或水经过螺旋桨）的不可压缩运动，通常由流体力学的欧拉方程来模拟。


\f{align*}
  \frac{\partial }{\partial t}\mathbf{v} + (\mathbf{v}\cdot\nabla)\mathbf{v}
  &=


  -\frac{1}{\rho}\nabla p + \mathbf{g}
  \qquad &\text{in } \mathbb{R}^n \backslash \Omega
  \\
  \nabla \cdot \mathbf{v}&=0
  &\text{in } \mathbb{R}^n\backslash\Omega
\f}其中流体密度 $\rho$ 和由外力引起的加速度 $\mathbf{g}$ 是给定的，速度 $\mathbf{v}$ 和压力 $p$ 是未知数。这里 $\Omega$ 是一个封闭的有界区域，代表流体在其周围运动的体。

上述方程可以从纳维-斯托克斯方程推导出来，假设与压力梯度、惯性力和外力的影响相比，粘度造成的影响可以忽略不计。这是步骤22中讨论的斯托克斯方程的相反情况，斯托克斯方程是主导粘度的极限情况，即速度非常小，惯性力可以被忽略掉。另一方面，由于假定的不可压缩性，该方程不适合于非常高速的气体流动，在那里必须考虑到压缩性和气体的状态方程，导致气体动力学的欧拉方程，一个双曲系统。

在本教程程序中，我们将只考虑没有外力的静止流动。

\f{align*}
  (\mathbf{v}\cdot\nabla)\mathbf{v}
  &=


  -\frac{1}{\rho}\nabla p
  \qquad &\text{in } \mathbb{R}^n \backslash \Omega
  \\
  \nabla \cdot \mathbf{v}&=0
  &\text{in } \mathbb{R}^n\backslash\Omega
\f}


欧拉方程的解的唯一性通过添加边界条件得到保证

\f[
  \label{eq:boundary-conditions}
  \begin{aligned}
    \mathbf{n}\cdot\mathbf{v}& = 0 \qquad && \text{ on } \partial\Omega \\
    \mathbf{v}& = \mathbf{v}_\infty && \text{ when } |\mathbf{x}| \to \infty,
  \end{aligned}
\f]

这就是说，身体在我们的坐标系中是静止的，不具有渗透性，而流体在无限远处具有（恒定）速度 $\mathbf{v}_\infty$ 。另一种观点是，我们的坐标系随着身体移动，而背景流体在无限远处处于静止状态。注意，我们将法线 $\mathbf{n}$ 定义为域 $\Omega$ 的<i>outer</i>法线，它与积分域的外法线相反。

对于静止和非静止的流动，求解过程从求解第二个方程中的速度开始，然后代入第一个方程，以找到压力。静止欧拉方程的求解通常是为了了解给定的（可能是复杂的）几何形状在系统上强制执行规定运动时的行为。

这个过程的第一步是将参照系从一个与身体一起运动的坐标系改变为一个身体在一个无限大的静止流体中运动的坐标系。这可以通过引入一个新的速度 $\mathbf{\tilde{v}}=\mathbf{v}-\mathbf{v}_\infty$ 来表示，对于这个速度，我们发现同样的方程成立（因为 $\nabla\cdot
\mathbf{v}_\infty=0$ ），我们有边界条件

\f[
  \label{eq:boundary-conditions-tilde}
  \begin{aligned}
    \mathbf{n}\cdot\mathbf{\tilde{v}}& = -\mathbf{n}\cdot\mathbf{v}_\infty \qquad && \text{ on } \partial\Omega \\
    \mathbf{\tilde{v}}& = 0 && \text{ when } |\mathbf{x}| \to \infty,
  \end{aligned}
\f]

如果我们假设流体是无旋转的，即 $\nabla \times
\mathbf{v}=0$ 中的 $\mathbb{R}^n\backslash\Omega$ ，我们可以用标量函数的梯度来表示速度，从而也可以表示扰动速度。

\f[
  \mathbf{\tilde{v}}=\nabla\phi,
\f] ，因此上述欧拉方程的第二部分可以改写为未知数的同质拉普拉斯方程  $\phi$  。

\f{align*}
\label{laplace}
\Delta\phi &= 0 \qquad &&\text{in}\ \mathbb{R}^n\backslash\Omega,
	   \\
	   \mathbf{n}\cdot\nabla\phi &= -\mathbf{n}\cdot\mathbf{v}_\infty
	   && \text{on}\ \partial\Omega
\f}而动量方程还原为伯努利方程，该方程将压力  $p$  表示为势的函数  $\phi$  。

\f[
\frac{p}{\rho} +\frac{1}{2} | \nabla \phi |^2 = 0 \in \Omega.
\f]

因此，我们可以通过解决势的拉普拉斯方程来解决这个问题。  我们回顾一下，下列函数，称为拉普拉斯方程的基本解。

\f[ \begin{aligned}
\label{eq:3} G(\mathbf{y}-\mathbf{x}) = &


-\frac{1}{2\pi}\ln|\mathbf{y}-\mathbf{x}| \qquad && \text{for } n=2 \\
G(\mathbf{y}-\mathbf{x}) = &
\frac{1}{4\pi}\frac{1}{|\mathbf{y}-\mathbf{x}|}&& \text{for } n=3,
\end{aligned}
\f]

在分布意义上满足方程的要求。

\f[


-\Delta_y G(\mathbf{y}-\mathbf{x}) = \delta(\mathbf{y}-\mathbf{x}),
\f]

其中导数是在变量 $\mathbf{y}$ 中完成的。通过使用通常的格林同一性，我们的问题可以只写在边界上  $\partial\Omega = \Gamma$  。我们回顾一下第二个格林同位数的一般定义。

\f[\label{green}
  \int_{\omega}
  (-\Delta u)v\,dx + \int_{\partial\omega} \frac{\partial u}{\partial \tilde{\mathbf{n}} }v \,ds
  =
  \int_{\omega}
  (-\Delta v)u\,dx + \int_{\partial\omega} u\frac{\partial v}{\partial \tilde{\mathbf{n}}} \,ds,
\f]

其中 $\tilde{\mathbf{n}}$ 是 $\omega$ 的表面的法线，从积分域 $\omega$ 向外指向。

在我们的例子中，积分域是 $\mathbb{R}^n\backslash\Omega$ ，其边界是 $ \Gamma_\infty \cup
\Gamma$ ，其中无穷大的 "边界 "被定义为

\f[
\Gamma_\infty \dealcoloneq \lim_{r\to\infty} \partial B_r(0).
\f]

在我们的程序中，法线被定义为<i>outer</i>到域 $\Omega$ ，也就是说，它们实际上是<i>inner</i>到积分域，在定义各种积分时，需要注意法线的正确符号，即用 $-\mathbf{n}$ 代替 $\tilde{\mathbf{n}}$ 。

如果我们把 $u$ 和 $v$ 分别与 $\phi$ 的解和拉普拉斯方程的基本解代入格林%同，只要 $\mathbf{x}$ 被选在 $\mathbb{R}^n\backslash\Omega$ 区域，就可以得到。

\f[
  \phi(\mathbf{x}) -
  \int_{\Gamma\cup\Gamma_\infty}\frac{\partial G(\mathbf{y}-\mathbf{x})}{\partial \mathbf{n}_y}\phi(\mathbf{y})\,ds_y
  =


  -\int_{\Gamma\cup\Gamma_\infty}G(\mathbf{y}-\mathbf{x})\frac{\partial \phi}{\partial \mathbf{n}_y}(\mathbf{y})\,ds_y
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega
\f]

其中法线现在指向<i>inward</i>的积分域。

请注意，在上述方程中，我们也有 $\Gamma_\infty$ 处的边界部分的积分。利用我们问题的边界条件，我们有 $\nabla \phi$ 在无限远处为零（这简化了右侧 $\Gamma_\infty$ 上的积分）。

左手边出现的 $\Gamma_\infty$ 上的积分可以通过观察 $\nabla\phi=0$ 来处理，这意味着 $\phi$ 在无穷远处必然是常数。我们把它的值定义为 $\phi_\infty$  。  要证明这一点是很容易的

\f[


-\int_{\Gamma_\infty} \frac{\partial G(\mathbf{y}-\mathbf{x})}
{\partial \mathbf{n}_y}\phi_\infty \,ds_y =
\lim_{r\to\infty} \int_{\partial B_r(0)} \frac{\mathbf{r}}{r} \cdot \nabla G(\mathbf{y}-\mathbf{x})
\phi_\infty \,ds_y = -\phi_\infty.
\f]

利用这一结果，我们可以利用所谓的单层和双层势能算子，只在边界上 $\Gamma$ 还原上述方程。

\f[\label{integral}
  \phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty


  -\left(S \frac{\partial \phi}{\partial n_y}\right)(\mathbf{x})
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega.
\f]

(这些算子的名称来自于它们分别描述了 $\mathbb{R}^n$ 中由于沿表面的单一薄片电荷和由于沿表面的双片电荷和反电荷而产生的电动势)。

在我们的例子中，我们知道边界上 $\phi$ 的纽曼值： $\mathbf{n}\cdot\nabla\phi = -\mathbf{n}\cdot\mathbf{v}_\infty$  。因此。

\f[
  \phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty +
   \left(S[\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
   \qquad \forall\mathbf{x} \in \mathbb{R}^n\backslash\Omega.
\f] 如果我们对上述方程的 $\mathbf{x}$ 采取趋向于 $\Gamma$ 的极限，利用众所周知的单层和双层算子的特性，我们得到一个正好在 $\Omega$ 的边界 $\Gamma$ 的方程。

\f[\label{SD}
  \alpha(\mathbf{x})\phi(\mathbf{x}) - (D\phi)(\mathbf{x}) = \phi_\infty +
  \left(S [\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
  \quad \mathbf{x}\in \partial\Omega,
\f]

这就是我们要找的边界积分方程（BIE），其中量 $\alpha(\mathbf{x})$ 是点 $\mathbf{x}$ 看到积分域 $\mathbb{R}^n\backslash\Omega$ 的角度或实体角的分数。

特别是，在边界 $\mathbf{x}$ 是可微的（即光滑）的点上，我们有 $\alpha(\mathbf{x})=\frac 12$ ，但在边界有角或边的点上，数值可能会更小或更大。

代入单层和双层运算符，我们得到。

\f[
  \alpha(\mathbf{x}) \phi(\mathbf{x})
  + \frac{1}{2\pi}\int_{\partial \Omega}  \frac{
  (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }{ |\mathbf{y}-\mathbf{x}|^2 }
  \phi(\mathbf{y}) \,ds_y
  = \phi_\infty


    -\frac{1}{2\pi}\int_{\partial \Omega}  \ln|\mathbf{y}-\mathbf{x}| \, \mathbf{n}\cdot\mathbf{v_\infty}\,ds_y
\f]为二维流动和

\f[
  \alpha(\mathbf{x}) \phi(\mathbf{x})
   + \frac{1}{4\pi}\int_{\partial \Omega} \frac{ (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }{ |\mathbf{y}-\mathbf{x}|^3 }\phi(\mathbf{y})\,ds_y
  = \phi_\infty +
  \frac{1}{4\pi}\int_{\partial \Omega} \frac{1}{|\mathbf{y}-\mathbf{x}|} \, \mathbf{n}\cdot\mathbf{v_\infty}\,ds_y
\f]适用于三维流动，其中基本解的法向导数被写成了便于计算的形式。在任何一种情况下， $\phi$ 都是完全在边界上提出的积分方程的解，因为 $\mathbf{x},\mathbf{y}\in\partial\Omega$  。

注意，点 $\mathbf{x}$ 看到域 $\Omega$ 的角度（在2D中）或实体角（在3D中） $\alpha(\mathbf{x})$ 的分数可以用双层势本身定义。

\f[
\alpha(\mathbf{x}) \dealcoloneq 1 -
\frac{1}{2(n-1)\pi}\int_{\partial \Omega} \frac{ (\mathbf{y}-\mathbf{x})\cdot\mathbf{n}_y  }
{ |\mathbf{y}-\mathbf{x}|^{n} }\phi(\mathbf{y})\,ds_y = 1+
\int_{\partial \Omega} \frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y.
\f]

如果我们考虑到这样一个事实，即纯诺伊曼问题的解在一个任意常数 $c$ 以内都是已知的，这意味着，如果我们将诺伊曼数据设为零，那么任何常数 $\phi = \phi_\infty$ 都将是一个解。在边界积分方程中插入常数解和诺伊曼边界条件，我们有

@f{align*}
\alpha\left(\mathbf{x}\right)\phi\left(\mathbf{x}\right)
&=\int_{\Omega}\phi\left(\mathbf{y}\right)\delta\left(\mathbf{y}-\mathbf{x}\right)\, dy\\
\Rightarrow
\alpha\left(\mathbf{x}\right)\phi_\infty
&=\phi_\infty\int_{\Gamma\cup\Gamma_\infty}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
=\phi_\infty\left[\int_{\Gamma_\infty}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
+\int_{\Gamma}\frac{ \partial G(\mathbf{y}-\mathbf{x}) }{\partial \mathbf{n}_y} \, ds_y
\right]


@f}

在 $\Gamma_\infty$ 上的积分是统一的，见上文，所以除以常数 $\phi_\infty$ 就得到了上面 $\alpha(\mathbf{x})$ 的明确表达。

虽然这个示例程序实际上只关注边界积分方程的求解，但在一个现实的设置中，我们仍然需要对速度进行求解。为此，请注意，我们刚刚计算了 $\phi(\mathbf{x})$ 的所有 $\mathbf{x}\in\partial\Omega$ 。在下一步中，我们可以计算（如果我们愿意，可以分析）所有 $\mathbb{R}^n\backslash\Omega$ 中的解 $\phi(\mathbf{x})$  。为此，回顾一下，我们有

\f[
  \phi(\mathbf{x})
  =
  \phi_\infty +
  (D\phi)(\mathbf{x})
  +
  \left(S[\mathbf{n}\cdot\mathbf{v}_\infty]\right)(\mathbf{x})
  \qquad \forall\mathbf{x}\in \mathbb{R}^n\backslash\Omega.
\f]，现在我们有了右手边的所有东西（ $S$ 和 $D$ 是我们可以评估的积分，边界上的法线速度已经给出，边界上的 $\phi$ 我们刚刚计算了）。最后，我们就可以恢复速度为  $\mathbf{\tilde v}=\nabla \phi$  。

注意，对 $\mathbf{x} \in
\Omega$ 的上述公式的评估结果应该是零，因为狄拉克三角 $\delta(\mathbf{x})$ 在域 $\mathbb{R}^n\backslash\Omega$ 的积分根据定义总是零。

作为最后的测试，让我们验证这个速度是否确实满足静止流场的动量平衡方程，即对于某个（未知）压力 $p$ 和一个给定的常数 $\rho$ ， $\mathbf{v}\cdot\nabla\mathbf{v} = -\frac 1\rho \nabla p$ 中是否 $\mathbf{v}=\mathbf{\tilde
v}+\mathbf{v}_\infty=\nabla\phi+\mathbf{v}_\infty$ 。换句话说，我们想验证上面所说的伯努利定律确实成立。为了证明这一点，我们用这个方程的左手边等同于

@f{align*}
  \mathbf{v}\cdot\nabla\mathbf{v}
  &=
  [(\nabla\phi+\mathbf{v}_\infty)\cdot\nabla] (\nabla\phi+\mathbf{v}_\infty)
  \\
  &=
  [(\nabla\phi+\mathbf{v}_\infty)\cdot\nabla] (\nabla\phi)


@f}

其中我们使用了 $\mathbf{v}_\infty$ 是常数。我们想把这个表达式写成某个东西的梯度（记住 $\rho$ 是一个常数）。如果我们单独考虑方程的组成部分，下一步会更方便（对出现两次的指数求和是隐含的）。

@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  (\partial_j\phi+v_{\infty,j}) \partial_j \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]


  -
  \partial_j [(\partial_j\phi+v_{\infty,j})] \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]


@f}

因为  $\partial_j \partial_j\phi = \Delta \phi = 0$  和  $\textrm{div}
\ \mathbf{v}_\infty=0$  。下一个。

@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_j [(\partial_j\phi+v_{\infty,j}) \partial_i\phi]
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  \partial_j [v_{\infty,j} \partial_i\phi]
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  \partial_j [v_{\infty,j}] \partial_i\phi
  +
  v_{\infty,j} \partial_j \partial_i\phi
  \\
  &=
  \partial_j [(\partial_j\phi) (\partial_i\phi)]
  +
  v_{\infty,j} \partial_j \partial_i\phi
  \\
  &=
  \partial_i \partial_j [(\partial_j\phi) \phi]


  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  +
  \partial_i [v_{\infty,j} \partial_j \phi]


  -
  \partial_i [v_{\infty,j}] \partial_j \phi


@f}

同样，最后一项消失了，因为 $\mathbf{v}_\infty$ 是常数，我们可以将第一项和第三项合并为一项。

@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_i (\partial_j [(\partial_j\phi) \phi + v_{\infty,j} \partial_j \phi])


  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  \\
  &=
  \partial_i [(\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi]


  -
  \partial_j [\partial_i (\partial_j\phi) \phi]


@f}



我们现在只需要对最后一项再做一下按摩。使用乘积规则，我们得到

@f{align*}
  \partial_j [\partial_i (\partial_j\phi) \phi]
  &=
  \partial_i [\partial_j \partial_j\phi] \phi
  +
  \partial_i [\partial_j \phi] (\partial_j \phi).


@f}

这些项中的第一个是零（因为，同样，对 $j$ 的求和得到 $\Delta\phi$ ，它是零）。最后一项可以写成 $\frac 12
\partial_i [(\partial_j\phi)(\partial_j\phi)]$ ，它是理想的梯度形式。因此，我们现在可以最终说明

@f{align*}
  [\mathbf{v}\cdot\nabla\mathbf{v}]_i
  &=
  \partial_i (\partial_j [(\partial_j\phi) \phi + v_{\infty,j} \partial_j \phi])


  -
  \partial_j [\partial_i (\partial_j\phi) \phi]
  \\
  &=
  \partial_i
  \left[
    (\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi


    -
    \frac 12 (\partial_j\phi)(\partial_j\phi)
  \right],
  \\
  &=
  \partial_i
  \left[
    \frac 12 (\partial_j\phi)(\partial_j \phi) + v_{\infty,j} \partial_j \phi
  \right],


@f}

或以矢量形式。

@f[
  \mathbf{v}\cdot\nabla\mathbf{v}
  =
  \nabla
  \left[
    \frac 12 \mathbf{\tilde v}^2
    + \mathbf{v}_{\infty} \cdot \mathbf{\tilde v}
  \right],


@f]

或者换句话说。

@f[
  p
  =


  -\rho
  \left[
    \frac 12 \mathbf{\tilde v}^2
    + \mathbf{v}_{\infty} \cdot \mathbf{\tilde v}
  \right]
  =


  -\rho
  \left[
    \frac 12 \mathbf{v}^2


    -
    \frac 12 \mathbf{v}_{\infty}^2
  \right]
  .


@f]

因为压力只确定到一个常数（它只在方程中出现一个梯度），一个同样有效的定义是

@f[
  p
  =


  -\frac 12 \rho \mathbf{v}^2
  .


@f]

这正是上面提到的伯努利定律。




<a name="Thenumericalapproximation"></a><h3>The numerical approximation</h3>


边界积分方程（BIE）的数值近似通常被称为边界元素法或面板法（后者主要用于计算流体力学界）。以下测试问题的目的是解决具有诺伊曼边界条件的拉普拉斯方程的积分表述，分别使用二维和三维空间的圆和球体，沿途说明了允许人们使用deal.II库处理边界元素问题几乎与有限元问题一样容易的特点。

为此，让 $\mathcal{T}_h = \bigcup_i K_i$ 成为流形 $\Gamma = \partial \Omega$ 的一个细分，如果 $n=2$ 则为 $M$ 线段，如果 $n=3$ 则为 $M$  四边形。我们将称每个单独的线段或四边形为<i>element</i>或<i>cell</i>，与周围空间的维度 $n$ 无关。我们将有限维空间 $V_h$ 定义为

\f[
  \label{eq:definition-Vh}
  V_h \dealcoloneq \{ v \in C^0(\Gamma) \text{ s.t. } v|_{K_i} \in \mathcal{Q}^1(K_i),
  \forall i\},
\f]的基函数 $\psi_i(\mathbf{x})$ ，我们将使用通常的FE_Q有限元，但这次它被定义在一个一维的流形上（我们通过使用第二个模板参数，通常默认为等于第一个；这里，我们将在一个 <code>dim</code> 维的空间中创建对象 <code>FE_Q@<dim-1,dim@></code> to indicate that we have <code>dim-1</code> 维单元）。一个 $\phi_h$ 的元素 $V_h$ 是由其系数的向量 $\boldsymbol{\phi}$ 唯一识别的，也就是说。

\f[
  \label{eq:definition-of-element}
  \phi_h(\mathbf{x}) \dealcoloneq \phi_i \psi_i(\mathbf{x}), \qquad
  \boldsymbol{\phi} \dealcoloneq \{ \phi_i \},
\f]，其中求和隐含在重复索引上。请注意，我们可以在这里使用不连续的元素&mdash；事实上，没有真正的理由使用连续的元素，因为积分表述并不意味着我们的试验函数有任何导数，所以连续性是不必要的，而且在文献中经常只使用片断常数元素。

<a name="Collocationboundaryelementmethod"></a><h3> Collocation boundary element method </h3>


到目前为止，最常见的边界积分方程的近似方法是使用基于碰撞的边界元素方法。

这种方法要求在一些与系统未知数数量相等的定位点上评估边界积分方程。这些点的选择是一个微妙的问题，需要仔细研究。假设这些点暂时是已知的，并称它们为 $\mathbf x_i$ 和 $i=0...n\_dofs$  。

那么问题就变成了。给定基准点 $\mathbf{v}_\infty$ ，在 $V_h$ 中找到一个函数 $\phi_h$ ，使得以下 $n\_dofs$ 方程得到满足。

\f{align*}
    \alpha(\mathbf{x}_i) \phi_h(\mathbf{x}_i)


    - \int_{\Gamma_y} \frac{ \partial G(\mathbf{y}-\mathbf{x}_i)}{\partial\mathbf{n}_y }
    \phi_h(\mathbf{y}) \,ds_y =
    \int_{\Gamma_y} G(\mathbf{y}-\mathbf{x}_i) \,
    \mathbf{n}_y\cdot\mathbf{v_\infty} \,ds_y
    ,
\f}

其中数量 $\alpha(\mathbf{x}_i)$ 是点 $\mathbf{x}_i$ 看到域 $\Omega$ 的（实体）角度的分数，如上所述，我们设定 $\phi_\infty$ 为零。  如果支持点 $\mathbf{x}_i$ 选择得当，那么问题可以写成以下线性系统。

\f[
\label{eq:linear-system}
(\mathbf{A}+\mathbf{N})\boldsymbol\phi = \mathbf{b},
\f]

其中

\f[
\begin{aligned}
\mathbf{A}_{ij}&=
\alpha(\mathbf{x}_i) \psi_j(\mathbf{x}_i)
= 1+\int_\Gamma
\frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}\,ds_y
\psi_j(\mathbf{x}_i)
\\
\mathbf{N}_{ij}&= - \int_\Gamma
  \frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}
  \psi_j(\mathbf{y}) \,ds_y
\\
\mathbf{b}_i&= \int_\Gamma
   G(\mathbf{y}-\mathbf{x}_i)  \, \mathbf{n}_y\cdot\mathbf{v_\infty}
   ds_y.
\end{aligned}
\f]

从线性代数的角度来看，可能的最佳选择是使矩阵 $\mathbf{A}+\mathbf{N}$ 成为最对角线的主导。一个自然的选择是选择 $\mathbf{x}_i$ 搭配点作为节点基函数 $\psi_i(\mathbf{x})$ 的支持点。在这种情况下， $\psi_j(\mathbf{x}_i)=\delta_{ij}$  ，因此，矩阵 $\mathbf{A}$ 是对角线，其条目为

\f[
  \mathbf{A}_{ii}
  =
  1+\int_\Gamma
  \frac{\partial G(\mathbf{y}-\mathbf{x}_i)}{\partial \mathbf{n}_y}\,ds_y
  =
  1-\sum_j N_{ij},
\f]，其中我们使用了 $\sum_j \psi_j(\mathbf{y})=1$ 作为通常的拉格朗日元素。有了这样的选择，矩阵 $\mathbf{A}$ 、 $\mathbf{N}$ 和右手边 $\mathbf{b}$ 的条目的计算需要对三角形 $\mathcal{T}_h$ 元素的奇异积分进行评估。在这些情况下，所有的积分都是在参考简单域上进行的，也就是说，我们假设 $\mathcal{T}_h$ 的每个元素 $K_i$ 可以表示为参考边界元素 $\hat K \dealcoloneq [0,1]^{n-1}$ 的线性（二维）或双线性（三维）变换，并且我们在从实数元素 $K_i$ 到参考元素 $\hat K$ 的变量改变后执行积分。

<a name="Treatingthesingularintegrals"></a><h3> Treating the singular integrals. </h3>


在二维空间，没有必要计算系统矩阵的对角线元素 $\mathbf{N}_{ii}$ ，因为即使分母在 $\mathbf{x}=\mathbf{y}$ 时归零，分子也总是零，因为 $\mathbf{n}_y$ 和 $(\mathbf{y}-\mathbf{x})$ 是正交的。]和 $(\mathbf{y}-\mathbf{x})$ 是正交的（在我们对 $\Omega$ 边界的多边形近似上），唯一的奇异积分出现在对 $\mathbf{b}_i$ 的第i个元素的计算上。

\f[
  \frac{1}{\pi}
  \int_{K_i}
  \ln|\mathbf{y}-\mathbf{x}_i| \, \mathbf{n}_y\cdot\mathbf{v_\infty} \,ds_y.
\f]

这可以通过QGaussLogR正交公式轻松处理。

同样，也可以使用QGaussOneOverR正交公式来进行三维空间的奇异积分。有兴趣的读者可以在其文档中找到关于这些正交规则如何工作的详细解释。

得到的矩阵 $\mathbf{A}+\mathbf{N}$ 是完整的。根据其大小，使用直接求解器或迭代求解器可能会很方便。为了这个例子代码的目的，我们选择只使用一个迭代求解器，而不提供任何预处理程序。

如果这是一个生产代码，而不是一个原理演示，有一些技术可以用来不存储完整的矩阵，而只存储那些大的和/或相关的条目。在边界元素方法的文献中，有大量的方法可以确定哪些元素是重要的，哪些是不重要的，从而使这些矩阵的表示方法明显稀疏，也有利于快速评估向量和矩阵之间的标量积。这不是本程序的目标，我们把它留给更复杂的实现。




<a name="Implementation"></a><h3>Implementation</h3>


实现起来相当直接。在以前的教程程序中都没有用到的主要一点是，deal.II中的大多数类不仅在维度上有模板，而且实际上在我们提出微分方程的流形的维度以及这个流形嵌入的空间的维度上也有模板。默认情况下，第二个模板参数等于第一个，这意味着我们要在二维空间的一个二维区域上求解。在这种情况下，要使用的三角化类是 <code>Triangulation@<2@></code> ，这相当于写成 <code>Triangulation@<2,2@></code> 的方式。

然而，事实并非如此：在目前的例子中，我们想在球体表面求解，这是一个嵌入三维空间的二维流形。因此，正确的类将是 <code>Triangulation@<2,3@></code> ，相应地，我们将使用 <code>DoFHandler@<2,3@></code> 作为DoF处理类， <code>FE_Q@<2,3@></code> 作为有限元。

关于如何处理生活在弯曲流形上的事物的一些进一步细节，可以在报告<a target="_top"
href="http://www.dealii.org/reports/codimension-one/desimone-heltai-manigrasso.pdf"><i>Tools
for the Solution of PDEs Defined on Curved Manifolds with the deal.II
Library</i><i>Tools
for the Solution of PDEs Defined on Curved Manifolds with the deal.II
Library</i> by A. DeSimone, L. Heltai, C. Manigrasso</a>中找到。此外，Step-38教程程序将我们在这里展示的内容扩展到了流形上提出的方程不是积分算子而实际上涉及导数的情况。




<a name="Testcase"></a><h3>Testcase</h3>


我们要解决的测试案例是一个圆形（2D）或球形（3D）的障碍物。这些几何体的网格将从当前目录下的文件中读入，然后一个SphericalManifold类型的对象将被附加到三角形上，以允许网格细化，尊重离散的初始网格背后的连续几何。

对于一个半径为 $a$ 的球体，以 $U$ 的速度在 $x$ 方向平移，势为

@f{align*}
\phi = -\frac{1}{2}U \left(\frac{a}{r}\right)3 r \cos\theta


@f}

例如，见J. N. Newman, <i>Marine Hydrodynamics</i>, 1977, pp.对于单位速度和半径，并限制 $(x,y,z)$ 位于球体表面， $\phi = -x/2$  。在测试问题中，流量为 $(1,1,1)$  ，因此在球体表面上适当的精确解是上述解与沿 $y$ 和 $z$ 轴的类似解的叠加，即 $\phi =
\frac{1}{2}(x + y + z)$  。


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
 * 程序一开始就包括了一堆include文件，我们将在程序的各个部分使用这些文件。其中大部分在以前的教程中已经讨论过了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/smartpointer.h> 
 * #include <deal.II/base/convergence_table.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/quadrature_selector.h> 
 * #include <deal.II/base/parsed_function.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/solver_control.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/precondition.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_in.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/grid/manifold_lib.h> 
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
 * 
 * @endcode
 * 
 * 这里有一些我们需要的C++标准头文件。
 * 

 * 
 * 
 * @code
 * #include <cmath> 
 * #include <iostream> 
 * #include <fstream> 
 * #include <string> 
 * 
 * @endcode
 * 
 * 这个序言的最后部分是将dealii命名空间中的所有内容导入到这个程序中的所有内容中。
 * 

 * 
 * 
 * @code
 * namespace Step34 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Singleanddoublelayeroperatorkernels"></a> 
 * <h3>Single and double layer operator kernels</h3>
 * 

 * 
 * 首先，让我们定义一下边界积分方程的机制。
 * 

 * 
 * 以下两个函数是单层和双层势能核的实际计算，即  $G$  和  $\nabla G$  。只有当矢量 $R = \mathbf{y}-\mathbf{x}$ 不同于零时，它们才是定义良好的。
 * 

 * 
 * 
 * @code
 *   namespace LaplaceKernel 
 *   { 
 *     template <int dim> 
 *     double single_layer(const Tensor<1, dim> &R) 
 *     { 
 *       switch (dim) 
 *         { 
 *           case 2: 
 *             return (-std::log(R.norm()) / (2 * numbers::PI)); 
 * 
 *           case 3: 
 *             return (1. / (R.norm() * 4 * numbers::PI)); 
 * 
 *           default: 
 *             Assert(false, ExcInternalError()); 
 *             return 0.; 
 *         } 
 *     } 
 * 
 *     template <int dim> 
 *     Tensor<1, dim> double_layer(const Tensor<1, dim> &R) 
 *     { 
 *       switch (dim) 
 *         { 
 *           case 2: 
 *             return R / (-2 * numbers::PI * R.norm_square()); 
 *           case 3: 
 *             return R / (-4 * numbers::PI * R.norm_square() * R.norm()); 
 * 
 *           default: 
 *             Assert(false, ExcInternalError()); 
 *             return Tensor<1, dim>(); 
 *         } 
 *     } 
 *   } // namespace LaplaceKernel 
 * @endcode
 * 
 * 
 * <a name="TheBEMProblemclass"></a> 
 * <h3>The BEMProblem class</h3>
 * 

 * 
 * 边界元素方法代码的结构与有限元素代码的结构非常相似，所以这个类的成员函数与其他大多数教程程序的成员函数一样。特别是，现在你应该熟悉从外部文件中读取参数，以及将不同的任务分割成不同的模块。这同样适用于边界元素方法，我们不会对其进行过多的评论，只是对其中的差异进行评论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BEMProblem 
 *   { 
 *   public: 
 *     BEMProblem(const unsigned int fe_degree      = 1, 
 *                const unsigned int mapping_degree = 1); 
 * 
 *     void run(); 
 * 
 *   private: 
 *     void read_parameters(const std::string &filename); 
 * 
 *     void read_domain(); 
 * 
 *     void refine_and_resize(); 
 * 
 * @endcode
 * 
 * 我们在这里发现的唯一真正不同的函数是装配程序。我们以最通用的方式编写了这个函数，以便能够方便地推广到高阶方法和不同的基本解（例如斯托克斯或麦克斯韦）。
 * 

 * 
 * 最明显的区别是，最终的矩阵是完整的，而且我们在通常的单元格循环内有一个嵌套的循环，访问所有自由度的支持点。 此外，当支持点位于我们所访问的单元内时，我们所执行的积分就会变成单数。
 * 

 * 
 * 实际的结果是，我们有两套正交公式、有限元值和临时存储，一套用于标准积分，另一套用于奇异积分，在必要时使用。
 * 

 * 
 * 
 * @code
 *     void assemble_system(); 
 * 
 * @endcode
 * 
 * 对于这个问题的解决有两种选择。第一个是使用直接求解器，第二个是使用迭代求解器。我们选择了第二种方案。
 * 

 * 
 * 我们组装的矩阵不是对称的，我们选择使用GMRES方法；然而为边界元素方法构建一个有效的预处理程序并不是一个简单的问题。这里我们使用一个非预处理的GMRES求解器。迭代求解器的选项，如公差、最大迭代次数等，都是通过参数文件选择的。
 * 

 * 
 * 
 * @code
 *     void solve_system(); 
 * 
 * @endcode
 * 
 * 一旦我们得到了解决方案，我们将计算计算出的势的 $L^2$ 误差，以及实体角的近似值的 $L^\infty$ 误差。我们使用的网格是平滑曲线的近似值，因此计算出的角的分量或实体角的对角线矩阵  $\alpha(\mathbf{x})$  应该一直等于  $\frac 12$  。在这个例程中，我们输出势的误差和计算角度的近似值的误差。注意，后者的误差实际上不是计算角度的误差，而是衡量我们对球体和圆的近似程度。
 * 

 * 
 * 对角度的计算做一些实验，对于较简单的几何形状，可以得到非常准确的结果。为了验证这一点，你可以在read_domain()方法中注释掉tria.set_manifold(1, manifold)一行，并检查程序生成的alpha。通过删除这个调用，每当细化网格时，新的节点将沿着构成粗略网格的直线放置，而不是被拉到我们真正想要近似的表面。在三维案例中，球体的粗网格是从一个立方体开始得到的，得到的字母值正好是面的节点上的 $\frac 12$ ，边的节点上的 $\frac 34$ 和顶点的8个节点上的 $\frac 78$ 。
 * 

 * 
 * 
 * @code
 *     void compute_errors(const unsigned int cycle); 
 * 
 * @endcode
 * 
 * 一旦我们在一维领域得到了一个解决方案，我们就想把它插值到空间的其他部分。这可以通过在compute_exterior_solution()函数中再次进行解与核的卷积来实现。
 * 

 * 
 * 我们想绘制速度变量，也就是势解的梯度。势解只在边界上是已知的，但我们使用与基本解的卷积在标准的二维连续有限元空间上进行插值。外推解的梯度图将给我们提供我们想要的速度。
 * 

 * 
 * 除了外域上的解，我们还在output_results()函数中输出域的边界上的解，当然了。
 * 

 * 
 * 
 * @code
 *     void compute_exterior_solution(); 
 * 
 *     void output_results(const unsigned int cycle); 
 * 
 * @endcode
 * 
 * 为了实现不受维度限制的编程，我们对这个单一的函数进行了专业化处理，以提取整合单元内部的奇异核所需的奇异正交公式。
 * 

 * 
 * 
 * @code
 *     const Quadrature<dim - 1> &get_singular_quadrature( 
 *       const typename DoFHandler<dim - 1, dim>::active_cell_iterator &cell, 
 *       const unsigned int index) const; 
 * 
 * @endcode
 * 
 * 通常的deal.II类可以通过指定问题的 "二维 "来用于边界元素方法。这是通过将Triangulation, FiniteElement和DoFHandler的可选第二模板参数设置为嵌入空间的维度来实现的。在我们的例子中，我们生成了嵌入在二维或三维空间的一维或二维网格。
 * 

 * 
 * 可选参数默认等于第一个参数，并产生我们在之前所有例子中看到的通常的有限元类。
 * 

 * 
 * 该类的构造方式是允许任意的域（通过高阶映射）和有限元空间的逼近顺序。有限元空间和映射的顺序可以在该类的构造函数中选择。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim - 1, dim> tria; 
 *     FE_Q<dim - 1, dim>          fe; 
 *     DoFHandler<dim - 1, dim>    dof_handler; 
 *     MappingQ<dim - 1, dim>      mapping; 
 * 
 * @endcode
 * 
 * 在BEM方法中，生成的矩阵是密集的。根据问题的大小，最终的系统可能通过直接的LU分解来解决，或者通过迭代方法来解决。在这个例子中，我们使用了一个无条件的GMRES方法。为BEM方法建立一个预处理程序是不容易的，我们在此不做处理。
 * 

 * 
 * 
 * @code
 *     FullMatrix<double> system_matrix; 
 *     Vector<double>     system_rhs; 
 * 
 * @endcode
 * 
 * 接下来的两个变量将表示解决方案 $\phi$ 以及一个向量，它将保存 $\alpha(\mathbf x)$ 的值（从一个点 $\mathbf x$ 可见的 $\Omega$ 的部分）在我们形状函数的支持点。
 * 

 * 
 * 
 * @code
 *     Vector<double> phi; 
 *     Vector<double> alpha; 
 * 
 * @endcode
 * 
 * 收敛表是用来输出精确解和计算出的字母的误差的。
 * 

 * 
 * 
 * @code
 *     ConvergenceTable convergence_table; 
 * 
 * @endcode
 * 
 * 下面的变量是我们通过参数文件来填充的。 本例中我们使用的新对象是 Functions::ParsedFunction 对象和QuadratureSelector对象。
 * 

 * 
 * Functions::ParsedFunction 类允许我们通过参数文件方便快捷地定义新的函数对象，自定义的定义可以非常复杂（关于所有可用的选项，见该类的文档）。
 * 

 * 
 * 我们将使用QuadratureSelector类来分配正交对象，该类允许我们根据一个识别字符串和公式本身的可能程度来生成正交公式。我们用它来允许自定义选择标准积分的正交公式，并定义奇异正交规则的顺序。
 * 

 * 
 * 我们还定义了几个参数，这些参数是在我们想把解决方案扩展到整个领域的情况下使用的。
 * 

 * 
 * 
 * @code
 *     Functions::ParsedFunction<dim> wind; 
 *     Functions::ParsedFunction<dim> exact_solution; 
 * 
 *     unsigned int                         singular_quadrature_order; 
 *     std::shared_ptr<Quadrature<dim - 1>> quadrature; 
 * 
 *     SolverControl solver_control; 
 * 
 *     unsigned int n_cycles; 
 *     unsigned int external_refinement; 
 * 
 *     bool run_in_this_dimension; 
 *     bool extend_solution; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="BEMProblemBEMProblemandBEMProblemread_parameters"></a> 
 * <h4>BEMProblem::BEMProblem and BEMProblem::read_parameters</h4>
 * 

 * 
 * 构造函数初始化各种对象的方式与有限元程序（如  step-4  或  step-6  ）中的方式基本相同。这里唯一的新成分是ParsedFunction对象，它在构造时需要说明组件的数量。
 * 

 * 
 * 对于精确解来说，向量分量的数量是1，而且不需要任何操作，因为1是ParsedFunction对象的默认值。然而，风需要指定dim组件。注意，在为 Functions::ParsedFunction, 的表达式声明参数文件中的条目时，我们需要明确指定分量的数量，因为函数 Functions::ParsedFunction::declare_parameters 是静态的，对分量的数量没有了解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BEMProblem<dim>::BEMProblem(const unsigned int fe_degree, 
 *                               const unsigned int mapping_degree) 
 *     : fe(fe_degree) 
 *     , dof_handler(tria) 
 *     , mapping(mapping_degree, true) 
 *     , wind(dim) 
 *     , singular_quadrature_order(5) 
 *     , n_cycles(4) 
 *     , external_refinement(5) 
 *     , run_in_this_dimension(true) 
 *     , extend_solution(true) 
 *   {} 
 * 
 *   template <int dim> 
 *   void BEMProblem<dim>::read_parameters(const std::string &filename) 
 *   { 
 *     deallog << std::endl 
 *             << "Parsing parameter file " << filename << std::endl 
 *             << "for a " << dim << " dimensional simulation. " << std::endl; 
 * 
 *     ParameterHandler prm; 
 * 
 *     prm.declare_entry("Number of cycles", "4", Patterns::Integer()); 
 *     prm.declare_entry("External refinement", "5", Patterns::Integer()); 
 *     prm.declare_entry("Extend solution on the -2,2 box", 
 *                       "true", 
 *                       Patterns::Bool()); 
 *     prm.declare_entry("Run 2d simulation", "true", Patterns::Bool()); 
 *     prm.declare_entry("Run 3d simulation", "true", Patterns::Bool()); 
 * 
 *     prm.enter_subsection("Quadrature rules"); 
 *     { 
 *       prm.declare_entry( 
 *         "Quadrature type", 
 *         "gauss", 
 *         Patterns::Selection( 
 *           QuadratureSelector<(dim - 1)>::get_quadrature_names())); 
 *       prm.declare_entry("Quadrature order", "4", Patterns::Integer()); 
 *       prm.declare_entry("Singular quadrature order", "5", Patterns::Integer()); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 对于二维和三维，我们将默认的输入数据设置为：解为  $x+y$  或  $x+y+z$  。实际计算出的解在无穷大时的数值为零。在这种情况下，这与精确解相吻合，不需要额外的修正，但是你应该注意，我们任意设置了 $\phi_\infty$ ，而我们传递给程序的精确解需要在无穷远处有相同的值，才能正确计算出误差。
 * 

 * 
 * Functions::ParsedFunction 对象的使用是非常直接的。 Functions::ParsedFunction::declare_parameters 函数需要一个额外的整数参数，指定给定函数的分量数量。它的默认值是1。当相应的 Functions::ParsedFunction::parse_parameters 方法被调用时，调用对象必须有与这里定义的相同数量的组件，否则会产生异常。
 * 

 * 
 * 在声明条目时，我们同时声明了二维和三维的函数。然而只有二维的最终被解析。这使得我们对二维和三维问题都只需要一个参数文件。
 * 

 * 
 * 注意，从数学的角度来看，边界上的风函数应该满足条件 $\int_{\partial\Omega} \mathbf{v}\cdot \mathbf{n} d \Gamma = 0$  ，这样问题才会有解。如果不满足这个条件，那么就找不到解，求解器也就不会收敛。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Wind function 2d"); 
 *     { 
 *       Functions::ParsedFunction<2>::declare_parameters(prm, 2); 
 *       prm.set("Function expression", "1; 1"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Wind function 3d"); 
 *     { 
 *       Functions::ParsedFunction<3>::declare_parameters(prm, 3); 
 *       prm.set("Function expression", "1; 1; 1"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Exact solution 2d"); 
 *     { 
 *       Functions::ParsedFunction<2>::declare_parameters(prm); 
 *       prm.set("Function expression", "x+y"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Exact solution 3d"); 
 *     { 
 *       Functions::ParsedFunction<3>::declare_parameters(prm); 
 *       prm.set("Function expression", "x+y+z"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 在求解器部分，我们设置所有的SolverControl参数。然后，该对象将在solve_system()函数中被送入GMRES求解器。
 * 

 * 
 * 
 * @code
 *     prm.enter_subsection("Solver"); 
 *     SolverControl::declare_parameters(prm); 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 在向ParameterHandler对象声明了所有这些参数后，让我们读取一个输入文件，该文件将为这些参数提供其值。然后我们继续从ParameterHandler对象中提取这些值。
 * 

 * 
 * 
 * @code
 *     prm.parse_input(filename); 
 * 
 *     n_cycles            = prm.get_integer("Number of cycles"); 
 *     external_refinement = prm.get_integer("External refinement"); 
 *     extend_solution     = prm.get_bool("Extend solution on the -2,2 box"); 
 * 
 *     prm.enter_subsection("Quadrature rules"); 
 *     { 
 *       quadrature = std::shared_ptr<Quadrature<dim - 1>>( 
 *         new QuadratureSelector<dim - 1>(prm.get("Quadrature type"), 
 *                                         prm.get_integer("Quadrature order"))); 
 *       singular_quadrature_order = prm.get_integer("Singular quadrature order"); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Wind function " + std::to_string(dim) + "d"); 
 *     { 
 *       wind.parse_parameters(prm); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Exact solution " + std::to_string(dim) + "d"); 
 *     { 
 *       exact_solution.parse_parameters(prm); 
 *     } 
 *     prm.leave_subsection(); 
 * 
 *     prm.enter_subsection("Solver"); 
 *     solver_control.parse_parameters(prm); 
 *     prm.leave_subsection(); 
 * 
 * @endcode
 * 
 * 最后，这里是另一个如何在独立维度编程中使用参数文件的例子。 如果我们想关闭两个模拟中的一个，我们可以通过设置相应的 "运行2D模拟 "或 "运行3D模拟 "标志为假来实现。
 * 

 * 
 * 
 * @code
 *     run_in_this_dimension = 
 *       prm.get_bool("Run " + std::to_string(dim) + "d simulation"); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemread_domain"></a> 
 * <h4>BEMProblem::read_domain</h4>
 * 

 * 
 * 边界元素法三角剖分与（dim-1）维三角剖分基本相同，不同之处在于顶点属于（dim）维空间。
 * 

 * 
 * deal.II中支持的一些网格格式默认使用三维点来描述网格。这些格式与deal.II的边界元素方法功能兼容。特别是我们可以使用UCD或GMSH格式。在这两种情况下，我们必须特别注意网格的方向，因为与标准有限元的情况不同，这里没有进行重新排序或兼容性检查。 所有的网格都被认为是有方向性的，因为它们被嵌入到一个高维空间中。参见GridIn和Triangulation的文档，以进一步了解三角结构中单元的方向。在我们的例子中，网格的法线是外在于2D的圆或3D的球体。
 * 

 * 
 * 对边界元素网格进行适当细化所需要的另一个细节是对网格所逼近的流形的准确描述。对于标准有限元网格的边界，我们已经多次看到了这一点（例如在 step-5 和 step-6 中），这里的原理和用法是一样的，只是SphericalManifold类需要一个额外的模板参数来指定嵌入空间维度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::read_domain() 
 *   { 
 *     const Point<dim>                      center = Point<dim>(); 
 *     const SphericalManifold<dim - 1, dim> manifold(center); 
 * 
 *     std::ifstream in; 
 *     switch (dim) 
 *       { 
 *         case 2: 
 *           in.open("coarse_circle.inp"); 
 *           break; 
 * 
 *         case 3: 
 *           in.open("coarse_sphere.inp"); 
 *           break; 
 * 
 *         default: 
 *           Assert(false, ExcNotImplemented()); 
 *       } 
 * 
 *     GridIn<dim - 1, dim> gi; 
 *     gi.attach_triangulation(tria); 
 *     gi.read_ucd(in); 
 * 
 *     tria.set_all_manifold_ids(1); 
 * 
 * @endcode
 * 
 * 对  Triangulation::set_manifold  的调用复制了流形（通过  Manifold::clone()),  所以我们不需要担心对  <code>manifold</code>  的无效指针。
 * 

 * 
 * 
 * @code
 *     tria.set_manifold(1, manifold); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemrefine_and_resize"></a> 
 * <h4>BEMProblem::refine_and_resize</h4>
 * 

 * 
 * 这个函数对网格进行全局细化，分配自由度，并调整矩阵和向量的大小。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::refine_and_resize() 
 *   { 
 *     tria.refine_global(1); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     const unsigned int n_dofs = dof_handler.n_dofs(); 
 * 
 *     system_matrix.reinit(n_dofs, n_dofs); 
 * 
 *     system_rhs.reinit(n_dofs); 
 *     phi.reinit(n_dofs); 
 *     alpha.reinit(n_dofs); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemassemble_system"></a> 
 * <h4>BEMProblem::assemble_system</h4>
 * 

 * 
 * 下面是这个程序的主要功能，组装与边界积分方程相对应的矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::assemble_system() 
 *   { 
 * 
 * @endcode
 * 
 * 首先我们用正交公式初始化一个FEValues对象，用于在非奇异单元中进行内核积分。这个正交公式是通过参数文件选择的，并且需要相当精确，因为我们要积分的函数不是多项式函数。
 * 

 * 
 * 
 * @code
 *     FEValues<dim - 1, dim> fe_v(mapping, 
 *                                 fe, 
 *                                 *quadrature, 
 *                                 update_values | update_normal_vectors | 
 *                                   update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int n_q_points = fe_v.n_quadrature_points; 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices( 
 *       fe.n_dofs_per_cell()); 
 * 
 *     std::vector<Vector<double>> cell_wind(n_q_points, Vector<double>(dim)); 
 *     double                      normal_wind; 
 * 
 * @endcode
 * 
 * 与有限元方法不同的是，如果我们使用拼合边界元方法，那么在每个装配循环中，我们只装配与一个自由度（与支撑点 $i$ 相关的自由度）和当前单元之间的耦合信息。这是用fe.dofs_per_cell元素的向量完成的，然后将其分配到全局行的矩阵中  $i$  。以下对象将持有这些信息。
 * 

 * 
 * 
 * @code
 *     Vector<double> local_matrix_row_i(fe.n_dofs_per_cell()); 
 * 
 * @endcode
 * 
 * 索引  $i$  运行在拼合点上，这是  $i$  第三个基函数的支持点，而  $j$  运行在内部积分点上。
 * 

 * 
 * 我们构建一个支持点的向量，它将用于局部积分。
 * 

 * 
 * 
 * @code
 *     std::vector<Point<dim>> support_points(dof_handler.n_dofs()); 
 *     DoFTools::map_dofs_to_support_points<dim - 1, dim>(mapping, 
 *                                                        dof_handler, 
 *                                                        support_points); 
 * 
 * @endcode
 * 
 * 这样做之后，我们就可以开始对所有单元进行积分循环，首先初始化FEValues对象，得到正交点的 $\mathbf{\tilde v}$ 的值（这个向量场应该是常数，但更通用也无妨）。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_v.reinit(cell); 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
 *         const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 
 *         wind.vector_value_list(q_points, cell_wind); 
 * 
 * @endcode
 * 
 * 然后我们在当前单元上形成所有自由度的积分（注意，这包括不在当前单元上的自由度，这与通常的有限元积分有偏差）。如果其中一个局部自由度与支持点 $i$ 相同，我们需要执行的积分是单数。因此，在循环的开始，我们检查是否是这种情况，并存储哪一个是奇异指数。
 * 

 * 
 * 
 * @code
 *         for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
 *           { 
 *             local_matrix_row_i = 0; 
 * 
 *             bool         is_singular    = false; 
 *             unsigned int singular_index = numbers::invalid_unsigned_int; 
 * 
 *             for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
 *               if (local_dof_indices[j] == i) 
 *                 { 
 *                   singular_index = j; 
 *                   is_singular    = true; 
 *                   break; 
 *                 } 
 * 
 * @endcode
 * 
 * 然后我们进行积分。如果指数 $i$ 不是局部自由度之一，我们只需将单层项加到右边，将双层项加到矩阵中。
 * 

 * 
 * 
 * @code
 *             if (is_singular == false) 
 *               { 
 *                 for (unsigned int q = 0; q < n_q_points; ++q) 
 *                   { 
 *                     normal_wind = 0; 
 *                     for (unsigned int d = 0; d < dim; ++d) 
 *                       normal_wind += normals[q][d] * cell_wind[q](d); 
 * 
 *                     const Tensor<1, dim> R = q_points[q] - support_points[i]; 
 * 
 *                     system_rhs(i) += (LaplaceKernel::single_layer(R) * 
 *                                       normal_wind * fe_v.JxW(q)); 
 * 
 *                     for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
 * 
 *                       local_matrix_row_i(j) -= 
 *                         ((LaplaceKernel::double_layer(R) * normals[q]) * 
 *                          fe_v.shape_value(j, q) * fe_v.JxW(q)); 
 *                   } 
 *               } 
 *             else 
 *               { 
 * 
 * @endcode
 * 
 * 现在我们处理更微妙的情况。如果我们在这里，这意味着在 $j$ 索引上运行的单元包含support_point[i]。在这种情况下，单层和双层势都是单数，它们需要特殊处理。            
 * 每当在给定单元内进行积分时，就会使用一个特殊的正交公式，允许人们对参考单元上的奇异权重进行任意函数的积分。            
 * 正确的正交公式由get_singular_quadrature函数选择，下面将详细说明。
 * 

 * 
 * 
 * @code
 *                 Assert(singular_index != numbers::invalid_unsigned_int, 
 *                        ExcInternalError()); 
 * 
 *                 const Quadrature<dim - 1> &singular_quadrature = 
 *                   get_singular_quadrature(cell, singular_index); 
 * 
 *                 FEValues<dim - 1, dim> fe_v_singular( 
 *                   mapping, 
 *                   fe, 
 *                   singular_quadrature, 
 *                   update_jacobians | update_values | update_normal_vectors | 
 *                     update_quadrature_points); 
 * 
 *                 fe_v_singular.reinit(cell); 
 * 
 *                 std::vector<Vector<double>> singular_cell_wind( 
 *                   singular_quadrature.size(), Vector<double>(dim)); 
 * 
 *                 const std::vector<Tensor<1, dim>> &singular_normals = 
 *                   fe_v_singular.get_normal_vectors(); 
 *                 const std::vector<Point<dim>> &singular_q_points = 
 *                   fe_v_singular.get_quadrature_points(); 
 * 
 *                 wind.vector_value_list(singular_q_points, singular_cell_wind); 
 * 
 *                 for (unsigned int q = 0; q < singular_quadrature.size(); ++q) 
 *                   { 
 *                     const Tensor<1, dim> R = 
 *                       singular_q_points[q] - support_points[i]; 
 *                     double normal_wind = 0; 
 *                     for (unsigned int d = 0; d < dim; ++d) 
 *                       normal_wind += 
 *                         (singular_cell_wind[q](d) * singular_normals[q][d]); 
 * 
 *                     system_rhs(i) += (LaplaceKernel::single_layer(R) * 
 *                                       normal_wind * fe_v_singular.JxW(q)); 
 * 
 *                     for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
 *                       { 
 *                         local_matrix_row_i(j) -= 
 *                           ((LaplaceKernel::double_layer(R) * 
 *                             singular_normals[q]) * 
 *                            fe_v_singular.shape_value(j, q) * 
 *                            fe_v_singular.JxW(q)); 
 *                       } 
 *                   } 
 *               } 
 * 
 * @endcode
 * 
 * 最后，我们需要将当前单元格的贡献添加到全局矩阵中。
 * 

 * 
 * 
 * @code
 *             for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j) 
 *               system_matrix(i, local_dof_indices[j]) += local_matrix_row_i(j); 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 积分运算符的第二部分是术语  $\alpha(\mathbf{x}_i) \phi_j(\mathbf{x}_i)$  。由于我们使用的是配位方案， $\phi_j(\mathbf{x}_i)=\delta_{ij}$  而相应的矩阵是一个对角线的矩阵，其条目等于 $\alpha(\mathbf{x}_i)$  。
 * 

 * 
 * 计算这个实体角的对角矩阵的一个快速方法是使用诺伊曼矩阵本身。只需将该矩阵与一个元素都等于-1的向量相乘，就可以得到阿尔法角或实体角的对角线矩阵（见介绍中的公式）。然后将这个结果加回到系统矩阵对象上，得到矩阵的最终形式。
 * 

 * 
 * 
 * @code
 *     Vector<double> ones(dof_handler.n_dofs()); 
 *     ones.add(-1.); 
 * 
 *     system_matrix.vmult(alpha, ones); 
 *     alpha.add(1); 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
 *       system_matrix(i, i) += alpha(i); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemsolve_system"></a> 
 * <h4>BEMProblem::solve_system</h4>
 * 

 * 
 * 下一个函数简单地解决了线性系统。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::solve_system() 
 *   { 
 *     SolverGMRES<Vector<double>> solver(solver_control); 
 *     solver.solve(system_matrix, phi, system_rhs, PreconditionIdentity()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemcompute_errors"></a> 
 * <h4>BEMProblem::compute_errors</h4>
 * 

 * 
 * 误差的计算在其他所有的例子程序中都是完全一样的，我们就不做过多的评论。请注意，在有限元方法中使用的方法在这里也可以使用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::compute_errors(const unsigned int cycle) 
 *   { 
 *     Vector<float> difference_per_cell(tria.n_active_cells()); 
 *     VectorTools::integrate_difference(mapping, 
 *                                       dof_handler, 
 *                                       phi, 
 *                                       exact_solution, 
 *                                       difference_per_cell, 
 *                                       QGauss<(dim - 1)>(2 * fe.degree + 1), 
 *                                       VectorTools::L2_norm); 
 *     const double L2_error = 
 *       VectorTools::compute_global_error(tria, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 * @endcode
 * 
 * 可以直接使用 Vector::linfty_norm() 函数来计算α向量的误差，因为在每个节点上，该值应该是 $\frac 12$  。然后，所有的误差都会被输出并附加到我们的ConvergenceTable对象中，以便以后计算收敛率。
 * 

 * 
 * 
 * @code
 *     Vector<double> difference_per_node(alpha); 
 *     difference_per_node.add(-.5); 
 * 
 *     const double       alpha_error    = difference_per_node.linfty_norm(); 
 *     const unsigned int n_active_cells = tria.n_active_cells(); 
 *     const unsigned int n_dofs         = dof_handler.n_dofs(); 
 * 
 *     deallog << "Cycle " << cycle << ':' << std::endl 
 *             << "   Number of active cells:       " << n_active_cells 
 *             << std::endl 
 *             << "   Number of degrees of freedom: " << n_dofs << std::endl; 
 * 
 *     convergence_table.add_value("cycle", cycle); 
 *     convergence_table.add_value("cells", n_active_cells); 
 *     convergence_table.add_value("dofs", n_dofs); 
 *     convergence_table.add_value("L2(phi)", L2_error); 
 *     convergence_table.add_value("Linfty(alpha)", alpha_error); 
 *   } 
 * 
 * @endcode
 * 
 * 奇异积分需要仔细选择正交规则。特别是deal.II库提供了为对数奇异性（QGaussLog, QGaussLogR）以及1/R奇异性（QGaussOneOverR）量身定制的正交规则。
 * 

 * 
 * 奇异积分通常是通过构建具有奇异权重的加权正交公式得到的，因此可以写成
 * 

 * 
 * \f[ \int_K f(x) s(x) dx = \sum_{i=1}^N w_i f(q_i) \f]
 * 

 * 
 * 其中 $s(x)$ 是一个给定的奇点，权重和正交点 $w_i,q_i$ 是精心选择的，以使上述公式对某类函数 $f(x)$ 是一个等式。
 * 

 * 
 * 在我们目前看到的所有有限元例子中，正交点本身的权重（即函数  $s(x)$  ），总是不断地等于1。 对于奇异积分，我们有两个选择：我们可以使用上面的定义，从积分中剔除奇异性（即用特殊的正交规则对 $f(x)$ 进行积分），或者我们可以要求正交规则用 $s(q_i)$ 对权重 $w_i$ 进行 "标准化"。
 * 

 * 
 * \f[ \int_K f(x) s(x) dx = \int_K g(x) dx = \sum_{i=1}^N \frac{w_i}{s(q_i)} g(q_i) \f]
 * 

 * 
 * 我们通过QGaussLogR和QGaussOneOverR的 @p factor_out_singularity 参数来使用这第二种选择。
 * 

 * 
 * 这些积分有些微妙，特别是在二维空间，由于从实数到参考单元的转换，积分的变量是以转换的行列式为尺度的。
 * 

 * 
 * 在二维空间中，这个过程不仅会导致一个因子作为常数出现在整个积分上，而且还会导致一个需要评估的额外积分。
 * 

 * 
 * \f[ \int_0^1 f(x)\ln(x/\alpha) dx = \int_0^1 f(x)\ln(x) dx - \int_0^1  f(x) \ln(\alpha) dx.  \f]
 * 

 * 
 * 这个过程由QGaussLogR类的构造函数来处理，它增加了额外的正交点和权重，以考虑到积分的第二部分。
 * 

 * 
 * 类似的推理应该在三维情况下进行，因为奇异正交是在参考单元的半径 $r$ 的逆上定制的，而我们的奇异函数生活在实空间，然而在三维情况下一切都更简单，因为奇异性与变换的行列式成线性比例。这使我们可以只建立一次奇异的二维正交规则，并在所有单元格中重复使用。
 * 

 * 
 * 在一维的奇异积分中，这是不可能的，因为我们需要知道正交的缩放参数，而这个参数并不是先验的。这里，正交规则本身也取决于当前单元格的大小。出于这个原因，有必要为每个单数积分创建一个新的正交。
 * 

 * 
 * 不同的正交规则是在get_singular_quadrature中建立的，它专门用于dim=2和dim=3，它们在assemble_system函数中被检索。作为参数给出的索引是奇异点所在的单位支持点的索引。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   const Quadrature<2> &BEMProblem<3>::get_singular_quadrature( 
 *     const DoFHandler<2, 3>::active_cell_iterator &, 
 *     const unsigned int index) const 
 *   { 
 *     Assert(index < fe.n_dofs_per_cell(), 
 *            ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 
 * 
 *     static std::vector<QGaussOneOverR<2>> quadratures; 
 *     if (quadratures.size() == 0) 
 *       for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i) 
 *         quadratures.emplace_back(singular_quadrature_order, 
 *                                  fe.get_unit_support_points()[i], 
 *                                  true); 
 *     return quadratures[index]; 
 *   } 
 * 
 *   template <> 
 *   const Quadrature<1> &BEMProblem<2>::get_singular_quadrature( 
 *     const DoFHandler<1, 2>::active_cell_iterator &cell, 
 *     const unsigned int                            index) const 
 *   { 
 *     Assert(index < fe.n_dofs_per_cell(), 
 *            ExcIndexRange(0, fe.n_dofs_per_cell(), index)); 
 * 
 *     static Quadrature<1> *q_pointer = nullptr; 
 *     if (q_pointer) 
 *       delete q_pointer; 
 * 
 *     q_pointer = new QGaussLogR<1>(singular_quadrature_order, 
 *                                   fe.get_unit_support_points()[index], 
 *                                   1. / cell->measure(), 
 *                                   true); 
 *     return (*q_pointer); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="BEMProblemcompute_exterior_solution"></a> 
 * <h4>BEMProblem::compute_exterior_solution</h4>
 * 

 * 
 * 我们还想知道一些关于外域中电势 $\phi$ 的值：毕竟我们考虑边界积分问题的动机是想知道外域中的速度!
 * 

 * 
 * 为此，我们在此假设边界元素域包含在盒子 $[-2,2]^{\text{dim}}$ 中，我们用与基本解的卷积来推算这个盒子内的实际解。这方面的公式在介绍中已经给出。
 * 

 * 
 * 整个空间的解的重建是在一个连续的、尺寸为dim的有限元网格上完成的。这些都是常用的，我们不做进一步评论。在函数的最后，我们再次以通常的方式输出这个外部解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::compute_exterior_solution() 
 *   { 
 *     Triangulation<dim> external_tria; 
 *     GridGenerator::hyper_cube(external_tria, -2, 2); 
 * 
 *     FE_Q<dim>       external_fe(1); 
 *     DoFHandler<dim> external_dh(external_tria); 
 *     Vector<double>  external_phi; 
 * 
 *     external_tria.refine_global(external_refinement); 
 *     external_dh.distribute_dofs(external_fe); 
 *     external_phi.reinit(external_dh.n_dofs()); 
 * 
 *     FEValues<dim - 1, dim> fe_v(mapping, 
 *                                 fe, 
 *                                 *quadrature, 
 *                                 update_values | update_normal_vectors | 
 *                                   update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int n_q_points = fe_v.n_quadrature_points; 
 * 
 *     std::vector<types::global_dof_index> dofs(fe.n_dofs_per_cell()); 
 * 
 *     std::vector<double>         local_phi(n_q_points); 
 *     std::vector<double>         normal_wind(n_q_points); 
 *     std::vector<Vector<double>> local_wind(n_q_points, Vector<double>(dim)); 
 * 
 *     std::vector<Point<dim>> external_support_points(external_dh.n_dofs()); 
 *     DoFTools::map_dofs_to_support_points<dim>(StaticMappingQ1<dim>::mapping, 
 *                                               external_dh, 
 *                                               external_support_points); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_v.reinit(cell); 
 * 
 *         const std::vector<Point<dim>> &q_points = fe_v.get_quadrature_points(); 
 *         const std::vector<Tensor<1, dim>> &normals = fe_v.get_normal_vectors(); 
 * 
 *         cell->get_dof_indices(dofs); 
 *         fe_v.get_function_values(phi, local_phi); 
 * 
 *         wind.vector_value_list(q_points, local_wind); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             normal_wind[q] = 0; 
 *             for (unsigned int d = 0; d < dim; ++d) 
 *               normal_wind[q] += normals[q][d] * local_wind[q](d); 
 *           } 
 * 
 *         for (unsigned int i = 0; i < external_dh.n_dofs(); ++i) 
 *           for (unsigned int q = 0; q < n_q_points; ++q) 
 *             { 
 *               const Tensor<1, dim> R = q_points[q] - external_support_points[i]; 
 * 
 *               external_phi(i) += 
 *                 ((LaplaceKernel::single_layer(R) * normal_wind[q] + 
 *                   (LaplaceKernel::double_layer(R) * normals[q]) * 
 *                     local_phi[q]) * 
 *                  fe_v.JxW(q)); 
 *             } 
 *       } 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(external_dh); 
 *     data_out.add_data_vector(external_phi, "external_phi"); 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = std::to_string(dim) + "d_external.vtk"; 
 *     std::ofstream     file(filename); 
 * 
 *  
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemoutput_results"></a> 
 * <h4>BEMProblem::output_results</h4>
 * 

 * 
 * 输出我们的计算结果是一个相当机械的任务。这个函数的所有组成部分之前已经讨论过了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::output_results(const unsigned int cycle) 
 *   { 
 *     DataOut<dim - 1, dim> dataout; 
 * 
 *     dataout.attach_dof_handler(dof_handler); 
 *     dataout.add_data_vector(phi, "phi", DataOut<dim - 1, dim>::type_dof_data); 
 *     dataout.add_data_vector(alpha, 
 *                             "alpha", 
 *                             DataOut<dim - 1, dim>::type_dof_data); 
 *     dataout.build_patches(mapping, 
 *                           mapping.get_degree(), 
 *                           DataOut<dim - 1, dim>::curved_inner_cells); 
 * 
 *     const std::string filename = std::to_string(dim) + "d_boundary_solution_" + 
 *                                  std::to_string(cycle) + ".vtk"; 
 *     std::ofstream file(filename); 
 * 
 *     dataout.write_vtk(file); 
 * 
 *     if (cycle == n_cycles - 1) 
 *       { 
 *         convergence_table.set_precision("L2(phi)", 3); 
 *         convergence_table.set_precision("Linfty(alpha)", 3); 
 * 
 *         convergence_table.set_scientific("L2(phi)", true); 
 *         convergence_table.set_scientific("Linfty(alpha)", true); 
 * 
 *         convergence_table.evaluate_convergence_rates(
 *           "L2(phi)", ConvergenceTable::reduction_rate_log2); 
 *         convergence_table.evaluate_convergence_rates( 
 *           "Linfty(alpha)", ConvergenceTable::reduction_rate_log2); 
 *         deallog << std::endl; 
 *         convergence_table.write_text(std::cout); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BEMProblemrun"></a> 
 * <h4>BEMProblem::run</h4>
 * 

 * 
 * 这是最主要的功能。它应该是不言自明的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void BEMProblem<dim>::run() 
 *   { 
 *     read_parameters("parameters.prm"); 
 * 
 *     if (run_in_this_dimension == false) 
 *       { 
 *         deallog << "Run in dimension " << dim 
 *                 << " explicitly disabled in parameter file. " << std::endl; 
 *         return; 
 *       } 
 * 
 *     read_domain(); 
 * 
 *     for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
 *       { 
 *         refine_and_resize(); 
 *         assemble_system(); 
 *         solve_system(); 
 *         compute_errors(cycle); 
 *         output_results(cycle); 
 *       } 
 * 
 *     if (extend_solution == true) 
 *       compute_exterior_solution(); 
 *   } 
 * } // namespace Step34 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main() function</h3>
 * 

 * 
 * 这是本程序的主要功能。它与以前所有的教程程序完全一样。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step34; 
 * 
 *       const unsigned int degree         = 1; 
 *       const unsigned int mapping_degree = 1; 
 * 
 *       deallog.depth_console(3); 
 *       BEMProblem<2> laplace_problem_2d(degree, mapping_degree); 
 *       laplace_problem_2d.run(); 
 * 
 *       BEMProblem<3> laplace_problem_3d(degree, mapping_degree); 
 *       laplace_problem_3d.run(); 
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
examples/step-34/doc/results.dox



<a name="Results"></a><h1>Results</h1>


我们使用以下 <code>parameters.prm</code> 文件（也可以在所有其他源文件所在的目录中找到）运行该程序。

@verbatim
# Listing of Parameters
# ---------------------
set Extend solution on the -2,2 box = true
set External refinement             = 5
set Number of cycles                = 4
set Run 2d simulation               = true
set Run 3d simulation               = true



subsection Exact solution 2d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =


  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = x+y   # default: 0


  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,t
end



subsection Exact solution 3d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =


  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = .5*(x+y+z)   # default: 0


  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,z,t
end



subsection Quadrature rules
  set Quadrature order          = 4
  set Quadrature type           = gauss
  set Singular quadrature order = 5
end



subsection Solver
  set Log frequency = 1
  set Log history   = false
  set Log result    = true
  set Max steps     = 100
  set Tolerance     = 1.e-10
end



subsection Wind function 2d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =


  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = 1; 1  # default: 0; 0


  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,t
end



subsection Wind function 3d
  # Any constant used inside the function which is not a variable name.
  set Function constants  =


  # Separate vector valued expressions by ';' as ',' is used internally by the
  # function parser.
  set Function expression = 1; 1; 1 # default: 0; 0; 0


  # The name of the variables as they will be used in the function, separated
  # by ','.
  set Variable names      = x,y,z,t
end
@endverbatim



当我们运行该程序时，屏幕上会打印出以下内容。

@verbatim
DEAL::
DEAL::Parsing parameter file parameters.prm
DEAL::for a 2 dimensional simulation.
DEAL:GMRES::Starting value 2.21576
DEAL:GMRES::Convergence step 1 value 2.37635e-13
DEAL::Cycle 0:
DEAL::   Number of active cells:       20
DEAL::   Number of degrees of freedom: 20
DEAL:GMRES::Starting value 3.15543
DEAL:GMRES::Convergence step 1 value 2.89310e-13
DEAL::Cycle 1:
DEAL::   Number of active cells:       40
DEAL::   Number of degrees of freedom: 40
DEAL:GMRES::Starting value 4.46977
DEAL:GMRES::Convergence step 1 value 3.11815e-13
DEAL::Cycle 2:
DEAL::   Number of active cells:       80
DEAL::   Number of degrees of freedom: 80
DEAL:GMRES::Starting value 6.32373
DEAL:GMRES::Convergence step 1 value 3.22474e-13
DEAL::Cycle 3:
DEAL::   Number of active cells:       160
DEAL::   Number of degrees of freedom: 160
DEAL::
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    20   20 4.465e-02    - 5.000e-02    -
    1    40   40 1.081e-02 2.05 2.500e-02 1.00
    2    80   80 2.644e-03 2.03 1.250e-02 1.00
    3   160  160 6.529e-04 2.02 6.250e-03 1.00
DEAL::
DEAL::Parsing parameter file parameters.prm
DEAL::for a 3 dimensional simulation.
DEAL:GMRES::Starting value 2.84666
DEAL:GMRES::Convergence step 3 value 8.68638e-18
DEAL::Cycle 0:
DEAL::   Number of active cells:       24
DEAL::   Number of degrees of freedom: 26
DEAL:GMRES::Starting value 6.34288
DEAL:GMRES::Convergence step 5 value 1.38740e-11
DEAL::Cycle 1:
DEAL::   Number of active cells:       96
DEAL::   Number of degrees of freedom: 98
DEAL:GMRES::Starting value 12.9780
DEAL:GMRES::Convergence step 5 value 3.29225e-11
DEAL::Cycle 2:
DEAL::   Number of active cells:       384
DEAL::   Number of degrees of freedom: 386
DEAL:GMRES::Starting value 26.0874
DEAL:GMRES::Convergence step 6 value 1.47271e-12
DEAL::Cycle 3:
DEAL::   Number of active cells:       1536
DEAL::   Number of degrees of freedom: 1538
DEAL::
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    24   26 3.437e-01    - 2.327e-01    -
    1    96   98 9.794e-02 1.81 1.239e-01 0.91
    2   384  386 2.417e-02 2.02 6.319e-02 0.97
    3  1536 1538 5.876e-03 2.04 3.176e-02 0.99
@endverbatim



从2d中的收敛表可以看出，如果我们选择足够精确的正交公式，那么我们得到的 $\alpha(\mathbf{x})$ 的误差应该正好是元素数量的倒数。用N段大小相等的圆近似产生一个有N个面的正多边形，其角度正好是 $\pi-\frac {2\pi}{N}$ ，因此我们的误差应该正好是 $\frac 12 - (\frac 12 -\frac 1N) = \frac 1N$  。事实上，这是一个很好的指标，表明我们正在以适当的方式进行奇异积分。

势的近似 $\phi$ 的误差主要是由于域的近似。通过使用高阶映射可以得到更好的近似值。

如果我们修改main()函数，将fe_degree和mapping_degree设置为2，并提高参数文件中正交公式的顺序，我们得到以下二维模拟的收敛表

@verbatim
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    20   40 5.414e-05    - 2.306e-04    -
    1    40   80 3.623e-06 3.90 1.737e-05 3.73
    2    80  160 2.690e-07 3.75 1.253e-05 0.47
    3   160  320 2.916e-08 3.21 7.670e-06 0.71
@endverbatim



和

@verbatim
cycle cells dofs    L2(phi)     Linfty(alpha)
    0    24   98 3.770e-03    - 8.956e-03    -
    1    96  386 1.804e-04 4.39 1.182e-03 2.92
    2   384 1538 9.557e-06 4.24 1.499e-04 2.98
    3  1536 6146 6.617e-07 3.85 1.892e-05 2.99
@endverbatim



三维的情况下。我们可以看到，高阶映射的收敛结果要好得多，这主要是由于曲线几何的分辨率更高。请注意，在自由度相同的情况下，例如在三维模拟中Q1案例的第3步和Q2案例的第2步，误差大约低三个数量级。

运行这些计算的结果是一堆输出文件，我们可以将其传递给我们选择的可视化程序。输出文件有两种：边界元素表面的势，以及扩展到内外域的势。在二维的情况下，这两个文件的组合看起来像

 <img src="https://www.dealii.org/images/steps/developer/step-34_2d.png" alt=""> 

而在三维空间中，我们首先显示的是表面上的电位，同时还有一个等高线图。

 <img src="https://www.dealii.org/images/steps/developer/step-34_3d.png" alt=""> 

然后是潜力的外部等高线图，不透明度设置为25%。

 <img src="https://www.dealii.org/images/steps/developer/step-34_3d-2.png" alt=""> 


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这是第一个考虑解决嵌入高维空间的曲面上定义的方程的教程程序。但这里讨论的方程相对简单，因为它只涉及一个积分算子，而不涉及在曲面上更难定义的导数。step-38教程程序考虑了这类问题并提供了必要的工具。

从实际角度来看，这里使用的边界元素方法（BEM）有两个瓶颈。首先是组装矩阵的成本是*二次方的未知数，即 ${\cal O}(N^2)$ ，其中 $N$ 是未知数的总数。通过查看`assemble_system()`函数可以看出这一点，它有这样的结构。

@code
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        ...


        for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
          ...
@endcode

这里，第一个循环走过了所有的单元格（一个系数 $N$ ），而内循环则贡献了另一个系数 $N$  。

这必须与*局部*微分算子的有限元方法进行对比。在那里，我们在所有单元上进行循环（ $N$ 的一个因子），并在每个单元上做一个与有多少单元或未知数无关的工作。这显然是一个瓶颈。

第二个瓶颈是系统矩阵是密集的（即是FullMatrix类型），因为每个自由度都与其他自由度相耦合。正如上面所指出的，仅仅*计算*这个带有 $N^2$ 非零项的矩阵必然需要至少 ${\cal O}(N^2)$ 次操作，但值得指出的是，仅仅做一个矩阵-向量乘积也需要这么多操作。如果用于解决线性系统的GMRES方法需要的迭代次数随着问题的大小而增长，就像通常的情况一样，那么解决线性系统需要的运算次数甚至比 ${\cal O}(N^2)$ 还要快。

"真实 "边界元素方法通过确定矩阵的哪些条目将是小的，因此可以被忽略的策略来解决这些问题（当然，代价是引入额外的误差）。这可以通过认识到矩阵项随着自由度 $i$ 和 $j$ 定义的位置之间的（物理）距离衰减而实现。这可以在快速多极法（FMM）等方法中得到利用，这些方法可以控制哪些矩阵项必须被存储和计算以达到一定的精度，并且--如果做得好的话--导致方法中线性系统的装配和求解都需要少于 ${\cal O}(N^2)$ 的操作。

实施这些方法显然为扩展目前的计划提供了机会。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-34.cc"
*/
