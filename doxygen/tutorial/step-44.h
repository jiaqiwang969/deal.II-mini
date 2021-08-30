/**
@page step_44 The step-44 tutorial program
This tutorial depends on step-18.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Listofreferences">List of references</a>
        <li><a href="#Notation"> Notation </a>
        <li><a href="#Kinematics">Kinematics</a>
        <li><a href="#Kinetics">Kinetics</a>
        <li><a href="#Pushforwardandpullbackoperators"> Push-forward and pull-back operators </a>
        <li><a href="#Hyperelasticmaterials">Hyperelastic materials</a>
      <ul>
        <li><a href="#NeoHookeanmaterials"> Neo-Hookean materials </a>
      </ul>
        <li><a href="#Elasticitytensors">Elasticity tensors</a>
        <li><a href="#Principleofstationarypotentialenergyandthethreefieldformulation">Principle of stationary potential energy and the three-field formulation</a>
        <li><a href="#Discretizationofgoverningequations"> Discretization of governing equations </a>
        <li><a href="#Thematerialclass"> The material class </a>
        <li><a href="#Numericalexample"> Numerical example </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Runtimeparameters">Run-time parameters</a>
      <ul>
        <li><a href="#FiniteElementsystem">Finite Element system</a>
        <li><a href="#Geometry">Geometry</a>
        <li><a href="#Materials">Materials</a>
        <li><a href="#Linearsolver">Linear solver</a>
        <li><a href="#Nonlinearsolver">Nonlinear solver</a>
        <li><a href="#Time">Time</a>
        <li><a href="#Allparameters">All parameters</a>
      </ul>
        <li><a href="#Timeclass">Time class</a>
        <li><a href="#CompressibleneoHookeanmaterialwithinathreefieldformulation">Compressible neo-Hookean material within a three-field formulation</a>
        <li><a href="#Quadraturepointhistory">Quadrature point history</a>
        <li><a href="#Quasistaticquasiincompressiblefinitestrainsolid">Quasi-static quasi-incompressible finite-strain solid</a>
        <li><a href="#ImplementationofthecodeSolidcodeclass">Implementation of the <code>Solid</code> class</a>
      <ul>
        <li><a href="#Publicinterface">Public interface</a>
      </ul>
        <li><a href="#Privateinterface">Private interface</a>
      <ul>
        <li><a href="#Threadingbuildingblocksstructures">Threading-building-blocks structures</a>
        <li><a href="#Solidmake_grid">Solid::make_grid</a>
        <li><a href="#Solidsystem_setup">Solid::system_setup</a>
        <li><a href="#Solidsetup_qph">Solid::setup_qph</a>
        <li><a href="#Solidupdate_qph_incremental">Solid::update_qph_incremental</a>
        <li><a href="#Solidsolve_nonlinear_timestep">Solid::solve_nonlinear_timestep</a>
        <li><a href="#Solidprint_conv_headerandSolidprint_conv_footer">Solid::print_conv_header and Solid::print_conv_footer</a>
        <li><a href="#Solidget_error_dilation">Solid::get_error_dilation</a>
        <li><a href="#Solidget_error_residual">Solid::get_error_residual</a>
        <li><a href="#Solidget_error_update">Solid::get_error_update</a>
        <li><a href="#Solidget_total_solution">Solid::get_total_solution</a>
        <li><a href="#Solidassemble_system">Solid::assemble_system</a>
        <li><a href="#Solidmake_constraints">Solid::make_constraints</a>
        <li><a href="#Solidassemble_scmathsfmathbfK_widetildeJwidetildeJ">Solid::assemble_sc}  解决整个块系统有点问题，因为对 $\mathsf{\mathbf{K}}_{ \widetilde{J} \widetilde{J}</a>
        <li><a href="#Solidsolve_linear_system">Solid::solve_linear_system</a>
        <li><a href="#Solidoutput_results">Solid::output_results</a>
      </ul>
        <li><a href="#Mainfunction">Main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-44/doc/intro.dox

 <br> 

<i>This program was contributed by Jean-Paul Pelteret and Andrew McBride.
<br>
This material is based upon work supported by  the German Science Foundation (Deutsche
Forschungsgemeinschaft, DFG), grant STE 544/39-1,  and the National Research Foundation of South Africa.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.439772,https://zenodo.org/badge/DOI/10.5281/zenodo.439772.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


本教程的主题是非线性固体力学。经典的单场方法（例如见步骤18）不能正确描述准不可压缩材料的响应。响应过于僵硬；这种现象被称为锁定。锁定问题可以通过各种替代策略来规避。其中一个策略是三场公式。在这里，它被用来模拟各向同性连续体的三维、完全非线性（几何和材料）响应。材料响应被近似为超弹性。此外，所采用的三场公式对准不可压缩和可压缩材料都有效。

本报告的目的是为使用deal.II处理非线性固体力学的问题提供基础。线性问题在步骤8中得到了解决。在第18步中部分考虑了几何非线性问题的非标准的、超弹性的形式：使用了线性化构成关系的速率形式，问题域随着运动的进行而变化。围绕非线性运动学的重要概念在理论和实施中都没有。然而，第18步确实描述了许多关键概念，以便在deal.II的框架内实现弹性。

我们从非线性运动学的速成课程开始。为了简单起见，我们将注意力限制在准静态问题上。此后，我们介绍了各种关键的应力测量，并描述了构成模型。然后，在解释用于管理材料的类的结构之前，我们详细描述了三场公式。然后介绍了例子问题的设置。

 @note  本教程是针对三维空间的弹性问题而开发的（并在介绍中进行了描述）。  虽然空间维度可以在main()例程中改变，但需要注意的是。  一般来说，二维弹性问题只是作为三维问题的理想化而存在。  也就是说，它们要么是平面应变，要么是平面应力。  这些选择中的任何一个的假设都需要被一致地施加。  更多信息请参见步骤8的说明。

<a name="Listofreferences"></a><h3>List of references</h3>


这里实施的三场公式是由Simo等人（1985）开创的，被称为混合雅各布-压力公式。重要的相关贡献包括Simo和Taylor（1991）以及Miehe（1994）的贡献。这里采用的符号在很大程度上借鉴了Holzapfel（2001）对非线性固体力学理论方面的出色概述。Hughes (2000)对与不可压缩弹性（小应变时）有关的问题作了很好的概述。

<ol>  <li>  J.C. Simo, R.L. Taylor and K.S. Pister (1985), Variational and projection methods for the volume constraint in finite deformation elasto-plasticity,  <em>  Computer Methods in Applied Mechanics and Engineering  </em>  , <strong> 51</strong>, 1-3, 177-208. 		DOI: <a href="http://doi.org/10.1016/0045-7825(85)90033-7">10.1016/0045-7825(85)90033-7</a>;  <li>  J.C. Simo and R.L. Taylor (1991), Quasi-incompressible finite elasticity in principal stretches.Continuum basis and numerical algorithms,  <em>  Computer Methods in Applied Mechanics and Engineering  </em>  , <strong> 85 </strong>, 3, 273-310. 		DOI: <a href="http://doi.org/10.1016/0045-7825(91)90100-K">10.1016/0045-7825(91)90100-K</a>;  <li>  C. Miehe (1994), Aspects of the formulation and finite element implementation of large strain isotropic elasticity  <em>  International Journal for Numerical Methods in Engineering  </em>  <strong> 37 /strong>, 12, 1981-2004. 		DOI: <a href="http://doi.org/10.1002/nme.1620371202">10.1002/nme.1620371202</a>;  <li>  G.A. Holzapfel (2001), Nonlinear Solid Mechanics.A Continuum Approach for Engineering, John Wiley & Sons. 		ISBN: 0-471-82304-X;  <li>  T.J.R. Hughes (2000), The Finite Element Method:线性静态和动态有限元分析》，多佛。 		ISBN: 978-0486411811  </ol> .

<ol>  <li>  J-P. V. Pelteret, D. Davydov, A. McBride, D. K. Vu, and P. Steinmann (2016), 在一个耦合问题中使用这种三场公式的例子记录在<ol>  <li>  J-P.V. Pelteret, D. Davydov, A. McBride, D. K. Vu, and P. Steinmann (2016), Computational electro-and magneto-elasticity for quasi-incompressible media immersed in free space,  <em>  International Journal for Numerical Methods in Engineering  </em>  。 		DOI: <a href="http://doi.org/10.1002/nme.5254">10.1002/nme.5254</a>  </ol>  。

<a name="Notation"></a><h3> Notation </h3>


我们可以把四阶张量看作是将二阶张量（矩阵）映射到自己身上的线性算子，其方式与矩阵将向量映射到向量上一样。有各种四阶单位张量，在即将到来的介绍中会用到。四阶单位张量 $\mathcal{I}$ 和 $\overline{\mathcal{I}}$ 定义如下

@f[
	\mathbf{A} = \mathcal{I}:\mathbf{A}
		\qquad \text{and} \qquad
	\mathbf{A}^T = \overline{\mathcal{I}}:\mathbf{A} \, .


@f]

注意  $\mathcal{I} \neq \overline{\mathcal{I}}^T$  。此外，我们通过以下方式定义对称和偏斜对称的四阶单位张量

@f[
	\mathcal{S} \dealcoloneq \dfrac{1}{2}[\mathcal{I} + \overline{\mathcal{I}}]
		\qquad \text{and} \qquad
	\mathcal{W} \dealcoloneq \dfrac{1}{2}[\mathcal{I} - \overline{\mathcal{I}}] \, ,


@f]

以致于

@f[
	\dfrac{1}{2}[\mathbf{A} + \mathbf{A}^T] = \mathcal{S}:\mathbf{A}
		\qquad \text{and} \qquad
	\dfrac{1}{2}[\mathbf{A} - \mathbf{A}^T] = \mathcal{W}:\mathbf{A} \, .


@f]

identity_tensor()返回的四阶  <code>SymmetricTensor</code>  是  $\mathcal{S}$  。




<a name="Kinematics"></a><h3>Kinematics</h3>


让时间域表示为 $\mathbb{T} = [0,T_{\textrm{end}}]$  ，其中 $t \in \mathbb{T}$ 和 $T_{\textrm{end}}$ 是总的问题持续时间。考虑一个连续体，在时间 $t=0$ 占据参考配置 $\Omega_0$ 。参考配置中的%粒子由位置矢量 $\mathbf{X}$ 识别。身体在后来的时间 $t>0$ 的配置被称为当前配置，表示为 $\Omega$ ，粒子由矢量 $\mathbf{x}$ 识别。参考配置和当前配置之间的非线性映射，表示为  $\boldsymbol{\varphi}$  ，作用如下。

@f[
	\mathbf{x} = \boldsymbol{\varphi}(\mathbf{X},t) \, .


@f]

粒子的位移的材料描述被定义为

@f[
	\mathbf{U}(\mathbf{X},t) = \mathbf{x}(\mathbf{X},t) - \mathbf{X} \, .


@f]



变形梯度 $\mathbf{F}$ 被定义为运动的材料梯度。

@f[
	\mathbf{F}(\mathbf{X},t)
		\dealcoloneq \dfrac{\partial \boldsymbol{\varphi}(\mathbf{X},t)}{\partial \mathbf{X}}
		= \textrm{Grad}\ \mathbf{x}(\mathbf{X},t)
		= \mathbf{I} + \textrm{Grad}\ \mathbf{U} \, .


@f]

变形梯度 $J(\mathbf{X},t) \dealcoloneq \textrm{det}\ \mathbf{F}(\mathbf{X},t) > 0$ 的行列式在参考配置和当前配置中映射出相应的体积元素，分别表示为 $\textrm{d}V$ 和 $\textrm{d}v$  ，为

@f[
	\textrm{d}v = J(\mathbf{X},t)\; \textrm{d}V \, .


@f]



就空间和材料坐标而言，变形的两个重要度量是左和右Cauchy-Green张量，分别表示为 $\mathbf{b} \dealcoloneq \mathbf{F}\mathbf{F}^T$ 和 $\mathbf{C} \dealcoloneq \mathbf{F}^T\mathbf{F}$  。它们都是对称的和正定的。

格林-拉格朗日应变张量的定义为

@f[
	\mathbf{E} \dealcoloneq \frac{1}{2}[\mathbf{C} - \mathbf{I} ]
		= \underbrace{\frac{1}{2}[\textrm{Grad}^T \mathbf{U} +	\textrm{Grad}\mathbf{U}]}_{\boldsymbol{\varepsilon}}
			+ \frac{1}{2}[\textrm{Grad}^T\ \mathbf{U}][\textrm{Grad}\ \mathbf{U}] \, .


@f]

如果假定变形为无限小，那么右边的第二项就可以忽略， $\boldsymbol{\varepsilon}$ （线性化的应变张量）是应变张量的唯一组成部分。从问题的设置来看，这个假设在步骤18中是不成立的，这使得在该教程程序中使用线性化的 $\boldsymbol{\varepsilon}$ 作为应变度量值得怀疑。

为了处理材料在受到体积和剪切型变形时表现出的不同响应，我们考虑将变形梯度 $\mathbf{F}$ 和左Cauchy-Green张量 $\mathbf{b}$ 分解为体积变化（体积）和体积保持（等效）部分。

@f[
	\mathbf{F}
		= (J^{1/3}\mathbf{I})\overline{\mathbf{F}}
	\qquad \text{and} \qquad
	\mathbf{b}
        = (J^{2/3}\mathbf{I})\overline{\mathbf{F}}\,\overline{\mathbf{F}}^T
		=  (J^{2/3}\mathbf{I})\overline{\mathbf{b}} \, .


@f]

显然， $\textrm{det}\ \mathbf{F} = \textrm{det}\ (J^{1/3}\mathbf{I}) = J$  。

空间速度场被表示为 $\mathbf{v}(\mathbf{x},t)$  。空间速度场相对于空间坐标的导数给出空间速度梯度  $\mathbf{l}(\mathbf{x},t)$  ，即

@f[
	\mathbf{l}(\mathbf{x},t)
		\dealcoloneq \dfrac{\partial \mathbf{v}(\mathbf{x},t)}{\partial \mathbf{x}}
		= \textrm{grad}\ \mathbf{v}(\mathbf{x},t) \, ,


@f]

其中 $\textrm{grad} \{\bullet \}
= \frac{\partial \{ \bullet \} }{ \partial \mathbf{x}}
= \frac{\partial \{ \bullet \} }{ \partial \mathbf{X}}\frac{\partial \mathbf{X} }{ \partial \mathbf{x}}
= \textrm{Grad} \{ \bullet \} \mathbf{F}^{-1}$  。




<a name="Kinetics"></a><h3>Kinetics</h3>


考奇应力定理将作用在当前构型 $\mathbf{t}$ 的无穷小表面元素上的考奇牵引力 $\mathrm{d}a$ 等同于考奇应力张量 $\boldsymbol{\sigma}$ （一个空间量）与表面的外向单位法线 $\mathbf{n}$ 的积，即

@f[
	\mathbf{t}(\mathbf{x},t, \mathbf{n}) = \boldsymbol{\sigma}\mathbf{n} \, .


@f]

Cauchy应力是对称的。同样，作用于参考构型 $\mathbf{T}$ 中的无穷小表面元素的第一皮奥拉-基尔霍夫牵引力 $\mathrm{d}A$ 是第一皮奥拉-基尔霍夫应力张量 $\mathbf{P}$ （两点张量）与表面的外向单位法线 $\mathbf{N}$ 的乘积，为

@f[
	\mathbf{T}(\mathbf{X},t, \mathbf{N}) = \mathbf{P}\mathbf{N} \, .


@f]

Cauchy牵引力 $\mathbf{t}$ 和第一个Piola-Kirchhoff牵引力 $\mathbf{T}$ 的关系为

@f[
	\mathbf{t}\mathrm{d}a = \mathbf{T}\mathrm{d}A \, .


@f]

这可以用<a href="http://en.wikipedia.org/wiki/Finite_strain_theory">Nanson's formula</a>来证明。

第一个Piola-Kirchhoff应力张量与Cauchy应力的关系为

@f[
	\mathbf{P} = J \boldsymbol{\sigma}\mathbf{F}^{-T} \, .


@f]

进一步的重要应力测量是（空间）基尔霍夫应力  $\boldsymbol{\tau} = J \boldsymbol{\sigma}$  和（参考）第二Piola-Kirchhoff应力  $\mathbf{S} = {\mathbf{F}}^{-1} \boldsymbol{\tau} {\mathbf{F}}^{-T}$  。




<a name="Pushforwardandpullbackoperators"></a><h3> Push-forward and pull-back operators </h3>


前推和后拉运算符允许人们在材料和空间设置之间转换各种措施。这里使用的应力测量是逆变的，而应变测量是协变的。

二阶协变张量 $(\bullet)^{\text{cov}}$ 的前推和后拉操作分别由以下方法给出。

@f[
	\chi_{*}(\bullet)^{\text{cov}} \dealcoloneq \mathbf{F}^{-T} (\bullet)^{\text{cov}} \mathbf{F}^{-1}
	\qquad \text{and} \qquad
	\chi^{-1}_{*}(\bullet)^{\text{cov}} \dealcoloneq \mathbf{F}^{T} (\bullet)^{\text{cov}} \mathbf{F} \, .


@f]



二阶禁忌张量 $(\bullet)^{\text{con}}$ 的前推和后拉操作分别由以下方法给出。

@f[
	\chi_{*}(\bullet)^{\text{con}} \dealcoloneq \mathbf{F} (\bullet)^{\text{con}} \mathbf{F}^T
	\qquad \text{and} \qquad
	\chi^{-1}_{*}(\bullet)^{\text{con}} \dealcoloneq \mathbf{F}^{-1} (\bullet)^{\text{con}} \mathbf{F}^{-T} \, .


@f]

例如  $\boldsymbol{\tau} = \chi_{*}(\mathbf{S})$  。




<a name="Hyperelasticmaterials"></a><h3>Hyperelastic materials</h3>


超弹性材料的响应受亥姆霍兹自由能函数 $\Psi = \Psi(\mathbf{F}) = \Psi(\mathbf{C}) = \Psi(\mathbf{b})$ 的制约，该函数作为应力的势能。例如，如果Helmholtz自由能取决于右Cauchy-Green张量 $\mathbf{C}$ ，那么各向同性的超弹性响应为

@f[
	\mathbf{S}
		= 2 \dfrac{\partial \Psi(\mathbf{C})}{\partial \mathbf{C}} \, .


@f]

如果亥姆霍兹自由能取决于左Cauchy-Green张量 $\mathbf{b}$ ，那么各向同性的超弹性响应为

@f[
	\boldsymbol{\tau}
		= 2 \dfrac{\partial \Psi(\mathbf{b})}{\partial \mathbf{b}} \mathbf{b}
		=  2 \mathbf{b} \dfrac{\partial \Psi(\mathbf{b})}{\partial \mathbf{b}} \, .


@f]



根据变形梯度的乘法分解，亥姆霍兹自由能可以分解为

@f[
	\Psi(\mathbf{b}) = \Psi_{\text{vol}}(J) + \Psi_{\text{iso}}(\overline{\mathbf{b}}) \, .


@f]

同样，基尔霍夫应力可以分解为体积部分和等效部分 $\boldsymbol{\tau} = \boldsymbol{\tau}_{\text{vol}} + \boldsymbol{\tau}_{\text{iso}}$ ，其中。

@f{align*}
	\boldsymbol{\tau}_{\text{vol}} &=
		2 \mathbf{b} \dfrac{\partial \Psi_{\textrm{vol}}(J)}{\partial \mathbf{b}}
		\\
		&= p J\mathbf{I} \, ,
		\\
	\boldsymbol{\tau}_{\text{iso}} &=
		2 \mathbf{b} \dfrac{\partial \Psi_{\textrm{iso}} (\overline{\mathbf{b}})}{\partial \mathbf{b}}
		\\
		&= \underbrace{( \mathcal{I} - \dfrac{1}{3} \mathbf{I} \otimes \mathbf{I})}_{\mathbb{P}} : \overline{\boldsymbol{\tau}} \, ,


@f}

其中 $p \dealcoloneq \dfrac{\partial \Psi_{\text{vol}}(J)}{\partial J}$ 是压力响应。   $\mathbb{P}$ 是投影张量，它提供了欧拉环境下的偏差算子。虚构的基尔霍夫应力张量 $\overline{\boldsymbol{\tau}}$ 被定义为

@f[
	\overline{\boldsymbol{\tau}}
		\dealcoloneq 2 \overline{\mathbf{b}} \dfrac{\partial \Psi_{\textrm{iso}}(\overline{\mathbf{b}})}{\partial \overline{\mathbf{b}}} \, .


@f]






 @note  上述定义的压力响应与固体力学中广泛使用的压力定义不同，即 $p = - 1/3 \textrm{tr} \boldsymbol{\sigma} = - 1/3 J^{-1} \textrm{tr} \boldsymbol{\tau}$  。这里 $p$ 是静水压力。我们在本教程中使用压力响应（尽管我们把它称为压力）。

<a name="NeoHookeanmaterials"></a><h4> Neo-Hookean materials </h4>


与可压缩<a href="http://en.wikipedia.org/wiki/Neo-Hookean_solid">neo-Hookean material</a>相对应的亥姆霍兹自由能由以下公式给出

@f[
    \Psi \equiv
        \underbrace{\kappa [ \mathcal{G}(J) ] }_{\Psi_{\textrm{vol}}(J)}
        + \underbrace{\bigl[c_1 [ \overline{I}_1 - 3] \bigr]}_{\Psi_{\text{iso}}(\overline{\mathbf{b}})} \, ,


@f]

其中 $\kappa \dealcoloneq \lambda + 2/3 \mu$ 是体积模量（ $\lambda$ 和 $\mu$ 是Lam&eacute; 参数）和 $\overline{I}_1 \dealcoloneq \textrm{tr}\ \overline{\mathbf{b}}$  。函数 $\mathcal{G}(J)$ 被要求是严格凸的，并满足 $\mathcal{G}(1) = 0$ 等条件，进一步的细节见Holzapfel（2001）。在这项工作中  $\mathcal{G} \dealcoloneq \frac{1}{4} [ J^2 - 1 - 2\textrm{ln}J ]$  .

不可压缩性对所有运动施加了等效约束  $J=1$  。对应于不可压缩的新胡克材料的亥姆霍兹自由能由以下公式给出

@f[
    \Psi \equiv
        \underbrace{\bigl[ c_1 [ I_1 - 3] \bigr] }_{\Psi_{\textrm{iso}}(\mathbf{b})} \, ,


@f]

其中  $ I_1 \dealcoloneq \textrm{tr}\mathbf{b} $  。因此，通过从可压缩自由能中去除体积分量并执行  $J=1$  得到不可压缩响应。




<a name="Elasticitytensors"></a><h3>Elasticity tensors</h3>


我们将使用Newton-Raphson策略来解决非线性边界值问题。因此，我们将需要将构成关系线性化。

材料描述中的四阶弹性张量定义为

@f[
	\mathfrak{C}
		= 2\dfrac{\partial \mathbf{S}(\mathbf{C})}{\partial \mathbf{C}}
		= 4\dfrac{\partial^2 \Psi(\mathbf{C})}{\partial \mathbf{C} \partial \mathbf{C}} \, .


@f]

空间描述 $\mathfrak{c}$ 中的四阶弹性张量由 $\mathfrak{C}$ 的推演得到，为

@f[
	\mathfrak{c} = J^{-1} \chi_{*}(\mathfrak{C})
		\qquad \text{and thus} \qquad
	J\mathfrak{c} = 4 \mathbf{b} \dfrac{\partial^2 \Psi(\mathbf{b})} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}	\, .


@f]

四阶弹性张量（对于超弹性材料）同时拥有主要和次要的对称性。

四阶空间弹性张量可以写成以下解耦形式。

@f[
	\mathfrak{c} = \mathfrak{c}_{\text{vol}} + \mathfrak{c}_{\text{iso}} \, ,


@f]

其中

@f{align*}
	J \mathfrak{c}_{\text{vol}}
		&= 4 \mathbf{b} \dfrac{\partial^2 \Psi_{\text{vol}}(J)} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}
		\\
		&= J[\widehat{p}\, \mathbf{I} \otimes \mathbf{I} - 2p \mathcal{I}]
			\qquad \text{where} \qquad
		\widehat{p} \dealcoloneq p + \dfrac{\textrm{d} p}{\textrm{d}J} \, ,
		\\
	J \mathfrak{c}_{\text{iso}}
		&=  4 \mathbf{b} \dfrac{\partial^2 \Psi_{\text{iso}}(\overline{\mathbf{b}})} {\partial \mathbf{b} \partial \mathbf{b}} \mathbf{b}
		\\
		&= \mathbb{P} : \mathfrak{\overline{c}} : \mathbb{P}
			+ \dfrac{2}{3}[\overline{\boldsymbol{\tau}}:\mathbf{I}]\mathbb{P}


			- \dfrac{2}{3}[ \mathbf{I}\otimes\boldsymbol{\tau}_{\text{iso}}
				+ \boldsymbol{\tau}_{\text{iso}} \otimes \mathbf{I} ] \, ,


@f}

其中空间描述中的虚构弹性张量 $\overline{\mathfrak{c}}$ 被定义为

@f[
	\overline{\mathfrak{c}}
		= 4 \overline{\mathbf{b}} \dfrac{ \partial^2 \Psi_{\textrm{iso}}(\overline{\mathbf{b}})} {\partial \overline{\mathbf{b}} \partial \overline{\mathbf{b}}} \overline{\mathbf{b}} \, .


@f]



<a name="Principleofstationarypotentialenergyandthethreefieldformulation"></a><h3>Principle of stationary potential energy and the three-field formulation</h3>


系统的总势能 $\Pi$ 是内部和外部势能之和，分别表示为 $\Pi_{\textrm{int}}$ 和 $\Pi_{\textrm{ext}}$  。我们希望通过最小化势能找到平衡配置。

如上所述，我们采用了三场的表述。我们用 $\mathbf{\Xi} \dealcoloneq \{ \mathbf{u}, \widetilde{p}, \widetilde{J} \}$ 表示主要未知数的集合。独立运动学变量 $\widetilde{J}$ 作为对 $J$ 的约束进入公式，由拉格朗日乘数 $\widetilde{p}$ （压力，我们将看到）强制执行。

这里使用的三场变分原理由以下公式给出

@f[
	\Pi(\mathbf{\Xi}) \dealcoloneq \int_\Omega \bigl[
		\Psi_{\textrm{vol}}(\widetilde{J})
		+ \widetilde{p}\,[J(\mathbf{u}) - \widetilde{J}]
		+ \Psi_{\textrm{iso}}(\overline{\mathbf{b}}(\mathbf{u}))
		\bigr] \textrm{d}v
	+ 	\Pi_{\textrm{ext}} \, ,


@f]

其中外部电势的定义为

@f[
	\Pi_{\textrm{ext}}
		= - \int_\Omega \mathbf{b}^\text{p} \cdot \mathbf{u}~\textrm{d}v


			- \int_{\partial \Omega_{\sigma}} \mathbf{t}^\text{p} \cdot \mathbf{u}~\textrm{d}a \, .


@f]

当前配置 $\partial \Omega$ 的边界由两部分组成： $\partial \Omega = \partial \Omega_{\mathbf{u}} \cup \partial \Omega_{\sigma}$  ，其中 $\partial \Omega_{\mathbf{u}} \cap \partial \Omega_{\boldsymbol{\sigma}} = \emptyset$  。规定的Cauchy牵引力，表示为  $\mathbf{t}^\text{p}$  ，被应用于  $ \partial \Omega_{\boldsymbol{\sigma}}$  ，而运动被规定在边界的其余部分  $\partial \Omega_{\mathbf{u}}$  。每单位电流体积的体力表示为  $\mathbf{b}^\text{p}$  。




势的静止性如下

@f{align*}
	R(\mathbf\Xi;\delta \mathbf{\Xi})
		&= D_{\delta \mathbf{\Xi}}\Pi(\mathbf{\Xi})
		\\
		&= \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \mathbf{u}} \cdot \delta \mathbf{u}
			+ \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \widetilde{p}} \delta \widetilde{p}
			+ \dfrac{\partial \Pi(\mathbf{\Xi})}{\partial \widetilde{J}} \delta \tilde{J}
			\\
		&= \int_{\Omega_0}  \left[
			\textrm{grad}\ \delta\mathbf{u} : [ \underbrace{[\widetilde{p} J \mathbf{I}]}_{\equiv \boldsymbol{\tau}_{\textrm{vol}}}
            +  \boldsymbol{\tau}_{\textrm{iso}}]
			+ \delta \widetilde{p}\, [ J(\mathbf{u}) - \widetilde{J}]
			+ \delta \widetilde{J}\left[ \dfrac{\textrm{d} \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}}


            -\widetilde{p}\right]
			\right]~\textrm{d}V
			\\
		&\quad - \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\textrm{d}V


			- \int_{\partial \Omega_{0,\boldsymbol{\sigma}}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\textrm{d}A
			\\
		&=0 \, ,


@f}

对于所有虚拟位移 $\delta \mathbf{u} \in H^1(\Omega)$ ，受 $\delta \mathbf{u} = \mathbf{0}$ 对 $\partial \Omega_{\mathbf{u}}$ 的约束，以及所有虚拟压力 $\delta \widetilde{p} \in L^2(\Omega)$ 和虚拟膨胀 $\delta \widetilde{J} \in L^2(\Omega)$ 。

人们应该注意到，在三个场的表述中 $\boldsymbol{\tau}_{\textrm{vol}} \equiv \widetilde{p} J \mathbf{I}$ ，体积基尔霍夫应力的定义和随后的体积正切与超弹性材料一节中给出的一般形式略有不同，其中 $\boldsymbol{\tau}_{\textrm{vol}} \equiv p J\mathbf{I}$ 。这是因为压力 $\widetilde{p}$ 现在是一个主要的场，而不是一个构成性的派生量。我们需要仔细区分主要场和从构成关系中得到的场。

 @note  虽然变量都是用空间量来表示的，但积分的领域是初始配置。这种方法被称为  <em>  总拉格朗日公式  </em>  。在步骤18中给出的方法，其积分域是当前配置，可以称为  <em>  更新的拉格朗日公式  </em>  。这两种方法的各种优点在文献中被广泛讨论。然而，应该指出的是，它们是等同的。


与残留物相对应的欧拉-拉格朗日方程为：。

@f{align*}
	&\textrm{div}\ \boldsymbol{\sigma} + \mathbf{b}^\text{p} = \mathbf{0} && \textrm{[equilibrium]}
		\\
	&J(\mathbf{u}) = \widetilde{J} 		&& \textrm{[dilatation]}
		\\
	&\widetilde{p} = \dfrac{\textrm{d} \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}} && \textrm{[pressure]} \, .


@f}

第一个方程是空间设置中的（准静态）平衡方程。第二个是约束条件  $J(\mathbf{u}) = \widetilde{J}$  。第三个是压力的定义  $\widetilde{p}$  。

 @note 下面的简化单场推导（ $\mathbf{u}$ 是唯一的主变量）使我们清楚地知道如何将积分的极限转化为参考域。

@f{align*}
\int_{\Omega}\delta \mathbf{u} \cdot [\textrm{div}\ \boldsymbol{\sigma} + \mathbf{b}^\text{p}]~\mathrm{d}v
&=
\int_{\Omega} [-\mathrm{grad}\delta \mathbf{u}:\boldsymbol{\sigma} + \delta \mathbf{u} \cdot\mathbf{b}^\text{p}]~\mathrm{d}v
  + \int_{\partial \Omega} \delta \mathbf{u} \cdot \mathbf{t}^\text{p}~\mathrm{d}a \\
&=


- \int_{\Omega_0} \mathrm{grad}\delta \mathbf{u}:\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot J\mathbf{b}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_0} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \\
&=


- \int_{\Omega_0} \mathrm{grad}\delta \mathbf{u}:\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_{0,\sigma}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \\
&=


- \int_{\Omega_0} [\mathrm{grad}\delta\mathbf{u}]^{\text{sym}} :\boldsymbol{\tau}~\mathrm{d}V
+ \int_{\Omega_0} \delta \mathbf{u} \cdot \mathbf{B}^\text{p}~\mathrm{d}V
 + \int_{\partial \Omega_{0,\sigma}} \delta \mathbf{u} \cdot \mathbf{T}^\text{p}~\mathrm{d}A \, ,


@f}

其中 $[\mathrm{grad}\delta\mathbf{u}]^{\text{sym}} = 1/2[ \mathrm{grad}\delta\mathbf{u} + [\mathrm{grad}\delta\mathbf{u}]^T] $  。

我们将使用迭代牛顿-拉弗森方法来解决非线性剩余方程  $R$  。为了简单起见，我们假设死荷载，即荷载不因变形而改变。

在  $t_{\textrm{n}-1}$  的已知状态和  $t_{\textrm{n}}$  的当前未知状态之间的数量变化被表示为  $\varDelta \{ \bullet \} = { \{ \bullet \} }^{\textrm{n}} - { \{ \bullet \} }^{\textrm{n-1}}$  。在当前迭代 $\textrm{i}$ 的数量值表示为  ${ \{ \bullet \} }^{\textrm{n}}_{\textrm{i}} = { \{ \bullet \} }_{\textrm{i}}$  。迭代  $\textrm{i}$  和  $\textrm{i}+1$  之间的增量变化被表示为  $d \{ \bullet \} \dealcoloneq \{ \bullet \}_{\textrm{i}+1} - \{ \bullet \}_{\textrm{i}}$  。

假设系统的状态在某个迭代中是已知的  $\textrm{i}$  。用牛顿-拉弗森方法求解的非线性治理方程的线性化近似值是：找到  $d \mathbf{\Xi}$  ，以便

@f[
	R(\mathbf{\Xi}_{\mathsf{i}+1}) =
		R(\mathbf{\Xi}_{\mathsf{i}})
		+ D^2_{d \mathbf{\Xi}, \delta \mathbf{\Xi}} \Pi(\mathbf{\Xi_{\mathsf{i}}}) \cdot d \mathbf{\Xi} \equiv 0 \, ,


@f]

然后设置  $\mathbf{\Xi}_{\textrm{i}+1} = \mathbf{\Xi}_{\textrm{i}}
+ d \mathbf{\Xi}$  。切线由以下公式给出

@f[
	D^2_{d \mathbf{\Xi}, \delta \mathbf{\Xi}} \Pi( \mathbf{\Xi}_{\mathsf{i}} )
		= D_{d \mathbf{\Xi}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})
		=: K(\mathbf{\Xi}_{\mathsf{i}}; d \mathbf{\Xi}, \delta \mathbf{\Xi}) \, .


@f]

因此。

@f{align*}
 	K(\mathbf{\Xi}_{\mathsf{i}}; d \mathbf{\Xi}, \delta \mathbf{\Xi})
 		&=
 			D_{d \mathbf{u}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi}) \cdot d \mathbf{u}
 			\\
 				&\quad +
 			 	D_{d \widetilde{p}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})  d \widetilde{p}
 			 \\
 			 	&\quad +
 			  D_{d \widetilde{J}} R( \mathbf{\Xi}_{\mathsf{i}}; \delta \mathbf{\Xi})  d \widetilde{J} \, ,


@f}

其中

@f{align*}
	D_{d \mathbf{u}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=
 	\int_{\Omega_0} \bigl[ \textrm{grad}\ \delta \mathbf{u} :
 			\textrm{grad}\ d \mathbf{u} [\boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}]
 			+ \textrm{grad}\ \delta \mathbf{u} :[
             \underbrace{[\widetilde{p}J[\mathbf{I}\otimes\mathbf{I} - 2 \mathcal{I}]}_{\equiv J\mathfrak{c}_{\textrm{vol}}} +
             J\mathfrak{c}_{\textrm{iso}}] :\textrm{grad} d \mathbf{u}
 		\bigr]~\textrm{d}V \, ,
 		\\
 	&\quad + \int_{\Omega_0} \delta \widetilde{p} J \mathbf{I} : \textrm{grad}\ d \mathbf{u} ~\textrm{d}V
 	\\
 	D_{d \widetilde{p}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=
 	\int_{\Omega_0} \textrm{grad}\ \delta \mathbf{u} : J \mathbf{I} d \widetilde{p} ~\textrm{d}V


 		-  \int_{\Omega_0} \delta \widetilde{J} d \widetilde{p}  ~\textrm{d}V \, ,
 	\\
 	D_{d \widetilde{J}} R( \mathbf{\Xi}; \delta \mathbf{\Xi})
 	&=  -\int_{\Omega_0} \delta \widetilde{p} d \widetilde{J}~\textrm{d}V
 	 + \int_{\Omega_0} \delta \widetilde{J}  \dfrac{\textrm{d}^2 \Psi_{\textrm{vol}}(\widetilde{J})}{\textrm{d} \widetilde{J}\textrm{d}\widetilde{J}} d \widetilde{J} ~\textrm{d}V \, .


@f}



注意，以下条款被称为几何应力和材料对切线矩阵的贡献。

@f{align*}
& \int_{\Omega_0} \textrm{grad}\ \delta \mathbf{u} :
 			\textrm{grad}\ d \mathbf{u} [\boldsymbol{\tau}_{\textrm{iso}} +  \boldsymbol{\tau}_{\textrm{vol}}]~\textrm{d}V
 			&& \quad {[\textrm{Geometrical stress}]} \, ,
 		\\
& \int_{\Omega_0} \textrm{grad} \delta \mathbf{u} :
 			[J\mathfrak{c}_{\textrm{vol}} + J\mathfrak{c}_{\textrm{iso}}] :\textrm{grad}\ d \mathbf{u}
 		~\textrm{d}V
 		&& \quad {[\textrm{Material}]} \, .


@f}






<a name="Discretizationofgoverningequations"></a><h3> Discretization of governing equations </h3>


这里使用的三场公式对准不可压缩材料是有效的，即在 $\nu \rightarrow 0.5$ （其中 $\nu$ 是<a
href="http://en.wikipedia.org/wiki/Poisson's_ratio">Poisson's ratio</a>）的地方，要很好地选择 $\mathbf{u},~\widetilde{p}$ 和 $\widetilde{J}$ 的插值场。通常情况下，选择 $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$ 。这里 $DGPM$ 是FE_DGPMonomial类。一个流行的选择是 $Q_1 \times DGPM_0 \times DGPM_0$ ，它被称为平均扩张法（见Hughes（2000）的直观讨论）。这个代码可以容纳 $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$ 的表述。不连续的近似允许 $\widetilde{p}$ 和 $\widetilde{J}$ 被浓缩出来，并恢复了基于位移的经典方法。

对于完全不可压缩的材料 $\nu = 0.5$ 和三场公式仍将表现出锁定行为。这可以通过在自由能中引入一个额外的约束条件来克服，其形式为  $\int_{\Omega_0} \Lambda [ \widetilde{J} - 1]~\textrm{d}V$  。这里 $\Lambda$ 是一个拉格朗日乘数，用于强制执行等时约束条件。进一步的细节见Miehe (1994)。

线性化的问题可以写成

@f[
	\mathbf{\mathsf{K}}( \mathbf{\Xi}_{\textrm{i}}) d\mathbf{\Xi}
	=
	\mathbf{ \mathsf{F}}(\mathbf{\Xi}_{\textrm{i}})


@f]

其中

@f{align*}
		\underbrace{\begin{bmatrix}
			\mathbf{\mathsf{K}}_{uu}	&	\mathbf{\mathsf{K}}_{u\widetilde{p}}	& \mathbf{0}
			\\
			\mathbf{\mathsf{K}}_{\widetilde{p}u}	&	\mathbf{0}	&	\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}
			\\
			\mathbf{0}	& 	\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}		& \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\end{bmatrix}}_{\mathbf{\mathsf{K}}(\mathbf{\Xi}_{\textrm{i}})}
		\underbrace{\begin{bmatrix}
			d \mathbf{\mathsf{u}}\\
            d \widetilde{\mathbf{\mathsf{p}}} \\
            d \widetilde{\mathbf{\mathsf{J}}}
		\end{bmatrix}}_{d \mathbf{\Xi}}
        =
        \underbrace{\begin{bmatrix}


			-\mathbf{\mathsf{R}}_{u}(\mathbf{u}_{\textrm{i}}) \\


            -\mathbf{\mathsf{R}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}}) \\


           -\mathbf{\mathsf{R}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
		\end{bmatrix}}_{ -\mathbf{\mathsf{R}}(\mathbf{\Xi}_{\textrm{i}}) }
=
        \underbrace{\begin{bmatrix}
			\mathbf{\mathsf{F}}_{u}(\mathbf{u}_{\textrm{i}}) \\
            \mathbf{\mathsf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}}) \\
           \mathbf{\mathsf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
		\end{bmatrix}}_{ \mathbf{\mathsf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .


@f}



在配方中没有压力和膨胀（主要）变量的导数存在。因此，压力和膨胀的不连续有限元插值产生了 $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$ 、 $\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}$ 和 $\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}$ 的块对角矩阵。因此，我们可以很容易地表达每个单元上的场 $\widetilde{p}$ 和 $\widetilde{J}$ ，只需倒置一个局部矩阵并乘以局部右手。然后我们可以将结果插入其余的方程中，并恢复一个经典的基于位移的方法。为了在元素水平上凝结出压力和膨胀的贡献，我们需要以下结果。

@f{align*}
		d \widetilde{\mathbf{\mathsf{p}}}
		& = \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \bigl[
			 \mathbf{\mathsf{F}}_{\widetilde{J}}


			 - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}} d \widetilde{\mathbf{\mathsf{J}}} \bigr]
			\\
		d \widetilde{\mathbf{\mathsf{J}}}
		& = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
			\mathbf{\mathsf{F}}_{\widetilde{p}}


			- \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
			\bigr]
		\\
		 \Rightarrow d \widetilde{\mathbf{\mathsf{p}}}
		&=  \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}


		- \underbrace{\bigl[\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathbf{\mathsf{K}}}}\bigl[ \mathbf{\mathsf{F}}_{\widetilde{p}}


 		- \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}} \bigr]


@f}

因此

@f[
		\underbrace{\bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
		}_{\mathbf{\mathsf{K}}_{\textrm{con}}} d \mathbf{\mathsf{u}}
		=
        \underbrace{
		\Bigl[
		\mathbf{\mathsf{F}}_{u}


			- \mathbf{\mathsf{K}}_{u\widetilde{p}} \bigl[ \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}


			- \overline{\mathbf{\mathsf{K}}}\mathbf{\mathsf{F}}_{\widetilde{p}} \bigr]
		\Bigr]}_{\mathbf{\mathsf{F}}_{\textrm{con}}}


@f]

其中

@f[
		\overline{\overline{\mathbf{\mathsf{K}}}} \dealcoloneq
			\mathbf{\mathsf{K}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{K}}} \mathbf{\mathsf{K}}_{\widetilde{p}u} \, .


@f]

请注意，由于 $\widetilde{p}$ 和 $\widetilde{J}$ 选择的是元素层面的不连续，所有需要反转的矩阵都是在元素层面定义的。

构建各种贡献的程序如下。

- 构建  $\mathbf{\mathsf{K}}$  。

- 形成  $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}$  的元素，并存储在  $\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$  中的  $\mathbf{\mathsf{K}}$  。

- 形成 $\overline{\overline{\mathbf{\mathsf{K}}}}$ 并添加到 $\mathbf{\mathsf{K}}_{uu}$ ，得到 $\mathbf{\mathsf{K}}_{\textrm{con}}$ 。

- 修改后的系统矩阵被称为  ${\mathbf{\mathsf{K}}}_{\textrm{store}}$  。   也就是@f[
        \mathbf{\mathsf{K}}_{\textrm{store}}
\dealcoloneq
        \begin{bmatrix}
			\mathbf{\mathsf{K}}_{\textrm{con}}	&	\mathbf{\mathsf{K}}_{u\widetilde{p}}	& \mathbf{0}
			\\
			\mathbf{\mathsf{K}}_{\widetilde{p}u}	&	\mathbf{0}	&	\mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
			\\
			\mathbf{0}	& 	\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}		& \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
		\end{bmatrix} \, .
  @f] 。






<a name="Thematerialclass"></a><h3> The material class </h3>


一个好的面向对象的材料类的设计将有利于本教程扩展到广泛的材料类型。在本教程中，我们只有一个名为Material_Compressible_Neo_Hook_Three_Field的材料类。理想情况下，这个类会派生自超弹性材料（HyperelasticMaterial），而超弹性材料会派生自基类Material。这里使用的三场性质的表述也使问题复杂化。

三场公式的亥姆霍兹自由能函数为  $\Psi = \Psi_\text{vol}(\widetilde{J}) + \Psi_\text{iso}(\overline{\mathbf{b}})$  。Kirchhoff应力的等效部分 ${\boldsymbol{\tau}}_{\text{iso}}(\overline{\mathbf{b}})$ 与使用超弹性材料的单场公式得到的相同。然而，自由能的体积部分现在是一个主要变量的函数  $\widetilde{J}$  。因此，对于三场公式来说，基尔霍夫应力 ${\boldsymbol{\tau}}_{\text{vol}}$ 的体积部分的构成反应（和正切）并不像单场公式那样由超弹性构成法给出。我们可以将术语 $\boldsymbol{\tau}_{\textrm{vol}} \equiv \widetilde{p} J \mathbf{I}$ 标记为体积基尔霍夫应力，但压力 $\widetilde{p}$ 不是由自由能得出的；它是一个主场。

为了有一个灵活的方法，我们决定Material_Compressible_Neo_Hook_Three_Field仍然能够计算并返回一个体积Kirchhoff应力和正切。为了做到这一点，我们选择在与正交点相关的Material_Compressible_Neo_Hook_Three_Field类中存储插值的主域 $\widetilde{p}$ 和 $\widetilde{J}$ 。这个决定应该在以后的阶段，当教程扩展到考虑其他材料时，再重新审视。




<a name="Numericalexample"></a><h3> Numerical example </h3>


这里考虑的数值例子是一个压缩下的几乎不可压缩的块。这个基准问题取自

- S. Reese, P. Wriggers, B.D. Reddy (2000), A new locking-free brick element technique for large deformation problems in elasticity,  <em>  Computers and Structures  </em>  , <strong> 75</strong>, 291-304.   DOI:<a href="http://doi.org/10.1016/S0045-7949(99)00137-6">10.1016/S0045-7949(99)00137-6</a>。

   <img src="https://www.dealii.org/images/steps/developer/step-44.setup.png" alt=""> 

该材料是具有<a href="http://en.wikipedia.org/wiki/Shear_modulus">shear modulus</a> $\mu = 80.194e6$ 和 $\nu = 0.4999$ 的准不可压缩的新胡克式。对于这样一个材料特性的选择，传统的单场 $Q_1$ 方法将锁定。也就是说，响应会过于僵硬。初始和最终配置显示在上面的图片中。利用对称性，我们只求解四分之一的几何体（即一个尺寸为 $0.001$ 的立方体）。域的上表面的内四分之一受到 $p_0$ 的载荷。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 我们首先包括所有必要的deal.II头文件和一些C++相关的文件。它们已经在以前的教程程序中详细讨论过了，所以你只需要参考过去的教程就可以了。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/base/point.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/symmetric_tensor.h> 
 * #include <deal.II/base/tensor.h> 
 * #include <deal.II/base/timer.h> 
 * #include <deal.II/base/work_stream.h> 
 * #include <deal.II/dofs/dof_renumbering.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * @endcode
 * 
 * 这个标头为我们提供了在正交点存储数据的功能
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_point_data.h> 
 * 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/grid/grid_in.h> 
 * #include <deal.II/grid/tria.h> 
 * 
 * #include <deal.II/fe/fe_dgp_monomial.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_system.h> 
 * #include <deal.II/fe/fe_tools.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/fe/mapping_q_eulerian.h> 
 * 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/precondition_selector.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/solver_selector.h> 
 * #include <deal.II/lac/sparse_direct.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * @endcode
 * 
 * 这里是使用LinearOperator类所需的头文件。这些头文件也都被方便地打包到一个头文件中，即<deal.II/lac/linear_operator_tools.h>，但为了透明起见，我们在此列出那些特别需要的头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/linear_operator.h> 
 * #include <deal.II/lac/packaged_operation.h> 
 * 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * 
 * @endcode
 * 
 * 在这两个标题中定义的是一些与有限应变弹性有关的操作。第一个将帮助我们计算一些运动量，第二个提供一些标准的张量定义。
 * 

 * 
 * 
 * @code
 * #include <deal.II/physics/elasticity/kinematics.h> 
 * #include <deal.II/physics/elasticity/standard_tensors.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * 
 * @endcode
 * 
 * 然后，我们将所有与本教程程序有关的东西都放入一个自己的命名空间，并将所有deal.II的函数和类名导入其中。
 * 

 * 
 * 
 * @code
 * namespace Step44 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Runtimeparameters"></a> 
 * <h3>Run-time parameters</h3>
 * 

 * 
 * 有几个参数可以在代码中设置，所以我们设置了一个ParameterHandler对象，在运行时读入选择。
 * 

 * 
 * 
 * @code
 *   namespace Parameters 
 *   { 
 * @endcode
 * 
 * 
 * <a name="FiniteElementsystem"></a> 
 * <h4>Finite Element system</h4>
 * 

 * 
 * 正如介绍中提到的，对于位移 $\mathbf{u}$ 应该使用不同的阶次插值，而不是压力 $\widetilde{p}$ 和膨胀 $\widetilde{J}$ 。 选择 $\widetilde{p}$ 和 $\widetilde{J}$ 作为元素级的不连续（常数）函数，导致了平均扩张方法。不连续的近似允许 $\widetilde{p}$ 和 $\widetilde{J}$ 被浓缩出来，并恢复了基于位移的经典方法。这里我们指定用于近似解的多项式阶数。正交阶数应作相应调整。
 * 

 * 
 * 
 * @code
 *     struct FESystem 
 *     { 
 *       unsigned int poly_degree; 
 *       unsigned int quad_order; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void FESystem::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Finite element system"); 
 *       { 
 *         prm.declare_entry("Polynomial degree", 
 *                           "2", 
 *                           Patterns::Integer(0), 
 *                           "Displacement system polynomial order"); 
 * 
 *         prm.declare_entry("Quadrature order", 
 *                           "3", 
 *                           Patterns::Integer(0), 
 *                           "Gauss quadrature order"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void FESystem::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Finite element system"); 
 *       { 
 *         poly_degree = prm.get_integer("Polynomial degree"); 
 *         quad_order  = prm.get_integer("Quadrature order"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Geometry"></a> 
 * <h4>Geometry</h4>
 * 

 * 
 * 对问题的几何形状和应用的载荷进行调整。 由于这里模拟的问题比较特殊，所以可以将载荷比例改变为特定的数值，以便与文献中给出的结果进行比较。
 * 

 * 
 * 
 * @code
 *     struct Geometry 
 *     { 
 *       unsigned int global_refinement; 
 *       double       scale; 
 *       double       p_p0; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Geometry::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Geometry"); 
 *       { 
 *         prm.declare_entry("Global refinement", 
 *                           "2", 
 *                           Patterns::Integer(0), 
 *                           "Global refinement level"); 
 * 
 *         prm.declare_entry("Grid scale", 
 *                           "1e-3", 
 *                           Patterns::Double(0.0), 
 *                           "Global grid scaling factor"); 
 * 
 *         prm.declare_entry("Pressure ratio p/p0", 
 *                           "100", 
 *                           Patterns::Selection("20|40|60|80|100"), 
 *                           "Ratio of applied pressure to reference pressure"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Geometry::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Geometry"); 
 *       { 
 *         global_refinement = prm.get_integer("Global refinement"); 
 *         scale             = prm.get_double("Grid scale"); 
 *         p_p0              = prm.get_double("Pressure ratio p/p0"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Materials"></a> 
 * <h4>Materials</h4>
 * 

 * 
 * 我们还需要新胡克材料的剪切模量 $ \mu $ 和泊松率 $ \nu $ 。
 * 

 * 
 * 
 * @code
 *     struct Materials 
 *     { 
 *       double nu; 
 *       double mu; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Materials::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Material properties"); 
 *       { 
 *         prm.declare_entry("Poisson's ratio", 
 *                           "0.4999", 
 *                           Patterns::Double(-1.0, 0.5), 
 *                           "Poisson's ratio"); 
 * 
 *         prm.declare_entry("Shear modulus", 
 *                           "80.194e6", 
 *                           Patterns::Double(), 
 *                           "Shear modulus"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Materials::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Material properties"); 
 *       { 
 *         nu = prm.get_double("Poisson's ratio"); 
 *         mu = prm.get_double("Shear modulus"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Linearsolver"></a> 
 * <h4>Linear solver</h4>
 * 

 * 
 * 接下来，我们选择求解器和预处理器的设置。 当牛顿增量中出现大的非线性运动时，使用有效的前置条件对于确保收敛性至关重要。
 * 

 * 
 * 
 * @code
 *     struct LinearSolver 
 *     { 
 *       std::string type_lin; 
 *       double      tol_lin; 
 *       double      max_iterations_lin; 
 *       bool        use_static_condensation; 
 *       std::string preconditioner_type; 
 *       double      preconditioner_relaxation; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void LinearSolver::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Linear solver"); 
 *       { 
 *         prm.declare_entry("Solver type", 
 *                           "CG", 
 *                           Patterns::Selection("CG|Direct"), 
 *                           "Type of solver used to solve the linear system"); 
 * 
 *         prm.declare_entry("Residual", 
 *                           "1e-6", 
 *                           Patterns::Double(0.0), 
 *                           "Linear solver residual (scaled by residual norm)"); 
 * 
 *         prm.declare_entry( 
 *           "Max iteration multiplier", 
 *           "1", 
 *           Patterns::Double(0.0), 
 *           "Linear solver iterations (multiples of the system matrix size)"); 
 * 
 *         prm.declare_entry("Use static condensation", 
 *                           "true", 
 *                           Patterns::Bool(), 
 *                           "Solve the full block system or a reduced problem"); 
 * 
 *         prm.declare_entry("Preconditioner type", 
 *                           "ssor", 
 *                           Patterns::Selection("jacobi|ssor"), 
 *                           "Type of preconditioner"); 
 * 
 *         prm.declare_entry("Preconditioner relaxation", 
 *                           "0.65", 
 *                           Patterns::Double(0.0), 
 *                           "Preconditioner relaxation value"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void LinearSolver::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Linear solver"); 
 *       { 
 *         type_lin                  = prm.get("Solver type"); 
 *         tol_lin                   = prm.get_double("Residual"); 
 *         max_iterations_lin        = prm.get_double("Max iteration multiplier"); 
 *         use_static_condensation   = prm.get_bool("Use static condensation"); 
 *         preconditioner_type       = prm.get("Preconditioner type"); 
 *         preconditioner_relaxation = prm.get_double("Preconditioner relaxation"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Nonlinearsolver"></a> 
 * <h4>Nonlinear solver</h4>
 * 

 * 
 * 采用牛顿-拉弗森方案来解决非线性治理方程组。 我们现在定义牛顿-拉弗森非线性求解器的公差和最大迭代次数。
 * 

 * 
 * 
 * @code
 *     struct NonlinearSolver 
 *     { 
 *       unsigned int max_iterations_NR; 
 *       double       tol_f; 
 *       double       tol_u; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void NonlinearSolver::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Nonlinear solver"); 
 *       { 
 *         prm.declare_entry("Max iterations Newton-Raphson", 
 *                           "10", 
 *                           Patterns::Integer(0), 
 *                           "Number of Newton-Raphson iterations allowed"); 
 * 
 *         prm.declare_entry("Tolerance force", 
 *                           "1.0e-9", 
 *                           Patterns::Double(0.0), 
 *                           "Force residual tolerance"); 
 * 
 *         prm.declare_entry("Tolerance displacement", 
 *                           "1.0e-6", 
 *                           Patterns::Double(0.0), 
 *                           "Displacement error tolerance"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void NonlinearSolver::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Nonlinear solver"); 
 *       { 
 *         max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson"); 
 *         tol_f             = prm.get_double("Tolerance force"); 
 *         tol_u             = prm.get_double("Tolerance displacement"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Time"></a> 
 * <h4>Time</h4>
 * 

 * 
 * 设置时间步长 $ \varDelta t $ 和模拟结束时间。
 * 

 * 
 * 
 * @code
 *     struct Time 
 *     { 
 *       double delta_t; 
 *       double end_time; 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     void Time::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Time"); 
 *       { 
 *         prm.declare_entry("End time", "1", Patterns::Double(), "End time"); 
 * 
 *         prm.declare_entry("Time step size", 
 *                           "0.1", 
 *                           Patterns::Double(), 
 *                           "Time step size"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * 
 *     void Time::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       prm.enter_subsection("Time"); 
 *       { 
 *         end_time = prm.get_double("End time"); 
 *         delta_t  = prm.get_double("Time step size"); 
 *       } 
 *       prm.leave_subsection(); 
 *     } 
 * @endcode
 * 
 * 
 * <a name="Allparameters"></a> 
 * <h4>All parameters</h4>
 * 

 * 
 * 最后，我们将上述所有的结构合并到一个容器中，这个容器可以容纳我们所有的运行时选择。
 * 

 * 
 * 
 * @code
 *     struct AllParameters : public FESystem, 
 *                            public Geometry, 
 *                            public Materials, 
 *                            public LinearSolver, 
 *                            public NonlinearSolver, 
 *                            public Time 
 * 
 *     { 
 *       AllParameters(const std::string &input_file); 
 * 
 *       static void declare_parameters(ParameterHandler &prm); 
 * 
 *       void parse_parameters(ParameterHandler &prm); 
 *     }; 
 * 
 *     AllParameters::AllParameters(const std::string &input_file) 
 *     { 
 *       ParameterHandler prm; 
 *       declare_parameters(prm); 
 *       prm.parse_input(input_file); 
 *       parse_parameters(prm); 
 *     } 
 * 
 *     void AllParameters::declare_parameters(ParameterHandler &prm) 
 *     { 
 *       FESystem::declare_parameters(prm); 
 *       Geometry::declare_parameters(prm); 
 *       Materials::declare_parameters(prm); 
 *       LinearSolver::declare_parameters(prm); 
 *       NonlinearSolver::declare_parameters(prm); 
 *       Time::declare_parameters(prm); 
 *     } 
 * 
 *     void AllParameters::parse_parameters(ParameterHandler &prm) 
 *     { 
 *       FESystem::parse_parameters(prm); 
 *       Geometry::parse_parameters(prm); 
 *       Materials::parse_parameters(prm); 
 *       LinearSolver::parse_parameters(prm); 
 *       NonlinearSolver::parse_parameters(prm); 
 *       Time::parse_parameters(prm); 
 *     } 
 *   } // namespace Parameters 
 * @endcode
 * 
 * 
 * <a name="Timeclass"></a> 
 * <h3>Time class</h3>
 * 

 * 
 * 一个简单的类来存储时间数据。它的功能是透明的，所以没有必要讨论。为了简单起见，我们假设一个恒定的时间步长。
 * 

 * 
 * 
 * @code
 *   class Time 
 *   { 
 *   public: 
 *     Time(const double time_end, const double delta_t) 
 *       : timestep(0) 
 *       , time_current(0.0) 
 *       , time_end(time_end) 
 *       , delta_t(delta_t) 
 *     {} 
 * 
 *     virtual ~Time() = default; 
 * 
 *     double current() const 
 *     { 
 *       return time_current; 
 *     } 
 *     double end() const 
 *     { 
 *       return time_end; 
 *     } 
 *     double get_delta_t() const 
 *     { 
 *       return delta_t; 
 *     } 
 *     unsigned int get_timestep() const 
 *     { 
 *       return timestep; 
 *     } 
 *     void increment() 
 *     { 
 *       time_current += delta_t; 
 *       ++timestep; 
 *     } 
 * 
 *   private: 
 *     unsigned int timestep; 
 *     double       time_current; 
 *     const double time_end; 
 *     const double delta_t; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="CompressibleneoHookeanmaterialwithinathreefieldformulation"></a> 
 * <h3>Compressible neo-Hookean material within a three-field formulation</h3>
 * 

 * 
 * 正如介绍中所讨论的，新胡克材料是一种超弹性材料。 整个领域被假定为由可压缩的新胡克材料组成。 这个类别定义了这种材料在三场公式中的行为。 可压缩的新胡克材料可以用应变能量函数（SEF）来描述  $ \Psi = \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(\widetilde{J})$ 
 * 

 * 
 * 等效响应由 $ \Psi_{\text{iso}}(\overline{\mathbf{b}}) = c_{1} [\overline{I}_{1} - 3] $ 给出，其中 $ c_{1} = \frac{\mu}{2} $ 和 $\overline{I}_{1}$ 是左或右等效Cauchy-Green变形张量的第一不变量。这就是 $\overline{I}_1 \dealcoloneq \textrm{tr}(\overline{\mathbf{b}})$  。在这个例子中，支配体积响应的SEF被定义为 $ \Psi_{\text{vol}}(\widetilde{J}) = \kappa \frac{1}{4} [ \widetilde{J}^2 - 1 - 2\textrm{ln}\; \widetilde{J} ]$  ，其中 $\kappa \dealcoloneq \lambda + 2/3 \mu$  是<a href="http:en.wikipedia.org/wiki/Bulk_modulus">bulk modulus</a>， $\lambda$ 是<a href="http:en.wikipedia.org/wiki/Lam%C3%A9_parameters">Lam&eacute;'s first parameter</a>。
 * 

 * 
 * 下面的类将被用来描述我们工作中的材料特征，并提供了一个中心点，如果要实现不同的材料模型，就需要对其进行修改。为了使其发挥作用，我们将在每个正交点存储一个这种类型的对象，并在每个对象中存储当前状态（由三个场的值或度量来表征），这样我们就可以围绕当前状态计算出线性化的弹性系数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Material_Compressible_Neo_Hook_Three_Field 
 *   { 
 *   public: 
 *     Material_Compressible_Neo_Hook_Three_Field(const double mu, const double nu) 
 *       : kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu))) 
 *       , c_1(mu / 2.0) 
 *       , det_F(1.0) 
 *       , p_tilde(0.0) 
 *       , J_tilde(1.0) 
 *       , b_bar(Physics::Elasticity::StandardTensors<dim>::I) 
 *     { 
 *       Assert(kappa > 0, ExcInternalError()); 
 *     } 
 * 
 * @endcode
 * 
 * 我们用基于  $F$  和压力  $\widetilde{p}$  以及膨胀  $\widetilde{J}$  的各种变形相关数据来更新材料模型，并在函数的最后包括一个内部一致性的物理检查。
 * 

 * 
 * 
 * @code
 *     void update_material_data(const Tensor<2, dim> &F, 
 *                               const double          p_tilde_in, 
 *                               const double          J_tilde_in) 
 *     { 
 *       det_F                      = determinant(F); 
 *       const Tensor<2, dim> F_bar = Physics::Elasticity::Kinematics::F_iso(F); 
 *       b_bar                      = Physics::Elasticity::Kinematics::b(F_bar); 
 *       p_tilde                    = p_tilde_in; 
 *       J_tilde                    = J_tilde_in; 
 * 
 *       Assert(det_F > 0, ExcInternalError()); 
 *     } 
 * 
 * @endcode
 * 
 * 第二个函数决定了基尔霍夫应力  $\boldsymbol{\tau} = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$  。
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau() 
 *     { 
 *       return get_tau_iso() + get_tau_vol(); 
 *     } 
 * 
 * @endcode
 * 
 * 空间设置中的四阶弹性张量 $\mathfrak{c}$ 由SEF $\Psi$ 计算为 $ J \mathfrak{c}_{ijkl} = F_{iA} F_{jB} \mathfrak{C}_{ABCD} F_{kC} F_{lD}$  其中 $ \mathfrak{C} = 4 \frac{\partial^2 \Psi(\mathbf{C})}{\partial \mathbf{C} \partial \mathbf{C}}$  
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc() const 
 *     { 
 *       return get_Jc_vol() + get_Jc_iso(); 
 *     } 
 * 
 * @endcode
 * 
 * 体积自由能相对于  $\widetilde{J}$  的导数，返回  $\frac{\partial \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$  。
 * 
 * @code
 *     double get_dPsi_vol_dJ() const 
 *     { 
 *       return (kappa / 2.0) * (J_tilde - 1.0 / J_tilde); 
 *     } 
 * 
 * @endcode
 * 
 * 体积自由能的二次导数，相对于  $\widetilde{J}$  。我们需要在切线中明确地进行以下计算，所以我们将其公开。 我们计算出  $\frac{\partial^2 \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial \widetilde{J}}$  。
 * 
 * @code
 *     double get_d2Psi_vol_dJ2() const 
 *     { 
 *       return ((kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde))); 
 *     } 
 * 
 * @endcode
 * 
 * 接下来的几个函数会返回各种数据，我们选择将其与材料一起存储。
 * 

 * 
 * 
 * @code
 *     double get_det_F() const 
 *     { 
 *       return det_F; 
 *     } 
 * 
 *     double get_p_tilde() const 
 *     { 
 *       return p_tilde; 
 *     } 
 * 
 *     double get_J_tilde() const 
 *     { 
 *       return J_tilde; 
 *     } 
 * 
 *   protected: 
 * 
 * @endcode
 * 
 * 定义构成模型参数  $\kappa$  （体积模量）和新胡克模型参数  $c_1$  。
 * 

 * 
 * 
 * @code
 *     const double kappa; 
 *     const double c_1; 
 * 
 * @endcode
 * 
 * 模型的具体数据，方便与材料一起存储。
 * 

 * 
 * 
 * @code
 *     double                  det_F; 
 *     double                  p_tilde; 
 *     double                  J_tilde; 
 *     SymmetricTensor<2, dim> b_bar; 
 * 
 * @endcode
 * 
 * 以下函数在内部用于确定上述一些公共函数的结果。第一个函数决定了体积基尔霍夫应力  $\boldsymbol{\tau}_{\textrm{vol}}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_vol() const 
 *     { 
 *       return p_tilde * det_F * Physics::Elasticity::StandardTensors<dim>::I; 
 *     } 
 * 
 * @endcode
 * 
 * 接下来，确定等效基尔霍夫应力  $\boldsymbol{\tau}_{\textrm{iso}} = \mathcal{P}:\overline{\boldsymbol{\tau}}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_iso() const 
 *     { 
 *       return Physics::Elasticity::StandardTensors<dim>::dev_P * get_tau_bar(); 
 *     } 
 * 
 * @endcode
 * 
 * 然后，确定虚构的基尔霍夫应力  $\overline{\boldsymbol{\tau}}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<2, dim> get_tau_bar() const 
 *     { 
 *       return 2.0 * c_1 * b_bar; 
 *     } 
 * 
 * @endcode
 * 
 * 计算切线的体积部分  $J \mathfrak{c}_\textrm{vol}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc_vol() const 
 *     { 
 *       return p_tilde * det_F * 
 *              (Physics::Elasticity::StandardTensors<dim>::IxI - 
 *               (2.0 * Physics::Elasticity::StandardTensors<dim>::S)); 
 *     } 
 * 
 * @endcode
 * 
 * 计算切线的等值部分  $J \mathfrak{c}_\textrm{iso}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<4, dim> get_Jc_iso() const 
 *     { 
 *       const SymmetricTensor<2, dim> tau_bar = get_tau_bar(); 
 *       const SymmetricTensor<2, dim> tau_iso = get_tau_iso(); 
 *       const SymmetricTensor<4, dim> tau_iso_x_I = 
 *         outer_product(tau_iso, Physics::Elasticity::StandardTensors<dim>::I); 
 *       const SymmetricTensor<4, dim> I_x_tau_iso = 
 *         outer_product(Physics::Elasticity::StandardTensors<dim>::I, tau_iso); 
 *       const SymmetricTensor<4, dim> c_bar = get_c_bar(); 
 * 
 *       return (2.0 / dim) * trace(tau_bar) * 
 *                Physics::Elasticity::StandardTensors<dim>::dev_P - 
 *              (2.0 / dim) * (tau_iso_x_I + I_x_tau_iso) + 
 *              Physics::Elasticity::StandardTensors<dim>::dev_P * c_bar * 
 *                Physics::Elasticity::StandardTensors<dim>::dev_P; 
 *     } 
 * 
 * @endcode
 * 
 * 计算虚构的弹性张量  $\overline{\mathfrak{c}}$  。对于所选择的材料模型，这只是零。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<4, dim> get_c_bar() const 
 *     { 
 *       return SymmetricTensor<4, dim>(); 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Quadraturepointhistory"></a> 
 * <h3>Quadrature point history</h3>
 * 

 * 
 * 正如在 step-18 中所看到的， <code> PointHistory </code> 类提供了一种在正交点存储数据的方法。 这里每个正交点都持有一个指向材料描述的指针。 因此，不同的材料模型可以用在域的不同区域。 在其他数据中，我们选择为正交点存储Kirchhoff应力 $\boldsymbol{\tau}$ 和正切 $J\mathfrak{c}$ 。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class PointHistory 
 *   { 
 *   public: 
 *     PointHistory() 
 *       : F_inv(Physics::Elasticity::StandardTensors<dim>::I) 
 *       , tau(SymmetricTensor<2, dim>()) 
 *       , d2Psi_vol_dJ2(0.0) 
 *       , dPsi_vol_dJ(0.0) 
 *       , Jc(SymmetricTensor<4, dim>()) 
 *     {} 
 * 
 *     virtual ~PointHistory() = default; 
 * 
 * @endcode
 * 
 * 第一个函数用于创建一个材料对象并正确初始化所有的张量。第二个函数根据当前的变形量 $\textrm{Grad}\mathbf{u}_{\textrm{n}}$ 、压力 $\widetilde{p}$ 和扩张 $\widetilde{J}$ 场值更新存储的数值和应力。
 * 

 * 
 * 
 * @code
 *     void setup_lqp(const Parameters::AllParameters &parameters) 
 *     { 
 *       material = 
 *         std::make_shared<Material_Compressible_Neo_Hook_Three_Field<dim>>( 
 *           parameters.mu, parameters.nu); 
 *       update_values(Tensor<2, dim>(), 0.0, 1.0); 
 *     } 
 * 
 * @endcode
 * 
 * 为此，我们从位移梯度 $\textrm{Grad}\ \mathbf{u}$ 中计算出变形梯度 $\mathbf{F}$  ，即 $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad}\ \mathbf{u}$ ，然后让与这个正交点相关的材料模型进行自我更新。当计算变形梯度时，我们必须注意与哪些数据类型进行比较 $\mathbf{I} + \textrm{Grad}\ \mathbf{u}$ ：由于 $I$ 有数据类型SymmetricTensor，只要写 <code>I + Grad_u_n</code> 就可以将第二个参数转换为对称张量，进行求和，然后将结果投给Tensor（即可能是非对称张量的类型）。然而，由于 <code>Grad_u_n</code> 在一般情况下是非对称的，转换为SymmetricTensor将会失败。我们可以通过先将 $I$ 转换为Tensor，然后像在非对称张量之间一样执行加法，来避免这种来回折腾。
 * 

 * 
 * 
 * @code
 *     void update_values(const Tensor<2, dim> &Grad_u_n, 
 *                        const double          p_tilde, 
 *                        const double          J_tilde) 
 *     { 
 *       const Tensor<2, dim> F = Physics::Elasticity::Kinematics::F(Grad_u_n); 
 *       material->update_material_data(F, p_tilde, J_tilde); 
 * 
 * @endcode
 * 
 * 材料已经更新，所以我们现在计算基尔霍夫应力 $\mathbf{\tau}$ ，切线 $J\mathfrak{c}$ 和体积自由能的一、二次导数。
 * 

 * 
 * 我们还存储了变形梯度的逆值，因为我们经常使用它。
 * 

 * 
 * 
 * @code
 *       F_inv         = invert(F); 
 *       tau           = material->get_tau(); 
 *       Jc            = material->get_Jc(); 
 *       dPsi_vol_dJ   = material->get_dPsi_vol_dJ(); 
 *       d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2(); 
 *     } 
 * 
 * @endcode
 * 
 * 我们提供一个接口来检索某些数据。 下面是运动学变量。
 * 

 * 
 * 
 * @code
 *     double get_J_tilde() const 
 *     { 
 *       return material->get_J_tilde(); 
 *     } 
 * 
 *     double get_det_F() const 
 *     { 
 *       return material->get_det_F(); 
 *     } 
 * 
 *     const Tensor<2, dim> &get_F_inv() const 
 *     { 
 *       return F_inv; 
 *     } 
 * 
 * @endcode
 * 
 * ...和动能变量。 这些在材料和全局切线矩阵以及残余装配操作中使用。
 * 

 * 
 * 
 * @code
 *     double get_p_tilde() const 
 *     { 
 *       return material->get_p_tilde(); 
 *     } 
 * 
 *     const SymmetricTensor<2, dim> &get_tau() const 
 *     { 
 *       return tau; 
 *     } 
 * 
 *     double get_dPsi_vol_dJ() const 
 *     { 
 *       return dPsi_vol_dJ; 
 *     } 
 * 
 *     double get_d2Psi_vol_dJ2() const 
 *     { 
 *       return d2Psi_vol_dJ2; 
 *     } 
 * 
 * @endcode
 * 
 * 最后是切线。
 * 

 * 
 * 
 * @code
 *     const SymmetricTensor<4, dim> &get_Jc() const 
 *     { 
 *       return Jc; 
 *     } 
 * 
 * @endcode
 * 
 * 在成员函数方面，这个类为它所代表的正交点存储了一个材料类型的副本，以备在域的不同区域使用不同的材料，以及变形梯度的逆值...
 * 

 * 
 * 
 * @code
 *   private: 
 *     std::shared_ptr<Material_Compressible_Neo_Hook_Three_Field<dim>> material; 
 * 
 *     Tensor<2, dim> F_inv; 
 * 
 * @endcode
 * 
 * ...... 和应力型变量以及切线  $J\mathfrak{c}$  。
 * 

 * 
 * 
 * @code
 *     SymmetricTensor<2, dim> tau; 
 *     double                  d2Psi_vol_dJ2; 
 *     double                  dPsi_vol_dJ; 
 * 
 *     SymmetricTensor<4, dim> Jc; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Quasistaticquasiincompressiblefinitestrainsolid"></a> 
 * <h3>Quasi-static quasi-incompressible finite-strain solid</h3>
 * 

 * 
 * Solid类是中心类，它代表了手头的问题。它遵循通常的方案，即它真正拥有的是一个构造函数、解构函数和一个 <code>run()</code> 函数，该函数将所有的工作分派给这个类的私有函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solid 
 *   { 
 *   public: 
 *     Solid(const std::string &input_file); 
 * 
 *     void run(); 
 * 
 *   private: 
 * 
 * @endcode
 * 
 * 在这个类的私有部分，我们首先向前声明一些对象，这些对象在使用WorkStream对象进行并行工作时使用（关于这方面的更多信息，请参见 @ref threads 模块）。
 * 

 * 
 * 我们声明这样的结构，用于正切（刚度）矩阵和右手边矢量的计算，静态冷凝，以及更新正交点。
 * 

 * 
 * 
 * @code
 *     struct PerTaskData_ASM; 
 *     struct ScratchData_ASM; 
 * 
 *     struct PerTaskData_SC; 
 *     struct ScratchData_SC; 
 * 
 *     struct PerTaskData_UQPH; 
 *     struct ScratchData_UQPH; 
 * 
 * @endcode
 * 
 * 我们从一个建立网格的成员函数开始收集。
 * 

 * 
 * 
 * @code
 *     void make_grid(); 
 * 
 * @endcode
 * 
 * 设置要解决的有限元系统。
 * 

 * 
 * 
 * @code
 *     void system_setup(); 
 * 
 *     void determine_component_extractors(); 
 * 
 * @endcode
 * 
 * 为增量位移场创建Dirichlet约束。
 * 

 * 
 * 
 * @code
 *     void make_constraints(const int it_nr); 
 * 
 * @endcode
 * 
 * 使用多线程的几个函数来组装系统和右手边的矩阵。它们中的每一个都是包装函数，一个是在WorkStream模型中对一个单元进行工作的执行函数，另一个是将对这一个单元的工作复制到代表它的全局对象中。
 * 

 * 
 * 
 * @code
 *     void assemble_system(); 
 * 
 *     void assemble_system_one_cell( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       ScratchData_ASM &                                     scratch, 
 *       PerTaskData_ASM &                                     data) const; 
 * 
 * @endcode
 * 
 * 还有类似的，执行全局静态冷凝。
 * 

 * 
 * 
 * @code
 *     void assemble_sc(); 
 * 
 *     void assemble_sc_one_cell( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       ScratchData_SC &                                      scratch, 
 *       PerTaskData_SC &                                      data); 
 * 
 *     void copy_local_to_global_sc(const PerTaskData_SC &data); 
 * 
 * @endcode
 * 
 * 创建并更新正交点。在这里，没有数据需要被复制到全局对象中，所以copy_local_to_global函数是空的。
 * 

 * 
 * 
 * @code
 *     void setup_qph(); 
 * 
 *     void update_qph_incremental(const BlockVector<double> &solution_delta); 
 * 
 *     void update_qph_incremental_one_cell( 
 *       const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *       ScratchData_UQPH &                                    scratch, 
 *       PerTaskData_UQPH &                                    data); 
 * 
 *     void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/) 
 *     {} 
 * 
 * @endcode
 * 
 * 用牛顿-拉弗森方法求解位移。我们把这个函数分成非线性循环和解决线性化的Newton-Raphson步骤的函数。
 * 

 * 
 * 
 * @code
 *     void solve_nonlinear_timestep(BlockVector<double> &solution_delta); 
 * 
 *     std::pair<unsigned int, double> 
 *     solve_linear_system(BlockVector<double> &newton_update); 
 * 
 * @endcode
 * 
 * 检索解决方案，以及后期处理和将数据写入文件。
 * 

 * 
 * 
 * @code
 *     BlockVector<double> 
 *     get_total_solution(const BlockVector<double> &solution_delta) const; 
 * 
 *     void output_results() const; 
 * 
 * @endcode
 * 
 * 最后是一些描述当前状态的成员变量。一个用于描述问题设置的参数集合...
 * 

 * 
 * 
 * @code
 *     Parameters::AllParameters parameters; 
 * 
 * @endcode
 * 
 * ...参考配置的体积...
 * 

 * 
 * 
 * @code
 *     double vol_reference; 
 * 
 * @endcode
 * 
 * ......以及对解决问题的几何形状的描述。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 * 
 * @endcode
 * 
 * 同时，记录当前时间和评估某些函数的时间
 * 

 * 
 * 
 * @code
 *     Time                time; 
 *     mutable TimerOutput timer; 
 * 
 * @endcode
 * 
 * 一个存储正交点信息的对象。与 step-18 不同，这里采用了deal.II的本地正交点数据管理器。
 * 

 * 
 * 
 * @code
 *     CellDataStorage<typename Triangulation<dim>::cell_iterator, 
 *                     PointHistory<dim>> 
 *       quadrature_point_history; 
 * 
 * @endcode
 * 
 * 对有限元系统的描述，包括位移多项式程度、自由度处理程序、每个单元的DoF数量以及用于从解向量中检索信息的提取器对象。
 * 

 * 
 * 
 * @code
 *     const unsigned int               degree; 
 *     const FESystem<dim>              fe; 
 *     DoFHandler<dim>                  dof_handler; 
 *     const unsigned int               dofs_per_cell; 
 *     const FEValuesExtractors::Vector u_fe; 
 *     const FEValuesExtractors::Scalar p_fe; 
 *     const FEValuesExtractors::Scalar J_fe; 
 * 
 * @endcode
 * 
 * 说明块系统是如何安排的。有3个块，第一个包含一个矢量DOF  $\mathbf{u}$  ，而另外两个描述标量DOF， $\widetilde{p}$  和  $\widetilde{J}$  。
 * 

 * 
 * 
 * @code
 *     static const unsigned int n_blocks          = 3; 
 *     static const unsigned int n_components      = dim + 2; 
 *     static const unsigned int first_u_component = 0; 
 *     static const unsigned int p_component       = dim; 
 *     static const unsigned int J_component       = dim + 1; 
 * 
 *     enum 
 *     { 
 *       u_dof = 0, 
 *       p_dof = 1, 
 *       J_dof = 2 
 *     }; 
 * 
 *     std::vector<types::global_dof_index> dofs_per_block; 
 *     std::vector<types::global_dof_index> element_indices_u; 
 *     std::vector<types::global_dof_index> element_indices_p; 
 *     std::vector<types::global_dof_index> element_indices_J; 
 * 
 * @endcode
 * 
 * 单元和面的高斯正交规则。单元和面的正交点的数量被记录下来。
 * 

 * 
 * 
 * @code
 *     const QGauss<dim>     qf_cell; 
 *     const QGauss<dim - 1> qf_face; 
 *     const unsigned int    n_q_points; 
 *     const unsigned int    n_q_points_f; 
 * 
 * @endcode
 * 
 * 用于存储收敛的解和右手边向量以及切线矩阵的对象。有一个AffineConstraints对象，用于跟踪约束条件。 我们利用了为块状系统设计的稀疏性模式。
 * 

 * 
 * 
 * @code
 *     AffineConstraints<double> constraints; 
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> tangent_matrix; 
 *     BlockVector<double>       system_rhs; 
 *     BlockVector<double>       solution_n; 
 * @endcode
 * 
 * 然后
 * 定义一些变量来存储规范，并更新规范和归一化系数。
 * 

 * 
 * 
 * @code
 *     struct Errors 
 *     { 
 *       Errors() 
 *         : norm(1.0) 
 *         , u(1.0) 
 *         , p(1.0) 
 *         , J(1.0) 
 *       {} 
 * 
 *       void reset() 
 *       { 
 *         norm = 1.0; 
 *         u    = 1.0; 
 *         p    = 1.0; 
 *         J    = 1.0; 
 *       } 
 *       void normalize(const Errors &rhs) 
 *       { 
 *         if (rhs.norm != 0.0) 
 *           norm /= rhs.norm; 
 *         if (rhs.u != 0.0) 
 *           u /= rhs.u; 
 *         if (rhs.p != 0.0) 
 *           p /= rhs.p; 
 *         if (rhs.J != 0.0) 
 *           J /= rhs.J; 
 *       } 
 * 
 *       double norm, u, p, J; 
 *     }; 
 * 
 *     Errors error_residual, error_residual_0, error_residual_norm, error_update, 
 *       error_update_0, error_update_norm; 
 * 
 * @endcode
 * 
 * 计算误差措施的方法
 * 

 * 
 * 
 * @code
 *     void get_error_residual(Errors &error_residual); 
 * 
 *     void get_error_update(const BlockVector<double> &newton_update, 
 *                           Errors &                   error_update); 
 * 
 *     std::pair<double, double> get_error_dilation() const; 
 * 
 * @endcode
 * 
 * 计算空间配置中的体积
 * 

 * 
 * 
 * @code
 *     double compute_vol_current() const; 
 * 
 * @endcode
 * 
 * 以悦目的方式向屏幕打印信息...
 * 

 * 
 * 
 * @code
 *     static void print_conv_header(); 
 * 
 *     void print_conv_footer(); 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeSolidcodeclass"></a> 
 * <h3>Implementation of the <code>Solid</code> class</h3>
 * 
 * <a name="Publicinterface"></a> 
 * <h4>Public interface</h4>
 * 

 * 
 * 我们使用从参数文件中提取的数据来初始化Solid类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   Solid<dim>::Solid(const std::string &input_file) 
 *     : parameters(input_file) 
 *     , vol_reference(0.) 
 *     , triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , time(parameters.end_time, parameters.delta_t) 
 *     , timer(std::cout, TimerOutput::summary, TimerOutput::wall_times) 
 *     , degree(parameters.poly_degree) 
 *     , 
 * 
 * @endcode
 * 
 * 有限元系统是由昏暗的连续位移DOF和不连续的压力和膨胀DOF组成。为了满足Babuska-Brezzi或LBB稳定性条件（见Hughes（2000）），我们设置了一个 $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$ 系统。  $Q_2 \times DGPM_1 \times DGPM_1$ 元素满足这个条件，而 $Q_1 \times DGPM_0 \times DGPM_0$ 元素不满足。然而，事实证明，后者还是表现出良好的收敛特性。
 * 

 * 
 * 
 * @code
 *     fe(FE_Q<dim>(parameters.poly_degree), 
 *        dim, // displacement 
 *        FE_DGPMonomial<dim>(parameters.poly_degree - 1), 
 *        1, // pressure 
 *        FE_DGPMonomial<dim>(parameters.poly_degree - 1), 
 *        1) 
 *     , // dilatation 
 *     dof_handler(triangulation) 
 *     , dofs_per_cell(fe.n_dofs_per_cell()) 
 *     , u_fe(first_u_component) 
 *     , p_fe(p_component) 
 *     , J_fe(J_component) 
 *     , dofs_per_block(n_blocks) 
 *     , qf_cell(parameters.quad_order) 
 *     , qf_face(parameters.quad_order) 
 *     , n_q_points(qf_cell.size()) 
 *     , n_q_points_f(qf_face.size()) 
 *   { 
 *     Assert(dim == 2 || dim == 3, 
 *            ExcMessage("This problem only works in 2 or 3 space dimensions.")); 
 *     determine_component_extractors(); 
 *   } 
 * 
 * @endcode
 * 
 * 在解决准静态问题时，时间成为一个加载参数，即我们随着时间线性增加加载量，使得这两个概念可以互换。我们选择用恒定的时间步长来线性递增时间。
 * 

 * 
 * 我们从预处理开始，设置初始扩张值，然后输出初始网格，然后开始模拟，开始第一次时间（和载荷）递增。
 * 

 * 
 * 在对初始解场施加约束 $\widetilde{J}=1$ 时，必须注意（或者至少要考虑一下）。该约束对应于未变形构型中变形梯度的行列式，也就是身份张量。我们使用FE_DGPMonomial基数来插值扩张场，因此我们不能简单地将相应的dof设置为unity，因为它们对应于单项式系数。因此，我们使用 VectorTools::project 函数来为我们做这项工作。 VectorTools::project 函数需要一个参数，表明悬挂节点的约束。我们在这个程序中没有 所以我们必须创建一个约束对象。在原始状态下，约束对象是没有排序的，必须先进行排序（使用 AffineConstraints::close 函数）才能使用。请看  step-21  以了解更多信息。我们只需要强制执行扩张的初始条件。为了做到这一点，我们使用ComponentSelectFunction，它作为一个掩码，将n_components的J_component设置为1。 这正是我们想要的。请看 step-20 中的用法，了解更多信息。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::run() 
 *   { 
 *     make_grid(); 
 *     system_setup(); 
 *     { 
 *       AffineConstraints<double> constraints; 
 *       constraints.close(); 
 * 
 *       const ComponentSelectFunction<dim> J_mask(J_component, n_components); 
 * 
 *       VectorTools::project( 
 *         dof_handler, constraints, QGauss<dim>(degree + 2), J_mask, solution_n); 
 *     } 
 *     output_results(); 
 *     time.increment(); 
 * 
 * @endcode
 * 
 * 然后我们宣布增量解决方案更新 $\varDelta \mathbf{\Xi} \dealcoloneq \{\varDelta \mathbf{u},\varDelta \widetilde{p}, \varDelta \widetilde{J} \}$ 并开始在时域上循环。
 * 

 * 
 * 在开始的时候，我们重置这个时间步长的解决方案更新...
 * 

 * 
 * 
 * @code
 *     BlockVector<double> solution_delta(dofs_per_block); 
 *     while (time.current() < time.end()) 
 *       { 
 *         solution_delta = 0.0; 
 * 
 * @endcode
 * 
 * ...求解当前时间步长并更新总解向量  $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} + \varDelta \mathbf{\Xi}$  ...
 * 

 * 
 * 
 * @code
 *         solve_nonlinear_timestep(solution_delta); 
 *         solution_n += solution_delta; 
 * 
 * @endcode
 * 
 * ...并在快乐地进入下一个时间步骤之前绘制结果。
 * 

 * 
 * 
 * @code
 *         output_results(); 
 *         time.increment(); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Privateinterface"></a> 
 * <h3>Private interface</h3>
 * 
 * <a name="Threadingbuildingblocksstructures"></a> 
 * <h4>Threading-building-blocks structures</h4>
 * 

 * 
 * 第一组私有成员函数与并行化有关。我们使用线程积木库（TBB）来执行尽可能多的计算密集型分布式任务。特别是，我们使用TBB组装正切矩阵和右手向量、静态凝结贡献，以及更新存储在正交点的数据。我们在这方面的主要工具是WorkStream类（更多信息见 @ref 线程模块）。
 * 

 * 
 * 首先我们要处理正切矩阵和右手边的装配结构。PerTaskData对象存储了本地对全局系统的贡献。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::PerTaskData_ASM 
 *   { 
 *     FullMatrix<double>                   cell_matrix; 
 *     Vector<double>                       cell_rhs; 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 * 
 *     PerTaskData_ASM(const unsigned int dofs_per_cell) 
 *       : cell_matrix(dofs_per_cell, dofs_per_cell) 
 *       , cell_rhs(dofs_per_cell) 
 *       , local_dof_indices(dofs_per_cell) 
 *     {} 
 * 
 *     void reset() 
 *     { 
 *       cell_matrix = 0.0; 
 *       cell_rhs    = 0.0; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 另一方面，ScratchData对象存储了较大的对象，如形状函数值数组（  <code>Nx</code>  ）和形状函数梯度和对称梯度向量，我们将在装配时使用。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::ScratchData_ASM 
 *   { 
 *     FEValues<dim>     fe_values; 
 *     FEFaceValues<dim> fe_face_values; 
 * 
 *     std::vector<std::vector<double>>                  Nx; 
 *     std::vector<std::vector<Tensor<2, dim>>>          grad_Nx; 
 *     std::vector<std::vector<SymmetricTensor<2, dim>>> symm_grad_Nx; 
 * 
 *     ScratchData_ASM(const FiniteElement<dim> &fe_cell, 
 *                     const QGauss<dim> &       qf_cell, 
 *                     const UpdateFlags         uf_cell, 
 *                     const QGauss<dim - 1> &   qf_face, 
 *                     const UpdateFlags         uf_face) 
 *       : fe_values(fe_cell, qf_cell, uf_cell) 
 *       , fe_face_values(fe_cell, qf_face, uf_face) 
 *       , Nx(qf_cell.size(), std::vector<double>(fe_cell.n_dofs_per_cell())) 
 *       , grad_Nx(qf_cell.size(), 
 *                 std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell())) 
 *       , symm_grad_Nx(qf_cell.size(), 
 *                      std::vector<SymmetricTensor<2, dim>>( 
 *                        fe_cell.n_dofs_per_cell())) 
 *     {} 
 * 
 *     ScratchData_ASM(const ScratchData_ASM &rhs) 
 *       : fe_values(rhs.fe_values.get_fe(), 
 *                   rhs.fe_values.get_quadrature(), 
 *                   rhs.fe_values.get_update_flags()) 
 *       , fe_face_values(rhs.fe_face_values.get_fe(), 
 *                        rhs.fe_face_values.get_quadrature(), 
 *                        rhs.fe_face_values.get_update_flags()) 
 *       , Nx(rhs.Nx) 
 *       , grad_Nx(rhs.grad_Nx) 
 *       , symm_grad_Nx(rhs.symm_grad_Nx) 
 *     {} 
 * 
 *     void reset() 
 *     { 
 *       const unsigned int n_q_points      = Nx.size(); 
 *       const unsigned int n_dofs_per_cell = Nx[0].size(); 
 *       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *         { 
 *           Assert(Nx[q_point].size() == n_dofs_per_cell, ExcInternalError()); 
 *           Assert(grad_Nx[q_point].size() == n_dofs_per_cell, 
 *                  ExcInternalError()); 
 *           Assert(symm_grad_Nx[q_point].size() == n_dofs_per_cell, 
 *                  ExcInternalError()); 
 *           for (unsigned int k = 0; k < n_dofs_per_cell; ++k) 
 *             { 
 *               Nx[q_point][k]           = 0.0; 
 *               grad_Nx[q_point][k]      = 0.0; 
 *               symm_grad_Nx[q_point][k] = 0.0; 
 *             } 
 *         } 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 然后我们定义结构来组装静态凝结的切线矩阵。回顾一下，我们希望解决一个基于位移的公式。由于 $\widetilde{p}$ 和 $\widetilde{J}$ 字段在元素层面上是不连续的，所以我们在元素层面上进行缩合。 由于这些操作是基于矩阵的，我们需要设置一些矩阵来存储一些切线矩阵子块的局部贡献。 我们把这些放在PerTaskData结构中。
 * 

 * 
 * 我们选择不在 <code>reset()</code> 函数中重置任何数据，因为矩阵提取和替换工具会处理这个问题。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::PerTaskData_SC 
 *   { 
 *     FullMatrix<double>                   cell_matrix; 
 *     std::vector<types::global_dof_index> local_dof_indices; 
 * 
 *     FullMatrix<double> k_orig; 
 *     FullMatrix<double> k_pu; 
 *     FullMatrix<double> k_pJ; 
 *     FullMatrix<double> k_JJ; 
 *     FullMatrix<double> k_pJ_inv; 
 *     FullMatrix<double> k_bbar; 
 *     FullMatrix<double> A; 
 *     FullMatrix<double> B; 
 *     FullMatrix<double> C; 
 * 
 *     PerTaskData_SC(const unsigned int dofs_per_cell, 
 *                    const unsigned int n_u, 
 *                    const unsigned int n_p, 
 *                    const unsigned int n_J) 
 *       : cell_matrix(dofs_per_cell, dofs_per_cell) 
 *       , local_dof_indices(dofs_per_cell) 
 *       , k_orig(dofs_per_cell, dofs_per_cell) 
 *       , k_pu(n_p, n_u) 
 *       , k_pJ(n_p, n_J) 
 *       , k_JJ(n_J, n_J) 
 *       , k_pJ_inv(n_p, n_J) 
 *       , k_bbar(n_u, n_u) 
 *       , A(n_J, n_u) 
 *       , B(n_J, n_u) 
 *       , C(n_p, n_u) 
 *     {} 
 * 
 *     void reset() 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * 我们希望在这里执行的操作的ScratchData对象是空的，因为我们不需要临时数据，但它仍然需要为当前deal.II中TBB的实现而定义。 所以我们为此创建了一个假的结构。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::ScratchData_SC 
 *   { 
 *     void reset() 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * 最后我们定义结构以协助更新正交点信息。与SC的装配过程类似，我们不需要PerTaskData对象（因为这里没有什么可存储的），但还是必须定义一个。请注意，这是因为对于我们这里的操作--更新正交点的数据--是纯粹的局部操作：我们在每个单元上做的事情在每个单元上都会被消耗掉，没有像使用WorkStream类时通常会有的全局聚合操作。我们仍然必须定义每个任务的数据结构，这表明WorkStream类可能不适合这种操作（原则上，我们可以简单地为每个单元使用 Threads::new_task 创建一个新的任务），但无论如何这样做也没有什么坏处。此外，如果一个正交点有不同的材料模型，需要不同程度的计算费用，那么这里使用的方法可能是有利的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::PerTaskData_UQPH 
 *   { 
 *     void reset() 
 *     {} 
 *   }; 
 * 
 * @endcode
 * 
 * ScratchData对象将被用来存储解向量的别名，这样我们就不必复制这个大的数据结构。然后我们定义一些向量来提取正交点的解值和梯度。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct Solid<dim>::ScratchData_UQPH 
 *   { 
 *     const BlockVector<double> &solution_total; 
 * 
 *     std::vector<Tensor<2, dim>> solution_grads_u_total; 
 *     std::vector<double>         solution_values_p_total; 
 *     std::vector<double>         solution_values_J_total; 
 * 
 *     FEValues<dim> fe_values; 
 * 
 *     ScratchData_UQPH(const FiniteElement<dim> & fe_cell, 
 *                      const QGauss<dim> &        qf_cell, 
 *                      const UpdateFlags          uf_cell, 
 *                      const BlockVector<double> &solution_total) 
 *       : solution_total(solution_total) 
 *       , solution_grads_u_total(qf_cell.size()) 
 *       , solution_values_p_total(qf_cell.size()) 
 *       , solution_values_J_total(qf_cell.size()) 
 *       , fe_values(fe_cell, qf_cell, uf_cell) 
 *     {} 
 * 
 *     ScratchData_UQPH(const ScratchData_UQPH &rhs) 
 *       : solution_total(rhs.solution_total) 
 *       , solution_grads_u_total(rhs.solution_grads_u_total) 
 *       , solution_values_p_total(rhs.solution_values_p_total) 
 *       , solution_values_J_total(rhs.solution_values_J_total) 
 *       , fe_values(rhs.fe_values.get_fe(), 
 *                   rhs.fe_values.get_quadrature(), 
 *                   rhs.fe_values.get_update_flags()) 
 *     {} 
 * 
 *     void reset() 
 *     { 
 *       const unsigned int n_q_points = solution_grads_u_total.size(); 
 *       for (unsigned int q = 0; q < n_q_points; ++q) 
 *         { 
 *           solution_grads_u_total[q]  = 0.0; 
 *           solution_values_p_total[q] = 0.0; 
 *           solution_values_J_total[q] = 0.0; 
 *         } 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Solidmake_grid"></a> 
 * <h4>Solid::make_grid</h4>
 * 

 * 
 * 进入第一个私有成员函数。在这里我们创建域的三角形，为此我们选择了按比例的立方体，每个面都有一个边界ID号。 对于缩进问题，网格必须至少被细化一次。
 * 

 * 
 * 然后，我们确定参考配置的体积，并将其打印出来进行比较。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::make_grid() 
 *   { 
 *     GridGenerator::hyper_rectangle( 
 *       triangulation, 
 *       (dim == 3 ? Point<dim>(0.0, 0.0, 0.0) : Point<dim>(0.0, 0.0)), 
 *       (dim == 3 ? Point<dim>(1.0, 1.0, 1.0) : Point<dim>(1.0, 1.0)), 
 *       true); 
 *     GridTools::scale(parameters.scale, triangulation); 
 *     triangulation.refine_global(std::max(1U, parameters.global_refinement)); 
 * 
 *     vol_reference = GridTools::volume(triangulation); 
 *     std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl; 
 * 
 * @endcode
 * 
 * 由于我们希望对顶面的一个补丁应用诺伊曼BC，我们必须找到域的这一部分的单元格面，并用一个明显的边界ID号来标记它们。 我们要找的面在+y面上，将得到边界ID 6（0到5已经在创建立方体域的六个面时使用了）。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         { 
 *           if (face->at_boundary() == true && 
 *               face->center()[1] == 1.0 * parameters.scale) 
 *             { 
 *               if (dim == 3) 
 *                 { 
 *                   if (face->center()[0] < 0.5 * parameters.scale && 
 *                       face->center()[2] < 0.5 * parameters.scale) 
 *                     face->set_boundary_id(6); 
 *                 } 
 *               else 
 *                 { 
 *                   if (face->center()[0] < 0.5 * parameters.scale) 
 *                     face->set_boundary_id(6); 
 *                 } 
 *             } 
 *         } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidsystem_setup"></a> 
 * <h4>Solid::system_setup</h4>
 * 

 * 
 * 接下来我们描述FE系统是如何设置的。 我们首先确定每块的分量数量。由于位移是一个矢量分量，所以前两个分量属于它，而后两个分量描述标量压力和扩张DOF。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::system_setup() 
 *   { 
 *     timer.enter_subsection("Setup system"); 
 * 
 *     std::vector<unsigned int> block_component(n_components, 
 *                                               u_dof); // Displacement 
 *     block_component[p_component] = p_dof;             // Pressure 
 *     block_component[J_component] = J_dof;             // Dilatation 
 * 
 * @endcode
 * 
 * 然后，DOF处理程序被初始化，我们以一种有效的方式对网格进行重新编号。我们还记录了每块DOF的数量。
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 *     DoFRenumbering::Cuthill_McKee(dof_handler); 
 *     DoFRenumbering::component_wise(dof_handler, block_component); 
 * 
 *     dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 * 
 *     std::cout << "Triangulation:" 
 *               << "\n\t Number of active cells: " 
 *               << triangulation.n_active_cells() 
 *               << "\n\t Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 * @endcode
 * 
 * 设置稀疏模式和切线矩阵
 * 

 * 
 * 
 * @code
 *     tangent_matrix.clear(); 
 *     { 
 *       const types::global_dof_index n_dofs_u = dofs_per_block[u_dof]; 
 *       const types::global_dof_index n_dofs_p = dofs_per_block[p_dof]; 
 *       const types::global_dof_index n_dofs_J = dofs_per_block[J_dof]; 
 * 
 *       BlockDynamicSparsityPattern dsp(n_blocks, n_blocks); 
 * 
 *       dsp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u); 
 *       dsp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p); 
 *       dsp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J); 
 * 
 *       dsp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u); 
 *       dsp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p); 
 *       dsp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J); 
 * 
 *       dsp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u); 
 *       dsp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p); 
 *       dsp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J); 
 *       dsp.collect_sizes(); 
 * 
 * @endcode
 * 
 * 全局系统矩阵最初具有以下结构 
 * @f{align*}
 * \underbrace{\begin{bmatrix}
 * \mathsf{\mathbf{K}}_{uu}  & \mathsf{\mathbf{K}}_{u\widetilde{p}} &
 * \mathbf{0}
 * \\ \mathsf{\mathbf{K}}_{\widetilde{p}u} & \mathbf{0} &
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}
 * \\ \mathbf{0} & \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}} &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * \end{bmatrix}}_{\mathsf{\mathbf{K}}(\mathbf{\Xi}_{\textrm{i}})}
 * \underbrace{\begin{bmatrix}
 * d \mathsf{u}
 * \\  d \widetilde{\mathsf{\mathbf{p}}}
 * \\  d \widetilde{\mathsf{\mathbf{J}}}
 * \end{bmatrix}}_{d \mathbf{\Xi}}
 * =
 * \underbrace{\begin{bmatrix}
 * \mathsf{\mathbf{F}}_{u}(\mathbf{u}_{\textrm{i}})
 * \\ \mathsf{\mathbf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}})
 * \\ \mathsf{\mathbf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
 * \end{bmatrix}}_{ \mathsf{\mathbf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .
 * @f}
 * 我们优化稀疏模式以反映这一结构，并防止为右对角块成分创建不必要的数据。
 * 

 * 
 * 
 * @code
 *       Table<2, DoFTools::Coupling> coupling(n_components, n_components); 
 *       for (unsigned int ii = 0; ii < n_components; ++ii) 
 *         for (unsigned int jj = 0; jj < n_components; ++jj) 
 *           if (((ii < p_component) && (jj == J_component)) || 
 *               ((ii == J_component) && (jj < p_component)) || 
 *               ((ii == p_component) && (jj == p_component))) 
 *             coupling[ii][jj] = DoFTools::none; 
 *           else 
 *             coupling[ii][jj] = DoFTools::always; 
 *       DoFTools::make_sparsity_pattern( 
 *         dof_handler, coupling, dsp, constraints, false); 
 *       sparsity_pattern.copy_from(dsp); 
 *     } 
 * 
 *     tangent_matrix.reinit(sparsity_pattern); 
 * 
 * @endcode
 * 
 * 然后，我们设置了存储向量
 * 

 * 
 * 
 * @code
 *     system_rhs.reinit(dofs_per_block); 
 *     system_rhs.collect_sizes(); 
 * 
 *     solution_n.reinit(dofs_per_block); 
 *     solution_n.collect_sizes(); 
 * 
 * @endcode
 * 
 * ...最后设置正交点历史。
 * 

 * 
 * 
 * @code
 *     setup_qph(); 
 * 
 *     timer.leave_subsection(); 
 *   } 
 * @endcode
 * 
 * 接下来我们从FE系统中计算出一些信息，描述哪些局部元素DOF连接到哪个块组件上。 这将在后面用于从全局矩阵中提取子块。
 * 

 * 
 * 本质上，我们所需要的就是让FES系统对象指出参考单元上的DOF连接到哪个块状部件上。 目前，插值域的设置是这样的：0表示位移DOF，1表示压力DOF，2表示膨胀DOF。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::determine_component_extractors() 
 *   { 
 *     element_indices_u.clear(); 
 *     element_indices_p.clear(); 
 *     element_indices_J.clear(); 
 * 
 *     for (unsigned int k = 0; k < fe.n_dofs_per_cell(); ++k) 
 *       { 
 *         const unsigned int k_group = fe.system_to_base_index(k).first.first; 
 *         if (k_group == u_dof) 
 *           element_indices_u.push_back(k); 
 *         else if (k_group == p_dof) 
 *           element_indices_p.push_back(k); 
 *         else if (k_group == J_dof) 
 *           element_indices_J.push_back(k); 
 *         else 
 *           { 
 *             Assert(k_group <= J_dof, ExcInternalError()); 
 *           } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::setup_qph}  用于存储正交信息的方法已经在  step-18  中描述。这里我们为SMP机器实现一个类似的设置。
 * 

 * 
 * 首先，实际的QPH数据对象被创建。这必须在网格被细化到最细的程度后才能完成。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::setup_qph() 
 *   { 
 *     std::cout << "    Setting up quadrature point data..." << std::endl; 
 * 
 *     quadrature_point_history.initialize(triangulation.begin_active(), 
 *                                         triangulation.end(), 
 *                                         n_q_points); 
 * 
 * @endcode
 * 
 * 接下来我们设置初始正交点数据。请注意，当检索正交点数据时，它将作为一个智能指针的向量返回。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       { 
 *         const std::vector<std::shared_ptr<PointHistory<dim>>> lqph = 
 *           quadrature_point_history.get_data(cell); 
 *         Assert(lqph.size() == n_q_points, ExcInternalError()); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           lqph[q_point]->setup_lqp(parameters); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::update_qph_incremental}  由于QP信息的更新经常发生，并且涉及一些昂贵的操作，我们定义了一个多线程的方法，将任务分布在一些CPU核心上。
 * 

 * 
 * 要开始这样做，首先我们需要获得这个牛顿增量时的总解，然后创建初始的从头开始的副本和复制数据对象。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta) 
 *   { 
 *     timer.enter_subsection("Update QPH data"); 
 *     std::cout << " UQPH " << std::flush; 
 * 
 *     const BlockVector<double> solution_total( 
 *       get_total_solution(solution_delta)); 
 * 
 *     const UpdateFlags uf_UQPH(update_values | update_gradients); 
 *     PerTaskData_UQPH  per_task_data_UQPH; 
 *     ScratchData_UQPH  scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total); 
 * 
 * @endcode
 * 
 * 然后，我们将它们和单格更新函数传递给WorkStream进行处理。
 * 

 * 
 * 
 * @code
 *     WorkStream::run(dof_handler.active_cell_iterators(), 
 *                     *this, 
 *                     &Solid::update_qph_incremental_one_cell, 
 *                     &Solid::copy_local_to_global_UQPH, 
 *                     scratch_data_UQPH, 
 *                     per_task_data_UQPH); 
 * 
 *     timer.leave_subsection(); 
 *   } 
 * 
 * @endcode
 * 
 * 现在我们描述一下我们如何从解决方案向量中提取数据，并将其传递给每个QP存储对象进行处理。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::update_qph_incremental_one_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     ScratchData_UQPH &                                    scratch, 
 *     PerTaskData_UQPH & /*data*/) 
 *   { 
 *     const std::vector<std::shared_ptr<PointHistory<dim>>> lqph = 
 *       quadrature_point_history.get_data(cell); 
 *     Assert(lqph.size() == n_q_points, ExcInternalError()); 
 * 
 *     Assert(scratch.solution_grads_u_total.size() == n_q_points, 
 *            ExcInternalError()); 
 *     Assert(scratch.solution_values_p_total.size() == n_q_points, 
 *            ExcInternalError()); 
 *     Assert(scratch.solution_values_J_total.size() == n_q_points, 
 *            ExcInternalError()); 
 * 
 *     scratch.reset(); 
 * 
 * @endcode
 * 
 * 我们首先需要找到当前单元内正交点的数值和梯度，然后利用位移梯度和总压力及扩张解数值更新每个局部QP。
 * 

 * 
 * 
 * @code
 *     scratch.fe_values.reinit(cell); 
 *     scratch.fe_values[u_fe].get_function_gradients( 
 *       scratch.solution_total, scratch.solution_grads_u_total); 
 *     scratch.fe_values[p_fe].get_function_values( 
 *       scratch.solution_total, scratch.solution_values_p_total); 
 *     scratch.fe_values[J_fe].get_function_values( 
 *       scratch.solution_total, scratch.solution_values_J_total); 
 * 
 *     for (const unsigned int q_point : 
 *          scratch.fe_values.quadrature_point_indices()) 
 *       lqph[q_point]->update_values(scratch.solution_grads_u_total[q_point], 
 *                                    scratch.solution_values_p_total[q_point], 
 *                                    scratch.solution_values_J_total[q_point]); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidsolve_nonlinear_timestep"></a> 
 * <h4>Solid::solve_nonlinear_timestep</h4>
 * 

 * 
 * 下一个函数是牛顿-拉弗逊方案的驱动方法。在它的顶部，我们创建一个新的向量来存储当前的牛顿更新步骤，重置错误存储对象并打印求解器头。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::solve_nonlinear_timestep(BlockVector<double> &solution_delta) 
 *   { 
 *     std::cout << std::endl 
 *               << "Timestep " << time.get_timestep() << " @ " << time.current() 
 *               << "s" << std::endl; 
 * 
 *     BlockVector<double> newton_update(dofs_per_block); 
 * 
 *     error_residual.reset(); 
 *     error_residual_0.reset(); 
 *     error_residual_norm.reset(); 
 *     error_update.reset(); 
 *     error_update_0.reset(); 
 *     error_update_norm.reset(); 
 * 
 *     print_conv_header(); 
 * 
 * @endcode
 * 
 * 我们现在进行一些牛顿迭代来迭代解决这个非线性问题。 由于问题是完全非线性的，而且我们使用的是完全牛顿方法，所以存储在切线矩阵和右手边向量中的数据是不能重复使用的，必须在每个牛顿步骤中清除。然后，我们最初建立线性系统并检查收敛性（并在第一次迭代中存储这个值）。rhs向量的无约束DOF持有失衡的力，并共同决定是否达到了平衡解。
 * 

 * 
 * 尽管对于这个特定的问题，我们可以在组合系统矩阵之前构建RHS向量，但为了扩展性，我们选择不这样做。分别组装RHS向量和系统矩阵的好处是，后者是一个昂贵的操作，我们可以通过在达到收敛时不组装切线矩阵来避免一个额外的组装过程。然而，这使得使用MPI并行化代码变得更加困难。此外，当把问题扩展到瞬态情况时，由于时间离散化和对速度和加速度场的约束应用，可能会对RHS产生额外的贡献。
 * 

 * 
 * 
 * @code
 *     unsigned int newton_iteration = 0; 
 *     for (; newton_iteration < parameters.max_iterations_NR; ++newton_iteration) 
 *       { 
 *         std::cout << " " << std::setw(2) << newton_iteration << " " 
 *                   << std::flush; 
 * 
 * @endcode
 * 
 * 我们构建线性系统，但暂不求解它（这一步应该比装配要贵得多）。
 * 

 * 
 * 
 * @code
 *         make_constraints(newton_iteration); 
 *         assemble_system(); 
 * 
 * @endcode
 * 
 * 我们现在可以确定归一化剩余误差，并检查解决方案的收敛性。
 * 

 * 
 * 
 * @code
 *         get_error_residual(error_residual); 
 *         if (newton_iteration == 0) 
 *           error_residual_0 = error_residual; 
 * 
 *         error_residual_norm = error_residual; 
 *         error_residual_norm.normalize(error_residual_0); 
 * 
 *         if (newton_iteration > 0 && error_update_norm.u <= parameters.tol_u && 
 *             error_residual_norm.u <= parameters.tol_f) 
 *           { 
 *             std::cout << " CONVERGED! " << std::endl; 
 *             print_conv_footer(); 
 * 
 *             break; 
 *           } 
 * 
 * @endcode
 * 
 * 如果我们决定要继续迭代，我们就解决线性化系统。
 * 

 * 
 * 
 * @code
 *         const std::pair<unsigned int, double> lin_solver_output = 
 *           solve_linear_system(newton_update); 
 * 
 * @endcode
 * 
 * 我们现在可以确定归一化的牛顿更新误差。
 * 

 * 
 * 
 * @code
 *         get_error_update(newton_update, error_update); 
 *         if (newton_iteration == 0) 
 *           error_update_0 = error_update; 
 * 
 *         error_update_norm = error_update; 
 *         error_update_norm.normalize(error_update_0); 
 * 
 * @endcode
 * 
 * 最后，由于我们隐含地接受了求解步骤，我们可以对当前时间步骤的求解增量进行实际更新，更新与这个新位移和应力状态有关的所有正交点信息，并继续迭代。
 * 

 * 
 * 
 * @code
 *         solution_delta += newton_update; 
 *         update_qph_incremental(solution_delta); 
 * 
 *         std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7) 
 *                   << std::scientific << lin_solver_output.first << "  " 
 *                   << lin_solver_output.second << "  " 
 *                   << error_residual_norm.norm << "  " << error_residual_norm.u 
 *                   << "  " << error_residual_norm.p << "  " 
 *                   << error_residual_norm.J << "  " << error_update_norm.norm 
 *                   << "  " << error_update_norm.u << "  " << error_update_norm.p 
 *                   << "  " << error_update_norm.J << "  " << std::endl; 
 *       } 
 * 
 * @endcode
 * 
 * 在最后，如果发现我们事实上做了比参数文件允许的更多的迭代，我们会引发一个异常，可以在main()函数中捕获。调用<code>AssertThrow(condition, exc_object)</code>实质上等同于<code>if (!cond) throw exc_object;</code>，但前一种形式在异常对象中填充了某些字段，以确定异常发生的位置（文件名和行号），使之更容易识别问题发生的位置。
 * 

 * 
 * 
 * @code
 *     AssertThrow(newton_iteration < parameters.max_iterations_NR, 
 *                 ExcMessage("No convergence in nonlinear solver!")); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidprint_conv_headerandSolidprint_conv_footer"></a> 
 * <h4>Solid::print_conv_header and Solid::print_conv_footer</h4>
 * 

 * 
 * 这个程序在一个漂亮的表格中打印出数据，这个表格在每次迭代的基础上被更新。接下来的两个函数设置了表头和表脚。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::print_conv_header() 
 *   { 
 *     static const unsigned int l_width = 150; 
 * 
 *     for (unsigned int i = 0; i < l_width; ++i) 
 *       std::cout << "_"; 
 *     std::cout << std::endl; 
 * 
 *     std::cout << "               SOLVER STEP               " 
 *               << " |  LIN_IT   LIN_RES    RES_NORM    " 
 *               << " RES_U     RES_P      RES_J     NU_NORM     " 
 *               << " NU_U       NU_P       NU_J " << std::endl; 
 * 
 *     for (unsigned int i = 0; i < l_width; ++i) 
 *       std::cout << "_"; 
 *     std::cout << std::endl; 
 *   } 
 * 
 *   template <int dim> 
 *   void Solid<dim>::print_conv_footer() 
 *   { 
 *     static const unsigned int l_width = 150; 
 * 
 *     for (unsigned int i = 0; i < l_width; ++i) 
 *       std::cout << "_"; 
 *     std::cout << std::endl; 
 * 
 *     const std::pair<double, double> error_dil = get_error_dilation(); 
 * 
 *     std::cout << "Relative errors:" << std::endl 
 *               << "Displacement:\t" << error_update.u / error_update_0.u 
 *               << std::endl 
 *               << "Force: \t\t" << error_residual.u / error_residual_0.u 
 *               << std::endl 
 *               << "Dilatation:\t" << error_dil.first << std::endl 
 *               << "v / V_0:\t" << error_dil.second * vol_reference << " / " 
 *               << vol_reference << " = " << error_dil.second << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_dilation"></a> 
 * <h4>Solid::get_error_dilation</h4>
 * 

 * 
 * 计算空间配置中的域的体积
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   double Solid<dim>::compute_vol_current() const 
 *   { 
 *     double vol_current = 0.0; 
 * 
 *     FEValues<dim> fe_values(fe, qf_cell, update_JxW_values); 
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 * 
 * @endcode
 * 
 * 与之前调用的不同，在这个例子中，正交点的数据是特别不可修改的，因为我们将只访问数据。我们通过将这个更新函数标记为常量来确保正确的get_data函数被调用。
 * 

 * 
 * 
 * @code
 *         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph = 
 *           quadrature_point_history.get_data(cell); 
 *         Assert(lqph.size() == n_q_points, ExcInternalError()); 
 * 
 *         for (const unsigned int q_point : fe_values.quadrature_point_indices()) 
 *           { 
 *             const double det_F_qp = lqph[q_point]->get_det_F(); 
 *             const double JxW      = fe_values.JxW(q_point); 
 * 
 *             vol_current += det_F_qp * JxW; 
 *           } 
 *       } 
 *     Assert(vol_current > 0.0, ExcInternalError()); 
 *     return vol_current; 
 *   } 
 * 
 * @endcode
 * 
 * 从 $L^2$ 的误差 $ \bigl[ \int_{\Omega_0} {[ J - \widetilde{J}]}^{2}\textrm{d}V \bigr]^{1/2}$ 中计算出扩张 $\widetilde{J}$ $J \dealcoloneq \textrm{det}\ \mathbf{F}$ 的吻合程度。我们还返回域的当前体积与参考体积的比率。这对于不可压缩介质来说是很有意义的，因为我们要检查等熵约束的执行情况。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<double, double> Solid<dim>::get_error_dilation() const 
 *   { 
 *     double dil_L2_error = 0.0; 
 * 
 *     FEValues<dim> fe_values(fe, qf_cell, update_JxW_values); 
 * 
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 * 
 *         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph = 
 *           quadrature_point_history.get_data(cell); 
 *         Assert(lqph.size() == n_q_points, ExcInternalError()); 
 * 
 *         for (const unsigned int q_point : fe_values.quadrature_point_indices()) 
 *           { 
 *             const double det_F_qp   = lqph[q_point]->get_det_F(); 
 *             const double J_tilde_qp = lqph[q_point]->get_J_tilde(); 
 *             const double the_error_qp_squared = 
 *               std::pow((det_F_qp - J_tilde_qp), 2); 
 *             const double JxW = fe_values.JxW(q_point); 
 * 
 *             dil_L2_error += the_error_qp_squared * JxW; 
 *           } 
 *       } 
 * 
 *     return std::make_pair(std::sqrt(dil_L2_error), 
 *                           compute_vol_current() / vol_reference); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_residual"></a> 
 * <h4>Solid::get_error_residual</h4>
 * 

 * 
 * 确定问题的真实残差误差。 也就是说，确定无约束自由度的残差误差。 注意，要做到这一点，我们需要忽略受约束的自由度，将这些向量分量的残差设置为零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::get_error_residual(Errors &error_residual) 
 *   { 
 *     BlockVector<double> error_res(dofs_per_block); 
 * 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
 *       if (!constraints.is_constrained(i)) 
 *         error_res(i) = system_rhs(i); 
 * 
 *     error_residual.norm = error_res.l2_norm(); 
 *     error_residual.u    = error_res.block(u_dof).l2_norm(); 
 *     error_residual.p    = error_res.block(p_dof).l2_norm(); 
 *     error_residual.J    = error_res.block(J_dof).l2_norm(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidget_error_update"></a> 
 * <h4>Solid::get_error_update</h4>
 * 

 * 
 * 确定问题的真实牛顿更新误差
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::get_error_update(const BlockVector<double> &newton_update, 
 *                                     Errors &                   error_update) 
 *   { 
 *     BlockVector<double> error_ud(dofs_per_block); 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
 *       if (!constraints.is_constrained(i)) 
 *         error_ud(i) = newton_update(i); 
 * 
 *     error_update.norm = error_ud.l2_norm(); 
 *     error_update.u    = error_ud.block(u_dof).l2_norm(); 
 *     error_update.p    = error_ud.block(p_dof).l2_norm(); 
 *     error_update.J    = error_ud.block(J_dof).l2_norm(); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="Solidget_total_solution"></a> 
 * <h4>Solid::get_total_solution</h4>
 * 

 * 
 * 这个函数提供了总解，它在任何牛顿步都有效。这是必须的，因为为了减少计算误差，总解只在时间步数结束时更新。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   BlockVector<double> Solid<dim>::get_total_solution( 
 *     const BlockVector<double> &solution_delta) const 
 *   { 
 *     BlockVector<double> solution_total(solution_n); 
 *     solution_total += solution_delta; 
 *     return solution_total; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Solidassemble_system"></a> 
 * <h4>Solid::assemble_system</h4>
 * 

 * 
 * 由于我们使用TBB进行装配，我们只需设置一份流程所需的数据结构，并将其与装配函数一起传递给WorkStream对象进行处理。请注意，我们必须确保在任何装配操作发生之前，矩阵和RHS向量被重置。此外，由于我们描述的是一个诺伊曼BC的问题，我们将需要面的法线，因此必须在面的更新标志中指定这个。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::assemble_system() 
 *   { 
 *     timer.enter_subsection("Assemble system"); 
 *     std::cout << " ASM_SYS " << std::flush; 
 * 
 *     tangent_matrix = 0.0; 
 *     system_rhs     = 0.0; 
 * 
 *     const UpdateFlags uf_cell(update_values | update_gradients | 
 *                               update_JxW_values); 
 *     const UpdateFlags uf_face(update_values | update_normal_vectors | 
 *                               update_JxW_values); 
 * 
 *     PerTaskData_ASM per_task_data(dofs_per_cell); 
 *     ScratchData_ASM scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face); 
 * 
 * @endcode
 * 
 * 这里用于向WorkStream类传递数据的语法在  step-13  中讨论。
 * 

 * 
 * 
 * @code
 *     WorkStream::run( 
 *       dof_handler.active_cell_iterators(), 
 *       [this](const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *              ScratchData_ASM &                                     scratch, 
 *              PerTaskData_ASM &                                     data) { 
 *         this->assemble_system_one_cell(cell, scratch, data); 
 *       }, 
 *       [this](const PerTaskData_ASM &data) { 
 *         this->constraints.distribute_local_to_global(data.cell_matrix, 
 *                                                      data.cell_rhs, 
 *                                                      data.local_dof_indices, 
 *                                                      tangent_matrix, 
 *                                                      system_rhs); 
 *       }, 
 *       scratch_data, 
 *       per_task_data); 
 * 
 *     timer.leave_subsection(); 
 *   } 
 * 
 * @endcode
 * 
 * 当然，我们仍然要定义如何组装单个单元的切线矩阵贡献。 我们首先需要重置和初始化一些从头开始的数据结构，并检索一些关于这个单元上DOF编号的基本信息。 我们可以预先计算单元的形状函数值和梯度。请注意，形状函数梯度是根据当前配置来定义的。 也就是  $\textrm{grad}\ \boldsymbol{\varphi} = \textrm{Grad}\ \boldsymbol{\varphi} \ \mathbf{F}^{-1}$  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::assemble_system_one_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     ScratchData_ASM &                                     scratch, 
 *     PerTaskData_ASM &                                     data) const 
 *   { 
 *     data.reset(); 
 *     scratch.reset(); 
 *     scratch.fe_values.reinit(cell); 
 *     cell->get_dof_indices(data.local_dof_indices); 
 * 
 *     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph = 
 *       quadrature_point_history.get_data(cell); 
 *     Assert(lqph.size() == n_q_points, ExcInternalError()); 
 * 
 *     for (const unsigned int q_point : 
 *          scratch.fe_values.quadrature_point_indices()) 
 *       { 
 *         const Tensor<2, dim> F_inv = lqph[q_point]->get_F_inv(); 
 *         for (const unsigned int k : scratch.fe_values.dof_indices()) 
 *           { 
 *             const unsigned int k_group = fe.system_to_base_index(k).first.first; 
 * 
 *             if (k_group == u_dof) 
 *               { 
 *                 scratch.grad_Nx[q_point][k] = 
 *                   scratch.fe_values[u_fe].gradient(k, q_point) * F_inv; 
 *                 scratch.symm_grad_Nx[q_point][k] = 
 *                   symmetrize(scratch.grad_Nx[q_point][k]); 
 *               } 
 *             else if (k_group == p_dof) 
 *               scratch.Nx[q_point][k] = 
 *                 scratch.fe_values[p_fe].value(k, q_point); 
 *             else if (k_group == J_dof) 
 *               scratch.Nx[q_point][k] = 
 *                 scratch.fe_values[J_fe].value(k, q_point); 
 *             else 
 *               Assert(k_group <= J_dof, ExcInternalError()); 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 现在我们建立本地单元刚度矩阵和RHS向量。由于全局和局部系统矩阵是对称的，我们可以利用这一特性，只建立局部矩阵的下半部分，并将其值复制到上半部分。 所以我们只组装一半的 $\mathsf{\mathbf{k}}_{uu}$  ,  $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{p}} = \mathbf{0}$  ,  $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{J}}$ 块，而整个 $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$  ,  $\mathsf{\mathbf{k}}_{u \widetilde{J}} = \mathbf{0}$  ,  $\mathsf{\mathbf{k}}_{u \widetilde{p}}$ 块被构建。
 * 

 * 
 * 在这样做的时候，我们首先从我们的正交历史对象中提取一些配置相关的变量，用于当前的正交点。
 * 

 * 
 * 
 * @code
 *     for (const unsigned int q_point : 
 *          scratch.fe_values.quadrature_point_indices()) 
 *       { 
 *         const SymmetricTensor<2, dim> tau     = lqph[q_point]->get_tau(); 
 *         const Tensor<2, dim>          tau_ns  = lqph[q_point]->get_tau(); 
 *         const SymmetricTensor<4, dim> Jc      = lqph[q_point]->get_Jc(); 
 *         const double                  det_F   = lqph[q_point]->get_det_F(); 
 *         const double                  p_tilde = lqph[q_point]->get_p_tilde(); 
 *         const double                  J_tilde = lqph[q_point]->get_J_tilde(); 
 *         const double dPsi_vol_dJ   = lqph[q_point]->get_dPsi_vol_dJ(); 
 *         const double d2Psi_vol_dJ2 = lqph[q_point]->get_d2Psi_vol_dJ2(); 
 *         const SymmetricTensor<2, dim> &I = 
 *           Physics::Elasticity::StandardTensors<dim>::I; 
 * 
 * @endcode
 * 
 * 这两个张量存储了一些预计算的数据。它们的用途将很快得到解释。
 * 

 * 
 * 
 * @code
 *         SymmetricTensor<2, dim> symm_grad_Nx_i_x_Jc; 
 *         Tensor<1, dim>          grad_Nx_i_comp_i_x_tau; 
 * 
 * @endcode
 * 
 * 接下来我们定义一些别名，使装配过程更容易操作。
 * 

 * 
 * 
 * @code
 *         const std::vector<double> &                 N = scratch.Nx[q_point]; 
 *         const std::vector<SymmetricTensor<2, dim>> &symm_grad_Nx = 
 *           scratch.symm_grad_Nx[q_point]; 
 *         const std::vector<Tensor<2, dim>> &grad_Nx = scratch.grad_Nx[q_point]; 
 *         const double                       JxW = scratch.fe_values.JxW(q_point); 
 * 
 *         for (const unsigned int i : scratch.fe_values.dof_indices()) 
 *           { 
 *             const unsigned int component_i = 
 *               fe.system_to_component_index(i).first; 
 *             const unsigned int i_group = fe.system_to_base_index(i).first.first; 
 * 
 * @endcode
 * 
 * 我们首先计算来自内力的贡献。 注意，根据rhs作为残差负数的定义，这些贡献被减去。
 * 

 * 
 * 
 * @code
 *             if (i_group == u_dof) 
 *               data.cell_rhs(i) -= (symm_grad_Nx[i] * tau) * JxW; 
 *             else if (i_group == p_dof) 
 *               data.cell_rhs(i) -= N[i] * (det_F - J_tilde) * JxW; 
 *             else if (i_group == J_dof) 
 *               data.cell_rhs(i) -= N[i] * (dPsi_vol_dJ - p_tilde) * JxW; 
 *             else 
 *               Assert(i_group <= J_dof, ExcInternalError()); 
 * 
 * @endcode
 * 
 * 在我们进入内循环之前，我们还有最后一次机会来引入一些优化。我们已经考虑到了系统的对称性，现在我们可以预先计算一些在内循环中反复应用的常用项。  我们在这里不会过分，而是将重点放在昂贵的操作上，即那些涉及等级4材料刚度张量和等级2应力张量的操作。    我们可以观察到的是，这两个张量都是以 "i "DoF为索引的形状函数梯度收缩的。这意味着，当我们在 "j "DoF上循环时，这个特殊的操作保持不变。出于这个原因，我们可以从内循环中提取这个操作，并节省许多操作，对于每个正交点和DoF索引 "i"，并在索引 "j "上重复，需要用等级4对称张量对等级2对称张量进行双重收缩，用等级2张量对等级1张量进行双重收缩。    在损失一些可读性的情况下，当使用模拟默认参数时，这个小变化将使对称系统的装配时间减少一半左右，并且随着h-细化水平的提高而变得更加显著。
 * 

 * 
 * 
 * @code
 *             if (i_group == u_dof) 
 *               { 
 *                 symm_grad_Nx_i_x_Jc    = symm_grad_Nx[i] * Jc; 
 *                 grad_Nx_i_comp_i_x_tau = grad_Nx[i][component_i] * tau_ns; 
 *               } 
 * 
 * @endcode
 * 
 * 现在我们准备计算正切矩阵的贡献。
 * 

 * 
 * 
 * @code
 *             for (const unsigned int j : 
 *                  scratch.fe_values.dof_indices_ending_at(i)) 
 *               { 
 *                 const unsigned int component_j = 
 *                   fe.system_to_component_index(j).first; 
 *                 const unsigned int j_group = 
 *                   fe.system_to_base_index(j).first.first; 
 * 
 * @endcode
 * 
 * 这就是 $\mathsf{\mathbf{k}}_{uu}$ 的贡献。它包括一个材料贡献和一个几何应力贡献，后者只沿局部矩阵对角线添加。
 * 

 * 
 * 
 * @code
 *                 if ((i_group == j_group) && (i_group == u_dof)) 
 *                   { 
 * 
 * @endcode
 * 
 * 材料贡献。
 * 

 * 
 * 
 * @code
 *                     data.cell_matrix(i, j) += symm_grad_Nx_i_x_Jc *  // 
 *                                               symm_grad_Nx[j] * JxW; // 
 * 
 * @endcode
 * 
 * 几何应力的贡献。
 * 

 * 
 * 
 * @code
 *                     if (component_i == component_j) 
 *                       data.cell_matrix(i, j) += 
 *                         grad_Nx_i_comp_i_x_tau * grad_Nx[j][component_j] * JxW; 
 *                   } 
 * 
 * @endcode
 * 
 * 接下来是 $\mathsf{\mathbf{k}}_{ \widetilde{p} u}$ 的贡献。
 * 

 * 
 * 
 * @code
 *                 else if ((i_group == p_dof) && (j_group == u_dof)) 
 *                   { 
 *                     data.cell_matrix(i, j) += N[i] * det_F *               // 
 *                                               (symm_grad_Nx[j] * I) * JxW; // 
 *                   } 
 * 
 * @endcode
 * 
 * 最后是  $\mathsf{\mathbf{k}}_{ \widetilde{J \widetilde{p}}$  和  $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$  的贡献。
 * 

 * 
 * 
 * @code
 *                 else if ((i_group == J_dof) && (j_group == p_dof)) 
 *                   data.cell_matrix(i, j) -= N[i] * N[j] * JxW; 
 *                 else if ((i_group == j_group) && (i_group == J_dof)) 
 *                   data.cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW; 
 *                 else 
 *                   Assert((i_group <= J_dof) && (j_group <= J_dof), 
 *                          ExcInternalError()); 
 *               } 
 *           } 
 *       } 
 * 
 * @endcode
 * 
 * 接下来，我们组装诺伊曼贡献。我们首先检查单元格面是否存在于施加了牵引力的边界上，如果是这样的话，就加入贡献。
 * 

 * 
 * 
 * @code
 *     for (const auto &face : cell->face_iterators()) 
 *       if (face->at_boundary() && face->boundary_id() == 6) 
 *         { 
 *           scratch.fe_face_values.reinit(cell, face); 
 * 
 *           for (const unsigned int f_q_point : 
 *                scratch.fe_face_values.quadrature_point_indices()) 
 *             { 
 *               const Tensor<1, dim> &N = 
 *                 scratch.fe_face_values.normal_vector(f_q_point); 
 * 
 * @endcode
 * 
 * 使用该正交点的面法线，我们指定参考配置中的牵引力。对于这个问题，在参考配置中应用了一个定义的压力。    假设施加的牵引力的方向不随领域的变形而变化。牵引力是用第一个Piola-Kirchhoff应力简单地定义的 $\mathbf{t} = \mathbf{P}\mathbf{N} = [p_0 \mathbf{I}] \mathbf{N} = p_0 \mathbf{N}$ 我们用时间变量来线性提升压力负荷。        请注意，我们在这里计算的对右手边向量的贡献只存在于向量的位移分量中。
 * 

 * 
 * 
 * @code
 *               static const double p0 = 
 *                 -4.0 / (parameters.scale * parameters.scale); 
 *               const double         time_ramp = (time.current() / time.end()); 
 *               const double         pressure  = p0 * parameters.p_p0 * time_ramp; 
 *               const Tensor<1, dim> traction  = pressure * N; 
 * 
 *               for (const unsigned int i : scratch.fe_values.dof_indices()) 
 *                 { 
 *                   const unsigned int i_group = 
 *                     fe.system_to_base_index(i).first.first; 
 * 
 *                   if (i_group == u_dof) 
 *                     { 
 *                       const unsigned int component_i = 
 *                         fe.system_to_component_index(i).first; 
 *                       const double Ni = 
 *                         scratch.fe_face_values.shape_value(i, f_q_point); 
 *                       const double JxW = scratch.fe_face_values.JxW(f_q_point); 
 * 
 *                       data.cell_rhs(i) += (Ni * traction[component_i]) * JxW; 
 *                     } 
 *                 } 
 *             } 
 *         } 
 * 
 * @endcode
 * 
 * 最后，我们需要将本地矩阵的下半部分复制到上半部分。
 * 

 * 
 * 
 * @code
 *     for (const unsigned int i : scratch.fe_values.dof_indices()) 
 *       for (const unsigned int j : 
 *            scratch.fe_values.dof_indices_starting_at(i + 1)) 
 *         data.cell_matrix(i, j) = data.cell_matrix(j, i); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::make_constraints}  这个问题的约束条件很容易描述。在这个特殊的例子中，边界值将被计算为牛顿算法的两次第一次迭代。一般来说，我们会在第2次迭代中建立非均质约束（也就是在后面的代码块中`apply_dirichlet_bc == true`时），在接下来的步骤中只建立相应的均质约束。虽然目前的例子只有同质约束，但以前的经验表明，一个常见的错误是在重构代码到特定用途时忘记添加额外的条件。这可能导致难以调试的错误。本着这种精神，我们选择让代码在每个牛顿步骤中执行什么操作方面更加啰嗦。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::make_constraints(const int it_nr) 
 *   { 
 * 
 * @endcode
 * 
 * 由于我们(a)处理的是牛顿迭代方法，(b)使用的是位移的增量公式，以及(c)将约束条件应用于增量位移场，所以对位移更新的任何非均质约束条件只应在第2次迭代时指定。由于该迭代后约束条件将得到完全满足，因此不需要做后续的贡献。
 * 

 * 
 * 
 * @code
 *     const bool apply_dirichlet_bc = (it_nr == 0); 
 * 
 * @endcode
 * 
 * 此外，
 * 在一个时间段内的第一次牛顿迭代之后，约束条件保持不变，只要不清除 @p constraints 对象，我们就不需要修改或重建它们。
 * 

 * 
 * 
 * @code
 *     if (it_nr > 1) 
 *       { 
 *         std::cout << " --- " << std::flush; 
 *         return; 
 *       } 
 * 
 *     std::cout << " CST " << std::flush; 
 * 
 *     if (apply_dirichlet_bc) 
 *       { 
 * 
 * @endcode
 * 
 * 在牛顿第2次迭代时，我们希望应用代表位移增量的边界条件的全套非均质和均质约束。因为一般来说，每个时间步长的约束条件可能是不同的，我们需要清除约束矩阵并完全重建它。一个例子是，如果一个表面正在加速，在这种情况下，每个时间步长的位移变化是不恒定的。
 * 

 * 
 * 
 * @code
 *         constraints.clear(); 
 * 
 * @endcode
 * 
 * 三维压痕问题的边界条件如下。在-x、-y和-z面（IDs 0,2,4）我们设置了一个对称条件，只允许平面运动，而+x和+z面（IDs 1,5）无牵引力。在这个设计好的问题中，+y面的一部分（ID 3）被设定为在x-和z-分量上没有运动。最后，如前所述，+y面的另一部分有一个施加的压力，但在x和z方向上也受到约束。
 * 

 * 
 * 在下文中，我们必须告诉函数插值的边界值应该约束解向量的哪些分量（也就是说，是x-、y-、z-位移还是它们的组合）。这是用ComponentMask对象完成的（见 @ref GlossComponentMask ），如果我们为有限元提供一个我们希望选择的分量的提取器对象，我们可以从有限元得到这些对象。为此，我们首先设置了这样的提取器对象，然后在生成相关构件掩码时使用它。
 * 

 * 
 * 
 * @code
 *         const FEValuesExtractors::Scalar x_displacement(0); 
 *         const FEValuesExtractors::Scalar y_displacement(1); 
 * 
 *         { 
 *           const int boundary_id = 0; 
 * 
 *           VectorTools::interpolate_boundary_values( 
 *             dof_handler, 
 *             boundary_id, 
 *             Functions::ZeroFunction<dim>(n_components), 
 *             constraints, 
 *             fe.component_mask(x_displacement)); 
 *         } 
 *         { 
 *           const int boundary_id = 2; 
 * 
 *           VectorTools::interpolate_boundary_values( 
 *             dof_handler, 
 *             boundary_id, 
 *             Functions::ZeroFunction<dim>(n_components), 
 *             constraints, 
 *             fe.component_mask(y_displacement)); 
 *         } 
 * 
 *         if (dim == 3) 
 *           { 
 *             const FEValuesExtractors::Scalar z_displacement(2); 
 * 
 *             { 
 *               const int boundary_id = 3; 
 * 
 *               VectorTools::interpolate_boundary_values( 
 *                 dof_handler, 
 *                 boundary_id, 
 *                 Functions::ZeroFunction<dim>(n_components), 
 *                 constraints, 
 *                 (fe.component_mask(x_displacement) | 
 *                  fe.component_mask(z_displacement))); 
 *             } 
 *             { 
 *               const int boundary_id = 4; 
 * 
 *               VectorTools::interpolate_boundary_values( 
 *                 dof_handler, 
 *                 boundary_id, 
 *                 Functions::ZeroFunction<dim>(n_components), 
 *                 constraints, 
 *                 fe.component_mask(z_displacement)); 
 *             } 
 * 
 *             { 
 *               const int boundary_id = 6; 
 * 
 *               VectorTools::interpolate_boundary_values( 
 *                 dof_handler, 
 *                 boundary_id, 
 *                 Functions::ZeroFunction<dim>(n_components), 
 *                 constraints, 
 *                 (fe.component_mask(x_displacement) | 
 *                  fe.component_mask(z_displacement))); 
 *             } 
 *           } 
 *         else 
 *           { 
 *             { 
 *               const int boundary_id = 3; 
 * 
 *               VectorTools::interpolate_boundary_values( 
 *                 dof_handler, 
 *                 boundary_id, 
 *                 Functions::ZeroFunction<dim>(n_components), 
 *                 constraints, 
 *                 (fe.component_mask(x_displacement))); 
 *             } 
 *             { 
 *               const int boundary_id = 6; 
 * 
 *               VectorTools::interpolate_boundary_values( 
 *                 dof_handler, 
 *                 boundary_id, 
 *                 Functions::ZeroFunction<dim>(n_components), 
 *                 constraints, 
 *                 (fe.component_mask(x_displacement))); 
 *             } 
 *           } 
 *       } 
 *     else 
 *       { 
 * 
 * @endcode
 * 
 * 由于所有的Dirichlet约束在牛顿第2次迭代后被完全满足，我们要确保对这些条目不做进一步的修改。这意味着我们要将所有非均质的Dirichlet约束转换成均质的约束。
 * 

 * 
 * 在这个例子中，这样做的程序是非常简单的，事实上，当只应用同质边界条件时，我们可以（也会）规避任何不必要的操作。在一个更普遍的问题中，我们应该注意悬挂节点和周期性约束，这也可能引入一些不均匀性。那么，为不同类型的约束保留不同的对象可能是有利的，一旦构建了同质Dirichlet约束，就将它们合并在一起。
 * 

 * 
 * 
 * @code
 *         if (constraints.has_inhomogeneities()) 
 *           { 
 * 
 * @endcode
 * 
 * 由于仿生约束是在上一次牛顿迭代中完成的，所以不能直接修改。所以我们需要将它们复制到另一个临时对象，并在那里进行修改。一旦我们完成了，我们将把它们转移回主 @p constraints 对象。
 * 

 * 
 * 
 * @code
 *             AffineConstraints<double> homogeneous_constraints(constraints); 
 *             for (unsigned int dof = 0; dof != dof_handler.n_dofs(); ++dof) 
 *               if (homogeneous_constraints.is_inhomogeneously_constrained(dof)) 
 *                 homogeneous_constraints.set_inhomogeneity(dof, 0.0); 
 * 
 *             constraints.clear(); 
 *             constraints.copy_from(homogeneous_constraints); 
 *           } 
 *       } 
 * 
 *     constraints.close(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::assemble_sc}  解决整个块系统有点问题，因为对 $\mathsf{\mathbf{K}}_{ \widetilde{J} \widetilde{J}}$ 块没有贡献，使其不可逆转（当使用迭代求解器时）。由于压力和扩张变量DOF是不连续的，我们可以将它们浓缩成一个较小的仅有位移的系统，然后我们将对其进行求解，随后进行后处理以检索出压力和扩张的解决方案。
 * 

 * 
 * 静态凝结过程可以在全局层面上进行，但我们需要其中一个块的逆向。然而，由于压力和扩张变量是不连续的，静态凝结（SC）操作也可以在每个单元的基础上进行，我们可以通过反转局部块来产生块对角线 $\mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}$ 块的逆。我们可以再次使用TBB来做这件事，因为每个操作都将是相互独立的。
 * 

 * 
 * 通过WorkStream类使用TBB，我们把每个元素的贡献集合起来形成 $ \mathsf{\mathbf{K}}_{\textrm{con}}= \bigl[ \mathsf{\mathbf{K}}_{uu} + \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]$ 。然后这些贡献被添加到全局刚度矩阵中。鉴于这样的描述，以下两个函数应该是清楚的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::assemble_sc() 
 *   { 
 *     timer.enter_subsection("Perform static condensation"); 
 *     std::cout << " ASM_SC " << std::flush; 
 * 
 *     PerTaskData_SC per_task_data(dofs_per_cell, 
 *                                  element_indices_u.size(), 
 *                                  element_indices_p.size(), 
 *                                  element_indices_J.size()); 
 *     ScratchData_SC scratch_data; 
 * 
 *     WorkStream::run(dof_handler.active_cell_iterators(), 
 *                     *this, 
 *                     &Solid::assemble_sc_one_cell, 
 *                     &Solid::copy_local_to_global_sc, 
 *                     scratch_data, 
 *                     per_task_data); 
 * 
 *     timer.leave_subsection(); 
 *   } 
 * 
 *   template <int dim> 
 *   void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC &data) 
 *   { 
 *     for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *       for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *         tangent_matrix.add(data.local_dof_indices[i], 
 *                            data.local_dof_indices[j], 
 *                            data.cell_matrix(i, j)); 
 *   } 
 * 
 * @endcode
 * 
 * 现在我们描述静态凝结过程。按照惯例，我们必须首先找出这个单元上的自由度有哪些全局数字，并重置一些数据结构。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::assemble_sc_one_cell( 
 *     const typename DoFHandler<dim>::active_cell_iterator &cell, 
 *     ScratchData_SC &                                      scratch, 
 *     PerTaskData_SC &                                      data) 
 *   { 
 *     data.reset(); 
 *     scratch.reset(); 
 *     cell->get_dof_indices(data.local_dof_indices); 
 * 
 * @endcode
 * 
 * 我们现在提取与当前单元相关的DFS对全局刚度矩阵的贡献。  $\widetilde{p}$ 和 $\widetilde{J}$ 插值的不连续性质意味着它们在全局水平上没有局部贡献的耦合。而 $\mathbf{u}$ 道夫则不是这样。 换句话说， $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$ 、 $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{p}}$ 和 $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$ ，当从全局刚度矩阵中提取时是元素贡献。 而 $\mathsf{\mathbf{k}}_{uu}$ 则不是这种情况。
 * 

 * 
 * 注：用小写的符号表示元素刚度矩阵。
 * 

 * 
 * 目前，与当前元素相关的dof矩阵（松散地表示为 $\mathsf{\mathbf{k}}$ ）是这样的。
 * @f{align*}
 * \begin{bmatrix}
 * \mathsf{\mathbf{k}}_{uu}  &  \mathsf{\mathbf{k}}_{u\widetilde{p}}
 * & \mathbf{0}
 * \\ \mathsf{\mathbf{k}}_{\widetilde{p}u} & \mathbf{0}  &
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}
 * \\ \mathbf{0}  &  \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}}  &
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}} \end{bmatrix}
 * @f}
 * 

 * 
 * 我们现在需要对其进行修改，使其显示为
 * @f{align*}
 * \begin{bmatrix}
 * \mathsf{\mathbf{k}}_{\textrm{con}}   &
 * \mathsf{\mathbf{k}}_{u\widetilde{p}}    & \mathbf{0}
 * \\ \mathsf{\mathbf{k}}_{\widetilde{p}u} & \mathbf{0} &
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \\ \mathbf{0} & \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}} &
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}} \end{bmatrix}
 * @f}
 * with $\mathsf{\mathbf{k}}_{\textrm{con}} = \bigl[
 * \mathsf{\mathbf{k}}_{uu} +\overline{\overline{\mathsf{\mathbf{k}}}}~
 * \bigr]$ where $               \overline{\overline{\mathsf{\mathbf{k}}}}
 * \dealcoloneq \mathsf{\mathbf{k}}_{u\widetilde{p}}
 * \overline{\mathsf{\mathbf{k}}} \mathsf{\mathbf{k}}_{\widetilde{p}u}
 * $
 * and
 * $
 * \overline{\mathsf{\mathbf{k}}} =
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{k}}_{\widetilde{J}\widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
 * $.
 * 

 * 
 * 在这一点上，我们需要注意到全局数据已经存在于 $\mathsf{\mathbf{K}}_{uu}$  ,  $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$  和  $\mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{p}}$  子块中。 因此，如果我们要修改它们，我们必须考虑到已经存在的数据（也就是说，如果需要的话，简单地添加到它或删除它）。 由于copy_local_to_global操作是一个 "+="操作，我们需要考虑到这一点
 * 

 * 
 * 特别是对于 $\mathsf{\mathbf{K}}_{uu}$ 块，这意味着从周围的单元格中加入了贡献，所以我们在操作这个块时需要小心。 我们不能直接擦除子块。
 * 

 * 
 * 我们将采用这种策略来获得我们想要的子块。
 * 

 * 
 * 

 * 
 * 

 * 
 * -  $ {\mathsf{\mathbf{k}}}_{\textrm{store}}$  : 由于我们不能访问 $\mathsf{\mathbf{k}}_{uu}$ ，但我们知道它的贡献被添加到全局 $\mathsf{\mathbf{K}}_{uu}$ 矩阵中，我们只想添加元素明智的静态凝结 $\overline{\overline{\mathsf{\mathbf{k}}}}$  。
 * 

 * 
 * 

 * 
 * 

 * 
 * -  $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$  : 类似地， $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ 存在于子块中。由于复制操作是一个+=操作，我们需要减去现有的 $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ 子矩阵，此外还需要 "添加 "我们想要替换它的东西。
 * 

 * 
 * 

 * 
 * 

 * 
 * -  $\mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}$  : 由于全局矩阵是对称的，这个块和上面那个块是一样的，我们可以简单地用 $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$ 来代替这个块。
 * 

 * 
 * 我们首先从系统矩阵中提取元素数据。因此，首先我们得到单元格的整个子块，然后提取 $\mathsf{\mathbf{k}}$ 作为与当前元素相关的道夫。
 * 

 * 
 * 
 * @code
 *     data.k_orig.extract_submatrix_from(tangent_matrix, 
 *                                        data.local_dof_indices, 
 *                                        data.local_dof_indices); 
 * 
 * @endcode
 * 
 * 接下来是 $\mathsf{\mathbf{k}}_{ \widetilde{p} u}$ 的局部矩阵 
 * $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$  和  $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$  的局部矩阵。
 * 

 * 
 * 
 * @code
 *     data.k_pu.extract_submatrix_from(data.k_orig, 
 *                                      element_indices_p, 
 *                                      element_indices_u); 
 *     data.k_pJ.extract_submatrix_from(data.k_orig, 
 *                                      element_indices_p, 
 *                                      element_indices_J); 
 *     data.k_JJ.extract_submatrix_from(data.k_orig, 
 *                                      element_indices_J, 
 *                                      element_indices_J); 
 * 
 * @endcode
 * 
 * 为了得到 $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ 的逆值，我们直接将其反转。 由于 $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ 是块状对角线，所以这个操作相对便宜。
 * 

 * 
 * 
 * @code
 *     data.k_pJ_inv.invert(data.k_pJ); 
 * 
 * @endcode
 * 
 * 现在，我们可以将凝结项添加到 $\mathsf{\mathbf{k}}_{uu}$ 块中，并将其放入单元格局部矩阵 
 * $ 
 * \mathsf{\mathbf{A}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $  中。
 * 

 * 
 * 
 * @code
 *     data.k_pJ_inv.mmult(data.A, data.k_pu); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{B}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $  
 * 
 * @code
 *     data.k_JJ.mmult(data.B, data.A); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{C}}
 * =
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $  
 * 
 * @code
 *     data.k_pJ_inv.Tmmult(data.C, data.B); 
 * @endcode
 * 
 * 
 * 
 * $
 * \overline{\overline{\mathsf{\mathbf{k}}}}
 * =
 * \mathsf{\mathbf{k}}_{u \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{k}}_{\widetilde{p} u}
 * $  
 * 
 * @code
 *     data.k_pu.Tmmult(data.k_bbar, data.C); 
 *     data.k_bbar.scatter_matrix_to(element_indices_u, 
 *                                   element_indices_u, 
 *                                   data.cell_matrix); 
 * 
 * @endcode
 * 
 * 接下来我们将 $\mathsf{\mathbf{k}}^{-1}_{ \widetilde{p} \widetilde{J}}$ 放在 $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$ 块中进行后处理。 再次注意，我们需要删除那里已经存在的贡献。
 * 

 * 
 * 
 * @code
 *     data.k_pJ_inv.add(-1.0, data.k_pJ); 
 *     data.k_pJ_inv.scatter_matrix_to(element_indices_p, 
 *                                     element_indices_J, 
 *                                     data.cell_matrix); 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::solve_linear_system}  我们现在拥有所有必要的组件，可以使用两种可能的方法之一来解决线性化系统。第一种是在元素层面上进行静态凝结，这需要对切线矩阵和RHS向量进行一些改动。另外，也可以通过在全局层面上进行凝结来解决全块系统。下面我们将实现这两种方法。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   std::pair<unsigned int, double> 
 *   Solid<dim>::solve_linear_system(BlockVector<double> &newton_update) 
 *   { 
 *     unsigned int lin_it  = 0; 
 *     double       lin_res = 0.0; 
 * 
 *     if (parameters.use_static_condensation == true) 
 *       { 
 * 
 * @endcode
 * 
 * 首先，这里是使用切线矩阵的（永久）增量的方法。对于下面的内容，回顾一下
 * @f{align*}
 * \mathsf{\mathbf{K}}_{\textrm{store}}
 * \dealcoloneq
 * \begin{bmatrix}
 * \mathsf{\mathbf{K}}_{\textrm{con}}      &
 * \mathsf{\mathbf{K}}_{u\widetilde{p}}    & \mathbf{0}
 * \\  \mathsf{\mathbf{K}}_{\widetilde{p}u}    &       \mathbf{0} &
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \\  \mathbf{0}      &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}                &
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}} \end{bmatrix} \, .
 * @f}
 * 和
 * @f{align*}
 * d \widetilde{\mathsf{\mathbf{p}}}
 * & =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}} \bigr]
 * \\ d \widetilde{\mathsf{\mathbf{J}}}
 * & =
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d
 * \mathsf{\mathbf{u}} \bigr]
 * \\ \Rightarrow d \widetilde{\mathsf{\mathbf{p}}}
 * &= \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \underbrace{\bigl[\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathsf{\mathbf{K}}}}\bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d
 * \mathsf{\mathbf{u}} \bigr]
 * @f}
 * ，从而
 * @f[
 * \underbrace{\bigl[ \mathsf{\mathbf{K}}_{uu} +
 * \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]
 * }_{\mathsf{\mathbf{K}}_{\textrm{con}}} d
 * \mathsf{\mathbf{u}}
 * =
 * \underbrace{
 * \Bigl[
 * \mathsf{\mathbf{F}}_{u}
 * - \mathsf{\mathbf{K}}_{u\widetilde{p}} \bigl[
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \overline{\mathsf{\mathbf{K}}}\mathsf{\mathbf{F}}_{\widetilde{p}}
 * \bigr]
 * \Bigr]}_{\mathsf{\mathbf{F}}_{\textrm{con}}}
 * @f]
 * 其中
 * @f[
 * \overline{\overline{\mathsf{\mathbf{K}}}} \dealcoloneq
 * \mathsf{\mathbf{K}}_{u\widetilde{p}}
 * \overline{\mathsf{\mathbf{K}}}
 * \mathsf{\mathbf{K}}_{\widetilde{p}u} \, .
 * @f]
 * 

 * 
 * 在顶部，我们分配了两个临时向量来帮助进行静态凝结，并分配了变量来存储线性求解器的迭代次数和（希望收敛的）残差。
 * 

 * 
 * 
 * @code
 *         BlockVector<double> A(dofs_per_block); 
 *         BlockVector<double> B(dofs_per_block); 
 * 
 * @endcode
 * 
 * 在这个函数的第一步，我们求解增量位移  $d\mathbf{u}$  。 为此，我们进行静态浓缩，使 
 * $\mathsf{\mathbf{K}}_{\textrm{con}}
 * = \bigl[ \mathsf{\mathbf{K}}_{uu} +
 * \overline{\overline{\mathsf{\mathbf{K}}}}~ \bigr]$ ，并将 $\mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}$ 放在原 $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$ 块中。也就是说，我们制作 $\mathsf{\mathbf{K}}_{\textrm{store}}$  。
 * 

 * 
 * 
 * @code
 *         { 
 *           assemble_sc(); 
 * 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $ 
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof) 
 *             .vmult(A.block(J_dof), system_rhs.block(p_dof)); 
 * 
 * @endcode
 * 
 * $
 * \mathsf{\mathbf{B}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $  
 * 
 * @code
 *           tangent_matrix.block(J_dof, J_dof) 
 *             .vmult(B.block(J_dof), A.block(J_dof)); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * $  
 * 
 * @code
 *           A.block(J_dof) = system_rhs.block(J_dof); 
 *           A.block(J_dof) -= B.block(J_dof); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $  
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof) 
 *             .Tvmult(A.block(p_dof), A.block(J_dof)); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{u}
 * =
 * \mathsf{\mathbf{K}}_{u \widetilde{p}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $  
 * 
 * @code
 *           tangent_matrix.block(u_dof, p_dof) 
 *             .vmult(A.block(u_dof), A.block(p_dof)); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{F}}_{\text{con}}
 * =
 * \mathsf{\mathbf{F}}_{u}
 * -
 * \mathsf{\mathbf{K}}_{u \widetilde{p}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * ]
 * $  
 * 
 * @code
 *           system_rhs.block(u_dof) -= A.block(u_dof); 
 * 
 *           timer.enter_subsection("Linear solver"); 
 *           std::cout << " SLV " << std::flush; 
 *           if (parameters.type_lin == "CG") 
 *             { 
 *               const auto solver_its = static_cast<unsigned int>( 
 *                 tangent_matrix.block(u_dof, u_dof).m() * 
 *                 parameters.max_iterations_lin); 
 *               const double tol_sol = 
 *                 parameters.tol_lin * system_rhs.block(u_dof).l2_norm(); 
 * 
 *               SolverControl solver_control(solver_its, tol_sol); 
 * 
 *               GrowingVectorMemory<Vector<double>> GVM; 
 *               SolverCG<Vector<double>> solver_CG(solver_control, GVM); 
 * 
 * @endcode
 * 
 * 我们默认选择了SSOR预处理程序，因为在单线程机器上，它似乎为这个问题提供了最快的求解器收敛特性。 然而，对于不同的问题规模，这可能不是真的。
 * 

 * 
 * 
 * @code
 *               PreconditionSelector<SparseMatrix<double>, Vector<double>> 
 *                 preconditioner(parameters.preconditioner_type, 
 *                                parameters.preconditioner_relaxation); 
 *               preconditioner.use_matrix(tangent_matrix.block(u_dof, u_dof)); 
 * 
 *               solver_CG.solve(tangent_matrix.block(u_dof, u_dof), 
 *                               newton_update.block(u_dof), 
 *                               system_rhs.block(u_dof), 
 *                               preconditioner); 
 * 
 *               lin_it  = solver_control.last_step(); 
 *               lin_res = solver_control.last_value(); 
 *             } 
 *           else if (parameters.type_lin == "Direct") 
 *             { 
 * 
 * @endcode
 * 
 * 否则，如果问题足够小，可以利用直接求解器。
 * 

 * 
 * 
 * @code
 *               SparseDirectUMFPACK A_direct; 
 *               A_direct.initialize(tangent_matrix.block(u_dof, u_dof)); 
 *               A_direct.vmult(newton_update.block(u_dof), 
 *                              system_rhs.block(u_dof)); 
 * 
 *               lin_it  = 1; 
 *               lin_res = 0.0; 
 *             } 
 *           else 
 *             Assert(false, ExcMessage("Linear solver type not implemented")); 
 * 
 *           timer.leave_subsection(); 
 *         } 
 * 
 * @endcode
 * 
 * 现在我们有了位移更新，将约束分配回牛顿更新。
 * 

 * 
 * 
 * @code
 *         constraints.distribute(newton_update); 
 * 
 *         timer.enter_subsection("Linear solver postprocessing"); 
 *         std::cout << " PP " << std::flush; 
 * 
 * @endcode
 * 
 * 解决位移问题后的下一步是进行后处理，从置换中得到扩张解。     
 * $
 * d \widetilde{\mathsf{\mathbf{J}}}
 * = \mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * - \mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * \bigr]
 * $  
 * 
 * @code
 *         { 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $  
 * 
 * @code
 *           tangent_matrix.block(p_dof, u_dof) 
 *             .vmult(A.block(p_dof), newton_update.block(u_dof)); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $  
 * 
 * @code
 *           A.block(p_dof) *= -1.0; 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{p}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * $  
 * 
 * @code
 *           A.block(p_dof) += system_rhs.block(p_dof); 
 * @endcode
 * 
 * 
 * 
 * $
 * d\mathsf{\mathbf{\widetilde{J}}}
 * =
 * \mathsf{\mathbf{K}}^{-1}_{\widetilde{p}\widetilde{J}}
 * [
 * \mathsf{\mathbf{F}}_{\widetilde{p}}
 * -\mathsf{\mathbf{K}}_{\widetilde{p}u} d \mathsf{\mathbf{u}}
 * ]
 * $  
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof) 
 *             .vmult(newton_update.block(J_dof), A.block(p_dof)); 
 *         } 
 * 
 * @endcode
 * 
 * 我们在此确保任何迪里希特约束都分布在更新的解决方案上。
 * 

 * 
 * 
 * @code
 *         constraints.distribute(newton_update); 
 * 
 * @endcode
 * 
 * 最后我们用代入法求解压力的更新。     
 * $
 * d \widetilde{\mathsf{\mathbf{p}}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * - \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * \bigr]
 * $  
 * 
 * @code
 *         { 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $  
 * 
 * @code
 *           tangent_matrix.block(J_dof, J_dof) 
 *             .vmult(A.block(J_dof), newton_update.block(J_dof)); 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * -\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $  
 * 
 * @code
 *           A.block(J_dof) *= -1.0; 
 * @endcode
 * 
 * 
 * 
 * $
 * \mathsf{\mathbf{A}}_{\widetilde{J}}
 * =
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * -
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * $  
 * 
 * @code
 *           A.block(J_dof) += system_rhs.block(J_dof); 
 * 
 * @endcode
 * 
 * 和
 * 最后....      
 * 

 * 
 * $
 * d \widetilde{\mathsf{\mathbf{p}}}
 * =
 * \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
 * \bigl[
 * \mathsf{\mathbf{F}}_{\widetilde{J}}
 * - \mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{J}}
 * d \widetilde{\mathsf{\mathbf{J}}}
 * \bigr]
 * $  
 * 
 * @code
 *           tangent_matrix.block(p_dof, J_dof) 
 *             .Tvmult(newton_update.block(p_dof), A.block(J_dof)); 
 *         } 
 * 
 * @endcode
 * 
 * 我们现在已经到了终点，所以我们将所有受限的道夫分配到牛顿更新中。
 * 

 * 
 * 
 * @code
 *         constraints.distribute(newton_update); 
 * 
 *         timer.leave_subsection(); 
 *       } 
 *     else 
 *       { 
 *         std::cout << " ------ " << std::flush; 
 * 
 *         timer.enter_subsection("Linear solver"); 
 *         std::cout << " SLV " << std::flush; 
 * 
 *         if (parameters.type_lin == "CG") 
 *           { 
 * 
 * @endcode
 * 
 * 在局部水平上手动凝结扩张和压力场，以及随后的后处理，需要花费相当大的努力才能实现。简而言之，我们必须产生逆矩阵 $\mathsf{\mathbf{K}}_{\widetilde{p}\widetilde{J}}^{-1}$ ，并将其永久写入全局切线矩阵中。然后我们对 $\mathsf{\mathbf{K}}_{uu}$ 进行永久修改，产生 $\mathsf{\mathbf{K}}_{\textrm{con}}$ 。这涉及到对切线矩阵的局部子块的提取和操作。在对位移进行求解后，对扩张和压力进行求解所需的各个矩阵-向量操作被仔细地执行。将这些众多的步骤与使用LinearOperator类提供的功能进行的更简单、更透明的实现形成对比。
 * 

 * 
 * 为了便于以后使用，我们为RHS向量中的块定义了一些别名
 * 

 * 
 * 
 * @code
 *             const Vector<double> &f_u = system_rhs.block(u_dof); 
 *             const Vector<double> &f_p = system_rhs.block(p_dof); 
 *             const Vector<double> &f_J = system_rhs.block(J_dof); 
 * 
 * @endcode
 * 
 * ... 对于牛顿更新向量中的块。
 * 

 * 
 * 
 * @code
 *             Vector<double> &d_u = newton_update.block(u_dof); 
 *             Vector<double> &d_p = newton_update.block(p_dof); 
 *             Vector<double> &d_J = newton_update.block(J_dof); 
 * 
 * @endcode
 * 
 * 我们将利用系统的对称性，所以不是所有的块都需要。
 * 

 * 
 * 
 * @code
 *             const auto K_uu = 
 *               linear_operator(tangent_matrix.block(u_dof, u_dof)); 
 *             const auto K_up = 
 *               linear_operator(tangent_matrix.block(u_dof, p_dof)); 
 *             const auto K_pu = 
 *               linear_operator(tangent_matrix.block(p_dof, u_dof)); 
 *             const auto K_Jp = 
 *               linear_operator(tangent_matrix.block(J_dof, p_dof)); 
 *             const auto K_JJ = 
 *               linear_operator(tangent_matrix.block(J_dof, J_dof)); 
 * 
 * @endcode
 * 
 * 然后我们构建一个LinearOperator，代表（方形块） $\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}$ 的逆。由于它是对角线的（或者，当使用高阶分解时，几乎是对角线的），所以雅可比预处理器是合适的。
 * 

 * 
 * 
 * @code
 *             PreconditionSelector<SparseMatrix<double>, Vector<double>> 
 *               preconditioner_K_Jp_inv("jacobi"); 
 *             preconditioner_K_Jp_inv.use_matrix( 
 *               tangent_matrix.block(J_dof, p_dof)); 
 *             ReductionControl solver_control_K_Jp_inv( 
 *               static_cast<unsigned int>(tangent_matrix.block(J_dof, p_dof).m() * 
 *                                         parameters.max_iterations_lin), 
 *               1.0e-30, 
 *               parameters.tol_lin); 
 *             SolverSelector<Vector<double>> solver_K_Jp_inv; 
 *             solver_K_Jp_inv.select("cg"); 
 *             solver_K_Jp_inv.set_control(solver_control_K_Jp_inv); 
 *             const auto K_Jp_inv = 
 *               inverse_operator(K_Jp, solver_K_Jp_inv, preconditioner_K_Jp_inv); 
 * 
 * @endcode
 * 
 * 现在我们可以构建 $\mathsf{\mathbf{K}}_{\widetilde{J}\widetilde{p}}^{-1}$ 的那个转置和一个线性算子，它代表了浓缩的操作 $\overline{\mathsf{\mathbf{K}}}$ 和 $\overline{\overline{\mathsf{\mathbf{K}}}}$ 以及最后的增强矩阵 $\mathsf{\mathbf{K}}_{\textrm{con}}$  。  请注意，schur_complement()算子在这里也能派上用场，但为了清楚起见，也为了展示线性求解方案的表述和实现之间的相似性，我们将手动执行这些操作。
 * 

 * 
 * 
 * @code
 *             const auto K_pJ_inv     = transpose_operator(K_Jp_inv); 
 *             const auto K_pp_bar     = K_Jp_inv * K_JJ * K_pJ_inv; 
 *             const auto K_uu_bar_bar = K_up * K_pp_bar * K_pu; 
 *             const auto K_uu_con     = K_uu + K_uu_bar_bar; 
 * 
 * @endcode
 * 
 * 最后，我们定义了一个增强刚度矩阵的逆运算，即  $\mathsf{\mathbf{K}}_{\textrm{con}}^{-1}$  。请注意，增强刚度矩阵的预处理程序与我们使用静态凝结的情况不同。在这种情况下，预处理程序是基于未修改的 $\mathsf{\mathbf{K}}_{uu}$ ，而在第一种方法中，我们实际上修改了这个子块的条目。然而，由于 $\mathsf{\mathbf{K}}_{\textrm{con}}$ 和 $\mathsf{\mathbf{K}}_{uu}$ 在同一空间操作，它对这个问题仍然足够。
 * 

 * 
 * 
 * @code
 *             PreconditionSelector<SparseMatrix<double>, Vector<double>> 
 *               preconditioner_K_con_inv(parameters.preconditioner_type, 
 *                                        parameters.preconditioner_relaxation); 
 *             preconditioner_K_con_inv.use_matrix( 
 *               tangent_matrix.block(u_dof, u_dof)); 
 *             ReductionControl solver_control_K_con_inv( 
 *               static_cast<unsigned int>(tangent_matrix.block(u_dof, u_dof).m() * 
 *                                         parameters.max_iterations_lin), 
 *               1.0e-30, 
 *               parameters.tol_lin); 
 *             SolverSelector<Vector<double>> solver_K_con_inv; 
 *             solver_K_con_inv.select("cg"); 
 *             solver_K_con_inv.set_control(solver_control_K_con_inv); 
 *             const auto K_uu_con_inv = 
 *               inverse_operator(K_uu_con, 
 *                                solver_K_con_inv, 
 *                                preconditioner_K_con_inv); 
 * 
 * @endcode
 * 
 * 现在我们可以对位移场进行求解了。  我们可以嵌套线性运算，结果立即写入牛顿更新向量中。  很明显，这个实现密切模仿了介绍中所说的推导。
 * 

 * 
 * 
 * @code
 *             d_u = 
 *               K_uu_con_inv * (f_u - K_up * (K_Jp_inv * f_J - K_pp_bar * f_p)); 
 * 
 *             timer.leave_subsection(); 
 * 
 * @endcode
 * 
 * 需要对扩张场和压力场进行后处理的操作，也同样容易表达。
 * 

 * 
 * 
 * @code
 *             timer.enter_subsection("Linear solver postprocessing"); 
 *             std::cout << " PP " << std::flush; 
 * 
 *             d_J = K_pJ_inv * (f_p - K_pu * d_u); 
 *             d_p = K_Jp_inv * (f_J - K_JJ * d_J); 
 * 
 *             lin_it  = solver_control_K_con_inv.last_step(); 
 *             lin_res = solver_control_K_con_inv.last_value(); 
 *           } 
 *         else if (parameters.type_lin == "Direct") 
 *           { 
 * 
 * @endcode
 * 
 * 用直接求解器求解全块系统。由于它是相对稳健的，它可能对因零 $\mathsf{\mathbf{K}}_{ \widetilde{J} \widetilde{J}}$ 块的存在而产生的问题免疫。
 * 

 * 
 * 
 * @code
 *             SparseDirectUMFPACK A_direct; 
 *             A_direct.initialize(tangent_matrix); 
 *             A_direct.vmult(newton_update, system_rhs); 
 * 
 *             lin_it  = 1; 
 *             lin_res = 0.0; 
 * 
 *             std::cout << " -- " << std::flush; 
 *           } 
 *         else 
 *           Assert(false, ExcMessage("Linear solver type not implemented")); 
 * 
 *         timer.leave_subsection(); 
 * 
 * @endcode
 * 
 * 最后，我们在这里再次确保任何Dirichlet约束都分布在更新的解决方案上。
 * 

 * 
 * 
 * @code
 *         constraints.distribute(newton_update); 
 *       } 
 * 
 *     return std::make_pair(lin_it, lin_res); 
 *   } 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect4{Solid::output_results}  这里我们介绍如何将结果写入文件，以便用ParaView或Visi来查看。该方法与以前的教程中的方法类似，因此将不作详细讨论。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void Solid<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     std::vector<std::string> solution_name(dim, "displacement"); 
 *     solution_name.emplace_back("pressure"); 
 *     solution_name.emplace_back("dilatation"); 
 * 
 *     DataOutBase::VtkFlags output_flags; 
 *     output_flags.write_higher_order_cells = true; 
 *     data_out.set_flags(output_flags); 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution_n, 
 *                              solution_name, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 * 
 * @endcode
 * 
 * 由于我们处理的是一个大的变形问题，如果能在一个位移的网格上显示结果就更好了!  与DataOut类相连的MappingQEulerian类提供了一个接口，通过该接口可以实现这一目的，而不需要我们自己物理地移动三角测量对象中的网格点。 我们首先需要将解决方案复制到一个临时矢量，然后创建欧拉映射。我们还向DataOut对象指定了多项式的度数，以便在使用高阶多项式时产生一个更精细的输出数据集。
 * 

 * 
 * 
 * @code
 *     Vector<double> soln(solution_n.size()); 
 *     for (unsigned int i = 0; i < soln.size(); ++i) 
 *       soln(i) = solution_n(i); 
 *     MappingQEulerian<dim> q_mapping(degree, dof_handler, soln); 
 *     data_out.build_patches(q_mapping, degree); 
 * 
 *     std::ofstream output("solution-" + std::to_string(dim) + "d-" + 
 *                          std::to_string(time.get_timestep()) + ".vtu"); 
 *     data_out.write_vtu(output); 
 *   } 
 * 
 * } // namespace Step44 
 * @endcode
 * 
 * 
 * <a name=""></a> 
 * @sect3{Main function}  最后我们提供了主要的驱动函数，它看起来与其他教程没有什么不同。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   using namespace Step44; 
 * 
 *   try 
 *     { 
 *       const unsigned int dim = 3; 
 *       Solid<dim>         solid("parameters.prm"); 
 *       solid.run(); 
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
examples/step-44/doc/results.dox



<a name="Results"></a><h1>Results</h1>


首先，我们提出了一系列3维结果与文献中的结果的比较（见Reese等人(2000)），以证明该程序按预期工作。

我们首先比较了 $Q_1-DGPM_0-DGPM_0$ 和 $Q_2-DGPM_1-DGPM_1$ 公式的网格细化的收敛性，如下图所总结的。块的上表面的中点的垂直位移被用来评估收敛性。对于不同的载荷参数 $p/p_0$ 值，两种方案都表现出良好的收敛特性。这些结果与文献中的结果一致。低阶公式通常高估了低层次细化的位移，而高阶插值方案则低估了位移，但程度较轻。这个基准，以及其他一系列没有在这里显示的基准，使我们相信代码在正常工作。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_convergence.png" alt="">
	<p align="center">
        Convergence of the $Q_1-DGPM_0-DGPM_0$ formulation in 3-d.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_convergence.png" alt="">
	<p align="center">
        Convergence of the $Q_2-DGPM_1-DGPM_1$ formulation in 3-d.
	</p>
    </td>
  </tr>
</table> 


下面是运行该问题产生的典型屏幕输出。所展示的特殊情况是 $Q_2-DGPM_1-DGPM_1$ 公式的情况。很明显，使用Newton-Raphson方法，可以得到二次收敛的解决方案。在所有的时间步长中，解的收敛是在5个牛顿增量内实现的。收敛后的位移的 $L_2$ -norm比几何尺度小几个数量级。

@code
Grid:
	 Reference volume: 1e-09
Triangulation:
	 Number of active cells: 64
	 Number of degrees of freedom: 2699
    Setting up quadrature point data...


Timestep 1 @ 0.1s
___________________________________________________________________________________________________________________________________________________________
                 SOLVER STEP                   |  LIN_IT   LIN_RES    RES_NORM     RES_U     RES_P      RES_J     NU_NORM      NU_U       NU_P       NU_J
___________________________________________________________________________________________________________________________________________________________
  0  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     786  2.118e-06  1.000e+00  1.000e+00  0.000e+00  0.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00
  1  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     552  1.031e-03  8.563e-02  8.563e-02  9.200e-13  3.929e-08  1.060e-01  3.816e-02  1.060e-01  1.060e-01
  2  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     667  5.602e-06  2.482e-03  2.482e-03  3.373e-15  2.982e-10  2.936e-03  2.053e-04  2.936e-03  2.936e-03
  3  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     856  6.469e-10  2.129e-06  2.129e-06  2.245e-19  1.244e-13  1.887e-06  7.289e-07  1.887e-06  1.887e-06
  4  ASM_R  CONVERGED!
___________________________________________________________________________________________________________________________________________________________
Relative errors:
Displacement:	7.289e-07
Force: 		2.451e-10
Dilatation:	1.353e-07
v / V_0:	1.000e-09 / 1.000e-09 = 1.000e+00



[...]


Timestep 10 @ 1.000e+00s
___________________________________________________________________________________________________________________________________________________________
                 SOLVER STEP                   |  LIN_IT   LIN_RES    RES_NORM     RES_U     RES_P      RES_J     NU_NORM      NU_U       NU_P       NU_J
___________________________________________________________________________________________________________________________________________________________
  0  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     874  2.358e-06  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00  1.000e+00
  1  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     658  2.942e-04  1.544e-01  1.544e-01  1.208e+13  1.855e+06  6.014e-02  7.398e-02  6.014e-02  6.014e-02
  2  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     790  2.206e-06  2.908e-03  2.908e-03  7.302e+10  2.067e+03  2.716e-03  1.433e-03  2.716e-03  2.717e-03
  3  ASM_R  ASM_K  CST  ASM_SC  SLV  PP  UQPH  |     893  2.374e-09  1.919e-06  1.919e-06  4.527e+07  4.100e+00  1.672e-06  6.842e-07  1.672e-06  1.672e-06
  4  ASM_R  CONVERGED!
___________________________________________________________________________________________________________________________________________________________
Relative errors:
Displacement:	6.842e-07
Force: 		8.995e-10
Dilatation:	1.528e-06
v / V_0:	1.000e-09 / 1.000e-09 = 1.000e+00
@endcode






使用定时器类，我们可以分辨出代码的哪些部分需要最高的计算费用。对于一个有大量自由度的案例（即高度精细化），下面给出了定时器的典型输出。本教程中的大部分代码都是基于Step-18和其他文章中描述、讨论和演示的优化而开发的。超过93%的时间花在线性求解器上，很明显，对于大型三维问题，可能有必要投资一个更好的求解器。SSOR预处理程序不是多线程的，但对于这类实体问题是有效的。研究使用另一种求解器，如通过Trilinos库提供的求解器，可能是有益的。




@code
+---------------------------------------------+------------+------------+
| Total wallclock time elapsed since start    | 9.874e+02s |            |
|                                             |            |            |
| Section                         | no. calls |  wall time | % of total |
+---------------------------------+-----------+------------+------------+
| Assemble system right-hand side |        53 | 1.727e+00s |  1.75e-01% |
| Assemble tangent matrix         |        43 | 2.707e+01s |  2.74e+00% |
| Linear solver                   |        43 | 9.248e+02s |  9.37e+01% |
| Linear solver postprocessing    |        43 | 2.743e-02s |  2.78e-03% |
| Perform static condensation     |        43 | 1.437e+01s |  1.46e+00% |
| Setup system                    |         1 | 3.897e-01s |  3.95e-02% |
| Update QPH data                 |        43 | 5.770e-01s |  5.84e-02% |
+---------------------------------+-----------+------------+------------+
@endcode




然后我们用ParaView对两种情况的结果进行了可视化。第一个是最粗的网格和最低阶插值方法。   $Q_1-DGPM_0-DGPM_0$  .第二种是在细化网格上使用 $Q_2-DGPM_1-DGPM_1$ 公式。位移的垂直分量、压力 $\widetilde{p}$ 和扩张 $\widetilde{J}$ 场显示如下。


对于第一种情况，很明显，粗略的空间离散化加上大位移导致了低质量的解决方案（加载比为 $p/p_0=80$ ）。此外，元素之间的压力差非常大。元素上的恒定压力场意味着大的压力梯度没有被捕获。然而，应该注意的是，即使在这种离散性差的情况下，在标准 $Q_1$ 位移公式中会出现的锁定现象也不会出现。块体顶面的跟踪节点的最终垂直位移仍在收敛解的12.5%以内。压力解决方案是非常粗略的，在相邻的单元之间有很大的跳跃。很明显，离施加的牵引力最近的体积经历了压缩，而域的外延则处于膨胀状态。膨胀解场和压力场明显相关，正的膨胀表示正压区域，负的表示压缩区域。正如介绍中所讨论的，压缩性压力有一个负号，而扩张性压力有一个正号。这源于体积应变能量函数的定义，与压力的物理现实的解释相反。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-displacement.png" alt="">
	<p align="center">
        Z-displacement solution for the 3-d problem.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-pressure.png" alt="">
	<p align="center">
        Discontinuous piece-wise constant pressure field.
	</p>
    </td>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q1-P0_gr_1_p_ratio_80-dilatation.png" alt="">
	<p align="center">
        Discontinuous piece-wise constant dilatation field.
	</p>
    </td>
  </tr>
</table> 

结合空间细化和高阶插值方案，产生了高质量的解决方案。三个网格细化加上 $Q_2-DGPM_1-DGPM_1$ 公式产生的结果清楚地抓住了问题的力学原理。牵引面的变形得到了很好的解决。我们现在可以观察到所施加的牵引力的实际范围，最大的力被施加在表面的中心点，导致最大的压缩。尽管领域中出现了很高的应变，特别是在施加牵引力的区域的边界，但解决方案仍然是准确的。压力场被捕捉到的细节比以前多得多。压缩和膨胀区域之间有明显的区别和过渡，压力场的线性近似允许在子元素尺度上对压力进行精细的可视化。然而，应该注意的是，压力场仍然是不连续的，可以在一个连续的网格上进行平滑处理，以达到后期处理的目的。




 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-displacement.png" alt="">
	<p align="center">
        Z-displacement solution for the 3-d problem.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-pressure.png" alt="">
	<p align="center">
        Discontinuous linear pressure field.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Q2-P1_gr_3_p_ratio_80-dilatation.png" alt="">
	<p align="center">
        Discontinuous linear dilatation field.
	</p>
    </td>
  </tr>
</table> 

这一简要的分析结果表明，三场公式能够有效地规避高度不可压缩介质的体积锁定。混合配方能够准确模拟近乎不可压缩的块体在压缩状态下的位移。命令行输出显示，在极度压缩下的体积变化导致泊松比为0.4999的体积变化小于0.01%。

在运行时间方面，对于类似的自由度数量， $Q_2-DGPM_1-DGPM_1$ 公式往往比 $Q_1-DGPM_0-DGPM_0$ 的计算成本更高（通过为低阶插值增加一个额外的网格细化级别产生）。下图显示了在一台4核（8线程）机器上连续运行的一批测试的情况。高阶方法计算时间的增加可能是由于高阶元素所需的带宽增加。如前所述，使用更好的求解器和预处理程序可以减轻使用高阶公式的费用。据观察，对于给定的问题，与单线程的SSOR预处理程序相比，使用多线程的Jacobi预处理程序可以减少72%的计算运行时间（在最坏的情况下是具有大量自由度的高阶公式）。然而，根据作者的经验，雅可比预处理方法可能不适合某些涉及替代构成模型的有限应变问题。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
     <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.Normalised_runtime.png" alt="">
	<p align="center">
        Runtime on a 4-core machine, normalised against the lowest grid resolution $Q_1-DGPM_0-DGPM_0$ solution that utilised a SSOR preconditioner.
	</p>
    </td>
  </tr>
</table> 


最后，下面展示了两个不同级别的网格细化的2维问题的位移解决方案的结果。很明显，由于二维模拟的额外约束，所产生的位移场虽然在质量上相似，但与三维模拟的情况不同。


 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.2d-gr_2.png" alt="">
	<p align="center">
        Y-displacement solution in 2-d for 2 global grid refinement levels.
	</p>
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-44.2d-gr_5.png" alt="">
	<p align="center">
        Y-displacement solution in 2-d for 5 global grid refinement levels.
	</p>
    </td>
  </tr>
</table> 

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


这项工作有许多明显的延伸。

- 首先，可以在自由能函数中加入一个额外的约束条件，以便在材料中强制执行高度的不可压缩性。一个额外的拉格朗日乘数将被引入，但这可以最容易地使用增强的拉格朗日乘数的原则来处理。这在  <em>  Simo和Taylor (1991)  </em>  中得到了证明。

- 这个模型中使用的构成关系是比较基本的。将材料类分成两个独立的类，一个处理体积响应，另一个处理等温线响应，并产生一个通用的材料类（即具有抽象的虚拟函数，派生类必须实现），允许增加更复杂的材料模型，这可能是有益的。这些模型可以包括其他超弹性材料、塑性和粘弹性材料以及其他材料。

- 该程序是为解决单节点多核机器上的问题而开发的。只要稍加努力，该程序就可以通过使用Petsc或Trilinos扩展到大规模的计算环境，使用的技术与step-40中演示的类似。这主要涉及对设置、装配、 <code>PointHistory</code> 和线性求解器例程的修改。

- 由于该程序假定为准静态平衡，为了研究惯性效应很重要的问题，例如涉及冲击的问题，有必要进行扩展以包括动态效应。

- 对于高度非线性问题，负载和解的限制程序可能是必要的。可以增加一个线搜索算法，将步长限制在牛顿增量内，以确保最佳收敛性。也可能需要使用负载限制方法，如Riks方法，来解决涉及几何不稳定性的不稳定问题，如屈曲和快穿。

- 许多物理问题涉及接触。有可能将物体间的摩擦或无摩擦接触的影响纳入这个程序。这将涉及到在自由能函数中增加一个额外的项，因此需要增加装配程序。我们还需要管理接触问题（检测和应力计算）本身。在自由能函数中增加惩罚项的一个替代方法是使用主动集方法，如步骤41中使用的方法。

- 使用LinearOperators的完整缩减程序已经被编码到线性求解器例程中。这也可以通过应用schur_complement()操作符来实现，以更自动化的方式缩减一个或多个字段。

- 最后，自适应网格细化，如步骤6和步骤18所示，可以提供额外的求解精度。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-44.cc"
*/
