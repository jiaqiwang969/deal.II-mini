/**
@page step_41 The step-41 tutorial program
This tutorial depends on step-15.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Introduction">Introduction</a>
        <li><a href="#Classicalformulation">Classical formulation</a>
        <li><a href="#Derivationofthevariationalinequality">Derivation of the variational inequality</a>
        <li><a href="#Formulationasasaddlepointproblem">Formulation as a saddle point problem</a>
        <li><a href="#ActiveSetmethodstosolvethesaddlepointproblem">Active Set methods to solve the saddle point problem</a>
        <li><a href="#Theprimaldualactivesetalgorithm">The primal-dual active set algorithm</a>
        <li><a href="#Implementation">Implementation</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeObstacleProblemcodeclasstemplate">The <code>ObstacleProblem</code> class template</a>
        <li><a href="#Righthandsideboundaryvaluesandtheobstacle">Right hand side, boundary values, and the obstacle</a>
        <li><a href="#ImplementationofthecodeObstacleProblemcodeclass">Implementation of the <code>ObstacleProblem</code> class</a>
      <ul>
        <li><a href="#ObstacleProblemObstacleProblem">ObstacleProblem::ObstacleProblem</a>
        <li><a href="#ObstacleProblemmake_grid">ObstacleProblem::make_grid</a>
        <li><a href="#ObstacleProblemsetup_system">ObstacleProblem::setup_system</a>
        <li><a href="#ObstacleProblemassemble_system">ObstacleProblem::assemble_system</a>
        <li><a href="#ObstacleProblemassemble_mass_matrix_diagonal">ObstacleProblem::assemble_mass_matrix_diagonal</a>
        <li><a href="#ObstacleProblemupdate_solution_and_constraints">ObstacleProblem::update_solution_and_constraints</a>
        <li><a href="#ObstacleProblemsolve">ObstacleProblem::solve</a>
        <li><a href="#ObstacleProblemoutput_results">ObstacleProblem::output_results</a>
        <li><a href="#ObstacleProblemrun">ObstacleProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-41/doc/intro.dox

 <br> 

<i>This program was contributed by Jörg Frohne (University of Siegen,
Germany) while on a long-term visit to Texas A&amp;M University.
<br>
This material is based upon work partly supported by ThyssenKrupp Steel Europe.
</i>


<a name="Intro"></a>

<a name="Introduction"></a><h3>Introduction</h3>


这个例子是基于二维的拉普拉斯方程，涉及的问题是，如果一个膜被一些外力偏转，但也被一个障碍物所限制，会发生什么。换句话说，想想一个弹性膜在边界处被夹在一个矩形框架上（我们选择 $\Omega =
\left[-1,1\right]^2$ ），由于重力作用而下垂。如果膜下有一个障碍物，阻止它达到平衡位置，如果重力是唯一存在的力，现在会发生什么？在目前的例子程序中，我们将考虑在膜下有一个楼梯的障碍物，重力推着膜。

这个问题通常被称为 "障碍问题"（也见<a
href="http://en.wikipedia.org/wiki/Obstacle_problem">this Wikipedia article</a>），它的结果是一个变分不等式，而不是变成弱形式的变分方程。下面我们将从经典的表述中推导出它，但在我们继续讨论数学问题之前，让我们展示一下我们在这个教程程序中要考虑的问题的解决方式，以获得一些我们应该期待的直觉。

 <table align="center" class="tutorial" cellspacing="3" cellpadding="3">
  <tr>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.png" alt="">
    </td>
    <td align="center">
        <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.png" alt="">
    </td>
  </tr>
</table> 

在这里，在左边，我们看到膜的位移。下面的障碍物的形状是清晰可见的。在右边，我们叠加了膜的哪些部分与障碍物接触。我们以后会把这组点称为 "活动集"，以表明这里有一个不等式约束在活动。




<a name="Classicalformulation"></a><h3>Classical formulation</h3>


该问题的经典表述具有以下形式。

@f{align*}


 -\textrm{div}\ \sigma &\geq f & &\quad\text{in } \Omega,\\
 \sigma &= \nabla u & &\quad\text{in } \Omega,\\
 u(\mathbf x) &= 0 & &\quad\text{on }\partial\Omega,\\
(-\Delta u - f)(u - g) &= 0 & &\quad\text{in } \Omega,\\
 u(\mathbf x) &\geq g(\mathbf x) & &\quad\text{in } \Omega


@f}

与 $u\in H^2(\Omega)$  。    $u$  是一个标量值函数，表示膜的垂直位移。第一个方程被称为平衡条件，有一个区域密度的力  $f$  。这里，我们将考虑这个力是重力。第二个方程被称为胡克定律，即应力 $\sigma$ 与位移 $u$ 的梯度成正比（比例常数，通常用 $E$ 表示，这里被设定为1，但不失一般性；如果它是常数，它可以被放入右边的函数）。在边界，我们有零迪里希特条件。很明显，前两个方程可以结合起来，得到 $-\Delta u \ge f$  。

直观地说，重力是向下作用的，所以 $f(\mathbf x)$ 是一个负函数（我们在这个程序中选择 $f=-10$ ）。那么，第一个条件意味着作用在膜上的总力是重力加上一些正值：即障碍物在它们两个接触的地方对膜施加的向上的力。这个额外的力有多大？我们还不知道（我们也不知道它实际作用的 "位置"），但它必须是使膜不穿透障碍物的。

上面的第四个等式和最后一个不等式构成了障碍条件，它必须在整个领域的每一点都成立。这两个条件中的后者意味着膜必须在任何地方都高于障碍物 $g(\mathbf x)$ 。倒数第二个方程，通常被称为 "互补条件"，说的是在膜不与障碍物接触的地方（即那些 $\mathbf x$ 的地方 $u(\mathbf x) - g(\mathbf x) \neq 0$ ），那么 $-\Delta u=f$ 在这些地方；换句话说，没有额外的力作用在那里，如预期的那样。另一方面，在 $u=g$ 的地方，我们可以有 $-\Delta u-f
\neq 0$ ，也就是说，可以有额外的力（尽管不一定要有：膜有可能只是接触而不是压住障碍物）。




<a name="Derivationofthevariationalinequality"></a><h3>Derivation of the variational inequality</h3>


获得障碍物问题的变量表述的一个明显方法是考虑总势能。

@f{equation*}
 E(u) \dealcoloneq \dfrac{1}{2}\int\limits_{\Omega} \nabla u \cdot \nabla u - \int\limits_{\Omega} fu.


@f}

我们必须找到以下最小化问题的解决方案 $u\in G$ 。

@f{equation*}
 E(u)\leq E(v)\quad \forall v\in G,


@f}

与可接受位移的凸集。

@f{equation*}
 G \dealcoloneq \lbrace v\in V: v\geq g \text{ a.e. in } \Omega\rbrace,\quad V\dealcoloneq H^1_0(\Omega).


@f}

这组数据照顾到了上述第三和第五个条件（边界值和互补条件）。

现在考虑 $E$ 的最小化器 $u\in G$ 和任何其他函数 $v\in
G$  。那么函数

@f{equation*}
 F(\varepsilon) \dealcoloneq E(u+\varepsilon(v-u)),\quad\varepsilon\in\left[0,1\right],


@f}

在 $\varepsilon = 0$ 处取最小值（因为 $u$ 是能量函数 $E(\cdot)$ 的最小值），因此，对于 $v$ 的任何选择， $F'(0)\geq 0$ 。请注意， $u+\varepsilon(v-u) = (1-\varepsilon)u+\varepsilon v\in G$  因为 $G$  的凸性。如果我们计算 $F'(\varepsilon)\vert_{\varepsilon=0}$ ，就可以得到我们要寻找的变异公式。

<i>Find a function $u\in G$ with</i>

@f{equation*}
 \left(\nabla u, \nabla(v-u)\right) \geq \left(f,v-u\right) \quad \forall v\in G.


@f}



这是变分不等式的典型形式，不仅仅是 $v$ 出现在双线性形式中，实际上还有 $v-u$  。原因是这样的：如果 $u$ 不受约束，那么我们可以在 $G$ 中找到测试函数 $v$ ，从而使 $v-u$ 可以有任何符号。通过选择测试函数 $v_1,v_2$ 使 $v_1-u = -(v_2-u)$ ，可以看出，只有当两边事实上相等时，不等式才能对 $v_1$ 和 $v_2$ 都成立，也就是说，我们得到一个变异的相等。

另一方面，如果 $u=g$ ，那么 $G$ 只允许测试函数 $v$ ，所以实际上 $v-u\ge 0$  。这意味着我们不能像上面那样用 $v-u$ 和 $-(v-u)$ 来测试这个方程，所以我们不能再得出两边实际上相等的结论。因此，这就模仿了我们上面讨论互补性条件的方式。




<a name="Formulationasasaddlepointproblem"></a><h3>Formulation as a saddle point problem</h3>


上面的变分不等式在工作中是很尴尬的。因此，我们想把它重新表述为一个等价的鞍点问题。我们引入拉格朗日乘子 $\lambda$ 和拉格朗日乘子 $K\subset V'$ 、 $V'$ 的凸锥 $V$ 、 $K \dealcoloneq \{\mu\in V': \langle\mu,v\rangle\geq 0,\quad \forall
v\in V, v \le 0 \}$ 的对偶空间，其中 $\langle\cdot,\cdot\rangle$ 表示 $V'$ 和 $V$  之间的对偶性。直观地说， $K$ 是所有 "非正函数 "的锥体，除了 $K\subset (H_0^1)'$ ，所以也包含了除正函数之外的其他对象。这就产生了。

<i>Find $u\in V$ and $\lambda\in K$ such that</i>

@f{align*}
 a(u,v) + b(v,\lambda) &= f(v),\quad &&v\in V\\
 b(u,\mu - \lambda) &\leq \langle g,\mu - \lambda\rangle,\quad&&\mu\in K,


@f}

<i>with</i>

@f{align*}
 a(u,v) &\dealcoloneq \left(\nabla u, \nabla v\right),\quad &&u,v\in V\\
 b(u,\mu) &\dealcoloneq \langle u,\mu\rangle,\quad &&u\in V,\quad\mu\in V'.


@f}

换句话说，我们可以把 $\lambda$ 看作是障碍物对膜施加的额外正力的负数。上面陈述的第二行中的不等式似乎只有错误的符号，因为我们在 $\lambda=0$ 的地方有 $\mu-\lambda<0$ ，鉴于 $K$ 的定义。

Glowinski, Lions and Tr&eacute;moli&egrave;res.中阐述了这个鞍点问题 $(u,\lambda)\in V\times K$ 的存在性和唯一性。Numerical Analysis of Variational Inequalities, North-Holland, 1981.




<a name="ActiveSetmethodstosolvethesaddlepointproblem"></a><h3>Active Set methods to solve the saddle point problem</h3>


有不同的方法来解决变量不等式。作为一种可能性，你可以把鞍点问题理解为一个带有不等式约束的凸二次方程序（QP）。

为了达到这个目的，让我们假设我们用相同的有限元空间来离散 $u$ 和 $\lambda$ ，例如通常的 $Q_k$ 空间。然后我们会得到方程

@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &[BU-G]_i \geq 0, \quad \Lambda_i \leq 0,\quad \Lambda_i[BU-G]_i = 0
\qquad \forall i.&


@f}

其中 $B$ 是所选有限元空间上的质量矩阵，上面的指数 $i$ 是针对位于域内部的自由度集合 $\cal S$ 中的所有自由度（我们在周边有迪里希条件）。然而，如果我们在组合产生这个质量矩阵的所有项时使用一个特殊的正交规则，即一个正交公式，其中正交点只位于定义了形状函数的插值点；因为除了一个形状函数外，所有的形状函数在这些位置都是零，所以我们得到一个对角线质量矩阵，具有

@f{align*}
  B_{ii} = \int_\Omega \varphi_i(\mathbf x)^2\ \textrm{d}x,
  \qquad
  B_{ij}=0 \ \text{for } i\neq j.


@f}

为了定义 $G$ ，我们使用与 $B$ 相同的技术。换句话说，我们定义

@f{align*}
  G_{i} = \int_\Omega g_h(x) \varphi_i(\mathbf x)\ \textrm{d}x,


@f}

其中 $g_h$ 是 $g$ 的一个合适的近似值。然后， $B_{ii}$ 和 $G_i$ 定义中的积分由梯形规则近似。有了这个，上面的方程可以重述为

@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &U_i-B_{ii}^{-1}G_i \ge 0, \quad \Lambda_i \leq 0,\quad \Lambda_i[U_i-B_{ii}^{-1}G_i] = 0
\qquad \forall i\in{\cal S}.&


@f}



现在我们为每个自由度 $i$ 定义函数

@f{equation*}
 C([BU]_i,\Lambda_i) \dealcoloneq -\Lambda_i + \min\lbrace 0, \Lambda_i + c([BU]_i - G_i) \rbrace,


@f}

在这个程序中，我们选择 $c>0$ 。这是一种惩罚参数，取决于问题本身，需要选择足够大的参数；例如，如果我们使用7个全局细化，使用当前程序对 $c = 1$ 没有收敛作用）。)

经过一番挠头，人们可以说服自己，上面的不等式可以等效地改写为

@f{equation*}
 C([BU]_i,\Lambda_i) = 0, \qquad \forall i\in{\cal S}.


@f}

我们在这里将使用的原始-双重主动集策略是一个迭代方案，它基于这个条件来预测下一个主动集和非主动集 $\mathcal{A}_k$ 和 $\mathcal{F}_k$ （即那些指数 $i$ 的互补集，对于这些指数 $U_i$ 要么等于要么不等于障碍物的值 $B^{-1}G$  ）。关于这种方法的更深入的处理，见Hintermueller, Ito, Kunisch:The primal-dual active set strategy as a semismooth newton method, SIAM J. OPTIM., 2003, Vol.13, No.3, pp.865-888.

<a name="Theprimaldualactivesetalgorithm"></a><h3>The primal-dual active set algorithm</h3>


初级-二级主动集方法的算法工作原理如下（注： $B = B^T$  ）。

1.初始化  $\mathcal{A}_k$  和  $\mathcal{F}_k$  ，使  $\mathcal{S}=\mathcal{A}_k\cup\mathcal{F}_k$  和  $\mathcal{A}_k\cap\mathcal{F}_k=\emptyset$  并设置  $k=1$  。2.找出满足@f{align*}
  AU^k + B\Lambda^k &= F,\\
  [BU^k]_i &= G_i\quad&&\forall i\in\mathcal{A}_k,\\
  \Lambda_i^k &= 0\quad&&\forall i\in\mathcal{F}_k.
 @f}的原始-双数对 $(U^k,\Lambda^k)$ 。

请注意，第二个和第三个条件意味着正好 $|S|$ 个未知数是固定的，第一个条件产生了确定 $U$ 和 $\Lambda$ 所需的剩余 $|S|$ 个方程。3.3. 用@f{equation*}
 \begin{split}
  \mathcal{A}_{k+1} \dealcoloneq \lbrace i\in\mathcal{S}:\Lambda^k_i + c([BU^k]_i - G_i)< 0\rbrace,\\
  \mathcal{F}_{k+1} \dealcoloneq \lbrace i\in\mathcal{S}:\Lambda^k_i + c([BU^k]_i - G_i)\geq 0\rbrace.
 \end{split}
 @f}定义新的活动和非活动集。

如果 $\mathcal{A}_{k+1}=\mathcal{A}_k$ （然后，显然也是 $\mathcal{F}_{k+1}=\mathcal{F}_k$ ），则停止，否则设置 $k=k+1$ 并转到步骤（2）。

该方法被称为 "原始-双重"，因为它同时使用原始变量（位移 $U$ ）以及双重变量（拉格朗日乘数 $\Lambda$ ）来确定下一个活动集。

在本节的最后，让我们补充两点意见。首先，对于任何满足这些条件的原始-双重对 $(U^k,\Lambda^k)$ ，我们可以区分以下几种情况。

1.   $\Lambda^k_i + c([BU^k]_i - G_i) < 0$  (i active)。     <br>  然后是 $[BU^k]_i<G_i$ 和 $\Lambda^k_i=0$ （渗透）或 $\Lambda^k_i<0$ 和 $[BU^k]_i=G_i$ （压载）。2.   $\Lambda^k_i + c([BU^k]_i - G_i)\geq 0$  (i不活动)。     <br>  然后是 $[BU^k]_i\geq G_i$ 和 $\Lambda^k_i=0$ （无接触）或 $\Lambda^k_i\geq0$ 和 $[BU^k]_i=G_i$ （无压迫负荷）。

第二，上面的方法在直觉上似乎是正确的，也是有用的，但有点临时性的。然而，它可以通过以下方式简明地推导出来。为此，请注意，我们要解决的是非线性系统

@f{eqnarray*}
 &A U + B\Lambda = F,&\\
 &C([BU-G]_i, \Lambda_i) = 0,
\qquad \forall i.&


@f}

我们可以通过始终围绕前一个迭代进行线性化（即应用牛顿方法）来迭代解决，但为此我们需要对不可微分的函数 $C(\cdot,\cdot)$ 进行线性化。也就是说，它是可微的，事实上我们有

@f{equation*}
 \dfrac{\partial}{\partial U^k_i}C([BU^k]_i,\Lambda^k_i) = \begin{cases}
                                   cB_{ii},& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)< 0\\
                                   0,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)\geq 0.
                                  \end{cases}


@f}



@f{equation*}
 \dfrac{\partial}{\partial\Lambda^k_i}C([BU^k]_i,\Lambda^k_i) = \begin{cases}
                                   0,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)< 0\\


                                   -1,& \text{if}\ \Lambda^k_i + c([BU^k]_i - G_i)\geq 0.
                                  \end{cases}


@f}

这表明一个半光滑的牛顿步骤，其形式为

@f{equation*}
 \begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & B_{\mathcal{F}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & 0 & B_{\mathcal{A}_k}\\
 0 & 0 & -Id_{\mathcal{F}_k} & 0\\
 0 & cB_{\mathcal{A}_k} & 0 & 0
\end{pmatrix}
\begin{pmatrix}
 \delta U^k_{\mathcal{F}_k}\\ \delta U^k_{\mathcal{A}_k}\\ \delta \Lambda^k_{\mathcal{F}_k}\\ \delta \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=


-\begin{pmatrix}
 (AU^k + \Lambda^k - F)_{\mathcal{F}_k}\\ (AU^k + \Lambda^k - F)_{\mathcal{A}_k}\\ -\Lambda^k_{\mathcal{F}_k}\\ c(B_{\mathcal{A}_k} U^k - G)_{\mathcal{A}_k}
\end{pmatrix},


@f}

其中，我们将矩阵  $A,B$  以及向量以自然的方式分成行和列，其索引属于活动集  ${\mathcal{A}_k}$  或非活动集  ${\mathcal{F}_k}$  。

我们也可以通过设置 $\delta U^k \dealcoloneq
U^{k+1} - U^k$ 和 $\delta \Lambda^k \dealcoloneq \Lambda^{k+1} - \Lambda^k$ 并将所有已知项带到右手边来解决我们感兴趣的变量，而不是求解更新 $\delta U, \delta \Lambda$  。这就得到了

@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & B_{\mathcal{F}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & 0 & B_{\mathcal{A}_k}\\
 0 & 0 & Id_{\mathcal{F}_k} & 0\\
 0 & B_{\mathcal{A}_k} & 0 & 0
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}\\ \Lambda^k_{\mathcal{F}_k}\\ \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}\\ 0\\ G_{\mathcal{A}_k}
\end{pmatrix}.


@f}

这些是上文描述基本算法时概述的方程式。

我们甚至可以进一步推动这一点。很容易看出，我们可以消除第三行和第三列，因为它意味着 $\Lambda_{\mathcal{F}_k} = 0$  。

@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} & 0\\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k} & B_{\mathcal{A}_k}\\
 0 & B_{\mathcal{A}_k} & 0
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}\\ \Lambda^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}\\ G_{\mathcal{A}_k}
\end{pmatrix}.


@f}

这表明，事实上我们只需要解决位于活动集上的拉格朗日乘数。通过考虑第二行，我们将通过以下方式恢复全部拉格朗日乘数向量

@f{equation*}
 \Lambda^k_S = B^{-1}\left(f_{\mathcal{S}} - A_{\mathcal{S}}U^k_{\mathcal{S}}\right).


@f}

由于第三行和 $B_{\mathcal{A}_k}$ 是一个对角线矩阵的事实，我们能够直接计算出 $U^k_{\mathcal{A}_k}=B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}$ 。因此，我们也可以把线性系统写成如下。

@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & 0\\
 0 & Id_{\mathcal{A}_k} \\
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k} - A_{\mathcal{F}_k\mathcal{A}_k}B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}
 \\
 B_{\mathcal{A}_k}^{-1}G_{\mathcal{A}_k}
\end{pmatrix}.


@f}

幸运的是，这种形式很容易得出：我们只需建立通常的拉普拉斯线性系统即可

@f{equation*}
\begin{pmatrix}
 A_{\mathcal{F}_k\mathcal{F}_k} & A_{\mathcal{F}_k\mathcal{A}_k} \\
 A_{\mathcal{A}_k\mathcal{F}_k} & A_{\mathcal{A}_k\mathcal{A}_k}
\end{pmatrix}
\begin{pmatrix}
 U^k_{\mathcal{F}_k}\\ U^k_{\mathcal{A}_k}
\end{pmatrix}
=
\begin{pmatrix}
 F_{\mathcal{F}_k}\\ F_{\mathcal{A}_k}
\end{pmatrix},


@f}

然后让AffineConstraints类消除所有受限自由度，即 $U^k_{\mathcal{A}_k}=B^{-1}_{\mathcal{A}_k}G_{\mathcal{A}_k}$ ，其方式与 $\mathcal{A}_k$ 中的自由度是Dirichlet数据一样。结果线性系统（上面的第二个到最后一个）是对称和正定的，我们用CG方法和Trilinos的AMG预处理程序来解决它。




<a name="Implementation"></a><h3>Implementation</h3>


本教程与第4步很相似。程序的总体结构遵循步骤4，但略有不同。

- 我们需要两个新的方法，  <code>assemble_mass_matrix_diagonal</code>  和  <code>update_solution_and_constraints</code>  。

- 我们需要新的成员变量来表示我们这里的约束。

- 我们改变求解器的预处理程序。


如果你想了解目前的计划，你可能想阅读一下步骤4。


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
 * 像往常一样，在开始的时候，我们把所有我们需要的头文件都包含在这里。除了为Trilinos库提供接口的各种文件外，没有什么意外。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/index_set.h> 
 * 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/trilinos_sparse_matrix.h> 
 * #include <deal.II/lac/trilinos_vector.h> 
 * #include <deal.II/lac/trilinos_precondition.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step41 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeObstacleProblemcodeclasstemplate"></a> 
 * <h3>The <code>ObstacleProblem</code> class template</h3>
 * 

 * 
 * 该类提供了描述障碍问题所需的所有函数和变量。它与我们在 step-4 中要做的事情很接近，所以相对简单。唯一真正的新组件是计算主动集合的update_solution_and_constraints函数和一些描述线性系统原始（无约束）形式所需的变量（ <code>complete_system_matrix</code> 和 <code>complete_system_rhs</code> ），以及主动集合本身和主动集合公式中用于缩放拉格朗日乘数的质量矩阵 $B$ 的对角线。其余的内容与 step-4 相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ObstacleProblem 
 *   { 
 *   public: 
 *     ObstacleProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void make_grid(); 
 *     void setup_system(); 
 *     void assemble_system(); 
 *     void 
 *          assemble_mass_matrix_diagonal(TrilinosWrappers::SparseMatrix &mass_matrix); 
 *     void update_solution_and_constraints(); 
 *     void solve(); 
 *     void output_results(const unsigned int iteration) const; 
 * 
 *     Triangulation<dim>        triangulation; 
 *     FE_Q<dim>                 fe; 
 *     DoFHandler<dim>           dof_handler; 
 *     AffineConstraints<double> constraints; 
 *     IndexSet                  active_set; 
 * 
 *     TrilinosWrappers::SparseMatrix system_matrix; 
 *     TrilinosWrappers::SparseMatrix complete_system_matrix; 
 * 
 *     TrilinosWrappers::MPI::Vector solution; 
 *     TrilinosWrappers::MPI::Vector system_rhs; 
 *     TrilinosWrappers::MPI::Vector complete_system_rhs; 
 *     TrilinosWrappers::MPI::Vector diagonal_of_mass_matrix; 
 *     TrilinosWrappers::MPI::Vector contact_force; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Righthandsideboundaryvaluesandtheobstacle"></a> 
 * <h3>Right hand side, boundary values, and the obstacle</h3>
 * 

 * 
 * 在下文中，我们定义了描述右侧函数、Dirichlet边界值以及作为 $\mathbf x$ 函数的障碍物高度的类。在这三种情况下，我们都从函数 @<dim@>, 派生出这些类，尽管在 <code>RightHandSide</code> 和 <code>Obstacle</code> 的情况下，这更多的是出于惯例而非必要，因为我们从未将此类对象传递给库。在任何情况下，鉴于我们选择了 $f=-10$  ,  $u|_{\partial\Omega}=0$  ...，右手和边界值类的定义是显而易见的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & /*p*/, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       AssertIndexRange(component, 1); 
 * 
 *       return -10; 
 *     } 
 *   }; 
 * 
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & /*p*/, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       AssertIndexRange(component, 1); 
 * 
 *       return 0; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 我们用一个级联的障碍物来描述障碍物的功能（想想看：楼梯的阶梯）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Obstacle : public Function<dim> 
 *   { 
 *   public: 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       (void)component; 
 *       Assert(component == 0, ExcIndexRange(component, 0, 1)); 
 * 
 *       if (p(0) < -0.5) 
 *         return -0.2; 
 *       else if (p(0) >= -0.5 && p(0) < 0.0) 
 *         return -0.4; 
 *       else if (p(0) >= 0.0 && p(0) < 0.5) 
 *         return -0.6; 
 *       else 
 *         return -0.8; 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeObstacleProblemcodeclass"></a> 
 * <h3>Implementation of the <code>ObstacleProblem</code> class</h3>
 * 
 * <a name="ObstacleProblemObstacleProblem"></a> 
 * <h4>ObstacleProblem::ObstacleProblem</h4>
 * 

 * 
 * 对每个看过前几个教程程序的人来说，构造函数是完全显而易见的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   ObstacleProblem<dim>::ObstacleProblem() 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemmake_grid"></a> 
 * <h4>ObstacleProblem::make_grid</h4>
 * 

 * 
 * 我们在二维的正方形 $[-1,1]\times [-1,1]$ 上解决我们的障碍物问题。因此这个函数只是设置了一个最简单的网格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::make_grid() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -1, 1); 
 *     triangulation.refine_global(7); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "Total number of cells: " << triangulation.n_cells() 
 *               << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemsetup_system"></a> 
 * <h4>ObstacleProblem::setup_system</h4>
 * 

 * 
 * 在这个值得注意的第一个函数中，我们设置了自由度处理程序，调整了向量和矩阵的大小，并处理了约束。最初，约束条件当然只是由边界值给出的，所以我们在函数的顶部对它们进行插值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::setup_system() 
 *   { 
 *     dof_handler.distribute_dofs(fe); 
 *     active_set.set_size(dof_handler.n_dofs()); 
 * 
 *     std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl 
 *               << std::endl; 
 * 
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              BoundaryValues<dim>(), 
 *                                              constraints); 
 *     constraints.close(); 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false); 
 * 
 *     system_matrix.reinit(dsp); 
 *     complete_system_matrix.reinit(dsp); 
 * 
 *     IndexSet solution_index_set = dof_handler.locally_owned_dofs(); 
 *     solution.reinit(solution_index_set, MPI_COMM_WORLD); 
 *     system_rhs.reinit(solution_index_set, MPI_COMM_WORLD); 
 *     complete_system_rhs.reinit(solution_index_set, MPI_COMM_WORLD); 
 *     contact_force.reinit(solution_index_set, MPI_COMM_WORLD); 
 * 
 * @endcode
 * 
 * 这里唯一要做的事情是计算 $B$ 矩阵中的因子，该矩阵用于缩放残差。正如在介绍中所讨论的，我们将使用一个小技巧来使这个质量矩阵成为对角线，在下文中，首先将所有这些计算成一个矩阵，然后提取对角线元素供以后使用。
 * 

 * 
 * 
 * @code
 *     TrilinosWrappers::SparseMatrix mass_matrix; 
 *     mass_matrix.reinit(dsp); 
 *     assemble_mass_matrix_diagonal(mass_matrix); 
 *     diagonal_of_mass_matrix.reinit(solution_index_set); 
 *     for (unsigned int j = 0; j < solution.size(); j++) 
 *       diagonal_of_mass_matrix(j) = mass_matrix.diag_element(j); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemassemble_system"></a> 
 * <h4>ObstacleProblem::assemble_system</h4>
 * 

 * 
 * 这个函数一次就把系统矩阵和右手边集合起来，并把约束条件（由于活动集以及来自边界值）应用到我们的系统中。否则，它在功能上等同于例如  step-4  中的相应函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::assemble_system() 
 *   { 
 *     std::cout << "   Assembling system..." << std::endl; 
 * 
 *     system_matrix = 0; 
 *     system_rhs    = 0; 
 * 
 *     const QGauss<dim>  quadrature_formula(fe.degree + 1); 
 *     RightHandSide<dim> right_hand_side; 
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
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         cell_matrix = 0; 
 *         cell_rhs    = 0; 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             { 
 *               for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                 cell_matrix(i, j) += 
 *                   (fe_values.shape_grad(i, q_point) * 
 *                    fe_values.shape_grad(j, q_point) * fe_values.JxW(q_point)); 
 * 
 *               cell_rhs(i) += 
 *                 (fe_values.shape_value(i, q_point) * 
 *                  right_hand_side.value(fe_values.quadrature_point(q_point)) * 
 *                  fe_values.JxW(q_point)); 
 *             } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         constraints.distribute_local_to_global(cell_matrix, 
 *                                                cell_rhs, 
 *                                                local_dof_indices, 
 *                                                system_matrix, 
 *                                                system_rhs, 
 *                                                true); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemassemble_mass_matrix_diagonal"></a> 
 * <h4>ObstacleProblem::assemble_mass_matrix_diagonal</h4>
 * 

 * 
 * 下一个函数用于计算对角线质量矩阵 $B$ ，用于在主动集方法中缩放变量。正如介绍中所讨论的，我们通过选择正交的梯形规则来获得质量矩阵的对角线。这样一来，我们就不再需要在正交点、指数 $i$ 和指数 $j$ 上进行三重循环，而是可以直接使用双重循环。考虑到我们在以前的许多教程程序中讨论过的内容，该函数的其余部分是显而易见的。
 * 

 * 
 * 注意在调用这个函数的时候，约束对象只包含边界值约束；因此我们在最后的复制-本地-全局步骤中不必注意保留矩阵项的值，这些项以后可能会受到活动集的约束。
 * 

 * 
 * 还需要注意的是，只有在我们拥有 $Q_1$ 元素的情况下，使用梯形规则的技巧才有效。对于更高阶的元素，我们需要使用一个正交公式，在有限元的所有支持点都有正交点。构建这样一个正交公式其实并不难，但不是这里的重点，所以我们只是在函数的顶部断言我们对有限元的隐含假设实际上得到了满足。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::assemble_mass_matrix_diagonal( 
 *     TrilinosWrappers::SparseMatrix &mass_matrix) 
 *   { 
 *     Assert(fe.degree == 1, ExcNotImplemented()); 
 * 
 *     const QTrapezoid<dim> quadrature_formula; 
 *     FEValues<dim>         fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *  
 *       { 
 *         fe_values.reinit(cell); 
 *         cell_matrix = 0; 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             cell_matrix(i, i) += 
 *               (fe_values.shape_value(i, q_point) * 
 *                fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         constraints.distribute_local_to_global(cell_matrix, 
 *                                                local_dof_indices, 
 *                                                mass_matrix); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemupdate_solution_and_constraints"></a> 
 * <h4>ObstacleProblem::update_solution_and_constraints</h4>
 * 

 * 
 * 在某种意义上，这是本程序的核心功能。 它更新了介绍中所讨论的受限自由度的活动集，并从中计算出一个AffineConstraints对象，然后可以用来在下一次迭代的解中消除受限自由度。同时，我们将解决方案的受限自由度设置为正确的值，即障碍物的高度。
 * 

 * 
 * 从根本上说，这个函数是相当简单的。我们必须在所有自由度上循环，并检查函数 $\Lambda^k_i + c([BU^k]_i - G_i) = \Lambda^k_i + cB_i(U^k_i - [g_h]_i)$ 的符号，因为在我们的例子中 $G_i = B_i[g_h]_i$  。为此，我们使用介绍中给出的公式，通过该公式我们可以计算出拉格朗日乘数，作为原始线性系统的残差（通过变量 <code>complete_system_matrix</code> and <code>complete_system_rhs</code> 给出。在这个函数的顶部，我们使用一个属于矩阵类的函数来计算这个残差。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::update_solution_and_constraints() 
 *   { 
 *     std::cout << "   Updating active set..." << std::endl; 
 * 
 *     const double penalty_parameter = 100.0; 
 * 
 *     TrilinosWrappers::MPI::Vector lambda( 
 *       complete_index_set(dof_handler.n_dofs())); 
 *     complete_system_matrix.residual(lambda, solution, complete_system_rhs); 
 * 
 * @endcode
 * 
 * 计算 contact_force[i] =
 * 

 * 
 * - lambda[i] * diagonal_of_mass_matrix[i]。
 * 

 * 
 * 
 * @code
 *     contact_force = lambda; 
 *     contact_force.scale(diagonal_of_mass_matrix); 
 *     contact_force *= -1; 
 * 
 * @endcode
 * 
 * 下一步是重置活动集和约束对象，并在所有自由度上开始循环。由于我们不能只是在解向量的所有元素上循环，所以这变得稍微复杂了一些，因为我们没有办法找出一个自由度与哪个位置相关；但是，我们需要这个位置来测试一个自由度的位移是大于还是小于这个位置的障碍物高度。
 * 

 * 
 * 我们通过在所有单元和定义在每个单元上的DoF上循环来解决这个问题。我们在这里使用一个 $Q_1$ 函数来描述位移，对于该函数，自由度总是位于单元格的顶点上；因此，我们可以通过询问顶点来获得每个自由度的索引及其位置。另一方面，这显然对高阶元素不起作用，因此我们添加了一个断言，确保我们只处理所有自由度都位于顶点的元素，以避免万一有人想玩增加解的多项式程度时用非功能性代码绊倒自己。
 * 

 * 
 * 循环单元格而不是自由度的代价是我们可能会多次遇到一些自由度，即每次我们访问与给定顶点相邻的一个单元格时。因此，我们必须跟踪我们已经接触过的顶点和尚未接触的顶点。我们通过使用一个标志数组来做到这一点  <code>dof_touched</code>  。
 * 

 * 
 * 
 * @code
 *     constraints.clear(); 
 *     active_set.clear(); 
 * 
 *     const Obstacle<dim> obstacle; 
 *     std::vector<bool>   dof_touched(dof_handler.n_dofs(), false); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       for (const auto v : cell->vertex_indices()) 
 *         { 
 *           Assert(dof_handler.get_fe().n_dofs_per_cell() == cell->n_vertices(), 
 *                  ExcNotImplemented()); 
 * 
 *           const unsigned int dof_index = cell->vertex_dof_index(v, 0); 
 * 
 *           if (dof_touched[dof_index] == false) 
 *             dof_touched[dof_index] = true; 
 *           else 
 *             continue; 
 * 
 * @endcode
 * 
 * 现在我们知道我们还没有触及这个DoF，让我们得到那里的位移函数的值以及障碍函数的值，并使用这个来决定当前DoF是否属于活动集。为此，我们使用上面和介绍中给出的函数。
 * 

 * 
 * 如果我们决定该DoF应该是活动集的一部分，我们将其索引添加到活动集中，在AffineConstraints对象中引入一个不均匀的平等约束，并将解的值重置为障碍物的高度。最后，系统的非接触部分的残差作为一个额外的控制（残差等于剩余的、未计算的力，在接触区之外应该为零），所以我们把残差向量的分量（即拉格朗日乘数lambda）清零，这些分量对应于身体接触的区域；在所有单元的循环结束时，残差将因此只包括非接触区的残差。我们在循环结束后输出这个残差的准则和活动集的大小。
 * 

 * 
 * 
 * @code
 *           const double obstacle_value = obstacle.value(cell->vertex(v)); 
 *           const double solution_value = solution(dof_index); 
 * 
 *           if (lambda(dof_index) + penalty_parameter * 
 *                                     diagonal_of_mass_matrix(dof_index) * 
 *                                     (solution_value - obstacle_value) < 
 *               0) 
 *             { 
 *               active_set.add_index(dof_index); 
 *               constraints.add_line(dof_index); 
 *               constraints.set_inhomogeneity(dof_index, obstacle_value); 
 * 
 *               solution(dof_index) = obstacle_value; 
 * 
 *               lambda(dof_index) = 0; 
 *             } 
 *         } 
 *     std::cout << "      Size of active set: " << active_set.n_elements() 
 *               << std::endl; 
 * 
 *     std::cout << "   Residual of the non-contact part of the system: " 
 *               << lambda.l2_norm() << std::endl; 
 * 
 * @endcode
 * 
 * 在最后一步中，我们将迄今为止从活动集合中得到的对DoF的约束加入到那些由Dirichlet边界值产生的约束中，并关闭约束对象。
 * 

 * 
 * 
 * @code
 *     VectorTools::interpolate_boundary_values(dof_handler, 
 *                                              0, 
 *                                              BoundaryValues<dim>(), 
 *                                              constraints); 
 *     constraints.close(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemsolve"></a> 
 * <h4>ObstacleProblem::solve</h4>
 * 

 * 
 * 关于求解函数，其实没有什么可说的。在牛顿方法的背景下，我们通常对非常高的精度不感兴趣（为什么要求一个高度精确的线性问题的解，而我们知道它只能给我们一个非线性问题的近似解），所以我们使用ReductionControl类，当达到一个绝对公差（为此我们选择 $10^{-12}$ ）或者当残差减少一定的系数（这里是 $10^{-3}$ ）时停止反复运算。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::solve() 
 *   { 
 *     std::cout << "   Solving system..." << std::endl; 
 * 
 *     ReductionControl                        reduction_control(100, 1e-12, 1e-3); 
 *     SolverCG<TrilinosWrappers::MPI::Vector> solver(reduction_control); 
 *     TrilinosWrappers::PreconditionAMG       precondition; 
 *     precondition.initialize(system_matrix); 
 * 
 *     solver.solve(system_matrix, solution, system_rhs, precondition); 
 *     constraints.distribute(solution); 
 * 
 *     std::cout << "      Error: " << reduction_control.initial_value() << " -> " 
 *               << reduction_control.last_value() << " in " 
 *               << reduction_control.last_step() << " CG iterations." 
 *               << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemoutput_results"></a> 
 * <h4>ObstacleProblem::output_results</h4>
 * 

 * 
 * 我们使用vtk-format进行输出。 该文件包含位移和活动集的数字表示。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::output_results(const unsigned int iteration) const 
 *   { 
 *     std::cout << "   Writing graphical output..." << std::endl; 
 * 
 *     TrilinosWrappers::MPI::Vector active_set_vector( 
 *       dof_handler.locally_owned_dofs(), MPI_COMM_WORLD); 
 *     for (const auto index : active_set) 
 *       active_set_vector[index] = 1.; 
 * 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "displacement"); 
 *     data_out.add_data_vector(active_set_vector, "active_set"); 
 *     data_out.add_data_vector(contact_force, "lambda"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output_vtk("output_" + 
 *                              Utilities::int_to_string(iteration, 3) + ".vtk"); 
 *     data_out.write_vtk(output_vtk); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ObstacleProblemrun"></a> 
 * <h4>ObstacleProblem::run</h4>
 * 

 * 
 * 这是一个对所有事情都有最高级别控制的函数。 它并不长，而且事实上相当直接：在主动集方法的每一次迭代中，我们都要组装线性系统，求解它，更新主动集并将解投射回可行集，然后输出结果。只要主动集在前一次迭代中没有变化，迭代就会终止。
 * 

 * 
 * 唯一比较棘手的部分是，我们必须在第一次迭代组装好线性系统（即矩阵和右手边）后保存它。原因是这是唯一一个我们可以在没有任何接触约束的情况下访问线性系统的步骤。我们需要这个来计算其他迭代中解决方案的残差，但是在其他迭代中，我们形成的线性系统中对应于约束自由度的行和列都被消除了，因此我们不能再访问原始方程的全部残差。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void ObstacleProblem<dim>::run() 
 *   { 
 *     make_grid(); 
 *     setup_system(); 
 * 
 *     IndexSet active_set_old(active_set); 
 *     for (unsigned int iteration = 0; iteration <= solution.size(); ++iteration) 
 *       { 
 *         std::cout << "Newton iteration " << iteration << std::endl; 
 * 
 *         assemble_system(); 
 * 
 *         if (iteration == 0) 
 *           { 
 *             complete_system_matrix.copy_from(system_matrix); 
 *             complete_system_rhs = system_rhs; 
 *           } 
 * 
 *         solve(); 
 *         update_solution_and_constraints(); 
 *         output_results(iteration); 
 * 
 *         if (active_set == active_set_old) 
 *           break; 
 * 
 *         active_set_old = active_set; 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step41 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 这就是主函数。它遵循所有其他主函数的模式。调用初始化MPI是因为我们在这个程序中建立线性求解器的Trilinos库需要它。
 * 

 * 
 * 
 * @code
 * int main(int argc, char *argv[]) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step41; 
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
 *                     "This program can only be run in serial, use ./step-41")); 
 * 
 *       ObstacleProblem<2> obstacle_problem; 
 *       obstacle_problem.run(); 
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
examples/step-41/doc/results.dox



<a name="Results"></a><h1>Results</h1>


运行该程序会产生这样的输出。

@code
Number of active cells: 16384
Total number of cells: 21845
Number of degrees of freedom: 16641


Newton iteration 0
   Assembling system...
   Solving system...
      Error: 0.310059 -> 5.16619e-05 in 5 CG iterations.
   Updating active set...
      Size of active set: 13164
   Residual of the non-contact part of the system: 1.61863e-05
   Writing graphical output...


Newton iteration 1
   Assembling system...
   Solving system...
      Error: 1.11987 -> 0.00109377 in 6 CG iterations.
   Updating active set...
      Size of active set: 12363
   Residual of the non-contact part of the system: 3.9373
   Writing graphical output...


...


Newton iteration 17
   Assembling system...
   Solving system...
      Error: 0.00713308 -> 2.29249e-06 in 4 CG iterations.
   Updating active set...
      Size of active set: 5399
   Residual of the non-contact part of the system: 0.000957525
   Writing graphical output...


Newton iteration 18
   Assembling system...
   Solving system...
      Error: 0.000957525 -> 2.8033e-07 in 4 CG iterations.
   Updating active set...
      Size of active set: 5399
   Residual of the non-contact part of the system: 2.8033e-07
   Writing graphical output...
@endcode



一旦活动集不再变化，迭代就会结束（此时它有5399个受限自由度）。代数前提条件显然工作得很好，因为我们只需要4-6次CG迭代来解决线性系统（尽管这也与我们对线性求解器的精度要求不高有很大关系）。

更具启示性的是看一连串的图形输出文件（每三步显示一次，最左边一栏是迭代的编号）。

 <table align="center">
  <tr>
    <td valign="top">
      0 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.00.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.00.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.00.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      3 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.03.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.03.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.03.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      6 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.06.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.06.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.06.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      9 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.09.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      12 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.12.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.12.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.12.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      15 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.15.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.15.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.15.png" alt="">
    </td>
  </tr>
  <tr>
    <td valign="top">
      18 &nbsp;
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.18.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.active-set.18.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.displacement.3d.18.png" alt="">
    </td>
  </tr>
</table> 

图片显示，在第一步中，解决方案（在没有任何约束条件的情况下被计算出来的）是如此的弯曲，以至于几乎每一个内部点都必须被反弹到阶梯函数上，产生一个不连续的解决方案。在活动集迭代的过程中，这种不切实际的膜的形状被平滑掉了，与最下层阶梯的接触消失了，解决方案也稳定下来。

除此以外，程序还输出拉格朗日乘数的值。请记住，这些是接触力，所以在接触集上只应该是正的，而在接触集之外是零。另一方面，如果一个拉格朗日乘数在活动集上是负的，那么这个自由度必须从活动集上删除。下面的图片显示了迭代1、9和18中的乘数，我们用红色和棕色表示正值，蓝色表示负值。

 <table align="center">
  <tr>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.01.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.09.png" alt="">
    </td>
    <td valign="top">
      <img src="https://www.dealii.org/images/steps/developer/step-41.forces.18.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      Iteration 1
    </td>
    <td align="center">
      Iteration 9
    </td>
    <td align="center">
      Iteration 18
    </td>
  </tr>
</table> 

很容易看出，正值在接触集的内部很好地收敛为适度的值，在台阶的边缘有很大的向上的力，正如人们所期望的那样（以支持那里的膜的大曲率）；在活动集的边缘，乘数最初是负的，导致集合缩小，直到在迭代18，不再有负的乘数，算法已经收敛了。




<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


与本教程的任何程序一样，有许多明显的扩展和实验的可能性。第一个很清楚：引入自适应性。接触问题是自适应网格的主要候选者，因为解决方案有沿着它不太规则的线（膜和障碍物之间建立接触的地方）和解决方案非常光滑的其他区域（或者，在目前的情况下，在它与障碍物接触的地方是恒定的）。在目前的程序中加入这一点应该不会造成太多困难，但要为此找到一个好的误差估计器并非易事。

一个更具挑战性的任务是扩展到3D。这里的问题不是简单地让一切都在三维中运行。相反，当一个三维物体变形并与一个障碍物接触时，障碍物并不像这里的情况那样在域内作为一个约束体的力量发挥作用。相反，接触力只作用于物体的边界。那么不等式就不在微分方程中，而实际上在（诺伊曼型）边界条件中，尽管这导致了一种类似的变分不等式。在数学上，这意味着拉格朗日乘数只存在于表面，当然，如果方便的话，它也可以通过零扩展到域中。在目前的程序中，人们不需要明确地形成和存储这个拉格朗日乘数。

对于三维案例来说，另一个有趣的问题是考虑有摩擦的接触问题。在几乎每个机械过程中，摩擦都有很大的影响。为了建模，我们必须考虑到接触面的切向应力。我们还必须注意到，摩擦给我们的问题增加了另一个非线性。

另一个不简单的修改是实现一个更复杂的构成法则，如非线性弹性或弹塑性材料行为。这里的困难在于如何处理通过非线性构成法产生的额外非线性。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-41.cc"
*/
