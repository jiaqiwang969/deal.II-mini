/**
@page step_22 The step-22 tutorial program
This tutorial depends on step-6, step-21.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Weakform">Weak form</a>
        <li><a href="#Boundaryconditions">Boundary conditions</a>
        <li><a href="#Discretization">Discretization</a>
        <li><a href="#Linearsolverandpreconditioningissues">Linear solver and preconditioning issues</a>
      <ul>
        <li><a href="#IsthishowoneshouldsolvetheStokesequations"> Is this how one should solve the Stokes equations? </a>
        <li><a href="#Anoteonthestructureofthelinearsystem"> A note on the structure of the linear system </a>
      </ul>
        <li><a href="#Thetestcase">The testcase</a>
        <li><a href="#Implementation">Implementation</a>
      <ul>
        <li><a href="#UsingimhomogeneousconstraintsforimplementingDirichletboundaryconditions">Using imhomogeneous constraints for implementing Dirichlet boundary conditions</a>
        <li><a href="#UsingAffineConstraintsforincreasingperformanceh4h4UsingAffineConstraintsforincreasingperformance">Using AffineConstraints for increasing performance</h4><h4>Using AffineConstraints for increasing performance</a>
        <li><a href="#Performanceoptimizations">Performance optimizations</a>
    </ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#Definingtheinnerpreconditionertype">Defining the inner preconditioner type</a>
        <li><a href="#ThecodeStokesProblemcodeclasstemplate">The <code>StokesProblem</code> class template</a>
        <li><a href="#Boundaryvaluesandrighthandside">Boundary values and right hand side</a>
        <li><a href="#Linearsolversandpreconditioners">Linear solvers and preconditioners</a>
      <ul>
        <li><a href="#ThecodeInverseMatrixcodeclasstemplate">The <code>InverseMatrix</code> class template</a>
        <li><a href="#ThecodeSchurComplementcodeclasstemplate">The <code>SchurComplement</code> class template</a>
      </ul>
        <li><a href="#StokesProblemclassimplementation">StokesProblem class implementation</a>
      <ul>
        <li><a href="#StokesProblemStokesProblem">StokesProblem::StokesProblem</a>
        <li><a href="#StokesProblemsetup_dofs">StokesProblem::setup_dofs</a>
        <li><a href="#StokesProblemassemble_system">StokesProblem::assemble_system</a>
        <li><a href="#StokesProblemsolve">StokesProblem::solve</a>
        <li><a href="#StokesProblemoutput_results">StokesProblem::output_results</a>
        <li><a href="#StokesProblemrefine_mesh">StokesProblem::refine_mesh</a>
        <li><a href="#StokesProblemrun">StokesProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Outputoftheprogramandgraphicalvisualization">Output of the program and graphical visualization</a>
      <ul>
        <li><a href="#2Dcalculations">2D calculations</a>
        <li><a href="#3Dcalculations">3D calculations</a>
      </ul>
        <li><a href="#Sparsitypattern">Sparsity pattern</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Improvedlinearsolverin3Dh4h4Improvedlinearsolverin3D">Improved linear solver in 3D</h4><h4>Improved linear solver in 3D</a>
      <ul>
        <li><a href="#BetterILUdecompositionbysmartreorderingh5h5BetterILUdecompositionbysmartreordering">Better ILU decomposition by smart reordering</h5><h5>Better ILU decomposition by smart reordering</a>
        <li><a href="#BetterpreconditionerfortheinnerCGsolver">Better preconditioner for the inner CG solver</a>
        <li><a href="#BlockSchurcomplementpreconditioner">Block Schur complement preconditioner</a>
        <li><a href="#Combiningtheblockpreconditionerandmultigrid">Combining the block preconditioner and multigrid</a>
        <li><a href="#Noblockmatricesandvectors">No block matrices and vectors</a>
      </ul>
        <li><a href="#Moreinterestingtestcases">More interesting testcases</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-22/doc/intro.dox

 <br> 

<i>This program was contributed by Martin Kronbichler and Wolfgang
Bangerth.
<br>
This material is based upon work partly supported by the National
Science Foundation under Award No. EAR-0426271 and The California Institute of
Technology. Any opinions, findings, and conclusions or recommendations
expressed in this publication are those of the author and do not
necessarily reflect the views of the National Science Foundation or of The
California Institute of Technology.
</i>




<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


这个程序处理斯托克斯方程组，其非维度形式如下。

@f{eqnarray*}


  -2\; \textrm{div}\; \varepsilon(\textbf{u}) + \nabla p &=& \textbf{f},
  \\


  -\textrm{div}\; \textbf{u} &=& 0,


@f}

其中 $\textbf u$ 表示流体的速度， $p$ 是其压力， $\textbf f$ 是外力， $\varepsilon(\textbf{u})= \nabla^s{\textbf{u}}= \frac 12 \left[
(\nabla \textbf{u}) + (\nabla \textbf{u})^T\right]$ 是对称梯度的第2级张量；它的分量定义是 $\varepsilon(\textbf{u})_{ij}=\frac
12\left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\right)$  。

斯托克斯方程描述了缓慢移动的粘性流体的稳态运动，如蜂蜜、地幔中的岩石，或其他惯性不起作用的情况。如果流体的运动速度足够快，与粘性摩擦力相比，惯性力很重要，那么斯托克斯方程就不再有效；考虑到惯性效应，就会产生非线性纳维-斯托克斯方程。然而，在这个教程程序中，我们将重点关注更简单的斯托克斯系统。

请注意，在推导更一般的可压缩纳维-斯托克斯方程时，扩散被建模为应力张量的发散。

@f{eqnarray*}
  \tau = - \mu (2\varepsilon(\textbf{u}) - \frac{2}{3}\nabla \cdot \textbf{u} I),


@f}

其中 $\mu$ 是流体的粘度。有了 $\mu=1$ 的假设（假设粘度恒定，通过除以 $\mu$ 使方程非立体化），并假设不可压缩性（ $\textrm{div}\; \textbf{u}=0$ ），我们就可以得出上面的公式。

@f{eqnarray*}
  \textrm{div}\; \tau = -2\textrm{div}\;\varepsilon(\textbf{u}).


@f}

一个不同的公式使用拉普拉斯算子（ $-\triangle \textbf{u}$ ）而不是对称梯度。这里一个很大的区别是，速度的不同分量并不耦合。如果你假定解的额外规律性 $\textbf{u}$ （第二部分导数存在并且是连续的），这些公式是等同的。

@f{eqnarray*}
  \textrm{div}\; \tau
  = -2\textrm{div}\;\varepsilon(\textbf{u})
  = -\triangle \textbf{u} + \nabla \cdot (\nabla\textbf{u})^T
  = -\triangle \textbf{u}.


@f}

这是因为 $i$ 中的第1个条目 $\nabla \cdot (\nabla\textbf{u})^T$ 是由以下内容给出的。

@f{eqnarray*}
[\nabla \cdot (\nabla\textbf{u})^T]_i
= \sum_j \frac{\partial}{\partial x_j} [(\nabla\textbf{u})^T]_{i,j}
= \sum_j \frac{\partial}{\partial x_j} [(\nabla\textbf{u})]_{j,i}
= \sum_j \frac{\partial}{\partial x_j} \frac{\partial}{\partial x_i} \textbf{u}_j
= \sum_j \frac{\partial}{\partial x_i} \frac{\partial}{\partial x_j} \textbf{u}_j
= \frac{\partial}{\partial x_i} \textrm{div}\; \textbf{u}
= 0.


@f}

如果你不能假设上述的规律性，或者你的粘度不是一个常数，那么这个等价关系就不再成立。因此，我们决定在本教程中坚持使用物理上更准确的对称张量公式。


为了得到良好的解决，我们将不得不在方程中加入边界条件。一旦我们讨论方程的微弱形式，哪些边界条件在这里是容易做到的就会变得清晰。

这里所涉及的方程属于矢量值问题的范畴。这个主题的顶层概述可以在 @ref vector_valued 模块中找到。




<a name="Weakform"></a><h3>Weak form</h3>


通过将其写成矢量形式，可以得到方程的弱形式为

@f{eqnarray*}
  \begin{pmatrix}
    {-2\; \textrm{div}\; \varepsilon(\textbf{u}) + \nabla p}
    \\
    {-\textrm{div}\; \textbf{u}}
  \end{pmatrix}
  =
  \begin{pmatrix}
  {\textbf{f}}
  \\
  0
  \end{pmatrix},


@f}

从左边开始与一个矢量值测试函数 $\phi = \begin{pmatrix}\textbf{v} \\ q\end{pmatrix}$ 形成点积，并在域 $\Omega$ 上进行积分，得到以下一组方程。

@f{eqnarray*}
  (\mathrm v,


   -2\; \textrm{div}\; \varepsilon(\textbf{u}) + \nabla p)_{\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  (\textbf{v}, \textbf{f})_\Omega,


@f}

这对所有的测试函数都必须成立  $\phi = \begin{pmatrix}\textbf{v}
\\ q\end{pmatrix}$  。

一般来说，一个好的经验法则是，如果一个人<i>can</i>减少公式中任何变量的导数，那么他<i>should</i>实际上是用部分积分来做。这是由<a
href="https://en.wikipedia.org/wiki/Partial_differential_equation">partial
differential equations</a>的理论引起的，特别是强和<a href="https://en.wikipedia.org/wiki/Weak_solution">weak
solutions</a>之间的区别）。我们已经为拉普拉斯方程做了这个工作，在那里我们对第二导数进行了分项积分，以获得在测试和试验函数上都只有一个导数的弱表述。

在当前情况下，我们对第二项进行分项积分。

@f{eqnarray*}
  (\textbf{v}, -2\; \textrm{div}\; \varepsilon(\textbf{u}))_{\Omega}


  - (\textrm{div}\; \textbf{v}, p)_{\Omega}
  + (\textbf{n}\cdot\textbf{v}, p)_{\partial\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  (\textbf{v}, \textbf{f})_\Omega.


@f}

同样地，我们对第一项进行分项积分，得到

@f{eqnarray*}
  (\nabla \textbf{v}, 2\; \varepsilon(\textbf{u}))_{\Omega}


  -
  (\textbf{n} \otimes \textbf{v}, 2\; \varepsilon(\textbf{u}))_{\partial\Omega}


  - (\textrm{div}\; \textbf{v}, p)_{\Omega}
  + (\textbf{n}\cdot\textbf{v}, p)_{\partial\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  (\textbf{v}, \textbf{f})_\Omega,


@f}

其中，两个张量值的量之间的标量乘积在此定义为

@f{eqnarray*}
  (\nabla \textbf{v}, 2\; \varepsilon(\textbf{u}))_{\Omega}
  =
  2 \int_\Omega \sum_{i,j=1}^d \frac{\partial v_j}{\partial x_i}
  \varepsilon(\textbf{u})_{ij} \ dx.


@f}

利用这一点，我们现在已经将对我们的变量的要求降低到 $\mathbf u,\mathbf v$ 的一阶导数，而对 $p,q$ 完全没有导数。

因为像 $\nabla\textbf{v}$ 这样的一般张量和 $\varepsilon(\textbf{u})$ 这样的对称张量之间的标量积等于两者的对称形式之间的标量积，我们也可以把上面的双线性形式写成如下。

@f{eqnarray*}
  (\varepsilon(\textbf{v}), 2\; \varepsilon(\textbf{u}))_{\Omega}


  -
  (\textbf{n} \otimes \textbf{v}, 2\; \varepsilon(\textbf{u}))_{\partial\Omega}


  - (\textrm{div}\; \textbf{v}, p)_{\Omega}
  + (\textbf{n}\cdot\textbf{v}, p)_{\partial\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  (\textbf{v}, \textbf{f})_\Omega,


@f}

我们将在下一节处理边界条款，但从域条款中已经可以看出

@f{eqnarray*}
  (\varepsilon(\textbf{v}), 2\; \varepsilon(\textbf{u}))_{\Omega}


  - (\textrm{div}\; \textbf{v}, p)_{\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}


@f}

的双线性形式，斯托克斯方程产生一个对称的双线性形式，并因此产生一个对称的（如果是不确定的）系统矩阵。




<a name="Boundaryconditions"></a><h3>Boundary conditions</h3>


 @dealiiVideoLecture{21.5}  （  @dealiiVideoLectureSeeAlso{21.55,21.6,21.65})  。

刚刚导出的弱形式立即为我们提供了施加边界条件的不同可能性。<ol>  <li>  迪里希特速度边界条件。在一个部分 $\Gamma_D\subset\partial\Omega$ ，我们可以对速度 $\textbf u$ 施加迪里希特条件。

    @f{eqnarray*}
        \textbf u = \textbf g_D \qquad\qquad \textrm{on}\ \Gamma_D.
    @f}

    因为测试函数 $\textbf{v}$ 来自解变量的切线空间，我们有 $\textbf{v}=0$ 对 $\Gamma_D$ ，因此有@f{eqnarray*}


      -(\textbf{n} \otimes \mathrm
        v, 2\; \varepsilon(\textbf{u}))_{\Gamma_D}
      +
      (\textbf{n}\cdot\textbf{v}, p)_{\Gamma_D}
      = 0.
    @f}。

    换句话说，像往常一样，强加的边界值并没有出现在弱形式中。

    值得注意的是，如果我们在整个边界上施加迪里希特边界值，那么压力就只能确定到一个常数。这方面的算法实现将使用类似于步骤11中的工具。

 <li>  诺伊曼型或自然边界条件。在边界的其余部分 $\Gamma_N=\partial\Omega\backslash\Gamma_D$ ，让我们把边界条款重新写成如下。     @f{eqnarray*}


      -(\textbf{n} \otimes \mathrm
        v, 2\; \varepsilon(\textbf{u}))_{\Gamma_N}
      +
      (\textbf{n}\cdot\textbf{v}, p)_{\Gamma_N}
      &=&
      \sum_{i,j=1}^d


      -(n_i v_j, 2\; \varepsilon(\textbf{u})_{ij})_{\Gamma_N}
      +
      \sum_{i=1}^d
      (n_i v_i, p)_{\Gamma_N}
      \\
      &=&
      \sum_{i,j=1}^d


      -(n_i v_j, 2\; \varepsilon(\textbf{u})_{ij})_{\Gamma_N}
      +
      \sum_{i,j=1}^d
      (n_i v_j, p \delta_{ij})_{\Gamma_N}
      \\
      &=&
      \sum_{i,j=1}^d
      (n_i v_j,p \delta_{ij} - 2\; \varepsilon(\textbf{u})_{ij})_{\Gamma_N}
      \\
      &=&
      (\textbf{n} \otimes \textbf{v},
      p \textbf{I} - 2\; \varepsilon(\textbf{u}))_{\Gamma_N}.
      \\
      &=&
      (\textbf{v},
       \textbf{n}\cdot [p \textbf{I} - 2\; \varepsilon(\textbf{u})])_{\Gamma_N}.
    @f}

    换句话说，在边界的诺伊曼部分，我们可以规定总应力的数值。     @f{eqnarray*}
      \textbf{n}\cdot [p \textbf{I} - 2\; \varepsilon(\textbf{u})]
      =
      \textbf g_N \qquad\qquad \textrm{on}\ \Gamma_N.
    @f}

    如果边界被细分为Dirichlet和Neumann部分  $\Gamma_D,\Gamma_N$  ，这就导致了以下弱形式。     @f{eqnarray*}
      (\varepsilon(\textbf{v}), 2\; \varepsilon(\textbf{u}))_{\Omega}


      - (\textrm{div}\; \textbf{v}, p)_{\Omega}


      -
      (q,\textrm{div}\; \textbf{u})_{\Omega}
      =
      (\textbf{v}, \textbf{f})_\Omega


      -
      (\textbf{v}, \textbf g_N)_{\Gamma_N}.
    @f}




 <li>  罗宾式边界条件。罗宾式边界条件是迪里切特和诺伊曼边界条件的混合物。它们将读作@f{eqnarray*}
      \textbf{n}\cdot [p \textbf{I} - 2\; \varepsilon(\textbf{u})]
      =
      \textbf S \textbf u \qquad\qquad \textrm{on}\ \Gamma_R,
    @f}。

    有一个等级2的张量（矩阵）  $\textbf S$  。相关的弱形式是@f{eqnarray*}
      (\varepsilon(\textbf{v}), 2\; \varepsilon(\textbf{u}))_{\Omega}


      - (\textrm{div}\; \textbf{v}, p)_{\Omega}


      -
      (q,\textrm{div}\; \textbf{u})_{\Omega}
      +
      (\textbf S \textbf u, \textbf{v})_{\Gamma_R}
      =
      (\textbf{v}, \textbf{f})_\Omega.
    @f} 。



 <li>  部分边界条件。可以通过只对速度的某些分量强制执行Dirichlet和Neumann边界条件来结合它们。例如，施加人工边界条件的一种方法是要求流动垂直于边界，即切向分量 $\textbf u_{\textbf t}=(\textbf
    1-\textbf n\otimes\textbf n)\textbf u$ 为零，从而约束速度的 <code>dim</code> -1分量。剩下的分量可以通过要求法向应力的法向分量为零来约束，产生以下一组边界条件。     @f{eqnarray*}
      \textbf u_{\textbf t} &=& 0,
      \\
      \textbf n \cdot \left(\textbf{n}\cdot [p \textbf{I} - 2\;
      \varepsilon(\textbf{u})] \right)
      &=&
      0.
    @f}



    另一种情况是当人们希望流动是<i>parallel</i>而不是垂直于边界时（在deal.II中， VectorTools::compute_no_normal_flux_constraints 函数可以为你这样做）。这种情况经常发生在自由边界的问题上（例如，在河流或湖泊的表面，如果流动的垂直力不足以使表面实际变形），或者如果边界对流体没有施加明显的摩擦力（例如，在地幔和地核的界面上，两种流体因密度不同而相遇，但它们的粘度都很小，不会对彼此产生很大的切向应力）。     在公式中，这意味着@f{eqnarray*}
      \textbf{n}\cdot\textbf u &=& 0,
      \\
      (\textbf 1-\textbf n\otimes\textbf n)
      \left(\textbf{n}\cdot [p \textbf{I} - 2\;
      \varepsilon(\textbf{u})] \right)
      &=&
      0,
    @f}

    第一个条件（需要强加）固定了速度的一个分量，第二个条件（将在弱形式下强制执行）固定了其余的两个分量。   </ol> 

尽管有这么多的可能性，我们在这个教程程序中只使用迪里希特和（同质）诺伊曼边界条件。




<a name="Discretization"></a><h3>Discretization</h3>


如上所述，在 $\Gamma_D$ 和 $\Gamma_N$ 上有迪里希特和诺伊曼边界条件的方程的弱形式是这样的：找到 $\textbf u\in \textbf V_g = \{\varphi \in H^1(\Omega)^d: \varphi_{\Gamma_D}=\textbf
g_D\}, p\in Q=L^2(\Omega)$ ，以便

@f{eqnarray*}
  (\varepsilon(\textbf{v}), 2\; \varepsilon(\textbf{u}))_{\Omega}


  - (\textrm{div}\; \textbf{v}, p)_{\Omega}


  -
  (q,\textrm{div}\; \textbf{u})_{\Omega}
  =
  (\textbf{v}, \textbf{f})_\Omega


  -
  (\textbf{v}, \textbf g_N)_{\Gamma_N}


@f}

为所有测试函数  $\textbf{v}\in \textbf V_0 = \{\varphi \in H^1(\Omega)^d: \varphi_{\Gamma_D}=0\},q\in
Q$  。

这些方程代表一个对称的<a
href="https://en.wikipedia.org/wiki/Ladyzhenskaya%E2%80%93Babu%C5%A1ka%E2%80%93Brezzi_condition">saddle
point problem</a>。众所周知，那么只有当我们寻找解决方案的函数空间必须满足某些条件时，解决方案才会存在，这些条件通常被称为Babuska-Brezzi或Ladyzhenskaya-Babuska-Brezzi（LBB）条件。上面的连续函数空间满足这些条件。然而，当我们将方程离散化，用有限维空间中的有限元函数取代连续变量和检验函数 $\textbf V_{g,h}\subset \textbf V_g,
Q_h\subset Q$ 时，我们必须确保 $\textbf V_h,Q_h$ 也满足LBB条件。这与我们在第20步中要做的事情类似。

对于斯托克斯方程，有许多可能的选择来确保有限元空间与LBB条件兼容。一个简单而准确的选择是 $\textbf u_h\in Q_{p+1}^d,
p_h\in Q_p$ ，即对速度使用比压力高一阶的元素。

这就导致了以下的离散问题：找到 $\textbf u_h,p_h$ ，以便于

@f{eqnarray*}
  (\varepsilon(\textbf{v}_h), 2\; \varepsilon(\textbf u_h))_{\Omega}


  - (\textrm{div}\; \textbf{v}_h, p_h)_{\Omega}


  -
  (q_h,\textrm{div}\; \textbf{u}_h)_{\Omega}
  =
  (\textbf{v}_h, \textbf{f})_\Omega


  -
  (\textbf{v}_h, \textbf g_N)_{\Gamma_N}


@f}

为所有测试函数  $\textbf{v}_h, q_h$  。组装与此问题相关的线性系统遵循  @ref step_20  "步骤-20"、步骤-21中使用的相同路线，并在  @ref
vector_valued  模块中详细解释。




<a name="Linearsolverandpreconditioningissues"></a><h3>Linear solver and preconditioning issues</h3>


离散方程的微弱形式自然导致了以下速度场和压力场的节点值的线性系统。

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{c}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{c}
    F \\ G
  \end{array}\right),


@f}

与第20步和第21步一样，我们将通过形成舒尔补数来解决这个方程组，也就是说，我们将首先找到 $P$ 的解。

@f{eqnarray*}
  BA^{-1}B^T P &=& BA^{-1} F - G, \\


@f}

然后

@f{eqnarray*}
  AU &=& F - B^TP.


@f}

我们这样做的方式与我们在以前的这些教程程序中做的差不多，也就是说，我们再次使用相同的类 <code>SchurComplement</code> 和 <code>InverseMatrix</code> 。然而，有两个显著的区别。

<ol>  <li>  首先，在混合拉普拉斯方程中，我们必须处理如何对舒尔补数 $B^TM^{-1}B$ 进行预处理的问题，它在谱上等同于压力空间上的拉普拉斯算子（因为 $B$ 代表梯度算子， $B^T$ 代表其邻接算子 $-\textrm{div}$ ，而 $M$ 代表身份（直到材料参数 $K^{-1}$ 为止），因此 $B^TM^{-1}B$ 类似于 $-\textrm{div} \mathbf 1 \nabla = -\Delta$ ）。因此，对于小的网格尺寸来说，矩阵的条件很差，我们不得不为Schur补数提出一个精心设计的预处理方案。

 <li>  其次，每次我们与 $B^TM^{-1}B$ 相乘时，我们必须用质量矩阵 $M$ 来解决。然而，这并不特别困难，因为质量矩阵总是有很好的条件的，所以使用CG和一点点预处理就能简单地反转。  换句话说， </ol> 的内部求解器的预处理很简单，而 $B^TM^{-1}B$ 的外部求解器的预处理很复杂。

在这里，情况几乎完全相反。差异源于这样一个事实，即舒尔补码的核心矩阵不是来自身份算子，而是来自拉普拉斯算子的一个变体， $-\textrm{div} \nabla^s$ （其中 $\nabla^s$ 是对称梯度），作用于一个矢量场。在对这个问题的研究中，我们主要遵循D. Silvester和A. Wathen的论文。"稳定的斯托克斯系统的快速迭代解第二部分。使用一般块状先决条件"。(SIAM J. Numer.Anal., 31 (1994), pp. 1352-1367)，可在线查阅<a
href="http://siamdl.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=SJNAAM000031000005001352000001&idtype=cvips&gifs=Yes" target="_top">here</a>。主要来说，舒尔补码的核心矩阵的差异有两个后果。

<ol>  <li>  首先，它使外部预处理变得简单：Schur补数对应于压力空间上的算子 $-\textrm{div} (-\textrm{div} \nabla^s)^{-1}
\nabla$ ；忘记我们处理的是对称梯度而不是常规梯度的事实，Schur补数类似于 $-\textrm{div} (-\textrm{div} \nabla)^{-1} \nabla =


-\textrm{div} (-\Delta)^{-1} \nabla$ ，即使在数学上不完全简明，在光谱上也等同于身份算子（一个启发式的论证是将算子换算成 $-\textrm{div}(-\Delta)^{-1} \nabla = -\textrm{div}\nabla(-\Delta)^{-1} =


-\Delta(-\Delta)^{-1} = \mathbf 1$ ）。事实证明，用CG方法直接解决这个Schur补数并不容易：在没有预处理的情况下，Schur补数矩阵的条件数取决于最大和最小单元的大小比，而且仍然需要50-100次CG迭代。然而，有一个简单的解决办法：用压力空间上的质量矩阵进行预处理，我们就可以减少到5-15次CG迭代，几乎不受网格结构的影响（看看这个程序的<a href="#Results">results section</a>，可以看到CG迭代的数量确实不会随着我们细化网格而改变）。

因此，除了我们已经有的东西之外，我们需要的是压力变量上的质量矩阵，我们将把它存储在一个单独的对象中。




 <li>  虽然与第20步讨论的混合拉普拉斯情况相比，外部预调节器变得简单了，但内部求解器的问题却变得更加复杂。在混合拉普拉斯离散化中，舒尔补数的形式为  $B^TM^{-1}B$  。因此，每当我们与舒尔补码相乘时，我们必须解决一个线性系统 $M_uz=y$ ；然而，这并不太复杂，因为压力空间上的质量矩阵 $M_u$ 是有条件的。


另一方面，对于我们这里考虑的斯托克斯方程，舒尔补码是 $BA^{-1}B^T$ ，其中矩阵 $A$ 与拉普拉斯算子有关（事实上，它是对应于双线性形式 $(\nabla^s \varphi_i, \nabla^s\varphi_j)$ 的矩阵）。因此，用 $A$ 求解要复杂得多：矩阵的条件很差，我们知道我们需要很多迭代，除非我们有一个非常好的预处理程序。更糟糕的是，我们每次与舒尔补码相乘时都要用 $A$ 求解，使用上述的预处理程序需要5-15次。

因为我们必须多次用 $A$ 求解，所以多花一次时间为这个矩阵创建一个好的预处理程序是值得的。所以我们要做的是：如果在2d中，我们使用终极预处理程序，即矩阵的直接稀疏LU分解。这是用SparseDirectUMFPACK类实现的，它使用UMFPACK直接求解器来计算分解。要使用它，你必须建立支持UMFPACK的deal.II（这是默认的）；参见<a href="../../readme.html#optional-software">ReadMe file</a>中的说明。有了它，内解器在一次迭代中就能收敛。

在2D中，我们可以做这种事情，因为即使是合理的大问题，也很少有超过100,000个未知数的，每行的非零项相对较少。此外，2D中矩阵的带宽是 ${\cal
O}(\sqrt{N})$ ，因此是中等的。对于这样的矩阵，稀疏因子可以在几秒钟内计算出来。作为参考，计算一个大小为 $N$ 、带宽为 $B$ 的矩阵的稀疏因子需要 ${\cal
O}(NB^2)$ 次操作。在2d中，这是 ${\cal O}(N^2)$ ；尽管这比例如组装线性系统的复杂度要高，后者需要 ${\cal
O}(N)$ ，但计算分解的常数非常小，直到我们达到非常大的未知数%，甚至更多，它才成为整个程序中的主导因素）。)

情况在3D中发生了变化，因为在那里我们很快就会有更多的未知数，而且矩阵的带宽（决定了稀疏LU因子中非零项的数量）是 ${\cal O}(N^{2/3})$ ，而且每行也有很多条目。这使得像UMFPACK这样的稀疏直接求解器的效率很低：只有在问题规模为10,000到100,000个未知数时，才能用合理的时间和内存资源计算稀疏分解。

在这种情况下，我们要做的是使用一个不完整的LU分解（ILU）作为预处理，而不是实际计算完整的LU因子。恰好，deal.II有一个类可以做到这一点。SparseILU。计算ILU所需要的时间只取决于稀疏矩阵中的非零项的数量（或者说我们愿意填入LU因子，如果这些非零项应该多于矩阵中的非零项），但与矩阵的带宽无关。因此，这也是一个可以有效地在三维中计算的操作。另一方面，根据定义，一个不完整的LU分解并不代表矩阵的精确逆  $A$  。因此，与稀疏直接求解器的预处理不同，用ILU进行预处理仍然需要一次以上的迭代。因此，当与舒尔补码相乘时，内解器将花费更多时间：这是一个不可避免的权衡。   </ol> 

在下面的程序中，我们将利用SparseILU和SparseDirectUMFPACK类具有非常相似的接口，可以互换使用的事实。我们所需要的是一个开关类，根据维度的不同，提供一个类型，即上述两个类中的任何一个。我们就是这样做的。

@code
template <int dim>
struct InnerPreconditioner;


template <>
struct InnerPreconditioner<2>
{
  using type = SparseDirectUMFPACK;
};


template <>
struct InnerPreconditioner<3>
{
  using type = SparseILU<double>;
};
@endcode



从这里开始，我们可以引用<code>typename  InnerPreconditioner@<dim@>::%type</code> 这个类型，并自动得到正确的预处理程序类。由于这两个类的接口相似，我们将能够在所有地方使用相同的语法来互换使用它们。




<a name="IsthishowoneshouldsolvetheStokesequations"></a><h4> Is this how one should solve the Stokes equations? </h4>


上面的讨论显示了解决由斯托克斯方程产生的线性系统的*种方式，由于辅导程序是教学工具，这是有意义的。但这是解决这个方程组的**方式吗？

这个问题的答案是否定的。上面已经指出了这个方法的主要瓶颈，就是我们必须反复求解Schur补数内的 $A$ 线性系统，由于我们没有一个好的Schur补数的预处理程序，这些求解就不得不经常发生。一个更好的方法是使用块分解，这是基于Silvester和Wathen  @cite SW94 的观察，在 @cite elman2005 中有更详细的解释。下面在本程序的结果部分的a <a href="#block-schur">block Schur
complementation preconditioner</a>一节中讨论了这种替代方法的实现。




<a name="Anoteonthestructureofthelinearsystem"></a><h4> A note on the structure of the linear system </h4>


以上，我们已经声称线性系统的形式是

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),


@f}

即，特别是在矩阵的右下方有一个零块。这样，我们就可以把舒尔补码写成 $S=B A^{-1} B^T$  。但这并不完全正确。

想一想，如果对某些压力变量有约束会怎样（见 @ref constraints "自由度的约束 "文件模块），例如因为我们使用自适应细化网格和连续压力有限元，所以会有悬挂节点。造成这种约束的另一个原因是压力的迪里希特边界条件。然后，AffineConstraints类，在将矩阵的局部贡献复制到全局线性系统时，会将对应于受限自由度的行和列清零，并在对角线上放一个正条目。为了简单起见，你可以认为这个条目是1，尽管实际上它是一个与其他矩阵条目相同数量级的值）。换句话说，右下角区块其实根本不是空的：它在对角线上有几个条目，每个受限的压力自由度都有一个，对我们要解决的线性系统的正确描述是，它的形式如下

@f{eqnarray*}
  \left(\begin{array}{cc}
    A & B^T \\ B & D_c
  \end{array}\right)
  \left(\begin{array}{cc}
    U \\ P
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    F \\ G
  \end{array}\right),


@f}

其中 $D_c$ 是零矩阵，除了受限自由度的正对角线项。那么正确的舒尔补码实际上是矩阵 $S = B A^{-1} B^T - D_c $ ，而不是上面说的那个。

思考这个问题使我们，首先，意识到所得到的舒尔补数现在是不确定的，因为 $B A^{-1} B^T$ 是对称的和正定的，而 $D_c$ 是正半定的，从前者减去后者可能不再是正定的。这很烦人，因为我们不能再对这个真正的舒尔补数采用共轭梯度法。也就是说，我们可以在 AffineConstraints::distribute_local_to_global() 中通过简单地将*负*值放在受限压力变量的对角线上来解决这个问题--因为我们实际上只是放了一些非零的东西来确保结果矩阵不是奇异的；我们真的不关心那个条目是正还是负。因此，如果 $D_c$ 的对角线上的条目是负的，那么 $S$ 将再次成为一个对称的正定矩阵。

但是，其次，下面的代码实际上并没有做这些事。它很高兴地用错误的Schur补码 $S = B A^{-1} B^T$ 来解决线性系统，完全忽略了这个问题。为什么会这样呢？为了理解为什么会这样，回顾一下，当把局部贡献写入全局矩阵时， AffineConstraints::distribute_local_to_global() 把对应于受限自由度的行和列置零。这意味着 $B$ 有一些零行， $B^T$ 零列。因此，如果要乘出 $S$ 的条目是什么，就会发现它的所有受限压力自由度的行和列都是零，包括对角线上的一个零。 $D_c$ 的非零条目将正好适合于这些零对角线位置，并确保 $S$ 是可逆的。不这样做，严格来说，意味着 $S$ 仍然是单数。它在非约束压力自由度子集上是对称和正定的，而在约束压力上只是一个零矩阵。为什么共轭梯度法对这个矩阵有效？因为 AffineConstraints::distribute_local_to_global() 也确保了对应于矩阵这些零行的右手边条目*也是零，也就是说，右手边是兼容的。

这意味着无论这些受限压力自由度的解向量的值是多少，这些行的残差总是为零，如果考虑到CG算法的内部操作，就永远不会对解向量产生任何更新。换句话说，CG算法只是*忽略*这些行，尽管矩阵是单数。这只是因为这些自由度与线性系统的其他部分完全解耦（因为整个行和相应的列都是零）。在求解过程结束时，求解向量中的受限压力值仍然和我们开始调用求解器时一模一样；当我们在CG求解器完成后调用 AffineConstraints::distribute() 时，它们最终被正确的值所覆盖。

这个讨论的结果是，大矩阵的右下角块为零的假设有点简化，但仅仅按照这个假设，实际上并没有导致任何值得解决的实际问题。




<a name="Thetestcase"></a><h3>The testcase</h3>


我们下面实现的域、右手边和边界条件与地球物理学中的一个问题有关：在那里，人们想计算大洋中裂缝下地球内部岩浆的流动场。裂缝是两个大陆板块非常缓慢地漂移开来的地方（每年最多几厘米），在地壳上留下一个裂缝，被下面的岩浆填充。在不试图完全现实的情况下，我们通过求解域 $\Omega=[-2,2]\times[0,1]\times[-1,0]$ 上的以下一组方程和边界条件来模拟这种情况。

@f{eqnarray*}


  -2\; \textrm{div}\; \varepsilon(\textbf{u}) + \nabla p &=& 0,
  \\


  -\textrm{div}\; \textbf{u} &=& 0,
  \\
  \mathbf u &=&   \left(\begin{array}{c}


    -1 \\ 0 \\0
  \end{array}\right)
  \qquad\qquad \textrm{at}\ z=0, x<0,
  \\
  \mathbf u &=&   \left(\begin{array}{c}
    +1 \\ 0 \\0
  \end{array}\right)
  \qquad\qquad \textrm{at}\ z=0, x>0,
  \\
  \mathbf u &=&   \left(\begin{array}{c}
    0 \\ 0 \\0
  \end{array}\right)
  \qquad\qquad \textrm{at}\ z=0, x=0,


@f}

并在其他地方使用自然边界条件 $\textbf{n}\cdot [p \textbf{I} - 2
\varepsilon(\textbf{u})] = 0$ 。换句话说，在顶面的左边部分，我们规定流体以速度 $-1$ 随大陆板向左移动，在顶面的右边部分向右移动，并在其他地方施加自然流动条件。如果我们在2d中，描述基本上是相同的，只是我们省略了上述所有矢量的第二部分。

正如在<a href="#Results">results section</a>中会变得很明显的那样，流场将从下面拉动材料，并将其移动到域的左右两端，这是预期的。速度边界条件的不连续性将在顶面中心产生一个压力奇点，将材料一直吸到顶面，以填补材料在此位置向外运动所留下的缺口。




<a name="Implementation"></a><h3>Implementation</h3>


<a name="UsingimhomogeneousconstraintsforimplementingDirichletboundaryconditions"></a><h4>Using imhomogeneous constraints for implementing Dirichlet boundary conditions</h4>


在之前的所有教程程序中，我们仅仅使用AffineConstraints对象来处理悬挂节点约束（步骤11除外）。然而，这个类也可以用来实现Dirichlet边界条件，正如我们将在这个程序中展示的，通过固定一些节点值  $x_i = b_i$  。注意，这些是不均匀约束，我们要特别注意一些。我们要实现的方法是，首先通过调用AffineConstraints对象读入边界值

@code
  VectorTools::interpolate_boundary_values (dof_handler,
                                            1,
                                            BoundaryValues<dim>(),
                                            constraints);
@endcode



非常类似于我们之前制作边界节点列表的方式（注意，我们只在边界标志为1的边界上设置Dirichlet条件）。然后边界值的实际应用由AffineConstraints对象直接处理，没有任何额外的干扰。

然后我们可以像以前一样进行，即通过填充矩阵，然后在约束对象上调用一个浓缩函数，其形式为

@code
  constraints.condense (system_matrix, system_rhs);
@endcode



请注意，我们在系统矩阵和系统右侧同时调用，因为解决不均匀约束需要对矩阵条目和右侧的知识。但出于效率的考虑，我们选择了另一种策略：所有收集在AffineConstraints对象中的约束都可以在将本地数据写入全局矩阵的同时得到解决，方法是使用调用

@code
  constraints.distribute_local_to_global (local_matrix, local_rhs,
                                          local_dof_indices,
                                          system_matrix, system_rhs);
@endcode



这个技术在step-27教程程序中进一步讨论。我们在这里需要知道的是，这个函数同时做了三件事：它把局部数据写入全局矩阵和右手边，它分布了悬挂的节点约束，另外还实现了（不均匀的）迪里切特边界条件。这很好，不是吗？

我们可以得出结论，AffineConstraints类提供了一个替代使用 MatrixTools::apply_boundary_values 来实现Dirichlet边界条件的方法。


<a name="constraint-matrix">


<a name="UsingAffineConstraintsforincreasingperformanceh4h4UsingAffineConstraintsforincreasingperformance"></a><h4>Using AffineConstraints for increasing performance</h4><h4>Using AffineConstraints for increasing performance</h4>

</a> 。

通常，稀疏矩阵包含大量的元素，当我们要开始线性求解时，这些元素实际上是零。这样的元素是在我们消除约束条件或实现Dirichlet条件时引入的，我们通常会删除受约束行和列中的所有条目，即把它们设置为零。对于本教程程序中所考虑的三维应用，存在于稀疏模式中但并不真正包含任何信息的那部分元素，可以达到矩阵中元素总数的四分之一。请记住，矩阵-向量乘积或预处理程序对稀疏矩阵的所有元素（甚至那些为零的元素）进行操作，这是我们在这里要避免的低效率。

直接解决约束自由度的一个好处是，我们可以避免在我们的稀疏矩阵中出现大部分要为零的条目&mdash；在矩阵构建过程中，我们不需要约束的条目（与传统算法相反，传统算法是先填充矩阵，之后才解决约束）。这将在形成矩阵-向量乘积时节省内存和时间。我们要做的是将约束信息传递给生成稀疏模式的函数，然后设置一个<tt>false</tt>参数，指定我们不打算使用约束条目。

@code
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern,
                                   constraints, false);
@endcode

顺便说一下，这个函数也避免了对稀疏模式的<tt>condense()</tt>函数的调用。




<a name="Performanceoptimizations"></a><h4>Performance optimizations</h4>


下面开发的程序已经看到了很多的TLC。我们在分析工具（主要是<a
href="http://www.valgrind.org/">valgrind</a>的cachegrind和callgrind工具，以及用于可视化的KDE<a
href="http://kcachegrind.sourceforge.net/">KCachegrind</a>程序）下一遍又一遍地运行它，看看瓶颈在哪里。这已经得到了回报：通过这种努力，当考虑到细化周期0到3的运行时间时，程序的速度已经提高了大约4倍，将CPU指令的总体执行数量从869,574,060,348减少到199,853,005,625。对于更高的细化水平，收益可能更大，因为一些不是 ${\cal O}(N)$ 的算法被取消了。

基本上，目前程序中有两种算法不随自由度数量的增加而线性扩展：自由度的重新编号（即 ${\cal O}(N \log N)$  ，以及线性求解器（即 ${\cal O}(N^{4/3})$  ）。至于第一个，虽然自由度的重新排序可能不是线性扩展，但它是整个算法中不可缺少的部分，因为它极大地提高了稀疏ILU的质量，很容易弥补计算重新编号的时间；证明这一点的图表和时间显示在DoFRenumbering命名空间的文档中，也强调了下面选择的Cuthill-McKee重新排序算法。

至于线性求解器：如上所述，我们在这里的实现使用了Schur补码公式。这不一定是非常好的选择，但展示了deal.II中的各种重要技术。关于哪种求解器最好的问题，在本程序的<a
href="#improved-solver">section on improved solvers in the results part</a>中再次进行了讨论，同时还有显示备选求解器的代码和对其结果的比较。

除此以外，在这个程序的创建过程中，许多其他的算法也被测试和改进。例如，在建立稀疏性模式时，我们最初使用了一个（现在已经不存在了）BlockCompressedSparsityPattern对象，每次增加一个元素；然而，它的数据结构对于我们在3d中离散化所产生的每行大量的非零条目来说适应性很差，导致了一个二次方的行为。替换deal.II中的内部算法，一次设置许多元素，并使用BlockCompressedSimpleSparsityPattern（截至2015年初，它又被BlockDynamicSparsityPattern取代）作为一个更好的适应性数据结构，消除了这个瓶颈，代价是内存消耗略高。同样，SparseILU类中的分解步骤的实现也非常低效，已经被一个快10倍的步骤所取代。甚至SparseILU的vmult函数也得到了改进，节省了大约20%的时间。小的改进在这里和那里被应用。此外，AffineConstraints对象被用来消除稀疏矩阵中很多最终将为零的条目，见<a href="#constraint-matrix">the section on using advanced
features of the AffineConstraints class</a>。

这里显示了在3D的细化周期0到3中，CPU指令在程序中各个不同的地方花费了多少的概况。

 <img src="https://www.dealii.org/images/steps/developer/step-22.profile-3.png" alt=""> 

可以看出，在这个细化级别，大约四分之三的指令数花在实际求解上（左边的 SparseILU::vmult 调用，中间的 SparseMatrix::vmult 调用用于舒尔补码求解，另一个方框代表在求解<i>U</i>时与稀疏ILU和稀疏矩阵的乘法）。大约五分之一的指令数用于矩阵装配和稀疏ILU计算（右下角的方框），其余的用于其他方面。由于 SparseILU::vmult 调用中的浮点运算通常比矩阵装配中的许多逻辑运算和查表要长得多，矩阵装配所占用的运行时间的比例实际上大大低于指令的比例，这在我们在结果部分的比较中会变得很明显。

对于更高的细化水平，代表求解器的方框以及右上角源自重排算法的蓝色方框将以牺牲程序的其他部分为代价而增长，因为它们不是线性扩展。在这个中等细化水平（3168个单元和93176个自由度），线性求解器已经占了大约四分之三的指令，这是一个很好的迹象，说明这个程序中使用的大多数算法都是经过良好调整的，加快程序速度的主要改进很可能不是来自手工优化的个别方面，而是通过改变求解器的算法。我们将在下面的结果讨论中也讨论这一点。

最后一点，作为参考，下图也显示了在优化这个程序的早期阶段，配置文件的样子。

 <img src="https://www.dealii.org/images/steps/developer/step-22.profile-3.original.png" alt=""> 

如上所述，这个版本的运行时间大约是第一个配置文件的四倍，其中SparseILU分解占用了大约30%的指令数，而操作早期低效版本的DynamicSparsityPattern大约占10%。这两个瓶颈后来都被完全消除了。


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
 * 像往常一样，我们从包括一些著名的文件开始。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/utilities.h> 
 * 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_tools.h> 
 * #include <deal.II/grid/grid_refinement.h> 
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
 * 
 * @endcode
 * 
 * 然后我们需要包括稀疏直接求解器UMFPACK的头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * @endcode
 * 
 * 这包括不完全LU因子化的库，它将被用作3D的预处理程序。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_ilu.h> 
 * 
 * @endcode
 * 
 * 这是C++语言。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * #include <memory> 
 * 
 * @endcode
 * 
 * 和所有的程序一样，名字空间dealii被包括在内。
 * 

 * 
 * 
 * @code
 * namespace Step22 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="Definingtheinnerpreconditionertype"></a> 
 * <h3>Defining the inner preconditioner type</h3>
 * 

 * 
 * 正如介绍中所解释的，我们将分别对两个和三个空间维度使用不同的预处理程序。我们通过使用空间维度作为模板参数来区分它们。关于模板的细节，请参见 step-4 。我们不打算在这里创建任何预处理对象，我们所做的只是创建一个持有确定预处理类的本地别名的类，这样我们就可以以独立于维度的方式编写我们的程序。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   struct InnerPreconditioner; 
 * 
 * @endcode
 * 
 * 在二维中，我们将使用一个稀疏的直接求解器作为预处理程序。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   struct InnerPreconditioner<2> 
 *   { 
 *     using type = SparseDirectUMFPACK; 
 *   }; 
 * 
 * @endcode
 * 
 * 还有三维的ILU预处理，由SparseILU调用。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   struct InnerPreconditioner<3> 
 *   { 
 *     using type = SparseILU<double>; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ThecodeStokesProblemcodeclasstemplate"></a> 
 * <h3>The <code>StokesProblem</code> class template</h3>
 * 

 * 
 * 这是对 step-20 的改编，所以主类和数据类型与那里使用的几乎相同。唯一不同的是，我们有一个额外的成员  <code>preconditioner_matrix</code>  ，用于预处理Schur补码，以及一个相应的稀疏模式  <code>preconditioner_sparsity_pattern</code>  。此外，我们没有依赖LinearOperator，而是实现了我们自己的InverseMatrix类。
 * 

 * 
 * 在这个例子中，我们还使用了自适应网格细化，其处理方式与  step-6  类似。根据介绍中的讨论，我们也将使用AffineConstraints对象来实现Dirichlet边界条件。因此，我们改变名称  <code>hanging_node_constraints</code> into <code>constraints</code>  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class StokesProblem 
 *   { 
 *   public: 
 *     StokesProblem(const unsigned int degree); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_dofs(); 
 *     void assemble_system(); 
 *     void solve(); 
 *     void output_results(const unsigned int refinement_cycle) const; 
 *     void refine_mesh(); 
 * 
 *     const unsigned int degree; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> system_matrix; 
 * 
 *     BlockSparsityPattern      preconditioner_sparsity_pattern; 
 *     BlockSparseMatrix<double> preconditioner_matrix; 
 * 
 *     BlockVector<double> solution; 
 *     BlockVector<double> system_rhs; 
 * 
 * @endcode
 * 
 * 这一条是新的：我们将使用一个所谓的共享指针结构来访问预处理程序。共享指针本质上只是指针的一种方便形式。几个共享指针可以指向同一个对象（就像普通的指针一样），但是当最后一个指向前提器对象的共享指针对象被删除时（例如共享指针对象超出了范围，它所在的类被销毁，或者指针被分配给了不同的前提器对象），那么指向的前提器对象也被销毁。这确保了我们不必手动跟踪有多少地方仍在引用一个前置条件器对象，它永远不会产生内存泄漏，也不会产生一个指向已被销毁对象的悬空指针。
 * 

 * 
 * 
 * @code
 *     std::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Boundaryvaluesandrighthandside"></a> 
 * <h3>Boundary values and right hand side</h3>
 * 

 * 
 * 与 step-20 和其他大多数例子程序一样，下一个任务是定义PDE的数据：对于斯托克斯问题，我们将在部分边界上使用自然边界值（即同质诺伊曼型），对于这些边界，我们不必做任何特殊处理（同质性意味着弱形式中的相应项只是零），而在边界的其余部分使用速度的边界条件（迪里希勒型），如介绍中所述。
 * 

 * 
 * 为了强制执行速度上的Dirichlet边界值，我们将像往常一样使用 VectorTools::interpolate_boundary_values 函数，这要求我们写一个具有与有限元一样多分量的函数对象。换句话说，我们必须在 $(u,p)$ -空间上定义函数，但在插值边界值时，我们要过滤掉压力分量。
 * 

 * 
 * 下面的函数对象是介绍中描述的边界值的表示。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class BoundaryValues : public Function<dim> 
 *   { 
 *   public: 
 *     BoundaryValues() 
 *       : Function<dim>(dim + 1) 
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
 *   double BoundaryValues<dim>::value(const Point<dim> & p, 
 *                                     const unsigned int component) const 
 *   { 
 *     Assert(component < this->n_components, 
 *            ExcIndexRange(component, 0, this->n_components)); 
 * 
 *     if (component == 0) 
 *       return (p[0] < 0 ? -1 : (p[0] > 0 ? 1 : 0)); 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   void BoundaryValues<dim>::vector_value(const Point<dim> &p, 
 *                                          Vector<double> &  values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = BoundaryValues<dim>::value(p, c); 
 *   } 
 * 
 * @endcode
 * 
 * 我们为右手边实现类似的函数，在目前的例子中，右手边只是零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class RightHandSide : public Function<dim> 
 *   { 
 *   public: 
 *     RightHandSide() 
 *       : Function<dim>(dim + 1) 
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
 *   double RightHandSide<dim>::value(const Point<dim> & /*p*/, 
 *                                    const unsigned int /*component*/) const 
 *   { 
 *     return 0; 
 *   } 
 * 
 *   template <int dim> 
 *   void RightHandSide<dim>::vector_value(const Point<dim> &p, 
 *                                         Vector<double> &  values) const 
 *   { 
 *     for (unsigned int c = 0; c < this->n_components; ++c) 
 *       values(c) = RightHandSide<dim>::value(p, c); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="Linearsolversandpreconditioners"></a> 
 * <h3>Linear solvers and preconditioners</h3>
 * 

 * 
 * 在介绍中广泛讨论了线性求解器和预处理器。在这里，我们创建将被使用的各自对象。
 * 

 * 
 * 
 * <a name=""></a> 
 * @sect4{The <code>InverseMatrix</code> class template}   <code>InverseMatrix</code> 类表示逆矩阵的数据结构。与 step-20 不同，我们用一个类来实现，而不是用辅助函数inverse_linear_operator()，我们将把这个类应用于不同种类的矩阵，这些矩阵需要不同的预处理程序（在 step-20 中，我们只对质量矩阵使用非同一性预处理程序）。矩阵和预处理器的类型通过模板参数传递给这个类，当创建 <code>InverseMatrix</code> 对象时，这些类型的矩阵和预处理器对象将被传递给构造器。成员函数 <code>vmult</code> 是通过解决一个线性系统得到的。
 * 

 * 
 * 
 * @code
 *   template <class MatrixType, class PreconditionerType> 
 *   class InverseMatrix : public Subscriptor 
 *   { 
 *   public: 
 *     InverseMatrix(const MatrixType &        m, 
 *                   const PreconditionerType &preconditioner); 
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const; 
 * 
 *   private: 
 *     const SmartPointer<const MatrixType>         matrix; 
 *     const SmartPointer<const PreconditionerType> preconditioner; 
 *   }; 
 * 
 *   template <class MatrixType, class PreconditionerType> 
 *   InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix( 
 *     const MatrixType &        m, 
 *     const PreconditionerType &preconditioner) 
 *     : matrix(&m) 
 *     , preconditioner(&preconditioner) 
 *   {} 
 * 
 * @endcode
 * 
 * 这就是 <code>vmult</code> 函数的实现。
 * 

 * 
 * 在这个类中，我们对解算器控制使用了一个相当大的容忍度。这样做的原因是，该函数被频繁使用，因此，任何使CG求解中的残差变小的额外努力都会使求解更加昂贵。请注意，我们不仅将该类作为Schur补码的预处理程序，而且在形成拉普拉斯矩阵的逆时也使用该类；因此，该类直接对解本身的精度负责，所以我们也不能选择太大的公差。
 * 

 * 
 * 
 * @code
 *   template <class MatrixType, class PreconditionerType> 
 *   void InverseMatrix<MatrixType, PreconditionerType>::vmult( 
 *     Vector<double> &      dst, 
 *     const Vector<double> &src) const 
 *   { 
 *     SolverControl            solver_control(src.size(), 1e-6 * src.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     dst = 0; 
 * 
 *     cg.solve(*matrix, dst, src, *preconditioner); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="ThecodeSchurComplementcodeclasstemplate"></a> 
 * <h4>The <code>SchurComplement</code> class template</h4>
 * 

 * 
 * 这个类实现了介绍中讨论的Schur补码。它与  step-20  相类似。 不过，我们现在用一个模板参数 <code>PreconditionerType</code> 来调用它，以便在指定逆矩阵类的各自类型时访问它。作为上述定义的结果，声明  <code>InverseMatrix</code>  现在包含了上述预处理类的第二个模板参数，这也影响到  <code>SmartPointer</code> object <code>m_inverse</code>  。
 * 

 * 
 * 
 * @code
 *   template <class PreconditionerType> 
 *   class SchurComplement : public Subscriptor 
 *   { 
 *   public: 
 *     SchurComplement( 
 *       const BlockSparseMatrix<double> &system_matrix, 
 *       const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse); 
 * 
 *     void vmult(Vector<double> &dst, const Vector<double> &src) const; 
 * 
 *   private: 
 *     const SmartPointer<const BlockSparseMatrix<double>> system_matrix; 
 *     const SmartPointer< 
 *       const InverseMatrix<SparseMatrix<double>, PreconditionerType>> 
 *       A_inverse; 
 * 
 *     mutable Vector<double> tmp1, tmp2; 
 *   }; 
 * 
 *   template <class PreconditionerType> 
 *   SchurComplement<PreconditionerType>::SchurComplement( 
 *     const BlockSparseMatrix<double> &system_matrix, 
 *     const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse) 
 *     : system_matrix(&system_matrix) 
 *     , A_inverse(&A_inverse) 
 *     , tmp1(system_matrix.block(0, 0).m()) 
 *     , tmp2(system_matrix.block(0, 0).m()) 
 *   {} 
 * 
 *   template <class PreconditionerType> 
 *   void 
 *   SchurComplement<PreconditionerType>::vmult(Vector<double> &      dst, 
 *                                              const Vector<double> &src) const 
 *   { 
 *     system_matrix->block(0, 1).vmult(tmp1, src); 
 *     A_inverse->vmult(tmp2, tmp1); 
 *     system_matrix->block(1, 0).vmult(dst, tmp2); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemclassimplementation"></a> 
 * <h3>StokesProblem class implementation</h3>
 * 
 * <a name="StokesProblemStokesProblem"></a> 
 * <h4>StokesProblem::StokesProblem</h4>
 * 

 * 
 * 这个类的构造函数看起来与  step-20  的构造函数非常相似。构造函数初始化了多项式程度、三角形、有限元系统和dof处理器的变量。矢量速度分量的基础多项式函数的阶数为 <code>degree+1</code> ，压力的阶数为 <code>degree</code> 。 这就得到了LBB稳定元对 $Q_{degree+1}^d\times Q_{degree}$ ，通常被称为泰勒-霍德元。
 * 

 * 
 * 请注意，我们用MeshSmoothing参数初始化三角形，这可以确保单元的细化是以PDE解的近似保持良好的方式进行的（如果网格过于非结构化就会出现问题），详情请参见 <code>Triangulation::MeshSmoothing</code> 的文档。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   StokesProblem<dim>::StokesProblem(const unsigned int degree) 
 *     : degree(degree) 
 *     , triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1) 
 *     , dof_handler(triangulation) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsetup_dofs"></a> 
 * <h4>StokesProblem::setup_dofs</h4>
 * 

 * 
 * 给定一个网格，该函数将自由度与之关联，并创建相应的矩阵和向量。在开始的时候，它还释放了指向预处理对象的指针（如果共享指针在此时指向任何东西的话），因为在这之后肯定不会再需要它了，在组装矩阵后必须重新计算，并将稀疏矩阵从其稀疏模式对象中解开。
 * 

 * 
 * 然后，我们继续分配自由度并重新编号。为了使ILU预处理程序（在3D中）有效地工作，重要的是以这样的方式列举自由度，以减少矩阵的带宽，或者也许更重要的是：以这样的方式使ILU尽可能地接近于真正的LU分解。另一方面，我们需要保留在  step-20  和  step-21  中已经看到的速度和压力的块状结构。这将分两步完成。首先，对所有的道次进行重新编号，以改善ILU，然后我们再一次按组件重新编号。由于 <code>DoFRenumbering::component_wise</code> 没有触及单个块内的重新编号，所以第一步的基本重新编号仍然存在。至于如何对自由度进行重新编号以提高ILU：deal.II有许多算法试图找到排序以提高ILU，或减少矩阵的带宽，或优化其他方面。DoFRenumbering命名空间显示了我们在本教程程序中基于这里讨论的测试案例而获得的几种算法的结果比较。在这里，我们将使用传统的Cuthill-McKee算法，该算法已经在之前的一些教程程序中使用。 在<a href="#improved-ilu">section on improved ILU</a>中我们将更详细地讨论这个问题。
 * 与以前的教程程序相比，
 * 还有一个变化。没有理由对 <code>dim</code> 的速度成分进行单独排序。事实上，与其先列举所有 $x$ -velocities，再列举所有 $y$ -velocities，等等，我们希望将所有速度放在一起，只在速度（所有分量）和压力之间分开。默认情况下， DoFRenumbering::component_wise 函数不是这样做的：它把每个矢量分量分开处理；我们要做的是把几个分量分成 "块"，并把这个块结构传递给该函数。因此，我们分配一个矢量 <code>block_component</code> ，有多少个元素就有多少个分量，描述所有的速度分量对应于块0，而压力分量将形成块1。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::setup_dofs() 
 *   { 
 *     A_preconditioner.reset(); 
 *     system_matrix.clear(); 
 *     preconditioner_matrix.clear(); 
 * 
 *     dof_handler.distribute_dofs(fe); 
 *     DoFRenumbering::Cuthill_McKee(dof_handler); 
 * 
 *     std::vector<unsigned int> block_component(dim + 1, 0); 
 *     block_component[dim] = 1; 
 *     DoFRenumbering::component_wise(dof_handler, block_component); 
 * 
 * @endcode
 * 
 * 现在是对Dirichlet边界条件的实现，在介绍中的讨论之后，这应该是很明显的。所有的变化是，这个函数已经出现在设置函数中，而我们习惯于在一些汇编例程中看到它。在我们设置网格的下面，我们将把施加Dirichlet边界条件的顶部边界与边界指标1联系起来。 我们必须将这个边界指标作为第二个参数传递给下面的插值函数。 不过，还有一件事。 描述Dirichlet条件的函数是为所有分量定义的，包括速度和压力。然而，Dirichlet条件只为速度而设置。 为此，我们使用一个只选择速度分量的ComponentMask。通过指定我们想要的特定分量，从有限元中获得该分量掩码。由于我们使用自适应细化网格，仿生约束对象需要首先填充由DoF处理程序生成的悬挂节点约束。注意这两个函数的顺序；我们首先计算悬挂节点约束，然后将边界值插入约束对象。这确保了我们在有悬挂节点的边界上尊重H<sup>1</sup>一致性（在三个空间维度上），悬挂节点需要支配Dirichlet边界值。
 * 

 * 
 * 
 * @code
 *     { 
 *       constraints.clear(); 
 * 
 *       FEValuesExtractors::Vector velocities(0); 
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                1, 
 *                                                BoundaryValues<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 *     } 
 * 
 *     constraints.close(); 
 * 
 * @endcode
 * 
 * 与 step-20 相类似，我们计算各个组件中的道夫。我们可以用与那里相同的方式来做，但我们想在我们已经用于重新编号的块结构上进行操作。函数  <code>DoFTools::count_dofs_per_fe_block</code>  的作用与  <code>DoFTools::count_dofs_per_fe_component</code>  相同，但现在通过  <code>block_component</code>  将速度和压力块分组。
 * 

 * 
 * 
 * @code
 *     const std::vector<types::global_dof_index> dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 *     const unsigned int n_u = dofs_per_block[0]; 
 *     const unsigned int n_p = dofs_per_block[1]; 
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (" << n_u << '+' << n_p << ')' << std::endl; 
 * 
 * @endcode
 * 
 * 下一个任务是为我们将创建的系统矩阵分配一个稀疏模式，为预处理矩阵分配一个稀疏模式。我们可以用与 step-20 相同的方式来做这件事，即通过 DoFTools::make_sparsity_pattern. 直接建立一个SparsityPattern类型的对象，但是，有一个重要的理由不这样做。在3D中，函数 DoFTools::max_couplings_between_dofs 对各个道夫之间的耦合产生了一个保守但相当大的数字，因此，最初为创建矩阵的稀疏模式提供的内存太多--实际上，对于中等大小的3D问题，初始稀疏模式甚至无法放入大多数系统的物理内存中，也请参见 step-18  中的讨论。相反，我们首先建立临时对象，使用不同的数据结构，不需要分配更多的内存，但不适合作为SparseMatrix或BlockSparseMatrix对象的基础；在第二步，我们将这些对象复制到BlockSparsityPattern类型的对象中。这完全类似于我们在  step-11  和  step-18  中已经做过的事情。特别是，我们利用了这样一个事实，即我们永远不会写入系统矩阵的 $(1,1)$ 块中，而且这是唯一需要填充的预处理矩阵块。
 * 

 * 
 * 所有这些都是在新范围内完成的，这意味着一旦信息被复制到  <code>sparsity_pattern</code>  ，  <code>dsp</code>  的内存将被释放。
 * 

 * 
 * 
 * @code
 *     { 
 *       BlockDynamicSparsityPattern dsp(2, 2); 
 * 
 *       dsp.block(0, 0).reinit(n_u, n_u); 
 *       dsp.block(1, 0).reinit(n_p, n_u); 
 *       dsp.block(0, 1).reinit(n_u, n_p); 
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
 *         dof_handler, coupling, dsp, constraints, false); 
 * 
 *  
 *     } 
 * 
 *     { 
 *       BlockDynamicSparsityPattern preconditioner_dsp(2, 2); 
 * 
 *       preconditioner_dsp.block(0, 0).reinit(n_u, n_u); 
 *       preconditioner_dsp.block(1, 0).reinit(n_p, n_u); 
 *       preconditioner_dsp.block(0, 1).reinit(n_u, n_p); 
 *       preconditioner_dsp.block(1, 1).reinit(n_p, n_p); 
 * 
 *       preconditioner_dsp.collect_sizes(); 
 * 
 *       Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1); 
 * 
 *       for (unsigned int c = 0; c < dim + 1; ++c) 
 *         for (unsigned int d = 0; d < dim + 1; ++d) 
 *           if (((c == dim) && (d == dim))) 
 *             preconditioner_coupling[c][d] = DoFTools::always; 
 *           else 
 *             preconditioner_coupling[c][d] = DoFTools::none; 
 * 
 *       DoFTools::make_sparsity_pattern(dof_handler, 
 *                                       preconditioner_coupling, 
 *                                       preconditioner_dsp, 
 *                                       constraints, 
 *                                       false); 
 * 
 *       preconditioner_sparsity_pattern.copy_from(preconditioner_dsp); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，与  step-20  中的方法类似，从块状结构中创建系统矩阵、前导矩阵、解决方案和右侧向量。
 * 

 * 
 * 
 * @code
 *     system_matrix.reinit(sparsity_pattern); 
 *     preconditioner_matrix.reinit(preconditioner_sparsity_pattern); 
 * 
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
 * <a name="StokesProblemassemble_system"></a> 
 * <h4>StokesProblem::assemble_system</h4>
 * 

 * 
 * 汇编过程遵循 step-20 和介绍中的讨论。我们使用众所周知的缩写来表示保存本单元自由度的局部矩阵、右手边和全局编号的数据结构。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::assemble_system() 
 *   { 
 *     system_matrix         = 0; 
 *     system_rhs            = 0; 
 *     preconditioner_matrix = 0; 
 * 
 *     QGauss<dim> quadrature_formula(degree + 2); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values | update_gradients); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 * 
 *     const unsigned int n_q_points = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> local_preconditioner_matrix(dofs_per_cell, 
 *                                                    dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const RightHandSide<dim>    right_hand_side; 
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 
 * 
 * @endcode
 * 
 * 接下来，我们需要两个对象，作为FEValues对象的提取器。它们的用途在  @ref  vector_valued 的报告中详细解释。
 * 

 * 
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 * @endcode
 * 
 * 作为对 step-20 和 step-21 的扩展，我们包括了一些优化，使这个特定问题的装配速度大大加快。这些改进是基于这样的观察：当我们像 step-20 那样做时，我们做了太多次的计算：对称梯度实际上在每个正交点有 <code>dofs_per_cell</code> 个不同的值，但是我们从FEValues对象中提取了 <code>dofs_per_cell*dofs_per_cell</code> 次。
 * 

 * 
 * - 在 <code>i</code> 的循环和 <code>j</code> 的内循环中。在3D中，这意味着评估它 $89^2=7921$ 次而不是 $89$ 次，这是一个不小的差别。
 * 

 * 
 * 所以我们在这里要做的是，在开始对单元上的道夫进行循环之前，在正交点得到一个秩-2张量的向量（类似的还有压力上的发散和基函数值）来避免这种重复计算。首先，我们创建各自的对象来保存这些值。然后，我们开始在所有单元上进行循环，并在正交点上进行循环，在那里我们首先提取这些值。我们在这里还实现了一个优化：本地矩阵（以及全局矩阵）将是对称的，因为所有涉及的操作都是相对于 $i$ 和 $j$ 对称的。这可以通过简单地运行内循环而不是 <code>dofs_per_cell</code>, but only up to <code>i</code> 来实现，即外循环的索引。
 * 

 * 
 * 
 * @code
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell); 
 *     std::vector<double>                  div_phi_u(dofs_per_cell); 
 *     std::vector<double>                  phi_p(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         local_matrix                = 0; 
 *         local_preconditioner_matrix = 0; 
 *         local_rhs                   = 0; 
 * 
 *         right_hand_side.vector_value_list(fe_values.get_quadrature_points(), 
 *                                           rhs_values); 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 symgrad_phi_u[k] = 
 *                   fe_values[velocities].symmetric_gradient(k, q); 
 *                 div_phi_u[k] = fe_values[velocities].divergence(k, q); 
 *                 phi_p[k]     = fe_values[pressure].value(k, q); 
 *               } 
 * 
 * @endcode
 * 
 * 最后是系统矩阵和我们用于预处理程序的矩阵的双线性形式。回顾一下，这两个的公式分别是
 * @f{align*}{
 * A_{ij} &= a(\varphi_i,\varphi_j)
 * \\     &= \underbrace{2(\varepsilon(\varphi_{i,\textbf{u}}),
 * \varepsilon(\varphi_{j,\textbf{u}}))_{\Omega}}
 * _{(1)}
 * \;
 * \underbrace{- (\textrm{div}\; \varphi_{i,\textbf{u}},
 * \varphi_{j,p})_{\Omega}}
 * _{(2)}
 * \;
 * \underbrace{- (\varphi_{i,p},
 * \textrm{div}\;
 * \varphi_{j,\textbf{u}})_{\Omega}}
 * _{(3)}
 * @f}
 * 和
 * @f{align*}{
 * M_{ij} &= \underbrace{(\varphi_{i,p},
 * \varphi_{j,p})_{\Omega}}
 * _{(4)},
 * @f} ， 
 * 其中 $\varphi_{i,\textbf{u}}$ 和 $\varphi_{i,p}$ 是 $i$ th形状函数的速度和压力成分。然后，上述各种术语在下面的实现中很容易识别。
 * 

 * 
 * 
 * @code
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               { 
 *                 for (unsigned int j = 0; j <= i; ++j) 
 *                   { 
 *                     local_matrix(i, j) += 
 *                       (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) // (1) 
 *                        - div_phi_u[i] * phi_p[j]                 // (2) 
 *                        - phi_p[i] * div_phi_u[j])                // (3) 
 *                       * fe_values.JxW(q);                        // * dx 
 * 
 *                     local_preconditioner_matrix(i, j) += 
 *                       (phi_p[i] * phi_p[j]) // (4) 
 *                       * fe_values.JxW(q);   // * dx 
 *                   } 
 * 
 * @endcode
 * 
 * 注意在上述（1）的实现中，`operator*`被重载用于对称张量，产生两个张量之间的标量乘积。            对于右手边，我们利用形状函数只在一个分量中不为零的事实（因为我们的元素是原始的）。 我们不是将代表形状函数i的dim+1值的张量与整个右手边的向量相乘，而是只看唯一的非零分量。函数 FiniteElement::system_to_component_index 将返回这个形状函数所处的分量（0=x速度，1=y速度，2=2d中的压力），我们用它来挑选出右手边向量的正确分量来相乘。
 * 

 * 
 * 
 * @code
 *                 const unsigned int component_i = 
 *                   fe.system_to_component_index(i).first; 
 *                 local_rhs(i) += (fe_values.shape_value(i, q)   // (phi_u_i(x_q) 
 *                                  * rhs_values[q](component_i)) // * f(x_q)) 
 *                                 * fe_values.JxW(q);            // * dx 
 *               } 
 *           } 
 * 
 * @endcode
 * 
 * 在我们将局部数据写入全局矩阵之前（同时使用AffineConstraints对象来应用Dirichlet边界条件并消除悬挂的节点约束，正如我们在介绍中讨论的那样），我们必须注意一件事。由于对称性，我们只建立了一半的局部矩阵，但我们要保存完整的矩阵，以便使用标准函数进行解算。这是通过翻转指数来实现的，以防我们指向本地矩阵的空部分。
 * 

 * 
 * 
 * @code
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
 *             { 
 *               local_matrix(i, j) = local_matrix(j, i); 
 *               local_preconditioner_matrix(i, j) = 
 *                 local_preconditioner_matrix(j, i); 
 *             } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global(local_matrix, 
 *                                                local_rhs, 
 *                                                local_dof_indices, 
 *                                                system_matrix, 
 *                                                system_rhs); 
 *         constraints.distribute_local_to_global(local_preconditioner_matrix, 
 *                                                local_dof_indices, 
 *                                                preconditioner_matrix); 
 *       } 
 * 
 * @endcode
 * 
 * 在我们要解决这个线性系统之前，我们为速度-速度矩阵生成一个预处理程序，即系统矩阵中的 <code>block(0,0)</code> 。如上所述，这取决于空间维度。由于 <code>InnerPreconditioner::type</code> 别名所描述的两个类具有相同的接口，因此无论我们想使用稀疏直接求解器还是ILU，都不需要做任何不同的事情。
 * 

 * 
 * 
 * @code
 *     std::cout << "   Computing preconditioner..." << std::endl << std::flush; 
 * 
 *     A_preconditioner = 
 *       std::make_shared<typename InnerPreconditioner<dim>::type>(); 
 *     A_preconditioner->initialize( 
 *       system_matrix.block(0, 0), 
 *       typename InnerPreconditioner<dim>::type::AdditionalData()); 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsolve"></a> 
 * <h4>StokesProblem::solve</h4>
 * 

 * 
 * 经过前面介绍中的讨论和各自类的定义， <code>solve</code> 函数的实现是相当直接的，其方式与 step-20 类似。首先，我们需要一个 <code>InverseMatrix</code> 类的对象，代表矩阵A的逆。正如在介绍中所描述的，在  <code>InnerPreconditioner::type</code>  类型的内部预处理器的帮助下，生成了逆。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::solve() 
 *   { 
 *     const InverseMatrix<SparseMatrix<double>, 
 *                         typename InnerPreconditioner<dim>::type> 
 *                    A_inverse(system_matrix.block(0, 0), *A_preconditioner); 
 *     Vector<double> tmp(solution.block(0).size()); 
 * 
 * @endcode
 * 
 * 这与  step-20  中的情况一样。我们生成 Schur 补数的右手边  $B A^{-1} F - G$  和一个代表各自线性运算的对象  $B A^{-1} B^T$  ，现在有一个模板参数表示预处理器
 * 

 * 
 * - 按照类的定义。
 * 

 * 
 * 
 * @code
 *     { 
 *       Vector<double> schur_rhs(solution.block(1).size()); 
 *       A_inverse.vmult(tmp, system_rhs.block(0)); 
 *       system_matrix.block(1, 0).vmult(schur_rhs, tmp); 
 *       schur_rhs -= system_rhs.block(1); 
 * 
 *       SchurComplement<typename InnerPreconditioner<dim>::type> schur_complement( 
 *         system_matrix, A_inverse); 
 * 
 * @endcode
 * 
 * 解算器调用的常规控制结构被创建...
 * 

 * 
 * 
 * @code
 *       SolverControl            solver_control(solution.block(1).size(), 
 *                                    1e-6 * schur_rhs.l2_norm()); 
 *       SolverCG<Vector<double>> cg(solver_control); 
 * 
 * @endcode
 * 
 * 现在是对舒尔补码的预处理。正如介绍中所解释的，预处理是由压力变量的质量矩阵来完成的。
 * 

 * 
 * 实际上，求解器需要有 $P^{-1}$ 形式的预处理，所以我们需要创建一个逆运算。我们再次使用一个 <code>InverseMatrix</code> 类的对象，它实现了求解器需要的 <code>vmult</code> 操作。 在这种情况下，我们必须对压力质量矩阵进行反转。正如在早期的教程程序中已经证明的那样，质量矩阵的反转是一个相当便宜和简单的操作（与拉普拉斯矩阵等相比）。带有ILU预处理的CG方法在5-10步内收敛，与网格大小无关。 这正是我们在这里所做的。我们选择另一个ILU预处理，并通过相应的模板参数将其带入InverseMatrix对象。 然后在逆矩阵的vmult操作中调用一个CG求解器。
 * 

 * 
 * 另一种方法是选择因子为1.2的SSOR预处理器，这种方法构建成本较低，但之后需要更多的迭代。它需要大约两倍的迭代次数，但其生成的成本几乎可以忽略不计。
 * 

 * 
 * 
 * @code
 *       SparseILU<double> preconditioner; 
 *       preconditioner.initialize(preconditioner_matrix.block(1, 1), 
 *                                 SparseILU<double>::AdditionalData()); 
 * 
 *       InverseMatrix<SparseMatrix<double>, SparseILU<double>> m_inverse( 
 *         preconditioner_matrix.block(1, 1), preconditioner); 
 * 
 * @endcode
 * 
 * 有了舒尔补码和高效的预处理程序，我们可以用通常的方法解决压力的相关方程（即解向量中的0块）。
 * 

 * 
 * 
 * @code
 *       cg.solve(schur_complement, solution.block(1), schur_rhs, m_inverse); 
 * 
 * @endcode
 * 
 * 在这第一个求解步骤之后，必须将悬挂的节点约束分布到求解中，以实现一致的压力场。
 * 

 * 
 * 
 * @code
 *       constraints.distribute(solution); 
 * 
 *       std::cout << "  " << solver_control.last_step() 
 *                 << " outer CG Schur complement iterations for pressure" 
 *                 << std::endl; 
 *     } 
 * 
 * @endcode
 * 
 * 和 step-20 一样，我们最后需要解速度方程，在这里我们插入压力方程的解。这只涉及我们已经知道的对象
 * 

 * 
 * 所以我们只需用 $p$ 乘以 $B^T$ ，减去右边的部分，再乘以 $A$ 的逆数。最后，我们需要分配悬挂节点的约束，以获得一个一致的流场。
 * 

 * 
 * 
 * @code
 *     { 
 *       system_matrix.block(0, 1).vmult(tmp, solution.block(1)); 
 *       tmp *= -1; 
 *       tmp += system_rhs.block(0); 
 * 
 *       A_inverse.vmult(solution.block(0), tmp); 
 * 
 *       constraints.distribute(solution); 
 *     } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemoutput_results"></a> 
 * <h4>StokesProblem::output_results</h4>
 * 

 * 
 * 下一个函数生成图形输出。在这个例子中，我们将使用VTK文件格式。 我们给问题中的各个变量附上名字： <code>velocity</code> to the <code>dim</code> 速度的组成部分和 <code>pressure</code> 压力的组成部分。
 * 

 * 
 * 并非所有的可视化程序都有能力将各个矢量分量组合成一个矢量来提供矢量图；特别是对于一些基于VTK的可视化程序来说，这一点是成立的。在这种情况下，在包含数据的文件中应该已经描述了组件的逻辑分组为矢量的情况。换句话说，我们需要做的是为我们的输出编写者提供一种方法，让他们知道有限元的哪些分量在逻辑上形成一个矢量（在 $d$ 空间维度上有 $d$ 分量），而不是让他们假设我们只是有一堆标量场。 这是用 <code>DataComponentInterpretation</code> 命名空间的成员实现的：和文件名一样，我们创建一个矢量，其中第一个 <code>dim</code> 分量指的是速度，并被赋予 DataComponentInterpretation::component_is_part_of_vector; 标签，我们最后推一个标签 DataComponentInterpretation::component_is_scalar 来描述压力变量的分组。
 * 

 * 
 * 然后函数的其余部分与  step-20  中的相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void 
 *   StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const 
 *   { 
 *     std::vector<std::string> solution_names(dim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 * 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output( 
 *       "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemrefine_mesh"></a> 
 * <h4>StokesProblem::refine_mesh</h4>
 * 

 * 
 * 这是 <code>StokesProblem</code> 类中最后一个有趣的函数。 正如它的名字所示，它获取问题的解决方案，并在需要时细化网格。其过程与 step-6 中的相应步骤相同，不同的是我们只根据压力的变化进行细化，也就是说，我们用ComponentMask类型的掩码对象调用Kelly误差估计器，选择我们感兴趣的压力的单一标量分量（我们通过指定我们想要的分量从有限元类中得到这样一个掩码）。此外，我们没有再次粗化网格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::refine_mesh() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 * 
 *     FEValuesExtractors::Scalar pressure(dim); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       solution, 
 *       estimated_error_per_cell, 
 *       fe.component_mask(pressure)); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.0); 
 *     triangulation.execute_coarsening_and_refinement(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemrun"></a> 
 * <h4>StokesProblem::run</h4>
 * 

 * 
 * 在斯托克斯类中的最后一步，像往常一样，是生成初始网格的函数，并按各自的顺序调用其他函数。
 * 

 * 
 * 我们从一个大小为 $4 \times 1$ （2D）或 $4 \times 1 \times 1$ （3D）的矩形开始，在 $R^2/R^3$ 中分别放置为 $(-2,2)\times(-1,0)$ 或 $(-2,2)\times(0,1)\times(-1,0)$  。在每个方向上以相等的网格大小开始是很自然的，所以我们在第一个坐标方向上将初始矩形细分四次。为了将创建网格所涉及的变量的范围限制在我们实际需要的范围内，我们将整个块放在一对大括号之间。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::run() 
 *   { 
 *     { 
 *       std::vector<unsigned int> subdivisions(dim, 1); 
 *       subdivisions[0] = 4; 
 * 
 *       const Point<dim> bottom_left = (dim == 2 ?                // 
 *                                         Point<dim>(-2, -1) :    // 2d case 
 *                                         Point<dim>(-2, 0, -1)); // 3d case 
 * 
 *       const Point<dim> top_right = (dim == 2 ?              // 
 *                                       Point<dim>(2, 0) :    // 2d case 
 *                                       Point<dim>(2, 1, 0)); // 3d case 
 * 
 *       GridGenerator::subdivided_hyper_rectangle(triangulation, 
 *                                                 subdivisions, 
 *                                                 bottom_left, 
 *                                                 top_right); 
 *     } 
 * 
 * @endcode
 * 
 * 边界指标1被设置为所有受Dirichlet边界条件约束的边界，即位于最后一个坐标方向上的0的面。详见上面的例子描述。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : triangulation.active_cell_iterators()) 
 *       for (const auto &face : cell->face_iterators()) 
 *         if (face->center()[dim - 1] == 0) 
 *           face->set_all_boundary_ids(1); 
 * 
 * @endcode
 * 
 * 然后，在第一次求解之前，我们应用一个初始细化。在3D中，会有更多的自由度，所以我们在那里细化得更少。
 * 

 * 
 * 
 * @code
 *     triangulation.refine_global(4 - dim); 
 * 
 * @endcode
 * 
 * 正如在 step-6 中第一次看到的那样，我们在不同的细化级别上循环细化（除了第一个循环），设置自由度和矩阵，组装，求解和创建输出。
 * 

 * 
 * 
 * @code
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 6; 
 *          ++refinement_cycle) 
 *       { 
 *         std::cout << "Refinement cycle " << refinement_cycle << std::endl; 
 * 
 *         if (refinement_cycle > 0) 
 *           refine_mesh(); 
 * 
 *         setup_dofs(); 
 * 
 *         std::cout << "   Assembling..." << std::endl << std::flush; 
 *         assemble_system(); 
 * 
 *         std::cout << "   Solving..." << std::flush; 
 *         solve(); 
 * 
 *         output_results(refinement_cycle); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step22 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 主函数与  step-20  中的相同。我们将元素度数作为参数传递，并在众所周知的模板槽中选择空间尺寸。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step22; 
 * 
 *       StokesProblem<2> flow_problem(1); 
 *       flow_problem.run(); 
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
examples/step-22/doc/results.dox

<a name="Results"></a>

<a name="Results"></a><h1>Results</h1>


<a name="Outputoftheprogramandgraphicalvisualization"></a><h3>Output of the program and graphical visualization</h3>


<a name="2Dcalculations"></a><h4>2D calculations</h4>


在 <code>main</code> 函数中空间维度设置为2的情况下运行程序，会产生以下输出（在 "释放模式 "下， @dealiiVideoLectureSeeAlso{18}): ）。

@code
examples/\step-22> make run
Refinement cycle 0
   Number of active cells: 64
   Number of degrees of freedom: 679 (594+85)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure


Refinement cycle 1
   Number of active cells: 160
   Number of degrees of freedom: 1683 (1482+201)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure


Refinement cycle 2
   Number of active cells: 376
   Number of degrees of freedom: 3813 (3370+443)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure


Refinement cycle 3
   Number of active cells: 880
   Number of degrees of freedom: 8723 (7722+1001)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure


Refinement cycle 4
   Number of active cells: 2008
   Number of degrees of freedom: 19383 (17186+2197)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure


Refinement cycle 5
   Number of active cells: 4288
   Number of degrees of freedom: 40855 (36250+4605)
   Assembling...
   Computing preconditioner...
   Solving...  11 outer CG Schur complement iterations for pressure
@endcode



上述整个计算在一台相当快的（以2015年的标准）机器上需要大约2秒。

我们立即看到的是，（外部）迭代的次数并没有随着我们细化网格而增加。这证实了介绍中的说法，即用质量矩阵对Schur补码进行预处理，确实可以得到一个与身份矩阵频谱等价的矩阵（即特征值上下受限，与网格大小或单元的相对大小无关）。换句话说，质量矩阵和Schur补码在光谱上是等价的。

在下面的图片中，我们展示了程序中前六个细化步骤的网格。  观察一下网格是如何在解迅速变化的区域被细化的。在上边界，我们的迪里希特边界条件在左半边是-1，右半边是1，所以在 $x=0$ 有一个突然的变化。同样地，在两个上角也有从Dirichlet到Neumann数据的变化，所以那里也需要细化。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-0.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-1.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-3.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-4.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.2d.mesh-5.png" alt="">
    </td>
  </tr>
</table> 

最后，下面是一个流场图。它显示了流体随着上层边界的移动而被来自下层的物质所取代。

 <img src="https://www.dealii.org/images/steps/developer/step-22.2d.solution.png" alt=""> 

该图使用了基于VTK的可视化程序的能力（在这种情况下是VisIt）来显示矢量数据；这是我们将使用中的有限元的速度分量声明为一组矢量分量，而不是本教程程序的 <code>StokesProblem@<dim@>::%output_results</code> 函数中的独立标量分量的结果。




<a name="3Dcalculations"></a><h4>3D calculations</h4>


在3D中，程序的屏幕输出看起来像这样。

@code
Refinement cycle 0
   Number of active cells: 32
   Number of degrees of freedom: 1356 (1275+81)
   Assembling...
   Computing preconditioner...
   Solving...  13 outer CG Schur complement iterations for pressure.


Refinement cycle 1
   Number of active cells: 144
   Number of degrees of freedom: 5088 (4827+261)
   Assembling...
   Computing preconditioner...
   Solving...  14 outer CG Schur complement iterations for pressure.


Refinement cycle 2
   Number of active cells: 704
   Number of degrees of freedom: 22406 (21351+1055)
   Assembling...
   Computing preconditioner...
   Solving...  14 outer CG Schur complement iterations for pressure.


Refinement cycle 3
   Number of active cells: 3168
   Number of degrees of freedom: 93176 (89043+4133)
   Assembling...
   Computing preconditioner...
   Solving...  15 outer CG Schur complement iterations for pressure.


Refinement cycle 4
   Number of active cells: 11456
   Number of degrees of freedom: 327808 (313659+14149)
   Assembling...
   Computing preconditioner...
   Solving...  15 outer CG Schur complement iterations for pressure.


Refinement cycle 5
   Number of active cells: 45056
   Number of degrees of freedom: 1254464 (1201371+53093)
   Assembling...
   Computing preconditioner...
   Solving...  14 outer CG Schur complement iterations for pressure.
@endcode



我们再次看到，随着我们对网格的细化，外迭代的次数并没有增加。然而，计算时间明显增加：对于上述每个迭代分别需要约0.14秒、0.63秒、4.8秒、35秒、2分33秒和13分12秒。这种运行时间的整体超线性（未知数的数量）增加是由于我们的内部求解器不是 ${\cal O}(N)$ ：一个简单的实验表明，随着我们不断细化网格，反演速度-速度块 $A$ 的ILU预处理的CG平均迭代次数会增加。

我们将解决如何可能改进我们的解算器<a
href="#improved-solver">below</a>的问题。

至于图形输出，在解决过程中产生的网格看起来如下。

 <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-0.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-1.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-2.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-3.png" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-4.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d.mesh-5.png" alt="">
    </td>
  </tr>
</table> 

同样，它们基本上显示了由边界条件引入的奇异点的位置。计算出的矢量场构成了一个有趣的图形。

 <img src="https://www.dealii.org/images/steps/developer/step-22.3d.solution.png" alt=""> 

这里显示的等值线也是压力变量的等值线，显示了在不连续的速度边界条件下的奇异现象。




<a name="Sparsitypattern"></a><h3>Sparsity pattern</h3>


正如在生成稀疏模式时解释的那样，在使用不完全LU分解等预处理程序时，必须牢记自由度的编号。使用刚度矩阵中的非零元素的分布可以最方便地看到这一点。

如果我们不对自由度的重新编号做任何特殊处理（即不使用 DoFRenumbering::Cuthill_McKee, ，而是使用 DoFRenumbering::component_wise 来确保自由度被适当地排序到矩阵和向量的相应块中），那么我们在二维的第一次自适应细化后得到以下图像。

 <img src="https://www.dealii.org/images/steps/developer/step-22.2d.sparsity-nor.png" alt=""> 

为了生成这样的图表，你必须在设置步骤的末尾插入一段类似以下的代码。

@code
  {
    std::ofstream out ("sparsity_pattern.gpl");
    sparsity_pattern.print_gnuplot(out);
  }
@endcode



可以清楚地看到，非零项几乎分布在整个矩阵中。  这使得用ILU进行预处理的效率很低。ILU生成的高斯消除（LU分解）没有填充元素，这意味着更多的暂定填充元素的遗漏将导致完整分解的近似度降低。

因此，在这个程序中，我们选择了一种更高级的元件重新编号的方法。  用 DoFRenumbering::Cuthill_McKee 进行重新编号，并将元件分组为速度和压力，得到以下输出。

 <img src="https://www.dealii.org/images/steps/developer/step-22.2d.sparsity-ren.png" alt=""> 

很明显，情况有了很大的改善。现在大多数元素都集中在矩阵中（0,0）块的对角线周围。其他区块也可以看到类似的效果。在这种情况下，ILU分解将更接近于完全的LU分解，这就提高了预处理程序的质量。值得注意的是，稀疏直接求解器UMFPACK在实际生成稀疏LU分解之前会对方程进行一些内部重新编号；这个过程导致了与我们从Cuthill-McKee算法中得到的模式非常相似）。

最后，我们想仔细看看三维的稀疏模式。我们只显示了矩阵的(0,0)块，还是在一次自适应细化之后。除了矩阵大小增加的事实之外，还可以看到矩阵中多了很多条目。此外，即使是优化后的重新编号，也会有相当数量的暂定填充元素。这说明了为什么UMFPACK在三维中不是一个好的选择--完整的分解需要许多新的条目，最终无法放入物理内存（RAM）。

 <img src="https://www.dealii.org/images/steps/developer/step-22.3d.sparsity_uu-ren.png" alt=""> 




<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="improved-solver">


<a name="Improvedlinearsolverin3Dh4h4Improvedlinearsolverin3D"></a><h4>Improved linear solver in 3D</h4><h4>Improved linear solver in 3D</h4>

</a> 。

我们在计算结果一节中看到，外迭代次数不取决于网格大小，从可扩展性的角度看，这是最佳的。然而，如上所述，这并不适用于整个求解器。在生成矩阵 $A$ 和质量矩阵 $M_p$ 的逆时，我们没有研究内部迭代的数量。当然，这在二维情况下是没有问题的，我们用直接求解器对 $A$ 进行预处理，逆矩阵结构的 <code>vmult</code> 操作将在一个单一的CG步骤中收敛，但这在三维情况下发生变化，我们只使用ILU预处理器。  在那里，反演 $A$ 所需的预处理CG步骤的数量随着网格的细化而增加，每个 <code>vmult</code> 操作平均涉及上述细化步骤中的大约14、23、36、59、75和101次内部CG迭代。(另一方面，应用反压力质量矩阵的迭代次数总是在5次左右，在二维和三维中都是如此)。  总而言之，大部分的工作都花在了用相同的矩阵 $A$ 反复解决线性系统上。让这看起来更糟糕的是，我们实际上反转了一个矩阵，其大小约为整个系统矩阵的95%，并代表了稀疏模式中85%的非零条目。因此，自然的问题是，在计算块状系统的解时，用矩阵 $A$ 解约15次的线性系统是否合理。

当然，答案是，我们可以用其他一些（大多数时候是更好的）方法来做。然而，必须指出的是，正如我们在早期教程中所看到的那样，眼前的这个不确定系统对线性代数的要求确实比标准椭圆问题高得多。如果与类似规模的椭圆问题相比，其改进仍然相当不令人满意。无论如何，我们将在下面介绍一些对线性求解器的改进，我们将在第31步程序中用额外的选项再次重新考虑这一讨论。

<a name="improved-ilu">


<a name="BetterILUdecompositionbysmartreorderingh5h5BetterILUdecompositionbysmartreordering"></a><h5>Better ILU decomposition by smart reordering</h5><h5>Better ILU decomposition by smart reordering</h5>

</a> 提高线性求解过程速度的第一个尝试是选择一个道夫重排序，使ILU更接近于全LU分解，这在代码中的注释中已经提到。DoFRenumbering命名空间比较了几种对斯托克斯方程的道夫重新编号的选择。关于计算时间的最佳结果是通过调用 DoFRenumbering::boost::king_ordering. 找到的，通过该程序，内部求解器需要的操作大大减少，例如，在第4周期对 $A$ 进行反演的内部CG迭代约62次，而标准Cuthill-McKee-算法的迭代约75次。另外，在第4周期，对于 <code>solve()</code> 的调用，计算时间从大约17分钟减少到11分钟。然而，King排序（以及一般由 DoFRenumbering::boost 命名空间提供的排序）有一个严重的缺点--它比构建中的交易版本使用更多的内存，因为它作用于抽象图而不是由三角化提供的几何图形。在目前的情况下，重新编号需要大约5倍的内存，这就产生了一个不可行的算法，在3D的最后一个周期有120万未知数。

<a name="BetterpreconditionerfortheinnerCGsolver"></a><h5>Better preconditioner for the inner CG solver</h5>另一个改善情况的想法是选择一个预处理程序，使(0,0)矩阵 $A$ 的CG在一个与网格无关的迭代次数中收敛，比如10到30次。我们在步骤16中已经看到了这样的候选方案：多网格。


<a name="BlockSchurcomplementpreconditioner"></a><h5>Block Schur complement preconditioner</h5> <a name="block-schur"></a> 即使对 $A$ 有一个好的预处理程序，我们仍然需要反复求解同一个线性系统（虽然右手边不同），以使舒尔补码的求解趋于一致。我们这里要讨论的方法是如何将内部迭代和外部迭代结合起来。如果我们坚持计算舒尔补码，就没有其他的可能性了。


另一种方法是一次性攻击块状系统，并使用近似的Schur补码作为有效的预处理程序。其思路如下。如果我们找到一个块状预处理器 $P$ ，使矩阵

@f{eqnarray*}
  P^{-1}\left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)


@f}

是简单的，那么使用该预处理程序的迭代求解器将在几次迭代中收敛。使用舒尔补码 $S = B A^{-1} B^T$ ，我们发现

@f{eqnarray*}
  P^{-1}
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
  \end{array}\right)


@f}

似乎是个不错的选择，因为

@f{eqnarray*}
  P^{-1}\left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    A^{-1} & 0 \\ S^{-1} B A^{-1} & -S^{-1}
  \end{array}\right)\cdot \left(\begin{array}{cc}
    A & B^T \\ B & 0
  \end{array}\right)
  =
  \left(\begin{array}{cc}
    I & A^{-1} B^T \\ 0 & I
  \end{array}\right).


@f}

这就是引言中提到的Silvester和Wathen的论文所采取的方法（不同的是Silvester和Wathen使用了右预处理）。在这种情况下，基于Krylov的迭代方法只有在应用 $A$ 和 $S$ 的精确求逆时才会一步收敛，因为所有的特征值都是1（而这种方法的迭代次数是由不同特征值的数量所决定的）。下面，我们将讨论如何为这个问题选择一个适当的求解器。首先，我们要仔细研究一下预处理程序的实现。

由于 $P$ 的目的只是作为一个预处理程序，我们将使用舒尔补码 $S$ 和矩阵 $A$ 的近似反向。因此，舒尔补集将由压力质量矩阵 $M_p$ 近似，我们使用 $A$ 的预处理器（周围没有反矩阵类）来近似 $A^{-1}$ 。

这里有一个实现块舒尔补码预处理的类。根据上面的推导，对块向量的 <code>vmult</code> 操作可以通过三个连续的操作来指定。

@code
template <class PreconditionerA, class PreconditionerMp>
class BlockSchurPreconditioner : public Subscriptor
{
  public:
    BlockSchurPreconditioner (const BlockSparseMatrix<double>         &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp>  &Mpinv,
          const PreconditionerA &Apreconditioner);


  void vmult (BlockVector<double>       &dst,
              const BlockVector<double> &src) const;


  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,
                       PreconditionerMp > > m_inverse;
    const PreconditionerA &a_preconditioner;


    mutable Vector<double> tmp;


};


template <class PreconditionerA, class PreconditionerMp>
BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::BlockSchurPreconditioner(
          const BlockSparseMatrix<double>                            &S,
          const InverseMatrix<SparseMatrix<double>,PreconditionerMp> &Mpinv,
          const PreconditionerA &Apreconditioner
          )
                :
                system_matrix           (&S),
                m_inverse               (&Mpinv),
                a_preconditioner        (Apreconditioner),
                tmp                     (S.block(1,1).m())
{}


        // Now the interesting function, the multiplication of
        // the preconditioner with a BlockVector.
template <class PreconditionerA, class PreconditionerMp>
void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
                                     BlockVector<double>       &dst,
                                     const BlockVector<double> &src) const
{
        // Form u_new = A^{-1} u
  a_preconditioner.vmult (dst.block(0), src.block(0));
        // Form tmp = - B u_new + p
        // (<code>SparseMatrix::residual</code>
        // does precisely this)
  system_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
        // Change sign in tmp
  tmp *= -1;
        // Multiply by approximate Schur complement
        // (i.e. a pressure mass matrix)
  m_inverse->vmult (dst.block(1), tmp);
}
@endcode



由于我们现在对整个区块系统采取行动，我们必须忍受一个缺点：我们需要对整个区块系统而不是较小的压力空间进行求解器迭代。

现在我们转向我们应该对块系统使用哪种求解器的问题。第一个观察结果是，所产生的预处理矩阵不能用CG求解，因为它既不是正定也不是对称的。

deal.II库实现了几个适合手头问题的求解器。一种选择是求解器 @ref SolverBicgstab "BiCGStab"，它被用于解决步骤9中的非对称平流问题。第二个选择，也就是我们要选择的，是 @ref SolverGMRES  "GMRES"（广义最小残差）。这两种方法都有其优点和缺点--在有些问题上，两种候选方法中的一种明显优于另一种，反之亦然。<a href="http://en.wikipedia.org/wiki/GMRES#Comparison_with_other_solvers">Wikipedia</a>关于GMRES方法的文章给出了一个比较的介绍。更全面和有根据的比较可以在J.W.Demmel的书中读到（Applied Numerical Linear Algebra, SIAM, 1997, section 6.6.6）。

对于我们用ILU预处理 $A$ 的具体问题，对于大的问题规模，我们当然需要在块系统上进行数百次迭代（我们不会打败CG！）。实际上，这不利于GMRES。在GMRES迭代过程中，Krylov向量的基础被陆续建立起来，并对这些向量进行一些操作。这个基础上的向量越多，需要的操作和内存就越多。操作的数量以 ${\cal O}(n + k^2)$ 的形式扩展，内存以 ${\cal O}(kn)$ 的形式扩展，其中 $k$ 是Krylov基础中的向量数量， $n$ 是（块）矩阵的大小。为了不让这些需求过度增长，deal.II将基的大小 $k$ 默认限制为30个向量。然后，重新建立基。这种GMRES方法的实现被称为GMRES(k)，默认为  $k=30$  。我们通过这一限制所获得的东西，即对操作和内存需求的约束，将被我们使用不完整的基础这一事实所补偿--这将增加所需的迭代次数。

另一方面，当需要多次迭代时，BiCGStab不会变慢（一次迭代只使用前面一个步骤的结果，而不是像GMRES那样使用所有的步骤）。除了BiCGStab由于需要两个矩阵-向量乘积（相比之下，CG或GMRES只需要一个），所以每一步的成本更高之外，还有一个主要原因使得BiCGStab不适合这个问题：预处理程序通过使用InverseMatrix类应用压力质量矩阵的逆。由于向量的逆矩阵应用只是以近似的方式进行（精确的逆太昂贵了），这也会影响求解器。在BiCGStab的情况下，由于这种扰动，Krylov向量将不会是正交的。虽然这对于少量的步骤（最多50步）来说是不关键的，但当这些扰动在迭代的粗放中增长到相当大的程度时，它就会破坏求解器的性能。

我们用BiCGStab做了一些实验，发现它在细化周期3之前比GMRES快（在3D中），但在周期4和5时变得非常慢（甚至比原来的Schur补码还慢），所以在这种情况下求解器是没有用的。为逆矩阵类选择一个更尖锐的容忍度（ <code>1e-10*src.l2_norm()</code> 而不是 <code>1e-6*src.l2_norm()</code> ）使BiCGStab在第4周期也表现良好，但没有改变在非常大的问题上的失败。

当然，GMRES也会受到近似求逆的影响，但它对正交性不那么敏感，而且对于大尺寸也能保持相对较好的性能，见下面的结果。

说到这里，我们转向用 $k=100$ 临时向量的GMRES实现求解器调用。

@code
      const SparseMatrix<double> &pressure_mass_matrix
        = preconditioner_matrix.block(1,1);
      SparseILU<double> pmass_preconditioner;
      pmass_preconditioner.initialize (pressure_mass_matrix,
        SparseILU<double>::AdditionalData());


      InverseMatrix<SparseMatrix<double>,SparseILU<double> >
        m_inverse (pressure_mass_matrix, pmass_preconditioner);


      BlockSchurPreconditioner<typename InnerPreconditioner<dim>::type,
                               SparseILU<double> >
        preconditioner (system_matrix, m_inverse, *A_preconditioner);


      SolverControl solver_control (system_matrix.m(),
                                    1e-6*system_rhs.l2_norm());
      GrowingVectorMemory<BlockVector<double> > vector_memory;
      SolverGMRES<BlockVector<double> >::AdditionalData gmres_data;
      gmres_data.max_n_tmp_vectors = 100;


      SolverGMRES<BlockVector<double> > gmres(solver_control, vector_memory,
                                              gmres_data);


      gmres.solve(system_matrix, solution, system_rhs,
                  preconditioner);


      constraints.distribute (solution);


      std::cout << " "
                << solver_control.last_step()
                << " block GMRES iterations";
@endcode



显然，人们需要添加include文件 @ref SolverGMRES "<lac/solver_gmres.h>"以使其运行。我们用BlockVector模板来调用求解器，以便使GMRES能够对块状向量和矩阵进行操作。还要注意的是，在我们将信息复制到另一个矩阵之后，我们需要将系统矩阵中的（1,1）块设置为零（我们将压力质量矩阵保存在那里，这不是问题的一部分）。

使用定时器类，我们收集了一些统计数据，将块状求解器的运行时间与上述问题实现中的运行时间进行比较。除了两个选项的解决方案，我们还检查了两个变体的解决方案是否接近（即这个求解器给出的解决方案确实与我们之前的解决方案相同），并计算矢量差的无穷大准则。

让我们先看看二维的结果。

@code
Refinement cycle 0
   Number of active cells: 64
   Number of degrees of freedom: 679 (594+85) [0.00162792 s]
   Assembling...  [0.00108981 s]
   Computing preconditioner... [0.0025959 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.00479603s ]
      Block Schur preconditioner: 12 GMRES iterations [0.00441718 s]
   l_infinity difference between solution vectors: 5.38258e-07


Refinement cycle 1
   Number of active cells: 160
   Number of degrees of freedom: 1683 (1482+201) [0.00345707 s]
   Assembling...  [0.00237417 s]
   Computing preconditioner... [0.00605702 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.0123992s ]
      Block Schur preconditioner: 12 GMRES iterations [0.011909 s]
   l_infinity difference between solution vectors: 1.74658e-05


Refinement cycle 2
   Number of active cells: 376
   Number of degrees of freedom: 3813 (3370+443) [0.00729299 s]
   Assembling...  [0.00529909 s]
   Computing preconditioner... [0.0167508 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.031672s ]
      Block Schur preconditioner: 12 GMRES iterations [0.029232 s]
   l_infinity difference between solution vectors: 7.81569e-06


Refinement cycle 3
   Number of active cells: 880
   Number of degrees of freedom: 8723 (7722+1001) [0.017709 s]
   Assembling...  [0.0126002 s]
   Computing preconditioner... [0.0435679 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.0971651s ]
      Block Schur preconditioner: 12 GMRES iterations [0.0992041 s]
   l_infinity difference between solution vectors: 1.87249e-05


Refinement cycle 4
   Number of active cells: 2008
   Number of degrees of freedom: 19383 (17186+2197) [0.039988 s]
   Assembling...  [0.028281 s]
   Computing preconditioner... [0.118314 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.252133s ]
      Block Schur preconditioner: 13 GMRES iterations [0.269125 s]
   l_infinity difference between solution vectors: 6.38657e-05


Refinement cycle 5
   Number of active cells: 4288
   Number of degrees of freedom: 40855 (36250+4605) [0.0880702 s]
   Assembling...  [0.0603511 s]
   Computing preconditioner... [0.278339 s]
   Solving...
      Schur complement: 11 outer CG iterations for p  [0.53846s ]
      Block Schur preconditioner: 13 GMRES iterations [0.578667 s]
   l_infinity difference between solution vectors: 0.000173363
@endcode



我们看到，块状舒尔补码预处理求解器和舒尔补码本身在求解时间上没有巨大差异。原因很简单：我们使用直接求解作为 $A$ 的预处理程序--所以我们不能指望通过避免内部迭代获得任何收益。我们看到，GMRES的迭代次数略有增加，但总的来说，这两种选择是相当相似的。

画面当然会发生三维变化。

@code
Refinement cycle 0
   Number of active cells: 32
   Number of degrees of freedom: 1356 (1275+81) [0.00845218 s]
   Assembling...  [0.019372 s]
   Computing preconditioner... [0.00712395 s]
   Solving...
      Schur complement: 13 outer CG iterations for p  [0.0320101s ]
      Block Schur preconditioner: 22 GMRES iterations [0.0048759 s]
   l_infinity difference between solution vectors: 2.15942e-05


Refinement cycle 1
   Number of active cells: 144
   Number of degrees of freedom: 5088 (4827+261) [0.0346942 s]
   Assembling...  [0.0857739 s]
   Computing preconditioner... [0.0465031 s]
   Solving...
      Schur complement: 14 outer CG iterations for p  [0.349258s ]
      Block Schur preconditioner: 35 GMRES iterations [0.048759 s]
   l_infinity difference between solution vectors: 1.77657e-05


Refinement cycle 2
   Number of active cells: 704
   Number of degrees of freedom: 22406 (21351+1055) [0.175669 s]
   Assembling...  [0.437447 s]
   Computing preconditioner... [0.286435 s]
   Solving...
      Schur complement: 14 outer CG iterations for p  [3.65519s ]
      Block Schur preconditioner: 63 GMRES iterations [0.497787 s]
   l_infinity difference between solution vectors: 5.08078e-05


Refinement cycle 3
   Number of active cells: 3168
   Number of degrees of freedom: 93176 (89043+4133) [0.790985 s]
   Assembling...  [1.97598 s]
   Computing preconditioner... [1.4325 s]
   Solving...
      Schur complement: 15 outer CG iterations for p  [29.9666s ]
      Block Schur preconditioner: 128 GMRES iterations [5.02645 s]
   l_infinity difference between solution vectors: 0.000119671


Refinement cycle 4
   Number of active cells: 11456
   Number of degrees of freedom: 327808 (313659+14149) [3.44995 s]
   Assembling...  [7.54772 s]
   Computing preconditioner... [5.46306 s]
   Solving...
      Schur complement: 15 outer CG iterations for p  [139.987s ]
      Block Schur preconditioner: 255 GMRES iterations [38.0946 s]
   l_infinity difference between solution vectors: 0.00020793


Refinement cycle 5
   Number of active cells: 45056
   Number of degrees of freedom: 1254464 (1201371+53093) [19.6795 s]
   Assembling...  [28.6586 s]
   Computing preconditioner... [22.401 s]
   Solving...
      Schur complement: 14 outer CG iterations for p  [796.767s ]
      Block Schur preconditioner: 524 GMRES iterations [355.597 s]
   l_infinity difference between solution vectors: 0.000501219
@endcode



在这里，块状预处理求解器明显优于Schur补数，但是网格点越多，优势就越小。这是因为GMRES(k)随着问题规模的扩大比CG更差，正如我们上面所讨论的。  尽管如此，对于中等规模的问题，3-6倍的改进是相当令人印象深刻的。




<a name="Combiningtheblockpreconditionerandmultigrid"></a><h5>Combining the block preconditioner and multigrid</h5>这个问题的终极线性求解器可以想象为 $A$ 的最佳预处理器（如多网格）和上述的块状预处理器的组合，这就是步骤31和步骤32教程程序（我们使用代数多网格方法）和步骤56（我们使用几何多网格方法）中所采取的方法。





<a name="Noblockmatricesandvectors"></a><h5>No block matrices and vectors</h5> 另一个可以考虑的可能性是不设置块状系统，而是一次性解决速度和压力系统。可以选择用UMFPACK直接求解（2D）或用ILU预处理的GMRES（3D）。这应该是很直接的尝试。





<a name="Moreinterestingtestcases"></a><h4>More interesting testcases</h4>


当然，这个程序也可以作为计算更有趣的情况下的流动的基础。编写这个程序的最初动机是希望它能成为一些地球物理流动问题的起点，例如大陆板块漂移分离的地方（例如洋中脊）下的岩浆运动。当然，在这种地方，几何形状比上面的例子更复杂，但要适应这种情况并不难。

例如，通过使用以下对边界值函数的修改

@code
template <int dim>
double
BoundaryValues<dim>::value (const Point<dim>  &p,
                            const unsigned int component) const
{
  Assert (component < this->n_components,
          ExcIndexRange (component, 0, this->n_components));


  const double x_offset = std::atan(p[1]*4)/3;


  if (component == 0)
    return (p[0] < x_offset ? -1 : (p[0] > x_offset ? 1 : 0));
  return 0;
}
@endcode

和以下方式生成网格作为域 $[-2,2]\times[-2,2]\times[-1,0]$ 。

@code
    std::vector<unsigned int> subdivisions (dim, 1);
    subdivisions[0] = 4;
    if (dim>2)
      subdivisions[1] = 4;


    const Point<dim> bottom_left = (dim == 2 ?
                                    Point<dim>(-2,-1) :
                                    Point<dim>(-2,-2,-1));
    const Point<dim> top_right   = (dim == 2 ?
                                    Point<dim>(2,0) :
                                    Point<dim>(2,2,0));


    GridGenerator::subdivided_hyper_rectangle (triangulation,
                                               subdivisions,
                                               bottom_left,
                                               top_right);
@endcode

那么我们就会得到断层线是弯曲的图像。   <table width="60%" align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d-extension.png" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-22.3d-grid-extension.png" alt="">
    </td>
  </tr>
</table> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-22.cc"
*/
