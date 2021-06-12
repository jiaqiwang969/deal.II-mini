/**
@page step_25 The step-25 tutorial program
This tutorial depends on step-24.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#Statementoftheproblem">Statement of the problem</a>
        <li><a href="#Weakformulationofthetimediscretizedequations">Weak formulation of the time-discretized equations</a>
        <li><a href="#Discretizationoftheweakformulationinspace">Discretization of the weak formulation in space</a>
        <li><a href="#Thetestcase">The test case</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefilesandglobalvariables">Include files and global variables</a>
        <li><a href="#ThecodeSineGordonProblemcodeclasstemplate">The <code>SineGordonProblem</code> class template</a>
        <li><a href="#Initialconditions">Initial conditions</a>
        <li><a href="#ImplementationofthecodeSineGordonProblemcodeclass">Implementation of the <code>SineGordonProblem</code> class</a>
      <ul>
        <li><a href="#SineGordonProblemSineGordonProblem">SineGordonProblem::SineGordonProblem</a>
        <li><a href="#SineGordonProblemmake_grid_and_dofs">SineGordonProblem::make_grid_and_dofs</a>
        <li><a href="#SineGordonProblemassemble_system">SineGordonProblem::assemble_system</a>
        <li><a href="#SineGordonProblemcompute_nl_term">SineGordonProblem::compute_nl_term</a>
        <li><a href="#SineGordonProblemcompute_nl_matrix">SineGordonProblem::compute_nl_matrix</a>
        <li><a href="#SineGordonProblemsolve">SineGordonProblem::solve</a>
        <li><a href="#SineGordonProblemoutput_results">SineGordonProblem::output_results</a>
        <li><a href="#SineGordonProblemrun">SineGordonProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#An11dSolution">An (1+1)-d Solution</a>
        <li><a href="#Afew21DSolutions">A few (2+1)D Solutions</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-25/doc/intro.dox

<a name="Introduction"></a><a name="Intro"></a> <h1>Introduction</h1>。


这个程序是由德克萨斯A&amp;M大学的Ivan Christov的一个学生项目发展而来的。这个程序的大部分工作都是由他完成的。

这个程序的目标是解决1、2或3空间维度的正弦戈登孤子方程。解决这个方程的动机是对二维和三维解的性质知之甚少，尽管一维情况已经被广泛研究。

确切地说，正弦-戈登方程的名称是对所谓的克莱因-戈登方程的双关语，它是薛定谔方程的相对论版本，适用于非零质量的粒子。这种相似性不仅仅是表面的，正弦-戈登方程已经被证明可以模拟一些统一场现象，如亚原子粒子的相互作用（见，例如，Perring &amp; Skyrme in Nuclear %Physics <b>31</b>）和超导体结中的约瑟夫森（量子）效应（见，例如<a
href="http://en.wikipedia.org/wiki/Long_Josephson_junction">http://en.wikipedia.org/wiki/Long_Josephson_junction</a>）。此外，从数学的角度来看，由于正弦戈登方程是 "完全可积分的"，它是使用反散射变换等通常方法研究的候选者。因此，多年来，人们发现了许多有趣的孤波，甚至是静止的正弦戈登方程解。在这些解中，粒子对应于局部特征。关于正弦戈登方程、反散射变换和其他寻找分析性孤子方程的方法的更多信息，读者应参考以下关于该主题的 "经典 "参考资料。G. L. Lamb的<i>Elements of Soliton
Theory</i>（第5章，第2节）和G. B. Whitham的<i>Linear and
Nonlinear Waves</i>（第17章，10-13节）。

 @note  我们将在第58步中介绍量子力学中一个单独的非线性方程，即非线性Schr&ouml;dinger方程。

<a name="Statementoftheproblem"></a><h3>Statement of the problem</h3> 我们希望解决的正弦戈登初始边界值问题（IBVP）由以下方程组成。


\f{eqnarray*}
  u_{tt}-\Delta u &=& -\sin(u) \quad\mbox{for}\quad (x,t) \in \Omega \times (t_0,t_f],\\
  {\mathbf n} \cdot \nabla u &=& 0 \quad\mbox{for}\quad (x,t) \in \partial\Omega
           \times (t_0,t_f],\\
  u(x,t_0) &=& u_0(x).
\f} 这是一个非线性方程，类似于我们在步骤23和步骤24中讨论的波浪方程。我们选择执行零诺伊曼边界条件，以使波从我们的域的边界上反射出去。然而，应该注意的是，迪里希特边界条件对这个问题并不合适。尽管正弦戈登方程的解是局部的，但只有在 $x=\pm\infty$ 处指定（迪里希特）边界条件才有意义，否则要么不存在解，要么只存在微不足道的解 $u=0$ 。

然而，上述方程的形式对于数值离散化来说并不理想。如果我们要直接准确地离散二阶时间导数，那么我们将需要一个大的模板（即需要在内存中保留几个时间步长），这可能变得很昂贵。因此，与我们在步骤23和步骤24中所做的完全类似，我们将二阶（时间上）正弦-戈登方程拆分为两个一阶（时间上）方程系统，我们称之为拆分，或速度公式。为此，通过设置 $v = u_t$ ，很容易看到正弦-戈登方程等同于

\f{eqnarray*}
  u_t - v &=& 0,\\
  v_t - \Delta u &=& -\sin(u).
\f}

现在，我们可以使用 $\theta$ 方法对分裂公式进行时间离散，该方法的模板只有两个时间步长。通过选择 $\theta\in [0,1]$ ，后者的离散化允许我们从一系列的方案中进行选择。特别是，如果我们选择 $\theta=0$ 或 $\theta=1$ ，我们可以分别得到一阶精确的显式或隐式欧拉方法。另一个重要的选择是 $\theta=\frac{1}{2}$ ，它给出了二阶精确的Crank-Nicolson方案。因此，上标 $n$ 表示在 $n^{\mathrm{th}}$ 时间步长的变量值，即在 $t=t_n \dealcoloneq n k$ ，其中 $k$ 是（固定）时间步长。因此，时间分解的正弦-戈登方程的拆分表述为

\f{eqnarray*}
  \frac{u^n - u^{n-1}}{k} - \left[\theta v^n + (1-\theta) v^{n-1}\right] &=& 0,\\
  \frac{v^n - v^{n-1}}{k} - \Delta\left[\theta u^n + (1-\theta) u^{n-1}\right]
  &=& -\sin\left[\theta u^n + (1-\theta) u^{n-1}\right].
\f}

我们可以通过一点代数来简化后者。从第一个方程中排除 $v^n$ 并重新排列，我们可以得到

\f{eqnarray*}
  \left[ 1-k^2\theta^2\Delta \right] u^n &=&
         \left[ 1+k^2\theta(1-\theta)\Delta\right] u^{n-1} + k v^{n-1}


         - k^2\theta\sin\left[\theta u^n + (1-\theta) u^{n-1}\right],\\
   v^n &=& v^{n-1} + k\Delta\left[ \theta u^n + (1-\theta) u^{n-1}\right]


         - k\sin\left[ \theta u^n + (1-\theta) u^{n-1} \right].
\f}

在这一点上，我们似乎可以直接进行空间离散化的方程。虽然这对第二个方程（在 $v^n$ 中是线性的）来说是正确的，但这对所有的 $\theta$ 来说是不可行的，因为上面的第一个方程是非线性的。因此，必须实现一个非线性求解器，然后将方程在空间中离散化并求解。

为此，我们可以使用牛顿方法。给定非线性方程 $F(u^n) = 0$ ，我们对 $u^n$ 产生连续的近似值，如下。

\f{eqnarray*}
  \mbox{ Find } \delta u^n_l \mbox{ s.t. } F'(u^n_l)\delta u^n_l = -F(u^n_l)
  \mbox{, set }  u^n_{l+1} = u^n_l + \delta u^n_l.
\f} 迭代可以用旧的时间步长进行初始化，即 $u^n_0 = u^{n-1}$  ，最终会产生分裂公式的第一个方程的解（见上文）。对于这里所考虑的正弦-戈登方程的时间离散化，我们有

\f{eqnarray*}
  F(u^n_l) &=&  \left[ 1-k^2\theta^2\Delta \right] u^n_l -
                 \left[ 1+k^2\theta(1-\theta)\Delta\right] u^{n-1} - k v^{n-1}
                 + k^2\theta\sin\left[\theta u^n_l + (1-\theta) u^{n-1}\right],\\
  F'(u^n_l) &=& 1-k^2\theta^2\Delta + k^2\theta^2\cos\left[\theta u^n_l
                        + (1-\theta) u^{n-1}\right].
\f} 注意，虽然 $F(u^n_l)$ 是一个函数，但 $F'(u^n_l)$ 是一个运算符。

<a name="Weakformulationofthetimediscretizedequations"></a><h3>Weak formulation of the time-discretized equations</h3> 事后看来，我们选择解和检验空间都是 $H^1(\Omega)$  。因此，乘以测试函数 $\varphi$ 并进行积分，我们得到每个时间步长的分裂公式（包括第一个方程的非线性求解器）的以下变分（或弱）公式。


\f{eqnarray*}
  &\mbox{ Find}& \delta u^n_l \in H^1(\Omega) \mbox{ s.t. }
  \left( F'(u^n_l)\delta u^n_l, \varphi \right)_{\Omega}
  = -\left(F(u^n_l), \varphi \right)_{\Omega} \;\forall\varphi\in H^1(\Omega),
  \mbox{ set } u^n_{l+1} = u^n_l + \delta u^n_l,\; u^n_0 = u^{n-1}.\\
  &\mbox{ Find}& v^n \in H^1(\Omega) \mbox{ s.t. }
  \left( v^n, \varphi \right)_{\Omega} = \left( v^{n-1}, \varphi \right)_{\Omega}


         - k\theta\left( \nabla u^n, \nabla\varphi \right)_{\Omega}


         - k (1-\theta)\left( \nabla u^{n-1}, \nabla\varphi \right)_{\Omega}


         - k\left(\sin\left[ \theta u^n + (1-\theta) u^{n-1} \right],
         \varphi \right)_{\Omega} \;\forall\varphi\in H^1(\Omega).
\f}注意，我们在涉及拉普拉斯算子的所有项上使用了分项积分和零诺伊曼边界条件。此外， $F(\cdot)$ 和 $F'(\cdot)$ 如上定义， $(\cdot,\cdot)_{\Omega}$ 表示域 $\Omega$ 上通常的 $L^2$ 内积，即 $(f,g)_{\Omega} = \int_\Omega fg
\,\mathrm{d}x$  。最后，请注意，第一个方程实际上是一个迭代程序的定义，所以在每个时间步骤中，它被多次解决，直到满足停止标准。

<a name="Discretizationoftheweakformulationinspace"></a><h3>Discretization of the weak formulation in space</h3>使用有限元方法，我们在空间中离散变量公式。为此，让 $V_h$ 成为一个具有节点基 $\{\varphi_1,\ldots,\varphi_N\}$ 的有限元空间（ $\mathrm{dim}\, V_h = N

< \infty$ ）。现在，我们可以用节点基来展开弱公式（见上文）中的所有函数。此后，我们将用大写字母表示一个函数的系数向量（在节点基中），用小写字母表示；例如， $u^n = \sum_{i=1}^N
U^n_i \varphi_i$  其中 $U^n \in {R}^N$  和 $u^n \in
H^1(\Omega)$  。因此，变量公式的有限维版本要求我们在每个时间步长中解决以下矩阵方程。

@f{eqnarray*}
  F_h'(U^{n,l})\delta U^{n,l} &=& -F_h(U^{n,l}), \qquad
        U^{n,l+1} = U^{n,l} + \delta U^{n,l}, \qquad U^{n,0} = U^{n-1}; \\
  MV^n &=& MV^{n-1} - k \theta AU^n -k (1-\theta) AU^{n-1} - k S(u^n,u^{n-1}).


@f}

以上，矩阵 $F_h'(\cdot)$ 和向量 $F_h(\cdot)$ 表示上面讨论的小工具的离散版本，即。

\f{eqnarray*}
  F_h(U^{n,l}) &=&  \left[ M+k^2\theta^2A \right] U^{n,l} -
                \left[ M-k^2\theta(1-\theta)A \right] U^{n-1} - k MV^{n-1}
                + k^2\theta S(u^n_l, u^{n-1}),\\
  F_h'(U^{n,l}) &=& M+k^2\theta^2A
                                + k^2\theta^2N(u^n_l,u^{n-1})
\f} 再次注意，上面的第一个矩阵方程实际上是迭代程序的定义，所以它被多次求解，直到满足停止标准。此外， $M$ 是质量矩阵，即 $M_{ij} = \left( \varphi_i,\varphi_j \right)_{\Omega}$  ， $A$ 是拉普拉斯矩阵，即 $A_{ij} = \left( \nabla \varphi_i, \nabla
\varphi_j \right)_{\Omega}$  ， $S$  是定义我们的辅助速度变量的方程中的非线性项，即 $S_j(f,g) = \left(
  \sin\left[ \theta f + (1-\theta) g\right], \varphi_j \right)_{\Omega}$  ， $N$  是 $F(\cdot)$  的雅各布矩阵中的非线性项，即 $N_{ij}(f,g) = \left( \cos\left[ \theta f + (1-\theta) g\right]\varphi_i,
  \varphi_j \right)_{\Omega}$  。

对于第一个方程，我们可以用什么解法？让我们来看看我们要反转的矩阵。

@f[
  (M+k^2\theta^2(A + N))_{ij} =
  \int_\Omega (1+k^2\theta^2 \cos \alpha)
  \varphi_i\varphi_j \; dx
  + k^2 \theta^2 \int_\Omega \nabla\varphi_i\nabla\varphi_j \; dx,


@f]

为一些 $\alpha$ ，取决于现在和以前的解决方案。首先，请注意，该矩阵是对称的。此外，如果时间步长 $k$ 足够小，即如果 $k\theta<1$ ，那么该矩阵也将是正定的。在下面的程序中，情况总是这样的，所以我们将使用共轭梯度法和SSOR法作为预处理。然而，我们应该记住，如果我们碰巧使用更大的时间步长，这将失败。幸运的是，在这种情况下，求解器将只是抛出一个异常，表明收敛失败，而不是默默地产生一个错误的结果。如果发生这种情况，那么我们可以简单地用能够处理不确定对称系统的方法取代CG方法。GMRES求解器通常是处理所有 "坏 "线性系统的标准方法，但它也是一个缓慢的方法。更好的方法可能是利用对称性的求解器，例如，SymmLQ，它也在deal.II中实现。

这个程序在步骤23和 @ref
step_24 "步骤24 "上使用了巧妙的优化。如果你仔细阅读上述公式，就会发现，速度 $V$ 只出现在与质量矩阵的乘积中。因此，在步骤23和步骤24中，我们有点浪费：在每个时间步骤中，我们会用质量矩阵求解一个线性系统，只是在下一个时间步骤中再次将该系统的解乘以 $M$ 。当然，这可以避免，我们在这个程序中也是这样做的。




<a name="Thetestcase"></a><h3>The test case</h3>


正弦戈登方程有几个分析解，包括一维和二维的。特别是，该程序如是计算一个具有单一的类似激波的初始条件的问题的解。  这个解是由Leibbrandt在Phys.Rev.Lett.中给出的。\b 41(7)中给出，并在 <code>ExactSolution</code> 类中实现。

应该注意的是，严格来说，这个闭式解只适用于无限空间的初值问题（而不是这里考虑的诺伊曼初界值问题）。然而，鉴于我们施加了零诺依曼边界条件，我们期望我们的初始边界值问题的解将接近无限空间初始值问题的解，如果波在我们领域的边界上的反射没有发生。在实践中，情况当然不是这样的，但我们至少可以假设是这样的。

二维解决方案中的常数 $\vartheta$ 和 $\lambda$ 以及三维解决方案中的 $\vartheta$ 、 $\phi$ 和 $\tau$ 被称为B&auml;cklund变换参数。它们控制诸如扭结的方向和陡度。为了测试代码与精确解的对比，我们应该选择这些参数，使扭结与网格保持一致。

我们在 <code>ExactSolution</code> 类中实现的解决方案是这些。   <ul>   <li>  在1D中：@f[
  u(x,t) =


  -4 \arctan\left[
     \frac{m}{\sqrt{1-m^2}}
     \frac{\sin\left(\sqrt{1-m^2}t+c_2\right)}
     {\cosh\left(mx+c_1\right)}
     \right],
  @f] 。

  其中我们选择  $m=\frac 12, c_1=c_2=0$  。

  在一维中，已知有更多有趣的分析解决方案。他们中的许多人被列在http://mathworld.wolfram.com/Sine-GordonEquation.html。

    <li>  在2D：@f[
    u(x,y,t) = 4 \arctan \left[a_0 e^{s\xi}\right],
  @f] 。

  其中 $\xi$ 被定义为@f[
    \xi = x \cos\vartheta + \sin(\vartheta) (y\cosh\lambda + t\sinh \lambda),
  @f]。

  而在这里我们选择 $\vartheta=\frac \pi 4, \lambda=a_0=s=1$ 。

    <li>  在三维：@f[
    u(x,y,z,t) = 4 \arctan \left[c_0 e^{s\xi}\right],
  @f] 。

  其中 $\xi$ 被定义为@f[
    \xi = x \cos\vartheta + y \sin \vartheta \cos\phi +
          \sin \vartheta \sin\phi (z\cosh\tau + t\sinh \tau),
  @f]。

  而在这里，我们选择 $\vartheta=\phi=\frac{\pi}{4}, \tau=c_1=s=1$  。   </ul> 


由于这使得玩起来更容易，用于设置&mdash;惊喜！&mdash;我们模拟的初始值的 <code>InitialValues</code> 类只是查询了描述初始时值的精确解的类，而不是重复实现解函数的努力。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefilesandglobalvariables"></a> 
 * <h3>Include files and global variables</h3>
 * 

 * 
 * 关于包含文件的解释，读者应该参考示例程序  step-1  到  step-4  。它们的标准顺序是  <code>base</code> -- <code>lac</code> -- <code>grid</code>  --  <code>dofs</code> -- <code>fe</code> -- <code>numerics</code>  （因为每一类大致都是建立在前面的基础上），然后是一些用于文件输入/输出和字符串流的C++头文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/utilities.h> 
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
 * 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后一步和以前所有的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step25 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeSineGordonProblemcodeclasstemplate"></a> 
 * <h3>The <code>SineGordonProblem</code> class template</h3>
 * 

 * 
 * 解决问题的整个算法被封装在这个类中。和以前的例子程序一样，这个类在声明时有一个模板参数，就是空间维度，这样我们就可以在一个、两个或三个空间维度上解决正弦-戈登方程。关于这个问题的独立于维度的类封装的更多信息，读者应该参考  step-3  和  step-4  。
 * 

 * 
 * 与 step-23 和 step-24 相比，在程序的总体结构中没有任何值得注意的地方（当然，在各种函数的内部运作中也有！）。最明显的区别是出现了两个新的函数 <code>compute_nl_term</code> 和 <code>compute_nl_matrix</code> ，计算系统矩阵的非线性贡献和第一个方程的右手边，正如在介绍中讨论的那样。此外，我们还必须有一个向量 <code>solution_update</code> ，它包含在每个牛顿步骤中对解向量的非线性更新。
 * 

 * 
 * 正如介绍中也提到的，我们在这个程序中不存储速度变量，而是质量矩阵乘以速度。这是在 <code>M_x_velocity</code> 变量中完成的（"x "是代表 "次数"）。
 * 

 * 
 * 最后， <code>output_timestep_skip</code> 变量存储了在生成图形输出前每次所需的时间步数。这一点在使用精细网格（因此时间步数较小）时非常重要，在这种情况下，我们会运行大量的时间步数，并创建大量的输出文件，这些文件中的解决方案在后续文件中看起来几乎是一样的。这只会堵塞我们的可视化程序，我们应该避免创建比我们真正感兴趣的更多的输出。因此，如果这个变量被设置为大于1的值 $n$ ，那么只有在每一个 $n$ 的时间步长时才会产生输出。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class SineGordonProblem 
 *   { 
 *   public: 
 *     SineGordonProblem(); 
 *     void run(); 
 * 
 *   private: 
 *     void         make_grid_and_dofs(); 
 *     void         assemble_system(); 
 *     void         compute_nl_term(const Vector<double> &old_data, 
 *                                  const Vector<double> &new_data, 
 *                                  Vector<double> &      nl_term) const; 
 *     void         compute_nl_matrix(const Vector<double> &old_data, 
 *                                    const Vector<double> &new_data, 
 *                                    SparseMatrix<double> &nl_matrix) const; 
 *     unsigned int solve(); 
 *     void         output_results(const unsigned int timestep_number) const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     SparsityPattern      sparsity_pattern; 
 *     SparseMatrix<double> system_matrix; 
 *     SparseMatrix<double> mass_matrix; 
 *     SparseMatrix<double> laplace_matrix; 
 * 
 *     const unsigned int n_global_refinements; 
 * 
 *     double       time; 
 *     const double final_time, time_step; 
 *     const double theta; 
 * 
 *     Vector<double> solution, solution_update, old_solution; 
 *     Vector<double> M_x_velocity; 
 *     Vector<double> system_rhs; 
 * 
 *     const unsigned int output_timestep_skip; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Initialconditions"></a> 
 * <h3>Initial conditions</h3>
 * 

 * 
 * 在下面两类中，我们首先实现了本程序介绍中提到的一维、二维和三维的精确解。如果想通过比较数值解和分析解来测试程序的准确性，那么这个时空解可能会有独立的意义（但是请注意，程序使用的是有限域，而这些是无界域的分析解）。例如，这可以用 VectorTools::integrate_difference 函数来完成。再次注意（正如在 step-23 中已经讨论过的），我们如何将时空函数描述为依赖于时间变量的空间函数，该变量可以使用FunctionTime基类的 FunctionTime::set_time() 和 FunctionTime::get_time() 成员函数进行设置和查询。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class ExactSolution : public Function<dim> 
 *   { 
 *   public: 
 *     ExactSolution(const unsigned int n_components = 1, const double time = 0.) 
 *       : Function<dim>(n_components, time) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> &p, 
 *                          const unsigned int /*component*/ = 0) const override 
 *     { 
 *       const double t = this->get_time(); 
 * 
 *       switch (dim) 
 *         { 
 *           case 1: 
 *             { 
 *               const double m  = 0.5; 
 *               const double c1 = 0.; 
 *               const double c2 = 0.; 
 *               return -4. * std::atan(m / std::sqrt(1. - m * m) * 
 *                                      std::sin(std::sqrt(1. - m * m) * t + c2) / 
 *                                      std::cosh(m * p[0] + c1)); 
 *             } 
 * 
 *           case 2: 
 *             { 
 *               const double theta  = numbers::PI / 4.; 
 *               const double lambda = 1.; 
 *               const double a0     = 1.; 
 *               const double s      = 1.; 
 *               const double arg    = p[0] * std::cos(theta) + 
 *                                  std::sin(theta) * (p[1] * std::cosh(lambda) + 
 *                                                     t * std::sinh(lambda)); 
 *               return 4. * std::atan(a0 * std::exp(s * arg)); 
 *             } 
 * 
 *           case 3: 
 *             { 
 *               const double theta = numbers::PI / 4; 
 *               const double phi   = numbers::PI / 4; 
 *               const double tau   = 1.; 
 *               const double c0    = 1.; 
 *               const double s     = 1.; 
 *               const double arg   = p[0] * std::cos(theta) + 
 *                                  p[1] * std::sin(theta) * std::cos(phi) + 
 *                                  std::sin(theta) * std::sin(phi) * 
 *                                    (p[2] * std::cosh(tau) + t * std::sinh(tau)); 
 *               return 4. * std::atan(c0 * std::exp(s * arg)); 
 *             } 
 * 
 *           default: 
 *             Assert(false, ExcNotImplemented()); 
 *             return -1e8; 
 *         } 
 *     } 
 *   }; 
 * 
 * @endcode
 * 
 * 在本节的第二部分，我们提供初始条件。我们很懒惰（也很谨慎），不想第二次实现与上面相同的函数。相反，如果我们被查询到初始条件，我们会创建一个对象 <code>ExactSolution</code> ，将其设置为正确的时间，并让它计算当时的精确解的任何值。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class InitialValues : public Function<dim> 
 *   { 
 *   public: 
 *     InitialValues(const unsigned int n_components = 1, const double time = 0.) 
 *       : Function<dim>(n_components, time) 
 *     {} 
 * 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override 
 *     { 
 *       return ExactSolution<dim>(1, this->get_time()).value(p, component); 
 *     } 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeSineGordonProblemcodeclass"></a> 
 * <h3>Implementation of the <code>SineGordonProblem</code> class</h3>
 * 

 * 
 * 让我们继续讨论主类的实现，因为它实现了介绍中概述的算法。
 * 

 * 
 * 
 * <a name="SineGordonProblemSineGordonProblem"></a> 
 * <h4>SineGordonProblem::SineGordonProblem</h4>
 * 

 * 
 * 这是 <code>SineGordonProblem</code> 类的构造函数。它指定了所需的有限元的多项式程度，关联了一个 <code>DoFHandler</code> to the <code>triangulation</code> 对象（就像在示例程序 step-3 和 step-4 中一样），初始化了当前或初始时间、最终时间、时间步长，以及时间步长方案的 $\theta$ 值。由于我们在这里计算的解是时间周期性的，所以开始时间的实际值并不重要，我们选择它是为了让我们在一个有趣的时间开始。
 * 

 * 
 * 请注意，如果我们选择显式欧拉时间步进方案（ $\theta = 0$ ），那么我们必须选择一个时间步进 $k \le h$ ，否则该方案不稳定，解中可能出现振荡。Crank-Nicolson方案（ $\theta = \frac{1}{2}$ ）和隐式Euler方案（ $\theta=1$ ）不存在这个缺陷，因为它们是无条件稳定的。然而，即使如此，时间步长也应选择在 $h$ 的数量级上，以获得一个好的解决方案。由于我们知道我们的网格是由矩形的均匀细分而来，我们可以很容易地计算出这个时间步长；如果我们有一个不同的域， step-24 中的技术使用 GridTools::minimal_cell_diameter 也是可以的。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   SineGordonProblem<dim>::SineGordonProblem() 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *     , n_global_refinements(6) 
 *     , time(-5.4414) 
 *     , final_time(2.7207) 
 *     , time_step(10 * 1. / std::pow(2., 1. * n_global_refinements)) 
 *     , theta(0.5) 
 *     , output_timestep_skip(1) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemmake_grid_and_dofs"></a> 
 * <h4>SineGordonProblem::make_grid_and_dofs</h4>
 * 

 * 
 * 这个函数创建了一个 <code>dim</code> 维度的矩形网格，并对其进行了多次细化。同时，一旦自由度被集合起来， <code>SineGordonProblem</code> 类的所有矩阵和向量成员都被初始化为相应的大小。像 step-24 一样，我们使用 <code>MatrixCreator</code> 函数来生成质量矩阵 $M$ 和拉普拉斯矩阵 $A$ ，并在程序的剩余时间里将它们存储在适当的变量中。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::make_grid_and_dofs() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -10, 10); 
 *     triangulation.refine_global(n_global_refinements); 
 * 
 *     std::cout << "   Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "   Total number of cells: " << triangulation.n_cells() 
 *               << std::endl; 
 * 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 *     DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs()); 
 *     DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *     sparsity_pattern.copy_from(dsp); 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 *     mass_matrix.reinit(sparsity_pattern); 
 *     laplace_matrix.reinit(sparsity_pattern); 
 * 
 *     MatrixCreator::create_mass_matrix(dof_handler, 
 *                                       QGauss<dim>(fe.degree + 1), 
 *                                       mass_matrix); 
 *     MatrixCreator::create_laplace_matrix(dof_handler, 
 *                                          QGauss<dim>(fe.degree + 1), 
 *                                          laplace_matrix); 
 * 
 *     solution.reinit(dof_handler.n_dofs()); 
 *     solution_update.reinit(dof_handler.n_dofs()); 
 *     old_solution.reinit(dof_handler.n_dofs()); 
 *     M_x_velocity.reinit(dof_handler.n_dofs()); 
 *     system_rhs.reinit(dof_handler.n_dofs()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemassemble_system"></a> 
 * <h4>SineGordonProblem::assemble_system</h4>
 * 

 * 
 * 这个函数为牛顿方法的每次迭代组装系统矩阵和右手向量。关于系统矩阵和右手边的明确公式，读者应该参考导论。
 * 

 * 
 * 请注意，在每个时间步长中，我们必须把对矩阵和右手边的各种贡献加起来。与 step-23 和 step-24 相比，这需要集合更多的项，因为它们取决于前一个时间步骤或前一个非线性步骤的解。我们使用函数 <code>compute_nl_matrix</code> 和 <code>compute_nl_term</code> 来做到这一点，而本函数提供了顶层逻辑。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::assemble_system() 
 *   { 
 * 
 * @endcode
 * 
 * 首先我们组装雅各布矩阵 $F'_h(U^{n,l})$  ，其中 $U^{n,l}$ 为方便起见被储存在向量 <code>solution</code> 中。
 * 

 * 
 * 
 * @code
 *     system_matrix.copy_from(mass_matrix); 
 *     system_matrix.add(std::pow(time_step * theta, 2), laplace_matrix); 
 * 
 *     SparseMatrix<double> tmp_matrix(sparsity_pattern); 
 *     compute_nl_matrix(old_solution, solution, tmp_matrix); 
 *     system_matrix.add(std::pow(time_step * theta, 2), tmp_matrix); 
 * 
 * @endcode
 * 
 * 接下来我们计算右手边的向量。这只是介绍中对 $-F_h(U^{n,l})$ 的描述所暗示的矩阵-向量的组合。
 * 

 * 
 * 
 * @code
 *     system_rhs = 0.; 
 * 
 *     Vector<double> tmp_vector(solution.size()); 
 * 
 *     mass_matrix.vmult(system_rhs, solution); 
 *     laplace_matrix.vmult(tmp_vector, solution); 
 *     system_rhs.add(std::pow(time_step * theta, 2), tmp_vector); 
 * 
 *     mass_matrix.vmult(tmp_vector, old_solution); 
 *     system_rhs.add(-1.0, tmp_vector); 
 *     laplace_matrix.vmult(tmp_vector, old_solution); 
 *     system_rhs.add(std::pow(time_step, 2) * theta * (1 - theta), tmp_vector); 
 * 
 *     system_rhs.add(-time_step, M_x_velocity); 
 * 
 *     compute_nl_term(old_solution, solution, tmp_vector); 
 *     system_rhs.add(std::pow(time_step, 2) * theta, tmp_vector); 
 * 
 *     system_rhs *= -1.; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemcompute_nl_term"></a> 
 * <h4>SineGordonProblem::compute_nl_term</h4>
 * 

 * 
 * 这个函数计算向量 $S(\cdot,\cdot)$ ，它出现在分裂公式的两个方程的非线性项中。这个函数不仅简化了这个项的重复计算，而且也是我们在时间步长为隐式时使用的非线性迭代求解器的基本组成部分（即 $\theta\ne 0$  ）。此外，我们必须允许该函数接收一个 "旧 "和一个 "新 "的解决方案作为输入。这些可能不是存储在 <code>old_solution</code> and <code>solution</code> 中的问题的实际解决方案，而只是我们线性化的两个函数。为了这个函数的目的，让我们在下面这个类的文档中分别调用前两个参数  $w_{\mathrm{old}}$  和  $w_{\mathrm{new}}$  。
 * 

 * 
 * 作为一个旁注，也许值得研究一下什么阶次的正交公式最适合这种类型的积分。由于 $\sin(\cdot)$ 不是一个多项式，可能没有正交公式可以准确地积分这些项。通常只需确保右手边的积分达到与离散化方案相同的精度即可，但通过选择更精确的正交公式，也许可以改善渐近收敛声明中的常数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::compute_nl_term(const Vector<double> &old_data, 
 *                                                const Vector<double> &new_data, 
 *                                                Vector<double> &nl_term) const 
 *   { 
 *     nl_term = 0; 
 *     const QGauss<dim> quadrature_formula(fe.degree + 1); 
 *     FEValues<dim>     fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_JxW_values | 
 *                               update_quadrature_points); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     Vector<double>                       local_nl_term(dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 *     std::vector<double>                  old_data_values(n_q_points); 
 *     std::vector<double>                  new_data_values(n_q_points); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         local_nl_term = 0; 
 * 
 * @endcode
 * 
 * 一旦我们将 <code>FEValues</code> 实例化重新初始化到当前单元格，我们就利用 <code>get_function_values</code> 例程来获取 "旧 "数据（大概在 $t=t_{n-1}$ ）和 "新 "数据（大概在 $t=t_n$ ）在所选正交公式节点的值。
 * 

 * 
 * 
 * @code
 *         fe_values.reinit(cell); 
 *         fe_values.get_function_values(old_data, old_data_values); 
 *         fe_values.get_function_values(new_data, new_data_values); 
 * 
 * @endcode
 * 
 * 现在，我们可以用所需的正交公式来评估  $\int_K \sin\left[\theta w_{\mathrm{new}} + (1-\theta) w_{\mathrm{old}}\right] \,\varphi_j\,\mathrm{d}x$  。
 * 

 * 
 * 
 * @code
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             local_nl_term(i) += 
 *               (std::sin(theta * new_data_values[q_point] + 
 *                         (1 - theta) * old_data_values[q_point]) * 
 *                fe_values.shape_value(i, q_point) * fe_values.JxW(q_point)); 
 * 
 * @endcode
 * 
 * 我们通过将各单元的积分对全局积分的贡献相加来得出结论。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           nl_term(local_dof_indices[i]) += local_nl_term(i); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemcompute_nl_matrix"></a> 
 * <h4>SineGordonProblem::compute_nl_matrix</h4>
 * 

 * 
 * 这是处理非线性方案的第二个函数。它计算矩阵  $N(\cdot,\cdot)$  ，它出现在  $F(\cdot)$  的雅各布项的非线性项中。正如 <code>compute_nl_term</code> 一样，我们必须让这个函数接收一个 "旧 "和一个 "新 "的解决方案作为输入，我们再次将其分别称为 $w_{\mathrm{old}}$ 和 $w_{\mathrm{new}}$ ，如下。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::compute_nl_matrix( 
 *     const Vector<double> &old_data, 
 *     const Vector<double> &new_data, 
 *     SparseMatrix<double> &nl_matrix) const 
 *   { 
 *     QGauss<dim>   quadrature_formula(fe.degree + 1); 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_JxW_values | 
 *                               update_quadrature_points); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> local_nl_matrix(dofs_per_cell, dofs_per_cell); 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 *     std::vector<double>                  old_data_values(n_q_points); 
 *     std::vector<double>                  new_data_values(n_q_points); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         local_nl_matrix = 0; 
 * 
 * @endcode
 * 
 * 同样，首先我们将我们的 <code>FEValues</code> 实例化重新初始化为当前单元。
 * 

 * 
 * 
 * @code
 *         fe_values.reinit(cell); 
 *         fe_values.get_function_values(old_data, old_data_values); 
 *         fe_values.get_function_values(new_data, new_data_values); 
 * 
 * @endcode
 * 
 * 然后，我们用所需的正交公式评估 $\int_K \cos\left[\theta w_{\mathrm{new}} + (1-\theta) w_{\mathrm{old}}\right]\, \varphi_i\, \varphi_j\,\mathrm{d}x$ 。
 * 

 * 
 * 
 * @code
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               local_nl_matrix(i, j) += 
 *                 (std::cos(theta * new_data_values[q_point] + 
 *                           (1 - theta) * old_data_values[q_point]) * 
 *                  fe_values.shape_value(i, q_point) * 
 *                  fe_values.shape_value(j, q_point) * fe_values.JxW(q_point)); 
 * 
 * @endcode
 * 
 * 最后，我们将各单元的积分对全局积分的贡献相加。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             nl_matrix.add(local_dof_indices[i], 
 *                           local_dof_indices[j], 
 *                           local_nl_matrix(i, j)); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemsolve"></a> 
 * <h4>SineGordonProblem::solve</h4>
 * 

 * 
 * 正如在介绍中所讨论的，这个函数在线性方程组上使用CG迭代求解器，该方程组是由牛顿方法的每个迭代的有限元空间离散化产生的，用于分割公式中的（非线性）第一个方程。该系统的解实际上是 $\delta U^{n,l}$ ，所以它被存储在 <code>solution_update</code> and used to update <code>solution</code> 的 <code>run</code> 函数中。
 * 

 * 
 * 注意，我们在求解前将解的更新值重新设置为零。这是没有必要的：迭代求解器可以从任何一点开始并收敛到正确的解。如果对线性系统的解有一个很好的估计，那么从这个向量开始可能是值得的，但是作为一个一般的观察，起点并不是很重要：它必须是一个非常非常好的猜测，以减少超过几个迭代的次数。事实证明，对于这个问题，使用之前的非线性更新作为起点实际上会损害收敛性并增加所需的迭代次数，所以我们简单地将其设置为零。
 * 

 * 
 * 该函数返回收敛到一个解决方案所需的迭代次数。这个数字以后将被用来在屏幕上生成输出，显示每次非线性迭代需要多少次迭代。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   unsigned int SineGordonProblem<dim>::solve() 
 *   { 
 *     SolverControl            solver_control(1000, 1e-12 * system_rhs.l2_norm()); 
 *     SolverCG<Vector<double>> cg(solver_control); 
 * 
 *     PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *     preconditioner.initialize(system_matrix, 1.2); 
 * 
 *     cg.solve(system_matrix, solution_update, system_rhs, preconditioner); 
 * 
 *     return solver_control.last_step(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemoutput_results"></a> 
 * <h4>SineGordonProblem::output_results</h4>
 * 

 * 
 * 这个函数将结果输出到一个文件。它与  step-23  和  step-24  中的相应函数基本相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::output_results( 
 *     const unsigned int timestep_number) const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(solution, "u"); 
 *     data_out.build_patches(); 
 * 
 *     const std::string filename = 
 *       "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu"; 
 *     DataOutBase::VtkFlags vtk_flags; 
 *     vtk_flags.compression_level = 
 *       DataOutBase::VtkFlags::ZlibCompressionLevel::best_speed; 
 *     data_out.set_flags(vtk_flags); 
 *     std::ofstream output(filename); 
 *     data_out.write_vtu(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="SineGordonProblemrun"></a> 
 * <h4>SineGordonProblem::run</h4>
 * 

 * 
 * 这个函数对一切都有最高级别的控制：它运行（外部）时间步长循环，（内部）非线性求解器循环，并在每个时间步长后输出解。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void SineGordonProblem<dim>::run() 
 *   { 
 *     make_grid_and_dofs(); 
 * 
 * @endcode
 * 
 * 为了确认初始条件，我们必须使用函数  $u_0(x)$  来计算  $U^0$  。为此，下面我们将创建一个 <code>InitialValues</code> 类型的对象；注意，当我们创建这个对象（它来自 <code>Function</code> 类）时，我们将其内部的时间变量设置为 $t_0$ ，以表明初始条件是在 $t=t_0$ 处评估的空间和时间的函数。
 * 

 * 
 * 然后我们通过使用 <code>VectorTools::project</code> 将 $u_0(x)$ 投影到网格上，产生 $U^0$ 。我们必须使用与 step-21 相同的悬挂节点约束结构： VectorTools::project 函数需要一个悬挂节点约束对象，但为了使用它，我们首先需要关闭它。
 * 

 * 
 * 
 * @code
 *     { 
 *       AffineConstraints<double> constraints; 
 *       constraints.close(); 
 *       VectorTools::project(dof_handler, 
 *                            constraints, 
 *                            QGauss<dim>(fe.degree + 1), 
 *                            InitialValues<dim>(1, time), 
 *                            solution); 
 *     } 
 * 
 * @endcode
 * 
 * 为了完整起见，我们像其他时间步长一样，将第2个时间步长输出到一个文件。
 * 

 * 
 * 
 * @code
 *     output_results(0); 
 * 
 * @endcode
 * 
 * 现在我们进行时间步进：在每个时间步进中，我们解决与问题的有限元离散化相对应的矩阵方程，然后根据我们在介绍中讨论的时间步进公式推进我们的解决方案。
 * 

 * 
 * 
 * @code
 *     unsigned int timestep_number = 1; 
 *     for (time += time_step; time <= final_time; 
 *          time += time_step, ++timestep_number) 
 *       { 
 *         old_solution = solution; 
 * 
 *         std::cout << std::endl 
 *                   << "Time step #" << timestep_number << "; " 
 *                   << "advancing to t = " << time << "." << std::endl; 
 * 
 * @endcode
 * 
 * 在每个时间步长的开始，我们必须通过牛顿方法求解拆分公式中的非线性方程---即先求解 $\delta U^{n,l}$ ，然后再计算 $U^{n,l+1}$ ，如此反复。这种非线性迭代的停止标准是： $\|F_h(U^{n,l})\|_2 \le 10^{-6} \|F_h(U^{n,0})\|_2$  。因此，我们需要记录第一次迭代中残差的规范。
 * 

 * 
 * 在每次迭代结束时，我们向控制台输出我们花了多少次线性求解器的迭代。当下面的循环完成后，我们有（一个近似的） $U^n$  。
 * 

 * 
 * 
 * @code
 *         double initial_rhs_norm = 0.; 
 *         bool   first_iteration  = true; 
 *         do 
 *           { 
 *             assemble_system(); 
 * 
 *             if (first_iteration == true) 
 *               initial_rhs_norm = system_rhs.l2_norm(); 
 * 
 *             const unsigned int n_iterations = solve(); 
 * 
 *             solution += solution_update; 
 * 
 *             if (first_iteration == true) 
 *               std::cout << "    " << n_iterations; 
 *             else 
 *               std::cout << '+' << n_iterations; 
 *             first_iteration = false; 
 *           } 
 *         while (system_rhs.l2_norm() > 1e-6 * initial_rhs_norm); 
 * 
 *         std::cout << " CG iterations per nonlinear step." << std::endl; 
 * 
 * @endcode
 * 
 * 在得到问题的第一个方程 $t=t_n$ 的解后，我们必须更新辅助速度变量  $V^n$  。然而，我们不计算和存储 $V^n$ ，因为它不是我们在问题中直接使用的数量。因此，为了简单起见，我们直接更新 $MV^n$ 。
 * 

 * 
 * 
 * @code
 *         Vector<double> tmp_vector(solution.size()); 
 *         laplace_matrix.vmult(tmp_vector, solution); 
 *         M_x_velocity.add(-time_step * theta, tmp_vector); 
 * 
 *         laplace_matrix.vmult(tmp_vector, old_solution); 
 *         M_x_velocity.add(-time_step * (1 - theta), tmp_vector); 
 * 
 *         compute_nl_term(old_solution, solution, tmp_vector); 
 *         M_x_velocity.add(-time_step, tmp_vector); 
 * 
 * @endcode
 * 
 * 很多时候，特别是对于细网格，我们必须选择相当小的时间步长，以使方案稳定。因此，有很多时间步长，在解的过程中 "没有什么有趣的事情发生"。为了提高整体效率--特别是加快程序速度和节省磁盘空间--我们每隔 <code>output_timestep_skip</code> 个时间步数才输出解。
 * 

 * 
 * 
 * @code
 *         if (timestep_number % output_timestep_skip == 0) 
 *           output_results(timestep_number); 
 *       } 
 *   } 
 * } // namespace Step25 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 这是该程序的主函数。它创建一个顶层类的对象并调用其主函数。如果在执行 <code>SineGordonProblem</code> 类的运行方法时抛出了异常，我们会在这里捕获并报告它们。关于异常的更多信息，读者应该参考  step-6  。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step25; 
 * 
 *       SineGordonProblem<1> sg_problem; 
 *       sg_problem.run(); 
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
examples/step-25/doc/results.dox



<a name="Results"></a><h1>Results</h1>显式欧拉时间步长方案（ $\theta=0$ ）对于我们希望解决的问题来说表现得很充分。不幸的是，由于稳定性问题，必须选择一个相当小的时间步长 ----  $k\sim h/10$ 似乎对我们进行的大多数模拟都有效。另一方面，Crank-Nicolson方案（ $\theta=\frac{1}{2}$ ）是无条件稳定的，而且（至少对于一维呼吸器的情况）我们可以选择大到 $25h$ 的时间步长而不对解决方案产生任何不良影响。隐式欧拉方案（ $\theta=1$ ）是 "指数阻尼 "的，所以它不是解决正弦戈登方程的好选择，因为它是保守的。然而， $\theta$ 方法提供的连续体中的一些阻尼方案对于消除边界效应引起的虚假振荡很有用。


在下面的模拟中，我们在一维的区间 $\Omega =
[-10,10]$ 和二维的正方形 $\Omega = [-10,10]\times [-10,10]$ 上解决正弦-戈登方程。在每种情况下，各自的网格被均匀地细化了6次，即 $h\sim
2^{-6}$  。

<a name="An11dSolution"></a><h3>An (1+1)-d Solution</h3> 我们讨论的第一个例子是正弦-戈登方程的所谓一维（静止的）呼吸器解。正如介绍中提到的，呼吸器有如下闭合形式的表达。


\f[
u_{\mathrm{breather}}(x,t) = -4\arctan \left(\frac{m}{\sqrt{1-m^2}} \frac{\sin\left(\sqrt{1-m^2}t +c_2\right)}{\cosh(mx+c_1)} \right),
\f] 其中 $c_1$  ,  $c_2$  和  $m<1$  是常数。在下面的模拟中，我们选择了  $c_1=0$  ,  $c_2=0$  ,  $m=0.5$  。此外，我们知道呼吸器的振荡周期是 $2\pi\sqrt{1-m^2}$ ，因此我们选择 $t_0=-5.4414$ 和 $t_f=2.7207$ ，这样我们可以观察到解决方案的三次振荡。然后，取 $u_0(x) = u_{\mathrm{breather}}(x,t_0)$ 、 $\theta=0$ 和 $k=h/10$ ，程序计算出以下解。

 <img src="https://www.dealii.org/images/steps/developer/step-25.1d-breather.gif" alt="Animation of the 1D stationary breather."> 

虽然程序中没有显示如何做，但另一种可视化(1+1)-d解决方案的方法是使用DataOutStack类产生的输出；它允许 "堆叠 "单个时间步骤的解决方案，因此我们从一维时间相关的解决方案中得到二维时空图。这就产生了下面的时空图，而不是上面的动画。

 <img src="https://www.dealii.org/images/steps/developer/step-25.1d-breather_stp.png" alt="A space-time plot of the 1D stationary breather."> 

此外，由于呼吸器是正弦-戈登方程的分析解，我们可以用它来验证我们的代码，尽管我们必须假设我们选择的诺伊曼边界条件引入的误差与数值误差相比是很小的。在这种假设下，可以使用 VectorTools::integrate_difference 函数来计算数值解和本程序的 <code>ExactSolution</code> 类描述的函数之间的差异。对于上面两幅图所示的模拟，每个时间步长的有限元解的 $L^2$ 误差的规范保持在 $10^{-2}$ 的数量级。因此，我们可以得出结论，数值方法在程序中得到了正确的实现。




<a name="Afew21DSolutions"></a><h3>A few (2+1)D Solutions</h3>


在文献中可以找到的正弦戈登方程在(2+1)D中的唯一分析解是所谓的Kink孤波。它有以下的闭合式表达。   @f[
    u(x,y,t) = 4 \arctan \left[a_0 e^{s\xi}\right]
  @f]

与@f[
    \xi = x \cos\vartheta + \sin(\vartheta) (y\cosh\lambda + t\sinh \lambda)
  @f]

其中 $a_0$ 、 $\vartheta$ 和 $\lambda$ 为常数。在下面的模拟中，我们选择了 $a_0=\lambda=1$  。请注意，如果 $\vartheta=\pi$ 是静止的，那么它将是一个很好的解决方案，我们可以用它来验证二维的程序，因为没有发生域边界的反射。

下面显示的模拟是用 $u_0(x) = u_{\mathrm{kink}}(x,t_0)$ 、 $\theta=\frac{1}{2}$ 、 $k=20h$ 、 $t_0=1$ 和 $t_f=500$ 进行的。每个时间步长的有限元解的误差的 $L^2$ 准则保持在 $10^{-2}$ 的数量级上，表明该程序在二维和一维中都能正常工作。不幸的是，这个解决方案并不十分有趣，不过为了完整起见，我们还是在下面附上了它的快照。

 <img src="https://www.dealii.org/images/steps/developer/step-25.2d-kink.png" alt="Stationary 2D kink."> 

现在我们已经在一维和二维中验证了代码，我们转到一个分析解未知的问题。

为此，我们围绕 $z$ 轴旋转上面讨论的扭结解：我们让 $\vartheta=\frac{\pi}{4}$  .后者的结果是一个不与网格对齐的孤波，所以反射立即发生在域的边界。对于下面所示的模拟，我们采取了 $u_0(x)=u_{\mathrm{kink}}(x,t_0)$ 、 $\theta=\frac{2}{3}$ 、 $k=20h$ 、 $t_0=0$ 和 $t_f=20$  。此外，我们不得不选择 $\theta=\frac{2}{3}$ ，因为对于任何 $\theta\le\frac{1}{2}$ 的边界都会产生振荡，这可能是由于方案而不是方程造成的，因此选择 $\theta$ 的值，在时间步进方案的 "指数阻尼 "频谱中，确保这些振荡不会被产生。

 <img src="https://www.dealii.org/images/steps/developer/step-25.2d-angled_kink.gif" alt="Animation of a moving 2D kink, at 45 degrees to the axes of the grid, showing boundary effects."> 

正弦-戈登方程的另一个有趣的解决方案（不能通过分析获得）可以通过使用两个一维呼吸器来构建以下可分离的二维初始条件来产生。

\f[
  u_0(x) =
  u_{\mathrm{pseudobreather}}(x,t_0) =
  16\arctan \left(
    \frac{m}{\sqrt{1-m^2}}
    \frac{\sin\left(\sqrt{1-m^2}t_0\right)}{\cosh(mx_1)} \right)
  \arctan \left(
    \frac{m}{\sqrt{1-m^2}}
    \frac{\sin\left(\sqrt{1-m^2}t_0\right)}{\cosh(mx_2)} \right),
\f] 其中 $x=(x_1,x_2)\in{R}^2$  ,  $m=0.5<1$  与我们上面讨论的一维情况一样。在下面的模拟中，我们选择了 $\theta=\frac{1}{2}$ 、 $k=10h$ 、 $t_0=-5.4414$ 和 $t_f=2.7207$ 。解决方案是相当有趣的



 <img src="https://www.dealii.org/images/steps/developer/step-25.2d-pseudobreather.gif" alt="Animation of a 2D pseudobreather."> 


<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


改变初始条件是有意义的。大多数选择不会导致保持局部的解决方案（在孤子界，这样的解决方案被称为 "静止的"，尽管解决方案确实随时间变化），而是导致方程的波状特征占主导地位的解决方案，并且波从局部初始条件的位置传播出去。例如，值得玩一玩 <code>InitialValues</code> 类，把对 <code>ExactSolution</code> 类的调用换成类似这样的函数。

@f[
  u_0(x,y) = \cos\left(\frac x2\right)\cos\left(\frac y2\right)


@f]

如果 $|x|,|y|\le \frac\pi 2$ ，和 $u_0(x,y)=0$ 在这个区域之外。

第二个方面是研究该方案是否是能量保持的。对于 @ref
step_23 "step-23 "中讨论的纯波方程，如果我们选择时间步进参数，使我们得到Crank-Nicolson方案，情况就是这样。我们可以在这里做类似的事情，注意到正弦-戈登解中的能量被定义为

@f[
  E(t) = \frac 12 \int_\Omega \left(\frac{\partial u}{\partial
  t}\right)^2
  + \left(\nabla u\right)^2 + 2 (1-\cos u) \; dx.


@f]

(我们在公式中使用 $1-\cos u$ 而不是 $-\cos u$ ，以确保对能量的所有贡献都是正的，从而使衰变的解在无界域上具有有限的能量。)

除此以外，还有两个明显的领域。

- 显然，自适应性（即时间自适应网格）会对这样的问题产生兴趣。它们的复杂性导致我们再次将其排除在本方案之外，尽管 @ref step_23 "step-23 "介绍中的一般评论仍然正确。

- 解决这个问题的更快方案。虽然今天的计算机已经足够快，可以在不长的时间内解决二维甚至三维静止的问题，但与时间相关的问题则是完全不同的一类问题。我们在步骤48中讨论了这个问题，我们展示了如何在不组装或倒置任何矩阵的情况下并行解决这个问题。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-25.cc"
*/
