/**
@page step_36 The step-36 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#EigenvaluesandDirichletboundaryconditions">Eigenvalues and Dirichlet boundary conditions</a>
        <li><a href="#Implementationdetails">Implementation details</a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeEigenvalueProblemcodeclasstemplate">The <code>EigenvalueProblem</code> class template</a>
        <li><a href="#ImplementationofthecodeEigenvalueProblemcodeclass">Implementation of the <code>EigenvalueProblem</code> class</a>
      <ul>
        <li><a href="#EigenvalueProblemEigenvalueProblem">EigenvalueProblem::EigenvalueProblem</a>
        <li><a href="#EigenvalueProblemmake_grid_and_dofs">EigenvalueProblem::make_grid_and_dofs</a>
        <li><a href="#EigenvalueProblemassemble_system">EigenvalueProblem::assemble_system</a>
        <li><a href="#EigenvalueProblemsolve">EigenvalueProblem::solve</a>
        <li><a href="#EigenvalueProblemoutput_results">EigenvalueProblem::output_results</a>
        <li><a href="#EigenvalueProblemrun">EigenvalueProblem::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Runningtheproblem">Running the problem</a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-36/doc/intro.dox

 <br> 

<i>This program was contributed by Toby D. Young and Wolfgang
Bangerth.  </i>

<a name="Preamble"></a>

<a name="Preamble"></a><h1>Preamble</h1>


在这个例子中，我们要解决的问题是一个特征谱问题。特征值问题出现在广泛的问题背景中，例如在计算腔体中的电磁驻波、鼓膜的振动模式或湖泊和河口的振荡中。最神秘的应用之一可能是量子力学中静止或准静止波函数的计算。后一种应用是我们在此想研究的，尽管本方案中概述的一般技术当然同样适用于上述其他应用。

Eigenspectrum问题的一般形式是

@f{align*}
	L \Psi &= \varepsilon \Psi \qquad &&\text{in}\ \Omega\quad,
	\\
	\Psi &= 0 &&\text{on}\ \partial\Omega\quad,


@f}

其中 $\Psi=\Psi(\mathbf x)$ 上的迪里切特边界条件也可以用诺伊曼或罗宾条件代替； $L$ 是一个算子，一般也包含微分算子。

在适当的条件下，上述方程有一组解 $\Psi_\ell,\varepsilon_\ell$  ,  $\ell\in {\cal I}$  ，其中 $\cal I$ 可以是一个有限的或无限的集合（在后一种情况下，它可能是一个离散的或有时至少是部分连续的集合）。在这两种情况下，让我们注意到，我们要计算的不再只是一个单一的解决方案，而是一组解决方案（各种特征函数和相应的特征值）。从数字上寻找这类特征值问题的所有特征值（特征函数）的问题是一个艰巨的挑战。事实上，如果集合 $\cal I$ 是无限的，这个挑战当然是难以解决的。  然而，大多数时候，我们实际上只对这些值（函数）的一小部分感兴趣；幸运的是，我们将在这个教程程序中使用的SLEPc库的接口允许我们选择特征谱的哪一部分和我们想要解决的多少个解决方案。

在这个程序中，我们使用的eigenspectrum求解器是由deal.II提供的类，围绕<a href="http://www.grycap.upv.es/slepc/" target="_top">SLEPc</a>库的线性代数实现；SLEPc本身建立在<a
href="http://www.mcs.anl.gov/petsc/" target="_top">PETSc</a>库的线性代数内容之上。

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


静止的量子力学的基本方程是薛定谔方程，它模拟粒子在外部势中的运动  $V(\mathbf x)$  。粒子由一个波函数 $\Psi(\mathbf x)$ 描述，它满足一个（非维度化）形式的关系

@f{align*} [-\Delta + V(\mathbf x)]
\Psi(\mathbf x) &= \varepsilon \Psi(\mathbf x) \qquad &&\text{in}\
\Omega\quad, \\ \Psi &= 0 &&\text{on}\ \partial\Omega\quad.


@f}

因此，这个粒子只能存在于一定数量的特征态中，这些特征态对应于作为该方程解的能量特征值 $\varepsilon_\ell$ 。量子力学的正统（哥本哈根）解释认为，如果一个粒子具有能量 $\varepsilon_\ell$ ，那么在位置 $\mathbf x$ 找到它的概率与 $|\Psi_\ell(\mathbf
x)|^2$ 成正比，其中 $\Psi_\ell$ 是对应于这个特征值的特征函数。

为了在数值上找到这个方程的解，即一组成对的特征值/特征函数，我们使用通常的有限元方法，将左边的方程与测试函数相乘，通过部分积分，在有限维空间中通过近似 $\Psi(\mathbf
x)\approx\Psi_h(\mathbf x)=\sum_{j}\phi_j(\mathbf x)\tilde\psi_j$ 寻找解，其中 $\tilde\psi$ 是一个扩展系数的矢量。然后，我们立即得出以下方程，将连续特征值问题离散化：@f[ \sum_j [(\nabla\phi_i,
\nabla\phi_j)+(V(\mathbf x)\phi_i,\phi_j)] \tilde{\psi}_j =
\varepsilon_h \sum_j (\phi_i, \phi_j) \tilde{\psi}_j\quad.  @f] 在矩阵和矢量符号中，这个方程然后读作。@f[ A
\tilde{\Psi} = \varepsilon_h M \tilde{\Psi} \quad, @f] 其中 $A$ 是由微分算子 $L$ 产生的刚度矩阵，而 $M$ 是质量矩阵。特征值问题的解决方案是一个特征谱  $\varepsilon_{h,\ell}$  ，以及相关的特征函数  $\Psi_\ell=\sum_j \phi_j\tilde{\psi}_j$  。




<a name="EigenvaluesandDirichletboundaryconditions"></a><h3>Eigenvalues and Dirichlet boundary conditions</h3>


在这个程序中，我们对波函数 $\Psi$ 使用了Dirichlet边界条件。从有限元代码的角度来看，这意味着只有内部的自由度是<i>freedom</i>的实数度：边界上的自由度不是自由的，但毕竟是被迫有一个零值。另一方面，有限元方法的力量和简单性主要来自于我们只是在每个单元上做同样的事情，而不需要过多考虑一个单元在哪里，它是否在一个不太精细的单元上，因此有一个悬挂的节点，或者与边界相邻。所有这些检查都会使有限元线性系统的组装变得难以忍受，甚至更难阅读。

因此，当然，当你用DoFHandler对象分配自由度时，你并不关心你列举的一些自由度是否处于迪里切特边界。他们都会得到数字。我们只需要在以后应用边界值的时候照顾到这些自由度。有两种基本的方法（要么使用 MatrixTools::apply_boundary_values() <i>after</i>装配线性系统，要么使用 AffineConstraints::distribute_local_to_global() <i>during</i>装配；更多信息见 @ref constraints "约束模块"），但两者的结果都一样：一个线性系统的总行数等于<i>all</i>自由度的数量，包括那些位于边界的自由度。然而，受迪里希特条件约束的自由度与线性系统的其他部分分开，方法是将相应的行和列清零，在对角线上放一个正条目，并在右手边放相应的迪里希特值。

如果你暂时假设我们将自由度重新编号，使迪里切特边界上的所有自由度排在最后，那么我们在求解一个有右手边的常规PDE时得到的线性系统将是这样的。

@f{align*}
  \begin{pmatrix}
    A_i & 0 \\ 0 & D_b
  \end{pmatrix}
  \begin{pmatrix}
    U_i \\ U_b
  \end{pmatrix}
  =
  \begin{pmatrix}
    F_i \\ F_b
  \end{pmatrix}.


@f}

这里，下标 $i$ 和 $b$ 分别对应于内部和边界自由度。内部自由度满足线性系统 $A_i U_i=F_i$ ，在内部产生正确的解，边界值由 $U_b = D_b^{-1} F_b$ 决定，其中 $D_b$ 是一个对角矩阵，由消除边界自由度的过程产生， $F_b$ 是以这样的方式选择的： $U_{b,j}=D_{b,jj}^{-1} F_{b,j}$ 对每个边界自由度 $j$ 都有正确的边界值。 对于好奇的人来说，矩阵 $D_b$ 的条目是将修改后的局部贡献加入全局矩阵的结果，对于局部矩阵，如果非零，对角线元素被设置为其绝对值；否则，它们被设置为对角线的绝对值的平均值。这个过程保证了 $D_b$ 的条目是正的，并且大小与对角线的其他条目相当，确保所产生的矩阵不会因为涉及到大小迥异的矩阵条目的舍入而产生不合理的精度损失。最终出现在对角线上的实际数值是很难预测的，你应该把它们当作任意的、不可预测的，但却是正的。)

对于 "常规 "的线性系统，这一切都导致了正确的解决方案。另一方面，对于特征值问题，这就不那么简单了。在那里，消除边界值会影响到我们在当前教程程序中要解决的矩阵 $A$ 和 $M$ 。消除边界值后，我们就会得到一个可以这样划分的特征值问题。

@f{align*}
  \begin{pmatrix}
    A_i & 0 \\ 0 & D_A
  \end{pmatrix}
  \begin{pmatrix}
    \tilde\Psi_i \\ \tilde\Psi_b
  \end{pmatrix}
  =
  \epsilon_h
  \begin{pmatrix}
    M_i & 0 \\ 0 & D_M
  \end{pmatrix}
  \begin{pmatrix}
    \tilde\Psi_i \\ \tilde\Psi_b
  \end{pmatrix}.


@f}

这种形式清楚地表明，有两组特征值：我们关心的那些，以及来自分离问题的虚假特征值。

@f[
  D_A \tilde \Psi_b = \epsilon_h D_M \Psi_b.


@f]

这些特征值是虚假的，因为它们是由一个只在边界节点上操作的特征值系统产生的--这些节点不是<i>freedom</i>的实数度。当然，由于两个矩阵 $D_A,D_M$ 是对角线，我们可以准确地量化这些虚假的特征值：它们是 $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$ （其中指数 $j$ 正好对应于受迪里切特边界值约束的自由度）。

那么，如何处理这些问题呢？第一部分是识别我们的特征值求解器何时找到其中一个。为此，程序通过计算所有受限自由度上表达式 $\varepsilon_{h,j}=D_{A,jj}/D_{M,jj}$ 的最小和最大，计算并打印出这些特征值所在的区间。在下面的程序中，这已经足够了：我们发现这个区间位于我们感兴趣并计算的最小特征值和相应的特征函数的集合之外，所以这里我们不需要做什么。

另一方面，我们可能会发现我们在这个程序中计算的一个特征值恰好在这个区间内，在这种情况下，我们不会立即知道它是一个虚假的还是一个真正的特征值。在这种情况下，我们可以在计算完两个矩阵后简单地缩放其中一个矩阵的对角线元素，从而将它们从特征谱的感兴趣的频率上移开。这可以通过使用以下代码来完成，确保所有假的特征值都正好等于 $1.234\cdot 10^5$  。

@code
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
      if (constraints.is_constrained(i))
        {
          stiffness_matrix.set(i, i, 1.234e5);
          mass_matrix.set(i, i, 1);
        }
@endcode

然而，这里没有采用这种策略，因为我们从程序中得到的虚假特征值恰好大于我们将计算的和感兴趣的最低的五个。




<a name="Implementationdetails"></a><h3>Implementation details</h3>


下面的程序实质上只是步骤4的一个稍加修改的版本。有所不同的是以下几点。

 <ul> 

 <li>  主类（名为  <code>EigenvalueProblem</code>  ）现在不再有单一的解向量，而是有一整套我们想要计算的各种特征函数的向量。此外， <code>main</code> 函数对这里的一切都有顶层控制，它通过 <code>SlepcInitialize</code> 和 <code>SlepFinalize</code> 同时初始化和最终确定SLEPc和PETSc的接口。 </li> 

 <li>  我们在步骤17和步骤18中使用PETSc矩阵和向量，因为这是SLEPc特征值求解器所要求的。 </li> 

 <li>  函数 <code>EigenvalueProblem::solve</code> 与教程中迄今为止所见的任何函数都完全不同，因为它不只是求解一个线性系统，而是实际求解特征值问题。它建立在SLEPc库上，更直接的是建立在类 SLEPcWrappers::SolverKrylovSchur.</li> 的deal.II SLEPc包装器上。

 <li>  我们使用ParameterHandler类来描述一些输入参数，如势的确切形式 $V({\mathbf
x})$  ，网格的全局细化步数，或我们要解决的特征值的数量。我们可以在这方面做得更多，但只限于在运行时选择一些实际的输入文件参数。为了看看在这方面可以做什么，看看 @ref step_29 "步骤-29 "和步骤-33。 </li> 

 <li>  我们使用FunctionParser类使潜在  $V(\mathbf
x)$  的运行时参数，可以在输入文件中指定为公式。 </li> 

 </ul> 

程序的其余部分以一种相当直接的方式从第4步开始。


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
 * 正如介绍中提到的，本程序基本上只是  step-4  的一个小修改版本。因此，以下大部分的include文件都是在那里使用的，或者至少是在以前的教程程序中已经使用的。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/function_parser.h> 
 * #include <deal.II/base/parameter_handler.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * 
 * @endcode
 * 
 * IndexSet用于设置每个  PETScWrappers::MPI::Vector:  的大小。
 * 
 * @code
 * #include <deal.II/base/index_set.h> 
 * 
 * @endcode
 * 
 * PETSc出现在这里是因为SLEPc依赖于这个库。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/petsc_sparse_matrix.h> 
 * #include <deal.II/lac/petsc_vector.h> 
 * 
 * @endcode
 * 
 * 然后我们需要实际导入SLEPc提供的求解器接口。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/slepc_solver.h> 
 * 
 * @endcode
 * 
 * 我们还需要一些标准的C++。
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后，和以前的程序一样，我们将所有的deal.II类和函数名导入到本程序中所有的名字空间中。
 * 

 * 
 * 
 * @code
 * namespace Step36 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeEigenvalueProblemcodeclasstemplate"></a> 
 * <h3>The <code>EigenvalueProblem</code> class template</h3>
 * 

 * 
 * 下面是主类模板的类声明。它看起来和在  step-4  中已经展示过的差不多了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class EigenvalueProblem 
 *   { 
 *   public: 
 *     EigenvalueProblem(const std::string &prm_file); 
 *     void run(); 
 * 
 *   private: 
 *     void         make_grid_and_dofs(); 
 *     void         assemble_system(); 
 *     unsigned int solve(); 
 *     void         output_results() const; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FE_Q<dim>          fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 * @endcode
 * 
 * 有了这些例外情况。对于我们的特征值问题，我们既需要左手边的刚度矩阵，也需要右手边的质量矩阵。我们还需要的不仅仅是一个解函数，而是一整套我们想要计算的特征函数，以及相应的特征值。
 * 

 * 
 * 
 * @code
 *     PETScWrappers::SparseMatrix             stiffness_matrix, mass_matrix; 
 *     std::vector<PETScWrappers::MPI::Vector> eigenfunctions; 
 *     std::vector<double>                     eigenvalues; 
 * 
 * @endcode
 * 
 * 然后，我们需要一个对象来存储几个运行时参数，我们将在输入文件中指定。
 * 

 * 
 * 
 * @code
 *     ParameterHandler parameters; 
 * 
 * @endcode
 * 
 * 最后，我们将有一个对象，包含对我们自由度的 "约束"。如果我们有自适应细化的网格（目前的程序中没有），这可能包括悬挂节点约束。这里，我们将存储边界节点的约束  $U_i=0$  。
 * 

 * 
 * 
 * @code
 *     AffineConstraints<double> constraints; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="ImplementationofthecodeEigenvalueProblemcodeclass"></a> 
 * <h3>Implementation of the <code>EigenvalueProblem</code> class</h3>
 * 
 * <a name="EigenvalueProblemEigenvalueProblem"></a> 
 * <h4>EigenvalueProblem::EigenvalueProblem</h4>
 * 

 * 
 * 首先是构造函数。主要的新部分是处理运行时的输入参数。我们需要首先声明它们的存在，然后从输入文件中读取它们的值，该文件的名称被指定为该函数的参数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   EigenvalueProblem<dim>::EigenvalueProblem(const std::string &prm_file) 
 *     : fe(1) 
 *     , dof_handler(triangulation) 
 *   { 
 * 
 * @endcode
 * 
 * TODO研究为什么获得正确的特征值退化所需的最小细化步骤数为6
 * 

 * 
 * 
 * @code
 *     parameters.declare_entry( 
 *       "Global mesh refinement steps", 
 *       "5", 
 *       Patterns::Integer(0, 20), 
 *       "The number of times the 1-cell coarse mesh should " 
 *       "be refined globally for our computations."); 
 *     parameters.declare_entry("Number of eigenvalues/eigenfunctions", 
 *                              "5", 
 *                              Patterns::Integer(0, 100), 
 *                              "The number of eigenvalues/eigenfunctions " 
 *                              "to be computed."); 
 *     parameters.declare_entry("Potential", 
 *                              "0", 
 *                              Patterns::Anything(), 
 *                              "A functional description of the potential."); 
 * 
 *     parameters.parse_input(prm_file); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemmake_grid_and_dofs"></a> 
 * <h4>EigenvalueProblem::make_grid_and_dofs</h4>
 * 

 * 
 * 下一个函数在域 $[-1,1]^d$ 上创建一个网格，根据输入文件的要求对其进行多次细化，然后给它附加一个DoFHandler，将矩阵和向量初始化为正确的大小。我们还建立了对应于边界值的约束  $u|_{\partial\Omega}=0$  。
 * 

 * 
 * 对于矩阵，我们使用PETSc包装器。这些包装器能够在非零条目被添加时分配必要的内存。这看起来效率很低：我们可以先计算稀疏模式，用它来初始化矩阵，然后在我们插入条目时，我们可以确定我们不需要重新分配内存和释放之前使用的内存。一种方法是使用这样的代码。用
 * <div class=CodeFragmentInTutorialComment>
 * @code
 *    DynamicSparsityPattern
 *       dsp (dof_handler.n_dofs(),
 *            dof_handler.n_dofs());
 *    DoFTools::make_sparsity_pattern (dof_handler, dsp);
 *    dsp.compress ();
 *    stiffness_matrix.reinit (dsp);
 *    mass_matrix.reinit (dsp);
 *  @endcode
 * </div>
 * 代替下面两个 <code>reinit()</code> 的刚度和质量矩阵的调用。
 * 

 * 
 * 不幸的是，这并不完全可行。上面的代码可能会导致在非零模式下的一些条目，我们只写零条目；最值得注意的是，对于那些属于边界节点的行和列的非对角线条目，这一点是成立的。这不应该是一个问题，但是不管什么原因，PETSc的ILU预处理程序（我们用来解决特征值求解器中的线性系统）不喜欢这些额外的条目，并以错误信息中止。
 * 

 * 
 * 在没有任何明显的方法来避免这种情况的情况下，我们干脆选择第二种最好的方法，即让PETSc在必要时分配内存。也就是说，由于这不是一个时间上的关键部分，这整个事件就不再重要了。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EigenvalueProblem<dim>::make_grid_and_dofs() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation, -1, 1); 
 *     triangulation.refine_global( 
 *       parameters.get_integer("Global mesh refinement steps")); 
 *     dof_handler.distribute_dofs(fe); 
 * 
 *     DoFTools::make_zero_boundary_constraints(dof_handler, constraints); 
 *     constraints.close(); 
 * 
 *     stiffness_matrix.reinit(dof_handler.n_dofs(), 
 *                             dof_handler.n_dofs(), 
 *                             dof_handler.max_couplings_between_dofs()); 
 *     mass_matrix.reinit(dof_handler.n_dofs(), 
 *                        dof_handler.n_dofs(), 
 *                        dof_handler.max_couplings_between_dofs()); 
 * 
 * @endcode
 * 
 * 下一步是处理特征谱的问题。在这种情况下，输出是特征值和特征函数，所以我们将特征函数和特征值列表的大小设置为与我们在输入文件中要求的一样大。当使用 PETScWrappers::MPI::Vector, 时，Vector是使用IndexSet初始化的。IndexSet不仅用于调整 PETScWrappers::MPI::Vector 的大小，而且还将 PETScWrappers::MPI::Vector 中的一个索引与一个自由度联系起来（更详细的解释见 step-40 ）。函数complete_index_set()创建了一个IndexSet，每个有效的索引都是这个集合的一部分。请注意，这个程序只能按顺序运行，如果并行使用，将抛出一个异常。
 * 

 * 
 * 
 * @code
 *     IndexSet eigenfunction_index_set = dof_handler.locally_owned_dofs(); 
 *     eigenfunctions.resize( 
 *       parameters.get_integer("Number of eigenvalues/eigenfunctions")); 
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
 *       eigenfunctions[i].reinit(eigenfunction_index_set, MPI_COMM_WORLD); 
 * 
 *     eigenvalues.resize(eigenfunctions.size()); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemassemble_system"></a> 
 * <h4>EigenvalueProblem::assemble_system</h4>
 * 

 * 
 * 在这里，我们从局部贡献 $A^K_{ij} = \int_K \nabla\varphi_i(\mathbf x) \cdot \nabla\varphi_j(\mathbf x) + V(\mathbf x)\varphi_i(\mathbf x)\varphi_j(\mathbf x)$ 和 $M^K_{ij} = \int_K \varphi_i(\mathbf x)\varphi_j(\mathbf x)$ 中分别组合出全局刚度和质量矩阵。如果你看过以前的教程程序，这个函数应该会很熟悉。唯一新的东西是使用我们从输入文件中得到的表达式，设置一个描述势 $V(\mathbf x)$ 的对象。然后我们需要在每个单元的正交点上评估这个对象。如果你见过如何评估函数对象（例如，见 step-5 中的系数），这里的代码也会显得相当熟悉。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EigenvalueProblem<dim>::assemble_system() 
 *   { 
 *     QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_gradients | 
 *                               update_quadrature_points | update_JxW_values); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell); 
 *     FullMatrix<double> cell_mass_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     FunctionParser<dim> potential; 
 *     potential.initialize(FunctionParser<dim>::default_variable_names(), 
 *                          parameters.get("Potential"), 
 *                          typename FunctionParser<dim>::ConstMap()); 
 * 
 *     std::vector<double> potential_values(n_q_points); 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         cell_stiffness_matrix = 0; 
 *         cell_mass_matrix      = 0; 
 * 
 *         potential.value_list(fe_values.get_quadrature_points(), 
 *                              potential_values); 
 * 
 *         for (unsigned int q_point = 0; q_point < n_q_points; ++q_point) 
 *           for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *             for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *               { 
 *                 cell_stiffness_matrix(i, j) +=           // 
 *                   (fe_values.shape_grad(i, q_point) *    // 
 *                      fe_values.shape_grad(j, q_point)    // 
 *                    +                                     // 
 *                    potential_values[q_point] *           // 
 *                      fe_values.shape_value(i, q_point) * // 
 *                      fe_values.shape_value(j, q_point)   // 
 *                    ) *                                   // 
 *                   fe_values.JxW(q_point);                // 
 * 
 *                 cell_mass_matrix(i, j) +=              // 
 *                   (fe_values.shape_value(i, q_point) * // 
 *                    fe_values.shape_value(j, q_point)   // 
 *                    ) *                                 // 
 *                   fe_values.JxW(q_point);              // 
 *               } 
 * 
 * @endcode
 * 
 * 现在我们有了本地矩阵的贡献，我们把它们转移到全局对象中，并处理好零边界约束。
 * 

 * 
 * 
 * @code
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         constraints.distribute_local_to_global(cell_stiffness_matrix, 
 *                                                local_dof_indices, 
 *                                                stiffness_matrix); 
 *         constraints.distribute_local_to_global(cell_mass_matrix, 
 *                                                local_dof_indices, 
 *                                                mass_matrix); 
 *       } 
 * 
 * @endcode
 * 
 * 在函数的最后，我们告诉PETSc，矩阵现在已经完全组装好了，稀疏矩阵表示法现在可以被压缩了，因为不会再添加任何条目。
 * 

 * 
 * 
 * @code
 *     stiffness_matrix.compress(VectorOperation::add); 
 *     mass_matrix.compress(VectorOperation::add); 
 * 
 * @endcode
 * 
 * 在离开函数之前，我们计算虚假的特征值，这些特征值是由零Dirichlet约束引入到系统中的。正如介绍中所讨论的，使用Dirichlet边界条件，加上位于域的边界的自由度仍然是我们所求解的线性系统的一部分，引入了一些虚假的特征值。下面，我们输出它们所处的区间，以确保我们在计算中出现时可以忽略它们。
 * 

 * 
 * 
 * @code
 *     double min_spurious_eigenvalue = std::numeric_limits<double>::max(), 
 *            max_spurious_eigenvalue = -std::numeric_limits<double>::max(); 
 * 
 *     for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i) 
 *       if (constraints.is_constrained(i)) 
 *         { 
 *           const double ev         = stiffness_matrix(i, i) / mass_matrix(i, i); 
 *           min_spurious_eigenvalue = std::min(min_spurious_eigenvalue, ev); 
 *           max_spurious_eigenvalue = std::max(max_spurious_eigenvalue, ev); 
 *         } 
 * 
 *     std::cout << "   Spurious eigenvalues are all in the interval " 
 *               << "[" << min_spurious_eigenvalue << "," 
 *               << max_spurious_eigenvalue << "]" << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemsolve"></a> 
 * <h4>EigenvalueProblem::solve</h4>
 * 

 * 
 * 这是该程序的关键新功能。现在系统已经设置好了，现在是实际解决问题的好时机：和其他例子一样，这是使用 "解决 "程序来完成的。从本质上讲，它的工作原理与其他程序一样：你设置一个SolverControl对象，描述我们要解决的线性系统的精度，然后我们选择我们想要的解算器类型。这里我们选择了SLEPc的Krylov-Schur求解器，对于这类问题来说，这是一个相当快速和强大的选择。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   unsigned int EigenvalueProblem<dim>::solve() 
 *   { 
 * 
 * @endcode
 * 
 * 我们从这里开始，就像我们通常做的那样，指定我们想要的收敛控制。
 * 

 * 
 * 
 * @code
 *     SolverControl                    solver_control(dof_handler.n_dofs(), 1e-9); 
 *     SLEPcWrappers::SolverKrylovSchur eigensolver(solver_control); 
 * 
 * @endcode
 * 
 * 在我们实际求解特征函数和-值之前，我们还必须选择哪一组特征值来求解。让我们选择那些实部最小的特征值和相应的特征函数（事实上，我们在这里解决的问题是对称的，所以特征值是纯实部的）。之后，我们就可以真正让SLEPc做它的工作了。
 * 

 * 
 * 
 * @code
 *     eigensolver.set_which_eigenpairs(EPS_SMALLEST_REAL); 
 * 
 *     eigensolver.set_problem_type(EPS_GHEP); 
 * 
 *     eigensolver.solve(stiffness_matrix, 
 *                       mass_matrix, 
 *                       eigenvalues, 
 *                       eigenfunctions, 
 *                       eigenfunctions.size()); 
 * 
 * @endcode
 * 
 * 上述调用的输出是一组向量和数值。在特征值问题中，特征函数只确定到一个常数，这个常数可以很随意地固定。由于对特征值问题的原点一无所知，SLEPc除了将特征向量归一到 $l_2$ （向量）准则外，没有其他选择。不幸的是，这个规范与我们从特征函数角度可能感兴趣的任何规范没有什么关系： $L_2(\Omega)$ 规范，或者也许是 $L_\infty(\Omega)$ 规范。
 * 

 * 
 * 让我们选择后者，重新划分特征函数的尺度，使其具有 $\|\phi_i(\mathbf x)\|_{L^\infty(\Omega)}=1$ 而不是 $\|\Phi\|_{l_2}=1$ （其中 $\phi_i$ 是 $i$ 第三个特征<i>function</i>， $\Phi_i$ 是相应的结点值矢量）。对于这里选择的 $Q_1$ 元素，我们知道函数 $\phi_i(\mathbf x)$ 的最大值是在其中一个节点达到的，所以 $\max_{\mathbf x}\phi_i(\mathbf x)=\max_j (\Phi_i)_j$ ，使得在 $L_\infty$ 准则下的归一化是微不足道的。请注意，如果我们选择 $Q_k$ 元素与 $k>1$ ，这就不容易了：在那里，一个函数的最大值不一定要在一个节点上达到，所以 $\max_{\mathbf x}\phi_i(\mathbf x)\ge\max_j (\Phi_i)_j$ （尽管平等通常几乎是真的）。
 * 

 * 
 * 
 * @code
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
 *       eigenfunctions[i] /= eigenfunctions[i].linfty_norm(); 
 * 
 * @endcode
 * 
 * 最后返回收敛所需的迭代次数。
 * 

 * 
 * 
 * @code
 *     return solver_control.last_step(); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemoutput_results"></a> 
 * <h4>EigenvalueProblem::output_results</h4>
 * 

 * 
 * 这是本程序的最后一个重要功能。它使用DataOut类来生成特征函数的图形输出，以便以后进行可视化。它的工作原理与其他许多教程中的程序一样。
 * 

 * 
 * 整个函数的集合被输出为一个单一的VTK文件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EigenvalueProblem<dim>::output_results() const 
 *   { 
 *     DataOut<dim> data_out; 
 * 
 *     data_out.attach_dof_handler(dof_handler); 
 * 
 *     for (unsigned int i = 0; i < eigenfunctions.size(); ++i) 
 *       data_out.add_data_vector(eigenfunctions[i], 
 *                                std::string("eigenfunction_") + 
 *                                  Utilities::int_to_string(i)); 
 * 
 * @endcode
 * 
 * 唯一值得讨论的可能是，由于势在输入文件中被指定为函数表达式，因此最好能将其与特征函数一起以图形形式表示。实现这一目的的过程相对简单：我们建立一个代表 $V(\mathbf x)$ 的对象，然后将这个连续函数插值到有限元空间。我们还将结果附加到DataOut对象上，以便进行可视化。
 * 

 * 
 * 
 * @code
 *     Vector<double> projected_potential(dof_handler.n_dofs()); 
 *     { 
 *       FunctionParser<dim> potential; 
 *       potential.initialize(FunctionParser<dim>::default_variable_names(), 
 *                            parameters.get("Potential"), 
 *                            typename FunctionParser<dim>::ConstMap()); 
 *       VectorTools::interpolate(dof_handler, potential, projected_potential); 
 *     } 
 *     data_out.add_data_vector(projected_potential, "interpolated_potential"); 
 * 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output("eigenvectors.vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="EigenvalueProblemrun"></a> 
 * <h4>EigenvalueProblem::run</h4>
 * 

 * 
 * 这是一个对一切都有顶层控制的函数。它几乎与  step-4  中的内容完全相同。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void EigenvalueProblem<dim>::run() 
 *   { 
 *     make_grid_and_dofs(); 
 * 
 *     std::cout << "   Number of active cells:       " 
 *               << triangulation.n_active_cells() << std::endl 
 *               << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << std::endl; 
 * 
 *     assemble_system(); 
 * 
 *     const unsigned int n_iterations = solve(); 
 *     std::cout << "   Solver converged in " << n_iterations << " iterations." 
 *               << std::endl; 
 * 
 *     output_results(); 
 * 
 *     std::cout << std::endl; 
 *     for (unsigned int i = 0; i < eigenvalues.size(); ++i) 
 *       std::cout << "      Eigenvalue " << i << " : " << eigenvalues[i] 
 *                 << std::endl; 
 *   } 
 * } // namespace Step36 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 
 * @code
 * int main(int argc, char **argv) 
 * { 
 *   try 
 *     { 
 *       using namespace dealii; 
 *       using namespace Step36; 
 * 
 *       Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
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
 *                     "This program can only be run in serial, use ./step-36")); 
 * 
 *       EigenvalueProblem<2> problem("step-36.prm"); 
 *       problem.run(); 
 *     } 
 * 
 * @endcode
 * 
 * 在这期间，我们一直在注意是否有任何异常应该被生成。如果是这样的话，我们就会惊慌失措...
 * 

 * 
 * 
 * @code
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
 * @endcode
 * 
 * 如果没有抛出异常，我们就告诉程序不要再胡闹了，乖乖地退出。
 * 

 * 
 * 
 * @code
 *   std::cout << std::endl << "   Job done." << std::endl; 
 * 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-36/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Runningtheproblem"></a><h3>Running the problem</h3>


该问题的输入由一个输入文件 <code>\step-36.prm</code> 设定参数，例如，该文件可以包含以下文本。

@code
set Global mesh refinement steps         = 5
set Number of eigenvalues/eigenfunctions = 5
set Potential                            = 0
@endcode



这里，域内电势为零，我们知道特征值由 $\lambda_{(mn)}=\frac{\pi^2}{4}(m^2+n^2)$ 给出，其中 $m,n\in{\mathbb N^+}$  。特征函数是正弦和余弦，在 $m$ 和 $n$ 方向的周期为 $x$ 和 $y$ 。这与我们的程序产生的输出相匹配。

@code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       1024
   Number of degrees of freedom: 1089
   Solver converged in 67 iterations.


      Eigenvalue 0 : 4.93877
      Eigenvalue 1 : 12.3707
      Eigenvalue 2 : 12.3707
      Eigenvalue 3 : 19.8027
      Eigenvalue 4 : 24.837


   Job done.  @endcode 这些特征值正是对应于 $(m,n)=(1,1)$ 、 $(1,2)$ 和 $(2,1)$ 、 $(2,2)$ 和 $(3,1)$ 等对。相应的特征方程的可视化看起来是这样的。

 <table width="80%">
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.0.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.1.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.2.png" alt=""></td>
</tr>
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.3.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.default.eigenfunction.4.png" alt=""></td>
<td></td>
</tr>
</table> 

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


在操场上玩几个游戏总是值得的!所以这里有几个建议。

 <ul> 

 <li> 上面使用的势（称为<i>infinite well</i>，因为它是一个由无限高的墙包围的平坦势）很有趣，因为它允许有分析上已知的解决方案。然而，除此之外，它是相当无聊的。也就是说，通过在输入文件中设置不同的势来玩弄这个势是很容易的。例如，让我们假设我们想在2d中使用以下势。

@f[
  V(x,y) = \left\{
       \begin{array}{ll}


         -100 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy>0
         \\


         -5 & \text{if}\ \sqrt{x^2+y^2}<\frac 34 \ \text{and}
                         \ xy\le 0
         \\
         0 & \text{otherwise}
      \end{array} \right.\quad.


@f]

换句话说，在半径为0.75的圆的两个扇面中，电位为-100，在另外两个扇面中为-5，而在圆外为零。我们可以通过在输入文件中使用以下内容来实现这一点。

@code
set Potential = if (x^2 + y^2 < 0.75^2, if (x*y > 0, -100, -5), 0)
@endcode

此外，如果我们还将网格细化程度提高一级，我们会得到以下结果。

@code
examples/\step-36> make run
============================ Running \step-36
   Number of active cells:       4096
   Number of degrees of freedom: 4225


   Eigenvalue 0 : -74.2562
   Eigenvalue 1 : -72.7322
   Eigenvalue 2 : -42.7406
   Eigenvalue 3 : -42.2232
   Eigenvalue 4 : -37.0744
@endcode



输出文件还包含一个内插的势的版本，看起来像这样（注意，正如预期的那样，最低的几个特征模式的概率密度 $|\Psi(\mathbf x)|^2$ 只有在势最低的地方才是显著的，即在势的内圈的右上角和左下角部分）。

 <img src="https://www.dealii.org/images/steps/developer/step-36.mod.potential.png" alt=""> 

前五个特征函数现在是这样的。

 <table width="80%">
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.0.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.1.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.2.png" alt=""></td>
</tr>
<tr>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.3.png" alt=""></td>
<td><img src="https://www.dealii.org/images/steps/developer/step-36.mod.eigenfunction.4.png" alt=""></td>
<td></td>
</tr>
</table> 

 <li> 在我们对问题的推导中，我们假设粒子被限制在一个域 $\Omega$ 中，并且在这个域的边界处，它的概率 $|\Psi|^2$ 为零。这相当于解决所有 ${\mathbb R}^d$ 上的特征值问题，并假设能量势只在 $\Omega$ 区域内是有限的，而在区域外是无限的。比较容易的是，在 $|\Psi(\mathbf x)|^2$ 的所有位置 $\mathbf x$ ， $V(\mathbf
x)=\infty$ 。那么问题来了，如果我们的势不是这种形式的，即没有一个势是无限的有界域，会发生什么？在这种情况下，可能值得只考虑一个非常大的边界域，其中 $V(\mathbf x)$ 即使不是无限的，也至少是非常大的。在这样的情况下玩一玩，探索一下当我们使计算区域越来越大时，频谱和特征函数如何变化。

 <li>  如果我们研究简单的谐波振荡器问题 $V(\mathbf x)=c|\mathbf x|^2$ 会怎样？这个势正是上一段所讨论的形式，具有超球面对称性。人们可能想用一个大的外半径的球面域，来近似于全空间问题（例如，通过引用 GridGenerator::hyper_ball). 

 <li>  上面的图显示了波函数  $\Psi(\mathbf x)$  ，但感兴趣的物理量实际上是粒子处于位置  $|\Psi(\mathbf x)|^2$  的概率密度。一些可视化程序可以从输入文件中的数据计算出衍生量，但我们也可以在创建输出文件时立即这样做。这样做的工具是DataPostprocessor类，可以和DataOut类一起使用。如何做到这一点的例子可以在步骤29和步骤33中找到。

 <li>  如果盒子里的粒子有%的内部自由度会怎样？例如，如果该粒子是一个自旋-  $1/2$  粒子？在这种情况下，我们可能要开始解决一个矢量值的问题，而不是。

 <li>  我们这里的deal.II库的实现使用PETScWrappers和SLEPcWrappers，适合在串行机器架构上运行。然而，对于更大的网格和更多的自由度，我们可能希望在并行架构上运行我们的应用程序。上述代码的并行实现在这里可能特别有用，因为广义的特征谱问题比前面大多数教程中考虑的标准问题的解决成本更高。幸运的是，修改上述程序使其符合MPI标准是一个相对简单的过程。关于如何做到这一点的简图可以在 @ref
step_17 "step-17 "中找到。

deal.II有与其中之一ARPACK（见<a
href="../../external-libs/arpack.html">the ARPACK configuration page</a>的设置说明）的接口，在ArpackSolver类中实现。下面是一个简短的、快速的概述，说明使用它需要改变什么，前提是你有一个工作的ARPACK安装，并且deal.II已经为它正确配置了（见deal.II <a href="../../readme.html" target="body">README</a>文件）。

首先，为了使用ARPACK接口，我们可以回到使用标准的deal.II矩阵和向量，所以我们首先替换PETSc和SLEPc头文件

@code
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/slepc_solver.h>
@endcode

与这些。

@code
#include <deal.II/lac/arpack_solver.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
@endcode

ARPACK允许复杂的特征值，所以我们还需要

@code
#include <complex>
@endcode



其次，我们在主类中切换回deal.II矩阵和向量定义。

@code
    SparsityPattern                     sparsity_pattern;
    SparseMatrix<double>                stiffness_matrix, mass_matrix;
    std::vector<Vector<double> >        eigenfunctions;
    std::vector<std::complex<double>>   eigenvalues;
@endcode

并按照惯例对它们进行初始化  <code>make_grid_and_dofs()</code>  。

@code
    sparsity_pattern.reinit (dof_handler.n_dofs(),
                             dof_handler.n_dofs(),
                             dof_handler.max_couplings_between_dofs());


    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    constraints.condense (sparsity_pattern);
    sparsity_pattern.compress();


    stiffness_matrix.reinit (sparsity_pattern);
    mass_matrix.reinit (sparsity_pattern);
@endcode



为了用ARPACK解决特征值问题，我们最后需要修改  <code>solve()</code>  。

@code
  template <int dim>
  unsigned int EigenvalueProblem<dim>::solve ()
  {
    SolverControl solver_control (dof_handler.n_dofs(), 1e-9);


    SparseDirectUMFPACK inverse;
    inverse.initialize (stiffness_matrix);


    const unsigned int num_arnoldi_vectors = 2*eigenvalues.size() + 2;
    ArpackSolver::AdditionalData additional_data(num_arnoldi_vectors);


    ArpackSolver eigensolver (solver_control, additional_data);
    eigensolver.solve (stiffness_matrix,
                       mass_matrix,
                       inverse,
                       eigenvalues,
                       eigenfunctions,
                       eigenvalues.size());


    for (unsigned int i=0; i<eigenfunctions.size(); ++i)
      eigenfunctions[i] /= eigenfunctions[i].linfty_norm ();


    return solver_control.last_step ();
  }
@endcode

请注意我们是如何使用精确分解（使用SparseDirectUMFPACK）作为ARPACK的预处理程序的。   </ul> 


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-36.cc"
*/
