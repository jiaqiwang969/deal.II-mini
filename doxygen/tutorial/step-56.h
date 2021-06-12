/**
@page step_56 The step-56 tutorial program
This tutorial depends on step-16, step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#StokesProblem"> Stokes Problem </a>
        <li><a href="#LinearSolverandPreconditioningIssues"> Linear Solver and Preconditioning Issues </a>
        <li><a href="#ReferenceSolution"> Reference Solution </a>
        <li><a href="#ComputingErrors"> Computing Errors </a>
        <li><a href="#DoFHandlers"> DoF Handlers </a>
        <li><a href="#DifferencesfromtheStep22tutorial"> Differences from the Step 22 tutorial </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#FunctionsforSolutionandRighthandside">Functions for Solution and Righthand side</a>
        <li><a href="#ASPECTBlockSchurPreconditioner">ASPECT BlockSchurPreconditioner</a>
        <li><a href="#TheStokesProblemclass">The StokesProblem class</a>
      <ul>
        <li><a href="#StokesProblemsetup_dofs">StokesProblem::setup_dofs</a>
        <li><a href="#StokesProblemassemble_system">StokesProblem::assemble_system</a>
        <li><a href="#StokesProblemassemble_multigrid">StokesProblem::assemble_multigrid</a>
        <li><a href="#StokesProblemsolve">StokesProblem::solve</a>
        <li><a href="#StokesProblemprocess_solution">StokesProblem::process_solution</a>
        <li><a href="#StokesProblemoutput_results">StokesProblem::output_results</a>
        <li><a href="#StokesProblemrun">StokesProblem::run</a>
      </ul>
        <li><a href="#Themainfunction">The main function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Errors"> Errors </a>
        <li><a href="#TimingResults"> Timing Results </a>
        <li><a href="#Possibilitiesforextensions"> Possibilities for extensions </a>
      <ul>
        <li><a href="#Checkhigherorderdiscretizations"> Check higher order discretizations </a>
        <li><a href="#Comparewithcheappreconditioner"> Compare with cheap preconditioner </a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-56/doc/intro.dox

<i>This program was contributed by Ryan Grove and Timo Heister.


This material is based upon work partially supported by National Science
Foundation grant DMS1522191 and the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-0949446 and The University of California-Davis.


The authors would like to thank the Isaac Newton Institute for
Mathematical Sciences, Cambridge, for support and hospitality during
the programme Melt in the Mantle where work on this tutorial was
undertaken. This work was supported by EPSRC grant no EP/K032208/1.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.400995,https://zenodo.org/badge/DOI/10.5281/zenodo.400995.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="StokesProblem"></a><h3> Stokes Problem </h3>


本教程的目的是为斯托克斯方程创建一个高效的线性求解器，并将其与其他方法进行比较。  在这里，我们将使用带有几何多栅的FGMRES作为预处理速度块，我们将在结果部分显示，这比step-22中使用的线性求解器（包括 "可能的扩展 "中描述的方案）从根本上来说是一种更好的方法。  从根本上说，这是因为只有多网格才有可能得到 $O(n)$ 的求解时间，其中 $n$ 是线性系统的未知数的数量。使用Timer类，我们收集一些统计数据来比较设置时间、求解时间和迭代次数。我们还计算了误差，以确保我们所实现的是正确的。

让  $u \in H_0^1 = \{ u \in H^1(\Omega), u|_{\partial \Omega} = 0 \}$  和  $p \in L_*^2 = \{ p \in L^2(\Omega), \int_\Omega p = 0
\}$  。斯托克斯方程的非维度形式如下。

@f{eqnarray*}


 - 2 \text{div} \frac {1}{2} \left[ (\nabla \textbf{u})
 + (\nabla \textbf{u})^T\right] + \nabla p & =& f \\


 - \nabla \cdot u &=& 0


@f}



请注意，我们使用的是变形张量，而不是 $\Delta u$ （关于两者之间的区别的详细描述可以在步骤22中找到，但总的来说，变形张量的物理性更强，也更昂贵）。

<a name="LinearSolverandPreconditioningIssues"></a><h3> Linear %Solver and Preconditioning Issues </h3>


离散方程的微弱形式自然导致了以下速度场和压力场的节点值的线性系统。

@f{eqnarray*}
\left(\begin{array}{cc} A & B^T \\ B & 0
\end{array}\right) \left(\begin{array}{c} U \\ P \end{array}\right) =
\left(\begin{array}{c} F \\ 0 \end{array}\right).


@f}



我们的目标是比较几种解决方法。  虽然step-22使用 "Schur补足法 "分两步解决线性系统，但我们本着step-22的 "结果 "部分所概述的方法的精神，使用FMGRES和一个有效的预处理器一次性地攻击块系统。其思路如下：如果我们找到一个块状预处理程序 $P$ ，使矩阵的

@f{eqnarray*}
\left(\begin{array}{cc} A & B^T \\ B & 0 \end{array}\right) P^{-1}


@f}



是简单的，那么使用该预处理的迭代求解器将在几次迭代后收敛。请注意，我们在这里做的是正确的预处理。  使用舒尔补码 $S=BA^{-1}B^T$ ，我们发现

@f{eqnarray*}
P^{-1} = \left(\begin{array}{cc} A & B^T \\ 0 &
 S \end{array}\right)^{-1}


@f}



是一个很好的选择。设 $\widetilde{A^{-1}}$ 是 $A^{-1}$ 的近似值， $\widetilde{S^{-1}}$ 是 $S^{-1}$ 的近似值，我们看到

@f{eqnarray*}
P^{-1} =
\left(\begin{array}{cc} A^{-1} & 0 \\ 0 & I \end{array}\right)
\left(\begin{array}{cc} I & B^T \\ 0 & -I \end{array}\right)
\left(\begin{array}{cc} I & 0 \\ 0 & S^{-1} \end{array}\right)
\approx
\left(\begin{array}{cc} \widetilde{A^{-1}} & 0 \\ 0 & I \end{array}\right)
\left(\begin{array}{cc} I & B^T \\ 0 & -I \end{array}\right)
\left(\begin{array}{cc} I & 0 \\ 0 & \widetilde{S^{-1}} \end{array}\right).
  @f}



由于 $P$ 的目的只是作为一个预处理程序，我们将在上式中使用右边的近似值。

正如步骤22所讨论的， $-M_p^{-1}=:\widetilde{S^{-1}} \approx
S^{-1}$ ，其中 $M_p$ 是压力质量矩阵，通过使用CG与ILU作为预处理程序近似求解， $\widetilde{A^{-1}}$ 是通过多种方法之一得到的：使用CG和ILU作为预处理程序求解线性系统，仅仅使用ILU的一次应用，使用CG和GMG（步骤16中描述的几何多网格）作为预处理程序求解线性系统，或者仅仅执行GMG的一个V-循环。

作为比较，我们也在整个系统上使用直接求解器UMFPACK来比较我们的结果，而不是FGMRES。  如果你想使用直接求解器（如UMFPACK），系统需要是可逆的。为了避免恒定压力给定的一维无效空间，我们将第一个压力未知数固定为零。这对迭代求解器来说是没有必要的。




<a name="ReferenceSolution"></a><h3> Reference Solution </h3>


测试问题是一个 "制造的解决方案"（详见步骤7），我们选择  $u=(u_1,u_2,u_3)=(2\sin (\pi x), - \pi y \cos
(\pi x),- \pi z \cos (\pi x))$  和  $p = \sin (\pi x)\cos (\pi y)\sin
(\pi z)$  。我们在域的整个边界上对速度应用迪里切特边界条件  $\Omega=[0,1]\times[0,1]\times[0,1]$  。为了执行边界条件，我们可以直接使用我们的参考解。

如果你在deal.II手册中查找创建一个从 <code>Function@<dim@></code> 派生的类所需要的东西，你会发现这个类有许多 @p virtual 函数，包括 Function::value(),   Function::vector_value(),   Function::value_list(),  等，所有这些都可以被重载。  deal.II的不同部分将需要这些特定函数中的不同部分。这在一开始会让人感到困惑，但幸运的是，你真正需要实现的只有 @p value(). 。Function类中的其他虚拟函数在里面有默认的实现，会默认调用你对 @p value 的实现。

注意，我们的参考解满足 $\nabla \cdot u = 0$  。此外，压力被选择为平均值为零。  对于步骤7的 "制造解决方案的方法"，我们需要找到 $\bf
f$ ，以便。

@f{align*}
{\bf f} =   - 2 \text{div} \frac {1}{2} \left[ (\nabla \textbf{u}) + (\nabla \textbf{u})^T\right] + \nabla p.


@f}



使用上面的参考解，我们得到。

@f{eqnarray*}
{\bf f} &=& (2 \pi^2 \sin (\pi x),- \pi^3 y \cos(\pi
x),- \pi^3 z \cos(\pi x))\\ & & + (\pi \cos(\pi x) \cos(\pi y)
\sin(\pi z) ,- \pi \sin(\pi y) \sin(\pi x) \sin(\pi z), \pi \cos(\pi
z) \sin(\pi x) \cos(\pi y)) @f}



<a name="ComputingErrors"></a><h3> Computing Errors </h3>


因为我们在线性系统中没有强制要求平均压力为零，所以我们需要在求解后对解决方案进行后处理。为了做到这一点，我们使用 VectorTools::compute_mean_value() 函数来计算压力的平均值，以从压力中减去它。




<a name="DoFHandlers"></a><h3> DoF Handlers </h3>


我们在这里实现几何多网格的方式只对速度变量（即上面描述的 $A$ 矩阵）进行执行，而不是压力。我们可以用不同的方法来实现这一点，包括将所有的粗网格操作视为作用于 $2\times
2$ 块系统，我们只考虑左上方的块。另外，我们也可以通过真正只考虑整个有限元离散化的速度部分的线性系统来实现。后者是我们在这里想要使用的方式。

为了实现这一点，我们需要能够提出这样的问题："我可以只拥有一个DoFHandler的一部分吗？"。在编写这个程序的时候，这是不可能的，所以为了满足我们的需求，我们只是为速度创建一个单独的、第二个DoFHandler。然后，我们只基于这个第二DoFHandler为多网格预处理程序建立线性系统，并简单地将第一块（整体）向量转移到整个第二DoFHandler的对应向量中。要做到这一点，我们必须保证两个DoFHandler对象中的（速度）自由度排序的<i>order</i>是相同的。事实上，首先在两个对象上分配自由度，然后在两个对象上使用相同的DoFRenumbering操作序列，就可以做到这一点。




<a name="DifferencesfromtheStep22tutorial"></a><h3> Differences from the Step 22 tutorial </h3>


第56步和第22步的主要区别是我们使用了块状求解器，而不是第22步中使用的Schur补码方法。这种方法的细节可以在步骤-22的 "可能的扩展 "部分的 "块状Schur补码预处理 "小节中找到。对于速度块的预处理，我们从<a href="https://aspect.geodynamics.org">ASPECT</a>中借用了一个叫做 @p BlockSchurPreconditioner 的类，该类可以选择求解 $A$ 的逆，或者只对其应用一个预处理扫频，这分别为我们提供了一种昂贵和便宜的方法。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * <a name="Includefiles"></a> 
 * <h3>Include files</h3>
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
 * #include <deal.II/lac/block_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_gmres.h> 
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
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * #include <deal.II/lac/sparse_ilu.h> 
 * #include <deal.II/grid/grid_out.h> 
 * 
 * @endcode
 * 
 * 我们需要包括以下文件来做计时。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/timer.h> 
 * 
 * @endcode
 * 
 * 这包括我们使用几何多网格所需的文件
 * 

 * 
 * 
 * @code
 * #include <deal.II/multigrid/multigrid.h> 
 * #include <deal.II/multigrid/mg_transfer.h> 
 * #include <deal.II/multigrid/mg_tools.h> 
 * #include <deal.II/multigrid/mg_coarse.h> 
 * #include <deal.II/multigrid/mg_smoother.h> 
 * #include <deal.II/multigrid/mg_matrix.h> 
 * 
 * #include <iostream> 
 * #include <fstream> 
 * 
 * namespace Step56 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 为了便于在所使用的不同求解器之间进行切换，我们声明了一个枚举，可以作为参数传递给主类的构造函数。
 * 

 * 
 * 
 * @code
 *   enum class SolverType 
 *   { 
 *     FGMRES_ILU, 
 *     FGMRES_GMG, 
 *     UMFPACK 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="FunctionsforSolutionandRighthandside"></a> 
 * <h3>Functions for Solution and Righthand side</h3>
 * 

 * 
 * Solution类用于定义边界条件和计算数值解的误差。请注意，我们需要定义数值和梯度，以便计算L2和H1误差。在这里，我们决定使用模板的特殊化来分离2D和3D的实现。
 * 

 * 
 * 请注意，前几个分量是速度分量，最后一个分量是压力。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class Solution : public Function<dim> 
 *   { 
 *   public: 
 *     Solution() 
 *       : Function<dim>(dim + 1) 
 *     {} 
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component = 0) const override; 
 *     virtual Tensor<1, dim> 
 *     gradient(const Point<dim> & p, 
 *              const unsigned int component = 0) const override; 
 *   }; 
 * 
 *   template <> 
 *   double Solution<2>::value(const Point<2> &   p, 
 *                             const unsigned int component) const 
 *   { 
 *     Assert(component <= 2 + 1, ExcIndexRange(component, 0, 2 + 1)); 
 * 
 *     using numbers::PI; 
 *     const double x = p(0); 
 *     const double y = p(1); 
 * 
 *     if (component == 0) 
 *       return sin(PI * x); 
 *     if (component == 1) 
 *       return -PI * y * cos(PI * x); 
 *     if (component == 2) 
 *       return sin(PI * x) * cos(PI * y); 
 * 
 *     return 0; 
 *   } 
 * 
 *   template <> 
 *   double Solution<3>::value(const Point<3> &   p, 
 *                             const unsigned int component) const 
 *   { 
 *     Assert(component <= 3 + 1, ExcIndexRange(component, 0, 3 + 1)); 
 * 
 *     using numbers::PI; 
 *     const double x = p(0); 
 *     const double y = p(1); 
 *     const double z = p(2); 
 * 
 *     if (component == 0) 
 *       return 2.0 * sin(PI * x); 
 *     if (component == 1) 
 *       return -PI * y * cos(PI * x); 
 *     if (component == 2) 
 *       return -PI * z * cos(PI * x); 
 *     if (component == 3) 
 *       return sin(PI * x) * cos(PI * y) * sin(PI * z); 
 * 
 *     return 0; 
 *   } 
 * 
 * @endcode
 * 
 * 注意，对于梯度，我们需要返回一个Tensor<1,dim>。
 * 

 * 
 * 
 * @code
 *   template <> 
 *   Tensor<1, 2> Solution<2>::gradient(const Point<2> &   p, 
 *                                      const unsigned int component) const 
 *   { 
 *     Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1)); 
 * 
 *     using numbers::PI; 
 *     const double x = p(0); 
 *     const double y = p(1); 
 * 
 *     Tensor<1, 2> return_value; 
 *     if (component == 0) 
 *       { 
 *         return_value[0] = PI * cos(PI * x); 
 *         return_value[1] = 0.0; 
 *       } 
 *     else if (component == 1) 
 *       { 
 *         return_value[0] = y * PI * PI * sin(PI * x); 
 *         return_value[1] = -PI * cos(PI * x); 
 *       } 
 *     else if (component == 2) 
 *       { 
 *         return_value[0] = PI * cos(PI * x) * cos(PI * y); 
 *         return_value[1] = -PI * sin(PI * x) * sin(PI * y); 
 *       } 
 * 
 *     return return_value; 
 *   } 
 * 
 *   template <> 
 *   Tensor<1, 3> Solution<3>::gradient(const Point<3> &   p, 
 *                                      const unsigned int component) const 
 *   { 
 *     Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1)); 
 * 
 *     using numbers::PI; 
 *     const double x = p(0); 
 *     const double y = p(1); 
 *     const double z = p(2); 
 * 
 *     Tensor<1, 3> return_value; 
 *     if (component == 0) 
 *       { 
 *         return_value[0] = 2 * PI * cos(PI * x); 
 *         return_value[1] = 0.0; 
 *         return_value[2] = 0.0; 
 *       } 
 *     else if (component == 1) 
 *       { 
 *         return_value[0] = y * PI * PI * sin(PI * x); 
 *         return_value[1] = -PI * cos(PI * x); 
 *         return_value[2] = 0.0; 
 *       } 
 *     else if (component == 2) 
 *       { 
 *         return_value[0] = z * PI * PI * sin(PI * x); 
 *         return_value[1] = 0.0; 
 *         return_value[2] = -PI * cos(PI * x); 
 *       } 
 *     else if (component == 3) 
 *       { 
 *         return_value[0] = PI * cos(PI * x) * cos(PI * y) * sin(PI * z); 
 *         return_value[1] = -PI * sin(PI * x) * sin(PI * y) * sin(PI * z); 
 *         return_value[2] = PI * sin(PI * x) * cos(PI * y) * cos(PI * z); 
 *       } 
 * 
 *     return return_value; 
 *   } 
 * 
 * @endcode
 * 
 * 实现  $f$  。更多信息请参见介绍。
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
 *   }; 
 * 
 *   template <> 
 *   double RightHandSide<2>::value(const Point<2> &   p, 
 *                                  const unsigned int component) const 
 *   { 
 *     Assert(component <= 2, ExcIndexRange(component, 0, 2 + 1)); 
 * 
 *     using numbers::PI; 
 *     double x = p(0); 
 *     double y = p(1); 
 *     if (component == 0) 
 *       return PI * PI * sin(PI * x) + PI * cos(PI * x) * cos(PI * y); 
 *     if (component == 1) 
 *       return -PI * PI * PI * y * cos(PI * x) - PI * sin(PI * y) * sin(PI * x); 
 *     if (component == 2) 
 *       return 0; 
 * 
 *     return 0; 
 *   } 
 * 
 *   template <> 
 *   double RightHandSide<3>::value(const Point<3> &   p, 
 *                                  const unsigned int component) const 
 *   { 
 *     Assert(component <= 3, ExcIndexRange(component, 0, 3 + 1)); 
 * 
 *     using numbers::PI; 
 *     double x = p(0); 
 *     double y = p(1); 
 *     double z = p(2); 
 *     if (component == 0) 
 *       return 2 * PI * PI * sin(PI * x) + 
 *              PI * cos(PI * x) * cos(PI * y) * sin(PI * z); 
 *     if (component == 1) 
 *       return -PI * PI * PI * y * cos(PI * x) + 
 *              PI * (-1) * sin(PI * y) * sin(PI * x) * sin(PI * z); 
 *     if (component == 2) 
 *       return -PI * PI * PI * z * cos(PI * x) + 
 *              PI * cos(PI * z) * sin(PI * x) * cos(PI * y); 
 *     if (component == 3) 
 *       return 0; 
 * 
 *     return 0; 
 *   } 
 * 
 * @endcode
 * 
 * 
 * <a name="ASPECTBlockSchurPreconditioner"></a> 
 * <h3>ASPECT BlockSchurPreconditioner</h3>
 * 

 * 
 * 在下文中，我们将实现一个预处理程序，它扩展了  step-22  的结果部分所讨论的想法。具体来说，我们1.使用一个上块三角的预处理器，因为我们想使用右预处理。2.可选择允许对速度块使用内部求解器，而不是使用单一的预处理程序。3.不使用InverseMatrix，而是明确地调用SolverCG。这种方法也用于ASPECT代码（见https:aspect.geodynamics.org），该代码在模拟地幔对流的背景下求解斯托克斯方程，该代码已被用于解决成千上万个处理器上的问题。
 * 

 * 
 * 构造函数中的bool标志 @p do_solve_A 允许我们对速度块应用一次预处理，或者使用内部迭代求解器来代替更精确的近似。
 * 

 * 
 * 注意我们是如何跟踪内部迭代的总和（预处理程序的应用）的。
 * 

 * 
 * 
 * @code
 *   template <class PreconditionerAType, class PreconditionerSType> 
 *   class BlockSchurPreconditioner : public Subscriptor 
 *   { 
 *   public: 
 *     BlockSchurPreconditioner( 
 *       const BlockSparseMatrix<double> &system_matrix, 
 *       const SparseMatrix<double> &     schur_complement_matrix, 
 *       const PreconditionerAType &      preconditioner_A, 
 *       const PreconditionerSType &      preconditioner_S, 
 *       const bool                       do_solve_A);
 * 
 *     void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const; 
 * 
 *     mutable unsigned int n_iterations_A; 
 *     mutable unsigned int n_iterations_S; 
 * 
 *   private: 
 *     const BlockSparseMatrix<double> &system_matrix; 
 *     const SparseMatrix<double> &     schur_complement_matrix; 
 *     const PreconditionerAType &      preconditioner_A; 
 *     const PreconditionerSType &      preconditioner_S; 
 * 
 *     const bool do_solve_A; 
 *   }; 
 * 
 *   template <class PreconditionerAType, class PreconditionerSType> 
 *   BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>:: 
 *     BlockSchurPreconditioner( 
 *       const BlockSparseMatrix<double> &system_matrix, 
 *       const SparseMatrix<double> &     schur_complement_matrix, 
 *       const PreconditionerAType &      preconditioner_A, 
 *       const PreconditionerSType &      preconditioner_S, 
 *       const bool                       do_solve_A) 
 *     : n_iterations_A(0)  
 *     , n_iterations_S(0) 
 *     , system_matrix(system_matrix) 
 *     , schur_complement_matrix(schur_complement_matrix) 
 *     , preconditioner_A(preconditioner_A) 
 *     , preconditioner_S(preconditioner_S) 
 *     , do_solve_A(do_solve_A) 
 *   {} 
 * 
 *   template <class PreconditionerAType, class PreconditionerSType> 
 *   void 
 *   BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::vmult( 
 *     BlockVector<double> &      dst, 
 *     const BlockVector<double> &src) const 
 *   { 
 *     Vector<double> utmp(src.block(0)); 
 * 
 * @endcode
 * 
 * 首先用S的近似值求解
 * 

 * 
 * 
 * @code
 *     { 
 *       SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm()); 
 *       SolverCG<Vector<double>> cg(solver_control); 
 * 
 *       dst.block(1) = 0.0; 
 *       cg.solve(schur_complement_matrix, 
 *                dst.block(1), 
 *                src.block(1), 
 *                preconditioner_S); 
 * 
 *       n_iterations_S += solver_control.last_step(); 
 *       dst.block(1) *= -1.0; 
 *     } 
 * 
 * @endcode
 * 
 * 第二，应用右上方的块（B^T
 * 

 * 
 * 
 * @code
 *     { 
 *       system_matrix.block(0, 1).vmult(utmp, dst.block(1)); 
 *       utmp *= -1.0; 
 *       utmp += src.block(0); 
 *     } 
 * 
 * @endcode
 * 
 * 最后，要么用左上角的块求解，要么只应用一个预设条件器扫频
 * 

 * 
 * 
 * @code
 *     if (do_solve_A == true) 
 *       { 
 *         SolverControl            solver_control(10000, utmp.l2_norm() * 1e-4); 
 *         SolverCG<Vector<double>> cg(solver_control); 
 * 
 *         dst.block(0) = 0.0; 
 *         cg.solve(system_matrix.block(0, 0), 
 *                  dst.block(0), 
 *                  utmp, 
 *                  preconditioner_A); 
 * 
 *         n_iterations_A += solver_control.last_step(); 
 *       } 
 *     else 
 *       { 
 *         preconditioner_A.vmult(dst.block(0), utmp); 
 *         n_iterations_A += 1; 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="TheStokesProblemclass"></a> 
 * <h3>The StokesProblem class</h3>
 * 

 * 
 * 这是该问题的主类。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class StokesProblem 
 *   { 
 *   public: 
 *     StokesProblem(const unsigned int pressure_degree, 
 *                   const SolverType   solver_type); 
 *     void run(); 
 * 
 *   private: 
 *     void setup_dofs(); 
 *     void assemble_system(); 
 *     void assemble_multigrid(); 
 *     void solve(); 
 *     void compute_errors(); 
 *     void output_results(const unsigned int refinement_cycle) const; 
 * 
 *     const unsigned int pressure_degree; 
 *     const SolverType   solver_type; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FESystem<dim>      velocity_fe; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 *     DoFHandler<dim>    velocity_dof_handler; 
 * 
 *     AffineConstraints<double> constraints; 
 * 
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> system_matrix; 
 *     SparseMatrix<double>      pressure_mass_matrix; 
 * 
 *     BlockVector<double> solution; 
 *     BlockVector<double> system_rhs; 
 * 
 *     MGLevelObject<SparsityPattern>      mg_sparsity_patterns; 
 *     MGLevelObject<SparseMatrix<double>> mg_matrices; 
 *     MGLevelObject<SparseMatrix<double>> mg_interface_matrices; 
 *     MGConstrainedDoFs                   mg_constrained_dofs; 
 * 
 *     TimerOutput computing_timer; 
 *   }; 
 * 
 *   template <int dim> 
 *   StokesProblem<dim>::StokesProblem(const unsigned int pressure_degree, 
 *                                     const SolverType   solver_type) 
 * 
 *     : pressure_degree(pressure_degree) 
 *     , solver_type(solver_type) 
 *     , triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , 
 * 
 * @endcode
 * 
 * 仅为速度的有限元。
 * 

 * 
 * 
 * @code
 *     velocity_fe(FE_Q<dim>(pressure_degree + 1), dim) 
 *     , 
 * 
 * @endcode
 * 
 * 整个系统的有限元。
 * 

 * 
 * 
 * @code
 *     fe(velocity_fe, 1, FE_Q<dim>(pressure_degree), 1) 
 *     , dof_handler(triangulation) 
 *     , velocity_dof_handler(triangulation) 
 *     , computing_timer(std::cout, TimerOutput::never, TimerOutput::wall_times) 
 *   {} 
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsetup_dofs"></a> 
 * <h4>StokesProblem::setup_dofs</h4>
 * 

 * 
 * 这个函数设置了DoFHandler、矩阵、向量和Multigrid结构（如果需要）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::setup_dofs() 
 *   { 
 *     TimerOutput::Scope scope(computing_timer, "Setup"); 
 * 
 *     system_matrix.clear(); 
 *     pressure_mass_matrix.clear(); 
 * 
 * @endcode
 * 
 * 主DoFHandler只需要活动的DoF，所以我们不在这里调用distribution_mg_dofs()
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 这个块结构将dim速度分量与压力分量（用于重新排序）分开。注意，我们有2个而不是像 step-22 中的dim+1块，因为我们的FESystem是嵌套的，dim速度分量作为一个块出现。
 * 

 * 
 * 
 * @code
 *     std::vector<unsigned int> block_component(2); 
 *     block_component[0] = 0; 
 *     block_component[1] = 1; 
 * 
 * @endcode
 * 
 * 速度从组件0开始。
 * 

 * 
 * 
 * @code
 *     const FEValuesExtractors::Vector velocities(0); 
 * 
 * @endcode
 * 
 * 如果我们应用重新排序来减少填充，
 * ILU的表现会更好。对于其他求解器来说，这样做并没有什么好处。
 * 

 * 
 * 
 * @code
 *     if (solver_type == SolverType::FGMRES_ILU) 
 *       { 
 *         TimerOutput::Scope ilu_specific(computing_timer, "(ILU specific)"); 
 *         DoFRenumbering::Cuthill_McKee(dof_handler); 
 *       } 
 * 
 * @endcode
 * 
 * 这确保了所有的速度DoFs在压力未知数之前被列举出来。这允许我们使用块来处理向量和矩阵，并允许我们为dof_handler和velocity_dof_handler获得相同的DoF编号。
 * 

 * 
 * 
 * @code
 *     DoFRenumbering::block_wise(dof_handler); 
 * 
 *     if (solver_type == SolverType::FGMRES_GMG) 
 *       { 
 *         TimerOutput::Scope multigrid_specific(computing_timer, 
 *                                               "(Multigrid specific)"); 
 *         TimerOutput::Scope setup_multigrid(computing_timer, 
 *                                            "Setup - Multigrid"); 
 * 
 * @endcode
 * 
 * 这将在一个单独的DoFHandler中分配速度空间的主动道夫和多网格道夫，如介绍中所述。
 * 

 * 
 * 
 * @code
 *         velocity_dof_handler.distribute_dofs(velocity_fe); 
 *         velocity_dof_handler.distribute_mg_dofs(); 
 * 
 * @endcode
 * 
 * 下面的代码块初始化了MGConstrainedDofs（使用速度的边界条件），以及每个层次的稀疏模式和矩阵。MGLevelObject<T>的resize()函数将破坏所有现有的包含对象。
 * 

 * 
 * 
 * @code
 *         std::set<types::boundary_id> zero_boundary_ids; 
 *         zero_boundary_ids.insert(0); 
 * 
 *         mg_constrained_dofs.clear(); 
 *         mg_constrained_dofs.initialize(velocity_dof_handler); 
 *         mg_constrained_dofs.make_zero_boundary_constraints(velocity_dof_handler, 
 *                                                            zero_boundary_ids); 
 *         const unsigned int n_levels = triangulation.n_levels(); 
 * 
 *         mg_interface_matrices.resize(0, n_levels - 1); 
 *         mg_matrices.resize(0, n_levels - 1); 
 *         mg_sparsity_patterns.resize(0, n_levels - 1); 
 * 
 *         for (unsigned int level = 0; level < n_levels; ++level) 
 *           { 
 *             DynamicSparsityPattern csp(velocity_dof_handler.n_dofs(level), 
 *                                        velocity_dof_handler.n_dofs(level)); 
 *             MGTools::make_sparsity_pattern(velocity_dof_handler, csp, level); 
 *             mg_sparsity_patterns[level].copy_from(csp); 
 * 
 *             mg_matrices[level].reinit(mg_sparsity_patterns[level]); 
 *             mg_interface_matrices[level].reinit(mg_sparsity_patterns[level]); 
 *           } 
 *       } 
 * 
 *     const std::vector<types::global_dof_index> dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 *     const unsigned int n_u = dofs_per_block[0]; 
 *     const unsigned int n_p = dofs_per_block[1]; 
 * 
 *     { 
 *       constraints.clear(); 
 * 
 * @endcode
 * 
 * 下面利用分量掩码对速度的边界值进行插值，这在矢量值dealii  step-20 教程中进一步说明。
 * 

 * 
 * 
 * @code
 *       DoFTools::make_hanging_node_constraints(dof_handler, constraints); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                Solution<dim>(), 
 *                                                constraints, 
 *                                                fe.component_mask(velocities)); 
 * 
 * @endcode
 * 
 * 正如在介绍中所讨论的，我们需要固定压力变量的一个自由度以确保问题的可解性。在这里，我们将第一个压力自由度标记为受限自由度，该自由度的索引为n_u。
 * 

 * 
 * 
 * @code
 *       if (solver_type == SolverType::UMFPACK) 
 *         constraints.add_line(n_u); 
 * 
 *       constraints.close(); 
 *     } 
 * 
 *     std::cout << "\tNumber of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "\tNumber of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (" << n_u << '+' << n_p << ')' << std::endl; 
 * 
 *     { 
 *       BlockDynamicSparsityPattern csp(dofs_per_block, dofs_per_block); 
 *       DoFTools::make_sparsity_pattern(dof_handler, csp, constraints, false); 
 *       sparsity_pattern.copy_from(csp); 
 *     } 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     solution.reinit(dofs_per_block); 
 *     system_rhs.reinit(dofs_per_block); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemassemble_system"></a> 
 * <h4>StokesProblem::assemble_system</h4>
 * 

 * 
 * 在这个函数中，系统矩阵被组装起来。我们在(1,1)块中组装压力质量矩阵（如果需要），并在此函数结束时将其移出此位置。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::assemble_system() 
 *   { 
 *     TimerOutput::Scope assemble(computing_timer, "Assemble"); 
 *     system_matrix = 0; 
 *     system_rhs    = 0; 
 * 
 * @endcode
 * 
 * 如果为真，我们将在(1,1)块中装配压力质量矩阵。
 * 

 * 
 * 
 * @code
 *     const bool assemble_pressure_mass_matrix = 
 *       (solver_type == SolverType::UMFPACK) ? false : true; 
 * 
 *     QGauss<dim> quadrature_formula(pressure_degree + 2); 
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
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const RightHandSide<dim>    right_hand_side; 
 *     std::vector<Vector<double>> rhs_values(n_q_points, Vector<double>(dim + 1)); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell); 
 *     std::vector<double>                  div_phi_u(dofs_per_cell); 
 *     std::vector<double>                  phi_p(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
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
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               { 
 *                 for (unsigned int j = 0; j <= i; ++j) 
 *                   { 
 *                     local_matrix(i, j) += 
 *                       (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) - 
 *                        div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] + 
 *                        (assemble_pressure_mass_matrix ? phi_p[i] * phi_p[j] : 
 *                                                         0)) * 
 *                       fe_values.JxW(q); 
 *                   } 
 * 
 *                 const unsigned int component_i = 
 *                   fe.system_to_component_index(i).first; 
 *                 local_rhs(i) += fe_values.shape_value(i, q) * 
 *                                 rhs_values[q](component_i) * fe_values.JxW(q); 
 *               } 
 *           } 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
 *             local_matrix(i, j) = local_matrix(j, i); 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 *         constraints.distribute_local_to_global(local_matrix, 
 *                                                local_rhs, 
 *                                                local_dof_indices, 
 *                                                system_matrix, 
 *                                                system_rhs); 
 *       } 
 * 
 *     if (solver_type != SolverType::UMFPACK) 
 *       { 
 *         pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1)); 
 *         pressure_mass_matrix.copy_from(system_matrix.block(1, 1)); 
 *         system_matrix.block(1, 1) = 0; 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemassemble_multigrid"></a> 
 * <h4>StokesProblem::assemble_multigrid</h4>
 * 

 * 
 * 在这里，与 step-16 中一样，我们有一个函数，用于组装多网格预处理程序所需的水平矩阵和界面矩阵。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::assemble_multigrid() 
 *   { 
 *     TimerOutput::Scope multigrid_specific(computing_timer, 
 *                                           "(Multigrid specific)"); 
 *     TimerOutput::Scope assemble_multigrid(computing_timer, 
 *                                           "Assemble Multigrid"); 
 * 
 *     mg_matrices = 0.; 
 * 
 *     QGauss<dim> quadrature_formula(pressure_degree + 2); 
 * 
 *     FEValues<dim> fe_values(velocity_fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values | update_gradients); 
 * 
 *     const unsigned int dofs_per_cell = velocity_fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 * 
 *     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell); 
 * 
 *     std::vector<AffineConstraints<double>> boundary_constraints( 
 *       triangulation.n_levels()); 
 *     std::vector<AffineConstraints<double>> boundary_interface_constraints( 
 *       triangulation.n_levels()); 
 *     for (unsigned int level = 0; level < triangulation.n_levels(); ++level) 
 *       { 
 *         boundary_constraints[level].add_lines( 
 *           mg_constrained_dofs.get_refinement_edge_indices(level)); 
 *         boundary_constraints[level].add_lines( 
 *           mg_constrained_dofs.get_boundary_indices(level)); 
 *         boundary_constraints[level].close(); 
 * 
 *         IndexSet idx = mg_constrained_dofs.get_refinement_edge_indices(level) & 
 *                        mg_constrained_dofs.get_boundary_indices(level); 
 * 
 *         boundary_interface_constraints[level].add_lines(idx); 
 *         boundary_interface_constraints[level].close(); 
 *       } 
 * 
 * @endcode
 * 
 * 这个迭代器会覆盖所有的单元格（不仅仅是活动的）。
 * 

 * 
 * 
 * @code
 *     for (const auto &cell : velocity_dof_handler.cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 *         cell_matrix = 0; 
 * 
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient(k, q); 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               for (unsigned int j = 0; j <= i; ++j) 
 *                 { 
 *                   cell_matrix(i, j) += 
 *                     (symgrad_phi_u[i] * symgrad_phi_u[j]) * fe_values.JxW(q); 
 *                 } 
 *           } 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = i + 1; j < dofs_per_cell; ++j) 
 *             cell_matrix(i, j) = cell_matrix(j, i); 
 * 
 *         cell->get_mg_dof_indices(local_dof_indices); 
 * 
 *         boundary_constraints[cell->level()].distribute_local_to_global( 
 *           cell_matrix, local_dof_indices, mg_matrices[cell->level()]); 
 * 
 *         for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *           for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *             if (!mg_constrained_dofs.at_refinement_edge(cell->level(), 
 *                                                         local_dof_indices[i]) || 
 *                 mg_constrained_dofs.at_refinement_edge(cell->level(), 
 *                                                        local_dof_indices[j])) 
 *               cell_matrix(i, j) = 0; 
 * 
 *         boundary_interface_constraints[cell->level()] 
 *           .distribute_local_to_global(cell_matrix, 
 *                                       local_dof_indices, 
 *                                       mg_interface_matrices[cell->level()]); 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemsolve"></a> 
 * <h4>StokesProblem::solve</h4>
 * 

 * 
 * 这个函数根据你想使用ILU或GMG作为预处理程序的情况进行不同的设置。 这两种方法共享相同的求解器（FGMRES），但需要初始化不同的预处理器。在这里，我们不仅为整个求解函数计时，还为预处理程序的设置以及求解本身分别计时。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::solve() 
 *   { 
 *     TimerOutput::Scope solve(computing_timer, "Solve"); 
 *     constraints.set_zero(solution); 
 * 
 *     if (solver_type == SolverType::UMFPACK) 
 *       { 
 *         computing_timer.enter_subsection("(UMFPACK specific)"); 
 *         computing_timer.enter_subsection("Solve - Initialize"); 
 * 
 *         SparseDirectUMFPACK A_direct; 
 *         A_direct.initialize(system_matrix); 
 * 
 *         computing_timer.leave_subsection(); 
 *         computing_timer.leave_subsection(); 
 * 
 *         { 
 *           TimerOutput::Scope solve_backslash(computing_timer, 
 *                                              "Solve - Backslash"); 
 *           A_direct.vmult(solution, system_rhs); 
 *         } 
 * 
 *         constraints.distribute(solution); 
 *         return; 
 *       } 
 * 
 * @endcode
 * 
 * 这里我们必须确保以 "足够好 "的精度求解残差
 * 

 * 
 * 
 * @code
 *     SolverControl solver_control(system_matrix.m(), 
 *                                  1e-10 * system_rhs.l2_norm()); 
 *     unsigned int  n_iterations_A; 
 *     unsigned int  n_iterations_S; 
 * 
 * @endcode
 * 
 * 这是用来传递我们是否要在预处理程序中解决A的问题。 我们可以把它改为false，看看是否还能收敛，如果能收敛，那么程序的运行速度是快是慢？
 * 

 * 
 * 
 * @code
 *     const bool use_expensive = true; 
 * 
 *     SolverFGMRES<BlockVector<double>> solver(solver_control); 
 * 
 *     if (solver_type == SolverType::FGMRES_ILU) 
 *       { 
 *         computing_timer.enter_subsection("(ILU specific)"); 
 *         computing_timer.enter_subsection("Solve - Set-up Preconditioner"); 
 * 
 *         std::cout << "   Computing preconditioner..." << std::endl 
 *                   << std::flush; 
 * 
 *         SparseILU<double> A_preconditioner; 
 *         A_preconditioner.initialize(system_matrix.block(0, 0)); 
 * 
 *         SparseILU<double> S_preconditioner; 
 *         S_preconditioner.initialize(pressure_mass_matrix); 
 * 
 *         const BlockSchurPreconditioner<SparseILU<double>, SparseILU<double>> 
 *           preconditioner(system_matrix, 
 *                          pressure_mass_matrix, 
 *                          A_preconditioner, 
 *                          S_preconditioner, 
 *                          use_expensive); 
 * 
 *         computing_timer.leave_subsection(); 
 *         computing_timer.leave_subsection(); 
 * 
 *         { 
 *           TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES"); 
 * 
 *           solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 *           n_iterations_A = preconditioner.n_iterations_A; 
 *           n_iterations_S = preconditioner.n_iterations_S; 
 *         } 
 *       } 
 *     else 
 *       { 
 *         computing_timer.enter_subsection("(Multigrid specific)"); 
 *         computing_timer.enter_subsection("Solve - Set-up Preconditioner"); 
 * 
 * @endcode
 * 
 * 在各级之间转移运算符
 * 

 * 
 * 
 * @code
 *         MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs); 
 *         mg_transfer.build(velocity_dof_handler); 
 * 
 * @endcode
 * 
 * 设置粗略的网格解算器
 * 

 * 
 * 
 * @code
 *         FullMatrix<double> coarse_matrix; 
 *         coarse_matrix.copy_from(mg_matrices[0]); 
 *         MGCoarseGridHouseholder<double, Vector<double>> coarse_grid_solver; 
 *         coarse_grid_solver.initialize(coarse_matrix); 
 * 
 *         using Smoother = PreconditionSOR<SparseMatrix<double>>; 
 *         mg::SmootherRelaxation<Smoother, Vector<double>> mg_smoother; 
 *         mg_smoother.initialize(mg_matrices); 
 *         mg_smoother.set_steps(2); 
 * 
 * @endcode
 * 
 * Multigrid作为CG的预处理程序时，需要是一个对称的运算器，所以平滑器必须是对称的
 * 

 * 
 * 
 * @code
 *         mg_smoother.set_symmetric(true); 
 * 
 *         mg::Matrix<Vector<double>> mg_matrix(mg_matrices); 
 *         mg::Matrix<Vector<double>> mg_interface_up(mg_interface_matrices); 
 *         mg::Matrix<Vector<double>> mg_interface_down(mg_interface_matrices); 
 * 
 * @endcode
 * 
 * 现在，我们准备设置V型循环算子和多级预处理程序。
 * 

 * 
 * 
 * @code
 *         Multigrid<Vector<double>> mg( 
 *           mg_matrix, coarse_grid_solver, mg_transfer, mg_smoother, mg_smoother); 
 *         mg.set_edge_matrices(mg_interface_down, mg_interface_up); 
 * 
 *         PreconditionMG<dim, Vector<double>, MGTransferPrebuilt<Vector<double>>> 
 *           A_Multigrid(velocity_dof_handler, mg, mg_transfer); 
 * 
 *         SparseILU<double> S_preconditioner; 
 *         S_preconditioner.initialize(pressure_mass_matrix, 
 *                                     SparseILU<double>::AdditionalData()); 
 * 
 *         const BlockSchurPreconditioner< 
 *           PreconditionMG<dim, 
 *                          Vector<double>, 
 *                          MGTransferPrebuilt<Vector<double>>>, 
 *           SparseILU<double>> 
 *           preconditioner(system_matrix, 
 *                          pressure_mass_matrix, 
 *                          A_Multigrid, 
 *                          S_preconditioner, 
 *                          use_expensive); 
 * 
 *         computing_timer.leave_subsection(); 
 *         computing_timer.leave_subsection(); 
 * 
 *         { 
 *           TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES"); 
 *           solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 *           n_iterations_A = preconditioner.n_iterations_A; 
 *           n_iterations_S = preconditioner.n_iterations_S; 
 *         } 
 *       } 
 * 
 *     constraints.distribute(solution); 
 * 
 *     std::cout 
 *       << std::endl 
 *       << "\tNumber of FGMRES iterations: " << solver_control.last_step() 
 *       << std::endl 
 *       << "\tTotal number of iterations used for approximation of A inverse: " 
 *       << n_iterations_A << std::endl 
 *       << "\tTotal number of iterations used for approximation of S inverse: " 
 *       << n_iterations_S << std::endl 
 *       << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemprocess_solution"></a> 
 * <h4>StokesProblem::process_solution</h4>
 * 

 * 
 * 这个函数计算出解决方案的L2和H1误差。为此，我们需要确保压力的平均值为零。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::compute_errors() 
 *   { 
 * 
 * @endcode
 * 
 * 计算平均压力 $\frac{1}{\Omega} \int_{\Omega} p(x) dx $ ，然后从每个压力系数中减去它。这将产生一个平均值为零的压力。这里我们利用了压力是分量 $dim$ 和有限元空间是结点的事实。
 * 

 * 
 * 
 * @code
 *     const double mean_pressure = VectorTools::compute_mean_value( 
 *       dof_handler, QGauss<dim>(pressure_degree + 2), solution, dim); 
 *     solution.block(1).add(-mean_pressure); 
 *     std::cout << "   Note: The mean value was adjusted by " << -mean_pressure 
 *               << std::endl; 
 * 
 *     const ComponentSelectFunction<dim> pressure_mask(dim, dim + 1); 
 *     const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), 
 *                                                      dim + 1); 
 * 
 *     Vector<float> difference_per_cell(triangulation.n_active_cells()); 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(pressure_degree + 2), 
 *                                       VectorTools::L2_norm, 
 *                                       &velocity_mask); 
 * 
 *     const double Velocity_L2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(pressure_degree + 2), 
 *                                       VectorTools::L2_norm, 
 *                                       &pressure_mask); 
 * 
 *     const double Pressure_L2_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::L2_norm); 
 * 
 *     VectorTools::integrate_difference(dof_handler, 
 *                                       solution, 
 *                                       Solution<dim>(), 
 *                                       difference_per_cell, 
 *                                       QGauss<dim>(pressure_degree + 2), 
 *                                       VectorTools::H1_norm, 
 *                                       &velocity_mask); 
 * 
 *     const double Velocity_H1_error = 
 *       VectorTools::compute_global_error(triangulation, 
 *                                         difference_per_cell, 
 *                                         VectorTools::H1_norm); 
 * 
 *     std::cout << std::endl 
 *               << "   Velocity L2 Error: " << Velocity_L2_error << std::endl 
 *               << "   Pressure L2 Error: " << Pressure_L2_error << std::endl 
 *               << "   Velocity H1 Error: " << Velocity_H1_error << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StokesProblemoutput_results"></a> 
 * <h4>StokesProblem::output_results</h4>
 * 

 * 
 * 这个函数生成图形输出，就像在  step-22  中所做的那样。
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
 * 
 * @endcode
 * 
 * 
 * <a name="StokesProblemrun"></a> 
 * <h4>StokesProblem::run</h4>
 * 

 * 
 * 斯托克斯类的最后一步是像往常一样，生成初始网格的函数，并按各自的顺序调用其他函数。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StokesProblem<dim>::run() 
 *   { 
 *     GridGenerator::hyper_cube(triangulation); 
 *     triangulation.refine_global(6 - dim); 
 * 
 *     if (solver_type == SolverType::FGMRES_ILU) 
 *       std::cout << "Now running with ILU" << std::endl; 
 *     else if (solver_type == SolverType::FGMRES_GMG) 
 *       std::cout << "Now running with Multigrid" << std::endl; 
 *     else 
 *       std::cout << "Now running with UMFPACK" << std::endl; 
 * 
 *     for (unsigned int refinement_cycle = 0; refinement_cycle < 3; 
 *          ++refinement_cycle) 
 *       { 
 *         std::cout << "Refinement cycle " << refinement_cycle << std::endl; 
 * 
 *         if (refinement_cycle > 0) 
 *           triangulation.refine_global(1); 
 * 
 *         std::cout << "   Set-up..." << std::endl; 
 *         setup_dofs(); 
 * 
 *         std::cout << "   Assembling..." << std::endl; 
 *         assemble_system(); 
 * 
 *         if (solver_type == SolverType::FGMRES_GMG) 
 *           { 
 *             std::cout << "   Assembling Multigrid..." << std::endl; 
 * 
 *             assemble_multigrid(); 
 *           } 
 * 
 *         std::cout << "   Solving..." << std::flush; 
 *         solve(); 
 * 
 *         compute_errors(); 
 * 
 *         output_results(refinement_cycle); 
 * 
 *         Utilities::System::MemoryStats mem; 
 *         Utilities::System::get_memory_stats(mem); 
 *         std::cout << "   VM Peak: " << mem.VmPeak << std::endl; 
 * 
 *         computing_timer.print_summary(); 
 *         computing_timer.reset(); 
 *       } 
 *   } 
 * } // namespace Step56 
 * @endcode
 * 
 * 
 * <a name="Themainfunction"></a> 
 * <h3>The main function</h3>
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step56; 
 * 
 *       const int degree = 1; 
 *       const int dim    = 3; 
 * 
 * @endcode
 * 
 * SolverType的选项。umfpack fgmres_ilu fgmres_gmg
 * 

 * 
 * 
 * @code
 *       StokesProblem<dim> flow_problem(degree, SolverType::FGMRES_GMG); 
 * 
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
 * @endcode
examples/step-56/doc/results.dox



<a name="Results"></a><h1>Results</h1>


<a name="Errors"></a><h3> Errors </h3>


我们首先运行代码，确认有限元解以混合有限元问题的误差分析所预测的正确速率收敛。考虑到足够平滑的精确解 $u$ 和 $p$ ，Taylor-Hood元素 $Q_k \times Q_{k-1}$ 的误差应该是

@f[
\| u -u_h \|_0 + h ( \| u- u_h\|_1 + \|p - p_h \|_0)
\leq C h^{k+1} ( \|u \|_{k+1} + \| p \|_k )


@f]



例如见Ern/Guermond《有限元的理论与实践》，第4.2.5节第195页。这确实是我们观察到的，以 $Q_2 \times Q_1$ 元素为例（这就是代码中的做法，但在 <code>main()</code> 中很容易改变）。

 <table align="center" class="doxtable">
  <tr>
    <th>&nbsp;</th>
    <th>L2 Velocity</th>
    <th>Reduction</th>
    <th>L2 Pressure</th>
    <th>Reduction</th>
    <th>H1 Velocity</th>
    <th>Reduction</th>
  </tr>
  <tr>
    <td>3D, 3 global refinements</td>
    <td>0.000670888</td>
    <td align="center">-</td>
    <td>0.0036533</td>
    <td align="center">-</td>
    <td>0.0414704</td>
    <td align="center">-</td>
  </tr>
  <tr>
    <td>3D, 4 global refinements</td>
    <td>8.38E-005</td>
    <td>8.0</td>
    <td>0.00088494</td>
    <td>4.1</td>
    <td>0.0103781</td>
    <td>4.0</td>
  </tr>
  <tr>
    <td>3D, 5 global refinements</td>
    <td>1.05E-005</td>
    <td>8.0</td>
    <td>0.000220253</td>
    <td>4.0</td>
    <td>0.00259519</td>
    <td>4.0</td>
</th>
  </tr>
</table> 

<a name="TimingResults"></a><h3> Timing Results </h3>


让我们比较一下使用UMFPACK的直接求解方法和两种方法，其中我们选择 $\widetilde {A^{-1}}=A^{-1}$ 和 $\widetilde{S^{-1}}=S^{-1}$ ，用CG求解 $A,S$ 的线性系统。然后CG的预处理程序是ILU或GMG。下表总结了求解器的迭代、时序和虚拟内存（VM）的峰值使用。

 <table align="center" class="doxtable">
<tr>
  <th></th>
  <th colspan="3">General</th>
  <th colspan="6">GMG</th>
  <th colspan="6">ILU</th>
  <th colspan="3">UMFPACK</th>
</tr>
<tr>
  <th></th>
  <th></th>
  <th colspan="2">Timings</th>
  <th colspan="2">Timings</th>
  <th colspan="3">Iterations</th>
  <th></th>
  <th colspan="2">Timings</th>
  <th colspan="3">Iterations</th>
  <th></th>
  <th colspan="2">Timings</th>
  <th></th>
</tr>
<tr>
  <th>Cycle</th>
  <th>DoFs</th>
  <th>Setup</th>
  <th>Assembly</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>Outer</th>
  <th>Inner (A)</th>
  <th>Inner (S)</th>
  <th>VM Peak</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>Outer</th>
  <th>Inner (A)</th>
  <th>Inner (S)</th>
  <th>VM Peak</th>
  <th>Setup</th>
  <th>Solve</th>
  <th>VM Peak</th>
</tr>
<tr>
  <td>0</td>
  <td>15468</td>
  <td>0.1s</td>
  <td>0.3s</td>
  <td>0.3s</td>
  <td>1.3s</td>
  <td>21</td>
  <td>67</td>
  <td>22</td>
  <td>4805</td>
  <td>0.3s</td>
  <td>0.6s</td>
  <td>21</td>
  <td>180</td>
  <td>22</td>
  <td>4783</td>
  <td>2.65s</td>
  <td>2.8s</td>
  <td>5054</td>
</tr>
<tr>
  <td>1</td>
  <td>112724</td>
  <td>1.0s</td>
  <td>2.4s</td>
  <td>2.6s</td>
  <td>14s</td>
  <td>21</td>
  <td>67</td>
  <td>22</td>
  <td>5441</td>
  <td>2.8s</td>
  <td>15.8s</td>
  <td>21</td>
  <td>320</td>
  <td>22</td>
  <td>5125</td>
  <td>236s</td>
  <td>237s</td>
  <td>11288</td>
</tr>
<tr>
  <td>2</td>
  <td>859812</td>
  <td>9.0s</td>
  <td>20s</td>
  <td>20s</td>
  <td>101s</td>
  <td>20</td>
  <td>65</td>
  <td>21</td>
  <td>10641</td>
  <td>27s</td>
  <td>268s</td>
  <td>21</td>
  <td>592</td>
  <td>22</td>
  <td>8307</td>
  <td>-</td>
  <td>-</td>
  <td>-</td>
</tr>
</table> 

从表中可以看出。

1.UMFPACK使用了大量的内存，特别是在3D中。另外，UMFPACK的计时并不随问题的大小而变化。

2.因为我们对 $A$ 和 $S$ 使用内部求解器，ILU和GMG需要相同数量的外部迭代。

3.对于ILU来说， $A$ 的（内部）迭代次数随着细化而增加，导致求解时间的线性扩展性较差。相比之下， $A$ 的内部迭代次数在GMG中保持不变，导致求解时间几乎完美的缩放。

4.GMG需要比ILU多一点的内存来存储电平和接口矩阵。

<a name="Possibilitiesforextensions"></a><h3> Possibilities for extensions </h3>


<a name="Checkhigherorderdiscretizations"></a><h4> Check higher order discretizations </h4>


用更高阶的稳定FE对进行实验，检查你是否观察到正确的收敛率。

<a name="Comparewithcheappreconditioner"></a><h4> Compare with cheap preconditioner </h4>


介绍中还概述了对整个系统进行预处理的另一种选择，即我们不选择上表中的 $\widetilde
{A^{-1}}=A^{-1}$ ，而是选择 $\widetilde{A^{-1}}$ ，只分别用GMG或ILU进行单一预处理的应用。

这实际上是在代码中实现的。目前，布尔值 <code>use_expensive</code> in <code>solve()</code> 被设置为 @p true.  上面提到的选项是通过设置为 @p false. 得到的。

你会发现，如果你用GMG这种方式，FGMRES的迭代次数在细化过程中保持不变。这意味着多重网格是最优的，并且与 $h$ 无关。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-56.cc"
*/
