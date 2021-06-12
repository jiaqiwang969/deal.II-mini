/**
@page step_57 The step-57 tutorial program
This tutorial depends on step-15, step-22.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
        <li><a href="#NavierStokesEquations"> Navier Stokes Equations </a>
        <li><a href="#LinearizationofNavierStokesEquations"> Linearization of Navier-Stokes Equations </a>
        <li><a href="#FindinganInitialGuess"> Finding an Initial Guess </a>
        <li><a href="#TheSolverandPreconditioner">The Solver and Preconditioner </a>
        <li><a href="#TestCase"> Test Case </a>
        <li><a href="#References"> References </a>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeNavierStokesProblemcodeclasstemplate">The <code>NavierStokesProblem</code> class template</a>
        <li><a href="#Boundaryvaluesandrighthandside">Boundary values and right hand side</a>
        <li><a href="#BlockSchurPreconditionerforNavierStokesequations">BlockSchurPreconditioner for Navier Stokes equations</a>
        <li><a href="#StationaryNavierStokesclassimplementation">StationaryNavierStokes class implementation</a>
      <ul>
        <li><a href="#StationaryNavierStokesStationaryNavierStokes">StationaryNavierStokes::StationaryNavierStokes</a>
        <li><a href="#StationaryNavierStokessetup_dofs">StationaryNavierStokes::setup_dofs</a>
        <li><a href="#StationaryNavierStokesinitialize_system">StationaryNavierStokes::initialize_system</a>
        <li><a href="#StationaryNavierStokesassemble">StationaryNavierStokes::assemble</a>
        <li><a href="#StationaryNavierStokessolve">StationaryNavierStokes::solve</a>
        <li><a href="#StationaryNavierStokesrefine_mesh">StationaryNavierStokes::refine_mesh</a>
        <li><a href="#StationaryNavierStokesdimnewton_iteration">StationaryNavierStokes<dim>::newton_iteration</a>
        <li><a href="#StationaryNavierStokescompute_initial_guess">StationaryNavierStokes::compute_initial_guess</a>
        <li><a href="#StationaryNavierStokesoutput_results">StationaryNavierStokes::output_results</a>
        <li><a href="#StationaryNavierStokesprocess_solution">StationaryNavierStokes::process_solution</a>
        <li><a href="#StationaryNavierStokesrun">StationaryNavierStokes::run</a>
      </ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
        <li><a href="#Testcase1LowReynoldsNumber"> Test case 1: Low Reynolds Number </a>
        <li><a href="#Testcase2HighReynoldsNumber"> Test case 2: High Reynolds Number </a>
        <li><a href="#Possibilitiesforextensions">Possibilities for extensions</a>
      <ul>
        <li><a href="#Comparetoothersolvers">Compare to other solvers</a>
        <li><a href="#3dcomputations">3d computations</a>
        <li><a href="#Parallelization">Parallelization</a>
    </ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-57/doc/intro.dox

 <br> 

<i>This program was contributed by Liang Zhao and Timo Heister.


This material is based upon work partially supported by National Science
Foundation grant DMS1522191 and the Computational Infrastructure in
Geodynamics initiative (CIG), through the National Science Foundation under
Award No. EAR-0949446 and The University of California-Davis.
</i>

 @dealiiTutorialDOI{10.5281/zenodo.484156,https://zenodo.org/badge/DOI/10.5281/zenodo.484156.svg} 

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


<a name="NavierStokesEquations"></a><h3> Navier Stokes Equations </h3>


在本教程中，我们展示了如何用牛顿方法求解不可压缩的纳维尔-斯托克斯方程（NSE）。我们在这里考虑的流动被假定是稳定的。在一个域 $\Omega \subset
\mathbb{R}^{d}$ ， $d=2,3$ ，具有片状光滑边界 $\partial \Omega$ ，和一个给定的力场 $\textbf{f}$ ，我们寻求一个速度场 $\textbf{u}$ 和压力场 $\textbf{p}$ ，满足以下条件

@f{eqnarray*}


- \nu \Delta\textbf{u} + (\textbf{u} \cdot \nabla)\textbf{u} + \nabla p &=& \textbf{f}\\


- \nabla \cdot \textbf{u} &=& 0.


@f}



与步骤22中讨论的斯托克斯方程不同，由于对流项的存在，NSE是一个非线性方程组  $(\textbf{u} \cdot
\nabla)\textbf{u}$  。计算数值解的第一步是将系统线性化，这将用牛顿方法完成。在步骤35中讨论了一个随时间变化的问题，其中系统是使用最后一个时间步骤的解来线性化的，没有必要进行非线性解。

<a name="LinearizationofNavierStokesEquations"></a><h3> Linearization of Navier-Stokes Equations </h3>


我们定义一个非线性函数，其根是NSE的解，即

@f{eqnarray*}
F(\mathbf{u}, p) =
  \begin{pmatrix}


    - \nu \Delta\mathbf{u} + (\mathbf{u} \cdot \nabla)\mathbf{u} + \nabla p - \mathbf{f} \\


    - \nabla \cdot \mathbf{u}
  \end{pmatrix}.


@f}



假设初始猜测足以保证牛顿迭代的收敛性，并表示 $\textbf{x} = (\textbf{u}, p)$ ，牛顿对向量函数的迭代可以定义为

@f{eqnarray*}
  \textbf{x}^{k+1} = \textbf{x}^{k} - (\nabla F(\textbf{x}^{k}))^{-1} F(\textbf{x}^{k}),


@f}



其中 $\textbf{x}^{k+1}$ 是步骤 $k+1$ 中的近似解， $\textbf{x}^{k}$ 代表上一步的解， $\nabla
F(\textbf{x}^{k})$ 是在 $\textbf{x}^{k}$ 处评估的雅各布矩阵。类似的迭代可以在步骤15中找到。

牛顿迭代公式意味着新的解决方案是通过在旧的解决方案上增加一个更新项而得到的。我们不是评估雅各布矩阵并取其倒数，而是将更新项视为一个整体，即

@f{eqnarray*}
  \delta \textbf{x}^{k} = - (\nabla F(\textbf{x}^{k}))^{-1} F(\textbf{x}^{k}),


@f}



其中  $\textbf{x}^{k+1}=\textbf{x}^{k}+\delta \textbf{x}^{k}$  。

我们可以通过求解系统找到更新项

@f{eqnarray*}
  \nabla F(\textbf{x}^{k}) \delta \textbf{x}^{k} = -F(\textbf{x}^{k}).


@f}



这里，前一个方程的左边代表 $F(\textbf{x})$ 沿 $\delta
\textbf{x}^{k}$ 在 $\textbf{x}^{k}$ 的方向梯度。根据定义，方向性梯度由以下公式给出

@f{eqnarray*}
  & &\nabla F(\mathbf{u}^{k}, p^{k}) (\delta \mathbf{u}^{k}, \delta p^{k}) \\
  \\
  &=& \lim_{\epsilon \to 0} \frac{1}{\epsilon}
      \left(
        F(\mathbf{u}^{k} + \epsilon \delta \mathbf{u}^{k},
          p^{k} + \epsilon \nabla \delta p^{k})


      - F(\mathbf{u}^{k}, p^{k})
      \right)\\
  \\
  &=& \lim_{\epsilon \to 0} \frac{1}{\epsilon}
      \begin{pmatrix}


        - \epsilon \nu \Delta \delta \mathbf{u}^{k}
        + \epsilon \mathbf{u}^{k} \cdot \nabla \delta \mathbf{u}^{k}
        + \epsilon \delta \mathbf{u}^{k} \cdot \nabla \mathbf{u}^{k}
        + \epsilon^{2} \delta \mathbf{u}^{k} \cdot \nabla \delta \mathbf{u}^{k}
        + \epsilon \nabla \delta p^{k}\\


        - \epsilon \nabla \cdot \delta \mathbf{u}^{k}
      \end{pmatrix} \\
  \\
  &=& \begin{pmatrix}


        - \nu \Delta \delta \mathbf{u}^{k}
        + \mathbf{u}^{k} \cdot \nabla \delta \mathbf{u}^{k}
        + \delta \mathbf{u}^{k} \cdot \nabla \mathbf{u}^{k}
        + \nabla \delta p^{k}\\


        - \nabla \cdot \delta \mathbf{u}^{k}
      \end{pmatrix}.


@f}



因此，我们得出了线性化系统。

@f{eqnarray*}


   -\nu \Delta \delta \mathbf{u}^{k}
  + \mathbf{u}^{k} \cdot \nabla \delta \mathbf{u}^{k}
  + \delta \mathbf{u}^{k} \cdot \nabla \mathbf{u}^{k}
  + \nabla \delta p^{k}
  = -F(\mathbf{x}^k), \\


   -\nabla \cdot\delta \mathbf{u}^{k}
  = \nabla \cdot \mathbf{u}^{k},


@f}



其中 $\textbf{u}^k$ 和 $p^k$ 是上一次迭代的解。此外，第二个方程的右边不是零，因为离散解不完全是无发散的（连续解是无发散的）。这里的右手边作为一个修正，导致速度的离散解沿着牛顿的迭代是无发散的。在这个线性系统中，唯一的未知数是更新项 $\delta \textbf{u}^{k}$ 和 $\delta p^{k}$ ，我们可以使用类似于步骤22中的策略（并以同样的方式推导出弱形式）。

现在，可以用牛顿迭代法来解决更新项。

<ol>  <li>  初始化。初始猜测  $u_0$  和  $p_0$  ，公差  $\tau$  ;  </li>   <li>  线性求解计算更新项  $\delta\textbf{u}^{k}$  和  $\delta p^k$  ;  </li>   <li>  更新近似。         $\textbf{u}^{k+1} = \textbf{u}^{k} + \delta\textbf{u}^{k}$  和  $p^{k+1} = p^{k} + \delta p^{k}$  ;  </li>   <li>  检查残差规范。   $E^{k+1} = \|F(\mathbf{u}^{k+1}, p^{k+1})\|$  :  <ul>   <li>  如果  $E^{k+1} \leq \tau$  , 停止。 </li>   <li>  如果  $E^{k+1} > \tau$  ，回到步骤2。 </li>   </ul>   </li>   </ol> 。

<a name="FindinganInitialGuess"></a><h3> Finding an Initial Guess </h3>


初始猜测需要足够接近解决方案，牛顿方法才能收敛；因此，找到一个好的起始值对非线性求解器来说至关重要。

当粘度 $\nu$ 较大时，通过解决带有粘度 $\nu$ 的斯托克斯方程可以得到一个好的初始猜测。虽然这取决于问题，但对于这里考虑的测试问题，这对 $\nu \geq 1/400$ 有效。

然而，如果粘度较小，对流项 $(\mathbf{u}\cdot\nabla)\mathbf{u}$ 将占主导地位，如测试案例2的 $1/7500$ 。  在这种情况下，我们使用延续法建立一系列的辅助NSE，其粘度接近于目标NSE中的粘度。相应地，我们创建了一个带有 $\{\nu_{i}\}$ 的序列，如果 $\nu_{n}= \nu$ 较小，我们接受两个具有粘度的NSE $\nu_{i}$ 和 $\nu_{i+1}$ 的解是接近的。  然后，我们使用带有粘度的NSE的解 $\nu_{i}$ 作为带有 $\nu_{i+1}$ 的NSE的初始猜测。这可以被认为是一个从斯托克斯方程到我们要解决的NSE的阶梯。

也就是说，我们首先解决一个斯托克斯问题

@f{eqnarray*}


  -\nu_{1} \Delta \textbf{u} + \nabla p &=& \textbf{f}\\


  -\nabla \cdot \textbf{u} &=& 0


@f}



以获得最初的猜测，即

@f{eqnarray*}


  -\nu_{1} \Delta \textbf{u} + (\textbf{u} \cdot \nabla)\textbf{u} + \nabla p &=& \textbf{f},\\


  -\nabla \cdot \textbf{u} &=& 0,


@f}



这也是延续法的初始猜测。这里 $\nu_{1}$ 相对较大，这样带粘度的斯托克斯问题的解 $\nu_{1}$ 就可以作为牛顿迭代中NSE的初始猜测。

那么，对

@f{eqnarray*}


  -\nu_{i} \Delta \textbf{u} + (\textbf{u} \cdot \nabla)\textbf{u} + \nabla p &=& \textbf{f},\\


  -\nabla \cdot \textbf{u} &=& 0.


@f}



的初始猜测。

@f{eqnarray*}


  -\nu_{i+1} \Delta \textbf{u} + (\textbf{u} \cdot \nabla)\textbf{u} + \nabla p &=& \textbf{f},\\


  -\nabla \cdot \textbf{u} &=& 0.


@f}



这个过程是用实验确定的粘度序列 $\{\nu_i\}$ 来重复的，这样最终的解决方案可以作为牛顿迭代的起始猜测。

<a name="TheSolverandPreconditioner"></a><h3>The %Solver and Preconditioner </h3>


在牛顿迭代的每一步，问题的结果是解决一个鞍点系统，其形式为

@f{eqnarray*}
    \begin{pmatrix}
      A & B^{T} \\
      B & 0
    \end{pmatrix}
    \begin{pmatrix}
      U \\
      P
    \end{pmatrix}
    =
    \begin{pmatrix}
      F \\
      0
    \end{pmatrix}.


@f}



这个系统矩阵的块状结构与步骤22中的相同。然而，左上角的矩阵 $A$ 不是对称的，因为有非线性项。我们可以不求解上述系统，而求解等效系统

@f{eqnarray*}
    \begin{pmatrix}
      A + \gamma B^TW^{-1}B & B^{T} \\
      B & 0
    \end{pmatrix}
    \begin{pmatrix}
      U \\
      P
    \end{pmatrix}
    =
    \begin{pmatrix}
      F \\
      0
    \end{pmatrix}


@f}



有一个参数  $\gamma$  和一个可逆矩阵  $W$  。这里 $\gamma B^TW^{-1}B$ 是Augmented Lagrangian项，详见[1]。

用 $G$ 表示新系统的系统矩阵，用 $b$ 表示右手边，我们用右预处理 $P^{-1}$ 反复求解，即 $GP^{-1}y = b$  ，其中

@f{eqnarray*}
P^{-1} =
  \begin{pmatrix}
    \tilde{A} & B^T \\
    0         & \tilde{S}
  \end{pmatrix}^{-1}


@f}



与 $\tilde{A} = A + \gamma B^TW^{-1}B$ ， $\tilde{S}$ 是相应的舒尔补数 $\tilde{S} = B^T \tilde{A}^{-1} B$  。我们让 $W = M_p$ 其中 $M_p$ 是压力质量矩阵，那么 $\tilde{S}^{-1}$ 可以近似为

@f{eqnarray*}
\tilde{S}^{-1} \approx -(\nu+\gamma)M_p^{-1}.


@f}



详见[1]。

我们将 $P^{-1}$ 分解为

@f{eqnarray*}
P^{-1} =
  \begin{pmatrix}
    \tilde{A}^{-1} & 0 \\
    0              & I
  \end{pmatrix}
  \begin{pmatrix}
    I & -B^T \\
    0 & I
  \end{pmatrix}
  \begin{pmatrix}
    I & 0 \\
    0 & \tilde{S}^{-1}
  \end{pmatrix}.


@f}



这里需要两个非精确求解器分别用于 $\tilde{A}^{-1}$ 和 $\tilde{S}^{-1}$ （见[1]）。由于压力质量矩阵是对称的和正定的，用ILU作为预处理程序的CG适合用于 $\tilde{S}^{-1}$ 。为了简单起见，我们对  $\tilde{A}^{-1}$  使用直接求解器UMFPACK。最后一个成分是与  $B^T$  的稀疏矩阵-向量乘积。我们不计算 $\tilde{A}$ 中增强拉格朗日项的矩阵乘积，而是集合Grad-Div稳定化 $(\nabla \cdot \phi _{i}, \nabla \cdot \phi _{j}) \approx (B^T
M_p^{-1}B)_{ij}$ ，如[2]中所解释。

<a name="TestCase"></a><h3> Test Case </h3>


我们使用盖子驱动的空腔流作为我们的测试案例，详情见[3]。计算域为单位平方，右手边为 $f=0$  。边界条件为

@f{eqnarray*}
  (u(x, y), v(x,y)) &=& (1,0) \qquad\qquad \textrm{if}\ y = 1 \\
  (u(x, y), v(x,y)) &=& (0,0) \qquad\qquad \textrm{otherwise}.


@f}



在解决这个问题时，误差由非线性误差（来自牛顿迭代）和离散化误差（取决于网格大小）组成。非线性部分随着牛顿的每次迭代而减少，离散化误差随着网格的细化而减少。在这个例子中，粗网的解被转移到连续的细网中，并作为初始猜测使用。因此，非线性误差总是低于牛顿迭代的容忍度，离散化误差也随着每次网格细化而减少。

在循环中，我们涉及三个求解器：一个用于 $\tilde{A}^{-1}$ ，一个用于 $M_p^{-1}$ 和一个用于 $Gx=b$  。前两个求解器在预处理程序中被调用，外部求解器给我们提供了更新项。总体收敛性由非线性残差控制；由于牛顿方法不需要精确的雅各布，我们采用FGMRES，外侧线性求解器的相对公差仅为1e-4。事实上，我们对这个系统采用了截断的牛顿解法。如步骤22所述，内部线性求解也不需要做得非常精确。这里我们使用CG，压力质量矩阵的相对公差为1e-6。正如预期的那样，我们仍然看到非线性残差收敛到了1e-14。另外，我们使用一个简单的线搜索算法来实现牛顿方法的全球化。

腔体参考值 $\mathrm{Re}=400$ 和 $\mathrm{Re}=7500$ 分别来自[4]和[5]，其中 $\mathrm{Re}$ 是雷诺数，可以定位在[8]。这里的粘度是由 $1/\mathrm{Re}$ 定义的。尽管我们仍然可以找到 $\mathrm{Re}=10000$ 的解决方案，而且参考文献中也包含了用于比较的结果，但我们在这里的讨论仅限于 $\mathrm{Re}=7500$  。这是因为从 $\mathrm{Re}=8000$ 附近开始，解不再是静止的，而是变成了周期性的，详见[7]。

<a name="References"></a><h3> References </h3> <醇>


    <li>  An Augmented Lagrangian-Based Approach to the Oseen Problem, M. Benzi and M. Olshanskii, SIAM J. SCI. COMPUT.COMPUT.2006  <li>  Efficient augmented Lagrangian-type preconditioning for the Oseen problem using Grad-Div stabilization, Timo Heister and Gerd Rapin  <li>  http://www.cfd-online.com/Wiki/Lid-driven_cavity_problem  <li>  High-Re solution for incompressible flow using the Navier-Stokes equations and a Multigrid Method, U. Ghia, K. N. Ghia, and C. T. Shin  <li>  E. Erturk, T.C. Corke and C. Gokcol  <li>  三维不可压缩Navier-Stokes方程的隐式加权ENO方案，Yang等人，1998  <li>  二维盖子驱动的空腔问题再探讨，C. Bruneau和M. Saad，2006  <li>  https://en.wikipedia.org/wiki/Reynolds_number  </ol> 


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
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/utilities.h> 
 * #include <deal.II/base/tensor.h> 
 * 
 * #include <deal.II/lac/block_vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/block_sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/solver_gmres.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/lac/affine_constraints.h> 
 * 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/grid_refinement.h> 
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
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * #include <deal.II/numerics/error_estimator.h> 
 * 
 * @endcode
 * 
 * 为了在网格之间传输解决方案，包括这个文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/numerics/solution_transfer.h> 
 * 
 * @endcode
 * 
 * 这个文件包括UMFPACK：直接求解器。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_direct.h> 
 * 
 * @endcode
 * 
 * 还有一个ILU预处理程序。
 * 

 * 
 * 
 * @code
 * #include <deal.II/lac/sparse_ilu.h> 
 * 
 * #include <fstream> 
 * #include <iostream> 
 * 
 * namespace Step57 
 * { 
 *   using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeNavierStokesProblemcodeclasstemplate"></a> 
 * <h3>The <code>NavierStokesProblem</code> class template</h3>
 * 

 * 
 * 该类管理介绍中描述的矩阵和向量：特别是，我们为当前的解决方案、当前的牛顿更新和直线搜索更新存储了一个BlockVector。 我们还存储了两个AffineConstraints对象：一个是强制执行Dirichlet边界条件的对象，另一个是将所有边界值设为0的对象。第一个约束解向量，第二个约束更新（也就是说，我们从不更新边界值，所以我们强制相关的更新向量值为零）。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   class StationaryNavierStokes 
 *   { 
 *   public: 
 *     StationaryNavierStokes(const unsigned int degree); 
 *     void run(const unsigned int refinement); 
 * 
 *   private: 
 *     void setup_dofs(); 
 * 
 *     void initialize_system(); 
 * 
 *     void assemble(const bool initial_step, const bool assemble_matrix); 
 * 
 *     void assemble_system(const bool initial_step); 
 * 
 *     void assemble_rhs(const bool initial_step); 
 * 
 *     void solve(const bool initial_step); 
 * 
 *     void refine_mesh(); 
 * 
 *     void process_solution(unsigned int refinement); 
 * 
 *     void output_results(const unsigned int refinement_cycle) const; 
 * 
 *     void newton_iteration(const double       tolerance, 
 *                           const unsigned int max_n_line_searches, 
 *                           const unsigned int max_n_refinements, 
 *                           const bool         is_initial_step, 
 *                           const bool         output_result); 
 * 
 *     void compute_initial_guess(double step_size); 
 * 
 *     double                               viscosity; 
 *     double                               gamma; 
 *     const unsigned int                   degree; 
 *     std::vector<types::global_dof_index> dofs_per_block; 
 * 
 *     Triangulation<dim> triangulation; 
 *     FESystem<dim>      fe; 
 *     DoFHandler<dim>    dof_handler; 
 * 
 *     AffineConstraints<double> zero_constraints; 
 *     AffineConstraints<double> nonzero_constraints; 
 * 
 *     BlockSparsityPattern      sparsity_pattern; 
 *     BlockSparseMatrix<double> system_matrix; 
 *     SparseMatrix<double>      pressure_mass_matrix; 
 * 
 *     BlockVector<double> present_solution; 
 *     BlockVector<double> newton_update; 
 *     BlockVector<double> system_rhs; 
 *     BlockVector<double> evaluation_point; 
 *   }; 
 * @endcode
 * 
 * 
 * <a name="Boundaryvaluesandrighthandside"></a> 
 * <h3>Boundary values and right hand side</h3>
 * 

 * 
 * 在这个问题中，我们设定沿空腔上表面的速度为1，其他三面墙的速度为0。右边的函数为零，所以我们在本教程中不需要设置右边的函数。边界函数的分量数为  <code>dim+1</code>  。我们最终将使用 VectorTools::interpolate_boundary_values 来设置边界值，这就要求边界值函数的分量数与解相同，即使没有全部使用。换个说法：为了让这个函数高兴，我们为压力定义了边界值，尽管我们实际上永远不会用到它们。
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
 *     virtual double value(const Point<dim> & p, 
 *                          const unsigned int component) const override; 
 *   }; 
 * 
 *   template <int dim> 
 *   double BoundaryValues<dim>::value(const Point<dim> & p, 
 *                                     const unsigned int component) const 
 *   { 
 *     Assert(component < this->n_components, 
 *            ExcIndexRange(component, 0, this->n_components)); 
 *     if (component == 0 && std::abs(p[dim - 1] - 1.0) < 1e-10) 
 *       return 1.0; 
 * 
 *     return 0; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="BlockSchurPreconditionerforNavierStokesequations"></a> 
 * <h3>BlockSchurPreconditioner for Navier Stokes equations</h3>
 * 

 * 
 * 正如介绍中所讨论的，Krylov迭代方法中的预处理器是作为一个矩阵-向量乘积算子实现的。在实践中，舒尔补码预处理器被分解为三个矩阵的乘积（如第一节所述）。第一个因素中的 $\tilde{A}^{-1}$ 涉及到对线性系统 $\tilde{A}x=b$ 的求解。在这里，为了简单起见，我们通过一个直接求解器来解决这个系统。第二个因素中涉及的计算是一个简单的矩阵-向量乘法。舒尔补码 $\tilde{S}$ 可以被压力质量矩阵很好地近似，其逆值可以通过不精确求解器得到。因为压力质量矩阵是对称和正定的，我们可以用CG来解决相应的线性系统。
 * 

 * 
 * 
 * @code
 *   template <class PreconditionerMp> 
 *   class BlockSchurPreconditioner : public Subscriptor 
 *   { 
 *   public: 
 *     BlockSchurPreconditioner(double                           gamma, 
 *                              double                           viscosity, 
 *                              const BlockSparseMatrix<double> &S, 
 *                              const SparseMatrix<double> &     P, 
 *                              const PreconditionerMp &         Mppreconditioner); 
 * 
 *     void vmult(BlockVector<double> &dst, const BlockVector<double> &src) const; 
 * 
 *   private: 
 *     const double                     gamma; 
 *     const double                     viscosity; 
 *     const BlockSparseMatrix<double> &stokes_matrix; 
 *     const SparseMatrix<double> &     pressure_mass_matrix; 
 *     const PreconditionerMp &         mp_preconditioner; 
 *     SparseDirectUMFPACK              A_inverse; 
 *   }; 
 * 
 * @endcode
 * 
 * 我们可以注意到，左上角的矩阵逆的初始化是在构造函数中完成的。如果是这样，那么预处理程序的每一次应用就不再需要计算矩阵因子了。
 * 

 * 
 * 
 * @code
 *   template <class PreconditionerMp> 
 *   BlockSchurPreconditioner<PreconditionerMp>::BlockSchurPreconditioner( 
 *     double                           gamma, 
 *     double                           viscosity, 
 *     const BlockSparseMatrix<double> &S, 
 *     const SparseMatrix<double> &     P, 
 *     const PreconditionerMp &         Mppreconditioner) 
 *     : gamma(gamma) 
 *     , viscosity(viscosity) 
 *     , stokes_matrix(S) 
 *     , pressure_mass_matrix(P) 
 *     , mp_preconditioner(Mppreconditioner) 
 *   { 
 *     A_inverse.initialize(stokes_matrix.block(0, 0)); 
 *   } 
 * 
 *   template <class PreconditionerMp> 
 *   void BlockSchurPreconditioner<PreconditionerMp>::vmult( 
 *     BlockVector<double> &      dst, 
 *     const BlockVector<double> &src) const 
 *   { 
 *     Vector<double> utmp(src.block(0)); 
 * 
 *     { 
 *       SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm()); 
 *       SolverCG<Vector<double>> cg(solver_control); 
 * 
 *       dst.block(1) = 0.0; 
 *       cg.solve(pressure_mass_matrix, 
 *                dst.block(1), 
 *                src.block(1), 
 *                mp_preconditioner); 
 *       dst.block(1) *= -(viscosity + gamma); 
 *     } 
 * 
 *     { 
 *       stokes_matrix.block(0, 1).vmult(utmp, dst.block(1)); 
 *       utmp *= -1.0; 
 *       utmp += src.block(0); 
 *     } 
 * 
 *     A_inverse.vmult(dst.block(0), utmp); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesclassimplementation"></a> 
 * <h3>StationaryNavierStokes class implementation</h3>
 * 
 * <a name="StationaryNavierStokesStationaryNavierStokes"></a> 
 * <h4>StationaryNavierStokes::StationaryNavierStokes</h4>
 * 

 * 
 * 该类的构造函数看起来与  step-22  中的构造函数非常相似。唯一的区别是粘度和增强的拉格朗日系数  <code>gamma</code>  。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   StationaryNavierStokes<dim>::StationaryNavierStokes(const unsigned int degree) 
 *     : viscosity(1.0 / 7500.0) 
 *     , gamma(1.0) 
 *     , degree(degree) 
 *     , triangulation(Triangulation<dim>::maximum_smoothing) 
 *     , fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1) 
 *     , dof_handler(triangulation) 
 *   {} 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokessetup_dofs"></a> 
 * <h4>StationaryNavierStokes::setup_dofs</h4>
 * 

 * 
 * 这个函数初始化DoFHandler，列举当前网格上的自由度和约束。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::setup_dofs() 
 *   { 
 *     system_matrix.clear(); 
 *     pressure_mass_matrix.clear(); 
 * 
 * @endcode
 * 
 * 第一步是将DoFs与给定的网格联系起来。
 * 

 * 
 * 
 * @code
 *     dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 我们对组件重新编号，使所有的速度DoF在压力DoF之前，以便能够将解向量分成两个块，在块预处理程序中分别访问。
 * 

 * 
 * 
 * @code
 *     std::vector<unsigned int> block_component(dim + 1, 0); 
 *     block_component[dim] = 1; 
 *     DoFRenumbering::component_wise(dof_handler, block_component); 
 * 
 *     dofs_per_block = 
 *       DoFTools::count_dofs_per_fe_block(dof_handler, block_component); 
 *     unsigned int dof_u = dofs_per_block[0]; 
 *     unsigned int dof_p = dofs_per_block[1]; 
 * 
 * @endcode
 * 
 * 在牛顿方案中，我们首先将边界条件应用于从初始步骤得到的解。为了确保边界条件在牛顿迭代过程中保持满足，在更新时使用零边界条件  $\delta u^k$  。因此我们设置了两个不同的约束对象。
 * 

 * 
 * 
 * @code
 *     FEValuesExtractors::Vector velocities(0); 
 *     { 
 *       nonzero_constraints.clear(); 
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, nonzero_constraints); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                BoundaryValues<dim>(), 
 *                                                nonzero_constraints, 
 *                                                fe.component_mask(velocities)); 
 *     } 
 *     nonzero_constraints.close(); 
 * 
 *     { 
 *       zero_constraints.clear(); 
 * 
 *       DoFTools::make_hanging_node_constraints(dof_handler, zero_constraints); 
 *       VectorTools::interpolate_boundary_values(dof_handler, 
 *                                                0, 
 *                                                Functions::ZeroFunction<dim>( 
 *                                                  dim + 1), 
 *                                                zero_constraints, 
 *                                                fe.component_mask(velocities)); 
 *     } 
 *     zero_constraints.close(); 
 * 
 *     std::cout << "Number of active cells: " << triangulation.n_active_cells() 
 *               << std::endl 
 *               << "Number of degrees of freedom: " << dof_handler.n_dofs() 
 *               << " (" << dof_u << " + " << dof_p << ')' << std::endl; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesinitialize_system"></a> 
 * <h4>StationaryNavierStokes::initialize_system</h4>
 * 

 * 
 * 在每个网格上，SparsityPattern和线性系统的大小是不同的。这个函数在网格细化后初始化它们。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::initialize_system() 
 *   { 
 *     { 
 *       BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block); 
 *       DoFTools::make_sparsity_pattern(dof_handler, dsp, nonzero_constraints); 
 *       sparsity_pattern.copy_from(dsp); 
 *     } 
 * 
 *     system_matrix.reinit(sparsity_pattern); 
 * 
 *     present_solution.reinit(dofs_per_block); 
 *     newton_update.reinit(dofs_per_block); 
 *     system_rhs.reinit(dofs_per_block); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesassemble"></a> 
 * <h4>StationaryNavierStokes::assemble</h4>
 * 

 * 
 * 这个函数建立了我们目前工作的系统矩阵和右手边。 @p initial_step 参数用于确定我们应用哪一组约束（初始步骤为非零，其他为零）。 @p assemble_matrix 参数分别决定了是组装整个系统还是只组装右手边的向量。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::assemble(const bool initial_step, 
 *                                              const bool assemble_matrix) 
 *   { 
 *     if (assemble_matrix) 
 *       system_matrix = 0; 
 * 
 *     system_rhs = 0; 
 * 
 *     QGauss<dim> quadrature_formula(degree + 2); 
 * 
 *     FEValues<dim> fe_values(fe, 
 *                             quadrature_formula, 
 *                             update_values | update_quadrature_points | 
 *                               update_JxW_values | update_gradients); 
 * 
 *     const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 *     const unsigned int n_q_points    = quadrature_formula.size(); 
 * 
 *     const FEValuesExtractors::Vector velocities(0); 
 *     const FEValuesExtractors::Scalar pressure(dim); 
 * 
 *     FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell); 
 *     Vector<double>     local_rhs(dofs_per_cell); 
 * 
 *     std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 对于线性化系统，我们为当前速度和梯度以及当前压力创建临时存储。在实践中，它们都是通过正交点的形状函数获得的。
 * 

 * 
 * 
 * @code
 *     std::vector<Tensor<1, dim>> present_velocity_values(n_q_points); 
 *     std::vector<Tensor<2, dim>> present_velocity_gradients(n_q_points); 
 *     std::vector<double>         present_pressure_values(n_q_points); 
 * 
 *     std::vector<double>         div_phi_u(dofs_per_cell); 
 *     std::vector<Tensor<1, dim>> phi_u(dofs_per_cell); 
 *     std::vector<Tensor<2, dim>> grad_phi_u(dofs_per_cell); 
 *     std::vector<double>         phi_p(dofs_per_cell); 
 * 
 *     for (const auto &cell : dof_handler.active_cell_iterators()) 
 *       { 
 *         fe_values.reinit(cell); 
 * 
 *         local_matrix = 0; 
 *         local_rhs    = 0; 
 * 
 *         fe_values[velocities].get_function_values(evaluation_point, 
 *                                                   present_velocity_values); 
 * 
 *         fe_values[velocities].get_function_gradients( 
 *           evaluation_point, present_velocity_gradients); 
 * 
 *         fe_values[pressure].get_function_values(evaluation_point, 
 *                                                 present_pressure_values); 
 * 
 * @endcode
 * 
 * 装配类似于  step-22  。一个以gamma为系数的附加项是增强拉格朗日（AL），它是通过grad-div稳定化组装的。 正如我们在介绍中所讨论的，系统矩阵的右下块应该为零。由于压力质量矩阵是在创建预处理程序时使用的，所以我们在这里组装它，然后在最后把它移到一个单独的SparseMatrix中（与 step-22 相同）。
 * 

 * 
 * 
 * @code
 *         for (unsigned int q = 0; q < n_q_points; ++q) 
 *           { 
 *             for (unsigned int k = 0; k < dofs_per_cell; ++k) 
 *               { 
 *                 div_phi_u[k]  = fe_values[velocities].divergence(k, q); 
 *                 grad_phi_u[k] = fe_values[velocities].gradient(k, q); 
 *                 phi_u[k]      = fe_values[velocities].value(k, q); 
 *                 phi_p[k]      = fe_values[pressure].value(k, q); 
 *               } 
 * 
 *             for (unsigned int i = 0; i < dofs_per_cell; ++i) 
 *               { 
 *                 if (assemble_matrix) 
 *                   { 
 *                     for (unsigned int j = 0; j < dofs_per_cell; ++j) 
 *                       { 
 *                         local_matrix(i, j) += 
 *                           (viscosity * 
 *                              scalar_product(grad_phi_u[j], grad_phi_u[i]) + 
 *                            present_velocity_gradients[q] * phi_u[j] * phi_u[i] + 
 *                            grad_phi_u[j] * present_velocity_values[q] * 
 *                              phi_u[i] - 
 *                            div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j] + 
 *                            gamma * div_phi_u[j] * div_phi_u[i] + 
 *                            phi_p[i] * phi_p[j]) * 
 *                           fe_values.JxW(q); 
 *                       } 
 *                   } 
 * 
 *                 double present_velocity_divergence = 
 *                   trace(present_velocity_gradients[q]); 
 *                 local_rhs(i) += 
 *                   (-viscosity * scalar_product(present_velocity_gradients[q], 
 *                                                grad_phi_u[i]) - 
 *                    present_velocity_gradients[q] * present_velocity_values[q] * 
 *                      phi_u[i] + 
 *                    present_pressure_values[q] * div_phi_u[i] + 
 *                    present_velocity_divergence * phi_p[i] - 
 *                    gamma * present_velocity_divergence * div_phi_u[i]) * 
 *                   fe_values.JxW(q); 
 *               } 
 *           } 
 * 
 *         cell->get_dof_indices(local_dof_indices); 
 * 
 *         const AffineConstraints<double> &constraints_used = 
 *           initial_step ? nonzero_constraints : zero_constraints; 
 * 
 *         if (assemble_matrix) 
 *           { 
 *             constraints_used.distribute_local_to_global(local_matrix, 
 *                                                         local_rhs, 
 *                                                         local_dof_indices, 
 *                                                         system_matrix, 
 *                                                         system_rhs); 
 *           } 
 *         else 
 *           { 
 *             constraints_used.distribute_local_to_global(local_rhs, 
 *                                                         local_dof_indices, 
 *                                                         system_rhs); 
 *           } 
 *       } 
 * 
 *     if (assemble_matrix) 
 *       { 
 * 
 * @endcode
 * 
 * 最后我们把压力质量矩阵移到一个单独的矩阵中。
 * 

 * 
 * 
 * @code
 *         pressure_mass_matrix.reinit(sparsity_pattern.block(1, 1)); 
 *         pressure_mass_matrix.copy_from(system_matrix.block(1, 1)); 
 * 
 * @endcode
 * 
 * 注意，将这个压力块设置为零并不等同于不在这个块中装配任何东西，因为这里的操作将（错误地）删除从压力作用力的悬挂节点约束中进来的对角线条目。这意味着，我们的整个系统矩阵将有完全为零的行。幸运的是，FGMRES处理这些行没有任何问题。
 * 

 * 
 * 
 * @code
 *         system_matrix.block(1, 1) = 0; 
 *       } 
 *   } 
 * 
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::assemble_system(const bool initial_step) 
 *   { 
 *     assemble(initial_step, true); 
 *   } 
 * 
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::assemble_rhs(const bool initial_step) 
 *   { 
 *     assemble(initial_step, false); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokessolve"></a> 
 * <h4>StationaryNavierStokes::solve</h4>
 * 

 * 
 * 在这个函数中，我们使用FGMRES和程序开始时定义的块状预处理程序来解决线性系统。我们在这一步得到的是解向量。如果这是初始步骤，解向量为我们提供了纳维尔-斯托克斯方程的初始猜测。对于初始步骤，非零约束被应用，以确保边界条件得到满足。在下面的步骤中，我们将求解牛顿更新，所以使用零约束。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::solve(const bool initial_step) 
 *   { 
 *     const AffineConstraints<double> &constraints_used = 
 *       initial_step ? nonzero_constraints : zero_constraints; 
 * 
 *     SolverControl solver_control(system_matrix.m(), 
 *                                  1e-4 * system_rhs.l2_norm(), 
 *                                  true); 
 * 
 *     SolverFGMRES<BlockVector<double>> gmres(solver_control); 
 *     SparseILU<double>                 pmass_preconditioner; 
 *     pmass_preconditioner.initialize(pressure_mass_matrix, 
 *                                     SparseILU<double>::AdditionalData()); 
 * 
 *     const BlockSchurPreconditioner<SparseILU<double>> preconditioner( 
 *       gamma, 
 *       viscosity, 
 *       system_matrix, 
 *       pressure_mass_matrix, 
 *       pmass_preconditioner); 
 * 
 *     gmres.solve(system_matrix, newton_update, system_rhs, preconditioner); 
 *     std::cout << "FGMRES steps: " << solver_control.last_step() << std::endl; 
 * 
 *     constraints_used.distribute(newton_update); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesrefine_mesh"></a> 
 * <h4>StationaryNavierStokes::refine_mesh</h4>
 * 

 * 
 * 在粗略的网格上找到一个好的初始猜测后，我们希望通过细化网格来减少误差。这里我们做了类似于 step-15 的自适应细化，只是我们只使用了速度上的Kelly估计器。我们还需要使用SolutionTransfer类将当前的解转移到下一个网格。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::refine_mesh() 
 *   { 
 *     Vector<float> estimated_error_per_cell(triangulation.n_active_cells()); 
 *     FEValuesExtractors::Vector velocity(0); 
 *     KellyErrorEstimator<dim>::estimate( 
 *       dof_handler, 
 *       QGauss<dim - 1>(degree + 1), 
 *       std::map<types::boundary_id, const Function<dim> *>(), 
 *       present_solution, 
 *       estimated_error_per_cell, 
 *       fe.component_mask(velocity)); 
 * 
 *     GridRefinement::refine_and_coarsen_fixed_number(triangulation, 
 *                                                     estimated_error_per_cell, 
 *                                                     0.3, 
 *                                                     0.0); 
 * 
 *     triangulation.prepare_coarsening_and_refinement(); 
 *     SolutionTransfer<dim, BlockVector<double>> solution_transfer(dof_handler); 
 *     solution_transfer.prepare_for_coarsening_and_refinement(present_solution); 
 *     triangulation.execute_coarsening_and_refinement(); 
 * 
 * @endcode
 * 
 * 首先，DoFHandler被设置，约束被生成。然后我们创建一个临时的BlockVector  <code>tmp</code>  ，其大小与新网格上的解决方案一致。
 * 

 * 
 * 
 * @code
 *     setup_dofs(); 
 * 
 *     BlockVector<double> tmp(dofs_per_block); 
 * 
 * @endcode
 * 
 * 将解决方案从粗网格转移到细网格，并对新转移的解决方案应用边界值约束。注意，present_solution仍然是对应于旧网格的一个向量。
 * 

 * 
 * 
 * @code
 *     solution_transfer.interpolate(present_solution, tmp); 
 *     nonzero_constraints.distribute(tmp); 
 * 
 * @endcode
 * 
 * 最后设置矩阵和向量，并将present_solution设置为插值后的数据。
 * 

 * 
 * 
 * @code
 *     initialize_system(); 
 *     present_solution = tmp; 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesdimnewton_iteration"></a> 
 * <h4>StationaryNavierStokes<dim>::newton_iteration</h4>
 * 

 * 
 * 这个函数实现了牛顿迭代，给定了公差、最大迭代次数和要做的网格细化次数。
 * 

 * 
 * 参数 <code>is_initial_step</code> 告诉我们是否需要 <code>setup_system</code> ，以及应该装配哪一部分，系统矩阵或右手边的矢量。如果我们做直线搜索，在最后一次迭代中检查残差准则时，右手边已经被组装起来了。因此，我们只需要在当前迭代中装配系统矩阵。最后一个参数 <code>output_result</code> 决定了是否应该产生图形输出。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::newton_iteration( 
 *     const double       tolerance, 
 *     const unsigned int max_n_line_searches, 
 *     const unsigned int max_n_refinements, 
 *     const bool         is_initial_step, 
 *     const bool         output_result) 
 *   { 
 *     bool first_step = is_initial_step; 
 * 
 *     for (unsigned int refinement_n = 0; refinement_n < max_n_refinements + 1; 
 *          ++refinement_n) 
 *       { 
 *         unsigned int line_search_n = 0; 
 *         double       last_res      = 1.0; 
 *         double       current_res   = 1.0; 
 *         std::cout << "grid refinements: " << refinement_n << std::endl 
 *                   << "viscosity: " << viscosity << std::endl; 
 * 
 *         while ((first_step || (current_res > tolerance)) && 
 *                line_search_n < max_n_line_searches) 
 *           { 
 *             if (first_step) 
 *               { 
 *                 setup_dofs(); 
 *                 initialize_system(); 
 *                 evaluation_point = present_solution; 
 *                 assemble_system(first_step); 
 *                 solve(first_step); 
 *                 present_solution = newton_update; 
 *                 nonzero_constraints.distribute(present_solution); 
 *                 first_step       = false; 
 *                 evaluation_point = present_solution; 
 *                 assemble_rhs(first_step); 
 *                 current_res = system_rhs.l2_norm(); 
 *  
 *  
 *  
 *  
 *  
 *  
 *                 evaluation_point = present_solution; 
 *                 assemble_system(first_step); 
 *                 solve(first_step); 
 * 
 * @endcode
 * 
 * 为了确保我们的解决方案越来越接近精确的解决方案，我们让解决方案用权重 <code>alpha</code> 更新，使新的残差小于上一步的残差，这是在下面的循环中完成。这与  step-15  中使用的线搜索算法相同。
 * 

 * 
 * 
 * @code
 *                 for (double alpha = 1.0; alpha > 1e-5; alpha *= 0.5) 
 *                   { 
 *                     evaluation_point = present_solution; 
 *                     evaluation_point.add(alpha, newton_update); 
 *                     nonzero_constraints.distribute(evaluation_point); 
 *                     assemble_rhs(first_step); 
 *                     current_res = system_rhs.l2_norm(); 
 *                     std::cout << "  alpha: " << std::setw(10) << alpha 
 *                               << std::setw(0) << "  residual: " << current_res 
 *                               << std::endl; 
 *                     if (current_res < last_res) 
 *                       break; 
 *                   } 
 *                 { 
 *                   present_solution = evaluation_point; 
 *                   std::cout << "  number of line searches: " << line_search_n 
 *                             << "  residual: " << current_res << std::endl; 
 *                   last_res = current_res; 
 *                 } 
 *                 ++line_search_n; 
 *               } 
 * 
 *             if (output_result) 
 *               { 
 *                 output_results(max_n_line_searches * refinement_n + 
 *                                line_search_n); 
 * 
 *                 if (current_res <= tolerance) 
 *                   process_solution(refinement_n); 
 *               } 
 *           } 
 * 
 *         if (refinement_n < max_n_refinements) 
 *           { 
 *             refine_mesh(); 
 *           } 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokescompute_initial_guess"></a> 
 * <h4>StationaryNavierStokes::compute_initial_guess</h4>
 * 

 * 
 * 这个函数将通过使用延续法为我们提供一个初始猜测，正如我们在介绍中讨论的那样。雷诺数被逐级增加 step- ，直到我们达到目标值。通过实验，斯托克斯的解足以成为雷诺数为1000的NSE的初始猜测，所以我们从这里开始。 为了确保前一个问题的解决方案与下一个问题足够接近，步长必须足够小。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::compute_initial_guess(double step_size) 
 *   { 
 *     const double target_Re = 1.0 / viscosity; 
 * 
 *     bool is_initial_step = true; 
 * 
 *     for (double Re = 1000.0; Re < target_Re; 
 *          Re        = std::min(Re + step_size, target_Re)) 
 *       { 
 *         viscosity = 1.0 / Re; 
 *         std::cout << "Searching for initial guess with Re = " << Re 
 *                   << std::endl; 
 *         newton_iteration(1e-12, 50, 0, is_initial_step, false); 
 *         is_initial_step = false; 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesoutput_results"></a> 
 * <h4>StationaryNavierStokes::output_results</h4>
 * 

 * 
 * 这个函数与 step-22 中的函数相同，只是我们为输出文件选择了一个同时包含雷诺数（即当前环境下的粘度的倒数）的名称。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::output_results( 
 *     const unsigned int output_index) const 
 *   { 
 *     std::vector<std::string> solution_names(dim, "velocity"); 
 *     solution_names.emplace_back("pressure"); 
 * 
 *     std::vector<DataComponentInterpretation::DataComponentInterpretation> 
 *       data_component_interpretation( 
 *         dim, DataComponentInterpretation::component_is_part_of_vector); 
 *     data_component_interpretation.push_back( 
 *       DataComponentInterpretation::component_is_scalar); 
 *     DataOut<dim> data_out; 
 *     data_out.attach_dof_handler(dof_handler); 
 *     data_out.add_data_vector(present_solution, 
 *                              solution_names, 
 *                              DataOut<dim>::type_dof_data, 
 *                              data_component_interpretation); 
 *     data_out.build_patches(); 
 * 
 *     std::ofstream output(std::to_string(1.0 / viscosity) + "-solution-" + 
 *                          Utilities::int_to_string(output_index, 4) + ".vtk"); 
 *     data_out.write_vtk(output); 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesprocess_solution"></a> 
 * <h4>StationaryNavierStokes::process_solution</h4>
 * 

 * 
 * 在我们的测试案例中，我们不知道分析解。该函数输出沿 $x=0.5$ 和 $0 \leq y \leq 1$ 的速度分量，以便与文献中的数据进行比较。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::process_solution(unsigned int refinement) 
 *   { 
 *     std::ofstream f(std::to_string(1.0 / viscosity) + "-line-" + 
 *                     std::to_string(refinement) + ".txt"); 
 *     f << "# y u_x u_y" << std::endl; 
 * 
 *     Point<dim> p; 
 *     p(0) = 0.5; 
 *     p(1) = 0.5; 
 * 
 *     f << std::scientific; 
 * 
 *     for (unsigned int i = 0; i <= 100; ++i) 
 *       { 
 *         p(dim - 1) = i / 100.0; 
 * 
 *         Vector<double> tmp_vector(dim + 1); 
 *         VectorTools::point_value(dof_handler, present_solution, p, tmp_vector); 
 *         f << p(dim - 1); 
 * 
 *         for (int j = 0; j < dim; j++) 
 *           f << " " << tmp_vector(j); 
 *         f << std::endl; 
 *       } 
 *   } 
 * @endcode
 * 
 * 
 * <a name="StationaryNavierStokesrun"></a> 
 * <h4>StationaryNavierStokes::run</h4>
 * 

 * 
 * 这是本程序的最后一步。在这一部分，我们分别生成网格和运行其他函数。最大细化度可以通过参数来设置。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void StationaryNavierStokes<dim>::run(const unsigned int refinement) 
 *   { 
 *     GridGenerator::hyper_cube(triangulation); 
 *     triangulation.refine_global(5); 
 * 
 *     const double Re = 1.0 / viscosity; 
 * 
 * @endcode
 * 
 * 如果粘度小于 $1/1000$ ，我们必须首先通过延续法搜索初始猜测。我们应该注意的是，搜索总是在初始网格上进行的，也就是这个程序中的 $8 \times 8$ 网格。之后，我们只需做与粘度大于 $1/1000$ 时相同的工作：运行牛顿迭代，细化网格，转移解决方案，并重复。
 * 

 * 
 * 
 * @code
 *     if (Re > 1000.0) 
 *       { 
 *         std::cout << "Searching for initial guess ..." << std::endl; 
 *         const double step_size = 2000.0; 
 *         compute_initial_guess(step_size); 
 *         std::cout << "Found initial guess." << std::endl; 
 *         std::cout << "Computing solution with target Re = " << Re << std::endl; 
 *         viscosity = 1.0 / Re; 
 *         newton_iteration(1e-12, 50, refinement, false, true); 
 *       } 
 *     else 
 *       { 
 * 
 * @endcode
 * 
 * 当粘度大于1/1000时，斯托克斯方程的解作为初始猜测已经足够好。如果是这样，我们就不需要用延续法来搜索初始猜测了。牛顿迭代可以直接开始。
 * 

 * 
 * 
 * @code
 *         newton_iteration(1e-12, 50, refinement, true, true); 
 *       } 
 *   } 
 * } // namespace Step57 
 * 
 * int main() 
 * { 
 *   try 
 *     { 
 *       using namespace Step57; 
 * 
 *       StationaryNavierStokes<2> flow(/* degree = */ 
 * 
 * 1); 
 *       flow.run(4); 
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
 *  
 *  
 *  
 *  
 *  
 *  
 *                 << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       std::cerr << "Unknown exception!" << std::endl 
 *                 << "Aborting!" << std::endl 
 *                 << "----------------------------------------------------" 
 *                 << std::endl; 
 *       return 1; 
 *     } 
 *   return 0; 
 * } 
 * 
 * 
 * @endcode
examples/step-57/doc/results.dox



<a name="Results"></a><h1>Results</h1>


现在我们用上面讨论的方法来解决有粘度的纳维尔-斯托克斯方程  $1/400$  和  $1/7500$  。

<a name="Testcase1LowReynoldsNumber"></a><h3> Test case 1: Low Reynolds Number </h3>


在第一个测试案例中，粘度被设定为 $1/400$  。正如我们在介绍中所讨论的，初始猜测是相应斯托克斯问题的解。在下面的表格中，显示了每个网格上的牛顿迭代的残差。表中的数据表明，牛顿迭代的收敛性是四等的。

 <table align="center" class="doxtable">
<tr>
    <th>$\mathrm{Re}=400$</th>
    <th colspan="2">Mesh0</th>
    <th colspan="2">Mesh1</th>
    <th colspan="2">Mesh2</th>
    <th colspan="2">Mesh3</th>
    <th colspan="2">Mesh4</th>
</tr>
<tr>
    <th>Newton iter   </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
</tr>
<tr>
  <td>1</td>
  <td>3.7112e-03</td>
  <td>5</td>
  <td>6.4189e-03</td>
  <td>3</td>
  <td>2.4338e-03</td>
  <td>3</td>
  <td>1.0570e-03</td>
  <td>3</td>
  <td>4.9499e-04</td>
  <td>3</td>
</tr>
<tr>
  <td>2</td>
  <td>7.0849e-04</td>
  <td>5</td>
  <td>9.9458e-04</td>
  <td>5</td>
  <td>1.1409e-04</td>
  <td>6</td>
  <td>1.3544e-05</td>
  <td>6</td>
  <td>1.4171e-06</td>
  <td>6</td>
</tr>
<tr>
  <td>3</td>
  <td>1.9980e-05</td>
  <td>5</td>
  <td>4.5007e-05</td>
  <td>5</td>
  <td>2.9020e-08</td>
  <td>5</td>
  <td>4.4021e-10</td>
  <td>6</td>
  <td>6.3435e-11</td>
  <td>6</td>
</tr>
<tr>
  <td>4</td>
  <td>2.3165e-09</td>
  <td>6</td>
  <td>1.6891e-07</td>
  <td>5</td>
  <td>1.2338e-14</td>
  <td>7</td>
  <td>1.8506e-14</td>
  <td>8</td>
  <td>8.8563e-15</td>
  <td>8</td>
</tr>
<tr>
  <td>5</td>
  <td>1.2585e-13</td>
  <td>7</td>
  <td>1.4520e-11</td>
  <td>6</td>
  <td>1.9044e-13</td>
  <td>8</td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
</tr>
<tr>
  <td>6</td>
  <td></td>
  <td></td>
  <td>1.3998e-15</td>
  <td>8</td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
</tr>
</table> 








下面的数字显示了生成网格的顺序。对于 $\mathrm{Re}=400$ 的情况，初始猜测是通过在 $8 \times 8$ 网格上求解斯托克斯得到的，并且网格是自适应细化的。在不同的网格之间，粗网格的解被内插到细网格中，作为初始猜测使用。

 <table align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Mesh0.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Mesh1.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Mesh2.png" width="232px" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Mesh3.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Mesh4.png" width="232px" alt="">
    </td>
  </tr>
</table> 

这张图片是用  $\mathrm{Re}=400$  的盖子驱动的空腔的图形流线结果。   <img src="https://www.dealii.org/images/steps/developer/step-57.Re400_Streamline.png" alt=""> 

然后将该方案与来自[4]的参考方案进行比较，参考方案的数据可以在文件 "ref_2d_ghia_u.txt "中找到。

 <img src="https://www.dealii.org/images/steps/developer/step-57.compare-Re400.svg" style="width:50%" alt=""> 

<a name="Testcase2HighReynoldsNumber"></a><h3> Test case 2: High Reynolds Number </h3>


牛顿迭代法需要一个良好的初始猜测。然而，当雷诺数很大时，非线性项占主导地位，因此，斯托克斯方程的解可能与精确解相去甚远。如果斯托克斯方程的解作为初始猜测，收敛性就会丧失。下图显示，非线性迭代被卡住了，在进一步的迭代中，残差不再减少。

 <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_loss_convergence.svg" style="width:50%" alt=""> 

因此，初始猜测必须通过延续法获得，这在导言中已经讨论过了。这里延续法的步长是 $|\nu_{i}-\nu_{i+1}|$ ，是2000，初始网格的大小是 $32 \times 32$ 。在得到一个初始猜测后，如前一个测试案例一样，对网格进行细化。下图显示，在每次细化时，牛顿迭代都有二次收敛性。为解决这个测试案例，共执行了52步牛顿迭代。

 <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_get_convergence.svg" style="width:50%" alt=""> 

我们还显示了每个网格上牛顿迭代的每一步的残差。表格中可以清楚地看到二次收敛的情况。

 <table align="center" class="doxtable">
  <tr>
    <th>$\mathrm{Re}=7500$</th>
    <th colspan="2">Mesh0</th>
    <th colspan="2">Mesh1</th>
    <th colspan="2">Mesh2</th>
    <th colspan="2">Mesh3</th>
    <th colspan="2">Mesh4</th>
  </tr>
  <tr>
    <th>Newton iter   </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
    <th>Residual      </th>
    <th>FGMRES        </th>
  </tr>
<tr>
  <td>1</td>
  <td>1.8922e-06</td>
  <td>6</td>
  <td>4.2506e-03</td>
  <td>3</td>
  <td>1.4299e-03</td>
  <td>3</td>
  <td>4.8793e-04</td>
  <td>2</td>
  <td>1.8998e-04</td>
  <td>2</td>
</tr>
<tr>
  <td>2</td>
  <td>3.1644e-09</td>
  <td>8</td>
  <td>1.3732e-03</td>
  <td>7</td>
  <td>4.1506e-04</td>
  <td>7</td>
  <td>9.1119e-05</td>
  <td>8</td>
  <td>1.3555e-05</td>
  <td>8</td>
</tr>
<tr>
  <td>3</td>
  <td>1.7611e-14</td>
  <td>9</td>
  <td>2.1946e-04</td>
  <td>6</td>
  <td>1.7881e-05</td>
  <td>6</td>
  <td>5.2678e-07</td>
  <td>7</td>
  <td>9.3739e-09</td>
  <td>7</td>
</tr>
<tr>
  <td>4</td>
  <td></td>
  <td></td>
  <td>8.8269e-06</td>
  <td>6</td>
  <td>6.8210e-09</td>
  <td>7</td>
  <td>2.2770e-11</td>
  <td>8</td>
  <td>1.2588e-13</td>
  <td>9</td>
</tr>
<tr>
  <td>5</td>
  <td></td>
  <td></td>
  <td>1.2974e-07</td>
  <td>7</td>
  <td>1.2515e-13</td>
  <td>9</td>
  <td>1.7801e-14</td>
  <td>1</td>
  <td></td>
  <td></td>
</tr>
<tr>
  <td>6</td>
  <td></td>
  <td></td>
  <td>4.4352e-11</td>
  <td>7</td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
</tr>
<tr>
  <td>7</td>
  <td></td>
  <td></td>
  <td>6.2863e-15</td>
  <td>9</td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
  <td></td>
</tr>
</table> 








生成的网格序列看起来像这样。   <table align="center">
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Mesh0.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Mesh1.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Mesh2.png" width="232px" alt="">
    </td>
  </tr>
  <tr>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Mesh3.png" width="232px" alt="">
    </td>
    <td align="center">
      <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Mesh4.png" width="232px" alt="">
    </td>
  </tr>
</table>  我们将我们的解决方案与[5]的参考解决方案进行比较。   <img src="https://www.dealii.org/images/steps/developer/step-57.compare-Re7500.svg" style="width:50%" alt="">  下图是图形结果。   <img src="https://www.dealii.org/images/steps/developer/step-57.Re7500_Streamline.png" alt=""> 

此外，误差由非线性误差和离散化误差组成，前者随着我们进行牛顿迭代而减少，后者则取决于网格大小。这就是为什么我们必须细化网格并在下一个更细的网格上重复牛顿迭代。从上表中，我们可以看到每个网格上的残差（非线性误差）都低于1e-12，但下图向我们展示了随后更细的网格上的解决方案之间的差异。

 <img src="https://www.dealii.org/images/steps/developer/step-57.converge-Re7500.svg" style="width:50%" alt=""> 

<a name="extensions"></a>

<a name="Possibilitiesforextensions"></a><h3>Possibilities for extensions</h3>


<a name="Comparetoothersolvers"></a><h4>Compare to other solvers</h4>


比较目前实现的线性求解器和仅仅使用UMFPACK求解整个线性系统是很容易的。你需要去除包含常压的空域，这在第56步中完成。更有趣的是与其他先进的预处理程序如PCD的比较。事实证明，这里的预处理器是非常有竞争力的，在论文[2]中可以看到。

下表显示了我们的迭代方法（FGMRES）与直接求解器（UMFPACK）之间在粘度设置为1/400的情况下对整个系统的计时结果。尽管我们在迭代求解器中的速度块使用了相同的直接求解器，但它的速度要快得多，消耗的内存也少。这在三维中会更加明显。

 <table align="center" class="doxtable">
<tr>
  <th>Refinement Cycle</th>
  <th>DoFs</th>
  <th>Iterative: Total/s (Setup/s)</th>
  <th>Direct: Total/s (Setup/s)</th>
</tr>
<tr>
  <td>5</td>
  <td>9539</td>
  <td>0.10 (0.06)</td>
  <td>0.13 (0.12)</td>
</tr>
<tr>
  <td>6</td>
  <td>37507</td>
  <td>0.58 (0.37)</td>
  <td>1.03 (0.97)</td>
</tr>
<tr>
  <td>7</td>
  <td>148739</td>
  <td>3.59 (2.73)</td>
  <td>7.78 (7.53)</td>
</tr>
<tr>
  <td>8</td>
  <td>592387</td>
  <td>29.17 (24.94)</td>
  <td>(>4GB RAM)</td>
</tr>
</table> 




<a name="3dcomputations"></a><h4>3d computations</h4>


该代码被设置为也可以在3D中运行。当然，参考值是不同的，例如，见[6]。在这个例子中，高分辨率的计算是无法做到的，因为速度块的直接求解器在三维中不能很好地工作。相反，需要一个基于代数或几何多网格的并行求解器。见下文。

<a name="Parallelization"></a><h4>Parallelization</h4>


对于较大的计算，特别是三维计算，有必要实现MPI并行求解器和预处理器。一个好的起点是step-55，它对斯托克斯方程的速度块使用代数多重网格。另一个选择是看一下<a href="https://www.dealii.org/code-gallery.html">deal.II code
gallery</a>中的代码列表，其中已经包含了并行的纳维-斯托克斯求解器。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-57.cc"
*/
