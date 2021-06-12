/**
@page step_5 The step-5 tutorial program
This tutorial depends on step-4.

@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
        <li><a href="#Includefiles">Include files</a>
        <li><a href="#ThecodeStep5codeclasstemplate">The <code>Step5</code> class template</a>
        <li><a href="#Workingwithnonconstantcoefficients">Working with nonconstant coefficients</a>
        <li><a href="#ThecodeStep5codeclassimplementation">The <code>Step5</code> class implementation</a>
      <ul>
        <li><a href="#Step5Step5">Step5::Step5</a>
        <li><a href="#Step5setup_system">Step5::setup_system</a>
        <li><a href="#Step5assemble_system">Step5::assemble_system</a>
        <li><a href="#Step5solve">Step5::solve</a>
        <li><a href="#Step5output_resultsandsettingoutputflags">Step5::output_results and setting output flags</a>
        <li><a href="#Step5run">Step5::run</a>
      </ul>
        <li><a href="#Thecodemaincodefunction">The <code>main</code> function</a>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-5/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>


 @dealiiVideoLecture{14} 

这个例子并没有显示出革命性的新东西，但它显示了许多比以前的例子更小的改进，也显示了许多通常可以在有限元程序中发现的小东西。其中包括。   <ul>   <li>  在连续细化的网格上进行计算。至少在数学科学中，在一个层次的网格上计算解是很常见的，以便对解的精度有一个感觉；如果你在一个网格上只有一个解，你通常无法猜测解的精度。此外，deal.II被设计用来支持自适应算法，其中在连续细化的网格上的迭代求解是算法的核心。虽然在这个例子中没有使用自适应网格，但这里为它们奠定了基础。     <li>  在实际应用中，领域经常被自动网格生成器细分为三角形。为了使用它们，从文件中读取粗大的网格是很重要的。在这个例子中，我们将读取一个UCD（非结构化单元数据）格式的粗网格。当这个程序在2000年左右首次编写时，UCD格式是AVS Explorer所使用的--这个程序在当时被合理地广泛使用，但现在已经不再重要了。尽管如此，该文件格式仍然存在，并且仍然被一些程序所理解）。     <li>  有限元程序通常会使用大量的计算时间，所以有时需要进行一些优化。我们将展示其中的一些。     <li>  另一方面，有限元程序往往是相当复杂的，所以调试是一个重要方面。我们通过使用断言来支持安全编程，断言在调试模式下检查参数和%内部状态的有效性，但在优化模式下被删除。(  @dealiiVideoLectureSeeAlso{18})   <li>  关于数学方面，我们展示了如何支持椭圆算子中的可变系数，以及如何对线性方程组使用预处理迭代求解器。   </ul> 

这里要解决的方程式如下。

@f{align*}


  -\nabla \cdot a(\mathbf x) \nabla u(\mathbf x) &= 1 \qquad\qquad & \text{in}\ \Omega,
  \\
  u &= 0 \qquad\qquad & \text{on}\ \partial\Omega.


@f}

如果 $a(\mathbf x)$ 是一个恒定的系数，这就只是泊松方程了。然而，如果它确实是空间可变的，它就是一个更复杂的方程（通常被称为 "扩展泊松方程"）。根据变量 $u$ 所指的内容，它可以模拟各种情况，具有广泛的适用性。

- 如果 $u$ 是电动势，那么 $-a\nabla u$ 是介质中的电流，系数 $a$ 是介质在任何特定点的电导率。在这种情况下，方程的右侧将是电源密度，通常为零或由局部的、类似德尔塔的函数组成）。

- 如果 $u$ 是薄膜的垂直挠度，那么 $a$ 将是对局部刚度的测量。这就是让我们解释下面结果部分所显示的图像的解释。

由于拉普拉斯/泊松方程出现在如此多的场合中，因此除了上面列出的两种解释外，还有许多其他解释。

当组装这个方程的线性系统时，我们需要弱的形式，这里的内容如下。

@f{align*}
  (a \nabla \varphi, \nabla u) &= (\varphi, 1) \qquad \qquad \forall \varphi.


@f}

 <code>assemble_system</code> 函数中的实现紧随其后。


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
 * 同样，前几个include文件已经知道了，所以我们不会对它们进行评论。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/function.h> 
 * #include <deal.II/base/logstream.h> 
 * #include <deal.II/lac/vector.h> 
 * #include <deal.II/lac/full_matrix.h> 
 * #include <deal.II/lac/sparse_matrix.h> 
 * #include <deal.II/lac/dynamic_sparsity_pattern.h> 
 * #include <deal.II/lac/solver_cg.h> 
 * #include <deal.II/lac/precondition.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/dofs/dof_tools.h> 
 * #include <deal.II/fe/fe_q.h> 
 * #include <deal.II/fe/fe_values.h> 
 * #include <deal.II/numerics/vector_tools.h> 
 * #include <deal.II/numerics/matrix_tools.h> 
 * #include <deal.II/numerics/data_out.h> 
 * 
 * @endcode
 * 
 * 这个是新的。我们想从磁盘上读取一个三角图，做这个的类在下面的文件中声明。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/grid_in.h> 
 * 
 * @endcode
 * 
 * 我们将使用一个圆形域，而描述其边界的对象来自这个文件。
 * 

 * 
 * 
 * @code
 * #include <deal.II/grid/manifold_lib.h> 
 * 
 * @endcode
 * 
 * 这是C++ ...
 * 

 * 
 * 
 * @code
 * #include <fstream> 
 * #include <iostream> 
 * 
 * @endcode
 * 
 * 最后，这在以前的教程程序中已经讨论过了。
 * 

 * 
 * 
 * @code
 * using namespace dealii; 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep5codeclasstemplate"></a> 
 * <h3>The <code>Step5</code> class template</h3>
 * 

 * 
 * 主类大部分和前面的例子一样。最明显的变化是删除了 <code>make_grid</code> 函数，因为现在创建网格是在 <code>run</code> 函数中完成的，其余功能都在 <code>setup_system</code> 中。除此以外，一切都和以前一样。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * class Step5 
 * { 
 * public: 
 *   Step5(); 
 *   void run(); 
 * 
 * private: 
 *   void setup_system(); 
 *   void assemble_system(); 
 *   void solve(); 
 *   void output_results(const unsigned int cycle) const; 
 * 
 *   Triangulation<dim> triangulation; 
 *   FE_Q<dim>          fe; 
 *   DoFHandler<dim>    dof_handler; 
 * 
 *   SparsityPattern      sparsity_pattern; 
 *   SparseMatrix<double> system_matrix; 
 * 
 *   Vector<double> solution; 
 *   Vector<double> system_rhs; 
 * }; 
 * @endcode
 * 
 * 
 * <a name="Workingwithnonconstantcoefficients"></a> 
 * <h3>Working with nonconstant coefficients</h3>
 * 

 * 
 * 在  step-4  中，我们展示了如何使用非恒定边界值和右手边。 在这个例子中，我们想在椭圆算子中使用一个可变系数来代替。由于我们有一个只取决于空间中的点的函数，我们可以做得更简单一些，使用一个普通的函数而不是继承自Function。
 * 

 * 
 * 这是对单点的系数函数的实现。如果与原点的距离小于0.5，我们让它返回20，否则返回1。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * double coefficient(const Point<dim> &p) 
 * { 
 *   if (p.square() < 0.5 * 0.5) 
 *     return 20; 
 *   else 
 *     return 1; 
 * } 
 * @endcode
 * 
 * 
 * <a name="ThecodeStep5codeclassimplementation"></a> 
 * <h3>The <code>Step5</code> class implementation</h3>
 * 
 * <a name="Step5Step5"></a> 
 * <h4>Step5::Step5</h4>
 * 

 * 
 * 这个函数和以前一样。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * Step5<dim>::Step5() 
 *   : fe(1) 
 *   , dof_handler(triangulation) 
 * {} 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5setup_system"></a> 
 * <h4>Step5::setup_system</h4>
 * 

 * 
 * 这是前面例子中的函数 <code>make_grid</code> ，减去了生成网格的部分。其他一切都没有变化。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step5<dim>::setup_system() 
 * { 
 *   dof_handler.distribute_dofs(fe); 
 * 
 *   std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() 
 *             << std::endl; 
 * 
 *   DynamicSparsityPattern dsp(dof_handler.n_dofs()); 
 *   DoFTools::make_sparsity_pattern(dof_handler, dsp); 
 *   sparsity_pattern.copy_from(dsp); 
 * 
 *   system_matrix.reinit(sparsity_pattern); 
 * 
 *   solution.reinit(dof_handler.n_dofs()); 
 *   system_rhs.reinit(dof_handler.n_dofs()); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5assemble_system"></a> 
 * <h4>Step5::assemble_system</h4>
 * 

 * 
 * 和前面的例子一样，这个函数在功能上没有太大变化，但仍有一些优化，我们将展示这些优化。对此，需要注意的是，如果使用高效的求解器（如预设条件的CG方法），组装矩阵和右手边会花费相当的时间，你应该考虑在某些地方使用一到两个优化。
 * 

 * 
 * 该函数的前几部分与之前完全没有变化。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step5<dim>::assemble_system() 
 * { 
 *   QGauss<dim> quadrature_formula(fe.degree + 1); 
 * 
 *   FEValues<dim> fe_values(fe, 
 *                           quadrature_formula, 
 *                           update_values | update_gradients | 
 *                             update_quadrature_points | update_JxW_values); 
 * 
 *   const unsigned int dofs_per_cell = fe.n_dofs_per_cell(); 
 * 
 *   FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell); 
 *   Vector<double>     cell_rhs(dofs_per_cell); 
 * 
 *   std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell); 
 * 
 * @endcode
 * 
 * 接下来是对所有单元的典型循环，以计算局部贡献，然后将它们转移到全局矩阵和向量中。与 step-4 相比，这部分的唯一变化是我们将使用上面定义的 <code>coefficient()</code> 函数来计算每个正交点的系数值。
 * 

 * 
 * 
 * @code
 *   for (const auto &cell : dof_handler.active_cell_iterators()) 
 *     { 
 *       cell_matrix = 0.; 
 *       cell_rhs    = 0.; 
 * 
 *       fe_values.reinit(cell); 
 * 
 *       for (const unsigned int q_index : fe_values.quadrature_point_indices()) 
 *         { 
 *           const double current_coefficient = 
 *             coefficient(fe_values.quadrature_point(q_index)); 
 *           for (const unsigned int i : fe_values.dof_indices()) 
 *             { 
 *               for (const unsigned int j : fe_values.dof_indices()) 
 *                 cell_matrix(i, j) += 
 *                   (current_coefficient *              // a(x_q) 
 *                    fe_values.shape_grad(i, q_index) * // grad phi_i(x_q) 
 *                    fe_values.shape_grad(j, q_index) * // grad phi_j(x_q) 
 *                    fe_values.JxW(q_index));           // dx 
 * 
 *               cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q) 
 *                               1.0 *                               // f(x_q) 
 *                               fe_values.JxW(q_index));            // dx 
 *             } 
 *         } 
 * 
 *       cell->get_dof_indices(local_dof_indices); 
 *       for (const unsigned int i : fe_values.dof_indices()) 
 *         { 
 *           for (const unsigned int j : fe_values.dof_indices()) 
 *             system_matrix.add(local_dof_indices[i], 
 *                               local_dof_indices[j], 
 *                               cell_matrix(i, j)); 
 * 
 *           system_rhs(local_dof_indices[i]) += cell_rhs(i); 
 *         } 
 *     } 
 * 
 * @endcode
 * 
 * 有了这样构建的矩阵，我们再次使用零边界值。
 * 

 * 
 * 
 * @code
 *   std::map<types::global_dof_index, double> boundary_values; 
 *   VectorTools::interpolate_boundary_values(dof_handler, 
 *                                            0, 
 *                                            Functions::ZeroFunction<dim>(), 
 *                                            boundary_values); 
 *   MatrixTools::apply_boundary_values(boundary_values, 
 *                                      system_matrix, 
 *                                      solution, 
 *                                      system_rhs); 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step5solve"></a> 
 * <h4>Step5::solve</h4>
 * 

 * 
 * 求解过程看起来又和前面的例子差不多。然而，我们现在将使用一个预设条件的共轭梯度算法。做出这种改变并不难。事实上，我们唯一需要改变的是，我们需要一个作为预处理程序的对象。我们将使用SSOR（对称连续过度放松），放松系数为1.2。为此， <code>SparseMatrix</code> 类有一个函数可以做一个SSOR步骤，我们需要把这个函数的地址和它应该作用的矩阵（也就是要反转的矩阵）以及松弛因子打包成一个对象。 <code>PreconditionSSOR</code> 类为我们做了这个。(  <code>PreconditionSSOR</code>  类需要一个模板参数，表示它应该工作的矩阵类型。默认值是 <code>SparseMatrix@<double@></code> ，这正是我们在这里需要的，所以我们只需坚持使用默认值，不在角括号中指定任何东西。)
 * 

 * 
 * 请注意，在目前的情况下，SSOR的表现并不比其他大多数预处理程序好多少（尽管比没有预处理好）。在下一个教程程序  step-6  的结果部分，将对不同的预处理进行简要比较。
 * 

 * 
 * 有了这个，函数的其余部分就很简单了：我们现在使用我们声明的预处理程序，而不是之前创建的 <code>PreconditionIdentity</code> 对象，CG求解器将为我们完成其余的工作。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step5<dim>::solve() 
 * { 
 *   SolverControl            solver_control(1000, 1e-12); 
 *   SolverCG<Vector<double>> solver(solver_control); 
 * 
 *   PreconditionSSOR<SparseMatrix<double>> preconditioner; 
 *   preconditioner.initialize(system_matrix, 1.2); 
 * 
 *   solver.solve(system_matrix, solution, system_rhs, preconditioner); 
 * 
 *   std::cout << "   " << solver_control.last_step() 
 *             << " CG iterations needed to obtain convergence." << std::endl; 
 * } 
 * @endcode
 * 
 * 
 * <a name="Step5output_resultsandsettingoutputflags"></a> 
 * <h4>Step5::output_results and setting output flags</h4>
 * 

 * 
 * 将输出写入文件的方法与上一个教程中的基本相同。唯一不同的是，我们现在需要为每个细化周期构建一个不同的文件名。
 * 

 * 
 * 这个函数以VTU格式写入输出，这是VTK格式的一个变种，因为它压缩了数据，所以需要更少的磁盘空间。当然，如果你希望使用一个不理解VTK或VTU的可视化程序，DataOut类还支持许多其他格式。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step5<dim>::output_results(const unsigned int cycle) const 
 * { 
 *   DataOut<dim> data_out; 
 * 
 *   data_out.attach_dof_handler(dof_handler); 
 *   data_out.add_data_vector(solution, "solution"); 
 * 
 *   data_out.build_patches(); 
 * 
 *   std::ofstream output("solution-" + std::to_string(cycle) + ".vtu"); 
 *   data_out.write_vtu(output); 
 * } 
 * 
 * @endcode
 * 
 * 
 * <a name="Step5run"></a> 
 * <h4>Step5::run</h4>
 * 

 * 
 * 在这个程序中，倒数第二件事是对 <code>run()</code> 函数的定义。与之前的程序不同，我们将在一个网格序列上进行计算，在每次迭代后都会进行全局细化。因此，该函数由6个周期的循环组成。在每个循环中，我们首先打印循环编号，然后决定如何处理网格。如果这不是第一个周期，我们就简单地对现有的网格进行一次全局精炼。然而，在运行这些循环之前，我们必须先生成一个网格。
 * 

 * 
 * 在前面的例子中，我们已经使用了 <code>GridGenerator</code> 类中的一些函数。在这里，我们想从一个存储单元的文件中读取网格，这个文件可能来自其他人，也可能是一个网格生成工具的产物。
 * 

 * 
 * 为了从文件中读取网格，我们生成一个数据类型为GridIn的对象，并将三角剖分与之相关联（也就是说，当我们要求它读取文件时，我们告诉它要填充我们的三角剖分对象）。然后我们打开相应的文件，用文件中的数据初始化三角剖分。
 * 

 * 
 * 
 * @code
 * template <int dim> 
 * void Step5<dim>::run() 
 * { 
 *   GridIn<dim> grid_in; 
 *   grid_in.attach_triangulation(triangulation); 
 *   std::ifstream input_file("circle-grid.inp"); 
 * 
 * @endcode
 * 
 * 我们现在想读取该文件。但是，输入文件只针对二维三角图，而这个函数是一个任意维度的模板。由于这只是一个演示程序，我们不会为不同的维度使用不同的输入文件，而是在不在二维的情况下迅速杀死整个程序。当然，由于下面的主函数假定我们是在二维空间工作，我们可以跳过这个检查，在这个版本的程序中，不会有任何不良影响。
 * 

 * 
 * 事实证明，90%以上的编程错误都是无效的函数参数，如无效的数组大小等，所以我们在整个deal.II中大量使用断言来捕捉此类错误。为此， <code>Assert</code> 宏是一个很好的选择，因为它确保作为第一个参数的条件是有效的，如果不是，就抛出一个异常（它的第二个参数），通常会终止程序，并给出错误发生的位置和原因的信息。关于 @p Assert 宏的具体作用，可以在 @ref Exceptions "异常文档模块 "中找到更详细的讨论）。这通常会大大减少发现编程错误的时间，我们发现断言是快速编程的宝贵手段。
 * 

 * 
 * 另一方面，如果你想做大的计算，所有这些检查（目前库中有超过10000个）不应该使程序太慢。为此， <code>Assert</code> 宏只在调试模式下使用，如果在优化模式下则扩展为零。因此，当你在小问题上测试你的程序并进行调试时，断言会告诉你问题出在哪里。一旦你的程序稳定了，你可以关闭调试，程序将在没有断言的情况下以最大速度运行你的实际计算。更准确地说：通过在优化模式下编译你的程序，关闭库中的所有检查（这些检查可以防止你用错误的参数调用函数，从数组中走出来，等等），通常可以使程序的运行速度提高四倍左右。即使优化后的程序性能更高，我们仍然建议在调试模式下开发，因为它允许库自动发现许多常见的编程错误。对于那些想尝试的人来说。从调试模式切换到优化模式的方法是用<code>make release</code>命令重新编译你的程序。现在 <code>make</code> 程序的输出应该向你表明，该程序现在是以优化模式编译的，以后也会被链接到已经为优化模式编译的库。为了切换回调试模式，只需用  <code>make debug</code>  命令重新编译。
 * 

 * 
 * 
 * @code
 *   Assert(dim == 2, ExcInternalError()); 
 * 
 * @endcode
 * 
 * ExcInternalError是一个全局定义的异常，每当出现严重的错误时就会抛出。通常，人们希望使用更具体的异常，特别是在这种情况下，如果 <code>dim</code> 不等于2，人们当然会尝试做其他事情，例如使用库函数创建一个网格。终止程序通常不是一个好主意，断言实际上只应该用于不应该发生的特殊情况，但由于程序员、用户或其他人的愚蠢而可能发生。上面的情况并不是对Assert的巧妙使用，但是再次强调：这是一个教程，也许值得展示一下什么是不应该做的，毕竟。
 * 

 * 
 * 所以，如果我们通过了断言，我们就知道dim==2，现在我们就可以真正地读取网格。它的格式是UCD（非结构化单元数据）（尽管惯例是使用UCD文件的后缀 <code>inp</code> ）。
 * 

 * 
 * 
 * @code
 *   grid_in.read_ucd(input_file); 
 * 
 * @endcode
 * 
 * 如果你想使用其他输入格式，你必须使用其他 <code>grid_in.read_xxx</code> 函数之一。(参见  <code>GridIn</code>  类的文档，以了解目前支持哪些输入格式)。
 * 

 * 
 * 文件中的网格描述了一个圆。因此，我们必须使用一个流形对象，告诉三角计算在细化网格时将边界上的新点放在哪里。与 step-1 不同的是，由于GridIn不知道域的边界是圆形的（与 GridGenerator::hyper_shell) 不同的是，我们必须在创建三角网格后明确地将流形附加到边界上，以便在细化网格时获得正确的结果。
 * 

 * 
 * 
 * @code
 *   const SphericalManifold<dim> boundary; 
 *   triangulation.set_all_manifold_ids_on_boundary(0); 
 *   triangulation.set_manifold(0, boundary); 
 * 
 *   for (unsigned int cycle = 0; cycle < 6; ++cycle) 
 *     { 
 *       std::cout << "Cycle " << cycle << ':' << std::endl; 
 * 
 *       if (cycle != 0) 
 *         triangulation.refine_global(1); 
 * 
 * @endcode
 * 
 * 现在我们有了一个确定的网格，我们写一些输出，做所有我们在前面的例子中已经看到的事情。
 * 

 * 
 * 
 * @code
 *       std::cout << "   Number of active cells: "  // 
 *                 << triangulation.n_active_cells() // 
 *                 << std::endl                      // 
 *                 << "   Total number of cells: "   // 
 *                 << triangulation.n_cells()        // 
 *                 << std::endl; 
 * 
 *       setup_system(); 
 *       assemble_system(); 
 *       solve(); 
 *       output_results(cycle); 
 *     } 
 * } 
 * @endcode
 * 
 * 
 * <a name="Thecodemaincodefunction"></a> 
 * <h3>The <code>main</code> function</h3>
 * 

 * 
 * 主函数看起来和前面的例子中的函数差不多，所以我们就不进一步评论了。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   Step5<2> laplace_problem_2d; 
 *   laplace_problem_2d.run(); 
 *   return 0; 
 * } 
 * 
 * @endcode
examples/step-5/doc/results.dox



<a name="Results"></a><h1>Results</h1>



下面是控制台的输出。

@code
Cycle 0:
   Number of active cells: 20
   Total number of cells: 20
   Number of degrees of freedom: 25
   13 CG iterations needed to obtain convergence.
Cycle 1:
   Number of active cells: 80
   Total number of cells: 100
   Number of degrees of freedom: 89
   18 CG iterations needed to obtain convergence.
Cycle 2:
   Number of active cells: 320
   Total number of cells: 420
   Number of degrees of freedom: 337
   29 CG iterations needed to obtain convergence.
Cycle 3:
   Number of active cells: 1280
   Total number of cells: 1700
   Number of degrees of freedom: 1313
   52 CG iterations needed to obtain convergence.
Cycle 4:
   Number of active cells: 5120
   Total number of cells: 6820
   Number of degrees of freedom: 5185
   95 CG iterations needed to obtain convergence.
Cycle 5:
   Number of active cells: 20480
   Total number of cells: 27300
   Number of degrees of freedom: 20609
   182 CG iterations needed to obtain convergence.
@endcode






在每个周期中，单元格的数量翻了两番，CG迭代的数量大约翻了一番。另外，在每个周期中，程序写出一个VTU格式的输出图形文件。它们被描述为以下内容。

 <table width="100%">
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-0-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-1-r9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-2-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-3-r9.2.png" alt="">
    </td>
  </tr>
  <tr>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-4-r9.2.png" alt="">
    </td>
    <td>
      <img src="https://www.dealii.org/images/steps/developer/step-5.solution-5-r9.2.png" alt="">
    </td>
  </tr>
</table> 




由于系数的可变性（那里的曲率减少的程度与系数增加的程度相同），溶液的顶部区域被压扁了。溶液的梯度沿着界面是不连续的，尽管这在上面的图片中不是很明显。我们将在下一个例子中更详细地研究这个问题。

图片还显示，这个程序计算出来的解在非常粗的网格上其实是相当错误的（它的大小是错误的）。这是因为没有任何数值方法能保证粗大网格上的解是特别准确的--但我们知道解<i>converges</i>是精确的解，事实上你可以看到从一个网格到下一个网格的解似乎在最后不再有太大的变化。


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-5.cc"
*/
