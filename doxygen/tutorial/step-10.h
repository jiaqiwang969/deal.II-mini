/**
@page step_10 The step-10 tutorial program
This tutorial depends on step-7.

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
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
examples/step-10/doc/intro.dox

<a name="Intro"></a>

<a name="Introduction"></a><h1>Introduction</h1>



这是一个相当短的例子，只显示了使用高阶映射的一些方面。我们所说的 <em> 映射 </em> 是指单元格（即单位线、正方形或立方体）与现实空间中的单元格之间的转换。在前面所有的例子中，我们都隐含地使用了线性或d-线性映射；你根本不会注意到这一点，因为如果你不做任何特别的事情，这就是发生的情况。然而，如果你的域有弯曲的边界，在有些情况下，边界的片状线性逼近（即通过直线段）是不够的，你希望你的计算域也是一个使用弯曲边界的真实域的逼近。如果边界近似使用分片的二次抛物线来近似真实边界，那么我们说这是二次或 $Q_2$ 的近似。如果我们使用成片的立体多项式的图形，那么这是一个 $Q_3$ 的近似，以此类推。




对于某些微分方程，已知如果精确域的边界是弯曲的，那么边界的片状线性近似，即 $Q_1$ 映射，是不够的。例如，使用 $C^1$ 元素的偏谐方程，或具有弯曲反射边界的域上的气体动力学的欧拉方程。在这些情况下，有必要使用高阶映射来计算积分。如果我们不使用这样的高阶映射，那么边界的逼近顺序将主导整个数值方案的收敛顺序，而不考虑域的内部离散化的收敛顺序。




我们没有用这些更复杂的例子来证明高阶映射的使用，而是只做了一个简单的计算：用两种不同的方法计算 $\pi=3.141592653589793238462643\ldots$ 的值。




第一种方法使用单位半径的圆的三角形近似，并在其上积分一个单位幅度的常数函数（ $f = 1$ ）。当然，如果领域是精确的单位圆，那么面积将是 $\pi$ ，但由于我们只使用片状多项式段的近似，我们积分的面积值并不完全是 $\pi$  。然而，众所周知，当我们细化三角形时， $Q_p$ 映射以 $h^{p+1}$ 的阶数逼近边界，其中 $h$ 是网格大小。我们将检查圆的计算面积值，以及它们在不同映射的网格细化下向 $\pi$ 收敛的情况。我们还将发现一种收敛行为，这种行为一开始是令人惊讶的，但有一个很好的解释。




第二种方法与此类似，但这次不使用三角形单位圆的面积，而是使用其周长。   $\pi$ 然后用周长的一半来近似，因为我们选择半径等于1。




 @note  本教程实质上展示了如何为积分选择一个特定的映射，方法是将一个特定的几何体附加到三角形上（例如在步骤1中已经完成），然后将一个映射参数传递给FEValues类，该类用于deal.II中的所有积分。我们选择的几何体是一个圆，deal.II已经有一个可以使用的类（SphericalManifold）。如果你想定义你自己的几何体，例如，因为它很复杂，不能用deal.II中已有的类来描述，你会想通过步骤53来阅读。


 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 以下第一个include文件现在可能已经众所周知，不需要进一步解释。
 * 

 * 
 * 
 * @code
 * #include <deal.II/base/quadrature_lib.h> 
 * #include <deal.II/base/convergence_table.h> 
 * #include <deal.II/grid/grid_generator.h> 
 * #include <deal.II/grid/manifold_lib.h> 
 * #include <deal.II/grid/tria.h> 
 * #include <deal.II/grid/grid_out.h> 
 * #include <deal.II/dofs/dof_handler.h> 
 * #include <deal.II/fe/fe_values.h> 
 * 
 * @endcode
 * 
 * 这个包含文件是新的。即使我们在本教程中不求解PDE，我们也要使用FE_Nothing类提供的自由度为零的假有限元。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/fe_nothing.h> 
 * 
 * @endcode
 * 
 * 下面的头文件也是新的：在其中，我们声明了MappingQ类，我们将使用该类来处理任意阶的多项式映射。
 * 

 * 
 * 
 * @code
 * #include <deal.II/fe/mapping_q.h> 
 * 
 * @endcode
 * 
 * 这又是一个C++的文件。
 * 

 * 
 * 
 * @code
 * #include <iostream> 
 * #include <fstream> 
 * #include <cmath> 
 * 
 * @endcode
 * 
 * 最后一步和以前的程序一样。
 * 

 * 
 * 
 * @code
 * namespace Step10 
 * { 
 *   using namespace dealii; 
 * 
 * @endcode
 * 
 * 现在，由于我们要计算 $\pi$ 的值，我们必须与一些东西进行比较。这些是 $\pi$ 的前几个数字，我们事先定义好，以便以后使用。由于我们想计算两个数字的差值，而这两个数字是相当精确的，计算出的 $\pi$ 的近似值的精度在一个双数变量可以容纳的数字范围内，所以我们宁可将参考值声明为 <code>long double</code> ，并给它增加一些数字。
 * 

 * 
 * 
 * @code
 *   const long double pi = 3.141592653589793238462643L; 
 * 
 * @endcode
 * 
 * 然后，第一个任务将是生成一些输出。由于这个程序非常小，我们在其中没有采用面向对象的技术，也没有声明类（当然，我们使用了库的面向对象的功能）。相反，我们只是将功能打包成独立的函数。我们使这些函数成为空间维数的模板，以符合使用deal.II时的通常做法，尽管我们只对两个空间维数使用这些函数，当试图对任何其他空间维数使用时，会出现异常。
 * 

 * 
 * 这些函数中的第一个只是生成一个圆的三角形（hyperball），并输出 $Q_p$ 的不同值的单元的映射。然后，我们细化一次网格，再做一次。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void gnuplot_output() 
 *   { 
 *     std::cout << "Output of grids into gnuplot files:" << std::endl 
 *               << "===================================" << std::endl; 
 * @endcode
 * 
 * 因此，
 * 首先生成一个圆的粗略三角剖分，并将一个合适的边界描述与之关联。默认情况下， GridGenerator::hyper_ball 将SphericalManifold附加到边界上（内部使用FlatManifold），所以我们简单地调用该函数并继续前进。
 * 

 * 
 * 
 * @code
 *     Triangulation<dim> triangulation; 
 *     GridGenerator::hyper_ball(triangulation); 
 * 
 * @endcode
 * 
 * 然后在当前网格上交替生成 $Q_1$ 、 $Q_2$ 和 $Q_3$ 映射的输出，以及（在循环体的末端）对网格进行一次全局细化。
 * 

 * 
 * 
 * @code
 *     for (unsigned int refinement = 0; refinement < 2; ++refinement) 
 *       { 
 *         std::cout << "Refinement level: " << refinement << std::endl; 
 * 
 *         std::string filename_base = "ball_" + std::to_string(refinement); 
 * 
 *         for (unsigned int degree = 1; degree < 4; ++degree) 
 *           { 
 *             std::cout << "Degree = " << degree << std::endl; 
 * 
 * @endcode
 * 
 * 为此，首先建立一个描述映射的对象。这是用MappingQ类来完成的，该类在构造函数中采用了它应使用的多项式程度作为参数。
 * 

 * 
 * 
 * @code
 *             const MappingQ<dim> mapping(degree); 
 * 
 * @endcode
 * 
 * 顺便提一下，对于一个片状线性映射，你可以给MappingQ的构造函数一个 <code>1</code> 的值，但也有一个MappingQ1类可以达到同样的效果。历史上，它以比MappingQ更简单的方式做了很多事情，但今天只是后者的一个包装。然而，如果你没有明确指定另一个映射，它仍然是库中许多地方隐含使用的类。
 * 

 * 
 * 为了真正用这个映射写出现在的网格，我们设置了一个对象，我们将用它来输出。我们将生成Gnuplot输出，它由一组描述映射的三角图的线条组成。默认情况下，三角剖分的每个面只画一条线，但由于我们想明确地看到映射的效果，所以我们想更详细地了解这些面。这可以通过传递给输出对象一个包含一些标志的结构来实现。在目前的情况下，由于Gnuplot只能画直线，我们在面孔上输出了一些额外的点，这样每个面孔就由30条小线来画，而不是只有一条。这足以让我们看到一条弯曲的线，而不是一组直线的印象。
 * 

 * 
 * 
 * @code
 *             GridOut               grid_out; 
 *             GridOutFlags::Gnuplot gnuplot_flags(false, 60); 
 *             grid_out.set_flags(gnuplot_flags); 
 * 
 * @endcode
 * 
 * 最后，生成一个文件名和一个用于输出的文件。
 * 

 * 
 * 
 * @code
 *             std::string filename = 
 *               filename_base + "_mapping_q_" + std::to_string(degree) + ".dat"; 
 *             std::ofstream gnuplot_file(filename); 
 * 
 * @endcode
 * 
 * 然后把三角图写到这个文件里。该函数的最后一个参数是一个指向映射对象的指针。这个参数有一个默认值，如果没有给出值，就会取一个简单的MappingQ1对象，我们在上面简单介绍过。这样就会在输出中产生一个真实边界的片状线性近似。
 * 

 * 
 * 
 * @code
 *             grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping); 
 *           } 
 *         std::cout << std::endl; 
 * 
 * @endcode
 * 
 * 在循环结束时，对网格进行全局细化。
 * 

 * 
 * 
 * @code
 *         triangulation.refine_global(); 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 现在我们进行代码的主要部分，即 $\pi$ 的近似。圆的面积当然是由 $\pi r^2$ 给出的，所以有一个半径为1的圆，面积代表的只是被搜索的数字。面积的数值计算是通过在整个计算域中积分值为1的常数函数来进行的，即通过计算面积 $\int_K 1 dx=\int_{\hat K} 1
 * \ \textrm{det}\ J(\hat x) d\hat x \approx \sum_i \textrm{det}
 * \ J(\hat x_i)w(\hat x_i)$ ，其中总和延伸到三角形中所有活动单元上的所有正交点， $w(x_i)$ 是正交点的重量 $x_i$ 。每个单元上的积分都是通过数字正交来逼近的，因此我们唯一需要的额外成分是建立一个FEValues对象，提供每个单元的相应`JxW`值。注意`JxW`是指<i>Jacobian determinant
 * times weight</i>的缩写；因为在数字正交中，两个因子总是出现在相同的地方，所以我们只提供合并的数量，而不是两个单独的数量）。我们注意到，在这里我们不会在其最初的目的中使用FEValues对象，即用于计算特定正交点上的特定有限元的基函数值。相反，我们只用它来获得正交点的 "JxW"，而不考虑我们将给FEValues对象的构造者的（假）有限元。给予FEValues对象的实际有限元根本不使用，所以我们可以给任何。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void compute_pi_by_area() 
 *   { 
 *     std::cout << "Computation of Pi by the area:" << std::endl 
 *               << "==============================" << std::endl; 
 * 
 * @endcode
 * 
 * 对于所有单元的数字正交，我们采用足够高的正交规则。我们选择8阶的QGauss（4点），以确保数字正交引起的误差比由于边界近似的阶数，即所采用的映射的阶数（最大6）要高。请注意，积分，雅各布行列式，不是一个多项式函数（相反，它是一个有理函数），所以我们不使用高斯正交来获得积分的精确值，就像在有限元计算中经常做的那样，但也可以使用任何类似阶数的正交公式来代替。
 * 

 * 
 * 
 * @code
 *     const QGauss<dim> quadrature(4); 
 * 
 * @endcode
 * 
 * 现在开始在多项式映射度=1...4的基础上进行循环。
 * 

 * 
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree) 
 *       { 
 *         std::cout << "Degree = " << degree << std::endl; 
 * 
 * @endcode
 * 
 * 首先生成三角形、边界和映射对象，正如已经看到的那样。
 * 

 * 
 * 
 * @code
 *         Triangulation<dim> triangulation; 
 *         GridGenerator::hyper_ball(triangulation); 
 * 
 *         const MappingQ<dim> mapping(degree); 
 * 
 * @endcode
 * 
 * 我们现在创建一个有限元。与其他的例子程序不同，我们实际上不需要用形状函数做任何计算；我们只需要FEValues对象的`JxW`值。因此，我们使用特殊的有限元类FE_Nothing，它的每个单元的自由度正好为零（顾名思义，每个单元的局部基础为空集）。FE_Nothing的一个比较典型的用法见  step-46  。
 * 

 * 
 * 
 * @code
 *         const FE_Nothing<dim> fe; 
 * 
 * @endcode
 * 
 * 同样地，我们需要创建一个DoFHandler对象。我们实际上并没有使用它，但是它将为我们提供`active_cell_iterators'，这是重新初始化三角形的每个单元上的FEValues对象所需要的。
 * 

 * 
 * 
 * @code
 *         DoFHandler<dim> dof_handler(triangulation); 
 * 
 * @endcode
 * 
 * 现在我们设置FEValues对象，向构造函数提供Mapping、假有限元和正交对象，以及要求只在正交点提供`JxW`值的更新标志。这告诉FEValues对象在调用 <code>reinit</code> 函数时不需要计算其他数量，从而节省计算时间。
 * 

 * 
 * 与之前的例子程序相比，FEValues对象的构造最重要的区别是，我们传递了一个映射对象作为第一个参数，它将被用于计算从单元到实数单元的映射。在以前的例子中，这个参数被省略了，结果是隐含地使用了MappingQ1类型的对象。
 * 

 * 
 * 
 * @code
 *         FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values); 
 * 
 * @endcode
 * 
 * 我们使用一个ConvergenceTable类的对象来存储所有重要的数据，如 $\pi$ 的近似值和与 $\pi$ 的真实值相比的误差。我们还将使用ConvergenceTable类提供的函数来计算 $\pi$ 的近似值的收敛率。
 * 

 * 
 * 
 * @code
 *         ConvergenceTable table; 
 * 
 * @endcode
 * 
 * 现在我们在三角形的几个细化步骤上循环。
 * 

 * 
 * 
 * @code
 *         for (unsigned int refinement = 0; refinement < 6; 
 *              ++refinement, triangulation.refine_global(1)) 
 *           { 
 * 
 * @endcode
 * 
 * 在这个循环中，我们首先将当前三角形的活动单元的数量添加到表格中。这个函数会自动创建一个上标为 "cells "的表格列，以防这个列之前没有被创建。
 * 

 * 
 * 
 * @code
 *             table.add_value("cells", triangulation.n_active_cells()); 
 * 
 * @endcode
 * 
 * 然后我们为虚拟有限元分配自由度。严格来说，在我们的特殊情况下，我们不需要这个函数的调用，但我们调用它是为了让DoFHandler高兴 -- 否则它将在下面的 FEValues::reinit 函数中抛出一个断言。
 * 

 * 
 * 
 * @code
 *             dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 我们将变量面积定义为 "长双"，就像我们之前为 "pi "变量所做的那样。
 * 

 * 
 * 
 * @code
 *             long double area = 0; 
 * 
 * @endcode
 * 
 * 现在我们循环所有的单元格，重新初始化每个单元格的FEValues对象，并将该单元格的所有`JxW`值加到`area`上......
 * 

 * 
 * 
 * @code
 *             for (const auto &cell : dof_handler.active_cell_iterators()) 
 *               { 
 *                 fe_values.reinit(cell); 
 *                 for (unsigned int i = 0; i < fe_values.n_quadrature_points; ++i) 
 *                   area += static_cast<long double>(fe_values.JxW(i)); 
 *               } 
 * 
 * @endcode
 * 
 * ...并将得到的区域值和错误存储在表中。我们需要静态转换为双数，因为没有实现add_value(string, long double)函数。请注意，这也涉及到第二个调用，因为 <code>std</code> 命名空间中的 <code>fabs</code> 函数在其参数类型上是重载的，所以存在一个获取并返回 <code>long double</code> 的版本，而全局命名空间中只有一个这样的函数被声明（获取并返回一个双数）。
 * 

 * 
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(area)); 
 *             table.add_value("error", static_cast<double>(std::fabs(area - pi))); 
 *           } 
 * 
 * @endcode
 * 
 * 我们想计算`error`列的收敛率。因此我们需要在调用`evaluate_all_convergence_rates`之前，将其他列从收敛率评估中省略。
 * 

 * 
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells"); 
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
 *         table.evaluate_all_convergence_rates( 
 *                                     ConvergenceTable::reduction_rate_log2); 
 * 
 * @endcode
 * 
 * 最后我们设置一些量的输出精度和科学模式...
 * 

 * 
 * 
 * @code
 *         table.set_precision("eval.pi", 16); 
 *         table.set_scientific("error", true); 
 * 
 * @endcode
 * 
 * ...并将整个表格写到  std::cout.  。
 * 
 * @code
 *         table.write_text(std::cout); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * 
 * @endcode
 * 
 * 下面的第二个函数也是计算 $\pi$ 的近似值，但这次是通过域的周长 $2\pi r$ 而不是面积。这个函数只是前一个函数的一个变体。因此，我们主要是给出不同之处的文件。
 * 

 * 
 * 
 * @code
 *   template <int dim> 
 *   void compute_pi_by_perimeter() 
 *   { 
 *     std::cout << "Computation of Pi by the perimeter:" << std::endl 
 *               << "===================================" << std::endl; 
 * 
 * @endcode
 * 
 * 我们采取同样的正交顺序，但这次是`dim-1`维正交，因为我们将在（边界）线上而不是在单元上积分。
 * 

 * 
 * 
 * @code
 *     const QGauss<dim - 1> quadrature(4); 
 * 
 * @endcode
 * 
 * 我们在所有度数上循环，创建三角形、边界、映射、假有限元和DoFHandler对象，如之前所见。
 * 

 * 
 * 
 * @code
 *     for (unsigned int degree = 1; degree < 5; ++degree) 
 *       { 
 *         std::cout << "Degree = " << degree << std::endl; 
 *         Triangulation<dim> triangulation; 
 *         GridGenerator::hyper_ball(triangulation); 
 * 
 *         const MappingQ<dim>   mapping(degree); 
 *         const FE_Nothing<dim> fe; 
 * 
 *         DoFHandler<dim> dof_handler(triangulation); 
 * 
 * @endcode
 * 
 * 然后我们创建一个FEFaceValues对象，而不是像前一个函数中的FEValues对象。同样，我们传递一个映射作为第一个参数。
 * 

 * 
 * 
 * @code
 *         FEFaceValues<dim> fe_face_values(mapping, 
 *                                          fe, 
 *                                          quadrature, 
 *                                          update_JxW_values); 
 *         ConvergenceTable  table; 
 * 
 *         for (unsigned int refinement = 0; refinement < 6; 
 *              ++refinement, triangulation.refine_global(1)) 
 *           { 
 *             table.add_value("cells", triangulation.n_active_cells()); 
 * 
 *             dof_handler.distribute_dofs(fe); 
 * 
 * @endcode
 * 
 * 现在我们在所有单元和每个单元的所有面上运行。只有边界面上的`JxW`值的贡献被添加到长双变量`周长`中。
 * 

 * 
 * 
 * @code
 *             long double perimeter = 0; 
 *             for (const auto &cell : dof_handler.active_cell_iterators()) 
 *               for (const auto &face : cell->face_iterators()) 
 *                 if (face->at_boundary()) 
 *                   { 
 * 
 * @endcode
 * 
 * 我们用单元格迭代器和面的编号重新启动FEFaceValues对象。
 * 

 * 
 * 
 * @code
 *                     fe_face_values.reinit(cell, face); 
 *                     for (unsigned int i = 0; 
 *                          i < fe_face_values.n_quadrature_points; 
 *                          ++i) 
 *                       perimeter += 
 *                         static_cast<long double>(fe_face_values.JxW(i)); 
 *                   } 
 * 
 * @endcode
 * 
 * 然后将评估后的数值存储在表中...
 * 

 * 
 * 
 * @code
 *             table.add_value("eval.pi", static_cast<double>(perimeter / 2.0L)); 
 *             table.add_value( 
 *               "error", static_cast<double>(std::fabs(perimeter / 2.0L - pi))); 
 *           } 
 * 
 * @endcode
 * 
 * ......然后像前一个函数那样结束这个函数。
 * 

 * 
 * 
 * @code
 *         table.omit_column_from_convergence_rate_evaluation("cells"); 
 *         table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
 *         table.evaluate_all_convergence_rates( 
 *           ConvergenceTable::reduction_rate_log2); 
 * 
 *         table.set_precision("eval.pi", 16); 
 *         table.set_scientific("error", true); 
 * 
 *         table.write_text(std::cout); 
 * 
 *         std::cout << std::endl; 
 *       } 
 *   } 
 * } // namespace Step10 
 * 
 * @endcode
 * 
 * 下面的主函数只是按照上述函数的出现顺序来调用它们。除此以外，它看起来就像以前的教程程序的主函数一样。
 * 

 * 
 * 
 * @code
 * int main() 
 * { 
 *   try 
 *     { 
 *       std::cout.precision(16); 
 * 
 *       const unsigned int dim = 2; 
 * 
 *       Step10::gnuplot_output<dim>(); 
 * 
 *       Step10::compute_pi_by_area<dim>(); 
 *       Step10::compute_pi_by_perimeter<dim>(); 
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
examples/step-10/doc/results.dox



<a name="Results"></a><h1>Results</h1>



该程序执行两个任务，第一个是生成映射域的可视化，第二个是通过所述的两种方法计算π。让我们先看一下生成的图形。它们是以Gnuplot格式生成的，可以用以下命令查看

@code
set style data lines
set size ratio -1
unset key
unset xtics
unset ytics
plot [-1:1][-1:1] "ball_0_mapping_q_1.dat" lw 4 lt rgb "black"
@endcode

或使用其他文件名之一。第二行确保生成的输出的长宽比实际上是1:1，也就是说，一个圆在你的屏幕上被画成一个圆，而不是一个椭圆。第三行关闭了图形中的键，因为那只会打印信息（文件名），而这些信息现在并不那么重要。同样地，第四行和第五行关闭了刻度线。然后生成具有特定线宽（"lw"，这里设置为4）和线型（"lt"，这里选择线应使用RGB颜色 "黑色"）的图。

下表显示了 $Q_1$ 、 $Q_2$ 和 $Q_3$ 映射的三角计算域，为原始粗网格（左）和一次均匀细化网格（右）。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q1.svg" alt="磁盘的五单元离散化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q1.svg" alt="磁盘的20单元离散化（即。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q2.svg" alt="五格离散化的圆盘，边缘为二次方。边界与实际的圆几乎没有区别。" width="400" height="400" > </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q2.svg" alt="带有二次方边缘的20个单元离散化。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_0_q3.svg" alt="带有三次方边缘的圆的五单元离散化。边界与实际的圆几乎没有区别。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_ball_1_q3.svg" alt="具有立方体边缘的20单元离散化。" width="400" height="400"> </div> </div>

这些图片显示了高阶映射的明显优势：它们在相当粗的网格上也能相当好地接近真实边界。为了进一步证明这一点，这里是使用 $Q_2$ 和 $Q_3$ 映射的粗网格的右上角四分之一圈的一部分，其中红色虚线标志着实际的圆。

<div class="twocolumn" style="width: 80%"> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q2.svg" alt="二次离散化的特写。二次插值和实际圆之间的距离很小。" width="400" height="400"> </div> <div> <img src="https://www.dealii.org/images/steps/developer/step_10_exact_vs_interpolate_q3.svg" alt="立方离散化的特写。立体插值和实际圆之间的距离非常小。" width="400" height="400"> </div> </div> </div>

很明显，二次映射很好地接近了边界，而对于三次映射来说，对于粗略的网格来说，近似域和真实域之间的差异已经很难看到了。你还可以看到，映射只在三角形的外部边界有一些变化。在内部，所有的线条仍然是由线性函数表示的，这就导致了只在边界的单元上进行额外的计算。因此，高阶映射通常不会比低阶映射明显地慢，因为额外的计算只在所有单元格的一个小子集上执行。




该程序的第二个目的是计算π的值，以达到良好的准确性。这是程序的这一部分的输出。

@code
Output of grids into gnuplot files:
===================================
Refinement level: 0
Degree = 1
Degree = 2
Degree = 3


Refinement level: 1
Degree = 1
Degree = 2
Degree = 3


Computation of Pi by the area:
==============================
Degree = 1
cells      eval.pi            error
    5 1.9999999999999993 1.1416e+00    -
   20 2.8284271247461890 3.1317e-01 1.87
   80 3.0614674589207174 8.0125e-02 1.97
  320 3.1214451522580511 2.0148e-02 1.99
 1280 3.1365484905459380 5.0442e-03 2.00
 5120 3.1403311569547516 1.2615e-03 2.00


Degree = 2
cells      eval.pi            error
    5 3.1045694996615860 3.7023e-02    -
   20 3.1391475703122267 2.4451e-03 3.92
   80 3.1414377167038290 1.5494e-04 3.98
  320 3.1415829366419006 9.7169e-06 4.00
 1280 3.1415920457576898 6.0783e-07 4.00
 5120 3.1415926155921117 3.7998e-08 4.00


Degree = 3
cells      eval.pi            error
    5 3.1410033851499288 5.8927e-04    -
   20 3.1415830393583839 9.6142e-06 5.94
   80 3.1415925017363797 1.5185e-07 5.98
  320 3.1415926512106696 2.3791e-09 6.00
 1280 3.1415926535525927 3.7200e-11 6.00
 5120 3.1415926535892100 5.8302e-13 6.00


Degree = 4
cells      eval.pi            error
    5 3.1415871927401131 5.4608e-06    -
   20 3.1415926314742428 2.2116e-08 7.95
   80 3.1415926535026202 8.7173e-11 7.99
  320 3.1415926535894498 3.4350e-13 7.99
 1280 3.1415926535897896 3.4671e-15 6.63
 5120 3.1415926535897909 2.4009e-15 0.53


Computation of Pi by the perimeter:
===================================
Degree = 1
cells      eval.pi            error
    5 2.8284271247461898 3.1317e-01    -
   20 3.0614674589207178 8.0125e-02 1.97
   80 3.1214451522580520 2.0148e-02 1.99
  320 3.1365484905459389 5.0442e-03 2.00
 1280 3.1403311569547525 1.2615e-03 2.00
 5120 3.1412772509327724 3.1540e-04 2.00


Degree = 2
cells      eval.pi            error
    5 3.1248930668550594 1.6700e-02    -
   20 3.1404050605605449 1.1876e-03 3.81
   80 3.1415157631807009 7.6890e-05 3.95
  320 3.1415878042798613 4.8493e-06 3.99
 1280 3.1415923498174534 3.0377e-07 4.00
 5120 3.1415926345931995 1.8997e-08 4.00


Degree = 3
cells      eval.pi            error
    5 3.1414940401456048 9.8613e-05    -
   20 3.1415913432549156 1.3103e-06 6.23
   80 3.1415926341726910 1.9417e-08 6.08
  320 3.1415926532906897 2.9910e-10 6.02
 1280 3.1415926535851355 4.6578e-12 6.00
 5120 3.1415926535897190 7.4216e-14 5.97


Degree = 4
cells      eval.pi            error
    5 3.1415921029432572 5.5065e-07     -
   20 3.1415926513737595 2.2160e-09  7.96
   80 3.1415926535810712 8.7222e-12  7.99
  320 3.1415926535897576 3.5525e-14  7.94
 1280 3.1415926535897936 4.6729e-16  6.25
 5120 3.1415926535897918 1.4929e-15 -1.68
@endcode






从输出的一个直接观察结果是，在所有情况下，数值都迅速收敛到 $\pi=3.141592653589793238462643$ 的真实值。请注意，对于 $Q_4$ 的映射，我们已经进入了四舍五入误差的制度，收敛率趋于平缓，这已经是相当多的了。然而，也请注意，对于 $Q_1$ 映射，即使在最细的网格上，精度也明显比 $Q_3$ 映射的粗网格上要差得多!




输出的最后一列显示了收敛顺序，以网格宽度的幂为单位  $h$  。在介绍中，我们曾说过  $Q_p$  映射的收敛顺序应该是  $h^{p+1}$  。然而，在所示的例子中，顺序是 $h^{2p}$  !这个起初令人惊讶的事实可以用 $Q_p$ 映射的特性来解释。在<i>p</i>阶时，它使用的支持点是基于<i>p</i>+1点的Gauss-Lobatto正交规则，以这样的方式选择支持点，使正交规则在2<i>p</i>阶时收敛。尽管这些点在这里只用于插值<i>p</i>阶多项式，但我们在数值评估积分时得到了超收敛效应，导致观察到的高阶收敛。这种效应在下面的出版物中也有详细讨论。A. Bonito, A. Demlow, and J. Owen:"拉普拉斯-贝特拉米算子的特征值和特征函数的有限元近似的先验误差估计"，提交，2018年）。)


 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-10.cc"
*/
