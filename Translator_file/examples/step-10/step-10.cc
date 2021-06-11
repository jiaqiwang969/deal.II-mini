CCTest_file/step-10.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2001 - 2021 by the deal.II authors 
 * 
 * This file is part of the deal.II library. 
 * 
 * The deal.II library is free software; you can use it, redistribute 
 * it, and/or modify it under the terms of the GNU Lesser General 
 * Public License as published by the Free Software Foundation; either 
 * version 2.1 of the License, or (at your option) any later version. 
 * The full text of the license can be found in the file LICENSE.md at 
 * the top level directory of deal.II. 
 * 
 * --------------------------------------------------------------------- 
 * 
 * Authors: Wolfgang Bangerth, Ralf Hartmann, University of Heidelberg, 2001 
 */ 



// 以下第一个include文件现在可能已经众所周知，不需要进一步解释。

#include <deal.II/base/quadrature_lib.h> 
#include <deal.II/base/convergence_table.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/dofs/dof_handler.h> 
#include <deal.II/fe/fe_values.h> 

// 这个包含文件是新的。即使我们在本教程中不求解PDE，我们也要使用FE_Nothing类提供的自由度为零的假有限元。

#include <deal.II/fe/fe_nothing.h> 

// 下面的头文件也是新的：在其中，我们声明了MappingQ类，我们将使用该类来处理任意阶的多项式映射。

#include <deal.II/fe/mapping_q.h> 

// 这又是一个C++的文件。

#include <iostream> 
#include <fstream> 
#include <cmath> 

// 最后一步和以前的程序一样。

namespace Step10 
{ 
  using namespace dealii; 

// 现在，由于我们要计算 $\pi$ 的值，我们必须与一些东西进行比较。这些是 $\pi$ 的前几个数字，我们事先定义好，以便以后使用。由于我们想计算两个数字的差值，而这两个数字是相当精确的，计算出的 $\pi$ 的近似值的精度在一个双数变量可以容纳的数字范围内，所以我们宁可将参考值声明为 <code>long double</code> ，并给它增加一些数字。

  const long double pi = 3.141592653589793238462643L; 

// 然后，第一个任务将是生成一些输出。由于这个程序非常小，我们在其中没有采用面向对象的技术，也没有声明类（当然，我们使用了库的面向对象的功能）。相反，我们只是将功能打包成独立的函数。我们使这些函数成为空间维数的模板，以符合使用deal.II时的通常做法，尽管我们只对两个空间维数使用这些函数，当试图对任何其他空间维数使用时，会出现异常。

// 这些函数中的第一个只是生成一个圆的三角形（hyperball），并输出 $Q_p$ 的不同值的单元的映射。然后，我们细化一次网格，再做一次。

  template <int dim> 
  void gnuplot_output() 
  { 
    std::cout << "Output of grids into gnuplot files:" << std::endl 
              << "===================================" << std::endl; 
//因此，
//首先生成一个圆的粗略三角剖分，并将一个合适的边界描述与之关联。默认情况下， GridGenerator::hyper_ball 将SphericalManifold附加到边界上（内部使用FlatManifold），所以我们简单地调用该函数并继续前进。

    Triangulation<dim> triangulation; 
    GridGenerator::hyper_ball(triangulation); 

// 然后在当前网格上交替生成 $Q_1$ 、 $Q_2$ 和 $Q_3$ 映射的输出，以及（在循环体的末端）对网格进行一次全局细化。

    for (unsigned int refinement = 0; refinement < 2; ++refinement) 
      { 
        std::cout << "Refinement level: " << refinement << std::endl; 

        std::string filename_base = "ball_" + std::to_string(refinement); 

        for (unsigned int degree = 1; degree < 4; ++degree) 
          { 
            std::cout << "Degree = " << degree << std::endl; 

// 为此，首先建立一个描述映射的对象。这是用MappingQ类来完成的，该类在构造函数中采用了它应使用的多项式程度作为参数。

            const MappingQ<dim> mapping(degree); 

// 顺便提一下，对于一个片状线性映射，你可以给MappingQ的构造函数一个 <code>1</code> 的值，但也有一个MappingQ1类可以达到同样的效果。历史上，它以比MappingQ更简单的方式做了很多事情，但今天只是后者的一个包装。然而，如果你没有明确指定另一个映射，它仍然是库中许多地方隐含使用的类。

// 为了真正用这个映射写出现在的网格，我们设置了一个对象，我们将用它来输出。我们将生成Gnuplot输出，它由一组描述映射的三角图的线条组成。默认情况下，三角剖分的每个面只画一条线，但由于我们想明确地看到映射的效果，所以我们想更详细地了解这些面。这可以通过传递给输出对象一个包含一些标志的结构来实现。在目前的情况下，由于Gnuplot只能画直线，我们在面孔上输出了一些额外的点，这样每个面孔就由30条小线来画，而不是只有一条。这足以让我们看到一条弯曲的线，而不是一组直线的印象。

            GridOut               grid_out; 
            GridOutFlags::Gnuplot gnuplot_flags(false, 60); 
            grid_out.set_flags(gnuplot_flags); 

// 最后，生成一个文件名和一个用于输出的文件。

            std::string filename = 
              filename_base + "_mapping_q_" + std::to_string(degree) + ".dat"; 
            std::ofstream gnuplot_file(filename); 

// 然后把三角图写到这个文件里。该函数的最后一个参数是一个指向映射对象的指针。这个参数有一个默认值，如果没有给出值，就会取一个简单的MappingQ1对象，我们在上面简单介绍过。这样就会在输出中产生一个真实边界的片状线性近似。

            grid_out.write_gnuplot(triangulation, gnuplot_file, &mapping); 
          } 
        std::cout << std::endl; 

// 在循环结束时，对网格进行全局细化。

        triangulation.refine_global(); 
      } 
  } 

// 现在我们进行代码的主要部分，即 $\pi$ 的近似。圆的面积当然是由 $\pi r^2$ 给出的，所以有一个半径为1的圆，面积代表的只是被搜索的数字。面积的数值计算是通过在整个计算域中积分值为1的常数函数来进行的，即通过计算面积 $\int_K 1 dx=\int_{\hat K} 1
// \ \textrm{det}\ J(\hat x) d\hat x \approx \sum_i \textrm{det}
// \ J(\hat x_i)w(\hat x_i)$ ，其中总和延伸到三角形中所有活动单元上的所有正交点， $w(x_i)$ 是正交点的重量 $x_i$ 。每个单元上的积分都是通过数字正交来逼近的，因此我们唯一需要的额外成分是建立一个FEValues对象，提供每个单元的相应`JxW`值。注意`JxW`是指<i>Jacobian determinant
// times weight</i>的缩写；因为在数字正交中，两个因子总是出现在相同的地方，所以我们只提供合并的数量，而不是两个单独的数量）。我们注意到，在这里我们不会在其最初的目的中使用FEValues对象，即用于计算特定正交点上的特定有限元的基函数值。相反，我们只用它来获得正交点的 "JxW"，而不考虑我们将给FEValues对象的构造者的（假）有限元。给予FEValues对象的实际有限元根本不使用，所以我们可以给任何。

  template <int dim> 
  void compute_pi_by_area() 
  { 
    std::cout << "Computation of Pi by the area:" << std::endl 
              << "==============================" << std::endl; 

// 对于所有单元的数字正交，我们采用足够高的正交规则。我们选择8阶的QGauss（4点），以确保数字正交引起的误差比由于边界近似的阶数，即所采用的映射的阶数（最大6）要高。请注意，积分，雅各布行列式，不是一个多项式函数（相反，它是一个有理函数），所以我们不使用高斯正交来获得积分的精确值，就像在有限元计算中经常做的那样，但也可以使用任何类似阶数的正交公式来代替。

    const QGauss<dim> quadrature(4); 

// 现在开始在多项式映射度=1...4的基础上进行循环。

    for (unsigned int degree = 1; degree < 5; ++degree) 
      { 
        std::cout << "Degree = " << degree << std::endl; 

// 首先生成三角形、边界和映射对象，正如已经看到的那样。

        Triangulation<dim> triangulation; 
        GridGenerator::hyper_ball(triangulation); 

        const MappingQ<dim> mapping(degree); 

// 我们现在创建一个有限元。与其他的例子程序不同，我们实际上不需要用形状函数做任何计算；我们只需要FEValues对象的`JxW`值。因此，我们使用特殊的有限元类FE_Nothing，它的每个单元的自由度正好为零（顾名思义，每个单元的局部基础为空集）。FE_Nothing的一个比较典型的用法见  step-46  。

        const FE_Nothing<dim> fe; 

// 同样地，我们需要创建一个DoFHandler对象。我们实际上并没有使用它，但是它将为我们提供`active_cell_iterators'，这是重新初始化三角形的每个单元上的FEValues对象所需要的。

        DoFHandler<dim> dof_handler(triangulation); 

// 现在我们设置FEValues对象，向构造函数提供Mapping、假有限元和正交对象，以及要求只在正交点提供`JxW`值的更新标志。这告诉FEValues对象在调用 <code>reinit</code> 函数时不需要计算其他数量，从而节省计算时间。

// 与之前的例子程序相比，FEValues对象的构造最重要的区别是，我们传递了一个映射对象作为第一个参数，它将被用于计算从单元到实数单元的映射。在以前的例子中，这个参数被省略了，结果是隐含地使用了MappingQ1类型的对象。

        FEValues<dim> fe_values(mapping, fe, quadrature, update_JxW_values); 

// 我们使用一个ConvergenceTable类的对象来存储所有重要的数据，如 $\pi$ 的近似值和与 $\pi$ 的真实值相比的误差。我们还将使用ConvergenceTable类提供的函数来计算 $\pi$ 的近似值的收敛率。

        ConvergenceTable table; 

// 现在我们在三角形的几个细化步骤上循环。

        for (unsigned int refinement = 0; refinement < 6; 
             ++refinement, triangulation.refine_global(1)) 
          { 

// 在这个循环中，我们首先将当前三角形的活动单元的数量添加到表格中。这个函数会自动创建一个上标为 "cells "的表格列，以防这个列之前没有被创建。

            table.add_value("cells", triangulation.n_active_cells()); 

// 然后我们为虚拟有限元分配自由度。严格来说，在我们的特殊情况下，我们不需要这个函数的调用，但我们调用它是为了让DoFHandler高兴 -- 否则它将在下面的 FEValues::reinit 函数中抛出一个断言。

            dof_handler.distribute_dofs(fe); 

// 我们将变量面积定义为 "长双"，就像我们之前为 "pi "变量所做的那样。

            long double area = 0; 

// 现在我们循环所有的单元格，重新初始化每个单元格的FEValues对象，并将该单元格的所有`JxW`值加到`area`上......

            for (const auto &cell : dof_handler.active_cell_iterators()) 
              { 
                fe_values.reinit(cell); 
                for (unsigned int i = 0; i < fe_values.n_quadrature_points; ++i) 
                  area += static_cast<long double>(fe_values.JxW(i)); 
              } 

// ...并将得到的区域值和错误存储在表中。我们需要静态转换为双数，因为没有实现add_value(string, long double)函数。请注意，这也涉及到第二个调用，因为 <code>std</code> 命名空间中的 <code>fabs</code> 函数在其参数类型上是重载的，所以存在一个获取并返回 <code>long double</code> 的版本，而全局命名空间中只有一个这样的函数被声明（获取并返回一个双数）。

            table.add_value("eval.pi", static_cast<double>(area)); 
            table.add_value("error", static_cast<double>(std::fabs(area - pi))); 
          } 

// 我们想计算`error`列的收敛率。因此我们需要在调用`evaluate_all_convergence_rates`之前，将其他列从收敛率评估中省略。

        table.omit_column_from_convergence_rate_evaluation("cells"); 
        table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
        table.evaluate_all_convergence_rates( 
                                    ConvergenceTable::reduction_rate_log2); 

// 最后我们设置一些量的输出精度和科学模式...

        table.set_precision("eval.pi", 16); 
        table.set_scientific("error", true); 

// ...并将整个表格写到  std::cout.  。
        table.write_text(std::cout); 

        std::cout << std::endl; 
      } 
  } 

// 下面的第二个函数也是计算 $\pi$ 的近似值，但这次是通过域的周长 $2\pi r$ 而不是面积。这个函数只是前一个函数的一个变体。因此，我们主要是给出不同之处的文件。

  template <int dim> 
  void compute_pi_by_perimeter() 
  { 
    std::cout << "Computation of Pi by the perimeter:" << std::endl 
              << "===================================" << std::endl; 

// 我们采取同样的正交顺序，但这次是`dim-1`维正交，因为我们将在（边界）线上而不是在单元上积分。

    const QGauss<dim - 1> quadrature(4); 

// 我们在所有度数上循环，创建三角形、边界、映射、假有限元和DoFHandler对象，如之前所见。

    for (unsigned int degree = 1; degree < 5; ++degree) 
      { 
        std::cout << "Degree = " << degree << std::endl; 
        Triangulation<dim> triangulation; 
        GridGenerator::hyper_ball(triangulation); 

        const MappingQ<dim>   mapping(degree); 
        const FE_Nothing<dim> fe; 

        DoFHandler<dim> dof_handler(triangulation); 

// 然后我们创建一个FEFaceValues对象，而不是像前一个函数中的FEValues对象。同样，我们传递一个映射作为第一个参数。

        FEFaceValues<dim> fe_face_values(mapping, 
                                         fe, 
                                         quadrature, 
                                         update_JxW_values); 
        ConvergenceTable  table; 

        for (unsigned int refinement = 0; refinement < 6; 
             ++refinement, triangulation.refine_global(1)) 
          { 
            table.add_value("cells", triangulation.n_active_cells()); 

            dof_handler.distribute_dofs(fe); 

// 现在我们在所有单元和每个单元的所有面上运行。只有边界面上的`JxW`值的贡献被添加到长双变量`周长`中。

            long double perimeter = 0; 
            for (const auto &cell : dof_handler.active_cell_iterators()) 
              for (const auto &face : cell->face_iterators()) 
                if (face->at_boundary()) 
                  { 

// 我们用单元格迭代器和面的编号重新启动FEFaceValues对象。

                    fe_face_values.reinit(cell, face); 
                    for (unsigned int i = 0; 
                         i < fe_face_values.n_quadrature_points; 
                         ++i) 
                      perimeter += 
                        static_cast<long double>(fe_face_values.JxW(i)); 
                  } 

// 然后将评估后的数值存储在表中...

            table.add_value("eval.pi", static_cast<double>(perimeter / 2.0L)); 
            table.add_value( 
              "error", static_cast<double>(std::fabs(perimeter / 2.0L - pi))); 
          } 

// ......然后像前一个函数那样结束这个函数。

        table.omit_column_from_convergence_rate_evaluation("cells"); 
        table.omit_column_from_convergence_rate_evaluation("eval.pi"); 
        table.evaluate_all_convergence_rates( 
          ConvergenceTable::reduction_rate_log2); 

        table.set_precision("eval.pi", 16); 
        table.set_scientific("error", true); 

        table.write_text(std::cout); 

        std::cout << std::endl; 
      } 
  } 
} // namespace Step10 

// 下面的主函数只是按照上述函数的出现顺序来调用它们。除此以外，它看起来就像以前的教程程序的主函数一样。

int main() 
{ 
  try 
    { 
      std::cout.precision(16); 

      const unsigned int dim = 2; 

      Step10::gnuplot_output<dim>(); 

      Step10::compute_pi_by_area<dim>(); 
      Step10::compute_pi_by_perimeter<dim>(); 
    } 
  catch (std::exception &exc) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Exception on processing: " << std::endl 
                << exc.what() << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 

      return 1; 
    } 
  catch (...) 
    { 
      std::cerr << std::endl 
                << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      std::cerr << "Unknown exception!" << std::endl 
                << "Aborting!" << std::endl 
                << "----------------------------------------------------" 
                << std::endl; 
      return 1; 
    } 

  return 0; 
} 


