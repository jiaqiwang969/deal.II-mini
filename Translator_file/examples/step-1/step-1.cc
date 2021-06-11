

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 1999 - 2021 by the deal.II authors 
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
 */ 


// @sect3{Include files}  

// 库中最基本的类是Triangulation类，它在这里声明。

#include <deal.II/grid/tria.h> 

// 这里有一些生成标准网格的函数。

#include <deal.II/grid/grid_generator.h> 

// 输出各种图形格式的网格。

#include <deal.II/grid/grid_out.h> 

// 这对于C++输出来说是需要的。

#include <iostream> 
#include <fstream> 

// 这是对 `std::sqrt` 和 `std::fabs` 函数声明的说明。

#include <cmath> 

//导入deal.II的最后一步是这样的。所有deal.II的函数和类都在一个命名空间 <code>dealii</code> 中，以确保它们不会与你可能想和deal.II一起使用的其他库的符号发生冲突。我们可以在使用这些函数和类时，在每个名字前加上 <code>dealii::</code> 的前缀，但这很快就会变得繁琐和令人厌烦。相反，我们只是简单地导入整个deal.II的名字空间，以供一般使用。

using namespace dealii; 
// @sect3{Creating the first mesh}  

// 在下面的第一个函数中，我们简单地使用单位方格作为域，并从中产生一个全局细化网格。

void first_grid() 
{ 

// 首先要做的是为二维域的三角化定义一个对象。

  Triangulation<2> triangulation; 

// 在这里和下面的许多情况下，类名后面的字符串"<2>"表示这是一个在两个空间维度上工作的对象。同样，也有一些三角形类的版本是在一个（"<1>"）和三个（"<3>"）空间维度上工作的。这种工作方式是通过一些模板魔法实现的，我们将在后面的示例程序中详细研究；在那里，我们也将看到如何以一种基本独立于维度的方式编写程序。

// 接下来，我们要用一个正方形领域的单个单元来填充三角结构。三角形被细化了四次，总共得到 $4^4=256$ 个单元。

  GridGenerator::hyper_cube(triangulation); 
  triangulation.refine_global(4); 

// 现在我们要将网格的图形表示写到输出文件中。deal.II的GridOut类可以用多种不同的输出格式来实现；在这里，我们选择可扩展矢量图（SVG）格式，你可以用你选择的网络浏览器来进行可视化。

  std::ofstream out("grid-1.svg"); 
  GridOut       grid_out; 
  grid_out.write_svg(triangulation, out); 
  std::cout << "Grid written to grid-1.svg" << std::endl; 
} 

//  @sect3{Creating the second mesh}  

// 下面第二个函数中的网格略微复杂一些，因为我们使用了一个环形域，并对结果进行了一次全局细化。

void second_grid() 
{ 

// 我们再次开始定义一个二维域的三角化对象。

  Triangulation<2> triangulation; 

// 然后我们用一个环形域来填充它。环的中心应是点(1,0)，内半径和外半径应是0.5和1。圆周单元的数量可以由这个函数自动调整，但我们选择在最后一个参数中明确设置为10。

  const Point<2> center(1, 0); 
  const double   inner_radius = 0.5, outer_radius = 1.0; 
  GridGenerator::hyper_shell( 
    triangulation, center, inner_radius, outer_radius, 10); 

// 默认情况下，三角测量假定所有边界都是直线，所有单元都是双线性四边形或三线性六边形，并且它们是由粗略网格（我们刚刚创建的）的单元定义的。除非我们做一些特别的事情，否则当需要引入新的点时，域被假定为由粗网格的直线划定，而新的点将简单地位于周围的中间。然而，在这里，我们知道领域是弯曲的，我们想让三角法根据底层的几何形状来放置新的点。幸运的是，一些优秀的灵魂实现了一个描述球状域的对象，而环是球状域的一个部分；它只需要环的中心，并自动计算出如何指示三角计算在哪里放置新的点。这在deal.II中的工作方式是，你用一个通常被称为 "流形指标 "的数字来标记你想要弯曲的三角形部分，然后告诉三角形在所有有这个流形指标的地方使用一个特定的 "流形对象"。具体如何操作在此并不重要（你可以在 step-53 和 @ref manifold 中阅读）。GridGenerator中的函数在大多数情况下为我们处理这个问题：它们将正确的流形附加到一个域上，这样当三角形被细化时，新的单元就会被放置在正确的位置上。在目前的情况下， GridGenerator::hyper_shell 为所有的单元格附加了一个球形流形：这将导致单元格在球面坐标的计算下被细化（因此新的单元格的边缘要么是径向的，要么是位于原点周围的同心圆）。

// 默认情况下（即对于手工创建的三角图或未调用GridGenerator函数（如 GridGenerator::hyper_shell 或 GridGenerator::hyper_ball), ），三角图的所有单元格和面都将其manifold_id设置为 numbers::flat_manifold_id, ，如果您想要一个产生直线边缘的流形，这是默认的，但您可以为个别单元格和面改变这个数字。在这种情况下，因此与数字0相关的曲面流形将不适用于那些流形指标为非零的部分，但其他流形描述对象可以与这些非零指标相关联。如果没有流形描述与特定的流形指标相关联，则暗示产生直角边缘的流形。(流形指标是一个略微复杂的话题；如果你对这里到底发生了什么感到困惑，你可能想看看 @ref GlossManifoldIndicator "关于这个话题的词汇表条目")。既然 GridGenerator::hyper_shell 选择的默认值是合理的，我们就不去管它。

// 为了演示如何在所有单元格上写一个循环，我们将分五个步骤向域的内圈细化网格。

  for (unsigned int step = 0; step < 5; ++step) 
    { 

// 接下来，我们需要对三角形的活动单元进行循环。你可以把三角形看作一个单元格的集合。如果它是一个数组，你只需要得到一个指针，用操作符`++`从一个元素递增到下一个元素。三角形的单元不是作为一个简单的数组来存储的，但是<i>iterator</i>的概念将指针的工作方式概括为任意的对象集合（更多信息见<a href= "http:en.wikipedia.org/wiki/Iterator#C.2B.2B">wikipedia</a>）。通常情况下，C++中的任何容器类型都会返回一个迭代器，指向集合的开始，方法称为`begin'，而迭代器则指向集合结束后的1，方法称为`end'。我们可以用操作符`++it`来增加一个迭代器`it`，用`*it`来解除引用以获得底层数据，并通过比较`it != collection.end()`来检查我们是否完成。

// 第二个重要的部分是我们只需要活动单元。活动单元是那些没有被进一步细化的单元，也是唯一可以被标记为进一步细化的单元。deal.II提供了迭代器类别，允许我们在<i>all</i>单元（包括活动单元的父单元）或只在活动单元上迭代。因为我们要的是后者，所以我们需要调用方法 Triangulation::active_cell_iterators().  。

//把所有这些放在一起，我们可以用
// @code{.cpp}
//      for (auto it = triangulation.active_cell_iterators().begin();
//           it != triangulation.active_cell_iterators().end();
//           ++it)
//        {
//          auto cell = *it;
//  //Then a miracle occurs...
//        }
//  @endcode
//  在一个三角形的所有活动单元上循环。 在这个循环的初始化器中，我们使用了`auto`关键字作为迭代器`it`的类型。`auto`关键字意味着被声明的对象的类型将从上下文中推断出来。当实际的类型名称很长，甚至可能是多余的时候，这个关键字很有用。如果你不确定类型是什么，想查一下结果支持什么操作，你可以去看方法的文档  Triangulation::active_cell_iterators().  在这个例子中，`it`的类型是  `Triangulation::active_cell_iterator`.  

// 虽然`auto`关键字可以让我们不用输入长长的数据类型名称，但我们仍然要输入大量冗余的关于开始和结束迭代器以及如何递增的声明。与其这样，我们不如使用<a href="http:en.cppreference.com/w/cpp/language/range-for">range-based for loops</a>，它将上面显示的所有语法包成一个更短的形式。

      for (auto &cell : triangulation.active_cell_iterators()) 
        { 
// @note  关于deal.II中使用的迭代器类的更多信息，见 @ref Iterators ，关于基于范围的for循环和`auto`关键字的更多信息，见 @ref CPP11 。

// 接下来，我们在单元格的所有顶点上循环。为此，我们查询一个顶点索引的迭代器（在2D中，这是一个包含元素`{0,1,2,3}`的数组，但是由于`cell->vertex_indices()`知道单元格所处的维度，因此返回的数组在所有维度上都是正确的，这使得无论我们在2D还是3D中运行这段代码都是正确的，也就是说，它实现了 "维度无关的编程" - 我们将在  step-4  中讨论一个重要部分）。

          for (const auto v : cell->vertex_indices()) 
            { 

// 如果这个单元格位于内边界，那么它至少有一个顶点必须位于内环上，因此与中心的径向距离正好是0.5，达到浮点精度。所以我们计算这个距离，如果我们发现一个顶点具有这个属性，我们就标记这个单元，以便以后进行细化。然后我们也可以打破所有顶点的循环，转到下一个单元。        因为离中心的距离是以浮点数计算的，所以我们必须期望我们所计算的东西只能精确到[round-off](https:en.wikipedia.org/wiki/Round-off_error)以内。因此，我们永远不能指望通过平等的方式来比较距离和内半径。诸如 "if (distance_from_center == inner_radius) "这样的语句将会失败，除非我们运气特别好。相反，我们需要以一定的容忍度进行比较，通常的方法是写成`if  (std::abs(distance_from_center  。

// - inner_radius) <= tolerance)`，其中`tolerance'是比四舍五入大的某个小数字。问题是如何选择它。我们可以直接选择，比如说，`1e-10'，但这只适合于我们比较的对象是大小为1的情况。如果我们创建了一个单元大小为`1e+10'的网格，那么`1e-10'将远远低于四舍五入，就像以前一样，只有在我们特别幸运的情况下，比较才会成功。相反，使公差*相对于被比较对象的典型 "比例 "几乎总是有用的。在这里，"尺度 "是指内半径，或者是细胞的直径。我们选择前者，并将公差设置为 $10^{-6}$ 倍环形物的内半径。

              const double distance_from_center = 
                center.distance(cell->vertex(v)); 

              if (std::fabs(distance_from_center - inner_radius) <= 
                  1e-6 * inner_radius) 
                { 
                  cell->set_refine_flag(); 
                  break; 
                } 
            } 
        } 

// 现在我们已经标记了所有我们想要细化的单元格，我们让三角化实际做这个细化。这样做的函数的名字很长，因为我们也可以标记单元格进行粗化，该函数一次完成粗化和细化。

      triangulation.execute_coarsening_and_refinement(); 
    } 

// 最后，在这五次细化迭代之后，我们要再次将得到的网格写入文件，同样是SVG格式。这和上面的工作一样。

  std::ofstream out("grid-2.svg"); 
  GridOut       grid_out; 
  grid_out.write_svg(triangulation, out); 

  std::cout << "Grid written to grid-2.svg" << std::endl; 
} 

//  @sect3{The main function}  

// 最后是主函数。这里没有什么可做的，只是调用两个子函数，产生两个网格。

int main() 
{ 
  first_grid(); 
  second_grid(); 
} 



