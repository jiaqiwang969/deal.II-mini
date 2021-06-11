CCTest_file/step-49.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2013 - 2021 by the deal.II authors 
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
 * Author: Timo Heister, Texas A&M University, 2013 
 */ 



// 这个教程程序很奇怪，与其他大多数步骤不同，介绍中已经提供了关于如何使用各种策略来生成网格的大部分信息。因此，这里没有什么需要评论的，我们在代码中穿插了相对较少的文字。从本质上讲，这里的代码只是提供了一个在介绍中已经描述过的内容的参考实现。

//  @sect3{Include files}  
#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_tools.h> 
#include <deal.II/grid/manifold_lib.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/grid_in.h> 

#include <iostream> 
#include <fstream> 

#include <map> 

using namespace dealii; 
// @sect3{Generating output for a given mesh}  

// 下面的函数为我们将在本程序的剩余部分中生成的任何网格生成一些输出。特别是，它生成了以下信息。



// - 一些关于这个网格所处的空间维数和它的单元数的一般信息。

// - 使用每个边界指标的边界面的数量，这样就可以和我们预期的情况进行比较。

// 最后，该函数将网格输出为VTU格式，可以方便地在Paraview或VisIt中进行可视化。

template <int dim> 
void print_mesh_info(const Triangulation<dim> &triangulation, 
                     const std::string &       filename) 
{ 
  std::cout << "Mesh info:" << std::endl 
            << " dimension: " << dim << std::endl 
            << " no. of cells: " << triangulation.n_active_cells() << std::endl; 

// 接下来循环所有单元格的所有面，找出每个边界指标的使用频率（请记住，如果你访问一个不存在的 std::map 对象的元素，它将被隐式创建并默认初始化--在当前情况下为零--然后我们再将其增加）。

  { 
    std::map<types::boundary_id, unsigned int> boundary_count; 
    for (const auto &face : triangulation.active_face_iterators()) 
      if (face->at_boundary()) 
        boundary_count[face->boundary_id()]++; 

    std::cout << " boundary indicators: "; 
    for (const std::pair<const types::boundary_id, unsigned int> &pair : 
         boundary_count) 
      { 
        std::cout << pair.first << "(" << pair.second << " times) "; 
      } 
    std::cout << std::endl; 
  } 

// 最后，产生一个网格的图形表示到一个输出文件。

  std::ofstream out(filename); 
  GridOut       grid_out; 
  grid_out.write_vtu(triangulation, out); 
  std::cout << " written to " << filename << std::endl << std::endl; 
} 
// @sect3{Main routines}  
// @sect4{grid_1: Loading a mesh generated by gmsh}  

// 在这第一个例子中，我们展示了如何加载我们在介绍中讨论过的如何生成的网格。这与 step-5 中加载网格的模式相同，尽管那里是以不同的文件格式（UCD而不是MSH）编写。

void grid_1() 
{ 
  Triangulation<2> triangulation; 

  GridIn<2> gridin; 
  gridin.attach_triangulation(triangulation); 
  std::ifstream f("example.msh"); 
  gridin.read_msh(f); 

  print_mesh_info(triangulation, "grid-1.vtu"); 
} 
// @sect4{grid_2: Merging triangulations}  

// 在这里，我们首先创建两个三角形，然后将它们合并成一个。 正如介绍中所讨论的，必须确保共同界面的顶点位于相同的坐标上。

void grid_2() 
{ 
  Triangulation<2> tria1; 
  GridGenerator::hyper_cube_with_cylindrical_hole(tria1, 0.25, 1.0); 

  Triangulation<2>          tria2; 
  std::vector<unsigned int> repetitions(2); 
  repetitions[0] = 3; 
  repetitions[1] = 2; 
  GridGenerator::subdivided_hyper_rectangle(tria2, 
                                            repetitions, 
                                            Point<2>(1.0, -1.0), 
                                            Point<2>(4.0, 1.0)); 

  Triangulation<2> triangulation; 
  GridGenerator::merge_triangulations(tria1, tria2, triangulation); 

  print_mesh_info(triangulation, "grid-2.vtu"); 
} 
// @sect4{grid_3: Moving vertices}  

// 在这个函数中，我们移动一个网格的顶点。这比人们通常想象的要简单：如果你用 <code>cell-@>vertex(i)</code> 询问一个单元格的 <code>i</code> 的顶点的坐标，它不只是提供这个顶点的位置，实际上是对存储这些坐标的位置的引用。然后我们可以修改存储在那里的值。

// 所以这就是我们在这个函数的第一部分所做的。我们创建一个几何形状为 $[-1,1]^2$ 的正方形，在原点处有一个半径为0.25的圆孔。然后我们在所有单元格和所有顶点上循环，如果一个顶点的 $y$ 坐标等于1，我们就把它向上移动0.5。

// 注意，这种程序通常不是这样工作的，因为通常会多次遇到相同的顶点，并且可能会多次移动它们。它在这里起作用是因为我们根据顶点的几何位置来选择要使用的顶点，而移动过一次的顶点在未来将无法通过这个测试。解决这个问题的一个更普遍的方法是保留一个 std::set ，即那些我们已经移动过的顶点索引（我们可以用 <code>cell-@>vertex_index(i)</code> 获得，并且只移动那些索引还不在这个集合中的顶点。

void grid_3() 
{ 
  Triangulation<2> triangulation; 
  GridGenerator::hyper_cube_with_cylindrical_hole(triangulation, 0.25, 1.0); 

  for (const auto &cell : triangulation.active_cell_iterators()) 
    { 
      for (const auto i : cell->vertex_indices()) 
        { 
          Point<2> &v = cell->vertex(i); 
          if (std::abs(v(1) - 1.0) < 1e-5) 
            v(1) += 0.5; 
        } 
    } 

// 在第二步，我们将对网格进行两次细化。为了正确做到这一点，我们应该沿着以原点为中心的圆的表面在内部边界上放置新的点。幸运的是， GridGenerator::hyper_cube_with_cylindrical_hole 已经在内部边界上附加了一个Manifold对象，所以我们不需要做任何事情，只需要细化网格（参见<a href="#Results">results section</a>中一个完全可行的例子，我们 <em> 做 </em> 附加一个Manifold对象）。

  triangulation.refine_global(2); 
  print_mesh_info(triangulation, "grid-3.vtu"); 
} 

// 如上图所示，做事有一个障碍。如果像这里所示的那样移动边界上的节点，由于内部的节点没有被移动，所以经常会出现内部的单元被严重扭曲的情况。在目前的情况下，这并不是一个很大的问题，因为当节点被移动时，网格并不包含任何内部节点--它是粗略的网格，而且恰好所有的顶点都在边界上。还有一种情况是，我们在这里的移动，与平均单元的大小相比，并没有太大影响。然而，有时我们确实想把顶点移动一段距离，在这种情况下，我们也需要移动内部节点。一个自动完成的方法是调用函数 GridTools::laplace_transform ，该函数接收一组转换后的顶点坐标并移动所有其他的顶点，使产生的网格在某种意义上有一个小的变形。

//  @sect4{grid_4: Demonstrating extrude_triangulation}  

// 这个例子从前面的函数中获取初始网格，并简单地将其挤压到第三空间维度。

void grid_4() 
{ 
  Triangulation<2> triangulation; 
  Triangulation<3> out; 
  GridGenerator::hyper_cube_with_cylindrical_hole(triangulation, 0.25, 1.0); 

  GridGenerator::extrude_triangulation(triangulation, 3, 2.0, out); 
  print_mesh_info(out, "grid-4.vtu"); 
} 
// @sect4{grid_5: Demonstrating GridTools::transform, part 1}  

// 这个例子和下一个例子首先创建一个网格，然后根据一个函数移动网格的每个节点，这个函数接收一个点并返回一个映射的点。在这个例子中，我们转换  $(x,y) \mapsto (x,y+\sin(\pi x/5))$  。

//  GridTools::transform()  需要一个三角形和一个参数，这个参数可以像一个函数一样被调用，接收一个点并返回一个点。有不同的方式来提供这样一个参数。它可以是一个函数的指针；它可以是一个具有`operator()`的类的对象；它可以是一个lambda函数；或者它可以是任何通过 <code>std::function@<Point@<2@>(const Point@<2@>)@></code> 对象描述的东西。

// 更现代的方法是使用一个接受一个点并返回一个点的lambda函数，这就是我们在下面所做的。

void grid_5() 
{ 
  Triangulation<2>          triangulation; 
  std::vector<unsigned int> repetitions(2); 
  repetitions[0] = 14; 
  repetitions[1] = 2; 
  GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                            repetitions, 
                                            Point<2>(0.0, 0.0), 
                                            Point<2>(10.0, 1.0)); 

  GridTools::transform( 
    [](const Point<2> &in) { 
      return Point<2>(in[0], in[1] + std::sin(numbers::PI * in[0] / 5.0)); 
    }, 
    triangulation); 
  print_mesh_info(triangulation, "grid-5.vtu"); 
} 

//  @sect4{grid_6: Demonstrating GridTools::transform, part 2}  

// 在第二个例子中，我们将使用映射  $(x,y) \mapsto (x,\tanh(2y)/\tanh(2))$  将原始网格中的点转换为新的网格。为了使事情更有趣，而不是像前面的例子那样在一个单一的函数中完成，我们在这里创建一个具有  <code>operator()</code>  的对象，这个对象将被  GridTools::transform.  所调用。当然，这个对象实际上可能要复杂得多：这个对象可能有成员变量，在计算顶点的新位置时起作用。

struct Grid6Func 
{ 
  double trans(const double y) const 
  { 
    return std::tanh(2 * y) / tanh(2); 
  } 

  Point<2> operator()(const Point<2> &in) const 
  { 
    return {in(0), trans(in(1))}; 
  } 
}; 

void grid_6() 
{ 
  Triangulation<2>          triangulation; 
  std::vector<unsigned int> repetitions(2); 
  repetitions[0] = repetitions[1] = 40; 
  GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                            repetitions, 
                                            Point<2>(0.0, 0.0), 
                                            Point<2>(1.0, 1.0)); 

  GridTools::transform(Grid6Func(), triangulation); 
  print_mesh_info(triangulation, "grid-6.vtu"); 
} 
// @sect4{grid_7: Demonstrating distort_random}  

// 在这最后一个例子中，我们创建了一个网格，然后通过随机扰动使其（内部）顶点变形。这不是你想在生产计算中做的事情（因为在具有 "良好形状 "单元的网格上的结果通常比在 GridTools::distort_random()), 产生的变形单元上的结果要好，但这是一个有用的工具，可以测试离散化和代码，确保它们不会因为网格恰好是均匀结构和支持超级收敛特性而意外地工作。

void grid_7() 
{ 
  Triangulation<2>          triangulation; 
  std::vector<unsigned int> repetitions(2); 
  repetitions[0] = repetitions[1] = 16; 
  GridGenerator::subdivided_hyper_rectangle(triangulation, 
                                            repetitions, 
                                            Point<2>(0.0, 0.0), 
                                            Point<2>(1.0, 1.0)); 

  GridTools::distort_random(0.3, triangulation, true); 
  print_mesh_info(triangulation, "grid-7.vtu"); 
} 
// @sect3{The main function}  

// 最后是主函数。这里没有什么可做的，只是调用我们上面写的所有各种函数。

int main() 
{ 
  try 
    { 
      grid_1(); 
      grid_2(); 
      grid_3(); 
      grid_4(); 
      grid_5(); 
      grid_6(); 
      grid_7(); 
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
} 


