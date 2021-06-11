CCTest_file/step-54.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2009 - 2021 by the deal.II authors 
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
 *  Authors: Andrea Mola, Luca Heltai, 2014 
 */ 


// @sect3{Include files}  

// 我们首先包括一堆我们将在程序的各个部分使用的文件。它们中的大多数已经在以前的教程中讨论过了。

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_in.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/numerics/data_out.h> 
#include <deal.II/numerics/vector_tools.h> 

// 这些是opencascade支持类和函数的头文件。注意，只有当你在编译deal.II库时支持OpenCASCADE，即在deal.II配置过程中调用 <code>-DDEAL_II_WITH_OPENCASCADE=ON</code> 和 <code>-DOPENCASCADE_DIR=/path/to/your/opencascade/installation</code> 时，这些将包含合理的数据。

#include <deal.II/opencascade/manifold_lib.h> 
#include <deal.II/opencascade/utilities.h> 

// 最后，几个C++标准头文件

#include <cmath> 
#include <iostream> 
#include <fstream> 
#include <string> 

// 我们将程序的其他部分隔离在自己的命名空间中

namespace Step54 
{ 
  using namespace dealii; 

//  @sect3{The TriangulationOnCAD class}  

// 这是主类。它真正做的是存储输入和输出文件的名称，以及一个三角图。然后，它提供了一个函数，可以从一个粗略的网格中生成这样一个三角形，使用介绍中所讨论的策略之一，并在类的顶部的枚举类型中列出。

// 这个类的成员函数与你在其他大多数教程程序中可以找到的类似，都是在模拟的网格设置阶段。

  class TriangulationOnCAD 
  { 
  public: 
    enum ProjectionType 
    { 
      NormalProjection       = 0, 
      DirectionalProjection  = 1, 
      NormalToMeshProjection = 2 
    }; 

    TriangulationOnCAD( 
      const std::string &  initial_mesh_filename, 
      const std::string &  cad_file_name, 
      const std::string &  output_filename, 
      const ProjectionType surface_projection_kind = NormalProjection); 

    void run(); 

  private: 
    void read_domain(); 

    void refine_mesh(); 

    void output_results(const unsigned int cycle); 

    Triangulation<2, 3> tria; 

    const std::string initial_mesh_filename; 
    const std::string cad_file_name; 
    const std::string output_filename; 

    const ProjectionType surface_projection_kind; 
  }; 
// @sect4{TriangulationOnCAD::TriangulationOnCAD}  

// TriangulationOnCAD类的构造函数非常简单。输入参数是输入和输出文件名的字符串，以及决定在网格细化循环中使用哪种曲面投影仪的枚举类型（详见下文）。

  TriangulationOnCAD::TriangulationOnCAD( 
    const std::string &  initial_mesh_filename, 
    const std::string &  cad_file_name, 
    const std::string &  output_filename, 
    const ProjectionType surface_projection_kind) 
    : initial_mesh_filename(initial_mesh_filename) 
    , cad_file_name(cad_file_name) 
    , output_filename(output_filename) 
    , surface_projection_kind(surface_projection_kind) 
  {} 
// @sect4{TriangulationOnCAD::read_domain}  

// 下面的函数代表了这个程序的核心。 在这个函数中，我们导入CAD形状，在此基础上生成并完善我们的三角测量。我们假设CAD曲面包含在 @p cad_file_name 文件中（我们在输入目录中提供了一个名为 "input/DTMB-5415_bulbous_bow.iges "的IGES文件例子，它代表了一艘船的球形船头）。几个凸和凹的高曲率区域的存在使得我们提供的几何体成为一个特别有意义的例子。

// 在导入船首表面后，我们提取了一些构成它的曲线和曲面，并使用它们来生成一组投影仪。这些投影仪定义了三角法在细化单元时必须遵循的规则，以定位每个新节点。

// 为了初始化Triangulation，就像在以前的教程中一样，我们导入一个以VTK格式保存的现有网格。在这里我们假设用户已经在外部生成了一个粗略的网格，与IGES的几何图形相匹配。在编写本教程的时候，deal.II库并不自动支持生成这样的网格，但是有一些工具可以从CAD文件开始为你提供合理的初始网格。在我们的例子中，导入的网格是由一个四边形单元组成的，其顶点被放置在CAD形状上。

// 在导入IGES几何体和初始网格后，我们将之前讨论过的投影仪分配给每个需要在CAD表面上进行细化的边和单元。

// 在本教程中，我们将测试介绍中所描述的三种不同的CAD表面投影器，并将分析每一种投影器所得到的结果。 如前所述，这些投影策略中的每一种都是在不同的类中实现的，这些类型的对象可以用 Triangulation::set_manifold 的方法分配给一个三角形。
//然后，
//下面的函数首先导入给定的CAD文件。函数的参数是一个包含所需文件名的字符串，以及一个比例因子。在这个例子中，比例因子被设置为1e-3，因为原始几何体是以毫米为单位的（这是大多数IGES文件的典型计量单位），而我们更喜欢以米为单位工作。 该函数的输出是一个OpenCASCADE通用拓扑形状类的对象，即一个 @p TopoDS_Shape. 。
  void TriangulationOnCAD::read_domain() 
  { 
    TopoDS_Shape bow_surface = OpenCASCADE::read_IGES(cad_file_name, 1e-3); 

// 每个CAD几何对象都被定义了一个公差，它表示其位置可能的不精确性。例如，顶点的公差 @p tol 表示它可以位于以标称位置为中心，半径为 @p tol. 的球体中的任何一点。

// 下面的方法是提取给定形状的公差，并使其变大一些以避免麻烦。

    const double tolerance = OpenCASCADE::get_shape_tolerance(bow_surface) * 5; 

// 我们现在要从通用形状中提取一组复合子形状。特别是，CAD文件的每个面都是由类型为 @p TopoDS_Wire, 的修剪曲线组成的，它是构成曲面边界的 @p TopoDS_Edges 的集合，以及曲面本身的NURBS描述。我们将使用线型投影仪将我们的三角形的边界与划定曲面的线联系起来。 为了提取所有的复合子形状，如线、壳或实体，我们求助于OpenCASCADE命名空间的一种方法。  OpenCASCADE::extract_compound_shapes 的输入是一个形状和一组空的 std::vectors 子形状，它将被填入在给定拓扑形状中发现的所有复合形状。

    std::vector<TopoDS_Compound>  compounds; 
    std::vector<TopoDS_CompSolid> compsolids; 
    std::vector<TopoDS_Solid>     solids; 
    std::vector<TopoDS_Shell>     shells; 
    std::vector<TopoDS_Wire>      wires; 

    OpenCASCADE::extract_compound_shapes( 
      bow_surface, compounds, compsolids, solids, shells, wires); 

// 接下来的几个步骤比较熟悉，允许我们从外部VTK文件中导入一个现有的网格，并将其转换为一个交易三角图。

    std::ifstream in; 

    in.open(initial_mesh_filename); 

    GridIn<2, 3> gi; 
    gi.attach_triangulation(tria); 
    gi.read_vtk(in); 

// 我们输出这个初始网格，将其保存为细化步骤0。

    output_results(0); 

// 导入的网格有一个位于三维空间中的单一的二维单元。我们现在要确保它是根据上面导入的CAD几何图形进行细化的。为此，我们得到一个单元的迭代器，并给它分配manifold_id 1（见  @ref GlossManifoldIndicator  "这个词汇表条目"）。我们还得到了一个指向其四个面的迭代器，并为每个面分配了manifold_id 2。

    Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active(); 
    cell->set_manifold_id(1); 

    for (const auto &face : cell->face_iterators()) 
      face->set_manifold_id(2); 

// 一旦CAD几何体和初始网格都被导入和消化，我们就用CAD的曲面和曲线来定义投影仪，并将它们分配给刚才指定的流形ID。

// 使用我们的CAD文件中的单线来定义第一个投影仪。 ArclengthProjectionLineManifold将确保位于导线上的每条网格边缘都被细化为一个位于导线上的点，并将其分割为两个位于边缘顶点之间的相等弧线。我们首先检查线的向量是否至少包含一个元素，然后为它创建一个Manifold对象。

// 一旦投影仪被创建，我们就把它分配给三角形的所有部分，manifold_id = 2。

    Assert( 
      wires.size() > 0, 
      ExcMessage( 
        "I could not find any wire in the CAD file you gave me. Bailing out.")); 

    OpenCASCADE::ArclengthProjectionLineManifold<2, 3> line_projector( 
      wires[0], tolerance); 

    tria.set_manifold(2, line_projector); 

// 根据构造函数的 @p surface_projection_kind 选项所指定的内容来创建表面投影仪。特别是，如果surface_projection_kind的值等于 @p NormalProjection, ，我们选择 OpenCASCADE::NormalProjectionManifold. ，那么新的网格点最初将在所考虑的单元格/边的arycenter处生成，然后沿其法线方向投影到CAD表面。 NormalProjectionManifold构造函数只需要一个形状和一个公差，然后我们把它分配给三角结构，用于所有具有id 1的流形的部分。

    switch (surface_projection_kind) 
      { 
        case NormalProjection: 
          { 
            OpenCASCADE::NormalProjectionManifold<2, 3> normal_projector( 
              bow_surface, tolerance); 
            tria.set_manifold(1, normal_projector); 

            break; 
          } 
// @p If  surface_projection_kind值为 @p DirectionalProjection,  我们选择 OpenCASCADE::DirectionalProjectionManifold 类。新的网格点将在所考虑的单元格/边的arycenter处初始生成，然后沿着 OpenCASCADE::DirectionalProjectionManifold 构造函数指定的方向投影到CAD表面上。在这个例子中，投影是沿着Y轴进行的。

        case DirectionalProjection: 
          { 
            OpenCASCADE::DirectionalProjectionManifold<2, 3> 
              directional_projector(bow_surface, 
                                    Point<3>(0.0, 1.0, 0.0), 
                                    tolerance); 
            tria.set_manifold(1, directional_projector); 

            break; 
          } 

// 作为第三个选项，如果 @p surface_projection_kind 的值是 @p NormalToMeshProjection, ，我们选择 OpenCASCADE::NormalToMeshProjectionManifold.  新的网格点将再次在所考虑的单元/边的arycenter处初始生成，然后沿着一个估计为网格法线方向的方向投影到CAD表面。 OpenCASCADE::NormalToMeshProjectionManifold  构造函数只需要一个形状（至少包含一个面）和一个公差。

        case NormalToMeshProjection: 
          { 
            OpenCASCADE::NormalToMeshProjectionManifold<2, 3> 
              normal_to_mesh_projector(bow_surface, tolerance); 
            tria.set_manifold(1, normal_to_mesh_projector); 

            break; 
          } 

// 最后，我们使用良好的软件清洁性，确保这真的涵盖了 @p case 语句的所有可能选项。如果我们得到任何其他的值，我们就直接中止程序。

        default: 
          AssertThrow(false, ExcInternalError()); 
      } 
  } 
// @sect4{TriangulationOnCAD::refine_mesh}  

// 这个函数是全局细化网格的。在其他教程中，它通常也会分配自由度，并调整矩阵和向量的大小。这里没有进行这些工作，因为我们没有在生成的三角形上运行任何模拟。

// 虽然这个函数看起来很简单，但这是我们对这个教程程序感兴趣的大部分工作的实际发生地点。特别是，在完善定义船体表面的四边形和直线时，Triangulation类将询问我们分配给处理单个流形ID的各种对象，以确定新顶点的位置。

  void TriangulationOnCAD::refine_mesh() 
  { 
    tria.refine_global(1); 
  } 

//  @sect4{TriangulationOnCAD::output_results}  

// 输出我们的计算结果是一个相当机械的任务。这个函数的所有组成部分之前已经讨论过了。

  void TriangulationOnCAD::output_results(const unsigned int cycle) 
  { 
    const std::string filename = 
      (output_filename + "_" + Utilities::int_to_string(cycle) + ".vtk"); 
    std::ofstream logfile(filename); 
    GridOut       grid_out; 
    grid_out.write_vtk(tria, logfile); 
  } 
// @sect4{TriangulationOnCAD::run}  

// 这是主函数。它应该是不言自明的。

  void TriangulationOnCAD::run() 
  { 
    read_domain(); 

    const unsigned int n_cycles = 5; 
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle) 
      { 
        refine_mesh(); 
        output_results(cycle + 1); 
      } 
  } 
} // namespace Step54 
// @sect3{The main() function}  

// 这是本程序的主要功能。它的基本结构与之前所有的教程程序一样，但通过新顶点放置的三种可能性来运行主类。

int main() 
{ 
  try 
    { 
      using namespace Step54; 

      const std::string in_mesh_filename = "input/initial_mesh_3d.vtk"; 
      const std::string cad_file_name    = "input/DTMB-5415_bulbous_bow.iges"; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in direction normal to CAD surface" 
                << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::string        out_mesh_filename = ("3d_mesh_normal_projection"); 
      TriangulationOnCAD tria_on_cad_norm(in_mesh_filename, 
                                          cad_file_name, 
                                          out_mesh_filename, 
                                          TriangulationOnCAD::NormalProjection); 
      tria_on_cad_norm.run(); 
      std::cout << "----------------------------------------------------------" 
                << std::endl;
      std::cout << std::endl; 
      std::cout << std::endl; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in y-axis direction" << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      out_mesh_filename = ("3d_mesh_directional_projection"); 
      TriangulationOnCAD tria_on_cad_dir( 
        in_mesh_filename, 
        cad_file_name, 
        out_mesh_filename, 
        TriangulationOnCAD::DirectionalProjection);  
      tria_on_cad_dir.run();  
      std::cout << "----------------------------------------------------------"  
                << std::endl; 
      std::cout << std::endl; 
      std::cout << std::endl; 

      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << "Testing projection in direction normal to mesh elements" 
                << std::endl; 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      out_mesh_filename = ("3d_mesh_normal_to_mesh_projection"); 
      TriangulationOnCAD tria_on_cad_norm_to_mesh( 
        in_mesh_filename, 
        cad_file_name, 
        out_mesh_filename, 
        TriangulationOnCAD::NormalToMeshProjection); 
      tria_on_cad_norm_to_mesh.run(); 
      std::cout << "----------------------------------------------------------" 
                << std::endl; 
      std::cout << std::endl; 
      std::cout << std::endl; 
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

