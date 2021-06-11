CCTest_file/step-53.cc

/* --------------------------------------------------------------------- 
 * 
 * Copyright (C) 2014 - 2020 by the deal.II authors 
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
 * Authors: Wolfgang Bangerth, Texas A&M University, 2014 
 *          Luca Heltai, SISSA, 2014 
 *          D. Sarah Stamps, MIT, 2014 
 */ 



// 让我们从这里需要的包含文件开始。显然，我们需要描述三角形的文件（  <code>tria.h</code>  ），以及允许我们创建和输出三角形的文件（  <code>grid_generator.h</code>  和  <code>grid_out.h</code>  ）。此外，我们需要声明Manifold和ChartManifold类的头文件，我们将需要这些类来描述几何体（ <code>manifold.h</code> ）。然后我们还需要以下头文件中的 GridTools::transform() 函数；这个函数的用途将在我们使用它时讨论。

#include <deal.II/grid/tria.h> 
#include <deal.II/grid/grid_generator.h> 
#include <deal.II/grid/grid_out.h> 
#include <deal.II/grid/manifold.h> 
#include <deal.II/grid/grid_tools.h> 

// 其余的包含文件与读取地形数据有关。正如介绍中所解释的，我们将从一个文件中读取它，然后使用下面头文件中第一个声明的 Functions::InterpolatedUniformGridData  类。因为数据很大，所以我们读取的文件是以gzip压缩数据的形式存储的，我们利用BOOST提供的一些功能来直接读取gzipped数据。

#include <deal.II/base/function_lib.h> 

#include <boost/iostreams/filtering_stream.hpp> 
#include <boost/iostreams/filter/gzip.hpp> 
#include <boost/iostreams/device/file.hpp> 

#include <fstream> 
#include <iostream> 
#include <memory> 

// 上事的最后部分是打开一个命名空间，把所有东西都放进去，然后把dealii命名空间导入其中。

namespace Step53 
{ 
  using namespace dealii; 
// @sect3{Describing topography: AfricaTopography}  

// 这个程序的第一个重要部分是描述地形 $h(\hat phi,\hat \theta)$ 作为经度和纬度的函数的类。正如在介绍中所讨论的那样，我们在这里将使我们的生活更容易一些，不以最普遍的方式来写这个类，而是只为我们在这里感兴趣的特定目的来写：插值从一个非常具体的数据文件中获得的数据，该文件包含了关于世界上一个特定地区的信息，我们知道该地区的范围。

// 该类的总体布局已经在上面讨论过了。下面是它的声明，包括我们在初始化 <code>topography_data</code> 成员变量时需要的三个静态成员函数。

  class AfricaTopography 
  { 
  public: 
    AfricaTopography(); 

    double value(const double lon, const double lat) const; 

  private: 
    const Functions::InterpolatedUniformGridData<2> topography_data; 

    static std::vector<double> get_data(); 
  }; 

// 让我们来看看这个类的实现。该类的有趣部分是构造函数和 <code>value()</code> 函数。前者初始化了 Functions::InterpolatedUniformGridData 成员变量，我们将使用这个构造函数，它要求我们传入我们要插值的二维数据集的端点（这里由区间 $[-6.983333, 11.98333]$ 给出）。 ]，使用介绍中讨论的切换端点的技巧，和 $[25, 35.983333]$ ，都是以度数给出的），数据被分割成的区间数（纬度方向379，经度方向219，总共 $380\times 220$ 个数据点），和一个包含数据的表对象。然后，数据的大小当然是 $380\times 220$ ，我们通过提供一个迭代器给下面 std::vector 函数返回的 <code>get_data()</code> 对象的83,600个元素中的第一个来初始化它。注意，我们在这里调用的所有成员函数都是静态的，因为(i)它们不访问类的任何成员变量，(ii)因为它们是在对象没有完全初始化的时候调用的。

  AfricaTopography::AfricaTopography() 
    : topography_data({{std::make_pair(-6.983333, 11.966667), 
                        std::make_pair(25, 35.95)}}, 
                      {{379, 219}}, 
                      Table<2, double>(380, 220, get_data().begin())) 
  {} 

  double AfricaTopography::value(const double lon, const double lat) const 
  { 
    return topography_data.value( 
      Point<2>(-lat * 180 / numbers::PI, lon * 180 / numbers::PI)); 
  } 

// 唯一一个更有意义的函数是 <code>get_data()</code> 函数。它返回一个临时向量，其中包含描述高度的所有83600个数据点，并从文件 <code>topography.txt.gz</code> 中读取。因为文件被gzip压缩了，所以我们不能直接通过类型为 std::ifstream, 的对象来读取它，但在BOOST库中有一些方便的方法（见http:www.boost.org），允许我们从压缩的文件中读取，而不用先在磁盘上解压缩。其结果是，基本上，只是另一个输入流，就所有的实际目的而言，看起来就像我们一直使用的那些输入流。

// 当读取数据时，我们读取三列数据，但忽略了前两列。最后一列的数据被附加到一个数组中，我们返回的数组将被复制到 <code>topography_data</code> 的表中，并被初始化。由于BOOST.iostreams库在输入文件不存在、不可读或不包含正确的数据行数时没有提供非常有用的异常，我们捕捉它可能产生的所有异常并创建我们自己的异常。为此，在 <code>catch</code> 子句中，我们让程序运行到一个 <code>AssertThrow(false, ...)</code> 语句中。由于条件总是假的，这总是会触发一个异常。换句话说，这相当于写了 <code>throw ExcMessage("...")</code> ，但它也填补了异常对象中的某些字段，这些字段以后会被打印在屏幕上，识别出发生异常的函数、文件和行。

  std::vector<double> AfricaTopography::get_data() 
  { 
    std::vector<double> data; 

// 创建一个流，我们从gzipped数据中读取

    boost::iostreams::filtering_istream in; 
    in.push(boost::iostreams::basic_gzip_decompressor<>()); 
    in.push(boost::iostreams::file_source("topography.txt.gz")); 

    for (unsigned int line = 0; line < 83600; ++line) 
      { 
        try 
          { 
            double lat, lon, elevation; 
            in >> lat >> lon >> elevation; 

            data.push_back(elevation); 
          } 
        catch (...) 
          { 
            AssertThrow(false, 
                        ExcMessage("Could not read all 83,600 data points " 
                                   "from the file <topography.txt.gz>!")); 
          } 
      } 

    return data; 
  } 
// @sect3{Describing the geometry: AfricaGeometry}  

// 下面的类是本程序的主类。它的结构已经在介绍中详细描述过了，不需要再多做介绍。

  class AfricaGeometry : public ChartManifold<3, 3> 
  { 
  public: 
    virtual Point<3> pull_back(const Point<3> &space_point) const override; 

    virtual Point<3> push_forward(const Point<3> &chart_point) const override; 

    virtual std::unique_ptr<Manifold<3, 3>> clone() const override; 

  private: 
    static const double R; 
    static const double ellipticity; 

    const AfricaTopography topography; 

    Point<3> push_forward_wgs84(const Point<3> &phi_theta_d) const; 
    Point<3> pull_back_wgs84(const Point<3> &x) const; 

    Point<3> push_forward_topo(const Point<3> &phi_theta_d_hat) const; 
    Point<3> pull_back_topo(const Point<3> &phi_theta_d) const; 
  }; 

  const double AfricaGeometry::R           = 6378137; 
  const double AfricaGeometry::ellipticity = 8.1819190842622e-2; 

// 如果你读过介绍，实现起来也是非常简单的。特别是，回拉和前推函数都只是WGS 84和地形图映射各自函数的串联。

  Point<3> AfricaGeometry::pull_back(const Point<3> &space_point) const 
  { 
    return pull_back_topo(pull_back_wgs84(space_point)); 
  } 

  Point<3> AfricaGeometry::push_forward(const Point<3> &chart_point) const 
  { 
    return push_forward_wgs84(push_forward_topo(chart_point)); 
  } 

// 下一个函数是Manifold基类的接口所要求的，它允许克隆AfricaGeometry类。注意，虽然该函数返回一个  `std::unique_ptr<Manifold<3,3>>`,  我们在内部创建了一个 `unique_ptr<AfricaGeometry>`。换句话说，这个库需要一个指向基类的指针，我们通过创建一个指向派生类的指针来提供这个指针。

  std::unique_ptr<Manifold<3, 3>> AfricaGeometry::clone() const 
  { 
    return std::make_unique<AfricaGeometry>(); 
  } 

// 下面的两个函数就定义了对应于地球WGS84参考形状的正向和反向变换。正向变换遵循介绍中所示的公式。反变换要复杂得多，至少不是直观的。它还存在一个问题，即它返回一个角度，在函数结束时，如果它应该从那里逃出来，我们需要将其夹回区间 $[0,2\pi]$ 。

  Point<3> AfricaGeometry::push_forward_wgs84(const Point<3> &phi_theta_d) const 
  { 
    const double phi   = phi_theta_d[0]; 
    const double theta = phi_theta_d[1]; 
    const double d     = phi_theta_d[2]; 

    const double R_bar = R / std::sqrt(1 - (ellipticity * ellipticity * 
                                            std::sin(theta) * std::sin(theta))); 

    return {(R_bar + d) * std::cos(phi) * std::cos(theta), 
            (R_bar + d) * std::sin(phi) * std::cos(theta), 
            ((1 - ellipticity * ellipticity) * R_bar + d) * std::sin(theta)}; 
  } 

  Point<3> AfricaGeometry::pull_back_wgs84(const Point<3> &x) const 
  { 
    const double b   = std::sqrt(R * R * (1 - ellipticity * ellipticity)); 
    const double ep  = std::sqrt((R * R - b * b) / (b * b)); 
    const double p   = std::sqrt(x(0) * x(0) + x(1) * x(1)); 
    const double th  = std::atan2(R * x(2), b * p); 
    const double phi = std::atan2(x(1), x(0)); 
    const double theta = 
      std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th), 3), 
                 (p - 
                  (ellipticity * ellipticity * R * std::pow(std::cos(th), 3)))); 
    const double R_bar = 
      R / (std::sqrt(1 - ellipticity * ellipticity * std::sin(theta) * 
                           std::sin(theta))); 
    const double R_plus_d = p / std::cos(theta); 

    Point<3> phi_theta_d; 
    if (phi < 0) 
      phi_theta_d[0] = phi + 2 * numbers::PI; 
    else if (phi > 2 * numbers::PI) 
      phi_theta_d[0] = phi - 2 * numbers::PI; 
    else 
      phi_theta_d[0] = phi; 
    phi_theta_d[1] = theta; 
    phi_theta_d[2] = R_plus_d - R_bar; 
    return phi_theta_d; 
  } 

// 与此相反，地形变换完全按照介绍中的描述进行。因此，没有什么可以补充的。

  Point<3> 
  AfricaGeometry::push_forward_topo(const Point<3> &phi_theta_d_hat) const 
  { 
    const double d_hat = phi_theta_d_hat[2]; 
    const double h = topography.value(phi_theta_d_hat[0], phi_theta_d_hat[1]); 
    const double d = d_hat + (d_hat + 500000) / 500000 * h; 
    return {phi_theta_d_hat[0], phi_theta_d_hat[1], d}; 
  } 

  Point<3> AfricaGeometry::pull_back_topo(const Point<3> &phi_theta_d) const 
  { 
    const double d     = phi_theta_d[2]; 
    const double h     = topography.value(phi_theta_d[0], phi_theta_d[1]); 
    const double d_hat = 500000 * (d - h) / (500000 + h); 
    return {phi_theta_d[0], phi_theta_d[1], d_hat}; 
  } 
// @sect3{Creating the mesh}  

// 在描述了几何体的属性之后，现在是处理用于离散它的网格的时候了。为此，我们为几何体和三角形创建对象，然后继续创建一个与参考域 $1\times 2\times 1$ 相对应的 $\hat U=[26,35]\times[-10,5]\times[-500000,0]$ 矩形网格。我们选择这个数目的细分，因为它导致了单元格大致上像立方体，而不是在某个方向上被拉伸。

// 当然，我们实际上对参考域的网格划分不感兴趣。我们感兴趣的是对真实域的网格划分。因此，我们将使用 GridTools::transform() 函数，它只是根据一个给定的变换来移动三角形的每个点。它想要的变换函数是一个将参考域中的一个点作为其单一参数的函数，并返回我们想要映射到的域中的相应位置。当然，这正是我们使用的几何学的前推函数。我们用一个lambda函数来包装它，以获得转换所需的那种函数对象。

  void run() 
  { 
    AfricaGeometry   geometry; 
    Triangulation<3> triangulation; 

    { 
      const Point<3> corner_points[2] = { 
        Point<3>(26 * numbers::PI / 180, -10 * numbers::PI / 180, -500000), 
        Point<3>(35 * numbers::PI / 180, 5 * numbers::PI / 180, 0)}; 
      std::vector<unsigned int> subdivisions(3); 
      subdivisions[0] = 1; 
      subdivisions[1] = 2; 
      subdivisions[2] = 1; 
      GridGenerator::subdivided_hyper_rectangle( 
        triangulation, subdivisions, corner_points[0], corner_points[1], true); 

      GridTools::transform( 
        [&geometry](const Point<3> &chart_point) { 
          return geometry.push_forward(chart_point); 
        }, 
        triangulation); 
    } 

// 下一步是向三角计算说明，在细化网格时，每当需要一个新的点时，都要使用我们的几何对象。我们通过告诉三角计算对所有流形指示器为零的物体使用我们的几何体，然后继续用流形指示器为零标记所有单元及其边界面和边。这确保了三角计算在每次需要新的顶点时都会参考我们的几何对象。由于流形指标是由母体继承给子体的，这也会在几个递归细化步骤之后发生。

    triangulation.set_manifold(0, geometry); 
    for (const auto &cell : triangulation.active_cell_iterators()) 
      cell->set_all_manifold_ids(0); 

// 最后一步是在最初的 $1\times 2\times 1$ 粗略网格之外细化该网格。我们可以在全局范围内细化若干次，但由于本教程程序的目的，我们实际上只对靠近表面的情况感兴趣，所以我们只是对所有在边界上有一个指标为5的面的单元进行6次细化。在我们上面使用的 GridGenerator::subdivided_hyper_rectangle() 函数的文档中查找，发现边界指标5对应于域的顶面（这就是上面调用 GridGenerator::subdivided_hyper_rectangle() 的最后一个 <code>true</code> 参数的含义：通过给每个边界分配一个独特的边界指标来给边界 "着色"）。

    for (unsigned int i = 0; i < 6; ++i) 
      { 
        for (const auto &cell : triangulation.active_cell_iterators()) 
          for (const auto &face : cell->face_iterators()) 
            if (face->boundary_id() == 5) 
              { 
                cell->set_refine_flag(); 
                break; 
              } 
        triangulation.execute_coarsening_and_refinement(); 

        std::cout << "Refinement step " << i + 1 << ": " 
                  << triangulation.n_active_cells() << " cells, " 
                  << GridTools::minimal_cell_diameter(triangulation) / 1000 
                  << "km minimal cell diameter" << std::endl; 
      } 

// 做完这一切，我们现在可以将网格输出到一个自己的文件中。

    const std::string filename = "mesh.vtu"; 
    std::ofstream     out(filename); 
    GridOut           grid_out; 
    grid_out.write_vtu(triangulation, out); 
  } 
} // namespace Step53 

//  @sect3{The main function}  

// 最后是主函数，它采用了从  step-6  开始的所有教程程序中使用的相同方案。这里没有什么可做的，只需要调用单一的  <code>run()</code>  函数。

int main() 
{ 
  try 
    { 
      Step53::run(); 
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


