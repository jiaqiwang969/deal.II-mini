//include/deal.II-translator/opencascade/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_occ_utilities_h
#  define dealii_occ_utilities_h

#  include <deal.II/base/config.h>

#  ifdef DEAL_II_WITH_OPENCASCADE

#    include <deal.II/base/point.h>

#    include <deal.II/fe/mapping_q1.h>

#    include <deal.II/grid/tria.h>

#    include <string>

// opencascade needs "HAVE_CONFIG_H" to be exported...
#    define HAVE_CONFIG_H
#    include <IFSelect_ReturnStatus.hxx>
#    include <TopoDS_CompSolid.hxx>
#    include <TopoDS_Compound.hxx>
#    include <TopoDS_Edge.hxx>
#    include <TopoDS_Face.hxx>
#    include <TopoDS_Shape.hxx>
#    include <TopoDS_Shell.hxx>
#    include <TopoDS_Solid.hxx>
#    include <TopoDS_Vertex.hxx>
#    include <TopoDS_Wire.hxx>
#    include <gp_Pnt.hxx>
#    undef HAVE_CONFIG_H



DEAL_II_NAMESPACE_OPEN

/**
 * 我们在这个命名空间中收集所有对OpenCASCADE实体进行操作的实用程序。OpenCASCADE把每个物体分成了一个拓扑描述和一个几何实体。基本的拓扑描述是一个TopoDS_Shape。TopoDS_Shape是轻型物体，可以被复制。最接近deal.II的类似物是一个TriaIterator。
 * OpenCASCADE的拓扑结构是参照STEP标准ISO-10303-42设计的。
 * 该结构是一个定向的单行图，父指其子，没有反向引用。抽象结构是以TopoDS包的C++类实现的。一个TopoDS_Shape是通过值来操作的，包含3个字段：位置、方向和一个myTShape句柄（TopoDS_TShape类型的）。根据OpenCASCADE的文档，myTShape和Location是用来在各种形状之间共享数据以节省内存。例如，一条属于两个面的边有相等的位置和myTShape字段，但有不同的方向（在一个面的上下文中是向前的，而在另一个面中是相反的）。
 * 有效的形状包括其他形状的集合、实体、面、边、顶点，等等。
 * 一旦有了拓扑描述，如果可以创建一个具体的几何对象，BRep类允许人们从一个形状中提取实际的几何信息。
 * 这是通过实现边界表示模型（来自BRep包）的人继承TopoDS包中的抽象拓扑学类来实现的。只有3种类型的拓扑对象具有几何表示法
 *
 * -顶点、边缘和面。
 * 每个TopoDS_Shape都可以被查询到它是什么类型的形状，而实际的几何对象，如曲面、曲线或点，可以用BRepTools提取。
 * 在这个命名空间中，我们提供了读取标准CAD文件并返回TopoDS_Shape的读取器和写入器，或者在给出TopoDS_Shape的情况下写入CAD文件。OpenCASCADE命名空间中的大多数函数都处理一种或另一种类型的TopoDS_Shape，并提供普通deal.II对象的接口，如Triangulation、Manifold等。
 * 注意，这些工具中的大多数只有在spacedim等于3时才有用，因为OpenCASCADE只在三维模式下工作。在某些情况下，它们也可以在二维模式下使用，而第三维将被设置为零。
 * 如果你希望在空间的维度为二时使用这些工具，那么请确保你的CAD文件实际上是平的，并且所有的Z坐标都等于零，否则你会得到许多异常。
 *
 *
 */
namespace OpenCASCADE
{
  /**
   * 计算一个形状的子对象。这个函数对于收集作为参数传递的TopoDS_Shape的信息很有用。它返回包含在给定形状中的面、边和顶点（唯一与实际几何图形相关的拓扑实体）的数量。
   *
   */
  std::tuple<unsigned int, unsigned int, unsigned int>
  count_elements(const TopoDS_Shape &shape);

  /**
   * 读取IGES文件并将其内容翻译成openCascade拓扑实体。选项scale_factor用于补偿在IGES文件和目标应用程序中使用的不同单位。IGES文件的标准单位是毫英寸。返回的对象是一个TopoDS_Shape，包含文件中的所有对象。
   *
   */
  TopoDS_Shape
  read_IGES(const std::string &filename, const double scale_factor = 1e-3);

  /**
   * 将给定的拓扑形状写入一个IGES文件中。
   *
   */
  void
  write_IGES(const TopoDS_Shape &shape, const std::string &filename);


  /**
   * 读取STL文件并将其内容翻译成openCascade拓扑实体。
   * 返回的对象是一个TopoDS_Shape，其中包含文件中的所有对象。
   *
   */
  TopoDS_Shape
  read_STL(const std::string &filename);

  /**
   * 将给定的拓扑形状写入一个STL文件中。为了做到这一点，该形状必须包含一个网格结构，该函数检查该形状的所有面是否有一个附加的网格，如果不是这样，它将自动进行网格化。我们注意到，OpenCASCADE的自动网格生成只考虑到形状和网格之间的几何相似性，要控制三角形的形状和规则性，你应该使用其他网格软件。两个参数`deflection`和`angular_deflection`选择创建的三角形相对于原始拓扑形状的精度。参数`sew_different_faces`给出了使用OpenCASCADE的Sewer的可能性，使用参数`sewer_tolerance`创建一个水密的封闭STL。参数`is_relative`指定距离是否是相对的，`in_parallel`如果执行应该是并行的。
   *
   */
  void
  write_STL(const TopoDS_Shape &shape,
            const std::string & filename,
            const double        deflection,
            const bool          sew_different_faces = false,
            const double        sewer_tolerance     = 1e-6,
            const bool          is_relative         = false,
            const double        angular_deflection  = 0.5,
            const bool          in_parallel         = false);


  /**
   * 读取STEP文件并将其内容翻译成openCascade拓扑实体。选项scale_factor用于补偿在STEP文件和目标应用程序中使用的不同单位。STEP文件的标准单位是毫英寸。返回的对象是一个TopoDS_Shape，包含文件中的所有对象。
   *
   */
  TopoDS_Shape
  read_STEP(const std::string &filename, const double scale_factor = 1e-3);


  /**
   * 将给定的拓扑形状写入一个STEP文件中。
   *
   */
  void
  write_STEP(const TopoDS_Shape &shape, const std::string &filename);

  /**
   * 该函数返回与该形状相关的公差。每个CAD几何对象都被定义了一个公差，它表示其位置可能的不精确性。例如，顶点的公差表示它可以位于以名义位置为中心、半径为tol的球体中的任何一点。在执行一个操作时，例如将一个点投影到一个曲面上（这反过来也有它的公差），我们必须记住，投影的精度将受到建立曲面的公差的限制。
   * 公差的计算是以组成形状的子形状中的最大公差为依据。
   *
   */
  double
  get_shape_tolerance(const TopoDS_Shape &shape);

  /**
   * 执行给定的拓扑形状与平面 $c_x x + c_y y + c_z z +c = 0$
   * 的交点。返回的拓扑形状将包含尽可能少的bsplines。如果相交产生一个空的形状，会产生一个异常。
   *
   */
  TopoDS_Shape
  intersect_plane(const TopoDS_Shape &in_shape,
                  const double        c_x,
                  const double        c_y,
                  const double        c_z,
                  const double        c,
                  const double        tolerance = 1e-7);

  /**
   * 尝试将给定的TopoDS_Shape中包含的所有边连接成一个TopoDS_Edge，包含尽可能少的BSP线。如果输入的形状包含面，它们将被这个函数忽略。如果包含的边不能被连接成一个单一的边，也就是说，它们形成了不相连的曲线，那么将抛出一个异常。
   *
   */
  TopoDS_Edge
  join_edges(const TopoDS_Shape &in_shape, const double tolerance = 1e-7);

  /**
   * 创建一条通过指定向量中的点的平滑BSpline曲线，并将其存储在返回的TopoDS_Shape中（TopoDS_Edge类型）。如果方向与零不同，这些点将根据它们与方向的标量乘积在内部重新排序，否则，它们将被作为传递的点使用。注意，如果算法需要，这个函数会改变输入的点。
   * 该类用于插值一条通过数组点的BsplineCurve，具有C2连续性。如果可选的参数
   * @p closed 被设置为
   * "true"，那么除了第一个点之外，曲线在所有的点上都是C2的（这里只给出C1的连续性），并且它将是一条封闭的曲线。
   * 该曲线保证与输入点的距离为 @p tolerance
   * 。如果算法不能生成这样的曲线，就会抛出一个异常。
   *
   */
  template <int dim>
  TopoDS_Edge
  interpolation_curve(std::vector<Point<dim>> &curve_points,
                      const Tensor<1, dim> &   direction = Tensor<1, dim>(),
                      const bool               closed    = false,
                      const double             tolerance = 1e-7);

  /**
   * 从TopoDS_Shape中提取所有子形状，并将结果存储到标准容器中。如果该形状不包含某种类型的形状，相应的容器将是空的。
   *
   */
  void
  extract_geometrical_shapes(const TopoDS_Shape &        shape,
                             std::vector<TopoDS_Face> &  faces,
                             std::vector<TopoDS_Edge> &  edges,
                             std::vector<TopoDS_Vertex> &vertices);

  /**
   * 从一个单一的面创建一个三角结构。这个类提取构成这个面的参数化曲面的第一个u和v参数，并创建一个三角形<2,spacedim>，包含反映这个面的单个粗略的单元。如果这个面不是一个修剪过的面，这个单元的顶点将与原始TopoDS_Face的TopoDS_Vertex顶点重合。然而，情况往往不是这样的，用户应该注意如何使用这个网格。
   * 如果你用Triangulation<2,2>调用这个函数，请确保输入的面的所有Z坐标都设置为零，否则你会得到一个异常。
   *
   */
  template <int spacedim>
  void
  create_triangulation(const TopoDS_Face &         face,
                       Triangulation<2, spacedim> &tria);


  /**
   * 给定一个Triangulation和一个可选的Mapping，创建一个平滑曲线的向量，对Triangulation的边界顶点的连接部分进行插值，并作为TopoDS_Edge对象的向量返回。
   * 这个函数构造封闭的Bspline曲线对象，通过三角形边界的所有顶点，每个顶点都有
   * $C^2$  连续性，除了第一个顶点，那里只保证 $C^1$
   * 连续性。
   * 返回的曲线是按照构成三角形边界的面的索引排序的，即第一条曲线是从索引最低的面开始提取的，以此类推。
   * @param[in]  triangulation 输入三角形  @param[in]  mapping
   * 可选的输入映射  @return  TopoDS_Edge对象的一个 std::vector
   * ，代表 "三角形 "边界的平滑插值。
   *
   */
  template <int spacedim>
  std::vector<TopoDS_Edge>
  create_curves_from_triangulation_boundary(
    const Triangulation<2, spacedim> &triangulation,
    const Mapping<2, spacedim> &      mapping =
      StaticMappingQ1<2, spacedim>::mapping);

  /**
   * 从TopoDS_Shape中提取所有复合形状，并将结果存储到标准容器中。如果形状不包含某种类型的化合物，相应的容器将是空的。
   *
   */
  void
  extract_compound_shapes(const TopoDS_Shape &           shape,
                          std::vector<TopoDS_Compound> & compounds,
                          std::vector<TopoDS_CompSolid> &compsolids,
                          std::vector<TopoDS_Solid> &    solids,
                          std::vector<TopoDS_Shell> &    shells,
                          std::vector<TopoDS_Wire> &     wires);

  /**
   * 将点 @p origin 投射到由 @p
   * in_shape给出的拓扑形状上，并返回投射的点，包含该点的子形状以及该点在结果形状中的参数u和v坐标。如果形状不是基本的，它的所有子形状都会被迭代，首先是面，然后是边，返回的形状是离点最近的一个
   * @p origin.
   * 如果返回的形状是一个边，那么只有u坐标被填上合理的信息，而v坐标被设置为零。
   * 这个函数返回一个包含投影点、形状、u坐标和v坐标的元组（只有当得到的形状是一个面时，v坐标才会与零不同）。
   *
   */
  template <int dim>
  std::tuple<Point<dim>, TopoDS_Shape, double, double>
  project_point_and_pull_back(const TopoDS_Shape &in_shape,
                              const Point<dim> &  origin,
                              const double        tolerance = 1e-7);

  /**
   * 返回点 @p origin 在 @p in_shape.
   * 给出的拓扑形状上的投影。如果该形状不是基本形状，它的所有子形状都会被迭代，首先是面，然后是边，返回的点是与
   * @p in_shape, 最接近的点，无论其类型如何。
   *
   */
  template <int dim>
  Point<dim>
  closest_point(const TopoDS_Shape &in_shape,
                const Point<dim> &  origin,
                const double        tolerance = 1e-7);

  /**
   * 给出一个基本形状 @p in_shape
   * 和形状中的参考坐标，返回实空间中的相应点。如果该形状是TopoDS_Edge，
   * @p v
   * 坐标被忽略。只有由函数project_point_and_pull_back()返回的边或面可以作为该函数的输入。如果不是这样，就会抛出一个异常。
   *
   */
  template <int dim>
  Point<dim>
  push_forward(const TopoDS_Shape &in_shape, const double u, const double v);


  /**
   * 给出一个TopoDS_Face  @p face
   * 和这个面内的参考坐标，返回实空间中的对应点，该点的表面法线，以及最小和最大曲率，作为一个元组。
   *
   */
  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  push_forward_and_differential_forms(const TopoDS_Face &face,
                                      const double       u,
                                      const double       v,
                                      const double       tolerance = 1e-7);


  /**
   * 获取与给定的拓扑形状最接近的点，以及该点的法线和最小及最大曲率。如果该形状不是基本的，它的所有子面（只有面）都会被迭代，先是面，然后只返回最近的点。如果
   * @p in_shape
   * 不包含至少一个面，这个函数将抛出一个异常。
   *
   */
  std::tuple<Point<3>, Tensor<1, 3>, double, double>
  closest_point_and_differential_forms(const TopoDS_Shape &in_shape,
                                       const Point<3> &    origin,
                                       const double        tolerance = 1e-7);


  /**
   * 沿着 @p 方向与给定的 @p origin
   * 点和给定的拓扑形状相交的一条线。如果有一个以上的交点，它将返回最近的一个。
   * 可选的 @p tolerance 参数是用来计算距离的。
   *
   */
  template <int dim>
  Point<dim>
  line_intersection(const TopoDS_Shape &  in_shape,
                    const Point<dim> &    origin,
                    const Tensor<1, dim> &direction,
                    const double          tolerance = 1e-7);


  /**
   * 将OpenCASCADE点转换成Point<spacedim>。
   * 容差参数用于检查OpenCASCADE点的非使用成分是否接近于零。如果不是这样，在调试模式下会抛出一个断言。
   *
   */
  template <int spacedim>
  Point<spacedim>
  point(const gp_Pnt &p, const double tolerance = 1e-10);


  /**
   * 将Point<3>转换成OpenCASCADE点。
   *
   */
  template <int spacedim>
  gp_Pnt
  point(const Point<spacedim> &p);


  /**
   * 根据两个点与方向的标量乘积对其进行排序。如果方向的规范是零，那么就使用lexicographical排序。可选的参数在比较对象时被用作相对公差。
   *
   */
  template <int dim>
  bool
  point_compare(const Point<dim> &    p1,
                const Point<dim> &    p2,
                const Tensor<1, dim> &direction = Tensor<1, dim>(),
                const double          tolerance = 1e-10);


  /**
   * 当作为参数指定的点不在给定的TopoDS_Shape的 @p tolerance
   * 之间时抛出异常。
   *
   */
  template <int dim>
  DeclException1(ExcPointNotOnManifold,
                 Point<dim>,
                 << "The point [ " << arg1 << " ] is not on the manifold.");

  /**
   * 当作为参数指定的点不能被投影到流形上时，会产生异常。
   *
   */
  template <int dim>
  DeclException1(ExcProjectionFailed,
                 Point<dim>,
                 << "Projection of point [ " << arg1 << " ] failed.");

  /**
   * 当内部的OpenCASCADE工具不能返回OK状态时被抛出。
   *
   */
  DeclException1(ExcOCCError,
                 IFSelect_ReturnStatus,
                 << "An OpenCASCADE routine failed with return status "
                 << arg1);

  /**
   * 试图对一个退化的边缘进行曲线操作。
   *
   */
  DeclException0(ExcEdgeIsDegenerate);

  /**
   * 试图对错误的形状类型进行操作。
   *
   */
  DeclException0(ExcUnsupportedShape);
} // namespace OpenCASCADE


DEAL_II_NAMESPACE_CLOSE

#  endif // DEAL_II_WITH_OPENCASCADE

#endif // dealii_occ_utilities_h
 /*----------------------------- occ_utilities.h -----------------------------*/ 


