//include/deal.II-translator/grid/grid_in_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_grid_in_h
#define dealii_grid_in_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>

#include <iostream>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int space_dim>
class Triangulation;
template <int dim>
struct CellData;
#endif

/**
 *
 *
 * @note
 * 由于deal.II只支持线、四边形和六面体网格，这个类中的函数只能读取完全由这些单元组成的网格。如果你绝对需要使用三角形或四面体的网格，那么你唯一的选择就是将网格转换成四边形和六面体。一个可以做到这一点的工具是tethex，可用<a
 * href="https://github.com/martemyev/tethex">here</a>。
 * 你读取的网格将形成 @p Triangulation
 * 对象的最粗层次。因此，它必须不包含悬空节点或其他形式的自适应细化，否则如果输入文件所代表的网格确实有这些节点，就会发生奇怪的事情。这是由于大多数网格描述格式并不存储单元间的邻接信息，所以网格读取函数必须重新生成这些信息。他们通过检查两个单元是否有一个共同的面来实现这一目的。如果在一个三角形中存在悬挂的节点，相邻的单元格就没有共同的（完整的）面，所以网格读取器得出结论，相邻的单元格沿这些面没有邻居，因此必须在边界上。实际上，域的内部裂缝就是这样被引入的。由于这种情况很难被发现（GridIn如何决定两个小单元的面与面重合的地方或者一个较大的单元实际上是一个与局部细化相关的悬空节点，还是确实是域中的一个裂缝？如果你的目标是保存并在以后再次读取经过自适应改进的三角图，那么这个类就不是你的解决方案；而是看一下PersistentTriangulation类。
 * 要读取网格数据，需要填充的三角图必须是空的。在调用这个类的函数时，输入文件可能只包含一维的线；二维的线和四边形；以及三维的线、四边形和六角形。所有其他的单元类型（例如，二维的三角形，三维的三角形或四面体）都会被拒绝。
 * 这里的 "维度
 * "指的是网格的维度；它可能被嵌入到一个更高的维度空间中，比如一个嵌入三维的球体二维表面上的网格，或者一个离散化三维线条的一维网格）。结果将是一个由输入文件中描述的单元组成的三角形，并尽可能正确地设置输入文件中描述的材料指标和边界指标。
 *
 *
 * @note
 * 你不能期望三角结构中的顶点和单元编号与输入文件中的一致。(根据我们对单元和顶点分别进行编号的事实，这一点已经很清楚了，而对于某些输入文件格式则不是这样的；有些格式也不要求连续编号，或者从零以外的索引开始编号)。
 *
 *  <h3>Supported input formats</h3>
 * 目前，支持以下输入格式。  <ul>   <li>   @p UCD  （非结构化单元数据）格式：这种格式用于网格输入以及数据输出。如果输入文件中有数据矢量，它们会被忽略，因为我们只对这一类的网格感兴趣。UCD格式要求顶点按以下顺序排列：在2d中
 *
 * @verbatim
 *    3-----2
 *    |     |
 *    |     |
 *    |     |
 *    0-----1
 * @endverbatim
 * 而在三维中
 *
 * @verbatim
 *       7-------6        7-------6
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   3   |       |    3-------2   |
 *   |   4-------5    |       |   5
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 * 注意，这种排序方式与deal.II的编号方式不同，见Triangulation类。
 * 关于UCD格式的确切描述可以在AVS
 * Explorer手册中找到（见http://www.avs.com）。  @p UCD
 * 格式可以通过read_ucd()函数读取。
 * <li>  <tt>DB网格</tt>格式：该格式被 @p BAMG
 * 网格生成器使用（见http://www-rocq.inria.fr/gamma/cdrom/www/bamg/eng.htm.
 * @p BAMG
 * 手册中关于该格式的文档非常不完整，所以我们实际上并没有解析输出的许多字段，因为我们不知道它们的含义，但是读取的数据足以按照网格生成器的意图建立起网格。这种格式可以通过read_dbmesh()函数来读取。
 * <li>   @p XDA
 * 格式：这是一个相当简单的格式，由MGF代码使用。我们没有确切的格式规范，但是阅读器可以读到几个例子文件。如果阅读器没有摸透你的文件，那么扩展它应该是相当简单的。
 * <li>  <tt>%Gmsh 1.0 mesh</tt>格式：这个格式被 @p %Gmsh
 * 网格生成器使用（见http://gmsh.info/）。 @p %Gmsh
 * 手册中的文档解释了如何生成与deal.II库兼容的网格（即四边形而非三角形）。为了使用这种格式，%Gmsh必须以旧的1.0格式来输出文件。这可以在输入文件中加入
 * "Mesh.MshFileVersion = 1 "一行来完成。
 * <li>  <tt>%Gmsh 2.0
 * mesh</tt>格式：这是上述格式的一个变种。read_msh()函数会自动判断一个输入文件是版本1还是版本2。
 * <li>  <tt>Tecplot</tt>格式：这个格式被 @p TECPLOT
 * 使用，通常作为不同应用程序之间数据交换的基础。注意，目前只支持ASCII格式，二进制数据不能被读取。
 * <li>  <tt>UNV</tt>格式：这种格式是由Salome网格生成器生成的，见http://www.salome-platform.org/ 。这里记录了 GridIn::read_unv 函数所支持的格式的部分。  <ul>   <li>  第2411节：http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2411  <li>  第2412节：http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2412  <li>  第2467节：http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2467  <li>  这个格式的所有部分，即使在我们的阅读器中可能不被支持，也可以在这里找到：http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse  </ul>  注意，Salome，比方在2D，只能在一个正好有4条边（或4块边界）的物体上制作一个四维网格。这意味着，如果你有一个更复杂的物体，并想用四边形网格，你将需要把这个物体分解成>=2个独立的对象。然后1）对这些独立的对象进行网格化，2）创建与这些独立对象相关的适当的单元组和/或面，3）建立一个复合网格，4）删除所有可能与该复合网格的一些内部面相关的数字。
 * <li>
 * <tt>VTK</tt>格式。VTK非结构化网格遗留文件阅读器生成器。该阅读器目前只能处理非结构化网格格式的二维和三维几何数据。一般遗留的vtk文件，包括非结构化网格格式的文档可以在这里找到：http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
 * VTK格式要求顶点按以下顺序排列：在2d中
 *
 * @verbatim
 *    3-----2
 *    |     |
 *    |     |
 *    |     |
 *    0-----1
 * @endverbatim
 * 在三维中
 *
 * @verbatim
 *       7-------6        7-------6
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   4   |       |    4-------5   |
 *   |   3-------2    |       |   2
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 *
 * <li>  <tt>Cubit</tt>格式:
 * deal.II目前不直接支持从Cubit导入。然而，Cubit可以用一个简单的插件以UCD格式导出，然后产生的UCD文件可以被这个类读取。插件的脚本可以在deal.II
 * wiki页面的<a
 * href="https://github.com/dealii/dealii/wiki/Mesh-Input-And-Output">Mesh
 * Input and Output</a>下找到。
 * 另外，Cubit可以生成ABAQUS文件，可以通过read_abaqus()函数读入。这可能是一个更好的选择，对于具有复杂的边界条件表面和多种材料的几何形状。
 *
 * - 目前不容易通过Cubit的python接口获得的信息。
 * </ul>
 * <h3>Structure of input grid data. The GridReordering class</h3>
 * 你有责任在单元格列表中使用正确的顶点编号，即对于1d的线，你必须首先给出坐标值较低的顶点，然后是坐标值较高的顶点。对于二维的四边形，
 * @p quad 列表中的顶点索引必须使顶点按逆时针方向编号。
 * 在二维空间中，出现了另一个困难，这与四边形的意义有关。一个四边形由四条线组成，这些线有一个方向，根据定义，这个方向如下。
 *
 * @verbatim
 * 3-->--2
 * |     |
 * ^     ^
 * |     |
 * 0-->--1
 * @endverbatim
 * 现在，两个相邻的单元格必须有一个顶点编号，使得公共边的方向是相同的。例如，以下两个四边形
 *
 * @verbatim
 * 3---4---5
 * |   |   |
 * 0---1---2
 * @endverbatim
 * 可以用顶点编号<tt>(0 1 4 3)</tt>和<tt>(1 2 5
 * 4)</tt>来描述，因为从两个单元格看，中间的线会得到方向<tt>1->4</tt>。
 * 编号<tt>(0 1 4 3)</tt>和<tt>(5 4 1
 * 2)</tt>将是不允许的，因为左边的四边形将给公共线以<tt>1->4</tt>的方向，而右边的四边形则希望使用<tt>4->1</tt>，导致模糊不清。Triangulation对象能够检测到这种特殊情况，可以通过将右侧四边形的指数旋转2来消除这种情况。然而，如果你给了顶点指数<tt>(4
 * 1 2
 * 5)</tt>，它就不知道该怎么做了，因为这样它就必须旋转一个元素或三个元素，决定采取哪种方式还没有实现。
 * 还有更多模棱两可的情况，如果不使用复杂的算法，三角法可能根本就不知道该怎么做。此外，类似的问题也存在于三个空间维度中，其中面和线的方向需要被照顾到。
 * 由于这个原因，这个类的<tt>read_*</tt>函数在读取各种输入格式的网格时，会调用GridReordering类，将定义单元格的顶点顺序变成满足三角计算类要求的顺序。如果你在通过这个类读取网格时遇到意外问题，请务必阅读该类的文档。
 *
 *  <h3>Dealing with distorted mesh cells</h3>
 * 对于每一个网格读取函数，最后一个调用总是指向 Triangulation::create_triangulation(). ，该函数检查它所创建的作为粗网格一部分的所有单元是否扭曲（这里的扭曲意味着从参考单元到实际单元的映射的Jacobian有一个非正的行列式，即该单元被挤压或扭曲；参见术语表中 @ref GlossDistorted "扭曲的单元 "
 * 条目）。如果它发现任何这样的单元格，它就会抛出一个异常。这个异常不会在当前类的网格阅读器函数中被捕获，因此会传播到调用它的函数中。在那里，如果你确定处理这样的单元格没有坏处的话，你可以捕捉并忽略这个异常。如果你不知道你的网格有这样的单元格，那么如果你忽略这个异常，你的结果可能最多只能是可疑的质量。
 *
 *
 *
 * @ingroup grid
 * @ingroup input  Pelteret 2015, Timo Heister 2015, Krzysztof Bzowski, 2015
 *
 */

template <int dim, int spacedim = dim>
class GridIn
{
public:
  /**
   * 可能的网格输入格式列表。这些值在调用函数read()时使用，以确定要调用的实际阅读器。
   *
   */
  enum Format
  {
    /// Use GridIn::default_format stored in this object
    Default,
    /// Use read_unv()
    unv,
    /// Use read_ucd()
    ucd,
    /// Use read_abaqus()
    abaqus,
    /// Use read_dbmesh()
    dbmesh,
    /// Use read_xda()
    xda,
    /// Use read_msh()
    msh,
    /// Use read_tecplot()
    tecplot,
    /// Use read_vtk()
    vtk,
    /// Use read_vtu()
    vtu,
    /// Use read_assimp()
    assimp,
    /// Use read_exodusii()
    exodusii,
  };

  /**
   * 构造函数。
   *
   */
  GridIn();

  /**
   * 构造函数。附上这个三角图，用网格数据来输入。
   *
   */
  GridIn(Triangulation<dim, spacedim> &tria);

  /**
   * 附上这个三角图，用网格数据来输入。
   *
   */
  void
  attach_triangulation(Triangulation<dim, spacedim> &tria);

  /**
   * 从给定的流中读取。如果没有给出格式，则使用
   * GridIn::Format::Default 。
   *
   */
  void
  read(std::istream &in, Format format = Default);

  /**
   * 打开由字符串给定的文件，并调用前面的函数read()。
   * 这个函数使用PathSearch机制来查找文件。使用的文件类是
   * <code>MESH</code>  。
   *
   */
  void
  read(const std::string &in, Format format = Default);

  /**
   * 从一个非结构化的vtk文件中读取网格数据。vtk文件可能包含以下VTK单元类型。VTK_HEXAHEDRON（12），VTK_TETRA（10），VTK_QUAD（9），VTK_TRIANGLE（5），以及VTK_LINE（3）。
   * 根据模板维度的不同，只接受上述的部分内容。
   * 特别是，在三个维度中，该函数希望文件包含
   *
   *
   *
   *
   * - VTK_HEXAHEDRON/VTK_TETRA细胞类型
   *
   * - VTK_QUAD/VTK_TRIANGLE单元类型，以指定可选的边界或内部四边形面
   *
   *
   *
   *
   *
   * - VTK_LINE单元格类型，用于指定可选的边界或内部边缘 在二维方面。
   *
   *
   *
   *
   *
   *
   * - VTK_QUAD/VTK_TRIANGLE单元类型
   *
   * - VTK_LINE单元格类型，用于指定可选的边界或内部边缘 在一个维度上
   *
   *
   *
   *
   *
   *
   * - VTK_LINE单元格类型 输入文件可以使用[VTK文件格式](http://www.vtk.org/VTK/img/file-formats.pdf)的CELL_DATA部分指定边界ID、材料ID和流形ID。    该函数解释了输入文件中包含的两种CELL_DATA类型。`SCALARS MaterialID`，用于指定单元的材料ID，或面和边的边界ID，以及`SCALARS ManifoldID`，可用于指定任何三角测量对象（单元、面或边）的流形ID。    配套的 GridOut::write_vtk 函数可用于编写与此方法兼容的VTK文件。
   * @ingroup simplex
   *
   */
  void
  read_vtk(std::istream &in);

  /**
   * 从非结构化的vtu文件中读取网格数据，由deal.II使用
   * GridOut::write_vtu(), 保存，标志
   * GridOutFlags::Vtu::serialize_triangulation 设置为真。
   * 注意，这个函数不支持读入任意的vtu文件，而只支持由deal.II自己编写的文件，使用函数
   * GridOut::write_vtu 并设置 GridOutFlags::Vtu::serialize_triangulation
   * 为true。    当这个标志被设置为 "true
   * "时，生成的vtu文件在一个xml部分中包含了三角测量，这个部分会被一般的vtu阅读器忽略。
   * 如果没有这个部分，就会产生一个异常。
   *
   */
  void
  read_vtu(std::istream &in);


  /**
   * 读取由Salome网格生成器生成的unv文件中的网格数据。数值数据被忽略。
   * 请注意在这个类的一般文档中对生成这种文件格式的评论。
   *
   */
  void
  read_unv(std::istream &in);

  /**
   * 从一个ucd文件中读取网格数据。数值数据被忽略。
   * 不可能使用ucd文件为同一个单元同时设置 boundary_id 和
   * manifold_id。但是可以使用标志apply_all_indicators_to_manifolds来决定文件中的指标是指流形（标志设置为真）还是指边界（标志设置为假）。如果该标志被设置，这些指标也会被用作单元格的流形标识。
   *
   */
  void
  read_ucd(std::istream &in,
           const bool    apply_all_indicators_to_manifolds = false);

  /**
   * 从Abaqus文件中读取网格数据。数值和构成数据被忽略。与ucd文件格式的情况一样，可以使用标志apply_all_indicators_to_manifolds来决定文件中的指标是指流形（标志设置为真）还是指边界（标志设置为假）。
   * @note
   * 目前这个网格阅读器的实现是次优的，因此对于大的网格来说可能会很慢。
   * @note  Cubit的使用提示。
   *
   *
   *
   *
   *
   *
   * - 可以在网格中定义多个材料ID。  这可以通过在预处理程序中指定块组来实现。
   *
   *
   *
   *
   * - 可以在网格中定义任意的表面边界。  这是通过在预处理程序中指定边集来实现的。特别是，边界不仅仅局限于表面（在三维中），单个元素面也可以被添加到边集中。当边界条件要应用在一个复杂的形状边界上时，这是非常有用的，因为仅仅使用 "面 "是难以定义的。类似的情况也可以在2d中完成。
   * @note  该文件格式的兼容性信息如下。
   *
   *
   *
   *
   *
   *
   * - 在Abaqus CAE 6.12中生成的文件已被验证可以正确导入，但较早（或较新）版本的Abaqus也可以生成有效的输入甲板。
   *
   *
   *
   *
   *
   *
   * - 使用Cubit 11.x, 12.x, 13.x, 14.x和15.x生成的文件是有效的,但只有在使用一组特定的导出步骤时。这些步骤如下。
   *
   *
   *
   *
   * - 点击右侧工具栏中的光盘图标，进入 "分析设置模式"。
   *
   *
   *
   *
   * - 在 "操作 "下选择 "导出网格"，点击右边工具栏上的必要图标。
   *
   * - 选择一个输出文件。在Cubit 11.0和12.0版本中，可能需要点击浏览按钮，在弹出的对话中输入。
   *
   *
   *
   *
   *
   *
   * - 选择要输出的尺寸。
   *
   *
   *
   *
   * - 勾选覆盖的方框。
   *
   *
   *
   *
   *
   *
   * - 如果使用Cubit v12.0以后,取消勾选 "使用Cubit ID的导出"。如果不勾选这个框，一个无效的文件将遇到错误。
   *
   *
   *
   *
   * - 点击应用。
   *
   */
  void
  read_abaqus(std::istream &in,
              const bool    apply_all_indicators_to_manifolds = false);

  /**
   * 从一个包含DB网格格式数据的文件中读取网格数据。
   *
   */
  void
  read_dbmesh(std::istream &in);

  /**
   * 从一个包含XDA格式数据的文件中读取网格数据。
   *
   */
  void
  read_xda(std::istream &in);

  /**
   * 从msh文件中读取网格数据，可以是该文件格式的第一版或第二版。%Gmsh格式的文件在http://www.gmsh.info/。
   * @note
   * deal.II的输入函数不区分换行和其他空白。因此，deal.II将能够读取比%Gmsh略微通用的格式的文件。
   * @ingroup simplex
   *
   */
  void
  read_msh(std::istream &in);

#ifdef DEAL_II_GMSH_WITH_API
  /**
   * 使用Gmsh
   * API读取网格数据。任何由Gmsh支持的文件都可以作为参数传递。格式是由文件名的扩展名推断出来的。
   * 该函数将非命名的物理id（gmsh格式<4.0）解释为材料或边界id（类似于其他read_msh()函数的情况）。如果你想指定非默认的流形或边界ID，你必须将所有需要非默认边界或流形ID的实体分组到命名的物理组中，其中名称是使用应用于
   * `std::map<std::string,  int>`的函数 Patterns::Tools::to_value()
   * 来解释的。键值可以是`MaterialID`（如果物理组指的是尺寸为`dim`的对象），`BoundaryID`（如果组指的是尺寸<`dim`的对象），或者`ManifoldID`。
   * 从Gmsh文档中，物理标签的格式遵循以下惯例。
   * @code
   * \$PhysicalNames // same as MSH version 2
   * numPhysicalNames(ASCII int)
   * dimension(ASCII int) physicalTag(ASCII int) "name"(127 characters max)
   * ...
   * \$EndPhysicalNames
   * @endcode
   * 例如，下面是网格文件的片段
   * @code
   * MeshFormat
   * 4.1 0 8
   * \$EndMeshFormat
   * \$PhysicalNames
   * 4
   * 1 1 "ManifoldID:0"
   * 1 2 "BoundaryID:
   *
   * -1, ManifoldID: 1"
   * 2 3 "ManifoldID: 1"
   * 2 4 "MaterialID: 2, ManifoldID: 1"
   * \$EndPhysicalNames
   * \$Entities
   * ...
   * @endcode
   * 指的是一个二维网格，其中。
   *
   *
   *
   *
   *
   *
   * - 维度1的部分边界具有物理标签1，流形ID 0
   *
   *
   *
   *
   *
   *
   * - 一些内部面（维度为1的线）有流形ID 1
   *
   *
   *
   *
   *
   * - 一些元素的流形ID为1（而材料ID等于默认值，即0）。
   *
   *
   *
   *
   *
   *
   * - 一些元素的流形id为1，材料id等于2 如果物理组没有命名，那么行为与其他read_msh()函数相同，即物理标签本身被解释为边界或材料id。
   * @ingroup simplex
   *
   */
  void
  read_msh(const std::string &filename);
#endif

  /**
   * 从一个包含tecplot
   * ASCII数据的文件中读取网格数据。这在没有安装任何tecplot的情况下也可以工作。
   *
   */
  void
  read_tecplot(std::istream &in);

  /**
   * 读取Assimp支持的文件，并从中生成一个三角图。
   * 如果你指定一个 @p mesh_index,
   * ，只有给定索引的网格会被提取出来，否则文件中存在的所有网格都会被用于生成三角图。
   * 这个函数只能用来读取二维网格（可能嵌入三维）。这是图形软件的标准，如blender或3D
   * studio max，这也是原始Assimp库的建立目的。我们把它
   * "弯曲
   * "到deal.II，以支持复杂的共维一维网格和复杂的二维网格。
   * 如果 @p remove_duplicates
   * 被设置为true（默认），那么重复的顶点将被移除，如果它们的距离低于
   * @p tol.
   * ，那么只有与给定维度和空间维度兼容的元素将被从网格中提取，并且只支持那些与deal.II兼容的元素。如果你设置了`ignore_unsupported_element_types`，所有其他的元素类型都会被这个算法直接忽略。例如，如果你的网格包含三角形和四边形的混合体，那么只有四边形会被提取出来。如果你混合了兼容和不兼容的元素类型，所得到的网格（在三角剖分对象中表示）可能没有任何意义。如果`ignore_unsupported_element_types`被设置为`false`，那么当遇到不支持的类型时，就会抛出一个异常。
   * @param  filename 要读取的文件  @param  mesh_index
   * 文件中网格的索引  @param  remove_duplicates 删除重复的顶点
   * @param  tol 删除顶点时使用的公差  @param
   * ignore_unsupported_element_types
   * 如果我们在解析过程中遇到不支持的类型，不要抛出异常。
   *
   */
  void
  read_assimp(const std::string &filename,
              const unsigned int mesh_index = numbers::invalid_unsigned_int,
              const bool         remove_duplicates                = true,
              const double       tol                              = 1e-12,
              const bool         ignore_unsupported_element_types = true);

  /**
   * 一个包含ExodusII提供的一些信息的结构，这些信息在Triangulation对象中没有直接表示。
   * @note
   * 这个结构的存在是为了与可能提供额外输出数据的read_exodusii的未来版本实现向前兼容，但现在它只有一个字段。
   *
   */
  struct ExodusIIData
  {
    /**
     * 一个向量，包含从deal.II边界id（或流形id）到所提供的ExodusII
     * sideset id的映射。
     *
     */
    std::vector<std::vector<int>> id_to_sideset_ids;
  };

  /**
   * 读取以ExodusII文件格式存储的网格。    ExodusII是一种功能丰富的文件格式，比本类支持的大多数其他网格格式支持更多的功能（如节点集、有限元场、质量保证数据等）。其中许多特征在deal.II中没有等效的表示，因此不被支持（例如，deal.II不直接将自由度分配给节点，所以以节点格式存储的数据不能被这个函数加载）。目前，只从输入文件中提取以下信息。      <ol>   <li>  块标识：元素的块标识被加载为其材料标识。 </li>   <li>  元素和顶点：存储在ExodusII文件中的核心几何信息填充到所附的三角测量对象中。高阶元素被自动截断为低阶元素，因为deal.II不支持这一功能（例如，在deal.II中没有相当于 <code>QUAD9</code> 的元素，因为所有四边形都有四个顶点，额外的几何信息被存储在Manifold或类似MappingQEulerian的东西）。 </li>   <li>  边集id：这些被解释为边界id或流形id（见下面关于输出值的说明）。如果你试图读取一个ExodusII文件，将一个sidet id分配给一个内部面的边界id，将会发生错误。 </li>   </ol>  对于非零二维的三角形，边集id不被翻译，因为这些三角形不支持边界id的设置。      @param  filename 要读取的文件的名称。      @param  apply_all_indicators_to_manifolds 布尔型，决定边集id是否应该被解释为流形id或边界id。默认值是<tt>false</tt>，即把所有sidet ids都当作边界ids。如果你的网格将sidet ids设置在内部面，那么就有必要将这个参数设置为 <code>true</code> ，然后做一些后处理来正确设置边界ids。      @return  这个函数返回一个结构体，包含一些由ExodusII文件存储的额外数据，这些数据不能被加载到Triangulation中。
   *
   * - 更多信息见ExodusIIData。    ExodusII中的一个单元面可以在任意数量的边集中（即，它可以有任意数量的边集ID
   *
   * - 然而，在deal.II中，一个边界单元面正好有一个边界id。所有不在边集中的边界面都被赋予（默认）边界ID  $0$  。这个函数然后将边集组合成唯一的集合，并给每个集合一个边界ID。比如说。考虑一个单四边形网格，其左侧没有边集ID，右侧有边集ID  $0$  和  $1$  ，其底部和顶部的边集ID为  $0$  。左侧面的边界id为 $0$ ，顶部和底部面的边界id为 $1$ ，右侧面的边界id为 $2$  。因此，在这种情况下，这个函数返回的向量将是  $\{\{\}, \{0\}, \{0, 1\}\}$  。
   *
   */
  ExodusIIData
  read_exodusii(const std::string &filename,
                const bool         apply_all_indicators_to_manifolds = false);

  /**
   * 以这种格式返回文件的标准后缀。
   *
   */
  static std::string
  default_suffix(const Format format);

  /**
   * 返回格式名称的枚举格式。
   *
   */
  static Format
  parse_format(const std::string &format_name);

  /**
   * 返回一个已实现的输入格式的列表。不同的名字用竖条符号（<tt>`|'</tt>）分开，如ParameterHandler类所使用的那样。
   *
   */
  static std::string
  get_format_names();

  /**
   * 异常情况
   *
   */
  DeclException1(ExcUnknownSectionType,
                 int,
                 << "The section type <" << arg1 << "> in an UNV "
                 << "input file is not implemented.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcUnknownElementType,
                 int,
                 << "The element type <" << arg1 << "> in an UNV "
                 << "input file is not implemented.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcUnknownIdentifier,
                 std::string,
                 << "The identifier <" << arg1 << "> as name of a "
                 << "part in an UCD input file is unknown or the "
                 << "respective input routine is not implemented."
                 << "(Maybe the space dimension of triangulation and "
                 << "input file do not match?");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNoTriangulationSelected,
                   "No Triangulation has been attached to this GridIn object "
                   "so that nothing can be filled during any read function "
                   "calls.  Please pass a reference to the Triangulation tria "
                   "to be  filled in the constructor GridIn(tria) or attach "
                   "it with the function call GridIn::attach_triangulation().");
  /**
   * 异常情况
   *
   */
  DeclException2(
    ExcInvalidVertexIndex,
    int,
    int,
    << "While creating cell " << arg1
    << ", you are referencing a vertex with index " << arg2
    << " but no vertex with this index has been described in the input file.");
  /**
   * 异常情况
   *
   */
  DeclException3(
    ExcInvalidVertexIndexGmsh,
    int,
    int,
    int,
    << "While creating cell " << arg1 << " (which is numbered as " << arg2
    << " in the input file), you are referencing a vertex with index " << arg3
    << " but no vertex with this index has been described in the input file.");
  /**
   * 异常情况
   *
   */
  DeclException0(ExcInvalidDBMeshFormat);
  /**
   * 异常情况
   *
   */
  DeclException1(ExcInvalidDBMESHInput,
                 std::string,
                 << "The string <" << arg1
                 << "> is not recognized at the present"
                 << " position of a DB Mesh file.");

  /**
   * 异常情况
   *
   */
  DeclException1(
    ExcDBMESHWrongDimension,
    int,
    << "The specified dimension " << arg1
    << " is not the same as that of the triangulation to be created.");

  DeclException1(ExcInvalidGMSHInput,
                 std::string,
                 << "The string <" << arg1
                 << "> is not recognized at the present"
                 << " position of a Gmsh Mesh file.");

  DeclException1(ExcGmshUnsupportedGeometry,
                 int,
                 << "The Element Identifier <" << arg1 << "> is not "
                 << "supported in the deal.II library when "
                 << "reading meshes in " << dim << " dimensions.\n"
                 << "Supported elements are: \n"
                 << "ELM-TYPE\n"
                 << "1 Line (2 nodes, 1 edge).\n"
                 << "3 Quadrilateral (4 nodes, 4 edges).\n"
                 << "5 Hexahedron (8 nodes, 12 edges, 6 faces) when in 3d.\n"
                 << "15 Point (1 node, ignored when read)");


  DeclException0(ExcGmshNoCellInformation);

protected:
  /**
   * 存储要用读入的数据输入的三角图的地址。
   *
   */
  SmartPointer<Triangulation<dim, spacedim>, GridIn<dim, spacedim>> tria;

  /**
   * 这个函数可以将<tt>read_*</tt>函数创建的原始单元数据对象以Gnuplot格式写入一个流。如果想看看实际创建的数据，这有时是很方便的，如果知道数据在某些方面不正确，但Triangulation类因为这些错误而拒绝生成三角图。特别是，这个类的输出写出了单元格的编号，以及每个单元格的面的方向。特别是需要后者的信息来验证单元格数据对象是否遵循单元格及其面的排序要求，即所有面需要有唯一的方向和相对于相邻单元格的指定方向（见该类和GridReordering类的文档）。
   * 这个函数的输出包括每条界于单元格的线的向量，表明它相对于这个单元格的方向，以及单元格的编号。整个输出的形式是，它可以被Gnuplot读入，并生成完整的图，而不需要用户进一步操作。
   *
   */
  static void
  debug_output_grid(const std::vector<CellData<dim>> &  cells,
                    const std::vector<Point<spacedim>> &vertices,
                    std::ostream &                      out);

private:
  /**
   * 跳过输入流中的空行，即不包含任何内容或只包含白色空间的行。
   *
   */
  static void
  skip_empty_lines(std::istream &in);

  /**
   * 跳过给定输入流目前所在点之后以指定字符（例如<tt>#</tt>）开始的注释行。调用此函数后，数据流处于注释行后第一行的开始位置，如果没有注释行，则处于与之前相同的位置。
   *
   */
  static void
  skip_comment_lines(std::istream &in, const char comment_start);

  /**
   * 这个函数做了一个讨厌的工作（由于非常宽松的约定和不同版本的tecplot格式），从一个tecplot头中提取重要的参数，包含在字符串
   * @p header.
   * 中，其他变量是输出变量，它们的值对函数执行没有影响。
   *
   */
  static void
  parse_tecplot_header(std::string &              header,
                       std::vector<unsigned int> &tecplot2deal,
                       unsigned int &             n_vars,
                       unsigned int &             n_vertices,
                       unsigned int &             n_cells,
                       std::vector<unsigned int> &IJK,
                       bool &                     structured,
                       bool &                     blocked);

  /**
   * 如果没有给出格式，read()使用的输入格式。
   *
   */
  Format default_format;
};

 /* -------------- declaration of explicit specializations ------------- */ 

#ifndef DOXYGEN

template <>
void
GridIn<2>::debug_output_grid(const std::vector<CellData<2>> &cells,
                             const std::vector<Point<2>> &   vertices,
                             std::ostream &                  out);


template <>
void
GridIn<2, 3>::debug_output_grid(const std::vector<CellData<2>> &cells,
                                const std::vector<Point<3>> &   vertices,
                                std::ostream &                  out);
template <>
void
GridIn<3>::debug_output_grid(const std::vector<CellData<3>> &cells,
                             const std::vector<Point<3>> &   vertices,
                             std::ostream &                  out);
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


