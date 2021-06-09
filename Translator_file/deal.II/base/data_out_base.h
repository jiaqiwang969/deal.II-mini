//include/deal.II-translator/base/data_out_base_0.txt
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

#ifndef dealii_data_out_base_h
#define dealii_data_out_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/numerics/data_component_interpretation.h>

// To be able to serialize XDMFEntry
#include <boost/serialization/map.hpp>

#include <limits>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

// Only include the Tecplot API header if the appropriate files
// were detected by configure
#ifdef DEAL_II_HAVE_TECPLOT
#  include <string.h>

#  include "TECIO.h"
#endif

#include <ostream>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class ParameterHandler;
class XDMFEntry;
#endif

/**
 * 这是一个基类，用于输出非常一般形式的网格上的数据。输出数据被期望为一组<tt>patches</tt>，并以可视化工具所期望的格式写入输出流。对于输出格式的列表，请查看枚举#OutputFormat。对于其中列出的每一种格式，这个类包含一个函数<tt>write_format</tt>，写出输出。关于某种格式的细节，请参考这些函数的文档。
 * <h3>Structure of the output data</h3>
 * 数据不是用deal.II网状结构写入的。相反，它依赖于一组由派生类（例如DataOut、DataOutStack、DataOutFaces、DataOutRotation或MatrixOut类）创建的<tt>patches</tt>。
 * 每个补丁描述一个网格的单一逻辑单元，可能会被细分若干次，以表示在这个单元上定义的高阶多项式。为此，一个补丁由一个<tt>dim</tt>维的规则网格组成，每个方向上的网格点数量相同。在最简单的情况下，它可能由单个网格单元的角点组成。对于这个局部网格的每个点，Patch包含任意数量的数据值，不过每个Patch上每个点的数据集数量必须相同。
 * 通过为不同的输出格式提供这个接口，可以简单地将这个类扩展到新的格式，而不需要依赖诸如实际的三角计算和数据矢量的处理。这些东西应该由派生类提供，它有一个用户可调用的接口，然后。
 * 在每个补丁中，数据是按照通常的词典顺序组织的，<i>x</i>运行最快，然后是<i>y</i>和<i>z</i>。节点是按照这个顺序存储的，单元也是如此。三维的每个单元都被存储，使其正面位于<i>xz</i>平面。为了提高这一概念的可理解性，以下两节是从本文档的前一版本中保留下来的。
 *
 *  <h4>Patches</h4>
 * 网格可以被认为是一个单元格的集合；如果你想在这样的网格上写出数据，你可以通过一次写一个单元格来实现。因此，这个类中的函数接收一个描述每个单元格上的数据的对象列表。每个单元的数据通常包括这个单元的顶点列表，以及每个顶点的数据值（例如，解决方案数据、错误信息等）的列表。
 * 然而，在某些情况下，单元格的这种接口过于局限。例如，你可能有高阶元素，只打印顶点的值是不够的。出于这个原因，我们不仅提供了只写顶点上的数据，而且还将数据组织成每个单元的张量积网格。参数<tt>n_subdivisions</tt>是为每个补丁单独给出的，它表示单元被分割输出的频率；例如，<tt>n_subdivisions==1</tt>产生的单元没有被分割，<tt>n_subdivisions==2</tt>将产生一个二维空间3乘3点的网格，三维空间3乘3点，<tt>n_subdivisions==3</tt>将产生4乘4（乘4）点，等等。这些点在补丁上的实际位置将通过多线性变换从为这个补丁给出的顶点计算出来。
 * 对于边界上的单元，可能会使用一个映射来计算内部点的位置。在这种情况下，坐标被存储在补丁内，因为它们不容易被恢复。
 * 鉴于这些注释，要在这个点的补丁上打印的实际数据由几个数据集组成，每个数据集在每个补丁点上都有一个值。例如在两个空间维度的<tt>n_subdivisions==2</tt>，每个数据集要提供9个值，由于补丁要被打印成张量积（或其向实空间单元的转化），其值要像<i>(x0,y0)
 * (x0,y1) (x0,y2) (x1,y0) (x1,y1) (x1,y2) (x2,y0) (x2,y1)
 * (x2,y2)</i>那样排序，即z坐标运行最快，然后是y坐标，然后是x（如果有那么多空间方向）。
 *
 *  <h4>Generalized patches</h4>
 * 一般来说，上面解释的补丁可能太受限制。例如，在一个三维计算中，如果人们对内部发生的事情不感兴趣，可能只想画出一个域的外表面。那么，在一个三维的世界中，应该绘制的对象是二维的。Patch类和相关的输出函数可以处理这些情况。因此，Patch类需要两个模板参数，第一个名为<tt>dim</tt>，表示对象的维度（在上面的例子中，这将是两个），而第二个名为<tt>spacedim</tt>，表示嵌入空间的维度（这将是三个）。一个补丁的角点具有空间的维度，而它们的数量则由补丁的维度决定。默认情况下，第二个模板参数的值与第一个相同，这将对应于输出一个单元，而不是一个面或其他东西。
 * <h3>DataOutBaseInterface</h3>
 * 这个命名空间的成员通常不会从用户代码中直接调用。相反，使用这里声明的函数的类通常是从DataOutInterface派生的。
 * 这个类的接口基本上由一个描述补丁的数据类型的声明和一堆函数组成，这些函数接收一个补丁列表，并将它们以一种格式或其他方式写入流。派生类的责任是提供这个补丁列表。除了补丁列表之外，还可以给每个数据集起一个名字。
 *
 *  <h3>Querying interface</h3>
 * 这个类还提供了几个函数（parse_output_format(),
 * get_output_format_names(),
 * default_suffix()），可以用来查询这个类支持哪些输出格式。这些函数提供了一个我们可以输出的所有格式的名称列表，解析一个字符串并返回一个表示每个格式的枚举，并提供一种方法将这个枚举的值转换为该名称的文件通常使用的后缀。使用这些函数，可以使应用程序完全摆脱对库目前允许输出的格式的了解；几个例子程序显示了如何做到这一点。
 * <h3>Output parameters</h3>
 * 所有的函数都有一个参数，是一个<tt>XFlags</tt>类型的结构，其中<tt>X</tt>是输出格式的名称。要知道目前支持哪些标志，请阅读不同结构的文档。
 * 注意，通常用于科学可视化程序的输出格式没有或只有很少的参数（除了一些兼容性标志），因为在那里，输出的实际外观是由可视化程序决定的，这个类产生的文件或多或少只存储原始数据。
 * 直接的输出格式，如Postscript或Povray，需要给予更多的参数，虽然，因为在那里，输出文件必须包含所有的观点、光源等细节。
 * <h3>Writing backends</h3>
 * 引入了一个抽象层，以方便其他可视化工具的编码后端。它适用于将信息分离成顶点字段、网格单元的连接信息字段和数据字段的数据格式。
 * 对于每一个字段，都实现了输出函数，即write_nodes()、write_cells()和write_data()。为了使用这些函数，必须编写一个特定格式的输出流，遵循DXStream、GmvStream、VtkStream等的例子，在.cc文件中实现。
 * 在这个框架中，一个新的输出格式的实现被简化为编写章节头和新的输出流类，用于编写单个网格对象。
 * <h3>Credits</h3>  <ul>
 * <li>  EPS输出基于Stefan Nauber对旧的DataOut类的早期实现
 * <li>  Thomas Richter的Povray输出。
 * <li>  Benjamin Shelton Kirk的Tecplot输出。
 * <li>  拉格朗日VTK输出 作者：Alexander Grayver
 * </ul>
 *
 *
 */
namespace DataOutBase
{
  /**
   * 描述一个<tt>dim</tt>空间维度的数据补丁的数据结构。    一个补丁由以下数据组成。    <ul>   <li>  角的#顶点， <li>  Patch在每个空间方向的单元数的#n_细分， <li>  附在每个顶点的#数据，以通常的lexicographic排序， <li>  关于#邻居的信息。    </ul>  有关其内容和目的的更多信息，请参见DataOutBase类的一般文档。 在二维的情况下，下图是一个<tt>n_subdivisions</tt>=4的例子，因为每个补丁内的（子）单元数等于<tt>2<sup>dim</sup></tt>。
   * @ingroup output
   *
   */
  template <int dim, int spacedim = dim>
  struct Patch
  {
    /**
     * 使<tt>spacedim</tt>模板参数可用。
     *
     */
    static const unsigned int space_dim = spacedim;

    /**
     * 一个补丁的角点。 如果
     * <code>points_are_available==false</code>
     * ，内部点是通过单元格到这些角点所指定的单元格的多线性变换来计算的。
     * 另一方面，如果 <code>points_are_available==true</code>
     * ，那么要产生输出的点的坐标被附加到 <code>data</code>
     * 表中的额外行中。
     * 点的顺序与三角测量中的单元相同。
     *
     */
    Point<spacedim> vertices[GeometryInfo<dim>::vertices_per_cell];

    /**
     * 当前补丁的邻居的补丁索引。这是为OpenDX格式提供的，该格式需要邻居信息来进行高级输出。
     *
     */
    std::array<unsigned int, GeometryInfo<dim>::faces_per_cell> neighbors;

    /**
     * 这个补丁的编号。由于我们不确定补丁是否总是以相同的顺序处理，我们最好存储这个。
     *
     */
    unsigned int patch_index;

    /**
     * 这个补丁要写的分区的数量。
     * <tt>1</tt>表示没有细分，<tt>2</tt>表示二分法，<tt>3</tt>三分法，等等。
     *
     */
    unsigned int n_subdivisions;

    /**
     * 数据向量。其格式如下。<tt>data(i,.)</tt>表示属于<tt>i</tt>个数据向量的数据。<tt>data.n_cols()</tt>因此等于输出点的数量；这个数字是<tt>（subdivisions+1）^{dim}</tt>。<tt>data.n_rows()</tt>等于数据向量的数量。就目前而言，一个数据向量等于一个标量，即使多个标量后来可能被解释为向量。
     * 在每一列中，<tt>data(.,j)</tt>是输出点<tt>j</tt>的数据值，其中<tt>j</tt>表示deal.II中通常的lexicographic排序。这也是<tt>QIterated</tt>类提供的点的顺序，当与<tt>QTrapezoid</tt>类作为子正交使用时。
     * 由于所有要打印的补丁的数据向量的数量通常是相同的，<tt>data.size()</tt>应该对所有提供的补丁产生相同的值。例外情况是设置了point_are_available的补丁，其中点的实际坐标被附加到'data'字段中，见point_are_available标志的文档。
     *
     */
    Table<2, float> data;

    /**
     * 表示内部补丁点的坐标（假设该补丁应该被进一步细分）是否被附加到
     * @p data 表中（  @p true)  或不附加（  @p false).
     * 后者是默认的，在这种情况下，该补丁内部点的位置是由补丁顶点的（双，三）线性内插计算的。
     * 这个选项的存在是因为补丁点可能是用Mapping（而不是通过线性插值）来评估的，因此必须存储在Patch结构中。
     *
     */
    bool points_are_available;

    /**
     * 这个补丁的基础单元的参考单元类型。
     *
     */
    ReferenceCell reference_cell;

    /**
     * 默认构造函数。设置#n_subdivisions为1，#points_are_available为假，#patch_index为#no_neighbor。
     *
     */
    Patch();

    /**
     * 比较目前的补丁与另一个补丁是否相同。这在我们的测试套件中的一些自动测试中使用。
     *
     */
    bool
    operator==(const Patch &patch) const;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
     * <tt>std::map</tt>) ）是困难的。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 将当前对象的内容与给定参数的内容交换。
     *
     */
    void
    swap(Patch<dim, spacedim> &other_patch);

    /**
     * 如果这个补丁的一侧没有邻居，则使用该值。
     *
     */
    static const unsigned int no_neighbor = numbers::invalid_unsigned_int;

    /**
     * @addtogroup  Exceptions  
     * @{ 
     *
     */

    /**
     * 例外情况
     *
     */
    DeclException2(
      ExcInvalidCombinationOfDimensions,
      int,
      int,
      << "It is not possible to have a structural dimension of " << arg1
      << " to be larger than the space dimension of the surrounding"
      << " space " << arg2);
    //@}
  };



  /**
   * 一般Patch<dim,spacedim>模板的特殊化，专门针对点的情况，即嵌入
   * @p spacedim 维空间的零维对象。
   * 目前的类与通用模板兼容，允许使用相同的函数以通用方式访问任意维度的补丁。然而，它使一些对零维补丁来说毫无意义的变量变成了
   * @p static
   * 变量，这些变量在整个程序中只存在一次，而不是每个补丁一次。具体来说，
   * @p neighbors 数组和 @p n_subdivisions
   * 成员变量就是这种情况，它们对零维补丁没有意义，因为点在其不存在的面上没有自然相邻，也不能合理地被细分。
   *
   */
  template <int spacedim>
  struct Patch<0, spacedim>
  {
    /**
     * 使<tt>spacedim</tt>模板参数可用。
     *
     */
    static const unsigned int space_dim = spacedim;

    /**
     * 一个补丁的角点。
     * 对于当前的零维补丁类，当然只有一个顶点。
     * 如果  <code>points_are_available==true</code>
     * ，那么要产生输出的点的坐标将作为附加行附在
     * <code>data</code>  表上。
     *
     */
    Point<spacedim> vertices[1];

    /**
     * 一个未使用的、 @p static
     * 的变量，它的存在只是为了允许从一般的代码中以通用方式访问。
     *
     */
    static unsigned int neighbors[1];

    /**
     * 这个补丁的编号。因为我们不确定补丁是否总是以相同的顺序处理，所以我们最好存储这个。
     *
     */
    unsigned int patch_index;

    /**
     * 这个补丁要写的分区的数量。
     * <tt>1</tt>表示没有细分，<tt>2</tt>表示二分法，<tt>3</tt>三分法，等等。
     * 由于细分对零维补丁没有意义，这个变量不被使用，它的存在只是为了允许从一般代码中以通用方式访问。
     *
     */
    static unsigned int n_subdivisions;

    /**
     * 数据向量。其格式如下。<tt>data(i,.)</tt>表示属于<tt>i</tt>个数据向量的数据。<tt>data.n_cols()</tt>因此等于输出点的数量；鉴于我们在点上产生输出，这个数字对于当前的类来说当然是1。<tt>data.n_rows()</tt>等于数据向量的数量。对于当前的目的，一个数据向量等于一个标量，即使多个标量后来可能被解释为向量。
     * 在每一列中，<tt>data(.,j)</tt>是输出点<tt>j</tt>的数据值；对于当前类，
     * @p j 只能是0。
     * 由于所有要打印的补丁的数据向量的数量通常是相同的，<tt>data.size()</tt>应该对所有提供的补丁产生相同的值。例外的情况是设置了point_are_available的补丁，其中点的实际坐标被附加到'data'字段中，见point_are_available标志的文档。
     *
     */
    Table<2, float> data;

    /**
     * 表示内部补丁点的坐标（假设该补丁应该被进一步细分）是否被附加到
     * @p data 表中（  @p true)  或不附加（  @p false).
     * 后者是默认的，在这种情况下，该补丁内部点的位置是由补丁顶点的（双，三）线性内插计算的。
     * 这个选项的存在是因为补丁点可能是用Mapping（而不是通过线性插值）来评估的，因此必须存储在Patch结构中。
     *
     */
    bool points_are_available;

    /**
     * 这个补丁的基础单元的参考单元类型。
     *
     */
    ReferenceCell reference_cell;

    /**
     * 默认构造函数。设置#points_are_available为false，#patch_index为#no_neighbor。
     *
     */
    Patch();

    /**
     * 比较目前的补丁与另一个补丁是否平等。这在我们的测试套件中的一些自动测试中使用。
     *
     */
    bool
    operator==(const Patch &patch) const;

    /**
     * 返回这个对象的内存消耗估计值，以字节为单位。这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
     * <tt>std::map</tt>) ）是困难的。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 将当前对象的内容与给定参数的内容交换。
     *
     */
    void swap(Patch<0, spacedim> &other_patch);

    /**
     * 如果这个补丁的一侧没有邻居，则使用该值。
     *
     */
    static const unsigned int no_neighbor = numbers::invalid_unsigned_int;

    /**
     * @addtogroup  Exceptions  
     * @{ 
     *
     */

    /**
     * 例外情况
     *
     */
    DeclException2(
      ExcInvalidCombinationOfDimensions,
      int,
      int,
      << "It is not possible to have a structural dimension of " << arg1
      << " to be larger than the space dimension of the surrounding"
      << " space " << arg2);
    //@}
  };


  /**
   * 描述不同输出标志之间共同功能的基类。    这是用
   * "奇怪的重复模板模式
   * "实现的；派生类使用自己的类型来填充类型名，以便<tt>memory_consumption</tt>正确工作。更多信息请参见维基百科关于该模式的页面。
   * @ingroup output
   *
   */
  template <typename FlagsType>
  struct OutputFlagsBase
  {
    /**
     * 用这个类提供的名称和类型声明所有标志，以便在输入文件中使用。
     * 这个方法什么都不做，但是子类可以覆盖这个方法，将字段添加到<tt>prm</tt>。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取 declare_parameters()
     * 中声明的参数，并为这种输出格式设置相应的标志。
     * 这个方法什么都不做，但是子类可以重写这个方法来向<tt>prm</tt>添加字段。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);

    /**
     * 返回这个对象的内存消耗估计值，单位是字节。这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
     * <tt>std::map</tt>) 是很困难的。
     *
     */
    std::size_t
    memory_consumption() const;
  };


  template <typename FlagsType>
  void
  OutputFlagsBase<FlagsType>::declare_parameters(ParameterHandler &)
  {}


  template <typename FlagsType>
  void
  OutputFlagsBase<FlagsType>::parse_parameters(const ParameterHandler &)
  {}


  template <typename FlagsType>
  std::size_t
  OutputFlagsBase<FlagsType>::memory_consumption() const
  {
    return sizeof(FlagsType);
  }


  /**
   * 控制OpenDX格式的输出细节的标志。
   * @ingroup output
   *
   */
  struct DXFlags : public OutputFlagsBase<DXFlags>
  {
    /**
     * 写入邻居信息。例如，如果OpenDX要计算积分曲线（流线），这个信息是必要的。如果它不存在，流线就会在单元格边界结束。
     *
     */
    bool write_neighbors;
    /**
     * 以二进制格式写入三角图的整数值。
     *
     */
    bool int_binary;
    /**
     * 以二进制格式写入坐标向量。
     *
     */
    bool coordinates_binary;

    /**
     * 以二进制格式写入数据向量。
     *
     */
    bool data_binary;

    /**
     * 将二进制坐标向量写成双数（64位）而不是浮点数（32位）。
     *
     */
    bool data_double;

    /**
     * 构造函数。
     *
     */
    DXFlags(const bool write_neighbors    = false,
            const bool int_binary         = false,
            const bool coordinates_binary = false,
            const bool data_binary        = false);

    /**
     * 用本类提供的名称和类型声明所有标志，以便在输入文件中使用。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取 declare_parameters()
     * 中声明的参数，并相应地设置该输出格式的标志。
     * 这样得到的标志会覆盖这个对象以前所有的内容。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * 控制AVS的UCD格式的输出细节的标志。
   * @ingroup output
   *
   */
  struct UcdFlags : public OutputFlagsBase<UcdFlags>
  {
    /**
     * 在文件的开头写一个注释，说明创建日期和其他一些数据。
     * 虽然UCD格式和AVS支持这一点，但其他一些程序对此感到困惑，所以默认情况是不写序言。然而，可以用这个标志写一个序言。
     * 默认值。  <code>false</code>  .
     *
     */
    bool write_preamble;

    /**
     * 构造函数。
     *
     */
    UcdFlags(const bool write_preamble = false);

    /**
     * 用本类提供的名称和类型声明所有标志，以便在输入文件中使用。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取 declare_parameters()
     * 中声明的参数，并相应地设置该输出格式的标志。
     * 这样得到的标志会覆盖这个对象以前所有的内容。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * 控制Gnuplot格式的输出细节的标志。
   * @ingroup output
   *
   */
  struct GnuplotFlags : public OutputFlagsBase<GnuplotFlags>
  {
    /**
     * 默认构造函数。用<tt>"x"</tt>, <tt>"y"</tt>,
     * 和<tt>"z"</tt>的默认值设置尺寸标签。
     *
     */
    GnuplotFlags();

    /**
     * 构造函数，为尺寸标签设置非默认值。
     *
     */
    GnuplotFlags(const std::vector<std::string> &space_dimension_labels);

    /**
     * 在每个空间维度上使用的标签。这些标签默认为<tt>"x"</tt>,
     * <tt>"y"</tt>,
     * 和<tt>"z"</tt>。标签被打印到Gnuplot文件中，周围有角括号。例如，如果空间维度是2，标签是<tt>"x"</tt>和<tt>"t"</tt>，那么相关的行将以下列内容开始
     * @verbatim
     * # <x> <t>
     * @endverbatim
     * 任何额外的标签都会被忽略。
     * 如果你自己指定这些标签，那么至少应该有<tt>spacedim</tt>标签，其中<tt>spacedim</tt>是输出数据的空间尺寸。
     *
     */
    std::vector<std::string> space_dimension_labels;

    /**
     * 返回这个对象的内存消耗估计值，单位是字节。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 当没有足够的指定维度标签时引发的异常情况。
     *
     */
    DeclExceptionMsg(ExcNotEnoughSpaceDimensionLabels,
                     "There should be at least one space dimension per spatial "
                     "dimension (extras are ignored).");
  };

  /**
   * 控制Povray格式的输出细节的标志。有几个标志已经实现，请看它们各自的文档。
   * @ingroup output
   *
   */
  struct PovrayFlags : public OutputFlagsBase<PovrayFlags>
  {
    /**
     * 正常矢量插值，如果设置为真，默认=假
     *
     */
    bool smooth;

    /**
     * 使用二次方补丁（b-splines）而不是三角形。 默认 =
     * false
     *
     */
    bool bicubic_patch;

    /**
     * 包括外部的
     * "data.inc"，包括摄像机、灯光和场景的纹理定义。 默认
     * = false
     *
     */
    bool external_data;

    /**
     * 构造函数。
     *
     */
    PovrayFlags(const bool smooth        = false,
                const bool bicubic_patch = false,
                const bool external_data = false);

    /**
     * 用本类提供的名称和类型声明所有标志，以便在输入文件中使用。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取 declare_parameters()
     * 中声明的参数，并相应地设置该输出格式的标志。
     * 这样得到的标志会覆盖这个对象以前所有的内容。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };


  /**
   * 控制封装的postscript格式输出细节的标志。
   * @ingroup output
   *
   */
  struct EpsFlags : public OutputFlagsBase<EpsFlags>
  {
    /**
     * 这表示将用于生成高度信息的数据向量的编号。默认情况下，如果有任何数据向量，就取第一个数据向量，即<tt>height_vector==0</tt>。如果没有数据向量，就不生成高度信息。
     *
     */
    unsigned int height_vector;

    /**
     * 将被用于给单元格着色的向量的编号。与#height_vector同样适用。
     *
     */
    unsigned int color_vector;

    /**
     * 枚举表示是否应该进行缩放，使给定的<tt>size</tt>等于结果图片的宽度或高度。
     *
     */
    enum SizeType
    {
      /// Scale to given width
      width,
      /// Scale to given height
      height
    };

    /**
     * 见上文。默认为<tt>width</tt>。
     *
     */
    SizeType size_type;

    /**
     * 以postscript单位给出的输出的宽度或高度
     * 这通常是以奇怪的单位1/72英寸给出。这是高度还是宽度，由标志<tt>size_type</tt>指定。
     * 默认值是300，这代表了大约10厘米的尺寸。
     *
     */
    unsigned int size;

    /**
     * 一行的宽度，以postscript为单位。默认值是0.5。
     *
     */
    double line_width;

    /**
     * 线条原点-查看器对Z轴的角度，单位是度。
     * 默认为Gnuplot默认的60。
     *
     */
    double azimut_angle;

    /**
     * 观察者的位置投射到x-y平面上，围绕z轴旋转的角度，从上面看时为正数。
     * 单位是度，零等于高于或低于负y轴的位置。
     * 默认值是Gnuplot的默认值30。
     * 下面是一个Gnuplot默认值为0的例子。
     * @verbatim
     *
     *        3________7
     *        /       /|
     *       /       / |
     *     2/______6/  |
     *     |   |   |   |
     * O-->  |   0___|___4
     *     |  /    |  /
     *     | /     | /
     *    1|/______5/
     *
     * @endverbatim
     *
     *
     */
    double turn_angle;

    /**
     * 与x轴和y轴相比，z轴要被拉伸的系数。这是为了补偿坐标值和溶液值的不同大小，并防止绘图看起来很不合适（如果溶液值比坐标值小得多，则完全没有高程，反之则是常见的
     * "极度山区"。        默认为<tt>1.0</tt>。
     *
     */
    double z_scaling;

    /**
     * 标志，决定是否要绘制单元格（或每个补丁的部分）的边界线。
     * 默认值。<tt>true</tt>。
     *
     */
    bool draw_mesh;

    /**
     * 标记是否要填充单元格边界线之间的区域。如果不这样做，就不进行隐藏线的清除，在这个粗略的实现中，是通过从后向前的顺序写单元格来完成的，从而将背景中的单元格隐藏在前景中的单元格中。
     * 如果这个标志是<tt>false</tt>，并且#draw_mesh也是<tt>false</tt>，就不会打印任何东西。
     * 如果这个标志是<tt>true</tt>，那么单元格将被画成一个数据集的颜色（如果#shade_cells是<tt>true</tt>），或者纯白色（如果#shade_cells是假的或者没有数据集）。
     * 默认为<tt>true</tt>。
     *
     */
    bool draw_cells;

    /**
     * 决定单元格是否应被#color_vector表示的数据集着色，或简单地涂成白色的标志。这个标志只有在<tt>#draw_cells==true</tt>时才有意义。着色是通过#color_function完成的。
     * 默认为<tt>true</tt>。
     *
     */
    bool shade_cells;

    /**
     * 保持RGB系统中三个颜色值的结构。
     *
     */
    struct RgbValues
    {
      float red;
      float green;
      float blue;

      /**
       * 如果三个颜色值所代表的颜色是灰度，即所有成分都相等，则返回<tt>true</tt>。
       *
       */
      bool
      is_grey() const;
    };

    /**
     * 一个函数指针类型的定义，接收一个值并返回RGB值的三色值。
     * 除了要计算颜色的实际值外，还要给出要着色的数据的最小和最大值。
     *
     */
    using ColorFunction = RgbValues (*)(const double value,
                                        const double min_value,
                                        const double max_value);

    /**
     * 这是一个指向用于给单元格着色的函数的指针。
     * 默认情况下，它指向静态函数default_color_function()，它是这个类的成员。
     *
     */
    ColorFunction color_function;


    /**
     * 默认的着色函数。这个函数做的是人们通常想要的。它将颜色从黑色（最低值）通过蓝色、绿色和红色转移到白色（最高值）。关于色阶的确切定义，请参考实现。
     * 这个函数最初是由Stefan Nauber编写的。
     *
     */
    static RgbValues
    default_color_function(const double value,
                           const double min_value,
                           const double max_value);

    /**
     * 这是一个替代性的颜色函数，产生黑色（最低值）和白色（最高值）之间的灰阶。你可以通过将#color_function变量设置为这个函数的地址来使用它。
     *
     */
    static RgbValues
    grey_scale_color_function(const double value,
                              const double min_value,
                              const double max_value);

    /**
     * 这是另外一个可供选择的颜色函数，产生介于白色（最低值）和黑色（最高值）之间的灰度，也就是说，刻度与之前的刻度是相反的。你可以通过将#color_function变量设置为这个函数的地址来使用它。
     *
     */
    static RgbValues
    reverse_grey_scale_color_function(const double value,
                                      const double min_value,
                                      const double max_value);

    /**
     * 构造函数。
     *
     */
    EpsFlags(const unsigned int  height_vector  = 0,
             const unsigned int  color_vector   = 0,
             const SizeType      size_type      = width,
             const unsigned int  size           = 300,
             const double        line_width     = 0.5,
             const double        azimut_angle   = 60,
             const double        turn_angle     = 30,
             const double        z_scaling      = 1.0,
             const bool          draw_mesh      = true,
             const bool          draw_cells     = true,
             const bool          shade_cells    = true,
             const ColorFunction color_function = &default_color_function);

    /**
     * 用本类提供的名称和类型声明所有标志，以便在输入文件中使用。
     * 对于着色，只提供本类中声明的颜色函数。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取 declare_parameters()
     * 中声明的参数，并为这个输出格式设置相应的标志。
     * 这样得到的标志会覆盖这个对象以前所有的内容。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * 控制GMV格式的输出细节的标志。目前没有实施任何标志。
   * @ingroup output
   *
   */
  struct GmvFlags : public OutputFlagsBase<GmvFlags>
  {};

  /**
   * 控制Tecplot格式输出细节的标志。
   * @ingroup output
   *
   */
  struct TecplotFlags : public OutputFlagsBase<TecplotFlags>
  {
    /**
     * Tecplot允许为区段分配名称。这个变量存储了这个名称。
     *
     */
    const char *zone_name;

    /**
     * 绞股蓝中每个区域的解决时间。该值必须为非负值，否则将不会被写入文件。如果是静态区，不要给它指定任何值。
     *
     */
    double solution_time;

    /**
     * 构造函数。
     *
     */
    TecplotFlags(const char * zone_name     = nullptr,
                 const double solution_time = -1.0);

    /**
     * 返回这个对象的内存消耗估计值，单位是字节。
     *
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * 控制VTK格式的输出细节的标志。
   * @ingroup output
   *
   */
  struct VtkFlags : public OutputFlagsBase<VtkFlags>
  {
    /**
     * 如果此文件是时间相关模拟的一部分，则为时间步长的时间。
     * 该变量的值将根据http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
     * 中提供的指示写入输出文件，除非它的默认值为
     * @verbatim std::numeric_limits<unsigned int>::min() @endverbatim.
     *
     */
    double time;

    /**
     * 如果该文件是时间相关模拟的一部分，则是时间步数，或者是非线性或其他迭代中的周期。
     * 该变量的值将根据http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
     * 中提供的指令写入输出文件，除非它的默认值为
     * @verbatim std::numeric_limits<unsigned int>::min() @endverbatim.
     *
     */
    unsigned int cycle;

    /**
     * 决定当前日期和时间是否应作为注释打印在文件的第二行的标志。
     * 默认为<tt>true</tt>。
     *
     */
    bool print_date_and_time;

    /**
     * 一个提供不同的zlib压缩级别的数据类型。这些直接映射到由zlib定义的常数。
     *
     */
    enum ZlibCompressionLevel
    {
      /**
       * 不使用任何压缩。
       *
       */
      no_compression,
      /**
       * 使用最快的可用压缩算法。
       *
       */
      best_speed,
      /**
       * 使用能产生最小的压缩文件的算法。这是默认的标志。
       *
       */
      best_compression,
      /**
       * 使用默认的压缩算法。这是在速度和文件大小之间的一个折衷。
       *
       */
      default_compression
    };

    /**
     * 决定运行zlib（如果有的话）的压缩级别的标志。默认是<tt>best_compression</tt>。
     *
     */
    ZlibCompressionLevel compression_level;

    /**
     * 决定是否将补丁写成线性单元或高阶拉格朗日单元的标志。
     * 默认为<tt>false</tt>。
     * @note
     * 写入对应于高阶多项式而不是简单的线性或双线性的数据的能力，是2017年12月VTK
     * 8.1.0中才引入的一个功能。你至少需要在2018年4月发布的Paraview
     * 5.5.0版本或类似的最新版本的VisIt才能使用这个功能（例如，2020年2月发布的VisIt
     * 3.1.1还不支持这个功能）。这些程序的旧版本在试图读取用这个标志设置为
     * "真
     * "生成的文件时，很可能会导致错误。这些程序的经验表明，这些错误信息很可能不太具有描述性，而且比较隐晦。
     *
     */
    bool write_higher_order_cells;

    /**
     * 构造器。
     *
     */
    VtkFlags(
      const double       time  = std::numeric_limits<double>::min(),
      const unsigned int cycle = std::numeric_limits<unsigned int>::min(),
      const bool         print_date_and_time              = true,
      const ZlibCompressionLevel compression_level        = best_compression,
      const bool                 write_higher_order_cells = false);
  };


  /**
   * 用于SVG输出的标志。
   * @ingroup output
   *
   */
  struct SvgFlags : public OutputFlagsBase<SvgFlags>
  {
    /**
     * 图像的高度，以SVG为单位。默认值为4000。
     *
     */
    unsigned int height;

    /**
     * 图像的宽度，以SVG为单位。如果留为零，宽度将由高度计算。
     *
     */
    unsigned int width;

    /**
     * 这表示将用于生成高度信息的数据向量的编号。默认情况下，如果有任何数据向量的话，将使用第一个数据向量，即：<tt>#height_vector==0</tt>。如果没有数据向量，就不产生高度信息。
     *
     */
    unsigned int height_vector;

    /**
     * 用于透视图的角度
     *
     */
    int azimuth_angle, polar_angle;

    unsigned int line_thickness;

    /**
     * 在绘制的区域周围画出5%的边距
     *
     */
    bool margin;

    /**
     * 绘制一个编码单元格着色的色条
     *
     */
    bool draw_colorbar;

    /**
     * 构造函数。
     *
     */
    SvgFlags(const unsigned int height_vector  = 0,
             const int          azimuth_angle  = 37,
             const int          polar_angle    = 45,
             const unsigned int line_thickness = 1,
             const bool         margin         = true,
             const bool         draw_colorbar  = true);
  };


  /**
   * 控制以deal.II中间格式输出的细节的标志。
   * 目前没有实现任何标志。
   * @ingroup output
   *
   */
  struct Deal_II_IntermediateFlags
    : public OutputFlagsBase<Deal_II_IntermediateFlags>
  {
    /**
     * 用于编写中间格式的当前文件格式版本的指标。我们并不试图向后兼容，所以这个数字只是用来验证我们所写的格式是当前读者和写者所理解的。
     *
     */
    static const unsigned int format_version;
  };

  /**
   * 控制DataOutFilter的标志。
   * @ingroup output
   *
   */
  struct DataOutFilterFlags
  {
    /**
     * 过滤重复的顶点和相关值。这将极大地减少输出数据的大小，但是如果数据对应的是不连续的领域，则会导致输出文件不能忠实地代表实际数据。特别是，沿着子域的边界，数据仍然是不连续的，而在子域内部则看起来是一个连续的场。
     *
     */
    bool filter_duplicate_vertices;

    /**
     * XDMF输出是否指的是HDF5文件。这将影响到输出的结构方式。
     *
     */
    bool xdmf_hdf5_output;

    /**
     * 构造函数。
     *
     */
    DataOutFilterFlags(const bool filter_duplicate_vertices = false,
                       const bool xdmf_hdf5_output          = false);

    /**
     * 用这个类提供的名称和类型来声明所有的标志，以便在输入文件中使用。
     *
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * 读取<tt>declare_parameters</tt>中声明的参数，并相应地设置该输出格式的标志。
     * 这样得到的标志会覆盖这个对象以前所有的内容。
     *
     */
    void
    parse_parameters(const ParameterHandler &prm);

    /**
     * 确定此对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * DataOutFilter提供了一种去除deal.II输出所产生的多余顶点和数值的方法。默认情况下，DataOutBase和建立在它之上的类在每个单元格的每个角上输出数据。这意味着数据对网格的每个顶点都要输出多次。这个方案的目的是支持不连续量的输出，要么是因为有限元空间是不连续的，要么是因为输出的量是由一个解场计算出来的，并且在各面之间是不连续的。
   * 这个类别是为了控制写入的数据量。
   * 如果写入文件的字段确实是不连续的，那么忠实地表示它们的唯一方法就是为每个顶点写入多个值（这通常是通过为同一个顶点写入多个节点位置并在这些节点上定义数据来实现的）。然而，对于精细的网格，人们不一定对输出场的精确表示感兴趣，因为输出场可能只有小的不连续。相反，每个顶点只输出一个值就足够了，这个值可以从任何相邻单元的顶点上定义的值中任意选择。
   *
   */
  class DataOutFilter
  {
  public:
    /**
     * 默认构造函数。
     *
     */
    DataOutFilter();

    /**
     * 带有一组给定标志的解构器。参见DataOutFilterFlags以了解可能的标志。
     *
     */
    DataOutFilter(const DataOutBase::DataOutFilterFlags &flags);

    /**
     * 将一个具有指定索引的点写入过滤后的数据集中。如果该点已经存在，并且我们正在过滤多余的值，所提供的索引将在内部指代另一个被记录的点。
     *
     */
    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &p);

    /**
     * 以内部重新排序的格式记录一个deal.II单元。
     *
     */
    template <int dim>
    void
    write_cell(const unsigned int index,
               const unsigned int start,
               const unsigned int d1,
               const unsigned int d2,
               const unsigned int d3);

    /**
     * 以内部重新排序的格式记录一个没有分区的单一deal.II单元（例如：单数）。
     *
     */
    void
    write_cell_single(const unsigned int index,
                      const unsigned int start,
                      const unsigned int n_points);

    /**
     * 过滤并记录一个数据集。如果在一个给定的顶点有多个值，并且冗余的值被删除，则任意选择一个作为记录的值。在未来，这可以扩展到在一个给定的顶点上的平均/最小/最大多个值。
     *
     */
    void
    write_data_set(const std::string &     name,
                   const unsigned int      dimension,
                   const unsigned int      set_num,
                   const Table<2, double> &data_vectors);

    /**
     * 调整大小，用所有过滤后的节点顶点填充一个向量，输出到文件。
     *
     */
    void
    fill_node_data(std::vector<double> &node_data) const;

    /**
     * 调整一个向量的大小，并填充所有过滤后的单元格顶点指数，以输出到文件中。
     *
     */
    void
    fill_cell_data(const unsigned int         local_node_offset,
                   std::vector<unsigned int> &cell_data) const;

    /**
     * 获取数据集编号所指示的数据集名称。
     *
     */
    std::string
    get_data_set_name(const unsigned int set_num) const;

    /**
     * 获取由数据集编号指示的数据集的维度。
     *
     */
    unsigned int
    get_data_set_dim(const unsigned int set_num) const;

    /**
     * 获取数据集编号所示数据集的原始双值数据。
     *
     */
    const double *
    get_data_set(const unsigned int set_num) const;

    /**
     * 返回这个DataOutFilter中的节点数。如果启用了过滤功能，这可能小于原始节点数。
     *
     */
    unsigned int
    n_nodes() const;

    /**
     * 返回这个DataOutFilter中被过滤的单元格的数量。单元没有被过滤，所以这将是原始的单元数。
     *
     */
    unsigned int
    n_cells() const;

    /**
     * 返回这个DataOutFilter中被过滤的数据集的数量。数据集没有被过滤，所以这将是原始数据集的数量。
     *
     */
    unsigned int
    n_data_sets() const;

    /**
     * 做基类继承的空函数。
     *
     */
    void
    flush_points();

    /**
     * 空函数来做基类的继承。
     *
     */
    void
    flush_cells();


  private:
    /**
     * 为Map3DPoint提供比较函数的空类。
     *
     */
    struct Point3Comp
    {
      bool
      operator()(const Point<3> &one, const Point<3> &two) const
      {
        /*下面的返回语句是以下代码的优化版本： for (unsigned int d=0; d<3; ++d){ if (one(d) < two(d)) return true; else if (one(d) > two(d) return false; } return false;        
*
*/

        return (one(0) < two(0) ||
                (!(two(0) < one(0)) &&
                 (one(1) < two(1) || (!(two(1) < one(1)) && one(2) < two(2)))));
      }
    };

    using Map3DPoint = std::multimap<Point<3>, unsigned int, Point3Comp>;

    /**
     * 用于指定过滤行为的标志。
     *
     */
    DataOutBase::DataOutFilterFlags flags;

    /**
     * 当前对象所代表的顶点所处的空间维度的数量。这对应于通常的
     * @p dim
     * 参数，但由于这个类不是以维度为模板的，所以我们需要在这里存储它。
     *
     */
    unsigned int node_dim;

    /**
     * 存储在 @ref
     * filtered_cells 中的单元格的数量。
     *
     */
    unsigned int num_cells;

    /**
     * 点的映射到内部索引。
     *
     */
    Map3DPoint existing_points;

    /**
     * 实际点的索引与内部点的索引的映射。
     *
     */
    std::map<unsigned int, unsigned int> filtered_points;

    /**
     * 将单元格映射到过滤后的点。
     *
     */
    std::map<unsigned int, unsigned int> filtered_cells;

    /**
     * 数据集名称。
     *
     */
    std::vector<std::string> data_set_names;

    /**
     * 数据集的尺寸。
     *
     */
    std::vector<unsigned int> data_set_dims;

    /**
     * 数据集的数据。
     *
     */
    std::vector<std::vector<double>> data_sets;

    /**
     * 记录基于内部重排的单元格顶点索引。
     *
     */
    void
    internal_add_cell(const unsigned int cell_index,
                      const unsigned int pt_index);
  };


  /**
   * 提供一个数据类型，指定目前支持的输出格式。
   *
   */
  enum OutputFormat
  {
    /**
     * 使用已经存储在对象中的格式。
     *
     */
    default_format,

    /**
     * 不写任何输出。
     *
     */
    none,

    /**
     * 为OpenDX输出。
     *
     */
    dx,

    /**
     * 为AVS提供UCD格式的输出。
     *
     */
    ucd,

    /**
     * 为Gnuplot工具的输出。
     *
     */
    gnuplot,

    /**
     * 为Povray光线跟踪器提供的输出。
     *
     */
    povray,

    /**
     * 在封装的PostScript中输出。
     *
     */
    eps,

    /**
     * 用于GMV的输出。
     *
     */
    gmv,

    /**
     * 为Tecplot提供文本格式的输出。
     *
     */
    tecplot,

    /**
     * 以二进制格式输出Tecplot。比文本格式更快、更小。
     * @deprecated 使用Tecplot的二进制输出已被弃用。
     *
     */
    tecplot_binary,

    /**
     * 以VTK格式输出。
     *
     */
    vtk,

    /**
     * 以VTK格式输出。
     *
     */
    vtu,

    /**
     * 以SVG格式输出。
     *
     */
    svg,

    /**
     * 以deal.II中间格式输出。
     *
     */
    deal_II_intermediate,

    /**
     * 以HDF5格式输出。
     *
     */
    hdf5
  };


  /**
   * 将给定的补丁列表以OpenDX格式写入输出流。
   *
   */
  template <int dim, int spacedim>
  void
  write_dx(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &            nonscalar_data_ranges,
    const DXFlags &flags,
    std::ostream & out);

  /**
   * 将给定的补丁列表以eps格式写到输出流中。
   * 以这种格式输出可以规避使用辅助图形程序将一些输出格式转换为图形格式。这样做的好处是输出简单而快速，缺点是你必须给出一大堆参数，这些参数决定了视线的方向、着色的模式、高度轴的缩放等等（当然，所有这些参数都有合理的默认值，你可能想改变它们）。
   * 这个函数只支持二维域的输出（即dim=2），垂直方向的数值取自数据矢量。
   * 基本上，输出包含了网格和它们之间的单元。你可以画其中之一，或者两者都画，或者不画，如果你真的对一个空的图片感兴趣的话。如果写出来，网格使用黑线。网格之间的单元格要么不打印（这将导致隐线去除的损失，即你可以
   * "看穿
   * "单元格到后面的线条），要么打印成白色（除了隐线去除外没有任何作用），要么使用其中一个数据向量（不必与计算高度信息的向量相同）和一个可定制的颜色函数进行着色。默认的颜色函数在黑色、蓝色、绿色、红色和白色之间选择颜色，选择的数据字段的值越来越大，用于着色。目前，每个单元格只显示一种颜色，该颜色取自单元格中心的数据字段的值；不使用单元格上颜色的双线性插值。
   * 默认情况下，视角的选择与GNUPLOT中的默认视角一样，即相对于正Z轴的角度为60度，相对于负Y轴的角度正向旋转30度（从上面看）。
   * 当然，你可以改变这些设置。
   * 编写EPS输出时，图片周围没有边界，也就是说，边界框在四面都靠近输出。坐标最多使用五位数来书写，以保持图片的合理尺寸。
   * 所有的参数以及它们的默认值都列在这个类的<tt>EpsFlags</tt>成员类的文档中。更多详细的信息请见那里。
   *
   */
  template <int spacedim>
  void
  write_eps(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string> &       data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const EpsFlags &flags,
    std::ostream &  out);

  /**
   * 这是一个与上面相同的函数，除了用于非二维的域。这个函数没有被实现（如果被调用将抛出一个错误），但被声明是为了允许与维度无关的程序。
   *
   */
  template <int dim, int spacedim>
  void
  write_eps(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const EpsFlags &flags,
    std::ostream &  out);


  /**
   * 将给定的补丁列表以GMV格式写到输出流中。
   * 数据以如下格式写入：节点被认为是斑块的点。在空间维度小于3的情况下，缺失的坐标会插入0。数据向量被写成节点或单元数据，对于第一种，数据空间被内插为（双，三）线性元素。
   *
   */
  template <int dim, int spacedim>
  void
  write_gmv(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const GmvFlags &flags,
    std::ostream &  out);

  /**
   * 将给定的补丁列表以gnuplot格式写到输出流中。
   * 然后可以通过启动<tt>gnuplot</tt>并输入命令来实现二维数据的可视化
   * @verbatim
   * set data style lines
   * splot "filename" using 1:2:n
   * @endverbatim
   * 本例假设显示的数据矢量的编号为<b>n-2</b>。
   * GNUPLOT格式不能直接处理非结构化网格上的数据。直接意味着你只给出顶点和其上的解值，而程序会构建自己的网格来表示这些数据。这只对二维的结构化张量积网格是可行的。然而，在一个文件内给出几个这样的补丁是可能的，这正是这个类的相应函数所做的：将每个单元格的数据写成一个数据补丁，至少在派生类传递的补丁代表单元格的情况下。请注意，补丁上的函数在补丁之间的界面上不需要是连续的，所以这个方法也适用于不连续的元素。还要注意的是，GNUPLOT可以对修补过的数据进行隐线去除。
   * 虽然这个讨论适用于两个空间维度，但在三维空间中则更为复杂。原因是我们仍然可以使用补丁，但当我们试图将它们可视化时就很困难了，因为如果我们使用切开数据（例如，通过使用x和z坐标，一个固定的y值和z方向的绘图函数值，那么补丁数据就不是GNUPLOT想要的补丁了。因此，我们使用了另一种方法，即把数据写在三维网格上，作为一连串的线，即每两个点与一个或多个数据集有关。
   * 因此，一个补丁的每个子单元有12条线。
   * 鉴于上述的线条，在Gnuplot中可以实现对这些数据的切割，就像这样。
   * @verbatim
   * set data style lines
   * splot [:][:][0:] "T" using 1:2:(\$3==.5 ? \$4 :
   *
   * -1)
   * @endverbatim
   * 这个命令在 $x$ -和 $y$ -方向上无限制地绘制数据，但在
   * $z$ -方向上只绘制那些在 $x$ - $y$
   * -平面以上的数据点（我们在这里假设一个正解，如果它有负值，你可能想减少下限）。此外，它只取z值（<tt>&3</tt>）等于0.5的数据点，即在<tt>z=0.5</tt>处切过域。对于这个平面上的数据点，第一个数据集（<tt>&4</tt>）的数据值在x-y平面上方的z方向上被提高；所有其他的点都被表示为<tt>-1</tt>的值，而不是数据向量的值，并且由于z绘图方向的下限，在第三对括号中给出，所以不被绘制。
   * 更复杂的切割是可能的，包括非线性的切割。然而，请注意，只有那些实际在切割面上的点才被绘制出来。
   *
   */
  template <int dim, int spacedim>
  void
  write_gnuplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &                 nonscalar_data_ranges,
    const GnuplotFlags &flags,
    std::ostream &      out);

  /**
   * 将给定的补丁列表写入Povray光线跟踪器的输出流中。    以这种格式输出会创建一个povray源文件，包括用povray 3.1渲染的标准相机和光源定义，目前，这种格式只支持二维数据的输出，第三个方向的数值取自数据矢量。    输出使用两个不同的povray-objects。      <ul>   <li>  <tt>BICUBIC_PATCH</tt> 一个<tt>bicubic_patch</tt>是一个3维的Bezier补丁。它由16个描述表面的点组成。4个角点被物体所接触，而其他12个点则拉动和拉伸补丁的形状。每个补丁上都会产生一个<tt>bicubic_patch</tt>。因此，细分的数量必须是3，以提供16个点的补丁。双三次元补丁并不精确，但能生成非常平滑的图像。      <li>  <tt>MESH</tt> 网格对象是用来存储大量的三角形的。补丁数据的每个正方形都被分割成一个左上角和一个右下角的三角形。如果细分的数量是3个，每个补丁就会产生32个三角形。    使用平滑标志povray在三角形上插值法线，模仿一个弯曲的表面 </ul> 。这个纹理必须在对象数据之前的某个地方声明。这可能是在一个外部数据文件中或在输出文件的开头。将<tt>external_data</tt>标志设置为false，一个标准的摄像机、灯光和纹理（按比例调整以适应场景）会被添加到输出文件中。设置为 "true"，一个包含文件 "data.inc "会被包括在内。这个文件不是由deal生成的，必须包括摄像机、灯光和纹理定义Tex。    你需要povray（>=3.0）来渲染这个场景。povray的最小选项是。
   * @verbatim
   * povray +I<inputfile> +W<horiz. size> +H<ver. size> +L<include path>
   * @endverbatim
   * 如果使用外部文件
   * "data.inc"，这个文件的路径必须包含在povray选项中。
   *
   */
  template <int dim, int spacedim>
  void
  write_povray(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &                nonscalar_data_ranges,
    const PovrayFlags &flags,
    std::ostream &     out);

  /**
   * 将给定的补丁列表以Tecplot
   * ASCII格式（FEBLOCK）写入输出流。
   * 更多信息请查阅Tecplot用户和参考手册。
   *
   */
  template <int dim, int spacedim>
  void
  write_tecplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &                 nonscalar_data_ranges,
    const TecplotFlags &flags,
    std::ostream &      out);

  /**
   * 将给定的补丁列表以AVS开发者指南（现在的AVS）中描述的UCD格式写到输出流中。由于目前格式的限制，只能输出基于节点的数据，这也是我们发明补丁概念的原因之一。为了编写高阶元素，你可以把它们分割成每个单元的几个子单元。
   * 然而，这些子单元也会被理解UCD格式的程序显示为不同的单元。
   * 我们没有利用提供模型数据的可能性，因为所有UCD程序都不支持这些数据。你可以在派生类中给出单元格数据，方法是将一个补丁上给定数据集的所有值设置为相同的值。
   *
   */
  template <int dim, int spacedim>
  void
  write_ucd(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const UcdFlags &flags,
    std::ostream &  out);

  /**
   * 将给定的补丁列表以VTK格式写到输出流中。数据是以传统的VTK格式写入的，而不是write_vtu()产生的基于XML的格式。
   * nonscalar_data_ranges参数表示输出中的组件范围，被认为是一个矢量，而不是简单的标量字段的集合。VTK的输出格式有特殊的规定，允许这些组件用一个名字来输出，而不是在可视化程序中把几个标量字段归为一个矢量。
   * @note
   * VTK是一种遗留格式，在很大程度上已经被VTU格式（VTK的XML结构版本）所取代了。特别是，VTU允许对数据进行压缩，因此导致大文件的文件大小要比VTK文件小得多。由于所有支持VTK的可视化程序也支持VTU，你应该考虑使用后者的文件格式，通过使用write_vtu()函数来代替。
   *
   */
  template <int dim, int spacedim>
  void
  write_vtk(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream &  out);


  /**
   * 将给定的补丁列表以VTU格式写入输出流中。数据是以基于XML的VTK格式写入的，而不是write_vtk()产生的传统格式。
   * nonscalar_data_ranges参数表示输出中组件的范围，被认为是一个矢量，而不是简单的标量字段的集合。VTK的输出格式有特殊的规定，允许这些组件用一个名字来输出，而不是在可视化程序中把几个标量字段归为一个矢量。
   * 一些可视化程序，如ParaView，可以读取几个独立的VTU文件来实现可视化的并行化。在这种情况下，你需要一个
   * <code>.pvtu</code> 文件来描述哪些VTU文件构成一个组。
   * DataOutInterface::write_pvtu_record()
   * 函数可以生成这样一个集中的记录。同样，
   * DataOutInterface::write_visit_record()
   * 对VisIt也有同样的作用（尽管VisIt从2.5.1版开始也可以读取
   * <code>pvtu</code>
   * 记录）。最后，对于与时间有关的问题，你可能还想看看
   * DataOutInterface::write_pvd_record()  这个函数的使用在  step-40
   * 中有解释。
   *
   */
  template <int dim, int spacedim>
  void
  write_vtu(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream &  out);

  /**
   * 这是为基于xml的vtu文件格式写头。这个例程与
   * DataOutInterface::write_vtu_footer() 和
   * DataOutInterface::write_vtu_main() 一起被 DataOutBase::write_vtu().
   * 内部使用。
   *
   */
  void
  write_vtu_header(std::ostream &out, const VtkFlags &flags);

  /**
   * 该函数为基于xml的vtu文件格式写入页脚。本例程与
   * DataOutInterface::write_vtu_header() 和
   * DataOutInterface::write_vtu_main() 一起被 DataOutBase::write_vtu().
   * 内部使用。
   *
   */
  void
  write_vtu_footer(std::ostream &out);

  /**
   * 该函数为基于xml的vtu文件格式写入主要部分。这个程序在内部与
   * DataOutInterface::write_vtu_header() 和
   * DataOutInterface::write_vtu_footer() 一起被 DataOutBase::write_vtu().
   * 使用。
   *
   */
  template <int dim, int spacedim>
  void
  write_vtu_main(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream &  out);

  /**
   * 一些可视化程序，如ParaView，可以读取几个独立的VTU文件，这些文件都是同一模拟的一部分，以实现可视化的并行。在这种情况下，你需要一个
   * <code>.pvtu</code> 文件来描述哪些VTU文件（例如，通过
   * DataOutInterface::write_vtu() 函数写入）构成一个组。
   * 当前的函数可以生成这样一个集中的记录。
   * 这个函数通常不会从用户空间自己调用，但你可能想通过
   * DataOutInterface::write_pvtu_record()
   * 调用它，因为DataOutInterface类可以访问你必须手工提供给当前函数的信息。
   * 在任何情况下，不管是直接调用这个函数还是通过
   * DataOutInterface::write_pvtu_record(),
   * 调用，这样编写的中央记录文件都包含一个（标量或矢量）字段的列表，描述哪些字段实际上可以在构成平行VTU文件集的各个文件中找到，以及这些文件的名称。这个函数通过第三和第四个参数获得字段的名称和类型；你可以用手来确定这些，但在实践中，这个函数最容易通过调用
   * DataOutInterfaces::write_pvtu_record(), 来调用，它通过调用
   * DataOutInterface::get_dataset_names() 和
   * DataOutInterface::get_nonscalar_data_ranges()
   * 函数确定最后两个参数。这个函数的第二个参数指定了构成并行集的文件名称。
   * @note  使用 DataOutBase::write_vtu() 和 DataOutInterface::write_vtu()
   * 来写每一块。还要注意，只有一个并行进程需要调用当前函数，列出所有并行进程写入的文件名。
   * @note  为了告诉Paraview将多个 <code>pvtu</code>
   * 文件组合在一起，每个文件描述一个与时间有关的模拟的一个时间步骤，请参阅
   * DataOutBase::write_pvd_record() 函数。
   * @note  旧版本的VisIt（2.5.1之前），不能读取
   * <code>pvtu</code>
   * 记录。然而，它可以读取由write_visit_record()函数写入的访问记录。
   *
   */
  void
  write_pvtu_record(
    std::ostream &                  out,
    const std::vector<std::string> &piece_names,
    const std::vector<std::string> &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &nonscalar_data_ranges);

  /**
   * 在ParaView中，可以将与时间有关的数据可视化，并将其标记为与时间有关的模拟的当前积分时间。为了使用这个功能，你需要一个
   * <code>.pvd</code>
   * 文件，描述哪个VTU或PVTU文件属于哪个时间段。这个函数写入一个提供这种映射的文件，也就是说，它需要一个对的列表，每个对表示一个特定的时间瞬时和包含这个时间瞬时的图形数据的相应文件。
   * 一个典型的用例，在计算与时间有关的解决方案的程序中，将是以下内容（
   * <code>time</code> and <code>time_step</code> 是类型为
   * <code>double</code> 和 <code>unsigned int</code>
   * 的类的成员变量；变量 <code>times_and_names</code> 的类型是
   * <code>std::vector@<std::pair@<double,std::string@> @></code> ）。
   * @code
   * template <int dim>
   * void MyEquation<dim>::output_results () const
   * {
   * DataOut<dim> data_out;
   *
   * data_out.attach_dof_handler(dof_handler);
   * data_out.add_data_vector(solution, "U");
   * data_out.build_patches();
   *
   * const std::string filename = "solution-" +
   *                              Utilities::int_to_string (timestep_n, 3) +
   *                              ".vtu";
   * std::ofstream output(filename);
   * data_out.write_vtu(output);
   *
   * times_and_names.emplace_back (time, filename);
   * std::ofstream pvd_output ("solution.pvd");
   * DataOutBase::write_pvd_record (pvd_output, times_and_names);
   * }
   * @endcode
   * @note  参见 DataOutInterface::write_vtu,
   * DataOutInterface::write_pvtu_record, 和
   * DataOutInterface::write_vtu_in_parallel
   * ，用于编写每个时间步长的解决方案。
   * @note
   * 每对文件的第二个元素，即储存每个时间的图形数据的文件，本身又可以是一个引用其他文件的文件。例如，它可以是一个
   * <code>.pvtu</code>
   * 文件的名称，该文件引用了一个并行计算的多个部分。
   *
   */
  void
  write_pvd_record(
    std::ostream &                                     out,
    const std::vector<std::pair<double, std::string>> &times_and_names);

  /**
   * 这个函数完全等同于write_pvtu_record()函数，但适用于旧版本的VisIt可视化程序和一个可视化图形（或仅一个时间步长）。关于这个函数的用途，见那里。
   * 这个函数在 "将数据输入VisIt "报告中的
   * "创建并行的主文件
   * "部分（第5.7节）有记录，可以在这里找到：
   * https://wci.llnl.gov/codes/visit/2.0.0/GettingDataIntoVisIt2.0.0.pdf
   *
   */
  void
  write_visit_record(std::ostream &                  out,
                     const std::vector<std::string> &piece_names);

  /**
   * 这个函数等同于上面的write_visit_record()，但用于多个时间步长。下面是一个如何使用该函数的例子。
   * @code
   * const unsigned int number_of_time_steps = 3;
   * std::vector<std::vector<std::string > > piece_names(number_of_time_steps);
   *
   * piece_names[0].emplace_back("subdomain_01.time_step_0.vtk");
   * piece_names[0].emplace_back("subdomain_02.time_step_0.vtk");
   *
   * piece_names[1].emplace_back("subdomain_01.time_step_1.vtk");
   * piece_names[1].emplace_back("subdomain_02.time_step_1.vtk");
   *
   * piece_names[2].emplace_back("subdomain_01.time_step_2.vtk");
   * piece_names[2].emplace_back("subdomain_02.time_step_2.vtk");
   *
   * std::ofstream visit_output ("solution.visit");
   *
   * DataOutBase::write_visit_record(visit_output, piece_names);
   * @endcode
   * 这个函数在 "将数据输入VisIt "报告的
   * "创建一个并行的主文件
   * "一节（第5.7节）中有记录，可以在这里找到：
   * https://wci.llnl.gov/codes/visit/2.0.0/GettingDataIntoVisIt2.0.0.pdf
   *
   */
  void
  write_visit_record(std::ostream &                               out,
                     const std::vector<std::vector<std::string>> &piece_names);

  /**
   * 这个函数等同于上面的write_visit_record()，但是对于多个时间步长，而且每个时间步长的时间的附加信息。下面是一个如何使用该函数的例子。
   * @code
   * const unsigned int number_of_time_steps = 3;
   * std::vector<std::pair<double,std::vector<std::string > > >
   * times_and_piece_names(number_of_time_steps);
   *
   * times_and_piece_names[0].first = 0.0;
   * times_and_piece_names[0].second.emplace_back("subdomain_01.time_step_0.vtk");
   * times_and_piece_names[0].second.emplace_back("subdomain_02.time_step_0.vtk");
   *
   * times_and_piece_names[1].first = 0.5;
   * times_and_piece_names[1].second.emplace_back("subdomain_01.time_step_1.vtk");
   * times_and_piece_names[1].second.emplace_back("subdomain_02.time_step_1.vtk");
   *
   * times_and_piece_names[2].first = 1.0;
   * times_and_piece_names[2].second.emplace_back("subdomain_01.time_step_2.vtk");
   * times_and_piece_names[2].second.emplace_back("subdomain_02.time_step_2.vtk");
   *
   * std::ofstream visit_output ("solution.visit");
   *
   * DataOutBase::write_visit_record(visit_output, times_and_piece_names);
   * @endcode
   * 这个函数在 "将数据输入VisIt "报告的
   * "为并行创建主文件
   * "一节（第5.7节）中有记录，可以在这里找到：
   * https://wci.llnl.gov/codes/visit/2.0.0/GettingDataIntoVisIt2.0.0.pdf
   *
   */
  void
  write_visit_record(
    std::ostream &out,
    const std::vector<std::pair<double, std::vector<std::string>>>
      &times_and_piece_names);

  /**
   * 将给定的补丁列表以SVG格式写入输出流中。
   * SVG（Scalable Vector
   * Graphics）是一种基于XML的矢量图像格式，由万维网联盟（W3C）开发和维护。该功能符合2011年8月16日发布的最新规范SVG
   * 1.1。通过设置或清除相应的标志（见SvgFlags结构），可以控制图形的输出。目前，这种格式只支持二维数据的输出，第三个方向的值取自数据矢量。
   * 对于输出，每个补丁被细分为四个三角形，然后被写成多边形，并用线性颜色梯度填充。补丁产生的颜色使顶点的数据值从指定的数据向量中可视化。可以画一个色条来编码着色。
   * @note
   * 这个函数到目前为止只在两个维度上实现，并为数据信息保留了一个附加维度。
   *
   */
  template <int spacedim>
  void
  write_svg(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string> &       data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &             nonscalar_data_ranges,
    const SvgFlags &flags,
    std::ostream &  out);

  /**
   * 将给定的补丁列表以deal.II中间格式写到输出流中。这不是任何其他图形程序所能理解的格式，而是直接转储deal.II所使用的内部中间格式。这种内部格式是由可以使用DataOutBase类产生输出的各种类产生的，例如从有限元求解中产生，然后在本类中转换为最终的图形格式。
   * 注意，中间格式就像它的名字一样：内部数据的直接表示。它不是标准化的，每当我们改变内部表示时就会改变。你只能期望使用用于写入的同一版本的deal.II来处理以这种格式写入的文件。
   * 我们提供写出这种中间格式的原因是，它可以使用DataOutReader类读回deal.II程序中，这至少在两种情况下是有帮助的。首先，这可以用来在以后生成任何其他目前能理解的图形格式的图形输出；这样，在运行时就不需要知道要求哪种输出格式，或者是否需要不同格式的多个输出文件。其次，与几乎所有其他图形格式相比，有可能合并几个包含中间格式数据的文件，并从中生成一个单一的输出文件，该文件可以再次采用中间格式或任何最终格式。后一种选择对并行程序最有帮助：正如
   * step-17
   * 示例程序所演示的，可以只让一个处理器为整个并行程序生成图形输出，但如果涉及许多处理器，这可能会变得效率极低，因为负载不再平衡。出路是让每个处理器为它的那块领域生成中间图形输出，而后将不同的文件合并成一个，这是一个比生成中间数据便宜得多的操作。
   * 中间格式的deal.II数据通常存储在以<tt>.d2</tt>结尾的文件中。
   *
   */
  template <int dim, int spacedim>
  void
  write_deal_II_intermediate(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &                              nonscalar_data_ranges,
    const Deal_II_IntermediateFlags &flags,
    std::ostream &                   out);

  /**
   * 将 @p data_filter
   * 中的数据写入一个包含网格和求解值的HDF5文件。
   *
   */
  template <int dim, int spacedim>
  void
  write_hdf5_parallel(const std::vector<Patch<dim, spacedim>> &patches,
                      const DataOutFilter &                    data_filter,
                      const std::string &                      filename,
                      const MPI_Comm &                         comm);

  /**
   * 将 @p data_filter 中的数据写入HDF5文件。如果 @p
   * write_mesh_file
   * 为假，网格数据将不被写入，解文件将只包含解的数值。如果
   * @p write_mesh_file
   * 为真，且文件名相同，生成的文件将同时包含网格数据和求解值。
   *
   */
  template <int dim, int spacedim>
  void
  write_hdf5_parallel(const std::vector<Patch<dim, spacedim>> &patches,
                      const DataOutFilter &                    data_filter,
                      const bool                               write_mesh_file,
                      const std::string &                      mesh_filename,
                      const std::string &solution_filename,
                      const MPI_Comm &   comm);

  /**
   * DataOutFilter是一种中间数据格式，可以减少将被写入文件的数据量。这个函数所填充的对象随后可以再次用于写入具体文件格式的数据；例如，见
   * DataOutBase::write_hdf5_parallel(). 。
   *
   */
  template <int dim, int spacedim>
  void
  write_filtered_data(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string> &         data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
      &            nonscalar_data_ranges,
    DataOutFilter &filtered_data);

  /**
   * 给出一个包含由write_deal_II_intermediate()写入的数据的输入流，确定调用该函数的<tt>dim</tt>和<tt>spacedim</tt>模板参数，并将它们作为一对值返回。
   * 注意，这个函数在流的当前位置吃了一些元素，因此改变了它。为了使用例如DataOutReader类从它那里读取数据，你可能希望将流重置到它以前的位置，或者关闭并重新打开它。
   *
   */
  std::pair<unsigned int, unsigned int>
  determine_intermediate_format_dimensions(std::istream &input);

  /**
   * 返回对应于给定字符串的OutputFormat值。如果该字符串与任何已知的格式不匹配，就会抛出一个异常。
   * 这个函数的主要目的是允许程序使用任何已实现的输出格式，而不需要在每次实现新的格式时扩展程序的分析器。
   * 要获得目前可用的格式名称的列表，例如，要把它交给ParameterHandler类，请使用函数get_output_format_names()。
   *
   */
  OutputFormat
  parse_output_format(const std::string &format_name);

  /**
   * 返回一个已实现的输出格式的列表。不同的名称由垂直条形符号（<tt>`|'</tt>）分开，正如ParameterHandler类所使用的那样。
   *
   */
  std::string
  get_output_format_names();

  /**
   * 提供一个函数，告诉我们一个给定输出格式的文件通常有哪个后缀。目前定义了以下格式。    <ul>   <li>  <tt>dx</tt>: <tt>.dx</tt>  <li>  <tt>ucd</tt>: <tt>.inp</tt>  <li>  <tt>gnuplot</tt>: <tt>.gnuplot</tt>  <li>  <tt>povray</tt>: <tt>.pov</tt>  <li>  <tt>eps</tt>: <tt>.eps</tt>  <li>  <tt>gmv</tt>: <tt>.gmv</tt>  <li>  <tt>tecplot</tt>: <tt>.dat</tt>  <li>  <tt>tecplot_binary</tt>: <tt>.plt</tt>  <li>  <tt>vtk</tt>: <tt>.vtk</tt>  <li>  <tt>vtu</tt>: <tt>.vtu</tt>  <li>  <tt>svg</tt>: <tt>.svg</tt>  <li>  <tt>deal_II_intermediate</tt>: <tt>.d2</tt>.     </ul>   @deprecated  使用Tecplot二进制输出已被废弃。
   *
   */
  std::string
  default_suffix(const OutputFormat output_format);

  /**
   * @addtogroup  Exceptions 
     * @{ 
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException2(ExcInvalidDatasetSize,
                 int,
                 int,
                 << "The number of points in this data set is " << arg1
                 << ", but we expected " << arg2
                 << " in each space direction.");
  /**
   * 一个输出函数没有收到任何用于写入的补丁。
   *
   */
  DeclExceptionMsg(ExcNoPatches,
                   "You are trying to write graphical data into a file, but "
                   "no data is available in the intermediate format that "
                   "the DataOutBase functions require. Did you forget to "
                   "call a function such as DataOut::build_patches()?");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcTecplotAPIError,
                   "The error code of one of the Tecplot functions was "
                   "not zero as expected.");
  /**
   * 异常情况
   *
   */
  DeclException1(ExcErrorOpeningTecplotFile,
                 char *,
                 << "There was an error opening Tecplot file " << arg1
                 << " for output.");

  //@}
} // namespace DataOutBase



/**
 * 这个类是DataOutBase命名空间中的函数的接口，正如其名字所暗示的那样。它没有提供太多的功能，只是提供了一种访问已实现的格式的方法和一种动态分配选择何种输出格式的方法。
 * 这个类被认为是实际生成输出数据的类的基类。它有两个抽象的虚拟函数，get_patches()和get_dataset_names()产生实际需要的数据。这些是唯一需要被派生类重载的函数。除此之外，它还有一个针对底层基类所支持的每种输出格式的函数，它使用这两个虚拟函数获得输出数据，并将它们传递给原始输出函数。
 * 这个类的目的主要有两个方面：支持存储标志，通过这些标志来控制不同输出格式的输出，并意味着以一种在运行时确定输出格式、标志和其他东西的方式来处理输出。除此之外，它还为上面简要讨论的派生类提供了抽象接口。
 *
 *  <h3>Output flags</h3>
 * 我们在这个类中处理标志的方式与<tt>GridOut</tt>类中使用的方式非常相似。关于为什么和如何处理的详细信息，以及一个编程的例子，我们参考该类的文档。
 * 基本上，这个类为底层的<tt>DataOutBase</tt>类所支持的每种输出格式存储了一组标志。只要使用<tt>write_*</tt>函数之一，就会用到这些标志。默认情况下，这些标志的值被设置为合理的启动，但如果你想改变它们，你可以创建一个结构，持有其中一种输出格式的标志，并使用这个类的<tt>set_flags</tt>函数来设置它，以确定该对象未来可能通过该输出格式产生的所有输出。
 * 关于不同输出函数支持哪些参数的信息，请参见<tt>DataOutBase</tt>类及其成员类的文档。
 *
 *  <h3>Run time selection of output parameters</h3>
 * 在上面描述的输出标志类中，为不同格式的输出定义了许多标志。为了使它们对输入文件处理程序类<tt>ParameterHandler</tt>可用，每个都有一个函数向参数处理程序声明这些标志，并从实际的输入文件中读回它们。为了避免在用户程序中必须为每个可用的输出格式和各自的标志类调用这些函数，目前的<tt>DataOutInterface</tt>类提供了一个函数<tt>declare_parameters</tt>，调用所有已知输出格式标志类的相应函数。每个此类格式的标志都被打包在输入文件的一个小节中。同样，还有一个函数<tt>parse_parameters</tt>，它读取这些参数并将它们存储在与此对象相关的标志中（见上文）。
 * 使用这些函数，你不必跟踪哪些格式是目前实现的。
 * 使用方法如下。
 *
 * @code
 * // within function declaring parameters:
 * prm.enter_subsection("Output format options");
 * DataOutInterface<dim>::declare_parameters(prm);
 * prm.leave_subsection();
 *
 * ...
 * // within function doing the output:
 * DataOut<dim> out;
 * prm.enter_subsection("Output format options");
 * out.parse_parameters(prm);
 * prm.leave_subsection();
 * @endcode
 * 注意，在本例中，使用了<tt>DataOut</tt>类。然而，任何从<tt>DataOutInterface</tt>派生的其他类都可以同样工作。
 *
 *  <h3>Run time selection of formats</h3>
 * 这个类，很像<tt>GridOut</tt>类，有一组函数提供支持的输出格式列表，一个<tt>enum</tt>表示所有这些，还有一个函数解析一个字符串，如果它是一个有效的输出格式的名称，则返回相应的<tt>enum</tt>值（实际上，这些函数是从基类继承的）。最后，有一个函数<tt>write</tt>，它接收这个<tt>enum</tt>的值，并根据这个值所选择的输出格式，分派给实际的<tt>write_*</tt>函数之一。
 * 提供不同输出格式名称的函数分别是：<tt>default_suffix</tt>,
 * <tt>parse_output_format</tt>, 和
 * <tt>get_output_format_names</tt>。它们使参数文件中输出格式的选择变得更加容易，尤其是独立于目前实现的格式。因此，每当一个新的格式被实施时，用户程序不需要改变。
 * 此外，这个类的对象有一个默认格式，可以通过参数文件的
 * "输出格式
 * "参数来设置。在一个程序中，这可以通过成员函数<tt>set_default_format</tt>来改变。使用这个默认格式，可以将格式选择完全留给参数文件。输出文件名的合适后缀可以通过不带参数的<tt>default_suffix</tt>获得。
 *
 *
 * @ingroup output
 *
 *
 */
template <int dim, int spacedim = dim>
class DataOutInterface
{
public:
  /**
   * 构建器。
   *
   */
  DataOutInterface();

  /**
   * 解构器。什么都不做，但由于这个类有虚拟函数，所以被声明为虚拟。
   *
   */
  virtual ~DataOutInterface() = default;

  /**
   * 通过get_patches()获取数据，并以OpenDX格式写入<tt>out</tt>。见
   * DataOutBase::write_dx. 。
   *
   */
  void
  write_dx(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以EPS格式写到<tt>out</tt>。参见
   * DataOutBase::write_eps.
   *
   */
  void
  write_eps(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据并以GMV格式写入<tt>out</tt>。参见
   * DataOutBase::write_gmv.
   *
   */
  void
  write_gmv(std::ostream &out) const;

  /**
   * 通过get_patches()获得数据并以GNUPLOT格式写入<tt>out</tt>。参见
   * DataOutBase::write_gnuplot.  。
   *
   */
  void
  write_gnuplot(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以POVRAY格式写到<tt>out</tt>。参见
   * DataOutBase::write_povray.
   *
   */
  void
  write_povray(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以Tecplot格式写到<tt>out</tt>。见
   * DataOutBase::write_tecplot.
   *
   */
  void
  write_tecplot(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以UCD格式写入<tt>out</tt>，用于AVS。参见
   * DataOutBase::write_ucd.
   *
   */
  void
  write_ucd(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以Vtk格式写到<tt>out</tt>。参见
   * DataOutBase::write_vtk.
   * @note
   * VTK是一种遗留格式，在很大程度上已经被VTU格式（VTK的XML结构版本）所取代了。特别是，VTU允许对数据进行压缩，因此导致大文件的文件大小要比VTK文件小得多。由于所有支持VTK的可视化程序也支持VTU，你应该考虑使用后者的文件格式，通过使用write_vtu()函数来代替。
   *
   */
  void
  write_vtk(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以Vtu（VTK的XML）格式写到<tt>out</tt>。参见
   * DataOutBase::write_vtu.
   * 一些可视化程序，如ParaView，可以读取几个独立的VTU文件以实现可视化的并行化。在这种情况下，你需要一个
   * <code>.pvtu</code> 文件来描述哪些VTU文件构成一个组。
   * DataOutInterface::write_pvtu_record()
   * 函数可以生成这样一个集中的记录。同样，
   * DataOutInterface::write_visit_record()
   * 对旧版本的VisIt也有同样的作用（尽管VisIt从2.5.1版本开始也可以读取
   * <code>pvtu</code> 记录）。最后，
   * DataOutInterface::write_pvd_record()
   * 可以用来将共同构成时间相关模拟的文件分组。
   *
   */
  void
  write_vtu(std::ostream &out) const;

  /**
   * 集体MPI调用，将所有参与节点（给定通信器中的节点）的解决方案写入共享文件系统上的一个压缩的.vtu文件。
   * 该通信器可以是计算所使用的通信器的一个子通信器。
   * 这个程序使用MPI
   * I/O来实现并行文件系统上的高性能。也可参见
   * DataOutInterface::write_vtu(). 。
   *
   */
  void
  write_vtu_in_parallel(const std::string &filename,
                        const MPI_Comm &   comm) const;

  /**
   * 一些可视化程序，如ParaView，可以读取几个独立的VTU文件，这些文件都是同一模拟的一部分，以实现可视化的并行化。在这种情况下，你需要一个
   * <code>.pvtu</code> 文件来描述哪些VTU文件（例如，通过
   * DataOutInterface::write_vtu() 函数写入）构成一个组。
   * 当前的函数可以生成这样一个集中记录。
   * 该函数生成的中央记录文件包含一个（标量或矢量）字段列表，描述哪些字段实际上可以在构成平行VTU文件组的各个文件中找到，以及这些文件的名称。这个函数通过本类的get_dataset_names()和get_nonscalar_data_ranges()函数获得字段的名称和类型。这个函数的第二个参数指定了构成平行集的文件名。
   * @note  使用 DataOutBase::write_vtu() 和 DataOutInterface::write_vtu()
   * 来写每一块。还要注意，只有一个并行进程需要调用当前函数，列出所有并行进程写入的文件名。
   * @note  这个函数的使用在  step-40  中解释。
   * @note  为了告诉Paraview将多个 <code>pvtu</code>
   * 文件组合在一起，每个文件描述一个与时间有关的仿真的一个时间步骤，请参见
   * DataOutBase::write_pvd_record() 函数。
   * @note  旧版本的VisIt（2.5.1之前），不能读取
   * <code>pvtu</code>
   * 记录。然而，它可以读取由write_visit_record()函数写入的访问记录。
   *
   */
  void
  write_pvtu_record(std::ostream &                  out,
                    const std::vector<std::string> &piece_names) const;

  /**
   * 这个函数并行地写入几个.vtu文件和一个.pvtu记录，并自动构建文件名。它是
   * DataOutInterface::write_vtu() 或
   * DataOutInterface::write_vtu_in_parallel(), 和
   * DataOutInterface::write_pvtu_record().
   * 的组合。例如，在10个进程中运行<code>
   * write_vtu_with_pvtu_record("output/", "solution", 3, comm, 4, 2)
   * </code>会生成这些文件
   * @code
   * output/solution_0003.0.vtu
   * output/solution_0003.1.vtu
   * output/solution_0003.pvtu
   * @endcode
   * 其中`.0.vtu`文件包含前一半进程的输出分组，而`.1.vtu`是其余一半进程的数据。
   * 一个指定的 @p directory 和一个 @p filename_without_extension
   * 构成文件名的第一部分。然后用 @p counter
   * 扩展文件名，标明当前的时间步数/迭代次数/等等，处理器ID，最后是.vtu/.pvtu结尾。由于要写入的时间步数取决于应用，在文件名中保留的数字可以作为参数
   * @p n_digits_for_counter,
   * 来指定，如果该参数保持默认值，则数字不加前导零
   * numbers::invalid_unsigned_int.
   * 如果需要一个以上的文件标识符（例如时间步数和求解器的迭代计数器），最后一个标识符作为
   * @p counter,
   * 使用，而所有其他标识符必须在调用该函数时加入 @p
   * filename_without_extension 。
   * 在并行设置中，每个时间步长通常要写几个文件。并行写入的文件数量取决于MPI进程的数量（见参数
   * @p mpi_communicator), 和默认值为0的指定数量 @p n_groups
   * 。其背景是VTU文件输出支持在并行文件系统上写入时，使用MPI
   * I/O将几个CPU的文件分组为给定数量的文件。 @p n_groups
   * 的默认值是0，意味着每个MPI等级将写入一个文件。1的值将生成一个包含整个域的解决方案的大文件，而更大的值将创建
   * @p n_groups 个文件（但不会超过MPI等级的数量）。
   * 请注意，只有一个处理器需要生成.pvtu文件，其中零号处理器被选择来承担这项工作。
   * 返回值是pvtu记录的集中文件的文件名。
   * @note  代码简单地结合了字符串 @p directory 和 @p
   * filename_without_extension, ，即用户必须确保 @p directory
   * 包含一个尾部字符，例如"/"，将目录和文件名分开。
   * @note
   * 如果要将输出写入当前工作目录，则使用空字符串""作为第一个参数。
   *
   */
  std::string
  write_vtu_with_pvtu_record(
    const std::string &directory,
    const std::string &filename_without_extension,
    const unsigned int counter,
    const MPI_Comm &   mpi_communicator,
    const unsigned int n_digits_for_counter = numbers::invalid_unsigned_int,
    const unsigned int n_groups             = 0) const;

  /**
   * 通过get_patches()获取数据，并以SVG格式写到<tt>out</tt>。参见
   * DataOutBase::write_svg. 。
   *
   */
  void
  write_svg(std::ostream &out) const;

  /**
   * 通过get_patches()获取数据，并以deal.II中间格式写到<tt>out</tt>中。参见
   * DataOutBase::write_deal_II_intermediate.
   * 注意，中间格式就像它的名字一样：内部数据的直接表示。它不是标准化的，每当我们改变内部表示时，它就会改变。你只能期望使用用于编写的相同版本的deal.II来处理以这种格式编写的文件。
   *
   */
  void
  write_deal_II_intermediate(std::ostream &out) const;

  /**
   * 基于data_filter中的数据创建一个XDMFEntry。这假设网格和求解数据被写到一个文件中。参见write_xdmf_file()中的使用实例。
   *
   */
  XDMFEntry
  create_xdmf_entry(const DataOutBase::DataOutFilter &data_filter,
                    const std::string &               h5_filename,
                    const double                      cur_time,
                    const MPI_Comm &                  comm) const;

  /**
   * 基于data_filter中的数据，创建一个XDMFEntry。这假设网格和解的数据被写入不同的文件。参见write_xdmf_file()中的使用实例。
   *
   */
  XDMFEntry
  create_xdmf_entry(const DataOutBase::DataOutFilter &data_filter,
                    const std::string &               h5_mesh_filename,
                    const std::string &               h5_solution_filename,
                    const double                      cur_time,
                    const MPI_Comm &                  comm) const;

  /**
   * 根据提供的XDMFEntry对象的向量，写一个XDMF文件。
   * 下面是一个如何用HDF5和DataOutFilter使用这个函数的例子。
   * @code
   * DataOutBase::DataOutFilterFlags flags(true, true);
   * DataOutBase::DataOutFilter data_filter(flags);
   * std::vector<XDMFEntry> xdmf_entries;
   * // Filter the data and store it in data_filter
   * data_out.write_filtered_data(data_filter);
   * // Write the filtered data to HDF5
   * data_out.write_hdf5_parallel(data_filter, "solution.h5", MPI_COMM_WORLD);
   * // Create an XDMF entry detailing the HDF5 file
   * auto new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
   *                                                "solution.h5",
   *                                                simulation_time,
   *                                                MPI_COMM_WORLD);
   * // Add the XDMF entry to the list
   * xdmf_entries.push_back(new_xdmf_entry);
   * // Create an XDMF file from all stored entries
   * data_out.write_xdmf_file(xdmf_entries, "solution.xdmf", MPI_COMM_WORLD);
   * @endcode
   *
   */
  void
  write_xdmf_file(const std::vector<XDMFEntry> &entries,
                  const std::string &           filename,
                  const MPI_Comm &              comm) const;

  /**
   * 将 @p data_filter
   * 中的数据写入一个单一的HDF5文件，包含网格和解的数值。下面是一个如何使用这个函数与DataOutFilter的例子。
   * @code
   * DataOutBase::DataOutFilterFlags flags(true, true);
   * DataOutBase::DataOutFilter data_filter(flags);
   * // Filter the data and store it in data_filter
   * data_out.write_filtered_data(data_filter);
   * // Write the filtered data to HDF5
   * data_out.write_hdf5_parallel(data_filter, "solution.h5", MPI_COMM_WORLD);
   * @endcode
   *
   *
   */
  void
  write_hdf5_parallel(const DataOutBase::DataOutFilter &data_filter,
                      const std::string &               filename,
                      const MPI_Comm &                  comm) const;

  /**
   * 将data_filter中的数据写到HDF5文件中。如果write_mesh_file为false，网格数据将不会被写入，而解文件将只包含解的数值。如果write_mesh_file为true，且文件名相同，则生成的文件将同时包含网格数据和求解值。
   *
   */
  void
  write_hdf5_parallel(const DataOutBase::DataOutFilter &data_filter,
                      const bool                        write_mesh_file,
                      const std::string &               mesh_filename,
                      const std::string &               solution_filename,
                      const MPI_Comm &                  comm) const;

  /**
   * DataOutFilter是一种中间数据格式，可以减少将被写入文件的数据量。这个函数所填充的对象随后可以再次用于写入具体文件格式的数据；例如，见
   * DataOutBase::write_hdf5_parallel(). 。
   *
   */
  void
  write_filtered_data(DataOutBase::DataOutFilter &filtered_data) const;


  /**
   * 根据给定的数据格式向<tt>out</tt>写入数据和网格。
   * 这个函数只是调用相应的<tt>write_*</tt>函数。如果没有要求输出格式，将写入<tt>default_format</tt>。
   * 如果没有提供格式，而默认格式是<tt>default_format</tt>，则会发生错误。
   *
   */
  void
  write(std::ostream &                  out,
        const DataOutBase::OutputFormat output_format =
          DataOutBase::default_format) const;

  /**
   * 设置默认格式。这里设置的值在任何时候都会被使用，要求输出格式为<tt>default_format</tt>。
   *
   */
  void
  set_default_format(const DataOutBase::OutputFormat default_format);


  /**
   * 设置用于输出的标志。这个方法希望<tt>flags</tt>是<tt>OutputFlagsBase</tt>的一个子类中的成员。
   *
   */
  template <typename FlagType>
  void
  set_flags(const FlagType &flags);


  /**
   * 一个函数，返回与基类中相应函数相同的字符串；唯一的例外是，如果省略了参数，则返回当前默认格式的值，即在调用此函数之前通过set_default_format()或parse_parameters()设置的格式的正确后缀。
   *
   */
  std::string
  default_suffix(const DataOutBase::OutputFormat output_format =
                   DataOutBase::default_format) const;

  /**
   * 通过在每个输出格式的参数文件中声明子段来声明所有输出格式的参数，并调用每个输出格式的标志类的相应<tt>declare_parameters</tt>函数。
   * 如果相应的格式不输出任何标志，那么某些声明的子段可能不包含条目。
   * 请注意，表示每个补丁的分区数量和输出格式的顶层参数没有被声明，因为它们只被传递给虚拟函数，而不被存储在这种类型的对象中。你必须自己声明它们。
   *
   */
  static void
  declare_parameters(ParameterHandler &prm);

  /**
   * 读取 declare_parameters()
   * 中声明的参数，并为输出格式设置相应的标志。
   * 这样得到的标志会覆盖之前所有默认构建的或由set_flags()函数设置的标志对象的内容。
   *
   */
  void
  parse_parameters(ParameterHandler &prm);

  /**
   * 返回这个对象的内存消耗估计值，单位是字节。
   * 这不是精确的（但通常会很接近），因为计算树的内存使用量（例如，
   * <tt>std::map</tt>) 是很困难的。
   *
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * 这是一个抽象函数，派生类通过这个函数将Patch结构形式的预处理数据（在基类DataOutBase中声明）传播给实际的输出函数。你需要重载这个函数，以便让输出函数知道它们应该打印什么。
   *
   */
  virtual const std::vector<DataOutBase::Patch<dim, spacedim>> &
  get_patches() const = 0;

  /**
   * 抽象的虚拟函数，基类的输出函数通过它获得数据集的名称。
   *
   */
  virtual std::vector<std::string>
  get_dataset_names() const = 0;

  /**
   * 该函数返回关于如何解释由一个以上数据集组成的输出文件的各个组成部分的信息。
   * 它返回一个索引对和相应的名称和类型的列表，表明输出的哪些成分被认为是矢量或张量值，而不仅仅是标量数据的集合。索引对是包括在内的；例如，如果我们有一个2d的斯托克斯问题，其分量是(u,v,p)，那么相应的矢量数据范围应该是(0,1)，返回的列表将只包括一个元组元素，如(0,1,
   * "速度",分量_是矢量的一部分)。
   * 由于一些派生类不知道非标量数据，这个函数有一个默认的实现，即简单地返回一个空字符串，意味着所有的数据都将被视为标量字段的集合。
   *
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const;

  /**
   * 验证get_dataset_names()和get_nonscalar_data_ranges()返回的数据集的名称是否有效。目前这包括检查名称是否被多次使用。如果遇到一个无效的状态，将在调试模式下触发一个Assert()。
   *
   */
  void
  validate_dataset_names() const;


  /**
   * 补丁的默认分区数。这是由parse_parameters()填充的，并且应该被派生类中的build_patches()遵守。
   *
   */
  unsigned int default_subdivisions;

private:
  /**
   * 标准的输出格式。
   * 如果输出格式default_format被要求，则使用此格式。它可以通过<tt>set_format</tt>函数或在参数文件中改变。
   *
   */
  DataOutBase::OutputFormat default_fmt;

  /**
   * OpenDX数据输出时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::DXFlags dx_flags;

  /**
   * 在输出UCD数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::UcdFlags ucd_flags;

  /**
   * 在输出GNUPLOT数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::GnuplotFlags gnuplot_flags;

  /**
   * 在输出POVRAY数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::PovrayFlags povray_flags;

  /**
   * 在一个空间维度上输出EPS数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::EpsFlags eps_flags;

  /**
   * 在一个空间维度上输出gmv数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::GmvFlags gmv_flags;

  /**
   * 在一个空间维度上输出Tecplot数据时使用的标志。可以通过<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::TecplotFlags tecplot_flags;

  /**
   * 在一个空间维度上输出vtk数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::VtkFlags vtk_flags;

  /**
   * 在一个空间维度上输出svg数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::SvgFlags svg_flags;

  /**
   * 在一个空间维度上输出deal.II中间数据时使用的标志。可以通过使用<tt>set_flags</tt>函数来改变。
   *
   */
  DataOutBase::Deal_II_IntermediateFlags deal_II_intermediate_flags;
};



/**
 * 一个用于读回以deal.II中间格式写入的数据的类，以便它能以任何其他支持的图形格式写入。这个类有两个主要用途。
 * 这个类的第一个用途是使应用程序可以推迟决定使用哪种图形格式，直到程序运行之后。数据以中间格式写入文件，以后可以将其转换为你希望的任何图形格式。这可能很有用，例如，如果你想把它转换为gnuplot格式以获得快速浏览，后来又想把它转换为OpenDX格式以获得高质量的数据版本。本类允许将这种中间格式读回程序中，并允许使用基类的相关函数将其写成任何其他支持的格式。
 * 第二种用途在并行程序中大多是有用的：与其让一个中心进程为整个程序生成图形输出，不如让每个进程为它所拥有的单元格生成图形数据，并以中间格式将其写入一个单独的文件中。稍后，所有这些中间文件可以被读回并合并在一起，这个过程与首先生成数据相比是很快的。使用中间格式主要是因为它允许单独的文件被合并，而一旦数据以任何支持的既定图形格式被写出来，这几乎是不可能的。
 * 这第二种使用情况在 step-18
 * 示例程序中做了一些详细解释。
 * 为了将数据读回这个对象，你必须知道写数据时使用的空间尺寸的模板参数。如果这种知识在编译时就可以得到，那么这就没有问题。然而，如果不是这样（比如在一个简单的格式转换器中），那么它就需要在运行时弄清楚，尽管编译器在编译时已经需要它。一个使用
 * DataOutBase::determine_intermediate_format_dimensions() 函数的方法。
 * 请注意，中间格式就像它的名字一样：内部数据的直接表示。它不是标准化的，每当我们改变内部表示时，它就会改变。你只能期望使用用于编写的同一版本的deal.II来处理以这种格式编写的文件。
 *
 *
 * @ingroup input output
 *
 *
 */
template <int dim, int spacedim = dim>
class DataOutReader : public DataOutInterface<dim, spacedim>
{
public:
  /**
   * 读取之前由 <tt>DataOutBase::write_deal_II_intermediate</tt>
   * 写入的补丁序列，并将其存储在当前对象中。这将覆盖之前的任何内容。
   *
   */
  void
  read(std::istream &in);

  /**
   * 这个函数可以用来将其他对象读取的补丁合并到本对象存储的补丁中。这有时很方便，例如，如果有一个领域分解算法，其中每个块由它自己的DoFHandler表示，但人们想同时输出所有块上的解决方案。另外，它也可用于并行程序，即每个进程只对其所占的单元产生输出，即使所有进程都能看到所有单元。
   * 为了使其发挥作用，本对象的输入文件和给定参数需要有相同数量的输出向量，并且它们需要在每个补丁中使用相同数量的细分。如果两个对象中的补丁在空间上重叠，输出结果可能会看起来相当有趣。
   * 如果你在合并补丁后为这个对象调用read()，之前的状态会被覆盖，合并的补丁会丢失。
   * 如果这个对象或另一个对象还没有设置任何补丁，这个函数将会失败。
   *
   */
  void
  merge(const DataOutReader<dim, spacedim> &other);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcIncompatibleDatasetNames,
                   "You are trying to merge two sets of patches for which the "
                   "declared names of the variables do not match.");
  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcIncompatiblePatchLists,
                   "You are trying to merge two sets of patches for which the "
                   "number of subdivisions or the number of vector components "
                   "do not match.");
  /**
   * 异常情况
   *
   */
  DeclException4(ExcIncompatibleDimensions,
                 int,
                 int,
                 int,
                 int,
                 << "Either the dimensions <" << arg1 << "> and <" << arg2
                 << "> or the space dimensions <" << arg3 << "> and <" << arg4
                 << "> do not match!");

protected:
  /**
   * 这是一个函数，该类通过该函数将预处理的数据以补丁结构（在基类DataOutBase中声明）的形式传播给实际的输出函数。
   * 它返回上次给read()函数的流所读取的补丁。
   *
   */
  virtual const std::vector<dealii::DataOutBase::Patch<dim, spacedim>> &
  get_patches() const override;

  /**
   * 抽象的虚拟函数，基类的输出函数通过它获得数据集的名称。
   * 返回上次读取文件时读取的变量名称。
   *
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * 该函数返回关于如何解释由一个以上数据集组成的输出文件的各个组成部分的信息。
   * 它返回一个索引对和相应的名称的列表，表明输出的哪些成分被认为是矢量值的，而不仅仅是标量数据的集合。索引对是包括在内的；例如，如果我们有一个2d的斯托克斯问题，其分量是(u,v,p)，那么相应的矢量数据范围应该是(0,1)，返回的列表将只包括一个元组的元素，如(0,1,
   * "速度") 。
   * 由于一些派生类不知道矢量数据，这个函数有一个默认的实现，只是返回一个空字符串，这意味着所有的数据都将被视为标量字段的集合。
   *
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const override;

private:
  /**
   * 保存补丁集以及输出变量名称的数组，所有这些都是我们从输入流中读取的。
   *
   */
  std::vector<dealii::DataOutBase::Patch<dim, spacedim>> patches;
  std::vector<std::string>                               dataset_names;

  /**
   * 关于输出域的某些成分是否要被视为矢量的信息。
   *
   */
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    nonscalar_data_ranges;
};



/**
 * 一个用于存储相关数据的类，在编写轻量级XDMF文件时使用。XDMF文件反过来指向存储实际模拟数据的重型数据文件（如HDF5）。这允许灵活安排数据，也允许将网格与点数据分开。
 *
 *
 */
class XDMFEntry
{
public:
  /**
   * 默认的构造函数，可以创建一个无效的对象。
   *
   */
  XDMFEntry();

  /**
   * 简化的构造函数，在  <code>solution_filename ==
   * mesh_filename</code>  , 和  <code>dim==spacedim</code>
   * 的情况下调用完整的构造函数。
   *
   */
  XDMFEntry(const std::string &filename,
            const double       time,
            const unsigned int nodes,
            const unsigned int cells,
            const unsigned int dim);

  /**
   * 简化的构造函数，在  <code>dim==spacedim</code>
   * 的情况下调用完整的构造函数。
   *
   */
  XDMFEntry(const std::string &mesh_filename,
            const std::string &solution_filename,
            const double       time,
            const unsigned int nodes,
            const unsigned int cells,
            const unsigned int dim);

  /**
   * 将所有成员设置为所提供的参数的构造函数。
   *
   */
  XDMFEntry(const std::string &mesh_filename,
            const std::string &solution_filename,
            const double       time,
            const unsigned int nodes,
            const unsigned int cells,
            const unsigned int dim,
            const unsigned int spacedim);

  /**
   * 记录一个属性和相关维度。
   *
   */
  void
  add_attribute(const std::string &attr_name, const unsigned int dimension);

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)读取或写入此对象的数据进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int  /*version*/ )
  {
    ar &valid &h5_sol_filename &h5_mesh_filename &entry_time &num_nodes
      &num_cells &dimension &space_dimension &attribute_dims;
  }

  /**
   * 获取与此条目相关的XDMF内容。
   * 如果该条目无效，则返回一个空字符串。      @deprecated
   * 使用重载，取一个`无符号int`和一个`const ReferenceCell
   * &`来代替。
   *
   */
  DEAL_II_DEPRECATED
  std::string
  get_xdmf_content(const unsigned int indent_level) const;

  /**
   * 获取与该条目相关的XDMF内容。
   * 如果该条目无效，则返回一个空字符串。
   *
   */
  std::string
  get_xdmf_content(const unsigned int   indent_level,
                   const ReferenceCell &reference_cell) const;

private:
  /**
   * 该条目是否有效并包含要写入的数据。
   *
   */
  bool valid;

  /**
   * 该条目引用的HDF5重数据解决方案文件的名称。
   *
   */
  std::string h5_sol_filename;

  /**
   * 此条目引用的HDF5网格文件的名称。
   *
   */
  std::string h5_mesh_filename;

  /**
   * 与此条目相关的模拟时间。
   *
   */
  double entry_time;

  /**
   * 数据节点的数量。
   *
   */
  unsigned int num_nodes;

  /**
   * 数据单元的数量。
   *
   */
  unsigned int num_cells;

  /**
   * 与数据相关的维度。
   *
   */
  unsigned int dimension;

  /**
   * 数据所处空间的维度。  请注意，维度<=空间维度。
   *
   */
  unsigned int space_dimension;

  /**
   * 与此条目相关的属性和它们的维度。
   *
   */
  std::map<std::string, unsigned int> attribute_dims;
};



 /* -------------------- inline functions ------------------- */ 

namespace DataOutBase
{
  inline bool
  EpsFlags::RgbValues::is_grey() const
  {
    return (red == green) && (red == blue);
  }


   /* -------------------- template functions ------------------- */ 

  /**
   * <tt>DataOutBase::Patch</tt>. 类型对象的输出操作符
   * 该操作符转储由补丁数据结构代表的中间图形格式。它以后可以被转换为一些图形程序的常规格式。
   *
   */
  template <int dim, int spacedim>
  std::ostream &
  operator<<(std::ostream &out, const Patch<dim, spacedim> &patch);



  /**
   * 类型对象的输入操作符  <tt>DataOutBase::Patch</tt>.
   * 该操作符读取由补丁数据结构代表的中间图形格式，使用操作符<<的格式写入。
   *
   */
  template <int dim, int spacedim>
  std::istream &
  operator>>(std::istream &in, Patch<dim, spacedim> &patch);
} // namespace DataOutBase


DEAL_II_NAMESPACE_CLOSE

#endif


