//include/deal.II-translator/grid/tria_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_tria_h
#define dealii_tria_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/point.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/tria_description.h>
#include <deal.II/grid/tria_iterator_selector.h>
#include <deal.II/grid/tria_levels.h>

#include <boost/serialization/map.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/unique_ptr.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/signals2.hpp>

#include <bitset>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <numeric>
#include <vector>


DEAL_II_NAMESPACE_OPEN

#ifdef signals
#  error \
    "The name 'signals' is already defined. You are most likely using the QT library \
and using the 'signals' keyword. You can either #include the Qt headers (or any conflicting headers) \
*after* the deal.II headers or you can define the 'QT_NO_KEYWORDS' macro and use the 'Q_SIGNALS' macro."
#endif

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Manifold;

template <int dim>
struct CellData;

struct SubCellData;

namespace TriangulationDescription
{
  template <int, int>
  struct Description;
}

namespace GridTools
{
  template <typename CellIterator>
  struct PeriodicFacePair;
}

template <int, int, int>
class TriaAccessor;
template <int spacedim>
class TriaAccessor<0, 1, spacedim>;
template <int, int, int>
class TriaAccessorBase;

namespace internal
{
  namespace TriangulationImplementation
  {
    class TriaFaces;

    class TriaObjects;

    template <int, int>
    class Policy;

    /**
     * 一个类的前向声明，我们把三角形类的大部分实现放在这个类中。参见.cc文件以了解更多信息。
     *
     */
    struct Implementation;
    struct ImplementationMixedMesh;
  } // namespace TriangulationImplementation

  namespace TriaAccessorImplementation
  {
    struct Implementation;
  }
} // namespace internal
#endif


 /*------------------------------------------------------------------------*/ 


namespace internal
{
  /**
   * 一个命名空间，用于三角化类和帮助器的内部类。
   *
   */
  namespace TriangulationImplementation
  {
    /**
     * 缓存类，用于存储三角剖分中各层内已使用和活动的元素（线或四边形等）的数量。这只是模板的声明，具体实例在下面。
     * 在过去，每当人们想要访问这些数字之一时，就必须在所有的线上执行一个循环，例如，计算元素，直到我们碰到终端迭代器。这很耗时，而且由于访问行数等是一个相当频繁的操作，这并不是一个最佳的解决方案。
     *
     */
    template <int dim>
    struct NumberCache
    {};

    /**
     * 缓存类，用于存储三角形的各层中已使用和活跃的元素（线或四边形等）的数量。这个特殊化存储了线的数量。
     * 在过去，每当人们想要访问这些数字之一时，就必须在所有的线上进行循环，例如，计算元素，直到我们碰到结束迭代器。这很耗时，而且由于访问行数等是一个相当频繁的操作，这并不是一个最佳的解决方案。
     *
     */
    template <>
    struct NumberCache<1>
    {
      /**
       * 我们使用过的对象的层数。
       *
       */
      unsigned int n_levels;

      /**
       * 整个三角测量中使用的线的数量。
       *
       */
      unsigned int n_lines;

      /**
       * 保存每层所用线数的数组。
       *
       */
      std::vector<unsigned int> n_lines_level;

      /**
       * 整个三角测量中的活动线数。
       *
       */
      unsigned int n_active_lines;

      /**
       * 保存每层活动线数的数组。
       *
       */
      std::vector<unsigned int> n_active_lines_level;

      /**
       * 构造函数。默认情况下，将数值设置为零。
       *
       */
      NumberCache();

      /**
       * 确定此对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };


    /**
     * 缓存类，用于存储三角形的各层中已使用和活动的元素（线或四边形等）的数量。这个特殊化存储四边形的数量。由于从基类NumberCache<1>的继承，线的数量也在这个类中。
     * 在过去，每当人们想要访问这些数字中的一个，就必须在所有的线上进行循环，例如，计算元素，直到我们碰到结束的迭代器。这很耗时，而且由于访问行数等是一个相当频繁的操作，这并不是一个最佳的解决方案。
     *
     */
    template <>
    struct NumberCache<2> : public NumberCache<1>
    {
      /**
       * 整个三角形中使用的四边形的数量。
       *
       */
      unsigned int n_quads;

      /**
       * 保存每层所用四边形数量的数组。
       *
       */
      std::vector<unsigned int> n_quads_level;

      /**
       * 整个三角结构中的活动四边形的数量。
       *
       */
      unsigned int n_active_quads;

      /**
       * 保存每层活动四边形数量的数组。
       *
       */
      std::vector<unsigned int> n_active_quads_level;

      /**
       * 构造函数。默认情况下，将数值设置为零。
       *
       */
      NumberCache();

      /**
       * 确定此对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };


    /**
     * 缓存类，用于存储三角图各层中已使用和活动的元素（线或四边形等）的数量。这个特殊化存储的是六边形的数量。由于从基类NumberCache<2>的继承，线和四边形的数量也在这个类中。
     * 在过去，每当人们想要访问这些数字之一时，就必须在所有的线上执行一个循环，例如，计算元素，直到我们碰到终点。这很耗时，而且由于访问行数等是一个相当频繁的操作，这并不是一个最佳的解决方案。
     *
     */
    template <>
    struct NumberCache<3> : public NumberCache<2>
    {
      /**
       * 整个三角测量中使用的六边形的数量。
       *
       */
      unsigned int n_hexes;

      /**
       * 保存每层使用的六边形数量的数组。
       *
       */
      std::vector<unsigned int> n_hexes_level;

      /**
       * 整个三角测量中的活动爻数。
       *
       */
      unsigned int n_active_hexes;

      /**
       * 保存每层活动六角的数量的数组。
       *
       */
      std::vector<unsigned int> n_active_hexes_level;

      /**
       * 构造函数。默认情况下，将数值设置为零。
       *
       */
      NumberCache();

      /**
       * 确定此对象的内存消耗（以字节为单位）的估计值。
       *
       */
      std::size_t
      memory_consumption() const;

      /**
       * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
       *
       */
      template <class Archive>
      void
      serialize(Archive &ar, const unsigned int version);
    };
  } // namespace TriangulationImplementation
} // namespace internal


 /*------------------------------------------------------------------------*/ 


/**
 * 三角形是一个单元的集合，这些单元共同覆盖了人们通常想要解决的偏微分方程的领域。这个域和覆盖它的网格代表了一个
 * @p dim 。
 *
 * - 扩张流形，并且生活在 @p spacedim 空间维度中，其中 @p dim 和 @p spacedim 是该类的模板参数。(如果没有指定 @p spacedim ，则采用默认值`spacedim=dim`)。
 * 因此，例如，一个 @p Triangulation<1,1>
 * 类型的对象（或者干脆是 @p  Triangulation<1>，因为默认为 @p
 * spacedim==dim
 * ）被用来表示和处理有限元方法中常用的一维三角形（因此，直线上的段）。另一方面，像
 * @p Triangulation<1,2> 或 @p Triangulation<2,3>
 * 这样的对象（与二维的曲线或三维的曲面有关）是人们想在边界元素方法中使用的对象。
 * 该类的名称主要是分层次的，并不意味着三角计算只能由三角形组成。相反，三角形由1d的线段组成（即，如果`dim==1`），以及由三维单元组成（如果`dim==3`）。此外，历史上，deal.II只支持二维的四边形（有四个顶点的单元：变形的矩形）和六面体（有六个边和八个顶点的单元，是变形的盒子），它们都不是三角形。换句话说，deal.II语言中的术语
 * "三角形 "是 "网格
 * "的同义词，应与它的语言来源分开理解。
 * 这个类被写成尽可能独立于维度（因此
 * dealii::internal::TriangulationImplementation::TriaLevel
 * 类的复杂结构），以允许代码共享，允许减少将一个维度的代码变化反映到其他维度的代码中的需要。尽管如此，一些函数是依赖于维度的，并且只存在针对不同维度的专门版本。
 * 这个类满足了 @ref ConceptMeshType "MeshType概念 "
 * 的要求。 <h3>Structure and iterators</h3>
 * Triangulation对象的实际数据结构是相当复杂的，如果试图直接对其进行操作的话，是相当不方便的，因为数据分布在相当多的数组和其他地方。然而，有足够强大的方法可以在不知道其确切关系的情况下对这些数据结构进行操作。deal.II使用类的局部别名（见下文）来使事情变得尽可能简单和不依赖维度。
 * Triangulation类提供了迭代器，可以在不知道用于描述单元的确切表示法的情况下，在所有单元上循环操作。更多信息见<tt>TriaIterator</tt>的文档。它们的名字是从迭代器类中导入的别名（从而使它们成为这个类的本地类型），具体如下。
 * <ul>   <li>  <tt>cell_iterator</tt>: 循环处理三角测量中使用的所有单元  <li>  <tt>active_cell_iterator</tt>: 循环处理所有活动单元  </ul>  。
 * 对于<tt>dim==1</tt>，这些迭代器被映射为如下。
 * @code
 *  using cell_iterator = line_iterator;
 *  using active_cell_iterator = active_line_iterator;
 * @endcode
 * 而对于 @p dim==2 我们有额外的面孔迭代器。
 * @code
 *  using cell_iterator = quad_iterator;
 *  using active_cell_iterator = active_quad_iterator;
 *
 *  using face_iterator = line_iterator;
 *  using active_face_iterator = active_line_iterator;
 * @endcode
 *
 * 通过使用单元格迭代器，你可以编写独立于空间维度的代码。这同样适用于子结构迭代器，子结构被定义为一个单元的面。单元的面在一维是一个顶点，在二维是一条线；但是，顶点的处理方式不同，因此线没有面。
 * Triangulation类提供了一些函数，如begin_active()，它给你一个通往第一个活动单元的迭代器。有相当多的函数返回迭代器。请看一下类的文档以获得一个概述。
 * 这些迭代器的使用与标准容器迭代器的使用类似。以下是Triangulation源代码中的一些例子（注意在最后两个例子中，模板参数
 * @p spacedim 被省略了，所以它采用默认值 <code>dim</code> ）。
 * <ul>   <li>   <em>  计算特定层次上的细胞数量  </em>  。
 * @code
 *    template <int dim, int spacedim>
 *    unsigned int
 *    Triangulation<dim, spacedim>::n_cells (const int level) const
 *    {
 *      int n=0;
 *      for (const auto &cell : cell_iterators_on_level(level))
 *        ++n;
 *      return n;
 *    }
 * @endcode
 * 另一种方法，即使用 <tt>std::distance</tt>, 将写成
 * @code
 *    template <int dim>
 *    unsigned int
 *    Triangulation<dim>::n_cells (const int level) const
 *    {
 *      int n=0;
 *      distance (begin(level),
 *                (level == levels.size()-1 ?
 *                 cell_iterator(end()) :
 *                 begin (level+1)),
 *                n);
 *      return n;
 *    }
 * @endcode
 * <li>   <em>  完善一个三角形的所有单元  </em>  。
 * @code
 *    template <int dim>
 *    void Triangulation<dim>::refine_global ()
 *    {
 *      for (const auto &cell : active_cell_iterators())
 *        cell->set_refine_flag ();
 *      execute_coarsening_and_refinement ();
 *    }
 * @endcode
 * </ul>
 *
 *  <h3>Usage</h3>
 * 三角形的使用主要是通过使用迭代器完成的。一个例子可能最能说明如何使用它。
 *
 * @code
 * int main ()
 * {
 * Triangulation<2> tria;
 *
 * // read in a coarse grid file
 *
 * // we want to log the refinement history
 * ofstream history ("mesh.history");
 *
 * // refine first cell
 * tria.begin_active()->set_refine_flag();
 * tria.save_refine_flags (history);
 * tria.execute_coarsening_and_refinement ();
 *
 * // refine first active cell on coarsest level
 * tria.begin_active()->set_refine_flag ();
 * tria.save_refine_flags (history);
 * tria.execute_coarsening_and_refinement ();
 *
 * Triangulation<2>::active_cell_iterator cell;
 * for (int i=0; i<17; ++i)
 *   {
 *     // refine the presently second last cell 17 times
 *     cell = tria.last_active(tria.n_levels()-1);
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * --cell;
 *     cell->set_refine_flag ();
 *     tria.save_refine_flags (history);
 *     tria.execute_coarsening_and_refinement ();
 *   };
 * // output the grid
 * ofstream out("grid.1");
 * GridOut::write_gnuplot (tria, out);
 * }
 * @endcode
 *
 *
 *  <h3>Creating a triangulation</h3>
 * 有几种创建三角的可能性。  <ul>   <li>  最常见的域，如超立方体（即直线、正方形、立方体等）、超球体（圆、球...）和其他一些更奇怪的域，如L形区域和高维泛化等，由GridGenerator类提供，它接受一个三角形，通过对所需域的划分来填充它。
 * <li>
 * 读入一个三角剖面。通过使用GridIn类的一个对象，你可以读入相当普遍的三角形。更多信息见那里。所提到的类使用下面直接描述的接口，将数据传输到三角测量中。
 * <li>
 * 明确创建一个三角形：你可以通过提供一个顶点列表和一个单元格列表来创建一个三角形。每个单元由一个向量组成，存储该单元在顶点列表中的顶点的索引。要看这是如何工作的，你可以看看
 * GridIn<dim>::read_*
 * 函数。要调用的适当函数是create_triangulation()。
 * 从只存储顶点信息的单元中创建本库所需的层次信息可能是一项相当复杂的任务。
 * 例如在二维中，我们必须在顶点之间创建线条（但只有一次，尽管有两个单元格将这两个顶点连接起来），我们还必须创建邻域信息。因此，被读入的网格不应该太大，读入细化的网格将是低效的（尽管从技术上讲，读入几个10,000或100,000个单元的网格是没有问题的；库可以处理这些问题）。除了性能方面，细化网格并不适合多网格算法，因为在最粗的层次上求解是昂贵的。在任何情况下，最明智的做法是尽可能粗略地读入网格，然后进行必要的细化步骤。
 * 你有责任保证单元有正确的方向。为了保证这一点，在保持单元格列表的输入向量中，每个单元格的顶点指数必须按规定的顺序排列，见GeometryInfo<dim>的文档。在一维中，第一个顶点索引必须是指具有较低坐标值的顶点。在二维和三维中，相应的条件不容易验证，也没有完全尝试这样做。如果你违反了这个条件，你可能会发现矩阵条目有错误的符号（顺时针方向的顶点编号，这将导致一个负的面积元素）或错误的矩阵元素（扭曲的四边形，即两个顶点互换；这将导致一个错误的面积元素）。
 * 有一些更微妙的条件必须施加在单元格内的顶点编号上。它们不仅适用于从UCD或任何其他输入文件中读取的数据，也适用于传递给create_triangulation()的数据。有关这方面的更多细节，请参见GridIn类的文档，最重要的是GridReordering类，它解释了许多问题和一种重新排序单元的算法，使它们满足上述条件。
 * <li>
 * 复制三角形：当在与时间相关的网格上计算或使用自适应细化时，你经常希望创建一个新的三角形，使其与另一个相同。
 * @p copy_triangulation函数为此提供了便利。
 * 它保证了两个三角形中的顶点、线或单元格的编号是相同的，并且如果两个迭代器在两个三角形上行走，它们会访问匹配的单元格，如果它们是平行递增的。可以设想在复制操作中实现清理，消除未使用的内存孔，重新连接分散的数据等等。原则上，这将是一个有用的操作，但是保证两个三角形的某些并行性似乎更重要，因为通常数据必须在网格之间传输。  </ul>
 * 最后，对于喜欢坏网格的人来说，还有一个特殊的函数：distort_random()。它将网格中的所有顶点按一个随机值左右移动，留下一个扭曲的网格。注意，你应该将这个函数应用于最终的网格，因为细化会使网格变得更加平滑。
 * 该函数将确保受限面的顶点（悬空节点）最终会出现在正确的位置，即在母线的另外两个顶点的中间，在更高的空间维度上也是如此（边界上的顶点不会被修正，所以不要在两个以上的空间维度上扭曲边界顶点，即在边界顶点可以成为悬空节点的维度上）。然而，应用该算法还有一个与单元格放置有关的缺点：一个单元格的子女不会像母单元格那样占据领域的同一区域。虽然这是边界上的单元的通常行为，但在使用多网格算法或将解从粗网格转移到细网格时，你可能会遇到麻烦。一般来说，只有当你只使用最精细的三角网格进行计算时，使用这个函数才是安全的。
 *
 *
 * <h3>Refinement and coarsening of a triangulation</h3>
 * 三角形的细化可以通过几种方式完成。最低级的方式是直接通过迭代器：让
 * @p i
 * 成为一个活动单元的迭代器（即指向的单元没有子代），然后函数调用<tt>i->set_refine_flag()</tt>标记相应的单元进行细化。标记非活动单元格会导致一个错误。
 * 在所有你想标记为细化的单元格之后，调用execute_coarsening_and_refinement()来实际执行细化。这个函数本身首先调用
 * @p prepare_coarsening_and_refinement
 * 函数来规范生成的三角形：由于两个相邻单元之间的面只能被细分一次（即两个相邻单元的层次最多只能相差一个；不可能一个单元被细化两次而相邻的单元没有被细化），一些额外的单元被标记为细化以平滑网格。这增加了结果单元的数量，但使网格更加规则，从而导致更好的近似特性，最重要的是，使数据结构和算法的处理更加容易。说实话，这主要是一个算法上的步骤，而不是有限元方法所需要的。
 * 要粗化一个网格，可以通过使用<tt>i->set_coarsen_flag</tt>和调用execute_coarsening_and_refinement()来实现上述同样的方法。
 * 先粗化，后细化的原因是，细化通常会增加一些额外的单元以保持三角形的规则，从而满足所有细化的要求，而粗化不会删除没有要求的单元；因此细化往往会恢复粗化的一些效果，而反之则不然。因此，所述的先粗化后细化的顺序通常会导致一个更接近预期的结果。
 * 通过迭代器 "手工
 * "标记单元进行细化是产生新网格的一种方法，特别是当你知道你在寻找什么样的网格时，比如你想让网格向边界连续细化或者总是在中心细化（参见示例程序，它们正是做这些事情）。然而，还有一些更高级的函数，它们更适合于在后验误差估计和自适应有限元的背景下自动生成分层网格。这些函数可以在GridRefinement类中找到。
 *
 *  <h3>Smoothing of a triangulation</h3>
 * 对于过于非结构化的网格，已经观察到一些近似特性的退化。因此，prepare_coarsening_and_refinement()被execute_coarsening_and_refinement()自动调用，可以对三角网格进行一些平滑处理。注意，网格平滑只对两个或多个空间维度进行，目前没有对一个空间维度的平滑。在下文中，让<tt>execute_*</tt>代表execute_coarsening_and_refinement（）。
 * 为了实现平滑化，Triangulation构造函数需要一个参数，指定每次调用<tt>execute_*</tt>时是否要对网格进行平滑化处理。默认情况下是不做这样的步骤，因为这将导致产生额外的单元格，这在所有情况下可能是不必要的。如果开启，调用<tt>execute_*</tt>的结果是标记额外的单元格进行细化，以避免出现上述的顶点。正则化和平滑三角形的算法将在下面的技术问题部分描述。这个参数必须给构造函数而不是给<tt>execute_*</tt>的原因是，如果你调用<tt>execute_*</tt>一次，没有平滑，一次就会导致算法问题，因为这样在某些细化步骤中需要细化两次。
 * 构造函数获取的参数是一个整数，它可以由定义在枚举#MeshSmoothing中的常数进行比特化组合（见那里的可能性）。
 *
 *
 * @note  虽然有可能将#MeshSmoothing中的所有标志传递给类型为
 * parallel::distributed::Triangulation,
 * 的对象，但如果它们需要了解不属于这个处理器的单元上的细化/粗化标志，则并不总是能够实现所有这些平滑选项。因此，对于其中的一些标志，并行三角形的最终单元数可能取决于它被分割成的处理器的数量。
 *
 *  <h3>Material and boundary information</h3>
 * 每个单元、面或边都存储了表示物体所属的材料或边界部分的信息。单元的材料ID通常用于识别哪些单元属于域的特定部分，例如，当你有不同的材料（钢铁、混凝土、木材）都属于同一个域时。然后，在组装双线性表格时，通常会查询与某个单元相关的材料ID，并使用它来确定（例如，通过表格查询，或一连串的if-else语句）该单元的正确材料系数是什么。另见 @ref GlossMaterialId "本词汇表条目"
 * 。 这个 material_id 可以在构建三角形时设置（通过 CellData
 * 数据结构），也可以在之后通过使用单元格迭代器设置。关于这个功能的典型使用，请看
 * step-28
 * 的教程程序。GridGenerator命名空间的函数通常将所有单元的材料ID设置为0。当通过GridIn类读取三角图时，不同的输入文件格式有不同的约定，但通常是明确指定材料ID，如果没有，则GridIn简单地将其设置为零。因为一个单元的材料是与域的特定区域相关的，所以材料ID在网格细化时由子单元从其父单元继承。
 * 低维对象上的边界指示器（这些对象没有材料ID）表示边界组件的数量。偏微分方程的弱表述可能在边界的不同部分有不同的边界条件。边界指标可以在创建矩阵或右侧向量时使用，以表示模型的这些不同部分（这种使用就像单元格的材料id）。边界指示器的范围可以从零到
 * numbers::internal_face_boundary_id-1. ，值
 * numbers::internal_face_boundary_id
 * 是保留的，用来表示没有边界指示器的内部线（在二维）和内部线和四边形（在三维）。这样一来，程序就可以很容易地确定这样的物体是否在边界上。材料指标的范围可以从零到
 * numbers::invalid_material_id-1. 。
 * 二维的线和三维的四边形在细化时将其边界指标继承给它们的子代。因此，你应该确保如果你有不同的边界部分，不同的部分被一个顶点（在二维）或一条线（在三维）分开，这样每个边界线或四边形都有一个唯一的边界指标。
 * 默认情况下（除非在创建三角图的过程中另有规定），边界的所有部分都有边界指标为零。作为一个历史遗留问题，这对于1d网格来说并不是真的。对于这些，最左边的顶点的边界指标为零，而最右边的顶点的边界指标为一。在这两种情况下，一个面的边界指示器都可以通过调用
 * <code>cell-@>face(1)-@>set_boundary_id(42);</code>  来改变。
 * @see   @ref GlossBoundaryIndicator  "关于边界指示器的词汇条目"
 *
 *  <h3>History of a triangulation</h3>
 * 可以从细化历史中重建网格，这些历史可以通过 @p
 * save_refine_flags  和  @p
 * load_refine_flags函数来存储和加载。通常情况下，代码会是这样的。
 * @code
 *   // open output file
 *   std::ofstream history("mesh.history");
 *   // do 10 refinement steps
 *   for (unsigned int step=0; step<10; ++step)
 *     {
 *       ...;
 *       // flag cells according to some criterion
 *       ...;
 *       tria.save_refine_flags (history);
 *       tria.execute_coarsening_and_refinement ();
 *     }
 * @endcode
 *
 * 如果你想从存储的信息中重新创建网格，你就写。
 * @code
 *   // open input file
 *   std::ifstream history("mesh.history");
 *   // do 10 refinement steps
 *   for (unsigned int step=0; step<10; ++step)
 *     {
 *       tria.load_refine_flags (history);
 *       tria.execute_coarsening_and_refinement ();
 *     }
 * @endcode
 *
 * 粗化和粗化标志也采用同样的方案。
 * 你可以在不同的细化信息集之间向输出文件写入其他信息，只要你在重新创建网格时读取这些信息。你应该确保将从保存的标志中创建的新三角图中的其他信息与旧三角图的信息相匹配，例如平滑水平；如果不是，从标志中实际创建的单元可能是其他单元，因为平滑增加了额外的单元，但它们的数量可能取决于平滑水平。
 * 实际上有两组<tt>save_*_flags</tt>和<tt>load_*_flags</tt>函数。一个以流为参数，从/到流中读/写信息，从而实现将标志存储到文件。另一个集合需要一个<tt>vector<bool></tt>类型的参数。这使得用户可以临时存储一些标志，例如，如果另一个函数需要它们，之后再恢复它们。
 *
 *  <h3>User flags and data</h3>
 * 三角形为用户标志提供了每行、四边形等的一个比特。这个字段可以像所有其他数据一样使用迭代器来访问。通常情况下，如果一个算法走过所有的单元并需要另一个单元，例如邻居，是否已经被处理过，那么这个用户标志就会被使用。参见 @ref GlossUserFlags "更多信息的词汇表"
 * 。
 * 还有一组用户数据，可以是<tt>无符号int</tt>或<tt>void</tt>，用于每一行、四边形等。你可以通过访问器类中<tt>用户数据</tt>下所列的函数访问这些数据。同样，见 @ref GlossUserData "词汇表的更多信息"
 * 。 这些用户索引或指针的值默认为 @p nullptr
 * 。请注意，这些指针在细化时不会被继承给子代。然而，在重新划分后，它们在所有单元上都是可用的，因为它们在以前的网格上被设置过。
 * 通常关于 @p void
 * 指针的类型安全缺失的警告在这里显然是适用的；类型的正确性等的责任完全在于指针的使用者。
 *
 *
 * @note
 * 用户指针和用户索引被存储在同一个地方。为了避免不必要的转换，Triangulation检查其中一个正在使用，并且不允许访问另一个，直到调用clear_user_data()。
 *
 *  <h3>Describing curved geometries</h3>
 * deal.II用继承自Manifold的类来实现所有的几何图形（弯曲的和其他的）；参见Manifold的文档， step-49
 * ，或 @ref manifold
 * 模块的例子和算法的完整描述。默认情况下，Triangulation中的所有单元都具有平坦的几何形状，也就是说，Triangulation中的所有线条都被假定为直线。如果一个单元的manifold_id不等于
 * numbers::flat_manifold_id
 * ，那么Triangulation会使用相关的Manifold对象对该单元进行计算（如单元细化）。下面是一个快速的例子，取自
 * GridGenerator::hyper_ball(), 的实现，它设置了一个极地网格。
 *
 *
 * @code
 * int main ()
 * {
 * Triangulation<2> triangulation;
 * const std::vector<Point<2>> vertices = {{-1.0,-1.0},
 *                                         {+1.0,-1.0},
 *                                         {-0.5,-0.5},
 *                                         {+0.5,-0.5},
 *                                         {-0.5,+0.5},
 *                                         {+1.0,+1.0},
 *                                         {-1.0,+1.0},
 *                                         {+1.0,+1.0}};
 * const std::vector<std::array<int,GeometryInfo<2>::vertices_per_cell>>
 *   cell_vertices = {{0, 1, 2, 3},
 *                    {0, 2, 6, 4},
 *                    {2, 3, 4, 5},
 *                    {1, 7, 3, 5},
 *                    {6, 4, 7, 5}};
 *
 * std::vector<CellData<2>> cells(cell_vertices.size(), CellData<2>());
 * for (unsigned int i=0; i<cell_vertices.size(); ++i)
 *   for (unsigned int j=0; j<GeometryInfo<2>::vertices_per_cell; ++j)
 *     cells[i].vertices[j] = cell_vertices[i][j];
 *
 * triangulation.create_triangulation (vertices, cells, SubCellData());
 * triangulation.set_all_manifold_ids_on_boundary(42);
 *
 * // set_manifold stores a copy of its second argument,
 * // so a temporary is okay
 * triangulation.set_manifold(42, PolarManifold<2>());
 * for (unsigned int i = 0; i < 4; ++i)
 *   {
 *     // refine all boundary cells
 *     for (const auto &cell : triangulation.active_cell_iterators())
 *       if (cell->at_boundary())
 *         cell->set_refine_flag();
 *
 *     triangulation.execute_coarsening_and_refinement();
 *   }
 * }
 * @endcode
 *
 * 这将设置一个网格，边界线将通过在极坐标中进行计算而被细化。当网格被细化时，与边界相邻的单元将使用这个新的线中点（以及其他三个中点和原来的单元顶点），用转折插值计算单元中点：这将使弯曲的边界以平滑的方式传播到内部。通过使用TransfiniteInterpolationManifold，可以生成一个更好的网格（在两个不同的Manifold描述之间对所有单元进行内插，而不是每次只插一个单元）；更多信息请参见该类的文档。
 * 你应该注意一个注意事项：如果你有凹形边界，你必须确保一个新的边界顶点不会位于要被细化的单元内太多。原因是中心顶点被放置在原单元的顶点、新面的中点和（在三维）新线的中点的加权平均点。因此，如果你的新边界顶点太靠近旧的四边形或六面体的中心，到中点顶点的距离就会变得太小，从而产生扭曲的单元。这个问题在 @ref GlossDistorted "扭曲的单元格 "
 * 中得到了广泛的讨论。 <h3>Getting notice when a triangulation
 * changes</h3>
 * 在有些情况下，一个对象希望知道每当一个三角形被细化、复制或以其他一些方式被修改时。当然，如果在你的用户代码中，每当你要细化三角测量时，你都要告诉每一个这样的对象，这可以实现，但这将变得很繁琐，而且容易出错。Triangulation类实现了一种更优雅的方式来实现这一点：信号。
 * 实质上，信号是一个对象（Triangulation类的一个成员），另一个对象可以连接到它。连接实质上是连接对象传递一个接受一定数量和种类的参数的函数对象。每当信号的拥有者想要表示某种事件时，它就会
 * "触发
 * "信号，这反过来意味着信号的所有连接都被触发：换句话说，函数对象被执行，可以采取必要的行动。
 * 作为一个简单的例子，下面的代码将在每次三角测量刚刚被完善时向输出打印一些东西。
 * @code
 *   void f()
 *   {
 *     std::cout << "Triangulation has been refined." << std::endl;
 *   }
 *
 *   void run ()
 *   {
 *     Triangulation<dim> triangulation;
 *     // fill it somehow
 *     triangulation.signals.post_refinement.connect (&f);
 *     triangulation.refine_global (2);
 *   }
 * @endcode
 * 这段代码将产生两次输出，每次精化周期一次。
 * 一个更有趣的应用是下面的，类似于FEValues类的作用。这个类存储了一个指向三角形的指针和一个指向最后处理的单元的迭代器（这样它就可以比较当前的单元和上一个单元，例如，如果新的单元是上一个单元的简单转换，那么就不需要重新计算雅各布矩阵）。然而，每当三角结构被修改时，以前处理的单元的迭代器就需要失效，因为它现在不再指向任何有用的单元（或者，至少，指向可能不一定与以前处理的单元相似的东西）。这段代码看起来是这样的（真正的代码有更多的错误检查，并且必须处理后续的单元格实际上可能属于不同的三角形的情况，但这对我们来说并不关心）。
 *
 * @code
 * template <int dim>
 * class FEValues
 * {
 * Triangulation<dim>::active_cell_iterator current_cell, previous_cell;
 * public:
 * void reinit (Triangulation<dim>::active_cell_iterator &cell);
 * void invalidate_previous_cell ();
 * };
 *
 * template <int dim>
 * void
 * FEValues<dim>::reinit (Triangulation<dim>::active_cell_iterator &cell)
 * {
 * if (previous_cell.status() != valid)
 *   {
 *     // previous_cell has not been set. set it now, and register with the
 *     // triangulation that we want to be informed about mesh refinement
 *     previous_cell = current_cell;
 *     previous_cell->get_triangulation().signals.post_refinement.connect(
 *       [this]()
 *       {
 *         this->invalidate_previous_cell();
 *       });
 *   }
 * else
 *  previous_cell = current_cell;
 *
 * current_cell = cell;
 * // ... do something with the cell...
 * }
 *
 * template <int dim>
 * void
 * FEValues<dim>::invalidate_previous_cell ()
 * {
 * previous_cell = Triangulation<dim>::active_cell_iterator();
 * }
 * @endcode
 * 这里，每当三角剖分被细化时，就会触发细化后的信号，调用附加在它身上的函数对象。这个函数对象是成员函数
 * <code>FEValues<dim>::invalidate_previous_cell</code>
 * ，我们将单一的参数（ <code>this</code>
 * 成员函数的指针，否则没有参数）绑定到FEValues对象的
 * <code>this</code>
 * 指针上。请注意，这里不需要拥有三角形和FEValues对象的代码在前者被完善时通知后者。(在实践中，该函数也希望连接到三角形所提供的其他一些信号，特别是创建和删除信号)。
 * Triangulation类有各种信号，表明三角形可以通过不同的行动来修改自己，并可能需要在其他地方采取后续行动。详情请参考
 * Triangulation::Signals 。 <h3>Serializing (loading or storing)
 * triangulations</h3>
 * 与deal.II中的许多其他类一样，Triangulation类可以使用BOOST的序列化设施将其内容流向一个档案。这样存储的数据以后可以再次从存档中检索，以恢复这个对象的内容。这个工具经常被用来将程序的状态保存到磁盘上，以便以后可能复活，通常是在长期运行的计算的检查点/重启策略的背景下，或者在不是很可靠的计算机上（例如，在非常大的集群上，个别节点偶尔会出现故障，然后导致整个MPI作业的崩溃）。
 * 由于技术原因，编写和恢复Triangulation对象并非易事。主要原因是，与许多其他对象不同，三角计算依赖于许多其他对象，它们存储指针或与之对接；例如，三角计算存储指向描述边界和流形的对象的指针，它们有存储指向其他对象的信号，以便它们能够被通知三角计算的变化（见本介绍中关于信号的部分）。由于这些对象是由用户空间拥有的（例如，用户可以创建一个自定义流形对象），它们可能无法被序列化。所以在这样的情况下，
 * boost::serialize
 * 可以存储一个对象的引用，而不是指针，但是在写的时候，这个引用永远不会被满足，因为所指向的对象没有被序列化。显然，在加载时，
 * boost::serialize
 * 将不知道让指针指向哪里，因为它从未得到重新创建最初指向的对象。
 * 由于这些原因，将三角图保存到档案中并不存储所有信息，而只是存储某些部分。更具体地说，被存储的信息是定义网格的所有信息，如顶点位置、顶点索引、顶点与单元的连接方式、边界指示器、子域ID、材料ID等。另一方面，以下信息不被存储。
 *
 *
 *
 *
 *
 * - 信号
 *
 * - 以前使用 Triangulation::set_manifold() 设置的Manifold对象的指针。
 * 另一方面，由于这些是通常在用户代码中设置的对象，它们通常可以很容易地在你重新加载三角图的那部分代码中再次设置。
 * 在某种意义上，这种序列化的方法意味着重新加载三角化更类似于调用
 * Triangulation::create_triangulation()
 * 函数，并在其中填充一些额外的内容，因为该函数也不接触属于这个三角化的信号和Manifold对象。为了保持这种类比，
 * Triangulation::load() 函数也触发了与
 * Triangulation::create_triangulation(). 相同种类的信号
 *
 *  <h3>Technical details</h3> <h4>%Algorithms for mesh regularization and
 * smoothing upon refinement</h4>
 * 我们选择了一个归纳的观点：由于在创建三角形时，所有的单元都在同一水平线上，所有关于共享一个共同面、边或顶点的单元的最大水平差异的正则性假设都是成立的。由于我们在网格历史的每一步中都使用了正则化和平滑化，所以在进一步细化时，这些假设也是成立的。
 * 正则化和平滑化是在 @p
 * 的prepare_coarsening_and_refinement函数中完成的，该函数在一开始就被
 * @p execute_coarsening_and_refinement调用。
 * 它通过查看旧网格和每个单元的细化标志来决定哪些额外的单元需要细化。
 * <ul>   <li>   <em>  正则化。 </em>  算法遍历所有单元，检查当前单元是否被标记为细化，以及当前单元的邻居是否比当前单元少细化一次。如果是，则标记该邻居进行细化。由于上面的归纳法，可能没有比现在的单元格少2级的邻居。
 * 这样标记为细化的邻居可能诱导出更多需要细化的单元。然而，这种需要额外细化的单元总是比现在的单元低一级，所以如果我们以相反的方式进行循环，从最高级别的单元开始，我们可以只对所有单元进行一次扫描。这样一来，我们就可以标记出更多低层的单元，但是如果这些单元需要更多的细化，那么在我们向后运行的循环中访问这些单元时，就会进行细化。
 * <li>   <em>  平滑化： </em>   <ul>   <li>   @p limit_level_difference_at_vertices:  首先建立一个列表，为每个顶点存储相邻单元所属的最高层。现在，由于我们在之前的细化步骤中也做了平滑处理，所以每个单元只能有顶点，其级别最多只能比当前单元的级别大一个。
 * 然而，如果我们为标记为细化的单元格存储级别加一，我们最终可能会发现单元格的顶点级别比该单元格的级别大两个。我们也需要细化这个单元，因此也需要更新其顶点的级别。这本身就可能导致需要细化的单元格，但这些单元格的层次较低，如上所述，这就是为什么我们可以只在一个循环中做各种额外标记的原因。
 * <li>   @p eliminate_unrefined_islands:
 * 对于每个单元，我们计算被细化或被标记为细化的邻居的数量。如果这个数字超过了没有被精炼和没有被标记为精炼的邻居的数量，那么当前单元格就被标记为精炼。由于这可能导致同一层次的单元也需要细化，我们将需要对所有单元进行额外的正则化和平滑化循环，直到没有任何变化。
 * <li>  <tt>eliminate_refined_*_islands</tt>:
 * 这个功能与上面的功能基本相同，但是用于粗化。如果一个单元格被标记为精简，或者它的所有子单元都是活动的，并且如果邻居的数量是活动的并且没有标记为精简，或者没有活动但所有子单元被标记为粗化的，那么这个单元的子单元被标记为粗化或者（如果这个单元被标记为精简）精简标记被清空。
 * 关于这两个版本的标志的区别，请看上面在本类描述的一般部分中关于网格平滑的部分。
 * 同样适用于上述情况：可能需要几个循环。  </ul>  </ul>  。
 * 正则化和平滑有一点互补性，当我们在一个标记为细化（正则化）的单元上或在一个未标记为细化的单元上时，我们会检查是否需要设置额外的细化标志。这使得可读的编程更容易。
 * 所有描述的算法只适用于一个以上的空间维度，因为对于一个维度没有任何限制。可能有必要对多网格算法进行一些平滑处理，但这必须在以后决定。
 *
 *  <h3>Warning</h3>
 * 似乎不可能通过迭代器的使用来保留三角形的 @p constness
 * 。因此，如果你声明指向 @p const
 * 三角形对象的指针，你应该清楚地知道你可能会不由自主地改变存储在三角形中的数据。
 *
 *
 * @ingroup grid aniso
 *
 */
template <int dim, int spacedim = dim>
class Triangulation : public Subscriptor
{
private:
  /**
   * 一个内部别名，使迭代器类的定义更简单。
   *
   */
  using IteratorSelector =
    dealii::internal::TriangulationImplementation::Iterators<dim, spacedim>;

public:
  /**
   * 声明一些网格平滑算法的符号名称。这些标志的含义在Triangulation类中有记载。
   *
   */
  enum MeshSmoothing
  {
    /**
     * 完全没有网格平滑，只是网格必须保持一个不规则。
     *
     */
    none = 0x0,
    /**
     *
     */
    limit_level_difference_at_vertices = 0x1,
    /**
     * 没有被细化的单个单元，被细化的单元所包围，通常也会导致局部的近似特性急剧下降。原因是未精化和精化单元之间的面的节点不是真实的自由度，而是带有约束。因此，没有额外自由度的补丁要比未精炼单元本身大得多。如果在传递给构造函数的参数中，#eliminate_unrefined_islands的位被设置，所有没有被标记为精炼的单元，但被比未精炼单元更多的精炼单元包围的单元都被标记为精炼。还没有被细化但被标记为细化的单元被计入细化邻居的数量。边界上的单元则完全不计算在内。根据这个定义，一个未精炼的岛也是一个被三个精炼的单元和一个未精炼的单元所包围的单元（在二维），或者一个被两个精炼的单元和一个未精炼的单元所包围的单元，并且在一侧处于边界上。因此，它并不是一个真正的岛屿，正如旗帜的名称所表明的那样。然而，到现在作者也没有想到更好的名字。
     *
     */
    eliminate_unrefined_islands = 0x2,
    /**
     * 补丁级别1的三角网格由补丁组成，也就是由被细化一次的单元组成。这个标志可以确保一个1级补丁的网格在粗化和细化之后仍然是1级补丁的。然而，在第一次调用
     * Triangulation::execute_coarsening_and_refinement()
     * 之前，用户有责任确保网格是属于补丁级的。最简单的方法是在创建三角网格后直接调用
     * global_refine(1)。
     * 由此可见，如果一个单元的至少一个子单元是或将被精炼，那么所有子单元都需要被精炼。如果设置了#patch_level_1标志，那么#eliminate_unrefined_islands、#eliminate_refined_inner_islands和#eliminate_refined_boundary_islands标志将被忽略，因为它们将被自动履行。
     *
     */
    patch_level_1 = 0x4,
    /**
     * 每个 @ref GlossCoarseMesh "粗略网格 "
     * 单元至少被精炼一次，也就是说，三角形可能在第1层有活动单元，但在第0层没有。这个标志可以确保一个具有最粗_级_1的网格在经过粗化和细化后仍然具有最粗_级_1。然而，用户有责任在第一次调用execute_coarsening_and_refinement之前，确保网格具有最粗的级别_1。最简单的方法是在创建三角网格后直接调用global_refine(1)。因此，第1层的活动单元可能不会被粗化。
     * 这个标志的主要用途是确保每个单元在每个坐标方向上至少有一个邻居（即每个单元至少有一个左或右，以及至少一个2d的上或下邻居）。这是某些计算单元间有限差异的算法的必要前提。DerivativeApproximation类是这些算法中的一种，它要求一个三角形是最粗级别的，除非所有单元在最粗级别的每个坐标方向上已经有至少一个邻居。
     *
     */
    coarsest_level_1 = 0x8,
    /**
     * 这个标志不包括在 @p maximum_smoothing.
     * 中，该标志与以下情况有关：考虑一个未精化的单元和一个精化的单元有一个共同的面，并且精化的单元沿共同面的一个子节点被标记为进一步精化的情况。在这种情况下，所得到的网格将沿着三角形的一条或多条边有一个以上的悬空节点，这种情况是不允许的。因此，为了进行细化，两个原始单元中较粗的单元也将被细化。
     * 然而，在许多情况下，以各向异性的方式细化两个原始单元中的较粗的单元就足够了，以避免在一条边上出现多个悬空顶点的情况。只做最小的各向异性细化可以节省单元和自由度。通过指定这个标志，库可以产生这些各向异性的细化。
     * 这个标志在默认情况下是不包含的，因为它可能导致各向异性细化的网格，即使没有单元被用户命令明确地进行各向异性细化。这个令人惊讶的事实可能会导致程序做错事，因为它们不是为各向异性网格可能发生的额外情况而编写的，见
     * step-30 的介绍中的讨论。
     *
     */
    allow_anisotropic_smoothing = 0x10,
    /**
     * 该算法寻求被细化或标记为细化的孤立单元。这个定义与#eliminate_unrefined_islands的定义不同，后者意味着一个岛屿被定义为一个被细化的单元，但其邻居中未被细化的单元多于被细化的单元。例如，在2D中，如果一个单元的邻居也被细化（或被细化但被标记为粗化），那么该单元的细化将被恢复。
     * 改变岛的定义的原因是，这个选项有点危险，因为如果你考虑一连串的细化单元（例如，沿着解决方案中的结点），两端的单元将被粗化，之后最外层的单元将需要被粗化。因此，只能对这样的细胞进行一次循环标记，以避免吃掉整个精炼细胞链（"链式反应"...）。
     * 这个算法也考虑到了那些实际上没有被细化但被标记为细化的单元。如果有必要的话，它会拿走细化标志。
     * 实际上这个标志有两个版本，#eliminate_refined_inner_islands和#eliminate_refined_boundary_islands。第一个消除了上面定义的位于域内部的岛屿，而第二个只消除了那些位于边界的单元的岛屿。之所以这样划分标志，是因为人们经常希望消除内部的岛屿，而边界上的岛屿却很有可能被需要，例如，当人们根据与边界积分相关的准则来细化网格，或者当人们有粗糙的边界数据时。
     *
     */
    eliminate_refined_inner_islands = 0x100,
    /**
     * 这个标志的结果与#eliminate_refined_inner_islands非常相似。见那里的文档。
     *
     */
    eliminate_refined_boundary_islands = 0x200,
    /**
     * 这个标志可以防止未精炼岛屿的出现。更详细地说。
     * 如果 "大部分邻居
     * "在该步骤后将被精炼，它将禁止粗化一个单元。
     *
     */
    do_not_produce_unrefined_islands = 0x400,

    /**
     * 这个标志总结了所有的平滑算法，这些算法可能在细化时通过标志一些更多的单元来进行细化。
     *
     */
    smoothing_on_refinement =
      (limit_level_difference_at_vertices | eliminate_unrefined_islands),
    /**
     * 这个标志总结了所有的平滑算法，这些算法可以在粗化时通过标记一些更多的单元格来进行粗化。
     *
     */
    smoothing_on_coarsening =
      (eliminate_refined_inner_islands | eliminate_refined_boundary_islands |
       do_not_produce_unrefined_islands),

    /**
     * 这个标志包括上述所有的（因此结合了所有实施的平滑算法），但各向异性的平滑算法除外。
     *
     */
    maximum_smoothing = 0xffff ^ allow_anisotropic_smoothing
  };

  /**
   * 一个别名，用于识别单元格迭代器。迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。
   * 当前的别名用于识别三角形中的单元。TriaIterator类的工作原理就像一个指针，当你解除引用时，会产生一个CellAccessor类型的对象。CellAccessor是一个标识三角形中单元格特定属性的类，但它派生（因此继承）于TriaAccessor，TriaAccessor描述了你可以对三角形中更一般的对象（线、面以及单元格）提出什么要求。
   * @ingroup Iterators
   *
   */
  using cell_iterator = TriaIterator<CellAccessor<dim, spacedim>>;

  /**
   * 和上面一样，允许在细化层次上使用 "MeshType概念"。
   *
   */
  using level_cell_iterator = cell_iterator;

  /**
   * @ingroup Iterators
   *
   */
  using active_cell_iterator = TriaActiveIterator<CellAccessor<dim, spacedim>>;

  /**
   * 一个别名，用于识别指向面的迭代器。  迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。
   * 当前的别名是识别三角形中的面。TriaIterator类的工作原理就像一个指向对象的指针，当你解除引用时，会产生一个TriaAccessor类型的对象，即可以用来查询面的几何属性的类，如它们的顶点、面积等。
   * @ingroup Iterators
   *
   */
  using face_iterator = TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>;

  /**
   * 一个别名，用于识别指向活动面的迭代器，即指向没有子节点的面。活动面必须是至少一个活动单元的面。
   * 除了 "活动 "的限定，这个别名与 @p face_iterator
   * 的别名相同。特别是，取消引用这两个别名都会产生相同的对象。
   * @ingroup Iterators
   *
   */
  using active_face_iterator =
    TriaActiveIterator<TriaAccessor<dim - 1, dim, spacedim>>;

  /**
   * 一个定义了迭代器类型的别名，用于迭代网格的顶点。 迭代器的概念在 @ref Iterators "迭代器文档模块 "
   * 中有详细的讨论。
   * @ingroup Iterators
   *
   */
  using vertex_iterator = TriaIterator<dealii::TriaAccessor<0, dim, spacedim>>;

  /**
   * 一个别名，定义了一个迭代器类型，用于迭代网格的顶点。 迭代器的概念在 @ref Iterators  "迭代器文档模块 "
   * 中有详细的讨论。    这个别名实际上与上面的 @p
   * vertex_iterator
   * 别名相同，因为网格中的所有顶点都是活动的（即，是活动单元的一个顶点）。
   * @ingroup Iterators
   *
   */
  using active_vertex_iterator =
    TriaActiveIterator<dealii::TriaAccessor<0, dim, spacedim>>;

  /**
   * 一个别名，定义了一个网格的（一维）线条的迭代器。在一维网格中，这些线是网格的单元，而在二维网格中，线是单元的面。
   * @ingroup Iterators
   *
   */
  using line_iterator = typename IteratorSelector::line_iterator;

  /**
   * 一个别名，允许在<i>active</i>线上迭代，即没有子节点的线的子集。在一维网格中，这些是网格的单元，而在二维网格中，线是单元的面。
   * 在二维或三维网格中，没有子节点的线（即活动线）是至少一个活动单元的一部分。每条这样的线还可能是与活动单元相邻的更粗的单元的线的子线。(这个较粗的邻居也是活动的)。
   * @ingroup Iterators
   *
   */
  using active_line_iterator = typename IteratorSelector::active_line_iterator;

  /**
   * 一个别名，定义了一个网格的（二维）四边形的迭代器。在二维网格中，这些是网格的单元，而在三维网格中，四边形是单元的面。
   * @ingroup Iterators
   *
   */
  using quad_iterator = typename IteratorSelector::quad_iterator;

  /**
   * 一个别名，允许在<i>active</i>四边形上进行迭代，即没有子节点的四边形子集。在二维网格中，这些是网格的单元，而在三维网格中，四边形是单元的面。
   * 在三维网格中，没有孩子的四边形（即活动四边形）是至少一个活动单元的面。每个这样的四边形还可能是与活动单元相邻的更粗的单元的四边形面的子。这个较粗的邻居也将是活动的）。
   * @ingroup Iterators
   *
   */
  using active_quad_iterator = typename IteratorSelector::active_quad_iterator;

  /**
   * 一个别名，定义了一个网格的（三维）六边形的迭代器。这个迭代器只有在三维网格中才有意义，在三维网格中六边形是网格的单元。
   * @ingroup Iterators
   *
   */
  using hex_iterator = typename IteratorSelector::hex_iterator;

  /**
   * 一个别名，允许在网格的<i>active</i>六边形上迭代。
   * 这个迭代器只有在三维网格中才有意义，在三维网格中六边形是网格的单元。因此，在这些三维网格中，这个迭代器等同于
   * @p active_cell_iterator  别名。
   * @ingroup Iterators
   *
   */
  using active_hex_iterator = typename IteratorSelector::active_hex_iterator;

  /**
   * 一个结构，被create_triangulation()函数用作异常对象，用于指示粗略网格单元中哪些单元是倒置的或严重扭曲的（见术语表中 @ref GlossDistorted "扭曲的单元 "
   * 条目）。
   * 这类对象会被create_triangulation()和execute_coarsening_and_refinement()函数抛出，如果要忽略这个条件，可以在用户代码中捕获它们。然而，请注意，只有在调用三角形类的构造函数时表明有必要进行这种检查，才会产生这种异常。
   * 如果从参考单元到实数单元的映射的Jacobian的行列式至少在一个顶点是负的，则该单元被称为<i>deformed</i>。这一计算是通过
   * GeometryInfo::jacobian_determinants_at_vertices 函数完成的。
   *
   */
  struct DistortedCellList : public dealii::ExceptionBase
  {
    /**
     * 解构器。空的，但为了异常规范而需要，因为基类有这个异常规范，而自动生成的析构器会因为成员对象而有一个不同的规范。
     *
     */
    virtual ~DistortedCellList() noexcept override;

    /**
     * 粗略网格单元中那些被变形或其子单元被变形的单元列表。
     *
     */
    std::list<typename Triangulation<dim, spacedim>::cell_iterator>
      distorted_cells;
  };

  /**
   * 使尺寸在函数模板中可用。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 使空间维度在函数模板中可用。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 创建一个空的三角结构。不要创建任何单元。      @param  smooth_grid 决定在网格细化时应强制执行的网格尺寸函数的平滑程度。      @param  check_for_distorted_cells 决定三角网格是否应该检查由create_triangulation()或execute_coarsening_and_refinement()创建的任何单元是否扭曲（见 @ref GlossDistorted  "扭曲的单元"
   * ）。
   * 如果设置，这两个函数在遇到扭曲的单元格时可能会抛出一个异常。
   *
   */
  Triangulation(const MeshSmoothing smooth_grid               = none,
                const bool          check_for_distorted_cells = false);

  /**
   * 复制构造函数。    你真的应该使用 @p copy_triangulation
   * 函数，所以这个构造函数被删除。这样做的原因是，我们可能想在集合中使用三角函数对象。然而，C++容器要求存储在其中的对象是可复制的，所以我们需要提供一个复制构造函数。另一方面，复制三角形是非常昂贵的，我们不希望这种对象被意外地复制，例如在编译器生成的临时对象中。通过定义一个复制构造函数但抛出一个错误，我们满足了容器的形式要求，但同时又不允许实际的复制。
   * 最后，通过这个异常，我们很容易找到需要修改代码以避免拷贝的地方。
   *
   */
  Triangulation(const Triangulation<dim, spacedim> &) = delete;

  /**
   * 移动构造器。
   * 通过窃取另一个三角结构的内部数据来创建一个新的三角结构。
   *
   */
  Triangulation(Triangulation<dim, spacedim> &&tria) noexcept;

  /**
   * 移动赋值运算符。
   *
   */
  Triangulation &
  operator=(Triangulation<dim, spacedim> &&tria) noexcept;

  /**
   * 删除对象和所有级别的层次结构。
   *
   */
  virtual ~Triangulation() override;

  /**
   * 通过删除所有数据，将这个三角结构重置为处女状态。
   * 注意，只有当这个对象不再存在任何订阅时，才允许进行这个操作，比如使用它的DoFHandler对象。
   *
   */
  virtual void
  clear();

  /**
   * 返回该三角函数所使用的MPI通信器。如果是一个串行Triangulation对象，将返回MPI_COMM_SELF。
   *
   */
  virtual MPI_Comm
  get_communicator() const;

  /**
   * 设置网格平滑度为 @p mesh_smoothing.
   * 这将覆盖给构造函数的MeshSmoothing。只有在三角结构为空时才允许调用这个函数。
   *
   */
  virtual void
  set_mesh_smoothing(const MeshSmoothing mesh_smoothing);

  /**
   * 返回遵守的网格平滑要求。
   *
   */
  virtual const MeshSmoothing &
  get_mesh_smoothing() const;

  /**
   * 将一个流形对象分配给三角剖面的某一部分。如果一个流形编号为
   * @p number
   * 的对象被细化，这个对象会被用来寻找新顶点的位置（参见
   * step-49
   * 的结果部分，对这个问题有更深入的讨论，并有例子）。
   * 它也被用于形状函数计算中的单元格向单元格的非线性（即：非Q1）转换。
   * 使用Manifold<dim,  spacedim>::clone() 创建一个 @p manifold_object
   * 的副本，并在内部存储。
   * 在非空三角的有效期内，有可能移除或替换一个Manifold对象。通常情况下，这是在第一次细化之前进行的，之后就很危险了。移除流形对象是通过reset_manifold()完成的。然后，该操作将之前给出的流形对象替换成一个直的流形近似值。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_manifold(const types::manifold_id       number,
               const Manifold<dim, spacedim> &manifold_object);

  /**
   * 重置三角形中具有给定 @p manifold_number
   * 的那些部分，以使用FlatManifold对象。这是一个非弯曲三角形的默认状态，并撤销由函数
   * Triangulation::set_manifold(). 分配的不同的Manifold对象。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  reset_manifold(const types::manifold_id manifold_number);

  /**
   * 重置三角形的所有部分，无论其manifold_id如何，都要使用FlatManifold对象。这将撤销函数对所有Manifold对象的分配
   * Triangulation::set_manifold().  。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  reset_all_manifolds();

  /**
   * 将所有单元格和面的manifold_id设置为给定参数。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_all_manifold_ids(const types::manifold_id number);

  /**
   * 将所有边界面的manifold_id设置为给定参数。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_all_manifold_ids_on_boundary(const types::manifold_id number);

  /**
   * 将所有边界面和边的manifold_id与给定的 boundary_id  @p b_id
   * 设置为给定的manifold_id  @p number.  。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_all_manifold_ids_on_boundary(const types::boundary_id b_id,
                                   const types::manifold_id number);

  /**
   * 返回一个对用于该三角测量的流形对象的常数引用。  @p
   * number 与set_manifold()中相同。
   * @note  如果找不到流形，则返回默认的平面流形。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  const Manifold<dim, spacedim> &
  get_manifold(const types::manifold_id number) const;

  /**
   * 返回一个向量，包含分配给该三角测量对象的活动单元的边界面的所有边界指标。注意，每个边界指标只被报告一次。返回向量的大小将代表不同指标的数量（大于或等于1）。
   * @ingroup boundary   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  virtual std::vector<types::boundary_id>
  get_boundary_ids() const;

  /**
   * 返回一个向量，该向量包含分配给该三角结构中活动单元格对象的所有流形指标。注意，每个流形指标只报告一次。返回向量的大小将代表不同指标的数量（大于或等于1）。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  virtual std::vector<types::manifold_id>
  get_manifold_ids() const;

  /**
   * 将 @p other_tria
   * 复制到这个三角区。这个操作并不便宜，所以你应该小心使用这个。我们没有把这个函数作为一个复制构造函数来实现，因为如果你可以在以后给它们赋值，那么维护三角形的集合就更容易了。
   * 请记住，这个函数也复制了之前由 @p set_manifold
   * 函数设置的边界描述符的指针。因此，你还必须保证描述边界的Manifold对象的寿命至少与复制的三角图一样长。
   * 这个三角形必须事先是空的。    该函数被制成 @p virtual
   * ，因为一些派生类可能想禁用或扩展该函数的功能。
   * @note
   * 调用这个函数会在other_tria上触发'copy'信号，也就是被复制的三角结构<i>from</i>。
   * 它还会触发当前三角形的'创建'信号。更多信息请参见通用文档中关于信号的部分。
   * @note
   * 信号的连接列表不会从旧的三角结构复制到新的三角结构，因为这些连接的建立是为了监视旧的三角结构如何变化，而不是它可能被复制到的任何三角结构如何变化。
   *
   */
  virtual void
  copy_triangulation(const Triangulation<dim, spacedim> &other_tria);

  /**
   *
   *
   *
   *
   * - 例如，它可能是行列式为零(表明你在一个单元中的边缘塌陷了)，但这是可以的，因为你并不打算在这个单元上进行积分。另一方面，变形的单元通常表明网格太粗，无法解决领域的几何问题，在这种情况下，忽略这个例外可能是不明智的。
   * @note  这个函数在  step-14  和  step-19  中使用。
   * @note  这个函数在做完它的工作后会触发 "创建
   * "信号。参见该类的一般文档中关于信号的部分。例如，作为其结果，所有连接到这个三角形的DoFHandler对象将通过
   * DoFHandler::reinit().  被重新初始化。
   * @note
   * 只有在dim==spacedim的情况下才会对扭曲的单元进行检查，否则，如果单元所描述的流形是扭曲的，那么单元就可以合法地被扭曲。
   *
   */
  virtual void
  create_triangulation(const std::vector<Point<spacedim>> &vertices,
                       const std::vector<CellData<dim>> &  cells,
                       const SubCellData &                 subcelldata);

  /**
   * 从提供的 TriangulationDescription::Description.
   * 中创建一个三角形。
   * @note
   * 如果需要流形，别忘了在调用此函数前用set_manifold()附加流形。
   * @note  命名空间 TriangulationDescription::Utilities 包含创建
   * TriangulationDescription::Description.   @param  construction_data
   * 这个过程需要的数据。
   *
   */
  virtual void
  create_triangulation(
    const TriangulationDescription::Description<dim, spacedim>
      &construction_data);

  /**
   * 仅用于向后兼容。这个函数按照5.2之前的deal.II版本的要求，以排序方式获取单元格数据，将其转换为新的（lexicographic）排序，并调用create_triangulation()。
   * @note
   * 该函数内部调用create_triangulation，因此可以抛出与其他函数相同的异常。
   *
   */
  DEAL_II_DEPRECATED
  virtual void
  create_triangulation_compatibility(
    const std::vector<Point<spacedim>> &vertices,
    const std::vector<CellData<dim>> &  cells,
    const SubCellData &                 subcelldata);

  /**
   * 还原或翻转dim<spacedim三角形的方向标志，见
   * @ref GlossDirectionFlag  。
   * 如果dim等于spacedim，这个函数会抛出一个异常。
   *
   */
  void
  flip_all_direction_flags();

  /**
   * @name  网格细化  
     * @{ 
   *
   */

  /**
   * 标记所有活动单元进行细化。
   * 这将细化所有级别的、尚未细化的单元格（也就是说，只有尚未有子代的单元格才会被细化）。这些单元格只被标记，不被细化，因此你有机会保存细化标记。
   *
   */
  void
  set_all_refine_flags();

  /**
   * 精炼所有单元格 @p times 次。换句话说，在每一次的 @p
   * times
   * 迭代中，循环所有单元格，并将每个单元格统一细化为
   * $2^\text{dim}$ 个子。在实践中，这个函数重复了以下操作
   * @p times
   * 次：调用set_all_refine_flags()，然后是execute_coarsening_and_refinement()。最终的结果是，单元格的数量增加了
   * $(2^\text{dim})^\text{times}=2^{\text{dim} \times \text{times}}$  倍。
   * 在这个循环中调用的execute_coarsening_and_refinement()函数如果创建了扭曲的单元格，可能会抛出一个异常（见其文档解释）。如果发生这种情况，这个异常将通过这个函数传播，在这种情况下，你可能不会得到实际的细化步骤数。
   * @note
   * 这个函数在做每个单独的细化循环之前和之后都会触发细化前和细化后的信号（即如果`times
   * >
   * 1'的话，会不止一次）。参见该类的一般文档中关于信号的部分。
   *
   */
  void
  refine_global(const unsigned int times = 1);

  /**
   * 对所有单元格进行给定次数的粗化。    在 @p times
   * 的每一次迭代中，所有单元格都将被标记为粗化。如果一个活动的单元格已经在最粗的层次上，它将被忽略。
   * @note
   * 该函数在做每个单独的粗化循环（即如果`时间>1'，则不止一次）之前和之后都会触发预精化信号。参见该类的一般文档中关于信号的部分。
   *
   */
  void
  coarsen_global(const unsigned int times = 1);

  /**
   * 同时执行三角形的细化和粗化。
   * 该函数将所有细化和粗化的标志重置为假。它将用户标志用于内部目的。因此，它们将被未定义的内容所覆盖。
   * 为了允许用户程序修复这些单元，如果需要的话，这个函数在完成所有其他工作后可能会抛出一个类型为DistortedCellList的异常，该异常包含一个已经被细化并且至少有一个子单元被扭曲的单元列表。如果没有单元格创建了扭曲的子代，该函数就不会创建这样的异常。
   * 请注意，为了实现对变形单元的检查，必须在创建三角形对象时指定
   * <code>check_for_distorted_cells</code> 标志。
   * 更多信息请参见通用文档。
   * @note
   * 这个函数在工作之前和之后都会触发精简前和精简后的信号。参见该类的一般文档中关于信号的部分。
   * @note  如果边界描述足够不规则，可能会发生一些由网格细化产生的子节点被扭曲（见 @ref GlossDistorted "扭曲的单元 "
   * 的广泛讨论）。
   * @note
   * 这个函数是<tt>虚拟的</tt>，以允许派生类插入钩子，如保存细化标志等（见例如PersistentTriangulation类）。
   *
   */
  virtual void
  execute_coarsening_and_refinement();

  /**
   * 既做细化和粗化的准备，也做网格的平滑。
   * 关于细化过程，它在<tt>dim>=2</tt>中固定细化的闭合（确保没有两个单元在细化水平上相差超过一个），等等。
   * 如果这个类的构造函数给出了相应的标志，它将执行一些网格平滑处理。
   * 该函数返回是否有额外的单元被标记为细化。
   * 更多关于细化后平滑的信息，请参见本类的通用文档。
   * 关于粗化部分，在准备实际的粗化步骤时，要对单元进行标记和去标记。这包括从可能不被删除的单元中删除粗化标志（例如，因为一个邻居比该单元更精细），做一些平滑处理，等等。
   * 其效果是，只有那些被标记为粗化的单元才会被实际粗化。这包括所有被标记的单元格都属于父单元格，其所有子单元格都被标记。
   * 该函数返回一些单元格的标记是否在这个过程中被改变。
   * 这个函数使用了用户标志，所以如果你事后还需要它们，请储存它们。
   *
   */
  virtual bool
  prepare_coarsening_and_refinement();

  /*  @} 。  
*
*/

  /**
   * @name  跟上一个三角形所发生的事情  
     * @{ 
   *
   */


  /**
   * 用于通知派生类中的函数，给定cell_iterator的单元格将如何变化。请注意，这可能与平行计算中cell_iterator中的refine_flag()和coarsen_flag()不同，因为细化约束是本机看不到的。
   *
   */
  enum CellStatus
  {
    /**
     * 单元不会被细化或粗化，可能会也可能不会移动到不同的处理器。
     *
     */
    CELL_PERSIST,
    /**
     * 该单元将被或被精简。
     *
     */
    CELL_REFINE,
    /**
     * 该单元的子单元将被或被粗化为该单元。
     *
     */
    CELL_COARSEN,
    /**
     * 无效状态。不会发生在用户身上。
     *
     */
    CELL_INVALID
  };

  /**
   * 一个用于累积下面的cell_weights槽函数的结果的结构。它接受一个迭代器范围，并返回值的总和。
   *
   */
  template <typename T>
  struct CellWeightSum
  {
    using result_type = T;

    template <typename InputIterator>
    T
    operator()(InputIterator first, InputIterator last) const
    {
      return std::accumulate(first, last, T());
    }
  };

  /**
   * 一个拥有 boost::signal
   * 对象的结构，用于三角形可以对自己做的一些操作。更多的信息和例子，请参考Triangulation类的一般文档中的
   * "三角形变化时获得通知 "部分。
   * 关于信号的文档，见http://www.boost.org/doc/libs/release/libs/signals2
   * 。
   *
   */
  struct Signals
  {
    /**
     * 只要调用 Triangulation::create_triangulation 或
     * Triangulation::copy_triangulation()
     * 就会触发这个信号。当通过 Triangulation::load().
     * 从档案中加载三角图时也会触发这个信号。
     *
     */
    boost::signals2::signal<void()> create;

    /**
     * 这个信号在 Triangulation::execute_coarsening_and_refinement()
     * 函数（其本身被其他函数如 Triangulation::refine_global()
     * 调用）开始执行时被触发。在这个信号被触发的时候，三角结构仍然没有改变。
     *
     */
    boost::signals2::signal<void()> pre_refinement;

    /**
     * 这个信号在 Triangulation::execute_coarsening_and_refinement()
     * 函数的执行结束时被触发，此时三角结构已经达到最终状态。
     *
     */
    boost::signals2::signal<void()> post_refinement;

    /**
     * 该信号在 GridTools::partition_triangulation() 和
     * GridTools::partition_triangulation_zorder()
     * 函数开始执行时被触发。在该信号被触发时，三角结构仍未改变。
     *
     */
    boost::signals2::signal<void()> pre_partition;

    /**
     * 当deal.II中的函数移动网格点时，这个信号会被触发，例如
     * GridTools::transform.  不幸的是，通过
     * <code>cell_iterator->vertex(v) = xxxx</code>
     * 修改用户代码中的顶点不能被这种方法检测到。
     *
     */
    boost::signals2::signal<void()> mesh_movement;

    /**
     * 这个信号对每个要被粗化的单元都会被触发。
     * @note
     * 该信号是以一组活动单元的直接父单元为参数触发的。该父单元的子单元随后将被粗化。
     *
     */
    boost::signals2::signal<void(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell)>
      pre_coarsening_on_cell;

    /**
     * 这个信号对每一个刚刚被粗化的单元格都会被触发。
     * @note  信号参数 @p cell
     * 对应于一组新创建的活动单元的直接父单元。
     *
     */
    boost::signals2::signal<void(
      const typename Triangulation<dim, spacedim>::cell_iterator &cell)>
      post_refinement_on_cell;

    /**
     * 每当拥有该信号的三角结构被另一个使用
     * Triangulation::copy_triangulation()
     * 的三角结构复制时，该信号就会被触发（即在<i>old</i>三角结构上被触发，但新的三角结构被作为参数传递）。
     *
     */
    boost::signals2::signal<void(
      const Triangulation<dim, spacedim> &destination_tria)>
      copy;

    /**
     * 只要调用 Triangulation::clear()
     * 函数，并在三角形的析构器中，这个信号就会被触发。当通过
     * Triangulation::load()
     * 从存档中加载三角函数时，该信号也会被触发，因为三角函数的先前内容首先被销毁。
     * 该信号在三角剖析的数据结构被销毁之前被触发。换句话说，连接到这个信号的函数可以最后看一下三角图，例如保存作为三角图的一部分存储的信息。
     *
     */
    boost::signals2::signal<void()> clear;

    /**
     * 这是一个总括性的信号，每当创建、细化后或清除信号被触发时都会被触发。实际上，它可以用来向连接到该信号的对象表明三角测量已经被改变，不管改变的具体原因是什么。
     * @note 单元级信号 @p pre_coarsening_on_cell 和 @p
     * post_refinement_on_cell不与此信号连接。
     *
     */
    boost::signals2::signal<void()> any_change;

    /**
     * 这个信号在每次自动或手动重新分区时对每个单元都会被触发。这个信号有些特殊，因为它只在分布式并行计算中被触发，而且只有当函数连接到它时才会被触发。它的目的是允许对领域进行加权重新划分，以平衡各进程的计算负载，而不是平衡单元的数量。任何连接的函数都需要一个指向单元格的迭代器，以及一个CellStatus参数，该参数表示该单元格是要被细化、粗化还是不被触及（更多信息请参见CellStatus枚举的文档）。该函数预计将返回一个无符号整数，它被解释为该单元的额外计算负荷。如果这个单元要被粗化，信号是为父单元调用的，你需要提供未来父单元的权重。如果这个单元将被精简，该函数应返回一个权重，该权重将被平均分配给当前单元的每个未来子单元。作为参考，每个单元格的总权重都要加上1000的值。这意味着信号返回值为1000（导致权重为2000）意味着一个进程处理这个特定单元的成本是两倍。如果有几个函数连接到这个信号，它们的返回值将被相加，以计算出最终的权重。
     * 这个函数在 step-68 中使用。
     *
     */
    boost::signals2::signal<unsigned int(const cell_iterator &,
                                         const CellStatus),
                            CellWeightSum<unsigned int>>
      cell_weight;

    /**
     * 这个信号在
     * parallel::distributed::Triangulation::execute_coarsening_and_refinement()
     * 函数开始执行时被触发（该函数本身被其他函数如
     * Triangulation::refine_global()
     * 调用）。在这个信号被触发的时候，三角结构仍然没有改变。这个信号与pre_refinement信号不同，因为在并行分布的情况下，pre_refinement信号被多次触发，没有办法区分最后的信号调用。
     *
     */
    boost::signals2::signal<void()> pre_distributed_refinement;

    /**
     * 这个信号是在执行
     * parallel::distributed::Triangulation::execute_coarsening_and_refinement()
     * 函数时触发的。在这个信号被触发的时候，p4est神谕已经被完善，单元格关系已经被更新。否则，三角结构是不变的，而且p4est谕言还没有被重新划分。
     *
     */
    boost::signals2::signal<void()> post_p4est_refinement;

    /**
     * 这个信号在
     * parallel::distributed::Triangulation::execute_coarsening_and_refinement()
     * 函数的执行结束时被触发，此时三角剖分已经达到最终状态。这个信号与post_refinement信号不同，因为在并行分布的情况下，post_refinement信号会被多次触发，而没有办法区分最后的信号调用。
     *
     */
    boost::signals2::signal<void()> post_distributed_refinement;

    /**
     * 这个信号在 parallel::distributed::Triangulation::repartition()
     * 函数的执行开始时被触发。在这个信号被触发的时候，三角关系仍然没有改变。
     * @note  parallel::distributed::Triangulation::repartition() 函数也被
     * parallel::distributed::Triangulation::load().
     * 调用。因此，pre_distributed_repartition信号将在pre_distributed_load信号之后被触发。
     *
     */
    boost::signals2::signal<void()> pre_distributed_repartition;

    /**
     * 这个信号在 parallel::distributed::Triangulation::repartition()
     * 函数的执行结束时被触发，此时三角化已经达到最终状态。
     *
     */
    boost::signals2::signal<void()> post_distributed_repartition;

    /**
     * 该信号在 parallel::distributed::Triangulation::save()
     * 函数开始执行时被触发。在这个信号被触发的时候，三角测量仍然没有变化。
     *
     */
    boost::signals2::signal<void()> pre_distributed_save;

    /**
     * 在 parallel::distributed::Triangulation::save()
     * 函数执行结束时，该信号被触发，此时三角结构已达到最终状态。
     *
     */
    boost::signals2::signal<void()> post_distributed_save;

    /**
     * 该信号在 parallel::distributed::Triangulation::load()
     * 函数开始执行时被触发。在这个信号被触发的时候，三角测量仍然没有变化。
     *
     */
    boost::signals2::signal<void()> pre_distributed_load;

    /**
     * 在 parallel::distributed::Triangulation::load()
     * 函数执行结束时，该信号被触发，此时三角结构已达到最终状态。
     *
     */
    boost::signals2::signal<void()> post_distributed_load;
  };

  /**
   * 三角形可以对自己进行的各种操作的信号。
   *
   */
  mutable Signals signals;

  /*  @} .   
*
*/

  /**
   * @name  一个三角形的历史  
     * @{ 
   *
   */

  /**
   * 将被标记为细化的单元格的地址保存到 @p 外。
   * 关于用法，请阅读该类的一般文档。
   *
   */
  void
  save_refine_flags(std::ostream &out) const;

  /**
   * 与上述相同，但将标志存储到一个位向量而不是文件中。
   *
   */
  void
  save_refine_flags(std::vector<bool> &v) const;

  /**
   * 读取由 @p save_refine_flags. 存储的信息。
   *
   */
  void
  load_refine_flags(std::istream &in);

  /**
   * 读取由 @p save_refine_flags. 存储的信息。
   *
   */
  void
  load_refine_flags(const std::vector<bool> &v);

  /**
   * 类似于 @p save_refine_flags. 的内容。
   *
   */
  void
  save_coarsen_flags(std::ostream &out) const;

  /**
   * 与上述相同，但将标志存储到一个位向量，而不是存储到一个文件。
   *
   */
  void
  save_coarsen_flags(std::vector<bool> &v) const;

  /**
   * 类似于 @p load_refine_flags. 。
   *
   */
  void
  load_coarsen_flags(std::istream &out);

  /**
   * 类似于 @p load_refine_flags. 。
   *
   */
  void
  load_coarsen_flags(const std::vector<bool> &v);

  /**
   * 返回该三角形是否经历过各向异性（而不是各向同性）的细化。
   *
   */
  bool
  get_anisotropic_refinement_flag() const;

  /*  @} .   
*
*/

  /**
   * @name  用户数据  @{
   *
   */

  /**
   * 清除所有的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  clear_user_flags();

  /**
   * 保存所有的用户标志。更多细节请参见该类的一般文档和
   * @p save_refine_flags  的文档。 也请参见  @ref GlossUserFlags  。
   *
   */
  void
  save_user_flags(std::ostream &out) const;

  /**
   * 和上面一样，但是把标志存储到一个位向量而不是文件中。
   * 如果有必要，输出向量会被调整大小。 也见  @ref
   * GlossUserFlags  。
   *
   */
  void
  save_user_flags(std::vector<bool> &v) const;

  /**
   * 读取由
   * @p save_user_flags.  储存的信息 也见  @ref GlossUserFlags  .
   *
   */
  void
  load_user_flags(std::istream &in);

  /**
   * 读取由 @p
   * save_user_flags. 存储的信息 也见 @ref GlossUserFlags  .
   *
   */
  void
  load_user_flags(const std::vector<bool> &v);

  /**
   * 清除所有行上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  clear_user_flags_line();

  /**
   * 保存行上的用户标志。 也见
   * @ref GlossUserFlags  .
   *
   */
  void
  save_user_flags_line(std::ostream &out) const;

  /**
   * 与上述相同，但将标志存储到一个位向量而不是文件中。
   * 如果有必要，输出向量会被调整大小。 也见  @ref
   * GlossUserFlags  。
   *
   */
  void
  save_user_flags_line(std::vector<bool> &v) const;

  /**
   * 加载位于行上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_line(std::istream &in);

  /**
   * 加载位于行上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_line(const std::vector<bool> &v);

  /**
   * 清除四边形上的所有用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  clear_user_flags_quad();

  /**
   * 保存四边形上的用户标志。 参见
   * @ref GlossUserFlags  。
   *
   */
  void
  save_user_flags_quad(std::ostream &out) const;

  /**
   * 与上述相同，但将标志存储到一个位向量而不是文件中。
   * 如果有必要，输出向量会被调整大小。 也见  @ref
   * GlossUserFlags  。
   *
   */
  void
  save_user_flags_quad(std::vector<bool> &v) const;

  /**
   * 加载位于四边形上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_quad(std::istream &in);

  /**
   * 加载位于四边形上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_quad(const std::vector<bool> &v);


  /**
   * 清除四边形上的所有用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  clear_user_flags_hex();

  /**
   * 保存赫兹上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  save_user_flags_hex(std::ostream &out) const;

  /**
   * 与上述相同，但将标志存储到一个位向量而不是文件中。
   * 如果有必要，输出向量会被调整大小。 也见  @ref
   * GlossUserFlags  。
   *
   */
  void
  save_user_flags_hex(std::vector<bool> &v) const;

  /**
   * 加载位于hexs上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_hex(std::istream &in);

  /**
   * 加载位于hexs上的用户标志。 也见
   * @ref GlossUserFlags  。
   *
   */
  void
  load_user_flags_hex(const std::vector<bool> &v);

  /**
   * 清除所有的用户指针和索引，并允许在下一次访问时使用这两者。 也见
   * @ref GlossUserData  。
   *
   */
  void
  clear_user_data();

  /**
   * 保存所有的用户索引。如果有必要，输出向量的大小会被调整。也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_indices(std::vector<unsigned int> &v) const;

  /**
   * 读取由save_user_indices()存储的信息。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_indices(const std::vector<unsigned int> &v);

  /**
   * 保存所有的用户指针。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_pointers(std::vector<void *> &v) const;

  /**
   * 读取由save_user_pointers()存储的信息。 也见
   * @ref GlossUserData  .
   *
   */
  void
  load_user_pointers(const std::vector<void *> &v);

  /**
   * 保存用户索引的行数。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_indices_line(std::vector<unsigned int> &v) const;

  /**
   * 加载位于行上的用户索引。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_indices_line(const std::vector<unsigned int> &v);

  /**
   * 保存位于四边形的用户索引。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_indices_quad(std::vector<unsigned int> &v) const;

  /**
   * 加载位于四边形上的用户索引。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_indices_quad(const std::vector<unsigned int> &v);

  /**
   * 保存位于十六进制的用户索引。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_indices_hex(std::vector<unsigned int> &v) const;

  /**
   * 加载位于十六进制上的用户索引。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_indices_hex(const std::vector<unsigned int> &v);
  /**
   * 保存位于行上的用户索引。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_pointers_line(std::vector<void *> &v) const;

  /**
   * 加载位于行上的用户指针。 也见
   * @ref GlossUserData  .
   *
   */
  void
  load_user_pointers_line(const std::vector<void *> &v);

  /**
   * 保存位于四边形上的用户指针。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_pointers_quad(std::vector<void *> &v) const;

  /**
   * 加载位于四边形上的用户指针。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_pointers_quad(const std::vector<void *> &v);

  /**
   * 保存位于十六进制的用户指针。如果有必要，输出向量会被调整大小。 也见
   * @ref GlossUserData  。
   *
   */
  void
  save_user_pointers_hex(std::vector<void *> &v) const;

  /**
   * 加载位于六角星上的用户指针。 也见
   * @ref GlossUserData  。
   *
   */
  void
  load_user_pointers_hex(const std::vector<void *> &v);

  /*  @}  。  
*
*/

  /**
   * @name  细胞迭代器函数  
     * @{ 
   *
   */

  /**
   * 迭代到第一层使用的单元格  @p level.  。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  cell_iterator
  begin(const unsigned int level = 0) const;

  /**
   * 迭代到第一层的活动单元  @p level.
   * 如果给定的层不包含任何活动单元（即这一层的所有单元都被进一步细化了，那么这个函数返回
   * <code>end_active(level)</code> ，这样的循环就会有
   * @code
   * for (const auto cell=tria.begin_active(level);
   *      cell!=tria.end_active(level);
   *      ++cell)
   *   {
   *     ...
   *   }
   * @endcode
   * 迭代次数为零，如果这一层没有活动单元，可能会有这样的预期。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  active_cell_iterator
  begin_active(const unsigned int level = 0) const;

  /**
   * 过了终点的迭代器；这个迭代器用于比较具有过了终点或开始前状态的迭代器。
   *
   */
  cell_iterator
  end() const;

  /**
   * 返回一个迭代器，它是第一个不在水平线上的迭代器。如果
   * @p level 是最后一个层次，那么这将返回<tt>end()</tt>。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  cell_iterator
  end(const unsigned int level) const;

  /**
   * 返回一个活动的迭代器，它是不在给定层次上的第一个活动迭代器。如果
   * @p level 是最后一个层次，那么这将返回<tt>end()</tt>。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  active_cell_iterator
  end_active(const unsigned int level) const;


  /**
   * 返回一个指向最后使用的单元格的迭代器。
   *
   */
  cell_iterator
  last() const;

  /**
   * 返回一个指向最后一个活动单元格的迭代器。
   *
   */
  active_cell_iterator
  last_active() const;

  /**
   * 返回一个迭代器，指向这个由独立的CellId对象构建的Triangulation对象的一个单元。
   * 如果给定的参数对应于这个三角形中的一个有效单元，对于当前处理器存储所有属于三角形的单元的顺序三角形，这个操作总是会成功。另一方面，如果这是一个平行三角形，那么当前的处理器可能实际上不知道这个单元。在这种情况下，对于本地相关的单元，这个操作会成功，但对于在当前处理器上不太精细的人工单元，这个操作可能不会成功。
   *
   */
  cell_iterator
  create_cell_iterator(const CellId &cell_id) const;

  /**
   * @name  细胞迭代器函数返回迭代器的范围
   *
   */

  /**
   * 返回一个迭代器范围，该范围包含构成这个三角形的所有单元格（无论是否激活）。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @return
   * 半开范围  <code>[this->begin(), this->end())</code>
   * @ingroup CPP11
   *
   */
  IteratorRange<cell_iterator>
  cell_iterators() const;

  /**
   * #Range-based_for_loop">here</a>）。一个例子是，如果没有基于范围的for循环，人们往往会写出如下的代码（假设我们的目标是在每个活动单元上设置用户标志）。
   * @code
   * Triangulation<dim> triangulation;
   * ...
   * typename Triangulation<dim>::active_cell_iterator
   *   cell = triangulation.begin_active(),
   *   endc = triangulation.end();
   * for (; cell!=endc; ++cell)
   *   cell->set_user_flag();
   * @endcode
   * 使用C++11的基于范围的for循环，现在这完全等同于下面的内容。
   * @code
   * Triangulation<dim> triangulation;
   * ...
   * for (const auto &cell : triangulation.active_cell_iterators())
   *   cell->set_user_flag();
   * @endcode
   * @return  半开范围<code>[this->begin_active(), this->end())</code>。
   * @ingroup CPP11
   *
   */
  IteratorRange<active_cell_iterator>
  active_cell_iterators() const;

  /**
   * 返回一个迭代器范围，该范围包含了构成该三角图给定层次的所有单元（无论是否激活）。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @param[in]  level
   * 这个三角结构的细化层次中的一个给定的级别。
   * @return  半开范围<code>[this->begin(level), this->end(level))</code>
   * @pre  level必须小于this->n_levels()。
   * @ingroup CPP11
   *
   */
  IteratorRange<cell_iterator>
  cell_iterators_on_level(const unsigned int level) const;

  /**
   * 返回一个迭代器范围，该范围包含所有构成该三角形的给定级别的活动单元。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * 参见active_cell_iterators()文档中的例子。      @param[in]  level
   * 这个三角形的细化层次中的一个给定级别。    @return
   * 半开范围<code>[this->begin_active(level), this->end(level))</code>
   * @pre  level必须小于this->n_levels()。
   * @ingroup CPP11
   *
   */
  IteratorRange<active_cell_iterator>
  active_cell_iterators_on_level(const unsigned int level) const;

  /*  @} 。  
*
*/

   /*-------------------------------------------------------------------------*/ 

  /**
   * @name  面对迭代器函数  
     * @{ 
   *
   */

  /**
   * 迭代到第一个使用的面。
   *
   */
  face_iterator
  begin_face() const;

  /**
   * 迭代到第一个活动面。
   *
   */
  active_face_iterator
  begin_active_face() const;

  /**
   * 超过终点的迭代器；这个迭代器用于比较超过终点或开始前状态的迭代器。
   *
   */
  face_iterator
  end_face() const;

  /**
   * 返回一个迭代器范围，其中包含构成这个三角形的所有活动面。这个函数是
   * Triangulation::active_cell_iterators(),
   * 的面的版本，允许人们写代码，例如。
   * @code
   * Triangulation<dim> triangulation;
   * ...
   * for (auto &face : triangulation.active_face_iterators())
   *   face->set_manifold_id(42);
   * @endcode
   * @return  半开范围<code>[this->begin_active_face(),
   * this->end_face())</code>。
   * @ingroup CPP11
   *
   */
  IteratorRange<active_face_iterator>
  active_face_iterators() const;

  /*  @} .   
*
*/

   /*-------------------------------------------------------------------------*/ 

  /**
   * @name  顶点迭代器函数  @{
   *
   */

  /**
   * 迭代到第一个使用的顶点。这个函数只有在dim不是一个的情况下才能使用。
   *
   */
  vertex_iterator
  begin_vertex() const;

  /**
   * 到第一个活动顶点的迭代器。因为所有的顶点都是活动的，begin_vertex()和begin_active_vertex()返回同一个顶点。这个函数只有在dim不是一个的情况下才可以使用。
   *
   */
  active_vertex_iterator
  begin_active_vertex() const;

  /**
   * 过去的迭代器；这个迭代器用于比较过去或开始前状态的迭代器。这个函数只有在dim不为一的情况下才能使用。
   *
   */
  vertex_iterator
  end_vertex() const;

  /*  @} 。  
*
*/

  /**
   * @name  关于三角形的信息 
     * @{ 
   *
   */

  /**
   * 在下文中，大多数函数都提供了两个版本，有和没有描述水平的参数。有这个参数的版本只适用于描述当前三角剖分的单元的对象。例如：在二维中不能调用<tt>n_lines(level)</tt>，只能调用<tt>n_lines()</tt>，因为线在二维中是面，因此没有层次。
   *
   */

  /**
   * 返回已使用的线的总数，无论是否激活。
   *
   */
  unsigned int
  n_lines() const;

  /**
   * 返回已使用的线的总数，在 @p level. 层上是否有效。
   *
   */
  unsigned int
  n_lines(const unsigned int level) const;

  /**
   * 返回有效线路的总数。
   *
   */
  unsigned int
  n_active_lines() const;

  /**
   * 返回活动线路的总数，在级别 @p level. 上。
   *
   */
  unsigned int
  n_active_lines(const unsigned int level) const;

  /**
   * 返回已使用的四边形总数，无论是否激活。
   *
   */
  unsigned int
  n_quads() const;

  /**
   * 返回已使用的四边形总数，活跃或不活跃，在级别 @p
   * level. 。
   *
   */
  unsigned int
  n_quads(const unsigned int level) const;

  /**
   * 返回活跃的四边形总数，活跃与否。
   *
   */
  unsigned int
  n_active_quads() const;

  /**
   * 返回活跃的四边形总数，活跃或不活跃的水平 @p level.
   * 。
   *
   */
  unsigned int
  n_active_quads(const unsigned int level) const;

  /**
   * 返回已使用的六面体的总数，活跃与否。
   *
   */
  unsigned int
  n_hexs() const;

  /**
   * 返回已使用的六面体总数，在 @p 层中是否有效。
   *
   */
  unsigned int
  n_hexs(const unsigned int level) const;

  /**
   * 返回使用中的六面体总数，无论是否激活。
   *
   */
  unsigned int
  n_active_hexs() const;

  /**
   * 返回活跃的六面体总数，在 @p 层中活跃或不活跃。
   *
   */
  unsigned int
  n_active_hexs(const unsigned int level) const;

  /**
   * 返回已使用的单元格总数，激活或未激活。
   * 在一个空间维度上映射到<tt>n_lines()</tt>，以此类推。
   *
   */
  unsigned int
  n_cells() const;

  /**
   * 返回已使用的单元格总数，无论是否激活，在级别 @p
   * level.
   * 中映射到<tt>n_lines(level)</tt>在一个空间维度上，等等。
   *
   */
  unsigned int
  n_cells(const unsigned int level) const;

  /**
   * 返回活动单元格的总数。在一个空间维度上映射到<tt>n_active_lines()</tt>，依此类推。
   *
   */
  unsigned int
  n_active_cells() const;

  /**
   * 返回活动单元格的总数。对于当前的类，这与n_active_cells()相同。然而，该函数可以在派生类中被重载（例如，在
   * parallel::distributed::Triangulation)
   * 中，它可以返回一个大于当前处理器上三角形对象报告的活动单元数的值。
   *
   */
  virtual types::global_cell_index
  n_global_active_cells() const;


  /**
   * 返回层面上的活动单元总数  @p level.
   * 在一个空间维度上映射为<tt>n_active_lines(level)</tt>，依此类推。
   *
   */
  unsigned int
  n_active_cells(const unsigned int level) const;

  /**
   * 返回已使用的面的总数，无论是否激活。
   * 在二维中，其结果等于n_lines()，在三维中等于n_quads()，而在一维中则等于使用的顶点数。
   *
   */
  unsigned int
  n_faces() const;

  /**
   * 返回活动面的总数。
   * 在2D中，其结果等于n_active_lines()，在3D中等于n_active_quads()，而在1D中等于使用的顶点数。
   *
   */
  unsigned int
  n_active_faces() const;

  /**
   * 返回这个三角形的层数。
   * @note
   * 在内部，三角形以层为单位存储数据，这个数据结构中的层数可能比人们想象的要多。
   *
   * - 例如，想象一下我们刚刚通过粗化最高层而得到的三角形，这样它就完全被删除了。这个层次并没有被删除，因为它很可能很快就会被下一个细化过程重新填充。  因此，如果你碰巧跑过原始单元格迭代器（作为这个类的用户你不能这样做，但在内部可以），那么层次结构中的对象数量就会大于最细化的单元格的层次加一。另一方面，由于这个类的用户很少关心这个问题，这个函数实际上只是返回最细化的活动单元的级别加一。(加一是因为在一个粗略的、未精炼的网格中，所有的单元都是零级的。
   *
   * - 使得级别的数量等于1）。)
   *
   */
  unsigned int
  n_levels() const;

  /**
   * 返回正在使用的层数。这个函数等同于串行三角网格的n_levels()，但是对于
   * parallel::distributed::Triangulation
   * 来说，它给出了所有处理器上n_levels()的最大值，因此可以比n_levels()大。
   *
   */
  virtual unsigned int
  n_global_levels() const;

  /**
   * 如果三角形有悬空节点，返回true。
   * 这个函数是虚拟的，因为其结果可以有不同的解释，这取决于三角剖分是否只存在于单个处理器上，或者像
   * parallel::distributed::Triangulation
   * 类中所做的那样，可能是分布式的（关于这个函数在并行情况下应该做什么，见那里的描述）。
   *
   */
  virtual bool
  has_hanging_nodes() const;

  /**
   * 返回顶点的总数。
   * 其中一些顶点可能没有被使用，这通常发生在粗化三角形时，一些顶点被丢弃，但我们不想对剩余的顶点重新编号，导致使用的顶点数量出现漏洞。
   * 你可以使用 @p n_used_vertices 函数获得已用顶点的数量。
   *
   */
  unsigned int
  n_vertices() const;

  /**
   * 返回一个对该三角结构中所有顶点的常数引用。注意，这个数组中不一定所有的顶点都被实际使用；例如，如果你粗化一个网格，那么一些顶点会被删除，但是它们在这个数组中的位置是不变的，因为顶点的索引只被分配一次。你可以通过函数get_used_vertices()来了解哪些顶点被实际使用。
   *
   */
  const std::vector<Point<spacedim>> &
  get_vertices() const;

  /**
   * 返回目前正在使用的顶点的数量，即至少属于一个被使用的元素。
   *
   */
  unsigned int
  n_used_vertices() const;

  /**
   * 如果具有此 @p index 的顶点被使用，返回 @p true 。
   *
   */
  bool
  vertex_used(const unsigned int index) const;

  /**
   * 返回一个对 @p bools
   * 数组的常数引用，表明顶点数组中的一个条目是否被使用。
   *
   */
  const std::vector<bool> &
  get_used_vertices() const;

  /**
   * 返回在一个共同顶点相遇的最大单元格数。由于这个数字在细化过程中是一个不变量，所以只考虑最粗层的单元。因此该操作是相当快的。该不变性只适用于最粗的三角形中足够多的单元格（例如，对于单个单元格会返回一个），因此在二维空间中会返回最小的4个，在三维空间中会返回8个，等等，这就是如果三角形被细化，会有多少单元格相遇。
   * 在一个空间维度上，会返回两个。
   *
   */
  unsigned int
  max_adjacent_cells() const;

  /**
   * 这个函数总是返回 @p invalid_subdomain_id
   * ，但是为了与派生的 @p parallel::distributed::Triangulation
   * 类兼容而存在。对于分布式并行三角计算，该函数返回由当前处理器拥有的那些单元的子域ID。
   *
   */
  virtual types::subdomain_id
  locally_owned_subdomain() const;

  /**
   * 返回一个对当前对象的引用。    这似乎不是很有用，但允许编写代码，可以访问任何满足 @ref ConceptMeshType "MeshType概念 "
   * 的底层三角形（这可能不仅是一个三角形，也可能是一个DoFHandler，例如）。
   *
   */
  Triangulation<dim, spacedim> &
  get_triangulation();

  /**
   * 返回一个对当前对象的引用。这是前一个函数的const-version。
   *
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;


  /*  @} .   
*
*/

  /**
   * @name  有关对象数量的内部信息  
     * @{ 
   *
   */

  /**
   * 已使用或未使用的总行数。
   * @note
   * 这个函数实际上是输出关于三角测量的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_lines() const;

  /**
   * 在给定级别上，已使用或未使用的行数。
   * @note
   * 这个函数实际上是输出关于三角形的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_lines(const unsigned int level) const;

  /**
   * 已使用或未使用的四边形总数。
   * @note
   * 这个函数实际上是在输出关于三角形的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_quads() const;

  /**
   * 在给定的水平上，已使用或未使用的四边形的数量。
   * @note
   * 这个函数实际上是输出关于三角测量的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_quads(const unsigned int level) const;

  /**
   * 在给定级别上，已使用或未使用的赫兹数。
   * @note
   * 这个函数实际上是输出关于三角测量的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_hexs(const unsigned int level) const;

  /**
   * 给定层次上已使用或未使用的单元格的数量。
   * @note
   * 这个函数实际上是输出关于三角形的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用于其他一些非常紧密地建立在它之上的类（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_cells(const unsigned int level) const;

  /**
   * 返回面的总数，无论是否使用。在2D中，这个结果等于n_raw_lines()，在3D中等于n_raw_quads()，而在1D中等于顶点的数量。
   * @note
   * 这个函数实际上是输出关于三角形的内部信息。它不应该在应用程序中使用。这个函数只是这个类的公共接口的一部分，因为它被用在其他一些非常紧密地建立在它之上的类中（特别是DoFHandler类）。
   *
   */
  unsigned int
  n_raw_faces() const;

  /*  @} 。  
*
*/

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为三角测量对象可能通过这个基类的指针被访问，即使实际对象是一个派生类。
   *
   */
  virtual std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
   * @note
   * 这个函数不保存<i>all</i>当前三角结构的成员变量。相反，只有某些种类的信息被保存。更多信息请参见该类的一般文档。
   *
   */
  template <class Archive>
  void
  save(Archive &ar, const unsigned int version) const;

  /**
   * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，从一个流中读取此对象的数据。
   * 扔掉之前的内容。
   * @note
   * 这个函数不会将当前三角结构的<i>all</i>成员变量重置为之前存储到存档的三角结构的成员变量。相反，只有某些种类的信息被加载。更多信息请参见该类的一般文档。
   * @note  这个函数调用 Triangulation::clear() 函数，因此触发了
   * "清除 "信号。在从存档中加载所有数据后，它又触发了
   * "创建
   * "信号。关于信号的更多信息，请参见该类的一般文档。
   *
   */
  template <class Archive>
  void
  load(Archive &ar, const unsigned int version);


  /**
   * 将这个函数的参数中给出的（粗略的）面对宣布为周期性的。这样就有可能获得跨越周期性边界的邻居。
   * 该向量可以由函数 GridTools::collect_periodic_faces.
   * 填充。关于周期性边界条件的更多信息，请参见
   * GridTools::collect_periodic_faces,
   * DoFTools::make_periodicity_constraints 和 step-45  。
   * @note
   * 在使用这个函数之前，必须先初始化三角结构，而且不能进行精化。
   *
   */
  virtual void
  add_periodicity(
    const std::vector<GridTools::PeriodicFacePair<cell_iterator>> &);

  /**
   * 返回 periodic_face_map。
   *
   */
  const std::map<
    std::pair<cell_iterator, unsigned int>,
    std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3>>> &
  get_periodic_face_map() const;

  /**
   * 返回填充有该三角图所用参考单元类型的向量。
   *
   */
  const std::vector<ReferenceCell> &
  get_reference_cells() const;

  /**
   * 指示三角形是否只由类似超立方体的单元组成，即线、四边形或六面体。
   *
   */
  bool
  all_reference_cells_are_hyper_cube() const;

#ifdef DOXYGEN
  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)，从一个流中写入和读取这个对象的数据，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &archive, const unsigned int version);
#else
  // This macro defines the serialize() method that is compatible with
  // the templated save() and load() method that have been implemented.
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  /**
   * @name  异常情况  @{
   *
   */

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcInvalidLevel,
                 int,
                 int,
                 << "You are requesting information from refinement level "
                 << arg1
                 << " of a triangulation, but this triangulation only has "
                 << arg2 << " refinement levels. The given level " << arg1
                 << " must be *less* than " << arg2 << ".");
  /**
   * 引发该异常的函数只能对空三角图进行操作，即没有网格单元的三角图。
   * @ingroup Exceptions
   *
   */
  DeclException2(
    ExcTriangulationNotEmpty,
    int,
    int,
    << "You are trying to perform an operation on a triangulation "
    << "that is only allowed if the triangulation is currently empty. "
    << "However, it currently stores " << arg1 << " vertices and has "
    << "cells on " << arg2 << " levels.");
  /**
   * 试图重新读取一个网格，发生了错误。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcGridReadError);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcFacesHaveNoLevel);
  /**
   * 你所访问的三角测量层是空的。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcEmptyLevel,
                 int,
                 << "You tried to do something on level " << arg1
                 << ", but this level is empty.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNonOrientableTriangulation);

  /**
   * 异常 请求的边界_id没有找到
   *
   */
  DeclException1(ExcBoundaryIdNotFound,
                 types::boundary_id,
                 << "The given boundary_id " << arg1
                 << " is not defined in this Triangulation!");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(
    ExcInconsistentCoarseningFlags,
    "A cell is flagged for coarsening, but either not all of its siblings "
    "are active or flagged for coarsening as well. Please clean up all "
    "coarsen flags on your triangulation via "
    "Triangulation::prepare_coarsening_and_refinement() beforehand!");

  /*  @}  * */

protected:
  /**
   * 异常
   *
   */

protected:
  /**
   * 在细化三角形的过程中做一些平滑处理。关于这方面的更多信息，请参见该类的通用文档。
   *
   */
  MeshSmoothing smooth_grid;

  /**
   * 矢量缓存给定三角形的所有参考单元类型（也是在分布式情况下）。
   *
   */
  std::vector<ReferenceCell> reference_cells;

  /**
   * 向给定的流写入一个bool向量，写入一个前缀和一个后缀的魔法数字。该向量以几乎二进制的格式写入，即bool标志被打包，但数据被写成ASCII文本。
   * 标志以二进制格式存储：对于每个 @p true, ，存储一个 @p
   * 1 位，否则存储一个 @p 0 位。
   * 这些位被存储为<tt>无符号字符</tt>，从而避免了字节数。它们被写入
   * @p out
   * 的纯文本中，因此平均每输入一个比特就有3.6个比特在输出中。其他信息（神奇的数字和输入向量的元素数）也是以纯文本形式存储。因此，该格式应该是跨平台兼容的。
   *
   */
  static void
  write_bool_vector(const unsigned int       magic_number1,
                    const std::vector<bool> &v,
                    const unsigned int       magic_number2,
                    std::ostream &           out);

  /**
   * 重新读取之前由 @p write_bool_vector
   * 编写的一个布尔矢量，并与神奇数字进行比较。
   *
   */
  static void
  read_bool_vector(const unsigned int magic_number1,
                   std::vector<bool> &v,
                   const unsigned int magic_number2,
                   std::istream &     in);

  /**
   * 从periodic_face_pairs_level_0重新创建周期性邻居的信息。
   *
   */
  void
  update_periodic_face_map();

  /**
   * 更新内部reference_cells向量。
   *
   */
  virtual void
  update_reference_cells();


private:
  /**
   * 与创建、细化和粗化相关的三角形具体任务的政策。
   *
   */
  std::unique_ptr<
    dealii::internal::TriangulationImplementation::Policy<dim, spacedim>>
    policy;

  /**
   * 如果调用add_periodicity()，该变量将给定的周期性面对存储在第0层，以便以后在识别多网格层次的鬼单元和设置period_face_map时访问。
   *
   */
  std::vector<GridTools::PeriodicFacePair<cell_iterator>>
    periodic_face_pairs_level_0;

  /**
   * 如果调用add_periodicity()，这个变量将存储活动的周期性面对。
   *
   */
  std::map<std::pair<cell_iterator, unsigned int>,
           std::pair<std::pair<cell_iterator, unsigned int>, std::bitset<3>>>
    periodic_face_map;

  /**
   * @name  内部使用的细胞迭代器函数  
     * @{ 
   *
   */

  /**
   * 声明一些原始迭代器的迭代器类型，即迭代器也可以迭代在以前的网格细化周期中被粗化掉的单元格列表中的孔。
   * 由于用户不应该访问我们存储数据的这些内部属性，所以这些迭代器类型被定为私有。
   *
   */
  using raw_cell_iterator = TriaRawIterator<CellAccessor<dim, spacedim>>;
  using raw_face_iterator =
    TriaRawIterator<TriaAccessor<dim - 1, dim, spacedim>>;
  using raw_vertex_iterator =
    TriaRawIterator<dealii::TriaAccessor<0, dim, spacedim>>;
  using raw_line_iterator = typename IteratorSelector::raw_line_iterator;
  using raw_quad_iterator = typename IteratorSelector::raw_quad_iterator;
  using raw_hex_iterator  = typename IteratorSelector::raw_hex_iterator;

  /**
   * 迭代器到第一个单元，无论是否使用，在关卡上  @p
   * level.
   * 如果一个关卡没有单元，将返回一个过去的迭代器。
   *
   */
  raw_cell_iterator
  begin_raw(const unsigned int level = 0) const;

  /**
   * 返回一个原始迭代器，该迭代器是不在关卡上的第一个迭代器。如果
   * @p 级别是最后一个级别，那么这将返回<tt>end()</tt>。
   *
   */
  raw_cell_iterator
  end_raw(const unsigned int level) const;

  /*  @} 。  
*
*/

  /**
   * @name  内部使用的行迭代器函数  @{
   *
   */

  /**
   * 到第一行的迭代器，无论是否使用，在级别上  @p level.
   * 如果一个级别没有行，将返回一个过去的迭代器。
   * 如果行没有单元格，即对于 @p dim>1 来说，必须给出 @p
   * level 参数。 当然，这也适用于上述所有其他函数。
   *
   */
  raw_line_iterator
  begin_raw_line(const unsigned int level = 0) const;

  /**
   * 迭代器到第 @p level. 层的第一个使用行。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  line_iterator
  begin_line(const unsigned int level = 0) const;

  /**
   * 迭代器到第一层的活动线  @p level.  。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  active_line_iterator
  begin_active_line(const unsigned int level = 0) const;

  /**
   * 过去的迭代器；这个迭代器用于比较具有过去或开始前状态的迭代器。
   *
   */
  line_iterator
  end_line() const;

  /*  @} .   
*
*/

  /**
   * @name  内部使用的四重迭代器函数  @{
   *
   */

  /**
   * 迭代到给定关卡上的第一个四边形，无论是否使用。如果一个级别没有四边形，则返回一个过去式的迭代器。
   * 如果四边形没有单元格，即对于 $dim>2$
   * 来说，不需要给出级别参数。
   * @note  给出的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  raw_quad_iterator
  begin_raw_quad(const unsigned int level = 0) const;

  /**
   * 迭代到第一层使用的四边形 @p level. 。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  quad_iterator
  begin_quad(const unsigned int level = 0) const;

  /**
   * 迭代到第一层的活动四边形  @p level.  。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的值。
   *
   */
  active_quad_iterator
  begin_active_quad(const unsigned int level = 0) const;

  /**
   * 过了终点的迭代器；这个迭代器用于比较具有过了终点或开始前状态的迭代器。
   *
   */
  quad_iterator
  end_quad() const;

  /*  @} .   
*
*/

  /**
   * @name  内部使用的十六进制迭代器函数  
     * @{ 
   *
   */

  /**
   * 迭代器到第一个十六进制，无论是否使用，在级别上  @p
   * level.
   * 如果一个级别没有十六进制，将返回一个过去的迭代器。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  raw_hex_iterator
  begin_raw_hex(const unsigned int level = 0) const;

  /**
   * 迭代器到第一层使用的六边形 @p level. 。
   * @note  给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  hex_iterator
  begin_hex(const unsigned int level = 0) const;

  /**
   * 迭代到第一层的活动六边形 @p level. 。
   * @note 给定的 @p level
   * 参数需要对应于三角形的一个层次，也就是说，应该小于n_levels()返回的值。另一方面，对于使用
   * parallel::distributed::Triangulation
   * 对象的并行计算，在全局网格的所有层次的单元上写循环通常是很方便的，即使三角形的<i>local</i>部分实际上没有高层次的单元。在这些情况下，如果
   * @p level
   * 参数小于n_global_levels()函数的返回值，则可以接受。如果给定的
   * @p level
   * 在n_levels()和n_global_levels()返回的值之间，那么在这个层次的三角形的本地部分不存在单元，该函数只是返回end()的结果。
   *
   */
  active_hex_iterator
  begin_active_hex(const unsigned int level = 0) const;

  /**
   * 过了终点的迭代器；这个迭代器用于比较具有过了终点或开始前状态的迭代器。
   *
   */
  hex_iterator
  end_hex() const;

  /*  @} .   
*
*/


  /**
   * (公共)函数clear()只有在三角图没有被其他用户订阅的情况下才会工作。现在，clear_despite_subscriptions()函数允许在有订阅的情况下清除三角图。
   * 请确保，在调用这个函数时，你知道你在做什么，因为它的使用只在非常罕见的情况下是合理的。例如，当订阅是针对最初的空三角图的，并且三角图对象希望在由于输入错误而抛出断言之前释放其内存（例如在create_triangulation()函数中）。
   *
   */
  void
  clear_despite_subscriptions();

  /**
   * 重置三角化策略。
   *
   */
  void
  reset_policy();

  /**
   * 对于所有的单元格，设置活动单元格指数，使活动单元格知道自己是第几个活动单元格，而其他所有的单元格都有一个无效的值。这个函数在网格创建、细化和序列化之后被调用。
   *
   */
  void
  reset_active_cell_indices();

  /**
   * 重置全局单元ID和全局级单元ID。
   *
   */
  void
  reset_global_cell_indices();

  /**
   * 重置单元格的顶点索引的缓存。
   *
   */
  void
  reset_cell_vertex_indices_cache();

  /**
   * 细化所有级别上的所有单元格，这些单元格之前被标记为细化。    注意，这个函数对<tt>dim=2,3</tt>使用<tt>line->user_flags</tt>，对<tt>dim=3</tt>使用<tt>quad->user_flags。    如果在创建此对象时指定了 <code>check_for_distorted_cells</code> 标志，该函数将返回一个产生了满足 @ref GlossDistorted "扭曲的单元格 "
   * 标准的单元格列表，在
   *
   */
  DistortedCellList
  execute_refinement();

  /**
   * 粗化所有被标记为粗化的单元格，或者说：删除那些所有子单元格都被标记为粗化的单元格的所有子单元格，并且其他一些约束条件也成立（见该类的一般文档）。
   *
   */
  void
  execute_coarsening();

  /**
   * 确保一个单元格的所有子单元格或没有子单元格被标记为粗化。
   *
   */
  void
  fix_coarsen_flags();

  /**
   * 将粗化单元格的唯一id转化为其索引。更多信息请参见术语表中的 @ref GlossCoarseCellId "粗放单元格ID "
   * 条目。
   * @note  对于串行和共享三角测量，id和index都是一样的。
   * 对于分布式三角形的设置，两者可能不同，因为id可能对应于一个全局id，而索引对应于一个局部id。
   * @param  coarse_cell_id 粗略单元的唯一id。    @return
   * 粗略单元在当前三角测量中的索引。
   *
   */
  virtual unsigned int
  coarse_cell_id_to_coarse_cell_index(
    const types::coarse_cell_id coarse_cell_id) const;


  /**
   * 将粗略单元的索引转换为其唯一的ID。更多信息请参见术语表中的 @ref GlossCoarseCellId "粗放单元ID "
   * 条目。
   * @note 见方法coarse_cell_id_to_coarse_cell_index()的说明。
   * @param  coarse_cell_index 粗放单元的索引。    @return
   * 粗体单元的Id。
   *
   */
  virtual types::coarse_cell_id
  coarse_cell_index_to_coarse_cell_id(
    const unsigned int coarse_cell_index) const;

  /**
   * 指向在不同层次上存储单元数据的对象的指针阵列。
   *
   */
  std::vector<
    std::unique_ptr<dealii::internal::TriangulationImplementation::TriaLevel>>
    levels;

  /**
   * 指向三角剖分的面的指针。在1d中，它不包含任何内容，在2D中，它包含关于线的数据，在3D中，包含四边形和线。
   * 所有这些都没有水平，因此被分开处理。
   *
   */
  std::unique_ptr<dealii::internal::TriangulationImplementation::TriaFaces>
    faces;


  /**
   * 该三角形的顶点数组。
   *
   */
  std::vector<Point<spacedim>> vertices;

  /**
   * 存储顶点使用的比特模式的数组。
   *
   */
  std::vector<bool> vertices_used;

  /**
   * 流形对象的集合。我们只存储不属于FlatManifold类型的对象。
   *
   */
  std::map<types::manifold_id, std::unique_ptr<const Manifold<dim, spacedim>>>
    manifold;

  /**
   * 表示是否进行了各向异性的细化的标志。
   *
   */
  bool anisotropic_refinement;


  /**
   * 一个决定我们是否在创建和细化网格时检查扭曲的单元的标志。
   *
   */
  const bool check_for_distorted_cells;

  /**
   * 用于保存线、四边形、六边形等数字的缓存。这些数字是在细化和粗化功能结束时设置的，可以使以后的访问更加快速。在过去，每当人们想要访问这些数字之一时，就必须在所有的行上执行一个循环，例如，计算元素，直到我们碰到终端迭代器。这很耗时，而且由于访问行数等是一个相当频繁的操作，这并不是一个最佳的解决方案。
   *
   */
  dealii::internal::TriangulationImplementation::NumberCache<dim> number_cache;

  /**
   * 一个将边界顶点的数量与边界指标联系起来的映射。这个字段只在1d中使用。我们有这个字段是因为我们在2D和更高的版本中用面来存储边界指示器信息，在那里我们有存储面的数据的结构空间，但在1D中没有这样的面的空间。
   * 这个字段被声明为指针是出于一个相当平凡的原因：这个类的所有其他可以被TriaAccessor层次结构修改的字段都是指针，所以这些访问器类存储了一个指向三角形的常数指针。如果这个字段（可以被
   * TriaAccessor::set_boundary_id)
   * 修改的字段不是指针，我们就不能再为TriaAccessor<0,1,spacedim>这样做。
   *
   */
  std::unique_ptr<std::map<unsigned int, types::boundary_id>>
    vertex_to_boundary_id_map_1d;


  /**
   * 一个将边界顶点的数量与流形指标联系起来的映射。这个字段只在1d中使用。我们之所以有这个字段，是因为在2d及以上版本中，我们将流形指示器信息与面孔一起存储，在这些结构中，我们有存储面孔数据的空间，但在1d中，没有为面孔提供这样的空间。
   * @note
   * 流形对象对于点来说是非常无用的，因为它们既没有被细化，也没有被映射到其内部。然而，我们允许为点存储流形ID，以便在独立于维度的程序中保持一致。
   * 这个字段被声明为指针，是出于一个相当平凡的原因：这个类中所有可以被TriaAccessor层次结构修改的其他字段都是指针，所以这些访问器类存储了一个指向三角的常量指针。如果这个字段（可以被
   * TriaAccessor::set_manifold_id)
   * 修改的字段不是指针，我们就不能再为TriaAccessor<0,1,spacedim>这样做。
   *
   */
  std::unique_ptr<std::map<unsigned int, types::manifold_id>>
    vertex_to_manifold_id_map_1d;

  // make a couple of classes friends
  template <int, int, int>
  friend class TriaAccessorBase;
  template <int, int, int>
  friend class TriaAccessor;
  friend class TriaAccessor<0, 1, spacedim>;

  friend class CellAccessor<dim, spacedim>;

  friend struct dealii::internal::TriaAccessorImplementation::Implementation;

  friend struct dealii::internal::TriangulationImplementation::Implementation;
  friend struct dealii::internal::TriangulationImplementation::
    ImplementationMixedMesh;

  friend class dealii::internal::TriangulationImplementation::TriaObjects;

  // explicitly check for sensible template arguments, but not on windows
  // because MSVC creates bogus warnings during normal compilation
#ifndef DEAL_II_MSVC
  static_assert(dim <= spacedim,
                "The dimension <dim> of a Triangulation must be less than or "
                "equal to the space dimension <spacedim> in which it lives.");
#endif
};


#ifndef DOXYGEN



namespace internal
{
  namespace TriangulationImplementation
  {
    template <class Archive>
    void
    NumberCache<1>::serialize(Archive &ar, const unsigned int)
    {
      ar &n_levels;
      ar &n_lines &n_lines_level;
      ar &n_active_lines &n_active_lines_level;
    }


    template <class Archive>
    void
    NumberCache<2>::serialize(Archive &ar, const unsigned int version)
    {
      this->NumberCache<1>::serialize(ar, version);

      ar &n_quads &n_quads_level;
      ar &n_active_quads &n_active_quads_level;
    }


    template <class Archive>
    void
    NumberCache<3>::serialize(Archive &ar, const unsigned int version)
    {
      this->NumberCache<2>::serialize(ar, version);

      ar &n_hexes &n_hexes_level;
      ar &n_active_hexes &n_active_hexes_level;
    }

  } // namespace TriangulationImplementation
} // namespace internal


template <int dim, int spacedim>
inline bool
Triangulation<dim, spacedim>::vertex_used(const unsigned int index) const
{
  AssertIndexRange(index, vertices_used.size());
  return vertices_used[index];
}



template <int dim, int spacedim>
inline unsigned int
Triangulation<dim, spacedim>::n_levels() const
{
  return number_cache.n_levels;
}

template <int dim, int spacedim>
inline unsigned int
Triangulation<dim, spacedim>::n_global_levels() const
{
  return number_cache.n_levels;
}


template <int dim, int spacedim>
inline unsigned int
Triangulation<dim, spacedim>::n_vertices() const
{
  return vertices.size();
}



template <int dim, int spacedim>
inline const std::vector<Point<spacedim>> &
Triangulation<dim, spacedim>::get_vertices() const
{
  return vertices;
}


template <int dim, int spacedim>
template <class Archive>
void
Triangulation<dim, spacedim>::save(Archive &ar, const unsigned int) const
{
  // as discussed in the documentation, do not store the signals as
  // well as boundary and manifold description but everything else
  ar &smooth_grid;

  unsigned int n_levels = levels.size();
  ar &         n_levels;
  for (const auto &level : levels)
    ar &level;

  // boost dereferences a nullptr when serializing a nullptr
  // at least up to 1.65.1. This causes problems with clang-5.
  // Therefore, work around it.
  bool faces_is_nullptr = (faces.get() == nullptr);
  ar & faces_is_nullptr;
  if (!faces_is_nullptr)
    ar &faces;

  ar &vertices;
  ar &vertices_used;

  ar &anisotropic_refinement;
  ar &number_cache;

  ar &check_for_distorted_cells;

  if (dim == 1)
    {
      ar &vertex_to_boundary_id_map_1d;
      ar &vertex_to_manifold_id_map_1d;
    }
}



template <int dim, int spacedim>
template <class Archive>
void
Triangulation<dim, spacedim>::load(Archive &ar, const unsigned int)
{
  // clear previous content. this also calls the respective signal
  clear();

  // as discussed in the documentation, do not store the signals as
  // well as boundary and manifold description but everything else
  ar &smooth_grid;

  unsigned int size;
  ar &         size;
  levels.resize(size);
  for (auto &level_ : levels)
    {
      std::unique_ptr<internal::TriangulationImplementation::TriaLevel> level;
      ar &                                                              level;
      level_ = std::move(level);
    }

  // Workaround for nullptr, see in save().
  bool faces_is_nullptr = true;
  ar & faces_is_nullptr;
  if (!faces_is_nullptr)
    ar &faces;

  ar &vertices;
  ar &vertices_used;

  ar &anisotropic_refinement;
  ar &number_cache;

  // the levels do not serialize the active_cell_indices because
  // they are easy enough to rebuild upon re-loading data. do
  // this here. don't forget to first resize the fields appropriately
  {
    for (auto &level : levels)
      {
        level->active_cell_indices.resize(level->refine_flags.size());
        level->global_active_cell_indices.resize(level->refine_flags.size());
        level->global_level_cell_indices.resize(level->refine_flags.size());
      }
    reset_cell_vertex_indices_cache();
    reset_active_cell_indices();
    reset_global_cell_indices();
  }

  reset_policy();

  bool my_check_for_distorted_cells;
  ar & my_check_for_distorted_cells;

  Assert(my_check_for_distorted_cells == check_for_distorted_cells,
         ExcMessage("The triangulation loaded into here must have the "
                    "same setting with regard to reporting distorted "
                    "cell as the one previously stored."));

  if (dim == 1)
    {
      ar &vertex_to_boundary_id_map_1d;
      ar &vertex_to_manifold_id_map_1d;
    }

  // trigger the create signal to indicate
  // that new content has been imported into
  // the triangulation
  signals.create();
}



template <int dim, int spacedim>
inline unsigned int
Triangulation<dim, spacedim>::coarse_cell_id_to_coarse_cell_index(
  const types::coarse_cell_id coarse_cell_id) const
{
  return coarse_cell_id;
}



template <int dim, int spacedim>
inline types::coarse_cell_id
Triangulation<dim, spacedim>::coarse_cell_index_to_coarse_cell_id(
  const unsigned int coarse_cell_index) const
{
  return coarse_cell_index;
}



 /* -------------- declaration of explicit specializations ------------- */ 

template <>
unsigned int
Triangulation<1, 1>::n_quads() const;
template <>
unsigned int
Triangulation<1, 1>::n_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 1>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<2, 2>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<3, 3>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<3, 3>::n_raw_quads() const;
template <>
unsigned int
Triangulation<1, 1>::n_active_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 1>::n_active_quads() const;
template <>
unsigned int
Triangulation<1, 1>::n_raw_hexs(const unsigned int level) const;
template <>
unsigned int
Triangulation<3, 3>::n_raw_hexs(const unsigned int level) const;
template <>
unsigned int
Triangulation<3, 3>::n_hexs() const;
template <>
unsigned int
Triangulation<3, 3>::n_active_hexs() const;
template <>
unsigned int
Triangulation<3, 3>::n_active_hexs(const unsigned int) const;
template <>
unsigned int
Triangulation<3, 3>::n_hexs(const unsigned int level) const;

template <>
unsigned int
Triangulation<1, 1>::max_adjacent_cells() const;


// -------------------------------------------------------------------
// -- Explicit specializations for codimension one grids


template <>
unsigned int
Triangulation<1, 2>::n_quads() const;
template <>
unsigned int
Triangulation<1, 2>::n_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 2>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<2, 3>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 2>::n_raw_hexs(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 2>::n_active_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 2>::n_active_quads() const;
template <>
unsigned int
Triangulation<1, 2>::max_adjacent_cells() const;

// -------------------------------------------------------------------
// -- Explicit specializations for codimension two grids

template <>
unsigned int
Triangulation<1, 3>::n_quads() const;
template <>
unsigned int
Triangulation<1, 3>::n_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 3>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<2, 3>::n_raw_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 3>::n_raw_hexs(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 3>::n_active_quads(const unsigned int level) const;
template <>
unsigned int
Triangulation<1, 3>::n_active_quads() const;
template <>
unsigned int
Triangulation<1, 3>::max_adjacent_cells() const;

template <>
bool
Triangulation<1, 1>::prepare_coarsening_and_refinement();
template <>
bool
Triangulation<1, 2>::prepare_coarsening_and_refinement();
template <>
bool
Triangulation<1, 3>::prepare_coarsening_and_refinement();


extern template class Triangulation<1, 1>;
extern template class Triangulation<1, 2>;
extern template class Triangulation<1, 3>;
extern template class Triangulation<2, 2>;
extern template class Triangulation<2, 3>;
extern template class Triangulation<3, 3>;

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

// Include tria_accessor.h here, so that it is possible for an end
// user to use the iterators of Triangulation<dim> directly without
// the need to include tria_accessor.h separately. (Otherwise the
// iterators are an 'opaque' or 'incomplete' type.)
#include <deal.II/grid/tria_accessor.h>

#endif


