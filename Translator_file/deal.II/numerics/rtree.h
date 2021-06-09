//include/deal.II-translator/numerics/rtree_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_numerics_rtree_h
#define dealii_numerics_rtree_h

#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/std_cxx20/iota_view.h>

#include <deal.II/boost_adaptors/bounding_box.h>
#include <deal.II/boost_adaptors/point.h>
#include <deal.II/boost_adaptors/segment.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/strategies.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <memory>


DEAL_II_NAMESPACE_OPEN

/**
 * boost::geometry::index::rtree
 * 类的一个封装器，实现了一个自平衡的空间索引（R-tree），能够使用不同的平衡算法来存储各种类型的值。
 * 来自[Wikipedia]（https://en.wikipedia.org/wiki/R-tree）。<blockquote>
 * R树是用于空间访问方法的树形数据结构，即用于索引多维信息，如地理坐标、矩形或多边形。R树是由Antonin
 * Guttman在1984年提出的，并在理论和应用方面都有很大的用途。R树在现实世界中的一个常见用途可能是存储空间对象，如餐厅位置或典型地图所构成的多边形：街道、建筑、湖泊轮廓、海岸线等，然后快速找到查询的答案，如
 * "找到我当前位置2公里内的所有博物馆"、"检索我位置2公里内的所有路段"（在导航系统中显示）或
 * "找到最近的加油站"（尽管没有考虑道路）。R-树还可以加速各种距离指标的近邻搜索，包括大圆距离。
 * 该数据结构的关键思想是将附近的物体分组，并在树的下一个较高层次中用它们的最小边界矩形来表示；R-树中的
 * "R
 * "是指矩形。由于所有的对象都在这个边界矩形内，一个不与边界矩形相交的查询也不能与任何包含的对象相交。在叶子层，每个矩形描述了一个单一的对象；在更高的层次上，聚集了越来越多的对象。这也可以看作是对数据集越来越粗的近似。
 * R-tree的关键困难是建立一个高效的树，一方面是平衡的（所以叶子节点在同一高度）另一方面是矩形不覆盖太多空隙，不重叠太多（这样在搜索过程中，需要处理的子树就少）。例如，插入元素以获得有效的树的原始想法是，总是插入到需要最少扩大其边界盒的子树中。一旦该页满了，数据就被分成两组，每组应覆盖最小的面积。大多数对R树的研究和改进都是为了改进树的构建方式，可以分为两个目标：从头开始构建一个高效的树（称为批量加载）和对现有的树进行修改（插入和删除）。</blockquote>
 * 一个RTree可以存储任何类型的 @p LeafType
 * ，只要可以提取一个RTree可以处理的 @p Indexable
 * ，并进行数值比较。一个 @p Indexable
 * 是适应于点（Point）、边界盒（BoundingBox）或段（Segment）概念的类型，对它来说，距离和相等的比较是可以实现的。deal.II的Point、Segment和BoundingBox类满足这个要求，但你可以混入任何
 * boost::geometry 接受为可索引的几何对象。
 * 特别是，给定一个 @p Indexable
 * 类型（例如一个点，一个BoundingBox，或一个段）， @p
 * LeafType 可以由 @p Indexable,   `std::pair<Indexable,  T>`，
 * `boost::tuple<Indexable,  ...>`或 `std::tuple<Indexable,
 * ...>`中的任何一个。 可选的参数 @p IndexType
 * 只在向树上逐一添加元素时使用。如果使用的是范围插入，那么树是使用打包算法建立的。
 * 如果想按顺序构建树，可以使用线性算法、二次算法和rstar算法。然而，这些都不是很有效，用户应该尽可能地使用打包算法。
 * 打包算法是一次性构建树，当你拥有所有的叶子时可以使用。
 * 这个类通常与两个辅助函数pack_rtree()中的一个结合使用，该函数接收一个容器或一系列迭代器，使用打包算法构建RTree。
 * 以下是一个使用实例。
 *
 *
 * @code
 * std::vector<Point<2>> points = generate_some_points();
 * auto tree = pack_rtree(points.begin(), points.end());
 * // or, equivalently:
 * // auto tree = pack_rtree(points);
 * @endcode
 *
 * 通过使用 [`boost::geometry::index`
 * 查询](https://www.boost.org/doc/libs/1_68_0/libs/geometry/doc/html/geometry/spatial_indexes/queries.html)来访问该树。例如，在用上面的片段构建了树之后，人们可以用下面的方式询问离某段最近的点。
 *
 *
 * @code
 * namespace bgi = boost::geometry::index;
 *
 * Segment<2> segment(Point<2>(0,0), Point<2>(1,1));
 *
 * std::vector<Point<2>> nearest;
 * tree.query(bgi::nearest(segment,3), std::back_inserter(intersection));
 * // Returns the 3 closest points to the Segment defined above.
 * @endcode
 *
 * 一般来说，树不需要存储实际的对象，只要它知道如何访问一个可索引类型的常量引用。这可以通过传递可选的模板参数
 * @p IndexableGetter,
 * 来实现，该模板参数可以提取给定类型的对象 @p LeafType.
 * 的一个可能的可索引类型的常量引用。
 * 作为一个例子，人们可以在一个容器中存储点，并且只在容器中创建一个索引树，使用下面定义的IndexableGetterFromIndices类，以及函数pack_rtree_of_indices()。
 *
 *
 */
template <typename LeafType,
          typename IndexType = boost::geometry::index::linear<16>,
          typename IndexableGetter =
            boost::geometry::index::indexable<LeafType>>
using RTree =
  boost::geometry::index::rtree<LeafType, IndexType, IndexableGetter>;

/**
 * 通过传递一个迭代器范围来构造正确的RTree对象。
 * 注意，参数的顺序与RTree类相反，因为我们可以从参数中自动推断出
 * @p LeafType 。
 *
 *
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename LeafTypeIterator,
          typename IndexableGetter = boost::geometry::index::indexable<
            typename LeafTypeIterator::value_type>>
RTree<typename LeafTypeIterator::value_type, IndexType, IndexableGetter>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end);

/**
 * 通过传递一个STL容器类型来构造一个RTree对象。这个函数在
 * step-70  中使用。
 * 注意模板参数的顺序与RTree类相反，因为我们可以从参数中自动推断出
 * @p LeafType ，而我们只需要在默认值不够的情况下指定 @p
 * IndexType 。
 *
 *
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename ContainerType,
          typename IndexableGetter = boost::geometry::index::indexable<
            typename ContainerType::value_type>>
RTree<typename ContainerType::value_type, IndexType, IndexableGetter>
pack_rtree(const ContainerType &container);

/**
 * 一个可以作为 @p IndexableGetter 模板参数的类，用于存储 @p
 * Container 类型中的条目索引的RTree。
 * 这个类可以作为代理，从兼容的容器中提取一个可索引的类型。比如说。
 *
 * @code
 * std::vector<std::pair<Point<dim>, double> > point_temperature = fill();
 * IndexableGetterFromIndices<decltype(point_temperature)>
 *  getter(point_temperature);
 *
 * const Point<dim> &p = getter(i); // returns point_temperature[i].first;
 * @endcode
 *
 * 这个类被pack_rtree_of_indices()函数用来构造一个RTree，其中的叶子是指向传递给这个类的容器的条目的指数。
 *
 *
 */
template <typename Container>
class IndexableGetterFromIndices
{
public:
  /**
   * 提升类型的别名，用于从兼容类型（对、图元等）中提取点、段或边界框。
   *
   */
  using IndexableGetter =
    typename boost::geometry::index::indexable<typename Container::value_type>;

  /**
   * 实际几何类型的别名。
   *
   */
  using result_type = typename IndexableGetter::result_type;

  /**
   * 索引类型的别名。
   *
   */
  using size_t = typename Container::size_type;

  /**
   * 构造函数。存储一个容器的常量引用。
   *
   */
  explicit IndexableGetterFromIndices(Container const &c)
    : container(c)
  {}

  /**
   * 实现了rtree类的 @p IndexableGetter 要求。
   *
   */
  result_type
  operator()(size_t i) const
  {
    return getter(container[i]);
  }

private:
  /**
   * 一个对容器的常量引用。
   *
   */
  const Container &container;

  /**
   * 一个getter的实例化，允许从容器 @p value_type
   * 和实际的可索引类型中轻松转换。
   *
   */
  IndexableGetter getter;
};

/**
 * 构建一个RTree对象，存储现有可索引类型的容器的索引。对容器的唯一要求是它支持0和容器大小之间的任何索引的operator[]（即，一个
 * std::vector, 或一个 std::array 可以，然而一个 std::map
 * 就不行了）。
 * 与pack_rtree()函数创建的对象不同，在这种情况下，我们不存储实际的几何类型，而只是存储它们的索引。
 *
 *
 * @code
 * namespace bgi = boost::geometry::index;
 * std::vector<Point<dim>> some_points = fill();
 * auto tree = pack_rtree(points); // the tree contains a copy of the points
 * auto index_tree = pack_rtree_of_indices(points); // the tree contains only
 *                                                // the indices of the
 *                                                // points
 * BoundingBox<dim> box = build_a_box();
 *
 * for(const auto &p: tree       | bgi::adaptors::queried(bgi::intersects(box)))
 * std::cout << "Point p: " << p << std::endl;
 *
 * for(const auto &i: index_tree | bgi::adaptors::queried(bgi::intersects(box)))
 * std::cout << "Point p: " << some_points[i] << std::endl;
 * @endcode
 *
 * 存储在树中的叶子是容器中相应条目的索引。对外部容器的引用被存储在内部，但请记住，如果你改变了容器，你应该重新建立树。
 *
 *
 */
template <typename IndexType = boost::geometry::index::linear<16>,
          typename ContainerType>
RTree<typename ContainerType::size_type,
      IndexType,
      IndexableGetterFromIndices<ContainerType>>
pack_rtree_of_indices(const ContainerType &container);

/**
 * 帮助结构，允许人们从RTree中提取一个层次作为BoundingBox对象的向量。
 * 这个结构实现了一个 boost::geometry::index::detail::rtree::visitor
 * 对象，允许人们访问任何现有的RTree对象，并返回与该树的特定目标层相关的边界框的向量。
 * 尽管有可能，直接使用这个结构是很麻烦的。建议通过辅助函数
 * extract_rtree_level() 来使用这个类。
 *
 *
 */
template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
struct ExtractLevelVisitor
  : public boost::geometry::index::detail::rtree::visitor<
      Value,
      typename Options::parameters_type,
      Box,
      Allocators,
      typename Options::node_tag,
      true>::type
{
  /**
   * 构建一个BoundingBox对象的向量 @p boxes ，对应于树的 @p
   * target_level 。
   *
   */
  inline ExtractLevelVisitor(
    Translator const & translator,
    const unsigned int target_level,
    std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes);

  /**
   * 一个标识树的内部节点的别名。
   *
   */
  using InternalNode =
    typename boost::geometry::index::detail::rtree::internal_node<
      Value,
      typename Options::parameters_type,
      Box,
      Allocators,
      typename Options::node_tag>::type;

  /**
   * 一个识别树的叶子的别名。
   *
   */
  using Leaf = typename boost::geometry::index::detail::rtree::leaf<
    Value,
    typename Options::parameters_type,
    Box,
    Allocators,
    typename Options::node_tag>::type;

  /**
   * 为InternalNode对象实现访问者接口。如果该节点属于 @p
   * target_leve, ，则填充边界框向量。
   *
   */
  inline void
  operator()(InternalNode const &node);

  /**
   * 执行Leaf对象的访问者接口。
   *
   */
  inline void
  operator()(Leaf const &);

  /**
   * 翻译器接口，是rtree的boost实现所需要的。
   *
   */
  Translator const &translator;

  /**
   * 存储我们当前访问的级别。
   *
   */
  size_t level;

  /**
   * 我们要从RTree对象中提取的级别。
   *
   */
  const size_t target_level;

  /**
   * 对BoundingBox对象的输入向量的引用。
   *
   */
  std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes;
};

/**
 * 给定一个RTree对象 @p rtree, 和一个目标级别 @p level,
 * ，返回一个BoundingBox对象的向量，其中包含使给定的 @p
 * level 的 @p rtree. 的所有边界框。
 * 这个函数是ExtractLevelVisitor类的一个方便的包装器。它被用于
 * step-70  。
 * 由于RTree对象是一棵平衡的树，你可以期望得到的向量的每个条目都包含大致相同数量的子对象，并最终包含相同数量的叶对象。如果你要求的级别在RTree中不存在，那么将返回最后的级别。
 * 这个函数的一个典型用法是在
 * parallel::distributed::Triangulation
 * 对象的背景下，人们希望构建一个活动进程的本地拥有的单元所覆盖的区域的粗略表示，并与其他进程交换这一信息。最精细的信息是由叶子提供的，在这种情况下，叶子是与三角形的本地所有单元相关的所有边界盒的集合。与所有参与的进程交换这些信息将违背并行计算的目的。然而，如果构建一个包含这些边界框的RT树（例如，通过调用
 * GridTools::Cache::get_cell_bounding_boxes_rtree()),
 * ，然后提取RT树的第一层之一，那么将只返回少量的BoundingBox对象，允许用户对域的几何形状及其在进程中的分布有一个非常有效的描述。
 * 下面的代码片断给出了一个使用实例。
 *
 * @code
 * parallel::distributed::Triangulation<2> tria(MPI_COMM_WORLD);
 * GridGenerator::hyper_ball(tria);
 * tria.refine_global(4);
 *
 * std::vector<BoundingBox<2>> all_boxes(tria.n_locally_owned_active_cells());
 * unsigned int                i = 0;
 * for (const auto &cell : tria.active_cell_iterators())
 * if (cell->is_locally_owned())
 *   all_boxes[i++] = cell->bounding_box();
 *
 * const auto tree  = pack_rtree(all_boxes);
 * const auto boxes = extract_rtree_level(tree, 1);
 * @endcode
 *
 * 当在三个进程上运行时，仅围绕本地拥有的单元格的BoundingBox对象的完整集合以及用这些盒子构建的rtree的第二层看起来就像下面的图片（每个进程一个图片）。
*  @image html rtree-process-0.png   @image html rtree-process-1.png   @image html rtree-process-2.png 。
 *
 */
template <typename Rtree>
inline std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>
extract_rtree_level(const Rtree &tree, const unsigned int level);



// Inline and template functions
#ifndef DOXYGEN

template <typename IndexType,
          typename LeafTypeIterator,
          typename IndexableGetter>
RTree<typename LeafTypeIterator::value_type, IndexType, IndexableGetter>
pack_rtree(const LeafTypeIterator &begin, const LeafTypeIterator &end)
{
  return RTree<typename LeafTypeIterator::value_type,
               IndexType,
               IndexableGetter>(begin, end);
}



template <typename IndexType, typename ContainerType, typename IndexableGetter>
RTree<typename ContainerType::value_type, IndexType, IndexableGetter>
pack_rtree(const ContainerType &container)
{
  return pack_rtree<IndexType, decltype(container.begin()), IndexableGetter>(
    container.begin(), container.end());
}



template <typename IndexType, typename ContainerType>
RTree<typename ContainerType::size_type,
      IndexType,
      IndexableGetterFromIndices<ContainerType>>
pack_rtree_of_indices(const ContainerType &container)
{
  std_cxx20::ranges::iota_view<typename ContainerType::size_type,
                               typename ContainerType::size_type>
    indices(0, container.size());
  return RTree<typename ContainerType::size_type,
               IndexType,
               IndexableGetterFromIndices<ContainerType>>(
    indices.begin(),
    indices.end(),
    IndexType(),
    IndexableGetterFromIndices<ContainerType>(container));
}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::
  ExtractLevelVisitor(
    const Translator & translator,
    const unsigned int target_level,
    std::vector<BoundingBox<boost::geometry::dimension<Box>::value>> &boxes)
  : translator(translator)
  , level(0)
  , target_level(target_level)
  , boxes(boxes)
{}



template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::
operator()(const ExtractLevelVisitor::InternalNode &node)
{
  using ElmentsType =
    typename boost::geometry::index::detail::rtree::elements_type<
      InternalNode>::type;

  const auto &elements = boost::geometry::index::detail::rtree::elements(node);

  if (level == target_level)
    {
      const auto offset = boxes.size();
      boxes.resize(offset + elements.size());

      unsigned int i = offset;
      for (typename ElmentsType::const_iterator it = elements.begin();
           it != elements.end();
           ++it)
        {
          boost::geometry::convert(it->first, boxes[i]);
          ++i;
        }
      return;
    }

  const size_t level_backup = level;
  ++level;

  for (typename ElmentsType::const_iterator it = elements.begin();
       it != elements.end();
       ++it)
    {
      boost::geometry::index::detail::rtree::apply_visitor(*this, *it->second);
    }

  level = level_backup;
}

template <typename Value,
          typename Options,
          typename Translator,
          typename Box,
          typename Allocators>
void
ExtractLevelVisitor<Value, Options, Translator, Box, Allocators>::
operator()(const ExtractLevelVisitor::Leaf &)
{}



template <typename Rtree>
inline std::vector<BoundingBox<
  boost::geometry::dimension<typename Rtree::indexable_type>::value>>
extract_rtree_level(const Rtree &tree, const unsigned int level)
{
  constexpr unsigned int dim =
    boost::geometry::dimension<typename Rtree::indexable_type>::value;

  using RtreeView =
    boost::geometry::index::detail::rtree::utilities::view<Rtree>;
  RtreeView rtv(tree);

  std::vector<BoundingBox<dim>> boxes;

  if (rtv.depth() == 0)
    {
      // The below algorithm does not work for `rtv.depth()==0`, which might
      // happen if the number entries in the tree is too small.
      // In this case, simply return a single bounding box.
      boxes.resize(1);
      boost::geometry::convert(tree.bounds(), boxes[0]);
    }
  else
    {
      const unsigned int target_level =
        std::min<unsigned int>(level, rtv.depth() - 1);

      ExtractLevelVisitor<typename RtreeView::value_type,
                          typename RtreeView::options_type,
                          typename RtreeView::translator_type,
                          typename RtreeView::box_type,
                          typename RtreeView::allocators_type>
        extract_level_visitor(rtv.translator(), target_level, boxes);
      rtv.apply_visitor(extract_level_visitor);
    }

  return boxes;
}



template <class Rtree>
unsigned int
n_levels(const Rtree &tree)
{
  boost::geometry::index::detail::rtree::utilities::view<Rtree> rtv(tree);
  return rtv.depth();
}



#endif

DEAL_II_NAMESPACE_CLOSE

#endif


