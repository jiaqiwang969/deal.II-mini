//include/deal.II-translator/grid/tria_accessor_0.txt
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

#ifndef dealii_tria_accessor_h
#define dealii_tria_accessor_h


#include <deal.II/base/config.h>

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/reference_cell.h>
#include <deal.II/grid/tria_iterator_base.h>
#include <deal.II/grid/tria_iterator_selector.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/container/small_vector.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <utility>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <typename Accessor>
class TriaRawIterator;
template <typename Accessor>
class TriaIterator;
template <typename Accessor>
class TriaActiveIterator;

namespace parallel
{
  template <int dim, int spacedim>
  class TriangulationBase;
}

template <int dim, int spacedim>
class Manifold;

template <int dim, int spacedim>
class Mapping;
#endif

namespace internal
{
  namespace TriangulationImplementation
  {
    class TriaObjects;
    struct Implementation;
    struct ImplementationMixedMesh;
  } // namespace TriangulationImplementation

  namespace TriaAccessorImplementation
  {
    struct Implementation;

    /**
     * 一个类型的实现，用它来存储访问器对象的级别。我们只在<tt>structdim
     * == dim</tt>的情况下需要它。
     * 否则，一个空对象就足够了。
     *
     */
    template <int structdim, int dim>
    struct PresentLevelType
    {
      struct type
      {
        /**
         * 默认构造函数。
         *
         */
        type() = default;

        /**
         * 虚构的构造函数。只允许零级。
         *
         */
        type(const int level)
        {
          Assert(level == 0, ExcInternalError());
          (void)level; // removes -Wunused-parameter warning in optimized mode
        }

        /**
         * 虚数转换操作符。返回零级。
         *
         */
        operator int() const
        {
          return 0;
        }

        void
        operator++() const
        {
          Assert(false, ExcInternalError());
        }

        void
        operator--() const
        {
          Assert(false, ExcInternalError());
        }
      };
    };


    /**
     * 实现一个类型，用它来存储访问器对象的级别。我们只在<tt>structdim
     * == dim</tt>的情况下需要它。
     * 否则，一个空对象就足够了。
     *
     */
    template <int dim>
    struct PresentLevelType<dim, dim>
    {
      using type = int;
    };
  } // namespace TriaAccessorImplementation
} // namespace internal
template <int structdim, int dim, int spacedim>
class TriaAccessor;
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>;
template <int spacedim>
class TriaAccessor<0, 1, spacedim>;

// note: the file tria_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates. originally, this was
// only done in debug mode, but led to cyclic reduction problems and
// so is now on by default.


/**
 * 一个包含访问器类所使用的异常类的命名空间。
 *
 *
 */
namespace TriaAccessorExceptions
{
  /**
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcCellNotUsed,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that point to valid "
                   "objects. These objects need not necessarily be active, "
                   "i.e., have no children, but they need to be part of a "
                   "triangulation. (The objects pointed to by an iterator "
                   "may -- after coarsening -- also be objects that used "
                   "to be part of a triangulation, but are now no longer "
                   "used. Their memory location may have been retained "
                   "for re-use upon the next mesh refinement, but is "
                   "currently unused.)");
  /**
   * 单元不是 @ref GlossActive 的 "活动 "
   * 单元，但它已经有了孩子。一些操作，如设置细化标志或访问自由度，只能在活动单元上实现。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcCellNotActive,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that point to 'active' "
                   "objects. 'Active' objects are those that do not have "
                   "children (in the case of cells), or that are part of "
                   "an active cell (in the case of faces or edges). However, "
                   "the object on which you are trying the current "
                   "operation is not 'active' in this sense.");
  /**
   * 试图访问一个事实上处于活动状态的单元的子单元。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcCellHasNoChildren,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that have children, "
                   "but the object on which you are trying the current "
                   "operation does not have any.");
  /**
   * 试图访问处于三角结构最粗层的单元的父级。
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcCellHasNoParent,
                   "The operation you are attempting can only be performed for "
                   "(cell, face, or edge) iterators that have a parent object, "
                   "but the object on which you are trying the current "
                   "operation does not have one -- i.e., it is on the "
                   "coarsest level of the triangulation.");
  /**
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcCantSetChildren,
                 int,
                 << "You can only set the child index if the cell does not "
                 << "currently have children registered; or you can clear it. "
                 << "The given index was " << arg1
                 << " (-1 means: clear children).");
  /**
   * @ingroup Exceptions
   *
   */
  template <typename AccessorType>
  DeclException1(ExcDereferenceInvalidObject,
                 AccessorType,
                 << "You tried to dereference an iterator for which this "
                 << "is not possible. More information on this iterator: "
                 << "index=" << arg1.index() << ", state="
                 << (arg1.state() == IteratorState::valid ?
                       "valid" :
                       (arg1.state() == IteratorState::past_the_end ?
                          "past_the_end" :
                          "invalid")));
  /**
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcCantCompareIterators,
                   "Iterators can only be compared if they point to the same "
                   "triangulation, or if neither of them are associated "
                   "with a triangulation.");
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNeighborIsCoarser);
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNeighborIsNotCoarser);
  /**
   * 你试图访问一个面孔的级别，但是面孔没有固有的级别。一个面的级别只能由相邻面的级别决定，这又意味着一个面可以有几个级别。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcFacesHaveNoLevel);
  /**
   * 你试图得到一个面的周期性邻居，而这个面并没有周期性邻居。关于这方面的更多信息，请参考 @ref GlossPeriodicConstraints "周期性边界的条目"
   * 。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNoPeriodicNeighbor);
  // TODO: Write documentation!
  /**
   * @ingroup Exceptions
   *
   */
  DeclException1(
    ExcSetOnlyEvenChildren,
    int,
    << "You can only set the child index of an even numbered child."
    << "The number of the child given was " << arg1 << ".");
} // namespace TriaAccessorExceptions


/**
 * TriaRawIterator和派生类所使用的访问器类的基类。
 * 该类只提供迭代器所需的基本功能（存储必要的数据成员，提供比较运算符等），但没有实际解除数据定义的功能。这是在派生类中完成的。
 * 在实现中，该类的行为在<tt>structdim==dim</tt>（网格的单元）和<tt>structdim&lt;dim</tt>（面和边）的情况下有所不同。对于后者，#present_level总是等于0，构造函数不可能在这里收到正值。对于单元格来说，任何级别都是可能的，但只有那些在Triangulation的级别范围内的级别是合理的。此外，函数objects()要么返回具有同一层次的所有单元的容器，要么返回具有该维度的所有对象的容器（<tt>structdim&lt;dim</tt>）。
 * 这个类的一些内部情况在
 * @ref IteratorAccessorInternals  中讨论。
 *
 *
 * @ingroup grid
 *
 * @ingroup Accessors
 *
 */
template <int structdim, int dim, int spacedim = dim>
class TriaAccessorBase
{
public:
  /**
   * 这个访问器所代表的对象所处空间的尺寸。
   * 例如，如果这个访问器代表一个四边形，它是四维空间中一个二维表面的一部分，那么这个值就是四。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 这个访问器所代表的事物的尺寸，是这个访问器的一部分。例如，如果这个访问器代表一条线，是六面体的一部分，那么这个值就是三。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 该访问器所代表的当前对象的尺寸。例如，如果它是线（不管它是四面体还是六面体的一部分，也不管我们处于什么维度），那么这个值就等于1。
   *
   */
  static const unsigned int structure_dimension = structdim;

  /**
   * 复制运算符。这些操作符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的情况下使用。据推测，这里的意图是将 @p b
   * 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个访问器；因此，这个操作对三角形的迭代器没有用。
   * 因此，这个操作被声明为删除，不能使用。
   *
   */
  void
  operator=(const TriaAccessorBase *) = delete;

protected:
  /**
   * 声明该访问器类期望从迭代器类中得到的数据类型。由于纯三角迭代器不需要额外的数据，这个数据类型是
   * @p void. 。
   *
   */
  using AccessorData = void;

  /**
   * 构造函数。受保护的，因此只能从好友类中调用。
   *
   */
  TriaAccessorBase(const Triangulation<dim, spacedim> *parent = nullptr,
                   const int                           level  = -1,
                   const int                           index  = -1,
                   const AccessorData *                       = nullptr);

  /**
   * 复制构造函数。用完全相同的数据创建一个对象。
   *
   */
  TriaAccessorBase(const TriaAccessorBase &);

  /**
   * 复制操作符。因为这只是从迭代器中调用，所以不要返回任何东西，因为迭代器会返回自己。
   * 这个方法是受保护的，因为它只能从迭代器类中调用。
   *
   */
  void
  copy_from(const TriaAccessorBase &);

  /**
   * 复制操作符。创建一个具有完全相同数据的对象。
   *
   */
  TriaAccessorBase &
  operator=(const TriaAccessorBase &);

  /**
   * 访问器的比较运算符。这个操作符在比较迭代器进入三角形的对象时使用，例如，当把它们放入一个
   * `std::map`.
   * 如果#structure_dimension小于#dimension，我们只是比较这样一个对象的索引，因为面和边没有层次。如果#structure_dimension等于#dimension，我们首先比较层次，只有在层次相等时才比较索引。
   *
   */
  bool
  operator<(const TriaAccessorBase &other) const;

protected:
  /**
   * 比较是否相等。
   *
   */
  bool
  operator==(const TriaAccessorBase &) const;

  /**
   * 比较不平等。
   *
   */
  bool
  operator!=(const TriaAccessorBase &) const;

  /**
   * @name  迭代器的推进
   *
   */
  /**
   * @{
   *
   */
  /**
   * 该操作符将迭代器推进到下一个元素。    仅适用于 @p
   * dim=1
   * 。如果有更多的元素，下一个元素就是这一层的下一个。如果现在的元素是这一层的最后一个，则访问下一层的第一个。
   *
   */
  void
  operator++();

  /**
   * 这个操作符将迭代器移动到前一个元素。    仅适用于 @p
   * dim=1
   * 。如果<tt>index>0</tt>，前一个元素是本层的前一个。如果现在的元素是这一层的第一个，则访问前一层的最后一个。
   *
   */
  void
  operator--();
  /**
   * @}
   *
   */

  /**
   * 访问具有相同维度的Triangulation的其他对象。
   *
   */
  dealii::internal::TriangulationImplementation::TriaObjects &
  objects() const;

public:
  /**
   * 用来从迭代器向访问器类传递参数的数据类型，无论这些参数的数量类型是什么，都是统一的。
   *
   */
  using LocalData = void *;

  /**
   * @name  迭代器地址和状态
   *
   */
  /**
   * @{
   *
   */

  /**
   * 对于单元格，该函数返回该单元格所处的网格层次。对于所有其他对象，该函数返回零。
   * @note  在一个三角形对象中，单元格由一对 <code>(level,
   * index)</code>
   * 唯一标识，前者是单元格的细化层次，后者是该细化层次中单元格的索引（前者是本函数的返回值）。
   * 因此，可能有多个处于不同细化级别的单元格，但在其级别中具有相同的索引。与此相反，如果当前对象对应的是一个面或边，那么该对象就只能通过其索引来识别，因为面和边没有细化级别。对于这些对象，当前函数总是返回0作为级别。
   *
   */
  int
  level() const;

  /**
   * 返回当前级别上指向的元素的索引。
   * 在Triangulation对象中，单元格由一对 <code>(level,
   * index)</code>
   * 唯一标识，前者是单元格的细化级别，后者是该细化级别中单元格的索引（后者是本函数的返回值）。因此，可能有多个处于不同细化级别的单元格，但在其级别中具有相同的索引。与此相反，如果当前对象对应于一个面或边，那么该对象只能通过其索引来识别，因为面和边没有细化级别。
   * @note
   * 这个函数返回的索引对象并不是每个层次上的连续数字集：从一个单元到另一个单元，一个层次中的一些索引可能没有使用。
   * @note  如果三角图实际上是 parallel::distributed::Triangulation
   * 类型的，那么指数就只相对于存储在当前处理器上的分布式三角图的那部分。换句话说，生活在不同处理器上的三角形分区中的单元可能具有相同的索引，即使它们指的是同一个单元，也可能具有不同的索引，即使它们指的是同一个单元（例如，如果一个单元被一个处理器拥有，但在另一个处理器上是一个幽灵单元）。
   *
   */
  int
  index() const;

  /**
   * 返回迭代器的状态。
   * 关于一个访问器可能处于的不同状态，请参考TriaRawIterator文档。
   *
   */
  IteratorState::IteratorStates
  state() const;

  /**
   * 返回该类所指向的对象所属于的三角形的引用。
   *
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * @}
   *
   */
protected:
  /**
   * 如果这是一个单元格（<tt>structdim==dim</tt>），则是水平。否则，包含零。
   *
   */
  typename dealii::internal::TriaAccessorImplementation::
    PresentLevelType<structdim, dim>::type present_level;

  /**
   * 用于存储目前所使用的级别上目前所指向的元素的索引。
   *
   */
  int present_index;

  /**
   * 指向我们行动的三角图的指针。
   *
   */
  const Triangulation<dim, spacedim> *tria;

private:
  template <typename Accessor>
  friend class TriaRawIterator;
  template <typename Accessor>
  friend class TriaIterator;
  template <typename Accessor>
  friend class TriaActiveIterator;
};



/**
 * 一个表示没有意义的迭代器的访问对象的类，如1d网格上的四维迭代器。
 * 这个类不能用来创建对象（事实上，如果试图这样做的话，它会抛出一个异常，但它有时允许以独立于维度的方式，以更简单的方式编写代码。例如，它允许编写独立于维度的四边形迭代器的代码
 *
 * - 即，也可以在1d中进行编译
 *
 * - 因为四元迭代器（通过当前的类）存在，并且在语法上是正确的。然而，你不能期望在1d中创建这些迭代器中的一个实际对象，这意味着你需要期望将使用四元迭代器的代码块包装成类似 <code>if (dim@>1)</code> 的东西。
 *
 * - 这反正是很有意义的。
 * 这个类提供了Accessor类与Iterator类交互所需的最小接口。然而，这只是为了语法上的正确性，这些函数除了产生错误之外，没有任何作用。
 *
 *
 * @ingroup Accessors
 *
 */
template <int structdim, int dim, int spacedim = dim>
class InvalidAccessor : public TriaAccessorBase<structdim, dim, spacedim>
{
public:
  /**
   * 从基类传播别名到这个类。
   *
   */
  using AccessorData =
    typename TriaAccessorBase<structdim, dim, spacedim>::AccessorData;

  /**
   * 构造器。
   * 这个类用于在给定维度中没有意义的迭代器，例如1D网格的四边形。因此，虽然这种对象的创建在语法上是有效的，但它们在语义上没有意义，当这种对象实际生成时，我们会产生一个异常。
   *
   */
  InvalidAccessor(const Triangulation<dim, spacedim> *parent     = nullptr,
                  const int                           level      = -1,
                  const int                           index      = -1,
                  const AccessorData *                local_data = nullptr);

  /**
   * 复制构造函数。
   * 这个类用于在给定维度中没有意义的迭代器，例如1D网格的四边形。因此，虽然这种对象的创建在语法上是有效的，但它们在语义上没有意义，当这种对象实际生成时，我们会产生一个异常。
   *
   */
  InvalidAccessor(const InvalidAccessor &);

  /**
   * 从其他访问器转换到当前无效的访问器。这当然也会导致运行时错误。
   *
   */
  template <typename OtherAccessor>
  InvalidAccessor(const OtherAccessor &);

  /**
   * 虚假的复制操作。
   *
   */
  void
  copy_from(const InvalidAccessor &);

  /**
   * 假的比较运算符。
   *
   */
  bool
  operator==(const InvalidAccessor &) const;
  bool
  operator!=(const InvalidAccessor &) const;

  /**
   * 虚构的运算符，以使事情得到编译。什么都不做。
   *
   */
  void
  operator++() const;
  void
  operator--() const;

  /**
   * 代表访问器是否指向已使用或未使用的对象的假函数。
   *
   */
  bool
  used() const;

  /**
   * 代表访问器是否指向有子对象的假函数。
   *
   */
  bool
  has_children() const;

  /**
   * 总是返回 numbers::flat_manifold_id. 的假函数。
   *
   */
  types::manifold_id
  manifold_id() const;

  /**
   * 总是返回 numbers::invalid_unsigned_int. 的假函数。
   *
   */
  unsigned int
  user_index() const;

  /**
   * 总是抛出的假函数。
   *
   */
  void
  set_user_index(const unsigned int p) const;

  /**
   * 总是抛出的假函数。
   *
   */
  void
  set_manifold_id(const types::manifold_id) const;

  /**
   * 抽取顶点的假函数。返回原点。
   *
   */
  Point<spacedim> &
  vertex(const unsigned int i) const;

  /**
   * 抽取线的假函数。返回一个默认构造的线段迭代器。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator
    line(const unsigned int i) const;

  /**
   * 抽取四边形的假函数。返回一个默认的四边形迭代器。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator
    quad(const unsigned int i) const;
};



/**
 * 一个提供对三角形中的对象的访问的类，如它的顶点、子对象、子女、几何信息等。这个类表示维度为
 * <code>structdim</code>
 * （即1代表线，2代表四边形，3代表六边形）的对象在维度为
 * <code>dim</code>
 * （即1代表线的三角结构，2代表四边形的三角结构，3代表六边形的三角结构）的三角结构中，该空间嵌入维度为
 * <code>spacedim</code>  ]（对于 <code>spacedim==dim</code>
 * ，三角形代表 $R^{dim}$ 中的一个域，对于
 * <code>spacedim@>dim</code>
 * ，三角形是嵌入高维空间中的流形）。 对于 @p structdim
 * 等于零的情况，即对于三角化的顶点，该类有一个特殊化。
 *
 *
 * @ingroup Accessors
 *
 */
template <int structdim, int dim, int spacedim>
class TriaAccessor : public TriaAccessorBase<structdim, dim, spacedim>
{
public:
  /**
   * 从基类传播别名到这个类。
   *
   */
  using AccessorData =
    typename TriaAccessorBase<structdim, dim, spacedim>::AccessorData;

  /**
   * 构造函数。
   *
   */
  TriaAccessor(const Triangulation<dim, spacedim> *parent     = nullptr,
               const int                           level      = -1,
               const int                           index      = -1,
               const AccessorData *                local_data = nullptr);

  /**
   * 复制构造函数不被删除，但复制的构造元素不应该被修改，也是对复制赋值运算符的注释。
   *
   */
  TriaAccessor(const TriaAccessor &) = default;

  /**
   * 移动构造函数。
   *
   */
  TriaAccessor(TriaAccessor &&) = default; // NOLINT

  /**
   * 转换构造器。这个构造器的存在是为了使某些构造在独立于维度的代码中写得更简单。例如，它允许将一个面的迭代器分配给一个线的迭代器，这个操作在2D中很有用，但在3D中没有任何意义。这里的构造函数是为了使代码符合C++的要求而存在的，但它会无条件地中止；换句话说，将一个面迭代器分配给一个线迭代器，最好放在一个if语句中，检查维度是否为2，并在3D中分配给一个四边形迭代器（如果没有这个构造函数，如果我们碰巧为2d编译，这个操作是非法的）。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 另一个对象之间的转换操作符，就像之前的那个一样，没有意义。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 复制运算符。这些运算符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的情况下使用。据推测，这里的意图是将 @p b
   * 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个存取器；因此，这个操作对三角形的迭代器没有用。
   * 因此，这个操作被声明为删除，不能使用。
   *
   */
  TriaAccessor &
  operator=(const TriaAccessor &) = delete;

  /**
   * 移动赋值运算符。允许移动。
   *
   */
  TriaAccessor &
  operator=(TriaAccessor &&) = default; // NOLINT

  /**
   * 默认的解构器。
   *
   */
  ~TriaAccessor() = default;

  /**
   * 测试元素是否被使用。 返回值是 @p true
   * 对于所有的迭代器都是正常的迭代器或活动的迭代器，只有原始迭代器可以返回
   * @p false.
   * 因为原始迭代器只在库的内部使用，你通常不需要这个函数。
   *
   */
  bool
  used() const;

  /**
   * @name  访问子对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 指向与此对象相连接的 @p ith 顶点的指针。如果
   * <code>dim=1</code> ，则抛出一个异常。
   *
   */
  TriaIterator<TriaAccessor<0, dim, spacedim>>
  vertex_iterator(const unsigned int i) const;

  /**
   * 返回当前对象的第i个顶点的全局索引。关于顶点编号的惯例在GeometryInfo类的文档中有所规定。
   * 请注意，返回值只是几何顶点的索引。
   * 它与可能与之相关的自由度无关。
   * 关于这一点，请参阅 @p DoFAccessor::vertex_dof_index 函数。
   * @note
   * 尽管有这个名字，这里返回的索引只是全局的，即它是特定于一个特定的三角形对象，或者，在三角形实际上是
   * parallel::distributed::Triangulation,
   * 类型的情况下，特定于存储在当前处理器上的分布式三角形的那一部分。
   *
   */
  unsigned int
  vertex_index(const unsigned int i) const;

  /**
   * 返回一个对 @p ith
   * 顶点的引用。该引用不是常量，也就是说，可以在赋值的左侧调用该函数，从而在三角剖分中移动单元的顶点。
   * 当然，这样做需要你确保顶点的新位置保持有用。
   *
   * - 例如，避免倒置或其他扭曲的情况（另见 @ref GlossDistorted "本词汇表条目"）。
   * @note
   * 当一个单元被细化时，它的子单元继承了它们与母单元共享的那些顶点的位置（加上为新的子单元创建的边、面和单元内部的新顶点的位置）。如果一个单元的顶点被移动，这意味着其子单元也将使用这些新的位置。
   * 另一方面，想象一下2D的情况，你有一个被精化的单元（有四个子单元），然后你移动连接所有四个子单元的中心顶点。如果你再次将这四个子单元粗化到母单元，那么移动的顶点的位置就会丢失，如果在以后的步骤中，你再次细化母单元，那么新的顶点将再次被放置在与第一次相同的位置上。
   *
   * - 即，不在你之前移动的位置。
   * @note  如果你有一个 parallel::distributed::Triangulation 对象，上述的行为是相关的。在那里，细化一个网格总是涉及到重新划分。换句话说，你可能在一个处理器上移动到不同位置的本地拥有的单元的顶点（见 @ref GlossLocallyOwnedCell "这个词汇表条目"
   * ），在网格细化时可能被移动到不同的处理器上（即使这些特定的单元没有被细化），这将根据他们之前拥有的粗略单元的位置重新创建他们的位置，而不是根据这些顶点在之前拥有它们的处理器的位置。换句话说，在并行计算中，你可能必须在每次网格细化后明确地移动节点，因为顶点的位置可能会也可能不会在伴随着网格细化的重新划分中被保留下来。
   *
   */
  Point<spacedim> &
  vertex(const unsigned int i) const;

  /**
   * 指向与此物体相邻的 @p ith 线的指针。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator
    line(const unsigned int i) const;

  /**
   * 围绕此对象的 @p ith 线的索引。
   * 只对<tt>structdim>1</tt>实现，否则会产生异常。
   *
   */
  unsigned int
  line_index(const unsigned int i) const;

  /**
   * 指向包围此对象的 @p ith 四边形的指针。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * 绑定此对象的 @p ith 四边形的索引。
   * 只对<tt>structdim>2</tt>实现，否则产生异常。
   *
   */
  unsigned int
  quad_index(const unsigned int i) const;
  /**
   * @}
   *
   */

  /**
   * @name  子对象的方向
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回索引为 @p face 的面的法线是否指向标准方向（ @p
   * true) ），或者是否相反（ @p false).
   * 哪个是标准方向，用GeometryInfo类来记录。在1d和2d中，它总是
   * @p true,
   * ，但在3d中它可能是不同的，见GeometryInfo类文档中的相应讨论。
   * 这个函数实际上只在库的内部使用，除非你绝对知道这是怎么回事。
   *
   */
  bool
  face_orientation(const unsigned int face) const;

  /**
   * 返回索引为 @p face 的面是否旋转了180度（ @p
   * 为真）或不旋转（ @p false). 在1d和2d中，这总是 @p false,
   * ，但在3d中可能不同，见GeometryInfo类文档中的相关讨论。
   * 这个函数实际上只在库的内部使用，除非你绝对知道这是怎么回事。
   *
   */
  bool
  face_flip(const unsigned int face) const;

  /**
   * 返回索引为 @p face 的面是否旋转了90度（ @p
   * 为真）或不旋转（ @p false). 在1d和2d中，这总是 @p false,
   * ，但在3d中可能不同，见GeometryInfo类文档中的相关讨论。
   * 这个函数实际上只在库的内部使用，除非你绝对知道这是怎么回事。
   *
   */
  bool
  face_rotation(const unsigned int face) const;

  /**
   * 返回索引为 @p line 的线条是否朝向标准方向。  @p true
   * 表示线的方向是从顶点0到顶点1，否则就是相反的方向。在1d和2d中，它总是
   * @p true,
   * ，但在3d中它可能是不同的，见GeometryInfo类文档中的相应讨论。
   * 这个函数实际上只在库的内部使用，除非你绝对知道这是怎么回事。
   * 这个函数查询 ReferenceCell::standard_vs_true_line_orientation().
   * 。
   *
   */
  bool
  line_orientation(const unsigned int line) const;
  /**
   * @}
   *
   */

  /**
   * @name 访问儿童
   *
   */
  /**
   * @{
   *
   */

  /**
   * 测试对象是否有子代。
   *
   */
  bool
  has_children() const;

  /**
   * 返回此对象的直系子女数。未精炼的单元格的子代数为零。
   *
   */
  unsigned int
  n_children() const;

  /**
   * @deprecated  使用n_active_descendants()代替。
   *
   */
  DEAL_II_DEPRECATED
  unsigned int
  number_of_children() const;

  /**
   * 计算并返回该对象的活动后裔的数量。例如，如果一个六边形的所有八个子单元都被进一步各向同性地精炼了一次，返回的数字将是64，而不是80。
   * 如果目前的单元没有被细化，则返回一个。
   * 如果把三角结构看成是一个森林，每棵树的根都是粗大的网格单元，节点都有后代（单元的子女），那么这个函数就会返回源自当前对象的子树中终端节点的数量；因此，如果当前对象没有被进一步细化，答案是1。
   *
   */
  unsigned int
  n_active_descendants() const;

  /**
   * 返回这个对象被精炼的次数。请注意，并不是所有的子对象都被精炼了那么多次（这就是为什么我们要在前面加上
   * @p max_),
   * 的原因，返回的数字是这个对象的任何分支中的最大精炼次数。
   * 例如，如果这个对象被精炼了，并且它的一个子对象正好再被精炼一次，那么<tt>max_refinement_depth</tt>应该返回2。
   * 如果这个对象没有被精炼（即它是活动的），那么返回值是0。
   *
   */
  unsigned int
  max_refinement_depth() const;

  /**
   * 返回一个指向 @p ith 子对象的迭代器。
   *
   */
  TriaIterator<TriaAccessor<structdim, dim, spacedim>>
  child(const unsigned int i) const;

  /**
   * 返回当前单元格上 @p child 的子数。这是
   * TriaAccessor::child(). 的反函数。
   *
   */
  unsigned int
  child_iterator_to_index(
    const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &child) const;

  /**
   * 返回该对象的迭代器，该对象与各向同性细化的第i个孩子相同。如果当前对象是各向同性细化的，那么返回的对象就是第i个子对象。如果当前对象是各向异性精炼的，那么返回的子对象实际上可能是该对象的孙子，或者根本就不存在（在这种情况下会产生一个异常）。
   *
   */
  TriaIterator<TriaAccessor<structdim, dim, spacedim>>
  isotropic_child(const unsigned int i) const;

  /**
   * 返回这个单元格的RefinementCase。
   *
   */
  RefinementCase<structdim>
  refinement_case() const;

  /**
   * @p ith
   * 子单元的索引。如果单元格的子代被访问，则子代的级别比本单元格的级别高一个。面的孩子没有级别。如果子单元不存在。
   *
   * - 被返回。
   *
   */
  int
  child_index(const unsigned int i) const;

  /**
   * @p ith
   * 各向同性的子女的索引。参见isotropic_child()函数对这个概念的定义。
   * 如果该子集不存在。
   *
   * - 被返回。
   *
   */
  int
  isotropic_child_index(const unsigned int i) const;
  /**
   * @}
   *
   */

  /**
   * @name  处理边界指标
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回此对象的边界指示器。    如果返回值是特殊值 numbers::internal_face_boundary_id, ，那么这个对象就在域的内部。      @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  types::boundary_id
  boundary_id() const;

  /**
   * 设置当前对象的边界指标。这与boundary_id()函数的情况相同。
   * 这个函数只设置当前对象本身的边界对象，而不是绑定它的指标。例如，在3D中，如果这个函数被调用到一个面，那么绑定这个面的4条边的边界指标保持不变。如果你想同时设置面和边的边界指标，请使用set_all_boundary_ids()函数。你可以在
   * step-49  的结果部分看到没有使用正确函数的结果。
   * @warning
   * 你不应该设置内部面（不在域的边界上的面）的边界指示器，也不应该将外部面的边界指示器设置为
   * numbers::internal_face_boundary_id
   * （这个值是为其他目的保留的）。如果边界单元的边界指示器为
   * numbers::internal_face_boundary_id
   * ，或者内部单元的边界指示器不是
   * numbers::internal_face_boundary_id.
   * ，算法可能无法工作或产生非常混乱的结果。
   * 不幸的是，当前对象没有办法找出它是否真的处于域的边界，因此无法确定你试图设置的值在当前情况下是否有意义。
   * @ingroup boundary   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  void
  set_boundary_id(const types::boundary_id) const;

  /**
   * 像set_boundary_id()那样做，但也要设置约束当前对象的边界指标。例如，在3D中，如果对一个面调用set_boundary_id()，那么绑定该面的4条边的边界指标保持不变。相反，如果你调用当前函数，面和边的边界指示器都被设置为给定值。
   * 如果你在3D中设置面的边界指示器（在2D中，该函数的作用与set_boundary_id()相同），并且你这样做是因为你想用一个弯曲的边界对象来表示与当前面相对应的那部分边界，那么这个函数就很有用。在这种情况下，Triangulation类需要弄清楚在网格细化时将新的顶点放在哪里，高阶Mapping对象也需要弄清楚曲线边界近似的新插值点应该在哪里。在这两种情况下，这两个类首先要确定边界面边缘的插值点，询问边界对象，然后再向边界对象询问对应于边界面内部的插值点。为了使其正常工作，仅仅设置了面的边界指示器是不够的，还需要设置约束面的边缘的边界指示器。这个函数一次就完成了所有这些工作。你可以在
   * step-49  的结果部分看到没有使用正确函数的结果。
   * @ingroup boundary   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  void
  set_all_boundary_ids(const types::boundary_id) const;

  /**
   * 返回这个对象是否在边界上。显然，这个函数的使用只适用于<tt>dim
   * @>structdim</tt>;
   * 然而，对于<tt>dim==structdim</tt>，一个对象是一个单元，CellAccessor类提供了另一种可能性来确定一个单元是否在边界。
   *
   */
  bool
  at_boundary() const;

  /**
   * 返回该对象所使用的流形对象的常数引用。
   * 正如 @ref manifold
   * 模块中所解释的，寻找合适的流形描述的过程涉及到查询流形或边界指标。更多信息见那里。
   *
   */
  const Manifold<dim, spacedim> &
  get_manifold() const;

  /**
   * @}
   *
   */

  /**
   * @name  处理流形指标的问题
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回此对象的歧管指标。    如果返回值是特殊值 numbers::flat_manifold_id, ，那么这个对象就与一个标准的笛卡尔流形描述相关。      @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  types::manifold_id
  manifold_id() const;

  /**
   * 设置流形指示器。 与<tt>manifold_id()</tt>函数同样适用。
   * 注意，它只设置当前对象本身的流形对象，而不是绑定它的指标，也不是其子对象的指标。例如，在3D中，如果对一个面调用这个函数，那么绑定该面的4条边的流形指标就不会改变。如果你想同时设置面、边和所有子节点的流形指标，请使用set_all_manifold_ids()函数。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_manifold_id(const types::manifold_id) const;

  /**
   * 像set_manifold_id()那样做，但也要设置绑定当前对象的流形指标。例如，在3D中，如果对一个面调用了set_manifold_id()，那么绑定该面的4条边的流形指标就不会改变。另一方面，面和边的流形指标都是用当前函数同时设置的。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_all_manifold_ids(const types::manifold_id) const;

  /**
   * @}
   *
   */


  /**
   * @name  用户数据
   *
   */
  /**
   * @{
   *
   */
  /**
   * 读取用户标志。更多信息见 @ref
   * GlossUserFlags 。
   *
   */
  bool
  user_flag_set() const;

  /**
   * 设置用户标志。参见 @ref
   * GlossUserFlags 了解更多信息。
   *
   */
  void
  set_user_flag() const;

  /**
   * 清除用户标志。参见
   * @ref GlossUserFlags  了解更多信息。
   *
   */
  void
  clear_user_flag() const;

  /**
   * 为这个和所有的子代设置用户标志。参见 @ref
   * GlossUserFlags 以了解更多信息。
   *
   */
  void
  recursively_set_user_flag() const;

  /**
   * 清除此项目和所有子项目的用户标志。参见 @ref
   * GlossUserFlags 以了解更多信息。
   *
   */
  void
  recursively_clear_user_flag() const;

  /**
   * 将用户数据重置为零，与指针或索引无关。更多信息见 @ref
   * GlossUserData 。
   *
   */
  void
  clear_user_data() const;

  /**
   * 将用户指针设置为 @p p. 。
   * @note
   * 用户指针和用户索引是相互排斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data() 。    参见 @ref GlossUserData
   * 以了解更多信息。
   *
   */
  void
  set_user_pointer(void *p) const;

  /**
   * 将用户指针重置为 @p
   * nullptr 指针。参见 @ref GlossUserData 以了解更多信息。
   *
   */
  void
  clear_user_pointer() const;

  /**
   * 访问用户指针的值。用户有责任确保该指针指向有用的东西。你应该使用新式的cast操作符来保持最低限度的类型安全，例如。
   * @note
   * 用户指针和用户索引是互斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data() 。<tt>A
   * a=static_cast<A*>(cell->user_pointer())；</tt>。    更多信息见
   * @ref GlossUserData 。
   *
   */
  void *
  user_pointer() const;

  /**
   * 将此对象及其所有子对象的用户指针设置为给定值。例如，如果某个子域的所有单元格，或者边界的某个部分的所有面都应该有用户指针指向描述这部分域或边界的对象，这就很有用。
   * 请注意，用户指针在网格细化过程中是不被继承的，所以在网格细化之后，可能会有单元格或面没有用户指针指向描述对象。在这种情况下，只需简单地循环所有具有此信息的最粗层次的元素，并使用此函数递归地设置三角形的所有更细层次的用户指针。
   * @note
   * 用户指针和用户索引是互斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data()  。    更多信息见 @ref
   * GlossUserData 。
   *
   */
  void
  recursively_set_user_pointer(void *p) const;

  /**
   * 清除这个对象的用户指针和它所有的子对象。这与recursively_set_user_pointer()函数的说法相同。更多信息见 @ref
   * GlossUserData 。
   *
   */
  void
  recursively_clear_user_pointer() const;

  /**
   * 将用户索引设置为 @p p. 。
   * @note
   * 用户指针和用户索引是相互排斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data() 。更多信息见 @ref GlossUserData
   * 。
   *
   */
  void
  set_user_index(const unsigned int p) const;

  /**
   * 将用户索引重置为0。更多信息见 @ref
   * GlossUserData 。
   *
   */
  void
  clear_user_index() const;

  /**
   * 访问用户索引的值。
   * @note
   * 用户指针和用户索引是相互排斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data() 。    更多信息见 @ref
   * GlossUserData 。
   *
   */
  unsigned int
  user_index() const;

  /**
   * 设置此对象及其所有子对象的用户索引。
   * 请注意，用户索引在网格细化中是不被继承的，所以在网格细化后，可能会有单元格或面的用户索引不符合预期。在这种情况下，只需循环查看所有拥有该信息的最粗层次的元素，并使用该函数递归设置三角结构中所有更细层次的用户索引。
   * @note
   * 用户指针和用户索引是相互排斥的。因此，你只能使用其中之一，除非你在中间调用
   * Triangulation::clear_user_data() 。    更多信息见 @ref
   * GlossUserData 。
   *
   */
  void
  recursively_set_user_index(const unsigned int p) const;

  /**
   * 清除这个对象的用户索引和它的所有后代。这与recursively_set_user_index()函数的说法相同。
   * 更多信息见 @ref GlossUserData 。
   *
   */
  void
  recursively_clear_user_index() const;
  /**
   * @}
   *
   */

  /**
   * @name  关于一个物体的几何信息
   *
   */
  /**
   *
   */

  /**
   * 物体的直径。
   * 一个物体的直径被计算为当前物体的最大对角线。如果这个物体是四边形，那么有两条这样的对角线，如果是六面体，那么有四条对角线连接
   * "对面
   * "的点。对于三角形和四面体，该函数只是返回最长边的长度。
   * 楔形和金字塔的情况就比较困难了。对于楔形，我们返回三个四边形面的最长对角线的长度或两个三角形面的最长边的长度。对于金字塔，同样的原则也适用。
   * 在所有这些情况下，这个 "直径
   * "的定义不一定是物体内部各点之间最大距离意义上的真正直径。事实上，我们经常可以构造出不是这样的物体，尽管这些物体与参考形状相比一般都有很大的变形。此外，对于可能使用高阶映射的物体，我们可能会有凸起的面，这也会给计算物体直径的精确表示带来麻烦。也就是说，上面使用的定义对于大多数计算来说是完全足够的。
   *
   */
  double
  diameter() const;

  /**
   * 返回一对Point和double，对应于物体的中心和合理的小包围球的半径。
   * 该函数实现了Ritter的O(n)算法，以获得一个围绕对象顶点的合理的小包围球。
   * 包围球的初始猜测是包含对象的最大对角线作为其直径的球。
   * 从这样的初始猜测开始，算法测试对象的所有顶点（除了最大对角线的顶点）在几何上是否在球内。
   * 如果发现任何顶点（v）在几何上不在球内，则通过移动球的中心和增加半径来构建一个新的球的迭代，以便在几何上包围先前的球和顶点（v）。当所有的顶点都在球的几何范围内时，该算法就结束了。
   * 如果一个顶点（v）在几何上位于球的某个迭代中，那么它将在球的后续迭代中继续如此（这是真实的，通过/a结构）。
   * @note  这个函数假定从参考单元开始的d-线性映射。    <a
   * href="http://geomalgorithms.com/a08-_containers.html">see
   * this</a>和[Ritter 1990]。
   *
   */
  std::pair<Point<spacedim>, double>
  enclosing_ball() const;

  /**
   * 返回包围该对象的最小边界框。
   * 注意，这个方法不知道你可能用来做计算的任何映射。如果你使用的是修改顶点位置的映射对象，比如MappingQEulerian，或者MappingFEField，那么你应该调用函数
   * Mapping::get_bounding_box() 代替。
   *
   */
  BoundingBox<spacedim>
  bounding_box() const;

  /**
   * 对象在给定轴线方向的长度，在本地坐标系中指定。关于本地轴的含义和列举，请参见GeometryInfo的文档。
   * 请注意，一个物体的 "长度
   * "可以有多种解释。在这里，我们选择它作为对象的任何一条边的最大长度，这些边与参考单元上选择的轴平行。
   *
   */
  double
  extent_in_direction(const unsigned int axis) const;

  /**
   * 返回任何两个顶点之间的最小距离。
   *
   */
  double
  minimum_vertex_distance() const;

  /**
   * 返回一个属于此对象所在的流形<dim,spacedim>的点，给定其在参考单元格
   * @p structdim
   * 上的参数化坐标。该函数可查询底层流形对象，并可用于获得该对象上任意点的精确几何位置。
   * 注意，参数 @p coordinates 是 <em> 参考单元 </em>
   * 上的坐标，以参考坐标给出。换句话说，该参数提供了不同顶点之间的权重。例如，对于线，调用这个函数的参数是Point<1>(.5)，相当于要求线的中心。
   *
   */
  Point<spacedim>
  intermediate_point(const Point<structdim> &coordinates) const;

  /**
   * 这个函数通过从参考 $d$ 维单元反转 $d$
   * 线性函数的仿射近似，计算从实数到单元格的快速近似转换。
   * 单位单元到实数单元映射的仿生近似是通过仿生函数与本对象的
   * $2^d$
   * 顶点的最小二乘法拟合找到的。对于任何有效的网格单元，其几何形状不是退化的，这个操作的结果是一个唯一的仿生映射。因此，对于所有给定的输入点，这个函数将返回一个有限的结果，即使在实际的双/三线或高阶映射的实际转换可能是单数的情况下。除了仅从顶点近似映射外，该函数还忽略了附加的流形描述。只有在从单元格到实际单元格的转换确实是仿生的情况下，结果才是准确的，比如在一维或二维/三维的笛卡尔和仿生（平行四边形）网格。
   * 对于单元格的精确变换，使用
   * Mapping::transform_real_to_unit_cell(). 。
   * @note  如果dim<spacedim，我们首先将p投影到平面上。
   *
   */
  Point<structdim>
  real_to_unit_cell_affine_approximation(const Point<spacedim> &point) const;

  /**
   * 物体的中心。对象的中心被定义为顶点位置的平均值，这也是
   * $Q_1$
   * 映射将映射参考单元的中心的地方。然而，你也可以要求这个函数返回与当前对象相关的底层流形对象所计算的顶点的平均值，方法是将可选参数
   * @p respect_manifold.  设置为 true
   * 流形通常会将顶点的坐标拉回到参考域（不一定是参考单元），在那里计算平均值，然后再将平均值点的坐标推到物理空间；结果点保证位于流形内，即使流形是弯曲的。
   * 当对象使用不同的流形描述作为其周围环境时，比如该TriaAccessor的部分边界对象使用非平面流形描述，但对象本身是平面的，
   * TriaAccessor::center()
   * 函数给出的结果可能不够准确，即使参数 @p
   * respect_manifold
   * 被设置为真。如果你发现这种情况，你可以通过将第二个附加参数
   * @p
   * interpolate_from_surrounding设置为真来进一步完善中心的计算。这将通过从所有边界对象的中心进行所谓的无限插值来计算中心的位置。对于一个二维物体，它在四条周围线中的每一条上加了
   * <code>1/2</code> 的权重，在四个顶点上加了 <code>-1/4</code>
   * 的权重。这相当于四个面的描述之间的线性插值，减去顶点的贡献，当通过与顶点相邻的两条线时，顶点的贡献被加了两次。在三维中，面的权重是
   * <code>1/2</code> ，线的权重是 <code>-1/4</code>
   * ，而顶点的权重是 <code>1/8</code>
   * 。为了进一步了解情况，还可以赋予TransfiniteInterpolationManifold类，该类不仅能够将这种有益的描述应用于单个单元，而且能够应用于粗略单元的所有子女。
   *
   */
  Point<spacedim>
  center(const bool respect_manifold             = false,
         const bool interpolate_from_surrounding = false) const;

  /**
   * 返回对象的arycenter（也叫中心点）。在 $D$ 空间维度中，维度为 $K$ 的对象的arycenter由@f[
   * \mathbf x_K = \frac{1}{|K|} \int_K \mathbf x \; \textrm{d}x
   * @f]定义的 $\mathbf x_K$ -维向量给出，其中对象的度量由@f[
   * |K| = \int_K \mathbf 1 \; \textrm{d}x. @f]给出。 这个函数假定
   * $K$ 由 $d$ -线性函数从参考 $d$
   * -维单元映射出来。然后，上面的积分可以被拉回到参考单元，并准确地进行评估（如果通过冗长的，并且与Center()函数相比，昂贵的计算）。
   *
   */
  Point<spacedim>
  barycenter() const;

  /**
   * 计算对象的dim-dimensional度量。对于二维空间中的二维单元，这等于其体积。另一方面，对于三维空间中的二维单元，或者如果当前指向的对象是三维空间中的三维单元的二维面，那么该函数就会计算该对象所占的面积。对于一个一维的物体，返回其长度。
   * 该函数只计算假定由（双/三）线性映射表示的单元、面或边的量度。换句话说，它只考虑绑定当前对象的顶点的位置，而不考虑对象的内部实际上如何被映射。在大多数简单的情况下，这正是你想要的。然而，对于不
   * "直
   * "的对象，例如，嵌入三维空间的二维单元作为弯曲域的三角化的一部分，三维单元的二维面不仅仅是平行四边形，或者位于域的边界的面不仅仅是由直线段或平面限定的，这个函数只计算由描述有关对象的真实几何的流形或边界对象定义的
   * "真实
   * "对象的（双/三）线性内插的二维度量。如果你想考虑
   * "真实
   * "几何，你将需要通过在物体上积分一个等于1的函数来计算这个度量，在应用正交之后，等于将你想用于积分的FEValues或FEFaceValues对象返回的JxW值相加。
   *
   */
  double
  measure() const;

  /**
   * 如果当前对象是给定参数的翻译，则返回真。
   * @note
   * 为了三角测量的目的，单元格、面等只由其顶点来描述。因此，当前函数只对顶点的位置进行比较。然而，对于许多实际应用来说，决定一个单元是否是另一个单元的翻译的不仅仅是顶点，还有单元如何从参考单元映射到它在现实空间中的位置。例如，如果我们使用高阶映射，那么不仅顶点必须是彼此的平移，而且还必须是沿边缘的点。因此，在这些问题中，应该问映射，而不是当前的函数，两个对象是否是彼此的平移。
   *
   */
  bool
  is_translation_of(
    const TriaIterator<TriaAccessor<structdim, dim, spacedim>> &o) const;

  /**
   * 当前对象的参考单元格类型。
   *
   */
  ReferenceCell
  reference_cell() const;

  /**
   * 顶点的数量。
   *
   */
  unsigned int
  n_vertices() const;

  /**
   * 行的数量。
   *
   */
  unsigned int
  n_lines() const;

  /**
   * 面的数量。
   * @note  只对单元格实现（dim==spacedim）。
   *
   */
  unsigned int
  n_faces() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_vertices()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_lines()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_faces()的所有索引的数组。
   * @note 只对单元格实现（dim==spacedim）。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  face_indices() const;

  /**
   * @}
   *
   */

private:
  /**
   * 与set_boundary_id类似，但没有检查内部面或无效的id。
   *
   */
  void
  set_boundary_id_internal(const types::boundary_id id) const;

  /**
   * 设置那些约束当前对象的对象的索引。例如，如果当前对象代表一个单元格，那么该参数表示绑定该单元格的面的索引。如果当前对象代表一条线，那么该参数表示绑定该线的顶点的指数。以此类推。
   *
   */
  void
  set_bounding_object_indices(
    const std::initializer_list<int> &new_indices) const;

  /**
   * 与上述相同，但对`无符号int`而言。
   *
   */
  void
  set_bounding_object_indices(
    const std::initializer_list<unsigned int> &new_indices) const;

  /**
   * 设置标志，表明 <code>line_orientation()</code> 将返回什么。
   * 只有在3D中才能设置面的线方向（即 <code>structdim==2 &&
   * dim==3</code>  ）。
   *
   */
  void
  set_line_orientation(const unsigned int line, const bool orientation) const;

  /**
   * 设置索引为 @p face 的四边形的法线是否指向标准方向（
   * @p true) ），或者是否相反（ @p false).
   * 哪个是标准方向，用GeometryInfo类记录。
   * 这个函数只用于库的内部使用。将这个标志设置为任何其他的值，而不是三角测量已经设置的值，必然会给你带来灾难。
   *
   */
  void
  set_face_orientation(const unsigned int face, const bool orientation) const;

  /**
   * 设置标志表明， <code>face_flip()</code> 将返回什么。
   * 只有在三维中才能设置单元格的面朝向（即
   * <code>structdim==3 && dim==3</code>  ）。
   *
   */
  void
  set_face_flip(const unsigned int face, const bool flip) const;

  /**
   * 设置标志，表明 <code>face_rotation()</code> 将返回什么。
   * 只有在3D中才能设置单元格的面朝向（即 <code>structdim==3
   * && dim==3</code>  ）。
   *
   */
  void
  set_face_rotation(const unsigned int face, const bool rotation) const;

  /**
   * 设置 @p used 标志。仅用于库的内部使用。
   *
   */
  void
  set_used_flag() const;

  /**
   * 清除 @p used 标志。仅供图书馆内部使用。
   *
   */
  void
  clear_used_flag() const;

  /**
   * 设置 @p RefinementCase<dim>
   * 这个TriaObject的精炼。对于<tt>structdim=1</tt>不定义，因为线条总是被细化成2条子线（各向同性细化）。
   * 如果你接触这个函数，你应该很清楚你在做什么。它是专门用于库的内部使用的。
   *
   */
  void
  set_refinement_case(const RefinementCase<structdim> &ref_case) const;

  /**
   * 清除这个TriaObject的RefinementCase<dim>，即重置为
   * RefinementCase<dim>::no_refinement.
   * 你应该很清楚你在做什么，如果你碰触这个函数。它是专门用于库的内部使用的。
   *
   */
  void
  clear_refinement_case() const;

  /**
   * 设置第1个孩子的索引。由于孩子们至少是成对出现的，我们只需要存储第二个孩子的索引，也就是偶数孩子的索引。请确保首先设置第i=0个孩子的索引。不允许对奇数的孩子调用这个函数。
   *
   */
  void
  set_children(const unsigned int i, const int index) const;

  /**
   * 清除子字段，即把它设置为一个值，表示这个单元格没有子字。
   *
   */
  void
  clear_children() const;

private:
  template <int, int>
  friend class Triangulation;

  friend struct dealii::internal::TriangulationImplementation::Implementation;
  friend struct dealii::internal::TriangulationImplementation::
    ImplementationMixedMesh;
  friend struct dealii::internal::TriaAccessorImplementation::Implementation;
};



/**
 * 这个类是<code>TriaAccessor<structdim, dim,
 * spacedim></code>的特殊化，用于 @p structdim
 * 为零的情况。该类代表维数为 <code>dim</code>
 * 的三角结构中的顶点（即1代表线的三角结构，2代表四边形的三角结构，3代表六边形的三角结构），该三角结构嵌入在维数为
 * <code>spacedim</code> 的空间中（对于 <code>spacedim==dim</code>
 * 的三角结构代表 ${\mathbb R}^\text{dim}$ 的域，对于
 * <code>spacedim@>dim</code>
 * 的三角结构是嵌入在高维空间的流形）。 对于 @p dim
 * 等于1的情况，即对于一维三角结构的顶点，这个类别还有一个进一步的特殊化，因为在这种情况下顶点也是面。
 *
 *
 * @ingroup Accessors
 *
 *
 */
template <int dim, int spacedim>
class TriaAccessor<0, dim, spacedim>
{
public:
  /**
   * 这个访问器所代表的对象所处的空间的维度。
   * 例如，如果这个访问器代表一个四边形，它是四维空间中一个二维表面的一部分，那么这个值就是四。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 这个访问器所代表的事物的尺寸，是这个访问器的一部分。例如，如果这个访问器代表一条线，是六面体的一部分，那么这个值就是三。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 这个访问器所代表的当前对象的尺寸。例如，如果它是线（不管它是四面体还是六面体的一部分，也不管我们处于什么维度），那么这个值就等于1。
   *
   */
  static const unsigned int structure_dimension = 0;

  /**
   * 指向内部数据的指针。
   *
   */
  using AccessorData = void;

  /**
   * 构造函数。第二个参数是我们指向的顶点的全局索引。
   *
   */
  TriaAccessor(const Triangulation<dim, spacedim> *tria,
               const unsigned int                  vertex_index);

  /**
   * 构造函数。这个构造函数的存在是为了保持与其他访问器类的接口兼容性。
   * @p index 可以用来设置我们所指向的顶点的全局索引。
   *
   */
  TriaAccessor(const Triangulation<dim, spacedim> *tria  = nullptr,
               const int                           level = 0,
               const int                           index = 0,
               const AccessorData *                      = nullptr);

  /**
   * 构造函数。不应该被调用，因此会产生一个错误。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 构造函数。不应该被调用，因此会产生一个错误。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 返回迭代器的状态。
   *
   */
  IteratorState::IteratorStates
  state() const;

  /**
   * 这个对象的水平。顶点没有级别，所以这个函数总是返回0。
   *
   */
  static int
  level();

  /**
   * 此对象的索引。返回这个对象所指向的顶点的全局索引。
   *
   */
  int
  index() const;

  /**
   * 返回这个类所指向的对象所属于的三角形的引用。
   *
   */
  const Triangulation<dim, spacedim> &
  get_triangulation() const;

  /**
   * @name  迭代器的进阶
   *
   */
  /**
   * @{
   *
   */
  /**
   * 该操作符将迭代器推进到下一个元素。
   *
   */
  void
  operator++();

  /**
   * 该操作符将迭代器移到上一个元素。
   *
   */
  void
  operator--();
  /**
   * 进行平等比较。
   *
   */
  bool
  operator==(const TriaAccessor &) const;

  /**
   * 不等式的比较。
   *
   */
  bool
  operator!=(const TriaAccessor &) const;

  /**
   * @}
   *
   */


  /**
   * @name  访问子对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回当前对象的第i个顶点的全局索引。如果 @p i
   * 为零，这将返回这个对象所指向的当前点的索引。否则，它会抛出一个异常。
   * 请注意，返回值只是几何顶点的索引。
   * 它与与之相关的可能的自由度没有关系。
   * 关于这一点，请参阅 @p DoFAccessor::vertex_dof_index 函数。
   * @note
   * 尽管有这个名字，这里返回的索引只是全局的，即它是特定于一个特定的三角形对象，或者，在三角形实际上是
   * parallel::distributed::Triangulation,
   * 类型的情况下，特定于存储在当前处理器上的分布式三角形的那一部分。
   *
   */
  unsigned int
  vertex_index(const unsigned int i = 0) const;

  /**
   * 返回一个对 @p ith
   * 顶点的引用。如果i是0，这将返回这个对象所指向的当前点的引用。否则，它会抛出一个异常。
   *
   */
  Point<spacedim> &
  vertex(const unsigned int i = 0) const;

  /**
   * 指向与此对象相邻的 @p ith
   * 线的指针。将指向一个无效的对象。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::line_iterator static line(const unsigned int);

  /**
   * 围绕此对象的 @p ith 行的行索引。抛出一个异常。
   *
   */
  static unsigned int
  line_index(const unsigned int i);

  /**
   * 指向包围此对象的 @p ith 四边形的指针。
   *
   */
  static typename dealii::internal::TriangulationImplementation::
    Iterators<dim, spacedim>::quad_iterator
    quad(const unsigned int i);

  /**
   * 绑定此对象的 @p ith 四边形的索引。抛出一个异常。
   *
   */
  static unsigned int
  quad_index(const unsigned int i);

  /**
   * @}
   *
   */


  /**
   * @name  关于一个物体的几何信息
   *
   */
  /**
   * @{
   *
   */

  /**
   * 物体的直径。这个函数总是返回零。
   *
   */
  double
  diameter() const;

  /**
   * 对象在给定轴线方向的长度，在本地坐标系中指定。有关本地轴的含义和列举，请参见GeometryInfo的文档。
   * 这个函数总是返回0。
   *
   */
  double
  extent_in_direction(const unsigned int axis) const;

  /**
   * 返回此对象的中心，当然，这与此对象所指的顶点的位置相吻合。参数
   * @p  respect_manifold和 @p interpolate_from_surrounding
   * 不被使用。它们的存在是为了提供与
   * <code>TriaAccessor<structdim,dim,spacedim></code>  相同的接口。
   *
   */
  Point<spacedim>
  center(const bool respect_manifold             = false,
         const bool interpolate_from_surrounding = false) const;

  /**
   * 计算对象的dim-dimensional度量。对于二维空间中的二维单元，这等于其体积。另一方面，对于三维空间中的二维单元，或者如果当前指向的对象是三维空间中的三维单元的二维面，那么该函数就会计算该对象所占的面积。对于一个一维的对象，返回其长度。对于一个零维的对象，返回零。
   *
   */
  double
  measure() const;
  /**
   * @}
   *
   */

  /**
   * @name  子对象的方向
   *
   */
  /**
   * @{
   *
   */

  /**
   * @brief Always return false
   *
   */
  static bool
  face_orientation(const unsigned int face);

  /**
   * @brief Always return false
   *
   */
  static bool
  face_flip(const unsigned int face);

  /**
   * @brief Always return false
   *
   */
  static bool
  face_rotation(const unsigned int face);

  /**
   * @brief Always return false
   *
   */
  static bool
  line_orientation(const unsigned int line);

  /**
   * @}
   *
   */

  /**
   * @name  访问儿童
   *
   */
  /**
   * @{
   *
   */

  /**
   * 测试该对象是否有子代。总是假的。
   *
   */
  static bool
  has_children();

  /**
   * 返回这个对象的直接子代数。这始终是零。
   *
   */
  static unsigned int
  n_children();

  /**
   * 计算并返回此对象的活动子孙的数量。  总是零。
   *
   */
  static unsigned int
  n_active_descendants();

  /**
   * @deprecated  使用n_active_descendants()代替。
   *
   */
  DEAL_II_DEPRECATED
  static unsigned int
  number_of_children();


  /**
   * 返回此对象被精炼的次数。总是0。
   *
   */
  static unsigned int
  max_refinement_depth();

  /**
   * @brief Return an invalid unsigned integer.
   *
   */
  static unsigned int
  child_iterator_to_index(const TriaIterator<TriaAccessor<0, dim, spacedim>> &);

  /**
   * @brief Return an invalid object.
   *
   */
  static TriaIterator<TriaAccessor<0, dim, spacedim>>
  child(const unsigned int);

  /**
   * @brief Return an invalid object.
   *
   */
  static TriaIterator<TriaAccessor<0, dim, spacedim>>
  isotropic_child(const unsigned int);

  /**
   * 总是不返回细化。
   *
   */
  static RefinementCase<0>
  refinement_case();

  /**
   * @brief Returns  -
   *
   */
  static int
  child_index(const unsigned int i);

  /**
   * @brief Returns  -
   *
   */
  static int
  isotropic_child_index(const unsigned int i);
  /**
   * @}
   *
   */

  /**
   * 返回这里指向的顶点是否被使用。
   *
   */
  bool
  used() const;

protected:
  /**
   * 复制操作符。因为这只是从迭代器中调用，所以不要返回任何东西，因为迭代器会返回自己。
   * 这个方法是受保护的，因为它只能从迭代器类中调用。
   *
   */
  void
  copy_from(const TriaAccessor &);

  /**
   * 访问器的比较运算符。这个操作符在比较迭代器进入三角形的对象时使用，例如在把它们放入
   * `std::map`.
   * 这个操作符简单地比较了当前对象所指向的顶点的全局索引。
   *
   */
  bool
  operator<(const TriaAccessor &other) const;

  /**
   * 指向我们所操作的三角形的指针。
   *
   */
  const Triangulation<dim, spacedim> *tria;

  /**
   * 这个对象所对应的顶点的全局索引。
   *
   */
  unsigned int global_vertex_index;

private:
  template <typename Accessor>
  friend class TriaRawIterator;
  template <typename Accessor>
  friend class TriaIterator;
  template <typename Accessor>
  friend class TriaActiveIterator;
};



/**
 * 这个类是<code>TriaAccessor<structdim, dim,
 * spacedim></code>的特化，用于 @p structdim 为零和 @p dim
 * 为一的情况。该类表示嵌入在维数为 <code>spacedim</code>
 * 的空间中的一维三角形的顶点（对于
 * <code>spacedim==dim==1</code> ，三角形表示 ${\mathbb R}^\text{dim}$
 * 中的一个域，对于 <code>spacedim@>dim==1</code>
 * ，三角形是嵌入在高维空间的流形）。
 * 目前TriaAccessor<0,dim,spacedim>类对一维三角形的顶点的专门化存在，因为在
 * @p dim  == 1的情况下，顶点也是面。
 *
 *
 * @ingroup Accessors
 *
 *
 */
template <int spacedim>
class TriaAccessor<0, 1, spacedim>
{
public:
  /**
   * 这个访问器所代表的对象所处空间的尺寸。
   * 例如，如果这个访问器代表一个四边形，它是四维空间中一个二维表面的一部分，那么这个值就是四。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 这个访问器所代表的事物的尺寸，是这个访问器的一部分。例如，如果这个访问器代表一条线，是六面体的一部分，那么这个值就是三。
   *
   */
  static const unsigned int dimension = 1;

  /**
   * 这个访问器所代表的当前对象的尺寸。例如，如果它是线（不管它是四面体还是六面体的一部分，也不管我们处于什么维度），那么这个值就等于1。
   *
   */
  static const unsigned int structure_dimension = 0;

  /**
   * 指向内部数据的指针。
   *
   */
  using AccessorData = void;

  /**
   * 这里代表的顶点是在域的左端、右端，还是在内部。
   *
   */
  enum VertexKind
  {
    /**
     * 左边的顶点。
     *
     */
    left_vertex,
    /**
     * 内部顶点。
     *
     */
    interior_vertex,
    /**
     * 右边的顶点。
     *
     */
    right_vertex
  };

  /**
   * 构造函数。
   * 由于没有从顶点到单元的映射，一个点的访问器对象没有办法弄清它是否在域的边界上。因此，第二个参数必须由生成这个访问器的对象传递。
   *
   * 例如，一个1d单元可以计算出它的左或右顶点是否在边界上。
   * 第三个参数是我们指向的顶点的全局索引。
   *
   */
  TriaAccessor(const Triangulation<1, spacedim> *tria,
               const VertexKind                  vertex_kind,
               const unsigned int                vertex_index);

  /**
   * 构造函数。这个构造函数的存在是为了保持与其他访问器类的界面兼容。然而，它在这里并没有做任何有用的事情，所以实际上可能不会被调用。
   *
   */
  TriaAccessor(const Triangulation<1, spacedim> *tria = nullptr,
               const int                              = 0,
               const int                              = 0,
               const AccessorData *                   = nullptr);

  /**
   * 构造函数。不应该被调用，因此会产生一个错误。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 构造函数。不应被调用，因此产生错误。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  TriaAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 复制操作符。因为这只是从迭代器中调用，所以不要返回任何东西，因为迭代器会返回自己。
   *
   */
  void
  copy_from(const TriaAccessor &);

  /**
   * 返回迭代器的状态。由于对点的迭代器不能被增加或减少，它的状态保持不变，特别是等于
   * IteratorState::valid.  。
   *
   */
  static IteratorState::IteratorStates
  state();

  /**
   * 这个对象的级别。顶点没有级别，所以这个函数总是返回0。
   *
   */
  static int
  level();

  /**
   * 此对象的索引。返回这个对象所指向的顶点的全局索引。
   *
   */
  int
  index() const;

  /**
   * 返回这个类所指向的对象所属于的三角形的引用。
   *
   */
  const Triangulation<1, spacedim> &
  get_triangulation() const;

  /**
   * @name  推进迭代器的工作
   *
   */
  /**
   * @{
   *
   */
  /**
   * 这个操作将迭代器推进到下一个元素。对于点，这个操作没有定义，所以你不能在点的迭代器上迭代。
   *
   */
  void
  operator++() const;

  /**
   * 该操作符将迭代器移到上一个元素。对于点来说，这个操作没有被定义，所以你不能对点的迭代器进行迭代。
   *
   */
  void
  operator--() const;
  /**
   * 进行平等比较。
   *
   */
  bool
  operator==(const TriaAccessor &) const;

  /**
   * 比较不等式。
   *
   */
  bool
  operator!=(const TriaAccessor &) const;

  /**
   * 访问器的比较运算器。当比较迭代器进入三角形的对象时，例如将其放入
   * `std::map`.
   * 中时，该操作符会被使用。该操作符只是比较当前对象所指向的顶点的全局索引。
   *
   */
  bool
  operator<(const TriaAccessor &other) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问子对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回当前对象的第i个顶点的全局索引。如果i为零，则返回该对象所指向的当前点的索引。否则，它会抛出一个异常。
   * 请注意，返回值只是几何顶点的索引。
   * 它与与之相关的可能的自由度没有关系。
   * 关于这一点，请参阅 @p DoFAccessor::vertex_dof_index 函数。
   * @note
   * 尽管有这个名字，这里返回的索引只是全局的，即它是特定于一个特定的三角形对象的，或者，在三角形实际上是
   * parallel::distributed::Triangulation,
   * 类型的情况下，特定于存储在当前处理器上的分布式三角形的那一部分。
   *
   */
  unsigned int
  vertex_index(const unsigned int i = 0) const;

  /**
   * 返回一个对 @p ith
   * 顶点的引用。如果i为零，这将返回该对象所指向的当前点的引用。否则，它会抛出一个异常。
   *
   */
  Point<spacedim> &
  vertex(const unsigned int i = 0) const;

  /**
   * 返回这个对象的中心，当然，这与这个对象所指的顶点的位置是一致的。
   *
   */
  Point<spacedim>
  center() const;

  /**
   * 指向限定此对象的 @p ith
   * 线的指针。将指向一个无效的对象。
   *
   */
  typename dealii::internal::TriangulationImplementation::
    Iterators<1, spacedim>::line_iterator static line(const unsigned int);

  /**
   * 围绕此对象的 @p ith 行的行索引。
   * 只对<tt>structdim>1</tt>实现，否则产生异常。
   *
   */
  static unsigned int
  line_index(const unsigned int i);

  /**
   * 指向包围此对象的 @p ith 四边形的指针。
   *
   */
  static typename dealii::internal::TriangulationImplementation::
    Iterators<1, spacedim>::quad_iterator
    quad(const unsigned int i);

  /**
   * 绑定此对象的 @p ith 四边形的索引。
   * 只对<tt>structdim>2</tt>实现，否则会产生异常。
   *
   */
  static unsigned int
  quad_index(const unsigned int i);

  /**
   * @}
   *
   */


  /**
   * 返回这个点是否在我们这里处理的一维三角的边界上。
   *
   */
  bool
  at_boundary() const;

  /**
   * 返回这个对象的边界指标。一维三角形的惯例是，左端顶点（可构建三角形的每条线段）的边界指示器为0，右端顶点的边界指示器为1，除非明确设置为不同的值。    如果返回值是特殊值 numbers::internal_face_boundary_id, ，那么这个对象是在域的内部。      @see   @ref GlossBoundaryIndicator  "关于边界指示器的词汇条目"
   *
   */
  types::boundary_id
  boundary_id() const;

  /**
   * 返回对该对象所使用的流形对象的常数引用。
   *
   */
  const Manifold<1, spacedim> &
  get_manifold() const;

  /**
   * 返回此对象的流形指标。      @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  types::manifold_id
  manifold_id() const;


  /**
   * @name  子对象的方向
   *
   */
  /**
   * @{
   *
   */

  /**
   * @brief Always return false
   *
   */
  static bool
  face_orientation(const unsigned int face);

  /**
   *
   */
  static bool
  face_flip(const unsigned int face);

  /**
   *
   */
  static bool
  face_rotation(const unsigned int face);

  /**
   * @brief Always return false
   *
   */
  static bool
  line_orientation(const unsigned int line);

  /**
   * @}
   *
   */

  /**
   * @name  访问儿童
   *
   */
  /**
   * @{
   *
   */

  /**
   * 测试该对象是否有子代。总是假的。
   *
   */
  static bool
  has_children();

  /**
   * 返回此对象的直系子代数。在维度0中，这始终是零。
   *
   */
  static unsigned int
  n_children();

  /**
   * 计算并返回此对象的活动子孙的数量。  总是零。
   *
   */
  static unsigned int
  n_active_descendants();

  /**
   * @deprecated  使用n_active_descendants()代替。
   *
   */
  DEAL_II_DEPRECATED
  static unsigned int
  number_of_children();


  /**
   * 返回此对象被精炼的次数。总是0。
   *
   */
  static unsigned int
  max_refinement_depth();

  /**
   * @brief Return an invalid unsigned integer.
   *
   */
  static unsigned int
  child_iterator_to_index(const TriaIterator<TriaAccessor<0, 1, spacedim>> &);

  /**
   * @brief Return an invalid object
   *
   */
  static TriaIterator<TriaAccessor<0, 1, spacedim>>
  child(const unsigned int);

  /**
   * @brief Return an invalid object
   *
   */
  static TriaIterator<TriaAccessor<0, 1, spacedim>>
  isotropic_child(const unsigned int);

  /**
   * 总是不返回细化。
   *
   */
  static RefinementCase<0>
  refinement_case();

  /**
   * @brief Returns
   *
   * -
   *
   */
  static int
  child_index(const unsigned int i);

  /**
   * @brief Returns  -
   *
   */
  static int
  isotropic_child_index(const unsigned int i);
  /**
   * @}
   *
   */

  /**
   * @name  处理边界指标
   *
   */
  /**
   * @{
   *
   */

  /**
   * 设置边界指示器。这与<tt>boundary_id()</tt>函数同样适用。
   * @warning
   * 你不应该设置内部面（不在域的边界上的面）的边界指示器，也不应该将外部面的边界指示器设置为
   * numbers::internal_face_boundary_id
   * （这个值是为其他目的保留的）。如果边界单元的边界指示器为
   * numbers::internal_face_boundary_id
   * ，或者内部单元的边界指示器不是
   * numbers::internal_face_boundary_id.
   * ，算法可能无法工作或产生非常混乱的结果。
   * 不幸的是，当前对象没有办法找出它是否真的处于域的边界，因此无法确定你试图设置的值在当前情况下是否有意义。
   * @ingroup boundary   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  void
  set_boundary_id(const types::boundary_id);

  /**
   * 设置该顶点的流形指标。到目前为止，这没有任何作用，因为流形只用于细化和映射对象，但顶点没有被细化，映射是微不足道的。这个函数在这里只是为了允许独立维度的编程。
   *
   */
  void
  set_manifold_id(const types::manifold_id);

  /**
   * 设置此对象和其所有低维子对象的边界指标。
   * 因为这个对象只代表一个顶点，所以没有低维对象，这个函数等同于用相同的参数调用set_boundary_id()。
   * @ingroup boundary   @see   @ref GlossBoundaryIndicator  "关于边界指标的词汇条目"
   *
   */
  void
  set_all_boundary_ids(const types::boundary_id);

  /**
   * 设置此对象及其所有低维子对象的流形指标。
   * 由于这个对象只代表一个顶点，所以没有低维对象，这个函数等同于用相同的参数调用set_manifold_id()。
   * @ingroup manifold   @see   @ref GlossManifoldIndicator  "关于流形指标的词汇条目"
   *
   */
  void
  set_all_manifold_ids(const types::manifold_id);
  /**
   * @}
   *
   */

  /**
   * 返回这里所指向的顶点是否被使用。
   *
   */
  bool
  used() const;

  /**
   * 参考当前对象的单元格类型。
   *
   */
  ReferenceCell
  reference_cell() const;

  /**
   * 顶点的数量。
   *
   */
  unsigned int
  n_vertices() const;

  /**
   * 行的数量。
   *
   */
  unsigned int
  n_lines() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_vertices()的所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  vertex_indices() const;

  /**
   * 返回一个对象，它可以被认为是一个包含从零到n_lines()所有索引的数组。
   *
   */
  std_cxx20::ranges::iota_view<unsigned int, unsigned int>
  line_indices() const;

protected:
  /**
   * 指向我们所操作的三角形的指针。
   *
   */
  const Triangulation<1, spacedim> *tria;

  /**
   * 这是一个左端、右端或内部顶点。这个信息是在创建时由单元格提供的。
   *
   */
  VertexKind vertex_kind;

  /**
   * 这个对象所对应的顶点的全局顶点索引。
   *
   */
  unsigned int global_vertex_index;
};



/**
 * 这个类允许访问一个单元：一维的线，二维的四边形，等等。
 * 以下指的是任何维度。
 * 该类允许访问一个<tt>单元格</tt>，在一维中是一条线，在二维中是一个四角形。单元比线或四边形本身有更多的功能，例如，它们可以被标记为细化，它们有邻居，它们有可能检查它们是否在边界上等等。该类提供了对所有这些数据的访问。
 *
 *
 * @ingroup grid
 * @ingroup Accessors
 *
 */
template <int dim, int spacedim = dim>
class CellAccessor : public TriaAccessor<dim, dim, spacedim>
{
public:
  /**
   * 将AccessorData类型传播到本类中。
   *
   */
  using AccessorData = typename TriaAccessor<dim, dim, spacedim>::AccessorData;

  /**
   * 定义这个容器的类型，是它的一部分。
   *
   */
  using Container = Triangulation<dim, spacedim>;

  /**
   * @name  构造函数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 构建器。
   *
   */
  CellAccessor(const Triangulation<dim, spacedim> *parent     = nullptr,
               const int                           level      = -1,
               const int                           index      = -1,
               const AccessorData *                local_data = nullptr);

  /**
   * 复制构造器。
   *
   */
  CellAccessor(const TriaAccessor<dim, dim, spacedim> &cell_accessor);

  /**
   * 转换构造器。这个构造器的存在是为了使某些构造在独立于维度的代码中写得更简单。例如，它允许将一个面的迭代器分配给一个线的迭代器，这个操作在2D中很有用，但在3D中没有任何意义。这里的构造函数是为了使代码符合C++的要求而存在的，但它会无条件地中止；换句话说，将面迭代器赋值给线迭代器最好放在一个if语句中，检查维度是否为2，并在3D中赋值给一个四维迭代器（如果没有这个构造函数，如果我们碰巧为2d编译，这个操作是非法的）。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  CellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 另一个对象之间的转换操作符，就像之前的那个一样，没有意义。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  CellAccessor(const TriaAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 复制构造器。
   *
   */
  CellAccessor(const CellAccessor<dim, spacedim> &) = default;

  /**
   * 移动构造函数。
   *
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  CellAccessor(CellAccessor<dim, spacedim> &&) = default;

  /**
   * 解除构造器。
   *
   */
  ~CellAccessor() = default;

  /**
   * 复制操作符。这些操作符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的情况下使用。据推测，这里的意图是将 @p b
   * 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个存取器；因此，这个操作对三角形的迭代器没有用。
   * 因此，这个操作被声明为删除，不能使用。
   *
   */
  CellAccessor<dim, spacedim> &
  operator=(const CellAccessor<dim, spacedim> &) = delete;

  /**
   * 移动赋值运算符。
   *
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  CellAccessor<dim, spacedim> &
  operator=(CellAccessor<dim, spacedim> &&) = default; // NOLINT

  /**
   * @}
   *
   */

  /**
   * @name  访问子对象和相邻对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回一个指向 @p ith
   * 子对象的指针。重载版本，返回一个更合理的迭代器类。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  child(const unsigned int i) const;

  /**
   * 返回该单元格所有子代的迭代器的数组。
   *
   */
  boost::container::small_vector<TriaIterator<CellAccessor<dim, spacedim>>,
                                 GeometryInfo<dim>::max_children_per_cell>
  child_iterators() const;

  /**
   * 返回此单元格的 @p ith 面的一个迭代器。
   *
   */
  TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>
  face(const unsigned int i) const;

  /**
   * 返回当前单元格上 @p face 的面数。这是 TriaAccessor::face().
   * 的反函数。
   *
   */
  unsigned int
  face_iterator_to_index(
    const TriaIterator<TriaAccessor<dim - 1, dim, spacedim>> &face) const;

  /**
   * 返回该单元格所有面的迭代器阵列。
   *
   */
  boost::container::small_vector<
    TriaIterator<TriaAccessor<dim - 1, dim, spacedim>>,
    GeometryInfo<dim>::faces_per_cell>
  face_iterators() const;

  /**
   * 返回此单元格的 @p ith 面的（全局）索引。
   * @note
   * 尽管有这个名字，这里返回的索引只是全局的，即它是特定于一个特定的三角形对象的，或者，在三角形实际上是
   * parallel::distributed::Triangulation,
   * 类型的情况下，特定于存储在当前处理器上的分布式三角形的那一部分。
   *
   */
  unsigned int
  face_index(const unsigned int i) const;

  /**
   * 返回一个单元格的迭代器，该单元格与给定面和子面编号上的当前单元格相邻。
   * 为了成功，目前的单元必须没有被进一步细化，而在给定面的邻居必须被进一步细化一次；然后返回的单元是该邻居的一个子单元。
   * 这个函数不能在1d中调用，因为在那里我们没有子面。
   * 这个函数在2d中的实现相当直接，首先确定当前单元格与邻居单元格的哪个面接壤（这就是
   * @p neighbor_of_neighbor 函数的作用），然后向 @p
   * GeometryInfo::child_cell_on_subface 询问子单元的索引。
   * 然而，在3D中情况更为复杂，因为面可能有不止一个方向，我们必须对这个单元和邻近的单元使用
   * @p face_orientation,  @p face_flip 和 @p face_rotation
   * ，以弄清我们想拥有哪个单元。
   * 这可能会导致令人惊讶的结果：如果我们坐在一个单元格上，并要求得到子面<tt>sf</tt>后面的单元格，那么这意味着我们考虑的是本单元格自然方向上的脸的子面。然而，如果从这个单元格看到的面有<tt>face_orientation()==false</tt>，那么将本单元格与邻近单元格的子面分开的子面不一定是本单元格的
   * @p sf-th 子面。之所以如此，是因为单元格上的 @p
   * subface_no
   * 对应于相对于本单元格的内在排序的子面，而面的迭代器的子女是相对于面的内在排序来计算的；这两个排序只有在面的方向是
   * @p true, 时才是相同的，否则就会颠倒。
   * 同样，<tt>face_flip()==true</tt>和<tt>face_rotation()==true()</tt>的影响也要考虑，这两种情况都表示非标准的脸。
   * 幸运的是，这只是非常少的关注，因为通常我们只希望在一个活动单元的给定面的所有更细的邻居上进行循环。只有在细化三角图的过程中，我们才希望为我们的子单元和邻居的子单元设置邻居信息。因为在这种情况下，我们可以尊重当前单元格的面的方向，在这个函数中，我们不尊重当前单元格的面的方向、面的翻转和面的旋转，即返回的邻居的孩子在关于给定面的内在排序的子面
   * @p subface 后面。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor_child_on_subface(const unsigned int face_no,
                            const unsigned int subface_no) const;

  /**
   返回编号为 @p face_no. 的面的另一侧的邻接单元的迭代器 如果邻接单元不存在，即如果当前对象的编号为 @p face_no 的面处于边界，那么将返回一个无效的迭代器。    因此，索引 @p face_no 必须小于n_faces（）。    一个单元格的邻居最多具有与这个单元格相同的级别。例如，考虑以下情况。    @image html limit_level_difference_at_vertices.png ""  在这里，如果你在右上角的单元格上，并要求得到它的左邻（根据GeometryInfo类中阐明的惯例，就是它的<i>zeroth</i>邻），那么你将得到左上角四个小单元格的母格。换句话说，你作为邻居得到的单元格与你现在所在的单元格（右上角的那个）具有相同的细化级别，而且它可能有子代。    另一方面，如果你在左上角四个小单元格中的右上角，并且你要求获得右邻（与索引 <code>face_no=1</code> 相关），那么你会得到右上角的大单元格，在这种情况下，它的细化级别较低，并且没有自己的孩子。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor(const unsigned int face_no) const;

  /**
   * 返回索引为 @p face_no.
   * 的面的另一侧的相邻单元格的单元格索引
   * 如果相邻单元格不存在，该函数返回
   *
   * - .     这个函数等同于<tt>cell->neighbor(face_no)->index()</tt>。  更多细节请参见 neighbor()。
   *
   */
  int
  neighbor_index(const unsigned int face_no) const;

  /**
   * 返回编号为 @p face_no. 的面的另一边的相邻单元的水平
   * 如果相邻单元不存在，该函数返回
   *
   * - .     这个函数等同于`cell->neighbor(face_no)->level()`。  更多细节请参见 neighbor()。
   *
   */
  int
  neighbor_level(const unsigned int face_no) const;

  /**
   * 返回这个单元格是<tt>cell->neighbor(face_no)</tt>的第几个邻居，即返回
   * @p other_face_no
   * ，使得<tt>cell->neighbor(face_no)->neighbor(other_face_no)==cell</tt>。如果你想知道如何从邻居回到现在的单元格，这个函数是正确的。
   * 请注意，这个操作只有在邻居不比当前单元格粗的情况下才有用。如果邻居更粗，这个函数会抛出一个异常。在这种情况下，请使用
   * @p neighbor_of_coarser_neighbor 函数。
   *
   */
  unsigned int
  neighbor_of_neighbor(const unsigned int face_no) const;

  /**
   * 返回，邻居是否比现在的单元更粗大。这在各向异性的细化中很重要，因为这一信息并不取决于单元格的级别。
   * 请注意，在各向异性的情况下，一个单元只能在给定的面比另一个单元更粗，而不是在一般的基础上。较细的单元的面包含在较粗的单元的相应面中，较细的面是较粗的面的子女或孙子。
   *
   */
  bool
  neighbor_is_coarser(const unsigned int face_no) const;

  /**
   * 这个函数是 @p neighbor_of_neighbor
   * 函数的泛化，适用于粗邻的情况。它返回一对数字，face_no和subface_no，如果邻居没有被细化，则具有以下属性。<tt>cell->neighbor(
   * neighbor)->neighbor_child_on_subface(face_no,
   * subface_no)==cell</tt>。在3D中，一个更粗的邻居仍然可以被细化。
   * 在这种情况下，subface_no表示与我们的脸有关的邻居脸的子索引。
   * <tt>cell->neighbor(neighbor)->face(face_no)->child(subface_no)==cell->face(neighbor)</tt>。
   * 在 step-30
   * 教程程序的介绍中讨论了3D中的这种情况以及它如何发生。
   * 这个函数对于<tt>dim==1</tt>来说是不可能的。
   *
   */
  std::pair<unsigned int, unsigned int>
  neighbor_of_coarser_neighbor(const unsigned int neighbor) const;

  /**
   * 这个函数是 @p neighbor_of_neighbor 和 @p
   * neighbor_of_coarser_neighbor
   * 函数的一个概括。它检查邻居是否更粗，并调用相应的函数。在这两种情况下，只有face_no被返回。
   *
   */
  unsigned int
  neighbor_face_no(const unsigned int neighbor) const;

  /**
   * 与DoFCellAccessor兼容的接口。总是返回 @p false. 。
   *
   */
  static bool
  is_level_cell();

  /**
   * @}
   *
   */
  /**
   * @name  处理周期性邻居的问题
   *
   */
  /**
   * @{
   *
   */
  /**
   * 如果单元格在其 @c
   * 第i个面有周期性邻居，该函数返回真，否则，返回值为假。
   *
   */
  bool
  has_periodic_neighbor(const unsigned int i) const;

  /**
   * 对于其 @c 第i个面位于周期性边界的单元格，见 @ref GlossPeriodicConstraints "周期性边界的条目"
   * ，该函数返回周期性边界另一侧的单元格的迭代器。如果在
   * @c 第i个面没有周期性边界，将抛出一个异常。
   * 为了避免遇到异常，在使用这个函数之前，请检查has_periodic_neighbor()对
   * @c 第i个面的结果。
   * periodic_neighbor()的行为与neighbor()类似，即返回的单元最多具有与当前单元相同的细化程度。在分布式网格上，通过调用
   * Triangulation::add_periodicity(),
   * ，我们可以确保周期性边界另一侧的元素在这个等级中以幽灵单元或局部拥有的单元存在。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  periodic_neighbor(const unsigned int i) const;

  /**
   * 对于 @c
   * 第i个面不在边界上的单元格，这个函数返回的结果与邻接()相同。如果
   * @c
   * 第i个面在一个周期性的边界上，这个函数返回与periodic_neighbor()相同的结果。如果上述两个条件都不满足，即
   * @c 第i个面在一个非周期性边界上，将抛出一个异常。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  neighbor_or_periodic_neighbor(const unsigned int i) const;

  /**
   * 返回单元格在给定面和子面编号处的周期性邻居的一个迭代器。使用这个函数的一般准则类似于函数
   * neighbor_child_on_subface()。这个函数的实现与periodic_neighbor_of_coarser_periodic_neighbor()一致。例如，假设我们坐在一个名为
   * @c cell1的单元格上，其在 @c
   * 第1个面后面的邻居是一个更粗的层次。让我们把这个更粗的邻居命名为
   * @c cell2。然后，通过调用
   * periodic_neighbor_of_coarser_periodic_neighbor()，从  @c
   * cell1，我们得到一个  @c  face_num 和一个  @c
   * subface_num。现在，如果我们从cell2调用periodic_neighbor_child_on_subface()，用上述face_num和subface_num，我们会得到一个前往
   * @c  cell1的迭代器。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  periodic_neighbor_child_on_subface(const unsigned int face_no,
                                     const unsigned int subface_no) const;

  /**
   * 这个函数是periodic_neighbor_of_periodic_neighbor()的泛化，用于那些有较粗的周期性邻居的单元。返回的一对数字可以在period_neighbor_child_on_subface()中使用，以回到当前单元。换句话说，对于有较粗的周期性邻居的单元格，下面的断言应该是真的：cell->periodic_neighbor(i)->periodic_neighbor_child_on_subface(face_no,
   * subface_no)==cell
   *
   */
  std::pair<unsigned int, unsigned int>
  periodic_neighbor_of_coarser_periodic_neighbor(const unsigned face_no) const;

  /**
   * 这个函数返回当前单元格的第 @c
   * 个面的周期性邻居的索引。如果在给定的面没有周期性邻居，返回值为
   *
   * - .
   *
   */
  int
  periodic_neighbor_index(const unsigned int i) const;

  /**
   * 该函数返回当前单元格的 @c
   * 第i个面的周期性邻居的水平。如果在给定的面没有周期性邻居，返回值为
   *
   * - .
   *
   */
  int
  periodic_neighbor_level(const unsigned int i) const;

  /**
   * 对于在其 @c
   * 第i个面有周期性邻居的单元格，该函数返回该周期性邻居的面数，以便当前单元格是该邻居的周期性邻居。换句话说，对于那些拥有与当前单元相同或更高细化程度的周期性邻居的单元，以下断言成立。
   * @c {cell->periodic_neighbor(i)->
   * periodic_neighbor(cell->periodic_neighbor_of_periodic_neighbor(i))==cell}。
   * 对于具有较粗的周期性邻居的单元格，应该使用periodic_neighbor_of_coarser_periodic_neighbor()和periodic_neighbor_child_on_subface()来回到当前单元格。
   *
   */
  unsigned int
  periodic_neighbor_of_periodic_neighbor(const unsigned int i) const;

  /**
   * 如果一个单元格在其 @c
   * 第i个面有一个周期性邻居，这个函数返回周期性邻居的面数，它与这个单元格的
   * @c 第i个面相连。
   *
   */
  unsigned int
  periodic_neighbor_face_no(const unsigned int i) const;

  /**
   * 如果周期性边界另一侧的元素更粗，该函数返回真，否则返回假。该实现允许该函数在各向异性细化的情况下工作。
   *
   */
  bool
  periodic_neighbor_is_coarser(const unsigned int i) const;

  /**
   * @}
   *
   */

  /**
   * @name  处理边界指标的问题
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回 @p ith
   * 的顶点或面（取决于维度）是否是边界的一部分。如果
   * @p ith 的邻居不存在，这就是真的。
   *
   */
  bool
  at_boundary(const unsigned int i) const;

  /**
   * 返回该单元格是否在边界上。在边界上的定义是有一个面在边界上。请注意，这并不包括四边形或六边形中只有一个顶点在边界上的情况，或者六边形中只有一条线在边界上，而所有面的内部都在域的内部。对于后一种情况，
   * @p has_boundary_lines函数才是正确的请求。
   *
   */
  bool
  at_boundary() const;

  /**
   * 这是对 @p at_boundary
   * 函数的一个轻微变化：对于一维和二维，它是等价的，对于三维，它返回六面体的12条线中是否至少有一条位于边界。当然，这包括整个面处于边界的情况，但也包括其他一些情况。
   *
   */
  bool
  has_boundary_lines() const;
  /**
   * @}
   *
   */

  /**
   * @name  处理细化指标的问题
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回此单元格被标记为细化的 @p RefinementCase<dim> 。
   * 这个函数的返回值可以与一个bool进行比较，以检查这个单元格是否被标记为任何种类的细化。例如，如果你之前为一个单元格调用了cell->set_refine_flag()，那么你将进入以下片段中的'if'块。
   * @code
   * if (cell->refine_flag_set())
   * {
   * // yes, this cell is marked for refinement.
   * }
   * @endcode
   *
   *
   */
  RefinementCase<dim>
  refine_flag_set() const;

  /**
   * 对指向的单元格进行标记，以便进行细化。这个函数只允许用于活动单元。保持
   * @p ref_case 的默认值将标记该单元为各向同性的细化。
   * 如果你选择各向异性的细化，例如通过传递一个标志
   * RefinementCase::cut_x,   RefinementCase::cut_y,   RefinementCase::cut_z,
   * 或这些标志的组合作为参数，那么请记住，X、Y或Z方向的细化发生在该单元的
   * <em> 局部 </em>
   * 坐标系。换句话说，这些标志决定了单元格的哪些边和面将被切割成新的边和面。另一方面，这个过程与单元在
   * <em> 全局 </em>
   * 坐标系中的方向无关，你不应该假定单元的局部坐标系在它所处空间的全局坐标系中的任何特定方向。
   *
   */
  void
  set_refine_flag(const RefinementCase<dim> ref_case =
                    RefinementCase<dim>::isotropic_refinement) const;

  /**
   * 清除细化标志。
   *
   */
  void
  clear_refine_flag() const;

  /**
   * 修改单元格的细化标志，以确保（至少）在面<tt>face_no</tt>的给定细化情况
   * @p face_refinement_case
   * ，考虑到面的方向、翻转和旋转。返回，是否必须修改细化标志。这个函数只允许用于活动单元。
   *
   */
  bool
  flag_for_face_refinement(
    const unsigned int             face_no,
    const RefinementCase<dim - 1> &face_refinement_case =
      RefinementCase<dim - 1>::isotropic_refinement) const;

  /**
   * 修改单元格的细化标志，确保行<tt>面_no</tt>将被细化。返回，是否必须修改细化标志。这个函数只允许用于活动单元。
   *
   */
  bool
  flag_for_line_refinement(const unsigned int line_no) const;

  /**
   * 返回面<tt>face_no</tt>的SubfaceCase。请注意，这与询问<tt>cell->face(face_no)->refinement_case()</tt>不一样，因为后者返回一个RefinementCase<dim-1>，因此只考虑一个（各向异性的）细化，而这个函数考虑完整的细化情况，包括可能对面的孩子进行细化。这个函数只能在2D和3D的活动单元中调用。
   *
   */
  dealii::internal::SubfaceCase<dim>
  subface_case(const unsigned int face_no) const;

  /**
   * 返回粗化标志是否被设置。
   *
   */
  bool
  coarsen_flag_set() const;

  /**
   * 指向粗化的单元格的标志。这个函数只允许用于活动单元。
   *
   */
  void
  set_coarsen_flag() const;

  /**
   * 清除粗化标志。
   *
   */
  void
  clear_coarsen_flag() const;
  /**
   * @}
   *
   */

  /**
   * @name  处理材料指标
   *
   */
  /**
   * @{
   *
   */

  /**
   *
   */
  types::material_id
  material_id() const;

  /**
   *
   */
  void
  set_material_id(const types::material_id new_material_id) const;

  /**
   * 将此单元格及其所有子单元格（以及孙子单元格，以此类推）的材质ID设置为给定值。    参见 @ref GlossMaterialId "词汇表 "
   * 以了解更多信息。
   *
   */
  void
  recursively_set_material_id(const types::material_id new_material_id) const;
  /**
   * @}
   *
   */

  /**
   * @name  处理子域指标的问题
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回这个单元格的子域指标。    参见 @ref GlossSubdomainId "词汇表 "
   * 以了解更多信息。
   * @note
   * 单元格的子域是一个只为活动单元格定义的属性，即没有被进一步细化的单元。因此，只有当它所指的单元格没有子域时，你才能调用这个函数。对于并行的多网格方法，知道哪个处理器拥有非活动单元也很重要，为此你可以调用level_subdomain_id()。
   *
   */
  types::subdomain_id
  subdomain_id() const;

  /**
   * 设置这个单元的子域id。    参见 @ref GlossSubdomainId "词汇表 "
   * 以了解更多信息。如果你使用
   * parallel::distributed::Triangulation 对象，则不应调用此函数。
   * @note
   * 单元的子域是一个只为活动单元定义的属性，即没有被进一步细化的单元。因此，只有当它所指的单元格没有子域时，你才能调用这个函数。对于并行的多网格方法，知道哪个处理器拥有非活动单元也很重要，为此你可以调用level_subdomain_id()。
   *
   */
  void
  set_subdomain_id(const types::subdomain_id new_subdomain_id) const;

  /**
   * 获取该单元的水平子域ID。这用于并行多网格，其中不仅全局网格（由活动单元组成）被划分到处理器中，而且还包括构成网格的递归细化单元的各个层次。换句话说，如果使用多网格层次结构，层次子域id是一个也为非活动单元定义的属性。
   *
   */
  types::subdomain_id
  level_subdomain_id() const;

  /**
   * 设置此单元格的水平子域ID。这用于并行多网格。
   *
   */
  void
  set_level_subdomain_id(
    const types::subdomain_id new_level_subdomain_id) const;


  /**
   * 设置此单元的子域ID（如果它是活动的）或其所有终端子单元（和孙子单元，等等，只要它们没有自己的子单元）的给定值。由于子域id是一个只为活跃的单元格（即没有自己的孩子）定义的概念，这个函数只为这个单元格的所有实际活跃的子和孙子设置子域id，跳过中间的子单元格。    更多信息请参见  @ref GlossSubdomainId  "词汇表"
   * 。如果你使用 parallel::distributed::Triangulation
   * 对象，就不应该调用这个函数，因为在那里，子域id是由你所在的处理器隐式定义的。
   *
   */
  void
  recursively_set_subdomain_id(
    const types::subdomain_id new_subdomain_id) const;
  /**
   * @}
   *
   */

  /**
   * 为当前单元格返回一个全局唯一的单元格索引，假设它不是人工的。如果该单元是串行三角形的一部分，其值与active_cell_index()相同。
   * 在平行三角剖分的情况下，本地拥有的单元在网格的每个子域中被连续地列举出来。这就保证了这个函数返回的索引可以作为总条目数为
   * Triangulation::n_globally_active_cells()
   * 的向量的索引，并且每个进程都存储一个连续的部分。
   * 如果这样一个单元数据向量已经被设置为
   * parallel::TriangulationBase::global_active_cell_index_partitioner(),
   * ，那么这个函数返回的索引就可以用来访问正确的向量条目。
   *
   */
  types::global_cell_index
  global_active_cell_index() const;

  /**
   * 为非人工水平单元返回一个全局唯一的索引。
   * @note 与global_active_cell_index()类似。
   *
   */
  types::global_cell_index
  global_level_cell_index() const;

  /**
   * @name  处理codim 1单元的方位问题。
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回此单元格的方向。
   * 关于这个标志的含义，见  @ref GlossDirectionFlag  。
   *
   */
  bool
  direction_flag() const;

  /**
   * 返回当前单元格是第几个活动单元格（假设当前单元格确实是活动的）。这很有用，例如，如果你要访问一个有多少个条目的向量的元素，就有多少个活动单元。这样的向量用于估计三角形的每个单元的误差，用于指定传递给GridRefinement中的函数的细化标准，以及用于生成单元的输出。
   * 如果当前单元格没有被激活，该函数会抛出一个异常。
   * @note
   * 如果这个函数被调用的三角形是
   * parallel::distributed::Triangulation,
   * 类型，那么活动单元可能是本地拥有的、幽灵单元或人工的（见
   * @ref GlossLocallyOwnedCell 、 @ref GlossGhostCell 和 @ref
   * GlossArtificialCell  ）。
   * 这个函数对所有这些单元进行计数，包括幽灵和人工活动单元。这意味着该函数返回的索引可以唯一地识别单个处理器上三角结构中的单元，但不能唯一地识别处理器之间共享的三角结构（部分）中的单元。如果你想识别跨处理器的活动单元，你需要考虑由
   * CellAccessor::id(). 返回的单元的CellId。
   *
   */
  unsigned int
  active_cell_index() const;

  /**
   * 返回该单元格的父单元格在父单元格所属的三角结构层次中的索引。父单元的层次当然要比本单元的层次低一个。如果父单元不存在（即，如果该对象处于网格层次结构的最粗层），将产生一个异常。
   *
   */
  int
  parent_index() const;

  /**
   * 返回一个到父对象的迭代器。如果父对象不存在（即，如果该对象处于网格层次结构的最粗层），将产生一个异常。
   *
   */
  TriaIterator<CellAccessor<dim, spacedim>>
  parent() const;

  /**
   * @}
   *
   */

  /**
   * @name  其他函数
   *
   */
  /**
   * @{
   *
   */

  /**
   *
   */
  bool
  is_active() const;

  /**
   * 返回该单元是否为当前处理器所拥有，或为其他处理器所拥有。如果应用于类型为 dealii::Triangulation, 的对象，该函数总是返回true，但如果三角形是类型为 parallel::distributed::Triangulation. ，则可能产生false。更多信息请参见 @ref GlossGhostCell "词汇表 "
   * 和 @ref distributed 模块。      @post
   * 返回值等于<code>!is_ghost() && !is_artificial()/code>。
   * @note
   * 一个细胞是否是幽灵细胞、人造的，或者是本地拥有的，或者是一个只与活动的细胞有关的属性。因此，只有当它所指的单元格没有孩子时，你才能调用这个函数。
   *
   */
  bool
  is_locally_owned() const;

  /**
   * 如果三角结构没有分布，或者level_subdomain_id()等于当前处理器的id，则返回true。
   *
   */
  bool
  is_locally_owned_on_level() const;

  /**
   * 返回这个单元在全局网格中是否存在，但是(i)被另一个处理器所拥有，也就是有一个不同于当前处理器所拥有的子域_id，以及(ii)与当前处理器拥有的单元相邻。    这个函数只有在使用的三角形是 parallel::distributed::Triangulation. 的情况下才有意义，在所有其他情况下，返回值总是错误的。    更多信息请参见 @ref GlossGhostCell  "词汇表 "
   * 和 @ref distributed  模块。      @post
   * 返回值等于<code>!is_locally_owned() && !is_artificial()/code>。
   * @note
   * 一个细胞是否是幽灵细胞、人造的，或者是本地拥有的，或者是一个只与活动的细胞有关的属性。因此，只有当它所指的单元格没有孩子时，你才能调用这个函数。
   *
   */
  bool
  is_ghost() const;

  /**
   * 返回这个单元格是否是人造的，即它不是当前处理器所拥有的单元格之一，而且也不在一个单元格的边界上。因此，它存在于网格中，以确保每个处理器拥有所有粗略的网格单元，并保持相邻单元的2：1比例，但它不是我们应该在当前处理器上工作的单元之一。特别是，不能保证这个单元实际上没有在其他处理器上进一步细化。    这个函数只有在使用的三角形是 parallel::distributed::Triangulation. 的情况下才有意义，在所有其他情况下，返回值总是假的。    参见 @ref GlossArtificialCell "词汇表 "
   * 和 @ref distributed 模块以了解更多信息。      @post
   * 返回值等于<code>!is_ghost() && !is_locally_owned()/code>。
   * @note
   * 一个单元是否是幽灵单元、人造的或本地拥有的是一个只与活动的单元有关的属性。因此，你只能在它所指的单元格没有孩子的情况下调用这个函数。
   *
   */
  bool
  is_artificial() const;

  /**
   * 测试点 @p p
   * 是否在这个单元格内。边界上的点被算作是在单元格内。
   * 请注意，这个函数假定单元格和实际单元格之间的映射是（双，三）线性的，也就是说，2D的面和3D的边都是直线。如果你有更高阶的变换，结果可能不同，因为一个点在实空间中是在单元内还是在单元外。
   * 在codim>0的情况下，首先将点投影到单元格所嵌入的流形上，然后检查这个投影是否在单元格内。
   *
   */
  bool
  point_inside(const Point<spacedim> &p) const;

  /**
   * 将此单元格的邻居 @p i 设置为 @p pointer. 所指向的单元格
   * 这个函数其实不应该公开（但由于各种原因需要公开，以便不使一长串函数成为朋友）：它修改了内部数据结构，可能会留下一些东西。请不要从应用程序代码中使用它。
   *
   */
  void
  set_neighbor(const unsigned int                               i,
               const TriaIterator<CellAccessor<dim, spacedim>> &pointer) const;

  /**
   * 返回当前单元格的唯一ID。这个ID是由层次结构中从粗父单元开始的路径构建的，并且在使用类型为
   * parallel::distributed::Triangulation.
   * 的对象的并行计算中正确工作。因此，这个函数在为单元（活动或不活动）提供一个唯一的标识符方面非常有用，也适用于并行三角计算。更多信息请参见CellId类的文档。
   * @note
   * 这个操作需要O(level)时间来计算。在大多数实际情况下，三角形的层数将取决于三角形中的单元格数量的对数。
   *
   */
  CellId
  id() const;

  using TriaAccessor<dim, dim, spacedim>::diameter;

  /**
   * 与 TriaAccessor::diameter() 相同，但也采取映射类。
   *
   */
  double
  diameter(const Mapping<dim, spacedim> &mapping) const;

  /**
   * @}
   *
   */


  /**
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcRefineCellNotActive);
  /**
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcCellFlaggedForRefinement);
  /**
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcCellFlaggedForCoarsening);

protected:
  /**
   * 这个函数假设邻居不比当前单元格粗。在这种情况下，它返回
   * neighbor_of_neighbor()值。然而，如果邻居更粗，这个函数返回一个
   * <code>invalid_unsigned_int</code>  。
   * 这个函数不供公众使用。请使用函数
   * neighbor_of_neighbor()来代替，如果对一个较粗的邻居进行调用，会抛出一个异常。如果邻居确实较粗（你可以通过neigher_is_coarser()函数知道这一点），那么应该调用neigher_of_coarser_neighbor()函数。如果你只想知道从邻居到现在的单元所需的
   * <code>face_no</code> ，那么只需使用
   * neighbor_face_no()函数，该函数可用于较粗的邻居和非较粗的邻居。
   *
   */
  unsigned int
  neighbor_of_neighbor_internal(const unsigned int neighbor) const;

  /**
   * 至于任何codim>0，我们可以使用类似的代码，C++不允许部分模板。我们使用这个辅助函数，然后从point_inside调用。
   *
   */
  template <int dim_, int spacedim_>
  bool
  point_inside_codim(const Point<spacedim_> &p) const;



private:
  /**
   * 设置一个单元格的活动单元格索引。这是在细化结束时进行的。
   *
   */
  void
  set_active_cell_index(const unsigned int active_cell_index) const;

  /**
   * 设置一个单元格的全局活动单元格索引。
   *
   */
  void
  set_global_active_cell_index(const types::global_cell_index index) const;

  /**
   * 为一个水平单元设置全局水平单元索引。
   *
   */
  void
  set_global_level_cell_index(const types::global_cell_index index) const;

  /**
   * 设置一个单元格的父级。
   *
   */
  void
  set_parent(const unsigned int parent_index);

  /**
   * 设置该单元格的方向。
   * 关于这个标志的含义，见  @ref GlossDirectionFlag  。
   *
   */
  void
  set_direction_flag(const bool new_direction_flag) const;

  template <int, int>
  friend class Triangulation;

  template <int, int>
  friend class parallel::TriangulationBase;

  friend struct dealii::internal::TriangulationImplementation::Implementation;
  friend struct dealii::internal::TriangulationImplementation::
    ImplementationMixedMesh;
};



 /* ----- declaration of explicit specializations and general templates ----- */ 


template <int structdim, int dim, int spacedim>
template <typename OtherAccessor>
InvalidAccessor<structdim, dim, spacedim>::InvalidAccessor(
  const OtherAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim, spacedim>::CellAccessor(
  const InvalidAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int structdim, int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
TriaAccessor<structdim, dim, spacedim>::TriaAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



template <int dim, int spacedim>
template <int structdim2, int dim2, int spacedim2>
CellAccessor<dim, spacedim>::CellAccessor(
  const TriaAccessor<structdim2, dim2, spacedim2> &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}


#ifndef DOXYGEN

template <>
bool
CellAccessor<1, 1>::point_inside(const Point<1> &) const;
template <>
bool
CellAccessor<2, 2>::point_inside(const Point<2> &) const;
template <>
bool
CellAccessor<3, 3>::point_inside(const Point<3> &) const;
template <>
bool
CellAccessor<1, 2>::point_inside(const Point<2> &) const;
template <>
bool
CellAccessor<1, 3>::point_inside(const Point<3> &) const;
template <>
bool
CellAccessor<2, 3>::point_inside(const Point<3> &) const;
// -------------------------------------------------------------------

template <>
void
TriaAccessor<3, 3, 3>::set_all_manifold_ids(const types::manifold_id) const;

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

// include more templates in debug and optimized mode
#include "tria_accessor.templates.h"

#endif


