//include/deal.II-translator/dofs/dof_accessor_0.txt
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

#ifndef dealii_dof_accessor_h
#define dealii_dof_accessor_h


#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_iterator_selector.h>

#include <deal.II/grid/tria_accessor.h>

#include <deal.II/hp/dof_handler.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/container/small_vector.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename number>
class FullMatrix;
template <typename number>
class Vector;
template <typename number>
class AffineConstraints;

template <typename Accessor>
class TriaRawIterator;

template <int, int>
class FiniteElement;

namespace internal
{
  namespace DoFCellAccessorImplementation
  {
    struct Implementation;
  }

  namespace DoFHandlerImplementation
  {
    struct Implementation;
    namespace Policy
    {
      struct Implementation;
    }
  } // namespace DoFHandlerImplementation

  namespace hp
  {
    namespace DoFHandlerImplementation
    {
      struct Implementation;
    }
  } // namespace hp
} // namespace internal
#endif

// note: the file dof_accessor.templates.h is included at the end of
// this file.  this includes a lot of templates and thus makes
// compilation slower, but at the same time allows for more aggressive
// inlining and thus faster code.


namespace internal
{
  namespace DoFAccessorImplementation
  {
    /**
     * 这是一个开关类，它只声明一个  @p alias.
     * 它是用来确定DoFAccessor类是从哪个类派生的。默认情况下，<tt>DoFAccessor
     * @<structdim,dim,spacedim@></tt>  派生自一般的<tt>Inheritance
     * @<structdim,dim,spacedim@></tt>  类中的别名，即<tt>TriaAccessor
     * @<structdim,dim,spacedim@></tt>,
     * ，但如果<tt>structdim==dim</tt>，则使用专门化<tt>Inheritance
     * @<dim,dim,spacedim@></tt>
     * ]被使用，它声明其本地类型为<tt>CellAccessor
     * @<dim,spacedim@></tt>.
     * 因此，如果所考虑的对象有完整的尺寸，即构成一个单元，则自动选择继承自CellAccessor。
     * @ingroup dofs
     * @ingroup Accessors
     *
     */
    template <int structdim, int dim, int spacedim>
    struct Inheritance
    {
      /**
       * @p alias. 的声明 更多信息请参见完整的文档。
       *
       */
      using BaseClass = dealii::TriaAccessor<structdim, dim, spacedim>;
    };


    /**
     * 这是一般模板的特殊化，用于对象有全尺寸的情况，即是一个单元。更多细节见一般模板。
     *
     */
    template <int dim, int spacedim>
    struct Inheritance<dim, dim, spacedim>
    {
      /**
       * @p alias. 的声明 更多信息请参见完整的文档。
       *
       */
      using BaseClass = dealii::CellAccessor<dim, spacedim>;
    };

    struct Implementation;
  } // namespace DoFAccessorImplementation
} // namespace internal


 /* -------------------------------------------------------------------------- */ 



/**
 * 一个可以访问存储在DoFHandler对象中的自由度的类。访问器用于访问与三角形的边、面和单元有关的数据。这个概念在
 * @ref Iterators  中有更详细的解释。
 * 该类主要遵循三角形库中声明的访问器库（TriaAccessor）所规定的路线。它使用户能够访问线、四边形或六边形的自由度。这个类的第一个模板参数决定了所考虑的对象的维度。1用于线，2用于四边形，3用于六边形。第二个参数表示我们应该在哪个类型的DoFHandler上工作。从第二个模板参数中，我们还可以推导出这个对象所指向的三角形的维度，以及它所嵌入的空间的维度。最后，模板参数
 * <code>level_dof_access</code>
 * 制约着函数get_active_or_mg_dof_indices()的行为。参见下面关于通用循环的部分。
 * <h3>Alias</h3>
 * 使用方法最好是通过DoFHandler类提供的各种迭代器的别名来实现，因为它们对类的命名和模板接口的变化更安全，同时也提供了更简单的输入方式（更少的复杂名称！）。
 * <h3>Generic loops and the third template argument</h3>
 * 许多循环看起来非常相似，无论它们是对三角结构的活动单元的活动道次进行操作，还是对单层或整个网格层次的水平道次进行操作。为了在这类循环中使用多态性，它们通过函数get_active_or_mg_dof_indices()访问自由度，该函数根据第三个模板参数改变行为。
 * 如果该参数为false，那么将访问活动单元的活动自由度。如果它是true，则使用水平道夫。DoFHandler有一些函数，例如begin()和begin_mg()，它们返回任一类型或其他类型。此外，它们可以相互转换，如果需要的话，因为它们访问的是相同的数据。
 * 建议在通用循环中使用函数get_active_or_mg_dof_indices()来代替get_dof_indices()或get_mg_dof_indices()。
 * <h3>Inheritance</h3>
 * 如果第一个模板参数给出的结构维度等于DoFHandler的维度（作为第二个模板参数给出），那么我们显然是在处理单元，而不是低维的对象。在这种情况下，继承自CellAccessor，以提供对该类所提供的所有细胞特定信息的访问。否则，也就是说，对于低维对象，继承自TriaAccessor。
 * 有一个DoFCellAccessor类，提供了与CellAccessor类等价的功能。
 * @tparam  structdim 访问器代表的对象的维度。例如，点的 @p
 * structdim 等于0，边的 @p structdim 等于1，等等。  @tparam  dim
 * 底层DoFHandler的维度。  @tparam  spacedim
 * 底层DoFHandler的空间尺寸。  @tparam  level_dof_access 如果 @p
 * false,
 * ，则访问器简单地表示DoFHandler中的单元、面或边，对于这些单元、面或边，自由度只存在于最精细的层面上。在这种情况下，有些操作是不允许的，比如询问非活动单元的自由度指数。另一方面，如果这个模板参数是
 * @p true,
 * ，那么访问器代表自由度多级层次中的一个对象。在这种情况下，访问<i>any</i>单元的DoF指数是可能的，并将返回<i>level</i>指数（对于活动单元，可能与<i>global</i>指数不同）。
 *
 *
 * @ingroup dofs
 *
 * @ingroup Accessors
 *
 */
template <int structdim, int dim, int spacedim, bool level_dof_access>
class DoFAccessor : public dealii::internal::DoFAccessorImplementation::
                      Inheritance<structdim, dim, spacedim>::BaseClass
{
public:
  /**
   * 一个静态变量，允许该类的用户发现第二个模板参数的值。
   *
   */
  static const unsigned int dimension = dim;

  /**
   * 一个静态变量，允许这个类的用户发现第三个模板参数的值。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 声明一个基类的别名，使访问一些异常类更加简单。
   *
   */
  using BaseClass = typename dealii::internal::DoFAccessorImplementation::
    Inheritance<structdim, dimension, space_dimension>::BaseClass;

  /**
   * 迭代器类所传递的数据类型。
   *
   */
  using AccessorData = DoFHandler<dimension, space_dimension>;

  /**
   * @name  构造函数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 默认构造函数。提供一个不能使用的访问器。
   *
   */
  DoFAccessor();

  /**
   * 构造函数，生成一个指向DoFHandler中特定单元或面或边的访问。
   * @param  tria 这个访问器所指向的三角结构。    @param  level
   * 指向对象的网格层次结构中的级别。例如，粗略的网格单元有零级，它们的子层有一级，以此类推。对于没有等级的面和边，这个参数会被忽略。
   * @param  index 指向指定细化层的对象的索引。    @param
   * dof_handler
   * 指向访问器应引用的DoFHandler对象的指针。当然，这个DoFHandler对象必须建立在与第一个参数中指定的相同的三角形上。
   *
   */
  DoFAccessor(const Triangulation<dim, spacedim> *tria,
              const int                           level,
              const int                           index,
              const DoFHandler<dim, spacedim> *   dof_handler);

  /**
   * 复制构造器。
   *
   */
  DoFAccessor(const DoFAccessor<structdim, dim, spacedim, level_dof_access> &) =
    default;

  /**
   * 移动构造函数。
   *
   */
  DoFAccessor(                                                    // NOLINT
    DoFAccessor<structdim, dim, spacedim, level_dof_access> &&) = // NOLINT
    default;                                                      // NOLINT

  /**
   * 解除构造器。
   *
   */
  ~DoFAccessor() = default;

  /**
   * 转换构造器。这个构造器的存在是为了使某些构造在独立于维度的代码中写得更简单。例如，它允许将一个面的迭代器分配给一个线的迭代器，这个操作在2D中很有用，但在3D中没有任何意义。这里的构造函数是为了使代码符合C++的要求而存在的，但它会无条件地中止；换句话说，将面迭代器赋值给线迭代器最好放在一个if语句中，检查维度是否为2，并在3D中赋值给一个四维迭代器（如果没有这个构造函数，如果我们碰巧为2d编译，这个操作是非法的）。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  DoFAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 另一个对象之间的转换操作符，就像之前的那个一样，没有意义。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  DoFAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * 复制构造函数允许切换级别访问和主动访问。
   *
   */
  template <bool level_dof_access2>
  DoFAccessor(const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &);

  /**
   * 拷贝操作符。这些操作符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的情况下使用。据推测，这里的意图是将 @p b
   * 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个访问器；因此，这个操作对DoF处理程序对象的迭代器没有用。
   * 因此，这个操作被声明为删除，不能使用。
   *
   */
  DoFAccessor<structdim, dim, spacedim, level_dof_access> &
  operator=(const DoFAccessor<structdim, dim, spacedim, level_dof_access> &da) =
    delete;

  /**
   * 移动赋值运算符。
   *
   */
  DoFAccessor<structdim, dim, spacedim, level_dof_access> &       // NOLINT
  operator=(                                                      // NOLINT
    DoFAccessor<structdim, dim, spacedim, level_dof_access> &&) = // NOLINT
    default;                                                      // NOLINT

  /**
   * @}
   *
   */

  /**
   * 返回一个我们正在使用的DoFHandler对象的句柄。
   *
   */
  const DoFHandler<dim, spacedim> &
  get_dof_handler() const;

  /**
   * 实现迭代器类所需的复制操作。
   *
   */
  template <bool level_dof_access2>
  void
  copy_from(const DoFAccessor<structdim, dim, spacedim, level_dof_access2> &a);

  /**
   * 迭代器类所使用的复制运算器。保留之前设置的dof处理程序，但设置TriaAccessor的对象坐标。
   *
   */
  void
  copy_from(const TriaAccessorBase<structdim, dim, spacedim> &da);

  /**
   * 告诉调用者get_active_or_mg_dof_indices()是访问活动的还是水平的道夫。
   *
   */
  static bool
  is_level_cell();

  /**
   * @name  访问子对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回一个指向 @p c-th 子对象的迭代器。
   *
   */
  TriaIterator<DoFAccessor<structdim, dim, spacedim, level_dof_access>>
  child(const unsigned int c) const;

  /**
   * 指向与此对象相界的 @p ith
   * 线的指针。如果当前对象本身是一条线，那么唯一有效的索引是
   * @p i 等于零，并且该函数返回一个指向自身的迭代器。
   *
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::line_iterator
    line(const unsigned int i) const;

  /**
   * 指向与此对象相邻的 @p ith
   * 四边形的指针。如果当前对象本身是一个四边形，那么唯一有效的索引是
   * @p i 等于零，并且该函数返回一个指向自身的迭代器。
   *
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<dim, spacedim, level_dof_access>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问此对象的DoF索引
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回位于此对象上的自由度的<i>global</i>指数，以有限元定义的标准排序（即顶点0上的自由度，顶点1上的自由度，等等，行0上的自由度，行1上的自由度，等等，quad 0上的自由度，等等）此函数仅对<i>active</i>对象可用（见 @ref GlossActive "此词汇条"
   * ）。
   * 单元需要是一个活跃的单元（而不是平行分布式计算中的人工）。
   * 向量在传递给这个函数之前必须有合适的大小。
   * 最后一个参数表示有限元素的索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的同一有限元。
   * 然而，当相关的DoFHandler对象启用了hp-capabilities，不同的有限元对象可以在不同的单元格上使用。因此，在两个单元之间的面以及顶点上，可能有两组自由度，相邻单元上使用的每个有限元都有一个。为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   * 对于单元，只有一个可能的有限元指数（即该单元的指数，由
   * <code>cell-@>active_fe_index</code> 返回。
   * 因此，派生的DoFCellAccessor类有一个该函数的重载版本，它以
   * <code>cell-@>active_fe_index</code> 为最后参数调用本函数。
   *
   */
  void
  get_dof_indices(std::vector<types::global_dof_index> &dof_indices,
                  const unsigned int                    fe_index =
                    DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * 返回当前对象上的自由度的全局多级指数，相对于多网格层次结构中的给定层次而言。指数是指该行所处层次的本地编号。
   *
   */
  void
  get_mg_dof_indices(const int                             level,
                     std::vector<types::global_dof_index> &dof_indices,
                     const unsigned int                    fe_index =
                       DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * 设置由get_mg_dof_indices返回的层次DoF指数。
   *
   */
  void
  set_mg_dof_indices(
    const int                                   level,
    const std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index = DoFHandler<dim, spacedim>::invalid_fe_index);

  /**
   * 与当前单元的 @p vertexth
   * 顶点相关的<i>i</i>度的全局DoF指数。
   * 最后一个参数表示的是有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类反正只支持所有单元上的相同有限元。
   * 然而，当hp-capabilities被启用时，不同的有限元对象可以被用于不同的单元。因此，在两个单元之间的面以及顶点上，可能有两组自由度，相邻单元上使用的每个有限元都有一个。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  types::global_dof_index
  vertex_dof_index(const unsigned int vertex,
                   const unsigned int i,
                   const unsigned int fe_index =
                     DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * 返回与 @p level. 层的 <code>vertex</code> 个顶点相关的
   * <code>i</code> 个自由度的全局DoF索引。
   * 也可以参见vertex_dof_index()。
   *
   */
  types::global_dof_index
  mg_vertex_dof_index(const int          level,
                      const unsigned int vertex,
                      const unsigned int i,
                      const unsigned int fe_index =
                        DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * 这个对象的第<i>i</i>个自由度的索引。
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于它的默认值，因为该类反正只支持所有单元上的同一个有限元。
   * 然而，当hp-capabilities被启用时，不同的有限元对象可以被用于不同的单元。因此，在两个单元之间的面以及顶点上，可能有两组自由度，相邻单元上使用的每个有限元都有一个。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果相匹配。
   * @note
   * 虽然get_dof_indices()函数返回一个数组，其中包含以某种方式存在于这个对象上的所有自由度的索引（即在这个对象的顶点、边或内部），但当前的dof_index()函数只考虑真正属于这个特定对象内部的自由度。换句话说，举个例子，如果当前对象指的是一个四边形（2D中的单元，3D中的面），并且与之相关的有限元是双线性的，那么get_dof_indices()会返回一个大小为4的数组，而dof_index()会产生一个异常，因为在面的内部没有定义度。
   *
   */
  types::global_dof_index
  dof_index(const unsigned int i,
            const unsigned int fe_index =
              DoFHandler<dim, spacedim>::invalid_fe_index) const;

  /**
   * 返回给定层面上的dof_index。也见dof_index。
   *
   */
  types::global_dof_index
  mg_dof_index(const int level, const unsigned int i) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问与此对象相关的有限元
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回在给定对象上激活的有限元的数量。
   * 当hp-capabilities被禁用时，答案当然总是1。
   * 然而，当hp-capabilities被启用时，情况就不是这样了。如果这是一个单元，答案当然是1。如果它是一个面，答案可能是1或2，取决于相邻的两个单元是否使用相同的有限元。如果它是一个3D的边，可能的返回值可能是1或大于这个值的任何其他值。
   *
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * 返回此对象上的 @p n-th
   * 活动FE索引。对于单元格和所有非hp-objects，只有一个活跃的FE索引，所以参数必须等于0。对于低维的hp-objects，有n_active_fe_indices()活动有限元，这个函数可以查询它们的指数。
   *
   */
  unsigned int
  nth_active_fe_index(const unsigned int n) const;

  /**
   * 返回此对象上的所有活动FE指数。
   * 返回的集合的大小等于此对象上活动的有限元的数量。
   *
   */
  std::set<unsigned int>
  get_active_fe_indices() const;

  /**
   * 如果具有给定索引的有限元在当前对象上处于活动状态，则返回真。当当前DoFHandler没有hp-能力时，当然只有当
   * @p fe_index 等于0时才是这种情况。对于单元格来说，如果
   * @p fe_index
   * 等于该单元格的active_fe_index()，则是这种情况。对于面和其他低维物体，可能有一个以上的
   * @p fe_index
   * 在任何给定的物体上是活跃的（见n_active_fe_indices()）。
   *
   */
  bool
  fe_index_is_active(const unsigned int fe_index) const;

  /**
   * 返回给定 @p fe_index. 在此对象上使用的有限元的引用 @p
   * fe_index 必须在此对象上使用，即
   * <code>fe_index_is_active(fe_index)</code> 必须返回true。
   *
   */
  const FiniteElement<dim, spacedim> &
  get_fe(const unsigned int fe_index) const;

  /**
   * @}
   *
   */

  /**
   * 子类的例外情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcInvalidObject,
                   "This accessor object has not been "
                   "associated with any DoFHandler object.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcVectorNotEmpty);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcVectorDoesNotMatch);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcMatrixDoesNotMatch);
  /**
   * 为一个应该是 @ref GlossActive "活动 "
   * 的单元格调用了一个函数，但它被细化了。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNotActive);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcCantCompareIterators);

protected:
  /**
   * 存储要访问的DoFHandler对象的地址。
   *
   */
  DoFHandler<dim, spacedim> *dof_handler;

public:
  /**
   * 比较是否相等。如果两个访问器引用的是同一个对象，则返回<tt>true</tt>。
   * 这个函数的模板参数允许对非常不同的对象进行比较。因此，其中一些被禁用。也就是说，如果两个对象的尺寸或dof处理程序不同，会产生一个异常。可以预见，这是一个不需要的比较。
   * 模板参数<tt>level_dof_access2</tt>被忽略了，这样，一个有级别访问的迭代器可以等于一个有活动自由度访问的迭代器。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator==(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * 比较不等式。操作符==()的布尔值不是。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator!=(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

protected:
  /**
   * 重置DoF处理程序指针。
   *
   */
  void
  set_dof_handler(DoFHandler<dim, spacedim> *dh);

  /**
   * 将此对象的<i>i</i>个自由度的索引设置为 @p 个索引。
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的同一有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  void
  set_dof_index(const unsigned int            i,
                const types::global_dof_index index,
                const unsigned int            fe_index =
                  DoFHandler<dim, spacedim>::invalid_fe_index) const;

  void
  set_mg_dof_index(const int                     level,
                   const unsigned int            i,
                   const types::global_dof_index index) const;

  /**
   * 将当前单元的 @p vertex-th
   * 顶点上的<i>i</i>度的全局索引设置为 @p index.
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的相同有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  void
  set_vertex_dof_index(const unsigned int            vertex,
                       const unsigned int            i,
                       const types::global_dof_index index,
                       const unsigned int            fe_index =
                         DoFHandler<dim, spacedim>::invalid_fe_index) const;

  void
  set_mg_vertex_dof_index(const int                     level,
                          const unsigned int            vertex,
                          const unsigned int            i,
                          const types::global_dof_index index,
                          const unsigned int            fe_index =
                            DoFHandler<dim, spacedim>::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;
  template <int, int, int, bool>
  friend class DoFAccessor;

private:
  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  template <int, int>
  friend class DoFHandler;

  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
  friend struct dealii::internal::DoFAccessorImplementation::Implementation;
};



/**
 * 一般DoFAccessor类模板的特化，用于零维对象（顶点）的情况，这些对象是空间维度上的一维单元的面。由于顶点的功能与一般的面不同，这个类做了一些与一般模板不同的事情，但界面看起来应该是一样的。
 *
 *
 */
template <int spacedim, bool level_dof_access>
class DoFAccessor<0, 1, spacedim, level_dof_access>
  : public TriaAccessor<0, 1, spacedim>
{
public:
  /**
   * 一个静态变量，允许这个类的用户发现第二个模板参数的值。
   *
   */
  static const unsigned int dimension = 1;

  /**
   * 一个静态变量，允许这个类的用户发现第三个模板参数的值。
   *
   */
  static const unsigned int space_dimension = spacedim;

  /**
   * 声明一个基类的别名，使访问一些异常类更加简单。
   *
   */
  using BaseClass = TriaAccessor<0, 1, spacedim>;

  /**
   * 迭代器类所传递的数据类型。
   *
   */
  using AccessorData = DoFHandler<1, spacedim>;

  /**
   * @name  构造函数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 默认构造函数。提供一个不能使用的访问器。
   *
   */
  DoFAccessor();

  /**
   * 如果这里的对象指的是一维三角形的一个顶点，即三角形的一个面，则使用构造函数。
   * 由于没有从顶点到单元的映射，一个点的访问器对象没有办法弄清它是否在域的边界上。因此，第二个参数必须由生成这个访问器的对象来传递
   *
   * 例如，一个1d单元可以计算出它的左或右顶点是否在边界上。
   * 第三个参数是我们指向的顶点的全局索引。
   * 第四个参数是一个指向DoFHandler对象的指针。
   * 这个迭代器只能为一维三角计算调用。
   *
   */
  DoFAccessor(
    const Triangulation<1, spacedim> *                      tria,
    const typename TriaAccessor<0, 1, spacedim>::VertexKind vertex_kind,
    const unsigned int                                      vertex_index,
    const DoFHandler<1, spacedim> *                         dof_handler);

  /**
   * 构造函数。这个构造函数的存在是为了保持与其他访问器类的接口兼容性。然而，它在这里并没有做任何有用的事情，所以实际上可能不会被调用。
   *
   */
  DoFAccessor(const Triangulation<1, spacedim> *,
              const int                                  = 0,
              const int                                  = 0,
              const DoFHandler<1, spacedim> *dof_handler = 0);

  /**
   * 转换构造函数。这个构造函数的存在是为了使某些构造在独立于维度的代码中写得更简单。例如，它允许将一个面的迭代器分配给一个线的迭代器，这个操作在2D中很有用，但在3D中没有任何意义。这里的构造函数是为了使代码符合C++的要求而存在的，但它会无条件地中止；换句话说，将面迭代器分配给线迭代器最好放在一个if语句中，检查维度是否为2，并在3D中分配给一个四边形迭代器（如果没有这个构造函数，如果我们碰巧为2d编译，这个操作是非法的）。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  DoFAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 另一个对象之间的转换操作符，就像之前的那个一样，没有意义。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  DoFAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * 复制构造器。
   *
   */
  DoFAccessor(const DoFAccessor<0, 1, spacedim, level_dof_access> &) = default;

  /**
   * 移动构造函数。
   *
   */
  // NOLINTNEXTLINE OSX does not compile with noexcept
  DoFAccessor(DoFAccessor<0, 1, spacedim, level_dof_access> &&) = default;

  /**
   * 解除构造器。
   *
   */
  ~DoFAccessor() = default;

  /**
   * 复制操作符。这些操作符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的上下文中使用。据推测，这里的意图是将
   * @p b 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个访问器；因此，这个操作对DoF处理程序对象的迭代器没有用。
   * 因此，这个操作被声明为删除，不能使用。
   *
   */
  DoFAccessor<0, 1, spacedim, level_dof_access> &
  operator=(const DoFAccessor<0, 1, spacedim, level_dof_access> &da) = delete;

  /**
   * 移动赋值运算符。
   *
   */
  DoFAccessor<0, 1, spacedim, level_dof_access> &operator      =(
    DoFAccessor<0, 1, spacedim, level_dof_access> &&) noexcept = default;

  /**
   * @}
   *
   */

  /**
   * 返回一个我们正在使用的DoFHandler对象的句柄。
   *
   */
  const DoFHandler<1, spacedim> &
  get_dof_handler() const;

  /**
   * 实现迭代器类所需的复制操作。
   *
   */
  template <bool level_dof_access2>
  void
  copy_from(const DoFAccessor<0, 1, spacedim, level_dof_access2> &a);

  /**
   * 迭代器类所使用的复制运算器。保留之前设置的dof处理程序，但设置TriaAccessor的对象坐标。
   *
   */
  void
  copy_from(const TriaAccessorBase<0, 1, spacedim> &da);

  /**
   * @name  访问子对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回一个无效的迭代器，其类型代表指向当前对象的一个子对象。该对象是无效的，因为点（由当前类代表）没有子代。
   *
   */
  TriaIterator<DoFAccessor<0, 1, spacedim, level_dof_access>>
  child(const unsigned int c) const;

  /**
   * 指向与此对象相界的 @p ith 线的指针。
   * 由于维度为1的网格没有四边形，这个方法只是抛出一个异常。
   *
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::line_iterator
    line(const unsigned int i) const;

  /**
   * 指向包围此对象的 @p ith 四边形的指针。
   * 由于维度为1的网格没有四边形，这个方法只是抛出一个异常。
   *
   */
  typename dealii::internal::DoFHandlerImplementation::
    Iterators<1, spacedim, level_dof_access>::quad_iterator
    quad(const unsigned int i) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问此对象的DoF指数
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回位于此物体上的自由度的<i>global</i>指数，其标准排序由有限元定义。这个函数只适用于<i>active</i>对象（见 @ref GlossActive "本词汇条"
   * ）。
   * 目前的顶点必须属于一个活动单元（而不是在并行分布式计算中的人工）。
   * 向量在传递给这个函数之前必须有合适的大小。
   * 最后一个参数表示有限元素的索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的同一有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   * 对于单元，只有一个可能的有限元指数（即该单元的指数，由
   * <code>cell-@>active_fe_index</code> 返回。
   * 因此，派生的DoFCellAccessor类有一个该函数的重载版本，它以
   * <code>cell-@>active_fe_index</code> 为最后参数调用本函数。
   *
   */
  void
  get_dof_indices(
    std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * 返回当前对象上的自由度的全局多级指数，相对于多网格层次结构中的给定层次而言。指数是指该行所处层次的本地编号。
   *
   */
  void
  get_mg_dof_indices(
    const int                             level,
    std::vector<types::global_dof_index> &dof_indices,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * 与当前单元的 @p vertexth
   * 顶点相关的<i>i</i>度的全局DoF索引。
   * 最后一个参数表示的是有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的同一有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  types::global_dof_index
  vertex_dof_index(
    const unsigned int vertex,
    const unsigned int i,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * 此对象的<i>i</i>个自由度的索引。
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于它的默认值，因为该类反正只支持所有单元上的同一个有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  types::global_dof_index
  dof_index(const unsigned int i,
            const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问与此对象相关的有限元
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回在给定对象上激活的有限元的数量。
   * 由于顶点没有存储计算所需的信息，这个方法只是引发一个异常，只是为了实现与尺寸无关的编程而存在。
   *
   */
  unsigned int
  n_active_fe_indices() const;

  /**
   * 返回此对象上的 @p n-th 活动FE索引。
   * 由于顶点没有存储计算所需的信息，这个方法只是引发一个异常，并且只是为了实现独立于维度的编程。
   *
   */
  unsigned int
  nth_active_fe_index(const unsigned int n) const;

  /**
   * 如果给定索引的有限元在当前对象上是活动的，则返回真。
   * 由于顶点没有存储计算所需的信息，这个方法只是引发一个异常，并且只存在于实现独立于维度的编程。
   *
   */
  bool
  fe_index_is_active(const unsigned int fe_index) const;

  /**
   * 返回给定 @p fe_index. 用于此对象的有限元的引用  @p
   * fe_index 必须用于此对象，即
   * <code>fe_index_is_active(fe_index)</code> 必须返回true。
   *
   */
  const FiniteElement<1, spacedim> &
  get_fe(const unsigned int fe_index) const;

  /**
   * @}
   *
   */

  /**
   * 子类的例外情况
   * @ingroup Exceptions
   *
   */
  DeclExceptionMsg(ExcInvalidObject,
                   "This accessor object has not been "
                   "associated with any DoFHandler object.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcVectorNotEmpty);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcVectorDoesNotMatch);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcMatrixDoesNotMatch);
  /**
   * 为一个应该是 @ref GlossActive "活动 "
   * 的单元格调用了一个函数，但它被细化了。
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcNotActive);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcCantCompareIterators);

protected:
  /**
   * 存储要访问的DoFHandler对象的地址。
   *
   */
  DoFHandler<1, spacedim> *dof_handler;

  /**
   * 比较是否相等。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator==(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * 比较不等式。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  bool
  operator!=(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &) const;

  /**
   * 重置DoF处理程序指针。
   *
   */
  void set_dof_handler(DoFHandler<1, spacedim> *dh);

  /**
   * 将此对象的<i>i</i>个自由度的索引设置为 @p 个索引。
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的同一个有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果相匹配。
   *
   */
  void
  set_dof_index(
    const unsigned int            i,
    const types::global_dof_index index,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  /**
   * 将当前单元的 @p vertex-th
   * 顶点上的<i>i</i>度的全局索引设置为 @p index.
   * 最后一个参数表示有限元索引。对于标准的::DoFHandler类，这个值必须等于其默认值，因为该类无论如何只支持所有单元上的相同有限元。
   * 然而，当相关的DoFHandler具有hp-capabilities时，不同的有限元对象可以在不同的单元上使用。在两个单元格之间的面上，以及顶点上，可能会有两组自由度，相邻单元格上使用的有限元各有一组。
   * 为了指定在哪一组自由度上工作，最后一个参数被用来区别对待。最后，如果这个函数是为一个单元对象调用的，那么只能有一个自由度集，而且fe_index必须与active_fe_index()的结果一致。
   *
   */
  void
  set_vertex_dof_index(
    const unsigned int            vertex,
    const unsigned int            i,
    const types::global_dof_index index,
    const unsigned int fe_index = AccessorData::invalid_fe_index) const;

  // Iterator classes need to be friends because they need to access
  // operator== and operator!=.
  template <typename>
  friend class TriaRawIterator;


  // Make the DoFHandler class a friend so that it can call the set_xxx()
  // functions.
  template <int, int>
  friend class DoFHandler;

  friend struct dealii::internal::DoFHandlerImplementation::Policy::
    Implementation;
  friend struct dealii::internal::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::hp::DoFHandlerImplementation::Implementation;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
};



 /* -------------------------------------------------------------------------- */ 


/**
 * 一个代表DoF访问器对象的类，用于表示不合理的迭代器，如1D网格上的四维迭代器。
 * 这个类不能用来创建对象（事实上，如果试图这样做的话，它会抛出一个异常，但它有时允许以独立于维度的方式，以更简单的方式编写代码。例如，它允许编写独立于维度的四边形迭代器的代码
 *
 * - 即，也可以在1d中进行编译
 *
 * - 因为四元迭代器（通过当前的类）存在，并且在语法上是正确的。然而，你不能期望在1d中创建这些迭代器中的一个实际对象，这意味着你需要期望将使用四元迭代器的代码块包装成类似 <code>if (dim@>1)</code> 的东西。
 *
 * - 反正这也是很有意义的。
 * 这个类提供了Accessor类与Iterator类交互所需的最小接口。然而，这只是为了语法上的正确性，这些函数除了产生错误之外，没有任何作用。
 *
 *
 * @ingroup Accessors
 *
 */
template <int structdim, int dim, int spacedim = dim>
class DoFInvalidAccessor : public InvalidAccessor<structdim, dim, spacedim>
{
public:
  /**
   * 从基类传播别名到这个类。
   *
   */
  using AccessorData =
    typename InvalidAccessor<structdim, dim, spacedim>::AccessorData;

  /**
   * 构造器。
   * 这个类用于在给定维度中没有意义的迭代器，例如1D网格的四边形。因此，虽然这种对象的创建在语法上是有效的，但它们在语义上没有意义，当这种对象实际生成时，我们会产生一个异常。
   *
   */
  DoFInvalidAccessor(const Triangulation<dim, spacedim> *parent     = 0,
                     const int                           level      = -1,
                     const int                           index      = -1,
                     const AccessorData *                local_data = 0);

  /**
   * 复制构造函数。
   * 这个类用于在给定维度中没有意义的迭代器，例如1D网格的四边形。因此，虽然这种对象的创建在语法上是有效的，但它们在语义上没有意义，当这种对象实际生成时，我们会产生一个异常。
   *
   */
  DoFInvalidAccessor(const DoFInvalidAccessor<structdim, dim, spacedim> &);

  /**
   * 从其他访问器转换到当前无效的访问器。这当然也会导致运行时错误。
   *
   */
  template <typename OtherAccessor>
  DoFInvalidAccessor(const OtherAccessor &);

  /**
   * 返回这个对象的<i>i</i>个自由度的索引到 @p index.
   * ，因为当前对象没有指向任何有用的东西，像这个类中的所有其他函数一样，这个函数只抛出一个异常。
   *
   */
  types::global_dof_index
  dof_index(const unsigned int i,
            const unsigned int fe_index =
              DoFHandler<dim, spacedim>::default_fe_index) const;

  /**
   * 将此对象的<i>i</i>个自由度的索引设置为 @p
   * 个索引。由于当前对象没有指向任何有用的东西，像这个类中的所有其他函数一样，这个函数只抛出一个异常。
   *
   */
  void
  set_dof_index(const unsigned int            i,
                const types::global_dof_index index,
                const unsigned int            fe_index =
                  DoFHandler<dim, spacedim>::invalid_fe_index) const;
};



 /* -------------------------------------------------------------------------- */ 


/**
 * 授予对单元格的自由度的访问权。
 * 注意，因为对于我们派生的类，即<tt>DoFAccessor<dim></tt>，两个模板参数是相等的，基类实际上是派生自CellAccessor，这使得这个类的函数对DoFCellAccessor类也可用。
 *
 *
 * @ingroup dofs
 *
 * @ingroup Accessors
 *
 */
template <int dimension_, int space_dimension_, bool level_dof_access>
class DoFCellAccessor : public DoFAccessor<dimension_,
                                           dimension_,
                                           space_dimension_,
                                           level_dof_access>
{
public:
  /**
   * 从DoFHandler中提取尺寸。
   *
   */
  static const unsigned int dim = dimension_;

  /**
   * 从DoFHandler中提取空间维度。
   *
   */
  static const unsigned int spacedim = space_dimension_;


  /**
   * 由迭代器类传递的数据类型。
   *
   */
  using AccessorData = DoFHandler<dimension_, space_dimension_>;

  /**
   * 声明基类的别名，使访问一些异常类更简单。
   *
   */
  using BaseClass =
    DoFAccessor<dimension_, dimension_, space_dimension_, level_dof_access>;

  /**
   * 定义这个容器的类型，是它的一部分。
   *
   */
  using Container = DoFHandler<dimension_, space_dimension_>;

  /**
   * 一个单元格的面的迭代器的类型。这就是face()函数的返回值。
   *
   */
  using face_iterator = TriaIterator<DoFAccessor<dimension_ - 1,
                                                 dimension_,
                                                 space_dimension_,
                                                 level_dof_access>>;

  /**
   * @name  构造器和初始化
   *
   */
  /**
   * @{
   *
   */

  /**
   * 构造函数
   *
   */
  DoFCellAccessor(const Triangulation<dimension_, space_dimension_> *tria,
                  const int                                          level,
                  const int                                          index,
                  const AccessorData *local_data);

  /**
   * 转换构造器。这个构造器的存在是为了使某些构造在独立于维度的代码中写得更简单。例如，它允许将一个面的迭代器分配给一个线的迭代器，这个操作在2D中很有用，但在3D中没有任何意义。这里的构造函数是为了使代码符合C++的要求而存在的，但它会无条件地中止；换句话说，将一个面迭代器分配给一个线迭代器，最好放在一个if语句中，检查维度是否为2，并在3D中分配给一个四边形迭代器（如果没有这个构造函数，如果我们碰巧为2d编译，这个操作是非法的）。
   *
   */
  template <int structdim2, int dim2, int spacedim2>
  DoFCellAccessor(const InvalidAccessor<structdim2, dim2, spacedim2> &);

  /**
   * 另一个对象之间的转换操作符，就像之前的那个一样，没有意义。
   *
   */
  template <int structdim2, int dim2, int spacedim2, bool level_dof_access2>
  explicit DoFCellAccessor(
    const DoFAccessor<structdim2, dim2, spacedim2, level_dof_access2> &);

  /**
   * 复制构造器。
   *
   */
  DoFCellAccessor(
    const DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &) =
    default;

  /**
   * 移动构造函数。
   *
   */
  DoFCellAccessor(                                                  // NOLINT
    DoFCellAccessor<dimension_, space_dimension_, level_dof_access> // NOLINT
      &&) = default;                                                // NOLINT

  /**
   * 解除构造函数
   *
   */
  ~DoFCellAccessor() = default;

  /**
   * 复制操作符。这些操作符通常在类似<tt>iterator
   * a,b;a=*b;</tt>的情况下使用。据推测，这里的意图是将 @p b
   * 所指向的对象复制到 @p a.
   * 所指向的对象。然而，取消引用迭代器的结果不是一个对象，而是一个访问器；因此，这个操作对于DoF处理程序对象上的迭代器是没用的。
   * 因此，这个操作被声明为删除，不能被使用。
   *
   */
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &
  operator=(
    const DoFCellAccessor<dimension_, space_dimension_, level_dof_access> &da) =
    delete;

  /**
   * 移动赋值运算符。
   *
   */
  DoFCellAccessor<dimension_, space_dimension_, level_dof_access> & // NOLINT
  operator=(                                                        // NOLINT
    DoFCellAccessor<dimension_, space_dimension_, level_dof_access> // NOLINT
      &&) = default;                                                // NOLINT

  /**
   * @}
   *
   */

  /**
   * 以DoF单元格迭代器的形式返回该单元格的父级。如果父对象不存在（即，如果该对象处于网格层次结构的最粗层），将产生一个异常。
   * 这个函数是需要的，因为基类CellAccessor的父函数返回一个没有访问DoF数据的三角形单元访问器。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  parent() const;

  /**
   * @name  访问子对象和相邻对象
   *
   */
  /**
   * @{
   *
   */

  /**
   * 将 @p ith
   * 邻居作为DoF单元迭代器返回。这个函数是需要的，因为基类的邻居函数返回一个没有访问DoF数据的单元格访问器。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor(const unsigned int i) const;

  /**
   * 返回 @p ith
   * 周期性邻居作为DoF单元的迭代器。这个函数是需要的，因为基类的邻居函数返回一个没有访问DoF数据的单元访问器。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  periodic_neighbor(const unsigned int i) const;

  /**
   * 返回 @p ith 邻居或周期性邻居作为DoF单元的迭代器。
   * 这个函数是需要的，因为基类的邻居函数返回一个没有访问DoF数据的单元格访问器。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor_or_periodic_neighbor(const unsigned int i) const;

  /**
   * 将 @p ith
   * 的子单元作为DoF单元迭代器返回。这个函数是需要的，因为基类的子函数返回一个没有访问DoF数据的单元格访问器。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  child(const unsigned int i) const;

  /**
   * 返回该单元格所有子节点的迭代器数组。
   *
   */
  boost::container::small_vector<
    TriaIterator<
      DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>,
    GeometryInfo<dimension_>::max_children_per_cell>
  child_iterators() const;

  /**
   * 返回此单元格的 @p ith 面的一个迭代器。
   * 这个函数返回一个一维的 <code>structdim == 0</code>
   * 的DoFAccessor，二维的 DoFAccessor::line ，以及三维的
   * DoFAccessor::quad 。
   *
   */
  face_iterator
  face(const unsigned int i) const;

  /**
   * 返回该单元格所有面的迭代器数组。
   *
   */
  boost::container::small_vector<face_iterator,
                                 GeometryInfo<dimension_>::faces_per_cell>
  face_iterators() const;

  /**
   * 返回基类中 @p neighbor_child_on_subface
   * 函数的结果，但将其转换为也可以访问DoF数据（基类中的函数只返回一个访问三角形数据的迭代器）。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  neighbor_child_on_subface(const unsigned int face_no,
                            const unsigned int subface_no) const;

  /**
   * 返回基类中 @p periodic_neighbor_child_on_subface
   * 函数的结果，但将其转换为也可以访问DoF数据（基类中的函数只返回一个可以访问三角测量数据的迭代器）。
   *
   */
  TriaIterator<DoFCellAccessor<dimension_, space_dimension_, level_dof_access>>
  periodic_neighbor_child_on_subface(const unsigned int face_no,
                                     const unsigned int subface_no) const;

  /**
   * @}
   *
   */

  /**
   * @name  从全局向量中提取数值
   *
   */
  /**
   * @{
   *
   */

  /**
   * 收集限制在此单元格的道夫上的给定向量的值，其标准排序为：顶点0上的道夫，顶点1上的道夫，等等，行0上的道夫，行1上的道夫，等等，四边0上的道夫，等等。换句话说，这个函数实现了一个[聚集操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在传递给这个函数之前，向量必须有合适的大小。这个函数只适用于活动单元的调用。
   * 输入的向量可以是<tt>Vector<float></tt>、Vector<double>或BlockVector<double>，或者是PETSc或Trilinos向量，如果deal.II被编译为支持这些库。调用者有责任保证存储在输入和输出向量中的数字类型是兼容的，并且具有相似的精度。
   *
   */
  template <class InputVector, typename number>
  void
  get_dof_values(const InputVector &values, Vector<number> &local_values) const;

  /**
   * 收集限制在该单元格的道夫上的给定向量的值，其标准排序为：顶点0的道夫，顶点1的道夫，等等，行0的道夫，行1的道夫，等等，四边0的道夫，等等。换句话说，这个函数实现了一个[聚集操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在传递给这个函数之前，向量必须有合适的大小。这个函数只适用于活动单元的调用。
   * 输入的向量可以是<tt>Vector<float></tt>、Vector<double>或BlockVector<double>，或者是PETSc或Trilinos向量，如果deal.II被编译为支持这些库。调用者有责任保证存储在输入和输出向量中的数字类型是兼容的，并且具有相似的精度。
   *
   */
  template <class InputVector, typename ForwardIterator>
  void
  get_dof_values(const InputVector &values,
                 ForwardIterator    local_values_begin,
                 ForwardIterator    local_values_end) const;

  /**
   * 收集限制在该单元格的道夫上的给定向量的值，其标准排序为：顶点0的道夫，顶点1的道夫，等等，行0的道夫，行1的道夫，等等，四边0的道夫，等等。换句话说，这个函数实现了一个[聚集操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在传递给这个函数之前，向量必须有合适的大小。这个函数只适用于活动单元的调用。
   * 输入的向量可以是<tt>Vector<float></tt>、Vector<double>或BlockVector<double>，或者是PETSc或Trilinos向量，如果deal.II被编译为支持这些库。调用者有责任保证存储在输入和输出向量中的数字类型是兼容的，并且具有相似的精度。作为参数传递给该函数的AffineConstraints对象可以确保在计算dof值时约束条件的分布是正确的。
   *
   */
  template <class InputVector, typename ForwardIterator>
  void
  get_dof_values(
    const AffineConstraints<typename InputVector::value_type> &constraints,
    const InputVector &                                        values,
    ForwardIterator local_values_begin,
    ForwardIterator local_values_end) const;

  /**
   * 这个函数是get_dof_values()的对应函数：它接收这个迭代器所指向的单元的自由度值的向量，并将这些值写入全局数据向量
   * @p values.
   * 换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 这个函数只对活动单元可调用。
   * 请注意，对于连续的有限元，调用这个函数也会影响到相邻单元的dof值。如果相邻单元的精细度低于当前单元，它还可能违反悬挂节点的连续性要求。这些要求没有被考虑到，必须由用户事后强制执行。
   * 向量在被传递给这个函数之前必须有合适的大小。
   * 输出的向量可以是Vector<float>, Vector<double>,
   * 或BlockVector<double>,
   * 或PETSc向量，如果deal.II被编译为支持这些库。调用者有责任保证存储在输入和输出向量中的数字类型是兼容的，并且具有相似的精度。
   *
   */
  template <class OutputVector, typename number>
  void
  set_dof_values(const Vector<number> &local_values,
                 OutputVector &        values) const;

  /**
   * 返回给定的有限元函数对当前单元的插值。在最简单的情况下，该单元是一个终端单元，即它没有子单元；那么，返回值是该单元上的节点值向量。你也可以通过
   * @p
   * get_dof_values函数获得所需的值。在另一种情况下，当单元有子节点时，我们使用有限元类提供的限制矩阵来计算从子节点到本单元的内插。
   * 如果单元是具有hp能力的DoFHandler的一部分，单元只有在活动时才有相关的有限元空间。然而，这个函数也应该提供有子单元的非活动单元的信息。因此，它带有第三个参数，可以在hp-context中使用，表示我们应该插值到的有限元空间。如果单元是活动的，这个函数就会从这个单元的
   * <code>values</code> 向量中获得有限元函数，并将其插值到
   * <code>fe_index</code>
   * 中的第1个元素所描述的空间，这个单元是DoFHandler中的一部分。如果该单元不是活动的，那么我们首先对其所有的终端子单元进行插值，然后将这个函数插值到所要求的单元，保持函数空间不变。
   * 假设两个输入向量事先已经有了合适的大小。
   * @note
   * 与get_dof_values()函数不同，这个函数只对单元格有效，而不是对线、四边形和六边形，因为插值目前只由有限元类为单元格提供。
   *
   */
  template <class InputVector, typename number>
  void
  get_interpolated_dof_values(
    const InputVector &values,
    Vector<number> &   interpolated_values,
    const unsigned int fe_index =
      DoFHandler<dimension_, space_dimension_>::invalid_fe_index) const;

  /**
   * 这个函数是get_interpolated_dof_values()的对应函数：你指定单元格上的dof值，这些值被内插到当前单元格的子单元，并设置在终端单元上。
   * 原则上，它的工作原理如下：如果这个对象指向的单元格是终端（即没有子单元），那么通过调用set_dof_values()函数在全局数据向量中设置dof值；否则，这些值被延长到每个子单元，并为每个子单元调用这个函数。
   * 使用get_interpolated_dof_values()和这个函数，你可以计算一个有限元函数对一个更粗的网格的内插，首先在粗网格的一个单元上得到内插的解，之后用这个函数重新分配。
   * 请注意，对于连续有限元，调用这个函数也会影响到相邻单元上的道夫值。如果相邻的单元比现在的单元细化程度低，或者它们的子单元比这个单元的子单元细化程度低，那么它也可能违反悬挂节点的连续性要求。这些要求没有得到照顾，必须由用户事后强制执行。
   * 如果单元格是具有hp能力的DoFHandler的一部分，单元格只有在活动时才有相关的有限元空间。然而，这个函数也应该在有子代的非活动单元上工作。
   * 因此，它带有第三个参数，可以在hp-上下文中使用，表示我们应该解释这个函数的输入矢量的有限元空间。如果单元格是活动的，这个函数就会将输入向量解释为与该单元格所属的DoFHandler相关的
   * <code>fe_index</code>
   * 的第三个元素所描述的空间元素，并将其插值到与该单元格相关的空间。另一方面，如果单元格不是活动的，那么我们首先使用给定的
   * <code>fe_index</code>
   * 从这个单元格到它的子单元格进行内插，直到我们在一个活动的单元格上结束，这时我们遵循本段开头的程序。
   * 假设两个向量事先已经有了合适的大小。
   * 这个函数依赖于一个单元的有限元空间对其子女的自然插值属性的存在，由有限元类的延长矩阵表示。对于某些元素，粗网格和细网格上的空间没有嵌套，在这种情况下，对子单元的插值不是相同的；关于在这种情况下延长矩阵所代表的内容，请参考各自的有限元类的文档。
   * @note
   * 与get_dof_values()函数不同，这个函数只对单元格有效，而不是对线、四边形和六边形，因为插值目前只由有限元类对单元格提供。
   *
   */
  template <class OutputVector, typename number>
  void
  set_dof_values_by_interpolation(
    const Vector<number> &local_values,
    OutputVector &        values,
    const unsigned int    fe_index =
      DoFHandler<dimension_, space_dimension_>::invalid_fe_index) const;

  /**
   * 通过将自由度的局部编号映射到全局编号并将局部值输入全局向量，将局部（基于单元）向量分配到全局向量。换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 这些元素被 <em>
   * 添加到全局向量的现有元素中，而不是直接设置，因为这通常是人们想要的。如果你需要处理约束条件，你可能还想看一下
   * AffineConstraints::distribute_local_to_global() 函数。
   *
   */
  template <typename number, typename OutputVector>
  void
  distribute_local_to_global(const Vector<number> &local_source,
                             OutputVector &        global_destination) const;

  /**
   * 通过将自由度的本地编号映射到全局编号并将本地值输入全局向量，将迭代器格式的本地（基于单元）向量分配给全局向量。
   * 换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 这些元素被 <em>
   * 添加到全局向量中的现有元素上，而不是直接设置，因为这通常是人们想要的。如果你需要处理约束条件，你可能还想看一下
   * AffineConstraints::distribute_local_to_global() 函数。
   *
   */
  template <typename ForwardIterator, typename OutputVector>
  void
  distribute_local_to_global(ForwardIterator local_source_begin,
                             ForwardIterator local_source_end,
                             OutputVector &  global_destination) const;

  /**
   * 通过将自由度的本地编号映射到全局编号并将本地值输入全局向量，将迭代器格式的本地（基于单元）向量分布到全局向量中。
   * 换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 这些元素被 <em> 添加到 </em>
   * 全局向量中的元素，而不是仅仅设置，因为这通常是人们想要的。此外，传递给这个函数的AffineConstraints对象确保了在这个过程中也消除了约束。
   *
   */
  template <typename ForwardIterator, typename OutputVector>
  void
  distribute_local_to_global(
    const AffineConstraints<typename OutputVector::value_type> &constraints,
    ForwardIterator local_source_begin,
    ForwardIterator local_source_end,
    OutputVector &  global_destination) const;

  /**
   * 这个函数与<tt>distribute_local_to_global(Vector,Vector)</tt>函数的作用基本相同，但是对矩阵而不是向量进行操作。如果矩阵类型是稀疏矩阵，那么它就应该在需要的地方有非零的入口槽。
   *
   */
  template <typename number, typename OutputMatrix>
  void
  distribute_local_to_global(const FullMatrix<number> &local_source,
                             OutputMatrix &global_destination) const;

  /**
   * 这个函数做了两个<tt>distribute_local_to_global</tt>函数与向量和矩阵参数的作用，但都是一次性的。
   *
   */
  template <typename number, typename OutputMatrix, typename OutputVector>
  void
  distribute_local_to_global(const FullMatrix<number> &local_matrix,
                             const Vector<number> &    local_vector,
                             OutputMatrix &            global_matrix,
                             OutputVector &            global_vector) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问此对象的DoF指数
   *
   */

  /**
   * @{
   *
   */

  /**
   * 获取此单元上的局部自由度的全局指数。
   * 如果这个对象访问了一个水平单元（由第三个模板参数或#is_level_cell表示），那么返回get_mg_dof_indices()的结果，否则返回get_dof_indices()。
   * 当调用begin_mg()时，你会得到一个level_cell_iterator，否则就是一个普通的level_cell_iterator。
   * 这种用法的例子是在DoFRenumbering的实现中。
   *
   */
  void
  get_active_or_mg_dof_indices(
    std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * 返回位于此对象上的自由度的<i>global</i>指数，其标准排序由有限元定义（即顶点0上的自由度，顶点1上的自由度，等等，行0上的自由度，行1上的自由度，等等，四边形0上的自由度，等等）该函数仅适用于<i>active</i>对象（见 @ref GlossActive  "本词汇条"
   * ）。      @param[out]  dof_indices
   * 指数将被写入的向量。在传递给这个函数之前，它必须有合适的大小（即
   * <code>fe.n_dofs_per_cell()</code>, <code>fe.dofs_per_face</code>  ，或
   * <code>fe.dofs_per_line</code>
   * ，取决于这个函数被调用的对象的种类）。
   * 这个函数重新实现了基类中的同一个函数。与基类中的函数相比，我们在这里不需要
   * <code>fe_index</code>
   * ，因为单元格上总是有一个唯一的有限元索引。
   * 这是一个要求单元格处于活动状态的函数。
   * 也请看get_active_or_mg_dof_indices()。
   * @note
   * 在本教程的许多地方和库的其他地方，这个函数的参数按惯例被称为
   * <code>local_dof_indices</code>
   * 。这个名字并不意味着表示自由度的<i>local</i>数量（总是在0和
   * <code>fe.n_dofs_per_cell()</code>
   * 之间），而是表示返回值是位于当前单元格上的那些自由度的<i>global</i>指数。
   *
   */
  void
  get_dof_indices(std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * 检索该单元上的自由度在与该单元的级别相关的级别向量中的全局指数。
   *
   */
  void
  get_mg_dof_indices(std::vector<types::global_dof_index> &dof_indices) const;

  /**
   * @}
   *
   */

  /**
   * @name  访问与此对象相关的有限元
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回用于此迭代器所指向的单元格的有限元。对于没有hp-capabilities的DoFHandler对象，这当然总是同一个元素，与我们目前所在的单元无关，但是对于hp-DoFHandler对象，这可能会在不同的单元之间变化。
   * @note
   * 由于自由度只存在于具有hp-capabilities的DoFHandler对象的活动单元上（即目前没有实现多层次的此类对象），在非活动单元上查询有限元是没有意义的，因为它们没有任何自由度的有限元空间与它们相关联。因此，这个函数在非活动单元上调用时将产生一个异常。
   *
   */
  const FiniteElement<dimension_, space_dimension_> &
  get_fe() const;

  /**
   * 返回用于该单元的有限元空间 hp::FECollection
   * 内的索引。这个函数只有在与当前单元格相关的DoFHandler对象启用了hp-capabilities时才有用。
   * @note
   * 由于自由度只存在于具有hp-capabilities的DoFHandler对象的活动单元上（即目前没有实现多层次的此类对象），在非活动单元上查询有限元是没有意义的，因为它们没有与之相关的有限元空间，没有任何自由度。因此，当在非活动单元上调用该函数时将产生一个异常。
   * @note  当使用平行网格时，无论是通过
   * parallel::shared::Triangulation 还是
   * parallel::distributed::Triangulation
   * 类，都只允许在本地拥有的或幽灵单元上调用此函数。没有关于人工单元的信息。
   * 此外， @p active_fe_index 信息只在调用 DoFHandler::set_fe() 和
   * DoFHandler::distribute_dofs().
   * 期间从一个处理器上的本地拥有的细胞交换到其他可能是幽灵细胞的处理器。
   * 请注意，如果你在调用这些函数之一后在一个细胞上调用set_active_fe_index()，那么这个信息将不会传播给其他可能将这个细胞作为幽灵细胞的处理器。更多信息请参见DoFHandler的文档。
   *
   */
  unsigned int
  active_fe_index() const;

  /**
   * 设置用于此单元的FiniteElement的索引。这决定了要使用
   * hp::FECollection
   * 中的哪个元素。这个函数只有在与当前单元相关的DoF处理程序对象启用了hp-功能时才有用。
   * @note
   * 由于自由度只存在于具有hp-能力的DoFHandler对象的活动单元上（即目前没有实现多层次的这种对象），在非活动单元上查询有限元是没有意义的，因为它们没有与之相关的有限元空间，没有任何自由度。因此，当在非活动单元上调用该函数时，将产生一个异常。
   * @note  当使用平行网格时，无论是通过
   * parallel::shared::Triangulation 还是
   * parallel::distributed::Triangulation
   * 类，都只允许在本地拥有的或幽灵单元上调用该函数。没有关于人工单元的信息。
   * 此外， @p active_fe_index 信息只在调用 DoFHandler::set_fe() 和
   * DoFHandler::distribute_dofs().
   * 期间，从一个处理器上的本地拥有的细胞交换到其他可能是幽灵细胞的处理器上。
   * 请注意，如果你在调用这些函数之一后，在一个细胞上调用set_active_fe_index()，那么这个信息将不会传播给其他可能将这个细胞作为幽灵细胞的处理器。更多信息请参见DoFHandler的文档。
   *
   */
  void
  set_active_fe_index(const unsigned int i) const;
  /**
   * @}
   *
   */

  /**
   * 将此单元的DoF指数设置为给定值。如果给定的DoF处理程序类存在DoF缓存，该函数会绕过DoF缓存。
   *
   */
  void
  set_dof_indices(const std::vector<types::global_dof_index> &dof_indices);

  /**
   * 将此单元格的Level DoF指数设置为给定值。
   *
   */
  void
  set_mg_dof_indices(const std::vector<types::global_dof_index> &dof_indices);

  /**
   * 更新缓存，我们在其中存储该单元的DoF指数。
   *
   */
  void
  update_cell_dof_indices_cache() const;

  /**
   * @name 处理细化指标
   *
   */
  /**
   * @{
   *
   */

  /**
   * 返回下次三角剖分被细化和粗化时将分配给该单元的有限元。如果没有通过set_future_fe_index()函数为该单元指定未来的有限元，则活动的有限元将保持不变，在这种情况下，将返回活动的有限元。
   * 对于没有启用hp-capabilities的DoFHandler，这当然总是同一个元素，与我们目前所在的单元无关，但是对于hp-
   * DoFHandler对象，这可能会在不同的单元之间变化。
   * @note
   * 由于自由度只存在于具有hp-capabilities的DoFHandler对象的活动单元上（即目前没有实现多层次的此类对象），在非活动单元上查询有限元是没有意义的，因为它们没有与之相关的有限元空间，没有任何自由度。因此，这个函数在非活动单元上调用时将产生一个异常。
   *
   */
  const FiniteElement<dimension_, space_dimension_> &
  get_future_fe() const;

  /**
   * 返回有限元的fe_index，该有限元将在下一次三角结构被细化和粗化时分配给该单元。如果没有通过set_future_fe_index()函数为该单元指定未来的有限元，则活动的单元将保持不变，在这种情况下，将返回活动有限元的fe_index。
   * @note
   * 由于自由度只存在于具有hp-capabilities的DoFHandler对象的活动单元上（即目前没有实现多层次的这种对象），在非活动单元上查询有限元是没有意义的，因为它们没有与之相关的有限元空间，没有任何自由度。因此，当在非活动单元上调用该函数时，将产生一个异常。
   * @note 当使用平行网格时，无论是通过
   * parallel::shared::Triangulation 还是
   * parallel::distributed::Triangulation
   * 类，只允许在本地拥有的单元上调用该函数。
   *
   */
  unsigned int
  future_fe_index() const;

  /**
   * 设置有限元的fe_index，该有限元将在下一次三角结构被细化和粗化时被分配给该单元。之前分配的未来有限元将被覆盖。
   * 参见future_fe_index()的注释，以了解对该功能的限制信息。
   *
   */
  void
  set_future_fe_index(const unsigned int i) const;

  /**
   * 返回是否已经设置了未来有限元。
   * 参见future_fe_index()的注释，以了解对该功能的限制信息。
   *
   */
  bool
  future_fe_index_set() const;

  /**
   * 撤销分配的未来有限元。因此，在下一次三角结构被细化和粗化时，活动的有限元将保持不变。
   * 参见future_fe_index()的说明，以了解对该功能的限制。
   *
   */
  void
  clear_future_fe_index() const;
  /**
   * @}
   *
   */

private:
  // Make the DoFHandler class a friend so that it can call the
  // update_cell_dof_indices_cache() function
  template <int, int>
  friend class DoFHandler;
  friend struct dealii::internal::DoFCellAccessorImplementation::Implementation;
};


template <int structdim, int dim, int spacedim, bool level_dof_access>
inline bool
DoFAccessor<structdim, dim, spacedim, level_dof_access>::is_level_cell()
{
  return level_dof_access;
}



template <int structdim, int dim, int spacedim>
template <typename OtherAccessor>
DoFInvalidAccessor<structdim, dim, spacedim>::DoFInvalidAccessor(
  const OtherAccessor &)
{
  Assert(false,
         ExcMessage("You are attempting an illegal conversion between "
                    "iterator/accessor types. The constructor you call "
                    "only exists to make certain template constructs "
                    "easier to write as dimension independent code but "
                    "the conversion is not valid in the current context."));
}



DEAL_II_NAMESPACE_CLOSE

// include more templates
#include "dof_accessor.templates.h"


#endif


