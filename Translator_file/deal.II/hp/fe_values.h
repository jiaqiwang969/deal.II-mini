//include/deal.II-translator/hp/fe_values_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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

#ifndef dealii_hp_fe_values_h
#define dealii_hp_fe_values_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <map>
#include <memory>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class FiniteElement;
#endif


namespace hp
{
  /**
   * <tt>hp::FE*Values</tt>
   * 类的基类，存储它们的共同数据。这个类的主要任务是提供一个表格，对于来自其相应集合对象的每一个有限元、映射和正交对象的组合，都有一个匹配的::FEValues,
   * ::FEFaceValues, 或::FESubfaceValues对象。
   * 然而，为了使事情更有效率，这些FE*Values对象只在请求时才被创建（懒惰分配）。如果需要的话，也可以通过相应的precalculate_fe_values()函数提前计算所有对象来绕过这一点。
   * 第一个模板参数表示我们所处的空间维度，第二个是我们所整合的对象的维度，即对于通常的
   * @p hp::FEValues
   * 来说，它等于第一个维度，而对于面部整合来说，它要少一个。第三个模板参数表示底层非hp-FE*Values的基础类型，即它可以是::FEValues，:FEFaceValues，或::FESubfaceValues。
   * @ingroup hp
   *
   */
  template <int dim, int q_dim, class FEValuesType>
  class FEValuesBase
  {
  public:
    /**
     * 构造函数。将这个类的字段设置为构造函数的参数所指示的值。
     *
     */
    FEValuesBase(
      const MappingCollection<dim, FEValuesType::space_dimension>
        &mapping_collection,
      const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
      const QCollection<q_dim> &                              q_collection,
      const UpdateFlags                                       update_flags);

    /**
     * 像上面的函数一样，但取一个正交集合的向量。
     * 对于 hp::FEFaceValues,
     * ，正交集合的第1项被解释为应用于第1个面的正交规则。
     *
     */
    FEValuesBase(
      const MappingCollection<dim, FEValuesType::space_dimension>
        &mapping_collection,
      const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
      const std::vector<QCollection<q_dim>> &                 q_collection,
      const UpdateFlags                                       update_flags);

    /**
     * 构造函数。这个构造函数与另一个构造函数相当，只是它使对象隐含地使用
     * $Q_1$ 映射（即MappingQGeneric（1）类型的对象）。
     *
     */
    FEValuesBase(
      const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
      const QCollection<q_dim> &                              q_collection,
      const UpdateFlags                                       update_flags);

    /**
     * 像上面的函数一样，但取一个矢量正交集合。    对于
     * hp::FEFaceValues,
     * ，正交集合的第i个条目被解释为应用于第i个面的面的正交规则。
     *
     */
    FEValuesBase(
      const FECollection<dim, FEValuesType::space_dimension> &fe_collection,
      const std::vector<QCollection<q_dim>> &                 q_collection,
      const UpdateFlags                                       update_flags);

    /**
     * 复制构造器。
     *
     */
    FEValuesBase(const FEValuesBase<dim, q_dim, FEValuesType> &other);

    /**
     * 复制操作符。虽然这种类型的对象可以被复制构建，但它们不能被复制，因此这个运算符被禁用。
     *
     */
    FEValuesBase &
    operator=(const FEValuesBase &) = delete;

    /**
     * 出于时间上的考虑，提前创建所有需要的FE*Values对象可能是有用的，而不是像本类中通常那样通过懒惰分配来计算它们。
     * 这个函数预先计算了与所提供的参数相对应的FE*Values对象。与同一索引相对应的所有向量条目的总和描述了一个FE*Values对象，与select_fe_values()类似。
     *
     */
    void
    precalculate_fe_values(const std::vector<unsigned int> &fe_indices,
                           const std::vector<unsigned int> &mapping_indices,
                           const std::vector<unsigned int> &q_indices);

    /**
     * 同上，面向最常用的 hp::FEValues
     * 对象，其中FE、正交和映射指数在每个单独的单元上都是相似的。
     * 为FECollection中的每个FE创建FE*Values对象，正交和映射分别对应于QuadratureCollection和MappingCollection中的相同索引。
     * 如果QuadratureCollection或MappingCollection只包含一个对象，它将用于所有FE*Values对象。
     *
     */
    void
    precalculate_fe_values();

    /**
     * 获取此处使用的有限元对象集合的引用。
     *
     */
    const FECollection<dim, FEValuesType::space_dimension> &
    get_fe_collection() const;

    /**
     * 获取此处使用的映射对象集合的引用。
     *
     */
    const MappingCollection<dim, FEValuesType::space_dimension> &
    get_mapping_collection() const;

    /**
     * 获取此处使用的正交对象集合的引用。
     *
     */
    const QCollection<q_dim> &
    get_quadrature_collection() const;

    /**
     * 获取底层更新标志。
     *
     */
    UpdateFlags
    get_update_flags() const;

    /**
     * 返回最后一次调用select_fe_values()所选择的 @p FEValues
     * 对象的引用。当你最后一次调用 <tt>hp::FE*Values</tt>
     * 类的 @p reinit 函数时，依次调用select_fe_values()。
     *
     */
    const FEValuesType &
    get_present_fe_values() const;

  protected:
    /**
     * 选择一个适合给定FE、正交和映射指数的FEValues对象。如果这样的对象还不存在，就创建一个。
     * 该函数返回一个可写的引用，这样派生类也可以重新引用()所选择的FEValues对象。
     *
     */
    FEValuesType &
    select_fe_values(const unsigned int fe_index,
                     const unsigned int mapping_index,
                     const unsigned int q_index);

  protected:
    /**
     * 一个指向要使用的有限元集合的指针。
     *
     */
    const SmartPointer<const FECollection<dim, FEValuesType::space_dimension>,
                       FEValuesBase<dim, q_dim, FEValuesType>>
      fe_collection;

    /**
     * 一个指向要使用的映射集合的指针。
     *
     */
    const SmartPointer<
      const MappingCollection<dim, FEValuesType::space_dimension>,
      FEValuesBase<dim, q_dim, FEValuesType>>
      mapping_collection;

    /**
     * 提供给构造函数的正交集合对象的副本。
     *
     */
    const QCollection<q_dim> q_collection;

    /**
     * 正交集合的矢量。对于 hp::FEFaceValues,
     * ，正交集的第1条被解释为应用于第1个面的正交规则。
     * 变量q_collection收集了每个正交集合的第一个正交规则。
     *
     */
    const std::vector<QCollection<q_dim>> q_collections;

  private:
    /**
     * 一个表格，我们在其中存储指向不同的有限元、映射和正交对象集合的fe_values对象的指针。
     * 第一个索引表示fe_collection中有限元的索引，第二个索引表示mapping集合中映射的索引，最后一个索引表示q_collection中正交公式的索引。
     * 最初，所有条目都是零指针，我们将在select_fe_values()或precalculate_fe_values()中根据需要懒散地分配它们。
     *
     */
    Table<3, std::unique_ptr<FEValuesType>> fe_values_table;

    /**
     * 指向上次调用select_fe_value()函数时选择的fe_values对象的一组索引。
     *
     */
    TableIndices<3> present_fe_values_index;

    /**
     * 给予构造函数的更新标志的值。
     *
     */
    const UpdateFlags update_flags;
  };

} // namespace hp


namespace hp
{
  /**
   * 相当于::FEValues类的一个hp值。参见 step-27
   * 教程程序中的使用实例。
   * 这个类的想法如下：当人们在hp-finite
   * element方法中装配矩阵时，不同的单元上可能有不同的有限元，因此人们也可能希望对不同的单元使用不同的正交公式。另一方面，::FEValues有效地处理了对单一有限元和正交对象所需的任何信息的预评估。这个类将这些概念结合起来：它提供了一个::FEValues对象的
   * "集合"。
   * 在构造时，人们传递的不是一个有限元和正交对象（以及可能的映射），而是整个类型
   * hp::FECollection 和 hp::QCollection. 的集合。
   * 后来，当人们坐在一个具体的单元上时，就会为这个特定的单元调用
   * reinit() 函数，就像对普通的 ::FEValues
   * 对象那样。不同的是，这一次，reinit()函数会查找该单元的活动FE索引，如果有必要的话，会创建一个::FEValues对象，在其集合中匹配具有该特定索引的有限元和正交公式，然后为当前单元重新初始化。然后可以使用get_present_fe_values()函数访问适合当前单元格的有限元和正交公式的::FEValues对象，人们可以像对待非hp-DoFHandler对象的任何::FEValues对象一样处理它。
   * reinit()函数有额外的参数，有默认值。如果没有指定，函数会从单元格的活动
   * FE 索引中获取进入  hp::FECollection,   hp::QCollection,  和
   * hp::MappingCollection
   * 对象的索引，如上所述。然而，人们也可以为当前单元格选择不同的索引。例如，通过指定不同的索引进入
   * hp::QCollection
   * 类，就不需要对正交集合中的正交对象进行排序，使其与FE集合中的有限元对象的顺序一一对应（尽管选择这样的顺序肯定很方便）。
   * 请注意::FEValues对象是即时创建的，也就是说，只有在需要时才创建。这确保了我们不会为每一个有限元、正交公式和映射的组合创建对象，而只是创建那些实际需要的对象。如果需要的话，也可以通过相应的
   * hp::FEValuesBase::precalculate_fe_values()
   * 函数提前计算所有对象来绕过这一点。
   * 这个类还没有实现在一维情况下的使用（<tt>spacedim != dim
   * </tt>）。
   * @ingroup hp hpcollection
   *
   */
  template <int dim, int spacedim = dim>
  class FEValues
    : public hp::FEValuesBase<dim, dim, dealii::FEValues<dim, spacedim>>
  {
  public:
    static const unsigned int dimension = dim;

    static const unsigned int space_dimension = spacedim;

    /**
     * 构造函数。用给定的参数初始化这个对象。
     *
     */
    FEValues(const MappingCollection<dim, spacedim> &mapping_collection,
             const FECollection<dim, spacedim> &     fe_collection,
             const QCollection<dim> &                q_collection,
             const UpdateFlags                       update_flags);


    /**
     * 构造函数。这个构造函数等同于另一个构造函数，除了它使对象隐含地使用
     * $Q_1$ 映射（即MappingQGeneric(1)类型的对象）。
     *
     */
    FEValues(const FECollection<dim, spacedim> &fe_collection,
             const QCollection<dim> &           q_collection,
             const UpdateFlags                  update_flags);


    /**
     * 重新初始化给定单元格的对象。
     * 调用后，你可以使用get_present_fe_values()函数得到一个与当前单元格对应的FEValues对象。
     * 对于这个FEValues对象，我们使用下面描述的附加参数来决定使用哪个有限元、映射和正交公式。它们的顺序是这样的：人们可能最想改变的参数排在前面。这些参数的规则如下。
     * 如果这个函数的 @p fe_index
     * 参数保持默认值，那么我们在传递给这个类的构造函数的
     * hp::FECollection 中使用该有限元，其索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出。因此，给这个对象的 hp::FECollection
     * 参数实际上应该与构建与本单元相关的DoFHandler所用的参数相同。另一方面，如果给这个参数一个值，它将覆盖
     * <code>cell-@>active_fe_index()</code>  的选择。        如果 @p
     * q_index
     * 参数被保留为默认值，那么我们使用传递给该类构造函数的
     * hp::QCollection 中的正交公式，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。在这种情况下，
     * hp::FECollection.
     * 中的每个有限元都应该有一个相应的正交公式。
     * 作为一种特殊情况，如果正交集合只包含一个元素（如果想在hp-discretization中对所有有限元使用同一个正交对象，即使这可能不是最有效的），那么这个单一的正交将被使用，除非为这个参数指定一个不同的值。另一方面，如果给这个参数一个值，它将覆盖
     * <code>cell-@>active_fe_index()</code>
     * 的选择或对单一正交的选择。        如果 @p mapping_index
     * 参数保持默认值，那么我们使用传递给该类构造函数的
     * hp::MappingCollection 中的映射对象，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。如上所述，如果映射集合只包含一个元素（如果想对hp-discretization中的所有有限元素使用
     * $Q_1$
     * 映射，这是一个常见的情况），那么这个单一的映射将被使用，除非为这个参数指定一个不同的值。
     *
     */
    template <bool lda>
    void
    reinit(const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);

    /**
     * 像前面的函数一样，但是对于非DoFHandler迭代器。这个函数存在的原因是，人们也可以对Triangulation对象使用
     * hp::FEValues 。        由于 <code>cell-@>active_fe_index()</code>
     * 对三角迭代器没有意义，这个函数从传递给这个对象构造器的相关构造中选择第零个有限元、映射和正交对象。唯一的例外是，如果你为这最后三个参数中的任何一个指定了一个不同于默认值的值。
     *
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);
  };



  /**
   * 这相当于 hp::FEValues
   * 类，但用于脸部整合，也就是说，它对 hp::FEValues
   * 来说就像:FEFaceValues对:FEValues一样。
   * 同样的评论适用于 hp::FEValues
   * 类的文档。然而，需要注意的是，在这里更常见的是，人们希望在reinit()函数中明确指定一个特定正交公式的索引。这是因为默认索引与当前函数上的有限元索引相对应。另一方面，在面的积分通常要用一个正交公式来进行，这个正交公式要根据面的两边使用的有限元来调整。如果我们将
   * hp::FECollection
   * 中的元素按多项式程度升序排列，并将这些有限元素与
   * hp::QCollection 中相应的正交公式相匹配，则 ]
   * 传递给构造函数，那么传递给
   * reinit()函数的正交指数通常应该是  <code>std::max
   * （单元格-  @>active_fe_index(),  邻居-  @>active_fe_index()</code>
   * 以确保选择的正交公式对  <em>  和  </em>
   * 有限元都足够精确。
   * @ingroup hp hpcollection
   *
   */
  template <int dim, int spacedim = dim>
  class FEFaceValues
    : public hp::FEValuesBase<dim, dim - 1, dealii::FEFaceValues<dim, spacedim>>
  {
  public:
    /**
     * 构造函数。用给定的参数初始化这个对象。
     *
     */
    FEFaceValues(const hp::MappingCollection<dim, spacedim> &mapping_collection,
                 const hp::FECollection<dim, spacedim> &     fe_collection,
                 const hp::QCollection<dim - 1> &            q_collection,
                 const UpdateFlags                           update_flags);

    /**
     * 像上面的函数一样，但取一个正交规则集合的向量。这允许为每个面分配一个不同的正交规则：集合的第1个条目被用作第1个面的正交规则。
     * 在集合只包含一个面的正交规则的情况下，这个正交规则被用于所有面。
     *
     */
    FEFaceValues(const hp::MappingCollection<dim, spacedim> &mapping_collection,
                 const hp::FECollection<dim, spacedim> &     fe_collection,
                 const std::vector<hp::QCollection<dim - 1>> &q_collections,
                 const UpdateFlags                            update_flags);


    /**
     * 构造函数。这个构造函数等同于其他的构造函数，只是它使对象隐含地使用
     * $Q_1$ 映射（即MappingQGeneric(1)类型的对象）。
     *
     */
    FEFaceValues(const hp::FECollection<dim, spacedim> &fe_collection,
                 const hp::QCollection<dim - 1> &       q_collection,
                 const UpdateFlags                      update_flags);

    /**
     * 像上面的函数一样，但取一个正交规则集合的向量。这允许为每个面分配一个不同的正交规则：集合的第1个条目被用作第1个面的正交规则。
     * 在集合只包含一个面的正交规则的情况下，这个正交规则被用于所有面。
     *
     */
    FEFaceValues(const hp::FECollection<dim, spacedim> &      fe_collection,
                 const std::vector<hp::QCollection<dim - 1>> &q_collections,
                 const UpdateFlags                            update_flags);

    /**
     * 重新初始化给定单元和面的对象。
     * 调用后，你可以使用get_present_fe_values()函数得到一个与当前单元格对应的FEFaceValues对象。
     * 对于这个FEFaceValues对象，我们使用下面描述的附加参数来决定使用哪个有限元、映射和正交公式。它们的顺序是，人们可能最想改变的参数排在前面。这些参数的规则如下。
     * 如果这个函数的 @p fe_index
     * 参数被保留为默认值，那么我们在传递给这个类的构造函数的
     * hp::FECollection 中使用该有限元，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出。因此，给这个对象的 hp::FECollection
     * 参数实际上应该与构建与本单元相关的DoFHandler所用的参数相同。另一方面，如果给这个参数一个值，它将覆盖
     * <code>cell-@>active_fe_index()</code> 的选择。        如果 @p
     * q_index
     * 参数被保留为默认值，那么我们使用传递给该类构造函数的
     * hp::QCollection 中的正交公式，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。在这种情况下，
     * hp::FECollection.
     * 中的每个有限元都应该有一个相应的正交公式。
     * 作为一种特殊情况，如果正交集合只包含一个元素（如果想在hp-discretization中对所有有限元使用同一个正交对象，即使这可能不是最有效的），那么这个单一的正交被使用，除非为这个参数指定一个不同的值。另一方面，如果给了这个参数一个值，它将覆盖
     * <code>cell-@>active_fe_index()</code>
     * 的选择或对单一正交的选择。        如果 @p mapping_index
     * 参数保持默认值，那么我们使用传递给该类构造函数的
     * hp::MappingCollection 中的映射对象，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。如上所述，如果映射集合只包含一个元素（如果想对hp-discretization中的所有有限元素使用
     * $Q_1$
     * 映射，这是一个常见的情况），那么这个单一的映射将被使用，除非为这个参数指定一个不同的值。
     *
     */
    template <bool lda>
    void
    reinit(const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell,
           const unsigned int                                       face_no,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);

    /**
     * 重新初始化给定单元和面的对象。
     * @note   @p face  必须是 @p cell's 面的迭代器之一。
     *
     */
    template <bool lda>
    void
    reinit(const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &   cell,
           const typename Triangulation<dim, spacedim>::face_iterator &face,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);

    /**
     * 和前面的函数一样，但是对于非DoFHandler迭代器。这个函数存在的原因是，人们也可以对Triangulation对象使用这个类。
     * 由于 <code>cell-@>active_fe_index()</code>
     * 对三角迭代器没有意义，这个函数从传递给这个对象的构造器的相关构造中选择第零个有限元、映射和正交对象。唯一的例外是，如果你为这最后三个参数中的任何一个指定了一个不同于默认值的值。
     *
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const unsigned int                                          face_no,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);

    /**
     * 为给定的单元格和面重新初始化该对象。
     * @note   @p face  必须是 @p cell's 面的迭代器之一。
     *
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const typename Triangulation<dim, spacedim>::face_iterator &face,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);
  };



  /**
   * 这个类对子表面实现了 hp::FEFaceValues 对面的作用。
   * 进一步的文档请见那里。
   * @ingroup hp hpcollection
   *
   */
  template <int dim, int spacedim = dim>
  class FESubfaceValues
    : public hp::
        FEValuesBase<dim, dim - 1, dealii::FESubfaceValues<dim, spacedim>>
  {
  public:
    /**
     * 构造函数。用给定的参数初始化这个对象。
     *
     */
    FESubfaceValues(
      const hp::MappingCollection<dim, spacedim> &mapping_collection,
      const hp::FECollection<dim, spacedim> &     fe_collection,
      const hp::QCollection<dim - 1> &            q_collection,
      const UpdateFlags                           update_flags);


    /**
     * 构造函数。这个构造函数等同于另一个构造函数，只是它使对象隐含地使用
     * $Q_1$ 映射（即MappingQGeneric(1)类型的对象）。
     *
     */
    FESubfaceValues(const hp::FECollection<dim, spacedim> &fe_collection,
                    const hp::QCollection<dim - 1> &       q_collection,
                    const UpdateFlags                      update_flags);

    /**
     * 为给定的单元、面和子面重新初始化对象。
     * 调用后，你可以使用get_present_fe_values()函数得到一个与当前单元格对应的FESubfaceValues对象。
     * 对于这个FESubfaceValues对象，我们使用下面描述的附加参数来决定使用哪个有限元、映射和正交公式。它们的顺序是这样的：人们可能最想改变的参数排在前面。这些参数的规则如下。
     * 如果 @p q_index
     * 参数保持默认值，那么我们使用传递给这个类的构造函数的
     * hp::QCollection 中的正交公式，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。在这种情况下，
     * hp::FECollection.
     * 中的每个有限元都应该有一个相应的正交公式。
     * 作为一种特殊情况，如果正交集合只包含一个元素（如果想在hp-discretization中对所有有限元使用同一个正交对象，即使这可能不是最有效的），那么这个单一的正交将被使用，除非为这个参数指定一个不同的值。另一方面，如果给了这个参数一个值，它将覆盖
     * <code>cell-@>active_fe_index()</code>
     * 的选择或对单一正交的选择。        如果 @p mapping_index
     * 参数保持默认值，那么我们使用传递给该类构造函数的
     * hp::MappingCollection 中的映射对象，索引由
     * <code>cell-@>active_fe_index()</code>
     * 给出，即与有限元的索引相同。如上所述，如果映射集合只包含一个元素（如果想对hp-discretization中的所有有限元素使用
     * $Q_1$
     * 映射，这种情况很常见），那么这个单一的映射将被使用，除非为这个参数指定一个不同的值。
     *
     */
    template <bool lda>
    void
    reinit(const TriaIterator<DoFCellAccessor<dim, spacedim, lda>> &cell,
           const unsigned int                                       face_no,
           const unsigned int                                       subface_no,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);

    /**
     * 像前面的函数一样，但是对于非DoFHandler迭代器。这个函数存在的原因是，人们也可以将这个类用于Triangulation对象。
     * 由于 <code>cell-@>active_fe_index()</code>
     * 对Triangulation迭代器没有意义，这个函数从传递给这个对象的构造器的相关构造中选择第零个有限元、映射和正交对象。唯一的例外是，如果你为这最后三个参数中的任何一个指定了一个不同于默认值的值。
     *
     */
    void
    reinit(const typename Triangulation<dim, spacedim>::cell_iterator &cell,
           const unsigned int                                          face_no,
           const unsigned int subface_no,
           const unsigned int q_index       = numbers::invalid_unsigned_int,
           const unsigned int mapping_index = numbers::invalid_unsigned_int,
           const unsigned int fe_index      = numbers::invalid_unsigned_int);
  };

} // namespace hp


// -------------- inline and template functions --------------

namespace hp
{
  template <int dim, int q_dim, class FEValuesType>
  inline const FEValuesType &
  FEValuesBase<dim, q_dim, FEValuesType>::get_present_fe_values() const
  {
    return *fe_values_table(present_fe_values_index);
  }



  template <int dim, int q_dim, class FEValuesType>
  inline const FECollection<dim, FEValuesType::space_dimension> &
  FEValuesBase<dim, q_dim, FEValuesType>::get_fe_collection() const
  {
    return *fe_collection;
  }



  template <int dim, int q_dim, class FEValuesType>
  inline const MappingCollection<dim, FEValuesType::space_dimension> &
  FEValuesBase<dim, q_dim, FEValuesType>::get_mapping_collection() const
  {
    return *mapping_collection;
  }



  template <int dim, int q_dim, class FEValuesType>
  inline const QCollection<q_dim> &
  FEValuesBase<dim, q_dim, FEValuesType>::get_quadrature_collection() const
  {
    return q_collection;
  }



  template <int dim, int q_dim, class FEValuesType>
  inline UpdateFlags
  FEValuesBase<dim, q_dim, FEValuesType>::get_update_flags() const
  {
    return update_flags;
  }
} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif


