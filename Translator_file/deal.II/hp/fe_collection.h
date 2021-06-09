//include/deal.II-translator/hp/fe_collection_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2021 by the deal.II authors
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

#ifndef dealii_fe_collection_h
#define dealii_fe_collection_h

#include <deal.II/base/config.h>

#include <deal.II/fe/component_mask.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values_extractors.h>

#include <deal.II/hp/collection.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * 该类作为DoFHandler中使用的有限元对象的集合。
   * 它实现了doxygen文档中描述的 @ref hpcollection
   * 模块中的概念。
   * 除了提供对集合元素的访问外，该类还提供对每个顶点、线等的最大自由度的访问，以便在使用与三角形单元相关的有限元时，在最坏的情况下分配尽可能多的内存。
   * 这个类还没有实现在一维情况下的使用（<tt>spacedim != dim
   * </tt>）。
   * @ingroup hp hpcollection
   *
   */
  template <int dim, int spacedim = dim>
  class FECollection : public Collection<FiniteElement<dim, spacedim>>
  {
  public:
    /**
     * 在hp-finite
     * element程序中考虑p-adaptivity时，需要建立一个有限元的层次结构，以确定细化的后续有限元和粗化的前期有限元。
     * 在这个结构中，我们提供了一个层次结构，默认情况下强加在所有FECollection对象上。
     *
     */
    struct DefaultHierarchy
    {
      /**
       * 返回 @p fe_collection. 中 @p fe_index 的后续索引 一旦到达
       * @p fe_collection
       * 的最后一个元素，在层次结构中就没有更高层次的元素，因此我们返回最后的值。
       *
       */
      static unsigned int
      next_index(const typename hp::FECollection<dim, spacedim> &fe_collection,
                 const unsigned int                              fe_index)
      {
        return ((fe_index + 1) < fe_collection.size()) ? fe_index + 1 :
                                                         fe_index;
      }

      /**
       * 返回 @p fe_collection. 中 @p fe_index 之前的索引 一旦到达
       * @p fe_collection
       * 的第一个元素，在层次结构中就没有较低层次的元素，因此我们返回第一个值。
       *
       */
      static unsigned int
      previous_index(
        const typename hp::FECollection<dim, spacedim> &fe_collection,
        const unsigned int                              fe_index)
      {
        (void)fe_collection;
        return (fe_index > 0) ? fe_index - 1 : fe_index;
      }
    };

    /**
     * 默认构造函数。导致一个空的集合，以后可以用push_back()来填充。建立一个有限元的层次结构，与它们在集合中的索引相对应。
     *
     */
    FECollection();

    /**
     * 转换构造函数。这个构造函数从一个单一的有限元创建一个FECollection。如果需要，可以用push_back()添加更多的有限元对象，尽管用同样的方式添加所有的映射可能会更清楚。
     *
     */
    explicit FECollection(const FiniteElement<dim, spacedim> &fe);

    /**
     * 构造函数。这个构造函数从传递给构造函数的一个或多个有限元对象创建一个FECollection。为了使这个调用有效，所有的参数都需要是派生自FiniteElement<dim,spacedim>类的类型。
     *
     */
    template <class... FETypes>
    explicit FECollection(const FETypes &... fes);

    /**
     * 构造函数。和上面一样，但是对于任何数量的元素。元素的指针以向量形式传递给这个构造函数。如上所述，参数所指向的有限元对象除了在内部创建副本外，实际上并不使用。因此，你可以在调用这个构造函数后立即再次删除这些指针。
     *
     */
    FECollection(const std::vector<const FiniteElement<dim, spacedim> *> &fes);

    /**
     * 拷贝构造函数。
     *
     */
    FECollection(const FECollection<dim, spacedim> &) = default;

    /**
     * 移动构造函数。
     * @note
     * 标准数据类型的实现可能会随着不同的库而改变，所以它们的移动成员可能会或不会被标记为非抛出。
     * 我们需要根据其成员变量显式地设置noexcept指定器，以便仍然获得性能优势（并满足clang-tidy）。
     *
     */
    FECollection(FECollection<dim, spacedim> &&) noexcept(
      std::is_nothrow_move_constructible<
        std::vector<std::shared_ptr<const FiniteElement<dim, spacedim>>>>::value
        &&std::is_nothrow_move_constructible<std::function<
          unsigned int(const typename hp::FECollection<dim, spacedim> &,
                       const unsigned int)>>::value) = default;

    /**
     * 移动赋值运算符。
     *
     */
    FECollection<dim, spacedim> &
    operator=(FECollection<dim, spacedim> &&) = default; // NOLINT

    /**
     * 等价比较运算符。所有存储的FiniteElement对象都按顺序进行比较。
     *
     */
    bool
    operator==(const FECollection<dim, spacedim> &fe_collection) const;

    /**
     * 非等价比较运算符。所有存储的FiniteElement对象按顺序进行比较。
     *
     */
    bool
    operator!=(const FECollection<dim, spacedim> &fe_collection) const;

    /**
     * 添加一个有限元素。这个函数生成一个给定元素的副本，即你可以做类似<tt>push_back(FE_Q<dim>(1));</tt>的事情。
     * 这个内部拷贝后来在整个集合被销毁时被这个对象销毁。
     * 当一个新的元素被添加时，它需要拥有与已经在集合中的所有其他元素相同数量的向量成分。
     *
     */
    void
    push_back(const FiniteElement<dim, spacedim> &new_fe);

    /**
     * 返回这个集合中的有限元素的向量分量的数量。 这个数字对于集合中的所有元素必须是相同的。        这个函数调用 FiniteElement::n_components. 更多信息见 @ref GlossComponent "术语表"
     * 。
     *
     */
    unsigned int
    n_components() const;

    /**
     * 返回这个集合中的有限元素的向量块的数量。虽然这个类保证了存储在其中的所有元素都有相同数量的向量分量，但是对于每个元素所组成的块数却没有这样的保证（一个元素的块数可能少于向量分量，更多信息见 @ref GlossBlock  "术语表"
     * ）。例如，你可能有一个FECollection对象，它存储了一份带有
     * <code>dim</code>
     * FE_Q对象的FESystem和一份FE_RaviartThomas元素的副本。两者都有
     * <code>dim</code> 向量成分，但前者有 <code>dim</code>
     * 块，后者只有一个。因此，如果所有元素的块数不一样，这个函数将抛出一个断言。如果它们相同，该函数返回的结果是
     * FiniteElement::n_blocks().  。
     *
     */
    unsigned int
    n_blocks() const;

    /**
     * 返回 FiniteElement::get_degree()
     * 在这个集合的所有元素上返回的值的最大值。
     *
     */
    unsigned int
    max_degree() const;

    /**
     * 返回此集合所有元素中每个顶点的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_vertex() const;

    /**
     * 返回这个集合的所有元素中每条线的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_line() const;

    /**
     * 返回这个集合中所有元素的每个四边形的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_quad() const;

    /**
     * 返回这个集合中所有元素中每个六度的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_hex() const;

    /**
     * 返回这个集合的所有元素中每个面的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_face() const;

    /**
     * 返回这个集合的所有元素中每个单元的最大自由度数。
     *
     */
    unsigned int
    max_dofs_per_cell() const;


    /**
     * 返回这个集合中的所有元素是否以新的方式实现了悬挂节点约束，这必须用于使元素
     * "hp-compatible"。如果不是这样，该函数返回false，这意味着FECollection中至少有一个元素不支持新的面孔接口约束。另一方面，如果这个方法确实返回true，这并不意味着hp-方法会起作用。
     * 这种行为与以下事实有关，即提供新式悬挂节点约束的FiniteElement类可能仍然没有为所有可能的情况提供这些约束。如果FE_Q和FE_RaviartThomas元素包含在FECollection中，并且都正确实现了get_face_interpolation_matrix方法，这个方法将返回true。
     * 但是get_face_interpolation_matrix可能仍然无法找到这两个元素之间的插值矩阵。
     *
     */
    bool
    hp_constraints_are_implemented() const;

    /**
     * 返回这个FECollection中支配所有与所提供的索引集相关的元素的有限元素的索引
     * @p fes.
     * 你可以在其各自的类文档或其继承的成员函数的实现中找到关于有限元素的支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元素可能支配也可能不支配它自己（例如FE_Nothing元素）。
     * 例如，如果一个FEC集合由`{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}`元素组成，我们要寻找支配这个集合中间元素的有限元素（即
     * @p fes
     * 为`{1,2}`），那么答案是`{FE_Q(1),FE_Q(2)`，因此这个函数将返回它们在FEC集合中的索引，即`{0,1}`。
     * @p codim
     * 参数描述了被调查的子空间的码率，并指定它受此比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    std::set<unsigned int>
    find_common_fes(const std::set<unsigned int> &fes,
                    const unsigned int            codim = 0) const;

    /**
     * 返回此FECollection中被与所提供的指数集相关的所有元素支配的有限元的指数
     * @p fes.
     * 你可以在其各自的类文件或其继承的成员函数的实现中找到关于有限元的支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元可能支配也可能不支配自己（例如，FE_Nothing元素）。
     * 例如，如果一个FEC集合由`{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}`元素组成，我们要寻找被这个集合的中间元素支配的有限元素（即。
     * @p fes
     * 是`{1,2}`，那么答案是`{FE_Q(3),FE_Q(4)`，因此这个函数将返回它们在FEC集合中的索引，即`{2,3}`。
     * @p codim
     * 参数描述了被调查的子空间的码率，并指定它受此比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    std::set<unsigned int>
    find_enclosing_fes(const std::set<unsigned int> &fes,
                       const unsigned int            codim = 0) const;

    /**
     * 返回所提供的指数集 @p fes
     * 中的一个有限元素的指数，该元素支配着这个非常集合中的所有其他元素。
     * 你可以在其各自的类文件中或在其继承的成员函数的实现中找到关于有限元的支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元可能支配也可能不支配自己（例如FE_Nothing元素）。
     * 如果这个集合正好由一个元素组成，我们认为它是支配性的，并返回其相应的索引。此外，如果函数根本无法找到一个有限元素，则返回
     * numbers::invalid_unsigned_int.
     * 例如，如果一个FEC集合由`{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}`元素组成，我们要在这个集合的中间元素中寻找主导的有限元素（即。
     * @p fes
     * 是`{1,2}`），那么答案是FE_Q(2)，因此这个函数将返回它在FEC集合中的索引，即`1`。
     * 当然，有可能存在不止一个元素支配着所有被选中的元素。例如，如果集合由`{FE_Q(1),FE_Q(1),FE_Q(2),FE_Q(2)}`组成，并且`fes`涵盖所有的索引，那么可以返回0或1。
     * 在这种情况下，该函数要么返回 "0"，要么返回
     * "1"，因为两者之间不存在平局。         @p codim
     * 参数描述了被调查的子空间的二维度，并指定它要接受这种比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    unsigned int
    find_dominating_fe(const std::set<unsigned int> &fes,
                       const unsigned int            codim = 0) const;

    /**
     * 返回所提供的指数集 @p fes
     * 中的一个有限元素的指数，该元素被这个非常集合的所有其他元素所支配。
     * 你可以在其各自的类文件中或其继承成员函数的实现中找到关于有限元的支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元可能支配自己，也可能不支配自己（例如FE_Nothing元素）。
     * 如果这个集合正好由一个元素组成，我们认为它是被支配的，并返回其相应的索引。此外，如果函数根本无法找到有限元，则返回
     * numbers::invalid_unsigned_int.
     * 例如，如果一个FEC集合由`{FE_Q(1),FE_Q(2),FE_Q(3),FE_Q(4)}`元素组成，我们在这个集合的中间元素中寻找支配的有限元（即。
     * @p fes
     * 是`{1,2}`），那么答案是FE_Q(3)，因此这个函数将返回它在FEC集合中的索引，即`2`。
     * 当然，有可能存在不止一个元素被所有选定的元素所支配。例如，如果集合由`{FE_Q(1),FE_Q(1),FE_Q(2),FE_Q(2)}`组成，并且`fes`覆盖了所有的索引，那么可以返回2或3。
     * 在这种情况下，该函数要么返回`2`，要么返回`3`，因为两者之间没有平局。
     * @p codim
     * 参数描述了被调查的子空间的维度，并指定它要接受这种比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    unsigned int
    find_dominated_fe(const std::set<unsigned int> &fes,
                      const unsigned int            codim = 0) const;

    /**
     * 返回所提供的指数集 @p fes
     * 中的一个有限元素的指数，该元素支配着这个非常集合中的所有其他元素。如果我们没有成功，我们就在整个集合上扩大搜索范围，挑选最不占优势的，也就是描述最大的有限元空间的元素，所提供的集合
     * @p fes 的所有有限元都是其中的一部分。
     * 你可以在其各自的类文件中或其继承的成员函数的实现中找到关于有限元的支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元可以支配也可以不支配自己（例如FE_Nothing元素）。
     * 如果这个集合正好由一个元素组成，我们认为它是被支配的，并返回其相应的索引。此外，如果函数根本无法找到一个有限元，则返回
     * numbers::invalid_unsigned_int.   @p codim
     * 参数描述了所研究的子空间的二维度，并指定它受此比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    unsigned int
    find_dominating_fe_extended(const std::set<unsigned int> &fes,
                                const unsigned int            codim = 0) const;

    /**
     * 返回所提供的指数集 @p fes
     * 中的一个有限元素的指数，该元素被这个非常集的所有其他元素所支配。如果我们没有成功，我们就在整个集合上扩展搜索，挑选出最被支配的，也就是描述最小的有限元空间的元素，该空间包括所提供集合的所有有限元
     * @p fes.
     * 你可以在其各自的类文档中或其继承成员函数的实现中找到关于有限元支配行为的信息
     * FiniteElement::compare_for_domination().
     * 考虑到一个有限元可能支配也可能不支配自己（例如，FE_Nothing元素）。
     * 如果这个集合正好由一个元素组成，我们就认为它是支配性的，并返回其相应的索引。此外，如果该函数根本无法找到一个有限元素，则返回
     * numbers::invalid_unsigned_int.   @p codim
     * 参数描述了被调查的子空间的二维度，并指定其受此比较。更多信息见
     * FiniteElement::compare_for_domination() 。
     *
     */
    unsigned int
    find_dominated_fe_extended(const std::set<unsigned int> &fes,
                               const unsigned int            codim = 0) const;

    /**
     * 确定有限元层次的集合函数，即一个函数 @p next
     * 返回给定的有限元后的索引，一个函数 @p prev
     * 返回前一个。        这两个函数都需要一个
     * hp::FECollection
     * 来传递有限元索引，在其基础上找到并返回新的索引。
     * @note
     * 传递和返回的索引都必须在这个集合的索引范围内有效，即在[0,
     * size()]内。
     *
     */
    void
    set_hierarchy(const std::function<unsigned int(
                    const typename hp::FECollection<dim, spacedim> &,
                    const unsigned int)> &next,
                  const std::function<unsigned int(
                    const typename hp::FECollection<dim, spacedim> &,
                    const unsigned int)> &prev);

    /**
     * 设置与集合中每个有限元的索引相对应的默认层次结构。
     * 这个默认层次是通过函数 DefaultHierarchy::next_index() 和
     * DefaultHierarchy::previous_index(). 建立的。
     *
     */
    void
    set_default_hierarchy();

    /**
     * 返回一个对应于注册层次结构的FE指数序列，以升序排列，即FE指数从低到高排序。
     * 用set_hierarchy()注册的一个自定义层次结构可以有多个FE指数序列。该函数将返回包含用户提供的索引
     * @p fe_index
     * 的序列，该索引可以位于序列内的任何位置。通过set_default_hierarchy()设置的默认层次结构，对应于升序的FE指数，只包含一个序列。
     * 例如，这个函数可以用来验证你所提供的层次结构是否涵盖了所有的元素，并符合所需的顺序。
     * 如果返回的容器的大小等于这个对象的元素数，则只存在一个FE指数序列，即
     * FECollection::size().  。
     *
     */
    std::vector<unsigned int>
    get_hierarchy_sequence(const unsigned int fe_index = 0) const;

    /**
     * %函数返回层次结构中给定的 @p fe_index
     * 之后的有限元的索引。        默认情况下，将返回 @p
     * fe_index 之后的索引。如果 @p fe_index
     * 已经对应于最后一个索引，将返回最后一个索引。
     * 可以通过成员函数set_hierachy()提供一个自定义的层次结构。
     *
     */
    unsigned int
    next_in_hierarchy(const unsigned int fe_index) const;

    /**
     * %函数返回层次结构中给定 @p fe_index
     * 之前的有限元的索引。        默认情况下，将返回 @p
     * fe_index 之前的索引。如果 @p fe_index
     * 已经对应于第一个索引，第一个索引将被返回。
     * 可以通过成员函数set_hierachy()提供一个自定义的层次结构。
     *
     */
    unsigned int
    previous_in_hierarchy(const unsigned int fe_index) const;

    /**
     * 返回一个分量掩码，其元素数量与此对象的向量分量相同，并且其中正好有一个分量是真实的，与给定的参数相对应。
     * @note  这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证了它是否从存储在这个FECollection中的每一个元素中得到了相同的结果。如果不是这样的话，它会抛出一个异常。
     * @param  标量
     * 一个代表该有限元的单一标量矢量分量的对象。
     * @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那一个。
     *
     */
    ComponentMask
    component_mask(const FEValuesExtractors::Scalar &scalar) const;

    /**
     * 返回一个分量掩码，其元素数与此对象的向量分量相同，其中与给定参数对应的
     * <code>dim</code> 分量为真。
     * @note 这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样，它会抛出一个异常。
     * @param  矢量
     * 一个表示该有限元的微弱矢量成分的对象。      @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的分量。
     *
     */
    ComponentMask
    component_mask(const FEValuesExtractors::Vector &vector) const;

    /**
     * 返回一个分量掩码，其元素数与此对象的向量分量相同，其中与给定参数对应的
     * <code>dim*(dim+1)/2</code> 分量为真。
     * @note 这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样，它会抛出一个异常。
     * @param  sym_tensor
     * 一个表示该有限元的dim*(dim+1)/2组件的对象，这些组件共同被解释为形成一个对称张量。
     * @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的分量。
     *
     */
    ComponentMask
    component_mask(
      const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
     * @note  这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样，它会抛出一个异常。
     * @param  block_mask 选择有限元单个块的掩码  @return
     * 选择那些与输入参数的选定块对应的组件的掩码。
     *
     */
    ComponentMask
    component_mask(const BlockMask &block_mask) const;

    /**
     * 返回一个块掩码，其元素数与此对象的块数相同，并且其中正好有一个与给定参数相对应的成分是真的。更多信息请参见  @ref GlossBlockMask  "术语表"
     * 。
     * @note
     * 这个函数只有在参数所引用的标量包含一个完整的块时才会成功。换句话说，例如，如果你传递了一个单一
     * $x$
     * 速度的提取器，并且这个对象代表一个FE_RaviartThomas对象，那么你选择的单一标量对象是一个更大的块的一部分，因此没有代表它的块屏蔽。然后该函数将产生一个异常。
     * @note  这个函数相当于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样的话，它会抛出一个异常。
     * @param  标量
     * 一个代表该有限元的单个标量向量分量的对象。
     * @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那一个。
     *
     */
    BlockMask
    block_mask(const FEValuesExtractors::Scalar &scalar) const;

    /**
     * 返回一个分量掩码，其元素数与此对象的向量分量相同，并且其中对应于给定参数的 <code>dim</code> 分量为真。更多信息见 @ref GlossBlockMask  "术语表"
     * 。
     * @note  这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证是否从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样的话，它会抛出一个异常。
     * @note  同样的注意事项适用于上述函数的版本。
     * 作为参数传递的提取器对象必须使其对应于完整的块，而不是分割这个元素的块。
     * @param  vector 一个表示该有限元的dim向量分量的对象。
     * @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的分量。
     *
     */
    BlockMask
    block_mask(const FEValuesExtractors::Vector &vector) const;

    /**
     * 返回一个分量掩码，其元素数与此对象的向量分量相同，其中与给定参数对应的 <code>dim*(dim+1)/2</code> 分量为真。更多信息见 @ref GlossBlockMask  "术语表"
     * 。
     * @note  同样的注意事项适用于上述函数的版本。
     * 作为参数传递的提取器对象必须使其对应于完整的块，而不是分割此元素的块。
     * @note  这个函数等同于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样，它会抛出一个异常。
     * @param  sym_tensor
     * 一个代表该有限元的dim*(dim+1)/2组件的对象，这些组件共同被解释为形成一个对称张量。
     * @return
     * 一个分量掩码，在所有分量中都是假的，除了与参数相对应的那些。
     *
     */
    BlockMask
    block_mask(const FEValuesExtractors::SymmetricTensor<2> &sym_tensor) const;

    /**
     * @note
     * 这个函数只有在参数所引用的组件包含完整的块时才会成功。换句话说，例如，如果你传递了一个单一
     * $x$
     * 速度的组件掩码，而这个对象代表一个FE_RaviartThomas对象，那么你选择的单一组件是一个更大的块的一部分，因此，没有代表它的块掩码。然后该函数将产生一个异常。
     * @note  这个函数相当于 FiniteElement::component_mask()
     * ，参数相同。它验证了它从存储在这个FECollection中的每一个元素中得到相同的结果。如果不是这样，它会抛出一个异常。
     * @param component_mask 选择有限元个别组件的掩码  @return
     * 选择那些与输入参数的选定块对应的掩码。
     *
     */
    BlockMask
    block_mask(const ComponentMask &component_mask) const;

    /**
     * @name  异常情况  @{
     *
     */

    /**
     * 异常情况
     * @ingroup Exceptions
     *
     */
    DeclException0(ExcNoFiniteElements);

    /**
     * @}
     *
     */

  private:
    /**
     * %函数返回层次结构中给定元素之后的有限元素的索引。
     *
     */
    std::function<unsigned int(const typename hp::FECollection<dim, spacedim> &,
                               const unsigned int)>
      hierarchy_next;

    /**
     * 返回在层次结构中给定的有限元之前的有限元的索引的函数。
     *
     */
    std::function<unsigned int(const typename hp::FECollection<dim, spacedim> &,
                               const unsigned int)>
      hierarchy_prev;
  };



   /* --------------- inline functions ------------------- */ 

  template <int dim, int spacedim>
  template <class... FETypes>
  FECollection<dim, spacedim>::FECollection(const FETypes &... fes)
  {
    static_assert(
      is_base_of_all<FiniteElement<dim, spacedim>, FETypes...>::value,
      "Not all of the input arguments of this function "
      "are derived from FiniteElement<dim,spacedim>!");

    // loop over all of the given arguments and add the finite elements to
    // this collection. Inlining the definition of fe_pointers causes internal
    // compiler errors on GCC 7.1.1 so we define it separately:
    const auto fe_pointers = {
      (static_cast<const FiniteElement<dim, spacedim> *>(&fes))...};
    for (const auto p : fe_pointers)
      push_back(*p);
  }



  template <int dim, int spacedim>
  inline unsigned int
  FECollection<dim, spacedim>::n_components() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    // note that there is no need
    // here to enforce that indeed
    // all elements have the same
    // number of components since we
    // have already done this when
    // adding a new element to the
    // collection.

    return this->operator[](0).n_components();
  }



  template <int dim, int spacedim>
  inline bool
  FECollection<dim, spacedim>::
  operator==(const FECollection<dim, spacedim> &fe_collection) const
  {
    const unsigned int n_elements = this->size();
    if (n_elements != fe_collection.size())
      return false;

    for (unsigned int i = 0; i < n_elements; ++i)
      if (!(this->operator[](i) == fe_collection[i]))
        return false;

    return true;
  }



  template <int dim, int spacedim>
  inline bool
  FECollection<dim, spacedim>::
  operator!=(const FECollection<dim, spacedim> &fe_collection) const
  {
    return !(*this == fe_collection);
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_degree() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).degree);

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_vertex() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_vertex());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_line() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_line());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_quad() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).max_dofs_per_quad());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_hex() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_hex());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_face() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).max_dofs_per_face());

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim, spacedim>::max_dofs_per_cell() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).n_dofs_per_cell());

    return max;
  }


  template <int dim, int spacedim>
  bool
  FECollection<dim, spacedim>::hp_constraints_are_implemented() const
  {
    Assert(this->size() > 0, ExcNoFiniteElements());

    for (unsigned int i = 0; i < this->size(); ++i)
      if (this->operator[](i).hp_constraints_are_implemented() == false)
        return false;

    return true;
  }


} // namespace hp

DEAL_II_NAMESPACE_CLOSE

#endif


