//include/deal.II-translator/hp/mapping_collection_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_mapping_collection_h
#define dealii_mapping_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/collection.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * 这个类实现了映射对象的集合，与 hp::FECollection
   * 实现有限元类集合的方式相同。
   * 它实现了doxygen文档中描述的 @ref hpcollection
   * 模块中的概念。
   * 尽管建议为hp-computation中使用的每个有限元种类提供一个适当的映射，但MappingCollection类实现了一个来自单一映射的转换构造器。
   * 因此，可以只为 hp::FEValues 类提供一个映射，而不是
   * hp::MappingCollection.
   * 类。这是为了方便用户，因为许多简单的几何形状不需要沿边界提供不同的映射来达到最佳收敛率。
   * 因此，提供一个单一的映射对象通常就足够了。参见
   * hp::FEValues 类中关于为给定单元选择映射的规则。
   * @ingroup hp hpcollection
   *
   */
  template <int dim, int spacedim = dim>
  class MappingCollection : public Collection<Mapping<dim, spacedim>>
  {
  public:
    /**
     * 默认构造函数。导致一个空的集合，以后可以用push_back()来填充。
     *
     */
    MappingCollection() = default;

    /**
     * 转换构造函数。这个构造函数从一个单一的映射中创建一个MappingCollection。如果需要的话，可以用push_back()添加更多的映射，尽管用同样的方式添加所有的映射可能会更清楚。
     *
     */
    explicit MappingCollection(const Mapping<dim, spacedim> &mapping);

    /**
     * 复制构造函数。
     *
     */
    MappingCollection(
      const MappingCollection<dim, spacedim> &mapping_collection);

    /**
     * 构造函数。这个构造函数从传递给构造函数的一个或多个映射对象创建一个MappingCollection。为了使这个调用有效，所有的参数都需要是从Mapping<dim,spacedim>类派生的类型。
     *
     */
    template <class... MappingTypes>
    explicit MappingCollection(const MappingTypes &... mappings);

    /**
     * 向MappingCollection添加一个新的映射。一般来说，你会希望对映射使用与你使用的
     * hp::FECollection 对象的元素相同的顺序。然而，与
     * hp::QCollection::push_back()
     * 函数讨论的相同的考虑因素也适用于当前的环境。
     * 这个类创建了一个给定映射对象的副本，也就是说，你可以做像<tt>push_back(MappingQ<dim>(3));</tt>的事情。这个内部拷贝后来在整个集合被销毁时被这个对象销毁。
     *
     */
    void
    push_back(const Mapping<dim, spacedim> &new_mapping);
  };


  /**
   * 库中的许多地方默认使用（bi-,tri-）线性映射，除非用户明确提供不同的映射来使用。在这些情况下，被调用的函数必须创建一个
   * $Q_1$
   * 映射对象，即一个MappingQGeneric(1)类型的对象。这是很昂贵的。在受影响的函数中创建这样的对象作为静态对象也是很昂贵的，因为静态对象在程序的整个生命周期中都不会被销毁，尽管它们只需要在代码第一次运行某个特定函数时创建一次。
   * 为了避免在整个库的这些上下文中创建（静态或动态）
   * $Q_1$ 映射对象，这个类定义了一个具有单一 $Q_1$
   * 映射对象的静态映射集合。然后这个集合可以在所有需要这种集合的地方使用。
   *
   */
  template <int dim, int spacedim = dim>
  struct StaticMappingQ1
  {
  public:
    /**
     * 公开可用的静态 $Q_1$  映射集合对象。
     *
     */
    static MappingCollection<dim, spacedim> mapping_collection;
  };


   /* --------------- inline functions ------------------- */ 

  template <int dim, int spacedim>
  template <class... MappingTypes>
  MappingCollection<dim, spacedim>::MappingCollection(
    const MappingTypes &... mappings)
  {
    static_assert(
      is_base_of_all<Mapping<dim, spacedim>, MappingTypes...>::value,
      "Not all of the input arguments of this function "
      "are derived from FiniteElement<dim,spacedim>!");

    // loop over all of the given arguments and add the mappings to
    // this collection. Inlining the definition of mapping_pointers causes
    // internal compiler errors on GCC 7.1.1 so we define it separately:
    const auto mapping_pointers = {
      (static_cast<const Mapping<dim, spacedim> *>(&mappings))...};
    for (const auto p : mapping_pointers)
      push_back(*p);
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif


