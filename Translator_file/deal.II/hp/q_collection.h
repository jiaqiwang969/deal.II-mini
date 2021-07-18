//include/deal.II-translator/hp/q_collection_0.txt
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

#ifndef dealii_q_collection_h
#define dealii_q_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/hp/collection.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * 这个类实现了一个正交对象的集合，与 hp::FECollection
   * 实现有限元类集合的方式相同。
   * 它实现了doxygen文档中描述的 @ref hpcollection
   * 模块中的概念。
   * @ingroup hp hpcollection
   *
   */
  template <int dim>
  class QCollection : public Collection<Quadrature<dim>>
  {
  public:
    /**
     * 默认构造函数。导致一个空的集合，以后可以用push_back()来填充。
     *
     */
    QCollection() = default;

    /**
     * 复制构造函数。
     *
     */
    template <int dim_in>
    QCollection(const QCollection<dim_in> &other);

    /**
     * 转换构造函数。这个构造函数从一个单一的正交法则创建一个QCollection。如果需要的话，可以用push_back()添加更多的正交公式，尽管用同样的方式添加所有的映射可能会更清楚。
     *
     */
    template <int dim_in>
    explicit QCollection(const Quadrature<dim_in> &quadrature);

    /**
     * 构造函数。这个构造函数从传递给构造函数的一个或多个正交对象创建一个QCollection。为了使这个调用有效，所有的参数需要是派生自类正交<dim>的类型。
     *
     */
    template <class... QTypes>
    explicit QCollection(const QTypes &... quadrature_objects);

    /**
     * 在QCollection中添加一个新的正交规则。在大多数情况下，你会希望按照元素被添加到
     * hp::FECollection
     * 中的相同顺序来添加正交规则，而这个正交规则集合就是为了这个正交规则。如果这样做，你将使用
     * hp::FECollection 和 hp::QCollection 对象的 hp::FEValues
     * 对象将自动选择相应的元素和正交公式。另一方面，在调用
     * hp::FEValues::reinit() 或 hp::FEFaceValues::reinit().
     * 时特别指定适当的索引，可以使用 hp::FECollection 和
     * hp::QCollection 对象中元素和正交公式的任意组合。
     * 在这些情况下， hp::FECollection 和 hp::QCollection
     * 对象的元素之间不需要有对应关系；在这种情况下它们甚至不需要有相同大小。
     * 顺便说一下，关于集合元素顺序的论点也可以针对
     * hp::MappingCollection 对象的元素提出。
     * 这个类创建了一个给定正交对象的副本，也就是说，你可以做类似<tt>push_back(QGauss<dim>(3));</tt>的事情。内部拷贝后来在整个集合被销毁时被这个对象销毁。
     *
     */
    template <int dim_in>
    void
    push_back(const Quadrature<dim_in> &new_quadrature);

    /**
     * 等价比较运算符。所有存储的正交对象都按顺序进行比较。
     *
     */
    bool
    operator==(const QCollection<dim> &q_collection) const;

    /**
     * 返回集合中所有元素的正交点的最大数量。这对初始化数组大多是有用的，以分配最大的内存量，当以后从这个集合中重新调整大小到一个特定的正交公式时，可能会用到这个内存量。
     *
     */
    unsigned int
    max_n_quadrature_points() const;

    /**
     * 例外
     *
     */
    DeclException0(ExcNoQuadrature);

  private:
    /**
     * 实数容器，用于存储指向不同正交对象的指针。
     *
     */
    std::vector<std::shared_ptr<const Quadrature<dim>>> quadratures;
  };



   /* --------------- inline functions ------------------- */ 

  template <int dim>
  template <int dim_in>
  QCollection<dim>::QCollection(const QCollection<dim_in> &other)
  {
    for (unsigned int i = 0; i < other.size(); ++i)
      push_back(other[i]);
  }



  template <int dim>
  template <class... QTypes>
  QCollection<dim>::QCollection(const QTypes &... quadrature_objects)
  {
    // loop over all of the given arguments and add the quadrature objects to
    // this collection. Inlining the definition of q_pointers causes internal
    // compiler errors on GCC 7.1.1 so we define it separately:
    if (is_base_of_all<Quadrature<dim>, QTypes...>::value)
      {
        const auto q_pointers = {
          (reinterpret_cast<const Quadrature<dim> *>(&quadrature_objects))...};
        for (const auto p : q_pointers)
          push_back(*p);
      }
    else if (is_base_of_all<Quadrature<1>, QTypes...>::value)
      {
        const auto q_pointers = {
          (reinterpret_cast<const Quadrature<1> *>(&quadrature_objects))...};
        for (const auto p : q_pointers)
          push_back(*p);
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }
  }



  template <int dim>
  inline unsigned int
  QCollection<dim>::max_n_quadrature_points() const
  {
    Assert(this->size() > 0,
           ExcMessage("You can't call this function for an empty collection"));

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).size());

    return max;
  }



  template <int dim>
  inline bool
  QCollection<dim>::operator==(const QCollection<dim> &q_collection) const
  {
    const unsigned int n_quadratures = this->size();
    if (n_quadratures != q_collection.size())
      return false;

    for (unsigned int i = 0; i < n_quadratures; ++i)
      if ((this->operator[](i) == q_collection[i]) == false)
        return false;

    return true;
  }



  template <int dim>
  template <int dim_in>
  inline QCollection<dim>::QCollection(const Quadrature<dim_in> &quadrature)
  {
    this->push_back(quadrature);
  }


  template <int dim>
  template <int dim_in>
  inline void
  QCollection<dim>::push_back(const Quadrature<dim_in> &new_quadrature)
  {
    Collection<Quadrature<dim>>::push_back(
      std::make_shared<const Quadrature<dim>>(new_quadrature));
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif


