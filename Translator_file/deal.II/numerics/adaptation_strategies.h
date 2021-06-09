//include/deal.II-translator/numerics/adaptation_strategies_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_adaptation_strategies_h
#define dealii_adaptation_strategies_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/grid/tria_accessor.h>

#include <algorithm>
#include <numeric>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * 当数据在适应过程中被转移时，决定如何处理从旧网格上的前单元被改变为新网格上的当前单元的数据并不简单。或者换句话说，数据应该如何存储在适应的网格上的单元中。
 * 在这个命名空间中，我们提供了一些应对这个问题的策略。这些策略可以传递给CellDataTransfer和
 * parallel::distributed::CellDataTransfer 构造函数。
 *
 *
 */
namespace AdaptationStrategies
{
  /**
   * 对于细化，所有的策略都接受父细胞及其相关数据。它们返回一个向量，包含父单元格将被细化到的每个单独的子单元格的数据。
   * 子女数据向量中的数值排序与调用 TriaAccessor::child_index.
   * 时的索引相对应。
   *
   */
  namespace Refinement
  {
    /**
     * 返回一个包含每个子单元的父单元数据副本的向量。        @f[
     * d_{K_c} = d_{K_p}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f]
     *
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    preserve(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
               &              parent,
             const value_type parent_value);

    /**
     * 返回一个包含父单元数据的向量，该向量被平均分配给所有子单元。        @f[
     * d_{K_c} = d_{K_p} / n_\text{children}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f] 这个策略保留了适应前后相应的全局数据向量的 $l_1$  -norm。
     *
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    split(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
            &              parent,
          const value_type parent_value);

    /**
     * 返回一个向量，其中包含父单元的平方数据被平均分配给所有子单元的平方。        @f[
     * d_{K_c}^2 = d_{K_p}^2 / n_\text{children}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f]该策略保留了适应前后相应全局数据向量的 $l_2$ -norm。
     *
     */
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
              &              parent,
            const value_type parent_value);
  } // namespace Refinement

  /**
   * 对于粗化，所有的策略都接受父单元和属于其前子单元的数据向量。它们返回将被分配给父单元的值。
   * 子代数据向量中的值的排序与调用 TriaAccessor::child_index.
   * 时的索引相对应。
   *
   */
  namespace Coarsening
  {
    /**
     * 检查所有子单元的数据是否匹配，并返回第一个子单元的值。        @f[
     * d_{K_p} = d_{K_c}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f]
     *
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    check_equality(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator
        &                            parent,
      const std::vector<value_type> &children_values);

    /**
     * 返回所有孩子的数据之和。        @f[
     * d_{K_p} = \sum d_{K_c}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f] 这个策略在适应前后都保留了相应的全局数据向量的 $l_1$  -norm。
     *
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    sum(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &                            parent,
        const std::vector<value_type> &children_values);

    /**
     * 返回所有子代数据的 $ l_2 $  -norm。        @f[
     * d_{K_p}^2 = \sum d_{K_c}^2
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f] 这个策略保留了适应前后相应全局数据向量的 $l_2$  -norm。
     *
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
              &                            parent,
            const std::vector<value_type> &children_values);

    /**
     * 返回所有孩子上的数据的平均值。        @f[
     * d_{K_p} = \sum d_{K_c} / n_\text{children}
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f]
     *
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    mean(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
           &                            parent,
         const std::vector<value_type> &children_values);

    /**
     * 返回所有孩子的数据的最大值。        @f[
     * d_{K_p} = \max \left( d_{K_c} \right)
     * \qquad
     * \forall K_c \text{ children of } K_p
     * @f]
     *
     */
    template <int dim, int spacedim, typename value_type>
    value_type
    max(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
          &                            parent,
        const std::vector<value_type> &children_values);
  } // namespace Coarsening
} // namespace AdaptationStrategies



 /* ---------------- template functions ---------------- */ 

#ifndef DOXYGEN

namespace AdaptationStrategies
{
  namespace Refinement
  {
    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    preserve(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
               &              parent,
             const value_type parent_value)
    {
      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(), parent_value);
    }



    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    split(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
            &              parent,
          const value_type parent_value)
    {
      static_assert(std::is_arithmetic<value_type>::value &&
                      !std::is_same<value_type, bool>::value,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(),
                                     parent_value / parent->n_children());
    }



    template <int dim, int spacedim, typename value_type>
    std::vector<value_type>
    l2_norm(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
              &              parent,
            const value_type parent_value)
    {
      static_assert(std::is_arithmetic<value_type>::value &&
                      !std::is_same<value_type, bool>::value,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(parent->n_children() > 0, ExcInternalError());
      return std::vector<value_type>(parent->n_children(),
                                     parent_value /
                                       std::sqrt(parent->n_children()));
    }
  } // namespace Refinement



  namespace Coarsening
  {
    template <int dim, int spacedim, typename value_type>
    value_type
    check_equality(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
      const std::vector<value_type> &children_values)
    {
      Assert(!children_values.empty(), ExcInternalError());

      const auto first_child = children_values.cbegin();
      for (auto other_child = first_child + 1;
           other_child != children_values.cend();
           ++other_child)
        Assert(*first_child == *other_child,
               ExcMessage(
                 "Values on cells that will be coarsened are not equal!"));

      return *first_child;
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    sum(const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
        const std::vector<value_type> &children_values)
    {
      static_assert(std::is_arithmetic<value_type>::value &&
                      !std::is_same<value_type, bool>::value,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(!children_values.empty(), ExcInternalError());
      return std::accumulate(children_values.cbegin(),
                             children_values.cend(),
                             static_cast<value_type>(0));
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    l2_norm(
      const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
      const std::vector<value_type> &children_values)
    {
      static_assert(std::is_arithmetic<value_type>::value &&
                      !std::is_same<value_type, bool>::value,
                    "The provided value_type may not meet the requirements "
                    "of this function.");

      Assert(!children_values.empty(), ExcInternalError());
      return std::sqrt(std::inner_product(children_values.cbegin(),
                                          children_values.cend(),
                                          children_values.cbegin(),
                                          static_cast<value_type>(0)));
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    mean(const typename dealii::Triangulation<dim, spacedim>::cell_iterator
           &                            parent,
         const std::vector<value_type> &children_values)
    {
      return sum<dim, spacedim, value_type>(parent, children_values) /
             children_values.size();
    }



    template <int dim, int spacedim, typename value_type>
    value_type
    max(const typename dealii::Triangulation<dim, spacedim>::cell_iterator &,
        const std::vector<value_type> &children_values)
    {
      Assert(!children_values.empty(), ExcInternalError());
      return *std::max_element(children_values.cbegin(),
                               children_values.cend());
    }
  } // namespace Coarsening
} // namespace AdaptationStrategies

#endif

DEAL_II_NAMESPACE_CLOSE

#endif  /* dealii_adaptation_strategies_h */ 


