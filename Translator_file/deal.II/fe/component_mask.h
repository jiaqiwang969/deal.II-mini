//include/deal.II-translator/fe/component_mask_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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

#ifndef dealii_fe_component_mask_h
#define dealii_fe_component_mask_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <algorithm>
#include <iosfwd>
#include <vector>

DEAL_II_NAMESPACE_OPEN



/**
 * 该类代表一个掩码，可用于选择有限元的单个矢量分量（另见 @ref GlossComponentMask "该术语条目"
 * ）。它通常有与有限元的矢量分量一样多的元素，人们可以使用
 * <code>operator[]</code> 来查询某个特定分量是否被选中。
 *
 *
 * @note  "掩码 "代表一个具有 @p true 和 @p false
 * 项的数据结构，通常用于启用或禁用某个特定向量分量的操作。根据这个定义，禁用的向量组件仍然存在
 *
 * - 它们只是没有被触及。因此，当你应用分量掩码来插值一个带有 $C$ 矢量分量的问题时（仅选择一个例子），描述边界值的输入参数仍然必须提供 $C$ 分量，即使掩码说我们只想把这些分量的一个子集插值到有限元空间。换句话说，分量掩码不代表<i>reduction</i>操作，它代表<i>selection</i>。
 * 这类对象在很多地方被使用，人们希望将操作限制在某个分量子集上，例如在 DoFTools::make_zero_boundary_values()
 * 或 VectorTools::interpolate_boundary_values(). 中。
 * 这些对象可以手工创建，或者更简单，要求有限元使用代码从某些选定的分量中生成一个分量掩码，例如我们创建一个只表示斯托克斯元速度分量的掩码（见
 * @ref vector_valued  ）。
 *
 * @code
 * // Q2 element for the velocities, Q1 element for the pressure
 * FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,
 *                          FE_Q<dim>(1), 1);
 * FEValuesExtractors::Scalar pressure(dim);
 * ComponentMask pressure_mask = stokes_fe.component_mask (pressure);
 * @endcode
 * 结果是一个分量掩码，在2d中，它的值是<code>[false, false,
 * true]</code>。同样地，使用
 *
 * @code
 * FEValuesExtractors::Vector velocities(0);
 * ComponentMask velocity_mask = stokes_fe.component_mask (velocities);
 * @endcode
 * 在2d中会产生一个 <code>[true, true, false]</code>
 * 的掩码。当然，在3D中，其结果将是 <code>[true, true, true,
 * false]</code>  。
 *
 *
 * @ingroup fe
 *
 * @ingroup vector_valued
 *
 */
class ComponentMask
{
public:
  /**
   * 初始化一个组件掩码。默认情况下，一个组件掩码代表一组<i>all</i>被选中的组件，也就是说，调用这个构造函数的结果是一个组件掩码，每当问到一个组件是否被选中时，总是返回
   * <code>true</code> 。
   *
   */
  ComponentMask() = default;

  /**
   * 用参数指定的一组选定的组件初始化这个类型的对象。
   * @param  component_mask 一个 <code>true/false</code>
   * 项的向量，决定有限元的哪些分量被选择。如果给定矢量的长度为零，则解释为<i>every</i>分量被选中的情况。
   *
   */
  ComponentMask(const std::vector<bool> &component_mask);

  /**
   * 用一定数量的元素初始化分量掩码，这些元素要么是真，要么是假。
   * @param  n_components 这个掩码的元素数量  @param  initializer
   * 这些元素中的每一个应该具有的值：要么是真，要么是假。
   *
   */
  ComponentMask(const unsigned int n_components, const bool initializer);

  /**
   * 将掩码中的某个条目设置为一个值。
   *
   */
  void
  set(const unsigned int index, const bool value);

  /**
   * 如果这个组件掩码已经被初始化为一个大小大于0的掩码，那么返回这个对象所代表的掩码的大小。
   * 另一方面，如果这个掩码已被初始化为一个空对象，代表一个对每个元素都是真的掩码（即，如果这个对象在调用
   * represents_the_all_selected_mask()时将返回true），那么返回0，因为没有明确的大小。
   *
   */
  unsigned int
  size() const;

  /**
   * 返回一个特定的组件是否被这个掩码所选择。如果这个掩码代表了选择<i>all
   * components</i>的对象的情况（例如，如果它是用默认的构造函数创建的，或者是从bool类型的空向量转换而来的），那么无论给定的参数是什么，这个函数都返回true。
   * @param  component_index
   * 该函数应返回该组件是否被选中的索引。如果这个对象代表一个掩码，其中所有组件总是被选中，那么这里允许任何索引。
   * 否则，给定的索引需要在零和该掩码所代表的组件数量之间。
   *
   */
  bool operator[](const unsigned int component_index) const;

  /**
   * 返回这个分量掩码是否正好代表 <code>n</code>
   * 分量的掩码。如果它被初始化为一个正好有 <code>n</code>
   * entries of type <code>bool</code> 的向量（在这种情况下， @p n
   * 必须等于size()的结果），或者它被初始化为一个空向量（或者使用默认构造函数），在这种情况下，它可以代表一个有任意数量组件的掩码，并且将总是说一个组件被选中，这就是真的。
   *
   */
  bool
  represents_n_components(const unsigned int n) const;

  /**
   * 返回被这个掩码选中的组件的数量。
   * 由于空的组件掩码代表每个组件都会返回 <code>true</code>
   * ，这个函数可能不知道组件掩码的真实大小，因此它需要一个参数来表示组件的总数量。
   * 如果该对象已经被初始化为一个非空的掩码（即，如果size()函数返回大于0的东西，或者等价地，如果respresent_the_all_selected_mask()返回false），那么该参数可以被省略，而size()的结果被取走。
   *
   */
  unsigned int
  n_selected_components(const unsigned int overall_number_of_components =
                          numbers::invalid_unsigned_int) const;

  /**
   * 返回第一个被选中的组件的索引。该参数存在的原因与n_selected_components()函数存在的原因相同。
   * 如果没有任何组件被选中，该函数会抛出一个异常。
   *
   */
  unsigned int
  first_selected_component(const unsigned int overall_number_of_components =
                             numbers::invalid_unsigned_int) const;

  /**
   * 如果这个掩码代表一个默认构建的掩码，对应于所有组件被选中的掩码，则返回true。如果为真，那么size()函数将返回0。
   *
   */
  bool
  represents_the_all_selected_mask() const;

  /**
   * 返回一个包含由当前对象选择的组件和作为参数传递的组件的联合体的组件掩码。
   *
   */
  ComponentMask
  operator|(const ComponentMask &mask) const;

  /**
   * 返回一个组件掩码，该掩码只包含那些在当前对象以及作为参数传递的对象中都被设置的元素。
   *
   */
  ComponentMask operator&(const ComponentMask &mask) const;

  /**
   * 返回这个对象和参数是否相同。
   *
   */
  bool
  operator==(const ComponentMask &mask) const;

  /**
   * 返回这个对象和参数是否不相同。
   *
   */
  bool
  operator!=(const ComponentMask &mask) const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcNoComponentSelected,
                   "The number of selected components in a mask "
                   "must be greater than zero.");

private:
  /**
   * 实际的组件掩码。
   *
   */
  std::vector<bool> component_mask;

  // make the output operator a friend so it can access
  // the component_mask array
  friend std::ostream &
  operator<<(std::ostream &out, const ComponentMask &mask);
};


/**
 * 将一个组件掩码写到输出流中。如果组件掩码代表所有的组件都被选中，而没有指定掩码的特定大小，那么它就会将字符串<code>[所有组件被选中]</code>写到流中。否则，它将以
 * <code>[true,true,true,false]</code> 这样的形式打印出组件掩码。
 * @param  out 要写入的流。  @param  mask 要写的掩码。  @return
 * 对第一个参数的引用。
 *
 *
 */
std::ostream &
operator<<(std::ostream &out, const ComponentMask &mask);

#ifndef DOXYGEN
// -------------------- inline functions ---------------------

inline ComponentMask::ComponentMask(const std::vector<bool> &component_mask)
  : component_mask(component_mask)
{}


inline ComponentMask::ComponentMask(const unsigned int n_components,
                                    const bool         initializer)
  : component_mask(n_components, initializer)
{}


inline unsigned int
ComponentMask::size() const
{
  return component_mask.size();
}


inline void
ComponentMask::set(const unsigned int index, const bool value)
{
  AssertIndexRange(index, component_mask.size());
  component_mask[index] = value;
}


inline bool ComponentMask::operator[](const unsigned int component_index) const
{
  // if the mask represents the all-component mask
  // then always return true
  if (component_mask.size() == 0)
    return true;
  else
    {
      // otherwise check the validity of the index and
      // return whatever is appropriate
      AssertIndexRange(component_index, component_mask.size());
      return component_mask[component_index];
    }
}


inline bool
ComponentMask::represents_n_components(const unsigned int n) const
{
  return ((component_mask.size() == 0) || (component_mask.size() == n));
}


inline unsigned int
ComponentMask::n_selected_components(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension(n, size());

  const unsigned int real_n = (n != numbers::invalid_unsigned_int ? n : size());
  if (component_mask.size() == 0)
    return real_n;
  else
    {
      AssertDimension(real_n, component_mask.size());
      return std::count_if(component_mask.begin(),
                           component_mask.end(),
                           [](const bool selected) { return selected; });
    }
}


inline unsigned int
ComponentMask::first_selected_component(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension(n, size());

  if (component_mask.size() == 0)
    return 0;
  else
    {
      for (unsigned int c = 0; c < component_mask.size(); ++c)
        if (component_mask[c] == true)
          return c;

      Assert(false, ExcMessage("No component is selected at all!"));
      return numbers::invalid_unsigned_int;
    }
}



inline bool
ComponentMask::represents_the_all_selected_mask() const
{
  return (component_mask.size() == 0);
}



inline ComponentMask
ComponentMask::operator|(const ComponentMask &mask) const
{
  // if one of the two masks denotes the all-component mask,
  // then return the other one
  if (component_mask.size() == 0)
    return mask;
  else if (mask.component_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(component_mask.size(), mask.component_mask.size());
      std::vector<bool> new_mask(component_mask.size());
      for (unsigned int i = 0; i < component_mask.size(); ++i)
        new_mask[i] = (component_mask[i] || mask.component_mask[i]);

      return new_mask;
    }
}


inline ComponentMask ComponentMask::operator&(const ComponentMask &mask) const
{
  // if one of the two masks denotes the all-component mask,
  // then return the other one
  if (component_mask.size() == 0)
    return mask;
  else if (mask.component_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(component_mask.size(), mask.component_mask.size());
      std::vector<bool> new_mask(component_mask.size());
      for (unsigned int i = 0; i < component_mask.size(); ++i)
        new_mask[i] = (component_mask[i] && mask.component_mask[i]);

      return new_mask;
    }
}


inline bool
ComponentMask::operator==(const ComponentMask &mask) const
{
  return component_mask == mask.component_mask;
}


inline bool
ComponentMask::operator!=(const ComponentMask &mask) const
{
  return component_mask != mask.component_mask;
}
#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif


