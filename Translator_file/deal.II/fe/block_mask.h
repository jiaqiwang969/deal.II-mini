//include/deal.II-translator/fe/block_mask_0.txt
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

#ifndef dealii_fe_block_mask_h
#define dealii_fe_block_mask_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <algorithm>
#include <iosfwd>
#include <vector>

DEAL_II_NAMESPACE_OPEN



/**
 * 该类代表一个掩码，可用于选择有限元的单个矢量块（另见 @ref GlossBlockMask "该术语条目"
 * ）。它通常有和有限元块一样多的元素，人们可以使用
 * <code>operator[]</code> 来查询某个特定块是否被选中。
 * 该类的语义与相关的ComponentMask类相同，即一个默认构造的掩码代表所有可能的块。关于这些语义的更多信息见那里。
 * 这类对象在很多地方被使用，人们希望将操作限制在某个块的子集上，例如在 DoFTools::extract_dofs.
 * 中，这些对象可以手工创建，或者更简单，要求有限元从某些选定的块中生成一个块掩码，使用这样的代码，我们创建一个掩码，只表示一个斯托克斯元的速度块（见
 * @ref vector_valued  ）。
 *
 * @code
 * // Q2 element for the velocities, Q1 element for the pressure
 * FESystem<dim> stokes_fe (FESystem<dim>(FE_Q<dim>(2), dim), 1,
 *                          FE_Q<dim>(1),                     1);
 * FEValuesExtractors::Scalar pressure(dim);
 * BlockMask pressure_mask = stokes_fe.block_mask (pressure);
 * @endcode
 * 注意，通过将速度元素包装成一个FESystem对象，我们确保整个元素只有两个块。其结果是一个块掩码，在2D和3D中都有<code>[false,
 * true]</code>的值。(将其与ComponentMask文档中讨论的相应组件掩码进行比较)。同样地，使用
 *
 * @code
 * FEValuesExtractors::Vector velocities(0);
 * BlockMask velocity_mask = stokes_fe.block_mask (velocities);
 * @endcode
 * 将导致2d和3d中的掩码 <code>[true, false]</code> 。
 *
 *
 * @ingroup fe
 *
 * @ingroup vector_valued
 *
 */
class BlockMask
{
public:
  /**
   * 初始化一个区块屏蔽。默认情况下，一个区块掩码代表一组被<i>all</i>选中的区块，也就是说，调用这个构造函数的结果是一个区块掩码，每当询问一个区块是否被选中时，总是返回
   * <code>true</code> 。
   *
   */
  BlockMask() = default;

  /**
   * 用参数指定的一组选定的块来初始化这个类型的对象。
   * @param  block_mask 一个 <code>true/false</code>
   * 项的向量，决定有限元的哪些块被选择。如果给定向量的长度为零，则解释为<i>every</i>块被选中的情况。
   *
   */
  BlockMask(const std::vector<bool> &block_mask);

  /**
   * 用一些全部为真或假的元素来初始化块掩码。      @param
   * n_blocks 这个掩码的元素数量  @param  initializer
   * 这些元素中的每一个应该具有的值：要么是真要么是假。
   *
   */
  BlockMask(const unsigned int n_blocks, const bool initializer);

  /**
   * 如果此块掩码已被初始化为大小大于0的掩码，那么返回此对象所代表的掩码的大小。另一方面，如果这个掩码已被初始化为一个空对象，代表一个对每个元素都是真的掩码（即，如果这个对象在调用
   * represents_the_all_selected_mask()时将返回true），那么返回0，因为没有明确的大小。
   *
   */
  unsigned int
  size() const;

  /**
   * 返回一个特定的块是否被这个掩码所选择。如果这个掩码代表了选择<i>all
   * blocks</i>的对象的情况（例如，如果它是使用默认的构造函数创建的，或者是从bool类型的空向量转换而来的），那么无论给定的参数是什么，这个函数都返回true。
   * @param  block_index
   * 该函数应返回是否选择该块的索引。如果这个对象代表一个掩码，其中所有的块总是被选中，那么这里允许任何索引。否则，给定的索引需要在零和这个掩码所代表的块数之间。
   *
   */
  bool operator[](const unsigned int block_index) const;

  /**
   * 返回此块掩码是否正好代表 <code>n</code>
   * 块的掩码。如果它被初始化为一个正好有 <code>n</code>
   * entries of type <code>bool</code> 的向量（在这种情况下， @p n
   * 必须等于size()的结果），或者它被初始化为一个空向量（或使用默认构造函数），在这种情况下，它可以代表一个有任意数量块的掩码，并且将总是说一个块被选中，这就是真的。
   *
   */
  bool
  represents_n_blocks(const unsigned int n) const;

  /**
   * 返回这个掩码所选择的区块的数量。
   * 由于空的区块掩码代表一个区块掩码，每个区块都会返回
   * <code>true</code>
   * ，这个函数可能不知道区块掩码的真实大小，因此它需要一个参数来表示区块的总数量。
   * 如果该对象已被初始化为一个非空的掩码（即，如果size()函数返回大于0的东西，或者等价地，如果respresent_the_all_selected_mask()返回false），那么该参数可以被省略，而size()的结果被取。
   *
   */
  unsigned int
  n_selected_blocks(const unsigned int overall_number_of_blocks =
                      numbers::invalid_unsigned_int) const;

  /**
   * 返回第一个被选中的块的索引。该参数存在的原因与n_selected_blocks()函数存在的原因相同。
   * 如果根本没有选择任何块，该函数会抛出一个异常。
   *
   */
  unsigned int
  first_selected_block(const unsigned int overall_number_of_blocks =
                         numbers::invalid_unsigned_int) const;

  /**
   * 如果这个掩码代表一个默认构造的掩码，对应于所有块都被选中的掩码，则返回true。如果为真，那么size()函数将返回0。
   *
   */
  bool
  represents_the_all_selected_mask() const;

  /**
   * 返回一个包含由当前对象选择的块和作为参数传递的块的联合体的块掩码。
   *
   */
  BlockMask
  operator|(const BlockMask &mask) const;

  /**
   * 返回一个块掩码，它只包含那些在当前对象和作为参数传递的对象中都被设置的元素。
   *
   */
  BlockMask operator&(const BlockMask &mask) const;

  /**
   * 返回这个对象和参数是否相同。
   *
   */
  bool
  operator==(const BlockMask &mask) const;

  /**
   * 返回这个对象和参数是否不相同。
   *
   */
  bool
  operator!=(const BlockMask &mask) const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 实际的块掩码。
   *
   */
  std::vector<bool> block_mask;

  // make the output operator a friend so it can access
  // the block_mask array
  friend std::ostream &
  operator<<(std::ostream &out, const BlockMask &mask);
};


/**
 * 将块掩码写到输出流中。如果该块掩码代表一个所有块都被选中而没有指定掩码的特定大小，那么它将把字符串
 * <code>[all blocks selected]</code> 写到流中。否则，它将以
 * <code>[true,true,true,false]</code> 的形式打印块掩码。
 * @param  out 要写入的流。  @param  mask 要写入的掩码。  @return
 * 对第一个参数的引用。
 *
 *
 */
std::ostream &
operator<<(std::ostream &out, const BlockMask &mask);


#ifndef DOXYGEN
// -------------------- inline functions ---------------------

inline BlockMask::BlockMask(const std::vector<bool> &block_mask)
  : block_mask(block_mask)
{}


inline BlockMask::BlockMask(const unsigned int n_blocks, const bool initializer)
  : block_mask(n_blocks, initializer)
{}


inline unsigned int
BlockMask::size() const
{
  return block_mask.size();
}


inline bool BlockMask::operator[](const unsigned int block_index) const
{
  // if the mask represents the all-block mask
  // then always return true
  if (block_mask.size() == 0)
    return true;
  else
    {
      // otherwise check the validity of the index and
      // return whatever is appropriate
      AssertIndexRange(block_index, block_mask.size());
      return block_mask[block_index];
    }
}


inline bool
BlockMask::represents_n_blocks(const unsigned int n) const
{
  return ((block_mask.size() == 0) || (block_mask.size() == n));
}


inline unsigned int
BlockMask::n_selected_blocks(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension(n, size());

  const unsigned int real_n = (n != numbers::invalid_unsigned_int ? n : size());
  if (block_mask.size() == 0)
    return real_n;
  else
    {
      AssertDimension(real_n, block_mask.size());
      return std::count_if(block_mask.begin(),
                           block_mask.end(),
                           [](const bool selected) { return selected; });
    }
}


inline unsigned int
BlockMask::first_selected_block(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension(n, size());

  if (block_mask.size() == 0)
    return 0;
  else
    {
      for (unsigned int c = 0; c < block_mask.size(); ++c)
        if (block_mask[c] == true)
          return c;

      Assert(false, ExcMessage("No block is selected at all!"));
      return numbers::invalid_unsigned_int;
    }
}



inline bool
BlockMask::represents_the_all_selected_mask() const
{
  return (block_mask.size() == 0);
}



inline BlockMask
BlockMask::operator|(const BlockMask &mask) const
{
  // if one of the two masks denotes the all-block mask,
  // then return the other one
  if (block_mask.size() == 0)
    return mask;
  else if (mask.block_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(block_mask.size(), mask.block_mask.size());
      std::vector<bool> new_mask(block_mask.size());
      for (unsigned int i = 0; i < block_mask.size(); ++i)
        new_mask[i] = (block_mask[i] || mask.block_mask[i]);

      return new_mask;
    }
}


inline BlockMask BlockMask::operator&(const BlockMask &mask) const
{
  // if one of the two masks denotes the all-block mask,
  // then return the other one
  if (block_mask.size() == 0)
    return mask;
  else if (mask.block_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(block_mask.size(), mask.block_mask.size());
      std::vector<bool> new_mask(block_mask.size());
      for (unsigned int i = 0; i < block_mask.size(); ++i)
        new_mask[i] = (block_mask[i] && mask.block_mask[i]);

      return new_mask;
    }
}


inline bool
BlockMask::operator==(const BlockMask &mask) const
{
  return block_mask == mask.block_mask;
}


inline bool
BlockMask::operator!=(const BlockMask &mask) const
{
  return block_mask != mask.block_mask;
}
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif


