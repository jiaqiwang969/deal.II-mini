//include/deal.II-translator/lac/block_indices_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2020 by the deal.II authors
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

#ifndef dealii_block_indices_h
#define dealii_block_indices_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/utilities.h>

#include <cstddef>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * BlockIndices表示一个指数范围（如一个矢量元素的有效指数范围
 * $[0,N)$ ），以及这一个范围如何被分解成更小但连续的
 * "块"（如一个解矢量的速度和压力部分）。特别是，它提供了在全局指数和一个块的指数<i>within</i>之间转换的能力。例如，该类用于BlockVector、BlockSparsityPattern和BlockMatrixBase类中。
 * 可以从这个类中获得的信息分为两组。首先，可以查询索引空间的全局大小（通过total_size()成员函数），以及块的数量和它们的大小（通过size()和block_size()函数）。
 * 其次，这个类管理全局索引到这个块内的局部索引的转换，以及反过来的转换。这是必要的，例如，当你在一个块向量中寻址一个全局元素，并想知道它在哪个块中，以及它在这个块中对应的索引。如果一个矩阵是由几个块组成的，你必须把全局的行和列索引翻译成局部的，这也很有用。
 *
 *
 * @ingroup data   @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
class BlockIndices : public Subscriptor
{
public:
  /**
   * 声明容器大小的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 默认构造函数。对零块进行初始化。
   *
   */
  BlockIndices();

  /**
   * 构造函数。初始化每个块 @p i
   * 中的条目数为<tt>block_sizes[i]</tt>。块的数量将是 @p
   * block_sizes的大小。
   *
   */
  BlockIndices(const std::vector<size_type> &block_sizes);

  /**
   * 移动构造函数。通过窃取另一个BlockIndices对象的内部数据来初始化一个新对象。
   *
   */
  BlockIndices(BlockIndices &&b) noexcept;

  /**
   * 复制构造函数。
   *
   */
  BlockIndices(const BlockIndices &) = default;

  /**
   * 用于具有相同大小块的结构的专门构造函数。
   *
   */
  explicit BlockIndices(const unsigned int n_blocks,
                        const size_type    block_size = 0);

  /**
   * 重新初始化块的数量，并给每个块分配相同数量的元素。
   *
   */
  void
  reinit(const unsigned int n_blocks, const size_type n_elements_per_block);

  /**
   * 根据给定的参数重新初始化每个块内的索引数。块的数量将被调整为<tt>block_sizes</tt>的大小，块
   * @p i 的大小被设置为<tt>block_sizes[i]</tt>。
   *
   */
  void
  reinit(const std::vector<size_type> &block_sizes);

  /**
   * 在块结构的末端添加另一个给定大小的块。
   *
   */
  void
  push_back(const size_type size);

  /**
   * @name  大小信息
   *
   */
  //@{

  /**
   * 索引字段中的块的数量。
   *
   */
  unsigned int
  size() const;

  /**
   * 返回所有块上累积的索引总数，也就是块向量空间的维度。
   *
   */
  size_type
  total_size() const;

  /**
   * @p ith 块的大小。
   *
   */
  size_type
  block_size(const unsigned int i) const;

  /**
   * 块大小的字符串表示。输出的形式是`[nb->b1,b2,b3|s]`，其中`nb`是n_blocks()，`s`是total_size()，`b1`等是block_size()对每个块返回的值。
   *
   */
  std::string
  to_string() const;

  //@}

  /**
   * @name  索引转换 本组函数假定一个对象，它是在按块排序后创建的，这样，每个块在对象中形成一组连续的索引。如果应用于其他对象，从这些函数得到的数字是没有意义的。
   *
   */
  //@{

  /**
   * 返回全局索引 @p
   * i的块和该块中的索引。这一对的第一个元素是块，第二个是其中的索引。
   *
   */
  std::pair<unsigned int, size_type>
  global_to_local(const size_type i) const;

  /**
   * 返回 @p block. 块中 @p index 的全局索引。
   *
   */
  size_type
  local_to_global(const unsigned int block, const size_type index) const;

  /**
   * 第1个区块的起始索引。
   *
   */
  size_type
  block_start(const unsigned int i) const;
  //@}

  /**
   * 复制操作者。
   *
   */
  BlockIndices &
  operator=(const BlockIndices &b);

  /**
   * 移动赋值运算符。通过转移BlockIndices对象的内容，将另一个BlockIndices对象移到当前对象上。
   *
   */
  BlockIndices &
  operator=(BlockIndices &&) noexcept;

  /**
   * 比较两个对象是否相同，即块的数量和所有块的大小是否相等。
   *
   */
  bool
  operator==(const BlockIndices &b) const;

  /**
   * 交换这两个对象的内容。
   *
   */
  void
  swap(BlockIndices &b);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * 块的数量。虽然这个值可以通过<tt>start_indices.size()-1</tt>获得，但我们对这个值进行了缓存以加快访问速度。
   *
   */
  unsigned int n_blocks;

  /**
   * 每个向量的全球起始索引。最后一个也是多余的值是总条目数。
   *
   */
  std::vector<size_type> start_indices;
};


/**
 * 记录BlockIndices的操作员。写下块的数量，每个块的大小和索引域的总大小。
 * @ref BlockIndices
 *
 *
 */
inline LogStream &
operator<<(LogStream &s, const BlockIndices &bi)
{
  const unsigned int n = bi.size();
  s << n << ":[";
  // Write first size without leading space
  if (n > 0)
    s << bi.block_size(0);
  // Write all other sizes
  for (unsigned int i = 1; i < n; ++i)
    s << ' ' << bi.block_size(i);
  s << "]->" << bi.total_size();
  return s;
}



 /* ---------------------- template and inline functions ------------------- */ 

inline void
BlockIndices::reinit(const unsigned int nb, const size_type block_size)
{
  n_blocks = nb;
  start_indices.resize(n_blocks + 1);
  for (size_type i = 0; i <= n_blocks; ++i)
    start_indices[i] = i * block_size;
}



inline void
BlockIndices::reinit(const std::vector<size_type> &block_sizes)
{
  if (start_indices.size() != block_sizes.size() + 1)
    {
      n_blocks = static_cast<unsigned int>(block_sizes.size());
      start_indices.resize(n_blocks + 1);
    }
  start_indices[0] = 0;
  for (size_type i = 1; i <= n_blocks; ++i)
    start_indices[i] = start_indices[i - 1] + block_sizes[i - 1];
}


inline BlockIndices::BlockIndices()
  : n_blocks(0)
  , start_indices(1, 0)
{}



inline BlockIndices::BlockIndices(const unsigned int n_blocks,
                                  const size_type    block_size)
  : n_blocks(n_blocks)
  , start_indices(n_blocks + 1)
{
  for (size_type i = 0; i <= n_blocks; ++i)
    start_indices[i] = i * block_size;
}



inline BlockIndices::BlockIndices(const std::vector<size_type> &block_sizes)
  : n_blocks(static_cast<unsigned int>(block_sizes.size()))
  , start_indices(block_sizes.size() + 1)
{
  reinit(block_sizes);
}



inline BlockIndices::BlockIndices(BlockIndices &&b) noexcept
  : n_blocks(b.n_blocks)
  , start_indices(std::move(b.start_indices))
{
  b.n_blocks      = 0;
  b.start_indices = std::vector<size_type>(1, 0);
}



inline void
BlockIndices::push_back(const size_type sz)
{
  start_indices.push_back(start_indices[n_blocks] + sz);
  ++n_blocks;
  AssertDimension(start_indices.size(), n_blocks + 1);
}


inline std::pair<unsigned int, BlockIndices::size_type>
BlockIndices::global_to_local(const size_type i) const
{
  AssertIndexRange(i, total_size());
  Assert(n_blocks > 0, ExcLowerRangeType<size_type>(i, size_type(1)));

  // start_indices[0] == 0 so we might as well start from the next one
  const auto it =
    --std::upper_bound(++start_indices.begin(), start_indices.end(), i);

  return {std::distance(start_indices.begin(), it), i - *it};
}


inline BlockIndices::size_type
BlockIndices::local_to_global(const unsigned int block,
                              const size_type    index) const
{
  AssertIndexRange(block, n_blocks);
  AssertIndexRange(index, start_indices[block + 1] - start_indices[block]);

  return start_indices[block] + index;
}


inline unsigned int
BlockIndices::size() const
{
  return n_blocks;
}



inline BlockIndices::size_type
BlockIndices::total_size() const
{
  if (n_blocks == 0)
    return 0;
  return start_indices[n_blocks];
}



inline BlockIndices::size_type
BlockIndices::block_size(const unsigned int block) const
{
  AssertIndexRange(block, n_blocks);
  return start_indices[block + 1] - start_indices[block];
}



inline std::string
BlockIndices::to_string() const
{
  std::string result = "[" + Utilities::int_to_string(n_blocks) + "->";
  for (unsigned int i = 0; i < n_blocks; ++i)
    {
      if (i > 0)
        result += ',';
      result += std::to_string(block_size(i));
    }
  result += "|" + std::to_string(total_size()) + ']';
  return result;
}



inline BlockIndices::size_type
BlockIndices::block_start(const unsigned int block) const
{
  AssertIndexRange(block, n_blocks);
  return start_indices[block];
}



inline BlockIndices &
BlockIndices::operator=(const BlockIndices &b)
{
  start_indices = b.start_indices;
  n_blocks      = b.n_blocks;
  return *this;
}



inline BlockIndices &
BlockIndices::operator=(BlockIndices &&b) noexcept
{
  start_indices = std::move(b.start_indices);
  n_blocks      = b.n_blocks;

  b.start_indices = std::vector<size_type>(1, 0);
  b.n_blocks      = 0;

  return *this;
}



inline bool
BlockIndices::operator==(const BlockIndices &b) const
{
  if (n_blocks != b.n_blocks)
    return false;

  for (size_type i = 0; i <= n_blocks; ++i)
    if (start_indices[i] != b.start_indices[i])
      return false;

  return true;
}



inline void
BlockIndices::swap(BlockIndices &b)
{
  std::swap(n_blocks, b.n_blocks);
  std::swap(start_indices, b.start_indices);
}



inline std::size_t
BlockIndices::memory_consumption() const
{
  return (sizeof(*this) + start_indices.size() * sizeof(start_indices[0]));
}



 /* ----------------- global functions ---------------------------- */ 


/**
 * 全局函数  @p swap
 * ，它重载了C++标准库中使用临时对象的默认实现。该函数简单地交换了两个对象的数据。
 * @relatesalso  BlockIndices
 *
 *
 */
inline void
swap(BlockIndices &u, BlockIndices &v)
{
  u.swap(v);
}



DEAL_II_NAMESPACE_CLOSE

#endif


