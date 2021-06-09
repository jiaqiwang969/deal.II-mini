//include/deal.II-translator/lac/block_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_block_vector_h
#define dealii_block_vector_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_type_traits.h>

#include <cstdio>
#include <vector>

DEAL_II_NAMESPACE_OPEN


// Forward declaration
#ifndef DOXYGEN
#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class BlockVector;
  }
} // namespace TrilinosWrappers
#  endif
#endif


/*!   @addtogroup  Vectors 
     * @{  

 
*
*/


/**
 * 一个基于deal.II向量的块向量的实现。虽然基类提供了大部分的接口，但这个类处理了向量的实际分配，并提供了底层向量类型的特定函数。
 *
 *
 * @note
 * 这个模板的实例提供给<tt>  @<float@>  和  @<double@></tt>;
 * 其他人可以在应用程序中生成（见手册中的 @ref
 * Instantiations 部分）。
 * @see   @ref GlossBlockLA  "块（线性代数）"
 *
 *
 */
template <typename Number>
class BlockVector : public BlockVectorBase<Vector<Number>>
{
public:
  /**
   * 对基类进行类型化定义，以便更简单地访问它自己的别名。
   *
   */
  using BaseClass = BlockVectorBase<Vector<Number>>;

  /**
   * 类型化底层向量的类型。
   *
   */
  using BlockType = typename BaseClass::BlockType;

  /**
   * 从基类中导入别名。
   *
   */
  using value_type      = typename BaseClass::value_type;
  using real_type       = typename BaseClass::real_type;
  using pointer         = typename BaseClass::pointer;
  using const_pointer   = typename BaseClass::const_pointer;
  using reference       = typename BaseClass::reference;
  using const_reference = typename BaseClass::const_reference;
  using size_type       = typename BaseClass::size_type;
  using iterator        = typename BaseClass::iterator;
  using const_iterator  = typename BaseClass::const_iterator;

  /**
   * 构造函数。有三种方法来使用这个构造函数。首先，没有任何参数，它生成一个没有块的对象。给定一个参数，它初始化<tt>n_blocks</tt>块，但这些块的大小为零。
   * 第三种变体最后将所有块初始化为相同的大小<tt>block_size</tt>。
   * 如果你打算使用不同大小的块，请参考下面的其他构造函数。
   *
   */
  explicit BlockVector(const unsigned int n_blocks   = 0,
                       const size_type    block_size = 0);

  /**
   * 复制构造函数。尺寸设置为 @p v,
   * 的尺寸，所有的组件都是从 @p v. 复制过来的。
   *
   */
  BlockVector(const BlockVector<Number> &V);


  /**
   * 移动构造函数。通过窃取给定参数向量的内部数据创建一个新的向量。
   *
   */
  BlockVector(BlockVector<Number> &&  /*v*/ ) noexcept = default;

  /**
   * 复制构造函数，获取另一种数据类型的BlockVector。如果没有从<tt>OtherNumber</tt>到<tt>Number</tt>的转换路径，这将失败。注意，当你复制到一个数据元素精度较低的BlockVector时，可能会失去精度。
   * 旧版本的gcc不尊重模板构造函数上的 @p explicit
   * 关键字。在这种情况下，很容易不小心写出效率很低的代码，因为编译器开始执行隐藏的转换。为了避免这种情况，如果我们在配置过程中检测到有损坏的编译器，这个功能就会被禁用。
   *
   */
  template <typename OtherNumber>
  explicit BlockVector(const BlockVector<OtherNumber> &v);

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 一个复制构造函数，获取一个（并行的）Trilinos块向量，并将其复制成deal.II自己的格式。
   *
   */
  BlockVector(const TrilinosWrappers::MPI::BlockVector &v);

#endif
  /**
   * 构造函数。设置块的数量为<tt>block_sizes.size()</tt>，并以<tt>block_sizes[i]</tt>零元素初始化每个块。
   *
   */
  BlockVector(const std::vector<size_type> &block_sizes);

  /**
   * 构造函数。将向量初始化为BlockIndices参数中的结构。
   *
   */
  BlockVector(const BlockIndices &block_indices);

  /**
   * 构造函数。设置块的数量为<tt>block_sizes.size()</tt>。
   * 用作为第二和第三参数的迭代器范围所指向的元素初始化向量。除了第一个参数外，这个构造函数与
   * <tt>std::vector</tt>
   * 类的相应构造函数完全类似，但需要第一个参数，以便知道如何将块向量细分为不同的块。
   *
   */
  template <typename InputIterator>
  BlockVector(const std::vector<size_type> &block_sizes,
              const InputIterator           first,
              const InputIterator           end);

  /**
   * 解除构造函数。清除内存
   *
   */
  ~BlockVector() override = default;

  /**
   * 在所有的子块上调用compress()函数。    这个功能只有在使用基于MPI的向量时才需要调用，为了兼容，在其他对象中也存在。    参见  @ref GlossCompress  "压缩分布式对象 "
   * 以获得更多信息。
   *
   */
  void
  compress(::dealii::VectorOperation::values operation =
             ::dealii::VectorOperation::unknown);

  /**
   * 返回`false`，因为这是一个串行块向量。
   * 只有在使用基于MPI的向量时才需要调用这个功能，并且为了兼容而存在于其他对象中。
   *
   */
  bool
  has_ghost_elements() const;

  /**
   * 复制操作：用给定的标量值填充向量的所有组件。
   *
   */
  BlockVector &
  operator=(const value_type s);

  /**
   * 对相同类型的参数进行复制操作。如果有必要，可以调整现在的向量的大小。
   *
   */
  BlockVector<Number> &
  operator=(const BlockVector<Number> &v);

  /**
   * 移动给定的向量。该操作符用给定参数向量的内容替换当前向量。
   *
   */
  BlockVector<Number> &
  operator=(BlockVector<Number> &&  /*v*/ ) = default; // NOLINT

  /**
   * 对不同类型的模板参数进行复制操作。如果有必要，可以调整当前向量的大小。
   *
   */
  template <class Number2>
  BlockVector<Number> &
  operator=(const BlockVector<Number2> &V);

  /**
   * 将一个常规向量复制到一个块向量中。
   *
   */
  BlockVector<Number> &
  operator=(const Vector<Number> &V);

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 从Trilinos块向量到deal.II块向量的复制构造函数。
   *
   */
  BlockVector<Number> &
  operator=(const TrilinosWrappers::MPI::BlockVector &V);
#endif

  /**
   * 重新初始化BlockVector，使其包含<tt>n_blocks</tt>每个大小为<tt>block_size</tt>的块。
   * 如果第二个参数是默认值，那么区块向量就会分配指定数量的区块，但让它们的大小为零。然后你需要重新初始化各个区块，并调用collect_sizes()来更新区块系统对其各个区块大小的认识。
   * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
   *
   */
  void
  reinit(const unsigned int n_blocks,
         const size_type    block_size           = 0,
         const bool         omit_zeroing_entries = false);

  /**
   * 重新初始化BlockVector，使其包含<tt>block_sizes.size()</tt>块。每个块都被重新初始化为<tt>block_sizes[i]</tt>的尺寸。
   * 如果块的数量与调用此函数前相同，所有的向量都保持不变，并且为每个向量调用
   * reinit()。
   * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
   * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新其缓存的向量大小。如果你在其中一个块上调用reinit()，那么对这个对象的后续操作可能会产生不可预测的结果，因为它们可能会被路由到错误的块上。
   *
   */
  void
  reinit(const std::vector<size_type> &block_sizes,
         const bool                    omit_zeroing_entries = false);

  /**
   * 重新初始化BlockVector以反映BlockIndices中发现的结构。
   * 如果块的数量与调用此函数前相同，所有的向量都保持不变，并且为每个向量调用
   * reinit()。
   * 如果<tt>omit_zeroing_entries==false</tt>，则向量被填充为零。
   *
   */
  void
  reinit(const BlockIndices &block_indices,
         const bool          omit_zeroing_entries = false);

  /**
   * 将维度改为向量<tt>V</tt>的维度。这与另一个reinit()函数同样适用。
   * <tt>V</tt>的元素不会被复制，也就是说，这个函数与调用<tt>reinit
   * (V.size(), omit_zeroing_entries)</tt>相同。
   * 注意，你必须调用这个（或其他reinit()函数）函数，而不是调用单个块的reinit()函数，以允许块向量更新它的向量大小缓存。如果你调用其中一个块的reinit()，那么这个对象的后续操作可能会产生不可预测的结果，因为它们可能被路由到错误的块。
   *
   */
  template <typename Number2>
  void
  reinit(const BlockVector<Number2> &V,
         const bool                  omit_zeroing_entries = false);

  /**
   * 将这个向量的每个元素乘以<tt>v</tt>的相应元素。
   *
   */
  template <class BlockVector2>
  void
  scale(const BlockVector2 &v);

  /**
   * 把这个向量的内容和另一个向量<tt>v</tt>交换。我们可以通过一个临时变量和复制数据元素来完成这个操作，但是这个函数的效率明显更高，因为它只交换了两个向量的数据指针，因此不需要分配临时存储空间和移动数据。
   * 这个函数类似于所有C++标准容器的swap()函数。此外，还有一个全局函数swap(u,v)，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
   *
   */
  void
  swap(BlockVector<Number> &v);

  /**
   * 打印到一个流。
   *
   */
  void
  print(std::ostream &     out,
        const unsigned int precision  = 3,
        const bool         scientific = true,
        const bool         across     = true) const;

  /**
   * 将向量全部写到流中。这是以二进制模式进行的，所以输出结果既不能被人类阅读，也不能（可能）被其他使用不同操作系统或数字格式的计算机阅读。
   *
   */
  void
  block_write(std::ostream &out) const;

  /**
   * 从一个文件中读取一个矢量en块。这是用上述函数的逆运算来完成的，所以它的速度相当快，因为位流没有被解释过。
   * 如果有必要，矢量会被调整大小。
   * 一个原始形式的错误检查被执行，它将识别最直截了当的尝试，将一些数据解释为存储在文件中的位流向量，但不会超过。
   *
   */
  void
  block_read(std::istream &in);

  /**
   * @addtogroup  Exceptions  @{
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException0(ExcIteratorRangeDoesNotMatchVectorSize);
  //@}
};

 /*@}*/ 

#ifndef DOXYGEN
 /*----------------------- Inline functions ----------------------------------*/ 



template <typename Number>
template <typename InputIterator>
BlockVector<Number>::BlockVector(const std::vector<size_type> &block_sizes,
                                 const InputIterator           first,
                                 const InputIterator           end)
{
  // first set sizes of blocks, but
  // don't initialize them as we will
  // copy elements soon
  (void)end;
  reinit(block_sizes, true);
  InputIterator start = first;
  for (size_type b = 0; b < block_sizes.size(); ++b)
    {
      InputIterator end = start;
      std::advance(end, static_cast<signed int>(block_sizes[b]));
      std::copy(start, end, this->block(b).begin());
      start = end;
    };
  Assert(start == end, ExcIteratorRangeDoesNotMatchVectorSize());
}



template <typename Number>
inline BlockVector<Number> &
BlockVector<Number>::operator=(const value_type s)
{
  AssertIsFinite(s);

  BaseClass::operator=(s);
  return *this;
}



template <typename Number>
inline BlockVector<Number> &
BlockVector<Number>::operator=(const BlockVector<Number> &v)
{
  reinit(v, true);
  BaseClass::operator=(v);
  return *this;
}



template <typename Number>
inline BlockVector<Number> &
BlockVector<Number>::operator=(const Vector<Number> &v)
{
  BaseClass::operator=(v);
  return *this;
}



template <typename Number>
template <typename Number2>
inline BlockVector<Number> &
BlockVector<Number>::operator=(const BlockVector<Number2> &v)
{
  reinit(v, true);
  BaseClass::operator=(v);
  return *this;
}

template <typename Number>
inline void
BlockVector<Number>::compress(::dealii::VectorOperation::values operation)
{
  for (size_type i = 0; i < this->n_blocks(); ++i)
    this->components[i].compress(operation);
}



template <typename Number>
inline bool
BlockVector<Number>::has_ghost_elements() const
{
  return false;
}



template <typename Number>
template <class BlockVector2>
void
BlockVector<Number>::scale(const BlockVector2 &v)
{
  BaseClass::scale(v);
}

#endif // DOXYGEN


/**
 * 全局函数，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  BlockVector
 *
 *
 */
template <typename Number>
inline void
swap(BlockVector<Number> &u, BlockVector<Number> &v)
{
  u.swap(v);
}


namespace internal
{
  namespace LinearOperatorImplementation
  {
    template <typename>
    class ReinitHelper;

    /**
     * linear_operator.h中内部使用的一个辅助类。对BlockVector<number>的特殊化。
     *
     */
    template <typename number>
    class ReinitHelper<BlockVector<number>>
    {
    public:
      template <typename Matrix>
      static void
      reinit_range_vector(const Matrix &       matrix,
                          BlockVector<number> &v,
                          bool                 omit_zeroing_entries)
      {
        v.reinit(matrix.get_row_indices(), omit_zeroing_entries);
      }

      template <typename Matrix>
      static void
      reinit_domain_vector(const Matrix &       matrix,
                           BlockVector<number> &v,
                           bool                 omit_zeroing_entries)
      {
        v.reinit(matrix.get_column_indices(), omit_zeroing_entries);
      }
    };

  } // namespace LinearOperatorImplementation
}  /* namespace internal */ 


/**
 * 将 dealii::BlockVector 声明为串行向量。
 *
 *
 */
template <typename Number>
struct is_serial_vector<BlockVector<Number>> : std::true_type
{};

DEAL_II_NAMESPACE_CLOSE

#endif


