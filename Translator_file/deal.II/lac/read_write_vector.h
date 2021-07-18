//include/deal.II-translator/lac/read_write_vector_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II authors
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

#ifndef dealii_read_write_vector_h
#define dealii_read_write_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/communication_pattern_base.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector_operation.h>

#include <cstdlib>
#include <cstring>
#include <iomanip>

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_epetra_communication_pattern.h>
#  include <deal.II/lac/trilinos_epetra_vector.h>
#  include <deal.II/lac/trilinos_tpetra_vector.h>

#  include <Epetra_MultiVector.h>

#endif

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
namespace LinearAlgebra
{
  template <typename>
  class Vector;
  namespace distributed
  {
    template <typename, typename>
    class Vector;
  } // namespace distributed
} // namespace LinearAlgebra

#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  namespace MPI
  {
    class Vector;
  }
} // namespace TrilinosWrappers
#  endif

#  ifdef DEAL_II_WITH_CUDA
namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    template <typename>
    class Vector;
  }
} // namespace LinearAlgebra
#  endif
#endif

namespace LinearAlgebra
{
  /*!   @addtogroup Vectors  
     * @{   
*
*/

  /**
   * ReadWriteVector旨在表示 ${\mathbb R}^N$
   * 中的向量，对于这些向量，它存储了所有或一个子集的元素。后一种情况在并行计算中很重要，因为
   * $N$
   * 可能大到没有一个处理器能真正解决向量的所有元素，但这也不是必须的：通常我们只需要存储住在本地拥有的单元上的自由度值，加上可能住在幽灵单元上的自由度。
   * 该类允许访问单个元素的读或写。
   * 然而，它不允许全局性的操作，如取法线。
   * ReadWriteVector可以用来读写从VectorSpaceVector派生的向量中的元素，如
   * TrilinosWrappers::MPI::Vector 和 PETScWrappers::MPI::Vector. <h3>Storing
   * elements</h3>，大多数情况下，人们将简单地使用这些自由度的全局数字从当前类的向量中读出或写入。这是用operator()()或operator[]()来完成的，它们调用global_to_local()将<i>global</i>的索引转化为<i>local</i>的。在这种情况下，很明显，人们只能访问当前对象确实存储的向量中的元素。
   * 然而，也可以按照当前对象所存储的元素的顺序来访问这些元素。换句话说，人们对用<i>global</i>的索引访问元素不感兴趣，而是使用一个枚举，只考虑实际存储的元素。local_element()函数为此提供了便利。为此，有必要知道<i>in
   * which
   * order</i>当前类存储的元素。所有连续范围的元素都按照每个范围的第一个索引的升序来存储。
   * 可以用函数 IndexSet::largest_range_starting_index()
   * 来获取最大范围的第一个索引。
   *
   */
  template <typename Number>
  class ReadWriteVector : public Subscriptor
  {
  public:
    /**
     * 声明所有容器中使用的标准类型。这些类型与<tt>C++</tt>标准库中的<tt>vector<...></tt>类中的类型平行。
     *
     */
    using value_type      = Number;
    using pointer         = value_type *;
    using const_pointer   = const value_type *;
    using iterator        = value_type *;
    using const_iterator  = const value_type *;
    using reference       = value_type &;
    using const_reference = const value_type &;
    using size_type       = types::global_dof_index;
    using real_type       = typename numbers::NumberTraits<Number>::real_type;

    /**
     * @name  1: 基本对象处理
     *
     */
    //@{
    /**
     * 空的构造函数。
     *
     */
    ReadWriteVector();

    /**
     * 复制构造函数。
     *
     */
    ReadWriteVector(const ReadWriteVector<Number> &in_vector);

    /**
     * 构建一个给定大小的向量，存储的元素的索引在[0,size]。
     *
     */
    explicit ReadWriteVector(const size_type size);

    /**
     * 构建一个向量，其存储元素的索引由IndexSet @p
     * locally_stored_indices. 给出。
     *
     */
    explicit ReadWriteVector(const IndexSet &locally_stored_indices);

    /**
     * 解构器。
     *
     */
    ~ReadWriteVector() override = default;

    /**
     * 将向量的全局大小设置为 @p size.
     * 存储元素的索引在[0,size]。        如果标志 @p
     * omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用它之前必须确保用合理的数据填充它）。
     *
     */
    virtual void
    reinit(const size_type size, const bool omit_zeroing_entries = false);

    /**
     * 使用与输入向量 @p in_vector
     * 相同的IndexSet，并为该向量分配内存。        如果标志
     * @p omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用前必须确保用合理的数据填充它）。
     *
     */
    template <typename Number2>
    void
    reinit(const ReadWriteVector<Number2> &in_vector,
           const bool                      omit_zeroing_entries = false);

    /**
     * 初始化向量。指数由 @p  locally_stored_indices指定。
     * 如果标志 @p omit_zeroing_entries
     * 被设置为false，内存将被初始化为0，否则内存将不被触动（用户在使用前必须确保用合理的数据填充它）。
     * local_stored_indices。
     *
     */
    virtual void
    reinit(const IndexSet &locally_stored_indices,
           const bool      omit_zeroing_entries = false);


#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_WITH_MPI
    /**
     * 通过提供对给定的ghosted或non-ghosted向量中所有本地可用条目的访问，初始化这个读写向量。
     * @note
     * 这个函数目前将参数中的值复制到ReadWriteVector中，所以这里的修改不会修改
     * @p trilinos_vec.
     * 这个函数主要是为了向后兼容而写的，以获得对库内重影
     * TrilinosWrappers::MPI::Vector 的元素访问。
     *
     */
    void
    reinit(const TrilinosWrappers::MPI::Vector &trilinos_vec);
#  endif
#endif

    /**
     * 对向量的每个元素应用漏斗 @p func
     * 。这个函式应该看起来像
     * @code
     * struct Functor
     * {
     * void operator() (Number &value);
     * };
     * @endcode
     * @note
     * 这个函数要求包含头文件read_write_vector.templates.h。
     *
     */
    template <typename Functor>
    void
    apply(const Functor &func);

    /**
     * 交换这个向量和另一个向量的内容  @p v.
     * 人们可以用一个临时变量和复制数据元素来完成这个操作，但是这个函数明显更有效率，因为它只交换了两个向量的数据指针，因此不需要分配临时存储和移动数据。
     * 这个函数类似于所有C++标准容器的 @p swap
     * 函数。此外，还有一个全局函数<tt>swap(u,v)</tt>，它简单地调用<tt>u.swap(v)</tt>，同样与标准函数相类似。
     *
     */
    void
    swap(ReadWriteVector<Number> &v);

    /**
     * 复制数据和输入向量的IndexSet  @p in_vector.  。
     *
     */
    ReadWriteVector<Number> &
    operator=(const ReadWriteVector<Number> &in_vector);

    /**
     * 复制数据和输入向量的IndexSet  @p in_vector.  。
     *
     */
    template <typename Number2>
    ReadWriteVector<Number> &
    operator=(const ReadWriteVector<Number2> &in_vector);

    /**
     * 将向量的所有元素设置为标量  @p s.  只有在 @p s
     * 等于零的情况下才允许此操作。
     *
     */
    ReadWriteVector<Number> &
    operator=(const Number s);

    /**
     * 从输入向量 @p vec.   VectorOperation::values   @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。
     * @note  参数 @p communication_pattern
     * 被忽略，因为我们在这里处理的是一个串行向量。
     *
     */
    void
    import(const dealii::Vector<Number> &vec,
           VectorOperation::values       operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});

    /**
     * 从输入向量 @p vec.   VectorOperation::values   @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。
     * @note 参数 @p communication_pattern
     * 被忽略，因为我们在这里处理的是一个串行向量。
     *
     */
    void
    import(const LinearAlgebra::Vector<Number> &vec,
           VectorOperation::values              operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});

    /**
     * 从输入向量 @p vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     *
     */
    template <typename MemorySpace>
    void
    import(const distributed::Vector<Number, MemorySpace> &vec,
           VectorOperation::values                         operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});

#ifdef DEAL_II_WITH_PETSC
    /**
     * 从输入向量 @p petsc_vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该添加到当前向量中或替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     *
     */
    void
    import(const PETScWrappers::MPI::Vector &petsc_vec,
           VectorOperation::values           operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});
#endif

#ifdef DEAL_II_WITH_TRILINOS
    /**
     * 从输入向量 @p trilinos_vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     * @note  @p trilinos_vec 不允许有幽灵条目。
     *
     */
    void
    import(const TrilinosWrappers::MPI::Vector &trilinos_vec,
           VectorOperation::values              operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});

#  ifdef DEAL_II_WITH_MPI
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * 从输入向量 @p tpetra_vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     *
     */
    void
    import(const TpetraWrappers::Vector<Number> &tpetra_vec,
           VectorOperation::values               operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});
#    endif

    /**
     * 从输入向量 @p epetra_vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。如果多次使用同一通信模式，可以使用最后一个参数。这可以用来提高性能。
     *
     */
    void
    import(const EpetraWrappers::Vector &epetra_vec,
           VectorOperation::values       operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});
#  endif
#endif

#ifdef DEAL_II_WITH_CUDA
    /**
     * 从输入向量 @p cuda_vec.  VectorOperation::values  @p operation
     * 中导入向量的IndexSet中存在的所有元素，用于决定 @p V
     * 中的元素是否应该被添加到当前向量中，或者替换当前元素。最后一个参数不使用。
     *
     */
    void
    import(const CUDAWrappers::Vector<Number> &cuda_vec,
           VectorOperation::values             operation,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern = {});
#endif

    /**
     * 这个函数返回的值表示这类对象所模拟的向量空间的维度。然而，当前类的对象实际上并不存储这个空间的向量的所有元素，而实际上可能只存储一个子集。存储的元素数量由n_elements()返回，并且小于或等于当前函数返回的数量。
     *
     */
    size_type
    size() const;

    /**
     * 这个函数返回存储的元素数。它小于或等于这种对象所模拟的向量空间的维度。这个维度是由size()返回的。
     * @deprecated  使用local_owned_size()代替。
     *
     */
    DEAL_II_DEPRECATED
    size_type
    n_elements() const;

    /**
     * 返回向量的本地大小，即本地拥有的索引数。
     *
     */
    size_type
    locally_owned_size() const;

    /**
     * 返回IndexSet，代表存储的元素的索引。
     *
     */
    const IndexSet &
    get_stored_elements() const;

    /**
     * 使 @p ReadWriteVector
     * 类有点像C++标准库的<tt>vector<>/tt>类，返回该向量的<i>locally
     * stored</i>元素的开始和结束的迭代器。
     *
     */
    iterator
    begin();

    /**
     * 返回本地存储的向量元素的开始的常数迭代器。
     *
     */
    const_iterator
    begin() const;

    /**
     * 返回一个迭代器，指向本地存储的条目数组结束后的元素。
     *
     */
    iterator
    end();

    /**
     * 返回一个恒定的迭代器，指向本地存储的条目数组结束后的元素。
     *
     */
    const_iterator
    end() const;
    //@}


    /**
     * @name  2: 数据访问
     *
     */
    //@{

    /**
     * 读取对应于 @p  global_index位置的数据。如果 @p
     * global_index
     * 不是由当前对象存储的，就会产生一个异常。
     *
     */
    Number
    operator()(const size_type global_index) const;

    /**
     * 对 @p  global_index对应位置的数据进行读写访问。如果 @p
     * global_index 没有被当前对象存储，则会产生一个异常。
     *
     */
    Number &
    operator()(const size_type global_index);

    /**
     * 读取对应于 @p  global_index位置的数据。如果 @p
     * global_index 没有被当前对象存储，就会产生一个异常。
     * 这个函数与operator()的作用相同。
     *
     */
    Number operator[](const size_type global_index) const;

    /**
     * 读取和写入对应于 @p global_index位置的数据。如果 @p
     * global_index
     * 不是由当前对象存储的，就会抛出一个异常。
     * 这个函数与operator()的作用相同。
     *
     */
    Number &operator[](const size_type global_index);

    /**
     * 与通过operator()获取向量中的单个元素不同，这个函数允许一次获取整个元素集。要读取的元素的索引在第一个参数中说明，相应的值在第二个参数中返回。
     * 如果当前的向量被称为 @p v,
     * ，那么这个函数就等同于代码
     * @code
     * for (unsigned int i=0; i<indices.size(); ++i)
     *   values[i] = v[indices[i]];
     * @endcode
     * @pre  @p indices 和 @p values 数组的大小必须是一致的。
     *
     */
    template <typename Number2>
    void
    extract_subvector_to(const std::vector<size_type> &indices,
                         std::vector<Number2> &        values) const;

    /**
     * 这个函数不是通过operator()获得向量的单个元素，而是允许一次获得整个元素集。与前一个函数不同的是，这个函数通过取消引用前两个参数提供的迭代器范围内的所有元素来获得元素的索引，并将向量的值放入通过取消引用从第三个参数指向的位置开始的迭代器范围获得的内存位置。
     * 如果当前的向量被称为 @p v,
     * ，那么这个函数就等同于代码
     * @code
     * ForwardIterator indices_p = indices_begin;
     * OutputIterator  values_p  = values_begin;
     * while (indices_p != indices_end)
     * {
     *  values_p = v[*indices_p];
     *   ++indices_p;
     *   ++values_p;
     * }
     * @endcode
     * @pre  必须能够写进从 @p values_begin
     * 开始的尽可能多的内存位置，因为有 @p indices_begin 和
     * @p indices_end. 之间的迭代器。
     *
     */
    template <typename ForwardIterator, typename OutputIterator>
    void
    extract_subvector_to(ForwardIterator       indices_begin,
                         const ForwardIterator indices_end,
                         OutputIterator        values_begin) const;

    /**
     * 读取访问 @p local_index. 指定的数据字段
     * 当你按照元素的存储顺序访问它们时，你有必要知道它们被存储在哪个位置。换句话说，你需要知道这个类所存储的元素的全局索引与这些全局元素的连续数组中的局部索引之间的映射。关于这一点，请看这个类的一般文档。
     * 性能。直接数组访问（快速）。
     *
     */
    Number
    local_element(const size_type local_index) const;

    /**
     * 对 @p local_index. 指定的数据域的读写访问
     * 当你按照元素的存储顺序访问它们时，你有必要知道它们是以何种方式存储的。换句话说，你需要知道这个类所存储的元素的全局索引与这些全局元素的连续数组中的局部索引之间的映射。关于这一点，请看这个类的一般文档。
     * 性能。直接数组访问（快速）。
     *
     */
    Number &
    local_element(const size_type local_index);
    //@}


    /**
     * @name  3: 修改向量
     *
     */
    //@{

    /**
     * 该函数将存储在 @p values 中的整组值添加到 @p indices.
     * 指定的向量成分中。
     *
     */
    template <typename Number2>
    void
    add(const std::vector<size_type> &indices,
        const std::vector<Number2> &  values);

    /**
     * 这个函数与前一个函数类似，但需要一个读写向量的值。
     *
     */
    template <typename Number2>
    void
    add(const std::vector<size_type> &  indices,
        const ReadWriteVector<Number2> &values);

    /**
     * 取一个<tt>n_elements</tt>连续存储的地址，并将其加入到向量中。处理所有其他两个<tt>add()</tt>函数没有涵盖的情况。
     *
     */
    template <typename Number2>
    void
    add(const size_type  n_elements,
        const size_type *indices,
        const Number2 *  values);

    /**
     * 将向量打印到输出流中  @p out.  。
     *
     */
    void
    print(std::ostream &     out,
          const unsigned int precision  = 3,
          const bool         scientific = true) const;

    /**
     * 以字节为单位返回这个类的内存消耗。
     *
     */
    std::size_t
    memory_consumption() const;
    //@}

  protected:
#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * 从输入向量中导入所有存在于向量IndexSet中的元素  @p
     * tpetra_vector.  这是一个辅助函数，不应该直接使用它。
     *
     */
    void
    import(
      const Tpetra::Vector<Number, int, types::global_dof_index> &tpetra_vector,
      const IndexSet &        locally_owned_elements,
      VectorOperation::values operation,
      const MPI_Comm &        mpi_comm,
      const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
        &communication_pattern);
#  endif

    /**
     * 从输入向量导入向量的IndexSet中的所有元素  @p
     * multivector.  这是一个辅助函数，不应该直接使用。
     *
     */
    void
    import(const Epetra_MultiVector &multivector,
           const IndexSet &          locally_owned_elements,
           VectorOperation::values   operation,
           const MPI_Comm &          mpi_comm,
           const std::shared_ptr<const Utilities::MPI::CommunicationPatternBase>
             &communication_pattern);
#endif

    /**
     * 返回 @p global_index. 的本地位置。
     *
     */
    unsigned int
    global_to_local(const types::global_dof_index global_index) const
    {
      // the following will throw an exception if the global_index is not
      // in the remaining_elements
      return static_cast<unsigned int>(
        stored_elements.index_within_set(global_index));
    }

    /**
     * 一个辅助函数，用于调整val数组的大小。
     *
     */
    void
    resize_val(const size_type new_allocated_size);

#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
    /**
     * 返回一个 TpetraWrappers::CommunicationPattern
     * ，并将其储存起来供将来使用。
     *
     */
    TpetraWrappers::CommunicationPattern
    create_tpetra_comm_pattern(const IndexSet &source_index_set,
                               const MPI_Comm &mpi_comm);
#  endif

    /**
     * 返回一个 EpetraWrappers::CommunicationPattern
     * ，并存储起来供将来使用。
     *
     */
    EpetraWrappers::CommunicationPattern
    create_epetra_comm_pattern(const IndexSet &source_index_set,
                               const MPI_Comm &mpi_comm);
#endif

    /**
     * 存储的元素的索引。
     *
     */
    IndexSet stored_elements;

    /**
     * 最后导入的向量的元素的索引集。
     *
     */
    IndexSet source_stored_elements;

    /**
     * Source_stored_elements IndexSet和当前向量之间的通信模式。
     *
     */
    std::shared_ptr<Utilities::MPI::CommunicationPatternBase> comm_pattern;

    /**
     * 指向这个向量的本地元素数组的指针。
     *
     */
    std::unique_ptr<Number[], decltype(std::free) *> values;

    /**
     * 对于带有TBB的并行循环，这个成员变量存储循环的亲和力信息。
     *
     */
    mutable std::shared_ptr<::dealii::parallel::internal::TBBPartitioner>
      thread_loop_partitioner;

    // Make all other ReadWriteVector types friends.
    template <typename Number2>
    friend class ReadWriteVector;

  private:
    /**
     * 这个类提供了一个围绕Functor的包装器，它作用于向量的单个元素。这对于使用
     * tbb::parallel_for 是必要的，因为它需要一个TBBForFunctor。
     *
     */
    template <typename Functor>
    class FunctorTemplate
    {
    public:
      /**
       * 构造器。取一个函数并存储它的副本。
       *
       */
      FunctorTemplate(ReadWriteVector<Number> &parent, const Functor &functor);

      /**
       * 用存储的functor的拷贝来评估元素。
       *
       */
      virtual void
      operator()(const size_type begin, const size_type end);

    private:
      /**
       * 对拥有FunctorTemplate的ReadWriteVector对象的别名。
       *
       */
      ReadWriteVector &parent;

      /**
       * 漏斗的副本。
       *
       */
      const Functor &functor;
    };
  };

   /*@}*/ 


   /*---------------------------- Inline functions ---------------------------*/ 

#ifndef DOXYGEN

  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector()
    : Subscriptor()
    , values(nullptr, free)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(0, true);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(
    const ReadWriteVector<Number> &v)
    : Subscriptor()
    , values(nullptr, free)
  {
    this->operator=(v);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(const size_type size)
    : Subscriptor()
    , values(nullptr, free)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(size, false);
  }



  template <typename Number>
  inline ReadWriteVector<Number>::ReadWriteVector(
    const IndexSet &locally_stored_indices)
    : Subscriptor()
    , values(nullptr, free)
  {
    // virtual functions called in constructors and destructors never use the
    // override in a derived class
    // for clarity be explicit on which function is called
    ReadWriteVector<Number>::reinit(locally_stored_indices);
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::size() const
  {
    return stored_elements.size();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::n_elements() const
  {
    return stored_elements.n_elements();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::size_type
  ReadWriteVector<Number>::locally_owned_size() const
  {
    return stored_elements.n_elements();
  }



  template <typename Number>
  inline const IndexSet &
  ReadWriteVector<Number>::get_stored_elements() const
  {
    return stored_elements;
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::begin()
  {
    return values.get();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::begin() const
  {
    return values.get();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::iterator
  ReadWriteVector<Number>::end()
  {
    return values.get() + this->n_elements();
  }



  template <typename Number>
  inline typename ReadWriteVector<Number>::const_iterator
  ReadWriteVector<Number>::end() const
  {
    return values.get() + this->n_elements();
  }



  template <typename Number>
  inline Number
  ReadWriteVector<Number>::operator()(const size_type global_index) const
  {
    return values[global_to_local(global_index)];
  }



  template <typename Number>
  inline Number &
  ReadWriteVector<Number>::operator()(const size_type global_index)
  {
    return values[global_to_local(global_index)];
  }



  template <typename Number>
  inline Number ReadWriteVector<Number>::
                operator[](const size_type global_index) const
  {
    return operator()(global_index);
  }



  template <typename Number>
  inline Number &ReadWriteVector<Number>::
                 operator[](const size_type global_index)
  {
    return operator()(global_index);
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::extract_subvector_to(
    const std::vector<size_type> &indices,
    std::vector<Number2> &        extracted_values) const
  {
    for (size_type i = 0; i < indices.size(); ++i)
      extracted_values[i] = operator()(indices[i]);
  }



  template <typename Number>
  template <typename ForwardIterator, typename OutputIterator>
  inline void
  ReadWriteVector<Number>::extract_subvector_to(
    ForwardIterator       indices_begin,
    const ForwardIterator indices_end,
    OutputIterator        values_begin) const
  {
    while (indices_begin != indices_end)
      {
        *values_begin = operator()(*indices_begin);
        indices_begin++;
        values_begin++;
      }
  }



  template <typename Number>
  inline Number
  ReadWriteVector<Number>::local_element(const size_type local_index) const
  {
    AssertIndexRange(local_index, this->n_elements());

    return values[local_index];
  }



  template <typename Number>
  inline Number &
  ReadWriteVector<Number>::local_element(const size_type local_index)
  {
    AssertIndexRange(local_index, this->n_elements());

    return values[local_index];
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const std::vector<size_type> &indices,
                               const std::vector<Number2> &  values)
  {
    AssertDimension(indices.size(), values.size());
    add(indices.size(), indices.data(), values.data());
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const std::vector<size_type> &  indices,
                               const ReadWriteVector<Number2> &values)
  {
    const size_type size = indices.size();
    for (size_type i = 0; i < size; ++i)
      {
        Assert(
          numbers::is_finite(values[i]),
          ExcMessage(
            "The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values[indices[i]];
      }
  }



  template <typename Number>
  template <typename Number2>
  inline void
  ReadWriteVector<Number>::add(const size_type  n_indices,
                               const size_type *indices,
                               const Number2 *  values_to_add)
  {
    for (size_type i = 0; i < n_indices; ++i)
      {
        Assert(
          numbers::is_finite(values[i]),
          ExcMessage(
            "The given value is not finite but either infinite or Not A Number (NaN)"));
        this->operator()(indices[i]) += values_to_add[i];
      }
  }



  template <typename Number>
  template <typename Functor>
  inline ReadWriteVector<Number>::FunctorTemplate<Functor>::FunctorTemplate(
    ReadWriteVector<Number> &parent,
    const Functor &          functor)
    : parent(parent)
    , functor(functor)
  {}



  template <typename Number>
  template <typename Functor>
  void
  ReadWriteVector<Number>::FunctorTemplate<Functor>::
  operator()(const size_type begin, const size_type end)
  {
    for (size_type i = begin; i < end; ++i)
      functor(parent.values[i]);
  }

#endif // ifndef DOXYGEN

} // end of namespace LinearAlgebra



/**
 * 全局函数 @p swap
 * ，它重载了C++标准库的默认实现，它使用一个临时对象。该函数简单地交换了两个向量的数据。
 * @relatesalso  Vectors
 *
 *
 */
template <typename Number>
inline void
swap(LinearAlgebra::ReadWriteVector<Number> &u,
     LinearAlgebra::ReadWriteVector<Number> &v)
{
  u.swap(v);
}


DEAL_II_NAMESPACE_CLOSE

#endif


