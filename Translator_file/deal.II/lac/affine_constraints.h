//include/deal.II-translator/lac/affine_constraints_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2021 by the deal.II authors
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

#ifndef dealii_affine_constraints_h
#define dealii_affine_constraints_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_local_storage.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_element_access.h>

#include <boost/range/iterator_range.hpp>

#include <set>
#include <type_traits>
#include <utility>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <typename>
class FullMatrix;
class SparsityPattern;
class DynamicSparsityPattern;
class BlockSparsityPattern;
class BlockDynamicSparsityPattern;
template <typename number>
class SparseMatrix;
template <typename number>
class BlockSparseMatrix;

namespace internal
{
  namespace AffineConstraints
  {
    using size_type = types::global_dof_index;

    /**
     * 这个结构包含了我们需要存储的关于每个全局条目（global_row）的所有信息：它们是由一些局部条目（local_row）或一些约束（constraint_position）直接获得的。这在用户代码中不直接使用，而是通过GlobalRowsFromLocal访问。
     * 这里执行的操作对应于将约束信息从全局自由度重塑为局部自由度（即单元格相关的DoF），同时将约束信息从压缩的行存储（每个被约束的局部dof都有一个与之相关的约束条目列表）转化为基于单元格相关DoF的压缩列存储（我们有一个全局自由度列表，对每个自由度我们有一个条目来源的局部行列表）。为了提高速度，我们另外存储一个条目是直接从本地自由度产生的，还是来自一个约束条件。
     *
     */
    struct Distributing
    {
      Distributing(const size_type global_row = numbers::invalid_size_type,
                   const size_type local_row  = numbers::invalid_size_type);

      Distributing(const Distributing &in);

      Distributing &
      operator=(const Distributing &in);

      bool
      operator<(const Distributing &in) const
      {
        return global_row < in.global_row;
      }

      size_type         global_row;
      size_type         local_row;
      mutable size_type constraint_position;
    };



    /**
     * 这个类代表了在局部层面上遇到的约束的缓存。其功能类似于
     * std::vector<std::vector<std::pair<uint,double>
     * >>，但经过调整，可以避免为每个条目频繁分配内存。数据被放入一个
     * std::vector<std::pair<uint,double>
     * >中，行的长度被固定为row_length。当这个结构被填满时，行数和行长都可以改变。在这种情况下，数据会被重新排列。这在用户代码中不直接使用，而是通过GlobalRowsFromLocal访问。
     *
     */
    template <typename number>
    struct DataCache
    {
      DataCache();

      void
      reinit();

      size_type
      insert_new_index(const std::pair<size_type, number> &pair);

      void
      append_index(const size_type                     index,
                   const std::pair<size_type, number> &pair);

      size_type
      get_size(const size_type index) const;

      const std::pair<size_type, number> *
      get_entry(const size_type index) const;

      size_type row_length;

      std::vector<std::pair<size_type, number>> data;

      std::vector<size_type> individual_size;
    };



    /**
     * 一个数据结构，收集所有来自本地贡献（单元）的全局行和它们的原点（直接/约束）。这基本上是一个由
     * "分布式
     * "结构组成的矢量，使用通过DataCache的访问。该结构提供了一些专门的排序和插入功能。
     * 在没有约束的情况下，这基本上是一个对`<uint,uint>`的列表，第一个索引是全局索引，第二个索引是本地索引。列表是按照全局索引排序的。
     * 在有约束的情况下，一个全局目标可能会得到一个贡献，因为它从一个受约束的目标获得数据。这意味着，除了直接的贡献之外，全局道夫还可能通过约束从局部道夫那里获得间接的贡献。
     * 这里执行的操作对应于将约束信息从全局自由度重塑为局部自由度（即单元格相关的自由度），同时将约束信息从压缩行存储（每个被约束的局部道夫都有一个与之相关的约束条目列表）转变为基于单元格相关自由度的压缩列存储（我们有一个全局自由度列表，对每个自由度我们有一个条目来自的局部行列表）。为了提高速度，我们另外存储一个条目是直接从本地自由度产生的，还是来自一个约束条件。
     *
     */
    template <typename number>
    class GlobalRowsFromLocal
    {
    public:
      /**
       * 构造函数。
       *
       */
      GlobalRowsFromLocal();

      void
      reinit(const size_type n_local_rows);

      void
      insert_index(const size_type global_row,
                   const size_type local_row,
                   const number    constraint_value);
      void
      sort();

      void
      print(std::ostream &os);

      /**
       * 返回结构中全局索引的数量。
       *
       */
      size_type
      size() const;

      /**
       * 返回与列表中counter_index-th条目相关的约束的数量。
       *
       */
      size_type
      size(const size_type counter_index) const;

      /**
       * 返回列表中counter_index-th条目的全局行。
       *
       */
      size_type
      global_row(const size_type counter_index) const;

      /**
       * 返回列表中counter_index-th条目的全局行。
       *
       */
      size_type &
      global_row(const size_type counter_index);

      /**
       * 返回与列表中counter_index-th条目相关的单元格矩阵中的本地行。对于受限制的行返回invalid_size_type。
       *
       */
      size_type
      local_row(const size_type counter_index) const;

      /**
       * 返回一个引用，而不是上面函数中的值。
       *
       */
      size_type &
      local_row(const size_type counter_index);

      /**
       * 返回单元格矩阵中的本地行，该行与index_in_constraint-th位置的列表中的counter_index-th条目相关。
       *
       */
      size_type
      local_row(const size_type counter_index,
                const size_type index_in_constraint) const;

      /**
       * 返回约束条件索引_in_constraint-th位置的列表中counter_index-th条目中的约束条件的值。
       *
       */
      number
      constraint_value(const size_type counter_index,
                       const size_type index_in_constraint) const;

      /**
       * 返回是否有一行有间接贡献（即至少有一个约束有非琐碎的ConstraintLine）。
       *
       */
      bool
      have_indirect_rows() const;

      /**
       * 追加一个被约束的条目。这意味着少了一个非实质性的行。
       *
       */
      void
      insert_constraint(const size_type constrained_local_dof);

      /**
       * 返回结构中受约束的道夫的数量。受约束因子并不直接对矩阵做出贡献，但是为了设置矩阵对角线和解决不均匀性，需要这些因子。
       *
       */
      size_type
      n_constraints() const;

      /**
       * 返回结构中具有不均匀性的受约束因子的数量。
       *
       */
      size_type
      n_inhomogeneities() const;

      /**
       * 这个函数告诉结构中的第1个约束是不均匀的。不均匀的约束对右手边有贡献，所以为了快速访问它们，把它们放在均匀的约束之前。
       *
       */
      void
      set_ith_constraint_inhomogeneous(const size_type i);

      /**
       * 检测到i号约束的本地行，当GlobalRowsToLocal被设置时，可以很容易地找到该行。
       *
       */
      size_type
      constraint_origin(size_type i) const;

      /**
       * 一个包含所有全局id和相应的局部id的向量，以及一个指向该数据的指针，我们在这里存储如何解决约束。
       *
       */
      std::vector<Distributing> total_row_indices;

    private:
      /**
       * 一个数据结构，保存约束条件的实际数据。
       *
       */
      DataCache<number> data_cache;

      /**
       * 一个说明有多少行的数字，约束条件不考虑。
       *
       */
      size_type n_active_rows;

      /**
       * 一个数字，表示具有不均匀约束的行的数量。
       *
       */
      size_type n_inhomogeneous_rows;
    };



    /**
     * 在调用distribution_local_to_global和add_entries_local_to_global时使用的抓取数据。为了避免频繁的内存分配，我们将这些数据保存在一个静态变量中，从一次调用到下一次。由于我们希望在矩阵中允许不同的数字类型，这是一个模板。
     * 由于每个线程都从ThreadLocalStorage中获得其私有版本的scratch数据，因此不会发生冲突的访问。为了使之有效，我们需要确保在distribution_local_to_global中没有任何调用本身可以产生任务。
     * 否则，我们可能会出现几个线程争夺数据的情况。
     * 只有通过一个访问器类才能访问从头开始的数据，这个访问器类处理访问，并标记数据为已用。
     *
     */
    template <typename number>
    struct ScratchData
    {
      /**
       * 构造函数，不做任何事情。
       *
       */
      ScratchData()
        : in_use(false)
      {}

      /**
       * 复制构造函数，不做任何事情
       *
       */
      ScratchData(const ScratchData &)
        : in_use(false)
      {}

      /**
       * 存储数据是否当前正在使用。
       *
       */
      bool in_use;

      /**
       * 用于列索引的临时数组
       *
       */
      std::vector<size_type> columns;

      /**
       * 列值的临时数组
       *
       */
      std::vector<number> values;

      /**
       * 用于块起始索引的临时数组
       *
       */
      std::vector<size_type> block_starts;

      /**
       * 用于向量索引的临时数组
       *
       */
      std::vector<size_type> vector_indices;

      /**
       * 向量值的临时数组
       *
       */
      std::vector<number> vector_values;

      /**
       * 用于重新排列行/列索引的数据数组。
       *
       */
      GlobalRowsFromLocal<number> global_rows;

      /**
       * 用于重新排序的行/列索引的数据数组。
       *
       */
      GlobalRowsFromLocal<number> global_columns;
    };
  } // namespace AffineConstraints
} // namespace internal

namespace internal
{
  namespace AffineConstraintsImplementation
  {
    template <class VectorType>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 VectorType &                                vec);

    template <class T>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 dealii::Vector<T> &                         vec);

    template <class T>
    void
    set_zero_all(const std::vector<types::global_dof_index> &cm,
                 dealii::BlockVector<T> &                    vec);
  } // namespace AffineConstraintsImplementation
} // namespace internal


template <typename number>
class AffineConstraints;
#endif

// TODO[WB]: We should have a function of the kind
//   AffineConstraints::add_constraint (const size_type constrained_dof,
//     const std::vector<std::pair<size_type, number> > &entries,
//     const number inhomogeneity = 0);
// rather than building up constraints piecemeal through add_line/add_entry
// etc. This would also eliminate the possibility of accidentally changing
// existing constraints into something pointless, see the discussion on the
// mailing list on "Tiny bug in interpolate_boundary_values" in Sept. 2010.

/**
 * 该类实现了对自由度的线性（可能是不均匀的）约束的处理。这种约束的概念和起源在 @ref
 * constraints
 * 模块中有广泛的描述。该类旨在处理相对于自由度总数而言数量有限的约束，例如百分之几到百分之三十；以及处理<i>M</i>其他自由度的线性组合，其中<i>M</i>也相对较小（最多不超过线性系统每行的平均条目数左右）。它是
 * <em> 而不是 </em> ，旨在描述全等级线性系统。
 * 在 @ref hp_paper "hp-paper "
 * 中详细描述了用于实现该类的算法。在 @ref constraints
 * 模块中也有大量的关于如何使用这个类的文档。
 *
 *  <h3>Description of constraints</h3>
 * 该类对象中的每一个 "行 "对应于一个受限自由度，行的编号为<i>i</i>，通过使用add_line()或add_lines()输入。这一行的条目是形式为(<i>j</i>, <i>a<sub>ij</sub></i>)的一对，通过add_entry()或add_entries()添加。该组织本质上是一个SparsityPattern，但只有几行包含非零元素，因此没有在其他行上浪费数据。对于通过上述机制添加的每一行，都会对形式为@f[
 * x_i = \sum_j a_{ij} x_j + b_i
 * @f]的约束自由度进行消除，其中<i>b<sub>i</sub></i>是可选的，由set_inhomogeneity()设置。因此，如果一个约束被表述为例如几个自由度的零平均值，则必须选择其中一个自由度来消除。
 * 请注意，约束条件在<i>x<sub>i</sub></i>中是线性的，而且约束条件中可能有一个常数（非均质）项。这正是我们需要的悬挂节点约束的形式，我们需要用其他自由度来约束一个自由度。还有其他可能的这种形式的条件，例如用于实现均值条件，正如在
 * step-11
 * 教程程序中所做的那样。该类的名称源于这些约束条件可以用矩阵形式表示为<b>X</b>
 * <i>x</i> =
 * <i>b</i>，然后该对象描述了矩阵<b>X</b>和向量<i>b</i>。创建/填充这种类型的对象最常用的方法是使用
 * DoFTools::make_hanging_node_constraints()
 * 函数。这些对象的使用首先在  step-6  中解释。
 * 本类型的对象被组织成行（row），但只有那些存在约束条件的行被存储。新的约束是通过使用add_line()函数添加新的行，然后使用add_entry()函数填充到给定的行中，或者使用add_entries()函数一次添加一个以上的条目。右边的元素，如果非零，可以使用set_inhomogeneity()函数进行设置。在所有约束条件被添加后，你需要调用close()，它可以压缩存储格式并对条目进行排序。
 *
 *
 * @note
 * 这个类实现的许多算法在  @ref hp_paper
 * 中讨论。这些算法也与<i>M. S. Shephard: Linear multipoint
 * constraints applied via transformation as part of a direct stiffness
 * assembly process. Int. J. Numer. Meth. Engrg., vol. 20 (1984), pp.
 * 2107-2112.</i>中所示的算法有关，不同的是，那里所示的算法完全消除了受限自由度，而我们通常将它们作为线性系统的一部分保留下来。
 *
 *
 * @ingroup dofs
 *
 * @ingroup constraints
 *
 */
template <typename number = double>
class AffineConstraints : public Subscriptor
{
public:
  /**
   * 声明容器尺寸的类型。
   *
   */
  using size_type = types::global_dof_index;

  /**
   * 一个枚举，描述了如果调用merge()函数所涉及的两个AffineConstraints对象恰好在相同的自由度上有约束，应该发生什么。
   *
   */
  enum MergeConflictBehavior
  {
    /**
     * 如果两个相关对象在同一自由度上有冲突的约束，则抛出一个异常。
     *
     */
    no_conflicts_allowed,

    /**
     * 在一个操作中 <code>cm1.merge(cm2)</code>, if <code>cm1</code> 和
     * <code>cm2</code> 对同一自由度有约束，取 <code>cm1</code>
     * 中的一个。
     *
     */
    left_object_wins,

    /**
     * 在一个操作中， <code>cm1.merge(cm2)</code>, if <code>cm1</code>
     * 和 <code>cm2</code> 对同一自由度有约束，从 <code>cm2</code>
     * 中选取一个。
     *
     */
    right_object_wins
  };

  /**
   * 构造函数。提供的IndexSet定义了在这个AffineConstraints容器内可能被约束的指数。在基于 parallel::distributed::Triangulation
   * 或 parallel::shared::Triangulation,
   * 的DoFHandler对象的计算中，应该使用本地相关的dofs集合（见
   * @ref GlossLocallyRelevantDof  ）。
   * 给定的IndexSet允许AffineConstraints容器通过不关心那些对当前处理器不重要的自由度来节省内存。另外，如果没有提供这样的IndexSet，将为<i>all</i>可能的索引创建内部数据结构，导致每个处理器的内存消耗与问题的<i>overall</i>大小成正比，而不仅仅是与当前处理器处理的整个问题的部分大小成正比。
   *
   */
  explicit AffineConstraints(const IndexSet &local_constraints = IndexSet());

  /**
   * 复制构造函数
   *
   */
  explicit AffineConstraints(const AffineConstraints &affine_constraints);

  /**
   * 移动构造函数
   *
   */
  AffineConstraints(AffineConstraints &&affine_constraints) noexcept =
    default; // NOLINT

  /**
   * 复制操作符。就像对许多其他大型对象一样，这个操作符被删除了，以避免它在一些地方被不小心使用，比如不小心将一个
   * @p AffineConstraints
   * 对象作为函数参数按值而不是按引用声明。
   * 然而，你可以使用copy_from()函数来明确地复制AffineConstraints对象。
   *
   */
  AffineConstraints &
  operator=(const AffineConstraints &) = delete;

  /**
   * 移动赋值运算符
   *
   */
  AffineConstraints &
  operator=(AffineConstraints &&affine_constraints) noexcept =
    default; // NOLINT

  /**
   * 将给定的对象复制到当前对象。
   * 这个函数的存在是因为 @p operator=() 被明确禁止。
   *
   */
  template <typename other_number>
  void
  copy_from(const AffineConstraints<other_number> &other);

  /**
   * clear()AffineConstraints对象，并提供一个带有可能被约束的线条的IndexSet。这个函数只在分布式情况下与提供不同的IndexSet有关。否则，这个例程就等同于调用clear()。详情见构造函数。
   *
   */
  void
  reinit(const IndexSet &local_constraints = IndexSet());

  /**
   * 确定我们是否可以为给定的 @p line_n.
   * 存储一个约束，这个例程只在分布式情况下有意义，并检查IndexSet是否允许存储这一行。如果不是在分布式情况下，总是返回true。
   *
   */
  bool
  can_store_line(const size_type line_n) const;

  /**
   * 返回描述本地相关行的索引集，如果有的话。请注意，如果没有给出本地行，这代表一个空的IndexSet，而否则它包含全局问题大小和本地范围。
   *
   */
  const IndexSet &
  get_local_lines() const;

  /**
   * 这个函数复制了 @p constraints_in
   * 的内容，其中的DoFs是IndexSet @p filter.
   * 中的元素，不存在于IndexSet中的元素被忽略。所有的DoFs将被转换到过滤器的本地索引空间，包括被约束的DoFs和这些条目所约束的其他DoFs。滤波器的局部索引空间是对作为滤波器元素的所有（全局）DoF的连续编号。
   * 例如，如果过滤器代表范围<tt>[10,20)</tt>，而约束对象 @p
   * constraints_in
   * 包括全局索引<tt>{7,13,14}</tt>，索引<tt>{3,4}</tt>被添加到调用的约束对象中（因为13和14是过滤器中的元素，元素13是索引中的第四元素，而14是第五元素）。
   * 这个函数提供了一个简单的方法，从一个完整的AffineConstraints中为矢量值问题中的某些矢量分量创建一个AffineConstraints，即从一个较大的AffineConstraints中提取一个对角线子块。该块是由IndexSet参数指定的。
   *
   */
  void
  add_selected_constraints(const AffineConstraints &constraints_in,
                           const IndexSet &         filter);

  /**
   * @name 添加约束  
     * @{ 
   *
   */

  /**
   * 在矩阵中添加一个新行。如果该行已经存在，那么该函数只是返回而不做任何事情。
   *
   */
  void
  add_line(const size_type line_n);

  /**
   * 为每个索引 <code>i</code> 调用第一个add_line()函数，其中
   * <code>lines[i]</code> 为真。
   * 这个函数的存在本质上是为了允许一次性添加几个形式为<i>x<sub>i</sub></i>=0的约束，其中应该添加这些约束的指数集<i>i</i>是由这个函数的参数给出。另一方面，就像重复调用单参数的add_line()函数一样，以后可以用add_entry()函数修改约束条件，以包括线性依赖，以及用set_inhomogeneity()修改不均匀性。
   *
   */
  void
  add_lines(const std::vector<bool> &lines);

  /**
   * 为参数中出现的每个索引 <code>i</code>
   * 调用第一个add_line()函数。
   * 这个函数的存在本质上是为了允许一次性添加几个形式为<i>x<sub>i</sub></i>=0的约束，其中应该添加这些约束的指数集<i>i</i>是由这个函数的参数给出。另一方面，就像重复调用单参数的add_line()函数一样，以后可以用add_entry()函数对约束条件进行修改，以包括线性依赖，以及用set_inhomogeneity()对不均匀性进行修改。
   *
   */
  void
  add_lines(const std::set<size_type> &lines);

  /**
   * 为参数中出现的每个索引 <code>i</code>
   * 调用第一个add_line()函数。
   * 这个函数的存在本质上是为了允许一次性添加几个形式为<i>x<sub>i</sub></i>=0的约束，其中应该添加这些约束的指数集<i>i</i>是由这个函数的参数给出。另一方面，就像重复调用单参数的add_line()函数一样，以后可以用add_entry()函数修改约束条件，以包括线性依赖，以及用set_inhomogeneity()修改不均匀性。
   *
   */
  void
  add_lines(const IndexSet &lines);

  /**
   * 在一个给定的行中添加一个条目。换句话说，这个函数为
   * $i$ 第1个自由度的约束条件增加了一个 $a_{ij} x_j$ 项。
   * 如果一个与这个函数调用所表示的条目具有相同的索引，那么这个函数只是简单地返回，条件是该条目的值是相同的。因此，两次输入一个约束条件并无大碍。
   * @param[in]  constrained_dof_index 被约束的自由度的索引  $i$
   * 。    @param[in]  column 被输入约束自由度 $i$
   * 的自由度的索引 $j$  。    @param[in]  权重 乘以 $x_j$
   * 的系数 $a_{ij}$  。
   *
   */
  void
  add_entry(const size_type constrained_dof_index,
            const size_type column,
            const number    weight);

  /**
   * 在一行约束条件中添加一整串条目，用成对的列索引和权重值表示。这个函数等同于多次调用前面的函数，但速度更快。
   *
   */
  void
  add_entries(
    const size_type                                  constrained_dof_index,
    const std::vector<std::pair<size_type, number>> &col_weight_pairs);

  /**
   * 为一个自由度的约束设置一个不均匀性。换句话说，它为自由度
   * $i$ 的约束添加一个常数 $b_i$
   * 。为了使其发挥作用，你需要先为给定的自由度调用add_line()。
   * @param[in]  constrained_dof_index 被约束的自由度的索引  $i$
   * 。    @param[in]  值 自由度上的约束的右手边值 $b_i$   $i$
   * 。
   *
   */
  void
  set_inhomogeneity(const size_type constrained_dof_index, const number value);

  /**
   * 关闭条目的填充。由于这种类型的矩阵的行通常是以任意顺序填充的，而且我们不想使用关联约束器来存储行，所以我们需要对行和行内的列进行排序，然后再使用矩阵。这可以通过这个函数完成。
   * 此外，零条目被丢弃，因为它们不被需要。
   * 关闭后，不再接受任何条目。如果该对象已经被关闭，那么该函数立即返回。
   * 这个函数也可以解决约束链的问题。例如，自由度13可能被约束为
   * $u_{13} = \frac{u_3}{2} + \frac{u_7}{2}$
   * ，而自由度7本身被约束为 $u_{7} = \frac{u_2}{2} +
   * \frac{u_4}{2}$  。然后，决议将是： $u_{13} = \frac{u_3}{2} +
   * \frac{u_2}{4} + \frac{u_4}{4}$
   * 。然而，请注意，这个约束图中的循环是不允许的，也就是说，例如
   * $u_4$ 本身不能直接或间接地再次约束到 $u_{13}$ 。
   *
   */
  void
  close();

  /**
   * 将作为参数的对象所代表的约束合并到这个对象所代表的约束中。两个对象都可能被关闭，也可能没有被关闭（通过之前调用它们的函数close()）。如果这个对象之前被关闭了，那么之后也会被关闭。
   * 然而，请注意，如果另一个参数是关闭的，那么合并的速度可能会大大加快。
   * 使用第二个参数的默认值，两个对象（这个对象和参数所代表的旧对象）中的约束可能并不是指相同的自由度，也就是说，一个对象中约束的自由度可能在第二个对象中没有约束。如果是这种情况，就会抛出一个异常。
   * 然而，这种行为可以通过为第二个参数提供一个不同的值来改变。
   * 默认情况下，不允许合并两个被初始化为不同IndexSet对象的AffineConstraints对象。
   * 这个行为可以通过适当地设置 @p allow_different_local_lines
   * 来改变。
   * 合并一个用IndexSet初始化的AffineConstraints和一个没有用IndexSet初始化的AffineConstraints还没有实现。
   *
   */
  void
  merge(
    const AffineConstraints &   other_constraints,
    const MergeConflictBehavior merge_conflict_behavior = no_conflicts_allowed,
    const bool                  allow_different_local_lines = false);

  /**
   * 将这个矩阵的所有条目向下移动 @p offset 行，向上移动
   * @p offset
   * 列。如果这个对象是用IndexSet初始化的，那么local_lines也会被移位。
   * 如果你正在构建块状矩阵，即所有的块都是由同一个DoFHandler对象构建的，也就是说，矩阵的大小大于自由度的数量，那么这个函数就很有用。由于几个矩阵的行和列对应着相同的自由度，你会生成几个约束对象，然后将它们移位，最后再将它们合并()在一起。
   *
   */
  void
  shift(const size_type offset);

  /**
   * 清除这个矩阵的所有条目。重置决定是否接受新条目的标志。
   * 这个函数也可以在对象为空或已经清空的情况下调用。
   *
   */
  void
  clear();

  /**
   * @}
   *
   */

  /**
   * @name 查询约束条件  
     * @{ 
   *
   */

  /**
   * 返回存储在该矩阵中的约束数量。
   *
   */
  size_type
  n_constraints() const;

  /**
   * 返回编号为 @p line_n 的自由度是否是一个约束度。
   * 请注意，如果之前调用了close()，那么这个函数会明显加快，因为此时受限自由度已经被排序，我们可以进行二进制搜索，而在调用close()之前，我们必须对所有条目进行线性搜索。
   *
   */
  bool
  is_constrained(const size_type line_n) const;

  /**
   * 返回该自由度是否被约束，以及是否被约束到只有一个权重为1的其他自由度。因此，该函数返回该自由度是否会被简单地消除，而恰恰是其他一个自由度。
   * 如果自由度完全不受约束，或者它受制于一个以上的其他自由度，或者它只受制于一个自由度，但其权重不同于1，则该函数返回
   * @p false 。
   *
   */
  bool
  is_identity_constrained(const size_type line_n) const;

  /**
   * 返回给定的两个自由度是否被一个平等约束所连接，该约束要么约束index1以便
   * <code>index1=index2</code> ，要么约束index2以便
   * <code>index2=index1</code>  。
   *
   */
  bool
  are_identity_constrained(const size_type line_n_1,
                           const size_type line_n_2) const;

  /**
   * 返回一个道夫所约束的其他道夫的最大数量。
   * 例如，在2d中，一个悬挂的节点只受制于它的两个邻居，所以返回值是2。然而，对于高阶元素和/或高维度，或其他类型的约束，这个数字就不明显了。
   * 这个名字表明，在系统矩阵内，对受约束节点的引用是间接指向它所受约束的节点的。
   *
   */
  size_type
  max_constraint_indirections() const;

  /**
   * 在道夫被约束的情况下，返回<tt>true</tt>，并且有一个非琐碎的不均匀值设置给道夫。
   *
   */
  bool
  is_inhomogeneously_constrained(const size_type index) const;

  /**
   * 如果AffineConstraints中的所有约束都是同质的，则返回<tt>false</tt>；如果至少有一个不均匀性，则返回<tt>true</tt>。
   *
   */
  bool
  has_inhomogeneities() const;

  /**
   * 如果一条线被约束，则返回一个指向条目向量的指针，如果道夫没有被约束，则返回一个零指针。
   *
   */
  const std::vector<std::pair<size_type, number>> *
  get_constraint_entries(const size_type line_n) const;

  /**
   * 返回存储在受约束道夫 @p
   * line_n中的不均匀性的值。不受约束的道夫也会返回一个零值。
   *
   */
  number
  get_inhomogeneity(const size_type line_n) const;

  /**
   * 将当前对象所代表的约束打印到给定的流中。    对于每个形式为@f[
   * x_{42} = 0.5 x_2 + 0.25 x_{14} + 2.75
   * @f]的约束，这个函数将写出一连串的行，看起来像这样。
   * @code
   * 42 2 : 0.5
   * 42 14 : 0.25
   * 42 : 2.75
   * @endcode
   * 最后一行只有在不均匀性（这里：2.75）非零时才会显示。
   * 像上面这样的线块对每个受限自由度都是重复的。
   *
   */
  void
  print(std::ostream &out) const;

  /**
   * 用'dot'格式写出约束的图形。dot'是一个程序，它可以接受一个节点列表，并产生一个约束自由度和它们所约束的自由度的图形表示。
   * 这个函数的输出可以作为'dot'程序的输入，该程序可以将图形转换成postscript、png、xfig和其他一些格式的图形表示。
   * 这个函数的存在主要是为了调试的目的。
   *
   */
  void
  write_dot(std::ostream &) const;

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 添加与给定向量中的指数相关的约束指数。
   * 在调用此函数后，索引向量包含初始元素和所有相关的约束指数。这个函数对这些元素进行排序，并抑制重复的元素。
   *
   */
  void
  resolve_indices(std::vector<types::global_dof_index> &indices) const;

  /**
   * @}
   *
   */

  /**
   * @name  在线性系统创建后消除其约束 
     * @{ 
   *
   */

  /**
   * 浓缩一个稀疏的模式。这个函数的名字模仿了我们用来凝结线性系统的函数的名字，但对于目前的语境来说，它有点名不副实。这是因为在线性系统的背景下，我们消除了线性系统的某些行和列，即我们
   * "减少 "或 "浓缩
   * "了线性系统。另一方面，在当前情况下，函数并不从稀疏模式中删除非零条目。相反，它将那些非零条目的位置添加到稀疏性模式中，这些非零条目以后将被用于浓缩线性系统的受限自由度的过程。
   * 由于该函数向稀疏性模式添加了新的非零条目，给定的稀疏性模式必须不被压缩。当前对象必须是封闭的。稀疏度模式在函数的最后被压缩。
   *
   */
  void
  condense(SparsityPattern &sparsity) const;

  /**
   * 与上面的函数相同，但浓缩了方形块的稀疏模式。
   *
   */
  void
  condense(BlockSparsityPattern &sparsity) const;

  /**
   * 与上面的函数相同，但浓缩了方形压缩的稀疏模式。
   *
   */
  void
  condense(DynamicSparsityPattern &sparsity) const;

  /**
   * 与上面的功能相同，但浓缩了方形压缩的稀疏模式。
   *
   */
  void
  condense(BlockDynamicSparsityPattern &sparsity) const;

  /**
   * 压缩一个给定的矩阵，即消除矩阵中对应于受限自由度的行和列。
   * 更详细的信息请参见该类的一般文档。
   *
   */
  void
  condense(SparseMatrix<number> &matrix) const;

  /**
   * 与上述函数相同，但浓缩了方形块稀疏矩阵。
   *
   */
  void
  condense(BlockSparseMatrix<number> &matrix) const;

  /**
   * 对给定的向量进行原地压缩。 @p VectorType
   * 可以是Vector<float>, Vector<number>, BlockVector<tt><...></tt>,
   * PETSc或Trilinos向量包装类，或任何其他具有相同接口的类型。注意，这个函数不考虑任何不均匀性，如果有任何不均匀性，会抛出一个异常。
   * 在这种情况下，请同时使用矩阵和矢量的函数。
   * @note
   * 这个函数对MPI向量不起作用。请使用带有两个向量参数的condense()来代替。
   *
   */
  template <class VectorType>
  void
  condense(VectorType &vec) const;

  /**
   * 该函数将 @p vec_ghosted 中的值复制并浓缩到 @p
   * 中输出。在串行代码中，它等同于调用condense
   * (vec)。如果以并行方式调用， @p vec_ghosted
   * 应该包含幽灵元素，而 @p output 不应该。
   *
   */
  template <class VectorType>
  void
  condense(const VectorType &vec_ghosted, VectorType &output) const;

  /**
   * 通过消除对应于受限自由度的线性系统的行和列，浓缩一个给定的矩阵和一个给定的向量。与矩阵相关的稀疏性模式需要被浓缩和压缩。
   * 这个函数是应用不均匀约束的合适选择。
   * 当前对象必须被关闭才能调用这个函数。
   * 更详细的信息请参见该类的一般文档。
   *
   */
  template <class VectorType>
  void
  condense(SparseMatrix<number> &matrix, VectorType &vector) const;

  /**
   * 与上述函数相同，但浓缩了方形块稀疏矩阵和向量。
   *
   */
  template <class BlockVectorType>
  void
  condense(BlockSparseMatrix<number> &matrix, BlockVectorType &vector) const;

  /**
   * 将一个向量中所有受约束的DoF的值设置为零。  @p
   * VectorType可以是Vector<float>, Vector<number>,
   * BlockVector<tt><...></tt>,
   * PETSc或Trilinos矢量包装类，或任何其他具有相同接口的类型。
   *
   */
  template <class VectorType>
  void
  set_zero(VectorType &vec) const;

  /**
   * @}
   *
   */

  /**
   * @name  在线性系统创建过程中消除其约束  
     * @{ 
   *
   */

  /**
   * 这个函数接收一个局部贡献的向量(  @p local_vector)
   * 对应于 @p
   * local_dof_indices中给出的自由度指数，并将它们分配到全局向量中。换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在大多数情况下，这些局部贡献将是在一个单元或一个单元的面上进行整合的结果。然而，只要
   * @p local_vector 和 @p
   * 的local_dof_indices有相同的元素数，这个函数就很乐意接受它所给的任何东西。
   * 与DoFAccessor类中的类似函数相比，这个函数也照顾到了约束，即如果
   * @p local_dof_indices
   * 中的一个元素属于一个受约束的节点，那么与其将 @p
   * local_vector 中的相应元素写入 @p
   * global_vector，不如将该元素分配到这个特定自由度受约束的global
   * vector中的条目。
   * 因此，通过使用这个函数将局部贡献分配给全局对象，就省去了在向量和矩阵完全集合后对缩合函数的调用。另一方面，由于这个原因，这个函数不仅写进了
   * @p local_dof_indices
   * 数组所列举的条目，而且还写进了（可能）其他必要的条目。
   * 请注意，这个函数将应用所有的约束，就像它们是同质的一样。为了正确设置非均质约束，请使用带有矩阵参数的类似函数，或者同时带有矩阵和矢量参数的函数。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局向量允许同时访问，并且在同一时间不访问具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者的现场进行确认。这个方法内部没有锁定机制来防止数据竞赛。
   * @param[in]  local_vector 本地贡献的向量。    @param[in]
   * local_dof_indices 本地贡献向量对应的本地自由度指数。
   * @param[out]  global_vector
   * 全局向量，所有局部贡献将被添加到该向量中。
   *
   */
  template <class InVector, class OutVector>
  void
  distribute_local_to_global(const InVector &              local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             OutVector &                   global_vector) const;

  /**
   * 这个函数接收一个与 @p
   * local_dof_indices中给出的自由度指数相对应的局部贡献向量(
   * @p local_vector)
   * ，并将其分配到全局向量中。换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在大多数情况下，这些局部贡献将是在一个单元或一个单元的面上进行整合的结果。然而，只要
   * @p local_vector 和 @p
   * 的local_dof_indices有相同的元素数，这个函数就很乐意接受它所给的任何东西。
   * 与DoFAccessor类中的类似函数相比，这个函数也照顾到了约束，即如果
   * @p local_dof_indices
   * 中的一个元素属于一个受约束的节点，那么与其将 @p
   * local_vector 中的相应元素写入 @p
   * global_vector，不如将该元素分配到这个特定自由度受约束的global
   * vector中的条目。
   * 因此，通过使用这个函数将局部贡献分配给全局对象，就省去了在向量和矩阵完全集合后对缩合函数的调用。另一方面，由于这个原因，该函数不仅写入
   * @p local_dof_indices
   * 数组所列举的条目，而且还（可能）写入其他必要的条目。这包括写进矩阵的对角线元素，如果相应的自由度受到限制。
   * 第四个参数<tt>local_matrix</tt>是为了在人们想只对向量应用不均匀约束时使用。这种情况可能是人们想在一个有不均匀约束的问题上装配一个右手边的向量，但全局矩阵已经在之前装配好了。一个典型的例子是一个时间步进算法，其中刚度矩阵被组装一次，而右手边在每个时间步进中被更新。然而，请注意，局部矩阵的列中的条目必须与已经写入全局矩阵的条目完全相同。否则，这个函数将不能正确处理不均匀性。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局向量允许同时访问，并且不同时访问具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者的现场进行确认。这个方法内部没有锁定机制来防止数据竞赛。
   *
   */
  template <typename VectorType>
  void
  distribute_local_to_global(const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             VectorType &                  global_vector,
                             const FullMatrix<number> &    local_matrix) const;

  /**
   * 与前一个函数相同，只是它使用了两个（可能）不同的索引集来正确处理局部矩阵由两个相邻元素组合计算时的不均匀性，例如DG中的边缘积分项。注意，在这两个元素具有不同多项式程度的情况下，局部矩阵是矩形的。
   * <tt>local_dof_indices_row</tt>是局部矩阵的行指数集合，<tt>local_dof_indices_col</tt>是局部矩阵的列指数集合。<tt>diagonal=false</tt>表示这两个索引集是否相等。
   * 如果两个索引集相等，<tt>diagonal</tt>必须设置为真，否则我们直接使用前面的函数。如果两个索引集不同（对角线=false），<tt>global_vector</tt>会被修改以处理不均匀性，但不会增加<tt>local_vector</tt>的条目。请注意，对于DG的内边缘积分不对右手边贡献任何数值。
   *
   */
  template <typename VectorType>
  void
  distribute_local_to_global(
    const Vector<number> &        local_vector,
    const std::vector<size_type> &local_dof_indices_row,
    const std::vector<size_type> &local_dof_indices_col,
    VectorType &                  global_vector,
    const FullMatrix<number> &    local_matrix,
    bool                          diagonal = false) const;

  /**
   * 在结果向量中输入一个单一的值，服从约束。
   *
   */
  template <class VectorType>
  void
  distribute_local_to_global(const size_type index,
                             const number    value,
                             VectorType &    global_vector) const;

  /**
   * 这个函数接收一个指向局部贡献向量的指针（ @p
   * local_vector），对应于 @p
   * local_dof_indices中给出的自由度指数，并将它们分配到全局向量中。换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在大多数情况下，这些局部贡献将是在一个单元或一个单元的面上进行整合的结果。然而，只要
   * @p
   * local_dof_indices中的条目表示合理的全局向量条目，这个函数就会对它所给出的任何东西感到满意。
   * 如果 @p local_dof_indices
   * 中的一个元素属于受限节点，那么与其将 @p
   * local_vector中的相应元素写入 @p global_vector,
   * ，不如将该元素分配到该特定自由度受限的全局向量中的条目。
   * 因此，通过使用这个函数将局部贡献分配给全局对象，可以省去在向量和矩阵完全集合后对缩合函数的调用。注意，这个函数完全忽略了不均匀约束。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局向量允许同时访问，并且不同时访问具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者的现场进行确认。这个方法里面没有锁定机制来防止数据竞赛。
   *
   */
  template <typename ForwardIteratorVec,
            typename ForwardIteratorInd,
            class VectorType>
  void
  distribute_local_to_global(ForwardIteratorVec local_vector_begin,
                             ForwardIteratorVec local_vector_end,
                             ForwardIteratorInd local_indices_begin,
                             VectorType &       global_vector) const;

  /**
   * 这个函数接收一个本地贡献的矩阵（  @p local_matrix)
   * 对应于 @p
   * local_dof_indices中给出的自由度指数，并将它们分配到全局矩阵中。换句话说，这个函数实现了一个[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 在大多数情况下，这些局部贡献将是在一个单元或一个单元的面上进行整合的结果。然而，只要
   * @p local_matrix 和 @p
   * 的local_dof_indices有相同的元素数，这个函数就很乐意接受它所给的任何东西。
   * 与DoFAccessor类中的类似函数相比，这个函数也照顾到了约束，即如果
   * @p local_dof_indices
   * 中的一个元素属于一个受约束的节点，那么与其将 @p
   * local_matrix 中的相应元素写入 @p
   * global_matrix，不如将该元素分配给这个特定自由度受约束的全局矩阵中的条目。
   * 通过这个方案，我们永远不会写进受约束自由度的行或列。为了确保得到的矩阵仍然可以被倒置，我们需要对对应于受约束节点的对角线元素做一些处理。因此，如果
   * @p
   * local_dof_indices中的一个自由度受到约束，我们就在矩阵中分配相应的条目，但同时也将局部矩阵的对角线条目的绝对值加到全局矩阵的相应条目上。
   * 假设离散算子是正定的，这就保证了对角线条目总是非零的、正的，并且与矩阵的其他条目具有相同的数量级。另一方面，在解决源问题时
   * $Au=f$
   * ，对角线元素的精确值并不重要，因为无论如何，各自自由度的值将被后面的distribution()调用覆盖。
   * @note
   * 上述程序在矩阵的频谱上增加了不可预见的人工特征值的数量。因此，在这种情况下，建议使用带有两个局部索引向量的等效函数。
   * 通过使用这个函数来分配对全局对象的局部贡献，可以省去在向量和矩阵完全集合后对缩合函数的调用。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局矩阵允许同时访问，并且访问的不是具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者的现场进行确认。这个方法内部没有锁定机制来防止数据竞赛。
   *
   */
  template <typename MatrixType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * 这个函数的作用与上面的函数几乎相同，但可以处理一般的矩形矩阵。实现这一点的主要区别是，受限行中的对角线条目不被触及，而被填充为任意值。
   * 由于对应于被消除的自由度的对角线项没有被设置，如果应用于正方形矩阵，其结果可能有一个零特征值。这在解决所产生的问题时必须加以考虑。对于解决源问题
   * $Au=f$
   * ，可以在建立矩阵后通过以下形式的代码来设置对角线条目
   * @code
   * for (unsigned int i=0;i<matrix.m();++i)
   *   if (constraints.is_constrained(i))
   *     matrix.diag_element(i) = 1.;
   * @endcode
   * 这里使用的1的值是任意的，但在Krylov空间方法的背景下是不关键的，因为它对应于一个不变的子空间。如果其他矩阵项的大小与机器精度相近，最好对其进行调整。
   * 对于解决特征值问题，这只会增加一个虚假的零特征值（其倍率可能大于1）。
   * 考虑到这一点，其他的就不用改了。
   *
   */
  template <typename MatrixType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &row_indices,
                             const std::vector<size_type> &col_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * 这个函数与上面的函数对一般矩形矩阵的作用几乎相同，但在行和列索引上使用不同的AffineConstraints对象。惯例是，行指数根据调用的AffineConstraints
   * <code>*this</code>
   * 进行约束，而列指数则根据给定的AffineConstraints
   * <code>column_affine_constraints</code>
   * 进行约束。这个函数允许处理这样的情况：矩阵的行和列由不同的函数空间表示，并有各自的指数列举，例如，在混合有限元问题中，有独立的DoFHandler对象，或者在多网格方法中不同层次之间的通量矩阵。
   * 与其他为行和列指数提供单独槽的方法一样，这种方法不为消除的自由度增加对角线条目。更详细的描述见那里。
   *
   */
  template <typename MatrixType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const std::vector<size_type> &row_indices,
                             const AffineConstraints &column_affine_constraints,
                             const std::vector<size_type> &column_indices,
                             MatrixType &                  global_matrix) const;

  /**
   * 这个函数根据调用AffineConstraints指定的约束条件，同时将元素写入矩阵和向量中。
   * 换句话说，它同时对矩阵和向量执行相应函数的[散点操作](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing))。
   * 这个函数也能正确处理不均匀约束。关于参数use_inhomogeneities_for_rhs见
   * @ref constraints 模块中的文档。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局矩阵和向量允许同时访问，并且在同一时间不访问具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者的现场进行确认。这个方法内部没有锁定机制来防止数据竞赛。
   *
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             bool use_inhomogeneities_for_rhs = false) const;

  /**
   * 做一个与distribut_local_to_global()函数类似的操作，该函数将写条目分配到受限自由度的矩阵中，只是这里我们不写进矩阵，而只是分配稀疏模式条目。    正如 @ref hp_paper "hp-paper "
   * 和 step-27
   * 中解释的那样，首先分配一个稀疏模式，然后再回来为那些由于消除受限自由度而被写入的矩阵条目分配额外的条目（使用
   * AffineConstraints::condense()
   * ），可能是一个非常昂贵的过程。
   * 更便宜的做法是立即分配这些条目，而不需要对稀疏模式对象进行第二次传递。这个函数正是这样做的。
   * 因为该函数只分配稀疏模式中的条目，所以它需要知道的是相互耦合的自由度。
   * 与前一个函数不同的是，没有写入实际值，所以第二个输入参数在这里是不必要的。
   * 该函数的第三个参数keep_constrained_entries决定了该函数是否要在稀疏性模式中分配条目，这些条目以后会在矩阵浓缩时被设置为零。如果矩阵是在无约束的情况下建立的，那么这些条目是必要的，只是后来被浓缩。如果使用该类的distribut_local_to_global()函数构建矩阵，则不需要这些条目，该函数在将本地矩阵复制到全局对象时，会立即分配条目。这个参数的默认值是true，意味着要分配少数条目，这些条目以后可能被设置为0。
   * 默认情况下，该函数将第一个参数中给出的所有索引对的条目添加到稀疏模式中（除非keep_constrained_entries为false）。然而，有时人们希望只添加所有这些对的一个子集。在这种情况下，可以使用最后一个参数，该参数指定了一个布尔掩码，哪些索引对应该被考虑。如果某对指数的掩码是假的，那么这对指数的稀疏性模式将不被添加，无论其中一个或两个指数是否对应于受限自由度。
   * 这个函数通常不是从用户代码中调用的，而是在传递AffineConstraints对象时在
   * DoFTools::make_sparsity_pattern() 函数中使用。
   * @note
   * 这个函数本身是线程安全的，也就是说，当几个线程同时调用它时，它也能正常工作。然而，只有当底层的全局稀疏性模式允许同时访问，并且在同一时间不访问具有相同全局索引的行时，该函数调用才是线程安全的。这需要从调用者那里得到确认。这个方法内部没有锁定机制来防止数据竞赛。
   *
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global(
    const std::vector<size_type> &local_dof_indices,
    SparsityPatternType &         sparsity_pattern,
    const bool                    keep_constrained_entries = true,
    const Table<2, bool> &        dof_mask = Table<2, bool>()) const;

  /**
   * 与另一个函数类似，但用于非二次疏散模式。
   *
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global(
    const std::vector<size_type> &row_indices,
    const std::vector<size_type> &col_indices,
    SparsityPatternType &         sparsity_pattern,
    const bool                    keep_constrained_entries = true,
    const Table<2, bool> &        dof_mask = Table<2, bool>()) const;

  /**
   * 这个函数从全局向量中导入数值（  @p global_vector)
   * 通过将约束条件应用于本地数值的向量，以迭代器格式表示。
   * 在大多数情况下，本地值将由单元格上的本地dof值来识别。然而，只要
   * @p
   * local_dof_indices中的条目表示合理的全局向量条目，这个函数就会对它所给的任何东西感到满意。
   * 如果 @p local_dof_indices
   * 中的一个元素属于受约束的节点，那么与其将 @p
   * global_vector中的相应元素写入 @p local_vector,
   * ，不如像各自的分发函数那样解决约束，也就是说，本地条目是由这个特定自由度受约束的全局条目构造的。
   * 与DoFAccessor类中的类似函数get_dof_values相比，这个函数不需要正确设置约束值（即调用distribution）。
   *
   */
  template <typename ForwardIteratorVec,
            typename ForwardIteratorInd,
            class VectorType>
  void
  get_dof_values(const VectorType & global_vector,
                 ForwardIteratorInd local_indices_begin,
                 ForwardIteratorVec local_vector_begin,
                 ForwardIteratorVec local_vector_end) const;

  /**
   * @}
   *
   */

  /**
   * @name  在求解线性系统后处理约束条件  
     * @{ 
   *
   */

  /**
   * 给定一个向量，将所有受约束的自由度设置为数值，使约束得到满足。例如，如果当前对象存储了约束条件
   * $x_3=\frac 12 x_1 + \frac 12 x_2$
   * ，那么这个函数将从给定的向量中读取 $x_1$ 和 $x_2$
   * 的值并根据这个约束条件设置元素 $x_3$
   * 。同样地，如果当前对象存储了约束条件 $x_{42}=208$
   * ，那么这个函数将把给定向量的第42个元素设置为208。
   * @note 如果这个函数是用一个平行向量 @p vec,
   * 调用的，那么这个向量不能包含鬼魂元素。
   *
   */
  template <class VectorType>
  void
  distribute(VectorType &vec) const;

  /**
   * @}
   *
   */

  /**
   * 该类代表AffineConstraints对象中的一个约束。
   *
   */
  struct ConstraintLine
  {
    /**
     * 一个数据类型，我们在其中存储构成一个约束的同质部分的条目列表。
     *
     */
    using Entries = std::vector<std::pair<size_type, number>>;

    /**
     * 这一行的全局DoF索引。由于只存储了很少的行，我们不能假设一个特定的顺序，必须明确地存储索引。
     *
     */
    size_type index;

    /**
     * 这一行中的条目的行号和值。
     * 对于我们使用向量而不是地图的原因及其后果，与对
     * AffineConstraints::lines. 所说的相同。
     *
     */
    Entries entries;

    /**
     * 不均匀性的值。
     *
     */
    number inhomogeneity;

    /**
     * 默认构造函数。
     *
     */
    ConstraintLine(const size_type &index         = numbers::invalid_dof_index,
                   const Entries &  entries       = {},
                   const number &   inhomogeneity = 0.0);

    /**
     * 复制构造器。
     *
     */
    template <typename ConstraintLineType>
    ConstraintLine(const ConstraintLineType &other);

    /**
     * 拷贝赋值。
     *
     */
    template <typename ConstraintLineType>
    ConstraintLine &
    operator=(const ConstraintLineType &other);

    /**
     * 这个运算符有点奇怪，而且不直观：它比较了两行的行号。我们需要这个来对行进行排序；事实上，我们可以用一个比较谓词来做这个。
     * 然而，这种方式更容易，尽管不直观，因为两行确实没有神赐的顺序关系。
     *
     */
    bool
    operator<(const ConstraintLine &) const;

    /**
     * 这个运算符同样也很奇怪：它检查两个操作数的行指数是否相等，而不考虑行的内容可能不同的事实。
     *
     */
    bool
    operator==(const ConstraintLine &) const;

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中写入和读取此对象的数据，以便进行序列化。
     *
     */
    template <class Archive>
    void
    serialize(Archive &ar, const unsigned int)
    {
      ar &index &entries &inhomogeneity;
    }
  };

  /**
   * 在LineRange容器中使用的迭代器类型的别名。
   *
   */
  using const_iterator = typename std::vector<ConstraintLine>::const_iterator;

  /**
   * get_lines()使用的返回类型的别名。
   *
   */
  using LineRange = boost::iterator_range<const_iterator>;

  /**
   * 返回一个包含（const）对存储在AffineConstraints容器中的所有线条的迭代器的范围对象。这样的范围对于初始化基于范围的for循环是很有用的，正如C++11所支持的那样。
   * @return  一个用于半开范围<code>[this->begin(),
   * this->end())</code>的行条目的范围对象。
   *
   */
  const LineRange
  get_lines() const;

  /**
   * 检查当前对象在分布式计算的所有处理器上是否一致。
   * 该方法检查所有处理器是否同意由 @p locally_active_dofs.
   * 给出的其本地行的约束，该方法是一个集体操作，只有当所有处理器都一致时才会返回
   * @p true 。    请提供由 Utilities::MPI::all_gather(MPI_Comm,
   * DoFHandler::locally_owned_dofs())
   * 返回的每个处理器拥有的DoF，作为 @p locally_owned_dofs
   * ，以及 DoFTools::extract_locally_active_dofs() 的结果，作为 @p
   * locally_active_dofs. 。
   * 前者用于确定特定DoF的所有权，而后者作为需要检查的行的集合。
   * 如果 @p verbose 被设置为 @p true,
   * ，额外的调试信息将被写入 std::cout. 。
   * @note
   * 这个方法交换了本地活动行的所有约束信息，因此对于大型计算来说很慢，可能只应该在调试模式下使用。我们不检查get_local_lines()返回的所有行，而只检查本地活动的行，因为我们允许处理器不知道一些本地相关行。
   * @return
   * 所有AffineConstraints对象是否一致。在所有处理器上返回相同的值。
   *
   */
  bool
  is_consistent_in_parallel(const std::vector<IndexSet> &locally_owned_dofs,
                            const IndexSet &             locally_active_dofs,
                            const MPI_Comm &             mpi_communicator,
                            const bool                   verbose = false) const;

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcMatrixIsClosed);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException0(ExcMatrixNotClosed);
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcLineInexistant,
                 size_type,
                 << "The specified line " << arg1 << " does not exist.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException4(ExcEntryAlreadyExists,
                 size_type,
                 size_type,
                 number,
                 number,
                 << "The entry for the indices " << arg1 << " and " << arg2
                 << " already exists, but the values " << arg3 << " (old) and "
                 << arg4 << " (new) differ "
                 << "by " << (arg4 - arg3) << ".");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcDoFConstrainedToConstrainedDoF,
                 int,
                 int,
                 << "You tried to constrain DoF " << arg1 << " to DoF " << arg2
                 << ", but that one is also constrained. This is not allowed!");
  /**
   * 异常情况。
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcDoFIsConstrainedFromBothObjects,
                 size_type,
                 << "Degree of freedom " << arg1
                 << " is constrained from both object in a merge operation.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcDoFIsConstrainedToConstrainedDoF,
                 size_type,
                 << "In the given argument a degree of freedom is constrained "
                 << "to another DoF with number " << arg1
                 << ", which however is constrained by this object. This is not"
                 << " allowed.");
  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException1(ExcRowNotStoredHere,
                 size_type,
                 << "The index set given to this constraints object indicates "
                 << "constraints for degree of freedom " << arg1
                 << " should not be stored by this object, but a constraint "
                 << "is being added.");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcColumnNotStoredHere,
                 size_type,
                 size_type,
                 << "The index set given to this constraints object indicates "
                 << "constraints using degree of freedom " << arg2
                 << " should not be stored by this object, but a constraint "
                 << "for degree of freedom " << arg1 << " uses it.");

  /**
   * 异常情况
   * @ingroup Exceptions
   *
   */
  DeclException2(ExcIncorrectConstraint,
                 int,
                 int,
                 << "While distributing the constraint for DoF " << arg1
                 << ", it turns out that one of the processors "
                 << "who own the " << arg2 << " degrees of freedom that x_"
                 << arg1 << " is constrained against does not know about "
                 << "the constraint on x_" << arg1
                 << ". Did you not initialize the AffineConstraints container "
                 << "with the appropriate locally_relevant set so "
                 << "that every processor who owns a DoF that constrains "
                 << "another DoF also knows about this constraint?");

  template <typename>
  friend class AffineConstraints;

private:
  /**
   * 存储矩阵的行数。
   * 条目通常是以任意顺序追加的，插入到向量中最好在最后进行，所以在所有条目插入后，顺序是不指定的。对条目的排序是在调用<tt>close()</tt>函数时进行的。
   * 我们可以不使用矢量，而是使用一个关联数组，像地图一样来存储这些行。然而，这将意味着一个更加碎片化的堆，因为它要分配许多小对象，而且还会使这个矩阵的使用变得更加缓慢。
   *
   */
  std::vector<ConstraintLine> lines;

  /**
   * 一个size_type的列表，包含受限自由度的ConstraintLine的位置，如果自由度不受限，则为
   * numbers::invalid_size_type 。 numbers::invalid_size_type
   * 的返回值因此返回是否有一个给定自由度指数的约束线。请注意，这个类对自由度的实际数量没有概念，所以如果我们检查某个自由度是否有约束线，那么这个向量实际上可能比我们检查的自由度的索引短。
   * 这个字段的存在是因为当添加一个新的约束线时，我们必须弄清楚它是否已经存在。以前，我们只是简单地在未排序的约束行列表中行走，直到我们碰到终点或找到它。如果N是约束条件的数量，这个算法是O(N)，这使得它在插入所有约束条件时是O(N^2)。对于有许多约束条件的大问题来说，这很容易占用总运行时间的5-10%。有了这个领域，我们可以节省这个时间，因为我们在O(1)时间内找到任何约束条件，或者得到它某个自由度不受约束。
   * 更糟糕的是，遍历现有的约束列表需要从内存的许多不同地方读取。因此，在大型3D应用中，add_line()函数在整个计算时间中表现得非常突出，主要是因为它产生了大量的高速缓存缺失。
   * 这也应该通过使用O(1)算法来访问这个数组的字段来解决。
   * 这个字段在其他一些情况下也很有用，例如，当人们需要随机访问约束条件时，就像所有在飞行中应用约束条件的函数一样，同时将单元贡献添加到向量和矩阵中。
   *
   */
  std::vector<size_type> lines_cache;

  /**
   * 这个IndexSet是用来限制保存在AffineConstraints中的行数的子集。这是必要的，因为在分布式计算中，lines_cache向量会变得太大。
   *
   */
  IndexSet local_lines;

  /**
   * 存储数组是否被排序。 如果是，就不能添加新的条目。
   *
   */
  bool sorted;

  mutable Threads::ThreadLocalStorage<
    internal::AffineConstraints::ScratchData<number>>
    scratch_data;

  /**
   * 内部函数，使用local_lines计算向量lines_cache中的行 @p
   * line_n 的索引。
   *
   */
  size_type
  calculate_line_index(const size_type line_n) const;

  /**
   * 这个函数实际上实现了标准（非块）矩阵的local_to_global函数。
   *
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             const bool use_inhomogeneities_for_rhs,
                             const std::integral_constant<bool, false>) const;

  /**
   * 这个函数实际上是为块状矩阵实现了local_to_global函数。
   *
   */
  template <typename MatrixType, typename VectorType>
  void
  distribute_local_to_global(const FullMatrix<number> &    local_matrix,
                             const Vector<number> &        local_vector,
                             const std::vector<size_type> &local_dof_indices,
                             MatrixType &                  global_matrix,
                             VectorType &                  global_vector,
                             const bool use_inhomogeneities_for_rhs,
                             const std::integral_constant<bool, true>) const;

  /**
   * 这个函数实际上是为标准（非块）稀疏类型实现了local_to_global函数。
   *
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global(const std::vector<size_type> &local_dof_indices,
                              SparsityPatternType &         sparsity_pattern,
                              const bool            keep_constrained_entries,
                              const Table<2, bool> &dof_mask,
                              const std::integral_constant<bool, false>) const;

  /**
   * 这个函数实际上是为块状稀疏度类型实现local_to_global函数。
   *
   */
  template <typename SparsityPatternType>
  void
  add_entries_local_to_global(const std::vector<size_type> &local_dof_indices,
                              SparsityPatternType &         sparsity_pattern,
                              const bool            keep_constrained_entries,
                              const Table<2, bool> &dof_mask,
                              const std::integral_constant<bool, true>) const;

  /**
   * 为distribution_local_to_global函数提供内部辅助函数。
   * 创建一个受影响的全局行列表用于分发，包括条目来自的本地行。该列表根据全局行的索引进行排序。
   *
   */
  void
  make_sorted_row_list(const std::vector<size_type> &local_dof_indices,
                       internal::AffineConstraints::GlobalRowsFromLocal<number>
                         &global_rows) const;

  /**
   * add_entries_local_to_global函数的内部辅助函数。
   * 创建一个受影响的行的列表，用于分发，没有任何额外的信息，否则与其他make_sorted_row_list()函数类似。
   *
   */
  void
  make_sorted_row_list(const std::vector<size_type> &local_dof_indices,
                       std::vector<size_type> &      active_dofs) const;

  /**
   * distribute_local_to_global函数的内部辅助函数。
   *
   */
  template <typename MatrixScalar, typename VectorScalar>
  typename ProductType<VectorScalar, MatrixScalar>::type
  resolve_vector_entry(
    const size_type                                                 i,
    const internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows,
    const Vector<VectorScalar> &    local_vector,
    const std::vector<size_type> &  local_dof_indices,
    const FullMatrix<MatrixScalar> &local_matrix) const;
};

 /* ---------------- template and inline functions ----------------- */ 

template <typename number>
inline AffineConstraints<number>::AffineConstraints(
  const IndexSet &local_constraints)
  : lines()
  , local_lines(local_constraints)
  , sorted(false)
{
  // make sure the IndexSet is compressed. Otherwise this can lead to crashes
  // that are hard to find (only happen in release mode).
  // see tests/mpi/affine_constraints_crash_01
  local_lines.compress();
}

template <typename number>
inline AffineConstraints<number>::AffineConstraints(
  const AffineConstraints &affine_constraints)
  : Subscriptor()
  , lines(affine_constraints.lines)
  , lines_cache(affine_constraints.lines_cache)
  , local_lines(affine_constraints.local_lines)
  , sorted(affine_constraints.sorted)
{}

template <typename number>
inline void
AffineConstraints<number>::add_line(const size_type line_n)
{
  Assert(sorted == false, ExcMatrixIsClosed());

  // the following can happen when we compute with distributed meshes and dof
  // handlers and we constrain a degree of freedom whose number we don't have
  // locally. if we don't abort here the program will try to allocate several
  // terabytes of memory to resize the various arrays below :-)
  Assert(line_n != numbers::invalid_size_type, ExcInternalError());
  const size_type line_index = calculate_line_index(line_n);

  // check whether line already exists; it may, in which case we can just quit
  if (is_constrained(line_n))
    return;

  // if necessary enlarge vector of existing entries for cache
  if (line_index >= lines_cache.size())
    lines_cache.resize(std::max(2 * static_cast<size_type>(lines_cache.size()),
                                line_index + 1),
                       numbers::invalid_size_type);

  // push a new line to the end of the list
  lines.emplace_back();
  lines.back().index         = line_n;
  lines.back().inhomogeneity = 0.;
  lines_cache[line_index]    = lines.size() - 1;
}



template <typename number>
inline void
AffineConstraints<number>::add_entry(const size_type constrained_dof_index,
                                     const size_type column,
                                     const number    weight)
{
  Assert(sorted == false, ExcMatrixIsClosed());
  Assert(constrained_dof_index != column,
         ExcMessage("Can't constrain a degree of freedom to itself"));

  // Ensure that the current line is present in the cache:
  const size_type line_index = calculate_line_index(constrained_dof_index);
  Assert(line_index < lines_cache.size(),
         ExcMessage("The current AffineConstraints does not contain the line "
                    "for the current entry. Call AffineConstraints::add_line "
                    "before calling this function."));

  // if in debug mode, check whether an entry for this column already exists
  // and if it's the same as the one entered at present
  //
  // in any case: exit the function if an entry for this column already
  // exists, since we don't want to enter it twice
  Assert(lines_cache[line_index] != numbers::invalid_size_type,
         ExcInternalError());
  Assert(!local_lines.size() || local_lines.is_element(column),
         ExcColumnNotStoredHere(constrained_dof_index, column));
  ConstraintLine *line_ptr = &lines[lines_cache[line_index]];
  Assert(line_ptr->index == constrained_dof_index, ExcInternalError());
  for (const auto &p : line_ptr->entries)
    if (p.first == column)
      {
        Assert(std::abs(p.second - weight) < 1.e-14,
               ExcEntryAlreadyExists(
                 constrained_dof_index, column, p.second, weight));
        return;
      }

  line_ptr->entries.emplace_back(column, weight);
}



template <typename number>
inline void
AffineConstraints<number>::set_inhomogeneity(
  const size_type constrained_dof_index,
  const number    value)
{
  const size_type line_index = calculate_line_index(constrained_dof_index);
  Assert(line_index < lines_cache.size() &&
           lines_cache[line_index] != numbers::invalid_size_type,
         ExcMessage("call add_line() before calling set_inhomogeneity()"));
  Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
  ConstraintLine *line_ptr = &lines[lines_cache[line_index]];
  line_ptr->inhomogeneity  = value;
}



template <typename number>
template <class VectorType>
inline void
AffineConstraints<number>::set_zero(VectorType &vec) const
{
  // since lines is a private member, we cannot pass it to the functions
  // above. therefore, copy the content which is cheap
  std::vector<size_type> constrained_lines(lines.size());
  for (unsigned int i = 0; i < lines.size(); ++i)
    constrained_lines[i] = lines[i].index;
  internal::AffineConstraintsImplementation::set_zero_all(constrained_lines,
                                                          vec);
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::n_constraints() const
{
  return lines.size();
}

template <typename number>
inline bool
AffineConstraints<number>::is_constrained(const size_type index) const
{
  const size_type line_index = calculate_line_index(index);
  return ((line_index < lines_cache.size()) &&
          (lines_cache[line_index] != numbers::invalid_size_type));
}

template <typename number>
inline bool
AffineConstraints<number>::is_inhomogeneously_constrained(
  const size_type line_n) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return false;
  else
    {
      Assert(lines_cache[line_index] < lines.size(), ExcInternalError());
      return !(lines[lines_cache[line_index]].inhomogeneity == number(0.));
    }
}

template <typename number>
inline const std::vector<std::pair<types::global_dof_index, number>> *
AffineConstraints<number>::get_constraint_entries(const size_type line_n) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return nullptr;
  else
    return &lines[lines_cache[line_index]].entries;
}

template <typename number>
inline number
AffineConstraints<number>::get_inhomogeneity(const size_type line_n) const
{
  // check whether the entry is constrained. could use is_constrained, but
  // that means computing the line index twice
  const size_type line_index = calculate_line_index(line_n);
  if (line_index >= lines_cache.size() ||
      lines_cache[line_index] == numbers::invalid_size_type)
    return 0;
  else
    return lines[lines_cache[line_index]].inhomogeneity;
}

template <typename number>
inline types::global_dof_index
AffineConstraints<number>::calculate_line_index(const size_type line_n) const
{
  // IndexSet is unused (serial case)
  if (!local_lines.size())
    return line_n;

  Assert(local_lines.is_element(line_n), ExcRowNotStoredHere(line_n));

  return local_lines.index_within_set(line_n);
}

template <typename number>
inline bool
AffineConstraints<number>::can_store_line(size_type line_n) const
{
  return !local_lines.size() || local_lines.is_element(line_n);
}

template <typename number>
inline const IndexSet &
AffineConstraints<number>::get_local_lines() const
{
  return local_lines;
}

template <typename number>
template <class VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const size_type index,
  const number    value,
  VectorType &    global_vector) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());

  if (is_constrained(index) == false)
    global_vector(index) += value;
  else
    {
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(index)]];
      for (size_type j = 0; j < position.entries.size(); ++j)
        global_vector(position.entries[j].first) +=
          value * position.entries[j].second;
    }
}

template <typename number>
template <typename ForwardIteratorVec,
          typename ForwardIteratorInd,
          class VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  ForwardIteratorVec local_vector_begin,
  ForwardIteratorVec local_vector_end,
  ForwardIteratorInd local_indices_begin,
  VectorType &       global_vector) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());
  for (; local_vector_begin != local_vector_end;
       ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        internal::ElementAccess<VectorType>::add(*local_vector_begin,
                                                 *local_indices_begin,
                                                 global_vector);
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          for (size_type j = 0; j < position.entries.size(); ++j)
            internal::ElementAccess<VectorType>::add(
              (*local_vector_begin) * position.entries[j].second,
              position.entries[j].first,
              global_vector);
        }
    }
}

template <typename number>
template <class InVector, class OutVector>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const InVector &              local_vector,
  const std::vector<size_type> &local_dof_indices,
  OutVector &                   global_vector) const
{
  Assert(local_vector.size() == local_dof_indices.size(),
         ExcDimensionMismatch(local_vector.size(), local_dof_indices.size()));
  distribute_local_to_global(local_vector.begin(),
                             local_vector.end(),
                             local_dof_indices.begin(),
                             global_vector);
}

template <typename number>
template <typename ForwardIteratorVec,
          typename ForwardIteratorInd,
          class VectorType>
inline void
AffineConstraints<number>::get_dof_values(
  const VectorType & global_vector,
  ForwardIteratorInd local_indices_begin,
  ForwardIteratorVec local_vector_begin,
  ForwardIteratorVec local_vector_end) const
{
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());
  for (; local_vector_begin != local_vector_end;
       ++local_vector_begin, ++local_indices_begin)
    {
      if (is_constrained(*local_indices_begin) == false)
        *local_vector_begin = global_vector(*local_indices_begin);
      else
        {
          const ConstraintLine &position =
            lines[lines_cache[calculate_line_index(*local_indices_begin)]];
          typename VectorType::value_type value = position.inhomogeneity;
          for (size_type j = 0; j < position.entries.size(); ++j)
            value += (global_vector(position.entries[j].first) *
                      position.entries[j].second);
          *local_vector_begin = value;
        }
    }
}

template <typename MatrixType>
class BlockMatrixBase;
template <typename SparsityPatternType>
class BlockSparsityPatternBase;
template <typename number>
class BlockSparseMatrixEZ;

namespace internal
{
  namespace AffineConstraints
  {
    /**
     * 一个 "特质
     * "类，可以用来确定一个给定的类型是否是块状矩阵类型。比如说。
     * @code
     * IsBlockMatrix<SparseMatrix<number> >::value
     * @endcode
     * 的值是`false`，而
     * @code
     * IsBlockMatrix<BlockSparseMatrix<number> >::value
     * @endcode
     * 是真的。这在模板上下文中有时很有用，我们可能想根据模板类型是表示普通的还是块状的矩阵类型来做不同的事情。          @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    template <typename MatrixType>
    struct IsBlockMatrix
    {
    private:
      /**
       * 如果该类派生于BlockMatrixBase，则重载返回true，这就是块状矩阵的作用（BlockSparseMatrixEZ除外）。
       *
       */
      template <typename T>
      static std::true_type
      check(const BlockMatrixBase<T> *);

      /**
       * BlockSparseMatrixEZ的重载，在编写这个类时，它是唯一不从BlockMatrixBase派生的块状矩阵。
       *
       */
      template <typename T>
      static std::true_type
      check(const BlockSparseMatrixEZ<T> *);

      /**
       * 为所有其他潜在的、当时显然不是块状矩阵的类型捕捉所有的东西。
       *
       */
      static std::false_type
      check(...);

    public:
      /**
       * 一个可静态计算的值，表示该类的模板参数是否是块状矩阵（事实上，该类型是否从BlockMatrixBase<T>派生，或者是其他块状矩阵类型之一）。
       *
       */
      static const bool value =
        std::is_same<decltype(check(std::declval<MatrixType *>())),
                     std::true_type>::value;
    };

    // instantiation of the static member
    template <typename MatrixType>
    const bool IsBlockMatrix<MatrixType>::value;


    /**
     * 一个可以用来确定给定类型是否是块状稀疏模式类型的类。在此，它与IsBlockMatrix类相匹配。          @see   @ref GlossBlockLA  "块（线性代数）"
     *
     */
    template <typename MatrixType>
    struct IsBlockSparsityPattern
    {
    private:
      /**
       * 如果该类派生自BlockSparsityPatternBase，则重载返回true，这就是块状稀疏模式的作用。
       *
       */
      template <typename T>
      static std::true_type
      check(const BlockSparsityPatternBase<T> *);

      /**
       * 捕捉所有其他潜在的类型，然后显然不是块状稀疏模式。
       *
       */
      static std::false_type
      check(...);

    public:
      /**
       * 一个可静态计算的值，表明该类的模板参数是否是块状稀疏模式（事实上，该类型是否来自BlockSparsityPatternBase<T>）。
       *
       */
      static const bool value =
        std::is_same<decltype(check(std::declval<MatrixType *>())),
                     std::true_type>::value;
    };

    // instantiation of the static member
    template <typename MatrixType>
    const bool IsBlockSparsityPattern<MatrixType>::value;

  } // namespace AffineConstraints
} // namespace internal



template <typename number>
template <typename other_number>
inline void
AffineConstraints<number>::copy_from(
  const AffineConstraints<other_number> &other)
{
  lines.clear();
  lines.insert(lines.begin(), other.lines.begin(), other.lines.end());
  lines_cache = other.lines_cache;
  local_lines = other.local_lines;
  sorted      = other.sorted;
}



template <typename number>
template <typename MatrixType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number> &    local_matrix,
  const std::vector<size_type> &local_dof_indices,
  MatrixType &                  global_matrix) const
{
  // create a dummy and hand on to the function actually implementing this
  // feature in the cm.templates.h file.
  Vector<typename MatrixType::value_type> dummy(0);
  distribute_local_to_global(
    local_matrix,
    dummy,
    local_dof_indices,
    global_matrix,
    dummy,
    false,
    std::integral_constant<
      bool,
      internal::AffineConstraints::IsBlockMatrix<MatrixType>::value>());
}



template <typename number>
template <typename MatrixType, typename VectorType>
inline void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number> &    local_matrix,
  const Vector<number> &        local_vector,
  const std::vector<size_type> &local_dof_indices,
  MatrixType &                  global_matrix,
  VectorType &                  global_vector,
  bool                          use_inhomogeneities_for_rhs) const
{
  // enter the internal function with the respective block information set,
  // the actual implementation follows in the cm.templates.h file.
  distribute_local_to_global(
    local_matrix,
    local_vector,
    local_dof_indices,
    global_matrix,
    global_vector,
    use_inhomogeneities_for_rhs,
    std::integral_constant<
      bool,
      internal::AffineConstraints::IsBlockMatrix<MatrixType>::value>());
}



template <typename number>
template <typename SparsityPatternType>
inline void
AffineConstraints<number>::add_entries_local_to_global(
  const std::vector<size_type> &local_dof_indices,
  SparsityPatternType &         sparsity_pattern,
  const bool                    keep_constrained_entries,
  const Table<2, bool> &        dof_mask) const
{
  // enter the internal function with the respective block information set,
  // the actual implementation follows in the cm.templates.h file.
  add_entries_local_to_global(
    local_dof_indices,
    sparsity_pattern,
    keep_constrained_entries,
    dof_mask,
    std::integral_constant<bool,
                           internal::AffineConstraints::IsBlockSparsityPattern<
                             SparsityPatternType>::value>());
}



template <typename number>
inline AffineConstraints<number>::ConstraintLine::ConstraintLine(
  const size_type &                                                  index,
  const typename AffineConstraints<number>::ConstraintLine::Entries &entries,
  const number &inhomogeneity)
  : index(index)
  , entries(entries)
  , inhomogeneity(inhomogeneity)
{}



template <typename number>
template <typename ConstraintLineType>
inline AffineConstraints<number>::ConstraintLine::ConstraintLine(
  const ConstraintLineType &other)
{
  this->index = other.index;

  entries.clear();
  entries.insert(entries.begin(), other.entries.begin(), other.entries.end());

  this->inhomogeneity = other.inhomogeneity;
}



template <typename number>
template <typename ConstraintLineType>
inline typename AffineConstraints<number>::ConstraintLine &
AffineConstraints<number>::ConstraintLine::
operator=(const ConstraintLineType &other)
{
  this->index = other.index;

  entries.clear();
  entries.insert(entries.begin(), other.entries.begin(), other.entries.end());

  this->inhomogeneity = other.inhomogeneity;

  return *this;
}

DEAL_II_NAMESPACE_CLOSE

#endif


