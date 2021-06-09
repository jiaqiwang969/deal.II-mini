//include/deal.II-translator/grid/filtered_iterator_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2020 by the deal.II authors
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

#ifndef dealii_filtered_iterator_h
#  define dealii_filtered_iterator_h


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/iterator_range.h>

#  include <deal.II/grid/tria_iterator_base.h>

#  include <memory>
#  include <set>
#  include <tuple>

DEAL_II_NAMESPACE_OPEN


/**
 * 在这个命名空间中，声明了一些可以作为FilteredIterator类中的过滤器的类。这些过滤器或者检查二进制信息（例如，
 * IteratorFilters::Active
 * 过滤器类检查指向的对象是否是活动的），或者通过与规定的值比较来检查有价值的信息（例如，LevelEqualTo过滤器类检查所考虑的迭代器指向的对象的级别是否等于构造时给过滤器的一个值。
 * 关于这些类的使用例子以及对过滤器的要求，请参见FilteredIterator类的一般描述。
 *
 *
 * @ingroup Iterators
 *
 *
 */
namespace IteratorFilters
{
  /**
   * 如果迭代器指向一个活动对象或迭代器过了终点，则评估为真的过滤器。
   * @ingroup Iterators
   *
   */
  class Active
  {
  public:
    /**
     * 评估迭代器，如果该对象是活动的或已过终点，则返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };

  /**
   * 如果迭代器指向一个设置了用户标志的对象或者迭代器超过了终点，则评估为真的过滤器。参见 @ref
   * GlossUserFlags 了解关于用户标志的信息。
   * @ingroup Iterators
   *
   */
  class UserFlagSet
  {
  public:
    /**
     * 评估迭代器，如果对象有设定的用户标志或过了终点，则返回true。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * 如果迭代器指向一个没有设置用户标志的对象或者迭代器超过了终点，则评估为真的过滤器。与前一个类别相反的过滤器。
   * @ingroup Iterators
   *
   */
  class UserFlagNotSet
  {
  public:
    /**
     * 评估迭代器，如果对象有一个未设置的用户标志或超过终点，则返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * 迭代器的过滤器，如果迭代器已经过了终点，或者指向的对象的级别等于给构造函数的一个值，则评估为真。
   * @ingroup Iterators
   *
   */
  class LevelEqualTo
  {
  public:
    /**
     * 构造函数。存储迭代器应被评估为真的水平。
     *
     */
    LevelEqualTo(const unsigned int level);

    /**
     * 评估操作符。如果指向的对象的级别等于存储的值，或者迭代器已过终点，则返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * 要与级别进行比较的存储值。
     *
     */
    const unsigned int level;
  };



  /**
   * 迭代器的过滤器，如果迭代器已经过了终点，或者所指向的对象的子域ID等于给构造函数的一个值，假设迭代器允许查询子域ID，则评价为真。
   * @ingroup Iterators
   *
   */
  class SubdomainEqualTo
  {
  public:
    /**
     * 构造函数。存储子域，迭代器应被评估为真。
     *
     */
    SubdomainEqualTo(const types::subdomain_id subdomain_id);

    /**
     * 评估操作符。如果指向的对象的子域等于存储的值或迭代器已过终点，则返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * 用于比较子域的存储值。
     *
     */
    const types::subdomain_id subdomain_id;
  };



  /**
   * 迭代器的过滤器，如果一个单元格为当前处理器所拥有，即如果它是一个 @ref GlossLocallyOwnedCell "本地拥有的单元格"
   * ，则评估为真。    该类在  step-32  中使用，与  @ref
   * distributed  模块的方法相关。
   * @ingroup Iterators
   *
   */
  class LocallyOwnedCell
  {
  public:
    /**
     * 评估运算符。如果单元格是本地拥有的，则返回true。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };



  /**
   * 迭代器的过滤器，如果单元格的水平子域id等于当前处理器id，则评估为真。
   * @ingroup Iterators
   *
   */
  class LocallyOwnedLevelCell
  {
  public:
    /**
     * 评价运算符。如果单元格的水平子域id等于当前处理器id，则返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * 迭代器的过滤器，如果指向的对象的迭代器等于给构造器的一个值或一组值，假设迭代器允许查询材料id，则评价为真。
   * @ingroup Iterators
   *
   */
  class MaterialIdEqualTo
  {
  public:
    /**
     * 构造器。存储迭代器应被评估为真的材料ID，并说明迭代器是否必须为本地所有。
     *
     */
    MaterialIdEqualTo(const types::material_id material_id,
                      const bool               only_locally_owned = false);

    /**
     * 构造函数。存储一个材料ID的集合，其迭代器应被评估为真，并说明迭代器是否必须是本地拥有的。
     *
     */
    MaterialIdEqualTo(const std::set<types::material_id> &material_ids,
                      const bool only_locally_owned = false);

    /**
     * 评价操作符。如果指向的对象的材料ID在存储的允许值集合内相等，则返回真，如果需要，如果单元格是本地拥有的。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * 用来比较材料ID的存储值。
     *
     */
    const std::set<types::material_id> material_ids;
    /**
     * 标志，说明是否只有本地拥有的单元格必须返回真。
     *
     */
    const bool only_locally_owned;
  };

  /**
   * 迭代器的过滤器，如果指向的对象的迭代器等于给构造器的一个值或一组值，假设迭代器允许查询活动的FE索引，则评价为真。
   * @ingroup Iterators
   *
   */
  class ActiveFEIndexEqualTo
  {
  public:
    /**
     * 构造函数。存储活动的FE索引，哪些迭代器应被评估为真，并说明迭代器是否必须为本地所有。
     *
     */
    ActiveFEIndexEqualTo(const unsigned int active_fe_index,
                         const bool         only_locally_owned = false);

    /**
     * 构造函数。存储一个活动FE指数的集合，迭代器必须被评估为真，如果迭代器必须是本地拥有的，则说明。
     *
     */
    ActiveFEIndexEqualTo(const std::set<unsigned int> &active_fe_indices,
                         const bool only_locally_owned = false);

    /**
     * 评估运算符。如果所指向的对象的活动FE索引在存储的允许值集合内相等，则返回真，如果需要的话，如果单元格是本地拥有的。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * 用于比较材料ID的存储值。
     *
     */
    const std::set<unsigned int> active_fe_indices;
    /**
     * 标志，说明是否只有本地拥有的单元格必须返回真。
     *
     */
    const bool only_locally_owned;
  };

  /**
   * 迭代器的过滤器，如果指向的对象的迭代器在边界上，则评估为真。
   * @ingroup Iterators
   *
   */
  class AtBoundary
  {
  public:
    /**
     * 评估迭代器，如果对象在边界上，返回真。
     *
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };
} // namespace IteratorFilters


/**
 * 这个类通过只迭代满足给定过滤器（称为 <em> 谓词 </em>
 * ，遵循C++标准库的符号）的元素，提供了对一系列三角形或DoFHandler迭代器的某种看法。一旦用谓词和迭代器的值进行初始化，如果调用操作符++或\--，过滤的迭代器会跳到满足谓词的下一个或上一个元素。位于两者之间但不满足谓词的中间迭代器值被跳过。因此，在某一类对象上写循环是非常简单的，不需要明确写出它们在每个循环迭代中必须满足的条件。如果函数被调用时有一对迭代器，表示它们将在一个范围内行动，通过选择一个过滤的迭代器而不是通常的迭代器，这尤其有帮助。
 * 该类在  step-32  中使用。
 *
 *  <h3>Predicates</h3>
 * 代表迭代器必须满足的条件的对象只需要提供一个允许调用评估操作符的接口，即
 * <code>bool operator() (const BaseIterator&)</code>
 * 。这包括函数指针以及实现<code>bool operator ()(const
 * BaseIterator&)</code>的类。然后，FilteredIterator将跳过所有该函数的返回值为
 * <code>false</code>  的对象。
 *
 *  一个简单的有效谓词的例子如下：给定函数
 *
 * @code
 * template <typename BIterator>
 * bool level_equal_to_3 (const BIterator& c)
 * {
 *   return (static_cast<unsigned int>(c->level()) == 3);
 * };
 * @endcode
 * 那么
 *
 * @code
 * &level_equal_to_3<typename Triangulation<dim>::active_cell_iterator>
 * @endcode
 * 是一个有效的谓词。 同样地，给定以下二元函数
 * @code
 * template <typename BIterator>
 * bool level_equal_to (const BIterator&     c,
 *                      const unsigned int level)
 * {
 *   return (static_cast<unsigned int>(c->level()) == level);
 * };
 * @endcode
 * 那么
 *
 * @code
 * [](const BIterator& c){ return level_equal_to<active_cell_iterator>(c, 3);}
 * @endcode
 * 是另一个有效的谓词（这里：一个函数，如果迭代器过了终点或者水平等于第二个参数，则返回真；这个第二个参数在创建lambda函数时被视为固定的）。
 * 最后，类可以是谓词。下面这个类就是一个。
 *
 * @code
 * class Active
 * {
 * public:
 *   template <class Iterator>
 *   bool operator () (const Iterator &i) const
 *   {
 *     return i->is_active();
 *   }
 * };
 * @endcode
 * 而这种类型的对象可以作为谓词使用。同样地，这个更复杂的也可以使用。
 *
 * @code
 * class SubdomainEqualTo
 * {
 * public:
 *   SubdomainEqualTo (const types::subdomain_id subdomain_id)
 *     : subdomain_id (subdomain_id)
 *   {};
 *
 *   template <class Iterator>
 *   bool operator () (const Iterator &i) const
 *   {
 *     return (i->subdomain_id() == subdomain_id);
 *   }
 *
 * private:
 *   const types::subdomain_id subdomain_id;
 * };
 * @endcode
 * 像 <code>SubdomainEqualTo(3)</code>
 * 这样的对象就可以作为谓词使用。
 * 因为每当一个谓词被评估时，都会检查被检查的迭代器是否真的有效（即没有超过终点），所以在谓词内部不需要对这种情况进行检查。
 * 很多过滤器类已经在IteratorFilters命名空间中实现了，但是按照上面的例子，编写不同的过滤器也很简单。
 *
 *  <h3>Initialization of filtered iterators</h3>
 * 过滤的迭代器在构造时被赋予一个谓词，这个谓词不能再被改变。如果这个谓词是作为一个模板参数给类的话，这种行为是可以预期的，但由于这将使过滤迭代器的声明成为一场噩梦，我们宁愿把谓词作为一个不可改变的实体给构造器。请注意，人们可以将一个具有一个谓词的过滤迭代器分配给另一个具有另一种类型的过滤迭代器；然而，这并不
 * <em>
 * 改变分配给迭代器的谓词，只有指示迭代器的指针被改变。
 * 如果一个被过滤的迭代器没有被分配一个底层（未被过滤的）迭代器类型的值，那么将采取默认值。然而，如果给构造函数一个值，该值必须超过终点，或者必须满足谓词。例如，如果谓词只在对象的级别等于3时才评估为真，那么
 * <code>tria.begin_active(3)</code> 将是一个有效的选择，而
 * <code>tria.begin()</code>
 * 则不是，因为后者也会返回非活动单元的迭代器，这些单元总是从级别0开始。
 * 由于人们经常只有一些迭代器，并希望将过滤后的迭代器设置为第一个满足谓词的迭代器（例如，第一个设置了用户标志的迭代器，或者第一个具有给定子域id的迭代器），因此有一些赋值函数#set_to_next_positive和#set_to_previous_positive，它们将满足谓词的下一个或上一个迭代器赋值，即。
 * 也就是说，它们沿着迭代器列表的任一方向追踪，直到找到一个匹配的迭代器（或过去的迭代器）。像
 * <code>operator=</code>
 * 一样，它们返回过滤后的迭代器的结果值。
 *
 *  <h3>Examples</h3>
 * 下面的调用计算了有设置用户标志的活动单元的数量。
 *
 * @code
 * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 *   begin (IteratorFilters::UserFlagSet()),
 *   end (IteratorFilters::UserFlagSet());
 * begin.set_to_next_positive(tria.begin_active());
 * end = tria.end();
 * n_flagged_cells = std::distance (begin, end);
 * @endcode
 * 请注意，通过 @p set_to_next_positive
 * 的调用，第一个有设置用户标志的单元被分配到 @p begin
 * 迭代器。对于结束迭代器来说，没有必要进行这样的调用，因为过去结束迭代器总是满足所有的谓词。
 * 同样的情况可以通过下面的片段来实现，虽然比较难读。
 *
 * @code
 * using FI =
 *   FilteredIterator<typename Triangulation<dim>::active_cell_iterator>;
 * n_flagged_cells =
 *   std::distance (
 *     FI(IteratorFilters::UserFlagSet()).set_to_next_positive(
 *       tria.begin_active()),
 *     FI(IteratorFilters::UserFlagSet(), tria.end()));
 * @endcode
 * 它依赖于这样一个事实：如果我们用一个给定的谓词创建一个未命名的过滤迭代器，但没有迭代器的值，并给它分配关于这个谓词的下一个正值，它就会返回自己，然后作为
 * @p std::distance
 * 函数的第一个参数。这个过程对于这里的这个函数的结束元素来说是没有必要的，因为过去的结束迭代器总是满足谓词，所以我们可以在构造函数中直接将这个值分配给过滤的迭代器。
 * 最后，下面的循环只在子域id等于3的单元格上组装矩阵。
 *
 * @code
 * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 * cell (IteratorFilters::SubdomainEqualTo(3)),
 * endc (IteratorFilters::SubdomainEqualTo(3), tria.end());
 * cell.set_to_next_positive (tria.begin_active());
 * for (; cell!=endc; ++cell)
 * assemble_local_matrix (cell);
 * @endcode
 *
 * 由于定义了过滤和未过滤的迭代器之间的比较，我们也可以让最后一个例子中的
 * @p endc 变量为 Triangulation::active_cell_iterator
 * 类型，因为它是不变的，其值不取决于过滤器。
 *
 *
 * @ingroup grid
 *
 * @ingroup Iterators
 *
 */
template <typename BaseIterator>
class FilteredIterator : public BaseIterator
{
public:
  /**
   * 对底层迭代器的访问器类型的类型定义。
   *
   */
  using AccessorType = typename BaseIterator::AccessorType;

  /**
   * 构造函数。将迭代器设置为默认状态并使用给定的谓词来过滤后续的赋值和迭代。
   *
   */
  template <typename Predicate>
  FilteredIterator(Predicate p);

  /**
   * 构造函数。使用给定的谓词进行过滤，用给定的值初始化迭代器。
   * 如果初始值 @p bi 不满足谓词 @p p
   * ，那么它就会被推进，直到我们碰到过去结束的迭代器，或者谓词被满足。例如，这允许写这样的代码
   * @code
   * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
   * cell (IteratorFilters::SubdomainEqualTo(13),
   *       triangulation.begin_active());
   * @endcode
   * 如果单元格 <code>triangulation.begin_active()</code>
   * 没有等于13的subdomain_id，那么迭代器将自动推进到第一个有的单元格。
   *
   */
  template <typename Predicate>
  FilteredIterator(Predicate p, const BaseIterator &bi);

  /**
   * 复制构造函数。复制给定参数的谓词和迭代器的值。
   *
   */
  FilteredIterator(const FilteredIterator &fi);

  /**
   * 赋值运算符。复制参数的迭代器值，但是正如在类的文档中所讨论的，参数的谓词并没有被复制。参数的迭代器值必须满足被赋值对象的谓词，这是在其构造时给出的。
   *
   */
  FilteredIterator &
  operator=(const FilteredIterator &fi);

  /**
   * 赋值运算符。复制参数的迭代器值，并保留该对象的谓词。给定的迭代器值必须满足分配到的对象的谓词，如在其构造时给出的。
   *
   */
  FilteredIterator &
  operator=(const BaseIterator &fi);

  /**
   * 从 @p bi
   * 开始搜索下一个满足此对象的谓词的迭代器，并将其分配给此对象。
   * 因为过滤后的迭代器会自动转换为底层的基本迭代器类型，所以你也可以给一个过滤后的迭代器作为这个函数的参数。
   *
   */
  FilteredIterator &
  set_to_next_positive(const BaseIterator &bi);

  /**
   * 如上所述，但是从 @p bi
   * 向后搜索满足此对象谓词的前一个迭代器，并将其分配给此对象。
   * 由于过滤后的迭代器会自动转换为底层的基本迭代器类型，你也可以给一个过滤后的迭代器作为这个函数的参数。
   *
   */
  FilteredIterator &
  set_to_previous_positive(const BaseIterator &bi);

  /**
   * 比较this和给定对象的基础迭代器值是否相等。
   * 我们不对谓词的平等性进行比较。
   *
   */
  bool
  operator==(const FilteredIterator &fi) const;

  /**
   * 比较这个对象和给定对象的基本迭代器值是否相等。
   * 此对象的谓词与此操作无关。
   *
   */
  bool
  operator==(const BaseIterator &fi) const;

  /**
   * 比较此对象和给定对象的基础迭代器值的不平等。
   * 我们不对谓词的平等性进行比较。
   *
   */
  bool
  operator!=(const FilteredIterator &fi) const;

  /**
   * 比较此对象的底层迭代器值与给定对象的不平等性。
   * 此对象的谓词与此操作无关。
   *
   */
  bool
  operator!=(const BaseIterator &fi) const;

  /**
   * 比较此对象和给定对象的基础迭代器值的排序。
   * 我们不比较谓词。
   *
   */
  bool
  operator<(const FilteredIterator &fi) const;

  /**
   * 比较这个对象和给定对象的底层迭代器值的排序。
   * 此对象的谓词与此操作无关。
   *
   */
  bool
  operator<(const BaseIterator &fi) const;

  /**
   * 前缀推进操作：移动到满足前提条件的下一个迭代器值，并返回新的迭代器值。
   *
   */
  FilteredIterator &
  operator++();

  /**
   * 后缀推进运算符：移动到满足谓词的下一个迭代器值，并返回旧的迭代器值。
   *
   */
  FilteredIterator
  operator++(int);

  /**
   * 前缀递减运算符：移动到满足谓词的前一个迭代器值，并返回新的迭代器值。
   *
   */
  FilteredIterator &
  operator--();

  /**
   * 后缀进位运算符：移动到满足谓词的前一个迭代器值，并返回旧的迭代器值。
   *
   */
  FilteredIterator
  operator--(int);

  /**
   * 异常情况。
   *
   */
  DeclException1(
    ExcInvalidElement,
    BaseIterator,
    << "The element " << arg1
    << " with which you want to compare or which you want to"
    << " assign from is invalid since it does not satisfy the predicate.");

private:
  /**
   * 用于封装谓词对象的基类。由于谓词可以有不同的类型，而且我们不想把这些类型编码到过滤迭代器类的模板参数列表中，所以我们使用一个带有抽象函数的基类和模板化的派生类，通过虚拟函数实现对实际谓词类型的使用。
   * @ingroup Iterators
   *
   */
  class PredicateBase
  {
  public:
    /**
     * 将析构器标记为虚拟，以允许通过指向基类的指针进行破坏。
     *
     */
    virtual ~PredicateBase() = default;

    /**
     * 抽象函数，在派生类中表示对给定迭代器的谓词的评估。
     *
     */
    virtual bool
    operator()(const BaseIterator &bi) const = 0;

    /**
     * 生成这个对象的一个副本，即这个指针的实际类型。
     *
     */
    virtual std::unique_ptr<PredicateBase>
    clone() const = 0;
  };


  /**
   * 上述抽象基类的实际实现。使用一个模板参数来表示谓词的实际类型，并存储它的副本。当虚拟函数被调用时，用存储的谓词副本评估给定的迭代器。
   * @ingroup Iterators
   *
   */
  template <typename Predicate>
  class PredicateTemplate : public PredicateBase
  {
  public:
    /**
     * 构造器。取一个谓词并存储它的一个副本。
     *
     */
    PredicateTemplate(const Predicate &predicate);

    /**
     * 用存储的谓词的副本来评估迭代器。
     *
     */
    virtual bool
    operator()(const BaseIterator &bi) const override;

    /**
     * 生成这个对象的一个副本，即这个指针的实际类型。
     *
     */
    virtual std::unique_ptr<PredicateBase>
    clone() const override;

  private:
    /**
     * 谓词的拷贝。
     *
     */
    const Predicate predicate;
  };

  /**
   * 指向一个对象的指针，该对象封装了给与构造函数的谓词的实际数据类型。
   *
   */
  std::unique_ptr<const PredicateBase> predicate;
};



/**
 * 给出基础迭代器和谓词，创建一个FilteredIterator类型的对象。
 * 这个函数使创建临时对象（例如作为函数参数）变得简单得多，因为人们不需要明确地用手指定基迭代器的类型。
 *
 * - 它是在这里自动推导出来的。
 * @relatesalso FilteredIterator
 *
 *
 */
template <typename BaseIterator, typename Predicate>
FilteredIterator<BaseIterator>
make_filtered_iterator(const BaseIterator &i, const Predicate &p)
{
  FilteredIterator<BaseIterator> fi(p);
  fi.set_to_next_positive(i);
  return fi;
}



namespace internal
{
  namespace FilteredIteratorImplementation
  {
    // The following classes create a nested sequence of
    // FilteredIterator<FilteredIterator<...<BaseIterator>...>> with as many
    // levels of FilteredIterator classes as there are elements in the TypeList
    // if the latter is given as a std::tuple<Args...>.
    template <typename BaseIterator, typename TypeList>
    struct NestFilteredIterators;

    template <typename BaseIterator, typename Predicate>
    struct NestFilteredIterators<BaseIterator, std::tuple<Predicate>>
    {
      using type = ::dealii::FilteredIterator<BaseIterator>;
    };

    template <typename BaseIterator, typename Predicate, typename... Targs>
    struct NestFilteredIterators<BaseIterator, std::tuple<Predicate, Targs...>>
    {
      using type = ::dealii::FilteredIterator<
        typename NestFilteredIterators<BaseIterator,
                                       std::tuple<Targs...>>::type>;
    };
  } // namespace FilteredIteratorImplementation
} // namespace internal



/**
 * 使用一个Predicate对给定的迭代器范围进行过滤。这允许替换。
 *
 * @code
 * DoFHandler<dim> dof_handler;
 * ...
 * for (const auto &cell : dof_handler.active_cell_iterators())
 *   {
 *     if (cell->is_locally_owned())
 *       {
 *         fe_values.reinit (cell);
 *         ...do the local integration on 'cell'...;
 *       }
 *   }
 * @endcode
 * 由:
 *
 * @code
 * DoFHandler<dim> dof_handler;
 * ...
 * const auto filtered_iterators_range =
 *   filter_iterators(dof_handler.active_cell_iterators(),
 *                    IteratorFilters::LocallyOwnedCell());
 * for (const auto &cell : filtered_iterators_range)
 *   {
 *     fe_values.reinit (cell);
 *     ...do the local integration on 'cell'...;
 *   }
 * @endcode
 *
 * @relatesalso FilteredIterator
 *
 * @ingroup CPP11
 *
 *
 */
template <typename BaseIterator, typename Predicate>
IteratorRange<FilteredIterator<BaseIterator>>
filter_iterators(IteratorRange<BaseIterator> i, const Predicate &p)
{
  FilteredIterator<BaseIterator> fi(p, *(i.begin()));
  FilteredIterator<BaseIterator> fi_end(p, *(i.end()));

  return IteratorRange<FilteredIterator<BaseIterator>>(fi, fi_end);
}



/**
 * 通过任意数量的Predicates过滤给定范围的迭代器。这允许替换。
 *
 * @code
 * DoFHandler<dim> dof_handler;
 * ...
 * for (const auto &cell : dof_handler.active_cell_iterators())
 *   {
 *     if (cell->is_locally_owned())
 *       {
 *         if (cell->at_boundary())
 *           {
 *             fe_values.reinit (cell);
 *             ...do the local integration on 'cell'...;
 *           }
 *       }
 *   }
 * @endcode
 * 由。
 *
 * @code
 * DoFHandler<dim> dof_handler;
 * ...
 * const auto filtered_iterators_range =
 *   filter_iterators(dof_handler.active_cell_iterators(),
 *                    IteratorFilters::LocallyOwnedCell(),
 *                    IteratorFilters::AtBoundary());
 * for (const auto &cell : filter_iterators_range)
 *   {
 *     fe_values.reinit (cell);
 *     ...do the local integration on 'cell'...;
 *   }
 * @endcode
 *
 * @relatesalso FilteredIterator
 *
 * @ingroup CPP11
 *
 *
 */
template <typename BaseIterator, typename Predicate, typename... Targs>
IteratorRange<
  typename internal::FilteredIteratorImplementation::
    NestFilteredIterators<BaseIterator, std::tuple<Predicate, Targs...>>::type>
filter_iterators(IteratorRange<BaseIterator> i,
                 const Predicate &           p,
                 const Targs... args)
{
  // Recursively create filtered iterators, one predicate at a time
  auto fi = filter_iterators(i, p);
  return filter_iterators(fi, args...);
}


 /* ------------------ Inline functions and templates ------------ */ 


template <typename BaseIterator>
template <typename Predicate>
inline FilteredIterator<BaseIterator>::FilteredIterator(Predicate p)
  : predicate(new PredicateTemplate<Predicate>(p))
{}



template <typename BaseIterator>
template <typename Predicate>
inline FilteredIterator<BaseIterator>::FilteredIterator(Predicate           p,
                                                        const BaseIterator &bi)
  : BaseIterator(bi)
  , predicate(new PredicateTemplate<Predicate>(p))
{
  if ((this->state() == IteratorState::valid) && !(*predicate)(*this))
    set_to_next_positive(bi);
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>::FilteredIterator(
  const FilteredIterator &fi)
  : // this construction looks strange, but without going through the
    // address of fi, GCC would not cast fi to the base class of type
    // BaseIterator but tries to go through constructing a new
    // BaseIterator with an Accessor.
  BaseIterator(*static_cast<const BaseIterator *>(&fi))
  , predicate(fi.predicate->clone())
{}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator=(const FilteredIterator &fi)
{
  // Using equivalent code to the one for 'operator=(const BaseIterator &bi)'
  // below, some compiler would not cast fi to the base class of type
  // BaseIterator but try to go through constructing a new Accessor from fi
  // which fails. Hence, we just use an explicit upcast and call the above-
  // mentioned method.
  const BaseIterator &bi      = fi;
  return              operator=(bi);
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator=(const BaseIterator &bi)
{
  Assert((bi.state() != IteratorState::valid) || (*predicate)(bi),
         ExcInvalidElement(bi));
  BaseIterator::operator=(bi);
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::set_to_next_positive(const BaseIterator &bi)
{
  BaseIterator::operator=(bi);
  while ((this->state() == IteratorState::valid) && (!(*predicate)(*this)))
    BaseIterator::operator++();

  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::set_to_previous_positive(const BaseIterator &bi)
{
  BaseIterator::operator=(bi);
  while ((this->state() == IteratorState::valid) && (!(*predicate)(*this)))
    BaseIterator::operator--();

  return *this;
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator==(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) ==
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator!=(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) !=
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator<(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) <
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator==(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) == bi);
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator!=(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) != bi);
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator<(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) < bi);
}


template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator++()
{
  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator++();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::operator++(int)
{
  const FilteredIterator old_state = *this;

  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator++();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return old_state;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator--()
{
  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator--();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::operator--(int)
{
  const FilteredIterator old_state = *this;

  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator--();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return old_state;
}



template <typename BaseIterator>
template <typename Predicate>
inline FilteredIterator<BaseIterator>::PredicateTemplate<
  Predicate>::PredicateTemplate(const Predicate &predicate)
  : predicate(predicate)
{}



template <typename BaseIterator>
template <typename Predicate>
bool
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::
operator()(const BaseIterator &bi) const
{
  return predicate(bi);
}



template <typename BaseIterator>
template <typename Predicate>
std::unique_ptr<typename FilteredIterator<BaseIterator>::PredicateBase>
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::clone() const
{
  return std::make_unique<PredicateTemplate>(predicate);
}



namespace IteratorFilters
{
  // ---------------- IteratorFilters::Active ---------

  template <class Iterator>
  inline bool
  Active::operator()(const Iterator &i) const
  {
    return i->is_active();
  }


  // ---------------- IteratorFilters::UserFlagSet ---------

  template <class Iterator>
  inline bool
  UserFlagSet::operator()(const Iterator &i) const
  {
    return (i->user_flag_set());
  }


  // ---------------- IteratorFilters::UserFlagNotSet ---------

  template <class Iterator>
  inline bool
  UserFlagNotSet::operator()(const Iterator &i) const
  {
    return (!i->user_flag_set());
  }


  // ---------------- IteratorFilters::LevelEqualTo ---------
  inline LevelEqualTo::LevelEqualTo(const unsigned int level)
    : level(level)
  {}



  template <class Iterator>
  inline bool
  LevelEqualTo::operator()(const Iterator &i) const
  {
    return (static_cast<unsigned int>(i->level()) == level);
  }



  // ---------------- IteratorFilters::SubdomainEqualTo ---------
  inline SubdomainEqualTo::SubdomainEqualTo(
    const types::subdomain_id subdomain_id)
    : subdomain_id(subdomain_id)
  {}



  template <class Iterator>
  inline bool
  SubdomainEqualTo::operator()(const Iterator &i) const
  {
    return (i->subdomain_id() == subdomain_id);
  }



  // ---------------- IteratorFilters::LocallyOwnedCell ---------

  template <class Iterator>
  inline bool
  LocallyOwnedCell::operator()(const Iterator &i) const
  {
    return (i->is_locally_owned());
  }


  // ---------------- IteratorFilters::LocallyOwnedLevelCell ---------

  template <class Iterator>
  inline bool
  LocallyOwnedLevelCell::operator()(const Iterator &i) const
  {
    return (i->is_locally_owned_on_level());
  }



  // ---------------- IteratorFilters::MaterialIdEqualTo ---------
  inline MaterialIdEqualTo::MaterialIdEqualTo(
    const types::material_id material_id,
    const bool               only_locally_owned)
    : material_ids{material_id}
    , only_locally_owned(only_locally_owned)
  {}



  inline MaterialIdEqualTo::MaterialIdEqualTo(
    const std::set<types::material_id> &material_ids,
    const bool                          only_locally_owned)
    : material_ids(material_ids)
    , only_locally_owned(only_locally_owned)
  {}



  template <class Iterator>
  inline bool
  MaterialIdEqualTo::operator()(const Iterator &i) const
  {
    return only_locally_owned == true ?
             (material_ids.find(i->material_id()) != material_ids.end() &&
              i->is_locally_owned()) :
             material_ids.find(i->material_id()) != material_ids.end();
  }



  // ---------------- IteratorFilters::ActiveFEIndexEqualTo ---------
  inline ActiveFEIndexEqualTo::ActiveFEIndexEqualTo(
    const unsigned int active_fe_index,
    const bool         only_locally_owned)
    : active_fe_indices{active_fe_index}
    , only_locally_owned(only_locally_owned)
  {}



  inline ActiveFEIndexEqualTo::ActiveFEIndexEqualTo(
    const std::set<unsigned int> &active_fe_indices,
    const bool                    only_locally_owned)
    : active_fe_indices(active_fe_indices)
    , only_locally_owned(only_locally_owned)
  {}



  template <class Iterator>
  inline bool
  ActiveFEIndexEqualTo::operator()(const Iterator &i) const
  {
    return only_locally_owned == true ?
             (active_fe_indices.find(i->active_fe_index()) !=
                active_fe_indices.end() &&
              i->is_locally_owned()) :
             active_fe_indices.find(i->active_fe_index()) !=
               active_fe_indices.end();
  }



  // ---------------- IteratorFilters::AtBoundary ---------

  template <class Iterator>
  inline bool
  AtBoundary::operator()(const Iterator &i) const
  {
    return (i->at_boundary());
  }
} // namespace IteratorFilters


DEAL_II_NAMESPACE_CLOSE

 /*------------------------- filtered_iterator.h ------------------------*/ 
#endif
 /*------------------------- filtered_iterator.h ------------------------*/ 


