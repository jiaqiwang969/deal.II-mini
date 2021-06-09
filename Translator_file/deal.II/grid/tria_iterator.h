//include/deal.II-translator/grid/tria_iterator_0.txt
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

#ifndef dealii_tria_iterator_h
#  define dealii_tria_iterator_h


 /*----------------------------   tria-iterator.h ---------------------------*/ 


#  include <deal.II/base/config.h>

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/point.h>

#  include <deal.II/grid/tria_iterator_base.h>

#  include <iterator>
#  include <ostream>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
template <int dim, int spacedim>
class Triangulation;
template <int, int, int>
class TriaAccessorBase;

template <typename>
class TriaIterator;
template <typename>
class TriaActiveIterator;
#  endif



// note: in non-debug mode, i.e. with optimizations, the file
// tria_iterator.templates.h is included at the end of this file.
// this includes a lot of templates and thus makes compilation
// slower, but at the same time allows for more aggressive
// inlining and thus faster code.


/**
 * 这个类实现了一个迭代器，类似于标准库中使用的那些迭代器。它满足了双向迭代器的要求。关于迭代器规格和用法的进一步细节，请参见C++文档。
 *
 *  除了标准接口之外，这个类的迭代器还提供了一个<tt>-
 * @></tt> 操作符，也就是说，你可以写这样的语句
 *
 * @code
 * cell->set_refine_flag ();
 * @endcode
 *
 * 每当要对所有的线、四边形、单元格等进行循环时，就会用到迭代器。这些循环可以像这样编码。
 *
 * @code
 * cell_iterator cell = tria.begin();
 * cell_iterator end  = tria.end();
 * for (; cell!=end; ++cell)
 * if (cell->at_boundary())
 *   cell->set_refine_flag();
 * @endcode
 *
 * 注意使用<tt>++cell</tt>而不是<tt>cell++</tt>，因为这不涉及暂存器和复制。建议在循环中使用一个固定值<tt>end</tt>而不是<tt>tria.end()</tt>，因为与普通的指针相比，这些迭代器的创建和复制相当昂贵。
 * 指向的对象是访问器，来自于TriaAccessorBase。哪种访问器由模板参数决定
 * <em>  访问器  </em>
 * 。这些访问器与其说是数据结构，不如说是提供访问存储在Triangulation或DoFHandler对象中的数据的函数集合。使用这些访问器，这些类的结构被隐藏在应用程序中。
 * <h3>Which iterator to use when</h3>
 * @attention
 * 应用程序很少使用TriaRawIterator，而是使用派生类TriaIterator或TriaActiveIterator中的一个。
 * <ul>   <li>  TriaRawIterator对象指向列表中的行、单元格等，无论它们是否被使用（在向量中，也存储<i>dead</i>对象，因为向量中的删除很昂贵，我们也不希望破坏向量中的编号所引起的排序）。因此不是所有的原始迭代器都指向有效的对象。
 * <li>
 * 派生类TriaIterator选择有效的单元，也就是在三角形层次结构中某个地方使用的单元。
 * <li>  TriaActiveIterator对象，只在活动单元上循环。  </ul>
 * <h3>Purpose</h3>
 * 迭代器并不比直接在数据结构上操作慢多少，因为它们执行的是你必须自己手工编码的循环。大多数迭代器和访问器函数是内联的。
 * 迭代器的主要功能在于<tt>++</tt>和<tt>--</tt>操作符。这些操作将迭代器向前或向后移动，就像它是一个指向数组的指针一样。在这里，这个操作并不容易，因为它可能包括跳过一些元素和三角化水平之间的过渡。这对用户来说是完全隐藏的，尽管你仍然可以创建一个指向任意元素的迭代器。
 * 实际上，来回移动迭代器的操作不是在迭代器类中完成的，而是在访问器类中完成的。因为这些是作为模板参数传递的，所以你可以在这里写出你自己的版本来增加更多的功能。
 * 此外，这里描述的迭代器满足了C++标准所规定的输入和双向迭代器的要求。因此，可以使用C++标准中算法部分的函数，例如，
 * <em>  count_if  </em>
 * （见Triangulation的文档中的一个例子）和其他几个。
 * <h3>Implementation</h3>
 * 迭代器类本身并没有什么功能。它只有在分配给Accessor（第二个模板参数）时才变得有用，Accessor才真正对数据进行访问。一个Accessor必须满足一些要求。
 * <ul>   <li>  它必须有两个名为<tt>present_level</tt>和<tt>present_index</tt>的成员，存储目前所指向的三角结构中的元素的地址。这些数据必须能被上面列出的所有三角形迭代器访问。
 * <li>
 * 它必须有一个构造函数，接收一个三角形*和两个无符号整数，表示初始级别和索引，以及一个取决于其类型的数据对象。
 * <li>
 * 对于TriaIterator和TriaActiveIterator类，它必须有一个成员函数<tt>bool
 * used()</tt>，对于后者有一个成员函数<tt>bool active()</tt>。
 * <li>  它必须有无效操作符<tt>++</tt>和<tt>--</tt>。
 * <li>  它必须声明一个本地别名<tt>AccessorData</tt>，说明访问器期望得到的数据类型作为第四个构造参数传递。通过声明本地数据类型，相应的迭代器类可以类型安全地强制该数据类型为它自己的构造器参数类型之一。如果一个存取器类不需要额外的数据，这个类型应是<tt>void</tt>。  </ul>
 * 那么迭代器就能够做它应该做的事情。所有必要的函数都在<tt>Accessor</tt>基类中实现，但你可以编写自己的版本（非虚拟的，因为我们使用模板）来增加功能。
 * 该库提供的访问器由两组组成，由它们是访问Triangulation对象还是DoFHandler对象的数据决定。它们分别来自TriaAccessor和DoFAccessor。每组都有专门的单元格（相对于面和线）的访问器，提供更多的功能，如访问邻居。
 * @attention
 * 似乎不可能通过迭代器的使用来保持三角形的不变性。因此，如果你声明指向<tt>const</tt>
 * triangulation对象的指针，你应该清楚地意识到你可能会不由自主地改变存储在三角结构中的数据。
 *
 *
 * @note
 * 关于有效和无效的迭代器的更多信息可以在TriaAccessorBase的文档中找到，其中检查和实现了迭代器状态。
 *
 *  <h3>Past-the-end iterators</h3>
 * 有一种past-the-end-pointers的表示方法，由成员变量的特殊值表示
 * @p present_level  和  @p present_index:  ]
 * 如果<tt>present_level>=0</tt>和<tt>present_index>=0</tt>，那么这个对象是有效的（当我们调查一个迭代器的状态时，没有检查三角形是否真的有那么多层或者在本层有那么多单元。然而，在许多迭代器被取消引用的地方，我们会进行这样的检查）；如果<tt>present_level==-1</tt>和<tt>present_index==-1</tt>，那么迭代器就指向了终点；在所有其他情况下，迭代器被视为无效。你可以通过调用<tt>state()</tt>函数来检查这一点。
 * 如果指向Triangulation对象的指针无效或为零，那么一个迭代器也是无效的。
 * 最后，如果 @p  present_level和 @p present_index
 * 所指向的元素没有被使用，即 @p used
 * 标志被设置为false，则迭代器是无效的。
 * 后面两个检查不在<tt>state()</tt>中进行，因为这两种情况应该只在通过
 * @p memcpy
 * 和类似的未初始化构造时发生（父三角形只能在构造时被设置）。如果一个迭代器通过空构造函数被构造为空，<tt>present_level==-2</tt>和<tt>present_index==-2</tt>。因此，无论如何，这个迭代器是无效的，不管三角指针的状态和所指向的元素的状态如何。
 * 过去的迭代器也可以用来比较一个迭代器和<i>before-the-start</i>的值，当向后运行时。指向一个向量两端过去的迭代器之间没有区别。
 * 通过只定义一个值为过端，并使所有其他的值无效，提供了第二个安全轨道：如果我们应该在迭代器被递增或递减时忘记了库中的检查，我们自动将迭代器从允许的状态
 * "过端 "转换为不允许的状态
 * "无效"，这增加了比过端迭代器早一些时间引发异常的机会。
 * @ref Triangulation
 * @ingroup grid
 *
 * @ingroup Iterators
 *
 */
template <typename Accessor>
class TriaRawIterator
{
public:
  /**
   * 声明Accessor的类型，以便在外部世界使用。这样其他函数就可以使用Accessor的类型而不知道具体的实现方式。
   *
   */
  using AccessorType = Accessor;

  /**
   * 默认构造函数。这个构造函数创建一个指向无效对象的迭代器。因此，该迭代器是不可用的。
   *
   */
  TriaRawIterator();

  /**
   * 复制构造函数。
   *
   */
  TriaRawIterator(const TriaRawIterator &);

  /**
   * 从给定的访问器中构造一个迭代器；给定的访问器不需要与本类的访问器是同一类型，但它需要是可转换的。
   * 通过这个构造函数，也可以构造派生迭代器的对象。
   * @code
   * DoFCellAccessor dof_accessor;
   * Triangulation::active_cell_iterator cell = dof_accessor;
   * @endcode
   *
   *
   */
  explicit TriaRawIterator(const Accessor &a);

  /**
   * 构造器。假设其他访问器类型可以转换为当前的类型。
   *
   */
  template <typename OtherAccessor>
  explicit TriaRawIterator(const OtherAccessor &a);

  /**
   * 适当的构造函数，用三角形、水平和指向的对象的索引进行初始化。最后一个参数是访问器类所声明的类型。
   *
   */
  TriaRawIterator(
    const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
    const int                                                            level,
    const int                                                            index,
    const typename AccessorType::AccessorData *local_data = nullptr);

  /**
   * 这是一个转换操作符（构造器），它接受另一个迭代器类型并复制数据；如果有一个从
   * @p OtherAccessor 类到这个对象的 @p Accessor
   * 类的转换路径，那么这个转换就能发挥作用。一个这样的路径是派生类到基类，例如可以用来从
   * DoFHandler::raw_cell_iterator, 中获得 Triangulation::raw_cell_iterator
   * ，因为DoFAccessor类是派生自TriaAccessorBase类。
   *
   */
  template <typename OtherAccessor>
  TriaRawIterator(const TriaRawIterator<OtherAccessor> &i);

  /**
   * 另一个转换操作符，我们使用来自TriaAccessorBase对象的Triangulation指针，同时根据Accessor的实际类型使用附加数据。
   *
   */
  TriaRawIterator(
    const TriaAccessorBase<Accessor::structure_dimension,
                           Accessor::dimension,
                           Accessor::space_dimension> &tria_accessor,
    const typename Accessor::AccessorData *            local_data);

  /**
   * 转换构造器。与上述相同，不同的是它从TriaIterator类（而不是TriaRawIterator）进行转换。
   *
   */
  template <typename OtherAccessor>
  TriaRawIterator(const TriaIterator<OtherAccessor> &i);

  /**
   * 转换构造函数。和上面一样，不同的是它转换自TriaActiveIterator类（而不是TriaRawIterator）。
   *
   */
  template <typename OtherAccessor>
  TriaRawIterator(const TriaActiveIterator<OtherAccessor> &i);

  /**
   * @name  解除引用
   *
   */
   /*@{*/ 
  /**
   * 解除引用操作符，返回对一个访问器的引用。因此，用法类似于<tt>(*i).index
   * ();</tt> 这个函数必须为不同的 @p
   * 指针进行明确的专业化处理，以允许一个
   * <tt>iterator<1,TriangulationLevel<1>::LinesData></tt>
   * 指向<tt>tria->lines.cell[index]</tt>，而对于高一维，则必须指向<tt>tria->quads.cell[index]</tt>。
   * 你不能解读无效的或超过终点的迭代器。
   *
   */
  const Accessor &operator*() const;

  /**
   * 去引用操作符，非 @p const 版本。
   *
   */
  Accessor &operator*();

  /**
   * 去引用操作符，返回指向的单元格的引用。使用方法类似于
   * <tt>i->index ();</tt> 有一个  @p const  和一个非  @p const
   * 版本。
   *
   */
  const Accessor *operator->() const;

  /**
   * 去引用操作符，非 @p const 版本。
   *
   */
  Accessor *operator->();


  /**
   * 为了能够将不同的访问器的末尾迭代器相互分配，我们需要一个访问函数，无论访问器的状态如何，都会返回访问器。
   * @warning
   * 这个函数不应该在应用程序中使用。它只用于库内的有限用途，而且会使调试工作变得更加困难。
   *
   */
  const Accessor &
  access_any() const;

   /*@}*/ 

  /**
   * 赋值运算符。
   *
   */
  TriaRawIterator &
  operator=(const TriaRawIterator &);

  /**
   * 等价比较。
   *
   */
  template <typename OtherAccessor = Accessor>
  typename std::enable_if<std::is_convertible<OtherAccessor, Accessor>::value,
                          bool>::type
  operator==(const TriaRawIterator<OtherAccessor> &) const;

  /**
   * 不等式比较。
   *
   */
  bool
  operator!=(const TriaRawIterator &) const;

  /**
   * 迭代器的排序关系。    这种关系试图对单元格进行总排序。    该关系定义如下。    对于 <tt>Accessor::structure_dimension < Accessor::dimension</tt>, 的对象，我们只需比较这种对象的索引。  排序是根据以下的层次结构进行的（在这个意义上，只有当前一个测试没有结果时，才会应用下一个测试）。      <ol>   <li>  过去结束的迭代器总是最后排序。两个过去-结束的迭代器排名相同，因此在这种情况下会返回false。 </li>   <li>  单元的级别。 </li>   <li>  级别内的单元格的索引。 </li>   </ol>  级别内的单元格。
   * @note  在一个 parallel::distributed::Triangulation
   * 中，不同的处理器之间的排序是不一致的，因为我们依靠index()，这很可能是不一样的。
   *
   */
  bool
  operator<(const TriaRawIterator &) const;

  /**
   * 另一个比较运算符，实现与#operator<相同的排序。
   *
   */
  bool
  operator>(const TriaRawIterator &) const;

  /**
   * @name  推进迭代器的发展
   *
   */
   /*@{*/ 
  /**
   * 前缀<tt>++</tt>操作符。<tt>++iterator</tt>。这个操作符将迭代器推进到下一个元素，并返回一个对<tt>*this</tt>的引用。
   *
   */
  TriaRawIterator &
  operator++();

  /**
   * 后缀<tt>++</tt>操作符。<tt>iterator++</tt>。该操作符将迭代器推进到下一个元素，但返回之前指向的元素的迭代器。
   * 由于这个操作涉及到一个临时和复制操作，而且对于一个指针来说，
   * @p iterator
   * 是一个相当大的对象，所以尽可能使用前缀操作符<tt>++iterator</tt>，特别是在for循环的头部（<tt>for
   * (; iterator!=end;
   * ++iterator)</tt>），因为在那里你通常不需要返回值。
   *
   */
  TriaRawIterator
  operator++(int);

  /**
   * 前缀  @p \--  操作符。  @p \--iterator.
   * 这个操作符将迭代器移动到前一个元素，并返回一个对<tt>*this</tt>的引用。
   *
   */
  TriaRawIterator &
  operator--();

  /**
   * 后缀  @p \--  操作符。  @p iterator\--.
   * 该操作符将迭代器移动到前一个元素，但返回之前指向的元素的迭代器。
   * 这与后缀操作符++的情况相同：如果可能的话，使用前缀操作符的形式来避免它，以避免使用临时变量。
   *
   */
  TriaRawIterator
  operator--(int);
   /*@}*/ 

  /**
   * 返回迭代器的状态。
   *
   */
  IteratorState::IteratorStates
  state() const;

  /**
   * 打印迭代器到一个流  <code>out</code>
   * 。格式是<tt>level.index</tt>。
   *
   */
  template <class StreamType>
  void
  print(StreamType &out) const;


  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 将该类标记为双向迭代器，并声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
   *
   */
  using iterator_category = std::bidirectional_iterator_tag;
  using value_type        = Accessor;
  using difference_type   = int;
  using pointer           = Accessor *;
  using reference         = Accessor &;

  /**
   * @name  异常情况
   *
   */
   /*@{*/ 
  /**
   * 有水平的TriaObjects的异常，即单元格。
   *
   */
  DeclException1(ExcDereferenceInvalidCell,
                 Accessor,
                 << "You tried to dereference a cell iterator for which this "
                 << "is not possible. More information on this iterator: "
                 << "level=" << arg1.level() << ", index=" << arg1.index()
                 << ", state="
                 << (arg1.state() == IteratorState::valid ?
                       "valid" :
                       (arg1.state() == IteratorState::past_the_end ?
                          "past_the_end" :
                          "invalid")));

  /**
   * 异常情况。
   *
   */
  DeclException1(ExcDereferenceInvalidObject,
                 Accessor,
                 << "You tried to dereference an iterator for which this "
                 << "is not possible. More information on this iterator: "
                 << "index=" << arg1.index() << ", state="
                 << (arg1.state() == IteratorState::valid ?
                       "valid" :
                       (arg1.state() == IteratorState::past_the_end ?
                          "past_the_end" :
                          "invalid")));

  /**
   * 异常情况
   *
   */
  DeclException0(ExcAdvanceInvalidObject);
  /**
   * 异常情况
   *
   */
  DeclException0(ExcInvalidComparison);

   /*@}*/ 
protected:
  /**
   * 持有真实数据的对象。
   *
   */
  Accessor accessor;


  // Make all other iterator class templates friends of this class. This is
  // necessary for the implementation of conversion constructors.
  //
  // In fact, we would not need them to be friends if they were for different
  // dimensions, but the compiler dislikes giving a fixed dimension and
  // variable accessor since then it says that would be a partial
  // specialization.
  template <typename SomeAccessor>
  friend class TriaRawIterator;
  template <typename SomeAccessor>
  friend class TriaIterator;
  template <typename SomeAccessor>
  friend class TriaActiveIterator;
};


/**
 * TriaRawIterator的这种特殊化只提供对 <em> 使用的 </em>
 * 线、四边形、单元格等的访问。
 *
 *
 * @ingroup grid
 *
 * @ingroup Iterators
 *
 */
template <typename Accessor>
class TriaIterator : public TriaRawIterator<Accessor>
{
public:
  /**
   * 默认构造函数。这个构造函数创建一个指向无效对象的迭代器。因此，该迭代器是不可用的。
   *
   */
  TriaIterator();

  /**
   * 复制构造函数。
   *
   */
  TriaIterator(const TriaIterator<Accessor> &);

  /**
   * 从可能指向非活动对象的迭代器的转换构造器（即，对于那些我们不能仅通过查看其类型就知道该对象被使用的对象）。
   * @pre
   * 传递给这个构造函数的参数必须是(i)一个过去的迭代器；或者(ii)它必须指向一个已使用的对象。所有其他情况都会导致异常。
   *
   */
  TriaIterator(const TriaRawIterator<Accessor> &);

  /**
   * 构造函数，用三角形、水平和指向的对象的索引进行初始化。最后一个参数是由访问器类声明的类型。
   * @pre
   * 传递给这个构造函数的参数必须是（i）一个过去的迭代器；或者（ii）它必须指向一个使用过的对象。所有其他情况都会导致异常。
   *
   */
  TriaIterator(
    const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
    const int                                                            level,
    const int                                                            index,
    const typename Accessor::AccessorData *local_data = nullptr);

  /**
   * 从一个可转换为Accessor类型的OtherAccessor的访问器构造。
   *
   */
  template <typename OtherAccessor>
  explicit TriaIterator(const OtherAccessor &a);

  /**
   * 这是一个转换操作符（构造器），它接收另一个迭代器类型并复制数据；如果有一个从
   * @p OtherAccessor 类到该对象的 @p Accessor
   * 类的转换路径，那么这种转换就能发挥作用。一个这样的路径是派生类到基类，例如可以用来从
   * DoFHandler::cell_iterator, 中获得 Triangulation::cell_iterator
   * ，因为DoFAccessor类是派生自TriaAccessorBase类。
   *
   */
  template <typename OtherAccessor>
  TriaIterator(const TriaIterator<OtherAccessor> &i);

  /**
   * 另一个转换操作符，我们使用TriaAccessorBase对象的三角函数的指针，同时根据Accessor的实际类型使用附加数据。
   *
   */
  TriaIterator(const TriaAccessorBase<Accessor::structure_dimension,
                                      Accessor::dimension,
                                      Accessor::space_dimension> &tria_accessor,
               const typename Accessor::AccessorData *            local_data);

  /**
   * 类似于上面的转换操作符，但做了一个检查迭代器是否指向一个使用的元素，这对原始迭代器来说是必要的。
   *
   */
  template <typename OtherAccessor>
  TriaIterator(const TriaRawIterator<OtherAccessor> &i);

  /**
   * 与上面那个类似的转换操作符，但用于从活动迭代器的转换。
   *
   */
  template <typename OtherAccessor>
  TriaIterator(const TriaActiveIterator<OtherAccessor> &i);

  /**
   * 赋值运算符。
   *
   */
  TriaIterator<Accessor> &
  operator=(const TriaIterator<Accessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个使用的元素时才有效。
   *
   */
  TriaIterator<Accessor> &
  operator=(const TriaRawIterator<Accessor> &);

  /**
   * 赋值运算符。要求，Accessor可以从OtherAccessor复制。
   *
   */
  template <class OtherAccessor>
  TriaIterator<Accessor> &
  operator=(const TriaIterator<OtherAccessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个使用的元素时才有效。要求，Accessor可以从OtherAccessor中复制。
   *
   */
  template <class OtherAccessor>
  TriaIterator<Accessor> &
  operator=(const TriaRawIterator<OtherAccessor> &);

  /**
   * @name  迭代器的推进
   *
   */
   /*@{*/ 
  /**
   * 前缀<tt>++</tt>操作符。<tt>++i</tt>。这个操作符将迭代器推进到下一个使用的元素，并返回一个对<tt>*this</tt>的引用。
   *
   */
  TriaIterator<Accessor> &
  operator++();

  /**
   * 后缀<tt>++</tt>操作符。<tt>i++</tt>。这个操作符将迭代器推进到下一个使用的元素，但返回一个迭代器到之前指向的元素。由于这涉及到一个临时和复制操作，而且对于一个指针来说，
   * @p active_iterator
   * 是一个相当大的对象，所以尽可能使用前缀操作符<tt>++i</tt>，特别是在for循环的头部（<tt>for
   * (; i!=end; ++i)</tt>），因为在那里你通常不需要返回值。
   *
   */
  TriaIterator<Accessor>
  operator++(int);

  /**
   * 前缀  @p \--  操作符。  @p \--i.
   * 这个操作符将迭代器推进到前一个使用的元素，并返回对<tt>*this</tt>的引用。
   *
   */
  TriaIterator<Accessor> &
  operator--();

  /**
   * Postfix  @p \--  操作符。  @p i\--.
   *
   */
  TriaIterator<Accessor>
  operator--(int);
   /*@}*/ 

  /**
   * 声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
   *
   */
  using iterator_category =
    typename TriaRawIterator<Accessor>::iterator_category;
  using value_type      = typename TriaRawIterator<Accessor>::value_type;
  using pointer         = typename TriaRawIterator<Accessor>::pointer;
  using reference       = typename TriaRawIterator<Accessor>::reference;
  using difference_type = typename TriaRawIterator<Accessor>::difference_type;

  /**
   * 异常情况
   *
   */
  DeclException0(ExcAssignmentOfUnusedObject);
};


/**
 * TriaIterator的这种特殊化只提供对 <em> 活动 </em>
 * 线、四边形、单元格等的访问。活动单元是指没有被细化的单元，因此是在最细级别上进行计算的单元。
 *
 *
 * @ingroup grid
 * @ingroup Iterators
 *
 */
template <typename Accessor>
class TriaActiveIterator : public TriaIterator<Accessor>
{
public:
  /**
   * 默认构造函数。这个构造函数创建一个指向无效对象的迭代器。因此，该迭代器是不可用的。
   *
   */
  TriaActiveIterator();

  /**
   * 复制构造函数。
   *
   */
  TriaActiveIterator(const TriaActiveIterator<Accessor> &);

  /**
   * 转换构造函数从指向潜在非活动对象的迭代器（或至少从其类型上看不出它是活动的）创建一个活动的迭代器。
   * @pre
   * 传递给这个构造函数的参数必须是(i)一个已过终点的迭代器；或者(ii)它必须指向一个活动对象。所有其他情况都会导致异常。
   *
   */
  TriaActiveIterator(const TriaRawIterator<Accessor> &);

  /**
   * 转换构造器从一个指向潜在非活动对象的迭代器中创建一个活动迭代器（或者至少从其类型上看不出它是活动的）。
   * @pre
   * 传递给这个构造函数的参数必须是(i)一个已过终点的迭代器；或者(ii)它必须指向一个活动对象。所有其他情况都会导致异常。
   *
   */
  TriaActiveIterator(const TriaIterator<Accessor> &);

  /**
   * 构造函数，用三角形、级别和指向的对象的索引进行初始化。最后一个参数是由当前迭代器使用的访问器类所声明的类型。
   * @pre
   * 传递给这个构造函数的参数必须是(i)一个已过终点的迭代器；或者(ii)它必须指向一个活动对象。所有其他的情况将导致一个异常。
   *
   */
  TriaActiveIterator(
    const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
    const int                                                            level,
    const int                                                            index,
    const typename Accessor::AccessorData *local_data = nullptr);

  /**
   * 这是一个转换操作符（构造函数），它接收另一个迭代器类型并复制数据；如果有一个从
   * @p OtherAccessor 类到该对象的 @p Accessor
   * 类的转换路径，那么这种转换就能发挥作用。一个这样的路径是派生类到基类，例如可以用来从
   * DoFHandler::active_cell_iterator, 中获得
   * Triangulation::active_cell_iterator
   * ，因为DoFAccessor类是派生自TriaAccessorBase类。
   *
   */
  template <typename OtherAccessor>
  TriaActiveIterator(const TriaActiveIterator<OtherAccessor> &i);

  /**
   * 另一个转换操作符，我们使用TriaAccessorBase对象的三角函数的指针，同时根据Accessor的实际类型使用附加数据。
   *
   */
  TriaActiveIterator(
    const TriaAccessorBase<Accessor::structure_dimension,
                           Accessor::dimension,
                           Accessor::space_dimension> &tria_accessor,
    const typename Accessor::AccessorData *            local_data);

  /**
   * 与上面的转换操作类似，但做了一个检查迭代器是否指向一个使用的元素，并且是活动的，这对原始迭代器来说是必要的。由于通常的迭代器也是原始迭代器，这个构造函数也适用于<tt>TriaIterator<OtherAccessor></tt>类型的参数。
   * @pre
   * 传递给这个构造函数的参数必须是(i)一个过去的迭代器；或者(ii)它必须指向一个活动对象。所有其他情况都会导致异常。
   *
   */
  template <typename OtherAccessor>
  TriaActiveIterator(const TriaRawIterator<OtherAccessor> &i);

  /**
   * 赋值运算符。
   *
   */
  TriaActiveIterator<Accessor> &
  operator=(const TriaActiveIterator<Accessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个活动元素时才有效。
   *
   */
  TriaActiveIterator<Accessor> &
  operator=(const TriaIterator<Accessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个活动元素或超过终点时才有效。
   *
   */
  TriaActiveIterator<Accessor> &
  operator=(const TriaRawIterator<Accessor> &);

  /**
   * 赋值运算符。要求，Accessor可以从OtherAccessor中复制。
   *
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator=(const TriaActiveIterator<OtherAccessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个活动元素或超过终点时才有效。要求，Accessor可以从OtherAccessor中复制。
   *
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator=(const TriaRawIterator<OtherAccessor> &);

  /**
   * 交叉赋值运算符。这个赋值只有在给定的迭代器指向一个活动元素时才有效。要求Accessor可以从OtherAccessor中复制。
   *
   */
  template <class OtherAccessor>
  TriaActiveIterator<Accessor> &
  operator=(const TriaIterator<OtherAccessor> &);

  /**
   * 前缀<tt>++</tt>操作符。<tt>++i</tt>。这个操作符将迭代器推进到下一个活动元素，并返回一个对<tt>*this</tt>的引用。
   *
   */
  TriaActiveIterator<Accessor> &
  operator++();

  /**
   * @name  迭代器的前进
   *
   */
   /*@{*/ 
  /**
   * 后缀<tt>++</tt>操作符。<tt>i++</tt>。该操作符将迭代器推进到下一个活动元素，但返回之前指向的元素的迭代器。由于这涉及到一个临时和复制操作，而且对于一个指针来说，
   * @p active_iterator
   * 是一个相当大的对象，所以尽可能使用前缀操作符<tt>++i</tt>，特别是在for循环的头部（<tt>for
   * (; i!=end; ++i)</tt>），因为在那里你通常不需要返回值。
   *
   */
  TriaActiveIterator<Accessor>
  operator++(int);

  /**
   * 前缀  @p \--  操作符。  @p \--i.
   * 这个操作符将迭代器推进到前一个活动元素，并返回一个对<tt>*this</tt>的引用。
   *
   */
  TriaActiveIterator<Accessor> &
  operator--();

  /**
   * Postfix  @p \--  操作符。  @p i\--.
   *
   */
  TriaActiveIterator<Accessor>
  operator--(int);
   /*@}*/ 

  /**
   * 声明一些别名，这些别名是迭代器的标准，被算法用来查询它们所工作的迭代器的具体内容。
   *
   */
  using iterator_category = typename TriaIterator<Accessor>::iterator_category;
  using value_type        = typename TriaIterator<Accessor>::value_type;
  using pointer           = typename TriaIterator<Accessor>::pointer;
  using reference         = typename TriaIterator<Accessor>::reference;
  using difference_type   = typename TriaIterator<Accessor>::difference_type;

  /**
   * 异常情况
   *
   */
  DeclException0(ExcAssignmentOfInactiveObject);
};


 /*----------------------- Inline functions -------------------*/ 


template <typename Accessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(const Accessor &a)
  : accessor(a)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(const OtherAccessor &a)
  : accessor(a)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const TriaRawIterator<OtherAccessor> &i)
  : accessor(i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const TriaIterator<OtherAccessor> &i)
  : accessor(i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const TriaActiveIterator<OtherAccessor> &i)
  : accessor(i.accessor)
{}



template <typename Accessor>
inline const Accessor &TriaRawIterator<Accessor>::operator*() const
{
  Assert(Accessor::structure_dimension != Accessor::dimension ||
           state() == IteratorState::valid,
         ExcDereferenceInvalidCell(accessor));
  Assert(Accessor::structure_dimension == Accessor::dimension ||
           state() == IteratorState::valid,
         ExcDereferenceInvalidObject(accessor));

  return accessor;
}



template <typename Accessor>
inline Accessor &TriaRawIterator<Accessor>::operator*()
{
  Assert(Accessor::structure_dimension != Accessor::dimension ||
           state() == IteratorState::valid,
         ExcDereferenceInvalidCell(accessor));
  Assert(Accessor::structure_dimension == Accessor::dimension ||
           state() == IteratorState::valid,
         ExcDereferenceInvalidObject(accessor));

  return accessor;
}



template <typename Accessor>
inline const Accessor &
TriaRawIterator<Accessor>::access_any() const
{
  return accessor;
}



template <typename Accessor>
inline const Accessor *TriaRawIterator<Accessor>::operator->() const
{
  return &(this->operator*());
}



template <typename Accessor>
inline Accessor *TriaRawIterator<Accessor>::operator->()
{
  return &(this->operator*());
}



template <typename Accessor>
inline IteratorState::IteratorStates
TriaRawIterator<Accessor>::state() const
{
  return accessor.state();
}



template <typename Accessor>
inline bool
TriaRawIterator<Accessor>::
operator<(const TriaRawIterator<Accessor> &other) const
{
  Assert(state() != IteratorState::invalid,
         ExcDereferenceInvalidObject(accessor));
  Assert(other.state() != IteratorState::invalid,
         ExcDereferenceInvalidObject(other.accessor));

  Assert(&accessor.get_triangulation() == &other.accessor.get_triangulation(),
         ExcInvalidComparison());

  // Deal with iterators past end
  if (state() == IteratorState::past_the_end)
    return false;
  if (other.state() == IteratorState::past_the_end)
    return true;

  return ((**this) < (*other));
}



template <typename Accessor>
inline bool
TriaRawIterator<Accessor>::
operator>(const TriaRawIterator<Accessor> &other) const
{
  return (other < *this);
}



template <typename Accessor>
inline TriaRawIterator<Accessor> &
TriaRawIterator<Accessor>::operator++()
{
  Assert(state() == IteratorState::valid, ExcAdvanceInvalidObject());

  ++accessor;
  return *this;
}



template <typename Accessor>
inline TriaRawIterator<Accessor> &
TriaRawIterator<Accessor>::operator--()
{
  Assert(state() == IteratorState::valid, ExcAdvanceInvalidObject());

  --accessor;
  return *this;
}



template <typename Accessor>
template <class StreamType>
inline void
TriaRawIterator<Accessor>::print(StreamType &out) const
{
  if (Accessor::structure_dimension == Accessor::dimension)
    out << accessor.level() << "." << accessor.index();
  else
    out << accessor.index();
}



template <typename Accessor>
inline std::size_t
TriaRawIterator<Accessor>::memory_consumption() const
{
  return sizeof(TriaRawIterator<Accessor>);
}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaIterator<Accessor>::TriaIterator(
  const TriaIterator<OtherAccessor> &i)
  : TriaRawIterator<Accessor>(i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaIterator<Accessor>::TriaIterator(
  const TriaActiveIterator<OtherAccessor> &i)
  : TriaRawIterator<Accessor>(i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaIterator<Accessor>::TriaIterator(
  const TriaRawIterator<OtherAccessor> &i)
  : TriaRawIterator<Accessor>(i.accessor)
{
#  ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#  endif
}

template <typename Accessor>
template <typename OtherAccessor>
TriaIterator<Accessor>::TriaIterator(const OtherAccessor &a)
  : TriaRawIterator<Accessor>(a)
{
#  ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#  endif
}

template <typename Accessor>
template <typename OtherAccessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaActiveIterator<OtherAccessor> &i)
  : TriaIterator<Accessor>(i.accessor)
{}



template <typename Accessor>
template <typename OtherAccessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaRawIterator<OtherAccessor> &i)
  : TriaIterator<Accessor>(i)
{
#  ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#  endif
}



/**
 * 打印这个迭代器指向的地址  @p out.
 * 这个地址由一对<tt>(level,index)</tt>给出，其中 @p index
 * 是相对于被指向的对象所在的层的索引。
 *
 *
 */
template <typename Accessor>
inline std::ostream &
operator<<(std::ostream &out, const TriaRawIterator<Accessor> &i)
{
  i.print(out);
  return out;
}



/**
 * 打印这个迭代器指向的地址  @p out.
 * 这个地址由一对<tt>(level,index)</tt>给出，其中 @p index
 * 是相对于被指向的对象所在级别的一个索引。
 *
 *
 */
template <typename Accessor>
inline std::ostream &
operator<<(std::ostream &out, const TriaIterator<Accessor> &i)
{
  i.print(out);
  return out;
}



/**
 * 打印这个迭代器指向的地址  @p out.
 * 这个地址由一对<tt>(level,index)</tt>给出，其中 @p index
 * 是相对于被指向的对象所在级别的一个索引。
 *
 *
 */
template <typename Accessor>
inline std::ostream &
operator<<(std::ostream &out, const TriaActiveIterator<Accessor> &i)
{
  i.print(out);
  return out;
}


DEAL_II_NAMESPACE_CLOSE


// if in optimized mode: include more templates
#  ifndef DEBUG
#    include "tria_iterator.templates.h"
#  endif


 /*----------------------------   tria-iterator.h ---------------------------*/ 
#endif
 /*----------------------------   tria-iterator.h ---------------------------*/ 


