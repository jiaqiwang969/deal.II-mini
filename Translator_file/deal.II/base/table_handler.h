//include/deal.II-translator/base/table_handler_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2021 by the deal.II authors
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

#ifndef dealii_table_handler_h
#define dealii_table_handler_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/serialization/map.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <fstream>
#include <map>
#include <ostream>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
class TableHandler;
#endif

namespace internal
{
  /**
   * 一个<tt>TableEntry</tt>存储一个表项的值。它可以是int、无符号int、
   * std::uint64_t, double或 std::string.
   * 类型。实质上，这个结构与 <code>boost::variant@<int,unsigned
   * int,std::uint64_t,double,std::string@></code>
   * 相同，但我们把这个对象包裹在一个结构中，我们可以为它写一个可以序列化的函数。这也是为什么该函数事实上不是
   * boost::any. 类型的原因。
   *
   */
  struct TableEntry
  {
  public:
    /**
     * 默认构造函数。
     *
     */
    TableEntry() = default;

    /**
     * 构造函数。用值初始化这个表元素  <code>t</code>  。
     *
     */
    template <typename T>
    TableEntry(const T &t);

    /**
     * 返回这个对象所存储的值。模板类型T必须是
     * <code>int,unsigned int,std::uint64_t,double,std::string</code>
     * 之一，它必须与这个TableEntry对象中最初存储的对象的数据类型相匹配。
     *
     */
    template <typename T>
    T
    get() const;

    /**
     * 如果数据已经以整数、无符号 integer,std::uint64_t,
     * 或双数形式存储在此对象中，则返回此对象的数值。
     * @return 双数
     *
     */
    double
    get_numeric_value() const;

    /**
     * 用给定的格式缓存包含的值并返回。来自列定义的给定参数被用于格式化。该值在内部以字符串形式缓存在cached_value中。
     * 如果列的格式发生变化，需要用这个例程使缓存失效。
     *
     */
    void
    cache_string(bool scientific, unsigned int precision) const;

    /**
     * 返回使用cache_string()缓存的值。这只是对cached_value的一个包装。
     *
     */
    const std::string &
    get_cached_string() const;


    /**
     * 返回一个TableEntry对象，该对象具有与存储值相同的数据类型，但其值是为该数据类型默认构建的。
     * 这用于在先前设置的列下面垫高。
     *
     */
    TableEntry
    get_default_constructed_copy() const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以达到序列化的目的。
     *
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
     *
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int version);
#else
    // This macro defines the serialize() method that is compatible with
    // the templated save() and load() method that have been implemented.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

  private:
    /**
     * 此对象所存储的数据类型的缩写。
     *
     */
    using value_type =
      boost::variant<int, unsigned int, std::uint64_t, double, std::string>;

    /**
     * 存储的值。
     *
     */
    value_type value;

    /**
     * 缓存当前值为一个字符串。
     *
     */
    mutable std::string cached_value;

    friend class dealii::TableHandler;
  };
} // namespace internal


/**
 * TableHandler存储任意值类型的TableEntries，并将表以文本或tex格式写入输出流。值的类型实际上可以从列到列，从行到行而变化。
 * <h3>Usage</h3> 最重要的函数是模板化的函数
 * <code>add_value(const std::string &key, const T value)</code>
 * ，如果这一列还不存在，就在表中添加一个名字为<tt>key</tt>的列，并在这一列中添加类型为<tt>T</tt>（必须是<tt>int</tt>,
 * <tt>unsigned int</tt>, <tt>double</tt>,  <tt>std::string</tt>)
 * 的给定值。
 * 表格完成后，有不同的输出可能性，例如，用write_tex()输出到latex文件，或用write_text()输出为文本。
 * 两个（或多个）列可以通过两次（或多次）调用add_column_to_supercolumn()合并成一个
 * "超级列"，见那里。此外，还有一个函数可以为每一列设置数字输出的精度，还有几个函数可以在tex模式下规定列的格式和标题。
 * 在 step-13 教程程序中也给出了对这个类的详细解释。
 *
 *  <h3>Example</h3>
 * 这是一个简单的例子，演示了这个类的用法。第一列包括数字
 * $i=1 \dots n$ ，第二列 $1^2 \dots n^2$ ，第三列
 * $\sqrt{1}\dots\sqrt{n}$
 * ，其中第二列和第三列被合并成一个超级列，超级键为<tt>平方和根</tt>。此外，第一列向右对齐（默认是<tt>居中</tt>），平方根的精度被设置为6（而不是默认的4）。
 *
 *
 * @code
 * TableHandler table;
 * for (unsigned int i = 1; i <= n; ++i)
 * {
 *   table.add_value("numbers", i);
 *   table.add_value("squares", i i);
 *   table.add_value("square roots", std::sqrt(i));
 * }
 * // merge the second and third column
 * table.add_column_to_supercolumn("squares", "squares and roots");
 * table.add_column_to_supercolumn("square roots", "squares and roots");
 *
 * // additional settings
 * table.set_tex_format("numbers", "r");
 * table.set_precision("square roots", 6);
 *
 * // output
 * std::ofstream out_file("number_table.tex");
 * table.write_tex(out_file);
 * out_file.close();
 * @endcode
 *
 *  <h3>Dealing with sparse data: auto-fill mode</h3>
 * 当生成输出时，TableHandler期望所有的列都有完全相同数量的元素，因此结果实际上是一个表。这假定在每个迭代（时间步骤、非线性迭代等）中，你都会填满每一列。另一方面，这不一定是你想做的。例如，可能计算非线性残差的函数每隔几个时间步数才被调用一次；或者，计算网格统计的函数只有在网格实际上被细化时才被调用。在这些情况下，对于某些列来说，add_value()函数的调用频率会比较低，因此该列的元素会比较少；此外，这些元素不会与包含本次迭代中产生的其他数据元素的行对齐。一个完全不同的情况是，表被填满，在以后的时间里，我们使用其中的数据来计算其他行的元素；ConvergenceTable类做了这样的事情。
 * 为了支持这两种情况，TableHandler类有一个叫做<i> auto-fill
 * mode</i>的属性。默认情况下，自动填充模式是关闭的，但它可以通过调用set_auto_fill_mode()来启用。如果自动填充模式被启用，我们使用以下算法。
 *
 *
 *
 * - 当调用  <code>add_value(key, value)</code>  时，我们计算对应于  <code>key</code>  的列中的元素数量。让我们把这个数字称为  $m$  。
 *
 *
 *
 * - 我们还确定其他列中元素的最大数量；称其为  $n$  。
 *
 *
 *
 * - 如果 $m < n-1$ ，那么我们将 $n-m-1$ 对象的副本 <code>T()</code> 添加到该列。这里， <code>T</code> 是给定的数据类型 <code>value</code>. For example, if <code>T</code> 是数字类型，那么 <code>T()</code> is the number zero; if <code>T</code> 就是 <code>std::string</code>, then <code>T()</code> 是空字符串 <code>""</code>  。
 *
 *
 *
 * - 将给定的值添加到该列中。
 * 用默认元素填充该列，确保添加后该列的条目数与其他最长的列一样多。换句话说，如果我们跳过了之前对给定键的add_value()的调用，那么填充将把默认值输入到这一列。
 * 如果你试图跳过为一个键添加值，如果为这个键添加元素是你在给定的迭代或时间步长中要做的第一件事，那么所描述的算法将失败，因为我们将填充到<i>previous</i>迭代或时间步长中最长的列的长度。你可能需要在你的程序中重新安排添加到这一列的顺序，在添加到将一直被添加到的一列之后；或者，你可能想在每次迭代开始时，在表格中添加迭代的编号，例如在第1列。
 * 在上面的案例中，我们总是在列<b>above</b>中填充被添加到一列的元素。然而，也有一种情况，我们必须要垫上<b>below</b>。也就是说，如果前面的一行已经用
 * TableHandler::add_value(),
 * 完全填充，后续的几行已经被部分填充，然后我们通过write_text()或write_tex()要求输出。在这种情况下，最后几行只填充了部分内容，需要在最后加入的元素下面进行填充。和以前一样，我们通过使用与该列最后一个元素相同类型的默认构造对象来实现。
 *
 *
 * @ingroup textoutput
 *
 *
 */
class TableHandler
{
public:
  /**
   * 一组选项，当用write_text()函数输出时，应该如何格式化一个表格。存在以下可能性。
   *
   *
   *
   *
   *
   *
   * -  <code>table_with_headers</code>  : 表格的格式是，内容在每一列的键下对齐，也就是说，键位于每一列的顶部。这适合于列数不多的表格，在这种情况下，整个表格可以显示在屏幕上。输出看起来像这样。
   * @code
   *   key1 key2 key3
   *   0    0    ""
   *   1    0    ""
   *   2    13   a
   *   1    0    ""
   * @endcode
   *
   *
   *
   *
   *
   *
   * -  <code>table_with_separate_column_description</code>  : 当有很多列，而表格整体不能显示在屏幕上时，这是一种更好的格式。在这里，列键首先被逐一列在自己的行上，并被编号以提高可读性。此外，这些描述行的前缀是'#'，以标明这些行是注释，供那些想把下面的表作为数据来读的程序使用，应该忽略这些描述行。GNUPLOT就是这样一个程序，它将自动忽略这些前缀的行。使用这个选项的输出结果看起来像这样。
   * @code
   *   # 1: key1
   *   # 2: key2
   *   # 3: key3
   *   0 0  ""
   *   1 0  ""
   *   2 13 a
   *   1 0  ""
   * @endcode
   *
   *
   *
   *
   *
   *
   * -  <code>simple_table_with_separate_column_description</code>  : 这种格式与 <code>table_with_separate_column_description</code> 非常相似，但它跳过了用额外的空白空间对齐列。这提高了大表的write_text()的性能。输出示例。
   * @code
   *   # 1: key1
   *   # 2: key2
   *   # 3: key3
   *   0 0 ""
   *   1 0 ""
   *   2 13 a
   *   1 0 ""
   * @endcode
   *
   *
   *
   *
   *
   * -  <code>org_mode_table</code>  : 输出为 org-mode (http://orgmode.org/) 表格格式。将org-mode表格转换为HTML/LaTeX/csv很容易。  输出示例。
   * @code
   * | key1 | key2 | key3 |
   * | 0    | 0    | ""   |
   * | 1    | 0    | ""   |
   * | 2    | 13   | a    |
   * | 1    | 0    | ""   |
   * @endcode
   *
   *
   */
  enum TextOutputFormat
  {
    /**
     * 打印带有页眉的表格。
     *
     */
    table_with_headers,
    /**
     * 为每一列标签单独打印表格。
     *
     */
    table_with_separate_column_description,
    /**
     * 像table_with_separate_column_description一样，但不对包含列标签的列进行对齐。
     *
     */
    simple_table_with_separate_column_description,
    /**
     * 以org模式的格式打印表格。
     *
     */
    org_mode_table
  };

  /**
   * 构造函数。
   *
   */
  TableHandler();


  /**
   * 通过给它一个名字来声明表格中的一个列的存在。
   * 正如在该类的文档中所讨论的，这通常是不必要的。
   *
   * - 只是通过add_value()函数为一个给定的列键添加一个值，也声明了该列。因此，这个函数只有在以下情况下才是必要的：你希望一个列也能显示出来，即使你从未在这一列的任何行中添加过条目；或者，更有可能的是，如果你想规定以后打印列的顺序，在条目被放进列之前以特定顺序声明列。    后者的目的也可以通过在程序中以任何顺序向表中添加条目，然后调用set_column_order()来实现。然而，这种方法需要知道
   *
   * - 在你的软件的一个中心位置
   *
   * - 软件的其他部分写入的所有列键，以及它们应该如何排序。这对小程序来说很容易做到，但对大型代码库来说可能不可行，因为代码库的部分内容只根据运行时的参数来执行。)
   *
   */
  void
  declare_column(const std::string &key);

  /**
   * 添加一个键为<tt>key</tt>的列（如果尚未存在），并将类型为<tt>T</tt>的值添加到该列。<tt>T</tt>类型的值必须可以转换为<code>int,
   * unsigned int, double,  std::uint64_t,   std::string</code>
   * 中的一种，否则将导致编译器错误。
   *
   */
  template <typename T>
  void
  add_value(const std::string &key, const T value);

  /**
   * 如果一行只被部分填满，那么就把该行中某一列中不存在元素的所有元素设置为空字符串。这类似于在介绍中描述的
   * "auto_fill_mode"，但更普遍，因为它允许你开始向新行的某一列中写入内容，而不必知道该列在前一行已经被写入。
   * 如果当前行的所有列都被写入了，那么这个函数根本就不做任何事情。换句话说，从概念上讲，这个函数
   * "完成 "了当前行，尽管它的用例是 "开始 "一个新行。
   *
   */
  void
  start_new_row();

  /**
   * 开启或关闭自动填充模式。关于自动填充模式的描述，请参见该类的一般文档。
   *
   */
  void
  set_auto_fill_mode(const bool state);

  /**
   * 创建一个超级列（如果还不存在的话）并将列纳入其中。
   * 列和超级列的键分别是<tt>key</tt>和<tt>superkey</tt>。
   * 要将两个列<tt>c1</tt>和<tt>c2</tt>合并到一个超级列<tt>sc</tt>，因此调用<tt>add_column_to_supercolumn(c1,sc)</tt>和<tt>add_column_to_supercolumn(c2,sc)</tt>。
   * 关于列的顺序，超级列取代了被添加到超级列的第一列。在超级列中，输出的顺序是按照添加到超级列中的列的顺序。
   *
   */
  void
  add_column_to_supercolumn(const std::string &key,
                            const std::string &superkey);

  /**
   * 改变表格中列和超级列的顺序。
   * <tt>new_order</tt>包括列和超列的键和超键，其顺序是用户希望它们被输出。如果包含了超键，子列的键就不需要在这个向量中明确提及。
   * 超级列内的子列的顺序是不可改变的，仍然是按照列被添加到超级列的顺序。
   * 这个函数也可以用来把列数太多的大表分解成小表。例如，你可以用前五列调用这个函数，然后调用其中一个<tt>write_*</tt>函数，然后用下五列调用这个函数，再一次<tt>write_*</tt>，如此循环。
   *
   */
  void
  set_column_order(const std::vector<std::string> &new_order);

  /**
   * 设置<tt>precision</tt>，例如，双数或浮点数变量的写入方式。<tt>precision</tt>与调用<tt>out<<setprecision(precision)</tt>时相同。
   *
   */
  void
  set_precision(const std::string &key, const unsigned int precision);

  /**
   * 设置<tt>scientific_flag</tt>。真表示科学，假表示定点记数法。
   *
   */
  void
  set_scientific(const std::string &key, const bool scientific);

  /**
   * 设置tex输出的<tt>key</tt>列的标题。你可能想选择与<tt>key</tt>不同的，如果它包含公式或类似结构。
   *
   */
  void
  set_tex_caption(const std::string &key, const std::string &tex_caption);

  /**
   * 设置整个<tt>table</tt>的tex标题，用于tex输出。
   *
   */
  void
  set_tex_table_caption(const std::string &table_caption);

  /**
   * 设置这个<tt>表的标签</tt>，用于文本输出。
   *
   */
  void
  set_tex_table_label(const std::string &table_label);

  /**
   * 设置超级列的标题<tt>superkey</tt>，用于特克斯输出。
   * 你可能想选择与<tt>superkey</tt>不同的，如果它包含公式或类似的结构。
   *
   */
  void
  set_tex_supercaption(const std::string &superkey,
                       const std::string &tex_supercaption);

  /**
   * 设置一个列的tex输出格式，例如：<tt>c</tt>, <tt>r</tt>,
   * <tt>l</tt>,
   * 或者<tt>p{3cm}</tt>。默认是<tt>c</tt>。另外，如果没有为某一列调用这个函数，默认值将被预设为<tt>c</tt>。
   *
   */
  void
  set_tex_format(const std::string &key, const std::string &format = "c");

  /**
   * 将表格作为格式化的文本写入给定的流中。文本的格式是这样的：它将数据表示为格式化的文本列。为了避免在自动读取这些表格时出现问题，例如用于后处理，如果该表的某个单元格中的条目是空的（即它是通过调用add_value()函数用空字符串创建的），那么该表的条目将被打印为
   * <code>""</code>  。
   * 第二个参数表示如何显示列键。更多信息请参见TextOutputFormat的描述。
   *
   */
  void
  write_text(std::ostream &         out,
             const TextOutputFormat format = table_with_headers) const;

  /**
   * 将表写成一个tex文件。如果 @p with_header
   * 被设置为false，那么就不使用 <code>\\documentclass{...}</code>,
   * <code>\\begin{document}</code> 和 <code>\\end{document}</code>
   * 。这样，该文件就可以用 <code>\\input{table_file}</code>
   * 这样的命令包含到一个现有的tex文件中。
   *
   */
  void
  write_tex(std::ostream &file, const bool with_header = true) const;

  /**
   * 清除表格的行，即在所有底层存储数据结构上调用clear()。
   *
   */
  void
  clear();

  /**
   * 删除所有在当前行添加的值。例如，当一个时间步长被拒绝，所有关于它的数据记录都需要被丢弃时，这就很有用。
   *
   */
  void
  clear_current_row();

  /**
   * 为了使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)进行序列化，将此对象的数据读入或写入一个流中。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * @addtogroup  Exceptions  @{
   *
   */

  /**
   * 异常情况
   *
   */
  DeclException1(ExcColumnNotExistent,
                 std::string,
                 << "Column <" << arg1 << "> does not exist.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcSuperColumnNotExistent,
                 std::string,
                 << "Supercolumn <" << arg1 << "> does not exist.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcColumnOrSuperColumnNotExistent,
                 std::string,
                 << "Column or supercolumn <" << arg1 << "> does not exist.");

  /**
   * 异常情况
   *
   */
  DeclException4(ExcWrongNumberOfDataEntries,
                 std::string,
                 int,
                 std::string,
                 int,
                 << "Column <" << arg1 << "> has " << arg2
                 << " rows, but Column <" << arg3 << "> has " << arg4
                 << " rows.");

  /**
   * 异常情况
   *
   */
  DeclException1(ExcUndefinedTexFormat,
                 std::string,
                 << "<" << arg1 << "> is not a tex column format. Use "
                 << "'l', 'c', or 'r' to indicate left, centered, or "
                 << "right aligned text.");
  //@}
protected:
  /**
   * 封装了描述一个表的一列所需的所有数据的结构。
   *
   */
  struct Column
  {
    /**
     * <tt>std::map</tt>. 所需的构造器。
     *
     */
    Column();

    /**
     * 构造函数。
     *
     */
    Column(const std::string &tex_caption);

    /**
     * 用默认构建的元素填充这一列，达到参数给定的行数。
     *
     */
    void
    pad_column_below(const unsigned int length);

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写到一个流中，以便进行序列化。
     *
     */
    template <class Archive>
    void
    save(Archive &ar, const unsigned int version) const;

    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从一个流中读取此对象的数据，以达到序列化的目的。
     *
     */
    template <class Archive>
    void
    load(Archive &ar, const unsigned int version);

#ifdef DOXYGEN
    /**
     * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)从流中写入和读取此对象的数据，以达到序列化的目的。
     *
     */
    template <class Archive>
    void
    serialize(Archive &archive, const unsigned int version);
#else
    // This macro defines the serialize() method that is compatible with
    // the templated save() and load() method that have been implemented.
    BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif


    /**
     * 使所有条目的字符串缓存无效并重新计算最大长度max_length。
     *
     */
    void
    invalidate_cache();

    /**
     * 本栏内的条目列表。值总是立即转换为字符串，以提供统一的查询方法。
     *
     */
    std::vector<internal::TableEntry> entries;

    /**
     * 在tex输出中该列的标题。 默认情况下，这是由
     * <tt>TableHandler::add_value(...)</tt>.
     * 给<tt>TableHandler</tt>的关键字符串，这可以通过调用
     * <tt>TableHandler::set_tex_caption(...)</tt>. 来改变。
     *
     */
    std::string tex_caption;

    /**
     * tex输出中的列格式。
     * 默认情况下，这是<tt>"c"</tt>，意味着
     * "居中"。这可以通过调用
     * <tt>TableHandler::set_tex_format(...)</tt> 来改变，用<tt>"c", "r",
     * "l"<tt>表示居中，右或左。
     *
     */

    std::string tex_format;

    /**
     * 双数或浮点数条目以该精度写入（由用户设置）。
     * 默认为4。
     *
     */
    unsigned int precision;

    /**
     * <tt>scientific</tt>=false表示固定点符号。
     *
     */
    bool scientific;

    /**
     * 可以被派生类用于任意用途的标志。
     * 特别是，ConvergenceTable类使用该标志来表示已经计算过收敛信息的列，或者根本不应该计算。
     *
     */
    unsigned int flag;

    /**
     * 这个条目缓存了这个表中所有条目的最大长度，以字符为单位。
     *
     */
    unsigned int max_length;
  };

  /**
   * 帮助函数，给出<tt>column_order</tt>中提到的所有列的键的向量，其中每个超级列的键被其子列的键所替代。
   * 这个函数隐含地检查数据的一致性。结果在<tt>sel_columns</tt>中返回。
   *
   */
  void
  get_selected_columns(std::vector<std::string> &sel_columns) const;

  /**
   * 内置函数，给出表中的行数，并检查每一列的行数是否相等。这个函数可以在写输出之前调用。
   *
   */
  unsigned int
  n_rows() const;

  /**
   * 一个变量，按照用户希望的顺序存储列和超列的键。默认情况下这是添加列的顺序。这个顺序可以通过set_column_order()来改变。
   *
   */
  std::vector<std::string> column_order;

  /**
   * 一个从列键到列的映射（不是超级列）。
   * 这个字段被声明为可变的，以便write_text()和write_tex()函数可以是常数，即使它们在'auto_fill_mode'开启的情况下可能会在下面填充列。
   *
   */
  mutable std::map<std::string, Column> columns;

  /**
   * 一个从每个超级列的键到其子列的键的正确顺序的映射。
   * 允许一个超级列与一个列有相同的键。
   * 注意，我们在这里不使用<tt>multimap</tt>，因为每个超级列键的顺序是相关的。
   *
   */
  std::map<std::string, std::vector<std::string>> supercolumns;

  /**
   * 一个从超级列键到超级列标题的映射，在tex输出中使用。
   * 默认情况下，这些只是超级列的键，但它们可以通过<tt>set_tex_supercaptions(...)</tt>改变。
   *
   */
  std::map<std::string, std::string> tex_supercaptions;

  /**
   * 表格本身的标题。
   *
   */
  std::string tex_table_caption;
  /**
   * 该表的标签。
   *
   */
  std::string tex_table_label;

  /**
   * 表示是否应该使用自动填充模式的标志。
   *
   */
  bool auto_fill_mode;
};


namespace internal
{
  template <typename T>
  TableEntry::TableEntry(const T &t)
    : value(t)
  {}


  template <typename T>
  T
  TableEntry::get() const
  {
    // we don't quite know the data type in 'value', but
    // it must be one of the ones in the type list of the
    // boost::variant. so if T is not in the list, or if
    // the data stored in the TableEntry is not of type
    // T, then we will get an exception that we can
    // catch and produce an error message
    try
      {
        return boost::get<T>(value);
      }
    catch (...)
      {
        Assert(false,
               ExcMessage(
                 "This TableEntry object does not store a datum of type T"));
        throw;
      }
  }



  template <class Archive>
  void
  TableEntry::save(Archive &ar, const unsigned int) const
  {
    // write first an identifier for the kind
    // of data stored and then the actual
    // data, in its correct data type
    if (const int *p = boost::get<int>(&value))
      {
        char c = 'i';
        ar &c &*p;
      }
    else if (const unsigned int *p = boost::get<unsigned int>(&value))
      {
        char c = 'u';
        ar &c &*p;
      }
    else if (const double *p = boost::get<double>(&value))
      {
        char c = 'd';
        ar &c &*p;
      }
    else if (const std::string *p = boost::get<std::string>(&value))
      {
        char c = 's';
        ar &c &*p;
      }
    else if (const std::uint64_t *p = boost::get<std::uint64_t>(&value))
      {
        char c = 'l';
        ar &c &*p;
      }
    else
      Assert(false, ExcInternalError());
  }



  template <class Archive>
  void
  TableEntry::load(Archive &ar, const unsigned int)
  {
    // following what we do in the save()
    // function, first read in the data type
    // as a one-character id, and then read
    // the data
    char c;
    ar & c;

    switch (c)
      {
        case 'i':
          {
            int val;
            ar &val;
            value = val;
            break;
          }

        case 'u':
          {
            unsigned int val;
            ar &         val;
            value = val;
            break;
          }

        case 'd':
          {
            double val;
            ar &   val;
            value = val;
            break;
          }

        case 's':
          {
            std::string val;
            ar &        val;
            value = val;
            break;
          }

        case 'l':
          {
            std::uint64_t val;
            ar &          val;
            value = val;
            break;
          }

        default:
          Assert(false, ExcInternalError());
      }
  }
} // namespace internal



template <typename T>
void
TableHandler::add_value(const std::string &key, const T value)
{
  // see if the column already exists
  if (columns.find(key) == columns.end())
    declare_column(key);

  if (auto_fill_mode == true)
    {
      // follow the algorithm given in the introduction to this class
      // of padding columns as necessary
      unsigned int max_col_length = 0;
      for (const auto &column : columns)
        max_col_length =
          std::max(max_col_length,
                   static_cast<unsigned int>(column.second.entries.size()));

      while (columns[key].entries.size() + 1 < max_col_length)
        {
          columns[key].entries.push_back(internal::TableEntry(T()));
          internal::TableEntry &entry = columns[key].entries.back();
          entry.cache_string(columns[key].scientific, columns[key].precision);
          columns[key].max_length =
            std::max(columns[key].max_length,
                     static_cast<unsigned int>(
                       entry.get_cached_string().length()));
        }
    }

  // now push the value given to this function
  columns[key].entries.push_back(internal::TableEntry(value));
  internal::TableEntry &entry = columns[key].entries.back();
  entry.cache_string(columns[key].scientific, columns[key].precision);
  columns[key].max_length =
    std::max(columns[key].max_length,
             static_cast<unsigned int>(entry.get_cached_string().length()));
}



template <class Archive>
void
TableHandler::Column::save(Archive &ar, const unsigned int  /*version*/ ) const
{
  ar &entries &tex_caption &tex_format &precision &scientific &flag &max_length;
}



template <class Archive>
void
TableHandler::Column::load(Archive &ar, const unsigned int  /*version*/ )
{
  ar &entries &tex_caption &tex_format &precision &scientific &flag &max_length;
  invalidate_cache();
}


template <class Archive>
void
TableHandler::serialize(Archive &ar, const unsigned int)
{
  ar &column_order &columns &supercolumns &tex_supercaptions &tex_table_caption
    &tex_table_label &auto_fill_mode;
}


DEAL_II_NAMESPACE_CLOSE

#endif


