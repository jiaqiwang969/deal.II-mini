//include/deal.II-translator/base/utilities_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2021 by the deal.II authors
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

#ifndef dealii_utilities_h
#define dealii_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <functional>
#include <string>
#include <tuple>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

#ifdef DEAL_II_WITH_TRILINOS
#  include <Epetra_Comm.h>
#  include <Epetra_Map.h>
#  include <Teuchos_Comm.hpp>
#  include <Teuchos_RCP.hpp>
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  else
#    include <Epetra_SerialComm.h>
#  endif
#endif

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/core/demangle.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/vector.hpp>

#ifdef DEAL_II_WITH_ZLIB
#  include <boost/iostreams/filter/gzip.hpp>
#endif

DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

DEAL_II_NAMESPACE_OPEN

// forward declare Point
#ifndef DOXYGEN
template <int dim, typename Number>
class Point;
#endif

/**
 * 一个用于实用功能的命名空间，这些功能不是特别针对有限元计算或数值程序的，但在编写应用程序时，在各种情况下都需要。
 *
 *
 * @ingroup utilities
 *
 *
 */
namespace Utilities
{
  /**
   * 返回一个形式为 "deal.II version x.y.z "的字符串，其中
   * "x.y.z
   * "标识了你正在使用的deal.II的版本。这个信息也是由DEAL_II_PACKAGE_NAME和DEAL_II_PACKAGE_VERSION预处理器变量提供的。
   *
   */
  std::string
  dealii_version_string();

  /**
   * 使用希尔伯特空间填充曲线给 @p points
   * 中的每个点分配一个索引。  为此，将确定 @p points
   * 的一个边界框，在此基础上计算它们的整数坐标。
   * 线性索引是以比特的昏暗集合的形式给出的，从高到低。
   * 这样做是为了保持沿每个轴的位深度的最大分辨率。请注意，这个dim-integer索引仍然可以很容易地用于分类和排序，例如使用整数组的lexicographic排序。
   * 希尔伯特曲线的深度（即每维的比特数）默认等于
   * <code>64</code>  .
   * @note
   * 这个函数也可以处理退化的情况，例如，当包围盒在其中一个维度的大小为零时。
   *
   */
  template <int dim, typename Number>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<Point<dim, Number>> &points,
    const int                              bits_per_dim = 64);

  /**
   * 与上述相同，但对于整数坐标的点。
   *
   */
  template <int dim>
  std::vector<std::array<std::uint64_t, dim>>
  inverse_Hilbert_space_filling_curve(
    const std::vector<std::array<std::uint64_t, dim>> &points,
    const int                                          bits_per_dim = 64);

  /**
   * 将 @p index 的每个元素中最不重要的 @p bits_per_dim
   * 位（从最后一位开始）打包成一个无符号整数。 @p index
   * 的最后一个元素将用于设置所得整数中的第一个 @p
   * bits_per_dim 位，倒数第二个元素用于设置下一个 @p
   * bits_per_dim
   * 位，等等。为了将所有的数据装入输出，下面应该保持
   * <code>bits\_per\_dim dim <= 64</code>  。
   * 该函数在调试和可视化由inverse_Hilbert_space_filling_curve()返回的指数时很有用。
   * @note
   * 为了比较inverse_Hilbert_space_filling_curve()返回的指数，没有必要使用这个函数，因为这可以通过
   * <code>std::lexicographical_compare()</code>  轻松完成。
   *
   */
  template <int dim>
  std::uint64_t
  pack_integers(const std::array<std::uint64_t, dim> &index,
                const int                             bits_per_dim);

  /**
   * 如果库被配置为ZLIB，那么这个函数会压缩输入字符串，并返回一个包含压缩后的输入的非零终止的字符串。
   * 如果库没有配置启用ZLIB，则返回的字符串与输入字符串相同。
   * @param[in]  input 要压缩的字符串  @return
   * 输入字符串的压缩版本
   *
   */
  std::string
  compress(const std::string &input);

  /**
   * 如果库中配置了ZLIB，那么这个函数假定输入的字符串已经用compress()函数进行了压缩，并返回原始的解压缩字符串。
   * 如果库没有配置启用ZLIB，则返回的字符串与输入字符串相同。
   * @param[in]  compressed_input
   * 一个压缩的字符串，由函数compress()返回  @return
   * 原始未压缩的字符串。
   *
   */
  std::string
  decompress(const std::string &compressed_input);

  /**
   * 将二进制输入编码为base64字符串。
   * Base64是一组二进制到文本的编码方案，通过将二进制数据转化为radix-64表示法，以ASCII字符串格式表示。Base64被设计用来在只可靠地支持文本内容的通道上携带以二进制格式存储的数据。它也被用来以独立于机器的方式存储二进制格式。
   * @param  binary_input
   * 一个字符向量，代表你的输入为二进制数据。    @return
   * 一个包含二进制输入的字符串，作为base64字符串。
   *
   */
  std::string
  encode_base64(const std::vector<unsigned char> &binary_input);

  /**
   * 将base64字符串解码为二进制输出。
   * 这是上面的encode_base64()函数的逆运算。      @param
   * base64_input 一个包含base64格式输入的字符串。    @return
   * 一个字符向量，表示你的输入为二进制数据。
   *
   */
  std::vector<unsigned char>
  decode_base64(const std::string &base64_input);

  /**
   * 将一个数字 @p value
   * 转换成一个字符串，有多少个数字就用多少个前导零填充。
   * 如果第二个参数保留为默认值，则数字不会被填充前导零。其结果与调用标准C++
   * `std::to_string` （或旧的C函数`itoa()`）的结果相同。
   * 这个函数需要一个`无符号int'作为参数。因此，如果你用`signed
   * int`（当然与`int`的类型相同）调用它，参数会被隐含地转换为无符号整数，负数可能不会像你希望的那样被打印出来。同样地，如果你用`long
   * int`调用函数，打印出来的结果可能会显示转换为`unsigned
   * int`时的溢出效果。
   * @note  不鼓励使用这个函数，用户应该使用
   * <code>Utilities::to_string()</code>
   * 代替。在目前的实现中，该函数只是简单地调用<code>to_string
   * @<unsigned  int  @>()</code>.  。
   *
   */
  std::string
  int_to_string(const unsigned int value,
                const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * 将一个数字 @p value 转换成一个字符串，其中有 @p digits
   * 个字符。字符串被填充了前导零，在可能的减号之后。
   * 因此，填充零的总数是 @p digits
   * 减去任何符号、小数点和 @p value. 的数字。
   * 如果第二个参数保持其默认值，数字就不会被填充前导零。结果与调用C++函数
   * `std::to_string()` （对于积分类型）或调用
   * `boost::lexical_cast()` （对于所有其他类型）的情况相同。
   *
   */
  template <typename number>
  std::string
  to_string(const number       value,
            const unsigned int digits = numbers::invalid_unsigned_int);

  /**
   * 确定需要多少个数字来表示最多和给定数字一样大的数字。
   *
   */
  unsigned int
  needed_digits(const unsigned int max_number);

  /**
   * 这个函数允许在 @p n_digits 的精度之后切断一个浮点数 @p
   * number ，即在科学浮点符号的小数点后 @p n_digits
   * 位之后切断。当解释为四舍五入操作时，这个函数减少了浮点数的绝对值，并且总是向零舍入，因为小数点位被简单地切断了。
   *
   */
  template <typename Number>
  Number
  truncate_to_n_digits(const Number number, const unsigned int n_digits);

  /**
   * 给定一个字符串，将其转换为一个整数。如果不可能，则抛出一个断言。
   *
   */
  int
  string_to_int(const std::string &s);

  /**
   * 返回一个描述对象尺寸的字符串。通常，deal.II库以及用户代码中的函数需要定义一个字符串，包含一些使用两个模板参数定义的对象的模板尺寸：dim（对象的拓扑尺寸）和spacedim（嵌入欧几里得空间的尺寸）。
   * 由于在所有deal.II类中，默认情况下spacedim等于dimension，上述字符串通常被缩减为"<dim>"，而不是"<dim,spacedim>"。
   * 如果dim等于spacedim，该函数返回一个包含 "dim
   * "的字符串，否则返回 "dim,spacedim"。
   *
   */
  std::string
  dim_string(const int dim, const int spacedim);

  /**
   * 给出一个字符串列表，将其转换为一个整数列表。如果不可能，则抛出一个断言。
   *
   */
  std::vector<int>
  string_to_int(const std::vector<std::string> &s);

  /**
   * 给定一个字符串，将其转换为一个双数。如果不可能，抛出一个断言。
   *
   */
  double
  string_to_double(const std::string &s);


  /**
   * 给出一个字符串的列表，将其转换为一个双数的列表。如果不可能，抛出一个断言。
   *
   */
  std::vector<double>
  string_to_double(const std::vector<std::string> &s);


  /**
   * 给出一个包含由 @p delimiter,
   * 分隔的文本的字符串，将其分割成其组成部分；对于每个组成部分，去除前面和后面的空格。分隔符的默认值是逗号，因此该函数可以分割逗号分隔的字符串列表。
   * 为了使来自表格的数据输入更简单，如果输入的字符串以定界符结束（后面可能有任意数量的空白），那么最后这个定界符将被忽略。比如说。
   * @code
   * Utilities::split_string_list("abc; def; ghi; ", ';');
   * @endcode
   * 产生的3个元素的输出 <code>{"abc","def","ghi"}</code>
   * 与你在输入时得到的相同
   * @code
   * Utilities::split_string_list("abc; def; ghi", ';');
   * @endcode
   * 或
   * @code
   * Utilities::split_string_list("abc; def; ghi;", ';');
   * @endcode
   * 作为这一规则的结果，像这样的调用
   * @code
   * Utilities::split_string_list(" ; ", ';');
   * @endcode
   * 产生了一个单元素的列表。由于修剪了空白，这个单元素是空字符串。
   * 这个函数可以消化定界符是空格的情况。在这种情况下，它返回字符串中的所有单词。结合上面的规则，这意味着
   * @code
   * Utilities::split_string_list("abc def ghi ", ' ');
   * @endcode
   * 尽管在字符串的末尾存在空格，*还是会产生上面的输出
   * <code>{"abc","def","ghi"}</code> 的3元素列表。此外。
   * @code
   * Utilities::split_string_list("      ", ' ');
   * @endcode
   * 无论字符串中有多少空格，都会产生一个空列表。
   *
   */
  std::vector<std::string>
  split_string_list(const std::string &s, const std::string &delimiter = ",");


  /**
   * split_string_list()的特殊化，用于分隔符为单个字符的情况。
   *
   */
  std::vector<std::string>
  split_string_list(const std::string &s, const char delimiter);


  /**
   * 取一个文本，通常是一个文档或其他东西，并尝试将其分割成最多
   * @p width 个字符宽的独立文本行，在文本中由 @p delimiter
   * 标记的位置进行分割。如果这是不可能的，则返回长于
   * @p width.
   * 的最短的几行，分隔符的默认值是一个空格字符。如果original_text包含换行符(/n)，字符串也会在这些位置被分割。
   *
   */
  std::vector<std::string>
  break_text_into_lines(const std::string &original_text,
                        const unsigned int width,
                        const char         delimiter = ' ');

  /**
   * 如果给定的模式字符串出现在字符串的第一个位置，则返回true。
   *
   */
  bool
  match_at_string_start(const std::string &name, const std::string &pattern);

  /**
   * 读取一个（有符号的）整数，从第二个参数指示的 @p
   * name
   * 中的位置开始，并将这个整数与它在字符串中所占的字符数一起作为一对返回。
   * 如果在指定的位置没有读到整数，则返回
   * (-1,numbers::invalid_unsigned_int)  。
   *
   */
  std::pair<int, unsigned int>
  get_integer_at_position(const std::string &name, const unsigned int position);

  /**
   * 返回一个字符串，其中 @p input 中所有出现的 @p from 都被
   * @p to. 替换。
   *
   */
  std::string
  replace_in_string(const std::string &input,
                    const std::string &from,
                    const std::string &to);

  /**
   * 返回一个字符串，删除 @p input
   * 开头和结尾的所有标准空白字符（包括'<tt>\t</tt>',
   * '<tt>\n</tt>', 和 '<tt>\r</tt>'）。
   *
   */
  std::string
  trim(const std::string &input);

  /**
   * 从以 @p a 为中心、标准差为 @p sigma.
   * 的归一化高斯概率分布中生成一个随机数，每次调用该函数时返回的数字都会不同。
   * 这个函数是可重入的，也就是说，它可以安全地从多个线程同时调用。此外，每个线程每次都会得到相同的数字序列。另一方面，如果你通过线程积木运行
   * Threads::Task
   * 对象，那么任务将被分配到大部分随机的线程中，并且可能在程序的不同运行中获得不同的随机数序列，因为之前的任务可能已经消耗了为你所在的线程生成的前几个随机数。如果这是一个问题，你需要在每次想从一个定义的点开始时创建自己的随机数生成器对象。
   * @note
   * 与系统函数rand()一样，这个函数在每次程序启动时都会产生相同的随机数序列。这是调试代码的一个重要特性，但它使我们无法真正验证代码的统计特性。对于`rand()`，你可以调用`srand()`来
   * "播种
   * "随机数发生器，以便在每次调用程序时得到不同的随机数序列。然而，这个函数不允许给随机数发生器播种。如果你需要这个，如上所述，请使用C++或BOOST设施之一。
   *
   */
  double
  generate_normal_random_number(const double a, const double sigma);

  /**
   * 返回变量类型的字符串描述  @p t.
   * 一般来说，C++使用混杂的名称来识别类型。这个函数使用
   * boost::core::demangle
   * 来返回一个人类可读的字符串，描述作为参数传递的变量的类型。
   *
   */
  template <class T>
  std::string
  type_to_string(const T &t);

  /**
   * 计算一个固定的幂，作为模板参数提供，是一个数字的计算。
   * 这个函数提供了一种有效的方法来计算诸如
   * <code>t^N</code> where <code>N</code>
   * 是一个已知的数字在编译时。    使用这个函数，如
   * <code>fixed_power@<dim@> (n)</code>  .
   *
   */
  template <int N, typename T>
  T
  fixed_power(const T t);

  /**
   * 替换 <code>std::pow</code>
   * ，允许对常量表达式参数进行编译时计算。 @p base
   * 必须是整数类型，指数 @p iexp 不能是负数。
   *
   */
  template <typename T>
  constexpr T
  pow(const T base, const int iexp)
  {
#if defined(DEBUG) && !defined(DEAL_II_CXX14_CONSTEXPR_BUG)
    // Up to __builtin_expect this is the same code as in the 'Assert' macro.
    // The call to __builtin_expect turns out to be problematic.
    if (!(iexp >= 0))
      ::dealii::deal_II_exceptions::internals::issue_error_noreturn(
        ::dealii::deal_II_exceptions::internals::abort_or_throw_on_exception,
        __FILE__,
        __LINE__,
        __PRETTY_FUNCTION__,
        "iexp>=0",
        "ExcMessage(\"The exponent must not be negative!\")",
        ExcMessage("The exponent must not be negative!"));
#endif
    // The "exponentiation by squaring" algorithm used below has to be
    // compressed to one statement due to C++11's restrictions on constexpr
    // functions. A more descriptive version would be:
    //
    // <code>
    // if (iexp <= 0)
    //   return 1;
    //
    // // avoid overflow of one additional recursion with pow(base * base, 0)
    // if (iexp == 1)
    //   return base;
    //
    // // if the current exponent is not divisible by two,
    // // we need to account for that.
    // const unsigned int prefactor = (iexp % 2 == 1) ? base : 1;
    //
    // // a^b = (a*a)^(b/2)      for b even
    // // a^b = a*(a*a)^((b-1)/2 for b odd
    // return prefactor * dealii::Utilities::pow(base*base, iexp/2);
    // </code>

    static_assert(std::is_integral<T>::value, "Only integral types supported");

    return iexp <= 0 ?
             1 :
             (iexp == 1 ? base :
                          (((iexp % 2 == 1) ? base : 1) *
                           dealii::Utilities::pow(base * base, iexp / 2)));
  }

  /**
   * 对 <tt>std::lower_bound</tt>
   * 的优化替换，用于在列索引的范围内搜索。对于目前的应用来说，执行时间大约减少了一半，部分原因是二进制搜索被小循环长度的线性搜索所取代。
   * 这个功能的另一个原因相当不明显：当使用GCC的libstdc++函数
   * std::lower_bound, 时，复杂性为O(log(N))，符合要求。
   * 然而，当使用GCC
   * libstdc++的调试版本时，正如我们在运行测试套件时所做的那样，那么
   * std::lower_bound 测试序列是否实际上是相对于枢轴 "值
   * "而言的分区（也就是说，本质上序列是按照二进制搜索的要求进行排序的）。
   * 然而，验证这一点意味着 std::lower_bound
   * 的复杂度跃升至O(N)；我们在下面调用这个函数O(N)次，使得整体复杂度为O(N*2)。其后果是，一些有大网格的测试完全跑出了测试的时间限制，并在libstdc++调试模式下失败。这个函数只是假设序列是排序的，而我们根本不做额外的检查。
   *
   */
  template <typename Iterator, typename T>
  Iterator
  lower_bound(Iterator first, Iterator last, const T &val);


  /**
   * 与上面的函数相同，但接受一个参数，用来比较迭代器所指向的对象序列的各个元素。
   *
   */
  template <typename Iterator, typename T, typename Comp>
  Iterator
  lower_bound(Iterator first, Iterator last, const T &val, const Comp comp);

  /**
   * 给定一个排列向量（即一个向量 $p_0\ldots p_{N-1}$
   * ，其中每个 $p_i\in [0,N)$ 和 $p_i\neq p_j$ 为 $i\neq j$
   * ），产生反向排列 $q_i=N-1-p_i$  。
   *
   */
  template <typename Integer>
  std::vector<Integer>
  reverse_permutation(const std::vector<Integer> &permutation);

  /**
   * 给出一个互换向量（即一个向量 $p_0\ldots p_{N-1}$
   * ，其中每个 $p_i\in [0,N)$ 和 $p_i\neq p_j$ 为 $i\neq j$
   * ），产生反互换 $q_0\ldots q_{N-1}$ ，以便 $q_{p_i}=p_{q_i}=i$
   * .
   *
   */
  template <typename Integer>
  std::vector<Integer>
  invert_permutation(const std::vector<Integer> &permutation);

  /**
   * 给定一个T类型的任意对象，使用 boost::serialization
   * 工具将该对象打包成一个字符向量，并将其追加到给定的缓冲区中。被添加到缓冲区的元素的数量将被返回。该对象可以使用下面的
   * Utilities::unpack 函数进行解包。
   * 如果库在编译时启用了ZLIB，那么输出缓冲区可以被压缩。这可以通过参数
   * @p allow_compression,
   * 来触发，并且只有在启用ZLIB时才有效。
   * 如果考虑用同一个缓冲区进行多次连续调用，出于性能考虑，建议确保缓冲区有足够的容量。
   *
   */
  template <typename T>
  size_t
  pack(const T &          object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression = true);

  /**
   * 使用上面提到的pack函数，只为给定对象创建并返回一个缓冲区。
   * 如果该库在编译时启用了ZLIB，那么输出缓冲区可以被压缩。这可以通过参数
   * @p allow_compression,
   * 来触发，并且只有在启用ZLIB时才有效。
   *
   */
  template <typename T>
  std::vector<char>
  pack(const T &object, const bool allow_compression = true);

  /**
   * 给定一个字符向量，通过调用函数 Utilities::pack,
   * 获得，在一个T类型的对象中恢复其内容。
   * 这个函数使用 boost::serialization
   * 实用程序从字符向量中解压缩对象，它是函数
   * Utilities::pack(). 的逆函数。 @p allow_compression
   * 参数表示要读取的缓冲区是否以前用ZLIB压缩过，并且只有在启用ZLIB时才有效。
   * @note  因为这个函数的参数不取决于模板类型  @p T,
   * 你必须在调用这个函数时手动指定模板参数。
   * @note
   * 如果你想打包()或解压()对象的数组，那么下面的方法可行。
   * @code
   *  double array[3] = {1,2,3};
   *  std::vector<char> buffer = Utilities::pack(array);
   * @endcode
   * 然而，反过来就不行了。
   * @code
   *  array = Utilities::unpack<double[3]>(buffer);
   * @endcode
   * 这是因为C++不允许函数返回数组。
   * 因此，有一个单独的解包()函数用于数组，见下文。
   *
   */
  template <typename T>
  T
  unpack(const std::vector<char> &buffer, const bool allow_compression = true);

  /**
   * 与上面的解包函数相同，但在给定的 std::vector<char>
   * 类型的打包缓冲区上（一部分）采取恒定的迭代器来代替。
   * @p allow_compression
   * 参数表示要读取的缓冲区是否以前用ZLIB压缩过，并且只有在启用ZLIB时才有效。
   *
   */
  template <typename T>
  T
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         const bool                               allow_compression = true);

  /**
   * 给定一个字符向量，通过调用函数 Utilities::pack,
   * 获得，将其内容还原为T型数组。    这个函数使用
   * boost::serialization
   * 实用程序从字符向量中解压缩对象，它是函数
   * Utilities::pack(). 的逆函数。 @p allow_compression
   * 参数表示要读取的缓冲区是否以前用ZLIB压缩过，只有当ZLIB被启用时才有效。
   * @note
   * 这个函数的存在是由于C++的一个怪癖。具体来说，如果你想打包()或解压()对象的数组，那么下面的方法就可以。
   * @code
   *  double array[3] = {1,2,3};
   *  std::vector<char> buffer = Utilities::pack(array);
   * @endcode
   * 然而，反过来就不行了。
   * @code
   *  array = Utilities::unpack<double[3]>(buffer);
   * @endcode
   * 这是因为C++不允许函数返回数组。
   * 因此，当前的函数允许写
   * @code
   *  Utilities::unpack(buffer, array);
   * @endcode
   * 注意，与其他unpack()函数不同，不需要明确指定模板参数，因为它们可以从第二个参数中推导出来。
   *
   */
  template <typename T, int N>
  void
  unpack(const std::vector<char> &buffer,
         T (&unpacked_object)[N],
         const bool allow_compression = true);

  /**
   * 和上面的解包函数一样，但是在给定的 std::vector<char>
   * 类型的打包缓冲区上（一部分）采取常数迭代器来代替。
   * @p allow_compression
   * 参数表示要读取的缓冲区是否以前用ZLIB压缩过，并且只有在启用ZLIB时才有效。
   *
   */
  template <typename T, int N>
  void
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         T (&unpacked_object)[N],
         const bool allow_compression = true);

  /**
   * 检查 @p number 中 @p n 位置的位是否被设置。
   *
   */
  bool
  get_bit(const unsigned char number, const unsigned int n);


  /**
   * 将 @p number 中位置 @p n 的位设置为值 @p x. 。
   *
   */
  void
  set_bit(unsigned char &number, const unsigned int n, const bool x);


  /**
   * 将一个类型为 `std::unique_ptr<From>` 的对象转换为类型为
   * `std::unique_ptr<To>`,
   * 的对象，这里假定我们可以使用`dynamic_cast`将指向`From`的指针转换为指向`To`的指针。
   *
   *
   *
   *
   *
   *
   * - 换句话说，我们假设`From'和`To'是通过一个类的层次结构连接起来的，而且所指向的对象实际上是一个同时包含`From'和`To'的类型。一个例子是如果`To`是从`From`派生出来的，或者反过来。    如果`dynamic_cast`不成功，该函数会抛出一个类型为 `std::bad_cast` 的异常。这和你在对象类型（但不是指针类型）之间的常规`dynamic_cast`不成功时得到的异常一样。    这个函数如何工作的例子如下。
   * @code
   * // A base class. Assume that it has virtual
   * // functions so that dynamic_cast can work.
   * class B
   * {
   *   ...
   * };
   *
   * // A derived class
   * class D : public B
   * {
   *   ...
   * };
   *
   * // A factory function
   * std::unique_ptr<B> create_object (...)
   * {
   *   ...
   * }
   *
   * void foo (...)
   * {
   *   std::unique_ptr<B> b = create_object (...);
   *
   *   // Assume that we know for some reason that the object above must
   *   // have created a D object but returned it as a std::unique_ptr<B>.
   *   // In order to access the D functionality, we need to cast the
   *   // pointer. Use the equivalent to dynamic_cast:
   *   std::unique_ptr<D> d = dynamic_unique_cast<D>(std::move(b));
   *
   *   // If the object really was a D, then 'd' now points to it. Note
   *   // also that in accordance with the semantics of std::unique_ptr,
   *   // it was necessary to std::move the 'b' object, and indeed 'b'
   *   // now no longer points to anything
   *
   * -- ownership has been
   *   // transferred to 'd'!
   * @endcode
   * @note  这个函数不会尝试转换 `std::unique_ptr`
   * 对象所存储的`Deleter`对象。因此，该函数只在deleter对象处于默认状态时起作用，即，如果它们是
   * `std::default_delete<To>` 和 `std::default_delete<From>`.
   * 类型的对象。
   *
   */
  template <typename To, typename From>
  std::unique_ptr<To>
  dynamic_unique_cast(std::unique_ptr<From> &&p);

  /**
   * 返回基本值。默认：返回输入。
   *
   */
  template <typename T>
  T &
  get_underlying_value(T &p);

  /**
   * 返回基本值。对 std::shared_ptr<T>. 的特殊化。
   *
   */
  template <typename T>
  T &
  get_underlying_value(std::shared_ptr<T> &p);

  /**
   * 返回底层值。对const  std::shared_ptr<T>. 的特化。
   *
   */
  template <typename T>
  T &
  get_underlying_value(const std::shared_ptr<T> &p);

  /**
   * 返回基本值。对 std::unique_ptr<T>. 的特化
   *
   */
  template <typename T>
  T &
  get_underlying_value(std::unique_ptr<T> &p);

  /**
   * 返回基本值。对const  std::unique_ptr<T>. 的特化。
   *
   */
  template <typename T>
  T &
  get_underlying_value(const std::unique_ptr<T> &p);

  /**
   * 一个用于探测系统属性的实用函数的命名空间。
   * @ingroup utilities
   *
   */
  namespace System
  {
    /**
     * 返回由 "uptime
     * "返回的CPU负载。注意，这个数字的解释取决于机器中实际的处理器数量。目前只在Linux上实现，使用/proc/loadavg伪文件，在其他系统上，我们只是返回0。
     *
     */
    double
    get_cpu_load();

    /**
     * 以字符串形式返回vectorization.h中DEAL_II_VECTORIZATION_WIDTH_IN_BITS描述的矢量化的指令集扩展。可能的返回值列表是。
     * <table> <tr> <td><tt>VECTORIZATION_LEVEL</tt></td> <td>Return
     * Value</td> <td>Width in bits</td> </tr> <tr> <td>0</td>
     * <td>disabled</td> <td>64</td> </tr> <tr> <td>1</td>
     * <td>SSE2/AltiVec</td> <td>128</td> </tr> <tr> <td>2</td> <td>AVX</td>
     * <td>256</td> </tr> <tr> <td>3</td> <td>AVX512</td> <td>512</td> </tr>
     * </table>
     *
     */
    const std::string
    get_current_vectorization_level();

    /**
     * 保存内存使用信息的结构，单位是KB。由get_memory_stats()使用。详见man
     * 5 proc entry /status。
     *
     */
    struct MemoryStats
    {
      /**
       * 虚拟内存的峰值大小，单位为 kB。
       *
       */
      unsigned long int VmPeak;

      /**
       * 当前虚拟内存大小，单位为 kB。
       *
       */
      unsigned long int VmSize;

      /**
       * 峰值常驻内存大小，单位为千字节。也被称为
       * "高水位线"(HWM)。
       *
       */
      unsigned long int VmHWM;

      /**
       * 当前常驻内存大小，单位为 kB。也被称为
       * "常驻组大小"(RSS)。
       *
       */
      unsigned long int VmRSS;
    };


    /**
     * 在 @p stats
     * 结构中填充有关该进程的内存消耗信息。这仅在Linux上实现。
     *
     */
    void
    get_memory_stats(MemoryStats &stats);


    /**
     * 返回这个进程所运行的主机的名称。
     *
     */
    std::string
    get_hostname();


    /**
     * 返回现在的时间为HH:MM:SS。
     *
     */
    std::string
    get_time();

    /**
     * 返回当前日期为YYY/MM/DD。MM和DD可以是一个或两个数字。
     *
     */
    std::string
    get_date();

    /**
     * 调用系统函数posix_memalign，如果没有，则调用一个替代函数，以分配具有某种最小对齐方式的内存。然后第一个参数将返回这个内存块的指针，以后可以通过标准的
     * <code>free</code> 调用来释放。          @param  memptr
     * 一个指针变量的地址，在这个调用之后将指向所分配的内存。
     * @param  alignment 内存块的最小对齐方式，字节数。
     * @param  size 要分配的内存块的大小，以字节为单位。
     * @note
     * 这个函数在内部检查错误代码，而不是把这个任务留给调用站点。
     *
     */
    void
    posix_memalign(void **memptr, std::size_t alignment, std::size_t size);
  } // namespace System


#ifdef DEAL_II_WITH_TRILINOS
  /**
   * 这个命名空间提供了一些用于初始化Trilinos对象的基本结构（例如，矩阵、向量和预调节器）。
   *
   */
  namespace Trilinos
  {
    /**
     * 返回创建Epetra_Maps所需的Trilinos Epetra_Comm对象。
     * 如果deal.II被配置为使用不支持MPI的编译器，那么生成的通信器将是一个串行的。否则，该通信器将对应于MPI_COMM_WORLD，即一个包含该MPI宇宙中所有进程的通信器。
     *
     */
    const Epetra_Comm &
    comm_world();

    /**
     * 返回一个创建Epetra_Maps所需的Trilinos Epetra_Comm对象。
     * 如果deal.II被配置为使用不支持MPI的编译器，那么产生的通信器将是一个串行通信器。否则，该通信器将对应于MPI_COMM_SELF，即一个只包括这一个处理器的通信器。
     *
     */
    const Epetra_Comm &
    comm_self();

    /**
     * 返回创建 Tpetra::Maps. 所需的 Teuchos::Comm 对象
     * 如果deal.II被配置为使用不支持MPI的编译器，那么产生的通信器将是一个串行通信器。否则，该通信器将对应于MPI_COMM_SELF，即一个只包括这一个处理器的通信器。
     *
     */
    const Teuchos::RCP<const Teuchos::Comm<int>> &
    tpetra_comm_self();

    /**
     * 给定一个通讯器，复制它。如果给定的通信器是串行的，这意味着只是返回它的一个副本。另一方面，如果它是%并行的，我们就复制底层的MPI_Comm对象：我们创建一个单独的MPI通信器，包含相同的处理器和相同的顺序，但有一个独立的标识符，与给定的通信器不同。该函数返回一个指向从Epetra_Comm派生的类的新对象的指针。这个函数的调用者需要承担这个函数的所有权。返回的对象应使用destroy_communicator()函数销毁。
     * 这个设施是用来分离通信流的。例如，一个程序可以简单地使用MPI_Comm_World来处理一切。但是很容易出现这样的情况：有时不是所有的处理器都参与到一个旨在实现全局的通信中来
     *
     * - 例如，如果我们在一个粗略的网格上组装一个矩阵，其单元数少于处理器的数量，一些处理器可能不会与其他处理器同步他们的矩阵，因为他们没有写进矩阵，因为他们没有拥有单元。这显然是一个错误。然而，如果这些处理器只是继续他们的工作，而下一个%并行操作恰好是对不同矩阵的同步，那么同步就可能成功
     *
     * - 偶然的，因为不同的处理器在谈论不同的矩阵。        如果我们对不同的矩阵使用不同的通信器，就可以避免这种情况，这就减少了本应分开的通信因为发生在同一通信器上而不被识别的可能性。此外，可以想象的是，一些MPI操作可以使用多个线程进行并行化，因为它们的通信器识别了相关的通信，而不是像只使用单个通信器的顺序程序那样识别它们的相对时间。
     *
     */
    Epetra_Comm *
    duplicate_communicator(const Epetra_Comm &communicator);

    /**
     * 给出一个由Diplicate_communicator()函数创建的Epetra通信器，销毁底层MPI通信器对象，并将Epetra_Comm对象重置为comm_self()的结果。
     * 当不再需要diplicate_communicator()的结果时，有必要调用这个函数。原因是在该函数中，我们首先创建一个新的MPI_Comm对象，然后围绕它创建一个Epetra_Comm。虽然我们可以负责销毁后者，但它并没有销毁通讯器，因为它只能假设它可能还被程序中的其他对象使用。因此，我们必须自己明确地销毁它。
     * 这个函数正是这样做的。因为这必须在Epetra_Comm对象仍然存在的情况下进行，所以它首先重置后者，然后销毁通信器对象。
     * @note
     * 如果你在一个不是由diplicate_communicator()创建的Epetra_Comm对象上调用这个函数，你很可能做错了什么。请不要这样做。
     *
     */
    void
    destroy_communicator(Epetra_Comm &communicator);

    /**
     * 返回给定的 @ref GlossMPICommunicator "communicator "
     * 对象中存在的MPI进程的数量。如果这是一个顺序作业（即程序根本没有使用MPI，或者使用MPI但只启动了一个MPI进程），那么通信器必然只涉及一个进程，该函数返回1。
     *
     */
    unsigned int
    get_n_mpi_processes(const Epetra_Comm &mpi_communicator);

    /**
     * 返回给定通信器所描述的进程空间中的当前MPI进程的编号。对于每个进程来说，这将是一个唯一的值，介于零和（小于）所有进程的数量（由get_n_mpi_processes()给出）之间。
     *
     */
    unsigned int
    get_this_mpi_process(const Epetra_Comm &mpi_communicator);

    /**
     * 给定一个Trilinos
     * Epetra地图，创建一个新的地图，其元素与处理器的细分相同，但使用给定的通信器对象而不是存储在第一个参数中的对象。实质上，这意味着我们创建一个地图，以同样的方式在相同的处理器之间进行通信，但使用一个单独的通道。
     * 这个函数通常与一个通过Diplicate_communicator()函数获得的通信器一起使用。
     *
     */
    Epetra_Map
    duplicate_map(const Epetra_BlockMap &map, const Epetra_Comm &comm);
  } // namespace Trilinos

#endif


} // namespace Utilities


// --------------------- inline functions

namespace Utilities
{
  template <int N, typename T>
  inline T
  fixed_power(const T x)
  {
    Assert(
      !std::is_integral<T>::value || (N >= 0),
      ExcMessage(
        "The non-type template parameter N must be a non-negative integer for integral type T"));

    if (N == 0)
      return T(1.);
    else if (N < 0)
      return T(1.) / fixed_power<-N>(x);
    else
      // Use exponentiation by squaring:
      return ((N % 2 == 1) ? x * fixed_power<N / 2>(x * x) :
                             fixed_power<N / 2>(x * x));
  }



  template <class T>
  inline std::string
  type_to_string(const T &t)
  {
    return boost::core::demangle(typeid(t).name());
  }



  template <typename Iterator, typename T>
  inline Iterator
  lower_bound(Iterator first, Iterator last, const T &val)
  {
    return Utilities::lower_bound(first, last, val, std::less<T>());
  }



  template <typename Iterator, typename T, typename Comp>
  inline Iterator
  lower_bound(Iterator first, Iterator last, const T &val, const Comp comp)
  {
    // verify that the two iterators are properly ordered. since
    // we need operator- for the iterator type anyway, do the
    // test as follows, rather than via 'last >= first'
    Assert(last - first >= 0,
           ExcMessage(
             "The given iterators do not satisfy the proper ordering."));

    unsigned int len = static_cast<unsigned int>(last - first);

    if (len == 0)
      return first;

    while (true)
      {
        // if length equals 8 or less,
        // then do a rolled out
        // search. use a switch without
        // breaks for that and roll-out
        // the loop somehow
        if (len < 8)
          {
            switch (len)
              {
                case 7:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 6:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 5:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 4:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 3:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 2:
                  if (!comp(*first, val))
                    return first;
                  ++first;
                  DEAL_II_FALLTHROUGH;
                case 1:
                  if (!comp(*first, val))
                    return first;
                  return first + 1;
                default:
                  // indices seem
                  // to not be
                  // sorted
                  // correctly!? or
                  // did len
                  // become==0
                  // somehow? that
                  // shouldn't have
                  // happened
                  Assert(false, ExcInternalError());
              }
          }



        const unsigned int half   = len >> 1;
        const Iterator     middle = first + half;

        // if the value is larger than
        // that pointed to by the
        // middle pointer, then the
        // insertion point must be
        // right of it
        if (comp(*middle, val))
          {
            first = middle + 1;
            len -= half + 1;
          }
        else
          len = half;
      }
  }


  // --------------------- non-inline functions

  template <typename T>
  size_t
  pack(const T &          object,
       std::vector<char> &dest_buffer,
       const bool         allow_compression)
  {
    std::size_t size = 0;

    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
#ifdef DEAL_II_HAVE_CXX17
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#else
    if (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#endif
      {
        (void)allow_compression;
        const std::size_t previous_size = dest_buffer.size();
        dest_buffer.resize(previous_size + sizeof(T));

        std::memcpy(dest_buffer.data() + previous_size, &object, sizeof(T));

        size = sizeof(T);
      }
    else
      {
        // use buffer as the target of a compressing
        // stream into which we serialize the current object
        const std::size_t previous_size = dest_buffer.size();
        {
          boost::iostreams::filtering_ostreambuf fosb;
#ifdef DEAL_II_WITH_ZLIB
          if (allow_compression)
            fosb.push(boost::iostreams::gzip_compressor());
#else
          (void)allow_compression;
#endif
          fosb.push(boost::iostreams::back_inserter(dest_buffer));

          boost::archive::binary_oarchive boa(fosb);
          boa << object;
          // the stream object has to be destroyed before the return statement
          // to ensure that all data has been written in the buffer
        }
        size = dest_buffer.size() - previous_size;
      }

    return size;
  }


  template <typename T>
  std::vector<char>
  pack(const T &object, const bool allow_compression)
  {
    std::vector<char> buffer;
    pack<T>(object, buffer, allow_compression);
    return buffer;
  }


  template <typename T>
  T
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         const bool                               allow_compression)
  {
    T object;

    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
#ifdef DEAL_II_HAVE_CXX17
    if constexpr (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#else
    if (std::is_trivially_copyable<T>() && sizeof(T) < 256)
#endif
      {
        (void)allow_compression;
        Assert(std::distance(cbegin, cend) == sizeof(T), ExcInternalError());
        std::memcpy(&object, &*cbegin, sizeof(T));
      }
    else
      {
        // decompress the buffer section into the object
        boost::iostreams::filtering_istreambuf fisb;
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          fisb.push(boost::iostreams::gzip_decompressor());
#else
        (void)allow_compression;
#endif
        fisb.push(boost::iostreams::array_source(&*cbegin, &*cend));

        boost::archive::binary_iarchive bia(fisb);
        bia >> object;
      }

    return object;
  }


  template <typename T>
  T
  unpack(const std::vector<char> &buffer, const bool allow_compression)
  {
    return unpack<T>(buffer.cbegin(), buffer.cend(), allow_compression);
  }


  template <typename T, int N>
  void
  unpack(const std::vector<char>::const_iterator &cbegin,
         const std::vector<char>::const_iterator &cend,
         T (&unpacked_object)[N],
         const bool allow_compression)
  {
    // see if the object is small and copyable via memcpy. if so, use
    // this fast path. otherwise, we have to go through the BOOST
    // serialization machinery
    if (std::is_trivially_copyable<T>() && sizeof(T) * N < 256)
      {
        Assert(std::distance(cbegin, cend) == sizeof(T) * N,
               ExcInternalError());
        std::memcpy(unpacked_object, &*cbegin, sizeof(T) * N);
      }
    else
      {
        // decompress the buffer section into the object
        boost::iostreams::filtering_istreambuf fisb;
#ifdef DEAL_II_WITH_ZLIB
        if (allow_compression)
          fisb.push(boost::iostreams::gzip_decompressor());
#else
        (void)allow_compression;
#endif
        fisb.push(boost::iostreams::array_source(&*cbegin, &*cend));

        boost::archive::binary_iarchive bia(fisb);
        bia >> unpacked_object;
      }
  }


  template <typename T, int N>
  void
  unpack(const std::vector<char> &buffer,
         T (&unpacked_object)[N],
         const bool allow_compression)
  {
    unpack<T, N>(buffer.cbegin(),
                 buffer.cend(),
                 unpacked_object,
                 allow_compression);
  }



  inline bool
  get_bit(const unsigned char number, const unsigned int n)
  {
    AssertIndexRange(n, 8);

    // source:
    // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
    // "Checking a bit"
    return (number >> n) & 1U;
  }



  inline void
  set_bit(unsigned char &number, const unsigned int n, const bool x)
  {
    AssertIndexRange(n, 8);

    // source:
    // https://stackoverflow.com/questions/47981/how-do-you-set-clear-and-toggle-a-single-bit
    // "Changing the nth bit to x"
    number ^= (-static_cast<unsigned char>(x) ^ number) & (1U << n);
  }



  template <typename To, typename From>
  inline std::unique_ptr<To>
  dynamic_unique_cast(std::unique_ptr<From> &&p)
  {
    // Let's see if we can cast from 'From' to 'To'. If so, do the cast,
    // and then release the pointer from the old
    // owner
    if (To *cast = dynamic_cast<To *>(p.get()))
      {
        std::unique_ptr<To> result(cast);
        p.release();
        return result;
      }
    else
      throw std::bad_cast();
  }



  template <typename T>
  inline T &
  get_underlying_value(T &p)
  {
    return p;
  }



  template <typename T>
  inline T &
  get_underlying_value(std::shared_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(const std::shared_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(std::unique_ptr<T> &p)
  {
    return *p;
  }



  template <typename T>
  inline T &
  get_underlying_value(const std::unique_ptr<T> &p)
  {
    return *p;
  }



  template <typename Integer>
  std::vector<Integer>
  reverse_permutation(const std::vector<Integer> &permutation)
  {
    const std::size_t n = permutation.size();

    std::vector<Integer> out(n);
    for (std::size_t i = 0; i < n; ++i)
      out[i] = n - 1 - permutation[i];

    return out;
  }



  template <typename Integer>
  std::vector<Integer>
  invert_permutation(const std::vector<Integer> &permutation)
  {
    const std::size_t n = permutation.size();

    std::vector<Integer> out(n, numbers::invalid_unsigned_int);

    for (std::size_t i = 0; i < n; ++i)
      {
        AssertIndexRange(permutation[i], n);
        out[permutation[i]] = i;
      }

    // check that we have actually reached
    // all indices
    for (std::size_t i = 0; i < n; ++i)
      Assert(out[i] != numbers::invalid_unsigned_int,
             ExcMessage("The given input permutation had duplicate entries!"));

    return out;
  }
} // namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#ifndef DOXYGEN
namespace boost
{
  namespace serialization
  {
    // Provides boost and c++11 with a way to serialize tuples and pairs
    // automatically.
    template <int N>
    struct Serialize
    {
      template <class Archive, typename... Args>
      static void
      serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        ar &std::get<N - 1>(t);
        Serialize<N - 1>::serialize(ar, t, version);
      }
    };

    template <>
    struct Serialize<0>
    {
      template <class Archive, typename... Args>
      static void
      serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
      {
        (void)ar;
        (void)t;
        (void)version;
      }
    };

    template <class Archive, typename... Args>
    void
    serialize(Archive &ar, std::tuple<Args...> &t, const unsigned int version)
    {
      Serialize<sizeof...(Args)>::serialize(ar, t, version);
    }
  } // namespace serialization
} // namespace boost
#endif

#endif


