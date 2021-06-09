//include/deal.II-translator/numerics/histogram_0.txt
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

#ifndef dealii_histogram_h
#define dealii_histogram_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN


/**
 * 这个类提供了一些生成2D和3D直方图的设施。它的使用方法是给它一个或几个数据集和一个如何将其中的数值范围分成区间的规则（例如线性间隔或对数间隔的区间）。然后，这些值被分类到不同的区间，每个区间的值的数量被储存起来，以便以后输出。如果只给了一个数据集，产生的直方图将是一个2D的，而如果给了一个以上的数据集，它将是一个3D的。对于一个以上的数据集，每个数据集都使用相同的区间，以便于比较。
 *
 *  <h3>Ways to generate the intervals</h3>
 * 目前，实现了以下的间隔方案。  <ul>   <li>  线性间隔。间隔在数据的最小值和最大值之间以恒定的步骤分布。  <li>  对数间隔。间隔在数值的最小值和最大值之间以恒定的步长分布。这个方案只有在数据只有正值的情况下才有用。负值和零值被排序到最左边的区间。  </ul>
 * 为了保持程序的可扩展性，你可以使用两个函数 @p
 * get_interval_spacing_names和 @p parse_interval_spacing,
 * ，它们总是给你一个目前支持的间隔格式的完整列表，并且能够生成
 * @p enum. 的相应值。
 * 如果你使用它们，你可以以这样的方式编写你的程序，即它只需要重新编译以使新增加的格式生效，而无需改变代码。
 *
 *  <h3>Output formats</h3> 目前，只支持GNUPLOT输出。
 *
 *
 *
 * @ingroup textoutput
 *
 *
 */
class Histogram
{
public:
  /**
   * 定义了几种安排间隔的方法。
   *
   */
  enum IntervalSpacing
  {
    /**
     * 以线性方式排列间隔。
     *
     */
    linear,
    /**
     * 以对数方式安排间隔。
     *
     */
    logarithmic
  };


  /**
   * 取几个值的列表，每个列表上产生一个直方图，然后将其一个个排列在后面。
   * 同时使用几个数据集可以更容易地进行比较，因为数据排序的区间对所有数据集都是一样的。
   * 直方图的排列方式是：<tt>值[i][j]<tt>的计算区间构成x区间，每个区间的值的数量是y区间（对于2d图）或z区间（对于3d图）。对于3D绘图，
   * @p y_values
   * 参数用来给每个数据集分配一个y方向的值，也就是生成的绘图中的深度坐标。对于2D绘图，
   * @p y_values 被忽略。
   * 如果你只给出一个数据集，即<tt>values.size()==1</tt>，那么得到的直方图将是一个2D的直方图。
   * @p n_intervals  表示数据将被排序的区间数； @p
   * interval_spacing
   * 表示计算区间边界的方式。关于这方面的更多信息，请参考通用文档。
   *
   */
  template <typename number>
  void
  evaluate(const std::vector<Vector<number>> &values,
           const std::vector<double> &        y_values,
           const unsigned int                 n_intervals,
           const IntervalSpacing              interval_spacing = linear);

  /**
   * 如果你只有一个数据集，这个函数只是上面那个函数的一个封装器。
   *
   */
  template <typename number>
  void
  evaluate(const Vector<number> &values,
           const unsigned int    n_intervals,
           const IntervalSpacing interval_spacing = linear);

  /**
   * 将 @p evaluate
   * 函数计算出的直方图以适合GNUPLOT程序的格式写入一个流。该函数可生成2D或3D直方图。
   *
   */
  void
  write_gnuplot(std::ostream &out) const;

  /**
   * 以字符串形式返回允许的区间间隔名称。目前是
   * "线性|logarithmic"。
   *
   */
  static std::string
  get_interval_spacing_names();

  /**
   * 获取一个包含上述函数返回的名称之一的字符串，并返回
   * @p IntervalSpacing. 的相应值
   * 如果该字符串不是有效的，则抛出一个错误。
   *
   */
  static IntervalSpacing
  parse_interval_spacing(const std::string &name);

  /**
   * 确定这个对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 异常情况。
   *
   */
  DeclExceptionMsg(ExcEmptyData,
                   "Your input argument to this function does not appear to "
                   "have any data in it.");
  /**
   * 异常情况。
   *
   */
  DeclException2(ExcIncompatibleArraySize,
                 int,
                 int,
                 << "The two array sizes " << arg1 << " and " << arg2
                 << " must match, but don't.");
  /**
   * 异常情况。
   *
   */
  DeclException1(ExcInvalidName,
                 std::string,
                 << "The given name <" << arg1
                 << "> does not match any of the known formats.");

private:
  /**
   * 表示一个区间的结构。
   *
   */
  struct Interval
  {
    /**
     * 构造函数。设置边界，并将这个区间的值的数量设置为零。
     *
     */
    Interval(const double left_point, const double right_point);

    /**
     * 确定这个对象的内存消耗（以字节为单位）的估计值。
     *
     */
    std::size_t
    memory_consumption() const;

    /**
     * 区间的左边界。
     *
     */
    double left_point;

    /**
     * 区间的右边界。
     *
     */
    double right_point;

    /**
     * 这个区间的值的数量。
     *
     */
    unsigned int content;
  };

  /**
   * "小于
   * "操作，通过对零和负值的排序找到最小的正值，使其大于最大的正数。
   * 用于寻找对数情况下区间间隔方案中最左边区间的下限。
   * 返回  @p true,
   * 如果（<tt>n1<n2</tt>，且（<tt>n1>0</tt>或<tt>n2<0</tt>）），或（n2<n1且n1>0且n2<=0）。这实际上是将所有的负数排序为比最大的正数大。
   *
   */
  template <typename number>
  static bool
  logarithmic_less(const number n1, const number n2);

  /**
   * 为给 @p 评估函数的每个数据集保存一组区间的向量。
   *
   */
  std::vector<std::vector<Interval>> intervals;

  /**
   * 3d直方图的深度轴的值。存储在 @p evaluate 函数中。
   *
   */
  std::vector<double> y_values;
};


DEAL_II_NAMESPACE_CLOSE

#endif


