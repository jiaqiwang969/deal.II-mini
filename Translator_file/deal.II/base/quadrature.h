//include/deal.II-translator/base/quadrature_0.txt
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

#ifndef dealii_quadrature_h
#define dealii_quadrature_h


#include <deal.II/base/config.h>

#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>

#include <array>
#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup Quadrature */ 
 /*@{*/ 

/**
 * 用于任意维度的正交公式的基类。该类存储单位线[0,1]、单位方格[0,1]x[0,1]等上的正交点和权重。
 * 有许多派生类，表示具体的积分公式。它们的名字以<tt>Q</tt>为前缀。更多细节请参考派生类的列表。
 * 更高维度的方案通常是一维公式的张量乘积，但请参考下面关于实现细节的部分。
 * 为了允许独立维度的编程，存在一个零维的正交公式。由于零维上的积分是在单点上的求值，这样一个公式的任何构造函数都初始化为一个权重为1的单点正交。对权重的访问是可能的，而对正交点的访问是不允许的，因为零维的点不包含任何信息。这些公式的主要目的是在QProjector中使用，它将从这些公式中创建一个有用的一维公式。
 * <h3>Mathematical background</h3>
 * 对于每个正交公式，我们用<tt>m</tt>表示，确切集成的多项式的最大程度。这个数字在每个公式的文档中给出。积分误差的阶数是<tt>m+1</tt>，也就是说，根据Bramble-
 * Hilbert定理，误差是到<tt>m+1</tt>的单元的大小。这个数字<tt>m</tt>可以在每个具体公式的文件中找到。对于最优公式QGauss，我们有
 * $m = 2N-1$
 * ，其中N是QGauss的构造器参数。张量积公式在每个空间方向上对度数为<tt>m</tt>的张量积多项式是精确的，但它们仍然只有<tt>m+1</tt>st
 * order。 <h3>Implementation details</h3>
 * 大多数在一个以上空间维度的积分公式是一个空间维度的正交公式的张量乘积，或者更普遍的是<tt>(dim-1)</tt>维度的公式与一个维度的公式的张量乘积。有一个特殊的构造函数可以从其他两个公式中生成一个正交公式。
 * 例如，QGauss  @<dim@>
 * 公式包括<tt>dim</tt>维度上的<i>N<sup>dim</sup></i>正交点，其中N是QGauss的构造参数。
 *
 *
 * @note
 * 这个模板的实例化提供给维度0、1、2和3（见 @ref
 * Instantiations 部分）。
 *
 *
 */
template <int dim>
class Quadrature : public Subscriptor
{
public:
  /**
   * 为作用于少一维的对象的正交定义一个别名。对于单元格来说，这将是一个面的四分法。
   *
   */
  using SubQuadrature = Quadrature<dim - 1>;

  /**
   * 构造函数。    这个构造函数被标记为显式，以避免像
   * <code>hp::QCollection@<dim@> q_collection(3)</code> 中的
   * <code>hp::QCollection@<dim@> q_collection(QGauss@<dim@>(3))</code>
   * 那样的非自愿事故。
   *
   */
  explicit Quadrature(const unsigned int n_quadrature_points = 0);

  /**
   * 将此正交公式构建为比现在少一维的公式与一维的公式的张量乘积。
   * 这个构造函数假设（并测试）常数函数被精确整合，即正交权重之和为1。
   * <tt>SubQuadrature<dim>::type</tt> 扩展为<tt>正交<dim-1></tt>。
   *
   */
  Quadrature(const SubQuadrature &, const Quadrature<1> &);

  /**
   * 将这个正交公式构建为一维公式的<tt>dim</tt>折张量积。
   * 假设一维规则中的点是按升序排列的，那么所产生的规则的点是按词典排序的，其中<i>x</i>运行最快。
   * 为了避免与1d中的复制构造器冲突，我们让参数在dim==1时为0d正交公式，在所有其他空间维度上为1d正交公式。
   * 这个构造函数并不要求常数函数被精确地整合。因此，如果一维公式是相对于一个加权函数定义的，它是合适的。
   *
   */
  explicit Quadrature(const Quadrature<dim != 1 ? 1 : 0> &quadrature_1d);

  /**
   * 复制构造函数。
   *
   */
  Quadrature(const Quadrature<dim> &q);

  /**
   * 移动构造函数。通过转移另一个正交对象的内部数据来构造一个新的正交对象。
   *
   */
  Quadrature(Quadrature<dim> &&) noexcept = default;

  /**
   * 从给定的正交点向量（实际上应该在单元格中）和相应的权重构造一个正交公式。
   * 你会希望权重之和为1，但这并不被检查。
   *
   */
  Quadrature(const std::vector<Point<dim>> &points,
             const std::vector<double> &    weights);

  /**
   * 从一个点的列表中构建一个假正交公式，权重设置为无穷大。因此，产生的对象并不是为了实际进行积分，而是与FEValues对象一起使用，以便找到一些点（此对象中的正交点）在实空间中转换后的单元上的位置。
   *
   */
  Quadrature(const std::vector<Point<dim>> &points);

  /**
   * 一个单点正交的构造函数。设置这个点的权重为1。
   *
   */
  Quadrature(const Point<dim> &point);

  /**
   * 虚拟解构器。
   *
   */
  virtual ~Quadrature() override = default;

  /**
   * 赋值运算符。复制#weights和#quadrature_points的内容以及大小。
   *
   */
  Quadrature &
  operator=(const Quadrature<dim> &);

  /**
   * 移动赋值运算符。将所有数据从另一个正交对象移动到这个对象。
   *
   */
  Quadrature &
  operator=(Quadrature<dim> &&) = default; // NOLINT

  /**
   * 测试两个正交点的相等性。
   *
   */
  bool
  operator==(const Quadrature<dim> &p) const;

  /**
   * 设置正交点和权重为参数中提供的值。
   *
   */
  void
  initialize(const std::vector<Point<dim>> &points,
             const std::vector<double> &    weights);

  /**
   * 正交点的数量。
   *
   */
  unsigned int
  size() const;

  /**
   * 返回<tt>i</tt>第1个正交点。
   *
   */
  const Point<dim> &
  point(const unsigned int i) const;

  /**
   * 返回对整个正交点数组的引用。
   *
   */
  const std::vector<Point<dim>> &
  get_points() const;

  /**
   * 返回<tt>i</tt>个正交点的权重。
   *
   */
  double
  weight(const unsigned int i) const;

  /**
   * 返回对整个权重数组的引用。
   *
   */
  const std::vector<double> &
  get_weights() const;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   *
   */
  std::size_t
  memory_consumption() const;

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据写入或读出到一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);

  /**
   * 如果正交对象是一维公式的张量乘积，并且正交点是按字母顺序排序的，则该函数返回真。
   *
   */
  bool
  is_tensor_product() const;

  /**
   * 如果正交公式是张量积，此函数返回 @p dim
   * 一维基础对象。  否则，不允许调用此函数。    对于 @p
   * dim 等于1的情况，我们不能将 std::array
   * 作为常量引用返回，必须通过值来返回。在这种情况下，数组将总是包含一个元素（
   * @p this).  ）。
   * @note  这个函数的实际返回类型是
   * @code
   * std::conditional<dim == 1,
   *                std::array<Quadrature<1>, dim>,
   *                const std::array<Quadrature<1>, dim> &>::type
   * @endcode
   * 在线文档中对该类型进行了缩写，以提高本页面的可读性。
   *
   */
#ifndef DOXYGEN
  typename std::conditional<dim == 1,
                            std::array<Quadrature<1>, dim>,
                            const std::array<Quadrature<1>, dim> &>::type
#else
  const std::array<Quadrature<1>, dim> &
#endif
  get_tensor_basis() const;

protected:
  /**
   * 正交点的列表。要由派生类的构造函数来填充。
   *
   */
  std::vector<Point<dim>> quadrature_points;

  /**
   * 正交点的权重列表。 将由派生类的构造函数填写。
   *
   */
  std::vector<double> weights;

  /**
   * 指示此对象是否代表正交公式是一维公式的张量乘积。
   * 如果dim==1或者调用了取自Quadrature<1>（可能还有Quadrature<dim-1>对象）的构造函数，则该标志被设置。这意味着正交点是按字母顺序排序的。
   *
   */
  bool is_tensor_product_flag;

  /**
   * 存储一维张量基础对象，以备此对象可由张量积表示。
   *
   */
  std::unique_ptr<std::array<Quadrature<1>, dim>> tensor_basis;
};


/**
 * 正交公式，实现参考单元上正交点的各向异性分布。为此，生成<tt>dim</tt>一维正交公式的张量积。
 *
 *
 * @note
 * 每个构造函数只能在与参数数相匹配的维度上使用。
 *
 *
 */
template <int dim>
class QAnisotropic : public Quadrature<dim>
{
public:
  /**
   * 一个一维公式的构造函数。这个构造函数只是复制给定的正交规则。
   *
   */
  QAnisotropic(const Quadrature<1> &qx);

  /**
   * 二维公式的构造函数。
   *
   */
  QAnisotropic(const Quadrature<1> &qx, const Quadrature<1> &qy);

  /**
   * 三维公式的构造函数。
   *
   */
  QAnisotropic(const Quadrature<1> &qx,
               const Quadrature<1> &qy,
               const Quadrature<1> &qz);
};


/**
 * 通过在每个方向上迭代另一个正交公式构建的正交公式。在一个以上的空间维度中，所得到的正交公式是通过在一个空间维度中建立各自迭代的正交公式的张量积而以通常的方式构建的。
 * 在一个空间维度上，给定的基础公式被复制和缩放到给定数量的长度为<tt>1/n_copies</tt>的子区间上。如果正交公式使用了单位区间的两个端点，那么在迭代后的正交公式内部会有两次使用的正交点；我们将它们合并为一个，其权重为最左和最右的正交点的权重之和。
 * 由于所有高于一维的维度都是由一维和<tt>dim-1</tt>维的正交公式的张量积建立起来的，所以给构造函数的参数需要是一个空间维度的正交公式，而不是<tt>dim</tt>维。
 * 本类的目的是提供一个低阶公式，其中误差常数可以通过增加正交点的数量来调整。这在整合单元格上的非微分函数时很有用。
 *
 *
 */
template <int dim>
class QIterated : public Quadrature<dim>
{
public:
  /**
   * 构造函数。在每个方向上迭代给定的正交公式<tt>n_copies</tt>次。
   *
   */
  QIterated(const Quadrature<1> &base_quadrature, const unsigned int n_copies);

  /**
   * 异常情况
   *
   */
  DeclExceptionMsg(ExcInvalidQuadratureFormula,
                   "The quadrature formula you provided cannot be used "
                   "as the basis for iteration.");
};



 /*@}*/ 

#ifndef DOXYGEN

// -------------------  inline and template functions ----------------


template <int dim>
inline unsigned int
Quadrature<dim>::size() const
{
  return weights.size();
}


template <int dim>
inline const Point<dim> &
Quadrature<dim>::point(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return quadrature_points[i];
}



template <int dim>
double
Quadrature<dim>::weight(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return weights[i];
}



template <int dim>
inline const std::vector<Point<dim>> &
Quadrature<dim>::get_points() const
{
  return quadrature_points;
}



template <int dim>
inline const std::vector<double> &
Quadrature<dim>::get_weights() const
{
  return weights;
}



template <int dim>
inline bool
Quadrature<dim>::is_tensor_product() const
{
  return is_tensor_product_flag;
}



template <int dim>
template <class Archive>
inline void
Quadrature<dim>::serialize(Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class.
  ar &static_cast<Subscriptor &>(*this);

  ar &quadrature_points &weights;
}



 /* -------------- declaration of explicit specializations ------------- */ 

template <>
Quadrature<0>::Quadrature(const unsigned int);
template <>
Quadrature<0>::Quadrature(const Quadrature<-1> &, const Quadrature<1> &);
template <>
Quadrature<0>::Quadrature(const Quadrature<1> &);
template <>
Quadrature<0>::Quadrature(const Point<0> &);

template <>
Quadrature<1>::Quadrature(const Quadrature<0> &, const Quadrature<1> &);

template <>
Quadrature<1>::Quadrature(const Quadrature<0> &);

template <>
QIterated<1>::QIterated(const Quadrature<1> &base_quadrature,
                        const unsigned int   n_copies);

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


