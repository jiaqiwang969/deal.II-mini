//include/deal.II-translator/base/point_0.txt
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

#ifndef dealii_point_h
#define dealii_point_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/tensor.h>

DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/geometries/point.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <cmath>

DEAL_II_NAMESPACE_OPEN

/**
 * 一个代表笛卡尔空间中一个点的类，其维度为  @p dim  。
 * 这个类的对象用于表示配备有<a
 * href="https://en.wikipedia.org/wiki/Cartesian_coordinate_system">Cartesian
 * coordinate
 * system</a>的向量空间的点（即锚定在原点的向量）。在其他用途中，它们被传递给在先验固定维度空间中对点进行操作的函数：与其使用<code>double
 * f(const double x)</code>和<code>double f(const double x, const double
 * y)</code>这样的函数，不如使用 <code>double f(const Point<dim>
 * &p)</code> ，因为它允许编写独立维度的代码。
 * deal.II特别使用Point对象来表示由笛卡尔坐标表示的点，也就是说，在
 * @p 二维空间维度中，一个点的特征是沿着由 @p dim
 * 相互正交的单位向量（称为
 * "坐标轴"）跨越的坐标系的轴的带符号距离。这种表示矢量的选择使得矢量的加法和缩放特别简单：人们只需要对每个坐标值进行加法或乘法。另一方面，当一个向量用其他类型的坐标系表示时（如<a
 * href="https://en.wikipedia.org/wiki/Spherical_coordinate_system">spherical
 * coordinate
 * systems</a>），添加或缩放向量就没有那么简单了。
 *
 *  <h3>What's a <code>Point@<dim@></code> and what is a
 * <code>Tensor@<1,dim@></code>?</h3> 点类派生于张量 @<1,dim@>
 * ，因此共享后者的成员函数和其他属性。事实上，它本身的附加函数相对较少（最明显的例外是计算空间中两点之间欧几里得距离的distance()函数），因此这两个类通常可以互换使用。
 * 尽管如此，还是有语义上的差异，使得我们在不同的、定义明确的语境中使用这些类。在deal.II中，我们用<tt>Point</tt>类来表示空间中的点，即表示
 * <em> 锚定在原点 </em>
 * 的向量（秩-1张量）。另一方面，锚定在其他地方的向量（因此在这个词的通常用法中不代表
 * <em> 点 </em> ）由张量 @<1,dim@>. 类型的对象表示。 ]
 * 特别是方向向量、法向量、梯度和两点之间的差值（即当你从一个点减去另一个点时得到的东西）：所有这些都由Tensor
 * @<1,dim@> 对象而不是Point @<dim@>. 表示。
 * 此外，点类只用于对象的坐标可以被认为拥有长度维度的地方。一个表示物体的重量、高度和成本的对象既不是点也不是张量（因为它缺乏坐标系旋转下的变换属性），因此不应该用这些类来表示。在这种情况下，使用一个大小为3的数组，或者使用
 * <code>std::array</code>
 * 类。另外，就像在矢量值函数的情况下，你可以使用矢量类型的对象或
 * <code>std::vector</code>  。
 *
 * @tparam  dim
 * 一个整数，表示一个点所在的空间的维度。当然，这等于确定一个点的坐标数。
 * @tparam  Number
 * 用于存储坐标值的数据类型。几乎在所有情况下，这都是默认的
 * @p double,
 * ，但在某些情况下，人们可能希望以不同的（而且总是标量的）类型来存储坐标。一个例子是一个区间类型，它可以存储一个坐标的值以及它的不确定性。另一个例子是一个允许自动微分的类型（例如，见
 * step-33
 * 中使用的Sacado类型），从而可以在传递一个坐标被存储在这种类型中的点对象时产生函数的解析（空间）导数。
 *
 *
 *
 * @ingroup geomprimitives
 *
 */
template <int dim, typename Number = double>
class Point : public Tensor<1, dim, Number>
{
public:
  /**
   * 标准构造函数。创建一个对应于原点的对象，即所有的坐标都设置为零。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV
  Point();

  /**
   * 将一个张量转换为一个点。
   *
   */
  explicit DEAL_II_CUDA_HOST_DEV
  Point(const Tensor<1, dim, Number> &);

  /**
   * 一维点的构造函数。这个函数只对<tt>dim==1</tt>实现，因为对<tt>dim!=1</tt>的点的使用被认为是不安全的，因为它将使点坐标的一些分量未被初始化。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  explicit DEAL_II_CUDA_HOST_DEV
  Point(const Number x);

  /**
   * 二维点的构造函数。这个函数只对<tt>dim==2</tt>实现，因为对于<tt>dim!=2</tt>的点来说，这个用法被认为是不安全的，因为它将使点坐标的某些分量未被初始化（如果dim>2）或不使用某些参数（如果dim<2）。
   * @note 这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV
  Point(const Number x, const Number y);

  /**
   * 三维点的构造函数。这个函数只对<tt>dim==3</tt>实现，因为对于<tt>dim!=3</tt>的点来说，这个用法被认为是不安全的，因为它将使点坐标的一些分量未被初始化（如果dim>3）或者不使用一些参数（如果dim<3）。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV
  Point(const Number x, const Number y, const Number z);

  /**
   * 将一个 boost::geometry::point 转换为一个 dealii::Point. 。
   *
   */
  template <std::size_t dummy_dim,
            typename std::enable_if<(dim == dummy_dim) && (dummy_dim != 0),
                                    int>::type = 0>
  Point(const boost::geometry::model::
          point<Number, dummy_dim, boost::geometry::cs::cartesian> &boost_pt);

  /**
   * 返回一个坐标方向<tt>i</tt>的单位向量，即除了<tt>i</tt>第1个坐标中的一个1之外，所有坐标中的向量都是0。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  static DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                               unit_vector(const unsigned int i);

  /**
   * 对<tt>index</tt>th坐标的读取访问。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Number
                        operator()(const unsigned int index) const;

  /**
   * 对<tt>index</tt>th坐标的读和写访问。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Number &
                        operator()(const unsigned int index);

  /**
   * 来自Tensor<1, dim,
   * Number>的赋值操作，其底层标量类型不同。这显然要求 @p
   * OtherNumber 类型可以转换为 @p Number. 。
   *
   */
  template <typename OtherNumber>
  Point<dim, Number> &
  operator=(const Tensor<1, dim, OtherNumber> &p);

  /**
   * @name  点的加法和减法。    @{
   *
   */

  /**
   * 将一个以张量<1,dim,Number>给出的偏移量添加到一个点上。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator+(const Tensor<1, dim, Number> &) const;

  /**
   * 减去两个点，即获得连接这两个点的向量。正如在这个类的文档中所讨论的，减去两个点的结果是一个锚定在两个点之一的向量（而不是原点），因此，结果是作为一个张量
   * @<1,dim@> 而不是作为一个点 @<dim@>. 返回。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Tensor<1, dim, Number>
                        operator-(const Point<dim, Number> &) const;

  /**
   * 从当前点减去一个差分向量（用张量 @<1,dim@>)
   * 表示）。这将产生另一个点，正如在这个类的文档中所讨论的，然后结果自然是作为一个点
   * @<dim@> 对象而不是作为一个张量 @<1,dim@>. 返回。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator-(const Tensor<1, dim, Number> &) const;

  /**
   * 相反的向量。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                        operator-() const;

  /**
   * @}
   *
   */

  /**
   * @name  点的乘法和缩放。点积。规范。    @{
   *
   */

  /**
   * 将当前点乘以一个系数。
   * @note  这个函数也可以在CUDA设备代码中使用。
   * @relatesalso  EnableIfScalar
   *
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>
  operator*(const OtherNumber) const;

  /**
   * 将当前点除以一个系数。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  template <typename OtherNumber>
  DEAL_II_CUDA_HOST_DEV Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>
  operator/(const OtherNumber) const;

  /**
   * 返回代表两点的向量的标量积。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV Number operator*(const Tensor<1, dim, Number> &p) const;

  /**
   * 返回该点向量与自身的标量乘积，即平方，或规范的平方。如果是复数类型，则相当于此点向量与自身的复共轭的收缩。
   * @note  这个函数等同于 Tensor<rank,dim,Number>::norm_square()
   * ，它返回弗罗本纽斯法线的平方。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  square() const;

  /**
   * 返回<tt>this</tt>点到<tt>p</tt>点的欧氏距离，即代表两点的向量之差的
   * $l_2$ 规范。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  distance(const Point<dim, Number> &p) const;

  /**
   * 返回<tt>this</tt>点到<tt>p</tt>点的平方欧氏距离。
   * @note  这个函数也可以在CUDA设备代码中使用。
   *
   */
  DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
  distance_square(const Point<dim, Number> &p) const;

  /**
   * @}
   *
   */

  /**
   * 使用[BOOST序列化库](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html)将此对象的数据读入或写入一个流中，以便进行序列化。
   *
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int version);
};

 /*--------------------------- Inline functions: Point -----------------------*/ 

#ifndef DOXYGEN

// At least clang-3.7 requires us to have a user-defined constructor
// and we can't use 'Point<dim,Number>::Point () = default' here.
template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point() // NOLINT
{}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Tensor<1, dim, Number> &t)
  : Tensor<1, dim, Number>(t)
{}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 1,
         ExcMessage(
           "You can only initialize Point<1> objects using the constructor "
           "that takes only one argument. Point<dim> objects with dim!=1 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the switch anyway so
  // as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  switch (dim)
    {
      case 1:
        this->values[0] = x;
        break;

      default:;
    }
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x, const Number y)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 2,
         ExcMessage(
           "You can only initialize Point<2> objects using the constructor "
           "that takes two arguments. Point<dim> objects with dim!=2 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the indirection anyway
  // so as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  this->values[0]                = x;
  this->values[y_index]          = y;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV
Point<dim, Number>::Point(const Number x, const Number y, const Number z)
{
#  ifndef __CUDA_ARCH__
  Assert(dim == 3,
         ExcMessage(
           "You can only initialize Point<3> objects using the constructor "
           "that takes three arguments. Point<dim> objects with dim!=3 "
           "require initialization with the constructor that takes 'dim' "
           "arguments."));
#  endif

  // we can only get here if we pass the assertion. use the indirection anyway
  // so as to avoid compiler warnings about uninitialized elements or writing
  // beyond the end of the 'values' array
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  constexpr unsigned int z_index = (dim < 3) ? 0 : 2;
  this->values[0]                = x;
  this->values[y_index]          = y;
  this->values[z_index]          = z;
}



template <int dim, typename Number>
template <
  std::size_t dummy_dim,
  typename std::enable_if<(dim == dummy_dim) && (dummy_dim != 0), int>::type>
inline Point<dim, Number>::Point(
  const boost::geometry::model::
    point<Number, dummy_dim, boost::geometry::cs::cartesian> &boost_pt)
{
  Assert(dim <= 3, ExcNotImplemented());
  this->values[0]                = boost::geometry::get<0>(boost_pt);
  constexpr unsigned int y_index = (dim < 2) ? 0 : 1;
  constexpr unsigned int z_index = (dim < 3) ? 0 : 2;

  if (dim >= 2)
    this->values[y_index] = boost::geometry::get<y_index>(boost_pt);

  if (dim >= 3)
    this->values[z_index] = boost::geometry::get<z_index>(boost_pt);
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
                             Point<dim, Number>::unit_vector(unsigned int i)
{
  Point<dim, Number> p;
  p[i] = 1.;
  return p;
}


template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number
Point<dim, Number>::operator()(const unsigned int index) const
{
#  ifndef __CUDA_ARCH__
  AssertIndexRange(static_cast<int>(index), dim);
#  endif
  return this->values[index];
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number &
Point<dim, Number>::operator()(const unsigned int index)
{
#  ifndef __CUDA_ARCH__
  AssertIndexRange(static_cast<int>(index), dim);
#  endif
  return this->values[index];
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE Point<dim, Number> &
Point<dim, Number>::operator=(const Tensor<1, dim, OtherNumber> &p)
{
  Tensor<1, dim, Number>::operator=(p);
  return *this;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator+(const Tensor<1, dim, Number> &p) const
{
  Point<dim, Number> tmp = *this;
  tmp += p;
  return tmp;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Tensor<1, dim, Number>
Point<dim, Number>::operator-(const Point<dim, Number> &p) const
{
  return (Tensor<1, dim, Number>(*this) -= p);
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator-(const Tensor<1, dim, Number> &p) const
{
  Point<dim, Number> tmp = *this;
  tmp -= p;
  return tmp;
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Point<dim, Number>
Point<dim, Number>::operator-() const
{
  Point<dim, Number> result;
  for (unsigned int i = 0; i < dim; ++i)
    result.values[i] = -this->values[i];
  return result;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
    Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
    Point<dim, Number>::operator*(const OtherNumber factor) const
{
  Point<dim, typename ProductType<Number, OtherNumber>::type> tmp;
  for (unsigned int i = 0; i < dim; ++i)
    tmp[i] = this->operator[](i) * factor;
  return tmp;
}



template <int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
  Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
  Point<dim, Number>::operator/(const OtherNumber factor) const
{
  const Tensor<1, dim, Number> &base_object = *this;
  return Point<
    dim,
    typename ProductType<Number,
                         typename EnableIfScalar<OtherNumber>::type>::type>(
    dealii::operator/(base_object, factor));
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV Number Point<dim, Number>::
                                    operator*(const Tensor<1, dim, Number> &p) const
{
  Number res = Number();
  for (unsigned int i = 0; i < dim; ++i)
    res += this->operator[](i) * p[i];
  return res;
}


template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::square() const
{
  return this->norm_square();
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::distance(const Point<dim, Number> &p) const
{
  return std::sqrt(distance_square(p));
}



template <int dim, typename Number>
inline DEAL_II_CUDA_HOST_DEV typename numbers::NumberTraits<Number>::real_type
Point<dim, Number>::distance_square(const Point<dim, Number> &p) const
{
  Number sum = internal::NumberType<Number>::value(0.0);
  for (unsigned int i = 0; i < dim; ++i)
    {
      const Number diff = static_cast<Number>(this->values[i]) - p(i);
      sum += numbers::NumberTraits<Number>::abs_square(diff);
    }

  return sum;
}



template <int dim, typename Number>
template <class Archive>
inline void
Point<dim, Number>::serialize(Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class
  ar &static_cast<Tensor<1, dim, Number> &>(*this);
}

#endif // DOXYGEN


 /*--------------------------- Global functions: Point -----------------------*/ 


/**
 * 用标量缩放点向量的全局运算符。
 *
 *
 * @note  这个函数也可以在CUDA设备代码中使用。
 * @relatesalso  点  @relatesalso  EnableIfScalar
 *
 *
 */
template <int dim, typename Number, typename OtherNumber>
inline DEAL_II_CUDA_HOST_DEV
  Point<dim,
        typename ProductType<Number,
                             typename EnableIfScalar<OtherNumber>::type>::type>
  operator*(const OtherNumber factor, const Point<dim, Number> &p)
{
  return p * factor;
}



/**
 * 点的输出运算符。连续地打印元素，中间有一个空格。
 * @relatesalso  点
 *
 *
 */
template <int dim, typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Point<dim, Number> &p)
{
  for (unsigned int i = 0; i < dim - 1; ++i)
    out << p[i] << ' ';
  out << p[dim - 1];

  return out;
}



/**
 * 点的输入运算符。连续地输入元素。  @relatesalso  点
 *
 *
 */
template <int dim, typename Number>
inline std::istream &
operator>>(std::istream &in, Point<dim, Number> &p)
{
  for (unsigned int i = 0; i < dim; ++i)
    in >> p[i];

  return in;
}


#ifndef DOXYGEN

/**
 * 维度为1的点的输出运算符。这是从一般模板中专门实现的，以避免编译器警告该循环是空的。
 *
 *
 */
template <typename Number>
inline std::ostream &
operator<<(std::ostream &out, const Point<1, Number> &p)
{
  out << p[0];

  return out;
}

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif


