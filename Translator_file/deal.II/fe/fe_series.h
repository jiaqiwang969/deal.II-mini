//include/deal.II-translator/fe/fe_series_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2021 by the deal.II authors
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

#ifndef dealii_fe_series_h
#define dealii_fe_series_h



#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/tensor.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/vector_tools_common.h>

#include <memory>
#include <string>
#include <vector>


DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup feall */ 
 /*@{*/ 


/**
 * 这个命名空间提供了计算参考元素上的解的扩展系列的函数。膨胀系数通常用于估计底层FiniteElement场的局部平滑性，以决定h-或p-适应性细化策略。
 *
 *
 */
namespace FESeries
{
  /**
   * 一个用于计算标量FE（或矢量FE的单一分量）场在参考元素上的傅里叶级数的扩展的类。  傅里叶级数的指数形式是基于指数函数集的完整性和赫米特正交性  $ \phi_{\bf k}({\bf x}) = \exp(2 \pi i\, {\bf k} \cdot {\bf x})$  。  例如在一维中，L2正交性条件为 @f[
   * \int_0^1 \phi_k(x) \phi_l^\ast(x) dx=\delta_{kl}.
   * @f] 注意  $ \phi_{\bf k} = \phi_{-\bf k}^\ast $  。    参考元素上的任意标量FE场可以在完整的正交指数基中展开为@f[
   * u({\bf x}) = \sum_{\bf k} c_{\bf k} \phi_{\bf k}({\bf x}).
   * @f] 从基的正交性属性来看，可以得出@f[
   * c_{\bf k} = \int_{[0,1]^d} u({\bf x}) \phi_{\bf k}^\ast ({\bf x}) d{\bf
   * x}\,. @f]
   * 正是这个复值展开系数，由这个类计算。请注意， $
   * u({\bf x}) = \sum_i u_i N_i({\bf x})$  ，其中 $ N_i({\bf x}) $
   * 是实值的FiniteElement形状函数。  因此 $ c_{\bf k} \equiv
   * c_{-\bf k}^\ast $ ，我们只需要计算 $ c_{\bf k} $ 的正指数 $
   * \bf k $  。
   *
   */
  template <int dim, int spacedim = dim>
  class Fourier : public Subscriptor
  {
  public:
    using CoefficientType = typename std::complex<double>;

    /**
     * 构造函数，初始化所有需要的数据结构。         @p
     * n_coefficients_per_direction 定义了每个方向上的系数数， @p
     * fe_collection 是将用于扩展的 hp::FECollection ， @p
     * q_collection 是用于整合每个FiniteElement在 @p fe_collection.
     * 的扩展的 hp::QCollection
     * 由于傅里叶扩展只能对标量场进行，这个类不能对矢量值的有限元操作，因此会抛出一个断言。然而，有限元场的每个分量可以分别被视为标量场，对其进行傅里叶展开也是可能的。为此，可选的参数
     * @p component 定义了每个有限元的哪个分量将被使用。
     * @p component
     * 的默认值只适用于标量FEs，在这种情况下，它表示唯一的分量将被分解。对于矢量FE，必须明确提供一个非默认值。
     *
     */
    Fourier(const std::vector<unsigned int> &      n_coefficients_per_direction,
            const hp::FECollection<dim, spacedim> &fe_collection,
            const hp::QCollection<dim> &           q_collection,
            const unsigned int component = numbers::invalid_unsigned_int);

    /**
     * 一个非默认的构造函数。 @p n_coefficients_per_direction
     * 定义了每个方向的模数， @p fe_collection
     * 是将使用扩展的 hp::FECollection ， @p q_collection
     * 是用于整合 @p fe_collection.  @deprecated
     * 中每个FiniteElement的扩展的 hp::QCollection
     * 使用不同的构造函数代替。
     *
     */
    DEAL_II_DEPRECATED
    Fourier(const unsigned int                     n_coefficients_per_direction,
            const hp::FECollection<dim, spacedim> &fe_collection,
            const hp::QCollection<dim> &           q_collection);

    /**
     * 计算 @p fourier_coefficients 中 @p local_dof_values
     * 给出的单元格矢量场，对应于FiniteElement的 @p
     * cell_active_fe_index  。
     *
     */
    template <typename Number>
    void
    calculate(const dealii::Vector<Number> &local_dof_values,
              const unsigned int            cell_active_fe_index,
              Table<dim, CoefficientType> & fourier_coefficients);

    /**
     * 返回与 @p index 相关的有限元在所提供的 hp::FECollection.
     * 中每个坐标方向的系数数。
     *
     */
    unsigned int
    get_n_coefficients_per_direction(const unsigned int index) const;

    /**
     * 计算所有的转换矩阵，将有限元解转移到系列展开表示。
     * 这些矩阵将通过调用calculate()按需生成，并存储起来供重复使用。通常情况下，这个操作会消耗大量的工作量。有了这个函数，所有的矩阵将被提前计算。
     * 这样，我们就可以将其昂贵的生成与实际应用分开。
     *
     */
    void
    precalculate_all_transformation_matrices();

    /**
     * 将此对象的所有转换矩阵写入一个流中，以便进行序列化。
     * 因为它的任何一个变换矩阵对于一个给定的场景只需要生成一次，所以通常的做法是提前调用precalculate_all_transformation_matrices()来确定它们，并通过序列化保留它们。
     *
     */
    template <class Archive>
    void
    save_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * 从一个流中读取所有的变换矩阵，并为这个对象恢复它们。
     *
     */
    template <class Archive>
    void
    load_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * 测试两个系列扩展对象的相等。
     *
     */
    bool
    operator==(const Fourier<dim, spacedim> &fourier) const;

  private:
    /**
     * 注册的每个有限元在每个方向上的系数数
     * hp::FECollection.  。
     *
     */
    const std::vector<unsigned int> n_coefficients_per_direction;

    /**
     * hp::FECollection 为其计算变换矩阵。
     *
     */
    SmartPointer<const hp::FECollection<dim, spacedim>> fe_collection;

    /**
     * hp::QCollection 用于计算转换矩阵。
     *
     */
    const hp::QCollection<dim> q_collection;

    /**
     * 角度频率  $ 2 \pi {\bf k} $  .
     *
     */
    Table<dim, Tensor<1, dim>> k_vectors;

    /**
     * 每个FiniteElement的变换矩阵。
     *
     */
    std::vector<FullMatrix<CoefficientType>> fourier_transform_matrices;

    /**
     * 用于存储未滚动系数的辅助向量。
     *
     */
    std::vector<CoefficientType> unrolled_coefficients;

    /**
     * 应该使用FiniteElement的哪个分量来计算展开。
     *
     */
    const unsigned int component;
  };



  /**
   * 一个用于计算标量FE（或矢量FE的单个分量）场在参考元素上扩展为一系列Legendre函数的类。    Legendre函数是Legendre微分方程@f[
   *  \frac{d}{dx}\left([1-x^2] \frac{d}{dx} P_n(x)\right) +
   *  n[n+1] P_n(x) = 0
   * @f]的解，可以用Rodrigues公式@f[
   * P_n(x) = \frac{1}{2^n n!} \frac{d^n}{dx^n}[x^2-1]^n.
   * @f]来表达，这些多项式相对于区间 $ L^2 $ @f[
   * \int_{-1}^1 P_m(x) P_n(x) = \frac{2}{2n + 1} \delta_{mn}
   * @f]的内积是正交的，是完整的。  在 $ L^2 $ 上的 $ [0;1] $ 的正交多项式家族可以通过@f[
   * \widetilde P_m = \sqrt{2} P_m(2x-1).
   * @f]构建。 参考元素 $ [0;1] $ 上的任意标量FE场可以在完全正交基础上展开为@f[
   * u(x) = \sum_{m} c_m \widetilde P_{m}(x).
   * @f] ]从基的正交特性可以看出，@f[
   * c_m = \frac{2m+1}{2} \int_0^1 u(x) \widetilde P_m(x) dx .
   * @f]该类计算系数 $ c_{\bf k} $ 使用 $ dim $
   * -维Legendre多项式，使用张量积规则从 $ \widetilde P_m(x) $
   * 构建。
   *
   */
  template <int dim, int spacedim = dim>
  class Legendre : public Subscriptor
  {
  public:
    using CoefficientType = double;

    /**
     * 构造函数，初始化所有需要的数据结构。         @p
     * n_coefficients_per_direction 定义了每个方向的系数数， @p
     * fe_collection 是将用于展开的 hp::FECollection ， @p
     * q_collection 是用于整合 @p fe_collection.
     * 中每个FiniteElement的展开的 hp::QCollection
     * 由于Legendre展开只能在标量场上进行，这个类不会对矢量值的有限元操作，因此会抛出一个断言。然而，有限元场的每个分量可以分别被视为标量场，对其进行Legendre扩展也是可能的。为此，可选的参数
     * @p component 定义了每个有限元的哪个分量将被使用。
     * 默认值 @p component
     * 只适用于标量FEs，在这种情况下，它表示唯一的分量将被分解。对于矢量FE，必须明确提供一个非默认值。
     *
     */
    Legendre(const std::vector<unsigned int> &n_coefficients_per_direction,
             const hp::FECollection<dim, spacedim> &fe_collection,
             const hp::QCollection<dim> &           q_collection,
             const unsigned int component = numbers::invalid_unsigned_int);

    /**
     * 一个非默认的构造函数。 @p size_in_each_direction
     * 定义了每个方向的系数数， @p fe_collection 是
     * hp::FECollection ，将对其进行扩展， @p q_collection
     * 是用于整合 @p fe_collection.  @deprecated
     * 中每个FiniteElement的扩展的 hp::QCollection
     * 使用一个不同的构造函数代替。
     *
     */
    DEAL_II_DEPRECATED
    Legendre(const unsigned int n_coefficients_per_direction,
             const hp::FECollection<dim, spacedim> &fe_collection,
             const hp::QCollection<dim> &           q_collection);

    /**
     * 计算 @p legendre_coefficients 中由 @p local_dof_values
     * 给出的单元向量场的 @p cell_active_fe_index
     * 对应的FiniteElement。
     *
     */
    template <typename Number>
    void
    calculate(const dealii::Vector<Number> &local_dof_values,
              const unsigned int            cell_active_fe_index,
              Table<dim, CoefficientType> & legendre_coefficients);

    /**
     * 返回与 @p index 相关的有限元在所提供的 hp::FECollection.
     * 中每个坐标方向的系数数。
     *
     */
    unsigned int
    get_n_coefficients_per_direction(const unsigned int index) const;

    /**
     * 计算所有的转换矩阵，将有限元解转移到系列展开表示。
     * 这些矩阵将通过调用calculate()按需生成，并存储起来供重复使用。通常情况下，这个操作会消耗大量的工作量。有了这个函数，所有的矩阵将被提前计算。
     * 这样，我们就可以将其昂贵的生成与实际应用分开。
     *
     */
    void
    precalculate_all_transformation_matrices();

    /**
     * 将此对象的所有转换矩阵写入一个流中，以便进行序列化。
     * 由于它的任何一个变换矩阵对于一个给定的场景只需要生成一次，通常的做法是提前调用precalculate_all_transformation_matrices()来确定它们，并通过序列化来保留它们。
     *
     */
    template <class Archive>
    void
    save_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * 从一个流中读取所有的变换矩阵，并为这个对象恢复它们。
     *
     */
    template <class Archive>
    void
    load_transformation_matrices(Archive &ar, const unsigned int version);

    /**
     * 测试两个系列扩展对象的相等。
     *
     */
    bool
    operator==(const Legendre<dim, spacedim> &legendre) const;

  private:
    /**
     * 注册的每个有限元在每个方向上的系数数
     * hp::FECollection.  。
     *
     */
    const std::vector<unsigned int> n_coefficients_per_direction;

    /**
     * hp::FECollection 将计算其变换矩阵。
     *
     */
    SmartPointer<const hp::FECollection<dim, spacedim>> fe_collection;

    /**
     * hp::QCollection  用于计算转换矩阵。
     *
     */
    const hp::QCollection<dim> q_collection;

    /**
     * 每个FiniteElement的变换矩阵。
     *
     */
    std::vector<FullMatrix<CoefficientType>> legendre_transform_matrices;

    /**
     * 用于存储解卷系数的辅助向量。
     *
     */
    std::vector<CoefficientType> unrolled_coefficients;

    /**
     * 应该使用FiniteElement的哪个分量来计算展开。
     *
     */
    const unsigned int component;
  };



  /**
   * 计算由 @p predicate 定义的 @p coefficients 的子集的 @p norm
   * 为常数。返回一对谓词值的向量和计算出的子集规范的向量。
   * @p predicate 应该返回一对 <code>bool</code>
   * 和<code>无符号int</code>。前者是一个标志，表明是否应该在计算中使用给定的TableIndices，而后者是指数的解卷值，根据这个指数将形成系数的子集。
   * 只有那些大于 @p smallest_abs_coefficient.
   * 的系数才会被考虑。
   * @note  只有以下 @p norm_type
   * 的值可以实现，并且在这种情况下是有意义的：均值、L1_norm、L2_norm、Linfty_norm。均值准则只适用于实值系数。
   *
   */
  template <int dim, typename CoefficientType>
  std::pair<std::vector<unsigned int>, std::vector<double>>
  process_coefficients(const Table<dim, CoefficientType> &coefficients,
                       const std::function<std::pair<bool, unsigned int>(
                         const TableIndices<dim> &)> &    predicate,
                       const VectorTools::NormType        norm_type,
                       const double smallest_abs_coefficient = 1e-10);

  /**
   * $y = k \, x + b$ 的线性回归最小平方拟合。
   * 输入向量的大小应该等于和大于1。返回的对子将包含
   * $k$ （第一）和 $b$ （第二）。
   *
   */
  std::pair<double, double>
  linear_regression(const std::vector<double> &x, const std::vector<double> &y);

} // namespace FESeries

 /*@}*/ 



#ifndef DOXYGEN

// -------------------  inline and template functions ----------------

namespace internal
{
  namespace FESeriesImplementation
  {
    template <int dim, typename CoefficientType>
    void
    fill_map_index(
      const Table<dim, CoefficientType> &coefficients,
      const TableIndices<dim> &          ind,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<dim> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>> &pred_to_values)
    {
      const std::pair<bool, unsigned int> pred_pair = predicate(ind);
      // don't add a value if predicate is false
      if (pred_pair.first == false)
        return;

      const unsigned int     pred_value  = pred_pair.second;
      const CoefficientType &coeff_value = coefficients(ind);
      // If pred_value is not in the pred_to_values map, the element will be
      // created. Otherwise a reference to the existing element is returned.
      pred_to_values[pred_value].push_back(coeff_value);
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<1, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<1> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>> &     pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); i++)
        {
          const TableIndices<1> ind(i);
          fill_map_index(coefficients, ind, predicate, pred_to_values);
        }
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<2, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<2> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>> &     pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); i++)
        for (unsigned int j = 0; j < coefficients.size(1); j++)
          {
            const TableIndices<2> ind(i, j);
            fill_map_index(coefficients, ind, predicate, pred_to_values);
          }
    }



    template <typename CoefficientType>
    void
    fill_map(
      const Table<3, CoefficientType> &coefficients,
      const std::function<
        std::pair<bool, unsigned int>(const TableIndices<3> &)> &predicate,
      std::map<unsigned int, std::vector<CoefficientType>> &     pred_to_values)
    {
      for (unsigned int i = 0; i < coefficients.size(0); i++)
        for (unsigned int j = 0; j < coefficients.size(1); j++)
          for (unsigned int k = 0; k < coefficients.size(2); k++)
            {
              const TableIndices<3> ind(i, j, k);
              fill_map_index(coefficients, ind, predicate, pred_to_values);
            }
    }



    template <typename Number>
    double
    complex_mean_value(const Number &value)
    {
      return value;
    }



    template <typename Number>
    double
    complex_mean_value(const std::complex<Number> &value)
    {
      AssertThrow(false,
                  ExcMessage(
                    "FESeries::process_coefficients() can not be used with "
                    "complex-valued coefficients and VectorTools::mean norm."));
      return std::abs(value);
    }
  } // namespace FESeriesImplementation
} // namespace internal



template <int dim, typename CoefficientType>
std::pair<std::vector<unsigned int>, std::vector<double>>
FESeries::process_coefficients(
  const Table<dim, CoefficientType> &coefficients,
  const std::function<std::pair<bool, unsigned int>(const TableIndices<dim> &)>
    &                         predicate,
  const VectorTools::NormType norm_type,
  const double                smallest_abs_coefficient)
{
  Assert(smallest_abs_coefficient >= 0.,
         ExcMessage("smallest_abs_coefficient should be non-negative."));

  std::vector<unsigned int> predicate_values;
  std::vector<double>       norm_values;

  // first, parse all table elements into a map of predicate values and
  // coefficients. We could have stored (predicate values ->TableIndicies) map,
  // but its processing would have been much harder later on.
  std::map<unsigned int, std::vector<CoefficientType>> pred_to_values;
  internal::FESeriesImplementation::fill_map(coefficients,
                                             predicate,
                                             pred_to_values);

  // now go through the map and populate the @p norm_values based on @p norm:
  for (const auto &pred_to_value : pred_to_values)
    {
      Vector<CoefficientType> values(pred_to_value.second.cbegin(),
                                     pred_to_value.second.cend());

      double norm_value = 0;
      switch (norm_type)
        {
          case VectorTools::L2_norm:
            {
              norm_value = values.l2_norm();
              break;
            }
          case VectorTools::L1_norm:
            {
              norm_value = values.l1_norm();
              break;
            }
          case VectorTools::Linfty_norm:
            {
              norm_value = values.linfty_norm();
              break;
            }
          case VectorTools::mean:
            {
              norm_value = internal::FESeriesImplementation::complex_mean_value(
                values.mean_value());
              break;
            }
          default:
            AssertThrow(false, ExcNotImplemented());
            break;
        }

      // will use all non-zero coefficients
      if (std::abs(norm_value) > smallest_abs_coefficient)
        {
          predicate_values.push_back(pred_to_value.first);
          norm_values.push_back(norm_value);
        }
    }

  return std::make_pair(predicate_values, norm_values);
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Fourier<dim, spacedim>::save_transformation_matrices(
  Archive &ar,
  const unsigned int  /*version*/ )
{
  // Store information about those resources which have been used to generate
  // the transformation matrices.
  // mode vector
  ar &n_coefficients_per_direction;

  // finite element collection
  unsigned int size = fe_collection->size();
  ar &         size;
  for (unsigned int i = 0; i < size; ++i)
    ar &(*fe_collection)[i].get_name();

  // quadrature collection
  size = q_collection.size();
  ar &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &q_collection[i];

  // Store the actual transform matrices.
  ar &fourier_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Fourier<dim, spacedim>::load_transformation_matrices(
  Archive &ar,
  const unsigned int  /*version*/ )
{
  // Check whether the currently registered resources are compatible with
  // the transformation matrices to load.
  // mode vector
  std::vector<unsigned int> compare_coefficients;
  ar &                      compare_coefficients;
  Assert(compare_coefficients == n_coefficients_per_direction,
         ExcMessage("A different number of coefficients vector has been used "
                    "to generate the transformation matrices you are about "
                    "to load!"));

  // finite element collection
  unsigned int size;
  ar &         size;
  AssertDimension(size, fe_collection->size());
  std::string name;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &name;
      Assert(name.compare((*fe_collection)[i].get_name()) == 0,
             ExcMessage("A different FECollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // quadrature collection
  ar &size;
  AssertDimension(size, q_collection.size());
  Quadrature<dim> quadrature;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &quadrature;
      Assert(quadrature == q_collection[i],
             ExcMessage("A different QCollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // Restore the transform matrices since all prerequisites are fulfilled.
  ar &fourier_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Legendre<dim, spacedim>::save_transformation_matrices(
  Archive &ar,
  const unsigned int  /*version*/ )
{
  // Store information about those resources which have been used to generate
  // the transformation matrices.
  // mode vector
  ar &n_coefficients_per_direction;

  // finite element collection
  unsigned int size = fe_collection->size();
  ar &         size;
  for (unsigned int i = 0; i < size; ++i)
    ar &(*fe_collection)[i].get_name();

  // quadrature collection
  size = q_collection.size();
  ar &size;
  for (unsigned int i = 0; i < size; ++i)
    ar &q_collection[i];

  // Store the actual transform matrices.
  ar &legendre_transform_matrices;
}



template <int dim, int spacedim>
template <class Archive>
inline void
FESeries::Legendre<dim, spacedim>::load_transformation_matrices(
  Archive &ar,
  const unsigned int  /*version*/ )
{
  // Check whether the currently registered resources are compatible with
  // the transformation matrices to load.
  // mode vector
  std::vector<unsigned int> compare_coefficients;
  ar &                      compare_coefficients;
  Assert(compare_coefficients == n_coefficients_per_direction,
         ExcMessage("A different number of coefficients vector has been used "
                    "to generate the transformation matrices you are about "
                    "to load!"));

  // finite element collection
  unsigned int size;
  ar &         size;
  AssertDimension(size, fe_collection->size());
  std::string name;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &name;
      Assert(name.compare((*fe_collection)[i].get_name()) == 0,
             ExcMessage("A different FECollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // quadrature collection
  ar &size;
  AssertDimension(size, q_collection.size());
  Quadrature<dim> quadrature;
  for (unsigned int i = 0; i < size; ++i)
    {
      ar &quadrature;
      Assert(quadrature == q_collection[i],
             ExcMessage("A different QCollection has been used to generate "
                        "the transformation matrices you are about to load!"));
    }

  // Restore the transform matrices since all prerequisites are fulfilled.
  ar &legendre_transform_matrices;
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_fe_series_h


