//include/deal.II-translator/matrix_free/tensor_product_kernels_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2021 by the deal.II authors
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


#ifndef dealii_matrix_free_tensor_product_kernels_h
#define dealii_matrix_free_tensor_product_kernels_h

#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/utilities.h>


DEAL_II_NAMESPACE_OPEN



namespace internal
{
  /**
   * 在这个命名空间中，实现了评估张量产品的评估器例程。
   *
   */
  enum EvaluatorVariant
  {
    /**
     * 不要使用比有限元的张量积结构更多的东西。
     *
     */
    evaluate_general,
    /**
     * 通过利用有限元的对称性进行评估：即通过利用形状函数和正交点的对称性跳过一些计算。
     *
     */
    evaluate_symmetric,
    /**
     * 利用对称性将算子分别应用于输入矢量的偶数和奇数部分：更多信息见EvaluatorTensorProduct专业化的文档。
     *
     */
    evaluate_evenodd,
    /**
     * 在Legendre和类似的多项式空间中使用对称性，其中偶数的形状函数围绕正交点的中心对称（考虑偶数多项式度数），奇数的形状函数围绕正交点的中心反对称（考虑奇数多项式度数）。这允许使用类似于偶数技术的策略，但不需要单独的系数数组。更多信息请参见EvaluatorTensorProduct专业化的文档。
     *
     */
    evaluate_symmetric_hierarchical
  };



  /**
   * 决定哪个数量应该通过张量积核计算。
   *
   */
  enum class EvaluatorQuantity
  {
    /**
     * 通过形状函数进行评估/积分。
     *
     */
    value,
    /**
     * 通过形状函数的梯度进行评估/积分。
     *
     */
    gradient,
    /**
     * 通过形状函数的Hessians进行评估/积分。
     *
     */
    hessian
  };



  /**
   * 通用的评估器框架，使用张量积形式对一般维度的给定形状数据进行估值。根据矩阵条目中的特定布局，这对应于通常的矩阵-矩阵乘积或包括一些对称性的矩阵-矩阵乘积。
   * @tparam  variant 用于创建模板特化的评估变量  @tparam  dim
   * 函数的尺寸  @tparam  n_rows
   * 变换矩阵中的行数，相当于通常张量收缩设置中的1d形状函数的数量
   * @tparam  ] n_columns
   * 变换矩阵中的列数，相当于通常张量收缩设置中的1d形状函数的数量
   * @tparam  Number 输入和输出数组的抽象数字类型  @tparam
   * Number2
   * 系数数组的抽象数字类型（默认为与输入/输出数组的类型相同）；必须用Number实现operator*才能有效
   *
   */
  template <EvaluatorVariant variant,
            int              dim,
            int              n_rows,
            int              n_columns,
            typename Number,
            typename Number2 = Number>
  struct EvaluatorTensorProduct
  {};



  /**
   * 使用基函数的张量积形式的任意维度的形状函数的内部评估器。
   * @tparam  dim 应用该类的空间维度  @tparam  n_rows
   * 变换矩阵中的行数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  n_columns
   * 变换矩阵中的列数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  ] Number 用于输入和输出数组的抽象数字类型
   * @tparam  Number2
   * 用于系数数组的抽象数字类型（默认为与输入/输出数组相同的类型）；必须用Number实现操作符*，并产生Number作为输出，才能成为有效类型。
   *
   */
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  struct EvaluatorTensorProduct<evaluate_general,
                                dim,
                                n_rows,
                                n_columns,
                                Number,
                                Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      Utilities::pow(n_rows, dim);
    static constexpr unsigned int n_columns_of_product =
      Utilities::pow(n_columns, dim);

    /**
     * 空的构造函数。什么都不做。在使用'value'和相关方法时要小心，因为它们需要用其他指针来填充。
     *
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {}

    /**
     * 构造函数，从ShapeInfo中获取数据
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            dummy1 = 0,
                           const unsigned int            dummy2 = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
    {
      // We can enter this function either for the apply() path that has
      // n_rows * n_columns entries or for the apply_face() path that only has
      // n_rows * 3 entries in the array. Since we cannot decide about the use
      // we must allow for both here.
      Assert(shape_values.size() == 0 ||
               shape_values.size() == n_rows * n_columns ||
               shape_values.size() == 3 * n_rows,
             ExcDimensionMismatch(shape_values.size(), n_rows * n_columns));
      Assert(shape_gradients.size() == 0 ||
               shape_gradients.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows * n_columns));
      Assert(shape_hessians.size() == 0 ||
               shape_hessians.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows * n_columns));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values(const Number in[], Number out[]) const
    {
      apply<direction, contract_over_rows, add>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients(const Number in[], Number out[]) const
    {
      apply<direction, contract_over_rows, add>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number in[], Number out[]) const
    {
      apply<direction, contract_over_rows, add>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_hessians, in, out);
    }

    /**
     * 这个函数沿着输入数组中张量数据的给定 @p direction
     * 应用张量积核，对应于一维条纹的乘法。这个函数允许
     * @p in 和 @p out 数组在n_rows ==
     * n_columns的情况下别名，也就是说，在 @p in 和 @p out
     * 指向相同地址的地方执行收缩是安全的。对于n_rows !=
     * n_columns的情况，一般来说输出是不正确的。
     * @tparam  被评估的方向  @tparam  contract_over_rows
     * 如果为真，张量收缩对给定的 @p shape_data
     * 数组中的行求和，否则对列求和  @tparam  ] add
     * 如果为真，结果将被添加到输出向量中，否则计算值将覆盖输出中的内容
     * @tparam  one_line
     * 如果为真，内核只沿着二维张量中的单一一维条纹应用，而不是像
     * @p false 情况下的全部n_rows^dim点。          @param  shape_data
     * 具有 @p n_rows 行和 @p n_columns
     * 列的变换矩阵，以行为主的格式存储  @param  in
     * 指向输入数据向量的起点  @param  out
     * 指向输出数据向量的起点
     *
     */
    template <int  direction,
              bool contract_over_rows,
              bool add,
              bool one_line = false>
    static void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number *                  in,
          Number *                        out);

    /**
     * 这个函数应用张量积操作，从单元格值中产生面值。与apply方法相反，这个方法假设与面的正交方向每个方向有n_rows自由度，而不是那些比当前应用的方向低的n_columns。换句话说，apply_face()必须在调用面内的任何插值之前被调用。
     * @tparam  face_direction 法向量的方向（0=x，1=y，等等）
     * @tparam  contract_onto_face
     * 如果为真，输入向量的大小为n_rows^dim，插值将被执行到n_rows^(dim-1)点。这是
     * FEFaceEvaluation::evaluate()
     * 调用中的一个典型场景。如果是假的，来自n_rows^(dim-1)点的数据被扩展到高维数据阵列的n_rows^dim点。在contract_onto_face==false的情况下，导数被加在一起
     * @tparam  add
     * 如果是true，结果被加到输出向量中，否则计算出的值会覆盖输出中的内容
     * @tparam  max_derivative
     * 设置应该被计算的导数的数量。0表示只有数值，1表示数值和第一导数，2表示第二导数。注意，所有的导数都要访问传递给类的构造函数的
     * @p shape_values 中的数据  @tparam  lex_faces
     * 设置面的评估点应该如何排序：按词典排序或按右手系统号排序（在三维中对方位1进行特殊处理）。默认情况下，右键系统号被启用，这只适用于3以内的尺寸。
     * @param  输入数据向量的地址  @param
     * 输出数据向量的地址
     *
     */
    template <int  face_direction,
              bool contract_onto_face,
              bool add,
              int  max_derivative,
              bool lex_faces = false>
    void
    apply_face(const Number *DEAL_II_RESTRICT in,
               Number *DEAL_II_RESTRICT out) const;

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int direction, bool contract_over_rows, bool add, bool one_line>
  inline void
  EvaluatorTensorProduct<evaluate_general,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::apply(const Number2 *DEAL_II_RESTRICT
                                                        shape_data,
                                         const Number * in,
                                         Number *       out)
  {
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage(
             "The given array shape_data must not be the null pointer!"));
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange(direction, dim);
    constexpr int mm = contract_over_rows ? n_rows : n_columns,
                  nn = contract_over_rows ? n_columns : n_rows;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            Number x[mm];
            for (int i = 0; i < mm; ++i)
              x[i] = in[stride * i];
            for (int col = 0; col < nn; ++col)
              {
                Number2 val0;
                if (contract_over_rows == true)
                  val0 = shape_data[col];
                else
                  val0 = shape_data[col * n_columns];
                Number res0 = val0 * x[0];
                for (int i = 1; i < mm; ++i)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_data[i * n_columns + col];
                    else
                      val0 = shape_data[col * n_columns + i];
                    res0 += val0 * x[i];
                  }
                if (add == false)
                  out[stride * col] = res0;
                else
                  out[stride * col] += res0;
              }

            if (one_line == false)
              {
                ++in;
                ++out;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm - 1);
            out += stride * (nn - 1);
          }
      }
  }



  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int  face_direction,
            bool contract_onto_face,
            bool add,
            int  max_derivative,
            bool lex_faces>
  inline void
  EvaluatorTensorProduct<evaluate_general,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::apply_face(const Number *DEAL_II_RESTRICT in,
                                              Number *DEAL_II_RESTRICT
                                                      out) const
  {
    Assert(dim > 0 && (lex_faces || dim < 4),
           ExcMessage("Only dim=1,2,3 supported"));
    static_assert(max_derivative >= 0 && max_derivative < 3,
                  "Only derivative orders 0-2 implemented");
    Assert(shape_values != nullptr,
           ExcMessage(
             "The given array shape_values must not be the null pointer."));

    constexpr int n_blocks1 =
      lex_faces ? dealii::Utilities::pow<unsigned int>(n_rows, face_direction) :
                  (dim > 1 ? n_rows : 1);
    constexpr int n_blocks2 =
      lex_faces ? dealii::Utilities::pow<unsigned int>(
                    n_rows, std::max(dim - face_direction - 1, 0)) :
                  (dim > 2 ? n_rows : 1);

    AssertIndexRange(face_direction, dim);
    constexpr int stride     = Utilities::pow(n_rows, face_direction);
    constexpr int out_stride = Utilities::pow(n_rows, dim - 1);
    const Number *DEAL_II_RESTRICT shape_values = this->shape_values;

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            if (contract_onto_face == true)
              {
                Number res0 = shape_values[0] * in[0];
                Number res1, res2;
                if (max_derivative > 0)
                  res1 = shape_values[n_rows] * in[0];
                if (max_derivative > 1)
                  res2 = shape_values[2 * n_rows] * in[0];
                for (int ind = 1; ind < n_rows; ++ind)
                  {
                    res0 += shape_values[ind] * in[stride * ind];
                    if (max_derivative > 0)
                      res1 += shape_values[ind + n_rows] * in[stride * ind];
                    if (max_derivative > 1)
                      res2 += shape_values[ind + 2 * n_rows] * in[stride * ind];
                  }
                if (add == false)
                  {
                    out[0] = res0;
                    if (max_derivative > 0)
                      out[out_stride] = res1;
                    if (max_derivative > 1)
                      out[2 * out_stride] = res2;
                  }
                else
                  {
                    out[0] += res0;
                    if (max_derivative > 0)
                      out[out_stride] += res1;
                    if (max_derivative > 1)
                      out[2 * out_stride] += res2;
                  }
              }
            else
              {
                for (int col = 0; col < n_rows; ++col)
                  {
                    if (add == false)
                      out[col * stride] = shape_values[col] * in[0];
                    else
                      out[col * stride] += shape_values[col] * in[0];
                    if (max_derivative > 0)
                      out[col * stride] +=
                        shape_values[col + n_rows] * in[out_stride];
                    if (max_derivative > 1)
                      out[col * stride] +=
                        shape_values[col + 2 * n_rows] * in[2 * out_stride];
                  }
              }

            if (lex_faces)
              {
                ++out;
                ++in;
              }
            else
              // increment: in regular case, just go to the next point in
              // x-direction. If we are at the end of one chunk in x-dir, need
              // to jump over to the next layer in z-direction
              switch (face_direction)
                {
                  case 0:
                    in += contract_onto_face ? n_rows : 1;
                    out += contract_onto_face ? 1 : n_rows;
                    break;
                  case 1:
                    ++in;
                    ++out;
                    // faces 2 and 3 in 3D use local coordinate system zx, which
                    // is the other way around compared to the tensor
                    // product. Need to take that into account.
                    if (dim == 3)
                      {
                        if (contract_onto_face)
                          out += n_rows - 1;
                        else
                          in += n_rows - 1;
                      }
                    break;
                  case 2:
                    ++in;
                    ++out;
                    break;
                  default:
                    Assert(false, ExcNotImplemented());
                }
          }
        if (lex_faces)
          {
            if (contract_onto_face)
              in += (dealii::Utilities::pow(n_rows, face_direction + 1) -
                     n_blocks1);
            else
              out += (dealii::Utilities::pow(n_rows, face_direction + 1) -
                      n_blocks1);
          }
        else if (face_direction == 1 && dim == 3)
          {
            // adjust for local coordinate system zx
            if (contract_onto_face)
              {
                in += n_rows * (n_rows - 1);
                out -= n_rows * n_rows - 1;
              }
            else
              {
                out += n_rows * (n_rows - 1);
                in -= n_rows * n_rows - 1;
              }
          }
      }
  }



  /**
   * 形状函数的内部评估器，使用基函数的张量积形式。与其他模板类相同，但没有利用模板参数和变量循环边界来代替。
   * @tparam  dim 应用该类的空间维度  @tparam  Number
   * 用于输入和输出数组的抽象数字类型  @tparam  Number2
   * 用于系数数组的抽象数字类型（默认为与输入/输出数组的类型相同）；必须用Number实现操作符*，并产生Number作为输出，才能成为有效类型
   *
   */
  template <int dim, typename Number, typename Number2>
  struct EvaluatorTensorProduct<evaluate_general, dim, 0, 0, Number, Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      numbers::invalid_unsigned_int;
    static constexpr unsigned int n_columns_of_product =
      numbers::invalid_unsigned_int;

    /**
     * 空的构造函数。什么都不做。在使用'values'和相关方法时要小心，因为它们需要用其他构造函数来填充。
     *
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
      , n_rows(numbers::invalid_unsigned_int)
      , n_columns(numbers::invalid_unsigned_int)
    {}

    /**
     * 构造函数，从ShapeInfo中获取数据
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            n_rows,
                           const unsigned int            n_columns)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
      , n_rows(n_rows)
      , n_columns(n_columns)
    {
      // We can enter this function either for the apply() path that has
      // n_rows * n_columns entries or for the apply_face() path that only has
      // n_rows * 3 entries in the array. Since we cannot decide about the use
      // we must allow for both here.
      Assert(shape_values.size() == 0 ||
               shape_values.size() == n_rows * n_columns ||
               shape_values.size() == n_rows * 3,
             ExcDimensionMismatch(shape_values.size(), n_rows * n_columns));
      Assert(shape_gradients.size() == 0 ||
               shape_gradients.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows * n_columns));
      Assert(shape_hessians.size() == 0 ||
               shape_hessians.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows * n_columns));
    }

    /**
     * 构造函数，从ShapeInfo中获取数据
     *
     */
    EvaluatorTensorProduct(const Number2 *    shape_values,
                           const Number2 *    shape_gradients,
                           const Number2 *    shape_hessians,
                           const unsigned int n_rows,
                           const unsigned int n_columns)
      : shape_values(shape_values)
      , shape_gradients(shape_gradients)
      , shape_hessians(shape_hessians)
      , n_rows(n_rows)
      , n_columns(n_columns)
    {}

    template <int direction, bool contract_over_rows, bool add>
    void
    values(const Number *in, Number *out) const
    {
      apply<direction, contract_over_rows, add>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients(const Number *in, Number *out) const
    {
      apply<direction, contract_over_rows, add>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number *in, Number *out) const
    {
      apply<direction, contract_over_rows, add>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, true>(shape_hessians, in, out);
    }

    template <int  direction,
              bool contract_over_rows,
              bool add,
              bool one_line = false>
    void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number *                  in,
          Number *                        out) const;

    template <int  face_direction,
              bool contract_onto_face,
              bool add,
              int  max_derivative,
              bool lex_faces = false>
    void
    apply_face(const Number *DEAL_II_RESTRICT in,
               Number *DEAL_II_RESTRICT out) const;

    const Number2 *    shape_values;
    const Number2 *    shape_gradients;
    const Number2 *    shape_hessians;
    const unsigned int n_rows;
    const unsigned int n_columns;
  };



  template <int dim, typename Number, typename Number2>
  template <int direction, bool contract_over_rows, bool add, bool one_line>
  inline void
  EvaluatorTensorProduct<evaluate_general, dim, 0, 0, Number, Number2>::apply(
    const Number2 *DEAL_II_RESTRICT shape_data,
    const Number *                  in,
    Number *                        out) const
  {
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    Assert(shape_data != nullptr,
           ExcMessage(
             "The given array shape_data must not be the null pointer!"));
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));
    AssertIndexRange(direction, dim);
    const int mm = contract_over_rows ? n_rows : n_columns,
              nn = contract_over_rows ? n_columns : n_rows;

    const int stride =
      direction == 0 ? 1 : Utilities::fixed_power<direction>(n_columns);
    const int n_blocks1 = one_line ? 1 : stride;
    const int n_blocks2 = direction >= dim - 1 ?
                            1 :
                            Utilities::fixed_power<dim - direction - 1>(n_rows);
    Assert(n_rows <= 128, ExcNotImplemented());

    // specialization for n_rows = 2 that manually unrolls the innermost loop
    // to make the operation perform better (not completely as good as the
    // templated one, but much better than the generic version down below,
    // because the loop over col can be more effectively unrolled by the
    // compiler)
    if (contract_over_rows && n_rows == 2)
      {
        const Number2 *shape_data_1 = shape_data + n_columns;
        for (int i2 = 0; i2 < n_blocks2; ++i2)
          {
            for (int i1 = 0; i1 < n_blocks1; ++i1)
              {
                const Number x0 = in[0], x1 = in[stride];
                for (int col = 0; col < nn; ++col)
                  {
                    const Number result =
                      shape_data[col] * x0 + shape_data_1[col] * x1;
                    if (add == false)
                      out[stride * col] = result;
                    else
                      out[stride * col] += result;
                  }

                if (one_line == false)
                  {
                    ++in;
                    ++out;
                  }
              }
            if (one_line == false)
              {
                in += stride * (mm - 1);
                out += stride * (nn - 1);
              }
          }
      }
    // specialization for n = 3
    else if (contract_over_rows && n_rows == 3)
      {
        const Number2 *shape_data_1 = shape_data + n_columns;
        const Number2 *shape_data_2 = shape_data + 2 * n_columns;
        for (int i2 = 0; i2 < n_blocks2; ++i2)
          {
            for (int i1 = 0; i1 < n_blocks1; ++i1)
              {
                const Number x0 = in[0], x1 = in[stride], x2 = in[2 * stride];
                for (int col = 0; col < nn; ++col)
                  {
                    const Number result = shape_data[col] * x0 +
                                          shape_data_1[col] * x1 +
                                          shape_data_2[col] * x2;
                    if (add == false)
                      out[stride * col] = result;
                    else
                      out[stride * col] += result;
                  }

                if (one_line == false)
                  {
                    ++in;
                    ++out;
                  }
              }
            if (one_line == false)
              {
                in += stride * (mm - 1);
                out += stride * (nn - 1);
              }
          }
      }
    // general loop for all other cases
    else
      for (int i2 = 0; i2 < n_blocks2; ++i2)
        {
          for (int i1 = 0; i1 < n_blocks1; ++i1)
            {
              Number x[129];
              for (int i = 0; i < mm; ++i)
                x[i] = in[stride * i];
              for (int col = 0; col < nn; ++col)
                {
                  Number2 val0;
                  if (contract_over_rows == true)
                    val0 = shape_data[col];
                  else
                    val0 = shape_data[col * n_columns];
                  Number res0 = val0 * x[0];
                  for (int i = 1; i < mm; ++i)
                    {
                      if (contract_over_rows == true)
                        val0 = shape_data[i * n_columns + col];
                      else
                        val0 = shape_data[col * n_columns + i];
                      res0 += val0 * x[i];
                    }
                  if (add == false)
                    out[stride * col] = res0;
                  else
                    out[stride * col] += res0;
                }

              if (one_line == false)
                {
                  ++in;
                  ++out;
                }
            }
          if (one_line == false)
            {
              in += stride * (mm - 1);
              out += stride * (nn - 1);
            }
        }
  }



  template <int dim, typename Number, typename Number2>
  template <int  face_direction,
            bool contract_onto_face,
            bool add,
            int  max_derivative,
            bool lex_faces>
  inline void
  EvaluatorTensorProduct<evaluate_general, dim, 0, 0, Number, Number2>::
    apply_face(const Number *DEAL_II_RESTRICT in,
               Number *DEAL_II_RESTRICT out) const
  {
    static_assert(lex_faces == false, "Not implemented yet.");

    Assert(shape_values != nullptr,
           ExcMessage(
             "The given array shape_data must not be the null pointer!"));
    static_assert(dim > 0 && dim < 4, "Only dim=1,2,3 supported");
    const int n_blocks1 = dim > 1 ? n_rows : 1;
    const int n_blocks2 = dim > 2 ? n_rows : 1;

    AssertIndexRange(face_direction, dim);
    const int stride =
      face_direction > 0 ? Utilities::fixed_power<face_direction>(n_rows) : 1;
    const int out_stride =
      dim > 1 ? Utilities::fixed_power<dim - 1>(n_rows) : 1;

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            if (contract_onto_face == true)
              {
                Number res0 = shape_values[0] * in[0];
                Number res1, res2;
                if (max_derivative > 0)
                  res1 = shape_values[n_rows] * in[0];
                if (max_derivative > 1)
                  res2 = shape_values[2 * n_rows] * in[0];
                for (unsigned int ind = 1; ind < n_rows; ++ind)
                  {
                    res0 += shape_values[ind] * in[stride * ind];
                    if (max_derivative > 0)
                      res1 += shape_values[ind + n_rows] * in[stride * ind];
                    if (max_derivative > 1)
                      res2 += shape_values[ind + 2 * n_rows] * in[stride * ind];
                  }
                if (add == false)
                  {
                    out[0] = res0;
                    if (max_derivative > 0)
                      out[out_stride] = res1;
                    if (max_derivative > 1)
                      out[2 * out_stride] = res2;
                  }
                else
                  {
                    out[0] += res0;
                    if (max_derivative > 0)
                      out[out_stride] += res1;
                    if (max_derivative > 1)
                      out[2 * out_stride] += res2;
                  }
              }
            else
              {
                for (unsigned int col = 0; col < n_rows; ++col)
                  {
                    if (add == false)
                      out[col * stride] = shape_values[col] * in[0];
                    else
                      out[col * stride] += shape_values[col] * in[0];
                    if (max_derivative > 0)
                      out[col * stride] +=
                        shape_values[col + n_rows] * in[out_stride];
                    if (max_derivative > 1)
                      out[col * stride] +=
                        shape_values[col + 2 * n_rows] * in[2 * out_stride];
                  }
              }

            // increment: in regular case, just go to the next point in
            // x-direction. If we are at the end of one chunk in x-dir, need
            // to jump over to the next layer in z-direction
            switch (face_direction)
              {
                case 0:
                  in += contract_onto_face ? n_rows : 1;
                  out += contract_onto_face ? 1 : n_rows;
                  break;
                case 1:
                  ++in;
                  ++out;
                  // faces 2 and 3 in 3D use local coordinate system zx, which
                  // is the other way around compared to the tensor
                  // product. Need to take that into account.
                  if (dim == 3)
                    {
                      if (contract_onto_face)
                        out += n_rows - 1;
                      else
                        in += n_rows - 1;
                    }
                  break;
                case 2:
                  ++in;
                  ++out;
                  break;
                default:
                  Assert(false, ExcNotImplemented());
              }
          }
        if (face_direction == 1 && dim == 3)
          {
            // adjust for local coordinate system zx
            if (contract_onto_face)
              {
                in += n_rows * (n_rows - 1);
                out -= n_rows * n_rows - 1;
              }
            else
              {
                out += n_rows * (n_rows - 1);
                in -= n_rows * n_rows - 1;
              }
          }
      }
  }



  /**
   * 使用基函数的张量积形式的1d-3d形状函数的内部评估器。该类专门针对
   * "对称
   * "有限元的基于张量积的元素的一般应用，即当形状函数关于0.5的对称性和正交点也是如此。
   * @tparam  dim 应用该类的空间维度  @tparam  n_rows
   * 变换矩阵中的行数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  n_columns
   * 变换矩阵中的列数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  ] Number 用于输入和输出数组的抽象数字类型
   * @tparam  Number2
   * 用于系数数组的抽象数字类型（默认为与输入/输出数组相同的类型）；必须用Number实现操作符*，并产生Number作为输出，才能成为有效类型。
   *
   */
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  struct EvaluatorTensorProduct<evaluate_symmetric,
                                dim,
                                n_rows,
                                n_columns,
                                Number,
                                Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      Utilities::pow(n_rows, dim);
    static constexpr unsigned int n_columns_of_product =
      Utilities::pow(n_columns, dim);

    /**
     * 构造函数，从ShapeInfo获取数据。
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            dummy1 = 0,
                           const unsigned int            dummy2 = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
    {
      Assert(shape_values.size() == 0 ||
               shape_values.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_values.size(), n_rows * n_columns));
      Assert(shape_gradients.size() == 0 ||
               shape_gradients.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_gradients.size(), n_rows * n_columns));
      Assert(shape_hessians.size() == 0 ||
               shape_hessians.size() == n_rows * n_columns,
             ExcDimensionMismatch(shape_hessians.size(), n_rows * n_columns));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values(const Number in[], Number out[]) const;

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients(const Number in[], Number out[]) const;

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number in[], Number out[]) const;

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  // In this case, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [ 0.687  0 -0.087 ]
  //        [ 0.4    1  0.4   ]
  //        [-0.087  0  0.687 ]
  // Q3 --> [ 0.66   0.003  0.002  0.049 ]
  //        [ 0.521  1.005 -0.01  -0.230 ]
  //        [-0.230 -0.01   1.005  0.521 ]
  //        [ 0.049  0.002  0.003  0.66  ]
  // Q4 --> [ 0.658  0.022  0 -0.007 -0.032 ]
  //        [ 0.608  1.059  0  0.039  0.176 ]
  //        [-0.409 -0.113  1 -0.113 -0.409 ]
  //        [ 0.176  0.039  0  1.059  0.608 ]
  //        [-0.032 -0.007  0  0.022  0.658 ]
  //
  // In these matrices, we want to use avoid computations involving zeros and
  // ones and in addition use the symmetry in entries to reduce the number of
  // read operations.
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline void
  EvaluatorTensorProduct<evaluate_symmetric,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::values(const Number in[], Number out[]) const
  {
    Assert(shape_values != nullptr, ExcNotInitialized());
    AssertIndexRange(direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns,
                  nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            for (int col = 0; col < n_cols; ++col)
              {
                Number2 val0, val1;
                Number  in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_values[col];
                    val1 = shape_values[nn - 1 - col];
                  }
                else
                  {
                    val0 = shape_values[col * n_columns];
                    val1 = shape_values[(col + 1) * n_columns - 1];
                  }
                if (mid > 0)
                  {
                    in0  = in[0];
                    in1  = in[stride * (mm - 1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_values[ind * n_columns + col];
                            val1 = shape_values[ind * n_columns + nn - 1 - col];
                          }
                        else
                          {
                            val0 = shape_values[col * n_columns + ind];
                            val1 =
                              shape_values[(col + 1) * n_columns - 1 - ind];
                          }
                        in0 = in[stride * ind];
                        in1 = in[stride * (mm - 1 - ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (contract_over_rows == true)
                  {
                    if (mm % 2 == 1)
                      {
                        val0 = shape_values[mid * n_columns + col];
                        in1  = val0 * in[stride * mid];
                        res0 += in1;
                        res1 += in1;
                      }
                  }
                else
                  {
                    if (mm % 2 == 1 && nn % 2 == 0)
                      {
                        val0 = shape_values[col * n_columns + mid];
                        in1  = val0 * in[stride * mid];
                        res0 += in1;
                        res1 += in1;
                      }
                  }
                if (add == false)
                  {
                    out[stride * col]            = res0;
                    out[stride * (nn - 1 - col)] = res1;
                  }
                else
                  {
                    out[stride * col] += res0;
                    out[stride * (nn - 1 - col)] += res1;
                  }
              }
            if (contract_over_rows == true && nn % 2 == 1 && mm % 2 == 1)
              {
                if (add == false)
                  out[stride * n_cols] = in[stride * mid];
                else
                  out[stride * n_cols] += in[stride * mid];
              }
            else if (contract_over_rows == true && nn % 2 == 1)
              {
                Number  res0;
                Number2 val0 = shape_values[n_cols];
                if (mid > 0)
                  {
                    res0 = val0 * (in[0] + in[stride * (mm - 1)]);
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        val0 = shape_values[ind * n_columns + n_cols];
                        res0 += val0 * (in[stride * ind] +
                                        in[stride * (mm - 1 - ind)]);
                      }
                  }
                else
                  res0 = Number();
                if (add == false)
                  out[stride * n_cols] = res0;
                else
                  out[stride * n_cols] += res0;
              }
            else if (contract_over_rows == false && nn % 2 == 1)
              {
                Number res0;
                if (mid > 0)
                  {
                    Number2 val0 = shape_values[n_cols * n_columns];
                    res0         = val0 * (in[0] + in[stride * (mm - 1)]);
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        val0       = shape_values[n_cols * n_columns + ind];
                        Number in1 = val0 * (in[stride * ind] +
                                             in[stride * (mm - 1 - ind)]);
                        res0 += in1;
                      }
                    if (mm % 2)
                      res0 += in[stride * mid];
                  }
                else
                  res0 = in[0];
                if (add == false)
                  out[stride * n_cols] = res0;
                else
                  out[stride * n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride * (mm - 1);
        out += stride * (nn - 1);
      }
  }



  // For the specialized loop used for the gradient computation in
  // here, the 1D shape values read (sorted lexicographically, rows
  // run over 1D dofs, columns over quadrature points):
  // Q2 --> [-2.549 -1  0.549 ]
  //        [ 3.098  0 -3.098 ]
  //        [-0.549  1  2.549 ]
  // Q3 --> [-4.315 -1.03  0.5  -0.44  ]
  //        [ 6.07  -1.44 -2.97  2.196 ]
  //        [-2.196  2.97  1.44 -6.07  ]
  //        [ 0.44  -0.5   1.03  4.315 ]
  // Q4 --> [-6.316 -1.3    0.333 -0.353  0.413 ]
  //        [10.111 -2.76  -2.667  2.066 -2.306 ]
  //        [-5.688  5.773  0     -5.773  5.688 ]
  //        [ 2.306 -2.066  2.667  2.76 -10.111 ]
  //        [-0.413  0.353 -0.333 -0.353  0.413 ]
  //
  // In these matrices, we want to use avoid computations involving
  // zeros and ones and in addition use the symmetry in entries to
  // reduce the number of read operations.
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline void
  EvaluatorTensorProduct<evaluate_symmetric,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::gradients(const Number in[],
                                             Number       out[]) const
  {
    Assert(shape_gradients != nullptr, ExcNotInitialized());
    AssertIndexRange(direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns,
                  nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            for (int col = 0; col < n_cols; ++col)
              {
                Number2 val0, val1;
                Number  in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_gradients[col];
                    val1 = shape_gradients[nn - 1 - col];
                  }
                else
                  {
                    val0 = shape_gradients[col * n_columns];
                    val1 = shape_gradients[(nn - col - 1) * n_columns];
                  }
                if (mid > 0)
                  {
                    in0  = in[0];
                    in1  = in[stride * (mm - 1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 -= val1 * in1;
                    res1 -= val0 * in1;
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_gradients[ind * n_columns + col];
                            val1 =
                              shape_gradients[ind * n_columns + nn - 1 - col];
                          }
                        else
                          {
                            val0 = shape_gradients[col * n_columns + ind];
                            val1 =
                              shape_gradients[(nn - col - 1) * n_columns + ind];
                          }
                        in0 = in[stride * ind];
                        in1 = in[stride * (mm - 1 - ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 -= val1 * in1;
                        res1 -= val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_gradients[mid * n_columns + col];
                    else
                      val0 = shape_gradients[col * n_columns + mid];
                    in1 = val0 * in[stride * mid];
                    res0 += in1;
                    res1 -= in1;
                  }
                if (add == false)
                  {
                    out[stride * col]            = res0;
                    out[stride * (nn - 1 - col)] = res1;
                  }
                else
                  {
                    out[stride * col] += res0;
                    out[stride * (nn - 1 - col)] += res1;
                  }
              }
            if (nn % 2 == 1)
              {
                Number2 val0;
                Number  res0;
                if (contract_over_rows == true)
                  val0 = shape_gradients[n_cols];
                else
                  val0 = shape_gradients[n_cols * n_columns];
                res0 = val0 * (in[0] - in[stride * (mm - 1)]);
                for (int ind = 1; ind < mid; ++ind)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_gradients[ind * n_columns + n_cols];
                    else
                      val0 = shape_gradients[n_cols * n_columns + ind];
                    Number in1 =
                      val0 * (in[stride * ind] - in[stride * (mm - 1 - ind)]);
                    res0 += in1;
                  }
                if (add == false)
                  out[stride * n_cols] = res0;
                else
                  out[stride * n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride * (mm - 1);
        out += stride * (nn - 1);
      }
  }



  // evaluates the given shape data in 1d-3d using the tensor product
  // form assuming the symmetries of unit cell shape hessians for
  // finite elements in FEEvaluation
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int direction, bool contract_over_rows, bool add>
  inline void
  EvaluatorTensorProduct<evaluate_symmetric,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::hessians(const Number in[],
                                            Number       out[]) const
  {
    Assert(shape_hessians != nullptr, ExcNotInitialized());
    AssertIndexRange(direction, dim);
    constexpr int mm     = contract_over_rows ? n_rows : n_columns;
    constexpr int nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            for (int col = 0; col < n_cols; ++col)
              {
                Number2 val0, val1;
                Number  in0, in1, res0, res1;
                if (contract_over_rows == true)
                  {
                    val0 = shape_hessians[col];
                    val1 = shape_hessians[nn - 1 - col];
                  }
                else
                  {
                    val0 = shape_hessians[col * n_columns];
                    val1 = shape_hessians[(col + 1) * n_columns - 1];
                  }
                if (mid > 0)
                  {
                    in0  = in[0];
                    in1  = in[stride * (mm - 1)];
                    res0 = val0 * in0;
                    res1 = val1 * in0;
                    res0 += val1 * in1;
                    res1 += val0 * in1;
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            val0 = shape_hessians[ind * n_columns + col];
                            val1 =
                              shape_hessians[ind * n_columns + nn - 1 - col];
                          }
                        else
                          {
                            val0 = shape_hessians[col * n_columns + ind];
                            val1 =
                              shape_hessians[(col + 1) * n_columns - 1 - ind];
                          }
                        in0 = in[stride * ind];
                        in1 = in[stride * (mm - 1 - ind)];
                        res0 += val0 * in0;
                        res1 += val1 * in0;
                        res0 += val1 * in1;
                        res1 += val0 * in1;
                      }
                  }
                else
                  res0 = res1 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_hessians[mid * n_columns + col];
                    else
                      val0 = shape_hessians[col * n_columns + mid];
                    in1 = val0 * in[stride * mid];
                    res0 += in1;
                    res1 += in1;
                  }
                if (add == false)
                  {
                    out[stride * col]            = res0;
                    out[stride * (nn - 1 - col)] = res1;
                  }
                else
                  {
                    out[stride * col] += res0;
                    out[stride * (nn - 1 - col)] += res1;
                  }
              }
            if (nn % 2 == 1)
              {
                Number2 val0;
                Number  res0;
                if (contract_over_rows == true)
                  val0 = shape_hessians[n_cols];
                else
                  val0 = shape_hessians[n_cols * n_columns];
                if (mid > 0)
                  {
                    res0 = val0 * (in[0] + in[stride * (mm - 1)]);
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          val0 = shape_hessians[ind * n_columns + n_cols];
                        else
                          val0 = shape_hessians[n_cols * n_columns + ind];
                        Number in1 = val0 * (in[stride * ind] +
                                             in[stride * (mm - 1 - ind)]);
                        res0 += in1;
                      }
                  }
                else
                  res0 = Number();
                if (mm % 2 == 1)
                  {
                    if (contract_over_rows == true)
                      val0 = shape_hessians[mid * n_columns + n_cols];
                    else
                      val0 = shape_hessians[n_cols * n_columns + mid];
                    res0 += val0 * in[stride * mid];
                  }
                if (add == false)
                  out[stride * n_cols] = res0;
                else
                  out[stride * n_cols] += res0;
              }

            ++in;
            ++out;
          }
        in += stride * (mm - 1);
        out += stride * (nn - 1);
      }
  }



  /**
   * 1d-3d形状函数的内部评估器，使用基函数的张量积形式。
   * 这个类对对称情况下的值、梯度和Hessians也用上述函数处理，实现了一种不同的方法。有可能将每个维度的成本从N^2降低到N^2/2，其中N是一维度的数量（形状矩阵中只有N^2/2个不同的条目，所以这是很合理的）。该方法是基于对输入向量的偶数和奇数部分分别应用算子的想法，因为在正交点上评估的形状函数是对称的。例如，在David
   * A.
   * Kopriva的《实施偏微分方程的频谱方法》一书中介绍了这种方法，Springer,
   * 2009，第3.5.3节（偶数-奇数分解）。尽管书中的实验说该方法在N<20的情况下效率不高，但在循环边界为编译时常量（模板）的情况下，它的效率更高。
   * @tparam  dim 应用该类的空间维度  @tparam  n_rows
   * 变换矩阵中的行数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  n_columns
   * 变换矩阵中的列数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  ] Number 用于输入和输出数组的抽象数字类型
   * @tparam  Number2
   * 用于系数数组的抽象数字类型（默认为与输入/输出数组相同的类型）；必须用Number实现操作符*，并产生Number作为输出，才能成为有效类型。
   *
   */
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  struct EvaluatorTensorProduct<evaluate_evenodd,
                                dim,
                                n_rows,
                                n_columns,
                                Number,
                                Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      Utilities::pow(n_rows, dim);
    static constexpr unsigned int n_columns_of_product =
      Utilities::pow(n_columns, dim);

    /**
     * 空的构造函数。什么都不做。在使用'values'和相关的方法时要小心，因为它们需要用另一个构造函数来填充，至少传入一个数组的值。
     *
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {}

    /**
     * 构造函数，从ShapeInfo中获取数据（使用存储在那里的偶数变体）。
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values)
      : shape_values(shape_values.begin())
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {
      AssertDimension(shape_values.size(), n_rows * ((n_columns + 1) / 2));
    }

    /**
     * 构造函数，从ShapeInfo中获取数据（使用存储在那里的偶数变体）。
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            dummy1 = 0,
                           const unsigned int            dummy2 = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
    {
      // In this function, we allow for dummy pointers if some of values,
      // gradients or hessians should not be computed
      if (!shape_values.empty())
        AssertDimension(shape_values.size(), n_rows * ((n_columns + 1) / 2));
      if (!shape_gradients.empty())
        AssertDimension(shape_gradients.size(), n_rows * ((n_columns + 1) / 2));
      if (!shape_hessians.empty())
        AssertDimension(shape_hessians.size(), n_rows * ((n_columns + 1) / 2));
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 1>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 2>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0, true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 1, true>(shape_gradients,
                                                         in,
                                                         out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 2, true>(shape_hessians,
                                                         in,
                                                         out);
    }

    /**
     * 这个函数沿着输入数组中张量数据的给定 @p direction
     * 应用张量积核，对应于一维条纹的乘法。这个函数允许
     * @p in 和 @p out 数组在n_rows ==
     * n_columns的情况下别名，也就是说，在 @p in 和 @p out
     * 指向相同地址的地方执行收缩是安全的。对于n_rows !=
     * n_columns的情况，只有当 @p one_line
     * 被设置为真时，输出才是正确的。          @tparam  方向
     * 被评估的方向  @tparam  contract_over_rows
     * 如果为真，张量收缩对给定的 @p shape_data
     * 数组中的行求和，否则对列求和  @tparam  add
     * 如果为真，结果被添加到输出向量中，否则计算值将覆盖输出中的内容
     * @tparam  ] type
     * 决定是否使用形状值（type=0）、形状梯度（type=1）或二阶导数（type=2，类似于type
     * 0，但没有两个额外的0条目）中出现的对称性  @tparam
     * one_line 如果为真，内核只沿着dim-dimensional
     * tensor中的单个1D条纹应用，而不是像 @p false
     * 情况中的全部n_rows^dim点。          @param  shape_data 具有
     * @p n_rows 行和 @p n_columns
     * 列的变换矩阵，以行为主的格式存储  @param  in
     * 指向输入数据向量的起点  @param  out
     * 指向输出数据向量的起点
     *
     */
    template <int  direction,
              bool contract_over_rows,
              bool add,
              int  type,
              bool one_line = false>
    static void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number *                  in,
          Number *                        out);

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int  direction,
            bool contract_over_rows,
            bool add,
            int  type,
            bool one_line>
  inline void
  EvaluatorTensorProduct<evaluate_evenodd,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::apply(const Number2 *DEAL_II_RESTRICT shapes,
                                         const Number *                  in,
                                         Number *                        out)
  {
    static_assert(type < 3, "Only three variants type=0,1,2 implemented");
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));

    // We cannot statically assert that direction is less than dim, so must do
    // an additional dynamic check
    AssertIndexRange(direction, dim);

    constexpr int nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int mm     = contract_over_rows ? n_rows : n_columns;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    constexpr int offset = (n_columns + 1) / 2;

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            Number xp[mid > 0 ? mid : 1], xm[mid > 0 ? mid : 1];
            for (int i = 0; i < mid; ++i)
              {
                if (contract_over_rows == true && type == 1)
                  {
                    xp[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                    xm[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                  }
                else
                  {
                    xp[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                    xm[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                  }
              }
            Number xmid = in[stride * mid];
            for (int col = 0; col < n_cols; ++col)
              {
                Number r0, r1;
                if (mid > 0)
                  {
                    if (contract_over_rows == true)
                      {
                        r0 = shapes[col] * xp[0];
                        r1 = shapes[(n_rows - 1) * offset + col] * xm[0];
                      }
                    else
                      {
                        r0 = shapes[col * offset] * xp[0];
                        r1 = shapes[(n_rows - 1 - col) * offset] * xm[0];
                      }
                    for (int ind = 1; ind < mid; ++ind)
                      {
                        if (contract_over_rows == true)
                          {
                            r0 += shapes[ind * offset + col] * xp[ind];
                            r1 += shapes[(n_rows - 1 - ind) * offset + col] *
                                  xm[ind];
                          }
                        else
                          {
                            r0 += shapes[col * offset + ind] * xp[ind];
                            r1 += shapes[(n_rows - 1 - col) * offset + ind] *
                                  xm[ind];
                          }
                      }
                  }
                else
                  r0 = r1 = Number();
                if (mm % 2 == 1 && contract_over_rows == true)
                  {
                    if (type == 1)
                      r1 += shapes[mid * offset + col] * xmid;
                    else
                      r0 += shapes[mid * offset + col] * xmid;
                  }
                else if (mm % 2 == 1 && (nn % 2 == 0 || type > 0 || mm == 3))
                  r0 += shapes[col * offset + mid] * xmid;

                if (add == false)
                  {
                    out[stride * col] = r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride * (nn - 1 - col)] = r1 - r0;
                    else
                      out[stride * (nn - 1 - col)] = r0 - r1;
                  }
                else
                  {
                    out[stride * col] += r0 + r1;
                    if (type == 1 && contract_over_rows == false)
                      out[stride * (nn - 1 - col)] += r1 - r0;
                    else
                      out[stride * (nn - 1 - col)] += r0 - r1;
                  }
              }
            if (type == 0 && contract_over_rows == true && nn % 2 == 1 &&
                mm % 2 == 1 && mm > 3)
              {
                if (add == false)
                  out[stride * n_cols] = shapes[mid * offset + n_cols] * xmid;
                else
                  out[stride * n_cols] += shapes[mid * offset + n_cols] * xmid;
              }
            else if (contract_over_rows == true && nn % 2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    r0 = shapes[n_cols] * xp[0];
                    for (int ind = 1; ind < mid; ++ind)
                      r0 += shapes[ind * offset + n_cols] * xp[ind];
                  }
                else
                  r0 = Number();
                if (type != 1 && mm % 2 == 1)
                  r0 += shapes[mid * offset + n_cols] * xmid;

                if (add == false)
                  out[stride * n_cols] = r0;
                else
                  out[stride * n_cols] += r0;
              }
            else if (contract_over_rows == false && nn % 2 == 1)
              {
                Number r0;
                if (mid > 0)
                  {
                    if (type == 1)
                      {
                        r0 = shapes[n_cols * offset] * xm[0];
                        for (int ind = 1; ind < mid; ++ind)
                          r0 += shapes[n_cols * offset + ind] * xm[ind];
                      }
                    else
                      {
                        r0 = shapes[n_cols * offset] * xp[0];
                        for (int ind = 1; ind < mid; ++ind)
                          r0 += shapes[n_cols * offset + ind] * xp[ind];
                      }
                  }
                else
                  r0 = Number();

                if ((type == 0 || type == 2) && mm % 2 == 1)
                  r0 += shapes[n_cols * offset + mid] * xmid;

                if (add == false)
                  out[stride * n_cols] = r0;
                else
                  out[stride * n_cols] += r0;
              }
            if (one_line == false)
              {
                in += 1;
                out += 1;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm - 1);
            out += stride * (nn - 1);
          }
      }
  }



  /**
   * 使用基函数的张量积形式的1d-3d形状函数的内部评估器。
   * 这个类实现了一种类似于偶数分解的方法，但具有不同类型的对称性。在这种情况下，我们假设单个形状函数已经显示了在正交点上的对称性，而不是在偶数情况下考虑的完整基础。特别是，我们假设形状函数的排序与Legendre基一样，偶数槽（数值阵列的行）中的形状函数是对称的，奇数槽是点对称的。与偶数分解一样，操作的数量是N^2/2，而不是N^2
   * FMAs（融合乘加），其中N是一维度数。区别在于输入和输出量的对称方式。
   * @tparam  dim 应用该类的空间维度  @tparam  n_rows
   * 变换矩阵中的行数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  n_columns
   * 变换矩阵中的列数，对应于通常张量收缩设置中的1d形状函数的数量
   * @tparam  ] Number 用于输入和输出数组的抽象数字类型
   * @tparam  Number2
   * 用于系数数组的抽象数字类型（默认为与输入/输出数组相同的类型）；必须用Number实现操作符*，并产生Number作为输出，才能成为有效类型。
   *
   */
  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  struct EvaluatorTensorProduct<evaluate_symmetric_hierarchical,
                                dim,
                                n_rows,
                                n_columns,
                                Number,
                                Number2>
  {
    static constexpr unsigned int n_rows_of_product =
      Utilities::pow(n_rows, dim);
    static constexpr unsigned int n_columns_of_product =
      Utilities::pow(n_columns, dim);

    /**
     * 空的构造函数。什么都不做。在使用'values'和相关的方法时要小心，因为它们需要用另一个构造函数来填充，至少要传入一个数组的值。
     *
     */
    EvaluatorTensorProduct()
      : shape_values(nullptr)
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {}

    /**
     * 构造函数，从ShapeInfo中获取数据（使用存储在那里的偶数变体）。
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number> &shape_values)
      : shape_values(shape_values.begin())
      , shape_gradients(nullptr)
      , shape_hessians(nullptr)
    {}

    /**
     * 构造函数，从ShapeInfo中获取数据（使用存储在那里的偶数变体）。
     *
     */
    EvaluatorTensorProduct(const AlignedVector<Number2> &shape_values,
                           const AlignedVector<Number2> &shape_gradients,
                           const AlignedVector<Number2> &shape_hessians,
                           const unsigned int            dummy1 = 0,
                           const unsigned int            dummy2 = 0)
      : shape_values(shape_values.begin())
      , shape_gradients(shape_gradients.begin())
      , shape_hessians(shape_hessians.begin())
    {
      (void)dummy1;
      (void)dummy2;
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 1>(shape_gradients, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0>(shape_hessians, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    values_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_values != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0, true>(shape_values, in, out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    gradients_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_gradients != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 1, true>(shape_gradients,
                                                         in,
                                                         out);
    }

    template <int direction, bool contract_over_rows, bool add>
    void
    hessians_one_line(const Number in[], Number out[]) const
    {
      Assert(shape_hessians != nullptr, ExcNotInitialized());
      apply<direction, contract_over_rows, add, 0, true>(shape_hessians,
                                                         in,
                                                         out);
    }

    /**
     * 这个函数沿着输入数组中张量数据的给定 @p direction
     * 应用张量积核，对应于一维条纹的乘法。这个函数允许
     * @p in 和 @p out 数组在n_rows ==
     * n_columns的情况下别名，也就是说，在 @p in 和 @p out
     * 指向相同地址的地方执行收缩是安全的。对于n_rows !=
     * n_columns的情况，只有当 @p one_line
     * 被设置为真时，输出才是正确的。          @tparam  方向
     * 被评估的方向  @tparam  contract_over_rows
     * 如果为真，张量收缩对给定的 @p shape_data
     * 数组中的行求和，否则对列求和  @tparam  add
     * 如果为真，结果被添加到输出向量中，否则计算值覆盖输出的内容
     * @tparam  ] type 决定评估是在 @p shape_data
     * 的偶数行（type=0）还是奇数行（type=1）中对称，以及在奇数行（type=0）或偶数行（type=1）中偏斜对称
     * @tparam  one_line
     * 如果为真，内核只沿着dim-dimensional张量中的单个1D条纹应用，而不是像
     * @p false 情况中的全部n_rows^dim点。          @param  shape_data
     * 具有 @p n_rows 行和 @p n_columns
     * 列的变换矩阵，以行为主的格式存储  @param  in
     * 指向输入数据向量的起点  @param  out
     * 指向输出数据向量的起点
     *
     */
    template <int  direction,
              bool contract_over_rows,
              bool add,
              int  type,
              bool one_line = false>
    static void
    apply(const Number2 *DEAL_II_RESTRICT shape_data,
          const Number *                  in,
          Number *                        out);

    const Number2 *shape_values;
    const Number2 *shape_gradients;
    const Number2 *shape_hessians;
  };



  template <int dim,
            int n_rows,
            int n_columns,
            typename Number,
            typename Number2>
  template <int  direction,
            bool contract_over_rows,
            bool add,
            int  type,
            bool one_line>
  inline void
  EvaluatorTensorProduct<evaluate_symmetric_hierarchical,
                         dim,
                         n_rows,
                         n_columns,
                         Number,
                         Number2>::apply(const Number2 *DEAL_II_RESTRICT shapes,
                                         const Number *                  in,
                                         Number *                        out)
  {
    static_assert(one_line == false || direction == dim - 1,
                  "Single-line evaluation only works for direction=dim-1.");
    static_assert(
      type == 0 || type == 1,
      "Only types 0 and 1 implemented for evaluate_symmetric_hierarchical.");
    Assert(dim == direction + 1 || one_line == true || n_rows == n_columns ||
             in != out,
           ExcMessage("In-place operation only supported for "
                      "n_rows==n_columns or single-line interpolation"));

    // We cannot statically assert that direction is less than dim, so must do
    // an additional dynamic check
    AssertIndexRange(direction, dim);

    constexpr int nn     = contract_over_rows ? n_columns : n_rows;
    constexpr int mm     = contract_over_rows ? n_rows : n_columns;
    constexpr int n_cols = nn / 2;
    constexpr int mid    = mm / 2;

    constexpr int stride    = Utilities::pow(n_columns, direction);
    constexpr int n_blocks1 = one_line ? 1 : stride;
    constexpr int n_blocks2 =
      Utilities::pow(n_rows, (direction >= dim) ? 0 : (dim - direction - 1));

    // this code may look very inefficient at first sight due to the many
    // different cases with if's at the innermost loop part, but all of the
    // conditionals can be evaluated at compile time because they are
    // templates, so the compiler should optimize everything away
    for (int i2 = 0; i2 < n_blocks2; ++i2)
      {
        for (int i1 = 0; i1 < n_blocks1; ++i1)
          {
            if (contract_over_rows)
              {
                Number x[mm];
                for (unsigned int i = 0; i < mm; ++i)
                  x[i] = in[stride * i];
                for (unsigned int col = 0; col < n_cols; ++col)
                  {
                    Number r0, r1;
                    if (mid > 0)
                      {
                        r0 = shapes[col] * x[0];
                        r1 = shapes[col + n_columns] * x[1];
                        for (unsigned int ind = 1; ind < mid; ++ind)
                          {
                            r0 +=
                              shapes[col + 2 * ind * n_columns] * x[2 * ind];
                            r1 += shapes[col + (2 * ind + 1) * n_columns] *
                                  x[2 * ind + 1];
                          }
                      }
                    else
                      r0 = r1 = Number();
                    if (mm % 2 == 1)
                      r0 += shapes[col + (mm - 1) * n_columns] * x[mm - 1];
                    if (add == false)
                      {
                        out[stride * col] = r0 + r1;
                        if (type == 1)
                          out[stride * (nn - 1 - col)] = r1 - r0;
                        else
                          out[stride * (nn - 1 - col)] = r0 - r1;
                      }
                    else
                      {
                        out[stride * col] += r0 + r1;
                        if (type == 1)
                          out[stride * (nn - 1 - col)] += r1 - r0;
                        else
                          out[stride * (nn - 1 - col)] += r0 - r1;
                      }
                  }
                if (nn % 2 == 1)
                  {
                    Number             r0;
                    const unsigned int shift = type == 1 ? 1 : 0;
                    if (mid > 0)
                      {
                        r0 = shapes[n_cols + shift * n_columns] * x[shift];
                        for (unsigned int ind = 1; ind < mid; ++ind)
                          r0 += shapes[n_cols + (2 * ind + shift) * n_columns] *
                                x[2 * ind + shift];
                      }
                    else
                      r0 = 0;
                    if (type != 1 && mm % 2 == 1)
                      r0 += shapes[n_cols + (mm - 1) * n_columns] * x[mm - 1];
                    if (add == false)
                      out[stride * n_cols] = r0;
                    else
                      out[stride * n_cols] += r0;
                  }
              }
            else
              {
                Number xp[mid + 1], xm[mid > 0 ? mid : 1];
                for (int i = 0; i < mid; ++i)
                  if (type == 0)
                    {
                      xp[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                      xm[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                    }
                  else
                    {
                      xp[i] = in[stride * i] - in[stride * (mm - 1 - i)];
                      xm[i] = in[stride * i] + in[stride * (mm - 1 - i)];
                    }
                if (mm % 2 == 1)
                  xp[mid] = in[stride * mid];
                for (unsigned int col = 0; col < n_cols; ++col)
                  {
                    Number r0, r1;
                    if (mid > 0)
                      {
                        r0 = shapes[2 * col * n_columns] * xp[0];
                        r1 = shapes[(2 * col + 1) * n_columns] * xm[0];
                        for (unsigned int ind = 1; ind < mid; ++ind)
                          {
                            r0 += shapes[2 * col * n_columns + ind] * xp[ind];
                            r1 +=
                              shapes[(2 * col + 1) * n_columns + ind] * xm[ind];
                          }
                      }
                    else
                      r0 = r1 = Number();
                    if (mm % 2 == 1)
                      {
                        if (type == 1)
                          r1 +=
                            shapes[(2 * col + 1) * n_columns + mid] * xp[mid];
                        else
                          r0 += shapes[2 * col * n_columns + mid] * xp[mid];
                      }
                    if (add == false)
                      {
                        out[stride * (2 * col)]     = r0;
                        out[stride * (2 * col + 1)] = r1;
                      }
                    else
                      {
                        out[stride * (2 * col)] += r0;
                        out[stride * (2 * col + 1)] += r1;
                      }
                  }
                if (nn % 2 == 1)
                  {
                    Number r0;
                    if (mid > 0)
                      {
                        r0 = shapes[(nn - 1) * n_columns] * xp[0];
                        for (unsigned int ind = 1; ind < mid; ++ind)
                          r0 += shapes[(nn - 1) * n_columns + ind] * xp[ind];
                      }
                    else
                      r0 = Number();
                    if (mm % 2 == 1 && type == 0)
                      r0 += shapes[(nn - 1) * n_columns + mid] * xp[mid];
                    if (add == false)
                      out[stride * (nn - 1)] = r0;
                    else
                      out[stride * (nn - 1)] += r0;
                  }
              }
            if (one_line == false)
              {
                in += 1;
                out += 1;
              }
          }
        if (one_line == false)
          {
            in += stride * (mm - 1);
            out += stride * (nn - 1);
          }
      }
  }



  /**
   * 在 evaluate_tensor_product_value_and_gradient 中避免使用 Tensor<1,
   * dim, Point<dim2> 的结构，因为点不能在 Tensor
   * 中使用。相反，这个结构的特殊化将点上传到一个Tensor<1,dim>。
   *
   */
  template <typename Number, typename Number2>
  struct ProductTypeNoPoint
  {
    using type = typename ProductType<Number, Number2>::type;
  };

  template <int dim, typename Number, typename Number2>
  struct ProductTypeNoPoint<Point<dim, Number>, Number2>
  {
    using type = typename ProductType<Tensor<1, dim, Number>, Number2>::type;
  };



  /**
   * 计算张量积形状函数  $\varphi_i$
   * 的多项式插值，给定系数向量  $u_i$  的形式
   * $u_h(\mathbf{x}) = \sum_{i=1}^{k^d} \varphi_i(\mathbf{x}) u_i$
   * 。形状函数 $\varphi_i(\mathbf{x}) =
   * \prod_{d=1}^{\text{dim}}\varphi_{i_d}^\text{1D}(x_d)$
   * 代表张量积。该函数返回一对，内插值为第一分量，参考坐标中的梯度为第二分量。注意，对于复合类型（例如，`values`字段开始一个Point<spacedim>参数），梯度的分量被排序为Tensor<1,
   * dim, Tensor<1,
   * spacedim>，导数为第一个索引；这是函数中通用参数的结果。
   * @param  poly 基础的一维多项式基  $\{\varphi^{1D}_{i_1}\}$
   * 以多项式的矢量形式给出。      @param  values
   * 多项式插值中类型为`Number`的扩展系数  $u_i$
   * 。这些系数可以是简单的 "双
   * "变量，但也可以是Point<spacedim>，如果它们定义了类型为
   * "Number2 "的算术运算。      @param  p
   * 在参考坐标中应该评估插值的位置。      @param  d_linear
   * 指定是否应该进行d-线性（一维的线性，二维的双线性，三维的三线性）插值，这允许解开循环并大大加快评估。
   * @param  renumber
   * 可选参数，用于指定系数向量中的重新编号，假设`values[renumber[i]]返回系数的lexicographic（张量积）条目。如果该向量为条目，则假设数值按词典排序。
   *
   */
  template <int dim, typename Number, typename Number2>
  inline std::pair<
    typename ProductTypeNoPoint<Number, Number2>::type,
    Tensor<1, dim, typename ProductTypeNoPoint<Number, Number2>::type>>
  evaluate_tensor_product_value_and_gradient(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const std::vector<Number> &                         values,
    const Point<dim, Number2> &                         p,
    const bool                                          d_linear = false,
    const std::vector<unsigned int> &                   renumber = {})
  {
    static_assert(dim >= 1 && dim <= 3, "Only dim=1,2,3 implemented");

    using Number3 = typename ProductTypeNoPoint<Number, Number2>::type;

    // use `int` type for this variable and the loops below to inform the
    // compiler that the loops below will never overflow, which allows it to
    // generate more optimized code for the variable loop bounds in the
    // present context
    const int n_shapes = poly.size();
    AssertDimension(Utilities::pow(n_shapes, dim), values.size());
    Assert(renumber.empty() || renumber.size() == values.size(),
           ExcDimensionMismatch(renumber.size(), values.size()));

    // shortcut for linear interpolation to speed up evaluation
    if (d_linear)
      {
        AssertDimension(poly.size(), 2);
        for (unsigned int i = 0; i < renumber.size(); ++i)
          AssertDimension(renumber[i], i);

        if (dim == 1)
          {
            Tensor<1, dim, Number3> derivative;
            derivative[0] = values[1] - values[0];
            return std::make_pair((1. - p[0]) * values[0] + p[0] * values[1],
                                  derivative);
          }
        else if (dim == 2)
          {
            const Number2           x0 = 1. - p[0], x1 = p[0];
            const Number3           tmp0   = x0 * values[0] + x1 * values[1];
            const Number3           tmp1   = x0 * values[2] + x1 * values[3];
            const Number3           mapped = (1. - p[1]) * tmp0 + p[1] * tmp1;
            Tensor<1, dim, Number3> derivative;
            derivative[0] = (1. - p[1]) * (values[1] - values[0]) +
                            p[1] * (values[3] - values[2]);
            derivative[1] = tmp1 - tmp0;
            return std::make_pair(mapped, derivative);
          }
        else if (dim == 3)
          {
            const Number2 x0 = 1. - p[0], x1 = p[0], y0 = 1. - p[1], y1 = p[1],
                          z0 = 1. - p[2], z1 = p[2];
            const Number3           tmp0   = x0 * values[0] + x1 * values[1];
            const Number3           tmp1   = x0 * values[2] + x1 * values[3];
            const Number3           tmpy0  = y0 * tmp0 + y1 * tmp1;
            const Number3           tmp2   = x0 * values[4] + x1 * values[5];
            const Number3           tmp3   = x0 * values[6] + x1 * values[7];
            const Number3           tmpy1  = y0 * tmp2 + y1 * tmp3;
            const Number3           mapped = z0 * tmpy0 + z1 * tmpy1;
            Tensor<1, dim, Number3> derivative;
            derivative[2] = tmpy1 - tmpy0;
            derivative[1] = z0 * (tmp1 - tmp0) + z1 * (tmp3 - tmp2);
            derivative[0] =
              z0 *
                (y0 * (values[1] - values[0]) + y1 * (values[3] - values[2])) +
              z1 *
                (y0 * (values[5] - values[4]) + y1 * (values[7] - values[6]));
            return std::make_pair(mapped, derivative);
          }
      }

    AssertIndexRange(n_shapes, 200);
    std::array<Number2, 2 * dim * 200> shapes;

    // Evaluate 1D polynomials and their derivatives
    for (unsigned int d = 0; d < dim; ++d)
      for (int i = 0; i < n_shapes; ++i)
        poly[i].value(p[d], 1, shapes.data() + 2 * (d * n_shapes + i));

    // Go through the tensor product of shape functions and interpolate
    // with optimal algorithm
    std::pair<Number3, Tensor<1, dim, Number3>> result = {};
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        Number3 value_y = {}, deriv_x = {}, deriv_y = {};
        for (int i1 = 0; i1 < (dim > 1 ? n_shapes : 1); ++i1)
          {
            // Interpolation + derivative x direction
            Number3 value = {}, deriv = {};

            // Distinguish the inner loop based on whether we have a
            // renumbering or not
            if (renumber.empty())
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                {
                  value += shapes[2 * i0] * values[i];
                  deriv += shapes[2 * i0 + 1] * values[i];
                }
            else
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                {
                  value += shapes[2 * i0] * values[renumber[i]];
                  deriv += shapes[2 * i0 + 1] * values[renumber[i]];
                }

            // Interpolation + derivative in y direction
            if (dim > 1)
              {
                value_y += value * shapes[2 * n_shapes + 2 * i1];
                deriv_x += deriv * shapes[2 * n_shapes + 2 * i1];
                deriv_y += value * shapes[2 * n_shapes + 2 * i1 + 1];
              }
            else
              {
                result.first     = value;
                result.second[0] = deriv;
              }
          }
        if (dim == 3)
          {
            // Interpolation + derivative in z direction
            result.first += value_y * shapes[4 * n_shapes + 2 * i2];
            result.second[0] += deriv_x * shapes[4 * n_shapes + 2 * i2];
            result.second[1] += deriv_y * shapes[4 * n_shapes + 2 * i2];
            result.second[2] += value_y * shapes[4 * n_shapes + 2 * i2 + 1];
          }
        else if (dim == 2)
          {
            result.first     = value_y;
            result.second[0] = deriv_x;
            result.second[1] = deriv_y;
          }
      }

    return result;
  }



  /**
   * 与evaluate_tensor_product_value_and_gradient()相同，但用于积分。
   *
   */
  template <int dim, typename Number, typename Number2>
  inline void
  integrate_add_tensor_product_value_and_gradient(
    const std::vector<Polynomials::Polynomial<double>> &poly,
    const Number2 &                                     value,
    const Tensor<1, dim, Number2> &                     gradient,
    const Point<dim, Number> &                          p,
    AlignedVector<Number2> &                            values,
    const std::vector<unsigned int> &                   renumber = {})
  {
    static_assert(dim >= 1 && dim <= 3, "Only dim=1,2,3 implemented");

    // as in evaluate, use `int` type to produce better code in this context
    const int n_shapes = poly.size();
    AssertDimension(Utilities::pow(n_shapes, dim), values.size());
    Assert(renumber.empty() || renumber.size() == values.size(),
           ExcDimensionMismatch(renumber.size(), values.size()));

    AssertIndexRange(n_shapes, 200);
    std::array<Number, 2 * dim * 200> shapes;

    // Evaluate 1D polynomials and their derivatives
    for (unsigned int d = 0; d < dim; ++d)
      for (int i = 0; i < n_shapes; ++i)
        poly[i].value(p[d], 1, shapes.data() + 2 * (d * n_shapes + i));

    // Implement the transpose of the function above
    for (int i2 = 0, i = 0; i2 < (dim > 2 ? n_shapes : 1); ++i2)
      {
        const Number2 test_value_z =
          dim > 2 ? (value * shapes[4 * n_shapes + 2 * i2] +
                     gradient[2] * shapes[4 * n_shapes + 2 * i2 + 1]) :
                    value;
        const Number2 test_grad_x =
          dim > 2 ? gradient[0] * shapes[4 * n_shapes + 2 * i2] : gradient[0];
        const Number2 test_grad_y =
          dim > 2 ? gradient[1] * shapes[4 * n_shapes + 2 * i2] :
                    (dim > 1 ? gradient[1] : Number2());
        for (int i1 = 0; i1 < (dim > 1 ? n_shapes : 1); ++i1)
          {
            const Number2 test_value_y =
              dim > 1 ? (test_value_z * shapes[2 * n_shapes + 2 * i1] +
                         test_grad_y * shapes[2 * n_shapes + 2 * i1 + 1]) :
                        test_value_z;
            const Number2 test_grad_xy =
              dim > 1 ? test_grad_x * shapes[2 * n_shapes + 2 * i1] :
                        test_grad_x;
            if (renumber.empty())
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                values[i] += shapes[2 * i0] * test_value_y +
                             shapes[2 * i0 + 1] * test_grad_xy;
            else
              for (int i0 = 0; i0 < n_shapes; ++i0, ++i)
                values[renumber[i]] += shapes[2 * i0] * test_value_y +
                                       shapes[2 * i0 + 1] * test_grad_xy;
          }
      }
  }


} // end of namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif


