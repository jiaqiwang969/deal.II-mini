//include/deal.II-translator/physics/notation_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#ifndef dealii_physics_notation_h
#define dealii_physics_notation_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <type_traits>

DEAL_II_NAMESPACE_OPEN


namespace Physics
{
  namespace Notation
  {
    /**
     * @brief A namespace with functions that assist in the conversion of
     * 向量和张量使用Kelvin记号和加权的压缩格式来回转换。
     * Kelvin和Voigt符号都采用相同的索引惯例。
     * 具体到空间维度3的情况，对于一个等级2的对称张量
     * $\mathbf{S}$ ，我们列举其张量成分@f[ \mathbf{S} \dealcoloneq
     * \left[ \begin{array}{ccc} S_{00}          & S_{01}          & S_{02} \\
     * S_{10} = S_{01} & S_{11}          & S_{12} \\ S_{20} = S_{02} & S_{21}
     * = S_{12} & S_{22} \end{array} \right] \quad \Rightarrow \quad \left[
     * \begin{array}{ccc} n = 0 & n = 5 & n = 4 \\ sym   & n = 1 & n = 3 \\
     * sym   & sym   & n = 2 \end{array} \right] ,
     * @f]，其中 $n$ 表示张量成分的开尔文索引，而对于一般等级2的张量 $\mathbf{T}$  ] @f[
     * \mathbf{T} \dealcoloneq \left[ \begin{array}{ccc} T_{00} & T_{01} &
     * T_{02} \\ T_{10} & T_{11} & T_{12} \\ T_{20} & T_{21} & T_{22}
     * \end{array}\right] \quad \Rightarrow \quad \left[ \begin{array}{ccc} n
     * = 0 & n = 5 & n = 4 \\ n = 6 & n = 1 & n = 3 \\ n = 7 & n = 8 & n = 2
     * \end{array}\right] ,
     * @f]，对于秩-1张量 $\mathbf{v}$ @f[
     * \mathbf{v} \dealcoloneq \left[ \begin{array}{c} v_{0} \\ v_{1} \\ v_{2}
     * \end{array}\right] \quad \Rightarrow \quad \left[ \begin{array}{c} n =
     * 0 \\ n = 1 \\ n = 2 \end{array}\right] . @f]
     * 总结一下，张量和开尔文指数之间的关系对于三维情况和类似的辨别二维情况概述在下面的表中。
     * <table> <tr> <th align="center"> Dimension 2 </th> <th align="center">
     * Dimension 3 </th> </tr> <tr> <td align="middle"> <table> <tr>
     * <th>Tensor index pairs</th> <th>Kelvin index</th> </tr> <tr> <td
     * align="center">00</td> <td align="center">0</td> </tr> <tr> <td
     * align="center">11</td> <td align="center">1</td> </tr> <tr> <td
     * align="center">12</td> <td align="center">2</td> </tr> <tr> <td
     * align="center">21</td> <td align="center">3</td> </tr> </table>  </td>
     * <td align="middle">  <table> <tr> <th>Tensor index pairs</th>
     * <th>Kelvin index</th> </tr> <tr> <td align="center">00</td> <td
     * align="center">0</td> </tr> <tr> <td align="center">11</td> <td
     * align="center">1</td> </tr> <tr> <td align="center">22</td> <td
     * align="center">2</td> </tr> <tr> <td align="center">12</td> <td
     * align="center">3</td> </tr> <tr> <td align="center">02</td> <td
     * align="center">4</td> </tr> <tr> <td align="center">01</td> <td
     * align="center">5</td> </tr> <tr> <td align="center">10</td> <td
     * align="center">6</td> </tr> <tr> <td align="center">20</td> <td
     * align="center">7</td> </tr> <tr> <td align="center">21</td> <td
     * align="center">8</td> </tr>
     * </table>  </td> </tr> </table>为了说明这个符号的目的，考虑等级2的对称张量 $\mathbf{S}$ 和 $\mathbf{E}$ ，它们之间的关系是 $\mathbf{S} = \cal{C} : \mathbf{E}$  ，其中算子 $\cal{C}$  是一个四阶对称的张量。相对于常用的Voigt符号，当 $\mathbf{S}$ 和 $\mathbf{E}$ 都是对称的时候，Kelvin（或Mandel）符号对内积 $\mathbf{S} : \mathbf{E}$ 的定义保持不变。一般来说，所有对称和一般张量的内积都是一样的，无论用什么符号表示。        为了实现这两个属性，即@f[
     * \mathbf{S} = \cal{C} : \mathbf{E}
     * \quad \Rightarrow   \quad
     * \tilde{\mathbf{S}} = \tilde{\cal{C}} \; \tilde{\mathbf{E}}
     * @f]和@f[
     * \mathbf{S} : \mathbf{E} \, \equiv \, \tilde{\mathbf{S}} \cdot
     * \tilde{\mathbf{E}} ,
     * @f]，它认为之前定义的对称张量的开尔文缩合等价物，由 $\tilde{\left(\bullet\right)}$ 表示，必须定义为@f[
     * \tilde{\mathbf{S}} = \left[ \begin{array}{c} S_{00} \\ S_{11} \\ S_{22}
     * \\ \sqrt{2} S_{12} \\ \sqrt{2} S_{02} \\ \sqrt{2} S_{01}
     * \end{array}\right] \quad \text{and} \quad \tilde{\mathbf{E}} = \left[
     * \begin{array}{c} E_{00} \\ E_{11} \\ E_{22} \\ \sqrt{2} E_{12} \\
     * \sqrt{2} E_{02} \\ \sqrt{2} E_{01} \end{array}\right] .
     * @f] ] 相应的、一致的浓缩四阶对称张量是@f[
     * \tilde{\cal{C}} = \left[ \begin{array}{cccccc} \tilde{\cal{C}}_{00} &
     * \tilde{\cal{C}}_{01} & \tilde{\cal{C}}_{02} & \tilde{\cal{C}}_{03} &
     * \tilde{\cal{C}}_{04} & \tilde{\cal{C}}_{05} \\ \tilde{\cal{C}}_{10} &
     * \tilde{\cal{C}}_{11} & \tilde{\cal{C}}_{12} & \tilde{\cal{C}}_{13} &
     * \tilde{\cal{C}}_{14} & \tilde{\cal{C}}_{15} \\ \tilde{\cal{C}}_{20} &
     * \tilde{\cal{C}}_{21} & \tilde{\cal{C}}_{22} & \tilde{\cal{C}}_{23} &
     * \tilde{\cal{C}}_{24} & \tilde{\cal{C}}_{25} \\ \tilde{\cal{C}}_{30} &
     * \tilde{\cal{C}}_{31} & \tilde{\cal{C}}_{32} & \tilde{\cal{C}}_{33} &
     * \tilde{\cal{C}}_{34} & \tilde{\cal{C}}_{35} \\ \tilde{\cal{C}}_{40} &
     * \tilde{\cal{C}}_{41} & \tilde{\cal{C}}_{42} & \tilde{\cal{C}}_{43} &
     * \tilde{\cal{C}}_{44} & \tilde{\cal{C}}_{45} \\ \tilde{\cal{C}}_{50} &
     * \tilde{\cal{C}}_{51} & \tilde{\cal{C}}_{52} & \tilde{\cal{C}}_{53} &
     * \tilde{\cal{C}}_{54} & \tilde{\cal{C}}_{55} \end{array}\right] \equiv
     * \left[ \begin{array}{cccccc} {\cal{C}}_{0000}           &
     * {\cal{C}}_{0011}          & {\cal{C}}_{0022}           & \sqrt{2}
     * {\cal{C}}_{0012}  & \sqrt{2} {\cal{C}}_{0002}  & \sqrt{2}
     * {\cal{C}}_{0001} \\ {\cal{C}}_{1100}           & {\cal{C}}_{1111}
     * & {\cal{C}}_{1122}           & \sqrt{2} {\cal{C}}_{1112}  & \sqrt{2}
     * {\cal{C}}_{1102}  & \sqrt{2} {\cal{C}}_{1101} \\ {\cal{C}}_{2200}
     * & {\cal{C}}_{2211}          & {\cal{C}}_{2222}           & \sqrt{2}
     * {\cal{C}}_{2212}  & \sqrt{2} {\cal{C}}_{2202}  & \sqrt{2}
     * {\cal{C}}_{2201} \\ \sqrt{2} {\cal{C}}_{1200}  & \sqrt{2}
     * {\cal{C}}_{1211} & \sqrt{2} {\cal{C}}_{1222}  & 2 {\cal{C}}_{1212}
     * & 2 {\cal{C}}_{1202} & 2 {\cal{C}}_{1201}        \\ \sqrt{2}
     * {\cal{C}}_{0200}  & \sqrt{2} {\cal{C}}_{0211} & \sqrt{2}
     * {\cal{C}}_{0222}  & 2 {\cal{C}}_{0212}         & 2 {\cal{C}}_{0202} & 2
     * {\cal{C}}_{0201}        \\ \sqrt{2} {\cal{C}}_{0100}  & \sqrt{2}
     * {\cal{C}}_{0111} & \sqrt{2} {\cal{C}}_{0122}  & 2 {\cal{C}}_{0112} & 2
     * {\cal{C}}_{0102}         & 2 {\cal{C}}_{0101} \end{array}\right] . @f]
     * 从FullMatrix $\tilde{\cal{C}}$
     * 的两个开尔文指数到等级4的SymmetricTensor $\cal{C}$
     * 的映射可以用上面的表格推断出来。
     * 一个重要的观察是，左侧张量 $\tilde{\mathbf{S}}$
     * 和右侧张量 $\tilde{\mathbf{E}}$
     * 都有相同的形式；这是一个在Voigt符号中不存在的属性。
     * 引入 $\tilde{\mathbf{S}}$ 、 $\tilde{\mathbf{E}}$ 和
     * $\tilde{\cal{C}}$
     * 的各种因素说明了张量的对称性。对其非对称对应物的开尔文描述不包括这些因素。
     * 一些有用的参考资料显示了这种符号的作用，其中包括。
     * @code{.bib}
     * @article{Nagel2016,
     * author  = {Nagel, T. and G{\"o}rke, U-J. and Moerman, K. and Kolditz,
     *            O.},
     * title   = {On advantages of the Kelvin mapping in finite element
     *            implementations of deformation processes},
     * journal = {Environmental Earth Sciences},
     * year    = {2016},
     * volume  = {75},
     * number  = {11},
     * pages   = {937}
     * }
     * @endcode
     * 和
     * @code{.bib}
     * @article{Dellinger1998,
     * author  = {Dellinger, J. and Vasicek, D. and Sondergeld, C.},
     * title   = {Kelvin notation for stabilizing elastic-constant inversion},
     * journal = {Revue de l'Institut Fran{\c{c}}ais du P{\'e}trole},
     * year    = {1998},
     * volume  = {53},
     * number  = {5},
     * pages   = {709--719},
     * url     = {http://sepwww.stanford.edu/oldsep/joe/Reprints/8IWSA.pdf},
     * }
     * @endcode
     * 以及在<a
     * href="https://en.wikipedia.org/wiki/Voigt_notation#Mandel_notation">this
     * wikipedia page</a>和<a
     * href="https://github.com/dealii/dealii/tree/master/tests/physics/notation-kelvin_02.cc">the
     * unit tests</a>上发现的在线参考。
     *
     */
    namespace Kelvin
    {
      /**
       * 输入矩阵的行数不正确。
       *
       */
      DeclException3(ExcNotationExcFullMatrixToTensorRowSize2,
                     int,
                     int,
                     int,
                     << "The number of rows in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << " or " << arg3
                     << ".");


      /**
       * 输入的矩阵有不正确的行数。
       *
       */
      DeclException4(ExcNotationExcFullMatrixToTensorRowSize3,
                     int,
                     int,
                     int,
                     int,
                     << "The number of rows in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << "," << arg3
                     << ", or " << arg4 << ".");


      /**
       * 输入的矩阵有不正确的列数。
       *
       */
      DeclException3(ExcNotationExcFullMatrixToTensorColSize2,
                     int,
                     int,
                     int,
                     << "The number of columns in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << " or " << arg3
                     << ".");


      /**
       * 输入的矩阵有不正确的列数。
       *
       */
      DeclException4(ExcNotationExcFullMatrixToTensorColSize3,
                     int,
                     int,
                     int,
                     int,
                     << "The number of columns in the input matrix is " << arg1
                     << ", but needs to be either " << arg2 << "," << arg3
                     << ", or " << arg4 << ".");


      /**
       * @name  正向操作。张量符号到开尔文符号
       *
       */
      //@{

      /**
       * 将一个标量值转换为其压缩的矢量等值。
       * 输出向量有一个条目。
       *
       */
      template <typename Number>
      Vector<Number>
      to_vector(const Number &s);


      /**
       * 将一个0级张量转换为其压缩向量等价物。
       * 输出向量有一个条目。
       *
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<0, dim, Number> &s);


      /**
       * 将一个秩-1张量转换为其压缩向量等价物。
       * 输出向量有 $dim$ 个条目。
       *
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<1, dim, Number> &v);


      /**
       * 将一个秩-2张量转换为其压缩向量等价物。
       * 输出向量有 Tensor<2,dim>::n_independent_components 个条目。
       *
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<2, dim, Number> &t);


      /**
       * 将一个等级2的对称张量转换为其压缩向量等价物。
       * 输出的向量有 SymmetricTensor<2,dim>::n_independent_components
       * 个条目。
       *
       */
      template <int dim, typename Number>
      Vector<Number>
      to_vector(const SymmetricTensor<2, dim, Number> &st);


      /**
       * 将一个标量值转换为它的压缩矩阵等值。
       * 输出的矩阵将有一行和一列。
       *
       */
      template <typename Number>
      FullMatrix<Number>
      to_matrix(const Number &s);


      /**
       * 将一个0级张量转换为其压缩矩阵等价物。
       * 输出的矩阵将有一行和一列。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<0, dim, Number> &s);


      /**
       * 将一个秩-1张量转换为其压缩矩阵等价物。
       * 输出的矩阵将有 $dim$ 行和一列。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<1, dim, Number> &v);


      /**
       * 将一个秩-2张量转换为其压缩矩阵等价物。
       * 输出的矩阵将有 $dim$ 行和 $dim$ 列。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<2, dim, Number> &t);


      /**
       * 将一个2级对称张量转换为其压缩矩阵等价物。
       * 输出矩阵将有 $dim$ 行和 $dim$
       * 列，其格式与非对称张量的等价函数相同。这是因为不可能将
       * SymmetricTensor<2,dim>::n_independent_components
       * 的唯一条目压缩到一个方形矩阵中。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<2, dim, Number> &st);

      /**
       * 将一个等级3的张量转换为其压缩的矩阵等价物。
       * 模板参数 @p SubTensor1 和 @p SubTensor2
       * 决定了解卷的发生方式，特别是如何解释秩-3张量的元素。
       * 因此，举例来说，在以下两个转换中
       * @code
       * Tensor<3,dim> r3_tnsr;      // All elements filled differently
       * Tensor<3,dim> r3_symm_tnsr; // Some elements filled symmetrically
       *
       * const FullMatrix<double> mtrx_1 =
       * Physics::Notation::to_matrix<dim,
       *                              Tensor<2,dim>,
       *                              Tensor<1,dim>*>(r3_tnsr);
       * const FullMatrix<double> mtrx_2 =
       * Physics::Notation::to_matrix<dim,
       *                              Tensor<1,dim>,
       *                              SymmetricTensor<2,dim>*>(r3_symm_tnsr);
       * @endcode
       * 矩阵 @p mtrx_1 将有 $dim \times dim$ 行和 $dim$
       * 列（即大小为 Tensor<2,dim>::n_independent_components  $\times$
       * Tensor<1,dim>::n_independent_components), ，而那些矩阵 @p
       * mtrx_2 将有 $dim$  ]行和 $(dim \times dim + dim)/2$
       * 列（即大小为 Tensor<1,dim>::n_independent_components  $\times$
       * SymmetricTensor<2,dim>::n_independent_components),
       * ，因为假定对应于第二和第三指数交替的条目是相等的。这就是说，
       * <code>r3_symm_tnsr[i][j][k] == r3_symm_tnsr[i][k][j]</code>  .
       *
       */
      template <int dim,
                typename SubTensor1 = Tensor<2, dim>,
                typename SubTensor2 = Tensor<1, dim>,
                typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<3, dim, Number> &t);


      /**
       * 将一个等级4的张量转换为其压缩矩阵等价物。
       * 输出的矩阵将有 Tensor<2,dim>::n_independent_components 行和
       * Tensor<2,dim>::n_independent_components 列。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<4, dim, Number> &t);


      /**
       * 将一个4级对称张量转换为其压缩矩阵等价物。
       * 输出的矩阵将有
       * SymmetricTensor<2,dim>::n_independent_components 行和
       * SymmetricTensor<2,dim>::n_independent_components 列。
       *
       */
      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<4, dim, Number> &st);

      //@}

      /**
       * @name  反向操作。开尔文符号到张量符号
       *
       */
      //@{

      /**
       * 将一个压缩的向量转换为其等效的标量值。
       *
       */
      template <typename Number>
      void
      to_tensor(const Vector<Number> &vec, Number &s);


      /**
       * 将一个压缩的向量转换为其等价的秩0张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<0, dim, Number> &s);


      /**
       * 将一个压缩的向量转换为其等效的秩-1张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<1, dim, Number> &v);


      /**
       * 将一个压缩的向量转换为它的等效秩-2张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<2, dim, Number> &t);


      /**
       * 将一个压缩的向量转换为其等效的秩-2对称张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, SymmetricTensor<2, dim, Number> &st);


      /**
       * 将一个压缩的矩阵转换为它的等效标量值。
       *
       */
      template <typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Number &s);


      /**
       * 将一个压缩矩阵转换为它的等效秩0张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<0, dim, Number> &s);


      /**
       * 将一个压缩的矩阵转换为其等效的秩-1张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<1, dim, Number> &v);


      /**
       * 将一个压缩矩阵转换为它的等效秩-2张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<2, dim, Number> &t);


      /**
       * 将一个压缩矩阵转换为其等效的秩-2对称张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<2, dim, Number> &st);


      /**
       * 将一个压缩矩阵转换为其等效的秩-3张量。
       * @note 基于矩阵 @p mtrx, 的大小， @p t
       * 的一些组件可以被解释为具有对称的对应关系。这与对应的to_matrix()函数的文档中解释的操作相反。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<3, dim, Number> &t);


      /**
       * 将一个压缩矩阵转换为其等价的等级4张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<4, dim, Number> &t);


      /**
       * 将一个压缩矩阵转换为其等效的秩-4对称张量。
       *
       */
      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<4, dim, Number> &st);


      /**
       * 一个通用的辅助函数，可以将一个压缩的向量转换为其等价物
       * @p TensorType.  。
       *
       */
      template <typename TensorType, typename Number>
      TensorType
      to_tensor(const Vector<Number> &vec);


      /**
       * 一个通用的辅助函数，将压缩矩阵转换为其等价物
       * @p TensorType.  。
       *
       */
      template <typename TensorType, typename Number>
      TensorType
      to_tensor(const FullMatrix<Number> &vec);
      //@}

    } // namespace Kelvin

  } // namespace Notation
} // namespace Physics


#ifndef DOXYGEN


// ------------------------- inline functions ------------------------


namespace Physics
{
  namespace Notation
  {
    namespace Kelvin
    {
      namespace internal
      {
        /**
         * 返回与一个压缩组件索引相关的张量指数 <code><row,
         * column></code> 。 @p symmetric 标志表示 @p component_n
         * 指数是否与具有对称项的张量相关。
         *
         */
        template <int dim>
        std::pair<unsigned int, unsigned int>
        indices_from_component(const unsigned int component_n,
                               const bool         symmetric);


        template <int dim>
        std::pair<unsigned int, unsigned int>
        indices_from_component(const unsigned int component_n, const bool)
        {
          AssertThrow(false, ExcNotImplemented());
          return std::make_pair(0u, 0u);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<1>(const unsigned int component_n, const bool)
        {
          AssertIndexRange(component_n, 1);

          return std::make_pair(0u, 0u);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<2>(const unsigned int component_n,
                                  const bool         symmetric)
        {
          if (symmetric == true)
            {
              Assert(
                (component_n < SymmetricTensor<2, 2>::n_independent_components),
                ExcIndexRange(component_n,
                              0,
                              SymmetricTensor<2, 2>::n_independent_components));
            }
          else
            {
              Assert((component_n < Tensor<2, 2>::n_independent_components),
                     ExcIndexRange(component_n,
                                   0,
                                   Tensor<2, 2>::n_independent_components));
            }

          static const unsigned int indices[4][2] = {{0, 0},
                                                     {1, 1},
                                                     {0, 1},
                                                     {1, 0}};
          return std::make_pair(indices[component_n][0],
                                indices[component_n][1]);
        }


        template <>
        inline std::pair<unsigned int, unsigned int>
        indices_from_component<3>(const unsigned int component_n,
                                  const bool         symmetric)
        {
          if (symmetric == true)
            {
              Assert(
                (component_n < SymmetricTensor<2, 3>::n_independent_components),
                ExcIndexRange(component_n,
                              0,
                              SymmetricTensor<2, 3>::n_independent_components));
            }
          else
            {
              Assert((component_n < Tensor<2, 3>::n_independent_components),
                     ExcIndexRange(component_n,
                                   0,
                                   Tensor<2, 3>::n_independent_components));
            }

          static const unsigned int indices[9][2] = {{0, 0},
                                                     {1, 1},
                                                     {2, 2},
                                                     {1, 2},
                                                     {0, 2},
                                                     {0, 1},
                                                     {1, 0},
                                                     {2, 0},
                                                     {2, 1}};
          return std::make_pair(indices[component_n][0],
                                indices[component_n][1]);
        }


        /**
         * 返回应用于浓缩向量中的条目的缩放因子。
         *
         */
        template <int dim>
        double
        vector_component_factor(const unsigned int component_i,
                                const bool         symmetric)
        {
          if (symmetric == false)
            return 1.0;

          if (component_i < dim)
            return 1.0;
          else
            return numbers::SQRT2;
        }


        /**
         * 返回应用于浓缩矩阵中的条目的缩放系数。
         *
         */
        template <int dim>
        double
        matrix_component_factor(const unsigned int component_i,
                                const unsigned int component_j,
                                const bool         symmetric)
        {
          if (symmetric == false)
            return 1.0;

          // This case check returns equivalent of this result:
          // internal::vector_component_factor<dim>(component_i,symmetric)*internal::vector_component_factor<dim>(component_j,symmetric);
          if (component_i < dim && component_j < dim)
            return 1.0;
          else if (component_i >= dim && component_j >= dim)
            return 2.0;
          else // ((component_i >= dim && component_j < dim) || (component_i <
               // dim && component_j >= dim))
            return numbers::SQRT2;
        }

      } // namespace internal


      template <typename Number>
      Vector<Number>
      to_vector(const Number &s)
      {
        Vector<Number>     out(1);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          out(r) = s;
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<0, dim, Number> &s)
      {
        return to_vector(s.operator const Number &());
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<1, dim, Number> &v)
      {
        Vector<Number>     out(v.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;
            out(r)               = v[i];
          }
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const Tensor<2, dim, Number> &t)
      {
        Vector<Number>     out(t.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;
            out(r)               = t[i][j];
          }
        return out;
      }


      template <int dim, typename Number>
      Vector<Number>
      to_vector(const SymmetricTensor<2, dim, Number> &st)
      {
        Vector<Number>     out(st.n_independent_components);
        const unsigned int n_rows = out.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, true);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            Assert(indices.second >= indices.first, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;

            const double factor =
              internal::vector_component_factor<dim>(r, true);

            out(r) = factor * st[i][j];
          }
        return out;
      }


      template <typename Number>
      FullMatrix<Number>
      to_matrix(const Number &s)
      {
        FullMatrix<Number> out(1, 1);
        out(0, 0) = s;
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<0, dim, Number> &s)
      {
        return to_matrix(s.operator const Number &());
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<1, dim, Number> &v)
      {
        FullMatrix<Number> out(v.n_independent_components, 1);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                Assert(c < 1, ExcInternalError());
                out(r, c) = v[i];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<2, dim, Number> &t)
      {
        FullMatrix<Number> out(dim, dim);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_i =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_i.first < dim, ExcInternalError());
            Assert(indices_i.second < dim, ExcInternalError());
            const unsigned int i = indices_i.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_j =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_j.first < dim, ExcInternalError());
                Assert(indices_j.second < dim, ExcInternalError());
                const unsigned int j = indices_j.second;

                out(r, c) = t[i][j];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<2, dim, Number> &st)
      {
        return to_matrix(Tensor<2, dim, Number>(st));
      }


      namespace internal
      {
        template <typename TensorType>
        struct is_rank_2_symmetric_tensor : std::false_type
        {};

        template <int dim, typename Number>
        struct is_rank_2_symmetric_tensor<SymmetricTensor<2, dim, Number>>
          : std::true_type
        {};
      } // namespace internal


      template <int dim,
                typename SubTensor1,
                typename SubTensor2,
                typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<3, dim, Number> &t)
      {
        static_assert(
          (SubTensor1::dimension == dim && SubTensor2::dimension == dim),
          "Sub-tensor spatial dimension is different from those of the input tensor.");

        static_assert(
          (SubTensor1::rank == 2 && SubTensor2::rank == 1) ||
            (SubTensor1::rank == 1 && SubTensor2::rank == 2),
          "Cannot build a rank 3 tensor from the given combination of sub-tensors.");

        FullMatrix<Number> out(SubTensor1::n_independent_components,
                               SubTensor2::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();

        if (SubTensor1::rank == 2 && SubTensor2::rank == 1)
          {
            const bool subtensor_is_rank_2_symmetric_tensor =
              internal::is_rank_2_symmetric_tensor<SubTensor1>::value;

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_ij =
                  internal::indices_from_component<dim>(
                    r, subtensor_is_rank_2_symmetric_tensor);
                Assert(indices_ij.first < dim, ExcInternalError());
                Assert(indices_ij.second < dim, ExcInternalError());
                if (subtensor_is_rank_2_symmetric_tensor)
                  {
                    Assert(indices_ij.second >= indices_ij.first,
                           ExcInternalError());
                  }
                const unsigned int i = indices_ij.first;
                const unsigned int j = indices_ij.second;

                const double factor = internal::vector_component_factor<dim>(
                  r, subtensor_is_rank_2_symmetric_tensor);

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_k =
                      internal::indices_from_component<dim>(c, false);
                    Assert(indices_k.first < dim, ExcInternalError());
                    const unsigned int k = indices_k.first;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      out(r, c) = factor * t[i][j][k];
                    else
                      out(r, c) = t[i][j][k];
                  }
              }
          }
        else if (SubTensor1::rank == 1 && SubTensor2::rank == 2)
          {
            const bool subtensor_is_rank_2_symmetric_tensor =
              internal::is_rank_2_symmetric_tensor<SubTensor2>::value;

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_k =
                  internal::indices_from_component<dim>(r, false);
                Assert(indices_k.first < dim, ExcInternalError());
                const unsigned int k = indices_k.first;

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_ij =
                      internal::indices_from_component<dim>(
                        c, subtensor_is_rank_2_symmetric_tensor);
                    Assert(indices_ij.first < dim, ExcInternalError());
                    Assert(indices_ij.second < dim, ExcInternalError());
                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        Assert(indices_ij.second >= indices_ij.first,
                               ExcInternalError());
                      }
                    const unsigned int i = indices_ij.first;
                    const unsigned int j = indices_ij.second;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        const double factor =
                          internal::vector_component_factor<dim>(
                            c, subtensor_is_rank_2_symmetric_tensor);
                        out(r, c) = factor * t[k][i][j];
                      }
                    else
                      out(r, c) = t[k][i][j];
                  }
              }
          }
        else
          {
            AssertThrow(false, ExcNotImplemented());
          }

        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const Tensor<4, dim, Number> &t)
      {
        FullMatrix<Number> out(
          Tensor<2, dim, Number>::n_independent_components,
          Tensor<2, dim, Number>::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                out(r, c) = t[i][j][k][l];
              }
          }
        return out;
      }


      template <int dim, typename Number>
      FullMatrix<Number>
      to_matrix(const SymmetricTensor<4, dim, Number> &st)
      {
        FullMatrix<Number> out(
          SymmetricTensor<2, dim, Number>::n_independent_components,
          SymmetricTensor<2, dim, Number>::n_independent_components);
        const unsigned int n_rows = out.m();
        const unsigned int n_cols = out.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, true);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            Assert(indices_ij.second >= indices_ij.first, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, true);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                Assert(indices_kl.second >= indices_kl.first,
                       ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                const double factor =
                  internal::matrix_component_factor<dim>(r, c, true);

                out(r, c) = factor * st[i][j][k][l];
              }
          }
        return out;
      }


      template <typename Number>
      void
      to_tensor(const Vector<Number> &vec, Number &s)
      {
        Assert(vec.size() == 1, ExcDimensionMismatch(vec.size(), 1));
        s = vec(0);
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<0, dim, Number> &s)
      {
        return to_tensor(vec, s.operator Number &());
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<1, dim, Number> &v)
      {
        Assert(vec.size() == v.n_independent_components,
               ExcDimensionMismatch(vec.size(), v.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            const unsigned int i = indices.first;
            v[i]                 = vec(r);
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, Tensor<2, dim, Number> &t)
      {
        Assert(vec.size() == t.n_independent_components,
               ExcDimensionMismatch(vec.size(), t.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;
            t[i][j]              = vec(r);
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const Vector<Number> &vec, SymmetricTensor<2, dim, Number> &st)
      {
        Assert(vec.size() == st.n_independent_components,
               ExcDimensionMismatch(vec.size(), st.n_independent_components));
        const unsigned int n_rows = vec.size();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, true);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second < dim, ExcInternalError());
            Assert(indices.second >= indices.first, ExcInternalError());
            const unsigned int i = indices.first;
            const unsigned int j = indices.second;

            const double inv_factor =
              1.0 / internal::vector_component_factor<dim>(r, true);

            st[i][j] = inv_factor * vec(r);
          }
      }


      template <typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Number &s)
      {
        Assert(mtrx.m() == 1, ExcDimensionMismatch(mtrx.m(), 1));
        Assert(mtrx.n() == 1, ExcDimensionMismatch(mtrx.n(), 1));
        Assert(mtrx.n_elements() == 1,
               ExcDimensionMismatch(mtrx.n_elements(), 1));
        s = mtrx(0, 0);
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<0, dim, Number> &s)
      {
        return to_tensor(mtrx, s.operator Number &());
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<1, dim, Number> &v)
      {
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == 1, ExcDimensionMismatch(mtrx.n(), 1));
        Assert(mtrx.n_elements() == v.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    v.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices =
              internal::indices_from_component<dim>(r, false);
            Assert(indices.first < dim, ExcInternalError());
            Assert(indices.second == 0, ExcInternalError());
            const unsigned int i = indices.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                Assert(c < 1, ExcInternalError());
                v[i] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<2, dim, Number> &t)
      {
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == dim, ExcDimensionMismatch(mtrx.n(), dim));
        Assert(mtrx.n_elements() == t.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    t.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_i =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_i.first < dim, ExcInternalError());
            Assert(indices_i.second < dim, ExcInternalError());
            const unsigned int i = indices_i.first;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_j =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_j.first < dim, ExcInternalError());
                Assert(indices_j.second < dim, ExcInternalError());
                const unsigned int j = indices_j.second;

                t[i][j] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<2, dim, Number> &st)
      {
        // Its impossible to fit the (dim^2 + dim)/2 entries into a square
        // matrix We therefore assume that its been converted to a standard
        // tensor format using to_matrix (SymmetricTensor<2,dim,Number>) at some
        // point...
        Assert(mtrx.m() == dim, ExcDimensionMismatch(mtrx.m(), dim));
        Assert(mtrx.n() == dim, ExcDimensionMismatch(mtrx.n(), dim));
        Assert((mtrx.n_elements() ==
                Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n_elements(),
                 Tensor<2, dim, Number>::n_independent_components));

        Tensor<2, dim, Number> tmp;
        to_tensor(mtrx, tmp);
        st = symmetrize(tmp);
        Assert((Tensor<2, dim, Number>(st) - tmp).norm() < 1e-12,
               ExcMessage(
                 "The entries stored inside the matrix were not symmetric"));
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<3, dim, Number> &t)
      {
        Assert((mtrx.m() == Tensor<1, dim, Number>::n_independent_components) ||
                 (mtrx.m() ==
                  Tensor<2, dim, Number>::n_independent_components) ||
                 (mtrx.m() ==
                  SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcNotationExcFullMatrixToTensorColSize3(
                 mtrx.m(),
                 Tensor<1, dim, Number>::n_independent_components,
                 Tensor<2, dim, Number>::n_independent_components,
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() == Tensor<1, dim, Number>::n_independent_components) ||
                 (mtrx.n() ==
                  Tensor<2, dim, Number>::n_independent_components) ||
                 (mtrx.n() ==
                  SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcNotationExcFullMatrixToTensorColSize3(
                 mtrx.n(),
                 Tensor<1, dim, Number>::n_independent_components,
                 Tensor<2, dim, Number>::n_independent_components,
                 SymmetricTensor<2, dim, Number>::n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        if (mtrx.n() == Tensor<1, dim, Number>::n_independent_components)
          {
            Assert(
              (mtrx.m() == Tensor<2, dim, Number>::n_independent_components) ||
                (mtrx.m() ==
                 SymmetricTensor<2, dim, Number>::n_independent_components),
              ExcNotationExcFullMatrixToTensorRowSize2(
                mtrx.m(),
                Tensor<2, dim, Number>::n_independent_components,
                SymmetricTensor<2, dim, Number>::n_independent_components));

            const bool subtensor_is_rank_2_symmetric_tensor =
              (mtrx.m() ==
               SymmetricTensor<2, dim, Number>::n_independent_components);

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_ij =
                  internal::indices_from_component<dim>(
                    r, subtensor_is_rank_2_symmetric_tensor);
                Assert(indices_ij.first < dim, ExcInternalError());
                Assert(indices_ij.second < dim, ExcInternalError());
                if (subtensor_is_rank_2_symmetric_tensor)
                  {
                    Assert(indices_ij.second >= indices_ij.first,
                           ExcInternalError());
                  }
                const unsigned int i = indices_ij.first;
                const unsigned int j = indices_ij.second;

                const double inv_factor =
                  1.0 / internal::vector_component_factor<dim>(
                          r, subtensor_is_rank_2_symmetric_tensor);

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_k =
                      internal::indices_from_component<dim>(c, false);
                    Assert(indices_k.first < dim, ExcInternalError());
                    const unsigned int k = indices_k.first;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        t[i][j][k] = inv_factor * mtrx(r, c);
                        t[j][i][k] = t[i][j][k];
                      }
                    else
                      t[i][j][k] = mtrx(r, c);
                  }
              }
          }
        else
          {
            Assert(
              (mtrx.m() == Tensor<1, dim, Number>::n_independent_components),
              ExcDimensionMismatch(
                mtrx.m(), Tensor<1, dim, Number>::n_independent_components));
            Assert(
              (mtrx.n() == Tensor<2, dim, Number>::n_independent_components) ||
                (mtrx.n() ==
                 SymmetricTensor<2, dim, Number>::n_independent_components),
              ExcNotationExcFullMatrixToTensorColSize2(
                mtrx.n(),
                Tensor<2, dim, Number>::n_independent_components,
                SymmetricTensor<2, dim, Number>::n_independent_components));

            const bool subtensor_is_rank_2_symmetric_tensor =
              (mtrx.n() ==
               SymmetricTensor<2, dim, Number>::n_independent_components);

            for (unsigned int r = 0; r < n_rows; ++r)
              {
                const std::pair<unsigned int, unsigned int> indices_k =
                  internal::indices_from_component<dim>(r, false);
                Assert(indices_k.first < dim, ExcInternalError());
                const unsigned int k = indices_k.first;

                for (unsigned int c = 0; c < n_cols; ++c)
                  {
                    const std::pair<unsigned int, unsigned int> indices_ij =
                      internal::indices_from_component<dim>(
                        c, subtensor_is_rank_2_symmetric_tensor);
                    Assert(indices_ij.first < dim, ExcInternalError());
                    Assert(indices_ij.second < dim, ExcInternalError());
                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        Assert(indices_ij.second >= indices_ij.first,
                               ExcInternalError());
                      }
                    const unsigned int i = indices_ij.first;
                    const unsigned int j = indices_ij.second;

                    if (subtensor_is_rank_2_symmetric_tensor)
                      {
                        const double inv_factor =
                          1.0 / internal::vector_component_factor<dim>(
                                  c, subtensor_is_rank_2_symmetric_tensor);
                        t[k][i][j] = inv_factor * mtrx(r, c);
                        t[k][j][i] = t[k][i][j];
                      }
                    else
                      t[k][i][j] = mtrx(r, c);
                  }
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &mtrx, Tensor<4, dim, Number> &t)
      {
        Assert((mtrx.m() == Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.m(), Tensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() == Tensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n(), Tensor<2, dim, Number>::n_independent_components));
        Assert(mtrx.n_elements() == t.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    t.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                t[i][j][k][l] = mtrx(r, c);
              }
          }
      }


      template <int dim, typename Number>
      void
      to_tensor(const FullMatrix<Number> &       mtrx,
                SymmetricTensor<4, dim, Number> &st)
      {
        Assert((mtrx.m() ==
                SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.m(),
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert((mtrx.n() ==
                SymmetricTensor<2, dim, Number>::n_independent_components),
               ExcDimensionMismatch(
                 mtrx.n(),
                 SymmetricTensor<2, dim, Number>::n_independent_components));
        Assert(mtrx.n_elements() == st.n_independent_components,
               ExcDimensionMismatch(mtrx.n_elements(),
                                    st.n_independent_components));

        const unsigned int n_rows = mtrx.m();
        const unsigned int n_cols = mtrx.n();
        for (unsigned int r = 0; r < n_rows; ++r)
          {
            const std::pair<unsigned int, unsigned int> indices_ij =
              internal::indices_from_component<dim>(r, false);
            Assert(indices_ij.first < dim, ExcInternalError());
            Assert(indices_ij.second < dim, ExcInternalError());
            const unsigned int i = indices_ij.first;
            const unsigned int j = indices_ij.second;

            for (unsigned int c = 0; c < n_cols; ++c)
              {
                const std::pair<unsigned int, unsigned int> indices_kl =
                  internal::indices_from_component<dim>(c, false);
                Assert(indices_kl.first < dim, ExcInternalError());
                Assert(indices_kl.second < dim, ExcInternalError());
                const unsigned int k = indices_kl.first;
                const unsigned int l = indices_kl.second;

                const double inv_factor =
                  1.0 / internal::matrix_component_factor<dim>(r, c, true);

                st[i][j][k][l] = inv_factor * mtrx(r, c);
              }
          }
      }


      template <typename TensorType, typename Number>
      inline TensorType
      to_tensor(const Vector<Number> &vec)
      {
        TensorType out;
        to_tensor(vec, out);
        return out;
      }


      template <typename TensorType, typename Number>
      inline TensorType
      to_tensor(const FullMatrix<Number> &mtrx)
      {
        TensorType out;
        to_tensor(mtrx, out);
        return out;
      }

    } // namespace Kelvin
  }   // namespace Notation
} // namespace Physics


#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif


