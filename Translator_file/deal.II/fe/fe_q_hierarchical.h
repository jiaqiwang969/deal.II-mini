//include/deal.II-translator/fe/fe_q_hierarchical_0.txt
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

#ifndef dealii_fe_q_hierarchical_h
#define dealii_fe_q_hierarchical_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 分层 @p Qp 形状函数的实现，产生连续的、程度为 @p p.
 * 的分片多项式的有限元空间。该类使用基于区间<tt>[0,1]</tt>上的分层基础
 * Polynomials::Hierarchical
 * 的张量积多项式实现，如果我们假设每个元素有一个单度，该张量积有限元适合构建一个
 * @p hp 张量积有限元。 该类的构造函数取该有限元的度数
 * @p p 。 该类没有实现对一维的情况（<tt>spacedim !=
 * dim</tt>）。 <h3>Implementation</h3>
 * 构造函数创建一个TensorProductPolynomials对象，其中包括度数为
 * @p Hierarchical 的多项式的张量积 @p p.  这个 @p
 * TensorProductPolynomials对象提供形状函数的所有值和导数。
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 * TensorProductPolynomials所代表的形状函数的原始排序是张量乘法的编号。然而，单元格上的形状函数被重新编号，从支持点在顶点的形状函数开始，然后是直线上的，四边形上的，最后是（对于三维）六边形上的。为了明确起见，这些编号列在下面。
 * <h4>Q1 elements</h4> $Q_1^H$
 * 元素是多项式一度，因此，与FE_Q类中的 $Q_1$
 * 元素完全相同。特别是，形状函数的定义完全相同。
 * <ul>   <li>  1D情况。
 * @verbatim
 *    0-------1
 * @endverbatim
 * <li>  二维情况。
 * @verbatim
 *    2-------3
 *    |       |
 *    |       |
 *    |       |
 *    0-------1
 * @endverbatim
 * <li>  3D情况。
 * @verbatim
 *       6-------7        6-------7
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   4   |       |    4-------5   |
 *   |   2-------3    |       |   3
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 *
 * 自由度的支持点的各自坐标值如下。  <ul>   <li>  形状函数0。<tt>[0, 0, 0]</tt>;  <li>  形状函数1：<tt>[1, 0, 0]</tt>;  <li>  形状函数2：<tt>[0, 1, 0]</tt>;  <li>  形状函数3：<tt>[1, 1, 0]</tt>;  <li>  ] 形状函数4：<tt>[0, 0, 1]</tt>;  <li>  形状函数5：<tt>[1, 0, 1]</tt>;  <li>  形状函数6：<tt>[0, 1, 1]</tt>;  <li>  形状函数7：<tt>[1, 1, 1]</tt>;  </ul>   </ul>
 * 在2d中，这些形状函数看起来如下。  <table> <tr> <td
 * align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_1^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_1^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q1/Q1H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_1^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_1^H$ element, shape function 3 </td> </tr> </table>
 *
 * <h4>Q2 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0---2---1
 * @endverbatim
 * <li>  二维情况。
 * @verbatim
 *    2---7---3
 *    |       |
 *    4   8   5
 *    |       |
 *    0---6---1
 * @endverbatim
 * <li>  3D情况。
 * @verbatim
 *       6--15---7        6--15---7
 *      /|       |       /       /|
 *    12 |       19     12      1319
 *    /  18      |     /       /  |
 *   4   |       |    4---14--5   |
 *   |   2---11--3    |       |   3
 *   |  /       /     |      17  /
 *  16 8       9     16       | 9
 *   |/       /       |       |/
 *   0---10--1        0---8---1
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *      /|       |       /       /|
 *     / |  23   |      /  25   / |
 *    /  |       |     /       /  |
 *     |       |
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*   |
 *   |20-------*    |       |21
 *   |  /       /     |   22  |  /
 *   | /  24   /      |       | /
 *   |/       /       |       |/
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 * -------*
 * @endverbatim
 * 中心顶点的编号为26。
 * 自由度的支持点的各自坐标值如下。  <ul>   <li>  形状函数0。<tt>[0, 0, 0]</tt>;  <li>  形状函数1：<tt>[1, 0, 0]</tt>;  <li>  形状函数2：<tt>[0, 1, 0]</tt>;  <li>  形状函数3：<tt>[1, 1, 0]</tt>;  <li>  ] 形状函数4：<tt>[0，0，1]</tt>;  <li>  形状函数5：<tt>[1，0，1]</tt>;  <li>  形状函数6：<tt>[0，1，1]</tt>;  <li>  ] 形状函数7：<tt>[1，1，1]</tt>;  <li>  形状函数8：<tt>[0，1/2，0]</tt>;  <li>  形状函数9：<tt>[1，1/2，0]</tt>;  <li>  形状函数10：<tt>[1/2，0，0]</tt>;  <li>  ] 形状函数11：<tt>[1/2, 1, 0]/tt>;  <li>  形状函数12：<tt>[0, 1/2, 1]/tt>;  <li>  形状函数13：<tt>[1, 1/2, 1]/tt>;  <li>  ] 形状函数14：<tt>[1/2, 0, 1]</tt>;  <li>  形状函数15：<tt>[1/2, 1, 1]</tt>;  <li>  形状函数16：<tt>[0, 0, 1/2]</tt>;  <li>  形状函数17：<tt>[1, 0, 1/2]</tt>;  <li>  ] 形状函数18：<tt>[0, 1, 1/2]/tt>;  <li>  形状函数19：<tt>[1, 1, 1/2]/tt>;  <li>  形状函数20：<tt>[0, 1/2, 1/2]/tt>;  <li>  ] 形状函数21：<tt>[1, 1/2, 1/2]/tt>;  <li>  形状函数22：<tt>[1/2, 0, 1/2]/tt>;  <li>  形状函数23：<tt>[1/2, 1/2]/tt>;  <li>  ] 形状函数24：<tt>[1/2, 1/2, 0]</tt>;  <li>  形状函数25：<tt>[1/2, 1/2, 1]</tt>;  <li>  形状函数26：<tt>[1/2, 1/2]</tt>;  </ul>   </ul>
 *
 * 在2d中，这些形状函数看起来如下（黑色平面对应于零；负的形状函数值可能不可见）。
 * <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_2^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_2^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q2/Q2H_shape0008.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_2^H$ element,
 * shape function 8 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 * <h4>Q3 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0--2--3--1
 * @endverbatim
 *
 * <li>  2D情况。
 * @verbatim
 *    2--10-11-3
 *    |        |
 *    5  14 15 7
 *    |        |
 *    4  12 13 6
 *    |        |
 *    0--8--9--1
 * @endverbatim
 * </ul>
 * 在2D中，这些形状函数看起来如下（黑色平面对应于零；负的形状函数值可能不可见）。
 * <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 8
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 10
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 12
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0014.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q3/Q3H_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_3^H$ element, shape function 14
 * </td>
 *
 * <td align="center"> $Q_3^H$ element, shape function 15 </td> </tr> </table>
 *
 * <h4>Q4 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0--2--3--4--1
 * @endverbatim
 * <li>  二维情况。
 * @verbatim
 *    2--13-14-15-3
 *    |           |
 *    6  22 23 24 9
 *    |           |
 *    5  19 20 21 8
 *    |           |
 *    4  16 17 18 7
 *    |           |
 *    0--10-11-12-1
 * @endverbatim
 * </ul>
 * 在2D中，这些形状函数看起来如下（黑色平面对应于零；负的形状函数值可能不可见）。
 * <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 0
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 2
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 4
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 6
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 8
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 10
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 12
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0014.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 14
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 15 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0016.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0017.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 16
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 17 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0018.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0019.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 18
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 19 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0020.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0021.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 20
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 21 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0022.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0023.png
 * </td> </tr> <tr> <td align="center"> $Q_4^H$ element, shape function 22
 * </td>
 *
 * <td align="center"> $Q_4^H$ element, shape function 23 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/hierarchical/Q4/Q4H_shape0024.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_4^H$ element,
 * shape function 24 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 *
 */
template <int dim>
class FE_Q_Hierarchical : public FE_Poly<dim>
{
public:
  /**
   * 度数为 @p p. 的张量乘积多项式的构造函数
   *
   */
  FE_Q_Hierarchical(const unsigned int p);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Q_Hierarchical<dim>(degree)</tt>，其中
   * @p dim 和 @p 度数由适当的值替换。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, dim>>
  clone() const override;

  /**
   * 如果形状函数 @p shape_index 在面 @p face_index.
   * 的某处有非零函数值，该函数返回 @p true, 。
   *
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * @name  支持hp的函数  
     * @{ 
   *
   */

  /**
   * 返回该元素是否以新的方式实现其悬挂节点约束，这必须用于使元素
   * "hp-compatible"。
   * 对于FE_Q_Hierarchical类，结果总是为真（与元素的程度无关），因为它实现了hp-capability所需的完整功能集。
   *
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * 返回从给定的有限元内插到现在的矩阵。只支持FE_Q_Hierarchical之间的插值。
   *
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim> &source,
                           FullMatrix<double> &      matrix) const override;

  /**
   * 网格之间的嵌入矩阵。只支持各向同性的细化。
   *
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * 如果在一个顶点上，有几个有限元被激活，hp-code首先为这些FEs的自由度分配不同的全局索引。然后调用这个函数来找出其中哪些应该得到相同的值，从而可以得到相同的全局自由度指数。
   * 因此，该函数返回当前有限元对象的自由度与 @p fe_other,
   * 的自由度之间的相同性列表，后者是对代表在该特定顶点上活动的其他有限元之一的有限元对象的引用。该函数计算两个有限元对象的哪些自由度是相等的，这两个自由度的编号都在零和两个有限元的n_dofs_per_vertex()的相应值之间。每一对的第一个索引表示本元素的一个顶点自由度，而第二个是另一个有限元素的相应索引。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 与上述相同，但对线而言。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(const FiniteElement<dim> &fe_other) const override;

  /**
   * 同上，但对面而言。
   *
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                         const unsigned int        face_no = 0) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination()
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim> &fe_other,
                         const unsigned int codim = 0) const override final;

   /*@}*/ 

  /**
   * 返回从一个元素的面插值到邻近元素的面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现从一个给定元素的内插，那么他们必须抛出一个类型为
   * <tt>FiniteElement<dim>::ExcInterpolationNotImplemented</tt>. 的异常。
   *
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim> &source,
                                FullMatrix<double> &      matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * 返回从一个元素的面内插到邻近元素的子面的矩阵。矩阵的大小是<tt>source.dofs_per_face</tt>乘以<tt>this->dofs_per_face</tt>。
   * 衍生元素将不得不实现这个函数。他们可能只为某些源有限元提供插值矩阵，例如那些来自同一家族的有限元。如果他们不实现给定元素的插值，那么他们必须抛出一个<tt>ExcInterpolationNotImplemented</tt>类型的异常。
   *
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim> &source,
    const unsigned int        subface,
    FullMatrix<double> &      matrix,
    const unsigned int        face_no = 0) const override;

  /**
   * 确定此对象的内存消耗（以字节为单位）的估计值。
   * 这个函数是虚拟的，因为有限元对象通常是通过指向其基类的指针来访问的，而不是类本身。
   *
   */
  virtual std::size_t
  memory_consumption() const override;

  /**
   * 对于度数为 @p sub_degree < @p degree,
   * 的有限元，我们返回一个向量，将度数为 @p sub_degree
   * 的FE上的编号映射到这个元素上的编号。
   *
   */
  std::vector<unsigned int>
  get_embedding_dofs(const unsigned int sub_degree) const;

  /**
   * 返回一个元素的常数模式的列表。对于这个元素，该列表由第一个顶点形状函数的真参数和其余的假参数组成。
   *
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

private:
  /**
   * 仅供内部使用。它的全称是 @p get_dofs_per_object_vector
   * 函数，它创建了 @p dofs_per_object
   * 向量，在构造函数内需要传递给 @p
   * FiniteElementData的构造函数。
   *
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * 连续有限元中自由度的编号是分层次的，也就是说，我们首先按照三角形定义的顶点顺序对顶点自由度进行编号，然后按照顺序和尊重线的方向对线自由度进行编号，然后是四边形上的自由度，等等。
   * 与1d分层多项式相关的道夫是以顶点为先（ $phi_0(x)=1-x$
   * 和 $phi_1(x)=x$ ），然后是线道夫（高阶多项式）。
   * 2d和3d等级多项式通过张量乘积来源于1d等级多项式。在下文中，由此产生的dofs的编号将用`fe_q_hierarchical
   * numbering`来表示。
   * 这个函数构建了一个表格，其中fe_q_hierarchical索引中的每个自由度在分级编号中会有。
   * 这个函数类似于 FETools::hierarchic_to_lexicographic_numbering()
   * 函数。然而，与上面定义的fe_q_hierarchical编号不同，lexicographic编号源于连续编号的自由度的张量积（就像对LagrangeEquidistant）。
   * 假设输出参数的大小已经符合正确的大小，等于有限元的自由度数。
   *
   */
  static std::vector<unsigned int>
  hierarchic_to_fe_q_hierarchical_numbering(const FiniteElementData<dim> &fe);

  /**
   * 这是一个类似于前一个函数的函数，但对面进行工作。
   *
   */
  static std::vector<unsigned int>
  face_fe_q_hierarchical_to_hierarchic_numbering(const unsigned int degree);

  /**
   * 初始化两个辅助字段，用于设置构造函数中的各种矩阵。
   *
   */
  void
  build_dofs_cell(std::vector<FullMatrix<double>> &dofs_cell,
                  std::vector<FullMatrix<double>> &dofs_subcell) const;

  /**
   * 初始化悬挂节点的约束矩阵。从构造函数中调用。
   *
   */
  void
  initialize_constraints(const std::vector<FullMatrix<double>> &dofs_subcell);

  /**
   * 初始化嵌入矩阵。从构造函数中调用。
   *
   */
  void
  initialize_embedding_and_restriction(
    const std::vector<FullMatrix<double>> &dofs_cell,
    const std::vector<FullMatrix<double>> &dofs_subcell);

  /**
   * 初始化FiniteElement类的 @p generalized_support_points 字段。
   * 从构造函数中调用。
   *
   */
  void
  initialize_generalized_support_points();

  /**
   * 初始化FiniteElement类的 @p generalized_face_support_points
   * 字段。从构造函数中调用。
   *
   */
  void
  initialize_generalized_face_support_points();

  /**
   * 在第一个面上从lexicographic到形状函数编号的映射。
   *
   */
  const std::vector<unsigned int> face_renumber;

  // Allow access from other dimensions. We need this since we want to call
  // the functions @p get_dpo_vector and @p
  // lexicographic_to_hierarchic_numbering for the faces of the finite element
  // of dimension dim+1.
  template <int dim1>
  friend class FE_Q_Hierarchical;
};

 /*@}*/ 

 /* -------------- declaration of explicit specializations ------------- */ 

template <>
void
FE_Q_Hierarchical<1>::initialize_generalized_face_support_points();

template <>
bool
FE_Q_Hierarchical<1>::has_support_on_face(const unsigned int,
                                          const unsigned int) const;

template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_fe_q_hierarchical_to_hierarchic_numbering(
  const unsigned int);

DEAL_II_NAMESPACE_CLOSE

#endif


