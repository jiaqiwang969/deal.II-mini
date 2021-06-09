//include/deal.II-translator/fe/fe_q_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2021 by the deal.II authors
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

#ifndef dealii_fe_q_h
#define dealii_fe_q_h

#include <deal.II/base/config.h>

#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_q_base.h>

DEAL_II_NAMESPACE_OPEN


 /*!@addtogroup fe */ 
 /*@{*/ 

/**
 * 标量拉格朗日有限元 @p Qp
 * 的实现，得到每个坐标方向上连续的、程度为 @p p
 * 的分片多项式的有限元空间。该类使用基于1D
 * Lagrange多项式的张量积多项式实现，具有等距（2度以内）、Gauss-Lobatto（从3度开始）或给定支持点。
 * 该类的标准构造函数取该有限元的度数 @p p
 * 。或者，它可以取一个正交公式 @p points
 * ，定义一个坐标方向上拉格朗日插值的支持点。
 * 关于<tt>spacedim</tt>模板参数的更多信息，请查阅FiniteElement或Triangulation的文档。
 * <h3>Implementation</h3>
 * 构造函数创建一个TensorProductPolynomials对象，其中包括度数为
 * @p LagrangeEquidistant 的多项式的张量积 @p p. 这个 @p
 * TensorProductPolynomials 对象提供形状函数的所有值和导数。
 * 在给出正交规则的情况下，构造函数创建一个TensorProductPolynomials对象，其中包括
 * @p Lagrange 多项式与 @p points. 的支持点的张量乘积。
 * 此外，构造函数还填充了 @p interface_constraints, 、 @p
 * 的延长（嵌入）和 @p restriction
 * 矩阵。这些只在一定程度上实现，对于非常高的多项式程度可能无法使用。
 * <h3>Unit support point distribution and conditioning of interpolation</h3>
 * 在构建多项式度数为1或2的FE_Q元素时，使用0和1（线性情况）或0、0.5和1（二次情况）的等距支持点。单位支持点或节点点<i>x<sub>i</sub></i>是那些<i>j</i>个拉格朗日多项式满足
 * $\delta_{ij}$ 属性的点，即一个多项式为1，其他都是0。
 * 对于更高的多项式度数，支持点默认为非流动性的，并选择为<tt>（度数+1）</tt>阶Gauss-Lobatto正交规则的支持点。这种点分布在任意多项式度数下产生条件良好的Lagrange插值。相比之下，基于等距点的多项式随着多项式度数的增加而变得越来越没有条件。在内插法中，这种效应被称为Runge现象。对于Galerkin方法，Runge现象通常在解的质量上不明显，而是在相关系统矩阵的条件数上。例如，10度的等距点的元素质量矩阵的条件数为2.6e6，而Gauss-Lobatto点的条件数约为400。
 * 一维的Gauss-Lobatto点包括单位区间的端点0和+1。内部点被移向端点，这使得靠近元素边界的点分布更加密集。
 * 如果与Gauss-Lobatto正交相结合，基于默认支持点的FE_Q可以得到对角线的质量矩阵。这种情况在
 * step-48
 * 中得到了证明。然而，这个元素可以通过通常的FEValues方法与任意的正交规则相结合，包括全高斯正交。在一般情况下，质量矩阵是非对角线的。
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 * TensorProductPolynomials所代表的形状函数的原始排序是张量乘法的编号。然而，单元格上的形状函数被重新编号，从支持点在顶点的形状函数开始，然后是直线上的，四边形上的，最后是（对于三维）六边形上的。为了明确起见，这些编号列在下面。
 * <h4>Q1 elements</h4>  <ul>   <li>  1D情况。
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
 * 形状函数的支持点的各自坐标值如下。  <ul>   <li>  形状函数0。<tt>[0, 0, 0]</tt>;  <li>  形状函数1：<tt>[1, 0, 0]</tt>;  <li>  形状函数2：<tt>[0, 1, 0]</tt>;  <li>  形状函数3：<tt>[1, 1, 0]</tt>;  <li>  ] 形状函数4：<tt>[0, 0, 1]</tt>;  <li>  形状函数5：<tt>[1, 0, 1]</tt>;  <li>  形状函数6：<tt>[0, 1, 1]</tt>;  <li>  形状函数7：<tt>[1, 1, 1]</tt>;  </ul>   </ul>
 * 在2d中，这些形状函数看起来如下。  <table> <tr> <td
 * align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_1$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_1$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q1/Q1_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_1$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_1$ element, shape function 3 </td> </tr> </table>
 *
 * <h4>Q2 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0---2---1
 * @endverbatim
 *
 * <li>  2D情况。
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
 *   0---10--1        0---10--1
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
 * 形状函数的支持点的各自坐标值如下。  <ul>   <li>  形状函数0。<tt>[0, 0, 0]</tt>;  <li>  形状函数1：<tt>[1, 0, 0]</tt>;  <li>  形状函数2：<tt>[0, 1, 0]</tt>;  <li>  形状函数3：<tt>[1, 1, 0]</tt>;  <li>  ] 形状函数4：<tt>[0，0，1]</tt>;  <li>  形状函数5：<tt>[1，0，1]</tt>;  <li>  形状函数6：<tt>[0，1，1]</tt>;  <li>  ] 形状函数7：<tt>[1，1，1]</tt>;  <li>  形状函数8：<tt>[0，1/2，0]</tt>;  <li>  形状函数9：<tt>[1，1/2，0]</tt>;  <li>  形状函数10：<tt>[1/2，0，0]</tt>;  <li>  ] 形状函数11：<tt>[1/2, 1, 0]/tt>;  <li>  形状函数12：<tt>[0, 1/2, 1]/tt>;  <li>  形状函数13：<tt>[1, 1/2, 1]/tt>;  <li>  ] 形状函数14：<tt>[1/2, 0, 1]</tt>;  <li>  形状函数15：<tt>[1/2, 1, 1]</tt>;  <li>  形状函数16：<tt>[0, 0, 1/2]</tt>;  <li>  形状函数17：<tt>[1, 0, 1/2]</tt>;  <li>  ] 形状函数18：<tt>[0, 1, 1/2]/tt>;  <li>  形状函数19：<tt>[1, 1, 1/2]/tt>;  <li>  形状函数20：<tt>[0, 1/2, 1/2]/tt>;  <li>  ] 形状函数21：<tt>[1, 1/2, 1/2]/tt>;  <li>  形状函数22：<tt>[1/2, 0, 1/2]/tt>;  <li>  形状函数23：<tt>[1/2, 1/2]/tt>;  <li>  ] 形状函数24：<tt>[1/2, 1/2, 0]</tt>;  <li>  形状函数25：<tt>[1/2, 1/2, 1]</tt>;  <li>  形状函数26：<tt>[1/2, 1/2]</tt>;  </ul>   </ul>
 *
 * 在2d中，这些形状函数看起来如下（黑色平面对应于零；负的形状函数值可能不可见）。
 * <table> <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_2$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_2$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q2/Q2_shape0008.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_2$ element,
 * shape function 8 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 * <h4>Q3 elements</h4>  <ul>   <li>  1D情况。
 * @verbatim
 *    0--2--3--1
 * @endverbatim
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
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 8 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 10 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 12 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0014.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q3/Q3_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_3$ element, shape function 14 </td>
 *
 * <td align="center"> $Q_3$ element, shape function 15 </td> </tr> </table>
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
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0000.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0001.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 0 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 1 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0002.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0003.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 2 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 3 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0004.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0005.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 4 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 5 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0006.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0007.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 6 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 7 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0008.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0009.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 8 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 9 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0010.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0011.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 10 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 11 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0012.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0013.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 12 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 13 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0014.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0015.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 14 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 15 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0016.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0017.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 16 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 17 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0018.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0019.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 18 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 19 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0020.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0021.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 20 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 21 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0022.png
 * </td>
 *
 * <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0023.png
 * </td> </tr> <tr> <td align="center"> $Q_4$ element, shape function 22 </td>
 *
 * <td align="center"> $Q_4$ element, shape function 23 </td> </tr>
 *
 * <tr> <td align="center">
 @image html http://www.dealii.org/images/shape-functions/lagrange/Q4/Q4_shape0024.png
 * </td>
 *
 * <td align="center"> </td> </tr> <tr> <td align="center"> $Q_4$ element,
 * shape function 24 </td>
 *
 * <td align="center"> </td> </tr> </table>
 *
 *
 */
template <int dim, int spacedim = dim>
class FE_Q : public FE_Q_Base<dim, spacedim>
{
public:
  /**
   * 基于Gauss-Lobatto支持（节点）点的度数 @p p
   * 的张量乘积多项式的构造器。对于度数为1和2的多项式，这些是通常的等距点。
   *
   */
  FE_Q(const unsigned int p);

  /**
   * 基于一维正交公式的支持点 @p points
   * 的张量乘积多项式的构造器。有限元的程度是<tt>points.size()-1</tt>。注意，第一个点必须是0，最后一个是1。构建<tt>FE_Q<dim>(QGaussLobatto<1>(fe_degree+1))</tt>等同于只指定多项式程度的构建器。对于选择<tt>fe_degree
   * >
   * 2</tt>的等距节点，构造<tt>FE_Q<dim>(QIterated<1>(QTrapezoid<1>(),fe_degree))</tt>。
   * 这个构造函数所创建的空间*与你调用`FE_Q<dim>(point.size()-1)'的情况相同，但*不同的是这个空间的基函数。
   * 不同的是这个空间的基函数。这在一些情况下是很有用的，人们希望通过整合这些形状函数的双线性形式来实现矩阵的某些属性。  例如，当计算单元格 $K$ 、@f[
   * M_{ij}^K = \int_K \varphi_i(\mathbf x) \varphi_j(\mathbf x) \; dx
   * @f]上的质量矩阵时，人们通常应用正交公式，并通过以下方式近似真实质量矩阵。  @f[
   * M_{ij}^K = \sum_q \varphi_i(\mathbf x_q) \varphi_j(\mathbf x_q) w_q,
   * @f] 其中正交点 $\mathbf x_q$ 和权重 $w_q$ 的位置取决于单元格 $K$  。如果用于定义这些点的正交公式 $\mathbf x_q$ 与传递给这个构造器的公式相同（或者，在更高维度上，由用于构造器的张量积产生），那么 $\varphi_i(\mathbf x_q) = \delta_{iq}$  ，矩阵还原为@f[
   * M_{ij}^K = \sum_q \delta_{iq} \delta_{jq} w_q = \delta_{ij} w_i,
   * @f]，即，对角线上有权重 $w_i$
   * 的斜向矩阵。这样的结构在使用显式时间步进方法时非常有用，因为，例如，在解决线性系统时，人们只需要在每个时间步进中反转对角线质量矩阵。
   *
   */
  FE_Q(const Quadrature<1> &points);

  /**
   * 返回一个唯一标识有限元的字符串。该类返回<tt>FE_Q<dim>(degree)</tt>，
   * @p dim 和 @p degree 用适当的值代替。
   *
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  /**
   * 在FiniteElement类中实现相应的函数。
   * 由于当前元素是插值的，所以节点值正好是支持点的值。此外，由于当前元素是标量的，支持点的值需要是长度为1的向量。
   *
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              nodal_values) const override;

  /**
   * @copydoc   FiniteElement::compare_for_domination() 。
   *
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;
};



 /*@}*/ 

DEAL_II_NAMESPACE_CLOSE

#endif


