//include/deal.II-translator/integrators/local_integrators_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2018 by the deal.II authors
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

#ifndef dealii_integrators_local_integrators_h
#define dealii_integrators_local_integrators_h

// This file only provides definition and documentation of the
// namespace LocalIntegrators. There is no necessity to include it
// anywhere in a C++ code. Only doxygen will make use of it.

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN


/**
 *
 * @brief Library of integrals over cells and faces
 * 这个命名空间包含了双线性形式、形式和误差估计的特定应用局部积分。它是一个函数的集合，被组织成专门用于某些应用的命名空间。例如，命名空间Laplace包含计算拉普拉斯算子的单元矩阵和单元残差的函数，以及Nitsche的弱边界条件或内部惩罚不连续Galerkin方法的函数。麦克斯韦命名空间对卷曲型问题也有同样的作用。
 * L2命名空间包含质量矩阵和<i>L<sup>2</sup></i>内积的函数。
 * <h3>Notational conventions</h3>
 * 在大多数情况下，这个命名空间中的函数的作用可以用一个单一的积分来描述。我们对单元格<i>Z</i>上的积分和面<i>F</i>上的积分加以区分。如果一个积分被表示为@f[
 * \int_Z u \otimes v \,dx,
 * @f]，它将产生以下结果，取决于操作的类型 <ul>   <li>  如果函数返回一个矩阵，位置<i>(i,j)</i>的条目将是测试函数<i>v<sub>i</sub></i> ]和试验函数<i>u<sub>j</sub></i>的综合积（注意指数的还原）  </li>   <li>  如果函数返回一个向量，那么在位置<i>i</i>的向量条目将是给定函数<i>u</i>与试验函数<i>v<sub>i</sub></i>的综合积。 </li>   <li>  如果函数返回一个数字，那么这个数字就是两个给定函数<i>u</i>和<i>v</i>的积分。  </ul>
 * 我们将使用常规的草书符号 $u$ 来表示标量，用粗体符号
 * $\mathbf u$
 * 表示向量。测试函数总是<i>v</i>，试用函数总是<i>u</i>。参数为希腊语，面法向量为
 * $\mathbf n = \mathbf n_1 =
 *
 * -\mathbf n_2$  。 <h3>Signature of functions</h3>
 * 这个命名空间的函数遵循一个通用的签名。在最简单的情况下，你有两个相关的函数
 *
 * @code
 * template <int dim>
 * void
 * cell_matrix (
 *   FullMatrix<double>& M,
 *   const FEValuesBase<dim>& fe,
 *   const double factor = 1.);
 *
 * template <int dim>
 * void
 * cell_residual (
 *   BlockVector<double>* v,
 *   const FEValuesBase<dim>& fe,
 *   const std::vector<Tensor<1,dim> >& input,
 *   const double factor = 1.);
 * @endcode
 *
 * 对于同一个算子通常有一对函数，函数<tt>cell_residual</tt>实现算子从有限元空间到其对偶的映射，函数<tt>cell_matrix</tt>生成对应于<tt>cell_residual</tt>的Frechet导数的双线性形式。
 * 这些函数的第一个参数是返回类型，它是 <ul>   <li>  FullMatrix&lt;double&gt; 用于矩阵  <li>  BlockVector&ltdouble&gt; 用于矢量  </ul>  。
 * 下一个参数是代表用于积分的有限元的FEValuesBase对象。如果积分运算符从一个有限元空间映射到另一个有限元空间的对偶（例如块系统中的对角线矩阵），那么首先要指定试验空间的FEValuesBase，之后是测试空间的FEValuesBase。
 * 这个列表后面是所需的数据集，顺序是 <ol>   <li>  来自有限元函数的数据向量  <li>  来自其他对象的数据向量  <li>  附加数据  <li>  一个与整个结果相乘的系数  </ol>  。
 * <h3>Usage</h3>
 * 本地积分器可以用在任何本来要实现本地积分循环的地方。下面的例子来自斯托克斯求解器的实现，使用 MeshWorker::Assembler::LocalBlocksToGlobalBlocks.  矩阵是  <ul>   <li>  0：速度的矢量拉普拉斯（这里有一个矢量元素）  <li>  1：发散矩阵  <li>  2：预处理器中使用的压力质量矩阵  </ul>  。
 * 有了这些矩阵， MeshWorker::loop() 调用的函数可以写成这样
 *
 * @code
 * using namespace dealii::LocalIntegrators;
 *
 * template <int dim>
 * void MatrixIntegrator<dim>::cell(MeshWorker::DoFInfo<dim>         &dinfo,
 *                                MeshWorker::IntegrationInfo<dim> &info)
 * {
 * Laplace::cell_matrix (dinfo.matrix(0,false).matrix,
 *                       info.fe_values(0));
 * Divergence::cell_matrix (dinfo.matrix(1,false).matrix,
 *                          info.fe_values(0),
 *                          info.fe_values(1));
 * L2::cell_matrix (dinfo.matrix(2,false).matrix,
 *                  info.fe_values(1));
 * }
 * @endcode
 * 参见 step-39 ，了解该代码的工作实例。
 *
 *
 * @ingroup Integrators
 *
 *
 */
namespace LocalIntegrators
{}

DEAL_II_NAMESPACE_CLOSE

#endif


