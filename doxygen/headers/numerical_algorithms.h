//include/deal.II-translator/A-headers/numerical_algorithms_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2015 by the deal.II authors
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


/**
 *    @defgroup numerics Numerical algorithms
 * 这个模块将一系列不同的类组合在一起，这些类通常在库中所有的基本三角形、DoFHandler和有限元类之上实现某种数值算法。它们之间一般是没有联系的。
 * 一些类，如DerivativeApproximation、KellyErrorEstimator和SolutionTransfer，作用于已经得到的解，并计算前两种情况下的派生量，或者帮助将一组向量从一个网格转移到另一个网格。
 * 命名空间MatrixCreator、MatrixTools和VectorTools提供了各种各样的服务，如创建拉普拉斯矩阵、将一个函数投影或内插到目前的有限元空间上，等等。
 * 与DoFTools和FETools函数的不同之处在于，它们对向量（即给定三角上的有限元函数空间的成员）进行工作，或者帮助创建它。另一方面，DoFTools函数只作用于给定的DoFHandler对象而不参考数据向量，FETools对象一般与有限元类一起工作，但同样没有任何相关的数据向量。
 *
 */


