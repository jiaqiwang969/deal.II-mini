//include/deal.II-translator/A-headers/solvers_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2014 by the deal.II authors
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
 *    @defgroup Solvers Linear solver classes
 * 为了正常工作，将矩阵和向量类作为模板参数的求解器要求这些类满足某种最小的接口，可以从求解器内部使用。对于迭代求解器，这个接口被定义在求解器类中。此外，求解器使用从SolverControl类派生出来的类的对象进行控制（例如其派生类ReductionControl），以确定最大的迭代次数或所需的公差。
 * 如果在配置过程中检测到（见ReadMe文件），一些稀疏的直接求解器也被支持。
 *
 * @ingroup LAC
 *
 */


