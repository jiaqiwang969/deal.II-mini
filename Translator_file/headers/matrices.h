//include/deal.II-translator/A-headers/matrices_0.txt
// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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
 *    @defgroup Matrices Matrix classes
 * deal.II带有许多不同的矩阵类，是为使用矩阵的各种目的而定制的。例如，有完整的矩阵、使用不同存储方案的稀疏矩阵、由单个块组成的矩阵，以及作为其他线性代数类接口实现的矩阵。在可能的情况下，所有这些实现都共享一个共同的接口，该接口至少包含编写迭代线性求解器所需的操作（见 @ref
 * Solvers ），但也包含从矩阵中读取和写入的元素访问。
 * @ref Preconditioners
 * 也是矩阵类，因为它们对向量进行线性运算。
 *
 * @ingroup LAC
 *
 */


/**
 *    @defgroup Matrix1 Basic matrices
 * 这些是由deal.II提供的实际矩阵类。可以在其中存储数值并检索它们。此外，它们还提供了线性求解器所需的全部接口（见 @ref
 * Solvers  ）。
 * 在这组矩阵中，有完整的矩阵、不同的稀疏矩阵和块状矩阵。此外，其他线性代数库（例如PETScWrappers）的接口中的一些类是矩阵。
 * 大多数deal.II稀疏矩阵类与它们的稀疏度模式分开，以使存储具有相同稀疏度模式的几个矩阵更有效率。更多信息见 @ref
 * Sparsity 。
 *
 * @ingroup Matrices
 *
 */


/**
 *    @defgroup Matrix2 Derived matrices
 * 这些矩阵是建立在基本矩阵之上的。它们使用 @ref ConceptMatrixType "MatrixType概念 "
 * 所定义的接口进行特殊操作。
 *
 * @ingroup Matrices
 *
 */


