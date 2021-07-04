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
 * @defgroup Solvers Linear solver classes
 *
 * This module groups iterative and direct solvers, eigenvalue solvers, and
 * some control classes. All these classes operate on objects of the
 * @ref Matrices "matrix" and @ref Vectors "vector classes" defined in deal.II.
 *
 * In order to work properly, solvers that take matrix and vector classes as
 * template arguments require that these classes satisfy a certain minimal
 * interface that can be used from inside the solver. For iterative solvers,
 * this interface is defined in the Solver class. In addition, solvers are
 * controlled using objects of classes that are derived from the SolverControl
 * class (for example its derived class ReductionControl), in order to
 * determine the maximal number of iterations or a desired tolerance.
 *
 * If detected during configuration (see the ReadMe file), some sparse direct
 * solvers are also supported.
 *
 * @ingroup LAC
 */
