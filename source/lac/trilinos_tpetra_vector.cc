// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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

#include <deal.II/lac/trilinos_tpetra_vector.templates.h>

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  ifdef DEAL_II_WITH_MPI

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  namespace TpetraWrappers
  {
    template class Vector<float>;
    template class Vector<double>;
#    ifdef DEAL_II_WITH_COMPLEX_VALUES
    template class Vector<std::complex<float>>;
    template class Vector<std::complex<double>>;
#    endif
  } // namespace TpetraWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#  endif
#endif
