// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2021 by the deal.II Authors
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

// a short (a few lines) description of what the program does

#include "../tests.h"

// all include files you need here

int main ()
{
  // Initialize deallog for test output.
  // This also reroutes deallog output to a file "output".
  initlog();

  // your testcode here:
  int i = 0;
  deallog << i << std::endl;

  return 0;
}
