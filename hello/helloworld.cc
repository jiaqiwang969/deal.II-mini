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

//#include "../tests.h"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <deal.II/base/logstream.h>
// all include files you need here


// Function to initialize deallog. Normally, it should be called at
// the beginning of main() like
//
// initlog();
//
// This will open the correct output file, divert log output there and
// switch off screen output. If screen output is desired, provide the
// optional first argument as 'true'.
std::string   deallogname;
std::ofstream deallogfile;

void
initlog(bool                          console = false,
        const std::ios_base::fmtflags flags   = std::ios::showpoint |
                                              std::ios::left)
{
  deallogname = "output";
  deallogfile.open(deallogname.c_str());
  deallog.attach(deallogfile, true, flags);
  deallog.depth_console(console ? 10 : 0);
}


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
