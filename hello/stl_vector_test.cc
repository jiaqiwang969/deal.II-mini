// ---------------------------------------------------------------------
//
// Copyright (C) 2012 - 2021 by the deal.II authors
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


// test for AlignedVector<unsigned int> which tests the basic stuff in the
// aligned vector

#include <deal.II/mystl/vector.h>

#include "../tests.h"




void 
test()
{
  int a[] = { 1,2,3,4,5 };
  mystl::vector<int> v1;
  mystl::vector<int> v2(10);
  mystl::vector<int> v3(10, 1);
  mystl::vector<int> v4(a, a + 5);
  mystl::vector<int> v5(v2);
  mystl::vector<int> v6(std::move(v2));
  mystl::vector<int> v7{ 1,2,3,4,5,6,7,8,9 };
  mystl::vector<int> v8, v9, v10;
  v8 = v3;
  v9 = std::move(v3);
  v10 = { 1,2,3,4,5,6,7,8,9 };

  deallog << "Constructor: ";
  for (unsigned int i = 0; i < v1.size(); ++i)
    deallog << v1[i] << " ";
  deallog << std::endl;

  v1.assign(8, 8);
  deallog << "Insertion: ";
  for (unsigned int i = 0; i < v1.size(); ++i)
    deallog << v1[i] << " ";
  deallog << std::endl;




//   FUN_AFTER(v1, v1.assign(8, 8));
//   FUN_AFTER(v1, v1.assign(a, a + 5));
//   FUN_AFTER(v1, v1.emplace(v1.begin(), 0));
//   FUN_AFTER(v1, v1.emplace_back(6));
//   FUN_AFTER(v1, v1.push_back(6));
//   FUN_AFTER(v1, v1.insert(v1.end(), 7));
//   FUN_AFTER(v1, v1.insert(v1.begin() + 3, 2, 3));
//   FUN_AFTER(v1, v1.insert(v1.begin(), a, a + 5));
//   FUN_AFTER(v1, v1.pop_back());
//   FUN_AFTER(v1, v1.erase(v1.begin()));
//   FUN_AFTER(v1, v1.erase(v1.begin(), v1.begin() + 2));
//   FUN_AFTER(v1, v1.reverse());
//   FUN_AFTER(v1, v1.swap(v4));
//   FUN_VALUE(*v1.begin());
//   FUN_VALUE(*(v1.end() - 1));
//   FUN_VALUE(*v1.rbegin());
//   FUN_VALUE(*(v1.rend() - 1));
//   FUN_VALUE(v1.front());
//   FUN_VALUE(v1.back());
//   FUN_VALUE(v1[0]);
//   FUN_VALUE(v1.at(1));
//   int* p = v1.data();
//   *p = 10;
//   *++p = 20;
//   p[1] = 30;
//   std::cout << " After change v1.data() :" << "\n";
//   COUT(v1);
//   std::cout << std::boolalpha;
//   FUN_VALUE(v1.empty());
//   std::cout << std::noboolalpha;
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.max_size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.resize(10));
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.shrink_to_fit());
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.resize(6, 6));
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.shrink_to_fit());
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.clear());
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.reserve(5));
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.reserve(20));
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   FUN_AFTER(v1, v1.shrink_to_fit());
//   FUN_VALUE(v1.size());
//   FUN_VALUE(v1.capacity());
//   PASSED;
// #if PERFORMANCE_TEST_ON
//   std::cout << "[--------------------- Performance Testing ---------------------]\n";
//   std::cout << "|---------------------|-------------|-------------|-------------|\n";
//   std::cout << "|      push_back      |";
// #if LARGER_TEST_DATA_ON
//   CON_TEST_P1(vector<int>, push_back, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
// #else
//   CON_TEST_P1(vector<int>, push_back, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
// #endif
//   std::cout << "\n";
//   std::cout << "|---------------------|-------------|-------------|-------------|\n";
//   PASSED;
// #endif
//   std::cout << "[----------------- End container test : vector -----------------]\n";
}
















// void
// test()
// {
//   using VEC = AlignedVector<unsigned int>;
//   VEC a(4);
//   deallog << "Constructor: ";
//   for (unsigned int i = 0; i < a.size(); ++i)
//     deallog << a[i] << " ";
//   deallog << std::endl;

//   a[2] = 1;
//   a.push_back(5);
//   a.push_back(42);

//   VEC b(a);
//   b.push_back(27);
//   a.insert_back(b.begin(), b.end());

//   deallog << "Insertion: ";
//   for (unsigned int i = 0; i < a.size(); ++i)
//     deallog << a[i] << " ";
//   deallog << std::endl;

//   a.resize(4);
//   deallog << "Shrinking: ";
//   for (unsigned int i = 0; i < a.size(); ++i)
//     deallog << a[i] << " ";
//   deallog << std::endl;

//   a.reserve(100);
//   deallog << "Reserve: ";
//   for (unsigned int i = 0; i < a.size(); ++i)
//     deallog << a[i] << " ";
//   deallog << std::endl;

//   a = b;
//   deallog << "Assignment: ";
//   for (unsigned int i = 0; i < a.size(); ++i)
//     deallog << a[i] << " ";
//   deallog << std::endl;

//   // check setting elements for large vectors
//   a.resize(0);
//   a.resize(100000, 1);
//   deallog << "Check large initialization: ";
//   for (unsigned int i = 0; i < 100000; ++i)
//     AssertDimension(a[i], 1);
//   deallog << "OK" << std::endl;

//   // check resize for large vectors
//   deallog << "Check large resize: ";
//   a.resize(200000, 2);
//   a.resize(400000);
//   for (unsigned int i = 0; i < 100000; ++i)
//     AssertDimension(a[i], 1);
//   for (unsigned int i = 100000; i < 200000; ++i)
//     AssertDimension(a[i], 2);
//   for (unsigned int i = 200000; i < 400000; ++i)
//     AssertDimension(a[i], 0);
//   deallog << "OK" << std::endl;
// }



int
main()
{
  initlog();

  test();
}