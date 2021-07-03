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
#include "../testfun.h"

void test()
{
        deallog << "[===============================================================]" << std::endl;
        deallog << "[----------------- Run container test : vector -----------------]" << std::endl;
        deallog << "[-------------------------- API test ---------------------------]" << std::endl;

        int a[] = {1, 2, 3, 4, 5};
        mystl::vector<int> v1;
        mystl::vector<int> v2(10);
        mystl::vector<int> v3(10, 1);
        mystl::vector<int> v4(a, a + 5);
        mystl::vector<int> v5(v2);
        mystl::vector<int> v6(std::move(v2));
        mystl::vector<int> v7{1, 2, 3, 4, 5, 6, 7, 8, 9};
        mystl::vector<int> v8, v9, v10;
        v8 = v3;
        v9 = std::move(v3);
        v10 = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        FUN_AFTER(v1, v1.assign(a, a + 5));
        FUN_AFTER(v1, v1.emplace(v1.begin(), 0));
        FUN_AFTER(v1, v1.emplace_back(6));
        FUN_AFTER(v1, v1.push_back(6));
        FUN_AFTER(v1, v1.insert(v1.end(), 7));
        FUN_AFTER(v1, v1.insert(v1.begin() + 3, 2, 3));
        FUN_AFTER(v1, v1.insert(v1.begin(), a, a + 5));
        FUN_AFTER(v1, v1.pop_back());
        FUN_AFTER(v1, v1.erase(v1.begin()));
        FUN_AFTER(v1, v1.erase(v1.begin(), v1.begin() + 2));
        FUN_AFTER(v1, v1.reverse());
        FUN_AFTER(v1, v1.swap(v4));
        FUN_VALUE(*v1.begin());
        FUN_VALUE(*(v1.end() - 1));
        FUN_VALUE(*v1.rbegin());
        FUN_VALUE(*(v1.rend() - 1));
        FUN_VALUE(v1.front());
        FUN_VALUE(v1.back());
        FUN_VALUE(v1[0]);
        FUN_VALUE(v1.at(1));
        int *p = v1.data();
        *p = 10;
        *++p = 20;
        p[1] = 30;
        deallog << " After change v1.data() :" << std::endl;
        COUT(v1);
        deallog << std::boolalpha;
        FUN_VALUE(v1.empty());
        deallog << std::noboolalpha;
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.max_size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.resize(10));
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.shrink_to_fit());
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.resize(6, 6));
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.shrink_to_fit());
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.clear());
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.reserve(5));
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.reserve(20));
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        FUN_AFTER(v1, v1.shrink_to_fit());
        FUN_VALUE(v1.size());
        FUN_VALUE(v1.capacity());
        deallog << "OK" << std::endl;

        // Performance Testing 暂时无法写入到output中，原因时间可变，影响ctest测试结果
        // 如果想看该测试结果到话，可以进入到 .debug 文件，运行可执行文件查看。
#if PERFORMANCE_TEST_ON
        std::cout << "[--------------------- Performance Testing ---------------------]\n";
        std::cout << "|---------------------|-------------|-------------|-------------|\n";
        std::cout << "|      push_back      |";
#if LARGER_TEST_DATA_ON
        CON_TEST_P1(vector<int>, push_back, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
#else
        CON_TEST_P1(vector<int>, push_back, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#endif
        std::cout << "\n";
        std::cout << "|---------------------|-------------|-------------|-------------|\n";
        PASSED;
#endif
        std::cout << "[----------------- End container test : vector -----------------]\n";
}

int main()
{
        initlog();

        test();
}
