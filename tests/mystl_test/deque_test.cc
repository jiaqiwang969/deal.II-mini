﻿
// deque test : 测试 deque 的接口和 push_front/push_back 的性能

#include <deque>

#include <deal.II/mystl/deque.h>
#include "../testfun.h"



void deque_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[----------------- Run container test : deque ------------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  int a[] = { 1,2,3,4,5 };
  mystl::deque<int> d1;
  mystl::deque<int> d2(5);
  mystl::deque<int> d3(5, 1);
  mystl::deque<int> d4(a, a + 5);
  mystl::deque<int> d5(d2);
  mystl::deque<int> d6(std::move(d2));
  mystl::deque<int> d7;
  d7 = d3;
  mystl::deque<int> d8;
  d8 = std::move(d3);
  mystl::deque<int> d9{ 1,2,3,4,5,6,7,8,9 };
  mystl::deque<int> d10;
  d10 = { 1,2,3,4,5,6,7,8,9 };

  FUN_AFTER(d1, d1.assign(5, 1));
  FUN_AFTER(d1, d1.assign(8, 8));
  FUN_AFTER(d1, d1.assign(a, a + 5));
  FUN_AFTER(d1, d1.assign({ 1,2,3,4,5 }));
  FUN_AFTER(d1, d1.insert(d1.end(), 6));
  FUN_AFTER(d1, d1.insert(d1.end() - 1, 2, 7));
  FUN_AFTER(d1, d1.insert(d1.begin(), a, a + 5));
  FUN_AFTER(d1, d1.erase(d1.begin()));
  FUN_AFTER(d1, d1.erase(d1.begin(), d1.begin() + 4));
  FUN_AFTER(d1, d1.emplace_back(8));
  FUN_AFTER(d1, d1.emplace_front(8));
  FUN_AFTER(d1, d1.emplace(d1.begin() + 1, 9));
  FUN_AFTER(d1, d1.push_front(1));
  FUN_AFTER(d1, d1.push_back(2));
  FUN_AFTER(d1, d1.pop_back());
  FUN_AFTER(d1, d1.pop_front());
  FUN_AFTER(d1, d1.shrink_to_fit());
  FUN_AFTER(d1, d1.resize(5));
  FUN_AFTER(d1, d1.resize(8, 8));
  FUN_AFTER(d1, d1.clear());
  FUN_AFTER(d1, d1.shrink_to_fit());
  FUN_AFTER(d1, d1.swap(d4));
  FUN_VALUE(*(d1.begin()));
  FUN_VALUE(*(d1.end() - 1));
  FUN_VALUE(*(d1.rbegin()));
  FUN_VALUE(*(d1.rend() - 1));
  FUN_VALUE(d1.front());
  FUN_VALUE(d1.back());
  FUN_VALUE(d1.at(1));
  FUN_VALUE(d1[2]);
  deallog << std::boolalpha;
  FUN_VALUE(d1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(d1.size());
  FUN_VALUE(d1.max_size());
  deallog << "OK" << std::endl;
#if PERFORMANCE_TEST_ON
  std::cout << "[--------------------- Performance Testing ---------------------]" << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|     push_front      |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(deque<int>, push_front, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
#else
  CON_TEST_P1(deque<int>, push_front, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|     push_back       |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(deque<int>, push_back, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
#else
  CON_TEST_P1(deque<int>, push_back, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  std::cout << "[----------------- End container test : deque ------------------]" << std::endl;
}




int
main()
{
  initlog();

  deque_test();
}