﻿
// 仅仅针对 sort, binary_search 做了性能测试

#include <algorithm>

#include <deal.II/mystl/algorithm.h>
#include "../testfun.h"

// 函数性能测试宏定义
#define FUN_TEST1(mode, fun, count)                                                     \
  do                                                                                    \
  {                                                                                     \
    std::string fun_name = #fun;                                                        \
    srand((int)time(0));                                                                \
    char buf[10];                                                                       \
    clock_t start, end;                                                                 \
    int *arr = new int[count];                                                          \
    for (size_t i = 0; i < count; ++i)                                                  \
      *(arr + i) = rand();                                                              \
    start = clock();                                                                    \
    mode::fun(arr, arr + count);                                                        \
    end = clock();                                                                      \
    int n = static_cast<int>(static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000); \
    std::snprintf(buf, sizeof(buf), "%d", n);                                           \
    std::string t = buf;                                                                \
    t += "ms   |";                                                                      \
    deallog << std::setw(WIDE) << t;                                                    \
    delete[] arr;                                                                       \
  } while (0)

#define FUN_TEST2(mode, fun, count)                                                     \
  do                                                                                    \
  {                                                                                     \
    std::string fun_name = #fun;                                                        \
    srand((int)time(0));                                                                \
    char buf[10];                                                                       \
    clock_t start, end;                                                                 \
    int *arr = new int[count];                                                          \
    for (size_t i = 0; i < count; ++i)                                                  \
      *(arr + i) = rand();                                                              \
    start = clock();                                                                    \
    for (size_t i = 0; i < count; ++i)                                                  \
      mode::fun(arr, arr + count, rand());                                              \
    end = clock();                                                                      \
    int n = static_cast<int>(static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000); \
    std::snprintf(buf, sizeof(buf), "%d", n);                                           \
    std::string t = buf;                                                                \
    t += "ms   |";                                                                      \
    deallog << std::setw(WIDE) << t;                                                    \
    delete[] arr;                                                                       \
  } while (0)

void binary_search_test()
{
  deallog << "[------------------- function : binary_search ------------------]" << std::endl;
  deallog << "| orders of magnitude |";
  TEST_LEN(LEN1, LEN2, LEN3, WIDE);
  deallog << "|         std         |";
  FUN_TEST2(std, binary_search, LEN1);
  FUN_TEST2(std, binary_search, LEN2);
  FUN_TEST2(std, binary_search, LEN3);
  deallog << std::endl
          << "|        mystl        |";
  FUN_TEST2(mystl, binary_search, LEN1);
  FUN_TEST2(mystl, binary_search, LEN2);
  FUN_TEST2(mystl, binary_search, LEN3);
  deallog << std::endl;
}

void sort_test()
{
  deallog << "[----------------------- function : sort -----------------------]" << std::endl;
  deallog << "| orders of magnitude |";
  TEST_LEN(LEN1, LEN2, LEN3, WIDE);
  deallog << "|         std         |";
  FUN_TEST1(std, sort, LEN1);
  FUN_TEST1(std, sort, LEN2);
  FUN_TEST1(std, sort, LEN3);
  deallog << std::endl
          << "|        mystl        |";
  FUN_TEST1(mystl, sort, LEN1);
  FUN_TEST1(mystl, sort, LEN2);
  FUN_TEST1(mystl, sort, LEN3);
  deallog << std::endl;
}

void algorithm_performance_test()
{

#if PERFORMANCE_TEST_ON
  std::cout << "[===============================================================]" << std::endl;
  std::cout << "[--------------- Run algorithm performance test ----------------]" << std::endl;
  sort_test();
  binary_search_test();
  std::cout << "[--------------- End algorithm performance test ----------------]" << std::endl;
  std::cout << "[===============================================================]" << std::endl;
#endif // PERFORMANCE_TEST_ON
}

int main()
{
  initlog();

  algorithm_performance_test();

}