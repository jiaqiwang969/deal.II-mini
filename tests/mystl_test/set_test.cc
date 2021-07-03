
// set test : 测试 set, multiset 的接口与它们 insert 的性能

#include <set>

#include <deal.II/mystl/set.h>
#include "../test.h"



void set_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[------------------ Run container test : set -------------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  int a[] = { 5,4,3,2,1 };
  mystl::set<int> s1;
  mystl::set<int, mystl::greater<int>> s2;
  mystl::set<int> s3(a, a + 5);
  mystl::set<int> s4(a, a + 5);
  mystl::set<int> s5(s3);
  mystl::set<int> s6(std::move(s3));
  mystl::set<int> s7;
  s7 = s4;
  mystl::set<int> s8;
  s8 = std::move(s4);
  mystl::set<int> s9{ 1,2,3,4,5 };
  mystl::set<int> s10;
  s10 = { 1,2,3,4,5 };

  for (int i = 5; i > 0; --i)
  {
    FUN_AFTER(s1, s1.emplace(i));
  }
  FUN_AFTER(s1, s1.emplace_hint(s1.begin(), 0));
  FUN_AFTER(s1, s1.erase(s1.begin()));
  FUN_AFTER(s1, s1.erase(0));
  FUN_AFTER(s1, s1.erase(1));
  FUN_AFTER(s1, s1.erase(s1.begin(), s1.end()));
  for (int i = 0; i < 5; ++i)
  {
    FUN_AFTER(s1, s1.insert(i));
  }
  FUN_AFTER(s1, s1.insert(a, a + 5));
  FUN_AFTER(s1, s1.insert(5));
  FUN_AFTER(s1, s1.insert(s1.end(), 5));
  FUN_VALUE(s1.count(5));
  FUN_VALUE(*s1.find(3));
  FUN_VALUE(*s1.lower_bound(3));
  FUN_VALUE(*s1.upper_bound(3));
  auto first = *s1.equal_range(3).first;
  auto second = *s1.equal_range(3).second;
  deallog << " s1.equal_range(3) : from " << first << " to " << second << std::endl;
  FUN_AFTER(s1, s1.erase(s1.begin()));
  FUN_AFTER(s1, s1.erase(1));
  FUN_AFTER(s1, s1.erase(s1.begin(), s1.find(3)));
  FUN_AFTER(s1, s1.clear());
  FUN_AFTER(s1, s1.swap(s5));
  FUN_VALUE(*s1.begin());
  FUN_VALUE(*s1.rbegin());
  deallog << std::boolalpha;
  FUN_VALUE(s1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(s1.size());
  FUN_VALUE(s1.max_size());
  deallog << "OK" << std::endl;
#if PERFORMANCE_TEST_ON
  std::cout << "[--------------------- Performance Testing ---------------------]" << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|       emplace       |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(set<int>, emplace, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#else
  CON_TEST_P1(set<int>, emplace, rand(), LEN1 _M, LEN2 _M, LEN3 _M);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  std::cout << "[------------------ End container test : set -------------------]" << std::endl;
}

void multiset_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[---------------- Run container test : multiset ----------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  int a[] = { 5,4,3,2,1 };
  mystl::multiset<int> s1;
  mystl::multiset<int, mystl::greater<int>> s2;
  mystl::multiset<int> s3(a, a + 5);
  mystl::multiset<int> s4(a, a + 5);
  mystl::multiset<int> s5(s3);
  mystl::multiset<int> s6(std::move(s3));
  mystl::multiset<int> s7;
  s7 = s4;
  mystl::multiset<int> s8;
  s8 = std::move(s4);
  mystl::multiset<int> s9{ 1,2,3,4,5 };
  mystl::multiset<int> s10;
  s10 = { 1,2,3,4,5 };

  for (int i = 5; i > 0; --i)
  {
    FUN_AFTER(s1, s1.emplace(i));
  }
  FUN_AFTER(s1, s1.emplace_hint(s1.begin(), 0));
  FUN_AFTER(s1, s1.erase(s1.begin()));
  FUN_AFTER(s1, s1.erase(0));
  FUN_AFTER(s1, s1.erase(1));
  FUN_AFTER(s1, s1.erase(s1.begin(), s1.end()));
  for (int i = 0; i < 5; ++i)
  {
    FUN_AFTER(s1, s1.insert(i));
  }
  FUN_AFTER(s1, s1.insert(a, a + 5));
  FUN_AFTER(s1, s1.insert(5));
  FUN_AFTER(s1, s1.insert(s1.end(), 5));
  FUN_VALUE(s1.count(5));
  FUN_VALUE(*s1.find(3));
  FUN_VALUE(*s1.lower_bound(3));
  FUN_VALUE(*s1.upper_bound(3));
  auto first = *s1.equal_range(3).first;
  auto second = *s1.equal_range(3).second;
  deallog << " s1.equal_range(3) : from " << first << " to " << second << std::endl;
  FUN_AFTER(s1, s1.erase(s1.begin()));
  FUN_AFTER(s1, s1.erase(1));
  FUN_AFTER(s1, s1.erase(s1.begin(), s1.find(3)));
  FUN_AFTER(s1, s1.clear());
  FUN_AFTER(s1, s1.swap(s5));
  FUN_VALUE(*s1.begin());
  FUN_VALUE(*s1.rbegin());
  deallog << std::boolalpha;
  FUN_VALUE(s1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(s1.size());
  FUN_VALUE(s1.max_size());
  PASSED;
#if PERFORMANCE_TEST_ON
  std::cout << "[--------------------- Performance Testing ---------------------]" << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|       emplace       |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(multiset<int>, emplace, rand(), LEN1 _M, LEN2 _M, LEN3 _M);
#else
  CON_TEST_P1(multiset<int>, emplace, rand(), LEN1 _S, LEN2 _S, LEN3 _S);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  std::cout << "[---------------- End container test : multiset ----------------]" << std::endl;
}




int
main()
{
  initlog();

  set_test();
  multiset_test();
}