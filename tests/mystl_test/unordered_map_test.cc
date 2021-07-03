﻿
// unordered_map test : 测试 unordered_map, unordered_multimap 的接口与它们 insert 的性能

#include <unordered_map>

#include <deal.II/mystl/unordered_map.h>
//#include "map_test.h"
#include "../testfun.h"

void unordered_map_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[-------------- Run container test : unordered_map -------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  mystl::vector<PAIR> v;
  for (int i = 0; i < 5; ++i)
    v.push_back(PAIR(5 - i, 5 - i));
  mystl::unordered_map<int, int> um1;
  mystl::unordered_map<int, int> um2(520);
  mystl::unordered_map<int, int> um3(520, mystl::hash<int>());
  mystl::unordered_map<int, int> um4(520, mystl::hash<int>(), mystl::equal_to<int>());
  mystl::unordered_map<int, int> um5(v.begin(), v.end());
  mystl::unordered_map<int, int> um6(v.begin(), v.end(), 100);
  mystl::unordered_map<int, int> um7(v.begin(), v.end(), 100, mystl::hash<int>());
  mystl::unordered_map<int, int> um8(v.begin(), v.end(), 100, mystl::hash<int>(), mystl::equal_to<int>());
  mystl::unordered_map<int, int> um9(um5);
  mystl::unordered_map<int, int> um10(std::move(um5));
  mystl::unordered_map<int, int> um11;
  um11 = um6;
  mystl::unordered_map<int, int> um12;
  um12 = std::move(um6);
  mystl::unordered_map<int, int> um13{PAIR(1, 1), PAIR(2, 3), PAIR(3, 3)};
  mystl::unordered_map<int, int> um14;
  um14 = {PAIR(1, 1), PAIR(2, 3), PAIR(3, 3)};

  MAP_FUN_AFTER(um1, um1.emplace(1, 1));
  MAP_FUN_AFTER(um1, um1.emplace_hint(um1.begin(), 1, 2));
  MAP_FUN_AFTER(um1, um1.insert(PAIR(2, 2)));
  MAP_FUN_AFTER(um1, um1.insert(um1.end(), PAIR(3, 3)));
  MAP_FUN_AFTER(um1, um1.insert(v.begin(), v.end()));
  MAP_FUN_AFTER(um1, um1.erase(um1.begin()));
  MAP_FUN_AFTER(um1, um1.erase(um1.begin(), um1.find(3)));
  MAP_FUN_AFTER(um1, um1.erase(1));
  deallog << std::boolalpha;
  FUN_VALUE(um1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.max_bucket_count());
  FUN_VALUE(um1.bucket(1));
  FUN_VALUE(um1.bucket_size(um1.bucket(5)));
  MAP_FUN_AFTER(um1, um1.clear());
  MAP_FUN_AFTER(um1, um1.swap(um7));
  MAP_VALUE(*um1.begin());
  FUN_VALUE(um1.at(1));
  FUN_VALUE(um1[1]);
  deallog << std::boolalpha;
  FUN_VALUE(um1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.max_size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.max_bucket_count());
  FUN_VALUE(um1.bucket(1));
  FUN_VALUE(um1.bucket_size(um1.bucket(1)));
  MAP_FUN_AFTER(um1, um1.reserve(1000));
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.bucket_size(1));
  FUN_VALUE(um1.bucket_size(2));
  FUN_VALUE(um1.bucket_size(3));
  MAP_FUN_AFTER(um1, um1.rehash(150));
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.count(1));
  MAP_VALUE(*um1.find(3));
  auto first = *um1.equal_range(3).first;
  auto second = *um1.equal_range(3).second;
  deallog << " um1.equal_range(3) : from <" << first.first << ", " << first.second
          << "> to <" << second.first << ", " << second.second << ">" << std::endl;
  FUN_VALUE(um1.load_factor());
  FUN_VALUE(um1.max_load_factor());
  MAP_FUN_AFTER(um1, um1.max_load_factor(1.5f));
  FUN_VALUE(um1.max_load_factor());
  deallog << "OK" << std::endl;
#if PERFORMANCE_TEST_ON
  deallog << "[--------------------- Performance Testing ---------------------]" << std::endl;
  deallog << "|---------------------|-------------|-------------|-------------|" << std::endl;
  deallog << "|       emplace       |";
#if LARGER_TEST_DATA_ON
  MAP_EMPLACE_TEST(unordered_map, LEN1 _M, LEN2 _M, LEN3 _M);
#else
  MAP_EMPLACE_TEST(unordered_map, LEN1 _S, LEN2 _S, LEN3 _S);
#endif
  deallog << std::endl;
  deallog << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  deallog << "[-------------- End container test : unordered_map -------------]" << std::endl;
}

void unordered_multimap_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[----------- Run container test : unordered_multimap -----------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  mystl::vector<PAIR> v;
  for (int i = 0; i < 5; ++i)
    v.push_back(PAIR(5 - i, 5 - i));
  mystl::unordered_multimap<int, int> um1;
  mystl::unordered_multimap<int, int> um2(520);
  mystl::unordered_multimap<int, int> um3(520, mystl::hash<int>());
  mystl::unordered_multimap<int, int> um4(520, mystl::hash<int>(), mystl::equal_to<int>());
  mystl::unordered_multimap<int, int> um5(v.begin(), v.end());
  mystl::unordered_multimap<int, int> um6(v.begin(), v.end(), 100);
  mystl::unordered_multimap<int, int> um7(v.begin(), v.end(), 100, mystl::hash<int>());
  mystl::unordered_multimap<int, int> um8(v.begin(), v.end(), 100, mystl::hash<int>(), mystl::equal_to<int>());
  mystl::unordered_multimap<int, int> um9(um5);
  mystl::unordered_multimap<int, int> um10(std::move(um5));
  mystl::unordered_multimap<int, int> um11;
  um11 = um6;
  mystl::unordered_multimap<int, int> um12;
  um12 = std::move(um6);
  mystl::unordered_multimap<int, int> um13{PAIR(1, 1), PAIR(2, 3), PAIR(3, 3)};
  mystl::unordered_multimap<int, int> um14;
  um14 = {PAIR(1, 1), PAIR(2, 3), PAIR(3, 3)};

  MAP_FUN_AFTER(um1, um1.emplace(1, 1));
  MAP_FUN_AFTER(um1, um1.emplace_hint(um1.begin(), 1, 2));
  MAP_FUN_AFTER(um1, um1.insert(PAIR(2, 2)));
  MAP_FUN_AFTER(um1, um1.insert(um1.end(), PAIR(3, 3)));
  MAP_FUN_AFTER(um1, um1.insert(v.begin(), v.end()));
  MAP_FUN_AFTER(um1, um1.erase(um1.begin()));
  MAP_FUN_AFTER(um1, um1.erase(um1.begin(), um1.find(3)));
  MAP_FUN_AFTER(um1, um1.erase(1));
  deallog << std::boolalpha;
  FUN_VALUE(um1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.max_bucket_count());
  FUN_VALUE(um1.bucket(1));
  FUN_VALUE(um1.bucket_size(um1.bucket(5)));
  MAP_FUN_AFTER(um1, um1.clear());
  MAP_FUN_AFTER(um1, um1.swap(um7));
  MAP_VALUE(*um1.begin());
  deallog << std::boolalpha;
  FUN_VALUE(um1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.max_size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.max_bucket_count());
  FUN_VALUE(um1.bucket(1));
  FUN_VALUE(um1.bucket_size(um1.bucket(1)));
  MAP_FUN_AFTER(um1, um1.reserve(1000));
  FUN_VALUE(um1.size());
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.bucket_size(1));
  FUN_VALUE(um1.bucket_size(2));
  FUN_VALUE(um1.bucket_size(3));
  MAP_FUN_AFTER(um1, um1.rehash(150));
  FUN_VALUE(um1.bucket_count());
  FUN_VALUE(um1.count(1));
  MAP_VALUE(*um1.find(3));
  auto first = *um1.equal_range(3).first;
  auto second = *um1.equal_range(3).second;
  deallog << " um1.equal_range(3) : from <" << first.first << ", " << first.second
          << "> to <" << second.first << ", " << second.second << ">" << std::endl;
  FUN_VALUE(um1.load_factor());
  FUN_VALUE(um1.max_load_factor());
  MAP_FUN_AFTER(um1, um1.max_load_factor(1.5f));
  FUN_VALUE(um1.max_load_factor());
  deallog << "OK" << std::endl;;
#if PERFORMANCE_TEST_ON
  std::cout << "[--------------------- Performance Testing ---------------------]" << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|       emplace       |";
#if LARGER_TEST_DATA_ON
  MAP_EMPLACE_TEST(unordered_multimap, LEN1 _M, LEN2 _M, LEN3 _M);
#else
  MAP_EMPLACE_TEST(unordered_multimap, LEN1 _S, LEN2 _S, LEN3 _S);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  std::cout << "[----------- End container test : unordered_multimap -----------]" << std::endl;
}

int main()
{
  initlog();

  unordered_map_test();
  unordered_multimap_test();
}