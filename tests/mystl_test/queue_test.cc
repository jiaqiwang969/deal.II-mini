

// queue test : 测试 queue, priority_queue 的接口和它们 push 的性能

#include <queue>

#include <deal.II/mystl/queue.h>
#include "../testfun.h"



void queue_print(mystl::queue<int> q)
{
  while (!q.empty())
  {
    deallog << " " << q.front();
    q.pop();
  }
  deallog << std::endl;
}

void p_queue_print(mystl::priority_queue<int> p)
{
  while (!p.empty())
  {
    deallog << " " << p.top();
    p.pop();
  }
  deallog << std::endl;
}

//  queue 的遍历输出
#define QUEUE_COUT(q) do {                       \
    std::string q_name = #q;                     \
    deallog << " " << q_name << " :";          \
    queue_print(q);                              \
} while(0)

// priority_queue 的遍历输出
#define P_QUEUE_COUT(p) do {                     \
    std::string p_name = #p;                     \
    deallog << " " << p_name << " :";          \
    p_queue_print(p);                            \
} while(0)

#define QUEUE_FUN_AFTER(con, fun) do {           \
  std::string fun_name = #fun;                   \
  deallog << " After " << fun_name << " :" << std::endl;  \
  fun;                                           \
  QUEUE_COUT(con);                               \
} while(0)

#define P_QUEUE_FUN_AFTER(con, fun) do {         \
  std::string fun_name = #fun;                   \
  deallog << " After " << fun_name << " :" << std::endl;  \
  fun;                                           \
  P_QUEUE_COUT(con);                             \
} while(0)

void queue_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[----------------- Run container test : queue ------------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  int a[] = { 1,2,3,4,5 };
  mystl::deque<int> d1(5);
  mystl::queue<int> q1;
  mystl::queue<int> q2(5);
  mystl::queue<int> q3(5, 1);
  mystl::queue<int> q4(a, a + 5);
  mystl::queue<int> q5(d1);
  mystl::queue<int> q6(std::move(d1));
  mystl::queue<int> q7(q2);
  mystl::queue<int> q8(std::move(q2));
  mystl::queue<int> q9;
  q9 = q3;
  mystl::queue<int> q10;
  q10 = std::move(q3);
  mystl::queue<int> q11{ 1,2,3,4,5 };
  mystl::queue<int> q12;
  q12 = { 1,2,3,4,5 };

  QUEUE_FUN_AFTER(q1, q1.push(1));
  QUEUE_FUN_AFTER(q1, q1.push(2));
  QUEUE_FUN_AFTER(q1, q1.push(3));
  QUEUE_FUN_AFTER(q1, q1.pop());
  QUEUE_FUN_AFTER(q1, q1.emplace(4));
  QUEUE_FUN_AFTER(q1, q1.emplace(5));
  deallog << std::boolalpha;
  FUN_VALUE(q1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(q1.size());
  FUN_VALUE(q1.front());
  FUN_VALUE(q1.back());
  while (!q1.empty())
  {
    QUEUE_FUN_AFTER(q1, q1.pop());
  }
  QUEUE_FUN_AFTER(q1, q1.swap(q4));
  QUEUE_FUN_AFTER(q1, q1.clear());
  PASSED;
#if PERFORMANCE_TEST_ON
  deallog << "[--------------------- Performance Testing ---------------------]" << std::endl;
  deallog << "|---------------------|-------------|-------------|-------------|" << std::endl;
  deallog << "|         push        |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(queue<int>, push, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
#else
  CON_TEST_P1(queue<int>, push, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#endif
  deallog << std::endl;
  deallog << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  deallog << "[----------------- End container test : queue ------------------]" << std::endl;
}

void priority_test()
{
  deallog << "[===============================================================]" << std::endl;
  deallog << "[------------- Run container test : priority_queue -------------]" << std::endl;
  deallog << "[-------------------------- API test ---------------------------]" << std::endl;
  int a[] = { 1,2,3,4,5 };
  mystl::vector<int> v1(5);
  mystl::priority_queue<int> p1;
  mystl::priority_queue<int> p2(5);
  mystl::priority_queue<int> p3(5, 1);
  mystl::priority_queue<int> p4(a, a + 5);
  mystl::priority_queue<int> p5(v1);
  mystl::priority_queue<int> p6(std::move(v1));
  mystl::priority_queue<int> p7(p2);
  mystl::priority_queue<int> p8(std::move(p2));
  mystl::priority_queue<int> p9;
  p9 = p3;
  mystl::priority_queue<int> p10;
  p10 = std::move(p3);
  mystl::priority_queue<int> p11{ 1,2,3,4,5 };
  mystl::priority_queue<int> p12;
  p12 = { 1,2,3,4,5 };

  P_QUEUE_FUN_AFTER(p1, p1.push(1));
  P_QUEUE_FUN_AFTER(p1, p1.push(5));
  P_QUEUE_FUN_AFTER(p1, p1.push(3));
  P_QUEUE_FUN_AFTER(p1, p1.pop());
  P_QUEUE_FUN_AFTER(p1, p1.emplace(7));
  P_QUEUE_FUN_AFTER(p1, p1.emplace(2));
  P_QUEUE_FUN_AFTER(p1, p1.emplace(8));
  deallog << std::boolalpha;
  FUN_VALUE(p1.empty());
  deallog << std::noboolalpha;
  FUN_VALUE(p1.size());
  FUN_VALUE(p1.top());
  while (!p1.empty())
  {
    P_QUEUE_FUN_AFTER(p1, p1.pop());
  }
  P_QUEUE_FUN_AFTER(p1, p1.swap(p4));
  P_QUEUE_FUN_AFTER(p1, p1.clear());
  PASSED;
#if PERFORMANCE_TEST_ON
  std::cout << "[--------------------- Performance Testing ---------------------]" << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  std::cout << "|         push        |";
#if LARGER_TEST_DATA_ON
  CON_TEST_P1(priority_queue<int>, push, rand(), LEN1 _LL, LEN2 _LL, LEN3 _LL);
#else
  CON_TEST_P1(priority_queue<int>, push, rand(), LEN1 _L, LEN2 _L, LEN3 _L);
#endif
  std::cout << std::endl;
  std::cout << "|---------------------|-------------|-------------|-------------|" << std::endl;
  PASSED;
#endif
  std::cout << "[------------- End container test : priority_queue -------------]" << std::endl;
}




int
main()
{
  initlog();

  queue_test();
}