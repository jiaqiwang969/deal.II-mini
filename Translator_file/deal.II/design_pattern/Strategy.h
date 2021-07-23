#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

/**
 * Strategy Design Pattern
 *
 * 意图: 让你定义一个算法系列，把它们每个都放到一个单独的类中，并使它们的对象可以互换。
 */

/**
 * 策略接口声明了一些算法的所有支持版本所共有的操作。
 * Context使用这个接口来调用Concrete Strategies定义的算法。
 */
class Strategy
{
public:
    virtual ~Strategy() {}
    virtual std::string DoAlgorithm(const std::vector<std::string> &data) const = 0;
};

/**
 * Context定义了客户感兴趣的接口。
 */

class Context
{
    /**
     * @var Strategy Context维护对策略对象之一的引用。Context不知道一个策略的具体类别。它应该通过策略接口与所有策略一起工作。
     */
private:
    Strategy *strategy_;
    /**
     * 通常，Context通过构造函数接受一个策略，但也提供一个setter来在运行时改变它。
     */
public:
    Context(Strategy *strategy = nullptr) : strategy_(strategy)
    {
    }
    ~Context()
    {
        delete this->strategy_;
    }
    /**
     * 通常，Context允许在运行时替换一个策略对象。
     */
    void set_strategy(Strategy *strategy)
    {
        delete this->strategy_;
        this->strategy_ = strategy;
    }
    /**
     *    Context将一些工作委托给Strategy对象，而不是自己实现+多个版本的算法。
     */
    void DoSomeBusinessLogic() const
    {
        // ...
        std::cout << "Context: Sorting data using the strategy (not sure how it'll do it)\n";
        std::string result = this->strategy_->DoAlgorithm(std::vector<std::string>{"a", "e", "c", "b", "d"});
        std::cout << result << "\n";
        // ...
    }
};

/**
 * 具体的策略在遵循基本策略接口的同时实现了该算法。该接口使它们在上下文中可以互换。
 */
class ConcreteStrategyA : public Strategy
{
public:
    std::string DoAlgorithm(const std::vector<std::string> &data) const override
    {
        std::string result;
        std::for_each(std::begin(data), std::end(data), [&result](const std::string &letter)
                      { result += letter; });
        std::sort(std::begin(result), std::end(result));

        return result;
    }
};
class ConcreteStrategyB : public Strategy
{
    std::string DoAlgorithm(const std::vector<std::string> &data) const override
    {
        std::string result;
        std::for_each(std::begin(data), std::end(data), [&result](const std::string &letter)
                      { result += letter; });
        std::sort(std::begin(result), std::end(result));
        for (int i = 0; i < result.size() / 2; i++)
        {
            std::swap(result[i], result[result.size() - i - 1]);
        }

        return result;
    }
};
/**
 * 客户端代码挑选一个具体的策略并将其传递给上下文。客户端应该意识到策略之间的差异，以便做出正确的选择。
 */

void ClientCode()
{
    Context *context = new Context(new ConcreteStrategyA);
    std::cout << "Client: Strategy is set to normal sorting.\n";
    context->DoSomeBusinessLogic();
    std::cout << "\n";
    std::cout << "Client: Strategy is set to reverse sorting.\n";
    context->set_strategy(new ConcreteStrategyB);
    context->DoSomeBusinessLogic();
    delete context;
}

// int main()
// {
//     ClientCode();
//     return 0;
// }