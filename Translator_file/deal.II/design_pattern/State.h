#include <iostream>
#include <typeinfo>
/**
 * State Design Pattern
 *
 *  意图: 当一个对象的内部状态改变时，让它改变其行为。它看起来就像该对象改变了它的类别。
 */

/**
 * 基状态类声明了所有具体状态应该实现的方法，并且还提供了一个与状态相关的上下文对象的反向引用。这个反向引用可以被状态用来将上下文过渡到另一个状态。
 */

class Context;

class State
{
    /**
   * @var Context
   */
protected:
    Context *context_;

public:
    virtual ~State()
    {
    }

    void set_context(Context *context)
    {
        this->context_ = context;
    }

    virtual void Handle1() = 0;
    virtual void Handle2() = 0;
};

/**
 * Context定义了客户感兴趣的接口。它还维护一个对State子类实例的引用，该实例代表Context的当前状态。
 */
class Context
{
    /**
   * @var State 对Context当前状态的引用。
   */
private:
    State *state_;

public:
    Context(State *state) : state_(nullptr)
    {
        this->TransitionTo(state);
    }
    ~Context()
    {
        delete state_;
    }
    /**
   * Context允许在运行时改变State对象。
   */
    void TransitionTo(State *state)
    {
        std::cout << "Context: Transition to " << typeid(*state).name() << ".\n";
        if (this->state_ != nullptr)
            delete this->state_;
        this->state_ = state;
        this->state_->set_context(this);
    }
    /**
   * Context将其部分行为委托给当前的State对象。
   */
    void Request1()
    {
        this->state_->Handle1();
    }
    void Request2()
    {
        this->state_->Handle2();
    }
};

/**
 * 具体状态实现各种行为，与Context的状态相关。
 */

class ConcreteStateA : public State
{
public:
    void Handle1() override;

    void Handle2() override
    {
        std::cout << "ConcreteStateA handles request2.\n";
    }
};

class ConcreteStateB : public State
{
public:
    void Handle1() override
    {
        std::cout << "ConcreteStateB handles request1.\n";
    }
    void Handle2() override
    {
        std::cout << "ConcreteStateB handles request2.\n";
        std::cout << "ConcreteStateB wants to change the state of the context.\n";
        this->context_->TransitionTo(new ConcreteStateA);
    }
};

void ConcreteStateA::Handle1()
{
    {
        std::cout << "ConcreteStateA handles request1.\n";
        std::cout << "ConcreteStateA wants to change the state of the context.\n";

        this->context_->TransitionTo(new ConcreteStateB);
    }
}

/**
 * The client code.
 */
void ClientCode()
{
    Context *context = new Context(new ConcreteStateA);
    context->Request1();
    context->Request2();
    delete context;
}

// int main() {
//   ClientCode();
//   return 0;
// }