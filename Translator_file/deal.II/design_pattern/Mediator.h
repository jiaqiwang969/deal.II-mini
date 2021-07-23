#include <iostream>
#include <string>
/**
 * Mediator Design Pattern
 *
 * 意图: 让你减少对象之间混乱的依赖关系。该模式限制了对象之间的直接通信，迫使它们只通过一个中介对象进行合作。
 */

/**
 * 调解器接口声明了一个被组件用来通知调解器各种事件的方法。调解器可以对这些事件做出反应，并将执行结果传递给其他组件。
 */
namespace DesignPattern
{

    namespace Mediator
    {

        class BaseComponent;
        class Mediator
        {
        public:
            virtual void Notify(BaseComponent *sender, std::string event) const = 0;
        };

        /**
 * 基础组件提供了在组件对象内存储调解人实例的基本功能。
 */
        class BaseComponent
        {
        protected:
            Mediator *mediator_;

        public:
            BaseComponent(Mediator *mediator = nullptr) : mediator_(mediator)
            {
            }
            void set_mediator(Mediator *mediator)
            {
                this->mediator_ = mediator;
            }
        };

        /**
 * 具体的组件实现各种功能。它们不依赖其他组件。它们也不依赖于任何具体的调解器类。
 */
        class Component1 : public BaseComponent
        {
        public:
            void DoA()
            {
                std::cout << "Component 1 does A.\n";
                this->mediator_->Notify(this, "A");
            }
            void DoB()
            {
                std::cout << "Component 1 does B.\n";
                this->mediator_->Notify(this, "B");
            }
        };

        class Component2 : public BaseComponent
        {
        public:
            void DoC()
            {
                std::cout << "Component 2 does C.\n";
                this->mediator_->Notify(this, "C");
            }
            void DoD()
            {
                std::cout << "Component 2 does D.\n";
                this->mediator_->Notify(this, "D");
            }
        };

        /**
 * 具体的调解人通过协调几个组件来实现合作行为。
 */
        class ConcreteMediator : public Mediator
        {
        private:
            Component1 *component1_;
            Component2 *component2_;

        public:
            ConcreteMediator(Component1 *c1, Component2 *c2) : component1_(c1), component2_(c2)
            {
                this->component1_->set_mediator(this);
                this->component2_->set_mediator(this);
            }
            void Notify(BaseComponent *sender, std::string event) const override
            {
                if (event == "A")
                {
                    std::cout << "Mediator reacts on A and triggers following operations:\n";
                    this->component2_->DoC();
                }
                if (event == "D")
                {
                    std::cout << "Mediator reacts on D and triggers following operations:\n";
                    this->component1_->DoB();
                    this->component2_->DoC();
                }
            }
        };

        /**
 * The client code.
 */

        void ClientCode()
        {
            Component1 *c1 = new Component1;
            Component2 *c2 = new Component2;
            ConcreteMediator *mediator = new ConcreteMediator(c1, c2);
            std::cout << "Client triggers operation A.\n";
            c1->DoA();
            std::cout << "\n";
            std::cout << "Client triggers operation D.\n";
            c2->DoD();

            delete c1;
            delete c2;
            delete mediator;
        }

    }
}

// int main() {
//   ClientCode();
//   return 0;
// }