#include <iostream>
#include <string>
/**
 * Mediator
 *
 * - 模式动机:
 * 
 *   让你减少对象之间混乱的依赖关系。该模式限制了对象之间的直接通信，迫使它们只通过一个中介对象进行合作。
 * 
 *   在用户与用户直接聊天的设计方案中，用户对象之间存在很强的关联性，将导致系统出现如下问题：
 * 
 *   系统结构复杂：对象之间存在大量的相互关联和调用，若有一个对象发生变化，则需要跟踪和该对象关联的其他所有对象，并进行适当处理。
 * 
 *   对象可重用性差：由于一个对象和其他对象具有很强的关联，若没有其他对象的支持，一个对象很难被另一个系统或模块重用，这些对象表现出来更像一个不可分割的整体，职责较为混乱。
 * 
 *   系统扩展性低：增加一个新的对象需要在原有相关对象上增加引用，增加新的引用关系也需要调整原有对象，系统耦合度很高，对象操作很不灵活，扩展性差。
 * 
 *   在面向对象的软件设计与开发过程中，根据“单一职责原则”，我们应该尽量将对象细化，使其只负责或呈现单一的职责。
 * 对于一个模块，可能由很多对象构成，而且这些对象之间可能存在相互的引用，为了减少对象两两之间复杂的引用关系，使之成为一个松耦合的系统，我们需要使用中介者模式，这就是中介者模式的模式动机。
 */

/**
 * - 模式定义: 
 * 
 *   中介者模式(Mediator Pattern)定义：用一个中介对象来封装一系列的对象交互，中介者使各对象不需要显式地相互引用，从而使其耦合松散，而且可以独立地改变它们之间的交互。中介者模式又称为调停者模式，它是一种对象行为型模式。
 *
 *   调解器接口声明了一个被组件用来通知调解器各种事件的方法。调解器可以对这些事件做出反应，并将执行结果传递给其他组件。
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
