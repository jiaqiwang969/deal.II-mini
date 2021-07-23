#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

/**
 * Memento
 *
 * - 模式动机:
 * 
 *   让你保存和恢复一个对象以前的状态，而不透露其实现的细节。
 */

/**
 * - 模式定义:
 * 
 *   Memento接口提供了一种检索Memento元数据的方法，例如创建日期或名称。然而，它并没有暴露发起者的状态。
 */
namespace DesignPattern
{

    namespace Memento
    {

        class Memento
        {
        public:
            virtual std::string GetName() const = 0;
            virtual std::string date() const = 0;
            virtual std::string state() const = 0;
        };

        /**
         * Concrete Memento包含用于存储发起人状态的基础设施。
         */
        class ConcreteMemento : public Memento
        {
        private:
            std::string state_;
            std::string date_;

        public:
            ConcreteMemento(std::string state) : state_(state)
            {
                this->state_ = state;
                std::time_t now = std::time(0);
                this->date_ = std::ctime(&now);
            }
            /**
             * 发起人在恢复其状态时使用这种方法。
             */
            std::string state() const override
            {
                return this->state_;
            }
            /**
             * 其余的方法是由看守所用来显示元数据。
             */
            std::string GetName() const override
            {
                return this->date_ + " / (" + this->state_.substr(0, 9) + "...)";
            }
            std::string date() const override
            {
                return this->date_;
            }
        };

        /**
         * 发起人持有一些可能随时间变化的重要状态。它还定义了一种将状态保存在备忘录内的方法和另一种将状态恢复的方法。
         */
        class Originator
        {
            /**
             * @var string
             * 为了简单起见，发起人的状态被储存在一个单一的变量内。
             */
        private:
            std::string state_;

            std::string GenerateRandomString(int length = 10)
            {
                const char alphanum[] =
                    "0123456789"
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                    "abcdefghijklmnopqrstuvwxyz";
                int stringLength = sizeof(alphanum) - 1;

                std::string random_string;
                for (int i = 0; i < length; i++)
                {
                    random_string += alphanum[std::rand() % stringLength];
                }
                return random_string;
            }

        public:
            Originator(std::string state) : state_(state)
            {
                std::cout << "Originator: My initial state is: " << this->state_ << "\n";
            }
            /**
             * 发起人的业务逻辑可能会影响其内部状态。因此，在通过save()方法启动业务逻辑的方法之前，客户应该备份状态。
             */
            void DoSomething()
            {
                std::cout << "Originator: I'm doing something important.\n";
                this->state_ = this->GenerateRandomString(30);
                std::cout << "Originator: and my state has changed to: " << this->state_ << "\n";
            }

            /**
             * 将当前状态保存在一个Memento内。
             */
            Memento *Save()
            {
                return new ConcreteMemento(this->state_);
            }
            /**
             * 从一个备忘录对象中恢复发起人的状态。
             */
            void Restore(Memento *memento)
            {
                this->state_ = memento->state();
                std::cout << "Originator: My state has changed to: " << this->state_ << "\n";
            }
        };

        /**
         * Caretaker并不依赖于具体的Memento类。因此，它不能访问存储在Memento中的发起人的状态。它通过基础的Memento接口与所有的纪念品一起工作。
         */
        class Caretaker
        {
            /**
             * @var Memento[]
             */
        private:
            std::vector<Memento *> mementos_;

            /**
             * @var Originator
             */
            Originator *originator_;

        public:
            Caretaker(Originator *originator) : originator_(originator)
            {
                this->originator_ = originator;
            }

            void Backup()
            {
                std::cout << "\nCaretaker: Saving Originator's state...\n";
                this->mementos_.push_back(this->originator_->Save());
            }
            void Undo()
            {
                if (!this->mementos_.size())
                {
                    return;
                }
                Memento *memento = this->mementos_.back();
                this->mementos_.pop_back();
                std::cout << "Caretaker: Restoring state to: " << memento->GetName() << "\n";
                try
                {
                    this->originator_->Restore(memento);
                }
                catch (...)
                {
                    this->Undo();
                }
            }
            void ShowHistory() const
            {
                std::cout << "Caretaker: Here's the list of mementos:\n";
                for (Memento *memento : this->mementos_)
                {
                    std::cout << memento->GetName() << "\n";
                }
            }
        };
        /**
         * Client code.
         */

        void ClientCode()
        {
            Originator *originator = new Originator("Super-duper-super-puper-super.");
            Caretaker *caretaker = new Caretaker(originator);
            caretaker->Backup();
            originator->DoSomething();
            caretaker->Backup();
            originator->DoSomething();
            caretaker->Backup();
            originator->DoSomething();
            std::cout << "\n";
            caretaker->ShowHistory();
            std::cout << "\nClient: Now, let's rollback!\n\n";
            caretaker->Undo();
            std::cout << "\nClient: Once more!\n\n";
            caretaker->Undo();

            delete originator;
            delete caretaker;
        }

    }
}
// int main() {
//   std::srand(static_cast<unsigned int>(std::time(NULL)));
//   ClientCode();
//   return 0;
// }
