#include <iostream>
#include <string>

/**
 * Command
 *
 * 模式动机:
 * 将一个请求转化为一个独立的对象，其中包含关于请求的所有信息。这种转换可以让你用不同的请求对方法进行参数化，延迟或排队执行请求，并支持可撤销操作。
 *
 * 在软件设计中，我们经常需要向某些对象发送请求，但是并不知道请求的接收者是谁，也不知道被请求的操作是哪个，我们只需在程序运行时指定具体的请求接收者即可，此时，可以使用命令模式来进行设计，使得请求发送者与请求接收者消除彼此之间的耦合，让对象之间的调用关系更加灵活。
 * 命令模式可以对发送者和接收者完全解耦，发送者与接收者之间没有直接引用关系，发送请求的对象只需要知道如何发送请求，而不必知道如何完成请求。这就是命令模式的模式动机。
 */
/**
 * 模式定义: 命令模式(Command
 * Pattern)：将一个请求封装为一个对象，从而使我们可用不同的请求对客户进行参数化；对请求排队或者记录请求日志，以及支持可撤销的操作。命令模式是一种对象行为型模式，其别名为动作(Action)模式或事务(Transaction)模式。
 *
 * 命令接口声明了一个执行命令的方法。
 */
namespace DesignPattern
{

    namespace Command
    {
        class Command
        {
        public:
            virtual ~Command()
            {
            }
            virtual void Execute() const = 0;
        };
        /**
         * 有些命令可以自行实现简单的操作。
         */
        class SimpleCommand : public Command
        {
        private:
            std::string pay_load_;

        public:
            explicit SimpleCommand(std::string pay_load) : pay_load_(pay_load)
            {
            }
            void Execute() const override
            {
                std::cout << "SimpleCommand: See, I can do simple things like printing (" << this->pay_load_ << ")\n";
            }
        };

        /**
         * Receiver类包含一些重要的业务逻辑。它们知道如何执行与执行请求相关的各种操作。事实上，任何类都可以作为一个接收者。
         */
        class Receiver
        {
        public:
            void DoSomething(const std::string &a)
            {
                std::cout << "Receiver: Working on (" << a << ".)\n";
            }
            void DoSomethingElse(const std::string &b)
            {
                std::cout << "Receiver: Also working on (" << b << ".)\n";
            }
        };

        /**
         * 然而，一些命令可以将更复杂的操作委托给其他对象，称为
         * "接收者"。
         */
        class ComplexCommand : public Command
        {
            /**
             * @var Receiver
             */
        private:
            Receiver *receiver_;
            /**
             * 上下文数据，为启动接收方的方法所需。
             */
            std::string a_;
            std::string b_;
            /**
             * 复杂的命令可以通过构造函数接受一个或几个接收器对象以及任何上下文数据。
             */
        public:
            ComplexCommand(Receiver *receiver, std::string a, std::string b) : receiver_(receiver), a_(a), b_(b)
            {
            }
            /**
             * 命令可以委托给接收器的任何方法。
             */
            void Execute() const override
            {
                std::cout << "ComplexCommand: Complex stuff should be done by a receiver object.\n";
                this->receiver_->DoSomething(this->a_);
                this->receiver_->DoSomethingElse(this->b_);
            }
        };

        /**
         * Invoker与一个或几个命令相关联。它向命令发送一个请求。
         */
        class Invoker
        {
            /**
             * @var Command
             */
        private:
            Command *on_start_;
            /**
             * @var Command
             */
            Command *on_finish_;
            /**
             * Initialize commands.
             */
        public:
            ~Invoker()
            {
                delete on_start_;
                delete on_finish_;
            }

            void SetOnStart(Command *command)
            {
                this->on_start_ = command;
            }
            void SetOnFinish(Command *command)
            {
                this->on_finish_ = command;
            }
            /**
             * Invoker不依赖于具体的命令或接收器类。Invoker通过执行一个命令，间接地将一个请求传递给接收器。
             */
            void DoSomethingImportant()
            {
                std::cout << "Invoker: Does anybody want something done before I begin?\n";
                if (this->on_start_)
                {
                    this->on_start_->Execute();
                }
                std::cout << "Invoker: ...doing something really important...\n";
                std::cout << "Invoker: Does anybody want something done after I finish?\n";
                if (this->on_finish_)
                {
                    this->on_finish_->Execute();
                }
            }
        };

    }
}
/**
 * 客户端代码可以用任何命令对调用者进行参数化。
 */

// int main() {
//   Invoker *invoker = new Invoker;
//   invoker->SetOnStart(new SimpleCommand("Say Hi!"));
//   Receiver *receiver = new Receiver;
//   invoker->SetOnFinish(new ComplexCommand(receiver, "Send email", "Save report"));
//   invoker->DoSomethingImportant();

//   delete invoker;
//   delete receiver;

//   return 0;
// }
