#include <iostream>
#include <string>
#include <vector>

namespace Pattern
{

    /**
     * ChainofResponsibility
     *
     * - 模式动机:
     *
     * 让你沿着一个处理程序链传递请求。当收到一个请求时，每个处理程序决定是处理该请求还是将其传递给链中的下一个处理程序。
     */
    /**
     * - 模式定义:
     *
     * 处理程序接口声明了一个用于建立处理程序链的方法。它还声明了一个执行请求的方法。
     */

    namespace ChainofResponsibility
    {

        class Handler
        {
        public:
            virtual Handler *SetNext(Handler *handler) = 0;
            virtual std::string Handle(std::string request) = 0;
        };
        /**
         * 默认的连锁行为可以在一个基础处理程序类里面实现。
         */
        class AbstractHandler : public Handler
        {
            /**
             * @var Handler
             */
        private:
            Handler *next_handler_;

        public:
            AbstractHandler() : next_handler_(nullptr)
            {
            }
            Handler *SetNext(Handler *handler) override
            {
                this->next_handler_ = handler;
                // 从这里返回一个处理程序将让我们以一种方便的方式连接处理程序，就像这样。
                // $monkey->setNext($squirrel)->setNext($dog);
                return handler;
            }
            std::string Handle(std::string request) override
            {
                if (this->next_handler_)
                {
                    return this->next_handler_->Handle(request);
                }

                return {};
            }
        };

        /**
         * 所有的具体处理程序要么处理一个请求，要么将其传递给链上的下一个处理程序。
         */
        class MonkeyHandler : public AbstractHandler
        {
        public:
            std::string Handle(std::string request) override
            {
                if (request == "Banana")
                {
                    return "Monkey: I'll eat the " + request + ".\n";
                }
                else
                {
                    return AbstractHandler::Handle(request);
                }
            }
        };
        class SquirrelHandler : public AbstractHandler
        {
        public:
            std::string Handle(std::string request) override
            {
                if (request == "Nut")
                {
                    return "Squirrel: I'll eat the " + request + ".\n";
                }
                else
                {
                    return AbstractHandler::Handle(request);
                }
            }
        };
        class DogHandler : public AbstractHandler
        {
        public:
            std::string Handle(std::string request) override
            {
                if (request == "MeatBall")
                {
                    return "Dog: I'll eat the " + request + ".\n";
                }
                else
                {
                    return AbstractHandler::Handle(request);
                }
            }
        };
        /**
         * 客户端代码通常适合于与单个处理程序一起工作。在大多数情况下，它甚至不知道处理程序是一个链的一部分。
         */
        void ClientCode(Handler &handler)
        {
            std::vector<std::string> food = {"Nut", "Banana", "Cup of coffee"};
            for (const std::string &f : food)
            {
                std::cout << "Client: Who wants a " << f << "?\n";
                const std::string result = handler.Handle(f);
                if (!result.empty())
                {
                    std::cout << "  " << result;
                }
                else
                {
                    std::cout << "  " << f << " was left untouched.\n";
                }
            }
        }

    }
}

/**
 * 客户端代码的另一部分构建了实际的链。
 */
// int main() {
//   MonkeyHandler *monkey = new MonkeyHandler;
//   SquirrelHandler *squirrel = new SquirrelHandler;
//   DogHandler *dog = new DogHandler;
//   monkey->SetNext(squirrel)->SetNext(dog);

//   /**
//    *  客户端应该能够向任何处理程序发送请求，而不仅仅是链中的第一个。
//    */
//   std::cout << "Chain: Monkey > Squirrel > Dog\n\n";
//   ClientCode(*monkey);
//   std::cout << "\n";
//   std::cout << "Subchain: Squirrel > Dog\n\n";
//   ClientCode(*squirrel);

//   delete monkey;
//   delete squirrel;
//   delete dog;

//   return 0;
// }
