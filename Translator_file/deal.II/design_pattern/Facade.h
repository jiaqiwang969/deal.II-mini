#include <iostream>
#include <string>

/**
 * Facade
 *
 * 模式动机:
 * 为一个库、一个框架或任何其他复杂的类集提供一个简化的接口。
 */

/**
 * 模式定义: 外观模式(Facade
 * Pattern)：外部与一个子系统的通信必须通过一个统一的外观对象进行，为子系统中的一组接口提供一个一致的界面，外观模式定义了一个高层接口，这个接口使得这一子系统更加容易使用。外观模式又称为门面模式，它是一种对象结构型模式。
 * 子系统可以直接接受来自门面或客户端的请求。在任何情况下，对于子系统来说，门面是另一个客户端，它不是子系统的一部分。
 */
namespace DesignPattern
{

    namespace Facade
    {

        class Subsystem1
        {
        public:
            std::string Operation1() const
            {
                return "Subsystem1: Ready!\n";
            }
            // ...
            std::string OperationN() const
            {
                return "Subsystem1: Go!\n";
            }
        };
        /**
         * 有些门面可以同时与多个子系统一起工作。
         */
        class Subsystem2
        {
        public:
            std::string Operation1() const
            {
                return "Subsystem2: Get ready!\n";
            }
            // ...
            std::string OperationZ() const
            {
                return "Subsystem2: Fire!\n";
            }
        };

        /**
         * Facade类为一个或几个子系统的复杂逻辑提供一个简单的接口。Facade将客户端的请求委托给子系统中适当的对象。Facade还负责管理它们的生命周期。所有这些都使客户端免受子系统不必要的复杂性的影响。
         */
        class Facade
        {
        protected:
            Subsystem1 *subsystem1_;
            Subsystem2 *subsystem2_;
            /**
             * 根据你的应用程序的需要，你可以向Facade提供现有的子系统对象，或者强迫Facade自己创建它们。
             */
        public:
            /**
             * 在这种情况下，我们将把内存所有权委托给Facade类。
             */
            Facade(
                Subsystem1 *subsystem1 = nullptr,
                Subsystem2 *subsystem2 = nullptr)
            {
                this->subsystem1_ = subsystem1 ?: new Subsystem1;
                this->subsystem2_ = subsystem2 ?: new Subsystem2;
            }
            ~Facade()
            {
                delete subsystem1_;
                delete subsystem2_;
            }
            /**
             * Facade的方法是通往子系统复杂功能的方便捷径。然而，客户只能获得子系统的一小部分功能。
             */
            std::string Operation()
            {
                std::string result = "Facade initializes subsystems:\n";
                result += this->subsystem1_->Operation1();
                result += this->subsystem2_->Operation1();
                result += "Facade orders subsystems to perform the action:\n";
                result += this->subsystem1_->OperationN();
                result += this->subsystem2_->OperationZ();
                return result;
            }
        };

        /**
         * 客户端代码通过Facade提供的简单接口与复杂的子系统一起工作。当facade管理子系统的生命周期时，客户端甚至可能不知道子系统的存在。这种方法可以让你把复杂度控制住。
         */
        void ClientCode(Facade *facade)
        {
            // ...
            std::cout << facade->Operation();
            // ...
        }
        /**
         * 客户端代码可能已经创建了一些子系统的对象。在这种情况下，用这些对象初始化Facade而不是让Facade创建新的实例可能是值得的。
         */

    }
}
