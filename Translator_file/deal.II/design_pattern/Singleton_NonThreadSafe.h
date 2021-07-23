#include <iostream>
#include <string>
#include <thread>

namespace DesignPattern
{

    namespace CreationalDesign
    {

        /**
         * Singleton
         *
         * - 模式动机:
         * 让你确保一个类只有一个实例，同时为这个实例提供一个全局访问点。
         *
         * 对于系统中的某些类来说，只有一个实例很重要，例如，一个系统中可以存在多个打印任务，但是只能有一个正在工作的任务；一个系统只能有一个窗口管理器或文件系统；一个系统只能有一个计时工具或ID（序号）生成器。
         *
         * 如何保证一个类只有一个实例并且这个实例易于被访问呢？定义一个全局变量可以确保对象随时都可以被访问，但不能防止我们实例化多个对象。
         *
         * 一个更好的解决办法是让类自身负责保存它的唯一实例。这个类可以保证没有其他实例被创建，并且它可以提供一个访问该实例的方法。这就是单例模式的模式动机。
         */
        /**
         *
         * - 模式定义:
         *
         * 单例模式(Singleton
         * Pattern)：单例模式确保某一个类只有一个实例，而且自行实例化并向整个系统提供这个实例，这个类称为单例类，它提供全局访问的方法。
         *
         * 单例模式的要点有三个：一是某个类只能有一个实例；二是它必须自行创建这个实例；三是它必须自行向整个系统提供这个实例。单例模式是一种对象创建型模式。单例模式又名单件模式或单态模式。
         *
         * Singleton类定义了
         * "GetInstance"方法，作为构造函数的替代方法，让客户反复访问该类的同一实例。
         */

        namespace Singleton
        {

            class Singleton
            {

                /**
                 * Singleton的构造函数应该始终是私有的，以防止用`new`操作符直接调用构造。
                 */

            protected:
                Singleton(const std::string value) : value_(value)
                {
                }

                static Singleton *singleton_;

                std::string value_;

            public:
                /**
                 * Singletons 不应该是可克隆的。
                 */
                Singleton(Singleton &other) = delete;
                /**
                 * Singletons 应该是不可转让的。
                 */
                void operator=(const Singleton &) = delete;
                /**
                 * 这是一个静态方法，控制对单子实例的访问。在第一次运行时，它创建一个单子对象并将其放入静态字段。在随后的运行中，它返回存储在静态字段中的客户现有对象。
                 */

                static Singleton *GetInstance(const std::string &value);
                /**
                 * 最后，任何单子都应该定义一些业务逻辑，这些逻辑可以在其实例上执行。
                 */
                void SomeBusinessLogic()
                {
                    // ...
                }

                std::string value() const
                {
                    return value_;
                }
            };

            Singleton *Singleton::singleton_ = nullptr;
            ;

            /**
             * 静态方法应该在类之外定义。
             */
            Singleton *Singleton::GetInstance(const std::string &value)
            {
                /**
                 * 这是一种更安全的创建实例的方式。 instance =
                 * new
                 * Singleton是危险的，以防两个实例线程想同时访问。
                 */
                if (singleton_ == nullptr)
                {
                    singleton_ = new Singleton(value);
                }
                return singleton_;
            }

            void ThreadFoo()
            {
                // 以下代码模拟了缓慢的初始化。
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                Singleton *singleton = Singleton::GetInstance("FOO");
                std::cout << singleton->value() << "\n";
            }

            void ThreadBar()
            {
                // 以下代码模拟了缓慢的初始化。
                std::this_thread::sleep_for(std::chrono::milliseconds(1000));
                Singleton *singleton = Singleton::GetInstance("BAR");
                std::cout << singleton->value() << "\n";
            }

        }
    }
}
