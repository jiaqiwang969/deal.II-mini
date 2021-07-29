
#include <iostream>
#include <mutex>
#include <thread>

namespace DesignPattern
{
    /**
     * Singleton
     *
     * - 模式动机:
     *
     * 让你确保一个类只有一个实例，同时为这个实例提供一个全局访问点。
     */
    /**
     * - 模式定义:
     *
     * Singleton类定义了
     * "GetInstance"方法，作为构造函数的替代方法，让客户反复访问该类的同一实例。
     */

    namespace Singleton
    {
        class Singleton
        {

            /**
             * Singleton的构造器/解构器应该总是私有的，以防止用`new`/`delete`操作符直接构造/解构调用。
             */
        private:
            static Singleton *pinstance_;
            static std::mutex mutex_;

        protected:
            Singleton(const std::string value) : value_(value)
            {
            }
            ~Singleton() {}
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
             * 最后，任何单子都应该定义一些业务逻辑，可以在其实例上执行。
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

        /**
         * 静态方法应该在类之外定义。
         */

        Singleton *Singleton::pinstance_{nullptr};
        std::mutex Singleton::mutex_;

        /**
         * 第一次调用GetInstance时，我们将锁定存储位置，然后再次确保变量为空，然后再设置值。
         * RU:
         */
        Singleton *Singleton::GetInstance(const std::string &value)
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (pinstance_ == nullptr)
            {
                pinstance_ = new Singleton(value);
            }
            return pinstance_;
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
