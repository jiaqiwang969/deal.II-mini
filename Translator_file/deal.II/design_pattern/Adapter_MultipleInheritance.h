#include <algorithm>
#include <iostream>
#include <string>
namespace DesignPattern
{

    /**
 * Adapter Design Pattern
 *
 * 意图: 提供一个统一的接口，允许具有不兼容接口的对象进行协作。
 */
    namespace Adapter
    {

        /**
 * 目标定义了客户端代码所使用的特定领域的接口。
 */
        class Target
        {
        public:
            virtual ~Target() = default;
            virtual std::string Request() const
            {
                return "Target: The default target's behavior.";
            }
        };

        /**
 * Adaptee包含一些有用的行为，但它的接口与现有的客户端代码不兼容。在客户端代码能够使用它之前，适应者需要进行一些调整。
 */
        class Adaptee
        {
        public:
            std::string SpecificRequest() const
            {
                return ".eetpadA eht fo roivaheb laicepS";
            }
        };

        /**
 * 适配器使用多重继承使适应者的接口与目标的接口兼容。
 */
        class Adapter : public Target, public Adaptee
        {
        public:
            Adapter() {}
            std::string Request() const override
            {
                std::string to_reverse = SpecificRequest();
                std::reverse(to_reverse.begin(), to_reverse.end());
                return "Adapter: (TRANSLATED) " + to_reverse;
            }
        };

        /**
 * 客户端代码支持所有遵循目标接口的类。
 */
        void ClientCode(const Target *target)
        {
            std::cout << target->Request();
        }
    }
}
