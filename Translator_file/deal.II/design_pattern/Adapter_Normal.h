#include <algorithm>
#include <iostream>
#include <string>

/**
 * Adapter
 *
 * - 模式动机:
 * 
 *   提供一个统一的接口，允许具有不兼容接口的对象进行协作。
 * 
 *   在软件开发中采用类似于电源适配器的设计和编码技巧被称为适配器模式。
 * 
 *   通常情况下，客户端可以通过目标类的接口访问它所提供的服务。有时，现有的类可以满足客户类的功能需要，但是它所提供的接口不一定是客户类所期望的，这可能是因为现有类中方法名与目标类中定义的方法名不一致等原因所导致的。
 * 
 *   在这种情况下，现有的接口需要转化为客户类期望的接口，这样保证了对现有类的重用。如果不进行这样的转化，客户类就不能利用现有类所提供的功能，适配器模式可以完成这样的转化。
 * 
 *   在适配器模式中可以定义一个包装类，包装不兼容接口的对象，这个包装类指的就是适配器(Adapter)，它所包装的对象就是适配者(Adaptee)，即被适配的类。
 * 
 *   适配器提供客户类需要的接口，适配器的实现就是把客户类的请求转化为对适配者的相应接口的调用。也就是说：当客户类调用适配器的方法时，在适配器类的内部将调用适配者类的方法，而这个过程对客户类是透明的，客户类并不直接访问适配者类。因此，适配器可以使由于接口不兼容而不能交互的类可以一起工作。这就是适配器模式的模式动机。
 */

/**
 * - 模式定义: 适配器模式(Adapter Pattern)
 * 
 *   将一个接口转换成客户希望的另一个接口，适配器模式使接口不兼容的那些类可以一起工作，其别名为包装器(Wrapper)。适配器模式既可以作为类结构型模式，也可以作为对象结构型模式。
 * 
 *   目标定义了客户端代码所使用的特定领域的接口。
 */
namespace DesignPattern
{
    namespace StructuralDesign
    {

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
         * 适配器使适应者的接口与目标的接口兼容。
         */
            class Adapter : public Target
            {
            private:
                Adaptee *adaptee_;

            public:
                Adapter(Adaptee *adaptee) : adaptee_(adaptee) {}
                std::string Request() const override
                {
                    std::string to_reverse = this->adaptee_->SpecificRequest();
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

}