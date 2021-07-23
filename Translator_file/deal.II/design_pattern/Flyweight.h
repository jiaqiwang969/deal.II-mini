#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

namespace DesignPattern
{
    namespace StructuralDesign
    {
        /**
         * Flyweight
         *
         * - 模式动机: 通过在多个对象之间共享状态的共同部分，而不是在每个对象中保留所有的数据，让你在可用的RAM数量中容纳更多的对象。
         *
         * 面向对象技术可以很好地解决一些灵活性或可扩展性问题，但在很多情况下需要在系统中增加类和对象的个数。当对象数量太多时，将导致运行代价过高，带来性能下降等问题。
         *
         * 享元模式正是为解决这一类问题而诞生的。享元模式通过共享技术实现相同或相似对象的重用。
         *
         * 在享元模式中可以共享的相同内容称为内部状态(IntrinsicState)，而那些需要外部环境来设置的不能共享的内容称为外部状态(Extrinsic
         * State)，由于区分了内部状态和外部状态，因此可以通过设置不同的外部状态使得相同的对象可以具有一些不同的特征，而相同的内部状态是可以共享的。
         *
         * 在享元模式中通常会出现工厂模式，需要创建一个享元工厂来负责维护一个享元池(Flyweight
         * Pool)用于存储具有相同内部状态的享元对象。
         *
         * 在享元模式中共享的是享元对象的内部状态，外部状态需要通过环境来设置。在实际使用中，能够共享的内部状态是有限的，因此享元对象一般都设计为较小的对象，它所包含的内部状态较少，这种对象也称为细粒度对象。享元模式的目的就是使用共享技术来实现大量细粒度对象的复用。
         *
         */

        /**
         * - 模式定义: 享元模式(Flyweight
         * Pattern)：运用共享技术有效地支持大量细粒度对象的复用。系统只使用少量的对象，而这些对象都很相似，状态变化很小，可以实现对象的多次复用。由于享元模式要求能够共享的对象必须是细粒度对象，因此它又称为轻量级模式，它是一种对象结构型模式。
         */

        namespace Flyweight
        {

            struct SharedState
            {
                std::string brand_;
                std::string model_;
                std::string color_;

                SharedState(const std::string &brand, const std::string &model, const std::string &color)
                    : brand_(brand), model_(model), color_(color)
                {
                }

                friend std::ostream &operator<<(std::ostream &os, const SharedState &ss)
                {
                    return os << "[ " << ss.brand_ << " , " << ss.model_ << " , " << ss.color_ << " ]";
                }
            };

            struct UniqueState
            {
                std::string owner_;
                std::string plates_;

                UniqueState(const std::string &owner, const std::string &plates)
                    : owner_(owner), plates_(plates)
                {
                }

                friend std::ostream &operator<<(std::ostream &os, const UniqueState &us)
                {
                    return os << "[ " << us.owner_ << " , " << us.plates_ << " ]";
                }
            };

            /**
             * Flyweight存储了属于多个真实商业实体的共同部分的状态（也叫内在状态）。Flyweight通过其方法参数接受其余的状态（外在状态，对每个实体都是独一无二的）。
             */
            class Flyweight
            {
            private:
                SharedState *shared_state_;

            public:
                Flyweight(const SharedState *shared_state) : shared_state_(new SharedState(*shared_state))
                {
                }
                Flyweight(const Flyweight &other) : shared_state_(new SharedState(*other.shared_state_))
                {
                }
                ~Flyweight()
                {
                    delete shared_state_;
                }
                SharedState *shared_state() const
                {
                    return shared_state_;
                }
                void Operation(const UniqueState &unique_state) const
                {
                    std::cout << "Flyweight: Displaying shared (" << *shared_state_ << ") and unique (" << unique_state << ") state.\n";
                }
            };
            /**
             * Flyweight Factory 创建并管理 Flyweight
             * 对象。它确保了飞翔权的正确共享。当客户端请求一个重量级对象时，工厂要么返回一个现有的实例，要么创建一个新的，如果它还不存在的话。
             */
            class FlyweightFactory
            {
                /**
                 * @var Flyweight[]
                 */
            private:
                std::unordered_map<std::string, Flyweight> flyweights_;
                /**
                 * 返回一个给定状态的Flyweight的字符串哈希值。
                 */
                std::string GetKey(const SharedState &ss) const
                {
                    return ss.brand_ + "_" + ss.model_ + "_" + ss.color_;
                }

            public:
                FlyweightFactory(std::initializer_list<SharedState> share_states)
                {
                    for (const SharedState &ss : share_states)
                    {
                        this->flyweights_.insert(std::make_pair<std::string, Flyweight>(this->GetKey(ss), Flyweight(&ss)));
                    }
                }

                /**
                 * 返回一个具有给定状态的现有Flyweight或创建一个新的。
                 */
                Flyweight GetFlyweight(const SharedState &shared_state)
                {
                    std::string key = this->GetKey(shared_state);
                    if (this->flyweights_.find(key) == this->flyweights_.end())
                    {
                        std::cout << "FlyweightFactory: Can't find a flyweight, creating new one.\n";
                        this->flyweights_.insert(std::make_pair(key, Flyweight(&shared_state)));
                    }
                    else
                    {
                        std::cout << "FlyweightFactory: Reusing existing flyweight.\n";
                    }
                    return this->flyweights_.at(key);
                }
                void ListFlyweights() const
                {
                    size_t count = this->flyweights_.size();
                    std::cout << "\nFlyweightFactory: I have " << count << " flyweights:\n";
                    for (std::pair<std::string, Flyweight> pair : this->flyweights_)
                    {
                        std::cout << pair.first << "\n";
                    }
                }
            };

            // ...
            void AddCarToPoliceDatabase(
                FlyweightFactory &ff, const std::string &plates, const std::string &owner,
                const std::string &brand, const std::string &model, const std::string &color)
            {
                std::cout << "\nClient: Adding a car to database.\n";
                const Flyweight &flyweight = ff.GetFlyweight({brand, model, color});
                // The client code either stores or calculates extrinsic state and passes it
                // to the flyweight's methods.
                flyweight.Operation({owner, plates});
            }

        }
    }
}

/**
 * 客户端代码通常会在应用程序的初始化阶段创建一堆预先填充的飞轮。
 */

// int main()
// {
//     FlyweightFactory *factory = new FlyweightFactory({{"Chevrolet", "Camaro2018", "pink"}, {"Mercedes Benz", "C300", "black"}, {"Mercedes Benz", "C500", "red"}, {"BMW", "M5", "red"}, {"BMW", "X6", "white"}});
//     factory->ListFlyweights();

//     AddCarToPoliceDatabase(*factory,
//                             "CL234IR",
//                             "James Doe",
//                             "BMW",
//                             "M5",
//                             "red");

//     AddCarToPoliceDatabase(*factory,
//                             "CL234IR",
//                             "James Doe",
//                             "BMW",
//                             "X1",
//                             "red");
//     factory->ListFlyweights();
//     delete factory;

//     return 0;
// }
