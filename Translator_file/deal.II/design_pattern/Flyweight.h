#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

/**
 * Flyweight Design Pattern
 *
 * 意图: 通过在多个对象之间共享状态的共同部分，而不是在每个对象中保留所有的数据，让你在可用的RAM数量中容纳更多的对象。
 */
namespace DesignPattern
{

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
 * Flyweight Factory 创建并管理 Flyweight 对象。它确保了飞翔权的正确共享。当客户端请求一个重量级对象时，工厂要么返回一个现有的实例，要么创建一个新的，如果它还不存在的话。
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