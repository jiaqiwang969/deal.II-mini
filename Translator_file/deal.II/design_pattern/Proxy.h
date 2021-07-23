#include <iostream>
/**
 * Proxy Design Pattern
 *
 * 意图: 为另一个对象提供一个代理或占位符，以控制对原始对象的访问或增加其他责任。
 */
/**
 * Subject接口声明了RealSubject和Proxy的共同操作。只要客户端使用这个接口与RealSubject一起工作，你就能给它传递一个代理而不是一个真正的主体。
 */
class Subject
{
public:
    virtual void Request() const = 0;
};
/**
 * RealSubject包含一些核心业务逻辑。通常，RealSubjects能够做一些有用的工作，这些工作也可能非常缓慢或敏感--例如，纠正输入数据。代理可以解决这些问题，而无需对RealSubject的代码进行任何修改。
 */
class RealSubject : public Subject
{
public:
    void Request() const override
    {
        std::cout << "RealSubject: Handling request.\n";
    }
};
/**
 * 代理人有一个与RealSubject相同的接口。
 */
class Proxy : public Subject
{
    /**
   * @var RealSubject
   */
private:
    RealSubject *real_subject_;

    bool CheckAccess() const
    {
        // 一些真正的检查应该在这里进行。
        std::cout << "Proxy: Checking access prior to firing a real request.\n";
        return true;
    }
    void LogAccess() const
    {
        std::cout << "Proxy: Logging the time of request.\n";
    }

    /**
   * 代理人维护对RealSubject类的一个对象的引用。它可以被懒惰地加载或由客户传递给代理。
   */
public:
    Proxy(RealSubject *real_subject) : real_subject_(new RealSubject(*real_subject))
    {
    }

    ~Proxy()
    {
        delete real_subject_;
    }
    /**
   *  代理模式最常见的应用是懒惰加载、缓存、控制访问、记录等。一个代理可以执行这些事情中的一个，然后根据结果，将执行结果传递给链接的RealSubject对象中的同一个方法。
   */
    void Request() const override
    {
        if (this->CheckAccess())
        {
            this->real_subject_->Request();
            this->LogAccess();
        }
    }
};
/**
 * 客户端代码应该通过主体接口与所有对象（包括主体和代理）一起工作，以便支持真实的主体和代理。然而，在现实生活中，客户端大多直接与他们的真实主体工作。在这种情况下，为了更容易实现该模式，你可以从真实主体的类中扩展你的代理。
 */
void ClientCode(const Subject &subject)
{
    // ...
    subject.Request();
    // ...
}

// int main() {
//   std::cout << "Client: Executing the client code with a real subject:\n";
//   RealSubject *real_subject = new RealSubject;
//   ClientCode(*real_subject);
//   std::cout << "\n";
//   std::cout << "Client: Executing the same client code with a proxy:\n";
//   Proxy *proxy = new Proxy(real_subject);
//   ClientCode(*proxy);

//   delete real_subject;
//   delete proxy;
//   return 0;
// }
