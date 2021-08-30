/**
@page step_0a The step-0a tutorial program
@htmlonly
<table class="tutorial" width="50%">
<tr><th colspan="2"><b><small>Table of contents</small></b></th></tr>
<tr><td width="50%" valign="top">
<ol>
  <li> <a href="#Intro" class=bold>Introduction</a>
    <ul>
    </ul>
  <li> <a href="#CommProg" class=bold>The commented program</a>
    <ul>
      </ul>
</ol></td><td width="50%" valign="top"><ol>
  <li value="3"> <a href="#Results" class=bold>Results</a>
    <ul>
    </ul>
  <li> <a href="#PlainProg" class=bold>The plain program</a>
</ol> </td> </tr> </table>
@endhtmlonly
<a name="23个设计模式"></a>

TODO: doxygen it!

参考: 
- 《图说设计模式》- https://github.com/me115/design_patterns
- 《深入设计模式》- Alexander Shvets
- 《设计模式》- 4人组
-  博士论文 - A Framework for Developing Finite Element Codes for MultiDisciplinary Applications 



<a name=""></a><h1> 一. 创建型模式 </h1>

创建型模式提供了创建对象的机制,能够提升已有代码的灵活性和可复用性。

<a name="1FactoryMethod"></a><h2> 1. 工厂方法 Factory Method </h2>

在父类中提供一个创建对象的接口以允许子类决定实例化对象的类型。

优缺点
- 你可以避免创建者和具体产品之间的紧密耦合。
- 单一职责原则。你可以将产品创建代码放在程序的单一位置, 从而使得代码更容易维护。
- 开闭原则。无需更改现有客户端代码,你就可以在程序中引入新的产品类型。
- 应用工厂方法模式需要引入许多新的子类,代码可能会因此变得更复杂。最好的情况是将该模式引入创建者类的现有层次结构中。


@code
#include <iostream>
#include <string>

/**
 * 工厂设计方法
 *
 * 意图：为在父类中创建对象提供一个接口，但允许子类改变将被创建的对象的类型。
 */

/**
 * 产品接口声明了所有具体产品必须 实现的操作。
 */

class Product {
 public:
  virtual ~Product() {}
  virtual std::string Operation() const = 0;
};

/**
 * 具体的产品提供了产品接口的各种实现。
 */
class ConcreteProduct1 : public Product {
 public:
  std::string Operation() const override {
    return "{Result of the ConcreteProduct1}";
  }
};
class ConcreteProduct2 : public Product {
 public:
  std::string Operation() const override {
    return "{Result of the ConcreteProduct2}";
  }
};

/**
 * Creator类声明了工厂方法，它应该返回一个产品类的对象。创造者的子类通常提供这个方法的实现。
 */

class Creator {
  /**
   * 请注意，创造者也可能提供一些工厂方法的默认实现。
   */
 public:
  virtual ~Creator(){};
  virtual Product* FactoryMethod() const = 0;
  /**
   * 还要注意的是，尽管它的名字叫 "创造者"，但它的主要责任不是创造产品。通常，它包含一些依赖于产品对象的核心业务逻辑，这些对象由工厂方法返回。子类可以通过重写工厂方法并从中返回不同类型的产品来间接地改变该业务逻辑。
   */

  std::string SomeOperation() const {
    // Call the factory method to create a Product object.
    Product* product = this->FactoryMethod();
    // Now, use the product.
    std::string result = "Creator: The same creator's code has just worked with " + product->Operation();
    delete product;
    return result;
  }
};

/**
 * 具体的创建者覆盖工厂方法，以改变生成的产品的类型。
 */
class ConcreteCreator1 : public Creator {
  /**
   *  请注意，方法的签名仍然使用抽象产品类型，尽管具体产品实际上是从该方法中返回的。这样，创造者就可以保持独立于具体的产品类型。
   */
 public:
  Product* FactoryMethod() const override {
    return new ConcreteProduct1();
  }
};

class ConcreteCreator2 : public Creator {
 public:
  Product* FactoryMethod() const override {
    return new ConcreteProduct2();
  }
};

/**
 * 客户端代码与一个具体的创造者的实例一起工作，尽管是通过它的基础接口。只要客户端继续通过基接口与创造者一起工作，你就可以把任何创造者的子类传递给它。
 */
void ClientCode(const Creator& creator) {
  // ...
  std::cout << "Client: I'm not aware of the creator's class, but it still works.\n"
            << creator.SomeOperation() << std::endl;
  // ...
}

/**
 * 应用程序根据配置或环境来挑选创建者的类型。
 */

int main() {
  std::cout << "App: Launched with the ConcreteCreator1.\n";
  Creator* creator = new ConcreteCreator1();
  ClientCode(*creator);
  std::cout << std::endl;
  std::cout << "App: Launched with the ConcreteCreator2.\n";
  Creator* creator2 = new ConcreteCreator2();
  ClientCode(*creator2);

  delete creator;
  delete creator2;
  return 0;
}
@endcode

<a name="2AbstractFactory"></a><h2> 2. 抽象工厂 Abstract Factory </h2>

让你能创建一系列相关的对象,而无需指定其具体类。
@code
#include <iostream>
#include <string>

/**
 * Abstract Factory
 *
 * - 模式动机: 让你产生相关对象的系列，而不需要指定它们的具体类别。
 */

/**
 * 一个产品系列的每个不同产品都应该有一个基本接口。该产品的所有变体都必须实现这个接口。
 */
class AbstractProductA {
 public:
  virtual ~AbstractProductA(){};
  virtual std::string UsefulFunctionA() const = 0;
};

/**
 * Concrete产品是由相应的Concrete工厂创造的。
 */
class ConcreteProductA1 : public AbstractProductA {
 public:
  std::string UsefulFunctionA() const override {
    return "The result of the product A1.";
  }
};

class ConcreteProductA2 : public AbstractProductA {
  std::string UsefulFunctionA() const override {
    return "The result of the product A2.";
  }
};

/**
 * 这里是另一个产品的基本接口。所有的产品都可以相互作用，但只有在相同的具体变体的产品之间才能进行适当的互动。
 */
class AbstractProductB {
  /**
   * 产品B能够做它自己的事情...
   */
 public:
  virtual ~AbstractProductB(){};
  virtual std::string UsefulFunctionB() const = 0;
  /**
   *   ...但它也可以与ProductA协作。
   * 抽象工厂确保它所创建的所有产品都具有相同的变体，因此是兼容的。
   */
  virtual std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const = 0;
};

/**
 * Concrete产品是由相应的Concrete工厂创造的。
 */
class ConcreteProductB1 : public AbstractProductB {
 public:
  std::string UsefulFunctionB() const override {
    return "The result of the product B1.";
  }
  /**
   * 变体 "产品B1 "只能够与变体 "产品A1 "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
   */
  std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override {
    const std::string result = collaborator.UsefulFunctionA();
    return "The result of the B1 collaborating with ( " + result + " )";
  }
};

class ConcreteProductB2 : public AbstractProductB {
 public:
  std::string UsefulFunctionB() const override {
    return "The result of the product B2.";
  }
  /**
   * 变体 "Product B2 "只能与变体 "Product A2 "正常工作。然而，它接受AbstractProductA的任何实例作为参数。
   */
  std::string AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override {
    const std::string result = collaborator.UsefulFunctionA();
    return "The result of the B2 collaborating with ( " + result + " )";
  }
};

/**
 * 抽象工厂接口声明了一组返回不同抽象产品的方法。这些产品被称为一个家族，并通过一个高级的主题或概念进行关联。一个家族的产品通常能够在它们之间进行协作。一个家族的产品可能有几个变体，但一个变体的产品与另一个变体的产品是不兼容的。
 */
class AbstractFactory {
 public:
  virtual AbstractProductA *CreateProductA() const = 0;
  virtual AbstractProductB *CreateProductB() const = 0;
};

/**
 * Concrete工厂生产属于单一变体的产品系列。工厂保证产生的产品是兼容的。请注意，具体工厂的方法的签名会返回一个抽象产品，而在方法内部会实例化一个具体产品。
 */
class ConcreteFactory1 : public AbstractFactory {
 public:
  AbstractProductA *CreateProductA() const override {
    return new ConcreteProductA1();
  }
  AbstractProductB *CreateProductB() const override {
    return new ConcreteProductB1();
  }
};

/**
 * 每个Concrete工厂都有一个相应的产品变体。
 */
class ConcreteFactory2 : public AbstractFactory {
 public:
  AbstractProductA *CreateProductA() const override {
    return new ConcreteProductA2();
  }
  AbstractProductB *CreateProductB() const override {
    return new ConcreteProductB2();
  }
};

/**
 * 客户端代码只通过抽象类型与工厂和产品一起工作。AbstractFactory和AbstractProduct。这让你可以将任何工厂或产品的子类传递给客户端代码而不破坏它。
 */

void ClientCode(const AbstractFactory &factory) {
  const AbstractProductA *product_a = factory.CreateProductA();
  const AbstractProductB *product_b = factory.CreateProductB();
  std::cout << product_b->UsefulFunctionB() << "\n";
  std::cout << product_b->AnotherUsefulFunctionB(*product_a) << "\n";
  delete product_a;
  delete product_b;
}

int main() {
  std::cout << "Client: Testing client code with the first factory type:\n";
  ConcreteFactory1 *f1 = new ConcreteFactory1();
  ClientCode(*f1);
  delete f1;
  std::cout << std::endl;
  std::cout << "Client: Testing the same client code with the second factory type:\n";
  ConcreteFactory2 *f2 = new ConcreteFactory2();
  ClientCode(*f2);
  delete f2;
  return 0;
}
@endcode


<a name="3Builder"></a><h2> 3. 生成器 Builder </h2>

使你能够分步骤创建复杂对象。该模式允许你使用相同的创建代码生成不同类型和形式的对象。
@code
#include <iostream>
#include <string>
#include <vector>

/**
 * Builder
 *
 * - 模式动机: 让你一步一步地构建复杂的对象。该模式允许你使用相同的构造代码产生不同类型和表现形式的对象。
 */

/**
 * 只有当你的产品相当复杂并需要大量配置时，使用构建器模式才有意义。
 * 与其他创建模式不同，不同的具体构建器可以产生不相关的产品。换句话说，各种构建器的结果可能并不总是遵循相同的接口。
 */

class Product1{
    public:
    std::vector<std::string> parts_;
    void ListParts()const{
        std::cout << "Product parts: ";
        for (size_t i=0;i<parts_.size();i++){
            if(parts_[i]== parts_.back()){
                std::cout << parts_[i];
            }else{
                std::cout << parts_[i] << ", ";
            }
        }
        std::cout << "\n\n"; 
    }
};


/**
 * 生成器接口指定了用于创建产品对象不同部分的方法。
 */
class Builder{
    public:
    virtual ~Builder(){}
    virtual void ProducePartA() const =0;
    virtual void ProducePartB() const =0;
    virtual void ProducePartC() const =0;
};
/**
 * 具体的Builder类遵循Builder接口并提供具体的构建步骤的实现。你的程序可能有几个不同的Builder，实现方式也不同。
 */
class ConcreteBuilder1 : public Builder{
    private:

    Product1* product;

    /**
     * 一个新的构建器实例应该包含一个空白的产品对象，它被用于进一步的装配。
     */
    public:

    ConcreteBuilder1(){
        this->Reset();
    }

    ~ConcreteBuilder1(){
        delete product;
    }

    void Reset(){
        this->product= new Product1();
    }
    /**
     * All production steps work with the same product instance.
     */

    void ProducePartA()const override{
        this->product->parts_.push_back("PartA1");
    }

    void ProducePartB()const override{
        this->product->parts_.push_back("PartB1");
    }

    void ProducePartC()const override{
        this->product->parts_.push_back("PartC1");
    }

    /**
     * Concrete建造商应该提供他们自己的方法来检索结果。这是因为各种类型的构建器可能会创建完全不同的产品，它们不遵循相同的接口。因此，这种方法不能在基本的Builder接口中声明（至少在静态类型的编程语言中）。请注意，PHP是一种动态类型的语言，这种方法可以在基础接口中出现。然而，为了清楚起见，我们不会在那里声明它。
     * 通常情况下，在将最终结果返回给客户之后，一个构建器实例被期望准备好开始生产另一个产品。这就是为什么通常的做法是在`getProduct`方法体的末尾调用重置方法。然而，这种行为并不是强制性的，你可以让你的构建器在处理之前的结果之前等待来自客户端代码的明确的重置调用。
     */

    /**
     * 请注意这里的内存所有权。一旦你调用GetProduct，这个函数的用户就有责任释放这段内存。这里可能是一个更好的选择，即使用智能指针来避免内存泄漏。
     */

    Product1* GetProduct() {
        Product1* result= this->product;
        this->Reset();
        return result;
    }
};

/**
 * 主任只负责按特定顺序执行建造步骤。当按照特定的顺序或配置生产产品时，它很有帮助。严格来说，Director类是可选的，因为客户可以直接控制建造者。
 */
class Director{
    /**
     * @var Builder
     */
    private:
    Builder* builder;
    /**
     * Director与客户代码传递给它的任何构建器实例一起工作。这样，客户代码可以改变新组装产品的最终类型。
     */

    public:

    void set_builder(Builder* builder){
        this->builder=builder;
    }

    /**
     * 主任可以使用相同的建造步骤建造多个产品变体。
     */

    void BuildMinimalViableProduct(){
        this->builder->ProducePartA();
    }
    
    void BuildFullFeaturedProduct(){
        this->builder->ProducePartA();
        this->builder->ProducePartB();
        this->builder->ProducePartC();
    }
};
/**
 * 客户端代码创建一个构建器对象，将其传递给主任，然后启动构建过程。最终结果从建造者对象中获取。
 */
/**
 * 为了简单起见，我使用了原始指针，但是你可能更喜欢在这里使用智能指针。
 */
void ClientCode(Director& director)
{
    ConcreteBuilder1* builder = new ConcreteBuilder1();
    director.set_builder(builder);
    std::cout << "Standard basic product:\n"; 
    director.BuildMinimalViableProduct();
    
    Product1* p= builder->GetProduct();
    p->ListParts();
    delete p;

    std::cout << "Standard full featured product:\n"; 
    director.BuildFullFeaturedProduct();

    p= builder->GetProduct();
    p->ListParts();
    delete p;

    // Remember, the Builder pattern can be used without a Director class.
    std::cout << "Custom product:\n";
    builder->ProducePartA();
    builder->ProducePartC();
    p=builder->GetProduct();
    p->ListParts();
    delete p;

    delete builder;
}

int main(){
    Director* director= new Director();
    ClientCode(*director);
    delete director;
    return 0;    
}

@endcode



<a name="4Prototype"></a><h2> 4. 原型 Prototype </h2>

让你能够复制已有对象,而又无需使代码依赖它们所属的类。
@code
#include <iostream>
#include <string>
#include <unordered_map>

using std::string;

/**
 * Prototype
 *  模式动机: 让你复制现有的对象而不使你的代码依赖于它们的类。
 */

enum Type {
  PROTOTYPE_1 = 0,
  PROTOTYPE_2
};

/**
 * 具有克隆能力的实例类。我们将看到不同类型的字段的值是如何被克隆的。
 */

class Prototype {
 protected:
  string prototype_name_;
  float prototype_field_;

 public:
  Prototype() {}
  Prototype(string prototype_name)
      : prototype_name_(prototype_name) {
  }
  virtual ~Prototype() {}
  virtual Prototype *Clone() const = 0;
  virtual void Method(float prototype_field) {
    this->prototype_field_ = prototype_field;
    std::cout << "Call Method from " << prototype_name_ << " with field : " << prototype_field << std::endl;
  }
};

/**
 * ConcretePrototype1是Prototype的一个子类，实现了Clone方法 在这个例子中，Prototype类的所有数据成员都在堆栈中。如果你的属性中有指针，例如：String* name_ ，你将需要实现Copy-Constructor以确保你有一个来自clone方法的深度拷贝。
 */

class ConcretePrototype1 : public Prototype {
 private:
  float concrete_prototype_field1_;

 public:
  ConcretePrototype1(string prototype_name, float concrete_prototype_field)
      : Prototype(prototype_name), concrete_prototype_field1_(concrete_prototype_field) {
  }

  /**
   * 注意，Clone方法返回一个指向新的ConcretePrototype1副本的指针。因此，客户端（调用Clone方法的人）有责任释放该内存。如果你有智能指针的知识，你可能更喜欢在这里使用unique_pointer。
   */
  Prototype *Clone() const override {
    return new ConcretePrototype1(*this);
  }
};

class ConcretePrototype2 : public Prototype {
 private:
  float concrete_prototype_field2_;

 public:
  ConcretePrototype2(string prototype_name, float concrete_prototype_field)
      : Prototype(prototype_name), concrete_prototype_field2_(concrete_prototype_field) {
  }
  Prototype *Clone() const override {
    return new ConcretePrototype2(*this);
  }
};

/**
 * 在PrototypeFactory中，你有两个具体的原型，每个具体的原型类都有一个，所以每次你想创建一个子弹，你可以使用现有的原型并克隆这些原型。
 */

class PrototypeFactory {
 private:
  std::unordered_map<Type, Prototype *, std::hash<int>> prototypes_;

 public:
  PrototypeFactory() {
    prototypes_[Type::PROTOTYPE_1] = new ConcretePrototype1("PROTOTYPE_1 ", 50.f);
    prototypes_[Type::PROTOTYPE_2] = new ConcretePrototype2("PROTOTYPE_2 ", 60.f);
  }

  /**
   * 要注意释放所有分配的内存。同样，如果你有智能指针，那么在这里使用它将会更好。
   */

  ~PrototypeFactory() {
    delete prototypes_[Type::PROTOTYPE_1];
    delete prototypes_[Type::PROTOTYPE_2];
  }

  /**
   * Notice here that you just need to specify the type of the prototype you
   * want and the method will create from the object with this type.
   */
  Prototype *CreatePrototype(Type type) {
    return prototypes_[type]->Clone();
  }
};

void Client(PrototypeFactory &prototype_factory) {
  std::cout << "Let's create a Prototype 1\n";

  Prototype *prototype = prototype_factory.CreatePrototype(Type::PROTOTYPE_1);
  prototype->Method(90);
  delete prototype;

  std::cout << "\n";

  std::cout << "Let's create a Prototype 2 \n";

  prototype = prototype_factory.CreatePrototype(Type::PROTOTYPE_2);
  prototype->Method(10);

  delete prototype;
}

int main() {
  PrototypeFactory *prototype_factory = new PrototypeFactory();
  Client(*prototype_factory);
  delete prototype_factory;

  return 0;
}
@endcode
<a name="5Singleton"></a><h2> 5. 单例 Singleton </h2>

让你能够保证一个类只有一个实例,并提供一个访问该实例的全局节点。

NonThreadSafe
@code 
#include <iostream>
#include <string>
#include <thread>

/**
 * Singleton
 *
 * - 模式动机: 让你确保一个类只有一个实例，同时为这个实例提供一个全局访问点。
 */
/**
 * Singleton类定义了 "GetInstance "方法，作为构造函数的替代方法，让客户反复访问该类的同一实例。
 */
class Singleton
{

    /**
     * Singleton的构造函数应该始终是私有的，以防止用`new`操作符直接调用构造。
     */

protected:
    Singleton(const std::string value): value_(value)
    {
    }

    static Singleton* singleton_;

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

    static Singleton *GetInstance(const std::string& value);
    /**
     * 最后，任何单子都应该定义一些业务逻辑，这些逻辑可以在其实例上执行。
     */
    void SomeBusinessLogic()
    {
        // ...
    }

    std::string value() const{
        return value_;
    } 
};

Singleton* Singleton::singleton_= nullptr;;

/**
 * 静态方法应该在类之外定义。
 */
Singleton *Singleton::GetInstance(const std::string& value)
{
    /**
     *  这是一种更安全的创建实例的方式。 instance = new Singleton是危险的，以防两个实例线程想同时访问。
     */
    if(singleton_==nullptr){
        singleton_ = new Singleton(value);
    }
    return singleton_;
}

void ThreadFoo(){
    // 以下代码模拟了缓慢的初始化。
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    Singleton* singleton = Singleton::GetInstance("FOO");
    std::cout << singleton->value() << "\n";
}

void ThreadBar(){
    // 以下代码模拟了缓慢的初始化。
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    Singleton* singleton = Singleton::GetInstance("BAR");
    std::cout << singleton->value() << "\n";
}


int main()
{
    std::cout <<"If you see the same value, then singleton was reused (yay!\n" <<
                "If you see different values, then 2 singletons were created (booo!!)\n\n" <<
                "RESULT:\n";   
    std::thread t1(ThreadFoo);
    std::thread t2(ThreadBar);
    t1.join();
    t2.join();

    return 0;
}
@endcode


ThreadSafe
@code
/**
 * 请记住，这只是一个说明性的小例子，在现实世界中，你可能会想到一些更多可能的问题。
 */

#include <iostream>
#include <mutex>
#include <thread>

/**
 * Singleton
 *
 * 意图。让你确保一个类只有一个实例，同时为这个实例提供一个全局访问点。
 */
/**
 * Singleton类定义了 "GetInstance "方法，作为构造函数的替代方法，让客户反复访问该类的同一实例。
 */
class Singleton
{

    /**
     * Singleton的构造器/解构器应该总是私有的，以防止用`new`/`delete`操作符直接构造/解构调用。
     */
private:
    static Singleton * pinstance_;
    static std::mutex mutex_;

protected:
    Singleton(const std::string value): value_(value)
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

    static Singleton *GetInstance(const std::string& value);
    /**
     * 最后，任何单子都应该定义一些业务逻辑，可以在其实例上执行。
     */
    void SomeBusinessLogic()
    {
        // ...
    }
    
    std::string value() const{
        return value_;
    } 
};

/**
 * 静态方法应该在类之外定义。
 */

Singleton* Singleton::pinstance_{nullptr};
std::mutex Singleton::mutex_;

/**
 * 第一次调用GetInstance时，我们将锁定存储位置，然后再次确保变量为空，然后再设置值。 RU:
 */
Singleton *Singleton::GetInstance(const std::string& value)
{
    std::lock_guard<std::mutex> lock(mutex_);
    if (pinstance_ == nullptr)
    {
        pinstance_ = new Singleton(value);
    }
    return pinstance_;
}

void ThreadFoo(){
    // 以下代码模拟了缓慢的初始化。
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    Singleton* singleton = Singleton::GetInstance("FOO");
    std::cout << singleton->value() << "\n";
}

void ThreadBar(){
    // 以下代码模拟了缓慢的初始化。
    std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    Singleton* singleton = Singleton::GetInstance("BAR");
    std::cout << singleton->value() << "\n";
}

int main()
{   
    std::cout <<"If you see the same value, then singleton was reused (yay!\n" <<
                "If you see different values, then 2 singletons were created (booo!!)\n\n" <<
                "RESULT:\n";   
    std::thread t1(ThreadFoo);
    std::thread t2(ThreadBar);
    t1.join();
    t2.join();
    
    return 0;
}
@endcode


<a name=""></a><h1> 二. 结构型模式 </h1>


<a name="6Adapter"></a><h2> 6. 适配器 Adapter </h2>

让接口不兼容的对象能够相互合作。

Normal
@code
#include <algorithm>
#include <iostream>
#include <string>

/**
 * Adapter
 *
 * - 模式动机: 提供一个统一的接口，允许具有不兼容接口的对象进行协作。
 */

/**
 * 目标定义了客户端代码所使用的特定领域的接口。
 */
class Target {
 public:
  virtual ~Target() = default;

  virtual std::string Request() const {
    return "Target: The default target's behavior.";
  }
};

/**
 * Adaptee包含一些有用的行为，但它的接口与现有的客户端代码不兼容。在客户端代码能够使用它之前，适应者需要进行一些调整。
 */
class Adaptee {
 public:
  std::string SpecificRequest() const {
    return ".eetpadA eht fo roivaheb laicepS";
  }
};

/**
 * 适配器使适应者的接口与目标的接口兼容。
 */
class Adapter : public Target {
 private:
  Adaptee *adaptee_;

 public:
  Adapter(Adaptee *adaptee) : adaptee_(adaptee) {}
  std::string Request() const override {
    std::string to_reverse = this->adaptee_->SpecificRequest();
    std::reverse(to_reverse.begin(), to_reverse.end());
    return "Adapter: (TRANSLATED) " + to_reverse;
  }
};

/**
 * 客户端代码支持所有遵循目标接口的类。
 */
void ClientCode(const Target *target) {
  std::cout << target->Request();
}

int main() {
  std::cout << "Client: I can work just fine with the Target objects:\n";
  Target *target = new Target;
  ClientCode(target);
  std::cout << "\n\n";
  Adaptee *adaptee = new Adaptee;
  std::cout << "Client: The Adaptee class has a weird interface. See, I don't understand it:\n";
  std::cout << "Adaptee: " << adaptee->SpecificRequest();
  std::cout << "\n\n";
  std::cout << "Client: But I can work with it via the Adapter:\n";
  Adapter *adapter = new Adapter(adaptee);
  ClientCode(adapter);
  std::cout << "\n";

  delete target;
  delete adaptee;
  delete adapter;

  return 0;
}

@endcode

MultipleInheritance
@code
#include <algorithm>
#include <iostream>
#include <string>

/**
 * Adapter
 *
 * - 模式动机: 提供一个统一的接口，允许具有不兼容接口的对象进行协作。
 */

/**
 * 目标定义了客户端代码所使用的特定领域的接口。
 */
class Target {
 public:
  virtual ~Target() = default;
  virtual std::string Request() const {
    return "Target: The default target's behavior.";
  }
};

/**
 * Adaptee包含一些有用的行为，但它的接口与现有的客户端代码不兼容。在客户端代码能够使用它之前，适应者需要进行一些调整。
 */
class Adaptee {
 public:
  std::string SpecificRequest() const {
    return ".eetpadA eht fo roivaheb laicepS";
  }
};

/**
 * 适配器使用多重继承使适应者的接口与目标的接口兼容。
 */
class Adapter : public Target, public Adaptee {
 public:
  Adapter() {}
  std::string Request() const override {
    std::string to_reverse = SpecificRequest();
    std::reverse(to_reverse.begin(), to_reverse.end());
    return "Adapter: (TRANSLATED) " + to_reverse;
  }
};

/**
 * 客户端代码支持所有遵循目标接口的类。
 */
void ClientCode(const Target *target) {
  std::cout << target->Request();
}

int main() {
  std::cout << "Client: I can work just fine with the Target objects:\n";
  Target *target = new Target;
  ClientCode(target);
  std::cout << "\n\n";
  Adaptee *adaptee = new Adaptee;
  std::cout << "Client: The Adaptee class has a weird interface. See, I don't understand it:\n";
  std::cout << "Adaptee: " << adaptee->SpecificRequest();
  std::cout << "\n\n";
  std::cout << "Client: But I can work with it via the Adapter:\n";
  Adapter *adapter = new Adapter;
  ClientCode(adapter);
  std::cout << "\n";

  delete target;
  delete adaptee;
  delete adapter;

  return 0;
}

@endcode

<a name="7Bridge"></a><h2> 7. 桥接 Bridge </h2>

可将一个大类或一系列紧密相关的类拆分为抽象和实现两个独立的层次结构,从而能在开发时分别使用。
@code
#include <iostream>
#include <string>

/**
 * Bridge
 *
 * - 模式动机: 让你把一个大类或一组密切相关的类分成两个独立的层次--抽象和实现--它们可以互相独立开发。
 *
 *               A
 *            /     \                        A         N
 *          Aa      Ab        ===>        /     \     / \
 *         / \     /  \                 Aa(N) Ab(N)  1   2
 *       Aa1 Aa2  Ab1 Ab2
 */

/**
 * 实现定义了所有实现类的接口。它不一定要与抽象的接口相匹配。事实上，这两个接口可以完全不同。通常情况下，实现接口只提供原始操作，而抽象接口在这些原始操作的基础上定义了更高层次的操作。
 */

class Implementation {
 public:
  virtual ~Implementation() {}
  virtual std::string OperationImplementation() const = 0;
};

/**
 * 每个具体的实现对应于一个特定的平台，并使用该平台的API实现实现接口。
 */
class ConcreteImplementationA : public Implementation {
 public:
  std::string OperationImplementation() const override {
    return "ConcreteImplementationA: Here's the result on the platform A.\n";
  }
};
class ConcreteImplementationB : public Implementation {
 public:
  std::string OperationImplementation() const override {
    return "ConcreteImplementationB: Here's the result on the platform B.\n";
  }
};

/**
 * 抽象（Abstraction）定义了两个类层次的 "控制 "部分的接口。它维护一个对实现层次结构的对象的引用，并将所有的实际工作委托给这个对象。
 */

class Abstraction {
  /**
   * @var Implementation
   */
 protected:
  Implementation* implementation_;

 public:
  Abstraction(Implementation* implementation) : implementation_(implementation) {
  }

  virtual ~Abstraction() {
  }

  virtual std::string Operation() const {
    return "Abstraction: Base operation with:\n" +
           this->implementation_->OperationImplementation();
  }
};
/**
 * 你可以扩展抽象而不改变实现类。
 */
class ExtendedAbstraction : public Abstraction {
 public:
  ExtendedAbstraction(Implementation* implementation) : Abstraction(implementation) {
  }
  std::string Operation() const override {
    return "ExtendedAbstraction: Extended operation with:\n" +
           this->implementation_->OperationImplementation();
  }
};

/**
 * 除了在初始化阶段，一个抽象对象与一个特定的实现对象相联系，客户机代码应该只依赖于抽象类。这样一来，客户端代码就可以支持任何抽象-实现的组合。
 */
void ClientCode(const Abstraction& abstraction) {
  // ...
  std::cout << abstraction.Operation();
  // ...
}
/**
 * 客户端代码应该能够与任何预先配置的抽象-实现组合一起工作。
 */

int main() {
  Implementation* implementation = new ConcreteImplementationA;
  Abstraction* abstraction = new Abstraction(implementation);
  ClientCode(*abstraction);
  std::cout << std::endl;
  delete implementation;
  delete abstraction;

  implementation = new ConcreteImplementationB;
  abstraction = new ExtendedAbstraction(implementation);
  ClientCode(*abstraction);

  delete implementation;
  delete abstraction;

  return 0;
}

@endcode
<a name="8Composite"></a><h2> 8. 组合 Composite </h2>

你可以使用它将对象组合成树状结构,并且能像使用独立对象一样使用它们。
@code

#include <algorithm>
#include <iostream>
#include <list>
#include <string>
/**
 * Composite
 *
 * - 模式动机: 让你把对象组成树状结构，然后像处理单个对象一样处理这些结构。
 */
/**
 * 基层的Component类声明了一个组合的简单和复杂对象的共同操作。
 */
class Component {
  /**
   * @var Component
   */
 protected:
  Component *parent_;
  /**
   *  作为选择，基础组件可以声明一个接口，用于设置和访问树状结构中组件的父级。它也可以为这些方法提供一些默认实现。
   */
 public:
  virtual ~Component() {}
  void SetParent(Component *parent) {
    this->parent_ = parent;
  }
  Component *GetParent() const {
    return this->parent_;
  }
  /**
   * 在某些情况下，在基础组件类中直接定义儿童管理操作将是有益的。这样，你就不需要向客户端代码暴露任何具体的组件类，即使是在对象树装配期间。缺点是，这些方法对于叶级组件来说将是空的。
   */
  virtual void Add(Component *component) {}
  virtual void Remove(Component *component) {}
  /**
   *   你可以提供一个方法，让客户端代码弄清一个组件是否可以产生子类。
   */
  virtual bool IsComposite() const {
    return false;
  }
  /**
   *   基组件可以实现一些默认的行为，或者将其留给具体的类（通过将包含该行为的方法声明为 "抽象的"）。
   */
  virtual std::string Operation() const = 0;
};
/**
 * 叶子类代表一个组合的终端对象。一个叶子不能有任何孩子。
 * 通常情况下，是叶子对象做实际工作，而复合对象只委托给它们的子构件。
 */
class Leaf : public Component {
 public:
  std::string Operation() const override {
    return "Leaf";
  }
};
/**
 * 复合类代表了可能有子类的复杂组件。通常情况下，复合对象将实际工作委托给它们的子代，然后将结果 "汇总"。
 */
class Composite : public Component {
  /**
   * @var \SplObjectStorage
   */
 protected:
  std::list<Component *> children_;

 public:
  /**
   *   一个复合对象可以将其他组件（包括简单的或复杂的）加入或移出其子列表。
   */
  void Add(Component *component) override {
    this->children_.push_back(component);
    component->SetParent(this);
  }
  /**
   * 请记住，这个方法删除了列表的指针，但没有释放内存，你应该手动操作，或者最好使用智能指针。
   */
  void Remove(Component *component) override {
    children_.remove(component);
    component->SetParent(nullptr);
  }
  bool IsComposite() const override {
    return true;
  }
  /**
   *  复合体以一种特殊的方式执行其主要逻辑。它递归地遍历其所有的子代，收集和计算其结果。由于复合体的子代将这些调用传递给他们的子代，以此类推，整个对象树被遍历的结果。
   */
  std::string Operation() const override {
    std::string result;
    for (const Component *c : children_) {
      if (c == children_.back()) {
        result += c->Operation();
      } else {
        result += c->Operation() + "+";
      }
    }
    return "Branch(" + result + ")";
  }
};
/**
 * The client code works with all of the components via the base interface.
 */
void ClientCode(Component *component) {
  // ...
  std::cout << "RESULT: " << component->Operation();
  // ...
}

/**
 * 由于儿童管理操作是在基础组件类中声明的，客户代码可以与任何组件一起工作，不管是简单的还是复杂的，都不需要依赖它们的具体类。
 */
void ClientCode2(Component *component1, Component *component2) {
  // ...
  if (component1->IsComposite()) {
    component1->Add(component2);
  }
  std::cout << "RESULT: " << component1->Operation();
  // ...
}

/**
 * 这样，客户端代码就可以支持简单的叶子组件...
 */

int main() {
  Component *simple = new Leaf;
  std::cout << "Client: I've got a simple component:\n";
  ClientCode(simple);
  std::cout << "\n\n";
  /**
   * ...以及复杂的复合材料。
   */

  Component *tree = new Composite;
  Component *branch1 = new Composite;

  Component *leaf_1 = new Leaf;
  Component *leaf_2 = new Leaf;
  Component *leaf_3 = new Leaf;
  branch1->Add(leaf_1);
  branch1->Add(leaf_2);
  Component *branch2 = new Composite;
  branch2->Add(leaf_3);
  tree->Add(branch1);
  tree->Add(branch2);
  std::cout << "Client: Now I've got a composite tree:\n";
  ClientCode(tree);
  std::cout << "\n\n";

  std::cout << "Client: I don't need to check the components classes even when managing the tree:\n";
  ClientCode2(tree, simple);
  std::cout << "\n";

  delete simple;
  delete tree;
  delete branch1;
  delete branch2;
  delete leaf_1;
  delete leaf_2;
  delete leaf_3;

  return 0;
}

@endcode
<a name="9Decorator"></a><h2> 9. 装饰 Decorator </h2>

允许你通过将对象放入包含行为的特殊封装对象中来为原对象绑定新的行为。
@code
#include <iostream>
#include <string>

/**
 * Decorator
 *
 * - 模式动机: 让你把新的行为附加到对象上，把这些对象放在包含行为的特殊包装对象中。
 */
/**
 * 基础组件接口定义了可由装饰器改变的操作。
 */
class Component {
 public:
  virtual ~Component() {}
  virtual std::string Operation() const = 0;
};
/**
 * 具体组件提供操作的默认实现。这些类可能有几种变化。
 */
class ConcreteComponent : public Component {
 public:
  std::string Operation() const override {
    return "ConcreteComponent";
  }
};
/**
 * 基础装饰器类遵循与其他组件相同的接口。这个类的主要目的是为所有具体的装饰器定义包装接口。包裹代码的默认实现可能包括一个用于存储被包裹组件的字段和初始化它的方法。
 */
class Decorator : public Component {
  /**
   * @var Component
   */
 protected:
  Component* component_;

 public:
  Decorator(Component* component) : component_(component) {
  }
  /**
   * The Decorator delegates all work to the wrapped component.
   */
  std::string Operation() const override {
    return this->component_->Operation();
  }
};
/**
 * Concrete Decorators call the wrapped object and alter its result in some way.
 */
class ConcreteDecoratorA : public Decorator {
  /**
   *  装饰器可以调用操作的父级实现，而不是直接调用被包装的对象。这种方法简化了装饰器类的扩展。
   */
 public:
  ConcreteDecoratorA(Component* component) : Decorator(component) {
  }
  std::string Operation() const override {
    return "ConcreteDecoratorA(" + Decorator::Operation() + ")";
  }
};
/**
 * 装饰器可以在调用包装对象之前或之后执行其行为。
 */
class ConcreteDecoratorB : public Decorator {
 public:
  ConcreteDecoratorB(Component* component) : Decorator(component) {
  }

  std::string Operation() const override {
    return "ConcreteDecoratorB(" + Decorator::Operation() + ")";
  }
};
/**
 * 客户端代码与所有使用组件接口的对象一起工作。这样，它就可以保持独立于它所使用的组件的具体类。
 */
void ClientCode(Component* component) {
  // ...
  std::cout << "RESULT: " << component->Operation();
  // ...
}

int main() {
  /**
   * 这样一来，客户端代码既可以支持简单的组件...
   */
  Component* simple = new ConcreteComponent;
  std::cout << "Client: I've got a simple component:\n";
  ClientCode(simple);
  std::cout << "\n\n";
  /**
   * ...以及装饰好的。
   *
   * 请注意，装饰器不仅可以包裹简单的组件，还可以包裹其他装饰器。
   */
  Component* decorator1 = new ConcreteDecoratorA(simple);
  Component* decorator2 = new ConcreteDecoratorB(decorator1);
  std::cout << "Client: Now I've got a decorated component:\n";
  ClientCode(decorator2);
  std::cout << "\n";

  delete simple;
  delete decorator1;
  delete decorator2;

  return 0;
}
@endcode

<a name="10Facade"></a><h2> 10. 外观 Facade </h2>

能为程序库、框架或其他复杂类提供一个简单的接口。
@code
#include <iostream>
#include <string>

/**
 * Facade
 *
 * - 模式动机: 为一个库、一个框架或任何其他复杂的类集提供一个简化的接口。
 */

/**
 * 子系统可以直接接受来自门面或客户端的请求。在任何情况下，对于子系统来说，门面是另一个客户端，它不是子系统的一部分。
 */
class Subsystem1 {
 public:
  std::string Operation1() const {
    return "Subsystem1: Ready!\n";
  }
  // ...
  std::string OperationN() const {
    return "Subsystem1: Go!\n";
  }
};
/**
 * 有些门面可以同时与多个子系统一起工作。
 */
class Subsystem2 {
 public:
  std::string Operation1() const {
    return "Subsystem2: Get ready!\n";
  }
  // ...
  std::string OperationZ() const {
    return "Subsystem2: Fire!\n";
  }
};

/**
 * Facade类为一个或几个子系统的复杂逻辑提供一个简单的接口。Facade将客户端的请求委托给子系统中适当的对象。Facade还负责管理它们的生命周期。所有这些都使客户端免受子系统不必要的复杂性的影响。
 */
class Facade {
 protected:
  Subsystem1 *subsystem1_;
  Subsystem2 *subsystem2_;
  /**
   *  根据你的应用程序的需要，你可以向Facade提供现有的子系统对象，或者强迫Facade自己创建它们。
   */
 public:
  /**
   * 在这种情况下，我们将把内存所有权委托给Facade类。
   */
  Facade(
      Subsystem1 *subsystem1 = nullptr,
      Subsystem2 *subsystem2 = nullptr) {
    this->subsystem1_ = subsystem1 ?: new Subsystem1;
    this->subsystem2_ = subsystem2 ?: new Subsystem2;
  }
  ~Facade() {
    delete subsystem1_;
    delete subsystem2_;
  }
  /**
   *   Facade的方法是通往子系统复杂功能的方便捷径。然而，客户只能获得子系统的一小部分功能。
   */
  std::string Operation() {
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
void ClientCode(Facade *facade) {
  // ...
  std::cout << facade->Operation();
  // ...
}
/**
 * 客户端代码可能已经创建了一些子系统的对象。在这种情况下，用这些对象初始化Facade而不是让Facade创建新的实例可能是值得的。
 */

int main() {
  Subsystem1 *subsystem1 = new Subsystem1;
  Subsystem2 *subsystem2 = new Subsystem2;
  Facade *facade = new Facade(subsystem1, subsystem2);
  ClientCode(facade);

  delete facade;

  return 0;
}

@endcode
<a name="11Flyweight"></a><h2> 11. 享元 Flyweight </h2>

摒弃了在每个对象中保存所有数据的方式,通过共享多个对象所共有的相同状态,让你能在有限的内存容量中载入更多对象。
@code
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

/**
 * Flyweight
 *
 * - 模式动机: 通过在多个对象之间共享状态的共同部分，而不是在每个对象中保留所有的数据，让你在可用的RAM数量中容纳更多的对象。
 */

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

/**
 * 客户端代码通常会在应用程序的初始化阶段创建一堆预先填充的飞轮。
 */

int main()
{
    FlyweightFactory *factory = new FlyweightFactory({{"Chevrolet", "Camaro2018", "pink"}, {"Mercedes Benz", "C300", "black"}, {"Mercedes Benz", "C500", "red"}, {"BMW", "M5", "red"}, {"BMW", "X6", "white"}});
    factory->ListFlyweights();

    AddCarToPoliceDatabase(*factory,
                            "CL234IR",
                            "James Doe",
                            "BMW",
                            "M5",
                            "red");

    AddCarToPoliceDatabase(*factory,
                            "CL234IR",
                            "James Doe",
                            "BMW",
                            "X1",
                            "red");
    factory->ListFlyweights();
    delete factory;

    return 0;
}

@endcode




<a name="12Proxy"></a><h2> 12. 代理 Proxy </h2>

让你能够提供对象的替代品或其占位符。代理控制着对于原对象的访问,并允许在将请求提交给对象前后进行一些处理。
@code
#include <iostream>
/**
 * Proxy
 *
 * - 模式动机: 为另一个对象提供一个代理或占位符，以控制对原始对象的访问或增加其他责任。
 */
/**
 * Subject接口声明了RealSubject和Proxy的共同操作。只要客户端使用这个接口与RealSubject一起工作，你就能给它传递一个代理而不是一个真正的主体。
 */
class Subject {
 public:
  virtual void Request() const = 0;
};
/**
 * RealSubject包含一些核心业务逻辑。通常，RealSubjects能够做一些有用的工作，这些工作也可能非常缓慢或敏感--例如，纠正输入数据。代理可以解决这些问题，而无需对RealSubject的代码进行任何修改。
 */
class RealSubject : public Subject {
 public:
  void Request() const override {
    std::cout << "RealSubject: Handling request.\n";
  }
};
/**
 * 代理人有一个与RealSubject相同的接口。
 */
class Proxy : public Subject {
  /**
   * @var RealSubject
   */
 private:
  RealSubject *real_subject_;

  bool CheckAccess() const {
    // 一些真正的检查应该在这里进行。
    std::cout << "Proxy: Checking access prior to firing a real request.\n";
    return true;
  }
  void LogAccess() const {
    std::cout << "Proxy: Logging the time of request.\n";
  }

  /**
   * 代理人维护对RealSubject类的一个对象的引用。它可以被懒惰地加载或由客户传递给代理。
   */
 public:
  Proxy(RealSubject *real_subject) : real_subject_(new RealSubject(*real_subject)) {
  }

  ~Proxy() {
    delete real_subject_;
  }
  /**
   *  代理模式最常见的应用是懒惰加载、缓存、控制访问、记录等。一个代理可以执行这些事情中的一个，然后根据结果，将执行结果传递给链接的RealSubject对象中的同一个方法。
   */
  void Request() const override {
    if (this->CheckAccess()) {
      this->real_subject_->Request();
      this->LogAccess();
    }
  }
};
/**
 * 客户端代码应该通过主体接口与所有对象（包括主体和代理）一起工作，以便支持真实的主体和代理。然而，在现实生活中，客户端大多直接与他们的真实主体工作。在这种情况下，为了更容易实现该模式，你可以从真实主体的类中扩展你的代理。
 */
void ClientCode(const Subject &subject) {
  // ...
  subject.Request();
  // ...
}

int main() {
  std::cout << "Client: Executing the client code with a real subject:\n";
  RealSubject *real_subject = new RealSubject;
  ClientCode(*real_subject);
  std::cout << "\n";
  std::cout << "Client: Executing the same client code with a proxy:\n";
  Proxy *proxy = new Proxy(real_subject);
  ClientCode(*proxy);

  delete real_subject;
  delete proxy;
  return 0;
}

@endcode



<a name=""></a><h1> 三. 行为模式 </h1>


<a name="13ChainofResponsibility"></a><h2> 13. 责任链 Chain of Responsibility </h2>

允许你将请求沿着处理者链进行发送。收到请求后,每个处理者均可对请求进行处理,或将其传递给链上的下个处理者。

@code
#include <iostream>
#include <string>
#include <vector>

/**
 * Chain of Responsibility
 *
 * - 模式动机: 让你沿着一个处理程序链传递请求。当收到一个请求时，每个处理程序决定是处理该请求还是将其传递给链中的下一个处理程序。
 */
/**
 * 处理程序接口声明了一个用于建立处理程序链的方法。它还声明了一个执行请求的方法。
 */
class Handler {
 public:
  virtual Handler *SetNext(Handler *handler) = 0;
  virtual std::string Handle(std::string request) = 0;
};
/**
 * 默认的连锁行为可以在一个基础处理程序类里面实现。
 */
class AbstractHandler : public Handler {
  /**
   * @var Handler
   */
 private:
  Handler *next_handler_;

 public:
  AbstractHandler() : next_handler_(nullptr) {
  }
  Handler *SetNext(Handler *handler) override {
    this->next_handler_ = handler;
    // 从这里返回一个处理程序将让我们以一种方便的方式连接处理程序，就像这样。
    // $monkey->setNext($squirrel)->setNext($dog);
    return handler;
  }
  std::string Handle(std::string request) override {
    if (this->next_handler_) {
      return this->next_handler_->Handle(request);
    }

    return {};
  }
};
@endcode

@code
/**
 * 所有的具体处理程序要么处理一个请求，要么将其传递给链上的下一个处理程序。
 */
class MonkeyHandler : public AbstractHandler {
 public:
  std::string Handle(std::string request) override {
    if (request == "Banana") {
      return "Monkey: I'll eat the " + request + ".\n";
    } else {
      return AbstractHandler::Handle(request);
    }
  }
};
class SquirrelHandler : public AbstractHandler {
 public:
  std::string Handle(std::string request) override {
    if (request == "Nut") {
      return "Squirrel: I'll eat the " + request + ".\n";
    } else {
      return AbstractHandler::Handle(request);
    }
  }
};
class DogHandler : public AbstractHandler {
 public:
  std::string Handle(std::string request) override {
    if (request == "MeatBall") {
      return "Dog: I'll eat the " + request + ".\n";
    } else {
      return AbstractHandler::Handle(request);
    }
  }
};
/**
 * 客户端代码通常适合于与单个处理程序一起工作。在大多数情况下，它甚至不知道处理程序是一个链的一部分。
 */
void ClientCode(Handler &handler) {
  std::vector<std::string> food = {"Nut", "Banana", "Cup of coffee"};
  for (const std::string &f : food) {
    std::cout << "Client: Who wants a " << f << "?\n";
    const std::string result = handler.Handle(f);
    if (!result.empty()) {
      std::cout << "  " << result;
    } else {
      std::cout << "  " << f << " was left untouched.\n";
    }
  }
}
/**
 * 客户端代码的另一部分构建了实际的链。
 */
int main() {
  MonkeyHandler *monkey = new MonkeyHandler;
  SquirrelHandler *squirrel = new SquirrelHandler;
  DogHandler *dog = new DogHandler;
  monkey->SetNext(squirrel)->SetNext(dog);

  /**
   *  客户端应该能够向任何处理程序发送请求，而不仅仅是链中的第一个。
   */
  std::cout << "Chain: Monkey > Squirrel > Dog\n\n";
  ClientCode(*monkey);
  std::cout << "\n";
  std::cout << "Subchain: Squirrel > Dog\n\n";
  ClientCode(*squirrel);

  delete monkey;
  delete squirrel;
  delete dog;

  return 0;
}
@endcode




<a name="14Command"></a><h2> 14. 命令 Command </h2>

它可将请求转换为一个包含与请求相关的所有信息的独立对象。
该转换让你能根据不同的请求将方法参数化、延迟请求执行或将其放入队列中,且能实现可撤销操作。

@code
#include <iostream>
#include <string>

/**
 * Command
 *
 * - 模式动机: 将一个请求转化为一个独立的对象，其中包含关于请求的所有信息。这种转换可以让你用不同的请求对方法进行参数化，延迟或排队执行请求，并支持可撤销操作。
 */
/**
 * 命令接口声明了一个执行命令的方法。
 */
class Command {
 public:
  virtual ~Command() {
  }
  virtual void Execute() const = 0;
};
/**
 * 有些命令可以自行实现简单的操作。
 */
class SimpleCommand : public Command {
 private:
  std::string pay_load_;

 public:
  explicit SimpleCommand(std::string pay_load) : pay_load_(pay_load) {
  }
  void Execute() const override {
    std::cout << "SimpleCommand: See, I can do simple things like printing (" << this->pay_load_ << ")\n";
  }
};

/**
 * Receiver类包含一些重要的业务逻辑。它们知道如何执行与执行请求相关的各种操作。事实上，任何类都可以作为一个接收者。
 */
class Receiver {
 public:
  void DoSomething(const std::string &a) {
    std::cout << "Receiver: Working on (" << a << ".)\n";
  }
  void DoSomethingElse(const std::string &b) {
    std::cout << "Receiver: Also working on (" << b << ".)\n";
  }
};

/**
 * 然而，一些命令可以将更复杂的操作委托给其他对象，称为 "接收者"。
 */
class ComplexCommand : public Command {
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
   *   复杂的命令可以通过构造函数接受一个或几个接收器对象以及任何上下文数据。
   */
 public:
  ComplexCommand(Receiver *receiver, std::string a, std::string b) : receiver_(receiver), a_(a), b_(b) {
  }
  /**
   * 命令可以委托给接收器的任何方法。
   */
  void Execute() const override {
    std::cout << "ComplexCommand: Complex stuff should be done by a receiver object.\n";
    this->receiver_->DoSomething(this->a_);
    this->receiver_->DoSomethingElse(this->b_);
  }
};

/**
 * Invoker与一个或几个命令相关联。它向命令发送一个请求。
 */
class Invoker {
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
  ~Invoker() {
    delete on_start_;
    delete on_finish_;
  }

  void SetOnStart(Command *command) {
    this->on_start_ = command;
  }
  void SetOnFinish(Command *command) {
    this->on_finish_ = command;
  }
  /**
   *  Invoker不依赖于具体的命令或接收器类。Invoker通过执行一个命令，间接地将一个请求传递给接收器。
   */
  void DoSomethingImportant() {
    std::cout << "Invoker: Does anybody want something done before I begin?\n";
    if (this->on_start_) {
      this->on_start_->Execute();
    }
    std::cout << "Invoker: ...doing something really important...\n";
    std::cout << "Invoker: Does anybody want something done after I finish?\n";
    if (this->on_finish_) {
      this->on_finish_->Execute();
    }
  }
};
/**
 * 客户端代码可以用任何命令对调用者进行参数化。
 */

int main() {
  Invoker *invoker = new Invoker;
  invoker->SetOnStart(new SimpleCommand("Say Hi!"));
  Receiver *receiver = new Receiver;
  invoker->SetOnFinish(new ComplexCommand(receiver, "Send email", "Save report"));
  invoker->DoSomethingImportant();

  delete invoker;
  delete receiver;

  return 0;
}

@endcode
<a name="15Iterator"></a><h2> 15. 迭代器 Iterator </h2>

让你能在不暴露集合底层表现形式(列表、栈和树等)的情况下遍历集合中所有的元素。
@code
/**
 * Iterator
 *
 * - 模式动机: 让你遍历一个集合的元素而不暴露它的底层表示（列表、堆栈、树等）
 */

#include <iostream>
#include <string>
#include <vector>

/**
 * C++有自己的迭代器实现，与标准库定义的不同泛型容器一起工作。
 */

template <typename T, typename U>
class Iterator {
 public:
  typedef typename std::vector<T>::iterator iter_type;
  Iterator(U *p_data, bool reverse = false) : m_p_data_(p_data) {
    m_it_ = m_p_data_->m_data_.begin();
  }

  void First() {
    m_it_ = m_p_data_->m_data_.begin();
  }

  void Next() {
    m_it_++;
  }

  bool IsDone() {
    return (m_it_ == m_p_data_->m_data_.end());
  }

  iter_type Current() {
    return m_it_;
  }

 private:
  U *m_p_data_;
  iter_type m_it_;
};

/**
 * 通用集合/容器提供一个或几个方法来检索新鲜的迭代器实例，与集合类兼容。
 */

template <class T>
class Container {
  friend class Iterator<T, Container>;

 public:
  void Add(T a) {
    m_data_.push_back(a);
  }

  Iterator<T, Container> *CreateIterator() {
    return new Iterator<T, Container>(this);
  }

 private:
  std::vector<T> m_data_;
};

class Data {
 public:
  Data(int a = 0) : m_data_(a) {}

  void set_data(int a) {
    m_data_ = a;
  }

  int data() {
    return m_data_;
  }

 private:
  int m_data_;
};

/**
 * 客户端代码可能知道也可能不知道具体的迭代器或集合类，对于这个实现来说，容器是通用的，所以你可以用一个int或一个自定义类来使用。
 */
void ClientCode() {
  std::cout << "________________Iterator with int______________________________________" << std::endl;
  Container<int> cont;

  for (int i = 0; i < 10; i++) {
    cont.Add(i);
  }

  Iterator<int, Container<int>> *it = cont.CreateIterator();
  for (it->First(); !it->IsDone(); it->Next()) {
    std::cout << *it->Current() << std::endl;
  }

  Container<Data> cont2;
  Data a(100), b(1000), c(10000);
  cont2.Add(a);
  cont2.Add(b);
  cont2.Add(c);

  std::cout << "________________Iterator with custom Class______________________________" << std::endl;
  Iterator<Data, Container<Data>> *it2 = cont2.CreateIterator();
  for (it2->First(); !it2->IsDone(); it2->Next()) {
    std::cout << it2->Current()->data() << std::endl;
  }
  delete it;
  delete it2;
}

int main() {
  ClientCode();
  return 0;
}

@endcode
<a name="16Mediator"></a><h2> 16. 中介者 Mediator </h2>

能让你减少对象之间混乱无序的依赖关系。该模式会限制对象之间的直接交互,迫使它们通过一个中介者对象进行合作。
@code

#include <iostream>
#include <string>
/**
 * Mediator
 *
 * - 模式动机: 让你减少对象之间混乱的依赖关系。该模式限制了对象之间的直接通信，迫使它们只通过一个中介对象进行合作。
 */

/**
 * 调解器接口声明了一个被组件用来通知调解器各种事件的方法。调解器可以对这些事件做出反应，并将执行结果传递给其他组件。
 */
class BaseComponent;
class Mediator {
 public:
  virtual void Notify(BaseComponent *sender, std::string event) const = 0;
};

/**
 * 基础组件提供了在组件对象内存储调解人实例的基本功能。
 */
class BaseComponent {
 protected:
  Mediator *mediator_;

 public:
  BaseComponent(Mediator *mediator = nullptr) : mediator_(mediator) {
  }
  void set_mediator(Mediator *mediator) {
    this->mediator_ = mediator;
  }
};

/**
 * 具体的组件实现各种功能。它们不依赖其他组件。它们也不依赖于任何具体的调解器类。
 */
class Component1 : public BaseComponent {
 public:
  void DoA() {
    std::cout << "Component 1 does A.\n";
    this->mediator_->Notify(this, "A");
  }
  void DoB() {
    std::cout << "Component 1 does B.\n";
    this->mediator_->Notify(this, "B");
  }
};

class Component2 : public BaseComponent {
 public:
  void DoC() {
    std::cout << "Component 2 does C.\n";
    this->mediator_->Notify(this, "C");
  }
  void DoD() {
    std::cout << "Component 2 does D.\n";
    this->mediator_->Notify(this, "D");
  }
};

/**
 * 具体的调解人通过协调几个组件来实现合作行为。
 */
class ConcreteMediator : public Mediator {
 private:
  Component1 *component1_;
  Component2 *component2_;

 public:
  ConcreteMediator(Component1 *c1, Component2 *c2) : component1_(c1), component2_(c2) {
    this->component1_->set_mediator(this);
    this->component2_->set_mediator(this);
  }
  void Notify(BaseComponent *sender, std::string event) const override {
    if (event == "A") {
      std::cout << "Mediator reacts on A and triggers following operations:\n";
      this->component2_->DoC();
    }
    if (event == "D") {
      std::cout << "Mediator reacts on D and triggers following operations:\n";
      this->component1_->DoB();
      this->component2_->DoC();
    }
  }
};

/**
 * The client code.
 */

void ClientCode() {
  Component1 *c1 = new Component1;
  Component2 *c2 = new Component2;
  ConcreteMediator *mediator = new ConcreteMediator(c1, c2);
  std::cout << "Client triggers operation A.\n";
  c1->DoA();
  std::cout << "\n";
  std::cout << "Client triggers operation D.\n";
  c2->DoD();

  delete c1;
  delete c2;
  delete mediator;
}

int main() {
  ClientCode();
  return 0;
}

@endcode
<a name="17Memento"></a><h2> 17. 备忘录 Memento </h2>

允许在不暴露对象实现细节的情况下保存和恢复对象之前的状态。
@code
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <string>
#include <vector>

/**
 * Memento
 *
 * - 模式动机: 让你保存和恢复一个对象以前的状态，而不透露其实现的细节。
 */

/**
 * Memento接口提供了一种检索Memento元数据的方法，例如创建日期或名称。然而，它并没有暴露发起者的状态。
 */
class Memento {
 public:
  virtual std::string GetName() const = 0;
  virtual std::string date() const = 0;
  virtual std::string state() const = 0;
};

/**
 * Concrete Memento包含用于存储发起人状态的基础设施。
 */
class ConcreteMemento : public Memento {
 private:
  std::string state_;
  std::string date_;

 public:
  ConcreteMemento(std::string state) : state_(state) {
    this->state_ = state;
    std::time_t now = std::time(0);
    this->date_ = std::ctime(&now);
  }
  /**
   * 发起人在恢复其状态时使用这种方法。
   */
  std::string state() const override {
    return this->state_;
  }
  /**
   * 其余的方法是由看守所用来显示元数据。
   */
  std::string GetName() const override {
    return this->date_ + " / (" + this->state_.substr(0, 9) + "...)";
  }
  std::string date() const override {
    return this->date_;
  }
};

/**
 * 发起人持有一些可能随时间变化的重要状态。它还定义了一种将状态保存在备忘录内的方法和另一种将状态恢复的方法。
 */
class Originator {
  /**
   * @var string 为了简单起见，发起人的状态被储存在一个单一的变量内。
   */
 private:
  std::string state_;

  std::string GenerateRandomString(int length = 10) {
    const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    int stringLength = sizeof(alphanum) - 1;

    std::string random_string;
    for (int i = 0; i < length; i++) {
      random_string += alphanum[std::rand() % stringLength];
    }
    return random_string;
  }

 public:
  Originator(std::string state) : state_(state) {
    std::cout << "Originator: My initial state is: " << this->state_ << "\n";
  }
  /**
   *   发起人的业务逻辑可能会影响其内部状态。因此，在通过save()方法启动业务逻辑的方法之前，客户应该备份状态。
   */
  void DoSomething() {
    std::cout << "Originator: I'm doing something important.\n";
    this->state_ = this->GenerateRandomString(30);
    std::cout << "Originator: and my state has changed to: " << this->state_ << "\n";
  }

  /**
   * 将当前状态保存在一个Memento内。
   */
  Memento *Save() {
    return new ConcreteMemento(this->state_);
  }
  /**
   * 从一个备忘录对象中恢复发起人的状态。
   */
  void Restore(Memento *memento) {
    this->state_ = memento->state();
    std::cout << "Originator: My state has changed to: " << this->state_ << "\n";
  }
};

/**
 * Caretaker并不依赖于具体的Memento类。因此，它不能访问存储在Memento中的发起人的状态。它通过基础的Memento接口与所有的纪念品一起工作。
 */
class Caretaker {
  /**
   * @var Memento[]
   */
 private:
  std::vector<Memento *> mementos_;

  /**
   * @var Originator
   */
  Originator *originator_;

 public:
  Caretaker(Originator *originator) : originator_(originator) {
    this->originator_ = originator;
  }

  void Backup() {
    std::cout << "\nCaretaker: Saving Originator's state...\n";
    this->mementos_.push_back(this->originator_->Save());
  }
  void Undo() {
    if (!this->mementos_.size()) {
      return;
    }
    Memento *memento = this->mementos_.back();
    this->mementos_.pop_back();
    std::cout << "Caretaker: Restoring state to: " << memento->GetName() << "\n";
    try {
      this->originator_->Restore(memento);
    } catch (...) {
      this->Undo();
    }
  }
  void ShowHistory() const {
    std::cout << "Caretaker: Here's the list of mementos:\n";
    for (Memento *memento : this->mementos_) {
      std::cout << memento->GetName() << "\n";
    }
  }
};
/**
 * Client code.
 */

void ClientCode() {
  Originator *originator = new Originator("Super-duper-super-puper-super.");
  Caretaker *caretaker = new Caretaker(originator);
  caretaker->Backup();
  originator->DoSomething();
  caretaker->Backup();
  originator->DoSomething();
  caretaker->Backup();
  originator->DoSomething();
  std::cout << "\n";
  caretaker->ShowHistory();
  std::cout << "\nClient: Now, let's rollback!\n\n";
  caretaker->Undo();
  std::cout << "\nClient: Once more!\n\n";
  caretaker->Undo();

  delete originator;
  delete caretaker;
}

int main() {
  std::srand(static_cast<unsigned int>(std::time(NULL)));
  ClientCode();
  return 0;
}
@endcode
<a name="18Observer"></a><h2> 18. 观察者 Observer </h2>

允许你定义一种订阅机制,可在对象事件发生时通知多个“观察”该对象的其他对象。
@code
/**
 * Observer
 *
 * - 模式动机: 让你定义一个订阅机制来通知多个对象关于他们正在观察的对象发生的任何事件。

 * 请注意，有很多不同的术语与这种模式有类似的含义。只要记住，主体也被称为发布者，而观察者通常被称为订阅者，反之亦然。另外，动词 "观察"、"倾听 "或 "跟踪 "通常也是同一个意思。
 */

#include <iostream>
#include <list>
#include <string>

class IObserver {
 public:
  virtual ~IObserver(){};
  virtual void Update(const std::string &message_from_subject) = 0;
};

class ISubject {
 public:
  virtual ~ISubject(){};
  virtual void Attach(IObserver *observer) = 0;
  virtual void Detach(IObserver *observer) = 0;
  virtual void Notify() = 0;
};

/**
 * 主体拥有一些重要的状态，并在状态改变时通知观察者。
 */

class Subject : public ISubject {
 public:
  virtual ~Subject() {
    std::cout << "Goodbye, I was the Subject.\n";
  }

  /**
   * 订阅管理方法。
   */
  void Attach(IObserver *observer) override {
    list_observer_.push_back(observer);
  }
  void Detach(IObserver *observer) override {
    list_observer_.remove(observer);
  }
  void Notify() override {
    std::list<IObserver *>::iterator iterator = list_observer_.begin();
    HowManyObserver();
    while (iterator != list_observer_.end()) {
      (*iterator)->Update(message_);
      ++iterator;
    }
  }

  void CreateMessage(std::string message = "Empty") {
    this->message_ = message;
    Notify();
  }
  void HowManyObserver() {
    std::cout << "There are " << list_observer_.size() << " observers in the list.\n";
  }

  /**
   * 通常情况下，订阅逻辑只是Subject真正能做的一小部分。主题通常持有一些重要的业务逻辑，每当有重要的事情即将发生（或发生后），就会触发一个通知方法。
   */
  void SomeBusinessLogic() {
    this->message_ = "change message message";
    Notify();
    std::cout << "I'm about to do some thing important\n";
  }

 private:
  std::list<IObserver *> list_observer_;
  std::string message_;
};

class Observer : public IObserver {
 public:
  Observer(Subject &subject) : subject_(subject) {
    this->subject_.Attach(this);
    std::cout << "Hi, I'm the Observer \"" << ++Observer::static_number_ << "\".\n";
    this->number_ = Observer::static_number_;
  }
  virtual ~Observer() {
    std::cout << "Goodbye, I was the Observer \"" << this->number_ << "\".\n";
  }

  void Update(const std::string &message_from_subject) override {
    message_from_subject_ = message_from_subject;
    PrintInfo();
  }
  void RemoveMeFromTheList() {
    subject_.Detach(this);
    std::cout << "Observer \"" << number_ << "\" removed from the list.\n";
  }
  void PrintInfo() {
    std::cout << "Observer \"" << this->number_ << "\": a new message is available --> " << this->message_from_subject_ << "\n";
  }

 private:
  std::string message_from_subject_;
  Subject &subject_;
  static int static_number_;
  int number_;
};

int Observer::static_number_ = 0;

void ClientCode() {
  Subject *subject = new Subject;
  Observer *observer1 = new Observer(*subject);
  Observer *observer2 = new Observer(*subject);
  Observer *observer3 = new Observer(*subject);
  Observer *observer4;
  Observer *observer5;

  subject->CreateMessage("Hello World! :D");
  observer3->RemoveMeFromTheList();

  subject->CreateMessage("The weather is hot today! :p");
  observer4 = new Observer(*subject);

  observer2->RemoveMeFromTheList();
  observer5 = new Observer(*subject);

  subject->CreateMessage("My new car is great! ;)");
  observer5->RemoveMeFromTheList();

  observer4->RemoveMeFromTheList();
  observer1->RemoveMeFromTheList();

  delete observer5;
  delete observer4;
  delete observer3;
  delete observer2;
  delete observer1;
  delete subject;
}

int main() {
  ClientCode();
  return 0;
}

@endcode
<a name="19State"></a><h2> 19. 状态 State </h2>

让你能在一个对象的内部状态变化时改变其行为,使其看上去就像改变了自身所属的类一样。
@code
#include <iostream>
#include <typeinfo>
/**
 * State
 *
 *  模式动机: 当一个对象的内部状态改变时，让它改变其行为。它看起来就像该对象改变了它的类别。
 */

/**
 * 基状态类声明了所有具体状态应该实现的方法，并且还提供了一个与状态相关的上下文对象的反向引用。这个反向引用可以被状态用来将上下文过渡到另一个状态。
 */

class Context;

class State {
  /**
   * @var Context
   */
 protected:
  Context *context_;

 public:
  virtual ~State() {
  }

  void set_context(Context *context) {
    this->context_ = context;
  }

  virtual void Handle1() = 0;
  virtual void Handle2() = 0;
};

/**
 * Context定义了客户感兴趣的接口。它还维护一个对State子类实例的引用，该实例代表Context的当前状态。
 */
class Context {
  /**
   * @var State 对Context当前状态的引用。
   */
 private:
  State *state_;

 public:
  Context(State *state) : state_(nullptr) {
    this->TransitionTo(state);
  }
  ~Context() {
    delete state_;
  }
  /**
   * Context允许在运行时改变State对象。
   */
  void TransitionTo(State *state) {
    std::cout << "Context: Transition to " << typeid(*state).name() << ".\n";
    if (this->state_ != nullptr)
      delete this->state_;
    this->state_ = state;
    this->state_->set_context(this);
  }
  /**
   * Context将其部分行为委托给当前的State对象。
   */
  void Request1() {
    this->state_->Handle1();
  }
  void Request2() {
    this->state_->Handle2();
  }
};

/**
 * 具体状态实现各种行为，与Context的状态相关。
 */

class ConcreteStateA : public State {
 public:
  void Handle1() override;

  void Handle2() override {
    std::cout << "ConcreteStateA handles request2.\n";
  }
};

class ConcreteStateB : public State {
 public:
  void Handle1() override {
    std::cout << "ConcreteStateB handles request1.\n";
  }
  void Handle2() override {
    std::cout << "ConcreteStateB handles request2.\n";
    std::cout << "ConcreteStateB wants to change the state of the context.\n";
    this->context_->TransitionTo(new ConcreteStateA);
  }
};

void ConcreteStateA::Handle1() {
  {
    std::cout << "ConcreteStateA handles request1.\n";
    std::cout << "ConcreteStateA wants to change the state of the context.\n";

    this->context_->TransitionTo(new ConcreteStateB);
  }
}

/**
 * The client code.
 */
void ClientCode() {
  Context *context = new Context(new ConcreteStateA);
  context->Request1();
  context->Request2();
  delete context;
}

int main() {
  ClientCode();
  return 0;
}

@endcode
<a name="20Strategy"></a><h2> 20. 策略 Strategy </h2>

能让你定义一系列算法,并将每种算法分别放入独立的类中, 以使算法的对象能够相互替换。
@code
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

/**
 * Strategy
 *
 * - 模式动机: 让你定义一个算法系列，把它们每个都放到一个单独的类中，并使它们的对象可以互换。
 */

/**
 * 策略接口声明了一些算法的所有支持版本所共有的操作。
 * Context使用这个接口来调用Concrete Strategies定义的算法。
 */
class Strategy
{
public:
    virtual ~Strategy() {}
    virtual std::string DoAlgorithm(const std::vector<std::string> &data) const = 0;
};

/**
 * Context定义了客户感兴趣的接口。
 */

class Context
{
    /**
     * @var Strategy Context维护对策略对象之一的引用。Context不知道一个策略的具体类别。它应该通过策略接口与所有策略一起工作。
     */
private:
    Strategy *strategy_;
    /**
     * 通常，Context通过构造函数接受一个策略，但也提供一个setter来在运行时改变它。
     */
public:
    Context(Strategy *strategy = nullptr) : strategy_(strategy)
    {
    }
    ~Context()
    {
        delete this->strategy_;
    }
    /**
     * 通常，Context允许在运行时替换一个策略对象。
     */
    void set_strategy(Strategy *strategy)
    {
        delete this->strategy_;
        this->strategy_ = strategy;
    }
    /**
     *    Context将一些工作委托给Strategy对象，而不是自己实现+多个版本的算法。
     */
    void DoSomeBusinessLogic() const
    {
        // ...
        std::cout << "Context: Sorting data using the strategy (not sure how it'll do it)\n";
        std::string result = this->strategy_->DoAlgorithm(std::vector<std::string>{"a", "e", "c", "b", "d"});
        std::cout << result << "\n";
        // ...
    }
};

/**
 * 具体的策略在遵循基本策略接口的同时实现了该算法。该接口使它们在上下文中可以互换。
 */
class ConcreteStrategyA : public Strategy
{
public:
    std::string DoAlgorithm(const std::vector<std::string> &data) const override
    {
        std::string result;
        std::for_each(std::begin(data), std::end(data), [&result](const std::string &letter) {
            result += letter;
        });
        std::sort(std::begin(result), std::end(result));

        return result;
    }
};
class ConcreteStrategyB : public Strategy
{
    std::string DoAlgorithm(const std::vector<std::string> &data) const override
    {
        std::string result;
        std::for_each(std::begin(data), std::end(data), [&result](const std::string &letter) {
            result += letter;
        });
        std::sort(std::begin(result), std::end(result));
        for (int i = 0; i < result.size() / 2; i++)
        {
            std::swap(result[i], result[result.size() - i - 1]);
        }

        return result;
    }
};
/**
 * 客户端代码挑选一个具体的策略并将其传递给上下文。客户端应该意识到策略之间的差异，以便做出正确的选择。
 */

void ClientCode()
{
    Context *context = new Context(new ConcreteStrategyA);
    std::cout << "Client: Strategy is set to normal sorting.\n";
    context->DoSomeBusinessLogic();
    std::cout << "\n";
    std::cout << "Client: Strategy is set to reverse sorting.\n";
    context->set_strategy(new ConcreteStrategyB);
    context->DoSomeBusinessLogic();
    delete context;
}

int main()
{
    ClientCode();
    return 0;
}
@endcode
<a name="21TemplateMethod"></a><h2> 21. 模板方法 Template Method </h2>

在超类中定义一个算法的框架,允许子类在不修改结构的情况下重写算法的特定步骤。
@code
#include <iostream>

/**
 * Template Method
 *
 * - 模式动机: 在超类中定义了算法的骨架，但允许子类在不改变算法结构的情况下覆盖算法的特定步骤。
 */
/**
 * 抽象类定义了一个模板方法，它包含了某种算法的骨架，由对（通常）抽象的原始操作的调用组成。
 * 具体的子类应该实现这些操作，但要保持模板方法本身的完整性。
 */
class AbstractClass {
  /**
   * 模板方法定义了一个算法的骨架。
   */
 public:
  void TemplateMethod() const {
    this->BaseOperation1();
    this->RequiredOperations1();
    this->BaseOperation2();
    this->Hook1();
    this->RequiredOperation2();
    this->BaseOperation3();
    this->Hook2();
  }
  /**
   * 这些操作已经有了实现。
   */
 protected:
  void BaseOperation1() const {
    std::cout << "AbstractClass says: I am doing the bulk of the work\n";
  }
  void BaseOperation2() const {
    std::cout << "AbstractClass says: But I let subclasses override some operations\n";
  }
  void BaseOperation3() const {
    std::cout << "AbstractClass says: But I am doing the bulk of the work anyway\n";
  }
  /**
   * 这些操作必须在子类中实现。
   */
  virtual void RequiredOperations1() const = 0;
  virtual void RequiredOperation2() const = 0;
  /**
   *  这些是 "钩子"。子类可以覆盖它们，但这不是强制性的，因为这些钩子已经有默认的（但空的）实现。钩子在算法的一些关键地方提供了额外的扩展点。
   */
  virtual void Hook1() const {}
  virtual void Hook2() const {}
};
/**
 * 具体类必须实现基类的所有抽象操作。他们也可以用一个默认的实现来覆盖一些操作。
 */
class ConcreteClass1 : public AbstractClass {
 protected:
  void RequiredOperations1() const override {
    std::cout << "ConcreteClass1 says: Implemented Operation1\n";
  }
  void RequiredOperation2() const override {
    std::cout << "ConcreteClass1 says: Implemented Operation2\n";
  }
};
/**
 * 通常情况下，具体类只覆盖基类的一部分操作。
 */
class ConcreteClass2 : public AbstractClass {
 protected:
  void RequiredOperations1() const override {
    std::cout << "ConcreteClass2 says: Implemented Operation1\n";
  }
  void RequiredOperation2() const override {
    std::cout << "ConcreteClass2 says: Implemented Operation2\n";
  }
  void Hook1() const override {
    std::cout << "ConcreteClass2 says: Overridden Hook1\n";
  }
};
/**
 * 客户端代码调用模板方法来执行该算法。客户端代码不必知道它所处理的对象的具体类别，只要它通过其基类的接口与对象一起工作即可。
 */
void ClientCode(AbstractClass *class_) {
  // ...
  class_->TemplateMethod();
  // ...
}

int main() {
  std::cout << "Same client code can work with different subclasses:\n";
  ConcreteClass1 *concreteClass1 = new ConcreteClass1;
  ClientCode(concreteClass1);
  std::cout << "\n";
  std::cout << "Same client code can work with different subclasses:\n";
  ConcreteClass2 *concreteClass2 = new ConcreteClass2;
  ClientCode(concreteClass2);
  delete concreteClass1;
  delete concreteClass2;
  return 0;
}

@endcode
<a name="22Visitor"></a><h2> 22. 访问者 Visitor </h2>

将算法与其所作用的对象隔离开来。
@code
#include <array>
#include <iostream>
#include <string>

/**
 * Visitor
 *
 * - 模式动机: 让你把算法与它们所操作的对象分开。
 */

/**
 * 访客接口声明了一组与组件类对应的访问方法。一个访问方法的签名允许访问者识别它所处理的组件的确切类别。
 */
class ConcreteComponentA;
class ConcreteComponentB;

class Visitor {
 public:
  virtual void VisitConcreteComponentA(const ConcreteComponentA *element) const = 0;
  virtual void VisitConcreteComponentB(const ConcreteComponentB *element) const = 0;
};

/**
 * 组件接口声明了一个`accept'方法，该方法应将基础访问者接口作为参数。
 */

class Component {
 public:
  virtual ~Component() {}
  virtual void Accept(Visitor *visitor) const = 0;
};

/**
 * 每个具体的组件都必须以这样的方式实现`接受'方法，即调用与组件类相对应的访问者的方法。
 */
class ConcreteComponentA : public Component {
  /**
   *  注意，我们正在调用`visitConcreteComponentA`，它与当前的类名相匹配。这样，我们让访问者知道它所工作的组件的类。
   */
 public:
  void Accept(Visitor *visitor) const override {
    visitor->VisitConcreteComponentA(this);
  }
  /**
   * 具体组件可能有一些特殊的方法，这些方法在其基类或接口中并不存在。访问者仍然能够使用这些方法，因为它知道该组件的具体类。
   */
  std::string ExclusiveMethodOfConcreteComponentA() const {
    return "A";
  }
};

class ConcreteComponentB : public Component {
  /**
   * Same here: visitConcreteComponentB => ConcreteComponentB
   */
 public:
  void Accept(Visitor *visitor) const override {
    visitor->VisitConcreteComponentB(this);
  }
  std::string SpecialMethodOfConcreteComponentB() const {
    return "B";
  }
};

/**
 * 具体的Visitor实现了同一算法的多个版本，它可以与所有的具体组件类一起工作。

 * 在与复杂的对象结构（如复合树）一起使用时，你可以体验到访问者模式的最大好处。在这种情况下，当在结构的各个对象上执行访问者的方法时，存储一些算法的中间状态可能会有帮助。
 */
class ConcreteVisitor1 : public Visitor {
 public:
  void VisitConcreteComponentA(const ConcreteComponentA *element) const override {
    std::cout << element->ExclusiveMethodOfConcreteComponentA() << " + ConcreteVisitor1\n";
  }

  void VisitConcreteComponentB(const ConcreteComponentB *element) const override {
    std::cout << element->SpecialMethodOfConcreteComponentB() << " + ConcreteVisitor1\n";
  }
};

class ConcreteVisitor2 : public Visitor {
 public:
  void VisitConcreteComponentA(const ConcreteComponentA *element) const override {
    std::cout << element->ExclusiveMethodOfConcreteComponentA() << " + ConcreteVisitor2\n";
  }
  void VisitConcreteComponentB(const ConcreteComponentB *element) const override {
    std::cout << element->SpecialMethodOfConcreteComponentB() << " + ConcreteVisitor2\n";
  }
};
/**
 * 客户端代码可以在任何元素的集合上运行访问者操作，而不需要弄清楚它们的具体类别。接受操作会引导对访问者对象中适当的操作的调用。
 */
void ClientCode(std::array<const Component *, 2> components, Visitor *visitor) {
  // ...
  for (const Component *comp : components) {
    comp->Accept(visitor);
  }
  // ...
}

int main() {
  std::array<const Component *, 2> components = {new ConcreteComponentA, new ConcreteComponentB};
  std::cout << "The client code works with all visitors via the base Visitor interface:\n";
  ConcreteVisitor1 *visitor1 = new ConcreteVisitor1;
  ClientCode(components, visitor1);
  std::cout << "\n";
  std::cout << "It allows the same client code to work with different types of visitors:\n";
  ConcreteVisitor2 *visitor2 = new ConcreteVisitor2;
  ClientCode(components, visitor2);

  for (const Component *comp : components) {
    delete comp;
  }
  delete visitor1;
  delete visitor2;

  return 0;
}

@endcode

 *
 *
 * <a name="CommProg"></a>
 * <h1> The commented program</h1>
 * 
 * 
 * @code
 * #include <iostream>
 * #include <string>
 * 
 * /**
 *  * Abstract Factory
 *  *
 *  * Intent: Lets you produce families of related objects without specifying their
 *  * concrete classes.
 *  */
 * 
 * /**
 *  * Each distinct product of a product family should have a base interface. All
 *  * variants of the product must implement this interface.
 *  */
 * class AbstractProductA
 * {
 * public:
 *   virtual ~AbstractProductA(){};
 *   virtual std::string UsefulFunctionA() const = 0;
 * };
 * 
 * /**
 *  * Concrete Products are created by corresponding Concrete Factories.
 *  */
 * class ConcreteProductA1 : public AbstractProductA
 * {
 * public:
 *   std::string UsefulFunctionA() const override
 *   {
 *     return "The result of the product A1.";
 *   }
 * };
 * 
 * class ConcreteProductA2 : public AbstractProductA
 * {
 *   std::string UsefulFunctionA() const override
 *   {
 *     return "The result of the product A2.";
 *   }
 * };
 * 
 * /**
 *  * Here's the the base interface of another product. All products can interact
 *  * with each other, but proper interaction is possible only between products of
 *  * the same concrete variant.
 *  */
 * class AbstractProductB
 * {
 *   /**
 *    * Product B is able to do its own thing...
 *    */
 * public:
 *   virtual ~AbstractProductB(){};
 *   virtual std::string UsefulFunctionB() const = 0;
 *   /**
 *    * ...but it also can collaborate with the ProductA.
 *    *
 *    * The Abstract Factory makes sure that all products it creates are of the
 *    * same variant and thus, compatible.
 *    */
 *   virtual std::string
 *   AnotherUsefulFunctionB(const AbstractProductA &collaborator) const = 0;
 * };
 * 
 * /**
 *  * Concrete Products are created by corresponding Concrete Factories.
 *  */
 * class ConcreteProductB1 : public AbstractProductB
 * {
 * public:
 *   std::string UsefulFunctionB() const override
 *   {
 *     return "The result of the product B1.";
 *   }
 *   /**
 *    * The variant, Product B1, is only able to work correctly with the variant,
 *    * Product A1. Nevertheless, it accepts any instance of AbstractProductA as an
 *    * argument.
 *    */
 *   std::string
 *   AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override
 *   {
 *     const std::string result = collaborator.UsefulFunctionA();
 *     return "The result of the B1 collaborating with ( " + result + " )";
 *   }
 * };
 * 
 * class ConcreteProductB2 : public AbstractProductB
 * {
 * public:
 *   std::string UsefulFunctionB() const override
 *   {
 *     return "The result of the product B2.";
 *   }
 *   /**
 *    * The variant, Product B2, is only able to work correctly with the variant,
 *    * Product A2. Nevertheless, it accepts any instance of AbstractProductA as an
 *    * argument.
 *    */
 *   std::string
 *   AnotherUsefulFunctionB(const AbstractProductA &collaborator) const override
 *   {
 *     const std::string result = collaborator.UsefulFunctionA();
 *     return "The result of the B2 collaborating with ( " + result + " )";
 *   }
 * };
 * 
 * /**
 *  * The Abstract Factory interface declares a set of methods that return
 *  * different abstract products. These products are called a family and are
 *  * related by a high-level theme or concept. Products of one family are usually
 *  * able to collaborate among themselves. A family of products may have several
 *  * variants, but the products of one variant are incompatible with products of
 *  * another.
 *  */
 * class AbstractFactory
 * {
 * public:
 *   virtual AbstractProductA *CreateProductA() const = 0;
 *   virtual AbstractProductB *CreateProductB() const = 0;
 * };
 * 
 * /**
 *  * Concrete Factories produce a family of products that belong to a single
 *  * variant. The factory guarantees that resulting products are compatible. Note
 *  * that signatures of the Concrete Factory's methods return an abstract product,
 *  * while inside the method a concrete product is instantiated.
 *  */
 * class ConcreteFactory1 : public AbstractFactory
 * {
 * public:
 *   AbstractProductA *CreateProductA() const override
 *   {
 *     return new ConcreteProductA1();
 *   }
 *   AbstractProductB *CreateProductB() const override
 *   {
 *     return new ConcreteProductB1();
 *   }
 * };
 * 
 * /**
 *  * Each Concrete Factory has a corresponding product variant.
 *  */
 * class ConcreteFactory2 : public AbstractFactory
 * {
 * public:
 *   AbstractProductA *CreateProductA() const override
 *   {
 *     return new ConcreteProductA2();
 *   }
 *   AbstractProductB *CreateProductB() const override
 *   {
 *     return new ConcreteProductB2();
 *   }
 * };
 * 
 * /**
 *  * The client code works with factories and products only through abstract
 *  * types: AbstractFactory and AbstractProduct. This lets you pass any factory or
 *  * product subclass to the client code without breaking it.
 *  */
 * 
 * void ClientCode(const AbstractFactory &factory)
 * {
 *   const AbstractProductA *product_a = factory.CreateProductA();
 *   const AbstractProductB *product_b = factory.CreateProductB();
 *   std::cout << product_b->UsefulFunctionB() << "\n";
 *   std::cout << product_b->AnotherUsefulFunctionB(*product_a) << "\n";
 *   delete product_a;
 *   delete product_b;
 * }
 * 
 * int main()
 * {
 *   std::cout << "Client: Testing client code with the first factory type:\n";
 *   ConcreteFactory1 *f1 = new ConcreteFactory1();
 *   ClientCode(*f1);
 *   delete f1;
 *   std::cout << std::endl;
 *   std::cout
 *     << "Client: Testing the same client code with the second factory type:\n";
 *   ConcreteFactory2 *f2 = new ConcreteFactory2();
 *   ClientCode(*f2);
 *   delete f2;
 *   return 0;
 * } * @endcode
<a name="Results"></a><h1>Results</h1>

用于测试代码的结果。 *
 *
<a name="PlainProg"></a>
<h1> The plain program</h1>
@include "step-0a.cc"
*/
