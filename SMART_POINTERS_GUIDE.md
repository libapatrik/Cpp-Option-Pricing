# Smart Pointers Guide: Learning from CppFM Codebase

## Table of Contents

1. [Introduction: Why Smart Pointers?](#introduction)
2. [std::unique_ptr](#stdunique_ptr)
3. [std::make_unique](#stdmake_unique)
4. [std::shared_ptr](#stdshared_ptr)
5. [std::make_shared](#stdmake_shared)
6. [std::weak_ptr](#stdweak_ptr)
7. [When to Use Which?](#when-to-use-which)
8. [Common Patterns in Your Codebase](#common-patterns)

---

## Introduction: Why Smart Pointers?

### The Problem with Raw Pointers

**Old C++ way (error-prone):**

```cpp
// Manual memory management - DANGEROUS!
ExtrapolationScheme* scheme = new FlatExtrapolation();
// ... use scheme ...
delete scheme;  // Must remember to delete!
// What if an exception occurs before delete? MEMORY LEAK!
```

**Problems:**

- Memory leaks if you forget `delete`
- Double deletion if you delete twice
- Dangling pointers if you delete and then use
- Exception safety issues

### The Solution: Smart Pointers

Smart pointers automatically manage memory using **RAII** (Resource Acquisition Is Initialization):

- **Acquire** memory when created
- **Release** memory when destroyed (automatically)
- **Exception-safe**: Even if exceptions occur, memory is cleaned up

---

## std::unique_ptr

### What is `std::unique_ptr`?

`std::unique_ptr` is a **smart pointer that owns and manages a single object** through a pointer. It ensures:

- **Exclusive ownership**: Only one `unique_ptr` can own an object
- **Automatic deletion**: Object is deleted when `unique_ptr` goes out of scope
- **No copying**: Cannot be copied (prevents double deletion)
- **Movable**: Can be moved (transfers ownership)

### Examples from Your Codebase

#### Example 1: Member Variable Storage

**Location:** `InterpolationSchemes.h` line 54

```cpp
class InterpolationScheme {
protected:
    std::unique_ptr<ExtrapolationScheme> _extrapolationScheme;
};
```

**What's happening:**

- `_extrapolationScheme` owns an `ExtrapolationScheme` object
- When `InterpolationScheme` is destroyed, `_extrapolationScheme` is automatically destroyed
- The owned `ExtrapolationScheme` is automatically deleted
- **No manual `delete` needed!**

#### Example 2: Creating Objects with `make_unique`

**Location:** `InterpolationSchemes.cpp` line 64

```cpp
LinearInterpolation::LinearInterpolation(const std::vector<double>& xData, 
                                         const std::vector<double>& yData)
{
    // Initialize with default QuadraticExtrapolation strategy
    _extrapolationScheme = std::make_unique<QuadraticExtrapolation>();
    _extrapolationScheme->initialize(*this);
}
```

**What's happening:**

- `std::make_unique<QuadraticExtrapolation>()` creates a new `QuadraticExtrapolation` object on the heap
- Returns a `std::unique_ptr<QuadraticExtrapolation>`
- Ownership is transferred to `_extrapolationScheme`
- Object is automatically deleted when `LinearInterpolation` is destroyed

#### Example 3: Transferring Ownership with `std::move`

**Location:** `InterpolationSchemes.cpp` line 12-19

```cpp
void InterpolationScheme::setExtrapolationScheme(std::unique_ptr<ExtrapolationScheme> scheme)
{
    if (!scheme) {
        throw std::invalid_argument("Extrapolation scheme cannot be null");
    }
    _extrapolationScheme = std::move(scheme);   // Transfer ownership
    _extrapolationScheme->initialize(*this);
}
```

**What's happening:**

- Function parameter takes `std::unique_ptr` by value (ownership transfer)
- `std::move(scheme)` transfers ownership from caller to `_extrapolationScheme`
- After `std::move`, `scheme` is empty (nullptr)
- **Cannot copy `unique_ptr`** - must use `std::move` to transfer

**Usage example:**

```cpp
// From main.cpp line 59-72
std::unique_ptr<ExtrapolationScheme> scheme = std::make_unique<FlatExtrapolation>();
auto interp = std::make_unique<LinearInterpolation>(xData, yData, std::move(scheme));
// After std::move(scheme), 'scheme' is now empty (nullptr)
```

#### Example 4: Returning `unique_ptr` from Functions

**Location:** `InterpolationSchemes.cpp` line 118-125

```cpp
std::unique_ptr<InterpolationScheme> LinearInterpolation::clone() const
{
    auto cloned = std::make_unique<LinearInterpolation>(_xData, _yData);
  
    // Clone extrapolation scheme 
    cloned->setExtrapolationScheme(_extrapolationScheme->clone());
  
    return cloned;  // Ownership transferred to caller
}
```

**What's happening:**

- `clone()` creates a new object and returns it
- Return value transfers ownership to caller
- **No copying** - ownership is moved
- Caller becomes the owner

**Usage:**

```cpp
auto original = std::make_unique<LinearInterpolation>(xData, yData);
auto cloned = original->clone();  // Ownership transferred to 'cloned'
// 'cloned' now owns the new object
```

### Key Operations with `unique_ptr`

```cpp
// 1. Create
auto ptr = std::make_unique<MyClass>(arg1, arg2);

// 2. Check if null
if (ptr) { /* ptr is not null */ }
if (!ptr) { /* ptr is null */ }

// 3. Dereference
MyClass& obj = *ptr;        // Get reference
MyClass* raw = ptr.get();   // Get raw pointer (don't delete!)

// 4. Reset (delete and set to null)
ptr.reset();                 // Delete object, set to nullptr
ptr.reset(new MyClass());    // Delete old, set to new

// 5. Release ownership (get raw pointer, unique_ptr becomes null)
MyClass* raw = ptr.release(); // You must delete raw manually!

// 6. Move (transfer ownership)
auto ptr2 = std::move(ptr);   // ptr is now null, ptr2 owns object
```

### Why `unique_ptr` Cannot Be Copied

```cpp
auto ptr1 = std::make_unique<MyClass>();
// auto ptr2 = ptr1;  // ERROR! Cannot copy unique_ptr
auto ptr2 = std::move(ptr1);  // OK! Transfer ownership
// Now ptr1 is null, ptr2 owns the object
```

**Reason:** Prevents double deletion. If you could copy, both would try to delete the same object.

---

## std::make_unique

### What is `std::make_unique`?

`std::make_unique` is a **factory function** that creates a `std::unique_ptr`. It was introduced in C++14.

### Syntax

```cpp
auto ptr = std::make_unique<Type>(arg1, arg2, ...);
```

### Why Use `make_unique`?

#### 1. Exception Safety

**BAD (potential leak):**

```cpp
void func() {
    process(std::unique_ptr<MyClass>(new MyClass()), 
            std::unique_ptr<MyClass>(new MyClass()));
    // If second new throws, first object leaks!
}
```

**GOOD (exception-safe):**

```cpp
void func() {
    process(std::make_unique<MyClass>(), 
            std::make_unique<MyClass>());
    // If second make_unique throws, first is automatically cleaned up
}
```

#### 2. Cleaner Code

**BAD:**

```cpp
std::unique_ptr<MyClass> ptr(new MyClass(arg1, arg2));
```

**GOOD:**

```cpp
auto ptr = std::make_unique<MyClass>(arg1, arg2);
```

#### 3. No Need for `new`/`delete`

You never write `new` or `delete` with `make_unique`!

### Examples from Your Codebase

**Location:** `InterpolationSchemes.cpp` line 120

```cpp
auto cloned = std::make_unique<LinearInterpolation>(_xData, _yData);
```

**What's happening:**

- Creates `LinearInterpolation` with constructor arguments `_xData, _yData`
- Returns `std::unique_ptr<LinearInterpolation>`
- Exception-safe and clean

**Location:** `main.cpp` line 61-66

```cpp
std::unique_ptr<ExtrapolationScheme> scheme;
if (choice == 1) {
    scheme = std::make_unique<FlatExtrapolation>();
} else if (choice == 2) {
    scheme = std::make_unique<LinearExtrapolation>();
} else {
    scheme = std::make_unique<QuadraticExtrapolation>();
}
```

**What's happening:**

- Creates different extrapolation schemes based on user choice
- All use `make_unique` for safety and clarity

---

## std::shared_ptr

### What is `std::shared_ptr`?

`std::shared_ptr` is a **smart pointer that shares ownership** of an object. Multiple `shared_ptr` instances can point to the same object.

**Key features:**

- **Shared ownership**: Multiple `shared_ptr` can own the same object
- **Reference counting**: Tracks how many `shared_ptr` point to the object
- **Automatic deletion**: Object is deleted when **last** `shared_ptr` is destroyed
- **Copyable**: Can be copied (increases reference count)
- **Thread-safe**: Reference counting is thread-safe (C++11)

### How Reference Counting Works

```cpp
{
    auto ptr1 = std::make_shared<MyClass>();  // Reference count = 1
    {
        auto ptr2 = ptr1;                     // Reference count = 2
        {
            auto ptr3 = ptr1;                 // Reference count = 3
        }  // ptr3 destroyed, reference count = 2
    }  // ptr2 destroyed, reference count = 1
}  // ptr1 destroyed, reference count = 0 -> object deleted!
```

### Example Usage (Not in Your Codebase, but Useful)

```cpp
class Node {
public:
    std::shared_ptr<Node> parent;
    std::shared_ptr<Node> child;
    // ...
};

void example() {
    auto parent = std::make_shared<Node>();
    auto child = std::make_shared<Node>();
  
    parent->child = child;      // parent owns child
    child->parent = parent;     // child owns parent
    // Circular reference! Both stay alive
  
    // When both go out of scope, both are deleted
}
```

### Operations with `shared_ptr`

```cpp
// 1. Create
auto ptr1 = std::make_shared<MyClass>(arg1, arg2);

// 2. Copy (increases reference count)
auto ptr2 = ptr1;  // Reference count = 2

// 3. Check reference count
size_t count = ptr1.use_count();  // Returns number of shared_ptr pointing to object

// 4. Reset
ptr1.reset();  // Decreases reference count, sets ptr1 to null

// 5. Check if unique
if (ptr1.unique()) { /* Only one shared_ptr owns object */ }
```

### When to Use `shared_ptr`

**Use `shared_ptr` when:**

- Multiple objects need to share ownership
- You don't know which object will outlive the others
- You need shared ownership semantics

**Example scenarios:**

- Shared resources (file handles, network connections)
- Observer patterns
- Caching systems
- Graph structures (with careful design to avoid cycles)

---

## std::make_shared

### What is `std::make_shared`?

`std::make_shared` is a **factory function** that creates a `std::shared_ptr`. Introduced in C++11.

### Syntax

```cpp
auto ptr = std::make_shared<Type>(arg1, arg2, ...);
```

### Why Use `make_shared`?

#### 1. More Efficient

`make_shared` allocates the object and control block (for reference counting) in **one allocation**, while `new` + `shared_ptr` constructor does **two allocations**.

**Less efficient:**

```cpp
std::shared_ptr<MyClass> ptr(new MyClass(arg1, arg2));  // Two allocations
```

**More efficient:**

```cpp
auto ptr = std::make_shared<MyClass>(arg1, arg2);  // One allocation
```

#### 2. Exception Safety

Same benefits as `make_unique`.

### Example Usage

```cpp
class Resource {
    // Expensive resource
};

class Manager {
    std::shared_ptr<Resource> _resource;
public:
    Manager() : _resource(std::make_shared<Resource>()) {}
  
    std::shared_ptr<Resource> getResource() const {
        return _resource;  // Returns shared_ptr (increases ref count)
    }
};

void example() {
    Manager manager;
    auto resource1 = manager.getResource();  // Ref count = 2
    auto resource2 = manager.getResource();  // Ref count = 3
  
    // When all go out of scope, Resource is deleted
}
```

---

## std::weak_ptr

### What is `std::weak_ptr`?

`std::weak_ptr` is a **non-owning smart pointer** that holds a reference to an object managed by `std::shared_ptr`. It doesn't affect the reference count.

**Key features:**

- **Non-owning**: Doesn't keep object alive
- **No reference count**: Doesn't increase `shared_ptr` reference count
- **Expires**: Can become invalid if object is deleted
- **Safe access**: Must convert to `shared_ptr` to access object

### Why Use `weak_ptr`?

**Problem:** Circular references with `shared_ptr`

```cpp
class Parent {
public:
    std::shared_ptr<Child> child;  // Parent owns Child
};

class Child {
public:
    std::shared_ptr<Parent> parent;  // Child owns Parent
    // CIRCULAR REFERENCE! Neither is ever deleted!
};

void bad_example() {
    auto parent = std::make_shared<Parent>();
    auto child = std::make_shared<Child>();
    parent->child = child;
    child->parent = parent;  // Circular reference!
    // When function ends, both still exist (ref count = 1 each)
    // MEMORY LEAK!
}
```

**Solution:** Use `weak_ptr` to break the cycle

```cpp
class Parent {
public:
    std::shared_ptr<Child> child;
};

class Child {
public:
    std::weak_ptr<Parent> parent;  // Weak reference - doesn't keep Parent alive
};

void good_example() {
    auto parent = std::make_shared<Parent>();
    auto child = std::make_shared<Child>();
    parent->child = child;
    child->parent = parent;  // Weak reference - no cycle!
    // When function ends, parent is deleted (ref count = 0)
    // Then child is deleted (ref count = 0)
}
```

### Operations with `weak_ptr`

```cpp
// 1. Create from shared_ptr
auto shared = std::make_shared<MyClass>();
std::weak_ptr<MyClass> weak = shared;  // Doesn't increase ref count

// 2. Check if expired
if (weak.expired()) { /* Object was deleted */ }

// 3. Lock (convert to shared_ptr)
if (auto locked = weak.lock()) {
    // Object still exists, locked is a shared_ptr
    // Can use locked safely
} else {
    // Object was deleted
}

// 4. Get shared_ptr (throws if expired)
try {
    auto shared = weak.lock();  // Returns shared_ptr or nullptr
    if (shared) {
        // Use shared
    }
} catch (...) {
    // Handle error
}
```

### Example: Observer Pattern

```cpp
class Subject {
    std::vector<std::weak_ptr<Observer>> observers;
public:
    void attach(std::weak_ptr<Observer> obs) {
        observers.push_back(obs);
    }
  
    void notify() {
        for (auto it = observers.begin(); it != observers.end();) {
            if (auto obs = it->lock()) {
                obs->update();  // Observer still exists
                ++it;
            } else {
                it = observers.erase(it);  // Observer was deleted
            }
        }
    }
};
```

---

## When to Use Which?

### Decision Tree

```
Do you need shared ownership?
├─ NO → Use std::unique_ptr
│   └─ Single owner, automatic cleanup
│
└─ YES → Use std::shared_ptr
    └─ Do you have circular references?
        ├─ NO → Use std::shared_ptr everywhere
        └─ YES → Use std::weak_ptr to break cycles
```

### Summary Table

| Smart Pointer       | Ownership | Copyable     | Use Case                             |
| ------------------- | --------- | ------------ | ------------------------------------ |
| `std::unique_ptr` | Exclusive | No (movable) | Single owner, automatic cleanup      |
| `std::shared_ptr` | Shared    | Yes          | Multiple owners, shared resources    |
| `std::weak_ptr`   | None      | Yes          | Break circular references, observers |

### Your Codebase Analysis

**Your codebase uses `std::unique_ptr` exclusively**, which is **excellent** because:

1. **Clear ownership**: Each object has one clear owner
2. **No circular references**: No risk of memory leaks from cycles
3. **Better performance**: No reference counting overhead
4. **Simpler design**: Easier to reason about

**When you might need `shared_ptr` in the future:**

- If you need to share an `InterpolationScheme` between multiple objects
- If you implement observer patterns
- If you have complex ownership relationships

**When you might need `weak_ptr`:**

- If you add `shared_ptr` and encounter circular references
- If you implement caching systems
- If you need to break ownership cycles

---

## Common Patterns in Your Codebase

### Pattern 1: Member Variable Storage

**Location:** `InterpolationSchemes.h` line 54

```cpp
class InterpolationScheme {
protected:
    std::unique_ptr<ExtrapolationScheme> _extrapolationScheme;
};
```

**Pattern:** Store owned objects as `unique_ptr` members
**Benefit:** Automatic cleanup, clear ownership

### Pattern 2: Factory Functions Returning `unique_ptr`

**Location:** `InterpolationSchemes.cpp` line 118

```cpp
std::unique_ptr<InterpolationScheme> LinearInterpolation::clone() const
{
    auto cloned = std::make_unique<LinearInterpolation>(_xData, _yData);
    cloned->setExtrapolationScheme(_extrapolationScheme->clone());
    return cloned;  // Ownership transferred
}
```

**Pattern:** Return `unique_ptr` from factory/clone functions
**Benefit:** Caller becomes owner, no memory leaks

### Pattern 3: Transferring Ownership in Constructors

**Location:** `InterpolationSchemes.cpp` line 67-73

```cpp
LinearInterpolation::LinearInterpolation(
    const std::vector<double>& xData, 
    const std::vector<double>& yData,
    std::unique_ptr<ExtrapolationScheme> extrapolationScheme)
    : LinearInterpolation(xData, yData)
{
    setExtrapolationScheme(std::move(extrapolationScheme));
}
```

**Pattern:** Constructor takes `unique_ptr` by value, transfers ownership
**Benefit:** Clear ownership transfer, exception-safe

### Pattern 4: Lazy Initialization with `mutable`

**Location:** `PDEs/Solver.h` line 82

```cpp
mutable std::unique_ptr<CubicSplineInterpolation> _interpolator;
```

**Pattern:** `mutable unique_ptr` for lazy initialization
**Benefit:** Can be initialized in `const` methods, created only when needed

**Usage:**

```cpp
const CubicSplineInterpolation& Solver::getInterpolator() const
{
    if (!_interpolator) {
        _interpolator = std::make_unique<CubicSplineInterpolation>(...);
    }
    return *_interpolator;
}
```

### Pattern 5: Container of `unique_ptr`

**Location:** `VolatilitySurface.h` line 87-88

```cpp
std::vector<std::unique_ptr<InterpolationScheme>> _smileInterpolators;
std::vector<std::unique_ptr<InterpolationScheme>> _termStructureInterpolators;
```

**Pattern:** Store polymorphic objects in containers
**Benefit:** Polymorphism, automatic cleanup, clear ownership

---

## Best Practices

### ✅ DO

1. **Use `make_unique` and `make_shared`** instead of `new`

   ```cpp
   auto ptr = std::make_unique<MyClass>();  // ✅
   // Not: std::unique_ptr<MyClass> ptr(new MyClass());  // ❌
   ```
2. **Use `unique_ptr` by default** - it's simpler and faster

   ```cpp
   std::unique_ptr<MyClass> ptr;  // ✅ Default choice
   ```
3. **Transfer ownership with `std::move`**

   ```cpp
   auto ptr2 = std::move(ptr1);  // ✅ Clear ownership transfer
   ```
4. **Return `unique_ptr` from functions** to transfer ownership

   ```cpp
   std::unique_ptr<MyClass> create() {
       return std::make_unique<MyClass>();
   }
   ```
5. **Use `weak_ptr` to break circular references**

   ```cpp
   std::weak_ptr<Parent> parent;  // ✅ Breaks cycle
   ```

### ❌ DON'T

1. **Don't mix raw pointers and smart pointers**

   ```cpp
   MyClass* raw = ptr.get();
   delete raw;  // ❌ DON'T! unique_ptr will delete it
   ```
2. **Don't create `shared_ptr` from raw pointer multiple times**

   ```cpp
   MyClass* raw = new MyClass();
   std::shared_ptr<MyClass> p1(raw);
   std::shared_ptr<MyClass> p2(raw);  // ❌ Double deletion!
   ```
3. **Don't use `shared_ptr` when `unique_ptr` is sufficient**

   ```cpp
   std::shared_ptr<MyClass> ptr;  // ❌ Unnecessary overhead
   std::unique_ptr<MyClass> ptr;  // ✅ Better choice
   ```
4. **Don't store raw pointers to objects managed by smart pointers**

   ```cpp
   MyClass* raw = ptr.get();
   // Later: ptr is destroyed, raw is dangling! ❌
   ```

---

## Memory Management Comparison

### Raw Pointers (Old Way)

```cpp
ExtrapolationScheme* scheme = new FlatExtrapolation();
// ... use scheme ...
delete scheme;  // Must remember!
// What if exception occurs? LEAK!
```

### Smart Pointers (Modern Way)

```cpp
auto scheme = std::make_unique<FlatExtrapolation>();
// ... use scheme ...
// Automatically deleted when out of scope!
// Exception-safe!
```

---

## Interview Quick Reference: Smart Pointers

### Essential Code Patterns

```cpp
// ============================================
// CREATION
// ============================================
auto uptr = std::make_unique<MyClass>(args);     // unique_ptr
auto sptr = std::make_shared<MyClass>(args);     // shared_ptr
std::weak_ptr<MyClass> wptr = sptr;            // weak_ptr

// ============================================
// USAGE
// ============================================
if (uptr) uptr->method();                        // unique_ptr
if (sptr) sptr->method();                        // shared_ptr
if (auto locked = wptr.lock()) locked->method(); // weak_ptr

// ============================================
// OWNERSHIP TRANSFER
// ============================================
auto uptr2 = std::move(uptr);                    // unique_ptr (move)
auto sptr2 = sptr;                                // shared_ptr (copy, ref count++)

// ============================================
// OPERATIONS
// ============================================
uptr.reset();                                     // Delete, set null
sptr.reset();                                     // Decrease ref count
MyClass* raw = uptr.get();                        // Get raw (don't delete!)
size_t count = sptr.use_count();                  // Reference count
if (wptr.expired()) { /* deleted */ }             // Check if valid
```

---

### Key Interview Answers

**Q: What is `std::unique_ptr`?**

> Exclusive ownership. One owner, auto-deletes. Cannot copy, can move. Use `make_unique`.

**Q: What is `std::shared_ptr`?**

> Shared ownership. Reference counting. Deleted when last `shared_ptr` destroyed. Copyable. Use `make_shared`.

**Q: What is `std::weak_ptr`?**

> Non-owning reference to `shared_ptr`. Doesn't affect ref count. Use `lock()` to access. Breaks circular references.

**Q: Why `make_unique`/`make_shared`?**

> Exception-safe. Single allocation (`make_shared`). Cleaner code. Never use `new`/`delete`.

**Q: When to use which?**

> Default: `unique_ptr`. Shared ownership: `shared_ptr`. Break cycles: `weak_ptr`.

---

### Common Patterns

```cpp
// Pattern 1: Member variable
class Container {
    std::unique_ptr<MyClass> _member;  // Auto-deleted
};

// Pattern 2: Factory function
std::unique_ptr<MyClass> create() {
    return std::make_unique<MyClass>();
}

// Pattern 3: Transfer in constructor
Container(std::unique_ptr<MyClass> ptr) 
    : _member(std::move(ptr)) {}

// Pattern 4: Break circular reference
class Child {
    std::weak_ptr<Parent> parent;  // Weak breaks cycle
    void use() {
        if (auto p = parent.lock()) p->doSomething();
    }
};
```

---

### Reference Counting Example

```cpp
auto sptr1 = std::make_shared<MyClass>();  // count = 1
auto sptr2 = sptr1;                        // count = 2
auto sptr3 = sptr1;                        // count = 3
sptr2.reset();                             // count = 2
sptr3.reset();                             // count = 1
sptr1.reset();                             // count = 0 → deleted
```

---

## Simplified Implementations (Interview Blueprints)

### 1. `std::unique_ptr` - Simplified Implementation

```cpp
template<typename T>
class unique_ptr {
public:
    // Constructor
    explicit unique_ptr(T* ptr = nullptr) : _rawPtr(ptr) {}
  
    // Destructor - auto cleanup
    ~unique_ptr() {
        delete _rawPtr;
    }
  
    // Delete copy (cannot copy unique_ptr)
    unique_ptr(const unique_ptr&) = delete;
    unique_ptr& operator=(const unique_ptr&) = delete;
  
    // Move constructor - transfer ownership
    unique_ptr(unique_ptr&& other) noexcept 
        : _rawPtr(other._rawPtr) {
        other._rawPtr = nullptr;  // Source becomes null
    }
  
    // Move assignment - transfer ownership
    unique_ptr& operator=(unique_ptr&& other) noexcept {
        if (this != &other) {
            delete _rawPtr;           // Clean up current resource
            _rawPtr = other._rawPtr;  // Steal ownership
            other._rawPtr = nullptr;  // Source becomes null
        }
        return *this;
    }
  
    // Access operators
    T& operator*() const { return *_rawPtr; }
    T* operator->() const { return _rawPtr; }
  
    // Check if valid
    explicit operator bool() const { return _rawPtr != nullptr; }
  
    // Get raw pointer (don't delete!)
    T* get() const { return _rawPtr; }
  
    // Release ownership (returns raw, sets to null)
    T* release() {
        T* temp = _rawPtr;
        _rawPtr = nullptr;
        return temp;
    }
  
    // Reset (delete and set to null)
    void reset(T* ptr = nullptr) {
        delete _rawPtr;
        _rawPtr = ptr;
    }

private:
    T* _rawPtr;
};

// Factory function (make_unique)
template<typename T, typename... Args>
unique_ptr<T> make_unique(Args&&... args) {
    return unique_ptr<T>(new T(std::forward<Args>(args)...));
}
```

**Key Points:**

- **Move semantics**: `&&` rvalue reference, `noexcept` for optimization
- **No copy**: Deleted copy constructor/assignment
- **Auto cleanup**: Destructor deletes pointer
- **Ownership transfer**: Move steals pointer, sets source to null

---

### 2. `std::shared_ptr` - Simplified Implementation

```cpp
template<typename T>
class shared_ptr {
private:
    T* _rawPtr;
    size_t* _refCount;  // Reference counter (shared between all copies)
  
    void cleanup() {
        if (_refCount) {
            (*_refCount)--;
            if (*_refCount == 0) {
                delete _rawPtr;
                delete _refCount;
            }
        }
    }

public:
    // Constructor
    explicit shared_ptr(T* ptr = nullptr) 
        : _rawPtr(ptr), _refCount(ptr ? new size_t(1) : nullptr) {}
  
    // Copy constructor - share ownership
    shared_ptr(const shared_ptr& other) 
        : _rawPtr(other._rawPtr), _refCount(other._refCount) {
        if (_refCount) (*_refCount)++;
    }
  
    // Copy assignment
    shared_ptr& operator=(const shared_ptr& other) {
        if (this != &other) {
            cleanup();  // Decrease old ref count
            _rawPtr = other._rawPtr;
            _refCount = other._refCount;
            if (_refCount) (*_refCount)++;
        }
        return *this;
    }
  
    // Destructor - decrease ref count, delete if last
    ~shared_ptr() {
        cleanup();
    }
  
    // Access operators
    T& operator*() const { return *_rawPtr; }
    T* operator->() const { return _rawPtr; }
  
    // Check if valid
    explicit operator bool() const { return _rawPtr != nullptr; }
  
    // Get raw pointer
    T* get() const { return _rawPtr; }
  
    // Get reference count
    size_t use_count() const { return _refCount ? *_refCount : 0; }
  
    // Check if unique owner
    bool unique() const { return use_count() == 1; }
  
    // Reset
    void reset(T* ptr = nullptr) {
        cleanup();
        _rawPtr = ptr;
        _refCount = ptr ? new size_t(1) : nullptr;
    }
};

// Factory function (make_shared)
template<typename T, typename... Args>
shared_ptr<T> make_shared(Args&&... args) {
    return shared_ptr<T>(new T(std::forward<Args>(args)...));
}
```

**Key Points:**

- **Reference counting**: Shared counter between all copies
- **Copyable**: Copy increases ref count
- **Auto cleanup**: Deletes when ref count reaches 0
- **Thread-safe**: Real `shared_ptr` uses atomic operations

---

### 3. `std::weak_ptr` - Simplified Implementation

```cpp
template<typename T>
class weak_ptr {
private:
    T* _rawPtr;
    size_t* _refCount;      // Points to shared_ptr's ref count (doesn't own)
    size_t* _weakCount;     // Weak reference counter
  
    void cleanup() {
        if (_weakCount) {
            (*_weakCount)--;
            // Control block deleted when ref_count = 0 AND weak_count = 0
        }
    }

public:
    // Constructor from shared_ptr
    weak_ptr(const shared_ptr<T>& sptr) 
        : _rawPtr(sptr.get()),
          _refCount(sptr._getRefCount()),  // Get ref count pointer
          _weakCount(new size_t(1)) {
        // Doesn't increment ref_count! Only tracks weak references
    }
  
    // Copy constructor
    weak_ptr(const weak_ptr& other) 
        : _rawPtr(other._rawPtr),
          _refCount(other._refCount),
          _weakCount(other._weakCount) {
        if (_weakCount) (*_weakCount)++;
    }
  
    // Destructor
    ~weak_ptr() {
        cleanup();
    }
  
    // Check if expired (object was deleted)
    bool expired() const {
        return !_refCount || *_refCount == 0;
    }
  
    // Lock - convert to shared_ptr
    shared_ptr<T> lock() const {
        if (expired()) {
            return shared_ptr<T>();  // Return null shared_ptr
        }
        // Create shared_ptr that shares same control block
        // Increases ref_count
        return shared_ptr<T>(_rawPtr, _refCount);
    }
  
    // Get reference count (of shared_ptrs)
    size_t use_count() const {
        return _refCount ? *_refCount : 0;
    }
};

// Key: weak_ptr doesn't keep object alive (ref_count = 0 → deleted)
//      But control block lives until weak_count also reaches 0
```

**Key Points:**

- **Non-owning**: Doesn't affect `shared_ptr` ref count
- **Expired check**: Checks if object still exists
- **Lock**: Converts to `shared_ptr` for access
- **Weak count**: Tracks weak references separately

---

### 4. Move Semantics Deep Dive

```cpp
// ============================================
// RVALUE REFERENCE: &&
// ============================================
// && means "rvalue reference" - binds to temporaries
unique_ptr(unique_ptr&& other)  // Move constructor
// std::move(x) converts x to rvalue, enabling move

// ============================================
// NOEXCEPT: Why it matters
// ============================================
// noexcept enables optimizations:
// - std::vector uses move instead of copy
// - Won't fall back to copy if move throws
unique_ptr(unique_ptr&& other) noexcept

// ============================================
// MOVE vs COPY
// ============================================
// Copy: Expensive, duplicates resources, both usable
unique_ptr(const unique_ptr&) = delete;  // Cannot copy!

// Move: Cheap, transfers ownership, source empty
unique_ptr(unique_ptr&& other) noexcept {
    _rawPtr = other._rawPtr;      // Steal pointer
    other._rawPtr = nullptr;      // Source becomes null
}

// ============================================
// SELF-ASSIGNMENT CHECK
// ============================================
unique_ptr& operator=(unique_ptr&& other) noexcept {
    if (this != &other) {  // Prevent moving to self
        delete _rawPtr;
        _rawPtr = other._rawPtr;
        other._rawPtr = nullptr;
    }
    return *this;
}
```

**Interview Answer:**

> "Move transfers ownership without copying. `&&` is rvalue reference, `noexcept` enables optimizations. Move constructor steals pointer and sets source to null. Self-assignment check prevents bugs."

---

### 5. Reference Counting Deep Dive

```cpp
// ============================================
// HOW REFERENCE COUNTING WORKS
// ============================================

// Create shared_ptr
shared_ptr<MyClass> p1(new MyClass());  
// _refCount = new size_t(1)  // count = 1

// Copy increases count
shared_ptr<MyClass> p2 = p1;  
// _refCount = p1._refCount   // Same counter!
// (*_refCount)++              // count = 2

// Another copy
shared_ptr<MyClass> p3 = p1;  
// (*_refCount)++              // count = 3

// Reset decreases count
p2.reset();  
// (*_refCount)--              // count = 2

// Destructor decreases count
// When p3 destroyed: count = 1
// When p1 destroyed: count = 0 → delete object

// ============================================
// CONTROL BLOCK (Real Implementation)
// ============================================
// Real shared_ptr uses control block:
struct ControlBlock {
    size_t ref_count;    // Shared references
    size_t weak_count;   // Weak references
    T* object;           // Actual object (or nullptr)
};

// make_shared allocates object + control block together
// More efficient than separate allocations
```

**Interview Answer:**

> "`shared_ptr` uses reference counting. All copies share same counter. Copy increments, destructor decrements. Deleted when count reaches 0. Control block stores ref count and object. `make_shared` allocates both together for efficiency."

---

## Summary

### Key Takeaways

1. **`std::unique_ptr`**: Exclusive ownership, use by default
2. **`std::make_unique`**: Safe factory function for `unique_ptr`
3. **`std::shared_ptr`**: Shared ownership, use when needed
4. **`std::make_shared`**: Efficient factory function for `shared_ptr`
5. **`std::weak_ptr`**: Break circular references, non-owning

### Your Codebase Uses

- ✅ `std::unique_ptr` - Excellent choice for your design
- ✅ `std::make_unique` - Safe and clean
- ❌ `std::shared_ptr` - Not needed (good!)
- ❌ `std::weak_ptr` - Not needed (good!)

**Your codebase demonstrates excellent modern C++ practices!**

---

## Further Reading

- [cppreference.com - unique_ptr](https://en.cppreference.com/w/cpp/memory/unique_ptr)
- [cppreference.com - shared_ptr](https://en.cppreference.com/w/cpp/memory/shared_ptr)
- [cppreference.com - weak_ptr](https://en.cppreference.com/w/cpp/memory/weak_ptr)
- [Herb Sutter&#39;s GotW #89: Smart Pointers](https://herbsutter.com/2013/06/05/gotw-89-solution-smart-pointers/)

---

*Last updated: Based on CppFM codebase analysis*
