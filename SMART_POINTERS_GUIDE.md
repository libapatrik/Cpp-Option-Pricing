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
9. [std::vector](#stdvector)
10. [operator() - Function Call Operator](#operator-function-call-operator)
11. [operator[] - Subscript Operator](#operator-subscript-operator)

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

### Container & Operator Summary Table

| Feature | Complexity | Bounds Check | Use Case |
| ------- | ---------- | ------------ | -------- |
| `std::vector::operator[]` | O(1) | No | Fast array access |
| `std::vector::at()` | O(1) | Yes (throws) | Safe array access |
| `std::vector::push_back()` | O(1) amortized | N/A | Add element |
| `operator()` | O(1) | N/A | Make objects callable |
| `operator[]` | O(1) typically | Depends on impl | Array-like access |

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

### Pattern 6: Implementing `operator()` for Callable Objects

```cpp
class Comparator {
    int _threshold;
public:
    Comparator(int threshold) : _threshold(threshold) {}
    
    bool operator()(int value) const {
        return value > _threshold;
    }
};

// Usage with STL
std::vector<int> vec = {1, 5, 10, 15};
Comparator gt5(5);
auto it = std::find_if(vec.begin(), vec.end(), gt5);
```

**Pattern:** Implement `operator()` to make objects callable
**Benefit:** Stateful predicates, reusable logic, STL compatibility

### Pattern 7: Implementing `operator[]` for Container Access

```cpp
class SafeArray {
    int* _data;
    size_t _size;
public:
    // Non-const version (modification)
    int& operator[](size_t index) {
        return _data[index];
    }
    
    // Const version (read-only)
    const int& operator[](size_t index) const {
        return _data[index];
    }
};
```

**Pattern:** Provide both const and non-const `operator[]`
**Benefit:** Array-like access, const correctness, intuitive interface

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

## Interview Quick Reference: Smart Pointers & Containers

### Essential Code Patterns

```cpp
// ============================================
// SMART POINTERS - CREATION
// ============================================
auto uptr = std::make_unique<MyClass>(args);     // unique_ptr
auto sptr = std::make_shared<MyClass>(args);     // shared_ptr
std::weak_ptr<MyClass> wptr = sptr;            // weak_ptr

// ============================================
// SMART POINTERS - USAGE
// ============================================
if (uptr) uptr->method();                        // unique_ptr
if (sptr) sptr->method();                        // shared_ptr
if (auto locked = wptr.lock()) locked->method(); // weak_ptr

// ============================================
// SMART POINTERS - OWNERSHIP TRANSFER
// ============================================
auto uptr2 = std::move(uptr);                    // unique_ptr (move)
auto sptr2 = sptr;                                // shared_ptr (copy, ref count++)

// ============================================
// SMART POINTERS - OPERATIONS
// ============================================
uptr.reset();                                     // Delete, set null
sptr.reset();                                     // Decrease ref count
MyClass* raw = uptr.get();                        // Get raw (don't delete!)
size_t count = sptr.use_count();                  // Reference count
if (wptr.expired()) { /* deleted */ }             // Check if valid

// ============================================
// VECTOR - OPERATIONS
// ============================================
std::vector<int> vec;
vec.push_back(42);                                // Add element
vec[0] = 10;                                      // Access/modify
vec.size();                                       // Get size
vec.capacity();                                   // Get capacity
vec.reserve(100);                                 // Reserve memory
vec.resize(50);                                   // Resize

// ============================================
// OPERATOR() - FUNCTOR
// ============================================
class Adder {
    int _value;
public:
    Adder(int v) : _value(v) {}
    int operator()(int x) const { return x + _value; }
};
Adder add5(5);
int result = add5(10);                           // Calls operator()

// ============================================
// OPERATOR[] - SUBSCRIPT
// ============================================
class Array {
    int* _data;
    size_t _size;
public:
    int& operator[](size_t i) { return _data[i]; }
    const int& operator[](size_t i) const { return _data[i]; }
};
Array arr(10);
arr[0] = 42;                                      // Modify
int val = arr[0];                                 // Read
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

**Q: How does `std::vector` work internally?**

> Dynamic array with capacity management. Stores pointer, size, capacity. Doubles capacity when full for amortized O(1) push_back. Uses placement new for construction, explicit destructor calls.

**Q: What is `operator()`?**

> Function call operator makes objects callable. Implement `ReturnType operator()(Params...)`. Can maintain state, overload for different signatures. Used with STL algorithms.

**Q: What is `operator[]`?**

> Subscript operator provides array-like access. Implement `T& operator[](Index)` and `const T& operator[](Index) const`. Usually no bounds checking (use `at()` for that). Should be O(1).

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

### 6. `std::vector` - Simplified Implementation

```cpp
template<typename T>
class vector {
private:
    T* _data;              // Pointer to dynamically allocated array
    size_t _size;          // Number of elements currently stored
    size_t _capacity;      // Total capacity of allocated memory

    void reallocate(size_t newCapacity) {
        T* newData = static_cast<T*>(::operator new(newCapacity * sizeof(T)));
        
        // Move existing elements to new memory
        size_t moveCount = _size < newCapacity ? _size : newCapacity;
        for (size_t i = 0; i < moveCount; ++i) {
            new (newData + i) T(std::move(_data[i]));  // Placement new + move
            _data[i].~T();  // Destroy old element
        }
        
        // Destroy remaining old elements if shrinking
        for (size_t i = moveCount; i < _size; ++i) {
            _data[i].~T();
        }
        
        ::operator delete(_data);
        _data = newData;
        _capacity = newCapacity;
        _size = moveCount;
    }

public:
    // Default constructor
    vector() : _data(nullptr), _size(0), _capacity(0) {}
    
    // Constructor with initial size
    explicit vector(size_t count) : _size(count), _capacity(count) {
        _data = static_cast<T*>(::operator new(_capacity * sizeof(T)));
        for (size_t i = 0; i < _size; ++i) {
            new (_data + i) T();  // Default construct each element
        }
    }
    
    // Copy constructor
    vector(const vector& other) : _size(other._size), _capacity(other._capacity) {
        _data = static_cast<T*>(::operator new(_capacity * sizeof(T)));
        for (size_t i = 0; i < _size; ++i) {
            new (_data + i) T(other._data[i]);  // Copy construct
        }
    }
    
    // Move constructor
    vector(vector&& other) noexcept 
        : _data(other._data), _size(other._size), _capacity(other._capacity) {
        other._data = nullptr;
        other._size = 0;
        other._capacity = 0;
    }
    
    // Destructor
    ~vector() {
        clear();
        ::operator delete(_data);
    }
    
    // Copy assignment
    vector& operator=(const vector& other) {
        if (this != &other) {
            clear();
            if (_capacity < other._size) {
                reallocate(other._size);
            }
            _size = other._size;
            for (size_t i = 0; i < _size; ++i) {
                new (_data + i) T(other._data[i]);
            }
        }
        return *this;
    }
    
    // Move assignment
    vector& operator=(vector&& other) noexcept {
        if (this != &other) {
            clear();
            ::operator delete(_data);
            _data = other._data;
            _size = other._size;
            _capacity = other._capacity;
            other._data = nullptr;
            other._size = 0;
            other._capacity = 0;
        }
        return *this;
    }
    
    // Subscript operator (non-const)
    T& operator[](size_t index) {
        return _data[index];
    }
    
    // Subscript operator (const)
    const T& operator[](size_t index) const {
        return _data[index];
    }
    
    // At with bounds checking
    T& at(size_t index) {
        if (index >= _size) {
            throw std::out_of_range("Index out of range");
        }
        return _data[index];
    }
    
    const T& at(size_t index) const {
        if (index >= _size) {
            throw std::out_of_range("Index out of range");
        }
        return _data[index];
    }
    
    // Size
    size_t size() const { return _size; }
    
    // Capacity
    size_t capacity() const { return _capacity; }
    
    // Empty
    bool empty() const { return _size == 0; }
    
    // Reserve
    void reserve(size_t newCapacity) {
        if (newCapacity > _capacity) {
            reallocate(newCapacity);
        }
    }
    
    // Resize
    void resize(size_t newSize) {
        if (newSize > _capacity) {
            reserve(newSize * 2);  // Common growth strategy
        }
        if (newSize > _size) {
            // Construct new elements
            for (size_t i = _size; i < newSize; ++i) {
                new (_data + i) T();
            }
        } else {
            // Destroy excess elements
            for (size_t i = newSize; i < _size; ++i) {
                _data[i].~T();
            }
        }
        _size = newSize;
    }
    
    // Push back
    void push_back(const T& value) {
        if (_size >= _capacity) {
            reserve(_capacity == 0 ? 1 : _capacity * 2);  // Double capacity
        }
        new (_data + _size) T(value);  // Copy construct
        ++_size;
    }
    
    void push_back(T&& value) {
        if (_size >= _capacity) {
            reserve(_capacity == 0 ? 1 : _capacity * 2);
        }
        new (_data + _size) T(std::move(value));  // Move construct
        ++_size;
    }
    
    // Pop back
    void pop_back() {
        if (_size > 0) {
            _data[_size - 1].~T();
            --_size;
        }
    }
    
    // Clear
    void clear() {
        for (size_t i = 0; i < _size; ++i) {
            _data[i].~T();
        }
        _size = 0;
    }
    
    // Front
    T& front() { return _data[0]; }
    const T& front() const { return _data[0]; }
    
    // Back
    T& back() { return _data[_size - 1]; }
    const T& back() const { return _data[_size - 1]; }
    
    // Data (get raw pointer)
    T* data() { return _data; }
    const T* data() const { return _data; }
};
```

**Key Points:**

- **Dynamic allocation**: Uses `operator new`/`operator delete` for raw memory
- **Placement new**: Constructs objects in-place using `new (ptr) T(...)`
- **Explicit destruction**: Calls `~T()` before deallocating
- **Growth strategy**: Doubles capacity when full (amortized O(1) push_back)
- **Move semantics**: Efficiently transfers ownership
- **Exception safety**: Proper cleanup in destructor

**Interview Answer:**

> "`vector` uses dynamic array with capacity management. Stores pointer, size, capacity. Doubles capacity when full for amortized O(1) push_back. Uses placement new for construction, explicit destructor calls. Move constructor transfers pointer ownership. Subscript operator provides O(1) access."

---

### 7. `operator()` - Function Call Operator Implementation

The function call operator `operator()` allows objects to be called like functions. These are called **functors** or **function objects**.

#### Basic Implementation

```cpp
class Adder {
private:
    int _value;
public:
    Adder(int value) : _value(value) {}
    
    // Function call operator
    int operator()(int x) const {
        return x + _value;
    }
};

// Usage
Adder add5(5);
int result = add5(10);  // Calls operator(), returns 15
```

#### Template Functor (Generic)

```cpp
template<typename T>
class Multiplier {
private:
    T _factor;
public:
    Multiplier(T factor) : _factor(factor) {}
    
    T operator()(T x) const {
        return x * _factor;
    }
};

// Usage
Multiplier<double> times2(2.0);
double result = times2(3.5);  // Returns 7.0
```

#### Multiple Parameters

```cpp
class Calculator {
public:
    // Overload operator() for different arities
    int operator()(int a, int b) const {
        return a + b;
    }
    
    int operator()(int a, int b, int c) const {
        return a + b + c;
    }
    
    double operator()(double a, double b) const {
        return a * b;
    }
};

// Usage
Calculator calc;
int sum = calc(1, 2);           // Calls first operator()
int sum3 = calc(1, 2, 3);       // Calls second operator()
double prod = calc(2.5, 3.0);   // Calls third operator()
```

#### Stateful Functor

```cpp
class Counter {
private:
    mutable int _count;  // mutable allows modification in const methods
public:
    Counter() : _count(0) {}
    
    int operator()() const {
        return ++_count;  // Increment and return
    }
    
    void reset() {
        _count = 0;
    }
    
    int getCount() const {
        return _count;
    }
};

// Usage
Counter counter;
int a = counter();  // Returns 1
int b = counter();  // Returns 2
int c = counter();  // Returns 3
```

#### Lambda Equivalent

```cpp
// Lambda (C++11)
auto add5 = [value = 5](int x) { return x + value; };
int result = add5(10);  // Returns 15

// Equivalent functor
class Add5 {
    int value = 5;
public:
    int operator()(int x) const { return x + value; }
};
```

#### Use Case: STL Algorithms

```cpp
class GreaterThan {
private:
    int _threshold;
public:
    GreaterThan(int threshold) : _threshold(threshold) {}
    
    bool operator()(int value) const {
        return value > _threshold;
    }
};

// Usage with STL algorithms
std::vector<int> vec = {1, 5, 10, 15, 20};
GreaterThan gt10(10);
auto it = std::find_if(vec.begin(), vec.end(), gt10);
// Finds first element > 10
```

#### Advanced: Generic Callable Wrapper

```cpp
template<typename T>
class Function {
private:
    T _func;
public:
    Function(T func) : _func(func) {}
    
    template<typename... Args>
    auto operator()(Args&&... args) -> decltype(_func(std::forward<Args>(args)...)) {
        return _func(std::forward<Args>(args)...);
    }
};

// Usage
Function<int(*)(int, int)> add([](int a, int b) { return a + b; });
int result = add(3, 4);  // Returns 7
```

**Key Points:**

- **Callable objects**: Makes objects behave like functions
- **State preservation**: Can maintain state between calls
- **Overloading**: Can have multiple `operator()` with different signatures
- **STL compatibility**: Works with algorithms expecting callables
- **Lambda alternative**: Lambdas are syntactic sugar for functors

**Interview Answer:**

> "`operator()` makes objects callable. Implement `ReturnType operator()(Params...)`. Can maintain state, overload for different signatures. Used with STL algorithms. Lambdas are syntactic sugar for functors with `operator()`."

---

### 8. `operator[]` - Subscript Operator Implementation

The subscript operator `operator[]` provides array-like access to container elements.

#### Basic Implementation

```cpp
class Array {
private:
    int* _data;
    size_t _size;
public:
    Array(size_t size) : _size(size) {
        _data = new int[_size];
    }
    
    ~Array() {
        delete[] _data;
    }
    
    // Non-const version (allows modification)
    int& operator[](size_t index) {
        return _data[index];
    }
    
    // Const version (read-only)
    const int& operator[](size_t index) const {
        return _data[index];
    }
    
    size_t size() const { return _size; }
};

// Usage
Array arr(10);
arr[0] = 42;        // Calls non-const operator[]
int val = arr[0];   // Calls const operator[]
```

#### With Bounds Checking

```cpp
class SafeArray {
private:
    int* _data;
    size_t _size;
public:
    SafeArray(size_t size) : _size(size) {
        _data = new int[_size];
    }
    
    ~SafeArray() {
        delete[] _data;
    }
    
    int& operator[](size_t index) {
        if (index >= _size) {
            throw std::out_of_range("Index out of bounds");
        }
        return _data[index];
    }
    
    const int& operator[](size_t index) const {
        if (index >= _size) {
            throw std::out_of_range("Index out of bounds");
        }
        return _data[index];
    }
};
```

#### 2D Array (Matrix)

```cpp
class Matrix {
private:
    double* _data;
    size_t _rows;
    size_t _cols;
    
    size_t index(size_t row, size_t col) const {
        return row * _cols + col;
    }
public:
    Matrix(size_t rows, size_t cols) : _rows(rows), _cols(cols) {
        _data = new double[_rows * _cols];
    }
    
    ~Matrix() {
        delete[] _data;
    }
    
    // Return proxy for row access
    class RowProxy {
    private:
        Matrix& _matrix;
        size_t _row;
    public:
        RowProxy(Matrix& matrix, size_t row) : _matrix(matrix), _row(row) {}
        
        double& operator[](size_t col) {
            return _matrix._data[_matrix.index(_row, col)];
        }
    };
    
    class ConstRowProxy {
    private:
        const Matrix& _matrix;
        size_t _row;
    public:
        ConstRowProxy(const Matrix& matrix, size_t row) 
            : _matrix(matrix), _row(row) {}
        
        const double& operator[](size_t col) const {
            return _matrix._data[_matrix.index(_row, col)];
        }
    };
    
    RowProxy operator[](size_t row) {
        return RowProxy(*this, row);
    }
    
    ConstRowProxy operator[](size_t row) const {
        return ConstRowProxy(*this, row);
    }
};

// Usage
Matrix mat(3, 4);
mat[1][2] = 3.14;        // Set element at row 1, col 2
double val = mat[1][2];  // Get element
```

#### Associative Container Style

```cpp
template<typename Key, typename Value>
class SimpleMap {
private:
    struct Pair {
        Key key;
        Value value;
    };
    std::vector<Pair> _pairs;
public:
    // Return reference (creates if doesn't exist)
    Value& operator[](const Key& key) {
        // Find existing pair
        for (auto& pair : _pairs) {
            if (pair.key == key) {
                return pair.value;
            }
        }
        // Not found, create new
        _pairs.push_back({key, Value()});
        return _pairs.back().value;
    }
    
    // Const version (read-only, throws if not found)
    const Value& operator[](const Key& key) const {
        for (const auto& pair : _pairs) {
            if (pair.key == key) {
                return pair.value;
            }
        }
        throw std::out_of_range("Key not found");
    }
};

// Usage
SimpleMap<std::string, int> map;
map["one"] = 1;        // Creates entry
int val = map["one"];  // Returns 1
```

#### String-like Container

```cpp
class MyString {
private:
    char* _data;
    size_t _length;
public:
    MyString(const char* str) {
        _length = strlen(str);
        _data = new char[_length + 1];
        strcpy(_data, str);
    }
    
    ~MyString() {
        delete[] _data;
    }
    
    // Character access
    char& operator[](size_t index) {
        return _data[index];
    }
    
    const char& operator[](size_t index) const {
        return _data[index];
    }
    
    size_t length() const { return _length; }
};

// Usage
MyString str("Hello");
str[0] = 'h';        // Modify first character
char c = str[1];     // Read second character
```

#### Key Differences: `operator[]` vs `at()`

```cpp
class Container {
    // operator[] - No bounds checking (faster)
    T& operator[](size_t index) {
        return _data[index];  // Undefined behavior if out of bounds
    }
    
    // at() - With bounds checking (safer)
    T& at(size_t index) {
        if (index >= _size) {
            throw std::out_of_range("Index out of range");
        }
        return _data[index];
    }
};
```

**Key Points:**

- **Two versions**: Non-const (modifiable) and const (read-only)
- **Return reference**: Usually returns `T&` for modification
- **No bounds checking**: `operator[]` typically doesn't check (use `at()` for safety)
- **Proxy pattern**: Can return proxy objects for multi-dimensional access
- **Performance**: `operator[]` should be O(1) for containers

**Interview Answer:**

> "`operator[]` provides array-like access. Implement `T& operator[](Index)` and `const T& operator[](Index) const`. Usually no bounds checking (use `at()` for that). Return reference for modification. Can use proxy pattern for multi-dimensional access. Should be O(1) for containers."

---

## Summary

### Key Takeaways

1. **`std::unique_ptr`**: Exclusive ownership, use by default
2. **`std::make_unique`**: Safe factory function for `unique_ptr`
3. **`std::shared_ptr`**: Shared ownership, use when needed
4. **`std::make_shared`**: Efficient factory function for `shared_ptr`
5. **`std::weak_ptr`**: Break circular references, non-owning
6. **`std::vector`**: Dynamic array with capacity management, amortized O(1) push_back
7. **`operator()`**: Makes objects callable (functors), maintains state
8. **`operator[]`**: Provides array-like access, typically O(1), no bounds checking

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
