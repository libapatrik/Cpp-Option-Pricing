// Smart Pointers

// unique_ptr 
// Wrapper around a raw pointer
#include <algorithm>
template<class T>
class unique_ptr
{
public:
	unique_ptr(T* rawPtr)
		: _rawPtr(new T(*rawPtr))
	{ }

	~unique_ptr()
	{
		delete _rawPtr;
	}

private:
	T* _rawPtr;
};



// ============================================================================
// Implement move operator
// ============================================================================
/**
 *  move operator
 *  Addresses the core problem of expensive copying
 *  The unique_ptr cannot be copied, because it would create another owner -> hence not unique anymore, but can be moved
 *  std::move(x) says "treat x as if it's a temporary, so it can be moved from."
 *  MOVING (cheap, transfers ownership)
 *  std::unique_ptr<House> myHouse = std::make_unique<House>();
 *  std::unique_ptr<House> yourHouse = std::move(myHouse);  
 *  Now: yourHouse owns the house, myHouse is empty (nullptr)

 * The move operator transfers ownership of resources without copying. 
 * The move constructor steals the pointer and sets the source to nullptr. 
 * The move assignment also handles self-assignment and cleans up existing resources first. Both are marked noexcept for performance and safety.

 * KEY: 
 * Q: "What's the difference from copy?"
 * A: "Copy duplicates resources (expensive). Move transfers ownership (cheap). Copy leaves both objects usable. Move leaves source empty."
 */
 class move_operator 
 {
public: 
	// constructor
	move_operator(move_operator&& other) noexcept // && rvalue reference, 
		: _rawPtr(std::move(other._rawPtr))
 	{
	}
	// Why noexcept? I declare and define a method that will not throw error.
	// Enables optimizations. Containers like std::vector will use move instead of copy when possible, and won't fall back to copy if move throws.

	// move assignment operator
	move_operator& operator=(move_operator&& other) noexcept   // & lvalue and && rvalue
	{
		if (this != &other) {         // self-assignment check; i.e. do not move to the same house again
			delete _rawPtr;           // clean the house
			_rawPtr = std::move(other._rawPtr);  // move the house i.e. transfer/steal ownership (no duplications)
			//other._rawPtr = nullptr;  // leave the other house empty; i.e. the source is empty but valid (i.e. valid pointer, but no ownership)
		}
		return *this;
	}

	// deletes
	move_operator(const move_operator& other) = delete;
	move_operator& operator=(const move_operator& other) = delete;
 };

// Learning about operators

// double operator()(double x) const // to handle routing to extrapolate() when needed
// {
//     auto [xMin, xMax] = getRange();
//     if (x < xMin || x > xMax) {
//         return _extrapolationScheme->extrapolate(x);
//     } else {
//         return interpolate(x);
//     }
// }
// // Class instance(param1, param2,...) // Constructor
// // LinearInterpolationScheme interp(params); // Constructor
// // interp(x) -> interp.operator()(x);
// operator[](size_t row, size_t column)
// Matrix A;
// A[2,3]

// operator=(double x)
// double a;
// double x;
// a = x -> a.operator=(x)