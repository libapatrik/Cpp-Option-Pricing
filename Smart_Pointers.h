// Smart Pointers

// unique_ptr 
// Wrapper around a raw pointer

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




