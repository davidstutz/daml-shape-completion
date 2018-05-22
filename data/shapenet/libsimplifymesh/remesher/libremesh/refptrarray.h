#ifndef REF_PTR_ARRAY_HEADER
#define REF_PTR_ARRAY_HEADER

#include "refptr.h"

template <class T, class E>
class RefPtrArray : public RefPtr<T>
{
  public:
    /* Ctor: Default one. */
    RefPtrArray (void)
    { }

    /* Ctor: From pointer. */
    explicit RefPtrArray (T* p) : RefPtr<T>(p)
    { }

    /* Ctor: Copy from other RefPtr. */
    RefPtrArray (RefPtrArray<T, E> const& src) : RefPtr<T>(src)
    { }

    /* Element access const operator. */
    E const& operator[] (std::size_t index) const
    { return (*this->ptr)[index]; }

    /* Element access operator. */
    E& operator[] (std::size_t index)
    { return (*this->ptr)[index]; }
};

#endif /* REF_PTR_ARRAY_HEADER */
