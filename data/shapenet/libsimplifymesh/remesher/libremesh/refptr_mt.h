/*
 *  A reference counting smart pointer for C++
 *  Copyright (c) 2005-2006 by Simon Fuhrmann
 *
 *  This code is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; version 2 dated June, 1991.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the application. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef REF_PTR_HEADER
#define REF_PTR_HEADER

#include <algorithm>
#include "thread.h"

template <class T>
class RefPtr
{
  /* Private declaration of class members. */
  protected:
    T* ptr;
    int* count;
    Mutex* mutex;

  /* Private definition of member methods. */
  private:
    void increment (void)
    {
      if (ptr == 0) return;
      MutexLock lock(*mutex);
      ++(*count);
    }

    void decrement (void)
    {
      if (ptr == 0) return;
      MutexLock lock(*mutex);
      if (--(*count) == 0)
      {
        lock.unlock();
        delete count;
        delete ptr;
        delete mutex;
      }
    }

  /* Public definition of member methods. */
  public:
    /* Ctor: Default one. */
    RefPtr (void) : ptr(0), count(0), mutex(0)
    { }

    /* Ctor: From pointer. */
    explicit RefPtr (T* p) : ptr(p)
    {
      count = (p == 0) ? 0 : new int(1);
      mutex = (p == 0) ? 0 : new Mutex;
    }

    /* Ctor: Copy from other RefPtr. */
    RefPtr (RefPtr<T> const& src) : ptr(src.ptr),
        count(src.count), mutex(src.mutex)
    { increment(); }

    /* Destructor. */
    ~RefPtr (void)
    { decrement(); }

    /* Assignment: From other RefPtr. */
    RefPtr<T>& operator= (RefPtr<T> const& rhs)
    {
      if (rhs.ptr == ptr) return *this;
      decrement();
      ptr = rhs.ptr;
      count = rhs.count;
      mutex = rhs.mutex;
      increment();
      return *this;
    }

    /* Assignment: From pointer. */
    RefPtr<T>& operator= (T* rhs)
    {
      if (rhs == ptr) return *this;
      decrement();
      ptr = rhs;
      count = (ptr == 0) ? 0 : new int(1);
      mutex = (ptr == 0) ? 0 : new Mutex;
      return *this;
    }

    /* Operations. */
    void reset (void)
    {
      decrement();
      ptr = 0; count = 0; mutex = 0;
    }

    void swap (RefPtr<T>& p)
    {
      std::swap(p.ptr, ptr);
      std::swap(p.count, count);
      std::swap(p.mutex, mutex);
    }

    /* Dereference. */
    T& operator* (void) const
    { return *ptr; }

    T* operator-> (void) const
    { return ptr; }

    /* Comparison. */
    bool operator== (T const* p) const
    { return p == ptr; }

    bool operator== (RefPtr<T> const& p) const
    { return p.ptr == ptr; }

    bool operator!= (T const* p) const
    { return p != ptr; }

    bool operator!= (RefPtr<T> const& p) const
    { return p.ptr != ptr; }

    bool operator< (RefPtr<T> const& rhs)
    { return ptr < rhs.ptr; }

    /* Information. */
    int use_count (void) const
    { return (count == 0) ? 0 : *count; }

    T* get (void) const
    { return ptr; }

  #ifndef NO_MEMBER_TEMPLATES

  /* Template friends for accessing diffrent RefPtr's. */
  private:
    template <class Y> friend class RefPtr;

  public:
    /* Ctor: From diffrent pointer. */
    template <class Y>
    explicit RefPtr (Y* p) : ptr(static_cast<T*>(p))
    {
      count = (p == 0) ? 0 : new int(1);
      mutex = (p == 0) ? 0 : new Mutex;
    }

    /* Ctor: Copy from diffrent RefPtr. */
    template <class Y>
    RefPtr (RefPtr<Y> const& src)
    {
      ptr = static_cast<T*>(src.ptr);
      count = src.count;
      mutex = src.mutex;
      increment();
    }

    /* Assignment: From diffrent RefPtr. */
    template <class Y>
    RefPtr<T>& operator= (RefPtr<Y> const& rhs)
    {
      if (rhs.ptr == ptr) return *this;
      decrement();
      ptr = static_cast<T*>(rhs.ptr);
      count = rhs.count;
      mutex = rhs.mutex;
      increment();
      return *this;
    }

    /* Assignment: From diffrent pointer. */
    template <class Y>
    RefPtr<T>& operator= (Y* rhs)
    {
      if (rhs == ptr) return *this;
      decrement();
      ptr = static_cast<T*>(rhs);
      count = (ptr == 0) ? 0 : new int(1);
      mutex = (ptr == 0) ? 0 : new Mutex;
      return *this;
    }

    /* Comparison with diffrent RefPtr type. */
    template <class Y>
    bool operator== (RefPtr<Y> const& p) const
    { return p.ptr == ptr; }

    template <class Y>
    bool operator!= (RefPtr<Y> const& p) const
    { return p.ptr != ptr; }

    template <class Y>
    bool operator< (RefPtr<Y> const& rhs) const
    { return ptr < rhs.ptr; }

  #endif /* NO_MEMBER_TEMPLATES */
};

#endif /* REF_PTR_HEADER */
