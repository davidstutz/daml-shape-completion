/*
 * This file is part of GtkEveMon.
 *
 * GtkEveMon is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Public License
 * along with GtkEveMon. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef WIN32_THREAD_HEADER
#define WIN32_THREAD_HEADER

#include <iostream>
#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>

#include "defines.h"

class Thread
{
  private:
    HANDLE handle;
    HANDLE cancel_event;
    bool cleanup;

    static void* stub (void* arg)
    { return ((Thread*)arg)->run(); }

  protected:
    virtual void* run (void) = 0;

  public:
    Thread (void)
    { this->cleanup = false; }

    virtual ~Thread (void)
    {
      if (this->cleanup)
      {
        CloseHandle(this->cancel_event);
        CloseHandle(this->handle);
      }
    }

    /* Creates a new thread and runs the run() method. */
    // FIXME: Abstract pthread_attr_t
    void pt_create (void)
    {
      this->cleanup = true;
      this->cancel_event = CreateEvent(NULL, FALSE, FALSE, NULL);
      this->handle = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)
          Thread::stub, (void*)this, 0, NULL);

      if (this->handle == NULL)
        std::cerr << "ERROR: CreateThread(): " << GetLastError() << std::endl;
    }

    /* Sends a cancelation request to the thread. */
    void pt_cancel (void)
    { SetEvent(this->cancel_event); }

    /* Blocks and waits for termination of the thread. */
    void pt_join (void)
    {
      this->cleanup = false;
      WaitForSingleObject(this->handle, INFINITE);
      CloseHandle(this->cancel_event);
      CloseHandle(this->handle);
    }
};

/* ---------------------------------------------------------------- */
// TODO Implement for Win32

class Mutex
{
  private:
    //pthread_mutex_t mutex;

  private:
    /* Don't allow copying of the mutex. */
    Mutex (Mutex const& UNUSED(rhs))
    { }

    /* Don't allow copying of the mutex. */
    Mutex& operator= (Mutex const& UNUSED(rhs))
    { return *this; }

  public:
    Mutex (void)
    { }
    //{ pthread_mutex_init(&this->mutex, 0); }

    ~Mutex (void)
    { }
    //{ pthread_mutex_destroy(&this->mutex); }

    /* Acquire ownership of the mutex. This will block if the mutex is
     * locked by another process until the mutex is released. */
    int lock (void)
    { return 0; }
    //{ return pthread_mutex_lock(&this->mutex); }

    /* This call acquires ownership of the mutex. If the mutex is
     * already locked, this call will NOT block and return with EBUSY. */
    int trylock (void)
    { return 0; }
    //{ return pthread_mutex_trylock(&this->mutex); }

    /* Release ownership of the mutex, and unlocks the mutex. */
    int unlock (void)
    { return 0; }
    //{ return pthread_mutex_unlock(&this->mutex); }
};

/* ---------------------------------------------------------------- */

class Semaphore
{
  private:
    HANDLE sem;
    volatile LONG value;

  private:
    /* Don't allow copying of semaphores. */
    Semaphore (Semaphore const& UNUSED(rhs))
    { }

    /* Don't allow copying of semaphores. */
    Semaphore& operator= (Semaphore const& UNUSED(rhs))
    { return *this; }

  public:
    /* Defaults to mutex. */
    Semaphore (unsigned int value = 1)
    {
      this->value = value;
      this->sem = CreateSemaphore(NULL, value, value, NULL);
    }

    ~Semaphore (void)
    {
      // TODO Destroy semaphore
    }

    /* DOWN the semaphore. */
    int wait (void)
    {
      int retval = WaitForSingleObject(this->sem, INFINITE);
      if (retval == WAIT_OBJECT_0)
      {
        InterlockedDecrement(&this->value);
        return 0;
      }
      return -1;
    }

    /* UP the semaphore. */
    int post (void)
    {
      InterlockedIncrement(&this->value);
      if (ReleaseSemaphore(this->sem, 1, NULL))
      {
        InterlockedDecrement(&this->value);
        return -1;
      }
      return 0;
    }

    /* Returns the current value. */
    int get_value (void)
    {
      return this->value;
    }
};

/* ---------------------------------------------------------------- */
// TODO Implement for Win32

/*
 * A read write lock using a shared read lock and an exclusive write lock.
 * Several readers can access critical data but writing to the data
 * is exclusive, i.e. no other writer or reader is allowed during write.
 */
class ReadWriteLock
{
  private:
    //pthread_rwlock_t rwlock;

  private:
    /* Don't allow copying of the lock. */
    ReadWriteLock (ReadWriteLock const& UNUSED(rhs))
    { }

    /* Don't allow copying of the lock. */
    ReadWriteLock& operator= (ReadWriteLock const& UNUSED(rhs))
    { return *this; }

  public:
    ReadWriteLock (void)
    { }
    //{ pthread_rwlock_init(&this->rwlock, 0); }

    ~ReadWriteLock (void)
    { }
    //{ pthread_rwlock_destroy(&this->rwlock); }

    int read_lock (void)
    { return 0; }
    //{ return pthread_rwlock_rdlock(&this->rwlock); }

    int write_lock (void)
    { return 0; }
    //{ return pthread_rwlock_wrlock(&this->rwlock); }

    int read_trylock (void)
    { return 0; }
    //{ return pthread_rwlock_tryrdlock(&this->rwlock); }

    int write_trylock (void)
    { return 0; }
    //{ return pthread_rwlock_trywrlock(&this->rwlock); }

    int unlock (void)
    { return 0; }
    //{ return pthread_rwlock_unlock(&this->rwlock); }
};


#endif /* WIN32_THREAD_HEADER */
