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

#ifndef POSIX_THREAD_HEADER
#define POSIX_THREAD_HEADER

#include <pthread.h>
#include <semaphore.h>

#include "defines.h"

class Thread
{
  private:
    pthread_t handle;
    bool cleanup;

    static void* stub (void* arg)
    { return ((Thread*)arg)->run(); }

  protected:
    virtual void* run (void) = 0;

  public:
    Thread (void)
    { this->cleanup = false; }

    virtual ~Thread (void)
    { if (this->cleanup) pthread_detach(this->handle); }

    /* Creates a new thread and runs the run() method. */
    void pt_create (const pthread_attr_t* p = 0)
    {
      this->cleanup = true;
      pthread_create(&this->handle, p, Thread::stub, (void*)this);
    }

    /* Sends a cancelation request to the thread. */
    void pt_cancel (void)
    { pthread_cancel(this->handle); }

    /* Blocks and waits for termination of the thread. */
    void pt_join (void)
    {
      this->cleanup = false;
      pthread_join(this->handle, 0);
    }
};

/* ---------------------------------------------------------------- */

class Mutex
{
  private:
    pthread_mutex_t mutex;

  private:
    /* Don't allow copying of the mutex. */
    Mutex (Mutex const& UNUSED(rhs))
    { }

    /* Don't allow copying of the mutex. */
    Mutex& operator= (Mutex const& UNUSED(rhs))
    { return *this; }

  public:
    Mutex (void)
    { pthread_mutex_init(&this->mutex, 0); }

    ~Mutex (void)
    { pthread_mutex_destroy(&this->mutex); }

    /* Acquire ownership of the mutex. This will block if the mutex is
     * locked by another process until the mutex is released. */
    int lock (void)
    { return pthread_mutex_lock(&this->mutex); }

    /* This call acquires ownership of the mutex. If the mutex is
     * already locked, this call will NOT block and return with EBUSY. */
    int trylock (void)
    { return pthread_mutex_trylock(&this->mutex); }

    /* Release ownership of the mutex, and unlocks the mutex. */
    int unlock (void)
    { return pthread_mutex_unlock(&this->mutex); }
};

/* ---------------------------------------------------------------- */

class Semaphore
{
  private:
    sem_t sem;

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
    { sem_init(&sem, 0, value); }

    Semaphore (unsigned int value, int pshared)
    { sem_init(&sem, pshared, value); }

    ~Semaphore (void)
    { sem_destroy(&sem); }

    /* DOWN the semaphore. */
    int wait (void)
    { return sem_wait(&sem); }

    /* UP the semaphore. */
    int post (void)
    { return sem_post(&sem); }

    /* Returns the current value. */
    int get_value (void)
    {
      int value;
      sem_getvalue(&sem, &value);
      return value;
    }
};

/* ---------------------------------------------------------------- */

/*
 * A read write lock using a shared read lock and an exclusive write lock.
 * Several readers can access critical data but writing to the data
 * is exclusive, i.e. no other writer or reader is allowed during write.
 */
class ReadWriteLock
{
  private:
    pthread_rwlock_t rwlock;

  private:
    /* Don't allow copying of the lock. */
    ReadWriteLock (ReadWriteLock const& UNUSED(rhs))
    { }

    /* Don't allow copying of the lock. */
    ReadWriteLock& operator= (ReadWriteLock const& UNUSED(rhs))
    { return *this; }

  public:
    ReadWriteLock (void)
    { pthread_rwlock_init(&this->rwlock, 0); }

    ~ReadWriteLock (void)
    { pthread_rwlock_destroy(&this->rwlock); }

    int read_lock (void)
    { return pthread_rwlock_rdlock(&this->rwlock); }

    int write_lock (void)
    { return pthread_rwlock_wrlock(&this->rwlock); }

    int read_trylock (void)
    { return pthread_rwlock_tryrdlock(&this->rwlock); }

    int write_trylock (void)
    { return pthread_rwlock_trywrlock(&this->rwlock); }

    int unlock (void)
    { return pthread_rwlock_unlock(&this->rwlock); }
};

#endif /* POSIX_THREAD_HEADER */
