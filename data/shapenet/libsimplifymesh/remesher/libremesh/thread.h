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

#ifndef THREAD_HEADER
#define THREAD_HEADER

#ifndef _WIN32
#   include <unistd.h>
#endif

#ifdef _WIN32
#   include "thread_win32.h"
#elif defined(_POSIX_THREADS)
#   include "thread_posix.h"
#else
#   error "Thread abstraction not available. PLEASE REPORT THIS!"
#endif

/* ---------------------------------------------------------------- */

/*
 * Exception safe lock for the mutex.
 */
class MutexLock
{
  private:
    Mutex* m;

  public:
    MutexLock (Mutex& mutex)
    { this->m = &mutex; this->m->lock(); }

    ~MutexLock (void)
    { if (this->m) this->m->unlock(); }

    void unlock (void)
    { this->m->unlock(); this->m = 0; }
};

/* ---------------------------------------------------------------- */

/*
 * Exception safe read lock for the reader/writer lock.
 */
class ReadLock
{
  private:
    ReadWriteLock* rwl;

  public:
    ReadLock (ReadWriteLock& rwlock)
    { this->rwl = &rwlock; this->rwl->read_lock(); }

    ~ReadLock (void)
    { if (this->rwl) this->rwl->unlock(); }

    void unlock (void)
    { this->rwl->unlock(); this->rwl = 0; }
};

/* ---------------------------------------------------------------- */

/*
 * Exception safe write lock for the reader/writer lock.
 */
class WriteLock
{
  private:
    ReadWriteLock* rwl;

  public:
    WriteLock (ReadWriteLock& rwlock)
    { this->rwl = &rwlock; this->rwl->write_lock(); }

    ~WriteLock (void)
    { if (this->rwl) this->rwl->unlock(); }

    void unlock (void)
    { this->rwl->unlock(); this->rwl = 0; }
};

#endif /* THREAD_HEADER */
