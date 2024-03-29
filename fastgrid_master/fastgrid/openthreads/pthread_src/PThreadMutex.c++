//
// OpenThread library, Copyright (C) 2002 - 2003  The Open Thread Group
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//

//
// PThreadMutex.c++ - C++ Mutex class built on top of posix threads.
// ~~~~~~~~~~~~~~~~
//

#include <unistd.h>
#include <pthread.h>
#include "../Mutex"
#include "PThreadMutexPrivateData.h"

using namespace OpenThreads;

//----------------------------------------------------------------------------
//
// Decription: Constructor
//
// Use: public.
//
Mutex::Mutex() {

    pthread_mutexattr_t mutex_attr;
    pthread_mutexattr_init( &mutex_attr );
    
    PThreadMutexPrivateData *pd = new PThreadMutexPrivateData();
	
#ifndef __linux__ // (not available until NPTL) [
    pthread_mutexattr_settype( &mutex_attr, PTHREAD_MUTEX_ERRORCHECK );
#endif // ] __linux__

#ifdef ALLOW_PRIORITY_SCHEDULING // [

#ifdef __sun // [
    pthread_mutexattr_setprotocol(&mutex_attr, PTHREAD_PRIO_NONE);
#endif // ] __sun

    //-------------------------------------------------------------------------
    // Initialization is a bit tricky, since we have to be able to be aware
    // that on many-to-many execution vehicle systems, we may run into
    // priority inversion deadlocks if a mutex is shared between threads
    // of differing priorities.  Systems that do this should provide the 
    // following protocol attributes to prevent deadlocks.  Check at runtime.
    //
    //  PRIO_INHERIT causes any thread locking the mutex to temporarily become
    //  the same priority as the highest thread also blocked on the mutex. 
    //  Although more expensive, this is the prefered method.
    //
    //  PRIO_PROTECT causes any thread locking the mutex to assume the priority
    //  specified by setprioceiling.  pthread_mutex_lock will fail if
    //  the priority ceiling is lower than the thread's priority.  Therefore,
    //  the priority ceiling must be set to the max priority in order to 
    //  garantee no deadlocks will occur.
    //
#if defined (_POSIX_THREAD_PRIO_INHERIT) || defined (_POSIX_THREAD_PRIO_PROTECT) // [

    if(sysconf(_POSIX_THREAD_PRIO_INHERIT)) {
	
	pthread_mutexattr_setprotocol(&mutex_attr, PTHREAD_PRIO_INHERIT);
	
    } else if (sysconf(_POSIX_THREAD_PRIO_PROTECT)) {

	int th_policy;
	struct sched_param th_param;
	pthread_getschedparam(pthread_self(), &th_policy, &th_param);

	pthread_mutexattr_setprotocol(&mutex_attr, PTHREAD_PRIO_PROTECT);
	pthread_mutexattr_setprioceiling(&mutex_attr, 
					 sched_get_priority_max(th_policy));

    }

#endif // ] Priority Scheduling.

#endif // ] ALLOW_PRIORITY_SCHEDULING

    pthread_mutex_init(&pd->mutex, &mutex_attr);
    _prvData = static_cast<void *>(pd);

}

//----------------------------------------------------------------------------
//
// Decription: Destructor
//
// Use: public.
//
Mutex::~Mutex() {

    PThreadMutexPrivateData *pd =
        static_cast<PThreadMutexPrivateData*>(_prvData);

    pthread_mutex_destroy(&pd->mutex);

    delete pd;

}

//----------------------------------------------------------------------------
//
// Decription: lock the mutex
//
// Use: public.
//
int Mutex::lock() {

    PThreadMutexPrivateData *pd =
        static_cast<PThreadMutexPrivateData*>(_prvData);

    return pthread_mutex_lock(&pd->mutex);

}

//----------------------------------------------------------------------------
//
// Decription: unlock the mutex
//
// Use: public.
//
int Mutex::unlock() {

    PThreadMutexPrivateData *pd =
        static_cast<PThreadMutexPrivateData*>(_prvData);

    return pthread_mutex_unlock(&pd->mutex);

}

//----------------------------------------------------------------------------
//
// Decription: test if the mutex may be locked
//
// Use: public.
//
int Mutex::trylock() {

    PThreadMutexPrivateData *pd =
        static_cast<PThreadMutexPrivateData*>(_prvData);

    return pthread_mutex_trylock(&pd->mutex);

}
