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
// PThreadBarrierPrivateData.h - private data structure for barrier
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#ifndef _PTHREADBARRIERPRIVATEDATA_H_
#define _PTHREADBARRIERPRIVATEDATA_H_

namespace OpenThreads {

#include <pthread.h>
#include "../Barrier"

class PThreadBarrierPrivateData {

    friend class Barrier;

private:

    PThreadBarrierPrivateData() {};
    
    virtual ~PThreadBarrierPrivateData() {};

    pthread_cond_t     cond;            // cv for waiters at barrier

    pthread_mutex_t    lock;            // mutex for waiters at barrier

    volatile int       maxcnt;          // number of threads to wait for

    volatile int       cnt;             // number of waiting threads

    volatile int       phase;           // flag to seperate two barriers

};

}

#endif //_PTHREADBARRIERPRIVATEDATA_H_
