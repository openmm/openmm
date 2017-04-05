#ifndef OPENMM_THREAD_POOL_H_
#define OPENMM_THREAD_POOL_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2017 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#define NOMINMAX
#include "windowsExport.h"
#include <functional>
#include <pthread.h>
#include <vector>

namespace OpenMM {

/**
 * A ThreadPool creates a set of worker threads that can be used to execute tasks in parallel.
 * After creating a ThreadPool, call execute() to start a task running then waitForThreads()
 * to block until all threads have finished.  You also can synchronize the threads in the middle
 * of the task by having them call syncThreads().  In this case, the parent thread should call
 * waitForThreads() an additional time; each call waits until all worker threads have reached the
 * next syncThreads(), and the final call waits until they exit from the Task's execute() method.
 * After calling waitForThreads() to block at a synchronization point, the parent thread should
 * call resumeThreads() to instruct the worker threads to resume.
 */
class OPENMM_EXPORT ThreadPool {
public:
    class Task;
    class ThreadData;
    /**
     * Create a ThreadPool.
     *
     * @param numThreads  the number of worker threads to create.  If this is 0 (the default), the
     *                    number of threads is set equal to the number of logical CPU cores available
     */
    ThreadPool(int numThreads=0);
    ~ThreadPool();
    /**
     * Get the number of worker threads in the pool.
     */
    int getNumThreads() const;
    /**
     * Execute a Task in parallel on the worker threads.
     */
    void execute(Task& task);
    /**
     * Execute a function in parallel on the worker threads.
     */
    void execute(std::function<void (ThreadPool&, int)> task);
    /**
     * This is called by the worker threads to block until all threads have reached the same point
     * and the master thread instructs them to continue by calling resumeThreads().
     */
    void syncThreads();
    /**
     * This is called by the master thread to wait until all threads have completed the Task.  Alternatively,
     * if the threads call syncThreads(), this blocks until all threads have reached the synchronization point.
     */
    void waitForThreads();
    /**
     * Instruct the threads to resume running after blocking at a synchronization point.
     */
    void resumeThreads();
private:
    bool isDeleted;
    int numThreads, waitCount;
    std::vector<pthread_t> thread;
    std::vector<ThreadData*> threadData;
    pthread_cond_t startCondition, endCondition;
    pthread_mutex_t lock;
    Task* currentTask;
    std::function<void (ThreadPool& pool, int)> currentFunction;
};

/**
 * This defines a task that can be executed in parallel by the worker threads.
 */
class OPENMM_EXPORT ThreadPool::Task {
public:
    /**
     * Execute the task on each thread.
     * 
     * @param pool         the ThreadPool being used to execute the task
     * @param threadIndex  the index of the thread invoking this method
     */
    virtual void execute(ThreadPool& pool, int threadIndex) = 0;
};

} // namespace OpenMM

#endif // OPENMM_THREAD_POOL_H_
