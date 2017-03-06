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

#include "openmm/internal/ThreadPool.h"
#include "openmm/internal/hardware.h"

using namespace std;

namespace OpenMM {

class ThreadPool::ThreadData {
public:
    ThreadData(ThreadPool& owner, int index) : owner(owner), index(index), isDeleted(false) {
    }
    void executeTask() {
        if (owner.currentTask != NULL)
            owner.currentTask->execute(owner, index);
        else
            owner.currentFunction(owner, index);
    }
    ThreadPool& owner;
    int index;
    bool isDeleted;
    Task* currentTask;
    function<void (ThreadPool& pool, int)> currentFunction;
};

static void* threadBody(void* args) {
    ThreadPool::ThreadData& data = *reinterpret_cast<ThreadPool::ThreadData*>(args);
    while (true) {
        // Wait for the signal to start running.
        
        data.owner.syncThreads();
        if (data.isDeleted)
            break;
        data.executeTask();
    }
    delete &data;
    return 0;
}

ThreadPool::ThreadPool(int numThreads) : currentTask(NULL) {
    if (numThreads <= 0)
        numThreads = getNumProcessors();
    this->numThreads = numThreads;
    pthread_cond_init(&startCondition, NULL);
    pthread_cond_init(&endCondition, NULL);
    pthread_mutex_init(&lock, NULL);
    thread.resize(numThreads);
    pthread_mutex_lock(&lock);
    waitCount = 0;
    for (int i = 0; i < numThreads; i++) {
        ThreadData* data = new ThreadData(*this, i);
        data->isDeleted = false;
        threadData.push_back(data);
        pthread_create(&thread[i], NULL, threadBody, data);
    }
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
}

ThreadPool::~ThreadPool() {
    for (auto data : threadData)
        data->isDeleted = true;
    pthread_mutex_lock(&lock);
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
    for (auto t : thread)
        pthread_join(t, NULL);
    pthread_mutex_destroy(&lock);
    pthread_cond_destroy(&startCondition);
    pthread_cond_destroy(&endCondition);
}

int ThreadPool::getNumThreads() const {
    return numThreads;
}

void ThreadPool::execute(Task& task) {
    currentTask = &task;
    resumeThreads();
}

void ThreadPool::execute(function<void (ThreadPool&, int)> task) {
    currentTask = NULL;
    currentFunction = task;
    resumeThreads();
}

void ThreadPool::syncThreads() {
    pthread_mutex_lock(&lock);
    waitCount++;
    pthread_cond_signal(&endCondition);
    pthread_cond_wait(&startCondition, &lock);
    pthread_mutex_unlock(&lock);
}

void ThreadPool::waitForThreads() {
    pthread_mutex_lock(&lock);
    while (waitCount < numThreads)
        pthread_cond_wait(&endCondition, &lock);
    pthread_mutex_unlock(&lock);
}

void ThreadPool::resumeThreads() {
    pthread_mutex_lock(&lock);
    waitCount = 0;
    pthread_cond_broadcast(&startCondition);
    pthread_mutex_unlock(&lock);
}

} // namespace OpenMM
