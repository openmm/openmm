/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2013-2025 Stanford University and the Authors.      *
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
    unique_lock<mutex> ul(lock);
    waitCount = 0;
    for (int i = 0; i < numThreads; i++) {
        ThreadData* data = new ThreadData(*this, i);
        data->isDeleted = false;
        threadData.push_back(data);
        threads.push_back(thread(threadBody, data));
    }
    while (waitCount < numThreads)
        endCondition.wait(ul);
}

ThreadPool::~ThreadPool() {
    for (auto data : threadData)
        data->isDeleted = true;
    lock.lock();
    startCondition.notify_all();
    lock.unlock();
    for (auto &t : threads)
        t.join();
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
    unique_lock<mutex> ul(lock);
    waitCount++;
    endCondition.notify_one();
    startCondition.wait(ul);
}

void ThreadPool::waitForThreads() {
    unique_lock<mutex> ul(lock);
    while (waitCount < numThreads)
        endCondition.wait(ul);
}

void ThreadPool::resumeThreads() {
    unique_lock<mutex> ul(lock);
    waitCount = 0;
    startCondition.notify_all();
}

} // namespace OpenMM
