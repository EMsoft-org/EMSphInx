/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019, William C. Lenthe                               *
 * All rights reserved.                                                *
 *                                                                     *
 * This package is free software; you can redistribute it and/or       *
 * modify it under the terms of the GNU General Public License as      *
 * published by the Free Software Foundation; either version 2 of the  *
 * License, or (at your option) any later version.                     *
 *                                                                     *
 * This program is distributed in the hope that it will be useful,     *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of      *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       *
 * GNU General Public License for more details.                        *
 *                                                                     *
 * You should have received a copy of the GNU General Public License   *
 * along with this program; if not, check the Free Software Foundation *
 * website: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>   *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _THREADPOOL_H_
#define _THREADPOOL_H_

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <chrono>//duration

class ThreadPool {
	private:
		//helper class to hold jobs with concurrent access
		struct TaskQueue {
			private:
				std::queue< std::function<void(const size_t)> > tasks;//queue of work items (argument-less / return-less functions to be executed)
				unsigned int                                    count;//counter for number of incomplete tasks (including tasks have been popped but not run)
				mutable std::condition_variable                 cCv  ;//condition variable for queue completing
				mutable std::mutex                              qMut ;//mutex to handle concurrent access to tasks
				mutable std::mutex                              iMut ;//mutex for condition variable and access to count

			public:
				TaskQueue() : count(0) {}//set number in queue to 0 at construction
				~TaskQueue() {clear();}//clear remaining tasks on destruction

				//@brief: add a work item to the queue
				//@param task: item to add to the queue
				void push(const std::function<void(const size_t)>& task);

				//@brief: do some work from the queue
				//@param i: index of thread doing work
				//@return: true/false if the queue was/wasn't empty (work wasn't/was done)
				bool process(const size_t i);

				//@brief: remove all unprocessed tasks from the queue (task that are already started won't be stopped)
				void clear();

				//@brief: check if the queue is empty
				//@return: true/false if the queue is/isn't empty 
				//@note: only check if there are pending tasks (that haven't been started), not if there are incomplete tasks running
				bool empty() const;

				//@brief: get the size of the queue
				//@return: size of the queue
				//@note: only counts pending tasks (that haven't been started), not incomplete tasks running
				unsigned int size() const;

				//@brief: wait for the queue to be empty and all running tasks to finish
				//@param timeout: duration to wait for queue to complete
				//@return: true/false if the queue was/wasn't complete before the timeout
				template <class Rep, class Period = std::ratio<1> >
				bool waitComplete(std::chrono::duration<Rep, Period> timeout) const;
		};

		std::vector<std::thread>        pool;//pool of worker threads
		bool                            live;//flag for if threads should continue looking for work
		TaskQueue                       tQue;//unprocessed work items
		std::condition_variable         wCv ;//condition variable to notify idle workers in pool
		std::mutex                      wMut;//mutex for wCv and write access to live

	public:
		//make noncopyable
		ThreadPool           (const ThreadPool &) = delete;
		ThreadPool& operator=(const ThreadPool &) = delete;
		
		//@brief: determine the number of threads in a default pool
		//@return: number of threads recommended by the system
		static size_t Concurrency() {return std::max<size_t>(1, std::thread::hardware_concurrency());}

		//@brief: construct a thread pool
		//@param threadCount: number of worker threads in pool (optional)
		ThreadPool(const size_t threadCount = Concurrency());

		//@brief: clean up thread pool (waits for all tasks to finish as written)
		~ThreadPool();

		//@brief: add a task to the work queue
		//@param task: work item to add to queue
		void schedule(const std::function<void(const size_t)>& task);

		//@brief: get number of worker threads in pool
		//@return: number of threads
		size_t size() const {return pool.size();}

		//@brief: get number of items in the queue
		//@return: number of items
		unsigned int items() const {return tQue.size();}

		//@brief: clear out any unfinished (and unstarted) work items in the queue
		void clear() {tQue.clear();}

		//@brief: blocking wait for all tasks in the queue to finish
		//@param timeout: duration to wait for tasks to complete
		//@return: true/false if the the tasks were / weren't completed before timeout
		template <class Rep, class Period = std::ratio<1> >
		bool waitAll(std::chrono::duration<Rep, Period> timeout) const {return tQue.waitComplete(timeout);}

		//@brief: blocking wait for all tasks in the queue to finish
		void waitAll() const {while(!waitAll(std::chrono::hours(1))){};}//wait indefinitely in 1 hour increments
};

/////////////////////////////////
// ThreadPool member functions //
/////////////////////////////////

#include <utility>//move

//@brief: construct a thread pool
//@param threadCount: number of worker threads in pool (optional)
ThreadPool::ThreadPool(const size_t threadCount) : live(true) {
	std::function<void(const size_t)> func = [this](const size_t i)->void{//every thread will run the same function:
		std::function<void(const size_t)> workItem;//the thread does work of type void()
		while(live) {//keep looking for work as long as the pool is live
			if(tQue.process(i)) {//try to do some work and check if the queue is empty (no work is available)
				std::unique_lock<std::mutex> lock(wMut);//lock condition variable mutex
				if(live && tQue.empty()) wCv.wait(lock);//wait for a work item to arrive (if one hasn't arrived meanwhile)
			}
		}
	};
	for(unsigned int i = 0; i < threadCount; i++) pool.emplace_back(std::bind(func,i));//loop over threads constructing
}

//@brief: clean up thread pool (waits for all tasks to finish as written)
ThreadPool::~ThreadPool() {
	waitAll();//wait for workers to finish waiting and in progress tasks (could call tQue.clear() before to remove remaining tasks)
	wMut.lock();//lock mutex so threads currently checking for work don't start sleeping before live is set to false
	live = false;//tell the workers to stop looking for new work
	wCv.notify_all();//wake up any workers waiting for a work item
	wMut.unlock();//release mutex so sleeping threads can wake up
	for(std::thread& t: pool) if(t.joinable()) t.join();//clean up workers
}

//@brief: add a task to the work queue
//@param task: work item to add to queue
void ThreadPool::schedule(const std::function<void(const size_t)>& task) {
	std::lock_guard<std::mutex> lock(wMut);//lock mutex so wCv isn't notified before a thread that just failed to get work goes to sleep
	tQue.push(task);//add task to queue
	wCv.notify_one();//notify a sleeping thread to wake up if needed
}

/////////////////////////////////
// TaskQueue member functions  //
/////////////////////////////////

//@brief: add a work item to the queue
//@param task: item to add to the queue
void ThreadPool::TaskQueue::push(const std::function<void(const size_t)>& task) {
	std::unique_lock<std::mutex> tLock(qMut, std::defer_lock);//build lock around queue mutex but don't acquire yet
	std::unique_lock<std::mutex> iLock(iMut, std::defer_lock);//build lock around count mutex but don't acquire yet
	std::lock(tLock, iLock);//simultaneously lock qMut and iMut (TODO: replace with std::scoped_lock(qMut, iMut) in c++17)
	tasks.push(task);//add task to queue
	++count;//increment count of items in queue
}

//@brief: do some work from the queue
//@param i: index of thread doing work
//@return: true/false if the queue was/wasn't empty (work wasn't/was done)
bool ThreadPool::TaskQueue::process(const size_t i) {
	std::function<void(const size_t)> task;//potential work to be done
	{//scope for lock_guard
		std::lock_guard<std::mutex> tLock(qMut);//lock queue mutex
		if(tasks.empty()) {//another thread got to the last item in the queue first
			return true;//return true if the queue was empty (work wasn't done)
		} else {//this thread was able to get a work item
			task = std::move(tasks.front());//get the work item
			tasks.pop();//remove the work item from the queue
		}
	}//free lock on queue before starting work
	task(i);//do the work item
	std::lock_guard<std::mutex> iLock(iMut);//lock counter mutex
	if(0 == --count) cCv.notify_all();//wake any threads waiting for the queue to be complete if needed
	return false;//return false if the queue wasn't empty (work was done)
}

//@brief: remove all unprocessed tasks from the queue (task that are already started won't be stopped)
void ThreadPool::TaskQueue::clear() {
	std::unique_lock<std::mutex> tLock(qMut, std::defer_lock);//build lock around queue mutex but don't acquire yet
	std::unique_lock<std::mutex> iLock(iMut, std::defer_lock);//build lock around count mutex but don't acquire yet
	std::lock(tLock, iLock);//simultaneously lock qMut and iMut (TODO: replace with std::scoped_lock(qMut, iMut) in c++17)
	while(!tasks.empty()) {
		tasks.pop();//remove next task
		--count;//decrement count of items in queue
	}
	if(0 == count) cCv.notify_all();//wake any threads waiting for the queue to be complete if needed
}

//@brief: check if the queue is empty
//@return: true/false if the queue is/isn't empty 
//@note: only check if there are pending tasks (that haven't been started), not if there are incomplete tasks running
bool ThreadPool::TaskQueue::empty() const {
	std::lock_guard<std::mutex> tLock(qMut);//lock queue mutex
	return tasks.empty();//return if there are pending tasks
}

//@brief: get the size of the queue
//@return: size of the queue
//@note: only counts pending tasks (that haven't been started), not incomplete tasks running
unsigned int ThreadPool::TaskQueue::size() const {
	std::lock_guard<std::mutex> tLock(iMut);//lock count mutex
	return count;//return if there are pending tasks
}

//@brief: wait for the queue to be empty and all running tasks to finish
//@param timeout: duration to wait for queue to complete
//@return: true/false if the queue was/wasn't complete before the timeout
template <class Rep, class Period>
bool ThreadPool::TaskQueue::waitComplete(std::chrono::duration<Rep, Period> timeout) const {
	std::unique_lock<std::mutex> iLock(iMut);//lock counter mutex
	if(count) return cCv.wait_for(iLock, timeout, [&](){return count == 0;});//wait for time to expire or queue to complete
	return true;//queue was already complete
}

#endif//_THREADPOOL_H_
