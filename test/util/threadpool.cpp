/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019-2019, De Graef Group, Carnegie Mellon University *
 * All rights reserved.                                                *
 *                                                                     *
 * Author: William C. Lenthe                                           *
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
 *                                                                     *
 * Interested in a commercial license? Contact:                        *
 *                                                                     *
 * Center for Technology Transfer and Enterprise Creation              *
 * 4615 Forbes Avenue, Suite 302                                       *
 * Pittsburgh, PA 15213                                                *
 *                                                                     *
 * phone. : 412.268.7393                                               *
 * email  : innovation@cmu.edu                                         *
 * website: https://www.cmu.edu/cttec/                                 *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>

//@brief : check if thead pool works
//@return: true / false if tests pass / fail
bool testThreadPool(std::ostream& os);

int main() {
	return testThreadPool(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include "util/threadpool.hpp"
#include "util/timer.hpp"

#include <cmath>
#include <algorithm>

//@brief : check if thead pool works
//@return: true / false if tests pass / fail
bool testThreadPool(std::ostream& os) {
	//build a thread pool
	ThreadPool pool;
	os << "checking if serial and parallel calculations are consistent\n";

	//make a function to do some work
	std::vector<int> buff(pool.size() * pool.size(), 0);
	auto func = [&](const size_t idx, const size_t tid) {buff[idx] = 1;};

	//make sure serial and calculation works
	for(size_t i = 0; i < buff.size(); i++) func(i, 0);
	if(!std::all_of(buff.cbegin(), buff.cend(), [](const int& v){return v == 1;})) {
		os << "got wrong answer from serial computation\n";
		return false;
	}
	std::fill(buff.begin(), buff.end(), 0);

	//make sure the parallel fill is correct
	for(size_t i = 0; i < buff.size(); i++) pool.schedule(std::bind(func, i, std::placeholders::_1));
	pool.waitAll();
	if(!std::all_of(buff.cbegin(), buff.cend(), [](const int& v){return v == 1;})) {
		os << "got wrong answer from parallel computation\n";
		return false;
	}
	std::fill(buff.begin(), buff.end(), 0);

	//make sure that thread index is working for passed functions
	os << "checking if thread index works correctly\n";
	std::vector< std::mutex > muts(pool.size());//1 mutex per thread
	auto funcMut = [&](const size_t tid) {
		for(size_t i = 0; i < pool.size(); i++) {
			if(!muts[tid].try_lock()) {//try to get this thread's mutex (there better not be any competition)
				os << "failed to lock thread mutex\n";
				os.flush();
				std::terminate();
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(100));//hold the lock for awhile
			muts[tid].unlock();//release the lock for the next loop
		}
	};
	for(size_t i = 0; i < buff.size(); i++) pool.schedule(funcMut);
	pool.waitAll();

	//make a slow function to burn time
	os << "checking speedup\n";
	auto funcWait = [&](const size_t tid) {std::this_thread::sleep_for(std::chrono::milliseconds(100));};//sleep for 0.1s
	Timer t;
	for(size_t i = 0; i < 10; i++) {//10 * 0.1 = 1s of work per thread
		for(size_t j = 0; j < pool.size(); j++) pool.schedule(funcWait);
	}
	pool.waitAll();
	double tWrk = t.poll();
	if(std::fabs(tWrk - 1.0) > 0.02) {//2% error
		os << "1s of work per thread took " << tWrk << "s\n";
		return false; 
	}

	//try the same thing with a single worker pool
	os << "checking destructor\n";
	{
		ThreadPool pool1(1);
		t.poll();
		pool1.schedule([&](const size_t tid){
			buff[0] = 0;
			std::this_thread::sleep_for(std::chrono::milliseconds(1000));
			buff[0] = 1;
		});
	}//test that destructor calls waitAll()
	tWrk = t.poll();
	if(1 != buff[0]) {
		os << "destructor didn't wait for thread completion\n";
		return false;
	}
	if(std::fabs(tWrk - 1.0) > 0.02) {//2% error
		os << "1s of work for 1 thread thread took " << tWrk << "s\n";
		return false; 
	}

	//make sure timeout works
	os << "checking timeout\n";
	pool.schedule([&](const size_t tid){std::this_thread::sleep_for(std::chrono::milliseconds(2000));});//2s of work
	t.poll();
	const bool done = pool.waitAll(std::chrono::milliseconds(500));//give the thread only 0.5s to finish
	tWrk = t.poll();
	if(done) {
		os << "too short waitout returned true\n";
		return false;
	}
	if(std::fabs(tWrk - 0.5) > 0.01) {//2% error
		os << "short waitout isn't accurate\n";
		return false;
	}
	pool.waitAll();

	//make sure task clearing works
	os << "checking task clearing\n";
	for(size_t i = 0; i < pool.size(); i++)//schedule short piece of work for every thread
		pool.schedule([&](const size_t tid){std::this_thread::sleep_for(std::chrono::milliseconds(600));});//0.5s of work
	for(size_t i = 0; i < pool.size(); i++)//schedule long piece of work for every thread
		pool.schedule([&](const size_t tid){std::this_thread::sleep_for(std::chrono::milliseconds(5000));});//5s of work
	std::this_thread::sleep_for(std::chrono::milliseconds(100));//give sleeping workers a chance to wake up
	if(pool.size() * 2 != pool.items()) {//make sure no workers have finished the short task and started the long task
		os << "unexpected queuing time\n";
		return false;
	}
	pool.clear();//remove the long work (shouldn't be started by any threads yet)
	t.poll();
	pool.waitAll();//wait for all the work to finish
	tWrk = t.poll();
	if(std::fabs(tWrk - 0.5) > 0.02) {//4% error
		os << "task clearing didn't remove unstarted work: " << tWrk << "\n";
		return false;
	}

	//if we made it this far all tests passed
	return true;
}
