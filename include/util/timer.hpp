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

#ifndef _TIMER_HPP_
#define _TIMER_HPP_

#include <chrono>
#include <ostream>

struct  Timer {
	std::chrono::high_resolution_clock::time_point tp;
	Timer() : tp(std::chrono::high_resolution_clock::now()) {}

	//@brief      : compute time since previous call (or construction)
	//@param reset: should the reference time be reset
	//@return     : time in seconds
	double poll(const bool reset = true) {
		std::chrono::high_resolution_clock::time_point now = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = now - tp;
		if(reset) set(now);
		return elapsed.count();
	}

	//@brief    : set the reference time
	//@param ref: reference time point
	void set(const std::chrono::high_resolution_clock::time_point ref = std::chrono::high_resolution_clock::now()) {tp = ref;}

	//@brief   : print a number of seconds formatted nicely as dy:hr:mn:sc
	//@param s : number of seconds
	//@param os: ostream to print to
	static void PrintSeconds(const double s, std::ostream& os) {
		//std::gmtime is creating a null pointer sometimes...
		// const time_t iSec = (time_t) s;//get number of complete seconds
		// std::tm* pTm = std::gmtime(&iSec);//format
		uint64_t tm_sec = (uint64_t) s;
		const uint64_t tm_year = tm_sec / 31536000ULL;
		tm_sec -= tm_year * 31536000ULL;
		const uint64_t tm_yday = tm_sec / 86400ULL;
		tm_sec -= tm_yday * 86400ULL;
		const uint64_t tm_hour = tm_sec / 3600ULL;
		tm_sec -= tm_hour * 3600ULL;
		const uint64_t tm_min = tm_sec / 60ULL;
		tm_sec -= tm_min  * 60ULL;
		std::tm t;
		t.tm_year = (int) tm_year;
		t.tm_yday = (int) tm_yday;
		t.tm_hour = (int) tm_hour;
		t.tm_min  = (int) tm_min ;
		t.tm_sec  = (int) tm_sec ;
		std::tm* pTm = &t;
		if       (pTm->tm_yday > 0) {
			os << pTm->tm_yday << ':';
			if(pTm->tm_hour < 10) os << '0';
			os << pTm->tm_hour << ':'; 
			if(pTm->tm_min  < 10) os << '0';
			os << pTm->tm_min  << ':'; 
			if(pTm->tm_sec  < 10) os << '0';
			os << pTm->tm_sec;
		} else if(pTm->tm_hour > 0) {
			os << pTm->tm_hour << ':'; 
			if(pTm->tm_min  < 10) os << '0';
			os << pTm->tm_min  << ':'; 
			if(pTm->tm_sec  < 10) os << '0';
			os << pTm->tm_sec;
		} else if(pTm->tm_min  > 0) {
			os << pTm->tm_min  << ':'; 
			if(pTm->tm_sec  < 10) os << '0';
			os << pTm->tm_sec;
		} else {
			os << pTm->tm_sec  << 's';
		}
	}
};

#endif//_TIMER_HPP_
