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

//@brief : check if timer works
//@return: true / false if tests pass / fail
bool testTimer(std::ostream& os);

int main() {
	return testTimer(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <thread>
#include <vector>
#include <numeric>
#include <cmath>

#include "util/timer.hpp"

//@brief : check if timer works
//@return: true / false if tests pass / fail
bool testTimer(std::ostream& os) {
	Timer t;//build a timer
	std::vector<double> times(3);//build space to store timing results (3 trials)
	for(size_t ms = 8; ms < 4096; ms *= 2) {
		//convert from ms to fractional seconds and print info
		const double dMs = double(ms) / 1000;
		os << "testing " << times.size() << " * " << dMs << "s: ";
		os.flush();

		//time sleep for a few runs
		for(double& v : times) {
			t.poll();
			std::this_thread::sleep_for(std::chrono::milliseconds(ms));
			v = t.poll();
		}
		double v = std::accumulate(times.cbegin(), times.cend(), 0.0) / times.size();

		//print result and check
		std::cout << v << '\n';
		if(std::fabs(dMs - v) * 1000 > 6) {
			os << "error outside of 4ms tolerance\n";
			return false;
		}
	}

	//if we made it this far all tests passed
	return true;
}
