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

#include <set>
#include <iostream>

//@brief  : compute the base 2 log of a number
//@param v: integer to compute base 2 log of
//@return : floor(log(v, 2))
//@note   : from stanford bit twiddles https://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
uint32_t log2(uint32_t v) {
	uint32_t r = (v > 0x0000FFFF) << 4; v >>= r;
	uint32_t s = (v > 0x000000FF) << 3; v >>= s; r |=  s      ;
	         s = (v > 0x0000000F) << 2; v >>= s; r |=  s      ;
	         s = (v > 0x00000003) << 1; v >>= s; r |=  s      ;
	                                             r |= (v >> 1);
	return r;
}

int main() {
	//specify the largest number we're interested and compute its base 2 log
	const uint32_t lb = 25 ;//lower bound on meaningful bandwidths
	const uint32_t ub = 2000;//upper bound on fft sizes
	const uint32_t l2 = log2(ub);//2^l2 is the largest value to consider

	//build initial set of small primes that FFTW reccomends
	std::set<uint32_t> factors;
	factors.insert(1);
	factors.insert(2);
	factors.insert(3);
	factors.insert(5);
	factors.insert(7);
	
	//build set of products of small primes (these are good fft sizes)
	for(uint32_t m = 0; m < l2; m++) {
		std::set<uint32_t> products;
		for(const uint32_t i : factors) {
			for(const uint32_t j : factors) {
				const uint32_t p = i * j;
				if(p <= ub) {
					products.insert(p);
				}
			}
		}
		factors = products;
	}

	//factors now contains all products of small primes <= ub, print valid bandwidths
	for(const uint32_t i : factors) {//loop over valid fft sizes
		if(1 == i % 2) {//check if there is a bandwidth such that 2 * bw - 1 == i
			const uint32_t bw = (i + 1) / 2;
			if(bw >= lb) std::cout << bw << '\n';
		}
	}

	return 0;
}
