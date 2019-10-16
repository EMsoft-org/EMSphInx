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

#include "util/fft.hpp"//for ffts in DiscreteSHT

int main(int argc, char *argv[]) {
	//parse arguments (file names / bandwidth)
	typedef double Real;
	if(2 != argc) {
		std::cout << "usage: " << argv[0] << "bandWidth\n";
		std::cout << "\tbandWidth  : max bandwidth to build FFTW wisdom fore\n";
		std::cout << "\t             some reasonable values are 63, 95, 158, 263\n";
		return EXIT_FAILURE;
	}
	const size_t bw = std::strtoul(argv[1], NULL, 0);

	//start by building wisdom for all spherical grid sizes up to the size for bw
	const size_t Nt = bw / 2 + 2;//number of rings in largest bandwidth grid
	for(size_t y = 0; y < Nt; y++) {
		std::cout << "\rplanning 1D fft for ring " << y+1 << "/" << Nt;
		std::cout.flush();
		fft::RealFFT<Real> plan(std::max<size_t>(1, 8 * y), fft::flag::Plan::Patient);
	}
	std::cout << '\n';

	//now build some fast sizes
	static const size_t FastSize[27] = {
		25, 32, 38, 41, 53, 63, 68, 74, 88, 95, 113, 122, 123, 158, 172, 188, 203, 221, 263, 284, 313, 338, 365, 368, 438, 473, 515
	};
	for(size_t i = 0; i < 27; i++) {
		if(FastSize[i] > bw) break;
		std::cout << "\rplanning 3D fft for BW = " << FastSize[i];
		std::cout.flush();
		fft::SepRealFFT3D<Real> plan(FastSize[i], fft::flag::Plan::Patient);
	}
	std::cout << '\n';

	return 0;
}
