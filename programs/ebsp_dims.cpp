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

#include <sstream>
#include <set>

#include "modality/ebsd/pattern.hpp"
#include "util/timer.hpp"

int main(int argc, char *argv[]) {
	const bool binToFloat = false;//should binned values be converted up to a float or kept as their current type

	//sanity check argument count
	if(2 != argc) {
		std::cout << "usage: " << argv[0] << " inputFile\n";
		std::cout << "\tinputFile  - pattern file to read (*.ebsp)\n";
		return EXIT_FAILURE;
	}

	//get some pattern info
	Timer t;
	emsphinx::ebsd::OxfordPatternFile pats(argv[1]);

	std::cout << "found " << pats.numPat() << " patterns in " << argv[1] << ":\n";
	std::cout << "\twidth : " << pats.width () << '\n';
	std::cout << "\thegiht: " << pats.height() << '\n';
	std::cout << "\ttype  : ";
	switch(pats.pixelType()) {
		case emsphinx::ImageSource::Bits::U8 : std::cout << "8 bit\n"  ; break;
		case emsphinx::ImageSource::Bits::U16: std::cout << "16 bit\n" ; break;
		case emsphinx::ImageSource::Bits::F32: std::cout << "float\n"  ; break;
		case emsphinx::ImageSource::Bits::UNK: // intentional fall through
		default: throw std::logic_error("unknown pixel type");
	}
	std::cout << "\tbytes : " << pats.imBytes() << "\n";
	double mb = double(pats.imBytes() * pats.numPat()) / (1024*1024);
	std::cout << "total size: ";
	if(mb > 1024) {
		std::cout << mb / 1024 << " GB\n";
	} else {
		std::cout << mb        << " MB\n";
	}

	//now loop over patterns accumulating all possible indices
	const size_t batchSize = 100;
	std::vector<double> vx, vy;
	std::set<double> sx, sy;

	//finally loop over patterns writing into file
	std::vector<char> buff(pats.imBytes() * batchSize);
	const size_t numBatch = (pats.numPat() + batchSize - 1) / batchSize;
	for(size_t i = 0; i <numBatch; i++) {
		//extract patterns and accumulate coordinates
		pats.extract(buff.data(), batchSize, &vx, &vy);
		sx.insert(vx.begin(), vx.end());
		sy.insert(vy.begin(), vy.end());

		//clear work space
		vx.clear();
		vy.clear();
	}
	std::cout << t.poll() << "s to read patterns\n";

	//print details of found dims
	std::cout << "found " << sx.size() << " x and " << sy.size() << " y coordinates\n\n";
	if(sx.size() * sy.size() != pats.numPat()) {
		std::cout << "X:";
		for(const double& x : sx) std::cout << ' ' << x; 
		std::cout << "\n\nY:";
		for(const double& y : sy) std::cout << ' ' << y; 
		std::cout << "\n";
	}

	return EXIT_SUCCESS;
}
