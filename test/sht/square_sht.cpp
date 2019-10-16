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

////////////////////////////////////////////////////////////////////////
//  test program for functions/classes in include/sht/square_sht.hpp  //
////////////////////////////////////////////////////////////////////////

#include <iostream>
#include "sht/square_sht.hpp"

namespace emsphinx {

	namespace square {
		//@brief       : test DiscreteSHT for self consistent behavior for a specific condition
		//@param bw    : bandwidth to test
		//@param lyt   : layout of square grid
		//@param maxEps: maximum allowable absolute error for a single harmonic coefficient
		//@param avgEps: maximum average error for a harmonic coefficient
		//@param os    : location to write error messages
		//@return      : true/false if test passes/fails
		template <typename Real> bool testSingleSHT(const size_t bw, Layout lyt, const Real maxEps, const Real avgEps, std::ostream& os);

		//@brief   : test DiscreteSHT for self consistent behavior over a range of reasonable conditions
		//@param os: location to write error messages
		//@return  : true/false if test passes/fails
		bool testSHT(std::ostream& os);
	}

}

int main() {
	//select location to write output
	std::ostream& os = std::cout;

	//do tests
	try {
		return emsphinx::square::testSHT(os) ? EXIT_SUCCESS : EXIT_FAILURE;
	} catch(std::exception& e) {
		os << "caught: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
}

#include <random>
#include <vector>
#include <limits>
#include <iomanip>

namespace emsphinx {

	namespace square {
		//@brief       : test DiscreteSHT for self consistent behavior
		//@param bw    : bandwidth to test
		//@param lyt   : layout of square grid
		//@param maxEps: maximum allowable absolute error for a single harmonic coefficient
		//@param avgEps: maximum average error for a harmonic coefficient
		//@param os    : location to write error messages
		//@return      : true/false if test passes/fails
		template <typename Real> bool testSingleSHT(const size_t bw, Layout lyt, const Real maxEps, const Real avgEps, std::ostream& os) {
			try {
				//build a random distribution
				const unsigned int seed = 0;//constant for deterministic behavior or time(NULL) for random
				std::mt19937_64 gen(seed);//64 bit mersenne twister
				std::uniform_real_distribution<Real> dist(-1,1);

				//generate a random spectra
				std::vector< std::complex<Real> > refSpec(bw * bw, std::complex<Real>(0));
				std::vector< std::complex<Real> > newSpec(bw * bw, std::complex<Real>(0));
				for(size_t m = 0; m < bw; m++) {
					for(size_t l = m; l < bw; l++) {
						if(m == 0)
							refSpec[m*bw+l] = std::complex<Real>(dist(gen), 0        );
						else
							refSpec[m*bw+l] = std::complex<Real>(dist(gen), dist(gen));
					}
				}

				//build discrete SHT calculators
				const size_t dimLam = 2 * bw + 1;
				const size_t dimLeg = bw + (bw % 2 == 0 ? 3 : 2);//only odd side lengths are supported
				const size_t dim = lyt == Layout::Legendre ? dimLeg : dimLam;
				DiscreteSHT<Real> sht(dim, bw, lyt);

				//do a round trip calculation
				std::vector<Real> func(dim * dim * 2);
				sht.synthesize(refSpec.data(), func.data());
				sht.analyze(func.data(), newSpec.data());

				//compare computed spectra
				Real avgErr = 0;
				for(size_t m = 0; m < bw; m++) {
					for(size_t l = 0; l < bw; l++) {
						std::complex<Real> delta = newSpec[m*bw+l] - refSpec[m*bw+l];
						Real rDelta = std::max(std::fabs(delta.real()), std::fabs(delta.imag()));
						if(rDelta > maxEps) {
							os << std::setprecision(std::numeric_limits<Real>::digits10);
							os << "original  : " << refSpec[m*bw+l] << '\n';
							os << "round trip: " << newSpec[m*bw+l] << '\n';
							os << "delta (%) : " << rDelta << '\n';
							os << "limit     : " << maxEps << '\n';
							os.flush();
							return false;
						}
						avgErr += rDelta;
					}
				}
				avgErr /= bw * bw;
				if(avgErr > avgEps) {
					os << "average error (" << avgEps << ") exceeds limit (" << avgEps << ")" << std::endl;
					return false;
				}
			} catch (std::exception& e) {
				os << e.what() << std::endl;
				return false;
			}
			return true;
		}


		//@brief   : test DiscreteSHT for self consistent behavior over a range of reasonable conditions
		//@param os: location to write error messages
		//@return  : true/false if test passes/fails
		bool testSHT(std::ostream& os) {
			//test double precision transform
			const double dMax = 0.005  ;//max absolute error
			const double dAvg = 0.00005;//max average error
			os << "testing round trip double transforms\n";
			os << "\tlegendre grid\n";
			for(size_t i = 4; i <= 384; i++) {
				if(0 == i % 8) os << '\t' << i << '\n';
				if(!square::testSingleSHT<double>(i, square::Layout::Legendre, dMax, dAvg, os)) {
					os << "bw == " << i << '\n';
					return false;
				}
			}
			os << "\tlambert grid\n";
			for(size_t i = 4; i <= 128; i++) {
				if(0 == i % 8) os << '\t' << i << '\n';
				if(!square::testSingleSHT<double>(i, square::Layout::Lambert , dMax, dAvg, os)) {
					os << "bw == " << i << '\n';
					return false;
				}
			}

		#ifdef EM_USE_F
			//test single precision transform
			const float fMax = 0.020  ;//max absolute error
			const float fAvg = 0.00020;//max average error
			os << "testing round trip float  transforms\n";
			os << "\tlegendre grid\n";
			for(size_t i = 4; i <= 192; i++) {
				if(0 == i % 8) os << '\t' << i << '\n';
				if(!square::testSingleSHT<float >(i, square::Layout::Legendre, fMax, fAvg, os)) {
					os << "bw == " << i << '\n';
					return false;
				}
			}
			os << "\tlambert grid\n";
			for(size_t i = 4; i <= 64; i++) {
				if(0 == i % 8) os << '\t' << i << '\n';
				if(!square::testSingleSHT<float >(i, square::Layout::Lambert , fMax, fAvg, os)) {
					os << "bw == " << i << '\n';
					return false;
				}
			}
		#endif

			//if we made it this far everything passed
			return true;
		}
	}
	
}

