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
//       test program for classes in include/sht/sht_xcorr.hpp        //
////////////////////////////////////////////////////////////////////////

//!!!!!!!!
// Tests are still missing for normalized cross correlation
//!!!!!!!!

#include <iostream>

namespace emsphinx {

	namespace sphere {
		//@brief     : test spherical cross correlation
		//@param bw  : bandwidth to test
		//@param mir : should a z mirror plane be used
		//@param nFld: degree of rotational symmtry about z axis
		//@param qu  : optional location to write applied rotation
		//@param qr  : optional location to write registered rotation
		//@return    : angle between applied and registered rotation
		//@note      : only tests unnormalized cross correlation
		template <typename Real> Real testCorr(const size_t bw, const bool mir, const size_t nFld, Real* qu = NULL, Real* qr = NULL);

		//@brief     : test normalized spherical cross correlation
		//@param bw  : bandwidth to test
		//@param mir : should a z mirror plane be used
		//@param nFld: degree of rotational symmtry about z axis
		//@param qu  : optional location to write applied rotation
		//@param qr  : optional location to write registered rotation
		//@return    : angle between applied and registered rotation
		//@note      : only tests unnormalized cross correlation
		template <typename Real> Real testNCorr(const size_t bw, const bool mir, const size_t nFld, Real* qu = NULL, Real* qr = NULL);

		//@brief   : run all spherical cross correlation test for a range of parameters
		//@param os: location to write errors
		//@return  : true / false if the tests pass/fail
		template <typename Real> bool runTests(std::ostream& os);
	}

}


int main() {
	try {
		return emsphinx::sphere::runTests<double>(std::cout) ? EXIT_SUCCESS : EXIT_FAILURE;
	} catch(std::exception& e) {
		std::cout << "caught: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
}

#include <random>
#include <limits>

#include "sht/wigner.hpp"
#include "sht/square_sht.hpp"
#include "sht/sht_xcorr.hpp"

#include "xtal/rotations.hpp"
#include "xtal/quaternion.hpp"
#include "idx/master.hpp"

namespace emsphinx {

	namespace sphere {

		//@brief     : generate a random spherical function
		//@param gen : random generator
		//@param dim : side length of spherical grid
		//@param mir : should the function have a mirror plane at the equator
		//@param nFld: order of rotational symmetry about z axis
		//@return    : random spherical function with specified size and symmetry
		template <typename Real>
		emsphinx::MasterPattern<Real> randomSphere(std::mt19937_64& gen, const size_t dim, const bool mir, const size_t nFld) {
			//start with a random north hemisphere
			std::uniform_real_distribution<Real> dist(-1,1);//build a random distribution
			emsphinx::MasterPattern<Real> mp(dim);//create an empty master pattern
			mp.lyt = square::Layout::Legendre;
			for(Real& v : mp.nh) v = dist(gen);//fill

			//apply mirror symmetry if needed
			if(mir) {
				mp.sh = mp.nh;//z mirror
			} else {
				for(Real& v : mp.sh) v = dist(gen);//build random south hemisphere
				mp.matchEquator();//make sure equators are the same
			}

			//apply rotational symmetry if needed
			if(nFld > 1) mp.makeNFold(nFld);

			return mp;
		}

		//@brief     : generate a random rotation
		//@param gen : random generator
		//@return    : random rotation
		template <typename Real>
		xtal::Quat<Real> randomRotation(std::mt19937_64& gen) {
			std::uniform_real_distribution<Real> dist(-1,1);//build a random distribution
			return xtal::Quat<Real>(std::fabs(dist(gen)), dist(gen), dist(gen), dist(gen)).normalize();
		}

		//@brief     : generate the spectra of a random pattern on the sphere and rotate it
		//@param bw  : bandwidth of spectra to generate
		//@param mir : should the spectra correspond to a function with z mirror symmetry
		//@param nFld: order of rotational symmetry about z axis in reference pattern
		//@param flm : location to write spectra of original unrotated (reference) function
		//@param gln : location to write spectra of rotated function
		//@return    : rotation applied to gln as quaternion
		template <typename Real>
		xtal::Quat<Real> randomPair(const size_t bw, const bool mir, const size_t nFld, std::vector< std::complex<Real> >& flm, std::vector< std::complex<Real> >& gln) {
			//build a random generator
			const unsigned int seed = 0;//constant for deterministic behavior or time(NULL) for random
			std::mt19937_64 gen(seed);//64 bit mersenne twister

			//generate a random function with the desired symmetry
			const size_t dim = bw + (bw % 2 == 0 ? 3 : 2);
			emsphinx::MasterPattern<Real> mp = randomSphere<Real>(gen, dim, mir, nFld);

			//compute the SHT of the master pattern
			emsphinx::MasterSpectra<Real> spec(mp, (uint16_t)bw, false);
			flm = std::vector< std::complex<Real> >(spec.data(), spec.data() + bw * bw);

			//exactly enforce rotational symmetry (since the real space implementation is a bit kludgey for # ring pts % nFld != 0)
			if(nFld > 1) {
				for(size_t m = 0; m < bw; m++) {
					if(0 != m % nFld) {
						std::fill(flm.begin() + m * bw, flm.begin() + (m+1) * bw, std::complex<Real>(0));
					}
				}
			}

			//apply a random rotation
			Real eu[3];
			xtal::Quat<Real> qu = randomRotation<Real>(gen);
			xtal::qu2zyz(qu.data(), eu);//convert to zyz euler angles
			gln.resize(flm.size());
			wigner::rotateHarmonics(bw, flm.data(), gln.data(), eu);
			return qu;
		}

		//@brief     : test spherical cross correlation
		//@param bw  : bandwidth to test
		//@param mir : should a z mirror plane be used
		//@param nFld: degree of rotational symmetry about z axis
		//@param qu  : optional location to write applied rotation
		//@param qr  : optional location to write registered rotation
		//@return    : angle between applied and registered rotation
		template <typename Real> Real testCorr(const size_t bw, const bool mir, const size_t nFld, Real* qu, Real* qr) {
			//generate a pair of rotated harmonics
			std::vector< std::complex<Real> > flm, gln;
			xtal::Quat<Real> u = randomPair(bw, mir, nFld, flm, gln);

			//compute the rotation of the peak cross correlation
			Real eu[3], r[4];
			sphere::Correlator<Real> s2(bw);//spherical cross correlation calculator
			Real xc = s2.correlate(flm.data(), gln.data(), mir, nFld, eu, true);
			xtal::zyz2qu(eu, r);//convert to a quaternion

			//save rotations if needed
			//these are conjugated to convert from crystal->sample to sample->crystal so that symmetry operators are correct
			if(NULL != qu) xtal::quat::conj(u.data(), qu);
			if(NULL != qr) xtal::quat::conj(r       , qr);

			//compare registered with applied rotations
			const Real dot = std::inner_product(r, r + 4, u.data(), Real(0));
			return 180.0 * std::acos(std::min(dot, Real(1))) / emsphinx::Constants<Real>::pi;
		}

		//@brief     : test spherical cross correlation
		//@param bw  : bandwidth to test
		//@param mir : should a z mirror plane be used
		//@param nFld: degree of rotational symmetry about z axis
		//@param qu  : optional location to write applied rotation
		//@param qr  : optional location to write registered rotation
		//@return    : angle between applied and registered rotation
		template <typename Real> Real testNCorr(const size_t bw, const bool mir, const size_t nFld, Real* qu, Real* qr) {
			//build a random generator
			const unsigned int seed = 0;//constant for deterministic behavior or time(NULL) for random
			std::mt19937_64 gen(seed);//64 bit mersenne twister

			//build a reasonable mask (something with minimal symmetry that has points in both hemispheres)
			const size_t dim = bw + (bw % 2 == 0 ? 3 : 2);
			std::vector<Real> mask(dim * dim * 2, 0), norms(dim * dim * 3);

			//get normal directions for each square legendre grid point and build mask limits
			square::legendre::normals(dim, norms.data());
			const Real tMin = (Real) -0.52359877559829887307710723054658L;// -30 degrees
			const Real tMax = (Real)  1.04719755119659774615421446109320L;//  60 degrees
			const Real zMax = (Real)  0.70710678118654752440084436210485L;// maximum negative z (45 degrees)

			//set mask points
			for(size_t j = 0; j < dim; j++) {
				for(size_t i = 0; i < dim; i++) {
					const size_t idx = j * dim + i;
					Real * n = norms.data() + 3 * idx;//get normal for this point (in NH)
					const Real theta = std::atan2(n[1], n[0]);//get angle
					if(theta >= tMin && theta <= tMax) {//inside of wedge
						mask[idx] = Real(1);//set north hemisphere for any lattitude
						if(n[2] < zMax) mask[idx+dim*dim] = Real(1);
					}
				}
			}

			//build a SHT calculator
			square::DiscreteSHT<Real> sht(dim, bw, square::Layout::Legendre);

			//generate a random function and do a round trip to remove frequencies above bandwidth
			std::vector< std::complex<Real> > flm(bw * bw), flm2(bw * bw), gln(bw * bw), mlm(bw * bw);
			emsphinx::MasterPattern<Real> mp = randomSphere<Real>(gen, dim, mir, nFld);
			sht.analyze(mp.nh.data(), mp.sh.data(), flm.data());
			sht.synthesize(flm.data(), mp.nh.data(), mp.sh.data());

			//now compute SHT of function squared
			for(Real& r : mp.nh) r *= r;
			for(Real& r : mp.sh) r *= r;
			sht.analyze(mp.nh.data(), mp.sh.data(), flm2.data());

			//next do the same for the mask
			sht.analyze(mask.data(), mlm.data());
			sht.synthesize(mlm.data(), mask.data());

			//now apply a random rotation to the original function
			Real eu[3];
			xtal::Quat<Real> u = randomRotation<Real>(gen);
			xtal::qu2zyz(u.data(), eu);//convert to zyz euler angles
			wigner::rotateHarmonics(bw, flm.data(), gln.data(), eu);

			//mask the rotated function
			std::vector<Real> work(dim * dim * 2);
			sht.synthesize(gln.data(), work.data());
			std::transform(work.cbegin(), work.cend(), mask.cbegin(), work.begin(), std::multiplies<Real>());
			sht.analyze(work.data(), gln.data());

			//build normalized spherical cross correlation calculator
			sphere::NormalizedCorrelator<Real> s2(bw, flm.data(), flm2.data(), mir, nFld, mlm.data());

			//compute the rotation of the peak normalized cross correlation
			Real r[4];
			Real xc = s2.correlate(gln.data(), eu, true);
			xtal::zyz2qu(eu, r);//convert to a quaternion

			//save rotations if needed
			//these are conjugated to convert from crystal->sample to sample->crystal so that symmetry operators are correct
			if(NULL != qu) xtal::quat::conj(u.data(), qu);
			if(NULL != qr) xtal::quat::conj(r       , qr);

			//compare registered with applied rotations
			const Real dot = std::inner_product(r, r + 4, u.data(), Real(0));
			// return 100;
			return 180.0 * std::acos(std::min(dot, Real(1))) / emsphinx::Constants<Real>::pi;
		}

		//@brief   : run all spherical cross correlation test for a range of parameters
		//@param os: location to write errors
		//@return  : true / false if the tests pass/fail
		template <typename Real> bool runTests(std::ostream& os) {
			//build tolerance and sizes
			Real eps = std::cbrt(std::numeric_limits<float>::epsilon());//this is ~6e-6 degrees for double and 0.005 degrees for float
			const std::vector<size_t> sizes = {
				53, 68, 88, 113, 123, 158,            //some fast sizes
				54, 55, 56, 57 , 58 , 59 , 60, 62, 64,//some sizes that need to be 0 padded for speed
			};

			//now loop over sizes doing testing
			xtal::Quat<Real> qu, qr;
			os << "testing symmetry free cross correlation" << std::endl;
			for(const size_t& i : sizes) {
				os << '\t' << i << '\n';
				const Real delta = testCorr<Real>(i, false, 1, qu.data(), qr.data());
				if(delta > eps) {
					os << "registered rotation is " << delta << " degrees away from applied\n";
					os << "\tbandwidth : " << i << '\n';
					os << "\tapplied   : " << qu << '\n';
					os << "\tregistered: " << qr << '\n';
					return false;
				}
			}

			//now loop over sizes doing normalized testing
			Real epsN = eps * 10;//allow additional error introduced by the masking
			os << "testing symmetry free normalized cross correlation" << std::endl;
			for(const size_t& i : sizes) {
				os << '\t' << i << '\n';
				const Real delta = testNCorr<Real>(i, false, 1, qu.data(), qr.data());
				if(delta > epsN) {
					os << "registered rotation is " << delta << " degrees away from applied\n";
					os << "\tbandwidth : " << i << '\n';
					os << "\tapplied   : " << qu << '\n';
					os << "\tregistered: " << qr << '\n';
					os << epsN << '\n';
					return false;
				}
			}

			//build combination of expected z mirror and nFold from basic crystallographic point groups
			std::vector<xtal::PointGroup> pgs = {
				//                          mir nFld
				xtal::PointGroup("112"  ),// F   2
				xtal::PointGroup("11m"  ),// T   1
				xtal::PointGroup("112/m"),// T   2
				xtal::PointGroup("3"    ),// F   3
				xtal::PointGroup("4"    ),// F   4
				xtal::PointGroup("4/m"  ),// T   4
				xtal::PointGroup("6"    ),// F   6
				xtal::PointGroup("6/m"  ),// T   6
			};

			//loosen tolerance some to accommodate symmetry
			eps = std::sqrt(eps) * 5;//this is ~0.012 and 0.35 degrees for double and float respectively

			//test each for a range of bandwidths
			for(const xtal::PointGroup& pg : pgs) {//loop over point groups
				os << "testing cross correlation for point group " << pg.name() << std::endl;
				for(size_t i = 53; i < 64; i++) {//smaller sample of faster sizes since there are a lot of tests in this loop
					Real delta = testCorr<Real>(i, pg.zMirror(), pg.zRot(), qu.data(), qr.data());//no symmetry accounted for
					if(delta > eps) {//make sure we don't have a symmetry operator
						xtal::Quat<Real> qw;
						pg.disoQu(qu.data(), qr.data(), qw.data());
						const Real deltaSym = 180.0 * std::acos(std::min(qw.w, Real(1))) / emsphinx::Constants<Real>::pi;
						delta = std::min(delta, deltaSym);
					}

					if(delta > eps) {
						os << "Registered rotation is " << delta << " degrees away from applied\n";
						os << "\tbandwidth  : " << i << '\n';
						os << "\tapplied    : " << qu << '\n';
						os << "\tregistered : " << qr << '\n';
						os << "\tpoint group: " << pg.name() << '\n';
						return false;
					}
				}
			}

			//normalized test each for a range of bandwidths
			for(const xtal::PointGroup& pg : pgs) {//loop over point groups
				os << "testing normalized cross correlation for point group " << pg.name() << std::endl;
				for(size_t i = 53; i < 64; i++) {//smaller sample of faster sizes since there are a lot of tests in this loop
					Real delta = testNCorr<Real>(i, pg.zMirror(), pg.zRot(), qu.data(), qr.data());//no symmetry accounted for
					if(delta > eps) {//make sure we don't have a symmetry operator
						xtal::Quat<Real> qw;
						pg.disoQu(qu.data(), qr.data(), qw.data());
						const Real deltaSym = 180.0 * std::acos(std::min(qw.w, Real(1))) / emsphinx::Constants<Real>::pi;
						delta = std::min(delta, deltaSym);
					}

					if(delta > eps) {
						os << "Registered rotation is " << delta << " degrees away from applied\n";
						os << "\tbandwidth  : " << i << '\n';
						os << "\tapplied    : " << qu << '\n';
						os << "\tregistered : " << qr << '\n';
						os << "\tpoint group: " << pg.name() << '\n';
						return false;
					}
				}
			}

			//if we made it this far everything passed
			return true;
		}
	}

}
