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

#include <iostream>

namespace xtal {
	//@brief    : test rotation conversions for 2 and 3 way self consistency, e.g. for 2 way make sure eu2ro, ro2eu give back the input
	//@param n  : euler grid density
	//@param do3: should 3 way tests be included (this can significantly increase computation time)
	//@param os : location to write errors
	//@return   : true / false if tests pass / fail
	template <typename Real> bool testRound(const size_t n, const bool do3, std::ostream& os);

	//@brief : run all rotation tests
	//@return: true / false if tests pass / fail
	template <typename Real> bool testRots(std::ostream& os);
}

int main() {
	std::ostream& os = std::cout;
	return xtal::testRots<float>(os) && xtal::testRots<double>(os) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <vector>
#include <cmath>
#include <algorithm>//transform
#include <numeric>//iota

#include "xtal/rotations.hpp"

namespace xtal {

	//@brief    : compute maximum of the element wise distance between 2 vectors
	//@param a  : first vector
	//@param b  : second vector
	//@param len: length of vectors
	//@return   : L1 distance
	template <typename Real> Real maxDelta(Real const * const a, Real const * const b, const size_t len) {
		Real maxDist = 0;
		for(size_t i = 0; i < len; i++) {
			const Real dist = std::fabs(a[i] - b[i]);
			if(dist > maxDist) maxDist = dist;
		}
		return maxDist;
	}

	//@brief    : maxDelta adapted for rodrigues vectors
	//@param a  : first vector
	//@param b  : second vector
	//@param dum: dummy parameter to match signature to maxDelta
	//@return   : ~L1 distance
	template <typename Real> Real roDelta(Real const * const a, Real const * const b, const size_t dum) {
		const Real dW = std::fabs( std::fabs( std::atan(a[3]) ) - std::fabs( std::atan(b[3]) ) ) * 2;//angular difference
		return std::max(maxDelta(a, b, 3), dW);//difference is max of angular and axis differences
	}

	//@brief    : maxDelta adapted for euler angles vectors
	//@param a  : first vector
	//@param b  : second vector
	//@param dum: dummy parameter to match signature to maxDelta
	//@return   : ~L1 distance
	template <typename Real> Real euDelta(Real const * const a, Real const * const b, const size_t dum) {
		Real qa[4], qb[4];
		eu2qu(a, qa);
		eu2qu(b, qb);
		return maxDelta(qa, qb, 4);
	}

	//@brief    : write a vector to an ostream
	//@param v  : vector
	//@param len: length of vector
	//@param os : ostream to write to
	//@return   : os
	template <typename Real> std::ostream& printVec(Real const * const v, const size_t len, std::ostream& os) {
		os << "(";
		os << v[0];
		for(size_t i = 1; i < len; i++) os << ", " << v[i];
		return os << ")";
	}

	//@brief    : test rotation conversions for 2 and 3 way self consistency, e.g. for 2 way make sure eu2ro, ro2eu give back the input
	//@param n  : euler grid density
	//@param do3: should 3 way tests be included (this can significantly increase computation time)
	//@param os : location to write errors
	//@return   : true / false if tests pass / fail
	template <typename Real> bool testRound(const size_t n, const bool do3, std::ostream& os) {
		std::string names[7] = {"eu", "om", "ax", "ro", "qu", "ho", "cu"};
		size_t      lens [7] = {  3,    9,    4,    4,    4,    3,    3 };

		//build a table of conversion functions
		typedef void (*conversionFunc)(Real const * const, Real * const);//signature for conversion function pointer
		conversionFunc conversion[7][7] = {
			//  2eu           2om           2ax           2ro           2qu           2ho           2cu   
			{        NULL, &eu2om<Real>, &eu2ax<Real>, &eu2ro<Real>, &eu2qu<Real>, &eu2ho<Real>, &eu2cu<Real>},// eu2
			{&om2eu<Real>,         NULL, &om2ax<Real>, &om2ro<Real>, &om2qu<Real>, &om2ho<Real>, &om2cu<Real>},// om2
			{&ax2eu<Real>, &ax2om<Real>,         NULL, &ax2ro<Real>, &ax2qu<Real>, &ax2ho<Real>, &ax2cu<Real>},// ax2
			{&ro2eu<Real>, &ro2om<Real>, &ro2ax<Real>,         NULL, &ro2qu<Real>, &ro2ho<Real>, &ro2cu<Real>},// ro2
			{&qu2eu<Real>, &qu2om<Real>, &qu2ax<Real>, &qu2ro<Real>,         NULL, &qu2ho<Real>, &qu2cu<Real>},// qu2
			{&ho2eu<Real>, &ho2om<Real>, &ho2ax<Real>, &ho2ro<Real>, &ho2qu<Real>,         NULL, &ho2cu<Real>},// ho2
			{&cu2eu<Real>, &cu2om<Real>, &cu2ax<Real>, &cu2ro<Real>, &cu2qu<Real>, &cu2ho<Real>,         NULL} // cu2
		};

		//build a table of comparison functions
		typedef Real (*comparisonFunc)(Real const * const, Real const * const, const size_t);//signature for comparison function pointer
		comparisonFunc comparison[7] = {
			&euDelta <Real>,
			&maxDelta<Real>,
			&maxDelta<Real>,
			&roDelta <Real>,
			&maxDelta<Real>,
			&maxDelta<Real>,
			&maxDelta<Real>
		};

		//build an evenly spaced list of euler angles
		const Real pi = static_cast<Real>(3.1415926535897932384626433832795L);
		std::vector<Real> phi((n - 1) * 2 + 1);//[0,2pi]
		std::iota(phi.begin(), phi.end(), Real(0));
		std::transform(phi.begin(), phi.end(), phi.begin(), [&pi, &n](const Real i){return pi * i / (n-1);});
		std::vector<Real> theta(phi.begin(), phi.begin() + n);////[0,pi]

		//2 way tests
		Real maxDiff = Real(0);
		size_t maxI = 0, maxJ = 0, maxK = 0, maxM = 0, maxN = 0, maxP = 0;
		Real eu[3], x[9], y[9], z[9], final[9];
		for(size_t i = 0; i < 7; i++) {
			for(size_t j = 0; j < 7; j++) {
				if(i == j) continue;
				for(size_t m = 0; m < phi.size(); m++) {
					eu[0] = phi[m];
					for(size_t n = 0; n < theta.size(); n++) {
						eu[1] = theta[n];
						for(size_t p = 0; p < phi.size(); p++) {
							eu[2] = phi[p];

							//get base orientation from eulers
							if(i == 0) std::copy(eu, eu+3, x);
							else conversion[0][i](eu, x);

							//do conversion
							conversion[i][j](x, y);//x2y
							conversion[j][i](y, final);//y2x

							//compute error
							Real diff = comparison[i](x, final, lens[i]);
							if(diff > maxDiff) {
								maxDiff = diff;
								maxI = i;
								maxJ = j;
								maxM = m;
								maxN = n;
								maxP = p;
							}
						}
					}
				}
			}
		}

		os << "max diff pairwise = " << names[maxI] << "2" << names[maxJ] << "-" << names[maxJ] << "2" << names[maxI] << ": " << maxDiff << "\n";

		//select threshold to pass, I choose cube root due to trig operators
		//cbrt(epsilon) is ~0.005 and 6e-6 for float/double respectively 
		const Real eps2 = std::cbrt(std::numeric_limits<Real>::epsilon());
		if(maxDiff > eps2) {
			os << "outside limit (" << eps2 << ")\n";

			//get worst euler angle
			eu[0] = phi  [maxM];
			eu[1] = theta[maxN];
			eu[2] = phi  [maxP];

			//convert to base orientation
			if(maxI == 0) std::copy(eu, eu+3, x);
			else conversion[0][maxI](eu, x);

			//do conversion
			conversion[maxI][maxJ](x, y);//x2y
			conversion[maxJ][maxI](y, final);//y2x

			os << "\teu = ";
				printVec(eu   , 3         , os) << '\n';
			os << "\t" << names[maxI] << " = ";
				printVec(x    , lens[maxI], os) << '\n';
			os << "\t" << names[maxJ] << " = ";
				printVec(y    , lens[maxJ], os) << '\n';
			os << "\t" << names[maxI] << " = ";
				printVec(final, lens[maxI], os) << '\n';
			return false;
		}
		
		if(do3) {
			//3 way tests
			maxDiff = Real(0);
			maxI = maxJ = maxK = maxM = maxN = maxP = 0;
			for(size_t i = 0; i < 7; i++) {
				for(size_t j = 0; j < 7; j++) {
					if(i == j) continue;
						for(size_t k = 0; k < 7; k++) {
						if(j == k || k == i) continue;
						for(size_t m = 0; m < phi.size(); m++) {
							eu[0] = phi[m];
							for(size_t n = 0; n < theta.size(); n++) {
								eu[1] = theta[n];
								for(size_t p = 0; p < phi.size(); p++) {
									eu[2] = phi[p];
									//get base orientation from eulers
									if(i == 0) std::copy(eu, eu+3, x);
									else conversion[0][i](eu, x);

									//do conversion
									conversion[i][j](x, y);//x2y
									conversion[j][k](y, z);//y2z
									conversion[k][i](z, final);//z2x

									//compute error
									Real diff = comparison[i](x, final, lens[i]);
									if(diff > maxDiff) {
										maxDiff = diff;
										maxI = i;
										maxJ = j;
										maxK = k;
										maxM = m;
										maxN = n;
										maxP = p;
									}
								}
							}
						}
					}
				}
			}

			os << "max diff three way = " << names[maxI] << "2" << names[maxK] << "-" << names[maxK] << "2" << names[maxJ]  << "-" << names[maxJ] << "2" << names[maxI] << ": " << maxDiff << "\n";

			//select threshold to pass, I choose cube root due to trig operators
			//cbrt(epsilon) is ~0.005 and 6e-6 for float/double respectively 
			const Real eps3 = std::cbrt(std::numeric_limits<Real>::epsilon());
			if(maxDiff > eps3) {
				os << "outside limit (" << eps3 << ")\n";

				//get worst euler angle
				eu[0] = phi  [maxM];
				eu[1] = theta[maxN];
				eu[2] = phi  [maxP];

				//convert to base orientation
				if(maxI == 0) std::copy(eu, eu+3, x);
				else conversion[0][maxI](eu, x);

				//do conversion
				conversion[maxI][maxJ](x, y);//x2y
				conversion[maxJ][maxK](y, z);//y2z
				conversion[maxK][maxI](z, final);//z2x

				os << "\teu = ";
					printVec(eu   , 3         , os) << '\n';
				os << "\t" << names[maxI] << " = ";
					printVec(x    , lens[maxI], os) << '\n';
				os << "\t" << names[maxJ] << " = ";
					printVec(y    , lens[maxJ], os) << '\n';
				os << "\t" << names[maxK] << " = ";
					printVec(z    , lens[maxK], os) << '\n';
				os << "\t" << names[maxI] << " = ";
					printVec(final, lens[maxI], os) << '\n';
				return false;
			}
		}

		//finally make sure that zyz vs zxz is correct
		const Real epsZYZ = std::numeric_limits<Real>::epsilon() * 10;
		for(size_t m = 0; m < phi.size(); m++) {
			eu[0] = phi[m];
			for(size_t n = 0; n < theta.size(); n++) {
				eu[1] = theta[n];
				for(size_t p = 0; p < phi.size(); p++) {
					eu[2] = phi[p];
					Real zyz[3] = {eu[0] - pi / 2, eu[1], eu[2] + pi / 2};

					//check zyz2qu
					eu2qu (eu , x);
					zyz2qu(zyz, y);

					Real diff = maxDelta<Real>(x, y, 4);
					if(diff > epsZYZ) {
						os << "zyz2qu inconsistent with eu2qu({alpha + pi/2, beta, gamma - pi/2}) for:\n";
						os << "\teu   : " << eu[0] << ' ' << eu[1] << ' ' << eu[2] << '\n';
						os << "\teu   : " << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << '\n';
						os << "\teu   : " << y[0] << ' ' << y[1] << ' ' << y[2] << ' ' << y[3] << '\n';
						os << "\tdelta: " << diff << '\n';
						return false;
					}

					//check qu2zyz
					qu2zyz(x, y);
					y[0] += pi / 2;
					y[2] -= pi / 2;
					diff = euDelta(eu, y, 3);
					if(diff > epsZYZ) {
						os << "qu2zyz inconsistent with qu2eu(){alpha - pi/2, beta, gamma + pi/2}) for:\n";
						os << "\teu   : " << eu[0] << ' ' << eu[1] << ' ' << eu[2] << '\n';
						os << "\teu   : " << x[0] << ' ' << x[1] << ' ' << x[2] << ' ' << x[3] << '\n';
						os << "\teu   : " << y[0] << ' ' << y[1] << ' ' << y[2] << ' ' << y[3] << '\n';
						os << "\tdelta: " << diff << '\n';
						return false;
					}
				}
			}
		}

		//if we made it this far everything passed
		return true;
	}

	//@brief : run all rotation tests
	//@return: true / false if tests pass / fail
	template <typename Real> bool testRots(std::ostream& os) {
		return testRound<Real>(25, false, os) && //high grid density for 2 way
		       testRound<Real>(15, true , os);   //lower grid density for 3 way
	}
}
