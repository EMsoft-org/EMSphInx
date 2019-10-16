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

namespace color {
	//@brief    : test color space conversions for 2 and 3 way self consistency, e.g. for 2 way make sure rgb2hsv, hvs2rgb gives back the input
	//@param n  : RGB grid density
	//@param do3: should 3 way tests be included (this can significantly increase computation time)
	//@param os : location to write errors
	//@return   : true / false if tests pass / fail
	template <typename Real> bool testRound(const size_t n, const bool do3, std::ostream& os);

	//@brief : run all color space tests (2 and 3 way for different grid sizes)
	//@return: true / false if tests pass / fail
	template <typename Real> bool testRots(std::ostream& os);
}

int main() {
	std::ostream& os = std::cout;
	return color::testRots<float>(os) && color::testRots<double>(os) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#include <vector>
#include <cmath>
#include <functional>
#include <string>
#include <algorithm>//transform
#include <numeric>//iota

#include "util/colorspace.hpp"

namespace color {

	//@brief    : compute maximum of the element wise distance between two 3d vectors
	//@param a  : first vector
	//@param b  : second vector
	//@return   : L1 distance
	template <typename Real> Real maxDelta(Real const * const a, Real const * const b) {
		Real maxDist = 0;
		for(size_t i = 0; i < 3; i++) {
			const Real dist = std::fabs(a[i] - b[i]);
			if(dist > maxDist) maxDist = dist;
		}
		return maxDist;
	}

	//@brief    : special comparison for hsl
	//@param a  : first hsl
	//@param b  : second hsl
	//@return   : L1 distance
	template <typename Real> Real compHsl(Real const * const a, Real const * const b) {
		//compute max chroma
		const Real cA = std::min(a[2], Real(1) - a[2]);
		const Real cB = std::min(b[2], Real(1) - b[2]);

		//start by comparing hue periodically
		Real maxDist = std::fabs(a[0] - b[0]);
		maxDist = std::min(maxDist, Real(1) - maxDist);
		maxDist *= std::max(a[1], b[1]) * std::max(cA, cB);//normalize by max saturation and chroma

		//now compare chroma normalized saturation
		maxDist = std::max(maxDist, std::fabs(a[1]*cA - b[1]*cB));

		//finally compare lightness
		return std::max(maxDist, std::fabs(a[2] - b[2]));
	}

	//@brief    : special comparison for hsv
	//@param a  : first hsv
	//@param b  : second hsv
	//@return   : L1 distance
	template <typename Real> Real compHsv(Real const * const a, Real const * const b) {
		//compute max chroma
		const Real cA = a[2];
		const Real cB = b[2];

		//start by comparing hue periodically
		Real maxDist = std::fabs(a[0] - b[0]);
		maxDist = std::min(maxDist, Real(1) - maxDist);
		maxDist *= std::max(a[1], b[1]);//normalize by max saturation

		//now compare chroma normalized saturation
		maxDist = std::max(maxDist, std::fabs(a[1]*cA - b[1]*cB));

		//finally compare value
		return std::max(maxDist, std::fabs(a[2] - b[2]));
	}

	//@brief    : test color space conversions for 2 and 3 way self consistency, e.g. for 2 way make sure rgb2hsv, hvs2rgb gives back the input
	//@param num: RGB grid density
	//@param do3: should 3 way tests be included (this can significantly increase computation time)
	//@param os : location to write errors
	//@return   : true / false if tests pass / fail
	template <typename Real> bool testRound(const size_t num, const bool do3, std::ostream& os) {
		std::string names[6] = {"rgb", "xyz", "luv", "lab", "hsv", "hsl"};

		//build illumination dependant pointers
		typedef std::function<void(Real const * const, Real * const)> conversionFunc;//signature for conversion function pointer
		Real const * ill = NULL;//use default illuminant
		conversionFunc rgbXluv = std::bind(rgb2luv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc rgbXlab = std::bind(rgb2lab<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc xyzXluv = std::bind(xyz2luv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc xyzXlab = std::bind(xyz2lab<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc luvXrgb = std::bind(luv2rgb<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc luvXxyz = std::bind(luv2xyz<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc luvXlab = std::bind(luv2lab<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc luvXhsv = std::bind(luv2hsv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc luvXhsl = std::bind(luv2hsl<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc labXrgb = std::bind(lab2rgb<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc labXxyz = std::bind(lab2xyz<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc labXluv = std::bind(lab2luv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc labXhsv = std::bind(lab2hsv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc labXhsl = std::bind(lab2hsl<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc hsvXluv = std::bind(hsv2luv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc hsvXlab = std::bind(hsv2lab<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc hslXluv = std::bind(hsl2luv<Real>, std::placeholders::_1, std::placeholders::_2, ill);
		conversionFunc hslXlab = std::bind(hsl2lab<Real>, std::placeholders::_1, std::placeholders::_2, ill);

		//build a table of conversion functions
		conversionFunc conversion[6][6] = {
			//   2rgb            2xyz            2luv            2lab            2hsv            2hsl
			{          NULL, &rgb2xyz<Real>,  rgbXluv      ,  rgbXlab      , &rgb2hsv<Real>, &rgb2hsl<Real>},// rgb2
			{&xyz2rgb<Real>,           NULL,  xyzXluv      ,  xyzXlab      , &xyz2hsv<Real>, &xyz2hsl<Real>},// xyz2
			{ luvXrgb      ,  luvXxyz      ,           NULL,  luvXlab      ,  luvXhsv      ,  luvXhsl      },// luv2
			{ labXrgb      ,  labXxyz      ,  labXluv      ,           NULL,  labXhsv      ,  labXhsl      },// lab2
			{&hsv2rgb<Real>, &hsv2xyz<Real>,  hsvXluv      ,  hsvXlab      ,           NULL, &hsv2hsl<Real>},// hsv2
			{&hsl2rgb<Real>, &hsl2xyz<Real>,  hslXluv      ,  hslXlab      , &hsl2hsv<Real>,           NULL},// hsl2
		};

		//build a table of comparison functions
		typedef Real (*comparisonFunc)(Real const * const, Real const * const);//signature for comparison function pointer
		comparisonFunc comparison[6] = {
			&maxDelta<Real>,
			&maxDelta<Real>,
			&maxDelta<Real>,
			&maxDelta<Real>,
			&compHsv <Real>,
			&compHsl <Real>,
		};

		//build a grid of RGB values
		std::vector<Real> lin(num);//[0,1]
		std::iota(lin.begin(), lin.end(), Real(0));
		std::for_each(lin.begin(), lin.end(), [num](Real& v){v /= (num-1);});

		//2 way tests
		Real maxDiff = Real(0);
		size_t maxI = 0, maxJ = 0, maxK = 0, maxM = 0, maxN = 0, maxP = 0;
		Real rgb[3], x[3], y[3], z[3], final[3];
		for(size_t i = 0; i < 6; i++) {
			for(size_t j = 0; j < 6; j++) {
				if(i == j) continue;
				for(size_t m = 0; m < num; m++) {
					rgb[0] = lin[m];
					for(size_t n = 0; n < num; n++) {
						rgb[1] = lin[n];
						for(size_t p = 0; p < num; p++) {
							rgb[2] = lin[p];

							//get base orientation from eulers
							if(i == 0) std::copy(rgb, rgb+3, x);
							else conversion[0][i](rgb, x);

							//do conversion
							conversion[i][j](x, y);//x2y
							conversion[j][i](y, final);//y2x

							//compute error
							Real diff = comparison[i](x, final);
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
			rgb[0] = lin[maxM];
			rgb[1] = lin[maxN];
			rgb[2] = lin[maxP];

			//convert to base orientation
			if(maxI == 0) std::copy(rgb, rgb+3, x);
			else conversion[0][maxI](rgb, x);

			//do conversion
			conversion[maxI][maxJ](x, y);//x2y
			conversion[maxJ][maxI](y, final);//y2x

			os << "\trgb = ("                   << rgb  [0] << '\t' << rgb  [1] << '\t' << rgb  [2] << ")\n";
			os << "\t" << names[maxI] << " = (" << x    [0] << '\t' << x    [1] << '\t' << x    [2] << ")\n";
			os << "\t" << names[maxJ] << " = (" << y    [0] << '\t' << y    [1] << '\t' << y    [2] << ")\n";
			os << "\t" << names[maxI] << " = (" << final[0] << '\t' << final[1] << '\t' << final[2] << ")\n";
			os << comparison[maxI](x, final) << '\n';
			return false;
		}

		if(do3) {
			//3 way tests
			maxDiff = Real(0);
			maxI = maxJ = maxK = maxM = maxN = maxP = 0;
			for(size_t i = 0; i < 6; i++) {
				for(size_t j = 0; j < 6; j++) {
					if(i == j) continue;
						for(size_t k = 0; k < 6; k++) {
						if(j == k || k == i) continue;
						for(size_t m = 0; m < num; m++) {
							rgb[0] = lin[m];
							for(size_t n = 0; n < num; n++) {
								rgb[1] = lin[n];
								for(size_t p = 0; p < num; p++) {
									rgb[2] = lin[p];
									//get base orientation from rgblers
									if(i == 0) std::copy(rgb, rgb+3, x);
									else conversion[0][i](rgb, x);

									//do conversion
									conversion[i][j](x, y);//x2y
									conversion[j][k](y, z);//y2z
									conversion[k][i](z, final);//z2x

									//compute error
									Real diff = comparison[i](x, final);
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
				rgb[0] = lin[maxM];
				rgb[1] = lin[maxN];
				rgb[2] = lin[maxP];

				//convert to base orientation
				if(maxI == 0) std::copy(rgb, rgb+3, x);
				else conversion[0][maxI](rgb, x);

				//do conversion
				conversion[maxI][maxJ](x, y);//x2y
				conversion[maxJ][maxK](y, z);//y2z
				conversion[maxK][maxI](z, final);//z2x

				os << "\trgb = ("                   << rgb  [0] << '\t' << rgb  [1] << '\t' << rgb  [2] << ")\n";
				os << "\t" << names[maxI] << " = (" << x    [0] << '\t' << x    [1] << '\t' << x    [2] << ")\n";
				os << "\t" << names[maxJ] << " = (" << y    [0] << '\t' << y    [1] << '\t' << y    [2] << ")\n";
				os << "\t" << names[maxK] << " = (" << z    [0] << '\t' << z    [1] << '\t' << z    [2] << ")\n";
				os << "\t" << names[maxI] << " = (" << final[0] << '\t' << final[1] << '\t' << final[2] << ")\n";
				return false;
			}
		}
		//if we made it this far everything passed
		return true;
	}

	//@brief : run all rotation tests
	//@return: true / false if tests pass / fail
	template <typename Real> bool testRots(std::ostream& os) {
		return testRound<Real>(51, false, os) && //high grid density for 2 way
		       testRound<Real>(31, true , os);   //lower grid density for 3 way
	}
}
