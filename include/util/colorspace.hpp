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


#ifndef _COLORSPACE_H_
#define _COLORSPACE_H_

namespace color {
	//@brief: color space conversion functions abc2ijk where abc/ijk are two of
	// -rgb, {r , g , b }: standard red, green, blue (sRGB)
	// -xyz, {X , Y , Z }: CIE 1931 XYZ color space ['master' (original) perceptually uniform color space]
	// -luv, {L*, u*, v*}: 1976 CIELUV color space [perceptually uniform space for computer displays]
	// -lab, {L*, a*, b*}: 1976 CIELab color space [perceptually uniform space for print]
	// -hsv, {h , s , v }: hue, saturation, value (cylindrical)
	// -hsl, {h , s , l }: hue, saturation, lightness (cylindrical)
	//@note: rgb, hsv, and hsl are restricted to the range [0,1] with values outside representing imaginary colors
	//@note: the range [0,1] is used for hue in hsv/hsl (not [0,2*pi] or [0,360])
	//@note: all conversions are available using the shortest path in network below
	//       .->hsv <---> rgb <---> xyz <---> luv
	//       `->hsl <-'                   `-> lab
	//       therefore the following direct transformations are implemented
	//        -rgb2xyz, rgb2hsv, rgb2hsl
	//        -xyz2rgb, xyz2luv, xyz2lab
	//        -luv2xyz, lab2xyz
	//        -hsv2rgb, hsv2hsl
	//        -hsl2rgb, hsl2hsv
	//       with indirect transformations requiring a combination of the above
	//@param abc: values to convert from abc color space
	//@param ijk: location to write values converted to ijk color space (can be the same as parameter abc)
	//@param ill: standard illuminant as xyz (only required for conversions involving xyz<->luv or xyz<->lab, defaults to CIE illuminant D65 for 2 degree observer)
	//@return: true/false if the values fall outside the ijk gamut for conversions to that pass through xyz2rgb (void for others)
	template <typename T> void rgb2xyz(T const * const rgb, T * const xyz                      );
	template <typename T> void rgb2luv(T const * const rgb, T * const luv, T const * ill = NULL);
	template <typename T> void rgb2lab(T const * const rgb, T * const lab, T const * ill = NULL);
	template <typename T> void rgb2hsv(T const * const rgb, T * const hsv                      );
	template <typename T> void rgb2hsl(T const * const rgb, T * const hsl                      );

	template <typename T> bool xyz2rgb(T const * const xyz, T * const rgb                      );
	template <typename T> void xyz2luv(T const * const xyz, T * const luv, T const * ill = NULL);
	template <typename T> void xyz2lab(T const * const xyz, T * const lab, T const * ill = NULL);
	template <typename T> bool xyz2hsv(T const * const xyz, T * const hsv                      );
	template <typename T> bool xyz2hsl(T const * const xyz, T * const hsl                      );

	template <typename T> bool luv2rgb(T const * const luv, T * const rgb, T const * ill = NULL);
	template <typename T> void luv2xyz(T const * const luv, T * const xyz, T const * ill = NULL);
	template <typename T> void luv2lab(T const * const luv, T * const lab, T const * ill = NULL);
	template <typename T> bool luv2hsv(T const * const luv, T * const hsv, T const * ill = NULL);
	template <typename T> bool luv2hsl(T const * const luv, T * const hsl, T const * ill = NULL);

	template <typename T> bool lab2rgb(T const * const lab, T * const rgb, T const * ill = NULL);
	template <typename T> void lab2xyz(T const * const lab, T * const xyz, T const * ill = NULL);
	template <typename T> void lab2luv(T const * const lab, T * const luv, T const * ill = NULL);
	template <typename T> bool lab2hsv(T const * const lab, T * const hsv, T const * ill = NULL);
	template <typename T> bool lab2hsl(T const * const lab, T * const hsl, T const * ill = NULL);

	template <typename T> void hsv2rgb(T const * const hsv, T * const rgb                      );
	template <typename T> void hsv2xyz(T const * const hsv, T * const xyz                      );
	template <typename T> void hsv2luv(T const * const hsv, T * const luv, T const * ill = NULL);
	template <typename T> void hsv2lab(T const * const hsv, T * const lab, T const * ill = NULL);
	template <typename T> void hsv2hsl(T const * const hsv, T * const hsl                      );

	template <typename T> void hsl2rgb(T const * const hsl, T * const rgb                      );
	template <typename T> void hsl2xyz(T const * const hsl, T * const xyz                      );
	template <typename T> void hsl2luv(T const * const hsl, T * const luv, T const * ill = NULL);
	template <typename T> void hsl2lab(T const * const hsl, T * const lab, T const * ill = NULL);
	template <typename T> void hsl2hsv(T const * const hsl, T * const hsv                      );
}

////////////////////////////////////////////////////////////////////////
//                       Implementation Details                       //
////////////////////////////////////////////////////////////////////////

#include <array>
#include <cmath>
#include <numeric>
#include <algorithm>

namespace color {
	namespace detail {
		//@brief    : invert a 3x3 matrix analytically
		//@param mat: matrix to invert in row major order
		template <typename T> std::array<T, 9> inv3x3(const std::array<T, 9>& mat) {
			const T det = mat[3*0+0] * mat[3*1+1] * mat[3*2+2] + mat[3*0+1] * mat[3*1+2] * mat[3*2+0] + mat[3*0+2] * mat[3*1+0] * mat[3*2+1]
						-(mat[3*0+0] * mat[3*1+2] * mat[3*2+1] + mat[3*0+1] * mat[3*1+0] * mat[3*2+2] + mat[3*0+2] * mat[3*1+1] * mat[3*2+0]);
			return std::array<T, 9> {
				(mat[3*1+1]*mat[3*2+2] - mat[3*1+2]*mat[3*2+1]) / det, (mat[3*0+2]*mat[3*2+1] - mat[3*0+1]*mat[3*2+2]) / det, (mat[3*0+1]*mat[3*1+2] - mat[3*0+2]*mat[3*1+1]) / det,
				(mat[3*1+2]*mat[3*2+0] - mat[3*1+0]*mat[3*2+2]) / det, (mat[3*0+0]*mat[3*2+2] - mat[3*0+2]*mat[3*2+0]) / det, (mat[3*0+2]*mat[3*1+0] - mat[3*0+0]*mat[3*1+2]) / det,
				(mat[3*1+0]*mat[3*2+1] - mat[3*1+1]*mat[3*2+0]) / det, (mat[3*0+1]*mat[3*2+0] - mat[3*0+0]*mat[3*2+1]) / det, (mat[3*0+0]*mat[3*1+1] - mat[3*0+1]*mat[3*1+0]) / det
			};
		}

		//@brief    : compute a matrix to convert from XYZ --> rgb
		//@param rgb: chromaticity of {red point, green point, blue point} as XYZ
		//@param w  : chromaticity of white point as XYZ
		template <typename T>
		std::array<T, 9> rgbMat(const T rgb[3][3], T const w[3]) {
			//build and invert 3x3 matrices to solve for rows of conversion matrix
			//matrix * (r, g, b) = {x, 0, 0}, {0, x, 0}, {0, x, 0} and matrix^-1 * {1,1,1} = w
			const T W[3] = {w[0] / w[1], T(1), w[2] / w[1]};
			const std::array<T, 9> invR = inv3x3<T>({  W[0]   ,   W[1]   ,   W[2]   , rgb[1][0], rgb[1][1], rgb[1][2], rgb[2][0], rgb[2][1], rgb[2][2]});
			const std::array<T, 9> invG = inv3x3<T>({rgb[0][0], rgb[0][1], rgb[0][2],   W[0]   ,   W[1]   ,   W[2]   , rgb[2][0], rgb[2][1], rgb[2][2]});
			const std::array<T, 9> invB = inv3x3<T>({rgb[0][0], rgb[0][1], rgb[0][2], rgb[1][0], rgb[1][1], rgb[1][2],   W[0]   ,   W[1]   ,   W[2]   });
			return {invR[3*0+0], invR[3*1+0], invR[3*2+0], invG[3*0+1], invG[3*1+1], invG[3*2+1], invB[3*0+2], invB[3*1+2], invB[3*2+2]};//assemble matrix
		}

		//@brief    : convert from hue, chroma, and minimum to rgb
		//@param h  : hue
		//@param c  : chroma
		//@param m  : minimum
		//@param rgb: location to write rgb (red, green, blue) values
		template <typename T> void hcm2rgb(const T h, const T c, const T m, T * const rgb) {
			const T h6 = h * 6;
			const T x = c * (T(1) - std::fabs(std::fmod(h6, T(2)) - T(1)));
			switch((size_t)h6) {
				case 6://intentional fall through
				case 0: rgb[0] = c+m; rgb[1] = x+m; rgb[2] =   m; return;
				case 1: rgb[0] = x+m; rgb[1] = c+m; rgb[2] =   m; return;
				case 2: rgb[0] =   m; rgb[1] = c+m; rgb[2] = x+m; return;
				case 3: rgb[0] =   m; rgb[1] = x+m; rgb[2] = c+m; return;
				case 4: rgb[0] = x+m; rgb[1] =   m; rgb[2] = c+m; return;
				case 5: rgb[0] = c+m; rgb[1] =   m; rgb[2] = x+m; return;
			}
		}

		//constants for standard illuminants and common rgb gamuts
		template <typename T>
		struct Standards {
			//standard illuminants as XYZ
			static const T A_2  [3], A_10  [3];//A for 2 and 10 degree observer
			static const T B_2  [3], B_10  [3];//B for 2 and 10 degree observer
			static const T C_2  [3], C_10  [3];//C for 2 and 10 degree observer
			static const T D50_2[3], D50_10[3];//D50 for 2 and 10 degree observer
			static const T D55_2[3], D55_10[3];//D55 for 2 and 10 degree observer
			static const T D65_2[3], D65_10[3];//D65 for 2 and 10 degree observer
			static const T D75_2[3], D75_10[3];//D75 for 2 and 10 degree observer
			static const T E    [3];

			//chromaticities of RGB standards ({{red point}, {green point}, {blue point}} as XYZ)
			//all use standard illuminant D65 and gamma 2.2 unless otherwise noted
			static const T sRGB    [3][3];//standard RGB (custom gamma)
			static const T cieRGB  [3][3];//CIE (1931) RGB (standard illuminant E)
			static const T appleRGB[3][3];//Apple RGB (1.8 gamma)
			static const T palRGB  [3][3];//PAL / SECAM RGB
			static const T ntscRGB [3][3];//NTSC / SMPTE C RGB
			static const T adobeRGB[3][3];//Adobe RGB (1998)

			//sRGB custom gamma constants
			static const T sA ;   //deviation of gamma correction coefficient from 1
			static const T sGamma;//gamma exponent
			static const T sPhi;  //scale factor for linear gamma region
			static const T sK0;   //cutoff for linear gamma correction in inverse direction

			//sRGB <--> XYZ conversion matrices
			static const std::array<T, 9> sRGBmat   ;//matrix to convert from XYZ  --> sRGB
			static const std::array<T, 9> sRGBmatInv;//matrix to convert from sRGB --> XYZ
		};

		//standard illuminants as xyz (normalized XYZ)
		template <typename T> const T Standards<T>::A_2   [3] = {T(0.44757), T(0.40745), T(0.14498)};
		template <typename T> const T Standards<T>::A_10  [3] = {T(0.45117), T(0.40594), T(0.14289)};
		template <typename T> const T Standards<T>::B_2   [3] = {T(0.34842), T(0.35161), T(0.29997)};
		template <typename T> const T Standards<T>::B_10  [3] = {T(0.34980), T(0.35270), T(0.29750)};
		template <typename T> const T Standards<T>::C_2   [3] = {T(0.31006), T(0.31616), T(0.37378)};
		template <typename T> const T Standards<T>::C_10  [3] = {T(0.31039), T(0.31905), T(0.37056)};
		template <typename T> const T Standards<T>::D50_2 [3] = {T(0.34567), T(0.35850), T(0.29583)};
		template <typename T> const T Standards<T>::D50_10[3] = {T(0.34773), T(0.35952), T(0.29275)};
		template <typename T> const T Standards<T>::D55_2 [3] = {T(0.33242), T(0.34743), T(0.32015)};
		template <typename T> const T Standards<T>::D55_10[3] = {T(0.33411), T(0.34877), T(0.31712)};
		template <typename T> const T Standards<T>::D65_2 [3] = {T(0.31271), T(0.32902), T(0.35827)};
		template <typename T> const T Standards<T>::D65_10[3] = {T(0.31382), T(0.33100), T(0.35518)};
		template <typename T> const T Standards<T>::D75_2 [3] = {T(0.29902), T(0.31485), T(0.38613)};
		template <typename T> const T Standards<T>::D75_10[3] = {T(0.29968), T(0.31740), T(0.38292)};
		template <typename T> const T Standards<T>::E     [3] = {T(1)/T(3) , T(1)/T(3) , T(1)/T(3) };

		//RGB chromaticities as xyz (normalized XYZ)
		template <typename T> const T Standards<T>::sRGB    [3][3] = {T(0.6400), T(0.3300), T(0.0300),
																	  T(0.3000), T(0.6000), T(0.1000),
																	  T(0.1500), T(0.0600), T(0.7900)};
		template <typename T> const T Standards<T>::cieRGB  [3][3] = {T(0.7347), T(0.2653), T(0.0000),
																	  T(0.2738), T(0.7174), T(0.0088),
																	  T(0.1666), T(0.0089), T(0.8245)};
		template <typename T> const T Standards<T>::appleRGB[3][3] = {T(0.6250), T(0.3400), T(0.0350),
																	  T(0.2800), T(0.5950), T(0.1250),
																	  T(0.1550), T(0.0700), T(0.7750)};
		template <typename T> const T Standards<T>::adobeRGB[3][3] = {T(0.6400), T(0.3300), T(0.0300),
																	  T(0.2100), T(0.7100), T(0.0800),
																	  T(0.1500), T(0.0600), T(0.7900)};
		template <typename T> const T Standards<T>::palRGB  [3][3] = {T(0.6400), T(0.3300), T(0.0300),
																	  T(0.2900), T(0.6000), T(0.1100),
																	  T(0.1500), T(0.0600), T(0.7900)};
		template <typename T> const T Standards<T>::ntscRGB [3][3] = {T(0.6300), T(0.3400), T(0.0300),
																	  T(0.3100), T(0.5950), T(0.0950),
																	  T(0.1550), T(0.0700), T(0.7750)};

		//sRGB gamma correction constants
		template <typename T> const T Standards<T>::sA     = T(0.055  );
		template <typename T> const T Standards<T>::sGamma = T(2.4    );
		template <typename T> const T Standards<T>::sPhi   = T(12.92  );
		template <typename T> const T Standards<T>::sK0    = T(0.04045);

		//sRGB conversion matrices
		template <typename T> const std::array<T, 9> Standards<T>::sRGBmat    = detail::rgbMat(Standards<T>::sRGB, Standards<T>::D65_2);
		template <typename T> const std::array<T, 9> Standards<T>::sRGBmatInv = detail::inv3x3(Standards<T>::sRGBmat);
	}

	////////////////////////////////////////////
	//  implementation of direct conversions  //
	////////////////////////////////////////////

	//@brief    : convert from XYZ to sRGB
	//@param xyz: XYZ (X, Y, Z) values to convert
	//@param rgb: location to write sRGB (red, green, blue) values
	//@return   : true/false if xyz falls outside/inside the sRGB color gamut
	template <typename T> bool xyz2rgb(T const * const xyz, T * const rgb) {
		using namespace detail;
		static const T gammaInv = T(1) / Standards<T>::sGamma;
		static const T k0Inv    = Standards<T>::sK0 / Standards<T>::sPhi;
		static const T a1       = T(1) + Standards<T>::sA;
		T work[3];
		for(size_t i = 0; i < 3; i++) work[i] = std::inner_product(xyz, xyz+3, Standards<T>::sRGBmat.begin() + 3*i, T(0));//XYZ -> linear rgb
		for(size_t i = 0; i < 3; i++) rgb[i] = work[i] <= k0Inv ? work[i] * Standards<T>::sPhi : a1 * std::pow(work[i], gammaInv) - Standards<T>::sA;//gamma correction

		//check if this value is outside the sRGB color gamut
		bool clamped = false;
		for(size_t i = 0; i < 3; i++) {
			if(std::signbit(rgb[i])) {
				rgb[i] = T(0);
				clamped = true;
			} else if(rgb[i] > T(1)) {
				rgb[i] = T(1);
				clamped = true;
			}
		}
		return clamped;
	}

	//@brief    : convert from sRGB to XYZ
	//@param rgb: sRGB (red, green, blue) values to convert
	//@param xyz: location to write XYZ (X, Y, Z) values
	template <typename T> void rgb2xyz(T const * const rgb, T * const xyz) {
		using namespace detail;
		static const T a1 = T(1) + Standards<T>::sA;
		T work[3];
		for(size_t i = 0; i < 3; i++) xyz[i] = rgb[i] <= Standards<T>::sK0 ? rgb[i] / Standards<T>::sPhi : std::pow((rgb[i]+Standards<T>::sA) / a1, Standards<T>::sGamma);//gamma correction
		for(size_t i = 0; i < 3; i++) work[i] = std::inner_product(xyz, xyz+3, Standards<T>::sRGBmatInv.begin() + 3*i, T(0));//XYZ -> linear rgb
		std::copy(work, work+3, xyz);
	}

	//@brief    : convert from XYZ to Lab
	//@param xyz: XYZ (X, Y, Z) values to convert
	//@param lab: location to write Lab (L*, a*, b*) values
	//@param ill: Lab illuminant as XYZ (or NULL to use illuminant D65 for a 2 degree observer)
	template <typename T> void xyz2lab(T const * const xyz, T * const lab, T const * ill) {
		//initialize constants once
		static const T delta = T(6) / 29;
		static const T d2 = delta * delta * 3;
		static const T d3 = delta * delta * delta;
		static const T k = T(4) / 29;

		//compute f(i/i_N)
		T fXYZ[3];
		T const * const illum = (NULL != ill) ? ill : detail::Standards<T>::D65_2;
		for(size_t i = 0; i < 3; i++) {
			const T t = xyz[i] / (illum[i] / illum[1]);//normalize with respect to white point
			fXYZ[i] = t > d3 ? std::cbrt(t) : t / d2 + k;//apply nonlinear scaling
		}

		//change basis and scale
		lab[0] = fXYZ[1] * 116 - 16;//L*
		lab[1] = (fXYZ[0] - fXYZ[1]) * 500;//a*
		lab[2] = (fXYZ[1] - fXYZ[2]) * 200;//b*
	}

	//@brief    : convert from Lab to xyz
	//@param lab: Lab (L*, a*, b*) values to convert
	//@param xyz: location to write XYZ (X, Y, Z) values
	//@param ill: Lab illuminant as XYZ (or NULL to use illuminant D65 for a 2 degree observer)
	template <typename T> void lab2xyz(T const * const lab, T * const xyz, T const * ill) {
		//initialize constants once
		static const T delta = T(6) / 29;
		static const T d2 = delta * delta * 3;
		static const T d3 = delta * delta * delta;
		static const T k = T(4) / 29;

		//change basis and scale
		const T Ln = (lab[0] + T(16)) / 116;
		xyz[0] = Ln + lab[1] / 500;
		xyz[1] = Ln;
		xyz[2] = Ln - lab[2] / 200;

		//compute inverse of f(i/i_N)
		T const * const illum = (NULL != ill) ? ill : detail::Standards<T>::D65_2;
		for(size_t i = 0; i < 3; i++) {
			xyz[i] = xyz[i] > delta ? xyz[i] * xyz[i] * xyz[i] : d2 * (xyz[i] - k);//remove nonlinear scaling
			xyz[i] *= illum[i] / illum[1];//remove white point normalization
		}
	}

	//@brief    : convert from XYZ to Luv
	//@param xyz: XYZ (X, Y, Z) values to convert
	//@param luv: location to write Luv (L*, u*, v*) values
	//@param ill: Lab illuminant as XYZ (or NULL to use illuminant D65 for a 2 degree observer)
	template <typename T> void xyz2luv(T const * const xyz, T * const luv, T const * ill) {
		//initialize constants once
		static const T d   = T(216) / 24389;//(6/29)^3
		static const T d_8 = T(8) / d;//(29/3)^3

		//compute normalized X, Y, and Z and denomentators and u/v
		T const * const illum = (NULL != ill) ? ill : detail::Standards<T>::D65_2;
		const T denn = (illum[0] + illum[1] * 15 + illum[2] * 3) / illum[1];
		const T den  =  xyz  [0] + xyz  [1] * 15 + xyz  [2] * 3;
		const bool zero = T(0) == den;
		const T u = (zero ? T(0) : xyz[0] / den) - (illum[0] / illum[1]) / denn;
		const T v = (zero ? T(0) : xyz[1] / den) - T(1)                 / denn;

		//compute Luv
		luv[0] = xyz[1] <= d ? xyz[1] * d_8 : std::cbrt(xyz[1]) * 116 - 16;
		luv[1] = luv[0] *  52 * u;
		luv[2] = luv[0] * 117 * v;
	}

	//@brief    : convert from Luv to xyz
	//@param luv: Luv (L*, u*, v*) values to convert
	//@param xyz: location to write XYZ (X, Y, Z) values
	//@param ill: Lab illuminant as XYZ (or NULL to use illuminant D65 for a 2 degree observer)
	template <typename T> void luv2xyz(T const * const luv, T * const xyz, T const * ill) {
		//handle L* = 0
		if(luv[0] == T(0)) {
			std::fill(xyz, xyz+3, T(0));
			return;
		}

		//compute u' and v'
		T const * const illum = (NULL != ill) ? ill : detail::Standards<T>::D65_2;
		const T denn = (illum[0] + illum[1] * 15 + illum[2] * 3) / illum[1];
		const T up = (luv[1] / 13 + luv[0] * (illum[0] / illum[1]) * 4 / denn) * 3;
		const T vp = (luv[2] / 13 + luv[0]                         * 9 / denn) * 4;
		const T Lp = (luv[0] + 16) / 116;

		//compute X, Y, and Z
		static const T d = T(27) / 24389;//(3/29)^3
		xyz[1] = luv[0] <= 8 ? luv[0] * d : Lp * Lp * Lp;
		xyz[2] = xyz[1] * (T(12) * luv[0] - up - vp * 5) / vp;
		xyz[0] = xyz[1] * (up * 3) / vp;
	}

	//@brief    : convert from hsv to rgb
	//@param hsv: hsv (hue, saturation, value) values to convert
	//@param rgb: location to write rgb (red, green, blue) values
	template <typename T> void hsv2rgb(T const * const hsv, T * const rgb) {
		const T c = hsv[1] * hsv[2];//chroma
		detail::hcm2rgb(hsv[0], c, hsv[2] - c  , rgb);
	}

	//@brief    : convert from hsl to rgb
	//@param hsl: hsl (hue, saturation, lightness) values to convert
	//@param rgb: location to write rgb (red, green, blue) values
	template <typename T> void hsl2rgb(T const * const hsl, T * const rgb) {
		const T c = (T(1) - std::fabs(hsl[2] * 2 - 1)) * hsl[1];//chroma
		detail::hcm2rgb(hsl[0], c, hsl[2] - c/2, rgb);
	}

	//@brief    : convert from hsl to hsv
	//@param hsl: hsl (hue, saturation, lightness) values to convert
	//@param hsv: location to write hsv (hue, saturation, value) values
	template <typename T> void hsl2hsv(T const * const hsl, T * const hsv) {
		const T s = hsl[1] * (hsl[2] < T(0.5) ? hsl[2] : T(1) - hsl[2]);
		hsv[0] = hsl[0];
		hsv[2] = s + hsl[2];
		hsv[1] = T(0) == s ? T(0) : T(2) * s / hsv[2];
	}

	//@brief    : convert from hsv to hsl
	//@param hsv: hsv (hue, saturation, lightness) values to convert
	//@param hsl: location to write hsl (hue, saturation, value) values
	template <typename T> void hsv2hsl(T const * const hsv, T * const hsl) {
		const T sv = hsv[1] * hsv[2];
		const T x = 2 * hsv[2] - sv;
		hsl[0] = hsv[0];
		hsl[2] = hsv[2] - sv / 2;
		hsl[1] = T(0) == sv ? T(0) : sv / (x < T(1) ? x : T(2) - x);
	}

	//@brief    : convert rgb to hsv
	//@param rgb: sRGB (red, green, blue) values to convert
	//@param hsv: location to write hsv (hue, saturation, value) values
	template <typename T> void rgb2hsv(T const * const rgb, T * const hsv) {
		auto minMax = std::minmax_element(rgb, rgb+3);
		const T vMax = *(minMax.second);
		if(T(0) == vMax) {//black
			std::fill(hsv, hsv + 3, T(0));
		} else {
			const T vMin = *(minMax.first );
			const T delta = vMax - vMin;
			if(T(0) == delta) {//gray
				hsv[0] = hsv[1] = T(0);
			} else {
				const size_t indMax = std::distance(rgb, minMax.second);
				hsv[0] = T(indMax) / 3 + (rgb[(indMax+1)%3] - rgb[(indMax+2)%3]) / (delta * 6);
				if(std::signbit(hsv[0])) hsv[0] += T(1);
				if(hsv[0] == T(1)) hsv[0] = T(0);// if instead of else if to handle -0
				hsv[1] = delta / vMax;
			}
			hsv[2] = vMax;
		}
	}

	//@brief    : convert rgb to hsl
	//@param rgb: sRGB (red, green, blue) values to convert
	//@param hsl: location to write hsl (hue, saturation, lightness) values
	template <typename T> void rgb2hsl(T const * const rgb, T * const hsl) {
		auto minMax = std::minmax_element(rgb, rgb+3);
		const T vMax = *(minMax.second);
		if(T(0) == vMax) {//black
			std::fill(hsl, hsl + 3, T(0));
		} else {
			const T vMin = *(minMax.first );
			const T delta = vMax - vMin;
			const T sigma = vMax + vMin;
			if(T(0) == delta) {//gray
				hsl[0] = hsl[1] = T(0);
			} else {
				const size_t indMax = std::distance(rgb, minMax.second);
				hsl[0] = T(indMax) / 3 + (rgb[(indMax+1)%3] - rgb[(indMax+2)%3]) / (delta * 6);//hue
				if(std::signbit(hsl[0])) hsl[0] += T(1);
				hsl[1] = delta / (sigma < T(1) ? sigma: T(2) - sigma);//saturation
			}
			hsl[2] = sigma / 2;//lightness
		}
	}

	////////////////////////////////////////////
	// implementation of indirect conversions //
	////////////////////////////////////////////

	//illuminant free cie spaces -> hsl/hsv (using cie spaces -> rgb)
	template <typename T> bool xyz2hsv(T const * const xyz, T * const hsv               ) {const bool b = xyz2rgb(xyz, hsv); rgb2hsv(hsv, hsv); return b;}//xyz->rgb->hsv
	template <typename T> bool xyz2hsl(T const * const xyz, T * const hsl               ) {const bool b = xyz2rgb(xyz, hsl); rgb2hsl(hsl, hsl); return b;}//xyz->rgb->hsl

	//illuminated cie spaces -> rgb/hsv/hsl (using cie spaces -> xyz -> rgb)
	template <typename T> bool luv2rgb(T const * const luv, T * const rgb, T const * ill) {luv2xyz(luv, rgb, ill); return xyz2rgb(rgb, rgb);}//luv->xyz->rgb
	template <typename T> bool luv2hsv(T const * const luv, T * const hsv, T const * ill) {luv2xyz(luv, hsv, ill); return xyz2hsv(hsv, hsv);}//luv->xyz->rgb->hsv
	template <typename T> bool luv2hsl(T const * const luv, T * const hsl, T const * ill) {luv2xyz(luv, hsl, ill); return xyz2hsl(hsl, hsl);}//luv->xyz->rgb->hsl
	template <typename T> bool lab2rgb(T const * const lab, T * const rgb, T const * ill) {lab2xyz(lab, rgb, ill); return xyz2rgb(rgb, rgb);}//lab->xyz->rgb
	template <typename T> bool lab2hsv(T const * const lab, T * const hsv, T const * ill) {lab2xyz(lab, hsv, ill); return xyz2hsv(hsv, hsv);}//lab->xyz->rgb->hsv
	template <typename T> bool lab2hsl(T const * const lab, T * const hsl, T const * ill) {lab2xyz(lab, hsl, ill); return xyz2hsl(hsl, hsl);}//lab->xyz->rgb->hsl

	//within cie spaces (through xyz)
	template <typename T> void luv2lab(T const * const luv, T * const lab, T const * ill) {luv2xyz(luv, lab, ill); xyz2lab(lab, lab, ill);}//luv->xyz->lab
	template <typename T> void lab2luv(T const * const lab, T * const luv, T const * ill) {lab2xyz(lab, luv, ill); xyz2luv(luv, luv, ill);}//lab->xyz->luv

	//rgb/hsv/hsl to cie spaces (all go through rgb -> xyz)
	template <typename T> void rgb2luv(T const * const rgb, T * const luv, T const * ill) {                   rgb2xyz(rgb, luv); xyz2luv(luv, luv, ill);}//     rgb->xyz->luv
	template <typename T> void rgb2lab(T const * const rgb, T * const lab, T const * ill) {                   rgb2xyz(rgb, lab); xyz2lab(lab, lab, ill);}//     rgb->xyz->lab
	template <typename T> void hsv2xyz(T const * const hsv, T * const xyz               ) {hsv2rgb(hsv, xyz); rgb2xyz(xyz, xyz);                        }//hsv->rgb->xyz
	template <typename T> void hsv2luv(T const * const hsv, T * const luv, T const * ill) {hsv2rgb(hsv, luv); rgb2xyz(luv, luv); xyz2luv(luv, luv, ill);}//hsv->rgb->xyz->luv
	template <typename T> void hsv2lab(T const * const hsv, T * const lab, T const * ill) {hsv2rgb(hsv, lab); rgb2xyz(lab, lab); xyz2lab(lab, lab, ill);}//hsv->rgb->xyz->lab
	template <typename T> void hsl2xyz(T const * const hsl, T * const xyz               ) {hsl2rgb(hsl, xyz); rgb2xyz(xyz, xyz);                        }//hsl->rgb->xyz
	template <typename T> void hsl2luv(T const * const hsl, T * const luv, T const * ill) {hsl2rgb(hsl, luv); rgb2xyz(luv, luv); xyz2luv(luv, luv, ill);}//hsl->rgb->xyz->luv
	template <typename T> void hsl2lab(T const * const hsl, T * const lab, T const * ill) {hsl2rgb(hsl, lab); rgb2xyz(lab, lab); xyz2lab(lab, lab, ill);}//hsl->rgb->xyz->lab
}

#endif//_COLORSPACE_H_
