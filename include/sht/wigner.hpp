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

#ifndef _WIGNER_H_
#define _WIGNER_H_

#include <complex>

namespace emsphinx {

	namespace wigner {
		//@remarks
		//everyone seems to use different variables for degree/order of wigner functions
		//wigner (lowercase) d functions use the notation of: Fukushima, Toshio. (2016). Numerical computation of Wigner's d-function of arbitrary high degree and orders by extending exponent of floating point numbers.
		//https://doi.org/10.13140/RG.2.2.31922.20160
		//that means I use j, k, and m for degree and order such that wigner (lowercase) d functions are:
		//  d^j_{k, m}(beta) = sqrt( ( (j+k)! * (j-k)! ) / ( (j+m)! * (j-m)! ) ) * pow(cos(beta / 2), k+m) * pow(sin(beta / 2), k-m) * P^{k-m, k+m}_{j-k}(cos(beta))
		//which is the the equation 1 in the Fukushima reference
		//P^{a, b}_{n}(x) is the Jacobi polynomial of degree n and orders a/b evaluated at x e.g. the mathematica function JacobiP[n, a, b, x]
		//there additionally several conventions for the wigner (uppercase) D function
		//I'll use the convention such that for passive ZYZ euler angles {alpha, beta, gamma) the SHT of a function can be rotated by aRot^l_k = \sum_{m=-1}^l a^l_m D^l_{k,m}(ZYZ)
		//
		//the Wigner (lowercase) d function is computed recursively so it is expensive to calculate a single value
		//table functions are provided to generate d^j_{k,m} efficiently for all values of j/k/m (to a desired max j)
		//the Wigner (uppercase) D function required complex exponentials and a (lowercase) d function so it is even more expensive
		//functions to evaluate a single upper/lowercase wigner value are provided primarily for sanity checks and comparision against other conventions
		//
		//if D^j_{k,m}(\alpha, \beta, \gamma) is needed for all j, k, and m it is best to
		//  -compute d^j_{k,m}(\beta) for all values (with a table function)
		//  -loop over m computing exp(I*m*alpha)
		//  -loop over k computing exp(I*k*gamma)
		//  -compute D^j_{k,m}(\alpha, \beta, \gamma) == d^j_{k,m}(\beta) * exp(I*m*alpha) * exp(I*k*gamma)

		//@brief   : compute Wigner (lowercase) d function (also called reduced wigner d)
		//@param j : degree in d^j_{k,m}(beta)
		//@param k : first order in d^j_{k,m}(beta)
		//@param m : second order in d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to D^j_{k,m}({0, beta, 0})
		//@note    : equivalent to the mathematica function WignerD[{j, k, m}, beta]
		template <typename Real> Real d(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB);

		//@brief  : compute Wigner (lowercase) d function (also called reduced wigner d) at pi/2
		//@param j: degree
		//@param k: first order
		//@param m: second order
		//@return : d^j_{k,m}(\frac{\pi}{2})
		//@note   : NAN when j < max(|k|, |m|)
		//@note   : equivalent to the mathematica function WignerD[{j, k, m}, Pi/2]
		template <typename Real> Real d(const int64_t j, const int64_t k, const int64_t m);

		//@brief  : compute symmetry of wigner (lowercase) d function @ pi/2
		//@param j: degree
		//@param k: first order
		//@param m: second order
		//@return : +/-1 such that dSign(j, k, m) * d(j, |k|, |m|) == d(j, k, m)
		int dSign(const int64_t j, const int64_t k, const int64_t m);

		//@brief   : compute Wigner (uppercase) D function
		//@param j : degree in d^j_{k,m}(beta)
		//@param k : first order in d^j_{k,m}(beta)
		//@param m : second order in d^j_{k,m}(beta)
		//@param eu: ZYZ euler angles
		//@return  : D^j_{k,m}(eu)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function WignerD[{j, k, m}, eu[2], eu[1], eu[0]]
		//@note    : for ZYZ euler angles this d^j_{k,m}(eu[1]) * exp(I*m*eu[0] + I*k*eu[2])
		template <typename Real> std::complex<Real> D(const int64_t j, const int64_t k, const int64_t m, Real const * const eu);

		//@brief      : compute a table of Wigner (lowercase) d functions at 0 <= beta <= pi
		//@param jMax : maximum order to compute table for
		//@param t    : cos(beta)
		//@param table: location to write d^j_{k,m}(beta) for all non negative j, k, m (must have space for jMax * jMax * jMax * 2 Reals)
		//@note       : the table has unused (an uninitialized space) for j < max(|k|, |m|)
		//@note       : d^j_{k,m}(   beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 0]
		//@note       : d^j_{k,m}(pi-beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 1]
		//@note       : negative k/m values can be found with the following symmetry relationships
		//              d^j_{-k,-m}( beta) = (-1)^(k-m) d^j_{k,m}(     beta)
		//              d^j_{ k,-m}( beta) = (-1)^(j+k) d^j_{k,m}(pi - beta)
		//              d^j_{-k, m}( beta) = (-1)^(j+m) d^j_{k,m}(pi - beta)
		template <typename Real> void dTable(const size_t jMax, const Real t, const bool nB, Real * const table);

		//@brief      : compute a table of Wigner (lowercase) d functions at 0 <= beta <= pi [with precomputed coefficient tables]
		//@param jMax : maximum order to compute table for
		//@param t    : cos(beta)
		//@param table: location to write d^j_{k,m}(beta) for all non negative j, k, m (must have space for jMax * jMax * jMax * 2 Reals)
		//@param pE   : precomputed table of e_km
		//@param pW   : precomputed table of w_jkm
		//@param pB   : precomputed table of b_jkm
		//@note       : the table has unused (an uninitialized space) for j < max(|k|, |m|)
		//@note       : d^j_{k,m}(   beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 0]
		//@note       : d^j_{k,m}(pi-beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 1]
		//@note       : negative k/m values can be found with the following symmetry relationships
		//              d^j_{-k,-m}( beta) = (-1)^(k-m) d^j_{k,m}(     beta)
		//              d^j_{ k,-m}( beta) = (-1)^(j+k) d^j_{k,m}(pi - beta)
		//              d^j_{-k, m}( beta) = (-1)^(j+m) d^j_{k,m}(pi - beta)
		template <typename Real> void dTablePre(const size_t jMax, const Real t, const bool nB, Real * const table, Real const * const pE, Real const * const pW, Real const * const pB);

		//@brief      : fill precomputed factor tables for dTablePre
		//@param jMax : maximum order to compute table for
		//@param pE   : location to write precomputed table of e_km
		//@param pW   : location to write precomputed table of w_jkm
		//@param pB   : location to write precomputed table of b_jkm
		template <typename Real> void dTablePreBuild(const size_t jMax, Real * const pE, Real * const pW, Real * const pB);

		//@brief      : compute a table of Wigner (lowercase) d functions at pi/2
		//@param jMax : max j
		//@param table: location to write d^j_{k,m}(\frac{\pi}{2}) for j = [0,jMax), k = [0,jMax), m = [0, jMax) (jMax^3 values)
		//@param trans: true/false to transpose k/m indices
		//@note       : d^j_{k,m} located at k * jMax * jMax + m * jMax + j (trans = false)
		//@note       : d^j_{k,m} located at m * jMax * jMax + k * jMax + j (trans = true )
		template <typename Real> void dTable(const size_t jMax, Real * const table, const bool trans = false);

		//@brief    : rotate the spherical harmonic transformation of a real function
		//@param bw : bandwidth of harmonic coefficients (max l exclusive)
		//@param alm: spherical harmonic coefficients to rotate with \hat{a}^l_{m} at alm[m * bw + j]
		//@param blm: location to write rotated spherical harmonic coefficients with \hat{b}^l_{m} at blm[m * bw + j]
		//@param zyz: rotation to apply as zyz euler angles
		//@note     : b^l_m = \sum_{n=-1}^l a^l_n D^l_{m,n}(qu)
		template <typename Real>
		void rotateHarmonics(const size_t bw, std::complex<Real> const * const alm, std::complex<Real> * const blm, Real const * const zyz);

		//@brief   : compute first derivative of Wigner (lowercase) d function (also called reduced wigner d) with respect to beta
		//@param j : degree in (d/dBeta) d^j_{k,m}(beta)
		//@param k : first order in (d/dBeta) d^j_{k,m}(beta)
		//@param m : second order in (d/dBeta) d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : (d/dBeta) d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function D[WignerD[{j, k, m}, beta], beta]
		//@note    : negative k/m values have the following symmetry relationships
		//           dPrime^j_{-k,-m}( beta) = (-1)^(k+m  ) dPrime^j_{k,m}(     beta)
		//           dPrime^j_{ k,-m}( beta) = (-1)^(j+k+1) dPrime^j_{k,m}(pi - beta)
		//           dPrime^j_{-k, m}( beta) = (-1)^(j+m+1) dPrime^j_{k,m}(pi - beta)
		template <typename Real> Real dPrime(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB);

		//@brief   : compute second derivative of Wigner (lowercase) d function (also called reduced wigner d) with respect to beta
		//@param j : degree in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param k : first order in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param m : second order in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : (d/dBeta)^2 d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function D[WignerD[{j, k, m}, beta], {beta, 2}]
		//@note    : negative k/m values have the following symmetry relationships
		//              dPrime2^j_{-k,-m}( beta) = (-1)^(k+m) dPrime2^j_{k,m}(     beta)
		//              dPrime2^j_{ k,-m}( beta) = (-1)^(j+k) dPrime2^j_{k,m}(pi - beta)
		//              dPrime2^j_{-k, m}( beta) = (-1)^(j+m) dPrime2^j_{k,m}(pi - beta)
		template <typename Real> Real dPrime2(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB);
	}
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <vector>

namespace emsphinx {

	namespace wigner {
		//simplified helper functions to compute Wigner (lowercase) d functions d^j_{k,m}(beta) for 0 <= beta <= pi and integer j,k,m
		//all functions use the notation of: Fukushima, Toshio. (2016). Numerical computation of Wigner's d-function of arbitrary high degree and orders by extending exponent of floating point numbers.
		//https://doi.org/10.13140/RG.2.2.31922.20160
		//blocks functions share the same basic structure but the factor of 2 has been removed since only whole (not half) integers are needed

		//@brief     : compute intermediate recursion coefficient u_{j,k,m} (equation 13)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@param t/tc: cos(beta) / 1 - cos(beta)
		//@return    : recursion coefficient
		//@note      : _0, _1, _2 for 0 <= beta < pi / 2, beta == pi / 2, and pi / 2 < beta <= pi / 2 respectively
		template <typename Real> inline Real    u_jkm_0(const int64_t j, const int64_t k, const int64_t m, const Real tc) {return -tc * ((j - 1) * j) - (k * m - (j - 1) * j);}
		                         inline int64_t u_jkm_1(const int64_t j, const int64_t k, const int64_t m               ) {return                     -  k * m               ;}
		template <typename Real> inline Real    u_jkm_2(const int64_t j, const int64_t k, const int64_t m, const Real t ) {return  t  * ((j - 1) * j) -  k * m               ;}

		//@brief     : compute intermediate recursion coefficient v_{j,k,m} (equation 14)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@return    : recursion coefficient
		template <typename Real> inline Real v_jkm  (const int64_t j, const int64_t k, const int64_t m               ) {return             std::sqrt( Real( (j+k-1) * (j-k-1) * (j+m-1) * (j-m-1) ) ) * (j  )  ;}

		//@brief     : compute intermediate recursion coefficient w_{j,k,m} (equation 15)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@return    : recursion coefficient
		//@note      : this function is most susceptible to integer overflow (particularly at k = m = 0) with a max k of only 215 for 32 bit ints (2^31-1)^(1/4), 64 bit integers buys up to ~55k which should be plenty
		template <typename Real> inline Real w_jkm  (const int64_t j, const int64_t k, const int64_t m               ) {return Real(1) / ( std::sqrt( Real( (j+k  ) * (j-k  ) * (j+m  ) * (j-m  ) ) ) * (j-1) );}

		//@brief     : compute intermediate recursion coefficient a_{j,k,m} (equation 11)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@param t/tc: cos(beta) / 1 - cos(beta)
		//@return    : recursion coefficient
		//@note      : _0, _1, _2 for 0 <= beta < pi / 2, beta == pi / 2, and pi / 2 < beta <= pi / 2 respectively
		template <typename Real> inline Real a_jkm_0(const int64_t j, const int64_t k, const int64_t m, const Real tc) {return w_jkm<Real>(j, k, m) * ( u_jkm_0<Real>(j, k, m, tc) * (2*j-1) );}
		template <typename Real> inline Real a_jkm_1(const int64_t j, const int64_t k, const int64_t m               ) {return w_jkm<Real>(j, k, m) * ( u_jkm_1      (j, k, m    ) * (2*j-1) );}
		template <typename Real> inline Real a_jkm_2(const int64_t j, const int64_t k, const int64_t m, const Real t ) {return w_jkm<Real>(j, k, m) * ( u_jkm_2<Real>(j, k, m, t ) * (2*j-1) );}

		//@brief     : compute intermediate recursion coefficient a_{j,k,m} (equation 11) [with precomputed w_jkm]
		//@param w   : w(j,k,m)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@param t/tc: cos(beta) / 1 - cos(beta)
		//@return    : recursion coefficient
		//@note      : _0, _1, _2 for 0 <= beta < pi / 2, beta == pi / 2, and pi / 2 < beta <= pi / 2 respectively
		template <typename Real> inline Real a_jkm_0_pre(const Real w, const int64_t j, const int64_t k, const int64_t m, const Real tc) {return w * ( u_jkm_0<Real>(j, k, m, tc) * (2*j-1) );}
		template <typename Real> inline Real a_jkm_1_pre(const Real w, const int64_t j, const int64_t k, const int64_t m               ) {return w * ( u_jkm_1      (j, k, m    ) * (2*j-1) );}
		template <typename Real> inline Real a_jkm_2_pre(const Real w, const int64_t j, const int64_t k, const int64_t m, const Real t ) {return w * ( u_jkm_2<Real>(j, k, m, t ) * (2*j-1) );}

		//@brief     : compute intermediate recursion coefficient a_{j,k,m} (equation 12)
		//@param j   : degree in d^j_{k,m}
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@return    : recursion coefficient
		template <typename Real> inline Real b_jkm  (const int64_t j, const int64_t k, const int64_t m               ) {return w_jkm<Real>(j, k, m) * v_jkm<Real>(j, k, m);}

		//@brief     : compute recursion seed coefficient u_{k,m} (equation 23)
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@return    : recursion seed
		template <typename Real> inline Real u_km_0 (             const int64_t k, const int64_t m, const Real tc) {return     (-tc * (k + 1) - (m - 1 - k));}
		template <typename Real> inline Real u_km_1 (             const int64_t k, const int64_t m               ) {return Real(              -  m         );}
		template <typename Real> inline Real u_km_2 (             const int64_t k, const int64_t m, const Real t ) {return     ( t  * (k + 1) -  m         );}

		//@brief     : compute recursion seed coefficient a_{k,m} (equation 22)
		//@param k   : first order in d^j_{k,m}
		//@param m   : second order in d^j_{k,m}
		//@param t/tc: cos(beta) / 1 - cos(beta)
		//@return    : recursion seed
		//@note      : _0, _1, _2 for 0 <= beta < pi / 2, beta == pi / 2, and pi / 2 < beta <= pi / 2 respectively
		template <typename Real> inline Real a_km_0 (             const int64_t k, const int64_t m, const Real tc) {return std::sqrt( Real( 2*k+1 ) / ( (k+m+1) * (k-m+1) ) ) * u_km_0<Real>(k, m, tc);}
		template <typename Real> inline Real a_km_1 (             const int64_t k, const int64_t m               ) {return std::sqrt( Real( 2*k+1 ) / ( (k+m+1) * (k-m+1) ) ) * u_km_1<Real>(k, m    );}
		template <typename Real> inline Real a_km_2 (             const int64_t k, const int64_t m, const Real t ) {return std::sqrt( Real( 2*k+1 ) / ( (k+m+1) * (k-m+1) ) ) * u_km_2<Real>(k, m, t );}

		//@brief: compute recursion seed coefficient e_{k,m} = \sqrt{\frac{(2k)!}{(k+m)!(k-m)!}} recursively (m <= k) (equation 21)
		//@param k: k in d^k_{k,m}
		//@param m: m in d^k_{k,m}
		//@return: e_km where d^k_{k,m} = 2^-k * e_km
		template <typename Real> inline Real e_km(const int64_t k, const int64_t m) {
			Real e_lm = 1;//e_mm;
			for(int64_t l = m+1; l <= k; l++) e_lm *= std::sqrt( Real( l*(2*l-1) ) / ( 2 * (l+m) * (l-m) ) ) * 2;
			return e_lm;//e_km
		}

		//@brief   : compute Wigner (lowercase) d function (also called reduced wigner d)
		//@param j : degree in d^j_{k,m}(beta)
		//@param k : first order in d^j_{k,m}(beta)
		//@param m : second order in d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to D^j_{k,m}({0, beta, 0})
		//@note    : equivalent to the mathematica function WignerD[{j, k, m}, beta]
		template <typename Real>
		Real d(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB) {
			//require 0 <= m <= k <= j and beta >= 0 (handle other cases with symmetry)
			if(nB) {                   //d^j_{ k, m}(-beta) =                d^j_{m,k}(     beta)
				return d(j, m, k, t, false);                //equation 5
			} else if(k < 0 && m < 0) {//d^j_{-k,-m}( beta) = (-1)^(   k- m) d^j_{k,m}(     beta)
				const int sign = (k-m) % 2 == 0 ? 1 : -1;
				return d<Real>(j, -k, -m,  t, false) * sign;//equation 6
			} else if(         m < 0) {//d^j_{ k,-m}( beta) = (-1)^(j+ k+2m) d^j_{k,m}(pi - beta)
				const int sign = (j+k) % 2 == 0 ? 1 : -1;
				return d<Real>(j,  k, -m, -t, false) * sign;//equation 7
			} else if(k < 0         ) {//d^j_{-k, m}( beta) = (-1)^(j+2k+3m) d^j_{k,m}(pi - beta)
				const int sign = (j+m) % 2 == 0 ? 1 : -1;
				return d<Real>(j, -k,  m, -t, false) * sign;//equation 8
			} else if(k     <  m    ) {//d^j_{ m, k}( beta) = (-1)^(   k- m) d^j_{k,m}(     beta)
				const int sign = (k-m) % 2 == 0 ? 1 : -1;
				return d<Real>(j,  m,  k,  t, false) * sign;//equation 9
			}

			if(j < k) return NAN;

			//determine if beta is < (0), > (2), or = (1) to pi/2
			const size_t type = t > 0 ? 0 : (t < 0 ? 2 : 1);
			const Real tc = Real(1) - t;

			//compute powers of cos/sin of beta / 2
			const Real c2 = std::sqrt( (Real(1) + t) / 2 );//cos(acos(t)) == cos(beta / 2), always positive since at this point 0 <= beta <= pi
			const Real s2 = std::sqrt( (Real(1) - t) / 2 );//sin(acos(t)) == sin(beta / 2), always positive since at this point 0 <= beta <= pi
			const Real cn = std::pow(c2, Real(k+m));//equation 20 for n = k+m
			const Real sn = std::pow(s2, Real(k-m));//equation 20 for n = k-m

			//compute first term for three term recursion 
			const Real d_kkm  = cn * sn * e_km<Real>(k, m);//equation 18, d^k_{k, m}(beta)
			if(j == k  ) return d_kkm;//if j == k we're done

			//compute second term for three term recursion 
			Real a_km;
			switch(type) {
				case 0: a_km = a_km_0<Real>(k, m, tc); break;//beta <  pi/2
				case 1: a_km = a_km_1<Real>(k, m    ); break;//beta == pi/2
				case 2: a_km = a_km_2<Real>(k, m, t ); break;//beta >  pi/2
			}
			const Real d_k1km = d_kkm * a_km;//equation 19, d^{k+1}_{k, m}(beta)
			if(j == k+1) return d_k1km;//if j == k + 1 we're done

			//recursively compute by degree to j
			Real d_ikm;
			Real d_i2km = d_kkm ;
			Real d_i1km = d_k1km;
			switch(type) {
				case 0://beta <  pi/2
					for(int64_t i = k + 2; i <= j; i++) {
						d_ikm = a_jkm_0<Real>(i, k, m, tc) * d_i1km - b_jkm<Real>(i, k, m) * d_i2km;//equation 10, d^i_{k, m}(beta)
						d_i2km = d_i1km;
						d_i1km = d_ikm ;
					}
					break;
				case 1://beta == pi/2
					for(int64_t i = k + 2; i <= j; i++) {
						d_ikm = a_jkm_1<Real>(i, k, m    ) * d_i1km - b_jkm<Real>(i, k, m) * d_i2km;//equation 10, d^i_{k, m}(beta)
						d_i2km = d_i1km;
						d_i1km = d_ikm ;
					}
					break;
				case 2://beta >  pi/2
					for(int64_t i = k + 2; i <= j; i++) {
						d_ikm = a_jkm_2<Real>(i, k, m, t ) * d_i1km - b_jkm<Real>(i, k, m) * d_i2km;//equation 10, d^i_{k, m}(beta)
						d_i2km = d_i1km;
						d_i1km = d_ikm ;
					}
					break;
			}
			return d_ikm;
		}

		//@brief  : compute Wigner (lowercase) d function (also called reduced wigner d) at pi/2
		//@param j: degree
		//@param k: first order
		//@param m: second order
		//@return : d^j_{k,m}(\frac{\pi}{2})
		//@note   : NAN when j < max(|k|, |m|)
		//@note   : equivalent to the mathematica function WignerD[{j, k, m}, Pi/2]
		template <typename Real> Real d(const int64_t j, const int64_t k, const int64_t m) {
			//require 0 <= m <= k <= j (handle with symmetry where possible)
			if(k < 0 && m < 0) {//d^j_{-k,-m} = (-1)^(   k- m) d^j_{k,m}
				return (    k-  m) % 2 == 0 ? d<Real>(j, -k, -m) : -d<Real>(j, -k, -m);
			} else if(m < 0) {//d^j_{ k,-m} = (-1)^(j+ k+2m) d^j_{k,m}
				return (j+  k+2*m) % 2 == 0 ? d<Real>(j,  k, -m) : -d<Real>(j,  k, -m);
			} else if(k < 0) {//d^j_{-k, m} = (-1)^(j+2k+3m) d^j_{k,m}
				return (j+2*k+3*m) % 2 == 0 ? d<Real>(j, -k,  m) : -d<Real>(j, -k,  m);
			} else if(k < m) {//d^j_{ m, k} = (-1)^(   k- m) d^j_{k,m}
				return (    k-  m) % 2 == 0 ? d<Real>(j,  m,  k) : -d<Real>(j,  m,  k);
			}
			if(j < k) return NAN;

			//compute first two terms for three term recursion 
			const Real d_kkm  = std::pow(Real(2), Real(-k)) * e_km<Real>(k, m);//equation 18, d^k_{k, m}(pi/2)
			if(j == k  ) return d_kkm;
			const Real d_k1km = d_kkm * a_km_1<Real>(k, m);//equation 19, d^{k+1}_{k, m}(pi/2)
			if(j == k+1) return d_k1km;

			//recursively compute
			Real d_ikm;
			Real d_i2km = d_kkm ;
			Real d_i1km = d_k1km;
			for(int64_t i = k + 2; i <= j; i++) {
				d_ikm = a_jkm_1<Real>(i, k, m) * d_i1km - b_jkm<Real>(i, k, m) * d_i2km;//equation 10, d^i_{k, m}(pi/2)
				d_i2km = d_i1km;
				d_i1km = d_ikm ;
			}
			return d_ikm;
		}

		//@brief  : compute symmetry of wigner (lowercase) d function @ pi/2
		//@param j: degree
		//@param k: first order
		//@param m: second order
		//@return : +/-1 such that dSign(j, k, m) * d(j, |k|, |m|) == d(j, k, m)
		int dSign(const int64_t j, const int64_t k, const int64_t m) {
			if(k < 0 && m < 0) {//d^j_{-k,-m} = (-1)^(   k- m) d^j_{k,m}
				return (k-m) % 2 == 0 ? 1 : -1;
			} else if(m < 0) {//d^j_{ k,-m} = (-1)^(j+ k+2m) d^j_{k,m}
				return (j+k) % 2 == 0 ? 1 : -1;
			} else if(k < 0) {//d^j_{-k, m} = (-1)^(j+2k+3m) d^j_{k,m}
				return (j+m) % 2 == 0 ? 1 : -1;
			}
			return 1;
		}

		//@brief   : compute Wigner (uppercase) D function
		//@param j : degree in d^j_{k,m}(beta)
		//@param k : first order in d^j_{k,m}(beta)
		//@param m : second order in d^j_{k,m}(beta)
		//@param eu: ZYZ euler angles
		//@return  : D^j_{k,m}(eu)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function WignerD[{j, k, m}, eu[2], eu[1], eu[0]]
		//@note    : for ZYZ euler angles this d^j_{k,m}(eu[1]) * exp(I*m*eu[0] + I*k*eu[2])
		template <typename Real> std::complex<Real> D(const int64_t j, const int64_t k, const int64_t m, Real const * const eu) {
			const Real sum = eu[0] * m + eu[2] * k;
			return std::complex<Real>(std::cos(sum), std::sin(sum)) * d(j, k, m, std::cos(eu[1]), std::signbit(eu[1]));
		}

		//@brief      : compute a table of Wigner (lowercase) d functions at 0 <= beta <= pi
		//@param jMax : maximum order to compute table for
		//@param t    : cos(beta)
		//@param table: location to write d^j_{k,m}(beta) for all non negative j, k, m (must have space for jMax * jMax * jMax * 2 Reals)
		//@note       : the table has unused (an uninitialized space) for j < max(|k|, |m|)
		//@note       : d^j_{k,m}(   beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 0]
		//@note       : d^j_{k,m}(pi-beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 1]
		//@note       : negative k/m values can be found with the following symmetry relationships
		//              d^j_{-k,-m}( beta) = (-1)^(k-m) d^j_{k,m}(     beta)
		//              d^j_{ k,-m}( beta) = (-1)^(j+k) d^j_{k,m}(pi - beta)
		//              d^j_{-k, m}( beta) = (-1)^(j+m) d^j_{k,m}(pi - beta)
		template <typename Real>
		void dTable(const size_t jMax, const Real t, const bool nB, Real * const table) {
			/*
			//naive implementation (has a bunch of redundant calculations)
			for(int k = 0; k < jMax; k++) {
				for(int m = 0; m < jMax; m++) {
					for(int j = 0; j < jMax; j++) {
						table[k * jMax * jMax * 2 + m * jMax * 2 + j * 2 + 0] = d<Real>(j, k, m,  t, false);//d^j_{k, m}(beta)
						table[k * jMax * jMax * 2 + m * jMax * 2 + j * 2 + 1] = d<Real>(j, k, m, -t, false);//d^j_{k, m}(pi - beta)
					}
				}
			}
			*/

			//determine which branch of function is needed and compute cos/sin of half angles
			const bool isType0 = !std::signbit(t);//is this type 0 or type 2? (type 1 will be grouped with type 0)
			const Real tc  = Real(1) - t;
			const Real tcN = Real(1) + t;//tc for -t
			const Real c2  = std::sqrt(tcN / 2);//cos(acos(t)) == cos(beta / 2), always positive since at this point 0 <= beta <= pi
			const Real s2  = std::sqrt(tc  / 2);//sin(acos(t)) == sin(beta / 2), always positive since at this point 0 <= beta <= pi
			
			//get function pointers to avoid repeated branch in loop
			Real (* const a_kmFunc  )(               const int64_t, const int64_t, const Real) = isType0 ? a_km_0 <Real> : a_km_2 <Real>;
			Real (* const a_kmFuncN )(               const int64_t, const int64_t, const Real) = isType0 ? a_km_2 <Real> : a_km_0 <Real>;
			Real (* const a_jkmFunc )(const int64_t, const int64_t, const int64_t, const Real) = isType0 ? a_jkm_0<Real> : a_jkm_2<Real>;
			Real (* const a_jkmFuncN)(const int64_t, const int64_t, const int64_t, const Real) = isType0 ? a_jkm_2<Real> : a_jkm_0<Real>;
			const Real t0 = isType0 ? tc : t  ;
			const Real tN = isType0 ? -t : tcN;

			//precompute integer powers of c2 and s2
			std::vector<Real> work(jMax * 4);
			Real* const pc2 = work.data();
			Real* const ps2 = work.data() + jMax * 2;
			for(size_t i = 0; i < jMax * 2; i++) {
				pc2[i] = std::pow(c2, Real(i));
				ps2[i] = std::pow(s2, Real(i));
			} 

			//fill table
			Real* pK = table;//table + k * jMax * jMAx * 2
			for(size_t k = 0; k < jMax; k++) {
				Real* pKM = pK                  ;//table + k * jMax * jMax * 2 + m * jMax * 2
				Real* pMK = table + k * jMax * 2;//table + m * jMax * jMax * 2 + k * jMax * 2
				for(size_t m = 0; m <= k; m++) {
					//determine sign change for swapping k/m
					int sign = 1;
					int signN = (k-m) % 2 == 0 ? 1 : -1;//symmetry from eq 9
					if(nB) std::swap(sign, signN);

					//compute powers of cos/sin of beta / 2
					const Real cn  = pc2[k+m];//std::pow(c2, k+m);//equation 20 for n = k+m for  t
					const Real sn  = ps2[k-m];//std::pow(s2, k-m);//equation 20 for n = k-m for  t
					const Real cnN = ps2[k+m];//std::pow(s2, k+m);//equation 20 for n = k+m for -t
					const Real snN = pc2[k-m];//std::pow(c2, k-m);//equation 20 for n = k-m for -t

					//compute first term for three term recursion 
					const Real ekm = e_km<Real>(k, m);
					const Real d_kkm  = cn  * sn  * ekm;//equation 18, d^k_{k, m}(beta) for  t
					const Real d_kkmN = cnN * snN * ekm;//equation 18, d^k_{k, m}(beta) for -t
					pKM[k * 2 + 0] = d_kkm  * sign ;
					pMK[k * 2 + 0] = d_kkm  * signN;//symmetry from eq 9
					pKM[k * 2 + 1] = d_kkmN * sign ;
					pMK[k * 2 + 1] = d_kkmN * signN;//symmetry from eq 9

					if(k+1 < jMax) {//if j == k we're done
						//compute second term for three term recursion 
						const Real a_km  = a_kmFunc (k, m, t0);
						const Real a_kmN = a_kmFuncN(k, m, tN);
						const Real d_k1km  = d_kkm  * a_km ;//equation 19, d^{k+1}_{k, m}(beta) for  t
						const Real d_k1kmN = d_kkmN * a_kmN;//equation 19, d^{k+1}_{k, m}(beta) for -t
						pKM[(k+1) * 2 + 0] = d_k1km  * sign ;
						pMK[(k+1) * 2 + 0] = d_k1km  * signN;//symmetry from eq 9
						pKM[(k+1) * 2 + 1] = d_k1kmN * sign ;
						pMK[(k+1) * 2 + 1] = d_k1kmN * signN;//symmetry from eq 9

						if(k+2 < jMax) {//if j == k + 1 we're done
							//recursively compute by degree to j
							Real d_ikm;
							Real d_i2km = d_kkm ;
							Real d_i1km = d_k1km;
							Real d_ikmN;
							Real d_i2kmN = d_kkmN ;
							Real d_i1kmN = d_k1kmN;

							for(size_t i = k + 2; i < jMax; i++) {
								d_ikm   = a_jkmFunc (i, k, m, t0) * d_i1km  - b_jkm<Real>(i, k, m) * d_i2km ;//equation 10, d^i_{k, m}(beta)
								d_i2km  = d_i1km;
								d_i1km  = d_ikm ;
								d_ikmN  = a_jkmFuncN(i, k, m, tN) * d_i1kmN - b_jkm<Real>(i, k, m) * d_i2kmN;//equation 10, d^i_{k, m}(beta)
								d_i2kmN = d_i1kmN;
								d_i1kmN = d_ikmN ;
								pKM[i * 2 + 0] = d_i1km  * sign ;
								pMK[i * 2 + 0] = d_i1km  * signN;
								pKM[i * 2 + 1] = d_i1kmN * sign ;
								pMK[i * 2 + 1] = d_i1kmN * signN;
							}
						}
					}

					//increment table pointers
					pKM += jMax        * 2;
					pMK += jMax * jMax * 2;
				}

				//increment table pointers
				pK += jMax * jMax * 2;
			}
		}

		//@brief      : compute a table of Wigner (lowercase) d functions at 0 <= beta <= pi [with precomputed coefficient tables]
		//@param jMax : maximum order to compute table for
		//@param t    : cos(beta)
		//@param table: location to write d^j_{k,m}(beta) for all non negative j, k, m (must have space for jMax * jMax * jMax * 2 Reals)
		//@param pE   : precomputed table of e_km
		//@param pW   : precomputed table of w_jkm
		//@param pB   : precomputed table of b_jkm
		//@note       : the table has unused (an uninitialized space) for j < max(|k|, |m|)
		//@note       : d^j_{k,m}(   beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 0]
		//@note       : d^j_{k,m}(pi-beta) at table[(k * jMax * jMax + m * jMax + j)*2 + 1]
		//@note       : negative k/m values can be found with the following symmetry relationships
		//              d^j_{-k,-m}( beta) = (-1)^(k-m) d^j_{k,m}(     beta)
		//              d^j_{ k,-m}( beta) = (-1)^(j+k) d^j_{k,m}(pi - beta)
		//              d^j_{-k, m}( beta) = (-1)^(j+m) d^j_{k,m}(pi - beta)
		template <typename Real>
		void dTablePre(const size_t jMax, const Real t, const bool nB, Real * const table, Real const * const pE, Real const * const pW, Real const * const pB) {
			//determine which branch of function is needed and compute cos/sin of half angles
			const bool isType0 = !std::signbit(t);//is this type 0 or type 2? (type 1 will be grouped with type 0)
			const Real tc  = Real(1) - t;
			const Real tcN = Real(1) + t;//tc for -t
			const Real c2  = std::sqrt(tcN / 2);//cos(acos(t)) == cos(beta / 2), always positive since at this point 0 <= beta <= pi
			const Real s2  = std::sqrt(tc  / 2);//sin(acos(t)) == sin(beta / 2), always positive since at this point 0 <= beta <= pi
			
			//get function pointers to avoid repeated branch in loop
			Real (* const a_kmFunc  )(            const int64_t, const int64_t,                const Real) = isType0 ? a_km_0     <Real> : a_km_2     <Real>;
			Real (* const a_kmFuncN )(            const int64_t, const int64_t,                const Real) = isType0 ? a_km_2     <Real> : a_km_0     <Real>;
			Real (* const a_jkmFunc )(const Real, const int64_t, const int64_t, const int64_t, const Real) = isType0 ? a_jkm_0_pre<Real> : a_jkm_2_pre<Real>;
			Real (* const a_jkmFuncN)(const Real, const int64_t, const int64_t, const int64_t, const Real) = isType0 ? a_jkm_2_pre<Real> : a_jkm_0_pre<Real>;
			const Real t0 = isType0 ? tc : t  ;
			const Real tN = isType0 ? -t : tcN;

			//precompute integer powers of c2 and s2
			std::vector<Real> work(jMax * 4);
			Real* const pc2 = work.data();
			Real* const ps2 = work.data() + jMax * 2;
			for(size_t i = 0; i < jMax * 2; i++) {
				pc2[i] = std::pow(c2, Real(i));
				ps2[i] = std::pow(s2, Real(i));
			} 

			//fill table
			Real* pK = table;//table + k * jMax * jMAx * 2
			for(size_t k = 0; k < jMax; k++) {
				Real* pKM = pK                  ;//table + k * jMax * jMax * 2 + m * jMax * 2
				Real* pMK = table + k * jMax * 2;//table + m * jMax * jMax * 2 + k * jMax * 2
				for(size_t m = 0; m <= k; m++) {
					//determine sign change for swapping k/m
					int sign = 1;
					int signN = (k-m) % 2 == 0 ? 1 : -1;//symmetry from eq 9
					if(nB) std::swap(sign, signN);

					//compute powers of cos/sin of beta / 2
					const Real cn  = pc2[k+m];//std::pow(c2, k+m);//equation 20 for n = k+m for  t
					const Real sn  = ps2[k-m];//std::pow(s2, k-m);//equation 20 for n = k-m for  t
					const Real cnN = ps2[k+m];//std::pow(s2, k+m);//equation 20 for n = k+m for -t
					const Real snN = pc2[k-m];//std::pow(c2, k-m);//equation 20 for n = k-m for -t

					//compute first term for three term recursion 
					const Real& ekm = pE[k * jMax + m];//e_km<Real>(k, m)
					const Real d_kkm  = cn  * sn  * ekm;//equation 18, d^k_{k, m}(beta) for  t
					const Real d_kkmN = cnN * snN * ekm;//equation 18, d^k_{k, m}(beta) for -t
					pKM[k * 2 + 0] = d_kkm  * sign ;
					pMK[k * 2 + 0] = d_kkm  * signN;//symmetry from eq 9
					pKM[k * 2 + 1] = d_kkmN * sign ;
					pMK[k * 2 + 1] = d_kkmN * signN;//symmetry from eq 9

					if(k+1 < jMax) {//if j == k we're done
						//compute second term for three term recursion 
						const Real a_km  = a_kmFunc (k, m, t0);
						const Real a_kmN = a_kmFuncN(k, m, tN);
						const Real d_k1km  = d_kkm  * a_km ;//equation 19, d^{k+1}_{k, m}(beta) for  t
						const Real d_k1kmN = d_kkmN * a_kmN;//equation 19, d^{k+1}_{k, m}(beta) for -t
						pKM[(k+1) * 2 + 0] = d_k1km  * sign ;
						pMK[(k+1) * 2 + 0] = d_k1km  * signN;//symmetry from eq 9
						pKM[(k+1) * 2 + 1] = d_k1kmN * sign ;
						pMK[(k+1) * 2 + 1] = d_k1kmN * signN;//symmetry from eq 9

						if(k+2 < jMax) {//if j == k + 1 we're done
							//recursively compute by degree to j
							Real d_ikm;
							Real d_i2km = d_kkm ;
							Real d_i1km = d_k1km;
							Real d_ikmN;
							Real d_i2kmN = d_kkmN ;
							Real d_i1kmN = d_k1kmN;

							for(size_t i = k + 2; i < jMax; i++) {
								const size_t idx = k * jMax * jMax + m * jMax + i;
								d_ikm   = a_jkmFunc (pW[idx], i, k, m, t0) * d_i1km  - pB[idx] * d_i2km ;//equation 10, d^i_{k, m}(beta), precomputed (a/b)_jkm<Real>(i, k, m)
								d_i2km  = d_i1km;
								d_i1km  = d_ikm ;
								d_ikmN  = a_jkmFuncN(pW[idx], i, k, m, tN) * d_i1kmN - pB[idx] * d_i2kmN;//equation 10, d^i_{k, m}(beta), precomputed (a/b)_jkm<Real>(i, k, m)
								d_i2kmN = d_i1kmN;
								d_i1kmN = d_ikmN ;
								pKM[i * 2 + 0] = d_i1km  * sign ;
								pMK[i * 2 + 0] = d_i1km  * signN;
								pKM[i * 2 + 1] = d_i1kmN * sign ;
								pMK[i * 2 + 1] = d_i1kmN * signN;
							}
						}
					}

					//increment table pointers
					pKM += jMax        * 2;
					pMK += jMax * jMax * 2;
				}

				//increment table pointers
				pK += jMax * jMax * 2;
			}
		}

		//@brief      : fill precomputed factor tables for dTablePre
		//@param jMax : maximum order to compute table for
		//@param pE   : location to write precomputed table of e_km
		//@param pW   : location to write precomputed table of w_jkm
		//@param pB   : location to write precomputed table of b_jkm
		template <typename Real>
		void dTablePreBuild(const size_t jMax, Real * const pE, Real * const pW, Real * const pB) {
			//fill table
			for(size_t k = 0; k < jMax; k++) {
				for(size_t m = 0; m <= k; m++) {
					pE[k * jMax + m] = e_km<Real>(k, m);
					for(size_t i = k + 2; i < jMax; i++) {
						const size_t idx = k * jMax * jMax + m * jMax + i;
						pW[idx] = w_jkm<Real>(i, k, m);
						pB[idx] = b_jkm<Real>(i, k, m);
					}
				}
			}
		}

		//@brief      : compute a table of Wigner (lowercase) d functions at pi/2
		//@param jMax : max j
		//@param table: location to write d^j_{k,m}(\frac{\pi}{2}) for j = [0,jMax), k = [0,jMax), m = [0, jMax) (jMax^3 values)
		//@param trans: true/false to transpose k/m indices
		//@note       : d^j_{k,m} located at k * jMax * jMax + m * jMax + j (trans = false)
		//@note       : d^j_{k,m} located at m * jMax * jMax + k * jMax + j (trans = true )
		template <typename Real> void dTable(const size_t jMax, Real * const table, const bool trans) {
			//naive implementation (has a bunch of redundant calculations)
			/*
			for(int k = 0; k < jMax; k++) {
				for(int m = 0; m < jMax; m++) {
					for(int j = 0; j < jMax; j++) {
						table[k * jMax * jMax + m * jMax + j] = d<Real>(j, k, m, 0, false);//k and am swapped for trans = true
					}
				}
			}
			*/

			//recursively compute wigner d function values for all d^j_{k,m} where k >= m >= 0
			//use symmetry to fill in table for m < 0 and m > k: /d^j_{ m, k} = (-1)^(  k- m) d^j_{k,m}
			for(int k = 0; k < jMax; k++) {
				for(int m = 0; m <= k; m++) {
					const bool km = (k-m) % 2 == 0;

					//compute d^k_{k,m} with closed form solution
					const Real d_kkm = std::pow(Real(2), -k) * e_km<Real>(k, m);//d^{k  }_{k,m}
					const Real d_kmk = km ? d_kkm : -d_kkm;//d^{k  }_{m,k}
					if(trans) {
						table[m * jMax * jMax +  k * jMax + k    ] = d_kkm;//save d^{k  }_{k,m}
						table[k * jMax * jMax +  m * jMax + k    ] = d_kmk;//save d^{k  }_{m,k}
					} else {
						table[k * jMax * jMax +  m * jMax + k    ] = d_kkm;//save d^{k  }_{k,m}
						table[m * jMax * jMax +  k * jMax + k    ] = d_kmk;//save d^{k  }_{m,k}
					}

					if(k + 1 == jMax) continue;//we don't need any higher order terms

					//compute d^{k+1}_{k,m} with recursion
					const Real d_k1km = d_kkm * a_km_1<Real>(k, m);//d^{k+1}_{k,m}
					const Real d_k1mk = km ? d_k1km : -d_k1km;//d^{k+1}_{m,k}
					if(trans) {
						table[m * jMax * jMax +  k * jMax + k + 1] = d_k1km;//save d^{k+1}_{k,m}
						table[k * jMax * jMax +  m * jMax + k + 1] = d_k1mk;//save d^{k+1}_{m,k}
					} else {
						table[k * jMax * jMax +  m * jMax + k + 1] = d_k1km;//save d^{k+1}_{k,m}
						table[m * jMax * jMax +  k * jMax + k + 1] = d_k1mk;//save d^{k+1}_{m,k}
					}
					if(k + 2 == jMax) continue;//we don't need any higher order terms

					//compute higher order terms with 3 term recursion
					Real djkm;           //d^{j  }_{k,m}
					Real d_j1km = d_k1km;//d^{j-1}_{k,m}
					Real d_j2km = d_kkm ;//d^{j-2}_{k,m}
					for(int j = k + 2; j < jMax; j++) {
						djkm = a_jkm_1<Real>(j, k, m) * d_j1km - b_jkm<Real>(j, k, m) * d_j2km;//d^{j}_{k,m}
						const Real djak = km             ? djkm : -djkm;//d^{j}_{m,k}
						if(trans) {
							table[m * jMax * jMax +  k * jMax + j] = djkm;//save d^{j}_{k,m}
							table[k * jMax * jMax +  m * jMax + j] = djak;//save d^{j}_{m,k}
						} else {
							table[k * jMax * jMax +  m * jMax + j] = djkm;//save d^{j}_{k,m}
							table[m * jMax * jMax +  k * jMax + j] = djak;//save d^{j}_{m,k}
						}
						d_j2km = d_j1km;
						d_j1km = djkm ;
					}
				}
			}
		}

		//@brief    : rotate the spherical harmonic transformation of a real function
		//@param bw : bandwidth of harmonic coefficients (max l exclusive)
		//@param alm: spherical harmonic coefficients to rotate with \hat{a}^l_{m} at alm[m * bw + j]
		//@param blm: location to write rotated spherical harmonic coefficients with \hat{b}^l_{m} at blm[m * bw + j]
		//@param zyz: rotation to apply as zyz euler angles
		//@note     : b^l_m = \sum_{n=-1}^l a^l_n D^l_{m,n}(qu)
		template <typename Real>
		void rotateHarmonics(const size_t bw, std::complex<Real> const * const alm, std::complex<Real> * const blm, Real const * const zyz) {
			//clear output array
			std::vector< std::complex<Real> > almRot(bw * bw, std::complex<Real>(0));//working array in case alm == blm

			//construct wigner (lowercase) d lookup table for beta
			std::vector<Real> dBeta(bw * bw * bw * 2);
			wigner::dTable(bw, std::cos(zyz[1]), std::signbit(zyz[1]), dBeta.data());//compute wigner (lowercase) d(beta) once

			//now that we can easily compute D^l_{m,n}(qu) do the actual summation
			for(size_t m = 0; m < bw; m++) {
				const std::complex<Real> expAlpha(std::cos(zyz[2] * m), std::sin(zyz[2] * m));//exp(I m gamma)
				for(size_t n = 0; n < bw; n++) {
					const std::complex<Real> expGamma (std::cos(zyz[0] * n), std::sin(zyz[0] * n));//exp(I n alpha)
					for(size_t j = std::max<size_t>(m, n); j < bw; j++) {
						const std::complex<Real> alGamma = alm[n * bw + j] * expGamma;//\hat{a}^l_{n}
						const Real rr = expAlpha.real() * alGamma.real();//ac
						const Real ri = expAlpha.real() * alGamma.imag();//ad
						const Real ir = expAlpha.imag() * alGamma.real();//bc
						const Real ii = expAlpha.imag() * alGamma.imag();//bd
						const std::complex<Real> vp(rr - ii, ir + ri);//expAlpha *      alGamma  = \hat{a}^l_{+n} * exp(I m gamma) * exp(I +n alpha)
						const std::complex<Real> vc(rr + ii, ir - ri);//expAlpha * conj(alGamma) = \hat{a}^l_{-n} * exp(I m gamma) * exp(I -n alpha) * (-1)^n using symmetry of real SHT
						const Real& dmn0 = dBeta[(m * bw * bw + n * bw + j) * 2 + 0];//d^j_{m,n}(     beta)
						const Real& dmn1 = dBeta[(m * bw * bw + n * bw + j) * 2 + 1];//d^j_{m,n}(pi - beta) since d^j_{m,-n}( beta) = (-1)^(j+m) d^j_{m,n}(pi - beta)
						almRot[m * bw + j] += vp * dmn0;//\hat{a}^l_{+n} * D^l_{m,+n}(zyz)
						if(n > 0) almRot[m * bw + j] += vc * (dmn1 * Real(0 == (j+m+n) % 2 ? 1 : -1));//\hat{a}^l_{-n} * D^l_{m,-n}(zyz)
					}
				}
			}
			std::copy(almRot.begin(), almRot.end(), blm);
		}

		//@brief   : compute first derivative of Wigner (lowercase) d function (also called reduced wigner d) with respect to beta
		//@param j : degree in (d/dBeta) d^j_{k,m}(beta)
		//@param k : first order in (d/dBeta) d^j_{k,m}(beta)
		//@param m : second order in (d/dBeta) d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : (d/dBeta) d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function D[WignerD[{j, k, m}, beta], beta]
		//@note    : negative k/m values have the following symmetry relationships
		//           dPrime^j_{-k,-m}( beta) = (-1)^(k+m  ) dPrime^j_{k,m}(     beta)
		//           dPrime^j_{ k,-m}( beta) = (-1)^(j+k+1) dPrime^j_{k,m}(pi - beta)
		//           dPrime^j_{-k, m}( beta) = (-1)^(j+m+1) dPrime^j_{k,m}(pi - beta)
		template <typename Real> Real dPrime(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB) {
			//compute prefactor (same for all j, k, and m for a given beta)
			const Real csc = Real(1) / std::sqrt(Real(1) - t * t) * (nB ? -1 : 1);//csc(beta), cot(beta) is csc * t

			//compute derivative
			const Real d0Term =              d(j, k  , m, t, nB) * (t * k - m) * csc;//d^j_{k,m}(beta) * (k*cot(beta) - m*csc(beta))
			const Real d1Term = j == k ? 0 : d(j, k+1, m, t, nB) * std::sqrt( Real( (j - k) * (j + k + 1) ) );//d^j_{k+1,m}(beta) * (j+k+1) * sqrt( (j-k) / (j+kl+1) )
			return d0Term - d1Term;
		}

		//@brief   : compute second derivative of Wigner (lowercase) d function (also called reduced wigner d) with respect to beta
		//@param j : degree in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param k : first order in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param m : second order in (d/dBeta)^2 d^j_{k,m}(beta)
		//@param t : cos(beta)
		//@param nB: true/false for negative/positive beta
		//@return  : (d/dBeta)^2 d^j_{k,m}(beta)
		//@note    : NAN when j < max(|k|, |m|)
		//@note    : equivalent to the mathematica function D[WignerD[{j, k, m}, beta], {beta, 2}]
		//@note    : negative k/m values have the following symmetry relationships
		//              dPrime2^j_{-k,-m}( beta) = (-1)^(k+m) dPrime2^j_{k,m}(     beta)
		//              dPrime2^j_{ k,-m}( beta) = (-1)^(j+k) dPrime2^j_{k,m}(pi - beta)
		//              dPrime2^j_{-k, m}( beta) = (-1)^(j+m) dPrime2^j_{k,m}(pi - beta)
		template <typename Real> Real dPrime2(const int64_t j, const int64_t k, const int64_t m, const Real t, const bool nB) {
			//compute prefactor (same for all j, k, and m for a given beta)
			const Real csc = Real(1) / std::sqrt(Real(1) - t * t) * (nB ? -1 : 1);//csc(beta), cot(beta) is csc * t

			//compute derivative prefactors
			const Real rjk  = std::sqrt( Real( (j - k    ) * (j + k + 1) ) );
			const Real d0Coef = ( t * t * k * k + t * m * (1 - 2 * k) + (m * m - k) ) * csc * csc;
			const Real d1Coef = rjk * (t * (1 + 2 * k) - 2 * m) * csc;
			const Real d2Coef = rjk * std::sqrt( Real( (j - k - 1) * (j + k + 2) ) );

			//compute derivative
			const Real d0Term =                d(j, k  , m, t, nB) * d0Coef;
			const Real d1Term = k   >= j ? 0 : d(j, k+1, m, t, nB) * d1Coef;
			const Real d2Term = k+1 >= j ? 0 : d(j, k+2, m, t, nB) * d2Coef;
			return d0Term - d1Term + d2Term;
		}
	}

}

#endif//_WIGNER_H_
