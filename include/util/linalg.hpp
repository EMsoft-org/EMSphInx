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

#ifndef _LINALG_H_
#define _LINALG_H_

#include <cstdlib>//size_t

// this header implements a minimal set of matrix decompositions preferring compactness and readability over speed
// it should only be used if you only need to do a few decompositions and/or the matrix size is small
// if you need more obscure decompositions and/or speed is crucial use a proper library (e.g. Eigen or LAPACK)

//@brief: functions to solve a linear system
//@note : these are convenience methods combining a decompse and backsolve call
//@note : here is the decision tree that matlab uses to solve dense systems of linear equations
          /*
            if(square) {
            	if([permuted] triangular) {
            		[p]tri
            	} else {
            		if(hermitian) {// conjugate transpose == self
            			if(diagonal is single sign)
            				try {cholesky} catch (...) {ldl}
            			else
            				ldl
            		} else {//not symmetric
            			if(upper hessenberg)
            				hess(not implemeneted)
            			else
            				lu
            		}
            	}
            } else {
            	qr
            }
          */
//@note: the following solvers aren't implemented here (but lu can be used for all of them with a modest speed penalty):
//          [p]tri - trivial
//          hess   - not a decomposition, if you want to implemented the solver see
//                   Henry, G. (1994). The shifted hessenberg system solve computation. Cornell Theory Center, Cornell University.
//          ldl    - this is tricky to implement (unlike cholesky stability requires [2x2] pivoting), see
//                   Ashcraft, C., R.G. Grimes, and J.G. Lewis. “Accurate Symmetric Indefinite Linear Equations Solvers.” SIAM J. Matrix Anal. Appl. Vol. 20. Number 2, 1998, pp. 513–561.
namespace solve {
	//@brief        : solve the linear system of equations A * x = b for x
	//@param a      : square matrix A in row major order, destroyed during function (over written with LU)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	template <typename Real> void lu(Real * const a, Real * const x, Real const * const b, const size_t n);

	//@brief        : solve the linear system of equations A * x = b for x
	//@param a      : square matrix A in row major order (only upper triangle is read, lower triangle is destroyed during function, overwritten with subdiagonal of L)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : matrix A must be symmetric positive-definite (real symmetric is a subset)
	template <typename Real> void cholesky(Real * const a, Real * const x, Real const * const b, const size_t n);
}

//@brief: functions to decompose a matrix
namespace decompose {
	//@brief        : perform in place LU decomposition
	//@param a      : square matrix A in row major order
	//@param p      : location to write pivots (compressed permutation matrix) [must hold at least n elements]
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : input is destroyed during function (over written with LU)
	template <typename Real> void lu(Real * const a, size_t * const p, const size_t n);

	//@brief        : perform semi in place cholesky decomposition
	//@param a      : square matrix A in row major order
	//@param d      : location to write diagonal of L
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@return       : true if matrix was negated to decompose (all negative diagonal)
	//@note         : only upper triangle is read
	//@note         : sub diagonal components are destroyed during function (overwritten with subdiagonal of L)
	template <typename Real> bool cholesky(Real * const a, Real * const d, const size_t n);

	//@brief        : perform in place qr decomposition
	//@param a      : rectangular matrix a A in row major order
	//@param m      : number of rows
	//@param n      : number of columns
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : upper triangle of A is overwritten with R
	//                remainder of A is overwritten with Q in a compressed format (sequence of givens rotations)
	//@note         : Givens rotation algorithm 5.2.2 in section 5.2.3 of Golub and Van Loan (1996) Matrix Computations
	template <typename Real> void qr(Real * const a, const size_t m, const size_t n);
}

//@brief: functions to solve a linear system of equations using an existing decomposition
//@note : may be preferred over solve:: functions if the same matrix A needs to be solved for multiple x/b values
namespace backsolve {
	//@brief        : solve the linear system of equations A * x = b for x using LU decomposition of A (from decompose::lu)
	//@param lu     : lu decomposition of A
	//@param p      : pivots (compressed permutation matrix)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : input is destroyed during function (over written with LU)
	template <typename Real> void lu(Real * const lu, size_t const * const p, Real * const x, Real const * const b, const size_t n);

	//@brief        : solve the linear system of equations A * x = b for x using cholesky decomposition of A (from decompose::cholesky)
	//@param ll     : subdiagonal of cholesky decomposition of A in row major order
	//@param d      : diagonal of cholesky decomposition of A
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@param neg    : true/false if A did / didn't need to be negated for decomposition
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : only upper triangle is read
	//@note         : sub diagonal components are destroyed during function (overwritten with subdiagonal of L)
	template <typename Real> void cholesky(Real const * const ll, Real const * const d, Real * const x, Real const * const b, const size_t n, const bool neg);
}

//@brief: additional functions to use the result of qr
namespace qr {
	//@brief        : compute the product of Q with another matrix y
	//@param a      : qr decomposition of a (from decompose::qr)
	//@param y      : matrix to multiply Q by
	//@param m      : number of rows in a (and y)
	//@param n      : number of columns in a
	//@param p      : number of columns in y
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : y is overwritten with Q * y
	//@note         : Q can be 'decompressed' by calling this function with y == Identity(m) ( && p == m)
	template <typename Real> void applyQ(Real const * const a, Real * const y, const size_t m, const size_t n, const size_t p);

	//@brief        : compute the product of Q^H (==Q^-1) with another matrix y
	//@param a      : qr decomposition of a (from decompose::qr)
	//@param y      : matrix to multiply Q^H by
	//@param m      : number of rows in a (and y)
	//@param n      : number of columns in a
	//@param p      : number of columns in y
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : y is overwritten with Q^H * y
	template <typename Real> void applyQH(Real const * const a, Real * const y, const size_t m, const size_t n, const size_t p);
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <limits>
#include <vector>
#include <complex>
#include <stdexcept>
#include <numeric>
#include <algorithm>
#include <functional>
#include <vector>
#include <type_traits>

//@brief: helper functions to handle real and complex valued decompositions with a single function
namespace detail {
	//helper to check if a type is std::complex<another type>
	template<class Real> struct is_complex                      : std::false_type {};//default to 'naked' real number
	template<class Real> struct is_complex<std::complex<Real> > : std::true_type  {};//handle std::complex<real number> specially, can get base type as typename Real::value_type

	//shortcuts for enable_if
	template<typename Real>	using IfReal = typename std::enable_if<!is_complex<Real>::value, int>::type;
	template<typename Real>	using IfCplx = typename std::enable_if< is_complex<Real>::value, int>::type;

	//helper to typedef real value from a real or complex template
	template <typename T, typename Enable=void> struct ReVal;
	template <typename R> struct ReVal<R, typename std::enable_if<!is_complex<R>::value >::type> {typedef          R             type;};//ReVal<Real>::type == Real
	template <typename C> struct ReVal<C, typename std::enable_if< is_complex<C>::value >::type> {typedef typename C::value_type type;};//ReVal<Cplx>::type == typename Cplx::value_type

	//helper to check if real part of number is negative
	template<class Real, IfReal<Real> = 0> inline          bool             signbit(const Real v) {return std::signbit(v       );}
	template<class Cplx, IfCplx<Cplx> = 0> inline          bool             signbit(const Cplx v) {return std::signbit(v.real());}

	//helper to get machine epsilon
	template<class Real, IfReal<Real> = 0> inline          Real             eps    ()             {return std::numeric_limits<         Real            >::epsilon();}
	template<class Cplx, IfCplx<Cplx> = 0> inline typename Cplx::value_type eps    ()             {return std::numeric_limits<typename Cplx::value_type>::epsilon();}

	//helper to return conjugate of real number os real (instead of complex)
	template<class Real, IfReal<Real> = 0> inline          Real             conj   (const Real v) {return           v ;}
	template<class Cplx, IfCplx<Cplx> = 0> inline          Cplx             conj   (const Cplx v) {return std::conj(v);}
}

namespace qr {
	//compute, compress, and decompress a givens rotation
	namespace givens {
		//@brief: compute the givens rotation to zero one element of a matrix
		//@param va: entry to rotate into
		//@param vb: entry to rotate from (will be zeroed)
		//@param c : location to write cosine term of givens rotation
		//@param s : location to write   sine term of givens rotation
		//@return  : true if a rotation is needed, false if |vb| is already ~0
		template <typename Real, detail::IfReal<Real> = 0>
		bool compute(const Real& va, const Real& vb, Real& c, Real& s) {
			//this method works for complex numbers as well if va/vb are replaces with |va|, |vb| but isn't ideal
			if(std::sqrt(std::norm(vb)) < detail::eps<Real>()) {//vb is already 0
				c = Real(1);
				s = Real(0);
				return false;
			} else {//vb is nonzero
				const Real h = std::hypot(va, vb);//vb is nonzero so no chance to divide by 0
				c = va / h;//cos of rotation to zero vb
				s = vb / h;//sin of rotation to zero vb
				return true;
			}
		}

		//@brief: compute the givens rotation to zero one element of a matrix
		//@param va: entry to rotate into
		//@param vb: entry to rotate from (will be zeroed)
		//@param c : location to write cosine term of givens rotation
		//@param s : location to write   sine term of givens rotation
		//@return  : true if a rotation is needed, false if |vb| is already ~0
		template <typename Cplx, detail::IfCplx<Cplx> = 0>
		bool compute(const Cplx& va, const Cplx& vb, Cplx& c, Cplx& s) {
			//a complex givens rotation to zero vb into va can be computed analogously to a real givens rotation:
			// h = sqrt( |va|^2 + |vb|^2 )
			// c = va / h
			// s = vb / h
			//if va and vb are considered as complex numbers va == ma * exp(I * pa) where ma/mb are the magnitudes and pa/pb the phases
			//then the givens rotation is essentially a combination of 3 rotations
			// c     = cos(theta) * exp(I * -pa)
			// s     = sin(theta) * exp(I * -pb)
			// theta = atan( |mb| / |ma| )
			//this is the easiest option but leaves us 3 degrees of freedom to store as 2 numbers for an in place decomposition
			//also the resulting givens matrix isn't hermetian (since c.imag() is non-zero)
			//the givens rotation can also be considered as a combination of 2 rotations e.g.
			// c   = cos(theta)
			// s   = sin(theta) * exp(I * phi)
			// phi = pa - pb
			//this also works for real numbers its just slightly more expensive (and confusing)
			const typename Cplx::value_type mb = std::sqrt(std::norm(vb));//get magnitude of element to zero
			if(mb < detail::eps<Cplx>()) {//vb is already 0 (theta == 0)
				c = Cplx(1);
				s = Cplx(0);
				return false;
			} else {//vb is nonzero
				const typename Cplx::value_type ma = std::sqrt(std::norm(va));//get magnitude of element to zero into
				if(ma < detail::eps<Cplx>()) {//don't divide by 0 (theta == pi/2)
					c = typename Cplx::value_type(0);//cos(theta)
					s = detail::conj(vb) / mb;
				} else {
					typename Cplx::value_type st(1);
					const typename Cplx::value_type ba = mb / ma;//we need to check that vb isn't already zero to avoid divide by 0 here
					c  = st / std::hypot(ba, st);//cos(theta)
					st = c.real() * ba;          //sin(theta)
					s = (va * detail::conj(vb)) * ( st / (ma * mb) );// cos/sin(phi) * sin(theta)
				}
				return true;
			}
		}

		//@brief  : compress a givens rotation into single number 
		//@param c: cos(theta)
		//@param s: sin(theta)
		//@return : 'compressed' representation of c and s
		//@note   : based on section 5.1.11 of Golub and Van Loan (1996) Matrix Computations
		template <typename Real, detail::IfReal<Real> = 0>
		Real compress(const Real& c, const Real& s) {
			if(std::fabs(c) < std::numeric_limits<Real>::epsilon()) {//don't divide by zero
				return Real(std::signbit(s) ? -1 : 1);//cos == 0 -> store sign of sine
			} else if(std::signbit(c) == std::signbit(s)) {//0 < theta < 90 or 180 < theta < 270
				return s / Real(2);//save sin / 2 -> [-0.5, 0.0), (0.0,0.5]
			} else {//90 < theta < 180 or 270 < theta < 360
				return Real(2) / c;//save 2 / c (-inf, -2], [2, inf)
			}
		}

		//@brief  : compress a givens rotation into single number 
		//@param c: cos(theta)
		//@param s: sin(theta) * exp(I * phi)
		//@return : 'compressed' representation of c and s
		//@note   : modified to handle complex givens rotations
		template <typename Cplx, detail::IfCplx<Cplx> = 0>
		Cplx compress(const Cplx& c, const Cplx& s) {
			const Cplx ePhi = s / std::sqrt(typename Cplx::value_type(1) - c.real() * c.real());//exp(I * phi) * sign( cos(theta) )
			return Cplx(c.real(), compress<typename Cplx::value_type>(ePhi.real(), ePhi.imag()));//leave s as is, compress cos/sin(phi) using real valued function
		}

		//@brief  : decompress a givens rotation back into a sin/cos pair
		//@param p: compressed rotation (output of compress)
		//@param c: location to write cos(theta)
		//@param s: location to write sin(theta)
		//@note   : based on section 5.1.11 of Golub and Van Loan (1996) Matrix Computations
		template <typename Real, detail::IfReal<Real> = 0>
		void decompress(const Real& p, Real& c, Real& s) {
			const Real ap = std::fabs(p);
			if(Real(1) == ap) {// +/-1 means we stored sine
				c = Real(0);
				s = p;
			} else if(ap < Real(1)) {//p is s/2 and signs match
				s = Real(2) * p;//extract sine
				c = std::copysign(std::sqrt(Real(1) - s * s), p);//mag from s^2+c^2==1 and sign from p (match)
			} else {//p is 2/c and signs don't match
				c = Real(2) / p;//extract cosine
				s = std::copysign(std::sqrt(Real(1) - c * c),-p);//mag from s^2+c^2==1 and sign from p (mismatch)
			}
		}

		//@brief  : decompress a givens rotation
		//@param p: compressed rotation (output of compress)
		//@param c: location to write cos(theta)
		//@param s: location to write sin(theta) * exp(I * phi)
		//@note   : modified to handle complex givens rotations
		template <typename Cplx, detail::IfCplx<Cplx> = 0>
		void decompress(const Cplx& p, Cplx& c, Cplx& s) {
			typename Cplx::value_type sr, si;
			c = Cplx(p.real(), typename Cplx::value_type(0));//c is already real valued
			decompress<typename Cplx::value_type>(p.imag(), sr, si);//extract exp(I * phi) * sign( sin(theta) )
			s = Cplx(sr, si) * std::sqrt(typename Cplx::value_type(1) - c.real() * c.real());//multiply by |sin(theta)|
		}
	}
}

namespace solve {
	//@brief        : solve the linear system of equations A * x = b for x
	//@param a      : square matrix A in row major order, destroyed during function (over written with LU)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	template <typename Real> void lu(Real * const a, Real * const x, Real const * const b, const size_t n) {
		std::vector<size_t> p(n);//allocate storage for pivots
		decompose::lu(a, p.data(), n);//decompose in place
		backsolve::lu(a, p.data(), x, b, n);//backsolve
	}

	//@brief        : solve the linear system of equations A * x = b for x
	//@param a      : square matrix A in row major order (only upper triangle is read, lower triangle is destroyed during function, overwritten with subdiagonal of L)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : matrix A must be symmetric positive-definite (real symmetric is a subset)
	template <typename Real> void cholesky(Real * const a, Real * const x, Real const * const b, const size_t n) {
		std::vector<Real> d(n);//allocate memory to hold diagonal (so we don't have to overwrite a's diagonal)
		const bool neg = decompose::cholesky(a, d.data(), n);//decompose (semi) in place
		backsolve::cholesky(a, d.data(), x, b, n, neg);//backsolve
	}
}

namespace decompose {
	//@brief        : perform in place LU decomposition
	//@param a      : square matrix A in row major order
	//@param p      : location to write pivots (compressed permutation matrix) [must hold at least n elements]
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : input is destroyed during function (over written with LU)
	template <typename Real> void lu(Real * const a, size_t * const p, const size_t n) {
		static_assert(std::is_floating_point<typename detail::ReVal<Real>::type>::value, "solve::lu must be tempalted on either a real or complex<real> type");
		//compute 1 / max(|row|) for implicitly scaled pivoting
		std::vector<Real> s(n, 1);//scaling factor
		for(size_t i=0; i<n; i++) s[i] /= std::real(*std::max_element(a+i*n, a+(i+1)*n, [](const Real& a,const Real& b){return std::abs(a)<std::abs(b);}));//scaling factor for each row is 1 / max(abs(row))
		std::iota(p, p + n, 0);//initially no permutation

		//loop down diagonal performing decomposition
		for(size_t i = 0; i < n; i++) {
			//determine pivot
			size_t iMax = i;//start with current row
			for(size_t j = i+1; j < n; j++)//loop over remaining rows
				if( std::norm( s[i] * a[i*n+i] ) < std::norm( s[j] * a[j*n+i] ) )//check if this row is bigger
					iMax =  j;//save biggest row index

			//don't divide by 0
			if(std::sqrt(std::norm( s[iMax] * a[iMax*n+i] )) < detail::eps<Real>())//if the biggest number is too small the matrix is singular
				throw std::runtime_error("singular matrix");//stop! (could handle but why bother if we can't back solve)

			//pivot if needed
			if(i != iMax) {
				std::swap_ranges(a+i*n, a+(i+1)*n, a+iMax*n);//swap row n with row iMax
				std::swap(s[i], s[iMax]);//swap row scalings
				std::swap(p[i], p[iMax]);//update permutation
			}

			//reduce
			for(size_t j = i+1; j < n; j++) {
				a[j*n+i] /= a[i*n+i];
				const Real& aji = a[j*n+i];//we'll need this value a bunch of times
				for(size_t k = i+1; k < n; k++) a[j*n+k] -= aji * a[i*n+k];
			}
		}
	}

	//@brief        : perform semi in place cholesky decomposition
	//@param a      : square matrix A in row major order
	//@param d      : location to write diagonal of L
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@return       : true if matrix was negated to decompose (all negative diagonal)
	//@note         : only upper triangle is read
	//@note         : sub diagonal components are destroyed during function (overwritten with subdiagonal of L)
	template <typename Real> bool cholesky(Real * const a, Real * const d, const size_t n) {
		static_assert(std::is_floating_point<typename detail::ReVal<Real>::type>::value, "cholesky must be templated on either a real or complex<real> type");
		if(n == 0) return false;//handle empty matrix
		const bool neg = detail::signbit(a[0]);//try negating matrix to handle negative definite matrices
		for(size_t i = 0; i < n; i++) {//loop over rows
			if(detail::signbit(a[n * i + i]) != neg) throw std::runtime_error("cholesky decomposition requires all positive or negative diagonal");//need to use LDL instead (which is why we've kept the diagonal intact)
			for(size_t j = i; j < n; j++) {//loop over columns
				Real sum = neg ? -detail::conj(a[n * i + j]) : detail::conj(a[n * i + j]);//get A_{i,j} (or negative if needed)
				for(size_t k = 0; k < i; k++) sum -= detail::conj(a[i * n + k]) * a[j * n + k];//accumulate A_{i,j} - \sum_{k=0}^{j-1} L_{i,k} * L^*_{j,k}
				if(i == j) {//we're computing L_{j,j} (save in d instead of overwriting elements of A)
					//sum will always be real here since A_{i,i} is real and L_{i,k} * L^*_{j,k} -> L_{j,k} * L^*_{j,k} for i == j
					if(std::real(sum) < detail::eps<Real>()) throw std::runtime_error("cholesky decomposition failed");//need to use LDL instead (which is why we've kept the diagonal intact)
					d[i] = std::sqrt(sum);//L_{j,j} = \sqrt{ A_{j,j} - \sum_{k=0}^{j-1} L_{j,k} * L_{j,k} }
				} else {//we're computing L_{i,j}
					a[j * n + i] = sum / d[i];//L_{i,j} = ( A_{i,j} - \sum_{k=0}^{j-1} L_{i,k} * L_{j,k} ) / L_{j,j}
				}
			}
		}
		return neg;
	}

	//@brief        : perform in place qr decomposition
	//@param a      : rectangular matrix a A in row major order
	//@param m      : number of rows
	//@param n      : number of columns
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : upper triangle of A is overwritten with R
	//                remainder of A is overwritten with Q in a compressed format (sequence of givens rotations)
	//@note         : Givens rotation algorithm 5.2.2 in section 5.2.3 of Golub and Van Loan (1996) Matrix Computations
	template <typename Real> void qr(Real * const a, const size_t m, const size_t n) {
		for(size_t i = 0; i < n; i++) {//loop over columns (left to right)
			for(size_t j = i+1; j < m; j++) {//loop down rows (top to bottom)
				const Real& vb = a[j*n+i];//this is the element we're zeroing
				const Real& va = a[i*n+i];//this is the element we're rotating with (I chose the diagonal but that is a bit arbitrary)
				Real s, c;//location to save givens rotation
				if(qr::givens::compute(va, vb, c, s)) {//get the givens rotation to zero out vb (and check if vb is already 0)
					for(size_t k = i; k < n; k++) {//loop over row applying rotation (skip columns that are already zeroed)
						const Real ak = a[i*n+k];//element from row we're rotating into
						const Real bk = a[j*n+k];//element from row rotating out of
						a[i*n+k] =              s  * bk +              c  * ak;//apply rotation
						a[j*n+k] = detail::conj(c) * bk - detail::conj(s) * ak;//apply rotation
					}
				}
				a[j*n+i] = qr::givens::compress(c, s);//save the rotation (accumulate elements of Q)
			}
		}
	}
}

//@brief: functions to solve a linear system of equations using an existing decomposition
//@note : may be preferred over solve:: functions if the same matrix A needs to be solved for multiple x/b values
namespace backsolve {
	//@brief        : solve the linear system of equations A * x = b for x using LU decomposition of A (from decompose::lu)
	//@param lu     : lu decomposition of A
	//@param p      : pivots (compressed permutation matrix)
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : input is destroyed during function (over written with LU)
	template <typename Real> void lu(Real * const lu, size_t const * const p, Real * const x, Real const * const b, const size_t n) {
		for(size_t i = 0; i < n; i++) x[i] = b[p[i]];//permute b
		for(size_t i = 0; i < n; i++) x[i] -= std::inner_product(x, x+i, lu+i*n, Real(0));//solve L y = b for y
		for(size_t i = n-1; i < n; i--) x[i] = (x[i] - std::inner_product(x+i+1, x+n, lu+i*n+i+1, Real(0))) / lu[i*n+i];//solve U x = y for x
	}

	//@brief        : solve the linear system of equations A * x = b for x using cholesky decomposition of A (from decompose::cholesky)
	//@param ll     : subdiagonal of cholesky decomposition of A in row major order
	//@param d      : diagonal of cholesky decomposition of A
	//@param x      : location to write x
	//@param b      : column vector b
	//@param n      : size of matrix
	//@param neg    : true/false if A did / didn't need to be negated for decomposition
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : only upper triangle is read
	//@note         : sub diagonal components are destroyed during function (overwritten with subdiagonal of L)
	template <typename Real> void cholesky(Real const * const ll, Real const * const d, Real * const x, Real const * const b, const size_t n, const bool neg) {
		for(size_t i = 0; i < n; i++) x[i] = (b[i] - std::inner_product(x, x + i, ll + i * n, Real(0))) / d[i];//solve L y = b for y with back substitution
		for(size_t i = n-1; i < n; i--) {//solve L^* x = y for x with forward substitution
			for(size_t j = n-1; j != i; --j) x[i] -= detail::conj(ll[j * n + i]) * x[j];
			x[i] /= d[i];
		}
		if(neg) std::transform(x, x + n, x, std::negate<Real>());//undo negation if needed
	}
}

//@brief: additional functions to use the result of qr
namespace qr {
	//@brief        : compute the product of Q with another matrix y
	//@param a      : qr decomposition of a (from decompose::qr)
	//@param y      : matrix to multiply Q by
	//@param m      : number of rows in a (and y)
	//@param n      : number of columns in a
	//@param p      : number of columns in y
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : y is overwritten with Q * y
	//@note         : Q can be 'decompressed' by calling this function with y == Identity(m) ( && p == m)
	template <typename Real> void applyQ(Real const * const a, Real * const y, const size_t m, const size_t n, const size_t p) {
		for(size_t i = n-1; i < n; i--) {//loop over columns (right to left) [reverse order from calculation]
			for(size_t j = m-1; j > i; j--) {//loop up rows (bottom to top) [reverse order from calculation]
				Real c, s;//storage for components of rotation
				qr::givens::decompress(a[j*n+i], c, s);//extract rotation that was applied
				for(size_t k = 0; k < p; k++) {//loop over columns of y
					const Real ak = y[j*p+k];//element from current  row [swapped from calculation]
					const Real bk = y[i*p+k];//element from diagonal row [swapped from calculation]
					y[j*p+k] = detail::conj(s) * bk + detail::conj(c) * ak;//apply rotation
					y[i*p+k] =              c  * bk -              s  * ak;//apply rotation
				}
			}
		}
	}

	//@brief        : compute the product of Q^H (==Q^-1) with another matrix y
	//@param a      : qr decomposition of a (from decompose::qr)
	//@param y      : matrix to multiply Q^H by
	//@param m      : number of rows in a (and y)
	//@param n      : number of columns in a
	//@param p      : number of columns in y
	//@template Real: matrix scalar type, must be a floating point type (real or std::complex)
	//@note         : y is overwritten with Q^H * y
	//@warning      : this is currently actually Q^T but for complex number the ivners of Q is Q^* (conjugate transpose)
	template <typename Real> void applyQH(Real const * const a, Real * const y, const size_t m, const size_t n, const size_t p) {
		for(size_t i = 0; i < n; i++) {//loop over columns (left to right) [same order as calculation]
			for(size_t j = i+1; j < m; j++) {//loop down rows (top to bottom) [same order as calculation]
				Real c, s;//storage for components of rotation
				qr::givens::decompress(a[j*n+i], c, s);//extract rotation that was applied
				for(size_t k = 0; k < p; k++) {//loop over columns of y
					const Real ak = y[(i-0)*p+k];//element from diagonal row [same as calculation]
					const Real bk = y[(j  )*p+k];//element from current  row [same as calculation]
					y[(i-0)*p+k] =              s  * bk +              c  * ak;//apply rotation
					y[(j  )*p+k] = detail::conj(c) * bk - detail::conj(s) * ak;//apply rotation
				}
			}
		}
	}
}

#endif//_LINALG_H_
