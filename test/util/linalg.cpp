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
#include <complex>

//@brief : test real valued linalg functions
//@return: true / false if tests pass / fail
template <typename Real> bool testLinAlg(std::ostream& os);

int main() {
	std::ostream& os = std::cout;
	try {
		os << "testing         float  \n";
		if(!testLinAlg<              float   >(os)) return EXIT_FAILURE;
		os << "testing         double \n";
		if(!testLinAlg<              double  >(os)) return EXIT_FAILURE;
	/*this doesn't currently pass (and isn't needed for EMSphInx)
		os << "testing complex<float >\n";
		if(!testLinAlg< std::complex<float > >(os)) return EXIT_FAILURE;
		os << "testing complex<double>\n";
		if(!testLinAlg< std::complex<double> >(os)) return EXIT_FAILURE;
	*/
	} catch (std::exception& e) {
		os << "caught: " << e.what() << '\n';
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

#include "util/linalg.hpp"

//short cut since we'll be using this alot
template <class T> using ReVal = typename detail::ReVal<T>::type;

#include <random>
#include <sstream>
#include <algorithm>
#include <complex>
#include <vector>
#include <iomanip>
#include <ios>

//@brief    : helper to build a random real/complex number
//@param mag: magnitude
//@param gen: random generator
//@return   : mag
template<typename Real, class URNG, detail::IfReal<Real> = 0>
Real rndPhs(const ReVal<Real>& mag, URNG& gen) {return mag;}

//@brief    : helper to build a random real/complex number
//@param mag: magnitude
//@param gen: random generator
//@return   : mag with random phase
template<typename Cplx, class URNG, detail::IfCplx<Cplx> = 0>
Cplx rndPhs(const ReVal<Cplx>& mag, URNG& gen) {
	static std::uniform_real_distribution<typename Cplx::value_type> dis(0, (ReVal<Cplx>)6.2831853071795864769252867665590);
	return std::polar(mag, dis(gen));
}

//lets make life a little easier by defining a bare bones matrix class
//this isn't super efficient (but it doesn't need to be)
template <typename Real> class Mat {
	size_t _m, _n;//rows, cols
	std::vector<Real> _a;//data in row major order

	public:
		//@brief  : construct an empty matrix (all zeros)
		//@param r: number of rows
		//@param c: number of columns
		Mat(const size_t r, const size_t c) : _m(r), _n(c), _a(r*c, Real(0)) {}

		//@brief: construct an empty matrix
		Mat() : Mat(0, 0) {}

		//@brief  : construct the identity matrix
		//@param n: side length
		static Mat Ident(const size_t n) {Mat a(n, n); for(size_t i = 0; i < n; i++) a[i][i] = Real(1); return a;}

		//@brief  : get access to the rth row of data
		//@param r: row to get
		//@return : pointer to row start
		Real const * operator[](const size_t r) const {return _a.data() + _n * r;}
		Real       * operator[](const size_t r)       {return _a.data() + _n * r;}

		//@brief : get the # of rows/columns in this matrix
		//@return: # of rows/columns
		size_t rows() const {return _m;}
		size_t cols() const {return _n;}

		//@brief: get underlying storage (row major order)
		Real const * data() const {return _a.data();}
		Real       * data()       {return _a.data();}

		//@brief : access to underlying storage iterators
		//@return: iterators for row major order storage
		typename std::vector<Real>::const_iterator cbegin() const {return _a.cbegin();} 
		typename std::vector<Real>::const_iterator cend  () const {return _a.cend  ();} 
		typename std::vector<Real>::      iterator  begin()       {return _a. begin();} 
		typename std::vector<Real>::      iterator  end  ()       {return _a. end  ();} 

		//@brief    : resize a matrix, padding with 0 if needed
		//@param r/c: new number of rows/columns
		void resize (const size_t r, const size_t c) {setRows(r); setCols(c);}
		void setRows(const size_t r) {_m = r; _a.resize(_m * _n, Real(0));}
		void setCols(const size_t c) {
			if(c < _n) {//shrinking
				for(size_t i = 1; i < _m; i++) std::copy(_a.begin() + i * _n, _a.begin() + i * _n + c, _a.begin() + i * c);//loop over rows repacking
				_a.resize(c * _m);//resize
			} else if(c > _n) {//expanding
				_a.resize(c * _m);//resize
				for(size_t i = _m; i-- > 0 ; ) {
					std::copy_backward(_a.begin() + i * _n, _a.begin() + i * _n + _n, _a.begin() + i * c + _n);//repack (this does nothing for i==0)
					std::fill(_a.begin() + i * c + _n, _a.begin() + i * c + c, Real(0));//fill new space with zeros
				}
			}
			_n = c;
		}

		//@brief  : multiply 2 matrices together
		//@param a: other matrix to multiply by
		//@return : this * a
		Mat operator*(const Mat& a) const {
			if(_n != a._m) throw std::runtime_error("matrix mutliplication size mismatch");//sanity check shapes
			const size_t r = _m;
			const size_t c = a._n;
			const size_t p = _n;
			Mat mat(r, c);
			for(size_t i = 0; i < r; i++) {//loop over output rows
				for(size_t j = 0; j < c; j++) {//loop over output columns
					Real v(0);//initialize value
					for(size_t k = 0; k < p; k++) v += _a[i*p+k] * a._a[k*c+j];//compute sum
					mat._a[i*c+j] = v;//save value
				}
			}
			return mat;
		}

		//@brief : transpose this matrix
		//@return: this^T
		Mat transpose() const {
			Mat<Real> res(_n, _m);
			for(size_t i = 0; i < _n; i++) {
				for(size_t j = 0; j < _m; j++) {
					res[i][j] = operator[](j)[i];
				}
			}
			return res;
		}

		//@brief : conjugate transpose this matrix
		//@return: this^*
		Mat conjugateTranspose() const {
			Mat<Real> res(_n, _m);
			for(size_t i = 0; i < _n; i++) {
				for(size_t j = 0; j < _m; j++) {
					res[i][j] = detail::conj(operator[](j)[i]);
				}
			}
			return res;
		}

		//@brief  : compare this matrix to another
		//@param a: matrix to compare to
		//@return : L2-norm sqrt ( ||this - a||^2 ) if the a is the same shape as this, largest possible Real otherwise
		ReVal<Real> compare(const Mat& a) {
			if(a._m != _m || a._n != _n) return std::numeric_limits<ReVal<Real> >::max();// matrices must have same size
			ReVal<Real> res(0);
			for(size_t i = 0; i < _a.size(); i++) res += std::norm(_a[i] - a._a[i]);
			return std::sqrt(res);
		}

		//@brief    : convert to string representation
		//@param os : ostream to print to
		//@param sci: should value be printed with scientific (instead of fixed) notation
		//@return   : os
		std::ostream& print(std::ostream& os, const bool sci = false) const {
			os << (sci ? std::scientific : std::fixed) << std::setprecision(4);
			for(int i = 0; i < _m; i++) {
				for(int j = 0; j < _n; j++) {
					os << std::setw(sci ? 13 : 10) << _a[i*_n+j];
				}
				os << '\n';
			}
			os << '\n';
			// os << std::defaultfloat; this isn't supported until fairly recent versions of gcc
			os.unsetf(std::ios_base::floatfield);
			return os;
		}

		friend std::ostream& operator<<(std::ostream& os, const Mat& m) {return m.print(os, false);}

};

//@brief  : compare 2 vectors
//@param a: first vector to compare
//@param b: second vector to compare
//@return : |a-b| / min(|a|, |b|)
template <typename Real>
ReVal<Real> vectorCompare(const std::vector<Real>& a, const std::vector<Real>& b) {
	if(a.size() != b.size()) return std::numeric_limits<ReVal<Real> >::max();// vectors must have same length
	ReVal<Real> ma(0), mb(0), mab(0);
	for(size_t i = 0; i < a.size(); i++) {
		ma  += std::norm(a[i]);
		mb  += std::norm(b[i]);
		mab += std::norm(a[i] - b[i]);
	}
	return std::sqrt(mab / std::max(ma, mb));
}

////////////////////////////////////////////////////////////////////////
//                      Random Matrix Generation                      //
////////////////////////////////////////////////////////////////////////

//@brief    : generate a random orthogonal matrix
//@param n  : side length
//@param gen: random generator
//@return   : nxn random matrix
//@reference: Stewart, G. W. (1980). The efficient generation of random orthogonal matrices with an application to condition estimators. SIAM Journal on Numerical Analysis, 17(3), 403-409.
template<typename Real, class URNG, detail::IfReal<Real> = 0>
Mat<Real> randomQ(const size_t n, URNG& gen) {
	Mat<Real> a = Mat<Real>::Ident(n);//start with identity matrix
	std::normal_distribution<ReVal<Real>> dis(ReVal<Real>(0), ReVal<Real>(1));//we'll also need a random distribution

	//to produce a random orthogonal (n+1) x (n+1) matrix:
	// -start with a random orthogonal  n x n matrix
	// -construct an (n+1)^2 matrix from the n^2 by augmenting with a 1 in the bottom right
	// -generate a normal random (n+1) unit vector
	// -apply a householder reflection to the (n+1)^2 matrix using the random vector
	//to build our random matrix we start from a random 1x1 orthogonal matrix (i.e. [1]) and work up to n x n
	std::vector<Real> x(n), row(n);//work space
	for(size_t i = 1; i < n; i++) {//loop over size we want to augment to
		//generate a random unit vector with normally distributed elements
		Real mag2 = 0;
		for(size_t j = 0; j <= i; j++) {
			const Real r = rndPhs<Real>(dis(gen), gen);//draw a random number
			mag2 += r * r;//accumulate magnitude
			x[j] = r;//save value
		}

		//apply house holder transform row by row such that a *= h
		for(size_t r = 0; r <= i; r++) {//loop over rows of a
			const Real vh = x[r] * ReVal<Real>(2) / mag2;//we'll need this value multiple times
			std::fill(row.begin(), row.begin() + i+1, Real(0));//fill output row with 0
			for(size_t c = 0; c <= i; c++) {//loop over cols of a
				Real& v = row[c];//get output value
				for(size_t k = 0; k <= i; k++) {//loop down rows of h
					const Real h = ReVal<Real>(c == k ? 1 : 0) - x[c] * x[k] * ReVal<Real>(2) / mag2;//compute value of householder matrix
					v += h * a[r][k];//accumulate matrix product
				}
			}
			std::copy(row.begin(), row.begin() + i+1, a[r]);//copy output back over a
		}
	}
	return a;
}

//@brief    : generate a random unitary matrix
//@param n  : side length
//@param gen: random generator
//@return   : nxn random matrix
//@reference: Mezzadri, F. (2007). How to generate random matrices from the classical compact groups. arXiv [arXiv:math-ph/0609050]
template<typename Cplx, class URNG, detail::IfCplx<Cplx> = 0>
Mat<Cplx> randomQ(const size_t n, URNG& gen) {
	//generate a random matrix
	Mat<Cplx> q  = Mat<Cplx>::Ident(n);//start with identity matrix
	Mat<Cplx> qr(q);//make a copy of the same size
	std::normal_distribution<ReVal<Cplx>> dis(ReVal<Cplx>(0), ReVal<Cplx>(1));//we'll also need a random distribution
	for(Cplx& v : qr) v = rndPhs<Cplx>(dis(gen), gen);//fill with random normally distributed numbers

	//now do the QR decomposition and extract Q
	decompose::qr(qr.data(), n, n);
	qr::applyQ(qr.data(), q.data(), n, n, n);

	//let D == diag(R) / element wise abs(diag(R))
	//simultaneously extract D and compute D * Q
	for(size_t j = 0; j < n; j++) {
		const Cplx d = qr[j][j] / std::sqrt(std::norm(qr[j][j]));//compute d/|d|
		for(size_t i = 0; i < n; i++) qr[j][i] = q[j][i] * d;//fill row with d * q
	}

	// random Q == q * d * q
	return q * qr;
}

//@brief    : generate a random square matrix with the specified eigen values and random eigen vectors
//@param d  : eigen values
//@return   : nxn random matrix
//@param gen: random generator
//@note     : a very slight modification would enable random rectangular matrices from specified singular values
template<typename Real, class URNG>
Mat<Real> randomEigVec(const std::vector<Real> d, URNG& gen) {
	//handle trivial case
	if(d.size() <= 1) {
		Mat<Real> res = Mat<Real>::Ident(d.size());
		if(1 == d.size()) res[0][0] = d[0];
		return res; 
	}

	//generate random eigen vectors and assmble
	Mat<Real> q = randomQ<Real>(d.size(), gen);
	Mat<Real> qt = q.transpose();
	for(size_t r = 0; r < d.size(); r++) std::transform(d.begin(), d.end(), q[r], q[r], std::multiplies<Real>());//q *= d;
	return q * qt;
}

//@brief    : generate random eigen values with the specified condition number and definite-ness
//@param n  : number of eigen values
//@param c  : condition number
//@param gen: random generator
//@param def: should the values be all real positive, all real negative, or mixed (and potentially complex) [+1, -1, or 0 respectively]
template<typename Real, class URNG>
std::vector<Real> randomEigVal(const size_t n, const ReVal<Real> c, URNG& gen, const int def = 0) {
	//handle trivial case
	if(n < 2) return std::vector<Real>(n, 1);

	//start by sanity checking condition number and generating eigen value bounds
	if(c < ReVal<Real>(1)) throw std::runtime_error("cannot generate matrix with condition number <1");
	const ReVal<Real> vMax = ReVal<Real>(1);// |largest eigen value| (this could be anything, but 1 is nice)
	const ReVal<Real> vMin = vMax / c;// |smallest eigen value|

	//now build our distributions
	static std::uniform_real_distribution<ReVal<Real>> disSgn(ReVal<Real>(-1), ReVal<Real>(1));//random sign generation for def == 0
	std::uniform_real_distribution<ReVal<Real>> dis(vMin, vMax);//magnitude generation

	//generate random eigen values by from the distribution
	std::vector<Real> diag(n);
	switch(def) {
		case  0://random eigen values
			for(Real& v : diag) v = rndPhs<Real>(std::copysign(dis(gen), disSgn(gen)), gen);//random value in range * random sign (w/ random phase if complex)
			break;

		case  1: for(Real& v : diag) v =  dis(gen); break;;// random value in range * random sign
		case -1: for(Real& v : diag) v = -dis(gen); break;//-random value in range * random sign
	}

	//make sure we achieve the desired condition number
	switch(def) {
		case  0:
			diag[0] = rndPhs<Real>(std::copysign(vMin, disSgn(gen)), gen);//set upper magnitude bound
			diag[1] = rndPhs<Real>(std::copysign(vMax, disSgn(gen)), gen);//set lower magnitude bound
			break;

		case  1:
			diag[0] =  vMin;
			diag[1] =  vMax;
			break;

		case -1:
			diag[0] = -vMin;
			diag[1] = -vMax;
			break;
	}

	//randomly shuffle our eigenvalues and return
	std::shuffle(diag.begin(), diag.end(), gen);
	return diag;
}

//@brief    : generate a random square matrix with the specified condition number
//@param n  : side length
//@param c  : condition number
//@param gen: random generator
//@return   : nxn random matrix
template<typename Real, class URNG>
Mat<Real> randomMatrix(const size_t n, const ReVal<Real> c, URNG& gen, const bool pd = false) {return randomEigVec(randomEigVal<Real>(n, c, gen, 0), gen);}

//@brief    : generate a random positive or negative definite square matrix with the specified condition number
//@param n  : side length
//@param c  : condition number
//@param pos: true/false for positive/negative definite
//@param gen: random generator
//@return   : nxn random matrix
template<typename Real, class URNG>
Mat<Real> randomDefinite(const size_t n, const ReVal<Real> c, URNG& gen, const bool pos = true) {return randomEigVec(randomEigVal<Real>(n, c, gen, pos ? 1 : -1), gen);}

////////////////////////////////////////////////////////////////////////
//                 Tests for Different Decompositions                 //
////////////////////////////////////////////////////////////////////////

//@brief    : test LU decomposition for a random matrix matrix
//@param os : location to write error messages
//@param n  : matrix size
//@param cn : condition number to test for
//@param gen: random generator e.g. std::mt19937
//@return   : true/false if test passed/failed
template<typename Real, class URNG>
bool testLU(std::ostream& os, const size_t n, const ReVal<Real> cn, URNG& gen) {
	//create a random matrix
	Mat<Real> a = randomMatrix<Real>(n, cn, gen);

	//do in place decomposition on a copy
	Mat<Real> lu(a);//copy to keep a
	std::vector<size_t> p(n);//allocate permutation matrix
	decompose::lu(lu.data(), p.data(), n);//in place decompose

	//extract L and U
	Mat<Real> l(lu), u(lu);
	for(size_t r = 0; r < n; r++) {
		for(size_t c = 0; c < r; c++) {
			u[r][c] = Real(0);
		}
		l[r][r] = Real(1);
		for(size_t c = r+1; c < n; c++) {
			l[r][c] = Real(0);
		}
	}

	//compute p * a
	Mat<Real> pa(n, n);
	for(size_t i = 0; i < p.size(); i++) {
		std::copy(a[p[i]], a[p[i]] + n, pa[i]);
	}

	//make sure that L * U == P * A
	Mat<Real> recon = l * u;
	ReVal<Real> delta = pa.compare(recon) / n;
	if(delta > detail::eps<Real>() * 50) {//quality of reconstruction shouldn't depend on condition number
		os << "L * U != A for :\n" << a << "L:\n" << l << "U:\n" << u;
		os << "L * U: " << recon << '\n';
		os << "P * A: " << pa << '\n';
		os << "error: " << delta << '\n';
		os << "condition number: " << cn << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}

	//solve a random system
	std::vector<Real> b(n), x(n), ax(n);
	std::uniform_real_distribution<ReVal<Real>> dis(ReVal<Real>(-1), ReVal<Real>(1));
	for(Real& v : b) v = rndPhs<Real>(dis(gen), gen);//build random b
	backsolve::lu(lu.data(), p.data(), x.data(), b.data(), n);//compute x
	for(size_t i = 0; i < n; i++) ax[i] = std::inner_product(x.begin(), x.end(), a[i], Real(0));//compute a * x
	
	//make sure than a * x == b
	delta = vectorCompare(ax, b) / n;
	if(delta > detail::eps<Real>() * std::max<ReVal<Real>>(cn, 50)) {//quality of back solution will depend on condition number (but give ourselves some breathing room for very small cn)
		os << "A * x != b for LU:\n" << a << '\n';
		os << "\nx :"; for(auto v : x ) os << ' ' << v;
		os << "\nb :"; for(auto v : b ) os << ' ' << v;
		os << "\nAx:"; for(auto v : ax) os << ' ' << v;
		os << "\nerror: " << delta << '\n';
		os << "condition number: " << cn << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}
	return true;
}

//@brief    : test cholesky decomposition for a random matrix matrix
//@param os : location to write error messages
//@param n  : matrix size
//@param cn : condition number to test for
//@param gen: random generator e.g. std::mt19937
//@return   : true/false if test passed/failed
template<typename Real, class URNG>
bool testCholesky(std::ostream& os, const size_t n, const ReVal<Real> cn, URNG& gen) {
	//create a random matrix positive/negative definite matrix
	Mat<Real> a = randomDefinite<Real>(n, cn, gen, (gen() % 2) == 0);

	//do in place decomposition on a copy
	Mat<Real> llt(a);//copy to keep a
	std::vector<Real> d(n);//allocate diagonal matrix
	const bool neg = decompose::cholesky(llt.data(), d.data(), n);//in place decompose

	//extract L
	Mat<Real> l(llt);
	for(size_t r = 0; r < n; r++) {
		l[r][r] = d[r];
		for(size_t c = r+1; c < n; c++) {
			l[r][c] = Real(0);
		}
	}

	//make sure that L * L^T == a
	Mat<Real> recon = l * l.transpose();
	if(neg) for(Real& v : recon) v = -v;//handle negative definite matrices
	ReVal<Real> delta = a.compare(recon) / n;
	if(delta > detail::eps<Real>() * 20) {
		os << "L * L^T != A for :\n" << a << "L:\n" << l;
		os << "L * L^T: " << recon << '\n';
		os << "error: " << delta << '\n';
		os << "condition number: " << cn << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}

	//solve a random system
	std::vector<Real> b(n), x(n), ax(n);
	std::uniform_real_distribution<ReVal<Real>> dis(ReVal<Real>(-1), ReVal<Real>(1));
	for(Real& v : b) v = rndPhs<Real>(dis(gen), gen);//build random b
	backsolve::cholesky(llt.data(), d.data(), x.data(), b.data(), n, neg);//compute x
	for(size_t i = 0; i < n; i++) ax[i] = std::inner_product(x.begin(), x.end(), a[i], Real(0));//compute a * x

	//make sure than a * x == b
	delta = vectorCompare(ax, b) / n;
	if(delta > detail::eps<Real>() * std::max<ReVal<Real>>(cn, 10)) {//quality of back solution will depend on condition number (but give ourselves some breathing room for very small cn)
		os << "A * x != b for cholesky:\n" << a << '\n';
		os << "\nx :"; for(auto v : x ) os << ' ' << v;
		os << "\nb :"; for(auto v : b ) os << ' ' << v;
		os << "\nAx:"; for(auto v : ax) os << ' ' << v;
		os << "\nerror: " << delta << '\n';
		os << "condition number: " << cn << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}
	return true;
}

//@brief    : test QR decomposition for a random mxn matrix
//@param os : location to write error messages
//@param m  : number of rows
//@param n  : number of columns
//@param gen: random generator e.g. std::mt19937
//@return   : true/false if test passed/failed
template<typename Real, class URNG>
bool testQR(std::ostream& os, const size_t m, const size_t n, URNG& gen) {
	//create a random matrix
	Mat<Real> a(m,n);
	std::uniform_real_distribution<ReVal<Real>> dis(ReVal<Real>(-1), ReVal<Real>(1));
	for(Real& v : a) v = rndPhs<Real>(dis(gen), gen);

	//copy and do the in place QR decomposition
	Mat<Real> qr(a);
	decompose::qr(qr.data(), m, n);

	//copy decomposition and clear out lower half (build R)
	Mat<Real> r(qr);
	for(size_t i = 1; i < m; i++)
		for(size_t j = 0; j < std::min(i,n); j++)
			r[i][j] = Real(0); 

	//reconstruct Q matrix
	Mat<Real> q = Mat<Real>::Ident(m);
	qr::applyQ(qr.data(), q.data(), m, n, m);

	//reconstruct Q^H matrix
	Mat<Real> qh = Mat<Real>::Ident(m);
	qr::applyQH(qr.data(), qh.data(), m, n, m);

	//compute q*r with matrix multiplication and q multiplication function
	Mat<Real> aReconMult = q*r;
	Mat<Real> aReconQ(r);
	qr::applyQ(qr.data(), aReconQ.data(), m, n, n);

	//compute r with applyQH
	Mat<Real> rReconQ(a);
	qr::applyQH(qr.data(), rReconQ.data(), m, n, n);

	//make sure qh is actually Q^H
	ReVal<Real> delta;
	delta = q.compare(qh.conjugateTranspose()) / std::max(m, n);
	if(delta > detail::eps<Real>()) {
		os << "Q != Q^H^H != decomposition of\n" << a;
		os << "Q:\n" << q << "Q^H:\n" << qh;


		os << "\nerror: " << delta << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}

	//make sure Q * Q^H is identity (check that Q is orthogonal/unitary)
	Mat<Real> qqh = q * qh;
	delta = qqh.compare(Mat<Real>::Ident(m)) / std::max(m, n);
	if(delta > detail::eps<Real>() * 10) {
		os << "Q * Q^H != I for decomposition of\n" << a;
		os << qqh;
		os << "Q:\n" << q << "Q^H:\n" << qh;


		os << "\nerror: " << delta << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		return false;
	}

	//make sure that Q R == A
	delta = a.compare(aReconMult) / std::max(m, n);
	ReVal<Real> eps = std::sqrt(detail::eps<Real>());
	if(delta > eps) {
		os << "Q * R != A for decomposition of\n" << a;
		os << "Q:\n" << q << "R:\n" << r << "Q*R\n" << aReconMult;


		os << "\nerror: " << delta << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		os << detail::eps<Real>() << '\n';
		return false;
	}

	// make sure that QR computed through applyQ(R) == A
	delta = a.compare(aReconQ) / std::max(m, n);
	if(delta > eps) {
		os << "applyQ(R) != A for decomposition of\n" << a;
		os << "Q:\n" << q << "R:\n" << r << "Q*R\n" << aReconQ;


		os << "\nerror: " << delta << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		os << detail::eps<Real>() << '\n';
		return false;
	}

	//make sure that R computed through applyQH(A) == R
	delta = r.compare(rReconQ) / std::max(m, n);
	if(delta > eps) {
		os << "applyQH(A) != R for decomposition of\n" << a;
		os << "Q^H:\n" << qh << "R:\n" << r << "Q^H*A\n" << rReconQ;


		os << "\nerror: " << delta << '\n';
		os << "eps: " << detail::eps<Real>() << '\n';
		os << detail::eps<Real>() << '\n';
		return false;
	}
	return true;
}

//@brief : test real valued linalg functions
//@return: true / false if tests pass / fail
template <typename Real>
bool testLinAlg(std::ostream& os) {
	//seed random number generation
	std::mt19937 gen(0);//deterministic behavior

	//build a range of difficulties (matrix condition numbers)
	std::vector<int> k;
	std::vector<ReVal<Real>> cond;
	k.push_back(std::numeric_limits< ReVal<Real> >::digits10);//start from the worst possible round off error (in base 10 sigfigs), this is 6 for float, 15 for double
	while(k.back() != 0) k.push_back(std::max<int>(0, k.back() - 3));//build k in multiples of 3 (somewhat arbitrary)
	for(const int& i : k) cond.push_back(std::pow(ReVal<Real>(10), i));//now convert from ~lost digits to condition number (10^k)

	//test QR decomposition for a range of sizes and condition numbers
	os << "\ttesting QR" << std::endl;
	for(size_t m = 1; m < 32; m++) {
		for(size_t n = 1; n < 32; n++) {
			if(!testQR<Real>(os, m, n, gen)) return false;
		}
	}

	//if QR passes we can use it to generate random unitary matrices

	//test LU decomposition for a range of sizes and condition numbers
	os << "\ttesting LU" << std::endl;
	for(size_t n = 1; n < 64; n++) {
		for(const ReVal<Real>& c : cond) {
			if(!testLU<Real>(os, n, c, gen)) return false;
		}
	}

	//test cholesky decomposition for a range of sizes and condition numbers
	os << "\ttesting Cholesky" << std::endl;
	cond.pop_back();//cholesky is less stable than lu, don't do the worst case
	for(size_t n = 1; n < 64; n++) {
		for(const ReVal<Real>& c : cond) {
			if(!testCholesky<Real>(os, n, c, gen)) return false;
		}
	}

	//if we made it this far all tests passed
	return true;
}
