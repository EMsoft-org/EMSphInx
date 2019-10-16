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

#ifndef _SQUARE_SHT_H_
#define _SQUARE_SHT_H_

#include <complex>
#include <vector>
#include <memory>

#include "util/fft.hpp"//for ffts in DiscreteSHT

//@brief: functions related to square <--> hemisphere mappings
//mappings must be the following criteria:
//  -rings of constant latitude points on the sphere map to square rings
//  -rings are symmtric across the equator
//  -points within a ring a equally spaced
//  -the sidelength of the square grid is odd
//  -the number of points in the rings (from pole -> equator) is 1, 8, 16, 24, 32, ... 
//  -the first point in each ring is at an azimuthal angle of 0
//  -the exterior of each square is the same and lies on the equator

namespace emsphinx {

	namespace square {

		enum class Layout {//types of square grids
			Lambert ,//equal area projection
			Legendre //equal lattitude rings at legendre roots (+ extra points at poles) for improved numerical stability
		};

		//helper class to perform discrete spherical harmonic transformations with the square lambert grid
		//@reference: Reinecke , M. (2011). Libpsht–algorithms for efficient spherical harmonic transforms. Astronomy & Astrophysics, 526, A108.
		//@reference: Schaeffer, N. (2013). Efficient spherical harmonic transforms aimed at pseudospectral numerical simulations. Geochemistry, Geophysics, Geosystems, 14(3), 751-758.
		//@note: grid requirements:
		//         -pixels are arranged on N_\phi \approx \sqrt{N_{pix}} iso-latitude rings w/ colatitude of ring y = \theta_y
		//         -within a ring pixels are equidistance in azimuthal angle and have identical weight (solid angle) w_y
		//         -the number of pixels in a given ring N_{\phi,y} can vary
		//         -the first pixel in a ring is at azimuthal angle \phi_{0,y}
		//       the square lambert grid satisfies these requirements:
		//         -\sqrt{ 2.5 * N_{pix} } rings for side length 3, rapidly approaching \sqrt{2 * N_{pix} } rings (within 1% at side length 9)
		//         -all grid points cover the same solid angle -> w_y = 1 / rings
		//         -\phi_{0,y} = 0 for odd side lengths (must be calculated for even side lengths)
		//         -for the reverse transformation (synthesis) all N_{\phi,y} are even
		//@note: optimizations implemented (see Schaeffer for details):
		//         -use of real values FFT -> factor of 2 savings on FFT
		//         -mirror (conjugate) symmetry of spherical harmonics -> factor of 2 on direct summations [only possible for ring positions that are symmetric across the equator]
		//         -polar optimization (implicitely via equal area grid)
		//@note: I've opted for on the fly (vs precomputed calculations) since the single thread performance was very similar and on the fly has much lower memory overhead
		//@note: single thread performance is comparable to SHTns compiled without SIMD enabled for the same number of rings, but non legendre root rings -> half the bandwidth
		//@note: complexity is ~n^2.7 for reasonable bandwidths (should probably be something like log(n)*n^2)
		template <typename Real>
		class DiscreteSHT {
			struct Constants;//helper struct to hold read only constants
			const std::shared_ptr<const Constants> shtLut      ;//read only values (can be shared across threads)
			fft::vector< std::complex<Real> >      cWrk1, cWrk2;//complex working arrays for a single ring
			fft::vector<              Real  >      rWrk1, rWrk2;//real working array for a single ring

			public:
				//@brief    : construct a spherical harmonic transformer for a given side length and max bandwidth
				//@param dim: side length of square projection to perform transformations for
				//@param mBw: maximum bandwidth to compute
				//@param lay: type of square projection to use
				DiscreteSHT(const size_t dim, const size_t mBw, const Layout lay);

				//@brief    : construct a spherical harmonic transformer for a given side length and max bandwidth
				//@param dim: side length of square projection to perform transformations for
				//@param mBw: maximum bandwidth to compute
				//@param cLt: cosines or ring latitudes (northern hemisphere only, mirrored across equator)
				//@param leg: true/false if the ring latitudes in cLt are legendre roots (+poles for odd sizes)
				DiscreteSHT(const size_t dim, const size_t mBw, Real const * const cLt, const bool leg);

				//@brief    : convince method for constructing a square legendre transformer for the specified bandwidth
				//@param dim: side length of square legendre projection to perform transformations for
				//@return   : transformer for specified sidelength, available bandwidth can be determined with maxBw()
				static DiscreteSHT Legendre(const size_t dim) {return DiscreteSHT(dim, dim - 2, Layout::Legendre);}//poles aren't used

				//@brief    : convince method for constructing a square lambert transformer for the specified bandwidth
				//@param dim: side length of square lambert projection to perform transformations for
				//@return   : transformer for specified sidelength, available bandwidth can be determined with maxBw()
				static DiscreteSHT Lambert(const size_t dim) {return DiscreteSHT(dim, (dim-1)/2, Layout::Lambert);}//half number of rings, equator has double cover

				//@brief    : compute spherical harmonic coefficients from a spherical function (forward transformation)
				//@param nh : value of function to analyze at each grid point for the north hemisphere (row major order)
				//@param sh : value of function to analyze at each grid point for the south hemisphere (row major order)
				//@param alm: location to write alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param bw : maximum bandwidth to compute (must be <= mBw argument from construction, 0 to use mBw from construction)
				//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
				//@note     : a^l_{-m} = std::conj((m % 2 == 0) ? a^l_{-m} : -a^l_{-m}) if they are needed
				void analyze(Real const * const nh, Real const * const sh, std::complex<Real> * const alm, const size_t bw, const size_t stM);

				//@brief    : compute spherical harmonic coefficients from a spherical function (forward transformation)
				//@param nh : value of function to analyze at each grid point for the north hemisphere (row major order)
				//@param sh : value of function to analyze at each grid point for the south hemisphere (row major order)
				//@param alm: location to write alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@note     : a^l_{-m} = std::conj((m % 2 == 0) ? a^l_{-m} : -a^l_{-m}) if they are needed
				void analyze(Real const * const nh, Real const * const sh, std::complex<Real> * const alm) {analyze(nh, sh, alm, 0, 0);}

				//@brief    : compute spherical harmonic coefficients from a spherical function (forward transformation)
				//@param pts: value of function to analyze at each grid point (row major order, northern followed by southern hemisphere)
				//@param alm: location to write alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param bw : maximum bandwidth to compute (must be <= mBw argument from construction, 0 to use mBw from construction)
				//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
				//@note     : a^l_{-m} = std::conj((m % 2 == 0) ? a^l_{-m} : -a^l_{-m}) if they are needed
				void analyze(Real const * const pts, std::complex<Real> * const alm, const size_t bw, const size_t stM) {analyze(pts, pts + shtLut->dim * shtLut->dim, alm, bw, stM);}

				//@brief    : compute spherical harmonic coefficients from a spherical function (forward transformation)
				//@param pts: value of function to analyze at each grid point (row major order, northern followed by southern hemisphere)
				//@param alm: location to write alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@note     : a^l_{-m} = std::conj((m % 2 == 0) ? a^l_{-m} : -a^l_{-m}) if they are needed
				void analyze(Real const * const pts, std::complex<Real> * const alm) {analyze(pts, alm, 0, 0);}

				//@brief    : compute spherical function from spherical harmonic coefficients (inverse transformation)
				//@param alm: alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param nh : location to write north hemisphere of spherical function (row major order)
				//@param sh : location to write north hemisphere of spherical function (row major order)
				//@param bw : maximum bandwidth to use in synthesis (must be <= mBw argument from construction, 0 to use mBw from construction)
				//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
				//@note     : only non-negative m values are used since SHT of real function is conjugate symmetric
				void synthesize(std::complex<Real> const * const alm, Real * const nh, Real * const sh, const size_t bw, const size_t stM);

				//@brief    : compute spherical function from spherical harmonic coefficients (inverse transformation)
				//@param alm: alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param nh : location to write north hemisphere of spherical function (row major order)
				//@param sh : location to write north hemisphere of spherical function (row major order)
				//@note     : only non-negative m values are used since SHT of real function is conjugate symmetric
				void synthesize(std::complex<Real> const * const alm, Real * const nh, Real * const sh) {synthesize(alm, nh, sh, 0, 0);}

				//@brief    : compute spherical function from spherical harmonic coefficients (inverse transformation)
				//@param alm: alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param pts: location to write spherical function (row major order, northern followed by southern hemisphere)
				//@param bw : maximum bandwidth to use in synthesis (must be <= mBw argument from construction, 0 to use mBw from construction)
				//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
				//@note     : only non-negative m values are used since SHT of real function is conjugate symmetric
				void synthesize(std::complex<Real> const * const alm, Real * const pts, const size_t bw, const size_t stM) {synthesize(alm, pts, pts + shtLut->dim * shtLut->dim, bw, stM);}

				//@brief    : compute spherical function from spherical harmonic coefficients (inverse transformation)
				//@param alm: alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
				//@param pts: location to write spherical function (row major order, northern followed by southern hemisphere)
				//@note     : only non-negative m values are used since SHT of real function is conjugate symmetric
				void synthesize(std::complex<Real> const * const alm, Real * const pts) {synthesize(alm, pts, 0, 0);}

				//@brief : get the maximum bandwidth that can be computed with analyze / used by synthesize
				//@return: maximum bandwidth
				size_t maxBw() const;

				//@brief : get the side length of input signals for analyze / output signals from synthesize
				//@return: side length of square projection
				size_t dim() const;
		};

		//these functions are specifically for the area preserving square <--> sphere mapping
		//@reference: Roşca, D. (2010). New uniform grids on the sphere. Astronomy & Astrophysics, 520, A63.
		namespace lambert {
			//@brief: square lambert projection from unit hemisphere to unit square
			//@param x: x coordinate on unit sphere
			//@param y: y coordinate on unit sphere
			//@param z: z coordinate on unit sphere
			//@param X: location to write x coordinate in unit square (0,1)
			//@param Y: location to write y coordinate in unit square (0,1)
			template <typename Real> void sphereToSquare(Real const& x, Real const& y, Real const& z, Real& X, Real& Y);

			//@brief: square lambert projection from unit square to unit hemisphere
			//@param X: x coordinate in unit square (0,1)
			//@param Y: y coordinate in unit square (0,1)
			//@param x: location to write x coordinate on unit sphere
			//@param y: location to write y coordinate on unit sphere
			//@param z: location to write z coordinate on unit sphere
			template <typename Real> void squareToSphere(Real const& X, Real const& Y, Real& x, Real& y, Real& z);

			//@brief    : compute cosine of the latitude of each ring in the northern hemisphere (including equator)
			//@note     : southern hemisphere cosines can be computed by symmetry with cos(lat[i]) = -cos(lat[# rings - i])
			//@param dim: side length of square lambert projection
			//@param lat: location to write cos(ring latitudes)
			template <typename Real> void cosLats(const size_t dim, Real * const lat);

			//@brief    : compute the real space coordinates of unprojected points in a square lambert grid (north hemisphere only)
			//@param dim: side length of square legendre grid (must be odd)
			//@param xyz: location to write real space coordiantes (dim * dim * 3)
			template <typename Real> void normals(const size_t dim, Real * const xyz);

			//@brief    : compute the solid angle correction of each grid point in the square lambert projection
			//@param dim: side length of square lambert projection to compute solid angles for
			//@param omg: location to write grid point solid angles (northern hemisphere only, row major order)
			//@reference: Mazonka, Oleg. "Solid angle of conical surfaces, polyhedral cones, and intersecting spherical caps." arXiv preprint arXiv:1205.1396 (2012).
			template <typename Real> void solidAngles(const size_t dim, Real * const omg);
		}

		//these functions are specifically for the grid with rings at legendre polynomial roots
		namespace legendre {
			//@brief    : compute cosine of the latitude of each ring in the northern hemisphere (including equator)
			//@note     : southern hemisphere cosines can be computed by symmetry with cos(lat[i]) = -cos(lat[# rings - i])
			//@param dim: side length of square legendre projection
			//@param lat: location to write cos(ring latitudes)
			//@reference: Barth, W., Martin, R. S., & Wilkinson, J. H. (1967). Calculation of the eigenvalues of a symmetric tridiagonal matrix by the method of bisection. Numerische Mathematik, 9(5), 386-393.
			//@method   : legendre roots are calculated as zeros of a symmetric tridiagonal matrix as described here - https://math.stackexchange.com/questions/12160/roots-of-legendre-polynomial/12209#12209
			template <typename Real> void roots(const size_t dim, Real * const lat);

			//@brief    : compute the real space coordinates of unprojected points in a square legendre grid (north hemisphere only)
			//@param dim: side length of square legendre grid (must be odd)
			//@param xyz: location to write real space coordiantes (dim * dim * 3)
			template <typename Real> void normals(const size_t dim, Real * const xyz);

			//@brief     : compute the 4 bounding indices of for a given unit direction in a the north hemisphere of square legendre grid
			//@param dim : side length of square legendre grid (must be odd)
			//@param zLat: z coordinates of each ring in the projection (e.g. from cosLats)
			//@param n   : direction to get bounding indices for
			//@param inds: location to write bounding indices
			template <typename Real> void boundingInds(const size_t dim, Real const * const zLat, Real const * const n, size_t * const inds);
		}

		//@brief     : extract a single ring from a row major pattern
		//@param dim : side length of square in pixels
		//@param ring: ring number to copy (0 for north pole, (dim-1)/2 for equator)
		//@param ptr : pointer to start of hemisphere to copy from (dim * dim array in row major order)
		//@param buff: location to write extracted ring
		//@return    : number of values copied
		template <typename T> size_t readRing(const size_t dim, const size_t ring, T const * const ptr, T * const buff);

		//@brief     : write a single ring into a row major pattern
		//@param dim : side length of square in pixels
		//@param ring: ring number to write (0 for north pole, (dim-1)/2 for equator)
		//@param ptr : pointer to start of hemisphere to write to (dim * dim array in row major order)
		//@param buff: location to read ring from
		//@return    : number of values copied
		template <typename T> size_t writeRing(const size_t dim, const size_t ring, T * const ptr, T const * const buff);

		//@brief    : compute quadrature weights for rings (w_y in equation 10 of Reinecke)
		//@param dim: side length of square lambert projection to compute weights for
		//@param lat: cosines of ring latitudes (symmetric across equator)
		//@param wgt: location to write weights for each row
		//@param skp: ring to exclude from weights (e.g. skip = 0 will exclude the poles)
		//@reference: https://doi.org/10.1111/j.1365-246X.1994.tb03995.x
		template <typename Real> void computeWeightsSkip(const size_t dim, Real const * const lat, Real * const wgt, const size_t skp);

		//@brief     : compute the cosines of ring latitudes for a given square grid
		//@param dim : side length of square projection to compute ring latitudes for
		//@param type: type of square projection to compute latitudes for
		//@return    : cosine(latitudes) from north pole -> equator
		template <typename Real> typename std::vector<Real> cosLats(const size_t dim, const Layout type);

		//@brief     : compute the normals of the northern hemisphere for a given square grid
		//@param dim : side length of square projection to compute ring latitudes for
		//@param type: type of square projection to compute latitudes for
		//@return    : {x,y,z} normals for north pole
		template <typename Real> typename std::vector<Real> normals(const size_t dim, const Layout type);

		//@brief     : compute the solid angle correction of each ring for a given square grid
		//@param dim : side length of square lambert projection to compute solid angles for
		//@param type: type of square projection to compute latitudes for
		//@return    : solide angle correction from north pole -> equator (actual pixel size / average pixel size)
		template <typename Real> typename std::vector<Real> solidAngles(const size_t dim, const Layout type);

		//@brief    : convert from a vecotrized index to a ring number
		//@param dim: side length of square lambert projection to compute index for
		//@param idx: vectorized index (y * dim + x) to compute ring number for [must be in north hemisphere i.e. < dim * dim]
		//@return   : ring number (e.g. for indexing into cosLats)
		size_t ringNum(const size_t dim, const size_t idx);
	}

}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>
#include <numeric>

#include "constants.hpp"//for computeWeights
#include "util/linalg.hpp"//for computeWeights

namespace emsphinx {
	
	namespace square {
		//helper struct to hold read only constants needed by DiscreteSHT (for sharing across threads)
		template <typename Real>
		struct DiscreteSHT<Real>::Constants {
			//these members could be shared across threads since they should be read only after construction
			const size_t        dim  ;//side length of square lambert projection
			const size_t        Nt   ;//number of pairs of equal latitude rings [(dim+1) / 2]
			const size_t        maxL ;//maximum bandwidth of square lambert projection (must be < Nt for arbitrary rings (Nt*2 for legendre rings))
			const size_t        Nw   ;//number of different types of weights [(dim-2) / 8 + 1]
			std::vector<Real  > wy   ;//weighting factor for each ring [Nt * Nw]
			std::vector<Real  > cosTy;//cosine of latitude of each ring [Nt]
			std::vector<Real  > amn  ;//precomputed a^m_n values for on the fly ylm calculation [maxL^2]
			std::vector<Real  > bmn  ;//precomputed b^m_n values for on the fly ylm calculation [maxL^2]
			std::vector< std::shared_ptr< fft::RealFFT<Real> > > ffts;//fft plans [Nt]

			//@brief  : construct a set of constants for a given side length and max bandwidth
			//@param d: side length of square lambert projection to perform transformations for
			//@param l: maximum bandwidth to compute
			//@param c: cosines or ring latitudes (northern hemisphere only, mirrored across equator)
			//@param b: true/false if the ring latitudes in c are legendre roots (+poles for odd sizes)
			Constants(const size_t d, const size_t l, Real const * const c, const bool b);
		};

		//@brief  : construct a set of constants for a given side length and max bandwidth
		//@param d: side length of square lambert projection to perform transformations for
		//@param l: maximum bandwidth to compute (exclusive)
		//@param c: cosines or ring latitudes (northern hemisphere only, mirrored across equator)
		//@param b: true/false if the ring latitudes in c are legendre roots (+poles for odd sizes)
		template <typename Real>
		DiscreteSHT<Real>::Constants::Constants(const size_t d, const size_t l, Real const * const c, const bool b) : dim(d), Nt( (dim+1) / 2 ), maxL(l), Nw((dim-2) / 4 + 1), wy(Nt * Nw), cosTy(c, c + Nt), amn(maxL*maxL), bmn(maxL*maxL), ffts(Nt) {
			//sanity check bandwidth and dimensions
			if(dim < 3) throw std::domain_error("square lambert side length must be at least 3");
			if(dim % 2 == 0) throw std::domain_error("only odd side lengths are supported");
			if(b) {
				const size_t limit = 2 * (Nt - 1) + (1 - dim % 2);
				if(maxL >= limit) throw std::domain_error("maximum bandwidth is side # rings - 1 (for legendre root latitudes)");
			} else {
				if(maxL >= Nt) throw std::domain_error("maximum bandwidth is side length / 2");
			}

			//compute amn and bmn values
			const Real k4p = Real(1) / (emsphinx::Constants<Real>::pi * 4);//1 / (4 pi): constant for a^m_m calculation
			Real kamm = 1;//\Pi_k=1^|m| \frac{2k+1}{2k} for m = 0: for a^m_m calculation
			for(size_t m = 0; m < maxL; m++) {
				//first compute a^m_m
				amn[m*maxL + m] = std::sqrt(kamm * k4p);//recursively compute a^m_m (Schaeffer equation 16)
				kamm *= Real(2*m+3) / (2*m+2);//update recursion for \Pi_k=1^|m| \frac{2k+1}{2k}
				if(m+1 == maxL) break;

				//now compute a^m_{m+1}
				const size_t m2 = m * m;//we'll need this value a bunch of times
				size_t n12 = m2;//(n-1)^2 for n = m+1
				size_t n2 = (m+1) * (m+1);//n^2 for n = m+1
				size_t n2m2 = n2 - m2;//n^2 - m^2 for n = m+1
				Real n12m2 = 0;//(n-1)^2 - m^2 for n = m+1 (not actually needed here but will be for recursion)
				amn[m*maxL + m+1] = std::sqrt( Real(4 * n2 - 1) / n2m2 );//a^m_n for n = m+1 (Schaeffer equation 17)
				
				//now compute remaining a^m_n and b^m_n values
				for(size_t n = m+2; n < maxL; n++) {//l is more commonly used but I'll stick with paper's notation
					n12 = n2;//use previously computed value of n^2 for (n-1)^2
					n12m2 = (Real)n2m2;//use previously computed values of n^2 - m^2 for (n-1)^2 - m^2
					n2 = n * n;//compute n^2
					n2m2 = n2 - m2;//compute n^2 - m^2
					amn[m*maxL + n] = std::sqrt( Real(4 * n2 - 1) / n2m2 );//a^m_n (Schaeffer equation 17)
					bmn[m*maxL + n] = std::sqrt( ( Real(2*n+1) / (2*n-3) ) * ( n12m2 / n2m2 ) );//b^m_n (Schaeffer equation 18)
				}
			}

			//compute ring weights (this is the most expensive step for larger sizes since the scalling is ~dim^4)
			if(b) {
				computeWeightsSkip(dim, cosTy.data(), wy.data(), 0);//skip poles
				for(size_t i = 1; i < Nw; i++) std::copy(wy.begin(), wy.begin() + Nt, wy.begin() + Nt * i);
			} else {
				for(size_t i = 0; i < Nw; i++) computeWeightsSkip(dim, cosTy.data(), wy.data() + i * Nt, i);
			}

			//build fft calculators
			for(size_t y = 0; y < Nt; y++) ffts[y] = std::make_shared< fft::RealFFT<Real> >(std::max<size_t>(1, 8 * y), fft::flag::Plan::Patient);//we'll be doing many transforms, take the time to find a fast path
		}

		//@brief    : construct a spherical harmonic transformer for a given side length and max bandwidth
		//@param dim: side length of square lambert projection to perform transformations for
		//@param mBw: maximum bandwidth to compute
		//@param lay: type of square projection to use
		template <typename Real>
		DiscreteSHT<Real>::DiscreteSHT(const size_t dim, const size_t mBw, const Layout lay) : DiscreteSHT(dim, mBw, cosLats<Real>(dim, lay).data(), lay == Layout::Legendre) {}

		//@brief    : construct a spherical harmonic transformer for a given side length and max bandwidth
		//@param dim: side length of square lambert projection to perform transformations for
		//@param mBw: maximum bandwidth to compute (exclusive)
		//@param cLt: cosines or ring latitudes (northern hemisphere only, mirrored across equator)
		//@param leg: true/false if the ring latitudes in cLt are legendre roots (+poles for odd sizes)
		template <typename Real>
		DiscreteSHT<Real>::DiscreteSHT(const size_t dim, const size_t mBw, Real const * const cLt, const bool leg) :
			shtLut(std::make_shared<const Constants>(dim, mBw, cLt, leg)),
			cWrk1 (std::max<size_t>(1, 4 * (dim-1))),
			cWrk2 (std::max<size_t>(1, 4 * (dim-1))),
			rWrk1 (std::max<size_t>(1, 4 * (dim-1))),
			rWrk2 (std::max<size_t>(1, 4 * (dim-1))) {}

		//@brief    : compute spherical harmonic coefficients from a spherical function (forward transformation)
		//@param nh : value of function to analyze at each grid point for the north hemisphere (row major order)
		//@param sh : value of function to analyze at each grid point for the south hemisphere (row major order)
		//@param alm: location to write alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
		//@param bw : maximum bandwidth to compute (must be <= mBw argument from construction, 0 to use mBw from construction)
		//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
		//@note     : a^l_{-m} = std::conj((m % 2 == 0) ? a^l_{-m} : -a^l_{-m}) if they are needed
		template <typename Real>
		void DiscreteSHT<Real>::analyze(Real const * const nh, Real const * const sh, std::complex<Real> * const alm, const size_t bw, const size_t stM) {
			//parse arguments, sanity check, and initialize output
			const size_t maxL = bw == 0 ? shtLut->maxL : bw;//what is the maximum bandwidth we need to compute
			if(maxL > shtLut->maxL) throw std::runtime_error("maximum bandwidth must be <= maximum from DiscreteSHT construction");
			const size_t stride = stM == 0 ? maxL : stM;//what is the stride of the output array
			if(stride < maxL) throw std::runtime_error("output array isn't big enough to store all coefficients");
			std::fill(alm, alm + stride * maxL, std::complex<Real>(0));//fill SHT with 0

			//compute SHT ring at a time
			for(size_t y = 0; y < shtLut->Nt; y++) {//loop over rings
				//copy current rings and compute G_{m,y} by fft (Reinecke equation 10) leveraging real symmetry
				const size_t Npy  = std::max<size_t>(1, 8 * y);//get number of points in this ring
				const size_t fftN = Npy / 2 + 1;//number of fft points: half from real -> conjugate symmetric, this includes problematic points (imag == 0 from conmjugate symmetry)
				readRing(shtLut->dim, y, nh, rWrk1.data());//copy current ring from row major order (north hemisphere)
				readRing(shtLut->dim, y, sh, rWrk2.data());//copy current ring from row major order (south hemisphere)
				shtLut->ffts[y]->forward(rWrk1.data(), cWrk1.data());//do ffts into complex working arrays (north hemisphere)
				shtLut->ffts[y]->forward(rWrk2.data(), cWrk2.data());//do ffts into complex working arrays (south hemisphere)

				//compute G_{m,y} +/- G_{m,Nt-1-y} since they are multiplied by symmetric values
				const size_t mLim = std::min<size_t>(maxL, fftN);//anything after l+1 isn't needed and anything after fftN is 0
				for(size_t m = 0; m < mLim; m++) {
					//get ring weight w_y
					//weights excluding only the closest problematic ring (missing complex value due to real even dft) are most stable + accurate
					//negate odd m values to correct for sign error in legendre polynomial calculation
					const Real wy = shtLut->wy[size_t(m/4) * shtLut->Nt + y] * (m % 2 == 1 ? -1 : 1);//mod 4 from rings having 8y points + real symmetry of dft

					//combine northern / southern hemisphere rings to take advantages of spherical harmonic symmetry
					const std::complex<Real> nPt = cWrk1[m] * wy;//northern hemisphere point
					const std::complex<Real> sPt = cWrk2[m] * wy;//southern hemisphere point
					cWrk1[m] = (nPt + sPt) * Real(0.5);//for even l+m (harmonics are     symmetric across equator)
					cWrk2[m] = (nPt - sPt) * Real(0.5);//for odd  l+m (harmonics are antisymmetric across equator)
				}

				//calculate seed values for on the fly legendre polynomial calculation
				const Real x = shtLut->cosTy[y];//cosine of ring latitude
				const Real r1x2 = std::sqrt(Real(1) - x * x);// (1-x)^\frac{1}{2}: constant for P^m_m calculation
				Real kpmm = 1;//(1 - x^2) ^ \frac{|m|}{2} for m = 0: for P^m_m calculation

				//accumulate a_{l,m} via direct summation (Reinecke equation 9) via on the fly summation
				Real const * amn = shtLut->amn.data();
				Real const * bmn = shtLut->bmn.data();
				std::complex<Real> * plm = alm;
				for(int m = 0; m < mLim; m++) {
					//get weighted values from fft
					const std::complex<Real>& gmyS = cWrk1[m];//for     symmetric modes
					const std::complex<Real>& gmyA = cWrk2[m];//for antisymmetric modes

					//first compute P^m_m
					Real pmn2 = amn[m] * kpmm;//recursively compute P^m_m (Schaeffer equation 13)
					kpmm *= r1x2;//update recursion for (1 - x^2) ^ \frac{|m|}{2}
					plm[m] += gmyS * pmn2;
					if(m+1 == maxL) break;//P^m_{m+1} doesn't exist

					//now compute P^m_{m+1}
					Real pmn1 = amn[m+1] * x * pmn2;//P^m_n for n = m+1 (Schaeffer equation 14)
					plm[m+1] += gmyA * pmn1;
					
					//now compute remaining P^m_n values
					for(int n = m+2; n < maxL; n++) {//l is more commonly used but I'll stick with paper's notation
						Real pmn = amn[n] * x * pmn1 - bmn[n] * pmn2;//P^m_n (Schaeffer equation 15)
						pmn2 = pmn1;//push back pmn values in recursion
						pmn1 = pmn; //push back pmn values in recursion
						plm[n] += ((n+m) % 2 == 0 ? gmyS : gmyA) * pmn;
					}

					//increment pointers to next m
					amn += shtLut->maxL;
					bmn += shtLut->maxL;
					plm += stride;
				}
			}
		}

		//@brief    : compute spherical function from spherical harmonic coefficients (inverse transformation)
		//@param alm: alm values bw * bw with (m,l) stored: (0,0), (0,1), (0,2), ..., (0,bw-1), (1,0), (1,1), (1,2), ..., (bw-1,0), (bw-1,1), (bw-1,bw-1)
		//@param nh : location to write north hemisphere of spherical function (row major order)
		//@param sh : location to write north hemisphere of spherical function (row major order)
		//@param bw : maximum bandwidth to use in synthesis (must be <= mBw argument from construction, 0 to use mBw from construction)
		//@param stM: stride between sequential m values of the same l in alm, i.e. a^l_m is at alm[stM * m + l], 0 to use bw
		//@note     : only non-negative m values are used since SHT of real function is conjugate symmetric
		template <typename Real>
		void DiscreteSHT<Real>::synthesize(std::complex<Real> const * const alm, Real * const nh, Real * const sh, const size_t bw, const size_t stM) {
			//parse arguments and sanity check
			const size_t maxL = bw == 0 ? shtLut->maxL : bw;//what is the maximum bandwidth we need to compute
			if(maxL > shtLut->maxL) throw std::runtime_error("maximum bandwidth must be <= maximum from DiscreteSHT construction");
			const size_t stride = stM == 0 ? maxL : stM;//what is the stride of the output array
			if(stride < maxL) throw std::runtime_error("input array isn't big enough to store all coefficients");

			//compute inverse SHT ring at a time
			for(size_t y = 0; y < shtLut->Nt; y++) {//loop over rings
				//compute F_{m,y} +- F_{m,Nt-1-y} leveraging symmetry of spherical harmonics across equator
				const size_t Npy  = std::max<size_t>(1, 8 * y);//get number of points in this ring
				const size_t fftN = Npy / 2 + 1;//number of fft points: half from real -> conjugate symmetric, this includes problematic points (imag == 0 from conmjugate symmetry)

				//calculate seed values for on the fly legendre polynomial calculation
				const Real x = shtLut->cosTy[y];//cosine of ring latitude
				const Real r1x2 = std::sqrt(Real(1) - x * x);// (1-x)^\frac{1}{2}: constant for P^m_m calculation
				Real kpmm = 1;//(1 - x^2) ^ \frac{|m|}{2} for m = 0: for P^m_m calculation

				//accumulate F_{m,y} +- F_{m,Nt-1-y} by direct summation (Reinecke equation 8)
				Real const * amn = shtLut->amn.data();
				Real const * bmn = shtLut->bmn.data();
				std::complex<Real> const * plm = alm;

				const size_t mLim = std::min<size_t>(maxL, fftN);//anything after l+1 isn't needed and anything after fftN is 0
				for(size_t m = 0; m < mLim; m++) {
					//get reference to even and odd sums and zero
					std::complex<Real>& fmyS = cWrk1[m];//for     symmetric modes
					std::complex<Real>& fmyA = cWrk2[m];//for antisymmetric modes
					fmyS = fmyA = Real(0);

					//first compute P^m_m
					Real pmn2 = amn[m] * kpmm;//recursively compute P^m_m (Schaeffer equation 13)
					kpmm *= r1x2;//update recursion for (1 - x^2) ^ \frac{|m|}{2}
					fmyS += plm[m] * pmn2;
					if(m+1 == maxL) break;//P^m_{m+1} doesn't exist

					//now compute P^m_{m+1}
					Real pmn1 = amn[m+1] * x * pmn2;//P^m_n for n = m+1 (Schaeffer equation 14)
					fmyA += plm[m+1] * pmn1;

					//now compute remaining P^m_n values
					for(size_t n = m+2; n < maxL; n++) {//l is more commonly used but I'll stick with paper's notation
						Real pmn = amn[n] * x * pmn1 - bmn[n] * pmn2;//P^m_n (Schaeffer equation 15)
						pmn2 = pmn1;//push back pmn values in recursion
						pmn1 = pmn; //push back pmn values in recursion
						((n+m) % 2 == 0 ? fmyS : fmyA) += plm[n] * pmn;
					}

					//increment pointers to next m
					amn += shtLut->maxL;
					bmn += shtLut->maxL;
					plm += stride;
				}

				//now convert from F_{m,y} +- F_{m,Nt-1-y} -> F_{m,y} and F_{m,Nt-1-y}
				for(size_t m = 0; m < mLim; m++) {
					//combine northern / southern hemisphere rings to take advantages of spherical harmonic symmetry
					//negate odd m values to correct for sign error in legendre polynomial calculation
					const std::complex<Real> sigma = cWrk1[m] * Real(m % 2 == 1 ? -1 : 1);//F_{m,y} + F_{m,Nt-1-y}
					const std::complex<Real> delta = cWrk2[m] * Real(m % 2 == 1 ? -1 : 1);//F_{m,y} - F_{m,Nt-1-y}
					cWrk1[m] = sigma + delta;//northern hemisphere point
					cWrk2[m] = sigma - delta;//southern hemisphere point
				}

				//fill in remaining fft with 0 (we don't have data for frequencies this high)
				if(fftN >= maxL) {
					std::fill(cWrk1.begin() + maxL, cWrk1.begin() + fftN + 1, std::complex<Real>(0));
					std::fill(cWrk2.begin() + maxL, cWrk2.begin() + fftN + 1, std::complex<Real>(0));
				}

				//do the inverse ffts of F_{m,y} and F_{m,Nt-1-y} (Reinecke equation 7) and copy to output
				shtLut->ffts[y]->inverse(cWrk1.data(), rWrk1.data());//do ffts into real working arrays (north hemisphere)
				shtLut->ffts[y]->inverse(cWrk2.data(), rWrk2.data());//do ffts into real working arrays (south hemisphere)
				writeRing(shtLut->dim, y, nh, rWrk1.data());//copy current ring to row major order (north hemisphere)
				writeRing(shtLut->dim, y, sh, rWrk2.data());//copy current ring to row major order (south hemisphere)
			}
		}

		//@brief : get the maximum bandwidth that can be computed with analyze / used by synthesize
		//@return: maximum bandwidth
		template <typename Real>
		size_t DiscreteSHT<Real>::maxBw() const {return shtLut->maxL;}

		//@brief : get the side length of input signals for analyze / output signals from synthesize
		//@return: side length of square
		template <typename Real>
		size_t DiscreteSHT<Real>::dim() const {return shtLut->dim;}

		namespace lambert {
			//@brief: square lambert projection from unit hemisphere to unit square
			//@param x: x coordinate on unit sphere
			//@param y: y coordinate on unit sphere
			//@param z: z coordinate on unit sphere
			//@param X: x coordinate in unit square (0,1)
			//@param Y: y coordinate in unit square (0,1)
			template <typename Real>
			void sphereToSquare(Real const& x, Real const& y, Real const& z, Real& X, Real& Y) {
				static const Real kPi_4 = Real(0.7853981633974483096156608458199);//pi/4
				const Real fZ = std::fabs(z);
				if(fZ == 1.0) {
					X = Y = Real(0.5);
				} else if(std::abs(y) <= std::abs(x)) {
					X = std::copysign(std::sqrt(Real(1) - fZ), x) * Real(0.5);//[-0.5, 0.5]
					Y = X * std::atan(y / x) / kPi_4 + Real(0.5);//[0, 1]
					X += Real(0.5);//[0, 1]
				} else {
					Y = std::copysign(std::sqrt(Real(1) - fZ), y) * Real(0.5);//[-0.5, 0.5]
					X = Y * std::atan(x / y) / kPi_4 + Real(0.5);//[0, 1]
					Y += Real(0.5);//[0, 1]
				}
			}

			//@brief: square lambert projection from unit square to unit hemisphere
			//@param X: x coordinate in unit square (0,1)
			//@param Y: y coordinate in unit square (0,1)
			//@param x: location to write x coordinate on unit sphere
			//@param y: location to write y coordinate on unit sphere
			//@param z: location to write z coordinate on unit sphere
			template <typename Real>
			void squareToSphere(Real const& X, Real const& Y, Real& x, Real& y, Real& z) {
				static const Real kPi_4 = Real(0.7853981633974483096156608458199);//pi/4
				const Real sX = Real(2) * X - 1;//[0,1] -> [-1, 1]
				const Real sY = Real(2) * Y - 1;//[0,1] -> [-1, 1]
				const Real aX = std::abs(sX);
				const Real aY = std::abs(sY);
				const Real vMax = std::max<Real>(aX, aY);
				if(vMax <= std::numeric_limits<Real>::epsilon()) {
					x = y = Real(0);
					z = Real(1);
				} else {
					if(vMax > Real(1) + std::numeric_limits<Real>::epsilon()) throw std::runtime_error("point doesn't lie in square (0,0) -> (1,1)");
					if(aX <= aY) {
						const Real q  = sY * std::sqrt(Real(2) - sY * sY);
						const Real qq = kPi_4 * sX / sY;
						x = q * std::sin(qq);
						y = q * std::cos(qq);
					} else {
						const Real q = sX * std::sqrt(Real(2) - sX * sX);
						const Real qq = kPi_4 * sY / sX;
						x = q * std::cos(qq);
						y = q * std::sin(qq);
					}
					z = Real(1) - vMax * vMax;
					const Real mag = std::sqrt(x*x + y*y + z*z);
					x /= mag; y /= mag; z /= mag;
				}
			}

			//@brief    : compute cosine of the latitude of each ring in the northern hemisphere (including equator)
			//@note     : southern hemisphere cosines can be computed by symmetry with cos(lat[i]) = -cos(lat[# rings - i])
			//@param dim: side length of square lambert projection
			//@param lat: location to write cos(ring latitudes)
			template <typename Real>
			void cosLats(const size_t dim, Real * const lat) {
				const size_t count = (dim + 1) / 2;//only need north hemisphere rings
				const bool even = (0 == dim % 2);//even and odd are slightly different
				const size_t denom = (dim - 1) * (dim - 1);//compute denominator of latitude cosines
				size_t numer = denom - (even ? 1 : 0);//odd starts @ pole, even just off
				size_t delta = even ? 8 : 4;//starting increment is similarly different
				for(size_t i = 0; i < count; i++) {//loop over rings computing cosine of ring latitude (this is the same for even and odd)
					lat[i] = Real(numer) / denom;//compute cos(latitude(ring i))
					numer -= delta;//update numerator
					delta += 8;//update change to numerator
				}
			}

			//@brief    : compute the real space coordinates of unprojected points in a square lambert grid (north hemisphere only)
			//@param dim: side length of square legendre grid (must be odd)
			//@param xyz: location to write real space coordiantes (dim * dim * 3)
			template <typename Real> void normals(const size_t dim, Real * const xyz) {
				//loop over northern hemisphere compute
				for(size_t j = 0; j < dim; j++) {
					const Real y = Real(j) / (dim - 1);//[0, 1]
					for(size_t i = 0; i < dim; i++) {
						const Real x = Real(i) / (dim - 1);//[0, 1]
						Real * const n = xyz + 3 * (dim * j + i);
						squareToSphere(x, y, n[0], n[1], n[2]);
					}
				}
			}

			//@brief    : compute the solid angle correction of each grid point in the square lambert projection
			//@param dim: side length of square lambert projection to compute solid angles for
			//@param omg: location to write grid point solid angles (northern hemisphere only, row major order)
			//@reference: Mazonka, Oleg. "Solid angle of conical surfaces, polyhedral cones, and intersecting spherical caps." arXiv preprint arXiv:1205.1396 (2012).
			template <typename Real>
			void solidAngles(const size_t dim, Real * const omg) {
				const bool even = 0 == dim % 2;//check if the grid has an odd or even side length
				const size_t totalPixels = dim * dim * 2 - 4 * (dim - 1);//count total number of pixels in square lambert sphere
				const Real invOmegaBar = Real(totalPixels) / (Real(M_PI) * 4);//determine the reciprocal of the average solid angle of pixels in the projection
				const size_t mid = dim / 2;//first point on / across from the poles
				const Real delta = Real(0.5) / (dim - 1);//spacing from center to corenrs of solid angle in square lambert space
				for(size_t y = mid; y < dim; y++) {//loop over rows skipping first half (mirror plane at Y = 0)
					const bool maxY = dim == y + 1;//check if we're on the equator
					const Real Y = Real(y) / (dim - 1);//fractional position in unit square
					const Real Ym = Y - delta;//compute bottom extent of pixel
					const Real Yp = Y + (maxY ? 0 : delta);//compute top extent of pixel (without crossing equator)
					for(size_t x = y; x < dim; x++) {//loop over columns skipping first half (mirror plane X = 0) and stopping at 45 degrees (mirror plane at X == Y)
						const bool maxX = dim == x + 1;//check if we're on the equator
						const Real X = Real(x) / (dim - 1);//fractional position in unit square
						const Real Xm = X - delta;//compute left extent of pixel
						const Real Xp = X + (maxX ? 0 : delta);//compute right extent of pixel (without crossing equator)

						//now compute the spherical coordinates of the corners of the pixels
						Real s[4][3];//corners of square lambert projection (clockwise order)
						squareToSphere(Xm, Ym, s[0][0], s[0][1], s[0][2]);//bottom left
						squareToSphere(Xp, Ym, s[1][0], s[1][1], s[1][2]);//bottom right
						squareToSphere(Xp, Yp, s[2][0], s[2][1], s[2][2]);//top right
						squareToSphere(Xm, Yp, s[3][0], s[3][1], s[3][2]);//top left

						//compute solid angle of this pixel (see @reference equation 25 for details)
						std::complex<Real> product = 1;
						for(size_t j = 0; j < 4; j++) {
							const size_t jm = (j + 3) % 4;//j-1
							const size_t jp = (j + 1) % 4;//j+1
							const Real aj = std::inner_product(s[jm], s[jm]+3, s[jp], Real(0));
							const Real bj = std::inner_product(s[jm], s[jm]+3, s[j ], Real(0));
							const Real cj = std::inner_product(s[j ], s[j ]+3, s[jp], Real(0));//==b[jp]
							Real cross[3] = {//s[j] X s[jp]
								s[j][1] * s[jp][2] - s[j][2] * s[jp][1],
								s[j][2] * s[jp][0] - s[j][0] * s[jp][2],
								s[j][0] * s[jp][1] - s[j][1] * s[jp][0]
							};
							const Real dj = std::inner_product(s[jm], s[jm]+3, cross, Real(0));
							product *= std::complex<Real>(bj * cj - aj, dj);
						}

						//compute the ratio of the actual pixel solid angle to the average pixel solid angle and mirror over X==Y plane
						int factor = 1;//most pixels are completely in the north hemisphere
						if(maxX) factor *= 2;//pixels on the right edge of the grid are half below the equator
						if(maxY) factor *= 2;//pixels on the top edge of the grid are half below the equator
						omg[dim * y + x] = omg[dim * x + y] = -std::atan2(product.imag(), product.real()) * factor * invOmegaBar;//area is arg(product)
					}

					//now that the +x half of the row has been filled in copy to the -x half
					std::reverse_copy(omg + dim * y + mid, omg + dim * y + dim, omg + dim * y);
				}

				//now that the entire +y half of the grid has been filled in copy to the -y half
				std::reverse_copy(omg + mid * dim, omg + dim * dim, omg);
			}
		}

		namespace legendre {
			//@brief    : compute cosine of the latitude of each ring in the northern hemisphere (including equator)
			//@note     : southern hemisphere cosines can be computed by symmetry with cos(lat[i]) = -cos(lat[# rings - i])
			//@param dim: side length of square legendre projection
			//@param lat: location to write cos(ring latitudes)
			//@reference: Barth, W., Martin, R. S., & Wilkinson, J. H. (1967). Calculation of the eigenvalues of a symmetric tridiagonal matrix by the method of bisection. Numerische Mathematik, 9(5), 386-393.
			//@method   : legendre roots are calculated as zeros of a symmetric tridiagonal matrix as described here - https://math.stackexchange.com/questions/12160/roots-of-legendre-polynomial/12209#12209
			template <typename Real>
			void roots(const size_t dim, Real * const lat) {
				//switch to variable names in paper for convenience
				const size_t n = dim;
				Real * const x = lat;

				//don't need to build diagonal (it is all zeros)
				//build subdiagonal of legendre matrix: i / sqrt(4i^2-1)
				std::vector<Real> b(n), beta(n);
				for(size_t i = 1; i < n; i++) {
					const Real den = Real(4*i*i-1);
					b   [i] = Real(i  ) / std::sqrt(den);//subdiagonal
					beta[i] = Real(i*i) /           den ;//subdiagonal^2
				}

				size_t z = 0;//total iterations to calculate all values
				const size_t limit = n * 32;//this should be enough for a 128 bit float (long double)

				const size_t m1 = n/2;//smallest eigen value to calculate (all larger values are calculated)
				Real eps1 = std::numeric_limits<Real>::epsilon();//target precision
				Real relfeh = std::numeric_limits<Real>::epsilon();//machine epsilon

				//reference pseudocode: calculation of xmin, xmax
				beta[1-1] = b[1-1] = 0;
				Real eps2 = relfeh;
				if(eps1 <= 0) eps1 = eps2;
				eps2 = eps1 / 2 + 7 * eps2;//maximum absolute error in eigenvalues

				//reference pseudocode: inner block
				std::vector<Real> wu(m1+1, Real(0));//array to hold lower bounds (upper bounds are stored in x), fill lower bounds with 0
				Real x0 = 1;//upper bound on current eigenvalue (start with largest)
				std::fill(x + 0, x + m1+1, Real(1));//fill upper bounds
				if(1 == n%2) x[m1] = 0;//explicitly specify middle eigenvalue for odd case
				//reference pseudocode: endi
				for(int k = 0; k < m1; k++) {//loop over eigen values in reverse order calculating (largest -> smallest)
					//get initial lower bound on eigen value k
					Real xu = 0;//start with lower bound on all eigen values
					for(int i = k; i < m1; i++) {//loop over previously computed eigenvalues
						if(xu < wu[i]) {//check if previous eigenvalue had a better lower bound
							xu = wu[i];//start from the better lower bound instead
							break;//reference pseudocode: goto contin
						}//reference pseudocode: end
					}//reference pseudocode: end i
					//reference pseudocode: contin
					if(x0 > x[k]) x0 = x[k];//initial upper bound on eigen value k
					while( x0 - xu > relfeh * (std::fabs(xu) + std::fabs(x0)) * 2 + eps1) {//keep iterating until we reach our precision goal
						//update estimate and check for failure to converge
						Real x1 = (xu + x0) / 2;//current estimate of eigen value is average of upper and lower bounds
						if(z++ > limit) throw std::runtime_error("too many iterations computing legendre roots");//count iterations required

						//reference pseudocode: Sturms sequence
						size_t a = n;
						Real q = 1;
						for(size_t i = 0; i < n; i++) {
							q = - ( x1 + (q != 0 ? beta[i] / q : std::fabs(b[i]) / relfeh) );
							if(std::signbit(q)) --a;
						}//reference pseudocode: end i

						//now update appropriate bound based on a
						if(a > k) {//update the lower bound
							if(a >= m1) {
								xu = wu[n-1-m1 ] = x1;//update lower bound of smallest eigen value to calculate
							} else {
								xu = wu[a-1] = x1;//update lower bound of eigenvalue a+1
								if(x[a] > x1) x[a] = x1;//also update corresponding upper bound if possible (all eigenvalues are distinct)
							}
						} else {//update the upper bound
							x0 = x1;
						}
						x[k] = (x0 + xu) / 2;
					}//reference pseudocode: end x1
				}//reference pseudocode: end k			
			}

			//@brief    : compute the real space coordinates of unprojected points in a square legendre grid (north hemisphere only)
			//@param dim: side length of square legendre grid (must be odd)
			//@param xyz: location to write real space coordiantes (dim * dim * 3)
			template <typename Real> void normals(const size_t dim, Real * const xyz) {
				static const Real kPi_4 = Real(0.7853981633974483096156608458199);//pi/4
				if(0 == dim % 2) throw std::runtime_error("only odd side lengths are supported");

				//compute ring latitudes once
				const int half = (int)(dim / 2);
				std::vector<Real> cosLats(half);
				roots(dim - 2, cosLats.data());

				//loop over northern hemisphere compute
				for(size_t j = 0; j < dim; j++) {
					const int rj = int(j) - half;//how many rings away from the pole are we in the x direction
					const size_t aj = (size_t)std::abs(rj);//|rj|
					const Real y = (Real(j) / (dim - 1)) * 2 - 1;//[-1,1]
					for(size_t i = 0; i < dim; i++) {
						const int ri = int(i) - half;//how many rings away from the pole are we in the x direction
						const size_t ai = std::abs(ri);//|ri|
						const Real x = (Real(j) / (dim - 1)) * 2 - 1;//[-1,1]
						const size_t ar = std::max(ai, aj);//how many rings away from the pole are we

						Real n[3];
						if(ar == 0) {//handle pole specially
							n[0] = n[1] = 0;
							n[2] = 1;
						} else {//not on pole
							const Real sX = Real(ri) / ar;
							const Real sY = Real(rj) / ar;
							Real x, y;
							if(ai <= aj) {
								const Real qq = kPi_4 * sX * sY;
								x = sY * std::sin(qq);
								y = sY * std::cos(qq);
							} else {
								const Real qq = kPi_4 * sY * sX;
								x = sX * std::cos(qq);
								y = sX * std::sin(qq);
							}
							const Real h = std::hypot(x, y);
							n[2] = cosLats[ar-1];//get z from lookup table
							n[0] = n[1] = std::sqrt(Real(1) - n[2] * n[2]);//sin(acos(z))
							n[0] *= x / h;//cos(atan2(y, x))
							n[1] *= y / h;//sin(atan2(y, x))
						}
						std::copy(n, n + 3, xyz + 3 * (dim * j + i));
					}
				}
			}

			//@brief     : compute the 4 bounding indices of for a given unit direction in a the north hemisphere of square legendre grid
			//@param dim : side length of square legendre grid (must be odd)
			//@param zLat: z coordinates of each ring in the projection (e.g. from cosLats)
			//@param n   : direction to get bounding indices for
			//@param inds: location to write bounding indices
			template <typename Real> void boundingInds(const size_t dim, Real const * const zLat, Real const * const n, size_t * const inds) {
				if(0 == dim % 2) throw std::runtime_error("only odd side lengths are supported");

				//get indices of bounding rings
				const size_t ub = std::distance(zLat, std::lower_bound(zLat, zLat + dim, n[2], std::greater<Real>()) );
				const size_t rN = ub-1;//index of ring with more northern latitude than n
				const size_t rS = ub == dim ? rN : ub;//index of ring with more southern latitude than n

				//get the index of the north pole
				const size_t idxPole = (dim * dim) / 2;//index of north pole
				const bool pole = rN == 0;//are we next to the north pole?

				//compute fractional progress around ring
				const Real ax = std::fabs(n[0]);
				const Real ay = std::fabs(n[1]);
				Real theta = ( std::atan(std::min(ax,ay) / std::max(ax,ay)) * 4 ) / emsphinx::Constants<Real>::pi;//[-0.5 -> 0.5]
				if(theta != theta) theta = 0;//ax == ay == 0;

				//initialize shift left/right or up/down in pixels
				size_t delta[4] = {1,1,1,1};

				//now handle octants of square grid separately
				const bool nx = std::signbit(n[0]);
				const bool ny = std::signbit(n[1]);
				bool sub;
				if(ay >= ax) {//[45,135] or [225,315]
					//get indices of (0, +/-y)
					if(ny) {//[225,315]
						inds[0] = inds[1] = idxPole - rN * dim;//index of 270 degree point in ring rN
						inds[2] = inds[3] = idxPole - rS * dim;//index of 270 degree point in ring rS
					} else {//[45,135]
						inds[0] = inds[1] = idxPole + rN * dim;//index of  90 degree point in ring rN
						inds[2] = inds[3] = idxPole + rS * dim;//index of  90 degree point in ring rS
					}
					sub = nx;
				} else {//[-45,45] or [135,225]
					//get indices of (+/-x, 0)
					if(nx) {//[135,225]
						inds[0] = inds[1] = idxPole - rN      ;//index of 180 degree point in ring rN
						inds[2] = inds[3] = idxPole - rS      ;//index of 180 degree point in ring rS
					} else {//[-45,45]
						inds[0] = inds[1] = idxPole + rN      ;//index of   0 degree point in ring rN
						inds[2] = inds[3] = idxPole + rS      ;//index of   0 degree point in ring rS
					}
					sub = ny;
					std::for_each(delta, delta+4, [dim](size_t& d){d *= dim;});//we need to do row shifts instead of column shifts
				}

				delta[0] *= (size_t)std::ceil (theta * rN);
				delta[1] *= (size_t)std::floor(theta * rN);
				delta[2] *= (size_t)std::ceil (theta * rS);
				delta[3] *= (size_t)std::floor(theta * rS);

				if(sub)
					for(size_t i = 0; i < 4; i++) inds[i] -= delta[i];
				else
					for(size_t i = 0; i < 4; i++) inds[i] += delta[i];
			}
		}

		//@brief     : extract a single ring from a row major pattern
		//@param dim : side length of square in pixels
		//@param ring: ring number to copy (0 for north pole, (dim-1)/2 for equator)
		//@param ptr : pointer to start of hemisphere to copy from (dim * dim array in row major order)
		//@param buff: location to write extracted ring
		//@return    : number of values copied
		template <typename T>
		size_t readRing(const size_t dim, const size_t ring, T const * const ptr, T * const buff) {
			//compute indicies of corners
			const size_t even = 1 - (dim % 2);//1 if even 0 if odd
			const size_t side = 2 * ring + 1 + even;//side length of ring to copy
			const size_t pole = (dim * (dim + even) ) / 2;//index of north pole (or closest point for even sizes)
			const size_t start = pole  +       ring      ;//first point to copy
			const size_t quad1 = start + dim * ring      ;//corner in quadrant 1
			const size_t quad2 = quad1 -       (side - 1);//corner in quadrant 2
			const size_t quad3 = quad2 - dim * (side - 1);//corner in quadrant 3
			const size_t quad4 = quad3 +       (side - 1);//corner in quadrant 4

			//compute offsets to quad2 and quad3 in the ring
			const size_t b1 = ring;
			const size_t b2 = b1 + side - 1;
			const size_t b3 = b2 + side - 1;
			const size_t b4 = b3 + side - 1;

			//copy data using access pattern that is sequential in the row major order
			//this should be the best option for caching behavior since the ring major buffer is always relatively small
			std::        copy(ptr + quad3, ptr + quad4 + 1, buff + b3);//copy from quadrant 3 corner to quadrant 4 corner
			for(size_t i = 1-even; i < ring; i++) {//copy from quadrant 3/4 corner to y == 0
				buff[b3-i-even] = ptr[quad3 + dim * (       i + even)];//-x (quadrant 3 -> y == 0)
				buff[b4+i+even] = ptr[quad4 + dim * (       i + even)];//+x (quadrant 4 -> y == 0)
			}
			for(size_t i = 0; i < ring; i++) {//copy from y == 0 to quadrant 2/1 corner
				buff[b2+ring-i] = ptr[quad3 + dim * (ring + i + even)];//-x (y == 0 -> quadrant 2)
				buff[        i] = ptr[quad4 + dim * (ring + i + even)];//+x (y == 0 -> quadrant 1)
			}
			std::reverse_copy(ptr + quad2, ptr + quad1 + 1, buff + b1);//copy from quadrant 1 corner to quadrant 2 corner
			if(0 == ring && 0 == even) return 1;
			return 4 * (side - 1);
		}

		//@brief     : write a single ring into a row major pattern
		//@param dim : side length of square in pixels
		//@param ring: ring number to write (0 for north pole, (dim-1)/2 for equator)
		//@param ptr : pointer to start of hemisphere to write to (dim * dim array in row major order)
		//@param buff: location to read ring from
		//@return    : number of values copied
		template <typename T>
		size_t writeRing(const size_t dim, const size_t ring, T * const ptr, T const * const buff) {
			//compute indicies of corners
			const size_t even = 1 - (dim % 2);//1 if even 0 if odd
			const size_t side = 2 * ring + 1 + even;//side length of ring to copy
			const size_t pole = (dim * (dim + even) ) / 2;//index of north pole (or closest point for even sizes)
			const size_t start = pole  +       ring      ;//first point to copy
			const size_t quad1 = start + dim * ring      ;//corner in quadrant 1
			const size_t quad2 = quad1 -       (side - 1);//corner in quadrant 2
			const size_t quad3 = quad2 - dim * (side - 1);//corner in quadrant 3
			const size_t quad4 = quad3 +       (side - 1);//corner in quadrant 4

			//compute offsets to quad2 and quad3 in the ring
			const size_t b1 = ring;
			const size_t b2 = b1 + side - 1;
			const size_t b3 = b2 + side - 1;
			const size_t b4 = b3 + side - 1;

			//copy data using access pattern that is sequential in the row major order
			//this should be the best option for caching behavior since the ring major buffer is always relatively small
			std::        copy(buff + b3, buff + b4 + 1, ptr + quad3);//copy from quadrant 3 corner to quadrant 4 corner
			for(size_t i = 1-even; i < ring; i++) {//copy from quadrant 3/4 corner to y == 0
				ptr[quad3 + dim * (       i + even)] = buff[b3-i-even];//-x (quadrant 3 -> y == 0)
				ptr[quad4 + dim * (       i + even)] = buff[b4+i+even];//+x (quadrant 4 -> y == 0)
			}
			for(size_t i = 0; i < ring; i++) {//copy from y == 0 to quadrant 2/1 corner
				ptr[quad3 + dim * (ring + i + even)] = buff[b2+ring-i];//-x (y == 0 -> quadrant 2)
				ptr[quad4 + dim * (ring + i + even)] = buff[        i];//+x (y == 0 -> quadrant 1)
			}
			std::reverse_copy(buff + b1, buff + b2 + 1, ptr + quad2);//copy from quadrant 1 corner to quadrant 2 corner
			if(0 == ring && 0 == even) return 1;
			return 4 * (side - 1);
		}

		//@brief    : compute quadrature weights for rings (w_y in equation 10 of Reinecke)
		//@param dim: side length of square lambert projection to compute weights for
		//@param lat: cosines of ring latitudes (symmetric across equator)
		//@param wgt: location to write weights for each row
		//@param skp: ring to exclude from weights (e.g. skip = 0 will exclude the poles)
		//@reference: https://doi.org/10.1111/j.1365-246X.1994.tb03995.x
		template <typename Real>
		void computeWeightsSkip(const size_t dim, Real const * const lat, Real * const wgt, const size_t skp) {
			//compute the cosine of the latitude of each ring directly
			const size_t num = dim;//compute the number of rings
			const size_t nMat = (num + 1) / 2 - 1;//only need north hemisphere (south by symmetry)

			//build matrix of the form a_ij = cos(2*i*latitude[j]) for Sneeuw equation (21)
			//I'll compute them using Chebyshev recursion: cos(ni) = T_n(cos(i)) with: T_0(x) = 1, T_1(x) = x, T_{n}(x) = 2x T_{n-1}(x) - T_{n-2}(x)
			std::vector<Real> a(nMat * nMat);//allocate matrix
			Real * const x = a.data() + nMat;//we'll need the second row x = cos(2*theta) for every subsequent row
			std::fill(a.begin(), a.begin() + nMat, Real(1));//cos(0*theta) along first row
			for(size_t i = 0  ; i < skp ; i++) x[i] = lat[i  ] * lat[i  ] * 2 - 1;//double angle formula for cos(x) -> cos(2x)
			for(size_t i = skp; i < nMat; i++) x[i] = lat[i+1] * lat[i+1] * 2 - 1;//double angle formula for cos(x) -> cos(2x)

			//now fill in remaining rows of matrix a
			for(size_t j = 2; j < nMat; j++) {//loop over remaining rows using Chebyshev recursion to fill in cos(2*n*i)
				Real const * const t2 = a.data() + (j-2) * nMat;//row of T_{n-2}
				Real const * const t1 = a.data() + (j-1) * nMat;//row of T_{n-1}
				Real       * const tn = a.data() +  j    * nMat;//row of T_{n  }
				for(size_t i = 0; i < nMat; i++) tn[i] = x[i] * t1[i] * 2 - t2[i];//T_{n}(x) = 2x T_{n-1}(x) - T_{n-2}(x)
			}
			//build column vector b and solve for weights (A * wgt = b)
			std::vector<Real> b(nMat);
			for(int i = 0; i < nMat; i++) b[i] = Real(-1) / (4 * i * i - 1);
			solve::lu(a.data(), wgt, b.data(), nMat);//solve A * wgt = b

			//compute the solid angle of a grid point
			const size_t gridPoints = dim * dim * 2 - (dim - 1) * 4;//determine total number of points on sphere (equator has double cover)
			const Real wn = emsphinx::Constants<Real>::pi2 * 4 / gridPoints;//solid angle of single point (solid angle of sphere / # grid points)
			const Real w0 = wn * Real( ( dim - 2 ) * dim + 2 );//initial scaling factor (doesn't account for different # pts in each ring)

			//rescale weights by solid angle and ratio of points to equatorial points and use symmetry to fill southern hemisphere weights
			const bool even = (0 == dim % 2);//even and odd are slightly different
			const size_t offset = even ? 4 : 0;
			const Real delta = std::accumulate(wgt, wgt + nMat, Real(0)) - 1;//should be 0 (weights should sum to 2, 1 for only northern hemisphere)
			if(delta > std::cbrt(std::numeric_limits<Real>::epsilon()) / 64 ) throw std::runtime_error("insufficient precision to accurately compute ring wieghts");
			std::copy(wgt + skp, wgt + nMat, wgt + skp + 1);//correct alignment for missing ring
			wgt[skp] = 0;//skipped ring has no weight
			wgt[0] *= w0 / (even ? 4 : 1);//handle first ring / pole point specially (4 points for even, 1 point for odd sizes)
			for(size_t i = 1; i < nMat+1; i++) wgt[i] *= w0 / (8 * i + offset);//scale to correct for different ring sizes to w0 / # points in this ring
		}

		//@brief     : compute the cosines of ring latitudes for a given square grid
		//@param dim : side length of square projection to compute ring latitudes for
		//@param type: type of square projection to compute latitudes for
		//@return    : cosine(latitudes) from north pole -> equator
		template <typename Real> typename std::vector<Real> cosLats(const size_t dim, const Layout type) {
			std::vector<Real> cLats((dim+1) / 2);//cosines of ring latitudes
			switch(type) {
				case Layout::Lambert:
					lambert::cosLats(dim, cLats.data());
					break;

				case Layout::Legendre:
					cLats[0] = 1;//manually specify pole
					legendre::roots(dim-2, cLats.data()+1);//fill in north ring latitude cosines
			}
			return cLats;
		}

		//@brief     : compute the normals of the northern hemisphere for a given square grid
		//@param dim : side length of square projection to compute ring latitudes for
		//@param type: type of square projection to compute latitudes for
		//@return    : {x,y,z} normals for north pole
		template <typename Real> typename std::vector<Real> normals(const size_t dim, const Layout type) {
			std::vector<Real> normals(dim * dim * 3);//cosines of ring latitudes
			switch(type) {
				case Layout ::Lambert:
					lambert ::normals(dim, normals.data());
					break;

				case Layout::Legendre:
					legendre::normals(dim, normals.data());
					break;
			}
			return normals;
		}

		//@brief     : compute the solid angle correction of each ring for a given square grid
		//@param dim : side length of square lambert projection to compute solid angles for
		//@param type: type of square projection to compute latitudes for
		//@return    : solide angle correction from north pole -> equator (actual pixel size / average pixel size)
		template <typename Real> typename std::vector<Real> solidAngles(const size_t dim, const Layout type) {
			//compute cosines of ring lattidues and average pixel size
			std::vector<Real> cLats = cosLats<Real>(dim, type);
			const size_t numPix = dim * dim * 2 - (dim - 1) * 4;
			// const Real avgPix = emsphinx::Constants<Real>::pi2 * Real(2) / numPix;//4 pi evenly split over pixels
			const Real avgPix = Real(2) / numPix;//2 evenly split over pixels (2 pi will cancel out)

			//convert form cosine lattidues of center to cosine latitudes of mid points
			std::vector<Real> cSplits(cLats.size() - 1);
			std::transform(cLats.cbegin(), cLats.cend() - 1, cLats.cbegin() + 1, cSplits.begin(), [](const Real& cA, const Real& cB) {
				return ( std::sqrt( (Real(1) + cA) * (Real(1) + cB) ) - std::sqrt( (Real(1) - cA) * (Real(1) - cB) ) ) / Real(2);//cos((a+b)/2) where a = acos(cA), b = acos(cB)
			});
			cSplits.push_back(-cSplits.back());//symmetric about equator 

			//convert from cosine latitudes of midpoints to spherical cap areas
			for(Real& c : cSplits) c = (Real(1) - c);// * emsphinx::Constants<Real>::pi2; (we're just going to divide by 2 pi later)

			//convert from spherical cap areas to area of ring
			std::adjacent_difference(cSplits.cbegin(), cSplits.cend(), cSplits.begin());

			//convert from ring area to pixel area and normalize by average pixel size
			size_t num = 8;
			if(0 == dim % 2) {//even [4,12,20]
				cSplits[0] /= avgPix * Real(4);
				num = 12;
			} else {//odd [1,8,16,24,...]
				cSplits[0] /= avgPix * Real(1);
			}
			for(size_t i = 1; i < cSplits.size(); i++) {
				cSplits[i] /= avgPix * Real(num);
				num += 8;
			}
			return cSplits;
		}

		//@brief    : convert from a vecotrized index to a ring number
		//@param dim: side length of square lambert projection to compute index for
		//@param idx: vectorized index (y * dim + x) to compute ring number for [must be in north hemisphere i.e. < dim * dim]
		//@return   : ring number (e.g. for indexing into cosLats)
		size_t ringNum(const size_t dim, const size_t idx) {
			//convert from vectorized to x/y indices
			const int y = int(idx / dim);
			const int x = int(idx - y * dim);

			//compute distance to center pixel
			int d2 = int(dim / 2);
			size_t dx = (size_t)std::abs(x - d2);//x distance from center
			size_t dy = (size_t)std::abs(y - d2);//y distance from center
			if(0 == dim % 2) {
				//even case is a little more complicated since the center isn't on a pixel
				++d2;//move to other corner of northern most ring
				dx = std::min<size_t>(dx, (size_t)std::abs(x - d2));//choose whichever x distance is smallest
				dy = std::min<size_t>(dy, (size_t)std::abs(y - d2));//choose whichever x distance is smallest
			}
			return std::max(dx, dy);//ring is maximum distance
		}
	}

}

#endif//_SQUARE_SHT_H_
