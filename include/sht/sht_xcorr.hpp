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

#ifndef _SHT_XCORR_H_
#define _SHT_XCORR_H_

#include <vector>
#include <complex>
#include <memory>

#include "util/fft.hpp"

namespace emsphinx {

	namespace sphere {
		//@brief: helper struct to compute the cross correlation of 2 spherical functions from their spherical harmonic transformations
		//@reference: Gutman, B., Wang, Y., Chan, T., Thompson, P. M., & Toga, A. W. (2008, October). Shape registration with spherical cross correlation. In 2nd MICCAI workshop on mathematical foundations of computational anatomy (pp. 56-67).
		//@note: the reference is for the generic case (complex functions) but restricting to real valued functions allow for savings from symmetry:
		//    \hat{f}(l, -m) = \hat{f}(l, m) * (-1)^m
		//  additionally since the decomposition of the wigner D function to 2 wigner d functions @ pi/2 introduces more symmetry
		//    d^j_{-k,-m} = (-1)^(   k- m) d^j_{k,m}
		//    d^j_{ k,-m} = (-1)^(j+ k+2m) d^j_{k,m}
		//    d^j_{-k, m} = (-1)^(j+2k+3m) d^j_{k,m}
		//    d^j_{ m, k} = (-1)^(   k- m) d^j_{k,m}
		//  finally since the cross correlation of the real functions is also real there is another factor of 2 savings
		//@note: variable names are consistent with Gutman et. al. (see eq 12 for details) except 'j' is used in place of 'l'
		template <typename Real>
		struct Correlator {
			//@brief          : construct a spherical correlator for a given bandwidth
			//@param bandWidth: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
			Correlator(const size_t bandWidth);

			//@brief    : compute the rotation of the maximum cross correlation between two spherical functions
			//@param flm: spherical harmonic coefficients for the first function
			//@param gln: spherical harmonic coefficients for the second function
			//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
			//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
			//@param eu : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param ref: true/false to use/not use real space refinement
			//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return   : maximum cross correlation
			Real correlate(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const eu, const bool ref = true, const Real eps = Real(0.01));

			//@brief     : compute a subpixel maxima in the cross correlation using interpolation
			//@param ind0: pixel to search from (should be near a local maxima)
			//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@return    : estimated cross correlation of peak
			Real interpPeak(size_t ind0, Real * const eu);
			
			//@brief     : compute a subpixel maxima in the cross correlation by either interpolation or refinement
			//@param flm : spherical harmonic coefficients for the first function
			//@param gln : spherical harmonic coefficients for the second function
			//@param fMr : true/false if there is/isn't a mirror plane in the first fuction
			//@param fNf : rotational symmetry about z axis in first function (1 for no rotational symmetry)
			//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return    : cross correlation of peak
			Real refinePeak(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const eu, const Real eps = Real(0.01));
			
			//@brief     : extract an (N+1)^3 neighborhood around a pixel using periodic boundary conditions
			//@template N: half window size of neighborhood to extract
			//@param idx : vectorized index to extract neighborhood around
			//@param nh  : 3D array of size (N+1)^3 to write neighborhood into
			template <size_t N> void extractNeighborhood(const size_t idx, Real (&nh)[N*2+1][N*2+1][N*2+1]) const;

			//@brief   : compute the index of the pixel closest to an orientation
			//@param eu: ZYZ euler angle to find nearest pixel of
			//@return  : vectorized pixel index
			size_t eulerIndex(Real const * const eu) const;

			//@brief    : compute the orientation of an index
			//@param idx: vectorized pixel index
			//@param eu : location to write ZYZ euler angles
			void indexEuler(const size_t idx, Real * const eu) const;

			//@brief : get the cross correlation from the previous call to correlate
			//@return: cross correlation grid
			const fft::vector<Real>& getXC() const {return xc;}

			//@brief    : copy the cross correlation from the previous call and convert to ZXZ euler angles with origin at 0
			//@param zxz: zxz euler angle grid to write cross correlation to
			//@note     : phi1 (rotation about Z) increments fastest, phi2 (rotation about Z'') middle, and Phi (rotation about X') slowest
			void extractBunge(Real * const zxz) const;

			//@brief : get the maximum bandwidth
			//@return: maximum bandwidth
			size_t getBw() const {return bw;}

			//@brief : get the cross correlation grid size
			//@return: side length of correlation cube (at least 2 * bw - 1, may be zero padded larger)
			size_t getCubeSize() const {return slP;}

			//@brief : get half the cross correlation grid size
			//@return: side length of correlation cube (at least 2 * bw - 1, may be zero padded larger)
			size_t getHalfSize() const {return bwP;}

			//@brief    : compute the first and second derivatives of the cross correlation at a single rotation
			//@param flm: spherical harmonic coefficients for the first function
			//@param gln: spherical harmonic coefficients for the second function
			//@param eu : rotation to compute derivatives of cross correlation for as ZYZ euler angle
			//@param jac: location to write jacobian of cross correlation {d/(d eu[0]), d/(d eu[1]), d/(d eu[2])}
			//@param hes: location to write hessian (3x3 matrix as 9 component vector) of cross correlation hes_ij = d/(d eu[i]) * d/(d eu[j])
			//@param mBW: maximum bandwidth to use in calculation (must be <= bw)
			//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
			//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
			//@param der: true/false to compute derivatives/only cross correlation
			//@return   : cross correlation for rotation eu
			Real derivatives(std::complex<Real> const * const flm, std::complex<Real> const * const gln, Real const * const eu, Real * const jac, Real * const hes, const size_t mBW, const bool fMr, const size_t fNf, const bool der = true);

			protected:

				//@brief    : compute the cross correlation between two spherical functions
				//@param flm: spherical harmonic coefficients for the first function
				//@param gln: spherical harmonic coefficients for the second function
				//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
				//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
				//@param pXc: location to write cross correlation
				void compute(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const pXc);

				//@brief : find the maximum cross correlation grid point
				//@return: index of maximum cross correlation from previous call to correlate
				size_t findPeak();

				struct Constants;//helper struct to hold read only constants
				const size_t                           bw   ;//maximum bandwidth to use (exclusive)
				const size_t                           sl   ;//side length of grid in euler space (2 * bandWidth - 1)
				const size_t                           slP  ;//zero padded side length of grid in euler space (2 * bandWidth - 1)
				const size_t                           bwP  ;//zero padded bandwidth (slP + 1) / 2
				const std::shared_ptr<const Constants> xcLut;//read only values (can be shared across threads)
				std::vector< std::complex<Real> >      fm   ;//2d lookup table to hold \hat{f}(j,m) * d^j_{m, k} for all j and m (for any given k)
				std::vector< std::complex<Real> >      gn   ;//1d lookup table to hold \bar{\hat{g}(j,n)} * d^j_{k,n} for all j (for any given k and n)
				fft::vector< std::complex<Real> >      fxc  ;//fft of cross correlation (in fftw multi dimensional real format)
				fft::vector<              Real  >      xc   ;//real space cross correlation (this is what we're after)
				std::vector<              Real  >      dBeta;//wigner (lowercase) d lookup table for arbitrary beta (for realspace refinement)
		};

		//@brief: intermediate class to impose condition that makes abstracting (un)normalized easier
		template <typename Real>
		struct PhaseCorrelator : protected Correlator<Real> {
			//@brief   : require bandwidth for construction
			//@param bw: bandwidth
			PhaseCorrelator(const size_t bw) : Correlator<Real>(bw) {}

			//@brief: default destructor (for unique_ptr)
			virtual ~PhaseCorrelator() = default;

			//@brief : get a copy of the stored spherical correlator
			//@return: unique pointer to copy of current correlator
			virtual std::unique_ptr<PhaseCorrelator> clone() const = 0;

			//@brief    : compute the rotation of the maximum (un)normalized cross correlation between two spherical functions
			//@param gln: spherical harmonic coefficients for the template function
			//@param eu : location to write rotation of maximum normalized cross correlation as ZYZ euler angle
			//@param ref: true/false to use/not use real space refinement
			//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return   : maximum (un)normalized cross correlation
			//@note     : flm, fMr, and fNf stored in derived class
			virtual Real correlate(std::complex<Real> const * const gln, Real * const eu, const bool ref = true, const Real eps = Real(0.01)) = 0;

			//@brief     : compute a subpixel maxima in the (un)normalized cross correlation by either interpolation or refinement
			//@param gln : spherical harmonic coefficients for the second function
			//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return    : cross correlation of peak
			virtual Real refinePeak(std::complex<Real> const * const gln, Real * const eu, const Real eps = Real(0.01)) = 0;
		};

		//@brief    : dervied class to compute un-normalized cross correlation
		template <typename Real>
		struct UnNormalizedCorrelator : public PhaseCorrelator<Real> {
			std::shared_ptr< std::vector< std::complex<Real> > > flm;//actual harmonics
			const bool                                           fMr;//z mirror flag
			const size_t                                         fNf;//n fold symmetry

			//@brief          : construct a spherical correlator for a given bandwidth
			//@param bandWidth: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
			//@param flm      : spherical harmonic coefficients for the reference function
			//@param fMr      : true/false if there is/isn't a mirror plane in the reference fuction
			//@param fNf      : rotational symmetry about z axis in reference function (1 for no rotational symmetry)
			UnNormalizedCorrelator(const size_t bandWidth, std::shared_ptr< std::vector< std::complex<Real> > > flm, const bool fMr, const size_t fNf) : PhaseCorrelator<Real>(bandWidth), flm(flm), fMr(fMr), fNf(fNf) {}

			//@brief    : compute the rotation of the maximum cross correlation between two spherical functions
			//@param gln: spherical harmonic coefficients for the template function
			//@param eu : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param ref: true/false to use/not use real space refinement
			//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return   : maximum
			Real correlate(std::complex<Real> const * const gln, Real * const eu, const bool ref = true, const Real eps = Real(0.01)) {return Correlator<Real>::correlate(flm->data(), gln, fMr, fNf, eu, ref, eps);}

			//@brief     : compute a subpixel maxima in the cross correlation by either interpolation or refinement
			//@param gln : spherical harmonic coefficients for the second function
			//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return    : cross correlation of peak
			Real refinePeak(std::complex<Real> const * const gln, Real * const eu, const Real eps = Real(0.01)) {return Correlator<Real>::refinePeak(flm->data(), gln, fMr, fNf, eu, eps);}

			//@brief : get a copy of the stored spherical correlator
			//@return: unique pointer to copy of current correlator
			std::unique_ptr<PhaseCorrelator<Real> > clone() const {return std::unique_ptr<PhaseCorrelator<Real> >(new UnNormalizedCorrelator(*this));}
		};

		//@brief    : dervied class to compute normalized cross correlation
		//@reference: Huhle, B., Schairer, T., and Strasser, W. (2009) Normalized Cross-Correlation Using SOFT
		//@note     : this class assumes the refernce function has full support
		template <typename Real>
		struct NormalizedCorrelator : public PhaseCorrelator<Real> {
			struct Constants;//helper struct to hold read only constants
			const std::shared_ptr<const Constants> ncLut;//read only values (can be shared across threads)

			//@brief          : construct a spherical correlator for a given bandwidth
			//@param bandWidth: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
			//@param flm      : spherical harmonic coefficients for the reference function
			//@param flm2     : spherical harmonic coefficients for the reference (function^2)
			//@param fMr      : true/false if there is/isn't a mirror plane in the reference fuction
			//@param fNf      : rotational symmetry about z axis in reference function (1 for no rotational symmetry)
			//@param mlm      : spherical harmonic coefficients for the mask function
			NormalizedCorrelator(const size_t bandWidth, std::complex<Real> const * const flm, std::complex<Real> const * const flm2, const bool fMr, const size_t fNf, std::complex<Real> const * const mlm);

			//@brief    : compute the rotation of the maximum normalized cross correlation between two spherical functions
			//@param gln: spherical harmonic coefficients for the template function
			//@param eu : location to write rotation of maximum normalized cross correlation as ZYZ euler angle
			//@param ref: true/false to use/not use real space refinement
			//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return   : maximum (semi)normalized cross correlation [still needs to be divided by the standard deviation of the pattern function]
			Real correlate(std::complex<Real> const * const gln, Real * const eu, const bool ref = true, const Real eps = Real(0.01));

			//@brief     : compute a subpixel maxima in the normalized cross correlation by either interpolation or refinement
			//@param gln : spherical harmonic coefficients for the second function
			//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
			//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
			//@return    : cross correlation of peak
			//@note      : to be fully robust this should use the chain rule to account for the effects of shifting the window (which it doesn't currently do)
			//             not using the chain rule assumes that the effect of the window is small near the peak
			Real refinePeak(std::complex<Real> const * const gln, Real * const eu, const Real eps = Real(0.01));

			//@brief : get a copy of the stored spherical correlator
			//@return: unique pointer to copy of current correlator
			std::unique_ptr<PhaseCorrelator<Real> > clone() const {return std::unique_ptr<PhaseCorrelator<Real> >(new NormalizedCorrelator(*this));}
		};
	}

}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cmath>

#include "constants.hpp"
#include "wigner.hpp"
#include "util/linalg.hpp"

namespace emsphinx {

	namespace sphere {
		//helper struct to hold read only constants needed by Correlator (for sharing across threads)
		template <typename Real>
		struct Correlator<Real>::Constants {
			//these members could be shared across threads since they should be read only after construction
			// const size_t                      bw  ;//maximum bandwidth to use (exclusive)
			// const size_t                      sl  ;//side length of grid in euler space (2 * bandWidth - 1)
			std::vector<              Real  > wigD;//lookup table for wigner d functions
			fft::SepRealFFT3D<Real>           plan;//fftw plan to compute Real cross correlation from DFT(cross correlation
			std::vector<              Real  > wigE;//lookup table for dTablePre
			std::vector<              Real  > wigW;//lookup table for dTablePre
			std::vector<              Real  > wigB;//lookup table for dTablePre

			//@brief   : construct a set of constants for a given max bandwidth
			//@param bw: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
			//@param fs: fft size (cube side length)
			Constants(const size_t bw, const size_t fs);
		};

		//helper struct to hold additional read only constants needed by NormalizedCorrelator (for sharing across threads)
		template <typename Real>
		struct NormalizedCorrelator<Real>::Constants {
			//these members could be shared across threads since they should be read only after construction
			std::vector<              Real  > rDen;//reciprocal of denominator for normalization over euler angle grid (computed via equation 8 in reference)
			fft::vector< std::complex<Real> > flm ;//spherical harmonic coefficients for the reference function
			fft::vector< std::complex<Real> > flm2;//spherical harmonic coefficients for the reference (function^2)
			fft::vector< std::complex<Real> > mlm ;//spherical harmonic coefficients for the mask function
			bool                              mr  ;//does flm have a mirror plane at the equator
			size_t                            nf  ;//rotational symmetry of flm about z axis

			//@brief     : construct a spherical correlator for a given bandwidth
			//@param cor : correlator object
			//@param flm : spherical harmonic coefficients for the reference function
			//@param flm2: spherical harmonic coefficients for the reference (function^2)
			//@param fMr : true/false if there is/isn't a mirror plane in the reference fuction
			//@param fNf : rotational symmetry about z axis in reference function (1 for no rotational symmetry)
			//@param mlm : spherical harmonic coefficients for the mask function
			Constants(NormalizedCorrelator<Real>* cor, std::complex<Real> const * const flm, std::complex<Real> const * const flm2, const bool fMr, const size_t fNf, std::complex<Real> const * const mlm);

			//@brief    : compute the normalization denominator at an arbitrary rotation
			//@param cor: correlator object
			//@param eu : rotation to compute normalization for
			//@return   : normalization
			Real denominator(Correlator<Real>& cor, Real * const eu) const;
		};

		namespace detail {
			//@brief   : compute product of ab * cd and ad * conj(cd) without duplicate flops
			//@param ab: first  complex number (a + bi)
			//@param cd: second complex number (c + di)
			//@param vp: location to write ab *      cd
			//@param vc: location to write ab * conj(cd)
			template <typename T>
			inline void conjMult(const std::complex<T>& ab, const std::complex<T>& cd, std::complex<T>& vp, std::complex<T>& vc);

			//@brief    : convert a vectorized index in the 3d cross correlation grid to individual components
			//@param idx: vectorized index
			//@param sl : side length of cube
			//@param knm: location to write indices (z,y,x)
			inline void extractInds(size_t idx, const size_t sl, size_t knm[3]);

			//@brief: interpolate subpixel peak location from a 3d voxel grid
			//@param p: neighborhood around peak
			//@param x: location to store subpixel maxima location within neighborhood (z, y, x from -1->1)
			//@param return: value of fit quadratic at maxima
			template <typename Real, typename PixelType>
			Real interpolateMaxima(PixelType p[3][3][3], Real x[3]);
		}

		//@brief   : construct a set of constants for a given max bandwidth
		//@param bw: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
		//@param fs: fft size (cube side length)
		template <typename Real>
		Correlator<Real>::Constants::Constants(const size_t bw, const size_t fs) : 
			wigD(bw * bw * bw),
			wigE(bw * bw     ),
			wigW(bw * bw * bw),
			wigB(bw * bw * bw),
			plan(fs, fft::flag::Plan::Patient)//the fft is the slowest part, find the best possible algorithm
		{
				wigner::dTable<Real>(bw, wigD.data(), true);
				wigner::dTablePreBuild<Real>(bw, wigE.data(), wigW.data(), wigB.data());
		}

		template <typename Real>
		Correlator<Real>::Correlator(const size_t bandWidth) :
			bw   (bandWidth                                 ),
			sl   (2 * bw - 1                                ),
			slP  ((size_t)fft::fastSize((uint32_t)sl)       ),
			bwP  ( (slP) / 2 + 1                            ),
			xcLut(std::make_shared<const Constants>(bw, slP)),
			fm   (bw * bw                                   ),
			gn   (bw                                        ),
			fxc  (slP * slP * bwP                           ),
			xc   (slP * slP * bwP                           ),
			dBeta(bw  * bw  * bw * 2                        ) {}

		//@brief: compute the cross correlation between two spherical functions
		//@param flm: spherical harmonic coefficients for the first function
		//@param gln: spherical harmonic coefficients for the second function
		//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
		//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
		//@param eu : location to write rotation of maximum cross correlation as ZYZ euler angle
		//@param ref: true/false to use/not use real space refinement
		//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
		//@return   : maximum cross correlation
		template <typename Real>
		Real Correlator<Real>::correlate(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const eu, const bool ref, const Real eps) {
			compute(flm, gln, fMr, fNf, xc.data());//compute cross correlation on euler grid
			const size_t ind0 = findPeak();//find maximum cross correlation
			const Real peak = interpPeak(ind0, eu);
			return ref ? refinePeak(flm, gln, fMr, fNf, eu, eps) : peak;//refine peak if needed
		}

		//@brief     : compute a subpixel maxima in the cross correlation using interpolation
		//@param ind0: pixel to search from (should be near a local maxima)
		//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
		//@return    : estimated cross correlation of peak
		template <typename Real>
		Real Correlator<Real>::interpPeak(size_t ind0, Real * const eu) {
			//convert from vectored index to k,n,m indices
			size_t knm[3];
			detail::extractInds(ind0, slP, knm);

			//compute indicies of neighborhood around peak with periodic boundary conditions
			const size_t N = 1;//half kernel size
			const size_t W = N * 2 + 1;//full kernel size
			Real xCorr[W][W][W];
			extractNeighborhood<N>(ind0, xCorr);

			//sub pixel interpolate peak using tri quadratic
			Real x[3] = {0,0,0};//initial guess at maximum voxel (z,y,x)
			Real peak0 = detail::interpolateMaxima(xCorr, x);
			const Real xMax = std::max(std::fabs(x[0]), std::max(std::fabs(x[1]), std::fabs(x[0])));
			if(xMax > Real(1)) {//dont step too far in case we're near a degeneracy
				std::fill(x, x + 3, Real(0));
				peak0 = xCorr[N][N][N];
			}

			//convert subpixel to ZYZ euler angle
			eu[0] = ((Real(knm[2]) + x[2])*4 - slP) * emsphinx::Constants<Real>::pi / (2 * slP);//alpha
			eu[1] = ((Real(knm[0]) + x[0])*2 - slP) * emsphinx::Constants<Real>::pi / (    slP);//beta
			eu[2] = ((Real(knm[1]) + x[1])*4 - slP) * emsphinx::Constants<Real>::pi / (2 * slP);//gamma
			return peak0;
		}

		//@brief     : compute a subpixel maxima in the cross correlation by either interpolation or refinement
		//@param flm : spherical harmonic coefficients for the first function
		//@param gln : spherical harmonic coefficients for the second function
		//@param fMr : true/false if there is/isn't a mirror plane in the first fuction
		//@param fNf : rotational symmetry about z axis in first function (1 for no rotational symmetry)
		//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
		//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
		//@return    : cross correlation of peak
		template <typename Real>
		Real Correlator<Real>::refinePeak(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const eu, const Real eps) {
			//initial setup
			Real eu0[3] = {eu[0], eu[1], eu[2]};//save a copy of the original orientation
			const Real   absEps  = eps * emsphinx::Constants<Real>::pi2 / slP     ;//compute stopping criterion is factor of grid resolution
			const Real   euEps   = std::sqrt(std::numeric_limits<Real>::epsilon());//epsilon for eu[1] being too close to 0 or pi
			const size_t maxIter = 15;//in tests of perfect rotations convergence generally takes at most 3 iterations (some edge cases near degeneracies are slower)
			Real hes[9], jac[3], step[3];//working space for hessian, jacobian, and euler update
			Real prevMag2 = emsphinx::Constants<Real>::pi2 * 3 / slP;//first step better not be more than 1 pixel in each direction

			//newton's method refinement
			try {
				Real peak;
				for(size_t iter = 1; iter <= maxIter; iter++) {
					//build hessian and derivatives
					peak = derivatives(flm, gln, eu, jac, hes, bw, fMr, fNf);

					//the hessian may be (nearly) non positive definate 
					try {//standard newton update computation
						if(std::isnan(hes[4])) throw std::runtime_error("degenerate");
						solve::cholesky(hes, step, jac, 3);//solve for step (this will fail for non definate matricies but thats good since we don't want a saddle point)
						const Real mag2 = step[0]*step[0] + step[1]*step[1] + step[2]*step[2];//compute step size
						if(mag2 > prevMag2) throw std::runtime_error("newton steps must always decrease in magnitude");//so far only a observed near degeneracies
						prevMag2 = mag2;//update step size
					} catch(...) {//there was a problem with the newton update
						//handle degeneracy near eu[1] = 0 and +/-pi
						if(std::isnan(hes[4])) {
							//if we're actually on the degeneracy the beta^2 derivative is undefined, do the 1x1 sub problem
							step[0] = jac[0] / hes[0];
							step[1] = 0;
							step[2] = 0;
						} else {//not degenerate but probably close
							//if we're just close to the degeneracy sovle 2x2 sub problem manually (3rd euler angle is false DoF)
							const Real det = hes[0] * hes[4] - hes[1] * hes[1];
							if(det < euEps) {
								if(std::fabs(det) < euEps) throw std::runtime_error("singular matrix during newton iteration");//don't divide by 0 this should be extremely rare
								if(det < euEps) throw std::runtime_error("converging to saddle during newton iteration");//it isn't clear how often this can happen (limiting the step size prevents occurance in all testing)
							}
							step[0] = (jac[0] * hes[4] - jac[1] * hes[1]) / det;
							step[1] = (jac[1] * hes[0] - jac[0] * hes[1]) / det;
							step[2] = 0;
						}
					}

					//note that we've computed the step, apply and check for convergence
					std::transform(eu, eu + 3, step, eu, std::minus<Real>());//apply step
					std::transform(step, step + 3, step, [](const Real& i){return std::fabs(i);});//convert steps to absolute value
					if(*std::max_element(step, step + 3) < absEps) break;//check for convergence
					if(iter == maxIter) throw std::runtime_error("failed to converge during cross correlation refinement");
				}
				return peak;

			} catch(...) {
				//there are still occasional failures due to the degeneracies behaving as saddle points
				std::copy(eu0, eu0 + 3, eu);//use original orientation
				return derivatives(flm, gln, eu, NULL, NULL, bw, fMr, fNf, false);//compute cross correlation at eu
			}
		}

		//@brief     : extract an (N+1)^3 neighborhood around a pixel using periodic boundary conditions
		//@template N: half window size of neighborhood to extract
		//@param idx : vectorized index to extract neighborhood around
		//@param nh  : 3D array of size (N+1)^3 to write neighborhood into
		template <typename Real>
		template <size_t N>
		void Correlator<Real>::extractNeighborhood(const size_t idx, Real (&nh)[N*2+1][N*2+1][N*2+1]) const {
			//convert from vectored index to k,n,m indices
			size_t knm[3];
			detail::extractInds(idx, slP, knm);

			//compute indicies of neighborhood around peak with periodic boundary conditions
			const size_t W = N * 2 + 1;//full kernel size
			size_t inds[3][W];
			for(size_t i = 0; i < 3; i++) {
				inds[i][N] = knm[i];
				for(size_t j = 0; j < N; j++) {
					//fill in -j indicies in direction i with periodic boundary conditions
					inds[i][N-1-j] = inds[i][N-j] == 0 ? slP - 1 : inds[i][N-j] - 1;

					//fill in +j indicies in direction i with periodic boundary conditions
					inds[i][N+1+j] = inds[i][N+j] + 1;
					if(inds[i][N+1+j] == slP) inds[i][N+1+j] = 0;
				}
			}

			//now use glide plane to bring top half points to bottom half
			for(size_t i = 0; i < W; i++) {
				if(inds[0][i] >= bwP) {//we're in the upper half of the euler cube
					inds[2][i] = inds[2][i] < bwP ? inds[2][i] + bwP - 1 : inds[2][i] - bwP;//glide alpha
					inds[1][i] = inds[1][i] < bwP ? inds[1][i] + bwP - 1 : inds[1][i] - bwP;//glide gamma
					inds[0][i] = slP - inds[0][i];//glide beta
				}
			}

			//extract cross correlation from 3x3x3 grid around peak
			for(size_t k = 0; k < W; k++) {
				for(size_t n = 0; n < W; n++) {
					for(size_t m = 0; m < W; m++) {
						nh[k][n][m] = xc[inds[0][k]*slP*slP + inds[1][n]*slP + inds[2][m]];
					}
				}
			}
		}

		//@brief   : compute the index of the pixel closest to an orientation
		//@param eu: ZYZ euler angle to find nearest pixel of
		//@return  : vectorized pixel index
		template <typename Real>
		size_t Correlator<Real>::eulerIndex(Real const * const eu) const {
			//convert from euler angle to fractional indices
			Real knmR[3] = {
				((eu[1] * (    slP)) / emsphinx::Constants<Real>::pi + slP ) / 2,//beta  -> k
				((eu[2] * (2 * slP)) / emsphinx::Constants<Real>::pi + slP ) / 4,//gamma -> n
				((eu[0] * (2 * slP)) / emsphinx::Constants<Real>::pi + slP ) / 4,//alpha -> m
			};

			//move beta -> [-pi,0] if needed
			const Real sl2 = Real(slP) / 2;
			if(knmR[0] > sl2) {//we're in the upper half of the euler cube
				knmR[0] = Real(slP) - knmR[0];//glide beta
				knmR[1] = std::fmod(knmR[1] + sl2, slP);;//glide gamma
				knmR[2] = std::fmod(knmR[2] + sl2, slP);;//glide gamma
			}

			//convert to nearest neighbors
			const size_t knm[3] = {
				(size_t) std::round(knmR[0]),
				(size_t) std::round(knmR[1]),
				(size_t) std::round(knmR[2]),
			};

			//vectorize
			return knm[0] * slP * slP + knm[1] * slP + knm[2];
		}

		//@brief    : compute the orientation of an index
		//@param idx: vectorized pixel index
		//@param eu : location to write ZYZ euler angles
		template <typename Real>
		void Correlator<Real>::indexEuler(const size_t idx, Real * const eu) const {
			//convert from vectorized to k,n,m
			size_t knm[3];
			detail::extractInds(idx, slP, knm);

			//convert from knm to alpha,beta, gamma
			eu[0] = (Real(knm[2])*4 - slP) * emsphinx::Constants<Real>::pi / (2 * slP);//alpha = 2 * pi * m / slP - pi / 2
			eu[1] = (Real(knm[0])*2 - slP) * emsphinx::Constants<Real>::pi / (    slP);//beta  = 2 * pi * k / slP - pi
			eu[2] = (Real(knm[1])*4 - slP) * emsphinx::Constants<Real>::pi / (2 * slP);//gamma = 2 * pi * n / slP - pi / 2
		}

		//@brief    : copy the cross correlation from the previous call and convert to ZXZ euler angles with origin at 0
		//@param zxz: zxz euler angle grid to write cross correlation to
		//@note     : phi1 (rotation about Z) increments fastest, phi2 (rotation about Z'') middle, and Phi (rotation about X') slowest
		template <typename Real>
		void Correlator<Real>::extractBunge(Real * const zxz) const {
			/*
			//this is extremely poorly implemented as written to keep things obvious
			//this hasn't been updated since zero padding was added
			for(size_t k = 0; k < bw; k++) {
				Real Phi = Real(k) / sl;//[0,0.5)
				for(size_t m = 0; m < sl; m++) {
					Real phi2 = Real(m) / sl;//[0,1)
					phi2 -= Real(0.25);//convert from zxz -> zyz: [-0.25,0.75)
					for(size_t n = 0; n < sl; n++) {
						Real phi1 = Real(n) / sl;//[0,1)
						phi1 += Real(0.25);//convert from zxz -> zyz: [0.25,1.25)

						//we now have fractional zyz euler angles, convert to fractional indices (correct for origin)
						double rm = phi1 + Real(0.25);//[0.5,1.5)
						double rk = Phi  + Real(0.5 );//[0.5,1.0)
						double rn = phi2 + Real(0.25);//[0.0,1.0)

						//now we have indices in the +beta half but we only have the -beta cube, move to bottom half using glide plane
						rm = rm - Real(0.5);//glide alpha: [ 0.0,1.0)
						rk = Real(1.0) - rk;//glide beta : ( 0.0,0.5]
						rn = rn - Real(0.5);//glide gamma: [-0.5,0.5)
						if(std::signbit(rn)) rn += Real(1.0);//[0.0,1.5)

						//finally convert to integer indices in our knm grid
						size_t ik = (size_t)std::round(rk * sl);//, sl-1);
						size_t in = (size_t)std::round(rn * sl);//, bw-1);//if rk == 0.5 we'll round up to bw but we can actually just mirror back instead of glide back
						size_t im = (size_t)std::round(rm * sl);//, sl-1);
						if(ik == bw) ik = bw - 1; 
						if(in == sl) in = 0; 
						if(im == sl) in = 0; 
						if(ik >= bw || in >= sl || im >= sl) continue;

						zxz[k * sl * sl + m * sl + n] = xc[ik * sl * sl + im * sl + in];
					}
				}
			}
			*/

			//slightly more difficult to understand but significantly more efficient
			//this isn't exactly correct since there is a half pixel shift (I picked the rounding direction that puts high correlation near the origin instead of the opposite corner)
			Real * pOut = zxz;//copy pointer to output start
			auto pIn = xc.cbegin() + bwP * slP * slP;//get pointer past last slice
			for(size_t k = 0; k < bwP; k++) {//loop over z slices
				pIn -= slP * slP;//move to previous slice start
				for(size_t m = 0; m < slP; m++) {//loop over rows
					const size_t im = ( slP - 1 - m + bwP + 1) % slP;//get index of row to copy
					auto iter = pIn + im * slP;//get pointer to input row start
					*pOut = *iter;
					std::reverse_copy(iter + 1, iter + slP, pOut + 1);//mirror copy into output row
					pOut += slP;//increment output
				}
			}
		}

		//@brief    : compute the cross correlation between two spherical functions
		//@param flm: spherical harmonic coefficients for the first function
		//@param gln: spherical harmonic coefficients for the second function
		//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
		//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
		//@param pXc: location to write cross correlation
		template <typename Real>
		void Correlator<Real>::compute(std::complex<Real> const * const flm, std::complex<Real> const * const gln, const bool fMr, const size_t fNf, Real * const pXc) {
			//naive implementation (no symmetry) to make summation clear
			//this also assumes no zero padding
			/*
			std::fill(fxc.begin(), fxc.end(), std::complex<Real>(0));
			const bool realFft = true;//true/false to use half sized fft format
			const size_t dm = realFft ? bw : sl;//length of fastest indexing dimension
			for(size_t ic = 0; ic < sl; ic++) {
				const int k = ic >= bw ? int(ic) - sl : ic;
				const size_t ak = std::abs(k);
				for(size_t ib = 0; ib < sl; ib++) {
					const int n = ib >= bw ? int(ib) - sl : ib;
					const size_t an = std::abs(n);
					const size_t maxKN = std::max(ak, an);
					for(size_t ia = 0; ia < dm; ia++) {
						const int m = ia >= bw ? int(ia) - sl : ia;
						const size_t am = std::abs(m);
						const size_t start = std::max<size_t>(am, maxKN);
						for(size_t j = start; j < bw; j++) {
							const Real dlkm = wigD[ak * bw * bw + am * bw + j] * wigner::dSign(j, k, m);//wigner::d<Real>(j, k, m);//this is for trans=false
							const Real dlnk = wigD[an * bw * bw + ak * bw + j] * wigner::dSign(j, n, k);//wigner::d<Real>(j, n, k);//this is for trans=false
							const std::complex<Real>& vflm = flm[am * bw + j];//\hat{f}^l_{|m|}
							const std::complex<Real>& vgln = gln[an * bw + j];//\hat{g}^l_{|n|}
							const std::complex<Real> f = std::signbit(m) ? std::conj(vflm) * Real(0 == am % 2 ? 1 : -1) : vflm;//symmetry of real SHT coefficients
							const std::complex<Real> g = std::signbit(n) ? std::conj(vgln) * Real(0 == an % 2 ? 1 : -1) : vgln;//symmetry of real SHT coefficients
							fxc[ic * sl * dm + ib * dm + ia] += f * std::conj(g) * dlkm * dlnk;
						}
					}
				}
			}
			*/

			//the above loop is conceptually simple but has many redundant calculations
			//the loop that follows is mathematically equivalent (for SHT of real functions only!!) much faster:
			// -use \hat{f}^l_{-m} = (-1)^m * \hat{f}^l_{m} for real valued functions (and the same for \hat{g}^l_{-n})
			// -build the fft of a real valued cross correlation
			// -precompute values of \hat{f}^l_{m} * d^l_{k,m}(\frac{\pi}{2}) [stored in fm]
			// -precompute values of \hat{g}^l_{n} * d^l_{n,k}(\frac{\pi}{2}) [stored in gn]
			// -eliminate redundant calculations from f * g and f * conj(g)

			//save some useful values
			const std::complex<Real> cz(0);
			const size_t mBw     = bw          ;//maximum bandwidth to use in calculation (better be <= bw but could be smaller for slight speedup, effectively a top hat filter size)
			const size_t flmFold = fNf         ;//rotational symmetry of flm about z axis (flm[m*bw+j] == 0 if m % flmFold != 0)
			const size_t glnFold = 1           ;//rotational symmetry of gln about z axis (gln[n*bw+j] == 0 if n % glnFold != 0)
			const bool   fMir    = fMr         ;//true/false if (flm[m*bw+j] == 0 if (m+j) % 2 != 0)
			const bool   gMir    = false       ;//true/false if (gln[n*bw+j] == 0 if (n+j) % 2 != 0)
			const bool   mirror  = fMir || gMir;//is there at least 1 mirror

			//precompute locations of systemic zeros
			//the storage for this should be moved to the constructor
			std::vector<uint_fast8_t> mFoldN0(bwP);
			for(size_t m = 0; m < mBw; m++) mFoldN0[m] = (uint_fast8_t) (0 != m % flmFold);//there are systemic zeros from rotational symmetry in flm

			//one value for each column for n % 2 == 0 and n % 2 == 1
			//true if cross correlation at n, m is nonzero, false if it is a systemic zero
			//using an integer array is much faster than bools becuase of the vector<bool> specialization
			std::vector<uint_fast8_t> nonZero0(bwP, false);//initialize with false for zero pad columns
			std::vector<uint_fast8_t> nonZero1(bwP, false);//initialize with false for zero pad columns
			const bool bMirror = fMir && gMir;//do both functions have a mirror
			for(size_t m = 0; m < mBw; m++) {
				const bool match0 = (m + 0) % 2 == 0;//check if parity of m matches parity of n == 0
				const bool match1 = (m + 1) % 2 == 0;//check if parity of m matches parity of n == 1
				const bool mir0   = bMirror && !match0;//there are systemic zeros if both functions have a mirror plane and a parity mismatch
				const bool mir1   = bMirror && !match1;//there are systemic zeros if both functions have a mirror plane and a parity mismatch
				const bool mFold  = 0 != m % flmFold;//there are systemic zeros from rotational symmetry in flm
				nonZero0[m] = !(mFold || mir0);//for n%2 == 0
				nonZero1[m] = !(mFold || mir1);//for n%2 == 1
			}

			//loop over planes
			std::complex<Real> v, vnc, vp, vc;
			std::complex<Real>* pK = fxc.data()                  ;//pointer to slice        k
			std::complex<Real>* nK = fxc.data() + slP * slP * bwP;//pointer to slice (slP - k)
			std::complex<Real> *pKpN, *pKnN, *nKpN, *nKnN;//row pointers
			for(size_t k = 0; k < mBw; k++) {//loop over non-zero z slices
				//precompute flm values * wigner d function
				//I transposed the wigner table from the simple (commented) loop to prioritize access for the inner loop
				Real const * pWig = xcLut->wigD.data() + k * bw;
				std::complex<Real>       * pFm  = fm .data();
				std::complex<Real> const * pFlm = flm       ;
				for(size_t m = 0; m < mBw; m++) {
					for(size_t j = std::max<size_t>(m, k); j < mBw; j++) pFm[j] = pFlm[j] * pWig[j];//f^j_m * wigner::d<Real>(j, k, m)
					pWig += bw * bw;
					pFm  += bw     ;
					pFlm += bw     ;
				}

				//compute pointers to row starts
				pKpN = pK            ;
				nKpN = nK            ;
				pKnN = pK + slP * bwP;
				nKnN = nK + slP * bwP;
				const bool posK = k > 0;

				//loop over rows
				pWig = xcLut->wigD.data() + k * bw * bw;
				std::complex<Real> const * pGln = gln;
				for(size_t n = 0; n < bwP; n++) {
					const bool posN = n > 0;
					const bool posKN = posK && posN;
					if(0 == n % glnFold && n < mBw) {//we haven't reached 0 padded rows and gln values are non-zero, compute dot product
						//precompute gln values * wigner d function
						const size_t maxKN = std::max<size_t>(k, n);
						for(size_t j = maxKN; j < mBw; j++) gn[j] = std::conj(pGln[j]) * pWig[j];//\hat{g}^j_n * wigner::d<Real>(j, n, k)

						//loop over columns
						std::vector<uint_fast8_t> const & nonZero = n % 2 == 0 ? nonZero0 : nonZero1;
						for(size_t m = 0; m < mBw; m++) {
							const size_t indMatch = (m + n) % 2;
							if(nonZero[m]) {//checking for systemic zeros from double mirror and parity mismatch here makes subsequent logic easy
								//build a pair of values as dot product
								const size_t m2 = m % 2;
								v = vnc = std::complex<Real>(0);
								size_t start = std::max<size_t>(m, maxKN);//first valid j
								if(mirror) {//there is a single mirror or double mirror with parity matching
									if(fMir && (start + m) % 2 != 0) ++start;//if fm[start] == 0 skip to next value (first value is zero)
									if(gMir && (start + n) % 2 != 0) ++start;//if gn[start] == 0 skip to next value (first value is zero) [for double mirrors no change here since parities match]
									const bool toggle = (start + m) % 2 == 0;//we don't need to toggle since we're incrementing by 2 but we still may need to negate the result
									for(size_t j = start; j < mBw; j+=2) {
										//do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
										detail::conjMult(fm[m * bw + j], gn[j], vp, vc);
										v   += vp;//pF[j] *           gn[j]
										vnc += vc;//pF[j] * std::conj(gn[j])
									}
									if(!toggle) vnc = -vnc;
								} else {
									bool toggle = (start + m) % 2 == 0;
									for(size_t j = start; j < mBw; j++) {
										//do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
										detail::conjMult(fm[m * bw + j], gn[j], vp, vc);
										v   +=          vp      ;//   pF[j] *           gn[j]
										vnc += toggle ? vc : -vc;//+/-pF[j] * std::conj(gn[j])
										toggle = !toggle;
									}
								}
								if(!(k % 2 == 0)) vnc = -vnc;//correct for computing negative vnc depending on j/m parity

								//fill in symmetric values using symmetry from: wigner d function, sht of real signal, sht of real pattern
								const bool match = 0 == indMatch;
								if(posKN) {
									pKpN    [m] =  v  ;//fxc( k,  n, m)
									nKnN    [m] =  vnc;//fxc(-k, -n, m)
									if(match) {
										nKpN[m] =  v  ;//fxc(-k,  n, m)
										pKnN[m] =  vnc;//fxc( k, -n, m)
									} else {
										nKpN[m] = -v  ;//fxc(-k,  n, m)
										pKnN[m] = -vnc;//fxc( k, -n, m)
									}
								} else {
									pKpN                  [m] =  v  ;//fxc( k,  n, m)
									if(match) {
										if     (posK) nKpN[m] =  v  ;//fxc(-k,  n, m)
										else if(posN) pKnN[m] =  vnc;//fxc( k, -n, m)
									} else {
										if     (posK) nKpN[m] = -v  ;//fxc(-k,  n, m)
										else if(posN) pKnN[m] = -vnc;//fxc( k, -n, m)
									}
								}
							} else {//in systemic zero
								if(posKN) {
									pKpN[m] = cz;//fxc( k,  n, m)
									nKpN[m] = cz;//fxc(-k,  n, m)
									pKnN[m] = cz;//fxc( k, -n, m)
									nKnN[m] = cz;//fxc(-k, -n, m)
								} else {
									              pKpN[m] = cz;//fxc( k,  n, m)
									if     (posK) nKpN[m] = cz;//fxc(-k,  n, m)
									else if(posN) pKnN[m] = cz;//fxc( k, -n, m)
								}
							}
						}
					} else {//we're in a row of systemic zeros or zero padding
						//fill in systemic zeros from rotational symmetry in gln
						std::fill          (pKpN, pKpN + bwP, cz);//fxc( k,  n, m)
						if(posKN) std::fill(nKnN, nKnN + bwP, cz);//fxc(-k, -n, m)
						if(posK ) std::fill(nKpN, nKpN + bwP, cz);//fxc(-k,  n, m)
						if(posN ) std::fill(pKnN, pKnN + bwP, cz);//fxc( k, -n, m)
					}

					//increment row pointers
					pKpN += bwP;
					nKpN += bwP;
					pKnN -= bwP;
					nKnN -= bwP;
					pWig += bw ;
					pGln += bw ;
				}

				//increment slice pointers
				pK += slP * bwP;
				nK -= slP * bwP;
			}

			//fill in zero pad slices
			if(mBw < bwP) std::fill(pK, nK + slP * bwP, cz);

			//compute cross correlation via fft
			xcLut->plan.inverse(fxc.data(), pXc, flmFold);//this skips systemic zeros and only computes the first half of euler space
		}

		//@brief : find the maximum cross correlation grid point
		//@return: index of maximum cross correlation from previous call to correlate
		template <typename Real>
		size_t Correlator<Real>::findPeak() {
			// return std::distance(xc.cbegin(), std::max_element(xc.cbegin(), xc.cbegin() + slP * slP * bwP));

			//profiler says this is noticeably faster:
			size_t iMax = 0;
			Real   vMax = xc.front();
			for(size_t i = 0; i < xc.size(); i++) {
				if(xc[i] > vMax) {
					vMax = xc[i];
					iMax = i;
				}
			}
			return iMax;
		}

		//@brief    : compute the first and second derivatives of the cross correlation at a single rotation
		//@param flm: spherical harmonic coefficients for the first function
		//@param gln: spherical harmonic coefficients for the second function
		//@param eu : rotation to compute derivatives of cross correlation for as ZYZ euler angle
		//@param jac: location to write jacobian of cross correlation {d/(d eu[0]), d/(d eu[1]), d/(d eu[2])}
		//@param hes: location to write hessian (3x3 matrix as 9 component vector) of cross correlation hes_ij = d/(d eu[i]) * d/(d eu[j])
		//@param mBW: maximum bandwidth to use in calculation (must be <= bw)
		//@param fMr: true/false if there is/isn't a mirror plane in the first fuction
		//@param fNf: rotational symmetry about z axis in first function (1 for no rotational symmetry)
		//@param der: true/false to compute derivatives/only cross correlation
		//@return   : cross correlation for rotation eu
		template <typename Real>
		Real Correlator<Real>::derivatives(std::complex<Real> const * const flm, std::complex<Real> const * const gln, Real const * const eu, Real * const jac, Real * const hes, const size_t mBW, const bool fMr, const size_t fNf, const bool der) {
			//initialze terms with 0
			Real wrk[10] = {0};//correlation, jacobian, hessian as 00, 11, 22, 01, 12, 20

			//bring middle euler angle to [-pi,pi] for wigner calculations
			Real beta = std::fmod(eu[1], emsphinx::Constants<Real>::pi2);
			if(beta > emsphinx::Constants<Real>::pi)
				beta -= emsphinx::Constants<Real>::pi2;
			else if(beta < -emsphinx::Constants<Real>::pi)
				beta += emsphinx::Constants<Real>::pi2;

			//compute sin/cos of alpha/gamma once (multiple angles)
			const Real sA = std::sin(eu[0]);//sin(alpha)
			const Real cA = std::cos(eu[0]);//cos(alpha)
			const Real sG = std::sin(eu[2]);//sin(gamma)
			const Real cG = std::cos(eu[2]);//cos(gamma)

			//precompute some values for on the fly wigner (uppercase) D calculation
			const Real t = std::cos(beta);
			const bool deg = std::fabs(std::fabs(t) - Real(1)) < std::numeric_limits<Real>::epsilon();//is beta nearly +/- pi?
			const bool nB = std::signbit(beta);
			const Real csc = Real(1) / std::sqrt(Real(1) - t * t) * (nB ? -1 : 1);//csc(beta), cot(beta) is csc * t
			// wigner::dTable   (mBW, t, nB, dBeta.data());//compute wigner (lowercase) d(beta) once
			wigner::dTablePre(mBW, t, nB, dBeta.data(), xcLut->wigE.data(), xcLut->wigW.data(), xcLut->wigB.data());//compute wigner (lowercase) d(beta) once

			//build symmetry information
			const size_t flmFold = fNf;//rotational symmetry of flm about z axis (flm[m*bw+j] == 0 if m % flmFold != 0)
			const size_t glnFold = 1  ;//rotational symmetry of gln about z axis (gln[n*bw+j] == 0 if n % glnFold != 0)
			const bool fMir = fMr  ;//true/false if (flm[m*bw+j] == 0 if (m+j) % 2 != 0)
			const bool gMir = false;//true/false if (gln[n*bw+j] == 0 if (n+j) % 2 != 0)
			const bool  mirror = fMir || gMir;//is there at least 1 mirror
			const bool bMirror = fMir && gMir;//do both functions have a mirror
			const int  dJ = fMir ? 2 : 1;

			////////////////////////////////////////
			//        loop over one order         //
			////////////////////////////////////////
			Real uA[3] = {0, sA * 2, 1};//recursion coefficients for chebyshev polynomial U_n(sin(alpha))
			Real tA[3] = {0, cA    , 1};//recursion coefficients for chebyshev polynomial T_n(cos(alpha))
			const int iBW = (int)mBW;//cast to integer once to silence signed/unsigned comparsion warnings
			for(int m = 0; m < iBW; m++) {
				////////////////////////////////////////
				// efficiently compute exp(I m alpha) //
				////////////////////////////////////////
				//update chebyshev recursion and use to compute exp(I m alpha)
				if(m < 2) {//use seed values for chebyshev recursion
					uA[0] = m == 0 ? 0 : sA;//sin(alpha * m)
					tA[0] = m == 0 ? 1 : cA;//cos(alpha * m)
				} else {//use chebyshev recursion
					//compute chebyshev polynomials(m, alpha) and multiple angle sin/cos
					uA[0] = sA * uA[1] * 2 - uA[2];//U_m(x) = 2 * x * U_{m-1}(x) - U_{m-2}(x)
					tA[0] = cA * tA[1] * 2 - tA[2];//T_m(x) = 2 * x * T_{m-1}(x) - T_{m-2}(x)
					const Real sm = m % 2 == 0 ? uA[1] * cA * (((m/2)-1) % 2 == 0 ? 1 : -1) : (uA[0] - sA * uA[1]) * (((m-1)/2) % 2 == 0 ? 1 : -1);//cos(alpha * m)
					const Real cm = tA[0];//cos(alpha * m)

					//update recursion and store values of sin/cos(alpha * m)
					uA[2] = uA[1]; uA[1] = uA[0];//update recursion coefficients for chebyshev polynomial of the first kind
					tA[2] = tA[1]; tA[1] = tA[0];//update recursion coefficients for chebyshev polynomial of the second kind
					uA[0] = sm;//store multiple sin value for subsequent access
					tA[0] = cm;//store multiple cos value for subsequent access
				}
				const bool mFold0 = 0 != m % flmFold;//there are systemic zeros from rotational symmetry in flm
				if(mFold0) continue;//flm[m * bw + j] == 0 so there is nothing to accumulate (but we still needed to update the multiple angle recursion)
				const std::complex<Real> expAlpha(tA[0], uA[0]);//exp(I m alpha) = cos(m * alpha) + I sin(m * alpha)

				////////////////////////////////////////
				//       loop over other order        //
				////////////////////////////////////////
				Real uG[3] = {0, sG * 2, 1};//recursion coefficients for chebyshev polynomial U_n(sin(gamma))
				Real tG[3] = {0, cG    , 1};//recursion coefficients for chebyshev polynomial T_n(cos(gamma))

				for(int n = 0; n < iBW; n++) {
					////////////////////////////////////////
					// efficiently compute exp(I n gamma) //
					////////////////////////////////////////
					//update chebyshev recursion and use to compute exp(I n gamma)
					if(n < 2) {//use seed values for chebyshev recursion
						uG[0] = n == 0 ? 0 : sG;//sin(gamma * n)
						tG[0] = n == 0 ? 1 : cG;//cos(gamma * n)
					} else {//use chebyshev recursion
						//compute chebyshev polynomials(n, gamma) and multiple angle sin/cos
						uG[0] = sG * uG[1] * 2 - uG[2];//U_n(x) = 2 * x * U_{n-1}(x) - U_{n-2}(x)
						tG[0] = cG * tG[1] * 2 - tG[2];//T_n(x) = 2 * x * T_{n-1}(x) - T_{n-2}(x)
						const Real sn = n % 2 == 0 ? uG[1] * cG * (((n/2)-1) % 2 == 0 ? 1 : -1) : (uG[0] - sG * uG[1]) * (((n-1)/2) % 2 == 0 ? 1 : -1);//cos(gamma * n)
						const Real cn = tG[0];//cos(gamma * n)

						//update recursion and store values of sin/cos(gamma * n)
						uG[2] = uG[1]; uG[1] = uG[0];//update recursion coefficients for chebyshev polynomial of the first kind
						tG[2] = tG[1]; tG[1] = tG[0];//update recursion coefficients for chebyshev polynomial of the second kind
						uG[0] = sn;//store multiple sin value for subsequent access
						tG[0] = cn;//store multiple cos value for subsequent access
					}
					const bool nFold0 = 0 != n % glnFold;//there are systemic zeros from rotational symmetry in flm
					if(nFold0) continue;//gln[n * bw + j] == 0 so there is nothing to accumulate (but we still needed to update the multiple angle recursion)
					const std::complex<Real> expGamma(tG[0], uG[0]);//exp(I n gamma) = cos(n * gamma) + I sin(n * gamma)

					//handle the case of 2 mirror planes with different parity
					const bool match = (m + n) % 2 == 0;//check if parity of m matches parity of n
					const bool mir0 = bMirror && !match;//there are systemic zeros if both functions have a mirror plane and a parity mismatch
					if(mir0) continue;//flm * gln = 0 for any l (alternating between flm and gln being 0)

					////////////////////////////////////////
					//  compute degree independent terms  //
					////////////////////////////////////////
					//compute exp(I * +m * Alpha) * exp(I * +/-n * Gamma) for on the fly wigner (uppercase) D calculation
					std::complex<Real> agP, agN;//exp(I * +m * Alpha) * exp(I * +/-n * Gamma)
					detail::conjMult(expAlpha, expGamma, agP, agN);
					const Real sign = (n+m)%2 == 0 ? 1 : -1;
					const Real sn = Real(n%2 == 0 ? 1 : -1);
					agP *= sign;
					agN *= sign * sn;

					//get loop start
					int start = std::max<int>(m, n);
					if(fMir && (start + m) % 2 != 0) ++start;//if fm[start] == 0 skip to next value (first value is zero)
					if(gMir && (start + n) % 2 != 0) ++start;//if gn[start] == 0 skip to next value (first value is zero) [for double mirrors no change here since parities match]
					
					if(der) {//we need derivatives
						//compute some prefactors for calculating derivatives of d^j_{m,n}(beta)
						const int mm = m * m;
						const int mn = m * n;
						const int nn = n * n;
						const Real coef2_0a  =   t * t * mm + (nn - m)             ;
						const Real coef2_0b  =   t * n * (1 - 2 * m)               ;
						const Real coef2_1a  =   t *     (1 + 2 * m)               ;
						const Real coef1_0PP = ( t * m      - n       ) * csc      ;//this is infinity if degenerate
						const Real coef1_0PN = ( t * m      + n       ) * csc      ;//this is infinity if degenerate
						const Real coef2_0PP = (coef2_0a    + coef2_0b) * csc * csc;//this is infinity if degenerate
						const Real coef2_0PN = (coef2_0a    - coef2_0b) * csc * csc;//this is infinity if degenerate
						const Real coef2_1PP = (coef2_1a    - 2 * n   ) * csc      ;//this is infinity if degenerate
						const Real coef2_1PN = (coef2_1a    + 2 * n   ) * csc      ;//this is infinity if degenerate

						////////////////////////////////////////
						//   loop over degrees accumulating   //
						////////////////////////////////////////
						for(int j = start; j < iBW; j+=dJ) {//increment by 1 for no mirror planes 2 if any are present
							//get wigner d^j_{m,+/-n} components
							const Real d0P    =                dBeta[((m  ) * mBW * mBW + n * mBW + j) * 2 + 0];//d^j_{m  ,n}(     beta)
							const Real d0N    =                dBeta[((m  ) * mBW * mBW + n * mBW + j) * 2 + 1];//d^j_{m  ,n}(pi - beta)
							const Real d0P_1  = m   >= j ? 0 : dBeta[((m+1) * mBW * mBW + n * mBW + j) * 2 + 0];//d^j_{m+1,n}(     beta)
							const Real d0N_1  = m   >= j ? 0 : dBeta[((m+1) * mBW * mBW + n * mBW + j) * 2 + 1];//d^j_{m+1,n}(pi - beta)
							const Real d0P_2  = m+1 >= j ? 0 : dBeta[((m+2) * mBW * mBW + n * mBW + j) * 2 + 0];//d^j_{m+2,n}(     beta)
							const Real d0N_2  = m+1 >= j ? 0 : dBeta[((m+2) * mBW * mBW + n * mBW + j) * 2 + 1];//d^j_{m+2,n}(pi - beta)

							//compute derivatives of d^j_{m,+/-n}(beta)
							const int  jm       = j - m;
							const Real rjm      =               std::sqrt( Real( (jm    ) * (j + m + 1) ) );
							const Real coef2_2  = 0 == jm ? 0 : std::sqrt( Real( (jm - 1) * (j + m + 2) ) ) * rjm;
							const Real d1P = d0P * coef1_0PP - d0P_1 * rjm                              ;//first  derivative of d^j_{+m,+n}(beta) w.r.t. beta
							const Real d1N = d0N * coef1_0PN + d0N_1 * rjm                              ;//first  derivative of d^j_{+m,-n}(beta) w.r.t. beta
							const Real d2P = d0P * coef2_0PP - d0P_1 * rjm * coef2_1PP + d0P_2 * coef2_2;//second derivative of d^j_{+m,+n}(beta) w.r.t. beta
							const Real d2N = d0N * coef2_0PN + d0N_1 * rjm * coef2_1PN + d0N_2 * coef2_2;//second derivative of d^j_{+m,-n}(beta) w.r.t. beta

							//compute f^l_m * g^l_n
							std::complex<Real> vp, vc;//\hat{f}^l_{+m} * hat{g}^l_{+n} and \hat{f}^l_{+m} * conj(hat{g}^l_{+n})
							detail::conjMult(flm[m * bw + j], gln[n * bw + j], vp, vc);//do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
							if((j+m)%2 != 0) vp = -vp;

							//compute components of cross correlation and partials
							const std::complex<Real> vcPP  = vc   * agP;//+n correlation term prefactor: \hat{f}^l_{+m} * conj(hat{g}^l_{+n}) * exp(I m alpha + I n gamma)
							const std::complex<Real> vcPP0 = vcPP * d0P;//+n correlation term: \hat{f}^l_{+m} * conj(hat{g}^l_{+n}) * D^l_{+m,+n}(alpha, beta, gamma)
							const std::complex<Real> vcPP1 = vcPP * d1P;//beta partial of +n correlation term
							const std::complex<Real> vpPN  = vp   * agN;//+n correlation term prefactor: \hat{f}^l_{+m} * conj(hat{g}^l_{-n}) * exp(I m alpha - I n gamma)
							const std::complex<Real> vpPN0 = vpPN * d0N;//-n correlation term: \hat{f}^l_{+m} * conj(hat{g}^l_{-n}) * D^l_{+m,-n}(alpha, beta, gamma)
							const std::complex<Real> vpPN1 = vpPN * d1N;//beta partial of -n correlation term

							//compute contributions to cross correlation, jacobian, and hessian from +m,+n (the real part of contributions from -m,-n are the same for real functions)
							const Real xc[10] = {
								vcPP0.real(),                                              //cross correlation
								vcPP0.imag() * -m , vcPP1.real()      , vcPP0.imag() * -n ,//alpha  , beta  , gamma   derivatives
								vcPP0.real() * -mm, vcPP .real() * d2P, vcPP0.real() * -nn,//alpha^2, beta^2, gamma^2 derivatives
								vcPP1.imag() * -m , vcPP1.imag() * -n , vcPP0.real() * -mn //alpha beta, beta gamma, gamma alpha derivatives
							};

							//compute contributions to cross correlation, jacobian, and hessian from +m,-n (the real part of contributions from -m,+n are the same for real functions)
							const Real xp[10] = {
								vpPN0.real(),                                              //cross correlation
								vpPN0.imag() * -m , vpPN1.real()      , vpPN0.imag() *  n ,//alpha  , beta  , gamma   derivatives
								vpPN0.real() * -mm, vpPN .real() * d2N, vpPN0.real() * -nn,//alpha^2, beta^2, gamma^2 derivatives
								vpPN1.imag() * -m , vpPN1.imag() *  n , vpPN0.real() *  mn //alpha beta, beta gamma, gamma alpha derivatives
							};

							//accumulate contributions
							std::transform(wrk, wrk + 10, xc, wrk, std::plus<Real>());//+m,+n
							if(n > 0) std::transform(wrk, wrk + 10, xp, wrk, std::plus<Real>());//+m,-n
							if(m > 0) {
								std::transform(wrk, wrk + 10, xp, wrk, std::plus<Real>());//-m,+n
								if(n > 0) std::transform(wrk, wrk + 10, xc, wrk, std::plus<Real>());//-m,-n
							}
						}
					} else {//we only need cross correlation
						for(int j = start; j < iBW; j+=dJ) {//increment by 1 for no mirror planes 2 if any are present
							//get wigner d^j_{m,+/-n} components
							const Real& d0P = dBeta[((m  ) * mBW * mBW + n * mBW + j) * 2 + 0];//d^j_{m  ,n}(     beta)
							const Real& d0N = dBeta[((m  ) * mBW * mBW + n * mBW + j) * 2 + 1];//d^j_{m  ,n}(pi - beta)

							//compute f^l_m * g^l_n
							std::complex<Real> vp, vc;//\hat{f}^l_{+m} * hat{g}^l_{+n} and \hat{f}^l_{+m} * conj(hat{g}^l_{+n})
							detail::conjMult(flm[m * bw + j], gln[n * bw + j], vp, vc);//do complex multiplication components by hand to eliminate duplicate flops from multiplying with conjugate
							if((j+m)%2 != 0) vp = -vp;

							//compute contributions to cross correlation
							const Real vcPP0 = vc.real() * agP.real() - vc.imag() * agP.imag();//(vc * agP).real()
							const Real vpPN0 = vp.real() * agN.real() - vp.imag() * agN.imag();//(vp * agN).real()
							const Real xc = vcPP0 * d0P;//from +m,+n
							const Real xp = vpPN0 * d0N;//from +m,-n

							//accumulate contributions
							wrk[0] += xc;//+m,+n
							if(n > 0) wrk[0] += xp;//+m,-n
							if(m > 0) {
								wrk[0] += xp;//-m,+n
								if(n > 0) wrk[0] += xc;//-m,-n
							}
						}
					}
				}
			}

			//copy result to outputs and return correlation
			if(der) {
				std::copy(wrk + 1, wrk + 4, jac);
				hes[0] = wrk[4+0]; hes[1] = wrk[4+3]; hes[2] = wrk[4+5];
				                   hes[4] = wrk[4+1]; hes[5] = wrk[4+4];
				                                      hes[8] = wrk[4+2];
				hes[3] = hes[1]; 
				hes[6] = hes[2]; hes[7] = hes[5]; 
			}
			return wrk[0];
		}

		//@brief          : construct a spherical correlator for a given bandwidth
		//@param bandWidth: maximum bandwidth of spherical harmonic to use in correlation (exclusive)
		//@param flm      : spherical harmonic coefficients for the reference function
		//@param flm2     : spherical harmonic coefficients for the reference (function^2)
		//@param fMr      : true/false if there is/isn't a mirror plane in the reference fuction
		//@param fNf      : rotational symmetry about z axis in reference function (1 for no rotational symmetry)
		//@param mlm      : spherical harmonic coefficients for the mask function
		template <typename Real>
		NormalizedCorrelator<Real>::NormalizedCorrelator(const size_t bandWidth, std::complex<Real> const * const flm, std::complex<Real> const * const flm2, const bool fMr, const size_t fNf, std::complex<Real> const * const mlm) :
		PhaseCorrelator<Real>::PhaseCorrelator(bandWidth),
		ncLut(std::make_shared<const Constants>(this, flm, flm2, fMr, fNf, mlm)) {}

		//@brief    : compute the rotation of the maximum normalized cross correlation between two spherical functions
		//@param gln: spherical harmonic coefficients for the template function
		//@param eu : location to write rotation of maximum normalized cross correlation as ZYZ euler angle
		//@param ref: true/false to use/not use real space refinement
		//@param eps: convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
		//@return   : maximum (semi)normalized cross correlation [still needs to be divided by the standard deviation of the pattern function]
		template <typename Real>
		Real NormalizedCorrelator<Real>::correlate(std::complex<Real> const * const gln, Real * const eu, const bool ref, const Real eps) {
			//compute cross correlation on euler grid
			Correlator<Real>::compute(ncLut->flm.data(), gln, ncLut->mr, ncLut->nf, Correlator<Real>::xc.data());

			//normalize and find peak in single pass
			size_t iMax = 0;
			Real   vMax = Correlator<Real>::xc.front() * ncLut->rDen.front();
			for(size_t i = 0; i < Correlator<Real>::xc.size(); i++) {
				Real& vI = Correlator<Real>::xc[i];//get value
				vI *= ncLut->rDen[i];//normalize
				if(vI > vMax) {//check for new max
					vMax = vI;
					iMax = i;
				}
			}

			//refine peak
			const Real peak = Correlator<Real>::interpPeak(iMax, eu);
			return ref ? NormalizedCorrelator<Real>::refinePeak(gln, eu, eps) : peak;
		}

		//@brief     : compute a subpixel maxima in the normalized cross correlation by either interpolation or refinement
		//@param gln : spherical harmonic coefficients for the second function
		//@param eu  : location to write rotation of maximum cross correlation as ZYZ euler angle
		//@param eps : convergence criterion for refinement (step size in fractional pixels where, absolute precision is ~eps * pi / bw radians)
		//@return    : cross correlation of peak
		//@note      : to be fully robust this should use the chain rule to account for the effects of shifting the window (which it doesn't currently do)
		//             not using the chain rule assumes that the effect of the window is small near the peak
		template <typename Real>
		Real NormalizedCorrelator<Real>::refinePeak(std::complex<Real> const * const gln, Real * const eu, const Real eps) {
			const Real cor = Correlator<Real>::refinePeak(ncLut->flm.data(), gln, ncLut->mr, ncLut->nf, eu, eps);//find unnormalized peak
			return cor / ncLut->denominator(*this, eu);//normalize peak
		}

		//@brief     : construct a spherical correlator for a given bandwidth
		//@param cor : correlator object
		//@param flm : spherical harmonic coefficients for the reference function
		//@param flm2: spherical harmonic coefficients for the reference (function^2)
		//@param fMr : true/false if there is/isn't a mirror plane in the reference fuction
		//@param fNf : rotational symmetry about z axis in reference function (1 for no rotational symmetry)
		//@param mlm : spherical harmonic coefficients for the mask function
		template <typename Real>
		NormalizedCorrelator<Real>::Constants::Constants(NormalizedCorrelator<Real>* cor, std::complex<Real> const * const flm, std::complex<Real> const * const flm2, const bool fMr, const size_t fNf, std::complex<Real> const * const mlm) :
			rDen(cor->getCubeSize() * cor->getCubeSize() * cor->getHalfSize()),//allocate space for denominator
			flm (flm , flm  + cor->bw * cor->bw                              ),//save reference spectra
			flm2(flm2, flm2 + cor->bw * cor->bw                              ),//save reference^2 spectra
			mlm (mlm , mlm  + cor->bw * cor->bw                              ),//save mask spectra
			mr  (fMr                                                         ),//save reference mirror
			nf  (fNf                                                         ) //save reference rotational symmetry
		{
			//first compute window function correlated with reference function
			cor->compute(flm, mlm, mr, nf, rDen.data());//store result in rDen array

			//next compute window function correlated with reference function^2
			cor->compute(flm2, mlm, mr, nf, cor->xc.data());//store result in working array of cor

			//finally compute the integral of the window function
			const Real s2m = mlm[0].real() * std::sqrt(emsphinx::Constants<Real>::pi * Real(4));//assumes the window function is a binary mask [0,4pi]

			//now we have all the parts needed to fill in the first part of the denominator (equation 8)
			std::transform(rDen.cbegin(), rDen.cend(), cor->xc.cbegin(), rDen.begin(), [s2m](const Real& mrf, const Real& mrf2){
				const Real fWbar = mrf / s2m;//equation 9
				return Real(1) / std::sqrt(mrf2 - Real(2) * fWbar * mrf + fWbar * fWbar * s2m);//compute reciprocal once to convert subsequent divisions into multiplications
			});
		}

		//@brief    : compute the normalization denominator at an arbitrary rotation
		//@param cor: correlator object
		//@param eu : rotation to compute normalization for
		//@return   : normalization
		template <typename Real>
		Real NormalizedCorrelator<Real>::Constants::denominator(Correlator<Real>& cor, Real * const eu) const {
			//first compute window function correlated with reference function
			const Real mrf = cor.derivatives(flm.data(), mlm.data(), eu, NULL, NULL, cor.getBw(), mr, nf, false);

			//next compute window function correlated with reference function^2
			const Real mrf2 = cor.derivatives(flm2.data(), mlm.data(), eu, NULL, NULL, cor.getBw(), mr, nf, false);

			//finally compute the integral of the window function
			// const Real s2m = mlm[0].real() / std::sqrt(emsphinx::Constants<Real>::pi * Real(4));//assumes the window function is a binary mask, {fractional from [0,1]}
			const Real s2m = mlm[0].real() * std::sqrt(emsphinx::Constants<Real>::pi * Real(4));//assumes the window function is a binary mask [0,4pi]

			//now we have all the parts needed to compute the denominator (equation 8)
			const Real fWbar = mrf / s2m;
			return std::sqrt(mrf2 - Real(2) * fWbar * mrf + fWbar * fWbar * s2m);
		}

		namespace detail {
			//@brief   : compute product of ab * cd and ad * conj(cd) without duplicate flops
			//@param ab: first  complex number (a + bi)
			//@param cd: second complex number (c + di)
			//@param vp: location to write ab *      cd
			//@param vc: location to write ab * conj(cd)
			template <typename T>
			inline void conjMult(const std::complex<T>& ab, const std::complex<T>& cd, std::complex<T>& vp, std::complex<T>& vc) {
				T rr = ab.real() * cd.real();//ac (real * real)
				T ii = ab.imag() * cd.imag();//bd (imag * imag)
				T ri = ab.real() * cd.imag();//ad (real * imag)
				T ir = ab.imag() * cd.real();//bc (imag * real)
				vp.real(rr - ii);
				vc.real(rr + ii);
				vp.imag(ir + ri);
				vc.imag(ir - ri);
			}

			//@brief    : convert a vectorized index in the 3d cross correlation grid to individual components
			//@param idx: vectorized index
			//@param sl : side length of cube
			//@param knm: location to write indices (z,y,x)
			inline void extractInds(size_t idx, const size_t sl, size_t knm[3]) {
				knm[0] = idx / (sl * sl);//k component of idx
				idx -= knm[0] * sl * sl;
				knm[1] = idx / sl;//n component of idx
				idx -= knm[1] * sl;
				knm[2] = idx;//m component of idx
			}

			//@brief: interpolate subpixel peak location from a 3d voxel grid
			//@param p: neighborhood around peak
			//@param x: location to store subpixel maxima location within neighborhood (z, y, x from -1->1)
			//@param return: value of fit quadratic at maxima
			template <typename Real, typename PixelType>
			Real interpolateMaxima(PixelType p[3][3][3], Real x[3]) {
				//compute the 27 biquadradic coefficients, f(x,y,z) = a_{kji} x^i y^j z^k
				//f(0, 0, 0) == a000
				const Real a000 = Real(p[1][1][1]);

				//f(1,0,0) = a000 + a100 + a200 && f(-1,0,0) = a000 - a100 + a200
				const Real a001 = Real(p[1][1][2] - p[1][1][0]) / 2;
				const Real a002 = Real(p[1][1][2] + p[1][1][0]) / 2 - a000;

				//same relationships for y and z
				const Real a010 = Real(p[1][2][1] - p[1][0][1]) / 2;
				const Real a020 = Real(p[1][2][1] + p[1][0][1]) / 2 - a000;
				const Real a100 = Real(p[2][1][1] - p[0][1][1]) / 2;
				const Real a200 = Real(p[2][1][1] + p[0][1][1]) / 2 - a000;
				
				//f( 1, 1,0) = a000 + a100 + a200 + a010 + a020 + a110 + a210 + a120 + a220
				//f( 1,-1,0) = a000 + a100 + a200 - a010 + a020 - a110 - a210 + a120 + a220
				//f(-1, 1,0) = a000 - a100 + a200 + a010 + a020 - a110 + a210 - a120 + a220
				//f(-1,-1,0) = a000 - a100 + a200 - a010 + a020 + a110 - a210 - a120 + a220
				// --> f( 1, 1,0) + f( 1,-1,0) + f(-1, 1,0) + f(-1,-1,0) = 4 * (a000 + a020 + a200 + a220)
				// --> f( 1, 1,0) - f( 1,-1,0) - f(-1, 1,0) + f(-1,-1,0) = 4 * a110
				// --> f( 1, 1,0) - f( 1,-1,0) + f(-1, 1,0) - f(-1,-1,0) = 4 * (a100 + a120)
				// --> f( 1, 1,0) + f( 1,-1,0) - f(-1, 1,0) - f(-1,-1,0) = 4 * (a010 + a210)
				const Real a022 = Real(p[1][2][2] + p[1][2][0] + p[1][0][2] + p[1][0][0]) / 4 - a000 - a020 - a002;
				const Real a011 = Real(p[1][2][2] - p[1][2][0] - p[1][0][2] + p[1][0][0]) / 4;
				const Real a012 = Real(p[1][2][2] + p[1][2][0] - p[1][0][2] - p[1][0][0]) / 4 - a010;
				const Real a021 = Real(p[1][2][2] - p[1][2][0] + p[1][0][2] - p[1][0][0]) / 4 - a001;

				//same relationships for yz and zx
				const Real a220 = Real(p[2][2][1] + p[2][0][1] + p[0][2][1] + p[0][0][1]) / 4 - a000 - a200 - a020;
				const Real a110 = Real(p[2][2][1] - p[2][0][1] - p[0][2][1] + p[0][0][1]) / 4;
				const Real a120 = Real(p[2][2][1] + p[2][0][1] - p[0][2][1] - p[0][0][1]) / 4 - a100;
				const Real a210 = Real(p[2][2][1] - p[2][0][1] + p[0][2][1] - p[0][0][1]) / 4 - a010;
				const Real a202 = Real(p[2][1][2] + p[0][1][2] + p[2][1][0] + p[0][1][0]) / 4 - a000 - a002 - a200;
				const Real a101 = Real(p[2][1][2] - p[0][1][2] - p[2][1][0] + p[0][1][0]) / 4;
				const Real a201 = Real(p[2][1][2] + p[0][1][2] - p[2][1][0] - p[0][1][0]) / 4 - a001;
				const Real a102 = Real(p[2][1][2] - p[0][1][2] + p[2][1][0] - p[0][1][0]) / 4 - a100;

				//similar relationships for corners
				const Real a222 = Real(p[2][2][2] + p[0][0][0] + p[0][2][2] + p[2][0][2] + p[2][2][0] + p[2][0][0] + p[0][2][0] + p[0][0][2]) / 8 - a000 - a200 - a020 - a002 - a022 - a202 - a220;
				const Real a211 = Real(p[2][2][2] + p[0][0][0] + p[0][2][2] - p[2][0][2] - p[2][2][0] + p[2][0][0] - p[0][2][0] - p[0][0][2]) / 8 - a011;
				const Real a121 = Real(p[2][2][2] + p[0][0][0] - p[0][2][2] + p[2][0][2] - p[2][2][0] - p[2][0][0] + p[0][2][0] - p[0][0][2]) / 8 - a101;
				const Real a112 = Real(p[2][2][2] + p[0][0][0] - p[0][2][2] - p[2][0][2] + p[2][2][0] - p[2][0][0] - p[0][2][0] + p[0][0][2]) / 8 - a110;
				const Real a111 = Real(p[2][2][2] - p[0][0][0] - p[0][2][2] - p[2][0][2] - p[2][2][0] + p[2][0][0] + p[0][2][0] + p[0][0][2]) / 8;
				const Real a122 = Real(p[2][2][2] - p[0][0][0] - p[0][2][2] + p[2][0][2] + p[2][2][0] + p[2][0][0] - p[0][2][0] - p[0][0][2]) / 8 - a100 - a120 - a102;
				const Real a212 = Real(p[2][2][2] - p[0][0][0] + p[0][2][2] - p[2][0][2] + p[2][2][0] - p[2][0][0] + p[0][2][0] - p[0][0][2]) / 8 - a010 - a012 - a210;
				const Real a221 = Real(p[2][2][2] - p[0][0][0] + p[0][2][2] + p[2][0][2] - p[2][2][0] - p[2][0][0] - p[0][2][0] + p[0][0][2]) / 8 - a001 - a201 - a021;

				//newton iterate to find maxima
				x[0] = x[1] = x[2] = 0;//initial guess at maximum voxel (z,y,x)
				const size_t maxIter = 25;
				const Real eps = std::sqrt(std::numeric_limits<Real>::epsilon());
				for(size_t i = 0; i < maxIter; i++) {
					//compute components of hessian matrix
					const Real xx = x[0] * x[0]; const Real yy = x[1] * x[1]; const Real zz = x[2] * x[2];
					const Real xy = x[0] * x[1]; const Real yz = x[1] * x[2]; const Real zx = x[2] * x[0];
					const Real h00 = (a200 + a210 * x[1] + a201 * x[2] + a220 * yy + a202 * zz + a211 * yz + a221 * yy * x[2] + a212 * x[1] * zz + a222 * yy * zz) * 2;
					const Real h11 = (a020 + a021 * x[2] + a120 * x[0] + a022 * zz + a220 * xx + a121 * zx + a122 * zz * x[0] + a221 * x[2] * xx + a222 * zz * xx) * 2;
					const Real h22 = (a002 + a102 * x[0] + a012 * x[1] + a202 * xx + a022 * yy + a112 * xy + a212 * xx * x[1] + a122 * x[0] * yy + a222 * xx * yy) * 2;
					const Real h01 = a110 + a111 * x[2] + a112 * zz + (a210 * x[0] + a120 * x[1] + a211 * zx + a121 * yz + a212 * x[0] * zz + a122 * x[1] * zz + (a220 * xy + a221 * xy * x[2] + a222 * xy * zz) * 2) * 2;
					const Real h12 = a011 + a111 * x[0] + a211 * xx + (a021 * x[1] + a012 * x[2] + a121 * xy + a112 * zx + a221 * x[1] * xx + a212 * x[2] * xx + (a022 * yz + a122 * yz * x[0] + a222 * yz * xx) * 2) * 2;
					const Real h02 = a101 + a111 * x[1] + a121 * yy + (a102 * x[2] + a201 * x[0] + a112 * yz + a211 * xy + a122 * x[2] * yy + a221 * x[0] * yy + (a202 * zx + a212 * zx * x[1] + a222 * zx * yy) * 2) * 2;
					
					//build inverse of hessian matrix
					const Real det = h00 * h11 * h22 - h00 * h12 * h12 - h11 * h02 * h02 - h22 * h01 * h01 + h01 * h12 * h02 * 2;
					const Real i00 = (h11 * h22 - h12 * h12) / det;
					const Real i11 = (h22 * h00 - h02 * h02) / det;
					const Real i22 = (h00 * h11 - h01 * h01) / det;
					const Real i01 = (h02 * h12 - h01 * h22) / det;
					const Real i12 = (h01 * h02 - h12 * h00) / det;
					const Real i02 = (h12 * h01 - h02 * h11) / det;

					//compute gradient
					const Real d0 = a100 + a110 * x[1] + a101 * x[2] + a120 * yy + a102 * zz + a111 * yz + a121 * yy * x[2] + a112 * x[1] * zz + a122 * yy * zz + x[0] * (a200 + a210 * x[1] + a201 * x[2] + a220 * yy + a202 * zz + a211 * yz + a221 * yy * x[2] + a212 * x[1] * zz + a222 * yy * zz) * 2;
					const Real d1 = a010 + a011 * x[2] + a110 * x[0] + a012 * zz + a210 * xx + a111 * zx + a112 * zz * x[0] + a211 * x[2] * xx + a212 * zz * xx + x[1] * (a020 + a021 * x[2] + a120 * x[0] + a022 * zz + a220 * xx + a121 * zx + a122 * zz * x[0] + a221 * x[2] * xx + a222 * zz * xx) * 2;
					const Real d2 = a001 + a101 * x[0] + a011 * x[1] + a201 * xx + a021 * yy + a111 * xy + a211 * xx * x[1] + a121 * x[0] * yy + a221 * xx * yy + x[2] * (a002 + a102 * x[0] + a012 * x[1] + a202 * xx + a022 * yy + a112 * xy + a212 * xx * x[1] + a122 * x[0] * yy + a222 * xx * yy) * 2;

					//update x
					const Real step[3] = {
						i00 * d0 + i01 * d1 + i02 * d2,
						i01 * d0 + i11 * d1 + i12 * d2,
						i02 * d0 + i12 * d1 + i22 * d2
					};
					x[0] -= step[0]; x[1] -= step[1]; x[2] -= step[2];

					//check for convergence
					const Real maxStep = std::max(std::max(std::fabs(step[0]), std::fabs(step[1])), std::fabs(step[2]));
					if(maxStep < eps) break;
					if(i+1 == maxIter) std::fill(x, x+3, Real(0));//don't interpolate if convergence wasn't reached
				}

				//compute interpolated value of maxima
				const Real xx = x[0] * x[0]; const Real yy = x[1] * x[1]; const Real zz = x[2] * x[2];
				const Real xy = x[0] * x[1]; const Real yz = x[1] * x[2]; const Real zx = x[2] * x[0];
				const Real vPeak = a000                    + a111 * x[0] * x[1] * x[2] + a222 * xx   * yy   * zz
				                 + a100 * x[0]             + a010 * x[1]               + a001 * x[2]
				                 + a200 * xx               + a020 * yy                 + a002 * zz
				                 + a110 * xy               + a011 * yz                 + a101 * zx
				                 + a120 * x[0] * yy        + a012 * x[1] * zz          + a201 * x[2] * xx
				                 + a210 * xx   * x[1]      + a021 * yy   * x[2]        + a102 * zz   * x[0]
				                 + a220 * xx   * yy        + a022 * yy   * zz          + a202 * zz   * xx
				                 + a112 * xy   * x[2]      + a211 * yz   * x[0]        + a121 * zx   * x[1]
				                 + a122 * x[0] * yy   * zz + a212 * xx   * x[1] * zz   + a221 * xx   * yy   * x[2];
				return vPeak;
			}
		}
	}

}

#endif//_SHT_XCORR_H_
