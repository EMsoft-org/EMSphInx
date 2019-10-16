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

#ifndef _SPHERE_INDEXER_H_
#define _SPHERE_INDEXER_H_

#include <memory>
#include <type_traits>
#include <vector>
#include <array>

#include "idx/base.hpp"
#include "idx/master.hpp"
#include "sht/square_sht.hpp"
#include "sht/sht_xcorr.hpp"
#include "xtal/quaternion.hpp"


namespace emsphinx {

	//@brief: a single indexing result
	template <typename Real>
	struct Result {
		Real corr ;//maximum cross correlation
		Real iq   ;//image quality metric
		int  phase;//phase of maximum cross correlation
		Real qu[4];//orientation of maximum cross correlation

		//@brief    : comparison operator
		//@param rhs: result to compare against
		//@return   : this < rhs (such that sort is in descending corr)
		bool operator<(const Result& rhs) const {return corr > rhs.corr;}//sort from largest to smallest cross correlation
	};

	//@brief: abstract indexer base class
	template <typename Real>	
	struct Indexer {
		//storage values values
		const size_t                                                    mBw ;//bandwidth to use for indexing
		const size_t                                                    dim ;//side length of square legendre grid for spherical harmonic transformer
		const std::array<Real, 4>                                       quNp;//quaternion to rotate detector such that is centered around the north pole (correction in orientation is then required)
		Real                                                            sum2;//sumsq of the back projected image in sph (integral over the window of sph^2 * dOmega)
		std::vector <std::vector< xtal::Quat<Real> > >                  pSym;//psuedosymmetric orientations to check for each phase	
		
		//larger work spaces
		std::vector<Real>                                               wrk ;//space to hold processed image
		std::vector<Real>                                               sph ;//space to hold image back projected onto sphere
		std::vector< std::complex<Real> >                               gln ;//space to hold SHT of image begin indexed

		//calculators
		std::unique_ptr< ImageProcessor<Real> >                         prc ;//image processor
		std::unique_ptr< BackProjector <Real> >                         prj ;//back projector
		square::DiscreteSHT<Real>                                       sht ;//spherical harmonic transform calculator
		std::vector< std::unique_ptr< sphere::PhaseCorrelator<Real> > > xc  ;//spherical correlator (for each phase)

		//@brief        : construct an indexer
		//@param bw     : bandwidth
		//@param imPrc  : image processor
		//@param backPrj: back projector
		//@param corr   : cross correlator
		Indexer(const size_t bw, std::unique_ptr< ImageProcessor<Real> > imPrc, std::unique_ptr< BackProjector<Real> > backPrj, const std::vector< std::unique_ptr< sphere::PhaseCorrelator<Real> > >& corrs);

		//@brief   : estimate a reasonable batch size
		//@param bw: bandwidth
		//@param nt: number of worker threads
		//@param np: number of images to index
		//@return  : a reasonable batch size (a number of images that should take on the order of 1s to index
		static size_t BatchEstimate(const size_t bw, const size_t nt, const size_t np);

		//@brief       : index a single image
		//@param    pat: pointer to image to index in row major order
		//@param    n  : number of results to keep (top n), user is responsible for making sure this is reasonable (extra points will be filled with an invalid phase)
		//@param    res: location to write results
		//@param    ref: true/false to use newton's method refinement/triquadratic interpolation
		//@template Pix: pixel type of image (will be cast to Real)
		//@note        : multiple results can be written (currently there is an additional result per psuedosymmetric orientation to check)
		template <typename Pix>
		void indexImage(Pix* pat, Result<Real>*const res, const size_t n, const bool ref);

		//@brief       : refine the orientation of a single image
		//@param    pat: pointer to image to refine in row major order
		//@param    res: location to write result (and read initial phase + orientation)
		//@template Pix: pixel type of image (will be cast to Real)
		template <typename Pix>
		void refineImage(Pix* pat, Result<Real>& res);

		//@brief : get a copy of the stored image processor
		//@return: unique pointer to copy of current processor
		std::unique_ptr<Indexer> clone() const {return std::unique_ptr<Indexer>(new Indexer(mBw, std::move(prc->clone()), std::move(prj->clone()), xc));}

		protected:

			//@brief   : process an image, unproject to the sphere, and compute spherical harmonic coefficients
			//@param im: image to process
			//@return  : image quality
			template <typename Pix>
			Real computeHarmonics(Pix* im);

			//@brief    : index the spherical function currently held in sph
			//@param p  : phase to correlate with
			//@param gln: spectra of image to index
			//@param ref: true/false to use newton's method refinement/triquadratic interpolation
			//@return   : correlation result (with ZYZ euler angle stored in qu[0,1,2])
			Result<Real> correlate(const size_t p, std::complex<Real> const * gln, const bool ref);

			//@brief    : refine the spherical function currently held in sph
			//@param p  : phase to refine with
			//@param gln: spectra of image to refine
			//@param eu : initial ZYZ euler angle to refine from
			//@return   : refinement result (with ZYZ euler angle stored in qu[0,1,2])
			Result<Real> refine(const size_t p, std::complex<Real> const * gln, Real const * const eu);
	};


}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

namespace emsphinx {

	////////////////////////////////////////////////////////////////////////
	//                              Indexer                               //
	////////////////////////////////////////////////////////////////////////

	//@brief        : construct an indexer
	//@param bw     : bandwidth
	//@param imPrc  : image processor
	//@param backPrj: back projector
	//@param corr   : cross correlator
	template <typename Real>
	Indexer<Real>::Indexer(const size_t bw, std::unique_ptr< ImageProcessor<Real> > imPrc, std::unique_ptr< BackProjector<Real> > backPrj, const std::vector< std::unique_ptr< sphere::PhaseCorrelator<Real> > >& corrs) :
		mBw (bw                                       ),//save bandwidth
		dim (mBw + (mBw % 2 == 0 ? 3 : 2)             ),//compute smallest odd dimension for bandwidth
		quNp(backPrj->northPoleQuat()                 ),//get back projection rotation
		wrk (imPrc->numPix()                          ),//allocate space for image processing
		sph (dim * dim * 2                            ),//allocate space for back projection
		gln (mBw  * mBw                               ),//allocate space for SHT of a single pattern
		prc (std::move(imPrc)                         ),//copy image processor
		prj (std::move(backPrj)                       ),//copy back projector
		sht (dim, dim - 2, square::Layout::Legendre   ) //build spherical harmonic transformer
		{
			pSym.reserve(corrs.size());
			xc  .reserve(corrs.size());
			for(const std::unique_ptr< sphere::PhaseCorrelator<Real> >& ptr : corrs) {
				pSym.push_back( std::vector< xtal::Quat<Real> >() );//place holder for psuedosymmetric equivalents
				xc  .push_back(ptr->clone());
			}
	}

	//@brief   : estimate a reasonable batch size
	//@param bw: bandwidth
	//@param nt: number of worker threads
	//@param np: number of patterns to index
	//@return  : a reasonable batch size (a number of patterns that should take on the order of 1s to index
	template <typename Real>
	size_t Indexer<Real>::BatchEstimate(const size_t bw, const size_t nt, const size_t np) {
		//first compute an estimate assuming a large number of patterns to index compared to the available thread count
		const double bw3 = double(bw * bw * bw);
		const double scl = bw3 * std::log(bw3);//complexity scaling is O(n^3 * ln(n^3))
		const double k = 1E-8;//scaling factor on my laptop (should be a good enough estimate for all modern computers)
		const double tPat = scl * k;//estimate time per pattern per thread
		const double pps = 1.0 / tPat;//number of patterns that can be indexed in 1 second on 1 thread
		size_t batchSize = std::max<size_t>(1, (size_t) (pps * 0.61803398874989484820458683436564));//start with a batch that takes ~1/golden ratio seconds (that way it won't sync with updates which could make early time estimates inaccurate)

		//next make sure that we don't have unused threads for small numbers of patterns
		const size_t numBatch = (size_t) std::ceil(double(np) / batchSize);//total number of batches
		if(numBatch < nt * nt) {//make sure there are enough batches for load balancing
			double newBatch = double(np) / (nt * nt);
			batchSize = (size_t) std::ceil(newBatch);
		}
		return batchSize;
	}

	//@brief       : index a single image
	//@param    pat: pointer to image to index in row major order
	//@param    n  : number of results to keep (top n), user is responsible for making sure this is reasonable (extra points will be filled with an invalid phase)
	//@param    res: location to write results
	//@param    ref: true/false to use newton's method refinement/triquadratic interpolation
	//@template Pix: pixel type of image (will be cast to Real)
	//@note        : multiple results can be written (currently there is an additional result per psuedosymmetric orientation to check)
	template <typename Real>
	template <typename Pix>
	void Indexer<Real>::indexImage(Pix* pat, Result<Real>*const res, const size_t n, const bool ref) {
		//fill output with empty result
		for(size_t i = 0; i < n; i++) {
			res[i].corr  =  0;//only keep something with a positive phase
			res[i].phase = -1;//start with 'bad' phase
			std::fill(res[i].qu, res[i].qu+4, Real(0));
		}
		const Real iq = computeHarmonics(pat);//do image processing, back project to sphere, and compute spherical harmonic transform

		//loop over phases
		xtal::Quat<Real> q0;
		Real eu[3];
		for(size_t p = 0; p < pSym.size(); p++) {
			//compute cross correlation with phase
			Result<Real> r = correlate(p, gln.data(), ref);
			r.iq = iq;
			r.phase = (int)p;

			//save the result if it is good enough
			size_t idx = std::distance(res, std::upper_bound(res, res+n, r));//determine where in sorted list of outputs this falls (should we keep it)
			if(idx < n) {//this is good enough to bump at least 1 result of the list
				for(size_t i = n-1; i != idx; i--) res[i] = res[i-1];//loop from back to front shifting results down one
				res[idx] = r;//insert this result into the list
			}

			//check pseudo symmetric misorientations
			xtal::zyz2qu(r.qu, q0.data());//convert to quaternion
			for(const xtal::Quat<Real>& q : pSym[p]) {
				//compute psuedo-symmetric orientation (order may be wrong here)
				//q0 is the quaternion to rotate from crystal -> sample
				//psuedo sym is in crystal frame
				//q may need to be conjugated depending on how they are input
				xtal::Quat<Real> qp = q0 * q;//first apply psuedo symmetry (q) then registered orientation (q0)
				xtal::qu2zyz(qp.data(), eu);//convert back to to euler angles

				//do refinement
				r = refine(p, gln.data(), eu);
				r.iq = iq;

				//save the result if it is good enough
				idx = std::distance(res, std::upper_bound(res, res+n, r));//determine where in sorted list of outputs this falls (should we keep it)
				if(idx < n) {//this is good enough to bump at least 1 result of the list
					for(size_t i = n-1; i != idx; i--) res[i] = res[i-1];//loop from back to front shifting results down one
					res[idx] = r;//insert this result into the list
				}
			}
		}

		//loop over results converting to quaterions
		for(size_t i = 0; i < n; i++) {
			xtal::zyz2qu(res[i].qu, res[i].qu);//convert to quaternion
			xtal::quat::mul(Indexer<Real>::quNp.data(), res[i].qu, res[i].qu);//correct for rotated detector frame
			for(size_t j = 1; j < 4; j++) res[i].qu[j] = -res[i].qu[j];//conjugate result (crystal->sample to sample->crystal)
		}
	}

	//@brief       : refine the orientation of a single image
	//@param    pat: pointer to image to refine in row major order
	//@param    res: location to write result (and read initial phase + orientation)
	//@template Pix: pixel type of image (will be cast to Real)
	template <typename Real>
	template <typename Pix>
	void Indexer<Real>::refineImage(Pix* pat, Result<Real>& res) {
		//compute spectra of pattern
		const Real iq = computeHarmonics(pat);//do image processing, back project to sphere, and compute spherical harmonic transform

		//unconjugate result
		xtal::Quat<Real> q0;
		xtal::quat::conj(res.qu, q0.data());

		//uncorrect for rotated detector frame
		xtal::Quat<Real> npC;
		xtal::quat::conj(Indexer<Real>::quNp.data(), npC.data());
		xtal::quat::mul (npC.data(), q0.data(), q0.data());

		//convert to euler angles
		Real eu[3];
		xtal::qu2zyz(q0.data(), eu);//convert to euler angles

		//do refinement
		refine(res.phase, gln.data(), eu);
		xtal::zyz2qu(eu, res.qu);

		//re-correct for rotated detector frame
		xtal::quat::mul(Indexer<Real>::quNp.data(), res.qu, res.qu);

		//re-conjugate result
		for(size_t i = 1; i < 4; i++) res.qu[i] = -res.qu[i];
		res.iq = iq;
	}

	//@brief   : process an image, unproject to the sphere, and compute spherical harmonic coefficients
	//@param im: image to process
	//@return  : image quality
	template <typename Real>
	template <typename Pix>
	Real Indexer<Real>::computeHarmonics(Pix* im) {
		Real iq;		
		prc->process(im, wrk.data());
		sum2 = prj->unproject(wrk.data(), sph.data(), &iq);
		sht.analyze(sph.data(), gln.data(), mBw, mBw);//compute SHT of function in sph
		return iq;
	}

	//@brief    : index the spherical function currently held in sph
	//@param p  : phase to correlate with
	//@param gln: spectra of image to index
	//@param ref: true/false to use newton's method refinement/triquadratic interpolation
	//@return   : correlation result (with ZYZ euler angle stored in qu[0,1,2])
	template <typename Real>
	Result<Real> Indexer<Real>::correlate(const size_t p, std::complex<Real> const * gln, const bool ref) {
		Result<Real> r;
		r.corr  = xc[p]->correlate(Indexer<Real>::gln.data(), r.qu, ref);
		r.phase = (int)p;
		return r;
	}

	//@brief    : refine the spherical function currently held in sph
	//@param p  : phase to refine with
	//@param gln: spectra of image to refine
	//@param eu : initial ZYZ euler angle to refine from
	//@return   : refinement result (with ZYZ euler angle stored in qu[0,1,2])
	template <typename Real>
	Result<Real> Indexer<Real>::refine(const size_t p, std::complex<Real> const * gln, Real const * const eu) {
		Result<Real> r;
		std::copy(eu, eu+3, r.qu);
		r.corr  = xc[p]->refinePeak(gln, r.qu);
		r.phase = (int)p;
		return r;
	}

}

#endif//_SPHERE_INDEXER_H_

