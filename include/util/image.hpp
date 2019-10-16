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

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include <memory>
#include <vector>

#include "fft.hpp"

namespace image {

	//@brief     : convert an image to 8 bit
	//@param im  : image to convert
	//@param nPix: number of pixels
	//@param p8  : location to write 8 bit image
	//@note      : intensities are rescaled so that [min,max] maps to [0,255]
	template <typename TPix>
	void to8Bit(TPix * const im, const size_t nPix, uint8_t * const p8);

	//@brief     : convert an image to 8 bit
	//@param im  : image to convert
	//@param nPix: number of pixels
	//@note      : intensities are rescaled so that [min,max] maps to [0,255]
	//@return    : rescaled im
	template <typename TPix>
	std::vector<uint8_t> to8Bit(TPix * const im, const size_t nPix);

	//@brief     : adaptive histogram equalization
	//@param im  : image to equalize
	//@param w   : width of image in pixels
	//@param h   : height of image in pixels
	//@param nx  : number of tiles (contrast regions) across image
	//@param ny  : number of tiles (contrast regions) down image
	//@param vMin: minimum histogram value (NAN to compute from input image)
	//@param vMax: maximum histogram value (NAN to compute from input image)
	//@reference : Pizer, S. M., Amburn, E. P., Austin, J. D., Cromartie, R., Geselowitz, A., Greer, T., ... & Zuiderveld, K. (1987). Adaptive histogram equalization and its variations. Computer vision, graphics, and image processing, 39(3), 355-368.
	//@note      : this currently uses 'mosaic' sampling (see fig 3a)
	//@note      : NAN pixels are ignored for equalization and then replaced with 0 afterwards
	template <typename TPix>
	void adHistEq(TPix * const im, const size_t w, const size_t h, const size_t nx, const size_t ny, double vMin = NAN, double vMax = NAN);

	//@brief     : compute the histogram of an image
	//@param im  : image to compute histogram of
	//@param w   : width of image in pixels
	//@param h   : height of image in pixels
	//@param bins: location to write values of histogram bins
	//@param cnts: number of pixels in each bin
	//@param nBin: number of histogram bins to use
	template <typename TPix>
	void hist(TPix * const im, const size_t w, const size_t h, TPix * bins, size_t * cnts, const size_t nBin = 256);

	//@brief     : compute the otsu threshold for an image from its hisogram
	//@param bins: location to write values of histogram bins
	//@param cnts: number of pixels in each bin
	//@param nBin: number of histogram bins to use
	//@return    : index of threshold value in bins
	template <typename TPix>
	size_t otsu(TPix * bins, size_t * cnts, const size_t nBin = 256);

	//@brief     : compute image qualty of a spectrum
	//@param dct: discrete cosine transform of image to compute IQ of
	//@param w   : width of image
	//@param h   : height of image
	//@return    : image quality
	template <typename Real>
	Real imageQuality(Real * dct, const size_t w, const size_t h);

	//@brief: helper for computing image quality
	template <typename Real>
	class ImageQualityCalc {
		fft::vector    <             Real  > wrk;//work space for frequency domain
		std::shared_ptr< fft::DCT2D <Real> > fwd;//plan for forward DCT
		public:
			const size_t w, h;//width/height of image to compute image quality of

			//@brief    : construct an image quality calculator
			//@param w  : width  of  images 
			//@param h  : height of  images 
			//@param flg: fft flag to use when creating DCT plans
			ImageQualityCalc(const size_t w, const size_t h, const fft::flag::Plan flg = fft::flag::Plan::Measure) : w(w), h(h), fwd(std::make_shared<fft::DCT2D <Real> >(w, h, flg, false)), wrk(w * h) {}

			//@brief   : compute image quality of an image
			//@param im: image to compute quality of
			//@return  : image quality
			double compute(Real const * im) {fwd->execute((Real*)im, wrk.data()); return imageQuality(wrk.data(), w, h);}
	};
	
	//helper for bilinearly interpolation
	template <typename Real>
	struct BiPix {
		size_t idx    ;//index of pixel in square sphere grid
		size_t inds[4];//index of 4 bounding pixels on detector
		Real   wgts[4];//relative weights of pixels to interpolate from

		//@brief    : bilinearly interpolate a pattern at the spherical pixel
		//@param pat: detector pattern to interpolate from
		//@return   : interpolated value
		Real interpolate(Real const * const pat) const;

		//@brief     : compute the indicies and coefficient to bilinearlly interpolate from an image
		//@param x   : fractional horizontal position in image [0,1]
		//@param y   : fractional vertical position in image [0,1]
		//@param w   : width of image in pixels
		//@param h   : height of image in pixels
		void bilinearCoeff(Real x, Real y, const size_t w, const size_t h);
	};

	//helper for resacling images
	template <typename Real>
	class Rescaler {
		fft::vector    <             Real  > work      ;//work space for frequency domain
		std::shared_ptr< fft::DCT2D <Real> > fwd , inv ;//plan for forward/inverse DCT
		public:
			const size_t   wIn , hIn ;//width/height of input  image
			const size_t   wOut, hOut;//width/height of output image

			//@brief: construct a scaler to rescale images from wI x hI -> wO x hO
			//@param wI : width  of input  images 
			//@param hI : height of input  images 
			//@param wO : width  of output images 
			//@param hO : height of output images 
			//@param flg: fft flag to use when creating DCT plans
			Rescaler(const size_t wI, const size_t hI, const size_t wO, const size_t hO, const fft::flag::Plan flg);

			//@brief: construct a scaler to rescale images from wI x hI -> (wI x hI) * scl
			//@param wI : width  of input  images 
			//@param hI : height of input  images 
			//@param scl: scale factor
			//@param flg: fft flag to use when creating DCT plans
			Rescaler(const size_t wI, const size_t hI, const Real scl, const fft::flag::Plan flg) : Rescaler(wI, hI, (size_t)std::round(scl * wI), (size_t)std::round(scl * hI), flg) {}

			//@brief    : rescale an image
			//@param in : input image to rescale (wIn * hIn)
			//@param out: location to write rescaled imaged (wOut * hOut), can be the same as in
			//@param zer: true/false to make average of the image zero
			//@param flt: width of high pass filter (0 to leave unfiltered) [dc value isn't affected]
			//@param iq : should the image quality be computed
			//@return   : iq ? image quality : 0
			Real scale(Real* in, Real* out, bool zer = false, const size_t flt = 0, const bool iq = false) {return scale(in, out, work.data(), zer, flt, iq);}

			//@breif: allocate space large enough for the working array
			fft::vector<Real> allocateWork() const {return fft::vector<Real>(work.size());}

			//@brief    : rescale an image using a user provided work array
			//@param in : input image to rescale (wIn * hIn)
			//@param out: location to write rescaled imaged (wOut * hOut), can be the same as in
			//@param wrk: work space (should be allocated via allocateWork())
			//@param zer: true/false to make average of the image zero
			//@param flt: width of high pass filter (0 to leave unfiltered) [dc value isn't affected]
			//@param iq : should the image quality be computed
			//@return   : iq ? image quality : 0
			Real scale(Real* in, Real* out, Real* wrk, bool zer = false, const size_t flt = 0, const bool iq = false) const;
	};

	//helper for subtracting backgrounds from non-rectangular images
	template <typename Real>
	class BlockRowBackground {
		std::shared_ptr< std::vector< std::pair<size_t, size_t> > > msk ;//range of pixels within circular mask for each row
		std::shared_ptr< std::vector<           size_t          > > yCnt;//number of pixels in each column
		const size_t                                                w, h;//width/height of image to subtract background
		std::shared_ptr< fft::DCT <Real> >                          xFwd;//plane for forward DCT of a single row
		std::shared_ptr< fft::DCT <Real> >                          yFwd;//plane for forward DCT of a single column
		std::shared_ptr< fft::DCT <Real> >                          xInv;//plane for inverse DCT of a single row
		std::shared_ptr< fft::DCT <Real> >                          yInv;//plane for inverse DCT of a single column
		fft::vector    <           Real  >                          wrk1;//work space
		fft::vector    <           Real  >                          wrk2;//work space
		fft::vector    <           Real  >                          wrk3;//work space

		//@brief      : private constructor to build a block row background subtractor from row blocks + size
		//@param rMsk : row blocks (height is size)
		//@param width: width of image
		BlockRowBackground(std::vector< std::pair<size_t, size_t> >& rMsk, const size_t width);

		public:
			//@brief  : construct a background subtractor for a circular mask inscribed in an image (centered)
			//@param w: width of image
			//@param h: height of image
			//@return : background subtractor
			static BlockRowBackground Circular(const size_t w, const size_t h);

			//@brief   : subtract x and y backgrounds from inside the block region
			//@param im: image to subtract background from
			//@param x : width of high pass filter for x direction (0 to leave unfiltered)
			//@param y : width of high pass filter for y direction (0 to leave unfiltered)
			void subtract(Real * const im, const size_t x, const size_t y);
	};
}

#include <algorithm>
#include <numeric>
#include <limits>
#include <type_traits>

#include "constants.hpp"

namespace image {
	////////////////////////////////////////////////////////////////////////
	//                               BiPix                                //
	////////////////////////////////////////////////////////////////////////

	//@brief     : convert an image to 8 bit
	//@param im  : image to convert
	//@param nPix: number of pixels
	//@param p8  : location to write 8 bit image
	//@note      : intensities are rescaled so that [min,max] maps to [0,255]
	template <typename TPix>
	void to8Bit(TPix * const im, const size_t nPix, uint8_t * const p8) {
		//handle trivial case first
		if(std::is_same<TPix, uint8_t>::value) {
			std::transform(im, im + nPix, p8, [](const TPix& v){return (uint8_t)v;});//casting needed to supress compiler warnings on nont same types
			return;
		}

		//compute minimum and maximum quality
		std::pair<TPix*, TPix*> minMax = std::minmax_element(im, im + nPix);
		const TPix delta = *minMax.first;

		//rescale to 8 bit range
		if(std::numeric_limits<TPix>::is_integer) {
			const double scale = double(255) / (*minMax.second - delta);
			std::transform(im, im + nPix, p8, [&](const TPix& v){return (uint8_t)std::round(scale * (v - delta));});
		} else {
			const TPix   scale = TPix  (255) / (*minMax.second - delta);
			std::transform(im, im + nPix, p8, [&](const TPix& v){return (uint8_t)std::round(scale * (v - delta));});
		}
	}

	//@brief     : convert an image to 8 bit
	//@param im  : image to convert
	//@param nPix: number of pixels
	//@note      : intensities are rescaled so that [min,max] maps to [0,255]
	//@return    : rescaled im
	template <typename TPix>
	std::vector<uint8_t> to8Bit(TPix * const im, const size_t nPix) {
		// if(std::is_same<TPix, uint8_t>::value) return std::vector<uint8_t>(im, im+nPix);//handle trivial case
		std::vector<uint8_t> ret(nPix);//allocate output space
		to8Bit(im, nPix, ret.data());//do conversion
		return ret;//return converted image
	}

	//@brief     : adaptive histogram equalization
	//@param im  : image to equalize
	//@param w   : width of image in pixels
	//@param h   : height of image in pixels
	//@param nx  : number of tiles (contrast regions) across image
	//@param ny  : number of tiles (contrast regions) down image
	//@param vMin: minimum histogram value (NAN to compute from input image)
	//@param vMax: maximum histogram value (NAN to compute from input image)
	//@reference : Pizer, S. M., Amburn, E. P., Austin, J. D., Cromartie, R., Geselowitz, A., Greer, T., ... & Zuiderveld, K. (1987). Adaptive histogram equalization and its variations. Computer vision, graphics, and image processing, 39(3), 355-368.
	//@note      : this currently uses 'mosaic' sampling (see fig 3a)
	//@note      : NAN pixels are ignored for equalization and then replaced with 0 afterwards
	template <typename TPix>
	void adHistEq(TPix * const im, const size_t w, const size_t h, const size_t nx, const size_t ny, double vMin, double vMax) {
		//first get histogram limits (minmax_element doesn't work with nans)
		// const std::pair<TPix*, TPix*> minMax = std::minmax_element(im, im + w * h);//get limits of image
		// const double vMin = double(*minMax.first );//darkest pixel value
		// const double vMax = double(*minMax.second);//brightest pixel value
		if(std::isnan(vMin) || std::isnan(vMax)) {
			const double vMin0 = vMin;//save initial vMin in case just vMax is NAN
			const double vMax0 = vMax;//save initial vMax in case just vMin is NAN (not currently possible)
			bool init = false;
			vMin = vMax = im[0];
			for(size_t i = 0; i < w * h; i++) {
				if(!std::isnan(im[i])) {
					if(!init) {
						vMin = vMax = im[i];
						init = true;
					} else {
						if(im[i] < vMin) vMin = im[i];
						if(im[i] > vMax) vMax = im[i];
					}
				}
			}
			if(!std::isnan(vMin0)) vMin = vMin0;
			if(!std::isnan(vMax0)) vMax = vMax0;
		}

		//create bins for histograms
		const size_t nBins = 256;//number of histogram bins
		const double range = vMax - vMin;//get range
		if(0.0 == range) return;//no contrast to equalize
		const double delta = range / nBins;//get histogram spacing
		std::vector<double> vBins(nBins);//allocate histogram
		for(size_t i = 0; i < nBins; i++) vBins[i] = vMin + delta * (i+1);//fill in bins
		vBins.back() = vMax;//make sure rounding errors don't put the last bin below the max

		//determine tile size
		const double tx = double(w) / nx;//x tile size in fractional pixels
		const double ty = double(h) / ny;//y tile size in fractional pixels
		const double hy = ty / 2;//y half tile size in fractional pixels

		//compute x tile assignments and fractional positions once
		std::vector<size_t> it(w);//which tile is each pixel in
		std::vector<size_t> il(w);//which grid point is to the left of each pixel
		std::vector<size_t> ir(w);//which grid point is to the right of each pixel
		std::vector<double> fx(w), cx(w);
		for(size_t i = 0; i < w; i++) {//loop over columns
			double f = double(i) / tx;//convert from pixels to fractional tiles
			it[i] = (size_t)f;//save the tile this pixel is in
			il[i] = (size_t)std::max<double>(f - 0.5, 0.0         );
			ir[i] = (size_t)std::min<double>(f + 0.5, double(nx-1));
			fx[i] = std::fmod(f - it[i] + 0.5, 1.0);//save progress between tile centers
			cx[i] = 1.0 - fx[i];//save complement of fx[i]
		}

		//allocate a histogram for each x tile
		std::vector<size_t> hist(nx * nBins);//histogram for current row of tiles
		std::vector<double> cdfPrev(nx * nBins);//cumulative histogram for previous row of tiles
		std::vector<double> cdfCur (nx * nBins);//cumulative histogram for current row of tiles
		std::vector<size_t> wrk(nBins);//work space for partial sum

		//loop over rows of tiles
		for(size_t k = 0; k < ny; k++) {
			//compute y coordinates of tile row bounds
			const double yS = ty *  k   ;//where does this tile row start
			const double yE = ty * (k+1);//where does the tile row end

			//convert y coordinates to nearest pixel
			const size_t jS = (size_t) std::round(yS);
			const size_t jE = (size_t) std::round(yE);

			//now loop over this tile computing the histogram for each x tile
			std::fill(hist.begin(), hist.end(), 0);//initialize histograms
			for(size_t j = jS; j < jE; j++) {//loop over rows within this tile
				for(size_t i = 0; i < w; i++) {//loop over width of image
					const TPix& v = im[j * w + i];//get pixel value
					if(!std::isnan(v)) {//don't count NaN pixels
						const auto iBin = std::lower_bound(vBins.cbegin(), vBins.cend(), v);//find the bin this pixel belongs to
						const size_t b = std::distance(vBins.cbegin(), iBin);//convert to bin index
						++hist[it[i] * nBins + b];//increment bin in correct histogram of current tile
					}
				}
			}

			//next convert from histograms to CDFs
			for(size_t i = 0; i < nx; i++) {
				// std::partial_sum(hist.cbegin() + i * nBins, hist.cbegin() + (i+1) * nBins, cdfCur.begin() + i * nBins);//raw cdf
				std::partial_sum(hist.cbegin() + i * nBins, hist.cbegin() + (i+1) * nBins, wrk.begin());//raw cdf
				std::transform(wrk.cbegin(), wrk.cend(), cdfCur.begin() + i * nBins, [](const size_t& sum){return (double)sum;});
				const double nrm = range / cdfCur[i * nBins + nBins - 1];//compute normalization to that CDF goes from [0, range]
				if(0 == cdfCur[i * nBins + nBins - 1]) {//don't divide by zero
					std::copy(vBins.begin(), vBins.end(), cdfCur.begin() + i * nBins);//entire window is nan, use linear ramp
				} else {//normal case, adjusth histogram
					std::for_each(cdfCur.begin() + i * nBins, cdfCur.begin() + (i+1) * nBins, [nrm, vMin](double& v){v = v * nrm + vMin;});//rescale histogram to [vMin, vMax]
				}
			}

			//handle first row of tiles specially
			if(k == 0) cdfPrev = cdfCur;//just make previous row the same as this one for first row

			//now loop over the parts of the tile we have the bounding points for and apply equalization
			const double yM     = yS + hy;//mid point of this tile row
			const double yMPrev = yS - hy;//mid point of previous tile row
			const size_t jMPrev = k   == 0  ? 0 : (size_t) std::round(yMPrev);//index of first row to compute values for
			const size_t jM     = k+1 == ny ? h : (size_t) std::round(yM    );//index of last row to compute values for
			for(size_t j = jMPrev; j < jM; j++) {//loop over rows (this is previous to current midpoint for normal rows but cliped/extened for first/last row)
				const double y = double(j);
				const double fy = std::min((double(j) - yMPrev) / ty, 1.0);//compute fractional progress between tiles (clamp to 1 for bottom half of last row)
				const double cy = 1.0 - fy;
				for(size_t i = 0; i < w; i++) {//loop over image columns
					TPix& v = im[j * w + i];//get pixel value
					if(!std::isnan(v)) {//don't change NaN pixels
						const auto iBin = std::lower_bound(vBins.cbegin(), vBins.cend(), v);//find the bin this pixel belongs to
						const size_t b = std::distance(vBins.cbegin(), iBin);//convert to bin index
						// const double grid[4] = {//get values of histogram equalization using 4 neighboring grid points
						double grid[4] = {//get values of histogram equalization using 4 neighboring grid points
							cdfPrev[nBins * il[i] + b],
							cdfPrev[nBins * ir[i] + b],
							cdfCur [nBins * il[i] + b],
							cdfCur [nBins * ir[i] + b]
						};
						//bilinearly interpolate grid points
						v = grid[0] * cy * cx[i]
						  + grid[1] * cy * fx[i]
						  + grid[2] * fy * cx[i]
						  + grid[3] * fy * fx[i];
					} else {
						v = 0;//replace nans with 0
					}
				}
			}
			cdfCur.swap(cdfPrev);
		}
	}

	//@brief     : compute the histogram of an image
	//@param im  : image to compute histogram of
	//@param w   : width of image in pixels
	//@param h   : height of image in pixels
	//@param bins: location to write values of histogram bins
	//@param cnts: number of pixels in each bin
	//@param nBin: number of histogram bins to use
	template <typename TPix>
	void hist(TPix * const im, const size_t w, const size_t h, TPix * bins, size_t * cnts, const size_t nBin) {
		//compute bins
		if(256 == nBin && std::is_same<TPix, std::uint8_t>::value) {
			//handle 8 bit images specially
			std::iota(bins, bins + 256, 0);
		} else {
			//get limits of image
			const std::pair<TPix*, TPix*> minMax = std::minmax_element(im, im + w * h);
			const double vMin = double(*minMax.first );//darkest pixel value
			const double vMax = double(*minMax.second);//brightest pixel value

			//convert to range
			const double range = vMax - vMin;//get range
			const double delta = range / nBin;//get histogram spacing
			for(size_t i = 0; i < nBin; i++) bins[i] = vMin + delta * (i+1);//fill in bins
			bins[nBin-1] = vMax;//make sure rounding errors don't put the last bin below the max
		}

		//now loop over image filling histogram
		const size_t nPix = w * h;
		std::fill(cnts, cnts + nBin, 0);
		for(size_t i = 0; i < nPix; i++) ++cnts[std::distance(bins, std::lower_bound(bins, bins + nBin, im[i]))];
	}

	//@brief     : compute the otsu threshold for an image from its hisogram
	//@param bins: location to write values of histogram bins
	//@param cnts: number of pixels in each bin
	//@param nBin: number of histogram bins to use
	//@return    : index of threshold value in bins
	template <typename TPix>
	size_t otsu(TPix * bins, size_t * cnts, const size_t nBin) {
		//compute number of pixels
		const size_t numPix = std::accumulate(cnts, cnts + nBin, size_t(0));

		//compute mean pixel
		double ut = std::inner_product(bins, bins + nBin, cnts, 0.0) / numPix;

		//now compute optimum threshold
		double wk   = double(cnts[0]) / numPix;//w(k)
		double uk   = wk * bins[0]            ;//mu(k)
		double sMax = 0                       ;//maximum sigma^2_b(k)
		size_t kMax = 0                       ;//index of sMak
		for(size_t k = 1; k < nBin - 1; k++) {
			const double pk = double(cnts[k]) / numPix;//convert from count to probability
			wk += pk          ;//update w(k)
			uk += pk * bins[k];//update u(k)`
			const double f = ut * wk - uk;//numerator prefactor
			const double s2 = f * f / (wk * (1.0 - wk));//compute merit sigma^2_b(k)
			if(s2 > sMax) {//is this a new best threshold?
				sMax = s2;//save best value
				kMax = k ;//save best index
			}
		}
		return kMax;
	}

	//@brief     : compute image qualty of a spectrum
	//@param dct: discrete cosine transform of image to compute IQ of
	//@param w   : width of image
	//@param h   : height of image
	//@return    : image quality
	template <typename Real>
	Real imageQuality(Real * dct, const size_t w, const size_t h) {
		Real vIq(0), sumP(0);
		const uint64_t numPix = w * h;//compute number of pixels
		uint64_t sumW(0);//this is big enough to hold the sum of the weights of at least a ~55000^2 image without overflow
		for(uint64_t j = 0; j < h; j++) {//loop over spectra rows
			uint64_t j2 = j * j;//compute y distance from origin^2
			for(uint64_t i = 0; i < w; i++) {//loop over spectra columns
				uint64_t r2 = j2 + i * i;//compute weighting for this position in spectrum (radius^2)
				const Real p = std::fabs(dct[j * w + i]);//get power
				vIq  += p * r2;//accumulate numerator (radius weighted power)
				sumP += p     ;//accumlate power
				sumW +=     r2;//accumulate weighting
			}
		}
		const Real den = Real(sumW) * sumP / numPix;
		if(sumP == Real(0)) vIq = 0;//if summed power is 0 there is no image (or only a DC value)
		else vIq = Real(1) - vIq / den;//normalize by power
		return vIq;
	}

	//@brief    : bilinearly interpolate a pattern at the spherical pixel
	//@param pat: detector pattern to interpolate from
	//@return   : interpolated value
	template <typename Real>
	Real BiPix<Real>::interpolate(Real const * const pat) const {
		Real v = 0;
		for(size_t i = 0; i < 4; i++) v += pat[inds[i]] * wgts[i];
		return v;
	}

	//@brief     : compute the indicies and coefficient to bilinearlly interpolate from an image
	//@param x   : fractional horizontal position in image [0,1]
	//@param y   : fractional vertical position in image [0,1]
	//@param w   : width of image in pixels
	//@param h   : height of image in pixels
	//@return    : BiPix holding appropriate weights and indices
	template <typename Real>
	void BiPix<Real>::bilinearCoeff(Real x, Real y, const size_t w, const size_t h) {
		//determine 4 bounding pixels
		x *= w-1; y *= h-1;
		const size_t indX0 = std::min((size_t)x, w-1);
		const size_t indY0 = std::min((size_t)y, h-1);
		const size_t indX1 = std::min(indX0+1, w - 1);//handle x == 1
		const size_t indY1 = std::min(indY0+1, h - 1);//handle y == 1

		//convert to vectorized indices
		inds[0] = indY0 * w + indX0;
		inds[1] = indY0 * w + indX1;
		inds[2] = indY1 * w + indX0;
		inds[3] = indY1 * w + indX1;

		//compute fractional progress between bounding pairs (and complements)
		Real wx1 = x - indX0;
		Real wy1 = y - indY0;
		Real wx0 = Real(1) - wx1;
		Real wy0 = Real(1) - wy1;

		//convert to pixel weights
		wgts[0] = wy0 * wx0;
		wgts[1] = wy0 * wx1;
		wgts[2] = wy1 * wx0;
		wgts[3] = wy1 * wx1;
	}

	////////////////////////////////////////////////////////////////////////
	//                              Rescaler                              //
	////////////////////////////////////////////////////////////////////////

	//@brief: construct a scaler to rescale images from wI x hI -> wO x hO
	//@param wI : width  of input  images 
	//@param hI : height of input  images 
	//@param wO : width  of output images 
	//@param hO : height of output images 
	//@param flg: fft flag to use when creating DCT plans
	template <typename Real>
	Rescaler<Real>::Rescaler(const size_t wI, const size_t hI, const size_t wO, const size_t hO, const fft::flag::Plan flg) : 
		wIn (wI),
		hIn (hI),
		wOut(wO),
		hOut(hO),
		fwd (std::make_shared<fft::DCT2D <Real> >(wI, hI, flg, false)),
		inv (std::make_shared<fft::DCT2D <Real> >(wO, hO, flg, true )),
		work(std::max(wI, wO) * std::max(hI, hO)) {}

	//@brief    : rescale an image
	//@param in : input image to rescale (wIn * hIn)
	//@param out: location to write rescaled imaged (wOut * hOut), can be the same as in
	//@param wrk: work space (should be allocated via allocateWork())
	//@param zer: true/false to make average of the image zero
	//@param flt: width of high pass filter (0 to leave unfiltered) [dc value isn't affected]
	//@param iq : should the image quality be computed
	//@return   : iq ? image quality : 0
	template <typename Real>
	Real Rescaler<Real>::scale(Real* in, Real* out, Real* wrk, bool zer, const size_t flt, const bool iq) const {
		//compute dct of input image (and image quality if needed)
		fwd->execute(in, wrk);//forward dct
		Real vIq = iq ? imageQuality(wrk, wIn, hIn) : Real(0);

		//truncate/pad to new size in frequency domain
		const size_t hCpy = std::min(hIn, hOut);//number of rows to copy
		if(wOut <= wIn) {//width constant or decreasing
			for(size_t j = 0; j < hCpy; j++) {//loop over rows to copy in order
				std::copy(wrk + j * wIn, wrk + j * wIn + wOut, wrk + j * wOut);//copy row by row
			}
		} else {//width is increasing
			for(size_t j = hCpy - 1; j < hCpy; j--) {//loop over rows to copy in reverse order
				for(size_t i = wIn - 1; i < wIn; i--) wrk[j*wOut + i] = wrk[j*wIn + i];
				std::fill(wrk + j * wOut + wIn, wrk + j * wOut + wOut, Real(0));//pad row to width 
			}
		}
		if(hOut > hIn) {//height is increasing
			std::fill(wrk + wOut * hIn, wrk + wOut * hOut, Real(0));//pad to height
		}

		//filter and compute inverse dct
		if(zer) wrk[0] = 0;//set dc value to 0 -> make mean value of image 0
		if(flt > 0) {
			//filter out long wavelengths to subtract non dc background
			for(size_t j = 0; j < flt; j++) {
				for(size_t i = 0; i < flt; i++) {
					const Real r = std::sqrt(Real(j*j+i*i));
					if(r <= flt && (i > 0 || j > 0)) {
						const Real factor = std::cos((r / (2 * flt) + Real(0.5)) * emsphinx::Constants<Real>::pi);
						wrk[j * wOut + i] *= factor * factor;
					}
				}
			}
		}
		inv->execute(wrk, out);//inverse dct
		return vIq;
	}

	//@brief      : private constructor to build a block row background subtractor from row blocks + size
	//@param rMsk : row blocks (height is size)
	//@param width: width of image
	template <typename Real>
	BlockRowBackground<Real>::BlockRowBackground(std::vector< std::pair<size_t, size_t> >& rMsk, const size_t width) : 
		msk (std::make_shared<std::vector< std::pair<size_t, size_t> > >(rMsk)),
		yCnt(std::make_shared<std::vector< size_t > >(width)),
		w   (width         ),
		h   (rMsk.size()   ),
		wrk1(std::max(w, h)),
		wrk2(std::max(w, h)),
		wrk3(std::max(w, h)),
		xFwd(std::make_shared<fft::DCT<Real> >(w, fft::flag::Plan::Measure, false)),
		yFwd(std::make_shared<fft::DCT<Real> >(h, fft::flag::Plan::Measure, false)),
		xInv(std::make_shared<fft::DCT<Real> >(w, fft::flag::Plan::Measure, true )),
		yInv(std::make_shared<fft::DCT<Real> >(h, fft::flag::Plan::Measure, true )) {
		//compute number of pixels in each column
		std::fill(yCnt->begin(), yCnt->end(), 0);
		for(size_t j = 0; j < h; j++) {
			const std::pair<size_t, size_t>& range = msk->operator[](j);
			for(size_t i = range.first; i < range.second; i++) {
				++yCnt->operator[](i);
			}
		}
	}


	//@brief  : construct a background subtractor for a circular mask inscribed in an image (centered)
	//@param w: width of image
	//@param h: height of image
	//@return : background subtractor
	template <typename Real>
	BlockRowBackground<Real> BlockRowBackground<Real>::Circular(const size_t w, const size_t h) {
		std::vector< std::pair<size_t, size_t> > mask(h);
		const size_t dim = std::min(w, h) - 1;
		const Real cutoff = Real(1) + Real(1) / (2 * dim);
		for(size_t j = 0; j < h; j++) {//loop over rows
			std::pair<size_t, size_t>& range = mask[j];
			range.first = w;
			range.second = w;
			const Real y = Real(j * 2) / dim - 1;//[0,1]
			for(int i = 0; i < w; i++) {
				const Real x = Real(i * 2) / dim - 1;//[0,1]
				const Real r = std::hypot(x, y);
				if(r <= cutoff) {//inside circle
					if(w == range.first) range.first = i;//if this is the first point inside the circle save the start
				} else if(w != range.first) {//this point is outside the circle but after at least one pixel inside the circle
					if(w == range.second) range.second = i;//save first point outside circle
				}
			}
		}
		return BlockRowBackground(mask, w);
	}

	//@brief   : subtract x and y backgrounds from inside the block region
	//@param im: image to subtract background from
	//@param x : width of high pass filter for x direction (0 to leave unfiltered)
	//@param y : width of high pass filter for y direction (0 to leave unfiltered)
	template <typename Real>
	void BlockRowBackground<Real>::subtract(Real * const im, const size_t x, const size_t y) {
		//initalize averages with zeros
		std::fill(wrk1.begin(), wrk1.begin() + h, 0);//average of rows
		std::fill(wrk2.begin(), wrk2.begin() + w, 0);//average of columns

		//loop over image accumulating region inside mask
		for(size_t j = 0; j < h; j++) {
			const std::pair<size_t, size_t>& range = msk->operator[](j);
			for(size_t i = range.first; i < range.second; i++) {
				wrk1[j] += im[j * w + i];//accumulate value of row
				wrk2[i] += im[j * w + i];//accumulate value of colummn
			}
			if(range.first != range.second) wrk1[j] /= range.second - range.first;//normalized to average (without dividing by 0)
		}
		for(size_t i = 0; i < w; i++) if(0 != yCnt->operator[](i)) wrk2[i] /= yCnt->operator[](i);//normalized to average (without dividing by 0)

		//now we have the average of each row/column within the masked region, convert to background with low pass filter
		xFwd->execute(wrk1.data(), wrk3.data());
		std::fill(wrk3.begin() + x, wrk3.end(), 0);
		xInv->execute(wrk3.data(), wrk1.data());
		yFwd->execute(wrk2.data(), wrk3.data());
		std::fill(wrk3.begin() + y, wrk3.end(), 0);
		if(x > 0) wrk3[0] = 0;//don't subtract DC value twice
		yInv->execute(wrk3.data(), wrk2.data());

		//correct for FFTW scaling
		for(size_t j = 0; j < h; j++) wrk1[j] /= h * 2;
		for(size_t i = 0; i < w; i++) wrk2[i] /= w * 2;

		//finally subtract background
		for(size_t j = 0; j < h; j++) {
			const std::pair<size_t, size_t>& range = msk->operator[](j);//get range of good pixels for this row
			for(size_t i = 0           ; i < range.first ; i++) im[j * w + i] = 0;//fill pxiels before mask with zero
			for(size_t i = range.first ; i < range.second; i++) im[j * w + i] -= wrk1[j] + wrk2[i];//subtract background
			for(size_t i = range.second; i < w           ; i++) im[j * w + i] = 0;//fill pixels after mask with zero
		}
	}
}

#endif//_IMAGE_H_
