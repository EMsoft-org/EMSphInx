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

#ifndef _AHE_H_
#define _AHE_H_

#include <vector>
#include <type_traits>

//@brief: helper to repeatedly do AHE for the same image size / tiles
//@note : this should theorhetically work for e.g. uint16_t or uint32_t, but be careful since a single 16 bit CDF is 512 kB and a single 32bit CDF is 32 GB
template <typename Real, typename TPix = uint8_t>
class AdaptiveHistogramEqualizer {
	static_assert(std::is_floating_point<Real>::value, "Real must be floating point type");
	static_assert(std::is_integral<TPix>::value, "TPix must be unsigned integer type");
	static_assert(std::is_unsigned<TPix>::value, "TPix must be unsigned integer type");
	static_assert(std::numeric_limits<size_t>::max() / std::numeric_limits<TPix>::max() > 64, "size_t not large enough to hold 64 TPix histograms");

	static const size_t HistBins = size_t(std::numeric_limits<TPix>::max()) + 1;//number of histogram bins

	//@brief: helper to hold tile bounds
	struct TileBounds {
		size_t iS, iE;//x start and end
		size_t jS, jE;//y start and end
	};

	//@brief: helper to hold pixel interpolation values
	struct InterpPair {
		size_t l, u;//lower and upper bounding tile indices
		Real c, f;//interpolation weights (sum to 1)
	};

	//these values are read only for equalization (can be shared across threads)
	std::vector< TileBounds > tiles;//bounds of tiles
	std::vector< InterpPair > jInds, iInds;//interpolation coefficients for each row/column

	//these members are working space for equalization (can't be shared across threads)
	std::vector<Real> cdfs;

	//@brief    : compute the patch histograms
	//@param im : image to compute histograms from (width and height are assumed to match prior call to setSize)
	//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
	void computeHist(TPix const* im, char const * const msk = NULL);

	public:
		//@brief   : set the image size and tile count for the equalizer
		//@param w : width of image in pixels
		//@param h : height of image in pixels
		//@param nx: number of x grid points
		//@param ny: number of y grid points
		//@note    : equalization complexity will scale as nx * ny for high grid densities
		void setSize(const size_t w, const size_t h, const size_t nx, const size_t ny);

		//@brief    : equalize an image in place using previously set conditions
		//@param im : image (width and height are assumed to match prior call to setSize)
		//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
		void equalize(TPix* im, char const * const msk = NULL);

		//@brief    : equalize an image out of place using previously set conditions
		//@param im : image (width and height are assumed to match prior call to setSize)
		//@param buf: location to write equalized image
		//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
		void equalize(TPix const * im, Real * buf, char const * const msk = NULL);
};

//@brief   : apply adaptive histogram equalization to a single image
//@param im: image to equalize
//@param w : width of image in pixels
//@param h : height of image in pixels
//@param nx: number of x grid points
//@param ny: number of y grid points
//@note    : this is a convenience function to build and use an AdaptiveHistogramEqualizer object so the performance won't be great
void adHistEq(uint8_t* im, const size_t w, const size_t h, const size_t nx, const size_t ny);

#include <cmath>
#include <algorithm>
#include <numeric>

//@brief   : set the image size and tile count for the equalizer
//@param w : width of image in pixels
//@param h : height of image in pixels
//@param nx: number of x grid points
//@param ny: number of y grid points
//@note    : equalization complexity will scale as nx * ny for high grid densities
template <typename Real, typename TPix>
void AdaptiveHistogramEqualizer<Real, TPix>::setSize(const size_t w, const size_t h, const size_t nx, const size_t ny) {
	tiles.resize(nx * ny);
	cdfs .resize(nx * ny * HistBins);

	//compute tile bounds		
	const Real tx = Real(w) / nx;//compute tile width in fractional pixels
	const Real ty = Real(h) / ny;//compute tile height in fractional tiles
	const Real hWdth = Real(0.5);//histogram calculation area, 1.0 for 50% overlap, 0.5 for mosaic (1.0 does background removal better, 0.5 does enhancement better if there isn't a background)
	std::vector<size_t> jMids(ny), iMids(nx);//tile midpoints
	for(size_t j = 0; j < ny; j++) {//loop over tile rows
		//compute top, middle, and bottom of tile in pixels
		Real midY = ty * j + ty / 2;
		Real minY = std::max(midY - ty * hWdth, Real(0));
		Real maxY = std::min(midY + ty * hWdth, Real(h));
		minY = std::round(minY); midY = std::round(midY); maxY = std::round(maxY);
		jMids[j] = (size_t)midY;
		for(size_t i = 0; i < nx; i++) {//loop over tile cols
			//compute left, middle, and right of tile in pixels
			Real midX = tx * i + tx / 2;
			Real minX = std::max(midX - tx * hWdth, Real(0));
			Real maxX = std::min(midX + tx * hWdth, Real(w));
			minX = std::round(minX); midX = std::round(midX); maxX = std::round(maxX);
			if(j == 0) iMids[i] = (size_t)midX;

			//save tile bounds
			TileBounds& t = tiles[j * nx + i];
			t.iS = (size_t)minX; t.iE = (size_t)maxX;
			t.jS = (size_t)minY; t.jE = (size_t)maxY;
		}
	}

	//compute y interpolation once
	jInds.resize(w);
	for(size_t j = 0; j < h; j++) {
		const size_t u = std::distance(jMids.cbegin(), std::upper_bound(jMids.cbegin(), jMids.cend(), j));
		if(jMids.size() == u) {//beyond last grid point
			jInds[j].l = jInds[j].u = jMids.size() - 1;
			jInds[j].c = jInds[j].f = Real(0.5);
		} else if(0 == u) {//before first grid point
			jInds[j].l = jInds[j].u = 0;
			jInds[j].c = jInds[j].f = Real(0.5);
		} else {//between 2 grid points, linear interpolate
			jInds[j].l = u - 1;
			jInds[j].u = u;
			jInds[j].f = Real(j - jMids[u - 1]) / (jMids[u] - jMids[u - 1]);
			jInds[j].c = Real(1) - jInds[j].f;
		}
		jInds[j].l *= nx * HistBins;//convert from tile to vectorized histograms index
		jInds[j].u *= nx * HistBins;//convert from tile to vectorized histograms index
	}

	//compute x interpolation once
	iInds.resize(h);
	for(size_t i = 0; i < w; i++) {
		const size_t u = std::distance(iMids.cbegin(), std::upper_bound(iMids.cbegin(), iMids.cend(), i));
		if(iMids.size() == u) {//beyond last grid point
			iInds[i].l = iInds[i].u = iMids.size() - 1;
			iInds[i].c = iInds[i].f = Real(0.5);
		} else if(0 == u) {//before first grid point
			iInds[i].l = iInds[i].u = 0;
			iInds[i].c = iInds[i].f = Real(0.5);
		} else {//between 2 grid points, linear interpolate
			iInds[i].l = u - 1;
			iInds[i].u = u;
			iInds[i].f = Real(i - iMids[u - 1]) / (iMids[u] - iMids[u - 1]);
			iInds[i].c = Real(1) - iInds[i].f;
		}
		iInds[i].l *= HistBins;//convert from tile to vectorized histograms index
		iInds[i].u *= HistBins;//convert from tile to vectorized histograms index
	}
}

//@brief    : compute the patch histograms
//@param im : image to compute histograms from (width and height are assumed to match prior call to setSize)
//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
template <typename Real, typename TPix>
void AdaptiveHistogramEqualizer<Real, TPix>::computeHist(TPix const* im, char const * const msk) {
	//compute CDF of each tile up front
	size_t hist[HistBins];
	for(size_t i = 0; i < tiles.size(); i++) {
		//compute histogram
		TileBounds& t = tiles[i];
		std::fill(hist, hist + HistBins, 0);
		if(NULL == msk) {
			for(size_t j = t.jS; j < t.jE; j++) {//loop over image rows of tile t
				for(size_t i = t.iS; i < t.iE; i++) ++hist[im[iInds.size() * j + i]];//loop over image columns of tile t
			}
		} else {
			for(size_t j = t.jS; j < t.jE; j++) {//loop over image rows of tile t
				for(size_t i = t.iS; i < t.iE; i++) {
					const size_t idx = iInds.size() * j + i;
					if(1 == msk[idx]) ++hist[im[idx]];//loop over image columns of tile t
				}
			}
			if(0 == *std::max_element(hist, hist + HistBins)) {//there were no good pixels
				std::fill(hist, hist + HistBins, 1);//fill with linear ramp (no adjustment)
			}
		}

		//convert to normalized cumulative distribution
		std::partial_sum(hist, hist + HistBins, hist);//unnormalized CDF
		const Real nrm = Real(HistBins-1) / hist[HistBins-1];//normalization so integral of CDF is pixel max
		std::transform(hist, hist + HistBins, cdfs.begin() + i * HistBins, [nrm](const size_t& v){return nrm * v;});
	}
}

//@brief    : equalize an image using previously set conditions
//@param im : image (width and height are assumed to match prior call to setSize)
//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
template <typename Real, typename TPix>
void AdaptiveHistogramEqualizer<Real, TPix>::equalize(TPix* im, char const * const msk) {
	computeHist(im, msk);

	//loop over pixels equalizing
	for(const InterpPair& j : jInds) {
		for(const InterpPair& i : iInds) {
			const TPix v = *im;
			*im++ = (TPix) ( cdfs[(j.l + i.l) + v] * j.c * i.c
			               + cdfs[(j.l + i.u) + v] * j.c * i.f
			               + cdfs[(j.u + i.l) + v] * j.f * i.c
			               + cdfs[(j.u + i.u) + v] * j.f * i.f
			               + Real(0.5));//add 0.5 so that casting to uint is rounding
		}
	}
}

//@brief    : equalize an image using previously set conditions
//@param im : image (width and height are assumed to match prior call to setSize)
//@param msk: mask of valid pixels (1/0 for valid/invalid) or NULL for all valid pixels
template <typename Real, typename TPix>
void AdaptiveHistogramEqualizer<Real, TPix>::equalize(TPix const * im, Real * buf, char const * const msk) {
	computeHist(im, msk);

	//loop over pixels equalizing
	for(const InterpPair& j : jInds) {
		for(const InterpPair& i : iInds) {
			const TPix v = *im++;
			*buf++ = cdfs[(j.l + i.l) + v] * j.c * i.c
			       + cdfs[(j.l + i.u) + v] * j.c * i.f
			       + cdfs[(j.u + i.l) + v] * j.f * i.c
			       + cdfs[(j.u + i.u) + v] * j.f * i.f;
		}
	}
}

//@brief   : apply adaptive histogram equalization to a single image
//@param im: image to equalize
//@param w : width of image in pixels
//@param h : height of image in pixels
//@param nx: number of x grid points
//@param ny: number of y grid points
//@note    : this is a convenience function to build and use an AdaptiveHistogramEqualizer object so the performance won't be great
void adHistEq(uint8_t* im, const size_t w, const size_t h, const size_t nx, const size_t ny) {
	AdaptiveHistogramEqualizer<double, uint8_t> eq;
	eq.setSize(w, h, nx, ny);
	eq.equalize(im);
}

#endif//_AHE_H_
