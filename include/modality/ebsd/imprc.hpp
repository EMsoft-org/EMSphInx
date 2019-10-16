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

#ifndef _EBSD_IM_PRC_
#define _EBSD_IM_PRC_

#include "util/gaussian.hpp"
#include "util/ahe.hpp"

namespace emsphinx {

	namespace ebsd {

		//@brief: class to hold pattern processing details
		template <typename Real>
		class PatternProcessor : public ImageProcessor<Real> {
			size_t                                   nPix ;
			bool                                     doBkg;
			bool                                     doAhe;
			bool                                     msk  ;//is there a mask (or are all pixels good)
			gaussian::BckgSub2D<Real>                bkg  ;
			AdaptiveHistogramEqualizer<Real,uint8_t> ahe  ;
			std::vector<uint8_t>                     work ;//work space for 8 bit conversion (this could be eliminated with some reordering of operations, but it doesn't seem worth it)

			public:
				//@brief: set processing parameters
				//@param w: image width in pixels
				//@param h: image height in pixels
				//@param r: circular mask radius in pixels (-1 for no mask, 0 to automatically size mask)
				//@param b: true to use background subtraction, false otherwise
				//@param n: adaptive histogram equalization nregions
				void setSize(const size_t w, const size_t h, const int r, const bool b, const size_t n);

				//@brief   : process an 8 bit image in place (background subtract and/or AHE)
				//@param im: image to process (modified)
				void process(uint8_t * const im);

				//@brief    : process an image out of place (background subtract and/or AHE)
				//@param im : image to process
				//@param buf: location to write floating point processed image
				template <typename TPix>
				void process(TPix const * const im, Real * const buf);
				void process(uint8_t  const * const im, Real * const buf) {return process<uint8_t >(im, buf);}
				void process(uint16_t const * const im, Real * const buf) {return process<uint16_t>(im, buf);}
				void process(float    const * const im, Real * const buf) {return process<float   >(im, buf);}

				//@brief : get background subtraction mask
				//@return: read only pointer to mask (1/0 for inside/outside mask)
				char const * getMask() const {return bkg.msk->data();}

				//@brief : get size of target image to process
				//@return: size of input image in pixels
				size_t numPix() const {return nPix;}

				//@brief : get a copy of the stored image processor
				//@return: unique pointer to copy of current processor
				std::unique_ptr<ImageProcessor<Real> > clone() const {return std::unique_ptr<PatternProcessor>(new PatternProcessor(*this));}
		};

	}//ebsd

}//emsphinx


#include <algorithm>

namespace emsphinx {

	namespace ebsd {

		//@param w: image width in pixels
		//@param h: image height in pixels
		//@param r: circular mask radius in pixels (-1 for no mask, 0 to automatically size mask)
		//@param b: true to use background subtraction, false otherwise
		//@param n: adaptive histogram equalization nregions
		template <typename Real>
		void PatternProcessor<Real>::setSize(const size_t w, const size_t h, const int r, const bool b, const size_t n) {
			nPix = w * h;

			doBkg = b    ;
			//always build background subtractor (for mask)
			if(-1 == r) {//no mask
				msk = false;
				bkg = gaussian::BckgSub2D<Real>(w, h);
			} else {//circular mask
				msk = true;
				if(r == 0)//automatically determine size
					bkg = gaussian::BckgSub2D<Real>::CircMask((int)w, (int)h);
				else//specify size
					bkg = gaussian::BckgSub2D<Real>::CircMask((int)w, (int)h, r);
			}

			//only build AHE calculator if needed
			doAhe = n > 0;
			if(doAhe) {
				ahe.setSize(w, h, n, n);
			} else {
				ahe = AdaptiveHistogramEqualizer<Real, uint8_t>();
			}

			//allocate 8 bit conversion work space if needed
			if(doBkg) {//background subtraction (+ potentially ahe)
				if(doAhe) {//convert background subtracted floating point image to 8 bit
					work.resize(w * h);
				} else {
					work.clear();
				}
			} else if(doAhe) {//only histogram equalization
				work.resize(w * h);
			} else {//no image processing
				work.clear();
			}


		}

		//@brief   : process an image in place (background subtract and/or AHE)
		//@param im: image to process (modified)
		template <typename Real>
		void PatternProcessor<Real>::process(uint8_t * const im) {
			if(doBkg) {
				bkg.fit(im);//determine background
				bkg.subtract(im);//in place background subtract
			}
			if(doAhe) {
				ahe.equalize(im, msk ? getMask() : NULL);
			}
		}

		//@brief    : process an image out of place (background subtract and/or AHE)
		//@param im : image to process
		//@param buf: location to write floating point processed image
		template <typename Real>
		template <typename TPix>
		void PatternProcessor<Real>::process(TPix const * const im, Real * const buf) {
			if(doBkg) {//background subtraction (+ potentially ahe)
				bkg.fit(im);//determine background
				bkg.subtract(im, buf);//in place background subtract

				if(doAhe) {//convert background subtracted floating point image to 8 bit
					std::pair<Real*, Real*> minMax = std::minmax_element(buf, buf + nPix);
					const Real vMin = *minMax.first;
					const Real factor = Real(255) / (*minMax.second - vMin);
					std::transform(buf, buf + nPix, work.data(), [vMin, factor](const Real& v){return (uint8_t)(factor * (v - vMin) + Real(0.5));});
					ahe.equalize(work.data(), buf, msk ? getMask() : NULL);
				}
			} else if(doAhe) {//only histogram equalization
				if(std::is_same<uint8_t, TPix>::value) {//already 8 bit
					ahe.equalize((uint8_t const *) im, buf, msk ? getMask() : NULL);
				} else {//need to convert input to 8 bit
					std::pair<TPix const*, TPix const*> minMax = std::minmax_element(im, im + nPix);
					const TPix vMin = *minMax.first;
					const Real factor = Real(255) / (*minMax.second - vMin);
					std::transform(im, im + nPix, work.data(), [vMin, factor](const TPix& v){return (uint8_t)(factor * (v - vMin) + Real(0.5));});
					ahe.equalize(work.data(), buf, msk ? getMask() : NULL);
				}
			} else {//no image processing
				std::transform(im, im + nPix, buf, [](const TPix& v){return (Real)v;});//just copy to output
			}
		}

	}//ebsd

}//emsphinx

#endif//_EBSD_IM_PRC_
