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

#ifndef _IDX_BASE_H_
#define _IDX_BASE_H_

#include <memory>
#include <type_traits>
#include <vector>
#include <array>

namespace emsphinx {

	//@brief: abstract base class to provide images to index
	class ImageSource {

		public:
			//possible pixel types
			enum class Bits {
				UNK = 0 ,//unknown
				U8  = 8 ,//uint8_t
				U16 = 16,//uint16_t
				F32 = 32 //float
			};

			//@brief : get the bytes per pixel
			//@return: bytes per pixel
			size_t pixBytes() const;

			//@brief : get the bytes per image
			//@return: bit depth of pixels
			size_t imBytes() const {return pixBytes() * numPix();}

			//@brief : get the pixel type of images
			//@return: bit depth of pixels
			Bits pixelType() const {return bits;}

			//@brief : get the width of images
			//@return: width of images in pixels
			size_t width() const {return w;}

			//@brief : get the width of images
			//@return: height of images in pixels
			size_t height() const {return h;}

			//@brief : get pixels per image
			//@return: pixels per image
			size_t numPix() const {return width() * height();}

			//@brief    : extract the next batch of images into a buffer
			//@param out: location to write extracted images
			//@param cnt: maximum number of images to extract
			//@return   : vector of the indices of each image extracted (e.g. {0,2,1,3} for the first 4 images but out of order)
			//@note     : implementations should be thread safe (only for other extract calls)
			virtual std::vector<size_t> extract(char * const out, const size_t cnt) const = 0;

		protected:
			Bits               bits;//pixel type
			size_t             w   ;//width of image
			size_t             h   ;//height of image
	};

	//@brief: abstract base class to preprocess images
	template <typename Real>
	class ImageProcessor {
		static_assert(std::is_floating_point<Real>::value, "Real must be a floating point type");

		public:

			//@brief    : process an image out of place (e.g. do background subtraction)
			//@param im : image to process
			//@param buf: location to write floating point processed image
			virtual void process(uint8_t  const * const im, Real * const buf) = 0;
			virtual void process(uint16_t const * const im, Real * const buf) = 0;
			virtual void process(float    const * const im, Real * const buf) = 0;

			//@brief : get size of target image to process
			//@return: size of input image in pixels
			virtual size_t numPix() const = 0;

			//@brief : get a copy of the stored image processor
			//@return: unique pointer to copy of current processor
			virtual std::unique_ptr<ImageProcessor> clone() const = 0;

			//@brief: default destructor (for unique_ptr)
			virtual ~ImageProcessor() = default;
	};

	//@brief: abstract bass class to encapsulate back projection of images to the sphere
	template <typename Real>
	class BackProjector {
		public:
			//@brief    : unproject an image from the detector to a square legendre grid
			//@param im : processed image to back project to unit sphere
			//@param sph: location to write back projected image
			//@param iq : location to write image image quality (NULL to skip computation)
			//@return   : sqrt(integral of im^2) (properly weighted on the sphere)
			virtual Real unproject(Real * const im, Real * const sph, Real * const iq = NULL) = 0;

			//@brief: compute quaternion to rotate back projection such that is centered around the north pole
			//@return: quaternion as wxyz
			virtual std::array<Real, 4> northPoleQuat() const {return std::array<Real, 4>({1, 0, 0, 0});}

			//@brief : get a copy of the stored back projector
			//@return: unique pointer to copy of current projector
			virtual std::unique_ptr<BackProjector> clone() const = 0;

			//@brief: default destructor (for unique_ptr)
			virtual ~BackProjector() = default;
	};
}

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

namespace emsphinx {

	//@brief : get the bytes per pixel
	//@return: bytes per pixel
	size_t ImageSource::pixBytes() const {
		switch(pixelType()) {
			case Bits::UNK: return 0;
			case Bits::U8 : return 1;
			case Bits::U16: return 2;
			case Bits::F32: return 4;
		}
		return 0;
	}
}

#endif//_IDX_BASE_H_

