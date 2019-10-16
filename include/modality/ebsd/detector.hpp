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

#ifndef _DETECTOR_H_
#define _DETECTOR_H_

#include <algorithm>
#include <memory>
#include <array>

#include "H5Cpp.h"

#include "util/fft.hpp"
#include "util/image.hpp"
#include "idx/base.hpp"

namespace emsphinx {

	namespace ebsd {
		//@brief: helper for describing the geometry of a detector
		template <typename Real>
		struct Geometry {
			//parameters to define the geometry of an EBSD detector
			Real   dTlt, dOmg;//detector tilt/rotation in degrees
			Real   sTlt, sOmg;//sample tilt/rotation in degrees degrees
			size_t w   , h   ;//camera width/height in pixels
			Real   pX  , pY  ;//pixel width/height in microns
			Real   cX  , cY  ;//x/y pattern center position relative to center of image (in pixels)
			Real   sDst      ;//scintillator distance in microns
			bool   circ      ;//true / false if a circular mask should/shouldn't be used
			bool   flip      ;//true / false if the images are vertically flipped (image coordinate system / cartesian coordinat system) [origin top left/origin bottom left]

			//@brief: construct an empty geometry
			Geometry() : dTlt(0), dOmg(0), sTlt(0), sOmg(0), w(0), h(0), pX(0), pY(0), cX(0), cY(0), sDst(0), circ(false), flip(false) {}

			//@brief    : set the sample tilt
			//@param sig: sample tilt in degrees (e.g. 70)
			void sampleTilt(const Real sig) {sTlt = sig;}

			//@brief    : set the camera tilt
			//@param sig: camera tilt in degrees (e.g. 15)
			void cameraTilt(const Real thetac) {dTlt = thetac;}

			//@brief      : set camera dimensions using EMsoft parameters
			//@param numsx: width  of detector in pixels
			//@param numsy: height of detector in pixels
			//@param delta: pixel size in microns
			void cameraSize(const size_t numsx, const size_t numsy, const Real delta) {w = numsx; h = numsy; pX = pY = delta;}

			//@brief    : set pattern center using EMsoft parameters
			//@param xpc: x pattern center in fractional pixels relative to center of image
			//@param ypc: y pattern center in fractional pixels relative to center of image
			//@param L  : scintillator distance in microns
			void patternCenter(const Real xpc, const Real ypc, const Real L) {cX = xpc; cY = ypc; sDst = L;}

			//@brief      : set pattern center using tsl parameters
			//@param xStar: tsl x^*
			//@param yStar: tsl y^*
			//@param zStar: tsl z^*
			//@param note : camera and pixel size need to be set first
			void patternCenterTSL(const Real xStar, const Real yStar, const Real zStar);

			//@brief      : set pattern center using oxford parameters
			//@param xStar: oxford x^*
			//@param yStar: oxford y^*
			//@param zStar: oxford z^*
			//@param note : camera and pixel size need to be set first
			void patternCenterOxford(const Real xStar, const Real yStar, const Real zStar);

			//@brief      : set pattern center using bruker parameters
			//@param xStar: bruker x^*
			//@param yStar: bruker y^*
			//@param zStar: bruker z^*
			//@param note : camera and pixel size need to be set first
			void patternCenterBruker(const Real xStar, const Real yStar, const Real zStar);

			//@brief      : build a detector using electron channeling pattern parameters
			//@param dim  : size (sidelength) of ECP in pixels
			//@param theta: maximum procession angle in degrees
			void ecp(const size_t dim, const Real theta);

			//@brief     : set if a circular mask should be used
			//@param mask: true/false if a circular mask should/shouldn't be used
			void maskPattern(const bool mask) {circ = mask;}

			//@brief         : set if detector images are vertically flipped
			//@param flipVert: true/false detector images are/aren't vertically flipped
			void flipPattern(const bool flipVert) {flip = flipVert;}

			//@brief  : rebin a detector (group pixels into nxn super pixel)
			//@param n: side length of block of pixels to merge
			//@return : binned geometry
			Geometry bin(const size_t n) const;

			//@brief    : bilinearlly interpolate a pixel value for the given direction
			//@param n  : unit direction to interpolate value for in the sample reference frame (pattern center is in +z direction)
			//@param pix: [optional] location to store interpolation pixel
			//@param flp: [optional] is the image to interpolate from flipped vertically
			//@return   : true / false if the direction lies inside / outside the detector
			bool interpolatePixel(Real const * const n, image::BiPix<Real> * const pix = NULL, const bool flp = false) const;

			//@brief  : get the sample direction the points toward the specified pixel
			//@param X: fractional x position on detector (relative to detector center)
			//@param Y: fractional y position on detector (relative to detector center)
			//@param n: location to write unit direction pointing to (X,Y) in sample reference frame
			void sampleDir(const Real X, const Real Y, Real * const n) const;

			//@brief        : approximate the solid angle captured by a detector
			//@param gridRes: side length of square lambert projection to grid over
			//@return       : fractional solid angle i.e [0,1] not [0,4*pi]
			Real solidAngle(size_t gridRes) const;

			//@brief      : create a new detector geometry by rescaling the pixels size while maintaining solid angle
			//@param scale: rescaling factor, must lie in [0, min(w, h)] with values < 1 corresponding to increasing pixel density
			//@note       : this is analogous to continuous camera binning with e.g. scale = 4.0 corresponding to 4x4 camera binning
			Geometry<Real> rescale(const Real scale) const;

			//@brief     : create a new detector geometry by rescaling the pixels size while maintaining solid angle
			//@param wNew: scaled detector width
			//@param hNew: scaled detector height
			Geometry<Real> rescale(size_t wNew, const size_t hNew) const;

			//@brief: compute quaternion to rotate detector such that is centered around the north pole
			//@return: quaternion as wxyz
			std::array<Real, 4> northPoleQuat() const;

			//@brief    : compute the rescale factor required to make average detector pixel size the same as an average pixel on a square spherical grid
			//@param dim: side length of square spherical grid to match pixel size to
			//@return   : geometry scale factor to make average pixel sizes the same
			Real scaleFactor(const size_t dim) const;

			//@brief    : read (partial) detector geometry info from an EMsoft master pattern file
			//@param grp: folder in h5 file to read from (e.g. "/NMLparameters/EBSDNameList/")
			//@note     : this doesn't get the sample tilt which is a monte carlo parameter
			void readEMsoft(H5::Group grp);
		};

		//helper for back projecting from an EBSD detector -> sphere
		template <typename Real>
		struct BackProjector : public emsphinx::BackProjector<Real> {
			struct Constants;//helper struct to hold read only constants
			const std::shared_ptr<const Constants> pLut;//read only values (can be shared across threads)
			std::vector<Real>                      rPat;//space for pattern as real
			fft::vector<Real>                      sWrk;//work space for pattern rescaler
			std::vector<Real>                      iVal;//space for intermediate interpolation result

			//@brief    : construct a back projector
			//@param geo: detector geometry to back project from
			//@param dim: side length of square legendre grid to back project onto
			//@param fct: additional scale factor for detector resizing
			//@param qu : quaternion to rotate detector location by
			BackProjector(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu = NULL);

			//@brief    : unproject a pattern from the detector to a square grid
			//@param pat: pattern to back project to unit sphere
			//@param sph: location to write back projected pattern
			//@param iq : location to write pattern image quality (NULL to skip computation)
			//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
			Real unproject(Real * const pat, Real * const sph, Real * const iq = NULL);

			//@brief : get a copy of the stored back projector
			//@return: unique pointer to copy of current projector
			std::unique_ptr<emsphinx::BackProjector<Real> > clone() const {return std::unique_ptr<BackProjector>(new BackProjector(*this));}

			//@brief    : build a mask of spherical pixels covered by the detector
			//@param sph: location to write mask
			void mask(Real * const sph);
		};

		//helper struct to hold read only constants needed by BackProjector (for sharing across threads)
		template <typename Real>
		struct BackProjector<Real>::Constants {
			std::vector< image::BiPix<Real> > iPts;//interpolation coefficients
			std::vector<Real>                 omeg;//solid angles of iPts (relative to average spherical pixel size)
			Real                              omgW;//sum of omeg (in window)
			Real                              omgS;//sum of omeg (over sphere)
			image::Rescaler<Real>             sclr;//pattern rescaler
			const bool                        flip;//do input patterns need to be vertically flipped

			//@brief    : construct a back projector
			//@param geo: detector geometry to back project from
			//@param dim: side length of square legendre grid to back project onto
			//@param fct: additional scale factor for detector resizing
			//@param qu : quaternion to rotate detector location by
			Constants(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu = NULL);
		};
	}//ebsd

}//emsphinx

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <stdexcept>
#include <algorithm>

#include "constants.hpp"
#include "sht/square_sht.hpp"
#include "sht/sht_xcorr.hpp"
#include "xtal/rotations.hpp"
#include "xtal/quaternion.hpp"

namespace emsphinx {

	namespace ebsd {

		////////////////////////////////////////////////////////////////////////
		//                              Geometry                              //
		////////////////////////////////////////////////////////////////////////

		//@brief      : set pattern center using tsl parameters
		//@param xStar: tsl x^*
		//@param yStar: tsl y^*
		//@param zStar: tsl z^*
		//@param note : camera and pixel size need to be set first
		template <typename Real>
		void Geometry<Real>::patternCenterTSL(const Real xStar, const Real yStar, const Real zStar) {
			cX = xStar * w - Real(w) / 2;
			cY = yStar * w - Real(h) / 2;
			sDst = zStar * w * pX;
		}

		//@brief      : set pattern center using oxford parameters
		//@param xStar: oxford x^*
		//@param yStar: oxford y^*
		//@param zStar: oxford z^*
		//@param note : camera and pixel size need to be set first
		template <typename Real>
		void Geometry<Real>::patternCenterOxford(const Real xStar, const Real yStar, const Real zStar) {
			cX = (xStar - Real(0.5)) * w;
			cY = (yStar - Real(0.5)) * h;
			sDst = zStar * w * pX;

		}

		//@brief      : set pattern center using bruker parameters
		//@param xStar: bruker x^*
		//@param yStar: bruker y^*
		//@param zStar: bruker z^*
		//@param note : camera and pixel size need to be set first
		template <typename Real>
		void Geometry<Real>::patternCenterBruker(const Real xStar, const Real yStar, const Real zStar) {
			cX = (xStar - Real(0.5)) * w;
			cY = (Real(0.5) - yStar) * h;
			sDst = zStar * h * pX;
		}

		//@brief      : build a detector using electron channeling pattern parameters
		//@param dim  : size (sidelength) of ECP in pixels
		//@param theta: maximum procession angle in degrees
		template <typename Real>
		void Geometry<Real>::ecp(const size_t dim, const Real theta) {
			sTlt = sOmg = 0;//no sample tilts
			dTlt = -90;//detector directly above sample
			dOmg = 0;//no detector rotation
			cX = cY = 0;//no pattern center adjustment
			circ = true ;//maximum procession angle -> circular patterns
			flip = false;//we shouldn't need to flip the pattern
			w = h = dim;//save effective detector size
			sDst = 10000;//select an arbitrary detector distance (I've used a common working distance)
			
			//compute effective pixel size
			const Real dTheta = theta * emsphinx::Constants<Real>::pi / (dim * 90);//what is the size of a single pixel in radians
			pX = pY = std::tan(dTheta) * sDst;//compute the size of a pixel compared to the arbitrary detector distance
		}

		//@brief  : rebin a detector (group pixels into nxn super pixel)
		//@param n: side length of block of pixels to merge
		//@return : binned geometry
		template <typename Real>
		Geometry<Real> Geometry<Real>::bin(const size_t n) const {
			if(n == 0) throw std::runtime_error("binning number must be at least 1");
			if(n == 1) return *this;
			if(w % n != 0 || h % n != 0) throw std::runtime_error("binning factor must be a divisor of both height and width");
			Geometry geo(*this);//copy everything
			Real bW = Real(w) / n;//compute fractional super pixel width
			Real bH = Real(h) / n;//compute fractional super pixel height
			geo.w = (size_t) std::ceil(bW);//adjust detector size
			geo.h = (size_t) std::ceil(bH);//adjust detector size
			geo.pX *= n;//adjust pixel width
			geo.pY *= n;//adjust pixel height
			geo.cX /= n;//adjust pixel width
			geo.cY /= n;//adjust pixel height
			/*
				//handling non exact divisor binning by converting pattern center to be relative to origin before scaling
				//this assumes that the first superpixel is at (0,0) and extra subpixels are cropped off
				const Real aX = cX + Real(0.5) * w;//pattern center x in fractional pixels relative to origin
				const Real aY = cY + Real(0.5) * h;//pattern center y in fractional pixels relative to origin
				geo.cX = (aX / n) - Real(0.5) * geo.w;//scale pattern center x relative to origin and the update back to relative to center
				geo.cY = (aY / n) - Real(0.5) * geo.h;//scale pattern center y relative to origin and the update back to relative to center
			*/
			return geo;
		}

		//@brief    : bilinearlly interpolate a pixel value for the given direction
		//@param n  : unit direction to interpolate value for in the sample reference frame (pattern center is in +z direction)
		//@param pix: [optional] location to store interpolation pixel
		//@param flp: [optional] is the image to interpolate from flipped vertically
		//@return   : true / false if the direction lies inside / outside the detector
		template <typename Real>
		bool Geometry<Real>::interpolatePixel(Real const * const n, image::BiPix<Real> * const pix, const bool flp) const {
			//sanity check
			if(std::signbit(n[2])) return false;//can't project through sample
			if(0 != dOmg || 0 != sOmg) throw std::runtime_error("omega tilt not yet supported");

			//compute angle between detector and sample normal
			static const Real degToRad = emsphinx::Constants<Real>::pi / 180;
			const Real alpha = (Real(90) - sTlt + dTlt) * degToRad;//angle between sample and detector is 90 with no tilts, increasing sample tilt decreases the angle, increasing camera tilt increases the angle
			if(std::fabs(alpha) > emsphinx::Constants<Real>::pi_2) throw std::runtime_error("pattern center not on detector");//beyond +/-90 the detector is behind the sample

			//compute sin / cos of angle and distance to detector
			const Real sA = std::sin(alpha);
			const Real cA = std::sqrt(Real(1) - sA * sA);//cos(alpha)
			const Real d =  sDst / (n[0] * sA + n[2] * cA);
			if(std::signbit(d)) return false;//negative distance

			//compute offset from pattern center in microns
			const Real x =                   n[1]  * d;
			const Real y = (sA * n[2] - cA * n[0]) * d;

			//convert from microns to fractional position (0->100%) and check if we're inside the detector
			const Real X = (cX + x / pX) / w + Real(0.5);
			      Real Y = (cY + y / pY) / h + Real(0.5);

			// check if we're inside the detector
			if(std::signbit(X) || std::signbit(Y) || X > 1 || Y > 1) return false;//outside of frame
			if(flp) Y = Real(1) - Y;

			//check against circular mask if needed
			if(circ) {
				const Real dX = (X - Real(0.5)) * w;//horizontal distance from center in pixels
				const Real dY = (Y - Real(0.5)) * h;//vertical distance from center in pixels
				const Real r = Real(std::min(w, h)) / 2;//radius of circular mask
				if(r * r < dX * dX + dY * dY) return false;
			}

			//if we made it this far we're in frame, compute relative contribution from neighboring pixels
			if(pix != NULL) pix->bilinearCoeff(X, Y, w, h);
			return true;
		}

		//@brief  : get the sample direction the points toward the specified pixel
		//@param X: fractional x position on detector (relative to detector center)
		//@param Y: fractional y position on detector (relative to detector center)
		//@param n: location to write unit direction pointing to (X,Y) in sample reference frame
		template <typename Real>
		void Geometry<Real>::sampleDir(const Real X, const Real Y, Real * const n) const {
			//compute angle between detector/sample and sin/cos
			const Real alpha = (Real(90) - sTlt + dTlt) * xtal::Constants<Real>::dg2rd;//angle between sample and detector is 90 with no tilts, increasing sample tilt decreases the angle, increasing camera tilt increases the angle
			const Real sA = std::sin(alpha);
			const Real cA = std::sqrt(Real(1) - sA * sA);//cos(alpha)

			//compute denomenator
			const Real fx = (cX - X * w) * pX;
			const Real fy = (cY - Y * h) * pY;
			const Real den = std::sqrt(sDst * sDst + fx * fx + fy * fy);

			//compute normal
			n[0] = (sDst * sA + fy * cA) / den;
			n[1] =            - fx       / den;
			n[2] = (sDst * cA - fy * sA) / den;
		}

		//@brief        : approximate the solid angle captured by a detector
		//@param gridRes: side length of square lambert projection to grid over
		//@return       : fractional solid angle i.e [0,1] not [0,4*pi]
		template <typename Real>
		Real Geometry<Real>::solidAngle(size_t gridRes) const {
			//loop over northern hemisphere
			Real XY[2], xyz[3];
			size_t goodPoints = 0;
			for(size_t j = 0; j <= gridRes; j++) {//loop over rows of square projection
				XY[1] = Real(j) / gridRes;
				for(size_t i = 0; i <= gridRes; i++) {//loop over columns of square projection
					XY[0] = Real(i) / gridRes;
					square::lambert::squareToSphere(XY[0], XY[1], xyz[0], xyz[1], xyz[2]);
					if(interpolatePixel(xyz)) ++goodPoints;//count this point if it falls on the detector
				}
			}
			const size_t totalPoints = (gridRes * gridRes + (gridRes - 2) * (gridRes - 2));//total number of points in both hemispheres (don't double count equators)
			return Real(goodPoints) / totalPoints;//count fraction of grid points on detector
		}

		//@brief      : create a new detector geometry by rescaling the pixels size while maintaining solid angle
		//@param scale: rescaling factor, must lie in [0, min(w, h)] with values < 1 corresponding to increasing pixel density
		//@note       : this is analogous to continuous camera binning with e.g. scale = 4.0 corresponding to 4x4 camera binning
		template <typename Real>
		Geometry<Real> Geometry<Real>::rescale(const Real scale) const {
			const size_t wNew = (size_t) std::round(Real(w) / scale);
			const size_t hNew = (size_t)std::round(Real(h) / scale);
			if(wNew == 0 || hNew == 0) throw std::runtime_error("cannot rescale detector to less than 1 pixel");
			return rescale(wNew, wNew);
		}

		//@brief     : create a new detector geometry by rescaling the pixels size while maintaining solid angle
		//@param wNew: scaled detector width
		//@param hNew: scaled detector height
		template <typename Real>
		Geometry<Real> Geometry<Real>::rescale(size_t wNew, const size_t hNew) const {
			//copy everything and update size
			Geometry<Real> geom(*this);
			geom.w = wNew;
			geom.h = hNew;

			//compute scale factor in x/y direction
			const Real sx = Real(w) / geom.w;
			const Real sy = Real(h) / geom.h;

			//rescale pixels so that the detector covers the same area
			geom.pX *= sx;
			geom.pY *= sy;

			//rescale the pattern center to account for the new pixel size
			geom.cX /= sx;
			geom.cY /= sy;
			return geom;
		}

		//@brief: compute quaternion to rotate detector such that is centered around the north pole
		//@return: quaternion as wxyz
		template <typename Real>
		std::array<Real, 4> Geometry<Real>::northPoleQuat() const {
			const Real theta = (Real(90) - sTlt + dTlt) * emsphinx::Constants<Real>::pi / 180;
			// return std::array<Real, 4>({std::cos(theta/2), 0, std::sin(theta/2) * emsphinx::pijk, 0});
			return std::array<Real, 4>({1, 0, 0, 0});
		}

		//@brief    : compute the rescale factor required to make average detector pixel size the same as an average pixel on a square spherical grid
		//@param dim: side length of square spherical grid to match pixel size to
		//@return   : geometry scale factor to make average pixel sizes the same
		template <typename Real>
		Real Geometry<Real>::scaleFactor(const size_t dim) const {
			const size_t sqrPix = dim * dim * 2 - (dim - 1) * 4;//number of pixels in square grid
			const size_t detPix = w * h;//number of pixels on detector
			return std::sqrt(solidAngle(501) * sqrPix / detPix);
		}

		//@brief    : read (partial) detector geometry info from an EMsoft master pattern file
		//@param grp: folder in h5 file to read from (e.g. "/NMLparameters/EBSDNameList/")
		//@note     : this doesn't get the sample tilt which is a monte carlo parameter
		template <typename Real>
		void Geometry<Real>::readEMsoft(H5::Group grp) {
			float thetac, delta, xpc, ypc, L;
			int32_t numsx, numsy, binning;
			std::string maskpattern;
			grp.openDataSet("thetac"     ).read(&thetac     , H5::PredType::NATIVE_FLOAT);
			grp.openDataSet("delta"      ).read(&delta      , H5::PredType::NATIVE_FLOAT);
			grp.openDataSet("xpc"        ).read(&xpc        , H5::PredType::NATIVE_FLOAT);
			grp.openDataSet("ypc"        ).read(&ypc        , H5::PredType::NATIVE_FLOAT);
			grp.openDataSet("L"          ).read(&L          , H5::PredType::NATIVE_FLOAT);
			grp.openDataSet("numsx"      ).read(&numsx      , H5::PredType::NATIVE_INT32);
			grp.openDataSet("numsy"      ).read(&numsy      , H5::PredType::NATIVE_INT32);
			grp.openDataSet("binning"    ).read(&binning    , H5::PredType::NATIVE_INT32);
			grp.openDataSet("maskpattern").read( maskpattern, H5::StrType(0, H5T_VARIABLE), H5::DataSpace(H5S_SCALAR));
			cameraTilt(thetac);
			cameraSize(numsx, numsy, delta);
			patternCenter(xpc, ypc, L);
			if(1 != binning) *this = bin(binning);
			if     (0 == maskpattern.compare("y")) maskPattern(true );
			else if(0 == maskpattern.compare("n")) maskPattern(false);
			else throw std::runtime_error("unknown maskpattern option '" + maskpattern + "'");
		}

		////////////////////////////////////////////////////////////////////////
		//                           BackProjector                            //
		////////////////////////////////////////////////////////////////////////

		template <typename Real>
		BackProjector<Real>::Constants::Constants(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu) :
			sclr(geo.w, geo.h, geo.scaleFactor(dim) * fct, fft::flag::Plan::Patient),//we're going to be doing lots of 
			flip(geo.flip)
		{
			//compute normals and solid angles of square legendre grid
			std::vector<Real> xyz(dim * dim * 3);
			square::legendre::normals(dim, xyz.data());
			std::vector<Real> omegaRing = square::solidAngles<Real>(dim, square::Layout::Legendre);

			//rescale detector
			Geometry<Real> g = geo.rescale(sclr.wOut, sclr.hOut);

			//determine weights north hemisphere
			image::BiPix<Real> p;
			for(size_t i = 0; i < xyz.size() / 3; i++) {//loop over square legendre grid points
				//get direction and apply rotatation
				Real n[3];
				if(NULL != qu) {
					xtal::quat::rotateVector(qu, xyz.data() + 3 * i, n);
				} else {
					std::copy(xyz.data() + 3 * i, xyz.data() + 3 * i + 3, n);
				}

				//save interpolation weights if it hits the detector
				if(g.interpolatePixel(n, &p, flip)) {
					p.idx = i;
					iPts.push_back(p);
					omeg.push_back(omegaRing[square::ringNum(dim, i)]);
				}
			}

			//determine weights south hemisphere
			for(size_t i = 0; i < xyz.size() / 3; i++) {//loop over square legendre grid points
				//split into x and y indices
				const size_t y = i / dim;
				if(y == 0 || y+1 == dim) continue;//dont' double count equator
				const size_t x = i - dim * y;
				if(x == 0 || x+1 == dim) continue;//dont' double count equator
				xyz[3*i+2] = -xyz[3*i+2];//move normal to southern hemisphere

				//get direction and apply rotatation
				Real n[3];
				if(NULL != qu) {
					xtal::quat::rotateVector(qu, xyz.data() + 3 * i, n);
				} else {
					std::copy(xyz.data() + 3 * i, xyz.data() + 3 * i + 3, n);
				}

				//save interpolation weights if it hits the detector
				if(g.interpolatePixel(n, &p, flip)) {
					p.idx = i;
					iPts.push_back(p);
					omeg.push_back(omegaRing[square::ringNum(dim, i)]);
				}
			}

			//accumulate solid angle sums
			omgW = std::accumulate(omeg.cbegin(), omeg.cend(), Real(0));
			omgS = 0;
			for(size_t j = 0; j < dim; j++) {
				for(size_t i = 0; i < dim; i++) {
					const Real& o = omegaRing[square::ringNum(dim, j * dim + i)];
					const bool eq = j == 0 || i == 0 || j == dim-1 || i == dim-1;//are we on the equator
					omgS += o;
					if(!eq) omgS += o;//don't double count equator
				}
			}
		}

		//@brief    : construct a back projector
		//@param geo: detector geometry to back project from
		//@param dim: side length of square legendre grid to back project onto
		//@param fct: additional scale factor for detector resizing
		//@param qu : quaternion to rotate detector location by
		template <typename Real>
		BackProjector<Real>::BackProjector(Geometry<Real> geo, const size_t dim, const Real fct, Real const * const qu) :
			pLut(std::make_shared<const Constants>(geo, dim, fct, qu)),
			rPat(std::max(geo.w * geo.h, pLut->sclr.wOut * pLut->sclr.hOut)),
			sWrk(pLut->sclr.allocateWork()),
			iVal(pLut->iPts.size()) {}

		//@brief    : unproject a pattern from the detector to a square grid
		//@param pat: pattern to back project to unit sphere
		//@param sph: location to write back projected pattern
		//@return   : sqrt(integral of pat^2) (properly weighted on the sphere), 0 if sum == false
		//@param iq : location to write pattern image quality (NULL to skip computation)
		template <typename Real>
		Real BackProjector<Real>::unproject(Real * const pat, Real * const sph, Real * const iq) {
			//rescale pattern
			Real vIq = pLut->sclr.scale(pat, rPat.data(), sWrk.data(), true, 0, NULL != iq);//rescale
			if(NULL != iq) *iq = vIq;//save iq value if needed

			//interpolate pattern intensities and compute weighted mean
			image::BiPix<Real> const * const pIpt = pLut->iPts.data();
			for(size_t i = 0; i < pLut->iPts.size(); i++) iVal[i] = pIpt[i].interpolate(rPat.data());//determine interpolated values
			const Real mean = std::inner_product(iVal.cbegin(), iVal.cend(), pLut->omeg.cbegin(), Real(0)) / pLut->omgW;//compute weighted average

			//make mean zero and compute weighted standard deviation
			Real stdev = 0;
			for(size_t i = 0; i < iVal.size(); i++) {
				iVal[i] -= mean;
				stdev += iVal[i] * iVal[i] * pLut->omeg[i];
			}
			stdev = std::sqrt(stdev / pLut->omgW);

			if(Real(0) == stdev) {//the back projected image had no contrast
				//copy to output grid making value uniformly 1
				for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = Real(1);
				return 0;
			} else {
				//make standard deviation 1
				for(size_t i = 0; i < iVal.size(); i++) iVal[i] /= stdev;

				//copy to output grid making mean 0
				for(size_t i = 0; i < pLut->iPts.size(); i++) sph[pIpt[i].idx] = iVal[i];

				//compute sum if needed
				static const Real var = std::sqrt(pLut->omgW / pLut->omgS * emsphinx::Constants<Real>::pi * Real(4));//standard deviation within window (since we normalized to 1)
				return var;
			}

		}

		//@brief    : build a mask of spherical pixels covered by the detector
		//@param sph: location to write mask
		template <typename Real>
		void BackProjector<Real>::mask(Real * const sph) {
			for(const image::BiPix<Real>& p: pLut->iPts) sph[p.idx] = 1;
		}
	}//ebsd

}//emsphinx

#endif//_DETECTOR_H_
