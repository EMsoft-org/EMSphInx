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

#ifndef _GAUSSIAN_H_
#define _GAUSSIAN_H_

#include <memory>
#include <vector>

namespace gaussian {
	//@brief a y-shifted gaussian function
	template <typename Real>
	struct Model {
		
		//model parameters
		Real a;//mean
		Real b;//2 * sigma^2
		Real c;//amplitude

		//@brief  : evauluate function f(x) = c * exp( -(x-a)^2 / b )
		//@prarm x: where to evaulate
		//@return : f(x)
		Real evaluate(const Real x) const;

		//@brief  : estimate the least squares fit gaussian model for a set of points
		//@param x   : x coordinates of points (or NULL to assume that x[i] == i)
		//@param y: y coordinates of points
		//@param n: number of points
		//@return : mean value of y
		Real estimate(Real const * const x, Real const * const y, size_t n);

		//@brief     : compute the least squares fit gaussian model for a set of points using current parameters as guess
		//@param x   : x coordinates of points (or NULL to assume that x[i] == i)
		//@param y   : y coordinates of points
		//@param n   : number of points
		//@param init: true to estimate initial values, false to use current values to initialize
		//@note      : minimizes sum( ( y - f(x) ) ^2 )
		//@return    : coefficient of determination (R^2)
		Real fit(Real const * const x, Real const * const y, size_t n, bool init = true);

	};

	//@brief: 2D Gaussian background subtraction
	template <typename Real>
	struct BckgSub2D {
		size_t m_w;//image width in pixels
		size_t m_h;//image height in pixels
		std::shared_ptr< std::vector<char> > msk;//image mask
		std::vector<Real> rWrk, cWrk;//work spaces for row / column
		gaussian::Model<Real> gx, gy;//x/y gaussian fits
		Real c;//combined x/y gaussian c

		//@brief  : construct a background subtracter with an arbitrary mask
		//@param w: image width in pixels
		//@param h: image width in pixels
		//@param m: pixel mask [assumed to be w * h with 0 for bad pixels 1 for good pixels (omitted to include all pixels)
		BckgSub2D(const size_t w, const size_t h, std::shared_ptr< std::vector<char> > m) : m_w(w), m_h(h), msk(m), rWrk(h), cWrk(w) {}
		BckgSub2D(const size_t w = 0, const size_t h = 0) : BckgSub2D(w, h, std::make_shared< std::vector<char> >(w * h, 1)) {}

		//@brief  : construct a circularly masked background subtracter
		//@param w: image width
		//@param h: image height
		//@param r: circle radius
		static BckgSub2D CircMask(const int w, const int h, const int r);
		static BckgSub2D CircMask(const int w, const int h) {return CircMask(w, h, std::min(w, h) / 2);}

		//@brief    : fit background from an image
		//@param im : image to subtract background from
		template <typename TPix>
		void fit(TPix * const im);

		//@brief   : in place subtract the 2d background from an image
		//@param im: image to subtract background from (in place)
		void subtract(uint8_t * const im);	
		void subtract(uint16_t* const im);	
		void subtract(float   * const im);	

		//@brief   : out of place subtract the 2d background from an image
		//@param im: image to subtract background from
		//@param bf: location to write background subtracted image
		template <typename TPix>
		void subtract(TPix const * const im, Real * const bf);
	};
}

#include <limits>
#include <cmath>
#include <algorithm>
#include <vector>

#include "linalg.hpp"

namespace gaussian {
	//@brief  : evauluate function f(x) = c * exp( -(x-a)^2 / b )
	//@prarm x: where to evaulate
	//@return : f(x)
	template <typename Real> Real Model<Real>::evaluate(const Real x) const {
		static_assert(std::is_floating_point<Real>::value, "gaussian model must be templated on floating point type");

		const Real dx = (x - a) ;
		return c * std::exp( -(dx * dx / b) );
	}

	//@brief  : estimate the least squares fit gaussian model for a set of points
	//@param x   : x coordinates of points (or NULL to assume that x[i] == i)
	//@param y: y coordinates of points
	//@param n: number of points
	//@return : mean value of y
	template <typename Real> Real Model<Real>::estimate(Real const * const x, Real const * const y, size_t n) {
		static_assert(std::is_floating_point<Real>::value, "gaussian model must be templated on floating point type");

		//estimate mean and scaling from max
		Real const * const yMax = std::max_element(y, y + n);
		const size_t iMax = std::distance(y, yMax);
		c = *yMax;
		a = NULL == x ? Real(iMax) : x[iMax];

		//estimate b with linear regression of ln( y[i] / c ) == -(x-a)^2 / b
		Real xy = 0, y2 = 0, yBar = 0;
		for(size_t i = 0; i < n; i++) {
			const Real yc = y[i] / c;
			yBar += y[i];
			if(yc > 0) {//don't compute log(-#)
				const Real dx  = a - (NULL == x ? Real(i) : x[i]);
				const Real yi = std::log( yc );
				y2 += yi * yi;
				xy -= yi * dx * dx;
			}
		}
		b = xy / y2;
		return yBar / n;
	}

	//@brief     : compute the least squares fit gaussian model for a set of points using current parameters as guess
	//@param x   : x coordinates of points (or NULL to assume that x[i] == i)
	//@param y   : y coordinates of points
	//@param n   : number of points
	//@param init: true to estimate initial values, false to use current values to initialize
	//@note      : minimizes sum( ( y - f(x) ) ^2 )
	//@return    : coefficient of determination (R^2)
	template <typename Real> Real Model<Real>::fit(Real const * const x, Real const * const y, size_t n, const bool init) {
		if(n < 3) throw std::runtime_error("not enough points for 3 degrees of freedom");
		Real yBar = 0;
		if(init) {
			yBar = estimate(x, y, n);
		} else {
			yBar = std::accumulate(y, y + n, Real(0)) / n;
		}

		//compute the total sum of squares once
		Real ssTot = 0;
		for(size_t i = 0; i < n; i++) {
			const Real dy = y[i] - yBar;
			ssTot += dy * dy;
		}

		//do gauss newton iteration
		const Real thr = (Real)0.0001;//convergence criterion (0.0001 => 0.01% relative step)
		const size_t maxIter = 50;
		Real ssPrev = 0, metricPrev = std::numeric_limits<Real>::max();
		for(size_t iter = 0; iter < maxIter; iter++) {
			//compute jacobian^T * jacobian and jacobian^T * residuals
			Real ss = 0;
			Real jTj[9] = {Real(0)}, jTr[3] = {Real(0)}, step[3];
			for(size_t i = 0; i < n; i++) {
				//compute derivatives of function w.r.t. x
				const Real dx  = a - (NULL == x ? Real(i) : x[i]);
				const Real dxb  = dx / b;
				const Real dfdc = std::exp( -dx * dxb );//df(x)/dc @ xi
				const Real fx   = c  * dfdc;//f(xi)
				const Real fxb  = fx * dxb ;
				const Real dfda = -fxb * Real(2);//df(x)/da @ xi
				const Real dfdb =  fxb * dxb;//df(x)/db @ xi
				const Real ri = y[i] - fx;//residual

				//accumulate
				jTj[0] += dfda * dfda; jTj[1] += dfda * dfdb; jTj[2] += dfda * dfdc;     jTr[0] += ri * dfda;
				                       jTj[4] += dfdb * dfdb; jTj[5] += dfdb * dfdc;     jTr[1] += ri * dfdb;
				                                              jTj[8] += dfdc * dfdc;     jTr[2] += ri * dfdc;
				ss += ri * ri;
			}

			//make symmetric and compute step
			jTj[3] = jTj[1];
			jTj[6] = jTj[2]; jTj[7] = jTj[5];
			solve::cholesky(jTj, step, jTr, 3);

			//apply step
			a += step[0];
			b += step[1];
			c += step[2];

			//check for convergence
			const Real metric = std::fabs((ssPrev - ss) / ss);
			if(metric >= metricPrev && metric < thr) return Real(1) - ss / ssTot;
			metricPrev = metric;
			ssPrev = ss;
		}
		throw std::runtime_error("failed to converged within allowed iterations");
	}

	//@brief  : construct a circularly masked background subtracter
	//@param w: image width
	//@param h: image height
	//@param r: circle radius
	template <typename Real>
	BckgSub2D<Real> BckgSub2D<Real>::CircMask(const int w, const int h, const int r) {
		const int w2 = w/2;
		const int h2 = h/2;
		const int rr = r*r;
		std::shared_ptr< std::vector<char> > mask = std::make_shared< std::vector<char> >(w * h, 0);
		for(int j = 0; j < h; j++) {
			const int y = j - h2;
			for(int i = 0; i < w; i++) {
				const int x = i - w2;
				if(x * x + y * y <= rr) mask->operator[](j * w + i) = 1;
			}
		}
		return BckgSub2D(w, h, mask);
	}

	//@brief    : fit background from an image
	//@param im : image to subtract background from
	template <typename Real>
	template <typename TPix>
	void BckgSub2D<Real>::fit(TPix * const im) {
		//start by computing row / column max
		std::fill(rWrk.begin(), rWrk.end(), im[0]);
		std::fill(cWrk.begin(), cWrk.end(), im[0]);
		for(size_t j = 0; j < m_h; j++) {
			for(size_t i = 0; i < m_w; i++) {
				if(1 == msk->operator[](j*m_w+i)) {
					const Real v = im[j*m_w+i];
					if(v > rWrk[j]) rWrk[j] = v;
					if(v > cWrk[i]) cWrk[i] = v;
				}
			}
		}

		//fit a guassian to row and column max
		try {
			gx.fit(NULL, cWrk.data()+1, m_w-2);
		} catch (...) {
			gx.a = Real(cWrk.size()) / 2;//mean in middle
			gx.b = std::numeric_limits<Real>::infinity();//infinite stddev (flat)
			gx.c = std::accumulate(cWrk.cbegin(), cWrk.cend(), Real(0)) / cWrk.size();//amplitude from average
		}
		try {
			gy.fit(NULL, rWrk.data()+1, m_h-2);
		} catch (...) {
			gy.a = Real(rWrk.size()) / 2;//mean in middle
			gy.b = std::numeric_limits<Real>::infinity();//infinite stddev (flat)
			gy.c = std::accumulate(rWrk.cbegin(), rWrk.cend(), Real(0)) / rWrk.size();//amplitude from average
		}

		//now re-estimate scaling using input image
		/*
		Real tk2 = 0, tky = 0;
		for(int j = 0; j < h; j++) {
			const Real dy = gy.a - j;
			const Real fy = std::exp( - (dy * dy ) / gy.b );
			for(int i = 0; i < w; i++) {
				if(1 == circMask[j * w + i]) {
					const Real dx = gx.a - i;
					const Real fx = std::exp( - (dx * dx ) / gx.b );
					const Real tk = fx * fy;
					tk2 += tk * tk  ;
					tky += tk * ptr[j*w+i];
				}
			}
		}
		c = tky / tk2;
		*/
		c = std::max(gx.c, gy.c);//this works essentially just as well but is much faster

		//compute x exponential once per column
		for(size_t i = 0; i < m_w; i++) {
			const Real dx = gx.a - i;
			cWrk[i] = std::exp( -(dx * dx / gx.b) ) * c;
		}
		for(size_t j = 0; j < m_h; j++) {
			const Real dy = gy.a - j;
			rWrk[j] = std::exp( -(dy * dy / gy.b) );
		}
	}

	//@brief    : subtract the 2d background from an image
	//@param im : image to subtract background from
	template <typename Real>
	void BckgSub2D<Real>::subtract(float * const im) {
		if(std::is_floating_point<Real>::value) {
			for(size_t j = 0; j < m_h; j++) {
				for(size_t i = 0; i < m_w; i++) {
					const size_t idx = j*m_w+i;
					im[idx] = (1 == msk->operator[](idx)) ? im[idx] - float(rWrk[j] * cWrk[i]) : 0.0;
				}
			}
		}
	}

	//@brief    : subtract the 2d background from an image
	//@param im : image to subtract background from
	template <typename Real>
	void BckgSub2D<Real>::subtract(uint8_t * const im) {
		const Real offset = c / 2;
		const uint8_t nVal = (uint8_t) (offset + Real(0.5));
		for(size_t j = 0; j < m_h; j++) {
			for(size_t i = 0; i < m_w; i++) {
				const size_t idx = j*m_w+i;
				if(1 == msk->operator[](idx)) {
					const int vNew = int(im[idx]) - (int)std::round(rWrk[j] * cWrk[i] - offset);//compute background
					if(vNew < 0) im[idx] = 0;
					else if (vNew > 0xFF) im[idx] = 0xFF;
					else im[idx] = (uint8_t) vNew;
				} else {
					im[idx] = nVal;
				}
			}
		}
	}

	//@brief    : subtract the 2d background from an image
	//@param im : image to subtract background from
	template <typename Real>
	void BckgSub2D<Real>::subtract(uint16_t * const im) {
		const Real offset = c / 2;
		const uint16_t nVal = (uint16_t) (offset + Real(0.5));
		for(size_t j = 0; j < m_h; j++) {
			for(size_t i = 0; i < m_w; i++) {
				const size_t idx = j*m_w+i;
				if(1 == msk->operator[](idx)) {
					const int vNew = int(im[idx]) - (int)std::round(rWrk[j] * cWrk[i] - offset);//compute background
					if(vNew < 0) im[idx] = 0;
					else if (vNew > 0xFFFF) im[idx] = 0xFFFF;
					else im[idx] = (uint16_t) vNew;
				} else {
					im[idx] = nVal;
				}
			}
		}
	}

	//@brief   : out of place subtract the 2d background from an image
	//@param im: image to subtract background from
	//@param bf: location to write background subtracted image
	template <typename Real>
	template <typename TPix>
	void BckgSub2D<Real>::subtract(TPix const * const im, Real * const bf) {
		for(size_t j = 0; j < m_h; j++) {
			for(size_t i = 0; i < m_w; i++) {
				const size_t idx = j*m_w+i;
				bf[idx] = (1 == msk->operator[](idx)) ? Real(im[idx]) - rWrk[j] * cWrk[i] : Real(0);
			}
		}
	}
}

#endif//_GAUSSIAN_H_
