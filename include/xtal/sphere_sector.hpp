/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                     *
 * Copyright (c) 2019, William C. Lenthe                               *
 * All rights reserved.                                                *
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
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef _sphere_sector_h_
#define _sphere_sector_h_

#include <string>
#include <vector>
#include <functional>

#include <iostream>

namespace xtal {
	namespace fs {
		//@brief   : reduce a direction to the fundamental sector
		//@param n : unit direction to reduce (magnitude is assumed to be 1) [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real> inline bool b1   (Real * const n);//-1
		template <typename Real> inline bool _211 (Real * const n);//211
		template <typename Real> inline bool _211r(Real * const n);//211 rotated 45 @ 2 (2 fold @ xy)
		template <typename Real> inline bool _121 (Real * const n);//121
		template <typename Real> inline bool _112 (Real * const n);//112
		template <typename Real> inline bool _m11 (Real * const n);//m11
		template <typename Real> inline bool _1m1 (Real * const n);//1m1
		template <typename Real> inline bool _11m (Real * const n);//11m
		template <typename Real> inline bool _12m1(Real * const n);//12/m1
		template <typename Real> inline bool _112m(Real * const n);//112/m
		template <typename Real> inline bool _222 (Real * const n);//222
		template <typename Real> inline bool _222r(Real * const n);//222r
		template <typename Real> inline bool  mm2 (Real * const n);//mm2
		template <typename Real> inline bool  mm2r(Real * const n);//mm2r
		template <typename Real> inline bool  mmm (Real * const n);//mmm
		template <typename Real> inline bool  mmmr(Real * const n);//mmmr
		template <typename Real> inline bool _4   (Real * const n);//4
		template <typename Real> inline bool b4   (Real * const n);//-4
		template <typename Real> inline bool _4m  (Real * const n);//4/m
		template <typename Real> inline bool _422 (Real * const n);//422
		template <typename Real> inline bool _4mm (Real * const n);//4mm
		template <typename Real> inline bool b42m (Real * const n);//-42m
		template <typename Real> inline bool b4m2 (Real * const n);//-4m2
		template <typename Real> inline bool _4mmm(Real * const n);//4/mmm
		template <typename Real> inline bool _3   (Real * const n);//3
		template <typename Real> inline bool _3r  (Real * const n);//3 rotated 30 degrees (use for other sectors)
		template <typename Real> inline bool b3   (Real * const n);//-3
		template <typename Real> inline bool _321 (Real * const n);//321
		template <typename Real> inline bool _312 (Real * const n);//312
		template <typename Real> inline bool _3m1 (Real * const n);//3m1
		template <typename Real> inline bool _31m (Real * const n);//31m
		template <typename Real> inline bool b3m1 (Real * const n);//-3m1
		template <typename Real> inline bool b31m (Real * const n);//-31m
		template <typename Real> inline bool _6   (Real * const n);//6
		template <typename Real> inline bool b6   (Real * const n);//-6
		template <typename Real> inline bool _6m  (Real * const n);//6/m
		template <typename Real> inline bool _622 (Real * const n);//622
		template <typename Real> inline bool _6mm (Real * const n);//6mm
		template <typename Real> inline bool b6m2 (Real * const n);//-6m2
		template <typename Real> inline bool b62m (Real * const n);//-62m
		template <typename Real> inline bool _6mmm(Real * const n);//6/mmm
		template <typename Real> inline bool _23  (Real * const n);//23
		template <typename Real> inline bool  mb3 (Real * const n);//m3
		template <typename Real> inline bool _432 (Real * const n);//432
		template <typename Real> inline bool b43m (Real * const n);//-43m
		template <typename Real> inline bool  mb3m(Real * const n);//m3m
	}

	//@brief: phenomenologically adjusted hsl2rgb
	//@param hsl: hue, saturation, lightness [0,1] (saturation is assumed to be 1)
	//@param rgb: location to write rgb [0,1]
	//@reference: Nolze, G., & Hielscher, R. (2016). Orientations–perfectly colored. Journal of Applied Crystallography, 49(5), 1786-1802.
	template <typename Real> void sph2rgb(Real const * const hsl, Real * const rgb);

	//spherical triangle -> color mappings
	// -Nolze, Gert and Hielscher Ralf. "Orientations Perfectly Colors." J. Appl. Crystallogr. 49.5 (2016): 1786-1802.
	namespace detail {
		//@brief: helper class for mapping a polygon on the sphere to the unit hemisphere
		template <size_t N, typename Real>
		struct SphericalPatch {
			//@brief      : construct a spherical triangle patch to map to the unit hemisphere
			//@param verts: vertices of spherical patch (maps to evenly spaced points to equator)
			//@param ctr  : center of spherical patch (maps to north pole)
			//@note       : verts in CCW order
			SphericalPatch(const Real verts[N][3], const Real ctr[3]) {build(verts, ctr);}

			//@brief    : convert a unit direction inside the patch to fractional polar coordinates on the hemisphere
			//@param n  : unit direction to color
			//@param tht: fractional azimuthal angle [0,1] maps to [0,2*pi]
			//@param phi: fractional polar angle [0,1] maps to [0,pi/2]
			void toHemi(Real const * const n, Real& tht, Real& phi) const;

			//@brief    : compute coloring for a unit direction in the fundamental sector
			//@param n  : unit direction in fundamental sector (undefined behavior for directions outside of sector)
			//@param rgb: location to write color [0,1]
			//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
			//@param wCn: white/black center (north/south hemisphere)
			//@param nTh: should theta be negated (to create enatiomorph / reverse color progression)
			void toColor(Real const * const n, Real * const rgb, std::function<void(Real const*const, Real*const)> h2r, const bool wCn = true, const bool nTh = false) const;

			protected:
				//@brief      : construct a spherical triangle patch to map to the unit hemisphere
				//@param verts: vertices of spherical patch (maps to evenly spaced points to equator)
				//@param ctr  : center of spherical patch (maps to north pole)
				//@param filF : fillet fraction, should be 0 for original paper, ~0.05 for perceptually uniform, must be in [0,0.5)
				//@note       : verts in CCW order
				void build(const Real verts[N][3], const Real ctr[3], const Real filF = Real(0.01));

				SphericalPatch() {}

			private:
				Real              rx          [N];//direction of red from center
				Real              ry          [N];//direction perpendicular to rx and center
				Real              center      [N];//centroid of spherical triangle
				Real              normals  [N][3];//normals of spherical triangle edges
				Real              cutoffs  [N *3];//cutoff angles (start/end of fillets and angles of vertices)
				Real              coeffs   [N][4];//fillet polynomial coefficients
				Real              cumAngles[ N+1];//cumulative angles of triangle vertices
				std::vector<Real> omega          ;//lookup table for nonlinear hue adjustment
		};

		//@brief: specialization for spherical triangles
		template < typename Real>
		struct SphericalTriangle : public SphericalPatch<3, Real> {
			//@brief       : construct a spherical triangle patch for IPF coloring
			//@param nRed  : unit direction to color red
			//@param nGreen: unit direction to color green
			//@param nBlue : unit direction to color blue
			//@note        : defines triangular patch in CCW order
			SphericalTriangle(Real const*const nRed, Real const*const nGreen, Real const*const nBlue);
		};

		//@brief: specialization for spherical wedges
		template < typename Real>
		struct SphericalWedge : public SphericalPatch<4, Real> {
			//@brief       : construct a spherical triangle patch for IPF coloring
			//@param nGreen: 2D unit direction to color green
			//@param nBlue : 2D unit direction to color blue
			//@note        : red is assumed to be at [0,0,1]
			SphericalWedge(Real const*const nGreen, Real const*const nBlue);
		};
	}
}

////////////////////////////////////////////////////////////////////////
//                       Implementation Details                       //
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>

#include "constants.hpp"

namespace xtal {

	////////////////////////////////////////////////////////////////////////
	//                   Fundamental Sector Reductions                    //
	////////////////////////////////////////////////////////////////////////

	namespace fs {
		//@brief   : rotate a 2d vector in place
		//@param xy: vector to rotate
		//@param c : cosine of rotation angle
		//@param s : sine of rotation angle
		template <typename Real> inline void rot2d(Real * const xy, const Real c, const Real s) {
			const Real xp = c * xy[0] - s * xy[1];
			const Real yp = c * xy[1] + s * xy[0];
			xy[0] = xp;
			xy[1] = yp;
		}

		//@brief   : mirror a 2d vector in place
		//@param xy: vector to mirror
		//@param c : cosine of mirror plane angle
		//@param s : sine of mirror plane angle
		template <typename Real> inline void mir2d(Real * const xy, const Real c, const Real s) {
			//compute nearest point to xy on mirror plane
			const Real d = c * xy[0] + s * xy[1];
			const Real pt[2] = {d * c, d * s};

			//reflect
			xy[0] = pt[0] - (xy[0] - pt[0]);
			xy[1] = pt[1] - (xy[1] - pt[1]);
		}

		//@brief  : move a point in the first octant to the maximum z rotation possible by 120 @ 111
		//@param n: vector to rotate [in place]
		//@return : true if z was already maximize, false otherwise
		template <typename Real> inline bool r111(Real * const n) {
			//get max and return if we're done
			const size_t idx = std::distance(n, std::max_element(n, n+3));
			if(n[2] == n[idx]) return true;
			std::rotate(n, n+1+idx, n+3);
			return false;
		}

		//@brief   : reduce a direction to the -1 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool b1   (Real * const n) {
			const bool ret = std::signbit(n[2]);
			if(ret) {
				n[0] = -n[0];
				n[1] = -n[1];
				n[2] = -n[2];
			}
			return !ret;
		}

		//@brief   : reduce a direction to the 211 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _211 (Real * const n) {
			const bool ret = std::signbit(n[2]);
			if(ret) {
				n[1] = -n[1];
				n[2] = -n[2];
			}
			return !ret;
		}

		//@brief   : reduce a direction to the 211 fundamental sector rotated 45 @ z
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _211r(Real * const n) {
			const bool ret = std::signbit(n[2]);
			if(ret) {
				std::swap(n[0], n[1]);
				n[2] = -n[2];
			}
			return !ret;
		}

		//@brief   : reduce a direction to the 121 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _121 (Real * const n) {
			const bool ret = std::signbit(n[2]);
			if(ret) {
				n[0] = -n[0];
				n[2] = -n[2];
			}
			return !ret;
		}

		//@brief   : reduce a direction to the 112 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _112 (Real * const n) {
			const bool ret = std::signbit(n[1]);
			if(ret) {
				n[0] = -n[0];
				n[1] = -n[1];
			}
			return !ret;
		}

		//@brief   : reduce a direction to the m11 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _m11 (Real * const n) {
			const bool ret = std::signbit(n[0]);
			if(ret) n[0] = -n[0];
			return !ret;
		}

		//@brief   : reduce a direction to the 1m1 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _1m1 (Real * const n) {
			const bool ret = std::signbit(n[1]);
			if(ret) n[1] = -n[1];
			return !ret;
		}

		//@brief   : reduce a direction to the 11m fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _11m (Real * const n) {
			const bool ret = std::signbit(n[2]);
			if(ret) n[2] = -n[2];
			return !ret;
		}

		//@brief   : reduce a direction to the 222r fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _222r(Real * const n) {
			bool ret = _211r(n);//first bring to +z hemisphere with rotation about xy
			if(n[1] < n[0]) {//are we below the y==x line
				_112(n);
				ret = false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the mm2r fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool mm2r(Real * const n) {
			const bool ret = _112(n);//first bring to +y hemisphere with rotation about z
			const Real ax = std::fabs(n[0]);
			if(ax > n[1]) {
				n[0] = std::copysign(n[1], n[0]);
				n[1] = ax;
				return false;
			}
			return ret;
		}
		
		//@brief   : reduce a direction to the 4 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _4   (Real * const n) {
			const bool ret = _112(n);//first move to +y with 2 fold
			if(std::signbit(n[0])) {//need to move to +x with 4 fold
				n[0] = -n[0];
				std::swap(n[0], n[1]);
				return false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the -4 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool b4   (Real * const n) {
			bool ret = true;
			if(std::signbit(n[2])) {//need to move to +z with -4
				n[0] = -n[0];
				std::swap(n[0], n[1]);
				n[2] = -n[2];
				ret = false;
			}
			if(std::signbit(n[1])) {//need to move to +y with 2
				n[0] = -n[0];
				n[1] = -n[1];
				ret = false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the 4mm fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _4mm (Real * const n) {
			const bool ret = _4(n);//first move to first quadrant with 4 fold
			if(n[1] > n[0]) {//use xy mirror if needed
				std::swap(n[0], n[1]);
				return false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the 3 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _3  (Real * const n) {
			//determine sector
			const Real t = std::fabs(n[1] / n[0]);//|tan(theta)|
			size_t idx = std::signbit(n[1]) ? 2 : 0;//0, 1, or 2 for first, second, or third sector, initialize with easy cases
			if(std::signbit(n[0]) && t <= Constants<Real>::r3) idx = 1;//check for second sector

			//apply rotation
			switch(idx) {
				case 0:                                       return true ;//in first sector, we're done
				case 1: rot2d(n, Real(-0.50), -Constants<Real>::r3_4); return false;//in second sector, rotate -120 @ z
				case 2: rot2d(n, Real(-0.50),  Constants<Real>::r3_4); return false;//in third sector, rotate 120 @ z
			}
			throw std::logic_error("unhandled 3 fold case");
			return false;
		}

		//@brief   : reduce a direction to the 3 fundamental sector rotated 30 degrees about z
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _3r (Real * const n) {
			//determine sector
			const Real t = std::fabs(n[1] / std::max(std::fabs(n[0]), std::numeric_limits<Real>::epsilon()));//|tan(theta)| if in +y hemisphere
			size_t idx = (std::signbit(n[0]) ? 1 : 2);//0, 1, or 2 for first, second, or third sector, initialize with easy cases
			if(!std::signbit(n[1]) && t >= Real(1) / Constants<Real>::r3) idx = 0;//check for first sector

			//apply rotation
			switch(idx) {
				case 0:                                       return true ;//in first sector, we're done
				case 1: rot2d(n, Real(-0.50), -Constants<Real>::r3_4); return false;//in second sector, rotate -120 @ z
				case 2: rot2d(n, Real(-0.50),  Constants<Real>::r3_4); return false;//in third sector, rotate 120 @ z
			}
			throw std::logic_error("unhandled 3 fold case");
			return false;
		}

		//@brief   : reduce a direction to the 31m fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _31m (Real * const n) {
			const bool ret = _3(n);//first reduce to 3 fold fz
			const Real t = n[1] / n[0];
			if(n[1] < std::numeric_limits<Real>::epsilon() || !(Real(0) <= t && t <= Constants<Real>::r3)) {//check if we're above 60 degrees (< eps instead of signbit to handle divide by 0)
				mir2d(n, Real(0.5), Constants<Real>::r3_4);
				return false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the 6 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _6   (Real * const n) {
			const bool ret = _112(n);//first bring to +y hemisphere

			//determine sector
			const Real t = std::fabs(n[1] / std::max(std::fabs(n[0]), std::numeric_limits<Real>::epsilon()));//|tan(theta)| if in +y hemisphere
			size_t idx = std::signbit(n[0]) ? 2 : 0;
			if(t > Constants<Real>::r3) idx = 1;//check for second sector

			//apply rotation
			switch(idx) {
				case 0:                                       return ret  ;//in first sector, we're done
				case 1: rot2d(n, Real( 0.50), -Constants<Real>::r3_4); return false;//in second sector, rotate -60 @ z
				case 2: rot2d(n, Real(-0.50), -Constants<Real>::r3_4); return false;//in third sector, rotate -120 @ z
			}
			throw std::logic_error("unhandled 6 fold case");
			return false;
		}

		//@brief   : reduce a direction to the 6mm fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _6mm (Real * const n) {
			const bool ret = _6(n);//start by reducing to 6 fold fs

			//now check if we're above 30 degrees
			const Real t = n[1] / n[0];//tan(theta)
			if(t > Real(1) / Constants<Real>::r3) {//above 30 degres
				mir2d(n, Constants<Real>::r3_4, Real(0.5));
				return false;
			}
			return ret;
		}

		//@brief   : reduce a direction to the 23 fundamental sector
		//@param n : unit direction to reduce [in place]
		//@return  : true/false if n was/wasn't already in the fundamental sector
		template <typename Real>
		inline bool _23 (Real * const n) {
			bool ret = _211(n);//first bring to +z hemisphere with rotation about y
			if(n[1] < n[0]) {//are we below the y==x line
				_112(n);
				ret = false;
			}

			//now handle 3 fold @ 111
			const Real v[5] = {n[0], n[1], -n[0], -n[0], n[2]};
			const size_t idx = std::distance(v, std::max_element(v, v+5));//find maximum potential z
			if(n[2] == v[idx]) return ret;

			//negate x/y if needed
			if(idx > 1) {
				n[0] = -n[0];
				n[1] = -n[1];
			}

			//reshuffle (rotation about 111) and bring back above the y==x line if needed
			std::rotate(n, n+1+(idx%2), n+3);
			if(n[1] < n[0]) _112(n);//are we below the y==x line
			return false;
		}

		template <typename Real> inline bool _12m1(Real * const n) {const bool tmp = _121 (n); return _1m1(n) && tmp;}//12/m1
		template <typename Real> inline bool _112m(Real * const n) {const bool tmp = _112 (n); return _11m(n) && tmp;}//112/m
		template <typename Real> inline bool _222 (Real * const n) {const bool tmp = _121 (n); return _112(n) && tmp;}//222
		template <typename Real> inline bool  mm2 (Real * const n) {const bool tmp = _m11 (n); return _1m1(n) && tmp;}//mm2
		template <typename Real> inline bool  mmm (Real * const n) {const bool tmp = mm2  (n); return _11m(n) && tmp;}//mmm
		template <typename Real> inline bool  mmmr(Real * const n) {const bool tmp = mm2r (n); return _11m(n) && tmp;}//mmmr
		template <typename Real> inline bool _4m  (Real * const n) {const bool tmp = _4   (n); return _11m(n) && tmp;}//4/m
		template <typename Real> inline bool _422 (Real * const n) {const bool tmp = _121 (n); return _4  (n) && tmp;}//422
		template <typename Real> inline bool b42m (Real * const n) {const bool tmp = _121 (n); return mm2r(n) && tmp;}//-42m
		template <typename Real> inline bool b4m2 (Real * const n) {const bool tmp = _211r(n); return mm2 (n) && tmp;}//-4m2
		template <typename Real> inline bool _4mmm(Real * const n) {const bool tmp = _11m (n); return _4mm(n) && tmp;}//4/mmm
		template <typename Real> inline bool b3   (Real * const n) {const bool tmp = b1   (n); return _3  (n) && tmp;}//-3
		template <typename Real> inline bool _321 (Real * const n) {const bool tmp = _211 (n); return _3  (n) && tmp;}//321
		template <typename Real> inline bool _312 (Real * const n) {const bool tmp = _121 (n); return _3r (n) && tmp;}//312
		template <typename Real> inline bool _3m1 (Real * const n) {const bool tmp = _3r  (n); return _m11(n) && tmp;}//3m1
		template <typename Real> inline bool b3m1 (Real * const n) {const bool tmp = _211 (n); return _3m1(n) && tmp;}//-3m1
		template <typename Real> inline bool b31m (Real * const n) {const bool tmp = _121 (n); return _31m(n) && tmp;}//-31m
		template <typename Real> inline bool b6   (Real * const n) {const bool tmp = _11m (n); return _3  (n) && tmp;}//-6
		template <typename Real> inline bool _6m  (Real * const n) {const bool tmp = _11m (n); return _6  (n) && tmp;}//6/m
		template <typename Real> inline bool _622 (Real * const n) {const bool tmp = _121 (n); return _6  (n) && tmp;}//622
		template <typename Real> inline bool b6m2 (Real * const n) {const bool tmp = _121 (n); return _3m1(n) && tmp;}//-6m2
		template <typename Real> inline bool b62m (Real * const n) {const bool tmp = _211 (n); return _31m(n) && tmp;}//-62m
		template <typename Real> inline bool _6mmm(Real * const n) {const bool tmp = _11m (n); return _6mm(n) && tmp;}//6/mmm
		template <typename Real> inline bool  mb3 (Real * const n) {const bool tmp = mmm  (n); return r111(n) && tmp;}//m3
		template <typename Real> inline bool _432 (Real * const n) {const bool tmp = _422 (n); return r111(n) && tmp;}//432
		template <typename Real> inline bool b43m (Real * const n) {const bool tmp = _23  (n); return mm2r(n) && tmp;}//-43m
		template <typename Real> inline bool  mb3m(Real * const n) {const bool tmp = mb3  (n); return _4mm(n) && tmp;}//m3m
	}

	//@brief: phenomenologically adjusted hsl2rgb
	//@param hsl: hue, saturation, lightness [0,1] (saturation is assumed to be 1)
	//@param rgb: location to write rgb [0,1]
	//@reference: Nolze, G., & Hielscher, R. (2016). Orientations–perfectly colored. Journal of Applied Crystallography, 49(5), 1786-1802.
	template <typename Real>
	void sph2rgb(Real const * const hsl, Real * const rgb) {
		//get lightness and saturation rescaling parameters
		const bool whiteCenter = hsl[2] >= Real(0.5);
		const Real yL = whiteCenter ? Real(0.25) : Real(0.50);
		const Real yS = whiteCenter ? Real(0.20) : Real(0.50);

		//constants for nonlinear hue adjustment
		static const Real iDen = Real( 0.570990316610288181236261564684297686279447800757942106831501845990856895);//1 / ( 1 + sqrt(2*pi) * std::erf( 5 * sqrt(2) / 3 ) * 0.3 );
		static const Real k1   = Real( 0.125331413731550025120788264240552262650349337030496915831496178817114683);// sqrt(pi/2) / 10
		static const Real k2   = Real(14.1421356237309504880168872420969807856967187537694807317667973799073248  );//10 * sqrt(2)
		static const Real k1_3 = Real( 0.333333333333333333333333333333333333333333333333333333333333333333333333);// 1/3
		static const Real k1_6 = Real( 0.166666666666666666666666666666666666666666666666666666666666666666666667);// 1/6
		static const Real pi_2 = Real(1.57079632679489661923132169163975144209858469968755291048747229615390820  );// pi/2

		//adjust hue gradient (A.5)
		const Real h3   = std::fmod(hsl[0], k1_3);
		const bool half = h3 > k1_6;
		const Real h6   = half ? k1_3 - h3 : h3;
		const Real hNew = (h6 + k1 * std::erf(k2 * h6)) * iDen;
		rgb[0] = hsl[0] - h3 + (half ? k1_3 - hNew : hNew);//save adjusted hue

		//adjust lightness gradient (A.9)
		const Real sP   = std::sin(hsl[2] * pi_2);
		const Real th   = yL * hsl[2] + (Real(1) - yL) * sP * sP;
		const Real gray = Real(1) - Real(2) * yS * std::fabs(th - Real(0.5));
		rgb[2] = (th - Real(0.5)) * gray + Real(0.5);//save adjusted lightness

		//adjust saturation gradient (A.10)
		rgb[1] = gray * ( Real(1) - std::fabs( Real(2) * th - Real(1) ) ) / ( Real(1) - std::fabs( Real(2) * hsl[2] - Real(1) ) );//save adjusted saturation
		if(std::isnan(rgb[1])) rgb[1] = Real(0);//correct for divide by 0

		//convert adjusted hsl to rgb
		const Real c = (Real(1) - std::fabs(rgb[2] * 2 - 1)) * rgb[1];//compute chroma
		const Real m = rgb[2] - c/2;//m
		const Real h = rgb[0] * 6;//hue [0,1] -> [0,6]
		const Real x = c * (Real(1) - std::fabs(std::fmod(h, 2) - 1));
		switch((size_t)h) {
			case 6://intentional fall through
			case 0: rgb[0] = c+m; rgb[1] = x+m; rgb[2] =   m; return;
			case 1: rgb[0] = x+m; rgb[1] = c+m; rgb[2] =   m; return;
			case 2: rgb[0] =   m; rgb[1] = c+m; rgb[2] = x+m; return;
			case 3: rgb[0] =   m; rgb[1] = x+m; rgb[2] = c+m; return;
			case 4: rgb[0] = x+m; rgb[1] =   m; rgb[2] = c+m; return;
			case 5: rgb[0] = c+m; rgb[1] =   m; rgb[2] = x+m; return;
		}
	}


	////////////////////////////////////////////////////////////////////////
	//                         Spherical Patches                          //
	////////////////////////////////////////////////////////////////////////

	namespace detail {
		namespace detail {
			//@brief   : 3d vector dot product
			//@param v1: first vector
			//@param v2: second vector
			//@return  : dot product
			template <typename Real>
			Real dot(Real const * const v1, Real const * const v2) {
				return std::inner_product(v1, v1+3, v2, Real(0));
			}
			
			//@brief  : 3d vector normalization
			//@param v: vector to normalize (in place)
			template <typename Real>
			void normalize(Real * const v) {
				const Real mag = std::sqrt(dot(v, v));
				std::transform(v, v+3, v, [mag](const Real& i){return i/mag;});
			}
			
			//@brief   : 3d vector cross product
			//@param v1: first vector
			//@param v2: second vector
			//@param x : location to write v1 cross v2
			template <typename Real>
			void cross(Real const * const v1, Real const * const v2, Real * const x) {
				const Real cross[3] = {
					x[0] = v1[1] * v2[2] - v1[2] * v2[1],
					x[1] = v1[2] * v2[0] - v1[0] * v2[2],
					x[2] = v1[0] * v2[1] - v1[1] * v2[0],
				};
				std::copy(cross, cross+3, x);
			}
		}

		//@brief       : construct a spherical triangle patch for IPF coloring
		//@param nRed  : unit direction to color red
		//@param nGreen: unit direction to color green
		//@param nBlue : unit direction to color blue
		//@note        : defines triangular patch in CCW order
		template <typename Real>
		SphericalTriangle<Real>::SphericalTriangle(Real const*const nRed, Real const*const nGreen, Real const*const nBlue) {
			//copy normals to array for convenience
			Real verts[3][3];
			std::copy(nRed  , nRed  +3, verts[0]);
			std::copy(nGreen, nGreen+3, verts[1]);
			std::copy(nBlue , nBlue +3, verts[2]);

			//first check if the three points are outside the hemisphere
			const Real det = verts[0][0] * verts[1][1] * verts[2][2]
			               + verts[0][1] * verts[1][2] * verts[2][0]
			               + verts[0][2] * verts[1][0] * verts[2][1]
			               - verts[0][0] * verts[1][2] * verts[2][1]
			               - verts[0][1] * verts[1][0] * verts[2][2]
			               - verts[0][2] * verts[1][1] * verts[2][0];
			if(det < std::numeric_limits<Real>::epsilon()) throw std::runtime_error("spherical triangle must be within single hemisphere");

			//compute center of spherical triangle (assumes triangle covers less than hemisphere)
			Real ctr[3];
			for(size_t i = 0; i < 3; i++) ctr[i] = verts[0][i] + verts[1][i] + verts[2][i];
			detail::normalize(ctr);

			//now build a general spherical patch
			SphericalPatch<3, Real>::build(verts, ctr);
		};


		//@brief       : construct a spherical triangle patch for IPF coloring
		//@param nGreen: 2D unit direction to color green
		//@param nBlue : 2D unit direction to color blue
		//@note        : red is assumed to be at [0,0,1]
		template <typename Real>
		SphericalWedge<Real>::SphericalWedge(Real const*const nGreen, Real const*const nBlue) {
			//copy normals to array for convenience
			Real verts[4][3];
			verts[0][0] = verts[0][1] = 0;         verts[0][2] = 1;
			std::copy(nGreen, nGreen+2, verts[1]); verts[1][2] = 0;
			verts[2][0] = verts[2][1] = 0;         verts[2][2] =-1;
			std::copy(nBlue , nBlue +2, verts[3]); verts[3][2] = 0;

			//compute center
			Real ctr[3] = {
				nGreen[0] + nBlue[0],
				nGreen[1] + nBlue[1],
				0
			};

			//handle special case of antipodal blue/green
			const Real dot = nGreen[0] * nBlue[0] + nGreen[1] * nBlue[1];
			if(std::fabs(dot + Real(1)) < std::numeric_limits<Real>::epsilon()) {
				ctr[0] = -nGreen[1];
				ctr[1] =  nGreen[0];
			}
			detail::normalize(ctr);

			//now build a general spherical patch
			SphericalPatch<4, Real>::build(verts, ctr);
		}

		//@brief      : construct a spherical triangle patch to map to the unit hemisphere
		//@param verts: vertices of spherical patch (maps to evenly spaced points to equator)
		//@param ctr  : center of spherical patch (maps to north pole)
		//@param filF : fillet fraction, should be 0 for original paper, ~0.05 for perceptually uniform, must be in [0,0.5)
		//@note       : verts in CCW order
		template <size_t N, typename Real>
		void SphericalPatch<N, Real>::build(const Real verts[N][3], const Real ctr[3], const Real filF) {
			//save center
			std::copy(ctr, ctr+3, center);

			//build orthogonal coordinate system for center -> each vertex
			Real vx[N][3], vy[N][3];
			for(size_t i =0; i < N; i++) {
				detail::cross(center, verts[i], vy[i]);
				detail::cross(vy[i] , center  , vx[i]);
				detail::normalize(vx[i]);
				detail::normalize(vy[i]);
			}
			std::copy(vx[0], vx[0]+3, rx);//red is the global x direction
			std::copy(vy[0], vy[0]+3, ry);//global y is perpindicular to x and patch center (z)

			//compute angles between successive verts
			Real angles[N];
			for(size_t i = 0; i < N; i++) angles[i] = std::acos(detail::dot(vx[i], vx[(i+1)%N]));
			std::partial_sum(angles, angles+N, cumAngles+1);

			//compute normals of circles defining edges of domain
			for(size_t i = 0; i < N; i++) {
				detail::cross(verts[i], verts[(i+1)%N], normals[i]);
				detail::normalize(normals[i]);
			}

			//compute cutoff angles for filleting
			Real deltas[N];
			std::fill(deltas, deltas + N, Constants<Real>::pi2 / N);//for now just evenly space deltas
			for(size_t i = 0; i < N; i++) {//loop over edges
				for(size_t j = 0; j < 3; j++) {//loop over verts
					const Real delta = 2 == j ? 0 : std::copysign(filF * deltas[i], 0 == j ? 1 : -1);
					cutoffs[i*3+j] = cumAngles[ i +(j==0 ? 0 : 1)] + delta;
				}
			}

			//numerically compute r and dr/dtheta at transition points between linear and filleted regions
			Real radii[2*N], dRadii[2*N];
			for(size_t i = 0; i < N; i++) {//loop over edges
				//compute angles where radius needs to be calculated
				const Real dT (0.1);//angular offset from vertex -> transition points (for numerical derivative calculation), ~1/2 degree
				const Real hf(0.01);//fractional distance numerical derivative calculation (should probably be <= 1)
				const Real thetas[6] = {
					cutoffs[3*i+0] - hf * dT,//symmetric points for derivative calculation
					cutoffs[3*i+0]          ,//first transition
					cutoffs[3*i+0] + hf * dT,//symmetric points for derivative calculation
					cutoffs[3*i+1] - hf * dT,//symmetric points for derivative calculation
					cutoffs[3*i+1]          ,//second transition
					cutoffs[3*i+1] + hf * dT //symmetric points for derivative calculation
				};

				//apply rotations and compute radius/angle at each point
				Real r[6];
				for(size_t j = 0; j < 6; j++) {
					//compute normal of circle at desired angle (ry rotated about center)
					const Real c = std::cos(thetas[j] / 2);
					const Real s = std::sin(thetas[j] / 2);

					//q * n (w == 0 since rotation axis is perpendicular to vector)
					const Real x = c * ry[0] + s * (center[1] * ry[2] - center[2] * ry[1]);
					const Real y = c * ry[1] + s * (center[2] * ry[0] - center[0] * ry[2]);
					const Real z = c * ry[2] + s * (center[0] * ry[1] - center[1] * ry[0]);

					const Real m[3] = {//q * n * q.conj() [normal of circle at desired angle]
						x * c + s * (z * center[1] - y * center[2]),
						y * c + s * (x * center[2] - z * center[0]),
						z * c + s * (y * center[0] - x * center[1]),
					};

					//now compute intersection of two unit circles at origin w/ normals v and normals[edge]
					const Real& nx = normals[i][0];
					const Real& ny = normals[i][1];
					const Real& nz = normals[i][2];
					const Real& mx = m[0];
					const Real& my = m[1];
					const Real& mz = m[2];
					const Real den = std::sqrt(  nx * nx * ( my * my + mz * mz ) + ny * ny * ( mz * mz + mx * mx ) + nz * nz * ( mx * mx + my * my ) - Real(2) * ( nz * nx * mz * mx + nx * ny * mx * my + ny * nz * my * mz ) );

					Real v[3] = {//intersection of two circles (point along edge i at angle thetas[j])
						(ny * mz - nz * my) / den,
						(nz * mx - nx * mz) / den,
						(nx * my - ny * mx) / den,
					};
					if(std::signbit(detail::dot(v, center))) std::transform(v, v+3, v, std::negate<Real>());//select intersection point closest to center
					r[j] = std::acos(detail::dot(v, center));//compute angle from center -> edge at this theta
				}

				//save radii and compute derivative
				radii [i*2+0] = r[1];
				radii [i*2+1] = r[4];
				dRadii[i*2+0] = (r[2] - r[0]) / (hf * dT * 2);
				dRadii[i*2+1] = (r[5] - r[3]) / (hf * dT * 2);
			}

			//compute polynomial coefficients to remove discontinuity in r
			for(size_t i = 0; i < N; i++) {//loop over edge
				const size_t j = (i+1)%N;//get index of next edge
				const Real v1 = radii [i*2+1];//value of radius at transition point in edge i (near edge j)
				const Real v2 = radii [j*2+0];//value of radius at transition point in edge j (near edge i)
				const Real m1 = dRadii[i*2+1] * filF * angles[i];//value of d(radius)/d(theta) at transition point (multiply by range of -1->0 to correct derivative scaling)
				const Real m2 = dRadii[j*2+0] * filF * angles[j];//value of d(radius)/d(theta) at transition point (multiply by range of  0->1 to correct derivative scaling)
				coeffs[i][0] = ( m1 + m2 + v1     - v2    ) / 4;
				coeffs[i][1] = (-m1 + m2                  ) / 4;
				coeffs[i][2] = (-m1 - m2 - v1 * 3 + v2 * 3) / 4;
				coeffs[i][3] = ( m1 - m2 + v1 * 2 + v2 * 2) / 4;
			}

			//build lookup table for nonlinear hue adjustment (to keep hue at vert i i/N)
			//this is a numerical lookup table to solve A.6

			//compute the fractional angle of each vertex
			Real rhoN[N];
			for(size_t i = 1; i < N; i++) {
				Real v[3];
				std::transform(verts[i], verts[i] + 3, center, v, std::minus<Real>());
				Real angle = std::atan2(detail::dot(ry, v), detail::dot(rx, v));
				if(std::signbit(angle)) angle+= Constants<Real>::pi2;
				rhoN[i-1] = angle / Constants<Real>::pi2;//convert angle w.r.t. red --> fractional
			}
			rhoN[N-1] = 1;

			//create evenly spaced list for angle from 0->1
			omega.resize(256 * N);//~256 point between successive verts
			std::vector<Real> irho(omega.size());
			irho.resize(omega.size());
			std::iota(irho.begin(), irho.end(), Real(0));
			std::for_each(irho.begin(), irho.end(), [&](Real&i){i /= Real(irho.size() - 1);});

			//compute the distance to the sector edge at each angle (in irho)
			omega[0] = 0;
			for(size_t i = 0; i < omega.size() - 1; i++) {
				//create vector normal to center at angle irho[i]
				Real n[3];
				Real s = std::sin(Constants<Real>::pi2 * irho[i]);
				Real c = std::cos(Constants<Real>::pi2 * irho[i]);
				std::transform(rx, rx+3, ry, n, [s, c](Real i, Real j){return i * s - j * c;});

				//determine which edge is closest and compute distance to edge
				Real normxn[3];
				const size_t j = std::distance(rhoN, std::upper_bound(rhoN, rhoN+N, irho[i]));//which edge is closest
				detail::cross(normals[j], n, normxn);
				const Real mag = std::sqrt(detail::dot(normxn, normxn));
				omega[i+1] = std::acos(detail::dot(normxn, center) / mag);
			}

			//get the offset to the vertices
			size_t idx[N+1] = {0};
			for(size_t i = 0; i < N; i++) idx[i+1] = std::distance(irho.begin(), std::upper_bound(irho.begin(), irho.end(), rhoN[i]));

			//normalize
			for(size_t i = 0; i < N; i++) {
				const Real sum = std::accumulate(omega.begin() + idx[i], omega.begin() + idx[i+1], Real(0)) * N;//we want the distance from vert i -> i+1 to cover 1/N
				std::for_each                   (omega.begin() + idx[i], omega.begin() + idx[i+1], [sum](Real& r){r /= sum;});
			}

			//integrate
			std::partial_sum(omega.begin(), omega.end(), omega.begin());//now we have our lookup table
		}

		//@brief    : convert a unit direction inside the patch to fractional polar coordinates on the hemisphere
		//@param n  : unit direction to color
		//@param tht: fractional azimuthal angle [0,1] maps to [0,2*pi]
		//@param phi: fractional polar angle [0,1] maps to [0,pi/2]
		template <size_t N, typename Real>
		void SphericalPatch<N, Real>::toHemi(Real const * const n, Real& tht, Real& phi) const {
			//compute angle with red direction
			Real v[3];
			const Real n0[3] = {n[0], n[1], n[2]};//in case input and output overlap
			std::transform(n0, n0 + 3, center, v, std::minus<Real>());
			Real angle = std::atan2(detail::dot(ry, v), detail::dot(rx, v));
			if(std::signbit(angle)) angle+= Constants<Real>::pi2;
			tht = angle / Constants<Real>::pi2;//convert angle w.r.t. red --> fractional

			//apply adaptive hue gradient
			tht *= omega.size() - 1;//rescale from [0,1] to lookup table size
			const size_t iOmg = (size_t)tht;//get index of lower bound in lookup table
			if(iOmg+1 < omega.size()) {//don't go out of bounds if tht == 1
				tht = omega[iOmg] + (omega[iOmg+1] - omega[iOmg]) * ((iOmg+1) - tht);//linearly interpolate from lookup table
			} else {
				tht = Real(1);
			}

			//compute polar angle
			const size_t idx = std::distance(cutoffs, std::lower_bound(cutoffs, cutoffs + 3 * N, angle));//determine which region this angle falls in
			const size_t i = idx / 3;//index of edge
			phi = std::acos(detail::dot(n0, center));//angle between center and point
			if(phi < std::numeric_limits<Real>::epsilon()) return;//avoid divide by zero issues

			//normalize polar angle
			switch(idx - i * 3) {
				case 1: {//in linear region, normalize angle by max possible angle
					Real nxc[3];
					detail::cross(n0        , center, nxc);//normal of arc through n/center
					detail::cross(normals[i], nxc   , v  );//intersection of two circles (edge of patch in direction tht)
					detail::normalize(v);
					phi /= std::acos(detail::dot(v, center)) * Real(2);//compute fractional progress ot edge
				} break;

				case 0: {//in first fillet
					const size_t j = (i+N-1)%N;//get i-1 with periodic boundary conditions
					Real x =  (angle - cumAngles[i  ]) / (cutoffs[idx  ] - cumAngles[i  ]);
					const Real den = coeffs[j][0] * x * x * x + coeffs[j][1] * x * x + coeffs[j][2] * x + coeffs[j][3];
					phi /= std::max(phi, den) * Real(2);//normalize, clipping at 1/2
				} break;

				case 2: {//in second fillet
					Real x = -(angle - cumAngles[i+1]) / (cutoffs[idx-1] - cumAngles[i+1]);
					const Real den = coeffs[i][0] * x * x * x + coeffs[i][1] * x * x + coeffs[i][2] * x + coeffs[i][3];
					phi /= std::max(phi, den) * Real(2);//normalize, clipping at 1/2
				} break;
			}
		}

		//@brief    : compute coloring for a unit direction in the fundamental sector
		//@param n  : unit direction in fundamental sector (undefined behavior for directions outside of sector)
		//@param rgb: location to write color [0,1]
		//@param h2r: hsl2rgb like coloring function to use with h, s, and l in [0,1] and output as [0,1] rgb: void(Real const * const hsl, Real * const rgb)
		//@param wCn: white/black center (north/south hemisphere)
		//@param nTh: should theta be negated (to create enatiomorph / reverse color progression)
		template <size_t N, typename Real>
		void SphericalPatch<N, Real>::toColor(Real const * const n, Real * const rgb, std::function<void(Real const*const, Real*const)> h2r, const bool wCn, const bool nTh) const {
			toHemi(n, rgb[0], rgb[2]);
			rgb[1] = Real(1);//fully satruated
			if(nTh) rgb[0] = Real(1) - rgb[0];
			if(wCn) rgb[2] = Real(1) - rgb[2];
			h2r(rgb, rgb);
		}
	}
}

#endif//_sphere_sector_h_
