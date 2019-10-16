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

#ifndef _ROTATION_H_
#define _ROTATION_H_

//orientation transform routines based on
// -Rowenhorst, David, et al. "Consistent Representations of and Conversions Between 3D Rotations." Model. Simul. Mater. Sci. Eng. 23.8 (2015): 083501.
// -Rosca, D., et al. "A New Method of Constructing a Grid in the Space of 3D rotations and its Applications to Texture Analysis." Model. Simul. Mater. Sci. Eng. 22.7 (2014): 075013.
// -fortran implementation of routines by Marc De Graef (https://github.com/marcdegraef/3Drotations)

//the following conventions are used:
// -passive rotations
// -quaternions as [w, x, y, z]
// -0 <= rotation angle <= pi
// -rotation axis in positive z hemisphere for rotations of pi (+y for z == 0, 100 for y == z == 0)
// -rotation axis = [0, 0, 1] for rotations of 0

#include <iostream>

namespace xtal {

	//@brief   : convert between 2 different orientation representations
	//@param xx: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param yy: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : xx and yy can be the same memory location
	//@note    : orientation representation abbreviations
	//           eu - ZXZ euler angles (Bunge convention, radians)
	//           om - orientation matrix (9 component vector in row major order)
	//           ax - axis angle (nx, ny, nz, w[radians])
	//           ro - rodrigues vector (nx, ny, nz, tan(w/2))
	//           qu - quaternion (w, x, y, z)
	//           ho - homochoric (x, y, z)
	//           cu - cubochoric (x, y, z)
	template<typename Real> void eu2om(Real const * const eu, Real * const om);//A.1
	template<typename Real> void eu2ax(Real const * const eu, Real * const ax);//A.2
	template<typename Real> void eu2ro(Real const * const eu, Real * const ro);//A.3
	template<typename Real> void eu2qu(Real const * const eu, Real * const qu);//A.4
	template<typename Real> void eu2ho(Real const * const eu, Real * const ho);//eu->ax->ho
	template<typename Real> void eu2cu(Real const * const eu, Real * const cu);//eu->ho->cu
	template<typename Real> void om2eu(Real const * const om, Real * const eu);//A.5
	template<typename Real> void om2ax(Real const * const om, Real * const ax);//A.6
	template<typename Real> void om2ro(Real const * const om, Real * const ro);//om->eu->ro
	template<typename Real> void om2qu(Real const * const om, Real * const qu);//A.7
	template<typename Real> void om2ho(Real const * const om, Real * const ho);//om->ax->ho
	template<typename Real> void om2cu(Real const * const om, Real * const cu);//om->ho->cu
	template<typename Real> void ax2eu(Real const * const ax, Real * const eu);//ax->om->eu
	template<typename Real> void ax2om(Real const * const ax, Real * const om);//A.8
	template<typename Real> void ax2ro(Real const * const ax, Real * const ro);//A.9
	template<typename Real> void ax2qu(Real const * const ax, Real * const qu);//A.10
	template<typename Real> void ax2ho(Real const * const ax, Real * const ho);//A.11
	template<typename Real> void ax2cu(Real const * const ax, Real * const cu);//ax->ho->cu
	template<typename Real> void ro2eu(Real const * const ro, Real * const eu);//ro->om->eu
	template<typename Real> void ro2om(Real const * const ro, Real * const om);//ro->ax->om
	template<typename Real> void ro2ax(Real const * const ro, Real * const ax);//A.12
	template<typename Real> void ro2qu(Real const * const ro, Real * const qu);//ro->ax->qu
	template<typename Real> void ro2ho(Real const * const ro, Real * const ho);//A.13
	template<typename Real> void ro2cu(Real const * const ro, Real * const cu);//ro->ho->cu
	template<typename Real> void qu2eu(Real const * const qu, Real * const eu);//A.13
	template<typename Real> void qu2om(Real const * const qu, Real * const om);//A.15
	template<typename Real> void qu2ax(Real const * const qu, Real * const ax);//A.16
	template<typename Real> void qu2ro(Real const * const qu, Real * const ro);//A.17
	template<typename Real> void qu2ho(Real const * const qu, Real * const ho);//A.18
	template<typename Real> void qu2cu(Real const * const qu, Real * const cu);//qu->ho->qu
	template<typename Real> void ho2eu(Real const * const ho, Real * const eu);//ho->ax->eu
	template<typename Real> void ho2om(Real const * const ho, Real * const om);//ho->ax->om
	template<typename Real> void ho2ax(Real const * const ho, Real * const ax);//A.19 [I use Newton's method instead of the table from the paper]
	template<typename Real> void ho2ro(Real const * const ho, Real * const ro);//ho->ax->ro
	template<typename Real> void ho2qu(Real const * const ho, Real * const qu);//ho->ax->qu
	template<typename Real> void ho2cu(Real const * const ho, Real * const qu);
	template<typename Real> void cu2eu(Real const * const cu, Real * const eu);//cu->ho->eu
	template<typename Real> void cu2om(Real const * const cu, Real * const om);//cu->ho->om
	template<typename Real> void cu2ax(Real const * const cu, Real * const ax);//cu->ho->ax
	template<typename Real> void cu2ro(Real const * const cu, Real * const ro);//cu->ho->ro
	template<typename Real> void cu2qu(Real const * const cu, Real * const qu);//cu->ho->qu
	template<typename Real> void cu2ho(Real const * const cu, Real * const qu);

	//@brief   : convert from an orientation matrix to another 
	//@param om: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param yy: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : orientation matrix as 3x3 array in row major order
	template<typename Real> void om2eu(Real const * const * const om, Real * const eu);
	template<typename Real> void om2ax(Real const * const * const om, Real * const ax);
	template<typename Real> void om2ro(Real const * const * const om, Real * const ro);
	template<typename Real> void om2qu(Real const * const * const om, Real * const qu);
	template<typename Real> void om2ho(Real const * const * const om, Real * const ho);
	template<typename Real> void om2cu(Real const * const * const om, Real * const cu);

	//@brief   : convert between 2 different orientation representations
	//@param xx: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param om: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : orientation matrix as 3x3 array in row major order
	template<typename Real> void eu2om(Real const * const eu, Real * const * const om);
	template<typename Real> void ax2om(Real const * const ax, Real * const * const om);
	template<typename Real> void ro2om(Real const * const ro, Real * const * const om);
	template<typename Real> void qu2om(Real const * const qu, Real * const * const om);
	template<typename Real> void ho2om(Real const * const ho, Real * const * const om);
	template<typename Real> void cu2om(Real const * const cu, Real * const * const om);

	//@brief   : convert ZYZ euler angles to quaternion
	//@param eu: euler angles to convert to quaternion (Z, Y', Z'')
	//@param qu: location to write quaternion as w, x, y, z
	//@note    : this is the standard wigner convention
	//@note    : equivalent to eu[0] -= pi/2 and eu[2] += pi/2 followed by eu2qu for ZXZ
	template <typename Real> void zyz2qu(Real const * const eu, Real * const qu);

	//@brief   : convert quaternion to ZYZ euler angles
	//@param qu: quaternion to convert to euler angles as w, x, y, z
	//@param eu: location to write euler angles (Z, Y', Z'')
	//@note    : this is the standard wigner convention
	//@note    : equivalent to qu2eu for ZXZ then eu[0] += pi/2 and eu[2] -= pi/2
	template<typename Real> void qu2zyz(Real const * const qu, Real * const eu);

	//@brief   : convert between 2 different orientation representations
	//@param xx: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param yy: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : xx and yy can be the same memory location
	//@note    : orientation representation abbreviations
	//           eu - ZXZ euler angles (Bunge convention, radians)
	//           om - orientation matrix (9 component vector in row major order)
	//           ax - axis angle (nx, ny, nz, w[radians])
	//           ro - rodrigues vector (nx, ny, nz, tan(w/2))
	//           qu - quaternion (w, x, y, z)
	//           ho - homochoric (x, y, z)
	//           cu - cubochoric (x, y, z)
	template<typename Real> void zyz2eu(Real const * const zyz, Real * const eu );//zyz[0] -= pi/2, zyz[2] += pi/2
	template<typename Real> void zyz2om(Real const * const zyz, Real * const om );//zyz->eu->om
	template<typename Real> void zyz2ax(Real const * const zyz, Real * const ax );//zyz->eu->ax
	template<typename Real> void zyz2ro(Real const * const zyz, Real * const ro );//zyz->eu->ro
	template<typename Real> void zyz2ho(Real const * const zyz, Real * const ho );//zyz->eu->ax->ho
	template<typename Real> void zyz2cu(Real const * const zyz, Real * const cu );//zyz->eu->ho->cu
	template<typename Real> void eu2zyz(Real const * const eu , Real * const zyz);//eu[0] += pi/2, eu[2] -= pi/2
	template<typename Real> void om2zyz(Real const * const om , Real * const zyz);//om->eu->zyz
	template<typename Real> void ax2zyz(Real const * const ax , Real * const zyz);//ax->om->eu->zyz
	template<typename Real> void ro2zyz(Real const * const ro , Real * const zyz);//ro->om->eu->zyz
	template<typename Real> void ho2zyz(Real const * const ho , Real * const zyz);//ho->ax->eu->zyz
	template<typename Real> void cu2zyz(Real const * const cu , Real * const zyz);//cu->ho->eu->zyz

	//enumeration of possible possible rotation representations 
	enum class Rotation {
		Unknown   ,//bad/invalid
		Euler     ,//ZXZ euler angles (Bunge convention, radians)
		Matrix    ,//orientation matrix (9 component vector in row major order)
		AxisAngle ,//axis angle (nx, ny, nz, w[radians])
		Rodrigues ,//rodrigues vector (nx, ny, nz, tan(w/2))
		Quaternion,//quaternion (w, x, y, z)
		Homochoric,//homochoric (x, y, z)
		Cubochoric,//cubochoric (x, y, z)
		EulerZYZ  ,//ZYZ euler angles
	};

	//@brief    : get the number of values in a given representation
	//@param rot: representatino to get length of
	//@return   : length of representation (e.g. 4 for quaternion)
	//@note     : returns 0 for Unknown
	size_t rotLen(const Rotation& rot);

	//@brief: typedef for conversino functions
	//@note : e.g. ConvFunc<double>::type
	template <typename Real> struct ConvFunc {typedef void(*type)(Real const * const, Real * const);};

	//@brief    : get a function pointer to convert between 2 representations
	//@param in : input representation
	//@param out: output representation
	//@return   : function pointer to convert from in => out (NULL if in or out is Unknown)
	template<typename Real> typename ConvFunc<Real>::type getConv(const Rotation& in, const Rotation& out);
}


//@brief    : print a string representation of the rotation type (e.g. 'eu' for Euler)
//@param os : ostream to write to
//@param rot: rotation type to write string representation of
//@return   : os
std::ostream& operator<<(std::ostream& os, const xtal::Rotation& rot);

//@brief    : parse a string representation of the rotation type (e.g. 'eu' for Euler)
//@param is : istream to parse from to write to
//@param rot: rotation type to store result in
//@return   : is
std::istream& operator>>(std::istream& is, xtal::Rotation& rot);

////////////////////////////////////////////////////////////////////////
//                          Implementations                           //
////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>

#include "constants.hpp"

namespace xtal {
	////////////////////////////////////////////////////////////////////////
	//                          Helper Functions                          //
	////////////////////////////////////////////////////////////////////////
	namespace detail {
		
		//@brief: helper class to hold commonly used floating point constants
		template <typename Real>
		struct RotConst {
			//make sure that the Real type has the appropriate attributes
			static_assert(std::is_floating_point<Real>::value, "Real must be a floating point type");//disallow integers

			//constants for ho <--> cu conversion
			static const Real k1;// = (Pi / 6) ^ (1 / 6);
			static const Real k2;// = Pi / 12;
			static const Real k3;// = sqrt(3 / Pi) * 2 ^ (3 / 4);
			static const Real k4;// = sqrt(6 / Pi);
			static const Real k5;// = sqrt(Pi / 24);
			static const Real i1;// = 1 / k1
			static const Real i2;// = 1 / k2
			static const Real i4;// = 1 / k4
		};

		template <typename Real> const Real RotConst<Real>::k1  = Real(0.897772786961286112895391569863263243784555740743629673879928327617808885);
		template <typename Real> const Real RotConst<Real>::k2  = Real(0.261799387799149436538553615273291907016430783281258818414578716025651367);
		template <typename Real> const Real RotConst<Real>::k3  = Real(1.64345640297250301250152526519400784863401933553126716815070358508401456 );
		template <typename Real> const Real RotConst<Real>::k4  = Real(1.38197659788534191706097858412755873255750948130497325394682941957806262 );
		template <typename Real> const Real RotConst<Real>::k5  = Real(0.361800627279133829681507313645397883936054447392272963488106163885184103);

		template <typename Real> const Real RotConst<Real>::i1  = Real(1.11386757821511262144383027557299082539112670911188377842237095371596563 );
		template <typename Real> const Real RotConst<Real>::i2  = Real(3.81971863420548805845321032094034468882703149777095476994401625741352314 );
		template <typename Real> const Real RotConst<Real>::i4  = Real(0.723601254558267659363014627290795767872108894784545926976212327770368205);

		//@brief   : helper function to consistently select rotation axis for rotations of pi
		//@param ax: axis to orient as {x, y, z}
		template <typename Real> void orientAxis(Real * const ax) {
			if(std::fabs(ax[2]) < Constants<Real>::rEps) {//z is 0 -> on equator
				ax[2] = 0;//make z actually zero
				if(std::fabs(ax[1]) < Constants<Real>::rEps) {//y and z are zero
					ax[1] = 0;//make y actually zero
					ax[0] = 1;//use 100 (not -100)
				} else {//y is nonzero
					const Real mag = std::copysign(std::sqrt(ax[0] * ax[0] + ax[1] * ax[1]), ax[1]);//renormalize since z was zerod
					ax[0] /= mag;//normalize and use +y half of equator
					ax[1] /= mag;//normalize and use +y half of equator
				}
			} else {//z is nonzero
				if(std::signbit(ax[2])) std::transform(ax, ax + 3, ax, std::negate<Real>());//use northern hemisphere
			}
		}

		//@brief  : inverse of ((3/4)*(x-sin(x)))^(2/3) via newton's method
		//@param y: ((3/4)*(x-sin(x)))^(2/3) )
		//@return : x
		//@note   : this is the magnitude of the homochoric vector^2
		template <typename Real>
		Real hoInv(Real y) {
			//sanity check input
			if(y < Real(0) || y >= Constants<Real>::hoR2 + Constants<Real>::thr) {
				std::stringstream ss;
				ss << "homochoric magnitude^2 " << y << " is outside of 0, (3*pi/4)^(2/3)]";
				throw std::domain_error(ss.str());
			}

			//handle edges
			if(y >= Constants<Real>::hoR2 - Constants<Real>::thr) return Constants<Real>::pi;
			else if(y <= Constants<Real>::thr) return Real(0);

			//newton iterate
			Real x = Real(2) * std::acos(Real(1) - y / Real(2));//initial guess from small taylor expansion
			y = std::sqrt(y);
			Real prevErr = Constants<Real>::inf;
			for(size_t i = 0; i < 16; i++) {
				const Real fx = std::pow(Real(3) * (x - std::sin(x)) / Real(4), Real(1) / Real(3));//compute forward value
				const Real delta = fx - y;//compute error
				const Real err = std::fabs(delta);
				if(0 == delta || err == prevErr) return x;//no error or flipping between +/- v
				x -= Real(4) * fx * fx * delta / (Real(1) - std::cos(x));//update
				if(err > prevErr) return x;//flipping between v / -2v
				prevErr = err;
			}

			//if we made it this far newton's method failed to converge
			std::stringstream ss;
			ss << "failed to invert ((3/4)*(x-sin(x)))^(2/3) for " << y;
			throw std::runtime_error(ss.str());
			return 0;
		}

		//@brief  : sextant type for cubochoric <-> homochoric transformation symmetry
		//@param v: vector to compute sextant of
		//@return : sextant type
		template<typename Real>
		size_t pyramidType(Real const * const v) {
			std::pair<Real const *, Real const *> minMax = std::minmax_element(v, v+3);
			const Real minEl = std::fabs(*(minMax.first));
			const Real maxEl =  *(minMax.second);
			return std::distance(v, minEl > maxEl ? minMax.first : minMax.second);
		}

		//@brief  : forward and reverse shuffling for cubochoric <-> homochoric transformation symmetry
		//@param v: vector to shuffle
		//@param p: pyramid type (from pyramidType(v))
		//@param u: true to unshuffle (undo a previous call of shufflePyramid(v, p))
		template<typename Real>
		void shufflePyramid(Real * const v, const size_t p, bool u = false) {
			if(p != 2) std::rotate(v, v + (u ? 2-p : p+1), v+3);
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                         Direct Conversions                         //
	////////////////////////////////////////////////////////////////////////

	//@brief   : convert from euler angle to orientation matrix
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@param om: orientation matrix as row major 3x3 array
	//@note    : Rowenhorst et.al. A.1
	template<typename Real> void eu2om(Real const * const eu, Real * const * const om) {
		//build matrix
		const Real c0 = std::cos(eu[0]);
		const Real c1 = std::cos(eu[1]);
		const Real c2 = std::cos(eu[2]);
		const Real s0 = std::sin(eu[0]);
		const Real s1 = std::sin(eu[1]);
		const Real s2 = std::sin(eu[2]);
		om[0][0] =  c0 * c2 - s0 * c1 * s2;
		om[0][1] =  s0 * c2 + c0 * c1 * s2;
		om[0][2] =  s1 * s2;
		om[1][0] = -c0 * s2 - s0 * c1 * c2;
		om[1][1] = -s0 * s2 + c0 * c1 * c2;
		om[1][2] =  s1 * c2;
		om[2][0] =  s0 * s1;
		om[2][1] = -c0 * s1;
		om[2][2] =  c1;

		//zero out tiny entries
		for(size_t i = 0; i < 3; i++) {
			for(size_t j = 0; j < 3; j++) {
				if(std::fabs(om[i][j]) <= Constants<Real>::thr) om[i][j] = Real(0);
			}
		}
	}

	//@brief   : convert from euler angle to axis angle axis angle
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@note    : Rowenhorst et.al. A.2
	template<typename Real> void eu2ax(Real const * const eu, Real * const ax) {
		//compute initial values and handle tiny rotations
		const Real t     = std::tan(eu[1] / Real(2));
		const Real sigma = (eu[0] + eu[2]) / Real(2);
		const Real s     = std::sin(sigma);
		const Real tau   = std::sqrt(t * t + s * s);
		if(tau <= Constants<Real>::thr * Real(2)) {
			std::fill(ax, ax+4, Real(0));
			ax[2] = Real(1);
			return;
		}

		//compute axis
		const Real delta = (eu[0] - eu[2]) / Real(2);
		Real alpha = std::fabs(sigma - Constants<Real>::pi / Real(2)) <= Constants<Real>::thr ? Constants<Real>::pi : Real(2) * std::atan(tau / std::cos(sigma));
		const Real k = Real(-pijk) / std::copysign(tau, alpha);
		ax[0] = k * std::cos(delta) * t;
		ax[1] = k * std::sin(delta) * t;
		ax[2] = k * std::sin(sigma);

		//normalize axis
		const Real mag = std::sqrt(std::inner_product(ax, ax+3, ax, Real(0)));
		std::transform(ax, ax+3, ax, [mag](const Real i){return i/mag;});

		//handle ambiguous case (rotation angle of pi)
		alpha = std::fabs(alpha);
		if(alpha + Constants<Real>::thr >= Constants<Real>::pi) {
			detail::orientAxis(ax);
			ax[3] = Constants<Real>::pi;
		} else if(alpha <= Constants<Real>::thr) {
			std::fill(ax, ax+4, Real(0));
			ax[3] = Real(1);
		} else {
			ax[3] = alpha;
		}
	}

	//@brief   : convert from euler angle to rodrigues
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@param ro: rodrigues vector as (nx, ny, nz, tan(w/2))
	//@note    : Rowenhorst et.al. A.3
	template<typename Real> void eu2ro(Real const * const eu, Real * const ro) {
		eu2ax(eu, ro);
		ro[3] = ro[3] == Constants<Real>::pi ? Constants<Real>::inf : std::tan(ro[3] / Real(2));
	}

	//@brief   : convert from euler angle to quaternion
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@param ro: quaternion as (w,x,y,z)
	//@note    : Rowenhorst et.al. A.4
	template<typename Real> void eu2qu(Real const * const eu, Real * const qu) {
		const Real c = std::cos(eu[1] / Real(2));
		const Real s = std::sin(eu[1] / Real(2));
		const Real sigma = (eu[0] + eu[2]) / Real(2);
		const Real delta = (eu[0] - eu[2]) / Real(2);
		qu[0] = c * std::cos(sigma);
		qu[1] = s * std::cos(delta) * Real(-pijk);
		qu[2] = s * std::sin(delta) * Real(-pijk);
		qu[3] = c * std::sin(sigma) * Real(-pijk);
		if(qu[0] < Real(0)) std::transform(qu, qu+4, qu, std::negate<Real>());

		//normalize
		Real mag = std::sqrt(std::inner_product(qu, qu+4, qu, Real(0)));
		std::transform(qu, qu+4, qu, [mag](const Real i){return i/mag;});

		//handle ambiguous case (rotation angle of pi)
		if(std::abs(qu[0]) <= Constants<Real>::thr) {
			detail::orientAxis(qu+1);
			qu[0] = Real(0);
		}
	}

	//@brief   : convert from orientation matrix to euler angle
	//@param om: orientation matrix as row major 3x3 array
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@note    : Rowenhorst et.al. A.5
	template<typename Real> void om2eu(Real const * const * const om, Real * const eu) {
		if(std::fabs(om[2][2]) >= Real(1) - Constants<Real>::thr) {
			if(om[2][2] > Real(0)) {
				eu[0] =  std::atan2( om[0][1], om[0][0]);//eu = [_, 0, _]
				eu[1] = Real(0);
			} else {
				eu[0] = -std::atan2(-om[0][1], om[0][0]);//eu = [_, pi, _]
				eu[1] = Constants<Real>::pi;
			}
			eu[2] = Real(0);
		} else {
			eu[1] = std::acos(om[2][2]);
			const Real zeta = Real(1) / std::sqrt(Real(1) - om[2][2] * om[2][2]);
			eu[0] = std::atan2(om[2][0] * zeta, -om[2][1] * zeta);
			eu[2] = std::atan2(om[0][2] * zeta,  om[1][2] * zeta);
		}
		std::for_each(eu, eu+3, [](Real& i){if(i < Real(0)) i += Constants<Real>::pi2;});
	}

	//@brief   : convert from orientation matrix to axis angle
	//@param om: orientation matrix as row major 3x3 array
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@note    : Rowenhorst et.al. A.6
	template<typename Real> void om2ax(Real const * const * const om, Real * const ax) {
		Real omega = (om[0][0] + om[1][1] + om[2][2] - Real(1)) / Real(2);
		if(omega >= Real(1) - Constants<Real>::thr) {
			std::fill(ax, ax+4, Real(0));
			ax[2] = Real(1);
			return;
		}

		//compute eigenvector for eigenvalue of 1 (cross product of 2 adjacent columns of A-y*I)
		const Real om00 = om[0][0] - Real(1);
		const Real om11 = om[1][1] - Real(1);
		const Real om22 = om[2][2] - Real(1);
		const Real vecs[3][3] = {
			{om[1][0]*om[2][1] - om[2][0]*  om11  , om[2][0]*om[0][1] -   om00  *om[2][1],   om00  *  om11   - om[1][0]*om[0][1]},
			{  om11  *  om22   - om[2][1]*om[1][2], om[2][1]*om[0][2] - om[0][1]*  om22  , om[0][1]*om[1][2] -   om11  *om[0][2]},
			{om[1][2]*om[2][0] -   om22  *om[1][0],   om22  *  om00   - om[0][2]*om[2][0], om[0][2]*om[1][0] - om[1][2]*  om00  }
		};

		//select vector with largest magnitude
		const Real mags[3] = {
			std::sqrt(std::inner_product(vecs[0], vecs[0]+3, vecs[0], Real(0))),
			std::sqrt(std::inner_product(vecs[1], vecs[1]+3, vecs[1], Real(0))),
			std::sqrt(std::inner_product(vecs[2], vecs[2]+3, vecs[2], Real(0)))
		};
		const size_t i = std::distance(mags, std::max_element(mags, mags+3));
		if(mags[i] <= Constants<Real>::thr) {
			std::fill(ax, ax+4, Real(0));
			ax[2] = Real(1);
			return;
		}
		const Real mag = mags[i];
		std::transform(vecs[i], vecs[i]+3, ax, [mag](const Real i){return i/mag;});

		//handle ambiguous case (rotation of pi)
		if(omega <= Constants<Real>::thr - Real(1)) {
			detail::orientAxis(ax);
			ax[3] = Constants<Real>::pi;
			return;
		}

		//check axis sign
		ax[0] = std::copysign(ax[0], Real(pijk) * (om[2][1] - om[1][2]));
		ax[1] = std::copysign(ax[1], Real(pijk) * (om[0][2] - om[2][0]));
		ax[2] = std::copysign(ax[2], Real(pijk) * (om[1][0] - om[0][1]));
		ax[3] = std::acos(omega);
	}

	//@brief   : convert from orientation matrix to quaternion
	//@param om: orientation matrix as row major 3x3 array
	//@param qu: quaternion as (w, x, y, z)
	//@note    : Rowenhorst et.al. A.7
	template<typename Real> void om2qu(Real const * const * const om, Real * const qu) {
		qu[0] = Real(1) + om[0][0] + om[1][1] + om[2][2];
		qu[1] = Real(1) + om[0][0] - om[1][1] - om[2][2];
		qu[2] = Real(1) - om[0][0] + om[1][1] - om[2][2];
		qu[3] = Real(1) - om[0][0] - om[1][1] + om[2][2];

		//handle ambiguous case (rotation of pi)
		if(std::fabs(qu[0]) <= Constants<Real>::thr) {
			om2ax(om, qu);
			ax2qu(qu, qu);
			return;
		}

		//handle rotation of 0
		if(qu[0] <= Constants<Real>::thr - Real(2)) {
			std::fill(qu+1, qu+4, Real(0));
			qu[0] = Real(1);
			return;
		}

		std::transform(qu, qu+4, qu, [](const Real i){return i <= Constants<Real>::thr ? Real(0) : Real(pijk) * std::sqrt(i) / Real(2);});
		if(Real(pijk) * om[1][2] > Real(pijk) * om[2][1]) qu[1] = -qu[1];
		if(Real(pijk) * om[2][0] > Real(pijk) * om[0][2]) qu[2] = -qu[2];
		if(Real(pijk) * om[0][1] > Real(pijk) * om[1][0]) qu[3] = -qu[3];

		//ensure rotation angle <= pi
		if(qu[0] < Real(0)) std::transform(qu, qu+4, qu, std::negate<Real>());

		//normalize
		Real mag = std::sqrt(std::inner_product(qu, qu+4, qu, Real(0)));
		std::transform(qu, qu+4, qu, [mag](const Real i){return i/mag;});
	}

	//@brief   : convert from axis angle to orientation matrix
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@param om: orientation matrix as row major 3x3 array
	//@note    : Rowenhorst et.al. A.8
	template<typename Real> void ax2om(Real const * const ax, Real * const * const om) {
		if(std::fabs(ax[3]) <= Constants<Real>::thr) {
			for(size_t i = 0; i < 3; i++) {
				std::fill(om[i], om[i]+3, Real(0));
				om[i][i] = Real(1);
			}
		} else if(std::fabs(ax[3]) >= Constants<Real>::pi - Constants<Real>::thr) {
			for(size_t i = 0; i < 3; i++)
				om[i][i] =Real(2) * ax[i] * ax[i] - Real(1);
			for(size_t i = 0; i < 3; i++) {
				const size_t j = (i+1)%3;
				const Real x = Real(2) * ax[i] * ax[j];
				om[i][j] = x;
				om[j][i] = x;
			}
		} else {
			const Real c = std::cos(ax[3]);
			const Real s = std::sin(ax[3]);
			const Real omc = Real(1) - c;
			for(size_t i = 0; i < 3; i++)
				om[i][i] = c + omc * ax[i] * ax[i];
			for(size_t i = 0; i < 3; i++) {
				const size_t j = (i+1)%3;
				const Real x = omc * ax[i] * ax[j];
				const Real y = Real(pijk) * s * ax[(i+2)%3];
				om[i][j] = x - y;
				om[j][i] = x + y;
			}
		}
	}

	//@brief   : convert from axis angle to rodrigues vector
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@param ro: rodrigues vector as (nx, ny, nz, tan(w/2))
	//@note    : Rowenhorst et.al. A.9
	template<typename Real> void ax2ro(Real const * const ax, Real * const ro) {
		if(std::fabs(ax[3]) <= Constants<Real>::thr) {
			std::fill(ro, ro+4, Real(0));
			ro[2] = Real(1);
		} else {
			std::copy(ax, ax+3, ro);
			ro[3] = std::fabs(ax[3] - Constants<Real>::pi) <= Constants<Real>::thr ? Constants<Real>::inf : std::tan(ax[3] / Real(2));
		}
	}

	//@brief   : convert from axis angle to quaternion
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@param qu: quaternion as (w, x, y, z)
	//@note    : Rowenhorst et.al. A.10
	template<typename Real> void ax2qu(Real const * const ax, Real * const qu) {
		if(std::fabs(ax[3]) <= Constants<Real>::thr) {
			std::fill(qu, qu+4, Real(0));
			qu[0] = Real(1);
		} else {
			if(ax == qu) std::rotate(qu, qu+3, qu+4);//rotate_copy doesn't work if source and destination are the same
			else std::rotate_copy(ax, ax+3, ax+4, qu);
			const Real s = std::sin(qu[0] / Real(2));
			qu[0] = std::cos(qu[0] / Real(2));
			std::for_each(qu+1, qu+4, [s](Real& i){i*=s;});

			//normalize
			Real mag = std::sqrt(std::inner_product(qu, qu+4, qu, Real(0)));
			std::for_each(qu, qu+4, [mag](Real& i){i/=mag;});
		}
	}

	//@brief   : convert from axis angle to homochoric vector
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@param ho: homochoric vector as (x, y, z)
	//@note    : Rowenhorst et.al. A.11
	template<typename Real> void ax2ho(Real const * const ax, Real * const ho) {
		if(std::fabs(ax[3]) <= Constants<Real>::thr) {
			std::fill(ho, ho+3, Real(0));
		} else {
			const Real k = std::pow(Real(3) / Real(4) * ( ax[3] - std::sin(ax[3]) ), Real(1) / Real(3));
			std::transform(ax, ax+3, ho, [k](const Real i){return i*k;});
		}
	}

	//@brief   : convert from rodrigues vector to axis angle
	//@param ro: rodrigues vector as (nx, ny, nz, tan(w/2))
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@note    : Rowenhorst et.al. A.12
	template<typename Real> void ro2ax(Real const * const ro, Real * const ax) {
		if(std::fabs(ro[3]) <= Constants<Real>::thr) {
			std::fill(ax, ax+4, Real(0));
			ax[2] = Real(1);
		} else {
			std::copy(ro, ro+3, ax);
			ax[3] = ro[3] < Constants<Real>::inf ? Real(2) * std::atan(ro[3]) : Constants<Real>::pi;
		}
	}

	//@brief   : convert from rodrigues vector to homochoric vector
	//@param ro: rodrigues vector as (nx, ny, nz, tan(w/2))
	//@param ho: homochoric vector as (x, y, z)
	//@note    : Rowenhorst et.al. A.13
	template<typename Real> void ro2ho(Real const * const ro, Real * const ho) {
		const Real t = ro[3] < Constants<Real>::inf ? Real(2) * std::atan(ro[3]) : Constants<Real>::pi;
		if(std::fabs(t) <= Constants<Real>::thr) {
			std::fill(ho, ho+3, Real(0));
		} else {
			const Real k = std::pow(Real(3) * ( t - std::sin(t) ) / Real(4), Real(1) / Real(3));
			std::transform(ro, ro+3, ho, [k](const Real i){return i*k;});
		}
	}

	//@brief   : convert from quaternion to euler angle
	//@param qu: quaternion as (w, x, y, z)
	//@param eu: ZXZ euler angles in radians (Bunge)
	//@note    : Rowenhorst et.al. A.13
	template<typename Real> void qu2eu(Real const * const qu, Real * const eu) {
		const Real q03 = qu[0] * qu[0] + qu[3] * qu[3];
		const Real q12 = qu[1] * qu[1] + qu[2] * qu[2];
		const Real chi = std::sqrt(q03 * q12);
		if(chi <= Constants<Real>::thr) {
			if(q12 <= Constants<Real>::thr){
				eu[0] = std::atan2(Real(-2*pijk) * qu[0] * qu[3], qu[0] * qu[0] - qu[3] * qu[3]);
				eu[1] = Real(0);
			} else {
				eu[0] = std::atan2(Real( 2     ) * qu[1] * qu[2], qu[1] * qu[1] - qu[2] * qu[2]);
				eu[1] = Constants<Real>::pi;
			}
			eu[2] = Real(0);
		} else {
			const Real y1 = qu[1] * qu[3]       ;//can divide by chi (but atan2 is magnitude independent)
			const Real y2 = qu[0] * qu[2] * pijk;//can divide by chi (but atan2 is magnitude independent)
			const Real x1 =-qu[0] * qu[1] * pijk;//can divide by chi (but atan2 is magnitude independent)
			const Real x2 = qu[2] * qu[3]       ;//can divide by chi (but atan2 is magnitude independent)
			eu[0] = std::atan2(y1 - y2, x1 - x2);
			eu[1] = std::atan2(Real(2) * chi, q03 - q12);
			eu[2] = std::atan2(y1 + y2, x1 + x2);
		}
		std::for_each(eu, eu+3, [](Real& i){if(i < Real(0)) i += Constants<Real>::pi2;});
	}

	//@brief   : convert from quaternion to orientation matrix
	//@param qu: quaternion as (w, x, y, z)
	//@param om: orientation matrix as row major 3x3 array
	//@note    : Rowenhorst et.al. A.15
	template<typename Real> void qu2om(Real const * const qu, Real * const * const om) {
		const Real qbar = qu[0] * qu[0] - qu[1] * qu[1] - qu[2] * qu[2] - qu[3] * qu[3];
		for(size_t i = 0; i < 3; i++) {
			om[i][i] = qbar + Real(2) * qu[i+1] * qu[i+1];

			const size_t j = (i+1)%3;
			const Real x = Real(2) * qu[i+1] * qu[j+1];
			const Real y = Real(2) * Real(pijk) * qu[0] * qu[(i+2)%3+1];
			om[i][j] = x - y;
			om[j][i] = x + y;
		}
	}

	//@brief   : convert from quaternion to axis angle
	//@param qu: quaternion as (w, x, y, z)
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@note    : Rowenhorst et.al. A.16
	template<typename Real> void qu2ax(Real const * const qu, Real * const ax) {
		const Real omega = Real(2) * std::acos(qu[0]);
		if(omega <= Constants<Real>::thr) {
			std::fill(ax, ax+4, Real(0));
			ax[2] = Real(1);
		} else {
			const Real s = std::copysign(Real(1) / std::sqrt(std::inner_product(qu+1, qu+4, qu+1, Real(0))), qu[0]);
			std::transform(qu+1, qu+4, ax, [s](const Real i){return i*s;});
			if(omega >= Constants<Real>::pi - Constants<Real>::thr) {
				detail::orientAxis(ax);
				ax[3] = Constants<Real>::pi;			
			} else {
				ax[3] = omega;
			}
		}
	}

	//@brief   : convert from quaternion to rodrigues vector
	//@param qu: quaternion as (w, x, y, z)
	//@param ro: rodrigues vector as (nx, ny, nz, tan(w/2))
	//@note    : Rowenhorst et.al. A.17
	template<typename Real> void qu2ro(Real const * const qu, Real * const ro) {
		if(qu[0] <= Constants<Real>::thr) {
			std::copy(qu+1, qu+4, ro);
			ro[3] = Constants<Real>::inf;
			return;
		}
		const Real s = std::sqrt(std::inner_product(qu+1, qu+4, qu+1, Real(0)));
		if(s <= Constants<Real>::thr) {
			std::fill(ro, ro+4, Real(0));
			ro[2] = Real(1);
		} else {
			std::transform(qu+1, qu+4, ro, [s](const Real i){return i/s;});
			ro[3] = std::tan(std::acos(qu[0]));
		}
	}

	//@brief   : convert from quaternion to homochoric vector
	//@param qu: quaternion as (w, x, y, z)
	//@param ho: homochoric vector as (x, y, z)
	//@note    : Rowenhorst et.al. A.18
	template<typename Real> void qu2ho(Real const * const qu, Real * const ho) {
		const Real omega = Real(2) * std::acos(qu[0]);
		if(std::fabs(omega) <= Constants<Real>::thr) {
			std::fill(ho, ho+3, Real(0));
		} else {
			const Real s = Real(1) / std::sqrt(std::inner_product(qu+1, qu+4, qu+1, Real(0)));
			const Real k = std::pow(Real(3) * ( omega - std::sin(omega) ) / Real(4), Real(1) / Real(3)) * s;
			std::transform(qu+1, qu+4, ho, [k](const Real i){return i*k;});
		}
	}

	//@brief   : convert from homochoric vector to axis angle
	//@param ho: homochoric vector as (x, y, z)
	//@param ax: axis angle pair as (nx, ny, nz, w [radians])
	//@note    : Rowenhorst et.al. A.19
	template<typename Real> void ho2ax(Real const * const ho, Real * const ax) {
		const Real mag2 = std::inner_product(ho, ho+3, ho, Real(0));
		if(mag2 <= Constants<Real>::thr) {
				std::fill(ax, ax+4, Real(0));
				ax[2] = Real(1);
		} else {
			ax[3] = detail::hoInv(mag2);
			const Real mag = std::sqrt(mag2);
			std::transform(ho, ho+3, ax, [mag](const Real i){return i/mag;});
			if(ax[3] >= Constants<Real>::pi - Constants<Real>::thr) {
				detail::orientAxis(ax);
				ax[3] = Constants<Real>::pi;
			}
		}
	}

	////////////////////////////////////////////////////////////////////////
	//                        Indirect Conversions                        //
	////////////////////////////////////////////////////////////////////////

	//no temporary storage required
	template<typename Real> void eu2cu(Real const * const         eu, Real * const         cu) {            eu2ho(eu, cu); ho2cu(cu, cu);}
	template<typename Real> void om2ro(Real const * const * const om, Real * const         ro) {            om2eu(om, ro); eu2ro(ro, ro);}
	template<typename Real> void om2cu(Real const * const * const om, Real * const         cu) {            om2ho(om, cu); ho2cu(cu, cu);}
	template<typename Real> void ax2cu(Real const * const         ax, Real * const         cu) {            ax2ho(ax, cu); ho2cu(cu, cu);}
	template<typename Real> void ro2qu(Real const * const         ro, Real * const         qu) {            ro2ax(ro, qu); ax2qu(qu, qu);}
	template<typename Real> void ro2cu(Real const * const         ro, Real * const         cu) {            ro2ho(ro, cu); ho2cu(cu, cu);}
	template<typename Real> void qu2cu(Real const * const         qu, Real * const         cu) {            qu2ho(qu, cu); ho2cu(cu, cu);}
	template<typename Real> void ho2ro(Real const * const         ho, Real * const         ro) {            ho2ax(ho, ro); ax2ro(ro, ro);}
	template<typename Real> void ho2qu(Real const * const         ho, Real * const         qu) {            ho2ax(ho, qu); ax2qu(qu, qu);}
	template<typename Real> void cu2eu(Real const * const         cu, Real * const         eu) {            cu2ho(cu, eu); ho2eu(eu, eu);}
	template<typename Real> void cu2qu(Real const * const         cu, Real * const         qu) {            cu2ho(cu, qu); ho2qu(qu, qu);}
	template<typename Real> void cu2ax(Real const * const         cu, Real * const         ax) {            cu2ho(cu, ax); ho2ax(ax, ax);}

	//temporary storage required
	template<typename Real> void eu2ho(Real const * const         eu, Real * const         ho) {Real ax[4]; eu2ax(eu, ax); ax2ho(ax, ho);}
	template<typename Real> void om2ho(Real const * const * const om, Real * const         ho) {Real ax[4]; om2ax(om, ax); ax2ho(ax, ho);}
	template<typename Real> void ax2eu(Real const * const         ax, Real * const         eu) {Real om[9]; ax2om(ax, om); om2eu(om, eu);}
	template<typename Real> void ro2eu(Real const * const         ro, Real * const         eu) {Real om[9]; ro2om(ro, om); om2eu(om, eu);}
	template<typename Real> void ro2om(Real const * const         ro, Real * const * const om) {Real ax[4]; ro2ax(ro, ax); ax2om(ax, om);}
	template<typename Real> void ho2eu(Real const * const         ho, Real * const         eu) {Real ax[4]; ho2ax(ho, ax); ax2eu(ax, eu);}
	template<typename Real> void ho2om(Real const * const         ho, Real * const * const om) {Real ax[4]; ho2ax(ho, ax); ax2om(ax, om);}
	template<typename Real> void cu2om(Real const * const         cu, Real * const * const om) {Real ho[3]; cu2ho(cu, ho); ho2om(ho, om);}
	template<typename Real> void cu2ro(Real const * const         cu, Real * const         ro) {Real ho[3]; cu2ho(cu, ho); ho2ro(ho, ro);}

	//zyz conversions
	template<typename Real> void zyz2om(Real const * const zyz, Real * const om) {zyz2eu(zyz, om); eu2om(om, om);}//zyz->eu->om
	template<typename Real> void zyz2ax(Real const * const zyz, Real * const ax) {zyz2eu(zyz, ax); eu2ax(ax, ax);}//zyz->eu->ax
	template<typename Real> void zyz2ro(Real const * const zyz, Real * const ro) {zyz2eu(zyz, ro); eu2ro(ro, ro);}//zyz->eu->ro
	template<typename Real> void zyz2ho(Real const * const zyz, Real * const ho) {zyz2eu(zyz, ho); eu2ho(ho, ho);}//zyz->eu->ax->ho
	template<typename Real> void zyz2cu(Real const * const zyz, Real * const cu) {zyz2eu(zyz, cu); eu2cu(cu, cu);}//zyz->eu->ho->cu

	template<typename Real> void om2zyz(Real const * const om, Real * const zyz) {om2eu(om, zyz); eu2zyz(zyz, zyz);}//om->eu->zyz
	template<typename Real> void ax2zyz(Real const * const ax, Real * const zyz) {ax2eu(ax, zyz); eu2zyz(zyz, zyz);}//ax->om->eu->zyz
	template<typename Real> void ro2zyz(Real const * const ro, Real * const zyz) {ro2eu(ro, zyz); eu2zyz(zyz, zyz);}//ro->om->eu->zyz
	template<typename Real> void ho2zyz(Real const * const ho, Real * const zyz) {ho2eu(ho, zyz); eu2zyz(zyz, zyz);}//ho->ax->eu->zyz
	template<typename Real> void cu2zyz(Real const * const cu, Real * const zyz) {cu2eu(cu, zyz); eu2zyz(zyz, zyz);}//cu->ho->eu->zyz


	////////////////////////////////////////////////////////////////////////
	//                       Vectorized OM Wrappers                       //
	////////////////////////////////////////////////////////////////////////

	//@brief   : convert from an orientation matrix to another 
	//@param om: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param yy: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : orientation matrix as 3x3 array in row major order
	// template<typename Real> void om2eu(Real const * const om, Real * const eu) {om2xx<Real>(&om2eu<Real>, om, eu);}
	// template<typename Real> void om2ax(Real const * const om, Real * const ax) {om2xx<Real>(&om2ax<Real>, om, ax);}
	// template<typename Real> void om2ro(Real const * const om, Real * const ro) {om2xx<Real>(&om2ro<Real>, om, ro);}
	// template<typename Real> void om2qu(Real const * const om, Real * const qu) {om2xx<Real>(&om2qu<Real>, om, qu);}
	// template<typename Real> void om2ho(Real const * const om, Real * const ho) {om2xx<Real>(&om2ho<Real>, om, ho);}
	// template<typename Real> void om2cu(Real const * const om, Real * const cu) {om2xx<Real>(&om2cu<Real>, om, cu);}
	
	template<typename Real> void om2eu(Real const * const om, Real * const eu) {Real const * const m[3] = {om, om+3, om+6}; om2eu(m, eu);}
	template<typename Real> void om2ax(Real const * const om, Real * const ax) {Real const * const m[3] = {om, om+3, om+6}; om2ax(m, ax);}
	template<typename Real> void om2ro(Real const * const om, Real * const ro) {Real const * const m[3] = {om, om+3, om+6}; om2ro(m, ro);}
	template<typename Real> void om2qu(Real const * const om, Real * const qu) {Real const * const m[3] = {om, om+3, om+6}; om2qu(m, qu);}
	template<typename Real> void om2ho(Real const * const om, Real * const ho) {Real const * const m[3] = {om, om+3, om+6}; om2ho(m, ho);}
	template<typename Real> void om2cu(Real const * const om, Real * const cu) {Real const * const m[3] = {om, om+3, om+6}; om2cu(m, cu);}
	
	//@brief   : convert between 2 different orientation representations
	//@param xx: input  orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@param om: output orientation (one of eu, om, ax, ro, qu, ho, or cu)
	//@note    : orientation matrix as 3x3 array in row major order
	template<typename Real> void eu2om(Real const * const eu, Real * const om) {Real * const m[3] = {om, om+3, om+6}; eu2om(eu, m);}
	template<typename Real> void ax2om(Real const * const ax, Real * const om) {Real * const m[3] = {om, om+3, om+6}; ax2om(ax, m);}
	template<typename Real> void ro2om(Real const * const ro, Real * const om) {Real * const m[3] = {om, om+3, om+6}; ro2om(ro, m);}
	template<typename Real> void qu2om(Real const * const qu, Real * const om) {Real * const m[3] = {om, om+3, om+6}; qu2om(qu, m);}
	template<typename Real> void ho2om(Real const * const ho, Real * const om) {Real * const m[3] = {om, om+3, om+6}; ho2om(ho, m);}
	template<typename Real> void cu2om(Real const * const cu, Real * const om) {Real * const m[3] = {om, om+3, om+6}; cu2om(cu, m);}

	////////////////////////////////////////////////////////////////////////
	//                       Cubochoric Conversions                       //
	////////////////////////////////////////////////////////////////////////

	//@brief   : convert from cubochoric vector to homochoric vector
	//@param cu: cubochoric vector as (x, y, z)
	//@param ho: homochoric vector as (x, y, z)
	//@note    : Rosca, D., et al
	template<typename Real> void cu2ho(Real const * const cu, Real * const ho) {
		//get pyramid type, shuffle coordinates to +z pyramid, and check bounds
		const size_t p = detail::pyramidType(cu);
		std::copy(cu, cu+3, ho);
		detail::shufflePyramid(ho, p);
		if(std::fabs(ho[2]) >= Constants<Real>::cuA / Real(2) + Constants<Real>::thr) {
			std::stringstream ss;
			ss << "cubochoric component " << ho[2] << " is outside of +/-pi^(2/3)";
			throw std::domain_error(ss.str());
		}

		//handle origin
		if(std::fabs(ho[2]) <= Constants<Real>::thr) {
			std::fill(ho, ho+3, Real(0));
			return;
		}

		//operation M1
		std::transform(ho, ho+3, ho, [](const Real i){return i*detail::RotConst<Real>::k1;});

		//operation M2
		bool swapped = std::fabs(ho[0]) > std::fabs(ho[1]);
		if(swapped) std::swap(ho[0], ho[1]);
		if(std::fabs(ho[1]) >= Constants<Real>::thr) {//skip points along z axis to avoid divide by zero
			const Real theta = detail::RotConst<Real>::k2 * ho[0] / ho[1];
			const Real k = detail::RotConst<Real>::k3 * ho[1] / std::sqrt(Constants<Real>::r2 - std::cos(theta));
			ho[0] = k * Constants<Real>::r2 * std::sin(theta);
			ho[1] = k * Constants<Real>::r2 * std::cos(theta) - k;
			if(swapped) std::swap(ho[0], ho[1]);
		} else {
			std::fill(ho, ho+2, Real(0));
		}

		//operation M3
		const Real ko = std::inner_product(ho, ho+2, ho, Real(0));
		const Real k = std::sqrt(Real(1) - Constants<Real>::pi * ko / (Real(24) * ho[2]*ho[2]));
		std::transform(ho, ho+2, ho, [k](const Real i){return i*k;});
		ho[2] = detail::RotConst<Real>::k4 * ho[2] - ko * detail::RotConst<Real>::k5 / ho[2];

		//unshuffle
		detail::shufflePyramid(ho, p, true);
	}
	//@brief   : convert from homochoric vector to cubochoric vector
	//@param ho: homochoric vector as (x, y, z)
	//@param cu: cubochoric vector as (x, y, z)
	//@note    : Rosca, D., et al
	template<typename Real> void ho2cu(Real const * const ho, Real * const cu) {
		//check bounds, get pyramid type, and shuffle coordinates to +z pyramid
		std::copy(ho, ho+3, cu);
		const Real rs = std::sqrt(std::inner_product(cu, cu+3, cu, Real(0)));
		if(rs >= Constants<Real>::hoR + Constants<Real>::thr) {
			std::stringstream ss;
			ss << "homochoric magnitude " << rs << " is outside of 0, (3*pi/4)^(2/3)]";
			throw std::domain_error(ss.str());
		}
		const size_t p = detail::pyramidType(cu);
		detail::shufflePyramid(cu, p);

		//handle origin
		if(rs <= Constants<Real>::thr) {
			std::fill(cu, cu+3, Real(0));
			return;
		}

		//invert operation M3
		const Real k = std::sqrt(Real(2) * rs / (rs + std::fabs(cu[2])));
		std::transform(cu, cu+2, cu, [k](const Real i){return i*k;});
		cu[2] = std::copysign(rs * detail::RotConst<Real>::i4, cu[2]);

		//invert operation M2
		Real x2 = cu[0] * cu[0];
		Real y2 = cu[1] * cu[1];
		const Real mag2 = x2 + y2;
		if(mag2 >= Constants<Real>::thr) {//skip points along z axis
			//only handle x <= y
			bool swapped = std::fabs(cu[0]) > std::fabs(cu[1]);
			if(swapped) {
				std::swap(cu[0], cu[1]);
				std::swap(x2, y2);
			}
			const Real x2y = mag2 + y2;
			const Real rx2y = std::sqrt(x2y);
			const Real kxy = std::sqrt( detail::RotConst<Real>::k2 * (x2y * mag2) / (x2y - std::fabs(cu[1]) * rx2y) );
			const Real ckx = (x2 + std::fabs(cu[1]) * rx2y) / Constants<Real>::r2 / mag2;
			std::transform(cu, cu+2, cu, [kxy](const Real i){return std::copysign(kxy, i);});
			if(ckx >= Real(1))
				cu[0] = Real(0);
			else if (ckx <= Real(-1))
				cu[0] *= Real(12);
			else
				cu[0] *= std::acos(ckx) * detail::RotConst<Real>::i2;
			if(swapped) std::swap(cu[0], cu[1]);
		} else {
			std::fill(cu, cu+2, Real(0));
		}

		//invert operation M1
		std::transform(cu, cu+3, cu, [](const Real i){return i*detail::RotConst<Real>::i1;});

		//unshuffle
		detail::shufflePyramid(cu, p, true);
	}

	////////////////////////////////////////////////////////////////////////
	//                   Nonstandard Euler Conversions                    //
	////////////////////////////////////////////////////////////////////////

	//@brief   : convert ZYZ euler angles to quaternion
	//@param eu: euler angles to convert to quaternion (Z, Y', Z'')
	//@param qu: location to write quaternion as w, x, y, z
	//@note    : this is the standard wigner convention
	//@note    : equivalent to eu[0] -= pi/2 and eu[2] += pi/2 followed by eu2qu for ZXZ
	//@note    : derived by simplifying [cos(a/2), 0, 0, sin(a/2)*pijk] * [cos(b/2), 0, sin(b/2)*pijk, 0] * [cos(y/2), 0, 0, sin(y/2)*pijk]
	template <typename Real> void zyz2qu(Real const * const eu, Real * const qu) {
		const Real c = std::cos(eu[1] / Real(2));
		const Real s = std::sin(eu[1] / Real(2));
		const Real sigma = (eu[2] + eu[0]) / Real(2);
		const Real delta = (eu[2] - eu[0]) / Real(2);
		qu[0] = c * std::cos(sigma);
		qu[1] = s * std::sin(delta) * Real(-pijk);
		qu[2] = s * std::cos(delta) * Real(-pijk);
		qu[3] = c * std::sin(sigma) * Real(-pijk);
		if(std::signbit(qu[0])) std::transform(qu, qu+4, qu, std::negate<Real>());//restrict rotation to [0,pi]

		//handle ambiguous case (rotation angle of pi)
		if(std::fabs(qu[0]) <= Constants<Real>::rEps) {//rotation is pi, restrict axis to +z hemisphere (+y semicircle on equator and 100 not -100)
			qu[0] = 0;//zero out w
			detail::orientAxis(qu+1);
		}
	}

	//@brief   : convert quaternion to ZYZ euler angles
	//@param qu: quaternion to convert to euler angles as w, x, y, z
	//@param eu: location to write euler angles (Z, Y', Z'')
	//@note    : this is the standard wigner convention
	//@note    : equivalent to qu2eu for ZXZ then eu[0] += pi/2 and eu[2] -= pi/2
	template<typename Real> void qu2zyz(Real const * const qu, Real * const eu) {
		const Real qu0 = qu[0] * Real(pijk);
		const Real q03 = qu[0] * qu[0] + qu[3] * qu[3];
		const Real q12 = qu[1] * qu[1] + qu[2] * qu[2];
		const Real chi = std::sqrt(q03 * q12);
		if(chi <= Constants<Real>::thr) {
			if(q12 <= Constants<Real>::thr){
				eu[0] = std::atan2(-Real(2) * qu0   * qu[3], qu[0] * qu[0] - qu[3] * qu[3]);
				eu[1] = Real(0);
			} else {
				eu[0] = std::atan2(-Real(2) * qu[1] * qu[2], qu[2] * qu[2] - qu[1] * qu[1]);
				eu[1] = Constants<Real>::pi;
			}
			eu[2] = Real(0);
		} else {
			const Real y1 = qu[2] * qu[3];//can divide by chi (but atan2 is magnitude independent)
			const Real y2 =-qu[1] * qu0  ;//can divide by chi (but atan2 is magnitude independent)
			const Real x1 =-qu[2] * qu0  ;//can divide by chi (but atan2 is magnitude independent)
			const Real x2 =-qu[1] * qu[3];//can divide by chi (but atan2 is magnitude independent)
			eu[0] = std::atan2(y1 - y2, x1 - x2);
			eu[1] = std::atan2(Real(2) * chi, q03 - q12);
			eu[2] = std::atan2(y1 + y2, x1 + x2);
		}
		std::for_each(eu, eu+3, [](Real& i){if(i < Real(0)) i += Constants<Real>::pi2;});
	}

	//@brief: convert from ZYZ to ZXZ (Bunge) euler angles
	//@param zyz: ZYZ euler angles
	//@param eu : where to write ZXZ euler angles
	template<typename Real> void zyz2eu(Real const * const zyz, Real * const eu ) {
		eu[0] = zyz[0] - Constants<Real>::pi_2;
		eu[1] = zyz[1]                        ;
		eu[2] = zyz[2] + Constants<Real>::pi_2;
	}

	//@brief: convert from ZYZ to ZXZ (Bunge) euler angles
	//@param eu : ZXZ euler angles
	//@param zyz: where to write ZYZ euler angles
	template<typename Real> void eu2zyz(Real const * const eu , Real * const zyz) {
		zyz[0] = eu[0] + Constants<Real>::pi_2;
		zyz[1] = eu[1]                        ;
		zyz[2] = eu[2] - Constants<Real>::pi_2;
	}

	//@brief    : get the number of values in a given representation
	//@param rot: representatino to get length of
	//@return   : length of representation (e.g. 4 for quaternion)
	//@note     : returns 0 for Unknown
	size_t rotLen(const Rotation& rot) {
		switch(rot) {
			case Rotation::Unknown   : return 0;
			case Rotation::Euler     : return 3;
			case Rotation::Matrix    : return 9;
			case Rotation::AxisAngle : return 4;
			case Rotation::Rodrigues : return 4;
			case Rotation::Quaternion: return 4;
			case Rotation::Homochoric: return 3;
			case Rotation::Cubochoric: return 3;
			case Rotation::EulerZYZ  : return 3;
		}
		throw std::logic_error("unhandled rotation type");
	}

	//@brief    : get a function pointer to convert between 2 representations
	//@param in : input representation
	//@param out: output representation
	//@return   : function pointer to convert from in => out (NULL if in or out is Unknown)
	template<typename Real> typename ConvFunc<Real>::type getConv(const Rotation& in, const Rotation& out) {
		//build functions to copy various sizes (for e.g. eu2eu)
		static const typename ConvFunc<Real>::type cpy3 = [](Real const * const in, Real * const out){std::copy(in, in+3, out);};
		static const typename ConvFunc<Real>::type cpy4 = [](Real const * const in, Real * const out){std::copy(in, in+4, out);};
		static const typename ConvFunc<Real>::type cpy9 = [](Real const * const in, Real * const out){std::copy(in, in+9, out);};

		switch(in) {
			case Rotation::Unknown   :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return NULL;
					case Rotation::Matrix    : return NULL;
					case Rotation::AxisAngle : return NULL;
					case Rotation::Rodrigues : return NULL;
					case Rotation::Quaternion: return NULL;
					case Rotation::Homochoric: return NULL;
					case Rotation::Cubochoric: return NULL;
					case Rotation::EulerZYZ  : return NULL;
				}

			case Rotation::Euler     :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return cpy3;
					case Rotation::Matrix    : return eu2om <Real>;
					case Rotation::AxisAngle : return eu2ax <Real>;
					case Rotation::Rodrigues : return eu2ro <Real>;
					case Rotation::Quaternion: return eu2qu <Real>;
					case Rotation::Homochoric: return eu2ho <Real>;
					case Rotation::Cubochoric: return eu2cu <Real>;
					case Rotation::EulerZYZ  : return eu2zyz<Real>;
				}

			case Rotation::Matrix    :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return om2eu <Real>;
					case Rotation::Matrix    : return cpy9;
					case Rotation::AxisAngle : return om2ax <Real>;
					case Rotation::Rodrigues : return om2ro <Real>;
					case Rotation::Quaternion: return om2qu <Real>;
					case Rotation::Homochoric: return om2ho <Real>;
					case Rotation::Cubochoric: return om2cu <Real>;
					case Rotation::EulerZYZ  : return om2zyz<Real>;
				}

			case Rotation::AxisAngle :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return ax2eu <Real>;
					case Rotation::Matrix    : return ax2om <Real>;
					case Rotation::AxisAngle : return cpy4;
					case Rotation::Rodrigues : return ax2ro <Real>;
					case Rotation::Quaternion: return ax2qu <Real>;
					case Rotation::Homochoric: return ax2ho <Real>;
					case Rotation::Cubochoric: return ax2cu <Real>;
					case Rotation::EulerZYZ  : return ax2zyz<Real>;
				}

			case Rotation::Rodrigues :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return ro2eu <Real>;
					case Rotation::Matrix    : return ro2om <Real>;
					case Rotation::AxisAngle : return ro2ax <Real>;
					case Rotation::Rodrigues : return cpy4;
					case Rotation::Quaternion: return ro2qu <Real>;
					case Rotation::Homochoric: return ro2ho <Real>;
					case Rotation::Cubochoric: return ro2cu <Real>;
					case Rotation::EulerZYZ  : return ro2zyz<Real>;
				}

			case Rotation::Quaternion:
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return qu2eu <Real>;
					case Rotation::Matrix    : return qu2om <Real>;
					case Rotation::AxisAngle : return qu2ax <Real>;
					case Rotation::Rodrigues : return qu2ro <Real>;
					case Rotation::Quaternion: return cpy4;
					case Rotation::Homochoric: return qu2ho <Real>;
					case Rotation::Cubochoric: return qu2cu <Real>;
					case Rotation::EulerZYZ  : return qu2zyz<Real>;
				}

			case Rotation::Homochoric:
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return ho2eu <Real>;
					case Rotation::Matrix    : return ho2om <Real>;
					case Rotation::AxisAngle : return ho2ax <Real>;
					case Rotation::Rodrigues : return ho2ro <Real>;
					case Rotation::Quaternion: return ho2qu <Real>;
					case Rotation::Homochoric: return cpy3;
					case Rotation::Cubochoric: return ho2cu <Real>;
					case Rotation::EulerZYZ  : return ho2zyz<Real>;
				}

			case Rotation::Cubochoric:
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return cu2eu <Real>;
					case Rotation::Matrix    : return cu2om <Real>;
					case Rotation::AxisAngle : return cu2ax <Real>;
					case Rotation::Rodrigues : return cu2ro <Real>;
					case Rotation::Quaternion: return cu2qu <Real>;
					case Rotation::Homochoric: return cu2ho <Real>;
					case Rotation::Cubochoric: return cpy3;
					case Rotation::EulerZYZ  : return cu2zyz<Real>;
				}

			case Rotation::EulerZYZ  :
				switch(out) {
					case Rotation::Unknown   : return NULL;
					case Rotation::Euler     : return zyz2eu <Real>;
					case Rotation::Matrix    : return zyz2om <Real>;
					case Rotation::AxisAngle : return zyz2ax <Real>;
					case Rotation::Rodrigues : return zyz2ro <Real>;
					case Rotation::Quaternion: return zyz2qu <Real>;
					case Rotation::Homochoric: return zyz2ho <Real>;
					case Rotation::Cubochoric: return zyz2cu <Real>;
					case Rotation::EulerZYZ  : return cpy3;
				}
		}
	}
}

//@brief    : print a string representation of the rotation type
//@param os : ostream to write to
//@param rot: rotation type to write string representation of
//@return   : os
std::ostream& operator<<(std::ostream& os, const xtal::Rotation& rot) {
	switch(rot) {
		case xtal::Rotation::Unknown   : return os << "unk";
		case xtal::Rotation::Euler     : return os << "eu";
		case xtal::Rotation::Matrix    : return os << "om";
		case xtal::Rotation::AxisAngle : return os << "ax";
		case xtal::Rotation::Rodrigues : return os << "ro";
		case xtal::Rotation::Quaternion: return os << "qu";
		case xtal::Rotation::Homochoric: return os << "ho";
		case xtal::Rotation::Cubochoric: return os << "cu";
		case xtal::Rotation::EulerZYZ  : return os << "zyz";
	}
	throw std::logic_error("unhandled rotation type");
}

//@brief    : parse a string representation of the rotation type (e.g. 'eu' for Euler)
//@param is : istream to parse from to write to
//@param rot: rotation type to store result in
//@return   : is
std::istream& operator>>(std::istream& is, xtal::Rotation& rot) {
	char nm[2];
	is >> nm[0] >> nm[1];
	if     ('e' == nm[0] && 'u' == nm[1]) rot = xtal::Rotation::Euler     ;
	else if('o' == nm[0] && 'm' == nm[1]) rot = xtal::Rotation::Matrix    ;
	else if('a' == nm[0] && 'x' == nm[1]) rot = xtal::Rotation::AxisAngle ;
	else if('r' == nm[0] && 'o' == nm[1]) rot = xtal::Rotation::Rodrigues ;
	else if('q' == nm[0] && 'u' == nm[1]) rot = xtal::Rotation::Quaternion;
	else if('h' == nm[0] && 'o' == nm[1]) rot = xtal::Rotation::Homochoric;
	else if('c' == nm[0] && 'u' == nm[1]) rot = xtal::Rotation::Cubochoric;
	else if('z' == nm[0] && 'y' == nm[1]) {
		is >> nm[1];
		if('z' == nm[1]) {
			rot = xtal::Rotation::EulerZYZ  ;
		}else {
			rot = xtal::Rotation::Unknown   ;
		}
	} else {
		rot = xtal::Rotation::Unknown;
	}
	return is;
}

#endif//_ROTATION_H_
